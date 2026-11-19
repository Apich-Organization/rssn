//! Checked native device objects and explicit normalized machine IR.
//! Physical map execution uses independent input/parameter/output references.
use super::native_backend_extended::TapeSnapshot;
use super::native_backend_extended::{
    self as execution,
};
use super::native_object_api::Dynamic;
use super::native_object_api::Handle;
use super::native_object_api::Input;
use super::native_object_api::Many;
use super::native_object_api::Maybe;
use super::native_object_api::NativeType;
use super::native_object_api::NewObject;
use super::native_object_api::object_binding;
use super::native_object_api::{
    self as api,
};
use crate::backend::device::CompiledKernel;
use crate::backend::device::Device;
use crate::backend::device::Dtype;
use crate::backend::device::KernelBuilder;
use crate::backend::device::KernelIr;
use crate::backend::device::KernelOp;
use crate::backend::device::ReduceOp;
use crate::backend::device::Storage;
use crate::backend::device::emit::EmitOptions;
use crate::backend::device::emit::Precision;
use crate::backend::device::emit::Target;
use crate::backend::device::emit::{
    self as native_emit,
};
use crate::backend::device::kernel::Columns;
use crate::backend::device::kernel::Scratch;
use crate::backend::device::kernel::{
    self as k,
};
use crate::backend::device::wgsl::KernelEmitter;
use crate::backend::device::wgsl::{
    self,
};
use crate::backend::device::{
    self as d,
};
use crate::backend::tape::Cmp;
use crate::graph::RuleError;
use crate::graph::rule::Installer;
use crate::kernel::OpError;
use crate::kernel::Value;
use crate::kernel::signature::DimExpr;
use crate::kernel::signature::TypeExpr;
use crate::kernel::value::TensorValue;
use crate::rules::numerical_abi::codecs::Output;
use std::sync::Arc;
const LIMIT: usize = 262_144;
struct DeviceObject(Arc<dyn Device>);
struct BufferObject(d::Buffer);
struct KernelObject {
    ir: KernelIr,
    frame: Option<TapeSnapshot>,
}
struct Builder(Option<KernelBuilder>);
struct Compiled {
    native: Arc<dyn CompiledKernel>,
    device: Arc<dyn Device>,
    frame: Option<TapeSnapshot>,
}
struct Registry(d::DeviceRegistry);
struct Op(KernelOp);
struct Options(EmitOptions);
struct DenseStorage(Vec<f64>);
impl Storage for DenseStorage {
    fn as_any(&self) -> &dyn std::any::Any {
        self
    }

    fn as_any_mut(&mut self) -> &mut dyn std::any::Any {
        self
    }
}
macro_rules! native {
    ($t:ty,$kind:literal,$bytes:expr,$body:expr) => {
        impl NativeType for $t {
            const KIND: &'static str = $kind;
            const MAX_BYTES: usize = $bytes;

            fn validate(&self) -> Result<(), OpError> {
                ($body)(self)
            }
        }
    };
}
mod extra;
mod features;
native!(
    DeviceObject,
    "rssn/native/device/v1",
    8_388_608,
    |_: &DeviceObject| Ok(())
);
native!(
    Registry,
    "rssn/native/device-registry/v1",
    16_777_216,
    |r: &Registry| {
        if r.0.list().len() <= 2 {
            Ok(())
        } else {
            Err(api::error(
                "owned device registry supports at most2 devices; global discovery lists all devices",
            ))
        }
    }
);
native!(
    BufferObject,
    "rssn/native/device-buffer/v1",
    2_097_152,
    |b: &BufferObject| {
        if b.0.len() <= LIMIT {
            Ok(())
        } else {
            Err(api::error("device buffer exceeds element budget"))
        }
    }
);
native!(
    KernelObject,
    "rssn/native/device-ir/v1",
    1_048_576,
    |b: &KernelObject| {
        budget(&b.ir)?;
        b.ir.validate().map_err(api::error)
    }
);
native!(
    Builder,
    "rssn/native/device-builder/v1",
    1_048_576,
    |b: &Builder| {
        if b.0.as_ref().is_some_and(|b| b.len() <= 4096) {
            Ok(())
        } else {
            Err(api::error(
                "builder consumed or instruction budget exceeded",
            ))
        }
    }
);
native!(
    Compiled,
    "rssn/native/device-kernel/v1",
    1_048_576,
    |b: &Compiled| {
        budget(b.native.ir())?;
        b.native.ir().validate().map_err(api::error)
    }
);
native!(Op, "rssn/native/device-op/v1", 128, |_: &Op| Ok(()));
native!(
    Options,
    "rssn/native/device-emitter-options/v1",
    4096,
    |o: &Options| {
        if o.0.entry.len() <= 128 {
            Ok(())
        } else {
            Err(api::error("entry identifier exceeds128bytes"))
        }
    }
);
fn budget(ir: &KernelIr) -> Result<(), OpError> {
    if ir.ops.len() > 4096 || ir.inputs > 64 || ir.params > 64 || ir.outputs.len() > 64 {
        return Err(api::error("device IR shape/work budget exceeded"));
    }
    Ok(())
}
fn work(
    ir: &KernelIr,
    rows: usize,
    lens: &[usize],
) -> Result<(), OpError> {
    budget(ir)?;
    let iterations = lens.iter().copied().max().unwrap_or(rows).max(1);
    let mut n = rows.checked_mul(ir.ops.len() + ir.outputs.len());
    if ir.fold.is_some() {
        n = n.and_then(|n| n.checked_mul(iterations));
    }
    if n.is_none_or(|n| n > 50_000_000)
        || rows.checked_mul(ir.outputs.len()).is_none_or(|n| n > LIMIT)
    {
        return Err(api::error("device launch work/allocation budget exceeded"));
    }
    Ok(())
}
fn check_columns(
    ir: &KernelIr,
    rows: usize,
    columns: &[&[f64]],
    params: &[f64],
) -> Result<(), OpError> {
    let lens = columns.iter().map(|c| c.len()).collect::<Vec<_>>();
    if lens.len() != ir.inputs || params.len() != ir.params {
        return Err(api::error("device argument arity mismatch"));
    }
    let data = columns
        .iter()
        .flat_map(|c| c.iter().copied())
        .collect::<Vec<_>>();
    ir.check_execution_work(
        rows,
        Columns {
            data: &data,
            offsets: &k::column_offsets(&lens),
            lens: &lens,
        },
        params,
        50_000_000,
    )
    .map_err(api::error)
}
fn check_buffer(
    device: &dyn Device,
    ir: &KernelIr,
    rows: usize,
    input: &d::Buffer,
    lens: &[usize],
    params: &[f64],
) -> Result<(), OpError> {
    if lens.len() != ir.inputs
        || params.len() != ir.params
        || lens.iter().try_fold(0usize, |n, l| n.checked_add(*l)) != Some(input.len())
    {
        return Err(api::error("device buffer/columns/parameter shape mismatch"));
    }
    let mut data = vec![0.0; input.len()];
    device.download(input, &mut data).map_err(api::error)?;
    ir.check_execution_work(
        rows,
        Columns {
            data: &data,
            offsets: &k::column_offsets(lens),
            lens,
        },
        params,
        50_000_000,
    )
    .map_err(api::error)
}
fn dtype(n: &str) -> Result<Dtype, OpError> {
    match n {
        | "F32" => Ok(Dtype::F32),
        | "F64" => Ok(Dtype::F64),
        | _ => Err(api::error("dtype must be F32/F64")),
    }
}
fn reduce(n: &str) -> Result<ReduceOp, OpError> {
    match n {
        | "Sum" => Ok(ReduceOp::Sum),
        | "Prod" => Ok(ReduceOp::Prod),
        | "Min" => Ok(ReduceOp::Min),
        | "Max" => Ok(ReduceOp::Max),
        | _ => Err(api::error("reduction must be Sum/Prod/Min/Max")),
    }
}
fn target(n: &str) -> Result<Target, OpError> {
    match n {
        | "Wgsl" => Ok(Target::Wgsl),
        | "Glsl" => Ok(Target::Glsl),
        | "OpenCl" => Ok(Target::OpenCl),
        | "Cuda" => Ok(Target::Cuda),
        | "Metal" => Ok(Target::Metal),
        | _ => Err(api::error("unknown source target")),
    }
}
fn cmp(n: &str) -> Result<Cmp, OpError> {
    match n {
        | "Lt" => Ok(Cmp::Lt),
        | "Le" => Ok(Cmp::Le),
        | "Gt" => Ok(Cmp::Gt),
        | "Ge" => Ok(Cmp::Ge),
        | "Eq" => Ok(Cmp::Eq),
        | "Ne" => Ok(Cmp::Ne),
        | _ => Err(api::error("unknown comparison")),
    }
}
fn ir(
    v: KernelIr,
    frame: Option<TapeSnapshot>,
) -> NewObject<KernelObject> {
    let dims = frame.as_ref().map_or_default(|f| {
        f.input_dimensions
            .iter()
            .chain(&f.parameter_dimensions)
            .chain(&f.output_dimensions)
            .cloned()
            .collect()
    });
    NewObject::new(KernelObject { ir: v, frame }, dims)
}
fn list(v: &Value) -> Result<&[Value], OpError> {
    let Value::Tensor(t) = v else {
        return Err(api::error("expected vector"));
    };
    if t.shape != [t.data.len()] || t.data.len() > LIMIT {
        return Err(api::error("invalid vector extent"));
    }
    Ok(&t.data)
}
fn physical(
    v: &Value,
    scales: &[f64],
    dimensions: &[DimExpr],
    columns: bool,
) -> Result<Vec<Vec<f64>>, OpError> {
    let v = list(v)?;
    if v.len() != scales.len() {
        return Err(api::error("physical columns/parameters arity mismatch"));
    }
    let mut cells = 0;
    v.iter()
        .zip(scales)
        .zip(dimensions)
        .map(|((v, s), dim)| {
            let inherited = if columns {
                let Value::Tensor(t) = v else {
                    return Err(api::error("physical column must be a vector"));
                };
                if matches!(t.elem, TypeExpr::Dyn(_)) {
                    None
                } else {
                    Some(api::dimension(v)?)
                }
            } else {
                None
            };
            let values = if columns {
                list(v)?
            } else {
                std::slice::from_ref(v)
            };
            cells += values.len();
            if cells > LIMIT {
                return Err(api::error("physical data allocation budget exceeded"));
            }
            if let Some(d) = &inherited {
                api::same_dimension(d, dim)?;
            }
            values
                .iter()
                .map(|v| {
                    let actual = if matches!(v, Value::Dyn(_) | Value::Quantity(_)) {
                        api::dimension(v)?
                    } else {
                        inherited.clone().unwrap_or(api::dimension(v)?)
                    };
                    api::same_dimension(&actual, dim)?;
                    api::scale_product(crate::rules::numerical_abi::numeric(v, 0)?, &[], &[*s])
                })
                .collect()
        })
        .collect()
}
fn physical_output(
    data: Vec<Vec<f64>>,
    f: &TapeSnapshot,
) -> Result<Dynamic, OpError> {
    if data.len() != f.output_scales.len() {
        return Err(api::error("output frame mismatch"));
    }
    let data = data
        .into_iter()
        .enumerate()
        .map(|(i, col)| {
            let dim = f.output_dimensions[i].clone();
            Ok(Value::Tensor(TensorValue {
                shape: vec![col.len()],
                elem: TypeExpr::Float.with_dim(dim.clone()),
                data: col
                    .into_iter()
                    .map(|x| {
                        api::quantity(
                            api::scale_product(x, &[f.output_scales[i]], &[])?,
                            dim.clone(),
                        )
                    })
                    .collect::<Result<_, _>>()?,
            }))
        })
        .collect::<Result<Vec<_>, OpError>>()?;
    Ok(Dynamic(Value::Tensor(TensorValue {
        shape: vec![data.len()],
        elem: api::dynamic_type(),
        data,
    })))
}
fn parse_op(
    kind: &str,
    args: &[usize],
    c: f64,
) -> Result<KernelOp, OpError> {
    let a = |i: usize| {
        args.get(i)
            .copied()
            .and_then(|v| u32::try_from(v).ok())
            .ok_or_else(|| api::error("missing/overflowed register"))
    };
    let index = |i: usize| a(i);
    Ok(match kind {
        | "Const" => KernelOp::Const(c),
        | "Input" => KernelOp::Input(index(0)?),
        | "Param" => KernelOp::Param(index(0)?),
        | "Add" => KernelOp::Add(a(0)?, a(1)?),
        | "Sub" => KernelOp::Sub(a(0)?, a(1)?),
        | "Mul" => KernelOp::Mul(a(0)?, a(1)?),
        | "Div" => KernelOp::Div(a(0)?, a(1)?),
        | "Neg" => KernelOp::Neg(a(0)?),
        | "PowF" => KernelOp::PowF(a(0)?, a(1)?),
        | "PowI" => {
            if c.fract() != 0.0 || c < f64::from(i32::MIN) || c > f64::from(i32::MAX) {
                return Err(api::error("integer exponent out of range"));
            }
            KernelOp::PowI(a(0)?, c as i32)
        },
        | "Atan2" => KernelOp::Atan2(a(0)?, a(1)?),
        | "Select" => KernelOp::Select(a(0)?, a(1)?, a(2)?),
        | "Gather" => KernelOp::Gather(index(0)?, a(1)?),
        | "Row" => KernelOp::Row,
        | "Rows" => KernelOp::Rows,
        | "Len" => KernelOp::Len(index(0)?),
        | "Inner" => KernelOp::Inner,
        | "IConst" => {
            if c.fract() != 0.0 || c < f64::from(i32::MIN) || c > f64::from(i32::MAX) {
                return Err(api::error("index constant must fit32bits"));
            }
            KernelOp::IConst(c as i64)
        },
        | "IAdd" => KernelOp::IAdd(a(0)?, a(1)?),
        | "ISub" => KernelOp::ISub(a(0)?, a(1)?),
        | "IMul" => KernelOp::IMul(a(0)?, a(1)?),
        | "IDiv" => KernelOp::IDiv(a(0)?, a(1)?),
        | "IMod" => KernelOp::IMod(a(0)?, a(1)?),
        | "IMin" => KernelOp::IMin(a(0)?, a(1)?),
        | "IMax" => KernelOp::IMax(a(0)?, a(1)?),
        | "ToNum" => KernelOp::ToNum(a(0)?),
        | "ToIdx" => KernelOp::ToIdx(a(0)?),
        | "Load" => KernelOp::Load(index(0)?, a(1)?),
        | "Acc" => KernelOp::Acc(index(0)?),
        | _ => {
            return Err(api::error(
                "unsupported opcode; use unary/comparison constructors for typed variants",
            ));
        },
    })
}
object_binding!(cpu,"num_native_device_cpu","crate::backend::device::cpu::CpuDevice::sequential","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,()->NewObject<DeviceObject>,|ctx,|{let _=ctx;Ok(NewObject::normalized(DeviceObject(Arc::new(d::cpu::CpuDevice::sequential()))))});
object_binding!(cpu_parallel,"num_native_device_cpu_parallel","crate::backend::device::cpu::CpuDevice::parallel","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,()->NewObject<DeviceObject>,|ctx,|{let _=ctx;Ok(NewObject::normalized(DeviceObject(Arc::new(d::cpu::CpuDevice::parallel()))))});
object_binding!(next_uid,"num_native_device_next_uid","crate::backend::device::next_uid","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,()->u64,|ctx,|{let _=ctx;Ok(d::next_uid())});
object_binding!(build_devices,"num_native_device_build_devices","crate::backend::device::build_devices","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,()->Many<String>,|ctx,|{let _=ctx;Ok(Many(d::build_devices()))});
object_binding!(uid,"num_native_device_uid","crate::backend::device::Device::uid","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<DeviceObject>)->u64,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.0.uid()))});
object_binding!(info,"num_native_device_info","crate::backend::device::Device::info","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<DeviceObject>)->String,|ctx,a|{let _=ctx;a.with(ctx,|a|serde_json::to_string(a.0.info()).map_err(api::error))});
object_binding!(kind,"num_native_device_kind","crate::backend::device::DeviceKind::tag","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<DeviceObject>)->String,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.0.info().kind.tag()))});
object_binding!(preferred_dtype,"num_native_device_preferred_dtype","crate::backend::device::preferred_dtype","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<DeviceObject>)->String,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(format!("{:?}",d::preferred_dtype(a.0.info()))))});
object_binding!(dtype_epsilon,"num_native_device_dtype_epsilon","crate::backend::device::dtype_epsilon","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<DeviceObject>,kind:String)->f64,|ctx,a,kind|{let _=ctx;a.with(ctx,|a|Ok(d::dtype_epsilon(a.0.info(),dtype(&kind)?)))});
object_binding!(synchronize,"num_native_device_synchronize","crate::backend::device::Device::synchronize","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<DeviceObject>)->bool,|ctx,a|{let _=ctx;a.with(ctx,|a|{a.0.synchronize().map_err(api::error)?;Ok(true)})});
object_binding!(registry_new,"num_native_device_registry_new","crate::backend::device::DeviceRegistry::new","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,()->NewObject<Registry>,|ctx,|{let _=ctx;Ok(NewObject::normalized(Registry(d::DeviceRegistry::new())))});
object_binding!(registry_empty,"num_native_device_registry_empty","crate::backend::device::DeviceRegistry::empty","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,()->NewObject<Registry>,|ctx,|{let _=ctx;Ok(NewObject::normalized(Registry(d::DeviceRegistry::empty())))});
object_binding!(registry_global,"num_native_device_registry_global","crate::backend::device::registry","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,()->String,|ctx,|{let _=ctx;serde_json::to_string(&d::registry().list()).map_err(api::error)});
object_binding!(registry_list,"num_native_device_registry_list","crate::backend::device::DeviceRegistry::list","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<Registry>)->String,|ctx,a|{let _=ctx;a.with(ctx,|a|serde_json::to_string(&a.0.list()).map_err(api::error))});
object_binding!(registry_register,"num_native_device_registry_register","crate::backend::device::DeviceRegistry::register","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<Registry>,device:Handle<DeviceObject>)->bool,|ctx,a,device|{let _=ctx;device.with(ctx,|d|a.with(ctx,|a|{if a.0.list().len()>=2{return Err(api::error("owned registry device count budget exceeded"));}a.0.register(d.0.clone()).map_err(api::error)?;Ok(true)}))});
object_binding!(registry_get,"num_native_device_registry_get","crate::backend::device::DeviceRegistry::get","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<Registry>,id:String)->Maybe<NewObject<DeviceObject>>,|ctx,a,id|{let _=ctx;a.with(ctx,|a|Ok(Maybe(a.0.get(&id).map(|d|NewObject::normalized(DeviceObject(d))))))});
object_binding!(registry_default,"num_native_device_registry_default","crate::backend::device::DeviceRegistry::default_device","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<Registry>)->Maybe<NewObject<DeviceObject>>,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(Maybe(a.0.default_device().map(|d|NewObject::normalized(DeviceObject(d))))))});
object_binding!(registry_first,"num_native_device_registry_first","crate::backend::device::DeviceRegistry::first_of","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<Registry>,kind:String)->Maybe<NewObject<DeviceObject>>,|ctx,a,kind|{let _=ctx;let kind=match kind.as_str(){"Cpu"=>d::DeviceKind::Cpu,"CpuParallel"=>d::DeviceKind::CpuParallel,_=>return Err(api::error("explicit CPU registry specialization expects Cpu/CpuParallel"))};a.with(ctx,|a|Ok(Maybe(a.0.first_of(&kind).map(|d|NewObject::normalized(DeviceObject(d))))))});
object_binding!(alloc,"num_native_device_alloc","crate::backend::device::Device::alloc","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<DeviceObject>,len:usize,kind:String)->NewObject<BufferObject>,|ctx,a,len,kind|{let _=ctx;if len>LIMIT{return Err(api::error("buffer allocation budget exceeded"));}a.with(ctx,|a|Ok(NewObject::normalized(BufferObject(a.0.alloc(len,dtype(&kind)?).map_err(api::error)?))))});
object_binding!(upload,"num_native_device_upload","crate::backend::device::Device::upload","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<DeviceObject>,values:Vec<f64>,kind:String)->NewObject<BufferObject>,|ctx,a,values,kind|{let _=ctx;if values.len()>LIMIT{return Err(api::error("buffer allocation budget exceeded"));}a.with(ctx,|a|Ok(NewObject::normalized(BufferObject(a.0.upload(&values,dtype(&kind)?).map_err(api::error)?))))});
object_binding!(download,"num_native_device_download","crate::backend::device::Device::download","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<DeviceObject>,buffer:Handle<BufferObject>)->Vec<f64>,|ctx,a,buffer|{let _=ctx;a.with(ctx,|a|buffer.with(ctx,|b|{let mut out=vec![0.0;b.0.len()];a.0.download(&b.0,&mut out).map_err(api::error)?;Ok(out)}))});
object_binding!(buffer_len,"num_native_device_buffer_len","crate::backend::device::Buffer::len","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(b:Handle<BufferObject>)->usize,|ctx,b|{let _=ctx;b.with(ctx,|b|Ok(b.0.len()))});
object_binding!(buffer_is_empty,"num_native_device_buffer_is_empty","crate::backend::device::Buffer::is_empty","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(b:Handle<BufferObject>)->bool,|ctx,b|{let _=ctx;b.with(ctx,|b|Ok(b.0.is_empty()))});
object_binding!(buffer_owner,"num_native_device_buffer_owner","crate::backend::device::Buffer::owner","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(b:Handle<BufferObject>)->u64,|ctx,b|{let _=ctx;b.with(ctx,|b|Ok(b.0.owner()))});
object_binding!(buffer_dtype,"num_native_device_buffer_dtype","crate::backend::device::Buffer::dtype","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(b:Handle<BufferObject>)->String,|ctx,b|{let _=ctx;b.with(ctx,|b|Ok(format!("{:?}",b.0.dtype())))});
object_binding!(copy,"num_native_device_copy","crate::backend::device::Device::copy","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<DeviceObject>,src:Handle<BufferObject>,src_at:usize,dst:Handle<BufferObject>,dst_at:usize,len:usize)->bool,|ctx,a,src,src_at,dst,dst_at,len|{let _=ctx;a.with(ctx,|a|src.with(ctx,|s|dst.with_mut(ctx,|t|{a.0.copy(&s.0,src_at,&mut t.0,dst_at,len).map_err(api::error)?;Ok(true)})))});
object_binding!(check_copy,"num_native_device_check_copy","crate::backend::device::check_copy","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<DeviceObject>,src:Handle<BufferObject>,src_at:usize,dst:Handle<BufferObject>,dst_at:usize,len:usize)->bool,|ctx,a,src,src_at,dst,dst_at,len|{let _=ctx;a.with(ctx,|a|src.with(ctx,|s|dst.with(ctx,|t|{d::check_copy(a.0.uid(),&s.0,src_at,&t.0,dst_at,len).map_err(api::error)?;Ok(true)})))});
object_binding!(reduce,"num_native_device_reduce","crate::backend::device::Device::reduce","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<DeviceObject>,op:String,b:Handle<BufferObject>)->f64,|ctx,a,op,b|{let _=ctx;a.with(ctx,|a|b.with(ctx,|b|a.0.reduce(reduce(&op)?,&b.0).map_err(api::error)))});
object_binding!(array_reduce,"num_native_device_array_reduce","crate::backend::device::array::reduce","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<DeviceObject>,op:String,b:Handle<BufferObject>)->f64,|ctx,a,op,b|{let _=ctx;a.with(ctx,|a|b.with(ctx,|b|d::array::reduce(&*a.0,reduce(&op)?,&b.0).map_err(api::error)))});
object_binding!(matmul,"num_native_device_matmul","crate::backend::device::Device::matmul","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<DeviceObject>,b:Handle<BufferObject>,c:Handle<BufferObject>,out:Handle<BufferObject>,m:usize,k:usize,n:usize)->bool,|ctx,a,b,c,out,m,k,n|{let _=ctx;if m.checked_mul(k).and_then(|v|v.checked_mul(n)).is_none_or(|v|v>50_000_000){return Err(api::error("matrix product work budget exceeded"));}a.with(ctx,|a|b.with(ctx,|b|c.with(ctx,|c|out.with_mut(ctx,|out|{a.0.matmul(&b.0,&c.0,&mut out.0,m,k,n).map_err(api::error)?;Ok(true)}))))});
object_binding!(array_matmul,"num_native_device_array_matmul","crate::backend::device::array::matmul","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<DeviceObject>,b:Handle<BufferObject>,c:Handle<BufferObject>,out:Handle<BufferObject>,m:usize,k:usize,n:usize)->bool,|ctx,a,b,c,out,m,k,n|{let _=ctx;if m.checked_mul(k).and_then(|v|v.checked_mul(n)).is_none_or(|v|v>50_000_000){return Err(api::error("matrix product work budget exceeded"));}a.with(ctx,|a|b.with(ctx,|b|c.with(ctx,|c|out.with_mut(ctx,|out|{d::array::matmul(&*a.0,&b.0,&c.0,&mut out.0,m,k,n).map_err(api::error)?;Ok(true)}))))});
object_binding!(kernel_from_tape,"num_native_device_kernel_from_tape","crate::backend::device::kernel::KernelIr::from_tape","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",false,(tape:Input)->NewObject<KernelObject>,|ctx,tape|{let _=ctx;let f=execution::tape_snapshot(ctx,&tape.0)?;let v=KernelIr::from_tape(&f.tape).map_err(api::error)?;Ok(ir(v,Some(f)))});
object_binding!(kernel_validate,"num_native_device_kernel_validate","crate::backend::device::kernel::KernelIr::validate","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<KernelObject>)->bool,|ctx,a|{let _=ctx;a.with(ctx,|a|{a.ir.validate().map_err(api::error)?;Ok(true)})});
object_binding!(kernel_is_elementwise,"num_native_device_kernel_is_elementwise","crate::backend::device::kernel::KernelIr::is_elementwise","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<KernelObject>)->bool,|ctx,a|{let _=ctx;a.with(ctx,|a|{Ok(a.ir.is_elementwise())})});
object_binding!(kernel_structural_hash,"num_native_device_kernel_structural_hash","crate::backend::device::kernel::KernelIr::structural_hash","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<KernelObject>)->u64,|ctx,a|{let _=ctx;a.with(ctx,|a|{Ok(a.ir.structural_hash())})});
object_binding!(kernel_analyse,"num_native_device_kernel_analyse","crate::backend::device::kernel::KernelIr::analyse","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<KernelObject>)->String,|ctx,a|{let _=ctx;a.with(ctx,|a|{let a=a.ir.analyse().map_err(api::error)?;Ok(format!("{a:?}"))})});
object_binding!(kernel_eval,"num_native_device_kernel_eval","crate::backend::device::kernel::KernelIr::eval_row","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<KernelObject>,x:Vec<f64>,p:Vec<f64>)->Vec<f64>,|ctx,a,x,p|{let _=ctx;a.with(ctx,|a|{if !a.ir.is_elementwise()||x.len()!=a.ir.inputs||p.len()!=a.ir.params{return Err(api::error("row evaluation expects elementwise IR and matching arguments"));}work(&a.ir,1,&[1])?;Ok(a.ir.eval_row(&x,&p))})});
object_binding!(map_columns,"num_native_device_map_columns","crate::backend::device::map_columns","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(device:Handle<DeviceObject>,a:Handle<KernelObject>,columns:Vec<Vec<f64>>,params:Vec<f64>)->Vec<Vec<f64>>,|ctx,device,a,columns,params|{let _=ctx;device.with(ctx,|d|a.with(ctx,|a|{let rows=columns.first().map_or(0,Vec::len);work(&a.ir,rows,&columns.iter().map(Vec::len).collect::<Vec<_>>())?;let cols=columns.iter().map(Vec::as_slice).collect::<Vec<_>>();check_columns(&a.ir,rows,&cols,&params)?;d::map_columns(&*d.0,&a.ir,&cols,&params).map_err(api::error)}))});
object_binding!(run_cols,"num_native_device_run_cols","crate::backend::device::array::run_cols","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(device:Handle<DeviceObject>,a:Handle<KernelObject>,columns:Vec<Vec<f64>>,params:Vec<f64>)->Vec<Vec<f64>>,|ctx,device,a,columns,params|{let _=ctx;device.with(ctx,|d|a.with(ctx,|a|{let rows=columns.first().map_or(0,Vec::len);work(&a.ir,rows,&columns.iter().map(Vec::len).collect::<Vec<_>>())?;let cols=columns.iter().map(Vec::as_slice).collect::<Vec<_>>();check_columns(&a.ir,rows,&cols,&params)?;d::array::run_cols(&*d.0,&a.ir,rows,&cols,&params).map_err(api::error)}))});
object_binding!(map_physical,"num_native_device_map_physical","","Execute the actual native device kernel with heterogeneous measured input columns, parameters and output units from a retained tape reference frame.",false,(device:Handle<DeviceObject>,a:Handle<KernelObject>,columns:Input,params:Input)->Dynamic,|ctx,device,a,columns,params|{let _=ctx;device.with(ctx,|d|a.with(ctx,|a|{let f=a.frame.as_ref().ok_or_else(||api::error("physical execution requires an IR imported from a framed tape"))?;let columns=physical(&columns.0,&f.input_scales,&f.input_dimensions,true)?;let params=physical(&params.0,&f.parameter_scales,&f.parameter_dimensions,false)?.into_iter().flatten().collect::<Vec<_>>();let rows=columns.first().map_or(0,Vec::len);work(&a.ir,rows,&columns.iter().map(Vec::len).collect::<Vec<_>>())?;let cols=columns.iter().map(Vec::as_slice).collect::<Vec<_>>();check_columns(&a.ir,rows,&cols,&params)?;physical_output(d::map_columns(&*d.0,&a.ir,&cols,&params).map_err(api::error)?,f)}))});
object_binding!(compile,"num_native_device_compile","crate::backend::device::Device::compile","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(device:Handle<DeviceObject>,a:Handle<KernelObject>)->NewObject<Compiled>,|ctx,device,a|{let _=ctx;device.with(ctx,|d|a.with(ctx,|a|Ok(NewObject::normalized(Compiled{native:d.0.compile(&a.ir).map_err(api::error)?,device:d.0.clone(),frame:a.frame.clone()}))))});
object_binding!(compiled_ir,"num_native_device_compiled_ir","crate::backend::device::CompiledKernel::ir","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<Compiled>)->NewObject<KernelObject>,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(ir(a.native.ir().clone(),a.frame.clone())))});
object_binding!(compiled_as_any,"num_native_device_compiled_as_any","crate::backend::device::CompiledKernel::as_any","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<Compiled>)->String,|ctx,a|{let _=ctx;a.with(ctx,|a|{let _=a.native.as_any();Ok(a.device.info().id.clone())})});
object_binding!(launch,"num_native_device_launch","crate::backend::device::Device::launch","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<Compiled>,rows:usize,input:Handle<BufferObject>,params:Vec<f64>,output:Handle<BufferObject>)->bool,|ctx,a,rows,input,params,output|{let _=ctx;a.with(ctx,|a|input.with(ctx,|input|output.with_mut(ctx,|output|{work(a.native.ir(),rows,&[rows])?;check_buffer(&*a.device,a.native.ir(),rows,&input.0,&vec![rows;a.native.ir().inputs],&params)?;a.device.launch(&*a.native,rows,&input.0,&params,&mut output.0).map_err(api::error)?;Ok(true)})))});
object_binding!(launch_cols,"num_native_device_launch_cols","crate::backend::device::Device::launch_cols","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<Compiled>,rows:usize,input:Handle<BufferObject>,lens:Vec<usize>,params:Vec<f64>,output:Handle<BufferObject>)->bool,|ctx,a,rows,input,lens,params,output|{let _=ctx;a.with(ctx,|a|input.with(ctx,|input|output.with_mut(ctx,|output|{work(a.native.ir(),rows,&lens)?;check_buffer(&*a.device,a.native.ir(),rows,&input.0,&lens,&params)?;a.device.launch_cols(&*a.native,rows,&input.0,&lens,&params,&mut output.0).map_err(api::error)?;Ok(true)})))});
object_binding!(check_launch,"num_native_device_check_launch","crate::backend::device::check_launch","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(device:Handle<DeviceObject>,a:Handle<KernelObject>,rows:usize,input:Handle<BufferObject>,output:Handle<BufferObject>)->bool,|ctx,device,a,rows,input,output|{let _=ctx;device.with(ctx,|d|a.with(ctx,|a|input.with(ctx,|input|output.with(ctx,|output|{d::check_launch(d.0.uid(),&a.ir,rows,&input.0,&output.0).map_err(api::error)?;Ok(true)}))))});
object_binding!(check_launch_cols,"num_native_device_check_launch_cols","crate::backend::device::check_launch_cols","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(device:Handle<DeviceObject>,a:Handle<KernelObject>,rows:usize,input:Handle<BufferObject>,output:Handle<BufferObject>,lens:Vec<usize>)->bool,|ctx,device,a,rows,input,output,lens|{let _=ctx;device.with(ctx,|d|a.with(ctx,|a|input.with(ctx,|input|output.with(ctx,|output|{d::check_launch_cols(d.0.uid(),&a.ir,rows,&input.0,&lens,&output.0).map_err(api::error)?;Ok(true)}))))});
object_binding!(pack_columns,"num_native_device_pack_columns","crate::backend::device::pack_columns","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(columns:Vec<Vec<f64>>)->Dynamic,|ctx,columns|{let _=ctx;let cols=columns.iter().map(Vec::as_slice).collect::<Vec<_>>();let(data,rows)=d::pack_columns(&cols).map_err(api::error)?;Ok(Dynamic(api::record(vec![("data",Output::encode(data)?),("rows",Value::int(rows))])))});
object_binding!(builder_new,"num_native_device_builder_new","crate::backend::device::kernel::KernelBuilder::new","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(inputs:usize,params:usize)->NewObject<Builder>,|ctx,inputs,params|{let _=ctx;if inputs>64||params>64{return Err(api::error("builder arity budget exceeded"));}Ok(NewObject::normalized(Builder(Some(KernelBuilder::new(inputs,params)))))});
object_binding!(builder_num,"num_native_device_builder_num","crate::backend::device::kernel::KernelBuilder::num","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<Builder>,value:f64)->u32,|ctx,a,value|{let _=ctx;a.with_mut(ctx,|a|{let b=a.0.as_mut().ok_or_else(||api::error("builder consumed"))?;if b.len()>=4096{return Err(api::error("builder instruction budget exceeded"));}Ok(b.num(value))})});
object_binding!(builder_idx,"num_native_device_builder_idx","crate::backend::device::kernel::KernelBuilder::idx","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<Builder>,value:i64)->u32,|ctx,a,value|{let _=ctx;a.with_mut(ctx,|a|{let b=a.0.as_mut().ok_or_else(||api::error("builder consumed"))?;if b.len()>=4096{return Err(api::error("builder instruction budget exceeded"));}Ok(b.idx(value))})});
object_binding!(builder_len,"num_native_device_builder_len","crate::backend::device::kernel::KernelBuilder::len","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<Builder>)->usize,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.0.as_ref().ok_or_else(||api::error("builder consumed"))?.len()))});
object_binding!(builder_is_empty,"num_native_device_builder_is_empty","crate::backend::device::kernel::KernelBuilder::is_empty","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<Builder>)->bool,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.0.as_ref().ok_or_else(||api::error("builder consumed"))?.is_empty()))});
object_binding!(op_new,"num_native_device_op_new","","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(kind:String,args:Vec<usize>,constant:f64)->NewObject<Op>,|ctx,kind,args,constant|{let _=ctx;if args.len()>3{return Err(api::error("opcode has at most3 operands"));}Ok(NewObject::normalized(Op(parse_op(&kind,&args,constant)?)))});
object_binding!(op_compare,"num_native_device_op_compare","","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(kind:String,a:u32,b:u32,index:bool)->NewObject<Op>,|ctx,kind,a,b,index|{let _=ctx;Ok(NewObject::normalized(Op(if index{KernelOp::ICmp(cmp(&kind)?,a,b)}else{KernelOp::Cmp(cmp(&kind)?,a,b)})))});
object_binding!(op_unary,"num_native_device_op_unary","","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(name:String,a:u32)->NewObject<Op>,|ctx,name,a|{let _=ctx;let unary=crate::rules::numerical_abi::native_backend_scalar::intrinsic(&name)?;Ok(NewObject::normalized(Op(KernelOp::Unary(unary,a))))});
object_binding!(builder_push,"num_native_device_builder_push","crate::backend::device::kernel::KernelBuilder::push","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<Builder>,op:Handle<Op>)->u32,|ctx,a,op|{let _=ctx;op.with(ctx,|op|a.with_mut(ctx,|a|{let b=a.0.as_mut().ok_or_else(||api::error("builder consumed"))?;if b.len()>=4096{return Err(api::error("builder instruction budget exceeded"));}Ok(b.push(op.0))}))});
object_binding!(builder_splice,"num_native_device_builder_splice","crate::backend::device::kernel::KernelBuilder::splice","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<Builder>,other:Handle<KernelObject>,inputs:Vec<u32>,params:Vec<u32>)->Vec<u32>,|ctx,a,other,inputs,params|{let _=ctx;other.with(ctx,|o|a.with_mut(ctx,|a|{let b=a.0.as_mut().ok_or_else(||api::error("builder consumed"))?;if b.len()+o.ir.ops.len()>4096||inputs.len()!=o.ir.inputs||params.len()!=o.ir.params||inputs.iter().chain(&params).any(|r|*r as usize>=b.len()){return Err(api::error("splice register/arity/instruction budget mismatch"));}Ok(b.splice(&o.ir,&inputs.into_iter().collect::<Vec<_>>(),&params.into_iter().collect::<Vec<_>>()).into_iter().collect())}))});
object_binding!(builder_finish,"num_native_device_builder_finish","crate::backend::device::kernel::KernelBuilder::finish","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<Builder>,outputs:Vec<u32>)->NewObject<KernelObject>,|ctx,a,outputs|{let _=ctx;a.with_mut(ctx,|a|{let b=a.0.take().ok_or_else(||api::error("builder consumed"))?;let out=b.finish(outputs.into_iter().collect()).map_err(api::error)?;Ok(ir(out,None))})});
object_binding!(builder_finish_fold,"num_native_device_builder_finish_fold","crate::backend::device::kernel::KernelBuilder::finish_fold","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<Builder>,outputs:Vec<u32>,count:u32,operations:Vec<String>,registers:Vec<u32>)->NewObject<KernelObject>,|ctx,a,outputs,count,operations,registers|{let _=ctx;if operations.len()!=registers.len()||registers.len()>64{return Err(api::error("fold accumulator arity/budget invalid"));}let fold=k::Fold{count,accs:operations.iter().zip(registers).map(|(op,r)|Ok((reduce(op)?,r))).collect::<Result<_,OpError>>()?};a.with_mut(ctx,|a|{let b=a.0.take().ok_or_else(||api::error("builder consumed"))?;let out=b.finish_fold(outputs.into_iter().collect(),Some(fold)).map_err(api::error)?;Ok(ir(out,None))})});
object_binding!(op_is_elementwise,"num_native_device_op_is_elementwise","crate::backend::device::kernel::KernelOp::is_elementwise","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<Op>)->bool,|ctx,a|{let _=ctx;a.with(ctx,|a|{Ok(a.0.is_elementwise())})});
object_binding!(op_ty,"num_native_device_op_ty","crate::backend::device::kernel::KernelOp::ty","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<Op>)->String,|ctx,a|{let _=ctx;a.with(ctx,|a|{Ok(format!("{:?}",a.0.ty()))})});
object_binding!(op_operand_ty,"num_native_device_op_operand_ty","crate::backend::device::kernel::KernelOp::operand_ty","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<Op>)->String,|ctx,a|{let _=ctx;a.with(ctx,|a|{Ok(format!("{:?}",a.0.operand_ty()))})});
object_binding!(op_operands,"num_native_device_op_operands","crate::backend::device::kernel::KernelOp::operands","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<Op>)->Vec<u32>,|ctx,a|{let _=ctx;a.with(ctx,|a|{let(r,n)=a.0.operands();Ok(r[..n].to_vec())})});
object_binding!(op_map,"num_native_device_op_map","crate::backend::device::kernel::KernelOp::map_operands","Explicit finite register-map specialization of the native generic operand mapper.",true,(a:Handle<Op>,mapping:Vec<u32>)->NewObject<Op>,|ctx,a,mapping|{let _=ctx;a.with(ctx,|a|{let(r,n)=a.0.operands();if r[..n].iter().any(|r|*r as usize>=mapping.len()){return Err(api::error("operand mapping does not cover source registers"));}Ok(NewObject::normalized(Op(a.0.map_operands(|r|mapping[r as usize]))))})});
object_binding!(options_new,"num_native_device_options_new","crate::backend::device::emit::EmitOptions::new","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(language:String)->NewObject<Options>,|ctx,language|{let _=ctx;Ok(NewObject::normalized(Options(EmitOptions::new(target(&language)?))))});
object_binding!(options_precision,"num_native_device_options_precision","crate::backend::device::emit::EmitOptions::precision","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<Options>,kind:String)->NewObject<Options>,|ctx,a,kind|{let _=ctx;a.with(ctx,|a|{let p=match kind.as_str(){"F32"=>Precision::F32,"F64"=>Precision::F64,"DoubleSingle"=>Precision::DoubleSingle,_=>return Err(api::error("unknown emitter precision"))};Ok(NewObject::normalized(Options(a.0.clone().precision(p))))})});
object_binding!(target_name,"num_native_device_target_name","crate::backend::device::emit::Target::name","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(language:String)->String,|ctx,language|{let _=ctx;Ok(target(&language)?.name().into())});
object_binding!(target_has_f64,"num_native_device_target_has_f64","crate::backend::device::emit::Target::has_f64","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(language:String)->bool,|ctx,language|{let _=ctx;Ok(target(&language)?.has_f64())});
object_binding!(emit,"num_native_device_emit","crate::backend::device::emit::emit","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<KernelObject>,options:Handle<Options>)->String,|ctx,a,options|{let _=ctx;a.with(ctx,|a|options.with(ctx,|o|native_emit::emit(&a.ir,&o.0).map_err(api::error)))});
object_binding!(emit_ds,"num_native_device_emit_ds","crate::backend::device::ds::emit_wgsl_ds","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<KernelObject>,options:Handle<Options>)->String,|ctx,a,options|{let _=ctx;a.with(ctx,|a|options.with(ctx,|o|d::ds::emit_wgsl_ds(&a.ir,&o.0).map_err(api::error)))});
object_binding!(emitter_target,"num_native_device_emitter_target","crate::backend::device::wgsl::KernelEmitter::target","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(options:Handle<Options>)->String,|ctx,options|{let _=ctx;options.with(ctx,|o|Ok(wgsl::SourceEmitter{options:o.0.clone()}.target().into()))});
object_binding!(emitter_emit,"num_native_device_emitter_emit","crate::backend::device::wgsl::KernelEmitter::emit","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<KernelObject>,options:Handle<Options>)->String,|ctx,a,options|{let _=ctx;a.with(ctx,|a|options.with(ctx,|o|wgsl::SourceEmitter{options:o.0.clone()}.emit(&a.ir).map_err(api::error)))});
object_binding!(emit_wgsl,"num_native_device_emit_wgsl","crate::backend::device::wgsl::emit_wgsl","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<KernelObject>,workgroup:u32,entry:String)->String,|ctx,a,workgroup,entry|{let _=ctx;if workgroup==0||workgroup>1024||entry.len()>128||entry.is_empty()||!entry.bytes().all(|c|c.is_ascii_alphanumeric()||c==b'_'){return Err(api::error("invalid WGSL entry/workgroup"));}a.with(ctx,|a|wgsl::emit_wgsl(&a.ir,&wgsl::WgslOptions{workgroup_size:workgroup,entry}).map_err(api::error))});
object_binding!(column_offsets,"num_native_device_column_offsets","crate::backend::device::kernel::column_offsets","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(lens:Vec<usize>)->Vec<usize>,|ctx,lens|{let _=ctx;if lens.iter().try_fold(0usize,|n,l|n.checked_add(*l)).is_none_or(|n|n>LIMIT){return Err(api::error("column size budget exceeded"));}Ok(k::column_offsets(&lens))});
object_binding!(num_to_idx,"num_native_device_num_to_idx","crate::backend::device::kernel::num_to_idx","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(x:f64)->i64,|ctx,x|{let _=ctx;Ok(k::num_to_idx(x))});
object_binding!(select_is_true,"num_native_device_select_is_true","crate::backend::device::kernel::select_is_true","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(x:f64)->bool,|ctx,x|{let _=ctx;Ok(k::select_is_true(x))});
object_binding!(gather_slot,"num_native_device_gather_slot","crate::backend::device::kernel::gather_slot","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(base:u32,offset:f64,inputs:usize)->Maybe<usize>,|ctx,base,offset,inputs|{let _=ctx;Ok(Maybe(k::gather_slot(base,offset,inputs)))});
object_binding!(cmp,"num_native_device_cmp","crate::backend::device::kernel::cmp_value","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(kind:String,a:f64,b:f64)->f64,|ctx,kind,a,b|{let _=ctx;Ok(k::cmp_value(cmp(&kind)?,a,b))});
object_binding!(literal,"num_native_device_literal","crate::backend::device::wgsl::f32_literal","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(x:f64)->String,|ctx,x|{let _=ctx;Ok(wgsl::f32_literal(x))});
object_binding!(ds_join,"num_native_device_ds_join","crate::backend::device::ds::join","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:f64,b:f64)->f64,|ctx,a,b|{let _=ctx;if [a,b].iter().any(|x|x.abs()>f64::from(f32::MAX)){return Err(api::error("double-single limbs must fit f32"));}Ok(d::ds::join([a as f32,b as f32]))});
object_binding!(ds_split,"num_native_device_ds_split","crate::backend::device::ds::split","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(x:f64)->Vec<f64>,|ctx,x|{let _=ctx;Ok(d::ds::split(x).into_iter().map(f64::from).collect())});
object_binding!(reduce_identity,"num_native_device_reduce_identity","crate::backend::device::ReduceOp::identity","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(op:String)->Dynamic,|ctx,op|{let _=ctx;Ok(Dynamic(Value::Float(reduce(&op)?.identity())))});
object_binding!(reduce_combine,"num_native_device_reduce_combine","crate::backend::device::ReduceOp::combine","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(op:String,a:f64,b:f64)->f64,|ctx,op,a,b|{let _=ctx;Ok(reduce(&op)?.combine(a,b))});
object_binding!(rk4,"num_native_device_rk4","crate::backend::device::kernel::rk4_step","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<KernelObject>,states:usize)->NewObject<KernelObject>,|ctx,a,states|{let _=ctx;a.with(ctx,|a|{if states>32||a.ir.ops.len()>128{return Err(api::error("step IR expansion budget exceeded"));}Ok(ir(k::rk4_step(&a.ir,states).map_err(api::error)?,None))})});
object_binding!(dopri,"num_native_device_dopri","crate::backend::device::array::dopri5_step","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<KernelObject>,states:usize)->NewObject<KernelObject>,|ctx,a,states|{let _=ctx;a.with(ctx,|a|{if states>32||a.ir.ops.len()>128{return Err(api::error("step IR expansion budget exceeded"));}Ok(ir(d::array::dopri5_step(&a.ir,states).map_err(api::error)?,None))})});
object_binding!(rhs_shape,"num_native_device_rhs_shape","crate::backend::device::kernel::rhs_shape","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<KernelObject>,states:usize)->bool,|ctx,a,states|{let _=ctx;a.with(ctx,|a|k::rhs_shape(&a.ir,states).map_err(api::error))});
object_binding!(rk4_alias,"num_native_device_rk4_alias","crate::backend::device::rk4_step","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<KernelObject>,states:usize)->NewObject<KernelObject>,|ctx,a,states|{let _=ctx;a.with(ctx,|a|{if states>32||a.ir.ops.len()>128{return Err(api::error("step IR expansion budget exceeded"));}Ok(ir(d::rk4_step(&a.ir,states).map_err(api::error)?,None))})});
object_binding!(matmul_kernel,"num_native_device_matmul_kernel","crate::backend::device::array::matmul_kernel","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,()->NewObject<KernelObject>,|ctx,|{let _=ctx;Ok(ir(d::array::matmul_kernel().map_err(api::error)?,None))});
object_binding!(spmv_kernel,"num_native_device_spmv_kernel","crate::backend::device::array::spmv_csr_kernel","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,()->NewObject<KernelObject>,|ctx,|{let _=ctx;Ok(ir(d::array::spmv_csr_kernel().map_err(api::error)?,None))});
object_binding!(reduce_pass,"num_native_device_reduce_pass","crate::backend::device::array::reduce_pass","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(op:String,inputs:usize,column:usize)->NewObject<KernelObject>,|ctx,op,inputs,column|{let _=ctx;if inputs>64||column>=inputs{return Err(api::error("reduction column/arity invalid"));}Ok(ir(d::array::reduce_pass(reduce(&op)?,inputs,column).map_err(api::error)?,None))});
object_binding!(map_reduce_pass,"num_native_device_map_reduce_pass","crate::backend::device::array::map_reduce_pass","Explicit checked native device adapter. Machine IR and buffers use normalized coordinates; use physical map for independently measured columns. Session, shape and bounded work are validated before native execution.",true,(a:Handle<KernelObject>,op:String)->NewObject<KernelObject>,|ctx,a,op|{let _=ctx;a.with(ctx,|a|Ok(ir(d::array::map_reduce_pass(&a.ir,reduce(&op)?).map_err(api::error)?,None)))});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    extra::install(i)?;
    features::install(i)?;
    api::install(
        i,
        &[
            &cpu::BINDING,
            &cpu_parallel::BINDING,
            &next_uid::BINDING,
            &build_devices::BINDING,
            &uid::BINDING,
            &info::BINDING,
            &kind::BINDING,
            &preferred_dtype::BINDING,
            &dtype_epsilon::BINDING,
            &synchronize::BINDING,
            &registry_new::BINDING,
            &registry_empty::BINDING,
            &registry_global::BINDING,
            &registry_list::BINDING,
            &registry_register::BINDING,
            &registry_get::BINDING,
            &registry_default::BINDING,
            &registry_first::BINDING,
            &alloc::BINDING,
            &upload::BINDING,
            &download::BINDING,
            &buffer_len::BINDING,
            &buffer_is_empty::BINDING,
            &buffer_owner::BINDING,
            &buffer_dtype::BINDING,
            &copy::BINDING,
            &check_copy::BINDING,
            &reduce::BINDING,
            &array_reduce::BINDING,
            &matmul::BINDING,
            &array_matmul::BINDING,
            &kernel_from_tape::BINDING,
            &kernel_validate::BINDING,
            &kernel_is_elementwise::BINDING,
            &kernel_structural_hash::BINDING,
            &kernel_analyse::BINDING,
            &kernel_eval::BINDING,
            &map_columns::BINDING,
            &run_cols::BINDING,
            &map_physical::BINDING,
            &compile::BINDING,
            &compiled_ir::BINDING,
            &compiled_as_any::BINDING,
            &launch::BINDING,
            &launch_cols::BINDING,
            &check_launch::BINDING,
            &check_launch_cols::BINDING,
            &pack_columns::BINDING,
            &builder_new::BINDING,
            &builder_num::BINDING,
            &builder_idx::BINDING,
            &builder_len::BINDING,
            &builder_is_empty::BINDING,
            &op_new::BINDING,
            &op_compare::BINDING,
            &op_unary::BINDING,
            &builder_push::BINDING,
            &builder_splice::BINDING,
            &builder_finish::BINDING,
            &builder_finish_fold::BINDING,
            &op_is_elementwise::BINDING,
            &op_ty::BINDING,
            &op_operand_ty::BINDING,
            &op_operands::BINDING,
            &op_map::BINDING,
            &options_new::BINDING,
            &options_precision::BINDING,
            &target_name::BINDING,
            &target_has_f64::BINDING,
            &emit::BINDING,
            &emit_ds::BINDING,
            &emitter_target::BINDING,
            &emitter_emit::BINDING,
            &emit_wgsl::BINDING,
            &column_offsets::BINDING,
            &num_to_idx::BINDING,
            &select_is_true::BINDING,
            &gather_slot::BINDING,
            &cmp::BINDING,
            &literal::BINDING,
            &ds_join::BINDING,
            &ds_split::BINDING,
            &reduce_identity::BINDING,
            &reduce_combine::BINDING,
            &rk4::BINDING,
            &dopri::BINDING,
            &rhs_shape::BINDING,
            &rk4_alias::BINDING,
            &matmul_kernel::BINDING,
            &spmv_kernel::BINDING,
            &reduce_pass::BINDING,
            &map_reduce_pass::BINDING,
        ],
        "numerical_vector",
    )
}
