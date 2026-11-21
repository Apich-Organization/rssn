//! Further concrete native device schemas; no raw pointers cross the ABI.
use super::BufferObject;
use super::Columns;
use super::DenseStorage;
use super::DeviceObject;
use super::Dynamic;
use super::Handle;
use super::Input;
use super::Installer;
use super::KernelObject;
use super::LIMIT;
use super::Maybe;
use super::NativeType;
use super::NewObject;
use super::OpError;
use super::RuleError;
use super::Scratch;
use super::TensorValue;
use super::TypeExpr;
use super::Value;
use super::api;
use super::check_columns;
use super::d;
use super::dtype;
use super::execution;
use super::ir;
use super::k;
use super::list;
use super::object_binding;
use super::reduce;
use super::work;
use crate::backend::device::array;
use crate::backend::device::plan;
use crate::rules::numerical_abi::codecs::Arg;
use crate::rules::numerical_abi::codecs::Output;
struct Analysis {
    native: k::Analysis,
    hash: u64,
}
struct Plan {
    native: plan::Plan,
    externals: usize,
    work: usize,
    nodes: usize,
}
struct StorageObject(Vec<f64>);
native!(
    Analysis,
    "rssn/native/device-analysis/v1",
    262144,
    |a: &Analysis| {
        if a.native.ty.len() <= 4096 {
            Ok(())
        } else {
            Err(api::error("analysis budget exceeded"))
        }
    }
);
native!(Plan, "rssn/native/device-plan/v1", 4_194_304, |a: &Plan| {
    if a.nodes <= 64 && a.work <= 4096 {
        Ok(())
    } else {
        Err(api::error("plan size budget exceeded"))
    }
});
native!(
    StorageObject,
    "rssn/native/device-storage/v1",
    2_097_152,
    |a: &StorageObject| {
        if a.0.len() <= LIMIT {
            Ok(())
        } else {
            Err(api::error("storage budget exceeded"))
        }
    }
);
fn source(v: &Value) -> Result<plan::Source, OpError> {
    let Value::Record(r) = v else {
        return Err(api::error("source must be {external} or {node,output}"));
    };
    if r.fields.len() == 1 {
        Ok(plan::Source::External(<usize as Arg>::decode(
            api::field(v, "external")?,
            0,
        )?))
    } else if r.fields.len() == 2 {
        Ok(plan::Source::Output(
            plan::NodeId(<usize as Arg>::decode(api::field(v, "node")?, 0)?),
            <usize as Arg>::decode(api::field(v, "output")?, 0)?,
        ))
    } else {
        Err(api::error("plan source fields invalid"))
    }
}
fn sources(v: &Value) -> Result<Vec<plan::Source>, OpError> {
    list(v)?.iter().map(source).collect()
}
fn boundary(v: &str) -> Result<array::Boundary, OpError> {
    match v {
        | "Zero" => Ok(array::Boundary::Zero),
        | "Periodic" => Ok(array::Boundary::Periodic),
        | "Clamp" => Ok(array::Boundary::Clamp),
        | _ => Err(api::error("boundary must be Zero/Periodic/Clamp")),
    }
}
object_binding!(analyse,"num_native_device_analyse_extra","crate::backend::device::kernel::KernelIr::analyse","Explicit native array/device object specialization. IR, buffer and generic storage values are normalized; physical programs retain their tape reference frame. Shapes and native execution work are bounded before launch.",true,(a:Handle<KernelObject>)->NewObject<Analysis>,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(NewObject::normalized(Analysis{native:a.ir.analyse().map_err(api::error)?,hash:a.ir.structural_hash()})))});
object_binding!(eval_at,"num_native_device_eval_at_extra","crate::backend::device::kernel::KernelIr::eval_at","Explicit native array/device object specialization. IR, buffer and generic storage values are normalized; physical programs retain their tape reference frame. Shapes and native execution work are bounded before launch.",true,(a:Handle<KernelObject>,analysis:Handle<Analysis>,row:usize,rows:usize,columns:Vec<Vec<f64>>,params:Vec<f64>)->Dynamic,|ctx,a,analysis,row,rows,columns,params|{let _=ctx;a.with(ctx,|a|analysis.with(ctx,|an|{if a.ir.structural_hash()!=an.hash||row>=rows{return Err(api::error("analysis belongs to another kernel or row exceeds launch"));}let lens=columns.iter().map(Vec::len).collect::<Vec<_>>();let refs=columns.iter().map(Vec::as_slice).collect::<Vec<_>>();work(&a.ir,rows,&lens)?;check_columns(&a.ir,rows,&refs,&params)?;let flat=columns.into_iter().flatten().collect::<Vec<_>>();let mut result=vec![0.0;a.ir.outputs.len()];a.ir.eval_at(&an.native,row,rows,Columns{data:&flat,offsets:&k::column_offsets(&lens),lens:&lens},&params,&mut Scratch::default(),&mut result);Ok(Dynamic(Value::Tensor(TensorValue{shape:vec![result.len()],elem:TypeExpr::Float,data:result.into_iter().map(Value::Float).collect()})))}))});
object_binding!(column_get,"num_native_device_column_get_extra","crate::backend::device::kernel::Columns<'_>::get","Explicit native array/device object specialization. IR, buffer and generic storage values are normalized; physical programs retain their tape reference frame. Shapes and native execution work are bounded before launch.",true,(columns:Vec<Vec<f64>>,column:usize,index:i64)->Dynamic,|ctx,columns,column,index|{let _=ctx;let lens=columns.iter().map(Vec::len).collect::<Vec<_>>();let flat=columns.into_iter().flatten().collect::<Vec<_>>();let value=Columns{data:&flat,offsets:&k::column_offsets(&lens),lens:&lens}.get(column,index);Ok(Dynamic(Value::Float(value)))});
object_binding!(to_tape,"num_native_device_to_tape_extra","crate::backend::device::kernel::KernelIr::to_tape","Explicit native array/device object specialization. IR, buffer and generic storage values are normalized; physical programs retain their tape reference frame. Shapes and native execution work are bounded before launch.",true,(a:Handle<KernelObject>)->Maybe<Dynamic>,|ctx,a|{let _=ctx;a.with(ctx,|a|{let native=a.ir.to_tape();Ok(Maybe(native.map(|tape|{let mut f=a.frame.as_ref().map_or_else(||execution::normalized_snapshot(tape.clone()),Clone::clone);f.tape=tape;execution::tape_value(ctx,f).map(Dynamic)}).transpose()?))})});
object_binding!(plan_new,"num_native_device_plan_new_extra","crate::backend::device::plan::Plan::new","Explicit native array/device object specialization. IR, buffer and generic storage values are normalized; physical programs retain their tape reference frame. Shapes and native execution work are bounded before launch.",true,(externals:usize)->NewObject<Plan>,|ctx,externals|{let _=ctx;if externals>64{return Err(api::error("plan external arity budget exceeded"));}Ok(NewObject::normalized(Plan{native:plan::Plan::new(externals),externals,work:0,nodes:0}))});
object_binding!(plan_add,"num_native_device_plan_add_extra","crate::backend::device::plan::Plan::add","Explicit native array/device object specialization. IR, buffer and generic storage values are normalized; physical programs retain their tape reference frame. Shapes and native execution work are bounded before launch.",true,(a:Handle<Plan>,ir:Handle<KernelObject>,inputs:Input,params:Vec<f64>)->usize,|ctx,a,ir,inputs,params|{let _=ctx;ir.with(ctx,|ir|a.with_mut(ctx,|a|{if !ir.ir.is_elementwise()||a.nodes>=64||a.work+ir.ir.ops.len()>4096{return Err(api::error("plan supports bounded elementwise native fusion; structured kernels use run_cols"));}let src=sources(&inputs.0)?;let id=a.native.add(ir.ir.clone(),src,params).map_err(api::error)?;a.work+=ir.ir.ops.len();a.nodes+=1;Ok(id.0)}))});
object_binding!(plan_result,"num_native_device_plan_result_extra","crate::backend::device::plan::Plan::result","Explicit native array/device object specialization. IR, buffer and generic storage values are normalized; physical programs retain their tape reference frame. Shapes and native execution work are bounded before launch.",true,(a:Handle<Plan>,result:Input)->bool,|ctx,a,result|{let _=ctx;a.with_mut(ctx,|a|{a.native.result(source(&result.0)?).map_err(api::error)?;Ok(true)})});
object_binding!(plan_fuse,"num_native_device_plan_fuse_extra","crate::backend::device::plan::Plan::fuse","Explicit native array/device object specialization. IR, buffer and generic storage values are normalized; physical programs retain their tape reference frame. Shapes and native execution work are bounded before launch.",true,(a:Handle<Plan>)->usize,|ctx,a|{let _=ctx;a.with_mut(ctx,|a|Ok(a.native.fuse()))});
object_binding!(plan_launches,"num_native_device_plan_launches_extra","crate::backend::device::plan::Plan::launches","Explicit native array/device object specialization. IR, buffer and generic storage values are normalized; physical programs retain their tape reference frame. Shapes and native execution work are bounded before launch.",true,(a:Handle<Plan>)->usize,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.native.launches()))});
object_binding!(plan_run,"num_native_device_plan_run_extra","crate::backend::device::plan::Plan::run","Explicit native array/device object specialization. IR, buffer and generic storage values are normalized; physical programs retain their tape reference frame. Shapes and native execution work are bounded before launch.",true,(a:Handle<Plan>,device:Handle<DeviceObject>,rows:usize,externals:Vec<Vec<f64>>)->Vec<Vec<f64>>,|ctx,a,device,rows,externals|{let _=ctx;a.with(ctx,|a|device.with(ctx,|d|{if externals.len()!=a.externals||rows.checked_mul(a.work).is_none_or(|n|n>50_000_000)||rows.checked_mul(64).is_none_or(|n|n>LIMIT){return Err(api::error("plan data/work budget exceeded"));}let cols=externals.iter().map(Vec::as_slice).collect::<Vec<_>>();a.native.run(&*d.0,rows,&cols).map_err(api::error)}))});
object_binding!(map_reduce,"num_native_device_map_reduce_extra","crate::backend::device::array::map_reduce","Explicit native array/device object specialization. IR, buffer and generic storage values are normalized; physical programs retain their tape reference frame. Shapes and native execution work are bounded before launch.",true,(device:Handle<DeviceObject>,a:Handle<KernelObject>,columns:Vec<Vec<f64>>,params:Vec<f64>,op:String)->Vec<f64>,|ctx,device,a,columns,params,op|{let _=ctx;device.with(ctx,|d|a.with(ctx,|a|{let rows=columns.first().map_or(0,Vec::len);let cols=columns.iter().map(Vec::as_slice).collect::<Vec<_>>();work(&a.ir,rows,&columns.iter().map(Vec::len).collect::<Vec<_>>())?;check_columns(&a.ir,rows,&cols,&params)?;array::map_reduce(&*d.0,&a.ir,&cols,&params,reduce(&op)?).map_err(api::error)}))});
object_binding!(reduce_column,"num_native_device_reduce_column_extra","crate::backend::device::array::reduce_column","Explicit native array/device object specialization. IR, buffer and generic storage values are normalized; physical programs retain their tape reference frame. Shapes and native execution work are bounded before launch.",true,(device:Handle<DeviceObject>,op:String,buffer:Handle<BufferObject>,lens:Vec<usize>,column:usize)->f64,|ctx,device,op,buffer,lens,column|{let _=ctx;device.with(ctx,|d|buffer.with(ctx,|b|{if lens.iter().try_fold(0usize,|n,l|n.checked_add(*l))!=Some(b.0.len()){return Err(api::error("reduction column lengths mismatch"));}array::reduce_column(&*d.0,reduce(&op)?,&b.0,&lens,column).map_err(api::error)}))});
object_binding!(pairwise_kernel,"num_native_device_pairwise_kernel_extra","crate::backend::device::array::pairwise_kernel","Explicit native array/device object specialization. IR, buffer and generic storage values are normalized; physical programs retain their tape reference frame. Shapes and native execution work are bounded before launch.",true,(a:Handle<KernelObject>,columns:usize,op:String,skip_self:bool)->NewObject<KernelObject>,|ctx,a,columns,op,skip_self|{let _=ctx;a.with(ctx,|a|{if columns>32{return Err(api::error("pairwise arity budget exceeded"));}Ok(ir(array::pairwise_kernel(&a.ir,columns,reduce(&op)?,skip_self).map_err(api::error)?,None))})});
object_binding!(stencil1d,"num_native_device_stencil1d_extra","crate::backend::device::array::stencil1d_kernel","Explicit native array/device object specialization. IR, buffer and generic storage values are normalized; physical programs retain their tape reference frame. Shapes and native execution work are bounded before launch.",true,(offsets:Vec<i64>,weights:Vec<f64>,edges:String)->NewObject<KernelObject>,|ctx,offsets,weights,edges|{let _=ctx;if offsets.len()!=weights.len()||offsets.len()>128{return Err(api::error("stencil arity budget exceeded"));}let points=offsets.into_iter().zip(weights).collect::<Vec<_>>();Ok(ir(array::stencil1d_kernel(&points,boundary(&edges)?).map_err(api::error)?,None))});
object_binding!(stencil2d,"num_native_device_stencil2d_extra","crate::backend::device::array::stencil2d_kernel","Explicit native array/device object specialization. IR, buffer and generic storage values are normalized; physical programs retain their tape reference frame. Shapes and native execution work are bounded before launch.",true,(offsets:Vec<[i64;2]>,weights:Vec<f64>,edges:String)->NewObject<KernelObject>,|ctx,offsets,weights,edges|{let _=ctx;if offsets.len()!=weights.len()||offsets.len()>128{return Err(api::error("stencil arity budget exceeded"));}let points=offsets.into_iter().zip(weights).map(|(p,w)|(p[0],p[1],w)).collect::<Vec<_>>();Ok(ir(array::stencil2d_kernel(&points,boundary(&edges)?).map_err(api::error)?,None))});
object_binding!(ensemble_rk4,"num_native_device_ensemble_rk4_extra","crate::backend::device::ensemble_rk4","Explicit native array/device object specialization. IR, buffer and generic storage values are normalized; physical programs retain their tape reference frame. Shapes and native execution work are bounded before launch.",true,(device:Handle<DeviceObject>,a:Handle<KernelObject>,params:Vec<f64>,initial:Vec<Vec<f64>>,t0:f64,t1:f64,steps:usize)->Vec<Vec<f64>>,|ctx,device,a,params,initial,t0,t1,steps|{let _=ctx;device.with(ctx,|d|a.with(ctx,|a|{let rows=initial.first().map_or(0,Vec::len);if !a.ir.is_elementwise()||steps==0||steps>65536||rows.checked_mul(steps).and_then(|n|n.checked_mul(a.ir.ops.len()+1)).is_none_or(|n|n>5_000_000){return Err(api::error("ensemble RK4 work budget exceeded"));}let y=initial.iter().map(Vec::as_slice).collect::<Vec<_>>();d::ensemble_rk4(&*d.0,&a.ir,&params,&y,t0,t1,steps).map_err(api::error)}))});
object_binding!(ensemble_dopri5,"num_native_device_ensemble_dopri5_extra","crate::backend::device::array::ensemble_dopri5","Explicit native array/device object specialization. IR, buffer and generic storage values are normalized; physical programs retain their tape reference frame. Shapes and native execution work are bounded before launch.",true,(device:Handle<DeviceObject>,a:Handle<KernelObject>,params:Vec<f64>,initial:Vec<Vec<f64>>,t0:f64,t1:f64,rtol:f64,atol:f64,first_step:f64,max_launches:usize)->Dynamic,|ctx,device,a,params,initial,t0,t1,rtol,atol,first_step,max_launches|{let _=ctx;device.with(ctx,|d|a.with(ctx,|a|{let rows=initial.first().map_or(0,Vec::len);if !a.ir.is_elementwise()||max_launches==0||max_launches>65536||rtol<=0.0||atol<=0.0||first_step<=0.0||rows.checked_mul(max_launches).and_then(|n|n.checked_mul(a.ir.ops.len()+1)).is_none_or(|n|n>5_000_000){return Err(api::error("adaptive ensemble controls/work budget invalid"));}let y=initial.iter().map(Vec::as_slice).collect::<Vec<_>>();let result=array::ensemble_dopri5(&*d.0,&a.ir,&params,&y,t0,t1,&array::Tolerance{rtol,atol,first_step:Some(first_step),max_launches}).map_err(api::error)?;Ok(Dynamic(api::record(vec![("values",Output::encode(result.y)?),("launches",Value::int(result.launches))])))}))});
object_binding!(storage_new,"num_native_device_storage_new_extra","","Explicit native array/device object specialization. IR, buffer and generic storage values are normalized; physical programs retain their tape reference frame. Shapes and native execution work are bounded before launch.",true,(values:Vec<f64>)->NewObject<StorageObject>,|ctx,values|{let _=ctx;Ok(NewObject::normalized(StorageObject(values)))});
object_binding!(buffer_new,"num_native_device_buffer_new_extra","crate::backend::device::Buffer::new","Explicit native array/device object specialization. IR, buffer and generic storage values are normalized; physical programs retain their tape reference frame. Shapes and native execution work are bounded before launch.",true,(owner:u64,kind:String,storage:Handle<StorageObject>)->NewObject<BufferObject>,|ctx,owner,kind,storage|{let _=ctx;storage.with(ctx,|s|{let dtype=dtype(&kind)?;Ok(NewObject::normalized(BufferObject(d::Buffer::new(owner,dtype,s.0.len(),Box::new(DenseStorage(s.0.clone()))))))})});
object_binding!(buffer_storage,"num_native_device_buffer_storage_extra","crate::backend::device::Buffer::storage","Explicit native array/device object specialization. IR, buffer and generic storage values are normalized; physical programs retain their tape reference frame. Shapes and native execution work are bounded before launch.",true,(b:Handle<BufferObject>)->Maybe<Vec<f64>>,|ctx,b|{let _=ctx;b.with(ctx,|b|Ok(Maybe(b.0.storage().as_any().downcast_ref::<DenseStorage>().map(|v|v.0.clone()))))});
object_binding!(buffer_storage_mut,"num_native_device_buffer_storage_mut_extra","crate::backend::device::Buffer::storage_mut","Explicit native array/device object specialization. IR, buffer and generic storage values are normalized; physical programs retain their tape reference frame. Shapes and native execution work are bounded before launch.",true,(b:Handle<BufferObject>)->Maybe<Vec<f64>>,|ctx,b|{let _=ctx;b.with_mut(ctx,|b|Ok(Maybe(b.0.storage_mut().as_any_mut().downcast_mut::<DenseStorage>().map(|v|v.0.clone()))))});
object_binding!(storage_as_any,"num_native_device_storage_as_any_extra","crate::backend::device::Storage::as_any","Explicit native array/device object specialization. IR, buffer and generic storage values are normalized; physical programs retain their tape reference frame. Shapes and native execution work are bounded before launch.",true,(b:Handle<BufferObject>)->Maybe<Dynamic>,|ctx,b|{let _=ctx;b.with(ctx,|b|Ok(Maybe(b.0.storage().as_any().downcast_ref::<DenseStorage>().map(|v|Output::encode(v.0.clone()).map(Dynamic)).transpose()?)))});
object_binding!(storage_as_any_mut,"num_native_device_storage_as_any_mut_extra","crate::backend::device::Storage::as_any_mut","Explicit native array/device object specialization. IR, buffer and generic storage values are normalized; physical programs retain their tape reference frame. Shapes and native execution work are bounded before launch.",true,(b:Handle<BufferObject>,values:Vec<f64>)->bool,|ctx,b,values|{let _=ctx;b.with_mut(ctx,|b|{let len=b.0.len();let v=b.0.storage_mut().as_any_mut().downcast_mut::<DenseStorage>().ok_or_else(||api::error("storage is not the explicit DenseStorage specialization"))?;if values.len()!=len{return Err(api::error("replacement storage length mismatch"));}v.0.copy_from_slice(&values);Ok(true)})});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &analyse::BINDING,
            &eval_at::BINDING,
            &column_get::BINDING,
            &to_tape::BINDING,
            &plan_new::BINDING,
            &plan_add::BINDING,
            &plan_result::BINDING,
            &plan_fuse::BINDING,
            &plan_launches::BINDING,
            &plan_run::BINDING,
            &map_reduce::BINDING,
            &reduce_column::BINDING,
            &pairwise_kernel::BINDING,
            &stencil1d::BINDING,
            &stencil2d::BINDING,
            &ensemble_rk4::BINDING,
            &ensemble_dopri5::BINDING,
            &storage_new::BINDING,
            &buffer_new::BINDING,
            &buffer_storage::BINDING,
            &buffer_storage_mut::BINDING,
            &storage_as_any::BINDING,
            &storage_as_any_mut::BINDING,
        ],
        "numerical_vector",
    )
}
