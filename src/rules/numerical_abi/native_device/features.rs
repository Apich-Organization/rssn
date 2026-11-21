//! Explicit feature-gated native device specializations.
#![cfg_attr(
    not(feature = "jit"),
    allow(unused_variables, dead_code, unused_imports)
)]
use super::DeviceObject;
use super::Dynamic;
use super::Handle;
use super::Installer;
use super::KernelIr;
use super::KernelObject;
use super::Many;
use super::Maybe;
use super::NativeType;
use super::NewObject;
use super::OpError;
use super::RuleError;
use super::api;
use super::budget;
#[cfg(feature = "jit")]
use super::check_columns;
use super::d;
use super::object_binding;
#[cfg(feature = "jit")]
use super::work;
#[cfg(feature = "jit")]
use crate::kernel::Value;
#[cfg(feature = "jit")]
use crate::rules::numerical_abi::codecs::Output;
#[cfg(feature = "gpu-wgpu")]
use std::sync::Arc;
struct NativeKernel {
    #[cfg(feature = "jit")]
    native: d::native::NativeKernel,
    ir: KernelIr,
}
struct DeviceBackend(d::DeviceBackend);
impl NativeType for NativeKernel {
    const KIND: &'static str = "rssn/native/device-native-kernel/v1";
    const MAX_BYTES: usize = 4_194_304;

    fn validate(&self) -> Result<(), OpError> {
        budget(&self.ir)
    }
}
impl NativeType for DeviceBackend {
    const KIND: &'static str = "rssn/native/device-backend/v1";
    const MAX_BYTES: usize = 4096;

    fn validate(&self) -> Result<(), OpError> {
        Ok(())
    }
}
macro_rules! feature_jit{($($body:tt)*)=>{{#[cfg(feature="jit")]{$($body)*}#[cfg(not(feature="jit"))]{Err(OpError::Unsupported{backend:"rssn".into(),feature:"requires RSSN jit feature".into()})}}};}
macro_rules! feature_gpu{($($body:tt)*)=>{{#[cfg(feature="gpu-wgpu")]{$($body)*}#[cfg(not(feature="gpu-wgpu"))]{Err(OpError::Unsupported{backend:"rssn".into(),feature:"requires RSSN gpu-wgpu feature".into()})}}};}
object_binding!(backend,"num_native_device_backend","crate::backend::device::DeviceBackend::new","Explicit native backend specialization. Feature-disabled JIT/GPU builds return Unsupported; machine representations are normalized and source, row, fold-work and allocation budgets are checked.",true,(device:Handle<DeviceObject>)->NewObject<DeviceBackend>,|ctx,device|{device.with(ctx,|d|Ok(NewObject::dependent(DeviceBackend(d::DeviceBackend::new(d.0.clone())),vec![],vec![device.value.clone()])))});
object_binding!(backend_name,"num_native_device_backend_name","crate::backend::Backend::name","Explicit native backend specialization. Feature-disabled JIT/GPU builds return Unsupported; machine representations are normalized and source, row, fold-work and allocation budgets are checked.",true,(a:Handle<DeviceBackend>)->String,|ctx,a|{a.with(ctx,|a|Ok(crate::backend::Backend::name(&a.0).into()))});
object_binding!(backend_eval,"num_native_device_backend_eval","crate::backend::Backend::compile","Explicit native backend specialization. Feature-disabled JIT/GPU builds return Unsupported; machine representations are normalized and source, row, fold-work and allocation budgets are checked.",true,(a:Handle<DeviceBackend>,expression:String,names:Vec<String>,point:Vec<f64>)->f64,|ctx,a,expression,names,point|{a.with(ctx,|a|{if expression.len()>4096||names.len()>32||names.len()!=point.len(){return Err(api::error("device-backend source/arity budget exceeded"));}let session=crate::api::Session::with_rules(&[crate::rules::elementary()]).map_err(api::error)?;let compiled=session.with_native_compilation(&[expression],&names,&[],|g,r,i,_|{let tape=crate::backend::tape::lower(g,r,i,&[])?;if tape.insts.len()>512{return Err(crate::backend::BackendError::Codegen("device backend instruction budget exceeded".into()));}crate::backend::Backend::compile(&a.0,g,r[0],i)}).map_err(api::error)?;Ok(compiled.call(&point))})});
object_binding!(native_compile,"num_native_device_native_compile","crate::backend::device::native::NativeKernel::compile","Explicit native backend specialization. Feature-disabled JIT/GPU builds return Unsupported; machine representations are normalized and source, row, fold-work and allocation budgets are checked.",true,(a:Handle<KernelObject>)->Maybe<NewObject<NativeKernel>>,|ctx,a|{feature_jit!{a.with(ctx,|a|{if a.ir.ops.len()>128{return Err(api::error("native device expansion budget exceeded"));}let analysis=a.ir.analyse().map_err(api::error)?;Ok(Maybe(d::native::NativeKernel::compile(&a.ir,&analysis).map(|native|NewObject::normalized(NativeKernel{native,ir:a.ir.clone()}))))})}});
object_binding!(native_lower,"num_native_device_native_lower","crate::backend::device::native::lower","Explicit native backend specialization. Feature-disabled JIT/GPU builds return Unsupported; machine representations are normalized and source, row, fold-work and allocation budgets are checked.",true,(a:Handle<KernelObject>)->Maybe<NewObject<crate::rules::numerical_abi::native_lir::Module>>,|ctx,a|{feature_jit!{a.with(ctx,|a|{if a.ir.ops.len()>128{return Err(api::error("native device expansion budget exceeded"));}let analysis=a.ir.analyse().map_err(api::error)?;Ok(Maybe(d::native::lower(&a.ir,&analysis).map(|m|NewObject::normalized(crate::rules::numerical_abi::native_lir::Module(m)))))})}});
object_binding!(native_run,"num_native_device_native_run","crate::backend::device::native::NativeKernel::run","Explicit native backend specialization. Feature-disabled JIT/GPU builds return Unsupported; machine representations are normalized and source, row, fold-work and allocation budgets are checked.",true,(a:Handle<NativeKernel>,columns:Vec<Vec<f64>>,params:Vec<f64>,row_lo:usize,total_rows:usize,count:usize)->Dynamic,|ctx,a,columns,params,row_lo,total_rows,count|{feature_jit!{a.with(ctx,|a|{if row_lo.checked_add(count).is_none_or(|n|n>total_rows){return Err(api::error("native row interval exceeds launch"));}let cols=columns.iter().map(Vec::as_slice).collect::<Vec<_>>();work(&a.ir,total_rows,&columns.iter().map(Vec::len).collect::<Vec<_>>())?;check_columns(&a.ir,total_rows,&cols,&params)?;let mut output=vec![vec![0.0;count];a.ir.outputs.len()];let mut out=output.iter_mut().map(Vec::as_mut_slice).collect::<Vec<_>>();let completed=a.native.run(&cols,&params,&mut out,row_lo,total_rows);Ok(Dynamic(api::record(vec![("completed",Value::Bool(completed)),("values",Output::encode(output)?)])))})}});
object_binding!(gpu_probe,"num_native_device_gpu_probe","crate::backend::device::gpu::probe","Explicit native backend specialization. Feature-disabled JIT/GPU builds return Unsupported; machine representations are normalized and source, row, fold-work and allocation budgets are checked.",true,()->Many<NewObject<DeviceObject>>,|ctx,|{let _=ctx;feature_gpu!{let devices=d::gpu::probe();if devices.len()>8{return Err(api::error("GPU discovery count budget exceeded"));}Ok(Many(devices.into_iter().map(|d|NewObject::normalized(DeviceObject(d))).collect()))}});
object_binding!(gpu_new,"num_native_device_gpu_new","crate::backend::device::gpu::WgpuDevice::new","Explicit native backend specialization. Feature-disabled JIT/GPU builds return Unsupported; machine representations are normalized and source, row, fold-work and allocation budgets are checked.",true,(ordinal:usize)->NewObject<DeviceObject>,|_ctx,ordinal|{let _=ordinal;feature_gpu!{if ordinal>=8{return Err(api::error("GPU ordinal budget exceeded"));}let instance=wgpu::Instance::new(&wgpu::InstanceDescriptor::default());let adapters=instance.enumerate_adapters(wgpu::Backends::all());let adapter=adapters.get(ordinal).ok_or_else(||api::error("GPU adapter ordinal unavailable"))?;let device=d::gpu::WgpuDevice::new(adapter,ordinal).map_err(api::error)?;Ok(NewObject::normalized(DeviceObject(Arc::new(device))))}});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &backend::BINDING,
            &backend_name::BINDING,
            &backend_eval::BINDING,
            &native_compile::BINDING,
            &native_lower::BINDING,
            &native_run::BINDING,
            &gpu_probe::BINDING,
            &gpu_new::BINDING,
        ],
        "numerical_vector",
    )
}
