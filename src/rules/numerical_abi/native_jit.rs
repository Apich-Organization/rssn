//! Explicit native JIT objects. Disabled builds report the required feature.
#![cfg_attr(
    not(feature = "jit"),
    allow(unused_variables, dead_code, unused_imports)
)]
use super::native_backend_extended::TapeSnapshot;
use super::native_backend_extended::{
    self as execution,
};
use super::native_object_api::Dynamic;
use super::native_object_api::Handle;
use super::native_object_api::Input;
use super::native_object_api::Maybe;
use super::native_object_api::NativeType;
use super::native_object_api::NewObject;
use super::native_object_api::object_binding;
use super::native_object_api::{
    self as api,
};
use crate::backend::CompiledMulti;
use crate::graph::RuleError;
use crate::graph::rule::Installer;
use crate::kernel::CallCtx;
use crate::kernel::OpError;
use crate::kernel::Value;
use std::sync::Arc;
struct BackendObject {
    #[cfg(feature = "jit")]
    native: crate::backend::jit::CraneliftBackend,
    compilations: usize,
}
struct Function {
    #[cfg(feature = "jit")]
    native: Arc<crate::backend::jit::JitFunction>,
    frame: TapeSnapshot,
}
struct Tail {
    #[cfg(feature = "jit")]
    native: crate::backend::jit::TailLoop,
    frame: TapeSnapshot,
}
struct Lir {
    #[cfg(feature = "jit")]
    native: crate::backend::lir_jit::LirModule,
}
macro_rules! native {
    ($t:ty,$kind:literal,$bytes:expr) => {
        impl NativeType for $t {
            const KIND: &'static str = $kind;
            const MAX_BYTES: usize = $bytes;

            fn validate(&self) -> Result<(), OpError> {
                Ok(())
            }
        }
    };
}
native!(BackendObject, "rssn/native/jit-backend/v1", 16_777_216);
native!(Function, "rssn/native/jit-function/v1", 4_194_304);
native!(Tail, "rssn/native/jit-tail-loop/v1", 4_194_304);
native!(Lir, "rssn/native/lir-jit/v1", 4_194_304);
macro_rules! jit{($($body:tt)*)=>{{#[cfg(feature="jit")]{$($body)*}#[cfg(not(feature="jit"))]{Err(OpError::Unsupported{backend:"rssn".into(),feature:"requires RSSN jit feature".into()})}}};}
fn tape(
    ctx: &CallCtx<'_>,
    v: &Value,
) -> Result<TapeSnapshot, OpError> {
    let s = execution::tape_snapshot(ctx, v)?;
    if s.tape.insts.len() > 512 {
        return Err(api::error(
            "native code generation budget exceeds512instructions",
        ));
    }
    Ok(s)
}
#[cfg(feature = "jit")]
fn straight_line(module: &crate::backend::lir::Module) -> Result<(), OpError> {
    use crate::backend::lir::Inst;
    use crate::backend::lir::Term;
    for f in &module.functions {
        for (at, b) in f.blocks.iter().enumerate() {
            let forward = |id: crate::backend::lir::BlockId| id.0 as usize > at;
            match &b.term {
                | Term::Jump(id, _) if !forward(*id) => {
                    return Err(api::error(
                        "bounded JIT specialization requires an acyclic CFG; use native interpreter for loops",
                    ));
                },
                | Term::Branch(_, (a, _), (b, _)) if !forward(*a) || !forward(*b) => {
                    return Err(api::error(
                        "bounded JIT specialization requires an acyclic CFG",
                    ));
                },
                | _ => {},
            }
            if b.insts.iter().any(|(_, i)| matches!(i, Inst::Call(..))) {
                return Err(api::error(
                    "bounded standalone JIT specialization excludes recursive/module calls; use native interpreter",
                ));
            }
        }
    }
    Ok(())
}
object_binding!(backend_new,"num_native_jit_backend_new","crate::backend::jit::CraneliftBackend::new","Actual feature-gated native JIT adapter. Requires RSSN jit; disabled builds return Unsupported. Compiled expression calls retain physical tape frames; low-level LIR has explicit normalized values and bounded acyclic control flow.",false,()->NewObject<BackendObject>,|ctx,|{let _=ctx;jit!{Ok(NewObject::normalized(BackendObject{native:crate::backend::jit::CraneliftBackend::new(),compilations:0}))}});
object_binding!(backend_with_key,"num_native_jit_backend_with_key","crate::backend::jit::CraneliftBackend::with_key","Actual feature-gated native JIT adapter. Requires RSSN jit; disabled builds return Unsupported. Compiled expression calls retain physical tape frames; low-level LIR has explicit normalized values and bounded acyclic control flow.",false,(key:u64)->NewObject<BackendObject>,|_ctx,key|{jit!{Ok(NewObject::normalized(BackendObject{native:crate::backend::jit::CraneliftBackend::with_key(key),compilations:0}))}});
object_binding!(backend_key_extra,"num_native_jit_backend_key_extra","crate::backend::jit::CraneliftBackend::key_extra","Actual feature-gated native JIT adapter. Requires RSSN jit; disabled builds return Unsupported. Compiled expression calls retain physical tape frames; low-level LIR has explicit normalized values and bounded acyclic control flow.",false,(a:Handle<BackendObject>)->u64,|ctx,a|{jit!{a.with(ctx,|a|Ok(a.native.key_extra()))}});
object_binding!(backend_cached,"num_native_jit_backend_cached","crate::backend::jit::CraneliftBackend::cached","Actual feature-gated native JIT adapter. Requires RSSN jit; disabled builds return Unsupported. Compiled expression calls retain physical tape frames; low-level LIR has explicit normalized values and bounded acyclic control flow.",false,(a:Handle<BackendObject>)->usize,|ctx,a|{jit!{a.with(ctx,|a|Ok(a.native.cached()))}});
object_binding!(backend_set_key,"num_native_jit_backend_set_key","crate::backend::jit::CraneliftBackend::set_key_extra","Actual feature-gated native JIT adapter. Requires RSSN jit; disabled builds return Unsupported. Compiled expression calls retain physical tape frames; low-level LIR has explicit normalized values and bounded acyclic control flow.",false,(a:Handle<BackendObject>,key:u64)->bool,|ctx,a,key|{jit!{a.with_mut(ctx,|a|{a.native.set_key_extra(key);Ok(true)})}});
object_binding!(backend_cache_key,"num_native_jit_backend_cache_key","crate::backend::jit::CraneliftBackend::cache_key","Actual feature-gated native JIT adapter. Requires RSSN jit; disabled builds return Unsupported. Compiled expression calls retain physical tape frames; low-level LIR has explicit normalized values and bounded acyclic control flow.",false,(a:Handle<BackendObject>,program:Input)->u64,|ctx,a,program|{jit!{{let s=tape(ctx,&program.0)?;a.with(ctx,|a|Ok(a.native.cache_key(&s.tape)))}}});
object_binding!(backend_compile,"num_native_jit_backend_compile","crate::backend::jit::CraneliftBackend::compile_tape","Actual feature-gated native JIT adapter. Requires RSSN jit; disabled builds return Unsupported. Compiled expression calls retain physical tape frames; low-level LIR has explicit normalized values and bounded acyclic control flow.",false,(a:Handle<BackendObject>,program:Input)->NewObject<Function>,|ctx,a,program|{jit!{{let s=tape(ctx,&program.0)?;a.with_mut(ctx,|a|{if a.compilations>=8{return Err(api::error("native backend compilation count budget exceeded"));}let native=a.native.compile_tape(&s.tape).map_err(api::error)?;a.compilations+=1;Ok(NewObject::new(Function{native,frame:s.clone()},execution::snapshot_dimensions(&s)))})}}});
object_binding!(compile,"num_native_jit_compile","crate::backend::jit::JitFunction::compile","Actual feature-gated native JIT adapter. Requires RSSN jit; disabled builds return Unsupported. Compiled expression calls retain physical tape frames; low-level LIR has explicit normalized values and bounded acyclic control flow.",false,(program:Input)->NewObject<Function>,|ctx,program|{jit!{{let s=tape(ctx,&program.0)?;let native=Arc::new(crate::backend::jit::JitFunction::compile(&s.tape).map_err(api::error)?);let d=execution::snapshot_dimensions(&s);Ok(NewObject::new(Function{native,frame:s},d))}}});
object_binding!(eval,"num_native_jit_eval","crate::backend::CompiledMulti::eval","Actual feature-gated native JIT adapter. Requires RSSN jit; disabled builds return Unsupported. Compiled expression calls retain physical tape frames; low-level LIR has explicit normalized values and bounded acyclic control flow.",false,(a:Handle<Function>,x:Input,p:Input)->Dynamic,|ctx,a,x,p|{jit!{a.with(ctx,|a|{let(x,p)=execution::snapshot_arguments(ctx,&a.frame,&x.0,&p.0)?;let mut out=vec![0.0;a.frame.output_scales.len()];a.native.eval(&x,&p,&mut out);execution::snapshot_outputs(out,&a.frame)})}});
object_binding!(target_key,"num_native_jit_target_key","crate::backend::jit::target_key","Actual feature-gated native JIT adapter. Requires RSSN jit; disabled builds return Unsupported. Compiled expression calls retain physical tape frames; low-level LIR has explicit normalized values and bounded acyclic control flow.",false,()->u64,|ctx,|{let _=ctx;jit!{Ok(crate::backend::jit::target_key())}});
object_binding!(is_compiled,"num_native_jit_is_compiled","crate::backend::jit::is_compiled","Actual feature-gated native JIT adapter. Requires RSSN jit; disabled builds return Unsupported. Compiled expression calls retain physical tape frames; low-level LIR has explicit normalized values and bounded acyclic control flow.",false,(a:Handle<Function>)->bool,|ctx,a|{jit!{a.with(ctx,|a|Ok(crate::backend::jit::is_compiled(&*a.native)))}});
object_binding!(tail_compile,"num_native_jit_tail_compile","crate::backend::jit::TailLoop::compile","Actual feature-gated native JIT adapter. Requires RSSN jit; disabled builds return Unsupported. Compiled expression calls retain physical tape frames; low-level LIR has explicit normalized values and bounded acyclic control flow.",false,(program:Input)->NewObject<Tail>,|ctx,program|{jit!{{let s=tape(ctx,&program.0)?;if s.output_scales.len()!=s.input_scales.len()+2||s.output_scales[..s.input_scales.len()]!=s.input_scales||s.output_dimensions[..s.input_scales.len()]!=s.input_dimensions||s.output_dimensions.last()!=Some(&crate::kernel::signature::DimExpr::Dimless){return Err(api::error("tail loop next-state references must exactly match input frame and final completion flag must be dimensionless"));}let native=crate::backend::jit::TailLoop::compile(&s.tape).map_err(api::error)?;let d=execution::snapshot_dimensions(&s);Ok(NewObject::new(Tail{native,frame:s},d))}}});
object_binding!(tail_run,"num_native_jit_tail_run","crate::backend::jit::TailLoop::run","Actual feature-gated native JIT adapter. Requires RSSN jit; disabled builds return Unsupported. Compiled expression calls retain physical tape frames; low-level LIR has explicit normalized values and bounded acyclic control flow.",false,(a:Handle<Tail>,x:Input,fuel:u64)->Maybe<Dynamic>,|ctx,a,x,fuel|{jit!{a.with(ctx,|a|{if fuel>100_000{return Err(api::error("tail jump fuel budget exceeded"));}let(x,_)=execution::snapshot_arguments(ctx,&a.frame,&x.0,&execution::empty_arguments())?;let result=a.native.run(&x,fuel);Ok(Maybe(result.map(|r|Ok::<_,OpError>(Dynamic(api::record(vec![("outputs",execution::snapshot_outputs(r.outputs,&a.frame)?.0),("remaining",Value::int(r.remaining))])))).transpose()?))})}});
object_binding!(lir_compile,"num_native_jit_lir_compile","crate::backend::lir_jit::LirModule::compile","Actual feature-gated native JIT adapter. Requires RSSN jit; disabled builds return Unsupported. Compiled expression calls retain physical tape frames; low-level LIR has explicit normalized values and bounded acyclic control flow.",false,(module:Handle<crate::rules::numerical_abi::native_lir::Module>)->NewObject<Lir>,|ctx,module|{jit!{module.with(ctx,|m|{straight_line(&m.0)?;Ok(NewObject::normalized(Lir{native:crate::backend::lir_jit::LirModule::compile(&m.0).map_err(api::error)?}))})}});
object_binding!(lir_module,"num_native_jit_lir_module","crate::backend::lir_jit::LirModule::module","Actual feature-gated native JIT adapter. Requires RSSN jit; disabled builds return Unsupported. Compiled expression calls retain physical tape frames; low-level LIR has explicit normalized values and bounded acyclic control flow.",false,(a:Handle<Lir>)->NewObject<crate::rules::numerical_abi::native_lir::Module>,|ctx,a|{jit!{a.with(ctx,|a|Ok(NewObject::normalized(crate::rules::numerical_abi::native_lir::Module(a.native.module().clone()))))}});
object_binding!(lir_call,"num_native_jit_lir_call","crate::backend::lir_jit::LirModule::call","Actual feature-gated native JIT adapter. Requires RSSN jit; disabled builds return Unsupported. Compiled expression calls retain physical tape frames; low-level LIR has explicit normalized values and bounded acyclic control flow.",false,(a:Handle<Lir>,function:u32,args:Input)->Dynamic,|ctx,a,function,args|{jit!{a.with(ctx,|a|{let mut owned=crate::rules::numerical_abi::native_lir::arguments(&args.0)?;let mut args=crate::rules::numerical_abi::native_lir::borrowed(&mut owned);let result=a.native.call(crate::backend::lir::FuncId(function),&mut args).map_err(api::error)?;Ok(Dynamic(api::record(vec![("outcome",crate::rules::numerical_abi::native_lir::outcome(&result)?),("arguments",crate::rules::numerical_abi::native_lir::args_value(&args)?)])))})}});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &backend_new::BINDING,
            &backend_with_key::BINDING,
            &backend_key_extra::BINDING,
            &backend_cached::BINDING,
            &backend_set_key::BINDING,
            &backend_cache_key::BINDING,
            &backend_compile::BINDING,
            &compile::BINDING,
            &eval::BINDING,
            &target_key::BINDING,
            &is_compiled::BINDING,
            &tail_compile::BINDING,
            &tail_run::BINDING,
            &lir_compile::BINDING,
            &lir_module::BINDING,
            &lir_call::BINDING,
        ],
        "numerical_vector",
    )
}
