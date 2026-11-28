//! Typed handles for the real native loop IR and its bounded interpreter.
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
use crate::backend::lir::BlockId;
use crate::backend::lir::FuncId;
use crate::backend::lir::Inst;
use crate::backend::lir::Reg;
use crate::backend::lir::Ty;
use crate::backend::lir::{
    self as l,
};
use crate::backend::tape::Cmp;
use crate::graph::RuleError;
use crate::graph::rule::Installer;
use crate::kernel::OpError;
use crate::kernel::Value;
use crate::kernel::value::TensorValue;
struct Builder {
    native: Option<l::FunctionBuilder>,
    regs: usize,
    blocks: usize,
    insts: usize,
}
struct Function(l::Function);
pub(super) struct Module(pub(super) l::Module);
struct Instruction(Inst);
struct Outcome(l::Outcome);
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
native!(
    Builder,
    "rssn/native/lir-builder/v1",
    1_048_576,
    |a: &Builder| {
        if a.native.is_some() && a.regs <= 4096 && a.blocks <= 256 && a.insts <= 4096 {
            Ok(())
        } else {
            Err(api::error("LIR builder consumed or size budget exceeded"))
        }
    }
);
native!(
    Function,
    "rssn/native/lir-function/v1",
    1_048_576,
    |a: &Function| { check_function(&a.0) }
);
native!(
    Module,
    "rssn/native/lir-module/v1",
    4_194_304,
    |a: &Module| {
        if a.0.functions.len() > 16 {
            return Err(api::error("LIR function count budget exceeded"));
        }
        for f in &a.0.functions {
            check_function(f)?;
        }
        a.0.verify().map_err(api::error)
    }
);
native!(
    Instruction,
    "rssn/native/lir-instruction/v1",
    4096,
    |_: &Instruction| Ok(())
);
native!(
    Outcome,
    "rssn/native/lir-outcome/v1",
    65536,
    |_: &Outcome| Ok(())
);
fn check_function(f: &l::Function) -> Result<(), OpError> {
    if f.regs.len() > 4096
        || f.blocks.len() > 256
        || f.blocks.iter().map(|b| b.insts.len()).sum::<usize>() > 4096
        || f.results.len() > 64
    {
        return Err(api::error("LIR function size budget exceeded"));
    }
    Ok(())
}
fn ty(s: &str) -> Result<Ty, OpError> {
    match s {
        | "F64" => Ok(Ty::F64),
        | "I64" => Ok(Ty::I64),
        | "Bool" => Ok(Ty::Bool),
        | "Slice" => Ok(Ty::Slice),
        | _ => Err(api::error("LIR type must be F64/I64/Bool/Slice")),
    }
}
fn compare(s: &str) -> Result<Cmp, OpError> {
    match s {
        | "Lt" => Ok(Cmp::Lt),
        | "Le" => Ok(Cmp::Le),
        | "Gt" => Ok(Cmp::Gt),
        | "Ge" => Ok(Cmp::Ge),
        | "Eq" => Ok(Cmp::Eq),
        | "Ne" => Ok(Cmp::Ne),
        | _ => Err(api::error("unknown LIR comparison")),
    }
}
fn builder(a: &mut Builder) -> Result<&mut l::FunctionBuilder, OpError> {
    if a.regs >= 4096 || a.insts >= 4096 {
        return Err(api::error("LIR instruction/register budget exceeded"));
    }
    a.native
        .as_mut()
        .ok_or_else(|| api::error("builder was consumed"))
}
fn operands(
    b: &l::FunctionBuilder,
    inst: &Inst,
) -> Result<(), OpError> {
    let regs: Vec<Reg> = match inst {
        | Inst::F(_) | Inst::I(_) | Inst::B(_) => vec![],
        | Inst::Copy(a)
        | Inst::FNeg(a)
        | Inst::FPowI(a, _)
        | Inst::FUnary(_, a)
        | Inst::INeg(a)
        | Inst::Not(a)
        | Inst::IToF(a)
        | Inst::FToI(a)
        | Inst::Len(a)
        | Inst::Extern1(_, a) => vec![*a],
        | Inst::Select(a, b, c) | Inst::Store(a, b, c) => vec![*a, *b, *c],
        | Inst::Call(_, args) | Inst::Eval(_, args) => args.clone(),
        | Inst::FAdd(a, b)
        | Inst::FSub(a, b)
        | Inst::FMul(a, b)
        | Inst::FDiv(a, b)
        | Inst::FPow(a, b)
        | Inst::FAtan2(a, b)
        | Inst::FMin(a, b)
        | Inst::FMax(a, b)
        | Inst::IAdd(a, b)
        | Inst::ISub(a, b)
        | Inst::IMul(a, b)
        | Inst::IDiv(a, b)
        | Inst::IRem(a, b)
        | Inst::IMin(a, b)
        | Inst::IMax(a, b)
        | Inst::FCmp(_, a, b)
        | Inst::ICmp(_, a, b)
        | Inst::And(a, b)
        | Inst::Or(a, b)
        | Inst::Load(a, b)
        | Inst::Extern2(_, a, b) => vec![*a, *b],
    };
    if regs.iter().any(|r| b.ty(*r).is_none()) {
        return Err(api::error("instruction references an unknown register"));
    }
    Ok(())
}
fn instruction(
    kind: &str,
    args: &[u32],
    float: f64,
    integer: i64,
    aux: &str,
) -> Result<Inst, OpError> {
    let r = |i: usize| {
        args.get(i)
            .copied()
            .map(Reg)
            .ok_or_else(|| api::error("instruction operand missing"))
    };
    Ok(match kind {
        | "F" => Inst::F(float),
        | "I" => Inst::I(integer),
        | "B" => Inst::B(integer != 0),
        | "Copy" => Inst::Copy(r(0)?),
        | "FNeg" => Inst::FNeg(r(0)?),
        | "INeg" => Inst::INeg(r(0)?),
        | "Not" => Inst::Not(r(0)?),
        | "IToF" => Inst::IToF(r(0)?),
        | "FToI" => Inst::FToI(r(0)?),
        | "Len" => Inst::Len(r(0)?),
        | "FAdd" => Inst::FAdd(r(0)?, r(1)?),
        | "FSub" => Inst::FSub(r(0)?, r(1)?),
        | "FMul" => Inst::FMul(r(0)?, r(1)?),
        | "FDiv" => Inst::FDiv(r(0)?, r(1)?),
        | "FPow" => Inst::FPow(r(0)?, r(1)?),
        | "FPowI" => Inst::FPowI(r(0)?, i32::try_from(integer).map_err(api::error)?),
        | "FUnary" => Inst::FUnary(super::native_backend_scalar::intrinsic(aux)?, r(0)?),
        | "FAtan2" => Inst::FAtan2(r(0)?, r(1)?),
        | "FMin" => Inst::FMin(r(0)?, r(1)?),
        | "FMax" => Inst::FMax(r(0)?, r(1)?),
        | "IAdd" => Inst::IAdd(r(0)?, r(1)?),
        | "ISub" => Inst::ISub(r(0)?, r(1)?),
        | "IMul" => Inst::IMul(r(0)?, r(1)?),
        | "IDiv" => Inst::IDiv(r(0)?, r(1)?),
        | "IRem" => Inst::IRem(r(0)?, r(1)?),
        | "IMin" => Inst::IMin(r(0)?, r(1)?),
        | "IMax" => Inst::IMax(r(0)?, r(1)?),
        | "FCmp" => Inst::FCmp(compare(aux)?, r(0)?, r(1)?),
        | "ICmp" => Inst::ICmp(compare(aux)?, r(0)?, r(1)?),
        | "And" => Inst::And(r(0)?, r(1)?),
        | "Or" => Inst::Or(r(0)?, r(1)?),
        | "Load" => Inst::Load(r(0)?, r(1)?),
        | "Select" => Inst::Select(r(0)?, r(1)?, r(2)?),
        | "Store" => Inst::Store(r(0)?, r(1)?, r(2)?),
        | "Call" => {
            Inst::Call(
                FuncId(u32::try_from(integer).map_err(api::error)?),
                args.iter().copied().map(Reg).collect(),
            )
        },
        | _ => {
            return Err(api::error(
                "unknown LIR instruction; raw external pointers are not imported",
            ));
        },
    })
}
fn encoded(v: &l::Value<'_>) -> Result<Value, OpError> {
    Ok(match v {
        | l::Value::F(x) => Value::Float(*x),
        | l::Value::I(x) => Value::int(*x),
        | l::Value::B(x) => Value::Bool(*x),
        | l::Value::Read(x) => {
            use super::codecs::Output;
            Output::encode(x.to_vec())?
        },
        | l::Value::Write(x) => {
            use super::codecs::Output;
            Output::encode(x.to_vec())?
        },
    })
}
pub(super) enum Owned {
    F(f64),
    I(i64),
    B(bool),
    S(Vec<f64>),
}
pub(super) fn arguments(v: &Value) -> Result<Vec<Owned>, OpError> {
    let Value::Tensor(t) = v else {
        return Err(api::error("LIR args must be a vector"));
    };
    if t.shape != [t.data.len()] || t.data.len() > 64 {
        return Err(api::error("LIR argument vector invalid"));
    }
    t.data
        .iter()
        .map(|v| {
            match v {
                | Value::Float(x) => Ok(Owned::F(*x)),
                | Value::Int(x) => {
                    use num_traits::ToPrimitive;
                    Ok(Owned::I(
                        x.to_i64()
                            .ok_or_else(|| api::error("LIR integer must fit i64"))?,
                    ))
                },
                | Value::Bool(x) => Ok(Owned::B(*x)),
                | Value::Tensor(t) if t.shape == [t.data.len()] => {
                    api::normalized(v)?;
                    Ok(Owned::S(
                        t.data
                            .iter()
                            .map(|v| super::numeric(v, 0))
                            .collect::<Result<_, _>>()?,
                    ))
                },
                | _ => {
                    Err(api::error(
                        "LIR arguments accept normalized Float/Int/Bool/slice values",
                    ))
                },
            }
        })
        .collect()
}
pub(super) fn borrowed(v: &mut [Owned]) -> Vec<l::Value<'_>> {
    v.iter_mut()
        .map(|v| {
            match v {
                | Owned::F(x) => l::Value::F(*x),
                | Owned::I(x) => l::Value::I(*x),
                | Owned::B(x) => l::Value::B(*x),
                | Owned::S(x) => l::Value::Write(x),
            }
        })
        .collect()
}
pub(super) fn args_value(v: &[l::Value<'_>]) -> Result<Value, OpError> {
    Ok(Value::Tensor(TensorValue {
        shape: vec![v.len()],
        elem: api::dynamic_type(),
        data: v.iter().map(encoded).collect::<Result<_, _>>()?,
    }))
}
pub(super) fn outcome(v: &l::Outcome) -> Result<Value, OpError> {
    Ok(match v {
        | l::Outcome::Returned(v) => {
            api::record(vec![
                ("status", Value::String("Returned".into())),
                ("values", args_value(v)?),
            ])
        },
        | l::Outcome::Bailed(code) => {
            api::record(vec![
                ("status", Value::String("Bailed".into())),
                ("code", Value::int(*code)),
            ])
        },
    })
}
object_binding!(builder_new,"num_native_lir_builder_new","crate::backend::lir::FunctionBuilder::new","Typed native LIR adapter. Machine instructions use explicit normalized IEEE values, checked register/block/type references and bounded native interpreter steps.",true,(name:String,params:Vec<String>,results:Vec<String>)->NewObject<Builder>,|ctx,name,params,results|{let _=ctx;if name.len()>128||params.len()>64||results.len()>64{return Err(api::error("LIR name/arity budget exceeded"));}let p=params.iter().map(|v|ty(v)).collect::<Result<Vec<_>,_>>()?;let r=results.iter().map(|v|ty(v)).collect::<Result<Vec<_>,_>>()?;Ok(NewObject::normalized(Builder{native:Some(l::FunctionBuilder::new(&name,&p,&r)),regs:p.len(),blocks:1,insts:0}))});
object_binding!(instruction,"num_native_lir_instruction","","Typed native LIR adapter. Machine instructions use explicit normalized IEEE values, checked register/block/type references and bounded native interpreter steps.",true,(kind:String,args:Vec<u32>,float:f64,integer:i64,aux:String)->NewObject<Instruction>,|ctx,kind,args,float,integer,aux|{let _=ctx;if args.len()>64{return Err(api::error("LIR instruction arity budget exceeded"));}Ok(NewObject::normalized(Instruction(instruction(&kind,&args,float,integer,&aux)?)))});
object_binding!(param,"num_native_lir_param","crate::backend::lir::FunctionBuilder::param","Typed native LIR adapter. Machine instructions use explicit normalized IEEE values, checked register/block/type references and bounded native interpreter steps.",true,(a:Handle<Builder>,index:usize)->u32,|ctx,a,index|{let _=ctx;a.with_mut(ctx,|a|{let b=builder(a)?;let r=b.param(index);if b.ty(r).is_none(){return Err(api::error("parameter out of range"));}Ok(r.0)})});
object_binding!(reg,"num_native_lir_reg","crate::backend::lir::FunctionBuilder::reg","Typed native LIR adapter. Machine instructions use explicit normalized IEEE values, checked register/block/type references and bounded native interpreter steps.",true,(a:Handle<Builder>,kind:String)->u32,|ctx,a,kind|{let _=ctx;a.with_mut(ctx,|a|{let t=ty(&kind)?;let r=builder(a)?.reg(t);a.regs+=1;Ok(r.0)})});
object_binding!(ty,"num_native_lir_ty","crate::backend::lir::FunctionBuilder::ty","Typed native LIR adapter. Machine instructions use explicit normalized IEEE values, checked register/block/type references and bounded native interpreter steps.",true,(a:Handle<Builder>,reg:u32)->Maybe<String>,|ctx,a,reg|{let _=ctx;a.with_mut(ctx,|a|Ok(Maybe(builder(a)?.ty(Reg(reg)).map(|t|format!("{t:?}")))))});
object_binding!(new_block,"num_native_lir_new_block","crate::backend::lir::FunctionBuilder::new_block","Typed native LIR adapter. Machine instructions use explicit normalized IEEE values, checked register/block/type references and bounded native interpreter steps.",true,(a:Handle<Builder>)->u32,|ctx,a|{let _=ctx;a.with_mut(ctx,|a|{if a.blocks>=256{return Err(api::error("LIR block budget exceeded"));}let r=builder(a)?.new_block();a.blocks+=1;Ok(r.0)})});
object_binding!(block_param,"num_native_lir_block_param","crate::backend::lir::FunctionBuilder::block_param","Typed native LIR adapter. Machine instructions use explicit normalized IEEE values, checked register/block/type references and bounded native interpreter steps.",true,(a:Handle<Builder>,block:u32,kind:String)->u32,|ctx,a,block,kind|{let _=ctx;a.with_mut(ctx,|a|{if block as usize>=a.blocks{return Err(api::error("unknown LIR block"));}let r=builder(a)?.block_param(BlockId(block),ty(&kind)?);a.regs+=1;Ok(r.0)})});
object_binding!(switch_to,"num_native_lir_switch_to","crate::backend::lir::FunctionBuilder::switch_to","Typed native LIR adapter. Machine instructions use explicit normalized IEEE values, checked register/block/type references and bounded native interpreter steps.",true,(a:Handle<Builder>,block:u32)->bool,|ctx,a,block|{let _=ctx;a.with_mut(ctx,|a|{if block as usize>=a.blocks{return Err(api::error("unknown LIR block"));}builder(a)?.switch_to(BlockId(block));Ok(true)})});
object_binding!(current,"num_native_lir_current","crate::backend::lir::FunctionBuilder::current","Typed native LIR adapter. Machine instructions use explicit normalized IEEE values, checked register/block/type references and bounded native interpreter steps.",true,(a:Handle<Builder>)->u32,|ctx,a|{let _=ctx;a.with_mut(ctx,|a|{let b=builder(a)?;Ok(b.current().0)})});
object_binding!(is_terminated,"num_native_lir_is_terminated","crate::backend::lir::FunctionBuilder::is_terminated","Typed native LIR adapter. Machine instructions use explicit normalized IEEE values, checked register/block/type references and bounded native interpreter steps.",true,(a:Handle<Builder>)->bool,|ctx,a|{let _=ctx;a.with_mut(ctx,|a|{let b=builder(a)?;Ok(b.is_terminated())})});
object_binding!(ins,"num_native_lir_ins","crate::backend::lir::FunctionBuilder::ins","Typed native LIR adapter. Machine instructions use explicit normalized IEEE values, checked register/block/type references and bounded native interpreter steps.",true,(a:Handle<Builder>,inst:Handle<Instruction>)->u32,|ctx,a,inst|{let _=ctx;inst.with(ctx,|inst|a.with_mut(ctx,|a|{let b=builder(a)?;operands(b,&inst.0)?;let r=b.ins(inst.0.clone());a.regs+=1;a.insts+=1;Ok(r.0)}))});
object_binding!(ins_typed,"num_native_lir_ins_typed","crate::backend::lir::FunctionBuilder::ins_typed","Typed native LIR adapter. Machine instructions use explicit normalized IEEE values, checked register/block/type references and bounded native interpreter steps.",true,(a:Handle<Builder>,inst:Handle<Instruction>,kind:String)->u32,|ctx,a,inst,kind|{let _=ctx;inst.with(ctx,|inst|a.with_mut(ctx,|a|{let b=builder(a)?;operands(b,&inst.0)?;let r=b.ins_typed(inst.0.clone(),ty(&kind)?);a.regs+=1;a.insts+=1;Ok(r.0)}))});
object_binding!(assign,"num_native_lir_assign","crate::backend::lir::FunctionBuilder::assign","Typed native LIR adapter. Machine instructions use explicit normalized IEEE values, checked register/block/type references and bounded native interpreter steps.",true,(a:Handle<Builder>,reg:u32,inst:Handle<Instruction>)->bool,|ctx,a,reg,inst|{let _=ctx;inst.with(ctx,|inst|a.with_mut(ctx,|a|{let b=builder(a)?;if b.ty(Reg(reg)).is_none(){return Err(api::error("assignment destination unknown"));}operands(b,&inst.0)?;b.assign(Reg(reg),inst.0.clone());a.insts+=1;Ok(true)}))});
object_binding!(jump,"num_native_lir_jump","crate::backend::lir::FunctionBuilder::jump","Typed native LIR adapter. Machine instructions use explicit normalized IEEE values, checked register/block/type references and bounded native interpreter steps.",true,(a:Handle<Builder>,block:u32,args:Vec<u32>)->bool,|ctx,a,block,args|{let _=ctx;a.with_mut(ctx,|a|{if block as usize>=a.blocks{return Err(api::error("unknown jump target"));}let b=builder(a)?;if args.iter().any(|r|b.ty(Reg(*r)).is_none()){return Err(api::error("jump operand unknown"));}b.jump(BlockId(block),&args.into_iter().map(Reg).collect::<Vec<_>>());Ok(true)})});
object_binding!(branch,"num_native_lir_branch","crate::backend::lir::FunctionBuilder::branch","Typed native LIR adapter. Machine instructions use explicit normalized IEEE values, checked register/block/type references and bounded native interpreter steps.",true,(a:Handle<Builder>,condition:u32,yes:u32,yes_args:Vec<u32>,no:u32,no_args:Vec<u32>)->bool,|ctx,a,condition,yes,yes_args,no,no_args|{let _=ctx;a.with_mut(ctx,|a|{if yes as usize>=a.blocks||no as usize>=a.blocks{return Err(api::error("unknown branch target"));}let b=builder(a)?;if b.ty(Reg(condition))!=Some(Ty::Bool)||yes_args.iter().chain(&no_args).any(|r|b.ty(Reg(*r)).is_none()){return Err(api::error("branch operands/type invalid"));}b.branch(Reg(condition),BlockId(yes),&yes_args.into_iter().map(Reg).collect::<Vec<_>>(),BlockId(no),&no_args.into_iter().map(Reg).collect::<Vec<_>>());Ok(true)})});
object_binding!(ret,"num_native_lir_ret","crate::backend::lir::FunctionBuilder::ret","Typed native LIR adapter. Machine instructions use explicit normalized IEEE values, checked register/block/type references and bounded native interpreter steps.",true,(a:Handle<Builder>,registers:Vec<u32>)->bool,|ctx,a,registers|{let _=ctx;a.with_mut(ctx,|a|{let b=builder(a)?;if registers.iter().any(|r|b.ty(Reg(*r)).is_none()){return Err(api::error("return operand unknown"));}b.ret(&registers.into_iter().map(Reg).collect::<Vec<_>>());Ok(true)})});
object_binding!(bail,"num_native_lir_bail","crate::backend::lir::FunctionBuilder::bail","Typed native LIR adapter. Machine instructions use explicit normalized IEEE values, checked register/block/type references and bounded native interpreter steps.",true,(a:Handle<Builder>,code:u32)->bool,|ctx,a,code|{let _=ctx;a.with_mut(ctx,|a|{builder(a)?.bail(code);Ok(true)})});
object_binding!(finish,"num_native_lir_finish","crate::backend::lir::FunctionBuilder::finish","Typed native LIR adapter. Machine instructions use explicit normalized IEEE values, checked register/block/type references and bounded native interpreter steps.",true,(a:Handle<Builder>)->NewObject<Function>,|ctx,a|{let _=ctx;a.with_mut(ctx,|a|{let b=a.native.take().ok_or_else(||api::error("builder consumed"))?;Ok(NewObject::normalized(Function(b.finish().map_err(api::error)?)))})});
object_binding!(function_verify,"num_native_lir_function_verify","crate::backend::lir::Module::verify","Typed native LIR adapter. Machine instructions use explicit normalized IEEE values, checked register/block/type references and bounded native interpreter steps.",true,(a:Handle<Function>)->bool,|ctx,a|{let _=ctx;a.with(ctx,|a|{l::Module{functions:vec![a.0.clone()]}.verify().map_err(api::error)?;Ok(true)})});
object_binding!(module,"num_native_lir_module","","Typed native LIR adapter. Machine instructions use explicit normalized IEEE values, checked register/block/type references and bounded native interpreter steps.",true,(functions:Input)->NewObject<Module>,|ctx,functions|{let _=ctx;let Value::Tensor(t)=functions.0 else{return Err(api::error("expected function handles"));};
if t.shape!=[t.data.len()]||t.data.len()>16{return Err(api::error("function count budget exceeded"));}let functions=t.data.iter().map(|v|api::get::<Function,_>(ctx,v,Function::KIND,|f|Ok(f.0.clone()))).collect::<Result<_,_>>()?;Ok(NewObject::normalized(Module(l::Module{functions})))});
object_binding!(module_verify,"num_native_lir_module_verify","crate::backend::lir::Module::verify","Typed native LIR adapter. Machine instructions use explicit normalized IEEE values, checked register/block/type references and bounded native interpreter steps.",true,(a:Handle<Module>)->bool,|ctx,a|{let _=ctx;a.with(ctx,|a|{a.0.verify().map_err(api::error)?;Ok(true)})});
object_binding!(written_params,"num_native_lir_written_params","crate::backend::lir::Module::written_params","Typed native LIR adapter. Machine instructions use explicit normalized IEEE values, checked register/block/type references and bounded native interpreter steps.",true,(a:Handle<Module>,function:u32)->Many<bool>,|ctx,a,function|{let _=ctx;a.with(ctx,|a|{if function as usize>=a.0.functions.len(){return Err(api::error("function out of range"));}Ok(Many(a.0.written_params(FuncId(function))))})});
object_binding!(check_args,"num_native_lir_check_args","crate::backend::lir::Module::check_args","Typed native LIR adapter. Machine instructions use explicit normalized IEEE values, checked register/block/type references and bounded native interpreter steps.",true,(a:Handle<Module>,function:u32,args:Input)->bool,|ctx,a,function,args|{let _=ctx;a.with(ctx,|a|{let mut owned=arguments(&args.0)?;let args=borrowed(&mut owned);a.0.check_args(FuncId(function),&args).map_err(api::error)?;Ok(true)})});
object_binding!(interpret,"num_native_lir_interpret","crate::backend::lir::Module::interpret","Typed native LIR adapter. Machine instructions use explicit normalized IEEE values, checked register/block/type references and bounded native interpreter steps.",true,(a:Handle<Module>,function:u32,args:Input,max_steps:u64)->Dynamic,|ctx,a,function,args,max_steps|{let _=ctx;a.with(ctx,|a|{if max_steps==0||max_steps>10_000_000{return Err(api::error("LIR step budget must be1..10000000"));}let mut owned=arguments(&args.0)?;let mut args=borrowed(&mut owned);let result=a.0.interpret(FuncId(function),&mut args,max_steps).map_err(api::error)?;Ok(Dynamic(api::record(vec![("outcome",outcome(&result)?),("arguments",args_value(&args)?)])))})});
object_binding!(outcome,"num_native_lir_outcome","crate::backend::lir::Module::interpret","Typed native LIR adapter. Machine instructions use explicit normalized IEEE values, checked register/block/type references and bounded native interpreter steps.",true,(a:Handle<Module>,function:u32,args:Input,max_steps:u64)->NewObject<Outcome>,|ctx,a,function,args,max_steps|{let _=ctx;a.with(ctx,|a|{if max_steps==0||max_steps>10_000_000{return Err(api::error("LIR step budget must be1..10000000"));}let mut owned=arguments(&args.0)?;let mut args=borrowed(&mut owned);Ok(NewObject::normalized(Outcome(a.0.interpret(FuncId(function),&mut args,max_steps).map_err(api::error)?)))})});
object_binding!(outcome_values,"num_native_lir_outcome_values","crate::backend::lir::Outcome::values","Typed native LIR adapter. Machine instructions use explicit normalized IEEE values, checked register/block/type references and bounded native interpreter steps.",true,(a:Handle<Outcome>)->Maybe<Dynamic>,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(Maybe(a.0.values().map(|v|args_value(v).map(Dynamic)).transpose()?)))});
object_binding!(value_ty,"num_native_lir_value_ty","crate::backend::lir::Value<'_>::ty","Typed native LIR adapter. Machine instructions use explicit normalized IEEE values, checked register/block/type references and bounded native interpreter steps.",true,(value:Input)->String,|ctx,value|{let _=ctx;let wrapped=Value::Tensor(TensorValue{shape:vec![1],elem:api::dynamic_type(),data:vec![value.0]});let mut owned=arguments(&wrapped)?;let values=borrowed(&mut owned);Ok(format!("{:?}",values[0].ty()))});
object_binding!(batch_from_tape,"num_native_lir_batch_from_tape","crate::backend::lir::batch_from_tape","Typed native LIR adapter. Machine instructions use explicit normalized IEEE values, checked register/block/type references and bounded native interpreter steps.",true,(tape:Input)->Maybe<NewObject<Function>>,|ctx,tape|{let _=ctx;let tape=super::super::native_backend_extended::tape_snapshot(ctx,&tape.0)?;if tape.tape.insts.len()>512{return Err(api::error("batch LIR expansion budget exceeded"));}Ok(Maybe(l::batch_from_tape(&tape.tape).map(|f|NewObject::normalized(Function(f)))))});
object_binding!(f64_to_i64,"num_native_lir_f64_to_i64","crate::backend::lir::f64_to_i64","Typed native LIR adapter. Machine instructions use explicit normalized IEEE values, checked register/block/type references and bounded native interpreter steps.",true,(x:f64)->Maybe<i64>,|ctx,x|{let _=ctx;Ok(Maybe(l::f64_to_i64(x)))});
object_binding!(fmin,"num_native_lir_fmin","crate::backend::lir::fmin","Typed native LIR adapter. Machine instructions use explicit normalized IEEE values, checked register/block/type references and bounded native interpreter steps.",true,(a:f64,b:f64)->f64,|ctx,a,b|{let _=ctx;Ok(l::fmin(a,b))});
object_binding!(fmax,"num_native_lir_fmax","crate::backend::lir::fmax","Typed native LIR adapter. Machine instructions use explicit normalized IEEE values, checked register/block/type references and bounded native interpreter steps.",true,(a:f64,b:f64)->f64,|ctx,a,b|{let _=ctx;Ok(l::fmax(a,b))});
object_binding!(fcmp,"num_native_lir_fcmp","crate::backend::lir::fcmp","Typed native LIR adapter. Machine instructions use explicit normalized IEEE values, checked register/block/type references and bounded native interpreter steps.",true,(kind:String,a:f64,b:f64)->bool,|ctx,kind,a,b|{let _=ctx;Ok(l::fcmp(compare(&kind)?,a,b))});
object_binding!(icmp,"num_native_lir_icmp","crate::backend::lir::icmp","Typed native LIR adapter. Machine instructions use explicit normalized IEEE values, checked register/block/type references and bounded native interpreter steps.",true,(kind:String,a:i64,b:i64)->bool,|ctx,kind,a,b|{let _=ctx;Ok(l::icmp(compare(&kind)?,a,b))});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &builder_new::BINDING,
            &instruction::BINDING,
            &param::BINDING,
            &reg::BINDING,
            &ty::BINDING,
            &new_block::BINDING,
            &block_param::BINDING,
            &switch_to::BINDING,
            &current::BINDING,
            &is_terminated::BINDING,
            &ins::BINDING,
            &ins_typed::BINDING,
            &assign::BINDING,
            &jump::BINDING,
            &branch::BINDING,
            &ret::BINDING,
            &bail::BINDING,
            &finish::BINDING,
            &function_verify::BINDING,
            &module::BINDING,
            &module_verify::BINDING,
            &written_params::BINDING,
            &check_args::BINDING,
            &interpret::BINDING,
            &outcome::BINDING,
            &outcome_values::BINDING,
            &value_ty::BINDING,
            &batch_from_tape::BINDING,
            &f64_to_i64::BINDING,
            &fmin::BINDING,
            &fmax::BINDING,
            &fcmp::BINDING,
            &icmp::BINDING,
        ],
        "numerical_vector",
    )
}
