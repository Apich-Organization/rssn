//! Explicit native Float, hyperdual and interval specializations of backend scalar APIs.
//! These are normalized IEEE primitives; physical tape execution uses reference frames.
use super::native_object_api::Dynamic;
use super::native_object_api::Handle;
use super::native_object_api::Input;
use super::native_object_api::NativeArg;
use super::native_object_api::NativeOutput;
use super::native_object_api::NativeType;
use super::native_object_api::NewObject;
use super::native_object_api::object_binding;
use super::native_object_api::{
    self as api,
};
use crate::backend::BackendError;
use crate::backend::Intrinsic;
use crate::backend::tape::Cmp;
use crate::backend::tape::FromExtern;
use crate::backend::tape::HyperDual;
use crate::backend::tape::Interval;
use crate::backend::tape::Scalar;
use crate::backend::tape::{
    self,
};
use crate::graph::Arity;
use crate::graph::RuleError;
use crate::graph::rule::Installer;
use crate::kernel::CallCtx;
use crate::kernel::OpError;
use crate::kernel::Value;
use crate::kernel::signature::TypeExpr;
fn raw(v: &Value) -> Option<f64> {
    match v {
        | Value::Float(x) => Some(*x),
        | Value::Dyn(d) => raw(&d.value),
        | _ => None,
    }
}
#[derive(Clone, Copy)]
struct Raw(f64);
impl NativeArg for Raw {
    fn ty() -> TypeExpr {
        TypeExpr::Float
    }

    fn decode(
        v: &Value,
        _: usize,
        _: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        api::normalized(v)?;
        Ok(Self(if let Some(x) = raw(v) {
            x
        } else {
            super::numeric(v, 0)?
        }))
    }
}
impl NativeOutput for Raw {
    fn ty() -> TypeExpr {
        TypeExpr::Float
    }

    fn encode(
        self,
        _: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        Ok(Value::Float(self.0))
    }
}
struct Function(Intrinsic);
impl NativeArg for Function {
    fn ty() -> TypeExpr {
        TypeExpr::String
    }

    fn decode(
        v: &Value,
        _: usize,
        _: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        let Value::String(v) = v else {
            return Err(api::error("expected intrinsic name"));
        };
        Ok(Self(match v.as_str() {
            | "Exp" => Intrinsic::Exp,
            | "Ln" => Intrinsic::Ln,
            | "Sin" => Intrinsic::Sin,
            | "Cos" => Intrinsic::Cos,
            | "Tan" => Intrinsic::Tan,
            | "Asin" => Intrinsic::Asin,
            | "Acos" => Intrinsic::Acos,
            | "Atan" => Intrinsic::Atan,
            | "Sinh" => Intrinsic::Sinh,
            | "Cosh" => Intrinsic::Cosh,
            | "Tanh" => Intrinsic::Tanh,
            | "Asinh" => Intrinsic::Asinh,
            | "Acosh" => Intrinsic::Acosh,
            | "Atanh" => Intrinsic::Atanh,
            | "Sqrt" => Intrinsic::Sqrt,
            | "Abs" => Intrinsic::Abs,
            | "Floor" => Intrinsic::Floor,
            | "Ceil" => Intrinsic::Ceil,
            | "Cbrt" => Intrinsic::Cbrt,
            | _ => return Err(api::error("unknown native intrinsic")),
        }))
    }
}
struct Comparison(Cmp);
impl NativeArg for Comparison {
    fn ty() -> TypeExpr {
        TypeExpr::String
    }

    fn decode(
        v: &Value,
        _: usize,
        _: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        let Value::String(v) = v else {
            return Err(api::error("expected comparison name"));
        };
        Ok(Self(match v.as_str() {
            | "Lt" => Cmp::Lt,
            | "Le" => Cmp::Le,
            | "Gt" => Cmp::Gt,
            | "Ge" => Cmp::Ge,
            | "Eq" => Cmp::Eq,
            | "Ne" => Cmp::Ne,
            | _ => return Err(api::error("unknown comparison")),
        }))
    }
}
impl NativeType for HyperDual {
    const KIND: &'static str = "rssn/native/backend/hyperdual/v1";
    const MAX_BYTES: usize = 128;

    fn validate(&self) -> Result<(), OpError> {
        Ok(())
    }
}
impl NativeType for Interval {
    const KIND: &'static str = "rssn/native/backend/interval/v1";
    const MAX_BYTES: usize = 128;

    fn validate(&self) -> Result<(), OpError> {
        if self.lo.is_nan() || self.hi.is_nan() || self.lo <= self.hi {
            Ok(())
        } else {
            Err(api::error("interval bounds reversed"))
        }
    }
}
fn eval_function(
    name: &str,
    n: usize,
) -> Result<crate::graph::op::EvalFn, OpError> {
    if name.len() > 128 {
        return Err(api::error("function name too long"));
    }
    let session =
        crate::api::Session::with_rules(&[crate::rules::elementary()]).map_err(api::error)?;
    session
        .with_native_compilation(&["0".into()], &[], &[], |g, _, _, _| {
            let (_, d) = g
                .ops()
                .iter()
                .find(|(_, d)| d.name.as_ref() == name)
                .ok_or_else(|| BackendError::NoSemantics(name.into()))?;
            let accepts = match d.arity {
                | Arity::Fixed(k) => usize::from(k) == n,
                | Arity::Variadic => true,
            };
            if !accepts {
                return Err(BackendError::Codegen(
                    "scalar function arity mismatch".into(),
                ));
            }
            d.eval.ok_or_else(|| BackendError::NoSemantics(name.into()))
        })
        .map_err(api::error)
}
object_binding!(hyperdual_new,"num_native_hyperdual_new","crate::backend::tape::HyperDual::new","Native normalized IEEE scalar operation. NaN/interval-entire sentinels retain their native meaning; no physical quantity units are erased.",true,(value:Raw,first:Raw,second:Raw,mixed:Raw)->NewObject<HyperDual>,|ctx,value,first,second,mixed|{let _=ctx;Ok(NewObject::normalized(HyperDual::new(value.0,first.0,second.0,mixed.0)))});
object_binding!(hyperdual_value,"num_native_hyperdual_value","crate::backend::tape::HyperDual::value","Native normalized IEEE scalar operation. NaN/interval-entire sentinels retain their native meaning; no physical quantity units are erased.",true,(value:Handle<HyperDual>)->Raw,|ctx,value|{let _=ctx;value.with(ctx,|v|Ok(Raw(v.value())))});
object_binding!(hyperdual_first,"num_native_hyperdual_first","crate::backend::tape::HyperDual::first","Native normalized IEEE scalar operation. NaN/interval-entire sentinels retain their native meaning; no physical quantity units are erased.",true,(value:Handle<HyperDual>)->Raw,|ctx,value|{let _=ctx;value.with(ctx,|v|Ok(Raw(v.first())))});
object_binding!(hyperdual_second,"num_native_hyperdual_second","crate::backend::tape::HyperDual::second","Native normalized IEEE scalar operation. NaN/interval-entire sentinels retain their native meaning; no physical quantity units are erased.",true,(value:Handle<HyperDual>)->Raw,|ctx,value|{let _=ctx;value.with(ctx,|v|Ok(Raw(v.second())))});
object_binding!(hyperdual_mixed,"num_native_hyperdual_mixed","crate::backend::tape::HyperDual::mixed","Native normalized IEEE scalar operation. NaN/interval-entire sentinels retain their native meaning; no physical quantity units are erased.",true,(value:Handle<HyperDual>)->Raw,|ctx,value|{let _=ctx;value.with(ctx,|v|Ok(Raw(v.mixed())))});
object_binding!(interval_point,"num_native_interval_point","crate::backend::tape::Interval::point","Native normalized IEEE scalar operation. NaN/interval-entire sentinels retain their native meaning; no physical quantity units are erased.",true,(x:Raw)->NewObject<Interval>,|ctx,x|{let _=ctx;Ok(NewObject::normalized(Interval::point(x.0)))});
object_binding!(interval_new,"num_native_interval_new","","Native normalized IEEE scalar operation. NaN/interval-entire sentinels retain their native meaning; no physical quantity units are erased.",true,(lo:Raw,hi:Raw)->NewObject<Interval>,|ctx,lo,hi|{let _=ctx;Ok(NewObject::normalized(Interval{lo:lo.0,hi:hi.0}))});
object_binding!(interval_bounds,"num_native_interval_bounds","","Native normalized IEEE scalar operation. NaN/interval-entire sentinels retain their native meaning; no physical quantity units are erased.",true,(value:Handle<Interval>)->Dynamic,|ctx,value|{let _=ctx;value.with(ctx,|v|Ok(Dynamic(api::record(vec![("lo",Value::Float(v.lo)),("hi",Value::Float(v.hi))]))))});
object_binding!(interval_contains,"num_native_interval_contains","crate::backend::tape::Interval::contains","Native normalized IEEE scalar operation. NaN/interval-entire sentinels retain their native meaning; no physical quantity units are erased.",true,(value:Handle<Interval>,x:Raw)->bool,|ctx,value,x|{let _=ctx;value.with(ctx,|v|Ok(v.contains(x.0)))});
object_binding!(interval_is_point,"num_native_interval_is_point","crate::backend::tape::Interval::is_point","Native normalized IEEE scalar operation. NaN/interval-entire sentinels retain their native meaning; no physical quantity units are erased.",true,(value:Handle<Interval>)->bool,|ctx,value|{let _=ctx;value.with(ctx,|v|Ok(v.is_point()))});
object_binding!(intrinsic_apply,"num_native_intrinsic_apply","crate::backend::tape::Intrinsic::apply","Native normalized IEEE scalar operation. NaN/interval-entire sentinels retain their native meaning; no physical quantity units are erased.",true,(function:Function,x:Raw)->Raw,|ctx,function,x|{let _=ctx;Ok(Raw(function.0.apply(x.0)))});
object_binding!(intrinsic_derivative,"num_native_intrinsic_derivative","crate::backend::tape::Intrinsic::derivative","Native normalized IEEE scalar operation. NaN/interval-entire sentinels retain their native meaning; no physical quantity units are erased.",true,(function:Function,x:Raw)->Raw,|ctx,function,x|{let _=ctx;Ok(Raw(function.0.derivative(x.0)))});
object_binding!(intrinsic_second_derivative,"num_native_intrinsic_second_derivative","crate::backend::tape::Intrinsic::second_derivative","Native normalized IEEE scalar operation. NaN/interval-entire sentinels retain their native meaning; no physical quantity units are erased.",true,(function:Function,x:Raw)->Raw,|ctx,function,x|{let _=ctx;Ok(Raw(function.0.second_derivative(x.0)))});
object_binding!(scalar_float_nan,"num_native_scalar_float_nan","crate::backend::tape::Scalar::nan","Native normalized IEEE scalar operation. NaN/interval-entire sentinels retain their native meaning; no physical quantity units are erased.",true,()->Raw,|ctx,|{let _=ctx;Ok(Raw(<f64 as Scalar>::nan()))});
object_binding!(scalar_float_constant,"num_native_scalar_float_constant","crate::backend::tape::Scalar::constant","Native normalized IEEE scalar operation. NaN/interval-entire sentinels retain their native meaning; no physical quantity units are erased.",true,(x:Raw)->Raw,|ctx,x|{let _=ctx;Ok(Raw(<f64 as Scalar>::constant(x.0)))});
object_binding!(scalar_float_to_f64,"num_native_scalar_float_to_f64","crate::backend::tape::Scalar::to_f64","Native normalized IEEE scalar operation. NaN/interval-entire sentinels retain their native meaning; no physical quantity units are erased.",true,(x:Raw)->Raw,|ctx,x|{let _=ctx;Ok(Raw(<f64 as Scalar>::to_f64(x.0)))});
object_binding!(scalar_float_add,"num_native_scalar_float_add","crate::backend::tape::Scalar::add","Native normalized IEEE scalar operation. NaN/interval-entire sentinels retain their native meaning; no physical quantity units are erased.",true,(a:Raw,b:Raw)->Raw,|ctx,a,b|{let _=ctx;Ok(Raw(<f64 as Scalar>::add(a.0,b.0)))});
object_binding!(scalar_float_sub,"num_native_scalar_float_sub","crate::backend::tape::Scalar::sub","Native normalized IEEE scalar operation. NaN/interval-entire sentinels retain their native meaning; no physical quantity units are erased.",true,(a:Raw,b:Raw)->Raw,|ctx,a,b|{let _=ctx;Ok(Raw(<f64 as Scalar>::sub(a.0,b.0)))});
object_binding!(scalar_float_mul,"num_native_scalar_float_mul","crate::backend::tape::Scalar::mul","Native normalized IEEE scalar operation. NaN/interval-entire sentinels retain their native meaning; no physical quantity units are erased.",true,(a:Raw,b:Raw)->Raw,|ctx,a,b|{let _=ctx;Ok(Raw(<f64 as Scalar>::mul(a.0,b.0)))});
object_binding!(scalar_float_div,"num_native_scalar_float_div","crate::backend::tape::Scalar::div","Native normalized IEEE scalar operation. NaN/interval-entire sentinels retain their native meaning; no physical quantity units are erased.",true,(a:Raw,b:Raw)->Raw,|ctx,a,b|{let _=ctx;Ok(Raw(<f64 as Scalar>::div(a.0,b.0)))});
object_binding!(scalar_float_neg,"num_native_scalar_float_neg","crate::backend::tape::Scalar::neg","Native normalized IEEE scalar operation. NaN/interval-entire sentinels retain their native meaning; no physical quantity units are erased.",true,(a:Raw)->Raw,|ctx,a|{let _=ctx;Ok(Raw(<f64 as Scalar>::neg(a.0)))});
object_binding!(scalar_float_powi,"num_native_scalar_float_powi","crate::backend::tape::Scalar::powi","Native normalized IEEE scalar operation. NaN/interval-entire sentinels retain their native meaning; no physical quantity units are erased.",true,(a:Raw,n:i32)->Raw,|ctx,a,n|{let _=ctx;Ok(Raw(<f64 as Scalar>::powi(a.0,n)))});
object_binding!(scalar_float_powf,"num_native_scalar_float_powf","crate::backend::tape::Scalar::powf","Native normalized IEEE scalar operation. NaN/interval-entire sentinels retain their native meaning; no physical quantity units are erased.",true,(a:Raw,b:Raw)->Raw,|ctx,a,b|{let _=ctx;Ok(Raw(<f64 as Scalar>::powf(a.0,b.0)))});
object_binding!(scalar_float_unary,"num_native_scalar_float_unary","crate::backend::tape::Scalar::unary","Native normalized IEEE scalar operation. NaN/interval-entire sentinels retain their native meaning; no physical quantity units are erased.",true,(a:Raw,function:Function)->Raw,|ctx,a,function|{let _=ctx;Ok(Raw(<f64 as Scalar>::unary(a.0,function.0)))});
object_binding!(scalar_float_atan2,"num_native_scalar_float_atan2","crate::backend::tape::Scalar::atan2","Native normalized IEEE scalar operation. NaN/interval-entire sentinels retain their native meaning; no physical quantity units are erased.",true,(a:Raw,b:Raw)->Raw,|ctx,a,b|{let _=ctx;Ok(Raw(<f64 as Scalar>::atan2(a.0,b.0)))});
object_binding!(scalar_float_compare,"num_native_scalar_float_compare","crate::backend::tape::Scalar::compare","Native normalized IEEE scalar operation. NaN/interval-entire sentinels retain their native meaning; no physical quantity units are erased.",true,(a:Raw,comparison:Comparison,b:Raw)->Raw,|ctx,a,comparison,b|{let _=ctx;Ok(Raw(<f64 as Scalar>::compare(a.0,comparison.0,b.0)))});
object_binding!(scalar_float_select,"num_native_scalar_float_select","crate::backend::tape::Scalar::select","Native normalized IEEE scalar operation. NaN/interval-entire sentinels retain their native meaning; no physical quantity units are erased.",true,(condition:Raw,a:Raw,b:Raw)->Raw,|ctx,condition,a,b|{let _=ctx;Ok(Raw(<f64 as Scalar>::select(condition.0,a.0,b.0)))});
object_binding!(scalar_float_call,"num_native_scalar_float_call","crate::backend::tape::Scalar::call","Native normalized IEEE scalar operation. NaN/interval-entire sentinels retain their native meaning; no physical quantity units are erased.",true,(function:String,values:Input)->Raw,|ctx,function,values|{let _=ctx;api::normalized(&values.0)?;let Value::Tensor(t)=values.0 else{return Err(api::error("expected scalar argument vector"));};if t.shape!=[t.data.len()]||t.data.len()>4096{return Err(api::error("scalar argument vector shape/work budget"));}let args=t.data.iter().map(|v|<Raw as NativeArg>::decode(v,0,ctx).map(|r|r.0)).collect::<Result<Vec<_>,_>>()?;let f=eval_function(&function,args.len())?;Ok(Raw(<f64 as Scalar>::call(f,&args)))});
extern "C" fn external_exp(x: f64) -> f64 {
    Intrinsic::Exp.apply(x)
}
extern "C" fn external_ln(x: f64) -> f64 {
    Intrinsic::Ln.apply(x)
}
extern "C" fn external_sin(x: f64) -> f64 {
    Intrinsic::Sin.apply(x)
}
extern "C" fn external_cos(x: f64) -> f64 {
    Intrinsic::Cos.apply(x)
}
extern "C" fn external_tan(x: f64) -> f64 {
    Intrinsic::Tan.apply(x)
}
extern "C" fn external_asin(x: f64) -> f64 {
    Intrinsic::Asin.apply(x)
}
extern "C" fn external_acos(x: f64) -> f64 {
    Intrinsic::Acos.apply(x)
}
extern "C" fn external_atan(x: f64) -> f64 {
    Intrinsic::Atan.apply(x)
}
extern "C" fn external_sinh(x: f64) -> f64 {
    Intrinsic::Sinh.apply(x)
}
extern "C" fn external_cosh(x: f64) -> f64 {
    Intrinsic::Cosh.apply(x)
}
extern "C" fn external_tanh(x: f64) -> f64 {
    Intrinsic::Tanh.apply(x)
}
extern "C" fn external_asinh(x: f64) -> f64 {
    Intrinsic::Asinh.apply(x)
}
extern "C" fn external_acosh(x: f64) -> f64 {
    Intrinsic::Acosh.apply(x)
}
extern "C" fn external_atanh(x: f64) -> f64 {
    Intrinsic::Atanh.apply(x)
}
extern "C" fn external_sqrt(x: f64) -> f64 {
    Intrinsic::Sqrt.apply(x)
}
extern "C" fn external_abs(x: f64) -> f64 {
    Intrinsic::Abs.apply(x)
}
extern "C" fn external_floor(x: f64) -> f64 {
    Intrinsic::Floor.apply(x)
}
extern "C" fn external_ceil(x: f64) -> f64 {
    Intrinsic::Ceil.apply(x)
}
extern "C" fn external_cbrt(x: f64) -> f64 {
    Intrinsic::Cbrt.apply(x)
}
fn extern1(f: Intrinsic) -> extern "C" fn(f64) -> f64 {
    match f {
        | Intrinsic::Exp => external_exp,
        | Intrinsic::Ln => external_ln,
        | Intrinsic::Sin => external_sin,
        | Intrinsic::Cos => external_cos,
        | Intrinsic::Tan => external_tan,
        | Intrinsic::Asin => external_asin,
        | Intrinsic::Acos => external_acos,
        | Intrinsic::Atan => external_atan,
        | Intrinsic::Sinh => external_sinh,
        | Intrinsic::Cosh => external_cosh,
        | Intrinsic::Tanh => external_tanh,
        | Intrinsic::Asinh => external_asinh,
        | Intrinsic::Acosh => external_acosh,
        | Intrinsic::Atanh => external_atanh,
        | Intrinsic::Sqrt => external_sqrt,
        | Intrinsic::Abs => external_abs,
        | Intrinsic::Floor => external_floor,
        | Intrinsic::Ceil => external_ceil,
        | Intrinsic::Cbrt => external_cbrt,
    }
}
extern "C" fn external_add(
    a: f64,
    b: f64,
) -> f64 {
    <f64 as Scalar>::add(a, b)
}
object_binding!(scalar_from_extern1,"num_native_scalar_from_extern1","crate::backend::tape::FromExtern::from_extern1","Native normalized IEEE scalar operation. NaN/interval-entire sentinels retain their native meaning; no physical quantity units are erased.",true,(function:Function,x:Raw)->Raw,|ctx,function,x|{let _=ctx;Ok(Raw(<f64 as FromExtern>::from_extern1(extern1(function.0),x.0)))});
object_binding!(scalar_from_extern2,"num_native_scalar_from_extern2","crate::backend::tape::FromExtern::from_extern2","Native normalized IEEE scalar operation. NaN/interval-entire sentinels retain their native meaning; no physical quantity units are erased.",true,(a:Raw,b:Raw)->Raw,|ctx,a,b|{let _=ctx;Ok(Raw(<f64 as FromExtern>::from_extern2(external_add,a.0,b.0)))});
object_binding!(tape_int_pow,"num_native_tape_int_pow","crate::backend::tape::int_pow","Native normalized IEEE scalar operation. NaN/interval-entire sentinels retain their native meaning; no physical quantity units are erased.",true,(x:Raw,n:i32)->Raw,|ctx,x,n|{let _=ctx;Ok(Raw(tape::int_pow(x.0,n)))});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &hyperdual_new::BINDING,
            &hyperdual_value::BINDING,
            &hyperdual_first::BINDING,
            &hyperdual_second::BINDING,
            &hyperdual_mixed::BINDING,
            &interval_point::BINDING,
            &interval_new::BINDING,
            &interval_bounds::BINDING,
            &interval_contains::BINDING,
            &interval_is_point::BINDING,
            &intrinsic_apply::BINDING,
            &intrinsic_derivative::BINDING,
            &intrinsic_second_derivative::BINDING,
            &scalar_float_nan::BINDING,
            &scalar_float_constant::BINDING,
            &scalar_float_to_f64::BINDING,
            &scalar_float_add::BINDING,
            &scalar_float_sub::BINDING,
            &scalar_float_mul::BINDING,
            &scalar_float_div::BINDING,
            &scalar_float_neg::BINDING,
            &scalar_float_powi::BINDING,
            &scalar_float_powf::BINDING,
            &scalar_float_unary::BINDING,
            &scalar_float_atan2::BINDING,
            &scalar_float_compare::BINDING,
            &scalar_float_select::BINDING,
            &scalar_float_call::BINDING,
            &scalar_from_extern1::BINDING,
            &scalar_from_extern2::BINDING,
            &tape_int_pow::BINDING,
        ],
        "numerical_vector",
    )
}

pub(super) fn intrinsic(name: &str) -> Result<Intrinsic, OpError> {
    let f = <Function as NativeArg>::decode(&Value::String(name.into()), 0, &CallCtx::new())?;
    Ok(f.0)
}
