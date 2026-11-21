// Stateless normalized dual arithmetic; physical references belong to the AD driver.
impl Arg for crate::kernels::autodiff::Dual {
    fn ty() -> TypeExpr {
        TypeExpr::Record(vec![
            ("val".into(), TypeExpr::Float),
            ("deriv".into(), TypeExpr::Float),
        ])
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        Ok(Self::new(
            <f64 as Arg>::decode(model_field(value, "val", &["val", "deriv"])?, index)?,
            <f64 as Arg>::decode(model_field(value, "deriv", &["val", "deriv"])?, index)?,
        ))
    }
}
record!(crate::kernels::autodiff::Dual { val: f64, deriv: f64 });
binding!(dual_constant,"num_dual_constant","numerical_vector","A dimensionless normalized dual constant; used inside exact forward-AD callbacks.","crate::kernels::autodiff::Dual::constant",DimensionPolicy::Dimensionless,(x:f64)->crate::kernels::autodiff::Dual,|x:f64|crate::kernels::autodiff::Dual::constant(x));
binding!(dual_variable,"num_dual_variable","numerical_vector","A normalized dual variable with tangent one.","crate::kernels::autodiff::Dual::var",DimensionPolicy::Dimensionless,(x:f64)->crate::kernels::autodiff::Dual,|x:f64|crate::kernels::autodiff::Dual::var(x));
binding!(dual_new,"num_dual_new","numerical_vector","Construct a normalized dual from finite primal and tangent values.","crate::kernels::autodiff::Dual::new",DimensionPolicy::Dimensionless,(x:f64,tangent:f64)->crate::kernels::autodiff::Dual,|x:f64,tangent:f64|crate::kernels::autodiff::Dual::new(x,tangent));
// Trait operators have no public inherent-function source-inventory path.
macro_rules! dual_binary {($module:ident,$name:literal,$path:literal,$op:tt)=>{binding!($module,$name,"numerical_vector","Native normalized dual arithmetic; tangent propagation stays in the RSSN Dual implementation.",$path,DimensionPolicy::Dimensionless,(x:crate::kernels::autodiff::Dual,y:crate::kernels::autodiff::Dual)->crate::kernels::autodiff::Dual,|x:crate::kernels::autodiff::Dual,y:crate::kernels::autodiff::Dual|x $op y);};}
dual_binary!(dual_add,"num_dual_add","",+);
dual_binary!(dual_sub,"num_dual_sub","",-);
dual_binary!(dual_mul,"num_dual_mul","",*);
dual_binary!(dual_div,"num_dual_div","",/);
macro_rules! dual_unary {($module:ident,$name:literal,$path:literal,$method:ident)=>{binding!($module,$name,"numerical_vector","Native normalized dual elementary function; rejects nonfinite primal/tangent results.",$path,DimensionPolicy::Dimensionless,(x:crate::kernels::autodiff::Dual)->crate::kernels::autodiff::Dual,|x:crate::kernels::autodiff::Dual|x.$method());};}
dual_unary!(
    dual_sin,
    "num_dual_sin",
    "crate::kernels::autodiff::Dual::sin",
    sin
);
dual_unary!(
    dual_cos,
    "num_dual_cos",
    "crate::kernels::autodiff::Dual::cos",
    cos
);
dual_unary!(
    dual_tan,
    "num_dual_tan",
    "crate::kernels::autodiff::Dual::tan",
    tan
);
dual_unary!(
    dual_exp,
    "num_dual_exp",
    "crate::kernels::autodiff::Dual::exp",
    exp
);
dual_unary!(
    dual_ln,
    "num_dual_ln",
    "crate::kernels::autodiff::Dual::ln",
    ln
);
dual_unary!(
    dual_sqrt,
    "num_dual_sqrt",
    "crate::kernels::autodiff::Dual::sqrt",
    sqrt
);
dual_unary!(
    dual_sinh,
    "num_dual_sinh",
    "crate::kernels::autodiff::Dual::sinh",
    sinh
);
dual_unary!(
    dual_cosh,
    "num_dual_cosh",
    "crate::kernels::autodiff::Dual::cosh",
    cosh
);
dual_unary!(
    dual_tanh,
    "num_dual_tanh",
    "crate::kernels::autodiff::Dual::tanh",
    tanh
);
binding!(dual_pow,"num_dual_pow","numerical_vector","Native normalized dual floating power with a dimensionless exponent.","crate::kernels::autodiff::Dual::powf",DimensionPolicy::Dimensionless,(x:crate::kernels::autodiff::Dual,p:f64)->crate::kernels::autodiff::Dual,|x:crate::kernels::autodiff::Dual,p:f64|x.powf(p));
dual_unary!(
    dual_abs,
    "num_dual_abs",
    "crate::kernels::autodiff::Dual::abs",
    abs
);
dual_unary!(
    dual_acos,
    "num_dual_acos",
    "crate::kernels::autodiff::Dual::acos",
    acos
);
dual_unary!(
    dual_asin,
    "num_dual_asin",
    "crate::kernels::autodiff::Dual::asin",
    asin
);
dual_unary!(
    dual_atan,
    "num_dual_atan",
    "crate::kernels::autodiff::Dual::atan",
    atan
);
binding!(dual_powi,"num_dual_powi","numerical_vector","Native normalized dual integer power.","crate::kernels::autodiff::Dual::powi",DimensionPolicy::Dimensionless,(x:crate::kernels::autodiff::Dual,p:i32)->crate::kernels::autodiff::Dual,|x:crate::kernels::autodiff::Dual,p:i32|x.powi(p));

pub(super) struct DualIntegerPower {
    val: f64,
    deriv: f64,
    second: f64,
}
record!(DualIntegerPower { val: f64, deriv: f64, second: f64 });
binding!(dual_powi_with_second_derivative,"num_dual_powi_with_second_derivative","numerical_vector","Evaluate the normalized integer power along val + t * deriv at t = 0. Returns the primal, first derivative and second directional derivative. The exponent is an Int in the signed 32-bit range. Nonfinite results are rejected.","crate::kernels::autodiff::Dual::powi_with_second_derivative",DimensionPolicy::Dimensionless,(x:crate::kernels::autodiff::Dual,p:i32)->super::DualIntegerPower,|x:crate::kernels::autodiff::Dual,p:i32|{
    let (power, second) = x.powi_with_second_derivative(p);
    super::DualIntegerPower { val: power.val, deriv: power.deriv, second }
});

pub(super) static DUAL_EXPORTS: &[&CodecExport] = &[
    &dual_abs::EXPORT,
    &dual_acos::EXPORT,
    &dual_asin::EXPORT,
    &dual_atan::EXPORT,
    &dual_powi::EXPORT,
    &dual_powi_with_second_derivative::EXPORT,
    &dual_constant::EXPORT,
    &dual_variable::EXPORT,
    &dual_new::EXPORT,
    &dual_add::EXPORT,
    &dual_sub::EXPORT,
    &dual_mul::EXPORT,
    &dual_div::EXPORT,
    &dual_sin::EXPORT,
    &dual_cos::EXPORT,
    &dual_tan::EXPORT,
    &dual_exp::EXPORT,
    &dual_ln::EXPORT,
    &dual_sqrt::EXPORT,
    &dual_sinh::EXPORT,
    &dual_cosh::EXPORT,
    &dual_tanh::EXPORT,
    &dual_pow::EXPORT,
];
