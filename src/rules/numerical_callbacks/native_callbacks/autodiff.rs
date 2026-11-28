//! Native forward dual sweeps and scoped tape-based reverse accumulation.
use super::CallCtx;
use super::Callbacks;
use super::DimExpr;
use super::EvalDyn;
use super::Installer;
use super::OpError;
use super::RuleError;
use super::Signature;
use super::TypeExpr;
use super::Value;
use super::bounded_size;
use super::count;
use super::equal;
use super::finite;
use super::function;
use super::ratio;
use super::real;
use super::record;
use super::register;
use super::scalar;
use super::scale;
use super::tensor_type;
use super::typed;
use super::var;
use super::vec_value;
use super::vector;
use crate::kernel::value::TensorValue;
use crate::kernels::autodiff::Dual;
use crate::kernels::autodiff::{
    self,
};
use crate::rules::numerical_abi::native_runtime::ad_var_type;
use crate::rules::numerical_abi::native_runtime::decode_ad_var;
use crate::rules::numerical_abi::native_runtime::encode_ad_var;
#[derive(Clone, Copy)]
enum Method {
    Derivative,
    Gradient,
    Jacobian,
    Reverse,
}
fn dual_type(vector: bool) -> TypeExpr {
    let component = if vector {
        tensor_type(1, DimExpr::Dimless)
    } else {
        TypeExpr::Float
    };
    TypeExpr::Record(vec![
        ("val".into(), component.clone()),
        ("deriv".into(), component),
    ])
}
fn dual_value(v: Dual) -> Value {
    record(vec![
        ("val", Value::Float(v.val)),
        ("deriv", Value::Float(v.deriv)),
    ])
}
fn field<'a>(
    v: &'a Value,
    name: &str,
) -> Result<&'a Value, OpError> {
    let Value::Record(r) = v else {
        return Err(super::domain("dual callback requires val/deriv record"));
    };
    if r.fields.len() != 2 {
        return Err(super::domain("dual record requires exactly val and deriv"));
    }
    let mut fields = r.fields.iter().filter(|(n, _)| n == name);
    let out = &fields
        .next()
        .ok_or_else(|| super::domain("missing dual component"))?
        .1;
    if fields.next().is_some() {
        return Err(super::domain("duplicate dual component"));
    }
    Ok(out)
}
fn decode_dual(
    v: &Value,
    ctx: &CallCtx<'_>,
) -> Result<Dual, OpError> {
    Ok(Dual::new(
        scalar(field(v, "val")?, &DimExpr::Dimless, 4, ctx)?,
        scalar(field(v, "deriv")?, &DimExpr::Dimless, 4, ctx)?,
    ))
}
fn dual_vector(v: &[Dual]) -> Value {
    record(vec![
        (
            "val",
            vec_value(
                &v.iter().map(|v| v.val).collect::<Vec<_>>(),
                &DimExpr::Dimless,
            ),
        ),
        (
            "deriv",
            vec_value(
                &v.iter().map(|v| v.deriv).collect::<Vec<_>>(),
                &DimExpr::Dimless,
            ),
        ),
    ])
}
impl EvalDyn for Method {
    fn signature(&self) -> Signature {
        let d = var("D");
        let e = var("E");
        let first = matches!(self, Self::Derivative);
        let x = if first {
            real(d.clone())
        } else {
            tensor_type(1, d.clone())
        };
        let derivative = DimExpr::Mul(vec![e.clone(), DimExpr::Inverse(Box::new(d.clone()))]);
        let mut args = vec![x, real(d), real(e.clone()), TypeExpr::Int];
        let callback = if matches!(self, Self::Reverse) {
            function(
                vec![TypeExpr::Tensor {
                    rank: crate::kernel::signature::NatExpr::Lit(1),
                    shape: None,
                    elem: Box::new(ad_var_type(DimExpr::Dimless)),
                }],
                ad_var_type(DimExpr::Dimless),
            )
        } else {
            function(
                vec![dual_type(!first)],
                dual_type(matches!(self, Self::Jacobian)),
            )
        };
        args.push(callback);
        let result = match self {
            | Self::Gradient => tensor_type(1, derivative),
            | Self::Jacobian => tensor_type(2, derivative),
            | Self::Derivative => {
                TypeExpr::Record(vec![
                    ("value".into(), real(e)),
                    ("derivative".into(), real(derivative)),
                    ("evaluations".into(), TypeExpr::Int),
                ])
            },
            | Self::Reverse => {
                TypeExpr::Record(vec![
                    ("value".into(), real(e)),
                    ("gradient".into(), tensor_type(1, derivative)),
                    ("evaluations".into(), TypeExpr::Int),
                ])
            },
        };
        typed(args, result, &["D", "E"])
    }

    // Keep each native method's validation, callback conversion and result units together.
    #[allow(clippy::too_many_lines)]
    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let expected = self.signature();
        super::arity(args, 5)?;
        let callbacks = Callbacks::new(ctx, args, &expected, 4, count(&args[3], 3)?)?;
        let d = if matches!(self, Self::Derivative) {
            super::super::scalar_dimension(&ctx.value_type(&args[0])?)?
        } else {
            let TypeExpr::Tensor { elem, .. } = ctx.value_type(&args[0])? else {
                return Err(super::domain("AD point must be a vector"));
            };
            super::super::scalar_dimension(&elem)?
        };
        let e = super::super::scalar_dimension(&ctx.value_type(&args[2])?)?;
        let xs = scale(&args[1], &d, 1, ctx)?;
        let ys = scale(&args[2], &e, 2, ctx)?;
        let x = if matches!(self, Self::Derivative) {
            vec![scalar(&args[0], &d, 0, ctx)? / xs]
        } else {
            vector(&args[0], &d, None, 0, ctx)?
                .into_iter()
                .map(|v| v / xs)
                .collect()
        };
        if x.is_empty() || x.iter().any(|x| !x.is_finite()) {
            return Err(super::domain(
                "AD requires a nonempty finite normalized point",
            ));
        }
        bounded_size(x.len(), x.len())?;
        let derivative = ratio(&e, &d)?;
        match self {
            | Self::Derivative => {
                let out = autodiff::forward_derivative(
                    &mut |v| {
                        callbacks
                            .call(0, &[dual_value(v)], |v| decode_dual(&v, ctx))
                            .unwrap_or(Dual::new(f64::NAN, f64::NAN))
                    },
                    x[0],
                );
                let used = callbacks.used.get();
                callbacks.finish()?;
                Ok(record(vec![
                    (
                        "value",
                        super::super::dimensional_float(finite(out.0 * ys)?, e),
                    ),
                    (
                        "derivative",
                        super::super::dimensional_float(finite(out.1 * (ys / xs))?, derivative),
                    ),
                    ("evaluations", Value::int(used)),
                ]))
            },
            | Self::Gradient => {
                let out = autodiff::forward_grad(
                    &mut |v| {
                        callbacks
                            .call(0, &[dual_vector(v)], |v| decode_dual(&v, ctx))
                            .unwrap_or(Dual::new(f64::NAN, f64::NAN))
                    },
                    &x,
                );
                callbacks.finish()?;
                let out = out.map_err(|e| super::domain(e.to_string()))?;
                let out = out
                    .into_iter()
                    .map(|v| finite(v * (ys / xs)))
                    .collect::<Result<Vec<_>, _>>()?;
                super::output_tensor(vec![x.len()], out, derivative)
            },
            | Self::Jacobian => {
                let out = autodiff::forward_jacobian(
                    &mut |v| {
                        callbacks
                            .call(0, &[dual_vector(v)], |v| {
                                let values =
                                    vector(field(&v, "val")?, &DimExpr::Dimless, None, 4, ctx)?;
                                let tangent = vector(
                                    field(&v, "deriv")?,
                                    &DimExpr::Dimless,
                                    Some(values.len()),
                                    4,
                                    ctx,
                                )?;
                                bounded_size(x.len(), values.len())?;
                                Ok(values
                                    .into_iter()
                                    .zip(tangent)
                                    .map(|(v, d)| Dual::new(v, d))
                                    .collect())
                            })
                            .unwrap_or_else(|| vec![Dual::new(f64::NAN, f64::NAN); x.len()])
                    },
                    &x,
                );
                callbacks.finish()?;
                let out = out.map_err(|e| super::domain(e.to_string()))?;
                let rows = out.len() / x.len();
                let out = out
                    .into_iter()
                    .map(|v| finite(v * (ys / xs)))
                    .collect::<Result<Vec<_>, _>>()?;
                super::output_tensor(vec![rows, x.len()], out, derivative)
            },
            | Self::Reverse => {
                let bank = ctx
                    .native_objects
                    .as_ref()
                    .ok_or_else(|| super::domain("reverse AD needs a session-owned object bank"))?;
                let _scope = bank.scoped()?;
                let _charge = bank.insert(
                    "rssn/native/ad/callback-charge/v1",
                    vec![],
                    (),
                    65_536 * 128,
                )?;
                let transport_error = std::cell::RefCell::new(None);
                let out = autodiff::reverse_grad(
                    &mut |vars| {
                        let input = vars
                            .iter()
                            .map(|v| encode_ad_var(ctx, v.clone(), DimExpr::Dimless))
                            .collect::<Result<Vec<_>, _>>();
                        let input = match input {
                            | Ok(input) => {
                                Value::Tensor(TensorValue {
                                    shape: vec![vars.len()],
                                    data: input,
                                    elem: ad_var_type(DimExpr::Dimless),
                                })
                            },
                            | Err(error) => {
                                *transport_error.borrow_mut() = Some(error);
                                return vars[0].clone();
                            },
                        };
                        callbacks
                            .call(0, &[input], |value| {
                                let (v, d) = decode_ad_var(ctx, &value)?;
                                equal(&d, &DimExpr::Dimless)?;
                                Ok(v)
                            })
                            .unwrap_or_else(|| vars[0].clone())
                    },
                    &x,
                );
                let used = callbacks.used.get();
                callbacks.finish()?;
                if let Some(error) = transport_error.into_inner() {
                    return Err(error);
                }
                let (value, gradient) = out.map_err(|e| super::domain(e.to_string()))?;
                let gradient = gradient
                    .into_iter()
                    .map(|v| finite(v * (ys / xs)))
                    .collect::<Result<Vec<_>, _>>()?;
                Ok(record(vec![
                    (
                        "value",
                        super::super::dimensional_float(finite(value * ys)?, e),
                    ),
                    (
                        "gradient",
                        super::output_tensor(vec![x.len()], gradient, derivative)?,
                    ),
                    ("evaluations", Value::int(used)),
                ]))
            },
        }
    }
}
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    for (name, path, m) in [
        (
            "num_dual_forward_derivative_scaled",
            "forward_derivative",
            Method::Derivative,
        ),
        (
            "num_dual_forward_gradient_scaled",
            "forward_grad",
            Method::Gradient,
        ),
        (
            "num_dual_forward_jacobian_scaled",
            "forward_jacobian",
            Method::Jacobian,
        ),
        (
            "num_native_reverse_gradient_scaled",
            "reverse_grad",
            Method::Reverse,
        ),
    ] {
        register(
            i,
            name,
            &format!("crate::kernels::autodiff::{path}"),
            "Exact native AD with normalized callback variables and independent physical references: point,positive coordinate_scale,positive value_scale,max_evaluations,callback. Forward callback scalar duals are {val,deriv}; vector duals are {val:vector,deriv:vector}; all dual arithmetic is dimensionless and uses native num_dual_* operators. Reverse callbacks use checked NativeAdVar handles on an ephemeral session tape; native operations enforce same-tape identity and a65536-node limit. Output derivatives restore E/D; no numerical differencing or opaque-callback finite-difference fallback.",
            5,
            m,
        )?;
    }
    Ok(())
}
