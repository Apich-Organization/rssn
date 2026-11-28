//! Native sequence and series algorithms on normalized physical values.
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
use super::count;
use super::finite;
use super::function;
use super::real;
use super::record;
use super::register;
use super::scalar;
use super::scalar_dimension;
use super::scale;
use super::typed;
use super::var;
use crate::kernels::convergence::LevinVariant;
use crate::kernels::convergence::{
    self,
};
use num_traits::ToPrimitive;
#[derive(Clone, Copy)]
enum Method {
    Series,
    Limit,
    LevinU,
    LevinT,
    LevinV,
}
impl Method {
    const fn levin(self) -> bool {
        matches!(self, Self::LevinU | Self::LevinT | Self::LevinV)
    }
}
impl EvalDyn for Method {
    fn signature(&self) -> Signature {
        let e = var("E");
        let mut params = if self.levin() {
            vec![TypeExpr::Int, real(e.clone()), TypeExpr::Int]
        } else {
            vec![
                TypeExpr::Int,
                real(e.clone()),
                real(e.clone()),
                TypeExpr::Int,
            ]
        };
        if matches!(self, Self::Series) {
            params.insert(0, TypeExpr::Int);
        }
        params.push(function(
            vec![if self.levin() {
                TypeExpr::Int
            } else {
                TypeExpr::Float
            }],
            real(e.clone()),
        ));
        let mut fields = vec![
            ("value".into(), real(e.clone())),
            ("evaluations".into(), TypeExpr::Int),
        ];
        if self.levin() {
            fields.push(("error".into(), real(e)));
        }
        typed(params, TypeExpr::Record(fields), &["E"])
    }

    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let sig = self.signature();
        super::arity(args, sig.params.len())?;
        let offset = usize::from(matches!(self, Self::Series));
        let maximum = if self.levin() {
            2
        } else {
            3 + offset
        };
        let callbacks = Callbacks::new(
            ctx,
            args,
            &sig,
            maximum + 1,
            count(&args[maximum], maximum)?,
        )?;
        let e = scalar_dimension(&callbacks.signatures[0].ret)?;
        super::equal(
            &scalar_dimension(&callbacks.signatures[0].params[0])?,
            &DimExpr::Dimless,
        )?;
        let terms = count(&args[offset], offset)?;
        let scale_at = if self.levin() {
            1
        } else {
            2 + offset
        };
        let ys = scale(&args[scale_at], &e, scale_at, ctx)?;
        let term = |n: f64| {
            callbacks
                .call(0, &[Value::Float(n)], |v| {
                    scalar(&v, &e, maximum + 1, ctx).map(|v| v / ys)
                })
                .unwrap_or(f64::NAN)
        };
        let value = if self.levin() {
            let variant = match self {
                | Self::LevinU => LevinVariant::U,
                | Self::LevinT => LevinVariant::T,
                | _ => LevinVariant::V,
            };
            let pair = convergence::levin_sum(
                |n| {
                    callbacks
                        .call(0, &[Value::int(n)], |v| {
                            scalar(&v, &e, maximum + 1, ctx).map(|v| v / ys)
                        })
                        .unwrap_or(f64::NAN)
                },
                terms,
                variant,
            );
            let used = callbacks.used.get();
            callbacks.finish()?;
            let (v, err) = pair.map_err(super::domain)?;
            return Ok(record(vec![
                (
                    "value",
                    super::super::dimensional_float(finite(v * ys)?, e.clone()),
                ),
                (
                    "error",
                    super::super::dimensional_float(finite(err * ys)?, e),
                ),
                ("evaluations", Value::int(used)),
            ]));
        } else {
            let tolerance = scale(&args[1 + offset], &e, 1 + offset, ctx)? / ys;
            if matches!(self, Self::Series) {
                let Value::Int(n) = &args[0] else {
                    return Err(super::domain("series start must be a nonnegative integer"));
                };
                let start = n
                    .to_usize()
                    .ok_or_else(|| super::domain("series start must fit usize"))?;
                Ok(convergence::sum_series_numerical(
                    term, start, terms, tolerance,
                ))
            } else {
                convergence::find_sequence_limit(term, terms, tolerance)
            }
        };
        let used = callbacks.used.get();
        callbacks.finish()?;
        let value = value.map_err(super::domain)?;
        Ok(record(vec![
            (
                "value",
                super::super::dimensional_float(finite(value * ys)?, e),
            ),
            ("evaluations", Value::int(used)),
        ]))
    }
}
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    for (name, path, m) in [
        (
            "num_native_sum_series",
            "sum_series_numerical",
            Method::Series,
        ),
        (
            "num_native_sequence_limit",
            "find_sequence_limit",
            Method::Limit,
        ),
        ("num_native_levin_u", "levin_sum", Method::LevinU),
        ("num_native_levin_t", "levin_sum", Method::LevinT),
        ("num_native_levin_v", "levin_sum", Method::LevinV),
    ] {
        register(
            i,
            name,
            &format!("crate::kernels::convergence::{path}"),
            "Physical native series/sequence convergence. sum_series: start,max_terms,tolerance,value_scale,max_evaluations,term(Float)->E. sequence_limit omits start. Levin: max_terms,value_scale,max_evaluations,term(Int)->E; result includes native error estimate. References/tolerances are positive physical E; callback failures propagate and budgets bound actual host calls.",
            if m.levin() {
                4
            } else if matches!(m, Method::Series) {
                6
            } else {
                5
            },
            m,
        )?;
    }
    Ok(())
}
