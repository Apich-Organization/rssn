//! Checked live callbacks for portable native Chebyshev approximation models.
use super::CallCtx;
use super::DimExpr;
use super::EvalDyn;
use super::Installer;
use super::OpError;
use super::RuleError;
use super::Signature;
use super::TypeExpr;
use super::Value;
use super::common::CallbackRun;
use super::common::{
    self,
};
use crate::kernel::signature::ArityInfo;
use crate::kernel::signature::TypeParam;
use crate::kernel::signature::TypeParamKind;
use crate::kernel::value::DynValue;
use crate::kernels::approx::ChebSeries;
use crate::rules::numerical_abi::native_models::cheb_type;
use crate::rules::numerical_abi::native_models::encode_cheb;
use std::cell::Cell;

struct Construct {
    adaptive: bool,
}
impl EvalDyn for Construct {
    fn signature(&self) -> Signature {
        let x = DimExpr::Var("X".into());
        let y = DimExpr::Var("Y".into());
        let scalar = |dim| TypeExpr::Float.with_dim(dim);
        let mut ret = cheb_type();
        if let TypeExpr::Record(fields) = &mut ret {
            for (name, ty) in fields {
                match name.as_str() {
                    | "coefficients" => {
                        if let TypeExpr::Tensor { elem, .. } = ty {
                            **elem = scalar(y.clone());
                        }
                    },
                    | "a" | "b" => *ty = scalar(x.clone()),
                    | _ => {},
                }
            }
        }
        let mut params = vec![scalar(x.clone()), scalar(x.clone())];
        if self.adaptive {
            params.extend([TypeExpr::Float, TypeExpr::Int]);
        } else {
            params.push(TypeExpr::Int);
        }
        params.push(TypeExpr::Int);
        params.push(TypeExpr::Function(vec![scalar(x)], Box::new(scalar(y))));
        let mut signature = common::effectful(Signature::new(params, ret));
        signature.type_params = ["X", "Y"]
            .into_iter()
            .map(|name| {
                TypeParam {
                    name: name.into(),
                    kind: TypeParamKind::Dim,
                    bounds: vec![],
                }
            })
            .collect();
        signature
    }

    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let signature = self.signature();
        let count = signature.params.len();
        if args.len() != count {
            return Err(OpError::Arity {
                expected: ArityInfo::fixed(count as u8),
                got: args.len(),
            });
        }
        let ci = count - 1;
        let (run, schema) = CallbackRun::new(ctx, &args[ci], ci, 1, signature.params[ci].clone())?;
        let x = super::scalar_dimension(&schema.params[0])?;
        let y = super::scalar_dimension(&schema.ret)?;
        super::check_dimension(ctx, &args[0], &x, 0)?;
        super::check_dimension(ctx, &args[1], &x, 1)?;
        let a = super::numeric(&args[0], 0)?;
        let b = super::numeric(&args[1], 1)?;
        if a >= b || !(b - a).is_finite() {
            return Err(OpError::Domain("Chebyshev interval must increase".into()));
        }
        let degree_index = if self.adaptive { 3 } else { 2 };
        let degree = common::positive_budget(&args[degree_index], degree_index, "degree")?;
        if degree > 64 || self.adaptive && (degree < 8 || !degree.is_power_of_two()) {
            return Err(OpError::Domain("degree must be 1..=64; adaptive max_degree must be a power of two from 8 through 64".into()));
        }
        let budget = common::positive_budget(&args[ci - 1], ci - 1, "callback evaluation budget")?;
        let exhausted = Cell::new(false);
        let callback = |coordinate: f64| {
            if run.failed() {
                return f64::NAN;
            }
            if run.evaluations() >= budget {
                exhausted.set(true);
                return f64::NAN;
            }
            run.evaluate(
                &[Value::Dyn(DynValue {
                    dim: x.clone(),
                    value: Box::new(Value::Float(coordinate)),
                })],
                |value| {
                    super::check_dimension(ctx, &value, &y, 0)?;
                    super::numeric(&value, 0)
                },
            )
            .unwrap_or(f64::NAN)
        };
        let series = if self.adaptive {
            super::check_dimension(ctx, &args[2], &DimExpr::Dimless, 2)?;
            let tolerance = super::numeric(&args[2], 2)?;
            if tolerance <= 0.0 {
                return Err(OpError::Domain(
                    "adaptive approximation tolerance must be positive".into(),
                ));
            }
            ChebSeries::adaptive(callback, a, b, tolerance, degree)
        } else {
            ChebSeries::from_fn(callback, a, b, degree)
        };
        run.finish()?;
        if exhausted.get() {
            return Err(OpError::Domain(
                "Chebyshev callback evaluation budget exhausted".into(),
            ));
        }
        let model = encode_cheb(series)?;
        let Value::Record(mut record) = model else {
            unreachable!()
        };
        for (name, value) in &mut record.fields {
            match name.as_str() {
                | "coefficients" => {
                    if let Value::Tensor(tensor) = value {
                        tensor.elem = TypeExpr::Float.with_dim(y.clone());
                    }
                },
                | "a" | "b" => {
                    let old = std::mem::replace(value, Value::Bool(false));
                    *value = Value::Dyn(DynValue {
                        dim: x.clone(),
                        value: Box::new(old),
                    });
                },
                | _ => {},
            }
        }
        Ok(Value::Record(record))
    }
}
pub(super) fn install(installer: &mut Installer<'_>) -> Result<(), RuleError> {
    for (adaptive, name, path) in [
        (
            false,
            "num_native_cheb_from_fn",
            "crate::kernels::approx::ChebSeries::from_fn",
        ),
        (
            true,
            "num_native_cheb_adaptive",
            "crate::kernels::approx::ChebSeries::adaptive",
        ),
    ] {
        let evaluator = Construct { adaptive };
        let arity = evaluator.signature().params.len() as u8;
        common::register(
            installer,
            name,
            path,
            "Construct a portable native Chebyshev series from a checked serial callback with independent physical coordinate and ordinate units. Bounds, degree, finite callback values, original callback failures and evaluation budgets are checked. Adaptive max_degree is a power of two from 8 through 64. Returns the same coefficient model accepted by native evaluation, derivative, integral and root methods.",
            arity,
            evaluator,
        )?;
    }
    Ok(())
}
