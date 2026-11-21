//! Explicit reference-scale adapters for measured nonlinear problems.
//!
//! Native solvers operate on nondimensional coordinates and residuals. This
//! makes mixed residual/step stopping criteria meaningful without assigning a
//! physical unit to a numerical control tolerance.
use super::CallCtx;
use super::DimExpr;
use super::EvalDyn;
use super::Installer;
use super::OpError;
use super::RecordValue;
use super::RuleError;
use super::Signature;
use super::TypeExpr;
use super::TypeParam;
use super::TypeParamKind;
use super::Value;
use super::check_dimension;
use super::common::CallbackRun;
use super::common::effectful;
use super::common::positive_budget;
use super::common::register;
use super::dimensional_float;
use super::dimensional_type;
use super::evolution::output_tensor;
use super::evolution::tensor_type;
use super::evolution::vector;
use super::evolution::vector_dimension;
use super::numeric;
use super::scalar_dimension;
use crate::kernels::optim;
use crate::kernels::rootfind;
use std::cell::Cell;
use std::cell::RefCell;

#[derive(Clone, Copy)]
enum Method {
    Bfgs,
    Lbfgs,
    Cg,
    Box,
    Trust,
    Broyden,
    Dogleg,
}
impl Method {
    const fn system(self) -> bool {
        matches!(self, Self::Broyden | Self::Dogleg)
    }

    const fn bounded(self) -> bool {
        matches!(self, Self::Box)
    }
}
fn control_type(
    d: DimExpr,
    e: DimExpr,
) -> TypeExpr {
    TypeExpr::Record(vec![
        ("coordinate_scale".into(), dimensional_type(d)),
        ("value_scale".into(), dimensional_type(e)),
        ("tolerance".into(), TypeExpr::Float),
        ("max_iterations".into(), TypeExpr::Int),
        ("max_evaluations".into(), TypeExpr::Int),
    ])
}
fn field<'a>(
    value: &'a Value,
    name: &str,
) -> Result<&'a Value, OpError> {
    let Value::Record(r) = value else {
        return Err(OpError::Domain("expected nonlinear control record".into()));
    };
    let matches: Vec<_> = r.fields.iter().filter(|(key, _)| key == name).collect();
    if matches.len() != 1 {
        return Err(OpError::Domain(format!(
            "missing or duplicate nonlinear control field {name}"
        )));
    }
    Ok(&matches[0].1)
}
fn function(
    input: TypeExpr,
    output: TypeExpr,
) -> TypeExpr {
    TypeExpr::Function(vec![input], Box::new(output))
}
fn equal(
    a: &DimExpr,
    b: &DimExpr,
) -> Result<(), OpError> {
    if super::canonical(a)? == super::canonical(b)? {
        Ok(())
    } else {
        Err(OpError::Domain(
            "nonlinear callback has incompatible dimensions".into(),
        ))
    }
}
fn normalized(
    value: f64,
    numerator: f64,
    denominator: f64,
) -> Result<f64, OpError> {
    if crate::kernels::exact_eq(value, 0.0) {
        return Ok(0.0);
    }
    let product = value * numerator;
    let mut result = product / denominator;
    if !result.is_finite() || crate::kernels::exact_eq(product, 0.0) {
        result = (value / denominator) * numerator;
        if !result.is_finite() || crate::kernels::exact_eq(result, 0.0) {
            result = (value.abs().ln() + numerator.ln() - denominator.ln())
                .exp()
                .copysign(value);
        }
    }
    if result.is_finite() {
        Ok(result)
    } else {
        Err(OpError::Domain(
            "reference scales produce a non-finite normalized coordinate/value".into(),
        ))
    }
}

impl EvalDyn for Method {
    fn signature(&self) -> Signature {
        let d = DimExpr::Var("D".into());
        let e = DimExpr::Var("E".into());
        let grad = DimExpr::Mul(vec![e.clone(), DimExpr::Inverse(Box::new(d.clone()))]);
        let x = tensor_type(1, d.clone());
        let y = tensor_type(1, e.clone());
        let mut args = vec![x.clone(), control_type(d, e.clone())];
        if self.bounded() {
            args.extend([x.clone(), x.clone()]);
        }
        args.push(function(
            x.clone(),
            if self.system() {
                y.clone()
            } else {
                dimensional_type(e.clone())
            },
        ));
        if !self.system() {
            args.push(function(x.clone(), tensor_type(1, grad)));
        }
        let mut fields = vec![
            ("x".into(), x),
            ("iterations".into(), TypeExpr::Int),
            ("evaluations".into(), TypeExpr::Int),
            ("converged".into(), TypeExpr::Bool),
        ];
        fields.push(if self.system() {
            ("residual".into(), y)
        } else {
            ("value".into(), dimensional_type(e))
        });
        let mut sig = effectful(Signature::new(args, TypeExpr::Record(fields)));
        sig.type_params = ["D", "E"]
            .into_iter()
            .map(|name| {
                TypeParam {
                    name: name.into(),
                    kind: TypeParamKind::Dim,
                    bounds: vec![],
                }
            })
            .collect();
        sig
    }

    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let expected = self.signature();
        if args.len() != expected.params.len() {
            return Err(OpError::Arity {
                expected: crate::kernel::signature::ArityInfo::fixed(expected.params.len() as u8),
                got: args.len(),
            });
        }
        let callback_index = if self.bounded() { 4 } else { 2 };
        let (f, fsig) = CallbackRun::new(
            ctx,
            &args[callback_index],
            callback_index,
            1,
            expected.params[callback_index].clone(),
        )?;
        let d = vector_dimension(&fsig.params[0], None)?;
        let e = if self.system() {
            vector_dimension(&fsig.ret, None)?
        } else {
            scalar_dimension(&fsig.ret)?
        };
        let gd = super::expression(&super::canonical(&DimExpr::Mul(vec![
            e.clone(),
            DimExpr::Inverse(Box::new(d.clone())),
        ]))?)?;
        let initial = vector(&args[0], &d, None, 0, ctx)?;
        let n = initial.len();
        if n == 0 {
            return Err(OpError::Domain("initial vector must be nonempty".into()));
        }
        let storage = if self.system() || matches!(self, Self::Bfgs | Self::Trust) {
            n.checked_mul(n).and_then(|v| v.checked_mul(4))
        } else {
            n.checked_mul(32)
        };
        if storage.is_none_or(|n| n > 8_388_608) {
            return Err(OpError::Domain(
                "nonlinear solver work arrays exceed the finite allocation budget".into(),
            ));
        }
        vector_dimension(&fsig.params[0], Some(n))?;
        if self.system() {
            vector_dimension(&fsig.ret, Some(n))?;
        }
        let controls = &args[1];
        let xval = field(controls, "coordinate_scale")?;
        let yval = field(controls, "value_scale")?;
        check_dimension(ctx, xval, &d, 1)?;
        check_dimension(ctx, yval, &e, 1)?;
        let xs = numeric(xval, 1)?;
        let ys = numeric(yval, 1)?;
        let tolval = field(controls, "tolerance")?;
        check_dimension(ctx, tolval, &DimExpr::Dimless, 1)?;
        let tol = numeric(tolval, 1)?;
        if xs <= 0.0 || ys <= 0.0 || tol <= 0.0 {
            return Err(OpError::Domain(
                "reference scales and tolerance must be positive".into(),
            ));
        }
        let max_iter = positive_budget(
            field(controls, "max_iterations")?,
            1,
            "nonlinear iteration budget",
        )?;
        let max_eval = positive_budget(
            field(controls, "max_evaluations")?,
            1,
            "nonlinear callback budget",
        )?;
        let used = Cell::new(0usize);
        let exhausted = Cell::new(false);
        let failed = Cell::new(false);
        let conversion_error = RefCell::new(None);
        let evaluate = |cb: &CallbackRun<'_, '_>,
                        z: &[f64],
                        decode: &dyn Fn(Value) -> Result<Vec<f64>, OpError>|
         -> Option<Vec<f64>> {
            if failed.get() {
                return None;
            }
            if used.get() >= max_eval {
                exhausted.set(true);
                failed.set(true);
                return None;
            }
            used.set(used.get() + 1);
            let x: Vec<_> = z.iter().map(|v| v * xs).collect();
            let input = match output_tensor(vec![n], x, d.clone()) {
                | Ok(v) => v,
                | Err(error) => {
                    *conversion_error.borrow_mut() = Some(error);
                    failed.set(true);
                    return None;
                },
            };
            let result = cb.evaluate(&[input], decode);
            if result.is_none() {
                failed.set(true);
            }
            result
        };
        let z: Vec<_> = initial
            .iter()
            .map(|x| normalized(*x, 1.0, xs))
            .collect::<Result<_, _>>()?;
        let result = if self.system() {
            let residual = |z: &[f64]| {
                evaluate(&f, z, &|value| {
                    vector(&value, &e, Some(n), callback_index, ctx)?
                        .into_iter()
                        .map(|v| normalized(v, 1.0, ys))
                        .collect()
                })
                .unwrap_or_else(|| vec![f64::NAN; n])
            };
            let solution = if matches!(self, Self::Broyden) {
                rootfind::broyden(residual, &z, tol, max_iter)
            } else {
                rootfind::dogleg(residual, &z, tol, max_iter)
            };
            let solution = solution
                .map_err(|err| OpError::NotReduced(format!("nonlinear system failed: {err:?}")));
            solution.map(|solution| {
                let residual = residual(&solution.x);
                (solution.x, solution.iterations, true, None, Some(residual))
            })
        } else {
            let (g, gsig) = CallbackRun::new(
                ctx,
                &args[callback_index + 1],
                callback_index + 1,
                1,
                expected.params[callback_index + 1].clone(),
            )?;
            equal(&vector_dimension(&gsig.params[0], Some(n))?, &d)?;
            equal(&vector_dimension(&gsig.ret, Some(n))?, &gd)?;
            let objective = |z: &[f64]| {
                evaluate(&f, z, &|value| {
                    check_dimension(ctx, &value, &e, callback_index)?;
                    Ok(vec![normalized(numeric(&value, callback_index)?, 1.0, ys)?])
                })
                .map_or(f64::NAN, |v| v[0])
            };
            let gradient = |z: &[f64]| {
                evaluate(&g, z, &|value| {
                    vector(&value, &gd, Some(n), callback_index + 1, ctx)?
                        .into_iter()
                        .map(|v| normalized(v, xs, ys))
                        .collect()
                })
                .unwrap_or_else(|| vec![f64::NAN; n])
            };
            let solution = match self {
                | Self::Bfgs => optim::bfgs(objective, Some(&gradient), &z, tol, max_iter),
                | Self::Lbfgs => optim::lbfgs(objective, Some(&gradient), &z, 10, tol, max_iter),
                | Self::Cg => optim::nonlinear_cg(objective, Some(&gradient), &z, tol, max_iter),
                | Self::Trust => {
                    optim::trust_region(
                        objective,
                        gradient,
                        None::<&fn(&[f64]) -> crate::kernels::dense::Mat>,
                        &z,
                        tol,
                        max_iter,
                    )
                },
                | Self::Box => {
                    let lower = vector(&args[2], &d, Some(n), 2, ctx)?;
                    let upper = vector(&args[3], &d, Some(n), 3, ctx)?;
                    if lower.iter().zip(&upper).any(|(a, b)| a > b) {
                        return Err(OpError::Domain("lower bound exceeds upper bound".into()));
                    }
                    let lo: Vec<_> = lower
                        .iter()
                        .map(|v| normalized(*v, 1.0, xs))
                        .collect::<Result<_, _>>()?;
                    let hi: Vec<_> = upper
                        .iter()
                        .map(|v| normalized(*v, 1.0, xs))
                        .collect::<Result<_, _>>()?;
                    optim::lbfgs_b(objective, Some(&gradient), &z, &lo, &hi, 10, tol, max_iter)
                },
                | Self::Broyden | Self::Dogleg => unreachable!(),
            };
            g.finish()?;
            Ok((
                solution.x,
                solution.iterations,
                solution.converged,
                Some(solution.fx),
                None,
            ))
        };
        f.finish()?;
        if let Some(error) = conversion_error.into_inner() {
            return Err(error);
        }
        if exhausted.get() {
            return Err(OpError::NotReduced(
                "nonlinear callback budget exhausted".into(),
            ));
        }
        let (x, iterations, converged, value, residual) = result?;
        let x = x.into_iter().map(|v| v * xs).collect();
        let mut fields = vec![
            ("x".into(), output_tensor(vec![n], x, d)?),
            ("iterations".into(), Value::int(iterations)),
            ("evaluations".into(), Value::int(used.get())),
            ("converged".into(), Value::Bool(converged)),
        ];
        if let Some(value) = value {
            let value = value * ys;
            if !value.is_finite() {
                return Err(OpError::Domain("nonfinite nonlinear objective".into()));
            }
            fields.push(("value".into(), dimensional_float(value, e)));
        } else if let Some(residual) = residual {
            fields.push((
                "residual".into(),
                output_tensor(vec![n], residual.into_iter().map(|v| v * ys).collect(), e)?,
            ));
        }
        Ok(Value::Record(RecordValue { fields }))
    }
}
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    for (name, path, method) in [
        (
            "num_bfgs_scaled",
            "crate::kernels::optim::bfgs",
            Method::Bfgs,
        ),
        (
            "num_lbfgs_scaled",
            "crate::kernels::optim::lbfgs",
            Method::Lbfgs,
        ),
        (
            "num_nonlinear_cg_scaled",
            "crate::kernels::optim::nonlinear_cg",
            Method::Cg,
        ),
        (
            "num_lbfgs_b_scaled",
            "crate::kernels::optim::lbfgs_b",
            Method::Box,
        ),
        (
            "num_trust_region_scaled",
            "crate::kernels::optim::trust_region",
            Method::Trust,
        ),
        (
            "num_broyden_scaled",
            "crate::kernels::rootfind::broyden",
            Method::Broyden,
        ),
        (
            "num_dogleg_scaled",
            "crate::kernels::rootfind::dogleg",
            Method::Dogleg,
        ),
    ] {
        let arity = method.signature().params.len() as u8;
        register(
            i,
            name,
            path,
            "Measured nonlinear problem with explicit coordinate_scale and value_scale, dimensionless tolerance, max_iterations and max_evaluations controls. Solves in normalized coordinates/residuals. Arguments: initial, controls, optional lower/upper, objective/residual, optional gradient. Gradient units are objective/coordinate. Homogeneous coordinate and residual vectors; callback failures propagate and budgets are shared.",
            arity,
            method,
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod normalization_tests {
    use super::normalized;
    #[test]
    fn reference_normalization_avoids_intermediate_overflow_and_underflow() {
        for value in [1e-300_f64, 1e300] {
            let scaled = normalized(value, value.abs(), value.abs()).unwrap();
            assert!((scaled / value - 1.0).abs() < 1e-12);
        }
        assert!(normalized(f64::MAX, 1.0, f64::MIN_POSITIVE).is_err());
    }
}
