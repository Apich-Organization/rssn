//! Finite-difference native operators with explicit coordinate/output references.
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
use super::register;
use super::scalar;
use super::scalar_dimension;
use super::scale;
use super::scaled_vec_value;
use super::tensor_type;
use super::typed;
use super::var;
use super::vector;
use super::vector_dimension;
use crate::kernels::calculus;
use num_traits::ToPrimitive;
#[derive(Clone, Copy)]
enum Method {
    Derivative,
    Richardson,
    Second,
    Partial,
    Gradient,
    Hessian,
    HessianDiagonal,
    HessianDiagonalWithSteps,
    Jacobian,
    Divergence,
    Curl,
    Laplacian,
}
impl Method {
    const fn scalar_input(self) -> bool {
        matches!(self, Self::Derivative | Self::Richardson | Self::Second)
    }

    const fn vector_output(self) -> bool {
        matches!(self, Self::Jacobian | Self::Divergence | Self::Curl)
    }

    const fn second(self) -> bool {
        matches!(
            self,
            Self::Second
                | Self::Hessian
                | Self::HessianDiagonal
                | Self::HessianDiagonalWithSteps
                | Self::Laplacian
        )
    }

    const fn extra(self) -> bool {
        self.uses_index() || matches!(self, Self::HessianDiagonalWithSteps)
    }

    const fn uses_index(self) -> bool {
        matches!(self, Self::Partial | Self::Jacobian)
    }

    const fn only_second_diagonal(self) -> bool {
        matches!(
            self,
            Self::HessianDiagonal | Self::HessianDiagonalWithSteps | Self::Laplacian
        )
    }
}
impl EvalDyn for Method {
    fn signature(&self) -> Signature {
        let d = var("D");
        let e = var("E");
        let derivative = DimExpr::Mul(vec![
            e.clone(),
            DimExpr::Inverse(Box::new(if self.second() {
                DimExpr::Pow(
                    Box::new(d.clone()),
                    crate::kernel::signature::Exp::Rat(crate::kernel::signature::Rational::new(
                        2, 1,
                    )),
                )
            } else {
                d.clone()
            })),
        ]);
        let x = if self.scalar_input() {
            real(d.clone())
        } else {
            tensor_type(1, d.clone())
        };
        let y = if self.vector_output() {
            tensor_type(1, e.clone())
        } else {
            real(e.clone())
        };
        let mut args = vec![x.clone(), real(d.clone()), real(e), TypeExpr::Int];
        if self.uses_index() {
            args.push(TypeExpr::Int);
        } else if matches!(self, Self::HessianDiagonalWithSteps) {
            args.push(tensor_type(1, d));
        }
        args.push(function(vec![x], y));
        let ret = match self {
            | Self::Gradient
            | Self::Curl
            | Self::HessianDiagonal
            | Self::HessianDiagonalWithSteps => tensor_type(1, derivative),
            | Self::Hessian | Self::Jacobian => tensor_type(2, derivative),
            | _ => real(derivative),
        };
        typed(args, ret, &["D", "E"])
    }

    // Keep each native method's validation, callback conversion and result units together.
    #[allow(clippy::too_many_lines)]
    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let expected = self.signature();
        super::arity(args, expected.params.len())?;
        let start = if self.extra() { 5 } else { 4 };
        let callbacks = Callbacks::new(ctx, args, &expected, start, count(&args[3], 3)?)?;
        let sig = &callbacks.signatures[0];
        let d = if self.scalar_input() {
            scalar_dimension(&sig.params[0])?
        } else {
            vector_dimension(&sig.params[0], None)?
        };
        let e = if self.vector_output() {
            vector_dimension(&sig.ret, None)?
        } else {
            scalar_dimension(&sig.ret)?
        };
        let xs = scale(&args[1], &d, 1, ctx)?;
        let ys = scale(&args[2], &e, 2, ctx)?;
        let point = if self.scalar_input() {
            vec![scalar(&args[0], &d, 0, ctx)? / xs]
        } else {
            vector(&args[0], &d, None, 0, ctx)?
                .into_iter()
                .map(|x| x / xs)
                .collect()
        };
        if point.is_empty() {
            return Err(super::domain("calculus point must be nonempty"));
        }
        let width = if self.only_second_diagonal() {
            1
        } else {
            point.len()
        };
        bounded_size(point.len(), width)?;
        if !self.scalar_input() {
            equal(&vector_dimension(&sig.params[0], Some(point.len()))?, &d)?;
        }
        let scalar_callback = |x: &[f64]| {
            let input = if self.scalar_input() {
                super::super::dimensional_float(x[0] * xs, d.clone())
            } else {
                scaled_vec_value(x, xs, &d)
            };
            callbacks
                .call(0, &[input], |v| scalar(&v, &e, start, ctx).map(|v| v / ys))
                .unwrap_or(f64::NAN)
        };
        let requested = if self.uses_index() {
            let Value::Int(n) = &args[4] else {
                return Err(super::domain("calculus index/count must be an integer"));
            };
            n.to_usize()
                .ok_or_else(|| super::domain("calculus index/count must fit usize"))?
        } else {
            0
        };
        let vector_callback = |x: &[f64], length: usize| {
            callbacks
                .call(0, &[scaled_vec_value(x, xs, &d)], |v| {
                    vector(&v, &e, Some(length), start, ctx)
                        .map(|v| v.into_iter().map(|v| v / ys).collect())
                })
                .unwrap_or_else(|| vec![f64::NAN; length])
        };
        let out = match self {
            | Self::Derivative => vec![calculus::derivative(|x| scalar_callback(&[x]), point[0])],
            | Self::Richardson => {
                vec![calculus::richardson_derivative(
                    &mut |x| scalar_callback(&[x]),
                    point[0],
                )]
            },
            | Self::Second => {
                vec![calculus::second_derivative(
                    &mut |x| scalar_callback(&[x]),
                    point[0],
                )]
            },
            | Self::Partial => {
                if requested >= point.len() {
                    return Err(super::domain("partial derivative index is out of range"));
                }
                vec![calculus::partial_derivative(
                    scalar_callback,
                    &point,
                    requested,
                )]
            },
            | Self::Gradient => calculus::gradient(scalar_callback, &point),
            | Self::Hessian => calculus::hessian(scalar_callback, &point),
            | Self::HessianDiagonal => calculus::hessian_diagonal(scalar_callback, &point),
            | Self::HessianDiagonalWithSteps => {
                let mut steps = vector(&args[4], &d, Some(point.len()), 4, ctx)?;
                for step in &mut steps {
                    *step /= xs;
                }
                calculus::hessian_diagonal_with_steps(scalar_callback, &point, &steps)
                    .map_err(|error| super::domain(error.to_string()))?
            },
            | Self::Laplacian => vec![calculus::laplacian(scalar_callback, &point)],
            | Self::Jacobian => {
                if requested == 0 {
                    return Err(super::domain("Jacobian output count must be positive"));
                }
                bounded_size(requested, point.len())?;
                equal(&vector_dimension(&sig.ret, Some(requested))?, &e)?;
                calculus::jacobian(
                    |x, y| y.copy_from_slice(&vector_callback(x, requested)),
                    &point,
                    requested,
                )
            },
            | Self::Divergence => {
                equal(&vector_dimension(&sig.ret, Some(point.len()))?, &e)?;
                vec![calculus::divergence(
                    |x| vector_callback(x, point.len()),
                    &point,
                )]
            },
            | Self::Curl => {
                if point.len() != 3 {
                    return Err(super::domain("curl requires a three-component point"));
                }
                equal(&vector_dimension(&sig.ret, Some(3))?, &e)?;
                calculus::curl(|x| vector_callback(x, 3), &point)
                    .ok_or_else(|| super::domain("curl requires three output components"))?
            },
        };
        callbacks.finish()?;
        let coefficient = if self.second() {
            (ys / xs) / xs
        } else {
            ys / xs
        };
        let data = out
            .into_iter()
            .map(|v| finite(v * coefficient))
            .collect::<Result<Vec<_>, _>>()?;
        let dim = if self.second() {
            ratio(&ratio(&e, &d)?, &d)?
        } else {
            ratio(&e, &d)?
        };
        match self {
            | Self::Gradient
            | Self::Curl
            | Self::HessianDiagonal
            | Self::HessianDiagonalWithSteps => super::output_tensor(vec![data.len()], data, dim),
            | Self::Hessian => super::output_tensor(vec![point.len(), point.len()], data, dim),
            | Self::Jacobian => super::output_tensor(vec![requested, point.len()], data, dim),
            | _ => Ok(super::super::dimensional_float(data[0], dim)),
        }
    }
}
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    for (name, path, method) in [
        ("derivative", "derivative", Method::Derivative),
        (
            "richardson_derivative",
            "richardson_derivative",
            Method::Richardson,
        ),
        ("second_derivative", "second_derivative", Method::Second),
        ("partial_derivative", "partial_derivative", Method::Partial),
        ("gradient", "gradient", Method::Gradient),
        ("hessian", "hessian", Method::Hessian),
        (
            "hessian_diagonal",
            "hessian_diagonal",
            Method::HessianDiagonal,
        ),
        (
            "hessian_diagonal_with_steps",
            "hessian_diagonal_with_steps",
            Method::HessianDiagonalWithSteps,
        ),
        ("jacobian", "jacobian", Method::Jacobian),
        ("divergence", "divergence", Method::Divergence),
        ("curl", "curl", Method::Curl),
        ("laplacian", "laplacian", Method::Laplacian),
    ] {
        let doc = match method {
            | Method::HessianDiagonal | Method::HessianDiagonalWithSteps => {
                "Evaluate the diagonal of the Hessian of a scalar callback at a measured point, returning one value per coordinate in output units divided by squared coordinate units. The common arguments are point, positive coordinate reference, positive output reference and evaluation budget. The with_steps variant takes a positive step per coordinate, in coordinate units, before the callback."
            },
            | _ => {
                "Native finite-difference calculus with physical reference scales: point,coordinate_scale,value_scale,max_evaluations,[partial index or Jacobian output count],callback. Independent D/E units; derivative E/D, Hessian/Laplacian E/D^2. References are positive. Shared callback errors and finite evaluation budgets are preserved."
            },
        };
        register(
            i,
            &format!("num_native_{name}_scaled"),
            &format!("crate::kernels::calculus::{path}"),
            doc,
            if method.extra() { 6 } else { 5 },
            method,
        )?;
    }
    Ok(())
}
