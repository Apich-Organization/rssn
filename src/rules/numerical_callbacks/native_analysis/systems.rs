//! Native system roots/Jacobians with normalized callbacks and physical references.
use super::super::CallCtx;
use super::super::DimExpr;
use super::super::EvalDyn;
use super::super::Installer;
use super::super::OpError;
use super::super::RuleError;
use super::super::Signature;
use super::super::TypeExpr;
use super::super::Value;
use super::super::common::register;
use super::super::evolution::output_tensor;
use super::super::evolution::tensor_type;
use super::super::evolution::vector;
use super::super::native_callbacks::Callbacks;
use super::super::native_callbacks::count;
use super::super::native_callbacks::domain;
use super::super::native_callbacks::finite;
use super::super::native_callbacks::function;
use super::super::native_callbacks::real;
use super::super::native_callbacks::record;
use super::super::native_callbacks::scalar;
use super::super::native_callbacks::scale;
use super::super::native_callbacks::typed;
use super::super::native_callbacks::var;
use super::super::native_callbacks::vec_value;
use crate::kernels::dense;
use crate::kernels::optim;
use crate::kernels::rootfind as root;
use crate::kernels::solve;
use crate::rules::numerical_abi::codecs::Arg;
#[derive(Clone, Copy)]
enum Method {
    Lm,
    NewtonFd,
    NewtonJac,
    Solve,
    Jacobian,
    TryJacobian,
    Gradient,
}
impl Method {
    const fn jacobian(self) -> bool {
        matches!(self, Self::Jacobian | Self::TryJacobian)
    }

    const fn gradient(self) -> bool {
        matches!(self, Self::Gradient)
    }

    const fn start(self) -> usize {
        if self.jacobian() {
            5
        } else if self.gradient() {
            4
        } else {
            6
        }
    }
}
impl EvalDyn for Method {
    fn signature(&self) -> Signature {
        let d = var("D");
        let e = var("E");
        let x = tensor_type(1, d.clone());
        let derivative = DimExpr::Mul(vec![e.clone(), DimExpr::Inverse(Box::new(d.clone()))]);
        let mut args = if self.jacobian() {
            vec![
                x.clone(),
                tensor_type(1, e.clone()),
                real(d),
                real(e.clone()),
                TypeExpr::Int,
            ]
        } else if self.gradient() {
            vec![x.clone(), real(d), real(e.clone()), TypeExpr::Int]
        } else {
            vec![
                x.clone(),
                real(d),
                real(e.clone()),
                TypeExpr::Float,
                TypeExpr::Int,
                TypeExpr::Int,
            ]
        };
        let normalized = tensor_type(1, DimExpr::Dimless);
        args.push(function(
            vec![normalized.clone()],
            if self.gradient() {
                TypeExpr::Float
            } else {
                normalized.clone()
            },
        ));
        if matches!(self, Self::NewtonJac) {
            args.push(function(vec![normalized], tensor_type(2, DimExpr::Dimless)));
        }
        let out = if self.jacobian() {
            tensor_type(2, derivative)
        } else if self.gradient() {
            tensor_type(1, derivative)
        } else if matches!(self, Self::Solve) {
            x
        } else {
            TypeExpr::Record(vec![
                ("x".into(), x),
                ("residual".into(), real(e)),
                ("iterations".into(), TypeExpr::Int),
            ])
        };
        typed(args, out, &["D", "E"])
    }

    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let signature = self.signature();
        let start = self.start();
        let budget = count(&args[start - 1], start - 1)?;
        let callbacks = Callbacks::new(ctx, args, &signature, start, budget)?;
        let TypeExpr::Tensor { elem, .. } = ctx.value_type(&args[0])? else {
            return Err(domain("system point must be tensor"));
        };
        let d = super::super::scalar_dimension(&elem)?;
        let scale_at = if self.jacobian() { 2 } else { 1 };
        let e = super::super::scalar_dimension(&ctx.value_type(&args[scale_at + 1])?)?;
        let xs = scale(&args[scale_at], &d, scale_at, ctx)?;
        let ys = scale(&args[scale_at + 1], &e, scale_at + 1, ctx)?;
        let x = vector(&args[0], &d, None, 0, ctx)?
            .into_iter()
            .map(|v| v / xs)
            .collect::<Vec<_>>();
        if x.is_empty() || x.len() > 256 || x.iter().any(|v| !v.is_finite()) {
            return Err(domain(
                "system dimension must be1..256 with finite normalized coordinates",
            ));
        }
        let f = |x: &[f64]| {
            callbacks
                .call(0, &[vec_value(x, &DimExpr::Dimless)], |v| {
                    vector(
                        &v,
                        &DimExpr::Dimless,
                        if self.jacobian() || matches!(self, Self::Lm) {
                            None
                        } else {
                            Some(x.len())
                        },
                        start,
                        ctx,
                    )
                })
                .unwrap_or_else(|| vec![f64::NAN; x.len()])
        };
        let derivative = DimExpr::Mul(vec![e.clone(), DimExpr::Inverse(Box::new(d.clone()))]);
        if self.gradient() {
            let g = optim::numerical_gradient(
                &|x| {
                    callbacks
                        .call(0, &[vec_value(x, &DimExpr::Dimless)], |v| {
                            scalar(&v, &DimExpr::Dimless, start, ctx)
                        })
                        .unwrap_or(f64::NAN)
                },
                &x,
            );
            callbacks.finish()?;
            return output_tensor(
                vec![x.len()],
                g.into_iter()
                    .map(|v| finite(v * (ys / xs)))
                    .collect::<Result<_, _>>()?,
                derivative,
            );
        }
        if self.jacobian() {
            let fx = vector(&args[1], &e, None, 1, ctx)?
                .into_iter()
                .map(|v| v / ys)
                .collect::<Vec<_>>();
            if fx.is_empty() || fx.len() > 256 {
                return Err(domain("Jacobian residual dimension must be1..256"));
            }
            let out = if matches!(self, Self::Jacobian) {
                Ok(root::numerical_jacobian(&f, &x, &fx))
            } else {
                root::try_numerical_jacobian(&f, &x, &fx)
            };
            callbacks.finish()?;
            let out = out.map_err(|e| domain(format!("{e:?}")))?;
            return output_tensor(
                vec![out.rows, out.cols],
                out.data
                    .into_iter()
                    .map(|v| finite(v * (ys / xs)))
                    .collect::<Result<_, _>>()?,
                derivative,
            );
        }
        let tol = scalar(&args[3], &DimExpr::Dimless, 3, ctx)?;
        if tol <= 0.0 {
            return Err(domain("system tolerance must be positive"));
        }
        let iterations = count(&args[4], 4)?;
        if matches!(self, Self::Solve) {
            let out = solve::solve_nonlinear_system(
                |x, out| {
                    let value = f(x);
                    if value.len() == out.len() {
                        out.copy_from_slice(&value);
                    } else {
                        out.fill(f64::NAN);
                    }
                },
                &x,
                tol,
                iterations,
            );
            callbacks.finish()?;
            let out = out.map_err(domain)?;
            return output_tensor(
                vec![out.len()],
                out.into_iter()
                    .map(|v| finite(v * xs))
                    .collect::<Result<_, _>>()?,
                d,
            );
        }
        let jac = |x: &[f64]| {
            callbacks
                .call(1, &[vec_value(x, &DimExpr::Dimless)], |v| {
                    let m = <dense::Mat as Arg>::decode(&v, start + 1)?;
                    if m.rows != x.len() || m.cols != x.len() {
                        return Err(domain(
                            "analytic system Jacobian must be square with state size",
                        ));
                    }
                    Ok(m)
                })
                .unwrap_or_else(|| dense::Mat::zeros(x.len(), x.len()))
        };
        let out = match self {
            | Self::Lm => root::levenberg_marquardt(f, &x, tol, iterations),
            | Self::NewtonFd => {
                root::newton_system(f, None::<&fn(&[f64]) -> dense::Mat>, &x, tol, iterations)
            },
            | Self::NewtonJac => root::newton_system(f, Some(&jac), &x, tol, iterations),
            | _ => return Err(domain("unknown native system method")),
        };
        callbacks.finish()?;
        let out = out.map_err(|e| domain(format!("{e:?}")))?;
        Ok(record(vec![
            (
                "x",
                output_tensor(
                    vec![out.x.len()],
                    out.x
                        .into_iter()
                        .map(|v| finite(v * xs))
                        .collect::<Result<_, _>>()?,
                    d,
                )?,
            ),
            (
                "residual",
                super::super::dimensional_float(finite(out.residual * ys)?, e),
            ),
            ("iterations", Value::int(out.iterations)),
        ]))
    }
}
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    for (name, path, method) in [
        ("lm", "rootfind::levenberg_marquardt", Method::Lm),
        ("newton_fd", "rootfind::newton_system", Method::NewtonFd),
        (
            "newton_jacobian",
            "rootfind::newton_system",
            Method::NewtonJac,
        ),
        (
            "solve_nonlinear",
            "solve::solve_nonlinear_system",
            Method::Solve,
        ),
        (
            "numerical_jacobian",
            "rootfind::numerical_jacobian",
            Method::Jacobian,
        ),
        (
            "try_numerical_jacobian",
            "rootfind::try_numerical_jacobian",
            Method::TryJacobian,
        ),
        (
            "optimization_gradient",
            "optim::numerical_gradient",
            Method::Gradient,
        ),
    ] {
        register(
            i,
            &format!("num_analysis_{name}_scaled"),
            &format!("crate::kernels::{path}"),
            "Genuine native method with physical point and positive coordinate/output references. Callbacks operate normalized vectors; residual/gradient/Jacobian results restore their proper E,D units. Solver args x,xref,yref,normalized tolerance,iterations,budget,callback(s); numerical Jacobian args x,physical fx,xref,yref,budget,f; gradient args x,xref,yref,budget,f. Shared callback budget and original errors are preserved.",
            (method.start() + 1 + usize::from(matches!(method, Method::NewtonJac))) as u8,
            method,
        )?;
    }
    Ok(())
}
