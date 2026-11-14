//! Explicit normalized callback protocols with positive physical reference scales.
#[path = "native_analysis/newton_cg.rs"]
mod newton_cg;
#[path = "native_analysis/systems.rs"]
mod systems;
#[path = "native_analysis/visitors.rs"]
mod visitors;
use super::CallCtx;
use super::DimExpr;
use super::EvalDyn;
use super::Installer;
use super::OpError;
use super::RuleError;
use super::Signature;
use super::TypeExpr;
use super::Value;
use super::common::register;
use super::evolution::tensor_type;
use super::evolution::vector;
use super::native_callbacks::Callbacks;
use super::native_callbacks::count;
use super::native_callbacks::domain;
use super::native_callbacks::finite;
use super::native_callbacks::function;
use super::native_callbacks::real;
use super::native_callbacks::record;
use super::native_callbacks::scalar;
use super::native_callbacks::scale;
use super::native_callbacks::typed;
use super::native_callbacks::var;
use crate::kernel::signature::NatExpr;
use crate::kernel::value::VariantValue;
use crate::kernels::integrate as i;
use crate::kernels::quadrature as q;
use crate::kernels::rootfind as r;
use crate::kernels::solve as s;
use num_traits::ToPrimitive;
#[derive(Clone, Copy)]
enum Mode {
    Integral,
    Root,
    Sum,
}
#[derive(Clone, Copy)]
struct Scalar {
    kind: usize,
    mode: Mode,
    np: usize,
    ni: usize,
    callbacks: usize,
}
fn rule(k: usize) -> Result<q::GkRule, OpError> {
    match k {
        | 0 => Ok(q::GkRule::G7K15),
        | 1 => Ok(q::GkRule::G10K21),
        | _ => Err(domain("Gauss-Kronrod rule is0(G7K15)or1(G10K21)")),
    }
}
const fn integral(v: q::Integral) -> (f64, Option<f64>, Option<bool>) {
    (v.value, Some(v.error), Some(v.converged))
}
impl EvalDyn for Scalar {
    fn signature(&self) -> Signature {
        let d = var("D");
        let e = var("E");
        let result = match self.mode {
            | Mode::Integral => DimExpr::Mul(vec![d.clone(), e.clone()]),
            | Mode::Root => d.clone(),
            | Mode::Sum => e.clone(),
        };
        let mut args = vec![
            tensor_type(1, DimExpr::Dimless),
            TypeExpr::Tensor {
                rank: NatExpr::Lit(1),
                shape: None,
                elem: Box::new(TypeExpr::Int),
            },
            real(d),
            real(e),
            TypeExpr::Int,
        ];
        for _ in 0..self.callbacks {
            args.push(function(vec![TypeExpr::Float], TypeExpr::Float));
        }
        typed(
            args,
            TypeExpr::Record(vec![
                ("value".into(), real(result.clone())),
                (
                    "error".into(),
                    TypeExpr::Variant(vec![
                        ("None".into(), vec![]),
                        ("Some".into(), vec![real(result)]),
                    ]),
                ),
                ("evaluations".into(), TypeExpr::Int),
                (
                    "converged".into(),
                    TypeExpr::Variant(vec![
                        ("None".into(), vec![]),
                        ("Some".into(), vec![TypeExpr::Bool]),
                    ]),
                ),
            ]),
            &["D", "E"],
        )
    }

    #[allow(clippy::too_many_lines)] // Native dispatch alternatives retain their exact, independently audited call signatures.
    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let signature = self.signature();
        let callbacks = Callbacks::new(ctx, args, &signature, 5, count(&args[4], 4)?)?;
        let p = vector(&args[0], &DimExpr::Dimless, Some(self.np), 0, ctx)?;
        let Value::Tensor(controls) = &args[1] else {
            return Err(domain("integer controls must be a tensor"));
        };
        if controls.shape != [self.ni] || controls.data.len() != self.ni {
            return Err(domain("native integer control count differs from protocol"));
        }
        let k = controls
            .data
            .iter()
            .map(|v| {
                let Value::Int(n) = v else {
                    return Err(domain("control must be exact Int"));
                };
                n.to_usize()
                    .filter(|n| *n <= 1_000_000)
                    .ok_or_else(|| domain("native integer control outside[0,1000000]"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        if self.kind == 10 && k[0] > 20 {
            return Err(domain("Romberg levels exceed bounded 20-level protocol"));
        }
        let d = super::scalar_dimension(&ctx.value_type(&args[2])?)?;
        let e = super::scalar_dimension(&ctx.value_type(&args[3])?)?;
        let xs = scale(&args[2], &d, 2, ctx)?;
        let ys = scale(&args[3], &e, 3, ctx)?;
        let mut f = |x| {
            callbacks
                .call(0, &[Value::Float(x)], |v| {
                    scalar(&v, &DimExpr::Dimless, 5, ctx)
                })
                .unwrap_or(f64::NAN)
        };
        let g = |x| {
            callbacks
                .call(1, &[Value::Float(x)], |v| {
                    scalar(&v, &DimExpr::Dimless, 6, ctx)
                })
                .unwrap_or(f64::NAN)
        };
        let native: Result<(f64, Option<f64>, Option<bool>), OpError> = match self.kind {
            | 0 => Ok(integral(q::integrate(f, p[0], p[1], p[2]))),
            | 1 => Ok(integral(q::exp_sinh(f, p[0], p[1]))),
            | 2 => Ok(integral(q::sinh_sinh(f, p[0]))),
            | 3 => Ok(integral(q::tanh_sinh(f, p[0], p[1], p[2]))),
            | 4 => Ok((q::filon(f, p[0], p[1], p[2], k[0] != 0, k[1]), None, None)),
            | 5 => {
                let (v, e, _) = q::gk_panel(&f, p[0], p[1], rule(k[0])?);
                Ok((v, Some(e), None))
            },
            | 6 => Ok((i::adaptive_quadrature(f, (p[0], p[1]), p[2]), None, None)),
            | 7 => {
                let r = i::gauss_kronrod(f, p[0], p[1], p[2], k[0]);
                Ok((r.value, Some(r.error), Some(r.error <= p[2])))
            },
            | 8 => {
                let r = i::gauss_kronrod_any(f, p[0], p[1], p[2], k[0]);
                Ok((r.value, Some(r.error), Some(r.error <= p[2])))
            },
            | 9 => Ok((i::gauss_legendre_quadrature(f, (p[0], p[1])), None, None)),
            | 10 => Ok((i::romberg_integration(f, (p[0], p[1]), k[0]), None, None)),
            | 11 => {
                i::simpson_rule(f, (p[0], p[1]), k[0])
                    .map(|v| (v, None, Some(true)))
                    .map_err(domain)
            },
            | 12 => Ok((i::trapezoidal_rule(f, (p[0], p[1]), k[0]), None, None)),
            | 13 => {
                r::newton_bracketed(f, g, p[0], p[1], p[2], k[0])
                    .map(|v| (v, None, Some(true)))
                    .map_err(|e| domain(format!("{e:?}")))
            },
            | 14 => {
                r::root_near(&mut f, p[0])
                    .map(|v| (v, None, Some(true)))
                    .map_err(|e| domain(format!("{e:?}")))
            },
            | 15 => {
                s::solve_root(f, p[0], p[1], k[0])
                    .map(|v| (v, None, Some(true)))
                    .map_err(|e| domain(format!("{e:?}")))
            },
            | 16 => {
                s::solve_root_bisection(f, (p[0], p[1]), p[2], k[0])
                    .map(|v| (v, None, Some(true)))
                    .map_err(|e| domain(format!("{e:?}")))
            },
            | 17 => {
                s::solve_root_newton(f, g, p[0], p[1], k[0])
                    .map(|v| (v, None, Some(true)))
                    .map_err(|e| domain(format!("{e:?}")))
            },
            | 18 => {
                Ok((
                    crate::kernels::series::sum_range(f, k[0] as i64, k[1] as i64),
                    None,
                    None,
                ))
            },
            | 19 => {
                let r = crate::kernels::series::sum_to_infinity(f, k[0] as i64, p[0], k[1]);
                Ok((r.value, Some(r.error), Some(r.error <= p[0])))
            },
            | _ => return Err(domain("unknown analysis scalar method")),
        };
        let used = callbacks.used.get();
        callbacks.finish()?;
        let (v, error, converged) = native?;
        let factor = match self.mode {
            | Mode::Integral => finite(xs * ys)?,
            | Mode::Root => xs,
            | Mode::Sum => ys,
        };
        let unit = match self.mode {
            | Mode::Integral => DimExpr::Mul(vec![d, e]),
            | Mode::Root => d,
            | Mode::Sum => e,
        };
        let value = super::dimensional_float(finite(v * factor)?, unit.clone());
        let error = error
            .filter(|e| e.is_finite())
            .map(|e| finite(e * factor).map(|v| super::dimensional_float(v, unit)))
            .transpose()?;
        Ok(record(vec![
            ("value", value),
            (
                "error",
                Value::Variant(VariantValue {
                    tag: if error.is_some() {
                        "Some"
                    } else {
                        "None"
                    }
                    .into(),
                    payload: error.into_iter().collect(),
                }),
            ),
            ("evaluations", Value::int(used)),
            (
                "converged",
                Value::Variant(VariantValue {
                    tag: if converged.is_some() {
                        "Some"
                    } else {
                        "None"
                    }
                    .into(),
                    payload: converged.map(Value::Bool).into_iter().collect(),
                }),
            ),
        ]))
    }
}
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    systems::install(i)?;
    newton_cg::install(i)?;
    visitors::install(i)?;
    register(
        i,
        "num_analysis_quadrature_integrate_scaled",
        "crate::kernels::quadrature::integrate",
        "Normalized native quadrature::integrate: parameters Float vector length3, exact nonnegative integer controls length0, positive coordinate reference D,positive value reference E,evaluation budget,then1 normalized Float callbacks. Results restore D*E for integrals,D for roots,E for sums. Error=None means native method supplies no finite error estimate; convergence flag retains native stopping status where available.",
        6,
        Scalar {
            kind: 0,
            mode: Mode::Integral,
            np: 3,
            ni: 0,
            callbacks: 1,
        },
    )?;
    register(
        i,
        "num_analysis_quadrature_exp_sinh_scaled",
        "crate::kernels::quadrature::exp_sinh",
        "Normalized native quadrature::exp_sinh: parameters Float vector length2, exact nonnegative integer controls length0, positive coordinate reference D,positive value reference E,evaluation budget,then1 normalized Float callbacks. Results restore D*E for integrals,D for roots,E for sums. Error=None means native method supplies no finite error estimate; convergence flag retains native stopping status where available.",
        6,
        Scalar {
            kind: 1,
            mode: Mode::Integral,
            np: 2,
            ni: 0,
            callbacks: 1,
        },
    )?;
    register(
        i,
        "num_analysis_quadrature_sinh_sinh_scaled",
        "crate::kernels::quadrature::sinh_sinh",
        "Normalized native quadrature::sinh_sinh: parameters Float vector length1, exact nonnegative integer controls length0, positive coordinate reference D,positive value reference E,evaluation budget,then1 normalized Float callbacks. Results restore D*E for integrals,D for roots,E for sums. Error=None means native method supplies no finite error estimate; convergence flag retains native stopping status where available.",
        6,
        Scalar {
            kind: 2,
            mode: Mode::Integral,
            np: 1,
            ni: 0,
            callbacks: 1,
        },
    )?;
    register(
        i,
        "num_analysis_quadrature_tanh_sinh_scaled",
        "crate::kernels::quadrature::tanh_sinh",
        "Normalized native quadrature::tanh_sinh: parameters Float vector length3, exact nonnegative integer controls length0, positive coordinate reference D,positive value reference E,evaluation budget,then1 normalized Float callbacks. Results restore D*E for integrals,D for roots,E for sums. Error=None means native method supplies no finite error estimate; convergence flag retains native stopping status where available.",
        6,
        Scalar {
            kind: 3,
            mode: Mode::Integral,
            np: 3,
            ni: 0,
            callbacks: 1,
        },
    )?;
    register(
        i,
        "num_analysis_quadrature_filon_scaled",
        "crate::kernels::quadrature::filon",
        "Normalized native quadrature::filon: parameters Float vector length3, exact nonnegative integer controls length2, positive coordinate reference D,positive value reference E,evaluation budget,then1 normalized Float callbacks. Results restore D*E for integrals,D for roots,E for sums. Error=None means native method supplies no finite error estimate; convergence flag retains native stopping status where available.",
        6,
        Scalar {
            kind: 4,
            mode: Mode::Integral,
            np: 3,
            ni: 2,
            callbacks: 1,
        },
    )?;
    register(
        i,
        "num_analysis_quadrature_gk_panel_scaled",
        "crate::kernels::quadrature::gk_panel",
        "Normalized native quadrature::gk_panel: parameters Float vector length2, exact nonnegative integer controls length1, positive coordinate reference D,positive value reference E,evaluation budget,then1 normalized Float callbacks. Results restore D*E for integrals,D for roots,E for sums. Error=None means native method supplies no finite error estimate; convergence flag retains native stopping status where available.",
        6,
        Scalar {
            kind: 5,
            mode: Mode::Integral,
            np: 2,
            ni: 1,
            callbacks: 1,
        },
    )?;
    register(
        i,
        "num_analysis_integrate_adaptive_scaled",
        "crate::kernels::integrate::adaptive_quadrature",
        "Normalized native integrate::adaptive_quadrature: parameters Float vector length3, exact nonnegative integer controls length0, positive coordinate reference D,positive value reference E,evaluation budget,then1 normalized Float callbacks. Results restore D*E for integrals,D for roots,E for sums. Error=None means native method supplies no finite error estimate; convergence flag retains native stopping status where available.",
        6,
        Scalar {
            kind: 6,
            mode: Mode::Integral,
            np: 3,
            ni: 0,
            callbacks: 1,
        },
    )?;
    register(
        i,
        "num_analysis_integrate_gauss_kronrod_scaled",
        "crate::kernels::integrate::gauss_kronrod",
        "Normalized native integrate::gauss_kronrod: parameters Float vector length3, exact nonnegative integer controls length1, positive coordinate reference D,positive value reference E,evaluation budget,then1 normalized Float callbacks. Results restore D*E for integrals,D for roots,E for sums. Error=None means native method supplies no finite error estimate; convergence flag retains native stopping status where available.",
        6,
        Scalar {
            kind: 7,
            mode: Mode::Integral,
            np: 3,
            ni: 1,
            callbacks: 1,
        },
    )?;
    register(
        i,
        "num_analysis_integrate_gauss_kronrod_any_scaled",
        "crate::kernels::integrate::gauss_kronrod_any",
        "Normalized native integrate::gauss_kronrod_any: parameters Float vector length3, exact nonnegative integer controls length1, positive coordinate reference D,positive value reference E,evaluation budget,then1 normalized Float callbacks. Results restore D*E for integrals,D for roots,E for sums. Error=None means native method supplies no finite error estimate; convergence flag retains native stopping status where available.",
        6,
        Scalar {
            kind: 8,
            mode: Mode::Integral,
            np: 3,
            ni: 1,
            callbacks: 1,
        },
    )?;
    register(
        i,
        "num_analysis_integrate_gauss_legendre_scaled",
        "crate::kernels::integrate::gauss_legendre_quadrature",
        "Normalized native integrate::gauss_legendre_quadrature: parameters Float vector length2, exact nonnegative integer controls length0, positive coordinate reference D,positive value reference E,evaluation budget,then1 normalized Float callbacks. Results restore D*E for integrals,D for roots,E for sums. Error=None means native method supplies no finite error estimate; convergence flag retains native stopping status where available.",
        6,
        Scalar {
            kind: 9,
            mode: Mode::Integral,
            np: 2,
            ni: 0,
            callbacks: 1,
        },
    )?;
    register(
        i,
        "num_analysis_integrate_romberg_scaled",
        "crate::kernels::integrate::romberg_integration",
        "Normalized native integrate::romberg_integration: parameters Float vector length2, exact nonnegative integer controls length1, positive coordinate reference D,positive value reference E,evaluation budget,then1 normalized Float callbacks. Results restore D*E for integrals,D for roots,E for sums. Error=None means native method supplies no finite error estimate; convergence flag retains native stopping status where available.",
        6,
        Scalar {
            kind: 10,
            mode: Mode::Integral,
            np: 2,
            ni: 1,
            callbacks: 1,
        },
    )?;
    register(
        i,
        "num_analysis_integrate_simpson_scaled",
        "crate::kernels::integrate::simpson_rule",
        "Normalized native integrate::simpson_rule: parameters Float vector length2, exact nonnegative integer controls length1, positive coordinate reference D,positive value reference E,evaluation budget,then1 normalized Float callbacks. Results restore D*E for integrals,D for roots,E for sums. Error=None means native method supplies no finite error estimate; convergence flag retains native stopping status where available.",
        6,
        Scalar {
            kind: 11,
            mode: Mode::Integral,
            np: 2,
            ni: 1,
            callbacks: 1,
        },
    )?;
    register(
        i,
        "num_analysis_integrate_trapezoidal_scaled",
        "crate::kernels::integrate::trapezoidal_rule",
        "Normalized native integrate::trapezoidal_rule: parameters Float vector length2, exact nonnegative integer controls length1, positive coordinate reference D,positive value reference E,evaluation budget,then1 normalized Float callbacks. Results restore D*E for integrals,D for roots,E for sums. Error=None means native method supplies no finite error estimate; convergence flag retains native stopping status where available.",
        6,
        Scalar {
            kind: 12,
            mode: Mode::Integral,
            np: 2,
            ni: 1,
            callbacks: 1,
        },
    )?;
    register(
        i,
        "num_analysis_newton_bracketed_scaled",
        "crate::kernels::rootfind::newton_bracketed",
        "Normalized native rootfind::newton_bracketed: parameters Float vector length3, exact nonnegative integer controls length1, positive coordinate reference D,positive value reference E,evaluation budget,then2 normalized Float callbacks. Results restore D*E for integrals,D for roots,E for sums. Error=None means native method supplies no finite error estimate; convergence flag retains native stopping status where available.",
        7,
        Scalar {
            kind: 13,
            mode: Mode::Root,
            np: 3,
            ni: 1,
            callbacks: 2,
        },
    )?;
    register(
        i,
        "num_analysis_root_near_scaled",
        "crate::kernels::rootfind::root_near",
        "Normalized native rootfind::root_near: parameters Float vector length1, exact nonnegative integer controls length0, positive coordinate reference D,positive value reference E,evaluation budget,then1 normalized Float callbacks. Results restore D*E for integrals,D for roots,E for sums. Error=None means native method supplies no finite error estimate; convergence flag retains native stopping status where available.",
        6,
        Scalar {
            kind: 14,
            mode: Mode::Root,
            np: 1,
            ni: 0,
            callbacks: 1,
        },
    )?;
    register(
        i,
        "num_analysis_solve_root_scaled",
        "crate::kernels::solve::solve_root",
        "Normalized native solve::solve_root: parameters Float vector length2, exact nonnegative integer controls length1, positive coordinate reference D,positive value reference E,evaluation budget,then1 normalized Float callbacks. Results restore D*E for integrals,D for roots,E for sums. Error=None means native method supplies no finite error estimate; convergence flag retains native stopping status where available.",
        6,
        Scalar {
            kind: 15,
            mode: Mode::Root,
            np: 2,
            ni: 1,
            callbacks: 1,
        },
    )?;
    register(
        i,
        "num_analysis_solve_root_bisection_scaled",
        "crate::kernels::solve::solve_root_bisection",
        "Normalized native solve::solve_root_bisection: parameters Float vector length3, exact nonnegative integer controls length1, positive coordinate reference D,positive value reference E,evaluation budget,then1 normalized Float callbacks. Results restore D*E for integrals,D for roots,E for sums. Error=None means native method supplies no finite error estimate; convergence flag retains native stopping status where available.",
        6,
        Scalar {
            kind: 16,
            mode: Mode::Root,
            np: 3,
            ni: 1,
            callbacks: 1,
        },
    )?;
    register(
        i,
        "num_analysis_solve_root_newton_scaled",
        "crate::kernels::solve::solve_root_newton",
        "Normalized native solve::solve_root_newton: parameters Float vector length2, exact nonnegative integer controls length1, positive coordinate reference D,positive value reference E,evaluation budget,then2 normalized Float callbacks. Results restore D*E for integrals,D for roots,E for sums. Error=None means native method supplies no finite error estimate; convergence flag retains native stopping status where available.",
        7,
        Scalar {
            kind: 17,
            mode: Mode::Root,
            np: 2,
            ni: 1,
            callbacks: 2,
        },
    )?;
    register(
        i,
        "num_analysis_series_sum_range_scaled",
        "crate::kernels::series::sum_range",
        "Normalized native series::sum_range: parameters Float vector length0, exact nonnegative integer controls length2, positive coordinate reference D,positive value reference E,evaluation budget,then1 normalized Float callbacks. Results restore D*E for integrals,D for roots,E for sums. Error=None means native method supplies no finite error estimate; convergence flag retains native stopping status where available.",
        6,
        Scalar {
            kind: 18,
            mode: Mode::Sum,
            np: 0,
            ni: 2,
            callbacks: 1,
        },
    )?;
    register(
        i,
        "num_analysis_series_sum_infinity_scaled",
        "crate::kernels::series::sum_to_infinity",
        "Normalized native series::sum_to_infinity: parameters Float vector length1, exact nonnegative integer controls length2, positive coordinate reference D,positive value reference E,evaluation budget,then1 normalized Float callbacks. Results restore D*E for integrals,D for roots,E for sums. Error=None means native method supplies no finite error estimate; convergence flag retains native stopping status where available.",
        6,
        Scalar {
            kind: 19,
            mode: Mode::Sum,
            np: 1,
            ni: 2,
            callbacks: 1,
        },
    )?;
    Ok(())
}
