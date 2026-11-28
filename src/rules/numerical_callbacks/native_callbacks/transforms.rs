//! Typed finite complex transforms and bounded inverse Laplace callbacks.
#![allow(
    clippy::arithmetic_side_effects,
    reason = "bounded control arithmetic and floating-point error estimates"
)]

use std::sync::Arc;

use num_complex::Complex64 as C;

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
use super::arity;
use super::complex::complex_type;
use super::complex::decode;
use super::complex::dimension;
use super::complex::value;
use super::count;
use super::domain;
use super::finite;
use super::function;
use super::numeric;
use super::positive;
use super::real;
use super::record;
use super::register;
use super::scalar;
use super::typed;
use super::var;
use crate::graph::OpDescriptor;
use crate::graph::dim::DimRule;
use crate::graph::op::Arity;
use crate::kernel::signature::BackendSupport;
use crate::kernel::signature::Stability;
use crate::kernels::complex_transforms as transforms;
use crate::kernels::complex_transforms::LaplaceContour;
use crate::kernels::contour_deform::ContourIntegral;
use crate::kernels::contour_deform::ContourQuadratureError;
use crate::kernels::quadrature::GkError;
use crate::rules::numerical_abi::codecs::Arg;

const DIMLESS: DimExpr = DimExpr::Dimless;
const BAD_SAMPLE: C = C::new(f64::NAN, f64::NAN);

#[derive(Clone, Copy)]
enum Transform {
    Mellin,
    InverseMellin,
    Hankel,
    Hilbert,
    Stieltjes,
    Laplace(LaplaceContour, bool),
    Stehfest,
}

#[derive(Clone, Copy)]
enum Algebraic {
    Borel,
    MellinBarnes,
    Meijer,
}

fn failed(error: ContourQuadratureError) -> OpError {
    match error {
        | ContourQuadratureError::Quadrature(GkError::PanelBudget) => OpError::BudgetExceeded,
        | ContourQuadratureError::Quadrature(GkError::Cancelled) => {
            crate::cancellation::Cancelled.into()
        },
        | error => domain(error.to_string()),
    }
}

fn estimate_type(
    d: DimExpr,
    is_real: bool,
) -> TypeExpr {
    let value = if is_real {
        real(d.clone())
    } else {
        complex_type(d.clone())
    };
    TypeExpr::Record(vec![
        ("value".into(), value),
        ("error".into(), real(d)),
        ("converged".into(), TypeExpr::Bool),
        ("evaluations".into(), TypeExpr::Int),
    ])
}

fn estimate_value(
    estimate: ContourIntegral,
    d: &DimExpr,
    is_real: bool,
) -> Result<Value, OpError> {
    let scalar = |v| super::super::dimensional_float(v, d.clone());
    let z = C::new(finite(estimate.value.re)?, finite(estimate.value.im)?);
    Ok(record(vec![
        (
            "value",
            if is_real {
                scalar(z.re)
            } else {
                value(z, d.clone())
            },
        ),
        ("error", scalar(finite(estimate.error)?)),
        ("converged", Value::Bool(estimate.converged)),
        ("evaluations", Value::int(estimate.evaluations)),
    ]))
}

fn compare(
    fine: C,
    coarse: C,
    tolerance: f64,
) -> ContourIntegral {
    let error = (fine - coarse).norm();
    ContourIntegral {
        value: fine,
        error,
        evaluations: 0,
        converged: error <= tolerance,
    }
}

impl Transform {
    const fn real_callback(self) -> bool {
        matches!(self, Self::Hilbert | Self::Stehfest)
    }

    const fn complex_input(self) -> bool {
        matches!(self, Self::InverseMellin | Self::Laplace(_, _))
    }

    const fn leading(self) -> usize {
        match self {
            | Self::InverseMellin | Self::Hankel => 3,
            | Self::Mellin | Self::Hilbert | Self::Stieltjes | Self::Laplace(_, true) => 2,
            | Self::Laplace(_, false) | Self::Stehfest => 1,
        }
    }
}

impl EvalDyn for Transform {
    fn signature(&self) -> Signature {
        let real_d = real(var("E"));
        let mut params = match self {
            | Self::Mellin | Self::Stieltjes => vec![TypeExpr::Complex, TypeExpr::Float],
            | Self::Hilbert | Self::Laplace(_, true) => vec![TypeExpr::Float; 2],
            | Self::InverseMellin => vec![TypeExpr::Complex, TypeExpr::Float, TypeExpr::Float],
            | Self::Hankel => vec![TypeExpr::Float; 3],
            | Self::Laplace(_, false) | Self::Stehfest => vec![TypeExpr::Float],
        };
        let input = if self.complex_input() {
            TypeExpr::Complex
        } else {
            TypeExpr::Float
        };
        let output = if self.real_callback() {
            real_d.clone()
        } else {
            complex_type(var("E"))
        };
        params.extend([
            real_d,
            TypeExpr::Int,
            TypeExpr::Int,
            function(vec![input], output),
        ]);
        typed(
            params,
            estimate_type(var("E"), self.real_callback()),
            &["E"],
        )
    }

    #[allow(
        clippy::indexing_slicing,
        reason = "the explicit signature checks arity before fixed arguments are read"
    )]
    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let signature = self.signature();
        arity(args, signature.params.len())?;
        let leading = self.leading();
        let work = count(&args[leading + 1], leading + 1)?;
        let budget = count(&args[leading + 2], leading + 2)?;
        let callback = leading + 3;
        let callbacks = Callbacks::new(ctx, args, &signature, callback, budget)?;
        let ret = &callbacks.signatures[0].ret;
        let dim = if self.real_callback() {
            super::scalar_dimension(ret)?
        } else {
            dimension(ret)?
        };
        let tolerance = scalar(&args[leading], &dim, leading, ctx)?;
        if tolerance <= 0.0 {
            return Err(domain("expected a positive absolute tolerance"));
        }
        let eval_real = |x| {
            callbacks
                .call(0, &[Value::Float(x)], |v| scalar(&v, &dim, callback, ctx))
                .unwrap_or(f64::NAN)
        };
        let eval_complex = |x| {
            callbacks
                .call(0, &[x], |v| decode(&v, &dim, callback, ctx))
                .unwrap_or(BAD_SAMPLE)
        };
        let z = |index: usize| decode(&args[index], &DIMLESS, index, ctx);
        let x = |index: usize| numeric(&args[index], index);
        let bound = |index: usize| positive(&args[index], index, ctx);
        let complex_at_real = |r| eval_complex(Value::Float(r));
        let complex_at_complex = |s| eval_complex(value(s, DIMLESS));
        let result = match self {
            | Self::Mellin => {
                transforms::mellin_checked(complex_at_real, z(0)?, bound(1)?, tolerance, work)
                    .map_err(failed)
            },
            | Self::InverseMellin => {
                transforms::inverse_mellin_checked(
                    complex_at_complex,
                    z(0)?,
                    x(1)?,
                    bound(2)?,
                    tolerance,
                    work,
                )
                .map_err(failed)
            },
            | Self::Hilbert => {
                transforms::hilbert_checked(eval_real, x(0)?, bound(1)?, tolerance, work)
                    .map_err(failed)
            },
            | Self::Hankel => {
                transforms::hankel_checked(
                    complex_at_real,
                    x(0)?,
                    x(1)?,
                    bound(2)?,
                    tolerance,
                    work,
                )
                .map_err(failed)
            },
            | Self::Stieltjes => {
                transforms::stieltjes_checked(complex_at_real, z(0)?, bound(1)?, tolerance, work)
                    .map_err(failed)
            },
            | Self::Laplace(contour, shifted) => {
                if !(8..=64).contains(&work) {
                    return Err(domain("Laplace contour resolution must be in 8..=64"));
                }
                let coarse = work / 2;
                if budget < contour.samples(work) + contour.samples(coarse) {
                    return Err(OpError::BudgetExceeded);
                }
                let sigma = if *shifted { x(1)? } else { 0.0 };
                let t = bound(0)?;
                let sample =
                    |n| transforms::invert_laplace(complex_at_complex, t, sigma, n, *contour);
                sample(work)
                    .and_then(|fine| sample(coarse).map(|coarse| compare(fine, coarse, tolerance)))
                    .map_err(domain)
            },
            | Self::Stehfest => {
                if !(4..=14).contains(&work) || !work.is_multiple_of(2) {
                    return Err(domain(
                        "Stehfest resolution must be an even integer in 4..=14",
                    ));
                }
                if budget < 2 * work - 2 {
                    return Err(OpError::BudgetExceeded);
                }
                let t = bound(0)?;
                let sample = |n| {
                    transforms::gaver_stehfest_inverse_laplace(eval_real, t, n)
                        .map(|r| C::new(r, 0.0))
                };
                sample(work)
                    .and_then(|fine| {
                        sample(work - 2).map(|coarse| compare(fine, coarse, tolerance))
                    })
                    .map_err(domain)
            },
        };
        let evaluations = callbacks.used.get();
        callbacks.finish()?;
        crate::cancellation::checkpoint()?;
        let mut estimate = result?;
        estimate.evaluations = evaluations;
        estimate_value(estimate, &dim, self.real_callback())
    }
}

fn complex_vector(
    input: &Value,
    index: usize,
    max: usize,
    dim: &DimExpr,
    ctx: &CallCtx<'_>,
) -> Result<Vec<C>, OpError> {
    let Value::Tensor(tensor) = input else {
        return Err(domain("expected a rank-one complex tensor"));
    };
    if tensor.shape != [tensor.data.len()] {
        return Err(domain("expected a rank-one complex tensor"));
    }
    if tensor.data.len() > max {
        return Err(OpError::BudgetExceeded);
    }
    tensor
        .data
        .iter()
        .map(|element| {
            // Tensor-level dimensions also annotate plain scalar storage.
            if element.type_expr() == TypeExpr::Complex {
                <C as Arg>::decode(element, index)
            } else {
                decode(element, dim, index, ctx)
            }
        })
        .collect()
}

impl EvalDyn for Algebraic {
    fn signature(&self) -> Signature {
        let vector = |dim| {
            TypeExpr::Tensor {
                rank: crate::kernel::signature::NatExpr::Lit(1),
                shape: None,
                elem: Box::new(complex_type(dim)),
            }
        };
        if matches!(self, Self::Borel) {
            let mut sig = Signature::new(
                vec![vector(var("E")), TypeExpr::Complex],
                complex_type(var("E")),
            );
            sig.type_params = typed(vec![], TypeExpr::Complex, &["E"]).type_params;
            return sig;
        }
        let mut params = if matches!(self, Self::MellinBarnes) {
            vec![vector(DIMLESS); 4]
        } else {
            vec![vector(DIMLESS); 2]
        };
        if matches!(self, Self::Meijer) {
            params.extend([TypeExpr::Int, TypeExpr::Int]);
        }
        params.extend([
            TypeExpr::Complex,
            TypeExpr::Float,
            TypeExpr::Float,
            TypeExpr::Float,
            TypeExpr::Int,
        ]);
        Signature::new(params, estimate_type(DIMLESS, false))
    }

    #[allow(
        clippy::indexing_slicing,
        reason = "the explicit signature checks arity before fixed arguments are read"
    )]
    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        arity(args, self.signature().params.len())?;
        if matches!(self, Self::Borel) {
            let TypeExpr::Tensor { elem, .. } = ctx.value_type(&args[0])? else {
                return Err(domain("expected a rank-one complex tensor"));
            };
            let dim = dimension(&elem)?;
            let coefficients = complex_vector(&args[0], 0, 1_000_000, &dim, ctx)?;
            let t = decode(&args[1], &DIMLESS, 1, ctx)?;
            let v = transforms::borel_transform_series(&coefficients, t);
            return Ok(value(C::new(finite(v.re)?, finite(v.im)?), dim));
        }
        let array = |i: usize| {
            let actual = ctx.value_type(&args[i])?;
            match actual {
                | TypeExpr::Tensor { elem, .. } if *elem == TypeExpr::Complex => {
                    complex_vector(&args[i], i, 32, &DIMLESS, ctx)
                },
                | _ => {
                    Err(domain(
                        "Mellin–Barnes parameters must be dimensionless complex vectors",
                    ))
                },
            }
        };
        let a = array(0)?;
        let b = array(1)?;
        let z = decode(&args[4], &DIMLESS, 4, ctx)?;
        let c = numeric(&args[5], 5)?;
        let limit = positive(&args[6], 6, ctx)?;
        let tol = positive(&args[7], 7, ctx)?;
        let panels = count(&args[8], 8)?;
        let estimate = match self {
            | Self::MellinBarnes => {
                transforms::MellinBarnes {
                    num_plus: &a,
                    num_minus: &b,
                    den_plus: &array(2)?,
                    den_minus: &array(3)?,
                }
                .integrate(z, c, limit, tol, panels)
            },
            | Self::Meijer => {
                let m = <usize as Arg>::decode(&args[2], 2)?;
                let n = <usize as Arg>::decode(&args[3], 3)?;
                transforms::meijer_checked(&a, &b, m, n, z, c, limit, tol, panels)
            },
            | Self::Borel => unreachable!(),
        }
        .map_err(failed)?;
        estimate_value(estimate, &DIMLESS, false)
    }
}

#[allow(
    clippy::single_call_fn,
    reason = "installation entry point for the transform family"
)]
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    for (name, source, doc, eval) in [
        (
            "mellin",
            "mellin_transform",
            "Truncated Mellin transform between exp(-limit) and exp(limit). The absolute error estimate covers finite quadrature only.",
            Transform::Mellin,
        ),
        (
            "inverse_mellin",
            "inverse_mellin_transform",
            "Truncated inverse Mellin transform on a caller-selected vertical line. The absolute error estimate excludes the omitted tail.",
            Transform::InverseMellin,
        ),
        (
            "hilbert",
            "hilbert_transform_real",
            "Finite Hilbert transform with a real density and singularity subtraction.",
            Transform::Hilbert,
        ),
        (
            "hankel",
            "hankel_transform",
            "Truncated Hankel transform of the given order and wavenumber, on [0, limit]. Coordinates and wavenumber are dimensionless.",
            Transform::Hankel,
        ),
        (
            "stieltjes",
            "stieltjes_transform",
            "Finite Stieltjes transform on [0, limit], with the spectral argument off the negative real cut.",
            Transform::Stieltjes,
        ),
        (
            "laplace_talbot",
            "bromwich_talbot",
            "Unshifted Talbot Laplace inversion with an empirical error from a second grid.",
            Transform::Laplace(LaplaceContour::Talbot, false),
        ),
        (
            "laplace_talbot_shifted",
            "bromwich_talbot_shifted",
            "Shifted Talbot Laplace inversion. The caller selects a shift to the right of relevant singularities.",
            Transform::Laplace(LaplaceContour::Talbot, true),
        ),
        (
            "laplace_parabola",
            "bromwich_weideman_parabola",
            "Unshifted parabolic Laplace inversion with an empirical error from a second grid.",
            Transform::Laplace(LaplaceContour::Parabola, false),
        ),
        (
            "laplace_parabola_shifted",
            "bromwich_weideman_parabola_shifted",
            "Shifted parabolic Laplace inversion. The caller selects a shift to the right of relevant singularities.",
            Transform::Laplace(LaplaceContour::Parabola, true),
        ),
        (
            "laplace_stehfest",
            "gaver_stehfest_inverse_laplace",
            "Real-axis Gaver–Stehfest inversion, comparing successive even orders to estimate error.",
            Transform::Stehfest,
        ),
    ] {
        register(
            i,
            &format!("num_transform_{name}_checked"),
            &format!("crate::kernels::complex_transforms::{source}"),
            doc,
            (eval.leading() + 4) as u8,
            eval,
        )?;
    }
    for (name, source, doc, eval) in [
        (
            "borel_series",
            "borel_transform_series",
            "Finite Borel polynomial with dimensionless argument and a shared coefficient dimension.",
            Algebraic::Borel,
        ),
        (
            "mellin_barnes_checked",
            "mellin_barnes_integral",
            "Finite Mellin–Barnes quadrature on a pole-free segment. The error estimate excludes contour tails.",
            Algebraic::MellinBarnes,
        ),
        (
            "meijer_g_checked",
            "meijer_g",
            "Finite Meijer G quadrature with a caller-selected pole-separating line and explicit termination information.",
            Algebraic::Meijer,
        ),
    ] {
        let signature = eval.signature();
        let op = i.op(OpDescriptor::new(
            &format!("num_transform_{name}"),
            Arity::Fixed(signature.arity() as u8),
        )
        .rust_path(&format!("crate::kernels::complex_transforms::{source}"))
        .category("numerical_special")
        .doc(doc)
        .stability(Stability::Experimental)
        .backend(BackendSupport {
            interpreter: true,
            jit: false,
            ffi: true,
            gpu: false,
        })
        .signature(signature)
        .eval_dyn(Arc::new(eval)))?;
        i.dim(op, DimRule::Container);
    }
    Ok(())
}
