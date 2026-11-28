//! Checked scalar Riemann–Hilbert transforms with measured densities.

use crate::kernel::value::ComplexValue;
use crate::kernels::contour_deform::ContourIntegral;
use crate::kernels::contour_deform::ContourQuadratureError;
use crate::kernels::quadrature::GkError;
use crate::kernels::riemann_hilbert::ScalarRhError;
use crate::kernels::riemann_hilbert::cauchy_transform_arc_checked;
use crate::kernels::riemann_hilbert::plemelj_boundary_values_checked;
use crate::kernels::riemann_hilbert::scalar_factorization_checked;
use num_complex::Complex64;

use super::super::super::dimensional_float;
use super::super::CallCtx;
use super::super::Callbacks;
use super::super::DimExpr;
use super::super::EvalDyn;
use super::super::Installer;
use super::super::OpError;
use super::super::RuleError;
use super::super::Signature;
use super::super::TypeExpr;
use super::super::Value;
use super::super::arity;
use super::super::complex::complex_type;
use super::super::complex::decode;
use super::super::complex::dimension;
use super::super::complex::value;
use super::super::count;
use super::super::domain;
use super::super::finite;
use super::super::function;
use super::super::numeric;
use super::super::real;
use super::super::record;
use super::super::register;
use super::super::scale;
use super::super::typed;
use super::super::var;

#[derive(Clone, Copy)]
enum ScalarMethod {
    Cauchy,
    Plemelj,
    Factor,
}

type ResultFields = Vec<(&'static str, Value)>;

fn fail(error: ScalarRhError) -> OpError {
    match error {
        | ScalarRhError::Quad(ContourQuadratureError::Quadrature(GkError::PanelBudget)) => {
            OpError::BudgetExceeded
        },
        | ScalarRhError::Quad(ContourQuadratureError::Quadrature(GkError::Cancelled)) => {
            crate::cancellation::Cancelled.into()
        },
        | error => domain(error.to_string()),
    }
}

fn estimated_fields(
    estimate: ContourIntegral,
    amplitude: f64,
    dim: &DimExpr,
) -> Result<ResultFields, OpError> {
    let physical = estimate.value * amplitude;
    let physical = Complex64::new(finite(physical.re)?, finite(physical.im)?);
    Ok(vec![
        ("value", value(physical, dim.clone())),
        (
            "error",
            dimensional_float(finite(estimate.error * amplitude)?, dim.clone()),
        ),
        ("converged", Value::Bool(estimate.converged)),
        ("evaluations", Value::int(estimate.evaluations)),
    ])
}

impl EvalDyn for ScalarMethod {
    fn signature(&self) -> Signature {
        let (params, mut fields, dims) = match self {
            | Self::Cauchy | Self::Plemelj => {
                let d = var("E");
                let complex = complex_type(d.clone());
                let first = if matches!(self, Self::Plemelj) {
                    TypeExpr::Float
                } else {
                    TypeExpr::Complex
                };
                let mut result = if matches!(self, Self::Plemelj) {
                    vec![
                        ("upper".into(), complex.clone()),
                        ("lower".into(), complex.clone()),
                    ]
                } else {
                    vec![("value".into(), complex.clone())]
                };
                result.push(("error".into(), real(d.clone())));
                (
                    vec![
                        first,
                        real(d.clone()),
                        real(d),
                        TypeExpr::Int,
                        TypeExpr::Int,
                        function(vec![TypeExpr::Float], complex),
                    ],
                    result,
                    &["E"] as &[_],
                )
            },
            | Self::Factor => {
                (
                    vec![
                        TypeExpr::Complex,
                        TypeExpr::Float,
                        TypeExpr::Int,
                        TypeExpr::Int,
                        TypeExpr::Int,
                        function(vec![TypeExpr::Float], TypeExpr::Complex),
                    ],
                    vec![
                        ("value".into(), TypeExpr::Complex),
                        ("error".into(), TypeExpr::Float),
                        ("logarithm".into(), TypeExpr::Complex),
                        ("log_error".into(), TypeExpr::Float),
                    ],
                    &[] as &[_],
                )
            },
        };
        fields.extend([
            ("converged".into(), TypeExpr::Bool),
            ("evaluations".into(), TypeExpr::Int),
        ]);
        typed(params, TypeExpr::Record(fields), dims)
    }

    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let expected = self.signature();
        arity(args, expected.params.len())?;
        let panels = count(&args[3], 3)?;
        let budget = count(&args[4], 4)?;
        let callbacks = Callbacks::new(ctx, args, &expected, 5, budget)?;
        let dim = dimension(&callbacks.signatures[0].ret)?;
        let amp = if matches!(self, Self::Factor) {
            1.0
        } else {
            scale(&args[1], &dim, 1, ctx)?
        };
        let sample = |t| {
            callbacks
                .call(0, &[Value::Float(t)], |value| {
                    let z = decode(&value, &dim, 5, ctx)? / amp;
                    Ok(Complex64::new(finite(z.re)?, finite(z.im)?))
                })
                .unwrap_or(Complex64::new(f64::NAN, f64::NAN))
        };
        let computed = match self {
            | Self::Cauchy => {
                let z = decode(&args[0], &DimExpr::Dimless, 0, ctx)?;
                let tolerance = scale(&args[2], &dim, 2, ctx)? / amp;
                cauchy_transform_arc_checked(&sample, z, tolerance, panels)
                    .map_err(fail)
                    .and_then(|estimate| estimated_fields(estimate, amp, &dim))
            },
            | Self::Plemelj => {
                let x = numeric(&args[0], 0)?;
                let tolerance = scale(&args[2], &dim, 2, ctx)? / amp;
                plemelj_boundary_values_checked(&sample, x, tolerance, panels)
                    .map_err(fail)
                    .and_then(|estimate| {
                        let [upper, lower] = [estimate.upper * amp, estimate.lower * amp];
                        if !upper.is_finite() || !lower.is_finite() {
                            return Err(domain(
                                "Plemelj boundary values are outside binary64 range",
                            ));
                        }
                        Ok(vec![
                            ("upper", value(upper, dim.clone())),
                            ("lower", value(lower, dim.clone())),
                            (
                                "error",
                                dimensional_float(finite(estimate.error * amp)?, dim.clone()),
                            ),
                            ("converged", Value::Bool(estimate.converged)),
                            ("evaluations", Value::int(estimate.evaluations)),
                        ])
                    })
            },
            | Self::Factor => {
                let z = decode(&args[0], &DimExpr::Dimless, 0, ctx)?;
                let tolerance = numeric(&args[1], 1)?;
                let phase_grid = count(&args[2], 2)?;
                if budget <= phase_grid {
                    return Err(OpError::BudgetExceeded);
                }
                scalar_factorization_checked(&sample, z, tolerance, phase_grid, panels)
                    .map_err(fail)
                    .and_then(|result| {
                        let mut fields = estimated_fields(result.estimate, 1.0, &DimExpr::Dimless)?;
                        fields.extend([
                            (
                                "logarithm",
                                Value::Complex(ComplexValue::Float {
                                    re: result.logarithm.re,
                                    im: result.logarithm.im,
                                }),
                            ),
                            ("log_error", Value::Float(result.log_error)),
                        ]);
                        Ok(fields)
                    })
            },
        };
        callbacks.finish()?;
        crate::cancellation::checkpoint()?;
        Ok(record(computed?))
    }
}

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    register(
        i,
        "num_riemann_hilbert_cauchy_checked",
        "crate::kernels::riemann_hilbert::cauchy_transform_arc",
        "Evaluates the Cauchy transform off [-1, 1] using checked adaptive quadrature. The density, output, amplitude scale and absolute tolerance share one dimension. The result reports the estimated error, convergence and actual host callbacks.",
        6,
        ScalarMethod::Cauchy,
    )?;
    register(
        i,
        "num_riemann_hilbert_plemelj_checked",
        "crate::kernels::riemann_hilbert::plemelj_boundary_values",
        "Evaluates both Plemelj boundary limits inside (-1, 1), using singularity subtraction and a shared error estimate. The density, output, amplitude scale and absolute tolerance share one dimension. Reports convergence and actual host callbacks.",
        6,
        ScalarMethod::Plemelj,
    )?;
    register(
        i,
        "num_riemann_hilbert_factor_checked",
        "crate::kernels::riemann_hilbert::scalar_rh_factorize_arc",
        "Factors a nonzero dimensionless jump on (-1, 1). The supplied mesh must resolve its continuous phase from the principal phase at -1. Reports the estimates and errors before and after exponentiation, convergence against the absolute value tolerance and actual host callbacks.",
        6,
        ScalarMethod::Factor,
    )
}
