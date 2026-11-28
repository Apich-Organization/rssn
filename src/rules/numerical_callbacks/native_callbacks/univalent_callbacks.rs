//! Bounded Schwarzian and unit-disk sampling with measured or normalized callbacks.
use crate::kernel::signature::NatExpr;
use crate::kernels::univalent as native;
use num_complex::Complex64;

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
use super::super::count;
use super::super::domain;
use super::super::function;
use super::super::ratio;
use super::super::real;
use super::super::register;
use super::super::scale;
use super::super::typed;
use super::super::var;
use super::super::vector;
use super::complex_type;
use super::decode;
use super::dimension;
use super::value;

const D0: DimExpr = DimExpr::Dimless;
const BAD: Complex64 = Complex64::new(f64::NAN, f64::NAN);

#[derive(Clone, Copy)]
enum Method {
    Schwarzian,
    Nehari,
    Starlikeness,
    Convexity,
}

fn inverse_square(dim: DimExpr) -> DimExpr {
    DimExpr::Inverse(Box::new(DimExpr::Mul(vec![dim.clone(), dim])))
}

fn bounded_angles(
    args: &[Value],
    arg: usize,
) -> Result<usize, OpError> {
    let angles = count(&args[arg], arg)?;
    if angles > 4096 {
        return Err(OpError::BudgetExceeded);
    }
    Ok(angles)
}

impl Method {
    fn sample(
        callbacks: &Callbacks<'_, '_>,
        input: &DimExpr,
        output: &DimExpr,
        value_scale: f64,
        index: usize,
        ctx: &CallCtx<'_>,
        z: Complex64,
    ) -> Complex64 {
        callbacks
            .call(index, &[value(z, input.clone())], |v| {
                decode(&v, output, index, ctx).map(|result| result / value_scale)
            })
            .unwrap_or(BAD)
    }

    fn evaluate(
        self,
        args: &[Value],
        callbacks: &Callbacks<'_, '_>,
        ctx: &CallCtx<'_>,
        budget: usize,
    ) -> Result<Value, OpError> {
        if matches!(self, Self::Schwarzian) {
            if budget < native::SCHWARZIAN_SAMPLES {
                return Err(OpError::BudgetExceeded);
            }
            let coordinate = dimension(&callbacks.signatures[0].params[0])?;
            let output = dimension(&callbacks.signatures[0].ret)?;
            let point = decode(&args[0], &coordinate, 0, ctx)?;
            let step = scale(&args[1], &coordinate, 1, ctx)?;
            let amplitude = scale(&args[2], &output, 2, ctx)?;
            let result = native::schwarzian_derivative_numeric(
                &|z| Self::sample(callbacks, &coordinate, &output, amplitude, 0, ctx, z),
                point,
                step,
            )
            .map_err(domain)?;
            let units = ratio(&D0, &DimExpr::Mul(vec![coordinate.clone(), coordinate]))?;
            return Ok(value(result, units));
        }
        let sample = |idx, z| Self::sample(callbacks, &D0, &D0, 1.0, idx, ctx, z);
        let angles = bounded_angles(args, 1)?;
        if matches!(self, Self::Nehari) {
            let Some(Value::Tensor(tensor)) = args.first() else {
                return Err(domain("Nehari sample radii need a rank-one tensor"));
            };
            if tensor
                .data
                .len()
                .saturating_mul(angles)
                .saturating_mul(native::SCHWARZIAN_SAMPLES)
                > budget
            {
                return Err(OpError::BudgetExceeded);
            }
            let radii = vector(&args[0], &D0, None, 0, ctx)?;
            let step = scale(&args[2], &D0, 2, ctx)?;
            let pass = native::check_nehari_univalence(&|z| sample(0, z), &radii, angles, step)
                .map_err(domain)?;
            return Ok(Value::Bool(pass));
        }
        if angles.saturating_mul(2 * native::RADIUS_ITERATIONS) > budget {
            return Err(OpError::BudgetExceeded);
        }
        let maximum = scale(&args[0], &D0, 0, ctx)?;
        if maximum > 1.0 {
            return Err(domain("sampled disk radius cannot exceed one"));
        }
        let first = |z| sample(0, z);
        let second = |z| sample(1, z);
        let radius = if matches!(self, Self::Starlikeness) {
            native::radius_of_starlikeness_numeric(&first, &second, maximum, angles)
        } else {
            native::radius_of_convexity_numeric(&first, &second, maximum, angles)
        };
        if radius.is_finite() {
            Ok(Value::Float(radius))
        } else {
            Err(domain(
                "sampled radius needs finite callback values and resolvable ratios",
            ))
        }
    }
}

impl EvalDyn for Method {
    fn signature(&self) -> Signature {
        let normalized = function(vec![TypeExpr::Complex], TypeExpr::Complex);
        match self {
            | Self::Schwarzian => {
                let d = var("D");
                let e = var("E");
                typed(
                    vec![
                        complex_type(d.clone()),
                        real(d.clone()),
                        real(e.clone()),
                        TypeExpr::Int,
                        function(vec![complex_type(d.clone())], complex_type(e)),
                    ],
                    complex_type(inverse_square(d)),
                    &["D", "E"],
                )
            },
            | Self::Nehari => {
                typed(
                    vec![
                        TypeExpr::Tensor {
                            rank: NatExpr::Lit(1),
                            shape: None,
                            elem: Box::new(TypeExpr::Float),
                        },
                        TypeExpr::Int,
                        TypeExpr::Float,
                        TypeExpr::Int,
                        normalized,
                    ],
                    TypeExpr::Bool,
                    &[],
                )
            },
            | Self::Starlikeness | Self::Convexity => {
                typed(
                    vec![
                        TypeExpr::Float,
                        TypeExpr::Int,
                        TypeExpr::Int,
                        normalized.clone(),
                        normalized,
                    ],
                    TypeExpr::Float,
                    &[],
                )
            },
        }
    }

    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let signature = self.signature();
        arity(args, signature.params.len())?;
        let budget_arg = if matches!(self, Self::Starlikeness | Self::Convexity) {
            2
        } else {
            3
        };
        let budget = count(&args[budget_arg], budget_arg)?;
        let callbacks = Callbacks::new(ctx, args, &signature, budget_arg + 1, budget)?;
        let result = self.evaluate(args, &callbacks, ctx, budget);
        callbacks.finish()?;
        crate::cancellation::checkpoint()?;
        result
    }
}

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    for (name, path, method, doc) in [
        (
            "num_univalent_schwarzian_checked",
            "crate::kernels::univalent::schwarzian_derivative_numeric",
            Method::Schwarzian,
            "Estimates the measured Schwarzian derivative from one 64-sample Cauchy circle. Pass center, radius, a positive callback-value scale, total callback budget, then the callback.",
        ),
        (
            "num_univalent_nehari_sampled",
            "crate::kernels::univalent::check_nehari_univalence",
            Method::Nehari,
            "Checks Nehari’s inequality on a finite angular/radial grid in the normalized unit disk. Pass sample radii, angles per radius, positive Cauchy radius, callback budget, then the callback.",
        ),
        (
            "num_univalent_starlikeness_radius_sampled",
            "crate::kernels::univalent::radius_of_starlikeness_numeric",
            Method::Starlikeness,
            "Bisects a finite angular diagnostic for starlikeness. Pass maximum normalized radius, angular count, callback budget, the map and its first derivative.",
        ),
        (
            "num_univalent_convexity_radius_sampled",
            "crate::kernels::univalent::radius_of_convexity_numeric",
            Method::Convexity,
            "Bisects a finite angular diagnostic for convexity. Pass maximum normalized radius, angular count, callback budget, the first and second derivatives.",
        ),
    ] {
        register(i, name, path, doc, 5, method)?;
    }
    Ok(())
}
