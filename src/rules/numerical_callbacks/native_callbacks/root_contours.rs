//! Bounded sampled moments and winding for a measured complex circle.

#[path = "root_iterations.rs"]
mod iterations;

use super::super::arity;
use super::super::domain;
use super::super::record;
use super::CallCtx;
use super::Callbacks;
use super::Complex64;
use super::DimExpr;
use super::EvalDyn;
use super::Installer;
use super::OpError;
use super::RuleError;
use super::Signature;
use super::TypeExpr;
use super::Value;
use super::complex_type;
use super::count;
use super::decode;
use super::dimension;
use super::equal;
use super::finite;
use super::function;
use super::ratio;
use super::real;
use super::register;
use super::scale;
use super::typed;
use super::value;
use super::var;
use crate::kernel::signature::NatExpr;
use crate::kernel::value::TensorValue;
use crate::kernels::complex_roots as native;
use crate::rules::numerical_abi::codecs::Arg;

#[derive(Clone, Copy)]
enum SampleMethod {
    Moments,
    Winding,
    Rouche,
}

const MAX_SAMPLES: usize = 8_192;
const MAX_MOMENT: usize = 32;

struct Circle<'a, 'ctx> {
    ctx: &'a CallCtx<'ctx>,
    callbacks: Callbacks<'a, 'ctx>,
    center: Complex64,
    radius: f64,
    axis: DimExpr,
    output: [(DimExpr, f64); 2],
}

impl<'a, 'ctx> Circle<'a, 'ctx> {
    fn new(
        args: &[Value],
        ctx: &'a CallCtx<'ctx>,
        expected: &Signature,
        samples: usize,
        budget: usize,
        method: SampleMethod,
    ) -> Result<Self, OpError> {
        let callbacks = Callbacks::new(ctx, args, expected, method.start(), budget)?;
        let axis = dimension(&callbacks.signatures[0].params[0])?;
        let amplitude_dim = dimension(&callbacks.signatures[0].ret)?;
        equal(&dimension(&callbacks.signatures[1].params[0])?, &axis)?;
        let center = decode(&args[0], &axis, 0, ctx)?;
        let radius = scale(&args[1], &axis, 1, ctx)?;
        let amplitude = scale(&args[2], &amplitude_dim, 2, ctx)?;
        native::validate_circle(center, radius, samples).map_err(domain)?;
        let (second_dim, second_scale) = if matches!(method, SampleMethod::Rouche) {
            (amplitude_dim.clone(), amplitude)
        } else {
            (ratio(&amplitude_dim, &axis)?, finite(amplitude / radius)?)
        };
        equal(&dimension(&callbacks.signatures[1].ret)?, &second_dim)?;
        if second_scale == 0.0 {
            return Err(domain("the normalized callback scale underflowed"));
        }
        Ok(Self {
            ctx,
            callbacks,
            center,
            radius,
            axis,
            output: [(amplitude_dim, amplitude), (second_dim, second_scale)],
        })
    }

    fn sample(
        &self,
        which: usize,
        point: Complex64,
    ) -> Complex64 {
        let point = self.center + self.radius * point;
        let bad_sample = Complex64::new(f64::NAN, f64::NAN);
        if !point.is_finite() {
            return bad_sample;
        }
        let (dimension, scale) = &self.output[which];
        self.callbacks
            .call(which, &[value(point, self.axis.clone())], |v| {
                let result = decode(&v, dimension, which, self.ctx)? / *scale;
                Ok(Complex64::new(finite(result.re)?, finite(result.im)?))
            })
            .unwrap_or(bad_sample)
    }
}

impl SampleMethod {
    const fn start(self) -> usize {
        if matches!(self, Self::Moments) {
            6
        } else {
            5
        }
    }

    fn sampled(
        self,
        circle: &Circle<'_, '_>,
        n: usize,
        order: usize,
    ) -> Result<Vec<(&'static str, Value)>, String> {
        let f = |z| circle.sample(0, z);
        let second = |z| circle.sample(1, z);
        let origin = Complex64::new(0.0, 0.0);
        match self {
            | Self::Moments => {
                let coarse = native::argument_principle_moments(f, second, origin, 1.0, order, n)?;
                let fine =
                    native::argument_principle_moments(f, second, origin, 1.0, order, 2 * n)?;
                let delta = coarse
                    .iter()
                    .zip(&fine)
                    .map(|(a, b)| (a - b).norm())
                    .fold(0.0_f64, f64::max);
                if !delta.is_finite() {
                    return Err("moment grid difference is nonfinite".into());
                }
                Ok(vec![
                    (
                        "moments",
                        Value::Tensor(TensorValue::of(
                            vec![fine.len()],
                            fine.into_iter()
                                .map(|z| value(z, DimExpr::Dimless))
                                .collect(),
                            TypeExpr::Complex,
                        )),
                    ),
                    ("sample_delta", Value::Float(delta)),
                ])
            },
            | Self::Winding => {
                let coarse = native::argument_principle_winding(f, second, origin, 1.0, n)?;
                let fine = native::argument_principle_winding(f, second, origin, 1.0, 2 * n)?;
                if coarse != fine {
                    return Err("coarse and refined circles disagree on winding".into());
                }
                Ok(vec![("winding", Value::int(fine))])
            },
            | Self::Rouche => {
                let checked = native::rouche_check_circle(f, second, origin, 1.0, n)?;
                Ok(vec![("satisfied_at_samples", Value::Bool(checked))])
            },
        }
    }
}

impl EvalDyn for SampleMethod {
    fn signature(&self) -> Signature {
        let d = var("D");
        let e = var("E");
        let point = complex_type(d.clone());
        let second = if matches!(self, Self::Rouche) {
            e.clone()
        } else {
            DimExpr::Mul(vec![e.clone(), DimExpr::Inverse(Box::new(d.clone()))])
        };
        let mut params = vec![point.clone(), real(d), real(e.clone()), TypeExpr::Int];
        if matches!(self, Self::Moments) {
            params.push(TypeExpr::Int);
        }
        params.extend([
            TypeExpr::Int,
            function(vec![point.clone()], complex_type(e)),
            function(vec![point], complex_type(second)),
        ]);
        let mut fields = match self {
            | Self::Moments => {
                vec![
                    (
                        "moments".into(),
                        TypeExpr::Tensor {
                            rank: NatExpr::Lit(1),
                            shape: None,
                            elem: Box::new(TypeExpr::Complex),
                        },
                    ),
                    ("sample_delta".into(), TypeExpr::Float),
                ]
            },
            | Self::Winding => vec![("winding".into(), TypeExpr::Int)],
            | Self::Rouche => vec![("satisfied_at_samples".into(), TypeExpr::Bool)],
        };
        fields.extend([
            ("samples".into(), TypeExpr::Int),
            ("evaluations".into(), TypeExpr::Int),
        ]);
        typed(params, TypeExpr::Record(fields), &["D", "E"])
    }

    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let expected = self.signature();
        arity(args, expected.params.len())?;
        let start = self.start();
        let n = count(&args[3], 3)?;
        if n < 8 {
            return Err(domain("the circle requires at least eight samples"));
        }
        let refined = if matches!(self, Self::Rouche) {
            n
        } else {
            n.saturating_mul(2)
        };
        let order = if matches!(self, Self::Moments) {
            usize::decode(&args[4], 4)?
        } else {
            0
        };
        let work = if matches!(self, Self::Rouche) {
            n.saturating_mul(2)
        } else {
            n.saturating_mul(6)
        };
        let budget = count(&args[start - 1], start - 1)?;
        if refined > MAX_SAMPLES || order > MAX_MOMENT || work > budget {
            return Err(OpError::BudgetExceeded);
        }
        let circle = Circle::new(args, ctx, &expected, n, budget, *self)?;
        let result = self.sampled(&circle, n, order);
        let evaluations = circle.callbacks.used.get();
        circle.callbacks.finish()?;
        crate::cancellation::checkpoint()?;
        let mut fields = result.map_err(domain)?;
        fields.extend([
            ("samples", Value::int(refined)),
            ("evaluations", Value::int(evaluations)),
        ]);
        Ok(record(fields))
    }
}

#[allow(
    clippy::single_call_fn,
    reason = "registration entrypoint for the measured root domain"
)]
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    iterations::install(i)?;
    register(
        i,
        "num_complex_root_circle_moments",
        "crate::kernels::complex_roots::argument_principle_moments",
        "Returns sampled logarithmic-derivative moments of (z-center)/radius, the largest coarse/fine difference, and sample/evaluation counts. Takes a measured center, radius, callback value scale, initial sample count, maximum moment, host budget, analytic function and derivative.",
        8,
        SampleMethod::Moments,
    )?;
    register(
        i,
        "num_complex_root_circle_winding",
        "crate::kernels::complex_roots::argument_principle_winding",
        "Returns a sampled winding estimate and actual counts if the coarse and refined complex integrals approach the same integer. Takes the measured circle, value scale, initial samples, host budget, and analytic or meromorphic function and derivative.",
        7,
        SampleMethod::Winding,
    )?;
    register(
        i,
        "num_complex_root_circle_rouche_samples",
        "crate::kernels::complex_roots::rouche_check_circle",
        "Reports |f|>|g| at each sampled point of a measured circle and returns actual counts. Takes the circle, shared callback value scale, sample count, host budget and two analytic functions. An entire-circle inequality requires a separate enclosure.",
        7,
        SampleMethod::Rouche,
    )
}
