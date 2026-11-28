//! Bounded, measured derivatives on a normalized Cauchy circle.

use crate::kernel::signature::NatExpr;
use crate::kernel::value::TensorValue;
use crate::kernels::complex::cauchy_fft_derivatives;
use crate::kernels::complex_roots::validate_circle;
use crate::rules::numerical_abi::codecs::Arg;

use super::super::super::dimensional_float;
use super::super::arity;
use super::super::domain;
use super::super::record;
use super::CallCtx;
use super::Callbacks;
use super::Complex64;
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
use super::finite;
use super::function;
use super::real;
use super::register;
use super::scale;
use super::typed;
use super::value;
use super::var;

struct Derivatives;

impl EvalDyn for Derivatives {
    fn signature(&self) -> Signature {
        let d = var("D");
        let e = var("E");
        let point = complex_type(d.clone());
        let output = complex_type(e.clone());
        typed(
            vec![
                point.clone(),
                real(d),
                real(e.clone()),
                real(e.clone()),
                TypeExpr::Int,
                TypeExpr::Int,
                TypeExpr::Int,
                function(vec![point], output.clone()),
            ],
            TypeExpr::Record(vec![
                (
                    "derivatives".into(),
                    TypeExpr::Tensor {
                        rank: NatExpr::Lit(1),
                        shape: None,
                        elem: Box::new(output),
                    },
                ),
                ("sample_error".into(), real(e)),
                ("converged".into(), TypeExpr::Bool),
                ("samples".into(), TypeExpr::Int),
                ("evaluations".into(), TypeExpr::Int),
            ]),
            &["D", "E"],
        )
    }

    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let expected = self.signature();
        arity(args, expected.params.len())?;
        let order = usize::decode(&args[4], 4)?;
        let samples = count(&args[5], 5)?;
        let budget = count(&args[6], 6)?;
        if samples < 8 || order >= samples {
            return Err(domain(
                "the Cauchy grid needs at least eight samples and more samples than the derivative order",
            ));
        }
        if order > 32 || samples > 4_096 || budget < 3 * samples {
            return Err(OpError::BudgetExceeded);
        }
        let callbacks = Callbacks::new(ctx, args, &expected, 7, budget)?;
        let d = dimension(&callbacks.signatures[0].params[0])?;
        let e = dimension(&callbacks.signatures[0].ret)?;
        let center = decode(&args[0], &d, 0, ctx)?;
        let radius = scale(&args[1], &d, 1, ctx)?;
        let amplitude = scale(&args[2], &e, 2, ctx)?;
        let tolerance = scale(&args[3], &e, 3, ctx)?;
        validate_circle(center, radius, samples).map_err(domain)?;
        let bad = Complex64::new(f64::NAN, f64::NAN);
        let sample = |u| {
            let point = center + u * radius;
            callbacks
                .call(0, &[value(point, d.clone())], |result| {
                    let sample = decode(&result, &e, 7, ctx)? / amplitude;
                    Ok(Complex64::new(finite(sample.re)?, finite(sample.im)?))
                })
                .unwrap_or(bad)
        };
        let origin = Complex64::new(0.0, 0.0);
        let coarse = cauchy_fft_derivatives(sample, origin, 1.0, order, samples);
        let fine = cauchy_fft_derivatives(sample, origin, 1.0, order, 2 * samples);
        let evaluations = callbacks.used.get();
        callbacks.finish()?;
        crate::cancellation::checkpoint()?;
        let error = coarse
            .iter()
            .zip(&fine)
            .try_fold(0.0_f64, |max, (&a, &b)| {
                finite((a - b).norm() * amplitude).map(|error| max.max(error))
            })?;
        let output = fine
            .into_iter()
            .map(|z| {
                let z = z * amplitude;
                Ok(value(
                    Complex64::new(finite(z.re)?, finite(z.im)?),
                    e.clone(),
                ))
            })
            .collect::<Result<Vec<_>, OpError>>()?;
        let derivatives = Value::Tensor(TensorValue::of(
            vec![output.len()],
            output,
            complex_type(e.clone()),
        ));
        Ok(record(vec![
            ("derivatives", derivatives),
            ("sample_error", dimensional_float(error, e)),
            ("converged", Value::Bool(error <= tolerance)),
            ("samples", Value::int(2 * samples)),
            ("evaluations", Value::int(evaluations)),
        ]))
    }
}

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    register(
        i,
        "num_complex_cauchy_fft_derivatives_checked",
        "crate::kernels::complex::cauchy_fft_derivatives",
        "Returns Cauchy FFT derivatives in the normalized local coordinate (z-center)/radius. All entries share the callback value dimension. A doubled grid estimates the largest absolute difference, and convergence compares this estimate with the requested absolute tolerance. Reports the refined sample count and actual callback evaluations.",
        8,
        Derivatives,
    )
}
