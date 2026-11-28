//! Laurent coefficients represented in one dimensionless local coordinate.

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
use super::function;
use super::real;
use super::record;
use super::register;
use super::scale;
use super::typed;
use super::var;
use crate::kernel::signature::NatExpr;
use crate::kernel::value::TensorValue;
use crate::rules::numerical_abi::codecs::Arg;
use num_complex::Complex64;

#[derive(Clone, Copy)]
struct Laurent;

impl EvalDyn for Laurent {
    fn signature(&self) -> Signature {
        let axis = var("D");
        let amplitude = var("E");
        let point = complex_type(axis.clone());
        let result = complex_type(amplitude);
        typed(
            vec![
                point.clone(),
                real(axis.clone()),
                TypeExpr::Int,
                TypeExpr::Int,
                TypeExpr::Int,
                TypeExpr::Int,
                function(vec![point], result.clone()),
            ],
            TypeExpr::Record(vec![
                (
                    "coefficients".into(),
                    TypeExpr::Tensor {
                        rank: NatExpr::Lit(1),
                        shape: None,
                        elem: Box::new(result),
                    },
                ),
                ("first_power".into(), TypeExpr::Int),
                ("radius".into(), real(axis)),
                ("evaluations".into(), TypeExpr::Int),
            ]),
            &["D", "E"],
        )
    }

    #[allow(
        clippy::indexing_slicing,
        clippy::arithmetic_side_effects,
        reason = "the fixed arity and bounded FFT order validate indices and sampling arithmetic"
    )]
    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        arity(args, 7)?;
        let negative = usize::decode(&args[2], 2)?;
        let positive = usize::decode(&args[3], 3)?;
        let requested = usize::decode(&args[4], 4)?;
        let budget = count(&args[5], 5)?;
        let min_samples = negative.saturating_add(positive).saturating_add(2).max(8);
        let samples = requested.max(min_samples);
        if negative > 256 || positive > 256 || samples > 65_536 || samples > budget {
            return Err(OpError::BudgetExceeded);
        }
        let callbacks = Callbacks::new(ctx, args, &self.signature(), 6, budget)?;
        let [signature] = callbacks.signatures.as_slice() else {
            return Err(domain("Laurent expansion expects one analytic callback"));
        };
        let [input] = signature.params.as_slice() else {
            return Err(domain("Laurent callback expects one coordinate"));
        };
        let axis = dimension(input)?;
        let amplitude = dimension(&signature.ret)?;
        let center = decode(&args[0], &axis, 0, ctx)?;
        let radius = scale(&args[1], &axis, 1, ctx)?;
        crate::kernels::complex_roots::validate_circle(center, radius, samples).map_err(domain)?;
        let output = crate::kernels::analytic_continuation::laurent_coeffs_fft(
            |unit| {
                callbacks
                    .call(
                        0,
                        &[value(center + radius * unit, axis.clone())],
                        |response| decode(&response, &amplitude, 6, ctx),
                    )
                    .unwrap_or(Complex64::new(f64::NAN, f64::NAN))
            },
            Complex64::new(0.0, 0.0),
            1.0,
            negative,
            positive,
            samples,
        );
        let evaluations = callbacks.used.get();
        callbacks.finish()?;
        crate::cancellation::checkpoint()?;
        let Some((first, _)) = output.first() else {
            return Err(OpError::NotReduced(
                "the sampled Laurent transform is nonfinite".into(),
            ));
        };
        let first = *first;
        let data = output
            .into_iter()
            .map(|(_, coefficient)| value(coefficient, amplitude.clone()))
            .collect::<Vec<_>>();
        let coefficients = Value::Tensor(TensorValue::of(
            vec![data.len()],
            data,
            complex_type(amplitude),
        ));
        let radius = if axis == DimExpr::Dimless {
            Value::Float(radius)
        } else {
            Value::Dyn(crate::kernel::value::DynValue {
                value: Box::new(Value::Float(radius)),
                dim: axis,
            })
        };
        Ok(record(vec![
            ("coefficients", coefficients),
            ("first_power", Value::int(first)),
            ("radius", radius),
            ("evaluations", Value::int(evaluations)),
        ]))
    }
}

#[allow(
    clippy::single_call_fn,
    reason = "registration entrypoint for the Laurent FFT operator"
)]
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    register(
        i,
        "num_continuation_laurent_fft",
        "crate::kernels::analytic_continuation::laurent_coeffs_fft",
        "Samples an analytic annulus along a measured circle and returns Laurent coefficients in the dimensionless coordinate (z-center)/radius. Coefficients share the callback amplitude units; first_power is the first exponent and evaluations counts actual callback calls.",
        7,
        Laurent,
    )
}
