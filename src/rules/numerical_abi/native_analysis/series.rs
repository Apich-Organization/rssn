//! Measured estimates from a supplied prefix of an infinite series.
use super::CallCtx;
use super::Installer;
use super::OpError;
use super::RuleError;
use super::native_object_api as api;
use crate::kernel::signature::DimExpr;
use crate::kernels::indefinite_sum::PolynomialAntidifference;
use crate::kernels::series::sum_terms_to_infinity;
use api::Dynamic;
use api::Handle;
use api::Input;
use api::Maybe;
use api::MeasuredVector;
use api::NativeOutput;
use api::NativeType;
use api::NewObject;
use api::object_binding;

const MAX_TERMS: usize = 262_143;
const MAX_ANTIDIFFERENCE_COEFFICIENTS: usize = 128;

object_binding!(
    sum_terms,
    "num_analysis_series_sum_terms",
    "crate::kernels::series::sum_terms_to_infinity",
    "Estimates an infinite sum from a finite measured prefix, stopping at the requested absolute tolerance, the end of the prefix, or the term budget. The record's value and heuristic error have the term dimension. Missing estimates are None, and terms reports the number consumed.",
    false,
    (terms: MeasuredVector<MAX_TERMS>, tolerance: Input, max_terms: usize) -> Dynamic,
    |ctx, terms, tolerance, max_terms| {
        api::same_dimension(&terms.dimension, &api::dimension(&tolerance.0)?)?;
        let tolerance = crate::rules::numerical_abi::numeric(&tolerance.0, 1)?;
        if tolerance < 0.0 || max_terms > MAX_TERMS {
            return Err(api::error("series needs a nonnegative tolerance and at most 262143 terms"));
        }
        let result = sum_terms_to_infinity(terms.values, tolerance, max_terms);
        let measured = |value: f64| {
            let payload = value.is_finite().then(|| api::quantity(value, terms.dimension.clone())).transpose()?;
            Maybe(payload.map(Dynamic)).encode(ctx)
        };
        Ok(Dynamic(api::record(vec![
            ("value", measured(result.value)?),
            ("error", measured(result.error)?),
            ("terms", crate::kernel::Value::int(result.terms)),
        ])))
    }
);

struct MeasuredAntidifference {
    polynomial: PolynomialAntidifference,
    origin: f64,
    step: f64,
}

impl NativeType for MeasuredAntidifference {
    const KIND: &'static str = "rssn/native/analysis/polynomial-antidifference/v1";
    const MAX_BYTES: usize = 2048;

    fn validate(&self) -> Result<(), OpError> {
        if self.polynomial.coefficients_are_finite() {
            Ok(())
        } else {
            Err(api::error(
                "polynomial antidifference coefficients exceed Float range",
            ))
        }
    }
}

fn coordinate(
    value: &Input,
    expected: &DimExpr,
    index: usize,
) -> Result<f64, OpError> {
    api::same_dimension(&api::dimension(&value.0)?, expected)?;
    super::numeric(&value.0, index)
}

fn scaled_delta(
    left: f64,
    right: f64,
    step: f64,
) -> Result<f64, OpError> {
    let difference = left - right;
    if difference.is_finite() {
        api::scale_product(difference, &[], &[step])
    } else {
        let left = api::scale_product(left, &[], &[step])?;
        let right = api::scale_product(right, &[], &[step])?;
        api::scale_product(left - right, &[], &[])
    }
}

fn antidifference_dimensions(
    plan: &Handle<MeasuredAntidifference>,
    ctx: &CallCtx<'_>,
) -> Result<[DimExpr; 2], OpError> {
    plan.dimensions(ctx)?.try_into().map_err(|_| {
        api::error("polynomial antidifference lacks its coordinate and value dimensions")
    })
}

object_binding!(
    antidifference_new,
    "num_analysis_polynomial_antidifference_new",
    "crate::kernels::indefinite_sum::PolynomialAntidifference::new",
    "Prepare the discrete antidifference for up to 128 ascending polynomial coefficients. Use normalized coordinate u = (x - origin) / step. The origin and positive step have the same dimension, and all coefficients have the common value dimension. The prepared function F is zero at the origin and F(x + step) - F(x) is the polynomial evaluated at u.",
    false,
    (coefficients: MeasuredVector<{MAX_ANTIDIFFERENCE_COEFFICIENTS}>, origin: Input, step: Input) -> NewObject<MeasuredAntidifference>,
    |_, coefficients, origin, step| {
        let (step, axis) = api::reference(&step.0)?;
        let origin = coordinate(&origin, &axis, 1)?;
        let polynomial = PolynomialAntidifference::new(&coefficients.values, 0.0);
        Ok(NewObject::new(MeasuredAntidifference { polynomial, origin, step }, vec![axis, coefficients.dimension]))
    }
);

object_binding!(
    antidifference_eval,
    "num_analysis_polynomial_antidifference_eval",
    "crate::kernels::indefinite_sum::PolynomialAntidifference::evaluate",
    "Evaluate the prepared polynomial antidifference at a measured coordinate. Its result keeps the polynomial coefficient dimension.",
    false,
    (plan: Handle<MeasuredAntidifference>, point: Input) -> Dynamic,
    |ctx, plan, point| {
        let [axis, values] = antidifference_dimensions(&plan, ctx)?;
        let point = coordinate(&point, &axis, 1)?;
        let value = plan.with(ctx, |plan| {
            let offset = scaled_delta(point, plan.origin, plan.step)?;
            Ok(plan.polynomial.evaluate(offset))
        })?;
        Ok(Dynamic(api::quantity(value, values)?))
    }
);

object_binding!(
    antidifference_delta,
    "num_analysis_polynomial_antidifference_difference",
    "crate::kernels::indefinite_sum::PolynomialAntidifference::difference",
    "Evaluate F(to) minus F(from) using a divided difference and the measured interval width. Both endpoints share the prepared coordinate dimension. The result keeps the polynomial coefficient dimension.",
    false,
    (plan: Handle<MeasuredAntidifference>, from: Input, to: Input) -> Dynamic,
    |ctx, plan, from, to| {
        let [axis, values] = antidifference_dimensions(&plan, ctx)?;
        let from = coordinate(&from, &axis, 1)?;
        let to = coordinate(&to, &axis, 2)?;
        let value = plan.with(ctx, |plan| {
            let start = scaled_delta(from, plan.origin, plan.step)?;
            let end = scaled_delta(to, plan.origin, plan.step)?;
            let width = scaled_delta(to, from, plan.step)?;
            Ok(plan.polynomial.difference_of_offsets(start, end, width))
        })?;
        Ok(Dynamic(api::quantity(value, values)?))
    }
);

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &sum_terms::BINDING,
            &antidifference_new::BINDING,
            &antidifference_eval::BINDING,
            &antidifference_delta::BINDING,
        ],
        "numerical_vector",
    )
}
