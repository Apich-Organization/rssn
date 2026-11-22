//! Survival analysis.
//!
//! | operator | value |
//! |---|---|
//! | `kaplan_meier(time, event, conf)` | `list(times, at_risk, events, survival, variance, lower, upper)`: the product-limit estimator with Greenwood variance and a complementary-log-log confidence band |
//! | `nelson_aalen(time, event)` | `list(times, hazard, cumulative_hazard)`: the Nelson–Aalen estimator of the cumulative hazard |
//! | `log_rank(time_a, event_a, time_b, event_b)` | `list(chi2, p_value, observed, expected)`: the two-sample log-rank test |
//! | `cox_fit(x, time, event, ties, tol)` | `list(beta, se, z, p, loglik, concordance)`: the Cox proportional-hazards fit |
//! | `cox_predict(x, beta)` | the relative hazard `exp(x . beta)` of each row |
//!
//! `event` is a list of `1` (the event was observed) and `0` (the observation
//! was censored); `x` is a list of covariate rows without an intercept. `ties`
//! is `breslow` or `efron` and decides how tied event times enter the partial
//! likelihood.
//!
//! A data set the estimator cannot use — no observed event, mismatched
//! lengths, an event flag that is neither `0` nor `1` — is a
//! [`DiagKind::Domain`] finding, not a curve.

use super::bridge::Data;
use super::bridge::Entry;
use super::bridge::bad;
use super::bridge::require;
use super::bridge::rows;
use super::bridge::scalar;
use super::bridge::vector;
use super::conic::numbers;
use crate::kernel::diag::DiagKind;
use crate::kernels::survival::Ties;

/// Every operator of this family.
pub(super) const ENTRIES: &[Entry] = &[
    Entry {
        name: "kaplan_meier",
        params: "vvd",
        eval: kaplan_meier,
        doc: "`kaplan_meier(time, event, conf)`: `list(times, at_risk, events, survival, variance, lower, upper)` of the product-limit estimator, with `event` a list of 1 (observed) and 0 (censored) and a `conf` confidence level.",
    },
    Entry {
        name: "nelson_aalen",
        params: "vv",
        eval: nelson_aalen,
        doc: "`nelson_aalen(time, event)`: `list(times, hazard, cumulative_hazard)` of the Nelson–Aalen estimator.",
    },
    Entry {
        name: "log_rank",
        params: "vvvv",
        eval: log_rank,
        doc: "`log_rank(time_a, event_a, time_b, event_b)`: `list(chi2, p_value, observed, expected)` of the two-sample log-rank test.",
    },
    Entry {
        name: "cox_fit",
        params: "tvvdd",
        eval: cox_fit,
        doc: "`cox_fit(x, time, event, ties, tol)`: `list(beta, se, z, p, loglik, concordance)` of the Cox proportional-hazards model; `x` is a list of covariate rows without an intercept and `ties` is `breslow` or `efron`.",
    },
    Entry {
        name: "cox_predict",
        params: "tv",
        eval: cox_predict,
        doc: "`cox_predict(x, beta)`: the relative hazard `exp(x . beta)` of each covariate row.",
    },
];

/// The event flags of `at`, rejecting anything that is not `0` or `1`.
fn events(
    args: &[Data],
    at: usize,
) -> Option<Vec<bool>> {
    let values = vector(args, at)?;
    let mut flags = Vec::with_capacity(values.len());
    for &value in values {
        if value == 1.0 {
            flags.push(true);
        } else if value == 0.0 {
            flags.push(false);
        } else {
            return bad(
                at,
                DiagKind::Domain,
                format_args!("an event flag must be 1 (observed) or 0 (censored), got {value}"),
            );
        }
    }
    Some(flags)
}

/// The observation times of `at`.
fn times(
    args: &[Data],
    at: usize,
) -> Option<&[f64]> {
    let values = vector(args, at)?;
    require(
        !values.is_empty(),
        Some(at),
        DiagKind::Domain,
        format_args!("the sample must not be empty"),
    )?;
    Some(values)
}

/// A list of counts as a numeric list.
fn counts(values: Vec<usize>) -> Data {
    Data::vector(
        values
            .into_iter()
            .map(|value| {
                #[allow(clippy::cast_precision_loss)]
                let scaled = value as f64;
                scaled
            })
            .collect(),
    )
}

fn kaplan_meier(args: &[Data]) -> Option<Data> {
    let time = times(args, 0)?;
    let event = events(args, 1)?;
    let conf = scalar(args, 2)?;
    require(
        conf > 0.0 && conf < 1.0,
        Some(2),
        DiagKind::Domain,
        format_args!("the confidence level must lie strictly between 0 and 1, got {conf}"),
    )?;
    require(
        time.len() == event.len(),
        Some(1),
        DiagKind::Shape,
        format_args!(
            "the times have {} entries but the flags have {}",
            time.len(),
            event.len()
        ),
    )?;
    match crate::kernels::survival::kaplan_meier(time, &event, conf) {
        | Ok(estimate) => {
            Some(Data::Tree(vec![
                Data::vector(estimate.times),
                counts(estimate.at_risk),
                counts(estimate.events),
                Data::vector(estimate.survival),
                Data::vector(estimate.variance),
                Data::vector(estimate.lower),
                Data::vector(estimate.upper),
            ]))
        },
        | Err(message) => bad(0, DiagKind::Domain, format_args!("{message}")),
    }
}

fn nelson_aalen(args: &[Data]) -> Option<Data> {
    let time = times(args, 0)?;
    let event = events(args, 1)?;
    require(
        time.len() == event.len(),
        Some(1),
        DiagKind::Shape,
        format_args!(
            "the times have {} entries but the flags have {}",
            time.len(),
            event.len()
        ),
    )?;
    match crate::kernels::survival::nelson_aalen(time, &event) {
        | Ok(estimate) => {
            Some(Data::Tree(vec![
                Data::vector(estimate.times),
                Data::vector(estimate.hazard),
                Data::vector(estimate.cumulative_hazard),
            ]))
        },
        | Err(message) => bad(0, DiagKind::Domain, format_args!("{message}")),
    }
}

fn log_rank(args: &[Data]) -> Option<Data> {
    let time_a = times(args, 0)?;
    let event_a = events(args, 1)?;
    let time_b = times(args, 2)?;
    let event_b = events(args, 3)?;
    match crate::kernels::survival::log_rank(time_a, &event_a, time_b, &event_b) {
        | Ok(result) => {
            Some(Data::Tree(vec![
                Data::Scalar(result.chi2),
                Data::Scalar(result.p_value),
                Data::vector(result.observed.to_vec()),
                Data::vector(result.expected.to_vec()),
            ]))
        },
        | Err(message) => bad(0, DiagKind::Domain, format_args!("{message}")),
    }
}

fn cox_fit(args: &[Data]) -> Option<Data> {
    let Some(x) = args.get(0).and_then(rows) else {
        return bad(
            0,
            DiagKind::Shape,
            format_args!("expected a list of covariate rows"),
        );
    };
    let time = times(args, 1)?;
    let event = events(args, 2)?;
    let ties_text = super::bridge::name(args, 3)?;
    let ties = match ties_text {
        | "breslow" => Ties::Breslow,
        | "efron" => Ties::Efron,
        | _ => {
            return bad(
                3,
                DiagKind::Domain,
                format_args!("`{ties_text}` is not a tie rule (breslow, efron)"),
            );
        },
    };
    let tolerance = scalar(args, 4)?;
    require(
        x.len() == time.len() && time.len() == event.len(),
        Some(0),
        DiagKind::Shape,
        format_args!(
            "the covariates have {} rows, the times {} and the flags {}",
            x.len(),
            time.len(),
            event.len()
        ),
    )?;
    require(
        tolerance > 0.0 && tolerance.is_finite(),
        Some(4),
        DiagKind::Domain,
        format_args!("the tolerance must be positive and finite, got {tolerance}"),
    )?;
    match crate::kernels::survival::cox_fit(&x, time, &event, ties, 100, tolerance) {
        | Ok(fit) => {
            Some(Data::Tree(vec![
                Data::vector(fit.beta),
                Data::vector(fit.se),
                Data::vector(fit.z),
                Data::vector(fit.p),
                Data::Scalar(fit.loglik),
                Data::Scalar(fit.concordance),
            ]))
        },
        | Err(message) => bad(0, DiagKind::Domain, format_args!("{message}")),
    }
}

fn cox_predict(args: &[Data]) -> Option<Data> {
    let Some(x) = args.get(0).and_then(rows) else {
        return bad(
            0,
            DiagKind::Shape,
            format_args!("expected a list of covariate rows"),
        );
    };
    let beta = numbers(args, 1)?;
    let Some(width) = x.first().map(Vec::len) else {
        return bad(0, DiagKind::Shape, format_args!("the covariates are empty"));
    };
    require(
        beta.len() == width,
        Some(1),
        DiagKind::Shape,
        format_args!(
            "the coefficients need {width} entries (one per covariate), got {}",
            beta.len()
        ),
    )?;
    let mut hazards = Vec::with_capacity(x.len());
    for row in &x {
        let mut score = 0.0;
        for (value, coefficient) in row.iter().zip(beta.iter()) {
            score += value * coefficient;
        }
        hazards.push(score.exp());
    }
    super::bridge::vec_out(hazards)
}
