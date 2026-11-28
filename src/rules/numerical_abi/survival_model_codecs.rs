// Portable Cox and Kaplan–Meier records shared by fitting and prediction.
use crate::kernels::survival::{CoxFit, KmEstimate, Ties};

record!(CoxFit {
    beta: Vec<f64>,
    se: Vec<f64>,
    z: Vec<f64>,
    p: Vec<f64>,
    loglik: f64,
    iterations: usize,
    converged: bool,
    concordance: f64,
});

record!(KmEstimate {
    times: Vec<f64>,
    at_risk: Vec<usize>,
    events: Vec<usize>,
    survival: Vec<f64>,
    variance: Vec<f64>,
    lower: Vec<f64>,
    upper: Vec<f64>,
});

impl Arg for CoxFit {
    fn ty() -> TypeExpr {
        <Self as Output>::ty()
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let record = native_fields(
            value,
            &[
                "beta",
                "se",
                "z",
                "p",
                "loglik",
                "iterations",
                "converged",
                "concordance",
            ],
        )?;
        let field = |name| native_field(record, name);
        Ok(Self {
            beta: Vec::<f64>::decode(field("beta")?, index)?,
            se: Vec::<f64>::decode(field("se")?, index)?,
            z: Vec::<f64>::decode(field("z")?, index)?,
            p: Vec::<f64>::decode(field("p")?, index)?,
            loglik: f64::decode(field("loglik")?, index)?,
            iterations: usize::decode(field("iterations")?, index)?,
            converged: bool::decode(field("converged")?, index)?,
            concordance: f64::decode(field("concordance")?, index)?,
        })
    }
}

impl Arg for KmEstimate {
    fn ty() -> TypeExpr {
        <Self as Output>::ty()
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let record = native_fields(
            value,
            &[
                "times", "at_risk", "events", "survival", "variance", "lower", "upper",
            ],
        )?;
        let field = |name| native_field(record, name);
        let estimate = Self {
            times: Vec::<f64>::decode(field("times")?, index)?,
            at_risk: Vec::<usize>::decode(field("at_risk")?, index)?,
            events: Vec::<usize>::decode(field("events")?, index)?,
            survival: Vec::<f64>::decode(field("survival")?, index)?,
            variance: Vec::<f64>::decode(field("variance")?, index)?,
            lower: Vec::<f64>::decode(field("lower")?, index)?,
            upper: Vec::<f64>::decode(field("upper")?, index)?,
        };
        let n = estimate.times.len();
        let matches = [
            estimate.at_risk.len(),
            estimate.events.len(),
            estimate.survival.len(),
            estimate.variance.len(),
            estimate.lower.len(),
            estimate.upper.len(),
        ]
        .iter()
        .all(|&length| length == n);
        if !matches
            || estimate.times.iter().any(|time| !time.is_finite())
            || estimate.times.windows(2).any(|pair| pair[0] >= pair[1])
            || estimate
                .survival
                .iter()
                .any(|value| !(0.0..=1.0).contains(value))
            || estimate.survival.windows(2).any(|pair| pair[0] < pair[1])
        {
            return Err(OpError::Domain(
                "Kaplan–Meier record requires equal-length columns, strictly increasing times and a nonincreasing survival probability".into(),
            ));
        }
        Ok(estimate)
    }
}

binding!(
    survival_fit_cox, "num_survival_fit_cox", "numerical_stats",
    "Fits a Cox proportional-hazards model for observed times, event flags and finite dimensionless covariate rows without an intercept. Select breslow or efron for tied events. Returns portable coefficients and fit diagnostics.",
    "crate::kernels::survival::cox::cox_fit", DimensionPolicy::UnitlessResult,
    (time: Vec<f64>, event: Vec<bool>, rows: Vec<Vec<f64>>, ties: String, max_iterations: usize, tolerance: f64)
        -> Result<super::CoxFit, OpError>,
    |time: Vec<f64>, event: Vec<bool>, rows: Vec<Vec<f64>>, ties: String, max_iterations: usize, tolerance: f64| {
        super::fit_cox_model(&time, &event, &rows, &ties, max_iterations, tolerance)
    }
);

binding!(
    survival_predict_cox, "num_survival_predict_cox", "numerical_stats",
    "Predicts finite relative hazards exp(row*beta) from a portable fitted Cox model and covariate rows with the fitted column count.",
    "crate::kernels::survival::cox::cox_predict", DimensionPolicy::Dimensionless,
    (fit: super::CoxFit, rows: Vec<Vec<f64>>) -> Result<Vec<f64>, OpError>,
    |fit: super::CoxFit, rows: Vec<Vec<f64>>| {
        if fit.beta.is_empty() || fit.beta.iter().any(|value| !value.is_finite())
            || rows.iter().flatten().any(|value| !value.is_finite()) {
            return Err(OpError::Domain("Cox prediction requires finite covariates and nonempty finite coefficients".into()));
        }
        crate::kernels::survival::cox_predict(&fit, &rows).map_err(OpError::Domain)
    }
);

binding!(
    survival_kaplan_meier, "num_survival_kaplan_meier", "numerical_stats",
    "Estimates the Kaplan–Meier curve and log-log confidence bounds. Event flags distinguish observed events from censoring. The estimate retains the input time dimension.",
    "crate::kernels::survival::kaplan_meier::kaplan_meier", DimensionPolicy::KaplanMeier,
    (time: Vec<f64>, event: Vec<bool>, confidence: f64) -> Result<super::KmEstimate, OpError>,
    |time: Vec<f64>, event: Vec<bool>, confidence: f64| {
        crate::kernels::survival::kaplan_meier(&time, &event, confidence).map_err(OpError::Domain)
    }
);

binding!(
    survival_median, "num_survival_median", "numerical_stats",
    "Returns the first measured time at or below one-half survival. Returns None if the fitted estimate never crosses one-half.",
    "crate::kernels::survival::kaplan_meier::median_survival", DimensionPolicy::MedianSurvival,
    (curve: super::KmEstimate) -> Option<f64>,
    |curve: super::KmEstimate| { crate::kernels::survival::median_survival(&curve) }
);

fn fit_cox_model(
    time: &[f64],
    event: &[bool],
    rows: &[Vec<f64>],
    ties: &str,
    max_iterations: usize,
    tolerance: f64,
) -> Result<CoxFit, OpError> {
    const MAX_WORK: usize = 50_000_000;
    let ties = match ties {
        | "breslow" => Ties::Breslow,
        | "efron" => Ties::Efron,
        | _ => return Err(OpError::Domain("Cox ties must be breslow or efron".into())),
    };
    let n = rows.len();
    let p = rows.first().map_or(1, Vec::len);
    let work = n.saturating_mul(
        p.saturating_mul(p)
            .saturating_mul(max_iterations)
            .saturating_add(n),
    );
    if work > MAX_WORK {
        return Err(OpError::Domain(
            "Cox fit exceeds the design or iteration work budget".into(),
        ));
    }
    crate::cancellation::checkpoint()?;
    crate::kernels::survival::cox_fit(rows, time, event, ties, max_iterations, tolerance)
        .map_err(OpError::Domain)
}
