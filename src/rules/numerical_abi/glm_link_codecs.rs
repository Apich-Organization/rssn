// The same link and family conventions used by the GLM fit operators.
use crate::kernels::glm;

fn glm_link_named(name: &str) -> Result<glm::Link, OpError> {
    glm::Link::named(name).ok_or_else(|| OpError::Domain(format!("unknown GLM link `{name}`")))
}

fn glm_mean_in_link_domain(
    kind: glm::Link,
    mean: f64,
) -> Result<(), OpError> {
    let supported = match kind {
        | glm::Link::Identity => true,
        | glm::Link::Log => mean > 0.0,
        | glm::Link::Inverse => mean != 0.0,
        | glm::Link::Logit
        | glm::Link::Probit
        | glm::Link::LogLog
        | glm::Link::ComplementaryLogLog => mean > 0.0 && mean < 1.0,
    };
    if !mean.is_finite() || !supported {
        return Err(OpError::Domain(format!(
            "{mean} lies outside the finite mean domain of {kind:?}"
        )));
    }
    Ok(())
}

binding!(
    glm_link_value,
    "num_glm_link",
    "numerical_stats",
    "Maps a dimensionless mean to the linear predictor. Choose identity, log, logit, probit, inverse, loglog, or cloglog. Means at a boundary with an infinite link are rejected.",
    "crate::kernels::glm::link::link",
    DimensionPolicy::Dimensionless,
    (name: String, mean: f64) -> Result<f64, OpError>,
    |name: String, mean: f64| {
        let kind = super::glm_link_named(&name)?;
        super::glm_mean_in_link_domain(kind, mean)?;
        Ok(super::glm::link::link(kind, mean))
    }
);

binding!(
    glm_inverse_link_value,
    "num_glm_inverse_link",
    "numerical_stats",
    "Maps a finite dimensionless linear predictor to its response mean. Uses the same named links as num_glm_link and preserves representable lower-tail probabilities.",
    "crate::kernels::glm::link::link_inverse",
    DimensionPolicy::Dimensionless,
    (name: String, predictor: f64) -> Result<f64, OpError>,
    |name: String, predictor: f64| {
        let kind = super::glm_link_named(&name)?;
        if !predictor.is_finite() {
            return Err(OpError::Domain("GLM predictor must be finite".into()));
        }
        let mean = super::glm::link::link_inverse(kind, predictor);
        if !mean.is_finite() {
            return Err(OpError::Domain("inverse GLM link has no finite mean at this predictor".into()));
        }
        Ok(mean)
    }
);

binding!(
    glm_link_derivative,
    "num_glm_link_derivative",
    "numerical_stats",
    "Returns the dimensionless derivative of a named GLM link at an admissible mean. The reciprocal link has a negative derivative.",
    "crate::kernels::glm::link::link_derivative",
    DimensionPolicy::Dimensionless,
    (name: String, mean: f64) -> Result<f64, OpError>,
    |name: String, mean: f64| {
        let kind = super::glm_link_named(&name)?;
        super::glm_mean_in_link_domain(kind, mean)?;
        let derivative = super::glm::link::link_derivative(kind, mean);
        if !derivative.is_finite() {
            return Err(OpError::Domain("GLM link derivative exceeds binary64 range".into()));
        }
        Ok(derivative)
    }
);

binding!(
    glm_variance,
    "num_glm_variance",
    "numerical_stats",
    "Returns the dimensionless variance function at a mean in the named family's support. Choose gaussian, poisson, binomial, gamma, or inverse_gaussian.",
    "crate::kernels::glm::family::variance",
    DimensionPolicy::Dimensionless,
    (name: String, mean: f64) -> Result<f64, OpError>,
    |name: String, mean: f64| {
        let family = super::glm::Family::named(&name)
            .ok_or_else(|| OpError::Domain(format!("unknown GLM family `{name}`")))?;
        super::glm::family::validate_response(family, std::slice::from_ref(&mean))
            .map_err(OpError::Domain)?;
        let variance = super::glm::family::variance(family, mean);
        if !variance.is_finite() {
            return Err(OpError::Domain("GLM variance exceeds binary64 range".into()));
        }
        Ok(variance)
    }
);

use crate::kernels::glm::GlmFit;

record!(GlmFit {
    beta: Vec<f64>,
    se: Vec<f64>,
    deviance: f64,
    null_deviance: f64,
    iterations: usize,
    converged: bool,
    aic: f64,
    dispersion: f64,
    fitted: Vec<f64>
});

struct GlmModel {
    estimation: GlmFit,
    link: String,
}
record!(GlmModel {
    estimation: GlmFit,
    link: String
});

impl Arg for GlmFit {
    fn ty() -> TypeExpr {
        <Self as Output>::ty()
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let fields = [
            "beta",
            "se",
            "deviance",
            "null_deviance",
            "iterations",
            "converged",
            "aic",
            "dispersion",
            "fitted",
        ];
        let field = |name| model_field(value, name, &fields);
        Ok(Self {
            beta: Vec::<f64>::decode(field("beta")?, index)?,
            se: Vec::<f64>::decode(field("se")?, index)?,
            deviance: f64::decode(field("deviance")?, index)?,
            null_deviance: f64::decode(field("null_deviance")?, index)?,
            iterations: usize::decode(field("iterations")?, index)?,
            converged: bool::decode(field("converged")?, index)?,
            aic: f64::decode(field("aic")?, index)?,
            dispersion: f64::decode(field("dispersion")?, index)?,
            fitted: Vec::<f64>::decode(field("fitted")?, index)?,
        })
    }
}

impl Arg for GlmModel {
    fn ty() -> TypeExpr {
        <Self as Output>::ty()
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let fields = ["estimation", "link"];
        Ok(Self {
            estimation: GlmFit::decode(model_field(value, "estimation", &fields)?, index)?,
            link: String::decode(model_field(value, "link", &fields)?, index)?,
        })
    }
}

const MAX_GLM_ROWS: usize = 100_000;
const MAX_GLM_FIT_WORK: usize = 50_000_000;

fn glm_design_size(rows: &[Vec<f64>]) -> usize {
    let columns = rows.first().map_or(1, |row| row.len().saturating_add(1));
    rows.len().saturating_mul(columns)
}

binding!(
    glm_fitted_model,
    "num_glm_fit",
    "numerical_stats",
    "Fits a generalized linear model. Design rows omit the intercept; optional weights must align with observations. Returns a portable named model with coefficients, convergence status and the inverse link used for prediction.",
    "crate::kernels::glm::fit::fit",
    DimensionPolicy::Dimensionless,
    (rows: Vec<Vec<f64>>, observed: Vec<f64>, family: String, link: String, weights: Option<Vec<f64>>, max_iterations: usize, tolerance: f64)
        -> Result<super::GlmModel, OpError>,
    |rows: Vec<Vec<f64>>, observed: Vec<f64>, family: String, link: String, weights: Option<Vec<f64>>, max_iterations: usize, tolerance: f64| {
        let family = super::glm::Family::named(&family)
            .ok_or_else(|| OpError::Domain(format!("unknown GLM family `{family}`")))?;
        let kind = super::glm_link_named(&link)?;
        let columns = rows.first().map_or(1, |row| row.len().saturating_add(1));
        let work = super::glm_design_size(&rows).saturating_mul(columns).saturating_mul(max_iterations);
        if rows.len() > super::MAX_GLM_ROWS || work > super::MAX_GLM_FIT_WORK {
            return Err(OpError::Domain("GLM fit exceeds the design or iteration work budget".into()));
        }
        let estimation = super::glm::fit(&rows, &observed, family, kind, weights.as_deref(), max_iterations, tolerance)
            .map_err(OpError::Domain)?;
        Ok(super::GlmModel { estimation, link })
    }
);
binding!(
    glm_model_predict,
    "num_glm_predict",
    "numerical_stats",
    "Predicts response means from a fitted model with the model's stored link. Design rows omit the intercept and must match the fitted slope count.",
    "crate::kernels::glm::predict::predict",
    DimensionPolicy::Dimensionless,
    (model: super::GlmModel, rows: Vec<Vec<f64>>) -> Result<Vec<f64>, OpError>,
    |model: super::GlmModel, rows: Vec<Vec<f64>>| {
        if rows.len() > super::MAX_GLM_ROWS || super::glm_design_size(&rows) > super::MAX_GLM_FIT_WORK {
            return Err(OpError::Domain("GLM prediction exceeds the design work budget".into()));
        }
        let kind = super::glm_link_named(&model.link)?;
        super::glm::predict(&model.estimation, &rows, kind).map_err(OpError::Domain)
    }
);
