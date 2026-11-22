//! Generalised linear models.
//!
//! | operator | value |
//! |---|---|
//! | `glm_fit(x, y, family, link, weights, tol)` | `list(beta, se, deviance, null_deviance, aic, iterations, fitted)`: the maximum-likelihood fit of a generalised linear model |
//! | `glm_predict(x, beta, link)` | the fitted mean of each row of `x` for the coefficients `beta` (intercept first) |
//!
//! `x` is a list of rows **without** an intercept column; `glm_fit` prepends
//! one. `family` is `gaussian`, `poisson`, `binomial`, `gamma` or
//! `inverse_gaussian`; `link` is `identity`, `log`, `logit`, `probit`,
//! `inverse`, `loglog` or `cloglog`. `weights` is a list of the same length as
//! `y`, or the empty list for equal weights. The fit runs iteratively
//! reweighted least squares with step halving and stops when the largest
//! coefficient update falls below `tol` (or after 100 iterations).
//!
//! A response the family cannot hold (a negative Poisson count, a binomial
//! proportion outside `[0, 1]`, a non-positive gamma response) is a
//! [`DiagKind::Domain`] finding, and the call stays a request rather than
//! returning a fit of a model the data do not support.

use super::bridge::Data;
use super::bridge::Entry;
use super::bridge::bad;
use super::bridge::name;
use super::bridge::require;
use super::bridge::rows;
use super::bridge::scalar;
use super::bridge::vector;
use crate::kernel::diag::DiagKind;
use crate::kernels::glm::Family;
use crate::kernels::glm::Link;

/// Most observations a fit accepts.
const MAX_ROWS: usize = 100_000;

/// Every operator of this family.
pub(super) const ENTRIES: &[Entry] = &[
    Entry {
        name: "glm_fit",
        params: "ttnnts",
        eval: glm_fit,
        doc: "`glm_fit(x, y, family, link, weights, tol)`: `list(beta, se, deviance, null_deviance, aic, iterations, fitted)` of a generalised linear model fitted by iteratively reweighted least squares; `x` is a list of rows without an intercept, `weights` is empty for equal weights.",
    },
    Entry {
        name: "glm_predict",
        params: "ttn",
        eval: glm_predict,
        doc: "`glm_predict(x, beta, link)`: the fitted mean of each row of `x` given the coefficients `beta` (intercept first) and the `link`.",
    },
];

/// The family named by `text`.
fn family_named(text: &str) -> Option<Family> {
    match text {
        | "gaussian" | "normal" => Some(Family::Gaussian),
        | "poisson" => Some(Family::Poisson),
        | "binomial" => Some(Family::Binomial),
        | "gamma" => Some(Family::Gamma),
        | "inverse_gaussian" | "inverse-gaussian" => Some(Family::InverseGaussian),
        | _ => None,
    }
}

/// The link named by `text`.
fn link_named(text: &str) -> Option<Link> {
    match text {
        | "identity" => Some(Link::Identity),
        | "log" => Some(Link::Log),
        | "logit" => Some(Link::Logit),
        | "probit" => Some(Link::Probit),
        | "inverse" => Some(Link::Inverse),
        | "loglog" | "log_log" => Some(Link::LogLog),
        | "cloglog" | "complementary_loglog" => Some(Link::ComplementaryLogLog),
        | _ => None,
    }
}

/// The design matrix at `at`.
fn design(
    args: &[Data],
    at: usize,
) -> Option<Vec<Vec<f64>>> {
    let Some(matrix) = args.get(at).and_then(rows) else {
        return bad(
            at,
            DiagKind::Shape,
            format_args!("expected a list of rows of numbers"),
        );
    };
    require(
        !matrix.is_empty() && matrix.len() <= MAX_ROWS,
        Some(at),
        DiagKind::Shape,
        format_args!("expected 1 to {MAX_ROWS} rows, got {}", matrix.len()),
    )?;
    Some(matrix)
}

fn glm_fit(args: &[Data]) -> Option<Data> {
    let x = design(args, 0)?;
    let y = vector(args, 1)?;
    let family_text = name(args, 2)?;
    let link_text = name(args, 3)?;
    let Some(family) = family_named(family_text) else {
        return bad(
            2,
            DiagKind::Domain,
            format_args!(
                "`{family_text}` is not a family (gaussian, poisson, binomial, gamma, \
                 inverse_gaussian)"
            ),
        );
    };
    let Some(link) = link_named(link_text) else {
        return bad(
            3,
            DiagKind::Domain,
            format_args!(
                "`{link_text}` is not a link (identity, log, logit, probit, inverse, loglog, \
                 cloglog)"
            ),
        );
    };
    let weights = match args.get(4).and_then(super::conic::numbers_of) {
        | Some(values) if values.is_empty() => None,
        | Some(values) => Some(values),
        | None => {
            return bad(
                4,
                DiagKind::Shape,
                format_args!("expected a list of weights, or the empty list for equal weights"),
            );
        },
    };
    let tolerance = scalar(args, 5)?;
    require(
        y.len() == x.len(),
        Some(1),
        DiagKind::Shape,
        format_args!(
            "the response has {} entries but the design has {} rows",
            y.len(),
            x.len()
        ),
    )?;
    let Some(width) = x.first().map(Vec::len) else {
        return bad(0, DiagKind::Shape, format_args!("the design is empty"));
    };
    require(
        x.iter().all(|row| row.len() == width),
        Some(0),
        DiagKind::Shape,
        format_args!("the design rows have different lengths"),
    )?;
    require(
        tolerance > 0.0 && tolerance.is_finite(),
        Some(5),
        DiagKind::Domain,
        format_args!("the tolerance must be positive and finite, got {tolerance}"),
    )?;
    if let Some(w) = weights.as_deref()
        && w.len() != y.len()
    {
        return bad(
            4,
            DiagKind::Shape,
            format_args!(
                "the weights have {} entries but the response has {}",
                w.len(),
                y.len()
            ),
        );
    }
    let fit =
        match crate::kernels::glm::fit(&x, y, family, link, weights.as_deref(), 100, tolerance) {
            | Ok(fit) => fit,
            | Err(message) => return bad(0, DiagKind::Domain, format_args!("{message}")),
        };
    #[allow(clippy::cast_precision_loss)]
    let iterations = fit.iterations as f64;
    Some(Data::Tree(vec![
        Data::vector(fit.beta),
        Data::vector(fit.se),
        Data::Scalar(fit.deviance),
        Data::Scalar(fit.null_deviance),
        Data::Scalar(fit.aic),
        Data::Scalar(iterations),
        Data::vector(fit.fitted),
    ]))
}

fn glm_predict(args: &[Data]) -> Option<Data> {
    let x = design(args, 0)?;
    let beta = vector(args, 1)?;
    let link_text = name(args, 2)?;
    let Some(link) = link_named(link_text) else {
        return bad(
            2,
            DiagKind::Domain,
            format_args!(
                "`{link_text}` is not a link (identity, log, logit, probit, inverse, loglog, \
                 cloglog)"
            ),
        );
    };
    let Some(width) = x.first().map(Vec::len) else {
        return bad(0, DiagKind::Shape, format_args!("the design is empty"));
    };
    require(
        beta.len() == width.saturating_add(1),
        Some(1),
        DiagKind::Shape,
        format_args!(
            "the coefficients need {} entries (intercept and one per column), got {}",
            width.saturating_add(1),
            beta.len()
        ),
    )?;
    let mut means = Vec::with_capacity(x.len());
    for row in &x {
        let mut eta = beta.first().copied().unwrap_or(0.0);
        for (value, coefficient) in row.iter().zip(beta.iter().skip(1)) {
            eta += value * coefficient;
        }
        means.push(crate::kernels::glm::link_inverse(link, eta));
    }
    super::bridge::vec_out(means)
}
