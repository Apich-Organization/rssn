//! Iteratively reweighted least squares for generalised linear models.

use crate::kernels::dense::Mat;
use crate::kernels::dense::cholesky_solve;
use crate::kernels::dense::lstsq;
use crate::kernels::dense::pinv;
use crate::kernels::glm::family::Family;
use crate::kernels::glm::family::clamp_mean;
use crate::kernels::glm::family::has_unit_dispersion;
use crate::kernels::glm::family::starting_mean;
use crate::kernels::glm::family::unit_deviance;
use crate::kernels::glm::family::validate_response;
use crate::kernels::glm::family::variance;
use crate::kernels::glm::link::Link;
use crate::kernels::glm::link::link;
use crate::kernels::glm::link::link_derivative;
use crate::kernels::glm::link::link_inverse;
use crate::kernels::special::ln_gamma_numerical;

/// A fitted generalised linear model.
#[derive(Debug, Clone, PartialEq)]
pub struct GlmFit {
    /// Coefficient vector with the intercept first, then one entry per column
    /// of the design matrix.
    pub beta: Vec<f64>,
    /// Standard errors of the coefficients, from the inverse Fisher
    /// information scaled by the dispersion estimate.
    pub se: Vec<f64>,
    /// Residual deviance of the fitted model.
    pub deviance: f64,
    /// Deviance of the intercept-only model.
    pub null_deviance: f64,
    /// Number of reweighting iterations performed.
    pub iterations: usize,
    /// Whether the coefficient update fell below the requested tolerance.
    pub converged: bool,
    /// Akaike information criterion.
    pub aic: f64,
    /// Dispersion estimate: `1` for the Poisson and binomial families, the
    /// Pearson estimate otherwise.
    pub dispersion: f64,
    /// Fitted mean for each observation.
    pub fitted: Vec<f64>,
}

/// Number of step-halving trials allowed before an update is abandoned.
const MAX_HALVINGS: usize = 40;

/// Smallest working weight, keeping the weighted design full rank.
const MIN_WORKING_WEIGHT: f64 = 1e-12;

/// Fits a generalised linear model by iteratively reweighted least squares.
///
/// `x` holds one row of covariates per observation and must **not** contain an
/// intercept column; `fit` prepends one, so the returned coefficient vector has
/// `x[0].len() + 1` entries with the intercept first. `weights` are prior
/// weights (number of trials for a binomial proportion, exposure for a Poisson
/// rate); when `None` every observation weighs one.
///
/// Convergence is declared when `max |Δβ| / (1 + |β|) < tol`. Each step is
/// halved up to `MAX_HALVINGS` times while the deviance would increase.
///
/// # Errors
/// Returns a lowercase message when the design and response disagree in
/// length, the design is empty or non-finite, a weight is not strictly
/// positive, the response leaves the family's support, `max_iter` is zero,
/// `tol` is not positive, there are fewer observations than parameters, or a
/// linear system in the iteration is singular.
pub fn fit(
    x: &[Vec<f64>],
    y: &[f64],
    family: Family,
    kind: Link,
    weights: Option<&[f64]>,
    max_iter: usize,
    tol: f64,
) -> Result<GlmFit, String> {
    let n = y.len();
    if n == 0 {
        return Err("response is empty".to_string());
    }
    if x.len() != n {
        return Err("design and response lengths differ".to_string());
    }
    let k = x.first().map_or(0, Vec::len);
    if x.iter().any(|row| row.len() != k) {
        return Err("design rows have unequal lengths".to_string());
    }
    if x.iter().any(|row| row.iter().any(|v| !v.is_finite())) {
        return Err("design contains a non-finite value".to_string());
    }
    validate_response(family, y)?;
    if max_iter == 0 {
        return Err("max_iter must be positive".to_string());
    }
    if tol <= 0.0 || !tol.is_finite() {
        return Err("tolerance must be positive and finite".to_string());
    }
    let p = k + 1;
    if n < p {
        return Err("fewer observations than parameters".to_string());
    }

    let prior: Vec<f64> = match weights {
        | Some(ws) => {
            if ws.len() != n {
                return Err("weight and response lengths differ".to_string());
            }
            if ws.iter().any(|v| !v.is_finite() || *v <= 0.0) {
                return Err("weights must be finite and positive".to_string());
            }
            ws.to_vec()
        },
        | None => vec![1.0; n],
    };

    // Starting means from the response, then one weighted least-squares step.
    let mut mu0 = crate::kernels::reuse::f64_scratch(n);
    mu0.resize(n, 0.0);
    let mut eta0 = crate::kernels::reuse::f64_scratch(n);
    eta0.resize(n, 0.0);
    for i in 0..n {
        let m = clamp_mean(family, starting_mean(family, y[i]));
        mu0[i] = m;
        eta0[i] = link(kind, m);
    }
    let mut work_weights = crate::kernels::reuse::f64_scratch(n);
    work_weights.resize(n, 0.0);
    let mut work_response = crate::kernels::reuse::f64_scratch(n);
    work_response.resize(n, 0.0);
    let mut wls_rhs = crate::kernels::reuse::f64_scratch(n);
    wls_rhs.resize(n, 0.0);
    working_arrays_into(
        kind,
        family,
        y,
        &mu0,
        &eta0,
        &prior,
        &mut work_weights,
        &mut work_response,
    );
    let mut beta = weighted_least_squares(x, &work_response, &work_weights, &mut wls_rhs)?;

    // A pathological starting point can overflow the inverse link; fall back
    // to the constant predictor given by the mean starting linear predictor.
    let initial = match evaluate(family, kind, x, &beta, y, &prior) {
        | Ok(state) => state,
        | Err(_) => {
            let num = crate::kernels::simd::dot_f64(&prior, &eta0);
            let den = crate::kernels::simd::sum_f64(&prior);
            let mut constant = vec![0.0; p];
            if den > 0.0 {
                constant[0] = num / den;
            }
            beta = constant;
            evaluate(family, kind, x, &beta, y, &prior)?
        },
    };
    let (mut mu, mut eta, mut deviance) = initial;

    let mut converged = false;
    let mut iterations = 0;

    for _ in 0..max_iter {
        iterations += 1;
        working_arrays_into(
            kind,
            family,
            y,
            &mu,
            &eta,
            &prior,
            &mut work_weights,
            &mut work_response,
        );
        let proposal = weighted_least_squares(x, &work_response, &work_weights, &mut wls_rhs)?;

        let mut step = 1.0;
        let mut accepted: Option<(Vec<f64>, Vec<f64>, Vec<f64>, f64)> = None;
        for _ in 0..MAX_HALVINGS {
            let candidate: Vec<f64> = beta
                .iter()
                .zip(proposal.iter())
                .map(|(b, q)| b + step * (q - b))
                .collect();
            if let Ok((m, e, d)) = evaluate(family, kind, x, &candidate, y, &prior) {
                if d <= deviance + 1e-10 * (1.0 + deviance.abs()) {
                    accepted = Some((candidate, m, e, d));
                    break;
                }
            }
            step *= 0.5;
        }

        let Some((candidate, m, e, d)) = accepted else {
            break;
        };

        let delta = candidate
            .iter()
            .zip(beta.iter())
            .map(|(new, old)| (new - old).abs() / (1.0 + old.abs()))
            .fold(0.0_f64, f64::max);

        beta = candidate;
        mu = m;
        eta = e;
        deviance = d;

        if delta < tol {
            converged = true;
            break;
        }
    }

    // Intercept-only model, used for the null deviance.
    let num = crate::kernels::simd::dot_f64(&prior, y);
    let den = crate::kernels::simd::sum_f64(&prior);
    let grand_mean = clamp_mean(family, num / den);
    let mut null_deviance = 0.0;
    for i in 0..n {
        null_deviance += prior[i] * unit_deviance(family, y[i], grand_mean);
    }

    let dispersion = if has_unit_dispersion(family) {
        1.0
    } else {
        let df = n.saturating_sub(p);
        if df == 0 {
            1.0
        } else {
            let mut pearson = 0.0;
            for i in 0..n {
                let v = variance(family, mu[i]);
                if v > 0.0 {
                    let r = y[i] - mu[i];
                    pearson += prior[i] * r * r / v;
                }
            }
            pearson / df as f64
        }
    };

    let se = coefficient_standard_errors(x, &prior, &mu, family, kind, dispersion)?;

    let aic = if has_unit_dispersion(family) {
        deviance + 2.0 * p as f64
    } else {
        neg2_loglik(family, y, &prior, dispersion, deviance) + 2.0 * p as f64
    };

    Ok(GlmFit {
        beta,
        se,
        deviance,
        null_deviance,
        iterations,
        converged,
        aic,
        dispersion,
        fitted: mu,
    })
}

/// Working weights `1 / (V(mu) g'(mu)²)` and working response
/// `eta + (y - mu) g'(mu)` for the current mean, written into caller scratch.
fn working_arrays_into(
    kind: Link,
    family: Family,
    y: &[f64],
    mu: &[f64],
    eta: &[f64],
    prior: &[f64],
    work_weights: &mut [f64],
    work_response: &mut [f64],
) {
    let n = y.len();
    for i in 0..n {
        let m = clamp_mean(family, mu[i]);
        let g_prime = link_derivative(kind, m);
        let v = variance(family, m);
        let curvature = (v * g_prime * g_prime).abs().max(f64::MIN_POSITIVE);
        work_weights[i] = (prior[i] / curvature).max(MIN_WORKING_WEIGHT);
        work_response[i] = eta[i] + (y[i] - m) * g_prime;
    }
}

/// Solves `min Σ w_i (z_i - x_i'β)²` by scaling the design by `√w` and using a
/// QR least-squares solve.
fn weighted_least_squares(
    x: &[Vec<f64>],
    z: &[f64],
    work_weights: &[f64],
    rhs: &mut [f64],
) -> Result<Vec<f64>, String> {
    let n = x.len();
    let k = x.first().map_or(0, Vec::len);
    let p = k + 1;
    let mut design = Mat::zeros(n, p);
    for i in 0..n {
        let row = x.get(i).ok_or("missing design row")?;
        let wi = work_weights.get(i).copied().ok_or("missing weight")?;
        let scale = wi.max(0.0).sqrt();
        design.set(i, 0, scale);
        for j in 0..row.len() {
            let value = *row.get(j).ok_or("missing design entry")?;
            design.set(i, j + 1, scale * value);
        }
        rhs[i] = scale * z.get(i).copied().ok_or("missing working response")?;
    }
    lstsq(&design, rhs).map_err(|err| format!("weighted least squares failed: {err}"))
}

/// Evaluates the linear predictor, mean and deviance for a coefficient vector.
fn evaluate(
    family: Family,
    kind: Link,
    x: &[Vec<f64>],
    beta: &[f64],
    y: &[f64],
    prior: &[f64],
) -> Result<(Vec<f64>, Vec<f64>, f64), String> {
    let n = y.len();
    let intercept = *beta.first().ok_or("coefficient vector is empty")?;
    let slopes = &beta[1..];
    let mut mu = vec![0.0; n];
    let mut eta = vec![0.0; n];
    let mut deviance = 0.0;
    for i in 0..n {
        let row = x.get(i).ok_or("missing design row")?;
        if row.len() > slopes.len() {
            return Err("coefficient index out of range".to_string());
        }
        let linear = intercept + crate::kernels::simd::dot_f64(row, slopes);
        let raw = link_inverse(kind, linear);
        if !raw.is_finite() {
            return Err("linear predictor produced a non-finite mean".to_string());
        }
        let m = clamp_mean(family, raw);
        eta[i] = linear;
        mu[i] = m;
        deviance += prior[i] * unit_deviance(family, y[i], m);
    }
    if !deviance.is_finite() {
        return Err("deviance is not finite".to_string());
    }
    Ok((mu, eta, deviance))
}

/// Diagonal of the inverse Fisher information, scaled by the dispersion.
fn coefficient_standard_errors(
    x: &[Vec<f64>],
    prior: &[f64],
    mu: &[f64],
    family: Family,
    kind: Link,
    dispersion: f64,
) -> Result<Vec<f64>, String> {
    let n = x.len();
    let k = x.first().map_or(0, Vec::len);
    let p = k + 1;
    let mut gram = Mat::zeros(p, p);
    let mut row = vec![0.0; p];
    for i in 0..n {
        let m = clamp_mean(family, mu[i]);
        let g_prime = link_derivative(kind, m);
        let v = variance(family, m);
        let curvature = (v * g_prime * g_prime).abs().max(f64::MIN_POSITIVE);
        let wi = prior[i] / curvature;
        row[0] = 1.0;
        if let Some(covariates) = x.get(i) {
            if let Some(slot) = row.get_mut(1..) {
                if slot.len() == covariates.len() {
                    slot.copy_from_slice(covariates);
                }
            }
        }
        for a in 0..p {
            for b in 0..p {
                gram.data[a * p + b] += wi * row[a] * row[b];
            }
        }
    }

    let inverse_diagonal = inverse_diagonal(&gram)?;
    let mut se = vec![0.0; p];
    for j in 0..p {
        let entry = dispersion * inverse_diagonal[j];
        se[j] = if entry > 0.0 {
            entry.sqrt()
        } else {
            0.0
        };
    }
    Ok(se)
}

/// Diagonal of the inverse of a symmetric positive semi-definite matrix,
/// falling back to the pseudo-inverse when the Cholesky factor does not exist.
fn inverse_diagonal(a: &Mat) -> Result<Vec<f64>, String> {
    let p = a.rows;
    let mut out = vec![0.0; p];
    let mut unit = vec![0.0; p];
    let mut cholesky_ok = true;
    for j in 0..p {
        unit.fill(0.0);
        unit[j] = 1.0;
        match cholesky_solve(a, &unit) {
            | Ok(column) => out[j] = column[j],
            | Err(_) => {
                cholesky_ok = false;
                break;
            },
        }
    }
    if cholesky_ok {
        return Ok(out);
    }
    let inverse = pinv(a).map_err(|_| "fisher information is singular".to_string())?;
    for j in 0..p {
        out[j] = inverse.at(j, j);
    }
    Ok(out)
}

/// `-2` times the maximised log-likelihood, written as the deviance plus the
/// saturated-model constant.
fn neg2_loglik(
    family: Family,
    y: &[f64],
    prior: &[f64],
    dispersion: f64,
    deviance: f64,
) -> f64 {
    if dispersion <= 0.0 || !dispersion.is_finite() {
        return deviance;
    }
    let mut saturated = 0.0;
    match family {
        | Family::Gaussian => {
            for i in 0..y.len() {
                let wi = prior[i].max(f64::MIN_POSITIVE);
                saturated += wi * ((2.0 * std::f64::consts::PI * dispersion / wi).ln() + 1.0);
            }
        },
        | Family::Gamma => {
            let shape = 1.0 / dispersion;
            let ln_gamma = ln_gamma_numerical(shape);
            let base = ln_gamma - shape * shape.ln() + shape;
            for i in 0..y.len() {
                let log_y = if y[i] > 0.0 { y[i].ln() } else { 0.0 };
                saturated += 2.0 * prior[i] * (base - (shape - 1.0) * log_y);
            }
        },
        | Family::InverseGaussian => {
            for i in 0..y.len() {
                let cube = y[i] * y[i] * y[i];
                let term = (2.0 * std::f64::consts::PI * dispersion * cube).ln();
                saturated += prior[i] * (term + 1.0 / dispersion);
            }
        },
        | Family::Poisson | Family::Binomial => {},
    }
    deviance / dispersion + saturated
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernels::dense::Mat;

    const TOL: f64 = 1e-6;

    #[test]
    fn gaussian_identity_matches_ordinary_least_squares() {
        let x = vec![vec![1.0], vec![2.0], vec![3.0], vec![4.0], vec![5.0]];
        let y = vec![2.1, 3.9, 6.2, 7.8, 10.1];

        let result =
            fit(&x, &y, Family::Gaussian, Link::Identity, None, 100, 1e-10).expect("gaussian fit");

        // Reference: least squares on the design with an explicit intercept.
        let mut design = Mat::zeros(y.len(), 2);
        for i in 0..y.len() {
            design.set(i, 0, 1.0);
            design.set(i, 1, x[i][0]);
        }
        let reference = lstsq(&design, &y).expect("reference least squares");

        assert!(
            (result.beta[0] - reference[0]).abs() < TOL,
            "{:?}",
            result.beta
        );
        assert!(
            (result.beta[1] - reference[1]).abs() < TOL,
            "{:?}",
            result.beta
        );
        assert!(result.converged);
        // Residuals are orthogonal to the design at the least-squares solution.
        let mut normal0 = 0.0;
        let mut normal1 = 0.0;
        for i in 0..y.len() {
            let r = y[i] - result.fitted[i];
            normal0 += r;
            normal1 += r * x[i][0];
        }
        assert!(normal0.abs() < 1e-8);
        assert!(normal1.abs() < 1e-8);
        assert!((result.dispersion - 0.0).abs() > 0.0);
    }

    #[test]
    fn poisson_log_satisfies_the_score_equations() {
        let x = vec![vec![0.0], vec![1.0], vec![2.0], vec![3.0], vec![4.0]];
        let y = vec![2.0, 3.0, 5.0, 8.0, 13.0];

        let result =
            fit(&x, &y, Family::Poisson, Link::Log, None, 200, 1e-12).expect("poisson fit");

        assert!(result.converged);
        assert_eq!(result.dispersion, 1.0);
        // Score equations: Σ (y - mu) = 0 and Σ x (y - mu) = 0.
        let mut score0 = 0.0;
        let mut score1 = 0.0;
        for i in 0..y.len() {
            let residual = y[i] - result.fitted[i];
            score0 += residual;
            score1 += x[i][0] * residual;
        }
        assert!(score0.abs() < 1e-8, "score0 = {score0}");
        assert!(score1.abs() < 1e-8, "score1 = {score1}");
        // The fitted means are the exponential of the linear predictor.
        for i in 0..y.len() {
            let eta = result.beta[0] + result.beta[1] * x[i][0];
            assert!((result.fitted[i] - eta.exp()).abs() < 1e-10);
        }
    }

    #[test]
    fn logistic_two_group_table_recovers_group_logits() {
        // Two groups of three Bernoulli trials with one success in the first
        // and two in the second; the fitted group probabilities are the sample
        // proportions, so the coefficients are the group logits.
        let x = vec![
            vec![0.0],
            vec![0.0],
            vec![0.0],
            vec![1.0],
            vec![1.0],
            vec![1.0],
        ];
        let y = vec![0.0, 0.0, 1.0, 0.0, 1.0, 1.0];

        let result =
            fit(&x, &y, Family::Binomial, Link::Logit, None, 200, 1e-12).expect("logistic fit");

        let p0: f64 = 1.0 / 3.0;
        let p1: f64 = 2.0 / 3.0;
        let expected0 = (p0 / (1.0 - p0)).ln();
        let expected1 = (p1 / (1.0 - p1)).ln() - expected0;

        assert!(
            (result.beta[0] - expected0).abs() < TOL,
            "{:?}",
            result.beta
        );
        assert!(
            (result.beta[1] - expected1).abs() < TOL,
            "{:?}",
            result.beta
        );
        assert!((result.fitted[0] - p0).abs() < TOL);
        assert!((result.fitted[5] - p1).abs() < TOL);

        // Residual deviance of the two-group table: two observations with
        // y = 0 and mu = 1/3, one with y = 1 and mu = 1/3, one with y = 0 and
        // mu = 2/3, and two with y = 1 and mu = 2/3.
        let expected_deviance = 8.0 * 1.5_f64.ln() + 4.0 * 3.0_f64.ln();
        assert!(
            (result.deviance - expected_deviance).abs() < 1e-8,
            "{} vs {}",
            result.deviance,
            expected_deviance
        );
    }

    #[test]
    fn binomial_deviance_matches_the_hand_computation() {
        // Intercept-only model: the design carries no covariates at all.
        let x = vec![Vec::new(), Vec::new()];
        let y = vec![0.2, 0.6];
        let w = vec![10.0, 10.0];

        let result =
            fit(&x, &y, Family::Binomial, Link::Logit, Some(&w), 200, 1e-12).expect("binomial fit");

        // Intercept-only logit gives mu = weighted mean = 0.4 for both groups.
        assert!((result.fitted[0] - 0.4).abs() < TOL);
        assert!((result.fitted[1] - 0.4).abs() < TOL);

        let mu = 0.4_f64;
        let expected = 20.0 * (0.2 * (0.2_f64 / mu).ln() + 0.8 * (0.8_f64 / 0.6).ln())
            + 20.0 * (0.6 * (0.6_f64 / mu).ln() + 0.4 * (0.4_f64 / 0.6).ln());
        assert!(
            (result.deviance - expected).abs() < 1e-9,
            "{} vs {}",
            result.deviance,
            expected
        );
        // Null deviance of an intercept-only model equals the residual one.
        assert!((result.null_deviance - result.deviance).abs() < 1e-9);
    }

    #[test]
    fn gamma_inverse_has_positive_fitted_values_and_finite_errors() {
        let x = vec![vec![0.0], vec![1.0], vec![2.0], vec![3.0], vec![4.0]];
        let y = vec![1.2, 2.1, 3.0, 4.4, 5.1];

        let result =
            fit(&x, &y, Family::Gamma, Link::Inverse, None, 200, 1e-12).expect("gamma fit");

        assert!(result.fitted.iter().all(|m| *m > 0.0 && m.is_finite()));
        assert!(result.se.iter().all(|s| s.is_finite() && *s >= 0.0));
        assert!(result.dispersion > 0.0);
        assert!(result.aic.is_finite());
    }

    #[test]
    fn probit_and_loglog_links_converge_on_binomial_data() {
        let x = vec![
            vec![-2.0],
            vec![-1.0],
            vec![-0.5],
            vec![0.5],
            vec![1.0],
            vec![2.0],
        ];
        let y = vec![0.05, 0.2, 0.35, 0.65, 0.8, 0.95];

        for kind in [Link::Probit, Link::LogLog, Link::ComplementaryLogLog] {
            let result =
                fit(&x, &y, Family::Binomial, kind, None, 300, 1e-12).expect("binomial fit");
            assert!(
                result.fitted.iter().all(|m| *m > 0.0 && *m < 1.0),
                "{kind:?}"
            );
            assert!(result.converged, "{kind:?}");
        }
    }

    #[test]
    fn inverse_gaussian_fit_is_finite() {
        let x = vec![vec![0.0], vec![1.0], vec![2.0], vec![3.0]];
        let y = vec![2.0, 2.5, 3.1, 3.9];
        let result = fit(&x, &y, Family::InverseGaussian, Link::Log, None, 200, 1e-10)
            .expect("inverse gaussian fit");
        assert!(result.fitted.iter().all(|m| m.is_finite() && *m > 0.0));
        assert!(result.beta.iter().all(|b| b.is_finite()));
    }

    #[test]
    fn standard_errors_match_the_gaussian_ols_formula() {
        let x = vec![vec![0.0], vec![1.0], vec![2.0], vec![3.0], vec![4.0]];
        let y = vec![1.0, 2.2, 2.8, 4.1, 5.2];
        let result = fit(&x, &y, Family::Gaussian, Link::Identity, None, 100, 1e-12).expect("fit");

        // Classic OLS standard errors for a simple linear regression.
        let n = y.len() as f64;
        let xbar = 2.0;
        let sxx: f64 = (0..5).map(|i| (x[i][0] - xbar) * (x[i][0] - xbar)).sum();
        let sigma2 = result.dispersion;
        let se_intercept = (sigma2 * (1.0 / n + xbar * xbar / sxx)).sqrt();
        let se_slope = (sigma2 / sxx).sqrt();
        assert!((result.se[0] - se_intercept).abs() < 1e-9);
        assert!((result.se[1] - se_slope).abs() < 1e-9);
    }

    #[test]
    fn invalid_inputs_are_rejected() {
        let empty: Vec<Vec<f64>> = Vec::new();
        assert!(
            fit(
                &empty,
                &[],
                Family::Gaussian,
                Link::Identity,
                None,
                10,
                1e-8
            )
            .is_err()
        );

        let x = vec![vec![1.0], vec![2.0]];
        assert!(fit(&x, &[1.0], Family::Gaussian, Link::Identity, None, 10, 1e-8).is_err());
        assert!(
            fit(&x, &[-1.0, 2.0], Family::Poisson, Link::Log, None, 10, 1e-8).is_err(),
            "negative poisson response must fail"
        );
        assert!(
            fit(&x, &[0.0, 2.0], Family::Gamma, Link::Log, None, 10, 1e-8).is_err(),
            "zero gamma response must fail"
        );
        assert!(
            fit(
                &x,
                &[1.5, 2.0],
                Family::Binomial,
                Link::Logit,
                None,
                10,
                1e-8
            )
            .is_err(),
            "binomial response above one must fail"
        );
        assert!(
            fit(
                &x,
                &[1.0, 2.0],
                Family::Gaussian,
                Link::Identity,
                Some(&[1.0]),
                10,
                1e-8
            )
            .is_err(),
            "weight length mismatch must fail"
        );
        assert!(
            fit(
                &x,
                &[1.0, 2.0],
                Family::Gaussian,
                Link::Identity,
                Some(&[1.0, 0.0]),
                10,
                1e-8
            )
            .is_err(),
            "zero weight must fail"
        );
        assert!(
            fit(
                &x,
                &[1.0, 2.0],
                Family::Gaussian,
                Link::Identity,
                None,
                0,
                1e-8
            )
            .is_err()
        );
        assert!(
            fit(
                &x,
                &[1.0, 2.0],
                Family::Gaussian,
                Link::Identity,
                None,
                10,
                0.0
            )
            .is_err()
        );

        let ragged = vec![vec![1.0, 2.0], vec![3.0]];
        assert!(
            fit(
                &ragged,
                &[1.0, 2.0],
                Family::Gaussian,
                Link::Identity,
                None,
                10,
                1e-8
            )
            .is_err()
        );
    }
}
