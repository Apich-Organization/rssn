//! Cox proportional-hazards regression.

use crate::kernels::dense::Mat;
use crate::kernels::dense::lu_factor;
use crate::kernels::dense::solve;
use crate::kernels::distributions::Distribution;
use crate::kernels::reuse::Scratch;
use crate::kernels::reuse::f64_scratch;
use crate::kernels::survival::common::group_end;
use crate::kernels::survival::common::time_order;
use std::ops::Range;

/// Handling of tied event times in the Cox partial likelihood.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Ties {
    /// Breslow approximation: all tied events share one risk-set term.
    Breslow,
    /// Efron approximation: tied events leave the risk set one at a time.
    Efron,
}

/// Maximum-likelihood fit of a Cox proportional-hazards model.
#[derive(Clone, Debug, PartialEq)]
pub struct CoxFit {
    /// Estimated log hazard ratios, one per covariate.
    pub beta: Vec<f64>,
    /// Standard errors from the inverse observed information.
    pub se: Vec<f64>,
    /// Wald statistics `beta / se`.
    pub z: Vec<f64>,
    /// Two-sided p-values of the Wald statistics.
    pub p: Vec<f64>,
    /// Log partial likelihood at the estimate.
    pub loglik: f64,
    /// Newton-Raphson iterations performed.
    pub iterations: usize,
    /// Whether the last accepted step fell below the tolerance.
    pub converged: bool,
    /// Harrell's concordance index over comparable pairs.
    pub concordance: f64,
}

/// One distinct observation time in descending order.
struct TimeGroup {
    rows: Range<usize>,
    events: usize,
}

struct RiskOrder {
    indices: Vec<usize>,
    groups: Vec<TimeGroup>,
}

impl RiskOrder {
    fn new(
        time: &[f64],
        event: &[bool],
    ) -> Self {
        let indices = time_order(time, true);
        let mut groups = Vec::new();
        let mut start = 0;
        while start < indices.len() {
            let end = group_end(time, &indices, start);
            let events = indices[start..end].iter().filter(|&&i| event[i]).count();
            groups.push(TimeGroup {
                rows: start..end,
                events,
            });
            start = end;
        }
        Self { indices, groups }
    }
}

/// Risk and observed-event weighted moments. The lower triangle suffices
/// because the information matrix is symmetric.
struct Moments {
    weight: f64,
    first: Scratch<'static, f64>,
    second: Scratch<'static, f64>,
}

impl Moments {
    fn new(dimensions: usize) -> Self {
        let cells = dimensions * dimensions;
        let mut first = f64_scratch(dimensions);
        let mut second = f64_scratch(cells);
        first.resize(dimensions, 0.0);
        second.resize(cells, 0.0);
        Self {
            weight: 0.0,
            first,
            second,
        }
    }

    fn reset(&mut self) {
        self.weight = 0.0;
        self.first.fill(0.0);
        self.second.fill(0.0);
    }

    fn rescale(
        &mut self,
        factor: f64,
    ) {
        self.weight *= factor;
        for value in self.first.iter_mut().chain(self.second.iter_mut()) {
            *value *= factor;
        }
    }

    fn add(
        &mut self,
        row: &[f64],
        weight: f64,
    ) {
        self.weight += weight;
        for (a, ((first, second), &value)) in self
            .first
            .iter_mut()
            .zip(self.second.chunks_exact_mut(row.len()))
            .zip(row)
            .enumerate()
        {
            let scaled = weight * value;
            *first += scaled;
            for (cell, &value) in second.iter_mut().zip(row).take(a + 1) {
                *cell += scaled * value;
            }
        }
    }
}

/// Dot product of two equally long slices.
fn dot(
    a: &[f64],
    b: &[f64],
) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// Checks the design matrix and survival data and returns the covariate count.
fn validate(
    x: &[Vec<f64>],
    time: &[f64],
    event: &[bool],
) -> Result<usize, String> {
    let n = time.len();
    if n == 0 {
        return Err("no observations".to_string());
    }
    if event.len() != n {
        return Err("time and event must have the same length".to_string());
    }
    if x.len() != n {
        return Err("the design matrix must have one row per observation".to_string());
    }
    if time.iter().any(|t| !t.is_finite()) {
        return Err("observation times must be finite".to_string());
    }
    let p = x.first().map_or(0, Vec::len);
    if p == 0 {
        return Err("the design matrix needs at least one covariate".to_string());
    }
    if x.iter().any(|row| row.len() != p) {
        return Err("the design matrix rows must have the same length".to_string());
    }
    if x.iter().flatten().any(|v| !v.is_finite()) {
        return Err("covariates must be finite".to_string());
    }
    if !event.iter().any(|&e| e) {
        return Err("all observations are censored".to_string());
    }
    for j in 0..p {
        let mean = x.iter().map(|row| row[j]).sum::<f64>() / n as f64;
        let variance = x
            .iter()
            .map(|row| (row[j] - mean) * (row[j] - mean))
            .sum::<f64>()
            / n as f64;
        if !variance.is_finite() || variance <= 0.0 {
            return Err("a covariate column has zero variance".to_string());
        }
    }
    Ok(p)
}

/// Log partial likelihood, gradient and observed information at `beta`.
///
/// The descending time sweep extends a scaled risk set once per observation.
/// At each event the numerator uses eta minus the current shift, as does the
/// denominator. This cancels the exponent scaling without adding it twice.
fn evaluate(
    beta: &[f64],
    x: &[Vec<f64>],
    event: &[bool],
    order: &RiskOrder,
    ties: Ties,
) -> (f64, Vec<f64>, Mat) {
    let p = beta.len();
    let mut eta = f64_scratch(x.len());
    eta.extend(x.iter().map(|row| dot(row, beta)));
    let mut loglik = 0.0;
    let mut grad = vec![0.0; p];
    let mut info = Mat::zeros(p, p);
    if eta.iter().any(|value| !value.is_finite()) {
        return (f64::NEG_INFINITY, grad, info);
    }

    let mut shift = f64::NEG_INFINITY;
    let mut risk = Moments::new(p);
    let mut removed = Moments::new(p);
    for group in &order.groups {
        let rows = &order.indices[group.rows.clone()];
        let next_shift = rows.iter().map(|&i| eta[i]).fold(shift, f64::max);
        if next_shift > shift {
            risk.rescale((shift - next_shift).exp());
            shift = next_shift;
        }
        removed.reset();
        for &i in rows {
            let scaled = eta[i] - shift;
            let weight = scaled.exp();
            let row = &x[i];
            risk.add(row, weight);
            if event[i] {
                loglik += scaled;
                for (gradient, &value) in grad.iter_mut().zip(row) {
                    *gradient += value;
                }
                if ties == Ties::Efron && group.events > 1 {
                    removed.add(row, weight);
                }
            }
        }
        if group.events == 0 {
            continue;
        }
        let (count, multiplicity) = match ties {
            | Ties::Breslow => (1, group.events as f64),
            | Ties::Efron => (group.events, 1.0),
        };
        for l in 0..count {
            let frac = l as f64 / group.events as f64;
            let denom = risk.weight - frac * removed.weight;
            if !denom.is_finite() || denom <= 0.0 {
                return (f64::NEG_INFINITY, grad, info);
            }
            loglik -= multiplicity * denom.ln();
            for a in 0..p {
                let mean_a = (risk.first[a] - frac * removed.first[a]) / denom;
                grad[a] -= multiplicity * mean_a;
                for b in 0..=a {
                    let cell = a * p + b;
                    let mean_b = (risk.first[b] - frac * removed.first[b]) / denom;
                    let cross = (risk.second[cell] - frac * removed.second[cell]) / denom;
                    info.data[cell] += multiplicity * (cross - mean_a * mean_b);
                }
            }
        }
    }
    for a in 0..p {
        for b in 0..a {
            info.data[b * p + a] = info.data[a * p + b];
        }
    }
    (loglik, grad, info)
}

/// Harrell's concordance index at the linear predictor `x . beta`.
///
/// A pair `(i, j)` with `event[i]` observed is comparable when `time[i] < time[j]`
/// or when `time[i] == time[j]` and observation `j` is censored (so `j`'s true
/// event time exceeds `time[i]`). Pairs of tied events are not comparable, and
/// tied risk scores contribute `0.5` to the concordant count.
#[allow(clippy::float_cmp)] // exact comparison identifies tied observation times
fn concordance(
    x: &[Vec<f64>],
    time: &[f64],
    event: &[bool],
    beta: &[f64],
) -> f64 {
    let n = time.len();
    let risk: Vec<f64> = x.iter().map(|row| dot(row, beta)).collect();
    let mut concordant = 0.0;
    let mut comparable = 0.0;
    for i in 0..n {
        if !event[i] {
            continue;
        }
        for j in 0..n {
            if i == j {
                continue;
            }
            let is_comparable = time[j] > time[i] || (time[j] == time[i] && !event[j]);
            if !is_comparable {
                continue;
            }
            comparable += 1.0;
            let tol = 16.0 * f64::EPSILON * (1.0 + risk[i].abs().max(risk[j].abs()));
            if (risk[i] - risk[j]).abs() <= tol {
                concordant += 0.5;
            } else if risk[i] > risk[j] {
                concordant += 1.0;
            }
        }
    }
    if comparable > 0.0 {
        concordant / comparable
    } else {
        0.0
    }
}

/// Fits a Cox proportional-hazards model by Newton-Raphson on the partial
/// likelihood, with step halving when a full step lowers the likelihood.
///
/// `x` holds one row of covariates per observation, `time` the observed times
/// and `event` whether the observation ended in the event. The iteration stops
/// when the largest absolute change of a coefficient falls below `tol` or
/// after `max_iter` iterations; the fit reports which happened through
/// [`CoxFit::converged`].
///
/// # Errors
/// Returns an error for empty or mismatched inputs, non-finite values, an empty
/// design matrix, a covariate column with zero variance, a sample in which
/// every observation is censored, a singular observed information matrix, or a
/// non-positive `max_iter` or `tol`.
pub fn cox_fit(
    x: &[Vec<f64>],
    time: &[f64],
    event: &[bool],
    ties: Ties,
    max_iter: usize,
    tol: f64,
) -> Result<CoxFit, String> {
    let p = validate(x, time, event)?;
    if max_iter == 0 {
        return Err("max_iter must be positive".to_string());
    }
    if !tol.is_finite() || tol <= 0.0 {
        return Err("tol must be positive and finite".to_string());
    }
    let order = RiskOrder::new(time, event);
    let mut beta = vec![0.0; p];
    let (mut loglik, mut grad, mut info) = evaluate(&beta, x, event, &order, ties);
    let mut converged = false;
    let mut iterations = 0usize;
    let mut candidate = beta.clone();
    for iter in 0..max_iter {
        iterations = iter + 1;
        let delta = solve(&info, &grad)
            .map_err(|_| "the observed information matrix is singular".to_string())?;
        let mut step = 1.0;
        let mut improved = false;
        for _ in 0..40 {
            for ((value, b), d) in candidate.iter_mut().zip(&beta).zip(&delta) {
                *value = b + step * d;
            }
            let (ll, next_grad, next_info) = evaluate(&candidate, x, event, &order, ties);
            if ll.is_finite() && ll >= loglik - 1e-12 {
                (loglik, grad, info) = (ll, next_grad, next_info);
                improved = true;
                break;
            }
            step *= 0.5;
        }
        if !improved {
            break;
        }
        let largest = delta.iter().map(|d| (step * d).abs()).fold(0.0, f64::max);
        std::mem::swap(&mut beta, &mut candidate);
        if largest < tol {
            converged = true;
            break;
        }
    }

    let factors =
        lu_factor(&info).map_err(|_| "the observed information matrix is singular".to_string())?;
    let inverse = factors.inverse();
    let mut se = vec![0.0; p];
    let mut z = vec![0.0; p];
    let mut p_value = vec![0.0; p];
    let normal = Distribution::Normal { mu: 0.0, sigma: 1.0 };
    for j in 0..p {
        let variance = inverse.at(j, j);
        if !variance.is_finite() || variance < 0.0 {
            return Err("the observed information matrix is not positive definite".to_string());
        }
        se[j] = variance.sqrt();
        if !se[j].is_finite() || se[j] <= 0.0 {
            return Err("a standard error is zero".to_string());
        }
        z[j] = beta[j] / se[j];
        p_value[j] = 2.0 * (1.0 - normal.cdf(z[j].abs()));
    }
    Ok(CoxFit {
        concordance: concordance(x, time, event, &beta),
        beta,
        se,
        z,
        p: p_value,
        loglik,
        iterations,
        converged,
    })
}

/// Relative hazards `exp(x . beta)` for new covariate rows.
///
/// # Errors
/// Returns an error when a row does not have one entry per fitted coefficient.
pub fn cox_predict(
    fit: &CoxFit,
    x: &[Vec<f64>],
) -> Result<Vec<f64>, String> {
    cox_predict_coefficients(&fit.beta, x)
}

/// Predicts relative hazards from fitted log-hazard coefficients.
///
/// # Errors
/// Fails when an observation has a different number of columns than the model.
pub fn cox_predict_coefficients(
    coefficients: &[f64],
    x: &[Vec<f64>],
) -> Result<Vec<f64>, String> {
    let mut out = Vec::with_capacity(x.len());
    for row in x {
        if row.len() != coefficients.len() {
            return Err("each row needs one entry per coefficient".to_string());
        }
        out.push(dot(row, coefficients).exp());
    }
    Ok(out)
}

#[cfg(test)]
mod tests {

    #[test]
    fn likelihood_uses_the_same_shift_in_the_numerator_and_risk_set() {
        use super::RiskOrder;
        use super::evaluate;
        let time = [1.0, 2.0, 2.0, 3.0, 4.0, 4.0, 6.0, 7.0, 7.0, 9.0];
        let events = [
            true, true, true, false, true, false, true, false, true, true,
        ];
        let x = vec![
            vec![0.0, 1.0],
            vec![1.0, 2.0],
            vec![2.0, 0.0],
            vec![0.0, -1.0],
            vec![1.0, 0.0],
            vec![2.0, 1.0],
            vec![0.0, 0.0],
            vec![1.0, -1.0],
            vec![2.0, 2.0],
            vec![1.0, 3.0],
        ];
        let shifted: Vec<_> = x.iter().map(|row| vec![row[0] + 2048.0, row[1]]).collect();
        let order = RiskOrder::new(&time, &events);
        for (mode, expected_beta, expected_likelihood, expected_initial) in [
            (
                Ties::Breslow,
                [-0.27517606604866644372059, 0.031559968384268934300931],
                -10.877790367190200294744,
                -11.744243680987878631446,
            ),
            (
                Ties::Efron,
                [-0.25949199911000958383773, 0.033591004225876807569548],
                -10.771430701428489793319,
                -11.605101320096243214325,
            ),
        ] {
            let finish = |rows: &[Vec<f64>]| evaluate(&[0.5, -0.25], rows, &events, &order, mode).0;
            assert!((finish(&x) - expected_initial).abs() < 2e-12);
            assert!((finish(&shifted) - expected_initial).abs() < 2e-12);
            let fitted = cox_fit(&x, &time, &events, mode, 30, 1e-10).expect("fit");
            assert!(fitted.converged);
            assert!((fitted.loglik - expected_likelihood).abs() < 2e-12);
            for (observed, expected) in fitted.beta.iter().zip(expected_beta) {
                assert!((observed - expected).abs() < 2e-10);
            }
        }
    }
    use super::Ties;
    use super::concordance;
    use super::cox_fit;
    use super::cox_predict;

    /// Eight rows, one binary covariate; the treated group fails earlier.
    fn binary_data() -> (Vec<Vec<f64>>, Vec<f64>, Vec<bool>) {
        let x = vec![
            vec![1.0],
            vec![1.0],
            vec![1.0],
            vec![1.0],
            vec![0.0],
            vec![0.0],
            vec![0.0],
            vec![0.0],
        ];
        let time = vec![1.0, 2.0, 4.0, 7.0, 3.0, 5.0, 6.0, 8.0];
        let event = vec![true, true, true, false, true, true, false, true];
        (x, time, event)
    }

    #[test]
    fn binary_covariate_gives_a_positive_effect() {
        let (x, time, event) = binary_data();
        let fit = cox_fit(&x, &time, &event, Ties::Breslow, 50, 1e-9).expect("fit");
        assert!(fit.converged, "the fit did not converge");
        assert!(fit.iterations < 20, "took {} iterations", fit.iterations);
        assert!(fit.beta[0] > 0.0, "beta = {}", fit.beta[0]);
        assert!(fit.se[0] > 0.0);
        assert!(fit.p[0] > 0.0 && fit.p[0] < 1.0);
        assert!(fit.loglik.is_finite());
        assert!(fit.concordance > 0.5, "c = {}", fit.concordance);
    }

    #[test]
    fn efron_and_breslow_agree_without_ties() {
        let (x, time, event) = binary_data();
        let breslow = cox_fit(&x, &time, &event, Ties::Breslow, 50, 1e-10).expect("fit");
        let efron = cox_fit(&x, &time, &event, Ties::Efron, 50, 1e-10).expect("fit");
        assert!((breslow.beta[0] - efron.beta[0]).abs() < 1e-8);
    }

    #[test]
    fn efron_and_breslow_differ_with_ties() {
        // Two tied events in the treated group at t = 2.
        let x = vec![
            vec![1.0],
            vec![1.0],
            vec![0.0],
            vec![0.0],
            vec![1.0],
            vec![1.0],
            vec![0.0],
            vec![0.0],
        ];
        let time = vec![2.0, 2.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
        let event = vec![true, true, true, true, false, true, true, false];
        let breslow = cox_fit(&x, &time, &event, Ties::Breslow, 50, 1e-10).expect("fit");
        let efron = cox_fit(&x, &time, &event, Ties::Efron, 50, 1e-10).expect("fit");
        assert!(breslow.beta[0].is_finite() && efron.beta[0].is_finite());
        assert!(
            (breslow.beta[0] - efron.beta[0]).abs() > 1e-6,
            "breslow {} efron {}",
            breslow.beta[0],
            efron.beta[0]
        );
    }

    #[test]
    fn predictions_and_validation() {
        let (x, time, event) = binary_data();
        let fit = cox_fit(&x, &time, &event, Ties::Efron, 50, 1e-9).expect("fit");
        let hazards = cox_predict(&fit, &[vec![0.0], vec![1.0]]).expect("predict");
        assert!((hazards[0] - 1.0).abs() < 1e-12);
        assert!((hazards[1] - fit.beta[0].exp()).abs() < 1e-12);
        assert!(cox_predict(&fit, &[vec![0.0, 1.0]]).is_err());

        assert!(cox_fit(&x, &time, &[false; 8], Ties::Efron, 50, 1e-9).is_err());
        let constant = vec![vec![1.0]; 8];
        assert!(cox_fit(&constant, &time, &event, Ties::Efron, 50, 1e-9).is_err());
        assert!(cox_fit(&x, &time, &event, Ties::Efron, 0, 1e-9).is_err());
        assert!(cox_fit(&x, &time, &event, Ties::Efron, 50, 0.0).is_err());
    }

    #[test]
    fn concordance_counts_tied_event_and_censored_pairs() {
        let time = vec![2.0, 2.0, 5.0];
        let event = vec![true, false, true];
        let beta = vec![1.0];

        // Event at t = 2 has higher risk than both the censored observation
        // tied at t = 2 and the later event at t = 5: 2/2 concordant.
        let x_concordant = vec![vec![2.0], vec![1.0], vec![0.0]];
        let c = concordance(&x_concordant, &time, &event, &beta);
        assert!((c - 1.0).abs() < 1e-15, "c = {c}");

        // Flipping the tied censored observation to have higher risk makes
        // (0, 1) discordant while (0, 2) stays concordant: 1/2 = 0.5.
        let x_mixed = vec![vec![2.0], vec![3.0], vec![0.0]];
        let c_mixed = concordance(&x_mixed, &time, &event, &beta);
        assert!((c_mixed - 0.5).abs() < 1e-15, "c = {c_mixed}");
    }

    #[test]
    fn concordance_excludes_tied_events_and_awards_half_for_risk_ties() {
        let x = vec![vec![1.0], vec![1.0], vec![1.0]];
        let time = vec![2.0, 2.0, 5.0];
        let event = vec![true, true, true];
        let beta = vec![1.0];
        // (0, 1) and (1, 0) are tied events (not comparable); (0, 2) and (1, 2)
        // are comparable with tied risk scores, earning 0.5 each.
        let c = concordance(&x, &time, &event, &beta);
        assert!((c - 0.5).abs() < 1e-15, "c = {c}");
    }

    #[test]
    fn concordance_returns_zero_when_no_pairs_are_comparable() {
        let x = vec![vec![1.0], vec![2.0]];
        let beta = vec![1.0];
        // Two tied events with no later observations: zero comparable pairs.
        assert_eq!(concordance(&x, &[2.0, 2.0], &[true, true], &beta), 0.0);
        // All observations censored: zero comparable pairs.
        assert_eq!(concordance(&x, &[1.0, 2.0], &[false, false], &beta), 0.0);
    }
}
