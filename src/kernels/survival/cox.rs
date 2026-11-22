//! Cox proportional-hazards regression.

use crate::kernels::dense::Mat;
use crate::kernels::dense::lu_factor;
use crate::kernels::dense::solve;
use crate::kernels::distributions::Distribution;

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

/// One event time: the indices of the tied events and of the risk set.
struct Block {
    events: Vec<usize>,
    risk: Vec<usize>,
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
        let column: Vec<f64> = x.iter().map(|row| row[j]).collect();
        let mean = column.iter().sum::<f64>() / n as f64;
        let variance = column.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / n as f64;
        if !variance.is_finite() || variance <= 0.0 {
            return Err("a covariate column has zero variance".to_string());
        }
    }
    Ok(p)
}

/// Builds one block per distinct event time.
#[allow(clippy::float_cmp)] // exact comparison groups tied observation times
fn blocks(
    time: &[f64],
    event: &[bool],
) -> Vec<Block> {
    let n = time.len();
    let mut times: Vec<f64> = time
        .iter()
        .copied()
        .zip(event.iter().copied())
        .filter(|pair| pair.1)
        .map(|pair| pair.0)
        .collect();
    times.sort_unstable_by(f64::total_cmp);
    times.dedup();
    times
        .iter()
        .map(|&t| {
            Block {
                events: (0..n).filter(|&i| event[i] && time[i] == t).collect(),
                risk: (0..n).filter(|&i| time[i] >= t).collect(),
            }
        })
        .collect()
}

/// Log partial likelihood, gradient and observed information at `beta`.
///
/// The exponentials are shifted by the largest linear predictor of each risk
/// set so that the accumulation cannot overflow; the shift cancels from the
/// ratios and is added back to the log likelihood.
fn evaluate(
    beta: &[f64],
    x: &[Vec<f64>],
    blocks: &[Block],
    p: usize,
    ties: Ties,
) -> (f64, Vec<f64>, Mat) {
    let eta: Vec<f64> = x.iter().map(|row| dot(row, beta)).collect();
    let mut loglik = 0.0;
    let mut grad = vec![0.0; p];
    let mut info = Mat::zeros(p, p);
    let mut s1 = vec![0.0; p];
    let mut s2 = vec![0.0; p * p];
    let mut d1 = vec![0.0; p];
    let mut d2 = vec![0.0; p * p];
    for block in blocks {
        let shift = block
            .risk
            .iter()
            .map(|&i| eta[i])
            .fold(f64::NEG_INFINITY, f64::max);
        let mut s0 = 0.0;
        s1.fill(0.0);
        s2.fill(0.0);
        d1.fill(0.0);
        d2.fill(0.0);
        let mut d0 = 0.0;
        for &i in &block.risk {
            let w = (eta[i] - shift).exp();
            s0 += w;
            let row = &x[i];
            for a in 0..p {
                s1[a] += w * row[a];
                for b in 0..p {
                    s2[a * p + b] += w * row[a] * row[b];
                }
            }
        }
        for &i in &block.events {
            let w = (eta[i] - shift).exp();
            d0 += w;
            let row = &x[i];
            for a in 0..p {
                d1[a] += w * row[a];
                for b in 0..p {
                    d2[a * p + b] += w * row[a] * row[b];
                }
            }
        }
        let count = block.events.len();
        for l in 0..count {
            let frac = match ties {
                | Ties::Breslow => 0.0,
                | Ties::Efron => l as f64 / count as f64,
            };
            let denom = s0 - frac * d0;
            if !denom.is_finite() || denom <= 0.0 {
                // Not reachable for valid input; force the step to be rejected.
                loglik = f64::NEG_INFINITY;
                continue;
            }
            loglik -= denom.ln();
            for a in 0..p {
                let num = s1[a] - frac * d1[a];
                grad[a] -= num / denom;
                for b in 0..p {
                    let num2 = s2[a * p + b] - frac * d2[a * p + b];
                    let value = info.at(a, b) + num2 / denom
                        - num * (s1[b] - frac * d1[b]) / (denom * denom);
                    info.set(a, b, value);
                }
            }
        }
        loglik += block.events.iter().map(|&i| eta[i]).sum::<f64>() + count as f64 * shift;
        for &i in &block.events {
            let row = &x[i];
            for a in 0..p {
                grad[a] += row[a];
            }
        }
    }
    (loglik, grad, info)
}

/// Harrell's concordance index at the linear predictor `x . beta`.
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
            if i == j || time[j] <= time[i] {
                continue;
            }
            comparable += 1.0;
            if risk[i] > risk[j] {
                concordant += 1.0;
            } else if (risk[i] - risk[j]).abs() < f64::EPSILON {
                concordant += 0.5;
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
    let blocks = blocks(time, event);
    let mut beta = vec![0.0; p];
    let (mut loglik, _, _) = evaluate(&beta, x, &blocks, p, ties);
    let mut converged = false;
    let mut iterations = 0usize;
    for iter in 0..max_iter {
        iterations = iter + 1;
        let (current, grad, info) = evaluate(&beta, x, &blocks, p, ties);
        loglik = current;
        let delta = solve(&info, &grad)
            .map_err(|_| "the observed information matrix is singular".to_string())?;
        let mut step = 1.0;
        let mut improved = false;
        let mut candidate = beta.clone();
        let mut candidate_loglik = current;
        for _ in 0..40 {
            for ((value, b), d) in candidate.iter_mut().zip(&beta).zip(&delta) {
                *value = b + step * d;
            }
            let (ll, _, _) = evaluate(&candidate, x, &blocks, p, ties);
            if ll.is_finite() && ll >= current - 1e-12 {
                candidate_loglik = ll;
                improved = true;
                break;
            }
            step *= 0.5;
        }
        if !improved {
            break;
        }
        let largest = delta.iter().map(|d| (step * d).abs()).fold(0.0, f64::max);
        beta = candidate;
        loglik = candidate_loglik;
        if largest < tol {
            converged = true;
            break;
        }
    }

    let (_, _, info) = evaluate(&beta, x, &blocks, p, ties);
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
    let p = fit.beta.len();
    let mut out = Vec::with_capacity(x.len());
    for row in x {
        if row.len() != p {
            return Err("each row needs one entry per coefficient".to_string());
        }
        out.push(dot(row, &fit.beta).exp());
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::Ties;
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
}
