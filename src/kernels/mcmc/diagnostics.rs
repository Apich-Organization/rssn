//! Convergence and mixing diagnostics for Markov chains.

use crate::kernels::mcmc::chain::mean_variance;

/// Sum of lagged products `(x_t - mean) (x_{t+lag} - mean)`.
fn lagged_product(
    chain: &[f64],
    mean: f64,
    lag: usize,
) -> f64 {
    let split = chain.len().saturating_sub(lag);
    chain
        .iter()
        .take(split)
        .zip(chain.iter().skip(lag))
        .map(|(a, b)| (a - mean) * (b - mean))
        .sum()
}

/// Sample autocorrelation of `chain` at lags `0..=max_lag`.
///
/// The estimate uses the usual biased denominator, the full sum of squared
/// deviations from the mean, so the sequence is positive semidefinite. Lag 0 is
/// always 1. Entries that cannot be estimated — a chain with fewer than two
/// draws, a lag of at least the chain length, or a constant chain — are `NaN`.
#[must_use]
pub fn autocorrelation(
    chain: &[f64],
    max_lag: usize,
) -> Vec<f64> {
    let n = chain.len();
    let mut out = vec![f64::NAN; max_lag + 1];
    if n < 2 {
        return out;
    }
    let mean = chain.iter().sum::<f64>() / n as f64;
    let total: f64 = chain.iter().map(|x| (x - mean) * (x - mean)).sum();
    if !total.is_finite() || total <= 0.0 {
        return out;
    }
    if let Some(first) = out.first_mut() {
        *first = 1.0;
    }
    for lag in 1..=max_lag {
        if lag >= n {
            break;
        }
        let mut numerator = 0.0;
        for t in 0..n - lag {
            numerator += (chain[t] - mean) * (chain[t + lag] - mean);
        }
        out[lag] = numerator / total;
    }
    out
}

/// Effective sample size by the initial positive sequence estimator.
///
/// The autocorrelations are summed in adjacent pairs and the sum is truncated
/// at the first non-positive pair, which keeps the estimate from being dragged
/// down by the noise of the tail. Lags are computed only until that truncation,
/// so a well-mixing chain costs two passes over the data. Returns the number of
/// independent draws with the same information content as the chain.
///
/// # Errors
/// Returns an error for fewer than four draws or a constant chain.
pub fn effective_sample_size(chain: &[f64]) -> Result<f64, String> {
    let n = chain.len();
    if n < 4 {
        return Err("at least four draws are required".to_string());
    }
    let (mean, variance) =
        mean_variance(chain).ok_or_else(|| "at least four draws are required".to_string())?;
    if !variance.is_finite() || variance <= 0.0 {
        return Err("the chain has zero variance".to_string());
    }
    let total = variance * (n - 1) as f64;
    let mut sum = 0.0;
    let mut lag = 1usize;
    while lag + 1 < n {
        // Pairs start at lag one: the leading 1 in the variance inflation
        // factor is the lag-zero autocorrelation itself.
        let pair =
            (lagged_product(chain, mean, lag) + lagged_product(chain, mean, lag + 1)) / total;
        if !pair.is_finite() || pair <= 0.0 {
            break;
        }
        sum += pair;
        lag += 2;
    }
    Ok(n as f64 / (1.0 + 2.0 * sum))
}

/// Split-Rhat potential scale reduction factor of several chains.
///
/// Every chain is cut in half, the two halves are treated as independent
/// chains, and the ratio of the pooled variance estimate to the mean
/// within-chain variance is returned. Values near one indicate that the chains
/// have mixed; values well above one indicate that they disagree.
///
/// # Errors
/// Returns an error for fewer than two chains, chains of different lengths,
/// chains shorter than four draws, or a zero within-chain variance.
pub fn rhat(chains: &[Vec<f64>]) -> Result<f64, String> {
    if chains.len() < 2 {
        return Err("at least two chains are required".to_string());
    }
    let n = chains.first().map_or(0, Vec::len);
    if chains.iter().any(|chain| chain.len() != n) {
        return Err("chains must have the same length".to_string());
    }
    let half = n / 2;
    if half < 2 {
        return Err("each chain needs at least four draws".to_string());
    }
    let mut means = Vec::with_capacity(2 * chains.len());
    let mut variances = Vec::with_capacity(2 * chains.len());
    for chain in chains {
        for segment in [&chain[..half], &chain[n - half..]] {
            let (mean, variance) =
                mean_variance(segment).ok_or_else(|| "a segment is too short".to_string())?;
            means.push(mean);
            variances.push(variance);
        }
    }
    let count = means.len() as f64;
    let within = variances.iter().sum::<f64>() / count;
    if !within.is_finite() || within <= 0.0 {
        return Err("the within-chain variance is zero".to_string());
    }
    let grand = means.iter().sum::<f64>() / count;
    let between = means.iter().map(|m| (m - grand) * (m - grand)).sum::<f64>() / (count - 1.0);
    let length = half as f64;
    // `between` is the variance of the segment means, which already carries the
    // factor `length` of the classical `B = n * Var(means)`.
    let pooled = (length - 1.0) / length * within + between;
    Ok((pooled / within).sqrt())
}

/// Geweke z-score comparing the first 10 % with the last 50 % of a chain.
///
/// The two segments are compared with a difference of means standardised by
/// their sample variances. Values far from zero in absolute value suggest that
/// the chain has not converged.
///
/// # Errors
/// Returns an error for fewer than twenty draws or a zero variance.
pub fn geweke(chain: &[f64]) -> Result<f64, String> {
    let n = chain.len();
    if n < 20 {
        return Err("at least twenty draws are required".to_string());
    }
    let first_len = n / 10;
    let last_len = n / 2;
    if first_len < 2 {
        return Err("at least twenty draws are required".to_string());
    }
    let first = chain.get(..first_len).unwrap_or(&[]);
    let last = chain.get(n - last_len..).unwrap_or(&[]);
    let (mean_first, variance_first) =
        mean_variance(first).ok_or_else(|| "the first segment is too short".to_string())?;
    let (mean_last, variance_last) =
        mean_variance(last).ok_or_else(|| "the last segment is too short".to_string())?;
    let denominator = (variance_first / first_len as f64 + variance_last / last_len as f64).sqrt();
    if !denominator.is_finite() || denominator <= 0.0 {
        return Err("the chain has zero variance".to_string());
    }
    Ok((mean_first - mean_last) / denominator)
}

#[cfg(test)]
mod tests {
    use super::autocorrelation;
    use super::effective_sample_size;
    use super::geweke;
    use super::rhat;
    use crate::kernels::mcmc::rng::Pcg32;

    /// Independent standard-normal draws.
    fn normal_draws(
        n: usize,
        seed: u64,
    ) -> Vec<f64> {
        let mut rng = Pcg32::seeded(seed);
        (0..n).map(|_| rng.normal()).collect()
    }

    #[test]
    fn independent_draws_have_negligible_autocorrelation() {
        // 20 000 draws put the standard error of each autocorrelation at 0.007,
        // so 0.03 is a four-sigma bound.
        let draws = normal_draws(20_000, 1);
        let acf = autocorrelation(&draws, 10);
        assert_eq!(acf.len(), 11);
        assert!((acf[0] - 1.0).abs() < 1e-12);
        for &value in &acf[1..] {
            assert!(value.abs() < 0.03, "autocorrelation {value}");
        }
        let empty = autocorrelation(&[], 3);
        assert_eq!(empty.len(), 4);
        assert!(empty.iter().all(|v| v.is_nan()));
    }

    #[test]
    fn effective_sample_size_of_independent_draws_is_the_draw_count() {
        // The initial positive sequence stops at the first negative pair, so
        // the estimate stays within a few percent of the draw count; 25 % is
        // the tolerance asked for by the interface contract.
        let draws = normal_draws(4_000, 2);
        let ess = effective_sample_size(&draws).expect("ess");
        let n = draws.len() as f64;
        assert!((ess - n).abs() < 0.25 * n, "ess = {ess}, n = {n}");
        assert!(effective_sample_size(&[1.0, 1.0, 1.0, 1.0]).is_err());
        assert!(effective_sample_size(&[1.0, 2.0]).is_err());
    }

    #[test]
    fn rhat_detects_identical_and_disagreeing_chains() {
        let draws = normal_draws(1_000, 3);
        let identical = rhat(&[draws.clone(), draws.clone()]).expect("rhat");
        assert!(
            (0.99..=1.01).contains(&identical),
            "identical chains give {identical}"
        );
        let shifted: Vec<f64> = draws.iter().map(|v| v + 3.0).collect();
        let disagreeing = rhat(&[draws, shifted]).expect("rhat");
        assert!(disagreeing > 1.2, "shifted chains give {disagreeing}");
        assert!(rhat(&[normal_draws(100, 4)]).is_err());
        assert!(rhat(&[normal_draws(100, 5), normal_draws(50, 6)]).is_err());
    }

    #[test]
    fn geweke_flags_a_shift_in_the_mean() {
        let stationary = normal_draws(4_000, 7);
        let z = geweke(&stationary).expect("z");
        assert!(z.abs() < 3.0, "stationary z = {z}");
        let drifting: Vec<f64> = stationary
            .iter()
            .enumerate()
            .map(|(i, &v)| if i < 2_000 { v } else { v + 5.0 })
            .collect();
        let shifted = geweke(&drifting).expect("z");
        assert!(shifted.abs() > 3.0, "drifting z = {shifted}");
        assert!(geweke(&[1.0, 2.0]).is_err());
    }
}
