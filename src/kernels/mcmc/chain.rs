//! Chain storage and shared option types for the samplers.

/// A Markov chain produced by one of the samplers.
#[derive(Clone, Debug, PartialEq)]
pub struct Chain {
    /// Post-burn-in, thinned draws; one inner vector per draw.
    pub draws: Vec<Vec<f64>>,
    /// Target log density at each stored draw.
    pub log_posterior: Vec<f64>,
    /// Number of accepted proposals, including those during burn-in.
    pub accepted: usize,
    /// Number of proposals attempted, including those during burn-in.
    pub proposed: usize,
}

/// Tuning parameters for the Metropolis-Hastings and Gibbs samplers.
#[derive(Clone, Debug, PartialEq)]
pub struct MhOptions {
    /// Per-coordinate standard deviation of the Gaussian proposal.
    pub proposal_sd: Vec<f64>,
    /// Number of initial iterations discarded.
    pub burn_in: usize,
    /// Keep every `thin`-th draw after burn-in.
    pub thin: usize,
    /// Seed of the internal generator.
    pub seed: u64,
}

/// Total number of iterations needed for `n` stored draws after burn-in.
pub(super) fn total_iterations(
    n: usize,
    burn_in: usize,
    thin: usize,
) -> Result<usize, String> {
    if n == 0 {
        return Err("the number of draws must be positive".to_string());
    }
    if thin == 0 {
        return Err("thin must be positive".to_string());
    }
    let kept = n
        .checked_mul(thin)
        .ok_or_else(|| "the draw budget overflows".to_string())?;
    burn_in
        .checked_add(kept)
        .ok_or_else(|| "the draw budget overflows".to_string())
}

/// Whether an iteration index belongs to the retained part of the chain.
pub(super) fn is_kept(
    iteration: usize,
    burn_in: usize,
    thin: usize,
) -> bool {
    iteration >= burn_in && (iteration - burn_in).is_multiple_of(thin)
}

/// Mean and unbiased variance of a slice, or `None` when it has fewer than two
/// entries.
pub(super) fn mean_variance(values: &[f64]) -> Option<(f64, f64)> {
    let n = values.len();
    if n < 2 {
        return None;
    }
    let mean = values.iter().sum::<f64>() / n as f64;
    let variance = values.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / (n - 1) as f64;
    Some((mean, variance))
}

#[cfg(test)]
mod tests {
    use super::is_kept;
    use super::mean_variance;
    use super::total_iterations;

    #[test]
    fn iteration_budget_counts_burn_in_and_thinning() {
        assert_eq!(total_iterations(10, 5, 2).expect("budget"), 25);
        assert!(total_iterations(0, 5, 2).is_err());
        assert!(total_iterations(10, 5, 0).is_err());
    }

    #[test]
    fn thinning_selects_every_nth_iteration_after_burn_in() {
        let kept: Vec<usize> = (0..12).filter(|&i| is_kept(i, 4, 3)).collect();
        assert_eq!(kept, vec![4, 7, 10]);
    }

    #[test]
    fn moments_use_the_unbiased_variance() {
        let (mean, variance) = mean_variance(&[1.0, 2.0, 3.0, 4.0]).expect("moments");
        assert!((mean - 2.5).abs() < 1e-15);
        assert!((variance - 5.0 / 3.0).abs() < 1e-15);
        assert!(mean_variance(&[1.0]).is_none());
    }
}
