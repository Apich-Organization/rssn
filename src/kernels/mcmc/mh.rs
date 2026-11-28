//! Metropolis-Hastings samplers with a symmetric Gaussian proposal.

use crate::kernels::mcmc::chain::Chain;
use crate::kernels::mcmc::chain::MhOptions;
use crate::kernels::mcmc::chain::is_kept;
use crate::kernels::mcmc::chain::total_iterations;
use crate::kernels::mcmc::rng::Pcg32;

/// Metropolis-Hastings sampler with a symmetric Gaussian random walk.
///
/// Each iteration proposes one joint move `x + proposal_sd * z` with `z`
/// standard normal per coordinate and accepts it with probability
/// `min(1, exp(log_target(new) - log_target(old)))`. Proposals whose target
/// density is not a number are rejected. The returned chain holds `n` draws
/// taken every `thin` iterations after `burn_in` discarded ones.
///
/// # Errors
/// Returns an error when the initial state or `proposal_sd` is empty or
/// mismatched, a proposal scale is not positive and finite, `n` or `thin` is
/// zero, the iteration budget overflows, or the target log density is not
/// finite at the initial state.
pub fn metropolis_hastings(
    log_target: &dyn Fn(&[f64]) -> f64,
    init: &[f64],
    n: usize,
    opts: &MhOptions,
) -> Result<Chain, String> {
    let total = total_iterations(n, opts.burn_in, opts.thin)?;
    if init.is_empty() {
        return Err("the initial state is empty".to_string());
    }
    if opts.proposal_sd.len() != init.len() {
        return Err("proposal_sd must have one entry per coordinate".to_string());
    }
    if opts.proposal_sd.iter().any(|s| !s.is_finite() || *s <= 0.0) {
        return Err("proposal_sd entries must be positive and finite".to_string());
    }
    let mut current_log = log_target(init);
    if !current_log.is_finite() {
        return Err("the target log density is not finite at the initial state".to_string());
    }

    let mut rng = Pcg32::seeded(opts.seed);
    let mut state = init.to_vec();
    let mut proposal = init.to_vec();
    let mut chain = Chain {
        draws: Vec::with_capacity(n),
        log_posterior: Vec::with_capacity(n),
        accepted: 0,
        proposed: 0,
    };
    for iteration in 0..total {
        for ((value, &x), &scale) in proposal.iter_mut().zip(&state).zip(&opts.proposal_sd) {
            *value = x + scale * rng.normal();
        }
        let proposed_log = log_target(&proposal);
        chain.proposed += 1;
        // A NaN difference fails the comparison and rejects the move; a
        // density of `-inf` is rejected because the exponential vanishes.
        let log_u = rng.next_f64().max(f64::MIN_POSITIVE).ln();
        if log_u < proposed_log - current_log {
            state.copy_from_slice(&proposal);
            current_log = proposed_log;
            chain.accepted += 1;
        }
        if is_kept(iteration, opts.burn_in, opts.thin) {
            chain.draws.push(state.clone());
            chain.log_posterior.push(current_log);
        }
    }
    Ok(chain)
}

/// Random-walk Metropolis sampler.
///
/// This is the same algorithm as [`metropolis_hastings`] with the per-coordinate
/// proposal scales taken as given: no adaptation of the step sizes is
/// performed, so the chain stays a time-homogeneous Markov chain.
///
/// # Errors
/// Returns the errors of [`metropolis_hastings`].
pub fn random_walk_metropolis(
    log_target: &dyn Fn(&[f64]) -> f64,
    init: &[f64],
    n: usize,
    opts: &MhOptions,
) -> Result<Chain, String> {
    metropolis_hastings(log_target, init, n, opts)
}

#[cfg(test)]
mod tests {
    use super::metropolis_hastings;
    use super::random_walk_metropolis;
    use crate::kernels::mcmc::chain::MhOptions;

    /// Log density of a standard normal.
    fn standard_normal(x: &[f64]) -> f64 {
        match x.first() {
            | Some(&v) => -0.5 * v * v,
            | None => f64::NEG_INFINITY,
        }
    }

    fn options() -> MhOptions {
        MhOptions {
            proposal_sd: vec![2.4],
            burn_in: 1_000,
            thin: 1,
            seed: 11,
        }
    }

    #[test]
    fn normal_target_recovers_its_moments() {
        // 40 000 draws with about 65 % acceptance: the standard error of the
        // mean is roughly 0.01 and of the variance roughly 0.015, so 0.03 and
        // 0.06 are comfortable tolerances.
        let chain = metropolis_hastings(&standard_normal, &[0.0], 40_000, &options()).expect("run");
        assert_eq!(chain.draws.len(), 40_000);
        assert_eq!(chain.log_posterior.len(), 40_000);
        assert_eq!(chain.proposed, 41_000);
        let n = chain.draws.len() as f64;
        let mean = chain.draws.iter().map(|d| d[0]).sum::<f64>() / n;
        let variance = chain.draws.iter().map(|d| d[0] * d[0]).sum::<f64>() / n - mean * mean;
        assert!(mean.abs() < 0.03, "mean = {mean}");
        assert!((variance - 1.0).abs() < 0.06, "variance = {variance}");
        let rate = chain.accepted as f64 / chain.proposed as f64;
        assert!(rate > 0.4 && rate < 0.9, "acceptance rate = {rate}");
    }

    #[test]
    fn random_walk_alias_matches_metropolis_hastings() {
        let direct = metropolis_hastings(&standard_normal, &[0.0], 500, &options()).expect("run");
        let alias = random_walk_metropolis(&standard_normal, &[0.0], 500, &options()).expect("run");
        assert_eq!(direct, alias);
    }

    #[test]
    fn invalid_options_are_rejected() {
        assert!(metropolis_hastings(&standard_normal, &[], 100, &options()).is_err());
        let mut bad = options();
        bad.proposal_sd = vec![0.0];
        assert!(metropolis_hastings(&standard_normal, &[0.0], 100, &bad).is_err());
        let mut bad = options();
        bad.proposal_sd = vec![1.0, 1.0];
        assert!(metropolis_hastings(&standard_normal, &[0.0], 100, &bad).is_err());
        let mut bad = options();
        bad.thin = 0;
        assert!(metropolis_hastings(&standard_normal, &[0.0], 100, &bad).is_err());
        let nowhere = |_: &[f64]| f64::NAN;
        assert!(metropolis_hastings(&nowhere, &[0.0], 100, &options()).is_err());
    }
}
