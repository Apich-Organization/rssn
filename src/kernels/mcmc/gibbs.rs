//! Gibbs sampling by inverse-CDF draws from the full conditionals.

use crate::kernels::mcmc::chain::Chain;
use crate::kernels::mcmc::chain::MhOptions;
use crate::kernels::mcmc::chain::is_kept;
use crate::kernels::mcmc::chain::total_iterations;
use crate::kernels::mcmc::rng::Pcg32;

/// Gibbs sampler driven by user-supplied full conditionals.
///
/// Each coordinate is refreshed in turn from its full conditional. A
/// conditional is a function of the current state and one uniform variate in
/// `(0, 1)`; it returns the new value of its coordinate, normally by inverting
/// the conditional CDF on a bracketed grid. One sweep over all coordinates is
/// one iteration. [`MhOptions::proposal_sd`] is not used, because Gibbs moves
/// have no proposal distribution.
///
/// The sampler never evaluates the joint density, so
/// [`Chain::log_posterior`] is filled with `NaN`.
///
/// # Errors
/// Returns an error when the number of conditionals does not match the state,
/// the state is empty, `n` or `thin` is zero, the iteration budget overflows,
/// or a conditional returns a non-finite value.
pub fn gibbs(
    conditionals: &[&dyn Fn(&[f64], f64) -> f64],
    init: &[f64],
    n: usize,
    opts: &MhOptions,
) -> Result<Chain, String> {
    let total = total_iterations(n, opts.burn_in, opts.thin)?;
    if init.is_empty() {
        return Err("the initial state is empty".to_string());
    }
    if conditionals.len() != init.len() {
        return Err("one conditional is required per coordinate".to_string());
    }
    let mut rng = Pcg32::seeded(opts.seed);
    let mut state = init.to_vec();
    let mut chain = Chain {
        draws: Vec::with_capacity(n),
        log_posterior: Vec::with_capacity(n),
        accepted: 0,
        proposed: 0,
    };
    for iteration in 0..total {
        for (coordinate, conditional) in conditionals.iter().enumerate() {
            let uniform = rng.next_f64();
            let value = conditional(&state, uniform);
            if !value.is_finite() {
                return Err("a conditional returned a non-finite value".to_string());
            }
            state[coordinate] = value;
        }
        chain.proposed += 1;
        chain.accepted += 1;
        if is_kept(iteration, opts.burn_in, opts.thin) {
            chain.draws.push(state.clone());
            chain.log_posterior.push(f64::NAN);
        }
    }
    Ok(chain)
}

#[cfg(test)]
mod tests {
    use super::gibbs;
    use crate::kernels::distributions::Distribution;
    use crate::kernels::mcmc::chain::MhOptions;

    /// Standard-normal quantile, guarded against the open ends of `(0, 1)`.
    fn normal_quantile(u: f64) -> f64 {
        Distribution::Normal { mu: 0.0, sigma: 1.0 }
            .quantile(u.clamp(1e-12, 1.0 - 1e-12))
            .unwrap_or(0.0)
    }

    #[test]
    fn bivariate_gaussian_conditionals_recover_the_means() {
        // X and Y are standard normal with correlation 0.8, so
        // X | Y = y ~ N(0.8 y, 0.36) and Y | X = x ~ N(0.8 x, 0.36).
        let conditional_x = |state: &[f64], u: f64| 0.8 * state[1] + 0.6 * normal_quantile(u);
        let conditional_y = |state: &[f64], u: f64| 0.8 * state[0] + 0.6 * normal_quantile(u);
        let conditionals: [&dyn Fn(&[f64], f64) -> f64; 2] = [&conditional_x, &conditional_y];
        let options = MhOptions {
            proposal_sd: Vec::new(),
            burn_in: 1_000,
            thin: 1,
            seed: 3,
        };
        // 40 000 sweeps; with a lag-one correlation near 0.8 the effective
        // sample size is about 4 500, so the mean has a standard error of
        // roughly 0.015 and 0.05 is a safe tolerance.
        let chain = gibbs(&conditionals, &[0.0, 0.0], 40_000, &options).expect("run");
        let n = chain.draws.len() as f64;
        let mean_x = chain.draws.iter().map(|d| d[0]).sum::<f64>() / n;
        let mean_y = chain.draws.iter().map(|d| d[1]).sum::<f64>() / n;
        assert!(mean_x.abs() < 0.05, "mean x = {mean_x}");
        assert!(mean_y.abs() < 0.05, "mean y = {mean_y}");
        assert!(chain.log_posterior.iter().all(|v| v.is_nan()));
    }

    #[test]
    fn invalid_inputs_are_rejected() {
        let conditional = |_: &[f64], _: f64| 0.0;
        let one: [&dyn Fn(&[f64], f64) -> f64; 1] = [&conditional];
        let options = MhOptions {
            proposal_sd: Vec::new(),
            burn_in: 0,
            thin: 1,
            seed: 1,
        };
        assert!(gibbs(&one, &[0.0, 0.0], 10, &options).is_err());
        assert!(gibbs(&one, &[], 10, &options).is_err());
        assert!(gibbs(&one, &[0.0], 0, &options).is_err());
        let infinite = |_: &[f64], _: f64| f64::INFINITY;
        let bad: [&dyn Fn(&[f64], f64) -> f64; 1] = [&infinite];
        assert!(gibbs(&bad, &[0.0], 10, &options).is_err());
    }
}
