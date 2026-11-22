//! Hamiltonian Monte Carlo with a standard-normal momentum.

use crate::kernels::mcmc::chain::Chain;
use crate::kernels::mcmc::chain::is_kept;
use crate::kernels::mcmc::chain::total_iterations;
use crate::kernels::mcmc::rng::Pcg32;

/// Tuning parameters for the Hamiltonian Monte Carlo sampler.
#[derive(Clone, Debug, PartialEq)]
pub struct HmcOptions {
    /// Leapfrog step size.
    pub step_size: f64,
    /// Number of leapfrog steps per proposal.
    pub leapfrog_steps: usize,
    /// Number of initial iterations discarded.
    pub burn_in: usize,
    /// Keep every `thin`-th draw after burn-in.
    pub thin: usize,
    /// Seed of the internal generator.
    pub seed: u64,
}

/// Hamiltonian Monte Carlo sampler.
///
/// A standard-normal momentum is drawn at every iteration, the Hamiltonian
/// dynamics of the target are integrated with the leapfrog scheme for
/// `leapfrog_steps` steps of length `step_size`, and the endpoint is accepted
/// with the Metropolis probability `min(1, exp(H_old - H_new))`. The returned
/// chain holds `n` draws taken every `thin` iterations after `burn_in`
/// discarded ones.
///
/// # Errors
/// Returns an error when the initial state is empty, the step size is not
/// positive and finite, the number of leapfrog steps, `n` or `thin` is zero,
/// the iteration budget overflows, the target log density is not finite at the
/// initial state, or the gradient has the wrong length or a non-finite entry.
pub fn hmc(
    log_target: &dyn Fn(&[f64]) -> f64,
    grad_log_target: &dyn Fn(&[f64]) -> Vec<f64>,
    init: &[f64],
    n: usize,
    opts: &HmcOptions,
) -> Result<Chain, String> {
    let total = total_iterations(n, opts.burn_in, opts.thin)?;
    if init.is_empty() {
        return Err("the initial state is empty".to_string());
    }
    if !opts.step_size.is_finite() || opts.step_size <= 0.0 {
        return Err("the step size must be positive and finite".to_string());
    }
    if opts.leapfrog_steps == 0 {
        return Err("at least one leapfrog step is required".to_string());
    }
    let dim = init.len();
    let gradient = |q: &[f64]| -> Result<Vec<f64>, String> {
        let g = grad_log_target(q);
        if g.len() != dim {
            return Err("the gradient has the wrong length".to_string());
        }
        if g.iter().any(|v| !v.is_finite()) {
            return Err("the gradient is not finite".to_string());
        }
        Ok(g)
    };
    let mut current_log = log_target(init);
    if !current_log.is_finite() {
        return Err("the target log density is not finite at the initial state".to_string());
    }
    let _ = gradient(init)?;

    let mut rng = Pcg32::seeded(opts.seed);
    let mut state = init.to_vec();
    let mut chain = Chain {
        draws: Vec::with_capacity(n),
        log_posterior: Vec::with_capacity(n),
        accepted: 0,
        proposed: 0,
    };
    let eps = opts.step_size;
    for iteration in 0..total {
        let mut momentum: Vec<f64> = (0..dim).map(|_| rng.normal()).collect();
        let mut position = state.clone();
        let mut grad = gradient(&position)?;
        let kinetic: f64 = 0.5 * momentum.iter().map(|p| p * p).sum::<f64>();
        let energy_old = -current_log + kinetic;
        for (p, g) in momentum.iter_mut().zip(&grad) {
            *p += 0.5 * eps * g;
        }
        for step in 0..opts.leapfrog_steps {
            for (q, p) in position.iter_mut().zip(&momentum) {
                *q += eps * p;
            }
            grad = gradient(&position)?;
            if step + 1 < opts.leapfrog_steps {
                for (p, g) in momentum.iter_mut().zip(&grad) {
                    *p += eps * g;
                }
            }
        }
        for (p, g) in momentum.iter_mut().zip(&grad) {
            *p += 0.5 * eps * g;
        }
        let proposed_log = log_target(&position);
        let kinetic_new: f64 = 0.5 * momentum.iter().map(|p| p * p).sum::<f64>();
        let energy_new = -proposed_log + kinetic_new;
        chain.proposed += 1;
        // A NaN energy difference fails the comparison and rejects the move.
        let log_u = rng.next_f64().max(f64::MIN_POSITIVE).ln();
        if log_u < energy_old - energy_new {
            state = position;
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

#[cfg(test)]
mod tests {
    use super::HmcOptions;
    use super::hmc;

    /// Unnormalised log density of `N(0, [[1, 0.8], [0.8, 1]])`.
    fn correlated_log(x: &[f64]) -> f64 {
        let (a, b) = (x[0], x[1]);
        -0.5 * (25.0 / 9.0 * a * a - 40.0 / 9.0 * a * b + 25.0 / 9.0 * b * b)
    }

    /// Gradient of [`correlated_log`].
    fn correlated_grad(x: &[f64]) -> Vec<f64> {
        let (a, b) = (x[0], x[1]);
        vec![
            -25.0 / 9.0 * a + 20.0 / 9.0 * b,
            20.0 / 9.0 * a - 25.0 / 9.0 * b,
        ]
    }

    fn options() -> HmcOptions {
        HmcOptions {
            step_size: 0.25,
            leapfrog_steps: 10,
            burn_in: 1_000,
            thin: 1,
            seed: 5,
        }
    }

    #[test]
    fn correlated_gaussian_recovers_its_covariance() {
        // 20 000 draws; the standard error of a covariance entry is about
        // 0.01 for this autocorrelation, so 0.05 is a comfortable tolerance.
        let chain = hmc(
            &correlated_log,
            &correlated_grad,
            &[0.0, 0.0],
            20_000,
            &options(),
        )
        .expect("run");
        assert_eq!(chain.draws.len(), 20_000);
        let n = chain.draws.len() as f64;
        let mean_a = chain.draws.iter().map(|d| d[0]).sum::<f64>() / n;
        let mean_b = chain.draws.iter().map(|d| d[1]).sum::<f64>() / n;
        let var_a = chain.draws.iter().map(|d| d[0] * d[0]).sum::<f64>() / n - mean_a * mean_a;
        let var_b = chain.draws.iter().map(|d| d[1] * d[1]).sum::<f64>() / n - mean_b * mean_b;
        let cov = chain.draws.iter().map(|d| d[0] * d[1]).sum::<f64>() / n - mean_a * mean_b;
        assert!((var_a - 1.0).abs() < 0.05, "var_a = {var_a}");
        assert!((var_b - 1.0).abs() < 0.05, "var_b = {var_b}");
        assert!((cov - 0.8).abs() < 0.05, "cov = {cov}");
        let rate = chain.accepted as f64 / chain.proposed as f64;
        assert!(rate > 0.6, "acceptance rate = {rate}");
    }

    #[test]
    fn invalid_gradients_and_options_are_rejected() {
        let wrong = |_: &[f64]| vec![1.0];
        assert!(hmc(&correlated_log, &wrong, &[0.0, 0.0], 10, &options()).is_err());
        let nan = |_: &[f64]| vec![f64::NAN, 0.0];
        assert!(hmc(&correlated_log, &nan, &[0.0, 0.0], 10, &options()).is_err());
        let mut bad = options();
        bad.step_size = 0.0;
        assert!(hmc(&correlated_log, &correlated_grad, &[0.0, 0.0], 10, &bad).is_err());
        let mut bad = options();
        bad.leapfrog_steps = 0;
        assert!(hmc(&correlated_log, &correlated_grad, &[0.0, 0.0], 10, &bad).is_err());
        assert!(hmc(&correlated_log, &correlated_grad, &[], 10, &options()).is_err());
    }
}
