//! MLMC estimates and variance retain finite payoff offsets.

use rssn::sim::stochastic::Method;
use rssn::sim::stochastic::MlmcOptions;
use rssn::sim::stochastic::mlmc;

fn level_zero() -> MlmcOptions {
    MlmcOptions {
        method: Method::EulerMaruyama,
        seed: 0,
        base_steps: 1,
        min_level: 0,
        max_level: 0,
        warmup_samples: 10,
        max_samples_per_level: 10,
        target_rmse: 1.0,
        ..MlmcOptions::default()
    }
}

#[test]
fn shifted_payoffs_keep_the_variance_of_the_observed_samples() {
    let mut payoffs = Vec::new();
    let result = mlmc(
        |_, _| Ok(vec![0.0]),
        |_, _| Ok(vec![1.0]),
        |state| {
            let value = (1_u64 << 42) as f64 + state[0];
            payoffs.push(value);
            Ok(value)
        },
        &[0.0],
        &level_zero(),
    )
    .unwrap();
    let centered: Vec<_> = payoffs.iter().map(|value| value - payoffs[0]).collect();
    let mean = centered.iter().sum::<f64>() / centered.len() as f64;
    let expected = centered
        .iter()
        .map(|value| (value - mean).powi(2))
        .sum::<f64>()
        / (centered.len() * (centered.len() - 1)) as f64;
    assert!(expected > 0.001);
    assert!(
        (result.variance / expected - 1.0).abs() < 1e-12,
        "reference {expected:.17e}, reported {:.17e}",
        result.variance
    );
    assert_eq!(result.estimate, payoffs[0] + mean);
}

#[test]
fn large_constant_payoff_keeps_a_finite_mean_and_zero_variance() {
    let result = mlmc(
        |_, _| Ok(vec![0.0]),
        |_, _| Ok(vec![0.0]),
        |_| Ok(1e308),
        &[0.0],
        &level_zero(),
    )
    .unwrap();
    assert_eq!(result.estimate, 1e308);
    assert_eq!(result.variance, 0.0);
    assert_eq!(result.std_error, 0.0);
}
