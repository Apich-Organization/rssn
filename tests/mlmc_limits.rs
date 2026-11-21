//! Sample budgets and cooperative cancellation in the MLMC solver.

use std::sync::Arc;
use std::sync::atomic::AtomicBool;
use std::sync::atomic::Ordering;

use rssn::cancellation;
use rssn::sim::stochastic::Method;
use rssn::sim::stochastic::MlmcOptions;
use rssn::sim::stochastic::SdeError;
use rssn::sim::stochastic::mlmc;

fn short_hierarchy() -> MlmcOptions {
    MlmcOptions {
        method: Method::EulerMaruyama,
        seed: 42,
        base_steps: 1,
        min_level: 0,
        max_level: 0,
        warmup_samples: 2,
        max_samples_per_level: 3,
        target_rmse: 1e-10,
        ..MlmcOptions::default()
    }
}

#[test]
fn a_short_warmup_is_preserved_and_adaptive_sampling_respects_the_cap() {
    let options = short_hierarchy();
    for (noise, expected) in [(0.0, 2), (1.0, 3)] {
        let mut calls = 0;
        let result = mlmc(
            |_, _| Ok(vec![0.0]),
            |_, _| Ok(vec![noise]),
            |state| {
                calls += 1;
                Ok(state[0])
            },
            &[0.0],
            &options,
        )
        .unwrap();
        assert_eq!(result.levels[0].samples, expected);
        assert_eq!(calls, expected);
        assert_eq!(result.total_cost, expected);
    }
}

#[test]
fn impossible_sample_bounds_fail_before_callbacks() {
    for (warmup_samples, max_samples_per_level) in [(1, 3), (3, 2)] {
        let options = MlmcOptions {
            warmup_samples,
            max_samples_per_level,
            ..short_hierarchy()
        };
        let result = mlmc(
            |_, _| panic!("invalid sample bounds ran a drift callback"),
            |_, _| panic!("invalid sample bounds ran a diffusion callback"),
            |_| panic!("invalid sample bounds ran a payoff"),
            &[0.0],
            &options,
        );
        assert!(matches!(result, Err(SdeError::Invalid(_))));
    }
}

#[test]
fn cancellation_stops_without_returning_a_partial_estimate() {
    let options = MlmcOptions {
        min_level: 1,
        max_level: 1,
        ..short_hierarchy()
    };
    // Three calls would finish level 0 and evaluate the first fine payoff
    // at level 1. Cancellation there must skip its paired coarse payoff.
    for cancel_after in [0, 1, 3] {
        let token = Arc::new(AtomicBool::new(cancel_after == 0));
        let _scope = cancellation::install(Arc::clone(&token));
        let mut calls = 0;
        let result = mlmc(
            |_, _| Ok(vec![0.0]),
            |_, _| Ok(vec![0.0]),
            |state| {
                calls += 1;
                if calls == cancel_after {
                    token.store(true, Ordering::Relaxed);
                }
                Ok(state[0])
            },
            &[0.0],
            &options,
        );
        assert_eq!(calls, cancel_after);
        assert!(matches!(result, Err(SdeError::StepLimit(reason)) if reason.contains("cancelled")));
    }
}
