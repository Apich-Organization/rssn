//! Seeded MLMC streams are isolated by level while fine and coarse paths stay coupled.

use rssn::sim::stochastic::Method;
use rssn::sim::stochastic::MlmcOptions;
use rssn::sim::stochastic::MlmcResult;
use rssn::sim::stochastic::mlmc;

fn hierarchy() -> MlmcOptions {
    MlmcOptions {
        method: Method::EulerMaruyama,
        seed: 0,
        base_steps: 1,
        min_level: 3,
        max_level: 3,
        warmup_samples: 10,
        max_samples_per_level: 10,
        target_rmse: 1.0,
        ..MlmcOptions::default()
    }
}

fn stream_trace(workers: usize) -> (MlmcResult, Vec<f64>) {
    rayon::ThreadPoolBuilder::new()
        .num_threads(workers)
        .build()
        .unwrap()
        .install(|| {
            let mut payoffs = Vec::new();
            let result = mlmc(
                |_, _| Ok(vec![0.0]),
                |t, _| Ok(vec![if t == 0.0 { 1.0 } else { 0.0 }]),
                |state| {
                    payoffs.push(state[0]);
                    Ok(state[0])
                },
                &[0.0],
                &hierarchy(),
            )
            .unwrap();
            (result, payoffs)
        })
}

fn float_bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

#[test]
fn cross_level_fine_paths_do_not_reuse_a_stream() {
    let (stats, payoffs) = stream_trace(1);
    assert_eq!(
        stats
            .levels
            .iter()
            .map(|stats| stats.samples)
            .collect::<Vec<_>>(),
        [10; 4]
    );
    assert_eq!(payoffs.len(), 70);
    // Level zero contributes ten fine payoffs. Levels 1 and 2 then
    // contribute twenty fine/coarse payoffs each. At level 3, the
    // third fine payoff is at offset 54. The first step sizes differ
    // by a factor of 4, so multiply the level 3 payoff by 2.
    let a = payoffs[10];
    let b = 2.0 * payoffs[54];
    assert_ne!(a.to_bits(), b.to_bits());
}

#[test]
fn caller_pool_size_does_not_change_seeded_paths_or_statistics() {
    let (reference, payoffs) = stream_trace(1);
    for workers in [2, 4, 8] {
        let (other, stream_payoffs) = stream_trace(workers);
        assert_eq!(reference, other);
        assert_eq!(float_bits(&payoffs), float_bits(&stream_payoffs));
    }
}

#[test]
fn same_level_fine_and_coarse_retain_the_brownian_coupling() {
    let mut payoffs = Vec::new();
    let result = mlmc(
        |_, _| Ok(vec![0.0]),
        |_, _| Ok(vec![1.0]),
        |state| {
            payoffs.push(state[0]);
            Ok(state[0])
        },
        &[0.0],
        &hierarchy(),
    )
    .unwrap();
    assert_eq!(result.levels.len(), 4);
    for pair in payoffs[10..].as_chunks::<2>().0 {
        assert!((pair[0] - pair[1]).abs() <= 2e-15, "{pair:?}");
    }
}

#[test]
fn first_seeded_terminal_solve_matches_the_full_path_and_callback_order() {
    use std::cell::RefCell;

    use rssn::kernels::sde;
    use rssn::sim::stochastic::ensemble_seed;
    use rssn::sim::stochastic::wiener;

    let initial = [0.5];
    let drift = |t: f64, y: &[f64]| 0.125 * t - 0.5 * y[0];
    let diffusion = |_: f64, y: &[f64]| 0.25 + 0.125 * y[0];
    for method in [Method::EulerMaruyama, Method::Milstein] {
        let options = MlmcOptions {
            method,
            t0: 0.125,
            t1: 0.625,
            base_steps: 4,
            min_level: 0,
            max_level: 0,
            warmup_samples: 2,
            max_samples_per_level: 2,
            ..MlmcOptions::default()
        };
        let trace = RefCell::new(Vec::new());
        let mut first = None;
        mlmc(
            |t, y| {
                trace.borrow_mut().push(('a', t.to_bits(), y[0].to_bits()));
                Ok(vec![drift(t, y)])
            },
            |t, y| {
                trace.borrow_mut().push(('g', t.to_bits(), y[0].to_bits()));
                Ok(vec![diffusion(t, y)])
            },
            |y| {
                first.get_or_insert_with(|| (y[0].to_bits(), trace.borrow().clone()));
                Ok(y[0])
            },
            &initial,
            &options,
        )
        .unwrap();

        let dt = (options.t1 - options.t0) / options.base_steps as f64;
        let increments = wiener(ensemble_seed(options.seed, 0), 1, options.base_steps, dt);
        let reference = RefCell::new(Vec::new());
        let record_a = |t: f64, y: &[f64], out: &mut [f64]| {
            reference
                .borrow_mut()
                .push(('a', t.to_bits(), y[0].to_bits()));
            out[0] = drift(t, y);
        };
        let record_g = |t: f64, y: &[f64], out: &mut [f64]| {
            reference
                .borrow_mut()
                .push(('g', t.to_bits(), y[0].to_bits()));
            out[0] = diffusion(t, y);
        };
        let states = match method {
            | Method::EulerMaruyama => {
                sde::euler_maruyama(record_a, record_g, 1, options.t0, dt, &initial, &increments)
            },
            | Method::Milstein => {
                sde::milstein(record_a, record_g, 1, options.t0, dt, &initial, &increments)
            },
            | _ => unreachable!(),
        }
        .unwrap();
        assert_eq!(states.len(), options.base_steps + 1);
        assert_eq!(
            first.unwrap(),
            (states.last().unwrap()[0].to_bits(), reference.into_inner()),
            "{method:?}"
        );
    }
}
