//! Directed-loop storage, block statistics and two-spin spectrum regressions.

use rssn::sim::sse::SseResult;
use rssn::sim::sse::SseXxz;
use rssn::sim::sse::chain_bonds;

fn bits(r: SseResult) -> [u64; 6] {
    [
        r.energy,
        r.energy_error,
        r.specific_heat,
        r.magnetisation_sq,
        r.staggered_sq,
        r.mean_ops,
    ]
    .map(f64::to_bits)
}

fn materialized(
    simulation: &mut SseXxz,
    therm: usize,
    sweeps: usize,
) -> SseResult {
    for _ in 0..therm {
        simulation.sweep();
    }
    let samples: Vec<_> = (0..sweeps)
        .map(|_| {
            simulation.sweep();
            simulation.sample()
        })
        .collect();
    let count = sweeps.max(1) as f64;
    let orders: Vec<_> = samples.iter().map(|s| s.n_ops as f64).collect();
    let n = orders.iter().fold(0.0, |sum, &x| sum + x) / count;
    let n2 = orders.iter().fold(0.0, |sum, &x| sum + x * x) / count;
    let blocks = sweeps.clamp(1, 20);
    let length = sweeps / blocks;
    let complete: Vec<_> = (0..blocks)
        .map(|b| orders[b * length..(b + 1) * length].iter().sum::<f64>() / length.max(1) as f64)
        .collect();
    let means = complete.iter().sum::<f64>() / blocks as f64;
    let variance =
        complete.iter().map(|x| (x - means).powi(2)).sum::<f64>() / (blocks.max(2) - 1) as f64;
    SseResult {
        energy: simulation.n_bonds() as f64 * simulation.constant() - n / simulation.beta,
        energy_error: (variance / blocks as f64).sqrt() / simulation.beta,
        specific_heat: n2 - n * n - n,
        magnetisation_sq: samples
            .iter()
            .fold(0.0, |sum, s| sum + s.magnetisation.powi(2))
            / count,
        staggered_sq: samples.iter().fold(0.0, |sum, s| sum + s.staggered.powi(2)) / count,
        mean_ops: n,
    }
}

#[test]
fn streamed_blocks_match_the_materialized_estimator_across_split_runs() {
    let init = SseXxz::new(8, chain_bonds(8, true), 1.25, 2.3, 2025).unwrap();
    let (mut fit, mut reference) = (init.clone(), init);
    for (i, count) in [0, 1, 19, 20, 21, 39, 40, 41, 103, 400]
        .into_iter()
        .enumerate()
    {
        // A new admissible anisotropy and temperature affect the next sweep.
        fit.jz = [1.25, 0.25, -0.75, -0.0, 0.0, 1.375][i % 6];
        fit.beta = 0.75 + i as f64 * 0.3;
        reference.jz = fit.jz;
        reference.beta = fit.beta;
        assert_eq!(
            bits(fit.run(6, count)),
            bits(materialized(&mut reference, 6, count))
        );
        fit.sweep();
        reference.sweep();
        let (a, b) = (fit.sample(), reference.sample());
        assert_eq!(a.n_ops, b.n_ops);
        assert_eq!(a.magnetisation.to_bits(), b.magnetisation.to_bits());
        assert_eq!(a.staggered.to_bits(), b.staggered.to_bits());
    }
}

#[test]
fn low_temperature_dimer_energy_and_spin_sectors_follow_the_closed_form_spectrum() {
    let jz = 0.8;
    let beta = 10.0;
    let levels = [jz * 0.25, jz * 0.25, -jz * 0.25 + 0.5, -jz * 0.25 - 0.5];
    let weights = levels.map(|e: f64| (-beta * e).exp());
    let z = weights.iter().sum::<f64>();
    let mean_e = levels.iter().zip(weights).map(|(e, w)| e * w).sum::<f64>() / z;
    let parallel = (weights[0] + weights[1]) / z;
    let mut sim = SseXxz::new(2, vec![(0, 1)], jz, beta, 29).unwrap();
    let result = sim.run(1_000, 45_000);
    assert!((result.energy - mean_e).abs() < 5.0 * result.energy_error + 0.006);
    assert!((result.magnetisation_sq - parallel).abs() < 0.001);
    assert!((result.staggered_sq + result.magnetisation_sq - 1.0).abs() < 2e-16);
}

#[test]
fn zero_measurement_sweeps_keeps_observables_empty() {
    let mut sim = SseXxz::new(4, chain_bonds(4, true), -0.3, 1.2, 13).unwrap();
    let result = sim.run(10, 0);
    let expected_shift = sim.n_bonds() as f64 * sim.constant();
    assert_eq!(bits(result), [expected_shift.to_bits(), 0, 0, 0, 0, 0]);
}

#[test]
fn constructor_rejects_nonfinite_settings_before_using_the_rng() {
    for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(SseXxz::new(2, vec![(0, 1)], invalid, 1.0, 9).is_err());
        assert!(SseXxz::new(2, vec![(0, 1)], 1.0, invalid, 9).is_err());
    }
}
