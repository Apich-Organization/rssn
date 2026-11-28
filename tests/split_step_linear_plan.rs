//! Tests for fixed-potential split-step propagation on periodic grids.

use num_complex::Complex64;
use rssn::sim::split_step::LinearSplitPlan;
use rssn::sim::split_step::SplitOrder;
use rssn::sim::split_step::SplitStep;
use std::f64::consts::TAU;

const ORDERS: [SplitOrder; 2] = [SplitOrder::Strang, SplitOrder::Yoshida4];

fn harmonic(
    n: &[usize],
    modes: &[isize],
) -> Vec<Complex64> {
    (0..n.iter().product())
        .map(|mut index| {
            let phase = n
                .iter()
                .zip(modes)
                .map(|(&extent, &mode)| {
                    let coordinate = index % extent;
                    index /= extent;
                    TAU * mode as f64 * coordinate as f64 / extent as f64
                })
                .sum();
            Complex64::from_polar(1.0, phase)
        })
        .collect()
}

fn max_error(
    actual: &[Complex64],
    expected: &[Complex64],
) -> f64 {
    assert_eq!(actual.len(), expected.len());
    assert!(actual.iter().all(|x| x.re.is_finite() && x.im.is_finite()));
    actual
        .iter()
        .zip(expected)
        .map(|(a, b)| (a - b).norm())
        .fold(0.0, f64::max)
}

#[test]
fn periodic_modes_follow_the_exact_dispersion_for_both_orders() {
    let cases: &[(&[usize], &[f64], &[isize])] = &[
        (&[1], &[2.5], &[0]),
        (&[7], &[3.7], &[-3]),
        (&[8], &[4.0], &[4]),
        (&[7, 6], &[1.4, 2.6], &[-1, 3]),
        (&[8, 5, 3], &[3.0, 2.0, 1.7], &[2, -2, 1]),
        (&[1, 5, 1], &[3.0, 2.1, 1.5], &[0, -2, 0]),
    ];
    let mass = 1.7;
    let constant = -0.35;
    let dt = 0.015;
    let steps = 17;
    for &(n, lengths, modes) in cases {
        let initial = harmonic(n, modes);
        let potential = vec![constant; initial.len()];
        let energy = modes
            .iter()
            .zip(lengths)
            .map(|(&m, &l)| (TAU * m as f64 / l).powi(2) / (2.0 * mass))
            .sum::<f64>()
            + constant;
        let phase = Complex64::from_polar(1.0, -energy * dt * steps as f64);
        let expected: Vec<_> = initial.iter().map(|&x| x * phase).collect();
        for order in ORDERS {
            let plan = LinearSplitPlan::new(n, lengths, mass, &potential, dt, order).unwrap();
            assert_eq!(plan.size(), initial.len());
            let actual = plan.evolve(&initial, steps).unwrap();
            let error = max_error(&actual, &expected);
            assert!(error < 1.5e-13, "{n:?}, {order:?}: {error}");
        }
    }
}

#[test]
fn fixed_nonuniform_potential_matches_legacy_and_reverses() {
    let solver = SplitStep::centered(&[16, 8], &[3.5, 6.9], 1.7).unwrap();
    let potential = |_t: f64, [x, y, _]: [f64; 3]| {
        -0.3 + 0.2 * (TAU * x / 3.5).cos() + 0.4 * (TAU * y / 6.9).sin()
    };
    let values: Vec<_> = (0..solver.grid.size())
        .map(|index| potential(0.0, solver.x(index)))
        .collect();
    let initial = solver.sample(|[x, y, _]| {
        Complex64::new(0.1 + (TAU * x / 3.5).cos(), 0.5 * (TAU * y / 6.9).sin())
    });
    let dt = 0.018;
    let steps = 31;
    let initial_norm = solver.norm2(&initial);
    for order in ORDERS {
        let plan = solver.prepare_linear(&values, dt, order).unwrap();
        let prepared = plan.evolve(&initial, steps).unwrap();
        let callback = solver.evolve(&initial, &potential, 0.0, 0.0, dt, steps, order);
        assert!(max_error(&prepared, &callback) < 1.0e-13);
        assert_eq!(prepared, plan.evolve(&initial, steps).unwrap());
        assert!((solver.norm2(&prepared) / initial_norm - 1.0).abs() < 5.0e-14);
        let inverse = solver.prepare_linear(&values, -dt, order).unwrap();
        assert!(max_error(&inverse.evolve(&prepared, steps).unwrap(), &initial) < 1.0e-13);
        // Evolving two portions with the reusable plan preserves the full trajectory.
        let partial = plan.evolve(&initial, 7).unwrap();
        let resumed = plan.evolve(&partial, steps - 7).unwrap();
        assert!(max_error(&prepared, &resumed) < 1.0e-13);
    }
}

#[test]
fn extreme_lengths_times_and_amplitudes_keep_finite_phase_and_norm() {
    for sign in [-1, 1] {
        let length = 2.0_f64.powi(550 * sign);
        let dt = 2.0_f64.powi(800 * sign);
        let mass = 2.0_f64.powi(-300 * sign);
        let magnitude = 2.0_f64.powi(-600 * sign);
        let normalized = harmonic(&[8], &[1]);
        let input: Vec<_> = normalized.iter().map(|&x| x * magnitude).collect();
        let step = SplitStep::centered(&[8], &[length], mass).unwrap();
        let expected_norm = 2.0_f64.powi(-650 * sign);
        assert!((step.norm2(&input) / expected_norm - 1.0).abs() < 3.0e-15);
        let potential = [0.125 / dt; 8];
        let expected_phase = Complex64::from_polar(1.0, -(0.5 * TAU * TAU + 0.125) * 3.0);
        let expected: Vec<_> = normalized.iter().map(|&x| x * expected_phase).collect();
        for order in ORDERS {
            let plan = step.prepare_linear(&potential, dt, order).unwrap();
            let actual = plan.evolve(&input, 3).unwrap();
            let scaled: Vec<_> = actual.iter().map(|x| x / magnitude).collect();
            assert!(max_error(&scaled, &expected) < 6.0e-14);
            assert!((step.norm2(&actual) / expected_norm - 1.0).abs() < 1.0e-14);
        }
    }

    // The integrated value is still finite when the cell volume itself overflows.
    let length = 2.0_f64.powi(550);
    let initial = [Complex64::new(2.0_f64.powi(-600), 0.0); 16];
    let solver = SplitStep::centered(&[4, 4], &[length; 2], 1.0).unwrap();
    assert_eq!(solver.norm2(&initial), 2.0_f64.powi(-100));
    assert_eq!(solver.norm2(&[Complex64::ZERO; 16]), 0.0);
}

#[test]
fn exact_identity_and_invalid_data() {
    let [order, _] = ORDERS;
    let enormous = Complex64::new(f64::MAX, -f64::MIN_POSITIVE);
    let input = [enormous, Complex64::new(-0.0, 0.0)];
    let zero = LinearSplitPlan::new(&[2], &[1.0], 1.0, &[f64::MAX; 2], -0.0, order).unwrap();
    assert_eq!(
        zero.evolve(&input, usize::MAX).unwrap()[1].re.to_bits(),
        (-0.0_f64).to_bits()
    );
    assert_eq!(zero.evolve(&input, usize::MAX).unwrap(), input);
    let plan = LinearSplitPlan::new(&[2], &[1.0], 1.0, &[0.0; 2], 1.0, order).unwrap();
    assert_eq!(plan.evolve(&input, 0).unwrap(), input);
    assert_eq!(
        plan.evolve(&[Complex64::ZERO; 2], 2).unwrap(),
        vec![Complex64::ZERO; 2]
    );
    assert!(plan.evolve(&[], 1).is_err());
    assert!(plan.evolve(&[Complex64::new(0.0, f64::NAN); 2], 0).is_err());
    for (n, lengths, mass) in [
        (vec![], vec![], 1.0),
        (vec![0], vec![1.0], 1.0),
        (vec![2], vec![], 1.0),
        (vec![2], vec![0.0], 1.0),
        (vec![2], vec![f64::INFINITY], 1.0),
        (vec![2], vec![1.0], 0.0),
        (vec![2], vec![1.0], f64::NAN),
        (vec![usize::MAX, 2], vec![1.0; 2], 1.0),
    ] {
        assert!(LinearSplitPlan::new(&n, &lengths, mass, &[0.0; 2], 1.0, order).is_err());
    }
    for (values, dt) in [
        (vec![0.0], 1.0),
        (vec![0.0, f64::NAN], 1.0),
        (vec![0.0; 2], f64::NAN),
    ] {
        assert!(LinearSplitPlan::new(&[2], &[1.0], 1.0, &values, dt, order).is_err());
    }
    assert!(LinearSplitPlan::new(&[2], &[1.0], 1.0, &[f64::MAX; 2], 4.0, order).is_err());
    assert!(LinearSplitPlan::new(&[2], &[f64::MIN_POSITIVE], 1.0, &[0.0; 2], 1.0, order).is_err());

    let quarter_period = 0.5 / TAU;
    let combine =
        LinearSplitPlan::new(&[2], &[1.0], 1.0, &[0.0; 2], quarter_period, order).unwrap();
    let limit = f64::MAX;
    assert!(
        combine
            .evolve(
                &[Complex64::new(limit, limit), Complex64::new(limit, -limit)],
                1
            )
            .is_err()
    );
}
