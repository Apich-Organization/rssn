//! Sampling, signed steps and workspace behavior of separable compositions.

use rssn::sim::symplectic::Scheme;
use rssn::sim::symplectic::integrate;
use rssn::sim::symplectic::integrate_sampled;
use rssn::sim::symplectic::step;

fn harmonic_force(
    q: &[f64],
    out: &mut [f64],
) {
    for (out, position) in out.iter_mut().zip(q) {
        *out = -position;
    }
}

#[test]
fn sampled_path_keeps_the_initial_cadence_and_off_cadence_final_state() {
    for h in [-0.1, 0.1] {
        let (mut q, mut p) = ([1.0, 2.0], [0.2, -0.5]);
        let sampled =
            integrate_sampled(Scheme::Yoshida6, &q, &p, h, 7, 3, &mut harmonic_force).unwrap();
        assert_eq!(
            sampled.iter().map(|sample| sample.step).collect::<Vec<_>>(),
            [0, 3, 6, 7]
        );
        let mut completed = 0;
        for sample in sampled {
            integrate(
                Scheme::Yoshida6,
                &mut q,
                &mut p,
                h,
                sample.step - completed,
                &mut harmonic_force,
            );
            assert_eq!(sample.q, q);
            assert_eq!(sample.p, p);
            completed = sample.step;
        }
    }
}

#[test]
fn sampled_verlet_has_bounded_long_term_energy_and_reverses() {
    let steps = 40_000;
    let h = 0.05;
    let path = integrate_sampled(
        Scheme::Verlet,
        &[1.0],
        &[0.0],
        h,
        steps,
        137,
        &mut harmonic_force,
    )
    .unwrap();
    assert_eq!(path.last().unwrap().step, steps);
    for state in &path {
        let energy = 0.5 * (state.q[0].powi(2) + state.p[0].powi(2));
        assert!((energy - 0.5).abs() < 3.2e-4);
    }
    let final_state = path.last().unwrap();
    let backward = integrate_sampled(
        Scheme::Verlet,
        &final_state.q,
        &final_state.p,
        -h,
        steps,
        steps,
        &mut harmonic_force,
    )
    .unwrap();
    let back = backward.last().unwrap();
    assert!((back.q[0] - 1.0).abs() < 1e-11);
    assert!(back.p[0].abs() < 1e-11);
}

#[test]
fn raw_integrate_restarts_the_acceleration_buffer_between_steps() {
    let mut trace = Vec::new();
    let stages = Scheme::Yoshida4.stages();
    let (mut q, mut p) = ([1.0], [0.25, 3.0]);
    integrate(Scheme::Yoshida4, &mut q, &mut p, 0.02, 2, &mut |_, a| {
        trace.push(a[0]);
        a[0] += 1.0;
    });
    let (mut qr, mut pr) = ([1.0], [0.25, 3.0]);
    for _ in 0..2 {
        step(&stages, &mut qr, &mut pr, 0.02, &mut |_, a| a[0] += 1.0);
    }
    assert_eq!(trace, [0.0, 1.0, 2.0, 0.0, 1.0, 2.0]);
    assert_eq!((q, p), (qr, pr));
    assert_eq!(p[1], 3.0);
}

#[test]
fn checked_sampling_rejects_invalid_inputs_before_evaluating_the_force() {
    let mut never = |_: &[f64], _: &mut [f64]| panic!("unexpected force evaluation");
    for (q, p, h, steps, cadence) in [
        (vec![], vec![], 0.1, 1, 1),
        (vec![1.0], vec![2.0, 3.0], 0.1, 1, 1),
        (vec![1.0], vec![2.0], f64::INFINITY, 1, 1),
        (vec![1.0], vec![f64::NAN], 0.1, 1, 1),
        (vec![1.0], vec![2.0], 0.1, 1, 0),
        (vec![1.0], vec![2.0], 0.1, usize::MAX, 1),
    ] {
        assert!(integrate_sampled(Scheme::Verlet, &q, &p, h, steps, cadence, &mut never).is_err());
    }
    assert!(integrate_sampled(Scheme::Verlet, &[1.0], &[4.0], f64::MAX, 1, 1, &mut never).is_err());
    let stationary =
        integrate_sampled(Scheme::Verlet, &[1.0], &[2.0], -0.0, 7, 3, &mut never).unwrap();
    assert_eq!(
        stationary.iter().map(|s| s.step).collect::<Vec<_>>(),
        [0, 3, 6, 7]
    );
    assert!(stationary.iter().all(|s| s.q == [1.0] && s.p == [2.0]));
    let initial = integrate_sampled(Scheme::Verlet, &[1.0], &[2.0], 0.1, 0, 3, &mut never).unwrap();
    assert_eq!(initial, [stationary[0].clone()]);
}

#[test]
fn incomplete_force_and_cancellation_on_the_last_step_return_errors() {
    let mut incomplete = |_: &[f64], out: &mut [f64]| out[0] = 0.0;
    assert!(
        integrate_sampled(
            Scheme::Verlet,
            &[1.0, 0.0],
            &[0.0, 1.0],
            0.1,
            2,
            1,
            &mut incomplete
        )
        .is_err()
    );

    use std::sync::Arc;
    use std::sync::atomic::AtomicBool;
    use std::sync::atomic::Ordering;
    let cancelled = Arc::new(AtomicBool::new(false));
    let _scope = rssn::cancellation::install(Arc::clone(&cancelled));
    let mut calls = 0;
    let mut force = |_: &[f64], out: &mut [f64]| {
        calls += 1;
        out.fill(0.0);
        cancelled.store(true, Ordering::Relaxed);
    };
    let error =
        integrate_sampled(Scheme::Yoshida8, &[1.0], &[0.0], 0.1, 1, 1, &mut force).unwrap_err();
    assert!(error.contains("cancelled"));
    assert_eq!(calls, 1);
}
