//! Scale, callback-order and symplectic checks for separable Hamiltonian steps.

use std::cell::RefCell;

use rssn::sim::integrators::SymplecticIntegratorMethod;
use rssn::sim::integrators::symplectic_step_with_workspace;
use rssn::sim::integrators::verlet_step;
use rssn::sim::integrators::yoshida4_step;

#[test]
fn constant_acceleration_retains_finite_kicks_and_drifts_at_extreme_scales() {
    for (scale, inverse) in [
        (2.0_f64.powi(700), 2.0_f64.powi(-700)),
        (2.0_f64.powi(-700), 2.0_f64.powi(700)),
    ] {
        for dt in [scale, -scale] {
            let acceleration = |_: &[f64]| vec![inverse];
            let (mut q, mut p) = ([0.0], [0.0]);
            verlet_step(&mut q, &mut p, scale, dt, &acceleration);
            // Constant acceleration is integrated exactly by Verlet. In this case h = ±m
            // and a = 1/m, so p = h and the signed drift cancels the sign of p/2.
            assert_eq!(q, [0.5 * scale]);
            assert_eq!(p, [dt]);

            let (mut q, mut p) = ([0.0], [0.0]);
            yoshida4_step(&mut q, &mut p, scale, dt, &acceleration);
            assert!((q[0] / scale - 0.5).abs() < 16.0 * f64::EPSILON);
            assert!((p[0] / dt - 1.0).abs() < 16.0 * f64::EPSILON);
        }
    }
    // Halving the least subnormal itself rounds to zero. The completed kick is finite.
    let (mut q, mut p) = ([0.0], [0.0]);
    verlet_step(
        &mut q,
        &mut p,
        2.0_f64.powi(1023),
        f64::from_bits(1),
        &|_| vec![2.0_f64.powi(51)],
    );
    assert_eq!(p, [1.0]);
    assert_eq!(q, [0.0]);
}

#[test]
fn free_drift_ignores_an_unrepresentable_intermediate_product() {
    for magnitude in [2.0_f64.powi(-700), 2.0_f64.powi(700)] {
        let (mut q, mut p) = ([0.0, 0.0], [magnitude, -magnitude]);
        verlet_step(&mut q, &mut p, magnitude, magnitude, &|_| vec![0.0; 2]);
        assert_eq!(q, [magnitude, -magnitude]);
        assert_eq!(p, [magnitude, -magnitude]);
    }

    let large = 2.0_f64.powi(700);
    let (mut q, mut p) = ([8.0, 0.0, 4.0], [1.0, large, 0.0]);
    verlet_step(&mut q, &mut p, large, large, &|_| vec![0.0; 3]);
    assert_eq!(q, [9.0, large, 4.0]);
    assert_eq!(p, [1.0, large, 0.0]);
}

#[test]
fn buffered_step_preserves_raw_callback_arguments_and_order() {
    for (method, calls) in [
        (SymplecticIntegratorMethod::Verlet, 2),
        (SymplecticIntegratorMethod::Yoshida4, 6),
    ] {
        let mass = 1.5;
        let h = -0.0625;
        let (mut q, mut p) = ([0.25, -0.5], [0.5, 0.75]);
        let (mut raw_q, mut raw_p) = (q, p);
        let (mut reference_q, mut reference_p) = (q, p);
        let mut reference_calls = Vec::new();
        let call = |positions: [f64; 2], call_count: usize| {
            std::array::from_fn::<_, 2, _>(|i| -positions[i] + f64::from(call_count as u32) * 0.125)
        };
        let mut reference_kick = |q, p: &mut [f64; 2], interval| {
            let a = call(q, reference_calls.len());
            reference_calls.push(q);
            for i in 0..2 {
                p[i] += 0.5 * interval * mass * a[i];
            }
        };
        let cube_root = 2.0_f64.cbrt();
        let outside = 1.0 / (2.0 - cube_root);
        let intervals = if method == SymplecticIntegratorMethod::Verlet {
            vec![h]
        } else {
            vec![outside * h, (-cube_root * outside) * h, outside * h]
        };
        for interval in intervals {
            reference_kick(reference_q, &mut reference_p, interval);
            for i in 0..2 {
                reference_q[i] += interval * reference_p[i] / mass;
            }
            reference_kick(reference_q, &mut reference_p, interval);
        }

        let raw_calls = RefCell::new(Vec::new());
        let raw = |q: &[f64]| {
            let q = <[f64; 2]>::try_from(q).unwrap();
            let mut visited = raw_calls.borrow_mut();
            let a = call(q, visited.len());
            visited.push(q);
            a.to_vec()
        };
        match method {
            | SymplecticIntegratorMethod::Verlet => {
                verlet_step(&mut raw_q, &mut raw_p, mass, h, &raw)
            },
            | SymplecticIntegratorMethod::Yoshida4 => {
                yoshida4_step(&mut raw_q, &mut raw_p, mass, h, &raw)
            },
        }
        let mut workspace = [0.0; 2];
        let workspace_ptr = workspace.as_ptr();
        let mut seen = Vec::new();
        symplectic_step_with_workspace(&mut q, &mut p, mass, h, method, &mut workspace, |q, a| {
            assert_eq!(a.as_ptr(), workspace_ptr);
            let q = <[f64; 2]>::try_from(q).unwrap();
            a.copy_from_slice(&call(q, seen.len()));
            seen.push(q);
        })
        .unwrap();
        assert_eq!(seen.len(), calls);
        assert_eq!(seen, reference_calls);
        assert_eq!(raw_calls.into_inner(), seen);
        assert_eq!((raw_q, raw_p), (reference_q, reference_p));
        assert_eq!((q, p), (reference_q, reference_p));
    }
}

#[test]
fn invalid_workspace_calls_leave_state_and_buffer_unchanged() {
    let (mut q, mut p) = ([1.0, -1.0], [0.0, -0.25]);
    let initial = (q, p);
    let mut workspace = [2.0, 3.0];
    let method = SymplecticIntegratorMethod::Verlet;
    let no_call = |_: &[f64], _: &mut [f64]| panic!("invalid calls must fail before the callback");
    for (mass, h) in [(0.0, 0.5), (f64::INFINITY, 0.5), (1.0, f64::NAN)] {
        assert!(
            symplectic_step_with_workspace(
                &mut q,
                &mut p,
                mass,
                h,
                method,
                &mut workspace,
                no_call
            )
            .is_err()
        );
    }
    assert!(
        symplectic_step_with_workspace(
            &mut q,
            &mut p[..1],
            1.0,
            0.5,
            method,
            &mut workspace,
            no_call
        )
        .is_err()
    );
    assert!(
        symplectic_step_with_workspace(
            &mut q,
            &mut p,
            1.0,
            0.5,
            method,
            &mut workspace[..1],
            no_call
        )
        .is_err()
    );
    assert!(
        symplectic_step_with_workspace(
            &mut q,
            &mut p,
            1.0,
            f64::MAX,
            SymplecticIntegratorMethod::Yoshida4,
            &mut workspace,
            no_call,
        )
        .is_err()
    );
    assert_eq!((q, p), initial);
    assert_eq!(workspace, [2.0, 3.0]);
}

#[test]
fn workspace_propagates_values_written_after_preflight() {
    let (mut q, mut p) = ([2.0], [1.0]);
    let mut work = [0.0];
    let mut calls = 0;
    symplectic_step_with_workspace(
        &mut q,
        &mut p,
        1.0,
        0.5,
        SymplecticIntegratorMethod::Verlet,
        &mut work,
        |_, a| {
            calls += 1;
            a[0] = if calls == 1 { 0.0 } else { f64::NAN };
        },
    )
    .unwrap();
    assert_eq!(calls, 2);
    assert_eq!(q, [2.5]);
    assert!(p[0].is_nan());
    assert!(work[0].is_nan());
}

#[test]
fn buffered_harmonic_step_is_reversible_with_bounded_energy() {
    for (method, tolerance) in [
        (SymplecticIntegratorMethod::Verlet, 2e-3),
        (SymplecticIntegratorMethod::Yoshida4, 2e-6),
    ] {
        let (mut q, mut p) = ([1.0], [0.0]);
        let mut workspace = [0.0];
        for h in [0.0625, -0.0625] {
            for _ in 0..1200 {
                symplectic_step_with_workspace(
                    &mut q,
                    &mut p,
                    2.0,
                    h,
                    method,
                    &mut workspace,
                    |q, a| {
                        a[0] = -q[0];
                    },
                )
                .unwrap();
                let energy = q[0] * q[0] + 0.25 * p[0] * p[0];
                assert!((energy - 1.0).abs() < tolerance);
            }
        }
        assert!((q[0] - 1.0).abs() < 5e-13);
        assert!(p[0].abs() < 5e-13);
    }
}
