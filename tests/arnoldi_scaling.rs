//! Analytic factorisations and numerical boundaries for plain Arnoldi.

use std::cell::Cell;

use rssn::kernels::eigs::arnoldi;

fn close(
    actual: f64,
    expected: f64,
    tolerance: f64,
) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:e}, expected {expected:e}, tolerance {tolerance:e}"
    );
}

#[test]
fn start_and_operator_scales_preserve_the_diagonal_factorisation() {
    let unit = std::f64::consts::FRAC_1_SQRT_2;
    for start_scale in [f64::from_bits(1), 1e-200, 1.0, 1e200, f64::MAX] {
        for magnitude in [2.0_f64.powi(-600), 1.0, 2.0_f64.powi(600)] {
            for sign in [1.0, -1.0] {
                let scale = sign * magnitude;
                let apply = |x: &[f64], y: &mut [f64]| {
                    y[0] = scale * x[0];
                    y[1] = (2.0 * scale) * x[1];
                };
                let result = arnoldi(apply, &[start_scale, start_scale], 2);
                let (v, h, beta) = &result;
                assert_eq!(*beta, 0.0);
                assert_eq!(v[2], [0.0; 2]);
                close(v[0][0], unit, 2e-15);
                close(v[0][1], unit, 2e-15);
                close(v[1][0], -sign * unit, 2e-15);
                close(v[1][1], sign * unit, 2e-15);
                close(h.at(0, 0) / magnitude, sign * 1.5, 3e-15);
                close(h.at(1, 1) / magnitude, sign * 1.5, 3e-15);
                close(h.at(1, 0) / magnitude, 0.5, 3e-15);
                close(h.at(0, 1) / magnitude, 0.5, 3e-15);
                assert_eq!(result, arnoldi(apply, &[start_scale, start_scale], 2));
            }
        }
    }
}

#[test]
fn subnormal_operator_retains_a_resolvable_projected_coupling() {
    let scale = f64::from_bits(128);
    let (v, h, beta) = arnoldi(
        |x, y| {
            y[0] = scale * x[0];
            y[1] = (2.0 * scale) * x[1];
        },
        &[1.0, 1.0],
        2,
    );
    // Allow callback arithmetic and the stored Hessenberg values to each
    // round to the nearest representable subnormal.
    close(h.at(1, 0) / scale, 0.5, 2.0 / 128.0);
    assert!(v[1][0] < -0.7 && v[1][1] > 0.7);
    assert_eq!(beta, 0.0);
}

#[test]
fn shift_operator_has_an_independent_residual_direction() {
    for scale in [2.0_f64.powi(-600), 1.0, 2.0_f64.powi(600)] {
        let (v, h, beta) = arnoldi(
            |x, y| {
                y.copy_from_slice(&[0.0, scale * x[0], -2.0 * scale * x[1], 3.0 * scale * x[2]]);
            },
            &[1.0, 0.0, 0.0, 0.0],
            2,
        );
        assert_eq!(
            v,
            [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, -1.0, 0.0],
            ]
        );
        assert_eq!(h.at(0, 0), 0.0);
        assert_eq!(h.at(0, 1), 0.0);
        assert_eq!(h.at(1, 1), 0.0);
        assert_eq!(h.at(1, 0), scale);
        assert_eq!(beta, 2.0 * scale);
    }
}

#[test]
fn near_breakdown_depends_on_the_relative_residual() {
    for scale in [2.0_f64.powi(-600), 1.0, 2.0_f64.powi(600)] {
        for (gap, expected_calls) in [(4e-15, 1), (4e-13, 2)] {
            let calls = Cell::new(0);
            let (v, _, beta) = arnoldi(
                |x, y| {
                    calls.set(calls.get() + 1);
                    y[0] = scale * x[0];
                    y[1] = (1.0 + gap) * scale * x[1];
                },
                &[1.0, 1.0],
                2,
            );
            assert_eq!(calls.get(), expected_calls);
            assert_eq!(v[1] == [0.0; 2], expected_calls == 1);
            assert_eq!(beta, 0.0);
        }
    }
}

#[test]
fn invariant_subspaces_zero_start_and_excess_directions_are_padded() {
    let (v, h, beta) = arnoldi(|_, _| unreachable!(), &[0.0; 2], 3);
    assert_eq!(v, [[0.0; 2]; 4]);
    assert!(h.data.iter().all(|&x| x == 0.0));
    assert_eq!(beta, 0.0);

    let (v, h, beta) = arnoldi(|_, _| unreachable!(), &[], 3);
    assert_eq!(v, [Vec::<f64>::new(), Vec::new(), Vec::new(), Vec::new()]);
    assert_eq!((h.rows, h.cols), (3, 3));
    assert_eq!(beta, 0.0);
    let (v, h, beta) = arnoldi(|_, _| unreachable!(), &[0.0, -f64::MAX], 0);
    assert_eq!(v, [[0.0, -1.0]]);
    assert_eq!((h.rows, h.cols, beta), (0, 0, 0.0));

    for scale in [0.0, 2.0_f64.powi(-600), 1.0, 2.0_f64.powi(600)] {
        let calls = Cell::new(0);
        let (v, h, beta) = arnoldi(
            |x, y| {
                calls.set(calls.get() + 1);
                for (target, source) in y.iter_mut().zip(x) {
                    *target = scale * source;
                }
            },
            &[1.0, -2.0, 2.0],
            8,
        );
        assert_eq!(calls.get(), 1);
        assert!(v[1..].iter().flatten().all(|&x| x == 0.0));
        close(h.at(0, 0), scale, 3e-15 * scale);
        assert!(h.data[1..].iter().all(|&x| x == 0.0));
        assert_eq!(beta, 0.0);
    }

    let calls = Cell::new(0);
    let (v, _, beta) = arnoldi(
        |x, y| {
            calls.set(calls.get() + 1);
            y.copy_from_slice(&[x[0], 2.0 * x[1]]);
        },
        &[1.0, 1.0],
        8,
    );
    assert_eq!(calls.get(), 2);
    assert_eq!(v.len(), 9);
    assert!(v[2..].iter().flatten().all(|&x| x == 0.0));
    assert_eq!(beta, 0.0);
}

#[test]
fn bad_input_reports_failure_without_erasing_completed_columns() {
    for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let (v, h, beta) = arnoldi(|_, _| unreachable!(), &[invalid, 1.0], 2);
        assert_eq!(v, [[0.0; 2]; 3]);
        assert!(h.at(0, 0).is_nan() && beta.is_nan());

        let calls = Cell::new(0);
        let (v, h, beta) = arnoldi(
            |x, y| {
                calls.set(calls.get() + 1);
                if calls.get() == 1 {
                    y.copy_from_slice(&[0.0, x[0], 0.0]);
                } else {
                    y.copy_from_slice(&[0.0, invalid, 0.0]);
                }
            },
            &[1.0, 0.0, 0.0],
            3,
        );
        assert_eq!(calls.get(), 2);
        assert_eq!(h.at(1, 0), 1.0);
        assert!(h.at(1, 1).is_nan() && beta.is_nan());
        assert_eq!(v, [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0; 3], [0.0; 3]]);
    }

    let (_, h, beta) = arnoldi(
        |x, y| {
            // This constant-entry matrix has an eigenvalue above f64::MAX,
            // even though A x has finite entries for this start.
            let entry = 0.625 * f64::MAX;
            y.fill(entry * x[0] + entry * x[1]);
        },
        &[1.0, 1.0],
        2,
    );
    assert!(h.at(0, 0).is_nan() && beta.is_nan());

    let (_, h, beta) = arnoldi(
        |x, y| {
            let entry = f64::MAX * x[0];
            y.copy_from_slice(&[0.0, entry, entry]);
        },
        &[1.0, 0.0, 0.0],
        1,
    );
    assert!(h.at(0, 0).is_nan() && beta.is_nan());
}
