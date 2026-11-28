//! Original-system verification and workspace boundaries for BiCGSTAB(l).

use std::cell::Cell;

use rssn::kernels::krylov::identity_precond;
use rssn::kernels::krylov_ext::bicgstab_l;

#[test]
fn exhausted_budget_does_not_certify_ordinary_improvement() {
    let calls = Cell::new(0);
    let apply = |x: &[f64], y: &mut [f64]| {
        calls.set(calls.get() + 1);
        y.copy_from_slice(x);
    };
    let initial = [0.5; 2];
    let result = bicgstab_l(
        apply,
        |_, _| unreachable!(),
        2,
        &[1.0; 2],
        Some(&initial),
        1e-12,
        0,
    );
    assert_eq!(calls.get(), 1);
    assert_eq!(result.x, initial);
    assert_eq!(result.iterations, 0);
    assert_eq!(result.residual, 0.5);
    assert!(!result.converged);

    let result = bicgstab_l(
        |x, y| y.copy_from_slice(&[x[0], 2.0 * x[1]]),
        identity_precond,
        1,
        &[1.0; 2],
        None,
        1e-12,
        1,
    );
    let exact = (1.0 - result.x[0]).hypot(1.0 - 2.0 * result.x[1]) / 2.0_f64.sqrt();
    assert_eq!(result.iterations, 1);
    assert!((result.residual - exact).abs() < 1e-16);
    assert!((1e-12..1.0).contains(&result.residual));
    assert!(!result.converged);
}

#[test]
fn actual_residuals_cover_extreme_scales_and_nonzero_subnormal_ratios() {
    for scale in [2.0_f64.powi(-600), 2.0_f64.powi(600)] {
        let result = bicgstab_l(
            identity_precond,
            identity_precond,
            2,
            &[scale; 2],
            None,
            1e-12,
            0,
        );
        assert_eq!(result.residual, 1.0);
        assert!(!result.converged);
    }
    let huge = 0.75 * f64::MAX;
    let result = bicgstab_l(
        identity_precond,
        identity_precond,
        2,
        &[huge],
        Some(&[-huge]),
        1e-12,
        0,
    );
    assert_eq!(result.residual, 2.0);
    assert!(!result.converged);

    let tiny = f64::from_bits(1);
    let result = bicgstab_l(
        identity_precond,
        identity_precond,
        1,
        &[huge, tiny],
        Some(&[huge, 0.0]),
        0.0,
        0,
    );
    assert_eq!(result.residual, tiny);
    assert!(!result.converged);
}

#[test]
fn right_preconditioned_nonsymmetric_solve_matches_manufactured_solution() {
    let expected = [1.0, -2.0, 3.0];
    let apply = |x: &[f64], out: &mut [f64]| {
        out.copy_from_slice(&[
            4.0 * x[0] - x[1] + 0.25 * x[2],
            2.0 * x[0] + 3.0 * x[1] - x[2],
            0.5 * x[1] + 2.5 * x[2],
        ]);
    };
    let mut b = [0.0; 3];
    apply(&expected, &mut b);
    for l in [0, 1, 2, 3, 4] {
        let result = bicgstab_l(
            apply,
            |r, z| z.copy_from_slice(&[r[0] / 4.0, r[1] / 3.0, r[2] / 2.5]),
            l,
            &b,
            Some(&[0.125, -0.25, 0.5]),
            1e-11,
            60,
        );
        let mut ax = [0.0; 3];
        apply(&result.x, &mut ax);
        let norm = |x: &[f64]| x.iter().fold(0.0_f64, |sum, &v| sum.hypot(v));
        let true_ratio = norm(&[b[0] - ax[0], b[1] - ax[1], b[2] - ax[2]]) / norm(&b);
        assert!((result.residual - true_ratio).abs() < 2e-16);
        assert_eq!(result.converged, true_ratio <= 1e-11, "l={l}");
        if l <= 3 {
            assert!(result.converged, "l={l}: {result:?}");
        }
        if result.converged {
            for (&actual, exact) in result.x.iter().zip(expected) {
                assert!((actual - exact).abs() < 2e-10, "l={l}");
            }
        }
    }
}

#[test]
fn empty_zero_rhs_and_invalid_controls_are_diagnosed() {
    let empty = bicgstab_l(
        |_, _| unreachable!(),
        |_, _| unreachable!(),
        4,
        &[],
        None,
        0.0,
        20,
    );
    assert!(empty.converged && empty.x.is_empty() && empty.residual == 0.0);

    let zero = bicgstab_l(
        identity_precond,
        identity_precond,
        2,
        &[0.0],
        Some(&[1.0]),
        1e-12,
        0,
    );
    assert!(!zero.converged && zero.residual.is_infinite());
    let solved = bicgstab_l(
        identity_precond,
        identity_precond,
        2,
        &[0.0],
        Some(&[1.0]),
        1e-12,
        5,
    );
    assert!(solved.converged && solved.residual == 0.0 && solved.x == [0.0]);

    for (initial, tolerance) in [
        (Some(&[1.0, 1.0][..]), 1e-12),
        (None, f64::NAN),
        (None, -1.0),
    ] {
        let result = bicgstab_l(
            |_, _| unreachable!(),
            |_, _| unreachable!(),
            2,
            &[1.0],
            initial,
            tolerance,
            10,
        );
        assert_eq!(result.x.len(), 1);
        assert!(!result.converged && result.residual.is_infinite() && result.iterations == 0);
    }
    let degree_overflow = bicgstab_l(
        identity_precond,
        identity_precond,
        usize::MAX,
        &[1.0],
        None,
        1e-12,
        10,
    );
    assert!(!degree_overflow.converged && degree_overflow.residual == 1.0);
}

#[test]
fn partial_operator_output_cannot_produce_a_convergence_claim() {
    let result = bicgstab_l(
        |x, y| y[0] = x[0],
        identity_precond,
        2,
        &[1.0; 2],
        None,
        1e-12,
        10,
    );
    assert!(result.residual.is_infinite() && !result.converged);
    let result = bicgstab_l(
        identity_precond,
        |x, y| y[0] = x[0],
        2,
        &[1.0; 2],
        None,
        1e-12,
        10,
    );
    assert_eq!(result.x, [0.0; 2]);
    assert_eq!(result.residual, 1.0);
    assert!(!result.converged);

    let calls = Cell::new(0);
    let result = bicgstab_l(
        |x, y| {
            calls.set(calls.get() + 1);
            if calls.get() == 2 {
                y[0] = x[0];
            } else {
                y.copy_from_slice(x);
            }
        },
        identity_precond,
        2,
        &[1.0; 2],
        None,
        1e-12,
        10,
    );
    assert_eq!(result.x, [0.0; 2]);
    assert_eq!(result.residual, 1.0);
    assert!(!result.converged);
}
