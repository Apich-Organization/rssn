//! Scale and failure contracts for matrix-free iterative solvers.
use rssn::kernels::krylov::IterResult;
use rssn::kernels::krylov::{
    self,
};
enum Method {
    Cg,
    Gmres(usize),
    Bicgstab,
}
fn solve(
    apply: impl Fn(&[f64], &mut [f64]),
    precond: impl Fn(&[f64], &mut [f64]),
    b: &[f64],
    x0: Option<&[f64]>,
    tol: f64,
    max_iter: usize,
    method: Method,
) -> IterResult {
    match method {
        | Method::Cg => krylov::cg(apply, precond, b, x0, tol, max_iter),
        | Method::Gmres(restart) => krylov::gmres(apply, precond, b, x0, restart, tol, max_iter),
        | Method::Bicgstab => krylov::bicgstab(apply, precond, b, x0, tol, max_iter),
    }
}
#[test]
fn diagonal_scale_sweep() {
    for method in [0, 1, 2] {
        for rhs in [1e-300, 1.0, 1e300] {
            for a in [1e-300, 1.0, 1e300] {
                for p in [1e-300, 1.0, 1e300] {
                    let expected: f64 = rhs / a;
                    if !expected.is_finite() || expected == 0.0 {
                        continue;
                    }
                    let mode = match method {
                        | 0 => Method::Cg,
                        | 1 => Method::Gmres(2),
                        | _ => Method::Bicgstab,
                    };
                    let result = solve(
                        |x, y| {
                            y[0] = a * x[0];
                            y[1] = a * x[1];
                        },
                        |r, z| {
                            z[0] = p * r[0];
                            z[1] = p * r[1];
                        },
                        &[rhs, -rhs],
                        None,
                        1e-12,
                        10,
                        mode,
                    );
                    assert!(
                        result.converged,
                        "method={method} b={rhs:e} a={a:e} p={p:e} {result:?}"
                    );
                    assert!((result.x[0] / expected - 1.0).abs() < 1e-12, "{result:?}");
                }
            }
        }
    }
}
fn mode(m: usize) -> Method {
    match m {
        | 0 => Method::Cg,
        | 1 => Method::Gmres(3),
        | _ => Method::Bicgstab,
    }
}
#[test]
fn coupled_matrix_scales_and_initial_guesses() {
    for method in 0..3 {
        for scale in [1e-300, 1.0, 1e300] {
            for guess in [None, Some([0.2 * scale, -0.4 * scale])] {
                let b = [4.0 * scale, 7.0 * scale];
                let result = solve(
                    |x, y| {
                        y[0] = 4.0 * x[0] + x[1];
                        y[1] = x[0] + 3.0 * x[1];
                    },
                    |r, z| {
                        z[0] = r[0] / 4.0;
                        z[1] = r[1] / 3.0;
                    },
                    &b,
                    guess.as_ref().map(|v| v.as_slice()),
                    1e-12,
                    20,
                    mode(method),
                );
                assert!(
                    result.converged,
                    "method{method} scale{scale:e} guess{guess:?} {result:?}"
                );
                assert!(
                    (result.x[0] / scale - 5.0 / 11.0).abs() < 1e-11,
                    "{result:?}"
                );
                assert!(
                    (result.x[1] / scale - 24.0 / 11.0).abs() < 1e-11,
                    "{result:?}"
                );
            }
        }
    }
}
#[test]
fn invalid_inputs_and_incomplete_operator_never_converge() {
    for method in 0..3 {
        for (b, guess, tol) in [
            (vec![1.0], Some(vec![]), 1e-8),
            (vec![f64::NAN], None, 1e-8),
            (vec![1.0], None, -1.0),
            (vec![1.0], None, f64::NAN),
        ] {
            let result = solve(
                |_, _| panic!("invalid input reached operator"),
                |_, _| panic!("invalid input reached preconditioner"),
                &b,
                guess.as_deref(),
                tol,
                10,
                mode(method),
            );
            assert!(!result.converged);
        }
        let result = solve(
            |_, y| {
                y[0] = 1.0;
            },
            |r, z| z.copy_from_slice(r),
            &[1.0, 2.0],
            None,
            1e-8,
            10,
            mode(method),
        );
        assert!(!result.converged, "{result:?}");
        let result = solve(
            |x, y| y.copy_from_slice(x),
            |_, _| {},
            &[1.0, 2.0],
            None,
            1e-8,
            10,
            mode(method),
        );
        assert!(!result.converged, "{result:?}");
    }
}
#[test]
fn zero_rhs_and_budgets_are_explicit() {
    for method in 0..3 {
        let result = solve(
            |_, _| panic!("A0 already known"),
            |_, _| panic!("already solved"),
            &[0.0, 0.0],
            None,
            0.0,
            0,
            mode(method),
        );
        assert!(result.converged);
        assert_eq!(result.iterations, 0);
        let result = solve(
            |x, y| y.copy_from_slice(x),
            |r, z| z.copy_from_slice(r),
            &[0.0, 0.0],
            Some(&[1.0, 2.0]),
            0.0,
            10,
            mode(method),
        );
        assert!(result.converged, "{result:?}");
        assert_eq!(result.residual, 0.0);
        let result = solve(
            |x, y| y.copy_from_slice(x),
            |r, z| z.copy_from_slice(r),
            &[1e-300, 2e-300],
            None,
            1e-12,
            0,
            mode(method),
        );
        assert!(!result.converged);
        assert_eq!(result.iterations, 0);
        assert_eq!(result.residual, 1.0);
    }
}

#[test]
fn solved_large_components_do_not_hide_unsolved_tiny_ones() {
    for method in 0..3 {
        let b = [1e300, 1e-300];
        let initial = [1e300, 0.0];
        let stopped = solve(
            |x, y| y.copy_from_slice(x),
            |r, z| z.copy_from_slice(r),
            &b,
            Some(&initial),
            0.0,
            0,
            mode(method),
        );
        assert!(!stopped.converged, "{stopped:?}");
        assert!(stopped.residual > 0.0);
        let solved = solve(
            |x, y| y.copy_from_slice(x),
            |r, z| z.copy_from_slice(r),
            &b,
            Some(&initial),
            0.0,
            4,
            mode(method),
        );
        assert!(solved.converged, "{solved:?}");
        assert_eq!(solved.x, b);
        assert_eq!(solved.residual, 0.0);
    }
}

#[test]
fn residual_refinement_recovers_components_lost_by_the_first_scale() {
    for method in 0..3 {
        for b in [[1e300, 1e-300], [0.0, 0.0]] {
            let initial = if b[0] == 0.0 {
                Some([1e300, 1e-300])
            } else {
                None
            };
            let result = solve(
                |x, y| y.copy_from_slice(x),
                |r, z| z.copy_from_slice(r),
                &b,
                initial.as_ref().map(|x| x.as_slice()),
                0.0,
                4,
                mode(method),
            );
            assert!(result.converged, "{result:?}");
            assert_eq!(result.x, b);
            assert!(result.iterations <= 4);
        }
    }
}

#[test]
fn uniformly_subnormal_systems_retain_finite_solutions() {
    for method in 0..3 {
        for scale in [1e-320, f64::from_bits(1), 1e308] {
            let result = solve(
                |x, y| {
                    for (i, value) in y.iter_mut().enumerate() {
                        *value = scale * x[i];
                    }
                },
                |r, z| z.copy_from_slice(r),
                &[scale, -scale],
                None,
                1e-12,
                8,
                mode(method),
            );
            assert!(result.converged, "{scale:e} {result:?}");
            assert!((result.x[0] - 1.0).abs() < 1e-12, "{result:?}");
            assert!((result.x[1] + 1.0).abs() < 1e-12, "{result:?}");
        }
    }
}

#[test]
fn coupled_huge_spd_systems_avoid_intermediate_product_overflow() {
    let matrix = krylov::Csr::from_triplets(
        2,
        2,
        &[(0, 0, 1e308), (0, 1, 9e307), (1, 0, 9e307), (1, 1, 1e308)],
    );
    let mut b = [0.0; 2];
    matrix.matvec(&[2.0, -2.0], &mut b);
    for method in 0..3 {
        let result = solve(
            |x, y| matrix.matvec(x, y),
            |r, z| z.copy_from_slice(r),
            &b,
            None,
            1e-12,
            12,
            mode(method),
        );
        assert!(result.converged, "{result:?}");
        assert!((result.x[0] - 2.0).abs() < 1e-12, "{result:?}");
        assert!((result.x[1] + 2.0).abs() < 1e-12, "{result:?}");
    }
}
