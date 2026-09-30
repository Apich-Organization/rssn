//! Linear / non-linear / scalar root solvers (ported from the old
//! `numerical_solve_test.rs`; closure-based API).

use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::kernels::matrix::Matrix;
use rssn::kernels::solve::{
    LinearSolution, solve_linear_system, solve_nonlinear_system, solve_root, solve_root_bisection,
    solve_root_newton,
};

fn cfg() -> ProptestConfig {
    ProptestConfig {
        rng_seed: RngSeed::Fixed(0x5EED),
        failure_persistence: None,
        ..ProptestConfig::default()
    }
}

fn mat_vec(a: &Matrix<f64>, x: &[f64]) -> Vec<f64> {
    (0..a.rows())
        .map(|i| (0..a.cols()).map(|j| a.get(i, j) * x[j]).sum())
        .collect()
}

#[test]
fn linear_unique() {
    // 2x + y = 5 ; x - y = 1  =>  (2, 1)
    let a = Matrix::new(2, 2, vec![2.0, 1.0, 1.0, -1.0]);
    match solve_linear_system(&a, &[5.0, 1.0]).unwrap_or_else(|e| panic!("{e}")) {
        LinearSolution::Unique(x) => {
            assert_approx_eq!(x[0], 2.0);
            assert_approx_eq!(x[1], 1.0);
        }
        other => panic!("expected unique solution, got {other:?}"),
    }
}

#[test]
fn linear_parametric() {
    // x + y + z = 3 twice (rank 1): a plane of solutions.
    let a = Matrix::new(2, 3, vec![1.0, 1.0, 1.0, 2.0, 2.0, 2.0]);
    match solve_linear_system(&a, &[3.0, 6.0]).unwrap_or_else(|e| panic!("{e}")) {
        LinearSolution::Parametric { particular, null_space_basis } => {
            let av = mat_vec(&a, &particular);
            assert_approx_eq!(av[0], 3.0);
            assert_approx_eq!(av[1], 6.0);
            assert_eq!(null_space_basis.cols(), 2, "null space of a rank-1 3-column matrix has dimension 2");
            for col in null_space_basis.get_cols() {
                for v in mat_vec(&a, &col) {
                    assert_approx_eq!(v, 0.0);
                }
            }
        }
        other => panic!("expected parametric solution, got {other:?}"),
    }
}

#[test]
fn linear_no_solution() {
    let a = Matrix::new(2, 2, vec![1.0, 1.0, 1.0, 1.0]);
    assert!(matches!(
        solve_linear_system(&a, &[2.0, 3.0]),
        Ok(LinearSolution::NoSolution)
    ));
}

#[test]
fn linear_dimension_mismatch_is_error() {
    let a = Matrix::new(2, 2, vec![1.0, 0.0, 0.0, 1.0]);
    assert!(solve_linear_system(&a, &[1.0, 2.0, 3.0]).is_err());
}

/// Regression case from the old suite (seed 650, 4x4).
#[test]
fn linear_regression_seed_650() {
    let (a, x) = lcg_system(4, 650);
    let b = mat_vec(&a, &x);
    match solve_linear_system(&a, &b).unwrap_or_else(|e| panic!("{e}")) {
        LinearSolution::Unique(sol) => {
            for (got, want) in mat_vec(&a, &sol).iter().zip(&b) {
                assert!((got - want).abs() < 1e-6, "A*sol = {got}, b = {want}");
            }
        }
        other => panic!("expected unique solution, got {other:?}"),
    }
}

/// Deterministic pseudo-random system used by the old property test.
fn lcg_system(n: usize, seed: u64) -> (Matrix<f64>, Vec<f64>) {
    let mut rng = seed;
    let mut next = || {
        rng = rng.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
        (rng % 100) as f64 / 10.0 - 5.0
    };
    let data: Vec<f64> = (0..n * n).map(|_| next()).collect();
    let x: Vec<f64> = (0..n).map(|_| next()).collect();
    (Matrix::new(n, n, data), x)
}

#[test]
fn nonlinear_system_circle_and_hyperbola() {
    // x^2 + y^2 = 4, x*y = 1
    let f = |x: &[f64], out: &mut [f64]| {
        out[0] = x[0] * x[0] + x[1] * x[1] - 4.0;
        out[1] = x[0] * x[1] - 1.0;
    };
    let r = solve_nonlinear_system(f, &[2.0, 0.5], 1e-12, 50).unwrap_or_else(|e| panic!("{e}"));
    assert_approx_eq!(r[0] * r[0] + r[1] * r[1], 4.0, 1e-9);
    assert_approx_eq!(r[0] * r[1], 1.0, 1e-9);
    // (x + y)^2 = 6, (x - y)^2 = 2 -> x = (sqrt6 + sqrt2)/2 for the branch near (2, 0.5)
    assert_approx_eq!(r[0], (6f64.sqrt() + 2f64.sqrt()) / 2.0, 1e-8);
}

#[test]
fn nonlinear_system_singular_jacobian_is_error() {
    // f(x) = x^2 at x = 0 has zero derivative.
    let r = solve_nonlinear_system(|x: &[f64], out: &mut [f64]| out[0] = x[0] * x[0] + 1.0, &[0.0], 1e-12, 10);
    assert!(r.is_err());
}

#[test]
fn newton_sqrt_two() {
    let r = solve_root_newton(|x| x * x - 2.0, |x| 2.0 * x, 1.0, 1e-14, 50).unwrap_or_else(|e| panic!("{e}"));
    assert_approx_eq!(r, 2f64.sqrt(), 1e-12);
}

#[test]
fn newton_flat_derivative_is_error() {
    assert!(solve_root_newton(|x| x * x - 2.0, |_| 0.0, 1.0, 1e-12, 10).is_err());
}

#[test]
fn bisection_cos_root() {
    let r = solve_root_bisection(f64::cos, (1.0, 2.0), 1e-12, 200).unwrap_or_else(|e| panic!("{e}"));
    assert_approx_eq!(r, std::f64::consts::FRAC_PI_2, 1e-9);
}

#[test]
fn bisection_requires_bracket() {
    assert!(solve_root_bisection(|x| x * x + 1.0, (-1.0, 1.0), 1e-12, 50).is_err());
}

#[test]
fn root_cubic_and_no_real_root() {
    let f = |x: f64| x * x * x - 2.0 * x - 5.0;
    let r = solve_root(f, 2.0, 1e-12, 100).unwrap_or_else(|e| panic!("{e}"));
    assert!(f(r).abs() < 1e-9);
    assert_approx_eq!(r, 2.094_551_481_542_326_7, 1e-9);
    assert!(solve_root(|x| x * x + 1.0, 0.3, 1e-12, 50).is_err());
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn prop_linear_system_solution_satisfies_equations(n in 2..5usize, seed in 0..1000u64) {
        let (a, x) = lcg_system(n, seed);
        let b = mat_vec(&a, &x);
        match solve_linear_system(&a, &b) {
            Ok(LinearSolution::Unique(sol)) => {
                let scale = b.iter().fold(1.0_f64, |m, v| m.max(v.abs()));
                for (got, want) in mat_vec(&a, &sol).iter().zip(&b) {
                    prop_assert!((got - want).abs() < 1e-6 * scale, "A*sol = {got}, b = {want}");
                }
            }
            Ok(LinearSolution::Parametric { particular, null_space_basis }) => {
                for (got, want) in mat_vec(&a, &particular).iter().zip(&b) {
                    prop_assert!((got - want).abs() < 1e-6, "particular: A*p = {got}, b = {want}");
                }
                for col in null_space_basis.get_cols() {
                    for v in mat_vec(&a, &col) {
                        prop_assert!(v.abs() < 1e-6, "null vector gives {v}");
                    }
                }
            }
            Ok(LinearSolution::NoSolution) => {}
            Err(e) => return Err(TestCaseError::fail(e)),
        }
    }

    #[test]
    fn prop_newton_finds_cube_roots(c in 0.5..50.0f64) {
        let r = solve_root_newton(|x| x * x * x - c, |x| 3.0 * x * x, c.max(1.0), 1e-13, 200)
            .map_err(TestCaseError::fail)?;
        prop_assert!((r - c.cbrt()).abs() < 1e-9);
        prop_assert!((r * r * r - c).abs() < 1e-9);
    }

    #[test]
    fn prop_bisection_root_has_small_residual(root in -5.0..5.0f64, slope in 0.5..4.0f64) {
        let f = move |x: f64| slope * (x - root);
        let r = solve_root_bisection(f, (root - 3.0, root + 4.0), 1e-12, 200).map_err(TestCaseError::fail)?;
        prop_assert!(f(r).abs() < 1e-9);
        prop_assert!((r - root).abs() < 1e-9);
    }

    #[test]
    fn prop_solve_root_returns_small_residual_for_odd_cubic(a in 0.5..3.0f64, shift in -4.0..4.0f64) {
        let f = move |x: f64| a * (x - shift).powi(3) + (x - shift);
        let r = solve_root(f, 0.0, 1e-12, 200).map_err(TestCaseError::fail)?;
        prop_assert!(f(r).abs() < 1e-8, "|f({r})| = {}", f(r).abs());
    }

    #[test]
    fn prop_nonlinear_system_recovers_planted_solution(x0 in 0.5..3.0f64, y0 in 0.5..3.0f64) {
        // F(x, y) = (x + y - s1, x*y - s2) has (x0, y0) as a root; start nearby.
        let (s1, s2) = (x0 + y0, x0 * y0);
        let f = move |x: &[f64], out: &mut [f64]| {
            out[0] = x[0] + x[1] - s1;
            out[1] = x[0] * x[1] - s2;
        };
        let r = solve_nonlinear_system(f, &[x0 + 0.05, y0 - 0.05], 1e-12, 100);
        // The Jacobian is singular when x0 == y0; skip that degenerate case.
        prop_assume!((x0 - y0).abs() > 0.2);
        let r = r.map_err(TestCaseError::fail)?;
        let mut res = [0.0; 2];
        f(&r, &mut res);
        prop_assert!(res[0].abs() < 1e-8 && res[1].abs() < 1e-8);
    }
}
