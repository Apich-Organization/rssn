//! Gradient-based optimisation drivers built on `argmin` (ported from
//! `numerical_optimize_test.rs`).

use argmin::core::{CostFunction, Gradient, State};
use ndarray::{Array1, Array2};
use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::kernels::optimize::{
    EquationOptimizer, LinearRegression, OptimizationConfig, ProblemType, Rastrigin, ResultAnalyzer,
    Rosenbrock, Sphere,
};

fn cfg() -> ProptestConfig {
    ProptestConfig {
        cases: 32,
        rng_seed: RngSeed::Fixed(0x5EED),
        failure_persistence: None,
        ..ProptestConfig::default()
    }
}

fn config(problem_type: ProblemType, max_iters: u64, tolerance: f64, dimension: usize) -> OptimizationConfig {
    OptimizationConfig { problem_type, max_iters, tolerance, dimension }
}

#[test]
fn default_config_values() {
    let c = OptimizationConfig::default();
    assert_eq!(c.max_iters, 1000);
    assert_eq!(c.dimension, 2);
    assert!((c.tolerance - 1e-6).abs() < 1e-18);
    assert!(matches!(c.problem_type, ProblemType::Rosenbrock));
}

#[test]
fn test_function_values_at_known_points() {
    let ros = Rosenbrock::default();
    assert_eq!(ros.cost(&Array1::from(vec![1.0, 1.0])).unwrap_or(f64::NAN), 0.0);
    assert_eq!(ros.cost(&Array1::from(vec![0.0, 0.0])).unwrap_or(f64::NAN), 1.0);
    let g = ros.gradient(&Array1::from(vec![1.0, 1.0])).unwrap_or_else(|e| panic!("{e}"));
    assert!(g.iter().all(|v| v.abs() < 1e-12), "gradient at the minimum: {g:?}");
    // f(-1.2, 1) = 4.84 + 100 * (1 - 1.44)^2 = 24.2
    let c = ros.cost(&Array1::from(vec![-1.2, 1.0])).unwrap_or(f64::NAN);
    assert!((c - 24.2).abs() < 1e-9, "cost = {c}");

    assert_eq!(Sphere.cost(&Array1::from(vec![1.0, 2.0, 3.0])).unwrap_or(f64::NAN), 14.0);
    assert_eq!(
        Sphere.gradient(&Array1::from(vec![1.0, -2.0])).unwrap_or_default(),
        Array1::from(vec![2.0, -4.0])
    );
    let ras = Rastrigin::default();
    assert!(ras.cost(&Array1::from(vec![0.0, 0.0])).unwrap_or(f64::NAN).abs() < 1e-12);
    assert!(ras.cost(&Array1::from(vec![1.0, 1.0])).unwrap_or(f64::NAN) > 1.0 - 1e-9);
}

#[test]
fn rosenbrock_bfgs_finds_the_minimum() {
    let cfg = config(ProblemType::Rosenbrock, 1000, 1e-8, 2);
    let res = EquationOptimizer::solve_with_bfgs(Rosenbrock::default(), Array1::from(vec![-1.2, 1.0]), &cfg)
        .unwrap_or_else(|e| panic!("{e}"));
    let best = res.state.get_best_param().cloned().unwrap_or_default();
    assert!(res.state.get_best_cost() < 1e-4, "cost {}", res.state.get_best_cost());
    assert!((best[0] - 1.0).abs() < 0.1, "x = {}", best[0]);
    assert!((best[1] - 1.0).abs() < 0.1, "y = {}", best[1]);
}

#[test]
fn linear_regression_recovers_intercept_and_slope() {
    // y = 2 + 3x ; the parameter vector is [intercept, slope].
    let x = Array2::from_shape_vec((5, 1), vec![1.0, 2.0, 3.0, 4.0, 5.0]).unwrap_or_else(|e| panic!("{e}"));
    let y = Array1::from(vec![5.0, 8.0, 11.0, 14.0, 17.0]);
    let problem = LinearRegression::new(x, y).unwrap_or_else(|e| panic!("{e}"));
    let cfg = config(ProblemType::Sphere, 1000, 1e-6, 2);
    let res = EquationOptimizer::solve_with_gradient_descent(problem, Array1::from(vec![0.0, 0.0]), &cfg)
        .unwrap_or_else(|e| panic!("{e}"));
    let best = res.state.get_best_param().cloned().unwrap_or_default();
    assert!((best[0] - 2.0).abs() < 0.1, "intercept = {}", best[0]);
    assert!((best[1] - 3.0).abs() < 0.1, "slope = {}", best[1]);
    assert!(res.state.get_best_cost() < 1e-3);
}

#[test]
fn linear_regression_rejects_mismatched_data() {
    let x = Array2::from_shape_vec((3, 1), vec![1.0, 2.0, 3.0]).unwrap_or_else(|e| panic!("{e}"));
    assert!(LinearRegression::new(x, Array1::from(vec![1.0, 2.0])).is_err());
}

#[test]
fn sphere_gradient_descent_reaches_the_origin() {
    let cfg = config(ProblemType::Sphere, 500, 1e-8, 3);
    let res = EquationOptimizer::solve_with_gradient_descent(Sphere, Array1::from(vec![2.0, -1.5, 3.0]), &cfg)
        .unwrap_or_else(|e| panic!("{e}"));
    assert!(res.state.get_best_cost() < 1e-6);
    let best = res.state.get_best_param().cloned().unwrap_or_default();
    assert!(best.iter().all(|v| v.abs() < 1e-3), "best {best:?}");
}

#[test]
fn convergence_summary_reflects_cost() {
    let cfg = config(ProblemType::Sphere, 500, 1e-8, 2);
    let res = EquationOptimizer::solve_with_bfgs(Sphere, Array1::from(vec![1.0, 1.0]), &cfg)
        .unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(ResultAnalyzer::analyze_convergence(&res.state), "Excellent convergence");
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn prop_sphere_converges_from_random_starts(x in -4.0..4.0f64, y in -4.0..4.0f64, z in -4.0..4.0f64) {
        let cfg = config(ProblemType::Sphere, 2000, 1e-4, 3);
        let res = EquationOptimizer::solve_with_gradient_descent(Sphere, Array1::from(vec![x, y, z]), &cfg)
            .map_err(|e| TestCaseError::fail(e.to_string()))?;
        prop_assert!(res.state.get_best_cost() < 0.1, "cost {} from {:?}", res.state.get_best_cost(), (x, y, z));
    }

    #[test]
    fn prop_bfgs_never_increases_cost_of_sphere(x in -5.0..5.0f64, y in -5.0..5.0f64) {
        let cfg = config(ProblemType::Sphere, 100, 1e-10, 2);
        let res = EquationOptimizer::solve_with_bfgs(Sphere, Array1::from(vec![x, y]), &cfg)
            .map_err(|e| TestCaseError::fail(e.to_string()))?;
        prop_assert!(res.state.get_best_cost() <= x * x + y * y + 1e-12);
    }
}
