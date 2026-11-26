//! Numerical optimization drivers built on `basin`.

use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::kernels::optimize::EquationOptimizer;
use rssn::kernels::optimize::LinearRegression;
use rssn::kernels::optimize::MinimizeOptions;
use rssn::kernels::optimize::OptimizationConfig;
use rssn::kernels::optimize::ProblemType;
use rssn::kernels::optimize::Rastrigin;
use rssn::kernels::optimize::ResultAnalyzer;
use rssn::kernels::optimize::Rosenbrock;
use rssn::kernels::optimize::SolverKind;
use rssn::kernels::optimize::Sphere;
use rssn::kernels::optimize::minimize;

fn cfg() -> ProptestConfig {
    ProptestConfig {
        cases: 16,
        rng_seed: RngSeed::Fixed(0x5EED),
        failure_persistence: None,
        ..ProptestConfig::default()
    }
}

fn config(
    problem_type: ProblemType,
    max_iters: u64,
    tolerance: f64,
    dimension: usize,
) -> OptimizationConfig {
    OptimizationConfig {
        problem_type,
        max_iters,
        tolerance,
        dimension,
    }
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
    assert!((Rosenbrock::cost(&[1.0, 1.0])).abs() < 1e-12);
    assert!((Rosenbrock::cost(&[0.0, 0.0]) - 1.0).abs() < 1e-12);
    // f(-1.2, 1) = 4.84 + 100 * (1 - 1.44)^2 = 24.2
    assert!((Rosenbrock::cost(&[-1.2, 1.0]) - 24.2).abs() < 1e-9);
    assert!((Sphere::cost(&[1.0, 2.0, 3.0]) - 14.0).abs() < 1e-12);
    assert!(Rastrigin::cost(&[0.0, 0.0]).abs() < 1e-12);
    assert!(Rastrigin::cost(&[1.0, 1.0]) > 1.0 - 1e-9);

    // The gradient of Rosenbrock vanishes at (1, 1).
    let problem = Rosenbrock::problem(2);
    let (_, gradient) = problem
        .value_and_gradient(&[1.0, 1.0])
        .unwrap_or_else(|e| panic!("{e}"));
    assert!(gradient.iter().all(|value| value.abs() < 1e-9));
}

#[test]
fn bfgs_converges_at_a_stationary_start_with_an_unreachable_target() {
    let cfg = config(ProblemType::Sphere, 50, -1.0, 2);
    let result = EquationOptimizer::solve_with_bfgs(&Sphere::problem(2), &[0.0, 0.0], &cfg)
        .unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(result.termination, "GradientTolerance");
    assert_eq!(result.iterations, 0);
    assert_eq!(result.x, vec![0.0, 0.0]);
    assert_eq!(result.cost, 0.0);
}

#[test]
fn bfgs_returns_the_start_at_a_zero_iteration_limit() {
    let cfg = config(ProblemType::Sphere, 0, 1e-8, 2);
    let initial = vec![1.0, -2.0];
    let result = EquationOptimizer::solve_with_bfgs(&Sphere::problem(2), &initial, &cfg)
        .unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(result.termination, "MaxIter");
    assert_eq!(result.iterations, 0);
    assert_eq!(result.x, initial);
    assert_eq!(result.cost, 5.0);
}

#[test]
fn bfgs_returns_the_start_when_its_cost_overflows_at_a_zero_iteration_limit() {
    let cfg = config(ProblemType::Sphere, 0, 1e-8, 1);
    let initial = vec![1e155];
    let result = EquationOptimizer::solve_with_bfgs(&Sphere::problem(1), &initial, &cfg)
        .unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(result.termination, "MaxIter");
    assert_eq!(result.iterations, 0);
    assert_eq!(result.x, initial);
    assert_eq!(result.cost, f64::INFINITY);
}

#[test]
fn minimize_result_round_trips_with_a_native_fingerprint() {
    let cfg = config(ProblemType::Sphere, 50, 1e-8, 2);
    let result = EquationOptimizer::solve_with_bfgs(&Sphere::problem(2), &[1.0, -2.0], &cfg)
        .unwrap_or_else(|e| panic!("{e}"));
    assert!(result.cost <= cfg.tolerance);
    assert!(result.cost_evals > 0);
    assert_eq!(result.termination, "TargetCost");

    let codec = bincode_next::config::standard()
        .with_deterministic_bincode()
        .with_fingerprint();
    let mut bytes =
        bincode_next::encode_to_vec(&result, codec).unwrap_or_else(|e| panic!("{e}"));
    let (decoded, read) =
        bincode_next::decode_from_slice(&bytes, codec).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(result, decoded);
    assert_eq!(read, bytes.len());

    bytes[0] ^= 1;
    let decoded: Result<(rssn::kernels::optimize::MinimizeResult, usize), _> =
        bincode_next::decode_from_slice(&bytes, codec);
    assert!(matches!(
        decoded,
        Err(bincode_next::error::DecodeError::SchemaMismatch { .. })
    ));
}

#[test]
fn rosenbrock_bfgs_finds_the_minimum() {
    let cfg = config(ProblemType::Rosenbrock, 1000, 1e-8, 2);
    let res = EquationOptimizer::solve_with_bfgs(&Rosenbrock::problem(2), &[-1.2, 1.0], &cfg)
        .unwrap_or_else(|e| panic!("{e}"));
    assert!(res.cost < 1e-4, "cost {}", res.cost);
    assert!((res.x[0] - 1.0).abs() < 0.1, "x = {}", res.x[0]);
    assert!((res.x[1] - 1.0).abs() < 0.1, "y = {}", res.x[1]);
}

#[test]
fn linear_regression_recovers_intercept_and_slope() {
    // y = 2 + 3x ; the parameter vector is [intercept, slope].
    let x = vec![vec![1.0], vec![2.0], vec![3.0], vec![4.0], vec![5.0]];
    let y = vec![5.0, 8.0, 11.0, 14.0, 17.0];
    let problem = LinearRegression::new(x, y)
        .unwrap_or_else(|e| panic!("{e}"))
        .problem();
    let cfg = config(ProblemType::Sphere, 1000, 1e-6, 2);
    let res = EquationOptimizer::solve_with_gradient_descent(&problem, &[0.0, 0.0], &cfg)
        .unwrap_or_else(|e| panic!("{e}"));
    assert!((res.x[0] - 2.0).abs() < 0.1, "intercept = {}", res.x[0]);
    assert!((res.x[1] - 3.0).abs() < 0.1, "slope = {}", res.x[1]);
    assert!(res.cost < 1e-3);
}

#[test]
fn linear_regression_rejects_mismatched_data() {
    assert!(LinearRegression::new(vec![vec![1.0], vec![2.0], vec![3.0]], vec![1.0, 2.0]).is_err());
    assert!(LinearRegression::new(Vec::new(), Vec::new()).is_err());
}

#[test]
fn sphere_gradient_descent_reaches_the_origin() {
    let cfg = config(ProblemType::Sphere, 500, 1e-8, 3);
    let res = EquationOptimizer::solve_with_gradient_descent(
        &Sphere::problem(3),
        &[2.0, -1.5, 3.0],
        &cfg,
    )
    .unwrap_or_else(|e| panic!("{e}"));
    assert!(res.cost < 1e-6);
    assert!(
        res.x.iter().all(|value| value.abs() < 1e-3),
        "best {:?}",
        res.x
    );
}

#[test]
fn convergence_summary_reflects_cost() {
    let cfg = config(ProblemType::Sphere, 500, 1e-8, 2);
    let res = EquationOptimizer::solve_with_bfgs(&Sphere::problem(2), &[1.0, 1.0], &cfg)
        .unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(
        ResultAnalyzer::analyze_convergence(&res),
        "Excellent convergence"
    );
}

#[test]
fn every_solver_kind_runs_on_a_quadratic() {
    // A bounded box around the origin keeps the global solvers honest.
    let bounds = (vec![-5.0, -5.0], vec![5.0, 5.0]);
    for kind in SolverKind::ALL {
        if matches!(kind, SolverKind::Brent | SolverKind::GoldenSection) {
            continue;
        }
        let problem = Sphere::problem(2);
        let options = MinimizeOptions {
            max_iter: 400,
            tolerance: 1e-6,
            seed: 7,
            bounds: Some(bounds.clone()),
            ..MinimizeOptions::default()
        };
        let result = minimize(&problem, &[2.5, -3.0], kind, &options)
            .unwrap_or_else(|e| panic!("{kind:?}: {e}"));
        assert!(
            result.cost <= 1.0,
            "{kind:?} left the cost at {} (x = {:?}, termination {}, {} iterations)",
            result.cost,
            result.x,
            result.termination,
            result.iterations
        );
        // `Direct` publishes its whole search as one initialisation step, so it
        // can legitimately report zero iterations; every solver must still
        // evaluate.
        assert!(
            result.cost_evals > 0,
            "{kind:?} reported no cost evaluations"
        );
    }
}

#[test]
fn scalar_solvers_minimise_one_variable() {
    // (x - 3)^2 has its minimum at 3.
    let problem = rssn::kernels::optimize::FnProblem::new(
        |x| (x[0] - 3.0).powi(2),
        Some(|x: &[f64]| vec![2.0 * (x[0] - 3.0)]),
        1,
        None,
    );
    for kind in [SolverKind::Brent, SolverKind::GoldenSection] {
        let result = minimize(&problem, &[0.0], kind, &MinimizeOptions::default())
            .unwrap_or_else(|e| panic!("{kind:?}: {e}"));
        assert!((result.x[0] - 3.0).abs() < 1e-3, "{kind:?}: {:?}", result.x);
    }
}

#[test]
fn invalid_requests_are_rejected() {
    let problem = Sphere::problem(2);
    let options = MinimizeOptions::default();
    // A bounded solver without bounds.
    assert!(minimize(&problem, &[1.0, 1.0], SolverKind::Bobyqa, &options).is_err());
    // A scalar solver with two variables.
    assert!(minimize(&problem, &[1.0, 1.0], SolverKind::Brent, &options).is_err());
    // An empty start.
    assert!(minimize(&problem, &[], SolverKind::Bfgs, &options).is_err());
    // A gradient solver on a gradient-free problem.
    let gradient_free = rssn::kernels::optimize::FnProblem::new(
        |x| Sphere::cost(x),
        None::<fn(&[f64]) -> Vec<f64>>,
        2,
        None,
    );
    assert!(minimize(&gradient_free, &[1.0, 1.0], SolverKind::Bfgs, &options).is_err());
}

#[test]
fn pso_stays_inside_its_box() {
    let cfg = config(ProblemType::Sphere, 200, 1e-6, 2);
    let bounds = (vec![-1.0, -1.0], vec![1.0, 1.0]);
    let res = EquationOptimizer::solve_with_pso(&Sphere::problem(2), bounds, &cfg)
        .unwrap_or_else(|e| panic!("{e}"));
    assert!(
        res.x.iter().all(|value| value.abs() <= 1.0 + 1e-9),
        "{:?}",
        res.x
    );
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn prop_sphere_converges_from_random_starts(x in -4.0..4.0f64, y in -4.0..4.0f64, z in -4.0..4.0f64) {
        let cfg = config(ProblemType::Sphere, 2000, 1e-4, 3);
        let res = EquationOptimizer::solve_with_gradient_descent(&Sphere::problem(3), &[x, y, z], &cfg)
            .unwrap_or_else(|e| panic!("{e}"));
        prop_assert!(res.cost < 1e-3, "cost {}", res.cost);
    }

    #[test]
    fn prop_bfgs_beats_the_start(x in -2.0..2.0f64, y in -2.0..2.0f64) {
        let cfg = config(ProblemType::Rosenbrock, 500, 1e-6, 2);
        let start = [x, y];
        let before = Rosenbrock::cost(&start);
        let res = EquationOptimizer::solve_with_bfgs(&Rosenbrock::problem(2), &start, &cfg)
            .unwrap_or_else(|e| panic!("{e}"));
        prop_assert!(res.cost <= before + 1e-9, "cost {} from {}", res.cost, before);
    }
}
