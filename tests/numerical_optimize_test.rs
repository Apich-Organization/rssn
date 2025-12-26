use argmin::core::Executor;
use argmin::core::IterState;
use argmin::core::State;
use argmin::solver::linesearch::MoreThuenteLineSearch;
use argmin::solver::quasinewton::BFGS;
use ndarray::Array1;
use ndarray::Array2;
use proptest::prelude::*;
use rssn::numerical::optimize::*;

type P = Array1<f64>;

type G = Array1<f64>;

type F = f64;

type H = Array2<f64>;

#[allow(dead_code)]

type BFGSState =
    IterState<P, G, (), H, (), F>;

#[test]
#[allow(unused_variables)]

fn test_rosenbrock_optimization() {

    let config = OptimizationConfig {
        problem_type:
            ProblemType::Rosenbrock,
        max_iters : 1000,
        tolerance : 1e-8,
        dimension : 2,
    };

    let problem = Rosenbrock::default();

    let initial_guess =
        Array1::from(vec![-1.2, 1.0]);

    let linesearch =
        MoreThuenteLineSearch::new();

    let solver = BFGS::new(linesearch);

    let dimension = initial_guess.len();

    let initial_inverse_hessian =
        Array2::eye(dimension);

    let initial_state : BFGSState =
        IterState::new()
            .param(initial_guess)
            .inv_hessian(
                initial_inverse_hessian,
            )
            .max_iters(config.max_iters)
            .target_cost(
                config.tolerance,
            );

    let result =
        Executor::new(problem, solver)
            .configure(
                |state : BFGSState| {

                    initial_state
                },
            )
            .run()
            .unwrap();

    let best_param = result
        .state
        .get_best_param()
        .unwrap();

    let best_cost = result
        .state
        .get_best_cost();

    assert!(best_cost < 1e-4);

    assert!(
        (best_param[0] - 1.0).abs()
            < 0.1
    );

    assert!(
        (best_param[1] - 1.0).abs()
            < 0.1
    );
}

#[test]

fn test_linear_regression() {

    let x = Array2::from_shape_vec(
        (5, 1),
        vec![
            1.0, 2.0, 3.0, 4.0, 5.0,
        ],
    )
    .unwrap();

    let y = Array1::from(vec![
        5.0, 8.0, 11.0, 14.0, 17.0,
    ]);

    let problem =
        match LinearRegression::new(
            x, y,
        ) {
            | Ok(p) => p,
            | Err(e) => {

                panic!(
                    "Failed to create \
                     LinearRegression \
                     problem: {}",
                    e
                )
            },
        };

    let config = OptimizationConfig {
        problem_type:
            ProblemType::Custom,
        max_iters : 1000,
        tolerance : 1e-6,
        dimension : 2,
    };

    let result = match EquationOptimizer::solve_with_gradient_descent(
        problem,
        Array1::from(vec![0.0, 0.0]),
        &config,
    ) {
        | Ok(r) => r,
        | Err(e) => {

            panic!(
                "Solver failed for linear regression: {}",
                e
            )
        },
    };

    let best_param = match result
        .state
        .get_best_param()
    {
        | Some(p) => p,
        | None => {

            panic!(
                "Best param should \
                 not be None after \
                 successful \
                 optimization"
            )
        },
    };

    assert!(
        (best_param[0] - 2.0).abs()
            < 0.5
    );

    assert!(
        (best_param[1] - 3.0).abs()
            < 0.5
    );
}

#[test]

fn test_sphere_function() {

    let config = OptimizationConfig {
        problem_type:
            ProblemType::Sphere,
        max_iters : 500,
        tolerance : 1e-8,
        dimension : 3,
    };

    let problem = Sphere;

    let initial_guess =
        Array1::from(vec![
            2.0, -1.5, 3.0,
        ]);

    let result = EquationOptimizer::solve_with_gradient_descent(
        problem,
        initial_guess,
        &config,
    )
    .unwrap();

    let best_cost = result
        .state
        .get_best_cost();

    assert!(best_cost < 1e-6);
}

proptest! {
    #[test]
    fn test_sphere_convergence_prop(x in -4.0..4.0f64, y in -4.0..4.0f64, z in -4.0..4.0f64) {
        let config = OptimizationConfig {
            problem_type: ProblemType::Sphere,
            max_iters: 2000,
            tolerance: 1e-4, // Relaxed tolerance for random starts with GD
            dimension: 3,
        };
        let problem = Sphere;
        let initial_guess = Array1::from(vec![x, y, z]);

        let result = EquationOptimizer::solve_with_gradient_descent(problem, initial_guess, &config);

        prop_assert!(result.is_ok());
        let res = result.unwrap();
        // Since GD can get stuck or be slow, check if cost decreased significantly or is low.
        // Sphere is convex, so it should be low.
        prop_assert!(res.state.get_best_cost() < 0.1, "Best cost {} is too high for start {:?}", res.state.get_best_cost(), vec![x, y, z]);
    }
}
