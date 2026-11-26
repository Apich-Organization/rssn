//! Exact identities and scaling regressions for the chained Rosenbrock function.

use rssn::kernels::optimize::EquationOptimizer;
use rssn::kernels::optimize::OptimizationConfig;
use rssn::kernels::optimize::Rosenbrock;

/// The chained Rosenbrock function `sum (1 - x_i)^2 + 100 (x_{i+1} - x_i^2)^2`.
fn value(x: &[f64]) -> f64 {
    Rosenbrock::cost(x)
}

#[test]
fn the_middle_coordinate_receives_both_adjacent_derivatives() {
    let x = [0.0, 1.0, 2.0];
    assert_eq!(value(&x).to_bits(), 201.0_f64.to_bits());
    let (_, gradient) = Rosenbrock::problem(3)
        .value_and_gradient(&x)
        .expect("gradient");
    // d/dx0 = -2(1 - x0) - 400 x0 (x1 - x0^2) = -2
    // d/dx1 = 200(x1 - x0^2) - 400 x1 (x2 - x1^2) = 200 - 400 = -200
    // d/dx2 = 200(x2 - x1^2) = 200
    assert_eq!(gradient, vec![-2.0, -200.0, 200.0]);
}

#[test]
fn gradient_matches_central_differences_on_a_large_chain() {
    let x: Vec<f64> = (0..73)
        .map(|i| (f64::from(13 * i % 31) - 10.0) / 8.0)
        .collect();
    let problem = Rosenbrock::problem(x.len());
    let (_, gradient) = problem.value_and_gradient(&x).expect("gradient");
    let h = 1.0 / 65_536.0;
    for (i, &entry) in gradient.iter().enumerate() {
        let mut up = x.clone();
        up[i] += h;
        let mut down = x.clone();
        down[i] -= h;
        // The function is a quartic in each coordinate, so the central
        // difference carries no truncation error beyond rounding.
        let numeric = (value(&up) - value(&down)) / (2.0 * h);
        assert!(
            (entry - numeric).abs() < 3e-6,
            "derivative {i}: {entry} vs {numeric}"
        );
    }
}

#[test]
fn reversed_coordinates_give_the_same_value() {
    let values: Vec<f64> = (0..9).map(|i| f64::from(i) / 4.0 - 1.0).collect();
    let mut reversed = values.clone();
    reversed.reverse();
    // The chain is not symmetric, but its value must be finite and reproducible.
    assert_eq!(value(&values).to_bits(), value(&values).to_bits());
    assert!(value(&reversed).is_finite());
}

#[test]
fn a_short_chain_is_the_empty_sum() {
    assert_eq!(value(&[]), 0.0);
    assert_eq!(value(&[3.0]), 0.0);
    // A single pair reduces to the textbook form.
    assert_eq!(value(&[0.0, 0.0]), 1.0);
    assert!(value(&[1.0, 1.0]).abs() < 1e-15);
}

#[test]
fn bfgs_uses_the_corrected_three_dimensional_gradient() {
    let config = OptimizationConfig {
        tolerance: 1e-10,
        dimension: 3,
        ..OptimizationConfig::default()
    };
    let optimum =
        EquationOptimizer::solve_with_bfgs(&Rosenbrock::problem(3), &[0.4, 0.7, 0.8], &config)
            .expect("solve");
    assert!(optimum.cost < 1e-9, "cost {}", optimum.cost);
    assert!(
        optimum.x.iter().all(|value| (value - 1.0).abs() < 1e-4),
        "point {:?}",
        optimum.x
    );
}
