//! Exact identities and scaling regressions for the chained Rosenbrock function.

use argmin::core::CostFunction;
use argmin::core::Gradient;
use argmin::core::State;
use ndarray::Array1;
use rssn::kernels::optimize::EquationOptimizer;
use rssn::kernels::optimize::OptimizationConfig;
use rssn::kernels::optimize::Rosenbrock;

#[test]
fn the_middle_coordinate_receives_both_adjacent_derivatives() {
    let objective = Rosenbrock::default();
    let x = Array1::from_vec(vec![0.0, 1.0, 2.0]);
    assert_eq!(objective.cost(&x).unwrap().to_bits(), 201.0_f64.to_bits());
    assert_eq!(
        objective.gradient(&x).unwrap(),
        Array1::from_vec(vec![-2.0, -200.0, 200.0])
    );
}

#[test]
fn strided_coordinates_match_independent_partial_derivatives() {
    let objective = Rosenbrock { a: -0.5, b: 3.0 };
    let coordinates = Array1::from_vec(vec![-0.75, 0.0, 0.5, 0.0, 1.25, 0.0, -1.5, 0.0, 0.25]);
    let x = coordinates.slice_move(ndarray::s![..;-2]);
    let evaluate = |x: &Array1<f64>| {
        x.iter()
            .zip(x.iter().skip(1))
            .map(|(&a, &b)| (objective.a - a).powi(2) + objective.b * (b - a * a).powi(2))
            .sum::<f64>()
    };
    let gradient = objective.gradient(&x).unwrap();
    assert_eq!(
        objective.cost(&x).unwrap().to_bits(),
        evaluate(&x).to_bits()
    );
    let h = 1.0 / 65_536.0;
    for (i, value) in gradient.iter().enumerate() {
        // A five-point stencil is exact on this quartic in a single coordinate.
        let sample = |offset: f64| {
            let mut shifted = x.clone();
            shifted[i] += offset * h;
            evaluate(&shifted)
        };
        let exact =
            (sample(-2.0) - 8.0 * sample(-1.0) + 8.0 * sample(1.0) - sample(2.0)) / (12.0 * h);
        assert!(
            (value - exact).abs() < 3e-9,
            "derivative {i}: {value} vs {exact}"
        );
    }
    for length in [0, 1] {
        assert!(objective.cost(&Array1::zeros(length)).is_err());
        assert!(objective.gradient(&Array1::zeros(length)).is_err());
    }
}

#[test]
fn large_contiguous_and_reversed_chains_match_coordinatewise_derivatives() {
    let objective = Rosenbrock { a: -0.25, b: -2.0 };
    let values = Array1::from_shape_fn(73, |i| ((13 * i % 31) as f64 - 10.0) / 8.0);
    for x in [values.clone(), values.slice_move(ndarray::s![..;-1])] {
        let expected = Array1::from_shape_fn(x.len(), |i| {
            let incoming = if i == 0 {
                0.0
            } else {
                2.0 * objective.b * (x[i] - x[i - 1] * x[i - 1])
            };
            let outgoing = if i + 1 == x.len() {
                0.0
            } else {
                2.0 * (x[i] - objective.a) - 4.0 * objective.b * x[i] * (x[i + 1] - x[i] * x[i])
            };
            incoming + outgoing
        });
        assert_eq!(objective.gradient(&x).unwrap(), expected);
    }
}

#[test]
fn finite_penalties_survive_both_raw_square_extremes() {
    for sign in [-1, 1] {
        let b = 2.0_f64.powi(-1000 * sign);
        let next = 2.0_f64.powi(600 * sign);
        let objective = Rosenbrock { a: 0.0, b };
        let x = Array1::from_vec(vec![0.0, next]);
        let expected = 2.0_f64.powi(200 * sign);
        assert_eq!(objective.cost(&x).unwrap().to_bits(), expected.to_bits());
        let right = 2.0_f64.powi(1 - 400 * sign);
        assert_eq!(
            objective.gradient(&x).unwrap(),
            Array1::from_vec(vec![0.0, right])
        );
    }
    let objective = Rosenbrock {
        a: 0.0,
        b: 2.0_f64.powi(100),
    };
    let x = Array1::from_vec(vec![0.0, 2.0_f64.powi(-545)]);
    assert_eq!(
        objective.cost(&x).unwrap().to_bits(),
        2.0_f64.powi(-990).to_bits()
    );
}

#[test]
fn weighted_derivatives_survive_unweighted_residual_overflow_and_underflow() {
    // Linear terms vanish because the left coordinate is a. The restoring
    // scales can be checked by adding the exact exponents in each monomial.
    for (a, b, expected_cost, left, right) in [
        (
            2.0_f64.powi(520),
            f64::from_bits(16),
            2.0_f64.powi(1010),
            2.0_f64.powi(492),
            -2.0_f64.powi(-29),
        ),
        (
            2.0_f64.powi(-560),
            2.0_f64.powi(1023),
            0.0,
            2.0_f64.powi(-655),
            -2.0_f64.powi(-96),
        ),
    ] {
        let objective = Rosenbrock { a, b };
        let x = Array1::from_vec(vec![a, 0.0]);
        assert_eq!(
            objective.cost(&x).unwrap().to_bits(),
            expected_cost.to_bits()
        );
        assert_eq!(
            objective.gradient(&x).unwrap(),
            Array1::from_vec(vec![left, right])
        );
    }
}

#[test]
fn a_large_weight_forms_finite_derivatives_and_zero_weight_skips_coupling() {
    let objective = Rosenbrock {
        a: 1.0,
        b: 2.0_f64.powi(1023),
    };
    let x = Array1::from_vec(vec![1.0, 1.0 + 2.0_f64.powi(-40)]);
    assert_eq!(
        objective.gradient(&x).unwrap(),
        Array1::from_vec(vec![-2.0_f64.powi(985), 2.0_f64.powi(984)])
    );

    let a = 2.0_f64.powi(600);
    let objective = Rosenbrock { a, b: 0.0 };
    let x = Array1::from_vec(vec![a, 0.0]);
    assert_eq!(objective.cost(&x).unwrap().to_bits(), 0);
    assert_eq!(objective.gradient(&x).unwrap(), Array1::<f64>::zeros(2));
}

#[test]
fn bfgs_uses_the_corrected_three_dimensional_gradient() {
    let config = OptimizationConfig {
        tolerance: 1e-10,
        dimension: 3,
        ..OptimizationConfig::default()
    };
    let optimum = EquationOptimizer::solve_with_bfgs(
        Rosenbrock::default(),
        Array1::from_vec(vec![0.4, 0.7, 0.8]),
        &config,
    )
    .unwrap();
    assert!(optimum.state.get_best_cost() < 1e-9);
    let point = optimum.state.get_best_param().unwrap();
    assert!(
        point.iter().all(|v| (v - 1.0).abs() < 1e-4),
        "point {point:?}"
    );
}
