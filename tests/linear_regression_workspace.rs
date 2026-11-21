//! Regression tests for the linear least-squares cost and gradient.

use argmin::core::CostFunction;
use argmin::core::Gradient;
use ndarray::Array1;
use ndarray::Array2;
use rssn::kernels::optimize::LinearRegression;

#[test]
fn cost_and_gradient_follow_logical_order_for_strided_inputs() {
    let x = Array2::from_shape_fn((5, 3), |(i, j)| (i as f64 - 1.0) * (j as f64 - 0.5));
    let x = x.slice_move(ndarray::s![..;-1, ..;-1]);
    let y = Array1::from_vec(vec![0.0, 1.5, -2.0, 1.0, 3.0]).slice_move(ndarray::s![..;-1]);
    let params = Array1::from_vec(vec![1.5, -2.0, 0.5, 1.0]).slice_move(ndarray::s![..;-1]);
    let model = LinearRegression::new(x, y).unwrap();

    let mut squares = 0.0;
    let mut gradient = [0.0; 4];
    for row in 0..5 {
        let mut residual = params[0] - model.y[row];
        for col in 0..3 {
            residual += params[col + 1] * model.x[(row, col)];
        }
        squares += residual * residual;
        gradient[0] += residual;
        for col in 0..3 {
            gradient[col + 1] += residual * model.x[(row, col)];
        }
    }
    let expected_cost = squares / 10.0;
    let expected_gradient = Array1::from_vec(gradient.map(|v| v / 5.0).to_vec());
    let (cost, slope) = model.value_and_gradient(&params).unwrap();
    assert_eq!(cost.to_bits(), expected_cost.to_bits());
    assert_eq!(slope, expected_gradient);
    assert_eq!(model.cost(&params).unwrap().to_bits(), cost.to_bits());
    assert_eq!(model.gradient(&params).unwrap(), slope);
}

#[test]
fn normalization_preserves_finite_cost_when_raw_squares_overflow_or_underflow() {
    let params = Array1::zeros(2);
    for target in [1e154_f64, 2.5e-162_f64] {
        let model =
            LinearRegression::new(Array2::zeros((4, 1)), Array1::from_elem(4, target)).unwrap();
        let expected = (0.5 * target) * target;
        let (cost, grad) = model.value_and_gradient(&params).unwrap();
        assert!(expected.is_finite() && expected > 0.0);
        assert_eq!(cost.to_bits(), expected.to_bits());
        assert_eq!(model.cost(&params).unwrap().to_bits(), cost.to_bits());
        assert_eq!(grad, Array1::from_vec(vec![-target, 0.0]));
    }
}

#[test]
fn normalization_preserves_representable_signed_gradient_means() {
    let x = Array2::from_shape_fn((4, 2), |(_, j)| if j == 0 { 1.0 } else { -0.5 });
    let model = LinearRegression::new(x, Array1::from_elem(4, -1e308)).unwrap();
    let params = Array1::zeros(3);
    let (cost, gradient) = model.value_and_gradient(&params).unwrap();
    assert!(cost.is_infinite() && cost.is_sign_positive());
    assert_eq!(gradient, Array1::from_vec(vec![1e308, 1e308, -5e307]));
    assert_eq!(model.gradient(&params).unwrap(), gradient);
}

#[test]
fn malformed_model_and_coefficient_shapes_report_errors() {
    let good = LinearRegression::new(Array2::zeros((3, 1)), Array1::zeros(3)).unwrap();
    for count in [0, 1, 3] {
        let coefficients = Array1::zeros(count);
        assert!(good.cost(&coefficients).is_err());
        assert!(good.gradient(&coefficients).is_err());
        assert!(good.value_and_gradient(&coefficients).is_err());
    }
    let malformed = LinearRegression {
        x: Array2::zeros((3, 1)),
        y: Array1::zeros(1),
    };
    let coefficients = Array1::zeros(2);
    assert!(malformed.cost(&coefficients).is_err());
    assert!(malformed.gradient(&coefficients).is_err());
    assert!(malformed.value_and_gradient(&coefficients).is_err());
}

#[test]
fn concurrent_reuse_overwrites_the_entire_residual_buffer() {
    let rows = 1024;
    let x = Array2::from_shape_fn((rows, 1), |(i, _)| if i % 2 == 0 { -1.0 } else { 1.0 });
    let y = Array1::from_iter(x.column(0).iter().map(|v| 2.0 + 3.0 * v));
    let model = LinearRegression::new(x, y).unwrap();

    std::thread::scope(|threads| {
        for phase in 0..4 {
            let model = &model;
            threads.spawn(move || {
                for i in 0..16 {
                    let delta = (phase * 16 + i) as f64 * 0.25;
                    let p = Array1::from_vec(vec![2.0 - delta, 3.0 + 2.0 * delta]);
                    let (cost, gradient) = model.value_and_gradient(&p).unwrap();
                    assert_eq!(cost.to_bits(), (2.5 * delta * delta).to_bits());
                    assert_eq!(gradient, Array1::from_vec(vec![-delta, 2.0 * delta]));
                }
            });
        }
    });
}

#[test]
fn contiguous_multi_feature_gradient_matches_central_differences() {
    let x = Array2::from_shape_fn((257, 11), |(i, j)| {
        ((i * 31 + j * 19) as f64 % 13.0 - 6.0) / 8.0
    });
    let y = Array1::from_shape_fn(257, |i| (i as f64 % 7.0) / 16.0);
    let parameters = Array1::from_shape_fn(12, |j| j as f64 / 16.0);
    let model = LinearRegression::new(x, y).unwrap();
    let (cost, gradient) = model.value_and_gradient(&parameters).unwrap();
    assert_eq!(model.cost(&parameters).unwrap().to_bits(), cost.to_bits());
    assert_eq!(model.gradient(&parameters).unwrap(), gradient);
    let delta = 1.0 / 1024.0;
    for (j, &entry) in gradient.iter().enumerate() {
        let mut trial = parameters.clone();
        trial[j] += delta;
        let upper = model.cost(&trial).unwrap();
        trial[j] -= 2.0 * delta;
        let lower = model.cost(&trial).unwrap();
        // The central difference of a quadratic has no truncation error.
        assert!(((upper - lower) / (2.0 * delta) - entry).abs() < 2e-13);
    }
}
