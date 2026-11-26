//! Regression tests for the linear least-squares cost and gradient.

use rssn::kernels::optimize::LinearRegression;

/// Builds a model from a row-major design matrix.
fn model(
    x: Vec<Vec<f64>>,
    y: Vec<f64>,
) -> LinearRegression {
    LinearRegression::new(x, y).expect("model")
}

#[test]
fn cost_and_gradient_follow_logical_order() {
    let design: Vec<Vec<f64>> = (0..5)
        .map(|i| {
            (0..3)
                .map(|j| (f64::from(i) - 1.0) * (f64::from(j) - 0.5))
                .collect()
        })
        .collect();
    let y = vec![0.0, 1.5, -2.0, 1.0, 3.0];
    let params = vec![1.5, -2.0, 0.5, 1.0];
    let model = model(design, y);

    let mut squares = 0.0;
    let mut gradient = [0.0; 4];
    for (row, &response) in model.x.iter().zip(&model.y) {
        let mut residual = params[0] - response;
        for (col, &value) in row.iter().enumerate() {
            residual += params[col + 1] * value;
        }
        squares += residual * residual;
        gradient[0] += residual;
        for (col, &value) in row.iter().enumerate() {
            gradient[col + 1] += residual * value;
        }
    }
    let count = model.y.len() as f64;
    let expected_cost = squares / (2.0 * count);
    let expected_gradient: Vec<f64> = gradient.iter().map(|value| value / count).collect();
    let (cost, slope) = model.value_and_gradient(&params).expect("value");
    assert_eq!(cost.to_bits(), expected_cost.to_bits());
    assert_eq!(slope, expected_gradient);
}

#[test]
fn malformed_model_and_coefficient_shapes_report_errors() {
    let good = model(vec![vec![0.0]; 3], vec![0.0; 3]);
    for count in [0, 1, 3] {
        let coefficients = vec![0.0; count];
        assert!(good.value_and_gradient(&coefficients).is_err());
    }
    assert!(LinearRegression::new(vec![vec![0.0; 1]; 3], vec![0.0; 1]).is_err());
    assert!(LinearRegression::new(Vec::new(), Vec::new()).is_err());
}

#[test]
fn concurrent_evaluations_agree() {
    let rows = 1024;
    let design: Vec<Vec<f64>> = (0..rows)
        .map(|i| vec![if i % 2 == 0 { -1.0 } else { 1.0 }])
        .collect();
    let y: Vec<f64> = design.iter().map(|row| 2.0 + 3.0 * row[0]).collect();
    let model = model(design, y);

    std::thread::scope(|threads| {
        for phase in 0..4 {
            let model = &model;
            threads.spawn(move || {
                for i in 0..16 {
                    let delta = f64::from(phase * 16 + i) * 0.25;
                    let p = vec![2.0 - delta, 3.0 + 2.0 * delta];
                    let (cost, gradient) = model.value_and_gradient(&p).expect("value");
                    // The residual is `delta * (2x - 1)` with x = -1 or 1, so
                    // the mean square is `5 delta^2` and the half mean square
                    // is `2.5 delta^2`. The mean residual is `-delta` and the
                    // mean residual times x is `2 delta`.
                    let expected = 2.5 * delta * delta;
                    assert_eq!(cost.to_bits(), expected.to_bits());
                    assert!((gradient[0] + delta).abs() < 1e-9, "{:?}", gradient);
                    assert!((gradient[1] - 2.0 * delta).abs() < 1e-9, "{:?}", gradient);
                }
            });
        }
    });
}

#[test]
fn multi_feature_gradient_matches_central_differences() {
    let design: Vec<Vec<f64>> = (0..257)
        .map(|i| {
            (0..11)
                .map(|j| (f64::from((i * 31 + j * 19) % 13) - 6.0) / 8.0)
                .collect()
        })
        .collect();
    let y: Vec<f64> = (0..257).map(|i| f64::from(i % 7) / 16.0).collect();
    let parameters: Vec<f64> = (0..12).map(|j| f64::from(j) / 16.0).collect();
    let model = model(design, y);
    let (cost, gradient) = model.value_and_gradient(&parameters).expect("value");
    assert_eq!(gradient.len(), 12);
    assert!(cost >= 0.0);

    let delta = 1.0 / 1024.0;
    for (j, &entry) in gradient.iter().enumerate() {
        let mut trial = parameters.clone();
        trial[j] += delta;
        let upper = model.value_and_gradient(&trial).expect("upper").0;
        trial[j] -= 2.0 * delta;
        let lower = model.value_and_gradient(&trial).expect("lower").0;
        // The central difference of a quadratic has no truncation error.
        let numeric = (upper - lower) / (2.0 * delta);
        assert!(
            (numeric - entry).abs() < 1e-8 * (1.0 + entry.abs()),
            "coefficient {j}: {entry} vs {numeric}"
        );
    }
}

#[test]
fn least_squares_recovers_the_generating_coefficients() {
    // y = 2 + 3 x with no noise; the normal equations recover it exactly.
    let design: Vec<Vec<f64>> = (1..=5).map(|x| vec![f64::from(x)]).collect();
    let y: Vec<f64> = (1..=5).map(|x| 2.0 + 3.0 * f64::from(x)).collect();
    let model = model(design, y);
    let (cost, gradient) = model.value_and_gradient(&[2.0, 3.0]).expect("value");
    assert!(cost < 1e-18, "cost {cost}");
    assert!(gradient.iter().all(|value| value.abs() < 1e-12));
}
