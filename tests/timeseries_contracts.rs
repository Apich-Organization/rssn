//! Regression cases for differential reconstruction and Kalman conditioning.

use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::Value;
use rssn::kernel::overload::matches_type;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::TensorValue;
use rssn::kernels::timeseries;

fn close(
    actual: f64,
    expected: f64,
) {
    assert!((actual - expected).abs() < 2e-12, "{actual} vs {expected}");
}

#[test]
fn integrates_every_difference_order_and_the_empty_tail() {
    let original = [1.0, 2.0, 4.0, 8.0, 16.0, 32.0];
    for order in 0..=original.len() {
        let differences = timeseries::difference(&original, order);
        assert_eq!(
            timeseries::integrate(&differences, &original[..order]),
            original
        );
    }
    assert!(timeseries::integrate(&[1.0], &[f64::NAN]).is_empty());
}

#[test]
fn a_constant_state_observed_twice_matches_gaussian_conditioning() {
    let y = [vec![2.0], vec![4.0]];
    let identity = [vec![1.0]];
    let zero = [vec![0.0]];
    let result = timeseries::kalman_filter(
        &y,
        &identity,
        &identity,
        &zero,
        &identity,
        &[0.0],
        &identity,
    )
    .unwrap();
    close(result.filtered[0][0], 1.0);
    close(result.filtered[1][0], 2.0);
    close(result.covariance[0][0][0], 0.5);
    close(result.covariance[1][0][0], 1.0 / 3.0);
    let smooth = timeseries::kalman_smooth(
        &y,
        &identity,
        &identity,
        &zero,
        &identity,
        &[0.0],
        &identity,
    )
    .unwrap();
    close(smooth[0][0], 2.0);
    close(smooth[1][0], 2.0);
}

#[test]
fn likelihood_stays_finite_when_the_covariance_determinant_is_unrepresentable() {
    let identity = [vec![1.0, 0.0], vec![0.0, 1.0]];
    let zero = [vec![0.0, 0.0], vec![0.0, 0.0]];
    for variance in [1e-200, 1e200] {
        let covariance = [vec![variance, 0.0], vec![0.0, variance]];
        let result = timeseries::kalman_filter(
            &[vec![0.0, 0.0]],
            &identity,
            &identity,
            &zero,
            &covariance,
            &[0.0, 0.0],
            &covariance,
        )
        .unwrap();
        close(
            result.loglik,
            -((2.0 * std::f64::consts::PI).ln() + (2.0 * variance).ln()),
        );
        close(result.covariance[0][0][0] / variance, 0.5);
    }
}

#[test]
fn a_positive_determinant_alone_does_not_make_an_innovation_covariance_valid() {
    let mut identity = vec![vec![0.0; 4]; 4];
    for (i, row) in identity.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    let zero = vec![vec![0.0; 4]; 4];
    let noise = [
        vec![1.0, 3.0, 0.0, 0.0],
        vec![3.0, 1.0, 0.0, 0.0],
        vec![0.0, 0.0, 1.0, 3.0],
        vec![0.0, 0.0, 3.0, 1.0],
    ];
    let error = timeseries::kalman_filter(
        &[vec![0.0; 4]],
        &identity,
        &identity,
        &zero,
        &noise,
        &[0.0; 4],
        &zero,
    )
    .unwrap_err();
    assert!(error.contains("positive semidefinite"));
}

#[test]
fn rejects_an_indefinite_noise_covariance_even_with_positive_innovation() {
    let identity = [
        vec![1.0, 0.0, 0.0],
        vec![0.0, 1.0, 0.0],
        vec![0.0, 0.0, 1.0],
    ];
    let zero = [vec![0.0; 3], vec![0.0; 3], vec![0.0; 3]];
    let p0 = [
        vec![10.0, 0.0, 0.0],
        vec![0.0, 10.0, 0.0],
        vec![0.0, 0.0, 10.0],
    ];
    // R has an eigenvalue of -0.5, although its 2-by-2 principal
    // minors are positive and adding P0 makes the innovation SPD.
    let r = [
        vec![1.0, 0.75, 0.75],
        vec![0.75, 1.0, -0.75],
        vec![0.75, -0.75, 1.0],
    ];
    let error = timeseries::kalman_filter(
        &[vec![0.0; 3]],
        &identity,
        &identity,
        &zero,
        &r,
        &[0.0; 3],
        &p0,
    )
    .unwrap_err();
    assert!(error.contains("noise covariance must be positive semidefinite"));
}

#[test]
fn accepts_a_rank_one_model_covariance_with_unequal_component_scales() {
    let identity = [vec![1.0, 0.0], vec![0.0, 1.0]];
    let zero = [vec![0.0, 0.0], vec![0.0, 0.0]];
    let p0 = [vec![1e200, 1.0], vec![1.0, 1e-200]];
    let result = timeseries::kalman_filter(
        &[vec![0.0; 2]],
        &identity,
        &identity,
        &zero,
        &identity,
        &[0.0; 2],
        &p0,
    )
    .unwrap();
    assert!(result.loglik.is_finite());
    close(result.covariance[0][0][0], 1.0);
}

#[test]
fn finite_half_quadratic_likelihood_does_not_overflow_when_squaring() {
    let identity = [vec![1.0]];
    let zero = [vec![0.0]];
    let observed = 1.5e154;
    let result = timeseries::kalman_filter(
        &[vec![observed]],
        &identity,
        &identity,
        &zero,
        &identity,
        &[0.0],
        &zero,
    )
    .unwrap();
    let expected = -(0.5 * observed) * observed;
    close(result.loglik / expected, 1.0);
}

#[test]
fn scalar_posterior_variance_stays_accurate_when_noise_is_small_relative_to_the_prior() {
    let identity = [vec![1.0]];
    let zero = [vec![0.0]];
    let prior = [vec![3e200]];
    let result = timeseries::kalman_filter(
        &[vec![0.0]],
        &identity,
        &identity,
        &zero,
        &identity,
        &[0.0],
        &prior,
    )
    .unwrap();
    close(result.covariance[0][0][0], 1.0);
}

#[test]
fn rank_one_observation_preserves_tiny_or_zero_noise_variance() {
    let identity = [vec![1.0, 0.0], vec![0.0, 1.0]];
    let zero = [vec![0.0, 0.0], vec![0.0, 0.0]];
    let prior = [vec![2.0_f64.powi(1020), 0.0], vec![0.0, 1.0]];
    for variance in [0.0, 1e-200] {
        let result = timeseries::kalman_filter(
            &[vec![0.0]],
            &identity,
            &[vec![1.0, 0.0]],
            &zero,
            &[vec![variance]],
            &[0.0; 2],
            &prior,
        )
        .unwrap();
        let covariance = &result.covariance[0];
        if variance == 0.0 {
            assert_eq!(covariance[0][0], 0.0);
        } else {
            close(covariance[0][0] / variance, 1.0);
        }
        assert_eq!(covariance[0][1], 0.0);
        assert_eq!(covariance[1], [0.0, 1.0]);
        assert!(result.loglik.is_finite());
    }
}

#[test]
fn finite_scalar_correction_survives_an_unrepresentable_gain() {
    let identity = [vec![1.0]];
    let zero = [vec![0.0]];
    let prior = [vec![1e308]];
    let sensitivity = 1e-309;
    let observed = 1e-200;
    let result = timeseries::kalman_filter(
        &[vec![observed]],
        &identity,
        &[vec![sensitivity]],
        &zero,
        &zero,
        &[0.0],
        &prior,
    )
    .unwrap();
    close(result.filtered[0][0] / (observed / sensitivity), 1.0);
    assert_eq!(result.covariance[0][0][0], 0.0);
}

fn measured(
    shape: &[usize],
    data: &[f64],
    dimension: DimExpr,
) -> Value {
    Value::Tensor(TensorValue {
        shape: shape.to_vec(),
        data: data.iter().copied().map(Value::Float).collect(),
        elem: TypeExpr::Float.with_dim(dimension),
    })
}

#[test]
fn measured_filter_preserves_state_units_and_uses_the_likelihood_reference() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let position = DimExpr::base("Length");
    let signal = DimExpr::base("Mass");
    let squared_position = position.clone().pow(2);
    let sensor = DimExpr::Mul(vec![signal.clone(), position.clone().pow(-1)]);
    let mut args = [
        measured(&[2, 1], &[2.0, 4.0], signal.clone()),
        measured(&[1, 1], &[1.0], DimExpr::Dimless),
        measured(&[1, 1], &[1.0], sensor),
        measured(&[1, 1], &[0.0], squared_position.clone()),
        measured(&[1, 1], &[1.0], signal.clone().pow(2)),
        measured(&[1], &[0.0], position.clone()),
        measured(&[1, 1], &[1.0], squared_position.clone()),
        Value::Dyn(DynValue {
            value: Box::new(Value::Float(1.0)),
            dim: signal.clone(),
        }),
    ];
    let call = |args: &[Value]| {
        session
            .call_value("num_timeseries_kalman_filter", args, &ctx)
            .unwrap()
    };
    let Value::Record(first) = call(&args) else {
        panic!("expected Kalman history")
    };
    let Some(Value::Tensor(cov)) = first.get("covariance") else {
        panic!("expected covariance history")
    };
    assert_eq!(cov.shape, [2, 1, 1]);
    assert!(matches_type(&TypeExpr::Float.with_dim(squared_position), &cov.elem).accepted());
    let Some(Value::Float(first_loglik)) = first.get("loglik") else {
        panic!("expected a scalar likelihood")
    };

    args[7] = Value::Dyn(DynValue {
        value: Box::new(Value::Float(10.0)),
        dim: signal,
    });
    let Value::Record(second) = call(&args) else {
        panic!("expected Kalman history")
    };
    assert_eq!(first.get("filtered"), second.get("filtered"));
    let Some(Value::Float(second_loglik)) = second.get("loglik") else {
        panic!("expected a scalar likelihood")
    };
    close(second_loglik - first_loglik, 2.0 * 10.0_f64.ln());

    args[7] = Value::Dyn(DynValue {
        value: Box::new(Value::Float(10.0)),
        dim: position,
    });
    assert!(
        session
            .call_value("num_timeseries_kalman_filter", &args, &ctx)
            .is_err()
    );
}

#[test]
fn symbolic_kalman_filter_accepts_a_rectangular_observation_matrix() {
    let session = Session::new();
    let request = session
        .parse(concat!(
            "kalman_filter(list(list(2.0), list(4.0)), ",
            "list(list(1.0, 0.0), list(0.0, 1.0)), list(list(1.0, 0.0)), ",
            "list(list(0.0, 0.0), list(0.0, 0.0)), list(list(1.0)))"
        ))
        .unwrap();
    let result = session
        .compute(request, &rssn::api::Config::default())
        .unwrap();
    let text = result.term.to_string();
    assert!(
        text.starts_with("list(") && !text.contains("kalman_filter"),
        "{text}"
    );
}

#[test]
fn rolling_moments_preserve_units_of_means_and_square_the_units_of_variances() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let length = DimExpr::base("Length");
    let samples = measured(&[3], &[1.0, 3.0, 5.0], length.clone());
    for width in [2, 4] {
        let args = [samples.clone(), Value::int(width), Value::int(0)];
        let Value::Record(result) = session
            .call_value("num_timeseries_rolling_moments", &args, &ctx)
            .unwrap()
        else {
            panic!("expected rolling moments");
        };
        for (field, dim, numbers) in [
            ("mean", length.clone(), [2.0, 4.0]),
            ("variance", length.clone().pow(2), [1.0, 1.0]),
        ] {
            let Some(Value::Tensor(tensor)) = result.get(field) else {
                panic!("expected a tensor for {field}");
            };
            let expected = TypeExpr::Float.with_dim(dim);
            assert!(matches_type(&expected, &tensor.elem).accepted());
            if width == 2 {
                assert_eq!(tensor.shape, [2]);
                assert_eq!(tensor.data.len(), numbers.len());
                for (element, number) in tensor.data.iter().zip(numbers) {
                    assert!(matches_type(&expected, &element.type_expr()).accepted());
                    let Value::Dyn(measured) = element else {
                        panic!("expected a measured value for {field}");
                    };
                    assert_eq!(*measured.value, Value::Float(number));
                }
            } else {
                assert_eq!(tensor.shape, [0]);
                assert!(tensor.data.is_empty());
            }
        }
    }
}
