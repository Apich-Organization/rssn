//! Scaled ordinary least-squares regression kernels and ABI.

use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::Exp;
use rssn::kernel::signature::Rational;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::TensorValue;
use rssn::kernels::stats::simple_linear_regression;

#[test]
fn linear_fit_preserves_slope_across_extreme_and_adjacent_samples() {
    for scale in [1e-200, 1.0, 1e200] {
        let (slope, intercept) =
            simple_linear_regression(&[(-scale, -2.0 * scale), (scale, 2.0 * scale)]);
        assert!((slope - 2.0).abs() <= f64::EPSILON);
        assert_eq!(intercept, 0.0);
    }

    let tiny = f64::from_bits(1);
    let (slope, intercept) = simple_linear_regression(&[(-1.0, -tiny), (1.0, tiny)]);
    assert_eq!(slope.to_bits(), tiny.to_bits());
    assert_eq!(intercept, 0.0);

    let (slope, intercept) = simple_linear_regression(&[(1.0, 0.0), (1.0_f64.next_up(), 1.0)]);
    assert_eq!(slope, 4_503_599_627_370_496.0);
    assert_eq!(intercept, -slope);

    let outer = 2.0_f64.powi(540);
    let inner = 2.0_f64.powi(-560);
    let response = 2.0_f64.powi(950);
    let (slope, intercept) = simple_linear_regression(&[
        (-outer, 0.0),
        (outer, 0.0),
        (-inner, -response),
        (inner, response),
    ]);
    assert_eq!(slope.to_bits(), 2.0_f64.powi(-690).to_bits());
    assert_eq!(intercept, 0.0);
}

#[test]
fn constant_and_non_finite_inputs_keep_documented_distinct_results() {
    let (slope, intercept) = simple_linear_regression(&[(f64::MAX, f64::MAX); 3]);
    assert!(slope.is_nan());
    assert_eq!(intercept, f64::MAX);
    for data in [vec![], vec![(f64::NAN, 1.0)], vec![(1.0, f64::INFINITY)]] {
        let (slope, intercept) = simple_linear_regression(&data);
        assert!(slope.is_nan() && intercept.is_nan());
    }
}

#[test]
fn two_measured_columns_keep_independent_slope_and_intercept_dimensions() {
    let session = Session::new();
    let time = DimExpr::Base("Time".into());
    let length = DimExpr::Base("Length".into());
    let vector = |dimension: DimExpr, samples: [f64; 3]| {
        Value::Tensor(TensorValue {
            shape: vec![samples.len()],
            data: samples
                .into_iter()
                .map(|sample| {
                    Value::Dyn(DynValue {
                        dim: dimension.clone(),
                        value: Box::new(Value::Float(sample)),
                    })
                })
                .collect(),
            elem: TypeExpr::Float.with_dim(dimension),
        })
    };
    let result = session
        .call_value(
            "num_stats_simple_linear_regression",
            &[
                vector(time.clone(), [-1e200, 0.0, 1e200]),
                vector(length.clone(), [-1e200, 1e200, 3e200]),
            ],
            &CallCtx::new(),
        )
        .unwrap();
    let Value::Record(record) = result else {
        panic!("expected measured slope and intercept");
    };
    for (name, expected, dimension) in [
        (
            "item0",
            2.0,
            DimExpr::Mul(vec![length.clone(), time.pow(Exp::Rat(Rational::int(-1)))]),
        ),
        ("item1", 1e200, length),
    ] {
        let Some(Value::Dyn(measured)) = record.get(name) else {
            panic!("expected a measured {name}");
        };
        assert_eq!(measured.dim, dimension);
        let actual = measured.value.as_f64().unwrap();
        assert!((actual / expected - 1.0).abs() <= 4.0 * f64::EPSILON);
    }
}
