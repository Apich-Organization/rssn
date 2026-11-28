//! The shared weights preserve scale and exact-node behavior.

use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::Value;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::TensorValue;
use rssn::kernels::polynomial::Polynomial;

fn vector(values: &[f64]) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![values.len()],
        data: values.iter().copied().map(Value::Float).collect(),
        elem: TypeExpr::Float,
    })
}

#[test]
fn interpolation_retains_distinct_small_nodes_and_overflowing_spans() {
    let step = f64::from_bits(1);
    let small = Polynomial::interpolate_barycentric_eval(
        &[0.0, step, 4.0 * step],
        &[0.0, 1.0, 4.0],
        2.0 * step,
    )
    .unwrap();
    assert!((small - 2.0).abs() <= 3e-13, "{small}");
    let wide =
        Polynomial::interpolate_barycentric_eval(&[-f64::MAX, f64::MAX], &[2.0, 4.0], 0.0).unwrap();
    assert_eq!(wide.to_bits(), 3.0_f64.to_bits());
    assert!(Polynomial::interpolate_barycentric_eval(&[0.0, -0.0], &[1.0, 2.0], 0.0).is_err());
}

#[test]
fn polynomials_reuse_the_correct_weights_for_multiple_ordinates() {
    let nodes = [-1.0, 0.0, 1.0];
    for index in -15..=15 {
        let point = f64::from(index) / 15.0;
        for (ordinates, expected) in [
            ([2.0, 1.0, 2.0], 1.0 + point * point),
            ([2.0, 4.0, 6.0], 4.0 + 2.0 * point),
        ] {
            let actual =
                Polynomial::interpolate_barycentric_eval(&nodes, &ordinates, point).unwrap();
            assert!((actual - expected).abs() <= 1e-14);
        }
    }
}

#[test]
fn numerical_abi_uses_relative_gcd_and_bounded_power() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let Value::Tensor(monic) = session
        .call_value(
            "num_real_polynomial_gcd",
            &[
                vector(&[1e-200, -3e-200, 2e-200]),
                vector(&[1e200, 2e200, -3e200]),
                Value::Float(1e-12),
            ],
            &ctx,
        )
        .unwrap()
    else {
        panic!("expected a monic coefficient vector");
    };
    assert_eq!(monic.shape, [2]);
    let [Value::Float(leading), Value::Float(constant)] = &monic.data[..] else {
        panic!("expected numerical coefficients");
    };
    assert_eq!(leading.to_bits(), 1.0_f64.to_bits());
    assert!((constant + 1.0).abs() <= 1e-12);
    let Value::Tensor(power) = session
        .call_value(
            "num_real_polynomial_power",
            &[vector(&[1.0, -1.0]), Value::int(3)],
            &ctx,
        )
        .unwrap()
    else {
        panic!("expected a coefficient vector");
    };
    assert_eq!(
        power.data,
        [
            Value::Float(1.0),
            Value::Float(-3.0),
            Value::Float(3.0),
            Value::Float(-1.0)
        ]
    );
    assert!(
        session
            .call_value(
                "num_real_polynomial_power",
                &[vector(&[1.0, -1.0]), Value::int(4097)],
                &ctx
            )
            .is_err()
    );
}
