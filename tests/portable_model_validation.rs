//! Host-provided descriptors and extreme physical scales cross the public ABI.
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::Exp;
use rssn::kernel::signature::Rational;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::OpaqueValue;
use rssn::kernel::value::RecordValue;
use rssn::kernel::value::TensorValue;

fn scalar(
    x: f64,
    dimension: &str,
) -> Value {
    Value::Dyn(DynValue {
        value: Box::new(Value::Float(x)),
        dim: DimExpr::Base(dimension.into()),
    })
}
fn tensor(
    shape: &[usize],
    values: &[f64],
    dimension: DimExpr,
) -> Value {
    Value::Tensor(TensorValue {
        shape: shape.to_vec(),
        data: values.iter().copied().map(Value::Float).collect(),
        elem: TypeExpr::Dim(Box::new(TypeExpr::Float), dimension),
    })
}
fn number(value: &Value) -> f64 {
    match value {
        | Value::Float(v) => *v,
        | Value::Dyn(v) => number(&v.value),
        | other => panic!("expected numeric result: {other:?}"),
    }
}
fn call(
    name: &str,
    arguments: &[Value],
) -> Value {
    Session::new()
        .call_value(name, arguments, &CallCtx::new())
        .unwrap_or_else(|error| panic!("{name}: {error:?}"))
}

#[test]
fn valid_portable_encoding_cannot_forge_a_panicking_spline() {
    // The envelope is valid and finite, but degree == coefficient count makes
    // the spline mathematically invalid. Validate this before de Boor indexing.
    let rows = [vec![0.0_f64, 0.0, 1.0], vec![1.0], vec![1.0]];
    let mut bytes = (rows.len() as u64).to_le_bytes().to_vec();
    for row in rows {
        bytes.extend_from_slice(&(row.len() as u64).to_le_bytes());
        for value in row {
            bytes.extend_from_slice(&value.to_le_bytes());
        }
    }
    let model = Value::Record(RecordValue {
        fields: vec![
            (
                "payload".into(),
                Value::Opaque(OpaqueValue {
                    kind: "rssn/bspline-defining-data/v1".into(),
                    bytes,
                }),
            ),
            ("coordinate_unit".into(), scalar(1.0, "Time")),
            ("value_unit".into(), scalar(1.0, "Length")),
        ],
    });
    assert!(matches!(
        Session::new().call_value(
            "num_bspline_query",
            &[model, scalar(0.5, "Time")],
            &CallCtx::new()
        ),
        Err(OpError::Domain(_))
    ));
}

#[test]
fn scaled_hessian_preserves_representable_result_when_scale_square_overflows() {
    let result = call(
        "num_ad_hessian_scaled",
        &[
            Value::String("x*x".into()),
            Value::Tensor(TensorValue {
                shape: vec![1],
                data: vec![Value::String("x".into())],
                elem: TypeExpr::String,
            }),
            tensor(&[1], &[1e200], DimExpr::Base("Length".into())),
            scalar(1e200, "Length"),
            scalar(1e300, "Time"),
        ],
    );
    let Value::Tensor(matrix) = result else {
        panic!("matrix expected")
    };
    let actual = number(&matrix.data[0]);
    assert!((actual / 2e-100 - 1.0).abs() < 1e-14, "{actual}");
}

#[test]
fn scaled_log_density_avoids_false_singular_covariance_from_scale_square_overflow() {
    let result = call(
        "num_multivariate_normal_log_pdf_scaled",
        &[
            tensor(&[1], &[0.0], DimExpr::Base("Length".into())),
            tensor(
                &[1, 1],
                &[1e300],
                DimExpr::Base("Length".into()).pow(Exp::Rat(Rational::new(2, 1))),
            ),
            tensor(&[1], &[0.0], DimExpr::Base("Length".into())),
            scalar(1e200, "Length"),
        ],
    );
    let expected = -0.5 * ((2.0 * std::f64::consts::PI).ln() + 1e-100_f64.ln());
    assert!((number(&result) - expected).abs() < 1e-12);
}
