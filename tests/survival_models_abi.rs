//! Exact measured-time survival records and reproducible Cox likelihood behavior.

use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::Value;
use rssn::kernel::manifest::SignatureSource;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::TensorValue;

fn vector(
    values: Vec<Value>,
    elem: TypeExpr,
) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![values.len()],
        data: values,
        elem,
    })
}

fn time_values(
    values: &[f64],
    measured: bool,
) -> Value {
    let dim = DimExpr::base("Time");
    let wrap = |value| {
        if measured {
            Value::Dyn(DynValue {
                value: Box::new(Value::Float(value)),
                dim: dim.clone(),
            })
        } else {
            Value::Float(value)
        }
    };
    let elem = if measured {
        TypeExpr::Dim(Box::new(TypeExpr::Float), dim.clone())
    } else {
        TypeExpr::Float
    };
    vector(values.iter().copied().map(wrap).collect(), elem)
}

fn flags(values: &[bool]) -> Value {
    vector(
        values.iter().copied().map(Value::Bool).collect(),
        TypeExpr::Bool,
    )
}

fn matrix(
    columns: usize,
    values: &[f64],
) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![values.len() / columns, columns],
        data: values.iter().copied().map(Value::Float).collect(),
        elem: TypeExpr::Float,
    })
}

#[test]
fn ties_and_median_round_trip_with_measured_time() {
    let session = Session::new();
    let context = CallCtx::new();
    for name in [
        "num_survival_fit_cox",
        "num_survival_predict_cox",
        "num_survival_kaplan_meier",
        "num_survival_median",
    ] {
        assert_eq!(
            session.manifest().op(name).unwrap().signature_source,
            SignatureSource::Declared
        );
    }
    for measured in [false, true] {
        let estimate = session
            .call_value(
                "num_survival_kaplan_meier",
                &[
                    time_values(&[2.0, 3.0, 5.0, 6.0, 8.0, 9.0], measured),
                    flags(&[true, false, true, true, false, false]),
                    Value::Float(0.95),
                ],
                &context,
            )
            .unwrap();
        let result = session
            .call_value("num_survival_median", &[estimate], &context)
            .unwrap();
        let Value::Variant(result) = result else {
            panic!("expected optional median")
        };
        assert_eq!(result.tag, "Some");
        let scalar = &result.payload[0];
        if measured {
            let Value::Dyn(scalar) = scalar else {
                panic!("expected measured result")
            };
            assert_eq!(scalar.dim, DimExpr::base("Time"));
            assert_eq!(*scalar.value, Value::Float(6.0));
        } else {
            assert_eq!(*scalar, Value::Float(6.0));
        }
    }
    let estimate = session
        .call_value(
            "num_survival_kaplan_meier",
            &[
                time_values(&[1.0, 2.0, 3.0], true),
                flags(&[true, false, false]),
                Value::Float(0.95),
            ],
            &context,
        )
        .unwrap();
    let result = session
        .call_value("num_survival_median", &[estimate], &context)
        .unwrap();
    let Value::Variant(result) = result else {
        panic!("expected optional median")
    };
    assert_eq!(result.tag, "None");
}

#[test]
fn portable_cox_fit_can_predict_and_rejects_wrong_columns() {
    let session = Session::new();
    let context = CallCtx::new();
    let fitted = session
        .call_value(
            "num_survival_fit_cox",
            &[
                time_values(&[1.0, 2.0, 4.0, 7.0, 3.0, 5.0, 6.0, 8.0], true),
                flags(&[true, true, true, false, true, true, false, true]),
                matrix(1, &[1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0]),
                Value::String("efron".into()),
                Value::int(20),
                Value::Float(1e-10),
            ],
            &context,
        )
        .unwrap();
    assert!(
        session
            .call_value(
                "num_survival_predict_cox",
                &[fitted.clone(), matrix(2, &[0.0, 1.0])],
                &context
            )
            .is_err()
    );
    let predicted = session
        .call_value(
            "num_survival_predict_cox",
            &[fitted, matrix(1, &[0.0, 1.0])],
            &context,
        )
        .unwrap();
    let Value::Tensor(hazards) = predicted else {
        panic!("expected hazard tensor")
    };
    assert_eq!(hazards.shape, vec![2]);
    assert_eq!(hazards.data[0], Value::Float(1.0));
    assert!(hazards.data[1].as_f64().unwrap() > 1.0);
}
