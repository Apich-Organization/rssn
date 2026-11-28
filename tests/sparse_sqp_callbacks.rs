//! Sparse SQP callback parity, work budgets and independent projections.

use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::NativeCallback;
use rssn::kernel::Value;
use rssn::kernel::signature::NatExpr;
use rssn::kernel::signature::Signature;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::RecordValue;
use rssn::kernel::value::TensorValue;
use std::sync::Arc;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering::Relaxed;

fn vector_type(elem: TypeExpr) -> TypeExpr {
    TypeExpr::Tensor {
        rank: NatExpr::Lit(1),
        shape: None,
        elem: Box::new(elem),
    }
}
fn record_type(fields: &[(&str, TypeExpr)]) -> TypeExpr {
    TypeExpr::Record(
        fields
            .iter()
            .map(|(name, ty)| ((*name).into(), ty.clone()))
            .collect(),
    )
}
fn record(fields: &[(&str, Value)]) -> Value {
    Value::Record(RecordValue {
        fields: fields
            .iter()
            .map(|(name, value)| ((*name).into(), value.clone()))
            .collect(),
    })
}
fn field<'a>(
    value: &'a Value,
    name: &str,
) -> &'a Value {
    let Value::Record(record) = value else {
        panic!("record")
    };
    record.get(name).unwrap()
}
fn tensor(
    data: Vec<Value>,
    elem: TypeExpr,
) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![data.len()],
        data,
        elem,
    })
}
fn floats(values: &[f64]) -> Value {
    tensor(
        values.iter().copied().map(Value::Float).collect(),
        TypeExpr::Float,
    )
}
fn indices(values: &[usize]) -> Value {
    tensor(
        values.iter().copied().map(Value::int).collect(),
        TypeExpr::Int,
    )
}
fn sparse(
    rows: &[usize],
    columns: &[usize],
    values: &[f64],
) -> Value {
    record(&[
        ("rows", indices(rows)),
        ("columns", indices(columns)),
        ("values", floats(values)),
    ])
}
fn sparse_type() -> TypeExpr {
    record_type(&[
        ("rows", vector_type(TypeExpr::Int)),
        ("columns", vector_type(TypeExpr::Int)),
        ("values", vector_type(TypeExpr::Float)),
    ])
}

#[test]
fn solver_shares_one_checked_evaluation_per_point_and_reports_callback_failures() {
    let mut ctx = CallCtx::new();
    let session = Session::new();
    let calls = Arc::new(AtomicUsize::new(0));
    let checked = Arc::clone(&calls);
    let vector = vector_type(TypeExpr::Float);
    let evaluate = ctx
        .add_callback(Arc::new(NativeCallback::new(
            move |args: &[Value]| {
                checked.fetch_add(1, Relaxed);
                let Value::Tensor(x) = &args[0] else {
                    panic!("vector expected")
                };
                let x = x.data[0].as_f64().unwrap();
                let d = x - 3.0;
                Ok(record(&[
                    ("objective", Value::Float(1e9 + d * d)),
                    ("gradient", floats(&[2.0 * d])),
                    ("constraints", floats(&[x - 1.0])),
                    ("jacobian", sparse(&[0], &[0], &[1.0])),
                ]))
            },
            Signature::new(
                vec![vector.clone()],
                record_type(&[
                    ("objective", TypeExpr::Float),
                    ("gradient", vector.clone()),
                    ("constraints", vector.clone()),
                    ("jacobian", sparse_type()),
                ]),
            ),
        )))
        .unwrap();
    let hessian = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |_| Ok(sparse(&[0], &[0], &[2.0])),
            Signature::new(
                vec![vector.clone(), vector.clone(), vector.clone()],
                sparse_type(),
            ),
        )))
        .unwrap();
    let broken = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |_| Ok(sparse(&[1], &[0], &[2.0])),
            Signature::new(vec![vector.clone(), vector.clone(), vector], sparse_type()),
        )))
        .unwrap();
    let shape = record(&[
        ("equalities", Value::int(0)),
        ("inequalities", Value::int(1)),
    ]);
    let args = |budget, h| {
        [
            floats(&[0.0]),
            shape.clone(),
            Value::Float(1e-10),
            Value::Float(10.0),
            Value::int(50),
            Value::int(budget),
            Value::Function(evaluate),
            Value::Function(h),
        ]
    };
    let result = session
        .call_value("num_sqp_solve", &args(100, hessian), &ctx)
        .unwrap();
    assert_eq!(
        *field(&result, "converged"),
        Value::Bool(true),
        "{result:?}"
    );
    let Value::Tensor(x) = field(&result, "x") else {
        panic!("solution")
    };
    assert!((x.data[0].as_f64().unwrap() - 1.0).abs() < 1e-10);
    assert!(calls.load(Relaxed) <= 3);
    let mut stationary = args(100, hessian);
    stationary[0] = floats(&[1.0]);
    let stationary = session
        .call_value("num_sqp_solve", &stationary, &ctx)
        .unwrap();
    assert_eq!(*field(&stationary, "converged"), Value::Bool(true));
    let Value::Tensor(lambda) = field(&stationary, "inequality_multipliers") else {
        panic!("multipliers")
    };
    assert!((lambda.data[0].as_f64().unwrap() - 4.0).abs() < 1e-12);
    assert!(
        session
            .call_value("num_sqp_solve", &args(1, hessian), &ctx)
            .unwrap_err()
            .to_string()
            .contains("budget")
    );
    assert!(
        session
            .call_value("num_sqp_solve", &args(100, broken), &ctx)
            .unwrap_err()
            .to_string()
            .contains("triplet")
    );
    assert_eq!(
        session
            .call_value(
                "num_sqp_oracle_equality_count",
                std::slice::from_ref(&shape),
                &ctx,
            )
            .unwrap(),
        Value::int(0)
    );
    assert_eq!(
        session
            .call_value(
                "num_sqp_oracle_inequality_count",
                std::slice::from_ref(&shape),
                &ctx,
            )
            .unwrap(),
        Value::int(1)
    );
    let query = [floats(&[0.0]), shape, Value::Function(evaluate)];
    let result = session
        .call_value("num_sqp_oracle_gradient", &query, &ctx)
        .unwrap();
    let Value::Tensor(t) = result else {
        panic!("gradient")
    };
    assert_eq!(t.data, vec![Value::Float(-6.0)]);
}
