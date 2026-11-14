//! Real compiled backends preserve retained sessions, physical references and checked lifetimes.
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::TensorValue;
fn call(
    session: &Session,
    name: &str,
    args: &[Value],
) -> Result<Value, OpError> {
    session.call_value(name, args, &CallCtx::new())
}
fn vector(
    values: Vec<Value>,
    elem: TypeExpr,
) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![values.len()],
        elem,
        data: values,
    })
}
fn measured(
    value: f64,
    dimension: DimExpr,
) -> Value {
    Value::Dyn(DynValue {
        value: Box::new(Value::Float(value)),
        dim: dimension,
    })
}
fn scalar(
    session: &Session,
    value: f64,
    dimension: DimExpr,
    reference: bool,
) -> Value {
    call(
        session,
        if reference {
            "num_backend_reference"
        } else {
            "num_backend_value"
        },
        &[measured(value, dimension)],
    )
    .unwrap()
}
fn compile(
    session: &Session,
    budget: usize,
) -> Value {
    let backend = call(session, "num_backend_interpreter", &[]).unwrap();
    let x = scalar(session, 2.0, DimExpr::Base("Length".into()), true);
    let t = scalar(session, 1.0, DimExpr::Base("Time".into()), true);
    let output = scalar(
        session,
        2.0,
        DimExpr::Mul(vec![
            DimExpr::Base("Length".into()),
            DimExpr::Inverse(Box::new(DimExpr::Base("Time".into()))),
        ]),
        true,
    );
    call(
        session,
        "num_backend_compile",
        &[
            backend,
            Value::String("x/t".into()),
            vector(
                vec![Value::String("x".into()), Value::String("t".into())],
                TypeExpr::String,
            ),
            vector(
                vec![x, t],
                TypeExpr::Opaque("rssn/native/backend-scalar/v1".into()),
            ),
            output,
            Value::int(budget),
        ],
    )
    .unwrap()
}
fn numeric(value: &Value) -> f64 {
    match value {
        | Value::Dyn(d) => numeric(&d.value),
        | _ => value.as_f64().unwrap(),
    }
}
#[test]
fn native_compiled_point_and_parallel_columns_use_real_normalized_program() {
    let session = Session::new();
    let program = compile(&session, 1_000_000);
    assert_eq!(
        call(&session, "num_backend_arity", &[program.clone()]).unwrap(),
        Value::int(2)
    );
    let x = scalar(&session, 4.0, DimExpr::Base("Length".into()), false);
    let t = scalar(&session, 2.0, DimExpr::Base("Time".into()), false);
    let point = vector(
        vec![x.clone(), t.clone()],
        TypeExpr::Opaque("rssn/native/backend-scalar/v1".into()),
    );
    let value = call(&session, "num_backend_call", &[program.clone(), point]).unwrap();
    assert_eq!(numeric(&value), 2.0);
    assert!(matches!(value, Value::Dyn(_)));
    let rows = 9000;
    let mut data = vec![x; rows];
    data.extend(vec![t; rows]);
    let columns = Value::Tensor(TensorValue {
        shape: vec![2, rows],
        elem: TypeExpr::Opaque("rssn/native/backend-scalar/v1".into()),
        data,
    });
    let sequential = call(
        &session,
        "num_backend_call_batch",
        &[program.clone(), columns.clone()],
    )
    .unwrap();
    let parallel = call(&session, "num_backend_call_batch_par", &[program, columns]).unwrap();
    assert_eq!(sequential, parallel);
    let Value::Tensor(tensor) = parallel else {
        panic!("tensor")
    };
    assert_eq!(tensor.shape, vec![rows]);
    assert!(tensor.data.iter().all(|v| numeric(v) == 2.0));
}
#[test]
fn native_compiled_rejects_wrong_units_budgets_stale_and_foreign_handles() {
    let session = Session::new();
    let program = compile(&session, 8);
    let wrong = scalar(&session, 4.0, DimExpr::Base("Time".into()), false);
    let time = scalar(&session, 2.0, DimExpr::Base("Time".into()), false);
    assert!(
        call(
            &session,
            "num_backend_call",
            &[
                program.clone(),
                vector(
                    vec![wrong.clone(), time.clone()],
                    TypeExpr::Opaque("rssn/native/backend-scalar/v1".into())
                )
            ]
        )
        .is_err()
    );
    let columns = Value::Tensor(TensorValue {
        shape: vec![2, 10],
        elem: TypeExpr::Opaque("rssn/native/backend-scalar/v1".into()),
        data: vec![time; 20],
    });
    assert!(
        call(
            &session,
            "num_backend_call_batch",
            &[program.clone(), columns]
        )
        .is_err()
    );
    assert!(call(&Session::new(), "num_backend_arity", &[program.clone()]).is_err());
    assert!(call(&session, "num_backend_reference", &[Value::Float(0.0)]).is_err());
    let mut forged = program.clone();
    if let Value::Opaque(v) = &mut forged {
        v.kind = "rssn/native/backend/v1".into();
    }
    assert!(call(&session, "num_backend_arity", &[forged]).is_err());
    call(&session, "num_backend_release", &[program.clone()]).unwrap();
    assert!(call(&session, "num_backend_arity", &[program]).is_err());
}
