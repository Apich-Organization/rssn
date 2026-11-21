//! Real compiled backends preserve retained sessions, physical references and checked lifetimes.
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::RepExpr;
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
fn native_column_gather_keeps_measured_values_and_checks_row_contracts() {
    let session = Session::new();
    let length = DimExpr::base("Length");
    let reference = measured(2.0, length.clone());
    let constructor = "num_native_tape_column_gather";
    let tape = call(&session, constructor, &[Value::int(2), reference.clone()]).unwrap();
    for invalid in [0, 63] {
        assert!(
            call(
                &session,
                constructor,
                &[Value::int(invalid), reference.clone()]
            )
            .is_err()
        );
    }
    let dynamic = |cells| {
        vector(
            cells,
            TypeExpr::Dyn(RepExpr {
                repr: "term".into(),
                dim: None,
            }),
        )
    };
    let distance = |m| measured(m, length.clone());
    let distance_handle = scalar(&session, 18.0, length.clone(), false);
    let source = dynamic(vec![distance(2.0), distance(4.0), distance(6.0)]);
    let other = dynamic(vec![
        distance(14.0),
        distance_handle.clone(),
        distance(22.0),
    ]);
    let offsets = |values: [f64; 3]| vector(values.map(Value::Float).into(), TypeExpr::Float);
    let columns = dynamic(vec![
        offsets([0.0, 1.9, 0.0]),
        source.clone(),
        other.clone(),
    ]);
    let eval = |session: &Session, rows, columns| {
        call(
            session,
            "num_native_tape_eval_columns",
            &[tape.clone(), Value::int(rows), columns, dynamic(Vec::new())],
        )
    };
    let Value::Tensor(outputs) = eval(&session, 3, columns.clone()).unwrap() else {
        panic!("expected output columns");
    };
    let [Value::Tensor(selected)] = outputs.data.as_slice() else {
        panic!("expected one selected column");
    };
    assert_eq!(selected.shape, [3]);
    assert_eq!(selected.elem, TypeExpr::Float.with_dim(length));
    assert_eq!(
        selected.data.iter().map(numeric).collect::<Vec<_>>(),
        [2.0, 18.0, 6.0]
    );
    let typed_source = vector(
        vec![Value::Float(2.0), Value::Float(4.0), Value::Float(6.0)],
        selected.elem.clone(),
    );
    let typed_other = vector(
        vec![Value::Float(14.0), distance_handle, Value::Float(22.0)],
        selected.elem.clone(),
    );
    let typed_columns = dynamic(vec![
        offsets([0.0, 1.9, 0.0]),
        typed_source.clone(),
        typed_other,
    ]);
    assert_eq!(
        eval(&session, 3, typed_columns).unwrap(),
        Value::Tensor(outputs.clone())
    );
    let inconsistent = vector(
        vec![measured(1.0, DimExpr::base("Time")); 3],
        selected.elem.clone(),
    );
    assert!(
        eval(
            &session,
            3,
            dynamic(vec![offsets([0.0; 3]), typed_source, inconsistent])
        )
        .is_err()
    );
    assert!(eval(&Session::new(), 3, columns.clone()).is_err());
    assert!(eval(&session, 2, columns).is_err());
    let bad_offset = dynamic(vec![offsets([0.0, 2.0, 1.0]), source.clone(), other]);
    assert!(eval(&session, 3, bad_offset).is_err());
    let wrong_unit = dynamic(vec![measured(1.0, DimExpr::base("Time")); 3]);
    let bad_source = dynamic(vec![offsets([0.0, 0.0, 0.0]), source, wrong_unit]);
    assert!(eval(&session, 3, bad_source).is_err());
    let empty = dynamic(vec![dynamic(Vec::new()); 3]);
    let Value::Tensor(outputs) = eval(&session, 0, empty).unwrap() else {
        panic!("expected empty measured output column");
    };
    assert!(
        matches!(&outputs.data[..], [Value::Tensor(t)] if t.shape == [0] && t.elem == selected.elem)
    );
}

#[test]
fn native_compiled_point_and_parallel_columns_use_real_normalized_program() {
    let session = Session::new();
    let program = compile(&session, 1_000_000);
    assert_eq!(
        call(
            &session,
            "num_backend_arity",
            std::slice::from_ref(&program)
        )
        .unwrap(),
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
    let length_with_cancelled_time = DimExpr::Mul(vec![
        DimExpr::Inverse(Box::new(DimExpr::base("Time"))),
        DimExpr::base("Length"),
        DimExpr::base("Time"),
    ]);
    let inputs = [x, scalar(&session, 6.0, length_with_cancelled_time, false)];
    let mut data = (0..rows)
        .map(|row| inputs[row % inputs.len()].clone())
        .collect::<Vec<_>>();
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
    for (row, value) in tensor.data.iter().enumerate() {
        assert_eq!(numeric(value).to_bits(), [2.0f64, 3.0][row % 2].to_bits());
    }

    let manifest = session.manifest();
    let operation = manifest
        .ops
        .iter()
        .find(|op| op.name == "num_backend_call_batch_par")
        .unwrap();
    assert!(operation.doc.contains("estimates the work per row"));
    assert!(
        operation
            .doc
            .contains("Cancellation is checked between blocks")
    );
}

#[test]
fn measured_compiled_columns_match_handle_batches_and_check_inputs() {
    let session = Session::new();
    let program = compile(&session, 1_000_000);
    let rows = 9000;
    let length = DimExpr::base("Length");
    let time = DimExpr::base("Time");
    let dyn_type = TypeExpr::Dyn(RepExpr {
        repr: "term".into(),
        dim: None,
    });
    let physical = |data, dim| vector(data, TypeExpr::Float.with_dim(dim));
    let x_column = physical(
        (0..rows).map(|i| Value::Float((i + 1) as f64)).collect(),
        length.clone(),
    );
    let t_column = physical(vec![Value::Float(2.0); rows], time.clone());
    let eval = |session: &Session, row_count, xs, ts| {
        call(
            session,
            "num_backend_call_batch_columns",
            &[
                program.clone(),
                Value::int(row_count),
                vector(vec![xs, ts], dyn_type.clone()),
            ],
        )
    };
    let Value::Tensor(result) = eval(&session, rows, x_column.clone(), t_column.clone()).unwrap()
    else {
        panic!("expected measured vector");
    };
    assert!(
        rssn::kernel::overload::matches_type(
            &TypeExpr::Float.with_dim(DimExpr::Mul(vec![length.clone(), time.clone().pow(-1)])),
            &result.elem
        )
        .accepted()
    );
    assert_eq!(result.shape, [rows]);
    for (i, v) in result.data.iter().enumerate() {
        assert_eq!(numeric(v), (i + 1) as f64 / 2.0);
    }

    let first_x = scalar(&session, 1.0, length.clone(), false);
    let first_t = scalar(&session, 2.0, time.clone(), false);
    let handles = Value::Tensor(TensorValue {
        shape: vec![2, 1],
        elem: TypeExpr::Opaque("rssn/native/backend-scalar/v1".into()),
        data: vec![first_x, first_t],
    });
    let first = call(
        &session,
        "num_backend_call_batch",
        &[program.clone(), handles],
    )
    .unwrap();
    assert_eq!(
        first,
        eval(
            &session,
            1,
            physical(vec![Value::Float(1.0)], length.clone()),
            physical(vec![Value::Float(2.0)], time.clone())
        )
        .unwrap()
    );
    assert!(eval(&session, rows - 1, x_column.clone(), t_column.clone()).is_err());
    assert!(eval(&session, rows, t_column.clone(), x_column.clone()).is_err());
    assert!(eval(&Session::new(), rows, x_column.clone(), t_column.clone()).is_err());
    let wrong_cell = physical(vec![measured(1.0, time.clone())], length.clone());
    assert!(
        eval(
            &session,
            1,
            wrong_cell,
            physical(vec![Value::Float(2.0)], time.clone())
        )
        .is_err()
    );
    let empty = eval(
        &session,
        0,
        physical(vec![], length),
        physical(vec![], time),
    )
    .unwrap();
    assert!(matches!(empty, Value::Tensor(t) if t.shape == [0] && t.elem == result.elem));
    let token = std::sync::Arc::new(std::sync::atomic::AtomicBool::new(true));
    let _scope = rssn::cancellation::install(token);
    assert!(matches!(
        eval(&session, rows, x_column, t_column),
        Err(OpError::HostDiagnostic { code: "E0812", .. })
    ));
}

#[test]
fn measured_columns_of_a_constant_program_use_the_explicit_row_extent() {
    let session = Session::new();
    let reference = scalar(&session, 1.0, DimExpr::Dimless, true);
    let program = call(
        &session,
        "num_backend_compile_current",
        &[
            Value::String("7".into()),
            vector(vec![], TypeExpr::String),
            vector(
                vec![],
                TypeExpr::Opaque("rssn/native/backend-scalar/v1".into()),
            ),
            reference,
            Value::int(64),
        ],
    )
    .unwrap();
    let cols = vector(
        vec![],
        TypeExpr::Dyn(RepExpr {
            repr: "term".into(),
            dim: None,
        }),
    );
    for rows in [0, 3] {
        assert_eq!(
            call(
                &session,
                "num_backend_call_batch_columns",
                &[program.clone(), Value::int(rows), cols.clone()]
            )
            .unwrap(),
            vector(
                vec![Value::Float(7.0); rows],
                TypeExpr::Float.with_dim(DimExpr::Dimless),
            ),
        );
    }
    assert!(
        call(
            &session,
            "num_backend_call_batch_columns",
            &[program, Value::int(16_385), cols]
        )
        .is_err()
    );
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
    assert!(
        call(
            &Session::new(),
            "num_backend_arity",
            std::slice::from_ref(&program)
        )
        .is_err()
    );
    assert!(call(&session, "num_backend_reference", &[Value::Float(0.0)]).is_err());
    let mut forged = program.clone();
    if let Value::Opaque(v) = &mut forged {
        v.kind = "rssn/native/backend/v1".into();
    }
    assert!(call(&session, "num_backend_arity", &[forged]).is_err());
    call(
        &session,
        "num_backend_release",
        std::slice::from_ref(&program),
    )
    .unwrap();
    assert!(call(&session, "num_backend_arity", &[program]).is_err());
}
