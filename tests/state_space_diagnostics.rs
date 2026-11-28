//! Checked discrete control diagnostics and their physical dimensions.

use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::overload::matches_type;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::RecordValue;
use rssn::kernel::value::TensorValue;

fn matrix(
    rows: usize,
    cols: usize,
    dim: DimExpr,
    samples: &[f64],
) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![rows, cols],
        elem: TypeExpr::Float.with_dim(dim),
        data: samples.iter().copied().map(Value::Float).collect(),
    })
}

fn model(
    a: Value,
    b: Value,
    c: Value,
    d: Value,
) -> Value {
    Value::Record(RecordValue {
        fields: vec![
            ("a".into(), a),
            ("b".into(), b),
            ("c".into(), c),
            ("d".into(), d),
            (
                "dt".into(),
                Value::Dyn(DynValue {
                    value: Box::new(Value::Float(0.25)),
                    dim: DimExpr::base("Time"),
                }),
            ),
        ],
    })
}

fn units() -> (DimExpr, DimExpr, DimExpr) {
    let b = DimExpr::base("Length");
    let c = DimExpr::base("Mass");
    let d = DimExpr::Mul(vec![b.clone(), c.clone()]);
    (b, c, d)
}

fn real(value: &Value) -> f64 {
    match value {
        | Value::Dyn(value) => real(&value.value),
        | Value::Float(value) => *value,
        | _ => panic!("expected a real value"),
    }
}

fn entries(
    value: &Value,
    shape: &[usize],
    dim: DimExpr,
) -> Vec<f64> {
    let Value::Tensor(array) = value else {
        panic!("expected a matrix")
    };
    assert_eq!(array.shape, shape);
    assert!(matches_type(&TypeExpr::Float.with_dim(dim), &array.elem).accepted());
    array.data.iter().map(real).collect()
}

fn call(
    session: &Session,
    operator: &str,
    model: &Value,
) -> Result<Value, OpError> {
    session.call_value(operator, std::slice::from_ref(model), &CallCtx::new())
}

fn diagonal(first: f64) -> Value {
    let (b, c, d) = units();
    model(
        matrix(2, 2, DimExpr::Dimless, &[first, 0.0, 0.0, -0.25]),
        matrix(2, 1, b, &[1.0, 2.0]),
        matrix(1, 2, c, &[3.0, 4.0]),
        matrix(1, 1, d, &[5.0]),
    )
}

#[test]
fn discrete_diagnostics_match_blocks_eigenvalues_and_units() {
    let diagonal_model = diagonal(0.5);
    let session = Session::new();
    for (name, count) in [("states", 2), ("inputs", 1), ("outputs", 1)] {
        assert_eq!(
            call(
                &session,
                &format!("num_control_discrete_n_{name}"),
                &diagonal_model
            )
            .unwrap(),
            Value::int(count)
        );
    }
    let (b, c, _) = units();
    for (name, dim, expected) in [
        ("controllability", b, [1.0, 0.5, 2.0, -0.5]),
        ("observability", c, [3.0, 4.0, 1.5, -1.0]),
    ] {
        let result = call(
            &session,
            &format!("num_control_discrete_{name}_matrix"),
            &diagonal_model,
        )
        .unwrap();
        assert_eq!(entries(&result, &[2, 2], dim), expected);
    }
    let Value::Tensor(poles) =
        call(&session, "num_control_discrete_poles", &diagonal_model).unwrap()
    else {
        panic!("expected complex poles")
    };
    assert_eq!(poles.shape, [2]);
    assert_eq!(poles.elem, TypeExpr::Complex);
    let mut real_parts: Vec<f64> = poles
        .data
        .iter()
        .map(|pole| {
            let Value::Complex(pole) = pole else {
                panic!("expected complex pole")
            };
            let (re, im) = pole.to_f64();
            assert_eq!(im, 0.0);
            re
        })
        .collect();
    real_parts.sort_by(f64::total_cmp);
    assert_eq!(real_parts, [-0.25, 0.5]);
    assert_eq!(
        call(&session, "num_control_discrete_is_stable", &diagonal_model).unwrap(),
        Value::Bool(true)
    );
    let (b, c, d) = units();
    let rotation = model(
        matrix(2, 2, DimExpr::Dimless, &[0.0, -0.5, 0.5, 0.0]),
        matrix(2, 1, b, &[1.0, 2.0]),
        matrix(1, 2, c, &[3.0, 4.0]),
        matrix(1, 1, d, &[5.0]),
    );
    let Value::Tensor(poles) = call(&session, "num_control_discrete_poles", &rotation).unwrap()
    else {
        panic!("expected complex poles")
    };
    let mut imaginary: Vec<f64> = poles
        .data
        .iter()
        .map(|pole| {
            let Value::Complex(pole) = pole else {
                panic!("expected complex pole")
            };
            let (re, im) = pole.to_f64();
            assert!(re.abs() < 1e-14);
            im
        })
        .collect();
    imaginary.sort_by(f64::total_cmp);
    assert!((imaginary[0] + 0.5).abs() < 1e-14);
    assert!((imaginary[1] - 0.5).abs() < 1e-14);
    assert_eq!(
        call(&session, "num_control_discrete_is_stable", &diagonal(1.0)).unwrap(),
        Value::Bool(false)
    );
}

#[test]
fn zero_state_model_reports_feedthrough_and_empty_matrices() {
    let (b, c, d) = units();
    let model = model(
        matrix(0, 0, DimExpr::Dimless, &[]),
        matrix(0, 1, b.clone(), &[]),
        matrix(1, 0, c.clone(), &[]),
        matrix(1, 1, d, &[5.0]),
    );
    let session = Session::new();
    assert_eq!(
        call(&session, "num_control_discrete_n_inputs", &model).unwrap(),
        Value::int(1)
    );
    assert_eq!(
        call(&session, "num_control_discrete_is_stable", &model).unwrap(),
        Value::Bool(true)
    );
    for (name, dim) in [("controllability", b), ("observability", c)] {
        let result = call(
            &session,
            &format!("num_control_discrete_{name}_matrix"),
            &model,
        )
        .unwrap();
        assert!(entries(&result, &[0, 0], dim).is_empty());
    }
}

#[test]
fn absent_channels_keep_the_nonempty_state_axis() {
    let (b, c, d) = units();
    let model = model(
        matrix(2, 2, DimExpr::Dimless, &[0.5, 0.0, 0.0, -0.25]),
        matrix(2, 0, b.clone(), &[]),
        matrix(0, 2, c.clone(), &[]),
        matrix(0, 0, d, &[]),
    );
    let session = Session::new();
    for (name, dim, shape) in [("controllability", b, [2, 0]), ("observability", c, [0, 2])] {
        let result = call(
            &session,
            &format!("num_control_discrete_{name}_matrix"),
            &model,
        )
        .unwrap();
        assert!(entries(&result, &shape, dim).is_empty());
    }
}

#[test]
fn checks_unit_shape_budget_and_eigen_failure() {
    let (b, c, d) = units();
    let session = Session::new();
    let bad_dimension = model(
        matrix(1, 1, DimExpr::base("Time").pow(-1), &[0.0]),
        matrix(1, 1, b.clone(), &[1.0]),
        matrix(1, 1, c.clone(), &[1.0]),
        matrix(1, 1, d.clone(), &[1.0]),
    );
    assert!(call(&session, "num_control_discrete_n_states", &bad_dimension).is_err());
    let bad_shape = model(
        matrix(1, 1, DimExpr::Dimless, &[0.0]),
        matrix(2, 1, b.clone(), &[1.0; 2]),
        matrix(1, 1, c.clone(), &[1.0]),
        matrix(1, 1, d.clone(), &[1.0]),
    );
    assert!(
        matches!(call(&session, "num_control_discrete_n_states", &bad_shape),
                     Err(OpError::Domain(message)) if message.contains("need shapes"))
    );
    let oversized = model(
        matrix(450, 450, DimExpr::Dimless, &vec![0.0; 450 * 450]),
        matrix(450, 1, b.clone(), &[1.0; 450]),
        matrix(1, 450, c.clone(), &[1.0; 450]),
        matrix(1, 1, d.clone(), &[1.0]),
    );
    assert!(
        matches!(call(&session, "num_control_discrete_poles", &oversized),
                     Err(OpError::Domain(message)) if message.contains("work limit"))
    );
    assert!(
        matches!(call(&session, "num_control_discrete_controllability_matrix", &oversized),
                     Err(OpError::Domain(message)) if message.contains("value limit"))
    );
    let unrepresentable = model(
        matrix(2, 2, DimExpr::Dimless, &[f64::MAX; 4]),
        matrix(2, 1, b, &[1.0; 2]),
        matrix(1, 2, c, &[1.0; 2]),
        matrix(1, 1, d, &[1.0]),
    );
    assert!(
        matches!(call(&session, "num_control_discrete_is_stable", &unrepresentable),
                     Err(OpError::Domain(message)) if message.contains("finite pole"))
    );
}
