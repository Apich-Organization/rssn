//! Checked PID and ERA trajectory contracts over the live value ABI.
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::DimRuleInfo;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::RecordValue;
use rssn::kernel::value::TensorValue;
use rssn::kernel::value::VariantValue;

fn tensor(
    shape: Vec<usize>,
    values: &[f64],
) -> Value {
    Value::Tensor(TensorValue {
        shape,
        data: values.iter().copied().map(Value::Float).collect(),
        elem: TypeExpr::Float,
    })
}
fn vector(values: &[f64]) -> Value {
    tensor(vec![values.len()], values)
}
fn record(fields: Vec<(&str, Value)>) -> Value {
    Value::Record(RecordValue {
        fields: fields
            .into_iter()
            .map(|(key, value)| (key.into(), value))
            .collect(),
    })
}
fn option(value: Option<Value>) -> Value {
    let (tag, payload) = value.map_or(("None", vec![]), |value| ("Some", vec![value]));
    Value::Variant(VariantValue {
        tag: tag.into(),
        payload,
    })
}
fn gains(
    kp: f64,
    ki: f64,
    kd: f64,
    filter: Option<f64>,
) -> Value {
    record(vec![
        ("kp", Value::Float(kp)),
        ("ki", Value::Float(ki)),
        ("kd", Value::Float(kd)),
        ("filter_time", option(filter.map(Value::Float))),
    ])
}
fn call(
    session: &Session,
    op: &str,
    args: &[Value],
) -> Result<Value, OpError> {
    session.call_value(&format!("num_control_{op}"), args, &CallCtx::new())
}

#[test]
fn positional_pid_filters_derivatives_and_prevents_negative_gain_windup() {
    let session = Session::new();
    let run = |pid, samples, limits| {
        call(
            &session,
            "simulate_pid",
            &[pid, samples, Value::Float(1.0), limits],
        )
    };
    assert_eq!(
        run(
            gains(0.0, 0.0, 1.0, Some(1.0)),
            vector(&[0.0, 1.0, 1.0]),
            option(None)
        )
        .unwrap(),
        vector(&[0.0, 0.5, 0.25]),
    );
    let limits = option(Some(record(vec![
        ("low", Value::Float(-1.0)),
        ("high", Value::Float(1.0)),
    ])));
    assert_eq!(
        run(
            gains(-2.0, -1.0, 0.0, None),
            vector(&[1.0, 1.0, 0.0]),
            limits
        )
        .unwrap(),
        vector(&[-1.0, -1.0, 0.0]),
    );
    let reversed = option(Some(record(vec![
        ("low", Value::Float(1.0)),
        ("high", Value::Float(-1.0)),
    ])));
    assert!(run(gains(1.0, 0.0, 0.0, None), vector(&[1.0]), reversed).is_err());
    assert!(
        run(
            gains(1.0, 0.0, 0.0, Some(0.0)),
            vector(&[1.0]),
            option(None)
        )
        .is_err()
    );
    assert!(
        run(
            gains(f64::MAX, 0.0, 0.0, None),
            vector(&[2.0]),
            option(None)
        )
        .is_err()
    );
}

#[test]
fn pid_plant_response_reuses_a_prepared_transfer_and_checks_time() {
    let session = Session::new();
    let plant = call(
        &session,
        "transfer_new_proper",
        &[vector(&[1.0]), vector(&[1.0, 1.0])],
    )
    .unwrap();
    let pid = gains(1.0, 0.0, 0.0, None);
    let time = [0.0, 0.01, 0.02];
    let response = call(
        &session,
        "pid_response",
        &[pid.clone(), plant.clone(), vector(&time), Value::Float(1.0)],
    )
    .unwrap();
    let Value::Tensor(response) = response else {
        panic!("expected PID samples")
    };
    assert_eq!(response.shape, [3]);
    assert_eq!(response.data[0], Value::Float(0.0));
    let second = response.data[1].as_f64().unwrap();
    assert!((second - (1.0 - (-0.01_f64).exp())).abs() < 1e-10);
    for t in [[0.0, 0.0], [-f64::MAX, f64::MAX]] {
        assert!(
            call(
                &session,
                "pid_response",
                &[pid.clone(), plant.clone(), vector(&t), Value::Float(1.0)]
            )
            .is_err()
        );
    }
    let measured = Value::Dyn(DynValue {
        value: Box::new(Value::Float(1.0)),
        dim: DimExpr::base("Length"),
    });
    assert!(
        call(
            &session,
            "pid_response",
            &[pid, plant, vector(&time), measured]
        )
        .is_err()
    );
}

fn era_model(
    session: &Session,
    scale: f64,
) -> Value {
    let u = [2.0, -1.0, 0.25, 0.0, 0.0, 0.0];
    let mut data = Vec::new();
    let mut state = 0.0;
    for value in u {
        data.extend([scale * state, -2.0 * scale * state]);
        state = 0.5 * state + value;
    }
    call(
        session,
        "era",
        &[tensor(vec![u.len(), 2], &data), vector(&u), Value::int(1)],
    )
    .unwrap()
}

#[test]
fn era_identifies_two_outputs_and_retains_empty_impulse_shape() {
    let session = Session::new();
    for scale in [1.0, 1e-15] {
        let fitted = era_model(&session, scale);
        let matrices = call(&session, "era_matrices", std::slice::from_ref(&fitted)).unwrap();
        let Value::Record(matrices) = matrices else {
            panic!("expected ERA matrices")
        };
        assert_eq!(
            matrices
                .fields
                .iter()
                .find(|(key, _)| key == "order")
                .unwrap()
                .1,
            Value::int(1)
        );
        assert_eq!(
            call(
                &session,
                "impulse_response",
                &[fitted.clone(), Value::int(0)]
            )
            .unwrap(),
            tensor(vec![0, 2], &[])
        );
        let body = call(&session, "impulse_response", &[fitted, Value::int(4)]).unwrap();
        let Value::Tensor(body) = body else {
            panic!("expected impulse samples")
        };
        assert_eq!(body.shape, [4, 2]);
        for (value, expected) in body
            .data
            .iter()
            .zip([0.0, 0.0, 1.0, -2.0, 0.5, -1.0, 0.25, -0.5])
        {
            assert!((value.as_f64().unwrap() / scale - expected).abs() < 1e-12);
        }
    }
}

#[test]
fn era_rejects_invalid_inputs_and_stale_or_foreign_handles() {
    let session = Session::new();
    let valid_y = tensor(vec![3, 1], &[0.0, 1.0, 0.5]);
    assert!(
        call(
            &session,
            "era",
            &[valid_y.clone(), vector(&[0.0, 1.0, 0.0]), Value::int(1)]
        )
        .is_err()
    );
    assert!(
        call(
            &session,
            "era",
            &[valid_y, vector(&[1.0, 0.0, 0.0]), Value::int(usize::MAX)]
        )
        .is_err()
    );
    assert!(
        call(
            &session,
            "era",
            &[tensor(vec![1_000_000, 0], &[]), vector(&[]), Value::int(1)]
        )
        .is_err()
    );
    let fitted = era_model(&session, 1.0);
    assert!(
        call(
            &Session::new(),
            "impulse_response",
            &[fitted.clone(), Value::int(1)]
        )
        .is_err()
    );
    assert!(
        call(
            &session,
            "impulse_response",
            &[fitted.clone(), Value::int(4097)]
        )
        .is_err()
    );
    assert_eq!(
        session
            .call_value(
                "num_native_release",
                std::slice::from_ref(&fitted),
                &CallCtx::new()
            )
            .unwrap(),
        Value::Bool(true)
    );
    assert!(call(&session, "impulse_response", &[fitted, Value::int(1)]).is_err());
}

#[test]
fn four_callable_sources_have_declared_dimensionless_overloads() {
    let manifest = Session::new().manifest();
    for (name, path) in [
        ("num_control_era", "crate::kernels::control::ident::era"),
        (
            "num_control_impulse_response",
            "crate::kernels::control::ident::impulse_response",
        ),
        (
            "num_control_pid_response",
            "crate::kernels::control::pid::pid_response",
        ),
        (
            "num_control_simulate_pid",
            "crate::kernels::control::pid::simulate_pid",
        ),
    ] {
        let op = manifest
            .ops
            .iter()
            .find(|entry| entry.name == name)
            .unwrap();
        assert_eq!(op.rust_path.as_deref(), Some(path));
        assert_eq!(op.dim_rule, Some(DimRuleInfo::Dimensionless));
        assert_eq!(
            op.signature_source,
            rssn::kernel::manifest::SignatureSource::Declared
        );
        assert_eq!(op.overloads.len(), 1);
    }
}
