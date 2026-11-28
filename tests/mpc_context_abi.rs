//! Bounded MPC preparation, shifted constraints and reusable rollout over the live ABI.
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::DimRuleInfo;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::RecordValue;
use rssn::kernel::value::TensorValue;
use rssn::kernel::value::VariantValue;
use std::sync::Arc;
use std::sync::atomic::AtomicBool;

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
fn matrix(value: f64) -> Value {
    tensor(vec![1, 1], &[value])
}
fn bound(value: Option<f64>) -> Value {
    let (tag, payload) = value.map_or(("None", vec![]), |value| ("Some", vec![vector(&[value])]));
    Value::Variant(VariantValue {
        tag: tag.into(),
        payload,
    })
}
fn options(
    horizon: usize,
    penalty: f64,
    umin: Option<f64>,
    umax: Option<f64>,
    xmax: Option<f64>,
) -> Value {
    Value::Record(RecordValue {
        fields: vec![
            ("horizon".into(), Value::int(horizon)),
            ("control_horizon".into(), Value::int(horizon)),
            ("q".into(), vector(&[1.0])),
            ("r".into(), vector(&[penalty])),
            ("u_min".into(), bound(umin)),
            ("u_max".into(), bound(umax)),
            ("x_min".into(), bound(None)),
            ("x_max".into(), bound(xmax)),
        ],
    })
}
fn call(
    session: &Session,
    name: &str,
    arguments: &[Value],
) -> Result<Value, OpError> {
    session.call_value(
        &format!("num_control_mpc_{name}"),
        arguments,
        &CallCtx::new(),
    )
}
fn prepare(
    session: &Session,
    options: Value,
) -> Result<Value, OpError> {
    call(
        session,
        "prepare",
        &[matrix(1.0), matrix(1.0), Value::Float(0.1), options],
    )
}
fn run(
    session: &Session,
    plan: Value,
    steps: usize,
) -> Result<Value, OpError> {
    call(
        session,
        "closed_loop",
        &[plan, vector(&[0.0]), vector(&[1.0]), Value::int(steps)],
    )
}
fn check_trajectory(
    value: Value,
    samples: &[f64],
) {
    let Value::Tensor(value) = value else {
        panic!("expected MPC trajectory")
    };
    assert_eq!(value.shape, [samples.len(), 1]);
    for (got, expected) in value.data.iter().zip(samples) {
        assert!((got.as_f64().unwrap() - expected).abs() < 1e-10);
    }
}

#[test]
fn prepared_rollout_matches_unconstrained_scalar_solution() {
    let session = Session::new();
    let plan = prepare(&session, options(1, 1.0, None, None, None)).unwrap();
    for _ in 0..2 {
        check_trajectory(
            run(&session, plan.clone(), 3).unwrap(),
            &[0.0, 0.5, 0.75, 0.875],
        );
    }
    check_trajectory(run(&session, plan, 0).unwrap(), &[0.0]);
}

#[test]
fn bounds_shift_after_every_closed_loop_state() {
    let session = Session::new();
    let limited_input = prepare(&session, options(1, 1.0, Some(-0.25), Some(0.25), None)).unwrap();
    check_trajectory(
        run(&session, limited_input, 4).unwrap(),
        &[0.0, 0.25, 0.5, 0.75, 0.875],
    );
    let limited_state = prepare(&session, options(1, 1.0, None, None, Some(0.1))).unwrap();
    check_trajectory(
        run(&session, limited_state, 3).unwrap(),
        &[0.0, 0.1, 0.1, 0.1],
    );
    let infeasible = prepare(&session, options(1, 1.0, Some(0.5), None, Some(0.1))).unwrap();
    assert!(run(&session, infeasible, 1).is_err());
}

#[test]
fn rejects_work_overflow_and_keeps_cancellation_and_release_visible() {
    let session = Session::new();
    assert!(prepare(&session, options(33, 1.0, None, None, None)).is_err());
    assert!(prepare(&session, options(1, 0.0, None, None, None)).is_err());
    assert!(prepare(&session, options(1, 1.0, Some(1.0), Some(-1.0), None)).is_err());
    assert!(
        call(
            &session,
            "prepare",
            &[
                matrix(f64::MAX),
                matrix(1.0),
                Value::Float(1.0),
                options(2, 1.0, None, None, None)
            ]
        )
        .is_err()
    );
    let plan = prepare(&session, options(1, 1.0, None, None, None)).unwrap();
    assert!(run(&session, plan.clone(), 257).is_err());
    let cancelled = rssn::cancellation::install(Arc::new(AtomicBool::new(true)));
    assert!(matches!(
        run(&session, plan.clone(), 1),
        Err(OpError::HostDiagnostic { code: "E0812", .. })
    ));
    drop(cancelled);
    assert_eq!(
        session
            .call_value(
                "num_native_release",
                std::slice::from_ref(&plan),
                &CallCtx::new()
            )
            .unwrap(),
        Value::Bool(true)
    );
    assert!(run(&session, plan, 1).is_err());
}

#[test]
fn trajectory_exports_the_canonical_checked_dimensionless_source() {
    let manifest = Session::new().manifest();
    let op = manifest
        .ops
        .iter()
        .find(|descriptor| descriptor.name == "num_control_mpc_closed_loop")
        .unwrap();
    assert_eq!(
        op.rust_path.as_deref(),
        Some("crate::kernels::control::mpc::mpc_closed_loop")
    );
    assert_eq!(op.dim_rule, Some(DimRuleInfo::Dimensionless));
    assert_eq!(
        op.signature_source,
        rssn::kernel::manifest::SignatureSource::Declared
    );
    assert_eq!(op.overloads.len(), 1);
}
