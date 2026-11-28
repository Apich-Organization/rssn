//! Prepared control models are reusable and enforce shape and sample time.
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::TensorValue;
use rssn::kernel::value::Value;

fn float_tensor(
    shape: &[usize],
    data: &[f64],
) -> Value {
    Value::Tensor(TensorValue {
        shape: shape.to_vec(),
        elem: TypeExpr::Float,
        data: data.iter().copied().map(Value::Float).collect(),
    })
}

#[test]
fn discrete_model_coordinates_are_reused_across_queries_and_simulation() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let call = |name, args: &[Value]| session.call_value(name, args, &ctx);
    let mk = |x: f64| float_tensor(&[1, 1], &[x]);
    let prepared = call(
        "num_control_state_space_discrete",
        &[mk(0.5), mk(1.0), mk(1.0), mk(0.0), Value::Float(0.125)],
    )
    .unwrap();
    for name in [
        "num_control_state_space_is_stable",
        "num_control_state_space_is_controllable",
        "num_control_state_space_is_observable",
    ] {
        assert_eq!(
            call(name, std::slice::from_ref(&prepared)).unwrap(),
            Value::Bool(true)
        );
    }
    let get_poles = || {
        call(
            "num_control_state_space_poles",
            std::slice::from_ref(&prepared),
        )
        .unwrap()
    };
    let Value::Tensor(poles) = get_poles() else {
        panic!("expected pole tensor")
    };
    assert_eq!(poles.shape, [1]);
    let Value::Complex(pole) = &poles.data[0] else {
        panic!("expected complex pole")
    };
    assert_eq!(pole.to_f64(), (0.5, 0.0));
    assert_eq!(get_poles(), Value::Tensor(poles));

    let sim = "num_control_state_space_simulate_with_dt";
    let mut args = [
        prepared.clone(),
        float_tensor(&[1], &[1.0]),
        float_tensor(&[3, 1], &[0.0, 0.0, 0.0]),
        Value::Float(0.125),
    ];
    assert_eq!(
        call(sim, &args).unwrap(),
        float_tensor(&[3, 1], &[1.0, 0.5, 0.25])
    );
    args[3] = Value::Float(0.25);
    assert!(call(sim, &args).is_err());
    assert!(
        Session::new()
            .call_value("num_control_state_space_poles", &[prepared], &ctx)
            .is_err()
    );
}

#[test]
fn both_continuous_entry_points_validate_shapes_and_have_identical_poles() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let mk = |x: f64| float_tensor(&[1, 1], &[x]);
    let args = [mk(-2.0), mk(0.0), mk(1.0), mk(0.0)];
    for name in [
        "num_control_state_space_new",
        "num_control_state_space_continuous",
    ] {
        let handle = session.call_value(name, &args, &ctx).unwrap();
        assert_eq!(
            session
                .call_value(
                    "num_control_state_space_is_stable",
                    std::slice::from_ref(&handle),
                    &ctx
                )
                .unwrap(),
            Value::Bool(true)
        );
        assert_eq!(
            session
                .call_value("num_control_state_space_is_controllable", &[handle], &ctx)
                .unwrap(),
            Value::Bool(false)
        );
        let invalid = [
            mk(-2.0),
            float_tensor(&[2, 1], &[1.0, 2.0]),
            mk(1.0),
            mk(0.0),
        ];
        assert!(session.call_value(name, &invalid, &ctx).is_err());
    }
}

#[test]
fn incomplete_spectra_do_not_report_legacy_stability() {
    use rssn::kernels::control::StateSpace;
    let mut model = StateSpace::continuous(
        vec![vec![-1.0]],
        vec![vec![1.0]],
        vec![vec![1.0]],
        vec![vec![0.0]],
    )
    .unwrap();
    model.a[0].clear();
    assert!(model.try_poles().is_err());
    assert!(model.poles().is_empty());
    assert!(!model.is_stable());
}
