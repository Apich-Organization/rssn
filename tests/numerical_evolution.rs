//! Checked host callbacks expose scalar roots and dimensioned vector IVPs.
use rssn::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::NativeCallback;
use rssn::kernel::OpError;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::NatExpr;
use rssn::kernel::signature::Signature;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::TensorValue;
use rssn::kernel::value::Value;
use std::sync::Arc;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering;

fn number(value: &Value) -> f64 {
    match value {
        | Value::Dyn(value) => number(&value.value),
        | Value::Quantity(value) => number(&value.value),
        | _ => value.as_f64().unwrap(),
    }
}
fn field<'a>(
    value: &'a Value,
    name: &str,
) -> &'a Value {
    let Value::Record(record) = value else {
        panic!("expected record")
    };
    &record.fields.iter().find(|(key, _)| key == name).unwrap().1
}
fn dimensional(
    value: f64,
    dim: DimExpr,
) -> Value {
    Value::Dyn(DynValue {
        value: Box::new(Value::Float(value)),
        dim,
    })
}
fn scalar(dim: DimExpr) -> TypeExpr {
    TypeExpr::Dim(Box::new(TypeExpr::Float), dim)
}
fn vector_type(dim: DimExpr) -> TypeExpr {
    TypeExpr::Tensor {
        rank: NatExpr::Lit(1),
        shape: None,
        elem: Box::new(scalar(dim)),
    }
}
fn vector(
    values: &[f64],
    dim: DimExpr,
) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![values.len()],
        data: values.iter().copied().map(Value::Float).collect(),
        elem: scalar(dim),
    })
}
fn time() -> DimExpr {
    DimExpr::Base("Time".into())
}
fn length() -> DimExpr {
    DimExpr::Base("Length".into())
}
fn rate() -> DimExpr {
    DimExpr::Mul(vec![length(), DimExpr::Inverse(Box::new(time()))])
}
fn methods() -> [&'static str; 4] {
    [
        "num_ode_dopri5",
        "num_ode_rosenbrock23",
        "num_ode_radau5",
        "num_ode_bdf",
    ]
}
fn ode_arguments() -> Vec<Value> {
    vec![
        dimensional(0.0, time()),
        dimensional(0.5, time()),
        vector(&[1.0, 2.0], length()),
        dimensional(1e-10, length()),
        Value::Float(1e-8),
        dimensional(0.1, time()),
        Value::int(10000),
    ]
}
fn ode_signature() -> Signature {
    Signature::new(
        vec![scalar(time()), vector_type(length())],
        vector_type(rate()),
    )
}

#[test]
fn every_vector_ivp_exposes_units_shapes_counters_and_pinned_dispatch() {
    let session = Session::new();
    for name in methods() {
        let calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&calls);
        let mut ctx = CallCtx::new();
        let handle = ctx
            .add_callback(Arc::new(NativeCallback::new(
                move |arguments: &[Value]| {
                    observed.fetch_add(1, Ordering::Relaxed);
                    assert!((0.0..=0.5).contains(&number(&arguments[0])));
                    let Value::Tensor(state) = &arguments[1] else {
                        panic!("state tensor")
                    };
                    assert_eq!(state.shape, vec![2]);
                    Ok(vector(&[3.0, 4.0], rate()))
                },
                ode_signature(),
            )))
            .unwrap();
        let args = ode_arguments();
        let mut combined = args.clone();
        combined.push(Value::Function(handle));
        let key = session
            .resolve_overload_with_context(name, &combined, &ctx)
            .unwrap();
        for pinned in [false, true] {
            calls.store(0, Ordering::Relaxed);
            let value = if pinned {
                session.call_key(&key, &combined, &ctx)
            } else {
                session.call_ho(name, &args, &[handle], &ctx)
            }
            .unwrap_or_else(|error| panic!("{name}: {error:?}"));
            let Value::Tensor(times) = field(&value, "times") else {
                panic!("times")
            };
            let Value::Tensor(states) = field(&value, "states") else {
                panic!("states")
            };
            assert_eq!(times.shape, vec![times.data.len()]);
            assert_eq!(number(times.data.last().unwrap()), 0.5);
            assert_eq!(states.shape, vec![times.data.len(), 2]);
            assert!(
                (number(&states.data[states.data.len() - 2]) - 2.5).abs() < 1e-8,
                "{name}"
            );
            assert!(
                (number(states.data.last().unwrap()) - 4.0).abs() < 1e-8,
                "{name}"
            );
            assert_eq!(times.elem, scalar(time()));
            assert_eq!(states.elem, scalar(length()));
            assert_eq!(
                number(field(&value, "evaluations")) as usize,
                calls.load(Ordering::Relaxed),
                "{name}"
            );
        }
    }
}

#[test]
fn fixed_width_callbacks_instantiate_the_operator_state_extent() {
    let session = Session::new();
    for name in methods() {
        let mut signature = ode_signature();
        signature.params[1] = scalar(length()).vector(2u64);
        signature.ret = scalar(rate()).vector(2u64);
        let mut ctx = CallCtx::new();
        let callback = ctx
            .add_callback(Arc::new(NativeCallback::new(
                |_: &[Value]| Ok(vector(&[0.0, 0.0], rate())),
                signature,
            )))
            .unwrap();
        let result = session
            .call_ho(name, &ode_arguments(), &[callback], &ctx)
            .unwrap_or_else(|error| panic!("{name}: {error:?}"));
        let Value::Tensor(states) = field(&result, "states") else {
            panic!("state tensor")
        };
        assert_eq!(states.shape[1], 2);
    }
}

#[test]
fn vector_derivative_dimensions_and_fixed_extents_are_checked_before_host_calls() {
    let session = Session::new();
    for name in methods() {
        for wrong_extent in [false, true] {
            let calls = Arc::new(AtomicUsize::new(0));
            let observed = Arc::clone(&calls);
            let mut signature = ode_signature();
            if wrong_extent {
                signature.ret = scalar(rate()).vector(3u64);
            } else {
                signature.ret = vector_type(length());
            }
            let mut ctx = CallCtx::new();
            let callback = ctx
                .add_callback(Arc::new(NativeCallback::new(
                    move |_: &[Value]| {
                        observed.fetch_add(1, Ordering::Relaxed);
                        Ok(vector(&[0.0, 0.0, 0.0], rate()))
                    },
                    signature,
                )))
                .unwrap();
            assert!(
                session
                    .call_ho(name, &ode_arguments(), &[callback], &ctx)
                    .is_err(),
                "{name}"
            );
            assert_eq!(calls.load(Ordering::Relaxed), 0, "{name}");
        }
    }
}

#[test]
fn malformed_vector_results_are_errors_and_original_host_failures_win() {
    let session = Session::new();
    for name in methods() {
        for malformed in [false, true] {
            let failures = Arc::new(AtomicUsize::new(0));
            let observed = Arc::clone(&failures);
            let mut ctx = CallCtx::new();
            let handle = ctx
                .add_callback(Arc::new(NativeCallback::new(
                    move |arguments: &[Value]| {
                        if number(&arguments[0]) > 0.1 {
                            observed.fetch_add(1, Ordering::Relaxed);
                            if malformed {
                                Ok(vector(&[1.0, 2.0, 3.0], rate()))
                            } else {
                                Err(OpError::HostDiagnostic {
                                    code: "E0802",
                                    message: "RHS fixture".into(),
                                })
                            }
                        } else {
                            Ok(vector(&[0.0, 0.0], rate()))
                        }
                    },
                    ode_signature(),
                )))
                .unwrap();
            let error = session
                .call_ho(name, &ode_arguments(), &[handle], &ctx)
                .unwrap_err();
            assert_eq!(
                failures.load(Ordering::Relaxed),
                1,
                "{name}: repeated a failed host callback"
            );
            if !malformed {
                assert!(
                    matches!(error, OpError::HostDiagnostic { code: "E0802", .. }),
                    "{name}: {error:?}"
                );
            }
        }
    }
}

#[test]
fn invalid_vector_problem_options_do_not_call_the_host() {
    let session = Session::new();
    for name in methods() {
        let calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&calls);
        let mut ctx = CallCtx::new();
        let callback = ctx
            .add_callback(Arc::new(NativeCallback::new(
                move |_: &[Value]| {
                    observed.fetch_add(1, Ordering::Relaxed);
                    Ok(vector(&[0.0, 0.0], rate()))
                },
                ode_signature(),
            )))
            .unwrap();
        for index in [3usize, 4, 5, 6] {
            let mut args = ode_arguments();
            args[index] = match index {
                | 3 => dimensional(-1.0, length()),
                | 4 => Value::Float(-1.0),
                | 5 => dimensional(0.0, time()),
                | _ => Value::int(0),
            };
            assert!(
                session.call_ho(name, &args, &[callback], &ctx).is_err(),
                "{name}"
            );
        }
        assert_eq!(calls.load(Ordering::Relaxed), 0, "{name}");
    }
}

#[test]
fn brent_root_has_independent_residual_units_and_no_extra_callback_evaluation() {
    let session = Session::new();
    let mut ctx = CallCtx::new();
    let calls = Arc::new(AtomicUsize::new(0));
    let observed = Arc::clone(&calls);
    let callback = ctx
        .add_callback(Arc::new(NativeCallback::new(
            move |arguments: &[Value]| {
                observed.fetch_add(1, Ordering::Relaxed);
                Ok(dimensional(number(&arguments[0]) - 2.0, time()))
            },
            Signature::new(vec![scalar(length())], scalar(time())),
        )))
        .unwrap();
    let arguments = vec![
        dimensional(0.0, length()),
        dimensional(4.0, length()),
        dimensional(1e-12, length()),
        Value::int(100),
    ];
    let mut all = arguments.clone();
    all.push(Value::Function(callback));
    let key = session
        .resolve_overload_with_context("num_brent", &all, &ctx)
        .unwrap();
    for pinned in [false, true] {
        calls.store(0, Ordering::Relaxed);
        let result = if pinned {
            session.call_key(&key, &all, &ctx)
        } else {
            session.call_ho("num_brent", &arguments, &[callback], &ctx)
        }
        .unwrap();
        assert!((number(field(&result, "root")) - 2.0).abs() < 1e-12);
        assert_eq!(
            ctx.value_type(field(&result, "root")).unwrap(),
            scalar(length())
        );
        assert_eq!(
            number(field(&result, "evaluations")) as usize,
            calls.load(Ordering::Relaxed)
        );
    }
}

#[test]
fn root_budget_bracket_and_first_error_contracts_are_explicit() {
    let session = Session::new();
    let mut ctx = CallCtx::new();
    let calls = Arc::new(AtomicUsize::new(0));
    let observed = Arc::clone(&calls);
    let callback = ctx
        .add_callback(Arc::new(NativeCallback::new(
            move |_: &[Value]| {
                observed.fetch_add(1, Ordering::Relaxed);
                Err(OpError::HostDiagnostic {
                    code: "E0802",
                    message: "root fixture".into(),
                })
            },
            Signature::new(vec![TypeExpr::Float], TypeExpr::Float),
        )))
        .unwrap();
    let mut args = vec![
        Value::Float(0.0),
        Value::Float(1.0),
        Value::Float(1e-8),
        Value::int(0),
    ];
    assert!(
        session
            .call_ho("num_brent", &args, &[callback], &ctx)
            .is_err()
    );
    assert_eq!(calls.load(Ordering::Relaxed), 0);
    args[3] = Value::int(10);
    assert!(matches!(
        session.call_ho("num_brent", &args, &[callback], &ctx),
        Err(OpError::HostDiagnostic { code: "E0802", .. })
    ));
    assert_eq!(calls.load(Ordering::Relaxed), 1);
    let constant = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |_: &[Value]| Ok(Value::Float(1.0)),
            Signature::new(vec![TypeExpr::Float], TypeExpr::Float),
        )))
        .unwrap();
    assert!(matches!(
        session.call_ho("num_brent", &args, &[constant], &ctx),
        Err(OpError::Domain(_))
    ));
}
