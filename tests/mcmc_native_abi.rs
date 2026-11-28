//! Checked deterministic streams and measured convergence diagnostics.
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::DimRuleInfo;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::TensorValue;
use rssn::kernels::mcmc::Pcg32;

fn call(
    session: &Session,
    op: &str,
    args: &[Value],
) -> Result<Value, OpError> {
    session.call_value(op, args, &CallCtx::new())
}

fn measured(
    data: &[f64],
    shape: &[usize],
    unit: &str,
) -> Value {
    let dim = DimExpr::base(unit);
    Value::Tensor(TensorValue {
        shape: shape.to_vec(),
        elem: TypeExpr::Float.with_dim(dim.clone()),
        data: data
            .iter()
            .map(|x| {
                Value::Dyn(DynValue {
                    dim: dim.clone(),
                    value: Box::new(Value::Float(*x)),
                })
            })
            .collect(),
    })
}

#[test]
fn exact_pcg_stream_matches_library_and_enforces_session_ownership() {
    let session = Session::new();
    let stream = call(
        &session,
        "num_mcmc_pcg32_new",
        &[Value::int(u64::MAX), Value::int(53)],
    )
    .unwrap();
    let mut expected = Pcg32::new(u64::MAX, 53);
    let args = std::slice::from_ref(&stream);
    for _ in 0..5 {
        assert_eq!(
            call(&session, "num_mcmc_pcg32_next_u32", args).unwrap(),
            Value::int(expected.next_u32())
        );
        assert_eq!(
            call(&session, "num_mcmc_pcg32_next_f64", args).unwrap(),
            Value::Float(expected.next_f64())
        );
        for _ in 0..2 {
            assert_eq!(
                call(&session, "num_mcmc_pcg32_normal", args).unwrap(),
                Value::Float(expected.normal())
            );
        }
    }
    assert!(call(&Session::new(), "num_mcmc_pcg32_normal", args).is_err());
    assert_eq!(
        call(&session, "num_native_release", args).unwrap(),
        Value::Bool(true)
    );
    assert!(call(&session, "num_mcmc_pcg32_normal", args).is_err());
    assert!(
        call(
            &session,
            "num_mcmc_pcg32_new",
            &[Value::int(12), Value::int(u64::MAX)]
        )
        .is_err()
    );
    let seeded = call(&session, "num_mcmc_pcg32_seeded", &[Value::int(8)]).unwrap();
    assert_eq!(
        call(&session, "num_mcmc_pcg32_next_u32", &[seeded]).unwrap(),
        Value::int(Pcg32::seeded(8).next_u32())
    );
}

#[test]
fn diagnostics_preserve_scale_invariance_and_return_dimensionless_values() {
    let session = Session::new();
    let samples = (0..24).map(|i| f64::from(i % 4)).collect::<Vec<_>>();
    let lengths = measured(&samples, &[24], "Length");
    let big = samples
        .iter()
        .map(|value| value * 1e270)
        .collect::<Vec<_>>();
    let times = measured(&big, &[24], "Time");
    for op in ["num_mcmc_effective_sample_size", "num_mcmc_geweke"] {
        let left = call(&session, op, std::slice::from_ref(&lengths))
            .unwrap()
            .as_f64()
            .unwrap();
        let right = call(&session, op, std::slice::from_ref(&times))
            .unwrap()
            .as_f64()
            .unwrap();
        assert!((left - right).abs() < 1e-12, "{op}: {left} vs {right}");
    }
    let Value::Tensor(acf) = call(
        &session,
        "num_mcmc_autocorrelation",
        &[lengths, Value::int(2)],
    )
    .unwrap() else {
        panic!("expected autocorrelation tensor");
    };
    assert_eq!(acf.shape, [3]);
    assert_eq!(acf.elem, TypeExpr::Float);
    assert_eq!(acf.data[0], Value::Float(1.0));
    let paired = measured(&[0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0], &[2, 4], "Length");
    let rhat = call(&session, "num_mcmc_split_rhat", &[paired])
        .unwrap()
        .as_f64()
        .unwrap();
    assert!((rhat - 0.5_f64.sqrt()).abs() < 1e-12);
}

#[test]
fn diagnostics_reject_mixed_units_and_bound_native_work() {
    let session = Session::new();
    let mut samples = measured(&[0.0, 1.0, 2.0, 3.0], &[4], "Length");
    if let Value::Tensor(ref mut tensor) = samples {
        tensor.data[3] = Value::Dyn(DynValue {
            dim: DimExpr::base("Time"),
            value: Box::new(Value::Float(3.0)),
        });
    }
    assert!(call(&session, "num_mcmc_effective_sample_size", &[samples]).is_err());
    let constant = measured(&[1.0; 4], &[4], "Length");
    assert!(
        call(
            &session,
            "num_mcmc_autocorrelation",
            &[constant, Value::int(2)]
        )
        .is_err()
    );
    let samples = measured(&vec![1.0; 4_097], &[4_097], "Length");
    assert!(
        call(
            &session,
            "num_mcmc_effective_sample_size",
            std::slice::from_ref(&samples)
        )
        .is_err()
    );
    assert!(
        call(
            &session,
            "num_mcmc_autocorrelation",
            &[samples, Value::int(4_096)]
        )
        .is_err()
    );
}

#[test]
fn mcmc_manifest_marks_exact_paths_and_dimensionless_diagnostics() {
    let session = Session::new();
    let manifest = session.manifest();
    let entries = manifest
        .ops
        .iter()
        .filter(|entry| entry.name.starts_with("num_mcmc_"))
        .collect::<Vec<_>>();
    assert_eq!(
        entries
            .iter()
            .filter(|entry| entry.rust_path.is_some())
            .count(),
        10
    );
    assert_eq!(entries.len(), 14);
    let entry = entries
        .iter()
        .find(|entry| entry.name == "num_mcmc_split_rhat")
        .unwrap();
    assert_eq!(entry.overloads[0].ret, TypeExpr::Float);
    assert_eq!(entry.dim_rule, Some(DimRuleInfo::Dimensionless));
    assert!(entries.iter().all(|entry| {
        entry.signature_source == rssn::kernel::manifest::SignatureSource::Declared
    }));
}

fn floats(values: &[f64]) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![values.len()],
        elem: TypeExpr::Float,
        data: values.iter().copied().map(Value::Float).collect(),
    })
}

fn vector_type() -> TypeExpr {
    use rssn::kernel::signature::NatExpr;
    TypeExpr::Tensor {
        rank: NatExpr::Lit(1),
        shape: None,
        elem: Box::new(TypeExpr::Float),
    }
}

fn config(draws: usize) -> Value {
    use rssn::kernel::value::RecordValue;
    Value::Record(RecordValue {
        fields: vec![
            ("draws".into(), Value::int(draws)),
            ("burn_in".into(), Value::int(10)),
            ("thin".into(), Value::int(2)),
            ("seed".into(), Value::int(42)),
        ],
    })
}

fn record_field<'a>(
    result: &'a Value,
    field: &str,
) -> &'a Value {
    let Value::Record(record) = result else {
        panic!("expected chain record")
    };
    &record
        .fields
        .iter()
        .find(|(name, _)| name == field)
        .unwrap()
        .1
}

#[test]
fn callback_samplers_replay_and_preserve_actual_draws() {
    use rssn::kernel::plugin::NativeCallback;
    use rssn::kernel::signature::Signature;
    use std::sync::Arc;
    let session = Session::new();
    let mut ctx = CallCtx::new();
    let log = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |args: &[Value]| {
                let Value::Tensor(state) = &args[0] else {
                    panic!("expected a state tensor")
                };
                Ok(Value::Float(
                    state
                        .data
                        .iter()
                        .map(|x| -0.5 * x.as_f64().unwrap().powi(2))
                        .sum(),
                ))
            },
            Signature::new(vec![vector_type()], TypeExpr::Float),
        )))
        .unwrap();
    let gradient = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |args: &[Value]| {
                let Value::Tensor(state) = &args[0] else {
                    panic!("expected a state tensor")
                };
                Ok(floats(
                    &state
                        .data
                        .iter()
                        .map(|x| -x.as_f64().unwrap())
                        .collect::<Vec<_>>(),
                ))
            },
            Signature::new(vec![vector_type()], vector_type()),
        )))
        .unwrap();
    let conditional = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |args: &[Value]| Ok(Value::Float(args[2].as_f64().unwrap())),
            Signature::new(
                vec![TypeExpr::Int, vector_type(), TypeExpr::Float],
                TypeExpr::Float,
            ),
        )))
        .unwrap();
    let invoke = |op, args: &[Value]| session.call_value(op, args, &ctx).unwrap();
    let initial = floats(&[0.0, 0.5]);
    let metropolis = [
        Value::Function(log),
        initial.clone(),
        floats(&[0.7, 0.7]),
        config(64),
    ];
    let mh = invoke("num_mcmc_metropolis_hastings", &metropolis);
    assert_eq!(mh, invoke("num_mcmc_random_walk_metropolis", &metropolis));
    assert_eq!(record_field(&mh, "proposed"), &Value::int(138));
    let hmc = [
        Value::Function(log),
        Value::Function(gradient),
        initial.clone(),
        Value::Float(0.3),
        Value::int(5),
        config(64),
    ];
    let run = invoke("num_mcmc_hmc", &hmc);
    assert_eq!(run, invoke("num_mcmc_hmc", &hmc));
    let Value::Tensor(draws) = record_field(&run, "draws") else {
        panic!("expected draw matrix")
    };
    assert_eq!(draws.shape, [64, 2]);
    let gibbs = invoke(
        "num_mcmc_gibbs",
        &[Value::Function(conditional), initial.clone(), config(64)],
    );
    let logged = invoke(
        "num_mcmc_gibbs_with_log_posterior",
        &[
            Value::Function(conditional),
            Value::Function(log),
            initial,
            config(64),
        ],
    );
    assert_eq!(
        record_field(&gibbs, "draws"),
        record_field(&logged, "draws")
    );
    assert_eq!(record_field(&gibbs, "log_posterior"), &floats(&[]));
    let Value::Tensor(logged) = record_field(&logged, "log_posterior") else {
        panic!("expected log vector")
    };
    assert_eq!(logged.shape, [64]);
}

#[test]
fn callback_failures_are_preserved_and_invalid_work_never_calls_host() {
    use rssn::kernel::plugin::NativeCallback;
    use rssn::kernel::signature::Signature;
    use std::sync::Arc;
    use std::sync::atomic::AtomicUsize;
    use std::sync::atomic::Ordering::Relaxed;
    let session = Session::new();
    let mut ctx = CallCtx::new();
    let counter = Arc::new(AtomicUsize::new(0));
    let callbacks = Arc::clone(&counter);
    let log = ctx
        .add_callback(Arc::new(NativeCallback::new(
            move |_: &[Value]| {
                if callbacks.fetch_add(1, Relaxed) == 0 {
                    Ok(Value::Float(0.0))
                } else {
                    Err(OpError::Host("precise callback failure".into()))
                }
            },
            Signature::new(vec![vector_type()], TypeExpr::Float),
        )))
        .unwrap();
    let run = |draws| {
        session.call_value(
            "num_mcmc_metropolis_hastings",
            &[
                Value::Function(log),
                floats(&[0.0]),
                floats(&[0.5]),
                config(draws),
            ],
            &ctx,
        )
    };
    assert!(matches!(run(20_000), Err(OpError::BudgetExceeded)));
    assert_eq!(counter.load(Relaxed), 0);
    let error = run(64).unwrap_err().to_string();
    assert!(error.contains("precise callback failure"), "{error}");
    assert_eq!(counter.load(Relaxed), 2);
}
