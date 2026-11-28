//! Prepared SISO transfer functions through the public ABI.
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::DimRuleInfo;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::ComplexValue;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::TensorValue;
use rssn::kernel::value::VariantValue;

fn vector(values: &[f64]) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![values.len()],
        data: values.iter().copied().map(Value::Float).collect(),
        elem: TypeExpr::Float,
    })
}

fn call(
    session: &Session,
    op: &str,
    args: &[Value],
) -> Result<Value, OpError> {
    session.call_value(&format!("num_control_transfer_{op}"), args, &CallCtx::new())
}

fn model(
    session: &Session,
    num: &[f64],
    den: &[f64],
) -> Value {
    call(session, "new", &[vector(num), vector(den)]).unwrap()
}

fn complex(
    re: f64,
    im: f64,
) -> Value {
    Value::Complex(ComplexValue::Float { re, im })
}

fn assert_close(
    got: Value,
    expected: (f64, f64),
) {
    let Value::Complex(got) = got else {
        panic!("complex expected, received {got:?}")
    };
    let (re, im) = got.to_f64();
    assert!((re - expected.0).abs() < 1e-12);
    assert!((im - expected.1).abs() < 1e-12);
}

#[test]
fn prepared_shape_stability_and_roots_are_checked() {
    let session = Session::new();
    let plan = model(&session, &[0.0, 2.0, 4.0], &[0.0, 1.0, 3.0, 2.0]);
    for _ in 0..2 {
        assert_eq!(
            call(&session, "order", std::slice::from_ref(&plan)).unwrap(),
            Value::int(2)
        );
        for predicate in ["is_proper", "is_strictly_proper", "is_stable"] {
            assert_eq!(
                call(&session, predicate, std::slice::from_ref(&plan)).unwrap(),
                Value::Bool(true)
            );
        }
        let zeros = call(&session, "zeros", std::slice::from_ref(&plan)).unwrap();
        let Value::Tensor(zeros) = zeros else {
            panic!("expected zeros")
        };
        assert_eq!(zeros.shape, [1]);
        assert_close(zeros.data[0].clone(), (-2.0, 0.0));
    }
    let failed_poles = model(&session, &[1.0], &[f64::MIN_POSITIVE, f64::MAX]);
    assert!(call(&session, "is_stable", &[failed_poles]).is_err());
    let unstable = model(&session, &[1.0], &[1.0, -1.0]);
    assert_eq!(
        call(&session, "is_stable", &[unstable]).unwrap(),
        Value::Bool(false)
    );
    let zero = model(&session, &[0.0], &[1.0]);
    assert_eq!(
        call(&session, "zeros", &[zero]).unwrap().type_expr(),
        TypeExpr::Tensor {
            rank: rssn::kernel::signature::NatExpr::Lit(1),
            shape: Some(vec![rssn::kernel::signature::NatExpr::Lit(0)]),
            elem: Box::new(TypeExpr::Complex),
        }
    );
}

#[test]
fn connections_match_frequency_algebra_and_reuse_handles() {
    let session = Session::new();
    let model = model(&session, &[1.0], &[1.0, 1.0]);
    let s = complex(0.0, 1.0);
    assert_close(
        call(&session, "evaluate", &[model.clone(), s.clone()]).unwrap(),
        (0.5, -0.5),
    );
    for (name, expected) in [
        ("series", (0.0, -0.5)),
        ("parallel", (1.0, -1.0)),
        ("feedback", (0.6, -0.2)),
    ] {
        let connected = call(&session, name, &[model.clone(), model.clone()]).unwrap();
        assert_close(
            call(&session, "evaluate", &[connected, s.clone()]).unwrap(),
            expected,
        );
    }
    assert!(call(&session, "evaluate", &[model, complex(-1.0, 0.0)]).is_err());
    let large = call(
        &session,
        "new_proper",
        &[vector(&[1e250]), vector(&[1e250])],
    )
    .unwrap();
    assert_close(call(&session, "evaluate", &[large, s]).unwrap(), (1.0, 0.0));
}

#[test]
fn margins_and_step_have_checked_numeric_contracts() {
    let session = Session::new();
    let model = model(&session, &[2.0], &[1.0, 1.0]);
    let margin = call(
        &session,
        "gain_margin_phase_margin",
        std::slice::from_ref(&model),
    )
    .unwrap();
    let Value::Record(margin) = margin else {
        panic!("expected margins")
    };
    for (field, value) in margin.fields {
        let Value::Variant(VariantValue { tag, payload }) = value else {
            panic!("expected option")
        };
        match (field.as_str(), tag.as_str(), payload.as_slice()) {
            | ("gain_db", "None", []) => {},
            | ("phase_degrees", "Some", [Value::Float(value)]) => {
                assert!((*value - 120.0).abs() < 1e-9)
            },
            | _ => panic!("unexpected {field} {tag} {payload:?}"),
        }
    }
    let times: Vec<_> = (0..=100).map(|i| f64::from(i) / 100.0).collect();
    let Value::Tensor(step) =
        call(&session, "step_response", &[model.clone(), vector(&times)]).unwrap()
    else {
        panic!("expected response")
    };
    let [.., Value::Float(last)] = step.data.as_slice() else {
        panic!("expected final sample")
    };
    assert!((*last - 2.0 * (1.0 - (-1.0_f64).exp())).abs() < 1e-8);
    for bad in [[0.0, 0.0], [2.0, 1.0], [-f64::MAX, f64::MAX]] {
        assert!(call(&session, "step_response", &[model.clone(), vector(&bad)]).is_err());
    }
    let marginal = call(&session, "new", &[vector(&[1.0]), vector(&[1.0, 0.0, 1.0])]).unwrap();
    assert!(call(&session, "gain_margin_phase_margin", &[marginal]).is_err());
}

#[test]
fn scope_workload_and_dimension_errors_are_explicit() {
    let session = Session::new();
    assert!(call(&session, "new", &[vector(&[1.0]), vector(&[0.0])]).is_err());
    assert!(call(&session, "new", &[vector(&[1.0; 130]), vector(&[1.0])]).is_err());
    assert!(
        call(
            &session,
            "new_proper",
            &[vector(&[1.0, 0.0]), vector(&[1.0])]
        )
        .is_err()
    );
    let improper = model(&session, &[1.0, 0.0], &[1.0]);
    assert!(call(&session, "step_response", &[improper, vector(&[0.0, 1.0])]).is_err());
    let plan = model(&session, &[1.0], &[1.0, 1.0]);
    let other = Session::new();
    assert!(call(&other, "order", std::slice::from_ref(&plan)).is_err());
    let measured = Value::Dyn(DynValue {
        value: Box::new(complex(1.0, 1.0)),
        dim: DimExpr::base("Time"),
    });
    assert!(call(&session, "evaluate", &[plan.clone(), measured]).is_err());
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
    assert!(call(&session, "order", &[plan]).is_err());
}

#[test]
fn each_callable_has_a_declared_normalized_transfer_signature() {
    let manifest = Session::new().manifest();
    let ops = manifest
        .ops
        .iter()
        .filter(|op| op.name.starts_with("num_control_transfer_"))
        .collect::<Vec<_>>();
    assert_eq!(ops.len(), 13);
    for op in ops {
        assert!(
            op.rust_path
                .as_deref()
                .unwrap()
                .starts_with("crate::kernels::control::transfer::")
        );
        assert_eq!(op.dim_rule, Some(DimRuleInfo::Dimensionless));
        assert_eq!(
            op.signature_source,
            rssn::kernel::manifest::SignatureSource::Declared
        );
        assert_eq!(op.overloads.len(), 1);
    }
}
