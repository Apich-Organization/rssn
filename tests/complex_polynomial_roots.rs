use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::ComplexValue;
use rssn::kernel::value::TensorValue;

fn complex(
    re: f64,
    im: f64,
) -> Value {
    Value::Complex(ComplexValue::Float { re, im })
}

fn field<'a>(
    result: &'a Value,
    name: &str,
) -> &'a Value {
    let Value::Record(record) = result else {
        panic!("expected result record")
    };
    record.get(name).unwrap()
}

fn tensor(value: &Value) -> &[Value] {
    let Value::Tensor(tensor) = value else {
        panic!("expected tensor")
    };
    &tensor.data
}

fn close(
    point: &Value,
    re: f64,
    im: f64,
    tol: f64,
) -> bool {
    let Value::Complex(point) = point else {
        panic!("expected complex root")
    };
    let (x, y) = point.to_f64();
    (x - re).hypot(y - im) < tol
}

fn call(
    session: &Session,
    method: &str,
    coefficients: &[f64],
    seed: Option<Value>,
    steps: u64,
) -> Result<Value, OpError> {
    let mut args = vec![Value::Tensor(TensorValue::of(
        vec![coefficients.len()],
        coefficients.iter().map(|&c| complex(c, 0.0)).collect(),
        TypeExpr::Complex,
    ))];
    args.extend(seed);
    args.extend([Value::int(steps), Value::Float(1e-12)]);
    session.call_value(
        &format!("num_complex_polynomial_{method}"),
        &args,
        &CallCtx::new(),
    )
}

#[test]
fn simultaneous_root_sets_verify_normalization_and_multiplicity() {
    let session = Session::new();
    for method in ["durand_kerner", "aberth_ehrlich"] {
        for scale in [1.0, 1e200, 1e-200] {
            let result = call(&session, method, &[-scale, 0.0, scale, 0.0], None, 96).unwrap();
            assert_eq!(field(&result, "degree"), &Value::int(2));
            let roots = tensor(field(&result, "roots"));
            assert!(roots.iter().any(|root| close(root, -1.0, 0.0, 1e-12)));
            assert!(roots.iter().any(|root| close(root, 1.0, 0.0, 1e-12)));
            assert!(field(&result, "coefficient_error").as_f64().unwrap() < 1e-12);
            assert!(
                tensor(field(&result, "backward_errors"))
                    .iter()
                    .all(|error| error.as_f64().unwrap() < 1e-12)
            );
        }
        let cubic = call(&session, method, &[-1.0, 0.0, 0.0, 1.0], None, 128).unwrap();
        let roots = tensor(field(&cubic, "roots"));
        assert!(roots.iter().any(|root| close(root, 1.0, 0.0, 1e-12)));
        for sign in [-1.0_f64, 1.0] {
            assert!(
                roots
                    .iter()
                    .any(|root| close(root, -0.5, sign * 0.75_f64.sqrt(), 1e-12))
            );
        }
        let repeated = call(&session, method, &[1.0, -2.0, 1.0], None, 256).unwrap();
        assert_eq!(field(&repeated, "degree"), &Value::int(2));
        assert!(
            tensor(field(&repeated, "roots"))
                .iter()
                .all(|root| close(root, 1.0, 0.0, 1e-6))
        );
        let repeated_zero = call(&session, method, &[0.0, 0.0, -1.0, 1.0], None, 1).unwrap();
        let roots = tensor(field(&repeated_zero, "roots"));
        assert_eq!(
            roots
                .iter()
                .filter(|root| close(root, 0.0, 0.0, 1e-12))
                .count(),
            2
        );
        assert!(roots.iter().any(|root| close(root, 1.0, 0.0, 1e-12)));
    }
}

#[test]
fn constants_empty_and_zero_polynomials_are_distinct() {
    let session = Session::new();
    for method in ["durand_kerner", "aberth_ehrlich"] {
        let constant = call(&session, method, &[2.0, 0.0], None, 1).unwrap();
        assert_eq!(field(&constant, "degree"), &Value::int(0));
        assert!(tensor(field(&constant, "roots")).is_empty());
        for zero in [&[] as &[f64], &[0.0, 0.0]] {
            assert!(matches!(
                call(&session, method, zero, None, 1),
                Err(OpError::Domain(_))
            ));
        }
        assert!(matches!(
            call(&session, method, &[-1.0, 0.0, 1.0], None, 0),
            Err(OpError::Domain(_))
        ));
        assert!(matches!(
            call(&session, method, &[-1.0, 0.0, 0.0, 1.0], None, 1),
            Err(OpError::NotReduced(_))
        ));
        assert!(matches!(
            call(&session, method, &[-1.0, 0.0, 1.0], None, 2_000_000),
            Err(OpError::BudgetExceeded)
        ));
    }
}

#[test]
fn laguerre_relative_residual_and_seed_validation() {
    let session = Session::new();
    let tiny = call(
        &session,
        "laguerre",
        &[-1e-200, 1e-200],
        Some(complex(4.0, 0.0)),
        1,
    )
    .unwrap();
    let Value::Variant(root) = tiny else {
        panic!("expected optional root")
    };
    assert_eq!(root.tag.as_str(), "Some");
    assert!(close(field(&root.payload[0], "root"), 1.0, 0.0, 1e-12));
    assert!(field(&root.payload[0], "backward_error").as_f64().unwrap() < 1e-12);
    let positive = call(
        &session,
        "laguerre",
        &[-1.0, 0.0, 1.0],
        Some(complex(2.0, 1.0)),
        32,
    )
    .unwrap();
    assert!(
        matches!(positive, Value::Variant(v) if v.tag=="Some" && close(field(&v.payload[0], "root"), 1.0, 0.0, 1e-10))
    );
    let constant = call(&session, "laguerre", &[2.0], Some(complex(4.0, 0.0)), 1).unwrap();
    assert!(matches!(constant, Value::Variant(v) if v.tag == "None"));
    assert!(matches!(
        call(&session, "laguerre", &[0.0], Some(complex(4.0, 0.0)), 1),
        Err(OpError::Domain(_))
    ));
    assert!(
        call(
            &session,
            "laguerre",
            &[2.0, 1.0],
            Some(complex(f64::INFINITY, 0.0)),
            1
        )
        .is_err()
    );
}
