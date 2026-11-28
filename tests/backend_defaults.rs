//! Default capability selection and portable numerical I/O through the value ABI.

use rssn::api::Session;
use rssn::backend;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::TensorValue;

fn call(
    session: &Session,
    name: &str,
    args: &[Value],
) -> Result<Value, OpError> {
    session.call_value(name, args, &CallCtx::new())
}

fn vector(
    elem: TypeExpr,
    data: Vec<Value>,
) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![data.len()],
        elem,
        data,
    })
}

fn string(value: &str) -> Value {
    Value::String(value.into())
}

#[test]
fn named_backends_execute_measured_expressions_and_auto_restores_scope() {
    let session = Session::new();
    let Value::Tensor(available) = call(&session, "num_backend_available", &[]).unwrap() else {
        panic!("backend selectors are reported as a string vector");
    };
    let auto = if available.data.contains(&string("tiered")) {
        "tiered"
    } else {
        "interpreter"
    };
    let current = || call(&session, "num_backend_current_name", &[]).unwrap();
    let select = |name| call(&session, "num_backend_select_default", &[string(name)]);
    assert_eq!(current(), string(auto));
    assert_eq!(select("cpu-parallel").unwrap(), string("cpu-parallel"));
    assert_eq!(current(), string("cpu-parallel"));
    assert!(matches!(
        select("gpu:missing:99"),
        Err(OpError::Unsupported { .. })
    ));
    assert_eq!(current(), string("cpu-parallel"));

    let scope_result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        backend::with_backend(backend::named("interpreter").unwrap(), || {
            assert_eq!(current(), string("interpreter"));
            panic!("scope test");
        });
    }));
    assert!(scope_result.is_err());
    assert_eq!(current(), string("cpu-parallel"));
    assert_eq!(select("auto").unwrap(), string(auto));
    assert_eq!(current(), string(auto));

    let reference = call(&session, "num_backend_reference", &[Value::Float(1.0)]).unwrap();
    let input = call(&session, "num_backend_value", &[Value::Float(2.5)]).unwrap();
    let operands = vector(
        TypeExpr::Opaque("rssn/native/backend-scalar/v1".into()),
        vec![input],
    );
    let references = vector(
        TypeExpr::Opaque("rssn/native/backend-scalar/v1".into()),
        vec![reference.clone()],
    );
    for selector in available.data {
        let Value::String(ref name) = selector else {
            panic!("invalid backend selector cell");
        };
        let named = call(
            &session,
            "num_backend_named",
            std::slice::from_ref(&selector),
        )
        .unwrap();
        assert_eq!(
            call(&session, "num_backend_name", std::slice::from_ref(&named)).unwrap(),
            if name == "auto" {
                string(auto)
            } else {
                selector.clone()
            }
        );
        let program = call(
            &session,
            "num_backend_compile",
            &[
                named.clone(),
                string("x*x+3"),
                vector(TypeExpr::String, vec![string("x")]),
                references.clone(),
                reference.clone(),
                Value::int(1_000),
            ],
        )
        .unwrap_or_else(|error| panic!("{name} compile: {error}"));
        assert_eq!(
            call(
                &session,
                "num_backend_call",
                &[program.clone(), operands.clone()]
            )
            .unwrap()
            .as_f64(),
            Some(9.25),
            "{name} arithmetic"
        );
        assert_eq!(
            call(&session, "num_backend_release", &[program]).unwrap(),
            Value::Bool(true)
        );
        assert_eq!(
            call(&session, "num_backend_release", &[named]).unwrap(),
            Value::Bool(true)
        );
    }
}

#[test]
fn portable_npy_reports_feature_availability_and_retains_ieee_values() {
    let session = Session::new();
    let encode = |shape: Vec<Value>, data: Vec<Value>| {
        call(
            &session,
            "num_npy_encode",
            &[vector(TypeExpr::Int, shape), vector(TypeExpr::Float, data)],
        )
    };
    let scalar = encode(vec![], vec![Value::Float(-0.0)]);
    if !cfg!(feature = "npy") {
        assert!(matches!(scalar, Err(OpError::Unsupported { .. })));
        return;
    }
    let bits = [
        0x8000_0000_0000_0000,
        0x7ff8_0000_0000_0001,
        0x7ff0_0000_0000_0000,
    ];
    for (shape, cases) in [
        (vec![], &bits[..1]),
        (vec![1, 3], bits.as_slice()),
        (vec![2, 0, 3], &bits[..0]),
    ] {
        let extents = shape.into_iter().map(Value::int).collect::<Vec<_>>();
        let encoded = encode(
            extents.clone(),
            cases
                .iter()
                .map(|b| Value::Float(f64::from_bits(*b)))
                .collect(),
        )
        .unwrap();
        let Value::String(ref hex) = encoded else {
            panic!("expected NPY hex string");
        };
        assert!(hex.starts_with("934e554d50590100"));
        let result = call(&session, "num_npy_decode", &[encoded]).unwrap();
        let Value::Record(decoded) = result else {
            panic!("expected NPY decoded record");
        };
        assert_eq!(decoded.get("shape"), Some(&vector(TypeExpr::Int, extents)));
        let Some(Value::Tensor(data)) = decoded.get("data") else {
            panic!("expected NPY data vector");
        };
        let actual = data
            .data
            .iter()
            .map(|v| v.as_f64().unwrap().to_bits())
            .collect::<Vec<_>>();
        assert_eq!(actual, cases);
    }
    assert!(encode(vec![Value::int(2)], vec![Value::Float(1.0)]).is_err());
    assert!(call(&session, "num_npy_decode", &[string("xx")]).is_err());
}
