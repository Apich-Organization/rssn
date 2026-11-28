//! Declared quadratic dynamics checks via the normal numerical registry.

use rssn::api::Session;
use rssn::kernel::Value;
use rssn::kernel::call::CallCtx;
use rssn::kernel::call::OpError;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::ComplexValue;
use rssn::kernel::value::TensorValue;

fn complex(
    re: f64,
    im: f64,
) -> Value {
    Value::Complex(ComplexValue::Float { re, im })
}
fn real(v: &Value) -> f64 {
    v.as_f64().unwrap()
}
fn parts(v: &Value) -> (f64, f64) {
    let Value::Complex(z) = v else {
        panic!("expected complex value")
    };
    z.to_f64()
}
fn field<'a>(
    v: &'a Value,
    name: &str,
) -> &'a Value {
    let Value::Record(r) = v else {
        panic!("expected result record")
    };
    r.get(name).unwrap()
}
fn array(v: &Value) -> &[Value] {
    let Value::Tensor(t) = v else {
        panic!("expected result tensor")
    };
    &t.data
}
fn call(
    s: &Session,
    name: &str,
    args: &[Value],
) -> Value {
    s.call_value(
        &format!("num_complex_dynamics_{name}"),
        args,
        &CallCtx::new(),
    )
    .unwrap()
}

#[test]
fn known_bounded_components_and_exterior_estimates() {
    let s = Session::new();
    for c in [complex(0.0, 0.0), complex(-1.0, 0.0)] {
        assert_eq!(
            call(&s, "in_main_components", std::slice::from_ref(&c)),
            Value::Bool(true)
        );
        let r = call(
            &s,
            "mandelbrot_escape",
            &[c.clone(), Value::int(32), Value::Float(2.0)],
        );
        assert_eq!(field(&r, "escaped"), &Value::Bool(false));
        assert_eq!(field(&r, "steps"), &Value::int(32));
        assert!(matches!(field(&r,"potential"),Value::Variant(v) if v.tag=="None"));
        let r = call(&s, "mandelbrot_distance", &[c, Value::int(32)]);
        assert!(matches!(r,Value::Variant(v) if v.tag=="None"));
    }
    for c in [complex(1.0, 1.0), complex(0.0, 1e200)] {
        assert_eq!(call(&s, "in_main_components", &[c]), Value::Bool(false));
    }
    let r = call(
        &s,
        "mandelbrot_escape",
        &[complex(1e200, 0.0), Value::int(4), Value::Float(2.0)],
    );
    assert_eq!(field(&r, "escaped"), &Value::Bool(true));
    assert_eq!(field(&r, "steps"), &Value::int(1));
    let Value::Variant(p) = field(&r, "potential") else {
        panic!("expected optional potential")
    };
    assert_eq!(p.tag, "Some");
    assert!((real(&p.payload[0]) - 1e200_f64.ln()).abs() < 1e-12);
    let r = call(
        &s,
        "julia_escape",
        &[
            complex(3.0, 0.0),
            complex(0.0, 0.0),
            Value::int(4),
            Value::Float(2.0),
        ],
    );
    assert_eq!(field(&r, "escaped"), &Value::Bool(true));
    assert_eq!(field(&r, "steps"), &Value::int(0));
    let d = call(
        &s,
        "mandelbrot_distance",
        &[complex(2.0, 0.0), Value::int(64)],
    );
    assert!(matches!(d,Value::Variant(v) if v.tag=="Some" && real(&v.payload[0])>0.0));
}

#[test]
fn known_orbits_primitive_cycles_and_four_multiplier_classes() {
    let s = Session::new();
    let r = call(
        &s,
        "quadratic_orbit",
        &[complex(0.0, 0.0), complex(-1.0, 0.0), Value::int(4)],
    );
    assert_eq!(
        array(&r).iter().map(parts).collect::<Vec<_>>(),
        [(0.0, 0.0), (-1.0, 0.0), (0.0, 0.0), (-1.0, 0.0), (0.0, 0.0)]
    );
    for (z, c, p, class, m) in [
        (0.0, -1.0, 2, "superattracting", 0.0),
        (0.25, 0.1875, 1, "attracting", 0.5),
        (0.5, 0.25, 1, "neutral", 1.0),
        (1.0, 0.0, 1, "repelling", 2.0),
    ] {
        let r = call(
            &s,
            "quadratic_multiplier",
            &[
                complex(z, 0.0),
                complex(c, 0.0),
                Value::int(p),
                Value::Float(1e-12),
            ],
        );
        assert_eq!(field(&r, "kind"), &Value::String(class.into()));
        assert_eq!(parts(field(&r, "multiplier")), (m, 0.0));
        assert_eq!(real(field(&r, "residual")), 0.0);
    }
    let seeds = Value::Tensor(TensorValue {
        shape: vec![4],
        elem: TypeExpr::Complex,
        data: vec![
            complex(0.0, 0.0),
            complex(1.0, 0.0),
            complex(-0.5, 0.866),
            complex(-0.5, -0.866),
        ],
    });
    let primitive = call(
        &s,
        "quadratic_periodic_points",
        &[
            complex(0.0, 0.0),
            Value::int(2),
            seeds,
            Value::int(20),
            Value::Float(1e-11),
        ],
    );
    let coords: Vec<_> = array(field(&primitive, "points"))
        .iter()
        .map(parts)
        .collect();
    assert_eq!(coords.len(), 2);
    assert!(coords.iter().all(|&(x, y)| {
        (x + 0.5).abs() < 1e-10 && (y.abs() - 3.0_f64.sqrt() * 0.5).abs() < 1e-10
    }));
    assert!(
        array(field(&primitive, "residuals"))
            .iter()
            .all(|v| real(v) < 1e-11)
    );
    let n = call(
        &s,
        "quadratic_lyapunov",
        &[
            complex(1.0, 0.0),
            complex(0.0, 0.0),
            Value::int(10),
            Value::int(2),
        ],
    );
    assert!((real(&n) - 2.0_f64.ln()).abs() < 1e-13);
}

#[test]
fn monomial_coordinate_ray_and_invalid_input() {
    let s = Session::new();
    for z in [complex(3.0, 2.0), complex(1e170, 1e170)] {
        let r = call(
            &s,
            "bottcher",
            &[z.clone(), complex(0.0, 0.0), Value::int(16)],
        );
        assert_eq!(parts(&r), parts(&z));
    }
    let r = call(
        &s,
        "external_ray",
        &[
            complex(0.0, 0.0),
            Value::Float(0.25),
            Value::Float(2.0_f64.ln()),
            Value::int(8),
        ],
    );
    let (x, y) = parts(field(&r, "point"));
    assert!(x.abs() < 1e-10 && (y - 2.0).abs() < 1e-10);
    assert!(real(field(&r, "residual")) < 1e-10);
    let ctx = CallCtx::new();
    assert!(matches!(
        s.call_value(
            "num_complex_dynamics_quadratic_orbit",
            &[complex(0.0, 0.0), complex(0.0, 0.0), Value::int(16_385)],
            &ctx
        ),
        Err(OpError::BudgetExceeded)
    ));
    assert!(
        s.call_value(
            "num_complex_dynamics_quadratic_multiplier",
            &[
                complex(0.3, 0.0),
                complex(0.0, 0.0),
                Value::int(1),
                Value::Float(1e-12)
            ],
            &ctx
        )
        .is_err()
    );
    assert!(
        s.call_value(
            "num_complex_dynamics_external_ray",
            &[
                complex(0.0, 0.0),
                Value::Float(0.25),
                Value::Float(400.0),
                Value::int(8)
            ],
            &ctx
        )
        .is_err()
    );
}
