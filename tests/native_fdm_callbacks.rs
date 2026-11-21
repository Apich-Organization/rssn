//! The generic native adaptive PDE controller uses the selected compiled step and event.
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::RecordValue;
use rssn::kernel::value::TensorValue;
use rssn::kernel::value::VariantValue;
fn call(
    s: &Session,
    n: &str,
    a: &[Value],
) -> Value {
    s.call_value(n, a, &CallCtx::new()).unwrap()
}
fn vector(
    v: Vec<Value>,
    elem: TypeExpr,
) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![v.len()],
        data: v,
        elem,
    })
}
fn floats(v: &[f64]) -> Value {
    vector(
        v.iter().copied().map(Value::Float).collect(),
        TypeExpr::Float,
    )
}
fn measured(
    x: f64,
    d: &str,
) -> Value {
    Value::Dyn(DynValue {
        value: Box::new(Value::Float(x)),
        dim: DimExpr::Base(d.into()),
    })
}
fn field<'a>(
    v: &'a Value,
    n: &str,
) -> &'a Value {
    let Value::Record(r) = v else {
        panic!("record")
    };
    &r.fields.iter().find(|(k, _)| k == n).unwrap().1
}
fn compile(
    s: &Session,
    expr: &str,
    step: bool,
) -> Value {
    let backend = call(s, "num_backend_interpreter", &[]);
    let mut refs = vec![measured(1.0, "Time")];
    if step {
        refs.push(measured(1.0, "Time"));
    }
    refs.push(measured(3.0, "Length"));
    let refs = refs
        .into_iter()
        .map(|v| call(s, "num_backend_reference", &[v]))
        .collect();
    let output = call(
        s,
        "num_backend_reference",
        &[if step {
            measured(3.0, "Length")
        } else {
            Value::Float(1.0)
        }],
    );
    let names = if step {
        vec!["t", "h", "u"]
    } else {
        vec!["t", "u"]
    };
    call(
        s,
        "num_backend_compile",
        &[
            backend,
            Value::String(expr.into()),
            vector(
                names.into_iter().map(|n| Value::String(n.into())).collect(),
                TypeExpr::String,
            ),
            vector(
                refs,
                TypeExpr::Opaque("rssn/native/backend-scalar/v1".into()),
            ),
            output,
            Value::int(100000),
        ],
    )
}
#[test]
fn generic_native_step_doubling_rejects_euler_and_localizes_coherent_event() {
    let s = Session::new();
    let frame = call(
        &s,
        "num_grid_frame_new",
        &[
            measured(2.0, "Length"),
            measured(2.0, "Time"),
            measured(3.0, "Length"),
        ],
    );
    let step = compile(&s, "u*(1-h)", true);
    let event = compile(&s, "u-0.5", false);
    let mut fields = vec![
        ("t0", 0.0),
        ("t1", 0.75),
        ("rtol", 1e-5),
        ("atol", 1e-9),
        ("h0", 0.5),
        ("h_min", 1e-8),
        ("h_max", 0.5),
    ]
    .into_iter()
    .map(|(n, v)| (n.into(), Value::Float(v)))
    .collect::<Vec<_>>();
    fields.extend(
        [
            ("max_steps", 1024),
            ("max_step_calls", 4096),
            ("max_rejections", 128),
        ]
        .into_iter()
        .map(|(n, v)| (n.into(), Value::int(v))),
    );
    let none = Value::Variant(VariantValue {
        tag: "None".into(),
        payload: vec![],
    });
    let out = call(
        &s,
        "num_grid_fdm_adaptive_evolve",
        &[
            frame,
            floats(&[1.0]),
            Value::int(1),
            floats(&[1e-9]),
            Value::Record(RecordValue { fields }),
            vector(
                vec![step],
                TypeExpr::Opaque("rssn/native/compiled/v1".into()),
            ),
            none,
            event,
            Value::int(-1),
            Value::int(10000),
        ],
    );
    assert!(field(&out, "rejected").as_f64().unwrap() > 0.0);
    let Value::Variant(event) = field(&out, "event_time") else {
        panic!("event")
    };
    assert_eq!(event.tag, "Some");
    let Value::Dyn(time) = &event.payload[0] else {
        panic!("measured event time")
    };
    let time = time.value.as_f64().unwrap();
    assert!((time - std::f64::consts::LN_2).abs() < 0.003);
    let Value::Tensor(times) = field(&out, "times") else {
        panic!("times")
    };
    assert!((times.data.last().unwrap().as_f64().unwrap() - time).abs() < 1e-12);
}
