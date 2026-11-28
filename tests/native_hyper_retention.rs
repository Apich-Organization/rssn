//! Native hyperelastic compiled program ownership and measured observer contracts.
use rssn::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::NativeCallback;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::overload::matches_type;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::NatExpr;
use rssn::kernel::signature::Signature;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::TensorValue;
use std::sync::Arc;
use std::sync::Mutex;

const COMPILED: &str = "rssn/native/compiled/v1";
const SCALAR: &str = "rssn/native/backend-scalar/v1";
const CONTACT: &str = "rssn/native/penalty-contact/v1";

fn call(
    session: &Session,
    name: &str,
    args: &[Value],
) -> Value {
    session
        .call_value(name, args, &CallCtx::new())
        .unwrap_or_else(|error| panic!("{name}: {error}"))
}

fn vector(
    data: Vec<Value>,
    elem: TypeExpr,
) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![data.len()],
        elem,
        data,
    })
}

fn floats(xs: &[f64]) -> Value {
    vector(
        xs.iter().copied().map(Value::Float).collect(),
        TypeExpr::Float,
    )
}

fn opaque(
    data: Vec<Value>,
    kind: &str,
) -> Value {
    vector(data, TypeExpr::Opaque(kind.into()))
}

fn measured(
    value: f64,
    dim: DimExpr,
) -> Value {
    Value::Dyn(DynValue {
        dim,
        value: Box::new(Value::Float(value)),
    })
}

fn length() -> DimExpr {
    DimExpr::base("Length")
}
fn time() -> DimExpr {
    DimExpr::base("Time")
}

fn refs(session: &Session) -> Value {
    call(
        session,
        "num_grid_mechanical_references",
        &[
            measured(2., length()),
            measured(4., time()),
            measured(3., DimExpr::base("Mass")),
            measured(1., DimExpr::base("Temperature")),
        ],
    )
}

fn program(
    session: &Session,
    expression: &str,
    value: f64,
    dimension: DimExpr,
) -> Value {
    let coordinate = call(session, "num_backend_reference", &[measured(1., length())]);
    let output = call(
        session,
        "num_backend_reference",
        &[measured(value, dimension)],
    );
    call(
        session,
        "num_backend_compile_current",
        &[
            Value::String(expression.into()),
            vector(
                ["x", "y", "z"].map(|s| Value::String(s.into())).to_vec(),
                TypeExpr::String,
            ),
            opaque(vec![coordinate; 3], SCALAR),
            output,
            Value::int(10_000),
        ],
    )
}

fn field<'a>(
    record: &'a Value,
    key: &str,
) -> &'a Value {
    let Value::Record(record) = record else {
        panic!("expected record: {record:?}")
    };
    record.get(key).unwrap()
}

fn number(value: &Value) -> f64 {
    match value {
        | Value::Dyn(d) => number(&d.value),
        | Value::Quantity(q) => number(&q.value),
        | _ => value.as_f64().unwrap(),
    }
}

fn at(
    value: &Value,
    index: usize,
) -> f64 {
    let Value::Tensor(tensor) = value else {
        panic!("expected tensor: {value:?}")
    };
    number(&tensor.data[index])
}

fn length_rate(order: i64) -> DimExpr {
    DimExpr::Mul(vec![length(), time().pow(order)])
}

fn vector_of(dim: DimExpr) -> TypeExpr {
    TypeExpr::Tensor {
        rank: NatExpr::Lit(1),
        shape: None,
        elem: Box::new(TypeExpr::Float.with_dim(dim)),
    }
}

fn state_signature() -> TypeExpr {
    TypeExpr::Record(vec![
        ("displacement".into(), vector_of(length())),
        ("velocity".into(), vector_of(length_rate(-1))),
        ("acceleration".into(), vector_of(length_rate(-2))),
        ("time".into(), TypeExpr::Float.with_dim(time())),
    ])
}

fn release(
    session: &Session,
    value: &Value,
) -> Result<Value, OpError> {
    session.call_value(
        "num_backend_release",
        std::slice::from_ref(value),
        &CallCtx::new(),
    )
}

#[test]
fn hyperelastic_programs_remain_charged_while_measured_states_share_the_model() {
    for name in ["auto", "interpreter", "cpu", "cpu-parallel", "cranelift"] {
        let Ok(backend) = rssn::backend::named(name) else {
            assert_eq!(name, "cranelift", "required backend is unavailable");
            continue;
        };
        rssn::backend::with_backend(backend, check_hyperelastic_program_retention);
    }
}

fn check_hyperelastic_program_retention() {
    let s = Session::new();
    let refs = refs(&s);
    let frame = call(
        &s,
        "num_grid_frame_new",
        &[
            measured(2., length()),
            measured(4., time()),
            measured(2., length()),
        ],
    );
    let mesh = call(
        &s,
        "num_grid_simplex_rectangle",
        &[
            frame,
            Value::int(1),
            Value::int(1),
            floats(&[0., 0.]),
            floats(&[1., 1.]),
        ],
    );
    let material = call(
        &s,
        "num_grid_hyper_material",
        &[
            refs.clone(),
            Value::String("NeoHookean".into()),
            Value::Float(10.),
            Value::Float(0.25),
            Value::Float(1.),
        ],
    );
    let bcs = call(
        &s,
        "num_grid_elastic_bcs",
        &[
            refs,
            vector(vec![], TypeExpr::Int),
            opaque(vec![], "rssn/native/elastic-boundary/v1"),
        ],
    );
    let force_dim = DimExpr::Mul(vec![
        DimExpr::base("Mass"),
        length().pow(-2),
        time().pow(-2),
    ]);
    let force_x = program(&s, "1", 0.75, force_dim.clone());
    let force_y = program(&s, "0", 1., force_dim.clone());
    let force_z = program(&s, "0", 1., force_dim);
    let problem = call(
        &s,
        "num_grid_hyper_problem",
        &[
            mesh,
            material,
            bcs,
            opaque(
                vec![force_x.clone(), force_y.clone(), force_z.clone()],
                COMPILED,
            ),
            opaque(vec![], CONTACT),
            Value::int(1024),
        ],
    );
    let rest = call(
        &s,
        "num_grid_hyper_rest",
        &[problem.clone(), floats(&[0.; 8])],
    );
    assert!(
        release(&s, &force_x)
            .unwrap_err()
            .to_string()
            .contains("dependent handles")
    );
    assert!(
        release(&s, &problem)
            .unwrap_err()
            .to_string()
            .contains("dependent handles")
    );

    let times = Arc::new(Mutex::new(Vec::new()));
    let observed = times.clone();
    let expected = state_signature();
    let mut ctx = CallCtx::new();
    let observer = ctx
        .add_callback(Arc::new(NativeCallback::new(
            move |args: &[Value]| {
                assert!(matches_type(&expected, &args[0].type_expr()).accepted());
                observed
                    .lock()
                    .unwrap()
                    .push(number(field(&args[0], "time")));
                Ok(Value::Bool(true))
            },
            Signature::new(vec![state_signature()], TypeExpr::Bool),
        )))
        .unwrap();
    let end = s
        .call_value(
            "num_grid_hyper_dynamic",
            &[
                rest.clone(),
                Value::String("Newmark".into()),
                Value::Float(0.25),
                Value::Float(0.5),
                Value::Float(1.),
                Value::Float(0.125),
                Value::int(2),
                Value::Float(1e-10),
                Value::Function(observer),
                Value::int(2),
            ],
            &ctx,
        )
        .unwrap();
    assert_eq!(*times.lock().unwrap(), [0.5, 1.]);
    let payload = call(
        &s,
        "num_grid_hyper_state_fields",
        std::slice::from_ref(&end),
    );
    assert_eq!(number(field(&payload, "time")), 1.);
    for (key, dim, expected_reference, expected_x) in [
        ("displacement", length(), 2., 1.),
        ("velocity", length_rate(-1), 0.5, 2.),
        ("acceleration", length_rate(-2), 0.125, 2.),
    ] {
        let stored = field(&payload, key);
        let reference = field(stored, "reference");
        assert!(matches_type(&TypeExpr::Float.with_dim(dim), &reference.type_expr()).accepted());
        assert_eq!(number(reference), expected_reference);
        let values = field(stored, "normalized");
        assert!(matches_type(&vector_of(DimExpr::Dimless), &values.type_expr()).accepted());
        for i in 0..4 {
            let x = at(values, 2 * i) * expected_reference;
            let y = at(values, 2 * i + 1) * expected_reference;
            assert!((x - expected_x).abs() < 1e-9, "{key} node {i}");
            assert!(y.abs() < 1e-9, "{key} node {i} y");
        }
    }
    assert_eq!(release(&s, &rest).unwrap(), Value::Bool(true));
    assert!(
        release(&s, &problem).is_err(),
        "the new state independently retains the shared model"
    );
    for handle in [&end, &problem, &force_x, &force_y, &force_z] {
        assert_eq!(release(&s, handle).unwrap(), Value::Bool(true));
    }
}

#[test]
fn penalty_contact_keeps_the_payload_cap_and_releases_parents_after_rejection() {
    let backend = rssn::backend::named("interpreter").unwrap();
    rssn::backend::with_backend(backend, check_penalty_contact_payload_cap);
}

fn check_penalty_contact_payload_cap() {
    let s = Session::new();
    let refs = refs(&s);
    let gap = program(&s, "1", 1., length());
    let normal = program(&s, "1", 1., DimExpr::Dimless);
    let contact = call(
        &s,
        "num_grid_hyper_contact",
        &[
            refs.clone(),
            Value::int(1),
            gap.clone(),
            opaque(vec![normal.clone(); 3], COMPILED),
            Value::Float(1.),
            Value::int(10),
        ],
    );
    assert!(
        release(&s, &normal).is_err(),
        "normal program is retained by the contact"
    );
    assert_eq!(release(&s, &contact).unwrap(), Value::Bool(true));

    let wide_expression = (1..=100)
        .map(|i| format!("sin({i}*x)"))
        .collect::<Vec<_>>()
        .join("+");
    let too_large = program(&s, &wide_expression, 1., DimExpr::Dimless);
    let error = s
        .call_value(
            "num_grid_hyper_contact",
            &[
                refs,
                Value::int(1),
                gap.clone(),
                opaque(vec![too_large.clone(); 3], COMPILED),
                Value::Float(1.),
                Value::int(10),
            ],
            &CallCtx::new(),
        )
        .unwrap_err();
    assert!(
        error
            .to_string()
            .contains("penalty contact compiled storage cap exceeded"),
        "{error}"
    );
    for handle in [&too_large, &gap, &normal] {
        assert_eq!(
            release(&s, handle).unwrap(),
            Value::Bool(true),
            "failed constructor kept a parent"
        );
    }
}
