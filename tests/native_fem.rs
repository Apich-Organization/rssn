//! Native P2 finite elements preserve diffusion/source/field references.
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::TensorValue;
fn call(
    s: &Session,
    n: &str,
    a: &[Value],
) -> Value {
    s.call_value(n, a, &CallCtx::new()).unwrap()
}
fn measured(
    x: f64,
    d: DimExpr,
) -> Value {
    Value::Dyn(DynValue {
        value: Box::new(Value::Float(x)),
        dim: d,
    })
}
fn base(n: &str) -> DimExpr {
    DimExpr::Base(n.into())
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
fn field<'a>(
    v: &'a Value,
    n: &str,
) -> &'a Value {
    let Value::Record(r) = v else {
        panic!("record")
    };
    &r.fields.iter().find(|(k, _)| k == n).unwrap().1
}
fn number(v: &Value) -> f64 {
    if let Value::Dyn(d) = v {
        number(&d.value)
    } else {
        v.as_f64().unwrap()
    }
}
fn compile(s: &Session) -> Value {
    let backend = call(s, "num_backend_interpreter", &[]);
    let references = [
        measured(1.0, base("Time")),
        measured(1.0, base("Length")),
        measured(1.0, base("Length")),
        measured(1.0, base("Length")),
    ]
    .into_iter()
    .map(|v| call(s, "num_backend_reference", &[v]))
    .collect();
    let output = call(
        s,
        "num_backend_reference",
        &[measured(
            1.5,
            DimExpr::Mul(vec![
                base("Length"),
                DimExpr::Inverse(Box::new(base("Time"))),
            ]),
        )],
    );
    call(
        s,
        "num_backend_compile",
        &[
            backend,
            Value::String("2".into()),
            vector(
                ["t", "x", "y", "z"]
                    .into_iter()
                    .map(|n| Value::String(n.into()))
                    .collect(),
                TypeExpr::String,
            ),
            vector(
                references,
                TypeExpr::Opaque("rssn/native/backend-scalar/v1".into()),
            ),
            output,
            Value::int(100000),
        ],
    )
}
#[test]
fn native_p2_static_solution_and_physical_evaluation_match_quadratic() {
    let s = Session::new();
    let frame = call(
        &s,
        "num_grid_frame_new",
        &[
            measured(2.0, base("Length")),
            measured(2.0, base("Time")),
            measured(3.0, base("Length")),
        ],
    );
    let boundary = call(
        &s,
        "num_grid_dirichlet",
        &[frame.clone(), Value::Float(0.0)],
    );
    let markers = call(&s, "num_grid_marker_with_default", &[boundary]);
    let source = compile(&s);
    let problem = call(
        &s,
        "num_grid_fem_constant_problem",
        &[
            markers,
            Value::Float(1.0),
            Value::Float(0.0),
            source.clone(),
            Value::int(10000),
        ],
    );
    assert!(
        s.call_value("num_backend_release", &[source], &CallCtx::new())
            .is_err()
    );
    let mesh = call(
        &s,
        "num_grid_simplex_interval",
        &[frame, Value::int(1), Value::Float(0.0), Value::Float(1.0)],
    );
    assert!(
        (number(&call(
            &s,
            "num_grid_simplex_measure",
            &[mesh.clone(), Value::int(0)]
        )) - 2.0)
            .abs()
            < 1e-12
    );
    let space = call(&s, "num_grid_fem_space_new", &[mesh, Value::int(2)]);
    let solution = call(
        &s,
        "num_grid_fem_solve_static",
        &[
            space.clone(),
            problem,
            Value::Float(0.0),
            Value::String("Auto".into()),
            Value::Float(1e-12),
        ],
    );
    let Value::Tensor(data) = field(&solution, "data") else {
        panic!("field")
    };
    let values = vector(data.data.clone(), TypeExpr::Float);
    let point = vector(
        vec![Value::Float(1.0), Value::Float(0.0), Value::Float(0.0)],
        TypeExpr::Float.with_dim(base("Length")),
    );
    let evaluated = call(&s, "num_grid_fem_eval", &[space, values, point]);
    let Value::Variant(v) = evaluated else {
        panic!("optional evaluation")
    };
    assert_eq!(v.tag, "Some");
    assert!((number(&v.payload[0]) - 0.75).abs() < 1e-12);
}
