//! Native grid callbacks retain coherent reference scales and their compiled parents.
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::TensorValue;
fn call(
    s: &Session,
    n: &str,
    a: &[Value],
) -> Result<Value, OpError> {
    s.call_value(n, a, &CallCtx::new())
}
fn vec(
    v: Vec<Value>,
    elem: TypeExpr,
) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![v.len()],
        data: v,
        elem,
    })
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
fn frame(s: &Session) -> Value {
    call(
        s,
        "num_grid_frame_new",
        &[
            measured(2.0, DimExpr::Base("Length".into())),
            measured(2.0, DimExpr::Base("Time".into())),
            measured(3.0, DimExpr::Base("Temperature".into())),
        ],
    )
    .unwrap()
}
fn compile(
    s: &Session,
    expr: &str,
    names: &[&str],
    refs: Vec<Value>,
    output: Value,
) -> Value {
    let backend = call(s, "num_backend_interpreter", &[]).unwrap();
    let refs = refs
        .into_iter()
        .map(|v| call(s, "num_backend_reference", &[v]).unwrap())
        .collect();
    let output = call(s, "num_backend_reference", &[output]).unwrap();
    call(
        s,
        "num_backend_compile",
        &[
            backend,
            Value::String(expr.into()),
            vec(
                names.iter().map(|s| Value::String((*s).into())).collect(),
                TypeExpr::String,
            ),
            vec(
                refs,
                TypeExpr::Opaque("rssn/native/backend-scalar/v1".into()),
            ),
            output,
            Value::int(10000),
        ],
    )
    .unwrap()
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
#[test]
fn retained_time_boundary_evaluates_coherent_seconds_and_holds_parent() {
    let s = Session::new();
    let f = frame(&s);
    let program = compile(
        &s,
        "1+t",
        &["t"],
        vec![measured(1.0, DimExpr::Base("Time".into()))],
        measured(3.0, DimExpr::Base("Temperature".into())),
    );
    let datum = call(
        &s,
        "num_grid_datum_time",
        &[
            f,
            program.clone(),
            Value::int(2),
            measured(3.0, DimExpr::Base("Temperature".into())),
        ],
    )
    .unwrap();
    assert!(call(&s, "num_backend_release", std::slice::from_ref(&program)).is_err());
    let point = vec(vec![Value::Float(0.0); 3], TypeExpr::Float);
    let value = call(
        &s,
        "num_grid_datum_eval",
        &[datum.clone(), Value::Float(0.5), point.clone()],
    )
    .unwrap();
    assert_eq!(number(&value), 6.0);
    assert!(matches!(value, Value::Dyn(_)));
    call(
        &s,
        "num_grid_datum_eval",
        &[datum.clone(), Value::Float(0.0), point.clone()],
    )
    .unwrap();
    assert!(
        call(
            &s,
            "num_grid_datum_eval",
            &[datum.clone(), Value::Float(0.0), point]
        )
        .is_err()
    );
    call(&s, "num_backend_release", &[datum]).unwrap();
    assert!(call(&s, "num_backend_release", &[program]).is_ok());
}
#[test]
fn native_helmholtz_source_keeps_physical_laplacian_units_and_exact_grid_profile() {
    let s = Session::new();
    let f = frame(&s);
    let b = call(&s, "num_grid_faces_dirichlet_zero", &[f]).unwrap();
    let g = call(
        &s,
        "num_grid_structured_unit",
        &[b, Value::int(1), Value::int(9)],
    )
    .unwrap();
    let forcing = DimExpr::Mul(vec![
        DimExpr::Base("Temperature".into()),
        DimExpr::Inverse(Box::new(DimExpr::Mul(vec![
            DimExpr::Base("Length".into()),
            DimExpr::Base("Length".into()),
        ]))),
    ]);
    let program = compile(
        &s,
        "1",
        &["x", "y", "z"],
        vec![measured(1.0, DimExpr::Base("Length".into())); 3],
        measured(0.75, forcing),
    );
    let out = call(
        &s,
        "num_grid_solve_helmholtz",
        &[
            g,
            Value::Float(0.0),
            program,
            Value::int(100),
            Value::String("Direct".into()),
            Value::Float(1e-12),
        ],
    )
    .unwrap();
    let Value::Tensor(data) = field(&out, "data") else {
        panic!("data")
    };
    assert!((number(&data.data[4]) - 0.125).abs() < 1e-12);
    assert_eq!(number(field(&out, "value_reference")), 3.0);
}
