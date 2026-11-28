//! Mixed translational/rotational finite-element references are not flattened to one unit.
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
    name: &str,
) -> Value {
    Value::Dyn(DynValue {
        value: Box::new(Value::Float(x)),
        dim: DimExpr::Base(name.into()),
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
fn vector(a: &[f64]) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![a.len()],
        elem: TypeExpr::Float,
        data: a.iter().copied().map(Value::Float).collect(),
    })
}
fn frame(s: &Session) -> Value {
    call(
        s,
        "num_grid_mechanical_references",
        &[
            measured(2.0, "Length"),
            measured(1.0, "Time"),
            measured(3.0, "Mass"),
            measured(10.0, "Temperature"),
        ],
    )
}
#[test]
fn beam_matrix_uses_energy_and_distinct_dof_references() {
    let s = Session::new();
    let f = frame(&s);
    let beam = call(
        &s,
        "num_grid_beam_new",
        &[
            f,
            Value::Float(2.0),
            Value::Float(5.0),
            Value::Float(3.0),
            Value::Float(4.0),
            Value::Float(0.0),
        ],
    );
    let k = call(&s, "num_grid_beam_local_stiffness_matrix", &[beam]);
    let Value::Tensor(data) = field(&k, "data") else {
        panic!("matrix")
    };
    let energy = number(field(&k, "quadratic_form_reference"));
    assert_eq!(energy, 12.0);
    let dofs = field(&k, "dof_references");
    assert_eq!(number(field(dofs, "dof0")), 2.0);
    assert_eq!(number(field(dofs, "dof2")), 1.0);
    assert!(matches!(field(dofs, "dof0"), Value::Dyn(_)));
    assert!(matches!(field(dofs, "dof2"), Value::Float(_)));
    assert!((number(&data.data[0]) * energy / 4.0 - 22.5).abs() < 1e-12);
    assert!((number(&data.data[14]) * energy - 480.0).abs() < 1e-12);
}
#[test]
fn principal_angle_keeps_radians_while_stresses_restore_pascal_reference() {
    let s = Session::new();
    let f = frame(&s);
    let p = call(
        &s,
        "num_grid_principal_stresses",
        &[f.clone(), vector(&[3.0, 1.0, 1.0])],
    );
    assert!((number(field(&p, "sigma1")) - (2.0 + 2.0f64.sqrt()) * 1.5).abs() < 1e-12);
    assert!((number(field(&p, "sigma2")) - (2.0 - 2.0f64.sqrt()) * 1.5).abs() < 1e-12);
    assert!((number(field(&p, "angle")) - std::f64::consts::FRAC_PI_8).abs() < 1e-12);
    assert!(matches!(field(&p, "angle"), Value::Float(_)));
    assert!(
        s.call_value(
            "num_grid_principal_stresses",
            &[f, vector(&[0.0; 6])],
            &CallCtx::new()
        )
        .is_err()
    );
}
