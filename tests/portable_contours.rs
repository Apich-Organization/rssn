//! Checked contour geometry and physical coordinates across the native registry.

use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::Value;
use rssn::kernel::overload::matches_type;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::ComplexValue;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::TensorValue;

fn measured(
    value: Value,
    dim: DimExpr,
) -> Value {
    Value::Dyn(DynValue {
        value: Box::new(value),
        dim,
    })
}

fn length(x: f64) -> Value {
    measured(Value::Float(x), DimExpr::base("Length"))
}

fn point(
    re: f64,
    im: f64,
) -> Value {
    measured(
        Value::Complex(ComplexValue::Float { re, im }),
        DimExpr::base("Length"),
    )
}

fn poles(values: &[f64]) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![values.len()],
        elem: TypeExpr::Float.with_dim(DimExpr::base("Length")),
        data: values.iter().copied().map(Value::Float).collect(),
    })
}

fn field<'a>(
    value: &'a Value,
    name: &str,
) -> &'a Value {
    let Value::Record(record) = value else {
        panic!("expected contour")
    };
    &record.fields.iter().find(|(key, _)| key == name).unwrap().1
}

fn tensor(
    value: &Value,
    size: usize,
    ty: TypeExpr,
) -> &[Value] {
    let Value::Tensor(tensor) = value else {
        panic!("expected contour vector")
    };
    assert_eq!(tensor.shape, [size]);
    assert!(matches_type(&ty, &tensor.elem).accepted());
    &tensor.data
}

fn payload(value: &Value) -> &Value {
    match value {
        | Value::Dyn(v) => payload(&v.value),
        | value => value,
    }
}

fn complex(value: &Value) -> (f64, f64) {
    let Value::Complex(c) = payload(value) else {
        panic!("expected complex coordinate")
    };
    c.to_f64()
}

fn call(
    session: &Session,
    name: &str,
    args: &[Value],
) -> Value {
    session
        .call_value(&format!("num_contour_{name}"), args, &CallCtx::new())
        .unwrap()
}

#[test]
fn measured_paths_preserve_direction_arc_metadata_and_open_state() {
    let session = Session::new();
    let coord = TypeExpr::Complex.with_dim(DimExpr::base("Length"));
    let c = call(&session, "bromwich", &[length(2.0), length(3.0)]);
    assert_eq!(field(&c, "closed"), &Value::Bool(false));
    assert_eq!(
        complex(&tensor(field(&c, "starts"), 1, coord.clone())[0]),
        (2.0, -3.0)
    );
    assert_eq!(
        complex(&tensor(field(&c, "ends"), 1, coord.clone())[0]),
        (2.0, 3.0)
    );
    assert!(tensor(field(&c, "arc_centers"), 0, coord.clone()).is_empty());
    let c = call(
        &session,
        "rectangle",
        &[length(-1.0), length(1.0), length(-2.0), length(2.0)],
    );
    assert_eq!(field(&c, "closed"), &Value::Bool(true));
    let starts = tensor(field(&c, "starts"), 4, coord.clone());
    assert_eq!(
        starts.iter().map(complex).collect::<Vec<_>>(),
        [(-1.0, -2.0), (1.0, -2.0), (1.0, 2.0), (-1.0, 2.0)]
    );
    let c = call(
        &session,
        "sector",
        &[length(2.0), Value::Float(-std::f64::consts::FRAC_PI_2)],
    );
    assert_eq!(
        tensor(field(&c, "arc_segments"), 1, TypeExpr::Int),
        &[Value::int(1)]
    );
    let last = &tensor(field(&c, "ends"), 3, coord)[1];
    let (re, im) = complex(last);
    assert!(re.abs() < 1e-14 && im == -2.0);
    assert_eq!(
        tensor(field(&c, "theta_end"), 1, TypeExpr::Float),
        &[Value::Float(-std::f64::consts::FRAC_PI_2)]
    );
}

#[test]
fn cut_and_indent_paths_keep_all_arcs_in_normalized_geometry() {
    let s = Session::new();
    let lens = TypeExpr::Float.with_dim(DimExpr::base("Length"));
    let key = call(
        &s,
        "keyhole",
        &[
            point(1.0, 1.0),
            length(1.0),
            length(3.0),
            Value::Float(-std::f64::consts::PI),
            Value::Float(0.15),
        ],
    );
    assert_eq!(
        tensor(field(&key, "arc_segments"), 2, TypeExpr::Int),
        &[Value::int(0), Value::int(2)]
    );
    assert_eq!(field(&key, "closed"), &Value::Bool(true));
    let h = call(
        &s,
        "hankel",
        &[
            length(1.0),
            length(3.0),
            Value::Float(-std::f64::consts::PI),
            Value::Float(0.15),
        ],
    );
    assert_eq!(field(&h, "closed"), &Value::Bool(false));
    assert_eq!(
        tensor(field(&h, "arc_segments"), 1, TypeExpr::Int),
        &[Value::int(1)]
    );
    let d = call(
        &s,
        "dogbone",
        &[point(1.0, 1.0), point(5.0, 2.0), length(0.5)],
    );
    assert_eq!(
        tensor(field(&d, "arc_segments"), 2, TypeExpr::Int),
        &[Value::int(1), Value::int(3)]
    );
    for scale in [1.0, 1e-80, 1e80] {
        let c = call(
            &s,
            "indented_semicircle",
            &[
                length(2.0 * scale),
                poles(&[0.5 * scale, -0.5 * scale]),
                length(0.125 * scale),
                Value::Bool(true),
            ],
        );
        let rs = tensor(field(&c, "arc_radii"), 3, lens.clone());
        let Value::Float(value) = payload(&rs[0]) else {
            panic!("expected real radius")
        };
        assert!((*value / scale - 0.125).abs() < 1e-13);
        assert_eq!(field(&c, "closed"), &Value::Bool(true));
    }
}

#[test]
fn invalid_geometry_and_dimension_mismatch_are_errors() {
    let s = Session::new();
    let ctx = CallCtx::new();
    assert!(
        s.call_value(
            "num_contour_bromwich",
            &[
                length(2.0),
                measured(Value::Float(1.0), DimExpr::base("Time"))
            ],
            &ctx
        )
        .is_err()
    );
    assert!(
        s.call_value(
            "num_contour_rectangle",
            &[length(2.0), length(1.0), length(-1.0), length(1.0)],
            &ctx
        )
        .is_err()
    );
    assert!(
        s.call_value(
            "num_contour_keyhole",
            &[
                point(0.0, 0.0),
                length(1.0),
                length(2.0),
                Value::Float(0.0),
                Value::Float(f64::EPSILON)
            ],
            &ctx
        )
        .is_err()
    );
    assert!(
        s.call_value(
            "num_contour_indented_semicircle",
            &[
                length(2.0),
                poles(&[1.99]),
                length(0.125),
                Value::Bool(true)
            ],
            &ctx
        )
        .is_err()
    );
    assert!(
        s.call_value(
            "num_contour_indented_semicircle",
            &[
                length(2.0),
                poles(&[0.1, 0.1]),
                length(0.125),
                Value::Bool(true)
            ],
            &ctx
        )
        .is_err()
    );
}
