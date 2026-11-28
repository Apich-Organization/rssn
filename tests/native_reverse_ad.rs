//! Stateful reverse AD uses native tapes and coherent physical dimensions.
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::value::DynValue;
fn call(
    s: &Session,
    name: &str,
    args: &[Value],
) -> Result<Value, OpError> {
    s.call_value(name, args, &CallCtx::new())
}
fn scalar(v: Value) -> f64 {
    match v {
        | Value::Dyn(d) => scalar(*d.value),
        | value => value.as_f64().unwrap(),
    }
}
#[test]
fn native_reverse_tape_propagates_physical_product_and_gradient_dimensions() {
    let s = Session::new();
    let tape = call(&s, "num_ad_tape_new", &[Value::int(32)]).unwrap();
    let physical = Value::Dyn(DynValue {
        value: Box::new(Value::Float(3.0)),
        dim: DimExpr::base("Length"),
    });
    let x = call(&s, "num_ad_tape_var", &[tape.clone(), physical]).unwrap();
    let square = call(&s, "num_ad_var_mul", &[x.clone(), x.clone()]).unwrap();
    let output = call(&s, "num_ad_var_value", std::slice::from_ref(&square)).unwrap();
    let Value::Dyn(d) = &output else {
        panic!("squared length witness")
    };
    assert_eq!(d.dim, DimExpr::base("Length").pow(2));
    assert_eq!(scalar(output), 9.0);
    let g = call(&s, "num_ad_var_backward", std::slice::from_ref(&square)).unwrap();
    let derivative = call(&s, "num_ad_var_grad_in", &[x.clone(), g.clone()]).unwrap();
    let Value::Dyn(d) = &derivative else {
        panic!("derivative length witness")
    };
    assert_eq!(d.dim, DimExpr::base("Length"));
    assert_eq!(scalar(derivative), 6.0);
    let later = call(&s, "num_ad_tape_var", &[tape.clone(), Value::Float(7.0)]).unwrap();
    assert_eq!(
        scalar(call(&s, "num_ad_var_grad_in", &[x.clone(), g.clone()]).unwrap()),
        6.0
    );
    assert_eq!(
        scalar(call(&s, "num_ad_var_grad_in", &[later.clone(), g.clone()]).unwrap()),
        0.0
    );
    call(&s, "num_native_release", &[later]).unwrap();
    assert!(call(&s, "num_ad_var_exp", std::slice::from_ref(&x)).is_err());
    assert!(call(&s, "num_native_release", std::slice::from_ref(&tape)).is_err());
    let foreign = Session::new();
    assert!(call(&foreign, "num_ad_var_value", std::slice::from_ref(&x)).is_err());
    for handle in [g, square, x, tape] {
        call(&s, "num_native_release", std::slice::from_ref(&handle)).unwrap();
        assert!(call(&s, "num_native_release", &[handle]).is_err());
    }
}
#[test]
fn native_reverse_rejects_foreign_tapes_unequal_dimensions_and_exhaustion() {
    let s = Session::new();
    let tape = call(&s, "num_ad_tape_new", &[Value::int(2)]).unwrap();
    let other = call(&s, "num_ad_tape_new", &[Value::int(8)]).unwrap();
    let x = call(&s, "num_ad_tape_var", &[tape.clone(), Value::Float(2.0)]).unwrap();
    let y = call(&s, "num_ad_tape_var", &[other.clone(), Value::Float(3.0)]).unwrap();
    assert!(call(&s, "num_ad_var_mul", &[x.clone(), y]).is_err());
    let z = call(
        &s,
        "num_ad_tape_var",
        &[
            other,
            Value::Dyn(DynValue {
                value: Box::new(Value::Float(3.0)),
                dim: DimExpr::base("Time"),
            }),
        ],
    )
    .unwrap();
    assert!(call(&s, "num_ad_var_add", &[x.clone(), z]).is_err());
    let square = call(&s, "num_ad_var_mul", &[x.clone(), x.clone()]).unwrap();
    assert_eq!(
        scalar(call(&s, "num_ad_var_value", &[square]).unwrap()),
        4.0
    );
    assert!(call(&s, "num_ad_var_sin", &[x]).is_err());
    assert!(call(&s, "num_ad_tape_var", &[tape, Value::Float(4.0)]).is_err());
}

#[test]
fn native_ad_division_keeps_representable_derivatives_at_extreme_scales() {
    use rssn::kernels::autodiff::Dual;
    let s = Session::new();
    for magnitude in [1e308, 1e-200] {
        let tape = call(&s, "num_ad_tape_new", &[Value::int(8)]).unwrap();
        let x = call(
            &s,
            "num_ad_tape_var",
            &[tape.clone(), Value::Float(magnitude)],
        )
        .unwrap();
        let y = call(&s, "num_ad_tape_var", &[tape, Value::Float(magnitude)]).unwrap();
        let q = call(&s, "num_ad_var_div", &[x.clone(), y.clone()]).unwrap();
        let g = call(&s, "num_ad_var_backward", &[q]).unwrap();
        let expected = -1.0 / magnitude;
        let got = scalar(call(&s, "num_ad_var_grad_in", &[y, g]).unwrap());
        assert!(
            (got / expected - 1.0).abs() < 1e-14,
            "{magnitude}: {got} != {expected}"
        );
        let tangent = Dual::constant(magnitude) / Dual::var(magnitude);
        assert!((tangent.deriv / expected - 1.0).abs() < 1e-14);
        let tangent = Dual::var(magnitude) / Dual::constant(magnitude);
        assert!((tangent.deriv * magnitude - 1.0).abs() < 1e-14);
    }
}

#[test]
fn rejected_native_ad_operations_do_not_consume_persistent_tape_nodes() {
    let s = Session::new();
    let tape = call(&s, "num_ad_tape_new", &[Value::int(2)]).unwrap();
    let x = call(&s, "num_ad_tape_var", &[tape.clone(), Value::Float(-1.0)]).unwrap();
    assert!(call(&s, "num_ad_var_sqrt", std::slice::from_ref(&x)).is_err());
    assert_eq!(
        call(&s, "num_ad_tape_len", std::slice::from_ref(&tape)).unwrap(),
        Value::int(1)
    );
    assert!(call(&s, "num_ad_var_sin", &[x]).is_ok());
    assert_eq!(call(&s, "num_ad_tape_len", &[tape]).unwrap(), Value::int(2));
    let full = Session::new();
    let tape = call(&full, "num_ad_tape_new", &[Value::int(5000)]).unwrap();
    for _ in 0..4095 {
        call(&full, "num_ad_tape_var", &[tape.clone(), Value::Float(1.0)]).unwrap();
    }
    assert!(call(&full, "num_ad_tape_var", &[tape.clone(), Value::Float(1.0)]).is_err());
    assert_eq!(
        call(&full, "num_ad_tape_len", &[tape]).unwrap(),
        Value::int(4095)
    );
}
