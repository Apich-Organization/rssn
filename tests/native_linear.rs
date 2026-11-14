//! Native tensor and block-sparse adapters preserve explicit reference scales.
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
fn ints(a: &[usize]) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![a.len()],
        elem: TypeExpr::Int,
        data: a.iter().map(|v| Value::int(*v)).collect(),
    })
}
fn floats(a: &[f64]) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![a.len()],
        elem: TypeExpr::Float,
        data: a.iter().copied().map(Value::Float).collect(),
    })
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
fn scalar(v: &Value) -> f64 {
    if let Value::Dyn(d) = v {
        scalar(&d.value)
    } else {
        v.as_f64().unwrap()
    }
}
fn array(
    s: &Session,
    shape: &[usize],
    data: &[f64],
    reference: Value,
) -> Value {
    call(
        s,
        "num_linear_array_new",
        &[ints(shape), floats(data), reference],
    )
    .unwrap()
}
#[test]
fn train_native_entry_roundtrip_and_inner_product_retain_units() {
    let s = Session::new();
    let a = array(&s, &[2, 2], &[1.0, 2.0, 2.0, 4.0], measured(3.0, "Length"));
    let train = call(
        &s,
        "num_linear_tt_svd",
        &[a.clone(), Value::Float(1e-12), Value::int(2)],
    )
    .unwrap();
    let entry = call(
        &s,
        "num_linear_train_entry",
        &[train.clone(), ints(&[1, 1])],
    )
    .unwrap();
    assert!((scalar(&entry) - 12.0).abs() < 1e-10);
    assert!(matches!(entry, Value::Dyn(_)));
    let restored = call(&s, "num_linear_train_to_dense", &[train.clone()]).unwrap();
    let raw = call(&s, "num_linear_array_data", &[restored]).unwrap();
    let Value::Tensor(t) = field(&raw, "data") else {
        panic!("tensor")
    };
    assert!((scalar(&t.data[3]) - 4.0).abs() < 1e-10);
    assert_eq!(scalar(field(&raw, "reference")), 3.0);
    let inner = call(&s, "num_linear_train_inner", &[train.clone(), train]).unwrap();
    assert!((scalar(&inner) - 225.0).abs() < 1e-8);
    let other = array(&s, &[2, 2], &[1.0; 4], measured(3.0, "Time"));
    let b = call(
        &s,
        "num_linear_tt_svd",
        &[other, Value::Float(1e-12), Value::int(2)],
    )
    .unwrap();
    let first = call(
        &s,
        "num_linear_tt_svd",
        &[a, Value::Float(1e-12), Value::int(2)],
    )
    .unwrap();
    assert!(call(&s, "num_linear_train_add", &[first, b]).is_err());
}
#[test]
fn bsr_native_transpose_product_and_input_guards_are_real() {
    let s = Session::new();
    let a = call(
        &s,
        "num_linear_csr_new",
        &[
            Value::int(2),
            Value::int(2),
            ints(&[0, 0, 1]),
            ints(&[0, 1, 1]),
            floats(&[2.0, 3.0, 4.0]),
            measured(2.0, "Time"),
        ],
    )
    .unwrap();
    let b = call(
        &s,
        "num_linear_bsr_from_csr",
        &[a, Value::int(1), Value::int(1)],
    )
    .unwrap();
    let x = array(&s, &[2], &[1.0, 2.0], measured(3.0, "Length"));
    let y = call(&s, "num_linear_bsr_matvec_transpose", &[b.clone(), x]).unwrap();
    let y = call(&s, "num_linear_array_data", &[y]).unwrap();
    let Value::Tensor(t) = field(&y, "data") else {
        panic!("tensor")
    };
    assert_eq!(scalar(&t.data[0]), 2.0);
    assert_eq!(scalar(&t.data[1]), 11.0);
    assert_eq!(scalar(field(&y, "reference")), 6.0);
    assert!(
        call(
            &s,
            "num_linear_bsr_from_csr",
            &[b, Value::int(1), Value::int(1)]
        )
        .is_err()
    );
    assert!(
        call(
            &s,
            "num_linear_array_new",
            &[ints(&[usize::MAX, 2]), floats(&[]), Value::Float(1.0)]
        )
        .is_err()
    );
    let foreign = Session::new();
    let local = array(&s, &[1], &[1.0], Value::Float(1.0));
    assert!(call(&foreign, "num_linear_array_norm", &[local]).is_err());
}
