//! Prepared weighted sampling and session-owned random streams.

use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::TensorValue;

const NEW: &str = "num_native_categorical_new";
const CATEGORIES: &str = "num_native_categorical_categories";
const SAMPLE: &str = "num_native_categorical_sample";
const FILL: &str = "num_native_categorical_fill";
const NEXT: &str = "num_native_rng_next_u64";

fn call(
    s: &Session,
    name: &str,
    args: &[Value],
) -> Result<Value, OpError> {
    s.call_value(name, args, &CallCtx::new())
}

fn weights(values: &[f64]) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![values.len()],
        elem: TypeExpr::Float,
        data: values.iter().copied().map(Value::Float).collect(),
    })
}

fn rng(
    s: &Session,
    seed: u64,
) -> Value {
    call(s, "num_native_rng_new", &[Value::int(seed)]).unwrap()
}

#[test]
fn batch_draws_match_scalar_sampling_and_exclude_zero_weight_categories() {
    let s = Session::new();
    let distribution = call(&s, NEW, &[weights(&[f64::MAX, 0.0, f64::MAX])]).unwrap();
    let batch_rng = rng(&s, 42);
    let scalar_rng = rng(&s, 42);
    assert_eq!(
        call(&s, CATEGORIES, std::slice::from_ref(&distribution)).unwrap(),
        Value::int(3)
    );
    let Value::Tensor(batch) = call(
        &s,
        FILL,
        &[distribution.clone(), batch_rng.clone(), Value::int(128)],
    )
    .unwrap() else {
        panic!("expected a tensor of indices");
    };
    assert_eq!(batch.shape, [128]);
    assert_eq!(batch.elem, TypeExpr::Int);
    assert!(batch.data.contains(&Value::int(0)));
    assert!(batch.data.contains(&Value::int(2)));
    for actual in batch.data {
        let expected = call(&s, SAMPLE, &[distribution.clone(), scalar_rng.clone()]).unwrap();
        assert!([Value::int(0), Value::int(2)].contains(&actual));
        assert_eq!(actual, expected);
    }
    assert_eq!(
        call(&s, NEXT, &[batch_rng]).unwrap(),
        call(&s, NEXT, &[scalar_rng]).unwrap()
    );
}

#[test]
fn invalid_and_degenerate_draws_leave_the_stream_unchanged() {
    let s = Session::new();
    let foreign = Session::new();
    for invalid in [
        weights(&[]),
        weights(&[0.0, 0.0]),
        weights(&[-1.0, 2.0]),
        weights(&[f64::NAN]),
    ] {
        assert!(call(&s, NEW, &[invalid]).is_err());
    }
    let mut measured = weights(&[2.0, 1.0]);
    let Value::Tensor(tensor) = &mut measured else {
        unreachable!()
    };
    let length = DimExpr::Base("Length".into());
    tensor.elem = TypeExpr::Float.with_dim(length.clone());
    tensor.data[0] = Value::Dyn(DynValue {
        dim: length,
        value: Box::new(Value::Float(2.0)),
    });
    assert!(call(&s, NEW, &[measured]).is_err());

    let one = call(&s, NEW, &[weights(&[0.0, 5.0, 0.0])]).unwrap();
    let tested = rng(&s, 42);
    let original = rng(&s, 42);
    assert!(call(&s, SAMPLE, &[one.clone(), rng(&foreign, 11)]).is_err());
    assert!(
        call(
            &s,
            FILL,
            &[one.clone(), tested.clone(), Value::int(262_144)]
        )
        .is_err()
    );
    let Value::Tensor(empty) =
        call(&s, FILL, &[one.clone(), tested.clone(), Value::int(0)]).unwrap()
    else {
        panic!("expected an empty tensor");
    };
    assert_eq!(empty.shape, [0]);
    assert_eq!(empty.elem, TypeExpr::Int);
    assert!(empty.data.is_empty());
    assert_eq!(
        call(&s, SAMPLE, &[one, tested.clone()]).unwrap(),
        Value::int(1)
    );
    assert_eq!(
        call(&s, NEXT, &[tested]).unwrap(),
        call(&s, NEXT, &[original]).unwrap()
    );
}
