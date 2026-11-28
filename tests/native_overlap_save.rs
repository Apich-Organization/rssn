//! Prepared overlap-save convolution across the authenticated measured ABI.

use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::TensorValue;

const NEW: &str = "num_native_overlap_save_plan_new";
const APPLY: &str = "num_native_overlap_save_plan_apply";
const BATCH: &str = "num_native_overlap_save_plan_apply_batch";

fn input(
    values: &[f64],
    dimension: DimExpr,
) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![values.len()],
        data: values.iter().copied().map(Value::Float).collect(),
        elem: TypeExpr::Float.with_dim(dimension),
    })
}

fn indices(values: &[usize]) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![values.len()],
        data: values.iter().copied().map(Value::int).collect(),
        elem: TypeExpr::Int,
    })
}

fn call(
    session: &Session,
    name: &str,
    values: &[Value],
) -> Result<Value, OpError> {
    session.call_value(name, values, &CallCtx::new())
}

fn check_vector(
    value: &Value,
    expected: &[f64],
    dimension: DimExpr,
) {
    let Value::Tensor(tensor) = value else {
        panic!("missing measured tensor")
    };
    assert_eq!(tensor.elem, TypeExpr::Float.with_dim(dimension));
    assert_eq!(tensor.shape, [expected.len()]);
    assert_eq!(tensor.data.len(), expected.len());
    for (actual, &expected) in tensor.data.iter().zip(expected) {
        let Value::Float(actual) = actual else {
            panic!("missing numeric sample")
        };
        assert!((actual - expected).abs() < 1e-13);
    }
}

#[test]
fn prepared_filter_preserves_product_units_and_ragged_rows() {
    let session = Session::new();
    let length = DimExpr::base("Length");
    let time = DimExpr::base("Time");
    let output = DimExpr::Mul(vec![length.clone(), time.clone()]);
    let plan = call(&session, NEW, &[input(&[2.0, -1.0], time)]).unwrap();

    let single = call(
        &session,
        APPLY,
        &[plan.clone(), input(&[1.0, 2.0], length.clone())],
    )
    .unwrap();
    check_vector(&single, &[2.0, 3.0, -2.0], output.clone());
    let vacant = call(&session, APPLY, &[plan.clone(), input(&[], length.clone())]).unwrap();
    check_vector(&vacant, &[], output.clone());

    let samples = input(&[1.0, 2.0, 3.0], length);
    let batch = call(
        &session,
        BATCH,
        &[plan.clone(), samples.clone(), indices(&[0, 2, 2, 3])],
    )
    .unwrap();
    let Value::Record(batch) = batch else {
        panic!("missing ragged output record")
    };
    assert_eq!(batch.get("offsets"), Some(&indices(&[0, 3, 3, 5])));
    check_vector(
        batch.get("values").unwrap(),
        &[2.0, 3.0, -2.0, 6.0, -3.0],
        output,
    );
    assert!(call(&session, BATCH, &[plan, samples, indices(&[0, 3, 2, 3])]).is_err());
}

#[test]
fn small_prepared_spectra_share_capacity_and_release_invalidates_the_handle() {
    let session = Session::new();
    let taps = input(&[1.0], DimExpr::Dimless);
    let plans: Vec<_> = (0..24)
        .map(|_| call(&session, NEW, std::slice::from_ref(&taps)).unwrap())
        .collect();
    let last = plans.last().unwrap();
    let samples = input(&[4.0], DimExpr::Dimless);
    let foreign = Session::new();
    assert!(call(&foreign, APPLY, &[last.clone(), samples.clone()]).is_err());
    assert!(call(&session, APPLY, &[last.clone(), samples.clone()]).is_ok());
    assert!(call(&session, "num_native_release", std::slice::from_ref(last)).is_ok());
    assert!(call(&session, APPLY, &[last.clone(), samples]).is_err());
}
