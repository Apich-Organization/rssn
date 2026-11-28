//! Round-trip of the portable fitted model through the public Session.

use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::Value;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::TensorValue;
use rssn::kernel::value::VariantValue;

fn tensor(
    shape: &[usize],
    values: &[f64],
) -> Value {
    Value::Tensor(TensorValue {
        shape: shape.to_vec(),
        data: values.iter().copied().map(Value::Float).collect(),
        elem: TypeExpr::Float,
    })
}

#[test]
fn fitted_model_retains_the_response_link_and_reports_bad_shapes() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let model = session
        .call_value(
            "num_glm_fit",
            &[
                tensor(&[4, 1], &[0.0, 1.0, 2.0, 3.0]),
                tensor(&[4], &[2.0, 4.0, 6.0, 8.0]),
                Value::String("gaussian".into()),
                Value::String("identity".into()),
                Value::Variant(VariantValue {
                    tag: "None".into(),
                    payload: vec![],
                }),
                Value::int(64),
                Value::Float(1e-12),
            ],
            &ctx,
        )
        .unwrap();
    let Value::Record(record) = &model else {
        panic!("expected a portable model record");
    };
    assert_eq!(record.get("link"), Some(&Value::String("identity".into())));
    let Value::Tensor(means) = session
        .call_value(
            "num_glm_predict",
            &[model.clone(), tensor(&[2, 1], &[4.0, 5.0])],
            &ctx,
        )
        .unwrap()
    else {
        panic!("expected means");
    };
    assert_eq!(means.shape, [2]);
    for (value, expected) in means.data.iter().zip([10.0, 12.0]) {
        assert!((value.as_f64().unwrap() - expected).abs() <= 1e-11);
    }
    assert!(
        session
            .call_value(
                "num_glm_predict",
                &[model, tensor(&[1, 2], &[2.0, 3.0])],
                &ctx
            )
            .is_err()
    );
}
