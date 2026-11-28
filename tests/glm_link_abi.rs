//! Exercise finite-domain GLM links through the portable RSSN value ABI.
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;

#[test]
fn typed_probability_link_round_trips_in_the_representable_tail() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let name = Value::String("cloglog".into());
    let result = session.call_value(
        "num_glm_inverse_link",
        &[name.clone(), Value::Float(-710.0)],
        &ctx,
    );
    let mean = result.expect("finite mean");
    assert!(mean.as_f64().expect("positive lower tail") > 0.0);
    let predictor = session
        .call_value("num_glm_link", &[name, mean], &ctx)
        .expect("finite predictor");
    assert!((predictor.as_f64().expect("predictor") + 710.0).abs() < 1e-12);
}

#[test]
fn variance_and_derivative_respect_named_support() {
    let session = Session::new();
    let ctx = CallCtx::new();
    assert_eq!(
        session.call_value(
            "num_glm_variance",
            &[Value::String("binomial".into()), Value::Float(0.25)],
            &ctx
        ),
        Ok(Value::Float(0.1875))
    );
    assert_eq!(
        session.call_value(
            "num_glm_link_derivative",
            &[Value::String("inverse".into()), Value::Float(2.0)],
            &ctx
        ),
        Ok(Value::Float(-0.25))
    );
    for (op, name, input) in [
        ("num_glm_link_derivative", "logit", -0.25),
        ("num_glm_link", "missing", 0.5),
        ("num_glm_inverse_link", "inverse", 0.0),
        ("num_glm_inverse_link", "logit", f64::INFINITY),
        ("num_glm_variance", "poisson", -1.0),
    ] {
        let result =
            session.call_value(op, &[Value::String(name.into()), Value::Float(input)], &ctx);
        assert!(
            matches!(result, Err(OpError::Domain(_))),
            "{op}: {result:?}"
        );
    }
}
