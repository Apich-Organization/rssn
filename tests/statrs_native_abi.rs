//! Checked session-owned statistical distributions and extreme categorical weights.
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::TensorValue;
use rssn::kernel::value::Value;
use rssn::kernels::distributions_ext::StatrsDistribution;

fn vector(values: &[f64]) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![values.len()],
        elem: TypeExpr::Float,
        data: values.iter().copied().map(Value::Float).collect(),
    })
}

#[test]
fn weighted_categorical_is_scale_invariant_when_weights_would_overflow() {
    let baseline = StatrsDistribution::from_name("categorical", &[1.0, 1.0, 0.5]).unwrap();
    let extreme =
        StatrsDistribution::from_name("categorical", &[f64::MAX, f64::MAX, 0.5 * f64::MAX])
            .unwrap();
    for point in [0.0, 1.0, 2.0] {
        assert_eq!(baseline.pdf(point), extreme.pdf(point));
        assert_eq!(baseline.cdf(point), extreme.cdf(point));
        assert_eq!(baseline.sf(point), extreme.sf(point));
    }
    assert!(StatrsDistribution::from_name("normal", &[0.0, 1.0, 2.0]).is_err());
    assert!(
        StatrsDistribution::from_name("binomial", &[18_446_744_073_709_551_616.0, 0.5]).is_err()
    );
}

#[test]
fn native_statistics_expose_tails_moments_and_checked_ownership() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let new = |name: &str, params: &[f64]| {
        session.call_value(
            "num_statrs_distribution_new",
            &[Value::String(name.into()), vector(params)],
            &ctx,
        )
    };
    let normal = new("normal", &[0.0, 1.0]).unwrap();
    let call = |suffix: &str, args: &[Value]| {
        session
            .call_value(&format!("num_statrs_distribution_{suffix}"), args, &ctx)
            .unwrap()
    };
    let args = [normal.clone(), Value::Float(9.0)];
    assert_eq!(call("cdf", &args), Value::Float(1.0));
    let Some(tail) = call("sf", &args).as_f64() else {
        panic!("expected real tail")
    };
    assert!(tail > 1e-20 && tail < 2e-19);
    assert_eq!(
        call("mean", std::slice::from_ref(&normal)),
        Value::Float(0.0)
    );
    assert_eq!(
        call("variance", std::slice::from_ref(&normal)),
        Value::Float(1.0)
    );
    let entropy = call("entropy", std::slice::from_ref(&normal))
        .as_f64()
        .unwrap();
    assert!((entropy - 1.418_938_533_204_672_7).abs() < 1e-14);
    assert_eq!(
        call("is_discrete", std::slice::from_ref(&normal)),
        Value::Bool(false)
    );
    let uniform = new("discrete-uniform", &[-2.0, 2.0]).unwrap();
    assert_eq!(
        call("is_discrete", std::slice::from_ref(&uniform)),
        Value::Bool(true)
    );
    assert_eq!(
        call("is_signed_discrete", std::slice::from_ref(&uniform)),
        Value::Bool(true)
    );
    assert!(new("normal", &[0.0, -1.0]).is_err());
    assert!(
        Session::new()
            .call_value("num_statrs_distribution_mean", &[normal], &ctx)
            .is_err()
    );

    let beta = new("beta", &[0.5, 1.0]).unwrap();
    assert_eq!(
        call("pdf", &[beta, Value::Float(0.0)]),
        Value::Float(f64::INFINITY)
    );
}
