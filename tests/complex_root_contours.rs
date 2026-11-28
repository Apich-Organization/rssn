use num_complex::Complex64;
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::NativeCallback;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::Signature;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::ComplexValue;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::QuantityValue;
use rssn::kernel::value::ScaleValue;
use rssn::kernel::value::UnitValue;
use std::sync::Arc;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering::Relaxed;

fn cx(z: Complex64) -> Value {
    Value::Complex(ComplexValue::Float { re: z.re, im: z.im })
}

fn measured(
    value: Value,
    dim: &DimExpr,
) -> Value {
    Value::Dyn(DynValue {
        value: Box::new(value),
        dim: dim.clone(),
    })
}

fn coordinate(value: &Value) -> Complex64 {
    match value {
        | Value::Complex(z) => {
            let (real, imaginary) = z.to_f64();
            Complex64::new(real, imaginary)
        },
        | Value::Dyn(v) => coordinate(&v.value),
        | _ => panic!("expected a coherent complex coordinate"),
    }
}

fn field<'a>(
    value: &'a Value,
    name: &str,
) -> &'a Value {
    let Value::Record(record) = value else {
        panic!("expected record")
    };
    record.get(name).unwrap()
}

fn scaled_output(
    z: Complex64,
    dim: &DimExpr,
) -> Value {
    Value::Quantity(QuantityValue {
        value: Box::new(cx(z * 1_000.0)),
        unit: UnitValue {
            name: "test milli-unit".into(),
            display: "milli-unit".into(),
            dimension: dim.clone(),
            scale: ScaleValue {
                numer: 1.into(),
                denom: 1000.into(),
                pi: 0,
            },
            offset: None,
            affine: false,
            aliases: vec![],
            prefixable: false,
            category: "test".into(),
            since: String::new(),
        },
    })
}

#[test]
fn measured_circle_normalizes_moments_and_accepts_scaled_callback_results() {
    let length = DimExpr::base("Length");
    let time = DimExpr::base("Time");
    let rate = DimExpr::Mul(vec![
        time.clone(),
        DimExpr::Inverse(Box::new(length.clone())),
    ]);
    let session = Session::new();
    let samples = Arc::new(AtomicUsize::new(0));
    let sampled = Arc::clone(&samples);
    let time_output = time.clone();
    let rate_output = rate.clone();
    let mut ctx = CallCtx::new();
    let f = ctx
        .add_callback(Arc::new(NativeCallback::new(
            move |args: &[Value]| {
                sampled.fetch_add(1, Relaxed);
                Ok(scaled_output(coordinate(&args[0]) - 1.5, &time_output))
            },
            Signature::new(
                vec![TypeExpr::Complex.with_dim(length.clone())],
                TypeExpr::Complex.with_dim(time.clone()),
            ),
        )))
        .unwrap();
    let df = ctx
        .add_callback(Arc::new(NativeCallback::new(
            move |_: &[Value]| Ok(measured(cx(Complex64::new(1.0, 0.0)), &rate_output)),
            Signature::new(
                vec![TypeExpr::Complex.with_dim(length.clone())],
                TypeExpr::Complex.with_dim(rate),
            ),
        )))
        .unwrap();
    let args = |budget| {
        vec![
            measured(cx(Complex64::new(1.0, 0.0)), &length),
            measured(Value::Float(2.0), &length),
            measured(Value::Float(2.0), &time),
            Value::int(48),
            Value::int(2),
            Value::int(budget),
        ]
    };
    assert!(matches!(
        session.call_ho(
            "num_complex_root_circle_moments",
            &args(287),
            &[f, df],
            &ctx
        ),
        Err(OpError::BudgetExceeded)
    ));
    assert_eq!(samples.load(Relaxed), 0);
    let result = session
        .call_ho(
            "num_complex_root_circle_moments",
            &args(288),
            &[f, df],
            &ctx,
        )
        .unwrap();
    let Value::Tensor(moments) = field(&result, "moments") else {
        panic!("expected moments")
    };
    assert_eq!(moments.shape, [3]);
    assert_eq!(moments.elem, TypeExpr::Complex);
    for (moment, expected) in moments.data.iter().zip([1.0, 0.25, 0.0625]) {
        assert!((coordinate(moment) - expected).norm() < 1e-12);
    }
    assert!(field(&result, "sample_delta").as_f64().unwrap() < 1e-12);
    assert_eq!(field(&result, "samples"), &Value::int(96));
    assert_eq!(field(&result, "evaluations"), &Value::int(288));
    assert_eq!(samples.load(Relaxed), 144);
}

#[test]
fn winding_rejects_fractional_alias_and_counts_poles_with_sign() {
    let session = Session::new();
    let mut ctx = CallCtx::new();
    let signature = Signature::new(vec![TypeExpr::Complex], TypeExpr::Complex);
    let f = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |args: &[Value]| Ok(cx(coordinate(&args[0]) - 0.9)),
            signature.clone(),
        )))
        .unwrap();
    let df = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |_: &[Value]| Ok(cx(Complex64::new(1.0, 0.0))),
            signature.clone(),
        )))
        .unwrap();
    let pole = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |args: &[Value]| Ok(cx(coordinate(&args[0]).inv())),
            signature.clone(),
        )))
        .unwrap();
    let pole_df = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |args: &[Value]| Ok(cx(-coordinate(&args[0]).inv().powu(2))),
            signature,
        )))
        .unwrap();
    let call = |n, left, right| {
        session.call_ho(
            "num_complex_root_circle_winding",
            &[
                cx(Complex64::new(0.0, 0.0)),
                Value::Float(1.0),
                Value::Float(1.0),
                Value::int(n),
                Value::int(1024),
            ],
            &[left, right],
            &ctx,
        )
    };
    assert!(
        call(8, f, df)
            .unwrap_err()
            .to_string()
            .contains("not near an integer")
    );
    assert_eq!(field(&call(64, f, df).unwrap(), "winding"), &Value::int(1));
    assert_eq!(
        field(&call(8, pole, pole_df).unwrap(), "winding"),
        &Value::int(-1)
    );
}

#[test]
fn rouche_keeps_sampled_semantics_and_host_failures() {
    let session = Session::new();
    let mut ctx = CallCtx::new();
    let signature = Signature::new(vec![TypeExpr::Complex], TypeExpr::Complex);
    let identity = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |args: &[Value]| Ok(args[0].clone()),
            signature.clone(),
        )))
        .unwrap();
    let half = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |_: &[Value]| Ok(cx(Complex64::new(0.5, 0.0))),
            signature.clone(),
        )))
        .unwrap();
    let nonfinite = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |_: &[Value]| Ok(cx(Complex64::new(f64::NAN, 0.0))),
            signature.clone(),
        )))
        .unwrap();
    let failed = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |_: &[Value]| Err(OpError::Host("sample failed".into())),
            signature,
        )))
        .unwrap();
    let call = |left, right| {
        session.call_ho(
            "num_complex_root_circle_rouche_samples",
            &[
                cx(Complex64::new(0.0, 0.0)),
                Value::Float(1.0),
                Value::Float(1.0),
                Value::int(32),
                Value::int(64),
            ],
            &[left, right],
            &ctx,
        )
    };
    assert_eq!(
        field(&call(identity, half).unwrap(), "satisfied_at_samples"),
        &Value::Bool(true)
    );
    assert_eq!(
        field(&call(half, identity).unwrap(), "satisfied_at_samples"),
        &Value::Bool(false)
    );
    assert!(call(nonfinite, identity).is_err());
    assert!(
        matches!(call(failed, identity), Err(OpError::Host(ref message)) if message=="sample failed")
    );
}

#[test]
fn raw_circle_algorithms_reject_nonfinite_samples_and_false_integer() {
    use rssn::kernels::complex_roots::argument_principle_moments;
    use rssn::kernels::complex_roots::argument_principle_winding;
    use rssn::kernels::complex_roots::rouche_check_circle;
    let one = Complex64::new(1.0, 0.0);
    let origin = Complex64::new(0.0, 0.0);
    let invalid = Complex64::new(f64::NAN, 0.0);
    assert!(argument_principle_winding(|z| z - 0.9, |_| one, origin, 1.0, 8).is_err());
    assert!(argument_principle_moments(|_| invalid, |_| one, origin, 1.0, 0, 8).is_err());
    assert!(rouche_check_circle(|_| invalid, |_| one, origin, 1.0, 8).is_err());
    // Both complex norms overflow while their ordering is representable.
    let large = Complex64::new(1.5e308, 1.5e308);
    let less_large = Complex64::new(1.4e308, 1.4e308);
    assert!(rouche_check_circle(|_| large, |_| less_large, origin, 1.0, 8).unwrap());
}
