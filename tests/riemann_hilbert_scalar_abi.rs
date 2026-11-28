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
use std::f64::consts::FRAC_PI_2;
use std::f64::consts::PI;
use std::f64::consts::TAU;
use std::sync::Arc;

fn cx(
    re: f64,
    im: f64,
) -> Value {
    Value::Complex(ComplexValue::Float { re, im })
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

fn complex(value: &Value) -> Complex64 {
    match value {
        | Value::Dyn(d) => complex(&d.value),
        | Value::Complex(z) => {
            let (re, im) = z.to_f64();
            Complex64::new(re, im)
        },
        | _ => panic!("expected a complex value"),
    }
}

fn field<'a>(
    value: &'a Value,
    key: &str,
) -> &'a Value {
    let Value::Record(record) = value else {
        panic!("expected status record")
    };
    record.get(key).unwrap()
}

fn as_scalar(value: &Value) -> f64 {
    if let Value::Dyn(d) = value {
        as_scalar(&d.value)
    } else {
        value.as_f64().unwrap()
    }
}

#[test]
fn scalar_cauchy_restores_measured_amplitude_and_actual_evaluations() {
    let session = Session::new();
    let time = DimExpr::base("Time");
    let unit = UnitValue {
        name: "ms".into(),
        display: "ms".into(),
        dimension: time.clone(),
        scale: ScaleValue {
            numer: 1.into(),
            denom: 1000.into(),
            pi: 0,
        },
        offset: None,
        affine: false,
        aliases: vec![],
        prefixable: false,
        category: "time".into(),
        since: String::new(),
    };
    let mut ctx = CallCtx::new();
    let callback = ctx
        .add_callback(Arc::new(NativeCallback::new(
            move |_: &[Value]| {
                Ok(Value::Quantity(QuantityValue {
                    value: Box::new(cx(3_000.0, 0.0)),
                    unit: unit.clone(),
                }))
            },
            Signature::new(
                vec![TypeExpr::Float],
                TypeExpr::Complex.with_dim(time.clone()),
            ),
        )))
        .unwrap();
    let evaluate = |point, budget| {
        session.call_ho(
            "num_riemann_hilbert_cauchy_checked",
            &[
                point,
                measured(Value::Float(2.0), &time),
                measured(Value::Float(1e-10), &time),
                Value::int(96),
                Value::int(budget),
            ],
            &[callback],
            &ctx,
        )
    };
    let answer = evaluate(cx(0.0, 1.0), 512).unwrap();
    assert_eq!(field(&answer, "converged"), &Value::Bool(true));
    assert!((complex(field(&answer, "value")) - 0.75).norm() < 1e-12);
    assert!(matches!(field(&answer,"value"), Value::Dyn(d) if d.dim == time));
    assert!(as_scalar(field(&answer, "error")) <= 1e-10);
    assert!(field(&answer, "evaluations").as_f64().unwrap() >= 15.0);
    assert!(matches!(
        evaluate(cx(0.0, 1.0), 1),
        Err(OpError::BudgetExceeded)
    ));
    assert!(evaluate(cx(0.0, 0.0), 512).is_err());
}

#[test]
fn plemelj_removes_the_central_singularity_and_distinguishes_panel_exhaustion() {
    let session = Session::new();
    let mut ctx = CallCtx::new();
    let signature = Signature::new(vec![TypeExpr::Float], TypeExpr::Complex);
    let linear = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |a: &[Value]| Ok(cx(a[0].as_f64().unwrap() + 1.0, 0.0)),
            signature.clone(),
        )))
        .unwrap();
    let oscillation = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |a: &[Value]| {
                let theta = 500.0 * a[0].as_f64().unwrap();
                Ok(cx(theta.cos(), theta.sin()))
            },
            signature,
        )))
        .unwrap();
    let exact = session
        .call_ho(
            "num_riemann_hilbert_plemelj_checked",
            &[
                Value::Float(0.0),
                Value::Float(1.0),
                Value::Float(1e-12),
                Value::int(2),
                Value::int(200),
            ],
            &[linear],
            &ctx,
        )
        .unwrap();
    assert_eq!(field(&exact, "converged"), &Value::Bool(true));
    let [upper, lower] = ["upper", "lower"].map(|key| complex(field(&exact, key)));
    assert!((upper - Complex64::new(0.5, -1.0 / PI)).norm() < 1e-12);
    assert!((lower - Complex64::new(-0.5, -1.0 / PI)).norm() < 1e-12);
    assert_eq!(field(&exact, "evaluations"), &Value::int(31));
    let capped = session
        .call_ho(
            "num_riemann_hilbert_cauchy_checked",
            &[
                cx(0.0, 1.0),
                Value::Float(1.0),
                Value::Float(1e-12),
                Value::int(1),
                Value::int(100),
            ],
            &[oscillation],
            &ctx,
        )
        .unwrap();
    assert_eq!(field(&capped, "converged"), &Value::Bool(false));
    assert!(complex(field(&capped, "value")).is_finite());
    assert!(as_scalar(field(&capped, "error")) > 1e-12);
}

#[test]
fn factorization_tracks_two_turns_and_checks_the_transformed_tolerance() {
    let session = Session::new();
    let mut ctx = CallCtx::new();
    let sig = Signature::new(vec![TypeExpr::Float], TypeExpr::Complex);
    let constant = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |_: &[Value]| Ok(cx(16.0, 0.0)),
            sig.clone(),
        )))
        .unwrap();
    let circle = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |v: &[Value]| {
                let arg = TAU * (v[0].as_f64().unwrap() + 1.0);
                Ok(cx(arg.cos(), arg.sin()))
            },
            sig.clone(),
        )))
        .unwrap();
    let zero = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |_: &[Value]| Ok(cx(0.0, 0.0)),
            sig,
        )))
        .unwrap();
    let call = |callback| {
        session.call_ho(
            "num_riemann_hilbert_factor_checked",
            &[
                cx(0.0, 1.0),
                Value::Float(1e-10),
                Value::int(32),
                Value::int(128),
                Value::int(2048),
            ],
            &[callback],
            &ctx,
        )
    };
    let actual = call(constant).unwrap();
    assert_eq!(field(&actual, "converged"), &Value::Bool(true));
    assert!((complex(field(&actual, "value")) - 2.0).norm() < 1e-12);
    assert!(as_scalar(field(&actual, "error")) <= 1e-10);
    let winding = call(circle).unwrap();
    let exact_log = Complex64::new(2.0 - FRAC_PI_2, FRAC_PI_2);
    assert!((complex(field(&winding, "logarithm")) - exact_log).norm() < 1e-10);
    assert!((complex(field(&winding, "value")) - exact_log.exp()).norm() < 1e-10);
    assert!(
        call(zero)
            .unwrap_err()
            .to_string()
            .contains("jump must be nonzero")
    );
}

#[test]
fn cauchy_fft_batch_preserves_value_units_and_exposes_sampling_aliases() {
    use rssn::kernel::value::TensorValue;
    let session = Session::new();
    let length = DimExpr::base("Length");
    let time = DimExpr::base("Time");
    let output = time.clone();
    let mut ctx = CallCtx::new();
    let cubic = ctx
        .add_callback(Arc::new(NativeCallback::new(
            move |args: &[Value]| {
                let z = complex(&args[0]).powu(3);
                Ok(measured(cx(z.re, z.im), &output))
            },
            Signature::new(
                vec![TypeExpr::Complex.with_dim(length.clone())],
                TypeExpr::Complex.with_dim(time.clone()),
            ),
        )))
        .unwrap();
    let aliquot = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |args: &[Value]| {
                let power = complex(&args[0]).powu(8);
                Ok(cx(power.re, power.im))
            },
            Signature::new(vec![TypeExpr::Complex], TypeExpr::Complex),
        )))
        .unwrap();
    let measured_args = [
        measured(cx(1.0, 0.0), &length),
        measured(Value::Float(2.0), &length),
        measured(Value::Float(1.0), &time),
        measured(Value::Float(1e-9), &time),
        Value::int(3),
        Value::int(16),
        Value::int(48),
    ];
    let found = session
        .call_ho(
            "num_complex_cauchy_fft_derivatives_checked",
            &measured_args,
            &[cubic],
            &ctx,
        )
        .unwrap();
    assert_eq!(field(&found, "converged"), &Value::Bool(true));
    assert_eq!(field(&found, "evaluations"), &Value::int(48));
    let Value::Tensor(TensorValue { data, elem, .. }) = field(&found, "derivatives") else {
        panic!("expected tensor")
    };
    assert!(
        rssn::kernel::overload::matches_type(&TypeExpr::Complex.with_dim(time), elem).accepted()
    );
    for (actual, wanted) in data.iter().zip([1.0, 6.0, 24.0, 48.0]) {
        assert!((complex(actual) - wanted).norm() < 1e-10);
    }
    let undersampled = session
        .call_ho(
            "num_complex_cauchy_fft_derivatives_checked",
            &[
                cx(0.0, 0.0),
                Value::Float(1.0),
                Value::Float(1.0),
                Value::Float(1e-10),
                Value::int(0),
                Value::int(8),
                Value::int(24),
            ],
            &[aliquot],
            &ctx,
        )
        .unwrap();
    assert_eq!(field(&undersampled, "converged"), &Value::Bool(false));
    assert!((as_scalar(field(&undersampled, "sample_error")) - 1.0).abs() < 1e-10);
}

#[test]
fn argument_principle_supports_only_the_zeroth_moment() {
    let session = Session::new();
    let mut ctx = CallCtx::new();
    let signature = Signature::new(vec![TypeExpr::Complex], TypeExpr::Complex);
    let f = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |args: &[Value]| Ok(args[0].clone()),
            signature.clone(),
        )))
        .unwrap();
    let df = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |_: &[Value]| Ok(cx(1.0, 0.0)),
            signature,
        )))
        .unwrap();
    let result = session
        .call_ho(
            "num_complex_root_circle_moments",
            &[
                cx(0.0, 0.0),
                Value::Float(1.0),
                Value::Float(1.0),
                Value::int(8),
                Value::int(0),
                Value::int(48),
            ],
            &[f, df],
            &ctx,
        )
        .unwrap();
    let Value::Tensor(tensor) = field(&result, "moments") else {
        panic!("expected a tensor of moments")
    };
    assert_eq!(tensor.shape, [1]);
    assert!((complex(&tensor.data[0]) - 1.0).norm() < 1e-12);
}
