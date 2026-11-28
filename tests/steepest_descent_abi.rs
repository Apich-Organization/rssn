use num_complex::Complex64 as C;
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
use rssn::kernels::steepest_descent as s;
use std::f64::consts::PI;
use std::sync::Arc;

fn complex(z: C) -> Value {
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

fn coordinate(value: &Value) -> C {
    match value {
        | Value::Complex(z) => {
            let (re, im) = z.to_f64();
            C::new(re, im)
        },
        | Value::Dyn(z) => coordinate(&z.value),
        | _ => panic!("expected complex coordinate"),
    }
}

#[test]
fn finds_only_verified_slopes_and_retains_tiny_nonzero_curvature() {
    assert!(
        s::find_saddle_point(
            |_| C::new(0.0, 0.0),
            |_| C::new(1.0, 0.0),
            |_| C::new(1e300, 0.0),
            C::new(1.0, 0.0),
            1e-12
        )
        .is_none()
    );
    let simple = s::find_saddle_point(
        |z| -0.5e-220 * z * z,
        |z| -1e-220 * z,
        |_| C::new(-1e-220, 0.0),
        C::new(0.0, 0.0),
        1e-240,
    )
    .unwrap();
    assert_eq!(simple.degeneracy, 2);
    let val = s::saddle_point_leading(|_| C::new(1.0, 0.0), &simple, 1e220);
    assert!((val.re - (2.0 * PI).sqrt()).abs() < 1e-13);
}

#[test]
fn vanishing_amplitude_keeps_the_mixed_cubic_correction() {
    let phi = |z: C| -0.5 * z.powi(2) + z.powi(3) / 6.0;
    let saddle = s::find_saddle_point(
        phi,
        |z| -z + 0.5 * z * z,
        |z| z - 1.0,
        C::new(0.0, 0.0),
        1e-12,
    )
    .unwrap();
    let lambda: f64 = 80.0;
    let expected = (2.0 * PI / lambda).sqrt() / (2.0 * lambda);
    let actual = s::saddle_point_asymptotic_order2(phi, |z| z, &saddle, lambda);
    assert!(
        (actual - expected).norm() < 1e-10 * expected,
        "{actual} vs {expected}"
    );
}

#[test]
fn tracing_respects_gaussian_orientation_and_validates_phase() {
    let phi = |z: C| 0.5 * z * z;
    let g = |_| C::new(1.0, 0.0);
    let saddle =
        s::find_saddle_point(phi, |z| z, |_| C::new(1.0, 0.0), C::new(0.1, 0.1), 1e-14).unwrap();
    let want = C::new(0.0, (2.0 * PI / 8.0).sqrt());
    let result = s::steepest_descent_integrate(phi, g, &saddle, 8.0, 4.0, 1e-12);
    assert!((result - want).norm() < 1e-12);
    let reversed = s::SaddlePoint {
        steepest_angle: saddle.steepest_angle + PI,
        ..saddle
    };
    assert!(
        (s::steepest_descent_integrate(phi, g, &reversed, 8.0, 4.0, 1e-12) + want).norm() < 1e-12
    );
    assert!(
        !s::steepest_descent_integrate(|z| phi(z) + 1.0, g, &saddle, 8.0, 4.0, 1e-12).is_finite()
    );
}

#[test]
fn curved_trace_recovers_the_transformed_gaussian() {
    let c = C::new(0.0, 0.12);
    let phase = |z: C| -0.5 * (z + c * z * z).powi(2);
    let amplitude = |z: C| 1.0 + 2.0 * c * z;
    let saddle = s::find_saddle_point(
        phase,
        |_| C::new(0.0, 0.0),
        |_| C::new(-1.0, 0.0),
        C::new(0.0, 0.0),
        1e-13,
    )
    .unwrap();
    let traced = s::steepest_descent_integrate(phase, amplitude, &saddle, 8.0, 4.0, 1e-12);
    assert!((traced - (2.0 * PI / 8.0).sqrt()).norm() < 1e-10);
}

#[test]
fn measured_saddle_record_round_trips_and_preserves_callback_errors() {
    let session = Session::new();
    let axis = DimExpr::base("Length");
    let amplitude = DimExpr::base("Time");
    let coord = TypeExpr::Complex.with_dim(axis.clone());
    let slope = DimExpr::Inverse(Box::new(axis.clone()));
    let curvature = DimExpr::Mul(vec![slope.clone(), slope.clone()]);
    let mut ctx = CallCtx::new();
    let phase = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |args: &[Value]| {
                let x = coordinate(&args[0]) - 1.0;
                Ok(complex(-0.5 * x * x))
            },
            Signature::new(vec![coord.clone()], TypeExpr::Complex),
        )))
        .unwrap();
    let d1 = slope.clone();
    let df = ctx
        .add_callback(Arc::new(NativeCallback::new(
            move |args: &[Value]| Ok(measured(complex(1.0 - coordinate(&args[0])), &d1)),
            Signature::new(vec![coord.clone()], TypeExpr::Complex.with_dim(slope)),
        )))
        .unwrap();
    let d2 = curvature.clone();
    let ddf = ctx
        .add_callback(Arc::new(NativeCallback::new(
            move |_: &[Value]| Ok(measured(complex(C::new(-1.0, 0.0)), &d2)),
            Signature::new(vec![coord.clone()], TypeExpr::Complex.with_dim(curvature)),
        )))
        .unwrap();
    let g_output = amplitude.clone();
    let g_signature = Signature::new(vec![coord], TypeExpr::Complex.with_dim(amplitude));
    let amplitude = ctx
        .add_callback(Arc::new(NativeCallback::new(
            move |_: &[Value]| Ok(measured(complex(C::new(3.0, 0.0)), &g_output)),
            g_signature.clone(),
        )))
        .unwrap();
    let failure = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |_: &[Value]| {
                Err(OpError::HostDiagnostic {
                    code: "E-test-saddle",
                    message: "failed amplitude".into(),
                })
            },
            g_signature,
        )))
        .unwrap();
    let find = |budget| {
        session.call_ho(
            "num_saddle_find",
            &[
                measured(complex(C::new(1.25, 0.1)), &axis),
                measured(Value::Float(2.0), &axis),
                Value::Float(1e-12),
                Value::int(budget),
            ],
            &[phase, df, ddf],
            &ctx,
        )
    };
    assert!(matches!(find(3), Err(OpError::BudgetExceeded)));
    let Value::Variant(found) = find(16).unwrap() else {
        panic!("missing saddle option")
    };
    assert_eq!(found.tag, "Some");
    let args = [
        found.payload[0].clone(),
        measured(Value::Float(2.0), &axis),
        Value::Float(8.0),
        Value::int(8),
    ];
    let lead = session
        .call_ho("num_saddle_leading", &args, &[amplitude], &ctx)
        .unwrap();
    assert!((coordinate(&lead).re - 3.0 * (2.0 * PI / 8.0).sqrt()).abs() < 1e-12);
    assert!(
        matches!(session.call_ho("num_saddle_leading", &args, &[failure], &ctx),
        Err(OpError::HostDiagnostic { code, .. }) if code == "E-test-saddle")
    );
}
