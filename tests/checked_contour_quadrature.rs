//! Complex integration result status across paths, folds and adaptive residues.

use std::cell::Cell;
use std::f64::consts::TAU;

use num_complex::Complex64;
use rssn::kernels::contour_deform::ContourQuadratureError;
use rssn::kernels::contour_deform::ContourSegment;
use rssn::kernels::contour_deform::cauchy_principal_value_checked;
use rssn::kernels::contour_deform::integrate_contour_checked;
use rssn::kernels::contour_deform::integrate_segment_checked;
use rssn::kernels::contour_deform::residue_at_infinity_checked;
use rssn::kernels::quadrature::GkError;

fn line(
    a: f64,
    b: f64,
) -> ContourSegment {
    ContourSegment::Line {
        start: a.into(),
        end: b.into(),
    }
}

#[test]
fn checked_contour_reuses_samples_and_preserves_budget_status() {
    let seen = Cell::new(0);
    let result = integrate_contour_checked(
        |z| {
            seen.set(seen.get() + 1);
            Complex64::new(z.re, 2.0 * z.re)
        },
        &[line(0.0, 0.5), line(0.5, 1.0)],
        1e-10,
        24,
    )
    .unwrap();
    assert!(result.converged);
    assert_eq!(result.evaluations, seen.get());
    assert_eq!(result.evaluations, 30);
    assert!((result.value - Complex64::new(0.5, 1.0)).norm() < 1e-12);
    let result = integrate_segment_checked(
        |z| (4000.0 * Complex64::i() * z).exp(),
        line(0.0, 1.0),
        1e-10,
        24,
    )
    .unwrap();
    assert!(result.value.is_finite());
    assert!(!result.converged && result.error > 1e-4);
    assert!(result.evaluations <= 2 * 705);
    assert_eq!(
        integrate_contour_checked(
            |_| Complex64::default(),
            &[line(0.0, 0.5), line(0.5, 1.0)],
            1e-10,
            1
        ),
        Err(ContourQuadratureError::Quadrature(GkError::PanelBudget))
    );
}

#[test]
fn checked_segments_reject_nonfinite_callbacks_and_malformed_geometry() {
    let calls = Cell::new(0);
    let output = integrate_segment_checked(
        |_| {
            calls.set(calls.get() + 1);
            Complex64::new(1.0, f64::NAN)
        },
        line(0.0, 1.0),
        1e-10,
        24,
    );
    assert_eq!(
        output,
        Err(ContourQuadratureError::Quadrature(
            GkError::NonFiniteEstimate
        ))
    );
    assert_eq!(calls.get(), 1);
    let arc = ContourSegment::Arc {
        center: Complex64::default(),
        radius: -1.0,
        theta0: 0.0,
        theta1: 1.0,
    };
    let output = integrate_segment_checked(|_| Complex64::new(1.0, 0.0), arc, 1e-10, 24);
    assert_eq!(output, Err(ContourQuadratureError::InvalidGeometry));
}

#[test]
fn adaptive_circle_eliminates_the_fixed_grid_residue_alias() {
    let result = residue_at_infinity_checked(|z| z.powu(31), 1.0, 1e-9, 256).unwrap();
    assert!(result.converged && result.value.norm() < 1e-10);
    let result =
        residue_at_infinity_checked(|z| Complex64::new(1.0, 0.0) / z, 3.0, 1e-9, 24).unwrap();
    assert!(result.converged && (result.value + 1.0).norm() < 1e-10);
    let arc = ContourSegment::Arc {
        center: Complex64::default(),
        radius: 3.0,
        theta0: TAU,
        theta1: 0.0,
    };
    let result =
        integrate_segment_checked(|z| Complex64::new(1.0, 0.0) / z, arc, 1e-9, 24).unwrap();
    assert!(result.converged && (result.value + Complex64::new(0.0, TAU)).norm() < 1e-10);
}

#[test]
fn principal_value_pairs_sides_counts_calls_and_supports_reversed_limits() {
    for (a, b, sign) in [(-1.0, 2.0, 1.0), (2.0, -1.0, -1.0)] {
        let calls = Cell::new(0);
        let result = cauchy_principal_value_checked(
            |x| {
                calls.set(calls.get() + 1);
                Complex64::new(1.0 / x, 0.0)
            },
            a,
            b,
            &[0.0],
            1e-10,
            64,
        )
        .unwrap();
        assert!(result.converged && (result.value.re - sign * 2.0_f64.ln()).abs() < 1e-10);
        assert_eq!(result.evaluations, calls.get());
        assert!(result.evaluations >= 45);
    }
    assert_eq!(
        cauchy_principal_value_checked(|_| Complex64::default(), -1.0, 1.0, &[0.0, 0.0], 1e-10, 24),
        Err(ContourQuadratureError::InvalidGeometry)
    );
}

#[test]
fn declared_contour_callbacks_round_trip_measured_path_and_report_failures() {
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

    let length = DimExpr::base("Length");
    let time = DimExpr::base("Time");
    let integral = DimExpr::Mul(vec![length.clone(), time.clone()]);
    let measured = |value, dim| {
        Value::Dyn(DynValue {
            value: Box::new(value),
            dim,
        })
    };
    let scalar = |v, d| measured(Value::Float(v), d);
    let session = Session::new();
    let path = session
        .call_value(
            "num_contour_bromwich",
            &[scalar(0.0, length.clone()), scalar(2.0, length.clone())],
            &CallCtx::new(),
        )
        .unwrap();
    let actual = Arc::new(AtomicUsize::new(0));
    let count = Arc::clone(&actual);
    let mut ctx = CallCtx::new();
    let milliseconds = UnitValue {
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
    let cb = ctx
        .add_callback(Arc::new(NativeCallback::new(
            move |_: &[Value]| {
                count.fetch_add(1, Relaxed);
                Ok(Value::Quantity(QuantityValue {
                    value: Box::new(Value::Complex(ComplexValue::Float { re: 3000.0, im: 0.0 })),
                    unit: milliseconds.clone(),
                }))
            },
            Signature::new(
                vec![TypeExpr::Complex.with_dim(length.clone())],
                TypeExpr::Complex.with_dim(time.clone()),
            ),
        )))
        .unwrap();
    let args = |budget| {
        [
            path.clone(),
            scalar(2.0, length.clone()),
            scalar(3.0, time.clone()),
            scalar(1e-9, integral.clone()),
            Value::int(24),
            Value::int(budget),
        ]
    };
    let output = session
        .call_ho(
            "num_contour_integrate_segment_checked",
            &args(900),
            &[cb],
            &ctx,
        )
        .unwrap();
    let Value::Record(output) = output else {
        panic!("expected integration status")
    };
    assert_eq!(output.get("converged"), Some(&Value::Bool(true)));
    assert_eq!(
        output.get("evaluations"),
        Some(&Value::int(actual.load(Relaxed)))
    );
    let Value::Dyn(value) = output.get("value").unwrap() else {
        panic!("expected measured integral")
    };
    let Value::Complex(estimate) = value.value.as_ref() else {
        panic!("expected complex integral")
    };
    let (re, im) = estimate.to_f64();
    assert!(re.abs() < 1e-12 && (im - 12.0).abs() < 1e-12);
    let err = session
        .call_ho("num_contour_integrate_checked", &args(1), &[cb], &ctx)
        .unwrap_err();
    assert_eq!(err, OpError::BudgetExceeded);
    let mut corrupt = args(900);
    let Value::Record(record) = &mut corrupt[0] else {
        panic!("expected path record")
    };
    let starts = &mut record
        .fields
        .iter_mut()
        .find(|(name, _)| name == "starts")
        .unwrap()
        .1;
    let Value::Tensor(starts) = starts else {
        panic!("expected coordinates")
    };
    starts.data[0] = measured(
        Value::Complex(ComplexValue::Float { re: 0.0, im: 1.0 }),
        time.clone(),
    );
    let before = actual.load(Relaxed);
    assert!(
        session
            .call_ho("num_contour_integrate_checked", &corrupt, &[cb], &ctx)
            .is_err()
    );
    assert_eq!(actual.load(Relaxed), before);

    let mut fail_ctx = CallCtx::new();
    let failure = fail_ctx
        .add_callback(Arc::new(NativeCallback::new(
            |_: &[Value]| Err(OpError::Host("user contour failure".into())),
            Signature::new(
                vec![TypeExpr::Complex.with_dim(length.clone())],
                TypeExpr::Complex.with_dim(time.clone()),
            ),
        )))
        .unwrap();
    let err = session
        .call_ho(
            "num_contour_integrate_checked",
            &args(900),
            &[failure],
            &fail_ctx,
        )
        .unwrap_err();
    assert!(matches!(err,OpError::Host(message) if message=="user contour failure"));
}

#[test]
fn cancellation_interrupts_the_panel_without_reporting_convergence() {
    use std::sync::Arc;
    use std::sync::atomic::AtomicBool;
    use std::sync::atomic::Ordering::Relaxed;

    let token = Arc::new(AtomicBool::new(false));
    let _scope = rssn::cancellation::install(Arc::clone(&token));
    let calls = Cell::new(0);
    let outcome = integrate_segment_checked(
        |_| {
            calls.set(calls.get() + 1);
            token.store(true, Relaxed);
            Complex64::new(1.0, 0.0)
        },
        line(0.0, 1.0),
        1e-10,
        24,
    );
    assert_eq!(
        outcome,
        Err(ContourQuadratureError::Quadrature(GkError::Cancelled))
    );
    assert_eq!(calls.get(), 1);
}
