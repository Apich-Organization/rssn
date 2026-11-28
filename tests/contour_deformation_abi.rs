//! Geometry and unit regressions for obstacle-aware contours.
use std::f64::consts::PI;
use std::f64::consts::TAU;
use std::sync::Arc;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering;

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
use rssn::kernel::value::TensorValue;
use rssn::kernels::contour_deform::BranchRay;
use rssn::kernels::contour_deform::ContourSegment;
use rssn::kernels::contour_deform::deform_path_around_obstacles;
use rssn::kernels::contour_deform::integrate_deformed;

#[test]
fn crossing_far_from_branch_point_takes_the_open_sheet_route() {
    for scale in [1e-130, 1.0, 1e130] {
        for angle in [0.0, 0.7, PI, -2.2] {
            let axis = Complex64::from_polar(1.0, angle);
            let origin = Complex64::new(0.4, -0.3) * scale;
            let ray = BranchRay { origin, angle };
            let a = origin + Complex64::new(3.0, -1.0) * axis * scale;
            let b = origin + Complex64::new(3.0, 1.0) * axis * scale;
            for (start, end, sign) in [(a, b, 1.0), (b, a, -1.0)] {
                let path = deform_path_around_obstacles(start, end, &[], &[ray], 0.05 * scale);
                assert!(
                    path.iter()
                        .any(|part| matches!(part, ContourSegment::Arc { .. }))
                );
                let integral = integrate_deformed(
                    &|z| Complex64::from(1.0) / (z - origin),
                    start,
                    end,
                    &[],
                    &[ray],
                    0.05 * scale,
                    1e-10,
                );
                let expected = sign * (2.0 * (1.0_f64 / 3.0).atan() - TAU);
                assert!(integral.re.abs() < 1e-9);
                assert!((integral.im - expected).abs() < 1e-9, "{angle}: {integral}");
            }
        }
    }
}

#[test]
fn pole_bypass_preserves_which_side_the_chord_passes_on() {
    for height in [-0.05_f64, 0.05] {
        let pole = Complex64::new(0.0, height);
        let a = Complex64::from(-1.0);
        let b = Complex64::from(1.0);
        let path = deform_path_around_obstacles(a, b, &[pole], &[], 0.1);
        let ContourSegment::Arc { center, radius, .. } = path[1] else {
            panic!("expected detour arc");
        };
        assert_eq!(center, pole);
        assert!(radius >= 0.1);
        let value = integrate_deformed(
            &|z| Complex64::from(1.0) / (z - pole),
            a,
            b,
            &[pole],
            &[],
            0.1,
            1e-10,
        );
        let expected = height.signum() * (PI - 2.0 * height.abs().atan());
        assert!((value.im - expected).abs() < 1e-10, "{value}");
        assert!(value.re.abs() < 1e-10);
    }
}

#[test]
fn unsafe_endpoints_and_mutually_blocked_rays_fail_closed() {
    let a = Complex64::new(1.0, -1.0);
    let b = Complex64::new(1.0, 1.0);
    let origin = Complex64::default();
    let positive = BranchRay { origin, angle: 0.0 };
    let negative = BranchRay { origin, angle: PI };
    assert!(deform_path_around_obstacles(a, b, &[], &[positive, negative], 0.1).is_empty());
    assert!(deform_path_around_obstacles(a, 1.0.into(), &[], &[positive], 0.1).is_empty());
    assert!(deform_path_around_obstacles(a, b, &[a], &[], 0.1).is_empty());
    assert!(deform_path_around_obstacles(a, b, &[], &[], f64::INFINITY).is_empty());
    let calls = AtomicUsize::new(0);
    let output = integrate_deformed(
        &|_| {
            calls.fetch_add(1, Ordering::Relaxed);
            Complex64::default()
        },
        a,
        b,
        &[],
        &[positive, negative],
        0.1,
        1e-8,
    );
    assert!(!output.is_finite());
    assert_eq!(calls.load(Ordering::Relaxed), 0);
}

#[test]
fn a_cut_tangent_to_a_pole_detour_is_rejected() {
    let start = Complex64::from(-2.0);
    let end = Complex64::from(2.0);
    let poles = [Complex64::default()];
    let tangent = BranchRay {
        origin: Complex64::new(-1.0, -0.5),
        angle: 0.0,
    };
    let separated = BranchRay {
        origin: Complex64::new(-1.0, -0.500001),
        ..tangent
    };
    assert!(deform_path_around_obstacles(start, end, &poles, &[tangent], 0.5).is_empty());
    assert!(!deform_path_around_obstacles(start, end, &poles, &[separated], 0.5).is_empty());
}

fn measured(
    value: Value,
    dim: DimExpr,
) -> Value {
    Value::Dyn(DynValue {
        value: Box::new(value),
        dim,
    })
}
fn complex(
    re: f64,
    im: f64,
) -> Value {
    Value::Complex(ComplexValue::Float { re, im })
}
fn tensor(
    elem: TypeExpr,
    data: Vec<Value>,
) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![data.len()],
        elem,
        data,
    })
}

#[test]
fn native_facades_preserve_units_and_do_not_sample_rejected_paths() {
    let length = DimExpr::base("Length");
    let time = DimExpr::base("Time");
    let integral = DimExpr::Mul(vec![length.clone(), time.clone()]);
    let position = TypeExpr::Complex.with_dim(length.clone());
    let inputs = [
        measured(complex(3.0, -1.0), length.clone()),
        measured(complex(3.0, 1.0), length.clone()),
        tensor(position.clone(), vec![]),
        tensor(position.clone(), vec![complex(0.0, 0.0)]),
        measured(Value::Float(0.05), length.clone()),
        tensor(TypeExpr::Float, vec![Value::Float(0.0)]),
    ];
    let session = Session::new();
    let path = session
        .call_value("num_contour_deform_checked", &inputs, &CallCtx::new())
        .unwrap();
    let Value::Record(path) = path else {
        panic!("expected path")
    };
    assert_eq!(path.get("closed"), Some(&Value::Bool(false)));
    let Value::Tensor(radii) = path.get("arc_radii").unwrap() else {
        panic!("expected measured radii")
    };
    assert_eq!(radii.elem, TypeExpr::Float.with_dim(length.clone()));
    assert_eq!(radii.shape, [1]);
    let mut args = inputs.to_vec();
    args.extend([
        measured(Value::Float(2.5), length.clone()),
        measured(Value::Float(7.0), time.clone()),
        measured(Value::Float(1e-9), integral),
        Value::int(64),
        Value::int(900),
    ]);
    let calls = Arc::new(AtomicUsize::new(0));
    let counter = Arc::clone(&calls);
    let codomain = time.clone();
    let mut ctx = CallCtx::new();
    let callback = ctx
        .add_callback(Arc::new(NativeCallback::new(
            move |_: &[Value]| {
                counter.fetch_add(1, Ordering::Relaxed);
                Ok(measured(complex(3.0, 0.0), codomain.clone()))
            },
            Signature::new(vec![position], TypeExpr::Complex.with_dim(time)),
        )))
        .unwrap();
    let output = session
        .call_ho(
            "num_contour_integrate_deformed_checked",
            &args,
            &[callback],
            &ctx,
        )
        .unwrap();
    let Value::Record(output) = output else {
        panic!("expected checked estimate")
    };
    let Value::Dyn(value) = output.get("value").unwrap() else {
        panic!("expected measured integral")
    };
    let Value::Complex(value) = value.value.as_ref() else {
        panic!("expected complex integral")
    };
    let (re, im) = value.to_f64();
    assert!(re.abs() < 1e-10 && (im - 6.0).abs() < 1e-10);
    assert_eq!(output.get("converged"), Some(&Value::Bool(true)));
    assert_eq!(
        output.get("evaluations"),
        Some(&Value::int(calls.load(Ordering::Relaxed)))
    );
    let before = calls.load(Ordering::Relaxed);
    let mut rejected = args.clone();
    rejected[1] = measured(complex(1.0, 0.0), length);
    assert!(
        session
            .call_ho(
                "num_contour_integrate_deformed_checked",
                &rejected,
                &[callback],
                &ctx
            )
            .is_err()
    );
    assert_eq!(calls.load(Ordering::Relaxed), before);
    args[10] = Value::int(1);
    assert_eq!(
        session.call_ho(
            "num_contour_integrate_deformed_checked",
            &args,
            &[callback],
            &ctx
        ),
        Err(OpError::BudgetExceeded)
    );
}
