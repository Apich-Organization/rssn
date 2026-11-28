//! Coordinate round trips, distance and validation across the hyperbolic models.

use num_complex::Complex64;
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::manifest::SignatureSource;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::DimRuleInfo;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::ComplexValue;
use rssn::kernel::value::DynValue;

const TO_DISK: &str = "num_uhp_to_disk";
const TO_UHP: &str = "num_disk_to_uhp";
const DISK_DISTANCE: &str = "num_poincare_disk_distance";
const UHP_DISTANCE: &str = "num_uhp_distance";

const fn complex(
    re: f64,
    im: f64,
) -> Value {
    Value::Complex(ComplexValue::Float { re, im })
}

fn decode(value: Value) -> Complex64 {
    let Value::Complex(value) = value else {
        panic!("expected a complex coordinate");
    };
    let (re, im) = value.to_f64();
    Complex64::new(re, im)
}

#[test]
fn cayley_charts_round_trip_and_preserve_distances() {
    let session = Session::new();
    let context = CallCtx::new();
    let p = complex(0.0, 1.0);
    let q = complex(0.0, 4.0);
    let a = session
        .call_value(TO_DISK, std::slice::from_ref(&p), &context)
        .unwrap();
    let b = session
        .call_value(TO_DISK, std::slice::from_ref(&q), &context)
        .unwrap();
    assert_eq!(decode(a.clone()), Complex64::new(0.0, 0.0));
    assert!((decode(b.clone()).re - 0.6).abs() <= 2.0 * f64::EPSILON);
    let round_trip = decode(
        session
            .call_value(TO_UHP, std::slice::from_ref(&b), &context)
            .unwrap(),
    );
    assert!((round_trip - Complex64::new(0.0, 4.0)).norm() <= 8.0 * f64::EPSILON);
    let half_plane = session
        .call_value(UHP_DISTANCE, &[p, q], &context)
        .unwrap()
        .as_f64()
        .unwrap();
    let disk = session
        .call_value(DISK_DISTANCE, &[a, b], &context)
        .unwrap()
        .as_f64()
        .unwrap();
    // The vertical geodesic has length log(y2/y1).
    assert!((half_plane - 4.0_f64.ln()).abs() <= 8.0 * f64::EPSILON);
    assert!((half_plane - disk).abs() <= 8.0 * f64::EPSILON);
}

#[test]
fn distances_keep_close_points_and_open_boundary_points_distinguishable() {
    let session = Session::new();
    let context = CallCtx::new();
    let distance = |name, a, b| {
        session
            .call_value(name, &[a, b], &context)
            .unwrap()
            .as_f64()
            .unwrap()
    };
    let close = distance(UHP_DISTANCE, complex(0.0, 1.0), complex(1e-12, 1.0));
    assert!((close / 1e-12 - 1.0).abs() <= 4.0 * f64::EPSILON);

    let radius = 1.0_f64.next_down();
    let near_boundary = distance(DISK_DISTANCE, complex(-radius, 0.0), complex(radius, 0.0));
    assert!((near_boundary - 4.0 * radius.atanh()).abs() < 4e-14);

    let low = f64::from_bits(1);
    let high = f64::MAX;
    let wide = distance(UHP_DISTANCE, complex(0.0, low), complex(0.0, high));
    assert!((wide - (high.ln() - low.ln())).abs() < 3e-12);
}

#[test]
fn nonaxial_disk_points_use_the_signed_radial_gap() {
    let session = Session::new();
    let context = CallCtx::new();
    let x = 0.956_304_755_963_035_4_f64;
    let y = 0.292_371_704_722_736_77;
    let interior = complex(x, y);
    let coincident = session
        .call_value(
            DISK_DISTANCE,
            &[interior.clone(), interior.clone()],
            &context,
        )
        .unwrap();
    assert_eq!(coincident.as_f64(), Some(0.0));
    let from_origin = session
        .call_value(
            DISK_DISTANCE,
            &[complex(0.0, 0.0), interior.clone()],
            &context,
        )
        .unwrap()
        .as_f64()
        .unwrap();
    assert!(from_origin.is_finite() && from_origin > 30.0);

    // The exact binary64 rational gap is positive even though the squared
    // modulus rounds to one. The inverse uses that gap in its numerator.
    let inverse = decode(session.call_value(TO_UHP, &[interior], &context).unwrap());
    let expected_im = 7.270_117_036_357_913e-16;
    assert!((inverse.im / expected_im - 1.0).abs() < 4e-15);
    let mapped = decode(
        session
            .call_value(TO_DISK, &[complex(inverse.re, inverse.im)], &context)
            .unwrap(),
    );
    assert!((mapped - Complex64::new(x, y)).norm() <= 2.0 * f64::EPSILON);

    // Stepping the real component upward crosses the circle in exact arithmetic.
    let exterior = complex(x.next_up(), y);
    assert!(matches!(
        session.call_value(TO_UHP, std::slice::from_ref(&exterior), &context),
        Err(OpError::Domain(_))
    ));
    assert!(matches!(
        session.call_value(DISK_DISTANCE, &[exterior, complex(0.0, 0.0)], &context),
        Err(OpError::Domain(_))
    ));
}

#[test]
fn chart_domains_reject_boundaries_and_measured_coordinates() {
    let session = Session::new();
    let context = CallCtx::new();
    for (op, args) in [
        (TO_DISK, vec![complex(0.0, 0.0)]),
        (TO_DISK, vec![complex(0.0, f64::MAX)]),
        (TO_UHP, vec![complex(1.0, 0.0)]),
        (DISK_DISTANCE, vec![complex(0.0, 0.0), complex(1.0, 0.0)]),
        (UHP_DISTANCE, vec![complex(0.0, 1.0), complex(0.0, -1.0)]),
    ] {
        assert!(
            matches!(
                session.call_value(op, &args, &context),
                Err(OpError::Domain(_))
            ),
            "{op}"
        );
    }
    let measured = Value::Dyn(DynValue {
        value: Box::new(complex(0.0, 1.0)),
        dim: DimExpr::base("Length"),
    });
    assert!(session.call_value(TO_DISK, &[measured], &context).is_err());
}

#[test]
fn hyperbolic_charts_have_declared_dimensionless_source_contracts() {
    let manifest = Session::new().manifest();
    for (op, path, arity, result) in [
        (TO_DISK, "uhp_to_disk", 1, TypeExpr::Complex),
        (TO_UHP, "disk_to_uhp", 1, TypeExpr::Complex),
        (DISK_DISTANCE, "poincare_disk_distance", 2, TypeExpr::Float),
        (UHP_DISTANCE, "uhp_distance", 2, TypeExpr::Float),
    ] {
        let op = manifest.op(op).unwrap();
        let path = format!("crate::kernels::complex_coords::{path}");
        assert_eq!(op.rust_path.as_deref(), Some(path.as_str()));
        assert_eq!(op.signature_source, SignatureSource::Declared);
        assert_eq!(op.dim_rule, Some(DimRuleInfo::Dimensionless));
        let [signature] = op.overloads.as_slice() else {
            panic!("expected one declared signature");
        };
        assert_eq!(signature.params, vec![TypeExpr::Complex; arity]);
        assert_eq!(signature.ret, result);
    }
}
