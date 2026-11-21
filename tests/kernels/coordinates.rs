//! Coordinate transformations and their Jacobians.

use rssn::kernels::coordinates::CoordinateSystem as C;
use rssn::kernels::coordinates::cartesian_to_spherical;
use rssn::kernels::coordinates::from_cartesian;
use rssn::kernels::coordinates::from_cartesian_pure;
use rssn::kernels::coordinates::jacobian_determinant;
use rssn::kernels::coordinates::numerical_jacobian;
use rssn::kernels::coordinates::spherical_to_cartesian;
use rssn::kernels::coordinates::to_cartesian;
use rssn::kernels::coordinates::to_cartesian_pure;
use rssn::kernels::coordinates::transform_point;
use rssn::kernels::coordinates::transform_point_pure;

fn det(m: &[Vec<f64>]) -> f64 {
    match m.len() {
        | 2 => m[0][0] * m[1][1] - m[0][1] * m[1][0],
        | 3 => {
            m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
        },
        | n => panic!("unsupported size {n}"),
    }
}

fn close(
    a: &[f64],
    b: &[f64],
    tol: f64,
) {
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(b) {
        assert!((x - y).abs() < tol, "{a:?} vs {b:?}");
    }
}

#[test]
fn known_points() {
    let p = transform_point(&[2.0, std::f64::consts::FRAC_PI_2], C::Polar, C::Cartesian)
        .unwrap_or_else(|e| panic!("{e}"));
    close(&p, &[0.0, 2.0], 1e-12);
    // theta from the z axis: theta = 0 is the north pole.
    let p = transform_point(&[3.0, 0.0, 1.0], C::Spherical, C::Cartesian)
        .unwrap_or_else(|e| panic!("{e}"));
    close(&p, &[0.0, 0.0, 3.0], 1e-12);
    let p = transform_point(&[0.0, 0.0, -2.0], C::Cartesian, C::Spherical)
        .unwrap_or_else(|e| panic!("{e}"));
    close(&p, &[2.0, std::f64::consts::PI, 0.0], 1e-12);
    let p = transform_point(&[0.0, 0.0, 0.0], C::Cartesian, C::Spherical)
        .unwrap_or_else(|e| panic!("{e}"));
    close(&p, &[0.0, 0.0, 0.0], 1e-15);
}

#[test]
fn round_trips() {
    let cases: [(C, &[f64]); 3] = [
        (C::Polar, &[1.7, 0.9]),
        (C::Cylindrical, &[1.7, -2.1, 0.4]),
        (C::Spherical, &[2.5, 1.1, -0.7]),
    ];
    for (system, p) in cases {
        let c = to_cartesian(p, system).unwrap_or_else(|e| panic!("{e}"));
        let back = from_cartesian(&c, system).unwrap_or_else(|e| panic!("{e}"));
        close(&back, p, 1e-12);
        let pure = transform_point_pure(
            &transform_point_pure(p, system, C::Cartesian),
            C::Cartesian,
            system,
        );
        close(&pure, p, 1e-12);
    }
    let s = transform_point(&[1.0, 2.0, 3.0], C::Cartesian, C::Spherical)
        .unwrap_or_else(|e| panic!("{e}"));
    let cyl = transform_point(&s, C::Spherical, C::Cylindrical).unwrap_or_else(|e| panic!("{e}"));
    let cart =
        transform_point(&cyl, C::Cylindrical, C::Cartesian).unwrap_or_else(|e| panic!("{e}"));
    close(&cart, &[1.0, 2.0, 3.0], 1e-12);
}

#[test]
fn spherical_polar_angle_resolves_directions_near_the_axis() {
    let point = [3e-16, 4e-16, 1.0];
    let actual = cartesian_to_spherical(point);
    let expected = [1.0, 5e-16, 4.0_f64.atan2(3.0)];
    assert_eq!(actual[0], expected[0]);
    assert!((actual[1] - expected[1]).abs() <= 1e-31);
    assert!((actual[2] - expected[2]).abs() <= f64::EPSILON);
    assert_eq!(from_cartesian_pure(&point, C::Spherical), actual);
    assert_eq!(
        transform_point(&point, C::Cartesian, C::Spherical).unwrap_or_else(|e| panic!("{e}")),
        actual
    );
    assert_eq!(
        transform_point_pure(&point, C::Cartesian, C::Spherical),
        actual
    );

    let round_trip = spherical_to_cartesian(actual);
    assert!((round_trip[0] - point[0]).abs() <= 1e-31);
    assert!((round_trip[1] - point[1]).abs() <= 1e-31);
    assert_eq!(round_trip[2], point[2]);
    assert_eq!(to_cartesian_pure(&actual, C::Spherical), round_trip);
    assert_eq!(
        transform_point(&actual, C::Spherical, C::Cartesian).unwrap_or_else(|e| panic!("{e}")),
        round_trip
    );
    assert_eq!(
        transform_point_pure(&actual, C::Spherical, C::Cartesian),
        round_trip
    );
}

#[test]
fn spherical_angles_survive_radius_overflow() {
    let max = f64::MAX;
    let [r, theta, phi] = cartesian_to_spherical([max, max, max]);
    assert_eq!(r, f64::INFINITY);
    assert!((theta - 0.955_316_618_124_509_3).abs() <= f64::EPSILON);
    assert_eq!(phi, std::f64::consts::FRAC_PI_4);

    // Only the full norm overflows.
    let [r, theta, phi] = cartesian_to_spherical([1.3e308, 0.0, 1.3e308]);
    assert_eq!(r, f64::INFINITY);
    assert_eq!(theta, std::f64::consts::FRAC_PI_4);
    assert_eq!(phi, 0.0);
}

#[test]
fn spherical_conversions_retain_subnormal_scale() {
    let min = f64::from_bits(1);
    let [r, theta, phi] = cartesian_to_spherical([min, min, min]);
    assert_eq!(r, 2.0 * min);
    assert!((theta - 0.955_316_618_124_509_3).abs() <= f64::EPSILON);
    assert_eq!(phi, std::f64::consts::FRAC_PI_4);

    // A subnormal transverse length can still give a normal polar angle.
    let [_, theta, _] = cartesian_to_spherical([min, min, 1e-300]);
    let expected = 6.987_143_370_513_132e-24;
    assert!((theta - expected).abs() <= 2e-15 * expected);

    // Each horizontal component is below half an ulp at this radius.
    assert_eq!(spherical_to_cartesian([min, 0.7, 0.7]), [0.0, 0.0, min]);
}

#[test]
fn spherical_signed_axis_and_nonfinite_conventions() {
    let [r, theta, phi] = cartesian_to_spherical([-0.0, -0.0, -0.0]);
    assert_eq!(r.to_bits(), 0);
    assert_eq!(theta.to_bits(), 0);
    assert_eq!(phi, -std::f64::consts::PI);
    let [r, theta, phi] = cartesian_to_spherical([1.0, 0.0, f64::INFINITY]);
    assert_eq!(r, f64::INFINITY);
    assert!(theta.is_nan());
    assert_eq!(phi, 0.0);
    let [r, theta, phi] = cartesian_to_spherical([f64::INFINITY, 0.0, 1.0]);
    assert_eq!(r, f64::INFINITY);
    assert_eq!(theta, std::f64::consts::FRAC_PI_2);
    assert_eq!(phi, 0.0);
    assert!(
        cartesian_to_spherical([f64::NAN, 0.0, 0.0])
            .iter()
            .all(|x| x.is_nan())
    );
    let direction = spherical_to_cartesian([2.0, 0.7, -1.0]);
    assert_eq!(
        spherical_to_cartesian([-2.0, 0.7, -1.0]),
        direction.map(|x| -x)
    );
}

#[test]
fn invalid_dimensions_are_errors() {
    assert!(transform_point(&[1.0, 2.0, 3.0], C::Polar, C::Cartesian).is_err());
    assert!(transform_point(&[1.0, 2.0], C::Cartesian, C::Spherical).is_err());
    assert!(transform_point(&[1.0, 2.0], C::Polar, C::Spherical).is_err());
    assert!(to_cartesian(&[1.0], C::Cartesian).is_err());
    assert!(from_cartesian(&[1.0, 2.0], C::Cylindrical).is_err());
}

#[test]
fn polar_jacobian_determinant_is_r() {
    let (r, t) = (2.3, 0.8);
    let j = numerical_jacobian(C::Polar, C::Cartesian, &[r, t]).unwrap_or_else(|e| panic!("{e}"));
    assert!((det(&j) - r).abs() < 1e-7);
    assert!((j[0][0] - t.cos()).abs() < 1e-8);
    assert!((j[1][1] - r * t.cos()).abs() < 1e-7);
}

#[test]
fn spherical_jacobian_determinant_is_r2_sin_theta() {
    let (r, t, p) = (1.9, 0.7, 2.2);
    let j = numerical_jacobian(C::Spherical, C::Cartesian, &[r, t, p])
        .unwrap_or_else(|e| panic!("{e}"));
    assert!((det(&j) - r * r * t.sin()).abs() < 1e-6, "{}", det(&j));
    // Cylindrical: determinant r.
    let j = numerical_jacobian(C::Cylindrical, C::Cartesian, &[r, p, 0.5])
        .unwrap_or_else(|e| panic!("{e}"));
    assert!((det(&j) - r).abs() < 1e-7);
}

#[test]
fn jacobian_of_inverse_map_is_inverse() {
    let cart = [0.6, -1.1, 0.9];
    let fwd =
        numerical_jacobian(C::Cartesian, C::Spherical, &cart).unwrap_or_else(|e| panic!("{e}"));
    let sph = transform_point(&cart, C::Cartesian, C::Spherical).unwrap_or_else(|e| panic!("{e}"));
    let back =
        numerical_jacobian(C::Spherical, C::Cartesian, &sph).unwrap_or_else(|e| panic!("{e}"));
    assert!((det(&fwd) * det(&back) - 1.0).abs() < 1e-6);
}

#[test]
fn jacobian_determinant_matches_the_known_volume_factors() {
    // polar -> Cartesian: |det J| = r.
    let j =
        jacobian_determinant(C::Polar, C::Cartesian, &[3.0, 0.7]).unwrap_or_else(|e| panic!("{e}"));
    assert!((j - 3.0).abs() < 1e-9, "{j}");
    // cylindrical -> Cartesian: |det J| = r.
    let j = jacobian_determinant(C::Cylindrical, C::Cartesian, &[4.0, -0.2, 1.5])
        .unwrap_or_else(|e| panic!("{e}"));
    assert!((j - 4.0).abs() < 1e-9, "{j}");
    // spherical -> Cartesian: |det J| = r^2 sin(theta).
    let (r, theta, phi) = (2.0, 1.1, -0.4);
    let j = jacobian_determinant(C::Spherical, C::Cartesian, &[r, theta, phi])
        .unwrap_or_else(|e| panic!("{e}"));
    let want = r * r * theta.sin();
    assert!((j - want).abs() < 1e-9 * want, "{j} vs {want}");
    // The inverse map carries the reciprocal determinant.
    let cart = transform_point(&[r, theta, phi], C::Spherical, C::Cartesian)
        .unwrap_or_else(|e| panic!("{e}"));
    let inv =
        jacobian_determinant(C::Cartesian, C::Spherical, &cart).unwrap_or_else(|e| panic!("{e}"));
    assert!((inv * j - 1.0).abs() < 1e-7, "{inv} * {j}");
    // Identity map: the central difference of `x +/- h` divided by `2h` is
    // one up to round-off, so the determinant is one to about 1e-10.
    let j = jacobian_determinant(C::Cartesian, C::Cartesian, &[1.0, 2.0, 3.0])
        .unwrap_or_else(|e| panic!("{e}"));
    assert!((j - 1.0).abs() < 1e-9, "{j}");
    let jac = numerical_jacobian(C::Spherical, C::Cartesian, &[r, theta, phi])
        .unwrap_or_else(|e| panic!("{e}"));
    let helper = jacobian_determinant(C::Spherical, C::Cartesian, &[r, theta, phi])
        .unwrap_or_else(|e| panic!("{e}"));
    assert!((det(&jac) - helper).abs() < 1e-12);
    // Mixing dimensions is an error.
    assert!(jacobian_determinant(C::Polar, C::Spherical, &[1.0, 0.5]).is_err());
}
