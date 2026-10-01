//! Coordinate transformations and their Jacobians.

use rssn::kernels::coordinates::{
    CoordinateSystem as C, from_cartesian, numerical_jacobian, to_cartesian, transform_point,
    transform_point_pure,
};

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

fn close(a: &[f64], b: &[f64], tol: f64) {
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
        let pure = transform_point_pure(&transform_point_pure(p, system, C::Cartesian), C::Cartesian, system);
        close(&pure, p, 1e-12);
    }
    let s = transform_point(&[1.0, 2.0, 3.0], C::Cartesian, C::Spherical)
        .unwrap_or_else(|e| panic!("{e}"));
    let cyl = transform_point(&s, C::Spherical, C::Cylindrical).unwrap_or_else(|e| panic!("{e}"));
    let cart = transform_point(&cyl, C::Cylindrical, C::Cartesian).unwrap_or_else(|e| panic!("{e}"));
    close(&cart, &[1.0, 2.0, 3.0], 1e-12);
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
    let fwd = numerical_jacobian(C::Cartesian, C::Spherical, &cart).unwrap_or_else(|e| panic!("{e}"));
    let sph = transform_point(&cart, C::Cartesian, C::Spherical).unwrap_or_else(|e| panic!("{e}"));
    let back = numerical_jacobian(C::Spherical, C::Cartesian, &sph).unwrap_or_else(|e| panic!("{e}"));
    assert!((det(&fwd) * det(&back) - 1.0).abs() < 1e-6);
}
