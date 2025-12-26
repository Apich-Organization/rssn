use assert_approx_eq::assert_approx_eq;
use rssn::prelude::numerical::*;
use rssn::prelude::*;
use rssn::symbolic::coordinates::CoordinateSystem;

#[test]

fn test_coord_transform_basic() {

    // Cartesian to Cylindrical: (1, 1, 1) -> (sqrt(2), 45°, 1)
    let p = vec![1.0, 1.0, 1.0];

    let res = numerical_transform_point_pure(
        &p,
        CoordinateSystem::Cartesian,
        CoordinateSystem::Cylindrical,
    )
    .unwrap();

    assert_approx_eq!(res[0], 2.0f64.sqrt());

    assert_approx_eq!(res[1], std::f64::consts::FRAC_PI_4);

    assert_approx_eq!(res[2], 1.0);
}

#[test]

fn test_coord_transform_spherical() {

    // Cartesian to Spherical: (1, 1, 1) -> (sqrt(3), theta=45°, phi=acos(1/sqrt(3)))
    let p = vec![1.0, 1.0, 1.0];

    let res = numerical_transform_point_pure(
        &p,
        CoordinateSystem::Cartesian,
        CoordinateSystem::Spherical,
    )
    .unwrap();

    assert_approx_eq!(res[0], 3.0f64.sqrt());

    assert_approx_eq!(res[1], std::f64::consts::FRAC_PI_4);

    assert_approx_eq!(res[2], (1.0 / 3.0f64.sqrt()).acos());
}

#[test]

fn test_coord_jacobian() {

    // Jacobian of Cartesian to Cylindrical at (1, 0, 0)
    // r = sqrt(x^2 + y^2), theta = atan2(y, x), z = z
    // dr/dx = x/r = 1, dr/dy = y/r = 0, dr/dz = 0
    // dtheta/dx = -y/r^2 = 0, dtheta/dy = x/r^2 = 1, dtheta/dz = 0
    // dz/dx = 0, dz/dy = 0, dz/dz = 1
    let p = vec![1.0, 0.0, 0.0];

    let jac = numerical_coordinates_jacobian(
        CoordinateSystem::Cartesian,
        CoordinateSystem::Cylindrical,
        &p,
    )
    .unwrap();

    assert_eq!(jac.rows(), 3);

    assert_eq!(jac.cols(), 3);

    // Expected Jacobian:
    // [[1, 0, 0],
    //  [0, 1, 0],
    //  [0, 0, 1]]
    assert_approx_eq!(jac.get(0, 0), 1.0);

    assert_approx_eq!(jac.get(0, 1), 0.0);

    assert_approx_eq!(jac.get(1, 0), 0.0);

    assert_approx_eq!(jac.get(1, 1), 1.0);
}

#[test]

fn test_coord_transform_symbolic_assisted() {

    // Test the non-pure version that uses symbolic rules.
    let p = vec![1.0, 0.0, 0.0];

    let res =
        numerical_transform_point(&p, CoordinateSystem::Cartesian, CoordinateSystem::Spherical)
            .unwrap();

    // (1, 0, 0) -> rho=1, theta=0, phi=acos(0)=pi/2
    assert_approx_eq!(res[0], 1.0);

    assert_approx_eq!(res[1], 0.0);

    assert_approx_eq!(res[2], std::f64::consts::FRAC_PI_2);
}
