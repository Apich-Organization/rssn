use assert_approx_eq::assert_approx_eq;
use rssn::numerical::differential_geometry::*;
use rssn::symbolic::coordinates::CoordinateSystem;

#[test]

fn test_metric_tensor_spherical() {

    // Spherical metric: diag(1, rho^2, rho^2 sin^2(phi))
    // x = rho sin(phi) cos(theta)
    // y = rho sin(phi) sin(theta)
    // z = rho cos(phi)
    // vars: rho, theta_sph, phi
    // g_11 = 1
    // g_22 = rho^2 sin^2(phi)
    // g_33 = rho^2
    // Wait, let's check CoordinateSystem::Spherical definition in coordinates.rs
    // sph_vars = [rho, theta_sph, phi]
    // rules = [rho sin(phi) cos(theta), rho sin(phi) sin(theta), rho cos(phi)]
    // h_rho = 1
    // h_theta = rho sin(phi)
    // h_phi = rho
    // g = diag(1, rho^2 sin^2(phi), rho^2)

    let point = vec![
        2.0,
        0.0,
        std::f64::consts::PI / 2.0,
    ]; // rho=2, theta=0, phi=pi/2
    let g = metric_tensor_at_point(
        CoordinateSystem::Spherical,
        &point,
    )
    .unwrap();

    assert_approx_eq!(
        g[0][0], 1.0, 1e-9
    );

    assert_approx_eq!(
        g[1][1], 4.0, 1e-9
    ); // 2^2 * sin^2(pi/2) = 4 * 1 = 4
    assert_approx_eq!(
        g[2][2], 4.0, 1e-9
    ); // 2^2 = 4
}

#[test]

fn test_christoffel_spherical() {

    let point = vec![
        1.0,
        0.0,
        std::f64::consts::PI / 2.0,
    ];

    let gamma = christoffel_symbols(
        CoordinateSystem::Spherical,
        &point,
    )
    .unwrap();

    // In spherical (rho, theta, phi):
    // Γ^rho_{theta theta} = -rho sin^2(phi) = -1 * 1 = -1
    // Γ^rho_{phi phi} = -rho = -1
    assert_approx_eq!(
        gamma[0][1][1],
        -1.0,
        1e-9
    );

    assert_approx_eq!(
        gamma[0][2][2],
        -1.0,
        1e-9
    );
}

#[test]

fn test_ricci_scalar_flat() {

    let point = vec![1.0, 2.0, 3.0];

    let r = ricci_scalar(
        CoordinateSystem::Cartesian,
        &point,
    )
    .unwrap();

    assert_approx_eq!(r, 0.0, 1e-9);
}
