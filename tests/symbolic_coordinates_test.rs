use rssn::symbolic::coordinates::*;
use rssn::symbolic::core::Expr;

#[test]

fn test_transform_point_cylindrical_to_cartesian_symbolic() {

    // Point in cylindrical: (r, theta, z)
    let r = Expr::new_variable("r");

    let theta = Expr::new_variable("theta");

    let z = Expr::new_variable("z");

    let cyl_point = vec![
        r.clone(),
        theta.clone(),
        z.clone(),
    ];

    let cart_point = transform_point(
        &cyl_point,
        CoordinateSystem::Cylindrical,
        CoordinateSystem::Cartesian,
    )
    .unwrap();

    // x = r*cos(theta), y = r*sin(theta), z = z
    assert_eq!(cart_point.len(), 3);

    println!("x = {:?}", cart_point[0]);

    println!("y = {:?}", cart_point[1]);

    println!("z = {:?}", cart_point[2]);
}

#[test]

fn test_transform_point_spherical_to_cartesian_symbolic() {

    // Point in spherical: (rho, theta, phi)
    let rho = Expr::new_variable("rho");

    let theta = Expr::new_variable("theta");

    let phi = Expr::new_variable("phi");

    let sph_point = vec![
        rho.clone(),
        theta.clone(),
        phi.clone(),
    ];

    let cart_point = transform_point(
        &sph_point,
        CoordinateSystem::Spherical,
        CoordinateSystem::Cartesian,
    )
    .unwrap();

    // x = rho*sin(phi)*cos(theta)
    // y = rho*sin(phi)*sin(theta)
    // z = rho*cos(phi)
    assert_eq!(cart_point.len(), 3);

    println!("x = {:?}", cart_point[0]);

    println!("y = {:?}", cart_point[1]);

    println!("z = {:?}", cart_point[2]);
}

#[test]

fn test_transform_expression_symbolic() {

    // Transform expression f(r, theta) from cylindrical to Cartesian
    // f = r^2
    let r = Expr::new_variable("r");

    let f = Expr::new_pow(r.clone(), Expr::new_constant(2.0));

    let f_cart = transform_expression(
        &f,
        CoordinateSystem::Cylindrical,
        CoordinateSystem::Cartesian,
    )
    .unwrap();

    // In Cartesian: r^2 = x^2 + y^2
    println!("f in Cartesian: {:?}", f_cart);
}

#[test]

fn test_get_metric_tensor_cylindrical() {

    // Metric tensor for cylindrical coordinates
    // g = diag(1, r^2, 1)
    let g = get_metric_tensor(CoordinateSystem::Cylindrical).unwrap();

    if let Expr::Matrix(rows) = g {

        assert_eq!(rows.len(), 3);

        assert_eq!(rows[0].len(), 3);

        println!("g[0][0] = {:?}", rows[0][0]);

        println!("g[1][1] = {:?}", rows[1][1]); // Should be r^2
        println!("g[2][2] = {:?}", rows[2][2]);
    } else {

        panic!("Expected matrix");
    }
}

#[test]

fn test_get_metric_tensor_spherical() {

    // Metric tensor for spherical coordinates
    // g = diag(1, rho^2*sin^2(phi), rho^2)
    let g = get_metric_tensor(CoordinateSystem::Spherical).unwrap();

    if let Expr::Matrix(rows) = g {

        assert_eq!(rows.len(), 3);

        assert_eq!(rows[0].len(), 3);

        println!("g[0][0] = {:?}", rows[0][0]); // Should be 1
        println!("g[1][1] = {:?}", rows[1][1]); // Should be rho^2*sin^2(phi)
        println!("g[2][2] = {:?}", rows[2][2]); // Should be rho^2
    } else {

        panic!("Expected matrix");
    }
}

#[test]

fn test_coordinate_system_identity() {

    // Test that transforming from a system to itself returns the same point
    let x = Expr::new_variable("x");

    let y = Expr::new_variable("y");

    let z = Expr::new_variable("z");

    let point = vec![
        x.clone(),
        y.clone(),
        z.clone(),
    ];

    let result = transform_point(
        &point,
        CoordinateSystem::Cartesian,
        CoordinateSystem::Cartesian,
    )
    .unwrap();

    assert_eq!(result.len(), 3);

    assert_eq!(result[0], x);

    assert_eq!(result[1], y);

    assert_eq!(result[2], z);
}
