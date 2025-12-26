use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;
use rssn::numerical::vector_calculus::*;
use rssn::symbolic::core::Expr;

#[test]

fn test_divergence_expr() {

    let x = Expr::new_variable("x");

    let y = Expr::new_variable("y");

    // F = [x^2, y^2] -> div F = 2x + 2y
    let f1 = Expr::new_pow(
        x.clone(),
        Expr::new_constant(2.0),
    );

    let f2 = Expr::new_pow(
        y.clone(),
        Expr::new_constant(2.0),
    );

    let div = divergence_expr(
        &[f1, f2],
        &["x", "y"],
        &[1.0, 2.0],
    )
    .unwrap();

    assert_approx_eq!(div, 6.0, 1e-5);
}

#[test]

fn test_curl_expr() {

    let x = Expr::new_variable("x");

    let y = Expr::new_variable("y");

    let _z = Expr::new_variable("z");

    // F = [y, -x, 0] -> curl F = [0, 0, -2]
    let f1 = y.clone();

    let f2 = Expr::new_mul(
        Expr::new_constant(-1.0),
        x.clone(),
    );

    let f3 = Expr::new_constant(0.0);

    let res = curl_expr(
        &[f1, f2, f3],
        &["x", "y", "z"],
        &[1.0, 1.0, 1.0],
    )
    .unwrap();

    assert_approx_eq!(
        res[0],
        0.0,
        1e-5
    );

    assert_approx_eq!(
        res[1],
        0.0,
        1e-5
    );

    assert_approx_eq!(
        res[2],
        -2.0,
        1e-5
    );
}

#[test]

fn test_laplacian() {

    let x = Expr::new_variable("x");

    let y = Expr::new_variable("y");

    // f = x^2 + 2y^2 -> laplacian = 2 + 4 = 6
    let f = Expr::new_add(
        Expr::new_pow(
            x,
            Expr::new_constant(2.0),
        ),
        Expr::new_mul(
            Expr::new_constant(2.0),
            Expr::new_pow(
                y,
                Expr::new_constant(2.0),
            ),
        ),
    );

    let lap = laplacian(
        &f,
        &["x", "y"],
        &[1.0, 1.0],
    )
    .unwrap();

    assert_approx_eq!(lap, 6.0, 1e-5);
}

#[test]

fn test_directional_derivative() {

    let x = Expr::new_variable("x");

    let y = Expr::new_variable("y");

    // f = x^2 + y^2, grad = [2x, 2y] = [2, 2] at (1,1)
    // direction = [1, 0] (unit vector already)
    // DD = [2, 2] . [1, 0] = 2
    let f = Expr::new_add(
        Expr::new_pow(
            x,
            Expr::new_constant(2.0),
        ),
        Expr::new_pow(
            y,
            Expr::new_constant(2.0),
        ),
    );

    let dd = directional_derivative(
        &f,
        &["x", "y"],
        &[1.0, 1.0],
        &[1.0, 0.0],
    )
    .unwrap();

    assert_approx_eq!(dd, 2.0, 1e-5);
}

proptest! {
    #[test]
    fn proptest_laplacian_linear(a in -10.0..10.0f64, b in -10.0..10.0f64) {
        let x = Expr::new_variable("x");
        let y = Expr::new_variable("y");
        // f = ax + by -> laplacian = 0
        let f = Expr::new_add(
            Expr::new_mul(Expr::new_constant(a), x),
            Expr::new_mul(Expr::new_constant(b), y)
        );
        let lap = laplacian(&f, &["x", "y"], &[1.0, 1.0]).unwrap();
        prop_assert!(lap.abs() < 1e-5);
    }
}
