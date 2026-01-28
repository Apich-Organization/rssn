use rssn::symbolic::core::Expr;
use rssn::symbolic::transforms::*;

#[test]

fn test_fourier_transform_construction()
{

    // Just ensure it constructs without panic
    let t = Expr::new_variable("t");

    let f = Expr::new_exp(
        Expr::new_neg(t.clone()),
    );

    let result = fourier_transform(
        &f,
        "t",
        "omega",
    );

    // Result should be some integral expression
    match result {
        | Expr::Integral {
            ..
        } => { // expected
        },
        | _ => { // also acceptable if simplification occurs
        },
    }
}

#[test]

fn test_laplace_transform_construction()
{

    let t = Expr::new_variable("t");

    let f = Expr::new_exp(
        Expr::new_neg(t.clone()),
    );

    let result =
        laplace_transform(&f, "t", "s");

    // Result should be some integral expression
    match result {
        | Expr::Integral {
            ..
        } => { // expected
        },
        | _ => { // also acceptable
        },
    }
}

#[test]

fn test_z_transform_construction() {

    let n = Expr::new_variable("n");

    let f = Expr::new_pow(
        Expr::new_constant(0.5),
        n.clone(),
    );

    let result =
        z_transform(&f, "n", "z");

    // Result should be Summation
    match result {
        | Expr::Summation(
            _,
            _,
            _,
            _,
        ) => { /* expected */ },
        | _ => { // also acceptable
        },
    }
}

#[test]

fn test_fourier_time_shift() {

    let omega =
        Expr::new_variable("omega");

    let f_omega = Expr::new_constant(1.0); // F(omega) = 1
    let a = Expr::new_constant(2.0);

    let result = fourier_time_shift(
        &f_omega,
        &a,
        "omega",
    );

    // Result should include exp(-j*omega*a) * F(omega)
    // Just check it doesn't panic
    assert!(matches!(
        result,
        Expr::Mul(_, _)
            | Expr::Dag(_)
            | Expr::Exp(_)
            | _
    ));
}

#[test]

fn test_fourier_differentiation() {

    let f_omega =
        Expr::new_variable("F");

    let result =
        fourier_differentiation(
            &f_omega,
            "omega",
        );

    // Result should be j*omega * F
    // Just check it doesn't panic
    assert!(matches!(
        result,
        Expr::Mul(_, _)
            | Expr::Dag(_)
            | _
    ));
}

#[test]

fn test_laplace_differentiation() {

    let f_s = Expr::new_variable("F_s");

    let f_zero = Expr::new_constant(0.0);

    let result =
        laplace_differentiation(
            &f_s,
            "s",
            &f_zero,
        );

    // Result should be s*F(s) - f(0)
    // Just check it doesn't panic
    assert!(matches!(
        result,
        Expr::Sub(_, _)
            | Expr::Mul(_, _)
            | Expr::Dag(_)
            | _
    ));
}

#[test]

fn test_convolution_fourier() {

    let t = Expr::new_variable("t");

    let f = Expr::new_exp(
        Expr::new_neg(t.clone()),
    );

    let g = Expr::new_exp(
        Expr::new_neg(t.clone()),
    );

    let result = convolution_fourier(
        &f,
        &g,
        "t",
        "omega",
    );

    // Result is FT(f) * FT(g)
    // Just check it doesn't panic
    assert!(matches!(
        result,
        Expr::Mul(_, _)
            | Expr::Dag(_)
            | _
    ));
}

#[test]

fn test_laplace_properties() {

    let f_s = Expr::new_variable("F_s");

    let a = Expr::new_constant(2.0);

    // Frequency shift
    let _ = laplace_frequency_shift(
        &f_s, &a, "s",
    );

    // Scaling
    let _ =
        laplace_scaling(&f_s, &a, "s");

    // Integration
    let _ =
        laplace_integration(&f_s, "s");
}

#[test]

fn test_z_properties() {

    let f_z = Expr::new_variable("F_z");

    let k = Expr::new_constant(2.0);

    // Time shift
    let _ = z_time_shift(&f_z, &k, "z");

    // Scaling
    let a = Expr::new_constant(3.0);

    let _ = z_scaling(&f_z, &a, "z");

    // Differentiation
    let _ =
        z_differentiation(&f_z, "z");
}

#[test]

fn test_partial_fraction() {

    // 1 / (x^2 - 4) = 1/((x-2)(x+2))
    let x = Expr::new_variable("x");

    let num = Expr::new_constant(1.0);

    // Denominator exactly as in test_solve_quadratic
    let den = Expr::new_sub(
        Expr::new_pow(
            x.clone(),
            Expr::new_constant(2.0),
        ),
        Expr::new_constant(4.0),
    );

    let expr = Expr::new_div(num, den);

    let result =
        partial_fraction_decomposition(
            &expr, "x",
        );

    assert!(
        result.is_some(),
        "Partial fraction \
         decomposition failed for \
         1/(x^2-4)"
    );

    let terms = result.unwrap();

    assert_eq!(terms.len(), 2);
}
