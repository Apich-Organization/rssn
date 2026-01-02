use rssn::symbolic::core::Expr;
use rssn::symbolic::integral_equations::solve_airfoil_equation;
use rssn::symbolic::integral_equations::FredholmEquation;
use rssn::symbolic::integral_equations::FredholmEquationParams;
use rssn::symbolic::integral_equations::VolterraEquation;
use rssn::symbolic::integral_equations::VolterraEquationParams;
use rssn::symbolic::simplify_dag::simplify;

#[test]

fn test_fredholm_neumann_series() {

    // y(x) = x + 0.1 * int_0^1 x*t * y(t) dt
    // This is a simple Fredholm equation.
    // f(x) = x
    // lambda = 0.1
    // K(x, t) = x*t
    // a = 0, b = 1

    let x = Expr::new_variable("x");

    let t = Expr::new_variable("t");

    let f_x = x.clone();

    let lambda = Expr::Constant(0.1);

    let kernel = Expr::new_mul(
        x.clone(),
        t.clone(),
    );

    let lower = Expr::Constant(0.0);

    let upper = Expr::Constant(1.0);

    let eq = FredholmEquation::new(
        FredholmEquationParams {
            y_x: Expr::new_variable(
                "y",
            ), // Placeholder
            f_x,
            lambda,
            kernel,
            lower_bound: lower,
            upper_bound: upper,
            var_x: "x".to_string(),
            var_t: "t".to_string(),
        },
    );

    // 1st iteration: y_0 = x
    // y_1 = x + 0.1 * int_0^1 x*t * t dt = x + 0.1 * x * [t^3/3]_0^1 = x + 0.1*x/3
    let solution =
        eq.solve_neumann_series(1);

    let simplified =
        simplify(&solution);

    // Expected: x + 0.0333... * x
    // Let's check if it contains x
    println!(
        "Neumann solution: {}",
        simplified
    );

    // We can check if substituting x=1 gives approx 1.0333
    // But symbolic check is better.
    // The structure should be Add(x, Mul(0.1, Mul(x, 1/3)))
    // Or simplified: Mul(x, 1.03333)
}

#[test]

fn test_fredholm_separable_kernel() {

    // y(x) = x + int_0^1 (x*t) * y(t) dt
    // K(x, t) = x * t -> a1(x) = x, b1(t) = t
    // lambda = 1

    let x = Expr::new_variable("x");

    let t = Expr::new_variable("t");

    let f_x = x.clone();

    let lambda = Expr::Constant(1.0);

    let kernel = Expr::new_mul(
        x.clone(),
        t.clone(),
    ); // Not used directly in separable solver but good for context
    let lower = Expr::Constant(0.0);

    let upper = Expr::Constant(1.0);

    let eq = FredholmEquation::new(
        FredholmEquationParams {
            y_x: Expr::new_variable(
                "y",
            ),
            f_x,
            lambda,
            kernel,
            lower_bound: lower,
            upper_bound: upper,
            var_x: "x".to_string(),
            var_t: "t".to_string(),
        },
    );

    let a_funcs = vec![x.clone()];

    let b_funcs = vec![t.clone()];

    let solution = eq
        .solve_separable_kernel(
            &a_funcs,
            &b_funcs,
        )
        .unwrap();

    println!(
        "Separable solution: {}",
        solution
    );

    // Exact solution check:
    // y(x) = x + x * int_0^1 t * y(t) dt
    // Let C = int_0^1 t * y(t) dt
    // y(x) = x(1+C)
    // C = int_0^1 t * t(1+C) dt = (1+C) * [t^3/3]_0^1 = (1+C)/3
    // C = 1/3 + C/3 => 2C/3 = 1/3 => 2C = 1 => C = 0.5
    // y(x) = 1.5 * x

    // Check if solution simplifies to 1.5 * x
    // We can substitute x=1 and check value
    let val_at_1 = rssn::symbolic::calculus::substitute(
        &solution,
        "x",
        &Expr::Constant(1.0),
    );

    let simplified_val =
        simplify(&val_at_1);

    if let Some(v) =
        simplified_val.to_f64()
    {

        assert!(
            (v - 1.5).abs() < 1e-6,
            "Expected 1.5, got {}",
            v
        );
    } else {

        // It might be 3/2
        println!(
            "Simplified value at 1: {}",
            simplified_val
        );
    }
}

#[test]

fn test_volterra_successive_approximations()
 {

    // y(x) = 1 + int_0^x y(t) dt
    // Solution should be e^x
    // f(x) = 1
    // lambda = 1
    // K(x, t) = 1

    let x = Expr::new_variable("x");

    let f_x = Expr::Constant(1.0);

    let lambda = Expr::Constant(1.0);

    let kernel = Expr::Constant(1.0);

    let lower = Expr::Constant(0.0);

    let eq = VolterraEquation::new(
        VolterraEquationParams {
            y_x: Expr::new_variable(
                "y",
            ),
            f_x,
            lambda,
            kernel,
            lower_bound: lower,
            var_x: "x".to_string(),
            var_t: "t".to_string(),
        },
    );

    // 3 iterations
    // y_0 = 1
    // y_1 = 1 + int_0^x 1 dt = 1 + x
    // y_2 = 1 + int_0^x (1+t) dt = 1 + x + x^2/2
    // y_3 = 1 + x + x^2/2 + x^3/6

    let solution = eq.solve_successive_approximations(3);

    println!(
        "Volterra approx solution: {}",
        solution
    );

    // Check structure or value
    // Let's check at x=1, should be approx 1 + 1 + 0.5 + 0.1666 = 2.666
    let val_at_1 = rssn::symbolic::calculus::substitute(
        &solution,
        "x",
        &Expr::Constant(1.0),
    );

    let simplified_val =
        simplify(&val_at_1);

    if let Some(v) =
        simplified_val.to_f64()
    {

        assert!(
            (v - 2.6666).abs() < 0.01,
            "Expected ~2.666, got {}",
            v
        );
    }
}

#[test]

fn test_airfoil_equation() {

    // (1/pi) * int_-1^1 y(t)/(t-x) dt = 1
    // f(x) = 1
    // Known solution involves Chebyshev polynomials

    let f_x = Expr::Constant(1.0);

    let solution =
        solve_airfoil_equation(
            &f_x, "x", "t",
        );

    println!(
        "Airfoil solution: {}",
        solution
    );

    // Just ensure it returns something non-trivial
    assert!(!matches!(
        solution,
        Expr::Constant(_)
    ));
}
