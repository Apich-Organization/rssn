//! Comprehensive test suite for the Equations & Dynamical Systems domain integrator in RSSN.
//!
//! Validates:
//! - Symbolic Ordinary Differential Equation solving (`ode(equation, func, var)`).
//! - Numerical ODE integration via Runge-Kutta 4 (`solve_ode_system_rk4`).
//! - Symbolic equation root solving (`solve(equation, target)`).
//! - Numerical equation root solving via Newton-Raphson / Bisection (`solve_root`).
//! - Unified Identity Transformation Operator pipeline `compute(expr, config)`.

use std::sync::Arc;
use rssn::compute::{compute, d, eq, ode, solve, var, ComputeConfig};
use rssn::symbolic::core::Expr;

#[test]
fn test_symbolic_ode_exponential_eq() {
    // dy/dx = y -> y = C * e^x (or equivalent canonical form)
    let y = var("y");
    let dy_dx = d(y.clone(), "x");
    let ode_eq = eq(dy_dx, y.clone());

    let ode_op = ode(ode_eq, "y", "x");
    let config = ComputeConfig::default()
        .with_ode()
        .target_symbolic();

    let res = compute(&ode_op, &config);
    println!("Symbolic ODE dy/dx = y result: {:?}", res);

    match res {
        Expr::Eq(lhs, rhs) => {
            assert_eq!(*lhs, var("y"));
            let rhs_str = format!("{rhs}");
            assert!(
                rhs_str.contains("exp") || rhs_str.contains('e'),
                "RHS should contain exponential term, got: {rhs_str}"
            );
        }
        other => panic!("Expected Eq(y, ...), got: {:?}", other),
    }
}

#[test]
fn test_symbolic_ode_exponential_sub() {
    // dy/dx - y = 0 -> y = C * e^x (or equivalent canonical form)
    let y = var("y");
    let dy_dx = d(y.clone(), "x");
    let ode_eq = Expr::new_sub(dy_dx, y.clone());

    let ode_op = ode(ode_eq, "y", "x");
    let config = ComputeConfig::default()
        .with_ode()
        .target_symbolic();

    let res = compute(&ode_op, &config);
    println!("Symbolic ODE (dy/dx - y = 0) result: {:?}", res);

    match res {
        Expr::Eq(lhs, rhs) => {
            assert_eq!(*lhs, var("y"));
            let rhs_str = format!("{rhs}");
            assert!(
                rhs_str.contains("exp") || rhs_str.contains('e'),
                "RHS should contain exponential term, got: {rhs_str}"
            );
        }
        other => panic!("Expected Eq(y, ...), got: {:?}", other),
    }
}

#[test]
fn test_symbolic_ode_first_order_linear() {
    // dy/dx + y = 0 -> y = C * e^-x
    let y = var("y");
    let dy_dx = d(y.clone(), "x");
    let ode_eq = Expr::new_add(dy_dx, y.clone());

    let ode_op = ode(ode_eq, "y", "x");
    let config = ComputeConfig::default()
        .with_ode()
        .target_symbolic();

    let res = compute(&ode_op, &config);
    println!("Symbolic ODE (dy/dx + y = 0) result: {:?}", res);

    match res {
        Expr::Eq(lhs, rhs) => {
            assert_eq!(*lhs, var("y"));
            let rhs_str = format!("{rhs}");
            assert!(
                rhs_str.contains("exp"),
                "RHS should contain exponential term, got: {rhs_str}"
            );
        }
        other => panic!("Expected Eq(y, ...), got: {:?}", other),
    }
}

#[test]
fn test_numerical_ode_rk4_exponential_with_bindings() {
    // dy/dx = y with y(0) = 1.0, evaluated at x = 1.0 -> y(1.0) = e^1 approx 2.71828
    let y = var("y");
    let dy_dx = d(y.clone(), "x");
    let ode_eq = eq(dy_dx, y.clone());

    let ode_op = ode(ode_eq, "y", "x");
    let config = ComputeConfig::default()
        .with_ode()
        .bind("y", 1.0)
        .bind("x", 1.0)
        .target_numerical(1e-6);

    let res = compute(&ode_op, &config);
    println!("Numerical ODE dy/dx = y at x=1 result: {:?}", res);

    match res {
        Expr::Constant(val) => {
            let expected = std::f64::consts::E;
            assert!(
                (val - expected).abs() < 1e-2,
                "Expected approx e ({expected}), got {val}"
            );
        }
        other => panic!("Expected Expr::Constant, got: {:?}", other),
    }
}

#[test]
fn test_numerical_ode_rk4_with_builder_methods() {
    // dy/dx = y with y(0) = 2.0, integrated on [0, 1] -> y(1) = 2 * e approx 5.43656
    let y = var("y");
    let dy_dx = d(y.clone(), "x");
    let ode_eq = eq(dy_dx, y.clone());

    let ode_op = ode(ode_eq, "y", "x");
    let config = ComputeConfig::default()
        .with_ode()
        .with_initial_condition("y", 2.0)
        .with_ode_range(0.0, 1.0)
        .with_ode_steps(100)
        .target_numerical(1e-6);

    let res = compute(&ode_op, &config);
    println!("Numerical ODE dy/dx = y with y0=2.0 result: {:?}", res);

    match res {
        Expr::Constant(val) => {
            let expected = 2.0 * std::f64::consts::E;
            assert!(
                (val - expected).abs() < 1e-2,
                "Expected approx 2e ({expected}), got {val}"
            );
        }
        other => panic!("Expected Expr::Constant, got: {:?}", other),
    }
}

#[test]
fn test_numerical_root_solving_quadratic() {
    // solve(x^2 - 2 = 0, x) with numerical target -> 1.4142...
    let x = var("x");
    let x_sq = Expr::new_pow(x.clone(), Expr::new_constant(2.0));
    let solve_eq = eq(
        Expr::new_sub(x_sq, Expr::new_constant(2.0)),
        Expr::new_constant(0.0),
    );

    let solve_op = solve(solve_eq, "x");
    let config = ComputeConfig::default()
        .with_ode()
        .target_numerical(1e-7);

    let res = compute(&solve_op, &config);
    println!("Numerical root solve(x^2 - 2 = 0, x) result: {:?}", res);

    match res {
        Expr::Constant(val) => {
            let expected = std::f64::consts::SQRT_2;
            assert!(
                (val - expected).abs() < 1e-4,
                "Expected approx sqrt(2) ({expected}), got {val}"
            );
        }
        other => panic!("Expected Expr::Constant, got: {:?}", other),
    }
}

#[test]
fn test_numerical_root_solving_quadratic_negative_root() {
    // solve(x^2 - 2 = 0, x) with initial guess x = -1.0 -> -1.4142...
    let x = var("x");
    let x_sq = Expr::new_pow(x.clone(), Expr::new_constant(2.0));
    let solve_eq = Expr::new_sub(x_sq, Expr::new_constant(2.0)); // implicitly = 0

    let solve_op = solve(solve_eq, "x");
    let config = ComputeConfig::default()
        .with_ode()
        .bind("x", -1.0)
        .target_numerical(1e-7);

    let res = compute(&solve_op, &config);
    println!("Numerical root solve(x^2 - 2 = 0, x, start=-1) result: {:?}", res);

    match res {
        Expr::Constant(val) => {
            let expected = -std::f64::consts::SQRT_2;
            assert!(
                (val - expected).abs() < 1e-4,
                "Expected approx -sqrt(2) ({expected}), got {val}"
            );
        }
        other => panic!("Expected Expr::Constant, got: {:?}", other),
    }
}

#[test]
fn test_numerical_root_solving_transcendental() {
    // solve(cos(x) - x = 0, x) -> Dottie number approx 0.739085...
    let x = var("x");
    let cos_x = Expr::Cos(Arc::new(x.clone()));
    let solve_eq = Expr::new_sub(cos_x, x.clone());

    let solve_op = solve(solve_eq, "x");
    let config = ComputeConfig::default()
        .with_ode()
        .target_numerical(1e-7);

    let res = compute(&solve_op, &config);
    println!("Numerical root solve(cos(x) - x = 0, x) result: {:?}", res);

    match res {
        Expr::Constant(val) => {
            let expected = 0.7390851332;
            assert!(
                (val - expected).abs() < 1e-4,
                "Expected approx {expected}, got {val}"
            );
        }
        other => panic!("Expected Expr::Constant, got: {:?}", other),
    }
}

#[test]
fn test_symbolic_root_solving_quadratic() {
    // solve(x^2 - 2 = 0, x) in symbolic mode -> [sqrt(2), -sqrt(2)]
    let x = var("x");
    let x_sq = Expr::new_pow(x.clone(), Expr::new_constant(2.0));
    let solve_eq = eq(
        Expr::new_sub(x_sq, Expr::new_constant(2.0)),
        Expr::new_constant(0.0),
    );

    let solve_op = solve(solve_eq, "x");
    let config = ComputeConfig::default()
        .with_ode()
        .target_symbolic();

    let res = compute(&solve_op, &config);
    println!("Symbolic root solve(x^2 - 2 = 0, x) result: {:?}", res);

    let res_str = format!("{res}");
    assert!(
        res_str.contains("sqrt") || res_str.contains("1.414"),
        "Symbolic solution should contain square root, got: {res_str}"
    );
}

#[test]
fn test_symbolic_ode_with_initial_conditions() {
    // dy/dx = y with y(0) = 1 -> y = exp(x)
    let y = var("y");
    let dy_dx = d(y.clone(), "x");
    let ode_eq = eq(dy_dx, y.clone());

    let ode_op = ode(ode_eq, "y", "x");
    let config = ComputeConfig::default()
        .with_ode()
        .with_initial_condition("y", 1.0)
        .target_symbolic();

    let res = compute(&ode_op, &config);
    println!("Symbolic ODE with IC y(0)=1 result: {:?}", res);

    match &res {
        Expr::Eq(lhs, rhs) => {
            assert_eq!(**lhs, var("y"));
            let rhs_str = format!("{rhs}");
            assert!(
                rhs_str.contains("exp"),
                "RHS should contain exp(x), got: {rhs_str}"
            );
        }
        Expr::Exp(_) => {
            let res_str = format!("{res}");
            assert!(
                res_str.contains("exp") && res_str.contains('x'),
                "Result should be exp(x), got: {res_str}"
            );
        }
        other => panic!("Expected Eq(y, ...) or Exp(...), got: {:?}", other),
    }
}
