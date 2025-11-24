// File: tests/symbolic/ode.rs

use rssn::symbolic::core::Expr;
use rssn::symbolic::ode::{
    solve_first_order_linear_ode, solve_ode, solve_ode_system, solve_riccati_ode,
    solve_separable_ode,
};
use rssn::symbolic::simplify_dag::simplify;
use std::sync::Arc;

fn var(name: &str) -> Expr {
    Expr::Variable(name.to_string())
}

fn c(val: f64) -> Expr {
    Expr::Constant(val)
}

#[test]
fn test_separable_ode() {
    // y' = y
    // dy/dx = y => dy/y = dx => ln|y| = x + C => y = C*e^x
    let y = var("y");
    let x = var("x");
    let y_prime = Expr::Derivative(Arc::new(y.clone()), "x".to_string());
    let eq = Expr::new_sub(y_prime, y.clone()); // y' - y = 0

    let sol = solve_separable_ode(&eq, "y", "x");
    assert!(sol.is_some());
    let sol_expr = sol.unwrap();
    println!("Separable Solution: {}", sol_expr);
    // Expected: ln(y) = x + C or similar implicit form
    // The solver returns implicit form: int(1/y) = int(1) + C
    // ln(y) = x + C
}

#[test]
fn test_first_order_linear() {
    // y' + y = x
    // Integrating factor e^x.
    // d(y e^x)/dx = x e^x
    // y e^x = x e^x - e^x + C
    // y = x - 1 + C e^-x
    let y = var("y");
    let x = var("x");
    let y_prime = Expr::Derivative(Arc::new(y.clone()), "x".to_string());
    let eq = Expr::new_sub(Expr::new_add(y_prime, y.clone()), x.clone());

    let sol = solve_first_order_linear_ode(&eq, "y", "x");
    assert!(sol.is_some());
    let sol_expr = sol.unwrap();
    println!("Linear Solution: {}", sol_expr);
}

#[test]
fn test_riccati_ode() {
    // y' = 1 + y^2
    // y1 = tan(x)
    let y = var("y");
    let x = var("x");
    let y_prime = Expr::Derivative(Arc::new(y.clone()), "x".to_string());

    // y' - (1 + y^2) = 0
    let rhs = Expr::new_add(c(1.0), Expr::new_pow(y.clone(), c(2.0)));
    let eq = Expr::new_sub(y_prime, rhs);

    let y1 = Expr::Tan(Arc::new(x.clone()));

    let sol = solve_riccati_ode(&eq, "y", "x", &y1);
    assert!(sol.is_some());
    let sol_expr = sol.unwrap();
    println!("Riccati Solution: {}", sol_expr);
}

#[test]
fn test_ode_system() {
    // Triangular system:
    // y' = x  => y = x^2/2 + C1
    // z' = y  => z = x^3/6 + C1*x + C2
    let y = var("y");
    let z = var("z");
    let x = var("x");

    let dy = Expr::Derivative(Arc::new(y.clone()), "x".to_string());
    let dz = Expr::Derivative(Arc::new(z.clone()), "x".to_string());

    let eq1 = Expr::Eq(Arc::new(dy), Arc::new(x.clone()));
    let eq2 = Expr::Eq(Arc::new(dz), Arc::new(y.clone()));

    let eqs = vec![eq1, eq2];
    let funcs = vec!["y", "z"];

    let solutions = solve_ode_system(&eqs, &funcs, "x");
    assert!(solutions.is_some());
    let sols = solutions.unwrap();
    assert_eq!(sols.len(), 2);
    println!("System Solutions: y={}, z={}", sols[0], sols[1]);
}

#[test]
fn test_solve_ode_dispatcher() {
    // Test that solve_ode correctly dispatches to linear solver
    // y' + y = 0 => y = C * e^-x
    let y = var("y");
    let x = var("x");
    let y_prime = Expr::Derivative(Arc::new(y.clone()), "x".to_string());
    let eq = Expr::new_add(y_prime, y.clone());

    let sol = solve_ode(&eq, "y", "x", None);
    println!("Dispatcher Solution: {}", sol);

    // Check if solution is valid (roughly)
    // It should be y = C * exp(-x)
    // The output might be implicit or explicit depending on solver picked.
    // Linear solver returns explicit y = ...
}
