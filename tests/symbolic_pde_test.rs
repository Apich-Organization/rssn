
use rssn::symbolic::core::Expr;
use rssn::symbolic::pde::{
    solve_pde_by_characteristics, solve_wave_equation_1d_dalembert,
    solve_heat_equation_1d, solve_laplace_equation_2d
};
use std::sync::Arc;

fn var(name: &str) -> Expr {
    Expr::Variable(name.to_string())
}

#[test]
fn test_method_of_characteristics() {
    // Solve a u_x + b u_y = c
    // Let's try: u_x + u_y = 1
    // a=1, b=1, c=1
    let u = var("u");
    
    let u_x = Expr::Derivative(Arc::new(u.clone()), "x".to_string());
    let u_y = Expr::Derivative(Arc::new(u.clone()), "y".to_string());
    
    // u_x + u_y - 1 = 0
    let lhs = Expr::new_add(u_x, u_y);
    let eq = Expr::new_sub(lhs, Expr::Constant(1.0));
    
    let sol = solve_pde_by_characteristics(&eq, "u", &["x", "y"]);
    assert!(sol.is_some());
    println!("Characteristics Solution: {}", sol.unwrap());
}

#[test]
fn test_wave_equation_dalembert() {
    // u_tt = c^2 * u_xx, with c = 2
    // u_tt - 4 * u_xx = 0
    let u = var("u");
    let _t = var("t");
    let _x = var("x");
    
    let u_t = Expr::Derivative(Arc::new(u.clone()), "t".to_string());
    let u_tt = Expr::Derivative(Arc::new(u_t), "t".to_string());
    
    let u_x = Expr::Derivative(Arc::new(u.clone()), "x".to_string());
    let u_xx = Expr::Derivative(Arc::new(u_x), "x".to_string());
    
    // u_tt - 4*u_xx = 0 (c = 2)
    let rhs = Expr::new_mul(Expr::Constant(4.0), u_xx);
    let eq = Expr::new_sub(u_tt, rhs);
    
    let sol = solve_wave_equation_1d_dalembert(&eq, "u", &["t", "x"]);
    assert!(sol.is_some());
    println!("D'Alembert Solution: {}", sol.unwrap());
}

#[test]
fn test_heat_equation_1d() {
    // u_t = α*u_xx, with α = 0.5
    // u_t - 0.5*u_xx = 0
    let u = var("u");
    
    let u_t = Expr::Derivative(Arc::new(u.clone()), "t".to_string());
    
    let u_x = Expr::Derivative(Arc::new(u.clone()), "x".to_string());
    let u_xx = Expr::Derivative(Arc::new(u_x), "x".to_string());
    
    // u_t - 0.5*u_xx = 0 (α = 0.5)
    let rhs = Expr::new_mul(Expr::Constant(0.5), u_xx);
    let eq = Expr::new_sub(u_t, rhs);
    
    let sol = solve_heat_equation_1d(&eq, "u", &["t", "x"]);
    assert!(sol.is_some());
    let solution = sol.unwrap();
    println!("Heat Equation Solution: {}", solution);
    
    // Verify it's an equation
    assert!(matches!(solution, Expr::Eq(_, _)));
}

#[test]
fn test_laplace_equation_2d() {
    // u_xx + u_yy = 0
    let u = var("u");
    
    let u_x = Expr::Derivative(Arc::new(u.clone()), "x".to_string());
    let u_xx = Expr::Derivative(Arc::new(u_x), "x".to_string());
    
    let u_y = Expr::Derivative(Arc::new(u.clone()), "y".to_string());
    let u_yy = Expr::Derivative(Arc::new(u_y), "y".to_string());
    
    // u_xx + u_yy = 0
    let eq = Expr::new_add(u_xx, u_yy);
    
    let sol = solve_laplace_equation_2d(&eq, "u", &["x", "y"]);
    assert!(sol.is_some());
    let solution = sol.unwrap();
    println!("Laplace Equation Solution: {}", solution);
    
    // Verify it's an equation
    assert!(matches!(solution, Expr::Eq(_, _)));
}

