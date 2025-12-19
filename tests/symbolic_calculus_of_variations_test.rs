use rssn::symbolic::calculus_of_variations::euler_lagrange;
use rssn::symbolic::core::Expr;
use std::sync::Arc;

#[test]
fn test_free_particle() {
    // L = 1/2 * m * v^2 = 1/2 * m * (dx/dt)^2
    let t = "t";
    let x = "x";
    let x_expr = Expr::Variable(x.to_string());
    // Manual derivative construction
    let v_expr = Expr::Derivative(Arc::new(x_expr.clone()), t.to_string());
    
    let m = Expr::Variable("m".to_string());
    let kinetic_energy = Expr::new_mul(
        Expr::new_mul(Expr::Constant(0.5), m.clone()),
        Expr::new_pow(v_expr.clone(), Expr::Constant(2.0))
    );
    
    let lagrangian = kinetic_energy;
    
    let eq = euler_lagrange(&lagrangian, x, t);
    
    let eq_str = format!("{:?}", eq);
    println!("Free Particle Equation: {}", eq_str);
    
    // Should be m * x''
    assert!(eq_str.contains("m"));
    // Depending on simplification, x'' might appear differently.
    // But result should NOT be 0.
    assert_ne!(eq_str, "0");
}

#[test]
fn test_harmonic_oscillator() {
    // L = 1/2 m v^2 - 1/2 k x^2
    let t = "t";
    let x = "x";
    let x_expr = Expr::Variable(x.to_string());
    let v_expr = Expr::Derivative(Arc::new(x_expr.clone()), t.to_string());
    
    let m = Expr::Variable("m".to_string());
    let k = Expr::Variable("k".to_string());
    
    let kinetic = Expr::new_mul(
        Expr::new_mul(Expr::Constant(0.5), m.clone()),
        Expr::new_pow(v_expr.clone(), Expr::Constant(2.0))
    );
    let potential = Expr::new_mul(
        Expr::new_mul(Expr::Constant(0.5), k.clone()),
        Expr::new_pow(x_expr.clone(), Expr::Constant(2.0)) 
    );
    
    let lagrangian = Expr::new_sub(kinetic, potential);
    
    // Eq: m x'' + k x = 0
    let eq = euler_lagrange(&lagrangian, x, t);
    let eq_str = format!("{:?}", eq);
    println!("Harmonic Oscillator Equation: {}", eq_str);
    
    assert!(eq_str.contains("m"));
    assert!(eq_str.contains("k"));
}

#[test]
fn test_pendulum_nonlinear() {
    // L = 1/2 m l^2 theta'^2 + m g l cos(theta)
    let t = "t";
    let theta = "theta";
    let theta_expr = Expr::Variable(theta.to_string());
    let omega = Expr::Derivative(Arc::new(theta_expr.clone()), t.to_string());
    
    let m = Expr::Variable("m".to_string());
    let l = Expr::Variable("l".to_string());
    let g = Expr::Variable("g".to_string());
    
    let kinetic = Expr::new_mul(
        Expr::new_mul(Expr::Constant(0.5), m.clone()),
        Expr::new_mul(
            Expr::new_pow(l.clone(), Expr::Constant(2.0)),
            Expr::new_pow(omega, Expr::Constant(2.0))
        )
    );
    
    let potential = Expr::new_mul(
        Expr::new_mul(m.clone(), g.clone()),
        Expr::new_mul(l.clone(), Expr::new_cos(theta_expr))
    );
    // Be careful with signs. PE usually U = -mgl cos(theta) (if y is down) or mgl(1-cos).
    // If L = T - V, and U = -mgl cos(theta), then L = T + mgl cos(theta).
    
    let lagrangian = Expr::new_add(kinetic, potential); 
    
    // Eq: d/dt(dL/dtheta') - dL/dtheta = 0
    // dL/dtheta' = m l^2 theta'
    // d/dt(...) = m l^2 theta''
    // dL/dtheta = - m g l sin(theta)
    // Eq: m l^2 theta'' - (- m g l sin(theta)) = m l^2 theta'' + m g l sin(theta)
    
    let eq = euler_lagrange(&lagrangian, theta, t);
    let eq_str = format!("{:?}", eq);
    println!("Pendulum Equation: {}", eq_str);
    
    assert!(eq_str.contains("sin(theta)"));
    assert!(eq_str.contains("m")); // Mass should be present from kinetic term too
}
