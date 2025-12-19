use rssn::symbolic::classical_mechanics::*;
use rssn::symbolic::core::Expr;
use rssn::symbolic::vector::Vector;

#[test]
fn test_kinetic_energy() {
    let m = Expr::new_variable("m");
    let v = Expr::new_variable("v");
    let t = kinetic_energy(&m, &v);
    let t_str = format!("{:?}", t);
    assert!(t_str.contains("0.5"));
    assert!(t_str.contains("m"));
    assert!(t_str.contains("v^(2)"));
}

#[test]
fn test_potential_energy_gravity_uniform() {
    let m = Expr::new_variable("m");
    let g = Expr::new_variable("g");
    let h = Expr::new_variable("h");
    let v = potential_energy_gravity_uniform(&m, &h, &g);
    let v_str = format!("{:?}", v);
    assert!(v_str.contains("m"));
    assert!(v_str.contains("g"));
    assert!(v_str.contains("h"));
}

#[test]
fn test_spring_potential() {
    let k = Expr::new_variable("k");
    let x = Expr::new_variable("x");
    let v = potential_energy_spring(&k, &x);
    let v_str = format!("{:?}", v);
    assert!(v_str.contains("0.5"));
    assert!(v_str.contains("k"));
    assert!(v_str.contains("x^(2)"));
}

#[test]
fn test_euler_lagrange_harmonic_oscillator() {
    let m = Expr::new_variable("m");
    let k = Expr::new_variable("k");
    let x = Expr::new_variable("x");
    let x_dot = Expr::new_variable("x_dot");
    
    let t = kinetic_energy(&m, &x_dot);
    let v = potential_energy_spring(&k, &x);
    let l = lagrangian(&t, &v);
    
    // d/dt (dL/dx_dot) - dL/dx = 0
    // dL/dx_dot = m * x_dot
    // d/dt(dL/dx_dot) = m * d/dt(x_dot)
    // dL/dx = -k * x
    // Equation: m * d/dt(x_dot) + k * x = 0
    
    let eq = euler_lagrange_equation(&l, "x", "x_dot", "t");
    let eq_str = format!("{:?}", eq);
    
    assert!(eq_str.contains("m"));
    assert!(eq_str.contains("d/dt(d/dt(x))"));
    assert!(eq_str.contains("k"));
    assert!(eq_str.contains("x"));
}

#[test]
fn test_torque_and_angular_momentum() {
    let rx = Expr::new_variable("rx");
    let ry = Expr::new_variable("ry");
    let rz = Expr::new_variable("rz");
    let r = Vector::new(rx, ry, rz);
    
    let fx = Expr::new_variable("fx");
    let fy = Expr::new_variable("fy");
    let fz = Expr::new_variable("fz");
    let f = Vector::new(fx, fy, fz);
    
    let tau = torque(&r, &f);
    // x component should be ry*fz - rz*fy
    let tau_x_str = format!("{:?}", tau.x);
    assert!(tau_x_str.contains("ry"));
    assert!(tau_x_str.contains("fz"));
    assert!(tau_x_str.contains("rz"));
    assert!(tau_x_str.contains("fy"));
}
