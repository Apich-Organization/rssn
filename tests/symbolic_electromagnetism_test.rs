use rssn::symbolic::core::Expr;
use rssn::symbolic::electromagnetism::*;
use rssn::symbolic::vector::Vector;

#[test]

fn test_lorentz_force() {

    let q = Expr::new_variable("q");

    let ex = Expr::new_variable("Ex");

    let ey = Expr::new_variable("Ey");

    let ez = Expr::new_variable("Ez");

    let e_field = Vector::new(ex, ey, ez);

    let vx = Expr::new_variable("vx");

    let vy = Expr::new_variable("vy");

    let vz = Expr::new_variable("vz");

    let velocity = Vector::new(vx, vy, vz);

    let bx = Expr::new_variable("bx");

    let by = Expr::new_variable("by");

    let bz = Expr::new_variable("bz");

    let b_field = Vector::new(bx, by, bz);

    let force = lorentz_force(
        &q,
        &e_field,
        &velocity,
        &b_field,
    );

    let fx_str = format!("{:?}", force.x);

    // Fx = q * (Ex + vy*bz - vz*by)
    assert!(fx_str.contains("q"));

    assert!(fx_str.contains("Ex"));

    assert!(fx_str.contains("vy"));

    assert!(fx_str.contains("bz"));
}

#[test]

fn test_maxwell_equations_gauss() {

    let ex = Expr::new_variable("x"); // E = [x, 0, 0] -> div E = 1
    let e_field = Vector::new(
        ex,
        Expr::Constant(0.0),
        Expr::Constant(0.0),
    );

    let b_field = Vector::new(
        Expr::Constant(0.0),
        Expr::Constant(0.0),
        Expr::Constant(0.0),
    );

    let rho = Expr::new_variable("rho");

    let j_field = Vector::new(
        Expr::Constant(0.0),
        Expr::Constant(0.0),
        Expr::Constant(0.0),
    );

    let maxwell = MaxwellEquations::new(
        &e_field,
        &b_field,
        &rho,
        &j_field,
        ("x", "y", "z"),
        "t",
    );

    let gauss_str = format!(
        "{:?}",
        maxwell.gauss_law_electric
    );

    // (div E - rho/epsilon_0) = (1 - rho/epsilon_0)
    assert!(gauss_str.contains("1"));

    assert!(gauss_str.contains("rho"));

    assert!(gauss_str.contains("epsilon_0"));
}

#[test]

fn test_energy_density() {

    let e_field = Vector::new(
        Expr::new_variable("Ex"),
        Expr::Constant(0.0),
        Expr::Constant(0.0),
    );

    let b_field = Vector::new(
        Expr::Constant(0.0),
        Expr::new_variable("By"),
        Expr::Constant(0.0),
    );

    let u = energy_density(&e_field, &b_field);

    let u_str = format!("{:?}", u);

    assert!(u_str.contains("epsilon_0"));

    assert!(u_str.contains("mu_0"));

    assert!(u_str.contains("Ex"));

    assert!(u_str.contains("By"));
}

#[test]

fn test_coulombs_law() {

    let q = Expr::new_variable("Q");

    let rx = Expr::new_variable("rx");

    let ry = Expr::new_variable("ry");

    let rz = Expr::new_variable("rz");

    let r = Vector::new(rx, ry, rz);

    let e = coulombs_law(&q, &r);

    let ex_str = format!("{:?}", e.x);

    assert!(ex_str.contains("Q"));

    assert!(ex_str.contains("rx"));

    assert!(ex_str.contains("epsilon_0"));

    assert!(ex_str.contains("Pi"));
}
