use rssn::symbolic::core::Expr;
use rssn::symbolic::quantum_field_theory::*;

#[test]
fn test_dirac_adjoint() {
    let psi = Expr::new_variable("psi");
    let bar = dirac_adjoint(&psi);
    assert!(bar.to_string().contains("psi"));
    assert!(bar.to_string().contains("gamma_0"));
}

#[test]
fn test_feynman_slash() {
    let v = Expr::new_variable("v_mu");
    let slash = feynman_slash(&v);
    assert!(slash.to_string().contains("gamma_mu"));
    assert!(slash.to_string().contains("v_mu"));
}

#[test]
fn test_scalar_field_lagrangian() {
    let phi = Expr::new_variable("phi");
    let m = Expr::new_variable("m");
    let l = scalar_field_lagrangian(&phi, &m);
    let l_str = l.to_string();
    assert!(l_str.contains("phi"));
    assert!(l_str.contains("m"));
    assert!(l_str.contains("partial_mu_phi"));
}

#[test]
fn test_qed_lagrangian() {
    let psi_bar = Expr::new_variable("psi_bar");
    let psi = Expr::new_variable("psi");
    let a_mu = Expr::new_variable("A_mu");
    let m = Expr::new_variable("m");
    let e = Expr::new_variable("e");
    let l = qed_lagrangian(&psi_bar, &psi, &a_mu, &m, &e);
    let l_str = l.to_string();
    assert!(l_str.contains("psi_bar"));
    assert!(l_str.contains("psi"));
    assert!(l_str.contains("gamma_mu"));
    assert!(l_str.contains("m"));
}

#[test]
fn test_propagator() {
    let p = Expr::new_variable("p");
    let m = Expr::new_variable("m");
    let prop_scalar = propagator(&p, &m, false);
    assert!(prop_scalar.to_string().contains("p"));
    assert!(prop_scalar.to_string().contains("m"));
    
    let prop_fermion = propagator(&p, &m, true);
    assert!(prop_fermion.to_string().contains("gamma_mu"));
}
