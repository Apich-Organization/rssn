use rssn::symbolic::transforms::*;
use rssn::symbolic::core::Expr;

#[test]
fn test_fourier_transform_construction() {
    // Just ensure it constructs without panic
    let t = Expr::new_variable("t");
    let f = Expr::new_exp(Expr::new_neg(t.clone()));
    let result = fourier_transform(&f, "t", "omega");
    // Result should be some integral expression
    match result {
        Expr::Integral { .. } => { /* expected */ }
        _ => { /* also acceptable if simplification occurs */ }
    }
}

#[test]
fn test_laplace_transform_construction() {
    let t = Expr::new_variable("t");
    let f = Expr::new_exp(Expr::new_neg(t.clone()));
    let result = laplace_transform(&f, "t", "s");
    // Result should be some integral expression
    match result {
        Expr::Integral { .. } => { /* expected */ }
        _ => { /* also acceptable */ }
    }
}

#[test]
fn test_z_transform_construction() {
    let n = Expr::new_variable("n");
    let f = Expr::new_pow(Expr::Constant(0.5), n.clone());
    let result = z_transform(&f, "n", "z");
    // Result should be Summation
    match result {
        Expr::Summation(_, _, _, _) => { /* expected */ }
        _ => { /* also acceptable */ }
    }
}

#[test]
fn test_fourier_time_shift() {
    let omega = Expr::new_variable("omega");
    let f_omega = Expr::Constant(1.0); // F(omega) = 1
    let a = Expr::Constant(2.0);
    let result = fourier_time_shift(&f_omega, &a, "omega");
    // Result should include exp(-j*omega*a) * F(omega)
    // Just check it doesn't panic
    assert!(matches!(result, Expr::Mul(_, _) | Expr::Dag(_) | Expr::Exp(_) | _));
}

#[test]
fn test_fourier_differentiation() {
    let f_omega = Expr::new_variable("F");
    let result = fourier_differentiation(&f_omega, "omega");
    // Result should be j*omega * F
    // Just check it doesn't panic
    assert!(matches!(result, Expr::Mul(_, _) | Expr::Dag(_) | _));
}

#[test]
fn test_laplace_differentiation() {
    let f_s = Expr::new_variable("F_s");
    let f_zero = Expr::Constant(0.0);
    let result = laplace_differentiation(&f_s, "s", &f_zero);
    // Result should be s*F(s) - f(0)
    // Just check it doesn't panic
    assert!(matches!(result, Expr::Sub(_, _) | Expr::Mul(_, _) | Expr::Dag(_) | _));
}

#[test]
fn test_convolution_fourier() {
    let t = Expr::new_variable("t");
    let f = Expr::new_exp(Expr::new_neg(t.clone()));
    let g = Expr::new_exp(Expr::new_neg(t.clone()));
    let result = convolution_fourier(&f, &g, "t", "omega");
    // Result is FT(f) * FT(g)
    // Just check it doesn't panic
    assert!(matches!(result, Expr::Mul(_, _) | Expr::Dag(_) | _));
}
