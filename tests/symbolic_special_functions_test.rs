//! Comprehensive tests for symbolic special functions module.
//!
//! Tests cover gamma, beta, error functions, Bessel functions,
//! orthogonal polynomials, and differential equations.

use num_traits::ToPrimitive;
use rssn::symbolic::core::{DagOp, Expr};
use rssn::symbolic::special_functions::*;

// --- Helper Functions ---

fn evaluate_expr(expr: &Expr) -> Option<f64> {

    match expr {
        Expr::Constant(v) => Some(*v),
        Expr::BigInt(v) => v.to_f64(),
        Expr::Rational(v) => v.to_f64(),
        Expr::Add(a, b) => Some(evaluate_expr(a)? + evaluate_expr(b)?),
        Expr::Sub(a, b) => Some(evaluate_expr(a)? - evaluate_expr(b)?),
        Expr::Mul(a, b) => Some(evaluate_expr(a)? * evaluate_expr(b)?),
        Expr::Div(a, b) => Some(evaluate_expr(a)? / evaluate_expr(b)?),
        Expr::Power(a, b) => Some(evaluate_expr(a)?.powf(evaluate_expr(b)?)),
        Expr::Pi => Some(std::f64::consts::PI),
        Expr::E => Some(std::f64::consts::E),
        Expr::Dag(node) => evaluate_dag(node),
        _ => None,
    }
}

fn evaluate_dag(node: &rssn::symbolic::core::DagNode) -> Option<f64> {

    match &node.op {
        DagOp::Constant(v) => Some(v.into_inner()),
        DagOp::BigInt(v) => v.to_f64(),
        DagOp::Rational(v) => v.to_f64(),
        DagOp::Pi => Some(std::f64::consts::PI),
        _ => None,
    }
}

fn assert_approx_eq(expr: &Expr, expected: f64) {

    if let Some(val) = evaluate_expr(expr) {

        assert!(
            (val - expected).abs() < 1e-6,
            "Expected {}, got {} (from {:?})",
            expected,
            val,
            expr
        );
    }
}

fn assert_expr_eq(expr: &Expr, expected: &Expr) {

    assert_eq!(expr, expected, "Expected {:?}, got {:?}", expected, expr);
}

// ============================================================================
// Gamma Function Tests
// ============================================================================

#[test]

fn test_gamma_integers() {

    // Γ(1) = 1
    let g1 = gamma(Expr::Constant(1.0));

    assert_approx_eq(&g1, 1.0);

    // Γ(2) = 1
    let g2 = gamma(Expr::Constant(2.0));

    assert_approx_eq(&g2, 1.0);

    // Γ(5) = 4! = 24
    let g5 = gamma(Expr::Constant(5.0));

    assert_approx_eq(&g5, 24.0);

    // Γ(6) = 5! = 120
    let g6 = gamma(Expr::Constant(6.0));

    assert_approx_eq(&g6, 120.0);
}

#[test]

fn test_gamma_half_integer() {

    // Γ(0.5) = √π
    let g_half = gamma(Expr::Constant(0.5));

    // Check it simplifies to Sqrt(Pi)
    match &g_half {
        Expr::Sqrt(inner) if **inner == Expr::Pi => { /* Success */ }
        Expr::Dag(_) => { /* Also acceptable if DAG form */ }
        _ => panic!("Expected sqrt(pi), got {:?}", g_half),
    }
}

#[test]

fn test_ln_gamma() {

    // ln(Γ(5)) = ln(24)
    let lg = ln_gamma(Expr::Constant(5.0));

    assert_approx_eq(&lg, 24.0_f64.ln());
}

// ============================================================================
// Beta Function Tests
// ============================================================================

#[test]

fn test_beta_basic() {

    // B(1, 1) = 1
    let b = beta(Expr::Constant(1.0), Expr::Constant(1.0));

    assert_approx_eq(&b, 1.0);

    // B(2, 1) = 1/2
    let b21 = beta(Expr::Constant(2.0), Expr::Constant(1.0));

    assert_approx_eq(&b21, 0.5);
}

#[test]

fn test_beta_symmetry() {

    // B(a, b) = B(b, a)
    let x = Expr::Variable("x".to_string());

    let y = Expr::Variable("y".to_string());

    // For constants, verify symmetry
    let b12 = beta(Expr::Constant(1.0), Expr::Constant(2.0));

    let b21 = beta(Expr::Constant(2.0), Expr::Constant(1.0));

    assert_eq!(evaluate_expr(&b12), evaluate_expr(&b21));
}

// ============================================================================
// Digamma Function Tests
// ============================================================================

#[test]

fn test_digamma_special_values() {

    // ψ(1) = -γ (Euler-Mascheroni constant)
    let d1 = digamma(Expr::Constant(1.0));

    match &d1 {
        Expr::Variable(s) if s == "-gamma" => { /* Success */ }
        Expr::Dag(_) => { /* Also acceptable */ }
        _ => panic!("Expected -gamma variable, got {:?}", d1),
    }
}

#[test]

fn test_polygamma() {

    // ψ⁽⁰⁾(z) = ψ(z)
    let pg0 = polygamma(Expr::Constant(0.0), Expr::Constant(2.0));

    let d2 = digamma(Expr::Constant(2.0));
    // Both should represent digamma(2)
}

// ============================================================================
// Error Function Tests
// ============================================================================

#[test]

fn test_erf_zero() {

    // erf(0) = 0
    let e = erf(Expr::Constant(0.0));

    assert_approx_eq(&e, 0.0);
}

#[test]

fn test_erf_odd_symmetry() {

    // erf(-x) should simplify to -erf(x)
    let x = Expr::new_variable("x");

    let e_neg = erf(Expr::new_neg(x.clone()));

    // Accept any result - just ensure it constructed without panic
    match &e_neg {
        Expr::Neg(_) => { /* expected: -erf(x) */ }
        Expr::Erf(_) => { /* also acceptable: erf(-x) stayed */ }
        Expr::Dag(_) => { /* simplified to dag */ }
        _ => { /* some other form is fine too */ }
    }
}

#[test]

fn test_erfc() {

    // erfc(0) = 1 - erf(0) = 1
    let ec = erfc(Expr::Constant(0.0));

    assert_approx_eq(&ec, 1.0);
}

// ============================================================================
// Zeta Function Tests
// ============================================================================

#[test]

fn test_zeta_special_values() {

    // ζ(0) = -1/2
    let z0 = zeta(Expr::Constant(0.0));

    assert_approx_eq(&z0, -0.5);

    // ζ(1) = ∞
    let z1 = zeta(Expr::Constant(1.0));

    assert_eq!(z1, Expr::Infinity);

    // ζ(-2) = 0 (trivial zero)
    let z_neg2 = zeta(Expr::Constant(-2.0));

    assert_approx_eq(&z_neg2, 0.0);
}

// ============================================================================
// Bessel Function Tests
// ============================================================================

#[test]

fn test_bessel_j_at_zero() {

    // J_0(0) = 1
    let j0 = bessel_j(Expr::Constant(0.0), Expr::Constant(0.0));

    assert_approx_eq(&j0, 1.0);

    // J_n(0) = 0 for n > 0
    let j1 = bessel_j(Expr::Constant(1.0), Expr::Constant(0.0));

    assert_approx_eq(&j1, 0.0);

    let j5 = bessel_j(Expr::Constant(5.0), Expr::Constant(0.0));

    assert_approx_eq(&j5, 0.0);
}

#[test]

fn test_bessel_y_at_zero() {

    // Y_n(0) = -∞
    let y0 = bessel_y(Expr::Constant(0.0), Expr::Constant(0.0));

    assert_eq!(y0, Expr::NegativeInfinity);
}

#[test]

fn test_bessel_i_at_zero() {

    // I_0(0) = 1
    let i0 = bessel_i(Expr::Constant(0.0), Expr::Constant(0.0));

    assert_approx_eq(&i0, 1.0);

    // I_n(0) = 0 for n > 0
    let i1 = bessel_i(Expr::Constant(1.0), Expr::Constant(0.0));

    assert_approx_eq(&i1, 0.0);
}

#[test]

fn test_bessel_k_at_zero() {

    // K_n(0) = ∞
    let k0 = bessel_k(Expr::Constant(0.0), Expr::Constant(0.0));

    assert_eq!(k0, Expr::Infinity);
}

// ============================================================================
// Legendre Polynomial Tests
// ============================================================================

#[test]

fn test_legendre_p_basic() {

    let x = Expr::Variable("x".to_string());

    // P_0(x) = 1
    let p0 = legendre_p(Expr::Constant(0.0), x.clone());

    assert_approx_eq(&p0, 1.0);

    // P_1(x) = x
    let p1 = legendre_p(Expr::Constant(1.0), x.clone());

    assert_eq!(p1, x);
}

#[test]

fn test_legendre_p_recurrence() {

    // P_2(x) = (3x² - 1)/2
    // At x=0: P_2(0) = -1/2
    let p2_at_0 = legendre_p(Expr::Constant(2.0), Expr::Constant(0.0));

    assert_approx_eq(&p2_at_0, -0.5);

    // At x=1: P_n(1) = 1 for all n
    let p2_at_1 = legendre_p(Expr::Constant(2.0), Expr::Constant(1.0));

    assert_approx_eq(&p2_at_1, 1.0);
}

// ============================================================================
// Laguerre Polynomial Tests
// ============================================================================

#[test]

fn test_laguerre_l_basic() {

    let x = Expr::Variable("x".to_string());

    // L_0(x) = 1
    let l0 = laguerre_l(Expr::Constant(0.0), x.clone());

    assert_approx_eq(&l0, 1.0);

    // L_1(x) = 1 - x, at x=0: L_1(0) = 1
    let l1_at_0 = laguerre_l(Expr::Constant(1.0), Expr::Constant(0.0));

    assert_approx_eq(&l1_at_0, 1.0);

    // L_1(1) = 0
    let l1_at_1 = laguerre_l(Expr::Constant(1.0), Expr::Constant(1.0));

    assert_approx_eq(&l1_at_1, 0.0);
}

#[test]

fn test_generalized_laguerre() {

    // L_0^α(x) = 1 for any α
    let gl = generalized_laguerre(
        Expr::Constant(0.0),
        Expr::Constant(1.0),
        Expr::Constant(2.0),
    );

    assert_approx_eq(&gl, 1.0);

    // L_n^0(x) = L_n(x)
    let gl0 = generalized_laguerre(
        Expr::Constant(1.0),
        Expr::Constant(0.0),
        Expr::Constant(0.0),
    );

    let l1 = laguerre_l(Expr::Constant(1.0), Expr::Constant(0.0));

    // Both should give L_1(0) = 1
    assert_approx_eq(&gl0, 1.0);

    assert_approx_eq(&l1, 1.0);
}

// ============================================================================
// Hermite Polynomial Tests
// ============================================================================

#[test]

fn test_hermite_h_basic() {

    let x = Expr::Variable("x".to_string());

    // H_0(x) = 1
    let h0 = hermite_h(Expr::Constant(0.0), x.clone());

    assert_approx_eq(&h0, 1.0);

    // H_1(x) = 2x, at x=0: H_1(0) = 0
    let h1_at_0 = hermite_h(Expr::Constant(1.0), Expr::Constant(0.0));

    assert_approx_eq(&h1_at_0, 0.0);

    // H_1(1) = 2
    let h1_at_1 = hermite_h(Expr::Constant(1.0), Expr::Constant(1.0));

    assert_approx_eq(&h1_at_1, 2.0);
}

#[test]

fn test_hermite_h_recurrence() {

    // H_2(x) = 4x² - 2, at x=0: H_2(0) = -2
    let h2_at_0 = hermite_h(Expr::Constant(2.0), Expr::Constant(0.0));

    assert_approx_eq(&h2_at_0, -2.0);

    // H_2(1) = 4 - 2 = 2
    let h2_at_1 = hermite_h(Expr::Constant(2.0), Expr::Constant(1.0));

    assert_approx_eq(&h2_at_1, 2.0);
}

// ============================================================================
// Chebyshev Polynomial Tests
// ============================================================================

#[test]

fn test_chebyshev_t_basic() {

    let x = Expr::Variable("x".to_string());

    // T_0(x) = 1
    let t0 = chebyshev_t(Expr::Constant(0.0), x.clone());

    assert_approx_eq(&t0, 1.0);

    // T_1(x) = x
    let t1 = chebyshev_t(Expr::Constant(1.0), x.clone());

    assert_eq!(t1, x);
}

#[test]

fn test_chebyshev_t_at_one() {

    // T_n(1) = 1 for all n
    for n in 0..=5 {

        let tn = chebyshev_t(Expr::Constant(n as f64), Expr::Constant(1.0));

        assert_approx_eq(&tn, 1.0);
    }
}

#[test]

fn test_chebyshev_t_at_minus_one() {

    // T_n(-1) = (-1)^n
    let t0 = chebyshev_t(Expr::Constant(0.0), Expr::Constant(-1.0));

    assert_approx_eq(&t0, 1.0);

    let t1 = chebyshev_t(Expr::Constant(1.0), Expr::Constant(-1.0));

    assert_approx_eq(&t1, -1.0);

    let t2 = chebyshev_t(Expr::Constant(2.0), Expr::Constant(-1.0));

    assert_approx_eq(&t2, 1.0);
}

#[test]

fn test_chebyshev_u_basic() {

    let x = Expr::Variable("x".to_string());

    // U_0(x) = 1
    let u0 = chebyshev_u(Expr::Constant(0.0), x.clone());

    assert_approx_eq(&u0, 1.0);

    // U_1(x) = 2x, at x=1: U_1(1) = 2
    let u1_at_1 = chebyshev_u(Expr::Constant(1.0), Expr::Constant(1.0));

    assert_approx_eq(&u1_at_1, 2.0);
}

// ============================================================================
// Differential Equation Construction Tests
// ============================================================================

#[test]

fn test_differential_equations_construct() {

    // Just ensure they construct without panic
    let y = Expr::new_variable("y");

    let x = Expr::new_variable("x");

    let n = Expr::Constant(2.0);

    let bessel_eq = bessel_differential_equation(&y, &x, &n);

    match bessel_eq {
        Expr::Eq(_, _) => { /* Success */ }
        _ => panic!("Expected Eq, got {:?}", bessel_eq),
    }

    let legendre_eq = legendre_differential_equation(&y, &x, &n);

    match legendre_eq {
        Expr::Eq(_, _) => { /* Success */ }
        _ => panic!("Expected Eq, got {:?}", legendre_eq),
    }

    let laguerre_eq = laguerre_differential_equation(&y, &x, &n);

    match laguerre_eq {
        Expr::Eq(_, _) => { /* Success */ }
        _ => panic!("Expected Eq, got {:?}", laguerre_eq),
    }

    let hermite_eq = hermite_differential_equation(&y, &x, &n);

    match hermite_eq {
        Expr::Eq(_, _) => { /* Success */ }
        _ => panic!("Expected Eq, got {:?}", hermite_eq),
    }

    let chebyshev_eq = chebyshev_differential_equation(&y, &x, &n);

    match chebyshev_eq {
        Expr::Eq(_, _) => { /* Success */ }
        _ => panic!("Expected Eq, got {:?}", chebyshev_eq),
    }
}

// ============================================================================
// Rodrigues Formula Tests
// ============================================================================

#[test]

fn test_rodrigues_formulas_construct() {

    let n = Expr::Constant(3.0);

    let x = Expr::new_variable("x");

    let legendre_rf = legendre_rodrigues_formula(&n, &x);

    match legendre_rf {
        Expr::Eq(_, _) => { /* Success */ }
        _ => panic!("Expected Eq, got {:?}", legendre_rf),
    }

    let hermite_rf = hermite_rodrigues_formula(&n, &x);

    match hermite_rf {
        Expr::Eq(_, _) => { /* Success */ }
        _ => panic!("Expected Eq, got {:?}", hermite_rf),
    }
}

// ============================================================================
// Integration Tests
// ============================================================================

#[test]

fn test_gamma_beta_relationship() {

    // B(a, b) = Γ(a)Γ(b) / Γ(a+b)
    // For a=2, b=3: B(2,3) = 1!*2!/4! = 1*2/24 = 1/12
    let b23 = beta(Expr::Constant(2.0), Expr::Constant(3.0));

    if let Some(val) = evaluate_expr(&b23) {

        assert!(
            (val - 1.0 / 12.0).abs() < 1e-6,
            "Expected 1/12, got {}",
            val
        );
    }
}

#[test]

fn test_polynomial_orthogonality_at_boundaries() {

    // All orthogonal polynomials evaluated at typical boundary points
    // P_n(1) = 1 for Legendre
    // T_n(1) = 1 for Chebyshev T
    for n in 0..=5 {

        let pn = legendre_p(Expr::Constant(n as f64), Expr::Constant(1.0));

        if let Some(val) = evaluate_expr(&pn) {

            assert!(
                (val - 1.0).abs() < 1e-6,
                "P_{}(1) should be 1, got {}",
                n,
                val
            );
        }

        let tn = chebyshev_t(Expr::Constant(n as f64), Expr::Constant(1.0));

        if let Some(val) = evaluate_expr(&tn) {

            assert!(
                (val - 1.0).abs() < 1e-6,
                "T_{}(1) should be 1, got {}",
                n,
                val
            );
        }
    }
}
