use rssn::symbolic::special::*;

fn assert_approx_eq(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-6, "Expected {}, got {}", b, a);
}

#[test]
fn test_gamma_numerical() {
    // Gamma(5) = 4! = 24
    assert_approx_eq(gamma_numerical(5.0), 24.0);
    // Gamma(0.5) = sqrt(pi)
    assert_approx_eq(gamma_numerical(0.5), std::f64::consts::PI.sqrt());
}

#[test]
fn test_ln_gamma_numerical() {
    // ln(Gamma(5)) = ln(24)
    assert_approx_eq(ln_gamma_numerical(5.0), 24.0f64.ln());
}

#[test]
fn test_beta_numerical() {
    // B(x, y) = Gamma(x)Gamma(y) / Gamma(x+y)
    // B(1, 1) = 1 * 1 / 1 = 1
    assert_approx_eq(beta_numerical(1.0, 1.0), 1.0);
    // B(2, 2) = 1 * 1 / Gamma(4) = 1 / 6
    assert_approx_eq(beta_numerical(2.0, 2.0), 1.0 / 6.0);
}

#[test]
fn test_erf_numerical() {
    // erf(0) = 0
    assert_approx_eq(erf_numerical(0.0), 0.0);
    // erf(inf) = 1 (approx)
    // erf(1) ~ 0.8427
    assert_approx_eq(erf_numerical(1.0), 0.84270079294971);
}

#[test]
fn test_erfc_numerical() {
    // erfc(0) = 1
    assert_approx_eq(erfc_numerical(0.0), 1.0);
    // erfc(1) = 1 - erf(1)
    assert_approx_eq(erfc_numerical(1.0), 1.0 - 0.84270079294971);
}
