//! Test suite for special functions module.

use rssn::symbolic::special::*;

fn assert_approx_eq(
    a: f64,
    b: f64,
) {

    assert!(
        (a - b).abs() < 1e-6,
        "Expected {}, got {}",
        b,
        a
    );
}

fn assert_approx_eq_rel(
    a: f64,
    b: f64,
    rel_tol: f64,
) {

    let diff = (a - b).abs();

    let max_val = a
        .abs()
        .max(b.abs())
        .max(1e-10);

    assert!(
        diff / max_val < rel_tol,
        "Expected {}, got {}, rel diff = {}",
        b,
        a,
        diff / max_val
    );
}

// ============================================================================
// Gamma and Related Functions Tests
// ============================================================================

#[test]

fn test_gamma_numerical() {

    // Gamma(5) = 4! = 24
    assert_approx_eq(
        gamma_numerical(5.0),
        24.0,
    );

    // Gamma(0.5) = sqrt(pi)
    assert_approx_eq(
        gamma_numerical(0.5),
        std::f64::consts::PI.sqrt(),
    );

    // Gamma(1) = 1
    assert_approx_eq(
        gamma_numerical(1.0),
        1.0,
    );
}

#[test]

fn test_ln_gamma_numerical() {

    // ln(Gamma(5)) = ln(24)
    assert_approx_eq(
        ln_gamma_numerical(5.0),
        24.0f64.ln(),
    );

    // ln(Gamma(1)) = 0
    assert_approx_eq(
        ln_gamma_numerical(1.0),
        0.0,
    );
}

#[test]

fn test_digamma_numerical() {

    // ψ(1) = -γ (Euler-Mascheroni constant)
    let euler_mascheroni = 0.5772156649015329;

    assert_approx_eq(
        digamma_numerical(1.0),
        -euler_mascheroni,
    );

    // ψ(2) = 1 - γ
    assert_approx_eq(
        digamma_numerical(2.0),
        1.0 - euler_mascheroni,
    );
}

#[test]

fn test_beta_numerical() {

    // B(x, y) = Gamma(x)Gamma(y) / Gamma(x+y)
    // B(1, 1) = 1 * 1 / 1 = 1
    assert_approx_eq(
        beta_numerical(1.0, 1.0),
        1.0,
    );

    // B(2, 2) = 1 * 1 / Gamma(4) = 1 / 6
    assert_approx_eq(
        beta_numerical(2.0, 2.0),
        1.0 / 6.0,
    );

    // B(a, b) = B(b, a) (symmetry)
    assert_approx_eq(
        beta_numerical(3.0, 5.0),
        beta_numerical(5.0, 3.0),
    );
}

#[test]

fn test_ln_beta_numerical() {

    // ln(B(2, 2)) = ln(1/6)
    assert_approx_eq(
        ln_beta_numerical(2.0, 2.0),
        (1.0 / 6.0_f64).ln(),
    );
}

#[test]

fn test_regularized_incomplete_beta() {

    // I_0(a, b) = 0
    assert_approx_eq(
        regularized_incomplete_beta(2.0, 3.0, 0.0),
        0.0,
    );

    // I_1(a, b) = 1
    assert_approx_eq(
        regularized_incomplete_beta(2.0, 3.0, 1.0),
        1.0,
    );

    // I_0.5(1, 1) = 0.5 (uniform distribution)
    assert_approx_eq(
        regularized_incomplete_beta(1.0, 1.0, 0.5),
        0.5,
    );
}

// ============================================================================
// Error Functions Tests
// ============================================================================

#[test]

fn test_erf_numerical() {

    // erf(0) = 0
    assert_approx_eq(
        erf_numerical(0.0),
        0.0,
    );

    // erf(1) ~ 0.8427
    assert_approx_eq(
        erf_numerical(1.0),
        0.84270079294971,
    );

    // erf(-x) = -erf(x) (odd function)
    assert_approx_eq(
        erf_numerical(-1.0),
        -erf_numerical(1.0),
    );
}

#[test]

fn test_erfc_numerical() {

    // erfc(0) = 1
    assert_approx_eq(
        erfc_numerical(0.0),
        1.0,
    );

    // erfc(1) = 1 - erf(1)
    assert_approx_eq(
        erfc_numerical(1.0),
        1.0 - 0.84270079294971,
    );

    // erf(x) + erfc(x) = 1
    let x = 1.5;

    assert_approx_eq(
        erf_numerical(x) + erfc_numerical(x),
        1.0,
    );
}

#[test]

fn test_inverse_erf() {

    // inverse_erf(erf(x)) = x
    let x = 0.5;

    let y = erf_numerical(x);

    assert_approx_eq(inverse_erf(y), x);

    // inverse_erf(0) = 0
    assert_approx_eq(
        inverse_erf(0.0),
        0.0,
    );
}

#[test]

fn test_inverse_erfc() {

    // inverse_erfc(erfc(x)) = x
    let x = 0.5;

    let y = erfc_numerical(x);

    assert_approx_eq(inverse_erfc(y), x);

    // inverse_erfc(1) = 0
    assert_approx_eq(
        inverse_erfc(1.0),
        0.0,
    );
}

// ============================================================================
// Combinatorial Functions Tests
// ============================================================================

#[test]

fn test_factorial() {

    assert_eq!(factorial(0), 1);

    assert_eq!(factorial(1), 1);

    assert_eq!(factorial(5), 120);

    assert_eq!(
        factorial(10),
        3628800
    );

    assert_eq!(
        factorial(20),
        2432902008176640000
    );
}

#[test]

fn test_double_factorial() {

    assert_eq!(
        double_factorial(0),
        1
    );

    assert_eq!(
        double_factorial(1),
        1
    );

    assert_eq!(
        double_factorial(5),
        15
    ); // 5 * 3 * 1
    assert_eq!(
        double_factorial(6),
        48
    ); // 6 * 4 * 2
    assert_eq!(
        double_factorial(7),
        105
    ); // 7 * 5 * 3 * 1
    assert_eq!(
        double_factorial(8),
        384
    ); // 8 * 6 * 4 * 2
}

#[test]

fn test_binomial() {

    assert_eq!(binomial(5, 0), 1);

    assert_eq!(binomial(5, 1), 5);

    assert_eq!(binomial(5, 2), 10);

    assert_eq!(binomial(5, 3), 10);

    assert_eq!(binomial(5, 4), 5);

    assert_eq!(binomial(5, 5), 1);

    assert_eq!(binomial(10, 3), 120);

    assert_eq!(
        binomial(20, 10),
        184756
    );

    // C(n, k) = C(n, n-k) (symmetry)
    assert_eq!(
        binomial(15, 6),
        binomial(15, 9)
    );

    // k > n returns 0
    assert_eq!(binomial(5, 6), 0);
}

#[test]

fn test_rising_factorial() {

    // (1)_n = n!
    assert_approx_eq(
        rising_factorial(1.0, 5),
        120.0,
    );

    // (3)_4 = 3 * 4 * 5 * 6 = 360
    assert_approx_eq(
        rising_factorial(3.0, 4),
        360.0,
    );

    // (x)_0 = 1
    assert_approx_eq(
        rising_factorial(5.0, 0),
        1.0,
    );
}

#[test]

fn test_falling_factorial() {

    // (5)_{(3)} = 5 * 4 * 3 = 60
    assert_approx_eq(
        falling_factorial(5.0, 3),
        60.0,
    );

    // (10)_{(4)} = 10 * 9 * 8 * 7 = 5040
    assert_approx_eq(
        falling_factorial(10.0, 4),
        5040.0,
    );

    // (x)_{(0)} = 1
    assert_approx_eq(
        falling_factorial(5.0, 0),
        1.0,
    );

    // (n)_{(n)} = n!
    assert_approx_eq(
        falling_factorial(5.0, 5),
        120.0,
    );
}

#[test]

fn test_ln_factorial() {

    assert_approx_eq(ln_factorial(0), 0.0);

    assert_approx_eq(ln_factorial(1), 0.0);

    assert_approx_eq(
        ln_factorial(5),
        120.0_f64.ln(),
    );

    assert_approx_eq(
        ln_factorial(10),
        3628800.0_f64.ln(),
    );
}

// ============================================================================
// Bessel Functions Tests
// ============================================================================

#[test]

fn test_bessel_j0() {

    // J_0(0) = 1
    assert_approx_eq(bessel_j0(0.0), 1.0);

    // J_0 has a zero near 2.405
    assert!(bessel_j0(2.405).abs() < 0.01);
}

#[test]

fn test_bessel_j1() {

    // J_1(0) = 0
    assert_approx_eq(bessel_j1(0.0), 0.0);

    // J_1 has a maximum near 1.84
    assert!(bessel_j1(1.84) > 0.5);
}

#[test]

fn test_bessel_y0() {

    // Y_0 is negative for small positive x
    assert!(bessel_y0(0.1) < 0.0);

    // Y_0 has a zero near 0.894
    assert!(bessel_y0(0.894).abs() < 0.01);
}

#[test]

fn test_bessel_y1() {

    // Y_1 is negative for small positive x
    assert!(bessel_y1(0.1) < 0.0);
}

#[test]

fn test_bessel_i0() {

    // I_0(0) = 1
    assert_approx_eq(bessel_i0(0.0), 1.0);

    // I_0(x) > 0 for all x
    assert!(bessel_i0(1.0) > 0.0);

    assert!(bessel_i0(-1.0) > 0.0);
}

#[test]

fn test_bessel_i1() {

    // I_1(0) = 0
    assert_approx_eq(bessel_i1(0.0), 0.0);

    // I_1(x) is odd: I_1(-x) = -I_1(x)
    assert_approx_eq(
        bessel_i1(-1.0),
        -bessel_i1(1.0),
    );
}

#[test]

fn test_bessel_k0() {

    // K_0(x) > 0 for x > 0
    assert!(bessel_k0(1.0) > 0.0);

    // K_0 decreases as x increases
    assert!(bessel_k0(1.0) > bessel_k0(2.0));
}

#[test]

fn test_bessel_k1() {

    // K_1(x) > 0 for x > 0
    assert!(bessel_k1(1.0) > 0.0);

    // K_1 decreases as x increases
    assert!(bessel_k1(1.0) > bessel_k1(2.0));
}

// ============================================================================
// Other Special Functions Tests
// ============================================================================

#[test]

fn test_sinc() {

    // sinc(0) = 1
    assert_approx_eq(sinc(0.0), 1.0);

    // sinc(n) = 0 for non-zero integer n
    assert_approx_eq(sinc(1.0), 0.0);

    assert_approx_eq(sinc(-1.0), 0.0);

    assert_approx_eq(sinc(2.0), 0.0);

    // sinc(0.5) = sin(pi/2) / (pi/2) = 2/pi
    assert_approx_eq(
        sinc(0.5),
        2.0 / std::f64::consts::PI,
    );
}

#[test]

fn test_zeta() {

    // ζ(2) = π²/6
    let expected = std::f64::consts::PI.powi(2) / 6.0;

    assert_approx_eq_rel(
        zeta(2.0),
        expected,
        1e-4,
    ); // Numerical approximation tolerance
       // ζ(4) = π⁴/90
    let expected4 = std::f64::consts::PI.powi(4) / 90.0;

    assert_approx_eq_rel(
        zeta(4.0),
        expected4,
        1e-4,
    );
}

#[test]

fn test_regularized_gamma_p() {

    // P(1, x) = 1 - e^(-x)
    let x = 2.0;

    assert_approx_eq(
        regularized_gamma_p(1.0, x),
        1.0 - (-x).exp(),
    );

    // P(a, x) for small x should be small
    let small_x = 0.1;

    let result = regularized_gamma_p(2.0, small_x);

    assert!(result > 0.0 && result < 0.1); // Should be small but positive
}

#[test]

fn test_regularized_gamma_q() {

    // Q(a, x) + P(a, x) = 1
    let a = 2.0;

    let x = 1.5;

    assert_approx_eq(
        regularized_gamma_p(a, x) + regularized_gamma_q(a, x),
        1.0,
    );

    // Q(1, x) = e^(-x)
    assert_approx_eq(
        regularized_gamma_q(1.0, x),
        (-x).exp(),
    );
}

// ============================================================================
// Integration Tests
// ============================================================================

#[test]

fn test_gamma_beta_relationship() {

    // B(a, b) = Γ(a)Γ(b) / Γ(a+b)
    let a = 3.0;

    let b = 4.0;

    let beta_direct = beta_numerical(a, b);

    let beta_from_gamma = gamma_numerical(a) * gamma_numerical(b) / gamma_numerical(a + b);

    assert_approx_eq(
        beta_direct,
        beta_from_gamma,
    );
}

#[test]

fn test_factorial_gamma_relationship() {

    // Γ(n+1) = n! for non-negative integers
    for n in 0..=10 {

        let gamma_val = gamma_numerical((n + 1) as f64);

        let fact_val = factorial(n as u64) as f64;

        assert_approx_eq(gamma_val, fact_val);
    }
}

#[test]

fn test_binomial_factorial_relationship() {

    // C(n, k) = n! / (k! * (n-k)!)
    let n = 10u64;

    let k = 4u64;

    let binomial_val = binomial(n, k);

    let factorial_val = factorial(n) / (factorial(k) * factorial(n - k));

    assert_eq!(
        binomial_val,
        factorial_val
    );
}

#[test]

fn test_erf_erfc_relationship() {

    // erf(x) + erfc(x) = 1 for various x
    for x in [
        0.0, 0.5, 1.0, 1.5, 2.0, -0.5, -1.0,
    ] {

        assert_approx_eq(
            erf_numerical(x) + erfc_numerical(x),
            1.0,
        );
    }
}

#[test]

fn test_bessel_recurrence() {

    // J_{n+1}(x) = (2n/x) * J_n(x) - J_{n-1}(x)
    // For n=1: J_2(x) = (2/x) * J_1(x) - J_0(x)
    let x = 2.0;

    let j0 = bessel_j0(x);

    let j1 = bessel_j1(x);

    let j2_from_recurrence = (2.0 / x) * j1 - j0;

    // J_2(2) ≈ 0.3528
    assert!((j2_from_recurrence - 0.3528).abs() < 0.01);
}
