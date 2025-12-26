use std::f64::consts::PI;

use rssn::numerical::special::*;

// ============================================================================
// Gamma Function Tests
// ============================================================================

#[test]

fn test_gamma() {

    // Γ(1) = 0! = 1
    assert!((gamma_numerical(1.0) - 1.0).abs() < 1e-10);

    // Γ(2) = 1! = 1
    assert!((gamma_numerical(2.0) - 1.0).abs() < 1e-10);

    // Γ(3) = 2! = 2
    assert!((gamma_numerical(3.0) - 2.0).abs() < 1e-10);

    // Γ(5) = 4! = 24
    assert!((gamma_numerical(5.0) - 24.0).abs() < 1e-10);

    // Γ(0.5) = √π
    assert!((gamma_numerical(0.5) - PI.sqrt()).abs() < 1e-10);
}

#[test]

fn test_ln_gamma() {

    // ln(Γ(5)) = ln(24) ≈ 3.178
    assert!((ln_gamma_numerical(5.0) - 24_f64.ln()).abs() < 1e-10);
}

#[test]

fn test_digamma() {

    // ψ(1) = -γ where γ is Euler's constant ≈ -0.5772
    let psi1 = digamma_numerical(1.0);

    assert!((psi1 - (-0.5772156649015329)).abs() < 1e-10);
}

#[test]

fn test_lower_incomplete_gamma() {

    // At x=0, γ(s,0) = 0
    assert!((lower_incomplete_gamma(2.0, 0.0) - 0.0).abs() < 1e-10);

    // γ(s, ∞) → Γ(s)
    let lic = lower_incomplete_gamma(2.0, 100.0);

    let full_gamma = gamma_numerical(2.0);

    assert!((lic - full_gamma).abs() < 1e-5);
}

// ============================================================================
// Beta Function Tests
// ============================================================================

#[test]

fn test_beta() {

    // B(1,1) = 1
    assert!((beta_numerical(1.0, 1.0) - 1.0).abs() < 1e-10);

    // B(a,b) = Γ(a)Γ(b)/Γ(a+b)
    let b23 = beta_numerical(2.0, 3.0);

    let expected = gamma_numerical(2.0) * gamma_numerical(3.0) / gamma_numerical(5.0);

    assert!((b23 - expected).abs() < 1e-10);
}

#[test]

fn test_regularized_beta() {

    // I_0(a,b) = 0
    assert!((regularized_beta(0.0, 2.0, 3.0) - 0.0).abs() < 1e-10);

    // I_1(a,b) = 1
    assert!((regularized_beta(1.0, 2.0, 3.0) - 1.0).abs() < 1e-10);

    // I_0.5(1,1) = 0.5 (uniform distribution)
    assert!((regularized_beta(0.5, 1.0, 1.0) - 0.5).abs() < 1e-10);
}

// ============================================================================
// Error Function Tests
// ============================================================================

#[test]

fn test_erf() {

    // erf(0) = 0
    assert!(erf_numerical(0.0).abs() < 1e-10);

    // erf(∞) → 1
    assert!((erf_numerical(10.0) - 1.0).abs() < 1e-10);

    // erf(-x) = -erf(x) (odd function)
    assert!((erf_numerical(-1.0) + erf_numerical(1.0)).abs() < 1e-10);
}

#[test]

fn test_erfc() {

    // erfc(x) = 1 - erf(x)
    let x = 1.5;

    assert!((erfc_numerical(x) - (1.0 - erf_numerical(x))).abs() < 1e-10);
}

#[test]

fn test_inverse_erf() {

    // erf⁻¹(erf(x)) = x
    let x = 0.5;

    let y = erf_numerical(x);

    assert!((inverse_erf_numerical(y) - x).abs() < 1e-10);
}

// ============================================================================
// Bessel Function Tests
// ============================================================================

#[test]

fn test_bessel_j0() {

    // J₀(0) = 1
    assert!((bessel_j0(0.0) - 1.0).abs() < 1e-10);

    // J₀ has a zero near x ≈ 2.4048
    assert!(bessel_j0(2.4048).abs() < 0.001);
}

#[test]

fn test_bessel_j1() {

    // J₁(0) = 0
    assert!(bessel_j1(0.0).abs() < 1e-10);
    // J₁ has an extremum at x ≈ 1.841
}

#[test]

fn test_bessel_y0() {

    // Y₀(x) is undefined for x < 0
    assert!(bessel_y0(-1.0).is_nan());

    // Y₀ has a zero near x ≈ 0.8936
    assert!(bessel_y0(0.8936).abs() < 0.01);
}

#[test]

fn test_bessel_i0() {

    // I₀(0) = 1
    assert!((bessel_i0(0.0) - 1.0).abs() < 1e-10);

    // I₀(x) > 1 for x > 0
    assert!(bessel_i0(1.0) > 1.0);
}

// ============================================================================
// Orthogonal Polynomial Tests
// ============================================================================

#[test]

fn test_legendre_p() {

    // P₀(x) = 1
    assert!((legendre_p(0, 0.5) - 1.0).abs() < 1e-10);

    // P₁(x) = x
    assert!((legendre_p(1, 0.5) - 0.5).abs() < 1e-10);

    // P₂(x) = (3x² - 1)/2
    let x = 0.5;

    let expected = (3.0 * x * x - 1.0) / 2.0;

    assert!((legendre_p(2, x) - expected).abs() < 1e-10);
}

#[test]

fn test_chebyshev_t() {

    // T₀(x) = 1
    assert!((chebyshev_t(0, 0.5) - 1.0).abs() < 1e-10);

    // T₁(x) = x
    assert!((chebyshev_t(1, 0.5) - 0.5).abs() < 1e-10);

    // T₂(x) = 2x² - 1
    let x = 0.5;

    assert!((chebyshev_t(2, x) - (2.0 * x * x - 1.0)).abs() < 1e-10);

    // Tₙ(cos(θ)) = cos(nθ)
    let theta = PI / 4.0;

    assert!((chebyshev_t(3, theta.cos()) - (3.0 * theta).cos()).abs() < 1e-10);
}

#[test]

fn test_hermite_h() {

    // H₀(x) = 1
    assert!((hermite_h(0, 1.0) - 1.0).abs() < 1e-10);

    // H₁(x) = 2x
    assert!((hermite_h(1, 1.0) - 2.0).abs() < 1e-10);

    // H₂(x) = 4x² - 2
    let x = 2.0;

    assert!((hermite_h(2, x) - (4.0 * x * x - 2.0)).abs() < 1e-10);
}

#[test]

fn test_laguerre_l() {

    // L₀(x) = 1
    assert!((laguerre_l(0, 1.0) - 1.0).abs() < 1e-10);

    // L₁(x) = 1 - x
    assert!((laguerre_l(1, 1.0) - 0.0).abs() < 1e-10);

    // L₂(x) = (x² - 4x + 2)/2
    let x = 2.0;

    let expected = (x * x - 4.0 * x + 2.0) / 2.0;

    assert!((laguerre_l(2, x) - expected).abs() < 1e-10);
}

// ============================================================================
// Other Special Function Tests
// ============================================================================

#[test]

fn test_factorial() {

    assert!((factorial(0) - 1.0).abs() < 1e-10);

    assert!((factorial(1) - 1.0).abs() < 1e-10);

    assert!((factorial(5) - 120.0).abs() < 1e-10);

    assert!((factorial(10) - 3628800.0).abs() < 1e-5);
}

#[test]

fn test_double_factorial() {

    // 0!! = 1, 1!! = 1
    assert!((double_factorial(0) - 1.0).abs() < 1e-10);

    assert!((double_factorial(1) - 1.0).abs() < 1e-10);

    // 5!! = 5*3*1 = 15
    assert!((double_factorial(5) - 15.0).abs() < 1e-10);

    // 6!! = 6*4*2 = 48
    assert!((double_factorial(6) - 48.0).abs() < 1e-10);
}

#[test]

fn test_binomial() {

    assert!((binomial(5, 0) - 1.0).abs() < 1e-10);

    assert!((binomial(5, 5) - 1.0).abs() < 1e-10);

    assert!((binomial(5, 2) - 10.0).abs() < 1e-10);

    assert!((binomial(10, 5) - 252.0).abs() < 1e-10);
}

#[test]

fn test_riemann_zeta() {

    // ζ(2) = π²/6
    println!(
        "riemann_zeta(2.0) = {}",
        riemann_zeta(2.0)
    );

    println!(
        "PI * PI / 6.0 = {}",
        PI * PI / 6.0
    );

    println!(
        "diff = {}",
        (riemann_zeta(2.0) - PI * PI / 6.0).abs()
    );

    println!("1e-3 = {}", 1e-3);

    assert!((riemann_zeta(2.0) - PI * PI / 6.0).abs() < 1e-3);

    // ζ(4) = π⁴/90
    println!(
        "riemann_zeta(4.0) = {}",
        riemann_zeta(4.0)
    );

    assert!((riemann_zeta(4.0) - PI.powi(4) / 90.0).abs() < 1e-5);
}

#[test]

fn test_sinc() {

    // sinc(0) = 1
    assert!((sinc(0.0) - 1.0).abs() < 1e-10);

    // sinc(1) = sin(π)/π = 0
    assert!(sinc(1.0).abs() < 1e-10);

    // sinc(0.5) = sin(π/2)/(π/2) = 2/π
    assert!((sinc(0.5) - 2.0 / PI).abs() < 1e-10);
}

#[test]

fn test_sigmoid() {

    // σ(0) = 0.5
    assert!((sigmoid(0.0) - 0.5).abs() < 1e-10);

    // σ(x) → 1 as x → ∞
    assert!((sigmoid(100.0) - 1.0).abs() < 1e-10);

    // σ(x) → 0 as x → -∞
    assert!(sigmoid(-100.0).abs() < 1e-10);

    // σ(-x) = 1 - σ(x)
    assert!((sigmoid(-2.0) - (1.0 - sigmoid(2.0))).abs() < 1e-10);
}

#[test]

fn test_softplus() {

    // softplus(0) = ln(2)
    assert!((softplus(0.0) - 2_f64.ln()).abs() < 1e-10);

    // softplus(x) ≈ x for large x
    assert!((softplus(100.0) - 100.0).abs() < 1e-10);
}

#[test]

fn test_logit() {

    // logit(0.5) = 0
    assert!(logit(0.5).abs() < 1e-10);

    // logit is inverse of sigmoid
    let p = 0.7;

    assert!((sigmoid(logit(p)) - p).abs() < 1e-10);
}

// ============================================================================
// Property-Based Tests
// ============================================================================

proptest::proptest! {
    /// Γ(x+1) = x * Γ(x) (recurrence relation)
    #[test]
    fn prop_gamma_recurrence(x in 1.0..10.0f64) {
        let lhs = gamma_numerical(x + 1.0);
        let rhs = x * gamma_numerical(x);
        proptest::prop_assert!((lhs - rhs).abs() < 1e-8 * lhs.abs());
    }

    /// erf + erfc = 1
    #[test]
    fn prop_erf_erfc_sum(x in -5.0..5.0f64) {
        let sum = erf_numerical(x) + erfc_numerical(x);
        proptest::prop_assert!((sum - 1.0).abs() < 1e-10);
    }

    /// sigmoid is bounded in (0, 1)
    #[test]
    fn prop_sigmoid_bounded(x in -100.0..100.0f64) {
        let s = sigmoid(x);
        println!("sigmoid({}) = {}", x, s);
        proptest::prop_assert!(s > 0.0 && s <= 1.0);
    }

    /// Legendre P_n(1) = 1 for all n
    #[test]
    fn prop_legendre_at_one(n in 0u32..20) {
        proptest::prop_assert!((legendre_p(n, 1.0) - 1.0).abs() < 1e-10);
    }

    /// Chebyshev T_n is bounded: |T_n(x)| ≤ 1 for |x| ≤ 1
    #[test]
    fn prop_chebyshev_bounded(n in 0u32..20, x in -1.0..1.0f64) {
        proptest::prop_assert!(chebyshev_t(n, x).abs() <= 1.0 + 1e-10);
    }

    /// Binomial symmetry: C(n,k) = C(n, n-k)
    #[test]
    fn prop_binomial_symmetry(n in 1u64..20, k in 0u64..20) {
        let k = k % (n + 1); // Ensure k <= n
        let lhs = binomial(n, k);
        let rhs = binomial(n, n - k);
        proptest::prop_assert!((lhs - rhs).abs() < 1e-10);
    }

    /// sinc is even: sinc(-x) = sinc(x)
    #[test]
    fn prop_sinc_even(x in 0.01..10.0f64) {
        proptest::prop_assert!((sinc(x) - sinc(-x)).abs() < 1e-10);
    }

    /// J₀(-x) = J₀(x) (even function)
    #[test]
    fn prop_bessel_j0_even(x in 0.0..10.0f64) {
        proptest::prop_assert!((bessel_j0(x) - bessel_j0(-x)).abs() < 1e-10);
    }

    /// J₁(-x) = -J₁(x) (odd function)
    #[test]
    fn prop_bessel_j1_odd(x in 0.01..10.0f64) {
        proptest::prop_assert!((bessel_j1(x) + bessel_j1(-x)).abs() < 1e-10);
    }
}
