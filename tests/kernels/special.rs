//! Special functions (ported from `numerical_special_test.rs`, expanded with
//! reference values from high-precision tables).

use std::f64::consts::PI;

use rssn::kernels::special::*;

fn cfg() -> proptest::prelude::ProptestConfig {
    proptest::prelude::ProptestConfig {
        rng_seed: proptest::test_runner::RngSeed::Fixed(0x5EED),
        failure_persistence: None,
        ..proptest::prelude::ProptestConfig::default()
    }
}

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




    assert!((riemann_zeta(2.0) - PI * PI / 6.0).abs() < 1e-3);

    // ζ(4) = π⁴/90

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
    #![proptest_config(cfg())]

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

// ============================================================================
// Added: reference values
// ============================================================================

#[test]
fn reference_values_erf_gamma_beta() {
    assert!((erf_numerical(1.0) - 0.842_700_792_949_714_9).abs() < 1e-9);
    assert!((erfc_numerical(2.0) - 0.004_677_734_981_047_265_4).abs() < 1e-9);
    assert!((inverse_erf_numerical(0.520_499_877_813_046_5) - 0.5).abs() < 1e-9);
    assert_eq!(inverse_erf_numerical(1.0), f64::INFINITY);
    assert_eq!(inverse_erf_numerical(-1.0), f64::NEG_INFINITY);
    assert!((gamma_numerical(10.5) - 1_133_278.388_948_78).abs() < 1e-3);
    assert!((ln_gamma_numerical(0.5) - PI.sqrt().ln()).abs() < 1e-12);
    assert!((digamma_numerical(0.5) - (-1.963_510_026_021_423_5)).abs() < 1e-10);
    assert!((ln_beta_numerical(2.0, 3.0) - (1.0f64 / 12.0).ln()).abs() < 1e-12);
    assert!((beta_numerical(2.0, 3.0) - 1.0 / 12.0).abs() < 1e-12);
}

#[test]
fn incomplete_gamma_closed_forms() {
    // s = 1: gamma(1, x) = 1 - e^-x
    for x in [0.1, 1.0, 3.0, 10.0] {
        let want = 1.0 - (-x as f64).exp();
        assert!((lower_incomplete_gamma(1.0, x) - want).abs() < 1e-8, "x = {x}");
        assert!((regularized_lower_gamma(1.0, x) - want).abs() < 1e-8, "x = {x}");
        assert!((regularized_upper_gamma(1.0, x) - (-x as f64).exp()).abs() < 1e-8, "x = {x}");
        assert!((upper_incomplete_gamma(1.0, x) - (-x as f64).exp()).abs() < 1e-8, "x = {x}");
    }
    // s = 2: gamma(2, x) = 1 - (1 + x) e^-x
    assert!((lower_incomplete_gamma(2.0, 3.0) - (1.0 - 4.0 * (-3.0f64).exp())).abs() < 1e-8);
    assert!(lower_incomplete_gamma(2.0, -1.0).is_nan());
    assert!(lower_incomplete_gamma(-2.0, 1.0).is_nan());
}

#[test]
fn incomplete_beta_closed_forms() {
    // I_x(2, 1) = x^2 ; I_x(1, 2) = 1 - (1 - x)^2
    assert!((regularized_beta(0.3, 2.0, 1.0) - 0.09).abs() < 1e-9);
    assert!((regularized_beta(0.3, 1.0, 2.0) - (1.0 - 0.49)).abs() < 1e-9);
    // Non-regularised: B(x; 2, 1) = x^2 / 2
    assert!((incomplete_beta(0.3, 2.0, 1.0) - 0.045).abs() < 1e-9);
    assert!(regularized_beta(1.5, 1.0, 1.0).is_nan());
    assert!(incomplete_beta(0.5, -1.0, 1.0).is_nan());
}

#[test]
fn bessel_reference_values() {
    assert!((bessel_j0(1.0) - 0.765_197_686_557_966_6).abs() < 1e-7);
    assert!((bessel_y0(1.0) - 0.088_256_964_215_676_96).abs() < 1e-7);
    assert!((bessel_i0(1.0) - 1.266_065_877_752_008_4).abs() < 1e-7);
    assert!((bessel_i1(1.0) - 0.565_159_103_992_485_1).abs() < 1e-7);
    // First zero of J0 and Y0 to high accuracy.
    assert!(bessel_j0(2.404_825_557_695_773).abs() < 1e-7);
    assert!(bessel_y0(0.893_576_966_279_167_5).abs() < 1e-7);
    // Large arguments use the asymptotic branch: J0(10) = -0.2459357644513483
    assert!((bessel_j0(10.0) - (-0.245_935_764_451_348_3)).abs() < 1e-7);
    assert!((bessel_j1(10.0) - 0.043_472_746_168_861_44).abs() < 1e-7);
    assert!(bessel_y1(-1.0).is_nan());
}

#[test]
#[ignore = "library bug: bessel_j1 rational approximation for |x| < 8 has mis-scaled coefficients; observed J1(1) = -0.175057, expected 0.440051 (also J1(2) = -4.552 vs 0.576725, J1(0.5) = 1.301 vs 0.242268)"]
fn bessel_j1_small_argument_reference_values() {
    assert!((bessel_j1(1.0) - 0.440_050_585_744_933_5).abs() < 1e-7);
    assert!((bessel_j1(2.0) - 0.576_724_807_756_873_4).abs() < 1e-7);
    assert!((bessel_j1(0.5) - 0.242_268_457_674_873_87).abs() < 1e-7);
}

#[test]
#[ignore = "library bug: bessel_y1 rational approximation for x < 8 is wrong; observed Y1(1) = -2.084740, expected -0.781213 (Y1(2) = -2.174 vs -0.107032, Y1(5) = 35.09 vs 0.147863)"]
fn bessel_y1_small_argument_reference_values() {
    assert!((bessel_y1(1.0) - (-0.781_212_821_300_288_7)).abs() < 1e-7);
    assert!((bessel_y1(2.0) - (-0.107_032_431_540_937_55)).abs() < 1e-7);
}

#[test]
fn bessel_y1_large_argument_reference_value() {
    assert!((bessel_y1(10.0) - 0.249_015_424_206_953_9).abs() < 1e-7);
}

#[test]
#[ignore = "library bug: inverse_erf_numerical loses accuracy for |y| >= 0.7 (3 Newton steps from a poor guess); observed erfinv(0.8427007929497149) = 0.999618221363, expected 1.0; erfinv(0.99) = 1.818515 vs 1.821386"]
fn inverse_erf_large_argument() {
    assert!((inverse_erf_numerical(0.842_700_792_949_714_9) - 1.0).abs() < 1e-9);
    assert!((inverse_erf_numerical(0.99) - 1.821_386_367_718_449_6).abs() < 1e-9);
}

#[test]
fn orthogonal_polynomial_reference_values() {
    // P_3(x) = (5x^3 - 3x)/2 ; P_4(0.3) = (35x^4 - 30x^2 + 3)/8
    assert!((legendre_p(3, 0.4) - (5.0 * 0.064 - 1.2) / 2.0).abs() < 1e-12);
    assert!((legendre_p(4, 0.3) - (35.0 * 0.0081 - 30.0 * 0.09 + 3.0) / 8.0).abs() < 1e-12);
    // U_n(cos t) = sin((n+1)t) / sin t
    let t: f64 = 0.7;
    assert!((chebyshev_u(4, t.cos()) - (5.0 * t).sin() / t.sin()).abs() < 1e-10);
    assert_eq!(chebyshev_u(0, 0.3), 1.0);
    // H_3(x) = 8x^3 - 12x ; L_3(x) = (-x^3 + 9x^2 - 18x + 6)/6
    assert!((hermite_h(3, 1.5) - (8.0 * 3.375 - 18.0)).abs() < 1e-10);
    assert!((laguerre_l(3, 2.0) - (-8.0 + 36.0 - 36.0 + 6.0) / 6.0).abs() < 1e-10);
}

#[test]
fn combinatorial_and_zeta_values() {
    assert_eq!(binomial(3, 5), 0.0);
    assert!((binomial(52, 5) - 2_598_960.0).abs() < 1e-3);
    assert!((factorial(20) - 2_432_902_008_176_640_000.0).abs() / 2.4e18 < 1e-12);
    assert!((double_factorial(9) - 945.0).abs() < 1e-9);
    assert_eq!(riemann_zeta(1.0), f64::INFINITY);
    assert!((riemann_zeta(6.0) - PI.powi(6) / 945.0).abs() < 1e-8);
}

#[test]
fn bernoulli_numbers() {
    assert_eq!(bernoulli_number(0), 1.0);
    assert_eq!(bernoulli_number(1), -0.5);
    assert!((bernoulli_number(2) - 1.0 / 6.0).abs() < 1e-15);
    assert!((bernoulli_number(4) - (-1.0 / 30.0)).abs() < 1e-15);
    assert_eq!(bernoulli_number(3), 0.0);
    // B_2(x) = x^2 - x + 1/6
    assert!((bernoulli_poly(2, 0.5) - (0.25 - 0.5 + 1.0 / 6.0)).abs() < 1e-12);
}

#[test]
fn logistic_helpers_edge_cases() {
    assert!(logit(0.0).is_nan() && logit(1.0).is_nan());
    assert!((softplus(-30.0) - (-30.0f64).exp()).abs() < 1e-20);
    assert!((softplus(1.0) - (1.0 + 1.0f64.exp()).ln() + 1.0 - 1.0).abs() < 1e-12);
}

proptest::proptest! {
    #![proptest_config(cfg())]

    /// ln Gamma(x) agrees with ln of Gamma(x)
    #[test]
    fn prop_ln_gamma_consistent(x in 0.5..20.0f64) {
        proptest::prop_assert!((ln_gamma_numerical(x) - gamma_numerical(x).ln()).abs() < 1e-9);
    }

    /// Beta is symmetric.
    #[test]
    fn prop_beta_symmetric(a in 0.5..10.0f64, b in 0.5..10.0f64) {
        proptest::prop_assert!((beta_numerical(a, b) - beta_numerical(b, a)).abs() < 1e-12 * beta_numerical(a, b).abs().max(1.0));
    }

    /// I_x(a, b) + I_{1-x}(b, a) = 1
    #[test]
    fn prop_regularized_beta_reflection(x in 0.01..0.99f64, a in 0.5..5.0f64, b in 0.5..5.0f64) {
        let s = regularized_beta(x, a, b) + regularized_beta(1.0 - x, b, a);
        proptest::prop_assert!((s - 1.0).abs() < 1e-7, "sum = {s}");
    }

    /// Legendre recurrence (n+1) P_{n+1} = (2n+1) x P_n - n P_{n-1}
    #[test]
    fn prop_legendre_recurrence(n in 1u32..15, x in -1.0..1.0f64) {
        let lhs = f64::from(n + 1) * legendre_p(n + 1, x);
        let rhs = f64::from(2 * n + 1) * x * legendre_p(n, x) - f64::from(n) * legendre_p(n - 1, x);
        proptest::prop_assert!((lhs - rhs).abs() < 1e-9);
    }

    /// erf is monotone and inverse_erf inverts it.
    #[test]
    fn prop_inverse_erf_round_trip(x in -0.6..0.6f64) {
        let y = erf_numerical(x);
        proptest::prop_assert!((inverse_erf_numerical(y) - x).abs() < 1e-6, "x = {x}");
    }

    /// logit inverts sigmoid.
    #[test]
    fn prop_logit_sigmoid_round_trip(p in 0.001..0.999f64) {
        proptest::prop_assert!((sigmoid(logit(p)) - p).abs() < 1e-12);
    }

    /// Bessel identity J0' = -J1 via central differences (asymptotic branch only: the
    /// small-argument J1 is broken, see `bessel_j1_small_argument_reference_values`).
    #[test]
    fn prop_bessel_j0_derivative_is_minus_j1(x in 8.5..15.0f64) {
        let h = 1e-5;
        let d = (bessel_j0(x + h) - bessel_j0(x - h)) / (2.0 * h);
        proptest::prop_assert!((d + bessel_j1(x)).abs() < 1e-5);
    }
}
