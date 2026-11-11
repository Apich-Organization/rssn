//! Real special functions at the edges: overflow thresholds, poles, tails,
//! argument reduction and the identities that tie the families together.

use std::f64::consts::PI;

use rssn::kernels::special::*;
use rssn::kernels::special_ext::*;

fn near(
    got: f64,
    want: f64,
    rel: f64,
) -> bool {
    if got.is_nan() || want.is_nan() {
        return got.is_nan() && want.is_nan();
    }
    got == want || (got - want).abs() <= rel * want.abs()
}

/// Collects failures so that one run reports every wrong value.
struct Check(Vec<String>);
impl Check {
    fn val(
        &mut self,
        name: &str,
        got: f64,
        want: f64,
        rel: f64,
    ) {
        if !near(got, want, rel) {
            self.0.push(format!("{name}: got {got:e}, want {want:e}"));
        }
    }

    fn finish(self) {
        assert!(
            self.0.is_empty(),
            "{} wrong values:\n{}",
            self.0.len(),
            self.0.join("\n")
        );
    }
}

#[test]
fn gamma_family_edges() {
    let mut c = Check(Vec::new());
    c.val(
        "gamma(-0.5)",
        gamma_numerical(-0.5),
        -3.544_907_701_811_031_8,
        1e-13,
    );
    c.val("gamma(0.5)", gamma_numerical(0.5), PI.sqrt(), 1e-14);
    c.val("gamma(1e-300)", gamma_numerical(1e-300), 1e300, 1e-13);
    // 170! is the largest factorial in f64: Gamma(171) must not overflow.
    c.val(
        "gamma(171)",
        gamma_numerical(171.0),
        7.257_415_615_307_994e306,
        1e-13,
    );
    c.val(
        "factorial(170)",
        factorial(170),
        7.257_415_615_307_994e306,
        1e-13,
    );
    c.val("gamma(172)", gamma_numerical(172.0), f64::INFINITY, 0.0);
    c.val("factorial(171)", factorial(171), f64::INFINITY, 0.0);
    c.val(
        "gamma(171.5)/gamma(170.5)",
        gamma_numerical(171.5) / gamma_numerical(170.5),
        170.5,
        1e-12,
    );
    c.val(
        "gamma(150.5)",
        gamma_numerical(150.5),
        (ln_gamma_numerical(150.5)).exp(),
        1e-11,
    );
    c.val("gamma(-1)", gamma_numerical(-1.0), f64::NAN, 0.0);
    c.val("gamma(-2)", gamma_numerical(-2.0), f64::NAN, 0.0);
    c.val("gamma(0)", gamma_numerical(0.0), f64::INFINITY, 0.0);
    c.val("gamma(-0)", gamma_numerical(-0.0), f64::NEG_INFINITY, 0.0);
    c.val("gamma(nan)", gamma_numerical(f64::NAN), f64::NAN, 0.0);
    c.val(
        "gamma(-inf)",
        gamma_numerical(f64::NEG_INFINITY),
        f64::NAN,
        0.0,
    );
    c.val(
        "gamma(inf)",
        gamma_numerical(f64::INFINITY),
        f64::INFINITY,
        0.0,
    );
    c.val("lgamma(1)", ln_gamma_numerical(1.0), 0.0, 0.0);
    c.val("lgamma(2)", ln_gamma_numerical(2.0), 0.0, 0.0);
    c.val("lgamma(3)", ln_gamma_numerical(3.0), 2f64.ln(), 1e-15);
    c.val(
        "lgamma(1e10)",
        ln_gamma_numerical(1e10),
        220_258_509_288.811_08,
        1e-13,
    );
    c.val(
        "lgamma(-0.5)",
        ln_gamma_numerical(-0.5),
        1.265_512_123_484_645_4,
        1e-13,
    );
    c.val(
        "digamma(1)",
        digamma_numerical(1.0),
        -0.577_215_664_901_532_9,
        1e-14,
    );
    c.val(
        "digamma(-0.5)",
        digamma_numerical(-0.5),
        0.036_489_973_978_576_52,
        1e-12,
    );
    c.val(
        "digamma(1e-10)",
        digamma_numerical(1e-10),
        -1e10 - 0.577_215_664_901_532_9,
        1e-12,
    );
    c.val("digamma(0)", digamma_numerical(0.0), f64::NEG_INFINITY, 0.0);
    c.val("digamma(-2)", digamma_numerical(-2.0), f64::NAN, 0.0);
    c.val(
        "polygamma(1,1)",
        polygamma_numerical(1, 1.0),
        PI * PI / 6.0,
        1e-13,
    );
    // psi'''(2.5) = 6 sum_k (2.5 + k)^-4, summed far enough
    let tri: f64 = 6.0
        * (0..2_000_000)
            .map(|k| (2.5 + f64::from(k)).powi(-4))
            .sum::<f64>();
    c.val("polygamma(3,2.5)", polygamma_numerical(3, 2.5), tri, 1e-9);
    // binomial coefficients: no NaN from overflowing gammas, no panic on u64::MAX
    c.val(
        "C(1000,500)",
        binomial(1000, 500),
        2.702_882_409_454_366e299,
        1e-12,
    );
    c.val(
        "C(10^6,3)",
        binomial(1_000_000, 3),
        166_666_166_667_000_000.0,
        1e-12,
    );
    c.val("C(u64max,1)", binomial(u64::MAX, 1), u64::MAX as f64, 1e-15);
    c.val(
        "C(u64max,u64max-1)",
        binomial(u64::MAX, u64::MAX - 1),
        u64::MAX as f64,
        1e-15,
    );
    c.val("C(5,7)", binomial(5, 7), 0.0, 0.0);
    c.val("C(52,5)", binomial(52, 5), 2_598_960.0, 1e-14);
    c.finish();
}

#[test]
fn error_function_family_edges() {
    let mut c = Check(Vec::new());
    c.val(
        "erf(1e-10)",
        erf_numerical(1e-10),
        1.128_379_167_095_512_6e-10,
        1e-14,
    );
    c.val(
        "erf(0.5)",
        erf_numerical(0.5),
        0.520_499_877_813_046_5,
        1e-14,
    );
    c.val(
        "erfc(5)",
        erfc_numerical(5.0),
        1.537_459_794_428_035e-12,
        1e-12,
    );
    c.val(
        "erfc(10)",
        erfc_numerical(10.0),
        2.088_487_583_762_545e-45,
        1e-12,
    );
    c.val(
        "erfc(26)",
        erfc_numerical(26.0),
        5.663_192_408_856_143e-296,
        1e-10,
    );
    c.val(
        "erfc(-5)",
        erfc_numerical(-5.0),
        2.0 - 1.537_459_794_428_035e-12,
        1e-15,
    );
    c.val("erf(inf)", erf_numerical(f64::INFINITY), 1.0, 0.0);
    c.val("erf(nan)", erf_numerical(f64::NAN), f64::NAN, 0.0);
    c.val(
        "ierf(.5)",
        inverse_erf_numerical(0.5),
        0.476_936_276_204_469_9,
        1e-13,
    );
    c.val("ierf(1)", inverse_erf_numerical(1.0), f64::INFINITY, 0.0);
    c.val(
        "ierf(-1)",
        inverse_erf_numerical(-1.0),
        f64::NEG_INFINITY,
        0.0,
    );
    c.val("ierf(2)", inverse_erf_numerical(2.0), f64::NAN, 0.0);
    c.val(
        "erf(ierf(1-1e-15))",
        erf_numerical(inverse_erf_numerical(1.0 - 1e-15)),
        1.0 - 1e-15,
        1e-13,
    );
    c.val("erfi(1)", erfi(1.0), 1.650_425_758_797_542_8, 1e-13);
    // inverse_erfc must stay accurate far into the tail (Newton on ln erfc).
    for p in [
        1.5, 1.0, 0.5, 1e-2, 1e-5, 1e-12, 1e-30, 1e-100, 1e-200, 1e-300,
    ] {
        let y = inverse_erfc(p);
        c.val(
            &format!("erfc(inverse_erfc({p:e}))"),
            erfc_numerical(y),
            p,
            1e-9,
        );
    }
    c.val("inverse_erfc(0)", inverse_erfc(0.0), f64::INFINITY, 0.0);
    c.val("inverse_erfc(2)", inverse_erfc(2.0), f64::NEG_INFINITY, 0.0);
    c.val("inverse_erfc(3)", inverse_erfc(3.0), f64::NAN, 0.0);
    c.val("inverse_erfc(1)", inverse_erfc(1.0), 0.0, 0.0);
    c.finish();
}

#[test]
fn incomplete_gamma_and_beta_tails() {
    let mut c = Check(Vec::new());
    c.val(
        "P(.5,2)",
        regularized_lower_gamma(0.5, 2.0),
        erf_numerical(2f64.sqrt()),
        1e-13,
    );
    c.val(
        "P(100,100)",
        regularized_lower_gamma(100.0, 100.0),
        0.513_298_798_279_148_6,
        1e-10,
    );
    c.val("P(a,0)", regularized_lower_gamma(2.0, 0.0), 0.0, 0.0);
    c.val(
        "P(a,inf)",
        regularized_lower_gamma(2.0, f64::INFINITY),
        1.0,
        0.0,
    );
    c.val(
        "Q(a,inf)",
        regularized_upper_gamma(2.0, f64::INFINITY),
        0.0,
        0.0,
    );
    c.val("P(0,1)", regularized_lower_gamma(0.0, 1.0), 1.0, 0.0);
    c.val("P(-1,1)", regularized_lower_gamma(-1.0, 1.0), f64::NAN, 0.0);
    c.val("P(1,-1)", regularized_lower_gamma(1.0, -1.0), f64::NAN, 0.0);
    c.val(
        "P(nan,1)",
        regularized_lower_gamma(f64::NAN, 1.0),
        f64::NAN,
        0.0,
    );
    // Large order: the mean of Gamma(a) is a, so P(a, a) is a little above 1/2.
    let p = regularized_lower_gamma(1e5, 1e5);
    assert!((p - 0.5).abs() < 1e-2 && p > 0.5, "P(1e5,1e5) = {p}");
    // Q for integer order is a finite Poisson sum, exactly: Q(n,x) = e^-x sum_{k<n} x^k/k!
    for (n, x) in [
        (1_u32, 0.5_f64),
        (5, 3.0),
        (20, 100.0),
        (20, 5.0),
        (50, 20.0),
        (100, 200.0),
    ] {
        let mut term = (-x).exp();
        let mut sum = term;
        for k in 1..n {
            term *= x / f64::from(k);
            sum += term;
        }
        c.val(
            &format!("Q({n},{x})"),
            regularized_upper_gamma(f64::from(n), x),
            sum,
            1e-11,
        );
        // Gamma(n, x) = (n - 1)! Q(n, x)
        c.val(
            &format!("Gamma({n},{x})"),
            upper_incomplete_gamma(f64::from(n), x),
            factorial(u64::from(n) - 1) * sum,
            1e-11,
        );
        // P + Q = 1 where neither is tiny
        let (p, q) = (
            regularized_lower_gamma(f64::from(n), x),
            regularized_upper_gamma(f64::from(n), x),
        );
        if p > 1e-3 && q > 1e-3 {
            c.val(&format!("P+Q({n},{x})"), p + q, 1.0, 1e-14);
        }
    }
    // Q(1/2, x) = erfc(sqrt x): the upper tail down to 1e-300
    for x in [0.1, 1.0, 10.0, 30.0, 100.0, 300.0] {
        c.val(
            &format!("Q(.5,{x})"),
            regularized_upper_gamma(0.5, x),
            erfc_numerical(x.sqrt()),
            1e-10,
        );
    }
    // Lower incomplete gamma of a large argument is not truncated by an
    // iteration cap any more: gamma(3, 500) = 2 - tiny.
    c.val(
        "gamma(3,500)",
        lower_incomplete_gamma(3.0, 500.0),
        2.0,
        1e-14,
    );
    c.val("gamma(2,0)", lower_incomplete_gamma(2.0, 0.0), 0.0, 0.0);
    // Gamma(0, x) = E1(x); Gamma(-1, x) = e^-x/x - E1(x); Gamma(-1/2, x)
    c.val(
        "Gamma(0,1)",
        upper_incomplete_gamma(0.0, 1.0),
        exp_integral_e1(1.0),
        1e-13,
    );
    c.val(
        "Gamma(-1,2)",
        upper_incomplete_gamma(-1.0, 2.0),
        (-2f64).exp() / 2.0 - exp_integral_e1(2.0),
        1e-12,
    );
    let x = 0.5_f64;
    c.val(
        "Gamma(-1/2,.5)",
        upper_incomplete_gamma(-0.5, x),
        2.0 * (-x).exp() / x.sqrt() - 2.0 * PI.sqrt() * erfc_numerical(x.sqrt()),
        1e-12,
    );
    // Regularised beta
    c.val("Ix(.5;2,3)", regularized_beta(0.5, 2.0, 3.0), 0.6875, 1e-14);
    c.val(
        "Ix(.5;1e3,1e3)",
        regularized_beta(0.5, 1e3, 1e3),
        0.5,
        1e-12,
    );
    c.val(
        "Ix(.5;1e-3,1e-3)",
        regularized_beta(0.5, 1e-3, 1e-3),
        0.5,
        1e-12,
    );
    c.val(
        "Ix(1e-10;.5,.5)",
        regularized_beta(1e-10, 0.5, 0.5),
        2.0 / PI * 1e-5f64.asin(),
        1e-9,
    );
    c.finish();
}

#[test]
fn bessel_edges() {
    let mut c = Check(Vec::new());
    c.val(
        "j0(100)",
        bessel_j0(100.0),
        0.019_985_850_304_223_122,
        1e-11,
    );
    c.val(
        "j1(100)",
        bessel_j1(100.0),
        -0.077_145_352_014_112_16,
        1e-11,
    );
    c.val("y0(1)", bessel_y0(1.0), 0.088_256_964_215_676_97, 1e-13);
    c.val("y0(0)", bessel_y0(0.0), f64::NEG_INFINITY, 0.0);
    c.val("y0(-1)", bessel_y0(-1.0), f64::NAN, 0.0);
    c.val("y1(1e-5)", bessel_y1(1e-5), -2.0 / PI / 1e-5, 1e-6);
    c.val(
        "j(.5,2)",
        bessel_j(0.5, 2.0),
        (2.0 / (PI * 2.0)).sqrt() * 2f64.sin(),
        1e-12,
    );
    c.val(
        "j(-.5,2)",
        bessel_j(-0.5, 2.0),
        (2.0 / (PI * 2.0)).sqrt() * 2f64.cos(),
        1e-12,
    );
    c.val(
        "j(2,-3)",
        bessel_j(2.0, -3.0),
        0.486_091_260_585_891_1,
        1e-11,
    );
    c.val(
        "y(.5,2)",
        bessel_y(0.5, 2.0),
        -(2.0 / (PI * 2.0)).sqrt() * 2f64.cos(),
        1e-12,
    );
    c.val(
        "i(.5,2)",
        bessel_i(0.5, 2.0),
        (2.0 / (PI * 2.0)).sqrt() * 2f64.sinh(),
        1e-12,
    );
    c.val("i0(800)", bessel_i0(800.0), f64::INFINITY, 0.0);
    c.val("i0(-3)", bessel_i0(-3.0), bessel_i0(3.0), 0.0);
    c.val("k0(.1)", bessel_k0(0.1), 2.427_069_024_702_017, 1e-13);
    c.val("k1(2)", bessel_k1(2.0), 0.139_865_881_816_522_43, 1e-13);
    c.val("k0(0)", bessel_k0(0.0), f64::INFINITY, 0.0);
    c.val("k0(800)", bessel_k0(800.0), 0.0, 0.0);
    c.val(
        "k(.5,2)",
        bessel_k(0.5, 2.0),
        (PI / 4.0).sqrt() * (-2f64).exp(),
        1e-12,
    );
    // J_nu(x) for x << nu is exponentially small and must keep its relative
    // accuracy: the ratio J_{nu+1}/J_nu tends to x / (2 (nu + 1)), and the
    // leading term (x/2)^nu / Gamma(nu + 1) fixes the magnitude.
    for (nu, x) in [(100.0, 1.0), (50.0, 2.0), (20.5, 0.5), (30.0, 10.0)] {
        let j = bessel_j(nu, x);
        assert!(j > 0.0, "J_{nu}({x}) = {j}");
        let ratio = bessel_j(nu + 1.0, x) / j;
        let lead = x / (2.0 * (nu + 1.0));
        c.val(
            &format!("J_{}/J_{nu}({x})", nu + 1.0),
            ratio,
            lead,
            4.0 * (x / nu).powi(2) + 1e-12,
        );
        let magnitude = (nu * (x / 2.0).ln() - ln_gamma_numerical(nu + 1.0)).exp();
        let correction = x * x / (4.0 * (nu + 1.0));
        if correction < 0.1 {
            c.val(
                &format!("|J_{nu}({x})|"),
                j,
                magnitude,
                1.5 * correction + 1e-12,
            );
        }
    }
    c.finish();
}

#[test]
fn airy_lambert_zeta_polylog() {
    let mut c = Check(Vec::new());
    c.val("ai(-10)", airy_ai(-10.0), 0.040_241_238_486_443_19, 1e-11);
    c.val("ai(10)", airy_ai(10.0), 1.104_753_255_289_868_7e-10, 1e-11);
    c.val("ai(0)", airy_ai(0.0), 0.355_028_053_887_817_2, 1e-14);
    c.val(
        "ai'(1)",
        airy_ai_prime(1.0),
        -0.159_147_441_296_793_28,
        1e-11,
    );
    c.val("bi(200)", airy_bi(200.0), f64::INFINITY, 0.0);
    c.val("ai(nan)", airy_ai(f64::NAN), f64::NAN, 0.0);
    // lambert_w(x, true) is W0, false is W-1
    c.val(
        "W0(1)",
        lambert_w(1.0, true),
        0.567_143_290_409_783_8,
        1e-14,
    );
    c.val("W0(-1/e)", lambert_w(-(-1f64).exp(), true), -1.0, 1e-7);
    c.val(
        "W-1(-.3)",
        lambert_w(-0.3, false),
        -1.781_337_023_421_627_7,
        1e-13,
    );
    c.val("W0(-1)", lambert_w(-1.0, true), f64::NAN, 0.0);
    c.val(
        "W0(inf)",
        lambert_w(f64::INFINITY, true),
        f64::INFINITY,
        0.0,
    );
    c.val("W0(0)", lambert_w(0.0, true), 0.0, 0.0);
    c.val("W0(1e-300)", lambert_w(1e-300, true), 1e-300, 1e-13);
    c.val("W-1(1)", lambert_w(1.0, false), f64::NAN, 0.0);
    let w = lambert_w(1e300, true);
    c.val("W0(1e300)e^W0", w * w.exp(), 1e300, 1e-12);
    let w = lambert_w(-1e-300, false);
    c.val("W-1(-1e-300)e^W", w * w.exp(), -1e-300, 1e-12);
    c.val("zeta(2)", riemann_zeta(2.0), PI * PI / 6.0, 1e-14);
    c.val("zeta(-1)", riemann_zeta(-1.0), -1.0 / 12.0, 1e-13);
    c.val("zeta(0)", riemann_zeta(0.0), -0.5, 1e-14);
    c.val("zeta(1)", riemann_zeta(1.0), f64::INFINITY, 0.0);
    c.val("zeta(-2)", riemann_zeta(-2.0), 0.0, 0.0);
    c.val("zeta(-3)", riemann_zeta(-3.0), 1.0 / 120.0, 1e-13);
    c.val(
        "zeta(-0.5)",
        riemann_zeta(-0.5),
        -0.207_886_224_977_354_5,
        1e-13,
    );
    c.val("zeta(1e3)", riemann_zeta(1000.0), 1.0, 1e-15);
    c.val("zeta(-200)", riemann_zeta(-200.0), 0.0, 0.0);
    // Near the pole: zeta(1 + e) = 1/e + gamma - gamma_1 e + ...
    for e in [1e-10_f64, 1e-6, 1e-3, 0.05, -1e-10, -1e-6, -1e-3, -0.05] {
        let e = (1.0 + e) - 1.0; // the exact offset of the double nearest 1 + e
        let want = 1.0 / e + 0.577_215_664_901_532_9 + 0.072_815_845_483_676_72 * e
            - 0.009_690_363_192_872_318 * e * e / 2.0
            - 0.002_053_834_420_303_346 * e.powi(3) / 6.0
            + 0.002_325_370_065_467_3 * e.powi(4) / 24.0;
        c.val(
            &format!("zeta(1+{e:e})"),
            riemann_zeta(1.0 + e),
            want,
            1e-12,
        );
    }
    // Continuity across the switch between the Laurent series and Hurwitz.
    let (below, above) = (
        riemann_zeta(1.0 + 0.1 - 1e-9),
        riemann_zeta(1.0 + 0.1 + 1e-9),
    );
    assert!(
        (below - above).abs() < 1e-5 * below.abs(),
        "{below} vs {above}"
    );
    c.val("hzeta(2,1)", hurwitz_zeta(2.0, 1.0), PI * PI / 6.0, 1e-13);
    c.val("hzeta(2,.5)", hurwitz_zeta(2.0, 0.5), PI * PI / 2.0, 1e-13);
    c.val("hzeta(1,1)", hurwitz_zeta(1.0, 1.0), f64::INFINITY, 0.0);
    c.val("Li2(.5)", polylog(2.0, 0.5), 0.582_240_526_465_012_5, 1e-13);
    c.val("Li2(1)", polylog(2.0, 1.0), PI * PI / 6.0, 1e-13);
    c.val("Li2(-1)", polylog(2.0, -1.0), -PI * PI / 12.0, 1e-13);
    c.val("Li1(.5)", polylog(1.0, 0.5), 2f64.ln(), 1e-14);
    c.val("Li3(1)", polylog(3.0, 1.0), 1.202_056_903_159_594_2, 1e-13);
    c.val("Li0(.5)", polylog(0.0, 0.5), 1.0, 1e-14);
    c.val("Li-1(.5)", polylog(-1.0, 0.5), 2.0, 1e-13);
    c.finish();
}

#[test]
fn elliptic_hypergeometric_and_integrals() {
    let mut c = Check(Vec::new());
    c.val("K(.5)", elliptic_k(0.5), 1.854_074_677_301_371_9, 1e-14);
    c.val("K(-1)", elliptic_k(-1.0), 1.311_028_777_146_06, 1e-14);
    c.val("K(1)", elliptic_k(1.0), f64::INFINITY, 0.0);
    c.val("K(2)", elliptic_k(2.0), f64::NAN, 0.0);
    c.val("E(.5)", elliptic_e(0.5), 1.350_643_881_047_675_5, 1e-14);
    c.val("E(-1)", elliptic_e(-1.0), 1.910_098_894_513_856, 1e-14);
    c.val("E(0)", elliptic_e(0.0), PI / 2.0, 1e-15);
    c.val(
        "F(pi/2,.5)",
        elliptic_f(PI / 2.0, 0.5),
        1.854_074_677_301_371_9,
        1e-13,
    );
    // Landen: F(phi + n pi) = F(phi) + 2 n K
    let (phi, m) = (0.7, 0.3);
    c.val(
        "F(phi+3pi)",
        elliptic_f(phi + 3.0 * PI, m),
        elliptic_f(phi, m) + 6.0 * elliptic_k(m),
        1e-12,
    );
    let (sn, cn, dn) = jacobi_elliptic(0.5, 0.5);
    c.val("sn", sn, 0.470_750_473_655_657_3, 1e-12);
    c.val("cn", cn, 0.882_266_394_890_440_3, 1e-12);
    c.val("dn", dn, 0.942_972_425_777_385_7, 1e-12);
    c.val(
        "2F1(1,1,2,-1)",
        hyp2f1(1.0, 1.0, 2.0, -1.0),
        2f64.ln(),
        1e-13,
    );
    c.val(
        "2F1(1,1,2,.5)",
        hyp2f1(1.0, 1.0, 2.0, 0.5),
        1.386_294_361_119_890_6,
        1e-13,
    );
    c.val(
        "2F1(.5,.5,1,.99)",
        hyp2f1(0.5, 0.5, 1.0, 0.99),
        2.0 / PI * elliptic_k(0.99),
        1e-10,
    );
    c.val(
        "2F1(.5,.5,1,-10)",
        hyp2f1(0.5, 0.5, 1.0, -10.0),
        2.0 / PI * elliptic_k(-10.0),
        1e-10,
    );
    c.val("2F1(1,1,3,1)", hyp2f1(1.0, 1.0, 3.0, 1.0), 2.0, 1e-10);
    c.val("2F1(1,1,1,.5)", hyp2f1(1.0, 1.0, 1.0, 0.5), 2.0, 1e-13);
    c.val(
        "2F1(1,1,2,1)",
        hyp2f1(1.0, 1.0, 2.0, 1.0),
        f64::INFINITY,
        0.0,
    );
    c.val("2F1(1,1,2,2)", hyp2f1(1.0, 1.0, 2.0, 2.0), f64::NAN, 0.0);
    c.val("2F1(1,1,0,.5)", hyp2f1(1.0, 1.0, 0.0, 0.5), f64::NAN, 0.0);
    c.val(
        "2F1(1,1,2,-1e6)",
        hyp2f1(1.0, 1.0, 2.0, -1e6),
        (1e6f64 + 1.0).ln() / 1e6,
        1e-10,
    );
    c.val(
        "1F1(1,2,3)",
        hyp1f1(1.0, 2.0, 3.0),
        (3f64.exp() - 1.0) / 3.0,
        1e-13,
    );
    c.val(
        "1F1(1,2,-3)",
        hyp1f1(1.0, 2.0, -3.0),
        (1.0 - (-3f64).exp()) / 3.0,
        1e-13,
    );
    c.val(
        "1F1(.5,1.5,-500)",
        hyp1f1(0.5, 1.5, -500.0),
        PI.sqrt() / (2.0 * 500f64.sqrt()),
        1e-12,
    );
    c.val("1F1(1,2,800)", hyp1f1(1.0, 2.0, 800.0), f64::INFINITY, 0.0);
    c.val(
        "E1(1)",
        exp_integral_e1(1.0),
        0.219_383_934_395_520_29,
        1e-13,
    );
    c.val(
        "E1(1e-300)",
        exp_integral_e1(1e-300),
        690.775_527_898_213_7 - 0.577_215_664_901_532_9,
        1e-12,
    );
    c.val("E1(0)", exp_integral_e1(0.0), f64::INFINITY, 0.0);
    c.val(
        "Ei(1)",
        exp_integral_ei(1.0),
        1.895_117_816_355_936_8,
        1e-13,
    );
    c.val(
        "Ei(-1)",
        exp_integral_ei(-1.0),
        -0.219_383_934_395_520_29,
        1e-13,
    );
    c.val("Ei(0)", exp_integral_ei(0.0), f64::NEG_INFINITY, 0.0);
    c.val("Si(1)", sine_integral(1.0), 0.946_083_070_367_183, 1e-13);
    c.val("Si(inf)", sine_integral(f64::INFINITY), PI / 2.0, 1e-15);
    c.val(
        "Ci(1)",
        cosine_integral(1.0),
        0.337_403_922_900_968_1,
        1e-13,
    );
    c.val("Ci(0)", cosine_integral(0.0), f64::NEG_INFINITY, 0.0);
    c.val("li(2)", log_integral(2.0), 1.045_163_780_117_492_8, 1e-13);
    c.val("li(1)", log_integral(1.0), f64::NEG_INFINITY, 0.0);
    c.finish();
}

#[test]
fn elementary_helpers_do_not_overflow_or_lose_digits() {
    let mut c = Check(Vec::new());
    c.val("sinc(0)", sinc(0.0), 1.0, 0.0);
    c.val("sinc(inf)", sinc(f64::INFINITY), 0.0, 0.0);
    c.val("logit(0)", logit(0.0), f64::NEG_INFINITY, 0.0);
    c.val("logit(1)", logit(1.0), f64::INFINITY, 0.0);
    c.val("logit(2)", logit(2.0), f64::NAN, 0.0);
    c.val("logit(-1)", logit(-1.0), f64::NAN, 0.0);
    c.val("logit(.5)", logit(0.5), 0.0, 0.0);
    c.val(
        "logit(1e-300)",
        logit(1e-300),
        -690.775_527_898_213_7,
        1e-14,
    );
    c.val("logit(sigmoid(3))", logit(sigmoid(3.0)), 3.0, 1e-14);
    c.val("sigmoid(-800)", sigmoid(-800.0), 0.0, 0.0);
    c.val("sigmoid(800)", sigmoid(800.0), 1.0, 0.0);
    c.val("softplus(800)", softplus(800.0), 800.0, 0.0);
    c.val("softplus(-800)", softplus(-800.0), 0.0, 0.0);
    c.val("softplus(-40)", softplus(-40.0), (-40f64).exp(), 1e-12);
    c.val(
        "bernoulli(20)",
        bernoulli_number(20),
        -529.124_242_424_242_4,
        1e-13,
    );
    c.val("bernoulli(5)", bernoulli_number(5), 0.0, 0.0);
    c.val("legendre(2,2)", legendre_p(2, 2.0), 5.5, 1e-15);
    c.finish();
    // Degenerate orders must not panic or hang.
    for n in [0_u32, 1, 2, 100, 1000, u32::MAX >> 16] {
        let _ = bernoulli_number(n);
        let _ = legendre_p(n.min(5000), 0.5);
        let _ = hermite_h(n.min(500), 10.0);
        let _ = laguerre_l(n.min(500), 2.0);
    }
    let _ = bernoulli_number(u32::MAX);
    let _ = double_factorial(300);
    let _ = ln_factorial(1e300);
}
