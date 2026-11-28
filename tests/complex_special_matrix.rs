//! Phase 4 integration tests: full complex special-functions matrix,
//! real-axis consistency, Wronskians, functional equations, ODEs, and Session API.

use std::f64::consts::PI;

use num_complex::Complex64;
use rssn::api::Config;
use rssn::api::Session;
use rssn::kernels::complex::cauchy_fft_derivatives;
use rssn::kernels::fresnel;
use rssn::kernels::special;
use rssn::kernels::special_complex;
use rssn::kernels::special_complex_ext as ext;
use rssn::kernels::special_ext;

fn near(
    got: Complex64,
    want: Complex64,
    tol: f64,
    label: &str,
) {
    let err = (got - want).norm();
    assert!(
        err <= tol * (1.0 + want.norm()),
        "{label}: got {got}, want {want}, err {err}"
    );
}

fn assert_carlson_rg_permutations(
    [x, y, z]: [Complex64; 3],
    expected: Complex64,
    tolerance: f64,
) {
    for [a, b, c] in [
        [x, y, z],
        [x, z, y],
        [y, x, z],
        [y, z, x],
        [z, x, y],
        [z, y, x],
    ] {
        near(
            ext::carlson_rg_complex(a, b, c),
            expected,
            tolerance,
            "Carlson R_G",
        );
    }
}

#[test]
fn carlson_rg_zero_parameters_follow_the_exact_limits() {
    let zero = Complex64::new(0.0, 0.0);
    let one = Complex64::new(1.0, 0.0);
    let z = Complex64::new(2.0, 3.0);
    for (argument, expected) in [
        ([one, one, zero], Complex64::new(PI / 4.0, 0.0)),
        ([z, z, zero], (PI / 4.0) * z.sqrt()),
        ([z, zero, zero], 0.5 * z.sqrt()),
        ([zero, zero, zero], zero),
    ] {
        assert_carlson_rg_permutations(argument, expected, 1e-15);
    }
    near(
        ext::elliptic_e_complex(zero),
        Complex64::new(PI / 2.0, 0.0),
        1e-15,
        "E(0)",
    );
    near(ext::elliptic_e_complex(one), one, 1e-15, "E(1)");
    for sign in [-1.0_f64, 1.0] {
        let cut = Complex64::new(-1.0, 0.0_f64.copysign(sign));
        assert_carlson_rg_permutations(
            [cut, cut, zero],
            Complex64::new(0.0, sign * PI / 4.0),
            1e-15,
        );
        assert_carlson_rg_permutations([cut, zero, zero], Complex64::new(0.0, sign * 0.5), 1e-15);
    }

    // The factor of one half acts on the root, after the subnormal input is scaled.
    let smallest = f64::from_bits(1);
    let result = ext::carlson_rg_complex(Complex64::new(0.0, smallest), zero, zero);
    let component = 0.5 * smallest.sqrt() * std::f64::consts::FRAC_1_SQRT_2;
    assert!((result.re / component - 1.0).abs() <= 2.0 * f64::EPSILON);
    assert!((result.im / component - 1.0).abs() <= 2.0 * f64::EPSILON);
}

#[test]
fn carlson_rg_permutations_match_the_defining_integral() {
    let zero = Complex64::new(0.0, 0.0);
    let one = Complex64::new(1.0, 0.0);
    // Independent quadrature of DLMF 19.16.2_5, with separate shifted square roots.
    for (args, expected) in [
        (
            [Complex64::new(1.0, 0.4), Complex64::new(3.0, -0.7), zero],
            Complex64::new(1.097_946_806_798_013_6, -0.022_302_029_297_672_13),
        ),
        (
            [Complex64::new(-1.0, 0.07), Complex64::new(2.0, 0.5), zero],
            Complex64::new(0.546_644_444_162_021_3, 0.332_204_453_516_196_76),
        ),
        (
            [Complex64::new(-1.0, 0.1), Complex64::new(-1.0, -0.1), zero],
            Complex64::new(-1.192_068_147_045_497_5, 0.0),
        ),
        (
            [Complex64::new(-1.0, 0.1), Complex64::new(-1.0, 0.1), one],
            Complex64::new(0.229_368_399_633_498_2, 0.567_058_829_970_660_7),
        ),
    ] {
        assert_carlson_rg_permutations(args, expected, 4e-13);
        assert_carlson_rg_permutations(args.map(|z| z.conj()), expected.conj(), 4e-13);
    }
}

#[test]
fn gamma_beta_incomplete_and_barnes_g_identities() {
    let z = Complex64::new(0.4, 0.7);
    let one = Complex64::new(1.0, 0.0);

    // Reflection: Gamma(z) * Gamma(1 - z) = pi / sin(pi * z)
    let lhs = special_complex::gamma_complex(z) * special_complex::gamma_complex(one - z);
    let rhs = Complex64::new(PI, 0.0) / (PI * z).sin();
    near(lhs, rhs, 1e-12, "gamma reflection");

    // Reciprocal gamma: rgamma(z) * Gamma(z) = 1, and rgamma(-3) = 0
    near(
        ext::rgamma_complex(z) * special_complex::gamma_complex(z),
        one,
        1e-12,
        "rgamma",
    );
    near(
        ext::rgamma_complex(Complex64::new(-3.0, 0.0)),
        Complex64::new(0.0, 0.0),
        1e-14,
        "rgamma(-3)",
    );

    // Beta: B(a, b) = Gamma(a)*Gamma(b)/Gamma(a+b)
    let a = Complex64::new(1.3, 0.4);
    let b = Complex64::new(2.1, -0.5);
    let b_want = special_complex::gamma_complex(a) * special_complex::gamma_complex(b)
        / special_complex::gamma_complex(a + b);
    near(ext::beta_complex(a, b), b_want, 1e-12, "beta");

    // Incomplete gamma: gamma(a, z) + Gamma(a, z) = Gamma(a)
    let z2 = Complex64::new(2.5, 1.1);
    let sum_g =
        ext::lower_incomplete_gamma_complex(a, z2) + ext::upper_incomplete_gamma_complex(a, z2);
    near(
        sum_g,
        special_complex::gamma_complex(a),
        1e-11,
        "incomplete gamma sum",
    );

    // Barnes G: G(z + 1) = Gamma(z) * G(z)
    near(
        ext::barnes_g_complex(z2 + one),
        special_complex::gamma_complex(z2) * ext::barnes_g_complex(z2),
        1e-11,
        "barnes G recurrence",
    );

    // Polygamma: psi'(z) + psi'(1 - z) = pi^2 / sin^2(pi*z)
    let s = (PI * z).sin();
    near(
        ext::polygamma_complex(1, z) + ext::polygamma_complex(1, one - z),
        Complex64::new(PI * PI, 0.0) / (s * s),
        1e-11,
        "trigamma reflection",
    );
}

#[test]
fn error_fresnel_and_exponential_integrals() {
    let z = Complex64::new(0.8, 0.6);
    // erf(z) + erfc(z) = 1
    near(
        fresnel::erf_complex(z) + ext::erfc_complex(z),
        Complex64::new(1.0, 0.0),
        1e-12,
        "erf + erfc",
    );
    // erfcx(z) = exp(z^2) * erfc(z)
    near(
        ext::erfcx_complex(z),
        (z * z).exp() * ext::erfc_complex(z),
        1e-12,
        "erfcx",
    );
    // Ei(i*x) = Ci(x) + i*(Si(x) + pi/2) for x > 0
    let x = 2.3;
    let ix = Complex64::new(0.0, x);
    let (si, ci) = ext::sin_cos_integral_complex(Complex64::new(x, 0.0));
    let ei_ix = ext::exp_integral_ei_complex(ix);
    let want = ci + Complex64::new(0.0, 1.0) * (si + 0.5 * PI);
    near(ei_ix, want, 1e-11, "Ei(ix) = Ci(x) + i(Si(x)+pi/2)");

    // Real-axis consistency for E1
    near(
        ext::exp_integral_e1_complex(Complex64::new(1.5, 0.0)),
        Complex64::new(special::exp_integral_e1(1.5), 0.0),
        1e-12,
        "E1 real axis",
    );
}

#[test]
fn bessel_hankel_airy_kelvin_and_wronskians() {
    let z = Complex64::new(2.2, 0.9);
    let nu = 1.3;
    // H1 + H2 = 2 * J
    near(
        ext::hankel_h1_complex(nu, z) + ext::hankel_h2_complex(nu, z),
        2.0 * special_complex::bessel_j_complex(nu, z),
        1e-12,
        "H1 + H2 = 2J",
    );

    // Bessel Wronskian via recurrence: J_{nu+1}(z)*Y_nu(z) - J_nu(z)*Y_{nu+1}(z) = 2/(pi*z)
    let wronsk = special_complex::bessel_j_complex(nu + 1.0, z) * ext::bessel_y_complex(nu, z)
        - special_complex::bessel_j_complex(nu, z) * ext::bessel_y_complex(nu + 1.0, z);
    near(
        wronsk,
        Complex64::new(2.0 / PI, 0.0) / z,
        1e-11,
        "J/Y Wronskian",
    );

    // Airy Wronskian: Ai(z)*Bi'(z) - Ai'(z)*Bi(z) = 1/pi
    for pt in [
        Complex64::new(0.5, 0.3),
        Complex64::new(-4.5, 1.2),
        Complex64::new(7.5, 2.0),
    ] {
        let (ai, aip) = ext::airy_ai_complex(pt);
        let (bi, bip) = ext::airy_bi_complex(pt);
        near(
            ai * bip - aip * bi,
            Complex64::new(1.0 / PI, 0.0),
            1e-10,
            "Airy Wronskian",
        );
    }

    // Airy ODE: Ai''(z) = z * Ai(z) via Cauchy-FFT
    let z0 = Complex64::new(1.2, -0.4);
    let derivs = cauchy_fft_derivatives(|w| ext::airy_ai_complex(w).0, z0, 0.4, 2, 32);
    near(derivs[2], z0 * derivs[0], 1e-10, "Airy ODE");

    // Real-axis consistency for Ai, Bi
    near(
        ext::airy_ai_complex(Complex64::new(1.5, 0.0)).0,
        Complex64::new(special::airy_ai(1.5), 0.0),
        1e-12,
        "Ai real axis",
    );
    near(
        ext::airy_bi_complex(Complex64::new(1.5, 0.0)).0,
        Complex64::new(special::airy_bi(1.5), 0.0),
        1e-12,
        "Bi real axis",
    );
}

#[test]
fn zeta_polylog_elliptic_and_modular_forms() {
    // First nontrivial zero of Riemann zeta: 0.5 + 14.134725141734693i
    let rho1 = Complex64::new(0.5, 14.134_725_141_734_693);
    near(
        ext::riemann_zeta_complex(rho1),
        Complex64::new(0.0, 0.0),
        1e-10,
        "zeta zero 1",
    );

    // Negative even integer trivial zeros: zeta(-2) = 0, zeta(-1) = -1/12
    near(
        ext::riemann_zeta_complex(Complex64::new(-2.0, 0.0)),
        Complex64::new(0.0, 0.0),
        1e-12,
        "zeta(-2)",
    );
    near(
        ext::riemann_zeta_complex(Complex64::new(-1.0, 0.0)),
        Complex64::new(-1.0 / 12.0, 0.0),
        1e-12,
        "zeta(-1)",
    );

    // Jacobi elliptic identities: sn^2 + cn^2 = 1, dn^2 + m*sn^2 = 1
    let u = Complex64::new(0.7, 0.4);
    let m = Complex64::new(0.3, 0.2);
    let (sn, cn, dn) = ext::jacobi_elliptic_complex(u, m);
    near(
        sn * sn + cn * cn,
        Complex64::new(1.0, 0.0),
        1e-12,
        "sn^2 + cn^2 = 1",
    );
    near(
        dn * dn + m * sn * sn,
        Complex64::new(1.0, 0.0),
        1e-12,
        "dn^2 + m*sn^2 = 1",
    );

    // Jacobi theta identity: theta_3^4 = theta_2^4 + theta_4^4
    let q = Complex64::new(0.1, 0.05);
    let t2 = ext::jacobi_theta_complex(2, Complex64::new(0.0, 0.0), q);
    let t3 = ext::jacobi_theta_complex(3, Complex64::new(0.0, 0.0), q);
    let t4 = ext::jacobi_theta_complex(4, Complex64::new(0.0, 0.0), q);
    near(
        t3.powi(4),
        t2.powi(4) + t4.powi(4),
        1e-12,
        "Jacobi theta identity",
    );

    // Weierstrass ODE: (wp')^2 = 4*wp^3 - g2*wp - g3
    let w1 = Complex64::new(1.0, 0.0);
    let w2 = Complex64::new(0.2, 1.1);
    let z_w = Complex64::new(0.3, 0.4);
    let wp = ext::weierstrass_p_complex(z_w, w1, w2);
    let wpp = ext::weierstrass_p_prime_complex(z_w, w1, w2);
    let (g2, g3) = ext::weierstrass_invariants_complex(w1, w2);
    near(
        wpp * wpp,
        4.0 * wp * wp * wp - g2 * wp - g3,
        1e-10,
        "Weierstrass ODE",
    );

    // Modular forms: j(i) = 1728, j(tau) = j(-1/tau)
    near(
        ext::klein_j_complex(Complex64::new(0.0, 1.0)),
        Complex64::new(1728.0, 0.0),
        1e-10,
        "j(i) = 1728",
    );
    let tau = Complex64::new(0.3, 1.2);
    near(
        ext::klein_j_complex(tau),
        ext::klein_j_complex(-Complex64::new(1.0, 0.0) / tau),
        1e-10,
        "j(-1/tau) = j(tau)",
    );
    let i_unit = Complex64::new(0.0, 1.0);
    near(
        ext::dedekind_eta_complex(-Complex64::new(1.0, 0.0) / tau),
        (-i_unit * tau).sqrt() * ext::dedekind_eta_complex(tau),
        1e-11,
        "Dedekind eta modular",
    );
}

#[test]
fn wright_omega_and_session_operators() {
    let z = Complex64::new(1.5, -2.1);
    let w = ext::wright_omega_complex(z);
    near(w + w.ln(), z, 1e-12, "wright omega defining equation");

    let session = Session::new();
    for expr in [
        "riemann_zeta_complex(2, 0)",
        "airy_ai_complex(1, 0)",
        "klein_j_complex(0, 1)",
        "wright_omega_complex(1, 1)",
        "jacobi_sn_complex(0.5, 0.2, 0.3, 0.1)",
    ] {
        let term = session.parse(expr).unwrap();
        let ans = session.compute(term, &Config::default()).unwrap();
        assert!(ans.reduced, "{expr} did not reduce");
    }
    let _ = special_ext::elliptic_k(0.5);
}

#[test]
fn complex_bessel_integer_series_reject_unrepresentable_orders() {
    let z = Complex64::new(0.0, std::f64::consts::PI);
    for order in [172.0, 1e8, -1e308, f64::INFINITY] {
        assert!(ext::bessel_y_complex(order, z).is_nan());
        assert!(ext::bessel_k_complex(order, z).is_nan());
        assert!(ext::bessel_ke_complex(order, z).is_nan());
    }
}

#[test]
fn adversarial_and_boundary_complex_special_functions() {
    // Ei(10) on positive real axis (E_1(-10) near cut)
    let ei10 = ext::exp_integral_ei_complex(Complex64::new(10.0, 0.0));
    assert!(
        (ei10.re - 2_492.228_976_241_877_7).abs() < 1e-8 && ei10.im.abs() < 1e-10,
        "Ei(10)={ei10}"
    );

    // Ci(-5) branch shift
    let (_, ci_m5) = ext::sin_cos_integral_complex(Complex64::new(-5.0, 0.0));
    let (_, ci_p5) = ext::sin_cos_integral_complex(Complex64::new(5.0, 0.0));
    near(
        ci_m5,
        ci_p5 + Complex64::new(0.0, std::f64::consts::PI),
        1e-12,
        "Ci(-5) = Ci(5) + i*pi",
    );

    // Elliptic F(pi | 0.5) = 2*K(0.5)
    let m = Complex64::new(0.5, 0.0);
    near(
        ext::elliptic_f_complex(Complex64::new(std::f64::consts::PI, 0.0), m),
        2.0 * ext::elliptic_k_complex(m),
        1e-12,
        "F(pi | m) = 2*K(m)",
    );

    // Parabolic cylinder D_2(z) = (z^2 - 1) * exp(-z^2/4) at integer nu=2 and large z=12
    let z_pc = Complex64::new(12.0, 0.5);
    let exact_d2 = (z_pc * z_pc - 1.0) * (-0.25 * z_pc * z_pc).exp();
    near(
        ext::parabolic_cylinder_d_complex(Complex64::new(2.0, 0.0), z_pc),
        exact_d2,
        1e-12 * exact_d2.norm().max(1.0),
        "D_2(12+0.5i)",
    );

    // Binomial(-3, 2) = 6
    near(
        ext::binomial_complex(Complex64::new(-3.0, 0.0), Complex64::new(2.0, 0.0)),
        Complex64::new(6.0, 0.0),
        1e-12,
        "C(-3, 2) = 6",
    );

    // Polygamma reflection for n=3 at z = -0.3 + 0.4i
    let z_pg = Complex64::new(-0.3, 0.4);
    let p3_l = ext::polygamma_complex(3, z_pg);
    let p3_r = ext::polygamma_complex(3, Complex64::new(1.0, 0.0) - z_pg);
    assert!(p3_l.is_finite() && p3_r.is_finite());

    // Large |Im(s)| Riemann zeta
    let z_high = ext::riemann_zeta_complex(Complex64::new(0.5, 60.0));
    assert!(z_high.is_finite() && z_high.norm() < 5.0);

    // Scaled Bessel K_0 at Re(z) = 800 (no exp overflow)
    let ke800 = ext::bessel_ke_complex(0.0, Complex64::new(800.0, 1.0));
    assert!(ke800.is_finite() && ke800.norm() > 0.0);

    // Mittag-Leffler E_{1,1}(z) = exp(z) and E_{2,1}(-z^2) = cos(z)
    let z_ml = Complex64::new(0.7, -0.4);
    near(
        ext::mittag_leffler_complex(Complex64::new(1.0, 0.0), Complex64::new(1.0, 0.0), z_ml),
        z_ml.exp(),
        1e-12,
        "E_{1,1}(z) = exp(z)",
    );
    near(
        ext::mittag_leffler_complex(
            Complex64::new(2.0, 0.0),
            Complex64::new(1.0, 0.0),
            -z_ml * z_ml,
        ),
        z_ml.cos(),
        1e-12,
        "E_{2,1}(-z^2) = cos(z)",
    );

    // ComplexBall special functions
    let b = rssn::kernels::complex_ball::ComplexBall::new(Complex64::new(1.5, 0.5), 1e-4).unwrap();
    assert!(
        b.erf()
            .unwrap()
            .contains(rssn::kernels::fresnel::erf_complex(b.mid))
    );
    assert!(
        b.gamma()
            .unwrap()
            .contains(rssn::kernels::special_complex::gamma_complex(b.mid))
    );
    assert!(
        b.bessel_j(0.0)
            .unwrap()
            .contains(rssn::kernels::special_complex::bessel_j_complex(0.0, b.mid))
    );

    // MpComplex special functions
    let (re_g, im_g) = rssn::kernels::mp::eval_complex_text("gamma(1 + I)", 25).unwrap();
    assert!(!re_g.is_empty() && !im_g.is_empty());
    let ((re_lo, re_hi), (im_lo, im_hi)) =
        rssn::kernels::mp::enclosure_complex_text("erf(0.5 + 0.5*I) + erfc(0.5 + 0.5*I)", 25)
            .unwrap();
    let re_lo_f: f64 = re_lo.parse().unwrap();
    let re_hi_f: f64 = re_hi.parse().unwrap();
    let im_lo_f: f64 = im_lo.parse().unwrap();
    let im_hi_f: f64 = im_hi.parse().unwrap();
    assert!(re_lo_f <= 1.0 && 1.0 <= re_hi_f && (re_hi_f - re_lo_f) < 1e-20);
    assert!(im_lo_f <= 0.0 && 0.0 <= im_hi_f && (im_hi_f - im_lo_f) < 1e-20);

    // Graph::eval_complex on symbolic special ops
    let mut g = rssn::graph::Graph::new();
    let _eng = rssn::graph::Engine::install(
        &mut g,
        &[
            rssn::rules::special::special(),
            rssn::rules::complex::complex(),
        ],
    )
    .unwrap();
    let t = g.parse("gamma(1 + I)").unwrap();
    let v = g
        .eval_complex(t, &std::collections::HashMap::new())
        .unwrap();
    near(
        v,
        rssn::kernels::special_complex::gamma_complex(Complex64::new(1.0, 1.0)),
        1e-12,
        "eval_complex gamma(1+I)",
    );
}
