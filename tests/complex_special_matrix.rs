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
