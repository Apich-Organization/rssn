//! Complex-argument gamma, digamma, Lambert W and Bessel J/I, checked against
//! published values, closed forms (half-integer orders) and functional
//! identities.

use std::f64::consts::PI;

use num_complex::Complex64;
use rssn::kernels::special::bessel_i;
use rssn::kernels::special::bessel_j;
use rssn::kernels::special::digamma_numerical;
use rssn::kernels::special::gamma_numerical;
use rssn::kernels::special::lambert_w;
use rssn::kernels::special_complex::*;

fn c(
    re: f64,
    im: f64,
) -> Complex64 {
    Complex64::new(re, im)
}

fn rel(
    got: Complex64,
    want: Complex64,
) -> f64 {
    (got - want).norm() / want.norm().max(1e-300)
}

fn grid() -> Vec<Complex64> {
    let mut v = Vec::new();
    for re in [-7.3, -2.5, -0.5, 0.01, 0.5, 1.0, 3.7, 12.5, 40.0] {
        for im in [-9.0, -1.0, -1e-3, 0.0, 0.25, 2.0, 15.0] {
            v.push(c(re, im));
        }
    }
    v
}

#[test]
fn gamma_published_values() {
    let g = gamma_complex(c(0.0, 1.0));
    assert!(
        rel(g, c(-0.154_949_828_301_810_68, -0.498_015_668_118_356)) < 1e-13,
        "{g}"
    );
    let g = gamma_complex(c(1.0, 1.0));
    assert!(
        rel(g, c(0.498_015_668_118_356, -0.154_949_828_301_810_68)) < 1e-13,
        "{g}"
    );
    // |Γ(1/2 + i t)|² = π / cosh(π t)
    for t in [0.1, 1.0, 3.0, 8.0] {
        let g = gamma_complex(c(0.5, t));
        let want = PI / (PI * t).cosh();
        assert!((g.norm_sqr() - want).abs() <= 1e-13 * want, "t = {t}");
    }
}

#[test]
fn gamma_agrees_with_the_real_routine_and_poles() {
    for x in [-3.5, -0.5, 0.1, 1.0, 2.5, 7.0, 25.5, 100.25] {
        let g = gamma_complex(c(x, 0.0));
        assert!(rel(g, c(gamma_numerical(x), 0.0)) < 1e-12, "x = {x}");
    }
    for n in 0..5 {
        assert!(gamma_complex(c(-f64::from(n), 0.0)).re.is_infinite());
    }
    assert!(gamma_complex(c(f64::NAN, 0.0)).re.is_nan());
    assert!(gamma_complex(c(0.0, f64::INFINITY)).re.is_nan());
}

#[test]
fn gamma_functional_equations() {
    for z in grid() {
        let (g, g1) = (gamma_complex(z), gamma_complex(z + 1.0));
        if !g.norm().is_finite() || g.norm() > 1e280 || g.norm() < 1e-280 {
            continue;
        }
        // Γ(z + 1) = z Γ(z)
        assert!(rel(g1, z * g) < 5e-13, "recurrence at {z}");
        // Γ(z) Γ(1 - z) = π / sin(π z)
        let refl = g * gamma_complex(1.0 - z);
        let want = PI / (PI * z).sin();
        if want.norm().is_finite() && want.norm() < 1e250 && refl.norm().is_finite() {
            assert!(rel(refl, want) < 5e-12, "reflection at {z}");
        }
        // conjugate symmetry
        assert!(
            rel(gamma_complex(z.conj()), g.conj()) < 1e-13,
            "conj at {z}"
        );
    }
}

#[test]
fn ln_gamma_is_the_analytic_branch() {
    for z in grid() {
        if z.im == 0.0 && z.re <= 0.0 {
            continue;
        }
        // ln Γ(z + 1) = ln Γ(z) + ln z with the principal ln z
        let lhs = ln_gamma_complex(z + 1.0);
        let rhs = ln_gamma_complex(z) + z.ln();
        assert!(
            (lhs - rhs).norm() < 1e-11 * (1.0 + lhs.norm()),
            "at {z}: {lhs} vs {rhs}"
        );
        // exp(ln Γ) = Γ
        let g = gamma_complex(z);
        if g.norm() < 1e250 && g.norm() > 1e-250 {
            assert!(
                rel(ln_gamma_complex(z).exp(), g) < 1e-9 * (1.0 + ln_gamma_complex(z).im.abs()),
                "exp at {z}"
            );
        }
    }
    // large arguments
    let z = c(1e5, 3e4);
    let l = ln_gamma_complex(z);
    let stirling = (z - 0.5) * z.ln() - z + 0.5 * (2.0 * PI).ln() + 1.0 / (12.0 * z);
    assert!((l - stirling).norm() < 1e-9 * l.norm());
    assert!(ln_gamma_complex(c(0.0, 0.0)).re.is_infinite());
}

#[test]
fn digamma_values_and_identities() {
    let psi_i = digamma_complex(c(0.0, 1.0));
    assert!(
        rel(psi_i, c(0.094_650_320_622_476_98, 2.076_674_047_468_581)) < 1e-12,
        "{psi_i}"
    );
    for x in [-3.5, -0.5, 0.1, 1.0, 2.5, 7.0, 50.0] {
        assert!(
            rel(digamma_complex(c(x, 0.0)), c(digamma_numerical(x), 0.0)) < 1e-12,
            "x = {x}: {} vs {}",
            digamma_complex(c(x, 0.0)),
            digamma_numerical(x)
        );
    }
    for z in grid() {
        if z.im == 0.0 && z.re <= 0.0 {
            continue;
        }
        // ψ(z + 1) = ψ(z) + 1 / z
        let lhs = digamma_complex(z + 1.0);
        let rhs = digamma_complex(z) + 1.0 / z;
        assert!(
            (lhs - rhs).norm() < 1e-11 * (1.0 + lhs.norm()),
            "recurrence at {z}"
        );
        // ψ(1 - z) - ψ(z) = π cot(π z)
        let d = digamma_complex(1.0 - z) - digamma_complex(z);
        let cot = PI / (PI * z).tan();
        if cot.norm() < 1e8 {
            assert!(
                (d - cot).norm() < 1e-10 * (1.0 + cot.norm()),
                "reflection at {z}"
            );
        }
    }
    assert!(digamma_complex(c(-2.0, 0.0)).re.is_nan());
    assert!(digamma_complex(c(f64::NAN, 1.0)).re.is_nan());
}

#[test]
fn lambert_w_branches() {
    let w = lambert_w_complex(c(1.0, 0.0), 0);
    assert!((w.re - 0.567_143_290_409_783_8).abs() < 1e-13, "{w}");
    assert!((lambert_w_complex(c(-(-1f64).exp(), 0.0), 0).re + 1.0).abs() < 1e-7);
    let w = lambert_w_complex(c(-0.3, 0.0), -1);
    assert!(
        (w.re + 1.781_337_023_421_627_6).abs() < 1e-12 && w.im.abs() < 1e-12,
        "{w}"
    );
    let w = lambert_w_complex(c(0.0, 1.0), 0);
    assert!(
        rel(w, c(0.374_699_020_737_117_5, 0.576_412_723_031_435_6)) < 1e-10,
        "{w}"
    );
    let w1 = lambert_w_complex(c(1.0, 0.0), 1);
    assert!(
        rel(w1, c(-1.533_913_319_793_574_5, 4.375_185_153_061_898)) < 1e-12,
        "{w1}"
    );
    assert!(lambert_w_complex(c(0.0, 0.0), 0).norm() == 0.0);
    assert_eq!(lambert_w_complex(c(0.0, 0.0), 2).re, f64::NEG_INFINITY);
    // the real routine agrees on the real branches
    for x in [-0.36, -0.2, -0.05, 0.0, 0.7, 3.0, 100.0, 1e6] {
        assert!(
            (lambert_w_complex(c(x, 0.0), 0).re - lambert_w(x, true)).abs()
                < 1e-9 * (1.0 + lambert_w(x, true).abs()),
            "W0({x})"
        );
    }
    // defining equation on a grid, all branches
    for z in grid() {
        for k in -3..=3 {
            let w = lambert_w_complex(z, k);
            let residual = (w * w.exp() - z).norm();
            assert!(
                residual <= 1e-11 * (1.0 + z.norm()),
                "k = {k}, z = {z}, w = {w}, residual {residual}"
            );
            // The imaginary part lies in the strip of branch k (loosely).
            let lo = (2.0 * f64::from(k) - 2.0) * PI;
            let hi = (2.0 * f64::from(k) + 2.0) * PI;
            assert!(
                w.im > lo - 1e-9 && w.im < hi + 1e-9,
                "strip k = {k}, z = {z}, w = {w}"
            );
        }
    }
}

#[test]
fn bessel_published_value() {
    let j = bessel_j_complex(0.0, c(1.0, 1.0));
    assert!(
        rel(j, c(0.937_608_476_806_029_3, -0.496_529_947_609_122_1)) < 1e-12,
        "{j}"
    );
    // J0(i y) = I0(y)
    let j = bessel_j_complex(0.0, c(0.0, 10.0));
    assert!(rel(j, c(2_815.716_628_466_254, 0.0)) < 1e-10, "{j}");
}

#[test]
fn bessel_half_integer_closed_forms() {
    // J_{1/2}(z) = sqrt(2/(πz)) sin z,   J_{-1/2}(z) = sqrt(2/(πz)) cos z,
    // J_{3/2}(z) = sqrt(2/(πz)) (sin z / z - cos z),
    // I_{1/2}(z) = sqrt(2/(πz)) sinh z,  I_{-1/2}(z) = sqrt(2/(πz)) cosh z.
    for z in [
        c(0.3, 0.1),
        c(2.0, -1.5),
        c(-1.5, 0.7),
        c(6.0, 6.0),
        c(13.0, 0.5),
        c(15.0, 3.0),
        c(30.0, 5.0),
        c(-40.0, 2.0),
        c(-25.0, -6.0),
        c(0.5, 60.0),
        c(5.0, -50.0),
        c(100.0, 0.0),
    ] {
        let p = (2.0 / (PI * z)).sqrt();
        let tol = if z.norm() > 11.0 && z.norm() < 16.0 {
            1e-9
        } else {
            1e-11
        };
        let scale = |want: Complex64| rel_scaled(want, z);
        assert!(
            rel(bessel_j_complex(0.5, z), p * z.sin()) < tol * scale(p * z.sin()),
            "J_1/2({z})"
        );
        assert!(
            rel(bessel_j_complex(-0.5, z), p * z.cos()) < tol * scale(p * z.cos()),
            "J_-1/2({z})"
        );
        let j32 = p * (z.sin() / z - z.cos());
        assert!(
            rel(bessel_j_complex(1.5, z), j32) < tol * scale(j32) * 10.0,
            "J_3/2({z})"
        );
        assert!(
            rel(bessel_i_complex(0.5, z), p * z.sinh()) < tol * scale(p * z.sinh()),
            "I_1/2({z})"
        );
        assert!(
            rel(bessel_i_complex(-0.5, z), p * z.cosh()) < tol * scale(p * z.cosh()),
            "I_-1/2({z})"
        );
    }
}

/// Cancellation allowance: when the answer is far smaller than the largest
/// intermediate (`e^{|Im z|}`) the relative error grows by that ratio.
fn rel_scaled(
    want: Complex64,
    z: Complex64,
) -> f64 {
    let big = z.im.abs().exp() / (z.norm()).sqrt();
    (big / want.norm().max(1e-300)).clamp(1.0, 1e6)
}

#[test]
fn bessel_agrees_with_the_real_routines() {
    for nu in [0.0, 1.0, 2.0, 0.3, 2.5, 7.0] {
        for x in [0.5, 3.0, 9.0, 13.9, 14.1, 20.0, 60.0, 150.0] {
            let want = bessel_j(nu, x);
            let got = bessel_j_complex(nu, c(x, 0.0));
            assert!(
                got.im.abs() <= 1e-9 * (1.0 + want.abs()),
                "J_{nu}({x}) imaginary part {got}"
            );
            assert!(
                (got.re - want).abs() <= 1e-8 * (1.0 + want.abs()) + 1e-9 * x.powf(-0.5),
                "J_{nu}({x}) = {got}, real {want}"
            );
            let want = bessel_i(nu, x);
            let got = bessel_i_complex(nu, c(x, 0.0));
            assert!(
                (got.re - want).abs() <= 1e-8 * want.abs(),
                "I_{nu}({x}) = {got}, real {want}"
            );
        }
    }
}

#[test]
fn bessel_recurrence_and_symmetries() {
    for z in [
        c(1.0, 2.0),
        c(-3.0, 1.0),
        c(8.0, -4.0),
        c(20.0, 3.0),
        c(0.2, 0.1),
    ] {
        let nu = 0.7;
        // J_{ν-1}(z) + J_{ν+1}(z) = (2ν / z) J_ν(z)
        let lhs = bessel_j_complex(nu - 1.0, z) + bessel_j_complex(nu + 1.0, z);
        let rhs = bessel_j_complex(nu, z) * (2.0 * nu / z);
        assert!(
            (lhs - rhs).norm() <= 1e-9 * (1.0 + lhs.norm()),
            "J recurrence at {z}"
        );
        // I_{ν-1}(z) - I_{ν+1}(z) = (2ν / z) I_ν(z)
        let lhs = bessel_i_complex(nu - 1.0, z) - bessel_i_complex(nu + 1.0, z);
        let rhs = bessel_i_complex(nu, z) * (2.0 * nu / z);
        assert!(
            (lhs - rhs).norm() <= 1e-9 * (1.0 + lhs.norm()),
            "I recurrence at {z}"
        );
        // J_{-n} = (-1)^n J_n,  I_{-n} = I_n
        for n in 1..4 {
            let s = if n % 2 == 0 { 1.0 } else { -1.0 };
            assert!(
                rel(
                    bessel_j_complex(-f64::from(n), z),
                    bessel_j_complex(f64::from(n), z) * s
                ) < 1e-12
            );
            assert!(
                rel(
                    bessel_i_complex(-f64::from(n), z),
                    bessel_i_complex(f64::from(n), z)
                ) < 1e-12
            );
        }
        // I_ν(z) = e^{-iνπ/2} J_ν(i z) in the right half plane
        if z.re > 0.0 {
            let lhs = bessel_i_complex(nu, z);
            let rhs = (c(0.0, -0.5 * nu * PI)).exp() * bessel_j_complex(nu, c(0.0, 1.0) * z);
            assert!(rel(lhs, rhs) < 1e-8, "I-J relation at {z}");
        }
    }
    // origin, NaN
    assert_eq!(bessel_j_complex(0.0, c(0.0, 0.0)), c(1.0, 0.0));
    assert_eq!(bessel_j_complex(2.0, c(0.0, 0.0)), c(0.0, 0.0));
    assert!(bessel_j_complex(-0.5, c(0.0, 0.0)).re.is_infinite());
    assert!(bessel_j_complex(f64::NAN, c(1.0, 0.0)).re.is_nan());
    assert!(bessel_i_complex(1.0, c(f64::INFINITY, 0.0)).re.is_nan());
}

#[test]
fn nothing_panics_on_degenerate_input() {
    let specials = [
        0.0,
        -0.0,
        1.0,
        -1.0,
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
        1e300,
        -1e300,
        5e-324,
        170.5,
        -170.5,
    ];
    for &re in &specials {
        for &im in &specials {
            let z = c(re, im);
            let _ = gamma_complex(z);
            let _ = ln_gamma_complex(z);
            let _ = digamma_complex(z);
            for k in [-1, 0, 1, i32::MAX, i32::MIN] {
                let _ = lambert_w_complex(z, k);
            }
            for &nu in &specials {
                let _ = bessel_j_complex(nu, z);
                let _ = bessel_i_complex(nu, z);
            }
        }
    }
}
