//! Carlson `R_C`/`R_J`, elliptic integrals of the third kind, Jacobi-function integrals.

use rssn::kernels::special_ext::carlson_rd;
use rssn::kernels::special_ext::elliptic_k;
use rssn::sim::elliptic::*;

/// Composite Gauss-Legendre (8 nodes) quadrature of `f` on `[0, b]` with `n` panels.
fn quad<F: Fn(f64) -> f64>(
    f: F,
    b: f64,
    n: usize,
) -> f64 {
    const X: [f64; 4] = [
        0.183_434_642_495_649_8,
        0.525_532_409_916_329,
        0.796_666_477_413_626_7,
        0.960_289_856_497_536_3,
    ];
    const W: [f64; 4] = [
        0.362_683_783_378_362,
        0.313_706_645_877_887_3,
        0.222_381_034_453_374_5,
        0.101_228_536_290_376_26,
    ];
    let h = b / n as f64;
    let mut s = 0.0;
    for i in 0..n {
        let c = (i as f64 + 0.5) * h;
        for k in 0..4 {
            s += W[k] * (f(c + 0.5 * h * X[k]) + f(c - 0.5 * h * X[k]));
        }
    }
    s * 0.5 * h
}

#[test]
fn carlson_rj_reduces_to_rd_and_rc() {
    for &(x, y, z) in &[
        (0.3_f64, 1.2, 2.5),
        (1.0, 2.0, 3.0),
        (0.0, 1.0, 4.0),
        (5.0, 0.1, 0.7),
    ] {
        let rj = carlson_rj(x, y, z, z);
        let rd = carlson_rd(x, y, z);
        assert!((rj - rd).abs() < 1e-13 * rd.abs(), "{rj} {rd}");
    }
    // R_C(x, y) = int ... ; R_C(0, 1) = pi/2, R_C(x, x) = 1/sqrt(x)
    assert!((carlson_rc(0.0, 1.0) - std::f64::consts::FRAC_PI_2).abs() < 1e-14);
    assert!((carlson_rc(2.0, 2.0) - 2.0_f64.powf(-0.5)).abs() < 1e-14);
    // R_C(x, y) = atanh-type closed form for y < x: ln((sqrt x + sqrt(x-y))/sqrt y)/sqrt(x-y)
    let (x, y) = (3.0_f64, 0.5_f64);
    let exact = ((x.sqrt() + (x - y).sqrt()) / y.sqrt()).ln() / (x - y).sqrt();
    assert!((carlson_rc(x, y) - exact).abs() < 1e-14);
}

#[test]
fn incomplete_third_kind_matches_quadrature() {
    for &m in &[0.0_f64, 0.3, 0.9, 0.999] {
        for &n in &[-3.0_f64, -0.5, 0.0, 0.4, 0.8] {
            for &phi in &[0.2_f64, 1.0, 1.5, 2.5, 4.0, 7.5, -3.3] {
                let f = |t: f64| {
                    1.0 / ((1.0 - n * t.sin().powi(2)) * (1.0 - m * t.sin().powi(2)).sqrt())
                };
                let exact = phi.signum() * quad(f, phi.abs(), 4000);
                let v = elliptic_pi(n, phi, m);
                assert!(
                    (v - exact).abs() < 1e-11 * (1.0 + exact.abs()),
                    "n={n} m={m} phi={phi}: {v} vs {exact}"
                );
            }
        }
    }
}

#[test]
fn complete_third_kind_special_values() {
    // Pi(0, m) = K(m); Pi(m, m) = E(m)/(1-m)
    for &m in &[0.1_f64, 0.5, 0.9] {
        assert!((elliptic_pi_complete(0.0, m) - elliptic_k(m)).abs() < 1e-13);
        let e = rssn::kernels::special_ext::elliptic_e(m);
        assert!((elliptic_pi_complete(m, m) - e / (1.0 - m)).abs() < 1e-12);
    }
    // Pi(n, 0) = pi / (2 sqrt(1 - n))
    for &n in &[-2.0_f64, 0.3, 0.9] {
        let exact = std::f64::consts::FRAC_PI_2 / (1.0 - n).sqrt();
        assert!((elliptic_pi_complete(n, 0.0) - exact).abs() < 1e-13 * exact);
    }
}

#[test]
fn jacobi_integrals_match_quadrature() {
    for &m in &[0.2_f64, 0.7, 0.95] {
        for &n in &[-1.5_f64, 0.3, 0.9] {
            for &phi in &[0.7_f64, 2.0, 4.5, 8.0] {
                let g = |t: f64| {
                    1.0 / ((1.0 - n * t.sin().powi(2)).powi(2) * (1.0 - m * t.sin().powi(2)).sqrt())
                };
                let exact = quad(g, phi, 6000);
                let v = jacobi_int_inv_sq(n, phi, m);
                assert!(
                    (v - exact).abs() < 1e-10 * (1.0 + exact.abs()),
                    "n={n} m={m} phi={phi}: {v} vs {exact}"
                );
            }
        }
        for &phi in &[0.7_f64, 2.0, 4.5] {
            let h = |t: f64| t.sin().powi(2) / (1.0 - m * t.sin().powi(2)).sqrt();
            let exact = quad(h, phi, 4000);
            assert!((jacobi_int_sn2(phi, m) - exact).abs() < 1e-11);
        }
    }
}
