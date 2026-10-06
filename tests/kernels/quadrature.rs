//! Tests for the advanced quadrature kernels.

use std::f64::consts::PI;

use rssn::kernels::quadrature::{
    GkRule, LowDiscrepancy, Sobol, clenshaw_curtis, cubature_nested, exp_sinh, filon,
    halton_point, integrate, integrate_gk, integrate_infinite, integrate_semi_infinite,
    monte_carlo, quasi_monte_carlo, romberg, sinh_sinh, tanh_sinh,
};

#[test]
fn gk_smooth_both_rules() {
    let exact = 1.0_f64.exp() - 1.0;
    for rule in [GkRule::G7K15, GkRule::G10K21] {
        let r = integrate_gk(f64::exp, 0.0, 1.0, rule, 1e-14, 1e-14, 100);
        assert!((r.value - exact).abs() < 1e-14);
        assert!(r.converged);
    }
}

#[test]
fn gk_polynomial_exact() {
    let r = integrate(|x| x.powi(13), 0.0, 1.0, 1e-12);
    assert!((r.value - 1.0 / 14.0).abs() < 1e-14);
}

#[test]
fn gk_endpoint_singularity() {
    let r = integrate_gk(|x| 1.0 / x.sqrt(), 0.0, 1.0, GkRule::G7K15, 1e-10, 1e-10, 2000);
    assert!((r.value - 2.0).abs() < 1e-7, "{r:?}");
    let r = integrate_gk(f64::ln, 0.0, 1.0, GkRule::G10K21, 1e-10, 1e-10, 2000);
    assert!((r.value + 1.0).abs() < 1e-8, "{r:?}");
}

#[test]
fn gk_infinite_ranges() {
    let r = integrate_semi_infinite(|x| (-x).exp(), 0.0, 1e-10);
    assert!((r.value - 1.0).abs() < 1e-9);
    let r = integrate_infinite(|x| (-x * x).exp(), 1e-10);
    assert!((r.value - PI.sqrt()).abs() < 1e-9);
}

#[test]
fn tanh_sinh_singular_and_smooth() {
    let r = tanh_sinh(|x| 1.0 / x.sqrt(), 0.0, 1.0, 1e-12);
    assert!((r.value - 2.0).abs() < 1e-11, "{r:?}");
    let r = tanh_sinh(f64::ln, 0.0, 1.0, 1e-12);
    assert!((r.value + 1.0).abs() < 1e-11, "{r:?}");
    let r = tanh_sinh(f64::cos, 0.0, PI / 2.0, 1e-12);
    assert!((r.value - 1.0).abs() < 1e-12);
    let r = tanh_sinh(|x| x, 1.0, 0.0, 1e-12);
    assert!((r.value + 0.5).abs() < 1e-12);
}

#[test]
fn exp_sinh_and_sinh_sinh() {
    let r = exp_sinh(|x| (-x).exp(), 0.0, 1e-12);
    assert!((r.value - 1.0).abs() < 1e-10, "{r:?}");
    let r = exp_sinh(|x| 1.0 / (1.0 + x * x), 0.0, 1e-12);
    assert!((r.value - PI / 2.0).abs() < 1e-9, "{r:?}");
    let r = sinh_sinh(|x| (-x * x).exp(), 1e-12);
    assert!((r.value - PI.sqrt()).abs() < 1e-11, "{r:?}");
}

#[test]
fn clenshaw_curtis_exactness() {
    let v = clenshaw_curtis(|x| x.powi(8) + 1.0, -1.0, 1.0, 8);
    assert!((v - (2.0 / 9.0 + 2.0)).abs() < 1e-14);
    let v = clenshaw_curtis(f64::exp, 0.0, 1.0, 24);
    assert!((v - (1.0_f64.exp() - 1.0)).abs() < 1e-14);
}

#[test]
fn romberg_reference() {
    let r = romberg(f64::sin, 0.0, PI, 1e-12, 14);
    assert!((r.value - 2.0).abs() < 1e-11);
    assert!(r.converged);
}

#[test]
fn filon_oscillatory() {
    let w: f64 = 50.0;
    let exact = w.sin() / w + (w.cos() - 1.0) / (w * w);
    let v = filon(|x| x, 0.0, 1.0, w, true, 10);
    assert!((v - exact).abs() < 1e-12, "{v} vs {exact}");
    let w: f64 = 100.0;
    let exact = w * (1.0 - PI.exp() * (w * PI).cos()) / (1.0 + w * w);
    let v = filon(f64::exp, 0.0, PI, w, false, 400);
    assert!((v - exact).abs() < 1e-8, "{v} vs {exact}");
    let v = filon(|x| x * x, 0.0, 1.0, 0.01, true, 4);
    assert!((v - 1.0 / 3.0).abs() < 1e-4);
}

#[test]
fn low_discrepancy_points() {
    assert_eq!(halton_point(1, 2), vec![0.5, 1.0 / 3.0]);
    assert_eq!(halton_point(2, 1), vec![0.25]);
    let mut s = Sobol::new(2);
    let p1 = s.next_point();
    let p2 = s.next_point();
    assert_eq!(p1, vec![0.5, 0.5]);
    assert_eq!(p2, vec![0.75, 0.25]);
}

#[test]
fn qmc_and_mc() {
    let f = |x: &[f64]| x.iter().map(|v| v * v).sum::<f64>();
    let lo = [0.0; 4];
    let hi = [1.0; 4];
    let exact = 4.0 / 3.0;
    let s = quasi_monte_carlo(f, &lo, &hi, 4096, LowDiscrepancy::Sobol);
    let h = quasi_monte_carlo(f, &lo, &hi, 4096, LowDiscrepancy::Halton);
    let m = monte_carlo(f, &lo, &hi, 4096, 42);
    assert!((s.value - exact).abs() < 1e-3, "{s:?}");
    assert!((h.value - exact).abs() < 1e-2, "{h:?}");
    assert!((m.value - exact).abs() < 5.0 * m.error + 1e-12);
    let m2 = monte_carlo(f, &lo, &hi, 4096, 42);
    assert_eq!(m.value.to_bits(), m2.value.to_bits());
}

#[test]
fn nested_cubature() {
    let v = cubature_nested(|x| (x[0] * x[1]).exp(), &[0.0, 0.0], &[1.0, 1.0], 1e-10);
    assert!((v.value - 1.317_902_151_454_404).abs() < 1e-8, "{v:?}");
    let v = cubature_nested(|x| x[0] + x[1] + x[2], &[0.0; 3], &[1.0; 3], 1e-10);
    assert!((v.value - 1.5).abs() < 1e-9);
}
