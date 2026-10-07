//! Tests of the completion layer: Struve functions and Fourier–Bessel
//! moments of odd powers, Robin problems with several non-positive modes,
//! modified Bessel modes of disks, Duhamel for massive waves, Euler–Cauchy
//! problems on `x < 0` and in several Euler variables.
//!
//! References are computed in the tests (own root scans and quadratures).

use super::dimension_tests::any;
use super::dimension_tests::close;
use super::dimension_tests::expand_integrals;
use super::dimension_tests::expand_sums;
use super::dimension_tests::rhs;
use super::dimension_tests::run;
use super::dimension_tests::series;
use super::dimension_tests::value;
use crate::kernels::special::bessel_i;
use crate::kernels::special::bessel_j;
use std::f64::consts::PI;

fn simpson(
    f: &dyn Fn(f64) -> f64,
    a: f64,
    b: f64,
) -> f64 {
    let n = 2000;
    let h = (b - a) / f64::from(n);
    let mut total = f(a) + f(b);
    for k in 1..n {
        total += f(a + f64::from(k) * h) * if k % 2 == 1 { 4.0 } else { 2.0 };
    }
    total * h / 3.0
}

/// The first `count` roots of `g` on `(0, top]` found by scanning and bisection.
fn roots(
    g: &dyn Fn(f64) -> f64,
    top: f64,
    count: usize,
) -> Vec<f64> {
    let step = 1e-3;
    let mut out = Vec::new();
    let mut x = 1e-4;
    while out.len() < count && x < top {
        let next = x + step;
        if g(x) * g(next) < 0.0 {
            let (mut lo, mut hi) = (x, next);
            for _ in 0..60 {
                let mid = f64::midpoint(lo, hi);
                if g(lo) * g(mid) <= 0.0 {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            out.push(f64::midpoint(lo, hi));
        }
        x = next;
    }
    out
}

#[test]
fn struve_functions() {
    let v = |src: &str| value(src, &[]);
    close(v("struve_h(0, 1)"), 0.568_656_627_048_287_9, 1e-12, "H0(1)");
    close(v("struve_h(1, 1)"), 0.198_457_336_201_943_8, 1e-12, "H1(1)");
    // d/dz [z^ν H_ν] = z^ν H_{ν-1} and d/dz [z^ν L_ν] = z^ν L_{ν-1}, across
    // the switch between the series and the integral representation.
    for x in [0.7_f64, 3.0, 7.95, 8.05, 12.0, 25.0] {
        let h = 1e-4;
        let f = |z: f64| z * super::numeric::struve_h(&[1.0, z]);
        let derivative = (f(x + h) - f(x - h)) / (2.0 * h);
        close(derivative, x * super::numeric::struve_h(&[0.0, x]), 1e-6, "(z H1)' = z H0");
        let g = |z: f64| z * super::numeric::struve_l(&[1.0, z]);
        let derivative = (g(x + h) - g(x - h)) / (2.0 * h);
        close(derivative / x.exp(), x * super::numeric::struve_l(&[0.0, x]) / x.exp(), 1e-6, "(z L1)' = z L0");
    }
    // H_{1/2}(x) = sqrt(2/(pi x)) (1 - cos x).
    for x in [0.5_f64, 3.0, 11.0] {
        close(super::numeric::struve_h(&[0.5, x]), (2.0 / (PI * x)).sqrt() * (1.0 - x.cos()), 1e-10, "H_1/2");
    }
}

#[test]
fn odd_power_fourier_bessel_moments_use_struve_functions() {
    let disk_heat = "diff(u(r, t), t) = diff(diff(u(r, t), r), r) + diff(u(r, t), r)/r";
    let s = any(&format!("pdsolve({disk_heat}, u(r, t), list(u(1, t) = 0, u(r, 0) = r))"));
    assert!(s.contains("struve_h") && !s.contains("defint"), "{s}");
    let (r, t) = (0.4, 0.05);
    let mut reference = 0.0;
    for m in 1..=30 {
        let j = value(&format!("bessel_zero(0, {m})"), &[]);
        let c = 2.0 * simpson(&|q| q * q * bessel_j(0.0, j * q), 0.0, 1.0) / bessel_j(1.0, j).powi(2);
        reference += c * (-j * j * t).exp() * bessel_j(0.0, j * r);
    }
    close(series(&s, &[("r", r), ("t", t)], 30), reference, 1e-8, &s);
    // Neumann data r^3 on a disk: three powers of the moment recursion.
    let s = any(&format!("pdsolve({disk_heat}, u(r, t), list(at(diff(u(r, t), r), r, 1) = 0, u(r, 0) = r^3))"));
    assert!(!s.contains("defint"), "{s}");
}

#[test]
fn robin_interval_with_two_negative_modes() {
    // X'(0) + 5 X(0) = 0 and X'(1) - 5 X(1) = 0: two negative eigenvalues.
    let s = any("pdsolve(diff(u(x, t), t) = diff(diff(u(x, t), x), x), u(x, t), list(at(diff(u(x, t), x), x, 0) + 5*u(0, t) = 0, at(diff(u(x, t), x), x, 1) - 5*u(1, t) = 0, u(x, 0) = x))");
    let (x, t) = (0.3, 0.02);
    // Reference: the modes cosh/sinh(κ(x - 1/2)), cos/sin(k(x - 1/2)) with own root scans.
    let mut modes: Vec<(f64, Box<dyn Fn(f64) -> f64>)> = Vec::new();
    for kappa in roots(&|k: f64| k * (k / 2.0).tanh() - 5.0, 20.0, 1) {
        modes.push((-kappa * kappa, Box::new(move |y: f64| (kappa * (y - 0.5)).cosh())));
    }
    for kappa in roots(&|k: f64| k - 5.0 * (k / 2.0).tanh(), 20.0, 1) {
        modes.push((-kappa * kappa, Box::new(move |y: f64| (kappa * (y - 0.5)).sinh())));
    }
    assert_eq!(modes.len(), 2, "two negative modes");
    for k in roots(&|k: f64| k * (k / 2.0).sin() + 5.0 * (k / 2.0).cos(), 130.0, 40) {
        modes.push((k * k, Box::new(move |y: f64| (k * (y - 0.5)).cos())));
    }
    for k in roots(&|k: f64| k * (k / 2.0).cos() - 5.0 * (k / 2.0).sin(), 130.0, 40) {
        modes.push((k * k, Box::new(move |y: f64| (k * (y - 0.5)).sin())));
    }
    let mut reference = 0.0;
    for (lambda, f) in &modes {
        let norm = simpson(&|y| f(y) * f(y), 0.0, 1.0);
        let c = simpson(&|y| y * f(y), 0.0, 1.0) / norm;
        reference += c * (-lambda * t).exp() * f(x);
    }
    let text = expand_integrals(&expand_sums(&rhs(&s), 60), 200);
    close(value(&text, &[("x", x), ("t", t)]), reference, 2e-4, &s);
}

#[test]
fn disk_with_a_negative_radial_mode() {
    let disk_heat = "diff(u(r, t), t) = diff(diff(u(r, t), r), r) + diff(u(r, t), r)/r";
    // u_r + h u = 0 at r = 1 with h = -3: one mode I_0(κ r), κ I_1(κ) = 3 I_0(κ).
    let s = any(&format!("pdsolve({disk_heat}, u(r, t), list(at(diff(u(r, t), r), r, 1) - 3*u(1, t) = 0, u(r, 0) = 1))"));
    assert!(s.contains("besseli"), "{s}");
    let (r, t) = (0.4, 0.02);
    let kappa = roots(&|k: f64| k * bessel_i(1.0, k) - 3.0 * bessel_i(0.0, k), 10.0, 1)[0];
    let mut reference = {
        let c = simpson(&|q| q * bessel_i(0.0, kappa * q), 0.0, 1.0) / simpson(&|q| q * bessel_i(0.0, kappa * q).powi(2), 0.0, 1.0);
        c * (kappa * kappa * t).exp() * bessel_i(0.0, kappa * r)
    };
    for k in roots(&|k: f64| -k * bessel_j(1.0, k) - 3.0 * bessel_j(0.0, k), 140.0, 40) {
        let c = simpson(&|q| q * bessel_j(0.0, k * q), 0.0, 1.0) / simpson(&|q| q * bessel_j(0.0, k * q).powi(2), 0.0, 1.0);
        reference += c * (-k * k * t).exp() * bessel_j(0.0, k * r);
    }
    let text = expand_integrals(&expand_sums(&rhs(&s), 40), 200);
    close(value(&text, &[("r", r), ("t", t)]), reference, 1e-5, &s);
}

#[test]
fn massive_wave_in_space_with_a_source() {
    let w3 = "diff(diff(u(x, y, z, t), t), t) = diff(diff(u(x, y, z, t), x), x) + diff(diff(u(x, y, z, t), y), y) + diff(diff(u(x, y, z, t), z), z)";
    // u_tt = Δu - u + 1 with zero data: u = 1 - cos t.
    let s = any(&format!("pdsolve({w3} - u(x, y, z, t) + 1, u(x, y, z, t), list(u(x, y, z, 0) = 0, at(diff(u(x, y, z, t), t), t, 0) = 0))"));
    let t = 0.9_f64;
    let text = expand_integrals(&rhs(&s), 160);
    close(value(&text, &[("x", 0.3), ("y", 0.5), ("z", 0.2), ("t", t)]), 1.0 - t.cos(), 2e-4, &s);
}

#[test]
fn euler_cauchy_on_the_negative_axis_and_in_two_variables() {
    let (x, t) = (-2.5_f64, 0.3_f64);
    let s = run("pdsolve(diff(u(x, t), t) = x^2*diff(diff(u(x, t), x), x) + x*diff(u(x, t), x), u(x, t), list(u(-1, t) = 0, u(-exp(pi), t) = 0, u(x, 0) = sin(ln(-x))))");
    close(value(&rhs(&s), &[("x", x), ("t", t)]), (-t).exp() * (-x).ln().sin(), 1e-12, &s);
    let eq = "diff(u(x, y, t), t) = x^2*diff(diff(u(x, y, t), x), x) + x*diff(u(x, y, t), x) + y^2*diff(diff(u(x, y, t), y), y) + y*diff(u(x, y, t), y)";
    let e = "exp(pi)";
    let s = run(&format!("pdsolve({eq}, u(x, y, t), list(u(1, y, t) = 0, u({e}, y, t) = 0, u(x, 1, t) = 0, u(x, {e}, t) = 0, u(x, y, 0) = sin(ln(x))*sin(ln(y))))"));
    close(value(&rhs(&s), &[("x", 2.5), ("y", 1.7), ("t", t)]), (-2.0 * t).exp() * 2.5_f64.ln().sin() * 1.7_f64.ln().sin(), 1e-12, &s);
}

/// The right-hand sides `(u, v)` of `list(u(x, t) = a, v(x, t) = b)`.
fn pair(text: &str) -> (String, String) {
    let inner = text.trim_start_matches("list(").strip_suffix(')').unwrap_or("");
    let (a, b) = inner.split_once(", v(x, t) = ").unwrap_or((inner, ""));
    (a.trim_start_matches("u(x, t) = ").to_owned(), b.to_owned())
}

#[test]
fn jordan_systems_with_sources_and_boundaries() {
    let numeric = |text: &str, x: f64, t: f64| crate::rules::testing::numeric(&crate::rules::standard(), text, &[("x", x), ("t", t)], 1e-9).0;
    // A source in the second chain element, whole line: v_t + v_x = cos x, u_t + u_x + v_x = 0.
    let s = any("pdsolve(list(diff(u(x, t), t) + diff(u(x, t), x) + diff(v(x, t), x) = 0, diff(v(x, t), t) + diff(v(x, t), x) = cos(x)), list(u(x, t), v(x, t)), list(u(x, 0) = 0, v(x, 0) = 0))");
    let (a, b) = pair(&s);
    let (x, t) = (0.7_f64, 0.4_f64);
    close(numeric(&b, x, t), x.sin() - (x - t).sin(), 1e-8, &s);
    close(numeric(&a, x, t), -(x.sin() - (x - t).sin()) + t * (x - t).cos(), 1e-8, &s);
    // A Jordan block with a derivative-free coupling on a bounded interval:
    // u_t = u_xx + v, v_t = v_xx with Dirichlet ends: u = t e^{-t} sin x for v(0) = sin x.
    let s = any("pdsolve(list(diff(u(x, t), t) = diff(diff(u(x, t), x), x) + v(x, t), diff(v(x, t), t) = diff(diff(v(x, t), x), x)), list(u(x, t), v(x, t)), list(u(x, 0) = 0, v(x, 0) = sin(x), u(0, t) = 0, u(pi, t) = 0, v(0, t) = 0, v(pi, t) = 0))");
    let (a, b) = pair(&s);
    close(numeric(&b, x, t), (-t).exp() * x.sin(), 1e-9, &s);
    close(numeric(&a, x, t), t * (-t).exp() * x.sin(), 1e-9, &s);
}
