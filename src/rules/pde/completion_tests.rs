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

type Eigenpair = (f64, Box<dyn Fn(f64) -> f64>);

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
    let mut modes: Vec<Eigenpair> = Vec::new();
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

#[test]
fn robin_half_line_and_exterior_ball_heat() {
    let numeric = |text: &str, at: &[(&str, f64)]| crate::rules::testing::numeric(&crate::rules::standard(), text, at, 1e-9).0;
    let heat = "diff(u(x, t), t) = diff(diff(u(x, t), x), x)";
    // u_x = 2 u at x = 0: u = e^{2x + 4t} is an exact solution.
    let s = any(&format!("pdsolve({heat}, u(x, t), list(at(diff(u(x, t), x), x, 0) - 2*u(0, t) = 0, u(x, 0) = exp(2*x)))"));
    let (x, t) = (0.5_f64, 0.3_f64);
    close(numeric(&rhs(&s), &[("x", x), ("t", t)]), (2.0 * x + 4.0 * t).exp(), 1e-6, &s);
    // Decaying data against the Green's function with its Robin integral computed by quadrature.
    let s = any(&format!("pdsolve({heat}, u(x, t), list(at(diff(u(x, t), x), x, 0) + 3*u(0, t) = 0, u(x, 0) = exp(-x)))"));
    let kernel = |y: f64| (-y * y / (4.0 * t)).exp() / (4.0 * PI * t).sqrt();
    let kappa = -3.0_f64;
    let green = |s_: f64| {
        let robin = simpson(&|z| (-kappa * z).exp() * kernel(x + s_ + z), 0.0, 12.0);
        kernel(x - s_) + kernel(x + s_) - 2.0 * kappa * robin
    };
    // κ < 0 makes e^{-κ z} grow, so integrate the combination stably: κ = -3, z up to 12 is enough for e^{3z - z²/(4t)}.
    let reference = simpson(&|s_| green(s_) * (-s_).exp(), 0.0, 14.0);
    close(numeric(&rhs(&s), &[("x", x), ("t", t)]), reference, 1e-5, &s);
    // Neumann with drift on the half-line: constant data stays constant.
    let s = any("pdsolve(diff(u(x, t), t) = diff(diff(u(x, t), x), x) + 2*diff(u(x, t), x), u(x, t), list(at(diff(u(x, t), x), x, 0) = 0, u(x, 0) = 1))");
    close(numeric(&rhs(&s), &[("x", x), ("t", t)]), 1.0, 1e-6, &s);
    // Outside a ball with the Robin condition u_r + 3 u = 0 at r = 1: u = e^{-2(r-1) + 4t}/r.
    let ball = "diff(u(r, t), t) = diff(diff(u(r, t), r), r) + 2*diff(u(r, t), r)/r";
    let s = any(&format!("pdsolve({ball}, u(r, t), list(at(diff(u(r, t), r), r, 1) + 3*u(1, t) = 0, u(oo, t) = 0, u(r, 0) = exp(-2*(r - 1))/r))"));
    let r = 1.8_f64;
    close(numeric(&rhs(&s), &[("r", r), ("t", 0.2)]), (-2.0 * (r - 1.0) + 0.8).exp() / r, 1e-6, &s);
}

#[test]
fn exterior_heat_in_the_plane() {
    let disk = "diff(u(r, t), t) = diff(diff(u(r, t), r), r) + diff(u(r, t), r)/r";
    let data = "(r - 1)*exp(-2*(r - 1))";
    let s = any(&format!("pdsolve({disk}, u(r, t), list(u(1, t) = 0, u(oo, t) = 0, u(r, 0) = {data}))"));
    assert!(s.contains("bessely(0,") && s.contains("defint"), "{s}");
    // Crank-Nicolson on [1, 13] as the independent reference.
    let (n, steps) = (1200_usize, 250_usize);
    let (t_end, h) = (0.5_f64, 12.0 / 1200.0);
    let dt = t_end / steps as f64;
    let radius = |i: usize| 1.0 + h * i as f64;
    let mut u: Vec<f64> = (0..=n).map(|i| (radius(i) - 1.0) * (-2.0 * (radius(i) - 1.0)).exp()).collect();
    for _ in 0..steps {
        let (mut lower, mut diag, mut upper, mut rhs_) = (vec![0.0; n + 1], vec![1.0; n + 1], vec![0.0; n + 1], vec![0.0; n + 1]);
        for i in 1..n {
            let r_i = radius(i);
            let (a, b) = (1.0 / (h * h) - 1.0 / (2.0 * h * r_i), 1.0 / (h * h) + 1.0 / (2.0 * h * r_i));
            let operator = a * u[i - 1] - 2.0 / (h * h) * u[i] + b * u[i + 1];
            rhs_[i] = u[i] + 0.5 * dt * operator;
            lower[i] = -0.5 * dt * a;
            diag[i] = 1.0 + dt / (h * h);
            upper[i] = -0.5 * dt * b;
        }
        // Thomas algorithm.
        for i in 1..=n {
            let m = lower[i] / diag[i - 1];
            diag[i] -= m * upper[i - 1];
            rhs_[i] -= m * rhs_[i - 1];
        }
        u[n] = rhs_[n] / diag[n];
        for i in (0..n).rev() {
            u[i] = (rhs_[i] - upper[i] * u[i + 1]) / diag[i];
        }
    }
    let r = 2.0_f64;
    let index = ((r - 1.0) / h).round() as usize;
    let reference = u[index];
    let text = s.replace(", rho, 1, oo)", ", rho, 1, 9)").replace(", k, 0, oo)", ", k, 0, 6)");
    let text = expand_integrals(&rhs(&text), 90);
    close(value(&text, &[("r", r), ("t", t_end)]), reference, 6e-3, &s);
}

/// The argument list of the numeric Sturm–Liouville operators for the
/// coefficients `p`, `q`, `w` on `[a, b]` with the end conditions `ends`.
fn sl_arguments(
    a: f64,
    b: f64,
    ends: [f64; 4],
    p: &dyn Fn(f64) -> f64,
    q: &dyn Fn(f64) -> f64,
    w: &dyn Fn(f64) -> f64,
) -> Vec<f64> {
    let n = 24_usize;
    let mut args = vec![a, b, ends[0], ends[1], ends[2], ends[3], n as f64];
    let xs: Vec<f64> = (0..=n).map(|k| f64::midpoint(a, b) + 0.5 * (b - a) * (PI * k as f64 / n as f64).cos()).collect();
    for f in [p, q, w] {
        args.extend(xs.iter().map(|&x| f(x)));
    }
    args
}

#[test]
fn numeric_sturm_liouville_operators() {
    use super::sl::sl_eigen;
    use super::sl::sl_mode;
    use super::sl::sl_weight;
    let with = |m: f64, rest: &[f64]| {
        let mut v = vec![m];
        v.extend_from_slice(rest);
        v
    };
    // p = w = e^x on [0, π] with Dirichlet ends: μ_m = m² + 1/4, y_m = e^{-x/2} sin(m x)/√(π/2).
    let args = sl_arguments(0.0, PI, [1.0, 0.0, 1.0, 0.0], &f64::exp, &|_| 0.0, &f64::exp);
    for m in 1..=6 {
        let mu = sl_eigen(&with(f64::from(m), &args));
        close(mu, f64::from(m * m) + 0.25, 1e-7, "eigenvalue of e^x");
    }
    for x in [0.4_f64, 1.3, 2.9] {
        let mut mode_args = vec![3.0, x];
        mode_args.extend_from_slice(&args);
        close(sl_mode(&mode_args), (-x / 2.0).exp() * (3.0 * x).sin() / (PI / 2.0).sqrt(), 1e-8, "eigenfunction");
        let mut weight_args = vec![x];
        weight_args.extend_from_slice(&args);
        close(sl_weight(&weight_args), x.exp(), 1e-12, "weight");
    }
    // Euler type: p = x, w = 1/x on [1, 2]: μ_m = (m π / ln 2)².
    let args = sl_arguments(1.0, 2.0, [1.0, 0.0, 1.0, 0.0], &|x| x, &|_| 0.0, &|x| 1.0 / x);
    for m in 1..=4 {
        let want = (f64::from(m) * PI / 2.0_f64.ln()).powi(2);
        close(sl_eigen(&with(f64::from(m), &args)) / want, 1.0, 1e-8, "Euler eigenvalue");
    }
    // Robin ends with a negative eigenvalue: y'(0) = 2 y(0), y'(1) + 3 y(1) = 0 (p0 = -2, q0 = 1, p1 = 3, q1 = 1).
    let args = sl_arguments(0.0, 1.0, [5.0, 1.0, 3.0, 1.0], &|_| 1.0, &|_| 0.0, &|_| 1.0);
    let kappa = super::numeric::sl_neg_root(&[5.0, 1.0, 3.0, 1.0, 1.0]);
    assert!(kappa.is_finite());
    close(sl_eigen(&with(1.0, &args)), -kappa * kappa, 1e-8, "negative Robin eigenvalue");
    for m in 1..=4 {
        let k = super::numeric::sl_root(&[5.0, 1.0, 3.0, 1.0, 1.0, f64::from(m)]);
        close(sl_eigen(&with(f64::from(m + 1), &args)), k * k, 1e-7, "positive Robin eigenvalue");
    }
    // With a potential: p = w = 1, q = 3 + x on [0, 1], Dirichlet: compare with the Rayleigh quotient of the mode.
    let q = |x: f64| 3.0 + x;
    let args = sl_arguments(0.0, 1.0, [1.0, 0.0, 1.0, 0.0], &|_| 1.0, &q, &|_| 1.0);
    let mu = sl_eigen(&with(1.0, &args));
    let mode = |x: f64| {
        let mut v = vec![1.0, x];
        v.extend_from_slice(&args);
        sl_mode(&v)
    };
    let h = 1e-4;
    let hr = 1e-5;
    let rayleigh = simpson(&|x| ((mode(x + hr) - mode(x - hr)) / (2.0 * hr)).powi(2) + q(x) * mode(x).powi(2), 2e-5, 1.0 - 2e-5);
    close(simpson(&|x| mode(x).powi(2), 0.0, 1.0), 1.0, 1e-6, "normalisation");
    close(rayleigh, mu, 2e-3, "Rayleigh quotient");
    let residual = |x: f64| {
        let second = (mode(x + h) - 2.0 * mode(x) + mode(x - h)) / h.powi(2);
        q(x) * mode(x) - second - mu * mode(x)
    };
    close(residual(0.37), 0.0, 5e-3, "residual of the eigenfunction");
}

#[test]
fn separation_of_variables_for_variable_coefficients() {
    // u_t = (1 + x)² u_xx + (1 + x) u_x on [0, 1] (Dirichlet): s = ln(1 + x), modes sin(m π s / ln 2).
    let eq = "diff(u(x, t), t) = (1 + x)^2*diff(diff(u(x, t), x), x) + (1 + x)*diff(u(x, t), x)";
    let l = 2.0_f64.ln();
    let data = "sin(pi*ln(1 + x)/ln(2)) + sin(3*pi*ln(1 + x)/ln(2))/2";
    let s = any(&format!("pdsolve({eq}, u(x, t), list(u(0, t) = 0, u(1, t) = 0, u(x, 0) = {data}))"));
    assert!(s.contains("sl_mode") && s.contains("sl_eigen"), "{s}");
    let (x, t) = (0.45_f64, 0.01_f64);
    let s_x = x.ln_1p();
    let exact = (-(PI / l).powi(2) * t).exp() * (PI * s_x / l).sin() + 0.5 * (-(3.0 * PI / l).powi(2) * t).exp() * (3.0 * PI * s_x / l).sin();
    let text = expand_integrals(&expand_sums(&rhs(&s), 5), 160);
    close(value(&text, &[("x", x), ("t", t)]), exact, 2e-4, &s);
    // A potential and general data: u_t = u_xx + x u_x - u on [0, 1] with a Robin end, residual of the series.
    let eq = "diff(u(x, t), t) = diff(diff(u(x, t), x), x) + x*diff(u(x, t), x) - u(x, t)";
    let s = any(&format!("pdsolve({eq}, u(x, t), list(u(0, t) = 0, at(diff(u(x, t), x), x, 1) + 2*u(1, t) = 0, u(x, 0) = x*(1 - x)))"));
    assert!(s.contains("sl_mode"), "{s}");
    let text = expand_integrals(&expand_sums(&rhs(&s), 8), 120);
    let u = |x: f64, t: f64| value(&text, &[("x", x), ("t", t)]);
    let (x0, t0, h) = (0.4, 0.05, 1e-3);
    let residual = (u(x0, t0 + h) - u(x0, t0 - h)) / (2.0 * h) - (u(x0 + h, t0) - 2.0 * u(x0, t0) + u(x0 - h, t0)) / (h * h) - x0 * (u(x0 + h, t0) - u(x0 - h, t0)) / (2.0 * h) + u(x0, t0);
    close(residual, 0.0, 5e-3, "residual of the series");
    close(u(1.0, t0) * 2.0 + (u(1.0 + h, t0) - u(1.0 - h, t0)) / (2.0 * h), 0.0, 2e-2, "Robin end");
}
