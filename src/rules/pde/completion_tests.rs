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
use crate::graph::Engine;
use crate::graph::Env;
use crate::graph::Graph;
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
        close(
            derivative,
            x * super::numeric::struve_h(&[0.0, x]),
            1e-6,
            "(z H1)' = z H0",
        );
        let g = |z: f64| z * super::numeric::struve_l(&[1.0, z]);
        let derivative = (g(x + h) - g(x - h)) / (2.0 * h);
        close(
            derivative / x.exp(),
            x * super::numeric::struve_l(&[0.0, x]) / x.exp(),
            1e-6,
            "(z L1)' = z L0",
        );
    }
    // H_{1/2}(x) = sqrt(2/(pi x)) (1 - cos x).
    for x in [0.5_f64, 3.0, 11.0] {
        close(
            super::numeric::struve_h(&[0.5, x]),
            (2.0 / (PI * x)).sqrt() * (1.0 - x.cos()),
            1e-10,
            "H_1/2",
        );
    }
}

#[test]
fn odd_power_fourier_bessel_moments_use_struve_functions() {
    let disk_heat = "diff(u(r, t), t) = diff(diff(u(r, t), r), r) + diff(u(r, t), r)/r";
    let s = any(&format!(
        "pdsolve({disk_heat}, u(r, t), list(u(1, t) = 0, u(r, 0) = r))"
    ));
    assert!(s.contains("struve_h") && !s.contains("defint"), "{s}");
    let (r, t) = (0.4, 0.05);
    let mut reference = 0.0;
    for m in 1..=30 {
        let j = value(&format!("bessel_zero(0, {m})"), &[]);
        let c =
            2.0 * simpson(&|q| q * q * bessel_j(0.0, j * q), 0.0, 1.0) / bessel_j(1.0, j).powi(2);
        reference += c * (-j * j * t).exp() * bessel_j(0.0, j * r);
    }
    close(series(&s, &[("r", r), ("t", t)], 30), reference, 1e-8, &s);
    // Neumann data r^3 on a disk: three powers of the moment recursion.
    let s = any(&format!(
        "pdsolve({disk_heat}, u(r, t), list(at(diff(u(r, t), r), r, 1) = 0, u(r, 0) = r^3))"
    ));
    assert!(!s.contains("defint"), "{s}");
}

#[test]
fn robin_interval_with_two_negative_modes() {
    // X'(0) + 5 X(0) = 0 and X'(1) - 5 X(1) = 0: two negative eigenvalues.
    let s = any(
        "pdsolve(diff(u(x, t), t) = diff(diff(u(x, t), x), x), u(x, t), list(at(diff(u(x, t), x), x, 0) + 5*u(0, t) = 0, at(diff(u(x, t), x), x, 1) - 5*u(1, t) = 0, u(x, 0) = x))",
    );
    let (x, t) = (0.3, 0.02);
    // Reference: the modes cosh/sinh(κ(x - 1/2)), cos/sin(k(x - 1/2)) with own root scans.
    let mut modes: Vec<Eigenpair> = Vec::new();
    for kappa in roots(&|k: f64| k * (k / 2.0).tanh() - 5.0, 20.0, 1) {
        modes.push((
            -kappa * kappa,
            Box::new(move |y: f64| (kappa * (y - 0.5)).cosh()),
        ));
    }
    for kappa in roots(&|k: f64| k - 5.0 * (k / 2.0).tanh(), 20.0, 1) {
        modes.push((
            -kappa * kappa,
            Box::new(move |y: f64| (kappa * (y - 0.5)).sinh()),
        ));
    }
    assert_eq!(modes.len(), 2, "two negative modes");
    for k in roots(
        &|k: f64| k * (k / 2.0).sin() + 5.0 * (k / 2.0).cos(),
        130.0,
        40,
    ) {
        modes.push((k * k, Box::new(move |y: f64| (k * (y - 0.5)).cos())));
    }
    for k in roots(
        &|k: f64| k * (k / 2.0).cos() - 5.0 * (k / 2.0).sin(),
        130.0,
        40,
    ) {
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
    let s = any(&format!(
        "pdsolve({disk_heat}, u(r, t), list(at(diff(u(r, t), r), r, 1) - 3*u(1, t) = 0, u(r, 0) = 1))"
    ));
    assert!(s.contains("besseli"), "{s}");
    let (r, t) = (0.4, 0.02);
    let kappa = roots(
        &|k: f64| k * bessel_i(1.0, k) - 3.0 * bessel_i(0.0, k),
        10.0,
        1,
    )[0];
    let mut reference = {
        let c = simpson(&|q| q * bessel_i(0.0, kappa * q), 0.0, 1.0)
            / simpson(&|q| q * bessel_i(0.0, kappa * q).powi(2), 0.0, 1.0);
        c * (kappa * kappa * t).exp() * bessel_i(0.0, kappa * r)
    };
    for k in roots(
        &|k: f64| -k * bessel_j(1.0, k) - 3.0 * bessel_j(0.0, k),
        140.0,
        40,
    ) {
        let c = simpson(&|q| q * bessel_j(0.0, k * q), 0.0, 1.0)
            / simpson(&|q| q * bessel_j(0.0, k * q).powi(2), 0.0, 1.0);
        reference += c * (-k * k * t).exp() * bessel_j(0.0, k * r);
    }
    let text = expand_integrals(&expand_sums(&rhs(&s), 40), 200);
    close(value(&text, &[("r", r), ("t", t)]), reference, 1e-5, &s);
}

#[test]
fn massive_wave_in_space_with_a_source() {
    let w3 = "diff(diff(u(x, y, z, t), t), t) = diff(diff(u(x, y, z, t), x), x) + diff(diff(u(x, y, z, t), y), y) + diff(diff(u(x, y, z, t), z), z)";
    // u_tt = Δu - u + 1 with zero data: u = 1 - cos t.
    let s = any(&format!(
        "pdsolve({w3} - u(x, y, z, t) + 1, u(x, y, z, t), list(u(x, y, z, 0) = 0, at(diff(u(x, y, z, t), t), t, 0) = 0))"
    ));
    let t = 0.9_f64;
    let text = expand_integrals(&rhs(&s), 160);
    close(
        value(&text, &[("x", 0.3), ("y", 0.5), ("z", 0.2), ("t", t)]),
        1.0 - t.cos(),
        2e-4,
        &s,
    );
}

#[test]
fn euler_cauchy_on_the_negative_axis_and_in_two_variables() {
    let (x, t) = (-2.5_f64, 0.3_f64);
    let s = run(
        "pdsolve(diff(u(x, t), t) = x^2*diff(diff(u(x, t), x), x) + x*diff(u(x, t), x), u(x, t), list(u(-1, t) = 0, u(-exp(pi), t) = 0, u(x, 0) = sin(ln(-x))))",
    );
    close(
        value(&rhs(&s), &[("x", x), ("t", t)]),
        (-t).exp() * (-x).ln().sin(),
        1e-12,
        &s,
    );
    let eq = "diff(u(x, y, t), t) = x^2*diff(diff(u(x, y, t), x), x) + x*diff(u(x, y, t), x) + y^2*diff(diff(u(x, y, t), y), y) + y*diff(u(x, y, t), y)";
    let e = "exp(pi)";
    let s = run(&format!(
        "pdsolve({eq}, u(x, y, t), list(u(1, y, t) = 0, u({e}, y, t) = 0, u(x, 1, t) = 0, u(x, {e}, t) = 0, u(x, y, 0) = sin(ln(x))*sin(ln(y))))"
    ));
    close(
        value(&rhs(&s), &[("x", 2.5), ("y", 1.7), ("t", t)]),
        (-2.0 * t).exp() * 2.5_f64.ln().sin() * 1.7_f64.ln().sin(),
        1e-12,
        &s,
    );
}

/// The right-hand sides `(u, v)` of `list(u(x, t) = a, v(x, t) = b)`.
fn pair(text: &str) -> (String, String) {
    let inner = text
        .trim_start_matches("list(")
        .strip_suffix(')')
        .unwrap_or("");
    let (a, b) = inner.split_once(", v(x, t) = ").unwrap_or((inner, ""));
    (a.trim_start_matches("u(x, t) = ").to_owned(), b.to_owned())
}

#[test]
fn jordan_systems_with_sources_and_boundaries() {
    let numeric = |text: &str, x: f64, t: f64| {
        crate::rules::testing::numeric(&crate::rules::standard(), text, &[("x", x), ("t", t)], 1e-9)
            .0
    };
    // A source in the second chain element, whole line: v_t + v_x = cos x, u_t + u_x + v_x = 0.
    let s = any(
        "pdsolve(list(diff(u(x, t), t) + diff(u(x, t), x) + diff(v(x, t), x) = 0, diff(v(x, t), t) + diff(v(x, t), x) = cos(x)), list(u(x, t), v(x, t)), list(u(x, 0) = 0, v(x, 0) = 0))",
    );
    let (a, b) = pair(&s);
    let (x, t) = (0.7_f64, 0.4_f64);
    close(numeric(&b, x, t), x.sin() - (x - t).sin(), 1e-8, &s);
    close(
        numeric(&a, x, t),
        -(x.sin() - (x - t).sin()) + t * (x - t).cos(),
        1e-8,
        &s,
    );
    // A Jordan block with a derivative-free coupling on a bounded interval:
    // u_t = u_xx + v, v_t = v_xx with Dirichlet ends: u = t e^{-t} sin x for v(0) = sin x.
    let s = any(
        "pdsolve(list(diff(u(x, t), t) = diff(diff(u(x, t), x), x) + v(x, t), diff(v(x, t), t) = diff(diff(v(x, t), x), x)), list(u(x, t), v(x, t)), list(u(x, 0) = 0, v(x, 0) = sin(x), u(0, t) = 0, u(pi, t) = 0, v(0, t) = 0, v(pi, t) = 0))",
    );
    let (a, b) = pair(&s);
    close(numeric(&b, x, t), (-t).exp() * x.sin(), 1e-9, &s);
    close(numeric(&a, x, t), t * (-t).exp() * x.sin(), 1e-9, &s);
}

#[test]
fn robin_half_line_and_exterior_ball_heat() {
    let numeric = |text: &str, at: &[(&str, f64)]| {
        crate::rules::testing::numeric(&crate::rules::standard(), text, at, 1e-9).0
    };
    let heat = "diff(u(x, t), t) = diff(diff(u(x, t), x), x)";
    // u_x = 2 u at x = 0: u = e^{2x + 4t} is an exact solution.
    let s = any(&format!(
        "pdsolve({heat}, u(x, t), list(at(diff(u(x, t), x), x, 0) - 2*u(0, t) = 0, u(x, 0) = exp(2*x)))"
    ));
    let (x, t) = (0.5_f64, 0.3_f64);
    close(
        numeric(&rhs(&s), &[("x", x), ("t", t)]),
        (2.0 * x + 4.0 * t).exp(),
        1e-6,
        &s,
    );
    // Decaying data against the Green's function with its Robin integral computed by quadrature.
    let s = any(&format!(
        "pdsolve({heat}, u(x, t), list(at(diff(u(x, t), x), x, 0) + 3*u(0, t) = 0, u(x, 0) = exp(-x)))"
    ));
    let kernel = |y: f64| (-y * y / (4.0 * t)).exp() / (4.0 * PI * t).sqrt();
    let kappa = -3.0_f64;
    let green = |s_: f64| {
        let robin = simpson(&|z| (-kappa * z).exp() * kernel(x + s_ + z), 0.0, 12.0);
        kernel(x - s_) + kernel(x + s_) - 2.0 * kappa * robin
    };
    // κ < 0 makes e^{-κ z} grow, so integrate the combination stably: κ = -3, z up to 12 is enough for e^{3z - z²/(4t)}.
    let reference = simpson(&|s_| green(s_) * (-s_).exp(), 0.0, 14.0);
    close(
        numeric(&rhs(&s), &[("x", x), ("t", t)]),
        reference,
        1e-5,
        &s,
    );
    // Neumann with drift on the half-line: constant data stays constant.
    let s = any(
        "pdsolve(diff(u(x, t), t) = diff(diff(u(x, t), x), x) + 2*diff(u(x, t), x), u(x, t), list(at(diff(u(x, t), x), x, 0) = 0, u(x, 0) = 1))",
    );
    close(numeric(&rhs(&s), &[("x", x), ("t", t)]), 1.0, 1e-6, &s);
    // Outside a ball with the Robin condition u_r + 3 u = 0 at r = 1: u = e^{-2(r-1) + 4t}/r.
    let ball = "diff(u(r, t), t) = diff(diff(u(r, t), r), r) + 2*diff(u(r, t), r)/r";
    let s = any(&format!(
        "pdsolve({ball}, u(r, t), list(at(diff(u(r, t), r), r, 1) + 3*u(1, t) = 0, u(oo, t) = 0, u(r, 0) = exp(-2*(r - 1))/r))"
    ));
    let r = 1.8_f64;
    close(
        numeric(&rhs(&s), &[("r", r), ("t", 0.2)]),
        (-2.0 * (r - 1.0) + 0.8).exp() / r,
        1e-6,
        &s,
    );
}

#[test]
fn exterior_heat_in_the_plane() {
    let disk = "diff(u(r, t), t) = diff(diff(u(r, t), r), r) + diff(u(r, t), r)/r";
    let data = "(r - 1)*exp(-2*(r - 1))";
    let s = any(&format!(
        "pdsolve({disk}, u(r, t), list(u(1, t) = 0, u(oo, t) = 0, u(r, 0) = {data}))"
    ));
    assert!(s.contains("bessely(0,") && s.contains("defint"), "{s}");
    // Crank-Nicolson on [1, 13] as the independent reference.
    let (n, steps) = (1200_usize, 250_usize);
    let (t_end, h) = (0.5_f64, 12.0 / 1200.0);
    let dt = t_end / steps as f64;
    let radius = |i: usize| 1.0 + h * i as f64;
    let mut u: Vec<f64> = (0..=n)
        .map(|i| (radius(i) - 1.0) * (-2.0 * (radius(i) - 1.0)).exp())
        .collect();
    for _ in 0..steps {
        let (mut lower, mut diag, mut upper, mut rhs_) = (
            vec![0.0; n + 1],
            vec![1.0; n + 1],
            vec![0.0; n + 1],
            vec![0.0; n + 1],
        );
        for i in 1..n {
            let r_i = radius(i);
            let (a, b) = (
                1.0 / (h * h) - 1.0 / (2.0 * h * r_i),
                1.0 / (h * h) + 1.0 / (2.0 * h * r_i),
            );
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
    let text = s
        .replace(", rho, 1, oo)", ", rho, 1, 9)")
        .replace(", k, 0, oo)", ", k, 0, 6)");
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
    let xs: Vec<f64> = (0..=n)
        .map(|k| f64::midpoint(a, b) + 0.5 * (b - a) * (PI * k as f64 / n as f64).cos())
        .collect();
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
    let args = sl_arguments(
        0.0,
        PI,
        [1.0, 0.0, 1.0, 0.0],
        &f64::exp,
        &|_| 0.0,
        &f64::exp,
    );
    for m in 1..=6 {
        let mu = sl_eigen(&with(f64::from(m), &args));
        close(mu, f64::from(m * m) + 0.25, 1e-7, "eigenvalue of e^x");
    }
    for x in [0.4_f64, 1.3, 2.9] {
        let mut mode_args = vec![3.0, x];
        mode_args.extend_from_slice(&args);
        close(
            sl_mode(&mode_args),
            (-x / 2.0).exp() * (3.0 * x).sin() / (PI / 2.0).sqrt(),
            1e-8,
            "eigenfunction",
        );
        let mut weight_args = vec![x];
        weight_args.extend_from_slice(&args);
        close(sl_weight(&weight_args), x.exp(), 1e-12, "weight");
    }
    // Euler type: p = x, w = 1/x on [1, 2]: μ_m = (m π / ln 2)².
    let args = sl_arguments(1.0, 2.0, [1.0, 0.0, 1.0, 0.0], &|x| x, &|_| 0.0, &|x| {
        1.0 / x
    });
    for m in 1..=4 {
        let want = (f64::from(m) * PI / 2.0_f64.ln()).powi(2);
        close(
            sl_eigen(&with(f64::from(m), &args)) / want,
            1.0,
            1e-8,
            "Euler eigenvalue",
        );
    }
    // Robin ends with a negative eigenvalue: y'(0) = 2 y(0), y'(1) + 3 y(1) = 0 (p0 = -2, q0 = 1, p1 = 3, q1 = 1).
    let args = sl_arguments(0.0, 1.0, [5.0, 1.0, 3.0, 1.0], &|_| 1.0, &|_| 0.0, &|_| 1.0);
    let kappa = super::numeric::sl_neg_root(&[5.0, 1.0, 3.0, 1.0, 1.0]);
    assert!(kappa.is_finite());
    close(
        sl_eigen(&with(1.0, &args)),
        -kappa * kappa,
        1e-8,
        "negative Robin eigenvalue",
    );
    for m in 1..=4 {
        let k = super::numeric::sl_root(&[5.0, 1.0, 3.0, 1.0, 1.0, f64::from(m)]);
        close(
            sl_eigen(&with(f64::from(m + 1), &args)),
            k * k,
            1e-7,
            "positive Robin eigenvalue",
        );
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
    let rayleigh = simpson(
        &|x| ((mode(x + hr) - mode(x - hr)) / (2.0 * hr)).powi(2) + q(x) * mode(x).powi(2),
        2e-5,
        1.0 - 2e-5,
    );
    close(
        simpson(&|x| mode(x).powi(2), 0.0, 1.0),
        1.0,
        1e-6,
        "normalisation",
    );
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
    let s = any(&format!(
        "pdsolve({eq}, u(x, t), list(u(0, t) = 0, u(1, t) = 0, u(x, 0) = {data}))"
    ));
    assert!(s.contains("sl_mode") && s.contains("sl_eigen"), "{s}");
    let (x, t) = (0.45_f64, 0.01_f64);
    let s_x = x.ln_1p();
    let exact = (-(PI / l).powi(2) * t).exp() * (PI * s_x / l).sin()
        + 0.5 * (-(3.0 * PI / l).powi(2) * t).exp() * (3.0 * PI * s_x / l).sin();
    let text = expand_integrals(&expand_sums(&rhs(&s), 5), 160);
    close(value(&text, &[("x", x), ("t", t)]), exact, 2e-4, &s);
    // A potential and general data: u_t = u_xx + x u_x - u on [0, 1] with a Robin end, residual of the series.
    let eq = "diff(u(x, t), t) = diff(diff(u(x, t), x), x) + x*diff(u(x, t), x) - u(x, t)";
    let s = any(&format!(
        "pdsolve({eq}, u(x, t), list(u(0, t) = 0, at(diff(u(x, t), x), x, 1) + 2*u(1, t) = 0, u(x, 0) = x*(1 - x)))"
    ));
    assert!(s.contains("sl_mode"), "{s}");
    let text = expand_integrals(&expand_sums(&rhs(&s), 8), 120);
    let u = |x: f64, t: f64| value(&text, &[("x", x), ("t", t)]);
    let (x0, t0, h) = (0.4, 0.05, 1e-3);
    let residual = (u(x0, t0 + h) - u(x0, t0 - h)) / (2.0 * h)
        - (u(x0 + h, t0) - 2.0 * u(x0, t0) + u(x0 - h, t0)) / (h * h)
        - x0 * (u(x0 + h, t0) - u(x0 - h, t0)) / (2.0 * h)
        + u(x0, t0);
    close(residual, 0.0, 5e-3, "residual of the series");
    close(
        u(1.0, t0) * 2.0 + (u(1.0 + h, t0) - u(1.0 - h, t0)) / (2.0 * h),
        0.0,
        2e-2,
        "Robin end",
    );
}

#[test]
fn classical_orthogonal_polynomial_operators() {
    let (x, t) = (0.7_f64, 0.4_f64);
    let numeric = |text: &str| {
        crate::rules::testing::numeric(&crate::rules::standard(), text, &[("x", x), ("t", t)], 1e-9)
            .0
    };
    // Hermite type (Ornstein-Uhlenbeck): u_t = u_xx - x u_x, x² = He_2 + 1.
    let s = any(
        "pdsolve(diff(u(x, t), t) = diff(diff(u(x, t), x), x) - x*diff(u(x, t), x), u(x, t), list(u(x, 0) = x^2))",
    );
    close(
        numeric(&rhs(&s)),
        1.0 + (x * x - 1.0) * (-2.0 * t).exp(),
        1e-8,
        &s,
    );
    // Laguerre type: u_t = x u_xx + (1 - x) u_x, x² = y₂ + 4 y₁ + 2.
    let s = any(
        "pdsolve(diff(u(x, t), t) = x*diff(diff(u(x, t), x), x) + (1 - x)*diff(u(x, t), x), u(x, t), list(u(x, 0) = x^2))",
    );
    close(
        numeric(&rhs(&s)),
        2.0 + 4.0 * (-t).exp() * (x - 1.0) + (-2.0 * t).exp() * (x * x - 4.0 * x + 2.0),
        1e-8,
        &s,
    );
    // Jacobi type on [0, 1]: u_t = x(1 - x) u_xx + (1/2 - 3x/2) u_x, x = 1/3 + y₁.
    let s = any(
        "pdsolve(diff(u(x, t), t) = x*(1 - x)*diff(diff(u(x, t), x), x) + (1/2 - 3*x/2)*diff(u(x, t), x), u(x, t), list(u(x, 0) = x))",
    );
    close(
        numeric(&rhs(&s)),
        1.0 / 3.0 + (x - 1.0 / 3.0) * (-1.5 * t).exp(),
        1e-8,
        &s,
    );
}

/// `v_tt = v_xx` on `x > 0` with `v_x + ρ v = 0` at 0, `v(x, 0) = f`, `v_t = 0`, by leapfrog.
fn robin_wave_reference(
    rho: f64,
    f: &dyn Fn(f64) -> f64,
    at: f64,
    t_end: f64,
) -> f64 {
    let (dx, dt) = (0.004_f64, 0.002_f64);
    let n = 4500_usize;
    let steps = (t_end / dt).round() as usize;
    let x = |i: usize| dx * i as f64;
    let mut previous: Vec<f64> = (0..=n).map(|i| f(x(i))).collect();
    let laplacian = |u: &[f64]| -> Vec<f64> {
        let mut out = vec![0.0; u.len()];
        for i in 0..u.len() {
            let left = if i == 0 {
                u[1] + 2.0 * dx * rho * u[0]
            } else {
                u[i - 1]
            };
            let right = if i + 1 == u.len() {
                0.0
            } else {
                u[i + 1]
            };
            out[i] = (left - 2.0 * u[i] + right) / (dx * dx);
        }
        out
    };
    let l0 = laplacian(&previous);
    let mut current: Vec<f64> = previous
        .iter()
        .zip(&l0)
        .map(|(u, l)| u + 0.5 * dt * dt * l)
        .collect();
    for _ in 1..steps {
        let l = laplacian(&current);
        let next: Vec<f64> = (0..=n)
            .map(|i| 2.0 * current[i] - previous[i] + dt * dt * l[i])
            .collect();
        previous = current;
        current = next;
    }
    current[(at / dx).round() as usize]
}

#[test]
fn robin_and_neumann_waves_on_half_lines_and_outside_balls() {
    let bump = |x: f64| (-4.0 * (x - 3.0) * (x - 3.0)).exp();
    let t = 1.2_f64;
    let wave = "diff(diff(u(x, t), t), t) = diff(diff(u(x, t), x), x)";
    let s = any(&format!(
        "pdsolve({wave}, u(x, t), list(at(diff(u(x, t), x), x, 0) + 7/10*u(0, t) = 0, u(x, 0) = exp(-4*(x - 3)^2), at(diff(u(x, t), t), t, 0) = 0))"
    ));
    assert!(s.contains("defint"), "{s}");
    let text = expand_integrals(&rhs(&s).replace(", oo)", ", 12)"), 500);
    close(
        value(&text, &[("x", 2.0), ("t", t)]),
        robin_wave_reference(0.7, &bump, 2.0, t),
        4e-3,
        &s,
    );
    // Outside the unit ball with a Neumann condition: v = r u, v_x - v = 0 (ρ = -1/R).
    let ball = "diff(diff(u(r, t), t), t) = diff(diff(u(r, t), r), r) + 2*diff(u(r, t), r)/r";
    let s = any(&format!(
        "pdsolve({ball}, u(r, t), list(at(diff(u(r, t), r), r, 1) = 0, u(oo, t) = 0, u(r, 0) = exp(-4*(r - 4)^2)/r, at(diff(u(r, t), t), t, 0) = 0))"
    ));
    let text = expand_integrals(&rhs(&s).replace(", oo)", ", 13)"), 500);
    close(
        value(&text, &[("r", 3.0), ("t", t)]),
        robin_wave_reference(-1.0, &bump, 2.0, t) / 3.0,
        4e-3,
        &s,
    );
}

#[test]
fn annulus_with_a_negative_radial_mode() {
    let heat = "diff(u(r, t), t) = diff(diff(u(r, t), r), r) + diff(u(r, t), r)/r";
    // u_r + 5u = 0 at r = 1 (wrong sign for the inner circle), u_r + u = 0 at r = 2.
    let s = any(&format!(
        "pdsolve({heat}, u(r, t), list(at(diff(u(r, t), r), r, 1) + 5*u(1, t) = 0, at(diff(u(r, t), r), r, 2) + u(2, t) = 0, u(r, 0) = 1))"
    ));
    assert!(s.contains("annulus_neg_root"), "{s}");
    let (r, t) = (1.4_f64, 0.02_f64);
    let (pa, pb) = (5.0_f64, 1.0_f64);
    // J_0 / Y_0 and I_0 / K_0 conditions p F + F' at the two radii.
    let cond_jy = |second: bool, p: f64, k: f64, radius: f64| {
        let z = k * radius;
        if second {
            p * crate::kernels::special::bessel_y0(z) - k * crate::kernels::special::bessel_y1(z)
        } else {
            p * crate::kernels::special::bessel_j0(z) - k * crate::kernels::special::bessel_j1(z)
        }
    };
    let cond_ik = |second: bool, p: f64, k: f64, radius: f64| {
        let z = k * radius;
        if second {
            p * crate::kernels::special::bessel_k0(z) - k * crate::kernels::special::bessel_k1(z)
        } else {
            p * crate::kernels::special::bessel_i0(z) + k * crate::kernels::special::bessel_i1(z)
        }
    };
    // Negative mode: F = cond_K(a) I_0 - cond_I(a) K_0 with cond_K(a) cond_I(b) = cond_I(a) cond_K(b).
    let g_neg = |k: f64| {
        cond_ik(true, pa, k, 1.0) * cond_ik(false, pb, k, 2.0)
            - cond_ik(false, pa, k, 1.0) * cond_ik(true, pb, k, 2.0)
    };
    let kappas = roots(&g_neg, 12.0, 1);
    assert_eq!(kappas.len(), 1, "one negative mode");
    let kappa = kappas[0];
    let neg_mode = move |q: f64| {
        cond_ik(true, pa, kappa, 1.0) * crate::kernels::special::bessel_i0(kappa * q)
            - cond_ik(false, pa, kappa, 1.0) * crate::kernels::special::bessel_k0(kappa * q)
    };
    let coefficient = |mode: &dyn Fn(f64) -> f64| {
        simpson(&|q| q * mode(q), 1.0, 2.0) / simpson(&|q| q * mode(q) * mode(q), 1.0, 2.0)
    };
    let mut reference = coefficient(&neg_mode) * (kappa * kappa * t).exp() * neg_mode(r);
    let g_pos = |k: f64| {
        cond_jy(true, pa, k, 1.0) * cond_jy(false, pb, k, 2.0)
            - cond_jy(false, pa, k, 1.0) * cond_jy(true, pb, k, 2.0)
    };
    for k in roots(&g_pos, 130.0, 40) {
        let mode = move |q: f64| {
            cond_jy(true, pa, k, 1.0) * crate::kernels::special::bessel_j0(k * q)
                - cond_jy(false, pa, k, 1.0) * crate::kernels::special::bessel_y0(k * q)
        };
        reference += coefficient(&mode) * (-k * k * t).exp() * mode(r);
    }
    let text = expand_integrals(&expand_sums(&rhs(&s), 40), 160);
    close(value(&text, &[("r", r), ("t", t)]), reference, 2e-4, &s);
}

#[test]
fn fourier_symbol_fallback_on_the_line() {
    use super::extension_tests::complex_value;
    let (x, t) = (0.8_f64, 0.3_f64);
    let bounds = |text: &str| text.replace("-oo", "-8").replace(", oo)", ", 8)");
    // Beam: u_tt + u_xxxx = 0, u(x, 0) = e^{-x²}, u_t(x, 0) = 0:  û = û₀ cos(k² t).
    let s = any(
        "pdsolve(diff(diff(u(x, t), t), t) + diff(diff(diff(diff(u(x, t), x), x), x), x) = 0, u(x, t), list(u(x, 0) = exp(-x^2), at(diff(u(x, t), t), t, 0) = 0))",
    );
    assert!(s.contains("defint"), "{s}");
    let g_hat = |k: f64| PI.sqrt() * (-k * k / 4.0).exp();
    let reference = simpson(
        &|k| g_hat(k) * (k * k * t).cos() * (k * x).cos(),
        -12.0,
        12.0,
    ) / (2.0 * PI);
    let text = expand_integrals(&bounds(&rhs(&s)), 110);
    let got = complex_value(&text, &[("x", x), ("t", t)]);
    close(got.re, reference, 2e-4, &s);
    // Third order (linearised KdV): u_t + u_xxx = 0, û = û₀ e^{i k³ t}.
    let s = any(
        "pdsolve(diff(u(x, t), t) + diff(diff(diff(u(x, t), x), x), x) = 0, u(x, t), list(u(x, 0) = exp(-x^2)))",
    );
    let reference =
        simpson(&|k| g_hat(k) * (k * x + k * k * k * t).cos(), -12.0, 12.0) / (2.0 * PI);
    let text = expand_integrals(&bounds(&rhs(&s)), 110);
    let got = complex_value(&text, &[("x", x), ("t", t)]);
    close(got.re, reference, 2e-4, &s);
}

#[test]
fn greens_functions() {
    let at = |s_x: f64, x: f64| [("x", x), ("s_x", s_x)];
    let (x, s0) = (0.9_f64, 2.1_f64);
    // L = d²/dx² on [0, π], Dirichlet: G = -x_<(π - x_>)/π.
    let s = any("pde_greens_function(diff(diff(u(x), x), x) = 0, u(x), list(u(0) = 0, u(pi) = 0))");
    let text = expand_sums(&rhs(&s), 800);
    close(value(&text, &at(s0, x)), -x * (PI - s0) / PI, 2e-3, &s);
    // L = d²/dx² - 1 with Neumann ends: G = -cosh(x_<) cosh(π - x_>)/sinh(π).
    let s = any(
        "pde_greens_function(diff(diff(u(x), x), x) - u(x) = 0, u(x), list(at(diff(u(x), x), x, 0) = 0, at(diff(u(x), x), x, pi) = 0))",
    );
    let text = expand_sums(&rhs(&s), 800);
    close(
        value(&text, &at(s0, x)),
        -x.cosh() * (PI - s0).cosh() / PI.sinh(),
        2e-3,
        &s,
    );
    // The heat kernel of the Dirichlet interval against the method of images.
    let s = any(
        "pde_greens_function(diff(u(x, t), t) = diff(diff(u(x, t), x), x), u(x, t), list(u(0, t) = 0, u(pi, t) = 0))",
    );
    let t = 0.3_f64;
    let k = |y: f64| (-y * y / (4.0 * t)).exp() / (4.0 * PI * t).sqrt();
    let images: f64 = (-6..=6)
        .map(|m| k(x - s0 + 2.0 * PI * f64::from(m)) - k(x + s0 + 2.0 * PI * f64::from(m)))
        .sum();
    let text = expand_sums(&rhs(&s), 40);
    close(
        value(&text, &[("x", x), ("s_x", s0), ("t", t)]),
        images,
        1e-9,
        &s,
    );
    // Free space in the plane and the Dirichlet half-plane by images.
    let lap = "diff(diff(u(x, y), x), x) + diff(diff(u(x, y), y), y) = 0";
    let s = any(&format!("pde_greens_function({lap}, u(x, y))"));
    let r = (0.3_f64 - 1.0).hypot(0.5_f64 - 0.2);
    close(
        value(
            &rhs(&s),
            &[("x", 0.3), ("y", 0.5), ("s_x", 1.0), ("s_y", 0.2)],
        ),
        r.ln() / (2.0 * PI),
        1e-12,
        &s,
    );
    let s = any(&format!(
        "pde_greens_function({lap}, u(x, y), list(u(x, 0) = 0))"
    ));
    let image = (0.3_f64 - 1.0).hypot(0.5_f64 + 0.2);
    close(
        value(
            &rhs(&s),
            &[("x", 0.3), ("y", 0.5), ("s_x", 1.0), ("s_y", 0.2)],
        ),
        (r.ln() - image.ln()) / (2.0 * PI),
        1e-12,
        &s,
    );
}

#[test]
fn laplace_in_a_ball_with_general_boundary_data() {
    let lap = "diff(diff(u(r, th, ph), r), r) + 2*diff(u(r, th, ph), r)/r + diff(diff(u(r, th, ph), th), th)/r^2 + cos(th)/sin(th)*diff(u(r, th, ph), th)/r^2 + diff(diff(u(r, th, ph), ph), ph)/(r^2*sin(th)^2)";
    let distance = |p: [f64; 3], q: [f64; 3]| {
        ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt()
    };
    let point = |r: f64, th: f64, ph: f64| {
        [
            r * th.sin() * ph.cos(),
            r * th.sin() * ph.sin(),
            r * th.cos(),
        ]
    };
    // Inside: the potential of a charge at (1, 0.8, 1.2) outside the unit ball.
    let outside = [1.0, 0.8, 1.2];
    let data = "1/((sin(th)*cos(ph) - 1)^2 + (sin(th)*sin(ph) - 4/5)^2 + (cos(th) - 6/5)^2)^(1/2)";
    let s = any(&format!(
        "pdsolve({lap} = 0, u(r, th, ph), list(u(1, th, ph) = {data}))"
    ));
    let (r, th, ph) = (0.55_f64, 1.0_f64, 0.7_f64);
    close(
        value(&rhs(&s), &[("r", r), ("th", th), ("ph", ph)]),
        1.0 / distance(point(r, th, ph), outside),
        1e-6,
        &s,
    );
    // Outside: a charge at (0.3, 0.2, 0.1) inside, decaying at infinity.
    let inside = [0.3, 0.2, 0.1];
    let data =
        "1/((sin(th)*cos(ph) - 3/10)^2 + (sin(th)*sin(ph) - 1/5)^2 + (cos(th) - 1/10)^2)^(1/2)";
    let s = any(&format!(
        "pdsolve({lap} = 0, u(r, th, ph), list(u(1, th, ph) = {data}, u(oo, th, ph) = 0))"
    ));
    let (r, th, ph) = (1.8_f64, 2.0_f64, 4.0_f64);
    close(
        value(&rhs(&s), &[("r", r), ("th", th), ("ph", ph)]),
        1.0 / distance(point(r, th, ph), inside),
        1e-6,
        &s,
    );
}

#[test]
fn burgers_hopf_cole_and_riemann_problems() {
    // Viscous Burgers u_t + u u_x = u_xx / 2: the stationary shock u = -tanh x.
    let s = any(
        "pdsolve(diff(u(x, t), t) + u(x, t)*diff(u(x, t), x) = 1/2*diff(diff(u(x, t), x), x), u(x, t), list(u(x, 0) = -tanh(x)))",
    );
    assert!(s.contains("defint"), "{s}");
    let text = expand_integrals(
        &rhs(&s).replace("-oo", "-14").replace(", oo)", ", 14)"),
        700,
    );
    let (x, t) = (0.7_f64, 0.8_f64);
    close(value(&text, &[("x", x), ("t", t)]), -x.tanh(), 2e-4, &s);
    // A moving frame and a scaling: v = 2u + 1 obeys the same law, so u_t + (2u + 1) u_x = u_xx/2 is the
    // stationary shock in the frame x - t: u = -(1 + tanh(x - t))/2 ... checked through the data u = -(tanh x + 1)/2.
    let s = any(
        "pdsolve(diff(u(x, t), t) + (2*u(x, t) + 1)*diff(u(x, t), x) = 1/2*diff(diff(u(x, t), x), x), u(x, t), list(u(x, 0) = -(tanh(x) + 1)/2))",
    );
    let text = expand_integrals(
        &rhs(&s).replace("-oo", "-14").replace(", oo)", ", 14)"),
        700,
    );
    let u_exact = |x: f64| -(x.tanh() + 1.0) / 2.0;
    // U = 2u + 1 = -tanh(x) is stationary, so u is too.
    close(value(&text, &[("x", x), ("t", t)]), u_exact(x), 2e-4, &s);
    // Riemann problems of Burgers' equation: shock s = (u_l + u_r)/2, fan u = x/t.
    let burgers = "diff(u(x, t), t) + u(x, t)*diff(u(x, t), x) = 0";
    let s = any(&format!(
        "pdsolve({burgers}, u(x, t), list(u(x, 0) = 1 - heaviside(x)))"
    ));
    close(value(&rhs(&s), &[("x", 0.4), ("t", 2.0)]), 1.0, 1e-12, &s);
    close(value(&rhs(&s), &[("x", 1.4), ("t", 2.0)]), 0.0, 1e-12, &s);
    let s = any(&format!(
        "pdsolve({burgers}, u(x, t), list(u(x, 0) = heaviside(x)))"
    ));
    close(value(&rhs(&s), &[("x", 0.5), ("t", 2.0)]), 0.25, 1e-12, &s);
    close(value(&rhs(&s), &[("x", -0.5), ("t", 2.0)]), 0.0, 1e-12, &s);
    close(value(&rhs(&s), &[("x", 2.5), ("t", 2.0)]), 1.0, 1e-12, &s);
    // Shock speed with a shifted jump and other states: u_l = 3, u_r = 1 at x₀ = 2 moves with speed 2.
    let s = any(&format!(
        "pdsolve({burgers}, u(x, t), list(u(x, 0) = 1 + 2*(1 - heaviside(x - 2))))"
    ));
    close(value(&rhs(&s), &[("x", 5.5), ("t", 2.0)]), 3.0, 1e-12, &s);
    close(value(&rhs(&s), &[("x", 6.5), ("t", 2.0)]), 1.0, 1e-12, &s);
    // Concave flux f = u (1 - u) (traffic): a jump up is a stationary shock, a jump down a fan u = (1 - ξ)/2.
    let traffic = "diff(u(x, t), t) + (1 - 2*u(x, t))*diff(u(x, t), x) = 0";
    let s = any(&format!(
        "pdsolve({traffic}, u(x, t), list(u(x, 0) = heaviside(x)))"
    ));
    close(value(&rhs(&s), &[("x", -0.3), ("t", 1.0)]), 0.0, 1e-12, &s);
    close(value(&rhs(&s), &[("x", 0.3), ("t", 1.0)]), 1.0, 1e-12, &s);
    let s = any(&format!(
        "pdsolve({traffic}, u(x, t), list(u(x, 0) = 1 - heaviside(x)))"
    ));
    close(value(&rhs(&s), &[("x", 0.5), ("t", 1.0)]), 0.25, 1e-12, &s);
    close(value(&rhs(&s), &[("x", -1.5), ("t", 1.0)]), 1.0, 1e-12, &s);
    close(value(&rhs(&s), &[("x", 1.5), ("t", 1.0)]), 0.0, 1e-12, &s);
}

#[test]
fn shock_and_rarefaction_interactions() {
    // Three constant states: the two outgoing waves of a convex flux collide
    // unless both are rarefactions, and the merged wave starts at the
    // collision point (Rankine-Hugoniot speeds, Oleinik selection).
    let burgers = "diff(u(x, t), t) + u(x, t)*diff(u(x, t), x) = 0";
    let solve = |data: &str| {
        any(&format!(
            "pdsolve({burgers}, u(x, t), list(u(x, 0) = {data}))"
        ))
    };
    let at = |s: &str, x: f64, t: f64| value(&rhs(s), &[("x", x), ("t", t)]);

    // (a) shock-shock: u = 3 | 2 | 0 at x = -1, 1. Speeds 5/2 and 1, collision
    // at T = 4/3, X = 7/3; afterwards one shock 3 -> 0 of speed 3/2.
    let s = solve("3 - heaviside(x + 1) - 2*heaviside(x - 1)");
    for (x, want) in [(0.0, 3.0), (1.7, 2.0), (2.2, 0.0)] {
        close(at(&s, x, 1.0), want, 1e-10, &s);
    }
    for (x, want) in [(3.0, 3.0), (3.6, 0.0)] {
        close(at(&s, x, 2.0), want, 1e-10, &s);
    }
    // the collision is at x = 7/3 when t = 4/3: just before it the middle
    // state still occupies a shrinking interval around 7/3
    close(at(&s, 7.0 / 3.0 - 0.1, 4.0 / 3.0 - 0.05), 2.0, 1e-10, &s);

    // (b) shock-rarefaction: u = 3 | 0 | 1. The shock (speed 3/2) overtakes
    // the fan (left edge 0) at T = 4/3, X = 1; afterwards one shock 3 -> 1
    // of speed 2.
    let s = solve("3 - 3*heaviside(x + 1) + heaviside(x - 1)");
    for (x, want) in [(0.0, 3.0), (0.7, 0.0), (1.5, 0.5), (2.5, 1.0)] {
        close(at(&s, x, 1.0), want, 1e-10, &s);
    }
    for (x, want) in [(2.0, 3.0), (2.6, 1.0)] {
        close(at(&s, x, 2.0), want, 1e-10, &s);
    }

    // (c) rarefaction-shock: u = 0 | 3 | 1. The fast fan edge (speed 3)
    // catches the shock (speed 2) at T = 2, X = 5; afterwards one rarefaction
    // 0 -> 1 centred at (5, 2).
    let s = solve("3*heaviside(x + 1) - 2*heaviside(x - 1)");
    for (x, want) in [(0.0, 1.0), (1.5, 2.5), (2.5, 3.0), (3.5, 1.0)] {
        close(at(&s, x, 1.0), want, 1e-10, &s);
    }
    for (x, want) in [(4.5, 0.0), (5.5, 0.5), (6.5, 1.0)] {
        close(at(&s, x, 3.0), want, 1e-10, &s);
    }

    // (d) rarefaction-rarefaction: u = 0 | 1 | 2; the fans touch and leave a
    // constant state between them, valid for every t.
    let s = solve("heaviside(x + 1) + heaviside(x - 1)");
    for (x, want) in [(-1.5, 0.0), (-0.5, 0.5), (1.0, 1.0), (2.5, 1.5), (3.5, 2.0)] {
        close(at(&s, x, 1.0), want, 1e-10, &s);
    }
    // the same data at a later time: the fans widen, the middle stays 1
    for (x, want) in [(-2.0, 0.0), (2.0, 1.0), (4.0, 1.5), (5.5, 2.0)] {
        close(at(&s, x, 2.0), want, 1e-10, &s);
    }

    // A concave flux (traffic f = u (1 - u)): u = 0 | 1/2 | 1 has two upward
    // shocks of speeds 1/2 and -1/2, which collide at T = 2, X = 0; the merged
    // jump 0 -> 1 is stationary.
    let traffic = "diff(u(x, t), t) + (1 - 2*u(x, t))*diff(u(x, t), x) = 0";
    let s = any(&format!(
        "pdsolve({traffic}, u(x, t), list(u(x, 0) = heaviside(x + 1)/2 + heaviside(x - 1)/2))"
    ));
    for (x, want) in [(-1.0, 0.0), (0.0, 0.5), (1.0, 1.0)] {
        close(at(&s, x, 1.0), want, 1e-10, &s);
    }
    for (x, want) in [(-0.5, 0.0), (0.5, 1.0)] {
        close(at(&s, x, 3.0), want, 1e-10, &s);
    }
}

#[test]
fn kdv_and_sine_gordon_solitons() {
    let v = |text: &str, x: f64, t: f64| value(text, &[("x", x), ("t", t)]);
    // One soliton: u = -2 κ² sech²(κ (x - 4 κ² t - x₀)).
    let s = run("pde_kdv_soliton(list(1/2), list(1/3))");
    let (x, t) = (0.9_f64, 0.4_f64);
    let kappa = 0.5_f64;
    let xi = kappa * (x - 4.0 * kappa * kappa * t - 1.0 / 3.0);
    close(
        v(&s, x, t),
        -2.0 * kappa * kappa / xi.cosh().powi(2),
        1e-10,
        &s,
    );
    // Two and three solitons solve u_t - 6 u u_x + u_xxx = 0 (finite differences).
    for (kappas, shifts) in [
        ("list(1, 3/2)", "list(0, 1/2)"),
        ("list(1/2, 1, 3/2)", "list(0, -1, 1/2)"),
    ] {
        let s = run(&format!("pde_kdv_soliton({kappas}, {shifts})"));
        let u = |x: f64, t: f64| v(&s, x, t);
        let (x, t, h) = (0.35_f64, 0.05_f64, 2e-3_f64);
        let u_t = (u(x, t + h) - u(x, t - h)) / (2.0 * h);
        let u_x = (u(x + h, t) - u(x - h, t)) / (2.0 * h);
        let u_xxx = (u(x + 2.0 * h, t) - 2.0 * u(x + h, t) + 2.0 * u(x - h, t) - u(x - 2.0 * h, t))
            / (2.0 * h.powi(3));
        close(
            u_t - 6.0 * u(x, t) * u_x + u_xxx,
            0.0,
            5e-3 * (1.0 + u_xxx.abs()),
            &s,
        );
    }
    // Sine-Gordon kinks in light-cone coordinates: u_xt = sin u.
    for (a, c) in [
        ("list(13/10)", "list(1/10)"),
        ("list(13/10, 7/10)", "list(1/10, -2/5)"),
        ("list(13/10, 7/10, 21/10)", "list(1/10, -2/5, 3/10)"),
    ] {
        let s = run(&format!("pde_sine_gordon_soliton({a}, {c})"));
        let u = |x: f64, t: f64| v(&s, x, t);
        let (x, t, h) = (0.2_f64, -0.3_f64, 1e-3_f64);
        let mixed =
            (u(x + h, t + h) - u(x + h, t - h) - u(x - h, t + h) + u(x - h, t - h)) / (4.0 * h * h);
        close(mixed, u(x, t).sin(), 2e-4, &s);
    }
}

#[test]
fn fourier_fallback_with_table_sources_and_two_dimensions() {
    use super::extension_tests::complex_value;
    let (x, y, t) = (0.8_f64, -0.4_f64, 0.3_f64);
    let bounds = |text: &str| text.replace("-oo", "-9").replace(", oo)", ", 9)");
    let g_hat = |k: f64| PI.sqrt() * (-k * k / 4.0).exp();
    let airy =
        |x: f64| simpson(&|k| g_hat(k) * (k * x + k * k * k * t).cos(), -12.0, 12.0) / (2.0 * PI);
    // Two dimensions: u_t + u_xxx + u_yyy = 0 with a Gaussian is the product of two one-dimensional solutions.
    let s = any(
        "pdsolve(diff(u(x, y, t), t) + diff(diff(diff(u(x, y, t), x), x), x) + diff(diff(diff(u(x, y, t), y), y), y) = 0, u(x, y, t), list(u(x, y, 0) = exp(-x^2 - y^2)))",
    );
    assert!(s.contains("defint") && !s.contains("fourier"), "{s}");
    let text = expand_integrals(&bounds(&rhs(&s)), 70);
    let got = complex_value(&text, &[("x", x), ("y", y), ("t", t)]);
    close(got.re, airy(x) * airy(y), 5e-4, &s);
    // A source: u_t + u_xxx = e^{-x²}, zero data: û = ĝ (e^{i k³ t} - 1)/(i k³).
    let s = any(
        "pdsolve(diff(u(x, t), t) + diff(diff(diff(u(x, t), x), x), x) = exp(-x^2), u(x, t), list(u(x, 0) = 0))",
    );
    let reference = simpson(
        &|k| {
            let k3 = k * k * k;
            let kernel = if k3.abs() < 1e-9 {
                t * (k * x).cos()
            } else {
                ((k3 * t).sin() * (k * x).cos() + ((k3 * t).cos() - 1.0) * (k * x).sin()) / k3
            };
            g_hat(k) * kernel
        },
        -12.0,
        12.0,
    ) / (2.0 * PI);
    let text = expand_integrals(&bounds(&rhs(&s)), 220);
    let got = complex_value(&text, &[("x", x), ("t", t)]);
    close(got.re, reference, 5e-4, &s);
}

#[test]
fn greens_functions_periodic_damped_and_heat_kernels() {
    use super::extension_tests::complex_value;
    let (x, s0, t) = (0.9_f64, 2.1_f64, 0.3_f64);
    let kernel = |y: f64| (-y * y / (4.0 * t)).exp() / (4.0 * PI * t).sqrt();
    // The periodic heat kernel on [0, 2π] by images.
    let s = any(
        "pde_greens_function(diff(u(x, t), t) = diff(diff(u(x, t), x), x), u(x, t), list(u(0, t) = u(2*pi, t)))",
    );
    let images: f64 = (-6..=6)
        .map(|m| kernel(x - s0 + 2.0 * PI * f64::from(m)))
        .sum();
    let text = expand_sums(&rhs(&s), 40);
    close(
        value(&text, &[("x", x), ("s_x", s0), ("t", t)]),
        images,
        1e-9,
        &s,
    );
    // Damped waves: u_tt + u_t = u_xx on [0, π] with Dirichlet ends: modes e^{-t/2} sin(ω t)/ω, ω² = n² - 1/4.
    let s = any(
        "pde_greens_function(diff(diff(u(x, t), t), t) + diff(u(x, t), t) = diff(diff(u(x, t), x), x), u(x, t), list(u(0, t) = 0, u(pi, t) = 0))",
    );
    let reference: f64 = (1..=40)
        .map(|n| {
            let nf = f64::from(n);
            let omega = (nf * nf - 0.25).sqrt();
            2.0 / PI * (nf * x).sin() * (nf * s0).sin() * (-t / 2.0).exp() * (omega * t).sin()
                / omega
        })
        .sum();
    let text = expand_sums(&rhs(&s), 40);
    close(
        complex_value(&text, &[("x", x), ("s_x", s0), ("t", t)]).re,
        reference,
        1e-9,
        &s,
    );
    // The heat kernel of the whole line and of the half-line with a Dirichlet end.
    let s = any("pde_greens_function(diff(u(x, t), t) = diff(diff(u(x, t), x), x), u(x, t))");
    close(
        value(&rhs(&s), &[("x", x), ("s_x", s0), ("t", t)]),
        kernel(x - s0),
        1e-12,
        &s,
    );
    let s = any(
        "pde_greens_function(diff(u(x, t), t) = diff(diff(u(x, t), x), x), u(x, t), list(u(0, t) = 0))",
    );
    close(
        value(&rhs(&s), &[("x", x), ("s_x", s0), ("t", t)]),
        kernel(x - s0) - kernel(x + s0),
        1e-12,
        &s,
    );
}

#[test]
fn annulus_with_a_zero_mode() {
    let heat = "diff(u(r, t), t) = diff(diff(u(r, t), r), r) + diff(u(r, t), r)/r";
    // Inner u_r - u = 0 at r = 1 and u_r - u/(2 (1 + ln 2)) = 0 at r = 2: F = ln r + 1 has eigenvalue 0.
    let s = any(&format!(
        "pdsolve({heat}, u(r, t), list(at(diff(u(r, t), r), r, 1) - u(1, t) = 0, at(diff(u(r, t), r), r, 2) - u(2, t)/(2*(1 + ln(2))) = 0, u(r, 0) = 1))"
    ));
    let (r, t) = (1.4_f64, 0.02_f64);
    let (pa, pb) = (-1.0_f64, -1.0 / (2.0 * (1.0 + 2.0_f64.ln())));
    let cond_jy = |second: bool, p: f64, k: f64, radius: f64| {
        let z = k * radius;
        if second {
            p * crate::kernels::special::bessel_y0(z) - k * crate::kernels::special::bessel_y1(z)
        } else {
            p * crate::kernels::special::bessel_j0(z) - k * crate::kernels::special::bessel_j1(z)
        }
    };
    let cond_ik = |second: bool, p: f64, k: f64, radius: f64| {
        let z = k * radius;
        if second {
            p * crate::kernels::special::bessel_k0(z) - k * crate::kernels::special::bessel_k1(z)
        } else {
            p * crate::kernels::special::bessel_i0(z) + k * crate::kernels::special::bessel_i1(z)
        }
    };
    let coefficient = |mode: &dyn Fn(f64) -> f64| {
        simpson(&|q| q * mode(q), 1.0, 2.0) / simpson(&|q| q * mode(q) * mode(q), 1.0, 2.0)
    };
    let zero_mode = |q: f64| q.ln() + 1.0;
    let mut reference = coefficient(&zero_mode) * zero_mode(r);
    let g_neg = |k: f64| {
        cond_ik(true, pa, k, 1.0) * cond_ik(false, pb, k, 2.0)
            - cond_ik(false, pa, k, 1.0) * cond_ik(true, pb, k, 2.0)
    };
    for kappa in roots(&g_neg, 12.0, 4) {
        let mode = move |q: f64| {
            cond_ik(true, pa, kappa, 1.0) * crate::kernels::special::bessel_i0(kappa * q)
                - cond_ik(false, pa, kappa, 1.0) * crate::kernels::special::bessel_k0(kappa * q)
        };
        reference += coefficient(&mode) * (kappa * kappa * t).exp() * mode(r);
    }
    let g_pos = |k: f64| {
        cond_jy(true, pa, k, 1.0) * cond_jy(false, pb, k, 2.0)
            - cond_jy(false, pa, k, 1.0) * cond_jy(true, pb, k, 2.0)
    };
    for k in roots(&g_pos, 130.0, 40) {
        let mode = move |q: f64| {
            cond_jy(true, pa, k, 1.0) * crate::kernels::special::bessel_j0(k * q)
                - cond_jy(false, pa, k, 1.0) * crate::kernels::special::bessel_y0(k * q)
        };
        reference += coefficient(&mode) * (-k * k * t).exp() * mode(r);
    }
    let text = expand_integrals(&expand_sums(&rhs(&s), 40), 160);
    close(value(&text, &[("r", r), ("t", t)]), reference, 3e-4, &s);
}

#[test]
fn hankel_transform_for_radial_evolution_on_whole_space() {
    let bounds = |text: &str| {
        text.replace(", rho, 0, oo)", ", rho, 0, 7)")
            .replace(", k, 0, oo)", ", k, 0, 8)")
    };
    // Plane: u_t = Δu with u(r, 0) = e^{-r²}: u = e^{-r²/(1+4t)}/(1+4t).
    let s = any(
        "pdsolve(diff(u(r, t), t) = diff(diff(u(r, t), r), r) + diff(u(r, t), r)/r, u(r, t), list(u(r, 0) = exp(-r^2)))",
    );
    assert!(s.contains("besselj(0,"), "{s}");
    let (r, t) = (0.7_f64, 0.5_f64);
    let text = expand_integrals(&bounds(&rhs(&s)), 100);
    close(
        value(&text, &[("r", r), ("t", t)]),
        (-r * r / (1.0 + 4.0 * t)).exp() / (1.0 + 4.0 * t),
        5e-4,
        &s,
    );
    // Space: u_tt = Δu with u(r, 0) = e^{-r²}, u_t = 0: u = [(r + t) e^{-(r+t)²} + (r - t) e^{-(r-t)²}]/(2 r).
    let s = any(
        "pdsolve(diff(diff(u(r, t), t), t) = diff(diff(u(r, t), r), r) + 2*diff(u(r, t), r)/r, u(r, t), list(u(r, 0) = exp(-r^2), at(diff(u(r, t), t), t, 0) = 0))",
    );
    assert!(s.contains("sin("), "{s}");
    let (r, t) = (1.1_f64, 0.6_f64);
    let text = expand_integrals(&bounds(&rhs(&s)), 100);
    let exact =
        ((r + t) * (-(r + t) * (r + t)).exp() + (r - t) * (-(r - t) * (r - t)).exp()) / (2.0 * r);
    close(value(&text, &[("r", r), ("t", t)]), exact, 1e-3, &s);
}

#[test]
fn numeric_modal_expansions_on_disks_and_balls() {
    use crate::kernels::special::bessel_j0;
    use crate::kernels::special::bessel_j1;
    let coefficient = |mode: &dyn Fn(f64) -> f64, f: &dyn Fn(f64) -> f64, weight: i32| {
        simpson(&|q| q.powi(weight) * f(q) * mode(q), 0.0, 1.0)
            / simpson(&|q| q.powi(weight) * mode(q) * mode(q), 0.0, 1.0)
    };
    // Disk heat with non-polynomial data e^{-r} (1 + r cos θ), Dirichlet at r = 1.
    let lap_disk = "diff(diff(u(r, th, t), r), r) + diff(u(r, th, t), r)/r + diff(diff(u(r, th, t), th), th)/r^2";
    let s = any(&format!(
        "pdsolve(diff(u(r, th, t), t) = {lap_disk}, u(r, th, t), list(u(1, th, t) = 0, u(r, th, 0) = exp(-r)*(1 + r*cos(th))))"
    ));
    assert!(!s.contains("defint"), "{s}");
    let (r, th, t) = (0.45_f64, 0.8_f64, 0.02_f64);
    let mut reference = 0.0;
    for k in roots(&bessel_j0, 60.0, 20) {
        let c = coefficient(&|q| bessel_j0(k * q), &|q| (-q).exp(), 1);
        reference += c * (-k * k * t).exp() * bessel_j0(k * r);
    }
    for k in roots(&bessel_j1, 60.0, 20) {
        let c = coefficient(&|q| bessel_j1(k * q), &|q| q * (-q).exp(), 1);
        reference += c * (-k * k * t).exp() * bessel_j1(k * r) * th.cos();
    }
    close(
        value(&rhs(&s), &[("r", r), ("th", th), ("t", t)]),
        reference,
        1e-6,
        &s,
    );
    // Disk wave, radial data (1 - r²)² cos r: u = Σ c cos(κ t) J_0(κ r).
    let s = any(
        "pdsolve(diff(diff(u(r, t), t), t) = diff(diff(u(r, t), r), r) + diff(u(r, t), r)/r, u(r, t), list(u(1, t) = 0, u(r, 0) = (1 - r^2)^2*cos(r), at(diff(u(r, t), t), t, 0) = 0))",
    );
    assert!(!s.contains("defint"), "{s}");
    let (r, t) = (0.35_f64, 0.3_f64);
    let mut reference = 0.0;
    for k in roots(&bessel_j0, 140.0, 40) {
        let c = coefficient(
            &|q| bessel_j0(k * q),
            &|q| (1.0 - q * q).powi(2) * q.cos(),
            1,
        );
        reference += c * (k * t).cos() * bessel_j0(k * r);
    }
    close(value(&rhs(&s), &[("r", r), ("t", t)]), reference, 2e-4, &s);
    // Axisymmetric ball heat, data e^{-r} (2 + r cos θ), Dirichlet at r = 1.
    let lap_ball = "diff(diff(u(r, th, t), r), r) + 2*diff(u(r, th, t), r)/r + diff(diff(u(r, th, t), th), th)/r^2 + cos(th)/sin(th)*diff(u(r, th, t), th)/r^2";
    let s = any(&format!(
        "pdsolve(diff(u(r, th, t), t) = {lap_ball}, u(r, th, t), list(u(1, th, t) = 0, u(r, th, 0) = exp(-r)*(2 + r*cos(th))))"
    ));
    assert!(!s.contains("defint"), "{s}");
    let (r, th, t) = (0.5_f64, 1.1_f64, 0.01_f64);
    let mut reference = 0.0;
    for k in 1..=40 {
        let kappa = PI * f64::from(k);
        let c = coefficient(
            &|q| {
                if q < 1e-12 {
                    kappa
                } else {
                    (kappa * q).sin() / q
                }
            },
            &|q| 2.0 * (-q).exp(),
            2,
        );
        reference += c * (-kappa * kappa * t).exp() * (kappa * r).sin() / r;
    }
    let j1 = |x: f64| {
        if x < 1e-12 {
            0.0
        } else {
            (PI / (2.0 * x)).sqrt() * bessel_j(1.5, x)
        }
    };
    for kappa in roots(&|x| bessel_j(1.5, x), 120.0, 30) {
        let c = coefficient(&|q| j1(kappa * q), &|q| q * (-q).exp(), 2);
        reference += c * (-kappa * kappa * t).exp() * j1(kappa * r) * th.cos();
    }
    close(
        value(&rhs(&s), &[("r", r), ("th", th), ("t", t)]),
        reference,
        1e-5,
        &s,
    );
}

/// `v_t = v_rr + (d/r) v_r - (λ/r²) v` on `[1, 13]` (Dirichlet at both ends) by Crank-Nicolson; the
/// value at `at`.
fn exterior_heat_reference(
    d: f64,
    lambda: f64,
    f: &dyn Fn(f64) -> f64,
    t_end: f64,
    at: f64,
) -> f64 {
    let (n, steps) = (1200_usize, 250_usize);
    let h = 12.0 / n as f64;
    let dt = t_end / steps as f64;
    let radius = |i: usize| 1.0 + h * i as f64;
    let mut u: Vec<f64> = (0..=n).map(|i| f(radius(i))).collect();
    u[0] = 0.0;
    for _ in 0..steps {
        let (mut lower, mut diag, mut upper, mut rhs_) = (
            vec![0.0; n + 1],
            vec![1.0; n + 1],
            vec![0.0; n + 1],
            vec![0.0; n + 1],
        );
        for i in 1..n {
            let r_i = radius(i);
            let (a, b) = (
                1.0 / (h * h) - d / (2.0 * h * r_i),
                1.0 / (h * h) + d / (2.0 * h * r_i),
            );
            let c = -2.0 / (h * h) - lambda / (r_i * r_i);
            let operator = a * u[i - 1] + c * u[i] + b * u[i + 1];
            rhs_[i] = u[i] + 0.5 * dt * operator;
            lower[i] = -0.5 * dt * a;
            diag[i] = 1.0 - 0.5 * dt * c;
            upper[i] = -0.5 * dt * b;
        }
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
    u[((at - 1.0) / h).round() as usize]
}

/// `v_tt = v_rr + (d/r) v_r - (λ/r²) v` on `[1, 21]` with `v = 0` at both ends, `v_t(0) = 0`, by leapfrog.
fn exterior_wave_reference(
    d: f64,
    lambda: f64,
    f: &dyn Fn(f64) -> f64,
    t_end: f64,
    at: f64,
) -> f64 {
    let (h, dt) = (0.005_f64, 0.0025_f64);
    let n = 4000_usize;
    let steps = (t_end / dt).round() as usize;
    let radius = |i: usize| 1.0 + h * i as f64;
    let acceleration = |u: &[f64]| -> Vec<f64> {
        let mut out = vec![0.0; u.len()];
        for i in 1..u.len() - 1 {
            let r_i = radius(i);
            out[i] = (u[i + 1] - 2.0 * u[i] + u[i - 1]) / (h * h)
                + d / r_i * (u[i + 1] - u[i - 1]) / (2.0 * h)
                - lambda / (r_i * r_i) * u[i];
        }
        out
    };
    let mut previous: Vec<f64> = (0..=n).map(|i| f(radius(i))).collect();
    previous[0] = 0.0;
    let a0 = acceleration(&previous);
    let mut current: Vec<f64> = previous
        .iter()
        .zip(&a0)
        .map(|(u, a)| u + 0.5 * dt * dt * a)
        .collect();
    for _ in 1..steps {
        let a = acceleration(&current);
        let next: Vec<f64> = (0..=n)
            .map(|i| 2.0 * current[i] - previous[i] + dt * dt * a[i])
            .collect();
        previous = current;
        current = next;
    }
    current[((at - 1.0) / h).round() as usize]
}

#[test]
fn weber_transforms_with_angular_dependence_and_exterior_waves() {
    let bounds = |text: &str| {
        text.replace(", rho, 1, oo)", ", rho, 1, 9)")
            .replace(", k, 0, oo)", ", k, 0, 8)")
    };
    // Heat outside the unit disk, first angular harmonic: u = v(r, t) cos θ.
    let lap = "diff(diff(u(r, th, t), r), r) + diff(u(r, th, t), r)/r + diff(diff(u(r, th, t), th), th)/r^2";
    let s = any(&format!(
        "pdsolve(diff(u(r, th, t), t) = {lap}, u(r, th, t), list(u(1, th, t) = 0, u(oo, th, t) = 0, u(r, th, 0) = exp(-(r - 4)^2)*cos(th)))"
    ));
    let (r, th, t) = (2.0_f64, 0.6_f64, 0.4_f64);
    let text = expand_integrals(&bounds(&rhs(&s)), 90);
    let reference =
        exterior_heat_reference(1.0, 1.0, &|q| (-(q - 4.0) * (q - 4.0)).exp(), t, r) * th.cos();
    close(
        value(&text, &[("r", r), ("th", th), ("t", t)]),
        reference,
        4e-3,
        &s,
    );
    // Waves outside the unit disk, radial: u_tt = Δu.
    let s = any(
        "pdsolve(diff(diff(u(r, t), t), t) = diff(diff(u(r, t), r), r) + diff(u(r, t), r)/r, u(r, t), list(u(1, t) = 0, u(oo, t) = 0, u(r, 0) = exp(-2*(r - 4)^2), at(diff(u(r, t), t), t, 0) = 0))",
    );
    let (r, t) = (3.0_f64, 1.2_f64);
    let text = expand_integrals(&bounds(&rhs(&s)), 90);
    let reference =
        exterior_wave_reference(1.0, 0.0, &|q| (-2.0 * (q - 4.0) * (q - 4.0)).exp(), t, r);
    close(value(&text, &[("r", r), ("t", t)]), reference, 6e-3, &s);
    // Heat outside the unit ball, degree one: u = r^{-1/2} F cos θ with the Bessel order 3/2.
    let lap = "diff(diff(u(r, th, t), r), r) + 2*diff(u(r, th, t), r)/r + diff(diff(u(r, th, t), th), th)/r^2 + cos(th)/sin(th)*diff(u(r, th, t), th)/r^2";
    let s = any(&format!(
        "pdsolve(diff(u(r, th, t), t) = {lap}, u(r, th, t), list(u(1, th, t) = 0, u(oo, th, t) = 0, u(r, th, 0) = exp(-(r - 4)^2)*cos(th)))"
    ));
    assert!(s.contains("legendre") || s.contains("cos(th)"), "{s}");
    let (r, th, t) = (2.0_f64, 0.6_f64, 0.4_f64);
    let text = expand_integrals(&bounds(&rhs(&s)), 90);
    let reference =
        exterior_heat_reference(2.0, 2.0, &|q| (-(q - 4.0) * (q - 4.0)).exp(), t, r) * th.cos();
    close(
        value(&text, &[("r", r), ("th", th), ("t", t)]),
        reference,
        4e-3,
        &s,
    );
}

#[test]
fn wave_kernels_mellin_and_nls_soliton() {
    use super::extension_tests::complex_value;
    // Impulse responses of the wave equation on whole space.
    let w1 = "diff(diff(u(x, t), t), t) = 4*diff(diff(u(x, t), x), x)";
    let s = any(&format!("pde_greens_function({w1}, u(x, t))"));
    let (x, s0, t) = (0.7_f64, 0.2_f64, 0.5_f64);
    close(
        value(&rhs(&s), &[("x", x), ("s_x", s0), ("t", t)]),
        0.25,
        1e-12,
        &s,
    );
    close(
        value(&rhs(&s), &[("x", 3.0), ("s_x", s0), ("t", t)]),
        0.0,
        1e-12,
        &s,
    );
    let s = any(&format!(
        "pde_greens_function({w1}, u(x, t), list(u(0, t) = 0))"
    ));
    close(
        value(&rhs(&s), &[("x", x), ("s_x", 0.3), ("t", 0.6)]),
        0.0,
        1e-12,
        &s,
    );
    close(
        value(&rhs(&s), &[("x", 0.25), ("s_x", 0.3), ("t", 0.1)]),
        0.25,
        1e-12,
        &s,
    );
    // The plane: ∫ G ds = t (the response to u_t(0) = 1), radial quadrature with r = c t sin φ.
    let w2 = "diff(diff(u(x, y, t), t), t) = diff(diff(u(x, y, t), x), x) + diff(diff(u(x, y, t), y), y)";
    let s = any(&format!("pde_greens_function({w2}, u(x, y, t))"));
    let t = 0.8_f64;
    let integral = simpson(
        &|phi| {
            let r = t * phi.sin();
            2.0 * PI
                * r
                * t
                * phi.cos()
                * value(
                    &rhs(&s),
                    &[("x", r), ("y", 0.0), ("s_x", 0.0), ("s_y", 0.0), ("t", t)],
                )
        },
        0.0,
        PI / 2.0 - 1e-6,
    );
    close(integral, t, 1e-3, &s);
    // Mellin: u_t = x² u_xx + x u_x on x > 0 is the heat equation in s = ln x.
    let s = any(
        "pdsolve(diff(u(x, t), t) = x^2*diff(diff(u(x, t), x), x) + x*diff(u(x, t), x), u(x, t), list(u(x, 0) = exp(-ln(x)^2)))",
    );
    let (x, t) = (1.7_f64, 0.3_f64);
    let text = expand_integrals(
        &rhs(&s).replace("-oo", "-12").replace(", oo)", ", 12)"),
        600,
    );
    let sx = x.ln();
    close(
        value(&text, &[("x", x), ("t", t)]),
        (-sx * sx / (1.0 + 4.0 * t)).exp() / (1.0 + 4.0 * t).sqrt(),
        1e-5,
        &s,
    );
    // The bright NLS soliton solves i ψ_t + ψ_xx/2 + |ψ|² ψ = 0.
    let s = run("pde_nls_soliton(list(3/2, 1/2, 1/5, 1/3))");
    let psi = |x: f64, t: f64| complex_value(&s, &[("x", x), ("t", t)]);
    let (x, t, h) = (0.4_f64, 0.3_f64, 1e-3_f64);
    let u_t = (psi(x, t + h) - psi(x, t - h)) / (2.0 * h);
    let u_xx = (psi(x + h, t) - 2.0 * psi(x, t) + psi(x - h, t)) / (h * h);
    let residual =
        num_complex::Complex64::new(0.0, 1.0) * u_t + 0.5 * u_xx + psi(x, t).norm_sqr() * psi(x, t);
    close(residual.norm(), 0.0, 5e-4, &s);
}

#[test]
fn riemann_problem_with_a_non_convex_flux() {
    // u_t + (u³)_x = 0: from -1 to 1 a shock of speed 3/4 to the tangent state 1/2 followed by a fan
    // u = sqrt(ξ/3); from 1 to -1 the mirror image.
    let eq = "diff(u(x, t), t) + 3*u(x, t)^2*diff(u(x, t), x) = 0";
    let s = any(&format!(
        "pdsolve({eq}, u(x, t), list(u(x, 0) = -1 + 2*heaviside(x)))"
    ));
    close(value(&rhs(&s), &[("x", 0.5), ("t", 1.0)]), -1.0, 1e-9, &s);
    close(
        value(&rhs(&s), &[("x", 1.5), ("t", 1.0)]),
        0.5_f64.sqrt(),
        1e-9,
        &s,
    );
    close(
        value(&rhs(&s), &[("x", 2.5), ("t", 1.0)]),
        (2.5_f64 / 3.0).sqrt(),
        1e-9,
        &s,
    );
    close(value(&rhs(&s), &[("x", 4.0), ("t", 1.0)]), 1.0, 1e-9, &s);
    let s = any(&format!(
        "pdsolve({eq}, u(x, t), list(u(x, 0) = 1 - 2*heaviside(x)))"
    ));
    close(value(&rhs(&s), &[("x", 0.5), ("t", 1.0)]), 1.0, 1e-9, &s);
    close(
        value(&rhs(&s), &[("x", 1.5), ("t", 1.0)]),
        -0.5_f64.sqrt(),
        1e-9,
        &s,
    );
    close(value(&rhs(&s), &[("x", 4.0), ("t", 1.0)]), -1.0, 1e-9, &s);
}

#[test]
fn mathieu_problems_by_hills_method() {
    use super::hill::hill_eigen;
    use super::hill::hill_mode;
    // Characteristic values of Mathieu's equation for q = 1 on a period 2π (Abramowitz–Stegun 20.1):
    // a₀, b₁, b₂ ... in ascending order: -0.4551386, -0.1102488, 1.8591081, 3.9170248, 4.3713010.
    let m = 256_usize;
    let mut args = vec![2.0 * PI, 24.0];
    args.extend((0..m).map(|j| 2.0 * (2.0 * 2.0 * PI * j as f64 / m as f64).cos()));
    let want = [
        -0.455_138_60,
        -0.110_248_82,
        1.859_108_07,
        3.917_024_77,
        4.371_300_98,
    ];
    for (k, w) in want.iter().enumerate() {
        let mut call = vec![(k + 1) as f64];
        call.extend_from_slice(&args);
        close(hill_eigen(&call), *w, 1e-6, "Mathieu characteristic value");
    }
    // The eigenfunctions are orthonormal on the period.
    let mode = |index: f64, x: f64| {
        let mut call = vec![index, x];
        call.extend_from_slice(&args);
        hill_mode(&call)
    };
    close(
        simpson(&|x| mode(1.0, x) * mode(1.0, x), 0.0, 2.0 * PI),
        1.0,
        1e-8,
        "norm",
    );
    close(
        simpson(&|x| mode(1.0, x) * mode(2.0, x), 0.0, 2.0 * PI),
        0.0,
        1e-8,
        "orthogonality",
    );
    // The heat equation u_t = u_xx - 2 cos(2x) u on a periodic interval: each term solves it (residual).
    let eq = "diff(u(x, t), t) = diff(diff(u(x, t), x), x) - 2*cos(2*x)*u(x, t)";
    let s = any(&format!(
        "pdsolve({eq}, u(x, t), list(u(0, t) = u(2*pi, t), u(x, 0) = 1 + sin(x)))"
    ));
    assert!(s.contains("hill_mode") && s.contains("hill_eigen"), "{s}");
    let text = expand_integrals(&expand_sums(&rhs(&s), 10), 130);
    let u = |x: f64, t: f64| value(&text, &[("x", x), ("t", t)]);
    let (x0, t0, h) = (1.0_f64, 0.5_f64, 2e-3_f64);
    let residual = (u(x0, t0 + h) - u(x0, t0 - h)) / (2.0 * h)
        - (u(x0 + h, t0) - 2.0 * u(x0, t0) + u(x0 - h, t0)) / (h * h)
        + 2.0 * (2.0 * x0).cos() * u(x0, t0);
    close(residual, 0.0, 2e-3, "residual");
    close(u(0.3, t0), u(0.3 + 2.0 * PI, t0), 1e-6, "periodic");
}

#[test]
fn tanh_method_finds_travelling_waves() {
    type Law = dyn Fn(f64, f64, f64, f64, f64) -> f64;
    // (equation, residual from u, u_t, u_x, u_xx, u_xxx)
    let cases: Vec<(&str, Box<Law>)> = vec![
        (
            "diff(u(x, t), t) + 6*u(x, t)*diff(u(x, t), x) + diff(diff(diff(u(x, t), x), x), x) = 0",
            Box::new(|u, ut, ux, _, uxxx| ut + 6.0 * u * ux + uxxx),
        ),
        (
            "diff(u(x, t), t) = diff(diff(u(x, t), x), x) + u(x, t)*(1 - u(x, t))",
            Box::new(|u, ut, _, uxx, _| ut - uxx - u * (1.0 - u)),
        ),
        (
            "diff(u(x, t), t) = diff(diff(u(x, t), x), x) + u(x, t) - u(x, t)^3",
            Box::new(|u, ut, _, uxx, _| ut - uxx - u + u * u * u),
        ),
        (
            "diff(u(x, t), t) + u(x, t)*diff(u(x, t), x) = 1/2*diff(diff(u(x, t), x), x)",
            Box::new(|u, ut, ux, uxx, _| ut + u * ux - 0.5 * uxx),
        ),
    ];
    for (equation, law) in cases {
        let s = any(&format!("pdsolve({equation}, u(x, t))"));
        assert!(s.contains("tanh"), "{s}");
        let text = rhs(&s);
        let u = |x: f64, t: f64| value(&text, &[("x", x), ("t", t)]);
        let (x, t, h) = (0.35_f64, 0.2_f64, 2e-3_f64);
        let ut = (u(x, t + h) - u(x, t - h)) / (2.0 * h);
        let ux = (u(x + h, t) - u(x - h, t)) / (2.0 * h);
        let uxx = (u(x + h, t) - 2.0 * u(x, t) + u(x - h, t)) / (h * h);
        let uxxx = (u(x + 2.0 * h, t) - 2.0 * u(x + h, t) + 2.0 * u(x - h, t) - u(x - 2.0 * h, t))
            / (2.0 * h.powi(3));
        close(
            law(u(x, t), ut, ux, uxx, uxxx),
            0.0,
            2e-3 * (1.0 + uxxx.abs()),
            &s,
        );
        // Non-trivial: the wave varies in space.
        assert!(
            (u(-3.0, t) - u(3.0, t)).abs() > 1e-3 || (u(0.0, t) - u(3.0, t)).abs() > 1e-3,
            "{s}"
        );
    }
}

#[test]
fn separable_variable_coefficients_in_two_variables() {
    // u_t = (1 + x)² u_xx + (1 + x) u_x + u_yy: an Euler mode times sin(y).
    let eq = "diff(u(x, y, t), t) = (1 + x)^2*diff(diff(u(x, y, t), x), x) + (1 + x)*diff(u(x, y, t), x) + diff(diff(u(x, y, t), y), y)";
    let l = 2.0_f64.ln();
    let data = "sin(pi*ln(1 + x)/ln(2))*sin(y)";
    let s = any(&format!(
        "pdsolve({eq}, u(x, y, t), list(u(0, y, t) = 0, u(1, y, t) = 0, u(x, 0, t) = 0, u(x, pi, t) = 0, u(x, y, 0) = {data}))"
    ));
    assert!(s.contains("sl_mode"), "{s}");
    let (x, y, t) = (0.45_f64, 1.1_f64, 0.01_f64);
    let exact = (-((PI / l).powi(2) + 1.0) * t).exp() * (PI * x.ln_1p() / l).sin() * y.sin();
    let text = expand_integrals(&expand_sums(&rhs(&s), 4), 120);
    close(
        value(&text, &[("x", x), ("y", y), ("t", t)]),
        exact,
        5e-4,
        &s,
    );
}

#[test]
fn reflectionless_darboux_and_sech_ansatz() {
    let v = |text: &str, x: f64, t: f64| value(text, &[("x", x), ("t", t)]);
    // Reflectionless data -N (N + 1) κ² sech²(κ x) at t = 0.
    for (n, kappa) in [(1_i32, 1.0_f64), (2, 1.0), (2, 0.7), (3, 1.0)] {
        let s = run(&format!("pde_kdv_reflectionless({n}, {kappa})"));
        for x in [-1.3_f64, 0.0, 0.4, 1.9] {
            let want = -f64::from(n * (n + 1)) * kappa * kappa / (kappa * x).cosh().powi(2);
            close(v(&s, x, 0.0), want, 1e-9, &s);
        }
    }
    // The N = 2 solution solves KdV.
    let s = run("pde_kdv_reflectionless(2, 1)");
    let u = |x: f64, t: f64| v(&s, x, t);
    let (x, t, h) = (0.35_f64, 0.05_f64, 2e-3_f64);
    let u_t = (u(x, t + h) - u(x, t - h)) / (2.0 * h);
    let u_x = (u(x + h, t) - u(x - h, t)) / (2.0 * h);
    let u_xxx = (u(x + 2.0 * h, t) - 2.0 * u(x + h, t) + 2.0 * u(x - h, t) - u(x - 2.0 * h, t))
        / (2.0 * h.powi(3));
    close(
        u_t - 6.0 * u(x, t) * u_x + u_xxx,
        0.0,
        5e-3 * (1.0 + u_xxx.abs()),
        &s,
    );
    // Darboux: the vacuum with psi = cosh(k x - 4 k^3 t) gives the one-soliton.
    let s = run("pde_darboux_kdv(0, cosh(x - 4*t))");
    let (x, t) = (0.8_f64, 0.3_f64);
    close(v(&s, x, t), -2.0 / (x - 4.0 * t).cosh().powi(2), 1e-9, &s);
    // Sech ansatz: the focusing mKdV soliton k sech(k (x - k² t)).
    let s = any(
        "pdsolve(diff(u(x, t), t) + 6*u(x, t)^2*diff(u(x, t), x) + diff(diff(diff(u(x, t), x), x), x) = 0, u(x, t))",
    );
    assert!(s.contains("cosh"), "{s}");
    let text = rhs(&s);
    let u = |x: f64, t: f64| value(&text, &[("x", x), ("t", t)]);
    let (x, t, h) = (0.35_f64, 0.2_f64, 2e-3_f64);
    let ut = (u(x, t + h) - u(x, t - h)) / (2.0 * h);
    let ux = (u(x + h, t) - u(x - h, t)) / (2.0 * h);
    let uxxx = (u(x + 2.0 * h, t) - 2.0 * u(x + h, t) + 2.0 * u(x - h, t) - u(x - 2.0 * h, t))
        / (2.0 * h.powi(3));
    close(
        ut + 6.0 * u(x, t) * u(x, t) * ux + uxxx,
        0.0,
        2e-3 * (1.0 + uxxx.abs()),
        &s,
    );
    assert!(u(0.0, t).abs() > 1e-3, "{s}");
}

#[test]
fn laplace_in_time_signalling_on_a_half_line() {
    let quiet = "u(x, 0) = 0, at(diff(u(x, t), t), t, 0) = 0";
    // Pure wave, Dirichlet signal: u = g(t - x/c).
    let s = any(&format!(
        "pdsolve(diff(diff(u(x, t), t), t) = 4*diff(diff(u(x, t), x), x), u(x, t), list(u(0, t) = sin(t), {quiet}))"
    ));
    let text = rhs(&s);
    close(
        value(&text, &[("x", 1.0), ("t", 2.0)]),
        (2.0_f64 - 0.5).sin(),
        1e-9,
        &s,
    );
    close(
        value(&text, &[("x", 3.0), ("t", 1.0)]),
        0.0,
        1e-9,
        "ahead of the front",
    );
    // Neumann signal: u = -c ∫_0^{t-a} g.
    let s = any(&format!(
        "pdsolve(diff(diff(u(x, t), t), t) = diff(diff(u(x, t), x), x), u(x, t), list(at(diff(u(x, t), x), x, 0) = cos(t), {quiet}))"
    ));
    close(
        value(&rhs(&s), &[("x", 1.0), ("t", 3.0)]),
        -(2.0_f64).sin(),
        1e-9,
        &s,
    );
    // Klein–Gordon, telegraph with damping, Klein–Gordon with Neumann data: residuals of the retarded integrals.
    let cases = [
        (
            "diff(diff(u(x, t), t), t) = diff(diff(u(x, t), x), x) - 4*u(x, t)",
            "u(0, t) = sin(t)",
            1.0_f64,
            0.0_f64,
            4.0_f64,
        ),
        (
            "diff(diff(u(x, t), t), t) + 2*diff(u(x, t), t) = diff(diff(u(x, t), x), x)",
            "u(0, t) = cos(t)",
            1.0,
            2.0,
            0.0,
        ),
        (
            "diff(diff(u(x, t), t), t) = diff(diff(u(x, t), x), x) - 3*u(x, t)",
            "at(diff(u(x, t), x), x, 0) = sin(t)",
            1.0,
            0.0,
            3.0,
        ),
        (
            "diff(diff(u(x, t), t), t) = diff(diff(u(x, t), x), x) + 2*u(x, t)",
            "u(0, t) = 1 + t",
            1.0,
            0.0,
            -2.0,
        ),
    ];
    for (equation, datum, c2, damping, mass) in cases {
        let s = any(&format!(
            "pdsolve({equation}, u(x, t), list({datum}, {quiet}))"
        ));
        let text = expand_integrals(&rhs(&s), 200);
        let u = |x: f64, t: f64| value(&text, &[("x", x), ("t", t)]);
        let (x0, t0, h) = (0.6_f64, 2.3_f64, 2e-3_f64);
        let utt = (u(x0, t0 + h) - 2.0 * u(x0, t0) + u(x0, t0 - h)) / (h * h);
        let ut = (u(x0, t0 + h) - u(x0, t0 - h)) / (2.0 * h);
        let uxx = (u(x0 + h, t0) - 2.0 * u(x0, t0) + u(x0 - h, t0)) / (h * h);
        close(
            utt + damping * ut + mass * u(x0, t0) - c2 * uxx,
            0.0,
            2e-2 * (1.0 + utt.abs()),
            &s,
        );
        assert!(u(x0, t0).abs() > 1e-3, "{s}");
    }
}

#[test]
fn poisson_integrals_on_disks_and_balls() {
    // Disk: g = cos² φ extends to 1/2 + (r/R)² cos(2θ)/2.
    let s = run("pde_poisson_integral(cos(phi)^2, 2, 2)");
    let text = expand_integrals(&s, 80);
    let (r, theta) = (1.1_f64, 0.7_f64);
    close(
        value(&text, &[("r", r), ("theta", theta)]),
        0.5 + (r / 2.0).powi(2) * (2.0 * theta).cos() / 2.0,
        1e-8,
        &s,
    );
    // Ball: g = cos² θ extends to 1/3 + (2/3)(r/R)² P₂(cos θ).
    let s = run("pde_poisson_integral(cos(theta)^2, 2, 3)");
    let text = expand_integrals(&s, 64);
    let (r, theta, phi) = (0.9_f64, 0.8_f64, 0.3_f64);
    let p2 = 0.5 * (3.0 * theta.cos().powi(2) - 1.0);
    close(
        value(&text, &[("r", r), ("theta", theta), ("phi", phi)]),
        1.0 / 3.0 + 2.0 / 3.0 * (r / 2.0).powi(2) * p2,
        3e-4,
        &s,
    );
}

#[test]
fn coulomb_bound_states_and_spheroidal_eigenvalues() {
    use super::coulomb::coulomb_eigen;
    use super::coulomb::coulomb_mode;
    use super::coulomb::spheroidal_eigen;
    use super::coulomb::spheroidal_mode;
    // Hydrogen-like spectrum and orthonormality.
    close(
        coulomb_eigen(&[2.0, 1.0, 3.0]),
        -4.0 / 9.0,
        1e-14,
        "eigenvalue",
    );
    for (l, n, m) in [
        (0.0, 1.0, 1.0),
        (0.0, 2.0, 2.0),
        (1.0, 2.0, 3.0),
        (1.0, 3.0, 3.0),
        (0.0, 3.0, 1.0),
    ] {
        let inner = simpson(
            &|r| coulomb_mode(&[1.5, l, n, r]) * coulomb_mode(&[1.5, l, m, r]),
            0.0,
            90.0,
        );
        let want = if (n - m).abs() < 1e-12 {
            1.0
        } else {
            0.0
        };
        close(inner, want, 3e-5, "orthonormality");
    }
    // The radial equation -u'' + (l(l+1)/r² - 2Z/r) u = E u.
    let (z, l, n, r, h) = (1.5_f64, 2.0_f64, 4.0_f64, 3.1_f64, 1e-4_f64);
    let u = |r: f64| coulomb_mode(&[z, l, n, r]);
    let second = (u(r + h) - 2.0 * u(r) + u(r - h)) / (h * h);
    close(
        -second + (l * (l + 1.0) / (r * r) - 2.0 * z / r) * u(r),
        coulomb_eigen(&[z, l, n]) * u(r),
        1e-5,
        "Coulomb equation",
    );
    // Evolution problems: data in the span of bound states.
    let bc = "u(0, t) = 0";
    let data = "2*r*exp(-r) + 0.5*r*(2 - r)*exp(-r/2)/sqrt(8)";
    let s = any(&format!(
        "pdsolve(diff(u(r, t), t) = diff(diff(u(r, t), r), r) + 2/r*u(r, t), u(r, t), list({bc}, u(r, 0) = {data}))"
    ));
    let (r, t) = (1.7_f64, 0.3_f64);
    let exact = 2.0 * r * (-r).exp() * t.exp()
        + 0.5 * (0.25 * t).exp() * r * (2.0 - r) * (-r / 2.0).exp() / 8.0_f64.sqrt();
    close(value(&rhs(&s), &[("r", r), ("t", t)]), exact, 1e-6, &s);
    let s = any(&format!(
        "pdsolve(diff(u(r, t), t) = diff(diff(u(r, t), r), r) + (2/r - 2/r^2)*u(r, t), u(r, t), list({bc}, u(r, 0) = r^2*exp(-r/2)))"
    ));
    close(
        value(&rhs(&s), &[("r", r), ("t", t)]),
        r * r * (-r / 2.0).exp() * (0.25 * t).exp(),
        1e-6,
        &s,
    );
    let s = any(&format!(
        "pdsolve(diff(diff(u(r, t), t), t) = diff(diff(u(r, t), r), r) + 2/r*u(r, t), u(r, t), list({bc}, u(r, 0) = 2*r*exp(-r), at(diff(u(r, t), t), t, 0) = 0))"
    ));
    close(
        value(&rhs(&s), &[("r", r), ("t", t)]),
        2.0 * r * (-r).exp() * t.cosh(),
        1e-6,
        &s,
    );
    // Data outside the bound-state span need the scattering states (see the continuum test).
    let s = any(&format!(
        "pdsolve(diff(u(r, t), t) = diff(diff(u(r, t), r), r) + 2/r*u(r, t), u(r, t), list({bc}, u(r, 0) = r*exp(-3*r)))"
    ));
    assert!(s.contains("coulomb_wave"), "{s}");
    // Spheroidal angular functions: gamma = 0 gives n(n+1); first-order perturbation theory.
    for n in 0..4 {
        close(
            spheroidal_eigen(&[0.0, f64::from(n), 0.0]),
            f64::from(n * (n + 1)),
            1e-10,
            "Legendre limit",
        );
    }
    let gamma = 1e-3_f64;
    for n in 0..4 {
        let nf = f64::from(n);
        let first = gamma * (2.0 * nf * (nf + 1.0) - 1.0) / ((2.0 * nf - 1.0) * (2.0 * nf + 3.0));
        close(
            spheroidal_eigen(&[0.0, nf, gamma]) - nf * (nf + 1.0),
            first,
            5e-6,
            "prolate perturbation",
        );
    }
    // The mode solves (1 - x²) y'' - 2 x y' + (λ - γ x² - m²/(1 - x²)) y = 0 and is normalised.
    let (m, n, gamma) = (1.0_f64, 3.0_f64, 2.0_f64);
    let lambda = spheroidal_eigen(&[m, n, gamma]);
    let y = |x: f64| spheroidal_mode(&[m, n, gamma, x]);
    let (x, h) = (0.3_f64, 1e-4_f64);
    let y2 = (y(x + h) - 2.0 * y(x) + y(x - h)) / (h * h);
    let y1 = (y(x + h) - y(x - h)) / (2.0 * h);
    close(
        (1.0 - x * x) * y2 - 2.0 * x * y1 + (lambda - gamma * x * x - m * m / (1.0 - x * x)) * y(x),
        0.0,
        5e-4,
        "spheroidal equation",
    );
    close(
        simpson(&|x| y(x) * y(x), -1.0, 1.0),
        1.0,
        1e-5,
        "spheroidal norm",
    );
    // Oblate case: negative gamma lowers the eigenvalues.
    assert!(spheroidal_eigen(&[0.0, 2.0, -3.0]) < 6.0);
}

#[test]
fn second_order_and_derivative_coupled_jordan_systems() {
    let numeric = |text: &str, x: f64, t: f64| {
        crate::rules::testing::numeric(&crate::rules::standard(), text, &[("x", x), ("t", t)], 1e-9)
            .0
    };
    // u_tt = u_xx + v_xx, v_tt = v_xx on the line: v = (g(x+t) + g(x-t))/2, u = t (g'(x+t) - g'(x-t))/4.
    let s = any(
        "pdsolve(list(diff(diff(u(x, t), t), t) = diff(diff(u(x, t), x), x) + diff(diff(v(x, t), x), x), diff(diff(v(x, t), t), t) = diff(diff(v(x, t), x), x)), list(u(x, t), v(x, t)), list(u(x, 0) = 0, at(diff(u(x, t), t), t, 0) = 0, v(x, 0) = exp(-x^2), at(diff(v(x, t), t), t, 0) = 0))",
    );
    let (a, b) = pair(&s);
    let (x, t) = (0.4_f64, 0.9_f64);
    let g = |z: f64| (-z * z).exp();
    let dg = |z: f64| -2.0 * z * (-z * z).exp();
    close(
        numeric(&b, x, t),
        f64::midpoint(g(x + t), g(x - t)),
        1e-8,
        &s,
    );
    close(
        numeric(&a, x, t),
        t * (dg(x + t) - dg(x - t)) / 4.0,
        1e-8,
        &s,
    );
    // A derivative coupling on a bounded interval: u_t = u_xx + v_x, v_t = v_xx, Dirichlet ends, v(x, 0) = sin 2x.
    let s = any(
        "pdsolve(list(diff(u(x, t), t) = diff(diff(u(x, t), x), x) + diff(v(x, t), x), diff(v(x, t), t) = diff(diff(v(x, t), x), x)), list(u(x, t), v(x, t)), list(u(x, 0) = 0, v(x, 0) = sin(2*x), u(0, t) = 0, u(pi, t) = 0, v(0, t) = 0, v(pi, t) = 0))",
    );
    let (a, b) = pair(&s);
    let (x, t) = (1.1_f64, 0.3_f64);
    close(
        numeric(&b, x, t),
        (-4.0 * t).exp() * (2.0 * x).sin(),
        1e-8,
        &s,
    );
    let reference: f64 = (1..=400)
        .filter(|&n| n != 2)
        .map(|n| {
            let nf = f64::from(n);
            let c = 2.0 / PI * simpson(&|y| 2.0 * (2.0 * y).cos() * (nf * y).sin(), 0.0, PI);
            c * ((-4.0 * t).exp() - (-nf * nf * t).exp()) / (nf * nf - 4.0) * (nf * x).sin()
        })
        .sum();
    let text = expand_integrals(&expand_sums(&a, 150), 60);
    close(value(&text, &[("x", x), ("t", t)]), reference, 2e-4, &s);
}

#[test]
fn numeric_gelfand_levitan_marchenko_for_general_data() {
    let evaluate = |datum: &str, x: f64, t: f64| {
        value(&run(&format!("pde_kdv_scattering({datum}, {x}, {t})")), &[])
    };
    // Reflectionless data: one soliton, u = -2 sech²(x - 4t).
    let (x, t) = (0.3_f64, 0.1_f64);
    close(
        evaluate("-2/cosh(x)^2", x, t),
        -2.0 / (x - 4.0 * t).cosh().powi(2),
        2e-5,
        "one soliton",
    );
    // Two solitons against the closed determinant formula.
    let closed = run("pde_kdv_reflectionless(2, 1)");
    close(
        evaluate("-6/cosh(x)^2", x, 0.05),
        value(&closed, &[("x", x), ("t", 0.05)]),
        2e-4,
        "two solitons",
    );
    // Data with a continuous spectrum: the datum is reproduced at t = 0 ...
    let gaussian = "-3*exp(-x^2)";
    for x in [-0.8_f64, 0.4, 1.3] {
        close(
            evaluate(gaussian, x, 0.0),
            -3.0 * (-x * x).exp(),
            2e-4,
            "initial datum",
        );
    }
    // An off-centre datum (the reflection coefficient carries the phase of the shift).
    for x in [-0.5_f64, 1.2] {
        close(
            evaluate("-3*exp(-(x - 0.7)^2)", x, 0.0),
            -3.0 * (-(x - 0.7) * (x - 0.7)).exp(),
            2e-4,
            "shifted datum",
        );
    }
    // ... and the evolved field solves u_t - 6 u u_x + u_xxx = 0.
    let (x, t, h) = (0.2_f64, 0.05_f64, 0.05_f64);
    let u = |x: f64, t: f64| evaluate(gaussian, x, t);
    let u_t = (u(x, t + 0.02) - u(x, t - 0.02)) / 0.04;
    let u_x = (u(x + h, t) - u(x - h, t)) / (2.0 * h);
    let u_xxx = (u(x + 2.0 * h, t) - 2.0 * u(x + h, t) + 2.0 * u(x - h, t) - u(x - 2.0 * h, t))
        / (2.0 * h.powi(3));
    let residual = u_t - 6.0 * u(x, t) * u_x + u_xxx;
    // The third difference carries an O(h²) error, so the tolerance is loose.
    close(residual, 0.0, 0.1 * (1.0 + u_xxx.abs()), "KdV residual");
}

#[test]
fn heat_outside_a_ball_without_axial_symmetry() {
    let bounds = |text: &str| {
        text.replace(", rho, 1, oo)", ", rho, 1, 9)")
            .replace(", k, 0, oo)", ", k, 0, 8)")
    };
    let lap = "diff(diff(u(r, th, ph, t), r), r) + 2*diff(u(r, th, ph, t), r)/r + diff(diff(u(r, th, ph, t), th), th)/r^2 + cos(th)/sin(th)*diff(u(r, th, ph, t), th)/r^2 + diff(diff(u(r, th, ph, t), ph), ph)/(r^2*sin(th)^2)";
    let radial = |q: f64| (-(q - 4.0) * (q - 4.0)).exp();
    let (r, th, ph, t) = (2.0_f64, 0.9_f64, 0.7_f64, 0.4_f64);
    // l = 1, m = 1: sin θ cos φ.
    let s = any(&format!(
        "pdsolve(diff(u(r, th, ph, t), t) = {lap}, u(r, th, ph, t), list(u(1, th, ph, t) = 0, u(oo, th, ph, t) = 0, u(r, th, ph, 0) = exp(-(r - 4)^2)*sin(th)*cos(ph)))"
    ));
    let text = expand_integrals(&bounds(&rhs(&s)), 90);
    let reference = exterior_heat_reference(2.0, 2.0, &radial, t, r) * th.sin() * ph.cos();
    close(
        value(&text, &[("r", r), ("th", th), ("ph", ph), ("t", t)]),
        reference,
        4e-3,
        &s,
    );
    // l = 2, m = 2: sin² θ sin 2φ.
    let s = any(&format!(
        "pdsolve(diff(u(r, th, ph, t), t) = {lap}, u(r, th, ph, t), list(u(1, th, ph, t) = 0, u(oo, th, ph, t) = 0, u(r, th, ph, 0) = exp(-(r - 4)^2)*sin(th)^2*sin(2*ph)))"
    ));
    let text = expand_integrals(&bounds(&rhs(&s)), 90);
    let reference =
        exterior_heat_reference(2.0, 6.0, &radial, t, r) * th.sin().powi(2) * (2.0 * ph).sin();
    close(
        value(&text, &[("r", r), ("th", th), ("ph", ph), ("t", t)]),
        reference,
        4e-3,
        &s,
    );
}

#[test]
fn dirichlet_greens_functions_of_disks_and_balls() {
    let w = |text: &str, point: &[(&str, f64)]| value(text, point);
    for (d, names) in [(2, vec!["x", "y"]), (3, vec!["x", "y", "z"])] {
        let s = run(&format!("pde_greens_ball(2, {d})"));
        let at = |x: &[f64], y: &[f64]| {
            let mut bindings: Vec<(String, f64)> = Vec::new();
            for (k, n) in names.iter().enumerate() {
                bindings.push(((*n).to_owned(), x[k]));
                bindings.push((format!("s_{n}"), y[k]));
            }
            let refs: Vec<(&str, f64)> = bindings.iter().map(|(a, b)| (a.as_str(), *b)).collect();
            w(&s, &refs)
        };
        let source: Vec<f64> = if d == 2 {
            vec![0.7, -0.4]
        } else {
            vec![0.7, -0.4, 0.5]
        };
        // Zero on the boundary |x| = R = 2.
        for angle in [0.3_f64, 1.9, 4.0] {
            let boundary: Vec<f64> = if d == 2 {
                vec![2.0 * angle.cos(), 2.0 * angle.sin()]
            } else {
                vec![2.0 * angle.cos() * 0.6, 2.0 * angle.sin() * 0.6, 2.0 * 0.8]
            };
            close(at(&boundary, &source), 0.0, 1e-12, &s);
        }
        // Symmetric, and harmonic away from the source.
        let x: Vec<f64> = if d == 2 {
            vec![-0.9, 0.8]
        } else {
            vec![-0.9, 0.8, -0.3]
        };
        close(at(&x, &source), at(&source, &x), 1e-12, "symmetry");
        let h = 1e-3_f64;
        let mut laplacian = -2.0 * f64::from(d) * at(&x, &source);
        for k in 0..usize::try_from(d).unwrap_or(2) {
            for sign in [-1.0, 1.0] {
                let mut moved = x.clone();
                moved[k] += sign * h;
                laplacian += at(&moved, &source);
            }
        }
        close(laplacian / (h * h), 0.0, 2e-3, "harmonic");
    }
    // Exterior points (|x|, |y| > R) vanish on the sphere as well.
    let s = run("pde_greens_ball(1, 3)");
    close(
        w(
            &s,
            &[
                ("x", 1.0),
                ("y", 0.0),
                ("z", 0.0),
                ("s_x", 3.0),
                ("s_y", 1.0),
                ("s_z", -2.0),
            ],
        ),
        0.0,
        1e-12,
        &s,
    );
}

#[test]
fn hodograph_linearises_diagonal_quasilinear_systems() {
    // r_t + r r_x = 0, s_t + s s_x = 0: x = r t + F(r) = s t + G(s), so t = (G(s) - F(r))/(r - s).
    let s = run("pde_hodograph(r, s)");
    let (equation, _) = s.split_once(" = 0").unwrap_or((&s, ""));
    let explicit = "((sin(s) - r^2)/(r - s))";
    let text = run(&equation.replace("t(r, s)", explicit));
    close(value(&text, &[("r", 0.7), ("s", -0.4)]), 0.0, 1e-9, &text);
    // A non-trivial pair is returned as an equation in the derivatives of t.
    let s = run("pde_hodograph(r + s, r*s)");
    assert!(s.contains("diff") && s.contains("= 0"), "{s}");
}

#[test]
fn waves_outside_a_ball_without_axial_symmetry() {
    let bounds = |text: &str| {
        text.replace(", rho, 1, oo)", ", rho, 1, 9)")
            .replace(", k, 0, oo)", ", k, 0, 8)")
    };
    let lap = "diff(diff(u(r, th, ph, t), r), r) + 2*diff(u(r, th, ph, t), r)/r + diff(diff(u(r, th, ph, t), th), th)/r^2 + cos(th)/sin(th)*diff(u(r, th, ph, t), th)/r^2 + diff(diff(u(r, th, ph, t), ph), ph)/(r^2*sin(th)^2)";
    let (r, th, ph, t) = (3.0_f64, 0.9_f64, 0.7_f64, 1.0_f64);
    let s = any(&format!(
        "pdsolve(diff(diff(u(r, th, ph, t), t), t) = {lap}, u(r, th, ph, t), list(u(1, th, ph, t) = 0, u(oo, th, ph, t) = 0, u(r, th, ph, 0) = exp(-2*(r - 4)^2)*sin(th)*cos(ph), at(diff(u(r, th, ph, t), t), t, 0) = 0))"
    ));
    let text = expand_integrals(&bounds(&rhs(&s)), 90);
    let reference =
        exterior_wave_reference(2.0, 2.0, &|q| (-2.0 * (q - 4.0) * (q - 4.0)).exp(), t, r)
            * th.sin()
            * ph.cos();
    close(
        value(&text, &[("r", r), ("th", th), ("ph", ph), ("t", t)]),
        reference,
        6e-3,
        &s,
    );
}

#[test]
fn laplace_route_with_initial_data_and_parabolic_robin_signals() {
    // Non-quiet start: the wave equation with a boundary signal and a pulse.
    let s = any(
        "pdsolve(diff(diff(u(x, t), t), t) = diff(diff(u(x, t), x), x), u(x, t), list(u(0, t) = sin(t), u(x, 0) = exp(-(x - 3)^2), at(diff(u(x, t), t), t, 0) = 0))",
    );
    let (x, t) = (1.2_f64, 2.0_f64);
    let odd = |z: f64| {
        if z >= 0.0 {
            (-(z - 3.0) * (z - 3.0)).exp()
        } else {
            -(-(-z - 3.0) * (-z - 3.0)).exp()
        }
    };
    let want = f64::midpoint(odd(x + t), odd(x - t)) + if t > x { (t - x).sin() } else { 0.0 };
    close(value(&rhs(&s), &[("x", x), ("t", t)]), want, 1e-9, &s);
    // Robin signals for the heat equation: residual, boundary condition and initial value.
    for (equation, reaction, h, datum, initial) in [
        (
            "diff(u(x, t), t) = diff(diff(u(x, t), x), x)",
            0.0_f64,
            2.0_f64,
            "1",
            None,
        ),
        (
            "diff(u(x, t), t) = diff(diff(u(x, t), x), x) - u(x, t)",
            -1.0,
            1.0,
            "sin(t)",
            Some("exp(-x)"),
        ),
    ] {
        let mut conditions = format!("at(diff(u(x, t), x), x, 0) + {h}*u(0, t) = {datum}");
        match initial {
            | Some(f) => conditions = format!("{conditions}, u(x, 0) = {f}"),
            | None => conditions.push_str(", u(x, 0) = 0"),
        }
        let s = any(&format!("pdsolve({equation}, u(x, t), list({conditions}))"));
        let text = expand_integrals(&rhs(&s).replace(", s, 0, oo)", ", s, 0, 25)"), 400);
        let u = |x: f64, t: f64| value(&text, &[("x", x), ("t", t)]);
        let (x0, t0, d) = (0.8_f64, 0.5_f64, 1e-3_f64);
        let residual = (u(x0, t0 + d) - u(x0, t0 - d)) / (2.0 * d)
            - (u(x0 + d, t0) - 2.0 * u(x0, t0) + u(x0 - d, t0)) / (d * d)
            - reaction * u(x0, t0);
        close(residual, 0.0, 3e-3, &s);
        let e = 0.03_f64;
        let g = if datum == "1" {
            1.0
        } else {
            t0.sin()
        };
        let slope = (-3.0 * u(0.0, t0) + 4.0 * u(e, t0) - u(2.0 * e, t0)) / (2.0 * e);
        close(slope + h * u(0.0, t0), g, 3e-2, "Robin condition");
    }
}

#[test]
fn crum_determinants_and_hodograph_solutions() {
    let v = |text: &str, x: f64, t: f64| value(text, &[("x", x), ("t", t)]);
    // Two-fold Darboux transformation of the vacuum: a two-soliton.
    let s = run("pde_crum_kdv(0, list(cosh(x - 4*t), sinh(2*x - 32*t)))");
    let u = |x: f64, t: f64| v(&s, x, t);
    let (x, t, h) = (0.35_f64, 0.05_f64, 2e-3_f64);
    let u_t = (u(x, t + h) - u(x, t - h)) / (2.0 * h);
    let u_x = (u(x + h, t) - u(x - h, t)) / (2.0 * h);
    let u_xxx = (u(x + 2.0 * h, t) - 2.0 * u(x + h, t) + 2.0 * u(x - h, t) - u(x - 2.0 * h, t))
        / (2.0 * h.powi(3));
    close(
        u_t - 6.0 * u(x, t) * u_x + u_xxx,
        0.0,
        5e-3 * (1.0 + u_xxx.abs()),
        &s,
    );
    assert!(u(x, t) < -0.1, "{s}");
    // One step agrees with the single Darboux transformation.
    let one = run("pde_crum_kdv(0, list(cosh(x - 4*t)))");
    close(
        v(&one, 0.8, 0.3),
        v(&run("pde_darboux_kdv(0, cosh(x - 4*t))"), 0.8, 0.3),
        1e-12,
        &one,
    );
    // The hodograph solution satisfies x = r t + F(r) = s t + G(s) and the linear equation.
    let s = run("pde_hodograph_solution(r^2, sin(s))");
    let (xs, ts) = {
        let inner = s
            .trim()
            .strip_prefix("list(")
            .and_then(|x| x.strip_suffix(')'))
            .unwrap_or("");
        let mut depth = 0_i32;
        let mut split = 0;
        for (i, c) in inner.char_indices() {
            match c {
                | '(' => depth += 1,
                | ')' => depth -= 1,
                | ',' if depth == 0 => split = i,
                | _ => {},
            }
        }
        (
            inner[..split].trim().to_owned(),
            inner[split + 1..].trim().to_owned(),
        )
    };
    let (r, sv) = (0.7_f64, -0.4_f64);
    let (x_value, t_value) = (
        value(&xs, &[("r", r), ("s", sv)]),
        value(&ts, &[("r", r), ("s", sv)]),
    );
    close(x_value, r * (t_value + r), 1e-12, &s);
    close(x_value, sv * t_value + sv.sin(), 1e-12, &s);
}

/// Splits `list(x, t)` printed by the engine.
fn split_pair(text: &str) -> (String, String) {
    let inner = text
        .trim()
        .strip_prefix("list(")
        .and_then(|x| x.strip_suffix(')'))
        .unwrap_or("");
    let mut depth = 0_i32;
    let mut split = 0;
    for (i, c) in inner.char_indices() {
        match c {
            | '(' => depth += 1,
            | ')' => depth -= 1,
            | ',' if depth == 0 => split = i,
            | _ => {},
        }
    }
    (
        inner[..split].trim().to_owned(),
        inner[split + 1..].trim().to_owned(),
    )
}

#[test]
fn hodograph_solutions_of_decoupled_affine_systems() {
    // r_t + (2 r + 1) r_x = 0, s_t + (3 s - 2) s_x = 0 with F = r^2, G = sin s:
    // x = (2 r + 1) t + r^2 = (3 s - 2) t + sin s.
    let s = run("pde_hodograph_solution(r^2, sin(s), 2*r + 1, 3*s - 2)");
    let (xs, ts) = split_pair(&s);
    let at = |text: &str, r: f64, sv: f64| value(text, &[("r", r), ("s", sv)]);
    let x = |r: f64, sv: f64| at(&xs, r, sv);
    let t = |r: f64, sv: f64| at(&ts, r, sv);
    for (r, sv) in [(0.7_f64, -0.4_f64), (1.9, 0.3), (-0.6, 2.2)] {
        let lambda = 2.0 * r + 1.0;
        let mu = 3.0 * sv - 2.0;
        close(x(r, sv), lambda * t(r, sv) + r * r, 1e-11, &s);
        close(x(r, sv), mu * t(r, sv) + sv.sin(), 1e-11, &s);
        // The hodograph equations x_s = lambda t_s and x_r = mu t_r (central differences).
        let h = 1e-5;
        let d = |f: &dyn Fn(f64, f64) -> f64, wrt_r: bool| {
            if wrt_r {
                (f(r + h, sv) - f(r - h, sv)) / (2.0 * h)
            } else {
                (f(r, sv + h) - f(r, sv - h)) / (2.0 * h)
            }
        };
        close(
            d(&x, false),
            lambda * d(&t, false),
            2e-6 * (1.0 + d(&x, false).abs()),
            "x_s = lambda t_s",
        );
        close(
            d(&x, true),
            mu * d(&t, true),
            2e-6 * (1.0 + d(&x, true).abs()),
            "x_r = mu t_r",
        );
    }
    // Non-affine decoupled speeds work the same way.
    let s = run("pde_hodograph_solution(r, s^2, r^2, cos(s))");
    let (xs, ts) = split_pair(&s);
    let (r, sv) = (0.8_f64, 0.5_f64);
    close(at(&xs, r, sv), r * r * at(&ts, r, sv) + r, 1e-11, &s);
    close(
        at(&xs, r, sv),
        sv.cos() * at(&ts, r, sv) + sv * sv,
        1e-11,
        &s,
    );
    // The two-argument form is unchanged (lambda = r, mu = s).
    let default = run("pde_hodograph_solution(r^2, sin(s))");
    let explicit = run("pde_hodograph_solution(r^2, sin(s), r, s)");
    assert_eq!(default, explicit);
    // Coupled systems outside the integer Euler-Poisson-Darboux family stay symbolic.
    for request in [
        "pde_hodograph_solution(r^2, sin(s), r + s, r + 2*s)",
        "pde_hodograph_solution(r^2, sin(s), 3*r + s, r + 3*s)",
        "pde_hodograph_solution(r^2, sin(s), 2*r + s, r + 2*s + 1)",
    ] {
        let text = any(request);
        assert!(text.contains("pde_hodograph_solution"), "{request}: {text}");
    }
}

#[test]
fn hodograph_solutions_of_coupled_affine_systems_of_integer_order() {
    // lambda = N k r + (N - 1) k s + c, mu = (N - 1) k r + N k s + c, k = 1/2, 1, 2:
    // the general solution of the E(N, N) hodograph equation with F, G arbitrary.
    let at = |text: &str, r: f64, sv: f64| value(text, &[("r", r), ("s", sv)]);
    for (n, k, c, speeds) in [
        (
            2_u32,
            1.0_f64,
            0.0_f64,
            "pde_hodograph_solution(r^3, sin(s), 2*r + s, r + 2*s)",
        ),
        (
            2,
            0.5,
            0.3,
            "pde_hodograph_solution(exp(r), s^2, r + s/2 + 3/10, r/2 + s + 3/10)",
        ),
        (
            3,
            1.0,
            -0.2,
            "pde_hodograph_solution(sin(r), cos(s), 3*r + 2*s - 1/5, 2*r + 3*s - 1/5)",
        ),
        (
            4,
            2.0,
            0.0,
            "pde_hodograph_solution(r^2, s^3, 8*r + 6*s, 6*r + 8*s)",
        ),
    ] {
        let s = run(speeds);
        let (xs, ts) = split_pair(&s);
        let x = |r: f64, sv: f64| at(&xs, r, sv);
        let t = |r: f64, sv: f64| at(&ts, r, sv);
        let nf = f64::from(n);
        let lambda = |r: f64, sv: f64| nf * k * r + (nf - 1.0) * k * sv + c;
        let mu = |r: f64, sv: f64| (nf - 1.0) * k * r + nf * k * sv + c;
        for (r, sv) in [(0.9_f64, -0.35_f64), (1.4, 0.6), (-0.5, 1.2)] {
            // The hodograph equations x_s = lambda t_s and x_r = mu t_r.
            let h = 2e-4;
            let ds = |f: &dyn Fn(f64, f64) -> f64| (f(r, sv + h) - f(r, sv - h)) / (2.0 * h);
            let dr = |f: &dyn Fn(f64, f64) -> f64| (f(r + h, sv) - f(r - h, sv)) / (2.0 * h);
            let scale = 1.0 + dr(&x).abs().max(ds(&x).abs());
            close(
                ds(&x),
                lambda(r, sv) * ds(&t),
                1e-5 * scale,
                &format!("N = {n}: x_s = lambda t_s at ({r}, {sv})"),
            );
            close(
                dr(&x),
                mu(r, sv) * dr(&t),
                1e-5 * scale,
                &format!("N = {n}: x_r = mu t_r at ({r}, {sv})"),
            );
            // The linear hodograph equation (lambda - mu) t_rs + lambda_r t_s - mu_s t_r = 0.
            let t_rs = (t(r + h, sv + h) - t(r + h, sv - h) - t(r - h, sv + h) + t(r - h, sv - h))
                / (4.0 * h * h);
            let residual = (lambda(r, sv) - mu(r, sv)) * t_rs + nf * k * ds(&t) - nf * k * dr(&t);
            close(
                residual,
                0.0,
                2e-3 * (1.0 + t_rs.abs() + ds(&t).abs() + dr(&t).abs()),
                &format!("N = {n}: EPD at ({r}, {sv})"),
            );
        }
    }
}

#[test]
fn robin_and_neumann_greens_functions() {
    let w = |text: &str, point: &[(&str, f64)]| value(text, point);
    // Half-space z > 0, ∂_n u + h u = 0 with ∂_n = -∂_z: u_z = h u at z = 0.
    for (d, names, source) in [
        (3, vec!["x", "y", "z"], vec![0.4, -0.3, 0.7]),
        (2, vec!["x", "y"], vec![0.4, 0.7]),
    ] {
        for h in [0.0_f64, 1.5] {
            let s = run(&format!("pde_greens_robin({d}, {h})"));
            let text = expand_integrals(&s.replace(", w, 0, oo)", ", w, 0, 30)"), 600);
            let at = |x: &[f64]| {
                let mut bindings: Vec<(String, f64)> = Vec::new();
                for (k, n) in names.iter().enumerate() {
                    bindings.push(((*n).to_owned(), x[k]));
                    bindings.push((format!("s_{n}"), source[k]));
                }
                let refs: Vec<(&str, f64)> =
                    bindings.iter().map(|(a, b)| (a.as_str(), *b)).collect();
                w(&text, &refs)
            };
            let last = names.len() - 1;
            let mut base: Vec<f64> = if d == 3 {
                vec![-0.9, 0.8, 0.0]
            } else {
                vec![-0.9, 0.0]
            };
            let e = 1e-4_f64;
            let g0 = at(&base);
            base[last] = e;
            let g1 = at(&base);
            close((g1 - g0) / e - h * g0, 0.0, 2e-3, &s);
            // Harmonic away from the source.
            let x: Vec<f64> = if d == 3 {
                vec![-0.9, 0.8, 0.6]
            } else {
                vec![-0.9, 0.6]
            };
            let step = 1e-2_f64;
            let mut laplacian = -2.0 * f64::from(d) * at(&x);
            for k in 0..names.len() {
                for sign in [-1.0, 1.0] {
                    let mut moved = x.clone();
                    moved[k] += sign * step;
                    laplacian += at(&moved);
                }
            }
            close(laplacian / (step * step), 0.0, 2e-2, "harmonic");
        }
    }
    // The interior Robin Green's function of the ball of radius 2.
    let s = run("pde_greens_ball_robin(2, 1)");
    let text = expand_sums(&s, 70);
    let source = [0.5_f64, -0.4, 0.3];
    let at = |x: [f64; 3]| {
        w(
            &text,
            &[
                ("x", x[0]),
                ("y", x[1]),
                ("z", x[2]),
                ("s_x", source[0]),
                ("s_y", source[1]),
                ("s_z", source[2]),
            ],
        )
    };
    let n = [0.6_f64, 0.0, 0.8];
    let delta = 1e-4_f64;
    let (outer, inner) = (
        at([2.0 * n[0], 2.0 * n[1], 2.0 * n[2]]),
        at([
            (2.0 - delta) * n[0],
            (2.0 - delta) * n[1],
            (2.0 - delta) * n[2],
        ]),
    );
    close((outer - inner) / delta + outer, 0.0, 2e-3, &s);
    let x = [-0.6_f64, 0.7, 0.2];
    let step = 1e-2_f64;
    let mut laplacian = -6.0 * at(x);
    for k in 0..3 {
        for sign in [-1.0, 1.0] {
            let mut moved = x;
            moved[k] += sign * step;
            laplacian += at(moved);
        }
    }
    close(laplacian / (step * step), 0.0, 2e-2, "harmonic ball");
    close(
        at(x),
        w(
            &text,
            &[
                ("x", source[0]),
                ("y", source[1]),
                ("z", source[2]),
                ("s_x", x[0]),
                ("s_y", x[1]),
                ("s_z", x[2]),
            ],
        ),
        1e-9,
        "symmetry",
    );
}

/// `text` compiled once: a closure evaluating it at real bindings. The
/// per-call `value` installs the rules and reparses every time, which is too
/// slow for the quadratures below.
fn compiled(text: &str) -> impl Fn(&[(&str, f64)]) -> f64 {
    let sets = crate::rules::standard();
    let mut graph = Graph::new();
    assert!(Engine::install(&mut graph, &sets).is_ok());
    let root = graph
        .parse(text)
        .unwrap_or_else(|e| panic!("cannot parse `{text}`: {e}"));
    let mut symbols = std::collections::HashMap::new();
    for name in ["x", "y", "z", "s_x", "s_y", "s_z"] {
        symbols.insert(name, graph.interner_mut().symbol(name));
    }
    move |at: &[(&str, f64)]| {
        let mut env = Env::numeric(0.0);
        for (name, value) in at {
            if let Some(&symbol) = symbols.get(*name) {
                env.bind(symbol, *value);
            }
        }
        graph.eval(root, &env).unwrap_or(f64::NAN)
    }
}

/// The second-order central Laplacian with step `h`, for a function of the
/// Cartesian coordinates `x`.
fn finite_laplacian(
    f: &dyn Fn(&[f64]) -> f64,
    x: &[f64],
    h: f64,
) -> f64 {
    let weight = if x.len() == 3 { -6.0 } else { -4.0 };
    let mut total = weight * f(x);
    for k in 0..x.len() {
        for sign in [-1.0, 1.0] {
            let mut moved = x.to_vec();
            moved[k] += sign * h;
            total += f(&moved);
        }
    }
    total / (h * h)
}

/// `finite_laplacian` with Richardson extrapolation (`O(h⁴)`).
fn richardson_laplacian(
    f: &dyn Fn(&[f64]) -> f64,
    x: &[f64],
) -> f64 {
    let coarse = finite_laplacian(f, x, 2e-3);
    let fine = finite_laplacian(f, x, 1e-3);
    (4.0 * fine - coarse) / 3.0
}

/// `∫_B f dx` over the disk (`d = 2`) or ball (`d = 3`) of radius `radius`,
/// for a source on the last axis at distance `rho`: Simpson in source-centred
/// polar/spherical coordinates with the radial substitution `r' = u² r'(θ)`
/// (which removes the endpoint singularity of the fundamental solution).
#[allow(clippy::single_call_fn)]
fn ball_mean(
    d: usize,
    radius: f64,
    rho: f64,
    f: &dyn Fn(&[f64]) -> f64,
) -> f64 {
    let n = if d == 3 { 400 } else { 800 };
    let theta_max = if d == 3 { PI } else { 2.0 * PI };
    let h_theta = theta_max / f64::from(n);
    let h_u = 1.0 / f64::from(n);
    let mut total = 0.0;
    for i in 0..=n {
        let theta = h_theta * f64::from(i);
        let w_theta = if i == 0 || i == n {
            1.0
        } else if i % 2 == 1 {
            4.0
        } else {
            2.0
        };
        let reach = -rho * theta.cos() + (radius * radius - rho * rho * theta.sin().powi(2)).sqrt();
        let mut inner = 0.0;
        for j in 0..=n {
            let u = h_u * f64::from(j);
            let w = if j == 0 || j == n {
                1.0
            } else if j % 2 == 1 {
                4.0
            } else {
                2.0
            };
            let value = if j == 0 {
                0.0
            } else {
                let mut x = vec![0.0; d];
                let t = u * u;
                x[0] = t * reach * theta.sin();
                x[d - 1] = rho + t * reach * theta.cos();
                f(&x)
            };
            let jacobian = if d == 3 {
                2.0 * u.powi(5) * reach.powi(3) * theta.sin()
            } else {
                2.0 * u.powi(3) * reach * reach
            };
            inner += w * jacobian * value;
        }
        total += w_theta * inner;
    }
    let factor = h_theta * h_u / 9.0 * if d == 3 { 2.0 * PI } else { 1.0 };
    total * factor
}

#[test]
fn generalised_neumann_greens_functions_of_disks_and_balls() {
    // The inverse of the Neumann Laplacian on zero-mean functions: symmetric,
    // ∂_r G = 0, ∫_B G = 0 and -ΔG = δ - 1/|B|. The auxiliary
    // T = G - (|x|² - ⟨|x|²⟩)/(2d|B|) carries the task's convention
    // ∂_r T = -1/|∂B| and -ΔT = δ.
    for d in [3_usize, 2] {
        let radius = 2.0_f64;
        let axes: Vec<(&str, &str)> = if d == 3 {
            vec![("x", "s_x"), ("y", "s_y"), ("z", "s_z")]
        } else {
            vec![("x", "s_x"), ("y", "s_y")]
        };
        let s = run(&format!("pde_greens_ball_neumann({radius}, {d})"));
        let g = compiled(&s);
        let at = |x: &[f64], q: &[f64]| {
            let mut bindings: Vec<(&str, f64)> = Vec::new();
            for (k, &(a, b)) in axes.iter().enumerate() {
                bindings.push((a, x[k]));
                bindings.push((b, q[k]));
            }
            g(&bindings)
        };
        let volume = if d == 3 {
            4.0 / 3.0 * PI * radius.powi(3)
        } else {
            PI * radius * radius
        };
        let boundary = if d == 3 {
            4.0 * PI * radius * radius
        } else {
            2.0 * PI * radius
        };
        let mean_square = if d == 3 {
            3.0 * radius * radius / 5.0
        } else {
            radius * radius / 2.0
        };
        let auxiliary = |x: &[f64], q: &[f64]| {
            let xx = x.iter().map(|v| v * v).sum::<f64>();
            at(x, q) - (xx - mean_square) / (2.0 * d as f64 * volume)
        };
        let source: Vec<f64> = if d == 3 {
            vec![0.0, 0.0, 0.9]
        } else {
            vec![0.0, 0.9]
        };
        let normals: Vec<Vec<f64>> = if d == 3 {
            vec![
                vec![1.0, 0.0, 0.0],
                vec![0.6, 0.0, 0.8],
                vec![0.0, 0.6, 0.8],
            ]
        } else {
            vec![vec![1.0, 0.0], vec![0.6, 0.8], vec![0.0, 1.0]]
        };
        let delta = 1e-4_f64;
        for normal in &normals {
            let outer: Vec<f64> = normal.iter().map(|v| (radius + delta) * v).collect();
            let inner: Vec<f64> = normal.iter().map(|v| (radius - delta) * v).collect();
            close(
                (at(&outer, &source) - at(&inner, &source)) / (2.0 * delta),
                0.0,
                1e-7,
                "∂_r G = 0",
            );
            let slope = (auxiliary(&outer, &source) - auxiliary(&inner, &source)) / (2.0 * delta);
            close(slope, -1.0 / boundary, 1e-6 / boundary, "∂_r T = -1/|∂B|");
        }
        // Exact symmetry.
        let first: Vec<f64> = if d == 3 {
            vec![0.4, 0.3, 0.2]
        } else {
            vec![0.4, 0.3]
        };
        let second: Vec<f64> = if d == 3 {
            vec![-0.5, 0.2, 0.7]
        } else {
            vec![-0.5, 0.2]
        };
        close(
            at(&first, &second),
            at(&second, &first),
            1e-10,
            "G(x, s) = G(s, x)",
        );
        close(
            at(&source, &first),
            at(&first, &source),
            1e-10,
            "G(x, s) = G(s, x)",
        );
        // Zero mean, by source-centred quadrature over the whole domain.
        let mean = ball_mean(d, radius, 0.9, &|x| at(x, &source));
        close(mean, 0.0, 1e-7, "∫_B G dx = 0");
        // -ΔG = 1/|B| (and -ΔT = 0) away from the source.
        let point: Vec<f64> = if d == 3 {
            vec![0.4, 0.3, 0.2]
        } else {
            vec![0.4, 0.3]
        };
        close(
            richardson_laplacian(&|x| at(x, &source), &point),
            1.0 / volume,
            1e-6,
            "ΔG = 1/|B|",
        );
        close(
            richardson_laplacian(&|x| auxiliary(x, &source), &point),
            0.0,
            1e-6,
            "ΔT = 0",
        );
    }
}

#[test]
fn exterior_robin_greens_functions_of_disks_and_balls() {
    // (∂_r + h)G = 0 with the outward normal, h > 0 and hR not an integer.
    for d in [3_usize, 2] {
        let radius = 1.0_f64;
        let h = 0.5_f64;
        let axes: Vec<(&str, &str)> = if d == 3 {
            vec![("x", "s_x"), ("y", "s_y"), ("z", "s_z")]
        } else {
            vec![("x", "s_x"), ("y", "s_y")]
        };
        let s = run(&format!("pde_greens_exterior_robin({radius}, 1/2, {d})"));
        let text = expand_sums(&s, 70);
        let g = compiled(&text);
        let at = |x: &[f64], q: &[f64]| {
            let mut bindings: Vec<(&str, f64)> = Vec::new();
            for (k, &(a, b)) in axes.iter().enumerate() {
                bindings.push((a, x[k]));
                bindings.push((b, q[k]));
            }
            g(&bindings)
        };
        let source: Vec<f64> = if d == 3 {
            vec![1.5, 0.3, 0.2]
        } else {
            vec![1.5, 0.3]
        };
        let normals: Vec<Vec<f64>> = if d == 3 {
            vec![
                vec![1.0, 0.0, 0.0],
                vec![0.6, 0.0, 0.8],
                vec![0.0, 0.6, 0.8],
            ]
        } else {
            vec![vec![1.0, 0.0], vec![0.6, 0.8], vec![0.0, 1.0]]
        };
        let delta = 1e-4_f64;
        for normal in &normals {
            let outer: Vec<f64> = normal.iter().map(|v| (radius + delta) * v).collect();
            let inner: Vec<f64> = normal.iter().map(|v| (radius - delta) * v).collect();
            let middle: Vec<f64> = normal.iter().map(|v| radius * v).collect();
            let slope = (at(&outer, &source) - at(&inner, &source)) / (2.0 * delta);
            close(
                slope + h * at(&middle, &source),
                0.0,
                1e-6,
                "∂_r G + h G = 0",
            );
        }
        // Exact symmetry (the Robin problem is self-adjoint).
        let first: Vec<f64> = if d == 3 {
            vec![1.2, -0.5, 0.7]
        } else {
            vec![1.2, -0.5]
        };
        let second: Vec<f64> = if d == 3 {
            vec![1.8, 0.9, -0.4]
        } else {
            vec![1.8, 0.9]
        };
        close(
            at(&first, &second),
            at(&second, &first),
            1e-10,
            "G(x, s) = G(s, x)",
        );
        close(
            at(&source, &first),
            at(&first, &source),
            1e-10,
            "G(x, s) = G(s, x)",
        );
        // Far field: 3D tends to the free-space fundamental solution plus the
        // leading 1/r correction; the 2D kernel is bounded.
        let far: Vec<f64> = if d == 3 {
            vec![0.0, 0.0, 200.0]
        } else {
            vec![0.0, 1e4]
        };
        let r = far.iter().map(|v| v * v).sum::<f64>().sqrt();
        let rho = source.iter().map(|v| v * v).sum::<f64>().sqrt();
        let distance = far
            .iter()
            .zip(&source)
            .map(|(a, b)| (a - b) * (a - b))
            .sum::<f64>()
            .sqrt();
        if d == 3 {
            let fundamental = 1.0 / (4.0 * PI * distance);
            let leading = h * radius * radius / (4.0 * PI * (1.0 - h * radius) * r * rho);
            close(
                at(&far, &source) - fundamental,
                leading,
                1e-5,
                "far-field 1/r correction",
            );
        } else {
            let limit = ((rho / radius).ln() - 1.0 / (radius * h)) / (2.0 * PI);
            close(at(&far, &source), limit, 1e-4, "bounded far field");
        }
        // h → ∞ against the independent exterior Dirichlet operator.
        let big = expand_sums(
            &run(&format!(
                "pde_greens_exterior_robin({radius}, 2000001/2, {d})"
            )),
            70,
        );
        let limit = compiled(&big);
        let dirichlet_text = run(&format!("pde_greens_ball({radius}, {d})"));
        let dirichlet = compiled(&dirichlet_text);
        for point in [vec![2.0, 0.5, 0.3], vec![-1.7, 1.1, 0.6]] {
            let point = &point[..d];
            close(
                limit(&bind(&axes, point, &source)),
                dirichlet(&bind(&axes, point, &source)),
                1e-5,
                "h → ∞ is Dirichlet",
            );
        }
    }
}

/// The bindings of the field point `x` and the source `q` for `axes`.
fn bind<'a>(
    axes: &[(&'a str, &'a str)],
    x: &[f64],
    q: &[f64],
) -> Vec<(&'a str, f64)> {
    let mut bindings = Vec::with_capacity(2 * axes.len());
    for (k, &(a, b)) in axes.iter().enumerate() {
        bindings.push((a, x[k]));
        bindings.push((b, q[k]));
    }
    bindings
}

#[test]
fn greens_kernels_represent_neumann_and_robin_solutions() {
    // Neumann ball of radius 2: g = cos θ on the sphere, u = z.
    let radius = 2.0_f64;
    let s = run("pde_greens_ball_neumann(2, 3)");
    let g = compiled(&s);
    let at = |x: &[f64], q: &[f64]| g(&bind(&[("x", "s_x"), ("y", "s_y"), ("z", "s_z")], x, q));
    let u = |r: f64| {
        let integrand = |theta: f64| {
            let source = [radius * theta.sin(), 0.0, radius * theta.cos()];
            let field = [0.0, 0.0, r];
            at(&field, &source) * theta.cos() * theta.sin()
        };
        simpson(&integrand, 0.0, PI) * 2.0 * PI * radius * radius
    };
    for r in [0.5_f64, 1.0, 1.5] {
        close(u(r), r, 1e-6, "Neumann ball u = z");
    }
    // pdsolve solves the same Neumann ball problem; its answer agrees with the
    // Green representation on the axis and with z off it.
    let lap3 = "diff(diff(u(r, th, ph), r), r) + 2*diff(u(r, th, ph), r)/r + diff(diff(u(r, th, ph), th), th)/r^2 + cos(th)/sin(th)*diff(u(r, th, ph), th)/r^2 + diff(diff(u(r, th, ph), ph), ph)/(r^2*sin(th)^2)";
    let solved = any(&format!(
        "pdsolve({lap3} = 0, u(r, th, ph), list(at(diff(u(r, th, ph), r), r, 1) = cos(th)))"
    ));
    close(
        value(&rhs(&solved), &[("r", 0.8), ("th", 0.7), ("ph", 0.3)]),
        0.8 * 0.7_f64.cos(),
        1e-12,
        "pdsolve Neumann ball",
    );
    // Neumann disk of radius 2: g = y/R on the circle, u = y.
    let s = run("pde_greens_ball_neumann(2, 2)");
    let g = compiled(&s);
    let at = |x: &[f64], q: &[f64]| g(&bind(&[("x", "s_x"), ("y", "s_y")], x, q));
    let u = |r: f64| {
        let integrand = |phi: f64| {
            let source = [radius * phi.cos(), radius * phi.sin()];
            let field = [0.0, r];
            at(&field, &source) * phi.sin()
        };
        simpson(&integrand, 0.0, 2.0 * PI) * radius
    };
    for r in [0.5_f64, 1.0, 1.5] {
        close(u(r), r, 1e-6, "Neumann disk u = y");
    }
    // Exterior Robin ball of radius 1, h = 1/2: g = cos θ, the decaying
    // dipole u = -(2/3) z/r³; the exterior single layer carries a minus sign.
    let s = expand_sums(&run("pde_greens_exterior_robin(1, 1/2, 3)"), 70);
    let g = compiled(&s);
    let at = |x: &[f64], q: &[f64]| g(&bind(&[("x", "s_x"), ("y", "s_y"), ("z", "s_z")], x, q));
    let u = |r: f64| {
        let integrand = |theta: f64| {
            let source = [theta.sin(), 0.0, theta.cos()];
            let field = [0.0, 0.0, r];
            at(&field, &source) * theta.cos() * theta.sin()
        };
        -simpson(&integrand, 0.0, PI) * 2.0 * PI
    };
    for r in [1.5_f64, 2.0, 3.0] {
        close(
            u(r),
            -2.0 / (3.0 * r * r),
            1e-6,
            "exterior Robin ball dipole",
        );
    }
    // Exterior Robin disk: g = sin φ, u = -2 sin φ/r.
    let s = expand_sums(&run("pde_greens_exterior_robin(1, 1/2, 2)"), 70);
    let g = compiled(&s);
    let at = |x: &[f64], q: &[f64]| g(&bind(&[("x", "s_x"), ("y", "s_y")], x, q));
    let u = |r: f64| {
        let integrand = |phi: f64| {
            let source = [phi.cos(), phi.sin()];
            let field = [0.0, r];
            at(&field, &source) * phi.sin()
        };
        -simpson(&integrand, 0.0, 2.0 * PI)
    };
    for r in [1.5_f64, 2.0, 3.0] {
        close(u(r), -2.0 / r, 1e-6, "exterior Robin disk");
    }
}

/// `v_t = v_rr + (2/r) v_r`-type radial problem `u_t = u_rr + (2Z/r - l(l+1)/r²) u` on `(0, 40)`
/// with Dirichlet ends by Crank-Nicolson; the value at `at`.
fn coulomb_heat_reference(
    z: f64,
    l: f64,
    f: &dyn Fn(f64) -> f64,
    t_end: f64,
    at: f64,
) -> f64 {
    let (n, steps) = (4000_usize, 100_usize);
    let h = 40.0 / n as f64;
    let dt = t_end / steps as f64;
    let mut u: Vec<f64> = (0..=n).map(|i| f(h * i as f64)).collect();
    u[0] = 0.0;
    u[n] = 0.0;
    for _ in 0..steps {
        let (mut lower, mut diag, mut upper, mut rhs_) = (
            vec![0.0; n + 1],
            vec![1.0; n + 1],
            vec![0.0; n + 1],
            vec![0.0; n + 1],
        );
        for i in 1..n {
            let r = h * i as f64;
            let c = -2.0 / (h * h) + 2.0 * z / r - l * (l + 1.0) / (r * r);
            let a = 1.0 / (h * h);
            rhs_[i] = u[i] + 0.5 * dt * (a * u[i - 1] + c * u[i] + a * u[i + 1]);
            lower[i] = -0.5 * dt * a;
            diag[i] = 1.0 - 0.5 * dt * c;
            upper[i] = -0.5 * dt * a;
        }
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
    let position = at / h;
    let i = position.floor() as usize;
    let w = position - i as f64;
    u[i] * (1.0 - w) + u[i + 1] * w
}

#[test]
fn coulomb_continuum_completes_the_expansion() {
    use crate::kernels::special_ext::coulomb_wave_value;
    // The scattering solution solves the radial equation and has the amplitude sqrt(2/pi).
    let (z, l, k) = (1.0_f64, 1_usize, 0.8_f64);
    let u = |r: f64| coulomb_wave_value(z, l, k, r);
    let (r, h) = (7.0_f64, 1e-3_f64);
    let second = (u(r + h) - 2.0 * u(r) + u(r - h)) / (h * h);
    close(
        second,
        (2.0 / (r * r) - 2.0 * z / r - k * k) * u(r),
        1e-4,
        "radial equation",
    );
    let r = 100.0_f64;
    let slope = (u(r + h) - u(r - h)) / (2.0 * h);
    let local = (k * k + 2.0 * z / r - 2.0 / (r * r)).sqrt();
    let amplitude = (local * u(r).powi(2) + slope * slope / local) / k;
    close(amplitude, 2.0 / PI, 2e-2, "amplitude");
    // Data outside the span of the bound states: the heat equation with a Coulomb potential.
    // The returned `defint`s are evaluated by nested Gauss–Legendre quadrature
    // (the new answer keeps the symbolic integrals but the adaptively chosen
    // cut-off can exceed the 16 of the old tabulated scheme).
    let s = any(
        "pdsolve(diff(u(r, t), t) = diff(diff(u(r, t), r), r) + 2/r*u(r, t), u(r, t), list(u(0, t) = 0, u(r, 0) = r*exp(-3*r)))",
    );
    assert!(s.contains("coulomb_wave"), "{s}");
    let (r, t) = (1.2_f64, 0.1_f64);
    let reference = coulomb_heat_reference(1.0, 0.0, &|q| q * (-3.0 * q).exp(), t, r);
    let got = super::extension_tests::coulomb_answer_value(&rhs(&s), &[("r", r), ("t", t)]);
    close(got, reference, 3e-2 * reference.abs().max(0.05), &s);
}

#[test]
fn spheroidal_helmholtz_and_laplace_problems() {
    let pro = "(xi^2 - 1)*diff(diff(u(xi, eta), xi), xi) + 2*xi*diff(u(xi, eta), xi) + (1 - eta^2)*diff(diff(u(xi, eta), eta), eta) - 2*eta*diff(u(xi, eta), eta)";
    // Laplace inside a prolate spheroid: P₂(ξ) P₂(η)/P₂(ξ₀).
    let s = any(&format!(
        "pdsolve({pro} = 0, u(xi, eta), list(u(2, eta) = (3*eta^2 - 1)/2))"
    ));
    let p2 = |x: f64| (3.0 * x * x - 1.0) / 2.0;
    let (xi, eta) = (1.4_f64, 0.3_f64);
    close(
        value(&rhs(&s), &[("xi", xi), ("eta", eta)]),
        p2(xi) * p2(eta) / p2(2.0),
        1e-5,
        &s,
    );
    // Helmholtz, γ = 4: the plane wave cos(2 ξ η) (z = ξ η d).
    let s = any(&format!(
        "pdsolve({pro} + 4*(xi^2 - eta^2)*u(xi, eta) = 0, u(xi, eta), list(u(3/2, eta) = cos(3*eta)))"
    ));
    let (xi, eta) = (1.2_f64, 0.4_f64);
    close(
        value(&rhs(&s), &[("xi", xi), ("eta", eta)]),
        (2.0 * xi * eta).cos(),
        1e-4,
        &s,
    );
    // Oblate coordinates: cos(2 ζ η) solves (ζ²+1) u_ζζ + … + 4 (ζ² + η²) u = 0.
    let obl = "(zeta^2 + 1)*diff(diff(u(zeta, eta), zeta), zeta) + 2*zeta*diff(u(zeta, eta), zeta) + (1 - eta^2)*diff(diff(u(zeta, eta), eta), eta) - 2*eta*diff(u(zeta, eta), eta)";
    let s = any(&format!(
        "pdsolve({obl} + 4*(zeta^2 + eta^2)*u(zeta, eta) = 0, u(zeta, eta), list(u(1, eta) = cos(2*eta)))"
    ));
    let (zeta, eta) = (0.6_f64, -0.5_f64);
    close(
        value(&rhs(&s), &[("zeta", zeta), ("eta", eta)]),
        (2.0 * zeta * eta).cos(),
        1e-4,
        &s,
    );
    // Azimuthal dependence: u = sqrt((ξ² - 1)(1 - η²)) cos φ is harmonic.
    let full = format!(
        "{} + (xi^2 - eta^2)/((xi^2 - 1)*(1 - eta^2))*diff(diff(u(xi, eta, phi), phi), phi)",
        pro.replace("u(xi, eta)", "u(xi, eta, phi)")
    );
    let s = any(&format!(
        "pdsolve({full} = 0, u(xi, eta, phi), list(u(2, eta, phi) = sqrt(3)*sqrt(1 - eta^2)*cos(phi)))"
    ));
    let (xi, eta, phi) = (1.5_f64, 0.3_f64, 0.7_f64);
    close(
        value(&rhs(&s), &[("xi", xi), ("eta", eta), ("phi", phi)]),
        ((xi * xi - 1.0) * (1.0 - eta * eta)).sqrt() * phi.cos(),
        1e-4,
        &s,
    );
}

#[test]
fn sources_outside_a_disk_by_duhamel() {
    let bounds = |text: &str| {
        text.replace(", rho, 1, oo)", ", rho, 1, 9)")
            .replace(", k, 0, oo)", ", k, 0, 8)")
    };
    let lap = "diff(diff(u(r, th, t), r), r) + diff(u(r, th, t), r)/r + diff(diff(u(r, th, t), th), th)/r^2";
    // Zero start, steady forcing exp(-(r - 3)²) cos θ.
    let s = any(&format!(
        "pdsolve(diff(u(r, th, t), t) = {lap} + exp(-(r - 3)^2)*cos(th), u(r, th, t), list(u(1, th, t) = 0, u(oo, th, t) = 0, u(r, th, 0) = 0))"
    ));
    let (r, th, t) = (2.0_f64, 0.6_f64, 0.3_f64);
    let text = expand_integrals(&bounds(&rhs(&s)), 90);
    let forcing = |q: f64| (-(q - 3.0) * (q - 3.0)).exp();
    // Duhamel: ∫_0^t U(t - τ)[forcing] dτ with U the homogeneous evolution (Simpson in τ).
    let intervals = 12_usize;
    let h = t / intervals as f64;
    let mut reference = 0.0;
    for j in 0..=intervals {
        let weight = if j == 0 || j == intervals {
            1.0
        } else if j % 2 == 1 {
            4.0
        } else {
            2.0
        };
        let elapsed = t - h * j as f64;
        let value = if elapsed <= 0.0 {
            forcing(r)
        } else {
            exterior_heat_reference(1.0, 1.0, &forcing, elapsed, r)
        };
        reference += weight * value;
    }
    reference *= h / 3.0 * th.cos();
    close(
        value(&text, &[("r", r), ("th", th), ("t", t)]),
        reference,
        4e-3 * reference.abs().max(0.02),
        &s,
    );
}

#[test]
fn time_varying_sources_outside_a_disk_by_duhamel() {
    let bounds = |text: &str| {
        text.replace(", rho, 1, oo)", ", rho, 1, 9)")
            .replace(", k, 0, oo)", ", k, 0, 8)")
    };
    let lap = "diff(diff(u(r, th, t), r), r) + diff(u(r, th, t), r)/r + diff(diff(u(r, th, t), th), th)/r^2";
    // A sum of separable terms g_i(t) s_i(r, θ): the first harmonic with
    // 1 + t² and the second with t. The Weber coefficient of each s_i is
    // formed once and only the time kernel is integrated against g_i.
    let s = any(&format!(
        "pdsolve(diff(u(r, th, t), t) = {lap} + exp(-(r - 3)^2)*(cos(th)*(1 + t^2) + sin(2*th)*t), u(r, th, t), list(u(1, th, t) = 0, u(oo, th, t) = 0, u(r, th, 0) = 0))"
    ));
    let (r, th, t) = (2.0_f64, 0.6_f64, 0.3_f64);
    let text = expand_integrals(&bounds(&rhs(&s)), 90);
    let forcing = |q: f64| (-(q - 3.0) * (q - 3.0)).exp();
    // Duhamel: ∫_0^t g(τ) U(t - τ)[s] dτ per harmonic (Simpson in τ).
    let intervals = 12_usize;
    let h = t / intervals as f64;
    let (mut first, mut second) = (0.0, 0.0);
    for j in 0..=intervals {
        let weight = if j == 0 || j == intervals {
            1.0
        } else if j % 2 == 1 {
            4.0
        } else {
            2.0
        };
        let (tau, elapsed) = (h * j as f64, t - h * j as f64);
        let (v1, v2) = if elapsed <= 0.0 {
            (forcing(r), forcing(r))
        } else {
            (
                exterior_heat_reference(1.0, 1.0, &forcing, elapsed, r),
                exterior_heat_reference(1.0, 4.0, &forcing, elapsed, r),
            )
        };
        first += weight * (1.0 + tau * tau) * v1;
        second += weight * tau * v2;
    }
    let reference = h / 3.0 * (first * th.cos() + second * (2.0 * th).sin());
    close(
        value(&text, &[("r", r), ("th", th), ("t", t)]),
        reference,
        4e-3 * reference.abs().max(0.02),
        &s,
    );
}

#[test]
fn exterior_time_varying_source_stays_cheap() {
    use crate::graph::Budget;
    use crate::graph::Saturate;
    // Regression for L6: the solution used to be a triple integral per mode
    // whose simplification took minutes (86 000 e-graph nodes and 9 800
    // unions for 1 + t; 92 000 and 13 000 for sin t); now it is the Weber
    // coefficient once per mode with a one-dimensional time kernel, and the
    // integrator declines Bessel integrands at once.
    let lap = "diff(diff(u(r, th, t), r), r) + diff(u(r, th, t), r)/r + diff(diff(u(r, th, t), th), th)/r^2";
    for source in [
        "exp(-(r - 3)^2)*cos(th)*(1 + t)",
        "exp(-(r - 3)^2)*cos(th)*sin(t)",
    ] {
        let src = format!(
            "pdsolve(diff(u(r, th, t), t) = {lap} + {source}, u(r, th, t), list(u(1, th, t) = 0, u(oo, th, t) = 0, u(r, th, 0) = 0))"
        );
        let sets = crate::rules::standard();
        let mut g = Graph::new();
        let engine = Engine::install(&mut g, &sets).unwrap();
        let root = g.parse(&src).unwrap();
        let report = engine.run(
            &mut g,
            &[root],
            &Env::symbolic(),
            &Saturate,
            &Budget::default(),
        );
        assert!(report.solved, "{source}: not solved");
        assert!(
            g.len() < 8_000 && report.unions < 2_000,
            "{source}: {} nodes, {} unions",
            g.len(),
            report.unions
        );
    }
}

// ---------------------------------------------------------------------
// Riemann's method for the hodograph equation of non-integer order.

/// A solution of `(x - y) u_xy + N (u_y - u_x) = 0` for every `N > 0`:
/// the Euler–Poisson integral `∫_0^1 f(y + τ (x - y)) τ^(N-1) (1-τ)^(N-1) dτ`.
fn epd_poisson(
    n: f64,
    f: &dyn Fn(f64) -> f64,
    x: f64,
    y: f64,
) -> f64 {
    crate::kernels::quadrature::integrate(
        |t| f(y + t * (x - y)) * t.powf(n - 1.0) * (1.0 - t).powf(n - 1.0),
        0.0,
        1.0,
        1e-13,
    )
    .value
}

#[test]
fn riemann_function_of_the_hodograph_equation() {
    use super::riemann::epd_riemann;
    use super::riemann::epd_riemann_gradient;
    let (r, s) = (2.0, 0.3);
    for n in [1.5_f64, 0.7, 2.5, 3.0] {
        // R = 1 at the field point, R = ((r - sigma)/(r - s))^N on rho = r and ((rho - s)/(r - s))^N on sigma = s
        close(
            epd_riemann(n, r, s, r, s),
            1.0,
            1e-13,
            "R at the field point",
        );
        for sigma in [-0.4, 0.0, 0.2] {
            close(
                epd_riemann(n, r, s, r, sigma),
                ((r - sigma) / (r - s)).powf(n),
                1e-12,
                "R on rho = r",
            );
        }
        for rho in [0.9, 1.4, 1.9] {
            close(
                epd_riemann(n, r, s, rho, s),
                ((rho - s) / (r - s)).powf(n),
                1e-12,
                "R on sigma = s",
            );
        }
        // the adjoint equation (rho - sigma)^2 R_rs - N (rho - sigma)(R_s - R_r) - 2 N R = 0
        let (rho, sigma) = (1.1, 0.5);
        let h = 1e-4;
        let f = |p: f64, q: f64| epd_riemann(n, r, s, p, q);
        let r_rs = (f(rho + h, sigma + h) - f(rho + h, sigma - h) - f(rho - h, sigma + h)
            + f(rho - h, sigma - h))
            / (4.0 * h * h);
        let (r_r, r_s) = epd_riemann_gradient(n, r, s, rho, sigma);
        close(
            r_r,
            (f(rho + h, sigma) - f(rho - h, sigma)) / (2.0 * h),
            1e-7,
            "R_rho",
        );
        close(
            r_s,
            (f(rho, sigma + h) - f(rho, sigma - h)) / (2.0 * h),
            1e-7,
            "R_sigma",
        );
        let d = rho - sigma;
        let residual = d * d * r_rs - n * d * (r_s - r_r) - 2.0 * n * f(rho, sigma);
        assert!(
            residual.abs() < 1e-5,
            "N = {n}: adjoint residual {residual}"
        );
    }
    // integer N: the Legendre polynomial, R = ((rho - sigma)/(r - s))^2 (1 - 2 z) for N = 2
    let (rho, sigma) = (1.2, 0.4);
    let z = (rho - r) * (sigma - s) / ((r - s) * (rho - sigma));
    close(
        epd_riemann(2.0, r, s, rho, sigma),
        ((rho - sigma) / (r - s)).powi(2) * (1.0 - 2.0 * z),
        1e-12,
        "N = 2",
    );
    // outside the hyperbolic region
    assert!(epd_riemann(1.5, 0.2, 0.9, 1.0, 0.5).is_nan());
}

#[test]
fn goursat_problem_of_the_hodograph_equation_for_non_integer_order() {
    use super::riemann::epd_goursat;
    let f = |v: f64| (0.9 * v).cos() + 0.3 * v * v;
    for n in [1.5_f64, 1.25, 2.5] {
        let (r0, s0) = (1.0, -0.5);
        let phi = |rho: f64| epd_poisson(n, &f, rho, s0);
        let psi = |sigma: f64| epd_poisson(n, &f, r0, sigma);
        for (r, s) in [(1.8, 0.4), (1.3, -0.2), (2.5, 0.9), (1.0, 0.7), (2.0, -0.5)] {
            let got = epd_goursat(n, (r0, s0), &phi, &psi, (r, s)).unwrap();
            let want = epd_poisson(n, &f, r, s);
            close(
                got,
                want,
                1e-8 * (1.0 + want.abs()),
                &format!("N = {n} at ({r}, {s})"),
            );
        }
    }
    // N = 2 against Darboux's closed form t = d_r d_s [(F(r) + G(s))/(r - s)] of pde_hodograph_solution
    let text = run("pde_hodograph_solution(sin(r), cos(s), 2*r + s, r + 2*s)");
    let (_, ts) = split_pair(&text);
    let darboux = |r: f64, s: f64| value(&ts, &[("r", r), ("s", s)]);
    let (r0, s0) = (1.0, -0.5);
    let phi = |rho: f64| darboux(rho, s0);
    let psi = |sigma: f64| darboux(r0, sigma);
    for (r, s) in [(1.8, 0.4), (1.3, -0.2)] {
        let got = epd_goursat(2.0, (r0, s0), &phi, &psi, (r, s)).unwrap();
        close(
            got,
            darboux(r, s),
            1e-8 * (1.0 + got.abs()),
            &format!("Darboux at ({r}, {s})"),
        );
    }
    // invalid geometry
    assert!(epd_goursat(1.5, (0.2, 0.0), &phi, &psi, (0.9, 0.5)).is_none());
    assert!(epd_goursat(1.5, (1.0, 0.0), &phi, &psi, (0.9, 0.5)).is_none());
}

#[test]
fn hodograph_riemann_function_in_hypergeometric_form() {
    use super::riemann::epd_riemann;
    let text = run("pde_hodograph_riemann(3/2)");
    assert!(!text.contains("pde_hodograph_riemann"), "{text}");
    for (r, s, rho, sigma) in [
        (2.0, 0.3, 1.1, 0.5),
        (3.0, -1.0, 1.7, 0.2),
        (1.5, 0.1, 1.5, 0.1),
    ] {
        let got = value(&text, &[("r", r), ("s", s), ("rho", rho), ("sigma", sigma)]);
        close(got, epd_riemann(1.5, r, s, rho, sigma), 1e-10, &text);
    }
    // an integer order gives a polynomial in z (terminating series): N = 3
    let text = run("pde_hodograph_riemann(3)");
    let got = value(
        &text,
        &[("r", 2.0), ("s", 0.3), ("rho", 1.1), ("sigma", 0.5)],
    );
    close(got, epd_riemann(3.0, 2.0, 0.3, 1.1, 0.5), 1e-10, &text);
}
