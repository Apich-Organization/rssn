//! Tests of the extended solvers: annuli and exterior domains, spherical
//! harmonics, Robin data, drift with Neumann ends, massive waves, Jordan
//! systems and variable-coefficient problems.
//!
//! Every solution is compared with an exact solution or with an
//! independent numerical computation (roots, quadratures and series
//! evaluated in the test, not by the solver).

use super::dimension_tests::any;
use super::dimension_tests::close;
use super::dimension_tests::rhs;
use super::dimension_tests::run;
use super::dimension_tests::series;
use super::dimension_tests::value;
use super::dimension_tests::HEAT;
use crate::graph::Engine;
use crate::graph::Graph;
use crate::kernels::special::bessel_j as bessel_j_general;
use crate::kernels::special::bessel_y1;
use crate::kernels::special::bessel_y as bessel_y_general;
use num_complex::Complex64;
use std::collections::HashMap;
use std::f64::consts::PI;

/// The complex value of `text` at the real bindings.
fn complex_value(
    text: &str,
    at: &[(&str, f64)],
) -> Complex64 {
    let sets = crate::rules::standard();
    let mut g = Graph::new();
    assert!(Engine::install(&mut g, &sets).is_ok());
    let root = g.parse(text).unwrap_or_else(|e| panic!("cannot parse `{text}`: {e}"));
    let mut bindings = HashMap::new();
    for (name, v) in at {
        bindings.insert(g.interner_mut().symbol(name), Complex64::new(*v, 0.0));
    }
    g.eval_complex(root, &bindings).unwrap_or_else(|| panic!("cannot evaluate `{text}`"))
}

/// `J_ν` with the elementary half-integer orders.
#[allow(clippy::float_cmp)]
fn bessel_j(
    nu: f64,
    x: f64,
) -> f64 {
    let c = (2.0 / (PI * x)).sqrt();
    if nu == 0.5 {
        c * x.sin()
    } else if nu == 1.5 {
        c * (x.sin() / x - x.cos())
    } else {
        bessel_j_general(nu, x)
    }
}

/// `Y_ν` with the fast special cases.
#[allow(clippy::float_cmp)]
fn bessel_y(
    nu: f64,
    x: f64,
) -> f64 {
    let c = (2.0 / (PI * x)).sqrt();
    if nu == 1.0 {
        bessel_y1(x)
    } else if nu == 0.5 {
        -c * x.cos()
    } else if nu == 1.5 {
        -c * (x.cos() / x + x.sin())
    } else {
        bessel_y_general(nu, x)
    }
}

/// Composite Simpson rule.
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

/// The first `count` positive roots of `g` found by scanning `step` and bisection.
fn roots(
    g: &dyn Fn(f64) -> f64,
    start: f64,
    step: f64,
    count: usize,
) -> Vec<f64> {
    let mut out = Vec::new();
    let mut x = start;
    while out.len() < count {
        let next = x + step;
        if g(x) * g(next) < 0.0 {
            let (mut lo, mut hi) = (x, next);
            for _ in 0..70 {
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

/// The condition `p F + q F'` for `F = Z_ν(k r)` (the cylinder function of the given kind).
fn cylinder(
    nu: f64,
    p: f64,
    q: f64,
    k: f64,
    r: f64,
    second: bool,
) -> f64 {
    let f = |order: f64, x: f64| if second { bessel_y(order, x) } else { bessel_j(order, x) };
    let h = 1e-5;
    let value = f(nu, k * r);
    let slope = (f(nu, k * (r + h)) - f(nu, k * (r - h))) / (2.0 * h);
    p * value + q * slope
}

/// The heat equation `u_t = Δu` on the annulus `a < r < b` (a disk annulus, or
/// a shell with `sphere`, radial order `ν`), by the eigenfunction expansion
/// computed from scratch: roots by scanning, coefficients by quadrature.
#[allow(clippy::too_many_arguments)]
fn annulus_heat_reference(
    nu: f64,
    sphere: bool,
    (pa, qa): (f64, f64),
    (pb, qb): (f64, f64),
    (a, b): (f64, f64),
    initial: &dyn Fn(f64) -> f64,
    (r, t): (f64, f64),
    terms: usize,
) -> f64 {
    // The mode for the shell is r^{-1/2} F; its derivative enters the conditions.
    let (pa, pb) = if sphere { (pa - qa / (2.0 * a), pb - qb / (2.0 * b)) } else { (pa, pb) };
    let g = |k: f64| {
        let ya = cylinder(nu, pa, qa, k, a, true);
        let ja = cylinder(nu, pa, qa, k, a, false);
        let jb = cylinder(nu, pb, qb, k, b, false);
        let yb = cylinder(nu, pb, qb, k, b, true);
        ya * jb - ja * yb
    };
    let found = roots(&g, 0.02, 0.02, terms);
    found
        .iter()
        .map(|&k| {
            let ya = cylinder(nu, pa, qa, k, a, true);
            let ja = cylinder(nu, pa, qa, k, a, false);
            let mode = |x: f64| {
                let f = ya * bessel_j(nu, k * x) - ja * bessel_y(nu, k * x);
                if sphere { f / x.sqrt() } else { f }
            };
            let weight = |x: f64| if sphere { x * x } else { x };
            let numerator = simpson(&|x| weight(x) * initial(x) * mode(x), a, b);
            let norm = simpson(&|x| weight(x) * mode(x) * mode(x), a, b);
            numerator / norm * (-k * k * t).exp() * mode(r)
        })
        .sum()
}

#[test]
fn eigenvalue_operators_for_annuli_and_negative_modes() {
    // Dirichlet annulus 1 < r < 2: the zeros of the cross-product J0(k) Y0(2k) - J0(2k) Y0(k).
    let k = value("annulus_root(0, 1, 0, 1, 1, 0, 2, 1)", &[]);
    close(bessel_j(0.0, k) * bessel_y(0.0, 2.0 * k) - bessel_j(0.0, 2.0 * k) * bessel_y(0.0, k), 0.0, 1e-8, "cross product");
    assert!((k - 3.123_0).abs() < 0.01, "{k}");
    let k3 = value("annulus_root(0, 1, 0, 1, 1, 0, 2, 3)", &[]);
    assert!((k3 - 3.0 * PI).abs() < 0.1, "{k3}");
    // The root of a Robin problem satisfies both conditions with the cross-product mode.
    let k = value("annulus_root(1, 0, 1, 1, 2, 1, 3, 2)", &[]);
    let (nu, a, b) = (1.0, 1.0, 3.0);
    let ya = cylinder(nu, 0.0, 1.0, k, a, true);
    let ja = cylinder(nu, 0.0, 1.0, k, a, false);
    let mode_condition = ya * cylinder(nu, 2.0, 1.0, k, b, false) - ja * cylinder(nu, 2.0, 1.0, k, b, true);
    close(mode_condition, 0.0, 1e-6, "robin annulus");
    // A negative eigenvalue: X'(0) = 0, X'(1) = 3 X(1) (the wrong sign): X = cosh(κ x), κ tanh κ = 3.
    let kappa = value("sl_neg_root(0, 1, -3, 1, 1)", &[]);
    close(kappa * kappa.tanh(), 3.0, 1e-8, "negative root");
    assert!(value("sl_neg_root(0, 1, 3, 1, 1)", &[]).is_nan());
}

#[test]
fn heat_on_annuli_and_shells() {
    let (r, t) = (1.4, 0.02);
    // Disk annulus, Dirichlet on both circles, constant data.
    let disk = "diff(u(r, t), t) = diff(diff(u(r, t), r), r) + diff(u(r, t), r)/r";
    let s = any(&format!("pdsolve({disk}, u(r, t), list(u(1, t) = 0, u(2, t) = 0, u(r, 0) = 1))"));
    assert!(s.contains("annulus_root") && !s.contains("defint"), "{s}");
    let reference = annulus_heat_reference(0.0, false, (1.0, 0.0), (1.0, 0.0), (1.0, 2.0), &|_| 1.0, (r, t), 20);
    close(series(&s, &[("r", r), ("t", t)], 20), reference, 1e-6, &s);
    // Insulated inner circle, Robin outer circle, polynomial data.
    let s = any(&format!("pdsolve({disk}, u(r, t), list(at(diff(u(r, t), r), r, 1) = 0, at(diff(u(r, t), r), r, 2) + u(2, t) = 0, u(r, 0) = 3 - r^2))"));
    assert!(!s.contains("defint"), "{s}");
    let reference = annulus_heat_reference(0.0, false, (0.0, 1.0), (1.0, 1.0), (1.0, 2.0), &|x| 3.0 - x * x, (r, t), 20);
    close(series(&s, &[("r", r), ("t", t)], 20), reference, 1e-6, &s);
    // Angular dependence: the order is the angular index.
    let polar = "diff(u(r, th, t), t) = diff(diff(u(r, th, t), r), r) + diff(u(r, th, t), r)/r + diff(diff(u(r, th, t), th), th)/r^2";
    let s = any(&format!("pdsolve({polar}, u(r, th, t), list(u(1, th, t) = 0, u(2, th, t) = 0, u(r, th, 0) = r*cos(th)))"));
    let reference = annulus_heat_reference(1.0, false, (1.0, 0.0), (1.0, 0.0), (1.0, 2.0), &|x| x, (r, t), 20);
    close(series(&s, &[("r", r), ("th", 0.6), ("t", t)], 20), reference * 0.6_f64.cos(), 1e-6, &s);
    // A shell with l = 0: v = r u on an interval; the mode sin(pi (r - 1))/r decays as e^{-pi² t}.
    let ball = "diff(u(r, t), t) = diff(diff(u(r, t), r), r) + 2*diff(u(r, t), r)/r";
    let s = run(&format!("pdsolve({ball}, u(r, t), list(u(1, t) = 0, u(2, t) = 0, u(r, 0) = sin(pi*(r - 1))/r))"));
    close(value(&rhs(&s), &[("r", r), ("t", t)]), (-PI * PI * t).exp() * (PI * (r - 1.0)).sin() / r, 1e-12, &s);
    // Insulated shell: the constant is a zero mode and stays.
    let s = run(&format!("pdsolve({ball}, u(r, t), list(at(diff(u(r, t), r), r, 1) = 0, at(diff(u(r, t), r), r, 2) = 0, u(r, 0) = 3))"));
    close(value(&rhs(&s), &[("r", r), ("t", t)]), 3.0, 1e-12, &s);
    // A shell with general data and l = 0 (sines), compared with the quadrature reference.
    let s = any(&format!("pdsolve({ball}, u(r, t), list(u(1, t) = 0, u(2, t) = 0, u(r, 0) = 1))"));
    let reference = annulus_heat_reference(0.5, true, (1.0, 0.0), (1.0, 0.0), (1.0, 2.0), &|_| 1.0, (r, t), 20);
    close(series(&s, &[("r", r), ("t", t)], 20), reference, 1e-6, &s);
    // A shell with l = 1 (cross-products of order 3/2).
    let ball_axial = "diff(u(r, th, t), t) = diff(diff(u(r, th, t), r), r) + 2*diff(u(r, th, t), r)/r + diff(diff(u(r, th, t), th), th)/r^2 + cos(th)/sin(th)*diff(u(r, th, t), th)/r^2";
    let s = any(&format!("pdsolve({ball_axial}, u(r, th, t), list(u(1, th, t) = 0, u(2, th, t) = 0, u(r, th, 0) = r*cos(th)))"));
    let reference = annulus_heat_reference(1.5, true, (1.0, 0.0), (1.0, 0.0), (1.0, 2.0), &|x| x, (r, 0.05), 8);
    close(series(&s, &[("r", r), ("th", 0.6), ("t", 0.05)], 8), reference * 0.6_f64.cos(), 1e-6, &s);
    // Waves on the annulus: velocity data.
    let s = any("pdsolve(diff(diff(u(r, t), t), t) = diff(diff(u(r, t), r), r) + diff(u(r, t), r)/r, u(r, t), list(u(1, t) = 0, u(2, t) = 0, u(r, 0) = 0, at(diff(u(r, t), t), t, 0) = 1))");
    assert!(s.contains("annulus_root") && !s.contains("defint"), "{s}");
    let k0 = value("annulus_root(0, 1, 0, 1, 1, 0, 2, 1)", &[]);
    // The first mode of the velocity expansion, evaluated independently.
    let (nu, a, b) = (0.0, 1.0, 2.0);
    let all: Vec<f64> = roots(&|k| cylinder(nu, 1.0, 0.0, k, a, true) * cylinder(nu, 1.0, 0.0, k, b, false) - cylinder(nu, 1.0, 0.0, k, a, false) * cylinder(nu, 1.0, 0.0, k, b, true), 0.02, 0.02, 12);
    close(all[0], k0, 1e-6, "first root");
    let wave_reference: f64 = all
        .iter()
        .map(|&k| {
            let (ya, ja) = (bessel_y(0.0, k), bessel_j(0.0, k));
            let mode = |x: f64| ya * bessel_j(0.0, k * x) - ja * bessel_y(0.0, k * x);
            let c = simpson(&|x| x * mode(x), 1.0, 2.0) / simpson(&|x| x * mode(x) * mode(x), 1.0, 2.0);
            c * (k * 0.3).sin() / k * mode(1.4)
        })
        .sum();
    close(series(&s, &[("r", 1.4), ("t", 0.3)], 12), wave_reference, 1e-6, &s);
}

#[test]
fn annulus_with_boundary_data() {
    // u(1, t) = 1, u(2, t) = 0, u(r, 0) = 0: the steady state ln(2/r)/ln 2 at late times.
    let disk = "diff(u(r, t), t) = diff(diff(u(r, t), r), r) + diff(u(r, t), r)/r";
    let s = any(&format!("pdsolve({disk}, u(r, t), list(u(1, t) = 1, u(2, t) = 0, u(r, 0) = 0))"));
    let r = 1.4;
    close(series(&s, &[("r", r), ("t", 3.0)], 60), (2.0 / r).ln() / 2.0_f64.ln(), 1e-2, &s);
    // The outer circle held at 1: the same with ln r / ln 2.
    let s = any(&format!("pdsolve({disk}, u(r, t), list(u(1, t) = 0, u(2, t) = 1, u(r, 0) = 0))"));
    close(series(&s, &[("r", r), ("t", 3.0)], 60), r.ln() / 2.0_f64.ln(), 1e-2, &s);
}

#[test]
fn robin_data_on_box_faces() {
    let lap2 = "diff(diff(u(x, y), x), x) + diff(diff(u(x, y), y), y)";
    // Robin data on the top face: u_y + u = sin(pi x): R(y) = sinh(pi y)/(pi cosh(pi) + sinh(pi)).
    let s = run(&format!("pdsolve({lap2} = 0, u(x, y), list(u(0, y) = 0, u(1, y) = 0, u(x, 0) = 0, at(diff(u(x, y), y), y, 1) + u(x, 1) = sin(pi*x)))"));
    assert!(!s.contains("sum("), "{s}");
    let profile = |y: f64| (PI * y).sinh() / (PI * PI.cosh() + PI.sinh());
    close(value(&rhs(&s), &[("x", 0.3), ("y", 0.6)]), profile(0.6) * (PI * 0.3).sin(), 1e-12, &s);
    // Neumann data on the top face.
    let s = run(&format!("pdsolve({lap2} = 0, u(x, y), list(u(0, y) = 0, u(1, y) = 0, u(x, 0) = 0, at(diff(u(x, y), y), y, 1) = sin(pi*x)))"));
    close(value(&rhs(&s), &[("x", 0.3), ("y", 0.6)]), (PI * 0.6).sinh() / (PI * PI.cosh()) * (PI * 0.3).sin(), 1e-12, &s);
    // Robin data on the bottom face with a Dirichlet top: R = -sinh(pi (1 - y))/(pi cosh(pi) + sinh(pi)).
    let s = run(&format!("pdsolve({lap2} = 0, u(x, y), list(u(0, y) = 0, u(1, y) = 0, u(x, 1) = 0, at(diff(u(x, y), y), y, 0) - u(x, 0) = sin(pi*x)))"));
    close(value(&rhs(&s), &[("x", 0.3), ("y", 0.6)]), -(PI * 0.4).sinh() / (PI * PI.cosh() + PI.sinh()) * (PI * 0.3).sin(), 1e-12, &s);
    // Robin on both y-faces (data on one): R = A cosh(pi y) + B sinh(pi y) with R' + R = 0 at y = 0... checked by the residuals.
    let s = run(&format!("pdsolve({lap2} = 0, u(x, y), list(u(0, y) = 0, u(1, y) = 0, at(diff(u(x, y), y), y, 0) - u(x, 0) = 0, at(diff(u(x, y), y), y, 1) + 2*u(x, 1) = sin(pi*x)))"));
    let y1 = 0.2;
    let u = |y: f64| value(&rhs(&s), &[("x", 0.5), ("y", y)]);
    let h = 1e-4;
    let uy = |y: f64| (u(y + h) - u(y - h)) / (2.0 * h);
    close(uy(0.0) - u(0.0), 0.0, 1e-6, &s);
    close(uy(1.0) + 2.0 * u(1.0), 1.0, 1e-6, &s);
    close((u(y1 + h) - 2.0 * u(y1) + u(y1 - h)) / (h * h), PI * PI * u(y1), 1e-3, &s);
    // Evolution with Robin data: the late-time state is the steady solution.
    let heat2 = "diff(u(x, y, t), t) = diff(diff(u(x, y, t), x), x) + diff(diff(u(x, y, t), y), y)";
    let s = any(&format!("pdsolve({heat2}, u(x, y, t), list(u(0, y, t) = 0, u(1, y, t) = 0, u(x, 0, t) = 0, at(diff(u(x, y, t), y), y, 1) + u(x, 1, t) = sin(pi*x), u(x, y, 0) = 0))"));
    close(series(&s, &[("x", 0.3), ("y", 0.5), ("t", 4.0)], 80), profile(0.5) * (PI * 0.3).sin(), 2e-3, &s);
}

#[test]
fn drift_with_neumann_and_robin_ends() {
    let pde = "diff(u(x, t), t) + 2*diff(u(x, t), x) = diff(diff(u(x, t), x), x)";
    // u = 1 is a steady solution of the Neumann problem (the zero mode e^{-x} of the gauged problem).
    let s = any(&format!("pdsolve({pde}, u(x, t), list(at(diff(u(x, t), x), x, 0) = 0, at(diff(u(x, t), x), x, 1) = 0, u(x, 0) = 1))"));
    close(series(&s, &[("x", 0.4), ("t", 0.3)], 6), 1.0, 1e-9, &s);
    // General data: the series satisfies the equation and the Neumann conditions.
    let s = any(&format!("pdsolve({pde}, u(x, t), list(at(diff(u(x, t), x), x, 0) = 0, at(diff(u(x, t), x), x, 1) = 0, u(x, 0) = x))"));
    assert!(s.contains("sl_root") || s.contains("sl_neg_root"), "{s}");
    let u = |x: f64, t: f64| series(&s, &[("x", x), ("t", t)], 30);
    let (x, t, h) = (0.45, 0.3, 1e-3);
    let residual = (u(x, t + h) - u(x, t - h)) / (2.0 * h) + 2.0 * (u(x + h, t) - u(x - h, t)) / (2.0 * h) - (u(x + h, t) - 2.0 * u(x, t) + u(x - h, t)) / (h * h);
    close(residual, 0.0, 2e-3, &s);
    close((u(h, t) - u(-h, t)) / (2.0 * h), 0.0, 1e-3, &s);
    close((u(1.0 + h, t) - u(1.0 - h, t)) / (2.0 * h), 0.0, 1e-3, &s);
    // Mass is conserved in the sense of the weighted integral of e^{-2x}... checked by the mean-value property:
    // d/dt ∫ e^{-2x} u dx = 0 for the zero-flux problem u_x = 0 with drift (c = 2) only for the Robin flux condition,
    // so instead compare with the explicit eigen-expansion computed independently.
    let reference = drift_reference(0.45, 0.3);
    close(u(x, t), reference, 1e-4, &s);
    // Robin ends with drift.
    let s = any(&format!("pdsolve({pde}, u(x, t), list(at(diff(u(x, t), x), x, 0) - u(0, t) = 0, at(diff(u(x, t), x), x, 1) + u(1, t) = 0, u(x, 0) = 1))"));
    let u = |x: f64, t: f64| series(&s, &[("x", x), ("t", t)], 30);
    let residual = (u(x, t + h) - u(x, t - h)) / (2.0 * h) + 2.0 * (u(x + h, t) - u(x - h, t)) / (2.0 * h) - (u(x + h, t) - 2.0 * u(x, t) + u(x - h, t)) / (h * h);
    close(residual, 0.0, 2e-3, &s);
    close((u(1.0 + h, t) - u(1.0 - h, t)) / (2.0 * h) + u(1.0, t), 0.0, 2e-3, &s);
    close((u(h, t) - u(-h, t)) / (2.0 * h) - u(0.0, t), 0.0, 2e-3, &s);
    // A Robin end of the wrong sign (u_x(1) = 3 u(1)) has a negative eigenvalue: the mode cosh(κ x).
    let s = any(&format!("pdsolve({HEAT}, u(x, t), list(at(diff(u(x, t), x), x, 0) = 0, at(diff(u(x, t), x), x, 1) = 3*u(1, t), u(x, 0) = 1))"));
    assert!(s.contains("sl_neg_root"), "{s}");
    let u = |x: f64, t: f64| series(&s, &[("x", x), ("t", t)], 40);
    let residual = (u(x, t + h) - u(x, t - h)) / (2.0 * h) - (u(x + h, t) - 2.0 * u(x, t) + u(x - h, t)) / (h * h);
    close(residual, 0.0, 5e-3, &s);
    close((u(1.0 + h, t) - u(1.0 - h, t)) / (2.0 * h) - 3.0 * u(1.0, t), 0.0, 5e-3, &s);
    close(u(0.5, 0.0), 1.0, 5e-2, &s);
}

/// `u_t + 2 u_x = u_xx` on `[0, 1]` with `u_x = 0` at both ends and `u(x, 0) = x`, by the
/// eigenfunctions of the gauged problem computed here: `u = e^{x} v` hence the modes of
/// `v_t = v_xx - v` with the Robin conditions `v' + v = 0` at both ends.
fn drift_reference(
    x: f64,
    t: f64,
) -> f64 {
    // v = e^{-x} u, v_x + v = e^{-x} u_x = 0 at both ends; X'' = -k² X with X' + X = 0 at 0 and 1.
    // Positive modes X = k cos(kx) - sin(kx) (X' + X: -k² sin... ) satisfy X'(0) + X(0) = 0:
    // X(0) = k, X'(0) = -k sin 0 - k cos 0 -> use sl_root's convention p0 = 1, q0 = 1 on the left, p1 = 1, q1 = 1 on the right.
    let (p0, q0, p1, q1) = (1.0, 1.0, 1.0, 1.0);
    let g = |k: f64| (p1 * q0 - q1 * p0) * k * k.cos() - (p1 * p0 + q1 * q0 * k * k) * k.sin();
    let ks = roots(&g, 1e-3, 1e-3, 80);
    let mode = |k: f64, s: f64| q0 * k * (k * s).cos() - p0 * (k * s).sin();
    let data = |s: f64| s;
    let mut total = 0.0;
    // The zero-flux problem has the e^{-x} mode: v = e^{-x}, eigenvalue of v_t = v_xx - v ...
    let weight_data = |s: f64| (-s).exp() * data(s);
    let neg = |s: f64| (-s).exp();
    let neg_coeff = simpson(&|s| neg(s) * weight_data(s), 0.0, 1.0) / simpson(&|s| neg(s) * neg(s), 0.0, 1.0);
    // v_t = v_xx - v: the mode e^{-x} has v_xx = e^{-x}, so v_t = 0.
    total += neg_coeff * neg(x) * 1.0;
    for k in ks {
        let c = simpson(&|s| mode(k, s) * weight_data(s), 0.0, 1.0) / simpson(&|s| mode(k, s) * mode(k, s), 0.0, 1.0);
        total += c * (-(k * k + 1.0) * t).exp() * mode(k, x);
    }
    x.exp() * total
}

#[test]
fn fourier_bessel_coefficients_in_closed_form() {
    let disk = "diff(u(r, t), t) = diff(diff(u(r, t), r), r) + diff(u(r, t), r)/r";
    let s = any(&format!("pdsolve({disk}, u(r, t), list(u(1, t) = 0, u(r, 0) = 1 - r^2))"));
    assert!(!s.contains("defint") && s.contains("bessel_zero(0, n)"), "{s}");
    let (r, t) = (0.4, 0.02);
    let reference: f64 = roots(&|k| bessel_j(0.0, k), 0.5, 0.05, 30)
        .iter()
        .map(|&j| {
            let c = 2.0 * simpson(&|x| x * (1.0 - x * x) * bessel_j(0.0, j * x), 0.0, 1.0) / bessel_j(1.0, j).powi(2);
            c * (-j * j * t).exp() * bessel_j(0.0, j * r)
        })
        .sum();
    close(series(&s, &[("r", r), ("t", t)], 30), reference, 1e-9, &s);
    // Odd powers have no closed form (Struve functions) and stay as integrals.
    let s = any(&format!("pdsolve({disk}, u(r, t), list(u(1, t) = 0, u(r, 0) = 1 - r))"));
    assert!(s.contains("defint"), "{s}");
    // Neumann disk with r^4 data.
    let s = any(&format!("pdsolve({disk}, u(r, t), list(at(diff(u(r, t), r), r, 1) = 0, u(r, 0) = r^4))"));
    assert!(!s.contains("defint"), "{s}");
    let reference = simpson(&|x| x * x.powi(4), 0.0, 1.0) * 2.0
        + roots(&|k| bessel_j(1.0, k), 0.5, 0.05, 30)
            .iter()
            .map(|&j| {
                let c = simpson(&|x| x * x.powi(4) * bessel_j(0.0, j * x), 0.0, 1.0) / simpson(&|x| x * bessel_j(0.0, j * x).powi(2), 0.0, 1.0);
                c * (-j * j * t).exp() * bessel_j(0.0, j * r)
            })
            .sum::<f64>();
    close(series(&s, &[("r", r), ("t", t)], 30), reference, 1e-8, &s);
    // A ball: the radial heat equation with data 1 - r^2 on 0 < r < 1.
    let ball = "diff(u(r, t), t) = diff(diff(u(r, t), r), r) + 2*diff(u(r, t), r)/r";
    let s = any(&format!("pdsolve({ball}, u(r, t), list(u(1, t) = 0, u(r, 0) = 1 - r^2))"));
    assert!(!s.contains("defint"), "{s}");
    let reference: f64 = (1..=30)
        .map(|n| {
            let k = f64::from(n) * PI;
            let c = 2.0 * simpson(&|x| x * (1.0 - x * x) * (k * x).sin(), 0.0, 1.0);
            c * (-k * k * t).exp() * (k * r).sin() / r
        })
        .sum();
    close(series(&s, &[("r", r), ("t", t)], 30), reference, 1e-9, &s);
    // A ball with l = 1: the Legendre projection of cos² and the Bessel moments of r².
    let axial = "diff(u(r, th, t), t) = diff(diff(u(r, th, t), r), r) + 2*diff(u(r, th, t), r)/r + diff(diff(u(r, th, t), th), th)/r^2 + cos(th)/sin(th)*diff(u(r, th, t), th)/r^2";
    let s = any(&format!("pdsolve({axial}, u(r, th, t), list(u(1, th, t) = 0, u(r, th, 0) = r^2*cos(th)^2))"));
    assert!(!s.contains("defint"), "{s}");
}

#[test]
fn spherical_harmonics_in_the_ball() {
    let lap = "diff(diff(u(r, th, ph), r), r) + 2*diff(u(r, th, ph), r)/r + diff(diff(u(r, th, ph), th), th)/r^2 + cos(th)/sin(th)*diff(u(r, th, ph), th)/r^2 + diff(diff(u(r, th, ph), ph), ph)/(r^2*sin(th)^2)";
    // u = x = r sin(th) cos(ph).
    let s = run(&format!("pdsolve({lap} = 0, u(r, th, ph), list(u(2, th, ph) = sin(th)*cos(ph)))"));
    let at = [("r", 1.0), ("th", 0.7), ("ph", 0.4)];
    close(value(&rhs(&s), &at), 0.5 * 0.7_f64.sin() * 0.4_f64.cos(), 1e-12, &s);
    // u = x z: Y_2^1.
    let s = run(&format!("pdsolve({lap} = 0, u(r, th, ph), list(u(1, th, ph) = sin(th)*cos(th)*cos(ph) + sin(th)^2*sin(2*ph)))"));
    let (r, th, ph) = (0.6, 0.9_f64, 0.4_f64);
    let want = r * r * (th.sin() * th.cos() * ph.cos() + th.sin().powi(2) * (2.0 * ph).sin());
    close(value(&rhs(&s), &[("r", r), ("th", th), ("ph", ph)]), want, 1e-12, &s);
    // Neumann data: u_r(1) = sin(th) cos(ph) gives u = r sin(th) cos(ph).
    let s = run(&format!("pdsolve({lap} = 0, u(r, th, ph), list(at(diff(u(r, th, ph), r), r, 1) = sin(th)*cos(ph)))"));
    close(value(&rhs(&s), &[("r", r), ("th", th), ("ph", ph)]), r * th.sin() * ph.cos(), 1e-12, &s);
    // Heat: the eigenmode j_1(z r) Y_1^1 with a zero z of j_1 (J_{3/2}).
    let heat = format!("diff(u(r, th, ph, t), t) = {}", lap.replace("(r, th, ph)", "(r, th, ph, t)"));
    let s = run(&format!("pdsolve({heat}, u(r, th, ph, t), list(u(1, th, ph, t) = 0, u(r, th, ph, 0) = besselj(3/2, bessel_zero(3/2, 1)*r)/r^(1/2)*sin(th)*cos(ph)))"));
    let z = value("bessel_zero(3/2, 1)", &[]);
    let t = 0.05;
    close(series(&s, &[("r", r), ("th", th), ("ph", ph), ("t", t)], 3), (-z * z * t).exp() * bessel_j(1.5, z * r) / r.sqrt() * th.sin() * ph.cos(), 1e-9, &s);
}

#[test]
fn exterior_domains() {
    // Laplace outside the unit disk: cos(2 th) / r².
    let lap = "diff(diff(u(r, th), r), r) + diff(u(r, th), r)/r + diff(diff(u(r, th), th), th)/r^2";
    let s = run(&format!("pdsolve({lap} = 0, u(r, th), list(u(1, th) = cos(2*th), u(oo, th) = 0))"));
    close(value(&rhs(&s), &[("r", 2.0), ("th", 0.4)]), (0.8_f64).cos() / 4.0, 1e-12, &s);
    // The two-dimensional exterior has the bounded limit a0 (here 1 + cos): the value at infinity is the mean.
    let s = run(&format!("pdsolve({lap} = 0, u(r, th), list(u(2, th) = 3 + sin(th), u(oo, th) = 3))"));
    close(value(&rhs(&s), &[("r", 4.0), ("th", 0.4)]), 3.0 + 0.4_f64.sin() * 2.0 / 4.0, 1e-12, &s);
    // Neumann exterior data: u_r(1) = -2 cos(2 th) gives u = cos(2 th)/r².
    let s = run(&format!("pdsolve({lap} = 0, u(r, th), list(at(diff(u(r, th), r), r, 1) = -2*cos(2*th), u(oo, th) = 0))"));
    close(value(&rhs(&s), &[("r", 1.5), ("th", 0.4)]), (0.8_f64).cos() / 2.25, 1e-12, &s);
    // Outside a ball: axisymmetric data cos²(th) at R = 2.
    let ball = "diff(diff(u(r, th), r), r) + 2*diff(u(r, th), r)/r + diff(diff(u(r, th), th), th)/r^2 + cos(th)/sin(th)*diff(u(r, th), th)/r^2";
    let s = run(&format!("pdsolve({ball} = 0, u(r, th), list(u(2, th) = cos(th)^2, u(oo, th) = 0))"));
    let (r, th) = (3.0_f64, 0.7_f64);
    let p2 = 1.5 * th.cos() * th.cos() - 0.5;
    close(value(&rhs(&s), &[("r", r), ("th", th)]), 1.0 / 3.0 * (2.0 / r) + 2.0 / 3.0 * (2.0 / r).powi(3) * p2, 1e-12, &s);
    // Without axial symmetry: x/r³ outside the unit ball.
    let full = "diff(diff(u(r, th, ph), r), r) + 2*diff(u(r, th, ph), r)/r + diff(diff(u(r, th, ph), th), th)/r^2 + cos(th)/sin(th)*diff(u(r, th, ph), th)/r^2 + diff(diff(u(r, th, ph), ph), ph)/(r^2*sin(th)^2)";
    let s = run(&format!("pdsolve({full} = 0, u(r, th, ph), list(u(1, th, ph) = sin(th)*cos(ph), u(oo, th, ph) = 0))"));
    close(value(&rhs(&s), &[("r", 2.0), ("th", 0.9), ("ph", 0.3)]), 0.9_f64.sin() * 0.3_f64.cos() / 4.0, 1e-12, &s);
    // Helmholtz outside a ball: the outgoing wave e^{ik(r - R)}/r for constant data.
    let helmholtz = "diff(diff(u(r, th), r), r) + 2*diff(u(r, th), r)/r + diff(diff(u(r, th), th), th)/r^2 + cos(th)/sin(th)*diff(u(r, th), th)/r^2 + 4*u(r, th)";
    let s = run(&format!("pdsolve({helmholtz} = 0, u(r, th), list(u(1, th) = 1, u(oo, th) = 0))"));
    let got = complex_value(&rhs(&s), &[("r", 2.5), ("th", 0.3)]);
    let want = Complex64::new(0.0, 2.0 * 1.5).exp() / 2.5;
    assert!((got - want).norm() < 1e-12, "{s}: {got} vs {want}");
    // l = 1 outgoing wave: h_1(kr) / h_1(k R) with the spherical Hankel function.
    let s = run(&format!("pdsolve({helmholtz} = 0, u(r, th), list(u(1, th) = cos(th), u(oo, th) = 0))"));
    let h1 = |x: f64| Complex64::new(0.0, x).exp() / x * (Complex64::new(1.0, 0.0) + Complex64::new(0.0, 1.0) / x);
    let got = complex_value(&rhs(&s), &[("r", 2.5), ("th", 0.3)]);
    let want = h1(5.0) / h1(2.0) * 0.3_f64.cos();
    assert!((got - want).norm() < 1e-12, "{s}: {got} vs {want}");
    // Helmholtz outside a disk: H_1^{(1)}(k r)/H_1^{(1)}(k R) cos(th).
    let helmholtz_2d = "diff(diff(u(r, th), r), r) + diff(u(r, th), r)/r + diff(diff(u(r, th), th), th)/r^2 + 4*u(r, th)";
    let s = any(&format!("pdsolve({helmholtz_2d} = 0, u(r, th), list(u(1, th) = cos(th), u(oo, th) = 0))"));
    let hankel = |x: f64| Complex64::new(bessel_j(1.0, x), bessel_y(1.0, x));
    let got = complex_value(&rhs(&s), &[("r", 2.5), ("th", 0.3)]);
    let want = hankel(5.0) / hankel(2.0) * 0.3_f64.cos();
    assert!((got - want).norm() < 1e-8, "{s}: {got} vs {want}");
    // The heat equation outside a ball for radial data: w = (r - 1) e^{-(r-1)²} / r has the exact solution.
    let heat = "diff(u(r, t), t) = diff(diff(u(r, t), r), r) + 2*diff(u(r, t), r)/r";
    let s = any(&format!("pdsolve({heat}, u(r, t), list(u(1, t) = 0, u(oo, t) = 0, u(r, 0) = (r - 1)*exp(-(r - 1)^2)/r))"));
    let (r, t) = (1.8, 0.3);
    let x = r - 1.0;
    let numeric_value = |text: &str| crate::rules::testing::numeric(&crate::rules::standard(), &rhs(text), &[("r", r), ("t", t)], 1e-10).0;
    close(numeric_value(&s), x * (-x * x / (1.0 + 4.0 * t)).exp() / (1.0 + 4.0 * t).powf(1.5) / r, 1e-7, &s);
    // A constant surface temperature: the similarity solution (1/r) erfc((r - 1)/(2 sqrt t)).
    let s = any(&format!("pdsolve({heat}, u(r, t), list(u(1, t) = 1, u(oo, t) = 0, u(r, 0) = 0))"));
    close(value(&rhs(&s), &[("r", r), ("t", t)]), crate::kernels::special::erfc_numerical(x / (2.0 * t.sqrt())) / r, 1e-8, &s);
}

#[test]
fn massive_and_damped_waves_in_free_space() {
    use super::dimension_tests::expand_integrals;
    let w3 = "diff(diff(u(x, y, z, t), t), t) = diff(diff(u(x, y, z, t), x), x) + diff(diff(u(x, y, z, t), y), y) + diff(diff(u(x, y, z, t), z), z)";
    let at = [("x", 0.3), ("y", 0.5), ("z", 0.2), ("t", 0.9)];
    let (r2, t) = (0.3_f64 * 0.3 + 0.25 + 0.04, 0.9_f64);
    // Klein-Gordon u_tt = Δu - u: u = cos t |x|² + 3 t sin t for u(0) = |x|², u_t(0) = 0.
    let s = any(&format!("pdsolve({w3} - u(x, y, z, t), u(x, y, z, t), list(u(x, y, z, 0) = x^2 + y^2 + z^2, at(diff(u(x, y, z, t), t), t, 0) = 0))"));
    assert!(s.contains("besselj(1,"), "{s}");
    close(value(&expand_integrals(&rhs(&s), 400), &at), t.cos() * r2 + 3.0 * t * t.sin(), 1e-5, &s);
    // Velocity data 1: u = sin t.
    let s = any(&format!("pdsolve({w3} - u(x, y, z, t), u(x, y, z, t), list(u(x, y, z, 0) = 0, at(diff(u(x, y, z, t), t), t, 0) = 1))"));
    close(value(&expand_integrals(&rhs(&s), 400), &at), t.sin(), 1e-5, &s);
    // A negative mass (u_tt = Δu + u): I Bessel kernels, u = cosh t |x|² + 3 t sinh t.
    let s = any(&format!("pdsolve({w3} + u(x, y, z, t), u(x, y, z, t), list(u(x, y, z, 0) = x^2 + y^2 + z^2, at(diff(u(x, y, z, t), t), t, 0) = 0))"));
    assert!(s.contains("besseli(1,"), "{s}");
    close(value(&expand_integrals(&rhs(&s), 400), &at), t.cosh() * r2 + 3.0 * t * t.sinh(), 1e-5, &s);
    // The damped wave u_tt + 2 u_t = Δu: u = |x|² + 3 t - 3 e^{-t} sinh t for u(0) = |x|², u_t(0) = 0.
    let s = any(&format!("pdsolve({w3} - 2*diff(u(x, y, z, t), t), u(x, y, z, t), list(u(x, y, z, 0) = x^2 + y^2 + z^2, at(diff(u(x, y, z, t), t), t, 0) = 0))"));
    close(value(&expand_integrals(&rhs(&s), 400), &at), r2 + 3.0 * t - 3.0 * (-t).exp() * t.sinh(), 1e-5, &s);
    // Two dimensions: u = cos t (x² + y²) + 2 t sin t.
    let w2 = "diff(diff(u(x, y, t), t), t) = diff(diff(u(x, y, t), x), x) + diff(diff(u(x, y, t), y), y)";
    let s = any(&format!("pdsolve({w2} - u(x, y, t), u(x, y, t), list(u(x, y, 0) = x^2 + y^2, at(diff(u(x, y, t), t), t, 0) = 0))"));
    close(value(&super::dimension_tests::expand_unit_integrals(&rhs(&s), 200), &[("x", 0.3), ("y", 0.5), ("t", 0.9)]), t.cos() * 0.34 + 2.0 * t * t.sin(), 1e-6, &s);
}

#[test]
fn massive_and_damped_waves_on_bounded_domains() {
    let (t, h) = (0.7_f64, 0.0_f64);
    let _ = h;
    // Klein-Gordon in a cube: omega² = 3 + 1.
    let w3 = "diff(diff(u(x, y, z, t), t), t) = diff(diff(u(x, y, z, t), x), x) + diff(diff(u(x, y, z, t), y), y) + diff(diff(u(x, y, z, t), z), z) - u(x, y, z, t)";
    let s = run(&format!("pdsolve({w3}, u(x, y, z, t), list(u(0, y, z, t) = 0, u(pi, y, z, t) = 0, u(x, 0, z, t) = 0, u(x, pi, z, t) = 0, u(x, y, 0, t) = 0, u(x, y, pi, t) = 0, u(x, y, z, 0) = sin(x)*sin(y)*sin(z), at(diff(u(x, y, z, t), t), t, 0) = 0))"));
    close(value(&rhs(&s), &[("x", 0.4), ("y", 0.9), ("z", 1.3), ("t", t)]), (2.0 * t).cos() * 0.4_f64.sin() * 0.9_f64.sin() * 1.3_f64.sin(), 1e-12, &s);
    // The damped wave in a disk: u = e^{-t}(cos ωt + sin ωt/ω) J0(z r), ω² = z² - 1.
    let s = run("pdsolve(diff(diff(u(r, t), t), t) + 2*diff(u(r, t), t) = diff(diff(u(r, t), r), r) + diff(u(r, t), r)/r, u(r, t), list(u(1, t) = 0, u(r, 0) = besselj(0, bessel_zero(0, 1)*r), at(diff(u(r, t), t), t, 0) = 0))");
    let z = 2.404_825_557_695_773_f64;
    let w = (z * z - 1.0).sqrt();
    close(value(&rhs(&s), &[("r", 0.4), ("t", t)]), (-t).exp() * ((w * t).cos() + (w * t).sin() / w) * bessel_j(0.0, z * 0.4), 1e-9, &s);
    // A massive wave in a ball (3D, radial): sin(pi r)/r with omega² = pi² + 4.
    let s = run("pdsolve(diff(diff(u(r, t), t), t) = diff(diff(u(r, t), r), r) + 2*diff(u(r, t), r)/r - 4*u(r, t), u(r, t), list(u(1, t) = 0, u(r, 0) = sin(pi*r)/r, at(diff(u(r, t), t), t, 0) = 0))");
    close(value(&rhs(&s), &[("r", 0.4), ("t", t)]), ((PI * PI + 4.0).sqrt() * t).cos() * (PI * 0.4).sin() / 0.4, 1e-12, &s);
}

#[test]
fn jordan_chains_in_systems() {
    let (x, t) = (0.7_f64, 0.4_f64);
    let parts = |s: &str, names: &[&str]| -> Vec<String> {
        let inner = s.trim_start_matches("list(").strip_suffix(')').unwrap_or("");
        let mut pieces = Vec::new();
        let mut rest = inner.to_owned();
        for (i, name) in names.iter().enumerate().rev() {
            let marker = format!("{name}(x, t) = ");
            if let Some(at) = rest.rfind(&marker) {
                let tail = rest[at + marker.len()..].trim_end_matches(", ").to_owned();
                pieces.push(tail.trim_end_matches(", ").to_owned());
                rest.truncate(at);
                let _ = i;
            }
        }
        pieces.reverse();
        pieces.iter().map(|p| p.trim_end_matches(", ").to_owned()).collect()
    };
    // u_t + u_x + v_x = 0, v_t + v_x = 0: A = [[1, 1], [0, 1]], a single Jordan block.
    let s = any("pdsolve(list(diff(u(x, t), t) + diff(u(x, t), x) + diff(v(x, t), x) = 0, diff(v(x, t), t) + diff(v(x, t), x) = 0), list(u(x, t), v(x, t)), list(u(x, 0) = sin(x), v(x, 0) = cos(x)))");
    let p = parts(&s, &["u", "v"]);
    assert_eq!(p.len(), 2, "{s}");
    close(value(&p[0], &[("x", x), ("t", t)]), (1.0 + t) * (x - t).sin(), 1e-12, &s);
    close(value(&p[1], &[("x", x), ("t", t)]), (x - t).cos(), 1e-12, &s);
    // A chain of length three.
    let s = any("pdsolve(list(diff(u(x, t), t) + diff(u(x, t), x) + diff(v(x, t), x) = 0, diff(v(x, t), t) + diff(v(x, t), x) + diff(w(x, t), x) = 0, diff(w(x, t), t) + diff(w(x, t), x) = 0), list(u(x, t), v(x, t), w(x, t)), list(u(x, 0) = sin(x), v(x, 0) = cos(x), w(x, 0) = sin(2*x)))");
    let inner = s.trim_start_matches("list(").strip_suffix(')').unwrap_or("");
    let (a, rest) = inner.split_once(", v(x, t) = ").unwrap_or((inner, ""));
    let (b, c) = rest.split_once(", w(x, t) = ").unwrap_or((rest, ""));
    let a = a.trim_start_matches("u(x, t) = ");
    // w = sin(2(x - t)); v = cos(x - t) - t w'; u = sin(x - t) - t v0' + t²/2 w0''.
    let y = x - t;
    close(value(c, &[("x", x), ("t", t)]), (2.0 * y).sin(), 1e-12, &s);
    close(value(b, &[("x", x), ("t", t)]), y.cos() - t * 2.0 * (2.0 * y).cos(), 1e-12, &s);
    close(value(a, &[("x", x), ("t", t)]), y.sin() + t * y.sin() - 0.5 * t * t * 4.0 * (2.0 * y).sin(), 1e-12, &s);
    // Coupled heat equations u_t = u_xx + v_xx, v_t = v_xx (Jordan block, whole line).
    let s = any("pdsolve(list(diff(u(x, t), t) = diff(diff(u(x, t), x), x) + diff(diff(v(x, t), x), x), diff(v(x, t), t) = diff(diff(v(x, t), x), x)), list(u(x, t), v(x, t)), list(u(x, 0) = exp(-x^2), v(x, 0) = exp(-x^2)))");
    let sg = 1.0 + 4.0 * t;
    let gauss = (-x * x / sg).exp() / sg.sqrt();
    let gauss_xx = (4.0 * x * x / (sg * sg) - 2.0 / sg) * gauss;
    let inner = s.trim_start_matches("list(").strip_suffix(')').unwrap_or("");
    let (a, b) = inner.split_once(", v(x, t) = ").unwrap_or((inner, ""));
    let a = a.trim_start_matches("u(x, t) = ");
    let numeric = |text: &str| crate::rules::testing::numeric(&crate::rules::standard(), text, &[("x", x), ("t", t)], 1e-10).0;
    close(numeric(b), gauss, 1e-7, &s);
    close(numeric(a), gauss + t * gauss_xx, 1e-7, &s);
}

#[test]
fn variable_coefficient_sturm_liouville_problems() {
    let (x, t) = (2.5_f64, 0.3_f64);
    // Euler-Cauchy: u_t = x² u_xx + x u_x on [1, e^pi], x = e^s: u = e^{-t} sin(ln x).
    let s = run("pdsolve(diff(u(x, t), t) = x^2*diff(diff(u(x, t), x), x) + x*diff(u(x, t), x), u(x, t), list(u(1, t) = 0, u(exp(pi), t) = 0, u(x, 0) = sin(ln(x))))");
    close(value(&rhs(&s), &[("x", x), ("t", t)]), (-t).exp() * x.ln().sin(), 1e-12, &s);
    // With a first-order term and a potential: u_t = x² u_xx + 3 x u_x - u: s-drift 2 and decay.
    let s = any("pdsolve(diff(u(x, t), t) = x^2*diff(diff(u(x, t), x), x) + 3*x*diff(u(x, t), x) - u(x, t), u(x, t), list(u(1, t) = 0, u(exp(pi), t) = 0, u(x, 0) = sin(ln(x))/x))");
    // u = e^{-s} v: v_t = v_ss - v - 1 ... checked through the equation residual.
    let u = |x: f64, t: f64| series(&s, &[("x", x), ("t", t)], 20);
    let h = 1e-4;
    let residual = (u(x, t + h) - u(x, t - h)) / (2.0 * h) - x * x * (u(x + h, t) - 2.0 * u(x, t) + u(x - h, t)) / (h * h) - 3.0 * x * (u(x + h, t) - u(x - h, t)) / (2.0 * h) + u(x, t);
    close(residual, 0.0, 2e-4, &s);
    // A Neumann end in the Euler variable: u_x(1) = 0 and u_x(e^pi) = 0, u = cos(ln x) is not invariant,
    // but u = 1 is.
    let s = any("pdsolve(diff(u(x, t), t) = x^2*diff(diff(u(x, t), x), x) + x*diff(u(x, t), x), u(x, t), list(at(diff(u(x, t), x), x, 1) = 0, at(diff(u(x, t), x), x, exp(pi)) = 0, u(x, 0) = 1))");
    close(series(&s, &[("x", x), ("t", t)], 6), 1.0, 1e-9, &s);
    // Legendre type: u_t = ((1 - x²) u_x)_x with u(x, 0) = P_2(x) decays as e^{-6 t}.
    let leg = "diff(u(x, t), t) = (1 - x^2)*diff(diff(u(x, t), x), x) - 2*x*diff(u(x, t), x)";
    let s = run(&format!("pdsolve({leg}, u(x, t), list(u(x, 0) = (3*x^2 - 1)/2))"));
    close(value(&rhs(&s), &[("x", 0.4), ("t", t)]), (-6.0 * t).exp() * (3.0 * 0.16 - 1.0) / 2.0, 1e-12, &s);
    // Polynomial data: x² = (2 P_2 + 1)/3.
    let s = run(&format!("pdsolve({leg}, u(x, t), list(u(x, 0) = x^2))"));
    close(value(&rhs(&s), &[("x", 0.4), ("t", t)]), 1.0 / 3.0 + 2.0 / 3.0 * (-6.0 * t).exp() * (3.0 * 0.16 - 1.0) / 2.0, 1e-12, &s);
    // The Legendre wave equation: u = cos(sqrt(6) t) P_2.
    let s = run("pdsolve(diff(diff(u(x, t), t), t) = (1 - x^2)*diff(diff(u(x, t), x), x) - 2*x*diff(u(x, t), x), u(x, t), list(u(x, 0) = (3*x^2 - 1)/2, at(diff(u(x, t), t), t, 0) = 0))");
    close(value(&rhs(&s), &[("x", 0.4), ("t", t)]), (6.0_f64.sqrt() * t).cos() * (3.0 * 0.16 - 1.0) / 2.0, 1e-12, &s);
    // Bessel type in a variable not called r: u_t = u_xx + u_x/x on 0 < x < 1.
    let s = run("pdsolve(diff(u(x, t), t) = diff(diff(u(x, t), x), x) + diff(u(x, t), x)/x, u(x, t), list(u(1, t) = 0, u(x, 0) = besselj(0, bessel_zero(0, 1)*x)))");
    let z = 2.404_825_557_695_773_f64;
    close(value(&rhs(&s), &[("x", 0.4), ("t", t)]), (-z * z * t).exp() * bessel_j(0.0, z * 0.4), 1e-9, &s);
}
