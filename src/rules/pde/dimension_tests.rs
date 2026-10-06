//! Tests of the solvers that treat one, two and three dimensions uniformly:
//! eigenfunction expansions on boxes, disks, cylinders and balls, heat and
//! wave kernels on whole space and half-spaces, half-space Green's
//! functions and linear systems.
//!
//! Closed forms are compared with the exact solutions at sample points;
//! series are cut off after a few dozen terms and compared with exact or
//! independently computed values.

use crate::rules::testing::eval;
use crate::rules::testing::reduce_with;
use crate::rules::testing::simplify;
use std::f64::consts::PI;

const HEAT: &str = "diff(u(x, t), t) = diff(diff(u(x, t), x), x)";
const HEAT_2D: &str = "diff(u(x, y, t), t) = diff(diff(u(x, y, t), x), x) + diff(diff(u(x, y, t), y), y)";
const LAPLACE_2D: &str = "diff(diff(u(x, y), x), x) + diff(diff(u(x, y), y), y)";
const WAVE: &str = "diff(diff(u(x, t), t), t) = diff(diff(u(x, t), x), x)";

fn run(src: &str) -> String {
    simplify(&crate::rules::standard(), src)
}

fn any(src: &str) -> String {
    reduce_with(&crate::rules::standard(), src, &[]).0
}

/// The right-hand side of `lhs = rhs`.
fn rhs(text: &str) -> String {
    text.split_once(" = ").map_or_else(|| text.to_owned(), |(_, r)| r.to_owned())
}

fn value(
    text: &str,
    at: &[(&str, f64)],
) -> f64 {
    eval(&crate::rules::standard(), text, at)
}

fn close(
    got: f64,
    want: f64,
    tolerance: f64,
    context: &str,
) {
    assert!((got - want).abs() <= tolerance, "{context}: {got} vs {want}");
}

/// Replaces every `sum(body, n, lo, oo)` by the explicit sum of its first
/// `terms` terms.
fn expand_sums(
    text: &str,
    terms: i64,
) -> String {
    let mut s = text.to_owned();
    while let Some(start) = s.rfind("sum(") {
        let open = start + 3;
        let mut depth = 0_i32;
        let mut end = open;
        for (i, c) in s.char_indices().skip(open) {
            match c {
                | '(' => depth += 1,
                | ')' => {
                    depth -= 1;
                    if depth == 0 {
                        end = i;
                        break;
                    }
                },
                | _ => {},
            }
        }
        let inner = s[open + 1..end].to_owned();
        let mut parts = Vec::new();
        let (mut depth, mut last) = (0_i32, 0);
        for (i, c) in inner.char_indices() {
            match c {
                | '(' => depth += 1,
                | ')' => depth -= 1,
                | ',' if depth == 0 => {
                    parts.push(inner[last..i].trim().to_owned());
                    last = i + 1;
                },
                | _ => {},
            }
        }
        parts.push(inner[last..].trim().to_owned());
        assert_eq!(parts.len(), 4, "malformed sum in {text}");
        let (body, var, lo, hi) = (&parts[0], &parts[1], parts[2].parse::<i64>().unwrap_or(1), &parts[3]);
        let hi = if hi == "oo" { terms } else { hi.parse::<i64>().unwrap_or(terms) };
        let pieces: Vec<String> = (lo..=hi).map(|k| replace_identifier(body, var, &format!("({k})"))).collect();
        s.replace_range(start..=end, &format!("({})", pieces.join(" + ")));
    }
    s
}

/// `text` with the identifier `name` replaced.
fn replace_identifier(
    text: &str,
    name: &str,
    with: &str,
) -> String {
    let chars: Vec<char> = text.chars().collect();
    let mut out = String::new();
    let mut i = 0;
    while i < chars.len() {
        if chars[i].is_alphabetic() || chars[i] == '_' {
            let mut j = i;
            while j < chars.len() && (chars[j].is_alphanumeric() || chars[j] == '_') {
                j += 1;
            }
            let word: String = chars[i..j].iter().collect();
            out.push_str(&if word == name { with.to_owned() } else { word });
            i = j;
        } else {
            out.push(chars[i]);
            i += 1;
        }
    }
    out
}

/// A series solution evaluated with `terms` terms.
fn series(
    text: &str,
    at: &[(&str, f64)],
    terms: i64,
) -> f64 {
    value(&expand_sums(&rhs(text), terms), at)
}

#[test]
fn numeric_eigenvalue_operators() {
    let v = |src: &str| value(src, &[]);
    close(v("bessel_zero(0, 1)"), 2.404_825_557_695_773, 1e-9, "j01");
    close(v("bessel_zero(0, 3)"), 8.653_727_912_911_012, 1e-8, "j03");
    close(v("bessel_zero(1, 2)"), 7.015_586_669_815_619, 1e-8, "j12");
    close(v("bessel_zero(1/2, 2)"), 2.0 * PI, 1e-8, "half order");
    close(v("bessel_root(0, 0, 1)"), 3.831_705_970_207_512, 1e-8, "j'01");
    close(v("bessel_root(0, 0, 0)"), 0.0, 0.0, "trivial root");
    close(v("sl_root(1, 0, 1, 0, 1, 2)"), 2.0 * PI, 1e-9, "dirichlet");
    close(v("sl_root(0, 1, 0, 1, pi, 3)"), 3.0, 1e-9, "neumann");
    // Robin: -k sin k + 2 cos k = 0.
    let k = v("sl_root(0, 1, 2, 1, 1, 1)");
    close(-k * k.sin() + 2.0 * k.cos(), 0.0, 1e-9, "robin residual");
    assert!(k > 0.0 && k < PI / 2.0, "{k}");
}

#[test]
fn intervals_with_every_boundary_type() {
    let (x, t) = (0.3, 0.05);
    let at = [("x", x), ("t", t)];
    // Dirichlet at both ends.
    let s = run(&format!("pdsolve({HEAT}, u(x, t), list(u(0, t) = 0, u(1, t) = 0, u(x, 0) = sin(pi*x) + 2*sin(3*pi*x)))"));
    let exact = (-PI * PI * t).exp() * (PI * x).sin() + 2.0 * (-9.0 * PI * PI * t).exp() * (3.0 * PI * x).sin();
    close(value(&rhs(&s), &at), exact, 1e-12, &s);
    // Neumann at 0, Dirichlet at 1: modes cos((n - 1/2) pi x).
    let s = run(&format!("pdsolve({HEAT}, u(x, t), list(at(diff(u(x, t), x), x, 0) = 0, u(1, t) = 0, u(x, 0) = cos(pi*x/2)))"));
    close(value(&rhs(&s), &at), (-PI * PI * t / 4.0).exp() * (PI * x / 2.0).cos(), 1e-12, &s);
    // Dirichlet at 0, Neumann at 1.
    let s = run(&format!("pdsolve({HEAT}, u(x, t), list(u(0, t) = 0, at(diff(u(x, t), x), x, 1) = 0, u(x, 0) = sin(3*pi*x/2)))"));
    close(value(&rhs(&s), &at), (-9.0 * PI * PI * t / 4.0).exp() * (1.5 * PI * x).sin(), 1e-12, &s);
    // A shifted interval [1, 3].
    let s = run(&format!("pdsolve({HEAT}, u(x, t), list(u(1, t) = 0, u(3, t) = 0, u(x, 0) = sin(pi*(x - 1)/2)))"));
    close(value(&rhs(&s), &[("x", 1.7), ("t", 0.2)]), (-PI * PI * 0.2 / 4.0).exp() * (PI * 0.7 / 2.0).sin(), 1e-12, &s);
    // Periodic: u(0, t) = u(2 pi, t).
    let s = run(&format!("pdsolve({HEAT}, u(x, t), list(u(0, t) = u(2*pi, t), u(x, 0) = 1 + cos(x) + sin(2*x)))"));
    close(value(&rhs(&s), &at), 1.0 + (-t).exp() * x.cos() + (-4.0 * t).exp() * (2.0 * x).sin(), 1e-12, &s);
}

/// Eigenvalue roots of `-k sin k + 2 cos k = 0` and the Robin series of
/// `u_t = u_xx`, `u_x(0) = 0`, `u_x(1) + 2 u(1) = 0`, `u(x, 0) = 1`.
fn robin_reference(
    x: f64,
    t: f64,
) -> f64 {
    let g = |k: f64| -k * k.sin() + 2.0 * k.cos();
    let mut total = 0.0;
    let mut k = 1e-6;
    let mut found = 0;
    while found < 60 {
        let next = k + 1e-3;
        if g(k) * g(next) < 0.0 {
            let (mut lo, mut hi) = (k, next);
            for _ in 0..60 {
                let mid = 0.5 * (lo + hi);
                if g(lo) * g(mid) <= 0.0 {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            let root = 0.5 * (lo + hi);
            let numerator = root.sin() / root;
            let norm = 0.5 + (2.0 * root).sin() / (4.0 * root);
            total += numerator / norm * (-root * root * t).exp() * (root * x).cos();
            found += 1;
        }
        k = next;
    }
    total
}

#[test]
fn robin_conditions() {
    let s = run(&format!(
        "pdsolve({HEAT}, u(x, t), list(at(diff(u(x, t), x), x, 0) = 0, at(diff(u(x, t), x), x, 1) + 2*u(1, t) = 0, u(x, 0) = 1))"
    ));
    assert!(s.contains("sl_root"), "{s}");
    for (x, t) in [(0.2, 0.05), (0.7, 0.2), (1.0, 0.1)] {
        close(series(&s, &[("x", x), ("t", t)], 12), robin_reference(x, t), 1e-6, &s);
    }
    // The Robin condition may also be written with the value on the right.
    let r = run(&format!(
        "pdsolve({HEAT}, u(x, t), list(at(diff(u(x, t), x), x, 0) = 0, at(diff(u(x, t), x), x, 1) = -2*u(1, t), u(x, 0) = 1))"
    ));
    close(series(&r, &[("x", 0.4), ("t", 0.1)], 12), robin_reference(0.4, 0.1), 1e-6, &r);
}

#[test]
fn wave_like_equations_on_intervals() {
    let (x, t) = (0.8, 0.9);
    let at = [("x", x), ("t", t)];
    // Dirichlet wave with speed 2.
    let s = run("pdsolve(diff(diff(u(x, t), t), t) = 4*diff(diff(u(x, t), x), x), u(x, t), list(u(0, t) = 0, u(pi, t) = 0, u(x, 0) = sin(3*x), at(diff(u(x, t), t), t, 0) = sin(x)))");
    close(value(&rhs(&s), &at), (3.0 * x).sin() * (6.0 * t).cos() + x.sin() * (2.0 * t).sin() / 2.0, 1e-12, &s);
    // Neumann ends: the constant mode moves with the initial velocity.
    let s = run(&format!("pdsolve({WAVE}, u(x, t), list(at(diff(u(x, t), x), x, 0) = 0, at(diff(u(x, t), x), x, pi) = 0, u(x, 0) = cos(x), at(diff(u(x, t), t), t, 0) = 1))"));
    close(value(&rhs(&s), &at), x.cos() * t.cos() + t, 1e-12, &s);
    // Telegraph equation u_tt + 2 u_t = u_xx.
    let s = run("pdsolve(diff(diff(u(x, t), t), t) + 2*diff(u(x, t), t) = diff(diff(u(x, t), x), x), u(x, t), list(u(0, t) = 0, u(pi, t) = 0, u(x, 0) = sin(4*x), at(diff(u(x, t), t), t, 0) = 0))");
    let w = 15.0_f64.sqrt();
    close(value(&rhs(&s), &at), (-t).exp() * (w * t).cos() * (4.0 * x).sin() + (-t).exp() * (w * t).sin() / w * (4.0 * x).sin(), 1e-12, &s);
    // Klein-Gordon.
    let s = run("pdsolve(diff(diff(u(x, t), t), t) = diff(diff(u(x, t), x), x) - u(x, t), u(x, t), list(u(0, t) = 0, u(pi, t) = 0, u(x, 0) = sin(2*x), at(diff(u(x, t), t), t, 0) = 0))");
    close(value(&rhs(&s), &at), 5.0_f64.sqrt().mul_add(t, 0.0).cos() * (2.0 * x).sin(), 1e-12, &s);
    // Free Schrödinger particle in a box, and the Euler-Bernoulli beam.
    assert_eq!(run("pdsolve(I*diff(u(x, t), t) = -diff(diff(u(x, t), x), x), u(x, t), list(u(0, t) = 0, u(pi, t) = 0, u(x, 0) = sin(x)))"), "u(x, t) = exp(-t*I)*sin(x)");
    let beam = run("pdsolve(diff(diff(u(x, t), t), t) + diff(diff(diff(diff(u(x, t), x), x), x), x) = 0, u(x, t), list(u(0, t) = 0, u(pi, t) = 0, at(diff(diff(u(x, t), x), x), x, 0) = 0, at(diff(diff(u(x, t), x), x), x, pi) = 0, u(x, 0) = sin(2*x), at(diff(u(x, t), t), t, 0) = 0))");
    close(value(&rhs(&beam), &at), (4.0 * t).cos() * (2.0 * x).sin(), 1e-12, &beam);
}

#[test]
fn sources_and_boundary_data_on_intervals() {
    let (x, t) = (0.4, 0.1);
    let at = [("x", x), ("t", t)];
    // A steady source in the mode sin(pi x).
    let s = run(&format!("pdsolve({HEAT} + sin(pi*x), u(x, t), list(u(0, t) = 0, u(1, t) = 0, u(x, 0) = 0))"));
    close(value(&rhs(&s), &at), (1.0 - (-PI * PI * t).exp()) * (PI * x).sin() / (PI * PI), 1e-12, &s);
    // Wave with a static source: u = sin x (1 - cos t).
    let s = run("pdsolve(diff(diff(u(x, t), t), t) = diff(diff(u(x, t), x), x) + sin(x), u(x, t), list(u(0, t) = 0, u(pi, t) = 0, u(x, 0) = 0, at(diff(u(x, t), t), t, 0) = 0))");
    close(value(&rhs(&s), &at), x.sin() * (1.0 - t.cos()), 1e-12, &s);
    // A source depending on time: Duhamel's integral.
    let s = any(&format!("pdsolve({HEAT} + t*sin(pi*x), u(x, t), list(u(0, t) = 0, u(1, t) = 0, u(x, 0) = 0))"));
    let want = (PI * x).sin() * (t / (PI * PI) - (1.0 - (-PI * PI * t).exp()) / PI.powi(4));
    if !s.contains("defint") {
        close(value(&rhs(&s), &at), want, 1e-10, &s);
    }
    // Constant boundary data: u = x + series.
    let s = run(&format!("pdsolve({HEAT}, u(x, t), list(u(0, t) = 0, u(1, t) = 1, u(x, 0) = 0))"));
    let exact: f64 = x + (1..=200)
        .map(|n| {
            let k = f64::from(n) * PI;
            2.0 * (-1.0_f64).powi(n) / k * (-k * k * t).exp() * (k * x).sin()
        })
        .sum::<f64>();
    close(series(&s, &at, 40), exact, 1e-6, &s);
    // Time-dependent boundary data: the series carries the data.
    let s = run(&format!("pdsolve({HEAT}, u(x, t), list(u(0, t) = sin(t), u(1, t) = 0, u(x, 0) = 0))"));
    assert!(s.contains("sum(") && s.contains("sin(t)"), "{s}");
    // Neumann data: the flux is prescribed.
    let s = run(&format!("pdsolve({HEAT}, u(x, t), list(at(diff(u(x, t), x), x, 0) = 0, at(diff(u(x, t), x), x, 1) = 1, u(x, 0) = 0))"));
    // u_t = u_xx with u_x(1) = 1: the mean grows linearly, u = t + x²/2 - 1/6 + series.
    let mean: f64 = series(&s, &[("x", 0.5), ("t", 0.3)], 60);
    let reference = 0.3 + 0.125 - 1.0 / 6.0 + (1..=200)
        .map(|n| {
            let k = f64::from(n) * PI;
            -2.0 * (-1.0_f64).powi(n) / (k * k) * (-k * k * 0.3).exp() * (k * 0.5).cos()
        })
        .sum::<f64>();
    close(mean, reference, 1e-6, &s);
}

#[test]
fn boxes_in_two_and_three_dimensions() {
    let (x, y, t) = (0.3, 0.6, 0.04);
    let at = [("x", x), ("y", y), ("t", t)];
    let box_2d = "u(0, y, t) = 0, u(1, y, t) = 0, u(x, 0, t) = 0, u(x, 1, t) = 0";
    // Heat: a single mode.
    let s = run(&format!("pdsolve({HEAT_2D}, u(x, y, t), list({box_2d}, u(x, y, 0) = sin(pi*x)*sin(2*pi*y)))"));
    close(value(&rhs(&s), &at), (-5.0 * PI * PI * t).exp() * (PI * x).sin() * (2.0 * PI * y).sin(), 1e-12, &s);
    // Heat with general data: a double series, compared with the
    // product of one-dimensional series.
    let s = run(&format!("pdsolve({HEAT_2D}, u(x, y, t), list({box_2d}, u(x, y, 0) = x*(1 - x)*y*(1 - y)))"));
    let coefficient = |n: i32| {
        let k = f64::from(n) * PI;
        2.0 * (2.0 - 2.0 * (-1.0_f64).powi(n)) / (k * k * k)
    };
    let one_d = |z: f64| -> f64 {
        (1..=41).step_by(2).map(|n| coefficient(n) * (-(f64::from(n) * PI).powi(2) * t).exp() * (f64::from(n) * PI * z).sin()).sum()
    };
    close(series(&s, &at, 14), one_d(x) * one_d(y), 1e-6, &s);
    // Mixed conditions: Neumann in y.
    let s = run(&format!("pdsolve({HEAT_2D}, u(x, y, t), list(u(0, y, t) = 0, u(1, y, t) = 0, at(diff(u(x, y, t), y), y, 0) = 0, at(diff(u(x, y, t), y), y, 1) = 0, u(x, y, 0) = sin(pi*x)*cos(pi*y)))"));
    close(value(&rhs(&s), &at), (-2.0 * PI * PI * t).exp() * (PI * x).sin() * (PI * y).cos(), 1e-12, &s);
    // Wave on a square.
    let s = run("pdsolve(diff(diff(u(x, y, t), t), t) = diff(diff(u(x, y, t), x), x) + diff(diff(u(x, y, t), y), y), u(x, y, t), list(u(0, y, t) = 0, u(pi, y, t) = 0, u(x, 0, t) = 0, u(x, pi, t) = 0, u(x, y, 0) = sin(x)*sin(2*y), at(diff(u(x, y, t), t), t, 0) = 0))");
    close(value(&rhs(&s), &at), 5.0_f64.sqrt().mul_add(t, 0.0).cos() * x.sin() * (2.0 * y).sin(), 1e-12, &s);
    // Heat in a cube.
    let s = run("pdsolve(diff(u(x, y, z, t), t) = diff(diff(u(x, y, z, t), x), x) + diff(diff(u(x, y, z, t), y), y) + diff(diff(u(x, y, z, t), z), z), u(x, y, z, t), list(u(0, y, z, t) = 0, u(1, y, z, t) = 0, u(x, 0, z, t) = 0, u(x, 1, z, t) = 0, u(x, y, 0, t) = 0, u(x, y, 1, t) = 0, u(x, y, z, 0) = sin(pi*x)*sin(pi*y)*sin(pi*z)))");
    close(value(&rhs(&s), &[("x", x), ("y", y), ("z", 0.5), ("t", t)]), (-3.0 * PI * PI * t).exp() * (PI * x).sin() * (PI * y).sin(), 1e-12, &s);
    // A source: Poisson's equation on the unit square, u(1/2, 1/2) = -0.0736713.
    let s = run(&format!("pdsolve({LAPLACE_2D} = 1, u(x, y), list(u(0, y) = 0, u(1, y) = 0, u(x, 0) = 0, u(x, 1) = 0))"));
    close(series(&s, &[("x", 0.5), ("y", 0.5)], 25), -0.073_671_353, 5e-4, &s);
    // Helmholtz with a resonant-free source.
    let s = run(&format!("pdsolve({LAPLACE_2D} + 5*u(x, y) = sin(x)*sin(y), u(x, y), list(u(0, y) = 0, u(pi, y) = 0, u(x, 0) = 0, u(x, pi) = 0))"));
    close(value(&rhs(&s), &[("x", x), ("y", y)]), x.sin() * y.sin() / 3.0, 1e-12, &s);
}

#[test]
fn laplace_with_data_on_every_side() {
    // u(0, y) = y, u(1, y) = 0, u(x, 0) = 0, u(x, 1) = sin(pi x).
    let s = run(&format!("pdsolve({LAPLACE_2D} = 0, u(x, y), list(u(0, y) = y, u(1, y) = 0, u(x, 0) = 0, u(x, 1) = sin(pi*x)))"));
    let (x, y) = (0.5, 0.5);
    let top = (PI * y).sinh() * (PI * x).sin() / PI.sinh();
    let left: f64 = (1..=200)
        .map(|n| {
            let k = f64::from(n) * PI;
            2.0 * (-1.0_f64).powi(n + 1) / k * (k * (1.0 - x)).sinh() / k.sinh() * (k * y).sin()
        })
        .sum();
    close(series(&s, &[("x", x), ("y", y)], 30), top + left, 2e-3, &s);
}

#[test]
fn disks_cylinders_and_balls() {
    let disk_heat = "diff(u(r, t), t) = diff(diff(u(r, t), r), r) + diff(u(r, t), r)/r";
    let (r, t) = (0.4, 0.1);
    // A Bessel eigenmode of the disk decays exponentially.
    let z = 2.404_825_557_695_773;
    let s = run(&format!("pdsolve({disk_heat}, u(r, t), list(u(1, t) = 0, u(r, 0) = besselj(0, bessel_zero(0, 1)*r)))"));
    close(value(&rhs(&s), &[("r", r), ("t", t)]), (-z * z * t).exp() * crate::kernels::special::bessel_j(0.0, z * r), 1e-9, &s);
    // Neumann condition: the eigenvalues are the zeros of J_1.
    let z = 3.831_705_970_207_512;
    let s = run(&format!("pdsolve({disk_heat}, u(r, t), list(at(diff(u(r, t), r), r, 1) = 0, u(r, 0) = besselj(0, bessel_root(0, 0, 1)*r)))"));
    close(value(&rhs(&s), &[("r", r), ("t", t)]), (-z * z * t).exp() * crate::kernels::special::bessel_j(0.0, z * r), 1e-9, &s);
    // General data: a Fourier-Bessel series (the coefficients stay as integrals).
    let s = any(&format!("pdsolve({disk_heat}, u(r, t), list(u(1, t) = 0, u(r, 0) = 1 - r^2))"));
    assert!(s.contains("sum(") && s.contains("bessel_zero(0, n)") && s.contains("besselj(1, bessel_zero(0, n))"), "{s}");
    // Angular dependence: J_1 cos(theta).
    let z = 3.831_705_970_207_512;
    let s = run("pdsolve(diff(u(r, th, t), t) = diff(diff(u(r, th, t), r), r) + diff(u(r, th, t), r)/r + diff(diff(u(r, th, t), th), th)/r^2, u(r, th, t), list(u(1, th, t) = 0, u(r, th, 0) = besselj(1, bessel_zero(1, 1)*r)*cos(th)))");
    let z1 = 3.831_705_970_207_512_f64;
    let _ = z;
    close(value(&rhs(&s), &[("r", r), ("th", 0.7), ("t", t)]), (-z1 * z1 * t).exp() * crate::kernels::special::bessel_j(1.0, z1 * r) * 0.7_f64.cos(), 1e-9, &s);
    // A ball with radial symmetry: u = e^{-t} sin(r)/r on r < pi.
    let s = run("pdsolve(diff(u(r, t), t) = diff(diff(u(r, t), r), r) + 2*diff(u(r, t), r)/r, u(r, t), list(u(pi, t) = 0, u(r, 0) = sin(r)/r))");
    close(value(&rhs(&s), &[("r", 1.3), ("t", t)]), (-t).exp() * 1.3_f64.sin() / 1.3, 1e-12, &s);
    // Poisson in a disk and in a ball with zero boundary values.
    let s = run("pdsolve(diff(diff(u(r, th), r), r) + diff(u(r, th), r)/r + diff(diff(u(r, th), th), th)/r^2 = 1, u(r, th), list(u(1, th) = 0))");
    close(value(&rhs(&s), &[("r", 0.5), ("th", 0.3)]), (0.25 - 1.0) / 4.0, 1e-12, &s);
    let s = run("pdsolve(diff(diff(u(r, th), r), r) + 2*diff(u(r, th), r)/r + diff(diff(u(r, th), th), th)/r^2 + cos(th)/sin(th)*diff(u(r, th), th)/r^2 = 1, u(r, th), list(u(1, th) = 0))");
    close(value(&rhs(&s), &[("r", 0.5), ("th", 0.3)]), (0.25 - 1.0) / 6.0, 1e-12, &s);
    // Harmonic functions with Dirichlet data, in a disk and in a ball.
    let s = run("pdsolve(diff(diff(u(r, th), r), r) + diff(u(r, th), r)/r + diff(diff(u(r, th), th), th)/r^2 = 0, u(r, th), list(u(2, th) = cos(3*th)))");
    close(value(&rhs(&s), &[("r", 1.0), ("th", 0.4)]), 0.125 * (1.2_f64).cos(), 1e-12, &s);
    let s = run("pdsolve(diff(diff(u(r, th), r), r) + 2*diff(u(r, th), r)/r + diff(diff(u(r, th), th), th)/r^2 + cos(th)/sin(th)*diff(u(r, th), th)/r^2 = 0, u(r, th), list(u(1, th) = cos(th)^2))");
    let c = 0.7_f64.cos();
    close(value(&rhs(&s), &[("r", 0.5), ("th", 0.7)]), 1.0 / 3.0 + 2.0 / 3.0 * 0.25 * (1.5 * c * c - 0.5), 1e-12, &s);
    // The ring 1 < r < 2.
    let s = run("pdsolve(diff(diff(u(r, th), r), r) + diff(u(r, th), r)/r + diff(diff(u(r, th), th), th)/r^2 = 0, u(r, th), list(u(1, th) = 0, u(2, th) = cos(th)))");
    close(value(&rhs(&s), &[("r", 1.5), ("th", 0.4)]), 2.0 / 3.0 * (1.5 - 1.0 / 1.5) * 0.4_f64.cos(), 1e-12, &s);
    // A solid cylinder: Laplace's equation with data on the top.
    let s = any("pdsolve(diff(diff(u(r, z), r), r) + diff(u(r, z), r)/r + diff(diff(u(r, z), z), z) = 0, u(r, z), list(u(1, z) = 0, u(r, 0) = 0, u(r, 1) = 1 - r^2))");
    assert!(s.contains("sum(") && s.contains("bessel_zero(0, n)") && s.contains("sin("), "{s}");
}

#[test]
fn heat_on_whole_space_and_half_spaces() {
    let (x, t) = (0.7, 0.3);
    // Constant data on the half-line: the similarity solution.
    let s = run(&format!("pdsolve({HEAT}, u(x, t), list(u(0, t) = 1, u(x, 0) = 0))"));
    close(value(&rhs(&s), &[("x", x), ("t", t)]), crate::kernels::special::erfc_numerical(x / (2.0 * t.sqrt())), 1e-9, &s);
    // Constant flux into the half-line.
    let s = run(&format!("pdsolve({HEAT}, u(x, t), list(at(diff(u(x, t), x), x, 0) = -1, u(x, 0) = 0))"));
    let exact = 2.0 * (t / PI).sqrt() * (-x * x / (4.0 * t)).exp() - x * crate::kernels::special::erfc_numerical(x / (2.0 * t.sqrt()));
    close(value(&rhs(&s), &[("x", x), ("t", t)]), exact, 1e-9, &s);
    // A quadrant: zero on both axes.
    let s = run(&format!("pdsolve({HEAT_2D}, u(x, y, t), list(u(x, y, 0) = exp(-x^2 - y^2), u(0, y, t) = 0, u(x, 0, t) = 0))"));
    close(value(&rhs(&s), &[("x", 0.0), ("y", 0.5), ("t", 0.2)]), 0.0, 1e-12, &s);
    close(value(&rhs(&s), &[("x", 0.8), ("y", 0.6), ("t", 1e-6)]), (-1.0_f64).exp(), 1e-3, &s);
    // Insulated plane in two dimensions: even images.
    let s = run(&format!("pdsolve({HEAT_2D}, u(x, y, t), list(u(x, y, 0) = exp(-x^2 - y^2), at(diff(u(x, y, t), x), x, 0) = 0))"));
    assert!(!s.contains("defint"), "{s}");
    // A source in two dimensions on the whole plane.
    close(value(&rhs(&run(&format!("pdsolve({HEAT_2D} + 1, u(x, y, t), list(u(x, y, 0) = 0))"))), &[("t", 0.4)]), 0.4, 1e-12, "heated plane");
    // Drift and reaction with a half-line and a Dirichlet condition.
    let s = any("pdsolve(diff(u(x, t), t) + 2*diff(u(x, t), x) = diff(diff(u(x, t), x), x) - u(x, t), u(x, t), list(u(x, 0) = exp(-x), u(0, t) = 0))");
    assert!(s.contains("defint") || s.contains("erf"), "{s}");
}

#[test]
fn waves_in_free_space_and_half_spaces() {
    let (x, t) = (0.9, 0.7);
    // d'Alembert with a source: u_tt = u_xx + x.
    let s = run("pdsolve(diff(diff(u(x, t), t), t) = diff(diff(u(x, t), x), x) + x, u(x, t), list(u(x, 0) = 0, at(diff(u(x, t), t), t, 0) = 0))");
    close(value(&rhs(&s), &[("x", x), ("t", t)]), 0.5 * t * t * x, 1e-12, &s);
    // Poisson's formula in two dimensions.
    let w2 = "diff(diff(u(x, y, t), t), t) = diff(diff(u(x, y, t), x), x) + diff(diff(u(x, y, t), y), y)";
    let at = [("x", x), ("y", 0.4), ("t", t)];
    let s = run(&format!("pdsolve({w2}, u(x, y, t), list(u(x, y, 0) = x^2 + y^2, at(diff(u(x, y, t), t), t, 0) = 0))"));
    close(value(&rhs(&s), &at), x * x + 0.16 + 2.0 * t * t, 1e-12, &s);
    let s = run(&format!("pdsolve({w2}, u(x, y, t), list(u(x, y, 0) = 0, at(diff(u(x, y, t), t), t, 0) = x*y))"));
    close(value(&rhs(&s), &at), t * x * 0.4, 1e-12, &s);
    // Kirchhoff with a source: u_tt = Δu + 6 gives 3 t².
    let w3 = "diff(diff(u(x, y, z, t), t), t) = diff(diff(u(x, y, z, t), x), x) + diff(diff(u(x, y, z, t), y), y) + diff(diff(u(x, y, z, t), z), z)";
    let s = run(&format!("pdsolve({w3} + 6, u(x, y, z, t), list(u(x, y, z, 0) = 0, at(diff(u(x, y, z, t), t), t, 0) = 0))"));
    close(value(&rhs(&s), &[("x", 0.1), ("y", 0.2), ("z", 0.3), ("t", t)]), 3.0 * t * t, 1e-12, &s);
    // A half-space with a fixed plane: the odd extension of z³ is z³.
    let s = run(&format!("pdsolve({w3}, u(x, y, z, t), list(u(x, y, z, 0) = z^3, at(diff(u(x, y, z, t), t), t, 0) = 0, u(x, y, 0, t) = 0))"));
    close(value(&rhs(&s), &[("x", 0.1), ("y", 0.2), ("z", 0.5), ("t", t)]), 0.125 + 3.0 * t * t * 0.5, 1e-12, &s);
    // Damped waves: telegraph equation on the line, with Bessel kernels.
    let s = any("pdsolve(diff(diff(u(x, t), t), t) + 2*diff(u(x, t), t) = diff(diff(u(x, t), x), x), u(x, t), list(u(x, 0) = exp(-x^2), at(diff(u(x, t), t), t, 0) = 0))");
    assert!(s.contains("besseli(0,") && s.contains("besseli(1,") && s.contains("exp(-t)"), "{s}");
    // Klein-Gordon on the line: J_0 kernel.
    let s = any("pdsolve(diff(diff(u(x, t), t), t) - diff(diff(u(x, t), x), x) + u(x, t) = 0, u(x, t), list(u(x, 0) = 0, at(diff(u(x, t), t), t, 0) = exp(-x^2)))");
    assert!(s.contains("besselj(0,"), "{s}");
}

#[test]
fn linear_systems() {
    // First-order hyperbolic system u_t + v_x = 0, v_t + u_x = 0.
    let s = run("pdsolve(list(diff(u(x, t), t) + diff(v(x, t), x) = 0, diff(v(x, t), t) + diff(u(x, t), x) = 0), list(u(x, t), v(x, t)))");
    assert!(s.starts_with("list(u(x, t) = ") && s.contains("F_1(x - t)") && s.contains("F_2(t + x)"), "{s}");
    let s = run("pdsolve(list(diff(u(x, t), t) + diff(v(x, t), x) = 0, diff(v(x, t), t) + diff(u(x, t), x) = 0), list(u(x, t), v(x, t)), list(u(x, 0) = sin(x), v(x, 0) = 0))");
    assert_eq!(s, "list(u(x, t) = 1/2*sin(t + x) + 1/2*sin(x - t), v(x, t) = 1/2*sin(x - t) - 1/2*sin(t + x))");
    // Coupled heat equations: decoupled by the eigenvectors of [[2, 1], [1, 2]].
    let s = run("pdsolve(list(diff(u(x, t), t) = 2*diff(diff(u(x, t), x), x) + diff(diff(v(x, t), x), x), diff(v(x, t), t) = diff(diff(u(x, t), x), x) + 2*diff(diff(v(x, t), x), x)), list(u(x, t), v(x, t)), list(u(0, t) = 0, u(pi, t) = 0, v(0, t) = 0, v(pi, t) = 0, u(x, 0) = sin(x), v(x, 0) = 0))");
    assert_eq!(s, "list(u(x, t) = 1/2*exp(-3*t)*sin(x) + 1/2*exp(-t)*sin(x), v(x, t) = 1/2*exp(-3*t)*sin(x) - 1/2*exp(-t)*sin(x))");
    // Coupled waves with a reaction coupling.
    let s = run("pdsolve(list(diff(diff(u(x, t), t), t) = 3*diff(diff(u(x, t), x), x) + diff(diff(v(x, t), x), x), diff(diff(v(x, t), t), t) = diff(diff(u(x, t), x), x) + 3*diff(diff(v(x, t), x), x)), list(u(x, t), v(x, t)), list(u(x, 0) = sin(x), v(x, 0) = 0, at(diff(u(x, t), t), t, 0) = 0, at(diff(v(x, t), t), t, 0) = 0))");
    let (x, t): (f64, f64) = (0.4, 0.9);
    let u = 0.25 * ((x - 2.0 * t).sin() + (x - 2.0_f64.sqrt() * t).sin() + (x + 2.0 * t).sin() + (x + 2.0_f64.sqrt() * t).sin());
    let inner = s.trim_start_matches("list(u(x, t) = ");
    let (u_text, _) = inner.split_once(", v(x, t) = ").unwrap_or((inner, ""));
    close(value(u_text, &[("x", x), ("t", t)]), u, 1e-12, &s);
    // Coupled heat equations in two dimensions on a square, with a coupling term.
    let s = run("pdsolve(list(diff(u(x, y, t), t) = diff(diff(u(x, y, t), x), x) + diff(diff(u(x, y, t), y), y) + v(x, y, t), diff(v(x, y, t), t) = diff(diff(v(x, y, t), x), x) + diff(diff(v(x, y, t), y), y) + u(x, y, t)), list(u(x, y, t), v(x, y, t)), list(u(0, y, t) = 0, u(pi, y, t) = 0, u(x, 0, t) = 0, u(x, pi, t) = 0, v(0, y, t) = 0, v(pi, y, t) = 0, v(x, 0, t) = 0, v(x, pi, t) = 0, u(x, y, 0) = sin(x)*sin(y), v(x, y, 0) = sin(x)*sin(y)))");
    assert!(s.starts_with("list(u(x, y, t) = "), "{s}");
    let (x, y, t) = (0.7, 0.9, 0.2);
    // u = v = e^{(1 - 2) t} sin x sin y.
    let u_text = s.trim_start_matches("list(u(x, y, t) = ").split(", v(x, y, t) = ").next().unwrap_or("");
    close(value(u_text, &[("x", x), ("y", y), ("t", t)]), (-t).exp() * x.sin() * y.sin(), 1e-12, &s);
}

#[test]
fn elliptic_problems_on_half_spaces() {
    // The upper half-plane with Dirichlet data: the Poisson kernel.
    let s = any(&format!("pdsolve({LAPLACE_2D} = 0, u(x, y), list(u(x, 0) = 1/(1 + x^2)))"));
    assert!(s.contains("defint(") && s.contains("pi"), "{s}");
    // Poisson's equation with a source and a grounded plane: images.
    let s = any(&format!("pdsolve({LAPLACE_2D} = exp(-x^2 - y^2), u(x, y), list(u(x, 0) = 0))"));
    assert!(s.contains("ln("), "{s}");
    // Three dimensions: Neumann data, the single layer.
    let lap3 = "diff(diff(u(x, y, z), x), x) + diff(diff(u(x, y, z), y), y) + diff(diff(u(x, y, z), z), z)";
    let s = any(&format!("pdsolve({lap3} = 0, u(x, y, z), list(at(diff(u(x, y, z), z), z, 0) = exp(-x^2 - y^2)))"));
    assert!(s.contains("defint(") && s.contains("pi"), "{s}");
}

#[test]
fn drift_on_intervals_and_plates() {
    // u_t + 2 u_x = u_xx on [0, pi] with zero ends: u = e^{x - 2t} sin x.
    let (x, t) = (0.8, 0.3);
    let s = run("pdsolve(diff(u(x, t), t) + 2*diff(u(x, t), x) = diff(diff(u(x, t), x), x), u(x, t), list(u(0, t) = 0, u(pi, t) = 0, u(x, 0) = exp(x)*sin(x)))");
    close(value(&rhs(&s), &[("x", x), ("t", t)]), (x - 2.0 * t).exp() * x.sin(), 1e-12, &s);
    // The Navier plate: biharmonic with simply supported edges.
    let d4 = "diff(diff(diff(diff(u(x, y), x), x), x), x) + 2*diff(diff(diff(diff(u(x, y), x), x), y), y) + diff(diff(diff(diff(u(x, y), y), y), y), y)";
    let edges = "u(0, y) = 0, u(pi, y) = 0, u(x, 0) = 0, u(x, pi) = 0, at(diff(diff(u(x, y), x), x), x, 0) = 0, at(diff(diff(u(x, y), x), x), x, pi) = 0, at(diff(diff(u(x, y), y), y), y, 0) = 0, at(diff(diff(u(x, y), y), y), y, pi) = 0";
    let s = run(&format!("pdsolve({d4} = sin(x)*sin(y), u(x, y), list({edges}))"));
    close(value(&rhs(&s), &[("x", 0.7), ("y", 1.1)]), 0.7_f64.sin() * 1.1_f64.sin() / 4.0, 1e-12, &s);
}

/// Replaces every `defint(body, var, lo, hi)` by the midpoint rule with
/// `points` nodes (the integrands here are analytic).
fn expand_integrals(
    text: &str,
    points: usize,
) -> String {
    let mut s = text.to_owned();
    while let Some(start) = s.rfind("defint(") {
        let open = start + 6;
        let mut depth = 0_i32;
        let mut end = open;
        for (i, c) in s.char_indices().skip(open) {
            match c {
                | '(' => depth += 1,
                | ')' => {
                    depth -= 1;
                    if depth == 0 {
                        end = i;
                        break;
                    }
                },
                | _ => {},
            }
        }
        let inner = s[open + 1..end].to_owned();
        let mut parts = Vec::new();
        let (mut depth, mut last) = (0_i32, 0);
        for (i, c) in inner.char_indices() {
            match c {
                | '(' => depth += 1,
                | ')' => depth -= 1,
                | ',' if depth == 0 => {
                    parts.push(inner[last..i].trim().to_owned());
                    last = i + 1;
                },
                | _ => {},
            }
        }
        parts.push(inner[last..].trim().to_owned());
        assert_eq!(parts.len(), 4, "malformed integral in {text}");
        let (body, var, lo, hi) = (&parts[0], &parts[1], &parts[2], &parts[3]);
        let pieces: Vec<String> = (0..points)
            .map(|k| {
                let fraction = (f64::from(u32::try_from(k).unwrap_or(0)) + 0.5) / f64::from(u32::try_from(points).unwrap_or(1));
                let node = format!("(({lo}) + (({hi}) - ({lo}))*{fraction})");
                replace_identifier(body, var, &node)
            })
            .collect();
        s.replace_range(start..=end, &format!("((({hi}) - ({lo}))/{points}*({}))", pieces.join(" + ")));
    }
    s
}

/// `∫_{x-t}^{x+t} f(s) K(R(s)) ds` by the midpoint rule, `R² = t² - (x - s)²`.
fn riemann(
    x: f64,
    t: f64,
    f: &dyn Fn(f64) -> f64,
    kernel: &dyn Fn(f64) -> f64,
) -> f64 {
    let n = 4000;
    let h = 2.0 * t / f64::from(n);
    (0..n)
        .map(|k| {
            let s = x - t + h * (f64::from(k) + 0.5);
            let r = (t * t - (x - s) * (x - s)).max(0.0).sqrt();
            f(s) * kernel(r) * h
        })
        .sum()
}

#[test]
fn klein_gordon_and_telegraph_on_the_line() {
    use crate::kernels::special::bessel_i;
    use crate::kernels::special::bessel_j;
    let gauss = |s: f64| (-s * s).exp();
    // Klein-Gordon u_tt = u_xx - u, u(x, 0) = 0, u_t(x, 0) = e^{-x²}.
    let kg = |x: f64, t: f64| 0.5 * riemann(x, t, &gauss, &|r| bessel_j(0.0, r));
    let (h, x, t) = (0.02, 0.4, 0.8);
    let residual = (kg(x, t + h) - 2.0 * kg(x, t) + kg(x, t - h)) / (h * h) - (kg(x + h, t) - 2.0 * kg(x, t) + kg(x - h, t)) / (h * h) + kg(x, t);
    assert!(residual.abs() < 2e-3, "the Riemann formula does not solve the equation: {residual}");
    let s = any("pdsolve(diff(diff(u(x, t), t), t) - diff(diff(u(x, t), x), x) + u(x, t) = 0, u(x, t), list(u(x, 0) = 0, at(diff(u(x, t), t), t, 0) = exp(-x^2)))");
    close(value(&expand_integrals(&rhs(&s), 96), &[("x", x), ("t", t)]), kg(x, t), 1e-5, &s);
    // Telegraph equation u_tt + 2 u_t = u_xx with u(x, 0) = e^{-x²}: u = e^{-t} w,
    // w_tt = w_xx + w, w(0) = f, w_t(0) = f.
    let w = |x: f64, t: f64| {
        0.5 * (gauss(x + t) + gauss(x - t))
            + 0.5 * riemann(x, t, &gauss, &|r| bessel_i(0.0, r))
            + 0.5 * t * riemann(x, t, &gauss, &|r| if r < 1e-9 { 0.5 } else { bessel_i(1.0, r) / r })
    };
    let u = |x: f64, t: f64| (-t).exp() * w(x, t);
    let residual = (u(x, t + h) - 2.0 * u(x, t) + u(x, t - h)) / (h * h) + 2.0 * (u(x, t + h) - u(x, t - h)) / (2.0 * h)
        - (u(x + h, t) - 2.0 * u(x, t) + u(x - h, t)) / (h * h);
    assert!(residual.abs() < 2e-3, "the damped Riemann formula does not solve the equation: {residual}");
    let s = any("pdsolve(diff(diff(u(x, t), t), t) + 2*diff(u(x, t), t) = diff(diff(u(x, t), x), x), u(x, t), list(u(x, 0) = exp(-x^2), at(diff(u(x, t), t), t, 0) = 0))");
    close(value(&expand_integrals(&rhs(&s), 96), &[("x", x), ("t", t)]), u(x, t), 1e-5, &s);
}

#[test]
fn robin_disks_neumann_balls_and_steady_faces() {
    use crate::kernels::special::bessel_j;
    let (r, t) = (0.4, 0.1);
    // Robin condition on a disk: z J_0'(z) + 2 J_0(z) = 0.
    let s = run("pdsolve(diff(u(r, t), t) = diff(diff(u(r, t), r), r) + diff(u(r, t), r)/r, u(r, t), list(at(diff(u(r, t), r), r, 1) + 2*u(1, t) = 0, u(r, 0) = besselj(0, bessel_root(0, 2, 1)*r)))");
    let z = value("bessel_root(0, 2, 1)", &[]);
    close(-z * bessel_j(1.0, z) + 2.0 * bessel_j(0.0, z), 0.0, 1e-8, "robin root");
    close(value(&rhs(&s), &[("r", r), ("t", t)]), (-z * z * t).exp() * bessel_j(0.0, z * r), 1e-9, &s);
    // Neumann condition on a ball: tan z = z.
    let s = run("pdsolve(diff(u(r, t), t) = diff(diff(u(r, t), r), r) + 2*diff(u(r, t), r)/r, u(r, t), list(at(diff(u(r, t), r), r, 1) = 0, u(r, 0) = sin(bessel_root(1/2, -1/2, 1)*r)/r))");
    let z = value("bessel_root(1/2, -1/2, 1)", &[]);
    close(z.tan() - z, 0.0, 1e-8, "ball root");
    close(value(&rhs(&s), &[("r", r), ("t", t)]), (-z * z * t).exp() * (z * r).sin() / r, 1e-9, &s);
    // Laplace's equation on a rectangle with an insulated edge: cosh profile.
    let s = run(&format!("pdsolve({LAPLACE_2D} = 0, u(x, y), list(u(0, y) = 0, u(1, y) = 0, at(diff(u(x, y), y), y, 0) = 0, u(x, 1) = sin(pi*x)))"));
    close(value(&rhs(&s), &[("x", 0.3), ("y", 0.6)]), (PI * 0.6).cosh() * (PI * 0.3).sin() / PI.cosh(), 1e-12, &s);
    // A harmonic function in a box with data on two opposite faces.
    let lap3 = "diff(diff(u(x, y, z), x), x) + diff(diff(u(x, y, z), y), y) + diff(diff(u(x, y, z), z), z)";
    let s = run(&format!(
        "pdsolve({lap3} = 0, u(x, y, z), list(u(0, y, z) = 0, u(1, y, z) = 0, u(x, 0, z) = 0, u(x, 1, z) = 0, u(x, y, 0) = sin(pi*x)*sin(pi*y), u(x, y, 1) = 2*sin(pi*x)*sin(pi*y)))"
    ));
    let k = 2.0_f64.sqrt() * PI;
    let profile = |z: f64| ((k * (1.0 - z)).sinh() + 2.0 * (k * z).sinh()) / k.sinh();
    close(value(&rhs(&s), &[("x", 0.3), ("y", 0.45), ("z", 0.7)]), profile(0.7) * (PI * 0.3).sin() * (PI * 0.45).sin(), 1e-12, &s);
}

#[test]
fn schroedinger_and_three_unknowns() {
    // The free particle on the half-line with a node at the wall.
    let s = any("pdsolve(I*diff(u(x, t), t) = -diff(diff(u(x, t), x), x), u(x, t), list(u(x, 0) = exp(-x^2), u(0, t) = 0))");
    assert!(s.contains("defint") || s.contains("erf"), "{s}");
    // A 3x3 first-order system: U_t + A U_x = 0, A = tridiagonal(1, 2, 1).
    let src = "pdsolve(list(diff(u(x, t), t) + 2*diff(u(x, t), x) + diff(v(x, t), x) = 0, diff(v(x, t), t) + diff(u(x, t), x) + 2*diff(v(x, t), x) + diff(w(x, t), x) = 0, diff(w(x, t), t) + diff(v(x, t), x) + 2*diff(w(x, t), x) = 0), list(u(x, t), v(x, t), w(x, t)), list(u(x, 0) = sin(x), v(x, 0) = 0, w(x, 0) = 0))";
    let s = run(src);
    let inner = s.trim_start_matches("list(").strip_suffix(')').unwrap_or("");
    let a: Vec<&str> = inner.split(", v(x, t) = ").collect();
    assert_eq!(a.len(), 2, "{s}");
    let b: Vec<&str> = a[1].split(", w(x, t) = ").collect();
    assert_eq!(b.len(), 2, "{s}");
    let parts = [a[0].trim_start_matches("u(x, t) = ").to_owned(), b[0].to_owned(), b[1].to_owned()];
    let at = |i: usize, x: f64, t: f64| value(&parts[i], &[("x", x), ("t", t)]);
    let (x, t, h) = (0.5, 0.3, 1e-4);
    let (u, v, w) = (|x, t| at(0, x, t), |x, t| at(1, x, t), |x, t| at(2, x, t));
    let dx = |f: &dyn Fn(f64, f64) -> f64| (f(x + h, t) - f(x - h, t)) / (2.0 * h);
    let dt = |f: &dyn Fn(f64, f64) -> f64| (f(x, t + h) - f(x, t - h)) / (2.0 * h);
    close(dt(&u) + 2.0 * dx(&u) + dx(&v), 0.0, 1e-6, &s);
    close(dt(&v) + dx(&u) + 2.0 * dx(&v) + dx(&w), 0.0, 1e-6, &s);
    close(dt(&w) + dx(&v) + 2.0 * dx(&w), 0.0, 1e-6, &s);
    close(u(x, 0.0), x.sin(), 1e-12, &s);
    close(v(x, 0.0), 0.0, 1e-12, &s);
}

#[test]
fn transport_with_data_and_step_data() {
    // u_t + u_x + 2 u_y + 3 u_z + u = 0 with a Gaussian: transported and damped.
    let s = run("pdsolve(diff(u(x, y, z, t), t) + diff(u(x, y, z, t), x) + 2*diff(u(x, y, z, t), y) + 3*diff(u(x, y, z, t), z) + u(x, y, z, t) = 0, u(x, y, z, t), list(u(x, y, z, 0) = exp(-x^2 - y^2 - z^2)))");
    let (x, y, z, t): (f64, f64, f64, f64) = (0.3, 0.5, 0.7, 0.4);
    let exact = (-t).exp() * (-(x - t).powi(2) - (y - 2.0 * t).powi(2) - (z - 3.0 * t).powi(2)).exp();
    close(value(&rhs(&s), &[("x", x), ("y", y), ("z", z), ("t", t)]), exact, 1e-12, &s);
    // A constant source along the characteristics.
    let s = run("pdsolve(diff(u(x, y, z, t), t) + diff(u(x, y, z, t), x) + diff(u(x, y, z, t), y) + diff(u(x, y, z, t), z) = 1, u(x, y, z, t), list(u(x, y, z, 0) = 0))");
    close(value(&rhs(&s), &[("x", x), ("y", y), ("z", z), ("t", t)]), t, 1e-12, &s);
    // Step data on the line: u = erfc(-x / 2√t) / 2.
    let s = any(&format!("pdsolve({HEAT}, u(x, t), list(u(x, 0) = heaviside(x)))"));
    if !s.contains("defint") {
        let want = 0.5 * crate::kernels::special::erfc_numerical(-0.3 / (2.0 * 0.5_f64.sqrt()));
        close(value(&rhs(&s), &[("x", 0.3), ("t", 0.5)]), want, 1e-8, &s);
    }
}

#[test]
fn rational_three_by_three_system() {
    // U_t + A U_x = 0 with the upper triangular A = [[1, 1, 0], [0, 2, 1], [0, 0, 3]].
    let src = "pdsolve(list(diff(u(x, t), t) + diff(u(x, t), x) + diff(v(x, t), x) = 0, diff(v(x, t), t) + 2*diff(v(x, t), x) + diff(w(x, t), x) = 0, diff(w(x, t), t) + 3*diff(w(x, t), x) = 0), list(u(x, t), v(x, t), w(x, t)), list(u(x, 0) = 0, v(x, 0) = 0, w(x, 0) = sin(x)))";
    let s = run(src);
    let inner = s.trim_start_matches("list(").strip_suffix(')').unwrap_or("");
    let a: Vec<&str> = inner.split(", v(x, t) = ").collect();
    assert_eq!(a.len(), 2, "{s}");
    let b: Vec<&str> = a[1].split(", w(x, t) = ").collect();
    assert_eq!(b.len(), 2, "{s}");
    let parts = [a[0].trim_start_matches("u(x, t) = ").to_owned(), b[0].to_owned(), b[1].to_owned()];
    let at = |i: usize, x: f64, t: f64| value(&parts[i], &[("x", x), ("t", t)]);
    let (x, t, h) = (0.5, 0.3, 1e-4);
    let dx = |i: usize| (at(i, x + h, t) - at(i, x - h, t)) / (2.0 * h);
    let dt = |i: usize| (at(i, x, t + h) - at(i, x, t - h)) / (2.0 * h);
    close(dt(0) + dx(0) + dx(1), 0.0, 1e-6, &s);
    close(dt(1) + 2.0 * dx(1) + dx(2), 0.0, 1e-6, &s);
    close(dt(2) + 3.0 * dx(2), 0.0, 1e-6, &s);
    close(at(2, x, 0.0), x.sin(), 1e-12, &s);
    close(at(0, x, 0.0), 0.0, 1e-12, &s);
}
