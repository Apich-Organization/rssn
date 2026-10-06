//! # Additional special functions
//!
//! Complete and incomplete elliptic integrals (Carlson symmetric forms and
//! the arithmetic-geometric mean), Gauss hypergeometric `2F1`, Kummer
//! `1F1`, Bessel functions of integer order by spectrally accurate
//! trapezoid quadrature, and Jacobi elliptic functions.
//!
//! Elliptic integrals use the *parameter* `m = k^2`.
#![allow(
    clippy::manual_midpoint,
    clippy::missing_const_for_fn,
    clippy::struct_field_names,
    clippy::or_fun_call,
    clippy::manual_swap,
    clippy::if_not_else,
    clippy::unnecessary_sort_by,
    clippy::while_float,
    clippy::too_long_first_doc_paragraph,
    clippy::cast_sign_loss,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::needless_pass_by_value,
    clippy::manual_map,
    clippy::unnecessary_map_or,
    clippy::suboptimal_flops,
    clippy::similar_names,
    clippy::unreadable_literal,
    clippy::excessive_precision,
    clippy::needless_range_loop,
    clippy::float_cmp,
    clippy::too_many_lines,
    clippy::cognitive_complexity,
    clippy::option_if_let_else,
    clippy::many_single_char_names,
    clippy::type_complexity,
    clippy::too_many_arguments,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_precision_loss,
    clippy::map_unwrap_or,
    clippy::cast_lossless
)]

use std::f64::consts::PI;

use crate::kernels::special::gamma_numerical;

/// Carlson's symmetric elliptic integral `R_F(x, y, z)` (arguments non-negative, at most one zero).
#[must_use]
pub fn carlson_rf(x: f64, y: f64, z: f64) -> f64 {
    let (mut x, mut y, mut z) = (x, y, z);
    let (mut delx, mut dely, mut delz, mut ave);
    loop {
        let (sx, sy, sz) = (x.sqrt(), y.sqrt(), z.sqrt());
        let lam = sx * (sy + sz) + sy * sz;
        x = 0.25 * (x + lam);
        y = 0.25 * (y + lam);
        z = 0.25 * (z + lam);
        ave = (x + y + z) / 3.0;
        delx = (ave - x) / ave;
        dely = (ave - y) / ave;
        delz = (ave - z) / ave;
        if delx.abs().max(dely.abs()).max(delz.abs()) < 1e-3 {
            break;
        }
    }
    let e2 = delx * dely - delz.powi(2);
    let e3 = delx * dely * delz;
    (1.0 + (e2 / 24.0 - 0.1 - 3.0 * e3 / 44.0) * e2 + e3 / 14.0) / ave.sqrt()
}

/// Carlson's `R_D(x, y, z)` (`z > 0`).
#[must_use]
pub fn carlson_rd(x: f64, y: f64, z: f64) -> f64 {
    let (mut x, mut y, mut z) = (x, y, z);
    let (mut sum, mut fac) = (0.0, 1.0);
    let (mut delx, mut dely, mut delz, mut ave);
    loop {
        let (sx, sy, sz) = (x.sqrt(), y.sqrt(), z.sqrt());
        let lam = sx * (sy + sz) + sy * sz;
        sum += fac / (sz * (z + lam));
        fac *= 0.25;
        x = 0.25 * (x + lam);
        y = 0.25 * (y + lam);
        z = 0.25 * (z + lam);
        ave = 0.2 * (x + y + 3.0 * z);
        delx = (ave - x) / ave;
        dely = (ave - y) / ave;
        delz = (ave - z) / ave;
        if delx.abs().max(dely.abs()).max(delz.abs()) < 1e-3 {
            break;
        }
    }
    let (c1, c2, c3, c4) = (3.0 / 14.0, 1.0 / 6.0, 9.0 / 22.0, 3.0 / 26.0);
    let (c5, c6) = (0.25 * c3, 1.5 * c4);
    let ea = delx * dely;
    let eb = delz * delz;
    let ec = ea - eb;
    let ed = ea - 6.0 * eb;
    let ee = ed + ec + ec;
    3.0 * sum
        + fac
            * (1.0 + ed * (-c1 + c5 * ed - c6 * delz * ee)
                + delz * (c2 * ee + delz * (-c3 * ec + delz * c4 * ea)))
            / (ave * ave.sqrt())
}

/// Complete elliptic integral of the first kind `K(m)`, `m < 1`.
#[must_use]
pub fn elliptic_k(m: f64) -> f64 {
    if m == 1.0 {
        return f64::INFINITY;
    }
    if m > 1.0 || m.is_nan() {
        return f64::NAN;
    }
    // AGM: K = pi / (2 agm(1, sqrt(1-m)))
    let (mut a, mut b) = (1.0, (1.0 - m).sqrt());
    for _ in 0..64 {
        let an = 0.5 * (a + b);
        b = (a * b).sqrt();
        a = an;
        if (a - b).abs() <= 1e-16 * a {
            break;
        }
    }
    PI / (2.0 * a)
}

/// Complete elliptic integral of the second kind `E(m)`, `m <= 1`.
#[must_use]
pub fn elliptic_e(m: f64) -> f64 {
    if m == 1.0 {
        return 1.0;
    }
    if m > 1.0 || m.is_nan() {
        return f64::NAN;
    }
    carlson_rf(0.0, 1.0 - m, 1.0) - m / 3.0 * carlson_rd(0.0, 1.0 - m, 1.0)
}

/// Incomplete elliptic integral of the first kind `F(phi, m)`.
#[must_use]
pub fn elliptic_f(phi: f64, m: f64) -> f64 {
    // reduce to |phi| <= pi/2 using F(phi + n pi) = F(phi) + 2 n K
    let n = (phi / PI).round();
    let r = phi - n * PI;
    let (s, c) = r.sin_cos();
    let v = s * carlson_rf(c * c, 1.0 - m * s * s, 1.0);
    v + 2.0 * n * elliptic_k(m)
}

/// Incomplete elliptic integral of the second kind `E(phi, m)`.
#[must_use]
pub fn elliptic_e_inc(phi: f64, m: f64) -> f64 {
    let n = (phi / PI).round();
    let r = phi - n * PI;
    let (s, c) = r.sin_cos();
    let q = 1.0 - m * s * s;
    let v = s * carlson_rf(c * c, q, 1.0) - m / 3.0 * s * s * s * carlson_rd(c * c, q, 1.0);
    v + 2.0 * n * elliptic_e(m)
}

/// Jacobi elliptic functions `(sn, cn, dn)` of `u` with parameter `m` in
/// `[0, 1]` (descending Landen / AGM method).
#[must_use]
pub fn jacobi_elliptic(u: f64, m: f64) -> (f64, f64, f64) {
    if m < 1e-15 {
        return (u.sin(), u.cos(), 1.0);
    }
    if m > 1.0 - 1e-15 {
        let t = u.tanh();
        let s = 1.0 / u.cosh();
        return (t, s, s);
    }
    let mut a = [0.0; 40];
    let mut c = [0.0; 40];
    a[0] = 1.0;
    let mut b = (1.0 - m).sqrt();
    c[0] = m.sqrt();
    let mut n = 0;
    while n < 39 && (c[n]).abs() > 1e-16 {
        a[n + 1] = 0.5 * (a[n] + b);
        c[n + 1] = 0.5 * (a[n] - b);
        b = (a[n] * b).sqrt();
        n += 1;
    }
    let mut phi = (2.0_f64).powi(n as i32) * a[n] * u;
    let mut i = n;
    while i > 0 {
        phi = 0.5 * (phi + (c[i] * phi.sin() / a[i]).asin());
        i -= 1;
    }
    let sn = phi.sin();
    let cn = phi.cos();
    (sn, cn, (1.0 - m * sn * sn).sqrt())
}

/// Gauss hypergeometric function `2F1(a, b; c; z)` for real `z < 1`.
///
/// Uses the power series for `|z| <= 0.9`, the Pfaff transformation for
/// negative `z`, and the `1 - z` connection formula (when the gamma
/// factors are finite) near `z = 1`. Returns NaN for `z >= 1`.
#[must_use]
pub fn hyp2f1(a: f64, b: f64, c: f64, z: f64) -> f64 {
    if z >= 1.0 || z.is_nan() {
        return f64::NAN;
    }
    let series = |a: f64, b: f64, c: f64, z: f64| -> f64 {
        let (mut term, mut sum) = (1.0, 1.0);
        for n in 0..200_000 {
            let nf = f64::from(n);
            term *= (a + nf) * (b + nf) / ((c + nf) * (nf + 1.0)) * z;
            sum += term;
            if term == 0.0 || term.abs() < 1e-17 * sum.abs() {
                break;
            }
        }
        sum
    };
    if z < -0.5 {
        let w = z / (z - 1.0);
        return (1.0 - z).powf(-a) * hyp2f1(a, c - b, c, w);
    }
    if z > 0.9 {
        let s = c - a - b;
        if (s - s.round()).abs() > 1e-9 {
            let v = gamma_numerical(c) * gamma_numerical(s)
                / (gamma_numerical(c - a) * gamma_numerical(c - b))
                * series(a, b, a + b - c + 1.0, 1.0 - z)
                + (1.0 - z).powf(s) * gamma_numerical(c) * gamma_numerical(-s)
                    / (gamma_numerical(a) * gamma_numerical(b))
                    * series(c - a, c - b, s + 1.0, 1.0 - z);
            if v.is_finite() {
                return v;
            }
        }
    }
    series(a, b, c, z)
}

/// Kummer confluent hypergeometric function `1F1(a; b; x)`.
#[must_use]
pub fn hyp1f1(a: f64, b: f64, x: f64) -> f64 {
    if x < 0.0 {
        return x.exp() * hyp1f1(b - a, b, -x);
    }
    let (mut term, mut sum) = (1.0, 1.0);
    for n in 0..500_000 {
        let nf = f64::from(n);
        term *= (a + nf) / ((b + nf) * (nf + 1.0)) * x;
        sum += term;
        if term == 0.0 || (term.abs() < 1e-17 * sum.abs() && nf > x) {
            break;
        }
    }
    sum
}

/// Bessel function of the first kind `J_n(x)` for integer order, via the
/// integral `(1/pi) * int_0^pi cos(n t - x sin t) dt` evaluated with the
/// trapezoid rule (geometrically convergent for periodic integrands).
#[must_use]
pub fn bessel_jn(n: i32, x: f64) -> f64 {
    let nn = n.unsigned_abs() as usize;
    let pts = 2 * (x.abs().ceil() as usize + nn + 40);
    let mut s = 0.0;
    for k in 0..=pts {
        let t = PI * k as f64 / pts as f64;
        let w = if k == 0 || k == pts { 0.5 } else { 1.0 };
        s += w * (f64::from(n) * t - x * t.sin()).cos();
    }
    s / pts as f64
}

/// Modified Bessel function of the first kind `I_n(x)` for integer order
/// (valid up to `|x| ~ 700`).
#[must_use]
pub fn bessel_in(n: i32, x: f64) -> f64 {
    let nn = n.unsigned_abs() as usize;
    let pts = 2 * (x.abs().ceil() as usize + nn + 40);
    let mut s = 0.0;
    for k in 0..=pts {
        let t = PI * k as f64 / pts as f64;
        let w = if k == 0 || k == pts { 0.5 } else { 1.0 };
        s += w * (x * t.cos()).exp() * (f64::from(n) * t).cos();
    }
    s / pts as f64
}
