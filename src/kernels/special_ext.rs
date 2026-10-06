//! # Additional special functions
//!
//! Complete and incomplete elliptic integrals (Carlson symmetric forms and
//! the arithmetic-geometric mean), Gauss hypergeometric `2F1`, Kummer
//! `1F1`, Bessel functions of integer order (series, Hankel asymptotics and
//! Miller backward recurrence) and Jacobi elliptic functions with exact
//! period reduction.
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

use crate::kernels::special::{digamma_numerical, gamma_numerical, ln_gamma_numerical};

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

/// Jacobi elliptic functions `(sn, cn, dn)` of `u` with parameter
/// `m = k^2` in `[0, 1]`; NaNs outside that range.
///
/// The argument is first reduced modulo the real period `4K(m)` and
/// folded with the quarter-period symmetries
/// (`sn(2K - v) = sn(v)`, `cn(2K - v) = -cn(v)`, `dn(2K - v) = dn(v)`) so
/// that the descending Landen / AGM iteration only ever sees an amplitude
/// in `[0, pi/2]`. The identities `sn^2 + cn^2 = 1` and
/// `dn^2 + m sn^2 = 1` therefore hold to machine precision for any
/// `|u|`, and the accuracy is limited only by the representation of `u`
/// itself (about `|u| * 1e-16`).
#[must_use]
pub fn jacobi_elliptic(u: f64, m: f64) -> (f64, f64, f64) {
    if !u.is_finite() || !(0.0..=1.0).contains(&m) {
        return (f64::NAN, f64::NAN, f64::NAN);
    }
    if m < 1e-15 {
        return (u.sin(), u.cos(), 1.0);
    }
    if m > 1.0 - 1e-15 {
        let t = u.tanh();
        let s = 1.0 / u.cosh();
        return (t, s, s);
    }
    let kk = elliptic_k(m);
    let period = 4.0 * kk;
    let r = u - period * (u / period).round();
    let a = r.abs();
    let (v, flip) = if a > kk { (2.0 * kk - a, true) } else { (a, false) };
    let (s, c) = jacobi_landen(v, m);
    let sn = if r < 0.0 { -s } else { s };
    let cn = if flip { -c } else { c };
    let dn = ((1.0 - m) + m * c * c).sqrt();
    (sn, cn, dn)
}

/// `(sin am, cos am)` for `u` in `[0, K]` by descending Landen iteration.
fn jacobi_landen(u: f64, m: f64) -> (f64, f64) {
    let mut a = [0.0; 40];
    let mut c = [0.0; 40];
    a[0] = 1.0;
    let mut b = (1.0 - m).sqrt();
    c[0] = m.sqrt();
    let mut n = 0;
    while n < 39 && c[n].abs() > 1e-16 {
        a[n + 1] = 0.5 * (a[n] + b);
        c[n + 1] = 0.5 * (a[n] - b);
        b = (a[n] * b).sqrt();
        n += 1;
    }
    let mut phi = (2.0_f64).powi(n as i32) * a[n] * u;
    let mut i = n;
    while i > 0 {
        let arg = (c[i] * phi.sin() / a[i]).clamp(-1.0, 1.0);
        phi = 0.5 * (phi + arg.asin());
        i -= 1;
    }
    phi.sin_cos()
}

fn is_nonpos_int(x: f64) -> bool {
    x <= 0.0 && (x - x.round()).abs() <= 1e-13 * x.abs().max(1.0)
}

/// `(ln|Gamma(x)|, sign)`; a pole gives `(inf, 1)`.
fn lgamma_signed(x: f64) -> (f64, f64) {
    if x > 0.0 {
        (ln_gamma_numerical(x), 1.0)
    } else if is_nonpos_int(x) {
        (f64::INFINITY, 1.0)
    } else {
        let s = (PI * x).sin();
        ((PI / s.abs()).ln() - ln_gamma_numerical(1.0 - x), s.signum())
    }
}

/// `prod Gamma(num) / prod Gamma(den)` evaluated through logarithms
/// (poles in the denominator give 0, poles in the numerator infinity).
fn gamma_ratio(num: &[f64], den: &[f64]) -> f64 {
    let mut ln = 0.0;
    let mut sign = 1.0;
    for &x in num {
        let (l, s) = lgamma_signed(x);
        if l.is_infinite() {
            return f64::INFINITY * sign;
        }
        ln += l;
        sign *= s;
    }
    for &x in den {
        let (l, s) = lgamma_signed(x);
        if l.is_infinite() {
            return 0.0;
        }
        ln -= l;
        sign *= s;
    }
    sign * ln.exp()
}

/// Plain power series of `2F1` (also the finite polynomial when a
/// numerator parameter is a non-positive integer).
fn hyp2f1_series(a: f64, b: f64, c: f64, z: f64) -> f64 {
    let (mut term, mut sum) = (1.0, 1.0);
    for n in 0..500_000 {
        let nf = f64::from(n);
        let ratio = (a + nf) * (b + nf) / ((c + nf) * (nf + 1.0)) * z;
        term *= ratio;
        sum += term;
        if term == 0.0 || (ratio.abs() < 1.0 && term.abs() < 1e-17 * sum.abs()) {
            break;
        }
    }
    sum
}

/// `2F1(a, b; a + b + m; 1 - y)` for an integer `m >= 0` and `0 < y < 1`
/// (Abramowitz-Stegun 15.3.10 and 15.3.11 with the logarithmic terms).
fn hyp2f1_log_case(a: f64, b: f64, m: u32, y: f64) -> f64 {
    let c = a + b + f64::from(m);
    let ly = y.ln();
    let euler = -0.577_215_664_901_532_9;
    if m == 0 {
        let pref = gamma_ratio(&[c], &[a, b]);
        let (mut psi_a, mut psi_b, mut psi_1) =
            (digamma_numerical(a), digamma_numerical(b), euler);
        let (mut t, mut yn) = (1.0, 1.0);
        let mut sum = 0.0;
        for n in 0..5000 {
            let nf = f64::from(n);
            let term = t * yn * (2.0 * psi_1 - psi_a - psi_b - ly);
            sum += term;
            if n > 2 && term.abs() < 1e-18 * sum.abs() {
                break;
            }
            t *= (a + nf) * (b + nf) / ((nf + 1.0) * (nf + 1.0));
            yn *= y;
            psi_a += 1.0 / (a + nf);
            psi_b += 1.0 / (b + nf);
            psi_1 += 1.0 / (nf + 1.0);
        }
        return pref * sum;
    }
    let mf = f64::from(m);
    // finite part: sum_{n<m} (a)_n (b)_n / (n! (1-m)_n) y^n
    let pref1 = gamma_ratio(&[mf, c], &[a + mf, b + mf]);
    let mut fin = 0.0;
    let mut t = 1.0;
    for n in 0..m {
        fin += t;
        let nf = f64::from(n);
        t *= (a + nf) * (b + nf) / ((nf + 1.0) * (nf + 1.0 - mf)) * y;
    }
    // infinite part
    let pref2 = gamma_ratio(&[c], &[a, b]) * if m.is_multiple_of(2) { 1.0 } else { -1.0 };
    let mut psi_n1 = euler; // psi(n+1)
    let mut psi_nm1 = digamma_numerical(mf + 1.0); // psi(n+m+1)
    let mut psi_am = digamma_numerical(a + mf);
    let mut psi_bm = digamma_numerical(b + mf);
    let mut coef = 1.0 / gamma_numerical(mf + 1.0);
    let mut yp = y.powi(m as i32);
    let mut inf = 0.0;
    for n in 0..5000 {
        let nf = f64::from(n);
        let term = coef * yp * (ly - psi_n1 - psi_nm1 + psi_am + psi_bm);
        inf += term;
        if n > 2 && term.abs() < 1e-18 * (inf.abs() + fin.abs()) {
            break;
        }
        coef *= (a + mf + nf) * (b + mf + nf) / ((nf + 1.0) * (nf + mf + 1.0));
        yp *= y;
        psi_n1 += 1.0 / (nf + 1.0);
        psi_nm1 += 1.0 / (nf + mf + 1.0);
        psi_am += 1.0 / (a + mf + nf);
        psi_bm += 1.0 / (b + mf + nf);
    }
    pref1 * fin - pref2 * inf
}

/// `2F1` for `0 <= x < 1`: power series for `x <= 0.6`, otherwise the
/// `1 - x` connection formula (with the logarithmic limit when `c - a - b`
/// is an integer).
fn hyp2f1_unit(a: f64, b: f64, c: f64, x: f64) -> f64 {
    if x <= 0.6 || is_nonpos_int(a) || is_nonpos_int(b) {
        return hyp2f1_series(a, b, c, x);
    }
    let y = 1.0 - x;
    let s = c - a - b;
    if is_nonpos_int(c - a) || is_nonpos_int(c - b) {
        return y.powf(s) * hyp2f1_series(c - a, c - b, c, x);
    }
    if (s - s.round()).abs() < 1e-9 {
        let mi = s.round();
        if mi < 0.0 {
            // Euler: 2F1(a,b;c;x) = (1-x)^(c-a-b) 2F1(c-a,c-b;c;x)
            return y.powf(mi) * hyp2f1_unit(c - a, c - b, c, x);
        }
        return hyp2f1_log_case(a, b, mi as u32, y);
    }
    let lead = gamma_ratio(&[c, s], &[c - a, c - b]);
    let tail = gamma_ratio(&[c, -s], &[a, b]);
    lead * hyp2f1_series(a, b, 1.0 - s, y)
        + y.powf(s) * tail * hyp2f1_series(c - a, c - b, s + 1.0, y)
}

/// Gauss hypergeometric function `2F1(a, b; c; z)` for real arguments.
///
/// * `z = 0` and terminating series (`a` or `b` a non-positive integer)
///   are evaluated exactly for any `z`;
/// * `0 <= z <= 0.6`: power series;
/// * `0.6 < z < 1`: the `1 - z` connection formula, including the
///   logarithmic cases when `c - a - b` is an integer (digamma
///   expansions of Abramowitz-Stegun 15.3.10-15.3.11) and Euler's
///   transformation for negative integer `c - a - b`;
/// * `z < 0`: Pfaff's transformation `z -> z / (z - 1)` followed by the
///   above, which also realises the `1 / (1 - z)` and `1 / z`
///   continuations for large `|z|`;
/// * `z = 1`: Gauss's summation theorem (infinite when divergent);
/// * `z > 1`: the value is complex in general, so NaN, except for
///   polynomials (terminating series, or `c - a`, `c - b` a non-positive
///   integer with integer `c - a - b`).
///
/// Accuracy is about `1e-13` generically; when `c - a - b` is within
/// `1e-9` of, but not equal to, an integer the connection formula loses
/// digits to cancellation (roughly `1e-16 / |c - a - b - n|` relative).
/// A non-positive integer `c` (other than a terminating series) gives NaN.
#[must_use]
pub fn hyp2f1(a: f64, b: f64, c: f64, z: f64) -> f64 {
    if a.is_nan() || b.is_nan() || c.is_nan() || z.is_nan() {
        return f64::NAN;
    }
    if z == 0.0 {
        return 1.0;
    }
    if is_nonpos_int(c) {
        let ok = |p: f64| is_nonpos_int(p) && p.round() > c.round();
        if !(ok(a) || ok(b)) {
            return f64::NAN;
        }
    }
    if is_nonpos_int(a) || is_nonpos_int(b) {
        return hyp2f1_series(a, b, c, z);
    }
    let s = c - a - b;
    if z > 1.0 {
        if (is_nonpos_int(c - a) || is_nonpos_int(c - b)) && (s - s.round()).abs() < 1e-13 {
            return (1.0 - z).powi(s.round() as i32) * hyp2f1_series(c - a, c - b, c, z);
        }
        return f64::NAN;
    }
    if z == 1.0 {
        if s > 0.0 {
            return gamma_ratio(&[c, s], &[c - a, c - b]);
        }
        let sign = gamma_ratio(&[c], &[a, b]).signum();
        return sign * f64::INFINITY;
    }
    if z < 0.0 {
        let w = z / (z - 1.0);
        return (1.0 - z).powf(-a) * hyp2f1_unit(a, c - b, c, w);
    }
    hyp2f1_unit(a, b, c, z)
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

/// Starting order of the Miller backward recurrence for order `n` at
/// argument `x` (the amplitude is negligible beyond it).
fn miller_start(n: usize, x: f64) -> usize {
    let m = (n as f64).max(x);
    let start = (m + 20.0 * m.cbrt() + 10.0).ceil() as usize;
    start + (start % 2)
}

/// Power series of `J_n(x)` (`sign = -1`) or `I_n(x)` (`sign = +1`).
fn bessel_series(n: usize, x: f64, sign: f64) -> f64 {
    let h = 0.5 * x;
    let mut term = h.powi(n as i32);
    for k in 1..=n {
        term /= k as f64;
    }
    let mut sum = term;
    let q = sign * h * h;
    for k in 1..2000 {
        let kf = k as f64;
        term *= q / (kf * (kf + n as f64));
        sum += term;
        if term.abs() < 1e-18 * sum.abs() {
            break;
        }
    }
    sum
}

/// Hankel asymptotic terms `t_k = prod_{j<=k} (mu - (2j-1)^2) / (j 8 x)`
/// summed as `(P, Q)` (alternating-sign split by parity); `None` when
/// the series does not reach machine precision without a large
/// cancellation. For the modified functions `P` is the full alternating
/// sum and `Q` is unused.
fn hankel_sums(n: usize, x: f64, modified: bool) -> Option<(f64, f64)> {
    let mu = 4.0 * (n as f64) * (n as f64);
    let (mut p, mut q) = (1.0, 0.0);
    let mut t = 1.0_f64;
    let mut maxt = 1.0_f64;
    for k in 1..400_i32 {
        let kf = f64::from(k);
        let tn = t * (mu - (2.0 * kf - 1.0).powi(2)) / (kf * 8.0 * x);
        if tn.abs() > t.abs() && kf > 2.0 * x {
            return None;
        }
        maxt = maxt.max(tn.abs());
        t = tn;
        if modified {
            p += if k % 2 == 0 { t } else { -t };
        } else if k % 2 == 0 {
            p += if (k / 2) % 2 == 0 { t } else { -t };
        } else {
            q += if ((k - 1) / 2) % 2 == 0 { t } else { -t };
        }
        if t.abs() < 1e-17 {
            return if maxt < 50.0 { Some((p, q)) } else { None };
        }
    }
    None
}

/// Hankel expansion of `J_n(x)` for large `x > 0`.
fn bessel_j_hankel(n: usize, x: f64) -> Option<f64> {
    let (p, q) = hankel_sums(n, x, false)?;
    // chi = x - (n/2 + 1/4) pi, with the constant reduced exactly mod 2 pi
    let phi0 = ((2 * n + 1) % 8) as f64 * std::f64::consts::FRAC_PI_4;
    let (sx, cx) = x.sin_cos();
    let (s0, c0) = phi0.sin_cos();
    let cos_chi = cx * c0 + sx * s0;
    let sin_chi = sx * c0 - cx * s0;
    Some((2.0 / (PI * x)).sqrt() * (p * cos_chi - q * sin_chi))
}

/// Miller backward recurrence for `J_n(x)`, `x > 0`, normalised by
/// `J_0 + 2 sum J_{2k} = 1`.
fn bessel_j_miller(n: usize, x: f64) -> f64 {
    let top = miller_start(n, x);
    let (mut bjp, mut bj) = (0.0_f64, 1e-30_f64);
    let mut sum = 0.0;
    let mut ans = if top == n { bj } else { 0.0 };
    for k in (1..=top).rev() {
        let bjm = 2.0 * k as f64 / x * bj - bjp;
        bjp = bj;
        bj = bjm;
        // bj is now J_{k-1}
        if bj.abs() > 1e250 {
            bj *= 1e-250;
            bjp *= 1e-250;
            ans *= 1e-250;
            sum *= 1e-250;
        }
        let idx = k - 1;
        if idx == n {
            ans = bj;
        }
        if idx == 0 {
            sum += bj;
        } else if idx % 2 == 0 {
            sum += 2.0 * bj;
        }
    }
    ans / sum
}

/// Bessel function of the first kind `J_n(x)` for integer order.
///
/// Small arguments (`x^2/4 <= n + 1`) use the power series; large
/// arguments use the Hankel asymptotic expansion (with the phase
/// reduced exactly, so the accuracy does not degrade with `x`, cost
/// `O(1)`); everything else uses Miller's backward recurrence started
/// at an asymptotically chosen order (`~ max(n, x) + 20 max(n, x)^(1/3)`)
/// and normalised by `J_0 + 2 sum J_{2k} = 1`.
#[must_use]
pub fn bessel_jn(n: i32, x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    let nn = n.unsigned_abs() as usize;
    let mut sign = if n < 0 && nn % 2 == 1 { -1.0 } else { 1.0 };
    if x < 0.0 && nn % 2 == 1 {
        sign = -sign;
    }
    let ax = x.abs();
    if ax.is_infinite() {
        return 0.0;
    }
    if ax == 0.0 {
        return if nn == 0 { 1.0 } else { 0.0 };
    }
    let v = if ax * ax / 4.0 <= (nn + 1) as f64 {
        bessel_series(nn, ax, -1.0)
    } else if ax >= 25.0 {
        bessel_j_hankel(nn, ax).unwrap_or_else(|| bessel_j_miller(nn, ax))
    } else {
        bessel_j_miller(nn, ax)
    };
    sign * v
}

/// Miller recurrence for `e^{-x} I_n(x)`, `x > 0`, normalised by
/// `I_0 + 2 sum_{k>=1} I_k = e^x`.
fn bessel_i_miller_scaled(n: usize, x: f64) -> f64 {
    let top = miller_start(n, x);
    let (mut ip, mut i) = (0.0_f64, 1e-30_f64);
    let mut sum = 0.0;
    let mut ans = if top == n { i } else { 0.0 };
    for k in (1..=top).rev() {
        let im = ip + 2.0 * k as f64 / x * i;
        ip = i;
        i = im;
        if i.abs() > 1e250 {
            i *= 1e-250;
            ip *= 1e-250;
            ans *= 1e-250;
            sum *= 1e-250;
        }
        let idx = k - 1;
        if idx == n {
            ans = i;
        }
        sum += if idx == 0 { i } else { 2.0 * i };
    }
    ans / sum
}

/// Exponentially scaled modified Bessel function `e^{-|x|} I_n(x)` for
/// integer order; finite for every finite `x` (unlike [`bessel_in`], which
/// overflows beyond `|x| ~ 709`).
#[must_use]
pub fn bessel_in_scaled(n: i32, x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    let nn = n.unsigned_abs() as usize;
    let sign = if x < 0.0 && nn % 2 == 1 { -1.0 } else { 1.0 };
    let ax = x.abs();
    if ax.is_infinite() {
        return 0.0;
    }
    if ax == 0.0 {
        return if nn == 0 { 1.0 } else { 0.0 };
    }
    let v = if ax * ax / 4.0 <= (nn + 1) as f64 {
        bessel_series(nn, ax, 1.0) * (-ax).exp()
    } else if ax >= 30.0 {
        match hankel_sums(nn, ax, true) {
            Some((p, _)) => p / (2.0 * PI * ax).sqrt(),
            None => bessel_i_miller_scaled(nn, ax),
        }
    } else {
        bessel_i_miller_scaled(nn, ax)
    };
    sign * v
}

/// Modified Bessel function of the first kind `I_n(x)` for integer
/// order, by the same strategy as [`bessel_jn`] (series, Hankel
/// expansion, Miller recurrence). Overflows to infinity beyond
/// `|x| ~ 709`; use [`bessel_in_scaled`] there.
#[must_use]
pub fn bessel_in(n: i32, x: f64) -> f64 {
    bessel_in_scaled(n, x) * x.abs().exp()
}
