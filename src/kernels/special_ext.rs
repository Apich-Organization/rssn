//! # Additional special functions
//!
//! Complete and incomplete elliptic integrals (Carlson symmetric forms and
//! the arithmetic-geometric mean), Gauss hypergeometric `2F1`, Kummer
//! `1F1`, Bessel functions of integer order (series, Hankel asymptotics and
//! Miller backward recurrence) and Jacobi elliptic functions with exact
//! period reduction.
//!
//! Elliptic integrals use the *parameter* `m = k^2`.
//!
//! The module also carries:
//!
//! * the Coulomb wave functions `F_l(η, ρ)`, `G_l(η, ρ)` and their
//!   derivatives, with the normalising constant `C_l(η)` and phase shift
//!   `σ_l(η)` ([`coulomb_f`], [`coulomb_g`], [`coulomb_f_prime`],
//!   [`coulomb_g_prime`], [`coulomb_c`], [`coulomb_sigma`]); `F` from the
//!   confluent series, `G` from `H⁺ = G + iF` as a Tricomi function
//!   ([`coulomb_g_from_u`], about `1e-13` against mpmath over `ρ` from 0.05 to
//!   15, `η` from -3 to 4, `l` to 10; the Steed continued fraction used before
//!   was off by up to `1e-6` and worse at small `ρ`) and the asymptotic
//!   expansion beyond `ρ = 24`;
//! * the polylogarithm `Li_s(z)` for real or complex `s` and `z`
//!   ([`polylog`], [`polylog_complex`]) and the Lerch transcendent
//!   `Φ(z, s, a)` ([`lerch`]);
//! * Appell's double hypergeometric function `F_1` ([`appell_f1`]);
//! * the Kelvin functions `ber_ν`, `bei_ν`, `ker_ν`, `kei_ν`
//!   ([`kelvin_ber`], [`kelvin_bei`], [`kelvin_ker`], [`kelvin_kei`]);
//! * the Anger and Weber functions `J_ν`, `E_ν` ([`anger_j`],
//!   [`weber_e`]).
//!
//! Complex parameters and `pFq`, Kummer `M` and Tricomi `U` live in
//! [`crate::kernels::hypergeometric`].

use std::f64::consts::PI;

use num_complex::Complex64;

use crate::kernels::dd::Cdd;
use crate::kernels::dd::Dd;
use crate::kernels::special::bernoulli_number;
use crate::kernels::special::digamma_numerical;
use crate::kernels::special::gamma_numerical;
use crate::kernels::special::ln_gamma_numerical;

/// Carlson's symmetric elliptic integral `R_F(x, y, z)` (arguments non-negative, at most one zero).
#[must_use]
pub fn carlson_rf(
    x: f64,
    y: f64,
    z: f64,
) -> f64 {
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
pub fn carlson_rd(
    x: f64,
    y: f64,
    z: f64,
) -> f64 {
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
            * (1.0
                + ed * (-c1 + c5 * ed - c6 * delz * ee)
                + delz * (c2 * ee + delz * (-c3 * ec + delz * c4 * ea)))
            / (ave * ave.sqrt())
}

/// Complete elliptic integral of the first kind `K(m)`, `m < 1`.
#[must_use]
pub fn elliptic_k(m: f64) -> f64 {
    if crate::kernels::exact_eq(m, 1.0) {
        return f64::INFINITY;
    }
    if m > 1.0 || m.is_nan() {
        return f64::NAN;
    }
    // AGM: K = pi / (2 agm(1, sqrt(1-m)))
    let (mut a, mut b) = (1.0, (1.0 - m).sqrt());
    for _ in 0..64 {
        let an = f64::midpoint(a, b);
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
    if crate::kernels::exact_eq(m, 1.0) {
        return 1.0;
    }
    if m > 1.0 || m.is_nan() {
        return f64::NAN;
    }
    carlson_rf(0.0, 1.0 - m, 1.0) - m / 3.0 * carlson_rd(0.0, 1.0 - m, 1.0)
}

/// Incomplete elliptic integral of the first kind `F(phi, m)`.
#[must_use]
pub fn elliptic_f(
    phi: f64,
    m: f64,
) -> f64 {
    // reduce to |phi| <= pi/2 using F(phi + n pi) = F(phi) + 2 n K
    let n = (phi / PI).round();
    let r = phi - n * PI;
    let (s, c) = r.sin_cos();
    let v = s * carlson_rf(c * c, 1.0 - m * s * s, 1.0);
    v + 2.0 * n * elliptic_k(m)
}

/// Incomplete elliptic integral of the second kind `E(phi, m)`.
#[must_use]
pub fn elliptic_e_inc(
    phi: f64,
    m: f64,
) -> f64 {
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
pub fn jacobi_elliptic(
    u: f64,
    m: f64,
) -> (f64, f64, f64) {
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
    let (v, flip) = if a > kk {
        (2.0 * kk - a, true)
    } else {
        (a, false)
    };
    let (s, c) = jacobi_landen(v, m);
    let sn = if r < 0.0 { -s } else { s };
    let cn = if flip { -c } else { c };
    let dn = ((1.0 - m) + m * c * c).sqrt();
    (sn, cn, dn)
}

/// `(sin am, cos am)` for `u` in `[0, K]` by descending Landen iteration.
fn jacobi_landen(
    u: f64,
    m: f64,
) -> (f64, f64) {
    let mut a = [0.0; 40];
    let mut c = [0.0; 40];
    a[0] = 1.0;
    let mut b = (1.0 - m).sqrt();
    c[0] = m.sqrt();
    let mut n = 0;
    while n < 39 && c[n].abs() > 1e-16 {
        a[n + 1] = f64::midpoint(a[n], b);
        c[n + 1] = 0.5 * (a[n] - b);
        b = (a[n] * b).sqrt();
        n += 1;
    }
    let mut phi = (2.0_f64).powi(n as i32) * a[n] * u;
    let mut i = n;
    while i > 0 {
        let arg = (c[i] * phi.sin() / a[i]).clamp(-1.0, 1.0);
        phi = f64::midpoint(phi, arg.asin());
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
        (
            (PI / s.abs()).ln() - ln_gamma_numerical(1.0 - x),
            s.signum(),
        )
    }
}

/// `prod Gamma(num) / prod Gamma(den)` evaluated through logarithms
/// (poles in the denominator give 0, poles in the numerator infinity).
fn gamma_ratio(
    num: &[f64],
    den: &[f64],
) -> f64 {
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
fn hyp2f1_series(
    a: f64,
    b: f64,
    c: f64,
    z: f64,
) -> f64 {
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
fn hyp2f1_log_case(
    a: f64,
    b: f64,
    m: u32,
    y: f64,
) -> f64 {
    let c = a + b + f64::from(m);
    let ly = y.ln();
    let euler = -std::f64::consts::EULER_GAMMA;
    if m == 0 {
        let pref = gamma_ratio(&[c], &[a, b]);
        let (mut psi_a, mut psi_b, mut psi_1) = (digamma_numerical(a), digamma_numerical(b), euler);
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
    let pref2 = gamma_ratio(&[c], &[a, b])
        * if m.is_multiple_of(2) {
            1.0
        } else {
            -1.0
        };
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
fn hyp2f1_unit(
    a: f64,
    b: f64,
    c: f64,
    x: f64,
) -> f64 {
    if x <= 0.6 || is_nonpos_int(a) || is_nonpos_int(b) {
        return hyp2f1_series(a, b, c, x);
    }
    let y = 1.0 - x;
    let s = c - a - b;
    if is_nonpos_int(c - a) || is_nonpos_int(c - b) {
        return y.powf(s) * hyp2f1_series(c - a, c - b, c, x);
    }
    let eps = (s - s.round()).abs();
    if eps > 0.0 && eps < 1e-4 {
        // The connection formula cancels catastrophically here; continue the
        // ODE solution analytically instead.
        let v = hyp2f1_ode(a, b, c, Complex64::new(x, 0.0));
        if v.re.is_finite() {
            return v.re;
        }
    }
    if eps == 0.0 {
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
/// `1e-4` of, but not equal to, an integer the evaluation switches to the
/// analytic continuation of the differential equation, which has no
/// cancellation (exact integers use the logarithmic limit formulas).
/// A non-positive integer `c` (other than a terminating series) gives NaN.
/// For `z > 1` use [`hyp2f1_complex`].
#[must_use]
pub fn hyp2f1(
    a: f64,
    b: f64,
    c: f64,
    z: f64,
) -> f64 {
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
    if crate::kernels::exact_eq(z, 1.0) {
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

/// Power series of `2F1` in a complex argument (`|z| < 1`).
fn hyp2f1_series_c(
    a: f64,
    b: f64,
    c: f64,
    z: Complex64,
) -> Complex64 {
    let mut term = Complex64::new(1.0, 0.0);
    let mut sum = term;
    for n in 0..2_000_000 {
        let nf = f64::from(n);
        let ratio = (a + nf) * (b + nf) / ((c + nf) * (nf + 1.0));
        term *= z * ratio;
        sum += term;
        if term.norm() == 0.0
            || (z.norm() < 1.0 && term.norm() < 1e-17 * sum.norm() && ratio.abs() * z.norm() < 1.0)
        {
            break;
        }
    }
    sum
}

/// Analytic continuation of `2F1(a, b; c; z)` by Taylor stepping of the
/// hypergeometric differential equation
/// `z (1 - z) w'' + [c - (a + b + 1) z] w' - a b w = 0` from `z = 1/2`
/// along a path in the plane cut along `[1, inf)` (for real `z > 1` the
/// path stays in the upper half plane, i.e. the value is the limit from
/// `Im z -> 0+`). Every step has a radius of at most 0.4 times the distance
/// to the nearest singular point (0 or 1), so the method is uniformly
/// accurate for `|z|` large, `z` near 1 and `c - a - b` close to an integer.
fn hyp2f1_ode(
    a: f64,
    b: f64,
    c: f64,
    z: Complex64,
) -> Complex64 {
    if z.norm() <= 0.5 {
        return hyp2f1_series_c(a, b, c, z);
    }
    let half = Complex64::new(0.5, 0.0);
    let f0 = hyp2f1_series_c(a, b, c, half);
    let d0 = hyp2f1_series_c(a + 1.0, b + 1.0, c + 1.0, half) * (a * b / c);
    let mut p = half;
    let (mut fv, mut dv) = (f0, d0);
    let side = if z.im < 0.0 { -1.0 } else { 1.0 };
    let waypoints = [Complex64::new(0.5, 0.5 * side), z];
    let ab = a * b;
    let apb1 = a + b + 1.0;
    for &target in &waypoints {
        let mut guard = 0;
        while (target - p).norm() > 0.0 && guard < 100_000 {
            guard += 1;
            let rho = p.norm().min((Complex64::new(1.0, 0.0) - p).norm());
            let rem = target - p;
            let len = rem.norm().min(0.4 * rho);
            let dir = rem / rem.norm();
            let xi = len / rho;
            // scaled Taylor coefficients v_k = u_k rho^k
            let p0 = p * (Complex64::new(1.0, 0.0) - p);
            let p1 = Complex64::new(1.0, 0.0) - p * 2.0;
            let q0 = Complex64::new(c, 0.0) - p * apb1;
            let mut v = vec![fv, dv * rho];
            let kmax = 80;
            let mut sum_f = Complex64::new(0.0, 0.0);
            let mut sum_d = Complex64::new(0.0, 0.0);
            let w = dir * xi; // x / rho = w t; use powers of w
            for k in 0..kmax {
                let kf = k as f64;
                let num = (p1 * kf + q0) * ((kf + 1.0) * rho) * v[k + 1]
                    - Complex64::new(kf * (kf - 1.0) + apb1 * kf + ab, 0.0) * (rho * rho) * v[k];
                v.push(-num / (p0 * ((kf + 2.0) * (kf + 1.0))));
            }
            let mut wp = Complex64::new(1.0, 0.0);
            for k in 0..v.len() {
                sum_f += v[k] * wp;
                if k + 1 < v.len() {
                    sum_d += v[k + 1] * wp * (k as f64 + 1.0);
                }
                wp *= w;
            }
            fv = sum_f;
            dv = sum_d / rho;
            p += dir * len;
            if (target - p).norm() < 1e-15 * (1.0 + target.norm()) {
                p = target;
                break;
            }
        }
    }
    fv
}

/// Gauss hypergeometric function `2F1(a, b; c; z)` for a complex argument
/// and real parameters, on the principal branch (cut along `[1, inf)`;
/// real `z > 1` gives the limit from the upper half plane).
///
/// Uses the series for `|z| <= 1/2` and analytic continuation of the
/// differential equation elsewhere, which covers `|z| > 1`, points near
/// `z = 1` and `c - a - b` near an integer without cancellation. Terminating
/// series are summed exactly. `z = 1` follows Gauss's theorem for real
/// parameters (infinite when divergent); a non-positive integer `c` without
/// termination gives NaN.
#[must_use]
pub fn hyp2f1_complex(
    a: f64,
    b: f64,
    c: f64,
    z: Complex64,
) -> Complex64 {
    let nan = Complex64::new(f64::NAN, f64::NAN);
    if a.is_nan() || b.is_nan() || c.is_nan() || z.re.is_nan() || z.im.is_nan() {
        return nan;
    }
    if z.norm() == 0.0 {
        return Complex64::new(1.0, 0.0);
    }
    if is_nonpos_int(a) || is_nonpos_int(b) {
        let ok = |p: f64| is_nonpos_int(p) && (!is_nonpos_int(c) || p.round() > c.round());
        if !(ok(a) || ok(b)) {
            return nan;
        }
        // polynomial
        let mut term = Complex64::new(1.0, 0.0);
        let mut sum = term;
        for n in 0..100_000 {
            let nf = f64::from(n);
            term *= z * ((a + nf) * (b + nf) / ((c + nf) * (nf + 1.0)));
            if term.norm() == 0.0 {
                break;
            }
            sum += term;
        }
        return sum;
    }
    if is_nonpos_int(c) {
        return nan;
    }
    if crate::kernels::exact_eq(z.im, 0.0) && crate::kernels::exact_eq(z.re, 1.0) {
        return Complex64::new(hyp2f1(a, b, c, 1.0), 0.0);
    }
    hyp2f1_ode(a, b, c, z)
}

/// Kummer confluent hypergeometric function `1F1(a; b; x)`.
#[must_use]
pub fn hyp1f1(
    a: f64,
    b: f64,
    x: f64,
) -> f64 {
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
fn miller_start(
    n: usize,
    x: f64,
) -> usize {
    let m = (n as f64).max(x);
    let start = (m + 20.0 * m.cbrt() + 10.0).ceil() as usize;
    start + (start % 2)
}

/// Power series of `J_n(x)` (`sign = -1`) or `I_n(x)` (`sign = +1`).
fn bessel_series(
    n: usize,
    x: f64,
    sign: f64,
) -> f64 {
    let h = 0.5 * x;
    let mut term = h.powi(n as i32);
    for k in 1..=n {
        term /= k as f64;
    }
    let mut sum = term;
    let q = sign * h * h;
    for k in 1..2000 {
        let kf = f64::from(k);
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
fn hankel_sums(
    n: usize,
    x: f64,
    modified: bool,
) -> Option<(f64, f64)> {
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
            q += if ((k - 1) / 2) % 2 == 0 {
                t
            } else {
                -t
            };
        }
        if t.abs() < 1e-17 {
            return if maxt < 50.0 {
                Some((p, q))
            } else {
                None
            };
        }
    }
    None
}

/// Hankel expansion of `J_n(x)` for large `x > 0`.
fn bessel_j_hankel(
    n: usize,
    x: f64,
) -> Option<f64> {
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
fn bessel_j_miller(
    n: usize,
    x: f64,
) -> f64 {
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
pub fn bessel_jn(
    n: i32,
    x: f64,
) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    let nn = n.unsigned_abs() as usize;
    let mut sign = if n < 0 && nn % 2 == 1 {
        -1.0
    } else {
        1.0
    };
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
fn bessel_i_miller_scaled(
    n: usize,
    x: f64,
) -> f64 {
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
pub fn bessel_in_scaled(
    n: i32,
    x: f64,
) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    let nn = n.unsigned_abs() as usize;
    let sign = if x < 0.0 && nn % 2 == 1 {
        -1.0
    } else {
        1.0
    };
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
            | Some((p, _)) => p / (2.0 * PI * ax).sqrt(),
            | None => bessel_i_miller_scaled(nn, ax),
        }
    } else {
        bessel_i_miller_scaled(nn, ax)
    };
    sign * v
}

/// Modified Bessel function of the first kind `I_n(x)` for integer order.
///
/// Modified Bessel function of the first kind `I_n(x)` for integer
/// order, by the same strategy as [`bessel_jn`] (series, Hankel
/// expansion, Miller recurrence). Overflows to infinity beyond
/// `|x| ~ 709`; use [`bessel_in_scaled`] there.
#[must_use]
pub fn bessel_in(
    n: i32,
    x: f64,
) -> f64 {
    bessel_in_scaled(n, x) * x.abs().exp()
}

// ============================================================================
// Coulomb wave functions
// ============================================================================

/// Bernoulli numbers `B_2, B_4, …`: the Stirling-series coefficients of
/// [`complex_ln_gamma`].
const STIRLING: [f64; 8] = [
    1.0 / 6.0,
    -1.0 / 30.0,
    1.0 / 42.0,
    -1.0 / 30.0,
    5.0 / 66.0,
    -691.0 / 2730.0,
    7.0 / 6.0,
    -3617.0 / 510.0,
];

/// Above this `ρ` the Coulomb functions use the large-`ρ` asymptotic
/// expansion, provided the classically allowed condition
/// `ρ ≥ 2|η| + 2` also holds.
const COULOMB_ASYMPTOTIC_MIN: f64 = 24.0;
/// The largest `ρ` at which the confluent series is used in the forbidden
/// region; beyond it neither representation is valid and the kernels return
/// `None`.
const COULOMB_SERIES_MAX: f64 = 34.0;

/// `(ln Γ(z), arg Γ(z))` for complex `z`, by Stirling's series after
/// shifting the argument to `Re z ≥ 16` (relative accuracy about `1e-15`).
///
/// Also used by [`polylog_complex`] through the complex Riemann zeta.
#[must_use]
pub(crate) fn complex_ln_gamma(z: Complex64) -> Complex64 {
    let mut z = z;
    let mut shift = Complex64::new(0.0, 0.0);
    while z.re < 16.0 {
        shift -= z.ln();
        z += Complex64::new(1.0, 0.0);
    }
    let mut tail = Complex64::new(0.0, 0.0);
    let mut power = z;
    let square = z * z;
    for (n, &coefficient) in STIRLING.iter().enumerate() {
        let k = (n + 1) as f64;
        tail += Complex64::new(coefficient / (2.0 * k * (2.0 * k - 1.0)), 0.0) / power;
        power *= square;
    }
    shift + (z - Complex64::new(0.5, 0.0)) * z.ln() - z
        + Complex64::new(0.5 * (2.0 * PI).ln(), 0.0)
        + tail
}

/// The Coulomb normalising constant `C_l(η)` (DLMF 33.2.5),
/// `C_l(η) = 2^l e^{-πη/2} |Γ(l + 1 + iη)| / Γ(2l + 2)`, evaluated in
/// logarithmic form so that large `|η|` neither overflows nor underflows
/// prematurely.
#[must_use]
pub fn coulomb_c(
    eta: f64,
    l: usize,
) -> f64 {
    if eta.is_nan() {
        return f64::NAN;
    }
    let lf = l as f64;
    let ln_abs = complex_ln_gamma(Complex64::new(lf + 1.0, eta)).re;
    let ln_den = complex_ln_gamma(Complex64::new(2.0 * lf + 2.0, 0.0)).re;
    let ln_c = lf * std::f64::consts::LN_2 - 0.5 * PI * eta + ln_abs - ln_den;
    if ln_c < -745.0 {
        0.0
    } else {
        ln_c.exp()
    }
}

/// The Coulomb phase shift `σ_l(η) = arg Γ(l + 1 + iη)` (DLMF 33.2.10).
#[must_use]
pub fn coulomb_sigma(
    eta: f64,
    l: usize,
) -> f64 {
    if eta.is_nan() {
        return f64::NAN;
    }
    complex_ln_gamma(Complex64::new(l as f64 + 1.0, eta)).im
}

/// The regular Coulomb function `F_l(η, ρ)` and its derivative from the
/// confluent-hypergeometric series (DLMF 33.2.4),
/// `F = C_l ρ^{l+1} Re[e^{-iρ} M(l + 1 - iη; 2l + 2; 2iρ)]`.
///
/// Below `ρ = 3` the sum is carried in `f64`; above it double-double
/// arithmetic is used, because the terms reach `e^{2ρ}` while the result is
/// of order one, so ordinary `f64` would lose about `2ρ / ln 10` digits to
/// cancellation.
#[allow(clippy::single_call_fn, clippy::float_cmp)]
pub(crate) fn coulomb_series(
    eta: f64,
    l: usize,
    rho: f64,
) -> (f64, f64) {
    let lf = l as f64;
    let c = coulomb_c(eta, l);
    if c == 0.0 || rho == 0.0 {
        return (0.0, 0.0);
    }
    if rho <= 3.0 {
        return coulomb_series_f64(eta, l, rho, c);
    }
    let a = Cdd::new(lf + 1.0, -eta);
    let b = Cdd::new(2.0 * lf + 2.0, 0.0);
    let z = Cdd::new(0.0, 2.0 * rho);
    let mut term = Cdd::new(1.0, 0.0);
    let mut sum = term;
    let mut moment = Cdd::new(0.0, 0.0);
    let mut k = 0_usize;
    loop {
        moment += Cdd::new(k as f64, 0.0) * term;
        k += 1;
        let kf = k as f64;
        term = term * (a + Cdd::new(kf - 1.0, 0.0)) / (b + Cdd::new(kf - 1.0, 0.0)) * z
            / Cdd::new(kf, 0.0);
        sum += term;
        if k > 4 && (term.norm() < 1e-30 * sum.norm() || k > 400) {
            break;
        }
    }
    let phase = Cdd::new(rho.cos(), -rho.sin());
    let rotated = phase * sum;
    let slope = phase * moment;
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let pref = c * rho.powi(i32::try_from(l).unwrap_or(0) + 1);
    let f = pref * rotated.re.to_f64();
    let derivative = (Dd::from_f64((lf + 1.0) / rho) * rotated.re
        + rotated.im
        + Dd::from_f64(1.0 / rho) * slope.re)
        .to_f64();
    (f, pref * derivative)
}

/// The same confluent series in ordinary `f64`, used where the cancellation
/// is mild (`ρ ≤ 3`).
#[allow(clippy::single_call_fn)]
fn coulomb_series_f64(
    eta: f64,
    l: usize,
    rho: f64,
    c: f64,
) -> (f64, f64) {
    let lf = l as f64;
    let (a, b) = (
        Complex64::new(lf + 1.0, -eta),
        Complex64::new(2.0 * lf + 2.0, 0.0),
    );
    let z = Complex64::new(0.0, 2.0 * rho);
    let mut term = Complex64::new(1.0, 0.0);
    let mut sum = term;
    let mut moment = Complex64::new(0.0, 0.0);
    let mut k = 0_usize;
    loop {
        moment += Complex64::new(k as f64, 0.0) * term;
        k += 1;
        let kf = k as f64;
        term *= (a + Complex64::new(kf - 1.0, 0.0)) / (b + Complex64::new(kf - 1.0, 0.0)) * z
            / Complex64::new(kf, 0.0);
        sum += term;
        if k > 4 && (term.norm() < 1e-30 * sum.norm() || k > 400) {
            break;
        }
    }
    let phase = Complex64::new(rho.cos(), -rho.sin());
    let rotated = phase * sum;
    let slope = phase * moment;
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let pref = c * rho.powi(i32::try_from(l).unwrap_or(0) + 1);
    let f = pref * rotated.re;
    let fp = pref * ((lf + 1.0) / rho * rotated.re + rotated.im + slope.re / rho);
    (f, fp)
}

/// `F`, `F'`, `G`, `G'` from the large-`ρ` asymptotic expansion
/// (DLMF 33.11.2-33.11.8) truncated at its smallest term. The phase is
/// `θ = ρ - η ln(2ρ) - lπ/2 + σ_l`, so `F ~ sin θ` and `G ~ cos θ`.
///
/// `None` when the expansion has not reached its smallest term before the
/// `f64` cancellation budget is exhausted.
#[allow(clippy::single_call_fn)]
pub(crate) fn coulomb_asymptotic(
    eta: f64,
    l: usize,
    rho: f64,
) -> Option<(f64, f64, f64, f64)> {
    let lf = l as f64;
    let (mut f, mut g) = (1.0, 0.0);
    let (mut f_hat, mut g_hat) = (0.0, 1.0 - eta / rho);
    let (mut sum_f, mut sum_g) = (1.0, 0.0);
    let (mut sum_f_hat, mut sum_g_hat) = (0.0, 1.0 - eta / rho);
    let mut smallest = f64::INFINITY;
    for k in 0..400_usize {
        let kf = k as f64;
        let lambda = eta * (2.0 * kf + 1.0) / (2.0 * rho * (kf + 1.0));
        let mu = (lf * (lf + 1.0) + eta * eta - kf * (kf + 1.0)) / (2.0 * rho * (kf + 1.0));
        let next_f = lambda * f - mu * g;
        let next_g = lambda * g + mu * f;
        let next_f_hat = lambda * f_hat - mu * g_hat - next_f / rho;
        let next_g_hat = lambda * g_hat + mu * f_hat - next_g / rho;
        let size = next_f
            .abs()
            .max(next_g.abs())
            .max(next_f_hat.abs())
            .max(next_g_hat.abs());
        if k > 4 && size > smallest {
            break;
        }
        sum_f += next_f;
        sum_g += next_g;
        sum_f_hat += next_f_hat;
        sum_g_hat += next_g_hat;
        smallest = smallest.min(size);
        f = next_f;
        g = next_g;
        f_hat = next_f_hat;
        g_hat = next_g_hat;
        if size < 1e-18 {
            break;
        }
    }
    let theta = rho - eta * (2.0 * rho).ln() - lf * 0.5 * PI + coulomb_sigma(eta, l);
    let (sine, cosine) = theta.sin_cos();
    let regular = sum_f * sine + sum_g * cosine;
    let irregular = sum_f * cosine - sum_g * sine;
    let regular_prime = sum_g_hat * cosine + sum_f_hat * sine;
    let irregular_prime = sum_f_hat * cosine - sum_g_hat * sine;
    (regular.is_finite()
        && irregular.is_finite()
        && regular_prime.is_finite()
        && irregular_prime.is_finite())
    .then_some((regular, regular_prime, irregular, irregular_prime))
}

/// Steed's continued fractions (DLMF 33.8.1, 33.8.2): the logarithmic
/// derivative `u = F'_l / F_l` and `p + iq = H⁺'_l / H⁺_l` with
/// `H⁺ = G + iF`. The caller checks `q > 0`, which certifies that the
/// fractions converged onto the physical solution.
#[allow(clippy::single_call_fn)]
fn steed_ratios(
    eta: f64,
    l: usize,
    rho: f64,
) -> Option<(f64, f64, f64)> {
    let lf = l as f64;
    let radius = |j: f64| (1.0 + eta * eta / (j * j)).sqrt();
    let shift = |j: f64| j / rho + eta / j;
    let term = |j: f64| shift(j) + shift(j + 1.0);
    let depth = (60.0 + 2.0 * rho + 2.0 * eta.abs() + 4.0 * lf).ceil() as usize;
    let start = lf + 1.0 + depth as f64;
    let mut denominator = term(start);
    let mut j = start;
    while j > lf + 1.0 + 1e-9 {
        j -= 1.0;
        let r = radius(j + 1.0);
        denominator = term(j) - r * r / denominator;
    }
    let u = shift(lf + 1.0) - radius(lf + 1.0).powi(2) / denominator;
    // H⁺'/H⁺ = c + (i/ρ) a b / C with C = D₀ + (a+1)(b+1)/(D₁ + ...),
    // D_k = 2(ρ - η + i(k+1)), evaluated by the modified Lentz algorithm.
    let a = Complex64::new(1.0 + lf, eta);
    let b = Complex64::new(-lf, eta);
    let dk = |k: usize| Complex64::new(2.0 * (rho - eta), 2.0 * (k as f64 + 1.0));
    let mut f = dk(0);
    let mut c = dk(0);
    let mut d = Complex64::new(0.0, 0.0);
    for k in 1..10_000_usize {
        let kf = k as f64;
        let numerator = (a + Complex64::new(kf, 0.0)) * (b + Complex64::new(kf, 0.0));
        let base = dk(k);
        d = base + numerator * d;
        if d.norm() < 1e-300 {
            d = Complex64::new(1e-300, 0.0);
        }
        d = Complex64::new(1.0, 0.0) / d;
        c = base + numerator / c;
        if c.norm() < 1e-300 {
            c = Complex64::new(1e-300, 0.0);
        }
        let delta = c * d;
        f *= delta;
        if (delta - Complex64::new(1.0, 0.0)).norm() < 1e-16 {
            break;
        }
    }
    let quotient = a * b / f;
    let p_iq = Complex64::new(0.0, 1.0 - eta / rho) + Complex64::new(0.0, 1.0 / rho) * quotient;
    let values = (u, p_iq.re, p_iq.im);
    (u.is_finite() && p_iq.re.is_finite() && p_iq.im.is_finite()).then_some(values)
}

/// The irregular Coulomb function `G_l(η, ρ)` and its derivative from the Tricomi function.
///
/// (DLMF 33.2.7): `H⁺ = G + iF = e^{iθ} (-2iρ)^{l+1+iη} U(l + 1 + iη, 2l + 2, -2iρ)` with `θ = ρ - η ln
/// 2ρ - lπ/2 + σ_l`, and `dU/dρ = 2i a U(a + 1, b + 1, -2iρ)` for the derivative.
///
/// Unlike Steed's continued fraction for `H⁺'/H⁺`, which converges slowly and loses digits at small
/// `ρ`, this needs no certification; its accuracy is that of
/// [`crate::kernels::hypergeometric::hyp_u_cplx`] on the imaginary axis. `None` for `ρ <= 0`, `l >
/// 64` or non-finite input.
#[must_use]
pub fn coulomb_g_from_u(
    eta: f64,
    l: usize,
    rho: f64,
) -> Option<(f64, f64)> {
    use crate::kernels::hypergeometric::hyp_u_cplx;
    if !eta.is_finite() || !rho.is_finite() || rho <= 0.0 || l > 64 {
        return None;
    }
    let lf = l as f64;
    let a = Complex64::new(lf + 1.0, eta);
    let b = Complex64::new(2.0 * lf + 2.0, 0.0);
    let z = Complex64::new(0.0, -2.0 * rho);
    let theta = rho - eta * (2.0 * rho).ln() - lf * 0.5 * PI + coulomb_sigma(eta, l);
    let ln_z = Complex64::new((2.0 * rho).ln(), -0.5 * PI);
    let prefactor = (Complex64::new(0.0, theta) + a * ln_z).exp();
    let u = hyp_u_cplx(a, b, z);
    let u_next = hyp_u_cplx(a + 1.0, b + 1.0, z);
    let h = prefactor * u;
    let hp = h * (Complex64::new(0.0, 1.0 - eta / rho) + a / rho)
        + prefactor * Complex64::new(0.0, 2.0) * a * u_next;
    (h.re.is_finite() && hp.re.is_finite()).then_some((h.re, hp.re))
}

/// The regular and irregular Coulomb functions and their derivatives.
///
/// `F_l(η, ρ)` is the solution regular at the origin and `G_l(η, ρ)` the
/// irregular solution normalised by the Wronskian `F' G - F G' = 1` (so
/// that `G ~ cos θ` at large `ρ`). The dimensionless variables are the
/// Sommerfeld parameter `η` and `ρ = k r`; `l ≥ 0` is the integer angular
/// momentum. The result is `None` outside the domain of the implemented
/// representations (`ρ ≤ 0`, non-finite arguments, `l > 64`, or
/// `24 < ρ ≤ 34` outside the classically allowed region), and also when
/// Steed's fraction fails to certify the physical solution.
///
/// References: DLMF §§33.2, 33.8 and 33.11; the series is evaluated in
/// double-double arithmetic beyond `ρ = 3` to defeat the `e^{2ρ}`
/// cancellation.
#[must_use]
#[allow(clippy::float_cmp)] // exact comparison against a sentinel / integer-valued input is intended
pub fn coulomb_fg(
    eta: f64,
    l: usize,
    rho: f64,
) -> Option<(f64, f64, f64, f64)> {
    if !eta.is_finite() || !rho.is_finite() || rho <= 0.0 || l > 64 {
        return None;
    }
    if rho >= COULOMB_ASYMPTOTIC_MIN && rho >= 2.0 * eta.abs() + 2.0 {
        return coulomb_asymptotic(eta, l, rho);
    }
    if rho > COULOMB_SERIES_MAX {
        return None;
    }
    let (f, fp) = coulomb_series(eta, l, rho);
    // G and G' from the Tricomi function (accurate to ~1e-13 on the whole
    // range); Steed's continued fraction, which loses digits at small rho
    // (relative errors of 1e-6 and worse for rho < 1), is the fallback.
    if let Some((g, gp)) = coulomb_g_from_u(eta, l, rho) {
        let values = (f, fp, g, gp);
        if f.is_finite() && fp.is_finite() {
            return Some(values);
        }
    }
    let (u, p, q) = steed_ratios(eta, l, rho)?;
    if !q.is_finite() || q <= 0.0 {
        return None;
    }
    let g = (u - p) * f / q;
    let gp = (p * (u - p) - q * q) * f / q;
    let values = (f, fp, g, gp);
    (f.is_finite() && fp.is_finite() && g.is_finite() && gp.is_finite()).then_some(values)
}

/// The regular Coulomb function `F_l(η, ρ)` (DLMF 33.2.4).
///
/// Domain and accuracy: see [`coulomb_fg`]; `None` where no representation
/// is available.
#[must_use]
pub fn coulomb_f(
    eta: f64,
    l: usize,
    rho: f64,
) -> Option<f64> {
    if rho == 0.0 && eta.is_finite() && l <= 64 {
        return Some(0.0); // F_l(eta, 0) = 0
    }
    if !eta.is_finite() || !rho.is_finite() || rho <= 0.0 || l > 64 {
        return None;
    }
    if rho >= COULOMB_ASYMPTOTIC_MIN && rho >= 2.0 * eta.abs() + 2.0 {
        return coulomb_asymptotic(eta, l, rho).map(|values| values.0);
    }
    if rho > COULOMB_SERIES_MAX {
        return None;
    }
    Some(coulomb_series(eta, l, rho).0)
}

/// The irregular Coulomb function `G_l(η, ρ)`, normalised by
/// `F' G - F G' = 1` (DLMF 33.2.6).
///
/// Domain and accuracy: see [`coulomb_fg`].
#[must_use]
pub fn coulomb_g(
    eta: f64,
    l: usize,
    rho: f64,
) -> Option<f64> {
    coulomb_fg(eta, l, rho).map(|values| values.2)
}

/// The derivative `dF_l/dρ` of the regular Coulomb function.
///
/// Domain and accuracy: see [`coulomb_fg`].
#[must_use]
pub fn coulomb_f_prime(
    eta: f64,
    l: usize,
    rho: f64,
) -> Option<f64> {
    if rho == 0.0 && eta.is_finite() && l <= 64 {
        // F'_0(eta, 0) = C_0(eta), F'_l(eta, 0) = 0 for l >= 1
        return Some(if l == 0 {
            coulomb_c(eta, 0)
        } else {
            0.0
        });
    }
    if !eta.is_finite() || !rho.is_finite() || rho <= 0.0 || l > 64 {
        return None;
    }
    if rho >= COULOMB_ASYMPTOTIC_MIN && rho >= 2.0 * eta.abs() + 2.0 {
        return coulomb_asymptotic(eta, l, rho).map(|values| values.1);
    }
    if rho > COULOMB_SERIES_MAX {
        return None;
    }
    Some(coulomb_series(eta, l, rho).1)
}

/// The derivative `dG_l/dρ` of the irregular Coulomb function.
///
/// Domain and accuracy: see [`coulomb_fg`].
#[must_use]
pub fn coulomb_g_prime(
    eta: f64,
    l: usize,
    rho: f64,
) -> Option<f64> {
    coulomb_fg(eta, l, rho).map(|values| values.3)
}

/// The energy-normalised scattering solution
/// `u_k(r) = √(2/π) F_l(η, k r)`, `η = -Z/k`, of
/// `u'' = (l(l+1)/r² - 2Z/r - k²) u` (used by the PDE solver).
///
/// `NaN` is returned where neither the series nor the asymptotic expansion
/// is accurate (`k ≤ 0`, `r < 0`, or the same domain limits as
/// [`coulomb_f`]).
#[must_use]
#[allow(clippy::float_cmp, clippy::single_call_fn)] // exact sentinel comparisons; one PDE call site
pub(crate) fn coulomb_wave_value(
    z: f64,
    l: usize,
    k: f64,
    r: f64,
) -> f64 {
    if !z.is_finite() || !k.is_finite() || k <= 0.0 || !r.is_finite() || r < 0.0 {
        return f64::NAN;
    }
    let rho = k * r;
    if rho == 0.0 {
        return 0.0;
    }
    coulomb_f(-z / k, l, rho).map_or(f64::NAN, |f| (2.0 / PI).sqrt() * f)
}

// ============================================================================
// Polylogarithm and Lerch transcendent
// ============================================================================

/// `B_{2j} / (2j)!` for `j = 1..=8`: the Euler-Maclaurin corrections of
/// [`complex_zeta`].
const ZETA_CORRECTIONS: [f64; 8] = [
    1.0 / 12.0,
    -1.0 / 720.0,
    1.0 / 30_240.0,
    -1.0 / 1_209_600.0,
    1.0 / 47_900_160.0,
    -691.0 / 1_307_674_368_000.0,
    1.0 / 74_724_249_600.0,
    -3617.0 / 10_670_622_842_880_000.0,
];

/// `base^exponent` with the principal branch of the logarithm.
#[allow(clippy::single_call_fn)]
fn complex_pow(
    base: Complex64,
    exponent: Complex64,
) -> Complex64 {
    (exponent * base.ln()).exp()
}

/// `-z` with a `+0.0` imaginary part on the negative real axis, so that
/// [`Complex64::ln`] takes the principal value `+iπ`. Rust's signed zero
/// would otherwise produce `-iπ` for a real `z > 0`, which is the wrong side
/// of the polylogarithm branch cut.
#[allow(clippy::float_cmp)]
fn negated_principal(z: Complex64) -> Complex64 {
    if z.im == 0.0 {
        Complex64::new(-z.re, 0.0)
    } else {
        Complex64::new(-z.re, -z.im)
    }
}

/// `k^(-s)` for a positive integer `k` and complex `s`.
fn complex_power_int(
    k: usize,
    s: Complex64,
) -> Complex64 {
    (-s * (k as f64).ln()).exp()
}

/// The Riemann zeta function for complex `s` on the principal branch.
///
/// Euler-Maclaurin summation for `Re s > 1`, Borwein's accelerated eta
/// series for `0 < Re s <= 1`, and the functional equation beyond (with
/// [`complex_ln_gamma`]). The pole at `s = 1` returns `+∞`.
#[allow(clippy::float_cmp)] // exact comparison against a sentinel / integer-valued input is intended
fn complex_zeta(s: Complex64) -> Complex64 {
    if s.re.is_nan() || s.im.is_nan() {
        return Complex64::new(f64::NAN, f64::NAN);
    }
    if (s - Complex64::new(1.0, 0.0)).norm() < 1e-14 {
        return Complex64::new(f64::INFINITY, 0.0);
    }
    // Non-positive integers: use the real implementation (the functional
    // equation is singular at s = 0).
    if s.im == 0.0 && s.re <= 0.0 && s.re.fract() == 0.0 {
        return Complex64::new(crate::kernels::special::riemann_zeta(s.re), 0.0);
    }
    if s.re > 1.0 {
        return zeta_euler_maclaurin(s);
    }
    if s.re > 0.0 {
        return zeta_via_eta_complex(s);
    }
    // ζ(s) = 2^s π^(s-1) sin(πs/2) Γ(1-s) ζ(1-s).
    let one_minus = Complex64::new(1.0, 0.0) - s;
    let reflected = complex_zeta(one_minus);
    let two_s = (s * std::f64::consts::LN_2).exp();
    let pi_s = ((s - Complex64::new(1.0, 0.0)) * PI.ln()).exp();
    let sine = (s * (0.5 * PI)).sin();
    let gamma = complex_ln_gamma(one_minus).exp();
    two_s * pi_s * sine * gamma * reflected
}

/// Euler-Maclaurin summation of `ζ(s)` for `Re s > 1`.
#[allow(clippy::single_call_fn)]
fn zeta_euler_maclaurin(s: Complex64) -> Complex64 {
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let terms = (32.0 + s.norm()).ceil().min(400.0) as usize;
    let mut sum = Complex64::new(0.0, 0.0);
    for k in 1..terms {
        sum += complex_power_int(k, s);
    }
    let ln_n = (terms as f64).ln();
    let n_pow = |exponent: Complex64| (exponent * ln_n).exp();
    let mut tail = n_pow(Complex64::new(1.0, 0.0) - s) / (s - Complex64::new(1.0, 0.0))
        + Complex64::new(0.5, 0.0) * n_pow(-s);
    let mut rising = s;
    let mut power = n_pow(-s - Complex64::new(1.0, 0.0));
    let inverse_n2 = 1.0 / ((terms * terms) as f64);
    for (index, &coefficient) in ZETA_CORRECTIONS.iter().enumerate() {
        tail += Complex64::new(coefficient, 0.0) * rising * power;
        #[allow(clippy::cast_precision_loss)]
        let j = (index + 1) as f64;
        rising *= (s + Complex64::new(2.0 * j - 1.0, 0.0)) * (s + Complex64::new(2.0 * j, 0.0));
        power *= Complex64::new(inverse_n2, 0.0);
    }
    sum + tail
}

/// Borwein's accelerated alternating series for the Dirichlet eta function,
/// continued to `ζ(s) = η(s) / (1 - 2^(1-s))` (`0 < Re s <= 1`).
#[allow(clippy::single_call_fn)]
fn zeta_via_eta_complex(s: Complex64) -> Complex64 {
    const N: usize = 40;
    #[allow(clippy::cast_precision_loss)]
    let n = N as f64;
    let mut d = [0.0_f64; N + 1];
    let mut term = 1.0_f64;
    let mut acc = 1.0_f64;
    d[0] = acc;
    for i in 0..N {
        #[allow(clippy::cast_precision_loss)]
        let fi = i as f64;
        term *= (n + fi) * 4.0 * (n - fi) / ((2.0 * fi + 1.0) * (2.0 * fi + 2.0));
        acc += term;
        d[i + 1] = acc;
    }
    let dn = d[N];
    let mut sum = Complex64::new(0.0, 0.0);
    for (k, &dk) in d.iter().enumerate().take(N) {
        let sign = if k % 2 == 0 { 1.0 } else { -1.0 };
        sum += Complex64::new(sign * (dk - dn), 0.0) * complex_power_int(k + 1, s);
    }
    let eta = -sum / Complex64::new(dn, 0.0);
    let denominator =
        Complex64::new(1.0, 0.0) - ((Complex64::new(1.0, 0.0) - s) * std::f64::consts::LN_2).exp();
    eta / denominator
}

/// `ζ(n)` for an integer `n`, using Bernoulli numbers for `n <= 0`.
#[allow(clippy::single_call_fn)]
fn zeta_at_integer(n: i64) -> f64 {
    if n >= 2 {
        return crate::kernels::special::riemann_zeta(n as f64);
    }
    if n == 1 {
        return f64::INFINITY;
    }
    if n == 0 {
        return -0.5;
    }
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let m = (-n) as u32;
    #[allow(clippy::cast_precision_loss)]
    let denominator = f64::from(m + 1);
    -bernoulli_number(m + 1) / denominator
}

/// Eulerian numbers `A(n, k)`, `k = 0..=n`, by
/// `A(n, k) = (n - k) A(n-1, k-1) + (k+1) A(n-1, k)`.
#[allow(clippy::single_call_fn)]
fn eulerian_numbers(n: usize) -> Vec<f64> {
    let mut row = vec![1.0];
    for size in 1..=n {
        let mut next = vec![0.0; size + 1];
        for (k, &value) in row.iter().enumerate() {
            #[allow(clippy::cast_precision_loss)]
            {
                next[k] += (k as f64 + 1.0) * value;
                next[k + 1] += (size as f64 - k as f64 - 1.0) * value;
            }
        }
        row = next;
    }
    row
}

/// `Li_{-n}(z) = z A_n(z) / (1 - z)^(n+1)` for `n >= 0`.
#[allow(clippy::single_call_fn)]
fn polylog_nonpositive_integer(
    n: usize,
    z: Complex64,
) -> Complex64 {
    let coefficients = eulerian_numbers(n);
    let mut polynomial = Complex64::new(0.0, 0.0);
    for &coefficient in coefficients.iter().rev() {
        polynomial = polynomial * z + Complex64::new(coefficient, 0.0);
    }
    let one_minus = Complex64::new(1.0, 0.0) - z;
    #[allow(clippy::cast_possible_wrap)]
    let power = one_minus.powi((n + 1) as i32);
    z * polynomial / power
}

/// Direct power series `Σ z^k / k^s` (`|z| <= 1/2`).
fn polylog_series(
    s: Complex64,
    z: Complex64,
) -> Complex64 {
    let mut term = z;
    let mut sum = z;
    for k in 2..100_000_usize {
        term *= z;
        let add = term * complex_power_int(k, s);
        sum += add;
        if add.norm() < 1e-17 * sum.norm() {
            break;
        }
    }
    sum
}

/// `B_n(x)` for a complex `x`, from the Bernoulli numbers.
#[allow(clippy::single_call_fn)]
fn bernoulli_polynomial(
    n: u32,
    x: Complex64,
) -> Complex64 {
    let mut sum = Complex64::new(0.0, 0.0);
    let mut binomial = 1.0_f64;
    let mut power = Complex64::new(1.0, 0.0);
    // B_n(x) = Σ_{j=0}^n C(n, j) B_{n-j} x^j.
    for j in 0..=n {
        sum += Complex64::new(binomial * bernoulli_number(n - j), 0.0) * power;
        if j < n {
            binomial *= f64::from(n - j) / f64::from(j + 1);
            power *= x;
        }
    }
    sum
}

/// `Li_n(z)` for a positive integer `n` and `|z| > 1/2`.
///
/// For `|z| <= 6/5` the logarithmic expansion
/// `Li_n(z) = μ^(n-1)/(n-1)! (H_{n-1} - ln(-μ)) + Σ_{k≠n-1} ζ(n-k) μ^k/k!`
/// (`μ = ln z`, `|μ| < 2π`) converges quickly; beyond it Kummer's relation
/// `Li_n(z) + (-1)^n Li_n(1/z) = -(2πi)^n/n! B_n(1/2 + ln(-z)/(2πi))`
/// reduces the argument to the unit disk.
#[allow(clippy::single_call_fn, clippy::float_cmp)]
fn polylog_positive_integer(
    n: u32,
    z: Complex64,
) -> Complex64 {
    let mu = z.ln();
    if z.norm() <= 1.2 {
        #[allow(clippy::cast_possible_wrap, clippy::cast_precision_loss)]
        let order = n as usize;
        let mut sum = Complex64::new(0.0, 0.0);
        let mut factorial = 1.0_f64;
        let mut mu_power = Complex64::new(1.0, 0.0);
        let mut term_norm = f64::INFINITY;
        for k in 0..=600_usize {
            if k != order - 1 {
                let index = i64::from(n) - i64::try_from(k).unwrap_or(i64::MAX);
                let zeta = zeta_at_integer(index);
                if zeta.is_infinite() {
                    return Complex64::new(f64::INFINITY, 0.0);
                }
                let term = Complex64::new(zeta * mu_power.re, zeta * mu_power.im) / factorial;
                sum += term;
                // A zero here is a trivial zero of ζ at a negative even
                // integer; it says nothing about convergence.
                term_norm = if zeta == 0.0 {
                    f64::INFINITY
                } else {
                    term.norm()
                };
            }
            #[allow(clippy::cast_precision_loss)]
            let next = (k + 1) as f64;
            factorial *= next;
            mu_power *= mu;
            if k > order + 5 && term_norm < 1e-18 * sum.norm() {
                break;
            }
        }
        #[allow(clippy::cast_precision_loss)]
        let harmonic: f64 = (1..n).map(|j| 1.0 / f64::from(j)).sum();
        #[allow(clippy::cast_possible_wrap)]
        let factorial_n = {
            let mut value = 1.0_f64;
            for j in 1..n {
                value *= f64::from(j);
            }
            value
        };
        let neg_mu = negated_principal(mu);
        let log_neg_mu = neg_mu.ln();
        #[allow(clippy::cast_possible_wrap, clippy::cast_precision_loss)]
        let log_term = mu.powi(n as i32 - 1)
            * Complex64::new(harmonic - log_neg_mu.re, -log_neg_mu.im)
            / factorial_n;
        return log_term + sum;
    }
    // Kummer's relation.
    let i_two_pi = Complex64::new(0.0, 2.0 * PI);
    let theta = Complex64::new(0.5, 0.0) + negated_principal(z).ln() / i_two_pi;
    let mut factorial = 1.0_f64;
    for k in 1..=n {
        factorial *= f64::from(k);
    }
    #[allow(clippy::cast_possible_wrap)]
    let leading = -(i_two_pi.powi(n as i32) * bernoulli_polynomial(n, theta) / factorial);
    let sign = if n.is_multiple_of(2) {
        1.0
    } else {
        -1.0
    };
    leading
        - Complex64::new(sign, 0.0)
            * polylog_series(
                Complex64::new(f64::from(n), 0.0),
                Complex64::new(1.0, 0.0) / z,
            )
}

/// The logarithmic expansion of `Li_s(z)` for non-integer `s`:
/// `Li_s(z) = Γ(1-s)(-ln z)^(s-1) + Σ_{k≥0} ζ(s-k) (ln z)^k / k!`
/// (DLMF 25.12.10), valid for `|ln z| < 2π` with the principal branches.
#[allow(clippy::single_call_fn)]
fn polylog_log_expansion(
    s: Complex64,
    z: Complex64,
) -> Complex64 {
    let mu = z.ln();
    let mut sum = complex_ln_gamma(Complex64::new(1.0, 0.0) - s).exp()
        * complex_pow(negated_principal(mu), s - Complex64::new(1.0, 0.0));
    let mut factorial = 1.0_f64;
    let mut mu_power = Complex64::new(1.0, 0.0);
    let mut last = f64::INFINITY;
    for k in 0..=600_usize {
        #[allow(clippy::cast_precision_loss)]
        let zeta = complex_zeta(s - Complex64::new(k as f64, 0.0));
        let term = zeta * mu_power / factorial;
        sum += term;
        if k > 10 && term.norm() > last && mu.norm() > 1.0 {
            break;
        }
        last = term.norm();
        if k > 10 && term.norm() < 1e-18 * sum.norm() {
            break;
        }
        #[allow(clippy::cast_precision_loss)]
        let next = (k + 1) as f64;
        factorial *= next;
        mu_power *= mu;
    }
    sum
}

/// The polylogarithm `Li_s(z)` for complex `s` and `z`, on the principal
/// branch (cut along `[1, ∞)`).
///
/// * `|z| <= 1/2`: the defining power series;
/// * non-positive integer `s = -n`: the rational function
///   `z A_n(z)/(1-z)^(n+1)` with Eulerian numbers `A_n`;
/// * positive integer `s`: the logarithmic (Bernoulli/Hurwitz) continuation
///   for `|z| <= 3/2` and Kummer's inversion relation beyond, so every `z` is
///   covered;
/// * other `s`: the logarithmic expansion about `z = 1` for `|ln z| < 4.5`
///   (the expansion is valid for `|ln z| < 2π`, but `ζ(s-k)` overflows
///   `f64` before it converges near that bound);
/// * `z = 1`: `ζ(s)` (the analytic continuation, as in `mpmath`).
///
/// Outside the stated domain (non-integer `s` with `|ln z| >= 4.5`) the
/// result is `NaN + NaN i`. Accuracy is about `1e-13` for moderate
/// parameters.
#[must_use]
#[allow(clippy::float_cmp)] // exact comparison against a sentinel / integer-valued input is intended
pub fn polylog_complex(
    s: Complex64,
    z: Complex64,
) -> Complex64 {
    let nan = Complex64::new(f64::NAN, f64::NAN);
    if s.re.is_nan() || s.im.is_nan() || z.re.is_nan() || z.im.is_nan() {
        return nan;
    }
    if z.norm() == 0.0 {
        return Complex64::new(0.0, 0.0);
    }
    if z == Complex64::new(1.0, 0.0) {
        return complex_zeta(s);
    }
    // The two elementary orders are exact and avoid any cancellation in
    // the logarithmic expansion.
    let one = Complex64::new(1.0, 0.0);
    if s == one {
        return -(one - z).ln();
    }
    if s == Complex64::new(0.0, 0.0) {
        return z / (one - z);
    }
    let integer = s.im.abs() < 1e-13 && s.re.fract() == 0.0;
    if integer && s.re <= 0.0 {
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let n = (-s.re) as usize;
        if n <= 48 {
            return polylog_nonpositive_integer(n, z);
        }
        return if z.norm() <= 0.5 {
            polylog_series(s, z)
        } else {
            nan
        };
    }
    if integer && s.re > 0.0 && s.re <= 60.0 {
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let n = s.re as u32;
        if z.norm() <= 0.5 {
            return polylog_series(s, z);
        }
        return polylog_positive_integer(n, z);
    }
    if z.norm() <= 0.5 {
        return polylog_series(s, z);
    }
    let mu = z.ln();
    // The logarithmic expansion needs `|ln z| < 2π`; in practice the term
    // `ζ(s-k)` would overflow `f64` long before the series converges for
    // `|ln z|` close to that bound, so the implemented range is `|ln z| <
    // 4.5` (integer `s` is handled by the Kummer relation for every `z`).
    if mu.norm() < 4.5 {
        return polylog_log_expansion(s, z);
    }
    nan
}

/// The polylogarithm `Li_s(z)` for real `s` and real `z`.
///
/// Real on `z <= 1` (including the negative axis); for `z > 1` the principal
/// value is complex in general and `NaN` is returned, except for
/// non-positive integer `s`, where `Li_{-n}(z)` is the real rational
/// function `z A_n(z)/(1-z)^(n+1)`. At `z = 1` the analytic continuation
/// `ζ(s)` is returned (and `+∞` at the pole `s = 1`).
#[must_use]
#[allow(clippy::float_cmp)] // exact comparison against a sentinel / integer-valued input is intended
pub fn polylog(
    s: f64,
    z: f64,
) -> f64 {
    if s.is_nan() || z.is_nan() {
        return f64::NAN;
    }
    if z == 0.0 {
        return 0.0;
    }
    if z == 1.0 {
        return if s == 1.0 {
            f64::INFINITY
        } else {
            crate::kernels::special::riemann_zeta(s)
        };
    }
    let value = polylog_complex(Complex64::new(s, 0.0), Complex64::new(z, 0.0));
    if z > 1.0 {
        return if s <= 0.0 && s.fract() == 0.0 && value.im.abs() < 1e-12 {
            value.re
        } else {
            f64::NAN
        };
    }
    if value.im.abs() <= 1e-10 * value.re.abs().max(1.0) {
        value.re
    } else {
        f64::NAN
    }
}

/// `Φ(z, -n, a)` for a non-positive integer order `-n`, from
/// `(z d/dz + a)^n 1/(1-z)`, the Stirling expansion
/// `(z d/dz)^j 1/(1-z) = Σ_m S(j, m) m! z^m/(1-z)^(m+1)` and the binomial
/// theorem for the commuting operators `z d/dz` and multiplication by `a`.
#[allow(clippy::single_call_fn)]
fn lerch_nonpositive_integer(
    z: f64,
    n: usize,
    a: f64,
) -> f64 {
    let mut stirling = vec![vec![0.0_f64; n + 1]; n + 1];
    stirling[0][0] = 1.0;
    for j in 1..=n {
        for m in 1..=j {
            #[allow(clippy::cast_precision_loss)]
            let mf = m as f64;
            stirling[j][m] = stirling[j - 1][m - 1] + mf * stirling[j - 1][m];
        }
    }
    let inverse = 1.0 / (1.0 - z);
    let mut sum = 0.0;
    let mut binomial = 1.0;
    let mut a_power = a.powi(i32::try_from(n).unwrap_or(0));
    for j in 0..=n {
        let mut inner = 0.0;
        let mut factorial = 1.0;
        let mut z_power = 1.0;
        let mut inverse_power = inverse;
        for m in 0..=j {
            #[allow(clippy::cast_precision_loss)]
            let mf = m as f64;
            inner += stirling[j][m] * factorial * z_power * inverse_power;
            factorial *= mf + 1.0;
            z_power *= z;
            inverse_power *= inverse;
        }
        sum += binomial * a_power * inner;
        binomial *= (n - j) as f64 / (j + 1) as f64;
        a_power /= a;
    }
    sum
}

/// The Lerch transcendent `Φ(z, s, a) = Σ_{k=0}^∞ z^k / (k+a)^s` for real
/// `s`, real `a > 0` and `z < 1` (the standard range).
///
/// For `|z| <= 9/10` the defining series is summed directly. For `s > 0` and
/// `|z| > 9/10` the Euler integral
/// `Φ = 1/Γ(s) ∫_0^∞ t^(s-1) e^(-a t)/(1 - z e^(-t)) dt`
/// is evaluated by double-exponential quadrature, which covers `z < 0`
/// (`|z| > 1` included) as well as `z ↑ 1`. Non-positive integer `s` use the
/// rational form `(z d/dz + a)^n 1/(1-z)`. At `z = 1` the analytic
/// continuation `ζ(s, a)` is returned; `z > 1` and non-positive non-integer
/// `s` with `|z| > 9/10` are outside the implemented range (`NaN`). Accuracy
/// is about `1e-13` for moderate parameters.
#[must_use]
#[allow(clippy::float_cmp)] // exact comparison against a sentinel / integer-valued input is intended
pub fn lerch(
    z: f64,
    s: f64,
    a: f64,
) -> f64 {
    if z.is_nan() || s.is_nan() || a.is_nan() || a <= 0.0 || z > 1.0 {
        return f64::NAN;
    }
    if z == 0.0 {
        return a.powf(-s);
    }
    if z == 1.0 {
        return if s == 1.0 {
            f64::INFINITY
        } else {
            crate::kernels::special::hurwitz_zeta(s, a)
        };
    }
    if z.abs() <= 0.9 {
        let mut sum = 0.0;
        let mut z_power = 1.0;
        for k in 0..100_000_usize {
            #[allow(clippy::cast_precision_loss)]
            let value = (k as f64 + a).powf(-s);
            let term = z_power * value;
            sum += term;
            if term.abs() < 1e-17 * sum.abs() {
                break;
            }
            z_power *= z;
        }
        return sum;
    }
    if s > 0.0 {
        // 1/Γ(s) ∫_0^∞ t^(s-1) e^(-a t) / (1 - z e^(-t)) dt, with the
        // power written through its logarithm so that large `t` cannot
        // overflow before the exponential decay sets in.
        let integrand = |t: f64| {
            if t <= 0.0 {
                return 0.0;
            }
            let log_value = (s - 1.0) * t.ln() - a * t;
            if log_value < -745.0 {
                return 0.0;
            }
            log_value.exp() / (1.0 - z * (-t).exp())
        };
        let quadrature = crate::kernels::quadrature::exp_sinh(integrand, 0.0, 1e-15);
        return quadrature.value / gamma_numerical(s);
    }
    if s <= 0.0 && s.fract() == 0.0 {
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let n = (-s) as usize;
        if n <= 64 {
            return lerch_nonpositive_integer(z, n, a);
        }
    }
    f64::NAN
}

// ============================================================================
// Appell's double hypergeometric function F1
// ============================================================================

/// `2F1(a, b; c; y)` by its power series, used for the inner sum of
/// [`appell_f1`]. The connection formula of [`hyp2f1`] becomes numerically
/// unstable when `a` and `c` are shifted far together, while this series is
/// uniformly convergent for the `|y| < 1` of the double series.
#[allow(clippy::single_call_fn)]
fn appell_inner(
    a: f64,
    b: f64,
    c: f64,
    y: f64,
) -> f64 {
    let mut term = 1.0;
    let mut sum = 1.0;
    for n in 1..100_000_usize {
        #[allow(clippy::cast_precision_loss)]
        let nf = n as f64;
        term *= (a + nf - 1.0) * (b + nf - 1.0) / ((c + nf - 1.0) * nf) * y;
        sum += term;
        if term.abs() < 1e-17 * sum.abs() {
            break;
        }
    }
    sum
}

/// Appell's first double hypergeometric function
/// `F_1(a; b_1, b_2; c; x, y)`
/// `= Σ_{m,n≥0} (a)_{m+n} (b_1)_m (b_2)_n / ((c)_{m+n} m! n!) x^m y^n`.
///
/// The double series converges for `|x| < 1`, `|y| < 1`; it is summed
/// directly when `max(|x|, |y|) <= 0.95` or when `a`/`b_1`/`b_2` is a
/// non-positive integer (the series then terminates). Otherwise Euler's
/// integral
/// `F_1 = Γ(c)/(Γ(a)Γ(c-a)) ∫_0^1 t^(a-1)(1-t)^(c-a-1)`
/// `× (1-xt)^(-b_1) (1-yt)^(-b_2) dt`
/// is evaluated by tanh-sinh quadrature, which requires `c > a > 0` and
/// `x, y < 1`. `NaN` outside those domains (`c <= 0`, `x >= 1`, `y >= 1`,
/// or no applicable representation). Accuracy is about `1e-13`.
#[must_use]
#[allow(clippy::float_cmp)] // exact comparison against a sentinel / integer-valued input is intended
pub fn appell_f1(
    a: f64,
    b1: f64,
    b2: f64,
    c: f64,
    x: f64,
    y: f64,
) -> f64 {
    if a.is_nan() || b1.is_nan() || b2.is_nan() || c.is_nan() || x.is_nan() || y.is_nan() {
        return f64::NAN;
    }
    if c <= 0.0 || x >= 1.0 || y >= 1.0 {
        return f64::NAN;
    }
    if a == 0.0 || (b1 == 0.0 && b2 == 0.0) || (x == 0.0 && y == 0.0) {
        return 1.0;
    }
    let terminating = |p: f64| p <= 0.0 && p.fract() == 0.0;
    let series = |a: f64, b1: f64, b2: f64, c: f64, x: f64, y: f64| -> f64 {
        let mut sum = 0.0;
        let mut term_x = 1.0;
        for m in 0..10_000_usize {
            #[allow(clippy::cast_precision_loss)]
            let mf = m as f64;
            let inner = appell_inner(a + mf, b2, c + mf, y);
            if !inner.is_finite() {
                return f64::NAN;
            }
            let contribution = term_x * inner;
            sum += contribution;
            if m > 5 && contribution.abs() < 1e-17 * sum.abs() {
                break;
            }
            term_x *= (a + mf) * (b1 + mf) / ((c + mf) * (mf + 1.0)) * x;
            if term_x == 0.0 {
                break;
            }
        }
        sum
    };
    if terminating(a) || terminating(b1) || terminating(b2) || x.abs().max(y.abs()) <= 0.95 {
        return series(a, b1, b2, c, x, y);
    }
    if c > a && a > 0.0 {
        let integrand = |t: f64| {
            if t <= 0.0 || t >= 1.0 {
                return 0.0;
            }
            t.powf(a - 1.0)
                * (1.0 - t).powf(c - a - 1.0)
                * (1.0 - x * t).powf(-b1)
                * (1.0 - y * t).powf(-b2)
        };
        let coefficient = gamma_numerical(c) / (gamma_numerical(a) * gamma_numerical(c - a));
        return coefficient
            * crate::kernels::quadrature::tanh_sinh(integrand, 0.0, 1.0, 1e-14).value;
    }
    f64::NAN
}

// ============================================================================
// Kelvin functions
// ============================================================================

/// Largest argument for which the Kelvin series are summed directly
/// (double-double arithmetic defeats the `e^(x/√2)` cancellation).
const KELVIN_SERIES_MAX: f64 = 30.0;
/// Largest argument for which `ker`/`kei` use the integral representation
/// (beyond it the asymptotic expansion is both faster and more accurate).
const KELVIN_INTEGRAL_MAX: f64 = 15.0;

/// Largest order of the Kelvin functions.
pub const KELVIN_MAX_ORDER: usize = 100;

/// `J_ν(x e^(3πi/4))` by its power series, in complex double-double
/// arithmetic (any integer order).
fn kelvin_j_series(
    nu: usize,
    x: f64,
) -> Complex64 {
    let a = x * std::f64::consts::FRAC_1_SQRT_2;
    let q = Cdd::new(0.0, 0.25 * x * x);
    // (z/2)^ν / ν! with z/2 = (-a + i a)/2.
    let half_z = Cdd::new(-0.5 * a, 0.5 * a);
    let mut term = Cdd::new(1.0, 0.0);
    for j in 1..=nu {
        #[allow(clippy::cast_precision_loss)]
        let jf = j as f64;
        term = term * half_z / Cdd::new(jf, 0.0);
    }
    let mut sum = term;
    for k in 1..4000_usize {
        #[allow(clippy::cast_precision_loss)]
        let kf = k as f64;
        term = term * q / Cdd::new(kf * (kf + nu as f64), 0.0);
        sum += term;
        if term.norm() < 1e-30 * sum.norm() {
            break;
        }
    }
    Complex64::new(sum.re.to_f64(), sum.im.to_f64())
}

/// `J_ν(z)` for `|z| >= 30` from the Hankel asymptotic expansion, with
/// `z = x e^(3πi/4)`.
fn kelvin_j_asymptotic(
    nu: usize,
    x: f64,
) -> Complex64 {
    let a = x * std::f64::consts::FRAC_1_SQRT_2;
    let z = Complex64::new(-a, a);
    let mu = 4.0 * (nu as f64) * (nu as f64);
    let mut term = Complex64::new(1.0, 0.0);
    let mut p = Complex64::new(1.0, 0.0);
    let mut q = Complex64::new(0.0, 0.0);
    let mut last = f64::INFINITY;
    for k in 1..400_i32 {
        let kf = f64::from(k);
        term *= Complex64::new(mu - (2.0 * kf - 1.0).powi(2), 0.0)
            / (Complex64::new(8.0 * kf, 0.0) * z);
        if term.norm() > last {
            break;
        }
        last = term.norm();
        match k % 4 {
            | 0 => p += term,
            | 1 => q += term,
            | 2 => p -= term,
            | _ => q -= term,
        }
        if term.norm() < 1e-18 {
            break;
        }
    }
    let chi = z - Complex64::new((0.5 * nu as f64 + 0.25) * PI, 0.0);
    (Complex64::new(2.0 / PI, 0.0) / z).sqrt() * (p * chi.cos() - q * chi.sin())
}

/// `J_ν(x e^(3πi/4))` for any order and `x > 30`: Miller's backward
/// recurrence `J_(k-1) = (2k/z) J_k - J_(k+1)` from an index beyond both `ν`
/// and `|z| = x` (where `J` is minimal), scaled to the Hankel expansion of
/// `J_0`, which is accurate there. (The expansion of `J_ν` itself needs
/// `x` well above `ν²`.)
fn kelvin_j_miller(
    nu: usize,
    x: f64,
) -> Complex64 {
    let a = x * std::f64::consts::FRAC_1_SQRT_2;
    let z = Complex64::new(-a, a);
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let start = nu.max(x.ceil() as usize) + 60;
    let mut y = vec![Complex64::new(0.0, 0.0); start + 2];
    y[start] = Complex64::new(1.0, 0.0);
    for k in (1..=start).rev() {
        #[allow(clippy::cast_precision_loss)]
        let factor = Complex64::new(2.0 * k as f64, 0.0) / z;
        y[k - 1] = factor * y[k] - y[k + 1];
        if y[k - 1].norm() > 1e200 {
            for v in &mut y[k - 1..] {
                *v *= 1e-200;
            }
        }
    }
    kelvin_j_asymptotic(0, x) * y[nu] / y[0]
}

/// `J_ν(x e^(3πi/4))` for any integer order.
fn kelvin_j(
    nu: usize,
    x: f64,
) -> Complex64 {
    if x <= KELVIN_SERIES_MAX {
        kelvin_j_series(nu, x)
    } else if nu <= 1 {
        kelvin_j_asymptotic(nu, x)
    } else {
        kelvin_j_miller(nu, x)
    }
}

/// `K_ν(x e^(iπ/4))` for `0 < x < 15` by
/// `∫_0^∞ e^(-z cosh t) cosh(νt) dt`, split into real and imaginary parts
/// and evaluated by adaptive Gauss-Kronrod quadrature.
#[allow(clippy::single_call_fn)]
fn kelvin_k_integral(
    nu: usize,
    x: f64,
) -> Complex64 {
    let a = x * std::f64::consts::FRAC_1_SQRT_2;
    let t_max = (45.0 / a).max(1.0).acosh();
    let envelope = |t: f64| {
        let phase = a * t.cosh();
        (-phase).exp() * (nu as f64 * t).cosh()
    };
    let real = crate::kernels::integrate::gauss_kronrod(
        |t| {
            let phase = a * t.cosh();
            envelope(t) * phase.cos()
        },
        0.0,
        t_max,
        1e-16,
        3000,
    )
    .value;
    let imaginary = crate::kernels::integrate::gauss_kronrod(
        |t| {
            let phase = a * t.cosh();
            -envelope(t) * phase.sin()
        },
        0.0,
        t_max,
        1e-16,
        3000,
    )
    .value;
    Complex64::new(real, imaginary)
}

/// `K_ν(z)` for `|z| >= 15`, `z = x e^(iπ/4)`, from
/// `K_ν(z) ~ sqrt(π/(2z)) e^(-z) Σ (μ-1)(μ-9)…(μ-(2k-1)²)/(k! (8z)^k)`.
#[allow(clippy::single_call_fn)]
fn kelvin_k_asymptotic(
    nu: usize,
    x: f64,
) -> Complex64 {
    let a = x * std::f64::consts::FRAC_1_SQRT_2;
    let z = Complex64::new(a, a);
    let mu = 4.0 * (nu as f64) * (nu as f64);
    let mut term = Complex64::new(1.0, 0.0);
    let mut sum = term;
    let mut last = f64::INFINITY;
    for k in 1..400_i32 {
        let kf = f64::from(k);
        term *= Complex64::new(mu - (2.0 * kf - 1.0).powi(2), 0.0)
            / (Complex64::new(8.0 * kf, 0.0) * z);
        if term.norm() > last {
            break;
        }
        last = term.norm();
        sum += term;
        if term.norm() < 1e-18 {
            break;
        }
    }
    (Complex64::new(PI / 2.0, 0.0) / z).sqrt() * (-z).exp() * sum
}

/// `K_ν(x e^(iπ/4))` for `ν = 0` or `1` and `0 < x`, switching from the
/// integral representation to the asymptotic expansion at `x = 15`.
fn kelvin_k(
    nu: usize,
    x: f64,
) -> Complex64 {
    if x < KELVIN_INTEGRAL_MAX {
        kelvin_k_integral(nu, x)
    } else {
        kelvin_k_asymptotic(nu, x)
    }
}

/// `K_ν(x e^(iπ/4))` for any integer order: the forward recurrence
/// `K_(n+1) = K_(n-1) + (2n/z) K_n` from `K_0` and `K_1` (stable, `K` being
/// the dominant solution).
fn kelvin_k_order(
    nu: usize,
    x: f64,
) -> Complex64 {
    if nu <= 1 {
        return kelvin_k(nu, x);
    }
    let a = x * std::f64::consts::FRAC_1_SQRT_2;
    let z = Complex64::new(a, a);
    let (mut below, mut current) = (kelvin_k(0, x), kelvin_k(1, x));
    for n in 1..nu {
        #[allow(clippy::cast_precision_loss)]
        let next = below + Complex64::new(2.0 * n as f64, 0.0) / z * current;
        below = current;
        current = next;
    }
    current
}

/// `e^(-iπν/2) K_ν(x e^(iπ/4)) = ker_ν(x) + i kei_ν(x)`.
fn kelvin_ker_kei(
    nu: usize,
    x: f64,
) -> Complex64 {
    let v = kelvin_k_order(nu, x);
    match nu % 4 {
        | 0 => v,
        | 1 => Complex64::new(v.im, -v.re),
        | 2 => -v,
        | _ => Complex64::new(-v.im, v.re),
    }
}

/// The limit of `ker_ν(x)` (`imaginary = false`) or `kei_ν(x)` as `x -> 0+`
/// for `ν >= 1`: the leading term of `K_ν(z) ~ (ν-1)!/2 (2/z)^ν` times
/// `e^(-iπν/2)` has the phase `e^(-3πiν/4)`; `NaN` when the component of
/// that phase vanishes (the limit is then finite and not computed).
fn kelvin_origin_limit(
    nu: usize,
    imaginary: bool,
) -> f64 {
    #[allow(clippy::cast_precision_loss)]
    let phase = -0.75 * PI * nu as f64;
    let component = if imaginary {
        phase.sin()
    } else {
        phase.cos()
    };
    if component.abs() < 1e-9 {
        f64::NAN
    } else {
        f64::INFINITY.copysign(component)
    }
}

/// The Kelvin function `ber_ν(x) = Re J_ν(x e^(3πi/4))` for an integer
/// order `ν <= 100` and `x >= 0`.
///
/// Power series in complex double-double arithmetic for `x <= 30` (the
/// terms reach `e^(x/√2)`, so ordinary `f64` loses all accuracy near the
/// switch), beyond it the Hankel asymptotic expansion for `ν = 0, 1` and
/// Miller's backward recurrence scaled to the expansion of `J_0` for higher
/// orders; `ker`/`kei` use the integral representation below `x = 15` and
/// their asymptotic expansion above for `ν = 0, 1`, and the forward
/// recurrence `K_(n+1) = K_(n-1) + (2n/z) K_n` for higher orders. `NaN` for
/// negative `x` or `ν > 100`.
#[must_use]
#[allow(clippy::float_cmp)] // exact comparison against a sentinel / integer-valued input is intended
pub fn kelvin_ber(
    nu: usize,
    x: f64,
) -> f64 {
    if x.is_nan() || x < 0.0 || nu > KELVIN_MAX_ORDER {
        return f64::NAN;
    }
    if x == 0.0 {
        return if nu == 0 { 1.0 } else { 0.0 };
    }
    kelvin_j(nu, x).re
}

/// The Kelvin function `bei_ν(x) = Im J_ν(x e^(3πi/4))`; domain and
/// accuracy as [`kelvin_ber`].
#[must_use]
#[allow(clippy::float_cmp)] // exact comparison against a sentinel / integer-valued input is intended
pub fn kelvin_bei(
    nu: usize,
    x: f64,
) -> f64 {
    if x.is_nan() || x < 0.0 || nu > KELVIN_MAX_ORDER {
        return f64::NAN;
    }
    if x == 0.0 {
        return 0.0;
    }
    kelvin_j(nu, x).im
}

/// The Kelvin function
/// `ker_ν(x) + i kei_ν(x) = e^(-iπν/2) K_ν(x e^(iπ/4))` for an integer order
/// `ν <= 100` and `x > 0`.
///
/// `ker_ν(x)` is its real part; see [`kelvin_ber`] for the algorithm. At
/// `x = 0` the value is `+inf` for `ν = 0` and, for `ν >= 1`, the infinity
/// of the sign of the leading term `cos(3πν/4)` (`NaN` where that vanishes
/// and the limit is finite).
#[must_use]
#[allow(clippy::float_cmp)] // exact comparison against a sentinel / integer-valued input is intended
pub fn kelvin_ker(
    nu: usize,
    x: f64,
) -> f64 {
    if x.is_nan() || x < 0.0 || nu > KELVIN_MAX_ORDER {
        return f64::NAN;
    }
    if x == 0.0 {
        return if nu == 0 {
            f64::INFINITY
        } else {
            kelvin_origin_limit(nu, false)
        };
    }
    kelvin_ker_kei(nu, x).re
}

/// The Kelvin function `kei_ν(x)`, the imaginary part of
/// `e^(-iπν/2) K_ν(x e^(iπ/4))`; domain and accuracy as [`kelvin_ker`].
#[must_use]
#[allow(clippy::float_cmp)] // exact comparison against a sentinel / integer-valued input is intended
pub fn kelvin_kei(
    nu: usize,
    x: f64,
) -> f64 {
    if x.is_nan() || x < 0.0 || nu > KELVIN_MAX_ORDER {
        return f64::NAN;
    }
    if x == 0.0 {
        // K_0(z) ~ -ln(z/2) - γ has the finite imaginary limit -π/4.
        return if nu == 0 {
            -0.25 * PI
        } else {
            kelvin_origin_limit(nu, true)
        };
    }
    kelvin_ker_kei(nu, x).im
}

// ============================================================================
// Anger and Weber functions
// ============================================================================

/// The Anger function
/// `J_ν(x) = (1/π) ∫_0^π cos(νθ - x sin θ) dθ` for real `ν` and `x`.
///
/// Evaluated by adaptive Gauss-Kronrod quadrature (DLMF §11.10). For
/// integer `ν = n` it coincides with the Bessel function `J_n(x)`; for
/// large `|x|` the quadrature cost grows with the number of oscillations.
#[must_use]
pub fn anger_j(
    nu: f64,
    x: f64,
) -> f64 {
    if nu.is_nan() || x.is_nan() {
        return f64::NAN;
    }
    let integrand = |theta: f64| (nu * theta - x * theta.sin()).cos();
    crate::kernels::integrate::gauss_kronrod(integrand, 0.0, PI, 1e-14, 4000).value / PI
}

/// The Weber function
/// `E_ν(x) = (1/π) ∫_0^π sin(νθ - x sin θ) dθ` for real `ν` and `x`.
///
/// Same method and domain notes as [`anger_j`].
#[must_use]
pub fn weber_e(
    nu: f64,
    x: f64,
) -> f64 {
    if nu.is_nan() || x.is_nan() {
        return f64::NAN;
    }
    let integrand = |theta: f64| (nu * theta - x * theta.sin()).sin();
    crate::kernels::integrate::gauss_kronrod(integrand, 0.0, PI, 1e-14, 4000).value / PI
}
