//! # Hypergeometric functions with complex parameters and arguments
//!
//! * [`hyp_pfq`]: the generalised hypergeometric function `pFq(a; b; z)` for
//!   complex parameter vectors and complex `z`;
//! * [`hyp2f1_cplx`], [`hyp1f1_cplx`] and [`hyp_u_cplx`]: Gauss `2F1`, Kummer
//!   `M(a, b, z) = 1F1` and Tricomi `U(a, b, z)` with every argument complex;
//! * [`gamma_c`] / [`ln_gamma_c`]: the complex gamma function (Lanczos).
//!
//! ## Method
//!
//! `pFq` satisfies the linear differential equation
//! `[theta (theta + b_1 - 1) ... (theta + b_q - 1) - z (theta + a_1) ... (theta + a_p)] w = 0`
//! with `theta = z d/dz`, of order `n = max(p, q + 1)`. For `|z| <= 1/2` the
//! defining series is summed. Otherwise the function is *continued
//! analytically* by Taylor stepping of this equation: from `z = 1/2` (where
//! the value and its first `n - 1` derivatives are summed from the shifted
//! series `pFq(a + k; b + k; z)`) along a polygonal path to `z`, each step
//! limited to `0.4` times the distance to the nearest singular point (`0`,
//! and `1` when `p = q + 1`). This gives uniform accuracy for `|z| > 1`, `z`
//! close to the singular point `1`, and for entire cases (`p <= q`) with large
//! `|z|`, without the cancellation of the direct series. The value on the
//! cut `[1, inf)` of `p = q + 1` is the limit from the half plane of `Im z`
//! (the upper one for real `z`).
//!
//! * `z = 1` with `p = q + 1`: Gauss's theorem for `2F1`; for higher `p`
//!   the slowly (algebraically) convergent series is accelerated with the
//!   Levin `u` transform of its partial sums.
//! * Terminating series (a numerator parameter equal to a non-positive
//!   integer) are summed exactly.
//! * `U(a, b, z)`: the asymptotic expansion for large `|z|`; for
//!   `Re z >= 0` the solution is continued *inwards* along a ray from a point
//!   where the expansion is accurate (the recessive solution grows in that
//!   direction, so the continuation is stable, also for integer `b` and for
//!   parameters of modulus 60 and more); for `Re z < 0` the value on the
//!   imaginary axis at the same modulus is continued along the circle
//!   `|z| = const` into the left half plane, where the competing solution
//!   `e^z z^{a-b}` is damped. No connection formula (and hence no
//!   cancellation near integer `b`) is used.
//!
//! Limitations: `pFq` and `M` with parameters of very large modulus (more than
//! about 50) lose accuracy; `b_j` equal to a non-positive integer (without termination) and
//! divergent cases (`p > q + 1`) give NaN.

use num_complex::Complex64;

use crate::kernels::exact_eq;

type C = Complex64;

const ONE: C = C::new(1.0, 0.0);
const ZERO: C = C::new(0.0, 0.0);

const fn nan() -> C {
    C::new(f64::NAN, f64::NAN)
}

fn is_nonpos_int_c(z: C) -> Option<u32> {
    if z.im.abs() < 1e-14 * z.re.abs().max(1.0) && z.re <= 0.0 && (z.re - z.re.round()).abs() < 1e-14 * z.re.abs().max(1.0)
    {
        Some((-z.re).round() as u32)
    } else {
        None
    }
}

// ------------------------------------------------------------ complex gamma

const LANCZOS: [f64; 9] = [
    0.999_999_999_999_809_9,
    676.520_368_121_885_1,
    -1_259.139_216_722_402_8,
    771.323_428_777_653_1,
    -176.615_029_162_140_6,
    12.507_343_278_686_905,
    -0.138_571_095_265_720_12,
    9.984_369_578_019_572e-6,
    1.505_632_735_149_311_6e-7,
];

/// Logarithm of the complex gamma function.
///
/// Lanczos approximation with relative accuracy about `1e-15`. For
/// `Re z < 1/2` the reflection formula is used, so the imaginary part is
/// determined only modulo `2 pi`.
#[must_use]
pub fn ln_gamma_c(z: C) -> C {
    if z.re < 0.5 {
        // ln Gamma(z) = ln(pi / sin(pi z)) - ln Gamma(1 - z)
        let s = (z * std::f64::consts::PI).sin();
        return (C::new(std::f64::consts::PI, 0.0) / s).ln() - ln_gamma_c(ONE - z);
    }
    let z = z - 1.0;
    let mut x = C::new(LANCZOS[0], 0.0);
    for (i, g) in LANCZOS.iter().enumerate().skip(1) {
        x += g / (z + i as f64);
    }
    let t = z + 7.5;
    0.5 * (2.0 * std::f64::consts::PI).ln() + (z + 0.5) * t.ln() - t + x.ln()
}

/// The complex gamma function; poles (non-positive integers) give NaN.
#[must_use]
pub fn gamma_c(z: C) -> C {
    if is_nonpos_int_c(z).is_some() {
        return nan();
    }
    if z.re < 0.5 {
        let s = (z * std::f64::consts::PI).sin();
        return C::new(std::f64::consts::PI, 0.0) / (s * gamma_c(ONE - z));
    }
    ln_gamma_c(z).exp()
}

// --------------------------------------------------------- Taylor ODE engine

struct HypOde {
    n: usize,
    a_coef: Vec<C>,
    b_coef: Vec<C>,
    singular_one: bool,
}

fn poly_mul_linear(p: &[C], shift: C) -> Vec<C> {
    // p(theta) * (theta + shift)
    let mut out = vec![ZERO; p.len() + 1];
    for (i, c) in p.iter().enumerate() {
        out[i + 1] += c;
        out[i] += c * shift;
    }
    out
}

fn stirling2(nmax: usize) -> Vec<Vec<f64>> {
    let mut s = vec![vec![0.0; nmax + 1]; nmax + 1];
    s[0][0] = 1.0;
    for m in 1..=nmax {
        for j in 1..=m {
            s[m][j] = j as f64 * s[m - 1][j] + s[m - 1][j - 1];
        }
    }
    s
}

fn binom(n: usize, k: usize) -> f64 {
    if k > n {
        return 0.0;
    }
    let mut r = 1.0;
    for i in 0..k {
        r *= (n - i) as f64 / (i + 1) as f64;
    }
    r
}

fn falling(k: usize, j: usize) -> f64 {
    if k < j {
        return 0.0;
    }
    let mut r = 1.0;
    for i in 0..j {
        r *= (k - i) as f64;
    }
    r
}

impl HypOde {
    fn new(a: &[C], b: &[C]) -> Self {
        let (p, q) = (a.len(), b.len());
        let n = p.max(q + 1);
        // P(theta) = theta prod (theta + b_j - 1), Q(theta) = prod (theta + a_i)
        let mut pp = vec![ZERO, ONE];
        for bj in b {
            pp = poly_mul_linear(&pp, bj - 1.0);
        }
        let mut qq = vec![ONE];
        for ai in a {
            qq = poly_mul_linear(&qq, *ai);
        }
        let st = stirling2(n);
        let mut a_coef = vec![ZERO; n + 1];
        let mut b_coef = vec![ZERO; n + 1];
        for j in 0..=n {
            for (m, c) in pp.iter().enumerate() {
                a_coef[j] += c * st[m][j];
            }
            for (m, c) in qq.iter().enumerate() {
                b_coef[j] -= c * st[m][j];
            }
        }
        Self { n, a_coef, b_coef, singular_one: p == q + 1 }
    }

    fn radius(&self, p: C) -> f64 {
        let r0 = p.norm();
        if self.singular_one { r0.min((ONE - p).norm()) } else { r0 }
    }

    /// Taylor coefficients (`c_d = w^(d)/d!`, `d < n`) at `p + dz` from those
    /// at `p`, with the scaling radius `rho`.
    fn step(&self, p: C, c: &[C], dz: C, rho: f64) -> Vec<C> {
        let n = self.n;
        let kmax = 130usize;
        // alpha_{j,i}: coefficients of (p + x)^j (A_j + B_j (p + x)) in x
        let mut beta = vec![vec![ZERO; n + 2]; n + 1];
        for j in 0..=n {
            let r0 = self.a_coef[j] + self.b_coef[j] * p;
            let r1 = self.b_coef[j];
            for i in 0..=j + 1 {
                let mut v = ZERO;
                if i <= j {
                    v += r0 * binom(j, i) * p.powi((j - i) as i32);
                }
                if i >= 1 && i - 1 <= j {
                    v += r1 * binom(j, i - 1) * p.powi((j + 1 - i) as i32);
                }
                beta[j][i] = v * rho.powi(i as i32 - j as i32);
            }
        }
        let mut ct = vec![ZERO; kmax + n + 1];
        for d in 0..n {
            ct[d] = c[d] * rho.powi(d as i32);
        }
        let lead = beta[n][0];
        for m in 0..kmax {
            let mut sum = ZERO;
            for j in 0..=n {
                for i in 0..=(j + 1).min(m + j) {
                    if j == n && i == 0 {
                        continue;
                    }
                    let k = m + j - i;
                    let ff = falling(k, j);
                    if ff != 0.0 {
                        sum += beta[j][i] * ct[k] * ff;
                    }
                }
            }
            ct[m + n] = -sum / (lead * falling(m + n, n));
        }
        let tau = dz / rho;
        let mut out = vec![ZERO; n];
        for d in 0..n {
            let mut s = ZERO;
            let mut tp = ONE; // tau^(k-d)
            for k in d..ct.len() {
                s += ct[k] * binom(k, d) * tp;
                tp *= tau;
            }
            out[d] = s / rho.powi(d as i32);
        }
        out
    }

    /// Walks from `p` (Taylor data `c`) through the waypoints.
    fn walk(&self, mut p: C, mut c: Vec<C>, waypoints: &[C]) -> Vec<C> {
        for &target in waypoints {
            let mut guard = 0;
            while (target - p).norm() > 0.0 && guard < 200_000 {
                guard += 1;
                let rho = self.radius(p);
                let rem = target - p;
                let len = rem.norm().min(0.4 * rho).min(8.0);
                let dir = rem / rem.norm();
                let dz = dir * len;
                c = self.step(p, &c, dz, rho);
                p += dz;
                if (target - p).norm() < 1e-15 * (1.0 + target.norm()) {
                    p = target;
                    break;
                }
            }
        }
        c
    }
}

/// Direct power series of `pFq` (valid for `|z| < 1` or entire cases).
fn pfq_series(a: &[C], b: &[C], z: C) -> C {
    let mut term = ONE;
    let mut sum = ONE;
    for k in 0..400_000usize {
        let kf = k as f64;
        let mut r = z / (kf + 1.0);
        for ai in a {
            r *= ai + kf;
        }
        for bj in b {
            r /= bj + kf;
        }
        term *= r;
        sum += term;
        let tn = term.norm();
        if tn == 0.0 || (r.norm() < 0.5 && tn < 1e-17 * sum.norm() && kf > a.len() as f64 + 2.0) {
            break;
        }
    }
    sum
}

/// Terminating series when a numerator parameter is a non-positive integer.
fn pfq_polynomial(a: &[C], b: &[C], z: C) -> Option<C> {
    let nterm = a.iter().filter_map(|x| is_nonpos_int_c(*x)).min()?;
    // the denominators must not vanish before the series terminates
    for bj in b {
        if let Some(m) = is_nonpos_int_c(*bj)
            && m < nterm
        {
            return Some(nan());
        }
        if let Some(m) = is_nonpos_int_c(*bj)
            && m == nterm
        {
            // 0/0 at the last term: terminates one term earlier by convention
            return Some(nan());
        }
    }
    let mut term = ONE;
    let mut sum = ONE;
    for k in 0..nterm as usize {
        let kf = k as f64;
        let mut r = z / (kf + 1.0);
        for ai in a {
            r *= ai + kf;
        }
        for bj in b {
            r /= bj + kf;
        }
        term *= r;
        sum += term;
    }
    Some(sum)
}

/// Levin `u` transform of a sequence of terms (complex), returning the
/// extrapolated sum of the first `terms.len()` terms' tail behaviour.
fn levin_u(terms: &[C], order: usize) -> C {
    let n = terms.len();
    let mut partial = Vec::with_capacity(n);
    let mut s = ZERO;
    for t in terms {
        s += t;
        partial.push(s);
    }
    let k = order.min(n - 1);
    let beta = 1.0;
    let (mut num, mut den) = (ZERO, ZERO);
    for j in 0..=k {
        let idx = j;
        let nn = beta + (idx + 1) as f64;
        let omega = terms[idx] * nn;
        let sign = if j % 2 == 0 { 1.0 } else { -1.0 };
        let w = sign * binom(k, j) * nn.powi(k as i32 - 1) / omega;
        num += w * partial[idx];
        den += w;
    }
    num / den
}

fn pfq_at_one(a: &[C], b: &[C]) -> C {
    let s: C = b.iter().sum::<C>() - a.iter().sum::<C>();
    if s.re <= 0.0 {
        return C::new(f64::INFINITY, 0.0);
    }
    if a.len() == 2 && b.len() == 1 {
        let (aa, bb, cc) = (a[0], a[1], b[0]);
        return gamma_c(cc) * gamma_c(cc - aa - bb) / (gamma_c(cc - aa) * gamma_c(cc - bb));
    }
    // terms t_k, accelerated by Levin's u transform from several starting points
    let nt = 48usize;
    let mut terms = Vec::with_capacity(nt);
    let mut t = ONE;
    for k in 0..nt {
        terms.push(t);
        let kf = k as f64;
        let mut r = ONE / (kf + 1.0);
        for ai in a {
            r *= ai + kf;
        }
        for bj in b {
            r /= bj + kf;
        }
        t *= r;
    }
    let mut best = levin_u(&terms, 8);
    let mut prev = best;
    let mut best_diff = f64::INFINITY;
    for order in 9..=18 {
        let v = levin_u(&terms, order);
        let d = (v - prev).norm();
        if d < best_diff {
            best_diff = d;
            best = v;
        }
        prev = v;
    }
    best
}

/// The generalised hypergeometric function `pFq(a; b; z)`.
///
/// See the module documentation for the algorithm. Returns NaN for
/// divergent series (`p > q + 1` without termination), a denominator
/// parameter at a non-positive integer, and `+inf` at `z = 1` when
/// `Re(sum b - sum a) <= 0`.
#[must_use]
pub fn hyp_pfq(a: &[C], b: &[C], z: C) -> C {
    if a.iter().chain(b).any(|x| x.re.is_nan() || x.im.is_nan()) || z.re.is_nan() || z.im.is_nan() {
        return nan();
    }
    if let Some(v) = pfq_polynomial(a, b, z) {
        return v;
    }
    if b.iter().any(|x| is_nonpos_int_c(*x).is_some()) {
        return nan();
    }
    if a.len() > b.len() + 1 {
        return nan();
    }
    if z.norm() == 0.0 {
        return ONE;
    }
    if z.norm() <= 0.5 {
        return pfq_series(a, b, z);
    }
    let ode = HypOde::new(a, b);
    if ode.singular_one && exact_eq(z.im, 0.0) && exact_eq(z.re, 1.0) {
        return pfq_at_one(a, b);
    }
    // Taylor data at z0 = 1/2
    let z0 = C::new(0.5, 0.0);
    let mut c = Vec::with_capacity(ode.n);
    let mut poch = ONE; // prod (a)_d / ((b)_d d!)
    for d in 0..ode.n {
        let df = d as f64;
        let a2: Vec<C> = a.iter().map(|x| x + df).collect();
        let b2: Vec<C> = b.iter().map(|x| x + df).collect();
        c.push(poch * pfq_series(&a2, &b2, z0));
        let mut r = ONE / (df + 1.0);
        for ai in a {
            r *= ai + df;
        }
        for bj in b {
            r /= bj + df;
        }
        poch *= r;
    }
    let side = if z.im < 0.0 { -1.0 } else { 1.0 };
    let waypoints = [C::new(0.5, 0.5 * side), z];
    let out = ode.walk(z0, c, &waypoints);
    out[0]
}

/// Gauss hypergeometric function `2F1(a, b; c; z)` for complex parameters
/// and argument (principal branch, cut along `[1, inf)`; the value on the
/// cut is the limit from `Im z -> 0+`).
#[must_use]
pub fn hyp2f1_cplx(a: C, b: C, c: C, z: C) -> C {
    hyp_pfq(&[a, b], &[c], z)
}

/// Kummer's confluent hypergeometric function `M(a, b, z) = 1F1(a; b; z)`.
///
/// All arguments are complex. For `Re z < 0` Kummer's transformation
/// `M(a, b, z) = e^z M(b - a, b, -z)` is applied so that the dominant
/// solution is the one being continued.
#[must_use]
pub fn hyp1f1_cplx(a: C, b: C, z: C) -> C {
    if z.re < 0.0 && is_nonpos_int_c(a).is_none() {
        return z.exp() * hyp_pfq(&[b - a], &[b], -z);
    }
    hyp_pfq(&[a], &[b], z)
}

/// Asymptotic series `z^{-a} sum_k (a)_k (a-b+1)_k / k! (-z)^{-k}`; `None` if
/// it does not reach `1e-15` relative accuracy.
fn u_asymptotic(a: C, b: C, z: C) -> Option<C> {
    let a1 = a - b + 1.0;
    let mut term = ONE;
    let mut sum = ONE;
    let mut last = f64::INFINITY;
    for k in 0..400usize {
        let kf = k as f64;
        term *= (a + kf) * (a1 + kf) / ((kf + 1.0) * (-z));
        let tn = term.norm();
        if tn > last && kf > 3.0 {
            return None;
        }
        last = tn;
        sum += term;
        if tn < 1e-17 * sum.norm() {
            return Some((-a * z.ln()).exp() * sum);
        }
        if tn == 0.0 {
            return Some((-a * z.ln()).exp() * sum);
        }
    }
    None
}

/// Tricomi's confluent hypergeometric function `U(a, b, z)` for complex
/// arguments on the principal branch (cut along the negative real axis;
/// the upper side is returned for `z < 0`).
#[must_use]
pub fn hyp_u_cplx(a: C, b: C, z: C) -> C {
    if a.re.is_nan() || b.re.is_nan() || z.re.is_nan() || a.im.is_nan() || b.im.is_nan() || z.im.is_nan() {
        return nan();
    }
    if z.norm() == 0.0 {
        return if b.re < 1.0 { gamma_c(ONE - b) / gamma_c(a - b + 1.0) } else { C::new(f64::INFINITY, 0.0) };
    }
    if let Some(n) = is_nonpos_int_c(a) {
        // U(-n, b, z) = (-1)^n (b)_n M(-n, b, z)
        let mut poch = ONE;
        for k in 0..n {
            poch *= b + f64::from(k);
        }
        let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
        return sign * poch * hyp_pfq(&[a], &[b], z);
    }
    if let Some(n) = is_nonpos_int_c(a - b + 1.0) {
        // U(a, b, z) = z^{1-b} U(a-b+1, 2-b, z)
        let b2 = C::new(2.0, 0.0) - b;
        let a2 = a - b + 1.0;
        let mut poch = ONE;
        for k in 0..n {
            poch *= b2 + f64::from(k);
        }
        let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
        return ((ONE - b) * z.ln()).exp() * sign * poch * hyp_pfq(&[a2], &[b2], z);
    }
    if let Some(v) = u_asymptotic(a, b, z) {
        return v;
    }
    if z.re >= 0.0 {
        // continue inwards from a ray point where the expansion is accurate
        let dir = z / z.norm();
        let mut r = 40.0_f64.max(z.norm());
        for _ in 0..12 {
            r *= 2.0;
            let zs = dir * r;
            if let (Some(u0), Some(u1)) = (u_asymptotic(a, b, zs), u_asymptotic(a + 1.0, b + 1.0, zs)) {
                // U' = -a U(a+1, b+1, z)
                let d = -a * u1;
                let ode = HypOde::new(&[a], &[b]);
                let c = ode.walk(zs, vec![u0, d], &[z]);
                return c[0];
            }
        }
        return nan();
    }
    // Re z < 0: start on the imaginary axis at the same modulus (reached by the
    // stable inward continuation above) and walk along the circle |z| = const
    // into the left half plane. Along that arc Re z decreases, so the
    // other solution of the equation, which behaves like e^z, is damped and the
    // algebraically behaving Tricomi function stays dominant.
    let rho = z.norm();
    let side = if z.im < 0.0 { -1.0 } else { 1.0 };
    let zb = C::new(0.0, side * rho);
    let u0 = hyp_u_cplx(a, b, zb);
    let u1 = hyp_u_cplx(a + 1.0, b + 1.0, zb);
    let phi_b = side * std::f64::consts::FRAC_PI_2;
    let phi = if z.im == 0.0 { std::f64::consts::PI } else { z.arg() };
    let steps = ((phi - phi_b).abs() / 0.15).ceil().max(1.0) as usize;
    let mut way: Vec<C> = (1..steps)
        .map(|k| C::from_polar(rho, phi_b + (phi - phi_b) * k as f64 / steps as f64))
        .collect();
    way.push(z);
    let ode = HypOde::new(&[a], &[b]);
    ode.walk(zb, vec![u0, -a * u1], &way)[0]
}
