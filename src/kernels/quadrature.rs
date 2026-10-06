//! # Advanced quadrature
//!
//! * adaptive Gauss-Kronrod (G7K15 and G10K21) with error control,
//!   semi-infinite and infinite ranges via variable transformations;
//! * double-exponential rules (tanh-sinh, exp-sinh, sinh-sinh) for endpoint
//!   singularities and infinite ranges;
//! * Clenshaw-Curtis rules;
//! * Filon-type quadrature for oscillatory integrands;
//! * multi-dimensional integration: nested adaptive Gauss-Kronrod, plain
//!   Monte Carlo and quasi-Monte Carlo with Halton and Sobol sequences;
//! * adaptive Genz-Malik cubature (embedded degree 7/5 rule,
//!   error-driven subdivision) on hyperrectangles.
//!
//! Everything is deterministic: random sampling uses a seeded generator.
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

use std::cell::RefCell;
use std::f64::consts::PI;

use crate::kernels::random::Rng;

/// Result of an integration: value, error estimate and cost.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Integral {
    /// Estimated value of the integral.
    pub value: f64,
    /// Estimated absolute error.
    pub error: f64,
    /// Number of integrand evaluations.
    pub evaluations: usize,
    /// Whether the requested tolerance was reached.
    pub converged: bool,
}

/// Gauss-Kronrod rule family.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GkRule {
    /// 7-point Gauss, 15-point Kronrod.
    G7K15,
    /// 10-point Gauss, 21-point Kronrod.
    G10K21,
}

const XGK15: [f64; 8] = [
    0.991455371120812639206854697526329,
    0.949107912342758524526189684047851,
    0.864864423359769072789712788640926,
    0.741531185599394439863864773280788,
    0.586087235467691130294144838258730,
    0.405845151377397166906606412076961,
    0.207784955007898467600689403773245,
    0.0,
];
const WGK15: [f64; 8] = [
    0.022935322010529224963732008058970,
    0.063092092629978553290700663189204,
    0.104790010322250183839876322541518,
    0.140653259715525918745189590510238,
    0.169004726639267902826583426598550,
    0.190350578064785409913256402421014,
    0.204432940075298892414161999234649,
    0.209482141084727828012999174891714,
];
const WG7: [f64; 4] = [
    0.129484966168869693270611432679082,
    0.279705391489276667901467771423780,
    0.381830050505118944950369775488975,
    0.417959183673469387755102040816327,
];

const XGK21: [f64; 11] = [
    0.995657163025808080735527280689003,
    0.973906528517171720077964012084452,
    0.930157491355708226001207180059508,
    0.865063366688984510732096688423493,
    0.780817726586416897063717578345042,
    0.679409568299024406234327365114874,
    0.562757134668604683339000099272694,
    0.433395394129247190799265943165784,
    0.294392862701460198131126603103866,
    0.148874338981631210884826001129720,
    0.0,
];
const WGK21: [f64; 11] = [
    0.011694638867371874278064396062192,
    0.032558162307964727478818972459390,
    0.054755896574351996031381300244580,
    0.075039674810919952767043140916190,
    0.093125454583697605535065465083366,
    0.109387158802297641899210590325805,
    0.123491976262065851077958109831074,
    0.134709217311473325928054001771707,
    0.142775938577060080797094273138717,
    0.147739104901338491374841515972068,
    0.149445554002916905664936468389821,
];
const WG10: [f64; 5] = [
    0.066671344308688137593568809893332,
    0.149451349150580593145776339657697,
    0.219086362515982043995534934228163,
    0.269266719309996355091226921569469,
    0.295524224714752870173892994651338,
];

/// One Gauss-Kronrod panel on `[a, b]`: returns `(kronrod, error, evals)`.
///
/// The error is the QUADPACK heuristic built from the Gauss/Kronrod
/// difference and the integral of `|f - mean|`.
pub fn gk_panel<F: Fn(f64) -> f64>(
    f: &F,
    a: f64,
    b: f64,
    rule: GkRule,
) -> (f64, f64, usize) {
    let (xgk, wgk): (&[f64], &[f64]) = match rule {
        GkRule::G7K15 => (&XGK15, &WGK15),
        GkRule::G10K21 => (&XGK21, &WGK21),
    };
    let n = xgk.len();
    let c = 0.5 * (a + b);
    let h = 0.5 * (b - a);
    let fc = f(c);
    let mut resk = wgk[n - 1] * fc;
    let mut resg = match rule {
        GkRule::G7K15 => WG7[3] * fc,
        GkRule::G10K21 => 0.0,
    };
    let mut resabs = wgk[n - 1] * fc.abs();
    let mut fv = vec![(0.0, 0.0); n - 1];
    for j in 0..n - 1 {
        let dx = h * xgk[j];
        let f1 = f(c - dx);
        let f2 = f(c + dx);
        fv[j] = (f1, f2);
        resk += wgk[j] * (f1 + f2);
        resabs += wgk[j] * (f1.abs() + f2.abs());
        if j % 2 == 1 {
            let wg = match rule {
                GkRule::G7K15 => WG7[j / 2],
                GkRule::G10K21 => WG10[j / 2],
            };
            resg += wg * (f1 + f2);
        }
    }
    let reskh = resk * 0.5;
    let mut resasc = wgk[n - 1] * (fc - reskh).abs();
    for j in 0..n - 1 {
        resasc += wgk[j] * ((fv[j].0 - reskh).abs() + (fv[j].1 - reskh).abs());
    }
    let result = resk * h;
    resabs *= h.abs();
    resasc *= h.abs();
    let mut err = ((resk - resg) * h).abs();
    if resasc != 0.0 && err != 0.0 {
        err = resasc * (1.0_f64).min((200.0 * err / resasc).powf(1.5));
    }
    if resabs > f64::MIN_POSITIVE / (50.0 * f64::EPSILON) {
        err = err.max(50.0 * f64::EPSILON * resabs);
    }
    (result, err, 2 * n - 1)
}

/// Adaptive Gauss-Kronrod integration of `f` over `[a, b]`.
///
/// The interval with the largest error estimate is bisected until the total
/// error is below `max(abs_tol, rel_tol * |value|)` or `max_panels` panels
/// exist. When plain bisection stagnates (typical of endpoint
/// singularities) the sequence of partial results is accelerated with
/// Wynn's epsilon algorithm (as in QUADPACK's QAGS) and the extrapolated
/// value is returned if its error estimate is smaller.
pub fn integrate_gk<F: Fn(f64) -> f64>(
    f: F,
    a: f64,
    b: f64,
    rule: GkRule,
    abs_tol: f64,
    rel_tol: f64,
    max_panels: usize,
) -> Integral {
    if a == b {
        return Integral { value: 0.0, error: 0.0, evaluations: 0, converged: true };
    }
    let (r0, e0, n0) = gk_panel(&f, a, b, rule);
    let mut panels = vec![(a, b, r0, e0)];
    let mut evals = n0;
    let mut total = r0;
    let mut total_err = e0;
    let mut partials = vec![r0];
    let tol = |v: f64| abs_tol.max(rel_tol * v.abs());
    while total_err > tol(total) && panels.len() < max_panels {
        let mut wi = 0;
        for (i, p) in panels.iter().enumerate() {
            if p.3 > panels[wi].3 {
                wi = i;
            }
        }
        let (pa, pb, _, _) = panels[wi];
        let mid = 0.5 * (pa + pb);
        if mid <= pa.min(pb) || mid >= pa.max(pb) {
            break;
        }
        let (r1, e1, n1) = gk_panel(&f, pa, mid, rule);
        let (r2, e2, n2) = gk_panel(&f, mid, pb, rule);
        evals += n1 + n2;
        panels[wi] = (pa, mid, r1, e1);
        panels.push((mid, pb, r2, e2));
        total = panels.iter().map(|p| p.2).sum();
        total_err = panels.iter().map(|p| p.3).sum();
        partials.push(total);
    }
    let converged = total_err <= tol(total);
    let mut out = Integral { value: total, error: total_err, evaluations: evals, converged };
    if !converged && partials.len() > 8 {
        let eps = crate::kernels::convergence::wynn_epsilon(&partials);
        if let Some(&last) = eps.last() {
            let est = (last - total).abs();
            if last.is_finite() && est < total_err {
                out.value = last;
                out.error = est.max(f64::EPSILON * last.abs());
                out.converged = out.error <= tol(last);
            }
        }
    }
    out
}

/// Adaptive G7K15 integration with a default panel budget (2000).
/// `tol` is the relative tolerance (absolute is `1e-3 * tol`).
pub fn integrate<F: Fn(f64) -> f64>(f: F, a: f64, b: f64, tol: f64) -> Integral {
    integrate_gk(f, a, b, GkRule::G7K15, tol * 1e-3, tol, 2000)
}

/// Integral over `[a, inf)` using `x = a + t/(1-t)`, `t in [0, 1)`.
pub fn integrate_semi_infinite<F: Fn(f64) -> f64>(f: F, a: f64, tol: f64) -> Integral {
    let g = |t: f64| {
        let s = 1.0 - t;
        if s <= 0.0 {
            return 0.0;
        }
        let v = f(a + t / s) / (s * s);
        if v.is_finite() { v } else { 0.0 }
    };
    integrate_gk(g, 0.0, 1.0, GkRule::G10K21, tol * 1e-3, tol, 2000)
}

/// Integral over `(-inf, inf)` using `x = t/(1-t^2)`, `t in (-1, 1)`.
pub fn integrate_infinite<F: Fn(f64) -> f64>(f: F, tol: f64) -> Integral {
    let g = |t: f64| {
        let s = 1.0 - t * t;
        if s <= 0.0 {
            return 0.0;
        }
        let v = f(t / s) * (1.0 + t * t) / (s * s);
        if v.is_finite() { v } else { 0.0 }
    };
    integrate_gk(g, -1.0, 1.0, GkRule::G10K21, tol * 1e-3, tol, 2000)
}

/// Generic double-exponential driver: `map(t) = (x, weight)`.
///
/// Evaluates trapezoid sums with step `2^-level` on `t in [-tmax, tmax]`,
/// reusing earlier levels, until two successive levels agree.
fn de_driver<F: Fn(f64) -> f64, M: Fn(f64) -> (f64, f64)>(
    f: &F,
    map: &M,
    tmax: f64,
    tol: f64,
    max_level: u32,
) -> Integral {
    let term = |t: f64| {
        let (x, w) = map(t);
        if !w.is_finite() || w == 0.0 {
            return 0.0;
        }
        let v = f(x) * w;
        if v.is_finite() { v } else { 0.0 }
    };
    let mut sum = term(0.0);
    let mut evals = 1;
    let mut h = 1.0;
    let mut k = 1.0;
    while k <= tmax {
        sum += term(k) + term(-k);
        evals += 2;
        k += 1.0;
    }
    let mut prev = sum * h;
    let mut result = Integral {
        value: prev,
        error: f64::INFINITY,
        evaluations: evals,
        converged: false,
    };
    for _ in 1..=max_level {
        h *= 0.5;
        let mut k = h;
        while k <= tmax {
            sum += term(k) + term(-k);
            evals += 2;
            k += 2.0 * h;
        }
        let cur = sum * h;
        let err = (cur - prev).abs();
        result = Integral {
            value: cur,
            error: err,
            evaluations: evals,
            converged: err <= tol * cur.abs() || err <= 1e-15,
        };
        if result.converged {
            break;
        }
        prev = cur;
    }
    result
}

/// Tanh-sinh quadrature on `[a, b]`.
///
/// Handles integrable endpoint singularities (`1/sqrt(x)`, `log x`, ...)
/// because the nodes cluster double-exponentially at the end points and the
/// integrand is never evaluated at them. `tol` is the relative tolerance.
pub fn tanh_sinh<F: Fn(f64) -> f64>(f: F, a: f64, b: f64, tol: f64) -> Integral {
    let half = 0.5 * (b - a);
    let map = |t: f64| {
        let u = 0.5 * PI * t.sinh();
        let ch = u.cosh();
        let w = 0.5 * PI * t.cosh() / (ch * ch) * half;
        // distance of the node from its end point, in units of `half`
        let d = 2.0 / ((2.0 * u.abs()).exp() + 1.0);
        if d <= 0.0 || !w.is_finite() {
            return (a, 0.0);
        }
        let x = if t >= 0.0 { b - half * d } else { a + half * d };
        (x, w)
    };
    de_driver(&f, &map, 6.5, tol, 10)
}

/// Exp-sinh quadrature on `[a, inf)`; suitable for algebraic or
/// exponential decay and endpoint singularities at `a`.
pub fn exp_sinh<F: Fn(f64) -> f64>(f: F, a: f64, tol: f64) -> Integral {
    let map = |t: f64| {
        let u = 0.5 * PI * t.sinh();
        let e = u.exp();
        (a + e, 0.5 * PI * t.cosh() * e)
    };
    de_driver(&f, &map, 6.5, tol, 10)
}

/// Sinh-sinh quadrature on `(-inf, inf)` for rapidly decaying integrands.
pub fn sinh_sinh<F: Fn(f64) -> f64>(f: F, tol: f64) -> Integral {
    let map = |t: f64| {
        let u = 0.5 * PI * t.sinh();
        (u.sinh(), 0.5 * PI * t.cosh() * u.cosh())
    };
    de_driver(&f, &map, 6.5, tol, 10)
}

/// Clenshaw-Curtis nodes and weights on `[-1, 1]` with `n + 1` points
/// (`n` is rounded up to an even number, at least 2).
pub fn clenshaw_curtis_rule(n: usize) -> (Vec<f64>, Vec<f64>) {
    let n = n.max(2).div_ceil(2) * 2;
    let nf = n as f64;
    let mut x = vec![0.0; n + 1];
    let mut w = vec![0.0; n + 1];
    for k in 0..=n {
        let th = PI * k as f64 / nf;
        x[k] = th.cos();
        let mut s = 0.0;
        for j in 1..=n / 2 {
            let b = if j == n / 2 { 1.0 } else { 2.0 };
            let jf = j as f64;
            s += b * (2.0 * jf * th).cos() / (4.0 * jf * jf - 1.0);
        }
        let c = if k == 0 || k == n { 1.0 } else { 2.0 };
        w[k] = c / nf * (1.0 - s);
    }
    (x, w)
}

/// Clenshaw-Curtis quadrature of `f` on `[a, b]` with `n + 1` points.
pub fn clenshaw_curtis<F: Fn(f64) -> f64>(f: F, a: f64, b: f64, n: usize) -> f64 {
    let (x, w) = clenshaw_curtis_rule(n);
    let h = 0.5 * (b - a);
    let m = 0.5 * (a + b);
    h * x.iter().zip(&w).map(|(&xi, &wi)| wi * f(m + h * xi)).sum::<f64>()
}

/// Romberg integration with Richardson extrapolation of the trapezoid rule.
/// The error is the last difference of diagonal entries.
pub fn romberg<F: Fn(f64) -> f64>(
    f: F,
    a: f64,
    b: f64,
    tol: f64,
    max_levels: usize,
) -> Integral {
    let m = max_levels.max(2);
    let mut r = vec![vec![0.0; m]; m];
    let mut h = b - a;
    r[0][0] = 0.5 * h * (f(a) + f(b));
    let mut evals = 2;
    let mut err = f64::INFINITY;
    let mut last = 0;
    for i in 1..m {
        h *= 0.5;
        let mut s = 0.0;
        let cnt = 1usize << (i - 1);
        for k in 0..cnt {
            s += f(a + (2 * k + 1) as f64 * h);
        }
        evals += cnt;
        r[i][0] = 0.5 * r[i - 1][0] + h * s;
        for j in 1..=i {
            let p = 4f64.powi(j as i32);
            r[i][j] = (p * r[i][j - 1] - r[i - 1][j - 1]) / (p - 1.0);
        }
        err = (r[i][i] - r[i - 1][i - 1]).abs();
        last = i;
        if err <= tol * r[i][i].abs().max(1e-300) && i >= 3 {
            break;
        }
    }
    let v = r[last][last];
    Integral {
        value: v,
        error: err,
        evaluations: evals,
        converged: err <= tol * v.abs().max(1e-300),
    }
}

/// Filon quadrature for `integral_a^b f(x) cos(w x) dx` (`cosine = true`)
/// or `sin(w x)` (`cosine = false`), using `2 * panels` sub-intervals.
///
/// Exact when `f` is quadratic on each double panel, so the accuracy does
/// not degrade as the frequency `w` grows.
pub fn filon<F: Fn(f64) -> f64>(
    f: F,
    a: f64,
    b: f64,
    omega: f64,
    cosine: bool,
    panels: usize,
) -> f64 {
    let n = panels.max(1);
    let h = (b - a) / (2 * n) as f64;
    let th = omega * h;
    let (alpha, beta, gamma) = if th.abs() < 0.1 {
        let t2 = th * th;
        (
            2.0 * th.powi(3) / 45.0 - 2.0 * th.powi(5) / 315.0 + 2.0 * th.powi(7) / 4725.0,
            2.0 / 3.0 + 2.0 * t2 / 15.0 + 4.0 * t2 * t2 / 105.0 + 2.0 * t2.powi(3) / 567.0,
            4.0 / 3.0 - 2.0 * t2 / 15.0 + t2 * t2 / 210.0 - t2.powi(3) / 11340.0,
        )
    } else {
        let (s, c) = th.sin_cos();
        let t3 = th.powi(3);
        (
            (th * th + th * s * c - 2.0 * s * s) / t3,
            2.0 * (th * (1.0 + c * c) - 2.0 * s * c) / t3,
            4.0 * (s - th * c) / t3,
        )
    };
    let g = |x: f64| if cosine { (omega * x).cos() } else { (omega * x).sin() };
    let (fa, fb) = (f(a), f(b));
    let mut even = 0.0;
    let mut odd = 0.0;
    for i in 0..=2 * n {
        let x = a + i as f64 * h;
        let v = f(x) * g(x);
        if i % 2 == 0 {
            let wgt = if i == 0 || i == 2 * n { 0.5 } else { 1.0 };
            even += wgt * v;
        } else {
            odd += v;
        }
    }
    let boundary = if cosine {
        fb * (omega * b).sin() - fa * (omega * a).sin()
    } else {
        -(fb * (omega * b).cos() - fa * (omega * a).cos())
    };
    h * (alpha * boundary + beta * even + gamma * odd)
}

/// Radical inverse of `i` in base `base`.
fn radical_inverse(mut i: u64, base: u64) -> f64 {
    let inv = 1.0 / base as f64;
    let mut f = inv;
    let mut r = 0.0;
    while i > 0 {
        r += (i % base) as f64 * f;
        i /= base;
        f *= inv;
    }
    r
}

const PRIMES: [u64; 16] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53];

/// The `index`-th point (`index >= 1`) of the Halton sequence in `dim <= 16` dimensions.
pub fn halton_point(index: u64, dim: usize) -> Vec<f64> {
    (0..dim.min(16)).map(|d| radical_inverse(index, PRIMES[d])).collect()
}

/// Generator of the Sobol sequence (up to 8 dimensions, Joe-Kuo directions).
#[derive(Debug, Clone)]
pub struct Sobol {
    dim: usize,
    dirs: Vec<[u32; 32]>,
    state: Vec<u32>,
    index: u32,
}

impl Sobol {
    /// Creates a generator for `dim` (clamped to 1..=8) dimensions.
    #[must_use]
    pub fn new(dim: usize) -> Self {
        // (degree s, polynomial coefficients a, initial m)
        const PARAMS: [(usize, u32, [u32; 5]); 7] = [
            (1, 0, [1, 0, 0, 0, 0]),
            (2, 1, [1, 3, 0, 0, 0]),
            (3, 1, [1, 3, 1, 0, 0]),
            (3, 2, [1, 1, 1, 0, 0]),
            (4, 1, [1, 1, 3, 3, 0]),
            (4, 4, [1, 3, 5, 13, 0]),
            (5, 2, [1, 1, 5, 5, 17]),
        ];
        let dim = dim.clamp(1, 8);
        let mut dirs = Vec::with_capacity(dim);
        let mut d0 = [0u32; 32];
        for (i, v) in d0.iter_mut().enumerate() {
            *v = 1u32 << (31 - i);
        }
        dirs.push(d0);
        for d in 1..dim {
            let (s, a, m0) = PARAMS[d - 1];
            let mut m = [0u32; 32];
            for i in 0..32 {
                if i < s {
                    m[i] = m0[i];
                } else {
                    let mut v = m[i - s] ^ (m[i - s] << s);
                    for k in 1..s {
                        if (a >> (s - 1 - k)) & 1 == 1 {
                            v ^= m[i - k] << k;
                        }
                    }
                    m[i] = v;
                }
            }
            let mut dv = [0u32; 32];
            for i in 0..32 {
                dv[i] = m[i] << (31 - i);
            }
            dirs.push(dv);
        }
        Self { dim, dirs, state: vec![0; dim], index: 0 }
    }

    /// Next point in `[0, 1)^dim` (the origin is skipped).
    pub fn next_point(&mut self) -> Vec<f64> {
        let c = (!self.index).trailing_zeros() as usize;
        self.index = self.index.wrapping_add(1);
        for d in 0..self.dim {
            self.state[d] ^= self.dirs[d][c.min(31)];
        }
        self.state.iter().map(|&s| f64::from(s) / 4_294_967_296.0).collect()
    }
}

/// Low-discrepancy sequence selector for [`quasi_monte_carlo`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LowDiscrepancy {
    /// Halton sequence (dimension up to 16).
    Halton,
    /// Sobol sequence (dimension up to 8).
    Sobol,
}

/// Quasi-Monte Carlo integral of `f` over the box `[lo, hi]`.
///
/// The error estimate is the difference to the estimate that uses only the
/// first half of the points.
pub fn quasi_monte_carlo<F: Fn(&[f64]) -> f64>(
    f: F,
    lo: &[f64],
    hi: &[f64],
    n: usize,
    seq: LowDiscrepancy,
) -> Integral {
    let dim = lo.len();
    let vol: f64 = lo.iter().zip(hi).map(|(a, b)| b - a).product();
    let mut sob = Sobol::new(dim);
    let mut sum = 0.0;
    let mut half = 0.0;
    let mut x = vec![0.0; dim];
    for i in 1..=n {
        let p = match seq {
            LowDiscrepancy::Halton => halton_point(i as u64, dim),
            LowDiscrepancy::Sobol => sob.next_point(),
        };
        for d in 0..dim {
            x[d] = lo[d] + (hi[d] - lo[d]) * p[d];
        }
        sum += f(&x);
        if i == n / 2 {
            half = sum;
        }
    }
    let value = vol * sum / n as f64;
    let vh = vol * half / (n / 2).max(1) as f64;
    Integral { value, error: (value - vh).abs(), evaluations: n, converged: true }
}

/// Plain Monte Carlo integral over `[lo, hi]` with a seeded generator;
/// the error is the standard error of the mean.
pub fn monte_carlo<F: Fn(&[f64]) -> f64>(
    f: F,
    lo: &[f64],
    hi: &[f64],
    n: usize,
    seed: u64,
) -> Integral {
    let dim = lo.len();
    let vol: f64 = lo.iter().zip(hi).map(|(a, b)| b - a).product();
    let mut rng = Rng::new(seed);
    let mut x = vec![0.0; dim];
    let (mut s, mut s2) = (0.0, 0.0);
    for _ in 0..n {
        for d in 0..dim {
            x[d] = lo[d] + (hi[d] - lo[d]) * rng.uniform();
        }
        let v = f(&x);
        s += v;
        s2 += v * v;
    }
    let nf = n as f64;
    let mean = s / nf;
    let var = (s2 / nf - mean.powi(2)).max(0.0);
    Integral { value: vol * mean, error: vol * (var / nf).sqrt(), evaluations: n, converged: true }
}

/// Nested adaptive Gauss-Kronrod cubature over the box `[lo, hi]`
/// (practical up to about four dimensions).
pub fn cubature_nested<F: Fn(&[f64]) -> f64>(
    f: F,
    lo: &[f64],
    hi: &[f64],
    tol: f64,
) -> Integral {
    fn rec<F: Fn(&[f64]) -> f64>(
        f: &F,
        lo: &[f64],
        hi: &[f64],
        level: usize,
        point: &RefCell<Vec<f64>>,
        count: &RefCell<usize>,
        tol: f64,
    ) -> f64 {
        if level == lo.len() {
            *count.borrow_mut() += 1;
            return f(&point.borrow());
        }
        let g = |x: f64| {
            point.borrow_mut()[level] = x;
            rec(f, lo, hi, level + 1, point, count, tol)
        };
        integrate_gk(g, lo[level], hi[level], GkRule::G7K15, tol * 1e-2, tol, 200).value
    }
    let point = RefCell::new(vec![0.0; lo.len()]);
    let count = RefCell::new(0usize);
    let value = rec(&f, lo, hi, 0, &point, &count, tol);
    let e = *count.borrow();
    Integral { value, error: tol * value.abs(), evaluations: e, converged: true }
}

/// A region of the Genz-Malik adaptive cubature.
struct GmRegion {
    center: Vec<f64>,
    half: Vec<f64>,
    value: f64,
    error: f64,
    axis: usize,
}

impl PartialEq for GmRegion {
    fn eq(&self, other: &Self) -> bool {
        self.error == other.error
    }
}
impl Eq for GmRegion {}
impl PartialOrd for GmRegion {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for GmRegion {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.error.partial_cmp(&other.error).unwrap_or(std::cmp::Ordering::Equal)
    }
}

/// Applies the embedded Genz-Malik degree-7/5 rule on one region and
/// returns the region with its value, error estimate and split axis.
fn gm_apply<F: Fn(&[f64]) -> f64>(f: &F, center: Vec<f64>, half: Vec<f64>) -> GmRegion {
    let n = center.len();
    let nf = n as f64;
    let l2 = (9.0_f64 / 70.0).sqrt();
    let l4 = (9.0_f64 / 10.0).sqrt();
    let l5 = (9.0_f64 / 19.0).sqrt();
    let w1 = (12824.0 - 9120.0 * nf + 400.0 * nf * nf) / 19683.0;
    let w2 = 980.0 / 6561.0;
    let w3 = (1820.0 - 400.0 * nf) / 19683.0;
    let w4 = 200.0 / 19683.0;
    let w5 = 6859.0 / 19683.0 / 2.0_f64.powi(n as i32);
    let e1 = (729.0 - 950.0 * nf + 50.0 * nf * nf) / 729.0;
    let e2 = 245.0 / 486.0;
    let e3 = (265.0 - 100.0 * nf) / 1458.0;
    let e4 = 25.0 / 729.0;
    let mut p = center.clone();
    let f0 = f(&p);
    let mut s2 = vec![0.0; n];
    let mut s3 = vec![0.0; n];
    for i in 0..n {
        let c = center[i];
        p[i] = c + l2 * half[i];
        let a = f(&p);
        p[i] = c - l2 * half[i];
        s2[i] = a + f(&p);
        p[i] = c + l4 * half[i];
        let a = f(&p);
        p[i] = c - l4 * half[i];
        s3[i] = a + f(&p);
        p[i] = c;
    }
    let mut s4 = 0.0;
    for i in 0..n {
        for j in i + 1..n {
            for (si, sj) in [(1.0, 1.0), (1.0, -1.0), (-1.0, 1.0), (-1.0, -1.0)] {
                p[i] = center[i] + si * l4 * half[i];
                p[j] = center[j] + sj * l4 * half[j];
                s4 += f(&p);
            }
            p[i] = center[i];
            p[j] = center[j];
        }
    }
    let mut s5 = 0.0;
    for mask in 0..(1usize << n) {
        for i in 0..n {
            let sgn = if (mask >> i) & 1 == 1 { 1.0 } else { -1.0 };
            p[i] = center[i] + sgn * l5 * half[i];
        }
        s5 += f(&p);
    }
    let vol: f64 = half.iter().map(|h| 2.0 * h).product();
    let t2: f64 = s2.iter().sum();
    let t3: f64 = s3.iter().sum();
    let i7 = vol * (w1 * f0 + w2 * t2 + w3 * t3 + w4 * s4 + w5 * s5);
    let i5 = vol * (e1 * f0 + e2 * t2 + e3 * t3 + e4 * s4);
    let ratio = (l4 * l4) / (l2 * l2);
    let mut axis = 0;
    let mut best = -1.0;
    for i in 0..n {
        let d = (s3[i] - 2.0 * f0 - ratio * (s2[i] - 2.0 * f0)).abs();
        if d > best * (1.0 + 1e-12) || (d >= best * (1.0 - 1e-12) && half[i] > half[axis]) {
            best = d;
            axis = i;
        }
    }
    GmRegion { center, half, value: i7, error: (i7 - i5).abs(), axis }
}

/// Number of integrand evaluations per Genz-Malik rule application.
fn gm_cost(n: usize) -> usize {
    1 + 4 * n + 2 * n * (n - 1) + (1usize << n)
}

/// Adaptive Genz-Malik cubature over the box `[lo, hi]`.
///
/// Each region is integrated with the degree-7 Genz-Malik rule
/// (`1 + 4n + 2n(n-1) + 2^n` points) whose embedded degree-5 rule gives
/// the error estimate; the region with the largest error is bisected
/// along the axis with the largest fourth divided difference. The box
/// is first split uniformly into a few sub-boxes so that narrow peaks
/// are not missed. Stops when the total error is below
/// `max(abs_tol, rel_tol * |value|)` or after `max_evals` integrand
/// evaluations (`converged == false`). Intended for 2 to about 8
/// dimensions.
pub fn cubature_genz_malik<F: Fn(&[f64]) -> f64>(
    f: F,
    lo: &[f64],
    hi: &[f64],
    abs_tol: f64,
    rel_tol: f64,
    max_evals: usize,
) -> Integral {
    let n = lo.len();
    if n == 0 || hi.len() != n || n > 20 {
        return Integral { value: f64::NAN, error: f64::NAN, evaluations: 0, converged: false };
    }
    let cost = gm_cost(n);
    let mut evals = 0usize;
    // initial uniform split
    let parts = ((64.0_f64).powf(1.0 / n as f64).floor() as usize).max(2);
    let mut heap = std::collections::BinaryHeap::new();
    let count = parts.pow(n as u32);
    for mut idx in 0..count {
        let mut c = vec![0.0; n];
        let mut h = vec![0.0; n];
        for i in 0..n {
            let k = idx % parts;
            idx /= parts;
            let w = (hi[i] - lo[i]) / parts as f64;
            c[i] = lo[i] + w * (k as f64 + 0.5);
            h[i] = 0.5 * w;
        }
        heap.push(gm_apply(&f, c, h));
        evals += cost;
    }
    loop {
        let value: f64 = heap.iter().map(|r| r.value).sum();
        let error: f64 = heap.iter().map(|r| r.error).sum();
        let target = abs_tol.max(rel_tol * value.abs());
        if error <= target {
            return Integral { value, error, evaluations: evals, converged: true };
        }
        if evals + 2 * cost > max_evals || !value.is_finite() {
            return Integral { value, error, evaluations: evals, converged: false };
        }
        let Some(r) = heap.pop() else {
            return Integral { value, error, evaluations: evals, converged: false };
        };
        let mut h = r.half.clone();
        h[r.axis] *= 0.5;
        let mut c1 = r.center.clone();
        let mut c2 = r.center.clone();
        c1[r.axis] -= h[r.axis];
        c2[r.axis] += h[r.axis];
        heap.push(gm_apply(&f, c1, h.clone()));
        heap.push(gm_apply(&f, c2, h));
        evals += 2 * cost;
    }
}
