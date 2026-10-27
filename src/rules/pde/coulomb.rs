//! Coulomb (hydrogen-like) radial problems and spheroidal eigenvalues.
//!
//! * **Coulomb endpoints.** `a u_t = b [u_rr + (e₀ + 2Z/r - ℓ(ℓ+1)/r²) u]`
//!   (also `a u_tt`) on the half-line `r > 0`: the operator
//!   `-∂² + ℓ(ℓ+1)/r² - 2Z/r` has the bound states
//!   `u_{nℓ}(r) = N r^{ℓ+1} e^{-Zr/n} L_{n-ℓ-1}^{2ℓ+1}(2Zr/n)` with
//!   `E_n = -Z²/n²` (exact, Whittaker/Laguerre type); the data are projected
//!   on the bound states by quadrature and the answer is returned when the
//!   projection reproduces the data (data in the span of the bound states,
//!   `n ≤ 8`); the continuum (scattering) states are *not* part of the
//!   expansion, so other data are refused.
//! * **Operators.** `coulomb_eigen(Z, ℓ, n)`, `coulomb_mode(Z, ℓ, n, r)`
//!   (orthonormal on `(0, ∞)`), and the spheroidal angular functions
//!   `spheroidal_eigen(m, n, γ)`, `spheroidal_mode(m, n, γ, x)` for
//!   `(1-x²) y'' - 2 x y' + (λ - γ x² - m²/(1-x²)) y = 0` on `[-1, 1]`
//!   (`γ = c² > 0` prolate, `< 0` oblate), computed in the orthonormal
//!   associated Legendre basis where the operator is `diag(j(j+1)) + γ X²`
//!   (Jacobi eigenvalue iteration, 60 basis functions).

use super::Conditions;
use super::Problem;
use super::hill::jacobi;
use super::modes::exact_value;
use super::spectral::evolution_time;
use super::util::call;
use super::util::defint;
use super::util::div;
use super::util::is_zero_number;
use super::util::sqrt;
use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Env;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpDescriptor;
use crate::graph::RuleError;
use crate::graph::rule::Installer;
use crate::rules::complex::build::add;
use crate::rules::complex::build::mul;
use crate::rules::complex::build::neg;
use crate::rules::complex::build::powi;
use crate::rules::complex::build::sub;
use std::cell::RefCell;
use std::rc::Rc;

/// The largest principal quantum number used in an expansion.
const MAX_STATE: usize = 14;
/// Basis size of the spheroidal matrix.
const BASIS: usize = 60;

fn factorial(n: usize) -> f64 {
    (1..=n).fold(1.0, |a, k| a * k as f64)
}

#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss, clippy::float_cmp)]
fn whole(v: f64) -> Option<usize> {
    ((0.0..=1e4).contains(&v) && v.fract() == 0.0).then_some(v as usize)
}

/// Generalised Laguerre polynomial `L_k^{α}(ρ)` by the three-term recurrence.
fn laguerre(
    k: usize,
    alpha: f64,
    rho: f64,
) -> f64 {
    let (mut previous, mut current) = (1.0, 1.0 + alpha - rho);
    if k == 0 {
        return previous;
    }
    for j in 1..k {
        let jf = j as f64;
        let next = ((2.0 * jf + 1.0 + alpha - rho) * current - (jf + alpha) * previous) / (jf + 1.0);
        previous = current;
        current = next;
    }
    current
}

/// The normalised radial function `u_{nℓ}(r) = r R_{nℓ}(r)` for charge `z`.
pub(super) fn coulomb_value(
    z: f64,
    l: usize,
    n: usize,
    r: f64,
) -> f64 {
    if n <= l {
        return f64::NAN;
    }
    let k = n - l - 1;
    let rho = 2.0 * z * r / n as f64;
    let norm = ((2.0 * z / n as f64).powi(3) * factorial(k) / (2.0 * n as f64 * factorial(n + l))).sqrt();
    r * norm * (-0.5 * rho).exp() * rho.powi(i32::try_from(l).unwrap_or(0)) * laguerre(k, 2.0 * l as f64 + 1.0, rho)
}

/// `coulomb_eigen(Z, ℓ, n)`: `-Z²/n²`.
pub(super) fn coulomb_eigen(args: &[f64]) -> f64 {
    let (Some(&z), Some(l), Some(n)) = (args.first(), args.get(1).and_then(|&v| whole(v)), args.get(2).and_then(|&v| whole(v))) else {
        return f64::NAN;
    };
    if n <= l || n == 0 { f64::NAN } else { -z * z / (n * n) as f64 }
}

/// `coulomb_mode(Z, ℓ, n, r)`.
pub(super) fn coulomb_mode(args: &[f64]) -> f64 {
    let (Some(&z), Some(l), Some(n), Some(&r)) = (args.first(), args.get(1).and_then(|&v| whole(v)), args.get(2).and_then(|&v| whole(v)), args.get(3)) else {
        return f64::NAN;
    };
    coulomb_value(z, l, n, r)
}

/// Eigenpairs of the spheroidal angular operator in the orthonormal
/// associated Legendre basis `P̃_{m+k}^m`.
fn spheroidal_uncached(
    m: usize,
    gamma: f64,
) -> (Vec<f64>, Vec<Vec<f64>>) {
    let a = |j: usize| -> f64 {
        let (jf, mf) = (j as f64, m as f64);
        if j <= m { 0.0 } else { ((jf * jf - mf * mf) / (4.0 * jf * jf - 1.0)).sqrt() }
    };
    let mut h = vec![vec![0.0; BASIS]; BASIS];
    for k in 0..BASIS {
        let j = m + k;
        h[k][k] = (j * (j + 1)) as f64 + gamma * (a(j + 1).powi(2) + a(j).powi(2));
        if k + 2 < BASIS {
            let v = gamma * a(j + 1) * a(j + 2);
            h[k][k + 2] = v;
            h[k + 2][k] = v;
        }
    }
    jacobi(h)
}

type SpheroidalPair = Rc<(Vec<f64>, Vec<Vec<f64>>)>;

thread_local! {
    static SPHEROIDAL: RefCell<Vec<((usize, u64), SpheroidalPair)>> = const { RefCell::new(Vec::new()) };
}

/// The cached eigenpairs `(values, vectors)` of order `m`.
pub(super) fn spheroidal_pair(
    m: usize,
    gamma: f64,
) -> SpheroidalPair {
    let key = (m, gamma.to_bits());
    if let Some(found) = SPHEROIDAL.with(|c| c.borrow().iter().find(|(k, _)| *k == key).map(|(_, v)| Rc::clone(v))) {
        return found;
    }
    let value = Rc::new(spheroidal_uncached(m, gamma));
    SPHEROIDAL.with(|c| {
        let mut cache = c.borrow_mut();
        if cache.len() >= 64 {
            cache.remove(0);
        }
        cache.push((key, Rc::clone(&value)));
    });
    value
}

/// The orthonormal angular function with expansion coefficients `d` (times
/// `sign`) at `x`.
pub(super) fn spheroidal_mode_value(
    m: usize,
    gamma: f64,
    x: f64,
    d: &[f64],
    sign: f64,
) -> f64 {
    let _ = gamma;
    // Orthonormal P̃_{m+k}^m(x) by the recurrence x P̃_j = a_{j+1} P̃_{j+1} + a_j P̃_{j-1}.
    let a = |j: usize| -> f64 {
        let (jf, mf) = (j as f64, m as f64);
        if j <= m { 0.0 } else { ((jf * jf - mf * mf) / (4.0 * jf * jf - 1.0)).sqrt() }
    };
    let mut start = (1.0 - x * x).max(0.0).powf(m as f64 / 2.0);
    let mut scale = 0.5_f64.sqrt();
    for i in 1..=m {
        scale *= (((2 * i + 1) as f64) / ((2 * i) as f64)).sqrt();
    }
    start *= scale;
    let mut total = d[0] * start;
    let (mut previous, mut current) = (0.0, start);
    for k in 1..d.len() {
        let j = m + k - 1;
        let next = (x * current - a(j) * previous) / a(j + 1);
        previous = current;
        current = next;
        total += d[k] * current;
    }
    sign * total
}

/// `spheroidal_eigen(m, n, γ)`, `n ≥ m`.
pub(super) fn spheroidal_eigen(args: &[f64]) -> f64 {
    let (Some(m), Some(n), Some(&gamma)) = (args.first().and_then(|&v| whole(v)), args.get(1).and_then(|&v| whole(v)), args.get(2)) else {
        return f64::NAN;
    };
    if n < m || n - m >= BASIS / 2 {
        return f64::NAN;
    }
    spheroidal_pair(m, gamma).0.get(n - m).copied().unwrap_or(f64::NAN)
}

/// `spheroidal_mode(m, n, γ, x)`: the angular function, normalised on
/// `[-1, 1]`, positive leading coefficient.
pub(super) fn spheroidal_mode(args: &[f64]) -> f64 {
    let (Some(m), Some(n), Some(&gamma), Some(&x)) = (args.first().and_then(|&v| whole(v)), args.get(1).and_then(|&v| whole(v)), args.get(2), args.get(3)) else {
        return f64::NAN;
    };
    if n < m || n - m >= BASIS / 2 || x.abs() > 1.0 {
        return f64::NAN;
    }
    let pair = spheroidal_pair(m, gamma);
    let Some(d) = pair.1.get(n - m) else {
        return f64::NAN;
    };
    let sign = if d[n - m] < 0.0 { -1.0 } else { 1.0 };
    spheroidal_mode_value(m, gamma, x, d, sign)
}

/// Registers the operators.
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    i.op(OpDescriptor::new("coulomb_eigen", Arity::Fixed(3)).eval(coulomb_eigen))?;
    i.op(OpDescriptor::new("coulomb_mode", Arity::Fixed(4)).eval(coulomb_mode))?;
    i.op(OpDescriptor::new("coulomb_wave", Arity::Fixed(4)).eval(coulomb_wave))?;
    i.op(OpDescriptor::new("spheroidal_eigen", Arity::Fixed(3)).eval(spheroidal_eigen))?;
    i.op(OpDescriptor::new("spheroidal_mode", Arity::Fixed(4)).eval(spheroidal_mode))?;
    Ok(())
}


// ----------------------------------------------------------------------
// Continuum (scattering) states
// ----------------------------------------------------------------------

/// Grid step and extent of the tabulated scattering solutions.
const WAVE_STEP: f64 = 0.01;
const WAVE_MAX: f64 = 120.0;
/// The upper limit of the continuum integral.
const K_MAX: f64 = 16.0;

type WaveTable = Rc<Vec<f64>>;

thread_local! {
    static WAVES: RefCell<Vec<([u64; 3], WaveTable)>> = const { RefCell::new(Vec::new()) };
}

/// The regular solution of `u'' = (ℓ(ℓ+1)/r² - 2Z/r - k²) u` on the grid
/// `r = i h`, normalised on the energy scale: `u → √(2/π) sin(kr + η ln 2kr + …)`
/// so that `∫ u_k u_{k'} dr = δ(k - k')` (`η = -Z/k`). The amplitude is read
/// off the WKB invariant `p u² + u'²/p`, `p² = k² + 2Z/r - ℓ(ℓ+1)/r²`, at the
/// end of the grid.
fn wave_table(
    z: f64,
    l: usize,
    k: f64,
) -> Option<WaveTable> {
    let key = [z.to_bits(), l as u64, k.to_bits()];
    if let Some(found) = WAVES.with(|w| w.borrow().iter().find(|(c, _)| *c == key).map(|(_, v)| Rc::clone(v))) {
        return Some(found);
    }
    if k <= 0.0 || z < 0.0 {
        return None;
    }
    let steps = (WAVE_MAX / WAVE_STEP) as usize;
    let lf = l as f64;
    let q = |r: f64| lf * (lf + 1.0) / (r * r) - 2.0 * z / r - k * k;
    // Frobenius start: u = Σ a_j r^{j+ℓ+1}.
    let series = |r: f64| -> (f64, f64) {
        let mut a = vec![1.0, -2.0 * z / ((lf + 1.0) * (lf + 2.0) - lf * (lf + 1.0))];
        for j in 2..10 {
            let jf = j as f64;
            let next = (-2.0 * z * a[j - 1] - k * k * a[j - 2]) / (jf * (jf + 2.0 * lf + 1.0));
            a.push(next);
        }
        let (mut u, mut du) = (0.0, 0.0);
        for (j, c) in a.iter().enumerate() {
            let e = j as f64 + lf + 1.0;
            u += c * r.powf(e);
            du += c * e * r.powf(e - 1.0);
        }
        (u, du)
    };
    let start = 5_usize;
    let mut table = vec![0.0; steps + 1];
    for (i, slot) in table.iter_mut().enumerate().take(start) {
        *slot = series(WAVE_STEP * i as f64).0;
    }
    let (mut u, mut p) = series(WAVE_STEP * start as f64);
    table[start] = u;
    let h = WAVE_STEP;
    for i in start..steps {
        let r = WAVE_STEP * i as f64;
        let (k1u, k1p) = (p, q(r) * u);
        let (k2u, k2p) = (p + 0.5 * h * k1p, q(r + 0.5 * h) * (u + 0.5 * h * k1u));
        let (k3u, k3p) = (p + 0.5 * h * k2p, q(r + 0.5 * h) * (u + 0.5 * h * k2u));
        let (k4u, k4p) = (p + h * k3p, q(r + h) * (u + h * k3u));
        u += h / 6.0 * (k1u + 2.0 * k2u + 2.0 * k3u + k4u);
        p += h / 6.0 * (k1p + 2.0 * k2p + 2.0 * k3p + k4p);
        table[i + 1] = u;
    }
    let r_end = WAVE_MAX;
    let local = (k * k + 2.0 * z / r_end - lf * (lf + 1.0) / (r_end * r_end)).max(1e-12).sqrt();
    let invariant = local * u * u + p * p / local;
    let amplitude = (invariant / k).sqrt();
    let scale = (2.0 / std::f64::consts::PI).sqrt() / amplitude;
    for v in &mut table {
        *v *= scale;
    }
    let value = Rc::new(table);
    WAVES.with(|w| {
        let mut cache = w.borrow_mut();
        if cache.len() >= 600 {
            cache.remove(0);
        }
        cache.push((key, Rc::clone(&value)));
    });
    Some(value)
}

/// The scattering solution `u_k(r)` of the Coulomb problem (cubic interpolation).
pub(super) fn coulomb_wave_value(
    z: f64,
    l: usize,
    k: f64,
    r: f64,
) -> f64 {
    let Some(table) = wave_table(z, l, k) else {
        return f64::NAN;
    };
    if !(0.0..=WAVE_MAX - 2.0 * WAVE_STEP).contains(&r) {
        return f64::NAN;
    }
    let position = r / WAVE_STEP;
    let i = (position.floor() as usize).clamp(1, table.len() - 3);
    let u = position - i as f64;
    let (p0, p1, p2, p3) = (table[i - 1], table[i], table[i + 1], table[i + 2]);
    let a = -u * (u - 1.0) * (u - 2.0) / 6.0;
    let b = (u + 1.0) * (u - 1.0) * (u - 2.0) / 2.0;
    let c = -(u + 1.0) * u * (u - 2.0) / 2.0;
    let d = (u + 1.0) * u * (u - 1.0) / 6.0;
    a * p0 + b * p1 + c * p2 + d * p3
}

/// `coulomb_wave(Z, ℓ, k, r)`.
pub(super) fn coulomb_wave(args: &[f64]) -> f64 {
    let (Some(&z), Some(l), Some(&k), Some(&r)) = (args.first(), args.get(1).and_then(|&v| whole(v)), args.get(2), args.get(3)) else {
        return f64::NAN;
    };
    coulomb_wave_value(z, l, k, r)
}

/// `∫_a^b f dr` by composite six-point Gauss–Legendre quadrature with panels of width `width`.
fn integrate_range(
    (a, b): (f64, f64),
    width: f64,
    f: &dyn Fn(f64) -> Option<f64>,
) -> Option<f64> {
    const NODES: [(f64, f64); 6] = [
        (-0.932_469_514_203_152, 0.171_324_492_379_170),
        (-0.661_209_386_466_265, 0.360_761_573_048_139),
        (-0.238_619_186_083_197, 0.467_913_934_572_691),
        (0.238_619_186_083_197, 0.467_913_934_572_691),
        (0.661_209_386_466_265, 0.360_761_573_048_139),
        (0.932_469_514_203_152, 0.171_324_492_379_170),
    ];
    let panels = (((b - a) / width).ceil() as usize).max(1);
    let h = (b - a) / panels as f64;
    let mut total = 0.0;
    for panel in 0..panels {
        let mid = a + (panel as f64 + 0.5) * h;
        for &(x, w) in &NODES {
            total += 0.5 * h * w * f(mid + 0.5 * h * x)?;
        }
    }
    Some(total)
}

/// Checks that bound states plus the scattering integral reproduce `data`
/// (to 3% of its size) and returns the upper radius of the data's support.
fn continuum_support(
    cx: &Cx<'_>,
    data: NodeId,
    symbol: crate::graph::SymbolId,
    (z, l): (f64, usize),
    bound: &[(usize, f64)],
) -> Option<f64> {
    let value = |r: f64| evaluate(cx, data, symbol, r);
    let scale = (1..=80).map(|k| value(0.25 * f64::from(k)).map_or(0.0, f64::abs)).fold(1e-12, f64::max);
    let mut rmax = 8.0;
    for k in 1..=((WAVE_MAX - 4.0) * 2.0) as usize {
        let r = 0.5 * k as f64;
        if value(r)?.abs() > 1e-9 * scale {
            rmax = r + 2.0;
        }
    }
    let rmax = rmax.min(WAVE_MAX - 4.0);
    let coefficient = |k: f64| integrate_range((0.0, rmax), 0.5, &|r| Some(value(r)? * coulomb_wave_value(z, l, k, r)));
    // The reconstruction at t = 0.
    let samples = [0.3, 0.8, 1.5, 2.5, 4.0, 6.0];
    let mut rebuilt = vec![0.0; samples.len()];
    let mut peak = 0.0_f64;
    let mut last = 0.0_f64;
    let panels = (K_MAX / 0.5) as usize;
    for panel in 0..panels {
        let (a, b) = (0.5 * panel as f64, 0.5 * (panel + 1) as f64);
        let mut panel_peak = 0.0_f64;
        // Gauss–Legendre nodes of the panel.
        for &(x, w) in &[
            (-0.932_469_514_203_152, 0.171_324_492_379_170),
            (-0.661_209_386_466_265, 0.360_761_573_048_139),
            (-0.238_619_186_083_197, 0.467_913_934_572_691),
            (0.238_619_186_083_197, 0.467_913_934_572_691),
            (0.661_209_386_466_265, 0.360_761_573_048_139),
            (0.932_469_514_203_152, 0.171_324_492_379_170),
        ] {
            let k = f64::midpoint(a, b) + 0.5 * (b - a) * x;
            let weight = 0.5 * (b - a) * w;
            let c = coefficient(k)?;
            panel_peak = panel_peak.max(c.abs());
            for (slot, &r) in rebuilt.iter_mut().zip(&samples) {
                *slot += weight * c * coulomb_wave_value(z, l, k, r);
            }
        }
        peak = peak.max(panel_peak);
        last = panel_peak;
    }
    if last > 5e-2 * peak.max(1e-12) {
        return None;
    }
    for (&r, &continuum) in samples.iter().zip(&rebuilt) {
        let bound_sum: f64 = bound.iter().map(|&(n, w)| w * coulomb_value(z, l, n, r)).sum();
        if (value(r)? - bound_sum - continuum).abs() > 3e-2 * scale {
            return None;
        }
    }
    Some(rmax)
}

// ----------------------------------------------------------------------
// The evolution solver
// ----------------------------------------------------------------------

/// `v` as a small fraction, if it is one.
fn snap(
    cx: &mut Cx<'_>,
    v: f64,
) -> Option<NodeId> {
    for q in 1..=64_i64 {
        #[allow(clippy::cast_precision_loss)]
        let scaled = (v * q as f64).round();
        #[allow(clippy::cast_precision_loss)]
        if (v * q as f64 - scaled).abs() < 1e-9 {
            #[allow(clippy::cast_possible_truncation)]
            return Some(cx.graph.num(Number::fraction(scaled as i64, q)?));
        }
    }
    None
}

fn evaluate(
    cx: &Cx<'_>,
    node: NodeId,
    symbol: crate::graph::SymbolId,
    r: f64,
) -> Option<f64> {
    let mut env = Env::numeric(0.0);
    env.bind(symbol, r);
    cx.graph.eval(node, &env).filter(|v| v.is_finite())
}

/// `∫_0^{upper} f dr`: panels of width 1/4 up to 16 (the oscillations of the
/// high Rydberg states near the nucleus), then of width 2.
fn integrate(
    upper: f64,
    f: &dyn Fn(f64) -> Option<f64>,
) -> Option<f64> {
    let near = integrate_range((0.0, upper.min(16.0)), 0.25, f)?;
    if upper <= 16.0 {
        return Some(near);
    }
    Some(near + integrate_range((16.0, upper), 2.0, f)?)
}

/// The Coulomb bound-state expansion of a half-line evolution problem.
#[allow(clippy::too_many_lines)]
pub(super) fn solve(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
) -> Option<NodeId> {
    if p.nonlinear || p.dimension() != 2 || !p.homogeneous(cx.graph) {
        return None;
    }
    let time = evolution_time(cx, p)?;
    let space = 1 - time;
    let (x, t) = (p.vars[space], p.vars[time]);
    let symbol = cx.graph.symbol_of(x)?;
    let zero_index = vec![0; 2];
    let (a2, a1) = (p.coefficient(cx.graph, &p.unit(time, 2)), p.coefficient(cx.graph, &p.unit(time, 1)));
    let (order, a) = if !cx.is_zero(a2) && cx.is_zero(a1) {
        (2_u32, a2)
    } else if cx.is_zero(a2) && !cx.is_zero(a1) {
        (1_u32, a1)
    } else {
        return None;
    };
    if !p.only(cx.graph, &[p.unit(time, order), zero_index.clone(), p.unit(space, 2)]) {
        return None;
    }
    let b = p.coefficient(cx.graph, &p.unit(space, 2));
    let c = p.coefficient(cx.graph, &zero_index);
    // The signs: a u_t + b u_rr + c u = 0 in the problem's residual form, so u_t = D (u_rr + g u), D = -b/a, g = c/b.
    let b_value = exact_value(cx, b)?;
    if b_value == 0.0 {
        return None;
    }
    let g = |cx: &Cx<'_>, r: f64| -> Option<f64> { Some(evaluate(cx, c, symbol, r)? / b_value) };
    let f_of = |cx: &Cx<'_>, r: f64| -> Option<f64> { Some(r * r * g(cx, r)?) };
    let (f1, f2, f3) = (f_of(cx, 1.0)?, f_of(cx, 2.0)?, f_of(cx, 3.0)?);
    // F(r) = e0 r² + β r + γ0.
    let e0 = f64::midpoint(f3 - 2.0 * f2, f1);
    let beta = f2 - f1 - 3.0 * e0;
    let gamma0 = f1 - e0 - beta;
    for r in [0.7, 5.3] {
        let model = e0 * r * r + beta * r + gamma0;
        if (f_of(cx, r)? - model).abs() > 1e-9 * (1.0 + model.abs()) {
            return None;
        }
    }
    let l_real = f64::midpoint(-1.0, (1.0 - 4.0 * gamma0).sqrt());
    let l = whole((l_real + 1e-9).floor()).filter(|&l| (l_real - l as f64).abs() < 1e-9)?;
    let z = beta / 2.0;
    if z <= 0.0 {
        return None;
    }
    let (z_node, e0_node) = (snap(cx, z)?, snap(cx, e0)?);
    // Conditions: data, and the regularity condition at the origin.
    let displacement = conditions.find(cx.graph, time, &zero_index, None)?;
    let velocity = conditions.find(cx.graph, time, &p.unit(time, 1), None);
    let mut counted = 1 + usize::from(velocity.is_some());
    for cond in &conditions.0 {
        if cond.on == space {
            if cond.robin.is_some() || cond.derivative != zero_index || !is_zero_number(cx.graph, cond.point) || !cx.is_zero(cond.value) {
                return None;
            }
            counted += 1;
        }
    }
    if counted != conditions.0.len() || (velocity.is_some() && order == 1) {
        return None;
    }
    // Projection of the data on the bound states.
    let upper = 5.0 * (MAX_STATE * MAX_STATE) as f64 / z + 60.0 / z;
    let project = |cx: &Cx<'_>, data: NodeId| -> Option<(Vec<(usize, f64)>, bool)> {
        let mut out = Vec::new();
        for n in (l + 1)..=MAX_STATE {
            let coefficient = integrate(upper, &|r| Some(evaluate(cx, data, symbol, r)? * coulomb_value(z, l, n, r)))?;
            out.push((n, coefficient));
        }
        let scale = (1..=40).map(|k| evaluate(cx, data, symbol, 0.5 * f64::from(k)).map_or(0.0, f64::abs)).fold(1e-12, f64::max);
        let mut reproduced = true;
        for k in 1..=30 {
            let r = 0.4 * f64::from(k);
            let rebuilt: f64 = out.iter().map(|&(n, w)| w * coulomb_value(z, l, n, r)).sum();
            if (evaluate(cx, data, symbol, r)? - rebuilt).abs() > 1e-6 * scale {
                reproduced = false;
            }
        }
        Some((out, reproduced))
    };
    let free = cx.graph.free_symbols(cx.graph.find(displacement.value)).to_vec();
    if free.iter().any(|&s| s != symbol) {
        return None;
    }
    let (position, position_ok) = project(cx, displacement.value)?;
    let (rate, rate_ok) = match velocity {
        | Some(v) if !cx.is_zero(v.value) => project(cx, v.value)?,
        | _ => (Vec::new(), true),
    };
    // Data outside the bound-state span need the scattering states.
    let mut continuum: Vec<(NodeId, f64, bool)> = Vec::new();
    if !position_ok {
        let support = continuum_support(cx, displacement.value, symbol, (z, l), &position)?;
        continuum.push((displacement.value, support, true));
    }
    if !rate_ok {
        let v = velocity?;
        let support = continuum_support(cx, v.value, symbol, (z, l), &rate)?;
        continuum.push((v.value, support, false));
    }
    // The exact mode and the time factors.
    let diffusivity = {
        let q = div(cx, b, a);
        let q = neg(cx.graph, q);
        cx.simplify(q)
    };
    let mut terms = Vec::new();
    let mode = |cx: &mut Cx<'_>, n: usize| -> Option<NodeId> {
        let k = n - l - 1;
        let nf = n as i64;
        let norm_sq = {
            let two_z = {
                let two = cx.graph.int(2);
                mul(cx.graph, &[two, z_node])
            };
            let ratio = {
                let n_node = cx.graph.int(nf);
                div(cx, two_z, n_node)
            };
            let cube = powi(cx.graph, ratio, 3);
            let top = cx.graph.num(Number::from(factorial(k)));
            let bottom = cx.graph.num(Number::from(2.0 * n as f64 * factorial(n + l)));
            let q = div(cx, top, bottom);
            mul(cx.graph, &[cube, q])
        };
        let norm = sqrt(cx, norm_sq)?;
        let n_node = cx.graph.int(nf);
        let kappa = div(cx, z_node, n_node);
        let rho = {
            let two = cx.graph.int(2);
            mul(cx.graph, &[two, kappa, x])
        };
        let decay = {
            let e = mul(cx.graph, &[kappa, x]);
            let e = neg(cx.graph, e);
            call(cx, "exp", &[e])?
        };
        let mut poly = Vec::new();
        for j in 0..=k {
            let binomial = factorial(k + 2 * l + 1) / (factorial(k - j) * factorial(j + 2 * l + 1));
            let sign = if j % 2 == 0 { 1.0 } else { -1.0 };
            let coefficient = cx.graph.num(Number::from(sign * binomial / factorial(j)));
            let power = powi(cx.graph, rho, i64::try_from(j).ok()?);
            poly.push(mul(cx.graph, &[coefficient, power]));
        }
        let laguerre_poly = add(cx.graph, &poly);
        let radial = powi(cx.graph, x, i64::try_from(l + 1).ok()?);
        Some(mul(cx.graph, &[norm, radial, decay, laguerre_poly]))
    };
    let amplitude_peak = position.iter().fold(1e-30_f64, |a, &(_, w)| a.max(w.abs()));
    let speed_peak = rate.iter().fold(1e-30_f64, |a, &(_, w)| a.max(w.abs()));
    let ts: Vec<usize> = position.iter().chain(rate.iter()).map(|&(n, _)| n).collect();
    let mut seen = Vec::new();
    for n in ts {
        if seen.contains(&n) {
            continue;
        }
        seen.push(n);
        let amplitude = position.iter().find(|&&(m, _)| m == n).map_or(0.0, |&(_, w)| w);
        let speed = rate.iter().find(|&&(m, _)| m == n).map_or(0.0, |&(_, w)| w);
        if amplitude.abs() < 1e-7 * amplitude_peak && speed.abs() < 1e-7 * speed_peak {
            continue;
        }
        // λ_n = e0 + Z²/n².
        let lambda = {
            let n_node = cx.graph.int(n as i64);
            let ratio = div(cx, z_node, n_node);
            let sq = powi(cx.graph, ratio, 2);
            add(cx.graph, &[e0_node, sq])
        };
        let mode_node = mode(cx, n)?;
        let factor = if order == 1 {
            let exponent = mul(cx.graph, &[diffusivity, lambda, t]);
            call(cx, "exp", &[exponent])?
        } else {
            // T'' = D λ T: ω² = -D λ.
            let d_value = exact_value(cx, diffusivity)?;
            let lambda_value = exact_value(cx, lambda)?;
            let omega2 = -d_value * lambda_value;
            let w2_node = {
                let q = mul(cx.graph, &[diffusivity, lambda]);
                neg(cx.graph, q)
            };
            let w2_node = cx.simplify(w2_node);
            let mut parts = Vec::new();
            if omega2 > 0.0 {
                let omega = sqrt(cx, w2_node)?;
                let wt = mul(cx.graph, &[omega, t]);
                if amplitude.abs() >= 1e-10 {
                    let a_node = cx.graph.num(Number::from(amplitude));
                    let cs = call(cx, "cos", &[wt])?;
                    parts.push(mul(cx.graph, &[a_node, cs]));
                }
                if speed.abs() >= 1e-10 {
                    let s_node = cx.graph.num(Number::from(speed));
                    let sn = call(cx, "sin", &[wt])?;
                    let q = div(cx, sn, omega);
                    parts.push(mul(cx.graph, &[s_node, q]));
                }
            } else if omega2 < 0.0 {
                let w2n = neg(cx.graph, w2_node);
                let omega = sqrt(cx, w2n)?;
                let wt = mul(cx.graph, &[omega, t]);
                if amplitude.abs() >= 1e-10 {
                    let a_node = cx.graph.num(Number::from(amplitude));
                    let cs = call(cx, "cosh", &[wt])?;
                    parts.push(mul(cx.graph, &[a_node, cs]));
                }
                if speed.abs() >= 1e-10 {
                    let s_node = cx.graph.num(Number::from(speed));
                    let sn = call(cx, "sinh", &[wt])?;
                    let q = div(cx, sn, omega);
                    parts.push(mul(cx.graph, &[s_node, q]));
                }
            } else {
                if amplitude.abs() >= 1e-10 {
                    parts.push(cx.graph.num(Number::from(amplitude)));
                }
                if speed.abs() >= 1e-10 {
                    let s_node = cx.graph.num(Number::from(speed));
                    parts.push(mul(cx.graph, &[s_node, t]));
                }
            }
            let combined = add(cx.graph, &parts);
            terms.push(mul(cx.graph, &[mode_node, combined]));
            continue;
        };
        let a_node = cx.graph.num(Number::from(amplitude));
        terms.push(mul(cx.graph, &[a_node, factor, mode_node]));
    }
    for (data, support, is_displacement) in continuum {
        if order == 2 {
            let d_value = exact_value(cx, diffusivity)?;
            if d_value <= 0.0 || e0 > 0.0 {
                return None;
            }
        }
        let (k, _) = super::dummy(cx, p, "k");
        let (rho, _) = super::dummy(cx, p, "rho");
        let (zero, l_node) = (cx.graph.int(0), cx.graph.int(i64::try_from(l).ok()?));
        let upper_k = cx.graph.int(K_MAX as i64);
        let upper_r = cx.graph.num(Number::from(support));
        let wave_at = |cx: &mut Cx<'_>, point: NodeId| call(cx, "coulomb_wave", &[z_node, l_node, k, point]);
        let (at_rho, at_x) = (wave_at(cx, rho)?, wave_at(cx, x)?);
        let f_rho = cx.graph.substitute(data, x, rho);
        let inner_body = mul(cx.graph, &[f_rho, at_rho]);
        let coefficient = defint(cx, inner_body, rho, zero, upper_r)?;
        let k2 = powi(cx.graph, k, 2);
        let factor = if order == 1 {
            let spread = sub(cx.graph, e0_node, k2);
            let exponent = mul(cx.graph, &[diffusivity, spread, t]);
            call(cx, "exp", &[exponent])?
        } else {
            let spread = sub(cx.graph, k2, e0_node);
            let w2 = mul(cx.graph, &[diffusivity, spread]);
            let omega = sqrt(cx, w2)?;
            let wt = mul(cx.graph, &[omega, t]);
            if is_displacement {
                call(cx, "cos", &[wt])?
            } else {
                let sn = call(cx, "sin", &[wt])?;
                div(cx, sn, omega)
            }
        };
        let body = mul(cx.graph, &[coefficient, factor, at_x]);
        terms.push(defint(cx, body, k, zero, upper_k)?);
    }
    if terms.is_empty() {
        return None;
    }
    let u = add(cx.graph, &terms);
    Some(cx.simplify(u))
}
