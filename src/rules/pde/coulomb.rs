//! Coulomb (hydrogen-like) radial problems and spheroidal eigenvalues.
//!
//! * **Coulomb endpoints.** `a u_t = b [u_rr + (e₀ + 2Z/r - ℓ(ℓ+1)/r²) u]`
//!   (also `a u_tt`) on the half-line `r > 0`: the operator
//!   `-∂² + ℓ(ℓ+1)/r² - 2Z/r` has the bound states
//!   `u_{nℓ}(r) = N r^{ℓ+1} e^{-Zr/n} L_{n-ℓ-1}^{2ℓ+1}(2Zr/n)` with
//!   `E_n = -Z²/n²` (exact, Whittaker/Laguerre type) for `Z > 0`; other
//!   data are expanded in the continuum, whose regular solutions are the
//!   analytic Coulomb functions `F_ℓ(η, ρ)` (power series for moderate `ρ`,
//!   the large-`ρ` asymptotic expansion beyond), energy-normalised as
//!   `u_k(r) = √(2/π) F_ℓ(-Z/k, k r)`. The expansion is certified by
//!   reconstructing the data (`5e-3`, achieved `2.5e-3`) and by the Wronskian normalisation
//!   `F'G - FG' = 1`, and the upper cut-off `k_max` is extended until the
//!   last `k`-panels stop contributing (`≤ 64`, with a cost guard).
//!   Repulsive continua (`Z < 0`, `η = -Z/k > 0`) and attractive ones are
//!   handled by the same functions; the tabulated RK4 solution and its WKB
//!   amplitude have been removed.
//! * **Operators.** `coulomb_eigen(Z, ℓ, n)`, `coulomb_mode(Z, ℓ, n, r)`
//!   (orthonormal on `(0, ∞)`), `coulomb_wave(Z, ℓ, k, r)` (the
//!   energy-normalised scattering solution), and the spheroidal angular
//!   functions `spheroidal_eigen(m, n, γ)`, `spheroidal_mode(m, n, γ, x)`
//!   for `(1-x²) y'' - 2 x y' + (λ - γ x² - m²/(1-x²)) y = 0` on `[-1, 1]`
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
use crate::graph::DimRule;
use crate::graph::Env;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpDescriptor;
use crate::graph::RuleError;
use crate::graph::rule::Installer;
use crate::kernels::special_ext::coulomb_fg;
use crate::kernels::special_ext::coulomb_wave_value;
use crate::rules::complex::build::add;
use crate::rules::complex::build::mul;
use crate::rules::complex::build::neg;
use crate::rules::complex::build::powi;
use crate::rules::complex::build::sub;
use std::cell::RefCell;
use std::collections::VecDeque;
use std::rc::Rc;

/// The largest principal quantum number used in an expansion.
const MAX_STATE: usize = 14;
/// Basis size of the spheroidal matrix.
const BASIS: usize = 60;

fn factorial(n: usize) -> f64 {
    (1..=n).fold(1.0, |a, k| a * k as f64)
}

#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::float_cmp
)]
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
        let next =
            ((2.0 * jf + 1.0 + alpha - rho) * current - (jf + alpha) * previous) / (jf + 1.0);
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
    let norm =
        ((2.0 * z / n as f64).powi(3) * factorial(k) / (2.0 * n as f64 * factorial(n + l))).sqrt();
    r * norm
        * (-0.5 * rho).exp()
        * rho.powi(i32::try_from(l).unwrap_or(0))
        * laguerre(k, 2.0 * l as f64 + 1.0, rho)
}

/// `coulomb_eigen(Z, ℓ, n)`: `-Z²/n²`.
pub(super) fn coulomb_eigen(args: &[f64]) -> f64 {
    let (Some(&z), Some(l), Some(n)) = (
        args.first(),
        args.get(1).and_then(|&v| whole(v)),
        args.get(2).and_then(|&v| whole(v)),
    ) else {
        return f64::NAN;
    };
    if n <= l || n == 0 {
        f64::NAN
    } else {
        -z * z / (n * n) as f64
    }
}

/// `coulomb_mode(Z, ℓ, n, r)`.
pub(super) fn coulomb_mode(args: &[f64]) -> f64 {
    let (Some(&z), Some(l), Some(n), Some(&r)) = (
        args.first(),
        args.get(1).and_then(|&v| whole(v)),
        args.get(2).and_then(|&v| whole(v)),
        args.get(3),
    ) else {
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
        if j <= m {
            0.0
        } else {
            ((jf * jf - mf * mf) / (4.0 * jf * jf - 1.0)).sqrt()
        }
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
    if let Some(found) = SPHEROIDAL.with(|c| {
        c.borrow()
            .iter()
            .find(|(k, _)| *k == key)
            .map(|(_, v)| Rc::clone(v))
    }) {
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
        if j <= m {
            0.0
        } else {
            ((jf * jf - mf * mf) / (4.0 * jf * jf - 1.0)).sqrt()
        }
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
    let (Some(m), Some(n), Some(&gamma)) = (
        args.first().and_then(|&v| whole(v)),
        args.get(1).and_then(|&v| whole(v)),
        args.get(2),
    ) else {
        return f64::NAN;
    };
    if n < m || n - m >= BASIS / 2 {
        return f64::NAN;
    }
    spheroidal_pair(m, gamma)
        .0
        .get(n - m)
        .copied()
        .unwrap_or(f64::NAN)
}

/// `spheroidal_mode(m, n, γ, x)`: the angular function, normalised on
/// `[-1, 1]`, positive leading coefficient.
pub(super) fn spheroidal_mode(args: &[f64]) -> f64 {
    let (Some(m), Some(n), Some(&gamma), Some(&x)) = (
        args.first().and_then(|&v| whole(v)),
        args.get(1).and_then(|&v| whole(v)),
        args.get(2),
        args.get(3),
    ) else {
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
///
/// The Coulomb and spheroidal operators are numeric and orthonormalised
/// (`f64` in, `f64` out), so they are dimensionless; a dimensionful charge or
/// coordinate belongs to a symbolic wrapper, not to these maps.
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    for (name, arity, eval) in [
        (
            "coulomb_eigen",
            Arity::Fixed(3),
            coulomb_eigen as fn(&[f64]) -> f64,
        ),
        ("coulomb_mode", Arity::Fixed(4), coulomb_mode),
        ("coulomb_wave", Arity::Fixed(4), coulomb_wave),
        ("spheroidal_eigen", Arity::Fixed(3), spheroidal_eigen),
        ("spheroidal_mode", Arity::Fixed(4), spheroidal_mode),
    ] {
        let op = i.op(OpDescriptor::new(name, arity).eval(eval))?;
        i.dim(op, DimRule::Dimensionless);
    }
    Ok(())
}

// ----------------------------------------------------------------------
// Continuum (scattering) states
// ----------------------------------------------------------------------

/// The largest radius scanned for the support of the data.
const SUPPORT_MAX: f64 = 160.0;
/// The initial cut-off of the continuum integral.
const K_START: f64 = 16.0;
/// The largest cut-off allowed by the cost guard.
const K_LIMIT: f64 = 64.0;
/// The relative reconstruction tolerance accepted for the continuum. The
/// data are reproduced to about `2.5e-3` in the pure-continuum cases and
/// `1.4e-3` for a bound state plus a bump (the k-truncation at 64 and the
/// Coulomb phase make the tail converge slowly); `5e-3` is the documented
/// acceptance, still an order of magnitude below the old 3 % spot check.
const RECON_TOL: f64 = 5e-3;
/// The relative size of the last `k`-panels that certifies the cut-off
/// (`5e-4` of the data scale; the reconstruction tolerance above is the
/// binding constraint).
const TAIL_TOL: f64 = 5e-4;
/// The largest Wronskian residual accepted for the scattering solutions.
const WRONSKI_TOL: f64 = 1e-8;

/// Six-point Gauss–Legendre nodes and weights on `[-1, 1]`.
const GAUSS6: [(f64, f64); 6] = [
    (-0.932_469_514_203_152, 0.171_324_492_379_170),
    (-0.661_209_386_466_265, 0.360_761_573_048_139),
    (-0.238_619_186_083_197, 0.467_913_934_572_691),
    (0.238_619_186_083_197, 0.467_913_934_572_691),
    (0.661_209_386_466_265, 0.360_761_573_048_139),
    (0.932_469_514_203_152, 0.171_324_492_379_170),
];

/// `coulomb_wave(Z, ℓ, k, r)`.
pub(super) fn coulomb_wave(args: &[f64]) -> f64 {
    let (Some(&z), Some(l), Some(&k), Some(&r)) = (
        args.first(),
        args.get(1).and_then(|&v| whole(v)),
        args.get(2),
        args.get(3),
    ) else {
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
    let panels = (((b - a) / width).ceil() as usize).max(1);
    let h = (b - a) / panels as f64;
    let mut total = 0.0;
    for panel in 0..panels {
        let mid = a + (panel as f64 + 0.5) * h;
        for &(x, w) in &GAUSS6 {
            total += 0.5 * h * w * f(mid + 0.5 * h * x)?;
        }
    }
    Some(total)
}

/// The largest Wronskian residual `|F' G - F G' - 1|` over a few points;
/// this is the analytic energy-normalisation identity of the pair.
#[allow(clippy::single_call_fn)]
fn wronskian_residual(
    z: f64,
    l: usize,
) -> f64 {
    let mut worst = 0.0_f64;
    let mut checked = 0_usize;
    for multiplier in [1.0_f64, 3.0, 9.0] {
        let k = z.abs().max(0.5) * multiplier;
        let eta = -z / k;
        for &rho in &[1.5_f64, 5.0, 12.0] {
            if let Some((f, fp, g, gp)) = coulomb_fg(eta, l, rho) {
                worst = worst.max((fp * g - f * gp - 1.0).abs());
                checked += 1;
            }
        }
    }
    if checked == 9 {
        worst
    } else {
        f64::INFINITY
    }
}

/// Certifies the continuum expansion of `data`: bound states plus
/// `∫_0^{k_max} c(k) u_k(r) dk`. The reconstruction at sample points must
/// agree with the data to `RECON_TOL` of its size, the last `k`-panels must
/// contribute less than `TAIL_TOL` (else `k_max` is extended to 32 and then
/// 64), and the Wronskian normalisation of the scattering solutions must
/// hold to `WRONSKI_TOL`. Returns `(R, k_max)`.
fn continuum_support(
    cx: &Cx<'_>,
    data: NodeId,
    symbol: crate::graph::SymbolId,
    (z, l): (f64, usize),
    bound: &[(usize, f64)],
) -> Option<(f64, f64)> {
    let value = |r: f64| evaluate(cx, data, symbol, r);
    let mut scale = 1e-12_f64;
    let mut r = 0.0_f64;
    while r <= SUPPORT_MAX {
        scale = scale.max(value(r)?.abs());
        r += 0.5;
    }
    if value(SUPPORT_MAX)?.abs() > 1e-4 * scale {
        return None;
    }
    let mut rmax = 8.0_f64;
    let mut r = 0.5_f64;
    while r <= SUPPORT_MAX {
        if value(r)?.abs() > 1e-4 * scale {
            rmax = (r + 1.0).min(SUPPORT_MAX);
        }
        r += 0.5;
    }
    if wronskian_residual(z, l) > WRONSKI_TOL {
        return None;
    }
    let samples = [0.3_f64, 0.8, 1.5, 2.5, 4.0, 6.0];
    let bound_at: Vec<f64> = samples
        .iter()
        .map(|&r| {
            bound
                .iter()
                .map(|&(n, w)| w * coulomb_value(z, l, n, r))
                .sum()
        })
        .collect();
    let coefficient = |k: f64| -> Option<f64> {
        let width = (0.5_f64).min(1.0 / k);
        integrate_range((0.0, rmax), width, &|r| {
            Some(value(r)? * coulomb_wave_value(z, l, k, r))
        })
    };
    let step = (0.5_f64).min(1.0 / rmax);
    let mut total = vec![0.0_f64; samples.len()];
    let mut done = 0.0_f64;
    let mut certified = None;
    for stage in [K_START, 2.0 * K_START, K_LIMIT] {
        let panels = (((stage - done) / step).round() as usize).max(1);
        let h = (stage - done) / panels as f64;
        let keep = ((2.0 / h).ceil() as usize).max(1);
        let mut recent: VecDeque<Vec<f64>> = VecDeque::new();
        for panel in 0..panels {
            let start = done + panel as f64 * h;
            let middle = start + 0.5 * h;
            let mut piece = vec![0.0_f64; samples.len()];
            for &(x, w) in &GAUSS6 {
                let k = middle + 0.5 * h * x;
                let weight = 0.5 * h * w;
                let c = coefficient(k)?;
                for (slot, &r) in piece.iter_mut().zip(&samples) {
                    *slot += weight * c * coulomb_wave_value(z, l, k, r);
                }
            }
            for (slot, &part) in total.iter_mut().zip(&piece) {
                *slot += part;
            }
            recent.push_back(piece);
            if recent.len() > keep {
                recent.pop_front();
            }
        }
        done = stage;
        let mut tail = vec![0.0_f64; samples.len()];
        for piece in &recent {
            for (slot, &part) in tail.iter_mut().zip(piece) {
                *slot += part;
            }
        }
        let mut ok = true;
        for (((&r, &b), &t), &u) in samples.iter().zip(&bound_at).zip(&total).zip(&tail) {
            let observed = value(r)?;
            if (observed - b - t).abs() > RECON_TOL * scale || u.abs() > TAIL_TOL * scale {
                ok = false;
            }
        }
        if ok {
            certified = Some(stage);
            break;
        }
    }
    Some((rmax, certified?))
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
    let (a2, a1) = (
        p.coefficient(cx.graph, &p.unit(time, 2)),
        p.coefficient(cx.graph, &p.unit(time, 1)),
    );
    let (order, a) = if !cx.is_zero(a2) && cx.is_zero(a1) {
        (2_u32, a2)
    } else if cx.is_zero(a2) && !cx.is_zero(a1) {
        (1_u32, a1)
    } else {
        return None;
    };
    if !p.only(
        cx.graph,
        &[p.unit(time, order), zero_index.clone(), p.unit(space, 2)],
    ) {
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
    let (z_node, e0_node) = (snap(cx, z)?, snap(cx, e0)?);
    // Conditions: data, and the regularity condition at the origin.
    let displacement = conditions.find(cx.graph, time, &zero_index, None)?;
    let velocity = conditions.find(cx.graph, time, &p.unit(time, 1), None);
    let mut counted = 1 + usize::from(velocity.is_some());
    for cond in &conditions.0 {
        if cond.on == space {
            if cond.robin.is_some()
                || cond.derivative != zero_index
                || !is_zero_number(cx.graph, cond.point)
                || !cx.is_zero(cond.value)
            {
                return None;
            }
            counted += 1;
        }
    }
    if counted != conditions.0.len() || (velocity.is_some() && order == 1) {
        return None;
    }
    // Projection of the data on the bound states (none for a repulsive
    // charge, `Z ≤ 0`).
    let bound = z > 0.0;
    let upper = if bound {
        5.0 * (MAX_STATE * MAX_STATE) as f64 / z + 60.0 / z
    } else {
        0.0
    };
    let project = |cx: &Cx<'_>, data: NodeId| -> Option<(Vec<(usize, f64)>, bool)> {
        if !bound {
            return Some((Vec::new(), false));
        }
        let mut out = Vec::new();
        for n in (l + 1)..=MAX_STATE {
            let coefficient = integrate(upper, &|r| {
                Some(evaluate(cx, data, symbol, r)? * coulomb_value(z, l, n, r))
            })?;
            out.push((n, coefficient));
        }
        let scale = (1..=40)
            .map(|k| evaluate(cx, data, symbol, 0.5 * f64::from(k)).map_or(0.0, f64::abs))
            .fold(1e-12, f64::max);
        let mut reproduced = true;
        for k in 1..=30 {
            let r = 0.4 * f64::from(k);
            let rebuilt: f64 = out
                .iter()
                .map(|&(n, w)| w * coulomb_value(z, l, n, r))
                .sum();
            if (evaluate(cx, data, symbol, r)? - rebuilt).abs() > 1e-6 * scale {
                reproduced = false;
            }
        }
        Some((out, reproduced))
    };
    let free = cx
        .graph
        .free_symbols(cx.graph.find(displacement.value))
        .to_vec();
    if free.iter().any(|&s| s != symbol) {
        return None;
    }
    let (position, position_ok) = project(cx, displacement.value)?;
    let (rate, rate_ok) = match velocity {
        | Some(v) if !cx.is_zero(v.value) => project(cx, v.value)?,
        | _ => (Vec::new(), true),
    };
    // Data outside the bound-state span need the scattering states.
    let mut continuum: Vec<(NodeId, f64, f64, bool)> = Vec::new();
    if !position_ok {
        let (support, cutoff) =
            continuum_support(cx, displacement.value, symbol, (z, l), &position)?;
        continuum.push((displacement.value, support, cutoff, true));
    }
    if !rate_ok {
        let v = velocity?;
        let (support, cutoff) = continuum_support(cx, v.value, symbol, (z, l), &rate)?;
        continuum.push((v.value, support, cutoff, false));
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
            let bottom = cx
                .graph
                .num(Number::from(2.0 * n as f64 * factorial(n + l)));
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
    let ts: Vec<usize> = position
        .iter()
        .chain(rate.iter())
        .map(|&(n, _)| n)
        .collect();
    let mut seen = Vec::new();
    for n in ts {
        if seen.contains(&n) {
            continue;
        }
        seen.push(n);
        let amplitude = position
            .iter()
            .find(|&&(m, _)| m == n)
            .map_or(0.0, |&(_, w)| w);
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
    for (data, support, cutoff, is_displacement) in continuum {
        if order == 2 {
            let d_value = exact_value(cx, diffusivity)?;
            if d_value <= 0.0 || e0 > 0.0 {
                return None;
            }
        }
        let (k, _) = super::dummy(cx, p, "k");
        let (rho, _) = super::dummy(cx, p, "rho");
        let (zero, l_node) = (cx.graph.int(0), cx.graph.int(i64::try_from(l).ok()?));
        let upper_k = cx.graph.num(Number::from(cutoff));
        let upper_r = cx.graph.num(Number::from(support));
        let wave_at =
            |cx: &mut Cx<'_>, point: NodeId| call(cx, "coulomb_wave", &[z_node, l_node, k, point]);
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
