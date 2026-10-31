//! General regular Sturm–Liouville problems with variable coefficients.
//!
//! `a(x) u_xx + b(x) u_x + c(x) u` on `[x₀, x₁]` (with `a` of one sign, `a`,
//! `b`, `c` smooth and free of other symbols) next to constant time
//! derivatives and Robin, Neumann or Dirichlet conditions with zero values
//! is brought to the self-adjoint form
//!
//! `-(p y')' + q y = μ w y`, `p = exp ∫ b̂/â`, `w = p/â`, `q = -p ĉ/â`
//!
//! (`â = κ a`, `κ = ±1`, so that `â > 0`; `S y = κ (â y'' + b̂ y' + ĉ y)`).
//! Closed forms (Liouville's transformation to normal form) exist only for
//! special coefficients, the Euler–Cauchy and Legendre families being
//! handled in `sturm`; for the general case the coefficients `p`, `q`, `w`
//! are sampled at Chebyshev points, and the eigenvalues and eigenfunctions
//! are numeric operators evaluated by the Prüfer transformation (shooting
//! on the angle `θ`, `y = ρ sin θ`, `p y' = ρ cos θ`) and Runge–Kutta
//! integration of the eigenfunction:
//!
//! * `sl_eigen(m, x₀, x₁, α₀, β₀, α₁, β₁, N, p…, q…, w…)`: the `m`-th
//!   eigenvalue `μ_m` (`m ≥ 1`, ascending) of the problem with
//!   `α₀ y + β₀ y' = 0` at `x₀` and `α₁ y + β₁ y' = 0` at `x₁`, the
//!   coefficients given by their `N + 1` values at the Chebyshev points
//!   `x_j = (x₀ + x₁)/2 + (x₁ - x₀)/2 cos(π j/N)`;
//! * `sl_mode(m, x, …)`: the eigenfunction, normalised by `∫ w y² = 1`
//!   with the sign fixed by `(y, p y')(x₀) ∝ (sin θ₀, cos θ₀)`;
//! * `sl_weight(x, x₀, x₁, …)`: the weight `w(x)` by barycentric
//!   interpolation.
//!
//! Separation of variables then yields `Σ c_m T_m(t) sl_mode(m, x, …)`,
//! the coefficients `c_m = ∫ w f y_m dx` being integrals of numeric
//! operators.

use super::Conditions;
use super::Problem;
use super::spectral::Axis;
use super::spectral::Engine;
use super::spectral::Item;
use super::spectral::Role;
use super::spectral::Shape;
use super::spectral::Symbol;
use super::spectral::Time;
use super::spectral::evolution_time;
use super::spectral::finish;
use super::spectral::interval;
use super::spectral::time_items;
use super::util::is_zero_number;
use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::DimRule;
use crate::graph::Env;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpDescriptor;
use crate::graph::RuleError;
use crate::graph::rule::Installer;
use std::cell::RefCell;
use std::f64::consts::PI;
use std::rc::Rc;

/// Number of Chebyshev intervals used to sample the coefficients.
const SAMPLES: usize = 24;

// ----------------------------------------------------------------------
// Numeric part
// ----------------------------------------------------------------------

/// Barycentric interpolation on Chebyshev points of the second kind.
struct Chebyshev {
    mid: f64,
    half: f64,
    nodes: Vec<f64>,
    values: Vec<f64>,
}

impl Chebyshev {
    fn new(
        a: f64,
        b: f64,
        values: Vec<f64>,
    ) -> Self {
        let n = values.len().saturating_sub(1).max(1);
        let nodes = (0..values.len())
            .map(|j| (PI * j as f64 / n as f64).cos())
            .collect();
        Self {
            mid: f64::midpoint(a, b),
            half: 0.5 * (b - a),
            nodes,
            values,
        }
    }

    fn eval(
        &self,
        x: f64,
    ) -> f64 {
        let t = (x - self.mid) / self.half;
        let last = self.values.len().saturating_sub(1);
        let (mut num, mut den) = (0.0, 0.0);
        for (j, (&tj, &f)) in self.nodes.iter().zip(&self.values).enumerate() {
            let d = t - tj;
            if d.abs() < 1e-15 {
                return f;
            }
            let sign = if j % 2 == 0 { 1.0 } else { -1.0 };
            let weight = if j == 0 || j == last {
                0.5 * sign
            } else {
                sign
            };
            let c = weight / d;
            num += c * f;
            den += c;
        }
        num / den
    }
}

/// A regular Sturm–Liouville problem `-(p y')' + q y = μ w y`.
struct Sl {
    a: f64,
    b: f64,
    ends: [f64; 4],
    p: Chebyshev,
    q: Chebyshev,
    w: Chebyshev,
}

/// The problem from `a, b, α₀, β₀, α₁, β₁, N, p…, q…, w…`.
#[allow(clippy::float_cmp)]
fn parse(args: &[f64]) -> Option<Sl> {
    let (&a, &b) = (args.first()?, args.get(1)?);
    let ends = [*args.get(2)?, *args.get(3)?, *args.get(4)?, *args.get(5)?];
    let n = *args.get(6)?;
    if !(2.0..=512.0).contains(&n) || n.fract() != 0.0 || b <= a || a.is_nan() || b.is_nan() {
        return None;
    }
    let count = n as usize + 1;
    let rest = args.get(7..)?;
    if rest.len() != 3 * count || rest.iter().any(|v| !v.is_finite()) {
        return None;
    }
    let (p, rest) = rest.split_at(count);
    let (q, w) = rest.split_at(count);
    let positive = p.iter().all(|&v| v > 0.0) && w.iter().all(|&v| v > 0.0);
    positive.then(|| {
        Sl {
            a,
            b,
            ends,
            p: Chebyshev::new(a, b, p.to_vec()),
            q: Chebyshev::new(a, b, q.to_vec()),
            w: Chebyshev::new(a, b, w.to_vec()),
        }
    })
}

impl Sl {
    /// The Prüfer angle at `x₀`, in `[0, π)`.
    fn start_angle(&self) -> f64 {
        let [alpha, beta, _, _] = self.ends;
        let t = (-beta / self.p.eval(self.a)).atan2(alpha).rem_euclid(PI);
        if PI - t < 1e-14 { 0.0 } else { t }
    }

    /// The target angle at `x₁` of the first eigenvalue, in `(0, π]`.
    fn end_angle(&self) -> f64 {
        let [_, _, alpha, beta] = self.ends;
        let t = (-beta / self.p.eval(self.b)).atan2(alpha).rem_euclid(PI);
        if t < 1e-14 { PI } else { t }
    }

    /// The number of integration steps for the `m`-th eigenvalue: a Weyl
    /// estimate `μ_m ≈ (m π / ∫ √(w/p))²` (bounded below by the potential)
    /// fixes the oscillation to be resolved.
    fn plan(
        &self,
        m: u32,
    ) -> usize {
        let samples = 400_usize;
        let h = (self.b - self.a) / samples as f64;
        let mut integral = 0.0;
        for k in 0..=samples {
            let x = self.a + h * k as f64;
            let weight = if k == 0 || k == samples {
                1.0
            } else if k % 2 == 1 {
                4.0
            } else {
                2.0
            };
            integral += weight * (self.w.eval(x) / self.p.eval(x)).sqrt();
        }
        integral *= h / 3.0;
        let (pmin, wmin) = (
            self.p.values.iter().copied().fold(f64::INFINITY, f64::min),
            self.w.values.iter().copied().fold(f64::INFINITY, f64::min),
        );
        let (wmax, qmax) = (
            self.w.values.iter().copied().fold(0.0, f64::max),
            self.q
                .values
                .iter()
                .copied()
                .fold(0.0, |a: f64, v| a.max(v.abs())),
        );
        let estimate = (f64::from(m) * PI / integral).powi(2).max(qmax / wmin);
        let frequency = (4.0 * estimate * wmax / pmin).sqrt() * (self.b - self.a);
        let steps = 3000.0 + 80.0 * frequency.powf(1.25);
        2 * (steps.min(600_000.0) as usize / 2)
    }

    /// The coefficients on the grid of `steps` intervals and their midpoints.
    fn grid(
        &self,
        steps: usize,
    ) -> Grid {
        let h = (self.b - self.a) / steps as f64;
        let points: Vec<f64> = (0..=2 * steps)
            .map(|k| self.a + 0.5 * h * k as f64)
            .collect();
        Grid {
            h,
            p: points.iter().map(|&x| self.p.eval(x)).collect(),
            q: points.iter().map(|&x| self.q.eval(x)).collect(),
            w: points.iter().map(|&x| self.w.eval(x)).collect(),
        }
    }

    /// The `m`-th eigenvalue (`m ≥ 1`) and the grid it was computed on.
    fn eigenvalue(
        &self,
        m: u32,
    ) -> Option<(f64, Grid)> {
        if m == 0 {
            return None;
        }
        let grid = self.grid(self.plan(m));
        let theta0 = self.start_angle();
        let target = self.end_angle() + f64::from(m - 1) * PI;
        let f = |mu: f64| grid.end_value(mu, theta0) - target;
        let mut hi = 1.0;
        for _ in 0..60 {
            if f(hi) >= 0.0 {
                break;
            }
            hi *= 2.0;
        }
        let mut lo = -1.0;
        for _ in 0..60 {
            if f(lo) <= 0.0 {
                break;
            }
            lo *= 2.0;
        }
        if f(hi) < 0.0 || f(lo) > 0.0 {
            return None;
        }
        for _ in 0..200 {
            let mid = f64::midpoint(lo, hi);
            if f(mid) < 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
            if hi - lo <= 1e-15 * mid.abs().max(1.0) {
                break;
            }
        }
        Some((f64::midpoint(lo, hi), grid))
    }
}

/// The coefficients at the grid points `x_k = a + k h/2`.
struct Grid {
    h: f64,
    p: Vec<f64>,
    q: Vec<f64>,
    w: Vec<f64>,
}

impl Grid {
    const fn steps(&self) -> usize {
        self.p.len().saturating_sub(1) / 2
    }

    /// `(p, q, w)` at the point number `k`.
    fn at(
        &self,
        k: usize,
    ) -> (f64, f64, f64) {
        (
            self.p.get(k).copied().unwrap_or(f64::NAN),
            self.q.get(k).copied().unwrap_or(f64::NAN),
            self.w.get(k).copied().unwrap_or(f64::NAN),
        )
    }

    /// `θ(x₁; μ)` for the Prüfer equation `θ' = cos²θ/p + (μ w - q) sin²θ`
    /// started at `theta0`, by the classical Runge–Kutta method.
    fn end_value(
        &self,
        mu: f64,
        theta0: f64,
    ) -> f64 {
        let rhs = |k: usize, theta: f64| {
            let (p, q, w) = self.at(k);
            let (s, c) = theta.sin_cos();
            c * c / p + (mu * w - q) * s * s
        };
        let h = self.h;
        let mut theta = theta0;
        for k in 0..self.steps() {
            let i = 2 * k;
            let k1 = rhs(i, theta);
            let k2 = rhs(i + 1, theta + 0.5 * h * k1);
            let k3 = rhs(i + 1, theta + 0.5 * h * k2);
            let k4 = rhs(i + 2, theta + h * k3);
            theta += h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        }
        theta
    }
}

/// An eigenfunction on a uniform grid (values and derivatives), normalised.
struct Mode {
    a: f64,
    h: f64,
    y: Vec<f64>,
    dy: Vec<f64>,
}

impl Mode {
    fn build(
        sl: &Sl,
        grid: &Grid,
        mu: f64,
    ) -> Self {
        let steps = grid.steps();
        let h = grid.h;
        let theta = sl.start_angle();
        let (mut y, mut v) = (theta.sin(), theta.cos());
        let rhs = |k: usize, y: f64, v: f64| {
            let (p, q, w) = grid.at(k);
            (v / p, (q - mu * w) * y)
        };
        let (mut ys, mut dys, mut ws) = (
            Vec::with_capacity(steps + 1),
            Vec::with_capacity(steps + 1),
            Vec::with_capacity(steps + 1),
        );
        for k in 0..=steps {
            let (p, _, w) = grid.at(2 * k);
            ys.push(y);
            dys.push(v / p);
            ws.push(w);
            if k == steps {
                break;
            }
            let i = 2 * k;
            let (a1, b1) = rhs(i, y, v);
            let (a2, b2) = rhs(i + 1, y + 0.5 * h * a1, v + 0.5 * h * b1);
            let (a3, b3) = rhs(i + 1, y + 0.5 * h * a2, v + 0.5 * h * b2);
            let (a4, b4) = rhs(i + 2, y + h * a3, v + h * b3);
            y += h / 6.0 * (a1 + 2.0 * a2 + 2.0 * a3 + a4);
            v += h / 6.0 * (b1 + 2.0 * b2 + 2.0 * b3 + b4);
        }
        // ∫ w y² by Simpson's rule.
        let square = |k: usize| {
            ws.get(k).copied().unwrap_or(0.0) * ys.get(k).copied().unwrap_or(0.0).powi(2)
        };
        let mut total = square(0) + square(steps);
        for k in 1..steps {
            total += square(k) * if k % 2 == 1 { 4.0 } else { 2.0 };
        }
        let norm = (total * h / 3.0).sqrt();
        let scale = if norm > 0.0 {
            1.0 / norm
        } else {
            f64::NAN
        };
        Self {
            a: sl.a,
            h,
            y: ys.iter().map(|v| v * scale).collect(),
            dy: dys.iter().map(|v| v * scale).collect(),
        }
    }

    /// Cubic Hermite interpolation.
    fn eval(
        &self,
        x: f64,
    ) -> f64 {
        let last = self.y.len().saturating_sub(2);
        let position = ((x - self.a) / self.h).max(0.0);
        let k = (position.floor() as usize).min(last);
        let t = position - k as f64;
        let (Some(&y0), Some(&y1), Some(&d0), Some(&d1)) = (
            self.y.get(k),
            self.y.get(k + 1),
            self.dy.get(k),
            self.dy.get(k + 1),
        ) else {
            return f64::NAN;
        };
        let (t2, t3) = (t * t, t * t * t);
        (2.0 * t3 - 3.0 * t2 + 1.0) * y0
            + (t3 - 2.0 * t2 + t) * self.h * d0
            + (-2.0 * t3 + 3.0 * t2) * y1
            + (t3 - t2) * self.h * d1
    }
}

type Cached = (Vec<u64>, Rc<(f64, Mode)>);

thread_local! {
    static CACHE: RefCell<Vec<Cached>> = const { RefCell::new(Vec::new()) };
}

/// The eigenvalue and eigenfunction of mode `m` (cached per problem).
fn solved(
    m: u32,
    args: &[f64],
) -> Option<Rc<(f64, Mode)>> {
    let mut key: Vec<u64> = args.iter().map(|v| v.to_bits()).collect();
    key.push(u64::from(m));
    let cached = CACHE.with(|c| {
        c.borrow()
            .iter()
            .find(|(k, _)| *k == key)
            .map(|(_, v)| Rc::clone(v))
    });
    if cached.is_some() {
        return cached;
    }
    let sl = parse(args)?;
    let (mu, grid) = sl.eigenvalue(m)?;
    let value = Rc::new((mu, Mode::build(&sl, &grid, mu)));
    CACHE.with(|c| {
        let mut cache = c.borrow_mut();
        if cache.len() >= 64 {
            cache.remove(0);
        }
        cache.push((key, Rc::clone(&value)));
    });
    Some(value)
}

#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::float_cmp
)]
fn index(m: f64) -> Option<u32> {
    if (1.0..=1e6).contains(&m) && m.fract() == 0.0 {
        Some(m as u32)
    } else {
        None
    }
}

/// `sl_eigen(m, a, b, α₀, β₀, α₁, β₁, N, p…, q…, w…)`.
pub(super) fn sl_eigen(args: &[f64]) -> f64 {
    let (Some(m), Some(rest)) = (args.first().and_then(|&m| index(m)), args.get(1..)) else {
        return f64::NAN;
    };
    solved(m, rest).map_or(f64::NAN, |s| s.0)
}

/// `sl_mode(m, x, a, b, α₀, β₀, α₁, β₁, N, p…, q…, w…)`.
pub(super) fn sl_mode(args: &[f64]) -> f64 {
    let (Some(m), Some(&x), Some(rest)) = (
        args.first().and_then(|&m| index(m)),
        args.get(1),
        args.get(2..),
    ) else {
        return f64::NAN;
    };
    solved(m, rest).map_or(f64::NAN, |s| s.1.eval(x))
}

/// `sl_weight(x, a, b, α₀, β₀, α₁, β₁, N, p…, q…, w…)`.
pub(super) fn sl_weight(args: &[f64]) -> f64 {
    let (Some(&x), Some(rest)) = (args.first(), args.get(1..)) else {
        return f64::NAN;
    };
    parse(rest).map_or(f64::NAN, |sl| sl.w.eval(x))
}

/// Registers the operators.
///
/// The Sturm-Liouville maps are numeric (`f64` coefficients at Chebyshev
/// points in, `f64` out), so they are dimensionless; the physical eigenvalue
/// `μ` of `-(p y')' + q y = μ w y` belongs to the symbolic expansion around
/// them.
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    for (name, eval) in [
        ("sl_eigen", sl_eigen as fn(&[f64]) -> f64),
        ("sl_mode", sl_mode),
        ("sl_weight", sl_weight),
    ] {
        let op = i.op(OpDescriptor::new(name, Arity::Variadic).eval(eval))?;
        i.dim(op, DimRule::Dimensionless);
    }
    Ok(())
}

// ----------------------------------------------------------------------
// The symbolic solver
// ----------------------------------------------------------------------

/// The value of `node` with the space variables set to `point`.
fn at_point(
    cx: &Cx<'_>,
    node: NodeId,
    symbols: &[crate::graph::SymbolId],
    point: &[f64],
) -> Option<f64> {
    let mut env = Env::numeric(0.0);
    for (&symbol, &value) in symbols.iter().zip(point) {
        env.bind(symbol, value);
    }
    cx.graph.eval(node, &env).filter(|v| v.is_finite())
}

/// One axis of a (separable) problem: its interval, numbers and the
/// coefficient nodes of `a u_xx + b u_x`.
struct AxisData {
    var: NodeId,
    start: NodeId,
    stop: NodeId,
    ends: [NodeId; 4],
    bounds: (f64, f64),
    a: NodeId,
    b: NodeId,
}

/// Separation of variables for `Σ_j [a_j(x_j) u_jj + b_j(x_j) u_j] + c(x) u`
/// with general smooth coefficients on boxes (one to three space variables,
/// `c` a sum of functions of single variables plus a constant), by the
/// numeric Sturm–Liouville operators of each axis.
#[allow(clippy::too_many_lines)]
pub(super) fn solve(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
) -> Option<NodeId> {
    if p.nonlinear {
        return None;
    }
    let n = p.dimension();
    let time_index = evolution_time(cx, p)?;
    let space: Vec<usize> = (0..n).filter(|&j| j != time_index).collect();
    if space.is_empty() || space.len() > 3 {
        return None;
    }
    let (a1, a2) = (
        p.coefficient(cx.graph, &p.unit(time_index, 1)),
        p.coefficient(cx.graph, &p.unit(time_index, 2)),
    );
    if !p.constant(cx.graph, a1) || !p.constant(cx.graph, a2) {
        return None;
    }
    for (index, c) in &p.linear {
        let known = *index == vec![0; n]
            || *index == p.unit(time_index, 1)
            || *index == p.unit(time_index, 2)
            || space
                .iter()
                .any(|&j| *index == p.unit(j, 1) || *index == p.unit(j, 2));
        if !known && !is_zero_number(cx.graph, *c) {
            return None;
        }
    }
    let value = |cx: &Cx<'_>, node: NodeId| super::modes::exact_value(cx, node);
    let cc = p.coefficient(cx.graph, &vec![0; n]);
    let mut data = Vec::new();
    for &j in &space {
        let ca = p.coefficient(cx.graph, &p.unit(j, 2));
        let cb = p.coefficient(cx.graph, &p.unit(j, 1));
        if cx.is_zero(ca) {
            return None;
        }
        // The coefficients depend on x_j alone.
        for &k in &space {
            if k == j {
                continue;
            }
            let other = cx.graph.symbol_of(p.vars[k])?;
            if cx.graph.depends_on(cx.graph.find(ca), other)
                || cx.graph.depends_on(cx.graph.find(cb), other)
            {
                return None;
            }
        }
        let iv = interval(cx, p, conditions, j)?;
        if iv.periodic || !cx.is_zero(iv.left.value) || !cx.is_zero(iv.right.value) {
            return None;
        }
        let (start, stop) = (value(cx, iv.start)?, value(cx, iv.stop)?);
        let ends = [
            value(cx, iv.left.p)?,
            value(cx, iv.left.q)?,
            value(cx, iv.right.p)?,
            value(cx, iv.right.q)?,
        ];
        if stop <= start || start.is_nan() || ends.iter().any(|v| !v.is_finite()) {
            return None;
        }
        data.push(AxisData {
            var: p.vars[j],
            start: iv.start,
            stop: iv.stop,
            ends: [iv.left.p, iv.left.q, iv.right.p, iv.right.q],
            bounds: (start, stop),
            a: ca,
            b: cb,
        });
    }
    let symbols: Vec<crate::graph::SymbolId> = space
        .iter()
        .map(|&j| cx.graph.symbol_of(p.vars[j]))
        .collect::<Option<_>>()?;
    let mids: Vec<f64> = data
        .iter()
        .map(|d| f64::midpoint(d.bounds.0, d.bounds.1))
        .collect();
    let dimension = space.len();
    // The zero-order term: c = C + Σ_j g_j(x_j) with g_j(mid_j) = 0.
    let base = at_point(cx, cc, &symbols, &mids)?;
    let single = dimension == 1;
    let leftover = if single { 0.0 } else { base };
    if !single {
        for trial in 0..3_u32 {
            let point: Vec<f64> = data
                .iter()
                .enumerate()
                .map(|(k, d)| {
                    d.bounds.0
                        + (d.bounds.1 - d.bounds.0)
                            * (0.17 + 0.29 * f64::from(trial) + 0.11 * k as f64)
                            % 1.0
                })
                .collect();
            let mut predicted = base;
            for k in 0..dimension {
                let mut shifted = mids.clone();
                shifted[k] = point[k];
                predicted += at_point(cx, cc, &symbols, &shifted)? - base;
            }
            let actual = at_point(cx, cc, &symbols, &point)?;
            if (predicted - actual).abs() > 1e-8 * (1.0 + actual.abs()) {
                return None;
            }
        }
    }
    let number = |cx: &mut Cx<'_>, v: f64| cx.graph.num(Number::from(v));
    let mut axes = Vec::new();
    let mut spatial: Vec<(Vec<u32>, NodeId)> = Vec::new();
    for (pos, d) in data.iter().enumerate() {
        let (start, stop) = d.bounds;
        let (mid, half) = (f64::midpoint(start, stop), 0.5 * (stop - start));
        let xs: Vec<f64> = (0..=SAMPLES)
            .map(|k| mid + half * (PI * k as f64 / SAMPLES as f64).cos())
            .collect();
        let along = |x_value: f64| -> Vec<f64> {
            let mut point = mids.clone();
            point[pos] = x_value;
            point
        };
        let kappa = if at_point(cx, d.a, &symbols, &along(mid))? > 0.0 {
            1.0
        } else {
            -1.0
        };
        let a_hat = |cx: &Cx<'_>, x_value: f64| {
            at_point(cx, d.a, &symbols, &along(x_value)).map(|v| kappa * v)
        };
        let b_hat = |cx: &Cx<'_>, x_value: f64| {
            at_point(cx, d.b, &symbols, &along(x_value)).map(|v| kappa * v)
        };
        let c_hat = |cx: &Cx<'_>, x_value: f64| -> Option<f64> {
            let full = at_point(cx, cc, &symbols, &along(x_value))?;
            Some(kappa * if single { full } else { full - base })
        };
        let (mut ps, mut qs, mut ws) = (Vec::new(), Vec::new(), Vec::new());
        for &xv in &xs {
            let a_value = a_hat(cx, xv)?;
            if a_value <= 0.0 {
                return None;
            }
            // ∫_{x₀}^{x} b̂/â by Simpson's rule.
            let steps = 400_usize;
            let h = (xv - start) / steps as f64;
            let mut integral = 0.0;
            if h > 0.0 {
                for k in 0..=steps {
                    let s = start + h * k as f64;
                    let g = b_hat(cx, s)? / a_hat(cx, s)?;
                    let weight = if k == 0 || k == steps {
                        1.0
                    } else if k % 2 == 1 {
                        4.0
                    } else {
                        2.0
                    };
                    integral += weight * g;
                }
                integral *= h / 3.0;
            }
            let p_value = integral.exp();
            ps.push(p_value);
            ws.push(p_value / a_value);
            qs.push(-p_value * c_hat(cx, xv)? / a_value);
        }
        let mut params = vec![
            d.start,
            d.stop,
            d.ends[0],
            d.ends[1],
            d.ends[2],
            d.ends[3],
            cx.graph.int(i64::try_from(SAMPLES).ok()?),
        ];
        for list in [&ps, &qs, &ws] {
            for &v in list {
                params.push(number(cx, v));
            }
        }
        let kappa_node = cx.graph.int(if kappa > 0.0 { 1 } else { -1 });
        let mut halves = vec![0; dimension];
        halves[pos] = 1;
        spatial.push((halves, kappa_node));
        axes.push(Axis {
            var: d.var,
            shape: Shape::General {
                start: d.start,
                stop: d.stop,
                params,
            },
        });
    }
    if leftover.abs() > 1e-14 {
        spatial.push((vec![0; dimension], number(cx, leftover)));
    }
    let order = if is_zero_number(cx.graph, a2) {
        u32::from(!is_zero_number(cx.graph, a1))
    } else {
        2
    };
    let time = Time {
        var: p.vars[time_index],
        order,
        a2,
        a1,
    };
    let mut items = time_items(cx, p, conditions, Some(time_index), order, n)?;
    if !p.homogeneous(cx.graph) {
        items.push(Item {
            role: Role::Source,
            h: p.source,
            special: None,
        });
    }
    let engine = Engine {
        axes,
        time: Some(time),
        symbol: Symbol::Cartesian(spatial),
        face: None,
        exterior: None,
    };
    finish(cx, p, conditions, &engine, &items, None, None)
}
