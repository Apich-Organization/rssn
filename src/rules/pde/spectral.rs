#![allow(dead_code)]
//! Eigenfunction expansions on separable domains.
//!
//! One engine serves every bounded domain whose Laplacian separates: the
//! domain is a list of *axes*, each with a family of eigenfunctions, a
//! weight and a norm, and the solution is the (nested) series
//! `Σ T_n(t) Φ_n(x)` with `Φ_n` the product of the axis eigenfunctions.
//! Axes are
//!
//! * an interval `[a, b]` with Dirichlet, Neumann, Robin or periodic ends
//!   (sines, cosines, `sl_root` eigenfunctions),
//! * a periodic angle (`cos nθ`, `sin nθ`),
//! * a polar angle with regularity at the poles (`P_l(cos θ)`),
//! * a radius `0 ≤ r ≤ R` in a disk (`J_ν(κ r)`) or a ball
//!   (`j_l(κ r)`), the order following the preceding angular axis,
//!
//! so that boxes in one to three dimensions (Cartesian), disks, annuli-free
//! cylinders, balls and spheres are all instances. The operator's symbol,
//! the time dependence (`a₂ T'' + a₁ T' + s T = −f`, solved in closed form
//! with Duhamel's integral for time-dependent forcing), nonhomogeneous
//! boundary data (as boundary terms of Green's second identity), sources
//! and initial data are handled once, here.
//!
//! The coefficients are computed exactly by projection; indices where the
//! symbolic projection differs from the projection with a literal index
//! (zero modes, resonances) are written out, and data that is a finite sum
//! of modes gives a finite sum.

use super::Condition;
use super::Conditions;
use super::Problem;
use super::agree;
use super::dummy;
use super::satisfies;
use super::util::call;
use super::util::contains_op;
use super::util::defint;
use super::util::div;
use super::util::fraction;
use super::util::infinity;
use super::util::is_zero_number;
use super::util::pi;
use super::util::sample;
use super::util::sample_env;
use super::util::series;
use super::util::sqrt;
use super::verified;
use crate::graph::Cx;
use crate::graph::Facts;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::rules::calculus::derivative;
use crate::rules::complex::build::add;
use crate::rules::complex::build::mul;
use crate::rules::complex::build::neg;
use crate::rules::complex::build::powi;
use crate::rules::complex::build::sub;

/// The kind of boundary condition at one end.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub(super) enum Kind {
    Dirichlet,
    Neumann,
    Robin,
}

/// `p u + q u_x = value` at `point`.
#[derive(Copy, Clone, Debug)]
pub(super) struct End {
    pub kind: Kind,
    pub point: NodeId,
    pub p: NodeId,
    pub q: NodeId,
    pub value: NodeId,
}

impl End {
    pub(super) fn from_condition(
        cx: &mut Cx<'_>,
        p: &Problem,
        c: &Condition,
    ) -> Option<Self> {
        let (zero, one) = (cx.graph.int(0), cx.graph.int(1));
        let first = p.unit(c.on, 1);
        let (kind, a, b) = if c.derivative.iter().all(|&d| d == 0) {
            (Kind::Dirichlet, one, zero)
        } else if c.derivative == first {
            match c.robin {
                | None => (Kind::Neumann, zero, one),
                | Some(h) => (Kind::Robin, h, one),
            }
        } else {
            return None;
        };
        Some(Self { kind, point: c.point, p: a, q: b, value: c.value })
    }
}

/// An interval `[start, stop]` with conditions at both ends, or a period.
#[derive(Clone, Debug)]
pub(super) struct Interval {
    pub var: NodeId,
    pub start: NodeId,
    pub stop: NodeId,
    pub length: NodeId,
    pub left: End,
    pub right: End,
    pub periodic: bool,
}

#[derive(Clone, Debug)]
pub(super) enum Shape {
    Interval(Interval),
    /// A periodic angle `[0, 2π)`.
    Angle,
    /// A polar angle `[0, π]`: Legendre functions of `cos θ`.
    Polar,
    /// A radius `[0, R]`; `sphere` for a ball (weight `r²`); `has_angle`
    /// when the Bessel order comes from the previous axis.
    Radial { radius: NodeId, end: End, sphere: bool, has_angle: bool },
}

#[derive(Clone, Debug)]
pub(super) struct Axis {
    pub var: NodeId,
    pub shape: Shape,
}

/// One family of eigenfunctions along an axis, for a given index.
#[derive(Clone, Debug)]
struct Family {
    mode: NodeId,
    /// The eigenvalue contributed to `-Δ`.
    eigen: NodeId,
    /// First index of the family.
    first: i64,
    /// How many leading indices are always written out.
    forced: i64,
    /// `∫ w φ² dx` when known in closed form for this index.
    norm: Option<NodeId>,
}

#[derive(Clone, Debug)]
struct Step {
    index: NodeId,
    family: Family,
}

/// Number of low indices compared between the symbolic and the literal
/// projection.
const WINDOW: i64 = 8;
/// Highest degree of a Legendre expansion.
const LEGENDRE_MAX: i64 = 10;

impl Axis {
    fn domain(
        &self,
        cx: &mut Cx<'_>,
    ) -> Option<(NodeId, NodeId)> {
        let zero = cx.graph.int(0);
        match &self.shape {
            | Shape::Interval(iv) => Some((iv.start, iv.stop)),
            | Shape::Angle => {
                let two = cx.graph.int(2);
                let p = pi(cx)?;
                Some((zero, mul(cx.graph, &[two, p])))
            },
            | Shape::Polar => Some((zero, pi(cx)?)),
            | Shape::Radial { radius, .. } => Some((zero, *radius)),
        }
    }

    fn weight(
        &self,
        cx: &mut Cx<'_>,
    ) -> Option<NodeId> {
        Some(match &self.shape {
            | Shape::Interval(_) | Shape::Angle => cx.graph.int(1),
            | Shape::Polar => call(cx, "sin", &[self.var])?,
            | Shape::Radial { sphere, .. } => powi(cx.graph, self.var, if *sphere { 2 } else { 1 }),
        })
    }

    /// Weight at a boundary point (for Green's identity).
    fn boundary_weight(
        &self,
        cx: &mut Cx<'_>,
        point: NodeId,
    ) -> NodeId {
        match &self.shape {
            | Shape::Radial { sphere, .. } => powi(cx.graph, point, if *sphere { 2 } else { 1 }),
            | _ => cx.graph.int(1),
        }
    }

    const fn finite_only(&self) -> bool {
        matches!(self.shape, Shape::Polar)
    }

    fn families(
        &self,
        cx: &mut Cx<'_>,
        path: &[Step],
        n: NodeId,
    ) -> Option<Vec<Family>> {
        let numeric = cx.graph.number_of(n).is_some();
        let zero = cx.graph.int(0);
        match &self.shape {
            | Shape::Interval(iv) => interval_families(cx, iv, n, numeric),
            | Shape::Angle => {
                let theta = self.var;
                let n_theta = mul(cx.graph, &[n, theta]);
                let cos = call(cx, "cos", &[n_theta])?;
                let sin = call(cx, "sin", &[n_theta])?;
                let norm = if numeric { None } else { Some(pi(cx)?) };
                Some(vec![
                    Family { mode: cos, eigen: zero, first: 0, forced: 1, norm },
                    Family { mode: sin, eigen: zero, first: 1, forced: 0, norm },
                ])
            },
            | Shape::Polar => {
                let c = call(cx, "cos", &[self.var])?;
                let mode = call(cx, "legendre", &[n, c])?;
                let two = cx.graph.int(2);
                let two_n = mul(cx.graph, &[two, n]);
                let one = cx.graph.int(1);
                let odd = add(cx.graph, &[two_n, one]);
                let norm = if numeric { None } else { Some(div(cx, two, odd)) };
                Some(vec![Family { mode, eigen: zero, first: 0, forced: 1, norm }])
            },
            | Shape::Radial { radius, end, sphere, has_angle } => {
                let order = if *has_angle { path.last()?.index } else { cx.graph.int(0) };
                radial_families(cx, self.var, *radius, end, *sphere, order, n, numeric)
            },
        }
    }
}

fn interval_families(
    cx: &mut Cx<'_>,
    iv: &Interval,
    n: NodeId,
    numeric: bool,
) -> Option<Vec<Family>> {
    let pi = pi(cx)?;
    let l = iv.length;
    let inverse_l = powi(cx.graph, l, -1);
    let xi = if is_zero_number(cx.graph, iv.start) { iv.var } else { sub(cx.graph, iv.var, iv.start) };
    let half_l = {
        let half = fraction(cx, 1, 2)?;
        mul(cx.graph, &[half, l])
    };
    let norm = if numeric { None } else { Some(half_l) };
    let trig = |cx: &mut Cx<'_>, name: &str, k: NodeId| -> Option<NodeId> {
        let kx = mul(cx.graph, &[k, xi]);
        call(cx, name, &[kx])
    };
    let eigen = |cx: &mut Cx<'_>, k: NodeId| {
        let e = powi(cx.graph, k, 2);
        cx.simplify(e)
    };
    let (left, right) = (iv.left.kind, iv.right.kind);
    if iv.periodic {
        let two = cx.graph.int(2);
        let k = mul(cx.graph, &[two, n, pi, inverse_l]);
        let e = eigen(cx, k);
        let (c, s) = (trig(cx, "cos", k)?, trig(cx, "sin", k)?);
        return Some(vec![
            Family { mode: c, eigen: e, first: 0, forced: 1, norm },
            Family { mode: s, eigen: e, first: 1, forced: 0, norm },
        ]);
    }
    let integer_k = mul(cx.graph, &[n, pi, inverse_l]);
    let half_integer_k = {
        let minus_half = fraction(cx, -1, 2)?;
        let shifted = add(cx.graph, &[n, minus_half]);
        mul(cx.graph, &[shifted, pi, inverse_l])
    };
    let one = |name: &str, k: NodeId, first: i64, forced: i64, cx: &mut Cx<'_>| -> Option<Vec<Family>> {
        let mode = trig(cx, name, k)?;
        let e = eigen(cx, k);
        Some(vec![Family { mode, eigen: e, first, forced, norm }])
    };
    match (left, right) {
        | (Kind::Dirichlet, Kind::Dirichlet) => one("sin", integer_k, 1, 0, cx),
        | (Kind::Neumann, Kind::Neumann) => one("cos", integer_k, 0, 1, cx),
        | (Kind::Dirichlet, Kind::Neumann) => one("sin", half_integer_k, 1, 0, cx),
        | (Kind::Neumann, Kind::Dirichlet) => one("cos", half_integer_k, 1, 0, cx),
        | _ => {
            // A Robin end: X = q0 k cos(kξ) − p0 sin(kξ).
            let k = call(cx, "sl_root", &[iv.left.p, iv.left.q, iv.right.p, iv.right.q, l, n])?;
            let cos = trig(cx, "cos", k)?;
            let sin = trig(cx, "sin", k)?;
            let a = mul(cx.graph, &[iv.left.q, k, cos]);
            let b = mul(cx.graph, &[iv.left.p, sin]);
            let mode = sub(cx.graph, a, b);
            let e = eigen(cx, k);
            Some(vec![Family { mode, eigen: e, first: 1, forced: 0, norm: None }])
        },
    }
}

#[allow(clippy::too_many_arguments)]
fn radial_families(
    cx: &mut Cx<'_>,
    r: NodeId,
    radius: NodeId,
    end: &End,
    sphere: bool,
    order: NodeId,
    n: NodeId,
    numeric: bool,
) -> Option<Vec<Family>> {
    let zero = cx.graph.int(0);
    let half = fraction(cx, 1, 2)?;
    // The Bessel order: n for a disk, l + 1/2 for a ball.
    let nu = if sphere { add(cx.graph, &[order, half]) } else { order };
    let nu = cx.simplify(nu);
    let order_value = cx.graph.number_of(nu).map(Number::to_f64);
    let neumann_like = is_zero_number(cx.graph, end.p);
    let constant_mode = neumann_like && order_value.is_some_and(|v| (v - if sphere { 0.5 } else { 0.0 }).abs() < 1e-12);
    if numeric && constant_mode && cx.graph.number_of(n).is_some_and(Number::is_zero) {
        let one = cx.graph.int(1);
        return Some(vec![Family { mode: one, eigen: zero, first: 0, forced: 1, norm: None }]);
    }
    // h R, with the sphere's shift.
    let hr = if neumann_like {
        zero
    } else {
        let ratio = div(cx, end.p, end.q);
        mul(cx.graph, &[ratio, radius])
    };
    let shifted = if sphere {
        let minus_half = fraction(cx, -1, 2)?;
        add(cx.graph, &[hr, minus_half])
    } else {
        hr
    };
    let root = match end.kind {
        | Kind::Dirichlet => call(cx, "bessel_zero", &[nu, n])?,
        | _ => call(cx, "bessel_root", &[nu, shifted, n])?,
    };
    let inverse_radius = powi(cx.graph, radius, -1);
    let kappa = mul(cx.graph, &[root, inverse_radius]);
    let eigen = powi(cx.graph, kappa, 2);
    let first = i64::from(!constant_mode);
    let kr = mul(cx.graph, &[kappa, r]);
    let sin_form = sphere && order_value.is_some_and(|v| (v - 0.5).abs() < 1e-12);
    let forced = i64::from(constant_mode);
    if sin_form {
        let s = call(cx, "sin", &[kr])?;
        let mode = div(cx, s, r);
        // ∫ sin²(κ r) dr = R/2 − sin(2 κ R)/(4 κ)
        let two_kr = {
            let two = cx.graph.int(2);
            mul(cx.graph, &[two, kappa, radius])
        };
        let s2 = call(cx, "sin", &[two_kr])?;
        let four = cx.graph.int(4);
        let denominator = mul(cx.graph, &[four, kappa]);
        let correction = div(cx, s2, denominator);
        let half_r = mul(cx.graph, &[half, radius]);
        let norm = sub(cx.graph, half_r, correction);
        let norm = if numeric { None } else { Some(norm) };
        return Some(vec![Family { mode, eigen, first, forced, norm }]);
    }
    let j = call(cx, "besselj", &[nu, kr])?;
    let mode = if sphere {
        let minus_half = fraction(cx, -1, 2)?;
        let power = crate::rules::complex::build::pow(cx.graph, r, minus_half);
        mul(cx.graph, &[power, j])
    } else {
        j
    };
    // ∫ r J_ν(κ r)² dr = R²/2 [J_ν'(z)² + (1 − ν²/z²) J_ν(z)²].
    let one = cx.graph.int(1);
    let (j_z, plus) = {
        let jz = call(cx, "besselj", &[nu, root])?;
        let up = add(cx.graph, &[nu, one]);
        (jz, call(cx, "besselj", &[up, root])?)
    };
    let r2_half = {
        let rr = powi(cx.graph, radius, 2);
        mul(cx.graph, &[half, rr])
    };
    let norm = if matches!(end.kind, Kind::Dirichlet) {
        let sq = powi(cx.graph, plus, 2);
        mul(cx.graph, &[r2_half, sq])
    } else {
        let down = {
            let m = neg(cx.graph, one);
            let lower = add(cx.graph, &[nu, m]);
            call(cx, "besselj", &[lower, root])?
        };
        let derivative_z = {
            let difference = sub(cx.graph, down, plus);
            mul(cx.graph, &[half, difference])
        };
        let d2 = powi(cx.graph, derivative_z, 2);
        let nu2 = powi(cx.graph, nu, 2);
        let z2 = powi(cx.graph, root, 2);
        let ratio = div(cx, nu2, z2);
        let factor = sub(cx.graph, one, ratio);
        let j2 = powi(cx.graph, j_z, 2);
        let second = mul(cx.graph, &[factor, j2]);
        let bracket = add(cx.graph, &[d2, second]);
        mul(cx.graph, &[r2_half, bracket])
    };
    let norm = if numeric { if constant_mode && cx.graph.number_of(n).is_some_and(Number::is_zero) { None } else { Some(norm) } } else { Some(norm) };
    Some(vec![Family { mode, eigen, first, forced, norm }])
}

// ----------------------------------------------------------------------
// The operator and the data
// ----------------------------------------------------------------------

/// The spatial symbol of the equation: how the operator acts on a mode,
/// as a function of the eigenvalues of the axes.
#[derive(Clone, Debug)]
pub(super) enum Symbol {
    /// `Σ c_α Π (−λ_j)^(α_j/2)` with `α_j / 2` per axis.
    Cartesian(Vec<(Vec<u32>, NodeId)>),
    /// `a0 − b Σ λ_j`.
    Laplacian { b: NodeId, a0: NodeId },
}

impl Symbol {
    fn value(
        &self,
        cx: &mut Cx<'_>,
        eigens: &[NodeId],
    ) -> NodeId {
        let total = match self {
            | Self::Cartesian(terms) => {
                let mut sum = Vec::new();
                for (halves, c) in terms {
                    let mut factors = vec![*c];
                    for (&h, &lambda) in halves.iter().zip(eigens) {
                        if h > 0 {
                            let minus = neg(cx.graph, lambda);
                            factors.push(powi(cx.graph, minus, i64::from(h)));
                        }
                    }
                    sum.push(mul(cx.graph, &factors));
                }
                add(cx.graph, &sum)
            },
            | Self::Laplacian { b, a0 } => {
                let lambdas = add(cx.graph, eigens);
                let bl = mul(cx.graph, &[*b, lambdas]);
                sub(cx.graph, *a0, bl)
            },
        };
        cx.simplify(total)
    }
}

#[derive(Copy, Clone, Debug)]
pub(super) struct Time {
    pub var: NodeId,
    pub order: u32,
    pub a2: NodeId,
    pub a1: NodeId,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub(super) enum Role {
    Displacement,
    Velocity,
    Source,
}

#[derive(Copy, Clone, Debug)]
pub(super) struct Special {
    pub axis: usize,
    pub kind: Kind,
    pub point: NodeId,
    pub weight: NodeId,
}

/// A function waiting to be projected: initial data, forcing, or boundary
/// data (projected by evaluation at its end on `special.axis`).
#[derive(Copy, Clone, Debug)]
pub(super) struct Item {
    pub role: Role,
    pub h: NodeId,
    pub special: Option<Special>,
}

pub(super) struct Engine {
    pub axes: Vec<Axis>,
    pub time: Option<Time>,
    pub symbol: Symbol,
}

fn zero_node(cx: &mut Cx<'_>) -> NodeId {
    cx.graph.int(0)
}

fn exp_of(
    cx: &mut Cx<'_>,
    x: NodeId,
) -> Option<NodeId> {
    call(cx, "exp", &[x])
}

/// The time dependence of one mode: `a₂ T'' + a₁ T' + s T + f = 0` with
/// `T(0) = c0`, `T'(0) = c1`; for steady problems `T = −f/s`.
#[allow(clippy::too_many_arguments, clippy::many_single_char_names)]
fn time_factor(
    cx: &mut Cx<'_>,
    p: &Problem,
    engine: &Engine,
    s: NodeId,
    c0: NodeId,
    c1: NodeId,
    f: NodeId,
) -> Option<NodeId> {
    let zero = zero_node(cx);
    let Some(time) = engine.time else {
        if cx.is_zero(f) {
            return Some(zero);
        }
        if cx.is_zero(s) {
            return None;
        }
        let ratio = div(cx, f, s);
        let t = neg(cx.graph, ratio);
        return Some(cx.simplify(t));
    };
    let t = time.var;
    let t_symbol = cx.graph.symbol_of(t)?;
    let forced = !cx.is_zero(f);
    let f_constant = !cx.graph.depends_on(cx.graph.find(f), t_symbol);
    let mut terms = Vec::new();
    match time.order {
        | 1 => {
            let rate = div(cx, s, time.a1);
            let minus_rate = neg(cx.graph, rate);
            let decay_exponent = mul(cx.graph, &[minus_rate, t]);
            let decay = exp_of(cx, decay_exponent)?;
            if !cx.is_zero(c0) {
                terms.push(mul(cx.graph, &[c0, decay]));
            }
            if forced {
                if f_constant {
                    if cx.is_zero(s) {
                        let ft = mul(cx.graph, &[f, t]);
                        let r = div(cx, ft, time.a1);
                        terms.push(neg(cx.graph, r));
                    } else {
                        let ratio = div(cx, f, s);
                        let one = cx.graph.int(1);
                        let gap = sub(cx.graph, one, decay);
                        let m = mul(cx.graph, &[ratio, gap]);
                        terms.push(neg(cx.graph, m));
                    }
                } else {
                    let (tau, _) = dummy(cx, p, "tau");
                    let f_tau = cx.graph.substitute(f, t, tau);
                    let elapsed = sub(cx.graph, t, tau);
                    let exponent = mul(cx.graph, &[minus_rate, elapsed]);
                    let kernel = exp_of(cx, exponent)?;
                    let body = mul(cx.graph, &[f_tau, kernel]);
                    let body = div(cx, body, time.a1);
                    let body = neg(cx.graph, body);
                    let zero_t = cx.graph.int(0);
                    terms.push(defint(cx, body, tau, zero_t, t)?);
                }
            }
        },
        | 2 => {
            let two = cx.graph.int(2);
            let two_a2 = mul(cx.graph, &[two, time.a2]);
            let gamma = div(cx, time.a1, two_a2);
            let gamma = cx.simplify(gamma);
            let ratio = div(cx, s, time.a2);
            let g2 = powi(cx.graph, gamma, 2);
            let omega2 = sub(cx.graph, ratio, g2);
            let omega2 = cx.simplify(omega2);
            let undamped = cx.is_zero(gamma);
            let damp = if undamped {
                cx.graph.int(1)
            } else {
                let minus_gamma = neg(cx.graph, gamma);
                let e = mul(cx.graph, &[minus_gamma, t]);
                exp_of(cx, e)?
            };
            if cx.is_zero(omega2) {
                // Critical: T = e^{-γt}(c0 + (c1 + γ c0) t) + forcing.
                if !undamped && forced {
                    return None;
                }
                let gc0 = mul(cx.graph, &[gamma, c0]);
                let slope = add(cx.graph, &[c1, gc0]);
                let line = {
                    let st = mul(cx.graph, &[slope, t]);
                    add(cx.graph, &[c0, st])
                };
                terms.push(mul(cx.graph, &[damp, line]));
                if forced {
                    if !f_constant {
                        let (tau, _) = dummy(cx, p, "tau");
                        let f_tau = cx.graph.substitute(f, t, tau);
                        let elapsed = sub(cx.graph, t, tau);
                        let body = mul(cx.graph, &[elapsed, f_tau]);
                        let body = div(cx, body, time.a2);
                        let body = neg(cx.graph, body);
                        let zero_t = cx.graph.int(0);
                        terms.push(defint(cx, body, tau, zero_t, t)?);
                    } else {
                        let t2 = powi(cx.graph, t, 2);
                        let half = fraction(cx, 1, 2)?;
                        let q = mul(cx.graph, &[half, f, t2]);
                        let q = div(cx, q, time.a2);
                        terms.push(neg(cx.graph, q));
                    }
                }
            } else {
                let omega = sqrt(cx, omega2)?;
                let omega_t = mul(cx.graph, &[omega, t]);
                let cos = call(cx, "cos", &[omega_t])?;
                let sin = call(cx, "sin", &[omega_t])?;
                let sin_over = div(cx, sin, omega);
                let c_part = mul(cx.graph, &[c0, cos]);
                let gc0 = mul(cx.graph, &[gamma, c0]);
                let slope = add(cx.graph, &[c1, gc0]);
                let s_part = mul(cx.graph, &[slope, sin_over]);
                let both = add(cx.graph, &[c_part, s_part]);
                terms.push(mul(cx.graph, &[damp, both]));
                if forced {
                    if f_constant {
                        if cx.is_zero(s) {
                            return None;
                        }
                        let ratio = div(cx, f, s);
                        let gs = mul(cx.graph, &[gamma, sin_over]);
                        let inside = add(cx.graph, &[cos, gs]);
                        let shape = mul(cx.graph, &[damp, inside]);
                        let one = cx.graph.int(1);
                        let gap = sub(cx.graph, one, shape);
                        let m = mul(cx.graph, &[ratio, gap]);
                        terms.push(neg(cx.graph, m));
                    } else {
                        let (tau, _) = dummy(cx, p, "tau");
                        let f_tau = cx.graph.substitute(f, t, tau);
                        let elapsed = sub(cx.graph, t, tau);
                        let wave_arg = mul(cx.graph, &[omega, elapsed]);
                        let wave = call(cx, "sin", &[wave_arg])?;
                        let damping = if undamped {
                            cx.graph.int(1)
                        } else {
                            let minus_gamma = neg(cx.graph, gamma);
                            let e = mul(cx.graph, &[minus_gamma, elapsed]);
                            exp_of(cx, e)?
                        };
                        let body = mul(cx.graph, &[damping, wave, f_tau]);
                        let scale = mul(cx.graph, &[time.a2, omega]);
                        let body = div(cx, body, scale);
                        let body = neg(cx.graph, body);
                        let zero_t = cx.graph.int(0);
                        terms.push(defint(cx, body, tau, zero_t, t)?);
                    }
                }
            }
        },
        | _ => return None,
    }
    let total = add(cx.graph, &terms);
    Some(cx.simplify(total))
}

// ----------------------------------------------------------------------
// The recursion over the axes
// ----------------------------------------------------------------------

fn is_zero_item(
    cx: &mut Cx<'_>,
    item: &Item,
) -> bool {
    cx.is_zero(item.h)
}

/// The projection of `item` on the family `fam` of axis number `level`.
fn project(
    cx: &mut Cx<'_>,
    engine: &Engine,
    level: usize,
    fam: &Family,
    item: &Item,
) -> Option<Item> {
    let axis = engine.axes.get(level)?;
    if is_zero_item(cx, item) {
        return Some(*item);
    }
    let norm = match fam.norm {
        | Some(n) => n,
        | None => {
            let (lo, hi) = axis.domain(cx)?;
            let w = axis.weight(cx)?;
            let square = powi(cx.graph, fam.mode, 2);
            let body = mul(cx.graph, &[w, square]);
            let integral = defint(cx, body, axis.var, lo, hi)?;
            cx.simplify(integral)
        },
    };
    let h = match item.special {
        | Some(sp) if sp.axis == level => {
            let basis = if sp.kind == Kind::Dirichlet {
                let d = derivative(cx.graph, fam.mode, axis.var)?;
                neg(cx.graph, d)
            } else {
                fam.mode
            };
            let at_end = cx.graph.substitute(basis, axis.var, sp.point);
            let factor = mul(cx.graph, &[sp.weight, at_end]);
            let c = mul(cx.graph, &[item.h, factor]);
            let c = div(cx, c, norm);
            return Some(Item { role: item.role, h: cx.simplify(c), special: None });
        },
        | _ => {
            let (lo, hi) = axis.domain(cx)?;
            let w = axis.weight(cx)?;
            let integrand = mul(cx.graph, &[item.h, fam.mode, w]);
            let integral = defint(cx, integrand, axis.var, lo, hi)?;
            div(cx, integral, norm)
        },
    };
    Some(Item { role: item.role, h: cx.simplify(h), special: item.special })
}

/// Whether the expression `h` (in the index `n`) vanishes, to rounding, at
/// the integers from `from` on.
fn vanishes_on_integers(
    cx: &mut Cx<'_>,
    h: NodeId,
    n: NodeId,
    from: i64,
) -> bool {
    let Some(n_symbol) = cx.graph.symbol_of(n) else {
        return false;
    };
    for k in 0..2_u32 {
        for m in from..from + 14 {
            let mut env = sample_env(cx.graph, h, k);
            #[allow(clippy::cast_precision_loss)]
            env.bind(n_symbol, m as f64);
            match cx.graph.eval(h, &env) {
                | Some(v) if v.is_finite() && v.abs() < 1e-9 => {},
                | _ => return false,
            }
        }
    }
    true
}

fn expand(
    cx: &mut Cx<'_>,
    p: &Problem,
    engine: &Engine,
    level: usize,
    items: &[Item],
    path: &mut Vec<Step>,
) -> Option<NodeId> {
    if level == engine.axes.len() {
        return leaf(cx, p, engine, items, path);
    }
    let axis = engine.axes.get(level)?;
    let stems = ["n", "m", "l", "q"];
    let (n_node, n_symbol) = dummy(cx, p, stems.get(level).copied().unwrap_or("q"));
    cx.graph.assume(n_symbol, Facts::INTEGER | Facts::NONNEGATIVE);
    let general_families = axis.families(cx, path, n_node)?;
    let zero = zero_node(cx);
    let mut terms = Vec::new();
    for (fi, gfam) in general_families.iter().enumerate() {
        let first = gfam.first;
        if axis.finite_only() {
            // Explicit modes up to a maximal degree; the data must be exhausted.
            let mut tail_zero = true;
            for m in first..=LEGENDRE_MAX {
                let number = cx.graph.int(m);
                let fam_m = axis.families(cx, path, number)?.get(fi)?.clone();
                let projected: Vec<Item> = items.iter().map(|it| project(cx, engine, level, &fam_m, it)).collect::<Option<_>>()?;
                if projected.iter().any(|it| contains_op(cx.graph, it.h, "defint")) {
                    return None;
                }
                let live = projected.iter().any(|it| !is_zero_item(cx, it));
                if m == LEGENDRE_MAX {
                    tail_zero = !live;
                }
                if live {
                    path.push(Step { index: number, family: fam_m });
                    let body = expand(cx, p, engine, level + 1, &projected, path);
                    path.pop();
                    terms.push(body?);
                }
            }
            if !tail_zero {
                return None;
            }
            continue;
        }
        let general: Vec<Item> = items.iter().map(|it| project(cx, engine, level, gfam, it)).collect::<Option<_>>()?;
        let unresolved = general.iter().any(|it| contains_op(cx.graph, it.h, "defint"));
        let mut general_zero = general.iter().all(|it| is_zero_item(cx, it));
        if !general_zero && !unresolved && general.iter().all(|it| vanishes_on_integers(cx, it.h, n_node, first + gfam.forced + WINDOW)) {
            general_zero = true;
        }
        // Literal indices.
        let mut exact: Vec<(i64, NodeId, Family, Vec<Item>)> = Vec::new();
        let mut last_exceptional: Option<i64> = None;
        let window = if unresolved { gfam.forced } else { WINDOW };
        for m in first..first + window {
            let number = cx.graph.int(m);
            let fam_m = axis.families(cx, path, number)?.get(fi)?.clone();
            let projected: Vec<Item> = items.iter().map(|it| project(cx, engine, level, &fam_m, it)).collect::<Option<_>>()?;
            let mut differs = m < first + gfam.forced;
            if !differs && !unresolved {
                for (g, e) in general.iter().zip(&projected) {
                    if general_zero {
                        let live = !is_zero_item(cx, e) && sample(cx.graph, e.h, 0).is_none_or(|v| v.abs() > 1e-9);
                        differs |= live;
                    } else {
                        let substituted = cx.graph.substitute(g.h, n_node, number);
                        let substituted = cx.simplify(substituted);
                        differs |= !agree(cx, substituted, e.h);
                    }
                }
            }
            if differs {
                last_exceptional = Some(m);
            }
            exact.push((m, number, fam_m, projected));
        }
        for (m, number, fam_m, projected) in exact {
            if last_exceptional.is_none_or(|last| m > last) {
                break;
            }
            if projected.iter().all(|it| is_zero_item(cx, it)) {
                continue;
            }
            path.push(Step { index: number, family: fam_m });
            let body = expand(cx, p, engine, level + 1, &projected, path);
            path.pop();
            terms.push(body?);
        }
        if !general_zero {
            let start = last_exceptional.map_or(first, |m| m + 1);
            path.push(Step { index: n_node, family: gfam.clone() });
            let body = expand(cx, p, engine, level + 1, &general, path);
            path.pop();
            let body = body?;
            if !cx.is_zero(body) {
                let start = cx.graph.int(start);
                let oo = infinity(cx)?;
                terms.push(series(cx, body, n_node, start, oo)?);
            }
        }
    }
    if terms.is_empty() {
        return Some(zero);
    }
    let total = add(cx.graph, &terms);
    Some(cx.simplify(total))
}

fn leaf(
    cx: &mut Cx<'_>,
    p: &Problem,
    engine: &Engine,
    items: &[Item],
    path: &[Step],
) -> Option<NodeId> {
    let eigens: Vec<NodeId> = path.iter().map(|s| s.family.eigen).collect();
    let s = engine.symbol.value(cx, &eigens);
    let zero = zero_node(cx);
    let (mut c0, mut c1) = (zero, zero);
    let mut forcing = Vec::new();
    for item in items {
        match item.role {
            | Role::Displacement => c0 = item.h,
            | Role::Velocity => c1 = item.h,
            | Role::Source => forcing.push(item.h),
        }
    }
    let f = add(cx.graph, &forcing);
    let f = cx.simplify(f);
    let t = time_factor(cx, p, engine, s, c0, c1, f)?;
    if cx.is_zero(t) {
        return Some(zero);
    }
    let mut factors = vec![t];
    factors.extend(path.iter().map(|s| s.family.mode));
    let term = mul(cx.graph, &factors);
    Some(cx.simplify(term))
}

/// Runs the expansion for the initial data and forcing in `items`.
fn run(
    cx: &mut Cx<'_>,
    p: &Problem,
    engine: &Engine,
    items: &[Item],
) -> Option<NodeId> {
    let mut path = Vec::new();
    expand(cx, p, engine, 0, items, &mut path)
}

// ----------------------------------------------------------------------
// Cartesian boxes
// ----------------------------------------------------------------------

/// Orders two end points: `(start, stop)` by numeric value.
fn order_points(
    cx: &mut Cx<'_>,
    a: NodeId,
    b: NodeId,
) -> Option<(NodeId, NodeId, bool)> {
    if cx.is_zero(a) {
        return Some((a, b, false));
    }
    if cx.is_zero(b) {
        return Some((b, a, true));
    }
    let (va, vb) = (sample(cx.graph, a, 0)?, sample(cx.graph, b, 0)?);
    if va <= vb { Some((a, b, false)) } else { Some((b, a, true)) }
}

/// The interval structure of variable `j`, from its two conditions.
pub(super) fn interval(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
    j: usize,
) -> Option<Interval> {
    let on: Vec<&Condition> = conditions.0.iter().filter(|c| c.on == j).collect();
    let var = p.vars[j];
    let zero_index = vec![0; p.dimension()];
    // Periodic: u(a) = u(b), optionally with the same for the derivative.
    let periodic = on.iter().find(|c| {
        c.derivative == zero_index
            && cx.graph.op(c.value) == crate::graph::op::core::APPLY
            && cx.graph.children(c.value).first() == Some(&p.function)
    });
    if let Some(c) = periodic {
        let args = cx.graph.children(c.value).to_vec();
        let other = *args.get(j + 1)?;
        let (start, stop, _) = order_points(cx, c.point, other)?;
        let length = sub(cx.graph, stop, start);
        let length = cx.simplify(length);
        let end = |point: NodeId, cx: &mut Cx<'_>| End {
            kind: Kind::Dirichlet,
            point,
            p: cx.graph.int(1),
            q: cx.graph.int(0),
            value: cx.graph.int(0),
        };
        let (left, right) = (end(start, cx), end(stop, cx));
        for c in &on {
            let ok = c.derivative == zero_index || c.derivative == p.unit(j, 1);
            if !ok {
                return None;
            }
        }
        return Some(Interval { var, start, stop, length, left, right, periodic: true });
    }
    let [a, b] = on.as_slice() else {
        return None;
    };
    let (start, stop, swapped) = order_points(cx, a.point, b.point)?;
    let (ca, cb) = if swapped { (*b, *a) } else { (*a, *b) };
    let left = End::from_condition(cx, p, ca)?;
    let right = End::from_condition(cx, p, cb)?;
    let length = sub(cx.graph, stop, start);
    let length = cx.simplify(length);
    Some(Interval { var, start, stop, length, left, right, periodic: false })
}

/// The time variable of an evolution equation, if the equation is one.
pub(super) fn evolution_time(
    cx: &mut Cx<'_>,
    p: &Problem,
) -> Option<usize> {
    let time = p.time(cx.graph);
    let named_t = cx.graph.symbol_of(p.vars[time]).is_some_and(|s| cx.graph.interner().symbol_name(s) == "t");
    let has = |order: u32| p.linear.iter().any(|(i, c)| i[time] == order && !is_zero_number(cx.graph, *c));
    if !has(1) && !has(2) {
        return None;
    }
    (named_t || (has(1) && !has(2))).then_some(time)
}

/// Whether the node contains a series or an unevaluated integral.
fn is_closed(
    cx: &Cx<'_>,
    node: NodeId,
) -> bool {
    !contains_op(cx.graph, node, "sum") && !contains_op(cx.graph, node, "defint")
}

/// Initial-value, boundary-value and eigenvalue problems for equations
/// with constant coefficients on boxes in any number of dimensions.
pub(super) fn solve_box(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
) -> Option<NodeId> {
    if p.nonlinear {
        return None;
    }
    let n = p.dimension();
    let time_index = evolution_time(cx, p);
    let space: Vec<usize> = (0..n).filter(|&j| Some(j) != time_index).collect();
    if space.is_empty() {
        return None;
    }
    // The operator: time coefficients, spatial terms.
    let (mut a2, mut a1) = (cx.graph.int(0), cx.graph.int(0));
    let mut spatial: Vec<(Vec<u32>, NodeId)> = Vec::new();
    for (index, c) in &p.linear {
        if is_zero_number(cx.graph, *c) {
            continue;
        }
        if !p.constant(cx.graph, *c) {
            return None;
        }
        if let Some(t) = time_index {
            if index[t] > 0 {
                if index.iter().enumerate().any(|(j, &d)| j != t && d != 0) {
                    return None;
                }
                match index[t] {
                    | 1 => a1 = *c,
                    | 2 => a2 = *c,
                    | _ => return None,
                }
                continue;
            }
        }
        if index.iter().any(|d| d % 2 != 0) {
            return None;
        }
        spatial.push((space.iter().map(|&j| index[j] / 2).collect(), *c));
    }
    if spatial.is_empty() {
        return None;
    }
    let order = if !is_zero_number(cx.graph, a2) {
        2
    } else if !is_zero_number(cx.graph, a1) {
        1
    } else {
        0
    };
    let time = time_index.map(|t| Time { var: p.vars[t], order, a2, a1 });
    // Axes.
    let mut axes = Vec::new();
    let mut intervals = Vec::new();
    for &j in &space {
        let iv = interval(cx, p, conditions, j)?;
        intervals.push(iv.clone());
        axes.push(Axis { var: p.vars[j], shape: Shape::Interval(iv) });
    }
    // Every condition is a boundary condition on an axis or an initial one.
    for c in &conditions.0 {
        if Some(c.on) == time_index {
            if !is_zero_number(cx.graph, c.point) {
                return None;
            }
            continue;
        }
        if !space.contains(&c.on) {
            return None;
        }
    }
    let mut items = time_items(cx, p, conditions, time_index, order, n)?;
    if !p.homogeneous(cx.graph) {
        items.push(Item { role: Role::Source, h: p.source, special: None });
    }
    // Boundary data as forcing (Green's second identity).
    for (level, iv) in intervals.iter().enumerate() {
        // Only a pure second derivative in this variable may occur.
        let mut cjj = None;
        for (halves, c) in &spatial {
            if halves.get(level).copied().unwrap_or(0) == 0 {
                continue;
            }
            let pure = halves.iter().enumerate().all(|(k, &h)| if k == level { h == 1 } else { h == 0 });
            if !pure || cjj.is_some() {
                cjj = Some(NodeId::NONE);
                break;
            }
            cjj = Some(*c);
        }
        let has_data = !iv.periodic && (!cx.is_zero(iv.left.value) || !cx.is_zero(iv.right.value));
        if !has_data {
            continue;
        }
        let cjj = cjj.filter(|c| *c != NodeId::NONE)?;
        items.extend(boundary_items(cx, iv, level, cjj)?);
    }
    let engine = Engine { axes, time, symbol: Symbol::Cartesian(spatial) };
    finish(cx, p, conditions, &engine, &items)
}

/// The initial-data items of an evolution problem (none for a steady
/// one); every condition on the time variable must be an initial one.
pub(super) fn time_items(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
    time_index: Option<usize>,
    order: u32,
    n: usize,
) -> Option<Vec<Item>> {
    let mut items = Vec::new();
    let Some(t) = time_index else {
        return Some(items);
    };
    let zero_index = vec![0; n];
    for c in conditions.0.iter().filter(|c| c.on == t) {
        if !is_zero_number(cx.graph, c.point) {
            return None;
        }
    }
    let disp = conditions.find(cx.graph, t, &zero_index, None);
    let vel = conditions.find(cx.graph, t, &p.unit(t, 1), None);
    if conditions.0.iter().filter(|c| c.on == t).count() != usize::from(disp.is_some()) + usize::from(vel.is_some()) {
        return None;
    }
    match (order, disp, vel) {
        | (1, Some(d), None) => items.push(Item { role: Role::Displacement, h: d.value, special: None }),
        | (2, d, v) if d.is_some() || v.is_some() => {
            if let Some(d) = d {
                items.push(Item { role: Role::Displacement, h: d.value, special: None });
            }
            if let Some(v) = v {
                items.push(Item { role: Role::Velocity, h: v.value, special: None });
            }
        },
        | _ => return None,
    }
    Some(items)
}

/// The forcing items of the boundary data at both ends of an interval
/// axis (`cjj` the coefficient of the second derivative there).
pub(super) fn boundary_items(
    cx: &mut Cx<'_>,
    iv: &Interval,
    level: usize,
    cjj: NodeId,
) -> Option<Vec<Item>> {
    let mut items = Vec::new();
    for (side, end) in [(-1_i64, iv.left), (1, iv.right)] {
        if iv.periodic || cx.is_zero(end.value) {
            continue;
        }
        let divisor = if end.kind == Kind::Dirichlet { end.p } else { end.q };
        let sign = cx.graph.int(side);
        let scaled = mul(cx.graph, &[sign, cjj, end.value]);
        let scaled = div(cx, scaled, divisor);
        let scaled = cx.simplify(scaled);
        items.push(Item {
            role: Role::Source,
            h: scaled,
            special: Some(Special { axis: level, kind: end.kind, point: end.point, weight: cx.graph.int(1) }),
        });
    }
    Some(items)
}

/// Runs the engine and checks closed-form results against the equation and
/// the conditions.
pub(super) fn finish(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
    engine: &Engine,
    items: &[Item],
) -> Option<NodeId> {
    let solution = run(cx, p, engine, items)?;
    if cx.is_zero(solution) {
        return None;
    }
    if is_closed(cx, solution) && (!verified(cx, p, solution) || !satisfies(cx, p, conditions, solution)) {
        return None;
    }
    Some(solution)
}
