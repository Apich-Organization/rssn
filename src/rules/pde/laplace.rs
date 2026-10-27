//! Laplace transform in time for signalling problems on a half-line.
//!
//! `a₂ u_tt + a₁ u_t + a₀ u + b u_xx = 0` on `x > 0` with data
//! `u(0, t) = g(t)` (Dirichlet) or `u_x(0, t) = g(t)` (Neumann), zero initial
//! displacement and velocity, and a wave leaving the boundary. With
//! `c² = -b/a₂`, `γ = a₁/(2a₂)`, `m² = a₀/a₂ - γ²` and `u = e^{-γt} w` the
//! transform in `t` turns the equation into `c² W'' = (s² + m²) W`, so that
//!
//! * Dirichlet: `W = G(s) e^{-x√(s²+m²)/c}`,
//! * Neumann:   `W = -c G(s) e^{-x√(s²+m²)/c}/√(s²+m²)`,
//!
//! and the inversion table (`L⁻¹[e^{-a√(s²+m²)}] = δ(t-a) - a m J₁(m√(t²-a²))/√(t²-a²)`,
//! `L⁻¹[e^{-a√(s²+m²)}/√(s²+m²)] = J₀(m√(t²-a²))`, `I` for `m² < 0`) with
//! `a = x/c` and the convolution theorem give the solution as a
//! retarded integral over the boundary datum.

use super::Conditions;
use super::Problem;
use super::dummy;
use super::spectral::evolution_time;
use super::util::call;
use super::util::defint;
use super::util::div;
use super::util::is_zero_number;
use super::util::sample;
use super::util::sqrt;
use crate::graph::Cx;
use crate::graph::Facts;
use crate::graph::NodeId;
use crate::rules::complex::build::add;
use crate::rules::complex::build::mul;
use crate::rules::complex::build::neg;
use crate::rules::complex::build::powi;
use crate::rules::complex::build::sub;

/// The retarded solution of the half-line signalling problem.
#[allow(clippy::too_many_lines)]
pub(super) fn signalling(
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
    let zero_index = vec![0; 2];
    let (a2, a1, a0, b) = (
        p.coefficient(cx.graph, &p.unit(time, 2)),
        p.coefficient(cx.graph, &p.unit(time, 1)),
        p.coefficient(cx.graph, &zero_index),
        p.coefficient(cx.graph, &p.unit(space, 2)),
    );
    if cx.is_zero(a2) || cx.is_zero(b) || ![a2, a1, a0, b].iter().all(|&c| p.constant(cx.graph, c)) {
        return None;
    }
    if !p.only(cx.graph, &[p.unit(time, 2), p.unit(time, 1), zero_index.clone(), p.unit(space, 2)]) {
        return None;
    }
    // The boundary datum and the quiet initial state.
    let mut boundary = None;
    let mut counted = 0;
    for c in &conditions.0 {
        if c.robin.is_some() {
            return None;
        }
        if c.on == space {
            if !is_zero_number(cx.graph, c.point) || boundary.is_some() {
                return None;
            }
            let neumann = if c.derivative == zero_index {
                false
            } else if c.derivative == p.unit(space, 1) {
                true
            } else {
                return None;
            };
            boundary = Some((neumann, c.value));
        } else if !(c.on == time && (c.derivative == zero_index || c.derivative == p.unit(time, 1)) && is_zero_number(cx.graph, c.point) && cx.is_zero(c.value)) {
            return None;
        }
        counted += 1;
    }
    if counted != conditions.0.len() {
        return None;
    }
    let (neumann, g) = boundary?;
    if cx.graph.depends_on(cx.graph.find(g), cx.graph.symbol_of(x)?) {
        return None;
    }
    let ratio = |cx: &mut Cx<'_>, num: NodeId, den: NodeId| {
        let q = div(cx, num, den);
        cx.simplify(q)
    };
    let c2 = {
        let q = ratio(cx, b, a2);
        let q = neg(cx.graph, q);
        cx.simplify(q)
    };
    if cx.graph.facts(c2).has(Facts::NEGATIVE) || sample(cx.graph, c2, 0).is_none_or(|v| v <= 0.0) {
        return None;
    }
    let two = cx.graph.int(2);
    let gamma = {
        let q = ratio(cx, a1, a2);
        let q = div(cx, q, two);
        cx.simplify(q)
    };
    let m2 = {
        let mass = ratio(cx, a0, a2);
        let g2 = powi(cx.graph, gamma, 2);
        let m = sub(cx.graph, mass, g2);
        cx.simplify(m)
    };
    let kind = if cx.is_zero(m2) {
        0
    } else if cx.graph.facts(m2).has(Facts::POSITIVE) || sample(cx.graph, m2, 0).is_some_and(|v| v > 0.0) {
        1
    } else if sample(cx.graph, m2, 0).is_some_and(|v| v < 0.0) {
        -1
    } else {
        return None;
    };
    let speed = sqrt(cx, c2)?;
    let delay = {
        let q = div(cx, x, speed);
        cx.simplify(q)
    };
    let (tau, _) = dummy(cx, p, "tau");
    let retarded = sub(cx.graph, t, delay);
    // The datum in the gauge variable: `e^{γ τ} g(τ)`.
    let gauge = |cx: &mut Cx<'_>, at: NodeId| -> Option<NodeId> {
        let value = cx.graph.substitute(g, t, at);
        if cx.is_zero(gamma) {
            return Some(value);
        }
        let exponent = mul(cx.graph, &[gamma, at]);
        let e = call(cx, "exp", &[exponent])?;
        Some(mul(cx.graph, &[e, value]))
    };
    let g_tau = gauge(cx, tau)?;
    // The memory kernel `K(t - τ)`.
    let lag = sub(cx.graph, t, tau);
    let kernel = if kind == 0 {
        None
    } else {
        let m = if kind > 0 { sqrt(cx, m2)? } else {
            let n = neg(cx.graph, m2);
            sqrt(cx, n)?
        };
        let gap = {
            let l2 = powi(cx.graph, lag, 2);
            let d2 = powi(cx.graph, delay, 2);
            sub(cx.graph, l2, d2)
        };
        let root = sqrt(cx, gap)?;
        let argument = mul(cx.graph, &[m, root]);
        let function = if kind > 0 { "besselj" } else { "besseli" };
        let order_node = cx.graph.int(i64::from(!neumann));
        let bessel = call(cx, function, &[order_node, argument])?;
        let k = if neumann {
            bessel
        } else {
            let scaled = mul(cx.graph, &[m, delay, bessel]);
            let k = div(cx, scaled, root);
            if kind > 0 { neg(cx.graph, k) } else { k }
        };
        Some(k)
    };
    let mut terms = Vec::new();
    if !neumann {
        // `g(t - a) e^{-γ a}`.
        let shifted = cx.graph.substitute(g, t, retarded);
        let direct = if cx.is_zero(gamma) {
            shifted
        } else {
            let exponent = {
                let ga = mul(cx.graph, &[gamma, delay]);
                neg(cx.graph, ga)
            };
            let e = call(cx, "exp", &[exponent])?;
            mul(cx.graph, &[e, shifted])
        };
        terms.push(direct);
    }
    if neumann || kernel.is_some() {
        let zero = cx.graph.int(0);
        let body = match kernel {
            | Some(k) => mul(cx.graph, &[g_tau, k]),
            | None => g_tau,
        };
        let integral = defint(cx, body, tau, zero, retarded)?;
        let decay = if cx.is_zero(gamma) {
            integral
        } else {
            let exponent = {
                let gt = mul(cx.graph, &[gamma, t]);
                neg(cx.graph, gt)
            };
            let e = call(cx, "exp", &[exponent])?;
            mul(cx.graph, &[e, integral])
        };
        if neumann {
            let minus = neg(cx.graph, speed);
            terms.push(mul(cx.graph, &[minus, decay]));
        } else {
            terms.push(decay);
        }
    }
    let total = add(cx.graph, &terms);
    let front = call(cx, "heaviside", &[retarded])?;
    let u = mul(cx.graph, &[front, total]);
    Some(cx.simplify(u))
}
