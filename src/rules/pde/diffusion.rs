//! Parabolic equations on whole space and on half-spaces, in any number of
//! dimensions.
//!
//! `a u_t = Σ_j (b_j u_{x_j x_j} + c_j u_{x_j}) + d u + f(x, t)` with
//! constant coefficients (`D_j = b_j / a`) on `R^d`, or on the half-space
//! `x_j > 0` for any subset of the variables with a homogeneous Dirichlet
//! or Neumann condition there:
//!
//! * the gauge `u = exp(Σ α_j x_j + β t) v`, `α_j = -c_j/(2 D_j a)`,
//!   `β = d/a - Σ c_j²/(4 D_j a)` removes drift and reaction;
//! * the heat kernel `Π K_j` gives the initial-value solution, with the
//!   method of images (`K(x - s) ∓ K(x + s)`) on every half-space axis;
//! * Duhamel's principle gives the source term;
//! * one half-space axis may carry boundary data `g(y, t)`: the
//!   Dirichlet kernel `x_n/(t - τ) Π K` (Neumann: `-2 D_n Π K`).
//!
//! Complex `D` (the free Schrödinger equation) is covered by the same
//! formulas. A constant-data half-line problem is solved in closed form
//! by the similarity variable `x / (2 √(D t))`.

use super::Conditions;
use super::Problem;
use super::classical::kernel;
use super::dummy;
use super::spectral::evolution_time;
use super::util::call;
use super::util::defint;
use super::util::div;
use super::util::infinity;
use super::util::is_zero_number;
use super::util::sqrt;
use super::verified;
use crate::graph::Cx;
use crate::graph::Facts;
use crate::graph::NodeId;
use crate::rules::complex::build::add;
use crate::rules::complex::build::mul;
use crate::rules::complex::build::neg;
use crate::rules::complex::build::powi;
use crate::rules::complex::build::sub;

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Domain {
    Whole,
    /// `x > 0` with a Dirichlet (`-1`) or Neumann (`+1`) condition.
    Half(i64),
}

struct Coefficients {
    time: usize,
    space: Vec<usize>,
    diffusivity: Vec<NodeId>,
    drift: Vec<NodeId>,
    reaction: NodeId,
    /// `-source / a`.
    forcing: NodeId,
}

fn coefficients(
    cx: &mut Cx<'_>,
    p: &Problem,
) -> Option<Coefficients> {
    if p.nonlinear {
        return None;
    }
    let time = evolution_time(cx, p)?;
    let n = p.dimension();
    let space: Vec<usize> = (0..n).filter(|&j| j != time).collect();
    let a = p.coefficient(cx.graph, &p.unit(time, 1));
    let a2 = p.coefficient(cx.graph, &p.unit(time, 2));
    if cx.is_zero(a) || !p.constant(cx.graph, a) || !is_zero_number(cx.graph, a2) {
        return None;
    }
    let ratio = |cx: &mut Cx<'_>, index: &[u32]| {
        let c = p.coefficient(cx.graph, index);
        let q = div(cx, c, a);
        let q = neg(cx.graph, q);
        cx.simplify(q)
    };
    let mut diffusivity = Vec::new();
    let mut drift = Vec::new();
    for &j in &space {
        let d = ratio(cx, &p.unit(j, 2));
        if cx.is_zero(d) || !p.constant(cx.graph, d) {
            return None;
        }
        // Real positive (or not decidable) diffusivity.
        let facts = cx.graph.facts(d);
        let value = cx.graph.eval(d, &crate::graph::Env::numeric(0.0));
        if facts.has(Facts::NEGATIVE) || value.is_some_and(|v| v <= 0.0) {
            return None;
        }
        diffusivity.push(d);
        let c = ratio(cx, &p.unit(j, 1));
        if !p.constant(cx.graph, c) {
            return None;
        }
        drift.push(c);
    }
    let reaction = ratio(cx, &vec![0; n]);
    if !p.constant(cx.graph, reaction) {
        return None;
    }
    // No other derivatives.
    for (index, c) in &p.linear {
        let known = *index == p.unit(time, 1)
            || *index == vec![0; n]
            || space.iter().any(|&j| *index == p.unit(j, 1) || *index == p.unit(j, 2));
        if !known && !is_zero_number(cx.graph, *c) {
            return None;
        }
    }
    let forcing = {
        let q = div(cx, p.source, a);
        let q = neg(cx.graph, q);
        cx.simplify(q)
    };
    Some(Coefficients { time, space, diffusivity, drift, reaction, forcing })
}

/// The Gauss–Weierstrass factor of one axis: `K(x - s)`, plus the image
/// `± K(x + s)` on a half-line.
fn axis_kernel(
    cx: &mut Cx<'_>,
    d: NodeId,
    x: NodeId,
    s: NodeId,
    tau: NodeId,
    domain: Domain,
) -> Option<NodeId> {
    let direct = sub(cx.graph, x, s);
    let k = kernel(cx, d, direct, tau)?;
    match domain {
        | Domain::Whole => Some(k),
        | Domain::Half(sign) => {
            let reflected = add(cx.graph, &[x, s]);
            let image = kernel(cx, d, reflected, tau)?;
            let sign = cx.graph.int(sign);
            let image = mul(cx.graph, &[sign, image]);
            Some(add(cx.graph, &[k, image]))
        },
    }
}

/// `∫ body ds_j` over each axis (innermost last), simplifying as it goes.
fn integrate_axes(
    cx: &mut Cx<'_>,
    mut body: NodeId,
    dummies: &[NodeId],
    domains: &[Domain],
    skip: Option<usize>,
) -> Option<NodeId> {
    let oo = infinity(cx)?;
    let minus_oo = neg(cx.graph, oo);
    for (j, (&s, &domain)) in dummies.iter().zip(domains).enumerate().rev() {
        if Some(j) == skip {
            continue;
        }
        let lower = if domain == Domain::Whole { minus_oo } else { cx.graph.int(0) };
        body = defint(cx, body, s, lower, oo)?;
        body = cx.simplify(body);
    }
    Some(body)
}

/// Initial-value and boundary-value problems on whole space and
/// half-spaces.
pub(super) fn parabolic(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
) -> Option<NodeId> {
    let c = coefficients(cx, p)?;
    let n = p.dimension();
    let zero_index = vec![0; n];
    let initial = conditions.find(cx.graph, c.time, &zero_index, None)?;
    if !is_zero_number(cx.graph, initial.point) {
        return None;
    }
    // Domains and boundary data.
    let mut domains = Vec::new();
    let mut data: Vec<(usize, NodeId)> = Vec::new();
    let mut counted = 1;
    for (axis, &j) in c.space.iter().enumerate() {
        let on: Vec<_> = conditions.0.iter().filter(|q| q.on == j).collect();
        counted += on.len();
        match on.as_slice() {
            | [] => domains.push(Domain::Whole),
            | [b] if b.robin.is_none() && is_zero_number(cx.graph, b.point) => {
                let dirichlet = b.derivative == zero_index;
                if !dirichlet && b.derivative != p.unit(j, 1) {
                    return None;
                }
                domains.push(Domain::Half(if dirichlet { -1 } else { 1 }));
                if !cx.is_zero(b.value) {
                    data.push((axis, b.value));
                }
            },
            | _ => return None,
        }
    }
    if conditions.0.len() != counted || data.len() > 1 {
        return None;
    }
    let t = p.vars[c.time];
    // Gauge.
    let two = cx.graph.int(2);
    let four = cx.graph.int(4);
    let mut alphas = Vec::new();
    let mut beta_terms = vec![c.reaction];
    for (axis, &d) in c.diffusivity.iter().enumerate() {
        let two_d = mul(cx.graph, &[two, d]);
        let q = div(cx, c.drift[axis], two_d);
        let alpha = neg(cx.graph, q);
        let alpha = cx.simplify(alpha);
        // Images need the gauge to be compatible with the condition:
        // Dirichlet always, Neumann only without drift.
        if domains[axis] == Domain::Half(1) && !cx.is_zero(alpha) {
            return None;
        }
        alphas.push(alpha);
        let c2 = powi(cx.graph, c.drift[axis], 2);
        let four_d = mul(cx.graph, &[four, d]);
        let q = div(cx, c2, four_d);
        beta_terms.push(neg(cx.graph, q));
    }
    let beta = add(cx.graph, &beta_terms);
    let beta = cx.simplify(beta);
    // Dummies of integration.
    let stems = ["s", "r", "q", "w", "z", "y"];
    let mut dummies = Vec::new();
    for axis in 0..c.space.len() {
        let stem = stems.get(axis).copied().unwrap_or("s");
        dummies.push(dummy(cx, p, stem).0);
    }
    let phase_at = |cx: &mut Cx<'_>, coords: &[NodeId], skip: Option<usize>| -> NodeId {
        let terms: Vec<NodeId> =
            coords.iter().zip(&alphas).enumerate().filter(|(j, _)| Some(*j) != skip).map(|(_, (&x, &a))| mul(cx.graph, &[a, x])).collect();
        add(cx.graph, &terms)
    };
    let xs: Vec<NodeId> = c.space.iter().map(|&j| p.vars[j]).collect();
    let substitute_space = |cx: &mut Cx<'_>, mut node: NodeId, to: &[NodeId], skip: Option<usize>| -> NodeId {
        for (j, (&x, &s)) in xs.iter().zip(to).enumerate() {
            if Some(j) != skip {
                node = cx.graph.substitute(node, x, s);
            }
        }
        node
    };
    let product = |cx: &mut Cx<'_>, d: &[NodeId], tau: NodeId, skip: Option<usize>| -> Option<NodeId> {
        let mut factors = Vec::new();
        for (axis, &x) in xs.iter().enumerate() {
            if Some(axis) == skip {
                continue;
            }
            factors.push(axis_kernel(cx, c.diffusivity[axis], x, d[axis], tau, domains[axis])?);
        }
        Some(mul(cx.graph, &factors))
    };
    let minus_one = cx.graph.int(-1);
    // Initial data.
    let f_at = substitute_space(cx, initial.value, &dummies, None);
    let phase_s = phase_at(cx, &dummies, None);
    let gauge_s = {
        let e = mul(cx.graph, &[minus_one, phase_s]);
        call(cx, "exp", &[e])?
    };
    let weight = product(cx, &dummies, t, None)?;
    let integrand = mul(cx.graph, &[weight, gauge_s, f_at]);
    let mut parts = Vec::new();
    if !cx.is_zero(initial.value) {
        parts.push(integrate_axes(cx, integrand, &dummies, &domains, None)?);
    }
    let zero = cx.graph.int(0);
    // Duhamel's integral for the source.
    if !cx.is_zero(c.forcing) {
        let (tau, _) = dummy(cx, p, "tau");
        let source = substitute_space(cx, c.forcing, &dummies, None);
        let source = cx.graph.substitute(source, t, tau);
        let elapsed = sub(cx.graph, t, tau);
        let weight = product(cx, &dummies, elapsed, None)?;
        let phase_s = phase_at(cx, &dummies, None);
        let b_tau = mul(cx.graph, &[beta, tau]);
        let gauge = {
            let e = add(cx.graph, &[phase_s, b_tau]);
            let e = mul(cx.graph, &[minus_one, e]);
            call(cx, "exp", &[e])?
        };
        let body = mul(cx.graph, &[weight, gauge, source]);
        let inner = integrate_axes(cx, body, &dummies, &domains, None)?;
        parts.push(defint(cx, inner, tau, zero, t)?);
    }
    // Boundary data on one half-space axis.
    if let Some(&(axis, g)) = data.first() {
        let sign = match domains[axis] {
            | Domain::Half(s) => s,
            | Domain::Whole => return None,
        };
        let (tau, _) = dummy(cx, p, "tau");
        let g_at = substitute_space(cx, g, &dummies, Some(axis));
        let g_at = cx.graph.substitute(g_at, t, tau);
        let elapsed = sub(cx.graph, t, tau);
        let weight = product(cx, &dummies, elapsed, Some(axis))?;
        let x_n = xs[axis];
        let k_n = kernel(cx, c.diffusivity[axis], x_n, elapsed)?;
        let phase_y = phase_at(cx, &dummies, Some(axis));
        let b_tau = mul(cx.graph, &[beta, tau]);
        let gauge = {
            let e = add(cx.graph, &[phase_y, b_tau]);
            let e = mul(cx.graph, &[minus_one, e]);
            call(cx, "exp", &[e])?
        };
        let normal = if sign < 0 {
            // Dirichlet: x_n / (t - τ) K_n.
            let factor = div(cx, x_n, elapsed);
            mul(cx.graph, &[factor, k_n])
        } else {
            // Neumann: -2 D_n K_n.
            let m = mul(cx.graph, &[minus_one, two, c.diffusivity[axis], k_n]);
            cx.simplify(m)
        };
        let body = mul(cx.graph, &[normal, weight, gauge, g_at]);
        let inner = integrate_axes(cx, body, &dummies, &domains, Some(axis))?;
        parts.push(defint(cx, inner, tau, zero, t)?);
    }
    let total = add(cx.graph, &parts);
    let phase_x = phase_at(cx, &xs, None);
    let b_t = mul(cx.graph, &[beta, t]);
    let exponent = add(cx.graph, &[phase_x, b_t]);
    let gauge = call(cx, "exp", &[exponent])?;
    let u = mul(cx.graph, &[gauge, total]);
    Some(cx.simplify(u))
}

/// Constant data on a half-line: `u(x, 0) = u₀`, `u(0, t) = g` (or
/// `u_x(0, t) = g`) for `u_t = D u_xx`: the similarity solution in
/// `ξ = x / (2 √(D t))`.
pub(super) fn similarity(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
) -> Option<NodeId> {
    let c = coefficients(cx, p)?;
    if c.space.len() != 1 || !cx.is_zero(c.forcing) || !cx.is_zero(c.reaction) || !cx.is_zero(c.drift[0]) {
        return None;
    }
    let n = p.dimension();
    let zero_index = vec![0; n];
    let (time, space) = (c.time, c.space[0]);
    let (x, t) = (p.vars[space], p.vars[time]);
    let initial = conditions.find(cx.graph, time, &zero_index, None)?;
    if conditions.0.len() != 2 || !is_zero_number(cx.graph, initial.point) {
        return None;
    }
    let free = |cx: &Cx<'_>, v: NodeId| p.vars.iter().all(|&q| cx.graph.symbol_of(q).is_none_or(|s| !cx.graph.depends_on(cx.graph.find(v), s)));
    let boundary = conditions.0.iter().find(|q| q.on == space)?;
    if !free(cx, initial.value) || !free(cx, boundary.value) || boundary.robin.is_some() || !is_zero_number(cx.graph, boundary.point) {
        return None;
    }
    let u0 = initial.value;
    let d = c.diffusivity[0];
    let dt = mul(cx.graph, &[d, t]);
    let root = sqrt(cx, dt)?;
    let two = cx.graph.int(2);
    let denominator = mul(cx.graph, &[two, root]);
    let xi = div(cx, x, denominator);
    let erfc = call(cx, "erfc", &[xi])?;
    let solution = if boundary.derivative == zero_index {
        // u₀ + (g - u₀) erfc(ξ)
        let gap = sub(cx.graph, boundary.value, u0);
        let m = mul(cx.graph, &[gap, erfc]);
        add(cx.graph, &[u0, m])
    } else if boundary.derivative == p.unit(space, 1) {
        // u₀ - 2 g √(D t / π) exp(-ξ²) + g x erfc(ξ)
        let pi = super::util::pi(cx)?;
        let ratio = div(cx, dt, pi);
        let r = sqrt(cx, ratio)?;
        let xi2 = powi(cx.graph, xi, 2);
        let minus_xi2 = neg(cx.graph, xi2);
        let gaussian = call(cx, "exp", &[minus_xi2])?;
        let minus_two = cx.graph.int(-2);
        let first = mul(cx.graph, &[minus_two, boundary.value, r, gaussian]);
        let second = mul(cx.graph, &[boundary.value, x, erfc]);
        add(cx.graph, &[u0, first, second])
    } else {
        return None;
    };
    let solution = cx.simplify(solution);
    verified(cx, p, solution).then_some(solution)
}
