//! Waves in free space and half-spaces, in one, two and three dimensions.
//!
//! `a u_tt + b u_t + c_s Σ_j u_{x_j x_j} + e u + f = 0` with initial data
//! `u(x, 0)`, `u_t(x, 0)`:
//!
//! * one dimension: d'Alembert's formula; with damping `b` (telegraph) and
//!   mass `e` (Klein–Gordon) Riemann's formula with `J_0` / `I_0` and the
//!   `J_1`-term for the displacement, after `u = e^{-γ t} w`;
//! * two dimensions: Poisson's formula, by Hadamard's method of descent from
//!   three dimensions: `∂_t (t M₂[f](ct)) + t M₂[g](ct)` with the weighted
//!   mean `M₂[h](r) = (1/2π) ∬ h(x + rρω) ρ/√(1-ρ²) dρ dφ`;
//! * three dimensions: Kirchhoff's formula with the spherical mean;
//! * sources by Duhamel's principle (the formula for the velocity data
//!   `F(·, τ)` at time `t - τ`);
//! * half-spaces `x_j > 0` with a homogeneous Dirichlet or Neumann
//!   condition: odd or even extension of data and source in `x_j`.

use super::Conditions;
use super::Problem;
use super::dummy;
use super::spectral::evolution_time;
use super::util::call;
use super::util::defint;
use super::util::div;
use super::util::fraction;
use super::util::is_zero_number;
use super::util::pi;
use super::util::sample;
use super::util::sqrt;
use crate::graph::Cx;
use crate::graph::Facts;
use crate::graph::NodeId;
use crate::rules::calculus::derivative;
use crate::rules::complex::build::add;
use crate::rules::complex::build::mul;
use crate::rules::complex::build::neg;
use crate::rules::complex::build::powi;
use crate::rules::complex::build::sub;
use crate::rules::poly::repr::Gens;
use crate::rules::poly::repr::Limits;
use crate::rules::poly::repr::from_term;
use crate::rules::poly::repr::to_term;

/// `+1` for a Neumann half-space axis, `-1` for Dirichlet, `0` for the
/// whole line.
type Parity = i64;

/// The extension of `h` across the half-space boundaries: even in an axis
/// with a Neumann condition, odd with a Dirichlet one.
fn extend(
    cx: &mut Cx<'_>,
    h: NodeId,
    xs: &[NodeId],
    parity: &[Parity],
) -> Option<NodeId> {
    let mut out = h;
    for (&x, &q) in xs.iter().zip(parity) {
        if q == 0 {
            continue;
        }
        let magnitude = call(cx, "abs", &[x])?;
        let folded = cx.graph.substitute(out, x, magnitude);
        out = if q < 0 {
            let s = call(cx, "sign", &[x])?;
            mul(cx.graph, &[s, folded])
        } else {
            folded
        };
    }
    Some(out)
}

/// The spherical mean (`d = 3`) or the weighted circular mean (`d = 2`) of
/// `h` over the ball of radius `radius` about `x`.
fn mean(
    cx: &mut Cx<'_>,
    p: &Problem,
    xs: &[NodeId],
    h: NodeId,
    radius: NodeId,
) -> Option<NodeId> {
    let d = xs.len();
    let pi = pi(cx)?;
    let zero = cx.graph.int(0);
    let one = cx.graph.int(1);
    let two = cx.graph.int(2);
    let two_pi = mul(cx.graph, &[two, pi]);
    let (phi, _) = dummy(cx, p, "phi");
    let (sin, cos) = (cx.graph.ops().lookup("sin")?, cx.graph.ops().lookup("cos")?);
    let sp = cx.graph.node(sin, &[phi]);
    let cp = cx.graph.node(cos, &[phi]);
    // Directions, scale factors and measure.
    let (directions, measure, outer_variable, outer_hi, normalisation): (Vec<NodeId>, NodeId, NodeId, NodeId, NodeId) = if d == 3 {
        let (theta, _) = dummy(cx, p, "theta");
        let st = cx.graph.node(sin, &[theta]);
        let ct = cx.graph.node(cos, &[theta]);
        let dirs = vec![mul(cx.graph, &[st, cp]), mul(cx.graph, &[st, sp]), ct];
        let four = cx.graph.int(4);
        let four_pi = mul(cx.graph, &[four, pi]);
        (dirs, st, theta, pi, powi(cx.graph, four_pi, -1))
    } else {
        let (rho, _) = dummy(cx, p, "rho");
        let dirs = vec![mul(cx.graph, &[rho, cp]), mul(cx.graph, &[rho, sp])];
        let root = {
            let r2 = powi(cx.graph, rho, 2);
            let gap = sub(cx.graph, one, r2);
            sqrt(cx, gap)?
        };
        let weight = div(cx, rho, root);
        (dirs, weight, rho, one, powi(cx.graph, two_pi, -1))
    };
    let mut shifted = h;
    let fresh: Vec<NodeId> = (0..d)
        .map(|_| {
            let s = cx.graph.interner_mut().fresh_symbol("w");
            cx.graph.symbol_node(s)
        })
        .collect();
    for (k, &x) in xs.iter().enumerate() {
        shifted = cx.graph.substitute(shifted, x, fresh[k]);
    }
    for (k, &x) in xs.iter().enumerate() {
        let step = mul(cx.graph, &[radius, directions[k]]);
        let moved = add(cx.graph, &[x, step]);
        shifted = cx.graph.substitute(shifted, fresh[k], moved);
    }
    // Expanded, the integrand is a sum of trigonometric monomials.
    let mut gens = Gens::default();
    if let Some(poly) = from_term(cx.graph, &mut gens, shifted, Limits { terms: 4096, exponent: 16 }) {
        shifted = to_term(cx.graph, &gens, &poly);
    }
    let body = mul(cx.graph, &[shifted, measure]);
    let inner = defint(cx, body, phi, zero, two_pi)?;
    let inner = cx.simplify(inner);
    let outer = defint(cx, inner, outer_variable, zero, outer_hi)?;
    let value = mul(cx.graph, &[normalisation, outer]);
    Some(cx.simplify(value))
}

/// Initial-value problems for the wave equation (with damping, mass and
/// sources) on whole space and half-spaces.
#[allow(clippy::too_many_lines)]
pub(super) fn wave(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
) -> Option<NodeId> {
    if p.nonlinear {
        return None;
    }
    let time = evolution_time(cx, p)?;
    let n = p.dimension();
    let space: Vec<usize> = (0..n).filter(|&j| j != time).collect();
    let d = space.len();
    if !(1..=3).contains(&d) {
        return None;
    }
    let a = p.coefficient(cx.graph, &p.unit(time, 2));
    if cx.is_zero(a) || !p.constant(cx.graph, a) {
        return None;
    }
    let ratio = |cx: &mut Cx<'_>, index: &[u32]| {
        let c = p.coefficient(cx.graph, index);
        let q = div(cx, c, a);
        cx.simplify(q)
    };
    // Isotropic: c² = -C/A.
    let mut c2 = None;
    for &j in &space {
        let q = ratio(cx, &p.unit(j, 2));
        let q = neg(cx.graph, q);
        let q = cx.simplify(q);
        if cx.is_zero(q) || !p.constant(cx.graph, q) {
            return None;
        }
        match c2 {
            | None => c2 = Some(q),
            | Some(c) => {
                let gap = sub(cx.graph, c, q);
                if !cx.is_zero(gap) {
                    return None;
                }
            },
        }
    }
    let c2 = c2?;
    let facts = cx.graph.facts(c2);
    if facts.has(Facts::NEGATIVE) || sample(cx.graph, c2, 0).is_some_and(|v| v <= 0.0) {
        return None;
    }
    let damping = ratio(cx, &p.unit(time, 1));
    let mass = ratio(cx, &vec![0; n]);
    for (index, c) in &p.linear {
        let known = *index == p.unit(time, 2)
            || *index == p.unit(time, 1)
            || *index == vec![0; n]
            || space.iter().any(|&j| *index == p.unit(j, 2));
        if !known && !is_zero_number(cx.graph, *c) {
            return None;
        }
    }
    if !p.constant(cx.graph, damping) || !p.constant(cx.graph, mass) {
        return None;
    }
    let half = fraction(cx, 1, 2)?;
    let two = cx.graph.int(2);
    let gamma = {
        let g = div(cx, damping, two);
        cx.simplify(g)
    };
    let m2 = {
        let g2 = powi(cx.graph, gamma, 2);
        let m = sub(cx.graph, mass, g2);
        cx.simplify(m)
    };
    let massless = cx.is_zero(m2) && cx.is_zero(gamma);
    if d > 1 && !massless {
        return None;
    }
    // Domains and conditions.
    let zero_index = vec![0; n];
    let displacement = conditions.find(cx.graph, time, &zero_index, None);
    let velocity = conditions.find(cx.graph, time, &p.unit(time, 1), None);
    if displacement.is_none() && velocity.is_none() {
        return None;
    }
    let mut parity = Vec::new();
    let mut counted = usize::from(displacement.is_some()) + usize::from(velocity.is_some());
    for &j in &space {
        let on: Vec<_> = conditions.0.iter().filter(|c| c.on == j).collect();
        counted += on.len();
        match on.as_slice() {
            | [] => parity.push(0),
            | [b] if b.robin.is_none() && is_zero_number(cx.graph, b.point) && cx.is_zero(b.value) => {
                if b.derivative == zero_index {
                    parity.push(-1);
                } else if b.derivative == p.unit(j, 1) {
                    parity.push(1);
                } else {
                    return None;
                }
            },
            | _ => return None,
        }
    }
    if conditions.0.len() != counted {
        return None;
    }
    for q in displacement.iter().chain(velocity.iter()) {
        if !is_zero_number(cx.graph, q.point) {
            return None;
        }
    }
    let t = p.vars[time];
    let xs: Vec<NodeId> = space.iter().map(|&j| p.vars[j]).collect();
    let zero = cx.graph.int(0);
    let c = sqrt(cx, c2)?;
    let ct = mul(cx.graph, &[c, t]);
    // Data of w = e^{γ t} u: w(0) = f, w_t(0) = g + γ f; forcing
    // e^{γ t} (-source / a).
    let f = displacement.map_or(zero, |q| q.value);
    let g = {
        let g0 = velocity.map_or(zero, |q| q.value);
        let gf = mul(cx.graph, &[gamma, f]);
        add(cx.graph, &[g0, gf])
    };
    let g = cx.simplify(g);
    let forcing = {
        let q = div(cx, p.source, a);
        let q = neg(cx.graph, q);
        let e = if cx.is_zero(gamma) {
            q
        } else {
            let gt = mul(cx.graph, &[gamma, t]);
            let e = call(cx, "exp", &[gt])?;
            mul(cx.graph, &[e, q])
        };
        cx.simplify(e)
    };
    let (f_e, g_e, forcing_e) = (extend(cx, f, &xs, &parity)?, extend(cx, g, &xs, &parity)?, extend(cx, forcing, &xs, &parity)?);
    let mut terms = Vec::new();
    if d == 1 {
        let x = xs[0];
        let (s, _) = dummy(cx, p, "s");
        let (tau, _) = dummy(cx, p, "tau");
        let right = add(cx.graph, &[x, ct]);
        let left = sub(cx.graph, x, ct);
        // Bessel kernels for the mass term.
        let kind = if cx.is_zero(m2) {
            0
        } else if cx.graph.facts(m2).has(Facts::POSITIVE) || sample(cx.graph, m2, 0).is_some_and(|v| v > 0.0) {
            1
        } else if cx.graph.facts(m2).has(Facts::NEGATIVE) || sample(cx.graph, m2, 0).is_some_and(|v| v < 0.0) {
            -1
        } else {
            return None;
        };
        let m_abs = if kind < 0 { neg(cx.graph, m2) } else { m2 };
        let m = sqrt(cx, m_abs)?;
        let mu = div(cx, m, c);
        // Radius function R(τ', s) = (c² τ'² - (x - s)²)^(1/2).
        let radius = |cx: &mut Cx<'_>, elapsed: NodeId| -> Option<NodeId> {
            let ce = mul(cx.graph, &[c, elapsed]);
            let a2 = powi(cx.graph, ce, 2);
            let diff = sub(cx.graph, x, s);
            let b2 = powi(cx.graph, diff, 2);
            let gap = sub(cx.graph, a2, b2);
            sqrt(cx, gap)
        };
        let bessel = |cx: &mut Cx<'_>, order: i64, r: NodeId| -> Option<NodeId> {
            let mur = mul(cx.graph, &[mu, r]);
            let o = cx.graph.int(order);
            call(cx, if kind > 0 { "besselj" } else { "besseli" }, &[o, mur])
        };
        let f_right = cx.graph.substitute(f_e, x, right);
        let f_left = cx.graph.substitute(f_e, x, left);
        let sum = add(cx.graph, &[f_right, f_left]);
        terms.push(mul(cx.graph, &[half, sum]));
        let two_c = mul(cx.graph, &[two, c]);
        let inverse_two_c = powi(cx.graph, two_c, -1);
        if !cx.is_zero(g) {
            let g_s = cx.graph.substitute(g_e, x, s);
            let r = radius(cx, t)?;
            let weight = if kind == 0 { cx.graph.int(1) } else { bessel(cx, 0, r)? };
            let body = mul(cx.graph, &[g_s, weight]);
            let integral = defint(cx, body, s, left, right)?;
            let integral = cx.simplify(integral);
            terms.push(mul(cx.graph, &[inverse_two_c, integral]));
        }
        if kind != 0 && !cx.is_zero(f) {
            let f_s = cx.graph.substitute(f_e, x, s);
            let r = radius(cx, t)?;
            let j1 = bessel(cx, 1, r)?;
            let ratio = div(cx, j1, r);
            let body = mul(cx.graph, &[f_s, ratio]);
            let integral = defint(cx, body, s, left, right)?;
            let integral = cx.simplify(integral);
            let tm = mul(cx.graph, &[half, t, m, integral]);
            terms.push(if kind > 0 { neg(cx.graph, tm) } else { tm });
        }
        if !cx.is_zero(forcing) {
            let elapsed = sub(cx.graph, t, tau);
            let reach = mul(cx.graph, &[c, elapsed]);
            let lo = sub(cx.graph, x, reach);
            let hi = add(cx.graph, &[x, reach]);
            let f_s = cx.graph.substitute(forcing_e, x, s);
            let f_s = cx.graph.substitute(f_s, t, tau);
            let r = radius(cx, elapsed)?;
            let weight = if kind == 0 { cx.graph.int(1) } else { bessel(cx, 0, r)? };
            let body = mul(cx.graph, &[f_s, weight]);
            let inner = defint(cx, body, s, lo, hi)?;
            let inner = cx.simplify(inner);
            let outer = defint(cx, inner, tau, zero, t)?;
            terms.push(mul(cx.graph, &[inverse_two_c, outer]));
        }
    } else {
        // Kirchhoff / Poisson: ∂_t (t M[f](ct)) + t M[g](ct) + Duhamel.
        if !cx.is_zero(f) {
            let m = mean(cx, p, &xs, f_e, ct)?;
            let tm = mul(cx.graph, &[t, m]);
            let tm = cx.simplify(tm);
            terms.push(derivative(cx.graph, tm, t)?);
        }
        if !cx.is_zero(g) {
            let m = mean(cx, p, &xs, g_e, ct)?;
            terms.push(mul(cx.graph, &[t, m]));
        }
        if !cx.is_zero(forcing) {
            let (tau, _) = dummy(cx, p, "tau");
            let elapsed = sub(cx.graph, t, tau);
            let radius = mul(cx.graph, &[c, elapsed]);
            let f_tau = cx.graph.substitute(forcing_e, t, tau);
            let m = mean(cx, p, &xs, f_tau, radius)?;
            let body = mul(cx.graph, &[elapsed, m]);
            terms.push(defint(cx, body, tau, zero, t)?);
        }
    }
    let total = add(cx.graph, &terms);
    let total = if cx.is_zero(gamma) {
        total
    } else {
        let minus_gamma = neg(cx.graph, gamma);
        let e = mul(cx.graph, &[minus_gamma, t]);
        let damp = call(cx, "exp", &[e])?;
        mul(cx.graph, &[damp, total])
    };
    Some(cx.simplify(total))
}
