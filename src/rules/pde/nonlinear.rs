//! Nonlinear scalar equations in one space variable.
//!
//! * **Viscous Burgers by the Hopf–Cole transformation.** For
//!   `A u_t + (a u + b) A u_x = ν A u_xx` with constants (the Burgers
//!   equation in a Galilean frame, scaled), `U = a u + b` solves
//!   `U_t + U U_x = ν U_xx`; `U = -2ν φ_x/φ` turns it into the heat equation
//!   `φ_t = ν φ_xx` with `φ(x, 0) = exp(-(1/2ν) ∫_0^x U₀)`, so that
//!   `U = ∫ ((x - s)/t) e^{-G/2ν} ds / ∫ e^{-G/2ν} ds`,
//!   `G = (x - s)²/(2t) + ∫_0^s U₀`.
//! * **Riemann problems of scalar conservation laws.** For
//!   `u_t + f'(u) u_x = 0` with the step data `u(x, 0) = u_l + (u_r - u_l)
//!   H(x - x₀)` and a convex or concave flux: the Rankine–Hugoniot shock
//!   `s = (f(u_r) - f(u_l))/(u_r - u_l)` when the Oleinik entropy condition
//!   holds (`u_l > u_r` for convex `f`, `u_l < u_r` for concave `f`), the
//!   rarefaction fan `f'(u) = (x - x₀)/t` otherwise.
//! * **Solitons.** `pde_kdv_soliton(list(κ…), list(x₀…))`: the
//!   reflectionless N-soliton solution of `u_t - 6 u u_x + u_xxx = 0` from
//!   the inverse scattering transform, `u = -2 ∂_x² ln det(I + A)`,
//!   `A_{mn} = 2 √(κ_m κ_n)/(κ_m + κ_n) e^{κ_m x₀_m + κ_n x₀_n -
//!   (κ_m + κ_n) x + 4 (κ_m³ + κ_n³) t}`;
//!   `pde_sine_gordon_soliton(list(a…), list(c…))`: the N-kink solution of
//!   `u_xt = sin u` from the Bäcklund transformation
//!   `u₁ = 4 arctan e^{a x + t/a + c}` of the vacuum and the Bianchi
//!   permutability `tan((u₁₂ - u₀)/4) = (a₁ + a₂)/(a₁ - a₂) tan((u₁ - u₂)/4)`.

use super::Conditions;
use super::Problem;
use super::derivative_by_node;
use super::dummy;
use super::modes::exact_value;
use super::util::call;
use super::util::contains_op;
use super::util::defint;
use super::util::div;
use super::util::infinity;
use crate::graph::Cx;
use crate::graph::NodeId;
use crate::rules::calculus::antiderivative;
use crate::rules::calculus::derivative;
use crate::rules::complex::build::add;
use crate::rules::complex::build::mul;
use crate::rules::complex::build::neg;
use crate::rules::complex::build::powi;
use crate::rules::complex::build::sub;
use crate::rules::solve::solve_for;

/// `A u_t + B(u) u_x + C u_xx`, read off the residual.
struct Law {
    /// `B/A`, a function of the unknown.
    speed: NodeId,
    /// `-C/A` if the equation has a second derivative.
    viscosity: Option<NodeId>,
    x: NodeId,
    t: NodeId,
}

fn decompose(
    cx: &mut Cx<'_>,
    p: &Problem,
) -> Option<Law> {
    if p.dimension() != 2 || !p.homogeneous(cx.graph) {
        return None;
    }
    let time = p.time(cx.graph);
    let space = 1 - time;
    let (u_t, u_x, u_xx) = (p.unit(time, 1), p.unit(space, 1), p.unit(space, 2));
    let jet = |index: &Vec<u32>| p.jets.iter().find(|(i, _)| i == index).map(|(_, n)| *n);
    let (t_node, x_node) = (jet(&u_t)?, jet(&u_x)?);
    let xx_node = jet(&u_xx);
    if p.jets.iter().any(|(i, _)| *i != u_t && *i != u_x && *i != u_xx && i.iter().any(|&d| d != 0)) {
        return None;
    }
    let a = derivative_by_node(cx, p.residual, t_node)?;
    let b = derivative_by_node(cx, p.residual, x_node)?;
    let c = match xx_node {
        | Some(n) => derivative_by_node(cx, p.residual, n)?,
        | None => cx.graph.int(0),
    };
    let u_symbol = cx.graph.symbol_of(p.function)?;
    let free_of_u = |cx: &Cx<'_>, n: NodeId| !cx.graph.depends_on(cx.graph.find(n), u_symbol);
    if cx.is_zero(a) || !p.constant(cx.graph, a) || !free_of_u(cx, a) || !p.constant(cx.graph, c) || !free_of_u(cx, c) {
        return None;
    }
    // The residual is exactly A u_t + B u_x + C u_xx.
    let mut parts = vec![mul(cx.graph, &[a, t_node]), mul(cx.graph, &[b, x_node])];
    if let Some(n) = xx_node {
        parts.push(mul(cx.graph, &[c, n]));
    }
    let rebuilt = add(cx.graph, &parts);
    let gap = sub(cx.graph, rebuilt, p.residual);
    if !cx.is_zero(gap) {
        return None;
    }
    let speed = div(cx, b, a);
    let speed = cx.simplify(speed);
    let viscosity = if cx.is_zero(c) {
        None
    } else {
        let q = div(cx, c, a);
        let q = neg(cx.graph, q);
        Some(cx.simplify(q))
    };
    Some(Law { speed, viscosity, x: p.vars[space], t: p.vars[time] })
}

/// The nonlinear solvers of this module.
pub(super) fn solve(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
) -> Option<NodeId> {
    let law = decompose(cx, p)?;
    let time = p.time(cx.graph);
    let zero_index = vec![0; 2];
    let initial = conditions.find(cx.graph, time, &zero_index, None)?;
    if conditions.0.len() != 1 || !cx.is_zero(initial.point) {
        return None;
    }
    match law.viscosity {
        | Some(nu) => hopf_cole(cx, p, &law, nu, initial.value),
        | None => riemann(cx, p, &law, initial.value),
    }
}

fn hopf_cole(
    cx: &mut Cx<'_>,
    p: &Problem,
    law: &Law,
    nu: NodeId,
    u0: NodeId,
) -> Option<NodeId> {
    if !p.constant(cx.graph, nu) || super::util::sample(cx.graph, nu, 0)? <= 0.0 {
        return None;
    }
    let u_symbol = cx.graph.symbol_of(p.function)?;
    // speed = a u + b.
    let a = derivative_by_node(cx, law.speed, p.unknown)?;
    if cx.graph.depends_on(cx.graph.find(a), u_symbol) || cx.is_zero(a) || !p.constant(cx.graph, a) {
        return None;
    }
    let au = mul(cx.graph, &[a, p.unknown]);
    let b = sub(cx.graph, law.speed, au);
    let b = cx.simplify(b);
    if cx.graph.depends_on(cx.graph.find(b), u_symbol) || !p.constant(cx.graph, b) {
        return None;
    }
    let (x, t) = (law.x, law.t);
    let (s, _) = dummy(cx, p, "s");
    // U₀ = a u₀ + b and F₀(s) = ∫_0^s U₀.
    let au0 = mul(cx.graph, &[a, u0]);
    let big_u0 = add(cx.graph, &[au0, b]);
    let big_u0 = cx.simplify(big_u0);
    let primitive = antiderivative(cx, big_u0, x)?;
    let zero = cx.graph.int(0);
    let at_s = cx.graph.substitute(primitive, x, s);
    let at_0 = cx.graph.substitute(primitive, x, zero);
    let f0 = sub(cx.graph, at_s, at_0);
    // G = (x - s)²/(2 t) + F₀(s), weight exp(-G/(2ν)).
    let gap = sub(cx.graph, x, s);
    let two = cx.graph.int(2);
    let two_t = mul(cx.graph, &[two, t]);
    let square = powi(cx.graph, gap, 2);
    let kinetic = div(cx, square, two_t);
    let g = add(cx.graph, &[kinetic, f0]);
    let two_nu = mul(cx.graph, &[two, nu]);
    let exponent = {
        let q = div(cx, g, two_nu);
        neg(cx.graph, q)
    };
    let weight = call(cx, "exp", &[exponent])?;
    let oo = infinity(cx)?;
    let minus_oo = neg(cx.graph, oo);
    let drift = div(cx, gap, t);
    let numerator_body = mul(cx.graph, &[drift, weight]);
    let numerator = defint(cx, numerator_body, s, minus_oo, oo)?;
    let denominator = defint(cx, weight, s, minus_oo, oo)?;
    let big_u = div(cx, numerator, denominator);
    let shifted = sub(cx.graph, big_u, b);
    let u = div(cx, shifted, a);
    Some(cx.simplify(u))
}

/// The heaviside atoms of `value`.
fn heavisides(
    cx: &Cx<'_>,
    value: NodeId,
) -> Vec<NodeId> {
    let Some(op) = cx.graph.ops().lookup("heaviside") else {
        return Vec::new();
    };
    let mut out = Vec::new();
    let mut stack = vec![value];
    while let Some(n) = stack.pop() {
        if cx.graph.op(n) == op {
            if !out.contains(&n) {
                out.push(n);
            }
            continue;
        }
        stack.extend_from_slice(cx.graph.children(n));
    }
    out
}

#[allow(clippy::too_many_lines)]
fn riemann(
    cx: &mut Cx<'_>,
    p: &Problem,
    law: &Law,
    data: NodeId,
) -> Option<NodeId> {
    let (x, t) = (law.x, law.t);
    let steps = heavisides(cx, data);
    let [step] = steps.as_slice() else {
        return None;
    };
    let argument = *cx.graph.children(*step).first()?;
    let slope = derivative(cx.graph, argument, x)?;
    let slope_value = exact_value(cx, slope)?;
    if slope_value == 0.0 {
        return None;
    }
    let zero = cx.graph.int(0);
    let at_zero = cx.graph.substitute(argument, x, zero);
    let x0 = {
        let q = div(cx, at_zero, slope);
        let q = neg(cx.graph, q);
        cx.simplify(q)
    };
    let (zero_value, one_value) = (cx.graph.int(0), cx.graph.int(1));
    let low = cx.graph.replace_subterm(data, *step, zero_value);
    let high = cx.graph.replace_subterm(data, *step, one_value);
    let (low, high) = (cx.simplify(low), cx.simplify(high));
    let (ul, ur) = if slope_value > 0.0 { (low, high) } else { (high, low) };
    let (vl, vr) = (exact_value(cx, ul)?, exact_value(cx, ur)?);
    if (vl - vr).abs() < 1e-14 {
        return Some(ul);
    }
    // The flux derivative as a function of a dummy w.
    let (w, _) = dummy(cx, p, "w");
    let speed_w = cx.graph.replace_subterm(law.speed, p.unknown, w);
    // Convexity of the flux between the states.
    let curvature = derivative(cx.graph, speed_w, w)?;
    let w_symbol = cx.graph.symbol_of(w)?;
    let mut sign = 0.0_f64;
    for k in 0..=8_u32 {
        let value = vl + (vr - vl) * f64::from(k) / 8.0;
        let mut env = crate::graph::Env::numeric(0.0);
        env.bind(w_symbol, value);
        let c = cx.graph.eval(curvature, &env).filter(|v| v.is_finite())?;
        if c.abs() < 1e-12 {
            return None;
        }
        if sign == 0.0 {
            sign = c.signum();
        } else if c.signum() != sign {
            return None;
        }
    }
    let shock = (sign > 0.0) == (vl > vr);
    let shifted = sub(cx.graph, x, x0);
    if shock {
        // s = (f(u_r) - f(u_l))/(u_r - u_l) with f' = speed.
        let flux = antiderivative(cx, speed_w, w)?;
        let at = |cx: &mut Cx<'_>, value: NodeId| cx.graph.substitute(flux, w, value);
        let (fl, fr) = (at(cx, ul), at(cx, ur));
        let jump_f = sub(cx.graph, fr, fl);
        let jump_u = sub(cx.graph, ur, ul);
        let s = div(cx, jump_f, jump_u);
        let st = mul(cx.graph, &[s, t]);
        let front = sub(cx.graph, shifted, st);
        let h = call(cx, "heaviside", &[front])?;
        let step_up = sub(cx.graph, ur, ul);
        let moved = mul(cx.graph, &[step_up, h]);
        let u = add(cx.graph, &[ul, moved]);
        return Some(cx.simplify(u));
    }
    // Rarefaction fan: f'(w) = (x - x₀)/t.
    let (a_speed, b_speed) = (cx.graph.substitute(speed_w, w, ul), cx.graph.substitute(speed_w, w, ur));
    let xi = div(cx, shifted, t);
    let equation = sub(cx.graph, speed_w, xi);
    let roots = solve_for(cx.graph, equation, w, 0)?;
    let &[fan] = roots.as_slice() else {
        return None;
    };
    let edge = |cx: &mut Cx<'_>, c: NodeId| {
        let ct = mul(cx.graph, &[c, t]);
        sub(cx.graph, shifted, ct)
    };
    let (left_edge, right_edge) = (edge(cx, a_speed), edge(cx, b_speed));
    let (h_left, h_right) = (call(cx, "heaviside", &[left_edge])?, call(cx, "heaviside", &[right_edge])?);
    let one = cx.graph.int(1);
    let not_left = sub(cx.graph, one, h_left);
    let not_right = sub(cx.graph, one, h_right);
    let t1 = mul(cx.graph, &[ul, not_left]);
    let t2 = mul(cx.graph, &[fan, h_left, not_right]);
    let t3 = mul(cx.graph, &[ur, h_right]);
    let u = add(cx.graph, &[t1, t2, t3]);
    Some(cx.simplify(u))
}

// ----------------------------------------------------------------------
// Solitons
// ----------------------------------------------------------------------

/// The determinant of a small matrix by cofactor expansion.
fn determinant(
    cx: &mut Cx<'_>,
    m: &[Vec<NodeId>],
) -> NodeId {
    let n = m.len();
    if n == 1 {
        return m[0][0];
    }
    let mut terms = Vec::new();
    for j in 0..n {
        let minor: Vec<Vec<NodeId>> = m.iter().skip(1).map(|row| row.iter().enumerate().filter(|(k, _)| *k != j).map(|(_, &e)| e).collect()).collect();
        let d = determinant(cx, &minor);
        let term = mul(cx.graph, &[m[0][j], d]);
        terms.push(if j % 2 == 0 { term } else { neg(cx.graph, term) });
    }
    add(cx.graph, &terms)
}

/// `pde_kdv_soliton(list(κ…), list(x₀…))`.
pub(super) fn kdv_soliton(
    cx: &mut Cx<'_>,
    kappas: &[NodeId],
    shifts: &[NodeId],
) -> Option<NodeId> {
    let n = kappas.len();
    if n == 0 || n > 4 || shifts.len() != n {
        return None;
    }
    let (x, t) = (cx.graph.sym("x"), cx.graph.sym("t"));
    let (two, four) = (cx.graph.int(2), cx.graph.int(4));
    let mut matrix = Vec::new();
    for m in 0..n {
        let mut row = Vec::new();
        for k in 0..n {
            let (km, kn) = (kappas[m], kappas[k]);
            let sum = add(cx.graph, &[km, kn]);
            let product = mul(cx.graph, &[km, kn]);
            let root = super::util::sqrt(cx, product)?;
            let amplitude = {
                let numerator = mul(cx.graph, &[two, root]);
                div(cx, numerator, sum)
            };
            let shift_m = mul(cx.graph, &[km, shifts[m]]);
            let shift_n = mul(cx.graph, &[kn, shifts[k]]);
            let spatial = {
                let s = mul(cx.graph, &[sum, x]);
                neg(cx.graph, s)
            };
            let (km3, kn3) = (powi(cx.graph, km, 3), powi(cx.graph, kn, 3));
            let cubes = add(cx.graph, &[km3, kn3]);
            let temporal = mul(cx.graph, &[four, cubes, t]);
            let exponent = add(cx.graph, &[shift_m, shift_n, spatial, temporal]);
            let e = call(cx, "exp", &[exponent])?;
            let entry = mul(cx.graph, &[amplitude, e]);
            let one = cx.graph.int(1);
            row.push(if m == k { add(cx.graph, &[one, entry]) } else { entry });
        }
        matrix.push(row);
    }
    let d = determinant(cx, &matrix);
    let log = call(cx, "ln", &[d])?;
    let first = derivative(cx.graph, log, x)?;
    let second = derivative(cx.graph, first, x)?;
    let minus_two = cx.graph.int(-2);
    let u = mul(cx.graph, &[minus_two, second]);
    Some(cx.simplify(u))
}

/// `4 arctan e^{a x + t/a + c}`.
fn kink(
    cx: &mut Cx<'_>,
    a: NodeId,
    c: NodeId,
) -> Option<NodeId> {
    let (x, t) = (cx.graph.sym("x"), cx.graph.sym("t"));
    let ax = mul(cx.graph, &[a, x]);
    let t_over_a = div(cx, t, a);
    let e = add(cx.graph, &[ax, t_over_a, c]);
    let exp = call(cx, "exp", &[e])?;
    let arctan = call(cx, "atan", &[exp])?;
    let four = cx.graph.int(4);
    Some(mul(cx.graph, &[four, arctan]))
}

/// The N-kink solution of `u_xt = sin u` by repeated Bianchi permutability.
fn bianchi(
    cx: &mut Cx<'_>,
    a: &[NodeId],
    c: &[NodeId],
) -> Option<NodeId> {
    match a.len() {
        | 0 => Some(cx.graph.int(0)),
        | 1 => kink(cx, a[0], c[0]),
        | n => {
            // Remove the last two parameters j, k.
            let (j, k) = (n - 2, n - 1);
            let rest_a: Vec<NodeId> = a.iter().take(j).copied().collect();
            let rest_c: Vec<NodeId> = c.iter().take(j).copied().collect();
            let with = |list: &[NodeId], extra: NodeId| {
                let mut v = list.to_vec();
                v.push(extra);
                v
            };
            let u_s = bianchi(cx, &rest_a, &rest_c)?;
            let u_j = bianchi(cx, &with(&rest_a, a[j]), &with(&rest_c, c[j]))?;
            let u_k = bianchi(cx, &with(&rest_a, a[k]), &with(&rest_c, c[k]))?;
            let (plus, minus) = (add(cx.graph, &[a[j], a[k]]), sub(cx.graph, a[j], a[k]));
            let ratio = div(cx, plus, minus);
            let diff = sub(cx.graph, u_j, u_k);
            let quarter = {
                let four = cx.graph.int(4);
                div(cx, diff, four)
            };
            let tan = call(cx, "tan", &[quarter])?;
            let inner = mul(cx.graph, &[ratio, tan]);
            let arctan = call(cx, "atan", &[inner])?;
            let four = cx.graph.int(4);
            let correction = mul(cx.graph, &[four, arctan]);
            Some(add(cx.graph, &[u_s, correction]))
        },
    }
}

/// `pde_sine_gordon_soliton(list(a…), list(c…))`.
pub(super) fn sine_gordon_soliton(
    cx: &mut Cx<'_>,
    a: &[NodeId],
    c: &[NodeId],
) -> Option<NodeId> {
    if a.is_empty() || a.len() > 3 || a.len() != c.len() {
        return None;
    }
    let u = bianchi(cx, a, c)?;
    if contains_op(cx.graph, u, "diff") {
        return None;
    }
    Some(cx.simplify(u))
}
