//! Lie point symmetries.
//!
//! A first-order equation `y' = ω(x, y)` is invariant under the flow of
//! `X = ξ ∂x + η ∂y` exactly when the linearised symmetry condition
//!
//! ```text
//! η_x + (η_y - ξ_x) ω - ξ_y ω² - ξ ω_x - η ω_y = 0
//! ```
//!
//! holds. Infinitesimals are sought as bivariate polynomials of total
//! degree up to two with unknown constant coefficients (translations,
//! scalings, rotations, projective maps and their combinations — the
//! symmetries behind separable, homogeneous, linear, Bernoulli, `F(ax+by)`
//! and linear-fractional equations, and many others). Clearing denominators
//! turns the condition into a polynomial identity in `x`, `y` and the other
//! functions occurring in `ω`, treated as independent indeterminates; its
//! coefficients give a homogeneous linear system whose null space is the
//! symmetry algebra found. Each non-trivial symmetry yields the
//! integrating factor `μ = 1/(η - ξ ω)`, which makes `dy - ω dx` exact, and
//! the solution is the level set of its potential.
//!
//! For higher-order equations the same machinery finds point symmetries
//! of `y^(n) = ω(x, y, …, y^(n-1))` from the prolonged condition; reduction
//! through canonical coordinates is used for the translation and scaling
//! symmetries handled in `reduce`.

use std::collections::BTreeMap;

use num_rational::BigRational;
use num_traits::One;
use num_traits::Zero;

use super::Problem;
use super::add;
use super::implicit;
use super::integrate;
use super::mul;
use super::sub;
use crate::graph::Cx;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::op::core;
use crate::rules::calculus::derivative;
use crate::rules::poly::repr::Gens;
use crate::rules::poly::repr::Limits;
use crate::rules::poly::repr::Mono;

/// Monomials `x^i y^j` with `i + j ≤ degree`.
fn monomials(degree: u32) -> Vec<(u32, u32)> {
    let mut out = Vec::new();
    for total in 0..=degree {
        for i in (0..=total).rev() {
            out.push((i, total - i));
        }
    }
    out
}

fn power(
    graph: &mut Graph,
    base: NodeId,
    e: u32,
) -> NodeId {
    match e {
        | 0 => graph.int(1),
        | 1 => base,
        | _ => {
            let e = graph.int(i64::from(e));
            graph.node(core::POW, &[base, e])
        },
    }
}

/// The null space of a matrix over `Q` (rows of equal length).
pub(super) fn null_space(
    mut rows: Vec<Vec<BigRational>>,
    columns: usize,
) -> Vec<Vec<BigRational>> {
    let mut pivots: Vec<usize> = Vec::new();
    let mut r = 0;
    for c in 0..columns {
        let Some(p) = (r..rows.len()).find(|&i| !rows[i][c].is_zero()) else {
            continue;
        };
        rows.swap(r, p);
        let lead = rows[r][c].clone();
        for v in &mut rows[r] {
            *v /= &lead;
        }
        for i in 0..rows.len() {
            if i != r && !rows[i][c].is_zero() {
                let factor = rows[i][c].clone();
                let pivot_row = rows[r].clone();
                for (v, p) in rows[i].iter_mut().zip(&pivot_row) {
                    *v -= &factor * p;
                }
            }
        }
        pivots.push(c);
        r += 1;
        if r == rows.len() {
            break;
        }
    }
    let free: Vec<usize> = (0..columns).filter(|c| !pivots.contains(c)).collect();
    free.iter()
        .map(|&f| {
            let mut v = vec![BigRational::zero(); columns];
            v[f] = BigRational::one();
            for (i, &p) in pivots.iter().enumerate() {
                v[p] = -rows[i][f].clone();
            }
            v
        })
        .collect()
}

/// Solutions of a linear homogeneous identity `expr = 0` in the
/// `unknowns`, required to hold identically in every other generator.
/// Returns a basis of the solution space, or `None` when `expr` is not a
/// rational expression linear in the unknowns.
pub(super) fn solve_identity(
    graph: &mut Graph,
    expr: NodeId,
    unknowns: &[NodeId],
) -> Option<Vec<Vec<BigRational>>> {
    let mut gens = Gens::default();
    let indices: Vec<u32> = unknowns.iter().map(|&u| gens.index(graph, u)).collect();
    let fraction = crate::rules::poly::ratio(graph, &mut gens, expr, Limits { terms: 20_000, exponent: 64 })?;
    // Group the numerator's terms by their part free of the unknowns.
    let mut equations: BTreeMap<Mono, Vec<BigRational>> = BTreeMap::new();
    for (mono, coeff) in fraction.numer.terms() {
        let mut rest = Mono::new();
        let mut unknown = None;
        for &(g, e) in mono {
            match indices.iter().position(|&i| i == g) {
                | Some(k) if e == 1 && unknown.is_none() => unknown = Some(k),
                | Some(_) => return None,
                | None => rest.push((g, e)),
            }
        }
        let k = unknown?;
        let row = equations.entry(rest).or_insert_with(|| vec![BigRational::zero(); unknowns.len()]);
        row[k] += coeff.to_rational()?;
    }
    Some(null_space(equations.into_values().collect(), unknowns.len()))
}

/// Point symmetries `(ξ, η)` of `y' = ω(x, y)` with polynomial
/// infinitesimals of degree at most `max_degree`.
pub(super) fn first_order_symmetries(
    cx: &mut Cx<'_>,
    x: NodeId,
    y: NodeId,
    omega: NodeId,
    max_degree: u32,
) -> Vec<(NodeId, NodeId)> {
    let mut found = Vec::new();
    for degree in 0..=max_degree {
        let Some(list) = symmetries_of_degree(cx, x, y, omega, degree) else {
            continue;
        };
        for pair in list {
            if !found.iter().any(|&(a, b)| (a, b) == pair) {
                found.push(pair);
            }
        }
        if !found.is_empty() {
            break;
        }
    }
    found
}

fn symmetries_of_degree(
    cx: &mut Cx<'_>,
    x: NodeId,
    y: NodeId,
    omega: NodeId,
    degree: u32,
) -> Option<Vec<(NodeId, NodeId)>> {
    let graph = &mut *cx.graph;
    let basis = monomials(degree);
    let mut unknowns = Vec::with_capacity(2 * basis.len());
    let mut xi_terms = Vec::new();
    let mut eta_terms = Vec::new();
    let mut monomial_terms = Vec::new();
    for &(i, j) in &basis {
        let (px, py) = (power(graph, x, i), power(graph, y, j));
        monomial_terms.push(mul(graph, &[px, py]));
    }
    for (side, terms) in [("a", &mut xi_terms), ("b", &mut eta_terms)] {
        for &m in &monomial_terms {
            let symbol = graph.interner_mut().fresh_symbol(side);
            let u = graph.symbol_node(symbol);
            unknowns.push(u);
            terms.push(mul(graph, &[u, m]));
        }
    }
    let xi = add(graph, &xi_terms);
    let eta = add(graph, &eta_terms);
    let (eta_x, eta_y) = (derivative(graph, eta, x)?, derivative(graph, eta, y)?);
    let (xi_x, xi_y) = (derivative(graph, xi, x)?, derivative(graph, xi, y)?);
    let (omega_x, omega_y) = (derivative(graph, omega, x)?, derivative(graph, omega, y)?);
    let omega_x = cx.simplify(omega_x);
    let omega_y = cx.simplify(omega_y);
    let graph = &mut *cx.graph;
    // η_x + (η_y - ξ_x) ω - ξ_y ω² - ξ ω_x - η ω_y
    let gap = sub(graph, eta_y, xi_x);
    let t1 = mul(graph, &[gap, omega]);
    let two = graph.int(2);
    let omega_sq = graph.node(core::POW, &[omega, two]);
    let minus_one = graph.int(-1);
    let t2 = mul(graph, &[minus_one, xi_y, omega_sq]);
    let t3 = mul(graph, &[minus_one, xi, omega_x]);
    let t4 = mul(graph, &[minus_one, eta, omega_y]);
    let condition = add(graph, &[eta_x, t1, t2, t3, t4]);
    let solutions = solve_identity(graph, condition, &unknowns)?;
    let mut out = Vec::new();
    for v in solutions {
        let with_values = |graph: &mut Graph, f: NodeId| -> NodeId {
            let mut result = f;
            for (u, value) in unknowns.iter().zip(&v) {
                let value = graph.num(Number::rat(value.clone()));
                result = graph.substitute(result, *u, value);
            }
            result
        };
        let (xi_v, eta_v) = (with_values(cx.graph, xi), with_values(cx.graph, eta));
        let xi_v = cx.simplify(xi_v);
        let eta_v = cx.simplify(eta_v);
        out.push((xi_v, eta_v));
    }
    Some(out)
}

/// Solves `y' = ω` with the help of a symmetry `(ξ, η)`: the integrating
/// factor `1/(η - ξ ω)` makes `dy - ω dx` exact.
pub(super) fn integrate_with_symmetry(
    cx: &mut Cx<'_>,
    problem: &mut Problem,
    omega: NodeId,
    xi: NodeId,
    eta: NodeId,
) -> Option<NodeId> {
    let (y, x) = (*problem.stand.first()?, problem.x);
    let graph = &mut *cx.graph;
    let characteristic = {
        let t = mul(graph, &[xi, omega]);
        sub(graph, eta, t)
    };
    let characteristic = cx.simplify(characteristic);
    if cx.is_zero(characteristic) {
        return None;
    }
    let minus_one = cx.graph.int(-1);
    let mu = cx.graph.node(core::POW, &[characteristic, minus_one]);
    let mu = cx.simplify(mu);
    // Φ_y = μ, Φ_x = -μ ω.
    let mu_omega = mul(cx.graph, &[minus_one, mu, omega]);
    let mu_omega = cx.simplify(mu_omega);
    let potential = potential(cx, mu, mu_omega, y, x).or_else(|| potential(cx, mu_omega, mu, x, y))?;
    let c = problem.constant(cx.graph);
    let relation = sub(cx.graph, potential, c);
    implicit(cx, problem, relation)
}

/// `Φ` with `Φ_u = f` and `Φ_v = g`: `∫ f du + ∫ (g - ∂_v ∫ f du) dv`.
pub(super) fn potential(
    cx: &mut Cx<'_>,
    f: NodeId,
    g: NodeId,
    u: NodeId,
    v: NodeId,
) -> Option<NodeId> {
    let u_symbol = cx.graph.symbol_of(u)?;
    let part = integrate(cx, f, u)?;
    let part_v = derivative(cx.graph, part, v)?;
    let rest = sub(cx.graph, g, part_v);
    let mut rest = cx.simplify(rest);
    if cx.graph.depends_on(cx.graph.find(rest), u_symbol) {
        // Free of u in value though not in form: specialise u.
        rest = specialise_constant(cx, rest, u)?;
    }
    let rest_integral = if cx.is_zero(rest) { cx.graph.int(0) } else { integrate(cx, rest, v)? };
    let total = add(cx.graph, &[part, rest_integral]);
    Some(cx.simplify(total))
}

/// `e` with `u` set to a constant, when `e` is numerically independent of
/// `u` (checked at several points of every free symbol).
pub(super) fn specialise_constant(
    cx: &mut Cx<'_>,
    e: NodeId,
    u: NodeId,
) -> Option<NodeId> {
    let u_symbol = cx.graph.symbol_of(u)?;
    let others: Vec<_> = cx.graph.free_symbols(cx.graph.find(e)).iter().copied().filter(|&s| s != u_symbol).collect();
    for base in [0.43, 1.31] {
        let mut reference = None;
        for at in [0.27, 0.91, 1.73, 2.6] {
            let mut env = crate::graph::Env::numeric(0.0);
            for (i, &s) in others.iter().enumerate() {
                env.bind(s, base + 0.29 * f64::from(u32::try_from(i).unwrap_or(0)));
            }
            env.bind(u_symbol, at);
            let v = cx.graph.eval(e, &env)?;
            if !v.is_finite() {
                return None;
            }
            match reference {
                | None => reference = Some(v),
                | Some(r) if (v - r).abs() <= 1e-9 * r.abs().max(1.0) => {},
                | Some(_) => return None,
            }
        }
    }
    for value in [0, 1, 2, -1, 3] {
        let at = cx.graph.int(value);
        let fixed = cx.graph.substitute(e, u, at);
        let fixed = cx.simplify(fixed);
        let defined = {
            let mut env = crate::graph::Env::numeric(0.0);
            for (i, &s) in others.iter().enumerate() {
                env.bind(s, 0.43 + 0.29 * f64::from(u32::try_from(i).unwrap_or(0)));
            }
            cx.graph.eval(fixed, &env).is_some_and(f64::is_finite)
        };
        if defined {
            return Some(fixed);
        }
    }
    None
}

/// The first-order equation `y' = ω` solved by Lie's method.
pub(super) fn first_order(
    cx: &mut Cx<'_>,
    problem: &mut Problem,
    omega: NodeId,
) -> Option<NodeId> {
    let (y, x) = (*problem.stand.first()?, problem.x);
    let symmetries = first_order_symmetries(cx, x, y, omega, 2);
    for (xi, eta) in symmetries {
        let saved = problem.constants;
        if let Some(found) = integrate_with_symmetry(cx, problem, omega, xi, eta) {
            return Some(found);
        }
        problem.constants = saved;
    }
    None
}

// ---------------------------------------------------------------------------
// Higher-order equations

/// `D f = f_x + Σ p_{k+1} ∂f/∂p_k` with `p_0 = y` (the `jets` are
/// `y, y', y'', …`, extended on demand).
fn total_x(
    cx: &mut Cx<'_>,
    f: NodeId,
    x: NodeId,
    jets: &mut Vec<NodeId>,
) -> Option<NodeId> {
    let mut terms = vec![derivative(cx.graph, f, x)?];
    let count = jets.len();
    for k in 0..count {
        let partial = derivative(cx.graph, f, jets[k])?;
        let partial = cx.simplify(partial);
        if cx.is_zero(partial) {
            continue;
        }
        if k + 1 == jets.len() {
            let symbol = cx.graph.interner_mut().fresh_symbol(&format!("y{}", k + 1));
            jets.push(cx.graph.symbol_node(symbol));
        }
        terms.push(mul(cx.graph, &[jets[k + 1], partial]));
    }
    let sum = add(cx.graph, &terms);
    Some(cx.simplify(sum))
}

/// Point symmetries `(ξ, η)` of `y^(n) = ω(x, y, …, y^(n-1))` with
/// bivariate polynomial infinitesimals of degree ≤ 2, from the prolonged
/// symmetry condition `η^(n) = ξ ω_x + η ω_y + Σ_k η^(k) ω_{y^(k)}` on
/// solutions, where `η^(k) = D^k Q + ξ y^(k+1)` and `Q = η - ξ y'`.
pub(super) fn higher_order_symmetries(
    cx: &mut Cx<'_>,
    problem: &Problem,
    omega: NodeId,
) -> Option<Vec<(NodeId, NodeId)>> {
    let n = problem.order();
    let x = problem.x;
    let mut jets: Vec<NodeId> = problem.stand[..n].to_vec();
    let y = jets[0];
    let basis = monomials(2);
    let mut unknowns = Vec::new();
    let mut xi_terms = Vec::new();
    let mut eta_terms = Vec::new();
    for &(i, j) in &basis {
        let (px, py) = (power(cx.graph, x, i), power(cx.graph, y, j));
        let m = mul(cx.graph, &[px, py]);
        for terms in [&mut xi_terms, &mut eta_terms] {
            let symbol = cx.graph.interner_mut().fresh_symbol("a");
            let u = cx.graph.symbol_node(symbol);
            unknowns.push(u);
            terms.push(mul(cx.graph, &[u, m]));
        }
    }
    let xi = add(cx.graph, &xi_terms);
    let eta = add(cx.graph, &eta_terms);
    let minus = cx.graph.int(-1);
    // Q = η - ξ y'
    let q = {
        let t = mul(cx.graph, &[minus, xi, jets[1]]);
        add(cx.graph, &[eta, t])
    };
    // η^(k) for k = 1..n.
    let mut prolonged = Vec::with_capacity(n + 1);
    let mut d = q;
    for k in 1..=n {
        d = total_x(cx, d, x, &mut jets)?;
        while jets.len() <= k + 1 {
            let symbol = cx.graph.interner_mut().fresh_symbol(&format!("y{}", jets.len()));
            jets.push(cx.graph.symbol_node(symbol));
        }
        let shift = mul(cx.graph, &[xi, jets[k + 1]]);
        prolonged.push(add(cx.graph, &[d, shift]));
    }
    // Condition, then on shell: y^(n) → ω, y^(n+1) → D ω.
    let (omega_x, omega_y) = (derivative(cx.graph, omega, x)?, derivative(cx.graph, omega, y)?);
    let mut rhs = vec![mul(cx.graph, &[xi, omega_x]), mul(cx.graph, &[eta, omega_y])];
    for k in 1..n {
        let partial = derivative(cx.graph, omega, jets[k])?;
        rhs.push(mul(cx.graph, &[prolonged[k - 1], partial]));
    }
    let rhs = add(cx.graph, &rhs);
    let mut condition = sub(cx.graph, prolonged[n - 1], rhs);
    let mut omega_jets: Vec<NodeId> = problem.stand[..n].to_vec();
    let d_omega = total_x(cx, omega, x, &mut omega_jets)?;
    // D ω may reference y^(n) (as jets[n]); replace it by ω as well.
    let d_omega = cx.graph.substitute(d_omega, omega_jets.get(n).copied().unwrap_or(jets[n]), omega);
    condition = cx.graph.substitute(condition, jets[n + 1], d_omega);
    condition = cx.graph.substitute(condition, jets[n], omega);
    let solutions = solve_identity(cx.graph, condition, &unknowns)?;
    let mut out = Vec::new();
    for v in solutions {
        let mut pair = [xi, eta];
        for f in &mut pair {
            for (u, value) in unknowns.iter().zip(&v) {
                let value = cx.graph.num(Number::rat(value.clone()));
                *f = cx.graph.substitute(*f, *u, value);
            }
            *f = cx.simplify(*f);
        }
        if !(cx.is_zero(pair[0]) && cx.is_zero(pair[1])) {
            out.push((pair[0], pair[1]));
        }
    }
    Some(out)
}

/// Canonical coordinates `(r, s)` of `X = ξ ∂x + η ∂y`: `X r = 0`,
/// `X s = 1`. `r` comes from the characteristic equation `dy/dx = η/ξ`
/// (solved by the ODE solver), `s` from `∫ dx/ξ` along characteristics (or
/// `∫ dy/η` when `ξ = 0`).
fn canonical_coordinates(
    cx: &mut Cx<'_>,
    x: NodeId,
    y: NodeId,
    xi: NodeId,
    eta: NodeId,
    depth: u32,
) -> Option<(NodeId, NodeId)> {
    if cx.is_zero(xi) {
        // r = x, s = ∫ dy / η(x, y).
        let minus = cx.graph.int(-1);
        let inverse = cx.graph.node(core::POW, &[eta, minus]);
        let s = integrate(cx, inverse, y)?;
        return Some((x, s));
    }
    // dY/dx = η(x, Y)/ξ(x, Y).
    let f_symbol = cx.graph.interner_mut().fresh_symbol("Y");
    let f = cx.graph.symbol_node(f_symbol);
    let big_y = cx.graph.node(core::APPLY, &[f, x]);
    let minus = cx.graph.int(-1);
    let inverse = cx.graph.node(core::POW, &[xi, minus]);
    let slope = mul(cx.graph, &[eta, inverse]);
    let slope = cx.graph.substitute(slope, y, big_y);
    let diff = cx.graph.ops().lookup("diff")?;
    let dy = cx.graph.node(diff, &[big_y, x]);
    let equation = cx.graph.node(core::EQ, &[dy, slope]);
    let (_, answer) = super::solve_equation(cx, equation, big_y, depth + 1)?;
    // The constant of the characteristic is the invariant r.
    let c1 = cx.graph.sym("C1");
    let relation = super::super::solve::as_expression(cx.graph, answer);
    let relation = cx.graph.replace_subterm(relation, big_y, y);
    let r = *super::solve_for(cx.graph, relation, c1, 0)?.first()?;
    let r = cx.simplify(r);
    // Along a characteristic y = Y(x; r): s = ∫ dx / ξ(x, Y(x; r)).
    let r_symbol = cx.graph.interner_mut().fresh_symbol("rho");
    let rho = cx.graph.symbol_node(r_symbol);
    let level = sub(cx.graph, r, rho);
    let along = *super::solve_for(cx.graph, level, y, 0)?.first()?;
    let xi_along = cx.graph.substitute(xi, y, along);
    let integrand = cx.graph.node(core::POW, &[xi_along, minus]);
    let integrand = cx.simplify(integrand);
    let s = integrate(cx, integrand, x)?;
    let s = cx.graph.substitute(s, rho, r);
    let s = cx.simplify(s);
    Some((r, s))
}

/// Second-order equations by one point symmetry: canonical coordinates
/// turn `y'' = ω` into a first-order equation for `v = ds/dr`.
pub(super) fn second_order(
    cx: &mut Cx<'_>,
    problem: &mut Problem,
    depth: u32,
) -> Option<NodeId> {
    if problem.order() != 2 || depth > 2 {
        return None;
    }
    let (x, y, dy, ddy) = (problem.x, problem.stand[0], problem.stand[1], problem.stand[2]);
    let omega = *super::solve_for(cx.graph, problem.expr, ddy, 0)?.first()?;
    let omega = cx.simplify(omega);
    let symmetries = higher_order_symmetries(cx, problem, omega)?;
    for (xi, eta) in symmetries {
        let saved = problem.constants;
        if let Some(found) = reduce_by_symmetry(cx, problem, x, y, dy, ddy, xi, eta, depth) {
            return Some(found);
        }
        problem.constants = saved;
    }
    None
}

#[allow(clippy::too_many_arguments)] // the coordinates and stand-ins of one problem
fn reduce_by_symmetry(
    cx: &mut Cx<'_>,
    problem: &mut Problem,
    x: NodeId,
    y: NodeId,
    dy: NodeId,
    ddy: NodeId,
    xi: NodeId,
    eta: NodeId,
    depth: u32,
) -> Option<NodeId> {
    let (r, s) = canonical_coordinates(cx, x, y, xi, eta, depth)?;
    // v = (s_x + s_y y')/(r_x + r_y y'), so y' = (s_x - v r_x)/(v r_y - s_y).
    let (rx, ry, sx, sy) = (
        derivative(cx.graph, r, x)?,
        derivative(cx.graph, r, y)?,
        derivative(cx.graph, s, x)?,
        derivative(cx.graph, s, y)?,
    );
    let v_symbol = cx.graph.interner_mut().fresh_symbol("v");
    let v = cx.graph.symbol_node(v_symbol);
    let w_symbol = cx.graph.interner_mut().fresh_symbol("w");
    let w = cx.graph.symbol_node(w_symbol);
    let minus = cx.graph.int(-1);
    let y1 = {
        let top = {
            let t = mul(cx.graph, &[minus, v, rx]);
            add(cx.graph, &[sx, t])
        };
        let bottom = {
            let t = mul(cx.graph, &[v, ry]);
            let u = mul(cx.graph, &[minus, sy]);
            add(cx.graph, &[t, u])
        };
        let inverse = cx.graph.node(core::POW, &[bottom, minus]);
        let q = mul(cx.graph, &[top, inverse]);
        cx.simplify(q)
    };
    // y'' = Y1_x + Y1_y Y1 + Y1_v w (r_x + r_y Y1)
    let y2 = {
        let a = derivative(cx.graph, y1, x)?;
        let b = derivative(cx.graph, y1, y)?;
        let c = derivative(cx.graph, y1, v)?;
        let speed = {
            let t = mul(cx.graph, &[ry, y1]);
            add(cx.graph, &[rx, t])
        };
        let bt = mul(cx.graph, &[b, y1]);
        let ct = mul(cx.graph, &[c, w, speed]);
        let sum = add(cx.graph, &[a, bt, ct]);
        cx.simplify(sum)
    };
    let mut reduced = cx.graph.substitute(problem.expr, ddy, y2);
    reduced = cx.graph.substitute(reduced, dy, y1);
    // Express x, y through (R, S) and fix S: the equation does not depend
    // on S (checked numerically at two values).
    let big_r_symbol = cx.graph.interner_mut().fresh_symbol("R");
    let big_r = cx.graph.symbol_node(big_r_symbol);
    let level_r = sub(cx.graph, r, big_r);
    let y_of = *super::solve_for(cx.graph, level_r, y, 0)?.first()?;
    let s_on = cx.graph.substitute(s, y, y_of);
    let in_rs = |cx: &mut Cx<'_>, s0: i64| -> Option<NodeId> {
        let target = cx.graph.int(s0);
        let level_s = sub(cx.graph, s_on, target);
        let x_of = *super::solve_for(cx.graph, level_s, x, 0)?.first()?;
        let y_at = cx.graph.substitute(y_of, x, x_of);
        let e = cx.graph.substitute(reduced, y, y_at);
        let e = cx.graph.substitute(e, x, x_of);
        Some(cx.simplify(e))
    };
    let (e0, e1) = (in_rs(cx, 0)?, in_rs(cx, 1)?);
    for point in [0.37, 0.81, 1.43] {
        let mut env = crate::graph::Env::numeric(0.0);
        env.bind(big_r_symbol, point);
        env.bind(v_symbol, 0.6 + point);
        env.bind(w_symbol, 0.3 + point);
        let others: Vec<_> = cx.graph.free_symbols(cx.graph.find(e0)).to_vec();
        for sym in others {
            if sym != big_r_symbol && sym != v_symbol && sym != w_symbol {
                env.bind(sym, 0.9);
            }
        }
        let (a, b) = (cx.graph.eval(e0, &env)?, cx.graph.eval(e1, &env)?);
        // Invariance up to a common non-zero factor.
        if (a.abs() < 1e-12) != (b.abs() < 1e-12) {
            return None;
        }
    }
    // First-order equation for V(R).
    let f_symbol = cx.graph.interner_mut().fresh_symbol("V");
    let f = cx.graph.symbol_node(f_symbol);
    let big_v = cx.graph.node(core::APPLY, &[f, big_r]);
    let diff = cx.graph.ops().lookup("diff")?;
    let dv = cx.graph.node(diff, &[big_v, big_r]);
    let first = cx.graph.substitute(e0, w, dv);
    let first = cx.graph.substitute(first, v, big_v);
    let zero = cx.graph.int(0);
    let equation = cx.graph.node(core::EQ, &[first, zero]);
    let (_, answer) = super::solve_equation(cx, equation, big_v, depth + 1)?;
    let &[lhs, slope] = cx.graph.children(answer) else {
        return None;
    };
    if lhs != big_v {
        return None;
    }
    problem.constants = problem.constants.max(constants_in(cx, slope));
    // s = ∫ V(R) dR + C, R = r(x, y).
    let integral = integrate(cx, slope, big_r)?;
    let c = problem.constant(cx.graph);
    let integral = cx.graph.substitute(integral, big_r, r);
    let rhs = add(cx.graph, &[integral, c]);
    let relation = sub(cx.graph, s, rhs);
    implicit(cx, problem, relation)
}

fn constants_in(
    cx: &Cx<'_>,
    node: NodeId,
) -> usize {
    cx.graph
        .free_symbols(cx.graph.find(node))
        .iter()
        .filter_map(|&s| cx.graph.interner().symbol_name(s).strip_prefix('C')?.parse::<usize>().ok())
        .max()
        .unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn null_space_of_small_systems() {
        let q = |v: i64| BigRational::from_integer(v.into());
        let basis = null_space(vec![vec![q(1), q(1), q(0)], vec![q(0), q(0), q(1)]], 3);
        assert_eq!(basis, vec![vec![q(-1), q(1), q(0)]]);
        assert!(null_space(vec![vec![q(1), q(0)], vec![q(0), q(1)]], 2).is_empty());
    }
}
