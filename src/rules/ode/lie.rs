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
fn potential(
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
fn specialise_constant(
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
