//! Systems of equations `solve(list(equations), list(unknowns))`.
//!
//! In order:
//!
//! 1. a square linear system by Cramer's rule;
//! 2. a polynomial system by a lexicographic Gröbner basis. Symbols other
//!    than the unknowns are parameters (ordered last), and when the system
//!    is underdetermined a subset of the unknowns is left *free*: the
//!    unknowns that remain are eliminated and expressed through the free
//!    ones, and a free unknown appears in the answer as itself. The basis
//!    is triangular, so the solution is read off by back-substitution;
//! 3. two polynomial equations in two unknowns by the resultant with
//!    respect to one of them;
//! 4. systems with transcendental parts by elimination: an equation that
//!    is linear in an unknown determines it, otherwise an unknown is
//!    solved for with the scalar solver and the branches are followed.
//!
//! Every solution tuple is substituted back into all the equations.

use super::as_expression;
use super::back_substitute;
use super::sample_envs;
use super::solve_for;
use super::solve_linear;
use crate::graph::op::core;
use crate::graph::Env;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::SymbolId;
use crate::rules::poly::best;
use crate::rules::poly::from_groebner;
use crate::rules::poly::groebner::groebner;
use crate::rules::poly::groebner::GroebnerLimits;
use crate::rules::poly::groebner::Order;
use crate::rules::poly::repr::from_term;
use crate::rules::poly::repr::to_term;
use crate::rules::poly::repr::Gens;
use crate::rules::poly::repr::Limits;
use crate::rules::poly::to_groebner;

/// Solution tuples of `equations = 0` for `unknowns`, verified.
pub(super) fn solve_system(
    graph: &mut Graph,
    equations: &[NodeId],
    unknowns: &[NodeId],
) -> Option<Vec<Vec<NodeId>>> {
    if unknowns.is_empty() {
        return None;
    }
    if equations.len() == unknowns.len() {
        if let Some(single) = solve_linear(graph, equations, unknowns) {
            return Some(vec![single]);
        }
    }
    let exprs: Vec<NodeId> = equations
        .iter()
        .map(|&e| {
            let expr = as_expression(graph, e);
            best(graph, expr)
        })
        .collect::<Option<_>>()?;
    let found = polynomial(graph, &exprs, unknowns)
        .or_else(|| resultant_pair(graph, &exprs, unknowns))
        .or_else(|| eliminate(graph, &exprs, unknowns, 0))?;
    let mut out: Vec<Vec<NodeId>> = Vec::new();
    for tuple in found {
        if !holds(graph, &exprs, unknowns, &tuple) {
            continue;
        }
        if !out.iter().any(|o| o.iter().zip(&tuple).all(|(&a, &b)| a == b || graph.same(a, b))) {
            out.push(tuple);
        }
    }
    Some(out)
}

/// Whether `tuple` makes every expression vanish: exactly (numerically)
/// when there are no free symbols, at several spot checks otherwise.
fn holds(
    graph: &mut Graph,
    exprs: &[NodeId],
    unknowns: &[NodeId],
    tuple: &[NodeId],
) -> bool {
    let mut substituted: Vec<NodeId> = Vec::with_capacity(exprs.len());
    for &e in exprs {
        let mut term = e;
        for (&u, &v) in unknowns.iter().zip(tuple) {
            term = graph.substitute(term, u, v);
        }
        substituted.push(term);
    }
    let mut all = substituted.clone();
    all.extend_from_slice(tuple);
    let has_free = substituted.iter().any(|&t| !graph.free_symbols(graph.find(t)).is_empty());
    if !has_free {
        return substituted.iter().all(|&t| match graph.eval(t, &Env::numeric(0.0)) {
            | Some(v) => v.abs() <= 1e-7 * (1.0 + tuple.iter().filter_map(|&n| graph.eval(n, &Env::numeric(0.0))).map(f64::abs).fold(0.0, f64::max).powi(3)),
            | None => true,
        });
    }
    for env in sample_envs(graph, &all, None) {
        let magnitude = tuple.iter().filter_map(|&n| graph.eval(n, &env)).filter(|v| v.is_finite()).map(f64::abs).fold(0.0, f64::max);
        for &t in &substituted {
            if let Some(v) = graph.eval(t, &env) {
                if v.is_finite() && v.abs() > 1e-7 * (1.0 + magnitude.powi(3)) {
                    return false;
                }
            }
        }
    }
    true
}

fn unknown_symbols(
    graph: &Graph,
    unknowns: &[NodeId],
) -> Option<Vec<SymbolId>> {
    unknowns.iter().map(|&u| graph.symbol_of(u)).collect()
}

/// The free symbols of `exprs` that are not unknowns.
fn parameters(
    graph: &mut Graph,
    exprs: &[NodeId],
    unknowns: &[NodeId],
) -> Option<Vec<NodeId>> {
    let symbols = unknown_symbols(graph, unknowns)?;
    let mut out: Vec<SymbolId> = Vec::new();
    for &e in exprs {
        for &s in graph.free_symbols(graph.find(e)) {
            if !symbols.contains(&s) && !out.contains(&s) {
                out.push(s);
            }
        }
    }
    Some(out.into_iter().map(|s| graph.symbol_node(s)).collect())
}

/// Subsets of `0..n` of size `k`, preferring late indices.
fn subsets(
    n: usize,
    k: usize,
) -> Vec<Vec<usize>> {
    let mut out: Vec<Vec<usize>> = Vec::new();
    for mask in 0_u32..(1 << n) {
        if mask.count_ones() as usize == k {
            out.push((0..n).filter(|&i| mask >> i & 1 == 1).collect());
        }
    }
    out.sort_by_key(|s| std::cmp::Reverse(s.clone()));
    out
}

/// A polynomial system through a lexicographic Gröbner basis, with free
/// unknowns when it is underdetermined.
fn polynomial(
    graph: &mut Graph,
    exprs: &[NodeId],
    unknowns: &[NodeId],
) -> Option<Vec<Vec<NodeId>>> {
    let n = unknowns.len();
    let params = parameters(graph, exprs, unknowns)?;
    if n > 6 || params.len() > 4 {
        return None;
    }
    let list = graph.node(core::LIST, exprs);
    for k in 0..n {
        for free in subsets(n, k) {
            let dependent: Vec<usize> = (0..n).filter(|i| !free.contains(i)).collect();
            let mut variables: Vec<NodeId> = dependent.iter().map(|&i| unknowns[i]).collect();
            variables.extend(free.iter().map(|&i| unknowns[i]));
            variables.extend_from_slice(&params);
            let vars = graph.node(core::LIST, &variables);
            let (generators, gens) = to_groebner(graph, list, vars, Order::Lex)?;
            let basis = groebner(&generators, Order::Lex, GroebnerLimits::default())?;
            // An inconsistent system.
            if basis.iter().any(|g| g.leading().is_some_and(|(m, _)| m.iter().all(|&e| e == 0))) {
                return Some(Vec::new());
            }
            let nd = dependent.len();
            let finite = (0..nd).all(|d| {
                basis.iter().any(|g| {
                    g.leading().is_some_and(|(m, _)| {
                        m.get(d).is_some_and(|&e| e > 0) && m.iter().take(nd).enumerate().all(|(j, &e)| j == d || e == 0)
                    })
                })
            });
            let constrained = basis.iter().any(|g| g.leading().is_some_and(|(m, _)| m.iter().take(nd).all(|&e| e == 0)));
            if !finite || constrained {
                continue;
            }
            let terms: Vec<NodeId> = basis.iter().map(|g| from_groebner(graph, &gens, g)).collect();
            let dep_nodes: Vec<NodeId> = dependent.iter().map(|&i| unknowns[i]).collect();
            let Some(solutions) = back_substitute(graph, &terms, &dep_nodes) else {
                continue;
            };
            let mut tuples = Vec::with_capacity(solutions.len());
            for solution in solutions {
                let mut tuple: Vec<NodeId> = unknowns.to_vec();
                for (&i, &value) in dependent.iter().zip(&solution) {
                    *tuple.get_mut(i)? = value;
                }
                tuples.push(tuple);
            }
            return Some(tuples);
        }
    }
    None
}

/// Two polynomial equations in two unknowns: the resultant with respect
/// to the second unknown gives the first, and the first equation the
/// second.
fn resultant_pair(
    graph: &mut Graph,
    exprs: &[NodeId],
    unknowns: &[NodeId],
) -> Option<Vec<Vec<NodeId>>> {
    let &[x, y] = unknowns else {
        return None;
    };
    let limits = Limits::default();
    let mut gens = Gens::default();
    let gx = gens.index(graph, x);
    let gy = gens.index(graph, y);
    let (symbol_x, symbol_y) = (graph.symbol_of(x)?, graph.symbol_of(y)?);
    let mut polys = Vec::new();
    for &e in exprs {
        let p = from_term(graph, &mut gens, e, limits)?;
        for g in p.support() {
            let node = gens.node(g)?;
            if g != gx && g != gy && (graph.depends_on(graph.find(node), symbol_x) || graph.depends_on(graph.find(node), symbol_y)) {
                return None;
            }
        }
        polys.push(p);
    }
    let [f, g, ..] = polys.as_slice() else {
        return None;
    };
    let r = super::elim::resultant(f, g, gy)?;
    if r.is_zero() || r.degree_in(gx) == 0 {
        return None;
    }
    let mut out = Vec::new();
    for x0 in super::polynomial_roots(graph, &gens, &r, gx)? {
        let first = graph.substitute(exprs[0], x, x0);
        let first = best(graph, first)?;
        for y0 in solve_for(graph, first, y, 0)? {
            out.push(vec![x0, y0]);
        }
    }
    Some(out)
}

fn size(
    graph: &Graph,
    node: NodeId,
) -> usize {
    let mut seen = std::collections::HashSet::new();
    let mut stack = vec![node];
    while let Some(n) = stack.pop() {
        if seen.insert(n) {
            stack.extend_from_slice(graph.children(n));
        }
    }
    seen.len()
}

/// `expr = c0 + c1 u` with `c0`, `c1` free of `u`: the value of `u`.
fn linear_value(
    graph: &mut Graph,
    expr: NodeId,
    u: NodeId,
) -> Option<NodeId> {
    let symbol = graph.symbol_of(u)?;
    let mut gens = Gens::default();
    let gu = gens.index(graph, u);
    let limits = Limits::default();
    let poly = from_term(graph, &mut gens, expr, limits)?;
    if poly.degree_in(gu) != 1 {
        return None;
    }
    for g in poly.support() {
        if g != gu && gens.node(g).is_some_and(|n| graph.depends_on(graph.find(n), symbol)) {
            return None;
        }
    }
    let parts = poly.coefficients_in(gu);
    let (c0, c1) = (to_term(graph, &gens, parts.first()?), to_term(graph, &gens, parts.get(1)?));
    let minus_one = graph.int(-1);
    let inverse = super::reciprocal(graph, c1);
    Some(super::product(graph, &[minus_one, c0, inverse]))
}

/// Elimination for systems that are not polynomial.
pub(super) fn eliminate(
    graph: &mut Graph,
    exprs: &[NodeId],
    unknowns: &[NodeId],
    depth: usize,
) -> Option<Vec<Vec<NodeId>>> {
    if depth > 6 {
        return None;
    }
    if exprs.is_empty() {
        // Nothing constrains the rest: they are free.
        return Some(vec![unknowns.to_vec()]);
    }
    let symbols = unknown_symbols(graph, unknowns)?;
    if unknowns.is_empty() {
        let ok = exprs.iter().all(|&e| {
            let v = graph.eval(e, &Env::numeric(0.0));
            v.is_none_or(|v| v.abs() < 1e-9)
        });
        return Some(if ok { vec![Vec::new()] } else { Vec::new() });
    }
    // (linear?, size, equation, unknown)
    let mut options: Vec<(bool, usize, usize, usize, Option<NodeId>)> = Vec::new();
    for (i, &e) in exprs.iter().enumerate() {
        let involved: Vec<usize> = (0..unknowns.len()).filter(|&j| graph.depends_on(graph.find(e), symbols[j])).collect();
        if involved.is_empty() {
            // A condition without unknowns.
            if let Some(v) = graph.eval(e, &Env::numeric(0.0)) {
                if v.abs() > 1e-9 {
                    return Some(Vec::new());
                }
            }
            continue;
        }
        for &j in &involved {
            let linear = linear_value(graph, e, unknowns[j]);
            options.push((linear.is_none(), size(graph, e), i, j, linear));
        }
    }
    options.sort_by_key(|o| (o.0, o.1, o.2, o.3));
    for (_, _, i, j, linear) in options {
        let (e, u) = (exprs[i], unknowns[j]);
        let values = match linear {
            | Some(v) => vec![v],
            | None => match solve_for(graph, e, u, 0) {
                | Some(values) => values,
                | None => continue,
            },
        };
        let rest_exprs: Vec<NodeId> = exprs.iter().enumerate().filter(|&(k, _)| k != i).map(|(_, &x)| x).collect();
        let rest_unknowns: Vec<NodeId> = unknowns.iter().enumerate().filter(|&(k, _)| k != j).map(|(_, &x)| x).collect();
        let mut out = Vec::new();
        let mut complete = true;
        for value in values {
            let substituted: Vec<NodeId> = rest_exprs.iter().map(|&r| graph.substitute(r, u, value)).collect();
            let substituted: Vec<NodeId> = substituted.iter().map(|&s| best(graph, s).unwrap_or(s)).collect();
            let Some(tails) = eliminate(graph, &substituted, &rest_unknowns, depth + 1) else {
                complete = false;
                break;
            };
            for tail in tails {
                let mut v = value;
                for (&other, &t) in rest_unknowns.iter().zip(&tail) {
                    v = graph.substitute(v, other, t);
                }
                let mut full = tail.clone();
                full.insert(j, v);
                out.push(full);
            }
        }
        if complete {
            return Some(out);
        }
    }
    None
}
