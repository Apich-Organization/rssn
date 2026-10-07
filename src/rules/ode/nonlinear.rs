//! Systems of ordinary differential equations beyond the constant-matrix
//! linear case.
//!
//! The system is first brought to first-order form `X' = F(t, X)`: every
//! derivative below the highest order of an unknown becomes a state
//! variable, and the highest derivatives are solved for (a linear system in
//! them, or one equation each). Then, in order:
//!
//! 1. **Constant linear systems** (including those that came from higher
//!    orders) by the elimination of [`super::systems`];
//! 2. **Triangular and decoupled systems**: a state whose equation involves
//!    only itself, `t` and states already solved is solved by the scalar
//!    solver and substituted into the rest;
//! 3. **Linear systems `X' = a(t) M X`** with a constant matrix `M`: the
//!    constant-coefficient system in `s = ∫ a dt`, then `s` replaced;
//! 4. **Two autonomous equations**: the time is eliminated, `dy/dx = g/f`
//!    (or `dx/dy = f/g`), and the first-order solver finds a first integral
//!    (Hamiltonian systems give the energy, Lotka–Volterra systems the
//!    implicit integral `d x - c ln x + b y - a ln y`);
//! 5. **Polynomial autonomous systems**: polynomial first integrals up to
//!    degree three by an ansatz with undetermined coefficients, each
//!    determined by an identity in the states.
//!
//! The answer is `list(x1(t) = …, …)` for an explicit solution and
//! `list(Φ(x1(t), …) = 0, …)` for first integrals (with the constants
//! `C1`, `C2`, … inside).

use super::add;
use super::mul;
use super::solve_equation;
use super::sub;
use super::systems;
use crate::graph::Cx;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::op::core;
use crate::rules::calculus::derivative;
use crate::rules::solve::as_expression;
use crate::rules::solve::solve_for;
use crate::rules::solve::solve_linear;

/// The first-order form of a system.
struct FirstOrder {
    /// The independent variable.
    t: NodeId,
    /// Original unknown functions.
    funcs: Vec<NodeId>,
    /// State functions `w_k(t)` (fresh undetermined functions) in order:
    /// for each unknown its derivatives of order `0..m`.
    states: Vec<NodeId>,
    /// State stand-in symbols, parallel to `states`.
    symbols: Vec<NodeId>,
    /// `X' = F`, parallel to `states`, in terms of the stand-ins and `t`.
    rates: Vec<NodeId>,
    /// For each original unknown the index of its order-0 state.
    primary: Vec<usize>,
}

fn occurs(
    cx: &Cx<'_>,
    term: NodeId,
    needle: NodeId,
) -> bool {
    let mut stack = vec![term];
    let mut seen = std::collections::HashSet::new();
    while let Some(n) = stack.pop() {
        if n == needle {
            return true;
        }
        if seen.insert(n) {
            stack.extend_from_slice(cx.graph.children(n));
        }
    }
    false
}

fn first_order_form(
    cx: &mut Cx<'_>,
    equations: &[NodeId],
    funcs: &[NodeId],
) -> Option<FirstOrder> {
    let diff = cx.graph.ops().lookup("diff")?;
    let t = *cx.graph.children(*funcs.first()?).get(1)?;
    let mut residuals: Vec<NodeId> = equations.iter().map(|&e| as_expression(cx.graph, e)).collect();
    if residuals.len() != funcs.len() {
        return None;
    }
    // The chains f, f', f'', … that occur.
    let mut orders: Vec<Vec<NodeId>> = Vec::new();
    for &f in funcs {
        let mut chain = vec![f];
        loop {
            let next = cx.graph.node(diff, &[*chain.last()?, t]);
            if chain.len() > 6 || !residuals.iter().any(|&e| occurs(cx, e, next)) {
                break;
            }
            chain.push(next);
        }
        if chain.len() < 2 {
            return None;
        }
        orders.push(chain);
    }
    let mut symbols = Vec::new();
    let mut state_funcs = Vec::new();
    let mut primary = Vec::new();
    let mut highest = Vec::new();
    let mut offsets = Vec::new();
    for chain in &orders {
        offsets.push(symbols.len());
        primary.push(symbols.len());
        for _ in 0..chain.len() - 1 {
            let s = cx.graph.interner_mut().fresh_symbol("s");
            symbols.push(cx.graph.symbol_node(s));
            state_funcs.push(super::reduce::fresh_function(cx.graph, "w", t));
        }
        let h = cx.graph.interner_mut().fresh_symbol("d");
        highest.push(cx.graph.symbol_node(h));
    }
    // Highest derivatives first, so that f does not hide f'.
    for (i, chain) in orders.iter().enumerate() {
        let m = chain.len();
        for k in (0..m).rev() {
            let replacement = if k + 1 == m { highest[i] } else { symbols[offsets[i] + k] };
            for r in &mut residuals {
                *r = cx.graph.replace_subterm(*r, chain[k], replacement);
            }
        }
    }
    let solved = match solve_linear(cx.graph, &residuals, &highest) {
        | Some(v) => v,
        | None => {
            let mut out = Vec::new();
            for (r, &h) in residuals.iter().zip(&highest) {
                out.push(*solve_for(cx.graph, *r, h, 0)?.first()?);
            }
            out
        },
    };
    let mut rates = vec![NodeId::NONE; symbols.len()];
    for (i, chain) in orders.iter().enumerate() {
        let order = chain.len() - 1;
        for k in 0..order {
            rates[offsets[i] + k] = if k + 1 < order { symbols[offsets[i] + k + 1] } else { solved[i] };
        }
    }
    Some(FirstOrder { t, funcs: funcs.to_vec(), states: state_funcs, symbols, rates, primary })
}

/// `rates` with the stand-in symbols replaced by the state functions.
fn in_functions(
    cx: &mut Cx<'_>,
    system: &FirstOrder,
    e: NodeId,
) -> NodeId {
    let mut out = e;
    for (&s, &w) in system.symbols.iter().zip(&system.states) {
        out = cx.graph.substitute(out, s, w);
    }
    out
}

/// Renames `C1, C2, …` in `e` to `C{offset+1}, …`.
fn shift_constants(
    cx: &mut Cx<'_>,
    e: NodeId,
    offset: usize,
) -> NodeId {
    if offset == 0 {
        return e;
    }
    let count = super::reduce::constants_in(cx.graph, e);
    let mut out = e;
    for k in (1..=count).rev() {
        let (from, to) = (cx.graph.sym(&format!("C{k}")), cx.graph.sym(&format!("C{}", k + offset)));
        out = cx.graph.substitute(out, from, to);
    }
    out
}

/// The answer for the original unknowns from state solutions.
fn assemble(
    cx: &mut Cx<'_>,
    system: &FirstOrder,
    values: &[NodeId],
) -> NodeId {
    let mut items = Vec::new();
    for (&f, &index) in system.funcs.iter().zip(&system.primary) {
        let v = cx.simplify(values[index]);
        items.push(cx.graph.node(core::EQ, &[f, v]));
    }
    cx.graph.node(core::LIST, &items)
}

fn depends_on_symbols(
    cx: &Cx<'_>,
    e: NodeId,
    symbols: &[NodeId],
) -> Vec<usize> {
    symbols
        .iter()
        .enumerate()
        .filter(|&(_, &s)| cx.graph.symbol_of(s).is_some_and(|sym| cx.graph.depends_on(cx.graph.find(e), sym)))
        .map(|(i, _)| i)
        .collect()
}

/// Triangular and decoupled systems.
fn triangular(
    cx: &mut Cx<'_>,
    system: &FirstOrder,
    depth: u32,
) -> Option<NodeId> {
    let n = system.symbols.len();
    let mut solved: Vec<Option<NodeId>> = vec![None; n];
    let mut constants = 0;
    for _ in 0..n {
        let pick = (0..n).find(|&i| {
            solved[i].is_none()
                && depends_on_symbols(cx, system.rates[i], &system.symbols).iter().all(|&j| j == i || solved[j].is_some())
        })?;
        let mut rhs = system.rates[pick];
        for (j, value) in solved.iter().enumerate() {
            if let Some(v) = *value {
                rhs = cx.graph.substitute(rhs, system.symbols[j], v);
            }
        }
        let own = system.symbols[pick];
        let w = system.states[pick];
        let rhs = cx.graph.substitute(rhs, own, w);
        // Constants already present must not collide with the new ones.
        let mut parked = Vec::new();
        let mut rhs = rhs;
        for k in 1..=constants {
            let from = cx.graph.sym(&format!("C{k}"));
            let to_symbol = cx.graph.interner_mut().fresh_symbol("K");
            let to = cx.graph.symbol_node(to_symbol);
            rhs = cx.graph.substitute(rhs, from, to);
            parked.push((to, from));
        }
        let diff = cx.graph.ops().lookup("diff")?;
        let lhs = cx.graph.node(diff, &[w, system.t]);
        let equation = cx.graph.node(core::EQ, &[lhs, rhs]);
        let (_, answer) = solve_equation(cx, equation, w, depth + 1)?;
        let &[left, value] = cx.graph.children(answer) else {
            return None;
        };
        if left != w {
            return None;
        }
        let mut value = shift_constants(cx, value, constants);
        for (to, from) in parked {
            value = cx.graph.substitute(value, to, from);
        }
        // The solved states may carry constants of their own.
        constants = constants.max(super::reduce::constants_in(cx.graph, value));
        solved[pick] = Some(value);
    }
    let values: Vec<NodeId> = solved.into_iter().collect::<Option<_>>()?;
    Some(assemble(cx, system, &values))
}

/// `X' = a(t) M X` through the constant system in `s = ∫ a dt`.
fn time_scaled_linear(
    cx: &mut Cx<'_>,
    system: &FirstOrder,
    depth: u32,
) -> Option<NodeId> {
    let n = system.symbols.len();
    let t = system.t;
    let t_symbol = cx.graph.symbol_of(t)?;
    // A_ij = ∂F_i/∂s_j must be free of the states.
    let mut a: Vec<Vec<NodeId>> = Vec::with_capacity(n);
    let mut reference: Option<NodeId> = None;
    for i in 0..n {
        let mut row = Vec::with_capacity(n);
        for j in 0..n {
            let d = derivative(cx.graph, system.rates[i], system.symbols[j])?;
            let d = cx.simplify(d);
            if !depends_on_symbols(cx, d, &system.symbols).is_empty() {
                return None;
            }
            if reference.is_none() && !cx.is_zero(d) {
                reference = Some(d);
            }
            row.push(d);
        }
        a.push(row);
        // Homogeneous: F_i must be exactly Σ A_ij s_j.
        let mut zeroed = system.rates[i];
        for &s in &system.symbols {
            let zero = cx.graph.int(0);
            zeroed = cx.graph.substitute(zeroed, s, zero);
        }
        let zeroed = cx.simplify(zeroed);
        if !cx.is_zero(zeroed) {
            return None;
        }
    }
    let scale = reference?;
    if !cx.graph.depends_on(cx.graph.find(scale), t_symbol) {
        return None;
    }
    let s_symbol = cx.graph.interner_mut().fresh_symbol("s");
    let s = cx.graph.symbol_node(s_symbol);
    let diff = cx.graph.ops().lookup("diff")?;
    let mut zs = Vec::new();
    for _ in 0..n {
        zs.push(super::reduce::fresh_function(cx.graph, "z", s));
    }
    let mut equations = Vec::new();
    for i in 0..n {
        let mut terms = Vec::new();
        for j in 0..n {
            let ratio = super::div(cx.graph, a[i][j], scale);
            let ratio = cx.simplify(ratio);
            if cx.graph.depends_on(cx.graph.find(ratio), t_symbol) {
                return None;
            }
            terms.push(mul(cx.graph, &[ratio, zs[j]]));
        }
        let rhs = add(cx.graph, &terms);
        let lhs = cx.graph.node(diff, &[zs[i], s]);
        equations.push(cx.graph.node(core::EQ, &[lhs, rhs]));
    }
    let eq_list = cx.graph.node(core::LIST, &equations);
    let fn_list = cx.graph.node(core::LIST, &zs);
    let _ = depth;
    let answer = systems::solve_system(cx, eq_list, fn_list)?;
    let integral = super::integrate(cx, scale, t)?;
    let items = cx.graph.children(answer).to_vec();
    let mut values = Vec::new();
    for (k, item) in items.iter().enumerate() {
        let &[_, v] = cx.graph.children(*item) else {
            return None;
        };
        let in_t = cx.graph.substitute(v, s, integral);
        let _ = k;
        values.push(in_t);
    }
    if values.len() != n {
        return None;
    }
    Some(assemble(cx, system, &values))
}

/// Two autonomous equations: a first integral from `dy/dx = g/f`.
fn autonomous_integral(
    cx: &mut Cx<'_>,
    system: &FirstOrder,
    depth: u32,
) -> Option<NodeId> {
    if system.symbols.len() != 2 || system.primary != [0, 1] {
        return None;
    }
    let t_symbol = cx.graph.symbol_of(system.t)?;
    let (f, g) = (system.rates[0], system.rates[1]);
    if cx.graph.depends_on(cx.graph.find(f), t_symbol) || cx.graph.depends_on(cx.graph.find(g), t_symbol) {
        return None;
    }
    let diff = cx.graph.ops().lookup("diff")?;
    for swap in [false, true] {
        let (a, b) = if swap { (1, 0) } else { (0, 1) };
        // Independent variable symbol X and function Y(X).
        let x_symbol = cx.graph.interner_mut().fresh_symbol("X");
        let x = cx.graph.symbol_node(x_symbol);
        let y_fn = super::reduce::fresh_function(cx.graph, "Y", x);
        let (num, den) = (system.rates[b], system.rates[a]);
        if cx.is_zero(den) {
            continue;
        }
        let mut num = cx.graph.substitute(num, system.symbols[a], x);
        let mut den = cx.graph.substitute(den, system.symbols[a], x);
        num = cx.graph.substitute(num, system.symbols[b], y_fn);
        den = cx.graph.substitute(den, system.symbols[b], y_fn);
        let slope = super::div(cx.graph, num, den);
        let slope = cx.simplify(slope);
        let lhs = cx.graph.node(diff, &[y_fn, x]);
        let equation = cx.graph.node(core::EQ, &[lhs, slope]);
        let Some((_, answer)) = solve_equation(cx, equation, y_fn, depth + 1) else {
            continue;
        };
        let &[left, right] = cx.graph.children(answer) else {
            continue;
        };
        // Back to the original functions.
        let (xa, yb) = (system.funcs[a], system.funcs[b]);
        let relation = if left == y_fn {
            // A relation between the two unknowns, written as `Φ = 0`.
            let v = cx.graph.substitute(right, x, xa);
            let difference = sub(cx.graph, yb, v);
            let zero = cx.graph.int(0);
            cx.graph.node(core::EQ, &[difference, zero])
        } else {
            let l = cx.graph.replace_subterm(left, y_fn, yb);
            let l = cx.graph.substitute(l, x, xa);
            let r = cx.graph.replace_subterm(right, y_fn, yb);
            let r = cx.graph.substitute(r, x, xa);
            cx.graph.node(core::EQ, &[l, r])
        };
        return Some(cx.graph.node(core::LIST, &[relation]));
    }
    None
}

/// Polynomial first integrals of an autonomous polynomial system.
fn polynomial_integrals(
    cx: &mut Cx<'_>,
    system: &FirstOrder,
) -> Option<NodeId> {
    let n = system.symbols.len();
    let t_symbol = cx.graph.symbol_of(system.t)?;
    if !(2..=4).contains(&n) {
        return None;
    }
    for &r in &system.rates {
        if cx.graph.depends_on(cx.graph.find(r), t_symbol) {
            return None;
        }
    }
    for degree in 1..=3_u32 {
        // Monomials of total degree <= degree (and at least 1).
        let mut monomials: Vec<Vec<u32>> = vec![vec![0; n]];
        for _ in 0..degree {
            let mut next = Vec::new();
            for m in &monomials {
                for i in 0..n {
                    let mut e = m.clone();
                    e[i] += 1;
                    if !next.contains(&e) {
                        next.push(e);
                    }
                }
            }
            for m in next {
                if !monomials.contains(&m) {
                    monomials.push(m);
                }
            }
        }
        monomials.retain(|m| m.iter().sum::<u32>() >= 1);
        let mut unknowns = Vec::new();
        let mut terms = Vec::new();
        let mut candidate_terms = Vec::new();
        for m in &monomials {
            let symbol = cx.graph.interner_mut().fresh_symbol("c");
            let c = cx.graph.symbol_node(symbol);
            unknowns.push(c);
            let mut factors = vec![c];
            for (i, &e) in m.iter().enumerate() {
                if e > 0 {
                    let power = cx.graph.int(i64::from(e));
                    factors.push(cx.graph.node(core::POW, &[system.symbols[i], power]));
                }
            }
            candidate_terms.push(mul(cx.graph, &factors));
        }
        let candidate = add(cx.graph, &candidate_terms);
        for i in 0..n {
            let d = derivative(cx.graph, candidate, system.symbols[i])?;
            terms.push(mul(cx.graph, &[d, system.rates[i]]));
        }
        let total = add(cx.graph, &terms);
        let basis = super::lie::solve_identity(cx.graph, total, &unknowns)?;
        if basis.is_empty() {
            continue;
        }
        let mut relations = Vec::new();
        for (k, v) in basis.iter().enumerate() {
            let mut h = candidate;
            for (u, value) in unknowns.iter().zip(v) {
                let value = cx.graph.num(Number::rat(value.clone()));
                h = cx.graph.substitute(h, *u, value);
            }
            let h = cx.simplify(h);
            if cx.is_zero(h) || relations.len() >= 3 {
                continue;
            }
            // In terms of the original unknowns.
            let mut in_funcs = h;
            for (&index, &f) in system.primary.iter().zip(&system.funcs) {
                in_funcs = cx.graph.substitute(in_funcs, system.symbols[index], f);
            }
            let c = cx.graph.sym(&format!("C{}", relations.len() + 1));
            let relation = sub(cx.graph, in_funcs, c);
            let zero = cx.graph.int(0);
            relations.push(cx.graph.node(core::EQ, &[relation, zero]));
            let _ = k;
        }
        if !relations.is_empty() {
            return Some(cx.graph.node(core::LIST, &relations));
        }
    }
    None
}

/// Entry point: systems that the constant-coefficient solver cannot do.
pub(super) fn solve(
    cx: &mut Cx<'_>,
    equations: NodeId,
    unknowns: NodeId,
    depth: u32,
) -> Option<NodeId> {
    if cx.graph.op(equations) != core::LIST || cx.graph.op(unknowns) != core::LIST {
        return None;
    }
    let eqs = cx.graph.children(equations).to_vec();
    let funcs = cx.graph.children(unknowns).to_vec();
    let system = first_order_form(cx, &eqs, &funcs)?;
    let n = system.symbols.len();
    // Constant linear systems, higher orders included.
    if let Some(found) = reduced_linear(cx, &system) {
        return Some(found);
    }
    if let Some(found) = triangular(cx, &system, depth) {
        return Some(found);
    }
    if let Some(found) = time_scaled_linear(cx, &system, depth) {
        return Some(found);
    }
    if n == 2
        && let Some(found) = autonomous_integral(cx, &system, depth) {
            return Some(found);
        }
    if system.primary.len() == n {
        return polynomial_integrals(cx, &system);
    }
    None
}

/// The first-order system as equations in the state functions, handed to
/// the constant-coefficient solver.
fn reduced_linear(
    cx: &mut Cx<'_>,
    system: &FirstOrder,
) -> Option<NodeId> {
    let diff = cx.graph.ops().lookup("diff")?;
    let mut equations = Vec::new();
    for (i, &w) in system.states.iter().enumerate() {
        let rhs = in_functions(cx, system, system.rates[i]);
        let lhs = cx.graph.node(diff, &[w, system.t]);
        equations.push(cx.graph.node(core::EQ, &[lhs, rhs]));
    }
    let eq_list = cx.graph.node(core::LIST, &equations);
    let fn_list = cx.graph.node(core::LIST, &system.states);
    let answer = systems::solve_system(cx, eq_list, fn_list)?;
    let items = cx.graph.children(answer).to_vec();
    let mut values = Vec::new();
    for item in items {
        let &[_, v] = cx.graph.children(item) else {
            return None;
        };
        values.push(v);
    }
    if values.len() != system.states.len() {
        return None;
    }
    Some(assemble(cx, system, &values))
}
