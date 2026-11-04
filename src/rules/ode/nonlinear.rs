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
//!
//! **Darboux polynomials and the Jacobi last multiplier** are separate report
//! operators for a planar autonomous polynomial system `x' = f(x, y)`,
//! `y' = g(x, y)`. A Darboux polynomial is a polynomial `P` with
//! `X(P) = k P` for a polynomial cofactor `k` of degree at most `deg X - 1`;
//! `P = 0` is then an invariant curve. The search fixes the leading
//! homogeneous form first: the linear factors of `x g_m - y f_m` (top
//! homogeneous parts of the rates) are the invariant lines of the leading
//! vector field, so the possible leading forms of a degree-`d` Darboux
//! polynomial are products of the irreducible factors of that form with
//! exponents adding up to `d` (factored over `Q`). The lower coefficients and
//! the cofactor then come from a linear descent by total degree: the
//! degree-`D` coefficient of `X(P) = k P` involves only the new homogeneous
//! part of `P`, the new part of `k` and already known quantities, never a
//! product of two unknowns. Each degree is an exact linear system over `Q`.
//! A Jacobi last multiplier is then `M = Π P_i^{λ_i}` with
//! `Σ λ_i k_i = -div X` solved over the rationals, and
//! `X(M) + (div X) M = 0` is checked numerically before it is returned.
//!
//! Limitations of the Darboux machinery: the search stops at the requested
//! total degree (four by default, six at most) and the cofactor degree
//! `deg X - 1`; a radial leading part (`x g_m = y f_m`) does not restrict the
//! leading forms and the search gives up; systems with symbolic parameters or
//! non-polynomial rates are out of scope; the last multiplier exists only when
//! the found Darboux polynomials span `-div X` over the rationals.

use num_traits::One;
use num_traits::Signed;
use num_traits::Zero;

use super::add;
use super::mul;
use super::solve_equation;
use super::sub;
use super::systems;
use crate::graph::Cx;
use crate::graph::Env;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::op::core;
use crate::rules::calculus::derivative;
use crate::rules::poly::repr::Gens;
use crate::rules::poly::repr::Limits;
use crate::rules::poly::repr::from_term;
use crate::rules::poly::univariate;
use crate::rules::poly::univariate::QPoly;
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
    let mut residuals: Vec<NodeId> = equations
        .iter()
        .map(|&e| as_expression(cx.graph, e))
        .collect();
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
            let replacement = if k + 1 == m {
                highest[i]
            } else {
                symbols[offsets[i] + k]
            };
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
            rates[offsets[i] + k] = if k + 1 < order {
                symbols[offsets[i] + k + 1]
            } else {
                solved[i]
            };
        }
    }
    Some(FirstOrder {
        t,
        funcs: funcs.to_vec(),
        states: state_funcs,
        symbols,
        rates,
        primary,
    })
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
        let (from, to) = (
            cx.graph.sym(&format!("C{k}")),
            cx.graph.sym(&format!("C{}", k + offset)),
        );
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
        .filter(|&(_, &s)| {
            cx.graph
                .symbol_of(s)
                .is_some_and(|sym| cx.graph.depends_on(cx.graph.find(e), sym))
        })
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
                && depends_on_symbols(cx, system.rates[i], &system.symbols)
                    .iter()
                    .all(|&j| j == i || solved[j].is_some())
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
    if cx.graph.depends_on(cx.graph.find(f), t_symbol)
        || cx.graph.depends_on(cx.graph.find(g), t_symbol)
    {
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
        && let Some(found) = autonomous_integral(cx, &system, depth)
    {
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
// ---------------------------------------------------------------------
// Darboux polynomials and the Jacobi last multiplier
// ---------------------------------------------------------------------
//
// For a planar autonomous polynomial system
//
//     x' = f(x, y),   y' = g(x, y),   X = f ∂x + g ∂y,
//
// a *Darboux polynomial* is a polynomial `P` with `X(P) = k P` for a
// polynomial *cofactor* `k` of degree at most `deg X - 1`; the curve
// `P = 0` is then invariant under the flow. A *Jacobi last multiplier* is a
// function `M` with `d ln M/dt = -div X`; when the Darboux polynomials and
// their cofactors are known, `M = Π P_i^{λ_i}` with `Σ λ_i k_i = -div X`
// solved over the rationals.
//
// The search has two steps. First the leading homogeneous form. The linear
// factors of `R(x, y) = x g_m - y f_m` (where `f_m, g_m` are the top
// homogeneous parts of degree `m = deg X`) are exactly the invariant lines
// of the leading vector field: for a line `ℓ = 0`, `X_m(ℓ)` is a multiple of
// `ℓ` precisely when `ℓ | R`. Over the rationals the possible leading forms
// of a Darboux polynomial of degree `d` are therefore products of the
// irreducible factors of `R` with exponents adding up to `d`. Second, with
// the leading form fixed, the remaining coefficients of `P` and of `k` are
// found by a *linear* descent by total degree. The coefficient of total
// degree `D` of `X(P) = k P` involves only the leading parts of the rates
// acting on the single unknown homogeneous part `P_{D-deg k}`, the known
// leading cofactor part acting on that same part, the one unknown cofactor
// part `k_{D-deg P}` acting on the known leading form, and quantities
// already determined — never a product of two unknowns. Each degree is thus
// an independent linear system over `Q`, solved exactly here.

/// Largest total degree of the first Darboux polynomial search, the default
/// of the operators.
const DARBOUX_DEGREE: u32 = 4;
/// Largest total degree an explicit argument of the operator may ask for.
const MAX_DARBOUX_DEGREE: u32 = 6;
/// Largest number of leading forms tried for one degree; beyond it the
/// enumeration of factor exponents is cut short (and the result may miss a
/// Darboux polynomial, never invent one).
const MAX_DARBOUX_CANDIDATES: usize = 96;
/// Points at which a constructed Jacobi last multiplier is checked.
const MULTIPLIER_CHECKS: [(f64, f64); 5] =
    [(1.3, 0.7), (2.1, 1.7), (0.4, 2.3), (1.7, 0.3), (3.1, 2.9)];

/// A polynomial in the two state variables over the rationals.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct Planar {
    /// `(exponent of x, exponent of y)` to a non-zero coefficient.
    terms: std::collections::BTreeMap<(u32, u32), num_rational::BigRational>,
}

type Rat = num_rational::BigRational;

fn rat(value: i64) -> Rat {
    Rat::from_integer(num_bigint::BigInt::from(value))
}

fn one() -> Rat {
    rat(1)
}

impl Planar {
    fn add_term(
        &mut self,
        exponent: (u32, u32),
        coefficient: Rat,
    ) {
        if coefficient.is_zero() {
            return;
        }
        match self.terms.get_mut(&exponent) {
            | Some(slot) => {
                *slot += coefficient;
                if slot.is_zero() {
                    self.terms.remove(&exponent);
                }
            },
            | None => {
                self.terms.insert(exponent, coefficient);
            },
        }
    }

    fn monomial(
        exponent: (u32, u32),
        coefficient: Rat,
    ) -> Self {
        let mut out = Self::default();
        out.add_term(exponent, coefficient);
        out
    }

    fn variable_x() -> Self {
        Self::monomial((1, 0), one())
    }

    fn variable_y() -> Self {
        Self::monomial((0, 1), one())
    }

    fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    fn degree(&self) -> Option<u32> {
        self.terms.keys().map(|&(x, y)| x.saturating_add(y)).max()
    }

    fn homogeneous_part(
        &self,
        degree: u32,
    ) -> Self {
        let mut out = Self::default();
        for (&(x, y), coefficient) in &self.terms {
            if x.saturating_add(y) == degree {
                out.add_term((x, y), coefficient.clone());
            }
        }
        out
    }

    fn add(
        &self,
        other: &Self,
    ) -> Self {
        let mut out = self.clone();
        for (&exponent, coefficient) in &other.terms {
            out.add_term(exponent, coefficient.clone());
        }
        out
    }

    fn sub(
        &self,
        other: &Self,
    ) -> Self {
        let mut out = self.clone();
        for (&exponent, coefficient) in &other.terms {
            out.add_term(exponent, -coefficient.clone());
        }
        out
    }

    fn scale(
        &self,
        factor: &Rat,
    ) -> Self {
        if factor.is_zero() {
            return Self::default();
        }
        let mut out = Self::default();
        for (&exponent, coefficient) in &self.terms {
            out.add_term(exponent, coefficient * factor);
        }
        out
    }

    fn mul(
        &self,
        other: &Self,
    ) -> Self {
        let mut out = Self::default();
        for (&(ax, ay), a) in &self.terms {
            for (&(bx, by), b) in &other.terms {
                out.add_term((ax + bx, ay + by), a * b);
            }
        }
        out
    }

    fn pow(
        &self,
        exponent: u32,
    ) -> Self {
        let mut out = Self::monomial((0, 0), one());
        for _ in 0..exponent {
            out = out.mul(self);
        }
        out
    }

    fn derivative_x(&self) -> Self {
        let mut out = Self::default();
        for (&(x, y), coefficient) in &self.terms {
            if x > 0 {
                out.add_term((x - 1, y), coefficient * rat(i64::from(x)));
            }
        }
        out
    }

    fn derivative_y(&self) -> Self {
        let mut out = Self::default();
        for (&(x, y), coefficient) in &self.terms {
            if y > 0 {
                out.add_term((x, y - 1), coefficient * rat(i64::from(y)));
            }
        }
        out
    }

    fn coefficient(
        &self,
        monomial: &Self,
    ) -> Rat {
        let Some((exponent, _)) = monomial.terms.iter().next() else {
            return Rat::zero();
        };
        self.terms.get(exponent).cloned().unwrap_or_else(Rat::zero)
    }

    /// The leading term in the graded lexicographic order, used by the
    /// bivariate division.
    fn leading(&self) -> Option<((u32, u32), Rat)> {
        let total = |m: (u32, u32)| m.0.saturating_add(m.1);
        self.terms
            .iter()
            .max_by(|a, b| {
                total(*a.0)
                    .cmp(&total(*b.0))
                    .then_with(|| a.0.0.cmp(&b.0.0))
                    .then_with(|| a.0.1.cmp(&b.0.1))
            })
            .map(|(&exponent, coefficient)| (exponent, coefficient.clone()))
    }

    /// The exact quotient `self / divisor`, or `None` when the division
    /// leaves a remainder.
    fn div_exact(
        &self,
        divisor: &Self,
    ) -> Option<Self> {
        let ((dx, dy), dc) = divisor.leading()?;
        let mut quotient = Self::default();
        let mut remainder = self.clone();
        let mut guard = 0;
        while !remainder.is_zero() {
            guard += 1;
            if guard > 20_000 {
                return None;
            }
            let ((rx, ry), rc) = remainder.leading()?;
            if rx < dx || ry < dy {
                return None;
            }
            let exponent = (rx - dx, ry - dy);
            let coefficient = rc / &dc;
            quotient.add_term(exponent, coefficient.clone());
            let term = Self::monomial(exponent, coefficient);
            remainder = remainder.sub(&divisor.mul(&term));
        }
        Some(quotient)
    }

    /// `self` scaled to primitive integer coefficients with a positive
    /// leading coefficient, and the rational scale applied.
    fn primitive_with_scale(&self) -> (Self, Rat) {
        if self.is_zero() {
            return (Self::default(), one());
        }
        let mut denominator = num_bigint::BigInt::from(1);
        for coefficient in self.terms.values() {
            denominator = num_integer::Integer::lcm(&denominator, coefficient.denom());
        }
        let mut integer: std::collections::BTreeMap<(u32, u32), num_bigint::BigInt> =
            std::collections::BTreeMap::new();
        for (&exponent, coefficient) in &self.terms {
            let scaled = coefficient * Rat::from_integer(denominator.clone());
            integer.insert(exponent, scaled.to_integer());
        }
        let mut content = num_bigint::BigInt::from(0);
        for value in integer.values() {
            content = num_integer::Integer::gcd(&content, value);
        }
        let mut scale = Rat::new(denominator, content.max(num_bigint::BigInt::from(1)));
        // The leading coefficient fixes the overall sign.
        let leading = self.leading().map_or_else(Rat::zero, |(_, c)| c);
        if leading.is_negative() {
            scale = -scale;
        }
        (self.scale(&scale), scale)
    }

    fn to_node(
        &self,
        graph: &mut Graph,
        x: NodeId,
        y: NodeId,
    ) -> NodeId {
        let mut terms = Vec::with_capacity(self.terms.len());
        for (&(ex, ey), coefficient) in &self.terms {
            let mut factors = Vec::new();
            if !coefficient.is_one() {
                factors.push(graph.num(Number::rat(coefficient.clone())));
            }
            if ex > 0 {
                factors.push(power_of(graph, x, ex));
            }
            if ey > 0 {
                factors.push(power_of(graph, y, ey));
            }
            terms.push(match factors.as_slice() {
                | [] => graph.int(1),
                | [only] => *only,
                | _ => graph.node(core::MUL, &factors),
            });
        }
        match terms.as_slice() {
            | [] => graph.int(0),
            | [only] => *only,
            | _ => graph.node(core::ADD, &terms),
        }
    }

    fn from_node(
        graph: &mut Graph,
        node: NodeId,
        x: NodeId,
        y: NodeId,
    ) -> Option<Self> {
        let mut gens = Gens::default();
        let gx = gens.index(graph, x);
        let gy = gens.index(graph, y);
        let poly = from_term(
            graph,
            &mut gens,
            node,
            Limits {
                terms: 20_000,
                exponent: 64,
            },
        )?;
        let mut out = Self::default();
        for (monomial, coefficient) in poly.terms() {
            let mut exponent = (0_u32, 0_u32);
            for &(generator, power) in monomial {
                if generator == gx {
                    exponent.0 = power;
                } else if generator == gy {
                    exponent.1 = power;
                } else {
                    return None;
                }
            }
            out.add_term(exponent, coefficient.to_rational()?);
        }
        Some(out)
    }
}

/// `base^exponent` for a positive literal exponent.
fn power_of(
    graph: &mut Graph,
    base: NodeId,
    exponent: u32,
) -> NodeId {
    match exponent {
        | 0 => graph.int(1),
        | 1 => base,
        | e => {
            let e = graph.int(i64::from(e));
            graph.node(core::POW, &[base, e])
        },
    }
}

/// `X(P) = f ∂x P + g ∂y P`.
fn derive(
    f: &Planar,
    g: &Planar,
    p: &Planar,
) -> Planar {
    f.mul(&p.derivative_x()).add(&g.mul(&p.derivative_y()))
}

/// The homogeneous factors of the invariant lines of the leading vector
/// field, from the factorization over `Q` of `R = x g_m - y f_m`.
///
/// `None` when the leading part is radial (`R = 0`): then every line is
/// invariant and the leading form does not restrict the search.
fn invariant_lines(r: &Planar) -> Option<Vec<Planar>> {
    if r.is_zero() {
        return None;
    }
    let x_power = r.terms.keys().map(|&(x, _)| x).min()?;
    let mut factors = Vec::new();
    let mut rest = r.clone();
    if x_power > 0 {
        let x = Planar::variable_x();
        factors.push(x.clone());
        rest = rest.div_exact(&x.pow(x_power))?;
    }
    // `rest(1, t)`: in a homogeneous polynomial the exponent of `y` alone
    // identifies the term, so this is a faithful univariate image.
    let mut image: QPoly = Vec::new();
    for (&(_, ey), coefficient) in &rest.terms {
        let index = usize::try_from(ey).ok()?;
        if image.len() <= index {
            image.resize(index.saturating_add(1), Rat::zero());
        }
        if let Some(slot) = image.get_mut(index) {
            *slot = coefficient.clone();
        }
    }
    let (_, irreducible) = univariate::factor(&image);
    for (factor, _multiplicity) in irreducible {
        let degree = u32::try_from(factor.len()).ok()?.saturating_sub(1);
        let mut homogeneous = Planar::default();
        for (index, coefficient) in factor.iter().enumerate() {
            let index = u32::try_from(index).ok()?;
            homogeneous.add_term(
                (degree.saturating_sub(index), index),
                Rat::from_integer(coefficient.clone()),
            );
        }
        factors.push(homogeneous);
    }
    Some(factors)
}

/// The leading forms of total degree `degree`: products of powers of the
/// invariant lines, in every combination of exponents.
fn leading_forms(
    factors: &[Planar],
    degree: u32,
    found: &mut Vec<Planar>,
) {
    fn walk(
        factors: &[Planar],
        index: usize,
        remaining: u32,
        current: &Planar,
        found: &mut Vec<Planar>,
    ) {
        if found.len() >= MAX_DARBOUX_CANDIDATES {
            return;
        }
        let Some(factor) = factors.get(index) else {
            if remaining == 0 {
                found.push(current.clone());
            }
            return;
        };
        let factor_degree = factor.degree().unwrap_or(0);
        if factor_degree == 0 {
            return;
        }
        let mut power = remaining / factor_degree;
        loop {
            let next = current.mul(&factor.pow(power));
            walk(
                factors,
                index.saturating_add(1),
                remaining - power * factor_degree,
                &next,
                found,
            );
            if found.len() >= MAX_DARBOUX_CANDIDATES || power == 0 {
                return;
            }
            power -= 1;
        }
    }
    walk(factors, 0, degree, &Planar::monomial((0, 0), one()), found);
}

/// The exact solution of an augmented rational linear system, free unknowns
/// set to zero, or `None` when it is inconsistent.
fn solve_rational(
    mut rows: Vec<Vec<Rat>>,
    columns: usize,
) -> Option<Vec<Rat>> {
    let mut pivots: Vec<usize> = Vec::new();
    let mut r = 0;
    for column in 0..columns {
        let Some(pivot) = (r..rows.len()).find(|&i| !rows[i][column].is_zero()) else {
            continue;
        };
        rows.swap(r, pivot);
        let lead = rows[r][column].clone();
        for value in &mut rows[r] {
            *value /= &lead;
        }
        for i in 0..rows.len() {
            if i == r || rows[i][column].is_zero() {
                continue;
            }
            let factor = rows[i][column].clone();
            let pivot_row = rows[r].clone();
            for (value, p) in rows[i].iter_mut().zip(&pivot_row) {
                *value -= &factor * p;
            }
        }
        pivots.push(column);
        r += 1;
        if r == rows.len() {
            break;
        }
    }
    if rows.iter().any(|row| {
        row.get(..columns)
            .is_some_and(|lhs| lhs.iter().all(Rat::is_zero))
            && row.get(columns).is_some_and(|rhs| !rhs.is_zero())
    }) {
        return None;
    }
    let mut solution = vec![Rat::zero(); columns];
    for (i, &column) in pivots.iter().enumerate() {
        solution[column] = rows[i][columns].clone();
    }
    Some(solution)
}

/// The lower-degree coefficients of a Darboux polynomial whose leading form
/// is `leading`, and its cofactor, by the linear descent by total degree.
fn darboux_from_leading(
    f: &Planar,
    g: &Planar,
    f_top: &Planar,
    g_top: &Planar,
    leading: &Planar,
    cofactor_degree: u32,
) -> Option<(Planar, Planar)> {
    let degree = leading.degree()?;
    let mut parts: std::collections::BTreeMap<u32, Planar> = std::collections::BTreeMap::new();
    parts.insert(degree, leading.clone());
    let mut cofactors: std::collections::BTreeMap<u32, Planar> = std::collections::BTreeMap::new();
    // The top-degree identity `X_m(P_d) = k_δ P_d` gives the leading
    // cofactor part.
    let top_image = derive(f_top, g_top, leading);
    let k_top = top_image.div_exact(leading)?;
    cofactors.insert(cofactor_degree, k_top.clone());
    let minus_one = rat(-1);
    for level in (0..degree.saturating_add(cofactor_degree)).rev() {
        let part_index = i64::from(level).checked_sub(i64::from(cofactor_degree))?;
        let cofactor_index = i64::from(level).checked_sub(i64::from(degree))?;
        let part_known = (0..=i64::from(degree)).contains(&part_index)
            && !parts.contains_key(&u32::try_from(part_index).ok()?);
        let cofactor_known = (0..=i64::from(cofactor_degree)).contains(&cofactor_index)
            && !cofactors.contains_key(&u32::try_from(cofactor_index).ok()?);
        let part_level = u32::try_from(part_index).ok().filter(|_| part_known);
        let cofactor_level = u32::try_from(cofactor_index)
            .ok()
            .filter(|_| cofactor_known);
        // LHS and both-known RHS terms, at this total degree.
        let mut known = Planar::default();
        for (&part_degree, part) in &parts {
            if i64::from(part_degree) > part_index {
                known = known.add(&derive(f, g, part).homogeneous_part(level));
            }
        }
        for (&cofactor_power, cofactor) in &cofactors {
            if i64::from(cofactor_power) <= cofactor_index {
                continue;
            }
            for (&part_degree, part) in &parts {
                if i64::from(part_degree) > part_index
                    && cofactor_power.saturating_add(part_degree) == level
                {
                    known = known.sub(&cofactor.mul(part));
                }
            }
        }
        let mut unknowns: Vec<Planar> = Vec::new();
        let mut is_part = Vec::new();
        if let Some(index) = part_level {
            for x_power in 0..=index {
                unknowns.push(Planar::monomial((x_power, index - x_power), one()));
                is_part.push(true);
            }
        }
        if let Some(index) = cofactor_level {
            for x_power in 0..=index {
                unknowns.push(Planar::monomial((x_power, index - x_power), one()));
                is_part.push(false);
            }
        }
        if unknowns.is_empty() {
            if !known.is_zero() {
                return None;
            }
            continue;
        }
        let mut rows = Vec::with_capacity(usize::try_from(level).ok()?.saturating_add(1));
        for x_power in 0..=level {
            let monomial = Planar::monomial((x_power, level - x_power), one());
            let mut row = Vec::with_capacity(unknowns.len().saturating_add(1));
            for (unknown, part) in unknowns.iter().zip(&is_part) {
                let column = if *part {
                    derive(f_top, g_top, unknown).sub(&k_top.mul(unknown))
                } else {
                    leading.mul(unknown).scale(&minus_one)
                };
                row.push(column.coefficient(&monomial));
            }
            row.push(-known.coefficient(&monomial));
            rows.push(row);
        }
        let solution = solve_rational(rows, unknowns.len())?;
        for (index, unknown) in unknowns.iter().enumerate() {
            let value = solution.get(index)?;
            let piece = unknown.scale(value);
            if is_part[index] {
                let slot = part_level?;
                let entry = parts.entry(slot).or_default();
                *entry = entry.add(&piece);
            } else {
                let slot = cofactor_level?;
                let entry = cofactors.entry(slot).or_default();
                *entry = entry.add(&piece);
            }
        }
    }
    let mut polynomial = Planar::default();
    for part in parts.values() {
        polynomial = polynomial.add(part);
    }
    let mut cofactor = Planar::default();
    for part in cofactors.values() {
        cofactor = cofactor.add(part);
    }
    // The descent is an exact construction, but a coding slip must not be
    // reported as a Darboux polynomial.
    (derive(f, g, &polynomial) == cofactor.mul(&polynomial)).then_some((polynomial, cofactor))
}

/// Every Darboux polynomial of total degree at most `max_degree`, with its
/// cofactor, up to the bounding normalisation (a non-zero scalar multiple of
/// a Darboux polynomial is the same curve).
fn darboux_polynomials(
    f: &Planar,
    g: &Planar,
    max_degree: u32,
) -> Vec<(Planar, Planar)> {
    let top_degree = f.degree().unwrap_or(0).max(g.degree().unwrap_or(0));
    if top_degree == 0 {
        return Vec::new();
    }
    let cofactor_degree = top_degree.saturating_sub(1);
    let f_top = f.homogeneous_part(top_degree);
    let g_top = g.homogeneous_part(top_degree);
    // `R = x g_m - y f_m`; its linear factors are the invariant lines.
    let r = Planar::variable_x()
        .mul(&g_top)
        .sub(&Planar::variable_y().mul(&f_top));
    let Some(factors) = invariant_lines(&r) else {
        return Vec::new();
    };
    let mut found = Vec::new();
    for degree in 1..=max_degree {
        let mut leading = Vec::new();
        leading_forms(&factors, degree, &mut leading);
        let mut seen: Vec<Planar> = Vec::new();
        for candidate in leading {
            let (candidate, scale) = candidate.primitive_with_scale();
            if candidate.is_zero() || seen.contains(&candidate) {
                continue;
            }
            seen.push(candidate.clone());
            let Some((polynomial, cofactor)) =
                darboux_from_leading(f, g, &f_top, &g_top, &candidate, cofactor_degree)
            else {
                continue;
            };
            // The leading form was normalised; the whole construction is
            // scaled consistently already, but keep the invariant explicit.
            let _ = scale;
            let (primitive, factor) = polynomial.primitive_with_scale();
            found.push((primitive, cofactor.scale(&factor)));
        }
    }
    found
}

/// The planar polynomial system of two first-order autonomous equations,
/// with the dependent variables as plain symbols.
struct PlanarSystem {
    f: Planar,
    g: Planar,
    x: NodeId,
    y: NodeId,
}

/// Reads `list(x' = f, y' = g)` and `list(x(t), y(t))` as a planar
/// autonomous polynomial system.
fn planar_system(
    cx: &mut Cx<'_>,
    equations: NodeId,
    unknowns: NodeId,
) -> Option<PlanarSystem> {
    if cx.graph.op(equations) != core::LIST || cx.graph.op(unknowns) != core::LIST {
        return None;
    }
    let eqs = cx.graph.children(equations).to_vec();
    let funcs = cx.graph.children(unknowns).to_vec();
    // The first-order form does the derivative bookkeeping: it matches each
    // `x'(t)` to its rate, so equations written in any order (and higher
    // derivatives, which are then refused) are handled by the same tested
    // parsing as the scalar and system solvers.
    let system = first_order_form(cx, &eqs, &funcs)?;
    if system.symbols.len() != 2 || system.primary != [0, 1] || system.rates.len() != 2 {
        return None;
    }
    let time = cx.graph.symbol_of(system.t)?;
    for rate in &system.rates {
        if cx.graph.depends_on(cx.graph.find(*rate), time) {
            return None;
        }
    }
    let x_state = *system.symbols.first()?;
    let y_state = *system.symbols.get(1)?;
    let f = Planar::from_node(cx.graph, system.rates[0], x_state, y_state)?;
    let g = Planar::from_node(cx.graph, system.rates[1], x_state, y_state)?;
    // The report uses the plain symbols of the given dependent variables.
    let x = *cx.graph.children(*system.funcs.first()?).first()?;
    let y = *cx.graph.children(*system.funcs.get(1)?).first()?;
    Some(PlanarSystem { f, g, x, y })
}
/// The total degree searched, from an optional third argument.
fn requested_degree(
    cx: &Cx<'_>,
    args: &[NodeId],
) -> u32 {
    args.get(2)
        .and_then(|&node| cx.graph.number_of(node))
        .map(Number::to_f64)
        .filter(|value| value.is_finite() && *value >= 1.0)
        .map_or(DARBOUX_DEGREE, |value| {
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let degree = value.floor() as u32;
            degree.min(MAX_DARBOUX_DEGREE)
        })
}

/// Kernel for `ode_darboux(eqs, vars)`.
pub(super) struct Darboux {
    pub(super) op: OpId,
}

impl Kernel for Darboux {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let args = cx.graph.children(node).to_vec();
        let [equations, unknowns, ..] = args.as_slice() else {
            return Outcome::Pass;
        };
        let degree = requested_degree(cx, &args);
        let Some(system) = planar_system(cx, *equations, *unknowns) else {
            return Outcome::Pass;
        };
        let found = darboux_polynomials(&system.f, &system.g, degree);
        let items: Vec<NodeId> = found
            .iter()
            .map(|(polynomial, cofactor)| {
                let p = polynomial.to_node(cx.graph, system.x, system.y);
                let k = cofactor.to_node(cx.graph, system.x, system.y);
                cx.graph.node(core::LIST, &[p, k])
            })
            .collect();
        let answer = cx.graph.node(core::LIST, &items);
        Outcome::Pinned(answer)
    }
}

/// Kernel for `ode_jacobi_last_multiplier(eqs, vars)`.
pub(super) struct JacobiLastMultiplier {
    pub(super) op: OpId,
}

impl Kernel for JacobiLastMultiplier {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let args = cx.graph.children(node).to_vec();
        let [equations, unknowns, ..] = args.as_slice() else {
            return Outcome::Pass;
        };
        let degree = requested_degree(cx, &args);
        let Some(system) = planar_system(cx, *equations, *unknowns) else {
            return Outcome::Pass;
        };
        last_multiplier(cx, &system, degree).map_or(Outcome::Pass, Outcome::Pinned)
    }
}

/// `M = Π P_i^{λ_i}` with `Σ λ_i k_i = -div X`, if the Darboux polynomials up
/// to `max_degree` determine one.
fn last_multiplier(
    cx: &mut Cx<'_>,
    system: &PlanarSystem,
    max_degree: u32,
) -> Option<NodeId> {
    let found = darboux_polynomials(&system.f, &system.g, max_degree);
    if found.is_empty() {
        return None;
    }
    let divergence = system.f.derivative_x().add(&system.g.derivative_y());
    let target = divergence.scale(&rat(-1));
    let mut monomials: Vec<(u32, u32)> = Vec::new();
    for (_, cofactor) in &found {
        monomials.extend(cofactor.terms.keys().copied());
    }
    monomials.extend(target.terms.keys().copied());
    monomials.sort_unstable();
    monomials.dedup();
    let mut rows: Vec<Vec<Rat>> = Vec::with_capacity(monomials.len());
    for monomial in monomials {
        let mut row = Vec::with_capacity(found.len().saturating_add(1));
        for (_, cofactor) in &found {
            row.push(
                cofactor
                    .terms
                    .get(&monomial)
                    .cloned()
                    .unwrap_or_else(Rat::zero),
            );
        }
        row.push(
            target
                .terms
                .get(&monomial)
                .cloned()
                .unwrap_or_else(Rat::zero),
        );
        rows.push(row);
    }
    let lambdas = solve_rational(rows, found.len())?;
    let mut factors = Vec::new();
    for ((polynomial, _), lambda) in found.iter().zip(&lambdas) {
        if lambda.is_zero() {
            continue;
        }
        let base = polynomial.to_node(cx.graph, system.x, system.y);
        let exponent = cx.graph.num(Number::rat(lambda.clone()));
        factors.push(cx.graph.node(core::POW, &[base, exponent]));
    }
    if factors.is_empty() {
        return None;
    }
    let raw = mul(cx.graph, &factors);
    let answer = cx.simplify(raw);
    multiplier_holds(cx, answer, system).then_some(answer)
}

/// Numeric check of `X(M) + (div X) M = 0` (i.e. `d ln M/dt = -div X`) at
/// several points where the Darboux polynomials do not vanish.
fn multiplier_holds(
    cx: &mut Cx<'_>,
    multiplier: NodeId,
    system: &PlanarSystem,
) -> bool {
    let divergence = system.f.derivative_x().add(&system.g.derivative_y());
    let divergence_node = divergence.to_node(cx.graph, system.x, system.y);
    let f_node = system.f.to_node(cx.graph, system.x, system.y);
    let g_node = system.g.to_node(cx.graph, system.x, system.y);
    let (Some(dx), Some(dy)) = (
        derivative(cx.graph, multiplier, system.x),
        derivative(cx.graph, multiplier, system.y),
    ) else {
        return false;
    };
    let df = mul(cx.graph, &[dx, f_node]);
    let dg = mul(cx.graph, &[dy, g_node]);
    let image = add(cx.graph, &[df, dg]);
    let loss = mul(cx.graph, &[divergence_node, multiplier]);
    let residual = add(cx.graph, &[image, loss]);
    let mut checked = 0;
    for (x_value, y_value) in MULTIPLIER_CHECKS {
        let mut env = Env::numeric(0.0);
        if let (Some(x_symbol), Some(y_symbol)) =
            (cx.graph.symbol_of(system.x), cx.graph.symbol_of(system.y))
        {
            env.bind(x_symbol, x_value);
            env.bind(y_symbol, y_value);
        }
        let (Some(value), Some(scale)) = (
            cx.graph.eval(residual, &env),
            cx.graph.eval(multiplier, &env),
        ) else {
            continue;
        };
        if !value.is_finite() || !scale.is_finite() {
            continue;
        }
        if value.abs() > 1e-9 * (1.0 + scale.abs()) {
            return false;
        }
        checked += 1;
    }
    checked >= 2
}

#[cfg(test)]
mod darboux_tests {
    use super::*;
    use crate::graph::Engine;
    use std::fmt::Write as _;

    fn system(
        graph: &mut Graph,
        f: &str,
        g: &str,
    ) -> (Planar, Planar, NodeId, NodeId) {
        let x = graph.sym("x");
        let y = graph.sym("y");
        let f_node = graph.parse(f).expect("f");
        let g_node = graph.parse(g).expect("g");
        (
            Planar::from_node(graph, f_node, x, y).expect("f is polynomial"),
            Planar::from_node(graph, g_node, x, y).expect("g is polynomial"),
            x,
            y,
        )
    }

    fn find(
        found: &[(Planar, Planar)],
        want: &Planar,
    ) -> Option<(Planar, Planar)> {
        found
            .iter()
            .find(|(polynomial, _)| {
                polynomial.primitive_with_scale().0 == want.primitive_with_scale().0
            })
            .cloned()
    }

    /// `X(P) = k P` exactly, for every polynomial found.
    fn identities_hold(
        f: &Planar,
        g: &Planar,
        found: &[(Planar, Planar)],
    ) -> bool {
        found
            .iter()
            .all(|(polynomial, cofactor)| derive(f, g, polynomial) == cofactor.mul(polynomial))
    }

    /// `d ln M/dt = -div X`, checked numerically at points where `M` is
    /// finite: `X(ln M) + div X = 0`.
    fn multiplier_holds(
        graph: &mut Graph,
        multiplier: NodeId,
        f: &Planar,
        g: &Planar,
        x: NodeId,
        y: NodeId,
    ) -> bool {
        let divergence = f.derivative_x().add(&g.derivative_y());
        let divergence_node = divergence.to_node(graph, x, y);
        let f_node = f.to_node(graph, x, y);
        let g_node = g.to_node(graph, x, y);
        let ln = graph.ops().lookup("ln").expect("ln");
        let log = graph.node(ln, &[multiplier]);
        let (Some(dx), Some(dy)) = (derivative(graph, log, x), derivative(graph, log, y)) else {
            return false;
        };
        let df = mul(graph, &[dx, f_node]);
        let dg = mul(graph, &[dy, g_node]);
        let image = add(graph, &[df, dg]);
        let residual = add(graph, &[image, divergence_node]);
        let mut checked = 0;
        for (x_value, y_value) in [(1.3_f64, 0.7_f64), (2.1, 1.7), (0.4, 2.3)] {
            let mut env = Env::numeric(0.0);
            let (Some(x_symbol), Some(y_symbol)) = (graph.symbol_of(x), graph.symbol_of(y)) else {
                return false;
            };
            env.bind(x_symbol, x_value);
            env.bind(y_symbol, y_value);
            let Some(value) = graph.eval(residual, &env) else {
                continue;
            };
            if !value.is_finite() || value.abs() > 1e-8 {
                return false;
            }
            checked += 1;
        }
        checked >= 2
    }

    #[test]
    fn lotka_volterra_darboux_and_multiplier() {
        let mut graph = Graph::new();
        let (f, g, x, y) = system(&mut graph, "x*(1 - y)", "y*(x - 1)");
        let found = darboux_polynomials(&f, &g, 4);
        assert!(
            identities_hold(&f, &g, &found),
            "a reported Darboux identity failed"
        );
        // x with cofactor 1 - y and y with cofactor x - 1.
        let one = Planar::monomial((0, 0), one());
        let (px, kx) = find(&found, &Planar::variable_x()).expect("x is Darboux");
        assert_eq!(kx, one.sub(&Planar::variable_y()), "cofactor of x");
        let (py, ky) = find(&found, &Planar::variable_y()).expect("y is Darboux");
        assert_eq!(ky, Planar::variable_x().sub(&one), "cofactor of y");
        let _ = (px, py);
        // The Jacobi last multiplier: d ln M/dt = -div X with div X = x - y.
        let engine = Engine::install(&mut graph, &crate::rules::standard()).expect("install");
        let env = Env::symbolic();
        let mut cx = Cx {
            graph: &mut graph,
            env: &env,
            engine: &engine,
        };
        let planar = PlanarSystem { f, g, x, y };
        let multiplier = last_multiplier(&mut cx, &planar, 4).expect("a multiplier exists");
        assert!(
            multiplier_holds(cx.graph, multiplier, &planar.f, &planar.g, x, y),
            "the multiplier of Lotka-Volterra does not satisfy d ln M/dt = -div X"
        );
        // The multiplier is finite and non-zero at a sample point; with the
        // rational solution chosen, 1/(x y) comes out for this system.
        let mut env = Env::numeric(0.0);
        if let (Some(x_symbol), Some(y_symbol)) = (cx.graph.symbol_of(x), cx.graph.symbol_of(y)) {
            env.bind(x_symbol, 2.0);
            env.bind(y_symbol, 3.0);
        }
        let value = cx.graph.eval(multiplier, &env).expect("evaluate M");
        assert!(value.is_finite() && value != 0.0, "M(2, 3) = {value}");
    }

    #[test]
    fn linear_system_multiplier() {
        let mut graph = Graph::new();
        let (f, g, x, y) = system(&mut graph, "x", "2*y");
        let found = darboux_polynomials(&f, &g, 4);
        assert!(identities_hold(&f, &g, &found));
        let one = Planar::monomial((0, 0), one());
        let (_, kx) = find(&found, &Planar::variable_x()).expect("x is Darboux");
        assert_eq!(kx, one, "cofactor of x is 1");
        let (_, ky) = find(&found, &Planar::variable_y()).expect("y is Darboux");
        assert_eq!(ky, one.scale(&rat(2)), "cofactor of y is 2");
        let engine = Engine::install(&mut graph, &crate::rules::standard()).expect("install");
        let env = Env::symbolic();
        let mut cx = Cx {
            graph: &mut graph,
            env: &env,
            engine: &engine,
        };
        let planar = PlanarSystem { f, g, x, y };
        let multiplier = last_multiplier(&mut cx, &planar, 4).expect("a multiplier exists");
        assert!(
            multiplier_holds(cx.graph, multiplier, &planar.f, &planar.g, x, y),
            "the multiplier of the linear system does not satisfy d ln M/dt = -div X"
        );
        let mut env = Env::numeric(0.0);
        if let (Some(x_symbol), Some(y_symbol)) = (cx.graph.symbol_of(x), cx.graph.symbol_of(y)) {
            env.bind(x_symbol, 2.0);
            env.bind(y_symbol, 3.0);
        }
        let value = cx.graph.eval(multiplier, &env).expect("evaluate M");
        assert!(value.is_finite() && value != 0.0, "M(2, 3) = {value}");
    }

    #[test]
    fn no_low_degree_darboux_polynomial() {
        // Van der Pol has no invariant algebraic curve of degree up to four;
        // the search must report nothing rather than a wrong curve.
        let mut graph = Graph::new();
        let (f, g, _, _) = system(&mut graph, "y", "(1 - x^2)*y - x");
        assert!(
            darboux_polynomials(&f, &g, 4).is_empty(),
            "Van der Pol should have no Darboux polynomial of degree <= 4"
        );
    }

    #[test]
    fn darboux_operators_are_registered() {
        // The operators are additive API: installing the rule set registers
        // them, and the parsing/report path runs end to end. The kernels are
        // called directly so the test does not depend on how the engine
        // schedules a report operator.
        let mut graph = Graph::new();
        let engine = Engine::install(&mut graph, &crate::rules::standard()).expect("install");
        let (Some(darboux_op), Some(multiplier_op)) = (
            graph.ops().lookup("ode_darboux"),
            graph.ops().lookup("ode_jacobi_last_multiplier"),
        ) else {
            panic!("the Darboux operators are not registered");
        };
        let request = graph
            .parse("ode_darboux(list(diff(x(t), t) = x(t)*(1 - y(t)), diff(y(t), t) = y(t)*(x(t) - 1)), list(x(t), y(t)))")
            .expect("parse");
        let env = Env::symbolic();
        let mut cx = Cx {
            graph: &mut graph,
            env: &env,
            engine: &engine,
        };
        let kernel = Darboux { op: darboux_op };
        let outcome = kernel.reduce(&mut cx, request);
        let Outcome::Pinned(answer) = outcome else {
            let args = cx.graph.children(request).to_vec();
            let mut report = format!("args={}", args.len());
            if let [equations, unknowns] = args.as_slice() {
                let eqs = cx.graph.children(*equations).to_vec();
                let _ = write!(
                    report,
                    " eqs_op={} eqs_len={}",
                    cx.graph.ops().get(cx.graph.op(*equations)).name,
                    eqs.len()
                );
                for &equation in &eqs {
                    let _ = write!(
                        report,
                        " eq_op={}",
                        cx.graph.ops().get(cx.graph.op(equation)).name
                    );
                }
                let funcs = cx.graph.children(*unknowns).to_vec();
                let _ = write!(report, " funcs={}", funcs.len());
                for &func in &funcs {
                    let children = cx.graph.children(func);
                    let names: Vec<&str> = children
                        .iter()
                        .map(|&child| cx.graph.ops().get(cx.graph.op(child)).name.as_ref())
                        .collect();
                    let _ = write!(
                        report,
                        " func_op={} children={names:?}",
                        cx.graph.ops().get(cx.graph.op(func)).name
                    );
                }
                let _ = write!(
                    report,
                    " first_order={}",
                    first_order_form(&mut cx, &eqs, &funcs).is_some()
                );
            }
            panic!("ode_darboux did not report the Darboux polynomials: {report}");
        };
        let text = cx.graph.display(answer);
        assert!(text.contains("list(list("), "{text}");
        let request = cx.graph
            .parse("ode_jacobi_last_multiplier(list(diff(x(t), t) = x(t)*(1 - y(t)), diff(y(t), t) = y(t)*(x(t) - 1)), list(x(t), y(t)))")
            .expect("parse");
        let kernel = JacobiLastMultiplier { op: multiplier_op };
        let Outcome::Pinned(answer) = kernel.reduce(&mut cx, request) else {
            panic!("ode_jacobi_last_multiplier did not report a multiplier");
        };
        let text = cx.graph.display(answer);
        assert!(
            !text.contains("ode_") && (text.contains('x') && text.contains('y')),
            "{text}"
        );
        // A system with no Darboux polynomial of the degree searched returns
        // an empty list.
        let request = cx.graph
            .parse("ode_darboux(list(diff(x(t), t) = y(t), diff(y(t), t) = (1 - x(t)^2)*y(t) - x(t)), list(x(t), y(t)))")
            .expect("parse");
        let kernel = Darboux { op: darboux_op };
        let Outcome::Pinned(answer) = kernel.reduce(&mut cx, request) else {
            panic!("ode_darboux did not report an (empty) list for Van der Pol");
        };
        assert_eq!(cx.graph.display(answer), "list()");
    }
}
