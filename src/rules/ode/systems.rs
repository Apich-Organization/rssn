//! Linear systems `X' = A X + F(t)` with a constant matrix `A`.
//!
//! `dsolve(list(diff(x(t), t) = …, …), list(x(t), y(t), …))`. The method
//! eliminates instead of diagonalising, so real solutions come out in real
//! form without complex arithmetic:
//!
//! 1. the characteristic polynomial `p(λ) = Σ c_k λ^k` of `A`
//!    (Faddeev–LeVerrier);
//! 2. by Cayley–Hamilton a component `x₁` satisfies the scalar equation
//!    `Σ c_k x₁^(k) = Σ c_k Σ_{j<k} e₁ᵀ A^(k-1-j) F^(j)`, solved by the scalar
//!    solver (constants `C1 … Cn`);
//! 3. the Krylov rows `e₁ᵀ A^k` determine the whole state from `x₁` and
//!    its derivatives: `K X = (x₁^(k) - Σ_{j<k} e₁ᵀ A^(k-1-j) F^(j))_k`;
//!    another component is used when `K` is singular.
//!
//! The answer is checked by substitution.

use super::add;
use super::mul;
use super::solve_equation;
use super::sub;
use crate::graph::Cx;
use crate::graph::Env;
use crate::graph::NodeId;
use crate::graph::op::core;
use crate::rules::calculus::derivative;
use crate::rules::solve::as_expression;
use crate::rules::solve::solve_linear;

type Matrix = Vec<Vec<NodeId>>;

fn mat_mul(
    cx: &mut Cx<'_>,
    a: &Matrix,
    b: &Matrix,
) -> Matrix {
    let n = a.len();
    let m = b.first().map_or(0, Vec::len);
    let mut out = vec![vec![NodeId::NONE; m]; n];
    for i in 0..n {
        for j in 0..m {
            let terms: Vec<NodeId> = (0..b.len()).map(|k| mul(cx.graph, &[a[i][k], b[k][j]])).collect();
            let sum = add(cx.graph, &terms);
            out[i][j] = cx.simplify(sum);
        }
    }
    out
}

pub(super) fn solve_system(
    cx: &mut Cx<'_>,
    equations: NodeId,
    unknowns: NodeId,
) -> Option<NodeId> {
    if cx.graph.op(equations) != core::LIST {
        return None;
    }
    let eqs = cx.graph.children(equations).to_vec();
    let funcs = cx.graph.children(unknowns).to_vec();
    let n = funcs.len();
    if n < 2 || eqs.len() != n || n > 4 {
        return None;
    }
    let diff = cx.graph.ops().lookup("diff")?;
    let t = *cx.graph.children(*funcs.first()?).get(1)?;
    cx.graph.symbol_of(t)?;
    // Stand-ins for x_i and x_i'.
    let mut states = Vec::with_capacity(n);
    let mut rates = Vec::with_capacity(n);
    for _ in 0..n {
        let s = cx.graph.interner_mut().fresh_symbol("s");
        states.push(cx.graph.symbol_node(s));
        let d = cx.graph.interner_mut().fresh_symbol("d");
        rates.push(cx.graph.symbol_node(d));
    }
    let mut residuals = Vec::with_capacity(n);
    for &e in &eqs {
        let mut r = as_expression(cx.graph, e);
        for i in 0..n {
            let d = cx.graph.node(diff, &[funcs[i], t]);
            r = cx.graph.replace_subterm(r, d, rates[i]);
            r = cx.graph.replace_subterm(r, funcs[i], states[i]);
        }
        residuals.push(cx.simplify(r));
    }
    let derivatives = solve_linear(cx.graph, &residuals, &rates)?;
    // A_ij = ∂f_i/∂s_j (constant), F_i = f_i at s = 0.
    let zero = cx.graph.int(0);
    let mut a: Matrix = vec![vec![zero; n]; n];
    let mut f = vec![zero; n];
    let t_symbol = cx.graph.symbol_of(t)?;
    for i in 0..n {
        for j in 0..n {
            let d = derivative(cx.graph, derivatives[i], states[j])?;
            let d = cx.simplify(d);
            let symbols = cx.graph.free_symbols(cx.graph.find(d)).to_vec();
            let state_symbols: Vec<_> = states.iter().filter_map(|&s| cx.graph.symbol_of(s)).collect();
            if symbols.iter().any(|s| *s == t_symbol || state_symbols.contains(s)) {
                return None;
            }
            a[i][j] = d;
        }
        let mut forcing = derivatives[i];
        for &s in &states {
            forcing = cx.graph.substitute(forcing, s, zero);
        }
        f[i] = cx.simplify(forcing);
    }
    // Characteristic polynomial by Faddeev–LeVerrier: c[n] = 1.
    let one = cx.graph.int(1);
    let identity: Matrix = (0..n).map(|i| (0..n).map(|j| if i == j { one } else { zero }).collect()).collect();
    let mut c = vec![zero; n + 1];
    c[n] = one;
    let mut m: Matrix = vec![vec![zero; n]; n];
    for k in 1..=n {
        let am = mat_mul(cx, &a, &m);
        m = am;
        for i in 0..n {
            let shifted = mul(cx.graph, &[c[n - k + 1], identity[i][i]]);
            let sum = add(cx.graph, &[m[i][i], shifted]);
            m[i][i] = cx.simplify(sum);
        }
        let am = mat_mul(cx, &a, &m);
        let trace = add(cx.graph, &(0..n).map(|i| am[i][i]).collect::<Vec<_>>());
        let scale = cx.graph.num(crate::graph::Number::fraction(-1, i64::try_from(k).ok()?)?);
        let value = mul(cx.graph, &[scale, trace]);
        c[n - k] = cx.simplify(value);
    }
    // Try each component as the pivot of the elimination.
    for pivot in 0..n {
        if let Some(answer) = eliminate(cx, &a, &f, &c, &funcs, t, pivot) {
            return Some(answer);
        }
    }
    None
}

fn eliminate(
    cx: &mut Cx<'_>,
    a: &Matrix,
    f: &[NodeId],
    c: &[NodeId],
    funcs: &[NodeId],
    t: NodeId,
    pivot: usize,
) -> Option<NodeId> {
    let n = funcs.len();
    let zero = cx.graph.int(0);
    let one = cx.graph.int(1);
    // Krylov rows e_pᵀ A^k, k = 0..n.
    let mut rows: Vec<Vec<NodeId>> = Vec::with_capacity(n + 1);
    rows.push((0..n).map(|j| if j == pivot { one } else { zero }).collect());
    for _ in 0..n {
        let last = vec![rows.last()?.clone()];
        let next = super::systems::mat_mul(cx, &last, a);
        rows.push(next.into_iter().next()?);
    }
    // Forcing corrections g_k = Σ_{j<k} e_pᵀ A^(k-1-j) F^(j).
    let mut f_derivatives = vec![f.to_vec()];
    for _ in 1..n {
        let last = f_derivatives.last()?.clone();
        let next: Option<Vec<NodeId>> = last.iter().map(|&v| derivative(cx.graph, v, t)).collect();
        f_derivatives.push(next?);
    }
    let mut corrections = Vec::with_capacity(n + 1);
    for k in 0..=n {
        let mut terms = Vec::new();
        for j in 0..k {
            let row = &rows[k - 1 - j];
            for (idx, &entry) in row.iter().enumerate() {
                let fj = *f_derivatives.get(j).and_then(|v| v.get(idx)).unwrap_or(&zero);
                terms.push(mul(cx.graph, &[entry, fj]));
            }
        }
        let sum = add(cx.graph, &terms);
        corrections.push(cx.simplify(sum));
    }
    // Scalar equation Σ c_k x^(k) = Σ c_k g_k for x = x_p.
    let diff = cx.graph.ops().lookup("diff")?;
    let x = funcs[pivot];
    let mut chain = vec![x];
    for _ in 0..n {
        let last = *chain.last()?;
        chain.push(cx.graph.node(diff, &[last, t]));
    }
    let lhs_terms: Vec<NodeId> = (0..=n).map(|k| mul(cx.graph, &[c[k], chain[k]])).collect();
    let rhs_terms: Vec<NodeId> = (0..=n).map(|k| mul(cx.graph, &[c[k], corrections[k]])).collect();
    let lhs = add(cx.graph, &lhs_terms);
    let rhs = add(cx.graph, &rhs_terms);
    let rhs = cx.simplify(rhs);
    let equation = cx.graph.node(core::EQ, &[lhs, rhs]);
    let (_, answer) = solve_equation(cx, equation, x, 1)?;
    let &[left, scalar] = cx.graph.children(answer) else {
        return None;
    };
    if left != x {
        return None;
    }
    // K X = (x^(k) - g_k)_{k<n}.
    let mut values = vec![scalar];
    for _ in 1..n {
        let last = *values.last()?;
        values.push(derivative(cx.graph, last, t)?);
    }
    let mut unknowns = Vec::with_capacity(n);
    for _ in 0..n {
        let s = cx.graph.interner_mut().fresh_symbol("X");
        unknowns.push(cx.graph.symbol_node(s));
    }
    let mut linear = Vec::with_capacity(n);
    for k in 0..n {
        let terms: Vec<NodeId> = (0..n).map(|j| mul(cx.graph, &[rows[k][j], unknowns[j]])).collect();
        let combination = add(cx.graph, &terms);
        let target = sub(cx.graph, values[k], corrections[k]);
        let equation = sub(cx.graph, combination, target);
        linear.push(cx.simplify(equation));
    }
    let solution = solve_linear(cx.graph, &linear, &unknowns)?;
    let mut items = Vec::with_capacity(n);
    for (i, &value) in solution.iter().enumerate() {
        let value = cx.simplify(value);
        items.push((funcs[i], value));
    }
    if !satisfies(cx, a, f, &items, t) {
        return None;
    }
    let equations: Vec<NodeId> = items.iter().map(|&(f, v)| cx.graph.node(core::EQ, &[f, v])).collect();
    Some(cx.graph.node(core::LIST, &equations))
}

/// `X' = A X + F` at a few points, with the constants at generic values.
fn satisfies(
    cx: &mut Cx<'_>,
    a: &Matrix,
    f: &[NodeId],
    items: &[(NodeId, NodeId)],
    t: NodeId,
) -> bool {
    let Some(t_symbol) = cx.graph.symbol_of(t) else {
        return false;
    };
    let n = items.len();
    let mut residuals = Vec::with_capacity(n);
    for i in 0..n {
        let Some(d) = derivative(cx.graph, items[i].1, t) else {
            return false;
        };
        let mut terms = vec![f[i]];
        for j in 0..n {
            terms.push(mul(cx.graph, &[a[i][j], items[j].1]));
        }
        let sum = add(cx.graph, &terms);
        residuals.push(sub(cx.graph, d, sum));
    }
    for point in [0.3, 0.9, 1.7] {
        for &r in &residuals {
            let mut env = Env::numeric(0.0);
            for &s in cx.graph.free_symbols(cx.graph.find(r)) {
                env.bind(s, if s == t_symbol { point } else { 0.5 + 0.13 * f64::from(s.raw() % 7) });
            }
            match cx.graph.eval(r, &env) {
                | Some(v) if v.is_finite() && v.abs() > 1e-6 => return false,
                | _ => {},
            }
        }
    }
    true
}
