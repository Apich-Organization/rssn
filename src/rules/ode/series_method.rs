//! Series solutions, transform methods for initial value problems, and
//! linear recurrences.
//!
//! * `ode_series(equation, y(x), x0, N)`: the solution as a polynomial in
//!   `x - x0` through degree `N` with `C1 … Cn` the free initial
//!   coefficients. The coefficients are fixed one at a time from the
//!   Taylor coefficients of the residual at `x0`, each step being linear
//!   in the newest coefficient — so non-linear equations with analytic
//!   right-hand sides work at ordinary points as well.
//! * Initial value problems with constant coefficients whose forcing the
//!   general method cannot integrate (`heaviside`, `dirac`, piecewise
//!   inputs) are solved with the Laplace transform: transform the
//!   equation, insert the initial values, solve for `Y(s)` and invert.
//! * `rsolve_z(equation, y(n), list(y(0) = a, …))`: linear recurrences with
//!   constant coefficients by the z-transform.

use super::Problem;
use super::add;
use super::mul;
use super::parse;
use super::sub;
use crate::graph::Cx;
use crate::graph::Env;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::op::core;
use crate::rules::calculus::derivative;
use crate::rules::solve::as_expression;

/// `ode_series(equation, y(x), x0, N)`.
pub(super) fn series_solution(
    cx: &mut Cx<'_>,
    equation: NodeId,
    unknown: NodeId,
    x0: NodeId,
    degree: NodeId,
) -> Option<NodeId> {
    let n_max = usize::try_from(cx.graph.number_of(degree)?.to_i64()?).ok().filter(|&n| n <= 16)?;
    let problem: Problem = parse(cx.graph, equation, unknown)?;
    let order = problem.order();
    if n_max < order {
        return None;
    }
    let x = problem.x;
    // y = Σ c_k (x - x0)^k with symbols c_k.
    let offset = sub(cx.graph, x, x0);
    let mut coefficients = Vec::with_capacity(n_max + 1);
    let mut terms = Vec::with_capacity(n_max + 1);
    for k in 0..=n_max {
        let symbol = cx.graph.interner_mut().fresh_symbol("c");
        let c = cx.graph.symbol_node(symbol);
        coefficients.push(c);
        let e = cx.graph.int(i64::try_from(k).ok()?);
        let power = cx.graph.node(core::POW, &[offset, e]);
        terms.push(mul(cx.graph, &[c, power]));
    }
    let y = add(cx.graph, &terms);
    let mut residual = problem.expr;
    let mut current = y;
    for (k, &stand) in problem.stand.iter().enumerate() {
        if k > 0 {
            current = derivative(cx.graph, current, x)?;
        }
        residual = cx.graph.substitute(residual, stand, current);
    }
    // The free coefficients become C1 … Cn.
    let mut values: Vec<Option<NodeId>> = vec![None; n_max + 1];
    for (k, value) in values.iter_mut().enumerate().take(order) {
        *value = Some(cx.graph.sym(&format!("C{}", k + 1)));
    }
    let mut factorial = 1_i64;
    let mut taylor = residual;
    for m in 0..=(n_max - order) {
        if m > 0 {
            taylor = derivative(cx.graph, taylor, x)?;
            factorial = factorial.checked_mul(i64::try_from(m).ok()?)?;
        }
        // The m-th Taylor coefficient of the residual with the known
        // coefficients inserted; it is linear in c_{m+order}.
        let mut at_point = cx.graph.substitute(taylor, x, x0);
        for (k, value) in values.iter().enumerate() {
            if let Some(v) = value {
                at_point = cx.graph.substitute(at_point, coefficients[k], *v);
            }
        }
        // Coefficients beyond the target do not enter this order.
        let zero = cx.graph.int(0);
        for &later in coefficients.iter().skip(m + order + 1) {
            at_point = cx.graph.substitute(at_point, later, zero);
        }
        let at_point = cx.simplify(at_point);
        let target = coefficients[m + order];
        let solutions = super::solve_for(cx.graph, at_point, target, 0)?;
        let value = cx.simplify(*solutions.first()?);
        values[m + order] = Some(value);
        let _ = factorial;
    }
    let mut series = Vec::with_capacity(n_max + 1);
    for (k, value) in values.into_iter().enumerate() {
        let value = value?;
        let e = cx.graph.int(i64::try_from(k).ok()?);
        let power = cx.graph.node(core::POW, &[offset, e]);
        series.push(mul(cx.graph, &[value, power]));
    }
    let sum = add(cx.graph, &series);
    let sum = cx.simplify(sum);
    Some(cx.graph.node(core::EQ, &[unknown, sum]))
}

/// The value given for `y^(k)(0)` in conditions `y(0) = a`,
/// `diff(y(t), t) = b`, … (derivative conditions at the point of the
/// condition before them); `None` unless every point is zero.
fn initial_values(
    cx: &Cx<'_>,
    unknown: NodeId,
    conditions: NodeId,
) -> Option<Vec<(usize, NodeId)>> {
    let diff = cx.graph.ops().lookup("diff")?;
    let function = *cx.graph.children(unknown).first()?;
    let mut out = Vec::new();
    for condition in cx.graph.children(conditions).to_vec() {
        let &[target, value] = cx.graph.children(condition) else {
            return None;
        };
        if cx.graph.op(condition) != core::EQ {
            return None;
        }
        match *cx.graph.children(target) {
            | [f, a] if cx.graph.op(target) == core::APPLY && f == function => {
                if !cx.graph.number_of(a).is_some_and(Number::is_zero) {
                    return None;
                }
                out.push((0, value));
            },
            | _ => {
                let mut order = 0;
                let mut inner = target;
                while cx.graph.op(inner) == diff {
                    inner = *cx.graph.children(inner).first()?;
                    order += 1;
                }
                if inner != unknown || order == 0 {
                    return None;
                }
                out.push((order, value));
            },
        }
    }
    Some(out)
}

/// An initial value problem at `t = 0` by the Laplace transform.
pub(super) fn laplace_ivp(
    cx: &mut Cx<'_>,
    equation: NodeId,
    unknown: NodeId,
    conditions: NodeId,
) -> Option<NodeId> {
    let (laplace, inverse, at, diff) = (
        cx.graph.ops().lookup("laplace")?,
        cx.graph.ops().lookup("inverse_laplace")?,
        cx.graph.ops().lookup("at")?,
        cx.graph.ops().lookup("diff")?,
    );
    let values = initial_values(cx, unknown, conditions)?;
    let t = *cx.graph.children(unknown).get(1)?;
    let s_symbol = cx.graph.interner_mut().fresh_symbol("s");
    let s = cx.graph.symbol_node(s_symbol);
    cx.graph.assume(s_symbol, crate::graph::Facts::POSITIVE);
    let expr = as_expression(cx.graph, equation);
    let pieces = transform_terms(cx, laplace, expr, t, s);
    let mut transformed = add(cx.graph, &pieces);
    // Initial values: at(diff^k(y(t)), t, 0) and y(0).
    let zero = cx.graph.int(0);
    let mut chain = vec![unknown];
    for _ in 0..values.iter().map(|v| v.0).max().unwrap_or(0) {
        let last = *chain.last()?;
        chain.push(cx.graph.node(diff, &[last, t]));
    }
    for &(k, value) in &values {
        let node = *chain.get(k)?;
        let marker = cx.graph.node(at, &[node, t, zero]);
        transformed = cx.graph.replace_subterm(transformed, marker, value);
        if k == 0 {
            let function = *cx.graph.children(unknown).first()?;
            let y0 = cx.graph.node(core::APPLY, &[function, zero]);
            transformed = cx.graph.replace_subterm(transformed, y0, value);
        }
    }
    let big_y = cx.graph.node(laplace, &[unknown, t, s]);
    let y_symbol = cx.graph.interner_mut().fresh_symbol("Y");
    let y_node = cx.graph.symbol_node(y_symbol);
    let algebraic = cx.graph.replace_subterm(transformed, big_y, y_node);
    let algebraic = cx.simplify(algebraic);
    let solved = *super::solve_for(cx.graph, algebraic, y_node, 0)?.first()?;
    let solved = cx.simplify(solved);
    let back = cx.graph.node(inverse, &[solved, s, t]);
    let back = cx.simplify(back);
    if contains_op(cx, back, &[laplace, inverse, at]) {
        return None;
    }
    Some(cx.graph.node(core::EQ, &[unknown, back]))
}

/// `ode_fourier(eq, y(x))`: the solution of a linear constant-coefficient
/// equation on the whole line that decays at `±∞`, by the Fourier
/// transform: `P(I w) Y(w) = F(w)`, `y = F⁻¹(F/P(I w))`.
pub(super) fn fourier_solution(
    cx: &mut Cx<'_>,
    equation: NodeId,
    unknown: NodeId,
) -> Option<NodeId> {
    let (fourier, inverse) = (
        cx.graph.ops().lookup("fourier")?,
        cx.graph.ops().lookup("inverse_fourier")?,
    );
    let x = *cx.graph.children(unknown).get(1)?;
    let w_symbol = cx.graph.interner_mut().fresh_symbol("w");
    let w = cx.graph.symbol_node(w_symbol);
    cx.graph.assume(w_symbol, crate::graph::Facts::REAL);
    let expr = as_expression(cx.graph, equation);
    let pieces = transform_terms(cx, fourier, expr, x, w);
    let transformed = add(cx.graph, &pieces);
    let big_y = cx.graph.node(fourier, &[unknown, x, w]);
    let y_symbol = cx.graph.interner_mut().fresh_symbol("Y");
    let y_node = cx.graph.symbol_node(y_symbol);
    let algebraic = cx.graph.replace_subterm(transformed, big_y, y_node);
    let algebraic = cx.simplify(algebraic);
    if contains_op(cx, algebraic, &[fourier]) {
        return None;
    }
    let solved = *super::solve_for(cx.graph, algebraic, y_node, 0)?.first()?;
    let solved = cx.simplify(solved);
    let back = cx.graph.node(inverse, &[solved, w, x]);
    let back = cx.simplify(back);
    if contains_op(cx, back, &[fourier, inverse]) {
        return None;
    }
    Some(cx.graph.node(core::EQ, &[unknown, back]))
}

/// `T(expr)` term by term (the transforms are linear), each piece
/// simplified on its own so nested runs stay small.
fn transform_terms(
    cx: &mut Cx<'_>,
    op: crate::graph::OpId,
    expr: NodeId,
    from: NodeId,
    to: NodeId,
) -> Vec<NodeId> {
    let expr = cx.simplify(expr);
    let expr = crate::rules::poly::expand_form(cx.graph, expr).unwrap_or(expr);
    let terms = if cx.graph.op(expr) == core::ADD { cx.graph.children(expr).to_vec() } else { vec![expr] };
    terms
        .into_iter()
        .map(|t| {
            let request = cx.graph.node(op, &[t, from, to]);
            cx.simplify(request)
        })
        .collect()
}

fn contains_op(
    cx: &Cx<'_>,
    node: NodeId,
    ops: &[crate::graph::OpId],
) -> bool {
    let mut stack = vec![node];
    while let Some(n) = stack.pop() {
        if ops.contains(&cx.graph.op(n)) {
            return true;
        }
        stack.extend_from_slice(cx.graph.children(n));
    }
    false
}

/// `rsolve(equation, y(n), conditions)`: linear constant-coefficient
/// recurrences by the z-transform.
pub(super) fn rsolve(
    cx: &mut Cx<'_>,
    equation: NodeId,
    unknown: NodeId,
    conditions: Option<NodeId>,
) -> Option<NodeId> {
    let (ztransform, inverse) = (cx.graph.ops().lookup("ztransform")?, cx.graph.ops().lookup("inverse_ztransform")?);
    let n = *cx.graph.children(unknown).get(1)?;
    let function = *cx.graph.children(unknown).first()?;
    cx.graph.symbol_of(n)?;
    let z_symbol = cx.graph.interner_mut().fresh_symbol("z");
    let z = cx.graph.symbol_node(z_symbol);
    let expr = as_expression(cx.graph, equation);
    // Shift so that the lowest index is y(n): y(n - 1) terms would need
    // values at negative indices.
    let pieces = transform_terms(cx, ztransform, expr, n, z);
    let mut transformed = add(cx.graph, &pieces);
    // Initial values y(0), y(1), …; without conditions they stay as
    // constants C1, C2, ….
    let mut given = Vec::new();
    if let Some(list) = conditions {
        for condition in cx.graph.children(list).to_vec() {
            let &[target, value] = cx.graph.children(condition) else {
                return None;
            };
            match *cx.graph.children(target) {
                | [f, k] if cx.graph.op(target) == core::APPLY && f == function => given.push((k, value)),
                | _ => return None,
            }
        }
    }
    for k in 0..8_i64 {
        let index = cx.graph.int(k);
        let y_k = cx.graph.node(core::APPLY, &[function, index]);
        let value = given
            .iter()
            .find(|(i, _)| cx.graph.number_of(*i).and_then(Number::to_i64) == Some(k))
            .map_or_else(|| cx.graph.sym(&format!("C{}", k + 1)), |&(_, v)| v);
        transformed = cx.graph.replace_subterm(transformed, y_k, value);
    }
    let big_x = cx.graph.node(ztransform, &[unknown, n, z]);
    let x_symbol = cx.graph.interner_mut().fresh_symbol("X");
    let x_node = cx.graph.symbol_node(x_symbol);
    let algebraic = cx.graph.replace_subterm(transformed, big_x, x_node);
    let algebraic = cx.simplify(algebraic);
    let solved = *super::solve_for(cx.graph, algebraic, x_node, 0)?.first()?;
    let solved = cx.simplify(solved);
    let back = cx.graph.node(inverse, &[solved, z, n]);
    let back = cx.simplify(back);
    if contains_op(cx, back, &[ztransform, inverse]) {
        return None;
    }
    // Check the recurrence at a few indices: every y(n + k) replaced by
    // the closed form at n + k.
    let mut check = expr;
    let mut stack = vec![expr];
    let mut applications = Vec::new();
    while let Some(node) = stack.pop() {
        if cx.graph.op(node) == core::APPLY && cx.graph.children(node).first() == Some(&function) {
            applications.push(node);
            continue;
        }
        stack.extend_from_slice(cx.graph.children(node));
    }
    for application in applications {
        let index = *cx.graph.children(application).get(1)?;
        let shifted = cx.graph.substitute(back, n, index);
        check = cx.graph.replace_subterm(check, application, shifted);
    }
    let n_symbol = cx.graph.symbol_of(n)?;
    for m in [3.0, 5.0, 8.0] {
        let mut env = Env::numeric(0.0);
        for &sym in cx.graph.free_symbols(cx.graph.find(check)) {
            env.bind(sym, if sym == n_symbol { m } else { 0.7 + 0.11 * f64::from(sym.raw() % 5) });
        }
        if let Some(v) = cx.graph.eval(check, &env)
            && v.is_finite() && v.abs() > 1e-7 {
                return None;
            }
    }
    Some(cx.graph.node(core::EQ, &[unknown, back]))
}
