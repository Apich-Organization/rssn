//! Linear systems of partial differential equations with constant
//! coefficients, by decoupling.
//!
//! `pdsolve(list(eq1, eq2, ...), list(u(x, t), v(x, t), ...), conditions)`
//! for a system
//!
//! `M ∂_t^o U + Σ_α C_α ∂^α U + S(x, t) = 0`
//!
//! with constant `n × n` matrices (`o = 1` or `2`, `M` invertible): the
//! matrices `A_α = -M⁻¹ C_α` must be simultaneously diagonalisable with
//! distinct eigenvalues for some combination, `P⁻¹ A_α P = diag(d_α)`. The
//! characteristic variables `W = P⁻¹ U` then satisfy `n` scalar equations
//! `∂_t^o w_i = Σ_α d_{α,i} ∂^α w_i + s_i`, solved by the scalar solvers
//! with the initial and boundary data transformed in the same way
//! (`w_i(data) = Σ_j (P⁻¹)_{ij} u_j(data)`), and `U = P W`. This covers
//!
//! * first-order hyperbolic systems `U_t + A U_x = S` (characteristic
//!   variables, `F_i(x - λ_i t)`),
//! * coupled heat, wave and Schrödinger equations with constant coupling
//!   matrices (`U_t = D U_xx + B U`, `U_tt = A U_xx`), in any dimension and
//!   on any domain the scalar solvers handle.
//!
//! The answer is `list(u(x, t) = …, v(x, t) = …)`.

use super::Condition;
use super::Conditions;
use super::Method;
use super::Problem;
use super::Index;
use super::solve;
use super::util::div;
use super::util::is_zero_number;
use crate::graph::Cx;
use crate::graph::NodeId;
use crate::graph::op::core;
use crate::rules::complex::build::add;
use crate::rules::complex::build::mul;
use crate::rules::complex::build::neg;
use crate::rules::complex::build::sub;
use crate::rules::poly::best;
use crate::rules::solve::solve_for;

type Matrix = Vec<Vec<NodeId>>;

fn zero_matrix(
    cx: &mut Cx<'_>,
    n: usize,
) -> Matrix {
    (0..n).map(|_| (0..n).map(|_| cx.graph.int(0)).collect()).collect()
}

fn minor(
    m: &Matrix,
    i: usize,
    j: usize,
) -> Matrix {
    m.iter()
        .enumerate()
        .filter(|(r, _)| *r != i)
        .map(|(_, row)| row.iter().enumerate().filter(|(c, _)| *c != j).map(|(_, &e)| e).collect())
        .collect()
}

fn determinant(
    cx: &mut Cx<'_>,
    m: &Matrix,
) -> NodeId {
    let n = m.len();
    if n == 0 {
        return cx.graph.int(1);
    }
    if n == 1 {
        return m[0][0];
    }
    let mut terms = Vec::new();
    for (j, &entry) in m[0].iter().enumerate() {
        if cx.is_zero(entry) {
            continue;
        }
        let sub_det = determinant(cx, &minor(m, 0, j));
        let signed = if j % 2 == 0 { sub_det } else { neg(cx.graph, sub_det) };
        terms.push(mul(cx.graph, &[entry, signed]));
    }
    let total = add(cx.graph, &terms);
    cx.simplify(total)
}

fn adjugate(
    cx: &mut Cx<'_>,
    m: &Matrix,
) -> Matrix {
    let n = m.len();
    let mut out = zero_matrix(cx, n);
    if n == 1 {
        out[0][0] = cx.graph.int(1);
        return out;
    }
    for i in 0..n {
        for j in 0..n {
            let d = determinant(cx, &minor(m, i, j));
            let signed = if (i + j) % 2 == 0 { d } else { neg(cx.graph, d) };
            // Transposed cofactors.
            if let Some(slot) = out.get_mut(j).and_then(|row| row.get_mut(i)) {
                *slot = signed;
            }
        }
    }
    out
}

fn multiply(
    cx: &mut Cx<'_>,
    a: &Matrix,
    b: &Matrix,
) -> Matrix {
    let n = a.len();
    let mut out = zero_matrix(cx, n);
    for i in 0..n {
        for j in 0..n {
            let terms: Vec<NodeId> = (0..n).map(|k| mul(cx.graph, &[a[i][k], b[k][j]])).collect();
            let total = add(cx.graph, &terms);
            out[i][j] = cx.simplify(total);
        }
    }
    out
}

fn apply(
    cx: &mut Cx<'_>,
    a: &Matrix,
    v: &[NodeId],
) -> Vec<NodeId> {
    a.iter()
        .map(|row| {
            let terms: Vec<NodeId> = row.iter().zip(v).map(|(&x, &y)| mul(cx.graph, &[x, y])).collect();
            let total = add(cx.graph, &terms);
            cx.simplify(total)
        })
        .collect()
}

fn inverse(
    cx: &mut Cx<'_>,
    m: &Matrix,
) -> Option<Matrix> {
    let d = determinant(cx, m);
    if cx.is_zero(d) {
        return None;
    }
    let adj = adjugate(cx, m);
    Some(
        adj.iter()
            .map(|row| {
                row.iter()
                    .map(|&e| {
                        let q = div(cx, e, d);
                        cx.simplify(q)
                    })
                    .collect()
            })
            .collect(),
    )
}

fn scale(
    cx: &mut Cx<'_>,
    m: &Matrix,
    k: NodeId,
) -> Matrix {
    m.iter()
        .map(|row| {
            row.iter()
                .map(|&e| {
                    let p = mul(cx.graph, &[k, e]);
                    cx.simplify(p)
                })
                .collect()
        })
        .collect()
}

fn is_zero_matrix(
    cx: &mut Cx<'_>,
    m: &Matrix,
) -> bool {
    m.iter().flatten().all(|&e| cx.is_zero(e))
}

/// Eigenvalues (distinct) and the matrix of eigenvectors of `g`.
fn eigen(
    cx: &mut Cx<'_>,
    g: &Matrix,
) -> Option<(Vec<NodeId>, Matrix)> {
    let n = g.len();
    let lambda_symbol = cx.graph.interner_mut().fresh_symbol("lambda");
    let lambda = cx.graph.symbol_node(lambda_symbol);
    let shifted = |cx: &mut Cx<'_>, value: NodeId| -> Matrix {
        let mut m = g.clone();
        for (i, row) in m.iter_mut().enumerate() {
            row[i] = sub(cx.graph, row[i], value);
        }
        m
    };
    let characteristic = {
        let m = shifted(cx, lambda);
        determinant(cx, &m)
    };
    let roots = solve_for(cx.graph, characteristic, lambda, 0)?;
    let roots: Vec<NodeId> = roots.into_iter().map(|r| cx.simplify(r)).collect();
    if roots.len() != n {
        return None;
    }
    for (i, &a) in roots.iter().enumerate() {
        for &b in &roots[i + 1..] {
            let gap = sub(cx.graph, a, b);
            if cx.is_zero(gap) {
                return None;
            }
        }
    }
    let mut columns: Vec<Vec<NodeId>> = Vec::new();
    for &root in &roots {
        let m = shifted(cx, root);
        let m: Matrix = m.iter().map(|row| row.iter().map(|&e| cx.simplify(e)).collect()).collect();
        let adj = adjugate(cx, &m);
        // A nonzero column of the adjugate spans the eigenspace.
        let column = (0..n).map(|j| adj.iter().map(|row| row[j]).collect::<Vec<NodeId>>()).find(|c| c.iter().any(|&e| !cx.is_zero(e)))?;
        columns.push(column);
    }
    let mut p = zero_matrix(cx, n);
    for (j, column) in columns.iter().enumerate() {
        for (i, &e) in column.iter().enumerate() {
            p[i][j] = e;
        }
    }
    Some((roots, p))
}

/// The jet `∂^α f`.
fn jet(
    cx: &mut Cx<'_>,
    f: NodeId,
    vars: &[NodeId],
    index: &[u32],
) -> Option<NodeId> {
    let diff = cx.graph.ops().lookup("diff")?;
    let mut out = f;
    for (k, &order) in index.iter().enumerate() {
        for _ in 0..order {
            out = cx.graph.node(diff, &[out, vars[k]]);
        }
    }
    Some(out)
}

/// Gives the arbitrary functions of component `i` their own names.
fn rename_functions(
    cx: &mut Cx<'_>,
    solution: NodeId,
    i: usize,
) -> NodeId {
    let mut heads = Vec::new();
    let mut stack = vec![solution];
    while let Some(n) = stack.pop() {
        if cx.graph.op(n) == core::APPLY {
            if let Some(&head) = cx.graph.children(n).first() {
                if cx.graph.symbol_of(head).is_some() && !heads.contains(&head) {
                    heads.push(head);
                }
            }
        }
        stack.extend_from_slice(cx.graph.children(n));
    }
    let mut out = solution;
    for head in heads {
        let Some(symbol) = cx.graph.symbol_of(head) else {
            continue;
        };
        let name = format!("{}_{}", cx.graph.interner().symbol_name(symbol), i + 1);
        let renamed = cx.graph.sym(&name);
        out = cx.graph.substitute(out, head, renamed);
    }
    out
}

/// Whether the node mentions any of the function symbols.
fn mentions(
    cx: &Cx<'_>,
    node: NodeId,
    functions: &[NodeId],
) -> bool {
    functions.iter().any(|&f| cx.graph.symbol_of(f).is_some_and(|s| cx.graph.depends_on(cx.graph.find(node), s)))
}

/// Solves a system; `None` if it is not of the supported kind.
#[allow(clippy::too_many_lines)]
pub(super) fn solve_system(
    cx: &mut Cx<'_>,
    equations: NodeId,
    unknowns: NodeId,
    conditions: Option<NodeId>,
) -> Option<NodeId> {
    let eq_list = best(cx.graph, equations)?;
    let eqs = cx.graph.children(eq_list).to_vec();
    let u_list = best(cx.graph, unknowns)?;
    let us = cx.graph.children(u_list).to_vec();
    let n = us.len();
    if n < 2 || eqs.len() != n {
        return None;
    }
    // One Problem per (equation, unknown).
    let mut problems: Vec<Vec<Problem>> = Vec::new();
    for &e in &eqs {
        let mut row = Vec::new();
        for &u in &us {
            row.push(Problem::parse(cx, e, u)?);
        }
        problems.push(row);
    }
    let vars = problems[0][0].vars.clone();
    let functions: Vec<NodeId> = problems[0].iter().map(|p| p.function).collect();
    for row in &problems {
        for (i, p) in row.iter().enumerate() {
            if p.nonlinear || p.vars.len() != vars.len() || p.vars.iter().zip(&vars).any(|(&a, &b)| !cx.graph.same(a, b)) || p.function != functions[i] {
                return None;
            }
        }
    }
    // Matrices per derivative index.
    let mut indices: Vec<Index> = Vec::new();
    for row in &problems {
        for p in row {
            for (index, c) in &p.linear {
                if is_zero_number(cx.graph, *c) {
                    continue;
                }
                if mentions(cx, *c, &functions) || !p.constant(cx.graph, *c) {
                    return None;
                }
                if !indices.contains(index) {
                    indices.push(index.clone());
                }
            }
        }
    }
    let mut matrices: Vec<(Index, Matrix)> = Vec::new();
    for index in &indices {
        let mut m = zero_matrix(cx, n);
        for k in 0..n {
            for i in 0..n {
                m[k][i] = problems[k][i].coefficient(cx.graph, index);
            }
        }
        matrices.push((index.clone(), m));
    }
    // The time variable and the order.
    let time = problems[0][0].time(cx.graph);
    let named_t = cx.graph.symbol_of(vars[time]).is_some_and(|s| cx.graph.interner().symbol_name(s) == "t");
    let has_order = |order: u32| indices.iter().any(|i| *i == problems[0][0].unit(time, order));
    let order = if has_order(2) {
        2
    } else if has_order(1) {
        1
    } else {
        return None;
    };
    if !named_t && order != 1 {
        return None;
    }
    let leading = problems[0][0].unit(time, order);
    let m_lead = matrices.iter().find(|(i, _)| *i == leading)?.1.clone();
    let m_inverse = inverse(cx, &m_lead)?;
    // A_α = -M⁻¹ C_α.
    let minus_one = cx.graph.int(-1);
    let mut a_matrices: Vec<(Index, Matrix)> = Vec::new();
    for (index, c) in &matrices {
        if *index == leading {
            continue;
        }
        if index[time] > 0 && !(order == 2 && *index == problems[0][0].unit(time, 1)) {
            return None;
        }
        let product = multiply(cx, &m_inverse, c);
        a_matrices.push((index.clone(), scale(cx, &product, minus_one)));
    }
    if a_matrices.is_empty() {
        return None;
    }
    // Source: the residuals with all unknowns set to zero.
    let zero = cx.graph.int(0);
    let mut source = Vec::new();
    for row in &problems {
        let mut residual = row[0].residual;
        for p in row {
            for &(_, node) in &p.jets {
                residual = cx.graph.replace_subterm(residual, node, zero);
            }
        }
        source.push(cx.simplify(residual));
    }
    let source = apply(cx, &m_inverse, &source);
    let source: Vec<NodeId> = source.into_iter().map(|s| neg(cx.graph, s)).collect();
    let source: Vec<NodeId> = source.into_iter().map(|s| cx.simplify(s)).collect();
    // A simultaneously diagonalising basis.
    let mut found: Option<(Matrix, Matrix, Vec<Vec<NodeId>>)> = None;
    let mut candidates: Vec<Matrix> = a_matrices.iter().map(|(_, m)| m.clone()).collect();
    if a_matrices.len() > 1 {
        let mut combination = zero_matrix(cx, n);
        for (w, (_, m)) in a_matrices.iter().enumerate() {
            let k = cx.graph.int(i64::try_from(w).ok()? + 1);
            let part = scale(cx, m, k);
            for i in 0..n {
                for j in 0..n {
                    let s = add(cx.graph, &[combination[i][j], part[i][j]]);
                    combination[i][j] = cx.simplify(s);
                }
            }
        }
        candidates.push(combination);
    }
    for g in candidates {
        if is_zero_matrix(cx, &g) {
            continue;
        }
        let Some((_, p)) = eigen(cx, &g) else {
            continue;
        };
        let Some(p_inverse) = inverse(cx, &p) else {
            continue;
        };
        let mut diagonals = Vec::new();
        let mut ok = true;
        for (_, a) in &a_matrices {
            let left = multiply(cx, &p_inverse, a);
            let d = multiply(cx, &left, &p);
            for (i, row) in d.iter().enumerate() {
                for (j, &entry) in row.iter().enumerate() {
                    if i != j && !cx.is_zero(entry) {
                        ok = false;
                    }
                }
            }
            diagonals.push((0..n).map(|i| d[i][i]).collect::<Vec<NodeId>>());
        }
        if ok {
            found = Some((p, p_inverse, diagonals));
            break;
        }
    }
    let (p_matrix, p_inverse, diagonals) = found?;
    // Conditions per component, grouped by shape.
    let shapes = match conditions {
        | Some(list) => group_conditions(cx, &problems, list)?,
        | None => Vec::new(),
    };
    let source_w = apply(cx, &p_inverse, &source);
    let mut solutions = Vec::new();
    for i in 0..n {
        let f_symbol = cx.graph.interner_mut().fresh_symbol("w");
        let f_node = cx.graph.symbol_node(f_symbol);
        let mut args = vec![f_node];
        args.extend(&vars);
        let unknown = cx.graph.node(core::APPLY, &args);
        // ∂_t^o w = Σ d ∂^α w + s.
        let lhs = jet(cx, unknown, &vars, &problems[0][0].unit(time, order))?;
        let mut rhs_terms = Vec::new();
        for ((index, _), diagonal) in a_matrices.iter().zip(&diagonals) {
            if cx.is_zero(diagonal[i]) {
                continue;
            }
            let j = jet(cx, unknown, &vars, index)?;
            rhs_terms.push(mul(cx.graph, &[diagonal[i], j]));
        }
        rhs_terms.push(source_w[i]);
        let rhs = add(cx.graph, &rhs_terms);
        let equation = cx.graph.node(core::EQ, &[lhs, rhs]);
        let scalar = Problem::parse(cx, equation, unknown)?;
        // Transformed conditions.
        let mut list = Vec::new();
        for shape in &shapes {
            let terms: Vec<NodeId> = (0..n).map(|j| mul(cx.graph, &[p_inverse[i][j], shape.values[j]])).collect();
            let total = add(cx.graph, &terms);
            let value = cx.simplify(total);
            list.push(Condition { on: shape.on, point: shape.point, derivative: shape.derivative.clone(), value, robin: None });
        }
        let found = solve(cx, &scalar, &Conditions(list), Method::Any)?;
        let &[_, rhs] = cx.graph.children(found) else {
            return None;
        };
        solutions.push(rename_functions(cx, rhs, i));
    }
    // U = P W.
    let mut result = Vec::new();
    for (j, &u) in us.iter().enumerate() {
        let terms: Vec<NodeId> = (0..n).map(|i| mul(cx.graph, &[p_matrix[j][i], solutions[i]])).collect();
        let total = add(cx.graph, &terms);
        let value = cx.simplify(total);
        result.push(cx.graph.node(core::EQ, &[u, value]));
    }
    Some(cx.graph.node(core::LIST, &result))
}

struct Shape {
    on: usize,
    point: NodeId,
    derivative: Index,
    values: Vec<NodeId>,
}

/// Parses the conditions of a system: each is `u_i(...) = value` for one
/// component; the same shape must be given for every component.
fn group_conditions(
    cx: &mut Cx<'_>,
    problems: &[Vec<Problem>],
    list: NodeId,
) -> Option<Vec<Shape>> {
    let n = problems.len();
    let list = best(cx.graph, list)?;
    if cx.graph.op(list) != core::LIST {
        return None;
    }
    let mut per_component: Vec<Vec<(usize, NodeId, Index, NodeId)>> = vec![Vec::new(); n];
    for condition in cx.graph.children(list).to_vec() {
        let &[target, value] = cx.graph.children(condition) else {
            return None;
        };
        if cx.graph.op(condition) != core::EQ {
            return None;
        }
        let mut placed = false;
        for (i, row) in problems.iter().enumerate() {
            if let Some((on, point, derivative)) = Conditions::atom_of(cx, &row[i], target) {
                per_component[i].push((on, point, derivative, value));
                placed = true;
                break;
            }
        }
        if !placed {
            return None;
        }
    }
    let first = per_component.first()?.clone();
    let mut shapes = Vec::new();
    for (on, point, derivative, value) in &first {
        let mut values = vec![*value];
        for component in per_component.iter().skip(1) {
            let found = component.iter().find(|(o, q, d, _)| o == on && d == derivative && cx.graph.same(*q, *point))?;
            values.push(found.3);
        }
        shapes.push(Shape { on: *on, point: *point, derivative: derivative.clone(), values });
    }
    if per_component.iter().any(|c| c.len() != first.len()) {
        return None;
    }
    Some(shapes)
}
