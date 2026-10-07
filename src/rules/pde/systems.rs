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
//! When no combination of the matrices has distinct eigenvalues, a common
//! Jordan basis `P⁻¹ A_α P = λ_α I + Σ_k μ_{α,k} N^k` (upper triangular
//! Toeplitz blocks, `N` the nilpotent shift) decouples the system into
//! triangular chains: for first-order systems without sources on whole
//! space, `w_a = Σ_j q_j(t, ∂) S_{a+j}` with `S` the scalar solutions and
//! `q_j` the coefficients of `exp(t Σ_k N^k M_k)` — polynomials in `t`.
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
use super::util::sample;
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
        if cx.graph.op(n) == core::APPLY
            && let Some(&head) = cx.graph.children(n).first()
                && cx.graph.symbol_of(head).is_some() && !heads.contains(&head) {
                    heads.push(head);
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
    for g in candidates.clone() {
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
                    // Radicals may hide a zero from the simplifier: test numerically too.
                    let numerically_zero = (0..2_u32).all(|k| sample(cx.graph, entry, k).is_some_and(|v| v.abs() < 1e-9));
                    if i != j && !cx.is_zero(entry) && !numerically_zero {
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
    let mut jordan: Option<JordanData> = None;
    if found.is_none() && order == 1 && source.iter().all(|&e| cx.is_zero(e)) {
        // No distinct eigenvalues: Jordan chains (polynomial factors in t).
        if let Some((p, p_inverse, diagonals, data)) = jordan_decouple(cx, &a_matrices, &candidates) {
            found = Some((p, p_inverse, diagonals));
            jordan = Some(data);
        }
    }
    let (p_matrix, p_inverse, diagonals) = found?;
    // Conditions per component, grouped by shape.
    let shapes = match conditions {
        | Some(list) => group_conditions(cx, &problems, list)?,
        | None => Vec::new(),
    };
    if jordan.is_some() && shapes.iter().any(|shape| shape.on != time) {
        // Derivatives of the solution do not preserve boundary conditions.
        return None;
    }
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
    if let Some(data) = &jordan {
        solutions = jordan_combine(cx, data, &vars, vars[time], &solutions)?;
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

// ----------------------------------------------------------------------
// Jordan chains
// ----------------------------------------------------------------------

/// Whether `e` is zero, to rounding if the simplifier cannot tell.
fn vanishes(
    cx: &mut Cx<'_>,
    e: NodeId,
) -> bool {
    cx.is_zero(e) || (0..2_u32).all(|k| sample(cx.graph, e, k).is_some_and(|v| v.abs() < 1e-9))
}

/// Gauss–Jordan elimination of the rows: the reduced rows and their pivot
/// columns.
#[allow(clippy::needless_range_loop)]
fn row_reduce(
    cx: &mut Cx<'_>,
    mut rows: Matrix,
    columns: usize,
) -> (Matrix, Vec<usize>) {
    let mut pivots = Vec::new();
    let mut top = 0;
    for col in 0..columns {
        let Some(found) = (top..rows.len()).find(|&r| !vanishes(cx, rows[r][col])) else {
            continue;
        };
        rows.swap(top, found);
        let pivot = rows[top][col];
        for c in 0..columns {
            let q = div(cx, rows[top][c], pivot);
            rows[top][c] = cx.simplify(q);
        }
        for r in 0..rows.len() {
            if r == top {
                continue;
            }
            let factor = rows[r][col];
            if vanishes(cx, factor) {
                continue;
            }
            for c in 0..columns {
                let t = mul(cx.graph, &[factor, rows[top][c]]);
                let d = sub(cx.graph, rows[r][c], t);
                rows[r][c] = cx.simplify(d);
            }
        }
        pivots.push(col);
        top += 1;
    }
    (rows, pivots)
}

/// A basis of the null space of `m`.
fn null_space(
    cx: &mut Cx<'_>,
    m: &Matrix,
) -> Vec<Vec<NodeId>> {
    let n = m.len();
    let (rows, pivots) = row_reduce(cx, m.clone(), n);
    let mut out = Vec::new();
    for free in (0..n).filter(|c| !pivots.contains(c)) {
        let mut v: Vec<NodeId> = (0..n).map(|_| cx.graph.int(0)).collect();
        v[free] = cx.graph.int(1);
        for (i, &pc) in pivots.iter().enumerate() {
            v[pc] = neg(cx.graph, rows[i][free]);
        }
        out.push(v);
    }
    out
}

/// The blocks of a Jordan form: for each, the first column, the size and
/// the eigenvalue.
type Blocks = Vec<(usize, usize, NodeId)>;

/// A Jordan basis of `g`: the matrix `P` whose columns are the chains
/// (`g p₁ = λ p₁`, `g p_{j+1} = λ p_{j+1} + p_j`) and the blocks.
fn jordan_basis(
    cx: &mut Cx<'_>,
    g: &Matrix,
) -> Option<(Matrix, Blocks)> {
    let n = g.len();
    let lambda_symbol = cx.graph.interner_mut().fresh_symbol("lambda");
    let lambda = cx.graph.symbol_node(lambda_symbol);
    let characteristic = {
        let mut m = g.clone();
        for (i, row) in m.iter_mut().enumerate() {
            row[i] = sub(cx.graph, row[i], lambda);
        }
        determinant(cx, &m)
    };
    let found = solve_for(cx.graph, characteristic, lambda, 0)?;
    let mut roots: Vec<NodeId> = Vec::new();
    for r in found {
        let r = cx.simplify(r);
        let mut duplicate = false;
        for &q in &roots {
            let gap = sub(cx.graph, r, q);
            duplicate |= vanishes(cx, gap);
        }
        if !duplicate {
            roots.push(r);
        }
    }
    let mut columns: Vec<Vec<NodeId>> = Vec::new();
    let mut blocks: Blocks = Vec::new();
    for &root in &roots {
        let mut shifted = g.clone();
        for (i, row) in shifted.iter_mut().enumerate() {
            row[i] = sub(cx.graph, row[i], root);
        }
        let shifted: Matrix = shifted.iter().map(|row| row.iter().map(|&e| cx.simplify(e)).collect()).collect();
        // Kernels of the powers.
        let mut kernels: Vec<Vec<Vec<NodeId>>> = vec![Vec::new()];
        let mut power = shifted.clone();
        for j in 1..=n {
            if j > 1 {
                power = multiply(cx, &power, &shifted);
            }
            let k = null_space(cx, &power);
            let done = kernels.last().is_some_and(|prev| prev.len() == k.len());
            if done {
                break;
            }
            kernels.push(k);
        }
        let height = kernels.len() - 1;
        let dimension: Vec<usize> = kernels.iter().map(Vec::len).collect();
        let mut chains: Vec<Vec<Vec<NodeId>>> = Vec::new();
        for j in (1..=height).rev() {
            let below = dimension[j - 1];
            let here = dimension[j];
            let above = if j < height { dimension[j + 1] } else { here };
            let count = (here - below).checked_sub(above - here)?;
            // Vectors already independent: the lower kernel and the images of longer chains.
            let mut spanned: Vec<Vec<NodeId>> = kernels[j - 1].clone();
            for chain in &chains {
                // chain[0] is the top vector, of height `chain.len()`; its
                // image of height j is chain[chain.len() - j].
                let idx = chain.len().checked_sub(j)?;
                spanned.push(chain.get(idx)?.clone());
            }
            let mut chosen = 0;
            for candidate in &kernels[j] {
                if chosen == count {
                    break;
                }
                let mut trial = spanned.clone();
                trial.push(candidate.clone());
                let (_, pivots) = row_reduce(cx, trial, n);
                let (_, base) = row_reduce(cx, spanned.clone(), n);
                if pivots.len() > base.len() {
                    spanned.push(candidate.clone());
                    // The chain v, M v, ..., M^{j-1} v.
                    let mut chain = vec![candidate.clone()];
                    for _ in 1..j {
                        let last = chain.last()?.clone();
                        let next = apply(cx, &shifted, &last);
                        chain.push(next);
                    }
                    chains.push(chain);
                    chosen += 1;
                }
            }
            if chosen != count {
                return None;
            }
        }
        for chain in chains {
            blocks.push((columns.len(), chain.len(), root));
            // Ascending: the eigenvector first.
            for v in chain.into_iter().rev() {
                columns.push(v);
            }
        }
    }
    if columns.len() != n {
        return None;
    }
    let mut p = zero_matrix(cx, n);
    for (j, column) in columns.iter().enumerate() {
        for (i, &e) in column.iter().enumerate() {
            p[i][j] = e;
        }
    }
    Some((p, blocks))
}

/// The data of a decoupling into Jordan blocks: for each block, the
/// superdiagonal coefficients `μ_k` (`k ≥ 1`) of every `A_α`.
struct JordanData {
    blocks: Blocks,
    /// `couplings[block][k - 1]`: `(index, μ)` for each derivative.
    couplings: Vec<Vec<Vec<(Index, NodeId)>>>,
}

/// Simultaneous Jordan decoupling: `P⁻¹ A_α P` must be block diagonal with
/// upper triangular Toeplitz blocks. The diagonals `λ_{α,block}` are
/// returned like the eigenvalues of a diagonalisation.
#[allow(clippy::type_complexity, clippy::needless_range_loop)]
fn jordan_decouple(
    cx: &mut Cx<'_>,
    a_matrices: &[(Index, Matrix)],
    candidates: &[Matrix],
) -> Option<(Matrix, Matrix, Vec<Vec<NodeId>>, JordanData)> {
    for g in candidates {
        if is_zero_matrix(cx, g) {
            continue;
        }
        let Some((p, blocks)) = jordan_basis(cx, g) else {
            continue;
        };
        let Some(p_inverse) = inverse(cx, &p) else {
            continue;
        };
        let n = p.len();
        let mut diagonals = Vec::new();
        let mut couplings: Vec<Vec<Vec<(Index, NodeId)>>> = blocks.iter().map(|&(_, size, _)| vec![Vec::new(); size.saturating_sub(1)]).collect();
        let mut ok = true;
        'matrices: for (index, a) in a_matrices {
            let left = multiply(cx, &p_inverse, a);
            let t = multiply(cx, &left, &p);
            let mut diagonal = vec![cx.graph.int(0); n];
            for (b, &(start, size, _)) in blocks.iter().enumerate() {
                for i in 0..n {
                    for j in 0..n {
                        let entry = t[i][j];
                        let in_i = i >= start && i < start + size;
                        let in_j = j >= start && j < start + size;
                        if in_i && in_j {
                            continue;
                        }
                        if (in_i || in_j) && !vanishes(cx, entry) {
                            ok = false;
                            break 'matrices;
                        }
                    }
                }
                // Toeplitz structure inside the block.
                let mut values: Vec<Option<NodeId>> = vec![None; size];
                for a_idx in 0..size {
                    for b_idx in 0..size {
                        let entry = t[start + a_idx][start + b_idx];
                        if b_idx < a_idx {
                            if !vanishes(cx, entry) {
                                ok = false;
                                break 'matrices;
                            }
                            continue;
                        }
                        let k = b_idx - a_idx;
                        match values[k] {
                            | None => values[k] = Some(entry),
                            | Some(first) => {
                                let gap = sub(cx.graph, first, entry);
                                if !vanishes(cx, gap) {
                                    ok = false;
                                    break 'matrices;
                                }
                            },
                        }
                    }
                }
                let lambda_value = values[0]?;
                for slot in diagonal.iter_mut().skip(start).take(size) {
                    *slot = lambda_value;
                }
                for k in 1..size {
                    let mu = values[k]?;
                    if !vanishes(cx, mu) {
                        couplings[b][k - 1].push((index.clone(), mu));
                    }
                }
            }
            diagonals.push(diagonal);
        }
        if ok && blocks.iter().any(|&(_, size, _)| size > 1) {
            return Some((p, p_inverse, diagonals, JordanData { blocks, couplings }));
        }
    }
    None
}

/// A constant-coefficient differential operator: `(index, coefficient)`.
type Operator = Vec<(Index, NodeId)>;

fn multiply_operators(
    cx: &mut Cx<'_>,
    a: &Operator,
    b: &Operator,
) -> Operator {
    let mut out: Operator = Vec::new();
    for (ia, ca) in a {
        for (ib, cb) in b {
            let index: Index = ia.iter().zip(ib).map(|(x, y)| x + y).collect();
            let c = mul(cx.graph, &[*ca, *cb]);
            match out.iter_mut().find(|(i, _)| *i == index) {
                | Some((_, existing)) => *existing = add(cx.graph, &[*existing, c]),
                | None => out.push((index, c)),
            }
        }
    }
    for (_, c) in &mut out {
        *c = cx.simplify(*c);
    }
    out
}

/// The solution `w_a = Σ_j q_j(t, ∂) S_{a+j}` of a Jordan block, with `q_j`
/// the coefficients of `exp(t Σ_k N^k M_k)`.
fn jordan_combine(
    cx: &mut Cx<'_>,
    data: &JordanData,
    vars: &[NodeId],
    time_var: NodeId,
    solutions: &[NodeId],
) -> Option<Vec<NodeId>> {
    let dimension = vars.len();
    let mut out = solutions.to_vec();
    let one = cx.graph.int(1);
    for (b, &(start, size, _)) in data.blocks.iter().enumerate() {
        if size == 1 {
            continue;
        }
        let identity: Operator = vec![(vec![0; dimension], one)];
        // x_k = t M_k.
        let x: Vec<Operator> = data.couplings[b]
            .iter()
            .map(|terms| terms.iter().map(|(i, mu)| (i.clone(), mul(cx.graph, &[time_var, *mu]))).collect())
            .collect();
        let mut q: Vec<Operator> = vec![identity];
        for j in 1..size {
            // j q_j = Σ_{k=1}^{j} k x_k q_{j-k}.
            let mut total: Operator = Vec::new();
            for k in 1..=j {
                let xk = x.get(k - 1)?;
                let product = multiply_operators(cx, xk, &q[j - k]);
                let scale = cx.graph.int(i64::try_from(k).ok()?);
                for (i, c) in product {
                    let c = mul(cx.graph, &[scale, c]);
                    match total.iter_mut().find(|(t, _)| *t == i) {
                        | Some((_, e)) => *e = add(cx.graph, &[*e, c]),
                        | None => total.push((i, c)),
                    }
                }
            }
            let inverse_j = {
                let jn = cx.graph.int(i64::try_from(j).ok()?);
                div(cx, one, jn)
            };
            for (_, c) in &mut total {
                *c = mul(cx.graph, &[inverse_j, *c]);
                *c = cx.simplify(*c);
            }
            q.push(total);
        }
        for a in 0..size {
            let mut terms = Vec::new();
            for j in 0..(size - a) {
                let source = solutions[start + a + j];
                for (index, c) in &q[j] {
                    let mut d = source;
                    for (k, &order) in index.iter().enumerate() {
                        for _ in 0..order {
                            d = crate::rules::calculus::derivative(cx.graph, d, vars[k])?;
                        }
                    }
                    terms.push(mul(cx.graph, &[*c, d]));
                }
            }
            let total = add(cx.graph, &terms);
            out[start + a] = cx.simplify(total);
        }
    }
    Some(out)
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
