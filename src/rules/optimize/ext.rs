//! Linear and quadratic programming, KKT points for inequality
//! constraints, convexity and Legendre transforms.
//!
//! | operator | value |
//! |---|---|
//! | `linprog(c, A, b)` | the exact minimum of `c·x` subject to `A x <= b`, `x >= 0` as `list(value, list(x1, ...))` (two-phase simplex over the rationals with Bland's rule); the symbols `infeasible` or `unbounded` otherwise |
//! | `linprog_max(c, A, b)` | the same for the maximum |
//! | `linprog_eq(c, A, b)` | the minimum of `c·x` subject to `A x = b`, `x >= 0` |
//! | `kkt_points(f, list(g1, ...), list(x, ...))` | `list(list(point, multipliers, value), ...)`: the Karush–Kuhn–Tucker points of `min f` subject to `g_i <= 0` (each `g_i` may also be written `le(a, b)` or `ge(a, b)`), found by solving the stationarity system for every active set; multipliers are `>= 0` and listed for every constraint, zero for inactive ones; sorted by value |
//! | `kkt_minimum(f, list(g1, ...), list(x, ...))` | `list(point, value)`, the KKT point of least value (the constrained minimum whenever one exists, e.g. on a compact feasible set) |
//! | `is_convex(f, list(x, ...))`, `is_concave(f, ...)` | `true` when every principal minor of the Hessian (of `-f` for concave) is provably non-negative; `false` when an exact sample point of the Hessian is not positive semidefinite; unreduced otherwise |
//! | `hessian_definiteness(f, list(x, ...), list(p, ...))` | one of `positive_definite`, `positive_semidefinite`, `negative_definite`, `negative_semidefinite`, `indefinite` for the Hessian at the point (exact when its entries are rational) |
//! | `qp_eq(Q, c, A, b)` | the minimiser and multipliers `list(x, lambda)` of `1/2 x'Qx + c'x` subject to `A x = b`, from the KKT linear system |
//! | `legendre_transform(f, x, p)` | `f*(p) = p x* - f(x*)` where `f'(x*) = p` |

use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::Zero;

use super::leading_minors;
use super::list_items;
use super::number;
use super::solve_system;
use super::substitute_point;
use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::RuleError;
use crate::graph::Tier;
use crate::rules::calculus::derivative;
use crate::rules::poly::best;
use crate::rules::solve::solve_for;

type Q = BigRational;

/// Largest number of inequality constraints `kkt_points` enumerates over.
const MAX_CONSTRAINTS: usize = 10;

// ----------------------------------------------------------------------
// Exact simplex
// ----------------------------------------------------------------------

/// Outcome of a linear program.
#[derive(Debug, PartialEq, Eq)]
enum Lp {
    Optimal(Q, Vec<Q>),
    Infeasible,
    Unbounded,
}

#[derive(Copy, Clone, PartialEq, Eq)]
enum Relation {
    AtMost,
    Equal,
}

struct Tableau {
    rows: Vec<Vec<Q>>,
    basis: Vec<usize>,
    objective: Vec<Q>,
}

impl Tableau {
    fn columns(&self) -> usize {
        self.objective.len() - 1
    }

    fn pivot(
        &mut self,
        row: usize,
        col: usize,
    ) {
        let pivot = self.rows[row][col].clone();
        for value in &mut self.rows[row] {
            *value = &*value / &pivot;
        }
        let pivot_row = self.rows[row].clone();
        for (r, other) in self.rows.iter_mut().enumerate() {
            if r != row && !other[col].is_zero() {
                let factor = other[col].clone();
                for (value, p) in other.iter_mut().zip(&pivot_row) {
                    *value -= &factor * p;
                }
            }
        }
        let factor = self.objective[col].clone();
        if !factor.is_zero() {
            for (value, p) in self.objective.iter_mut().zip(&pivot_row) {
                *value -= &factor * p;
            }
        }
        self.basis[row] = col;
    }

    /// Simplex iterations (Bland's rule) over the first `allowed` columns;
    /// `false` when unbounded.
    fn optimise(
        &mut self,
        allowed: usize,
    ) -> bool {
        loop {
            let Some(col) = (0..allowed).find(|&j| self.objective[j].is_negative()) else {
                return true;
            };
            let rhs = self.columns();
            let mut best: Option<(usize, Q)> = None;
            for (r, row) in self.rows.iter().enumerate() {
                if row[col].is_positive() {
                    let ratio = &row[rhs] / &row[col];
                    let better = match &best {
                        | None => true,
                        | Some((at, current)) => {
                            ratio < *current || (ratio == *current && self.basis[r] < self.basis[*at])
                        },
                    };
                    if better {
                        best = Some((r, ratio));
                    }
                }
            }
            let Some((row, _)) = best else {
                return false;
            };
            self.pivot(row, col);
        }
    }

    /// Resets the objective row to the costs `c` (zero beyond them),
    /// expressed in the current basis.
    fn price(
        &mut self,
        costs: &[Q],
    ) {
        let width = self.objective.len();
        let mut objective: Vec<Q> = (0..width).map(|j| costs.get(j).cloned().unwrap_or_else(Q::zero)).collect();
        objective[width - 1] = Q::zero();
        for (r, &b) in self.basis.iter().enumerate() {
            let cost = costs.get(b).cloned().unwrap_or_else(Q::zero);
            if !cost.is_zero() {
                for (value, entry) in objective.iter_mut().zip(&self.rows[r]) {
                    *value -= &cost * entry;
                }
            }
        }
        self.objective = objective;
    }
}

/// Minimises `c·x` over `x >= 0` subject to the rows `a x (<=|=) b`.
fn simplex(
    c: &[Q],
    a: &[Vec<Q>],
    b: &[Q],
    relations: &[Relation],
) -> Lp {
    let (n, m) = (c.len(), a.len());
    // Columns: variables, slacks, artificials.
    let slack_rows: Vec<usize> = (0..m).filter(|&i| relations[i] == Relation::AtMost).collect();
    let slack_start = n;
    let artificial_start = n + slack_rows.len();
    let mut rows: Vec<Vec<Q>> = Vec::with_capacity(m);
    let mut basis = Vec::with_capacity(m);
    let mut artificial_rows = Vec::new();
    for i in 0..m {
        let flip = b[i].is_negative();
        let sign = if flip { -Q::one() } else { Q::one() };
        let mut row: Vec<Q> = a[i].iter().map(|v| v * &sign).collect();
        row.resize(artificial_start, Q::zero());
        if let Some(k) = slack_rows.iter().position(|&r| r == i) {
            row[slack_start + k] = sign.clone();
        }
        let needs_artificial = relations[i] == Relation::Equal || flip;
        if needs_artificial {
            artificial_rows.push(i);
        } else if let Some(k) = slack_rows.iter().position(|&r| r == i) {
            basis.push(slack_start + k);
        }
        row.push(&b[i] * &sign);
        rows.push(row);
    }
    let total = artificial_start + artificial_rows.len();
    for row in &mut rows {
        let rhs = row.pop().unwrap_or_else(Q::zero);
        row.resize(total, Q::zero());
        row.push(rhs);
    }
    let mut basis_full = vec![0; m];
    let mut next_basic = basis.into_iter();
    for (i, slot) in basis_full.iter_mut().enumerate() {
        if let Some(k) = artificial_rows.iter().position(|&r| r == i) {
            rows[i][artificial_start + k] = Q::one();
            *slot = artificial_start + k;
        } else {
            *slot = next_basic.next().unwrap_or(0);
        }
    }
    let mut tableau = Tableau { rows, basis: basis_full, objective: vec![Q::zero(); total + 1] };
    if !artificial_rows.is_empty() {
        let mut costs = vec![Q::zero(); total];
        for cost in costs.iter_mut().skip(artificial_start) {
            *cost = Q::one();
        }
        tableau.price(&costs);
        tableau.optimise(total);
        if !tableau.objective[total].is_zero() {
            return Lp::Infeasible;
        }
        // Drive remaining artificial variables out of the basis.
        let mut r = 0;
        while r < tableau.rows.len() {
            if tableau.basis[r] >= artificial_start {
                if let Some(col) = (0..artificial_start).find(|&j| !tableau.rows[r][j].is_zero()) {
                    tableau.pivot(r, col);
                } else {
                    tableau.rows.remove(r);
                    tableau.basis.remove(r);
                    continue;
                }
            }
            r += 1;
        }
        for row in &mut tableau.rows {
            let rhs = row.pop().unwrap_or_else(Q::zero);
            row.truncate(artificial_start);
            row.push(rhs);
        }
        tableau.objective.truncate(artificial_start);
        tableau.objective.push(Q::zero());
    }
    let width = tableau.objective.len() - 1;
    tableau.price(c);
    if !tableau.optimise(width) {
        return Lp::Unbounded;
    }
    let mut point = vec![Q::zero(); n];
    for (r, &b) in tableau.basis.iter().enumerate() {
        if b < n {
            point[b] = tableau.rows[r][width].clone();
        }
    }
    let value = c.iter().zip(&point).map(|(ci, xi)| ci * xi).fold(Q::zero(), |acc, v| acc + v);
    Lp::Optimal(value, point)
}

fn rational(
    graph: &mut Graph,
    node: NodeId,
) -> Option<Q> {
    let node = best(graph, node)?;
    graph.number_of(node).and_then(Number::to_rational)
}

fn rational_vector(
    graph: &mut Graph,
    node: NodeId,
) -> Option<Vec<Q>> {
    let node = best(graph, node)?;
    if graph.op(node) != core::LIST {
        return None;
    }
    graph.children(node).to_vec().into_iter().map(|c| rational(graph, c)).collect()
}

fn rational_matrix(
    graph: &mut Graph,
    node: NodeId,
) -> Option<Vec<Vec<Q>>> {
    let node = best(graph, node)?;
    if graph.op(node) != core::LIST {
        return None;
    }
    let rows: Option<Vec<Vec<Q>>> =
        graph.children(node).to_vec().into_iter().map(|r| rational_vector(graph, r)).collect();
    let rows = rows?;
    let width = rows.first().map_or(0, Vec::len);
    rows.iter().all(|r| r.len() == width).then_some(rows)
}

fn rational_node(
    graph: &mut Graph,
    value: &Q,
) -> NodeId {
    graph.num(Number::rat(value.clone()))
}

fn linprog(
    cx: &mut Cx<'_>,
    args: &[NodeId],
    maximise: bool,
    equality: bool,
) -> Option<NodeId> {
    let &[c, a, b] = args else { return None };
    let mut c = rational_vector(cx.graph, c)?;
    let a = rational_matrix(cx.graph, a)?;
    let b = rational_vector(cx.graph, b)?;
    if a.len() != b.len() || a.iter().any(|row| row.len() != c.len()) {
        return None;
    }
    if maximise {
        for v in &mut c {
            *v = -v.clone();
        }
    }
    let relation = if equality { Relation::Equal } else { Relation::AtMost };
    let relations = vec![relation; a.len()];
    Some(match simplex(&c, &a, &b, &relations) {
        | Lp::Optimal(value, point) => {
            let value = if maximise { -value } else { value };
            let value = rational_node(cx.graph, &value);
            let items: Vec<NodeId> = point.iter().map(|p| rational_node(cx.graph, p)).collect();
            let point = cx.graph.node(core::LIST, &items);
            cx.graph.node(core::LIST, &[value, point])
        },
        | Lp::Infeasible => cx.graph.sym("infeasible"),
        | Lp::Unbounded => cx.graph.sym("unbounded"),
    })
}

// ----------------------------------------------------------------------
// Definiteness of rational symmetric matrices
// ----------------------------------------------------------------------

fn rational_determinant(m: &[Vec<Q>]) -> Q {
    let n = m.len();
    let mut work: Vec<Vec<Q>> = m.to_vec();
    let mut det = Q::one();
    for col in 0..n {
        let Some(pivot) = (col..n).find(|&r| !work[r][col].is_zero()) else {
            return Q::zero();
        };
        if pivot != col {
            work.swap(pivot, col);
            det = -det;
        }
        det *= work[col][col].clone();
        for r in col + 1..n {
            let factor = &work[r][col] / &work[col][col];
            for k in col..n {
                let delta = &factor * &work[col][k];
                work[r][k] -= delta;
            }
        }
    }
    det
}

fn principal_minors_nonnegative(m: &[Vec<Q>]) -> bool {
    let n = m.len();
    (1_u32..1 << n).all(|mask| {
        let picked: Vec<usize> = (0..n).filter(|&k| mask >> k & 1 == 1).collect();
        let sub: Vec<Vec<Q>> = picked.iter().map(|&r| picked.iter().map(|&c| m[r][c].clone()).collect()).collect();
        !rational_determinant(&sub).is_negative()
    })
}

fn definiteness(m: &[Vec<Q>]) -> &'static str {
    let n = m.len();
    let negated: Vec<Vec<Q>> = m.iter().map(|r| r.iter().map(|v| -v.clone()).collect()).collect();
    let positive_leading =
        (1..=n).all(|k| rational_determinant(&m[..k].iter().map(|r| r[..k].to_vec()).collect::<Vec<_>>()).is_positive());
    let negative_leading = (1..=n)
        .all(|k| rational_determinant(&negated[..k].iter().map(|r| r[..k].to_vec()).collect::<Vec<_>>()).is_positive());
    if positive_leading {
        "positive_definite"
    } else if negative_leading {
        "negative_definite"
    } else if principal_minors_nonnegative(m) {
        "positive_semidefinite"
    } else if principal_minors_nonnegative(&negated) {
        "negative_semidefinite"
    } else {
        "indefinite"
    }
}

/// The Hessian of `f` with respect to `vars`, simplified.
fn hessian(
    cx: &mut Cx<'_>,
    f: NodeId,
    vars: &[NodeId],
) -> Option<Vec<Vec<NodeId>>> {
    let mut gradient = Vec::with_capacity(vars.len());
    for &v in vars {
        let d = derivative(cx.graph, f, v)?;
        gradient.push(cx.simplify(d));
    }
    let mut h = Vec::with_capacity(vars.len());
    for &g in &gradient {
        let mut row = Vec::with_capacity(vars.len());
        for &v in vars {
            let d = derivative(cx.graph, g, v)?;
            row.push(cx.simplify(d));
        }
        h.push(row);
    }
    Some(h)
}

fn rational_hessian_at(
    cx: &mut Cx<'_>,
    h: &[Vec<NodeId>],
    vars: &[NodeId],
    point: &[NodeId],
) -> Option<Vec<Vec<Q>>> {
    let mut out = Vec::with_capacity(h.len());
    for row in h {
        let mut values = Vec::with_capacity(row.len());
        for &entry in row {
            let at = substitute_point(cx.graph, entry, vars, point);
            let at = cx.simplify(at);
            values.push(cx.graph.number_of(at).and_then(Number::to_rational)?);
        }
        out.push(values);
    }
    Some(out)
}

/// Whether `ge(expression, 0)` is decided true by the engine.
fn provably_nonnegative(
    cx: &mut Cx<'_>,
    expression: NodeId,
) -> bool {
    if let Some(q) = cx.graph.number_of(expression).and_then(Number::to_rational) {
        return !q.is_negative();
    }
    let Some(ge) = cx.graph.ops().lookup("ge") else {
        return false;
    };
    let zero = cx.graph.int(0);
    let Some(request) = cx.graph.try_node(ge, &[expression, zero]) else {
        return false;
    };
    let decided = cx.simplify(request);
    matches!(cx.graph.payload(decided), Some(crate::graph::Payload::Bool(true)))
}

fn truth(
    graph: &mut Graph,
    value: bool,
) -> NodeId {
    graph.lit(crate::graph::Payload::Bool(value))
}

fn convexity(
    cx: &mut Cx<'_>,
    args: &[NodeId],
    concave: bool,
) -> Option<NodeId> {
    let &[f, vars] = args else { return None };
    let f = best(cx.graph, f)?;
    let vars = list_items(cx.graph, vars)?;
    if vars.iter().any(|&v| cx.graph.symbol_of(v).is_none()) || vars.len() > 6 {
        return None;
    }
    let mut h = hessian(cx, f, &vars)?;
    if concave {
        for row in &mut h {
            for entry in row.iter_mut() {
                let minus_one = cx.graph.int(-1);
                let negated = cx.graph.node(core::MUL, &[minus_one, *entry]);
                *entry = cx.simplify(negated);
            }
        }
    }
    let n = vars.len();
    // Exact witnesses against convexity at rational sample points.
    let samples = [-2_i64, -1, 0, 1, 2, 3];
    for code in 0..samples.len().pow(u32::try_from(n).ok()?).min(400) {
        let mut rest = code;
        let point: Vec<NodeId> = (0..n)
            .map(|_| {
                let v = samples[rest % samples.len()];
                rest /= samples.len();
                cx.graph.int(v)
            })
            .collect();
        if let Some(values) = rational_hessian_at(cx, &h, &vars, &point) {
            if !principal_minors_nonnegative(&values) {
                return Some(truth(cx.graph, false));
            }
        }
    }
    // Proof: every principal minor provably non-negative.
    for mask in 1_u32..1 << n {
        let picked: Vec<usize> = (0..n).filter(|&k| mask >> k & 1 == 1).collect();
        let sub: Vec<Vec<NodeId>> =
            picked.iter().map(|&r| picked.iter().map(|&c| h[r][c]).collect()).collect();
        let minor = symbolic_determinant(cx, &sub)?;
        if !provably_nonnegative(cx, minor) {
            return None;
        }
    }
    Some(truth(cx.graph, true))
}

/// Determinant of a small symbolic matrix by cofactor expansion.
fn symbolic_determinant(
    cx: &mut Cx<'_>,
    m: &[Vec<NodeId>],
) -> Option<NodeId> {
    let rows = crate::rules::linalg::matrix_term(cx.graph, m);
    let det = cx.graph.ops().lookup("det")?;
    let request = cx.graph.try_node(det, &[rows])?;
    Some(cx.simplify(request))
}

fn hessian_definiteness(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Option<NodeId> {
    let &[f, vars, point] = args else { return None };
    let f = best(cx.graph, f)?;
    let vars = list_items(cx.graph, vars)?;
    let point = list_items(cx.graph, point)?;
    if vars.len() != point.len() || vars.is_empty() || vars.len() > 8 {
        return None;
    }
    let h = hessian(cx, f, &vars)?;
    let values = rational_hessian_at(cx, &h, &vars, &point)?;
    Some(cx.graph.sym(definiteness(&values)))
}

// ----------------------------------------------------------------------
// KKT
// ----------------------------------------------------------------------

/// `g <= 0` normal form of an inequality.
fn inequality(
    cx: &mut Cx<'_>,
    node: NodeId,
) -> Option<NodeId> {
    let node = best(cx.graph, node)?;
    eprintln!("DEBUG ineq {}", cx.graph.display(node));
    let mut negated = false;
    let mut current = node;
    let mut name = cx.graph.ops().get(cx.graph.op(current)).name.to_string();
    while name == "not" {
        let &[inner] = cx.graph.children(current) else { return None };
        current = inner;
        negated = !negated;
        name = cx.graph.ops().get(cx.graph.op(current)).name.to_string();
    }
    let kids = cx.graph.children(current).to_vec();
    let minus = |cx: &mut Cx<'_>, a: NodeId, b: NodeId| {
        let minus_one = cx.graph.int(-1);
        let negated = cx.graph.node(core::MUL, &[minus_one, b]);
        let sum = cx.graph.node(core::ADD, &[a, negated]);
        cx.simplify(sum)
    };
    // `lower` means `a <= b` (so `a - b <= 0`); `not` swaps the sides.
    let (lower, a, b) = match (name.as_str(), kids.as_slice()) {
        | ("le" | "lt", &[a, b]) => (true, a, b),
        | ("ge" | "gt", &[a, b]) => (false, a, b),
        | _ if negated => return None,
        | _ => return Some(node),
    };
    let lower = lower != negated;
    Some(if lower { minus(cx, a, b) } else { minus(cx, b, a) })
}

fn subsets(
    m: usize,
    max_size: usize,
) -> Vec<Vec<usize>> {
    let mut out: Vec<Vec<usize>> = (0_u32..1 << m)
        .map(|mask| (0..m).filter(|&k| mask >> k & 1 == 1).collect::<Vec<_>>())
        .filter(|s| s.len() <= max_size)
        .collect();
    out.sort_by_key(Vec::len);
    out
}

struct KktPoint {
    point: Vec<NodeId>,
    multipliers: Vec<NodeId>,
    value: NodeId,
    numeric: Option<f64>,
}

fn kkt(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Option<Vec<KktPoint>> {
    let &[f, constraints, vars] = args else { return None };
    let f = best(cx.graph, f)?;
    let vars = list_items(cx.graph, vars)?;
    let mut gs = Vec::new();
    for g in list_items(cx.graph, constraints)? {
        gs.push(inequality(cx, g)?);
    }
    if vars.iter().any(|&v| cx.graph.symbol_of(v).is_none()) || gs.len() > MAX_CONSTRAINTS {
        return None;
    }
    let mut grad_f = Vec::new();
    for &v in &vars {
        let d = derivative(cx.graph, f, v)?;
        grad_f.push(cx.simplify(d));
    }
    let mut grad_g = Vec::new();
    for &g in &gs {
        let mut row = Vec::new();
        for &v in &vars {
            let d = derivative(cx.graph, g, v)?;
            row.push(cx.simplify(d));
        }
        grad_g.push(row);
    }
    let mut found: Vec<KktPoint> = Vec::new();
    let mut seen: Vec<String> = Vec::new();
    for active in subsets(gs.len(), vars.len()) {
        let mus: Vec<NodeId> = active
            .iter()
            .map(|&k| {
                let s = cx.graph.interner_mut().fresh_symbol(&format!("mu{}", k + 1));
                cx.graph.symbol_node(s)
            })
            .collect();
        let mut equations = Vec::new();
        for (j, &gf) in grad_f.iter().enumerate() {
            let mut terms = vec![gf];
            for (&mu, &k) in mus.iter().zip(&active) {
                terms.push(cx.graph.node(core::MUL, &[mu, grad_g[k][j]]));
            }
            let sum = cx.graph.node(core::ADD, &terms);
            equations.push(cx.simplify(sum));
        }
        for &k in &active {
            equations.push(gs[k]);
        }
        let mut unknowns = vars.clone();
        unknowns.extend(&mus);
        let Some(solutions) = solve_system(cx, &equations, &unknowns) else {
            continue;
        };
        for solution in solutions {
            let (point, active_mus) = solution.split_at(vars.len());
            // Multipliers must be non-negative.
            let mut ok = true;
            for &mu in active_mus {
                match number(cx, mu) {
                    | Some(v) if v >= -1e-9 => {},
                    | _ => ok = false,
                }
            }
            if !ok {
                continue;
            }
            // Feasibility of the inactive constraints.
            for (k, &g) in gs.iter().enumerate() {
                if active.contains(&k) {
                    continue;
                }
                let at = substitute_point(cx.graph, g, &vars, point);
                match number(cx, at) {
                    | Some(v) if v <= 1e-9 => {},
                    | _ => ok = false,
                }
            }
            if !ok {
                continue;
            }
            let key = point.iter().map(|&p| cx.graph.display(p)).collect::<Vec<_>>().join(",");
            if seen.contains(&key) {
                continue;
            }
            seen.push(key);
            let zero = cx.graph.int(0);
            let mut multipliers = vec![zero; gs.len()];
            for (&k, &mu) in active.iter().zip(active_mus) {
                multipliers[k] = mu;
            }
            let value = substitute_point(cx.graph, f, &vars, point);
            let value = cx.simplify(value);
            let numeric = number(cx, value);
            found.push(KktPoint { point: point.to_vec(), multipliers, value, numeric });
        }
    }
    found.sort_by(|a, b| match (a.numeric, b.numeric) {
        | (Some(x), Some(y)) => x.total_cmp(&y),
        | _ => std::cmp::Ordering::Equal,
    });
    Some(found)
}

fn kkt_list(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Option<NodeId> {
    let points = kkt(cx, args)?;
    let mut items = Vec::with_capacity(points.len());
    for p in points {
        let point = cx.graph.node(core::LIST, &p.point);
        let mus = cx.graph.node(core::LIST, &p.multipliers);
        items.push(cx.graph.node(core::LIST, &[point, mus, p.value]));
    }
    Some(cx.graph.node(core::LIST, &items))
}

fn kkt_minimum(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Option<NodeId> {
    let points = kkt(cx, args)?;
    if points.iter().any(|p| p.numeric.is_none()) {
        return None;
    }
    let first = points.into_iter().next()?;
    let point = cx.graph.node(core::LIST, &first.point);
    Some(cx.graph.node(core::LIST, &[point, first.value]))
}

// ----------------------------------------------------------------------
// Quadratic programming and Legendre transform
// ----------------------------------------------------------------------

fn request(
    graph: &mut Graph,
    name: &str,
    args: &[NodeId],
) -> Option<NodeId> {
    let op = graph.ops().lookup(name)?;
    graph.try_node(op, args)
}

fn qp_eq(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Option<NodeId> {
    let &[q, c, a, b] = args else { return None };
    let (q, c, a, b) = (matrix_of(cx.graph, q)?, vector_of(cx.graph, c)?, matrix_of(cx.graph, a)?, vector_of(cx.graph, b)?);
    let (n, m) = (c.len(), b.len());
    if q.len() != n || q.iter().any(|r| r.len() != n) || a.len() != m || a.iter().any(|r| r.len() != n) {
        return None;
    }
    let zero = cx.graph.int(0);
    let mut kkt_rows = Vec::with_capacity(n + m);
    for (i, row) in q.iter().enumerate() {
        let mut full = row.clone();
        full.extend(a.iter().map(|ar| ar[i]));
        kkt_rows.push(full);
    }
    for row in &a {
        let mut full = row.clone();
        full.extend(std::iter::repeat_n(zero, m));
        kkt_rows.push(full);
    }
    let mut rhs = Vec::with_capacity(n + m);
    for &ci in &c {
        let minus_one = cx.graph.int(-1);
        rhs.push(cx.graph.node(core::MUL, &[minus_one, ci]));
    }
    rhs.extend(&b);
    let matrix = crate::rules::linalg::matrix_term(cx.graph, &kkt_rows);
    let rhs = cx.graph.node(core::LIST, &rhs);
    let system = request(cx.graph, "linsolve", &[matrix, rhs])?;
    let solved = cx.simplify(system);
    if cx.graph.op(solved) != core::LIST {
        return None;
    }
    let items = cx.graph.children(solved).to_vec();
    if items.len() != n + m {
        return None;
    }
    let (x, lambda) = items.split_at(n);
    let (x, lambda) = (cx.graph.node(core::LIST, x), cx.graph.node(core::LIST, lambda));
    Some(cx.graph.node(core::LIST, &[x, lambda]))
}

fn vector_of(
    graph: &mut Graph,
    node: NodeId,
) -> Option<Vec<NodeId>> {
    let node = best(graph, node)?;
    (graph.op(node) == core::LIST).then(|| graph.children(node).to_vec())
}

fn matrix_of(
    graph: &mut Graph,
    node: NodeId,
) -> Option<Vec<Vec<NodeId>>> {
    let rows = vector_of(graph, node)?;
    rows.into_iter().map(|r| vector_of(graph, r)).collect()
}

fn legendre_transform(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Option<NodeId> {
    let &[f, x, p] = args else { return None };
    let f = best(cx.graph, f)?;
    let derivative_of_f = derivative(cx.graph, f, x)?;
    let equation = {
        let minus_one = cx.graph.int(-1);
        let negated = cx.graph.node(core::MUL, &[minus_one, p]);
        let sum = cx.graph.node(core::ADD, &[derivative_of_f, negated]);
        cx.simplify(sum)
    };
    let roots = solve_for(cx.graph, equation, x, 0)?;
    let &[root] = roots.as_slice() else { return None };
    let at = cx.graph.substitute(f, x, root);
    let minus_one = cx.graph.int(-1);
    let negated = cx.graph.node(core::MUL, &[minus_one, at]);
    let px = cx.graph.node(core::MUL, &[p, root]);
    let sum = cx.graph.node(core::ADD, &[px, negated]);
    Some(cx.simplify(sum))
}

// ----------------------------------------------------------------------
// Kernel
// ----------------------------------------------------------------------

#[derive(Copy, Clone)]
enum Kind {
    Linprog { maximise: bool, equality: bool },
    KktPoints,
    KktMinimum,
    Convex { concave: bool },
    Definiteness,
    QpEq,
    Legendre,
}

struct Extension {
    op: OpId,
    kind: Kind,
}

impl Kernel for Extension {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let args = cx.graph.children(node).to_vec();
        let result = match self.kind {
            | Kind::Linprog { maximise, equality } => linprog(cx, &args, maximise, equality),
            | Kind::KktPoints => kkt_list(cx, &args),
            | Kind::KktMinimum => kkt_minimum(cx, &args),
            | Kind::Convex { concave } => convexity(cx, &args, concave),
            | Kind::Definiteness => hessian_definiteness(cx, &args),
            | Kind::QpEq => qp_eq(cx, &args),
            | Kind::Legendre => legendre_transform(cx, &args),
        };
        result.map_or(Outcome::Pass, Outcome::Equal)
    }
}

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    for (name, arity, kind) in [
        ("linprog", 3, Kind::Linprog { maximise: false, equality: false }),
        ("linprog_max", 3, Kind::Linprog { maximise: true, equality: false }),
        ("linprog_eq", 3, Kind::Linprog { maximise: false, equality: true }),
        ("kkt_points", 3, Kind::KktPoints),
        ("kkt_minimum", 3, Kind::KktMinimum),
        ("is_convex", 2, Kind::Convex { concave: false }),
        ("is_concave", 2, Kind::Convex { concave: true }),
        ("hessian_definiteness", 3, Kind::Definiteness),
        ("qp_eq", 4, Kind::QpEq),
        ("legendre_transform", 3, Kind::Legendre),
    ] {
        let op = i.op(OpDescriptor::new(name, Arity::Fixed(arity)).flags(OpFlags::HEAVY).cost(100))?;
        i.kernel(&format!("optimize/{name}"), Tier::Reduce, Extension { op, kind });
    }
    Ok(())
}

/// Sign of the constrained second variation on the tangent space of the
/// constraints: Hessian of the Lagrangian projected onto the null space of
/// the constraint Jacobian, classified by its leading minors.
pub(super) fn projected_kind(
    cx: &mut Cx<'_>,
    lagrangian: NodeId,
    constraints: &[NodeId],
    all: &[NodeId],
    solution: &[NodeId],
) -> Option<&'static str> {
    let n = all.len().checked_sub(constraints.len())?;
    let vars = all.get(..n)?;
    let at = |cx: &mut Cx<'_>, term: NodeId| -> Option<f64> {
        let point = substitute_point(cx.graph, term, all, solution);
        number(cx, point)
    };
    let mut jacobian = Vec::with_capacity(constraints.len());
    for &g in constraints {
        let mut row = Vec::with_capacity(n);
        for &v in vars {
            let d = derivative(cx.graph, g, v)?;
            row.push(at(cx, d)?);
        }
        jacobian.push(row);
    }
    let mut h = Vec::with_capacity(n);
    for &vi in vars {
        let li = derivative(cx.graph, lagrangian, vi)?;
        let mut row = Vec::with_capacity(n);
        for &vj in vars {
            let d = derivative(cx.graph, li, vj)?;
            row.push(at(cx, d)?);
        }
        h.push(row);
    }
    let basis = null_basis(&jacobian, n);
    if basis.is_empty() {
        return Some("isolated");
    }
    let d = basis.len();
    let mut reduced = vec![vec![0.0; d]; d];
    for (a, za) in basis.iter().enumerate() {
        for (b, zb) in basis.iter().enumerate() {
            let mut total = 0.0;
            for (i, hi) in h.iter().enumerate() {
                for (j, &hij) in hi.iter().enumerate() {
                    total += za[i] * hij * zb[j];
                }
            }
            reduced[a][b] = total;
        }
    }
    let minors = leading_minors(&reduced);
    let tiny = |v: f64| v.abs() < 1e-9;
    if minors.iter().any(|&m| tiny(m)) {
        return None;
    }
    if minors.iter().all(|&m| m > 0.0) {
        Some("local_min")
    } else if minors.iter().enumerate().all(|(k, &m)| (m < 0.0) == (k % 2 == 0)) {
        Some("local_max")
    } else {
        Some("saddle")
    }
}

/// An orthonormal basis of the null space of the rows of `jacobian`.
fn null_basis(
    jacobian: &[Vec<f64>],
    n: usize,
) -> Vec<Vec<f64>> {
    // Row reduce, then read off the free columns.
    let mut m: Vec<Vec<f64>> = jacobian.to_vec();
    let mut pivots = Vec::new();
    let mut row = 0;
    for col in 0..n {
        if row >= m.len() {
            break;
        }
        let Some(p) = (row..m.len()).max_by(|&a, &b| m[a][col].abs().total_cmp(&m[b][col].abs())) else {
            break;
        };
        if m[p][col].abs() < 1e-10 {
            continue;
        }
        m.swap(p, row);
        let pivot = m[row][col];
        for value in &mut m[row] {
            *value /= pivot;
        }
        let pivot_row = m[row].clone();
        for (r, other) in m.iter_mut().enumerate() {
            if r != row {
                let factor = other[col];
                for (value, p) in other.iter_mut().zip(&pivot_row) {
                    *value -= factor * p;
                }
            }
        }
        pivots.push(col);
        row += 1;
    }
    let mut basis: Vec<Vec<f64>> = Vec::new();
    for free in (0..n).filter(|c| !pivots.contains(c)) {
        let mut v = vec![0.0; n];
        v[free] = 1.0;
        for (r, &p) in pivots.iter().enumerate() {
            v[p] = -m[r][free];
        }
        // Gram–Schmidt against the previous basis vectors.
        for b in &basis {
            let dot: f64 = v.iter().zip(b).map(|(x, y)| x * y).sum();
            for (x, y) in v.iter_mut().zip(b) {
                *x -= dot * y;
            }
        }
        let length = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        if length > 1e-12 {
            basis.push(v.into_iter().map(|x| x / length).collect());
        }
    }
    basis
}

#[cfg(test)]
mod tests {
    use crate::rules::optimize::optimize;
    use crate::rules::testing::simplify;

    fn run(src: &str) -> String {
        simplify(&[optimize()], src)
    }

    #[test]
    fn exact_linear_programs() {
        // max 3x + 2y, x + y <= 4, x + 3y <= 6, x, y >= 0: optimum 12 at (4, 0).
        assert_eq!(run("linprog_max(list(3, 2), list(list(1, 1), list(1, 3)), list(4, 6))"), "list(12, list(4, 0))");
        // min -x - y with the same constraints: -4.
        assert_eq!(run("linprog(list(-1, -1), list(list(1, 1), list(1, 3)), list(4, 6))"), "list(-4, list(4, 0))");
        // Fractional optimum.
        assert_eq!(run("linprog_max(list(1, 1), list(list(2, 1), list(1, 3)), list(4, 5))"), "list(13/5, list(7/5, 6/5))");
        assert_eq!(run("linprog_eq(list(1, 2), list(list(1, 1)), list(3))"), "list(3, list(3, 0))");
        // x >= 2 and x <= 1 is infeasible; -x is unbounded.
        assert_eq!(run("linprog(list(1), list(list(-1), list(1)), list(-2, 1))"), "infeasible");
        assert_eq!(run("linprog(list(-1), list(list(-1)), list(0))"), "unbounded");
    }

    #[test]
    fn kkt_conditions() {
        // min x^2 + y^2 subject to x + y >= 2: the point (1, 1) with multiplier 2.
        assert_eq!(
            run("kkt_points(x^2 + y^2, list(ge(2, x + y)), list(x, y))"),
            "list(list(list(0, 0), list(0), 0))"
        );
        assert_eq!(
            run("kkt_points(x^2 + y^2, list(le(2, x + y)), list(x, y))"),
            "list(list(list(1, 1), list(2), 2))"
        );
        assert_eq!(run("kkt_minimum(x^2 + y^2, list(le(2, x + y)), list(x, y))"), "list(list(1, 1), 2)");
        // min (x - 2)^2 subject to x <= 1: boundary point.
        assert_eq!(run("kkt_minimum((x - 2)^2, list(x - 1), list(x))"), "list(list(1), 1)");
        // Inactive constraint: the unconstrained minimum survives.
        assert_eq!(run("kkt_minimum((x - 2)^2, list(x - 5), list(x))"), "list(list(2), 0)");
    }

    #[test]
    fn convexity_and_definiteness() {
        assert_eq!(run("is_convex(x^2 + y^2 + x*y, list(x, y))"), "true");
        assert_eq!(run("is_convex(x^2 - y^2, list(x, y))"), "false");
        assert_eq!(run("is_concave(-x^2 - y^2, list(x, y))"), "true");
        assert_eq!(run("is_convex(x^3, list(x))"), "false");
        assert_eq!(run("hessian_definiteness(x^2 + y^2, list(x, y), list(0, 0))"), "positive_definite");
        assert_eq!(run("hessian_definiteness(x*y, list(x, y), list(0, 0))"), "indefinite");
        assert_eq!(run("hessian_definiteness(x^2, list(x, y), list(1, 1))"), "positive_semidefinite");
        assert_eq!(run("hessian_definiteness(-x^2 - y^2, list(x, y), list(1, 1))"), "negative_definite");
    }

    #[test]
    fn quadratic_programs_and_transforms() {
        // min 1/2 (x^2 + y^2) subject to x + y = 2: (1, 1) with multiplier -1.
        assert_eq!(
            run("qp_eq(list(list(1, 0), list(0, 1)), list(0, 0), list(list(1, 1)), list(2))"),
            "list(list(1, 1), list(-1))"
        );
        assert_eq!(run("legendre_transform(x^2/2, x, p)"), "1/2*p^2");
        assert_eq!(run("legendre_transform(x^2, x, p)"), "1/4*p^2");
    }

    #[test]
    fn constrained_extrema_in_three_variables() {
        // Extremes of x + y + z on the unit sphere: +-sqrt(3), classified by the projected Hessian.
        let text = run("find_constrained_extrema(x + y + z, list(x^2 + y^2 + z^2 - 1), list(x, y, z))");
        assert!(text.contains("local_max") && text.contains("local_min"), "{text}");
    }
}
