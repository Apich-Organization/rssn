//! Linear algebra and vector calculus.
//!
//! A vector is `list(a, b, ...)`; a matrix is a list of rows,
//! `list(list(a, b), list(c, d))`. Matrix products are their own operator,
//! `matmul`, because `*` is commutative.
//!
//! Exact matrices (every entry a rational number) are handled with exact
//! rational arithmetic; symbolic matrices with terms, simplifying and
//! zero-testing pivots through the engine. Numeric decompositions that have
//! no exact counterpart (singular values, eigenvectors of symmetric float
//! matrices) use the numeric kernels.
//!
//! | operator | value |
//! |---|---|
//! | `matmul(A, B)`, `madd(A, B)`, `smul(c, A)`, `transpose(A)`, `identity(n)`, `zeros(m, n)` | matrix arithmetic |
//! | `dims(A)` | `list(rows, columns)` |
//! | `det(A)`, `trace(A)`, `rank(A)`, `inverse(A)`, `rref(A)`, `nullspace(A)` | the usual |
//! | `linsolve(A, b)` | the solution of `A x = b` (`list()` if inconsistent; free unknowns become parameters `t1, t2, ...`) |
//! | `charpoly(A, l)`, `eigenvals(A)`, `eigenvects(A)` | characteristic polynomial; eigenvalues (with multiplicity); `list(list(value, multiplicity, list(basis...)), ...)` |
//! | `lu(A)`, `qr(A)`, `svd(A)` | `list(P, L, U)` with `P A = L U`; `list(Q, R)`; numeric `list(U, S, V)` |
//! | `dot`, `cross`, `norm`, `normalize`, `angle`, `project`, `outer` | vector operations |
//! | `grad(f, vars)`, `div(F, vars)`, `curl(F, vars)`, `laplacian(f, vars)`, `jacobian(F, vars)`, `hessian(f, vars)`, `directional(f, vars, v)` | vector calculus |
//! | `line_integral(f, curve, t, a, b)`, `line_integral_vec(F, curve, t, a, b)`, `surface_integral(f, surface, u, v, ua, ub, va, vb)`, `volume_integral(f, vars, bounds)` | integrals over parametrised curves, surfaces and boxes, reduced to definite integrals |

use num_rational::BigRational;
use num_traits::One;
use num_traits::Zero;

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
use crate::graph::RuleSet;
use crate::graph::Tier;
use crate::kernels::matrix::FaerDecompositionResult;
use crate::kernels::matrix::FaerDecompositionType;
use crate::kernels::matrix::Matrix;
use crate::rules::poly::best;
use crate::rules::poly::repr::from_term;
use crate::rules::poly::repr::to_term;
use crate::rules::poly::repr::Gens;
use crate::rules::poly::repr::Limits;
use crate::rules::poly::repr::Poly;

use super::calculus::calculus;
use super::calculus::derivative;
use super::solve::solve;
use super::solve::solve_for;

/// The linear algebra rule set.
#[must_use]
pub fn linalg() -> RuleSet {
    RuleSet::new("linalg", install).needs(calculus()).needs(solve())
}

/// Which operation a [`LinalgKernel`] performs.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Request {
    Matmul,
    Madd,
    Smul,
    Transpose,
    Identity,
    Zeros,
    Dims,
    Det,
    Trace,
    Rank,
    Inverse,
    Rref,
    Nullspace,
    Linsolve,
    Charpoly,
    Eigenvals,
    Eigenvects,
    Lu,
    Qr,
    Svd,
    Dot,
    Cross,
    Norm,
    Normalize,
    Angle,
    Project,
    Outer,
    Grad,
    Div,
    Curl,
    Laplacian,
    Jacobian,
    Hessian,
    Directional,
    LineIntegral,
    LineIntegralVec,
    SurfaceIntegral,
    VolumeIntegral,
}

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let table: [(&str, Arity, Request); 38] = [
        ("matmul", Arity::Variadic, Request::Matmul),
        ("madd", Arity::Variadic, Request::Madd),
        ("smul", Arity::Fixed(2), Request::Smul),
        ("transpose", Arity::Fixed(1), Request::Transpose),
        ("identity", Arity::Fixed(1), Request::Identity),
        ("zeros", Arity::Fixed(2), Request::Zeros),
        ("dims", Arity::Fixed(1), Request::Dims),
        ("det", Arity::Fixed(1), Request::Det),
        ("trace", Arity::Fixed(1), Request::Trace),
        ("rank", Arity::Fixed(1), Request::Rank),
        ("inverse", Arity::Fixed(1), Request::Inverse),
        ("rref", Arity::Fixed(1), Request::Rref),
        ("nullspace", Arity::Fixed(1), Request::Nullspace),
        ("linsolve", Arity::Fixed(2), Request::Linsolve),
        ("charpoly", Arity::Fixed(2), Request::Charpoly),
        ("eigenvals", Arity::Fixed(1), Request::Eigenvals),
        ("eigenvects", Arity::Fixed(1), Request::Eigenvects),
        ("lu", Arity::Fixed(1), Request::Lu),
        ("qr", Arity::Fixed(1), Request::Qr),
        ("svd", Arity::Fixed(1), Request::Svd),
        ("dot", Arity::Fixed(2), Request::Dot),
        ("cross", Arity::Fixed(2), Request::Cross),
        ("norm", Arity::Fixed(1), Request::Norm),
        ("normalize", Arity::Fixed(1), Request::Normalize),
        ("angle", Arity::Fixed(2), Request::Angle),
        ("project", Arity::Fixed(2), Request::Project),
        ("outer", Arity::Fixed(2), Request::Outer),
        ("grad", Arity::Fixed(2), Request::Grad),
        ("div", Arity::Fixed(2), Request::Div),
        ("curl", Arity::Fixed(2), Request::Curl),
        ("laplacian", Arity::Fixed(2), Request::Laplacian),
        ("jacobian", Arity::Fixed(2), Request::Jacobian),
        ("hessian", Arity::Fixed(2), Request::Hessian),
        ("directional", Arity::Fixed(3), Request::Directional),
        ("line_integral", Arity::Fixed(5), Request::LineIntegral),
        ("line_integral_vec", Arity::Fixed(5), Request::LineIntegralVec),
        ("surface_integral", Arity::Fixed(8), Request::SurfaceIntegral),
        ("volume_integral", Arity::Fixed(3), Request::VolumeIntegral),
    ];
    for (name, arity, request) in table {
        let op = i.op(OpDescriptor::new(name, arity).flags(OpFlags::HEAVY).cost(100))?;
        i.kernel(&format!("linalg/{name}"), Tier::Reduce, LinalgKernel { op, request });
    }
    i.kernel("linalg/broadcast", Tier::Normalize, Broadcast);
    Ok(())
}

/// Arithmetic on lists is elementwise: a sum of lists of one length is
/// the list of sums, and a product with exactly one list factor scales
/// every element. Nested lists (matrices, tensors) follow, one level per
/// application.
struct Broadcast;

impl Broadcast {
    fn as_list(
        graph: &Graph,
        node: NodeId,
    ) -> Option<NodeId> {
        graph.enodes(graph.find(node)).find(|&n| graph.op(n) == core::LIST)
    }
}

impl Kernel for Broadcast {
    fn ops(&self) -> Vec<OpId> {
        vec![core::ADD, core::MUL]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let graph = &mut *cx.graph;
        let children = graph.children(node).to_vec();
        let lists: Vec<Option<NodeId>> = children.iter().map(|&c| Self::as_list(graph, c)).collect();
        if graph.op(node) == core::ADD {
            let Some(lists) = lists.into_iter().collect::<Option<Vec<NodeId>>>() else {
                return Outcome::Pass;
            };
            let len = graph.children(lists[0]).len();
            if lists.iter().any(|&l| graph.children(l).len() != len) {
                return Outcome::Pass;
            }
            let items: Vec<NodeId> = (0..len)
                .map(|k| {
                    let terms: Vec<NodeId> = lists.iter().map(|&l| graph.children(l)[k]).collect();
                    add(graph, &terms)
                })
                .collect();
            return Outcome::Equal(list(graph, &items));
        }
        let mut found = lists.iter().enumerate().filter_map(|(k, l)| l.map(|l| (k, l)));
        let (Some((at, l)), None) = (found.next(), found.next()) else {
            return Outcome::Pass;
        };
        let others: Vec<NodeId> = children.iter().enumerate().filter(|&(k, _)| k != at).map(|(_, &c)| c).collect();
        let items: Vec<NodeId> = graph
            .children(l)
            .to_vec()
            .into_iter()
            .map(|item| {
                let mut factors = others.clone();
                factors.push(item);
                mul(graph, &factors)
            })
            .collect();
        Outcome::Equal(list(graph, &items))
    }

    fn revisit(&self) -> bool {
        true
    }
}

// ----------------------------------------------------------------------
// Reading and writing matrices
// ----------------------------------------------------------------------

/// Entries of a vector `list(...)`, after taking the best form of each.
fn vector(
    graph: &mut Graph,
    node: NodeId,
) -> Option<Vec<NodeId>> {
    let node = best(graph, node)?;
    if graph.op(node) != core::LIST {
        return None;
    }
    Some(graph.children(node).to_vec())
}

/// Rows of a matrix `list(list(...), ...)`; all rows must have the same
/// length.
fn matrix(
    graph: &mut Graph,
    node: NodeId,
) -> Option<Vec<Vec<NodeId>>> {
    let rows = vector(graph, node)?;
    let mut out = Vec::with_capacity(rows.len());
    for row in rows {
        if graph.op(row) != core::LIST {
            return None;
        }
        out.push(graph.children(row).to_vec());
    }
    let width = out.first().map_or(0, Vec::len);
    out.iter().all(|r| r.len() == width).then_some(out)
}

/// A matrix or a vector operand.
enum Operand {
    Matrix(Vec<Vec<NodeId>>),
    Vector(Vec<NodeId>),
}

/// Reads a list of lists as a matrix and a list of scalars as a vector.
fn operand(
    graph: &mut Graph,
    node: NodeId,
) -> Option<Operand> {
    if let Some(m) = matrix(graph, node) {
        if !m.is_empty() {
            return Some(Operand::Matrix(m));
        }
    }
    let items = vector(graph, node)?;
    let mut scalars = Vec::with_capacity(items.len());
    for item in items {
        let item = best(graph, item)?;
        if graph.op(item) == core::LIST {
            return None;
        }
        scalars.push(item);
    }
    Some(Operand::Vector(scalars))
}

fn list(
    graph: &mut Graph,
    items: &[NodeId],
) -> NodeId {
    graph.node(core::LIST, items)
}

fn matrix_term(
    graph: &mut Graph,
    rows: &[Vec<NodeId>],
) -> NodeId {
    let row_nodes: Vec<NodeId> = rows.iter().map(|r| list(graph, r)).collect();
    list(graph, &row_nodes)
}

fn mul(
    graph: &mut Graph,
    factors: &[NodeId],
) -> NodeId {
    match factors {
        | [] => graph.int(1),
        | [only] => *only,
        | _ => graph.node(core::MUL, factors),
    }
}

fn add(
    graph: &mut Graph,
    terms: &[NodeId],
) -> NodeId {
    match terms {
        | [] => graph.int(0),
        | [only] => *only,
        | _ => graph.node(core::ADD, terms),
    }
}

fn neg(
    graph: &mut Graph,
    node: NodeId,
) -> NodeId {
    let minus_one = graph.int(-1);
    mul(graph, &[minus_one, node])
}

fn reciprocal(
    graph: &mut Graph,
    node: NodeId,
) -> NodeId {
    let minus_one = graph.int(-1);
    graph.node(core::POW, &[node, minus_one])
}

fn sqrt(
    graph: &mut Graph,
    node: NodeId,
) -> Option<NodeId> {
    let half = graph.num(Number::fraction(1, 2)?);
    Some(graph.node(core::POW, &[node, half]))
}

/// Whether the best form of `node` divides by something.
fn has_division(
    graph: &mut Graph,
    node: NodeId,
) -> bool {
    let Some(term) = best(graph, node) else {
        return false;
    };
    let mut stack = vec![term];
    let mut seen = std::collections::HashSet::new();
    while let Some(n) = stack.pop() {
        if !seen.insert(n) {
            continue;
        }
        if graph.op(n) == core::POW {
            if let Some(e) = graph.children(n).get(1).and_then(|&e| graph.number_of(e)) {
                if e.to_f64() < 0.0 {
                    return true;
                }
            }
        }
        stack.extend_from_slice(graph.children(n));
    }
    false
}

/// Brings entries that divide over one denominator; other entries are
/// left for the engine.
fn normalise_entries(
    cx: &mut Cx<'_>,
    m: &mut [Vec<NodeId>],
) {
    for row in m {
        for entry in row.iter_mut() {
            if has_division(cx.graph, *entry) {
                *entry = normal(cx, *entry);
            }
        }
    }
}

/// Simplifies `node`, bringing rational functions over one denominator
/// (with common univariate factors cancelled) so that entries of symbolic
/// matrices stay in a normal form and zero tests are reliable.
fn normal(
    cx: &mut Cx<'_>,
    node: NodeId,
) -> NodeId {
    let simplified = cx.simplify(node);
    if !has_division(cx.graph, simplified) {
        return simplified;
    }
    match request(cx.graph, "cancel", &[simplified]) {
        | Some(cancel) => {
            let together = cx.simplify(cancel);
            if cx.graph.op(together) == cx.graph.op(cancel) { simplified } else { together }
        },
        | None => simplified,
    }
}

// ----------------------------------------------------------------------
// Fields: exact rationals and simplified terms
// ----------------------------------------------------------------------

/// The arithmetic Gaussian elimination needs.
trait Field {
    type E: Clone;
    fn zero(&mut self) -> Self::E;
    fn one(&mut self) -> Self::E;
    fn sub(
        &mut self,
        a: &Self::E,
        b: &Self::E,
    ) -> Self::E;
    fn mul(
        &mut self,
        a: &Self::E,
        b: &Self::E,
    ) -> Self::E;
    fn div(
        &mut self,
        a: &Self::E,
        b: &Self::E,
    ) -> Option<Self::E>;
    fn is_zero(
        &mut self,
        a: &Self::E,
    ) -> bool;
}

struct Rationals;

impl Field for Rationals {
    type E = BigRational;

    fn zero(&mut self) -> BigRational {
        BigRational::zero()
    }

    fn one(&mut self) -> BigRational {
        BigRational::one()
    }

    fn sub(
        &mut self,
        a: &BigRational,
        b: &BigRational,
    ) -> BigRational {
        a - b
    }

    fn mul(
        &mut self,
        a: &BigRational,
        b: &BigRational,
    ) -> BigRational {
        a * b
    }

    fn div(
        &mut self,
        a: &BigRational,
        b: &BigRational,
    ) -> Option<BigRational> {
        (!b.is_zero()).then(|| a / b)
    }

    fn is_zero(
        &mut self,
        a: &BigRational,
    ) -> bool {
        a.is_zero()
    }
}

struct Terms<'c, 'a> {
    cx: &'c mut Cx<'a>,
}

impl Field for Terms<'_, '_> {
    type E = NodeId;

    fn zero(&mut self) -> NodeId {
        self.cx.graph.int(0)
    }

    fn one(&mut self) -> NodeId {
        self.cx.graph.int(1)
    }

    fn sub(
        &mut self,
        a: &NodeId,
        b: &NodeId,
    ) -> NodeId {
        let negated = neg(self.cx.graph, *b);
        let sum = add(self.cx.graph, &[*a, negated]);
        normal(self.cx, sum)
    }

    fn mul(
        &mut self,
        a: &NodeId,
        b: &NodeId,
    ) -> NodeId {
        let product = mul(self.cx.graph, &[*a, *b]);
        normal(self.cx, product)
    }

    fn div(
        &mut self,
        a: &NodeId,
        b: &NodeId,
    ) -> Option<NodeId> {
        if self.is_zero(b) {
            return None;
        }
        let inverse = reciprocal(self.cx.graph, *b);
        let quotient = mul(self.cx.graph, &[*a, inverse]);
        Some(normal(self.cx, quotient))
    }

    fn is_zero(
        &mut self,
        a: &NodeId,
    ) -> bool {
        if let Some(n) = self.cx.graph.number_of(*a) {
            return n.is_zero();
        }
        let value = normal(self.cx, *a);
        self.cx.graph.number_of(value).is_some_and(Number::is_zero)
    }
}

/// Reduced row echelon form in place; returns the pivot columns and the
/// determinant factor (product of pivots, with the sign of the row swaps),
/// which is the determinant when the matrix is square and of full rank.
fn rref<F: Field>(
    field: &mut F,
    m: &mut [Vec<F::E>],
) -> Option<(Vec<usize>, F::E)> {
    let rows = m.len();
    let cols = m.first().map_or(0, Vec::len);
    let mut pivots = Vec::new();
    let mut factor = field.one();
    let mut row = 0;
    for col in 0..cols {
        if row >= rows {
            break;
        }
        let Some(pivot_row) = (row..rows).find(|&r| !field.is_zero(&m[r][col])) else {
            continue;
        };
        if pivot_row != row {
            m.swap(pivot_row, row);
            let minus_one = field.zero();
            let one = field.one();
            let negative = field.sub(&minus_one, &one);
            factor = field.mul(&factor, &negative);
        }
        let pivot = m[row][col].clone();
        factor = field.mul(&factor, &pivot);
        for c in 0..cols {
            m[row][c] = field.div(&m[row][c], &pivot)?;
        }
        for r in 0..rows {
            if r == row || field.is_zero(&m[r][col]) {
                continue;
            }
            let scale = m[r][col].clone();
            for c in 0..cols {
                let delta = field.mul(&scale, &m[row][c]);
                m[r][c] = field.sub(&m[r][c], &delta);
            }
        }
        pivots.push(col);
        row += 1;
    }
    Some((pivots, factor))
}

/// The entries as exact rationals, if all of them are.
fn rational_matrix(
    graph: &Graph,
    m: &[Vec<NodeId>],
) -> Option<Vec<Vec<BigRational>>> {
    m.iter().map(|row| row.iter().map(|&e| graph.number_of(e).and_then(Number::to_rational)).collect()).collect()
}

fn rational_node(
    graph: &mut Graph,
    value: &BigRational,
) -> NodeId {
    graph.num(Number::rat(value.clone()))
}

/// Runs `rref` over the best field for the matrix and returns the reduced
/// matrix as terms, the pivots, and the determinant factor.
fn reduce(
    cx: &mut Cx<'_>,
    m: &[Vec<NodeId>],
) -> Option<(Vec<Vec<NodeId>>, Vec<usize>, NodeId)> {
    if let Some(mut exact) = rational_matrix(cx.graph, m) {
        let (pivots, factor) = rref(&mut Rationals, &mut exact)?;
        let reduced = exact.iter().map(|row| row.iter().map(|v| rational_node(cx.graph, v)).collect()).collect();
        let factor = rational_node(cx.graph, &factor);
        return Some((reduced, pivots, factor));
    }
    let mut work: Vec<Vec<NodeId>> = m.iter().map(|row| row.iter().map(|&e| cx.simplify(e)).collect()).collect();
    let mut field = Terms { cx };
    let (pivots, factor) = rref(&mut field, &mut work)?;
    Some((work, pivots, factor))
}

fn determinant(
    cx: &mut Cx<'_>,
    m: &[Vec<NodeId>],
) -> Option<NodeId> {
    let n = m.len();
    if m.iter().any(|r| r.len() != n) {
        return None;
    }
    if n == 0 {
        return Some(cx.graph.int(1));
    }
    if rational_matrix(cx.graph, m).is_none() && n <= 10 {
        // Symbolic: Laplace expansion over memoised minors is division
        // free, so the result is a polynomial in the entries.
        let value = laplace(cx.graph, m);
        return Some(normal(cx, value));
    }
    let (_, pivots, factor) = reduce(cx, m)?;
    Some(if pivots.len() == n { factor } else { cx.graph.int(0) })
}

/// The determinant by expansion along successive rows, memoising the minor
/// of every set of remaining columns: `O(n 2^n)` products.
fn laplace(
    graph: &mut Graph,
    m: &[Vec<NodeId>],
) -> NodeId {
    let n = m.len();
    // minors[mask] = det of rows n-|mask|.. restricted to the columns in mask
    let mut minors: std::collections::HashMap<u32, NodeId> = std::collections::HashMap::new();
    minors.insert(0, graph.int(1));
    let full: u32 = (1_u32 << n) - 1;
    let mut masks: Vec<u32> = (1..=full).collect();
    masks.sort_by_key(|mask| mask.count_ones());
    for mask in masks {
        let row = n - mask.count_ones() as usize;
        let mut terms = Vec::new();
        let mut sign_even = true;
        for col in 0..n {
            if mask & (1 << col) == 0 {
                continue;
            }
            let rest = minors[&(mask & !(1 << col))];
            let term = mul(graph, &[m[row][col], rest]);
            terms.push(if sign_even { term } else { neg(graph, term) });
            sign_even = !sign_even;
        }
        let value = add(graph, &terms);
        minors.insert(mask, value);
    }
    minors[&full]
}

/// The inverse of a square matrix of terms (exact or symbolic), for other
/// rule sets.
pub(crate) fn invert(
    cx: &mut Cx<'_>,
    m: &[Vec<NodeId>],
) -> Option<Vec<Vec<NodeId>>> {
    inverse(cx, m)
}

/// [`normal`] for other rule sets: simplified, over one denominator.
pub(crate) fn normal_form(
    cx: &mut Cx<'_>,
    node: NodeId,
) -> NodeId {
    normal(cx, node)
}

fn identity_rows(
    graph: &mut Graph,
    n: usize,
) -> Vec<Vec<NodeId>> {
    let (zero, one) = (graph.int(0), graph.int(1));
    (0..n).map(|i| (0..n).map(|j| if i == j { one } else { zero }).collect()).collect()
}

fn inverse(
    cx: &mut Cx<'_>,
    m: &[Vec<NodeId>],
) -> Option<Vec<Vec<NodeId>>> {
    let n = m.len();
    if m.iter().any(|r| r.len() != n) {
        return None;
    }
    if rational_matrix(cx.graph, m).is_none() && n <= 6 {
        // Symbolic: adj(A) / det(A) keeps entries as quotients of
        // polynomials in the entries.
        let det = determinant(cx, m)?;
        if (Terms { cx }).is_zero(&det) {
            return None;
        }
        let inverse_det = reciprocal(cx.graph, det);
        let mut out = vec![vec![cx.graph.int(0); n]; n];
        for (i, row) in m.iter().enumerate() {
            for j in 0..row.len() {
                let minor: Vec<Vec<NodeId>> = m
                    .iter()
                    .enumerate()
                    .filter(|&(r, _)| r != i)
                    .map(|(_, r)| r.iter().enumerate().filter(|&(c, _)| c != j).map(|(_, &e)| e).collect())
                    .collect();
                let cofactor = laplace(cx.graph, &minor);
                let cofactor = if (i + j) % 2 == 0 { cofactor } else { neg(cx.graph, cofactor) };
                let entry = mul(cx.graph, &[cofactor, inverse_det]);
                out[j][i] = normal(cx, entry);
            }
        }
        return Some(out);
    }
    let id = identity_rows(cx.graph, n);
    let augmented: Vec<Vec<NodeId>> = m.iter().zip(&id).map(|(r, i)| r.iter().chain(i).copied().collect()).collect();
    let (reduced, pivots, _) = reduce(cx, &augmented)?;
    if pivots.len() < n || pivots.iter().take(n).enumerate().any(|(k, &p)| p != k) {
        return None;
    }
    Some(reduced.into_iter().map(|row| row[n..].to_vec()).collect())
}

/// A basis of the null space.
fn nullspace(
    cx: &mut Cx<'_>,
    m: &[Vec<NodeId>],
) -> Option<Vec<Vec<NodeId>>> {
    let cols = m.first().map_or(0, Vec::len);
    let (reduced, pivots, _) = reduce(cx, m)?;
    let mut basis = Vec::new();
    for free in (0..cols).filter(|c| !pivots.contains(c)) {
        let mut v = vec![cx.graph.int(0); cols];
        v[free] = cx.graph.int(1);
        for (row, &p) in pivots.iter().enumerate() {
            let entry = reduced[row][free];
            v[p] = neg(cx.graph, entry);
            v[p] = cx.simplify(v[p]);
        }
        basis.push(v);
    }
    Some(basis)
}

fn transpose(m: &[Vec<NodeId>]) -> Vec<Vec<NodeId>> {
    let cols = m.first().map_or(0, Vec::len);
    (0..cols).map(|c| m.iter().map(|row| row[c]).collect()).collect()
}

fn matmul(
    graph: &mut Graph,
    a: &[Vec<NodeId>],
    b: &[Vec<NodeId>],
) -> Option<Vec<Vec<NodeId>>> {
    let inner = a.first().map_or(0, Vec::len);
    if inner != b.len() {
        return None;
    }
    let cols = b.first().map_or(0, Vec::len);
    let mut out = Vec::with_capacity(a.len());
    for row in a {
        let mut new_row = Vec::with_capacity(cols);
        for c in 0..cols {
            let terms: Vec<NodeId> = (0..inner).map(|k| mul(graph, &[row[k], b[k][c]])).collect();
            new_row.push(add(graph, &terms));
        }
        out.push(new_row);
    }
    Some(out)
}

/// Characteristic polynomial `det(l I - A)` as an expanded polynomial in
/// `l`, by the Faddeev–LeVerrier recurrence (no zero tests needed).
fn charpoly(
    cx: &mut Cx<'_>,
    m: &[Vec<NodeId>],
    l: NodeId,
) -> Option<NodeId> {
    let n = m.len();
    if m.iter().any(|r| r.len() != n) {
        return None;
    }
    // c_n = 1, M_0 = 0; M_k = A M_{k-1} + c_{n-k+1} I, c_{n-k} = -tr(A M_k)/k
    let mut coefficients = vec![cx.graph.int(0); n + 1];
    coefficients[n] = cx.graph.int(1);
    let zero_rows: Vec<Vec<NodeId>> = (0..n).map(|_| vec![cx.graph.int(0); n]).collect();
    let mut previous = zero_rows;
    for k in 1..=n {
        let product = matmul(cx.graph, m, &previous)?;
        let c = coefficients[n - k + 1];
        let mut current = product;
        for (i, row) in current.iter_mut().enumerate() {
            row[i] = add(cx.graph, &[row[i], c]);
        }
        for row in &mut current {
            for entry in row.iter_mut() {
                *entry = cx.simplify(*entry);
            }
        }
        let am = matmul(cx.graph, m, &current)?;
        let diagonal: Vec<NodeId> = (0..n).map(|i| am[i][i]).collect();
        let trace = add(cx.graph, &diagonal);
        let scale = cx.graph.num(Number::fraction(-1, i64::try_from(k).ok()?)?);
        let value = mul(cx.graph, &[scale, trace]);
        coefficients[n - k] = cx.simplify(value);
        previous = current;
    }
    let mut terms = Vec::with_capacity(n + 1);
    for (power, &c) in coefficients.iter().enumerate() {
        let term = match power {
            | 0 => c,
            | 1 => mul(cx.graph, &[c, l]),
            | p => {
                let e = cx.graph.int(i64::try_from(p).ok()?);
                let raised = cx.graph.node(core::POW, &[l, e]);
                mul(cx.graph, &[c, raised])
            },
        };
        terms.push(term);
    }
    let sum = add(cx.graph, &terms);
    let mut gens = Gens::default();
    gens.index(cx.graph, l);
    let poly = from_term(cx.graph, &mut gens, sum, Limits::default())?;
    let expanded = to_term(cx.graph, &gens, &poly);
    Some(cx.simplify(expanded))
}

/// Eigenvalues with multiplicities, from the roots of the characteristic
/// polynomial.
fn eigenvalues(
    cx: &mut Cx<'_>,
    m: &[Vec<NodeId>],
) -> Option<Vec<(NodeId, usize)>> {
    let fresh = cx.graph.interner_mut().fresh_symbol("lambda");
    let l = cx.graph.symbol_node(fresh);
    let p = charpoly(cx, m, l)?;
    let roots = solve_for(cx.graph, p, l, 0)?;
    // Multiplicity: the order of vanishing of p at the root.
    let mut gens = Gens::default();
    let gl = gens.index(cx.graph, l);
    let poly = from_term(cx.graph, &mut gens, p, Limits::default())?;
    let mut out = Vec::with_capacity(roots.len());
    for root in roots {
        let mut multiplicity = 0;
        let mut current: Poly = poly.clone();
        loop {
            let value = to_term(cx.graph, &gens, &current);
            let at = cx.graph.substitute(value, l, root);
            if !cx.is_zero(at) || multiplicity >= m.len() {
                break;
            }
            multiplicity += 1;
            current = current.derivative(gl);
        }
        out.push((root, multiplicity.max(1)));
    }
    let total: usize = out.iter().map(|(_, k)| k).sum();
    // Some roots have no closed form (or are complex): incomplete.
    (total == m.len()).then_some(out)
}

fn eigenvectors(
    cx: &mut Cx<'_>,
    m: &[Vec<NodeId>],
) -> Option<NodeId> {
    let values = eigenvalues(cx, m)?;
    let n = m.len();
    let mut entries = Vec::with_capacity(values.len());
    for (value, multiplicity) in values {
        let shifted: Vec<Vec<NodeId>> = (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| {
                        if i == j {
                            let negated = neg(cx.graph, value);
                            add(cx.graph, &[m[i][j], negated])
                        } else {
                            m[i][j]
                        }
                    })
                    .collect()
            })
            .collect();
        let basis = nullspace(cx, &shifted)?;
        let vectors: Vec<NodeId> = basis.iter().map(|v| list(cx.graph, v)).collect();
        let basis_node = list(cx.graph, &vectors);
        let k = cx.graph.int(i64::try_from(multiplicity).ok()?);
        entries.push(list(cx.graph, &[value, k, basis_node]));
    }
    Some(list(cx.graph, &entries))
}

/// `P A = L U` with partial pivoting (pivot on the first non-zero entry).
fn lu(
    cx: &mut Cx<'_>,
    m: &[Vec<NodeId>],
) -> Option<NodeId> {
    let n = m.len();
    if m.iter().any(|r| r.len() != n) {
        return None;
    }
    let mut u: Vec<Vec<NodeId>> = m.iter().map(|r| r.iter().map(|&e| cx.simplify(e)).collect()).collect();
    let mut l = identity_rows(cx.graph, n);
    let mut perm: Vec<usize> = (0..n).collect();
    let mut field = Terms { cx };
    for col in 0..n {
        let pivot_row = (col..n).find(|&r| !field.is_zero(&u[r][col]))?;
        if pivot_row != col {
            u.swap(pivot_row, col);
            perm.swap(pivot_row, col);
            for c in 0..col {
                let (a, b) = (l[pivot_row][c], l[col][c]);
                l[pivot_row][c] = b;
                l[col][c] = a;
            }
        }
        for r in col + 1..n {
            let factor = field.div(&u[r][col], &u[col][col])?;
            l[r][col] = factor;
            for c in col..n {
                let delta = field.mul(&factor, &u[col][c]);
                u[r][c] = field.sub(&u[r][c], &delta);
            }
        }
    }
    let graph = &mut *field.cx.graph;
    let (zero, one) = (graph.int(0), graph.int(1));
    let p: Vec<Vec<NodeId>> =
        perm.iter().map(|&k| (0..n).map(|j| if j == k { one } else { zero }).collect()).collect();
    let (p, l, u) = (matrix_term(graph, &p), matrix_term(graph, &l), matrix_term(graph, &u));
    Some(list(graph, &[p, l, u]))
}

/// `A = Q R` by Gram–Schmidt on the columns (full column rank).
fn qr(
    cx: &mut Cx<'_>,
    m: &[Vec<NodeId>],
) -> Option<NodeId> {
    let columns = transpose(m);
    let n = columns.len();
    let mut q: Vec<Vec<NodeId>> = Vec::with_capacity(n);
    let mut r = vec![vec![cx.graph.int(0); n]; n];
    for (j, column) in columns.iter().enumerate() {
        let mut v = column.clone();
        for (i, qi) in q.iter().enumerate() {
            let projection = dot(cx.graph, qi, column);
            let projection = cx.simplify(projection);
            r[i][j] = projection;
            for (vk, &qk) in v.iter_mut().zip(qi) {
                let delta = mul(cx.graph, &[projection, qk]);
                let delta = neg(cx.graph, delta);
                *vk = add(cx.graph, &[*vk, delta]);
                *vk = cx.simplify(*vk);
            }
        }
        let squared = dot(cx.graph, &v, &v);
        let squared = cx.simplify(squared);
        if cx.is_zero(squared) {
            return None;
        }
        let length = sqrt(cx.graph, squared)?;
        let length = cx.simplify(length);
        r[j][j] = length;
        let inverse = reciprocal(cx.graph, length);
        let normalised: Vec<NodeId> = v
            .iter()
            .map(|&vk| {
                let scaled = mul(cx.graph, &[vk, inverse]);
                cx.simplify(scaled)
            })
            .collect();
        q.push(normalised);
    }
    let q_rows = transpose(&q);
    let (q_node, r_node) = (matrix_term(cx.graph, &q_rows), matrix_term(cx.graph, &r));
    Some(list(cx.graph, &[q_node, r_node]))
}

/// Numeric singular value decomposition.
fn svd(
    graph: &mut Graph,
    m: &[Vec<NodeId>],
) -> Option<NodeId> {
    let rows = m.len();
    let cols = m.first().map_or(0, Vec::len);
    let mut data = Vec::with_capacity(rows * cols);
    for row in m {
        for &e in row {
            data.push(graph.number_of(e)?.to_f64());
        }
    }
    let matrix = Matrix::new(rows, cols, data);
    let FaerDecompositionResult::Svd { u, s, v } = matrix.decompose(FaerDecompositionType::Svd)? else {
        return None;
    };
    let to_term = |graph: &mut Graph, mat: &Matrix<f64>| {
        let rows: Vec<Vec<NodeId>> =
            (0..mat.rows()).map(|i| (0..mat.cols()).map(|j| graph.float(*mat.get(i, j))).collect()).collect();
        matrix_term(graph, &rows)
    };
    let (u, v) = (to_term(graph, &u), to_term(graph, &v));
    let values: Vec<NodeId> = s.iter().map(|&x| graph.float(x)).collect();
    let s = list(graph, &values);
    Some(list(graph, &[u, s, v]))
}

fn dot(
    graph: &mut Graph,
    a: &[NodeId],
    b: &[NodeId],
) -> NodeId {
    let terms: Vec<NodeId> = a.iter().zip(b).map(|(&x, &y)| mul(graph, &[x, y])).collect();
    add(graph, &terms)
}

/// Builds a request node `name(args...)` for another rule set.
fn request(
    graph: &mut Graph,
    name: &str,
    args: &[NodeId],
) -> Option<NodeId> {
    let op = graph.ops().lookup(name)?;
    graph.try_node(op, args)
}

fn derivative_simplified(
    cx: &mut Cx<'_>,
    f: NodeId,
    x: NodeId,
) -> Option<NodeId> {
    let term = best(cx.graph, f)?;
    let d = derivative(cx.graph, term, x)?;
    Some(cx.simplify(d))
}

/// The symbolic length of the derivative of a parametrised curve.
fn speed(
    cx: &mut Cx<'_>,
    curve: &[NodeId],
    t: NodeId,
) -> Option<(Vec<NodeId>, NodeId)> {
    let mut tangent = Vec::with_capacity(curve.len());
    for &c in curve {
        tangent.push(derivative_simplified(cx, c, t)?);
    }
    let squared = dot(cx.graph, &tangent, &tangent);
    let squared = cx.simplify(squared);
    let length = sqrt(cx.graph, squared)?;
    Some((tangent, cx.simplify(length)))
}

struct LinalgKernel {
    op: OpId,
    request: Request,
}

impl Kernel for LinalgKernel {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let args = cx.graph.children(node).to_vec();
        self.compute(cx, &args).map_or(Outcome::Pass, Outcome::Equal)
    }

    fn revisit(&self) -> bool {
        true
    }
}

impl LinalgKernel {
    #[allow(clippy::too_many_lines)]
    fn compute(
        &self,
        cx: &mut Cx<'_>,
        args: &[NodeId],
    ) -> Option<NodeId> {
        let arg = |i: usize| args.get(i).copied();
        match self.request {
            | Request::Matmul => {
                // A vector is a row on the left of a product and a column
                // anywhere else; the result of a product with a vector is
                // a vector again.
                let last = args.len().saturating_sub(1);
                let mut flat = false;
                let mut acc = Vec::new();
                for (k, &next) in args.iter().enumerate() {
                    let b = match operand(cx.graph, next)? {
                        | Operand::Matrix(m) => m,
                        | Operand::Vector(v) if k == 0 && last > 0 => {
                            flat = true;
                            vec![v]
                        },
                        | Operand::Vector(v) => {
                            flat = true;
                            v.into_iter().map(|e| vec![e]).collect()
                        },
                    };
                    acc = if k == 0 { b } else { matmul(cx.graph, &acc, &b)? };
                }
                normalise_entries(cx, &mut acc);
                if flat && (acc.len() == 1 || acc.iter().all(|r| r.len() == 1)) {
                    let items: Vec<NodeId> = acc.into_iter().flatten().collect();
                    return Some(list(cx.graph, &items));
                }
                Some(matrix_term(cx.graph, &acc))
            },
            | Request::Madd if operand(cx.graph, arg(0)?).is_some_and(|o| matches!(o, Operand::Vector(_))) => {
                let mut acc = vector(cx.graph, arg(0)?)?;
                for &next in args.get(1..)? {
                    let Operand::Vector(b) = operand(cx.graph, next)? else {
                        return None;
                    };
                    if b.len() != acc.len() {
                        return None;
                    }
                    for (x, y) in acc.iter_mut().zip(b) {
                        *x = add(cx.graph, &[*x, y]);
                    }
                }
                let mut rows = vec![acc];
                normalise_entries(cx, &mut rows);
                Some(list(cx.graph, &rows[0]))
            },
            | Request::Smul if operand(cx.graph, arg(1)?).is_some_and(|o| matches!(o, Operand::Vector(_))) => {
                let c = arg(0)?;
                let v = vector(cx.graph, arg(1)?)?;
                let scaled: Vec<NodeId> = v.iter().map(|&e| mul(cx.graph, &[c, e])).collect();
                Some(list(cx.graph, &scaled))
            },
            | Request::Madd => {
                let mut acc = matrix(cx.graph, arg(0)?)?;
                for &next in args.get(1..)? {
                    let b = matrix(cx.graph, next)?;
                    if b.len() != acc.len() || b.first().map(Vec::len) != acc.first().map(Vec::len) {
                        return None;
                    }
                    for (row, other) in acc.iter_mut().zip(&b) {
                        for (x, &y) in row.iter_mut().zip(other) {
                            *x = add(cx.graph, &[*x, y]);
                        }
                    }
                }
                normalise_entries(cx, &mut acc);
                Some(matrix_term(cx.graph, &acc))
            },
            | Request::Smul => {
                let c = arg(0)?;
                let m = matrix(cx.graph, arg(1)?)?;
                let scaled: Vec<Vec<NodeId>> =
                    m.iter().map(|row| row.iter().map(|&e| mul(cx.graph, &[c, e])).collect()).collect();
                Some(matrix_term(cx.graph, &scaled))
            },
            | Request::Transpose => {
                let m = matrix(cx.graph, arg(0)?)?;
                Some(matrix_term(cx.graph, &transpose(&m)))
            },
            | Request::Identity => {
                let n = usize::try_from(cx.graph.number_of(arg(0)?)?.to_i64()?).ok().filter(|&n| n <= 256)?;
                let rows = identity_rows(cx.graph, n);
                Some(matrix_term(cx.graph, &rows))
            },
            | Request::Zeros => {
                let size = |g: &Graph, n: NodeId| usize::try_from(g.number_of(n)?.to_i64()?).ok().filter(|&n| n <= 256);
                let (r, c) = (size(cx.graph, arg(0)?)?, size(cx.graph, arg(1)?)?);
                let zero = cx.graph.int(0);
                Some(matrix_term(cx.graph, &vec![vec![zero; c]; r]))
            },
            | Request::Dims => {
                let m = matrix(cx.graph, arg(0)?)?;
                let rows = cx.graph.int(i64::try_from(m.len()).ok()?);
                let cols = cx.graph.int(i64::try_from(m.first().map_or(0, Vec::len)).ok()?);
                Some(list(cx.graph, &[rows, cols]))
            },
            | Request::Det => {
                let m = matrix(cx.graph, arg(0)?)?;
                determinant(cx, &m)
            },
            | Request::Trace => {
                let m = matrix(cx.graph, arg(0)?)?;
                if m.iter().any(|r| r.len() != m.len()) {
                    return None;
                }
                let diagonal: Vec<NodeId> = (0..m.len()).map(|i| m[i][i]).collect();
                Some(add(cx.graph, &diagonal))
            },
            | Request::Rank => {
                let m = matrix(cx.graph, arg(0)?)?;
                let (_, pivots, _) = reduce(cx, &m)?;
                Some(cx.graph.int(i64::try_from(pivots.len()).ok()?))
            },
            | Request::Inverse => {
                let m = matrix(cx.graph, arg(0)?)?;
                let inv = inverse(cx, &m)?;
                Some(matrix_term(cx.graph, &inv))
            },
            | Request::Rref => {
                let m = matrix(cx.graph, arg(0)?)?;
                let (reduced, _, _) = reduce(cx, &m)?;
                Some(matrix_term(cx.graph, &reduced))
            },
            | Request::Nullspace => {
                let m = matrix(cx.graph, arg(0)?)?;
                let basis = nullspace(cx, &m)?;
                let vectors: Vec<NodeId> = basis.iter().map(|v| list(cx.graph, v)).collect();
                Some(list(cx.graph, &vectors))
            },
            | Request::Linsolve => {
                let a = matrix(cx.graph, arg(0)?)?;
                let b = vector(cx.graph, arg(1)?)?;
                if b.len() != a.len() {
                    return None;
                }
                let cols = a.first().map_or(0, Vec::len);
                let augmented: Vec<Vec<NodeId>> =
                    a.iter().zip(&b).map(|(row, &bi)| row.iter().copied().chain([bi]).collect()).collect();
                let (reduced, pivots, _) = reduce(cx, &augmented)?;
                if pivots.contains(&cols) {
                    // A pivot in the right-hand column: inconsistent.
                    return Some(list(cx.graph, &[]));
                }
                // Free columns become parameters `t1, t2, ...`.
                let mut solution = vec![cx.graph.int(0); cols];
                let mut parameter = 0;
                for col in (0..cols).filter(|c| !pivots.contains(c)) {
                    parameter += 1;
                    solution[col] = cx.graph.sym(&format!("t{parameter}"));
                }
                for (row, &p) in pivots.iter().enumerate() {
                    let mut terms = vec![reduced[row][cols]];
                    for col in (0..cols).filter(|c| !pivots.contains(c)) {
                        let product = mul(cx.graph, &[reduced[row][col], solution[col]]);
                        terms.push(neg(cx.graph, product));
                    }
                    let value = add(cx.graph, &terms);
                    solution[p] = cx.simplify(value);
                }
                Some(list(cx.graph, &solution))
            },
            | Request::Charpoly => {
                let m = matrix(cx.graph, arg(0)?)?;
                charpoly(cx, &m, arg(1)?)
            },
            | Request::Eigenvals => {
                let m = matrix(cx.graph, arg(0)?)?;
                let values = eigenvalues(cx, &m)?;
                let mut flat = Vec::new();
                for (value, multiplicity) in values {
                    flat.extend(std::iter::repeat_n(value, multiplicity));
                }
                Some(list(cx.graph, &flat))
            },
            | Request::Eigenvects => {
                let m = matrix(cx.graph, arg(0)?)?;
                eigenvectors(cx, &m)
            },
            | Request::Lu => {
                let m = matrix(cx.graph, arg(0)?)?;
                lu(cx, &m)
            },
            | Request::Qr => {
                let m = matrix(cx.graph, arg(0)?)?;
                qr(cx, &m)
            },
            | Request::Svd => {
                let m = matrix(cx.graph, arg(0)?)?;
                svd(cx.graph, &m)
            },
            | Request::Dot => {
                let (a, b) = (vector(cx.graph, arg(0)?)?, vector(cx.graph, arg(1)?)?);
                (a.len() == b.len()).then(|| dot(cx.graph, &a, &b))
            },
            | Request::Cross => {
                let (a, b) = (vector(cx.graph, arg(0)?)?, vector(cx.graph, arg(1)?)?);
                let (&[a1, a2, a3], &[b1, b2, b3]) = (a.as_slice(), b.as_slice()) else {
                    return None;
                };
                let component = |g: &mut Graph, p: NodeId, q: NodeId, r: NodeId, s: NodeId| {
                    let first = mul(g, &[p, q]);
                    let second = mul(g, &[r, s]);
                    let second = neg(g, second);
                    add(g, &[first, second])
                };
                let c1 = component(cx.graph, a2, b3, a3, b2);
                let c2 = component(cx.graph, a3, b1, a1, b3);
                let c3 = component(cx.graph, a1, b2, a2, b1);
                Some(list(cx.graph, &[c1, c2, c3]))
            },
            | Request::Norm => {
                let v = vector(cx.graph, arg(0)?)?;
                let squared = dot(cx.graph, &v, &v);
                sqrt(cx.graph, squared)
            },
            | Request::Normalize => {
                let v = vector(cx.graph, arg(0)?)?;
                let squared = dot(cx.graph, &v, &v);
                let length = sqrt(cx.graph, squared)?;
                let inverse = reciprocal(cx.graph, length);
                let scaled: Vec<NodeId> = v.iter().map(|&e| mul(cx.graph, &[e, inverse])).collect();
                Some(list(cx.graph, &scaled))
            },
            | Request::Angle => {
                let (a, b) = (vector(cx.graph, arg(0)?)?, vector(cx.graph, arg(1)?)?);
                let ab = dot(cx.graph, &a, &b);
                let aa = dot(cx.graph, &a, &a);
                let bb = dot(cx.graph, &b, &b);
                let product = mul(cx.graph, &[aa, bb]);
                let length = sqrt(cx.graph, product)?;
                let inverse = reciprocal(cx.graph, length);
                let cosine = mul(cx.graph, &[ab, inverse]);
                request(cx.graph, "acos", &[cosine])
            },
            | Request::Project => {
                // projection of a onto b
                let (a, b) = (vector(cx.graph, arg(0)?)?, vector(cx.graph, arg(1)?)?);
                let ab = dot(cx.graph, &a, &b);
                let bb = dot(cx.graph, &b, &b);
                let inverse = reciprocal(cx.graph, bb);
                let scale = mul(cx.graph, &[ab, inverse]);
                let scaled: Vec<NodeId> = b.iter().map(|&e| mul(cx.graph, &[scale, e])).collect();
                Some(list(cx.graph, &scaled))
            },
            | Request::Outer => {
                let (a, b) = (vector(cx.graph, arg(0)?)?, vector(cx.graph, arg(1)?)?);
                let rows: Vec<Vec<NodeId>> =
                    a.iter().map(|&x| b.iter().map(|&y| mul(cx.graph, &[x, y])).collect()).collect();
                Some(matrix_term(cx.graph, &rows))
            },
            | Request::Grad => {
                let vars = vector(cx.graph, arg(1)?)?;
                let f = arg(0)?;
                let mut components = Vec::with_capacity(vars.len());
                for &x in &vars {
                    components.push(derivative_simplified(cx, f, x)?);
                }
                Some(list(cx.graph, &components))
            },
            | Request::Div => {
                let (field, vars) = (vector(cx.graph, arg(0)?)?, vector(cx.graph, arg(1)?)?);
                if field.len() != vars.len() {
                    return None;
                }
                let mut terms = Vec::with_capacity(vars.len());
                for (&f, &x) in field.iter().zip(&vars) {
                    terms.push(derivative_simplified(cx, f, x)?);
                }
                Some(add(cx.graph, &terms))
            },
            | Request::Curl => {
                let (field, vars) = (vector(cx.graph, arg(0)?)?, vector(cx.graph, arg(1)?)?);
                let (&[fx, fy, fz], &[x, y, z]) = (field.as_slice(), vars.as_slice()) else {
                    return None;
                };
                let d = |cx: &mut Cx<'_>, f: NodeId, v: NodeId| derivative_simplified(cx, f, v);
                let pairs = [(fz, y, fy, z), (fx, z, fz, x), (fy, x, fx, y)];
                let mut components = Vec::with_capacity(3);
                for (p, u, q, w) in pairs {
                    let first = d(cx, p, u)?;
                    let second = d(cx, q, w)?;
                    let second = neg(cx.graph, second);
                    components.push(add(cx.graph, &[first, second]));
                }
                Some(list(cx.graph, &components))
            },
            | Request::Laplacian => {
                let vars = vector(cx.graph, arg(1)?)?;
                let f = arg(0)?;
                let mut terms = Vec::with_capacity(vars.len());
                for &x in &vars {
                    let first = derivative_simplified(cx, f, x)?;
                    terms.push(derivative_simplified(cx, first, x)?);
                }
                Some(add(cx.graph, &terms))
            },
            | Request::Jacobian => {
                let (field, vars) = (vector(cx.graph, arg(0)?)?, vector(cx.graph, arg(1)?)?);
                let mut rows = Vec::with_capacity(field.len());
                for &f in &field {
                    let mut row = Vec::with_capacity(vars.len());
                    for &x in &vars {
                        row.push(derivative_simplified(cx, f, x)?);
                    }
                    rows.push(row);
                }
                Some(matrix_term(cx.graph, &rows))
            },
            | Request::Hessian => {
                let vars = vector(cx.graph, arg(1)?)?;
                let f = arg(0)?;
                let mut rows = Vec::with_capacity(vars.len());
                for &x in &vars {
                    let first = derivative_simplified(cx, f, x)?;
                    let mut row = Vec::with_capacity(vars.len());
                    for &y in &vars {
                        row.push(derivative_simplified(cx, first, y)?);
                    }
                    rows.push(row);
                }
                Some(matrix_term(cx.graph, &rows))
            },
            | Request::Directional => {
                let vars = vector(cx.graph, arg(1)?)?;
                let direction = vector(cx.graph, arg(2)?)?;
                if direction.len() != vars.len() {
                    return None;
                }
                let f = arg(0)?;
                let mut gradient = Vec::with_capacity(vars.len());
                for &x in &vars {
                    gradient.push(derivative_simplified(cx, f, x)?);
                }
                let squared = dot(cx.graph, &direction, &direction);
                let length = sqrt(cx.graph, squared)?;
                let inverse = reciprocal(cx.graph, length);
                let projection = dot(cx.graph, &gradient, &direction);
                Some(mul(cx.graph, &[projection, inverse]))
            },
            | Request::LineIntegral | Request::LineIntegralVec => {
                // f(r(t)) |r'(t)| dt, or F(r(t)) . r'(t) dt
                let (integrand, curve_node, t, a, b) = (arg(0)?, arg(1)?, arg(2)?, arg(3)?, arg(4)?);
                let curve = vector(cx.graph, curve_node)?;
                let (tangent, length) = speed(cx, &curve, t)?;
                let coordinates = coordinate_symbols(cx.graph, curve.len());
                let substitute = |graph: &mut Graph, mut term: NodeId| {
                    for (&symbol, &value) in coordinates.iter().zip(&curve) {
                        term = graph.substitute(term, symbol, value);
                    }
                    term
                };
                let body = if self.request == Request::LineIntegral {
                    let along = substitute(cx.graph, integrand);
                    mul(cx.graph, &[along, length])
                } else {
                    let field = vector(cx.graph, integrand)?;
                    let along: Vec<NodeId> = field.iter().map(|&f| substitute(cx.graph, f)).collect();
                    dot(cx.graph, &along, &tangent)
                };
                let body = cx.simplify(body);
                request(cx.graph, "defint", &[body, t, a, b])
            },
            | Request::SurfaceIntegral => {
                // f(r(u, v)) |r_u x r_v| du dv
                let surface = vector(cx.graph, arg(1)?)?;
                let (u, v) = (arg(2)?, arg(3)?);
                if surface.len() != 3 {
                    return None;
                }
                let mut ru = Vec::with_capacity(3);
                let mut rv = Vec::with_capacity(3);
                for &c in &surface {
                    ru.push(derivative_simplified(cx, c, u)?);
                    rv.push(derivative_simplified(cx, c, v)?);
                }
                let ru_node = list(cx.graph, &ru);
                let rv_node = list(cx.graph, &rv);
                let normal = Self { op: self.op, request: Request::Cross }.compute(cx, &[ru_node, rv_node])?;
                let normal = vector(cx.graph, normal)?;
                let squared = dot(cx.graph, &normal, &normal);
                let squared = cx.simplify(squared);
                let area = sqrt(cx.graph, squared)?;
                let coordinates = coordinate_symbols(cx.graph, 3);
                let mut along = arg(0)?;
                for (&symbol, &value) in coordinates.iter().zip(&surface) {
                    along = cx.graph.substitute(along, symbol, value);
                }
                let body = mul(cx.graph, &[along, area]);
                let body = cx.simplify(body);
                let inner = request(cx.graph, "defint", &[body, u, arg(4)?, arg(5)?])?;
                request(cx.graph, "defint", &[inner, v, arg(6)?, arg(7)?])
            },
            | Request::VolumeIntegral => {
                // volume_integral(f, list(x, y, z), list(list(x0, x1), ...))
                let vars = vector(cx.graph, arg(1)?)?;
                let bounds = matrix(cx.graph, arg(2)?)?;
                if bounds.len() != vars.len() || bounds.iter().any(|b| b.len() != 2) {
                    return None;
                }
                let mut body = arg(0)?;
                for (&x, limits) in vars.iter().zip(&bounds) {
                    body = request(cx.graph, "defint", &[body, x, limits[0], limits[1]])?;
                }
                Some(body)
            },
        }
    }
}

/// The coordinate symbols `x, y, z` (or `x1..xn` beyond three) that
/// integrands of curve and surface integrals are written in.
fn coordinate_symbols(
    graph: &mut Graph,
    n: usize,
) -> Vec<NodeId> {
    if n <= 3 {
        ["x", "y", "z"].iter().take(n).map(|name| graph.sym(name)).collect()
    } else {
        (1..=n).map(|k| graph.sym(&format!("x{k}"))).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rules::testing::eval;
    use crate::rules::testing::numeric;
    use crate::rules::testing::reduce_with;
    use crate::rules::testing::simplify;

    fn run(src: &str) -> String {
        simplify(&[linalg()], src)
    }

    #[test]
    fn arithmetic() {
        assert_eq!(run("matmul(list(list(1, 2), list(3, 4)), list(list(0, 1), list(1, 0)))"), "list(list(2, 1), list(4, 3))");
        assert_eq!(run("madd(list(list(1, 2)), list(list(a, b)))"), "list(list(a + 1, b + 2))");
        assert_eq!(run("smul(2, list(list(1, x)))"), "list(list(2, 2*x))");
        assert_eq!(run("transpose(list(list(1, 2, 3), list(4, 5, 6)))"), "list(list(1, 4), list(2, 5), list(3, 6))");
        assert_eq!(run("identity(2)"), "list(list(1, 0), list(0, 1))");
        assert_eq!(run("trace(list(list(a, b), list(c, d)))"), "a + d");
        let (text, reduced) = reduce_with(&[linalg()], "matmul(list(list(1, 2)), list(list(1, 2)))", &[]);
        assert!(!reduced, "shapes do not match: {text}");
    }

    #[test]
    fn determinants_inverses_and_rank() {
        assert_eq!(run("det(list(list(1, 2), list(3, 4)))"), "-2");
        assert_eq!(run("det(list(list(a, b), list(c, d)))"), "a*d - b*c");
        assert_eq!(run("det(list(list(2, 0, 1), list(1, 3, 2), list(1, 1, 2)))"), "6");
        assert_eq!(run("det(list(list(1, 2), list(2, 4)))"), "0");
        assert_eq!(run("inverse(list(list(1, 2), list(3, 4)))"), "list(list(-2, 1), list(3/2, -1/2))");
        assert_eq!(run("inverse(list(list(a, 0), list(0, b)))"), "list(list(1/a, 0), list(0, 1/b))");
        let (text, reduced) = reduce_with(&[linalg()], "inverse(list(list(1, 2), list(2, 4)))", &[]);
        assert!(!reduced, "a singular matrix has no inverse: {text}");
        assert_eq!(run("rank(list(list(1, 2, 3), list(2, 4, 6), list(1, 0, 1)))"), "2");
        assert_eq!(run("rref(list(list(1, 2, 3), list(4, 5, 6)))"), "list(list(1, 0, -1), list(0, 1, 2))");
        assert_eq!(run("nullspace(list(list(1, 2, 3), list(4, 5, 6)))"), "list(list(1, -2, 1))");
    }

    #[test]
    fn inverse_times_matrix_is_identity() {
        let a = "list(list(2, 1, 1), list(1, 3, 2), list(1, 0, 0))";
        assert_eq!(run(&format!("matmul(inverse({a}), {a})")), "list(list(1, 0, 0), list(0, 1, 0), list(0, 0, 1))");
        let s = "list(list(a, 1), list(1, a))";
        assert_eq!(run(&format!("matmul({s}, inverse({s}))")), "list(list(1, 0), list(0, 1))");
    }

    #[test]
    fn linear_systems() {
        assert_eq!(run("linsolve(list(list(2, 1), list(1, 3)), list(3, 5))"), "list(4/5, 7/5)");
        assert_eq!(run("linsolve(list(list(1, 1), list(1, 1)), list(1, 2))"), "list()");
        assert_eq!(run("linsolve(list(list(a, 0), list(0, b)), list(1, 1))"), "list(1/a, 1/b)");
        assert_eq!(run("linsolve(list(list(1, 1, 1)), list(3))"), "list(3 - t1 - t2, t1, t2)");
        assert_eq!(run("dims(list(list(1, 2, 3), list(4, 5, 6)))"), "list(2, 3)");
    }

    #[test]
    fn eigen() {
        assert_eq!(run("charpoly(list(list(1, 2), list(3, 4)), l)"), "l^2 - 5*l - 2");
        assert_eq!(run("eigenvals(list(list(2, 0), list(0, 3)))"), "list(2, 3)");
        assert_eq!(run("eigenvals(list(list(2, 1), list(1, 2)))"), "list(1, 3)");
        assert_eq!(run("eigenvals(list(list(1, 1), list(0, 1)))"), "list(1, 1)");
        assert_eq!(run("eigenvects(list(list(2, 1), list(1, 2)))"), "list(list(1, 1, list(list(-1, 1))), list(3, 1, list(list(1, 1))))");
        // Symbolic entries.
        assert_eq!(run("charpoly(list(list(a, b), list(b, a)), l)"), "a^2 - 2*a*l - b^2 + l^2");
    }

    #[test]
    fn decompositions() {
        // P A = L U must hold.
        let a = "list(list(0, 2, 1), list(1, 1, 0), list(2, 1, 3))";
        let lu = run(&format!("lu({a})"));
        let check = run(&format!("madd(matmul(item0, {a}), smul(-1, matmul(item1, item2)))").replace("item0", &part(&lu, 0)).replace("item1", &part(&lu, 1)).replace("item2", &part(&lu, 2)));
        assert_eq!(check, "list(list(0, 0, 0), list(0, 0, 0), list(0, 0, 0))", "{lu}");
        // A = Q R with orthonormal Q.
        let qr = run("qr(list(list(3, 1), list(4, 2)))");
        let product = run(&format!("matmul({}, {})", part(&qr, 0), part(&qr, 1)));
        assert_eq!(product, "list(list(3, 1), list(4, 2))", "{qr}");
        // Numeric SVD: singular values of diag(3, -2).
        let svd = run("svd(list(list(3.0, 0.0), list(0.0, -2.0)))");
        let values = part(&svd, 1);
        assert!(values.contains('3') && values.contains('2'), "{svd}");
    }

    /// The `index`-th top-level item of a `list(...)` text.
    fn part(
        text: &str,
        index: usize,
    ) -> String {
        let inner = text.strip_prefix("list(").and_then(|s| s.strip_suffix(')')).unwrap_or("");
        let mut depth = 0_i32;
        let mut start = 0;
        let mut items = Vec::new();
        for (i, ch) in inner.char_indices() {
            match ch {
                | '(' => depth += 1,
                | ')' => depth -= 1,
                | ',' if depth == 0 => {
                    items.push(inner[start..i].trim().to_owned());
                    start = i + 1;
                },
                | _ => {},
            }
        }
        items.push(inner[start..].trim().to_owned());
        items.get(index).cloned().unwrap_or_default()
    }

    #[test]
    fn vectors() {
        assert_eq!(run("dot(list(1, 2, 3), list(4, 5, 6))"), "32");
        assert_eq!(run("cross(list(1, 0, 0), list(0, 1, 0))"), "list(0, 0, 1)");
        assert_eq!(run("norm(list(3, 4))"), "5");
        assert_eq!(run("normalize(list(3, 4))"), "list(3/5, 4/5)");
        assert_eq!(run("angle(list(1, 0), list(0, 1))"), "1/2*pi");
        assert_eq!(run("project(list(1, 1), list(2, 0))"), "list(1, 0)");
        assert_eq!(run("outer(list(1, 2), list(3, 4))"), "list(list(3, 4), list(6, 8))");
    }

    #[test]
    fn vector_calculus() {
        assert_eq!(run("grad(x^2*y + z, list(x, y, z))"), "list(2*x*y, x^2, 1)");
        assert_eq!(run("div(list(x, y^2, z*x), list(x, y, z))"), "x + 2*y + 1");
        assert_eq!(run("curl(list(-y, x, 0), list(x, y, z))"), "list(0, 0, 2)");
        assert_eq!(run("laplacian(x^2 + y^2 + z^2, list(x, y, z))"), "6");
        assert_eq!(run("jacobian(list(x*y, x + y), list(x, y))"), "list(list(y, x), list(1, 1))");
        assert_eq!(run("hessian(x^2*y, list(x, y))"), "list(list(2*y, 2*x), list(2*x, 0))");
        assert_eq!(run("directional(x^2 + y^2, list(x, y), list(3, 4))"), "6/5*x + 8/5*y");
        // div curl = 0 and curl grad = 0 for any field.
        assert_eq!(run("div(curl(list(x*y*z, sin(x)*y, exp(z)*x), list(x, y, z)), list(x, y, z))"), "0");
        assert_eq!(run("curl(grad(x^2*sin(y)*z, list(x, y, z)), list(x, y, z))"), "list(0, 0, 0)");
    }

    #[test]
    fn integrals_over_curves_surfaces_and_volumes() {
        // Arc length of the unit circle.
        let (value, _) = numeric(&[linalg()], "line_integral(1, list(cos(t), sin(t)), t, 0, 2*pi)", &[], 1e-10);
        assert!((value - 2.0 * std::f64::consts::PI).abs() < 1e-9, "{value}");
        // Work of F = (-y, x) around the unit circle: 2 pi.
        assert_eq!(run("line_integral_vec(list(-y, x), list(cos(t), sin(t)), t, 0, 2*pi)"), "2*pi");
        // Area of the unit sphere.
        let (value, _) = numeric(
            &[linalg()],
            "surface_integral(1, list(sin(u)*cos(v), sin(u)*sin(v), cos(u)), u, v, 0, pi, 0, 2*pi)",
            &[],
            1e-8,
        );
        assert!((value - 4.0 * std::f64::consts::PI).abs() < 1e-6, "{value}");
        assert_eq!(run("volume_integral(x*y*z, list(x, y, z), list(list(0, 1), list(0, 2), list(0, 3)))"), "9/2");
        let _ = eval;
    }
}
