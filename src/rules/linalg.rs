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
//! | `lu(A)`, `qr(A)`, `svd(A)` | `list(P, L, U)` with `P A = L U`; `list(Q, R)`; `list(U, S, V)` (numeric for float entries, exact otherwise) |
//! | `dot`, `cross`, `norm`, `normalize`, `angle`, `project`, `outer` | vector operations |
//! | `grad(f, vars)`, `div(F, vars)`, `curl(F, vars)`, `laplacian(f, vars)`, `jacobian(F, vars)`, `hessian(f, vars)`, `directional(f, vars, v)` | vector calculus |
//! | `kron`, `diag`, `vandermonde`, `hilbert`, `companion`, `adjugate`, `cofactors`, `colspace`, `rowspace`, `left_nullspace`, `pinv`, `lstsq`, `minpoly`, `matpoly`, `cayley_hamilton`, `jordan`, `jordan_blocks`, `matfun`, `matexp`, `mpow`, `matsqrt`, `spectral`, `ldl`, `cholesky`, `orthogonalize`, `orthonormalize`, `smith`, `invariant_factors`, `hermite_form`, `is_symmetric`, `is_orthogonal`, `is_positive_definite` | exact matrix functions and factorisations, documented in the table below |
//! | `line_integral(f, curve, t, a, b)`, `line_integral_vec(F, curve, t, a, b)`, `surface_integral(f, surface, u, v, ua, ub, va, vb)`, `volume_integral(f, vars, bounds)` | integrals over parametrised curves, surfaces and boxes, reduced to definite integrals |


use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Zero;

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
use crate::graph::dim::DimCtx;
use crate::graph::dim::DimError;
use crate::graph::dim::DimRule;
use crate::graph::dim::DimTerm;
use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::kernels::matrix::FaerDecompositionResult;
use crate::kernels::matrix::FaerDecompositionType;
use crate::kernels::matrix::Matrix;
use crate::rules::poly::best;
use crate::rules::poly::repr::Gens;
use crate::rules::poly::repr::Limits;
use crate::rules::poly::repr::Poly;
use crate::rules::poly::repr::from_term;
use crate::rules::poly::repr::to_term;

use super::calculus::calculus;
use super::calculus::derivative;
use super::solve::solve;
use super::solve::solve_for;

/// The linear algebra rule set.
#[must_use]
pub fn linalg() -> RuleSet {
    RuleSet::new("linalg", install)
        .needs(calculus())
        .needs(solve())
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

/// The dimension signature of a linear-algebra request.
///
/// Matrix and vector arguments have the dimension of their entries (the
/// `List` rule unifies the leaves), so `matmul`/`dot`/`kron` multiply while
/// `madd` adds. `det` is `entry^n`, `inverse` is `1/entry`, `linsolve` is
/// `rhs / matrix`, and the vector-calculus requests differentiate with
/// respect to their variable argument.
fn request_rule(request: Request) -> DimRule {
    match request {
        | Request::Matmul | Request::Smul | Request::Dot | Request::Cross | Request::Outer => {
            DimRule::Product
        },
        | Request::Madd => DimRule::Sum,
        | Request::Transpose
        | Request::Trace
        | Request::Rref
        | Request::Norm
        | Request::Project => DimRule::Same(0),
        | Request::Identity
        | Request::Zeros
        | Request::Nullspace
        | Request::Eigenvals
        | Request::Eigenvects
        | Request::Lu
        | Request::Qr
        | Request::Svd => DimRule::List,
        | Request::Dims | Request::Rank | Request::Normalize => {
            DimRule::Custom(dimensionless_result)
        },
        | Request::Det => DimRule::Custom(det_dimension),
        | Request::Inverse => DimRule::Inverse,
        | Request::Linsolve => DimRule::Custom(quotient_dimension),
        | Request::Charpoly => DimRule::NonMetric,
        | Request::Angle => DimRule::Angle,
        | Request::Grad
        | Request::Div
        | Request::Curl
        | Request::Jacobian
        | Request::Directional => DimRule::Diff { f: 0, x: 1 },
        | Request::Laplacian | Request::Hessian => DimRule::DiffN { f: 0, x: 1, n: 2 },
        | Request::LineIntegral | Request::LineIntegralVec => DimRule::Integral { f: 0, x: 2 },
        | Request::SurfaceIntegral => DimRule::Custom(surface_integral_dimension),
        | Request::VolumeIntegral => DimRule::Custom(volume_integral_dimension),
    }
}

/// `dims`, `rank` and `normalize` are pure numbers or directions: the
/// result is dimensionless whatever the entries are.
#[allow(clippy::missing_const_for_fn, clippy::unnecessary_wraps)]
fn dimensionless_result(
    _ctx: &DimCtx<'_>,
    _args: &[DimTerm],
) -> Result<DimTerm, DimError> {
    Ok(DimTerm::dimensionless())
}

/// `det(A)` for an `n x n` matrix is `entry^n`; an unknown shape (a symbolic
/// matrix argument) falls back to `entry`.
#[allow(clippy::unnecessary_wraps)]
fn det_dimension(
    ctx: &DimCtx<'_>,
    args: &[DimTerm],
) -> Result<DimTerm, DimError> {
    let Some(entry) = args.first() else {
        return Ok(ctx.any());
    };
    let order = matrix_order(ctx).unwrap_or(1);
    Ok(entry.pow(&BigRational::from_integer(BigInt::from(order))))
}

/// `linsolve(A, b)` is `b / A`, not `A / b`: `Quotient` would be backwards.
#[allow(clippy::unnecessary_wraps)]
fn quotient_dimension(
    _ctx: &DimCtx<'_>,
    args: &[DimTerm],
) -> Result<DimTerm, DimError> {
    Ok(match (args.first(), args.get(1)) {
        | (Some(matrix), Some(rhs)) => rhs.mul(&matrix.inverse()),
        | _ => DimTerm::dimensionless(),
    })
}

#[allow(clippy::unnecessary_wraps)]
fn surface_integral_dimension(
    _ctx: &DimCtx<'_>,
    args: &[DimTerm],
) -> Result<DimTerm, DimError> {
    let field = args.first().cloned().unwrap_or_else(DimTerm::dimensionless);
    let first = args.get(2).cloned().unwrap_or_else(DimTerm::dimensionless);
    let second = args.get(3).cloned().unwrap_or_else(DimTerm::dimensionless);
    Ok(field.mul(&first).mul(&second))
}

#[allow(clippy::unnecessary_wraps)]
fn volume_integral_dimension(
    _ctx: &DimCtx<'_>,
    args: &[DimTerm],
) -> Result<DimTerm, DimError> {
    let field = args.first().cloned().unwrap_or_else(DimTerm::dimensionless);
    let volume = args
        .get(1)
        .map_or_else(DimTerm::dimensionless, |variables| {
            variables.pow(&BigRational::from_integer(BigInt::from(3)))
        });
    Ok(field.mul(&volume))
}

/// The order of the literal matrix argument of a `det`, when it is one.
fn matrix_order(ctx: &DimCtx<'_>) -> Option<usize> {
    let matrix = ctx.graph.children(ctx.node).first().copied()?;
    if ctx.graph.op(matrix) != core::LIST {
        return None;
    }
    let rows = ctx.graph.children(matrix);
    let first_row = rows.first().copied()?;
    (ctx.graph.op(first_row) == core::LIST).then_some(rows.len())
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
        (
            "line_integral_vec",
            Arity::Fixed(5),
            Request::LineIntegralVec,
        ),
        (
            "surface_integral",
            Arity::Fixed(8),
            Request::SurfaceIntegral,
        ),
        ("volume_integral", Arity::Fixed(3), Request::VolumeIntegral),
    ];
    for (name, arity, request) in table {
        let op = i.op(OpDescriptor::new(name, arity)
            .flags(OpFlags::HEAVY)
            .cost(100))?;
        i.dim(op, request_rule(request));
        i.kernel(
            &format!("linalg/{name}"),
            Tier::Reduce,
            LinalgKernel { op, request },
        );
    }
    i.kernel("linalg/broadcast", Tier::Normalize, Broadcast);
    install_ext(i)?;
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
        graph
            .enodes(graph.find(node))
            .find(|&n| graph.op(n) == core::LIST)
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
        let lists: Vec<Option<NodeId>> =
            children.iter().map(|&c| Self::as_list(graph, c)).collect();
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
        let mut found = lists
            .iter()
            .enumerate()
            .filter_map(|(k, l)| l.map(|l| (k, l)));
        let (Some((at, l)), None) = (found.next(), found.next()) else {
            return Outcome::Pass;
        };
        let others: Vec<NodeId> = children
            .iter()
            .enumerate()
            .filter(|&(k, _)| k != at)
            .map(|(_, &c)| c)
            .collect();
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
    if let Some(m) = matrix(graph, node)
        && !m.is_empty()
    {
        return Some(Operand::Matrix(m));
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

pub(crate) fn matrix_term(
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
        if graph.op(n) == core::POW
            && let Some(e) = graph.children(n).get(1).and_then(|&e| graph.number_of(e))
            && e.to_f64() < 0.0
        {
            return true;
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
            if cx.graph.op(together) == cx.graph.op(cancel) {
                simplified
            } else {
                together
            }
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
        #[allow(clippy::needless_range_loop)]
        // index is used for more than one array / arithmetic; iterator form would not be clearer
        for c in 0..cols {
            m[row][c] = field.div(&m[row][c], &pivot)?;
        }
        for r in 0..rows {
            if r == row || field.is_zero(&m[r][col]) {
                continue;
            }
            let scale = m[r][col].clone();
            #[allow(clippy::needless_range_loop)]
            // index is used for more than one array / arithmetic; iterator form would not be clearer
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
    m.iter()
        .map(|row| {
            row.iter()
                .map(|&e| graph.number_of(e).and_then(Number::to_rational))
                .collect()
        })
        .collect()
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
        let reduced = exact
            .iter()
            .map(|row| row.iter().map(|v| rational_node(cx.graph, v)).collect())
            .collect();
        let factor = rational_node(cx.graph, &factor);
        return Some((reduced, pivots, factor));
    }
    let mut work: Vec<Vec<NodeId>> = m
        .iter()
        .map(|row| row.iter().map(|&e| cx.simplify(e)).collect())
        .collect();
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
    Some(if pivots.len() == n {
        factor
    } else {
        cx.graph.int(0)
    })
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
            terms.push(if sign_even {
                term
            } else {
                neg(graph, term)
            });
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
    (0..n)
        .map(|i| (0..n).map(|j| if i == j { one } else { zero }).collect())
        .collect()
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
            #[allow(clippy::needless_range_loop)]
            // index is used for more than one array / arithmetic; iterator form would not be clearer
            for j in 0..row.len() {
                let minor: Vec<Vec<NodeId>> = m
                    .iter()
                    .enumerate()
                    .filter(|&(r, _)| r != i)
                    .map(|(_, r)| {
                        r.iter()
                            .enumerate()
                            .filter(|&(c, _)| c != j)
                            .map(|(_, &e)| e)
                            .collect()
                    })
                    .collect();
                let cofactor = laplace(cx.graph, &minor);
                let cofactor = if (i + j) % 2 == 0 {
                    cofactor
                } else {
                    neg(cx.graph, cofactor)
                };
                let entry = mul(cx.graph, &[cofactor, inverse_det]);
                out[j][i] = normal(cx, entry);
            }
        }
        return Some(out);
    }
    let id = identity_rows(cx.graph, n);
    let augmented: Vec<Vec<NodeId>> = m
        .iter()
        .zip(&id)
        .map(|(r, i)| r.iter().chain(i).copied().collect())
        .collect();
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
    (0..cols)
        .map(|c| m.iter().map(|row| row[c]).collect())
        .collect()
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
        #[allow(clippy::needless_range_loop)]
        // index is used for more than one array / arithmetic; iterator form would not be clearer
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

/// Eigenvalues by the shifted QR algorithm, for a matrix of numbers with
/// at least one float entry (or in a numeric run) whose characteristic
/// polynomial has no closed-form roots. Complex pairs become `re + im*I`.
fn numeric_eigenvalues(
    cx: &mut Cx<'_>,
    m: &[Vec<NodeId>],
) -> Option<NodeId> {
    let mut floats = cx.env.numeric;
    let mut rows = Vec::with_capacity(m.len());
    for row in m {
        let mut values = Vec::with_capacity(row.len());
        for &e in row {
            let number = cx.graph.number_of(e)?;
            floats |= matches!(number, Number::Float(_));
            values.push(number.to_f64());
        }
        rows.push(values);
    }
    if !floats {
        return None;
    }
    let a = crate::kernels::dense::Mat::from_rows(&rows).ok()?;
    let mut values = crate::kernels::dense::eigenvalues(&a).ok()?;
    values.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.total_cmp(&y.1)));
    let unit = cx.graph.ops().lookup("I");
    let mut terms = Vec::with_capacity(values.len());
    for (re, im) in values {
        let real = cx.graph.float(re);
        if im.abs() < f64::MIN_POSITIVE {
            terms.push(real);
        } else {
            let unit = match unit {
                | Some(op) => cx.graph.node(op, &[]),
                | None => {
                    let minus_one = cx.graph.int(-1);
                    sqrt(cx.graph, minus_one)?
                },
            };
            let imag = cx.graph.float(im);
            let imag = mul(cx.graph, &[imag, unit]);
            terms.push(add(cx.graph, &[real, imag]));
        }
    }
    Some(list(cx.graph, &terms))
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
#[allow(clippy::tuple_array_conversions)] // false positive: the tuple is a destructuring of separate values, not a conversion
fn lu(
    cx: &mut Cx<'_>,
    m: &[Vec<NodeId>],
) -> Option<NodeId> {
    let n = m.len();
    if m.iter().any(|r| r.len() != n) {
        return None;
    }
    let mut u: Vec<Vec<NodeId>> = m
        .iter()
        .map(|r| r.iter().map(|&e| cx.simplify(e)).collect())
        .collect();
    let mut l = identity_rows(cx.graph, n);
    let mut perm: Vec<usize> = (0..n).collect();
    let mut field = Terms { cx };
    for col in 0..n {
        let pivot_row = (col..n).find(|&r| !field.is_zero(&u[r][col]))?;
        if pivot_row != col {
            u.swap(pivot_row, col);
            perm.swap(pivot_row, col);
            #[allow(clippy::needless_range_loop)]
            // index is used for more than one array / arithmetic; iterator form would not be clearer
            for c in 0..col {
                let (a, b) = (l[pivot_row][c], l[col][c]);
                l[pivot_row][c] = b;
                l[col][c] = a;
            }
        }
        for r in col + 1..n {
            let factor = field.div(&u[r][col], &u[col][col])?;
            l[r][col] = factor;
            #[allow(clippy::needless_range_loop)]
            // index is used for more than one array / arithmetic; iterator form would not be clearer
            for c in col..n {
                let delta = field.mul(&factor, &u[col][c]);
                u[r][c] = field.sub(&u[r][c], &delta);
            }
        }
    }
    let graph = &mut *field.cx.graph;
    let (zero, one) = (graph.int(0), graph.int(1));
    let p: Vec<Vec<NodeId>> = perm
        .iter()
        .map(|&k| (0..n).map(|j| if j == k { one } else { zero }).collect())
        .collect();
    let (p, l, u) = (
        matrix_term(graph, &p),
        matrix_term(graph, &l),
        matrix_term(graph, &u),
    );
    Some(list(graph, &[p, l, u]))
}

/// `A = Q R` by Gram–Schmidt on the columns (full column rank).
#[allow(clippy::tuple_array_conversions)] // false positive: the tuple is a destructuring of separate values, not a conversion
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
/// Exact singular value decomposition `list(U, S, V)`: `S` holds the
/// square roots of the eigenvalues of `AᵀA` (largest first when they can
/// be ordered), `V` its normalised eigenvectors and `U` the columns
/// `A v / σ` for the non-zero singular values.
fn symbolic_svd(
    cx: &mut Cx<'_>,
    m: &[Vec<NodeId>],
) -> Option<NodeId> {
    let rows = m.len();
    let cols = m.first().map_or(0, Vec::len);
    // AᵀA
    let mut ata = vec![vec![NodeId::NONE; cols]; cols];
    for i in 0..cols {
        for j in 0..cols {
            let terms: Vec<NodeId> = (0..rows)
                .map(|k| cx.graph.node(core::MUL, &[m[k][i], m[k][j]]))
                .collect();
            let sum = cx.graph.node(core::ADD, &terms);
            ata[i][j] = cx.simplify(sum);
        }
    }
    let vectors = eigenvectors(cx, &ata)?;
    let mut triples: Vec<(f64, NodeId, Vec<NodeId>)> = Vec::new();
    for entry in cx.graph.children(vectors).to_vec() {
        let &[value, _, basis] = cx.graph.children(entry) else {
            return None;
        };
        for v in cx.graph.children(basis).to_vec() {
            let comps = cx.graph.children(v).to_vec();
            if comps.len() != cols {
                return None;
            }
            // Normalise.
            let squares: Vec<NodeId> = comps
                .iter()
                .map(|&c| {
                    let two = cx.graph.int(2);
                    cx.graph.node(core::POW, &[c, two])
                })
                .collect();
            let norm2 = cx.graph.node(core::ADD, &squares);
            let half = cx.graph.num(Number::fraction(-1, 2)?);
            let inv_norm = cx.graph.node(core::POW, &[norm2, half]);
            let unit: Vec<NodeId> = comps
                .iter()
                .map(|&c| {
                    let p = cx.graph.node(core::MUL, &[c, inv_norm]);
                    cx.simplify(p)
                })
                .collect();
            let order = cx
                .graph
                .eval(value, &crate::graph::Env::numeric(0.0))
                .unwrap_or(0.0);
            triples.push((order, value, unit));
        }
    }
    triples.sort_by(|a, b| b.0.total_cmp(&a.0));
    let half = cx.graph.num(Number::fraction(1, 2)?);
    let mut sigmas = Vec::new();
    let mut v_columns = Vec::new();
    let mut u_columns = Vec::new();
    for (_, value, unit) in &triples {
        let sigma = cx.graph.node(core::POW, &[*value, half]);
        let sigma = cx.simplify(sigma);
        sigmas.push(sigma);
        v_columns.push(unit.clone());
        if !cx.is_zero(sigma) {
            let minus = cx.graph.int(-1);
            let inv = cx.graph.node(core::POW, &[sigma, minus]);
            let column: Vec<NodeId> = (0..rows)
                .map(|i| {
                    let terms: Vec<NodeId> = (0..cols)
                        .map(|j| cx.graph.node(core::MUL, &[m[i][j], unit[j]]))
                        .collect();
                    let sum = cx.graph.node(core::ADD, &terms);
                    let scaled = cx.graph.node(core::MUL, &[sum, inv]);
                    cx.simplify(scaled)
                })
                .collect();
            u_columns.push(column);
        }
    }
    let transpose = |cols: &[Vec<NodeId>]| -> Vec<Vec<NodeId>> {
        let n = cols.first().map_or(0, Vec::len);
        (0..n)
            .map(|i| cols.iter().map(|c| c[i]).collect())
            .collect()
    };
    let u = matrix_term(cx.graph, &transpose(&u_columns));
    let v = matrix_term(cx.graph, &transpose(&v_columns));
    let s = list(cx.graph, &sigmas);
    Some(list(cx.graph, &[u, s, v]))
}

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
    let FaerDecompositionResult::Svd { u, s, v } = matrix.decompose(FaerDecompositionType::Svd)?
    else {
        return None;
    };
    let to_term = |graph: &mut Graph, mat: &Matrix<f64>| {
        let rows: Vec<Vec<NodeId>> = (0..mat.rows())
            .map(|i| {
                (0..mat.cols())
                    .map(|j| graph.float(*mat.get(i, j)))
                    .collect()
            })
            .collect();
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
    let terms: Vec<NodeId> = a
        .iter()
        .zip(b)
        .map(|(&x, &y)| mul(graph, &[x, y]))
        .collect();
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
        self.compute(cx, &args)
            .map_or(Outcome::Pass, Outcome::Equal)
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
                    acc = if k == 0 {
                        b
                    } else {
                        matmul(cx.graph, &acc, &b)?
                    };
                }
                normalise_entries(cx, &mut acc);
                if flat && (acc.len() == 1 || acc.iter().all(|r| r.len() == 1)) {
                    let items: Vec<NodeId> = acc.into_iter().flatten().collect();
                    return Some(list(cx.graph, &items));
                }
                Some(matrix_term(cx.graph, &acc))
            },
            | Request::Madd
                if operand(cx.graph, arg(0)?).is_some_and(|o| matches!(o, Operand::Vector(_))) =>
            {
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
            | Request::Smul
                if operand(cx.graph, arg(1)?).is_some_and(|o| matches!(o, Operand::Vector(_))) =>
            {
                let c = arg(0)?;
                let v = vector(cx.graph, arg(1)?)?;
                let scaled: Vec<NodeId> = v.iter().map(|&e| mul(cx.graph, &[c, e])).collect();
                Some(list(cx.graph, &scaled))
            },
            | Request::Madd => {
                let mut acc = matrix(cx.graph, arg(0)?)?;
                for &next in args.get(1..)? {
                    let b = matrix(cx.graph, next)?;
                    if b.len() != acc.len() || b.first().map(Vec::len) != acc.first().map(Vec::len)
                    {
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
                let scaled: Vec<Vec<NodeId>> = m
                    .iter()
                    .map(|row| row.iter().map(|&e| mul(cx.graph, &[c, e])).collect())
                    .collect();
                Some(matrix_term(cx.graph, &scaled))
            },
            | Request::Transpose => {
                let m = matrix(cx.graph, arg(0)?)?;
                Some(matrix_term(cx.graph, &transpose(&m)))
            },
            | Request::Identity => {
                let n = usize::try_from(cx.graph.number_of(arg(0)?)?.to_i64()?)
                    .ok()
                    .filter(|&n| n <= 256)?;
                let rows = identity_rows(cx.graph, n);
                Some(matrix_term(cx.graph, &rows))
            },
            | Request::Zeros => {
                let size = |g: &Graph, n: NodeId| {
                    usize::try_from(g.number_of(n)?.to_i64()?)
                        .ok()
                        .filter(|&n| n <= 256)
                };
                let (r, c) = (size(cx.graph, arg(0)?)?, size(cx.graph, arg(1)?)?);
                let zero = cx.graph.int(0);
                Some(matrix_term(cx.graph, &vec![vec![zero; c]; r]))
            },
            | Request::Dims => {
                let m = matrix(cx.graph, arg(0)?)?;
                let rows = cx.graph.int(i64::try_from(m.len()).ok()?);
                let cols = cx
                    .graph
                    .int(i64::try_from(m.first().map_or(0, Vec::len)).ok()?);
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
                let augmented: Vec<Vec<NodeId>> = a
                    .iter()
                    .zip(&b)
                    .map(|(row, &bi)| row.iter().copied().chain([bi]).collect())
                    .collect();
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
                let Some(values) = eigenvalues(cx, &m) else {
                    return numeric_eigenvalues(cx, &m);
                };
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
                svd(cx.graph, &m).or_else(|| symbolic_svd(cx, &m))
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
                let rows: Vec<Vec<NodeId>> = a
                    .iter()
                    .map(|&x| b.iter().map(|&y| mul(cx.graph, &[x, y])).collect())
                    .collect();
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
                let (integrand, curve_node, t, a, b) =
                    (arg(0)?, arg(1)?, arg(2)?, arg(3)?, arg(4)?);
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
                    let along: Vec<NodeId> =
                        field.iter().map(|&f| substitute(cx.graph, f)).collect();
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
                let normal = Self {
                    op: self.op,
                    request: Request::Cross,
                }
                .compute(cx, &[ru_node, rv_node])?;
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
        ["x", "y", "z"]
            .iter()
            .take(n)
            .map(|name| graph.sym(name))
            .collect()
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
        assert_eq!(
            run("matmul(list(list(1, 2), list(3, 4)), list(list(0, 1), list(1, 0)))"),
            "list(list(2, 1), list(4, 3))"
        );
        assert_eq!(
            run("madd(list(list(1, 2)), list(list(a, b)))"),
            "list(list(a + 1, b + 2))"
        );
        assert_eq!(run("smul(2, list(list(1, x)))"), "list(list(2, 2*x))");
        assert_eq!(
            run("transpose(list(list(1, 2, 3), list(4, 5, 6)))"),
            "list(list(1, 4), list(2, 5), list(3, 6))"
        );
        assert_eq!(run("identity(2)"), "list(list(1, 0), list(0, 1))");
        assert_eq!(run("trace(list(list(a, b), list(c, d)))"), "a + d");
        let (text, reduced) = reduce_with(
            &[linalg()],
            "matmul(list(list(1, 2)), list(list(1, 2)))",
            &[],
        );
        assert!(!reduced, "shapes do not match: {text}");
    }

    #[test]
    fn determinants_inverses_and_rank() {
        assert_eq!(run("det(list(list(1, 2), list(3, 4)))"), "-2");
        assert_eq!(run("det(list(list(a, b), list(c, d)))"), "a*d - b*c");
        assert_eq!(
            run("det(list(list(2, 0, 1), list(1, 3, 2), list(1, 1, 2)))"),
            "6"
        );
        assert_eq!(run("det(list(list(1, 2), list(2, 4)))"), "0");
        assert_eq!(
            run("inverse(list(list(1, 2), list(3, 4)))"),
            "list(list(-2, 1), list(3/2, -1/2))"
        );
        assert_eq!(
            run("inverse(list(list(a, 0), list(0, b)))"),
            "list(list(1/a, 0), list(0, 1/b))"
        );
        let (text, reduced) =
            reduce_with(&[linalg()], "inverse(list(list(1, 2), list(2, 4)))", &[]);
        assert!(!reduced, "a singular matrix has no inverse: {text}");
        assert_eq!(
            run("rank(list(list(1, 2, 3), list(2, 4, 6), list(1, 0, 1)))"),
            "2"
        );
        assert_eq!(
            run("rref(list(list(1, 2, 3), list(4, 5, 6)))"),
            "list(list(1, 0, -1), list(0, 1, 2))"
        );
        assert_eq!(
            run("nullspace(list(list(1, 2, 3), list(4, 5, 6)))"),
            "list(list(1, -2, 1))"
        );
    }

    #[test]
    fn inverse_times_matrix_is_identity() {
        let a = "list(list(2, 1, 1), list(1, 3, 2), list(1, 0, 0))";
        assert_eq!(
            run(&format!("matmul(inverse({a}), {a})")),
            "list(list(1, 0, 0), list(0, 1, 0), list(0, 0, 1))"
        );
        let s = "list(list(a, 1), list(1, a))";
        assert_eq!(
            run(&format!("matmul({s}, inverse({s}))")),
            "list(list(1, 0), list(0, 1))"
        );
    }

    #[test]
    fn linear_systems() {
        assert_eq!(
            run("linsolve(list(list(2, 1), list(1, 3)), list(3, 5))"),
            "list(4/5, 7/5)"
        );
        assert_eq!(
            run("linsolve(list(list(1, 1), list(1, 1)), list(1, 2))"),
            "list()"
        );
        assert_eq!(
            run("linsolve(list(list(a, 0), list(0, b)), list(1, 1))"),
            "list(1/a, 1/b)"
        );
        assert_eq!(
            run("linsolve(list(list(1, 1, 1)), list(3))"),
            "list(3 - t1 - t2, t1, t2)"
        );
        assert_eq!(
            run("dims(list(list(1, 2, 3), list(4, 5, 6)))"),
            "list(2, 3)"
        );
    }

    #[test]
    fn eigen() {
        assert_eq!(
            run("charpoly(list(list(1, 2), list(3, 4)), l)"),
            "l^2 - 5*l - 2"
        );
        assert_eq!(run("eigenvals(list(list(2, 0), list(0, 3)))"), "list(2, 3)");
        assert_eq!(run("eigenvals(list(list(2, 1), list(1, 2)))"), "list(1, 3)");
        assert_eq!(run("eigenvals(list(list(1, 1), list(0, 1)))"), "list(1, 1)");
        // x^5 - x - 1.5 has no closed-form roots: the companion matrix with
        // a float entry falls back to the QR algorithm.
        let quintic = run(
            "eigenvals(list(list(0, 0, 0, 0, 1.5), list(1, 0, 0, 0, 1), list(0, 1, 0, 0, 0), list(0, 0, 1, 0, 0), list(0, 0, 0, 1, 0)))",
        );
        assert!(
            quintic.starts_with("list(") && quintic.contains("1.22171074511"),
            "{quintic}"
        );
        assert_eq!(quintic.matches("(-1)^(1/2)").count(), 4, "{quintic}");
        assert_eq!(
            run("eigenvects(list(list(2, 1), list(1, 2)))"),
            "list(list(1, 1, list(list(-1, 1))), list(3, 1, list(list(1, 1))))"
        );
        // Symbolic entries.
        assert_eq!(
            run("charpoly(list(list(a, b), list(b, a)), l)"),
            "a^2 - 2*a*l - b^2 + l^2"
        );
    }

    #[test]
    fn decompositions() {
        // P A = L U must hold.
        let a = "list(list(0, 2, 1), list(1, 1, 0), list(2, 1, 3))";
        let lu = run(&format!("lu({a})"));
        let check = run(
            &format!("madd(matmul(item0, {a}), smul(-1, matmul(item1, item2)))")
                .replace("item0", &part(&lu, 0))
                .replace("item1", &part(&lu, 1))
                .replace("item2", &part(&lu, 2)),
        );
        assert_eq!(
            check, "list(list(0, 0, 0), list(0, 0, 0), list(0, 0, 0))",
            "{lu}"
        );
        // A = Q R with orthonormal Q.
        let qr = run("qr(list(list(3, 1), list(4, 2)))");
        let product = run(&format!("matmul({}, {})", part(&qr, 0), part(&qr, 1)));
        assert_eq!(product, "list(list(3, 1), list(4, 2))", "{qr}");
        // Numeric SVD: singular values of diag(3, -2).
        let svd = run("svd(list(list(3.0, 0.0), list(0.0, -2.0)))");
        let values = part(&svd, 1);
        assert!(values.contains('3') && values.contains('2'), "{svd}");
        // Exact input: the symbolic SVD keeps exact singular values.
        let exact = run("svd(list(list(3, 0), list(0, -2)))");
        let values = part(&exact, 1);
        assert!(
            !exact.contains('.') && values.contains('3') && values.contains('2'),
            "{exact}"
        );
        assert_eq!(values, "list(3, 2)", "{exact}");
    }

    /// The `index`-th top-level item of a `list(...)` text.
    fn part(
        text: &str,
        index: usize,
    ) -> String {
        let inner = text
            .strip_prefix("list(")
            .and_then(|s| s.strip_suffix(')'))
            .unwrap_or("");
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
        assert_eq!(
            run("outer(list(1, 2), list(3, 4))"),
            "list(list(3, 4), list(6, 8))"
        );
    }

    #[test]
    fn vector_calculus() {
        assert_eq!(run("grad(x^2*y + z, list(x, y, z))"), "list(2*x*y, x^2, 1)");
        assert_eq!(run("div(list(x, y^2, z*x), list(x, y, z))"), "x + 2*y + 1");
        assert_eq!(run("curl(list(-y, x, 0), list(x, y, z))"), "list(0, 0, 2)");
        assert_eq!(run("laplacian(x^2 + y^2 + z^2, list(x, y, z))"), "6");
        assert_eq!(
            run("jacobian(list(x*y, x + y), list(x, y))"),
            "list(list(y, x), list(1, 1))"
        );
        assert_eq!(
            run("hessian(x^2*y, list(x, y))"),
            "list(list(2*y, 2*x), list(2*x, 0))"
        );
        assert_eq!(
            run("directional(x^2 + y^2, list(x, y), list(3, 4))"),
            "6/5*x + 8/5*y"
        );
        // div curl = 0 and curl grad = 0 for any field.
        assert_eq!(
            run("div(curl(list(x*y*z, sin(x)*y, exp(z)*x), list(x, y, z)), list(x, y, z))"),
            "0"
        );
        assert_eq!(
            run("curl(grad(x^2*sin(y)*z, list(x, y, z)), list(x, y, z))"),
            "list(0, 0, 0)"
        );
    }

    #[test]
    fn integrals_over_curves_surfaces_and_volumes() {
        // Arc length of the unit circle.
        let (value, _) = numeric(
            &[linalg()],
            "line_integral(1, list(cos(t), sin(t)), t, 0, 2*pi)",
            &[],
            1e-10,
        );
        assert!((value - 2.0 * std::f64::consts::PI).abs() < 1e-9, "{value}");
        // Work of F = (-y, x) around the unit circle: 2 pi.
        assert_eq!(
            run("line_integral_vec(list(-y, x), list(cos(t), sin(t)), t, 0, 2*pi)"),
            "2*pi"
        );
        // Area of the unit sphere.
        let (value, _) = numeric(
            &[linalg()],
            "surface_integral(1, list(sin(u)*cos(v), sin(u)*sin(v), cos(u)), u, v, 0, pi, 0, 2*pi)",
            &[],
            1e-8,
        );
        assert!((value - 4.0 * std::f64::consts::PI).abs() < 1e-6, "{value}");
        assert_eq!(
            run("volume_integral(x*y*z, list(x, y, z), list(list(0, 1), list(0, 2), list(0, 3)))"),
            "9/2"
        );
        let _ = eval;
    }
}

// Exact matrix functions and factorisations: Jordan normal form, matrix
// exponentials, powers and arbitrary functions of a matrix, minimal
// polynomial, Cayley–Hamilton, pseudo-inverse and least squares, LDL and
// Cholesky, Smith and Hermite normal forms over the integers, spaces of a
// matrix, Kronecker products and special matrices.
//
// Entries may be exact rationals or symbolic terms; the algorithms use the
// same zero-testing arithmetic as `rref` and `det`, so the results are
// exact for generic symbolic entries.
//
// | operator | value |
// |---|---|
// | `kron(A, B)` | the Kronecker product |
// | `diag(list(d1, ...))`, `vandermonde(list(x1, ...))`, `hilbert(n)`, `companion(p, x)` | special matrices (`companion` of a monic-normalised polynomial) |
// | `adjugate(A)`, `cofactors(A)` | the adjugate (classical adjoint) and the matrix of cofactors |
// | `colspace(A)`, `rowspace(A)`, `left_nullspace(A)` | bases, as lists of vectors |
// | `pinv(A)` | the Moore–Penrose pseudo-inverse, from a full-rank factorisation `A = C R` |
// | `lstsq(A, b)` | the minimum-norm least-squares solution `pinv(A) b` |
// | `minpoly(A, l)` | the minimal polynomial in `l` (from the first linear dependence among the powers of `A`) |
// | `matpoly(p, x, A)` | the polynomial `p(x)` evaluated at the matrix `A` (Horner) |
// | `cayley_hamilton(A)` | the characteristic polynomial evaluated at `A`: the zero matrix |
// | `jordan(A)` | `list(P, J)` with `A = P J P^-1` and `J` in Jordan normal form (eigenvalues that have a closed form) |
// | `jordan_blocks(A)` | `list(list(eigenvalue, size), ...)` |
// | `matfun(f, x, A)` | `f(A)` for an expression `f` in `x`: `P f(J) P^-1` with the derivative entries `f^(k)(λ)/k!` in each Jordan block |
// | `matexp(A, t)` | `exp(A t)` (the matrix exponential), exact via the Jordan form |
// | `mpow(A, n)` | the matrix power: repeated squaring for an integer `n` (negative via the inverse), the Jordan form for a symbolic `n` |
// | `matsqrt(A)` | a square root through `matfun(sqrt(x), x, A)` |
// | `spectral(A)` | `list(list(eigenvalue, projector), ...)` for a diagonalisable matrix: `A = Σ λ P_λ` |
// | `ldl(A)` | `list(L, D)` with `A = L diag(D) Lᵀ` for a symmetric matrix with non-vanishing pivots |
// | `cholesky(A)` | `L` with `A = L Lᵀ` (positive pivots; for symbolic pivots they are assumed positive) |
// | `orthogonalize(vs)`, `orthonormalize(vs)` | Gram–Schmidt on a list of vectors (dependent vectors dropped) |
// | `smith(A)` | `list(U, D, V)` with `U A V = D` the Smith normal form of an integer matrix, `U` and `V` unimodular |
// | `invariant_factors(A)` | the non-zero diagonal of the Smith normal form |
// | `hermite_form(A)` | `list(H, U)` with `U A = H` the row-style Hermite normal form of an integer matrix |
// | `is_symmetric(A)`, `is_orthogonal(A)`, `is_positive_definite(A)` | truth values (the last for exact rational matrices, by Sylvester's criterion) |

use num_traits::Signed;


type Mat = Vec<Vec<NodeId>>;

/// Largest matrix the Jordan machinery accepts.
const MAX_JORDAN: usize = 12;

fn square(m: &[Vec<NodeId>]) -> bool {
    m.iter().all(|r| r.len() == m.len()) && !m.is_empty()
}

fn zero_matrix(
    graph: &mut Graph,
    rows: usize,
    cols: usize,
) -> Mat {
    let zero = graph.int(0);
    vec![vec![zero; cols]; rows]
}

/// `a b` with simplified entries.
fn product(
    cx: &mut Cx<'_>,
    a: &[Vec<NodeId>],
    b: &[Vec<NodeId>],
) -> Option<Mat> {
    let mut out = matmul(cx.graph, a, b)?;
    for row in &mut out {
        for entry in row.iter_mut() {
            *entry = normal(cx, *entry);
        }
    }
    Some(out)
}

fn sum(
    cx: &mut Cx<'_>,
    a: &[Vec<NodeId>],
    b: &[Vec<NodeId>],
    sign: i64,
) -> Option<Mat> {
    if a.len() != b.len() {
        return None;
    }
    let factor = cx.graph.int(sign);
    let mut out = Vec::with_capacity(a.len());
    for (ra, rb) in a.iter().zip(b) {
        if ra.len() != rb.len() {
            return None;
        }
        let mut row = Vec::with_capacity(ra.len());
        for (&x, &y) in ra.iter().zip(rb) {
            let scaled = mul(cx.graph, &[factor, y]);
            let total = add(cx.graph, &[x, scaled]);
            row.push(normal(cx, total));
        }
        out.push(row);
    }
    Some(out)
}

fn scale(
    cx: &mut Cx<'_>,
    factor: NodeId,
    a: &[Vec<NodeId>],
) -> Mat {
    a.iter()
        .map(|row| {
            row.iter()
                .map(|&e| {
                    let p = mul(cx.graph, &[factor, e]);
                    normal(cx, p)
                })
                .collect()
        })
        .collect()
}

fn shifted(
    cx: &mut Cx<'_>,
    m: &[Vec<NodeId>],
    value: NodeId,
) -> Mat {
    let minus = neg(cx.graph, value);
    let mut out = m.to_vec();
    for (i, row) in out.iter_mut().enumerate() {
        if let Some(entry) = row.get_mut(i) {
            let total = add(cx.graph, &[*entry, minus]);
            *entry = normal(cx, total);
        }
    }
    out
}

fn rank_of(
    cx: &mut Cx<'_>,
    rows: &[Vec<NodeId>],
) -> Option<usize> {
    if rows.is_empty() {
        return Some(0);
    }
    let (_, pivots, _) = reduce(cx, rows)?;
    Some(pivots.len())
}

fn apply(
    cx: &mut Cx<'_>,
    m: &[Vec<NodeId>],
    v: &[NodeId],
) -> Vec<NodeId> {
    m.iter()
        .map(|row| {
            let terms: Vec<NodeId> = row
                .iter()
                .zip(v)
                .map(|(&a, &b)| mul(cx.graph, &[a, b]))
                .collect();
            let total = add(cx.graph, &terms);
            normal(cx, total)
        })
        .collect()
}

// ----------------------------------------------------------------------
// Jordan normal form
// ----------------------------------------------------------------------

struct Block {
    value: NodeId,
    /// Chain vectors `p_1, ..., p_s` with `(A - λ) p_1 = 0`, `(A - λ) p_k = p_{k-1}`.
    chain: Vec<Vec<NodeId>>,
}

fn jordan_blocks(
    cx: &mut Cx<'_>,
    m: &[Vec<NodeId>],
) -> Option<Vec<Block>> {
    let n = m.len();
    if !square(m) || n > MAX_JORDAN {
        return None;
    }
    let values = eigenvalues(cx, m)?;
    let mut blocks = Vec::new();
    for (value, multiplicity) in values {
        let shift = shifted(cx, m, value);
        // Kernels of N^k until the dimension reaches the multiplicity.
        let mut kernels: Vec<Vec<Vec<NodeId>>> = Vec::new();
        let mut power = identity_rows(cx.graph, n);
        let mut dims = vec![0_usize];
        loop {
            power = product(cx, &shift, &power)?;
            let basis = nullspace(cx, &power)?;
            let dim = basis.len();
            kernels.push(basis);
            dims.push(dim);
            if dim >= multiplicity || kernels.len() > multiplicity {
                break;
            }
        }
        if dims.last().copied() != Some(multiplicity) {
            return None;
        }
        let top_size = kernels.len();
        // at_least[k] = number of blocks of size >= k (1-based; at_least[top+1] = 0).
        let mut at_least = vec![0_usize; top_size + 2];
        for k in 1..=top_size {
            at_least[k] = dims[k] - dims[k - 1];
        }
        let mut tops: Vec<(usize, Vec<NodeId>)> = Vec::new();
        for size in (1..=top_size).rev() {
            let exact = at_least[size] - at_least[size + 1];
            if exact == 0 {
                continue;
            }
            // The span already accounted for inside ker N^size.
            let mut current: Vec<Vec<NodeId>> = if size >= 2 {
                kernels[size - 2].clone()
            } else {
                Vec::new()
            };
            for (bigger, top) in &tops {
                let mut image = top.clone();
                for _ in 0..(bigger - size) {
                    image = apply(cx, &shift, &image);
                }
                current.push(image);
            }
            let mut rank = rank_of(cx, &current)?;
            let mut chosen = 0;
            for candidate in kernels[size - 1].clone() {
                if chosen == exact {
                    break;
                }
                current.push(candidate.clone());
                let new_rank = rank_of(cx, &current)?;
                if new_rank > rank {
                    rank = new_rank;
                    tops.push((size, candidate));
                    chosen += 1;
                } else {
                    current.pop();
                }
            }
            if chosen != exact {
                return None;
            }
        }
        for (size, top) in tops {
            let mut chain = Vec::with_capacity(size);
            let mut v = top;
            for _ in 0..size {
                chain.push(v.clone());
                v = apply(cx, &shift, &v);
            }
            chain.reverse();
            blocks.push(Block { value, chain });
        }
    }
    Some(blocks)
}

/// `(P, blocks)`: the columns of `P` are the chains, block after block.
fn jordan_basis(
    cx: &mut Cx<'_>,
    m: &[Vec<NodeId>],
) -> Option<(Mat, Vec<Block>)> {
    let blocks = jordan_blocks(cx, m)?;
    let n = m.len();
    let mut columns: Vec<Vec<NodeId>> = Vec::with_capacity(n);
    for block in &blocks {
        columns.extend(block.chain.iter().cloned());
    }
    if columns.len() != n {
        return None;
    }
    Some((transpose(&columns), blocks))
}

fn jordan_matrix(
    cx: &mut Cx<'_>,
    n: usize,
    blocks: &[Block],
) -> Mat {
    let mut j = zero_matrix(cx.graph, n, n);
    let one = cx.graph.int(1);
    let mut at = 0;
    for block in blocks {
        let size = block.chain.len();
        for k in 0..size {
            j[at + k][at + k] = block.value;
            if k + 1 < size {
                j[at + k][at + k + 1] = one;
            }
        }
        at += size;
    }
    j
}

fn jordan(
    cx: &mut Cx<'_>,
    m: &[Vec<NodeId>],
) -> Option<NodeId> {
    let (p, blocks) = jordan_basis(cx, m)?;
    let j = jordan_matrix(cx, m.len(), &blocks);
    let (p, j) = (matrix_term(cx.graph, &p), matrix_term(cx.graph, &j));
    Some(list(cx.graph, &[p, j]))
}

fn jordan_block_list(
    cx: &mut Cx<'_>,
    m: &[Vec<NodeId>],
) -> Option<NodeId> {
    let blocks = jordan_blocks(cx, m)?;
    let mut items = Vec::new();
    for block in &blocks {
        let size = cx.graph.int(i64::try_from(block.chain.len()).ok()?);
        items.push(list(cx.graph, &[block.value, size]));
    }
    Some(list(cx.graph, &items))
}

/// `f(A)` where `f` is an expression in the symbol `x`.
fn matrix_function(
    cx: &mut Cx<'_>,
    f: NodeId,
    x: NodeId,
    m: &[Vec<NodeId>],
) -> Option<Mat> {
    let n = m.len();
    let (p, blocks) = jordan_basis(cx, m)?;
    let inverse_p = inverse(cx, &p)?;
    // Derivatives f^(k)(x) / k! as needed.
    let largest = blocks.iter().map(|b| b.chain.len()).max()?;
    let mut derivatives = vec![f];
    for _ in 1..largest {
        let last = *derivatives.last()?;
        let d = derivative(cx.graph, last, x)?;
        derivatives.push(cx.simplify(d));
    }
    let mut middle = zero_matrix(cx.graph, n, n);
    let mut at = 0;
    for block in &blocks {
        let size = block.chain.len();
        let mut factorial = BigInt::one();
        let mut entries = Vec::with_capacity(size);
        for (k, &d) in derivatives.iter().enumerate().take(size) {
            if k > 0 {
                factorial *= k;
            }
            let value = cx.graph.substitute(d, x, block.value);
            let divisor = cx.graph.num(Number::rat(BigRational::new(
                BigInt::one(),
                factorial.clone(),
            )));
            let term = mul(cx.graph, &[divisor, value]);
            entries.push(normal(cx, term));
        }
        for row in 0..size {
            middle[at + row][at + row..at + size].copy_from_slice(&entries[..size - row]);
        }
        at += size;
    }
    let left = product(cx, &p, &middle)?;
    product(cx, &left, &inverse_p)
}

fn fresh(cx: &mut Cx<'_>) -> NodeId {
    let s = cx.graph.interner_mut().fresh_symbol("x");
    cx.graph.symbol_node(s)
}

fn matrix_exponential(
    cx: &mut Cx<'_>,
    m: &[Vec<NodeId>],
    t: NodeId,
) -> Option<Mat> {
    let x = fresh(cx);
    let argument = mul(cx.graph, &[x, t]);
    let f = request(cx.graph, "exp", &[argument])?;
    matrix_function(cx, f, x, m)
}

fn integer_power(
    cx: &mut Cx<'_>,
    m: &[Vec<NodeId>],
    exponent: u64,
) -> Option<Mat> {
    let mut result = identity_rows(cx.graph, m.len());
    let mut base = m.to_vec();
    let mut e = exponent;
    while e > 0 {
        if e & 1 == 1 {
            result = product(cx, &result, &base)?;
        }
        e >>= 1;
        if e > 0 {
            base = product(cx, &base, &base)?;
        }
    }
    Some(result)
}

fn matrix_power(
    cx: &mut Cx<'_>,
    m: &[Vec<NodeId>],
    exponent: NodeId,
) -> Option<Mat> {
    if let Some(e) = cx.graph.number_of(exponent).and_then(Number::to_i64) {
        let magnitude = Some(e.unsigned_abs()).filter(|&v| v <= 100_000)?;
        if e >= 0 {
            return integer_power(cx, m, magnitude);
        }
        let inv = inverse(cx, m)?;
        return integer_power(cx, &inv, magnitude);
    }
    let x = fresh(cx);
    let f = cx.graph.node(core::POW, &[x, exponent]);
    matrix_function(cx, f, x, m)
}

fn spectral(
    cx: &mut Cx<'_>,
    m: &[Vec<NodeId>],
) -> Option<NodeId> {
    let n = m.len();
    let values = eigenvalues(cx, m)?;
    let mut items = Vec::new();
    for (value, multiplicity) in values {
        let shift = shifted(cx, m, value);
        let basis = nullspace(cx, &shift)?;
        if basis.len() != multiplicity {
            // Defective: no spectral decomposition into projectors.
            return None;
        }
        // P = V (Vᵀ V)⁻¹ Vᵀ with the eigenvectors as columns of V.
        let v = transpose(&basis);
        let vt = basis.clone();
        let gram = product(cx, &vt, &v)?;
        let gram_inverse = inverse(cx, &gram)?;
        let left = product(cx, &v, &gram_inverse)?;
        let projector = product(cx, &left, &vt)?;
        debug_assert_eq!(projector.len(), n);
        let projector = matrix_term(cx.graph, &projector);
        items.push(list(cx.graph, &[value, projector]));
    }
    Some(list(cx.graph, &items))
}

// ----------------------------------------------------------------------
// Polynomials of matrices
// ----------------------------------------------------------------------

fn flatten(m: &[Vec<NodeId>]) -> Vec<NodeId> {
    m.iter().flatten().copied().collect()
}

fn minimal_polynomial(
    cx: &mut Cx<'_>,
    m: &[Vec<NodeId>],
    l: NodeId,
) -> Option<NodeId> {
    let n = m.len();
    if !square(m) {
        return None;
    }
    let mut powers: Vec<Mat> = vec![identity_rows(cx.graph, n)];
    for k in 1..=n {
        let next = product(cx, m, &powers[k - 1])?;
        powers.push(next);
        // Columns are the flattened powers I, A, ..., A^k.
        let columns: Vec<Vec<NodeId>> = powers.iter().map(|p| flatten(p)).collect();
        let rows = transpose(&columns);
        let basis = nullspace(cx, &rows)?;
        if let Some(dependency) = basis.first() {
            let lead = *dependency.last()?;
            let mut terms = Vec::new();
            for (power, &c) in dependency.iter().enumerate() {
                let inverse_lead = reciprocal(cx.graph, lead);
                let coefficient = mul(cx.graph, &[c, inverse_lead]);
                let coefficient = normal(cx, coefficient);
                let exponent = cx.graph.int(i64::try_from(power).ok()?);
                let raised = cx.graph.node(core::POW, &[l, exponent]);
                terms.push(mul(cx.graph, &[coefficient, raised]));
            }
            let total = add(cx.graph, &terms);
            return Some(cx.simplify(total));
        }
    }
    None
}

/// `p(A)` for a polynomial `p` in `x`, by Horner's rule on the
/// coefficients of `p`.
fn matrix_polynomial(
    cx: &mut Cx<'_>,
    p: NodeId,
    x: NodeId,
    m: &[Vec<NodeId>],
) -> Option<Mat> {
    let n = m.len();
    if !square(m) {
        return None;
    }
    let mut gens = crate::rules::poly::repr::Gens::default();
    let gx = gens.index(cx.graph, x);
    let poly = crate::rules::poly::repr::from_term(
        cx.graph,
        &mut gens,
        p,
        crate::rules::poly::repr::Limits::default(),
    )?;
    // Coefficients by repeated differentiation at zero: c_k = p^(k)(0)/k!.
    let mut coefficients = Vec::new();
    let mut current = poly;
    let zero = cx.graph.int(0);
    let mut factorial = BigInt::one();
    let mut order = 0_u32;
    while !current.is_zero() {
        if order > 0 {
            factorial *= order;
        }
        let term = crate::rules::poly::repr::to_term(cx.graph, &gens, &current);
        let at_zero = cx.graph.substitute(term, x, zero);
        let divisor = cx.graph.num(Number::rat(BigRational::new(
            BigInt::one(),
            factorial.clone(),
        )));
        let value = mul(cx.graph, &[divisor, at_zero]);
        coefficients.push(cx.simplify(value));
        current = current.derivative(gx);
        order += 1;
    }
    let identity = identity_rows(cx.graph, n);
    let mut acc = zero_matrix(cx.graph, n, n);
    for &c in coefficients.iter().rev() {
        acc = product(cx, &acc, m)?;
        let c_identity = scale(cx, c, &identity);
        acc = sum(cx, &acc, &c_identity, 1)?;
    }
    Some(acc)
}

// ----------------------------------------------------------------------
// Pseudo-inverse and friends
// ----------------------------------------------------------------------

fn pseudo_inverse(
    cx: &mut Cx<'_>,
    m: &[Vec<NodeId>],
) -> Option<Mat> {
    let rows = m.len();
    let cols = m.first().map_or(0, Vec::len);
    let (reduced, pivots, _) = reduce(cx, m)?;
    if pivots.is_empty() {
        return Some(zero_matrix(cx.graph, cols, rows));
    }
    // A = C R with C the pivot columns and R the non-zero rows of rref(A).
    let c: Mat = (0..rows)
        .map(|i| pivots.iter().map(|&p| m[i][p]).collect())
        .collect();
    let r: Mat = reduced.into_iter().take(pivots.len()).collect();
    let (ct, rt) = (transpose(&c), transpose(&r));
    let rrt = product(cx, &r, &rt)?;
    let ctc = product(cx, &ct, &c)?;
    let (rrt_inverse, ctc_inverse) = (inverse(cx, &rrt)?, inverse(cx, &ctc)?);
    let left = product(cx, &rt, &rrt_inverse)?;
    let middle = product(cx, &left, &ctc_inverse)?;
    product(cx, &middle, &ct)
}

fn colspace(
    cx: &mut Cx<'_>,
    m: &[Vec<NodeId>],
) -> Option<Vec<Vec<NodeId>>> {
    let (_, pivots, _) = reduce(cx, m)?;
    Some(
        pivots
            .iter()
            .map(|&p| m.iter().map(|row| row[p]).collect())
            .collect(),
    )
}

fn rowspace(
    cx: &mut Cx<'_>,
    m: &[Vec<NodeId>],
) -> Option<Vec<Vec<NodeId>>> {
    let (reduced, pivots, _) = reduce(cx, m)?;
    Some(reduced.into_iter().take(pivots.len()).collect())
}

fn ldl(
    cx: &mut Cx<'_>,
    m: &[Vec<NodeId>],
) -> Option<(Mat, Vec<NodeId>)> {
    let n = m.len();
    if !square(m) {
        return None;
    }
    let mut l = identity_rows(cx.graph, n);
    let mut d: Vec<NodeId> = Vec::with_capacity(n);
    for j in 0..n {
        // d_j = a_jj - Σ_k l_jk² d_k
        let mut value = m[j][j];
        for (k, &dk) in d.iter().enumerate() {
            let square = mul(cx.graph, &[l[j][k], l[j][k], dk]);
            let negated = neg(cx.graph, square);
            value = add(cx.graph, &[value, negated]);
        }
        let value = normal(cx, value);
        if (Terms { cx }).is_zero(&value) {
            return None;
        }
        d.push(value);
        for i in j + 1..n {
            let mut entry = m[i][j];
            for (k, &dk) in d.iter().enumerate().take(j) {
                let t = mul(cx.graph, &[l[i][k], l[j][k], dk]);
                let negated = neg(cx.graph, t);
                entry = add(cx.graph, &[entry, negated]);
            }
            let inverse_d = reciprocal(cx.graph, value);
            let scaled = mul(cx.graph, &[entry, inverse_d]);
            l[i][j] = normal(cx, scaled);
        }
    }
    Some((l, d))
}

fn is_symmetric(
    cx: &mut Cx<'_>,
    m: &[Vec<NodeId>],
) -> bool {
    if !square(m) {
        return false;
    }
    for (i, row) in m.iter().enumerate() {
        for (j, &entry) in row.iter().enumerate().skip(i + 1) {
            let negated = neg(cx.graph, m[j][i]);
            let difference = add(cx.graph, &[entry, negated]);
            if !cx.is_zero(difference) {
                return false;
            }
        }
    }
    true
}

fn cholesky(
    cx: &mut Cx<'_>,
    m: &[Vec<NodeId>],
) -> Option<Mat> {
    if !is_symmetric(cx, m) {
        return None;
    }
    let (mut l, d) = ldl(cx, m)?;
    for (j, &dj) in d.iter().enumerate() {
        if let Some(q) = cx.graph.number_of(dj).and_then(Number::to_rational)
            && !q.is_positive()
        {
            return None;
        }
        let root = sqrt(cx.graph, dj)?;
        let root = cx.simplify(root);
        for row in &mut l {
            let scaled = mul(cx.graph, &[row[j], root]);
            row[j] = normal(cx, scaled);
        }
    }
    Some(l)
}

/// Gram–Schmidt on vectors (orthogonal, or orthonormal), dropping the
/// dependent ones.
fn gram_schmidt(
    cx: &mut Cx<'_>,
    vectors: &[Vec<NodeId>],
    normalise: bool,
) -> Option<Vec<Vec<NodeId>>> {
    let mut basis: Vec<Vec<NodeId>> = Vec::new();
    let mut squares: Vec<NodeId> = Vec::new();
    for v in vectors {
        let mut w = v.clone();
        for (b, &bb) in basis.iter().zip(&squares) {
            let numerator = dot_ext(cx, v, b);
            let inverse_bb = reciprocal(cx.graph, bb);
            let coefficient = mul(cx.graph, &[numerator, inverse_bb]);
            let coefficient = normal(cx, coefficient);
            for (wi, &bi) in w.iter_mut().zip(b) {
                let t = mul(cx.graph, &[coefficient, bi]);
                let negated = neg(cx.graph, t);
                let total = add(cx.graph, &[*wi, negated]);
                *wi = normal(cx, total);
            }
        }
        let length_squared = dot_ext(cx, &w, &w);
        if cx.is_zero(length_squared) {
            continue;
        }
        if normalise {
            let root = sqrt(cx.graph, length_squared)?;
            let root = cx.simplify(root);
            let inverse_root = reciprocal(cx.graph, root);
            let unit: Vec<NodeId> = w
                .iter()
                .map(|&wi| {
                    let p = mul(cx.graph, &[wi, inverse_root]);
                    normal(cx, p)
                })
                .collect();
            let one = cx.graph.int(1);
            basis.push(unit);
            squares.push(one);
        } else {
            basis.push(w);
            squares.push(length_squared);
        }
    }
    Some(basis)
}

fn dot_ext(
    cx: &mut Cx<'_>,
    a: &[NodeId],
    b: &[NodeId],
) -> NodeId {
    let terms: Vec<NodeId> = a
        .iter()
        .zip(b)
        .map(|(&x, &y)| mul(cx.graph, &[x, y]))
        .collect();
    let total = add(cx.graph, &terms);
    normal(cx, total)
}

// ----------------------------------------------------------------------
// Integer normal forms
// ----------------------------------------------------------------------

type Ints = Vec<Vec<BigInt>>;

fn integer_matrix(
    graph: &mut Graph,
    node: NodeId,
) -> Option<Ints> {
    let m = matrix(graph, node)?;
    m.iter()
        .map(|row| {
            row.iter()
                .map(|&e| {
                    match graph.number_of(e)? {
                        | Number::Int(v) => Some(v.clone()),
                        | _ => None,
                    }
                })
                .collect()
        })
        .collect()
}

fn integer_identity(n: usize) -> Ints {
    (0..n)
        .map(|i| (0..n).map(|j| BigInt::from(i32::from(i == j))).collect())
        .collect()
}

fn integers_to_term(
    graph: &mut Graph,
    m: &Ints,
) -> NodeId {
    let rows: Mat = m
        .iter()
        .map(|r| {
            r.iter()
                .map(|v| graph.num(Number::Int(v.clone())))
                .collect()
        })
        .collect();
    matrix_term(graph, &rows)
}

fn row_sub(
    m: &mut Ints,
    target: usize,
    source: usize,
    factor: &BigInt,
) {
    let source_row = m[source].clone();
    for (value, s) in m[target].iter_mut().zip(&source_row) {
        *value -= factor * s;
    }
}

fn col_sub(
    m: &mut Ints,
    target: usize,
    source: usize,
    factor: &BigInt,
) {
    for row in m.iter_mut() {
        let s = row[source].clone();
        row[target] -= factor * s;
    }
}

fn floor_div(
    a: &BigInt,
    b: &BigInt,
) -> BigInt {
    num_integer::Integer::div_floor(a, b)
}

/// Smith normal form `U A V = D` with unimodular `U`, `V`.
fn smith(a: &Ints) -> (Ints, Ints, Ints) {
    let rows = a.len();
    let cols = a.first().map_or(0, Vec::len);
    let mut d = a.clone();
    let mut u = integer_identity(rows);
    let mut v = integer_identity(cols);
    let limit = rows.min(cols);
    for t in 0..limit {
        loop {
            // Smallest non-zero entry in the trailing block.
            let mut best: Option<(usize, usize)> = None;
            for i in t..rows {
                for j in t..cols {
                    if !d[i][j].is_zero()
                        && best.is_none_or(|(bi, bj)| d[i][j].abs() < d[bi][bj].abs())
                    {
                        best = Some((i, j));
                    }
                }
            }
            let Some((pi, pj)) = best else { break };
            d.swap(t, pi);
            u.swap(t, pi);
            for row in &mut d {
                row.swap(t, pj);
            }
            for row in &mut v {
                row.swap(t, pj);
            }
            let pivot = d[t][t].clone();
            let mut clean = true;
            for i in t + 1..rows {
                if !d[i][t].is_zero() {
                    let q = floor_div(&d[i][t], &pivot);
                    row_sub(&mut d, i, t, &q);
                    row_sub(&mut u, i, t, &q);
                    if !d[i][t].is_zero() {
                        clean = false;
                    }
                }
            }
            for j in t + 1..cols {
                if !d[t][j].is_zero() {
                    let q = floor_div(&d[t][j], &pivot);
                    col_sub(&mut d, j, t, &q);
                    col_sub(&mut v, j, t, &q);
                    if !d[t][j].is_zero() {
                        clean = false;
                    }
                }
            }
            if !clean {
                continue;
            }
            // Divisibility of the rest by the pivot.
            let mut offender = None;
            'search: for (i, row) in d.iter().enumerate().skip(t + 1) {
                for entry in row.iter().skip(t + 1) {
                    if !(entry % &pivot).is_zero() {
                        offender = Some(i);
                        break 'search;
                    }
                }
            }
            if let Some(i) = offender {
                // Add row i to row t and reduce again.
                row_sub(&mut d, t, i, &-BigInt::one());
                row_sub(&mut u, t, i, &-BigInt::one());
                continue;
            }
            break;
        }
        if d[t][t].is_negative() {
            for value in &mut d[t] {
                *value = -value.clone();
            }
            for value in &mut u[t] {
                *value = -value.clone();
            }
        }
    }
    (u, d, v)
}

/// Row-style Hermite normal form `U A = H`.
fn hermite(a: &Ints) -> (Ints, Ints) {
    let rows = a.len();
    let cols = a.first().map_or(0, Vec::len);
    let mut h = a.clone();
    let mut u = integer_identity(rows);
    let mut lead = 0;
    let mut pivots = Vec::new();
    for c in 0..cols {
        if lead >= rows {
            break;
        }
        // Euclid on the rows from `lead` down in column c.
        loop {
            let mut nonzero: Vec<usize> = (lead..rows).filter(|&r| !h[r][c].is_zero()).collect();
            if nonzero.is_empty() {
                break;
            }
            nonzero.sort_by(|&x, &y| h[x][c].abs().cmp(&h[y][c].abs()));
            let smallest = nonzero[0];
            h.swap(lead, smallest);
            u.swap(lead, smallest);
            let pivot = h[lead][c].clone();
            let mut done = true;
            for r in lead + 1..rows {
                if !h[r][c].is_zero() {
                    let q = floor_div(&h[r][c], &pivot);
                    row_sub(&mut h, r, lead, &q);
                    row_sub(&mut u, r, lead, &q);
                    if !h[r][c].is_zero() {
                        done = false;
                    }
                }
            }
            if done {
                break;
            }
        }
        if h[lead][c].is_zero() {
            continue;
        }
        if h[lead][c].is_negative() {
            for value in &mut h[lead] {
                *value = -value.clone();
            }
            for value in &mut u[lead] {
                *value = -value.clone();
            }
        }
        pivots.push((lead, c));
        lead += 1;
    }
    // Reduce the entries above each pivot.
    for &(r, c) in &pivots {
        let pivot = h[r][c].clone();
        for above in 0..r {
            let q = floor_div(&h[above][c], &pivot);
            if !q.is_zero() {
                row_sub(&mut h, above, r, &q);
                row_sub(&mut u, above, r, &q);
            }
        }
    }
    (h, u)
}

// ----------------------------------------------------------------------
// Miscellaneous
// ----------------------------------------------------------------------

fn kron(
    cx: &mut Cx<'_>,
    a: &[Vec<NodeId>],
    b: &[Vec<NodeId>],
) -> Mat {
    let (ar, ac) = (a.len(), a.first().map_or(0, Vec::len));
    let (br, bc) = (b.len(), b.first().map_or(0, Vec::len));
    let mut out = Vec::with_capacity(ar * br);
    for a_row in a {
        for b_row in b {
            let mut row = Vec::with_capacity(ac * bc);
            for &x in a_row {
                for &y in b_row {
                    let p = mul(cx.graph, &[x, y]);
                    row.push(cx.simplify(p));
                }
            }
            out.push(row);
        }
    }
    out
}

fn minor(
    m: &[Vec<NodeId>],
    row: usize,
    col: usize,
) -> Mat {
    m.iter()
        .enumerate()
        .filter(|&(r, _)| r != row)
        .map(|(_, r)| {
            r.iter()
                .enumerate()
                .filter(|&(c, _)| c != col)
                .map(|(_, &e)| e)
                .collect()
        })
        .collect()
}

fn cofactors(
    cx: &mut Cx<'_>,
    m: &[Vec<NodeId>],
) -> Option<Mat> {
    let n = m.len();
    if !square(m) {
        return None;
    }
    if n == 1 {
        return Some(vec![vec![cx.graph.int(1)]]);
    }
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let mut row = Vec::with_capacity(n);
        for j in 0..n {
            let sub = minor(m, i, j);
            let value = determinant(cx, &sub)?;
            let value = if (i + j) % 2 == 0 {
                value
            } else {
                neg(cx.graph, value)
            };
            row.push(normal(cx, value));
        }
        out.push(row);
    }
    Some(out)
}

fn companion(
    cx: &mut Cx<'_>,
    p: NodeId,
    x: NodeId,
) -> Option<Mat> {
    let mut gens = crate::rules::poly::repr::Gens::default();
    let gx = gens.index(cx.graph, x);
    let poly = crate::rules::poly::repr::from_term(
        cx.graph,
        &mut gens,
        p,
        crate::rules::poly::repr::Limits::default(),
    )?;
    let mut degree = 0_usize;
    let mut current = poly;
    let mut coefficients = Vec::new();
    let zero = cx.graph.int(0);
    let mut factorial = BigInt::one();
    while !current.is_zero() {
        let term = crate::rules::poly::repr::to_term(cx.graph, &gens, &current);
        let at_zero = cx.graph.substitute(term, x, zero);
        let divisor = cx.graph.num(Number::rat(BigRational::new(
            BigInt::one(),
            factorial.clone(),
        )));
        let value = mul(cx.graph, &[divisor, at_zero]);
        coefficients.push(cx.simplify(value));
        degree += 1;
        factorial *= degree;
        current = current.derivative(gx);
    }
    let n = degree.checked_sub(1)?;
    if n == 0 {
        return None;
    }
    let lead = *coefficients.last()?;
    let mut out = zero_matrix(cx.graph, n, n);
    let one = cx.graph.int(1);
    for i in 1..n {
        out[i][i - 1] = one;
    }
    for (i, &c) in coefficients.iter().take(n).enumerate() {
        let inverse_lead = reciprocal(cx.graph, lead);
        let quotient = mul(cx.graph, &[c, inverse_lead]);
        let negated = neg(cx.graph, quotient);
        out[i][n - 1] = normal(cx, negated);
    }
    Some(out)
}

fn positive_definite(
    graph: &Graph,
    m: &[Vec<NodeId>],
) -> Option<bool> {
    if !square(m) {
        return None;
    }
    let q: Vec<Vec<BigRational>> = m
        .iter()
        .map(|r| {
            r.iter()
                .map(|&e| graph.number_of(e).and_then(Number::to_rational))
                .collect()
        })
        .collect::<Option<_>>()?;
    for (i, row) in q.iter().enumerate() {
        for (j, entry) in row.iter().enumerate().take(i) {
            if *entry != q[j][i] {
                return Some(false);
            }
        }
    }
    // Gaussian elimination without pivoting: all pivots positive.
    let n = q.len();
    let mut work = q;
    for col in 0..n {
        if !work[col][col].is_positive() {
            return Some(false);
        }
        let pivot_row = work[col].clone();
        for row in work.iter_mut().skip(col + 1) {
            let factor = &row[col] / &pivot_row[col];
            for (value, p) in row.iter_mut().zip(&pivot_row).skip(col) {
                *value -= &factor * p;
            }
        }
    }
    Some(true)
}

// ----------------------------------------------------------------------
// Kernel
// ----------------------------------------------------------------------

#[derive(Copy, Clone)]
enum Kind {
    Kron,
    Diag,
    Vandermonde,
    Hilbert,
    Companion,
    Adjugate,
    Cofactors,
    Colspace,
    Rowspace,
    LeftNullspace,
    Pinv,
    Lstsq,
    Minpoly,
    Matpoly,
    CayleyHamilton,
    Jordan,
    JordanBlocks,
    Matfun,
    Matexp,
    Mpow,
    Matsqrt,
    Spectral,
    Ldl,
    Cholesky,
    Orthogonalize,
    Orthonormalize,
    Smith,
    InvariantFactors,
    Hermite,
    IsSymmetric,
    IsOrthogonal,
    IsPositiveDefinite,
}

struct Extension {
    op: OpId,
    kind: Kind,
}

fn vectors_of(
    graph: &mut Graph,
    node: NodeId,
) -> Option<Vec<Vec<NodeId>>> {
    let items = vector(graph, node)?;
    items.into_iter().map(|v| vector(graph, v)).collect()
}

impl Extension {
    #[allow(clippy::too_many_lines)]
    fn compute(
        &self,
        cx: &mut Cx<'_>,
        args: &[NodeId],
    ) -> Option<NodeId> {
        let arg = |i: usize| args.get(i).copied();
        match self.kind {
            | Kind::Kron => {
                let (a, b) = (matrix(cx.graph, arg(0)?)?, matrix(cx.graph, arg(1)?)?);
                let k = kron(cx, &a, &b);
                Some(matrix_term(cx.graph, &k))
            },
            | Kind::Diag => {
                let d = vector(cx.graph, arg(0)?)?;
                let mut m = zero_matrix(cx.graph, d.len(), d.len());
                for (i, &v) in d.iter().enumerate() {
                    m[i][i] = v;
                }
                Some(matrix_term(cx.graph, &m))
            },
            | Kind::Vandermonde => {
                let xs = vector(cx.graph, arg(0)?)?;
                let mut rows = Vec::with_capacity(xs.len());
                for &x in &xs {
                    let mut row = Vec::with_capacity(xs.len());
                    for k in 0..xs.len() {
                        let e = cx.graph.int(i64::try_from(k).ok()?);
                        let power = cx.graph.node(core::POW, &[x, e]);
                        row.push(cx.simplify(power));
                    }
                    rows.push(row);
                }
                Some(matrix_term(cx.graph, &rows))
            },
            | Kind::Hilbert => {
                let n = usize::try_from(cx.graph.number_of(arg(0)?)?.to_i64()?)
                    .ok()
                    .filter(|&n| (1..=64).contains(&n))?;
                let mut rows = Vec::with_capacity(n);
                for i in 0..n {
                    let mut row = Vec::with_capacity(n);
                    for j in 0..n {
                        let denominator = i64::try_from(i + j + 1).ok()?;
                        row.push(cx.graph.num(Number::fraction(1, denominator)?));
                    }
                    rows.push(row);
                }
                Some(matrix_term(cx.graph, &rows))
            },
            | Kind::Companion => {
                let p = best(cx.graph, arg(0)?)?;
                let c = companion(cx, p, arg(1)?)?;
                Some(matrix_term(cx.graph, &c))
            },
            | Kind::Cofactors | Kind::Adjugate => {
                let m = matrix(cx.graph, arg(0)?)?;
                let c = cofactors(cx, &m)?;
                let c = if matches!(self.kind, Kind::Adjugate) {
                    transpose(&c)
                } else {
                    c
                };
                Some(matrix_term(cx.graph, &c))
            },
            | Kind::Colspace => {
                let m = matrix(cx.graph, arg(0)?)?;
                let basis = colspace(cx, &m)?;
                let items: Vec<NodeId> = basis.iter().map(|v| list(cx.graph, v)).collect();
                Some(list(cx.graph, &items))
            },
            | Kind::Rowspace => {
                let m = matrix(cx.graph, arg(0)?)?;
                let basis = rowspace(cx, &m)?;
                let items: Vec<NodeId> = basis.iter().map(|v| list(cx.graph, v)).collect();
                Some(list(cx.graph, &items))
            },
            | Kind::LeftNullspace => {
                let m = matrix(cx.graph, arg(0)?)?;
                let basis = nullspace(cx, &transpose(&m))?;
                let items: Vec<NodeId> = basis.iter().map(|v| list(cx.graph, v)).collect();
                Some(list(cx.graph, &items))
            },
            | Kind::Pinv => {
                let m = matrix(cx.graph, arg(0)?)?;
                let p = pseudo_inverse(cx, &m)?;
                Some(matrix_term(cx.graph, &p))
            },
            | Kind::Lstsq => {
                let m = matrix(cx.graph, arg(0)?)?;
                let b = vector(cx.graph, arg(1)?)?;
                if b.len() != m.len() {
                    return None;
                }
                let p = pseudo_inverse(cx, &m)?;
                let x = apply(cx, &p, &b);
                Some(list(cx.graph, &x))
            },
            | Kind::Minpoly => {
                let m = matrix(cx.graph, arg(0)?)?;
                minimal_polynomial(cx, &m, arg(1)?)
            },
            | Kind::Matpoly => {
                let m = matrix(cx.graph, arg(2)?)?;
                let p = best(cx.graph, arg(0)?)?;
                let result = matrix_polynomial(cx, p, arg(1)?, &m)?;
                Some(matrix_term(cx.graph, &result))
            },
            | Kind::CayleyHamilton => {
                let m = matrix(cx.graph, arg(0)?)?;
                let l = fresh(cx);
                let p = charpoly(cx, &m, l)?;
                let result = matrix_polynomial(cx, p, l, &m)?;
                Some(matrix_term(cx.graph, &result))
            },
            | Kind::Jordan => {
                let m = matrix(cx.graph, arg(0)?)?;
                jordan(cx, &m)
            },
            | Kind::JordanBlocks => {
                let m = matrix(cx.graph, arg(0)?)?;
                jordan_block_list(cx, &m)
            },
            | Kind::Matfun => {
                let m = matrix(cx.graph, arg(2)?)?;
                let f = best(cx.graph, arg(0)?)?;
                let result = matrix_function(cx, f, arg(1)?, &m)?;
                Some(matrix_term(cx.graph, &result))
            },
            | Kind::Matexp => {
                let m = matrix(cx.graph, arg(0)?)?;
                let result = matrix_exponential(cx, &m, arg(1)?)?;
                Some(matrix_term(cx.graph, &result))
            },
            | Kind::Mpow => {
                let m = matrix(cx.graph, arg(0)?)?;
                let result = matrix_power(cx, &m, arg(1)?)?;
                Some(matrix_term(cx.graph, &result))
            },
            | Kind::Matsqrt => {
                let m = matrix(cx.graph, arg(0)?)?;
                let x = fresh(cx);
                let f = sqrt(cx.graph, x)?;
                let result = matrix_function(cx, f, x, &m)?;
                Some(matrix_term(cx.graph, &result))
            },
            | Kind::Spectral => {
                let m = matrix(cx.graph, arg(0)?)?;
                spectral(cx, &m)
            },
            | Kind::Ldl => {
                let m = matrix(cx.graph, arg(0)?)?;
                if !is_symmetric(cx, &m) {
                    return None;
                }
                let (l, d) = ldl(cx, &m)?;
                let (l, d) = (matrix_term(cx.graph, &l), list(cx.graph, &d));
                Some(list(cx.graph, &[l, d]))
            },
            | Kind::Cholesky => {
                let m = matrix(cx.graph, arg(0)?)?;
                let l = cholesky(cx, &m)?;
                Some(matrix_term(cx.graph, &l))
            },
            | Kind::Orthogonalize | Kind::Orthonormalize => {
                let vs = vectors_of(cx.graph, arg(0)?)?;
                let basis = gram_schmidt(cx, &vs, matches!(self.kind, Kind::Orthonormalize))?;
                let items: Vec<NodeId> = basis.iter().map(|v| list(cx.graph, v)).collect();
                Some(list(cx.graph, &items))
            },
            | Kind::Smith => {
                let a = integer_matrix(cx.graph, arg(0)?)?;
                let (u, d, v) = smith(&a);
                let (u, d, v) = (
                    integers_to_term(cx.graph, &u),
                    integers_to_term(cx.graph, &d),
                    integers_to_term(cx.graph, &v),
                );
                Some(list(cx.graph, &[u, d, v]))
            },
            | Kind::InvariantFactors => {
                let a = integer_matrix(cx.graph, arg(0)?)?;
                let (_, d, _) = smith(&a);
                let cols = d.first().map_or(0, Vec::len);
                let items: Vec<NodeId> = (0..d.len().min(cols))
                    .filter(|&i| !d[i][i].is_zero())
                    .map(|i| cx.graph.num(Number::Int(d[i][i].clone())))
                    .collect();
                Some(list(cx.graph, &items))
            },
            | Kind::Hermite => {
                let a = integer_matrix(cx.graph, arg(0)?)?;
                let (h, u) = hermite(&a);
                let (h, u) = (
                    integers_to_term(cx.graph, &h),
                    integers_to_term(cx.graph, &u),
                );
                Some(list(cx.graph, &[h, u]))
            },
            | Kind::IsSymmetric => {
                let m = matrix(cx.graph, arg(0)?)?;
                let value = is_symmetric(cx, &m);
                Some(cx.graph.lit(crate::graph::Payload::Bool(value)))
            },
            | Kind::IsOrthogonal => {
                let m = matrix(cx.graph, arg(0)?)?;
                if !square(&m) {
                    return Some(cx.graph.lit(crate::graph::Payload::Bool(false)));
                }
                let gram = product(cx, &transpose(&m), &m)?;
                let identity = identity_rows(cx.graph, m.len());
                let difference = sum(cx, &gram, &identity, -1)?;
                let mut value = true;
                for entry in difference.iter().flatten() {
                    if !cx.is_zero(*entry) {
                        value = false;
                    }
                }
                Some(cx.graph.lit(crate::graph::Payload::Bool(value)))
            },
            | Kind::IsPositiveDefinite => {
                let m = matrix(cx.graph, arg(0)?)?;
                let value = positive_definite(cx.graph, &m)?;
                Some(cx.graph.lit(crate::graph::Payload::Bool(value)))
            },
        }
    }
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
        self.compute(cx, &args)
            .map_or(Outcome::Pass, Outcome::Equal)
    }

    fn revisit(&self) -> bool {
        true
    }
}

/// The dimension signature of a matrix extension.
fn extension_rule(kind: Kind) -> DimRule {
    match kind {
        | Kind::Kron => DimRule::Product,
        | Kind::Diag
        | Kind::Hilbert
        | Kind::Adjugate
        | Kind::Cofactors
        | Kind::Colspace
        | Kind::Rowspace
        | Kind::LeftNullspace
        | Kind::Jordan
        | Kind::JordanBlocks
        | Kind::Matfun
        | Kind::Matexp
        | Kind::Mpow
        | Kind::Matsqrt
        | Kind::Spectral
        | Kind::Ldl
        | Kind::Cholesky
        | Kind::Orthogonalize
        | Kind::Orthonormalize
        | Kind::Smith
        | Kind::InvariantFactors
        | Kind::Hermite => DimRule::List,
        | Kind::Vandermonde
        | Kind::Companion
        | Kind::Matpoly
        | Kind::Minpoly
        | Kind::CayleyHamilton
        | Kind::IsSymmetric
        | Kind::IsOrthogonal
        | Kind::IsPositiveDefinite => DimRule::NonMetric,
        | Kind::Pinv => DimRule::Inverse,
        | Kind::Lstsq => DimRule::Custom(quotient_dimension),
    }
}

fn install_ext(i: &mut Installer<'_>) -> Result<(), RuleError> {
    for (name, arity, kind) in [
        ("kron", 2, Kind::Kron),
        ("diag", 1, Kind::Diag),
        ("vandermonde", 1, Kind::Vandermonde),
        ("hilbert", 1, Kind::Hilbert),
        ("companion", 2, Kind::Companion),
        ("adjugate", 1, Kind::Adjugate),
        ("cofactors", 1, Kind::Cofactors),
        ("colspace", 1, Kind::Colspace),
        ("rowspace", 1, Kind::Rowspace),
        ("left_nullspace", 1, Kind::LeftNullspace),
        ("pinv", 1, Kind::Pinv),
        ("lstsq", 2, Kind::Lstsq),
        ("minpoly", 2, Kind::Minpoly),
        ("matpoly", 3, Kind::Matpoly),
        ("cayley_hamilton", 1, Kind::CayleyHamilton),
        ("jordan", 1, Kind::Jordan),
        ("jordan_blocks", 1, Kind::JordanBlocks),
        ("matfun", 3, Kind::Matfun),
        ("matexp", 2, Kind::Matexp),
        ("mpow", 2, Kind::Mpow),
        ("matsqrt", 1, Kind::Matsqrt),
        ("spectral", 1, Kind::Spectral),
        ("ldl", 1, Kind::Ldl),
        ("cholesky", 1, Kind::Cholesky),
        ("orthogonalize", 1, Kind::Orthogonalize),
        ("orthonormalize", 1, Kind::Orthonormalize),
        ("smith", 1, Kind::Smith),
        ("invariant_factors", 1, Kind::InvariantFactors),
        ("hermite_form", 1, Kind::Hermite),
        ("is_symmetric", 1, Kind::IsSymmetric),
        ("is_orthogonal", 1, Kind::IsOrthogonal),
        ("is_positive_definite", 1, Kind::IsPositiveDefinite),
    ] {
        let op = i.op(OpDescriptor::new(name, Arity::Fixed(arity))
            .flags(OpFlags::HEAVY)
            .cost(100))?;
        i.dim(op, extension_rule(kind));
        i.kernel(
            &format!("linalg/{name}"),
            Tier::Reduce,
            Extension { op, kind },
        );
    }
    Ok(())
}

#[cfg(test)]
mod ext_tests {
    use crate::rules::linalg::linalg;
    use crate::rules::testing::simplify;

    fn run(src: &str) -> String {
        simplify(&[linalg()], src)
    }

    const A: &str = "list(list(2, 1), list(0, 2))";
    const B: &str = "list(list(4, 1, 0), list(0, 4, 1), list(0, 0, 4))";

    #[test]
    fn products_and_special_matrices() {
        assert_eq!(
            run("kron(list(list(1, 2), list(3, 4)), list(list(0, 1), list(1, 0)))"),
            "list(list(0, 1, 0, 2), list(1, 0, 2, 0), list(0, 3, 0, 4), list(3, 0, 4, 0))"
        );
        assert_eq!(
            run("diag(list(1, a, 3))"),
            "list(list(1, 0, 0), list(0, a, 0), list(0, 0, 3))"
        );
        assert_eq!(
            run("vandermonde(list(1, 2, 3))"),
            "list(list(1, 1, 1), list(1, 2, 4), list(1, 3, 9))"
        );
        assert_eq!(run("hilbert(2)"), "list(list(1, 1/2), list(1/2, 1/3))");
        assert_eq!(
            run("companion(x^3 - 2*x^2 + 3*x - 4, x)"),
            "list(list(0, 0, 4), list(1, 0, -3), list(0, 1, 2))"
        );
        assert_eq!(run("det(vandermonde(list(a, b)))"), "b - a");
    }

    #[test]
    fn adjugate_and_spaces() {
        assert_eq!(
            run("adjugate(list(list(1, 2), list(3, 4)))"),
            "list(list(4, -2), list(-3, 1))"
        );
        assert_eq!(
            run("cofactors(list(list(1, 2), list(3, 4)))"),
            "list(list(4, -3), list(-2, 1))"
        );
        assert_eq!(
            run("colspace(list(list(1, 2, 3), list(2, 4, 6), list(1, 0, 1)))"),
            "list(list(1, 2, 1), list(2, 4, 0))"
        );
        assert_eq!(
            run("rowspace(list(list(1, 2, 3), list(4, 5, 6)))"),
            "list(list(1, 0, -1), list(0, 1, 2))"
        );
        assert_eq!(
            run("left_nullspace(list(list(1, 2), list(2, 4)))"),
            "list(list(-2, 1))"
        );
    }

    #[test]
    fn pseudo_inverse_and_least_squares() {
        assert_eq!(
            run("pinv(list(list(1, 0), list(0, 2), list(0, 0)))"),
            "list(list(1, 0, 0), list(0, 1/2, 0))"
        );
        // Rank deficient.
        assert_eq!(
            run("pinv(list(list(1, 1), list(1, 1)))"),
            "list(list(1/4, 1/4), list(1/4, 1/4))"
        );
        // The line fit y = a + b x through (0,1), (1,2), (2,2): a = 7/6, b = 1/2.
        assert_eq!(
            run("lstsq(list(list(1, 0), list(1, 1), list(1, 2)), list(1, 2, 2))"),
            "list(7/6, 1/2)"
        );
    }

    #[test]
    fn minimal_and_matrix_polynomials() {
        assert_eq!(run(&format!("minpoly({A}, l)")), "l^2 - 4*l + 4");
        assert_eq!(run("minpoly(list(list(2, 0), list(0, 2)), l)"), "l - 2");
        assert_eq!(
            run(&format!("cayley_hamilton({B})")),
            "list(list(0, 0, 0), list(0, 0, 0), list(0, 0, 0))"
        );
        assert_eq!(
            run("cayley_hamilton(list(list(a, b), list(c, d)))"),
            "list(list(0, 0), list(0, 0))"
        );
        assert_eq!(
            run(&format!("matpoly(x^2 - 4*x + 4, x, {A})")),
            "list(list(0, 0), list(0, 0))"
        );
        assert_eq!(
            run("matpoly(x^2 + 1, x, list(list(0, -1), list(1, 0)))"),
            "list(list(0, 0), list(0, 0))"
        );
    }

    #[test]
    fn jordan_form_and_matrix_functions() {
        assert_eq!(run(&format!("jordan_blocks({B})")), "list(list(4, 3))");
        assert_eq!(
            run("jordan_blocks(list(list(2, 0), list(0, 2)))"),
            "list(list(2, 1), list(2, 1))"
        );
        // A P = P J.
        let pj = run(&format!("jordan({B})"));
        assert!(pj.starts_with("list(list("), "{pj}");
        let p = tests_part(&pj, 0);
        let j = tests_part(&pj, 1);
        assert_eq!(
            run(&format!(
                "madd(matmul({B}, {p}), smul(-1, matmul({p}, {j})))"
            )),
            "list(list(0, 0, 0), list(0, 0, 0), list(0, 0, 0))"
        );
        // exp(A t) for a Jordan block: e^(2t) (I + t N).
        assert_eq!(
            run(&format!("matexp({A}, t)")),
            "list(list(exp(2*t), t*exp(2*t)), list(0, exp(2*t)))"
        );
        assert_eq!(
            run("matexp(list(list(1, 0), list(0, 2)), t)"),
            "list(list(exp(t), 0), list(0, exp(2*t)))"
        );
        assert_eq!(
            run(&format!("mpow({A}, 5)")),
            "list(list(32, 80), list(0, 32))"
        );
        assert_eq!(
            run(&format!("mpow({A}, -1)")),
            "list(list(1/2, -1/4), list(0, 1/2))"
        );
        assert_eq!(
            run(&format!("mpow({A}, n)")),
            "list(list(2^n, n*2^(n - 1)), list(0, 2^n))"
        );
        assert_eq!(
            run("matsqrt(list(list(4, 0), list(0, 9)))"),
            "list(list(2, 0), list(0, 3))"
        );
        assert_eq!(
            run(&format!("matfun(x^2 + 1, x, {A})")),
            "list(list(5, 4), list(0, 5))"
        );
    }

    #[test]
    fn spectral_decomposition() {
        assert_eq!(
            run("spectral(list(list(2, 1), list(1, 2)))"),
            "list(list(1, list(list(1/2, -1/2), list(-1/2, 1/2))), list(3, list(list(1/2, 1/2), list(1/2, 1/2))))"
        );
    }

    #[test]
    fn symmetric_factorisations() {
        assert_eq!(
            run("ldl(list(list(4, 2), list(2, 3)))"),
            "list(list(list(1, 0), list(1/2, 1)), list(4, 2))"
        );
        assert_eq!(
            run("cholesky(list(list(4, 2), list(2, 3)))"),
            "list(list(2, 0), list(1, 2^(1/2)))"
        );
        let (text, reduced) = crate::rules::testing::reduce_with(
            &[linalg()],
            "cholesky(list(list(1, 2), list(2, 1)))",
            &[],
        );
        assert!(!reduced, "not positive definite: {text}");
        assert_eq!(
            run("orthogonalize(list(list(1, 1, 0), list(1, 0, 1), list(0, 1, 1)))"),
            "list(list(1, 1, 0), list(1/2, -1/2, 1), list(-2/3, 2/3, 2/3))"
        );
        assert_eq!(
            run("orthonormalize(list(list(3, 4), list(1, 0)))"),
            "list(list(3/5, 4/5), list(4/5, -3/5))"
        );
        assert_eq!(run("is_symmetric(list(list(1, 2), list(2, 5)))"), "true");
        assert_eq!(run("is_orthogonal(list(list(0, 1), list(-1, 0)))"), "true");
        assert_eq!(
            run("is_positive_definite(list(list(2, 1), list(1, 2)))"),
            "true"
        );
        assert_eq!(
            run("is_positive_definite(list(list(1, 2), list(2, 1)))"),
            "false"
        );
    }

    #[test]
    fn integer_normal_forms() {
        let a = "list(list(2, 4, 4), list(-6, 6, 12), list(10, 4, 16))";
        assert_eq!(run(&format!("invariant_factors({a})")), "list(2, 2, 156)");
        let smith = run(&format!("smith({a})"));
        let (u, d, v) = (
            tests_part(&smith, 0),
            tests_part(&smith, 1),
            tests_part(&smith, 2),
        );
        assert_eq!(run(&format!("matmul({u}, {a}, {v})")), d);
        assert_eq!(run("det(list(list(1, 0), list(0, 1)))"), "1");
        let hermite =
            run("hermite_form(list(list(2, 3, 6, 2), list(5, 6, 1, 6), list(8, 3, 1, 1)))");
        let (h, u) = (tests_part(&hermite, 0), tests_part(&hermite, 1));
        assert_eq!(
            run(&format!(
                "matmul({u}, list(list(2, 3, 6, 2), list(5, 6, 1, 6), list(8, 3, 1, 1)))"
            )),
            h
        );
        assert_eq!(
            h,
            "list(list(1, 0, 50, -11), list(0, 3, 28, -2), list(0, 0, 61, -13))"
        );
        assert_eq!(run("det(list(list(1, 2), list(3, 4)))"), "-2");
    }

    /// The `index`-th top-level item of a `list(...)` text.
    fn tests_part(
        text: &str,
        index: usize,
    ) -> String {
        let inner = text
            .strip_prefix("list(")
            .and_then(|s| s.strip_suffix(')'))
            .unwrap_or("");
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
}
