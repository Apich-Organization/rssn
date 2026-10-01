//! Lie algebras: brackets, structure constants, adjoint and Killing forms,
//! the exponential map and the Baker–Campbell–Hausdorff series.
//!
//! An element of a *matrix* Lie algebra is a square matrix
//! `list(list(..), ..)`; an element of a Lie algebra of *vector fields* is
//! the list of its component expressions, and every operator that takes a
//! basis then takes the list of coordinate symbols as a last argument.
//! The bracket is the commutator `[A, B] = A B - B A` for matrices and
//! `[X, Y]^i = X^j ∂_j Y^i - Y^j ∂_j X^i` for fields.
//!
//! | operator | value |
//! |---|---|
//! | `lie_bracket(A, B)`, `lie_bracket(X, Y, coords)` | the bracket |
//! | `structure_constants(basis[, coords])` | `c[i][j][k]` with `[e_i, e_j] = Σ_k c_ijk e_k` |
//! | `adjoint_rep(X, basis[, coords])` | matrix of `ad_X` in the basis: column `j` holds the coordinates of `[X, e_j]` |
//! | `adjoint_group(g, X)` | `Ad_g X = g X g⁻¹` |
//! | `killing_form(basis[, coords])` | `K_ab = tr(ad_a ad_b)` |
//! | `check_jacobi(basis[, coords])` | `true` when `[x,[y,z]] + [y,[z,x]] + [z,[x,y]] = 0` for all basis triples |
//! | `is_lie_algebra(basis[, coords])` | closure of the span under the bracket, and Jacobi |
//! | `commutator_table(basis[, coords])` | `list(list([e_i, e_j]))` |
//! | `exp_map(A)` | the matrix exponential in closed form (see below) |
//! | `exp_map(A, n)` | the Taylor polynomial `Σ_{k≤n} A^k / k!` |
//! | `baker_campbell_hausdorff(X, Y, n[, coords])` | `log(e^X e^Y)` up to order `n ≤ 4` |
//! | `so3_basis()`, `su2_basis()`, `sl2_basis()` | standard bases |
//!
//! # Closed forms of `exp_map`
//!
//! In order: a nilpotent matrix (finite series); a diagonal matrix; a 2×2
//! matrix (`e^A = e^m (cosh s + sinh(s)/s · B)` with `B = A - m I`, `m` the
//! mean eigenvalue and `s² = B²`, which becomes `cos`/`sin` when `s²` is
//! negative); a 3×3 antisymmetric matrix (Rodrigues' formula); a
//! diagonalizable matrix whose eigenvalues are known (`P e^D P⁻¹`).
//! Anything else is left unreduced.
//!
//! # Conventions of the standard bases
//!
//! * `so3_basis() = (Lx, Ly, Lz)`, the infinitesimal rotations about the
//!   axes, with `[Li, Lj] = ε_ijk Lk`.
//! * `su2_basis() = (σx, σy, σz) / (2i)`, i.e. `-iσ/2`, so that the
//!   structure constants are again `ε_ijk`. (The legacy implementation used
//!   `+iσ/2`, which flips the sign of every structure constant.)
//! * `sl2_basis() = (h, e, f)` with `[h, e] = 2e`, `[h, f] = -2f`,
//!   `[e, f] = h`.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::Signed;

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
use crate::rules::poly::best;

use super::calculus::derivative;
use super::complex::complex;
use super::elementary::elementary;
use super::linalg::linalg;
use super::logic::logic;

/// The Lie algebra rule set.
#[must_use]
pub fn lie() -> RuleSet {
    RuleSet::new("lie", install).needs(linalg()).needs(complex()).needs(elementary()).needs(logic())
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Request {
    Bracket,
    StructureConstants,
    AdjointRep,
    AdjointGroup,
    Killing,
    Jacobi,
    IsLieAlgebra,
    CommutatorTable,
    ExpMap,
    Bch,
    So3,
    Su2,
    Sl2,
}

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    for (name, arity, request) in [
        ("lie_bracket", Arity::Variadic, Request::Bracket),
        ("structure_constants", Arity::Variadic, Request::StructureConstants),
        ("adjoint_rep", Arity::Variadic, Request::AdjointRep),
        ("adjoint_group", Arity::Fixed(2), Request::AdjointGroup),
        ("killing_form", Arity::Variadic, Request::Killing),
        ("check_jacobi", Arity::Variadic, Request::Jacobi),
        ("is_lie_algebra", Arity::Variadic, Request::IsLieAlgebra),
        ("commutator_table", Arity::Variadic, Request::CommutatorTable),
        ("exp_map", Arity::Variadic, Request::ExpMap),
        ("baker_campbell_hausdorff", Arity::Variadic, Request::Bch),
        ("so3_basis", Arity::Fixed(0), Request::So3),
        ("su2_basis", Arity::Fixed(0), Request::Su2),
        ("sl2_basis", Arity::Fixed(0), Request::Sl2),
    ] {
        let op = i.op(OpDescriptor::new(name, arity).flags(OpFlags::HEAVY).cost(100))?;
        i.kernel(&format!("lie/{name}"), Tier::Reduce, LieKernel { op, request });
    }
    Ok(())
}

// ----------------------------------------------------------------------
// Term helpers (shared with the geometric algebra rule set)
// ----------------------------------------------------------------------

/// `name(args...)` for an installed operator.
pub(crate) fn call(
    graph: &mut Graph,
    name: &str,
    args: &[NodeId],
) -> Option<NodeId> {
    let op = graph.ops().lookup(name)?;
    graph.try_node(op, args)
}

/// Product of `factors` (no node for fewer than two).
pub(crate) fn mul(
    graph: &mut Graph,
    factors: &[NodeId],
) -> NodeId {
    match factors {
        | [] => graph.int(1),
        | [only] => *only,
        | _ => graph.node(core::MUL, factors),
    }
}

/// Sum of `terms` (no node for fewer than two).
pub(crate) fn add(
    graph: &mut Graph,
    terms: &[NodeId],
) -> NodeId {
    match terms {
        | [] => graph.int(0),
        | [only] => *only,
        | _ => graph.node(core::ADD, terms),
    }
}

/// `x^e` for an integer exponent.
pub(crate) fn pow_int(
    graph: &mut Graph,
    x: NodeId,
    e: i64,
) -> NodeId {
    let e = graph.int(e);
    graph.node(core::POW, &[x, e])
}

/// Whether `node` is the literal zero.
pub(crate) fn is_zero_literal(
    graph: &Graph,
    node: NodeId,
) -> bool {
    graph.number_of(node).is_some_and(Number::is_zero)
}

fn rational_sqrt(r: &BigRational) -> Option<Number> {
    if r.is_negative() {
        return None;
    }
    let n = r.numer().sqrt();
    let d = r.denom().sqrt();
    (&n * &n == *r.numer() && &d * &d == *r.denom()).then(|| Number::rat(BigRational::new(n, d)))
}

/// A `w` with `w² = t` when one is evident from the form of `t`: a
/// perfect-square rational, an even power, a product of such. The sign of
/// `w` is not specified.
pub(crate) fn exact_sqrt(
    graph: &mut Graph,
    t: NodeId,
) -> Option<NodeId> {
    if let Some(n) = graph.number_of(t) {
        let root = rational_sqrt(&n.to_rational()?)?;
        return Some(graph.num(root));
    }
    let kids = graph.children(t).to_vec();
    match graph.op(t) {
        | core::POW => {
            let &[base, exponent] = kids.as_slice() else {
                return None;
            };
            let e = graph.number_of(exponent)?.to_rational()?;
            let two = BigRational::from_integer(BigInt::from(2));
            let half = &e / &two;
            if !e.is_integer() || !half.is_integer() || !e.is_positive() {
                return None;
            }
            if half == BigRational::from_integer(BigInt::from(1)) {
                return Some(base);
            }
            let h = graph.num(Number::rat(half));
            Some(graph.node(core::POW, &[base, h]))
        },
        | core::MUL => {
            let roots: Vec<NodeId> = kids.iter().map(|&k| exact_sqrt(graph, k)).collect::<Option<_>>()?;
            Some(mul(graph, &roots))
        },
        | _ => None,
    }
}

/// `(c, s)` with `exp(X) = c + s X` for every `X` with `X² = delta`
/// (a scalar): `cos`/`sin` when `delta < 0`, `cosh`/`sinh` when `delta > 0`.
pub(crate) fn exp_even(
    cx: &mut Cx<'_>,
    delta: NodeId,
) -> Option<(NodeId, NodeId)> {
    let delta = cx.simplify(delta);
    if is_zero_literal(cx.graph, delta) {
        return Some((cx.graph.int(1), cx.graph.int(1)));
    }
    let minus_one = cx.graph.int(-1);
    let negated = mul(cx.graph, &[minus_one, delta]);
    let negated = cx.simplify(negated);
    let half = cx.graph.num(Number::fraction(1, 2)?);
    let (trig, w) = if let Some(w) = exact_sqrt(cx.graph, negated) {
        (true, w)
    } else if let Some(w) = exact_sqrt(cx.graph, delta) {
        (false, w)
    } else if cx.graph.number_of(delta).is_some_and(Number::is_negative) {
        (true, cx.graph.node(core::POW, &[negated, half]))
    } else {
        (false, cx.graph.node(core::POW, &[delta, half]))
    };
    let (cos, sin) = if trig { ("cos", "sin") } else { ("cosh", "sinh") };
    let c = call(cx.graph, cos, &[w])?;
    let sn = call(cx.graph, sin, &[w])?;
    let w_inv = pow_int(cx.graph, w, -1);
    let s = mul(cx.graph, &[sn, w_inv]);
    Some((cx.simplify(c), cx.simplify(s)))
}

// ----------------------------------------------------------------------
// Elements
// ----------------------------------------------------------------------

/// An element of a Lie algebra, flattened: an `n × n` matrix (row-major)
/// or the components of a vector field.
#[derive(Clone, Debug)]
struct Elem {
    size: Option<usize>,
    data: Vec<NodeId>,
}

type Mat = Vec<Vec<NodeId>>;

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
    let nodes: Vec<NodeId> = rows.iter().map(|r| list(graph, r)).collect();
    list(graph, &nodes)
}

fn elem_term(
    graph: &mut Graph,
    e: &Elem,
) -> NodeId {
    match e.size {
        | Some(n) => {
            let rows: Mat = e.data.chunks(n).map(<[NodeId]>::to_vec).collect();
            matrix_term(graph, &rows)
        },
        | None => list(graph, &e.data),
    }
}

/// Entries of a `list(..)` after taking the best form of the term.
fn items(
    graph: &mut Graph,
    node: NodeId,
) -> Option<Vec<NodeId>> {
    let node = best(graph, node)?;
    (graph.op(node) == core::LIST).then(|| graph.children(node).to_vec())
}

fn matrix(
    graph: &mut Graph,
    node: NodeId,
) -> Option<Mat> {
    let rows = items(graph, node)?;
    let mut out = Vec::with_capacity(rows.len());
    for row in rows {
        out.push(items(graph, row)?);
    }
    let width = out.first().map_or(0, Vec::len);
    (width > 0 && out.iter().all(|r| r.len() == width)).then_some(out)
}

fn read_elem(
    graph: &mut Graph,
    node: NodeId,
    coords: Option<&[NodeId]>,
) -> Option<Elem> {
    if let Some(c) = coords {
        let data = items(graph, node)?;
        let scalar = data.iter().all(|&d| graph.op(d) != core::LIST);
        return (scalar && data.len() == c.len()).then_some(Elem { size: None, data });
    }
    let m = matrix(graph, node)?;
    (m.len() == m[0].len()).then(|| Elem { size: Some(m.len()), data: m.into_iter().flatten().collect() })
}

fn read_basis(
    graph: &mut Graph,
    node: NodeId,
    coords: Option<&[NodeId]>,
) -> Option<Vec<Elem>> {
    let nodes = items(graph, node)?;
    let basis: Vec<Elem> = nodes.into_iter().map(|n| read_elem(graph, n, coords)).collect::<Option<_>>()?;
    let first = basis.first()?;
    basis.iter().all(|e| e.size == first.size && e.data.len() == first.data.len()).then_some(basis)
}

fn read_coords(
    graph: &mut Graph,
    node: NodeId,
) -> Option<Vec<NodeId>> {
    let c = items(graph, node)?;
    c.iter().all(|&x| graph.symbol_of(x).is_some()).then_some(c)
}

/// Splits `args` into the first `fixed` arguments and the optional
/// coordinate list that may follow them.
fn with_coords(
    graph: &mut Graph,
    args: &[NodeId],
    fixed: usize,
) -> Option<(Vec<NodeId>, Option<Vec<NodeId>>)> {
    if args.len() == fixed {
        Some((args.to_vec(), None))
    } else if args.len() == fixed + 1 {
        Some((args[..fixed].to_vec(), Some(read_coords(graph, args[fixed])?)))
    } else {
        None
    }
}

/// Whether `t` is provably zero: it simplifies, possibly after expansion,
/// to the literal zero.
fn vanishes(
    cx: &mut Cx<'_>,
    t: NodeId,
) -> bool {
    let s = cx.simplify(t);
    if is_zero_literal(cx.graph, s) {
        return true;
    }
    if cx.graph.number_of(s).is_some() {
        return false;
    }
    call(cx.graph, "expand", &[s]).is_some_and(|e| {
        let e = cx.simplify(e);
        is_zero_literal(cx.graph, e)
    })
}

fn matprod(
    cx: &mut Cx<'_>,
    a: &[NodeId],
    b: &[NodeId],
    n: usize,
) -> Vec<NodeId> {
    let mut out = Vec::with_capacity(n * n);
    for i in 0..n {
        for j in 0..n {
            let terms: Vec<NodeId> = (0..n).map(|k| mul(cx.graph, &[a[i * n + k], b[k * n + j]])).collect();
            let sum = add(cx.graph, &terms);
            out.push(cx.simplify(sum));
        }
    }
    out
}

/// `sum coef_k * e_k`, entry by entry.
fn combine(
    cx: &mut Cx<'_>,
    parts: &[(NodeId, &Elem)],
) -> Option<Elem> {
    let first = parts.first()?.1;
    let mut data = Vec::with_capacity(first.data.len());
    for at in 0..first.data.len() {
        let terms: Vec<NodeId> = parts.iter().map(|&(c, e)| mul(cx.graph, &[c, e.data[at]])).collect();
        let sum = add(cx.graph, &terms);
        data.push(cx.simplify(sum));
    }
    Some(Elem { size: first.size, data })
}

fn bracket(
    cx: &mut Cx<'_>,
    a: &Elem,
    b: &Elem,
    coords: Option<&[NodeId]>,
) -> Option<Elem> {
    if a.size != b.size || a.data.len() != b.data.len() {
        return None;
    }
    let mut data = Vec::with_capacity(a.data.len());
    match (coords, a.size) {
        | (None, Some(n)) => {
            let ab = matprod(cx, &a.data, &b.data, n);
            let ba = matprod(cx, &b.data, &a.data, n);
            let minus_one = cx.graph.int(-1);
            for (x, y) in ab.into_iter().zip(ba) {
                let neg_y = mul(cx.graph, &[minus_one, y]);
                let diff = add(cx.graph, &[x, neg_y]);
                data.push(cx.simplify(diff));
            }
        },
        | (Some(c), None) => {
            let minus_one = cx.graph.int(-1);
            for i in 0..c.len() {
                let mut terms = Vec::new();
                for (j, &x) in c.iter().enumerate() {
                    let db = derivative(cx.graph, b.data[i], x)?;
                    let da = derivative(cx.graph, a.data[i], x)?;
                    terms.push(mul(cx.graph, &[a.data[j], db]));
                    terms.push(mul(cx.graph, &[minus_one, b.data[j], da]));
                }
                let sum = add(cx.graph, &terms);
                data.push(cx.simplify(sum));
            }
        },
        | _ => return None,
    }
    Some(Elem { size: a.size, data })
}

fn is_zero_elem(
    cx: &mut Cx<'_>,
    e: &Elem,
) -> bool {
    e.data.iter().all(|&d| vanishes(cx, d))
}

/// Number of sample points used to solve for constant coefficients.
const SAMPLES: usize = 10;

/// The coordinates of each target in the basis (free coordinates of a
/// dependent basis are zero), with the rank of the basis; `None` when some
/// target is outside the span.
fn express(
    cx: &mut Cx<'_>,
    basis: &[Elem],
    targets: &[Elem],
    coords: Option<&[NodeId]>,
) -> Option<(Vec<Vec<NodeId>>, usize)> {
    let d = basis.len();
    let len = basis.first()?.data.len();
    // Vector fields: the coefficients are constants, so the equations are
    // collected at several sample points of the coordinates; the result is
    // verified symbolically below.
    let mut m: Mat = Vec::new();
    for k in 0..if coords.is_some() { SAMPLES } else { 1 } {
        for r in 0..len {
            let mut row = Vec::with_capacity(d + targets.len());
            for e in basis.iter().chain(targets) {
                let mut t = e.data[r];
                for (j, &x) in coords.unwrap_or_default().iter().enumerate() {
                    let v = cx.graph.int(i64::try_from((k * 7 + j * 13 + k * j * 3) % 17).unwrap_or(0) - 8);
                    t = cx.graph.substitute(t, x, v);
                }
                row.push(if coords.is_some() { cx.simplify(t) } else { t });
            }
            m.push(row);
        }
    }
    let len = m.len();
    let minus_one = cx.graph.int(-1);
    let mut pivots = Vec::new();
    let mut r = 0;
    for c in 0..d {
        let Some(p) = (r..len).find(|&p| !is_zero_literal(cx.graph, m[p][c])) else {
            continue;
        };
        m.swap(r, p);
        let inv = pow_int(cx.graph, m[r][c], -1);
        for entry in &mut m[r] {
            let scaled = mul(cx.graph, &[*entry, inv]);
            *entry = cx.simplify(scaled);
        }
        for q in 0..len {
            if q == r || is_zero_literal(cx.graph, m[q][c]) {
                continue;
            }
            let f = m[q][c];
            for k in 0..m[q].len() {
                let t = mul(cx.graph, &[minus_one, f, m[r][k]]);
                let sum = add(cx.graph, &[m[q][k], t]);
                m[q][k] = cx.simplify(sum);
            }
        }
        pivots.push(c);
        r += 1;
    }
    for row in &m[r..] {
        if !row[d..].iter().all(|&e| vanishes(cx, e)) {
            return None;
        }
    }
    let zero = cx.graph.int(0);
    let mut out = Vec::with_capacity(targets.len());
    for j in 0..targets.len() {
        let mut v = vec![zero; d];
        for (row, &c) in pivots.iter().enumerate() {
            v[c] = m[row][d + j];
        }
        out.push(v);
    }
    if coords.is_some() {
        for (v, target) in out.iter().zip(targets) {
            let parts: Vec<(NodeId, &Elem)> = v.iter().copied().zip(basis).collect();
            let minus_one = cx.graph.int(-1);
            let mut all = parts;
            all.push((minus_one, target));
            let residual = combine(cx, &all)?;
            if !is_zero_elem(cx, &residual) {
                return None;
            }
        }
    }
    Some((out, pivots.len()))
}

/// `[e_i, e_j]` expressed in the basis: `c[i][j][k]`. Requires a basis.
fn structure(
    cx: &mut Cx<'_>,
    basis: &[Elem],
    coords: Option<&[NodeId]>,
) -> Option<Vec<Vec<Vec<NodeId>>>> {
    let n = basis.len();
    let mut targets = Vec::with_capacity(n * n);
    for a in basis {
        for b in basis {
            targets.push(bracket(cx, a, b, coords)?);
        }
    }
    let (coefficients, rank) = express(cx, basis, &targets, coords)?;
    if rank != n {
        return None;
    }
    Some(coefficients.chunks(n).map(<[Vec<NodeId>]>::to_vec).collect())
}

/// `ad_X` in the basis: entry `(i, j)` is the `i`th coordinate of `[X, e_j]`.
fn adjoint(
    cx: &mut Cx<'_>,
    x: &Elem,
    basis: &[Elem],
    coords: Option<&[NodeId]>,
) -> Option<Mat> {
    let targets: Vec<Elem> = basis.iter().map(|b| bracket(cx, x, b, coords)).collect::<Option<_>>()?;
    let (columns, rank) = express(cx, basis, &targets, coords)?;
    if rank != basis.len() {
        return None;
    }
    Some((0..basis.len()).map(|i| columns.iter().map(|col| col[i]).collect()).collect())
}

fn jacobi(
    cx: &mut Cx<'_>,
    basis: &[Elem],
    coords: Option<&[NodeId]>,
) -> Option<bool> {
    let n = basis.len();
    for i in 0..n {
        for j in i + 1..n {
            for k in j + 1..n {
                let (x, y, z) = (&basis[i], &basis[j], &basis[k]);
                let mut parts = Vec::new();
                for (a, b, c) in [(x, y, z), (y, z, x), (z, x, y)] {
                    let inner = bracket(cx, b, c, coords)?;
                    parts.push(bracket(cx, a, &inner, coords)?);
                }
                let one = cx.graph.int(1);
                let sum = combine(cx, &[(one, &parts[0]), (one, &parts[1]), (one, &parts[2])])?;
                if !is_zero_elem(cx, &sum) {
                    return Some(false);
                }
            }
        }
    }
    Some(true)
}

// ----------------------------------------------------------------------
// The exponential map
// ----------------------------------------------------------------------

fn exp_of(
    cx: &mut Cx<'_>,
    x: NodeId,
) -> Option<NodeId> {
    let e = call(cx.graph, "exp", &[x])?;
    Some(cx.simplify(e))
}

fn identity(
    cx: &mut Cx<'_>,
    n: usize,
) -> Vec<NodeId> {
    (0..n * n).map(|k| cx.graph.int(i64::from(k / n == k % n))).collect()
}

/// `sum coef_k * data_k` with `data_k` flattened `n × n` matrices.
fn series(
    cx: &mut Cx<'_>,
    coefs: &[NodeId],
    powers: &[Vec<NodeId>],
) -> Vec<NodeId> {
    (0..powers[0].len())
        .map(|at| {
            let terms: Vec<NodeId> = coefs.iter().zip(powers).map(|(&c, p)| mul(cx.graph, &[c, p[at]])).collect();
            let sum = add(cx.graph, &terms);
            cx.simplify(sum)
        })
        .collect()
}

fn factorial_inverse(
    graph: &mut Graph,
    k: u32,
) -> NodeId {
    let f: BigInt = (1..=k).map(BigInt::from).product();
    graph.num(Number::rat(BigRational::new(BigInt::from(1), f)))
}

fn taylor(
    cx: &mut Cx<'_>,
    a: &[NodeId],
    n: usize,
    order: u32,
) -> Vec<NodeId> {
    let mut powers = vec![identity(cx, n)];
    for _ in 0..order {
        let last = powers.last().cloned().unwrap_or_default();
        powers.push(matprod(cx, &last, a, n));
    }
    let coefs: Vec<NodeId> = (0..=order).map(|k| factorial_inverse(cx.graph, k)).collect();
    series(cx, &coefs, &powers)
}

fn nilpotent_exp(
    cx: &mut Cx<'_>,
    a: &[NodeId],
    n: usize,
) -> Option<Vec<NodeId>> {
    let mut p = a.to_vec();
    for k in 1..=n {
        if p.iter().all(|&e| is_zero_literal(cx.graph, e)) {
            let order = u32::try_from(k).ok()?;
            return Some(taylor(cx, a, n, order.saturating_sub(1)));
        }
        p = matprod(cx, &p, a, n);
    }
    None
}

fn diagonal_exp(
    cx: &mut Cx<'_>,
    a: &[NodeId],
    n: usize,
) -> Option<Vec<NodeId>> {
    let zero = cx.graph.int(0);
    let mut out = vec![zero; n * n];
    for r in 0..n {
        for c in 0..n {
            if r != c && !is_zero_literal(cx.graph, a[r * n + c]) {
                return None;
            }
        }
        out[r * n + r] = exp_of(cx, a[r * n + r])?;
    }
    Some(out)
}

fn two_by_two_exp(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<Vec<NodeId>> {
    let half = cx.graph.num(Number::fraction(1, 2)?);
    let minus_one = cx.graph.int(-1);
    let sum = add(cx.graph, &[a[0], a[3]]);
    let m = mul(cx.graph, &[half, sum]);
    let m = cx.simplify(m);
    let shift = |cx: &mut Cx<'_>, entry: NodeId, on_diagonal: bool| {
        let t = if on_diagonal {
            let neg_m = mul(cx.graph, &[minus_one, m]);
            add(cx.graph, &[entry, neg_m])
        } else {
            entry
        };
        cx.simplify(t)
    };
    let b = [shift(cx, a[0], true), shift(cx, a[1], false), shift(cx, a[2], false), shift(cx, a[3], true)];
    let b2 = matprod(cx, &b, &b, 2);
    // B² = delta I for a traceless 2x2 matrix.
    let (c, s) = exp_even(cx, b2[0])?;
    let scale = exp_of(cx, m)?;
    let id = identity(cx, 2);
    let mut out = Vec::with_capacity(4);
    for k in 0..4 {
        let ci = mul(cx.graph, &[c, id[k]]);
        let sb = mul(cx.graph, &[s, b[k]]);
        let t = add(cx.graph, &[ci, sb]);
        let t = mul(cx.graph, &[scale, t]);
        out.push(cx.simplify(t));
    }
    Some(out)
}

/// Rodrigues' formula for `A = [w]_x`.
fn rodrigues(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<Vec<NodeId>> {
    for r in 0..3 {
        for c in 0..3 {
            let sum = add(cx.graph, &[a[r * 3 + c], a[c * 3 + r]]);
            if !vanishes(cx, sum) {
                return None;
            }
        }
    }
    let squares: Vec<NodeId> = [a[7], a[2], a[3]].iter().map(|&x| pow_int(cx.graph, x, 2)).collect();
    let theta2 = add(cx.graph, &squares);
    let theta2 = cx.simplify(theta2);
    // A³ = -θ² A, hence e^A = I + sin(θ)/θ A + (1 - cos θ)/θ² A².
    let a2 = matprod(cx, a, a, 3);
    let half = cx.graph.num(Number::fraction(1, 2)?);
    let w = match exact_sqrt(cx.graph, theta2) {
        | Some(w) => w,
        | None => cx.graph.node(core::POW, &[theta2, half]),
    };
    let sin = call(cx.graph, "sin", &[w])?;
    let cos = call(cx.graph, "cos", &[w])?;
    let w_inv = pow_int(cx.graph, w, -1);
    let c1 = mul(cx.graph, &[sin, w_inv]);
    let one = cx.graph.int(1);
    let minus_one = cx.graph.int(-1);
    let neg_cos = mul(cx.graph, &[minus_one, cos]);
    let one_minus_cos = add(cx.graph, &[one, neg_cos]);
    let c2 = mul(cx.graph, &[one_minus_cos, w_inv, w_inv]);
    let (c1, c2) = (cx.simplify(c1), cx.simplify(c2));
    let id = identity(cx, 3);
    Some(series(cx, &[one, c1, c2], &[id, a.to_vec(), a2]))
}

fn eigen_exp(
    cx: &mut Cx<'_>,
    rows: &Mat,
) -> Option<Vec<NodeId>> {
    let n = rows.len();
    let a_term = matrix_term(cx.graph, rows);
    let request = call(cx.graph, "eigenvects", &[a_term])?;
    let result = cx.simplify(request);
    let entries = items(cx.graph, result)?;
    let mut columns: Vec<Vec<NodeId>> = Vec::new();
    let mut values = Vec::new();
    for entry in entries {
        let &[value, _, basis] = cx.graph.children(entry) else {
            return None;
        };
        for v in items(cx.graph, basis)? {
            columns.push(items(cx.graph, v)?);
            values.push(value);
        }
    }
    if columns.len() != n || columns.iter().any(|c| c.len() != n) {
        return None;
    }
    let p: Mat = (0..n).map(|r| columns.iter().map(|c| c[r]).collect()).collect();
    let p_term = matrix_term(cx.graph, &p);
    let inverse = call(cx.graph, "inverse", &[p_term])?;
    let inverse = cx.simplify(inverse);
    let p_inv: Vec<NodeId> = matrix(cx.graph, inverse)?.into_iter().flatten().collect();
    let zero = cx.graph.int(0);
    let mut d = vec![zero; n * n];
    for (k, &v) in values.iter().enumerate() {
        d[k * n + k] = exp_of(cx, v)?;
    }
    let flat: Vec<NodeId> = p.into_iter().flatten().collect();
    let pd = matprod(cx, &flat, &d, n);
    Some(matprod(cx, &pd, &p_inv, n))
}

fn exp_closed(
    cx: &mut Cx<'_>,
    rows: &Mat,
) -> Option<Vec<NodeId>> {
    let n = rows.len();
    let a: Vec<NodeId> = rows.iter().flatten().copied().collect();
    if let Some(e) = nilpotent_exp(cx, &a, n) {
        return Some(e);
    }
    if let Some(e) = diagonal_exp(cx, &a, n) {
        return Some(e);
    }
    if n == 2 {
        return two_by_two_exp(cx, &a);
    }
    if n == 3 {
        if let Some(e) = rodrigues(cx, &a) {
            return Some(e);
        }
    }
    eigen_exp(cx, rows)
}

fn bch(
    cx: &mut Cx<'_>,
    x: &Elem,
    y: &Elem,
    order: i64,
    coords: Option<&[NodeId]>,
) -> Option<Elem> {
    let one = cx.graph.int(1);
    let frac = |cx: &mut Cx<'_>, n: i64, d: i64| Number::fraction(n, d).map(|f| cx.graph.num(f));
    let xy = bracket(cx, x, y, coords)?;
    let mut parts: Vec<(NodeId, Elem)> = vec![(one, x.clone()), (one, y.clone())];
    if order >= 2 {
        parts.push((frac(cx, 1, 2)?, xy.clone()));
    }
    if order >= 3 {
        let xxy = bracket(cx, x, &xy, coords)?;
        let yx = bracket(cx, y, x, coords)?;
        let yyx = bracket(cx, y, &yx, coords)?;
        let twelfth = frac(cx, 1, 12)?;
        parts.push((twelfth, xxy.clone()));
        parts.push((twelfth, yyx));
        if order >= 4 {
            let yxxy = bracket(cx, y, &xxy, coords)?;
            parts.push((frac(cx, -1, 24)?, yxxy));
        }
    }
    let refs: Vec<(NodeId, &Elem)> = parts.iter().map(|(c, e)| (*c, e)).collect();
    combine(cx, &refs)
}

// ----------------------------------------------------------------------
// The kernel
// ----------------------------------------------------------------------

struct LieKernel {
    op: OpId,
    request: Request,
}

impl Kernel for LieKernel {
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

fn truth(
    graph: &mut Graph,
    verdict: bool,
) -> Option<NodeId> {
    let op = graph.ops().lookup(if verdict { "true" } else { "false" })?;
    Some(graph.node(op, &[]))
}

const SO3: &str = "list(list(list(0, 0, 0), list(0, 0, -1), list(0, 1, 0)), list(list(0, 0, 1), list(0, 0, 0), list(-1, 0, 0)), list(list(0, -1, 0), list(1, 0, 0), list(0, 0, 0)))";
const SU2: &str = "list(list(list(0, -I/2), list(-I/2, 0)), list(list(0, -1/2), list(1/2, 0)), list(list(-I/2, 0), list(0, I/2)))";
const SL2: &str = "list(list(list(1, 0), list(0, -1)), list(list(0, 1), list(0, 0)), list(list(0, 0), list(1, 0)))";

impl LieKernel {
    #[allow(clippy::too_many_lines)]
    fn compute(
        &self,
        cx: &mut Cx<'_>,
        args: &[NodeId],
    ) -> Option<NodeId> {
        match self.request {
            | Request::So3 | Request::Su2 | Request::Sl2 => {
                if !args.is_empty() {
                    return None;
                }
                let src = match self.request {
                    | Request::So3 => SO3,
                    | Request::Su2 => SU2,
                    | _ => SL2,
                };
                let term = cx.graph.parse(src).ok()?;
                Some(cx.simplify(term))
            },
            | Request::Bracket => {
                let (head, coords) = with_coords(cx.graph, args, 2)?;
                let a = read_elem(cx.graph, head[0], coords.as_deref())?;
                let b = read_elem(cx.graph, head[1], coords.as_deref())?;
                let c = bracket(cx, &a, &b, coords.as_deref())?;
                Some(elem_term(cx.graph, &c))
            },
            | Request::AdjointGroup => {
                let &[g, x] = args else {
                    return None;
                };
                let g_rows = matrix(cx.graph, g)?;
                let n = g_rows.len();
                let x = read_elem(cx.graph, x, None)?;
                if x.size != Some(n) || g_rows.len() != g_rows[0].len() {
                    return None;
                }
                let g_term = matrix_term(cx.graph, &g_rows);
                let inverse = call(cx.graph, "inverse", &[g_term])?;
                let inverse = cx.simplify(inverse);
                let g_inv: Vec<NodeId> = matrix(cx.graph, inverse)?.into_iter().flatten().collect();
                let g_flat: Vec<NodeId> = g_rows.into_iter().flatten().collect();
                let gx = matprod(cx, &g_flat, &x.data, n);
                let data = matprod(cx, &gx, &g_inv, n);
                Some(elem_term(cx.graph, &Elem { size: Some(n), data }))
            },
            | Request::ExpMap => {
                let rows = matrix(cx.graph, *args.first()?)?;
                if rows.len() != rows[0].len() {
                    return None;
                }
                let n = rows.len();
                let data = match args {
                    | [_] => exp_closed(cx, &rows)?,
                    | &[_, order] => {
                        let order = u32::try_from(cx.graph.number_of(order)?.to_i64()?).ok()?;
                        if order > 64 {
                            return None;
                        }
                        let a: Vec<NodeId> = rows.into_iter().flatten().collect();
                        taylor(cx, &a, n, order)
                    },
                    | _ => return None,
                };
                Some(elem_term(cx.graph, &Elem { size: Some(n), data }))
            },
            | Request::Bch => {
                let (head, coords) = with_coords(cx.graph, args, 3)?;
                let x = read_elem(cx.graph, head[0], coords.as_deref())?;
                let y = read_elem(cx.graph, head[1], coords.as_deref())?;
                let order = cx.graph.number_of(head[2])?.to_i64()?;
                if !(1..=4).contains(&order) {
                    return None;
                }
                let r = bch(cx, &x, &y, order, coords.as_deref())?;
                Some(elem_term(cx.graph, &r))
            },
            | Request::AdjointRep => {
                let (head, coords) = with_coords(cx.graph, args, 2)?;
                let x = read_elem(cx.graph, head[0], coords.as_deref())?;
                let basis = read_basis(cx.graph, head[1], coords.as_deref())?;
                let m = adjoint(cx, &x, &basis, coords.as_deref())?;
                Some(matrix_term(cx.graph, &m))
            },
            | Request::StructureConstants | Request::Killing | Request::CommutatorTable | Request::Jacobi | Request::IsLieAlgebra => {
                let (head, coords) = with_coords(cx.graph, args, 1)?;
                let coords = coords.as_deref();
                let basis = read_basis(cx.graph, head[0], coords)?;
                self.on_basis(cx, &basis, coords)
            },
        }
    }

    fn on_basis(
        &self,
        cx: &mut Cx<'_>,
        basis: &[Elem],
        coords: Option<&[NodeId]>,
    ) -> Option<NodeId> {
        let n = basis.len();
        match self.request {
            | Request::StructureConstants => {
                let c = structure(cx, basis, coords)?;
                let planes: Vec<NodeId> = c
                    .iter()
                    .map(|plane| {
                        let rows: Vec<NodeId> = plane.iter().map(|v| list(cx.graph, v)).collect();
                        list(cx.graph, &rows)
                    })
                    .collect();
                Some(list(cx.graph, &planes))
            },
            | Request::Killing => {
                let c = structure(cx, basis, coords)?;
                // ad_a[k][j] = c[a][j][k]; K_ab = sum_jk ad_a[k][j] ad_b[j][k].
                let rows: Vec<Vec<NodeId>> = (0..n)
                    .map(|a| {
                        (0..n)
                            .map(|b| {
                                let mut terms = Vec::with_capacity(n * n);
                                for (j, row) in c[a].iter().enumerate() {
                                    for (k, &x) in row.iter().enumerate() {
                                        terms.push(mul(cx.graph, &[x, c[b][k][j]]));
                                    }
                                }
                                let sum = add(cx.graph, &terms);
                                cx.simplify(sum)
                            })
                            .collect()
                    })
                    .collect();
                Some(matrix_term(cx.graph, &rows))
            },
            | Request::CommutatorTable => {
                let mut rows = Vec::with_capacity(n);
                for a in basis {
                    let mut row = Vec::with_capacity(n);
                    for b in basis {
                        let c = bracket(cx, a, b, coords)?;
                        row.push(elem_term(cx.graph, &c));
                    }
                    rows.push(list(cx.graph, &row));
                }
                Some(list(cx.graph, &rows))
            },
            | Request::Jacobi => {
                let verdict = jacobi(cx, basis, coords)?;
                truth(cx.graph, verdict)
            },
            | _ => {
                let mut targets = Vec::with_capacity(n * n);
                for a in basis {
                    for b in basis {
                        targets.push(bracket(cx, a, b, coords)?);
                    }
                }
                let closed = express(cx, basis, &targets, coords).is_some();
                let verdict = closed && jacobi(cx, basis, coords)?;
                truth(cx.graph, verdict)
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rules::testing::eval;
    use crate::rules::testing::reduce_with;
    use crate::rules::testing::simplify;

    fn run(src: &str) -> String {
        simplify(&[lie()], src)
    }

    const E: &str = "list(list(0, 1), list(0, 0))";
    const F: &str = "list(list(0, 0), list(1, 0))";

    #[test]
    fn standard_bases() {
        assert_eq!(run("so3_basis()"), "list(list(list(0, 0, 0), list(0, 0, -1), list(0, 1, 0)), list(list(0, 0, 1), list(0, 0, 0), list(-1, 0, 0)), list(list(0, -1, 0), list(1, 0, 0), list(0, 0, 0)))");
        assert_eq!(run("su2_basis()"), "list(list(list(0, -1/2*I), list(-1/2*I, 0)), list(list(0, -1/2), list(1/2, 0)), list(list(-1/2*I, 0), list(0, 1/2*I)))");
        assert_eq!(run("sl2_basis()"), "list(list(list(1, 0), list(0, -1)), list(list(0, 1), list(0, 0)), list(list(0, 0), list(1, 0)))");
    }

    #[test]
    fn matrix_bracket() {
        assert_eq!(run(&format!("lie_bracket({E}, {F})")), "list(list(1, 0), list(0, -1))");
        assert_eq!(run("lie_bracket(list(list(a, b), list(c, d)), list(list(a, b), list(c, d)))"), "list(list(0, 0), list(0, 0))");
        assert_eq!(run("lie_bracket(list(list(0, x), list(0, 0)), list(list(0, 0), list(y, 0)))"), "list(list(x*y, 0), list(0, -x*y))");
        // not square: left alone
        assert!(!reduce_with(&[lie()], "lie_bracket(list(list(1, 2)), list(list(1, 2)))", &[]).1);
    }

    #[test]
    fn vector_field_bracket() {
        // [-y d_x + x d_y... ] rotation field against d_y
        assert_eq!(run("lie_bracket(list(y, -x), list(0, 1), list(x, y))"), "list(-1, 0)");
        // [x d_x, d_x] = -d_x
        assert_eq!(run("lie_bracket(list(x, 0), list(1, 0), list(x, y))"), "list(-1, 0)");
        // rotations of R^3 close into so(3)
        let basis = "list(list(0, -z, y), list(z, 0, -x), list(-y, x, 0))";
        assert_eq!(run(&format!("is_lie_algebra({basis}, list(x, y, z))")), "true");
        assert_eq!(run(&format!("check_jacobi({basis}, list(x, y, z))")), "true");
        // [Lx, Ly] = -Lz for these fields, so the constants are minus those of so3_basis().
        let c = run(&format!("structure_constants({basis}, list(x, y, z))"));
        assert_eq!(entry(&c, &[0, 1]), "list(0, 0, -1)");
        assert_eq!(entry(&c, &[1, 2]), "list(-1, 0, 0)");
        // d_x, x d_x, x^2 d_x close (sl(2) acting on the line); x d_x, x^2 d_x, x^3 d_x do not
        assert_eq!(run("is_lie_algebra(list(list(1), list(x), list(x^2)), list(x))"), "true");
        assert_eq!(run("is_lie_algebra(list(list(x), list(x^2), list(x^3)), list(x))"), "false");
    }

    #[test]
    fn so3_structure_constants_are_levi_civita() {
        let c = run("structure_constants(so3_basis())");
        let epsilon = |i: usize, j: usize, k: usize| -> i32 {
            match (i, j, k) {
                | (0, 1, 2) | (1, 2, 0) | (2, 0, 1) => 1,
                | (0, 2, 1) | (2, 1, 0) | (1, 0, 2) => -1,
                | _ => 0,
            }
        };
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    let v = eval(&[lie()], &format!("1*{}", entry(&c, &[i, j, k])), &[]);
                    assert!((v - f64::from(epsilon(i, j, k))).abs() < 1e-12, "c[{i}][{j}][{k}] = {v}");
                }
            }
        }
        // su(2) in the -i sigma/2 convention has the same constants.
        assert_eq!(run("structure_constants(su2_basis())"), c);
    }

    /// The entry at `path` of a printed nested `list(...)`.
    fn entry(text: &str, path: &[usize]) -> String {
        let mut cur = text.to_owned();
        for &p in path {
            let inner = cur.strip_prefix("list(").and_then(|s| s.strip_suffix(')')).unwrap_or(&cur);
            let mut depth = 0;
            let mut parts = vec![String::new()];
            for ch in inner.chars() {
                match ch {
                    | '(' => depth += 1,
                    | ')' => depth -= 1,
                    | _ => {},
                }
                if ch == ',' && depth == 0 {
                    parts.push(String::new());
                } else if let Some(last) = parts.last_mut() {
                    last.push(ch);
                }
            }
            cur = parts[p].trim().to_owned();
        }
        cur
    }

    #[test]
    fn sl2_constants_killing_and_adjoint() {
        assert_eq!(run("structure_constants(sl2_basis())"), "list(list(list(0, 0, 0), list(0, 2, 0), list(0, 0, -2)), list(list(0, -2, 0), list(0, 0, 0), list(1, 0, 0)), list(list(0, 0, 2), list(-1, 0, 0), list(0, 0, 0)))");
        assert_eq!(run("killing_form(so3_basis())"), "list(list(-2, 0, 0), list(0, -2, 0), list(0, 0, -2))");
        assert_eq!(run("killing_form(sl2_basis())"), "list(list(8, 0, 0), list(0, 0, 4), list(0, 4, 0))");
        // ad_h = diag(0, 2, -2) in the basis (h, e, f)
        assert_eq!(run("adjoint_rep(list(list(1, 0), list(0, -1)), sl2_basis())"), "list(list(0, 0, 0), list(0, 2, 0), list(0, 0, -2))");
        // ad_e: [e, h] = -2e, [e, e] = 0, [e, f] = h
        assert_eq!(run(&format!("adjoint_rep({E}, sl2_basis())")), "list(list(0, 0, 1), list(-2, 0, 0), list(0, 0, 0))");
    }

    #[test]
    fn jacobi_and_closure() {
        assert_eq!(run("check_jacobi(so3_basis())"), "true");
        assert_eq!(run("check_jacobi(su2_basis())"), "true");
        assert_eq!(run("check_jacobi(sl2_basis())"), "true");
        assert_eq!(run("is_lie_algebra(so3_basis())"), "true");
        assert_eq!(run("is_lie_algebra(su2_basis())"), "true");
        assert_eq!(run("is_lie_algebra(sl2_basis())"), "true");
        // e, f alone do not close ([e, f] = h)
        assert_eq!(run(&format!("is_lie_algebra(list({E}, {F}))")), "false");
        // strictly upper triangular 3x3 (Heisenberg) closes
        assert_eq!(run("is_lie_algebra(list(list(list(0,1,0),list(0,0,0),list(0,0,0)), list(list(0,0,0),list(0,0,1),list(0,0,0)), list(list(0,0,1),list(0,0,0),list(0,0,0))))"), "true");
        // a symbolic algebra: span of diagonal matrices is abelian
        assert_eq!(run("check_jacobi(list(list(list(a, 0), list(0, b)), list(list(c, 0), list(0, d))))"), "true");
    }

    #[test]
    fn commutator_table_of_so3() {
        let t = run("commutator_table(so3_basis())");
        assert_eq!(entry(&t, &[0, 0]), "list(list(0, 0, 0), list(0, 0, 0), list(0, 0, 0))");
        // [Lx, Ly] = Lz
        assert_eq!(entry(&t, &[0, 1]), "list(list(0, -1, 0), list(1, 0, 0), list(0, 0, 0))");
        assert_eq!(entry(&t, &[1, 0]), "list(list(0, 1, 0), list(-1, 0, 0), list(0, 0, 0))");
    }

    #[test]
    fn adjoint_of_the_group() {
        // g E g^-1 for g = [[1, 1], [0, 1]]
        assert_eq!(run("adjoint_group(list(list(1, 1), list(0, 1)), list(list(0, 0), list(1, 0)))"), "list(list(1, -1), list(1, -1))");
        assert_eq!(run("adjoint_group(list(list(2, 0), list(0, 1)), list(list(0, 1), list(0, 0)))"), "list(list(0, 2), list(0, 0))");
    }

    #[test]
    fn exponential_of_so3_generators_is_a_rotation() {
        assert_eq!(run("exp_map(list(list(0, 0, 0), list(0, 0, -t), list(0, t, 0)))"), "list(list(1, 0, 0), list(0, cos(t), -sin(t)), list(0, sin(t), cos(t)))");
        assert_eq!(run("exp_map(list(list(0, -t), list(t, 0)))"), "list(list(cos(t), -sin(t)), list(sin(t), cos(t)))");
        // rotation about the axis (1, 2, 2) by angle 3: trace 1 + 2 cos 3, R R^T = I
        let a = "list(list(0, -2, 2), list(2, 0, -1), list(-2, 1, 0))";
        let tr = run(&format!("trace(exp_map({a}))"));
        assert!((eval(&[lie()], &tr, &[]) - (1.0 + 2.0 * 3.0_f64.cos())).abs() < 1e-12, "{tr}");
        let rrt = run(&format!("matmul(exp_map({a}), transpose(exp_map({a})))"));
        for (i, j) in [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2)] {
            let v = eval(&[lie()], &entry(&rrt, &[i, j]), &[]);
            assert!((v - f64::from(u8::from(i == j))).abs() < 1e-12, "{rrt}");
        }
        let det = run(&format!("det(exp_map({a}))"));
        assert!((eval(&[lie()], &det, &[]) - 1.0).abs() < 1e-12, "{det}");
    }

    #[test]
    fn exponential_closed_forms() {
        assert_eq!(run(&format!("exp_map({E})")), "list(list(1, 1), list(0, 1))");
        assert_eq!(run("exp_map(list(list(0, a, b), list(0, 0, c), list(0, 0, 0)))"), run("exp_map(list(list(0, a, b), list(0, 0, c), list(0, 0, 0)))"));
        assert_eq!(run("exp_map(list(list(a, 0), list(0, b)))"), "list(list(exp(a), 0), list(0, exp(b)))");
        assert_eq!(run("exp_map(list(list(0, 0), list(0, 0)))"), "list(list(1, 0), list(0, 1))");
        // hyperbolic: exp of [[0, t], [t, 0]]
        assert_eq!(run("exp_map(list(list(0, t), list(t, 0)))"), "list(list(cosh(t), sinh(t)), list(sinh(t), cosh(t)))");
        // a 2x2 with a double eigenvalue: exp([[1, 1], [0, 1]]) = e [[1, 1], [0, 1]]
        assert_eq!(run("exp_map(list(list(1, 1), list(0, 1)))"), "list(list(E, E), list(0, E))");
        // diagonalizable 3x3 (triangular, distinct eigenvalues)
        let m = run("exp_map(list(list(1, 2, 0), list(0, 3, 0), list(0, 0, 5)))");
        assert_eq!(entry(&m, &[2, 2]), "exp(5)");
        // the Taylor polynomial
        assert_eq!(run("exp_map(list(list(0, 1), list(0, 0)), 3)"), "list(list(1, 1), list(0, 1))");
        assert_eq!(run("exp_map(list(list(0, -1), list(1, 0)), 2)"), "list(list(1/2, -1), list(1, 1/2))");
        // trace of exp = exp of the trace for su(2) elements -> det = 1
        let d = run("det(exp_map(list(list(a, b), list(c, -a))))");
        let v = eval(&[lie()], &d, &[("a", 0.3), ("b", 0.2), ("c", -0.4)]);
        assert!((v - 1.0).abs() < 1e-12, "{d}");
        // not closed-form: left unreduced
        assert!(!reduce_with(&[lie()], "exp_map(list(list(0, 1, 0, 0), list(0, 0, 1, 0), list(0, 0, 0, 1), list(1, 0, 0, 0)))", &[]).1);
    }

    #[test]
    fn baker_campbell_hausdorff_terms() {
        // so(3): X = Lx, Y = Ly: order 1 = X + Y, order 2 adds Lz / 2
        let x = "list(list(0, 0, 0), list(0, 0, -1), list(0, 1, 0))";
        let y = "list(list(0, 0, 1), list(0, 0, 0), list(-1, 0, 0))";
        assert_eq!(run(&format!("baker_campbell_hausdorff({x}, {y}, 1)")), "list(list(0, 0, 1), list(0, 0, -1), list(-1, 1, 0))");
        assert_eq!(run(&format!("baker_campbell_hausdorff({x}, {y}, 2)")), "list(list(0, -1/2, 1), list(1/2, 0, -1), list(-1, 1, 0))");
        // sl(2): X = e, Y = f
        assert_eq!(run(&format!("baker_campbell_hausdorff({E}, {F}, 2)")), "list(list(1/2, 1), list(1, -1/2))");
        assert_eq!(run(&format!("baker_campbell_hausdorff({E}, {F}, 3)")), "list(list(1/2, 5/6), list(5/6, -1/2))");
        assert_eq!(run(&format!("baker_campbell_hausdorff({E}, {F}, 4)")), "list(list(5/12, 5/6), list(5/6, -5/12))");
        // commuting elements: X + Y at every order
        assert_eq!(run("baker_campbell_hausdorff(list(list(a, 0), list(0, b)), list(list(c, 0), list(0, d)), 4)"), "list(list(a + c, 0), list(0, b + d))");
        // Heisenberg group: the series stops after order 2
        let hx = "list(list(0, 1, 0), list(0, 0, 0), list(0, 0, 0))";
        let hy = "list(list(0, 0, 0), list(0, 0, 1), list(0, 0, 0))";
        assert_eq!(run(&format!("matmul(exp_map({hx}), exp_map({hy}))")), run(&format!("exp_map(baker_campbell_hausdorff({hx}, {hy}, 3))")));
        // fields: [d_x, x d_y] = d_y, so X + Y + [X, Y] / 2
        assert_eq!(run("baker_campbell_hausdorff(list(1, 0), list(0, x), 2, list(x, y))"), "list(1, x + 1/2)");
        // order outside 1..=4 is not handled
        assert!(!reduce_with(&[lie()], &format!("baker_campbell_hausdorff({E}, {F}, 5)"), &[]).1);
    }
}
