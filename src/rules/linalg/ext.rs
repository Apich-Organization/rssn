//! Exact matrix functions and factorisations: Jordan normal form, matrix
//! exponentials, powers and arbitrary functions of a matrix, minimal
//! polynomial, Cayley–Hamilton, pseudo-inverse and least squares, LDL and
//! Cholesky, Smith and Hermite normal forms over the integers, spaces of a
//! matrix, Kronecker products and special matrices.
//!
//! Entries may be exact rationals or symbolic terms; the algorithms use the
//! same zero-testing arithmetic as `rref` and `det`, so the results are
//! exact for generic symbolic entries.
//!
//! | operator | value |
//! |---|---|
//! | `kron(A, B)` | the Kronecker product |
//! | `diag(list(d1, ...))`, `vandermonde(list(x1, ...))`, `hilbert(n)`, `companion(p, x)` | special matrices (`companion` of a monic-normalised polynomial) |
//! | `adjugate(A)`, `cofactors(A)` | the adjugate (classical adjoint) and the matrix of cofactors |
//! | `colspace(A)`, `rowspace(A)`, `left_nullspace(A)` | bases, as lists of vectors |
//! | `pinv(A)` | the Moore–Penrose pseudo-inverse, from a full-rank factorisation `A = C R` |
//! | `lstsq(A, b)` | the minimum-norm least-squares solution `pinv(A) b` |
//! | `minpoly(A, l)` | the minimal polynomial in `l` (from the first linear dependence among the powers of `A`) |
//! | `matpoly(p, x, A)` | the polynomial `p(x)` evaluated at the matrix `A` (Horner) |
//! | `cayley_hamilton(A)` | the characteristic polynomial evaluated at `A`: the zero matrix |
//! | `jordan(A)` | `list(P, J)` with `A = P J P^-1` and `J` in Jordan normal form (eigenvalues that have a closed form) |
//! | `jordan_blocks(A)` | `list(list(eigenvalue, size), ...)` |
//! | `matfun(f, x, A)` | `f(A)` for an expression `f` in `x`: `P f(J) P^-1` with the derivative entries `f^(k)(λ)/k!` in each Jordan block |
//! | `matexp(A, t)` | `exp(A t)` (the matrix exponential), exact via the Jordan form |
//! | `mpow(A, n)` | the matrix power: repeated squaring for an integer `n` (negative via the inverse), the Jordan form for a symbolic `n` |
//! | `matsqrt(A)` | a square root through `matfun(sqrt(x), x, A)` |
//! | `spectral(A)` | `list(list(eigenvalue, projector), ...)` for a diagonalisable matrix: `A = Σ λ P_λ` |
//! | `ldl(A)` | `list(L, D)` with `A = L diag(D) Lᵀ` for a symmetric matrix with non-vanishing pivots |
//! | `cholesky(A)` | `L` with `A = L Lᵀ` (positive pivots; for symbolic pivots they are assumed positive) |
//! | `orthogonalize(vs)`, `orthonormalize(vs)` | Gram–Schmidt on a list of vectors (dependent vectors dropped) |
//! | `smith(A)` | `list(U, D, V)` with `U A V = D` the Smith normal form of an integer matrix, `U` and `V` unimodular |
//! | `invariant_factors(A)` | the non-zero diagonal of the Smith normal form |
//! | `hermite_form(A)` | `list(H, U)` with `U A = H` the row-style Hermite normal form of an integer matrix |
//! | `is_symmetric(A)`, `is_orthogonal(A)`, `is_positive_definite(A)` | truth values (the last for exact rational matrices, by Sylvester's criterion) |

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::Zero;

use super::add;
use super::eigenvalues;
use super::identity_rows;
use super::inverse;
use super::list;
use super::matrix;
use super::matrix_term;
use super::mul;
use super::neg;
use super::normal;
use super::nullspace;
use super::reciprocal;
use super::reduce;
use super::request;
use super::sqrt;
use super::transpose;
use super::vector;
use super::Field;
use super::Terms;
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
    let mut out = super::matmul(cx.graph, a, b)?;
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
            let terms: Vec<NodeId> = row.iter().zip(v).map(|(&a, &b)| mul(cx.graph, &[a, b])).collect();
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
            let mut current: Vec<Vec<NodeId>> = if size >= 2 { kernels[size - 2].clone() } else { Vec::new() };
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
            let divisor = cx.graph.num(Number::rat(BigRational::new(BigInt::one(), factorial.clone())));
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
    let poly = crate::rules::poly::repr::from_term(cx.graph, &mut gens, p, crate::rules::poly::repr::Limits::default())?;
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
        let divisor = cx.graph.num(Number::rat(BigRational::new(BigInt::one(), factorial.clone())));
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
    let c: Mat = (0..rows).map(|i| pivots.iter().map(|&p| m[i][p]).collect()).collect();
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
    Some(pivots.iter().map(|&p| m.iter().map(|row| row[p]).collect()).collect())
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
            && !q.is_positive() {
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
            let numerator = dot(cx, v, b);
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
        let length_squared = dot(cx, &w, &w);
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

fn dot(
    cx: &mut Cx<'_>,
    a: &[NodeId],
    b: &[NodeId],
) -> NodeId {
    let terms: Vec<NodeId> = a.iter().zip(b).map(|(&x, &y)| mul(cx.graph, &[x, y])).collect();
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
                .map(|&e| match graph.number_of(e)? {
                    | Number::Int(v) => Some(v.clone()),
                    | _ => None,
                })
                .collect()
        })
        .collect()
}

fn integer_identity(n: usize) -> Ints {
    (0..n).map(|i| (0..n).map(|j| BigInt::from(i32::from(i == j))).collect()).collect()
}

fn integers_to_term(
    graph: &mut Graph,
    m: &Ints,
) -> NodeId {
    let rows: Mat = m.iter().map(|r| r.iter().map(|v| graph.num(Number::Int(v.clone()))).collect()).collect();
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
                    if !d[i][j].is_zero() && best.is_none_or(|(bi, bj)| d[i][j].abs() < d[bi][bj].abs()) {
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
        .map(|(_, r)| r.iter().enumerate().filter(|&(c, _)| c != col).map(|(_, &e)| e).collect())
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
            let value = super::determinant(cx, &sub)?;
            let value = if (i + j) % 2 == 0 { value } else { neg(cx.graph, value) };
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
    let poly = crate::rules::poly::repr::from_term(cx.graph, &mut gens, p, crate::rules::poly::repr::Limits::default())?;
    let mut degree = 0_usize;
    let mut current = poly;
    let mut coefficients = Vec::new();
    let zero = cx.graph.int(0);
    let mut factorial = BigInt::one();
    while !current.is_zero() {
        let term = crate::rules::poly::repr::to_term(cx.graph, &gens, &current);
        let at_zero = cx.graph.substitute(term, x, zero);
        let divisor = cx.graph.num(Number::rat(BigRational::new(BigInt::one(), factorial.clone())));
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
    let q: Vec<Vec<BigRational>> =
        m.iter().map(|r| r.iter().map(|&e| graph.number_of(e).and_then(Number::to_rational)).collect()).collect::<Option<_>>()?;
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
                let n = usize::try_from(cx.graph.number_of(arg(0)?)?.to_i64()?).ok().filter(|&n| (1..=64).contains(&n))?;
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
                let c = if matches!(self.kind, Kind::Adjugate) { transpose(&c) } else { c };
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
                let p = super::charpoly(cx, &m, l)?;
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
                let (u, d, v) =
                    (integers_to_term(cx.graph, &u), integers_to_term(cx.graph, &d), integers_to_term(cx.graph, &v));
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
                let (h, u) = (integers_to_term(cx.graph, &h), integers_to_term(cx.graph, &u));
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
        self.compute(cx, &args).map_or(Outcome::Pass, Outcome::Equal)
    }

    fn revisit(&self) -> bool {
        true
    }
}

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
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
        let op = i.op(OpDescriptor::new(name, Arity::Fixed(arity)).flags(OpFlags::HEAVY).cost(100))?;
        i.kernel(&format!("linalg/{name}"), Tier::Reduce, Extension { op, kind });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
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
        assert_eq!(run("diag(list(1, a, 3))"), "list(list(1, 0, 0), list(0, a, 0), list(0, 0, 3))");
        assert_eq!(run("vandermonde(list(1, 2, 3))"), "list(list(1, 1, 1), list(1, 2, 4), list(1, 3, 9))");
        assert_eq!(run("hilbert(2)"), "list(list(1, 1/2), list(1/2, 1/3))");
        assert_eq!(run("companion(x^3 - 2*x^2 + 3*x - 4, x)"), "list(list(0, 0, 4), list(1, 0, -3), list(0, 1, 2))");
        assert_eq!(run("det(vandermonde(list(a, b)))"), "b - a");
    }

    #[test]
    fn adjugate_and_spaces() {
        assert_eq!(run("adjugate(list(list(1, 2), list(3, 4)))"), "list(list(4, -2), list(-3, 1))");
        assert_eq!(run("cofactors(list(list(1, 2), list(3, 4)))"), "list(list(4, -3), list(-2, 1))");
        assert_eq!(run("colspace(list(list(1, 2, 3), list(2, 4, 6), list(1, 0, 1)))"), "list(list(1, 2, 1), list(2, 4, 0))");
        assert_eq!(run("rowspace(list(list(1, 2, 3), list(4, 5, 6)))"), "list(list(1, 0, -1), list(0, 1, 2))");
        assert_eq!(run("left_nullspace(list(list(1, 2), list(2, 4)))"), "list(list(-2, 1))");
    }

    #[test]
    fn pseudo_inverse_and_least_squares() {
        assert_eq!(run("pinv(list(list(1, 0), list(0, 2), list(0, 0)))"), "list(list(1, 0, 0), list(0, 1/2, 0))");
        // Rank deficient.
        assert_eq!(run("pinv(list(list(1, 1), list(1, 1)))"), "list(list(1/4, 1/4), list(1/4, 1/4))");
        // The line fit y = a + b x through (0,1), (1,2), (2,2): a = 7/6, b = 1/2.
        assert_eq!(run("lstsq(list(list(1, 0), list(1, 1), list(1, 2)), list(1, 2, 2))"), "list(7/6, 1/2)");
    }

    #[test]
    fn minimal_and_matrix_polynomials() {
        assert_eq!(run(&format!("minpoly({A}, l)")), "l^2 - 4*l + 4");
        assert_eq!(run("minpoly(list(list(2, 0), list(0, 2)), l)"), "l - 2");
        assert_eq!(run(&format!("cayley_hamilton({B})")), "list(list(0, 0, 0), list(0, 0, 0), list(0, 0, 0))");
        assert_eq!(run("cayley_hamilton(list(list(a, b), list(c, d)))"), "list(list(0, 0), list(0, 0))");
        assert_eq!(run(&format!("matpoly(x^2 - 4*x + 4, x, {A})")), "list(list(0, 0), list(0, 0))");
        assert_eq!(run("matpoly(x^2 + 1, x, list(list(0, -1), list(1, 0)))"), "list(list(0, 0), list(0, 0))");
    }

    #[test]
    fn jordan_form_and_matrix_functions() {
        assert_eq!(run(&format!("jordan_blocks({B})")), "list(list(4, 3))");
        assert_eq!(run("jordan_blocks(list(list(2, 0), list(0, 2)))"), "list(list(2, 1), list(2, 1))");
        // A P = P J.
        let pj = run(&format!("jordan({B})"));
        assert!(pj.starts_with("list(list("), "{pj}");
        let p = tests_part(&pj, 0);
        let j = tests_part(&pj, 1);
        assert_eq!(run(&format!("madd(matmul({B}, {p}), smul(-1, matmul({p}, {j})))")), "list(list(0, 0, 0), list(0, 0, 0), list(0, 0, 0))");
        // exp(A t) for a Jordan block: e^(2t) (I + t N).
        assert_eq!(run(&format!("matexp({A}, t)")), "list(list(exp(2*t), t*exp(2*t)), list(0, exp(2*t)))");
        assert_eq!(run("matexp(list(list(1, 0), list(0, 2)), t)"), "list(list(exp(t), 0), list(0, exp(2*t)))");
        assert_eq!(run(&format!("mpow({A}, 5)")), "list(list(32, 80), list(0, 32))");
        assert_eq!(run(&format!("mpow({A}, -1)")), "list(list(1/2, -1/4), list(0, 1/2))");
        assert_eq!(run(&format!("mpow({A}, n)")), "list(list(2^n, n*2^(n - 1)), list(0, 2^n))");
        assert_eq!(run("matsqrt(list(list(4, 0), list(0, 9)))"), "list(list(2, 0), list(0, 3))");
        assert_eq!(run(&format!("matfun(x^2 + 1, x, {A})")), "list(list(5, 4), list(0, 5))");
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
        assert_eq!(run("cholesky(list(list(4, 2), list(2, 3)))"), "list(list(2, 0), list(1, 2^(1/2)))");
        let (text, reduced) = crate::rules::testing::reduce_with(&[linalg()], "cholesky(list(list(1, 2), list(2, 1)))", &[]);
        assert!(!reduced, "not positive definite: {text}");
        assert_eq!(
            run("orthogonalize(list(list(1, 1, 0), list(1, 0, 1), list(0, 1, 1)))"),
            "list(list(1, 1, 0), list(1/2, -1/2, 1), list(-2/3, 2/3, 2/3))"
        );
        assert_eq!(run("orthonormalize(list(list(3, 4), list(1, 0)))"), "list(list(3/5, 4/5), list(4/5, -3/5))");
        assert_eq!(run("is_symmetric(list(list(1, 2), list(2, 5)))"), "true");
        assert_eq!(run("is_orthogonal(list(list(0, 1), list(-1, 0)))"), "true");
        assert_eq!(run("is_positive_definite(list(list(2, 1), list(1, 2)))"), "true");
        assert_eq!(run("is_positive_definite(list(list(1, 2), list(2, 1)))"), "false");
    }

    #[test]
    fn integer_normal_forms() {
        let a = "list(list(2, 4, 4), list(-6, 6, 12), list(10, 4, 16))";
        assert_eq!(run(&format!("invariant_factors({a})")), "list(2, 2, 156)");
        let smith = run(&format!("smith({a})"));
        let (u, d, v) = (tests_part(&smith, 0), tests_part(&smith, 1), tests_part(&smith, 2));
        assert_eq!(run(&format!("matmul({u}, {a}, {v})")), d);
        assert_eq!(run("det(list(list(1, 0), list(0, 1)))"), "1");
        let hermite = run("hermite_form(list(list(2, 3, 6, 2), list(5, 6, 1, 6), list(8, 3, 1, 1)))");
        let (h, u) = (tests_part(&hermite, 0), tests_part(&hermite, 1));
        assert_eq!(run(&format!("matmul({u}, list(list(2, 3, 6, 2), list(5, 6, 1, 6), list(8, 3, 1, 1)))")), h);
        assert_eq!(h, "list(list(1, 0, 50, -11), list(0, 3, 28, -2), list(0, 0, 61, -13))");
        assert_eq!(run("det(list(list(1, 2), list(3, 4)))"), "-2");
    }

    /// The `index`-th top-level item of a `list(...)` text.
    fn tests_part(
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
}
