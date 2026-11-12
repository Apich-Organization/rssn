//! Exact sparse linear algebra over an arbitrary coefficient domain.
//!
//! The dense exact routines ([`crate::kernels::qlinalg`], the symbolic
//! `linalg` rule set) cost `O(n^3)` coefficient operations whatever the
//! matrix looks like, and with exact coefficients each operation grows. A
//! sparse matrix of the sizes that arise from discretisations, incidence
//! structures or large structured models needs different algorithms, and
//! they come in three families, all here:
//!
//! 1. **Elimination with sparsity-preserving pivoting** ([`echelon`]), generic
//!    over a [`Domain`]: the rationals, a prime field, or *symbolic terms*
//!    (the `sparse_symbolic` rule set supplies a domain whose elements are
//!    expression-graph nodes, so the same code eliminates a matrix with
//!    polynomial entries). The pivot is chosen in a sparsest remaining row,
//!    on its column with the fewest remaining entries (a Markowitz count
//!    restricted to that row), which keeps fill-in low without a symbolic
//!    analysis phase.
//! 2. **Modular methods** for rational matrices, which never let a
//!    coefficient grow during elimination: [`det_rational`] combines
//!    determinants modulo word-sized primes by the Chinese remainder theorem
//!    up to Hadamard's bound, and [`solve_dixon`] solves `A x = b` by Dixon's
//!    `p`-adic lifting: one sparse factorisation modulo `p`, then one sparse
//!    solve and one sparse product per `p`-adic digit, and a rational
//!    reconstruction at the end. Both results are exact, not probabilistic:
//!    the determinant by the bound, the solution by substitution.
//! 3. **Black-box methods** over a prime field, which only multiply by the
//!    matrix and therefore create no fill-in at all: [`minpoly_wiedemann`]
//!    and [`solve_wiedemann`] (Wiedemann's algorithm with the
//!    Berlekamp-Massey recurrence solver).
//!
//! Matrices are [`SparseMat`]: one sorted list of `(column, value)` per row,
//! without explicit zeros. Indices are zero-based here; the term-level
//! operators use the kernel's one-based convention and convert.
//!
//! Degenerate input never panics: an index outside the declared shape is an
//! error from the constructor, and shape mismatches return `None`.

#![allow(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    clippy::many_single_char_names,
    clippy::module_name_repetitions,
    clippy::single_call_fn,
    clippy::missing_docs_in_private_items,
    reason = "index arithmetic over rows whose bounds are established by the constructors"
)]

use std::collections::BTreeMap;
use std::fmt::Debug;

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;

/// The arithmetic exact elimination needs.
///
/// Methods take `&mut self` because a domain may have state: the symbolic
/// domain builds nodes in an expression graph and runs a simplifier for its
/// zero test.
pub trait Domain {
    /// An element.
    type E: Clone + Debug;
    /// The additive identity.
    fn zero(&mut self) -> Self::E;
    /// The multiplicative identity.
    fn one(&mut self) -> Self::E;
    /// Whether `a` is zero. Must be exact: elimination divides by anything
    /// this reports as non-zero.
    fn is_zero(
        &mut self,
        a: &Self::E,
    ) -> bool;
    /// `a + b`.
    fn add(
        &mut self,
        a: &Self::E,
        b: &Self::E,
    ) -> Self::E;
    /// `a - b`.
    fn sub(
        &mut self,
        a: &Self::E,
        b: &Self::E,
    ) -> Self::E;
    /// `a * b`.
    fn mul(
        &mut self,
        a: &Self::E,
        b: &Self::E,
    ) -> Self::E;
    /// `-a`.
    fn neg(
        &mut self,
        a: &Self::E,
    ) -> Self::E {
        let zero = self.zero();
        self.sub(&zero, a)
    }
    /// `a / b`, or `None` when `b` is zero.
    fn div(
        &mut self,
        a: &Self::E,
        b: &Self::E,
    ) -> Option<Self::E>;
}

/// The rational numbers.
#[derive(Clone, Copy, Debug, Default)]
pub struct Rationals;

impl Domain for Rationals {
    type E = BigRational;

    fn zero(&mut self) -> BigRational {
        BigRational::zero()
    }

    fn one(&mut self) -> BigRational {
        BigRational::one()
    }

    fn is_zero(
        &mut self,
        a: &BigRational,
    ) -> bool {
        a.is_zero()
    }

    fn add(
        &mut self,
        a: &BigRational,
        b: &BigRational,
    ) -> BigRational {
        a + b
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
}

/// The field of integers modulo a prime below `2^63`; elements are reduced
/// residues.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Gfp(pub u64);

impl Gfp {
    /// `a * b mod p`.
    #[inline]
    #[must_use]
    pub fn mul_mod(
        self,
        a: u64,
        b: u64,
    ) -> u64 {
        ((u128::from(a) * u128::from(b)) % u128::from(self.0)) as u64
    }

    /// `a^e mod p`.
    #[must_use]
    pub fn pow_mod(
        self,
        mut a: u64,
        mut e: u64,
    ) -> u64 {
        let mut acc = 1 % self.0;
        a %= self.0;
        while e > 0 {
            if e & 1 == 1 {
                acc = self.mul_mod(acc, a);
            }
            a = self.mul_mod(a, a);
            e >>= 1;
        }
        acc
    }

    /// The inverse of `a` (Fermat), or `None` for zero.
    #[must_use]
    pub fn inv(
        self,
        a: u64,
    ) -> Option<u64> {
        (a % self.0 != 0).then(|| self.pow_mod(a, self.0 - 2))
    }

    /// The residue of an integer.
    #[must_use]
    pub fn reduce(
        self,
        a: &BigInt,
    ) -> u64 {
        a.mod_floor(&BigInt::from(self.0)).to_u64().unwrap_or(0)
    }
}

impl Domain for Gfp {
    type E = u64;

    fn zero(&mut self) -> u64 {
        0
    }

    fn one(&mut self) -> u64 {
        1 % self.0
    }

    fn is_zero(
        &mut self,
        a: &u64,
    ) -> bool {
        *a == 0
    }

    fn add(
        &mut self,
        a: &u64,
        b: &u64,
    ) -> u64 {
        ((u128::from(*a) + u128::from(*b)) % u128::from(self.0)) as u64
    }

    fn sub(
        &mut self,
        a: &u64,
        b: &u64,
    ) -> u64 {
        ((u128::from(*a) + u128::from(self.0) - u128::from(*b % self.0)) % u128::from(self.0))
            as u64
    }

    fn mul(
        &mut self,
        a: &u64,
        b: &u64,
    ) -> u64 {
        self.mul_mod(*a, *b)
    }

    fn div(
        &mut self,
        a: &u64,
        b: &u64,
    ) -> Option<u64> {
        self.inv(*b).map(|i| self.mul_mod(*a, i))
    }
}

/// A sparse matrix: per row, the non-zero entries sorted by column.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SparseMat<E> {
    /// Number of rows.
    pub rows: usize,
    /// Number of columns.
    pub cols: usize,
    /// `data[i]`: the `(column, value)` pairs of row `i`, columns ascending.
    pub data: Vec<Vec<(usize, E)>>,
}

impl<E: Clone + Debug> SparseMat<E> {
    /// The zero matrix.
    #[must_use]
    pub fn zeros(
        rows: usize,
        cols: usize,
    ) -> Self {
        Self {
            rows,
            cols,
            data: vec![Vec::new(); rows],
        }
    }

    /// Builds a matrix from `(row, column, value)` triplets. Repeated
    /// positions are added and zeros dropped.
    ///
    /// # Errors
    /// When a triplet lies outside the shape.
    pub fn from_triplets<D: Domain<E = E>>(
        dom: &mut D,
        rows: usize,
        cols: usize,
        triplets: impl IntoIterator<Item = (usize, usize, E)>,
    ) -> Result<Self, String> {
        let mut maps: Vec<BTreeMap<usize, E>> = vec![BTreeMap::new(); rows];
        for (i, j, v) in triplets {
            if i >= rows || j >= cols {
                return Err(format!(
                    "entry ({i}, {j}) is outside a {rows} x {cols} matrix"
                ));
            }
            let merged = match maps[i].get(&j) {
                | Some(old) => dom.add(old, &v),
                | None => v,
            };
            maps[i].insert(j, merged);
        }
        let mut data = Vec::with_capacity(rows);
        for map in maps {
            let mut row = Vec::with_capacity(map.len());
            for (j, v) in map {
                if !dom.is_zero(&v) {
                    row.push((j, v));
                }
            }
            data.push(row);
        }
        Ok(Self { rows, cols, data })
    }

    /// A sparse copy of a dense matrix.
    #[must_use]
    pub fn from_dense<D: Domain<E = E>>(
        dom: &mut D,
        dense: &[Vec<E>],
    ) -> Self {
        let cols = dense.iter().map(Vec::len).max().unwrap_or(0);
        let data = dense
            .iter()
            .map(|row| {
                row.iter()
                    .enumerate()
                    .filter(|(_, v)| !dom.is_zero(v))
                    .map(|(j, v)| (j, v.clone()))
                    .collect()
            })
            .collect();
        Self {
            rows: dense.len(),
            cols,
            data,
        }
    }

    /// The dense form.
    #[must_use]
    pub fn to_dense<D: Domain<E = E>>(
        &self,
        dom: &mut D,
    ) -> Vec<Vec<E>> {
        self.data
            .iter()
            .map(|row| {
                let mut out = vec![dom.zero(); self.cols];
                for (j, v) in row {
                    out[*j] = v.clone();
                }
                out
            })
            .collect()
    }

    /// Number of stored entries.
    #[must_use]
    pub fn nnz(&self) -> usize {
        self.data.iter().map(Vec::len).sum()
    }

    /// The entries as `(row, column, value)`, in row-major order.
    pub fn triplets(&self) -> impl Iterator<Item = (usize, usize, &E)> + '_ {
        self.data
            .iter()
            .enumerate()
            .flat_map(|(i, row)| row.iter().map(move |(j, v)| (i, *j, v)))
    }

    /// The transpose.
    #[must_use]
    pub fn transpose(&self) -> Self {
        let mut data: Vec<Vec<(usize, E)>> = vec![Vec::new(); self.cols];
        for (i, row) in self.data.iter().enumerate() {
            for (j, v) in row {
                data[*j].push((i, v.clone()));
            }
        }
        Self {
            rows: self.cols,
            cols: self.rows,
            data,
        }
    }

    /// `self + other`, or `None` when the shapes differ.
    #[must_use]
    pub fn add<D: Domain<E = E>>(
        &self,
        dom: &mut D,
        other: &Self,
    ) -> Option<Self> {
        if (self.rows, self.cols) != (other.rows, other.cols) {
            return None;
        }
        let one = dom.one();
        let data = self
            .data
            .iter()
            .zip(&other.data)
            .map(|(a, b)| axpy_rows(dom, a, &one, b))
            .collect();
        Some(Self {
            rows: self.rows,
            cols: self.cols,
            data,
        })
    }

    /// `c * self`.
    #[must_use]
    pub fn scale<D: Domain<E = E>>(
        &self,
        dom: &mut D,
        c: &E,
    ) -> Self {
        let data = self
            .data
            .iter()
            .map(|row| {
                let mut out = Vec::with_capacity(row.len());
                for (j, v) in row {
                    let scaled = dom.mul(c, v);
                    if !dom.is_zero(&scaled) {
                        out.push((*j, scaled));
                    }
                }
                out
            })
            .collect();
        Self {
            rows: self.rows,
            cols: self.cols,
            data,
        }
    }

    /// `self * other`, or `None` when the inner dimensions differ. Each
    /// result row is accumulated from the rows of `other` it touches
    /// (Gustavson's algorithm), so the cost follows the number of products
    /// that are actually non-zero.
    #[must_use]
    pub fn matmul<D: Domain<E = E>>(
        &self,
        dom: &mut D,
        other: &Self,
    ) -> Option<Self> {
        if self.cols != other.rows {
            return None;
        }
        let mut data = Vec::with_capacity(self.rows);
        for row in &self.data {
            let mut acc: BTreeMap<usize, E> = BTreeMap::new();
            for (k, a) in row {
                for (j, b) in &other.data[*k] {
                    let term = dom.mul(a, b);
                    let merged = match acc.get(j) {
                        | Some(old) => dom.add(old, &term),
                        | None => term,
                    };
                    acc.insert(*j, merged);
                }
            }
            data.push(
                acc.into_iter()
                    .filter(|(_, v)| !dom.is_zero(v))
                    .collect::<Vec<_>>(),
            );
        }
        Some(Self {
            rows: self.rows,
            cols: other.cols,
            data,
        })
    }

    /// `self * x`, or `None` when the length of `x` is not the column count.
    #[must_use]
    pub fn matvec<D: Domain<E = E>>(
        &self,
        dom: &mut D,
        x: &[E],
    ) -> Option<Vec<E>> {
        if x.len() != self.cols {
            return None;
        }
        Some(
            self.data
                .iter()
                .map(|row| {
                    let mut acc = dom.zero();
                    for (j, v) in row {
                        let term = dom.mul(v, &x[*j]);
                        acc = dom.add(&acc, &term);
                    }
                    acc
                })
                .collect(),
        )
    }
}

/// `a + c * b` for two sorted sparse rows.
fn axpy_rows<D: Domain>(
    dom: &mut D,
    a: &[(usize, D::E)],
    c: &D::E,
    b: &[(usize, D::E)],
) -> Vec<(usize, D::E)> {
    let mut out = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0, 0);
    while i < a.len() || j < b.len() {
        let take_a = j >= b.len() || (i < a.len() && a[i].0 < b[j].0);
        let take_b = i >= a.len() || (j < b.len() && b[j].0 < a[i].0);
        if take_a {
            out.push(a[i].clone());
            i += 1;
        } else if take_b {
            let v = dom.mul(c, &b[j].1);
            if !dom.is_zero(&v) {
                out.push((b[j].0, v));
            }
            j += 1;
        } else {
            let scaled = dom.mul(c, &b[j].1);
            let v = dom.add(&a[i].1, &scaled);
            if !dom.is_zero(&v) {
                out.push((a[i].0, v));
            }
            i += 1;
            j += 1;
        }
    }
    out
}

fn entry<E>(
    row: &[(usize, E)],
    col: usize,
) -> Option<&E> {
    row.binary_search_by(|(j, _)| j.cmp(&col))
        .ok()
        .map(|at| &row[at].1)
}

/// The result of [`echelon`]: a row echelon form up to row and column
/// permutations.
#[derive(Clone, Debug)]
pub struct Echelon<E> {
    /// Number of columns of the eliminated matrix.
    pub cols: usize,
    /// `(row, column)` of each pivot, in elimination order.
    pub pivots: Vec<(usize, usize)>,
    /// The pivot rows after elimination, in the same order. Pivot row `k` is
    /// zero in the pivot columns of every earlier pivot.
    pub rows: Vec<Vec<(usize, E)>>,
    /// The eliminations performed, in order: `(target row, pivot index,
    /// factor)` meaning `row[target] -= factor * pivot_row`. Replaying them
    /// on a right-hand side reproduces the elimination.
    pub ops: Vec<(usize, usize, E)>,
    /// The rows that hold no pivot, as they ended (all zero in the pivotable
    /// columns; possibly non-zero beyond them).
    pub rest: Vec<(usize, Vec<(usize, E)>)>,
}

impl<E: Clone + Debug> Echelon<E> {
    /// The rank over the pivotable columns.
    #[must_use]
    pub fn rank(&self) -> usize {
        self.pivots.len()
    }

    /// The sign of the permutation that brings the pivots to the diagonal
    /// (`true` for odd), for a square matrix of full rank.
    #[must_use]
    pub fn odd_permutation(&self) -> bool {
        let rows: Vec<usize> = self.pivots.iter().map(|p| p.0).collect();
        let cols: Vec<usize> = self.pivots.iter().map(|p| p.1).collect();
        parity(&rows) ^ parity(&cols)
    }
}

/// Whether a sequence of distinct indices is an odd permutation of its
/// sorted order (counted by cycles of the ranking).
fn parity(sequence: &[usize]) -> bool {
    let mut order: Vec<usize> = (0..sequence.len()).collect();
    order.sort_by_key(|&k| sequence[k]);
    // `order[r]` is the position holding the r-th smallest value.
    let mut seen = vec![false; order.len()];
    let mut odd = false;
    for start in 0..order.len() {
        if seen[start] {
            continue;
        }
        let mut len = 0;
        let mut at = start;
        while !seen[at] {
            seen[at] = true;
            at = order[at];
            len += 1;
        }
        if len % 2 == 0 {
            odd = !odd;
        }
    }
    odd
}

/// Gaussian elimination with sparsity-preserving pivoting.
///
/// Pivots are only taken in columns below `pivot_cols`, which lets a caller
/// append right-hand sides as extra columns. `record` keeps the list of row
/// operations (needed to solve for further right-hand sides later).
pub fn echelon<D: Domain>(
    dom: &mut D,
    matrix: &SparseMat<D::E>,
    pivot_cols: usize,
    record: bool,
) -> Echelon<D::E> {
    let pivot_cols = pivot_cols.min(matrix.cols);
    let mut rows: Vec<Vec<(usize, D::E)>> = matrix.data.clone();
    let mut active = vec![true; matrix.rows];
    // Rows that may hold an entry in each column. Entries are added as fill
    // appears and never removed, so a listed row is re-checked before use.
    let mut col_rows: Vec<Vec<usize>> = vec![Vec::new(); matrix.cols];
    let mut col_count = vec![0_usize; matrix.cols];
    for (i, row) in rows.iter().enumerate() {
        for (j, _) in row {
            col_rows[*j].push(i);
            col_count[*j] += 1;
        }
    }
    let width = |row: &[(usize, D::E)]| row.iter().take_while(|(j, _)| *j < pivot_cols).count();
    let mut out = Echelon {
        cols: matrix.cols,
        pivots: Vec::new(),
        rows: Vec::new(),
        ops: Vec::new(),
        rest: Vec::new(),
    };
    loop {
        // The sparsest active row that still has a pivotable entry.
        let mut best: Option<(usize, usize)> = None;
        for (i, row) in rows.iter().enumerate() {
            if !active[i] {
                continue;
            }
            let w = width(row);
            if w > 0 && best.is_none_or(|(_, bw)| w < bw) {
                best = Some((i, w));
                if w == 1 {
                    break;
                }
            }
        }
        let Some((pr, _)) = best else {
            break;
        };
        // In that row, the column with the fewest remaining entries.
        let pc = rows[pr]
            .iter()
            .take_while(|(j, _)| *j < pivot_cols)
            .map(|(j, _)| *j)
            .min_by_key(|&j| (col_count[j], j))
            .unwrap_or(0);
        let Some(pivot) = entry(&rows[pr], pc).cloned() else {
            break;
        };
        active[pr] = false;
        let index = out.pivots.len();
        let pivot_row = rows[pr].clone();
        for (j, _) in &pivot_row {
            col_count[*j] = col_count[*j].saturating_sub(1);
        }
        let targets = std::mem::take(&mut col_rows[pc]);
        for target in targets {
            if !active[target] {
                continue;
            }
            let Some(value) = entry(&rows[target], pc).cloned() else {
                continue;
            };
            let Some(factor) = dom.div(&value, &pivot) else {
                continue;
            };
            let minus = dom.neg(&factor);
            let before: Vec<usize> = rows[target].iter().map(|(j, _)| *j).collect();
            let mut updated = axpy_rows(dom, &rows[target], &minus, &pivot_row);
            // The pivot column cancels exactly by construction; a domain
            // whose normal form misses that must not leave a ghost entry.
            updated.retain(|(j, _)| *j != pc);
            for (j, _) in &updated {
                if before.binary_search(j).is_err() {
                    col_rows[*j].push(target);
                    col_count[*j] += 1;
                }
            }
            for j in &before {
                if updated.binary_search_by(|(c, _)| c.cmp(j)).is_err() {
                    col_count[*j] = col_count[*j].saturating_sub(1);
                }
            }
            rows[target] = updated;
            if record {
                out.ops.push((target, index, factor));
            }
        }
        out.pivots.push((pr, pc));
        out.rows.push(pivot_row);
    }
    for (i, row) in rows.into_iter().enumerate() {
        if active[i] {
            out.rest.push((i, row));
        }
    }
    out
}

/// The rank of a matrix.
pub fn rank<D: Domain>(
    dom: &mut D,
    matrix: &SparseMat<D::E>,
) -> usize {
    echelon(dom, matrix, matrix.cols, false).rank()
}

/// The determinant of a square matrix, or `None` when it is not square.
pub fn det<D: Domain>(
    dom: &mut D,
    matrix: &SparseMat<D::E>,
) -> Option<D::E> {
    if matrix.rows != matrix.cols {
        return None;
    }
    let e = echelon(dom, matrix, matrix.cols, false);
    if e.rank() < matrix.rows {
        return Some(dom.zero());
    }
    let mut product = dom.one();
    for (row, &(_, pc)) in e.rows.iter().zip(&e.pivots) {
        let pivot = entry(row, pc)?;
        product = dom.mul(&product, pivot);
    }
    Some(if e.odd_permutation() {
        dom.neg(&product)
    } else {
        product
    })
}

/// Back-substitution on the pivot rows: the values of the pivot variables
/// given the values `x` of every other column and the right-hand sides
/// `rhs[k]` of the pivot rows.
fn back_substitute<D: Domain>(
    dom: &mut D,
    e: &Echelon<D::E>,
    x: &mut [D::E],
    rhs: &[D::E],
    limit: usize,
) -> Option<()> {
    for k in (0..e.pivots.len()).rev() {
        let pc = e.pivots[k].1;
        let mut acc = rhs[k].clone();
        let mut pivot = None;
        for (j, v) in &e.rows[k] {
            if *j == pc {
                pivot = Some(v.clone());
            } else if *j < limit {
                let term = dom.mul(v, &x[*j]);
                acc = dom.sub(&acc, &term);
            }
        }
        x[pc] = dom.div(&acc, &pivot?)?;
    }
    Some(())
}

/// What [`solve`] found.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Solution<E> {
    /// The system has no solution.
    Inconsistent,
    /// A particular solution and a basis of the null space; the general
    /// solution is the particular one plus any combination of the basis.
    Family {
        /// One solution (free variables set to zero).
        particular: Vec<E>,
        /// A basis of the null space, empty when the solution is unique.
        nullspace: Vec<Vec<E>>,
    },
}

/// Solves `A x = b`, or returns `None` when `b` has the wrong length.
pub fn solve<D: Domain>(
    dom: &mut D,
    matrix: &SparseMat<D::E>,
    b: &[D::E],
) -> Option<Solution<D::E>> {
    if b.len() != matrix.rows {
        return None;
    }
    let n = matrix.cols;
    let mut augmented = matrix.clone();
    augmented.cols = n + 1;
    for (row, value) in augmented.data.iter_mut().zip(b) {
        if !dom.is_zero(value) {
            row.push((n, value.clone()));
        }
    }
    let e = echelon(dom, &augmented, n, false);
    if e.rest.iter().any(|(_, row)| !row.is_empty()) {
        return Some(Solution::Inconsistent);
    }
    let zero = dom.zero();
    let rhs: Vec<D::E> = e
        .rows
        .iter()
        .map(|row| entry(row, n).cloned().unwrap_or_else(|| zero.clone()))
        .collect();
    let mut particular = vec![zero.clone(); n];
    back_substitute(dom, &e, &mut particular, &rhs, n)?;
    let nullspace = null_basis(dom, &e, n)?;
    Some(Solution::Family {
        particular,
        nullspace,
    })
}

/// A basis of the null space from an echelon form over `n` columns.
fn null_basis<D: Domain>(
    dom: &mut D,
    e: &Echelon<D::E>,
    n: usize,
) -> Option<Vec<Vec<D::E>>> {
    let mut is_pivot = vec![false; n];
    for &(_, pc) in &e.pivots {
        if pc < n {
            is_pivot[pc] = true;
        }
    }
    let zero = dom.zero();
    let zeros = vec![zero.clone(); e.pivots.len()];
    let mut basis = Vec::new();
    for free in (0..n).filter(|&j| !is_pivot[j]) {
        let mut x = vec![zero.clone(); n];
        x[free] = dom.one();
        back_substitute(dom, e, &mut x, &zeros, n)?;
        basis.push(x);
    }
    Some(basis)
}

/// A basis of the null space `{x : A x = 0}`.
pub fn nullspace<D: Domain>(
    dom: &mut D,
    matrix: &SparseMat<D::E>,
) -> Vec<Vec<D::E>> {
    let e = echelon(dom, matrix, matrix.cols, false);
    null_basis(dom, &e, matrix.cols).unwrap_or_default()
}

// ---------------------------------------------------------------------------
// Modular methods for rational matrices
// ---------------------------------------------------------------------------

/// Deterministic Miller-Rabin for 64-bit integers.
#[must_use]
pub fn is_prime_u64(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for p in [2_u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n % p == 0 {
            return n == p;
        }
    }
    let field = Gfp(n);
    let (mut d, mut s) = (n - 1, 0_u32);
    while d % 2 == 0 {
        d /= 2;
        s += 1;
    }
    'witness: for a in [2_u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = field.pow_mod(a, d);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 1..s {
            x = field.mul_mod(x, x);
            if x == n - 1 {
                continue 'witness;
            }
        }
        return false;
    }
    true
}

/// The primes below `2^62`, largest first: the moduli of the modular
/// methods. Large enough that an unlucky prime (one dividing a determinant)
/// is practically never met, small enough for `u128` products.
pub fn word_primes() -> impl Iterator<Item = u64> {
    let mut n = 1_u64 << 62;
    std::iter::from_fn(move || {
        loop {
            n -= 1;
            if is_prime_u64(n) {
                return Some(n);
            }
        }
    })
}

/// A rational matrix scaled row by row to integer entries: the integer
/// matrix and the factor each row was multiplied by.
fn integer_rows(matrix: &SparseMat<BigRational>) -> (SparseMat<BigInt>, Vec<BigInt>) {
    let mut scales = Vec::with_capacity(matrix.rows);
    let mut data = Vec::with_capacity(matrix.rows);
    for row in &matrix.data {
        let scale = row
            .iter()
            .fold(BigInt::one(), |acc, (_, v)| acc.lcm(v.denom()));
        data.push(
            row.iter()
                .map(|(j, v)| (*j, (v * BigRational::from(scale.clone())).to_integer()))
                .collect(),
        );
        scales.push(scale);
    }
    (
        SparseMat {
            rows: matrix.rows,
            cols: matrix.cols,
            data,
        },
        scales,
    )
}

fn reduce_matrix(
    matrix: &SparseMat<BigInt>,
    field: Gfp,
) -> SparseMat<u64> {
    SparseMat {
        rows: matrix.rows,
        cols: matrix.cols,
        data: matrix
            .data
            .iter()
            .map(|row| {
                row.iter()
                    .map(|(j, v)| (*j, field.reduce(v)))
                    .filter(|(_, v)| *v != 0)
                    .collect()
            })
            .collect(),
    }
}

/// An upper bound of `sqrt(sum of squares)` of a vector of integers, at
/// least one.
fn norm_bound<'a>(values: impl Iterator<Item = &'a BigInt>) -> BigInt {
    let squares: BigInt = values.map(|v| v * v).sum();
    (squares.sqrt() + BigInt::one()).max(BigInt::one())
}

/// The determinant of a square rational matrix, exactly, by modular
/// determinants combined with the Chinese remainder theorem up to
/// Hadamard's bound. `None` when the matrix is not square.
#[must_use]
pub fn det_rational(matrix: &SparseMat<BigRational>) -> Option<BigRational> {
    if matrix.rows != matrix.cols {
        return None;
    }
    if matrix.rows == 0 {
        return Some(BigRational::one());
    }
    let (integer, scales) = integer_rows(matrix);
    if integer.data.iter().any(Vec::is_empty) {
        return Some(BigRational::zero());
    }
    // |det| <= product of the row norms.
    let bound: BigInt = integer
        .data
        .iter()
        .map(|row| norm_bound(row.iter().map(|(_, v)| v)))
        .product();
    let target = bound * BigInt::from(2) + BigInt::one();
    let mut modulus = BigInt::one();
    let mut residue = BigInt::zero();
    for p in word_primes() {
        let mut field = Gfp(p);
        let reduced = reduce_matrix(&integer, field);
        let d = det(&mut field, &reduced)?;
        // Chinese remainder: residue + modulus * t = d (mod p).
        let have = field.reduce(&residue);
        let diff = field.sub(&d, &have);
        let inv = field.inv(field.reduce(&modulus))?;
        let t = field.mul_mod(diff, inv);
        residue += &modulus * BigInt::from(t);
        modulus *= BigInt::from(p);
        if modulus > target {
            break;
        }
    }
    // The symmetric representative.
    if &residue * BigInt::from(2) > modulus {
        residue -= &modulus;
    }
    let scale: BigInt = scales.iter().product();
    Some(BigRational::new(residue, scale))
}

/// Rational reconstruction: the fraction `n / d` with `n = a d (mod m)` and
/// `|n|, d <= sqrt(m / 2)`, if there is one.
#[must_use]
pub fn rational_reconstruct(
    a: &BigInt,
    m: &BigInt,
) -> Option<BigRational> {
    let bound = (m / BigInt::from(2)).sqrt();
    let (mut r0, mut r1) = (m.clone(), a.mod_floor(m));
    let (mut t0, mut t1) = (BigInt::zero(), BigInt::one());
    while r1 > bound {
        let q = &r0 / &r1;
        let r2 = &r0 - &q * &r1;
        let t2 = &t0 - &q * &t1;
        r0 = r1;
        r1 = r2;
        t0 = t1;
        t1 = t2;
    }
    if t1.is_zero() || t1.abs() > bound || !r1.gcd(&t1).is_one() {
        return None;
    }
    Some(BigRational::new(r1, t1))
}

/// A sparse factorisation modulo a prime, kept for repeated solves.
struct ModularFactor {
    field: Gfp,
    echelon: Echelon<u64>,
    n: usize,
}

impl ModularFactor {
    fn new(
        matrix: &SparseMat<u64>,
        field: Gfp,
    ) -> Option<Self> {
        let mut dom = field;
        let echelon = echelon(&mut dom, matrix, matrix.cols, true);
        (echelon.rank() == matrix.rows).then_some(Self {
            field,
            echelon,
            n: matrix.cols,
        })
    }

    /// `A^-1 b mod p`.
    fn solve(
        &self,
        b: &[u64],
    ) -> Option<Vec<u64>> {
        let mut dom = self.field;
        let mut r = b.to_vec();
        for &(target, index, factor) in &self.echelon.ops {
            let source = r[self.echelon.pivots[index].0];
            let scaled = dom.mul(&factor, &source);
            r[target] = dom.sub(&r[target], &scaled);
        }
        let rhs: Vec<u64> = self.echelon.pivots.iter().map(|&(pr, _)| r[pr]).collect();
        let mut x = vec![0_u64; self.n];
        back_substitute(&mut dom, &self.echelon, &mut x, &rhs, self.n)?;
        Some(x)
    }
}

/// Solves a square non-singular rational system `A x = b` exactly by
/// Dixon's `p`-adic lifting.
///
/// The matrix is factored once modulo a word-sized prime; every further
/// `p`-adic digit of the solution costs one sparse solve modulo `p` and one
/// sparse integer product, so no rational number ever appears during the
/// elimination. The lifted residue is turned back into fractions by rational
/// reconstruction and the result is checked by substitution, which makes the
/// answer exact rather than "correct up to a bound".
///
/// Returns `None` when the matrix is not square, `b` has the wrong length, or
/// the matrix is singular modulo several primes in a row (which, for these
/// primes, means singular); use [`solve`] to analyse a singular system.
#[must_use]
pub fn solve_dixon(
    matrix: &SparseMat<BigRational>,
    b: &[BigRational],
) -> Option<Vec<BigRational>> {
    let n = matrix.rows;
    if n != matrix.cols || b.len() != n {
        return None;
    }
    if n == 0 {
        return Some(Vec::new());
    }
    let (integer, scales) = integer_rows(matrix);
    // Scale each right-hand side with its row, then clear what is left.
    let scaled: Vec<BigRational> = b
        .iter()
        .zip(&scales)
        .map(|(v, s)| v * BigRational::from(s.clone()))
        .collect();
    let common = scaled
        .iter()
        .fold(BigInt::one(), |acc, v| acc.lcm(v.denom()));
    let rhs: Vec<BigInt> = scaled
        .iter()
        .map(|v| (v * BigRational::from(common.clone())).to_integer())
        .collect();
    let (field, factor) = word_primes().take(4).find_map(|p| {
        let field = Gfp(p);
        ModularFactor::new(&reduce_matrix(&integer, field), field).map(|f| (field, f))
    })?;
    let p = BigInt::from(field.0);
    // Cramer's rule with Hadamard's inequality bounds numerators and
    // denominators of the solution of the integer system.
    let transposed = integer.transpose();
    let hadamard: BigInt = transposed
        .data
        .iter()
        .map(|column| norm_bound(column.iter().map(|(_, v)| v)))
        .product();
    let limit = &hadamard * &hadamard * norm_bound(rhs.iter()) * BigInt::from(4);
    let mut ints = BigIntegers;
    let mut residual = rhs.clone();
    let mut lifted = vec![BigInt::zero(); n];
    let mut power = BigInt::one();
    let mut next_try = 4_usize;
    let mut step = 0_usize;
    loop {
        let digits = factor.solve(&residual.iter().map(|v| field.reduce(v)).collect::<Vec<_>>())?;
        let digit_ints: Vec<BigInt> = digits.iter().map(|&d| BigInt::from(d)).collect();
        let product = integer.matvec(&mut ints, &digit_ints)?;
        for ((r, ax), (x, d)) in residual
            .iter_mut()
            .zip(&product)
            .zip(lifted.iter_mut().zip(&digit_ints))
        {
            // Exact: r = A d (mod p) by construction.
            *r = (&*r - ax) / &p;
            *x += d * &power;
        }
        power *= &p;
        step += 1;
        let done = power > limit;
        // Small solutions are found long before the bound; try early, at
        // geometrically spaced steps, and always at the bound.
        if step == next_try || done || residual.iter().all(Zero::is_zero) {
            next_try = next_try.saturating_mul(2);
            let candidate: Option<Vec<BigRational>> = lifted
                .iter()
                .map(|v| rational_reconstruct(v, &power))
                .collect();
            if let Some(x) = candidate
                && verifies(&integer, &rhs, &x)
            {
                let unscale = BigRational::from(common.clone());
                return Some(x.into_iter().map(|v| v / &unscale).collect());
            }
            if done {
                return None;
            }
        }
    }
}

/// Whether `A x = b` holds exactly for an integer matrix and vector.
fn verifies(
    matrix: &SparseMat<BigInt>,
    b: &[BigInt],
    x: &[BigRational],
) -> bool {
    matrix.data.iter().zip(b).all(|(row, want)| {
        let sum: BigRational = row
            .iter()
            .map(|(j, v)| BigRational::from(v.clone()) * &x[*j])
            .sum();
        sum == BigRational::from(want.clone())
    })
}

/// The integers as a [`Domain`] for products (division is exact division
/// and fails otherwise; elimination does not use this domain).
struct BigIntegers;

impl Domain for BigIntegers {
    type E = BigInt;

    fn zero(&mut self) -> BigInt {
        BigInt::zero()
    }

    fn one(&mut self) -> BigInt {
        BigInt::one()
    }

    fn is_zero(
        &mut self,
        a: &BigInt,
    ) -> bool {
        a.is_zero()
    }

    fn add(
        &mut self,
        a: &BigInt,
        b: &BigInt,
    ) -> BigInt {
        a + b
    }

    fn sub(
        &mut self,
        a: &BigInt,
        b: &BigInt,
    ) -> BigInt {
        a - b
    }

    fn mul(
        &mut self,
        a: &BigInt,
        b: &BigInt,
    ) -> BigInt {
        a * b
    }

    fn div(
        &mut self,
        a: &BigInt,
        b: &BigInt,
    ) -> Option<BigInt> {
        (!b.is_zero() && (a % b).is_zero()).then(|| a / b)
    }
}

// ---------------------------------------------------------------------------
// Black-box methods over a prime field
// ---------------------------------------------------------------------------

/// Berlekamp-Massey: the shortest linear recurrence of a sequence over
/// `GF(p)`, as the monic connection polynomial `c` (lowest degree first)
/// with `sum_k c[k] * s[i + k] = 0` for every `i`.
#[must_use]
pub fn berlekamp_massey(
    sequence: &[u64],
    field: Gfp,
) -> Vec<u64> {
    let mut dom = field;
    // c(x) = 1 + c1 x + ...: s[i] + c1 s[i-1] + ... = 0.
    let mut c: Vec<u64> = vec![1];
    let mut prev: Vec<u64> = vec![1];
    let mut len = 0_usize;
    let mut shift = 1_usize;
    let mut last = 1_u64;
    for i in 0..sequence.len() {
        let mut d = sequence[i] % field.0;
        for k in 1..=len.min(c.len() - 1) {
            let term = dom.mul(&c[k], &sequence[i - k]);
            d = dom.add(&d, &term);
        }
        if d == 0 {
            shift += 1;
            continue;
        }
        let Some(coef) = dom.div(&d, &last) else {
            break;
        };
        let old = c.clone();
        if c.len() < prev.len() + shift {
            c.resize(prev.len() + shift, 0);
        }
        for (k, v) in prev.iter().enumerate() {
            let term = dom.mul(&coef, v);
            c[k + shift] = dom.sub(&c[k + shift], &term);
        }
        if 2 * len <= i {
            len = i + 1 - len;
            prev = old;
            last = d;
            shift = 1;
        } else {
            shift += 1;
        }
    }
    c.resize(len + 1, 0);
    // Reverse into the form sum_k m[k] s[i + k] = 0, monic at the top.
    c.reverse();
    c
}

/// A small deterministic generator for the random projections.
struct Lcg(u64);

impl Lcg {
    fn next(
        &mut self,
        field: Gfp,
    ) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0 % field.0
    }
}

/// The minimal polynomial of the sequence `u^T A^i v`, `i = 0, 1, ...`.
fn projected_minpoly(
    matrix: &SparseMat<u64>,
    field: Gfp,
    u: &[u64],
    v: &[u64],
) -> Option<Vec<u64>> {
    let mut dom = field;
    let n = matrix.rows;
    let mut w = v.to_vec();
    let mut sequence = Vec::with_capacity(2 * n + 2);
    for _ in 0..2 * n + 2 {
        let mut dot = 0_u64;
        for (a, b) in u.iter().zip(&w) {
            let term = dom.mul(a, b);
            dot = dom.add(&dot, &term);
        }
        sequence.push(dot);
        w = matrix.matvec(&mut dom, &w)?;
    }
    Some(berlekamp_massey(&sequence, field))
}

/// The minimal polynomial of a square matrix over `GF(p)` by Wiedemann's
/// algorithm, monic, lowest degree first.
///
/// The algorithm multiplies by the matrix `2n` times per trial and stores
/// two vectors: no entry of the matrix is ever changed, so there is no
/// fill-in. One random projection yields a divisor of the minimal
/// polynomial, and the right one with probability at least `1 - 2n/p`; the
/// result here is the least common multiple over `trials` projections
/// (seeded, so the function is deterministic). `None` when the matrix is not
/// square.
#[must_use]
pub fn minpoly_wiedemann(
    matrix: &SparseMat<u64>,
    field: Gfp,
    trials: usize,
    seed: u64,
) -> Option<Vec<u64>> {
    let n = matrix.rows;
    if n != matrix.cols {
        return None;
    }
    let mut rng = Lcg(seed | 1);
    let mut best: Vec<u64> = vec![1];
    for _ in 0..trials.max(1) {
        let u: Vec<u64> = (0..n).map(|_| rng.next(field)).collect();
        let v: Vec<u64> = (0..n).map(|_| rng.next(field)).collect();
        let m = projected_minpoly(matrix, field, &u, &v)?;
        best = poly_lcm(&best, &m, field);
    }
    Some(best)
}

/// Solves a square non-singular system `A x = b` over `GF(p)` by Wiedemann's
/// algorithm, or returns `None` when the projection does not determine a
/// solution (singular matrix, or an unlucky projection in each of the
/// `trials`). A returned vector always satisfies `A x = b`.
#[must_use]
pub fn solve_wiedemann(
    matrix: &SparseMat<u64>,
    b: &[u64],
    field: Gfp,
    trials: usize,
    seed: u64,
) -> Option<Vec<u64>> {
    let n = matrix.rows;
    if n != matrix.cols || b.len() != n {
        return None;
    }
    let mut dom = field;
    let mut rng = Lcg(seed | 1);
    for _ in 0..trials.max(1) {
        let u: Vec<u64> = (0..n).map(|_| rng.next(field)).collect();
        let m = projected_minpoly(matrix, field, &u, b)?;
        // m(A) b = 0 with m(0) != 0 gives x = -(1/m0) sum_{k>=1} m_k A^(k-1) b.
        let Some(m0_inv) = m.first().and_then(|&c| field.inv(c)) else {
            continue;
        };
        let mut x = vec![0_u64; n];
        let mut w = b.to_vec();
        for &coefficient in m.iter().skip(1) {
            for (xi, wi) in x.iter_mut().zip(&w) {
                let term = dom.mul(&coefficient, wi);
                *xi = dom.add(xi, &term);
            }
            w = matrix.matvec(&mut dom, &w)?;
        }
        let scale = dom.neg(&m0_inv);
        for xi in &mut x {
            *xi = dom.mul(xi, &scale);
        }
        if matrix.matvec(&mut dom, &x)? == b.iter().map(|v| v % field.0).collect::<Vec<_>>() {
            return Some(x);
        }
    }
    None
}

fn poly_trim(p: &mut Vec<u64>) {
    while p.len() > 1 && p.last() == Some(&0) {
        p.pop();
    }
}

/// Remainder of `a` by `b` over `GF(p)` (lowest degree first).
fn poly_rem(
    a: &[u64],
    b: &[u64],
    field: Gfp,
) -> Vec<u64> {
    let mut dom = field;
    let mut r = a.to_vec();
    poly_trim(&mut r);
    let Some(lead_inv) = b.last().and_then(|&c| field.inv(c)) else {
        return r;
    };
    while r.len() >= b.len() && !(r.len() == 1 && r[0] == 0) {
        let top = *r.last().unwrap_or(&0);
        let factor = dom.mul(&top, &lead_inv);
        let offset = r.len() - b.len();
        for (k, c) in b.iter().enumerate() {
            let term = dom.mul(&factor, c);
            r[offset + k] = dom.sub(&r[offset + k], &term);
        }
        r.pop();
        if r.is_empty() {
            r.push(0);
        }
        poly_trim(&mut r);
    }
    r
}

fn poly_gcd(
    a: &[u64],
    b: &[u64],
    field: Gfp,
) -> Vec<u64> {
    let (mut x, mut y) = (a.to_vec(), b.to_vec());
    poly_trim(&mut x);
    poly_trim(&mut y);
    while !(y.len() == 1 && y[0] == 0) {
        let r = poly_rem(&x, &y, field);
        x = y;
        y = r;
    }
    // Monic.
    if let Some(inv) = x.last().and_then(|&c| field.inv(c)) {
        for c in &mut x {
            *c = field.mul_mod(*c, inv);
        }
    }
    x
}

/// Exact quotient `a / b` over `GF(p)`.
fn poly_div_exact(
    a: &[u64],
    b: &[u64],
    field: Gfp,
) -> Vec<u64> {
    let mut dom = field;
    let mut r = a.to_vec();
    let Some(lead_inv) = b.last().and_then(|&c| field.inv(c)) else {
        return vec![0];
    };
    if r.len() < b.len() {
        return vec![0];
    }
    let mut q = vec![0_u64; r.len() - b.len() + 1];
    for offset in (0..q.len()).rev() {
        let top = r[offset + b.len() - 1];
        let factor = dom.mul(&top, &lead_inv);
        q[offset] = factor;
        for (k, c) in b.iter().enumerate() {
            let term = dom.mul(&factor, c);
            r[offset + k] = dom.sub(&r[offset + k], &term);
        }
    }
    q
}

fn poly_mul(
    a: &[u64],
    b: &[u64],
    field: Gfp,
) -> Vec<u64> {
    let mut dom = field;
    let mut out = vec![0_u64; a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            let term = dom.mul(x, y);
            out[i + j] = dom.add(&out[i + j], &term);
        }
    }
    out
}

fn poly_lcm(
    a: &[u64],
    b: &[u64],
    field: Gfp,
) -> Vec<u64> {
    let g = poly_gcd(a, b, field);
    let mut l = poly_div_exact(&poly_mul(a, b, field), &g, field);
    poly_trim(&mut l);
    if let Some(inv) = l.last().and_then(|&c| field.inv(c)) {
        for c in &mut l {
            *c = field.mul_mod(*c, inv);
        }
    }
    l
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

    use super::*;
    use crate::kernels::qlinalg;

    fn q(n: i64) -> BigRational {
        BigRational::from(BigInt::from(n))
    }

    fn frac(
        n: i64,
        d: i64,
    ) -> BigRational {
        BigRational::new(BigInt::from(n), BigInt::from(d))
    }

    /// A deterministic sparse rational matrix with about `per_row` entries
    /// per row and a dominant diagonal (so that it is non-singular).
    fn sample(
        n: usize,
        per_row: usize,
        seed: u64,
    ) -> SparseMat<BigRational> {
        let mut s = seed | 1;
        let mut next = move || {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            s
        };
        let mut triplets = Vec::new();
        for i in 0..n {
            triplets.push((i, i, q(20 + (next() % 7) as i64)));
            for _ in 0..per_row {
                let j = (next() % n as u64) as usize;
                let num = (next() % 9) as i64 - 4;
                let den = 1 + (next() % 3) as i64;
                triplets.push((i, j, frac(num, den)));
            }
        }
        SparseMat::from_triplets(&mut Rationals, n, n, triplets).unwrap()
    }

    #[test]
    fn construction_merges_drops_zeros_and_checks_bounds() {
        let m = SparseMat::from_triplets(
            &mut Rationals,
            2,
            3,
            vec![(0, 2, q(1)), (0, 0, q(2)), (0, 2, q(-1)), (1, 1, q(5))],
        )
        .unwrap();
        assert_eq!(m.data[0], vec![(0, q(2))]);
        assert_eq!(m.nnz(), 2);
        assert!(SparseMat::from_triplets(&mut Rationals, 2, 2, vec![(2, 0, q(1))]).is_err());
        let dense = m.to_dense(&mut Rationals);
        assert_eq!(SparseMat::from_dense(&mut Rationals, &dense), m);
        assert_eq!(m.transpose().transpose(), m);
    }

    #[test]
    fn products_and_sums_agree_with_the_dense_routines() {
        let a = sample(12, 3, 1);
        let b = sample(12, 2, 2);
        let mut r = Rationals;
        let product = a.matmul(&mut r, &b).unwrap().to_dense(&mut r);
        assert_eq!(
            product,
            qlinalg::matmul(&a.to_dense(&mut r), &b.to_dense(&mut r))
        );
        let sum = a.add(&mut r, &b).unwrap().to_dense(&mut r);
        for (i, row) in sum.iter().enumerate() {
            for (j, v) in row.iter().enumerate() {
                assert_eq!(*v, &a.to_dense(&mut r)[i][j] + &b.to_dense(&mut r)[i][j]);
            }
        }
        assert!(a.matmul(&mut r, &SparseMat::zeros(5, 5)).is_none());
        assert_eq!(a.scale(&mut r, &q(0)).nnz(), 0);
    }

    #[test]
    fn determinant_three_ways() {
        for (n, seed) in [(1, 3), (2, 4), (9, 5), (25, 6)] {
            let a = sample(n, 3, seed);
            let dense = qlinalg::det(&a.to_dense(&mut Rationals));
            assert_eq!(
                det(&mut Rationals, &a).unwrap(),
                dense,
                "elimination, n = {n}"
            );
            assert_eq!(det_rational(&a).unwrap(), dense, "modular, n = {n}");
        }
        // Singular, and a permutation matrix with an odd permutation.
        let singular = SparseMat::from_triplets(
            &mut Rationals,
            3,
            3,
            vec![(0, 0, q(1)), (1, 0, q(2)), (2, 2, q(3))],
        )
        .unwrap();
        assert_eq!(det(&mut Rationals, &singular).unwrap(), q(0));
        assert_eq!(det_rational(&singular).unwrap(), q(0));
        let swap = SparseMat::from_triplets(
            &mut Rationals,
            3,
            3,
            vec![(0, 1, q(1)), (1, 0, q(1)), (2, 2, q(1))],
        )
        .unwrap();
        assert_eq!(det(&mut Rationals, &swap).unwrap(), q(-1));
        assert_eq!(det_rational(&swap).unwrap(), q(-1));
        assert!(det_rational(&SparseMat::zeros(2, 3)).is_none());
        assert_eq!(det_rational(&SparseMat::zeros(0, 0)).unwrap(), q(1));
    }

    #[test]
    fn dixon_lifting_solves_exactly() {
        for (n, seed) in [(1, 7), (6, 8), (30, 9)] {
            let a = sample(n, 3, seed);
            let b: Vec<BigRational> = (0..n)
                .map(|i| frac(i as i64 * 3 - 5, 1 + (i % 4) as i64))
                .collect();
            let x = solve_dixon(&a, &b).expect("non-singular");
            assert_eq!(a.matvec(&mut Rationals, &x).unwrap(), b, "n = {n}");
            // The generic eliminator finds the same unique solution.
            match solve(&mut Rationals, &a, &b).unwrap() {
                | Solution::Family {
                    particular,
                    nullspace,
                } => {
                    assert_eq!(particular, x);
                    assert!(nullspace.is_empty());
                },
                | Solution::Inconsistent => panic!("a solvable system"),
            }
        }
        // A solution with large numerators and denominators: the Hilbert matrix.
        let n = 8;
        let hilbert = SparseMat::from_dense(
            &mut Rationals,
            &(0..n)
                .map(|i| (0..n).map(|j| frac(1, (i + j + 1) as i64)).collect())
                .collect::<Vec<Vec<BigRational>>>(),
        );
        let ones = vec![q(1); n];
        let x = solve_dixon(&hilbert, &ones).unwrap();
        assert_eq!(hilbert.matvec(&mut Rationals, &x).unwrap(), ones);
        assert_eq!(x[0], q(-8));
        // Singular and misshapen systems are refused.
        let singular =
            SparseMat::from_triplets(&mut Rationals, 2, 2, vec![(0, 0, q(1)), (1, 0, q(1))])
                .unwrap();
        assert!(solve_dixon(&singular, &[q(1), q(1)]).is_none());
        assert!(solve_dixon(&sample(3, 1, 1), &[q(1)]).is_none());
    }

    #[test]
    fn underdetermined_and_inconsistent_systems() {
        // x + y + z = 6, y - z = 0: one free variable.
        let a = SparseMat::from_triplets(
            &mut Rationals,
            2,
            3,
            vec![
                (0, 0, q(1)),
                (0, 1, q(1)),
                (0, 2, q(1)),
                (1, 1, q(1)),
                (1, 2, q(-1)),
            ],
        )
        .unwrap();
        let b = vec![q(6), q(0)];
        let Solution::Family {
            particular,
            nullspace,
        } = solve(&mut Rationals, &a, &b).unwrap()
        else {
            panic!("consistent");
        };
        assert_eq!(a.matvec(&mut Rationals, &particular).unwrap(), b);
        assert_eq!(nullspace.len(), 1);
        assert_eq!(
            a.matvec(&mut Rationals, &nullspace[0]).unwrap(),
            vec![q(0), q(0)]
        );
        assert_eq!(rank(&mut Rationals, &a), 2);
        assert_eq!(super::nullspace(&mut Rationals, &a).len(), 1);
        // x = 1 and x = 2.
        let clash =
            SparseMat::from_triplets(&mut Rationals, 2, 1, vec![(0, 0, q(1)), (1, 0, q(1))])
                .unwrap();
        assert_eq!(
            solve(&mut Rationals, &clash, &[q(1), q(2)]).unwrap(),
            Solution::Inconsistent
        );
        assert!(solve(&mut Rationals, &clash, &[q(1)]).is_none());
    }

    #[test]
    fn pivoting_avoids_the_fill_of_an_arrow_matrix() {
        // Dense first row and column, diagonal otherwise. Eliminating in
        // natural order fills the whole matrix; a sparsity-aware order
        // creates no fill at all.
        let n = 60;
        let mut triplets = Vec::new();
        for i in 0..n {
            triplets.push((i, i, q(i as i64 + 2)));
            if i > 0 {
                triplets.push((0, i, q(1)));
                triplets.push((i, 0, q(1)));
            }
        }
        let arrow = SparseMat::from_triplets(&mut Rationals, n, n, triplets).unwrap();
        let e = echelon(&mut Rationals, &arrow, n, false);
        assert_eq!(e.rank(), n);
        let stored: usize = e.rows.iter().map(Vec::len).sum();
        assert!(stored <= arrow.nnz(), "{stored} entries after elimination");
        assert_eq!(
            det(&mut Rationals, &arrow).unwrap(),
            qlinalg::det(&arrow.to_dense(&mut Rationals))
        );
    }

    #[test]
    fn prime_field_arithmetic_and_primes() {
        assert!(is_prime_u64(2) && is_prime_u64(1_000_000_007) && !is_prime_u64(1));
        assert!(!is_prime_u64(3_215_031_751));
        let p = word_primes().next().unwrap();
        assert!(p < (1 << 62) && is_prime_u64(p));
        let field = Gfp(101);
        let mut dom = field;
        assert_eq!(dom.sub(&3, &5), 99);
        assert_eq!(dom.div(&1, &2), Some(51));
        assert_eq!(dom.div(&1, &0), None);
        assert_eq!(field.reduce(&BigInt::from(-1)), 100);
    }

    #[test]
    fn rational_reconstruction_recovers_small_fractions() {
        let p = BigInt::from(1_000_000_007_u64);
        let inv7 = BigInt::from(7).modpow(&(&p - BigInt::from(2)), &p);
        assert_eq!(
            rational_reconstruct(&(BigInt::from(3) * inv7).mod_floor(&p), &p),
            Some(frac(3, 7))
        );
        assert_eq!(rational_reconstruct(&BigInt::from(5), &p), Some(q(5)));
        let inv9 = BigInt::from(9).modpow(&(&p - BigInt::from(2)), &p);
        assert_eq!(
            rational_reconstruct(&(BigInt::from(-22) * inv9).mod_floor(&p), &p),
            Some(frac(-22, 9))
        );
        // A fraction too large for the modulus is not invented.
        let inv = BigInt::from(99_991).modpow(&(&p - BigInt::from(2)), &p);
        assert_ne!(
            rational_reconstruct(&(BigInt::from(99_989) * inv).mod_floor(&p), &p),
            Some(frac(99_989, 99_991))
        );
    }

    #[test]
    fn berlekamp_massey_finds_fibonacci() {
        let field = Gfp(1_000_000_007);
        let mut s = vec![0_u64, 1];
        for i in 2..20 {
            s.push((s[i - 1] + s[i - 2]) % field.0);
        }
        // s[i+2] - s[i+1] - s[i] = 0.
        assert_eq!(
            berlekamp_massey(&s, field),
            vec![field.0 - 1, field.0 - 1, 1]
        );
        assert_eq!(berlekamp_massey(&[0, 0, 0], field), vec![1]);
    }

    #[test]
    fn wiedemann_agrees_with_elimination() {
        let field = Gfp(word_primes().next().unwrap());
        let rational = sample(40, 3, 11);
        let (integer, _) = integer_rows(&rational);
        let a = reduce_matrix(&integer, field);
        let b: Vec<u64> = (0..40).map(|i| (i * i + 1) as u64).collect();
        let x = solve_wiedemann(&a, &b, field, 3, 42).expect("non-singular");
        let mut dom = field;
        assert_eq!(a.matvec(&mut dom, &x).unwrap(), b);
        // The same solution as the factorisation gives (it is unique).
        let direct = ModularFactor::new(&a, field).unwrap().solve(&b).unwrap();
        assert_eq!(x, direct);
        // The polynomial found annihilates the matrix: m(A) e = 0 for a
        // vector with every coordinate set.
        let m = minpoly_wiedemann(&a, field, 2, 7).unwrap();
        assert!(m.len() <= 41 && m.last() == Some(&1));
        let mut power = vec![1_u64; 40];
        let mut total = vec![0_u64; 40];
        for &c in &m {
            for (t, v) in total.iter_mut().zip(&power) {
                let term = dom.mul(&c, v);
                *t = dom.add(t, &term);
            }
            power = a.matvec(&mut dom, &power).unwrap();
        }
        assert!(total.iter().all(|&v| v == 0));
        // A tridiagonal matrix with a non-zero subdiagonal has a cyclic
        // vector, so its minimal polynomial is its characteristic polynomial:
        // full degree, and the constant term is (-1)^n det.
        let n = 12;
        let mut triplets = Vec::new();
        for i in 0..n {
            triplets.push((i, i, (i + 2) as u64));
            if i + 1 < n {
                triplets.push((i, i + 1, 1));
                triplets.push((i + 1, i, 1));
            }
        }
        let jacobi = SparseMat::from_triplets(&mut dom, n, n, triplets).unwrap();
        let m = minpoly_wiedemann(&jacobi, field, 3, 7).unwrap();
        assert_eq!(m.len(), n + 1);
        assert_eq!(m[0], det(&mut dom, &jacobi).unwrap(), "(-1)^12 det");
        // A projection (A^2 = A) has minimal polynomial x^2 - x.
        let projection = SparseMat::from_triplets(
            &mut dom,
            4,
            4,
            vec![(0, 0, 1), (1, 1, 1), (2, 0, 5), (3, 1, 9)],
        )
        .unwrap();
        assert_eq!(
            minpoly_wiedemann(&projection, field, 4, 3).unwrap(),
            vec![0, field.0 - 1, 1]
        );
        assert!(solve_wiedemann(&projection, &[1, 2, 3, 4], field, 3, 5).is_none());
        assert!(minpoly_wiedemann(&SparseMat::zeros(2, 3), field, 1, 1).is_none());
    }

    #[test]
    fn empty_and_degenerate_matrices() {
        let empty: SparseMat<BigRational> = SparseMat::zeros(0, 0);
        assert_eq!(rank(&mut Rationals, &empty), 0);
        assert_eq!(det(&mut Rationals, &empty).unwrap(), q(1));
        assert_eq!(solve_dixon(&empty, &[]).unwrap(), Vec::<BigRational>::new());
        let zero: SparseMat<BigRational> = SparseMat::zeros(3, 4);
        assert_eq!(rank(&mut Rationals, &zero), 0);
        assert_eq!(nullspace(&mut Rationals, &zero).len(), 4);
        assert!(det(&mut Rationals, &zero).is_none());
    }
}
