//! Exact sparse linear algebra.
//!
//! Sparsity-preserving elimination over any domain, modular determinants,
//! Dixon lifting and Wiedemann's algorithm.
//!
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

use std::alloc::Layout;
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

/// The field of integers modulo a prime, with the reduction context for that
/// prime attached.
///
/// Elements are reduced residues in `[0, p)`. The context carries a Barrett
/// reducer, so a multiplication inside an elimination is a multiply, a shift
/// and at most two subtractions instead of a 128-bit hardware division (30–100
/// cycles). The previous representation was the bare modulus, which forced
/// every product through `u128` division.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Gfp {
    /// The prime modulus.
    pub p: u64,
    reducer: Barrett,
}

/// Division-free Montgomery arithmetic for odd moduli.
///
/// The implementation is [`crate::kernels::modular::Montgomery`]; the name is
/// kept because the sparse elimination code was written against it.
pub use crate::kernels::modular::Montgomery as MontGfp;

use crate::kernels::modular::Barrett;

/// A Montgomery context for a modulus the row-accumulation kernel can handle.
///
/// An inner product adds up to sixteen products of residues below `p` in a
/// `u128` before reducing once, which stays below `2^128` only while `p < 2^62`.
/// A larger modulus keeps the plain reduction path: the shared `Montgomery`
/// type has no such bound, so the restriction belongs to this caller and is
/// applied here rather than assumed of the context.
fn mont_for_rows(p: u64) -> Option<MontGfp> {
    (p < (1 << 62)).then(|| MontGfp::new(p)).flatten()
}

/// Division-free Barrett arithmetic for any modulus.
///
/// The implementation is [`crate::kernels::modular::Barrett`].
pub use crate::kernels::modular::Barrett as BarrettGfp;

/// Fast modular inverse of `a` modulo `p` by the extended Euclidean algorithm.
///
/// The implementation is [`crate::kernels::modular::mod_inv_u64`]; the name is
/// kept because the sparse elimination code was written against it.
pub use crate::kernels::modular::mod_inv_u64 as fast_inv_mod;

impl Gfp {
    /// The field of integers modulo `p`.
    ///
    /// A modulus below two is not a field; such a value keeps the historical
    /// behaviour of mapping every residue to zero, which the operations below
    /// check before touching the reducer.
    #[must_use]
    pub const fn new(p: u64) -> Self {
        // `Barrett::new_prime` needs a modulus of at least two. The sentinel is
        // never consulted, because every operation guards on `p <= 1` first.
        let modulus = if p < 2 { 2 } else { p };
        Self {
            p,
            reducer: Barrett::new_prime(modulus),
        }
    }

    /// The prime modulus.
    #[inline]
    #[must_use]
    pub const fn prime(self) -> u64 {
        self.p
    }

    /// The division-free reducer for this modulus.
    #[inline]
    #[must_use]
    pub const fn reducer(self) -> Barrett {
        self.reducer
    }

    /// `a * b mod p`.
    #[inline]
    #[must_use]
    pub const fn mul_mod(
        self,
        a: u64,
        b: u64,
    ) -> u64 {
        if self.p <= 1 {
            0
        } else {
            self.reducer.mul(a, b)
        }
    }

    /// `(a + b) mod p`.
    #[inline]
    #[must_use]
    pub const fn add_mod(
        self,
        a: u64,
        b: u64,
    ) -> u64 {
        if self.p <= 1 {
            0
        } else {
            self.reducer.add(a, b)
        }
    }

    /// `(a - b) mod p`.
    #[inline]
    #[must_use]
    pub const fn sub_mod(
        self,
        a: u64,
        b: u64,
    ) -> u64 {
        if self.p <= 1 {
            0
        } else {
            self.reducer.sub(a, b)
        }
    }

    /// `a^e mod p`.
    #[must_use]
    pub const fn pow_mod(
        self,
        a: u64,
        e: u64,
    ) -> u64 {
        if self.p <= 1 {
            return 0;
        }
        self.reducer.pow(a, e)
    }

    /// The inverse of `a` modulo `p`, or `None` for zero.
    #[inline]
    #[must_use]
    pub const fn inv(
        self,
        a: u64,
    ) -> Option<u64> {
        if self.p <= 1 {
            return None;
        }
        self.reducer.inv(a)
    }

    /// The residue of an integer.
    #[must_use]
    pub fn reduce(
        self,
        a: &BigInt,
    ) -> u64 {
        if self.p <= 1 {
            return 0;
        }
        if let Some(v) = a.to_i64() {
            let r = v.unsigned_abs() % self.p;
            return if v < 0 && r != 0 {
                self.p - r
            } else {
                r
            };
        }
        if let Some(v) = a.to_i128() {
            let r = (v.unsigned_abs() % u128::from(self.p)) as u64;
            return if v < 0 && r != 0 {
                self.p - r
            } else {
                r
            };
        }
        let mut rem = (a % BigInt::from(self.p)).to_i128().unwrap_or(0);
        if rem < 0 {
            rem += i128::from(self.p);
        }
        rem as u64
    }

    /// Batch modular inversion of non-zero elements in `values` modulo `p`
    /// using Montgomery's trick. The output buffer holds prefix products
    /// until the reverse pass replaces them with the corresponding inverses.
    ///
    /// Returns `None` if `p <= 1` or any element is divisible by `p`.
    #[must_use]
    pub fn batch_inv(
        self,
        values: &[u64],
    ) -> Option<Vec<u64>> {
        if self.p <= 1 {
            return None;
        }
        if values.is_empty() {
            return Some(Vec::new());
        }
        // These products already have canonical operands. Reuse the wide reduction
        // directly while keeping the Barrett residue contract at the input boundary.
        let product = |a, b| self.reducer.reduce_u128(u128::from(a) * u128::from(b));
        let mut out = Vec::with_capacity(values.len());
        let mut running = 1_u64;
        for &value in values {
            let residue = self.reducer.reduce_u64(value);
            if residue == 0 {
                return None;
            }
            out.push(running);
            running = product(running, residue);
        }
        let mut inverse = self.inv(running)?;
        for (prefix, &value) in out.iter_mut().zip(values).rev() {
            let residue = self.reducer.reduce_u64(value);
            *prefix = product(inverse, *prefix);
            inverse = product(inverse, residue);
        }
        Some(out)
    }
}

impl Domain for Gfp {
    type E = u64;

    #[inline]
    fn zero(&mut self) -> u64 {
        0
    }

    #[inline]
    fn one(&mut self) -> u64 {
        1 % self.p
    }

    #[inline]
    fn is_zero(
        &mut self,
        a: &u64,
    ) -> bool {
        *a == 0
    }

    #[inline]
    fn add(
        &mut self,
        a: &u64,
        b: &u64,
    ) -> u64 {
        let p = self.p;
        let ar = if *a >= p { *a % p } else { *a };
        let br = if *b >= p { *b % p } else { *b };
        let (sum, over) = ar.overflowing_add(br);
        if over || sum >= p {
            sum.wrapping_sub(p)
        } else {
            sum
        }
    }

    #[inline]
    fn sub(
        &mut self,
        a: &u64,
        b: &u64,
    ) -> u64 {
        let p = self.p;
        let ar = if *a >= p { *a % p } else { *a };
        let br = if *b >= p { *b % p } else { *b };
        if ar >= br {
            ar - br
        } else {
            p - (br - ar)
        }
    }

    #[inline]
    fn mul(
        &mut self,
        a: &u64,
        b: &u64,
    ) -> u64 {
        self.mul_mod(*a, *b)
    }

    #[inline]
    fn neg(
        &mut self,
        a: &u64,
    ) -> u64 {
        let r = if *a >= self.p {
            *a % self.p
        } else {
            *a
        };
        if r == 0 { 0 } else { self.p - r }
    }

    #[inline]
    fn div(
        &mut self,
        a: &u64,
        b: &u64,
    ) -> Option<u64> {
        self.inv(*b).map(|i| self.mul_mod(*a, i))
    }
}

impl Domain for Barrett {
    type E = u64;

    #[inline]
    fn zero(&mut self) -> u64 {
        0
    }

    #[inline]
    fn one(&mut self) -> u64 {
        1 % self.modulus()
    }

    #[inline]
    fn is_zero(
        &mut self,
        a: &u64,
    ) -> bool {
        *a == 0
    }

    #[inline]
    fn add(
        &mut self,
        a: &u64,
        b: &u64,
    ) -> u64 {
        Self::add(*self, *a, *b)
    }

    #[inline]
    fn sub(
        &mut self,
        a: &u64,
        b: &u64,
    ) -> u64 {
        Self::sub(*self, *a, *b)
    }

    #[inline]
    fn mul(
        &mut self,
        a: &u64,
        b: &u64,
    ) -> u64 {
        Self::mul(*self, *a, *b)
    }

    #[inline]
    fn neg(
        &mut self,
        a: &u64,
    ) -> u64 {
        Self::neg(*self, *a)
    }

    #[inline]
    fn div(
        &mut self,
        a: &u64,
        b: &u64,
    ) -> Option<u64> {
        Self::div(*self, *a, *b)
    }
}

/// Reduces a sparse rational matrix directly to `GF(p)` using Montgomery's
/// trick batch inversion over all non-trivial denominators, leasing its
/// workspace from [`crate::kernels::reuse::global_pool`].
///
/// Returns `None` if the modulus is below two or any entry has a denominator
/// divisible by it.
#[must_use]
pub fn reduce_rational_matrix_mod_p(
    matrix: &SparseMat<BigRational>,
    field: Gfp,
) -> Option<SparseMat<u64>> {
    if field.p <= 1 {
        return None;
    }
    let nnz = matrix.nnz();
    if nnz == 0 {
        return Some(SparseMat::zeros(matrix.rows, matrix.cols));
    }
    // Residues can occupy all 64 bits, independently of the host's index width.
    let pool = crate::kernels::reuse::global_pool();
    let mut denoms = pool.scratch_i64(nnz);
    let mut numer_residues = pool.scratch_i64(nnz);
    let mut prefix = pool.scratch_i64(nnz);
    denoms.clear();
    numer_residues.clear();
    prefix.clear();

    let mut running = 1_u64;
    let mut has_nontrivial = false;
    for row in &matrix.data {
        for (_, q) in row {
            let num_r = field.reduce(q.numer());
            numer_residues.push(num_r.cast_signed());
            if q.is_integer() {
                denoms.push(1);
                prefix.push(0);
            } else {
                let den_r = field.reduce(q.denom());
                if den_r == 0 {
                    return None;
                }
                if den_r == 1 {
                    denoms.push(1);
                    prefix.push(0);
                } else {
                    has_nontrivial = true;
                    prefix.push(running.cast_signed());
                    running = field.mul_mod(running, den_r);
                    denoms.push(den_r.cast_signed());
                }
            }
        }
    }

    if has_nontrivial {
        let mut inv_acc = field.inv(running)?;
        for idx in (0..nnz).rev() {
            let d = denoms[idx].cast_unsigned();
            if d > 1 {
                let d_inv = field.mul_mod(inv_acc, prefix[idx].cast_unsigned());
                inv_acc = field.mul_mod(inv_acc, d);
                numer_residues[idx] = field
                    .mul_mod(numer_residues[idx].cast_unsigned(), d_inv)
                    .cast_signed();
            }
        }
    }

    let mut flat_idx = 0_usize;
    let mut data = Vec::with_capacity(matrix.rows);
    for row in &matrix.data {
        let mut out_row = Vec::with_capacity(row.len());
        for &(col, _) in row {
            let val = numer_residues[flat_idx].cast_unsigned();
            flat_idx += 1;
            if val != 0 {
                out_row.push((col, val));
            }
        }
        data.push(out_row);
    }

    Some(SparseMat {
        rows: matrix.rows,
        cols: matrix.cols,
        data,
    })
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
        let mut raw_rows: Vec<Vec<(usize, E)>> = vec![Vec::new(); rows];
        for (i, j, v) in triplets {
            if i >= rows || j >= cols {
                return Err(format!(
                    "entry ({i}, {j}) is outside a {rows} x {cols} matrix"
                ));
            }
            raw_rows[i].push((j, v));
        }
        let mut data = Vec::with_capacity(rows);
        for mut row in raw_rows {
            if row.is_empty() {
                data.push(row);
                continue;
            }
            // Stable sort by column preserves input combination order for duplicate positions.
            row.sort_by_key(|&(j, _)| j);
            let mut merged_row: Vec<(usize, E)> = Vec::with_capacity(row.len());
            for (j, v) in row {
                if let Some((last_col, last_val)) = merged_row.last_mut()
                    && *last_col == j
                {
                    *last_val = dom.add(last_val, &v);
                } else {
                    merged_row.push((j, v));
                }
            }
            merged_row.retain(|(_, v)| !dom.is_zero(v));
            data.push(merged_row);
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
    /// (Gustavson's algorithm) using a sparse dense-indexed accumulator with
    /// pooled index scratch, so the cost follows the number of products that
    /// are actually non-zero without `BTreeMap` node allocations.
    #[must_use]
    pub fn matmul<D: Domain<E = E>>(
        &self,
        dom: &mut D,
        other: &Self,
    ) -> Option<Self> {
        if self.cols != other.rows {
            return None;
        }
        let n_cols = other.cols;
        let pool = crate::kernels::reuse::global_pool();
        let mut marker = pool.scratch_usize(n_cols);
        marker.resize(n_cols, usize::MAX);
        let mut touched = pool.scratch_usize(n_cols.min(64));
        let mut values: Vec<E> = Vec::new();
        let mut data = Vec::with_capacity(self.rows);

        for row in &self.data {
            touched.clear();
            values.clear();
            for (k, a) in row {
                for (j, b) in &other.data[*k] {
                    let term = dom.mul(a, b);
                    let pos = marker[*j];
                    if pos == usize::MAX {
                        let idx = values.len();
                        marker[*j] = idx;
                        touched.push(*j);
                        values.push(term);
                    } else {
                        values[pos] = dom.add(&values[pos], &term);
                    }
                }
            }
            touched.sort_unstable();
            let mut out_row = Vec::with_capacity(touched.len());
            for &col in touched.iter() {
                let pos = marker[col];
                marker[col] = usize::MAX;
                let v = &values[pos];
                if !dom.is_zero(v) {
                    out_row.push((col, v.clone()));
                }
            }
            data.push(out_row);
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

    /// Computes and caches a symmetric Reverse Cuthill-McKee bandwidth-reducing
    /// ordering for this matrix's sparsity pattern via
    /// [`crate::kernels::reuse::ComputeReusePool::sparse_ordering_plan`].
    #[must_use]
    pub fn rcm_ordering_plan(&self) -> std::sync::Arc<crate::kernels::reuse::SparseOrderingPlan> {
        let nnz = self.nnz();
        let pool = crate::kernels::reuse::global_pool();
        let mut row_ptrs = pool.scratch_usize(self.rows.saturating_add(1));
        let mut col_indices = pool.scratch_usize(nnz);
        row_ptrs.clear();
        col_indices.clear();
        row_ptrs.push(0);
        for row in &self.data {
            for &(col, _) in row {
                col_indices.push(col);
            }
            row_ptrs.push(col_indices.len());
        }
        pool.sparse_ordering_plan(&row_ptrs, &col_indices)
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

/// Intrusive bucketed priority queue of active rows keyed by `(width, row_index)`,
/// replacing `BTreeSet<(usize, usize)>` with `O(1)` bucket updates and tie-breaking
/// in ascending row order.
struct RowDegreeQueue {
    /// `buckets[w]`: sorted list of active row indices with pivotable width `w`.
    buckets: Vec<Vec<usize>>,
    /// Current recorded width for each row (`0` if inactive or zero width).
    row_width: Vec<usize>,
    /// Lower bound on the minimum non-empty width bucket (`>= 1`).
    min_width: usize,
    /// Number of active rows currently stored in buckets `>= 1`.
    active_count: usize,
}

impl RowDegreeQueue {
    fn new(
        rows: usize,
        max_width: usize,
        initial_widths: &[usize],
    ) -> Self {
        let mut buckets = vec![Vec::new(); max_width.saturating_add(1)];
        let mut row_width = vec![0_usize; rows];
        let mut min_width = max_width.saturating_add(1);
        let mut active_count = 0_usize;
        // Iterating i in 0..rows naturally populates each bucket in ascending row order!
        for (i, &w) in initial_widths.iter().enumerate() {
            row_width[i] = w;
            if w > 0 {
                buckets[w].push(i);
                if w < min_width {
                    min_width = w;
                }
                active_count += 1;
            }
        }
        Self {
            buckets,
            row_width,
            min_width,
            active_count,
        }
    }

    #[inline]
    fn pop_first(&mut self) -> Option<(usize, usize)> {
        if self.active_count == 0 {
            return None;
        }
        while self.min_width < self.buckets.len() && self.buckets[self.min_width].is_empty() {
            self.min_width += 1;
        }
        if self.min_width >= self.buckets.len() {
            return None;
        }
        let w = self.min_width;
        let bucket = &mut self.buckets[w];
        let pr = bucket.remove(0);
        self.row_width[pr] = 0;
        self.active_count -= 1;
        Some((w, pr))
    }

    #[inline]
    fn update_width(
        &mut self,
        row: usize,
        new_width: usize,
    ) {
        let old_width = self.row_width[row];
        if old_width == new_width {
            return;
        }
        if old_width > 0 {
            let bucket = &mut self.buckets[old_width];
            if let Ok(pos) = bucket.binary_search(&row) {
                bucket.remove(pos);
                self.active_count -= 1;
            }
        }
        self.row_width[row] = new_width;
        if new_width > 0 {
            if new_width >= self.buckets.len() {
                self.buckets.resize(new_width.saturating_add(1), Vec::new());
            }
            let bucket = &mut self.buckets[new_width];
            let pos = bucket.partition_point(|&r| r < row);
            bucket.insert(pos, row);
            if new_width < self.min_width {
                self.min_width = new_width;
            }
            self.active_count += 1;
        }
    }
}

struct EliminationWorkspace<E> {
    // A row can be listed more than once if a column disappears and fills in
    // again. Active rows are re-checked before each attempted elimination.
    candidates: Vec<Vec<usize>>,
    counts: Vec<usize>,
    scratch: Vec<(usize, E)>,
    pivot_able_widths: Vec<usize>,
    pivot_cols: usize,
}

impl<E> EliminationWorkspace<E> {
    fn new(
        rows: &[Vec<(usize, E)>],
        cols: usize,
        pivot_cols: usize,
    ) -> Self {
        let mut candidates = vec![Vec::new(); cols];
        let mut counts = vec![0; cols];
        let mut pivot_able_widths = Vec::with_capacity(rows.len());
        for (i, row) in rows.iter().enumerate() {
            let mut w = 0_usize;
            for (j, _) in row {
                candidates[*j].push(i);
                counts[*j] += 1;
                if *j < pivot_cols {
                    w += 1;
                }
            }
            pivot_able_widths.push(w);
        }
        Self {
            candidates,
            counts,
            scratch: Vec::new(),
            pivot_able_widths,
            pivot_cols,
        }
    }

    fn update<D: Domain<E = E>>(
        &mut self,
        dom: &mut D,
        target: usize,
        row: &mut Vec<(usize, E)>,
        minus_factor: &E,
        pivot_row: &[(usize, E)],
        pivot_col: usize,
    ) -> usize {
        self.scratch.reserve(row.len() + pivot_row.len());
        let mut left = row.drain(..).peekable();
        let mut new_pivot_width = 0_usize;
        let pivot_cols = self.pivot_cols;
        for &(col, ref pivot_value) in pivot_row {
            while let Some(entry) = left.next_if(|(j, _)| *j < col) {
                if entry.0 < pivot_cols {
                    new_pivot_width += 1;
                }
                self.scratch.push(entry);
            }
            let term = dom.mul(minus_factor, pivot_value);
            let prior = left.next_if(|(j, _)| *j == col);
            let filled = prior.is_none();
            let value = match prior {
                | Some((_, old)) => dom.add(&old, &term),
                | None => term,
            };
            // The pivot vanishes algebraically even when a symbolic domain
            // cannot see the cancellation in its current normal form.
            if dom.is_zero(&value) || col == pivot_col {
                if !filled {
                    self.counts[col] = self.counts[col].saturating_sub(1);
                }
            } else {
                if filled {
                    self.candidates[col].push(target);
                    self.counts[col] += 1;
                }
                if col < pivot_cols {
                    new_pivot_width += 1;
                }
                self.scratch.push((col, value));
            }
        }
        for entry in left {
            if entry.0 < pivot_cols {
                new_pivot_width += 1;
            }
            self.scratch.push(entry);
        }
        self.pivot_able_widths[target] = new_pivot_width;
        // Draining leaves the original allocation available for the next row.
        std::mem::swap(row, &mut self.scratch);
        new_pivot_width
    }
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
    pub const fn rank(&self) -> usize {
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
    let n = sequence.len();
    if n <= 1 {
        return false;
    }
    let pool = crate::kernels::reuse::global_pool();
    let mut order = pool.scratch_usize(n);
    order.clear();
    order.extend(0..n);
    order.sort_by_key(|&k| sequence[k]);
    // `order[r]` is the position holding the r-th smallest value.
    let mut seen = vec![false; n];
    let mut odd = false;
    for start in 0..n {
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

fn scan_pivot_row(
    widths: &[usize],
    active: &[bool],
) -> Option<(usize, usize)> {
    let mut best = None;
    for (i, &w) in widths.iter().enumerate() {
        if !active[i] {
            continue;
        }
        if w > 0 && best.is_none_or(|(bw, _)| w < bw) {
            best = Some((w, i));
            if w == 1 {
                break;
            }
        }
    }
    best
}

/// Gaussian elimination with sparsity-preserving pivoting.
///
/// Pivots are only taken in columns below `pivot_cols`, which lets a caller
/// append right-hand sides as extra columns. `record` keeps the list of row
/// operations (needed to solve for further right-hand sides later).
/// Pivot ties are broken by the original row and column indices.
pub fn echelon<D: Domain>(
    dom: &mut D,
    matrix: &SparseMat<D::E>,
    pivot_cols: usize,
    record: bool,
) -> Echelon<D::E> {
    echelon_owned(dom, matrix.clone(), pivot_cols, record)
}

fn echelon_owned<D: Domain>(
    dom: &mut D,
    matrix: SparseMat<D::E>,
    pivot_cols: usize,
    record: bool,
) -> Echelon<D::E> {
    const INDEXED_PIVOT_MIN_ROWS: usize = 64;
    let pivot_cols = pivot_cols.min(matrix.cols);
    let mut rows = matrix.data;
    let mut active = vec![true; matrix.rows];
    let mut workspace = EliminationWorkspace::new(&rows, matrix.cols, pivot_cols);
    let mut queue = (rows.len() >= INDEXED_PIVOT_MIN_ROWS)
        .then(|| RowDegreeQueue::new(rows.len(), pivot_cols, &workspace.pivot_able_widths));
    let mut out = Echelon {
        cols: matrix.cols,
        pivots: Vec::with_capacity(matrix.rows.min(pivot_cols)),
        rows: Vec::with_capacity(matrix.rows.min(pivot_cols)),
        ops: Vec::new(),
        rest: Vec::new(),
    };
    loop {
        let best = queue.as_mut().map_or_else(
            || scan_pivot_row(&workspace.pivot_able_widths, &active),
            RowDegreeQueue::pop_first,
        );
        let Some((pivot_width, pr)) = best else {
            break;
        };
        // In that row, the column with the fewest remaining entries.
        let mut pc = rows[pr][0].0;
        let mut best_key = (workspace.counts[pc], pc);
        let mut pivot_index = 0_usize;
        for (idx, &(j, _)) in rows[pr][..pivot_width].iter().enumerate().skip(1) {
            let key = (workspace.counts[j], j);
            if key < best_key {
                best_key = key;
                pc = j;
                pivot_index = idx;
            }
        }
        active[pr] = false;
        workspace.pivot_able_widths[pr] = 0;
        let index = out.pivots.len();
        let pivot_row = std::mem::take(&mut rows[pr]);
        let pivot = &pivot_row[pivot_index].1;
        for (j, _) in &pivot_row {
            workspace.counts[*j] = workspace.counts[*j].saturating_sub(1);
        }
        let mut targets = std::mem::take(&mut workspace.candidates[pc]);
        for &target in &targets {
            if !active[target] {
                continue;
            }
            let Some(value) = entry(&rows[target], pc) else {
                continue;
            };
            let Some(factor) = dom.div(value, pivot) else {
                continue;
            };
            let minus = dom.neg(&factor);
            let updated_width =
                workspace.update(dom, target, &mut rows[target], &minus, &pivot_row, pc);
            if let Some(q) = &mut queue {
                q.update_width(target, updated_width);
            }
            if record {
                out.ops.push((target, index, factor));
            }
        }
        targets.clear();
        workspace.candidates[pc] = targets;
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

/// Reduced row echelon form (Gauss-Jordan form) of a sparse matrix over `dom`.
///
/// Scales every pivot row so its pivot entry is `1` and eliminates each pivot
/// column from all earlier pivot rows as well as subsequent rows. Returns the
/// reduced [`SparseMat<D::E>`] with pivot rows placed in canonical order at the
/// top followed by any zero rows, together with the sorted list of pivot columns.
pub fn rref<D: Domain>(
    dom: &mut D,
    matrix: &SparseMat<D::E>,
) -> (SparseMat<D::E>, Vec<usize>) {
    let mut e = echelon(dom, matrix, matrix.cols, false);
    let rank = e.rank();
    if rank == 0 {
        return (SparseMat::zeros(matrix.rows, matrix.cols), Vec::new());
    }
    // Order pivot rows by ascending pivot column.
    let mut perm: Vec<usize> = (0..rank).collect();
    perm.sort_unstable_by_key(|&k| e.pivots[k].1);

    let mut piv_cols: Vec<usize> = Vec::with_capacity(rank);
    let mut piv_rows: Vec<Vec<(usize, D::E)>> = Vec::with_capacity(rank);
    for k in perm {
        let pc = e.pivots[k].1;
        let row = std::mem::take(&mut e.rows[k]);
        piv_cols.push(pc);
        piv_rows.push(row);
    }

    // Normalize each pivot row so its pivot entry is 1.
    for k in 0..rank {
        let pc = piv_cols[k];
        if let Some(piv_val) = entry(&piv_rows[k], pc).cloned()
            && let Some(inv) = {
                // `one` takes the domain mutably, so it is materialised before
                // the borrow `div` needs.
                let one = dom.one();
                dom.div(&one, &piv_val)
            }
        {
            for (col, val) in &mut piv_rows[k] {
                if *col == pc {
                    *val = dom.one();
                } else {
                    *val = dom.mul(val, &inv);
                }
            }
        }
    }

    // Eliminate upward from k = rank - 1 down to 0.
    for k in (0..rank).rev() {
        let pc = piv_cols[k];
        let pivot_row = piv_rows[k].clone();
        for target in 0..k {
            if let Some(factor) = entry(&piv_rows[target], pc).cloned() {
                let minus = dom.neg(&factor);
                let mut merged = axpy_rows(dom, &piv_rows[target], &minus, &pivot_row);
                // Ensure the pivot column itself is removed even in unsimplified symbolic domains.
                merged.retain(|(col, _)| *col != pc);
                piv_rows[target] = merged;
            }
        }
    }

    let mut data = piv_rows;
    data.resize(matrix.rows, Vec::new());
    (
        SparseMat {
            rows: matrix.rows,
            cols: matrix.cols,
            data,
        },
        piv_cols,
    )
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

/// Solves `A x = b`.
///
/// Returns `None` when the matrix and right-hand side row counts disagree
/// or the augmented column count cannot be represented by the workspace.
pub fn solve<D: Domain>(
    dom: &mut D,
    matrix: &SparseMat<D::E>,
    b: &[D::E],
) -> Option<Solution<D::E>> {
    if b.len() != matrix.rows || matrix.data.len() != matrix.rows {
        return None;
    }
    let n = matrix.cols;
    let augmented_cols = n.checked_add(1)?;
    Layout::array::<Vec<usize>>(augmented_cols).ok()?;
    Layout::array::<D::E>(n).ok()?;
    let mut augmented = matrix.clone();
    augmented.cols = augmented_cols;
    for (row, value) in augmented.data.iter_mut().zip(b) {
        if !dom.is_zero(value) {
            row.push((n, value.clone()));
        }
    }
    let e = echelon_owned(dom, augmented, n, false);
    if e.rest.iter().any(|(_, row)| !row.is_empty()) {
        return Some(Solution::Inconsistent);
    }
    let zero = dom.zero();
    let rhs: Vec<D::E> = e
        .rows
        .iter()
        .map(|row| entry(row, n).cloned().unwrap_or_else(|| zero.clone()))
        .collect();
    let mut particular = vec![zero; n];
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
pub const fn is_prime_u64(n: u64) -> bool {
    super::number_theory::is_prime_miller_rabin(n)
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
///
/// Evaluates word-sized prime moduli in parallel via [`crate::kernels::par::map_indexed`]
/// when [`crate::execution::should_split_work`] indicates enough work, and folds
/// the modular determinants in deterministic prime order via incremental CRT.
#[must_use]
pub fn det_rational(matrix: &SparseMat<BigRational>) -> Option<BigRational> {
    if matrix.rows != matrix.cols {
        return None;
    }
    if matrix.rows == 0 {
        return Some(BigRational::one());
    }
    let (integer, scales) = integer_rows(matrix);
    let det = det_over_integers(&integer)?;
    let scale: BigInt = scales.iter().product();
    Some(BigRational::new(det, scale))
}

/// Exact modular determinant of a square integer matrix.
fn det_over_integers(integer: &SparseMat<BigInt>) -> Option<BigInt> {
    if integer.data.iter().any(Vec::is_empty) {
        return Some(BigInt::zero());
    }
    // |det| <= product of the row norms.
    let bound: BigInt = integer
        .data
        .iter()
        .map(|row| norm_bound(row.iter().map(|(_, v)| v)))
        .product();
    let target = bound * BigInt::from(2) + BigInt::one();

    // Pre-collect the exact list of word primes needed to exceed `target`.
    let mut primes = Vec::new();
    let mut prod_check = BigInt::one();
    for p in word_primes() {
        primes.push(p);
        prod_check *= BigInt::from(p);
        if prod_check > target {
            break;
        }
    }

    let num_primes = primes.len();
    let nnz = integer.nnz().max(integer.rows);
    let per_prime_work = nnz.saturating_mul(integer.rows.max(1));
    let residues: Vec<Option<u64>> =
        if crate::execution::should_split_work(num_primes, per_prime_work) {
            crate::kernels::par::map_indexed(num_primes, 2, |idx| {
                let p = primes[idx];
                let mut field = Gfp::new(p);
                let reduced = reduce_matrix(integer, field);
                det(&mut field, &reduced)
            })
        } else {
            primes
                .iter()
                .map(|&p| {
                    let mut field = Gfp::new(p);
                    let reduced = reduce_matrix(integer, field);
                    det(&mut field, &reduced)
                })
                .collect()
        };

    let mut modulus = BigInt::one();
    let mut residue = BigInt::zero();
    for (&p, maybe_d) in primes.iter().zip(residues) {
        let d = maybe_d?;
        let field = Gfp::new(p);
        let mut dom = field;
        // Chinese remainder: residue + modulus * t = d (mod p).
        let have = field.reduce(&residue);
        let diff = dom.sub(&d, &have);
        let inv = field.inv(field.reduce(&modulus))?;
        let t = field.mul_mod(diff, inv);
        residue += &modulus * BigInt::from(t);
        modulus *= BigInt::from(p);
    }
    // The symmetric representative.
    if &residue * BigInt::from(2) > modulus {
        residue -= &modulus;
    }
    Some(residue)
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
///
/// Precomputes the inverse of every pivot in the echelon form and uses
/// Montgomery/Barrett arithmetic plus pooled scratch buffers for each solve step.
struct ModularFactor {
    field: Gfp,
    echelon: Echelon<u64>,
    pivot_invs: Vec<u64>,
    n: usize,
}

impl ModularFactor {
    fn new(
        matrix: &SparseMat<u64>,
        field: Gfp,
    ) -> Option<Self> {
        let mut dom = field;
        let echelon = echelon(&mut dom, matrix, matrix.cols, true);
        if echelon.rank() != matrix.rows {
            return None;
        }
        let mut pivot_vals = Vec::with_capacity(echelon.pivots.len());
        for (row, &(_, pc)) in echelon.rows.iter().zip(&echelon.pivots) {
            pivot_vals.push(*entry(row, pc)?);
        }
        let pivot_invs = field.batch_inv(&pivot_vals)?;
        Some(Self {
            field,
            echelon,
            pivot_invs,
            n: matrix.cols,
        })
    }

    /// Apply the row operations to `rhs` and overwrite `solution` with `A^-1 b mod p`.
    fn solve(
        &self,
        rhs: &mut [u64],
        solution: &mut [u64],
    ) {
        let mut dom = self.field;
        for &(target, index, factor) in &self.echelon.ops {
            let source = rhs[self.echelon.pivots[index].0];
            if source != 0 {
                let scaled = dom.mul(&factor, &source);
                let cur = rhs[target];
                rhs[target] = dom.sub(&cur, &scaled);
            }
        }
        solution.fill(0);
        for k in (0..self.echelon.pivots.len()).rev() {
            let (pr, pc) = self.echelon.pivots[k];
            let mut acc = rhs[pr];
            for &(j, v) in &self.echelon.rows[k] {
                if j != pc && j < self.n {
                    let xj = solution[j];
                    if xj != 0 {
                        let term = dom.mul(&v, &xj);
                        acc = dom.sub(&acc, &term);
                    }
                }
            }
            solution[pc] = dom.mul(&acc, &self.pivot_invs[k]);
        }
    }
}

/// Flat CSR representation of an integer sparse matrix for fast `A * d` products
/// with `u64` digit vectors during Dixon `p`-adic lifting.
struct IntegerCsr {
    rows: usize,
    row_ptr: Vec<usize>,
    col_idx: Vec<usize>,
    values: Vec<BigInt>,
    small_i64: Option<Vec<i64>>,
}

impl IntegerCsr {
    fn from_sparse(matrix: &SparseMat<BigInt>) -> Self {
        let nnz = matrix.nnz();
        let mut row_ptr = Vec::with_capacity(matrix.rows.saturating_add(1));
        let mut col_idx = Vec::with_capacity(nnz);
        let mut values = Vec::with_capacity(nnz);
        let mut small_vals = Vec::with_capacity(nnz);
        let mut all_small = true;
        row_ptr.push(0);
        for row in &matrix.data {
            for (col, val) in row {
                col_idx.push(*col);
                if all_small {
                    if let Some(v64) = val.to_i64() {
                        small_vals.push(v64);
                    } else {
                        all_small = false;
                    }
                }
                values.push(val.clone());
            }
            row_ptr.push(col_idx.len());
        }
        Self {
            rows: matrix.rows,
            row_ptr,
            col_idx,
            values,
            small_i64: all_small.then_some(small_vals),
        }
    }

    fn mul_digits(
        &self,
        digits: &[u64],
        out: &mut Vec<BigInt>,
    ) {
        let n = self.rows;
        let nnz = self.col_idx.len();
        let grain = if crate::execution::should_split_work(n, (nnz / n.max(1)).max(16)) {
            64
        } else {
            usize::MAX
        };
        if let Some(small) = &self.small_i64 {
            crate::kernels::par::map_indexed_into(n, grain, out, |i| {
                let start = self.row_ptr[i];
                let end = self.row_ptr[i + 1];
                let mut acc = 0_i128;
                for k in start..end {
                    let d = digits[self.col_idx[k]];
                    if d != 0 {
                        acc += i128::from(small[k]) * i128::from(d);
                    }
                }
                BigInt::from(acc)
            });
        } else {
            crate::kernels::par::map_indexed_into(n, grain, out, |i| {
                let start = self.row_ptr[i];
                let end = self.row_ptr[i + 1];
                let mut acc = BigInt::zero();
                for k in start..end {
                    let d = digits[self.col_idx[k]];
                    if d != 0 {
                        acc += &self.values[k] * BigInt::from(d);
                    }
                }
                acc
            });
        }
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
    if b.len() != matrix.rows {
        return None;
    }
    let mut sols = solve_dixon_multi(matrix, std::slice::from_ref(&b.to_vec()))?;
    sols.pop()
}

/// Solves a square nonsingular rational system with multiple right-hand sides.
///
/// Dixon’s `p`-adic lifting shares one modular factorization of `A` and one flat CSR workspace
/// across all columns of `rhs`.
///
/// Returns `None` when `matrix` is not square, any vector in `rhs` does not
/// have length `matrix.rows`, or `matrix` is singular.
#[must_use]
pub fn solve_dixon_multi(
    matrix: &SparseMat<BigRational>,
    rhs: &[Vec<BigRational>],
) -> Option<Vec<Vec<BigRational>>> {
    let n = matrix.rows;
    if n != matrix.cols || rhs.iter().any(|b| b.len() != n) {
        return None;
    }
    if rhs.is_empty() {
        return Some(Vec::new());
    }
    if n == 0 {
        return Some(vec![Vec::new(); rhs.len()]);
    }
    let (integer, scales) = integer_rows(matrix);
    let factor_at = |prime| {
        let field = Gfp::new(prime);
        ModularFactor::new(&reduce_matrix(&integer, field), field).map(|factor| (field, factor))
    };
    let (field, factor) = word_primes().take(4).find_map(factor_at).or_else(|| {
        // A nonsingular rational matrix can be singular at all trial primes.
        // Its exact determinant distinguishes that case from singularity.
        let determinant = det_over_integers(&integer)?;
        if determinant.is_zero() {
            return None;
        }
        let compatible = word_primes()
            .skip(4)
            .find(|&prime| Gfp::new(prime).reduce(&determinant) != 0)?;
        factor_at(compatible)
    })?;
    let p = BigInt::from(field.p);

    // Cramer's rule with Hadamard's inequality bounds numerators and
    // denominators of the solution of the integer system.
    let transposed = integer.transpose();
    let hadamard: BigInt = transposed
        .data
        .iter()
        .map(|column| norm_bound(column.iter().map(|(_, v)| v)))
        .product();
    let hadamard_sq_4 = &hadamard * &hadamard * BigInt::from(4);
    let csr = IntegerCsr::from_sparse(&integer);

    let lift_one = |b: &Vec<BigRational>| -> Option<Vec<BigRational>> {
        let scaled: Vec<BigRational> = b
            .iter()
            .zip(&scales)
            .map(|(v, s)| v * BigRational::from(s.clone()))
            .collect();
        let common = scaled
            .iter()
            .fold(BigInt::one(), |acc, v| acc.lcm(v.denom()));
        let int_rhs: Vec<BigInt> = scaled
            .iter()
            .map(|v| (v * BigRational::from(common.clone())).to_integer())
            .collect();
        let limit = &hadamard_sq_4 * norm_bound(int_rhs.iter());

        let mut modular_rhs = vec![0_u64; n];
        let mut digits = vec![0_u64; n];
        let mut residual = int_rhs.clone();
        let mut lifted = vec![BigInt::zero(); n];
        let mut product = Vec::with_capacity(n);
        let mut power = BigInt::one();
        let mut next_try = 4_usize;
        let mut step = 0_usize;

        loop {
            for (residue, value) in modular_rhs.iter_mut().zip(&residual) {
                *residue = field.reduce(value);
            }
            factor.solve(&mut modular_rhs, &mut digits);
            csr.mul_digits(&digits, &mut product);
            for ((r, ax), (x, &d)) in residual
                .iter_mut()
                .zip(&product)
                .zip(lifted.iter_mut().zip(&digits))
            {
                // Exact: r = A d (mod p) by construction.
                *r = (&*r - ax) / &p;
                if d != 0 {
                    *x += BigInt::from(d) * &power;
                }
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
                    && verifies(&integer, &int_rhs, &x)
                {
                    let unscale = BigRational::from(common);
                    return Some(x.into_iter().map(|v| v / &unscale).collect());
                }
                if done {
                    return None;
                }
            }
        }
    };

    let k = rhs.len();
    let per_rhs_work = integer.nnz().max(n).saturating_mul(8);
    if crate::execution::should_split_work(k, per_rhs_work) {
        crate::kernels::par::try_map_indexed(k, 2, |idx| lift_one(&rhs[idx]).ok_or(())).ok()
    } else {
        rhs.iter().map(lift_one).collect()
    }
}

/// Flattened Compressed Sparse Row (CSR) matrix in Montgomery domain for
/// division-free Wiedemann Krylov sequences.
///
/// Each row product accumulates `u128` products in Montgomery representation
/// and performs a single `MontGfp::redc` at the end of the row (chunked every
/// 16 entries so `u128` never overflows for `p < 2^62`).
struct MontgomeryCsr {
    rows: usize,
    cols: usize,
    row_ptr: Vec<usize>,
    col_idx: Vec<usize>,
    val_mont: Vec<u64>,
    mont: Option<MontGfp>,
    field: Gfp,
}

impl MontgomeryCsr {
    fn new(
        matrix: &SparseMat<u64>,
        field: Gfp,
    ) -> Self {
        let nnz = matrix.nnz();
        let mut row_ptr = Vec::with_capacity(matrix.rows.saturating_add(1));
        let mut col_idx = Vec::with_capacity(nnz);
        let mut val_mont = Vec::with_capacity(nnz);
        let mont = mont_for_rows(field.p);
        row_ptr.push(0);
        for row in &matrix.data {
            for &(c, v) in row {
                let vr = v % field.p;
                if vr != 0 {
                    col_idx.push(c);
                    val_mont.push(mont.map_or(vr, |m| m.to_mont(vr)));
                }
            }
            row_ptr.push(col_idx.len());
        }
        Self {
            rows: matrix.rows,
            cols: matrix.cols,
            row_ptr,
            col_idx,
            val_mont,
            mont,
            field,
        }
    }

    #[inline]
    fn row_dot(
        &self,
        row: usize,
        x: &[usize],
    ) -> usize {
        let start = self.row_ptr[row];
        let end = self.row_ptr[row + 1];
        if start == end {
            return 0;
        }
        let cols = &self.col_idx[start..end];
        let vals = &self.val_mont[start..end];
        self.mont.map_or_else(
            || {
                let p = u128::from(self.field.p);
                let mut sum = 0_u128;
                for (&c, &v) in cols.iter().zip(vals) {
                    sum = (sum + u128::from(v) * (x[c] as u128)) % p;
                }
                (sum as u64) as usize
            },
            |m| {
                // For p < 2^62, each product is < 2^124, so up to 16 products fit in u128 (< 2^128).
                // Chunking by 16 and taking redc or modulo p*2^64 avoids all overflow and hardware division!
                let mut total = 0_u64;
                for (c_chunk, v_chunk) in cols.chunks(16).zip(vals.chunks(16)) {
                    let mut sum = 0_u128;
                    for (&c, &v) in c_chunk.iter().zip(v_chunk) {
                        sum += u128::from(v) * (x[c] as u128);
                    }
                    // Split 128-bit sum into low 64 bits and high 64 bits:
                    // redc(lo) + hi mod p = sum * R^{-1} mod p.
                    let lo = u128::from(sum as u64);
                    let hi = (sum >> 64) as u64;
                    let red_lo = m.redc(lo);
                    let red_hi = if hi >= m.modulus() {
                        hi % m.modulus()
                    } else {
                        hi
                    };
                    let chunk_val = m.add_mont(red_lo, red_hi);
                    total = m.add_mont(total, chunk_val);
                }
                total as usize
            },
        )
    }

    fn spmv_into(
        &self,
        x: &[usize],
        y: &mut Vec<usize>,
    ) {
        let n = self.rows;
        let nnz = self.col_idx.len();
        let work_per_row = (nnz / n.max(1)).max(4);
        let grain = if crate::execution::should_split_work(n, work_per_row) {
            128
        } else {
            usize::MAX
        };
        crate::kernels::par::map_indexed_into(n, grain, y, |i| self.row_dot(i, x));
    }

    #[inline]
    fn dot(
        &self,
        u: &[u64],
        w: &[usize],
    ) -> u64 {
        self.mont.map_or_else(
            || {
                let mut dom = self.field;
                let mut acc = 0_u64;
                for (&a, &b) in u.iter().zip(w) {
                    let term = dom.mul(&a, &(b as u64));
                    acc = dom.add(&acc, &term);
                }
                acc
            },
            |m| {
                // u is in standard domain, w is instandard domain (< p).
                // Delayed u128 accumulation with chunks of 16, single reduction per chunk!
                let mut total = 0_u64;
                for (u_chunk, w_chunk) in u.chunks(16).zip(w.chunks(16)) {
                    let mut sum = 0_u128;
                    for (&ui, &wi) in u_chunk.iter().zip(w_chunk) {
                        sum += u128::from(ui) * (wi as u128);
                    }
                    let lo = u128::from(sum as u64);
                    let hi = (sum >> 64) as u64;
                    let red_lo = m.redc(lo);
                    let red_hi = if hi >= m.modulus() {
                        hi % m.modulus()
                    } else {
                        hi
                    };
                    // chunk_val = sum * R^{-1} mod p; we convert back via mul_mont(..., r2_mod_p) at the end!
                    let chunk_val = m.add_mont(red_lo, red_hi);
                    total = m.add_mont(total, chunk_val);
                }
                m.mul_mont(total, m.r2_mod())
            },
        )
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
    let mont = mont_for_rows(field.p);
    // c(x) = 1 + c1 x + ...: s[i] + c1 s[i-1] + ... = 0.
    let mut c: Vec<u64> = vec![1];
    let mut prev: Vec<u64> = vec![1];
    let mut len = 0_usize;
    let mut shift = 1_usize;
    let mut last = 1_u64;
    for i in 0..sequence.len() {
        let mut d = sequence[i] % field.p;
        let bound = len.min(c.len() - 1);
        if let Some(m) = mont
            && bound > 0
        {
            // Delayed u128 accumulation in chunks of 16, single reduction per chunk!
            let mut total = 0_u64;
            let mut k = 1_usize;
            while k <= bound {
                let end = (k + 16).min(bound + 1);
                let mut sum = 0_u128;
                for idx in k..end {
                    sum += u128::from(c[idx]) * u128::from(sequence[i - idx] % field.p);
                }
                let lo = u128::from(sum as u64);
                let hi = (sum >> 64) as u64;
                let red_lo = m.redc(lo);
                let red_hi = if hi >= m.modulus() {
                    hi % m.modulus()
                } else {
                    hi
                };
                let chunk_val = m.add_mont(red_lo, red_hi);
                total = m.add_mont(total, chunk_val);
                k = end;
            }
            let inner = m.mul_mont(total, m.r2_mod());
            d = m.add_mont(d, inner);
        } else {
            for k in 1..=bound {
                let term = dom.mul(&c[k], &sequence[i - k]);
                d = dom.add(&d, &term);
            }
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
        if let Some(m) = mont {
            let coef_mont = m.to_mont(coef);
            for (k, &v) in prev.iter().enumerate() {
                let term = m.mul_mont(coef_mont, v);
                c[k + shift] = m.sub_mont(c[k + shift], term);
            }
        } else {
            for (k, v) in prev.iter().enumerate() {
                let term = dom.mul(&coef, v);
                c[k + shift] = dom.sub(&c[k + shift], &term);
            }
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
    const fn next(
        &mut self,
        field: Gfp,
    ) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0 % field.p
    }
}

/// The minimal polynomial of the sequence `u^T A^i v`, `i = 0, 1, ...`,
/// using the flattened Montgomery CSR representation and pooled Krylov buffers.
fn projected_minpoly_csr(
    csr: &MontgomeryCsr,
    u: &[u64],
    v: &[u64],
) -> Option<Vec<u64>> {
    if u.len() != csr.rows || v.len() != csr.cols {
        return None;
    }
    let n = csr.rows;
    let pool = crate::kernels::reuse::global_pool();
    let mut w = pool.scratch_usize(n);
    let mut next_w = pool.scratch_usize(n);
    w.clear();
    w.extend(v.iter().map(|&val| (val % csr.field.p) as usize));
    let steps = 2 * n + 2;
    let mut sequence = Vec::with_capacity(steps);
    for step in 0..steps {
        sequence.push(csr.dot(u, &w));
        if step + 1 < steps {
            csr.spmv_into(&w, &mut next_w);
            std::mem::swap(&mut *w, &mut *next_w);
        }
    }
    Some(berlekamp_massey(&sequence, csr.field))
}

/// The minimal polynomial of a square matrix over `GF(p)` by Wiedemann's
/// algorithm, monic, lowest degree first.
///
/// The algorithm flattens `matrix` once into a Montgomery-domain CSR
/// workspace, multiplies by the matrix `2n` times per trial with delayed
/// `u128` row accumulation, and leases Krylov vectors from
/// [`crate::kernels::reuse::global_pool`]: no entry of the matrix is ever
/// changed, so there is no fill-in. One random projection yields a divisor
/// of the minimal polynomial, and the right one with probability at least
/// `1 - 2n/p`; the result here is the least common multiple over `trials`
/// projections (seeded, so the function is deterministic). `None` when the
/// matrix is not square.
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
    let csr = MontgomeryCsr::new(matrix, field);
    let mut rng = Lcg(seed | 1);
    let mut best: Vec<u64> = vec![1];
    for _ in 0..trials.max(1) {
        let u: Vec<u64> = (0..n).map(|_| rng.next(field)).collect();
        let v: Vec<u64> = (0..n).map(|_| rng.next(field)).collect();
        let m = projected_minpoly_csr(&csr, &u, &v)?;
        best = poly_lcm(&best, &m, field);
    }
    Some(best)
}

/// Solves a square non-singular system `A x = b` over `GF(p)` by Wiedemann's
/// algorithm.
///
/// Flattens `matrix` once into a Montgomery-domain CSR workspace and leases
/// Krylov vectors from [`crate::kernels::reuse::global_pool`].
/// Returns `None` when the projection does not determine a solution
/// (singular matrix, or an unlucky projection in each of the `trials`). A
/// returned vector always satisfies `A x = b`.
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
    let csr = MontgomeryCsr::new(matrix, field);
    let pool = crate::kernels::reuse::global_pool();
    let mut w = pool.scratch_usize(n);
    let mut next_w = pool.scratch_usize(n);
    let mut dom = field;
    let mut rng = Lcg(seed | 1);
    let b_reduced: Vec<u64> = b.iter().map(|&v| v % field.p).collect();

    for _ in 0..trials.max(1) {
        let u: Vec<u64> = (0..n).map(|_| rng.next(field)).collect();
        let m = projected_minpoly_csr(&csr, &u, &b_reduced)?;
        // m(A) b = 0 with m(0) != 0 gives x = -(1/m0) sum_{k>=1} m_k A^(k-1) b.
        let Some(m0_inv) = m.first().and_then(|&c| field.inv(c)) else {
            continue;
        };
        let mut x = vec![0_u64; n];
        w.clear();
        w.extend(b_reduced.iter().map(|&v| v as usize));
        let coeffs = &m[1..];
        for (step, &coefficient) in coeffs.iter().enumerate() {
            if coefficient != 0 {
                if let Some(mont) = csr.mont {
                    let c_mont = mont.to_mont(coefficient);
                    for (xi, &wi) in x.iter_mut().zip(w.iter()) {
                        let term = mont.mul_mont(c_mont, wi as u64);
                        *xi = mont.add_mont(*xi, term);
                    }
                } else {
                    for (xi, &wi) in x.iter_mut().zip(w.iter()) {
                        let term = dom.mul(&coefficient, &(wi as u64));
                        *xi = dom.add(xi, &term);
                    }
                }
            }
            if step + 1 < coeffs.len() {
                csr.spmv_into(&w, &mut next_w);
                std::mem::swap(&mut *w, &mut *next_w);
            }
        }
        let scale = dom.neg(&m0_inv);
        if let Some(mont) = csr.mont {
            let scale_mont = mont.to_mont(scale);
            for xi in &mut x {
                *xi = mont.mul_mont(*xi, scale_mont);
            }
        } else {
            for xi in &mut x {
                *xi = dom.mul(xi, &scale);
            }
        }
        w.clear();
        w.extend(x.iter().map(|&v| v as usize));
        csr.spmv_into(&w, &mut next_w);
        if next_w
            .iter()
            .zip(&b_reduced)
            .all(|(&got, &want)| (got as u64) == want)
        {
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
                    assert_eq!(nullspace, Vec::<Vec<BigRational>>::new());
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
        let field = Gfp::new(101);
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
        let field = Gfp::new(1_000_000_007);
        let mut s = vec![0_u64, 1];
        for i in 2..20 {
            s.push((s[i - 1] + s[i - 2]) % field.p);
        }
        // s[i+2] - s[i+1] - s[i] = 0.
        assert_eq!(
            berlekamp_massey(&s, field),
            vec![field.p - 1, field.p - 1, 1]
        );
        assert_eq!(berlekamp_massey(&[0, 0, 0], field), vec![1]);
    }

    #[test]
    fn wiedemann_agrees_with_elimination() {
        let field = Gfp::new(word_primes().next().unwrap());
        let rational = sample(40, 3, 11);
        let (integer, _) = integer_rows(&rational);
        let a = reduce_matrix(&integer, field);
        let b: Vec<u64> = (0..40).map(|i| (i * i + 1) as u64).collect();
        let x = solve_wiedemann(&a, &b, field, 3, 42).expect("non-singular");
        let mut dom = field;
        assert_eq!(a.matvec(&mut dom, &x).unwrap(), b);
        // The same solution as the factorisation gives (it is unique).
        let mut rhs = b;
        let mut direct = vec![0; rhs.len()];
        ModularFactor::new(&a, field)
            .unwrap()
            .solve(&mut rhs, &mut direct);
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
            vec![0, field.p - 1, 1]
        );
        assert!(solve_wiedemann(&projection, &[1, 2, 3, 4], field, 3, 5).is_none());
        assert!(minpoly_wiedemann(&SparseMat::zeros(2, 3), field, 1, 1).is_none());
    }

    /// The field arithmetic must agree with `u128` arithmetic for every prime
    /// the engine accepts, including the ones near `2^62` where a careless
    /// Barrett estimate is off by more than one.
    #[test]
    fn gfp_arithmetic_matches_the_u128_reference() {
        for p in [
            2_u64,
            3,
            97,
            1_000_003,
            1_000_000_007,
            (1 << 31) - 1,
            (1 << 61) - 1,
            (1 << 62) - 57,
        ] {
            let field = Gfp::new(p);
            assert_eq!(field.prime(), p);
            let mut s = 0x2545_F491_4F6C_DD1D_u64;
            let mut next = move || {
                s ^= s << 13;
                s ^= s >> 7;
                s ^= s << 17;
                s
            };
            for _ in 0..200 {
                let a = next() % p;
                let b = next() % p;
                assert_eq!(
                    u128::from(field.mul_mod(a, b)),
                    u128::from(a) * u128::from(b) % u128::from(p),
                    "p = {p}, {a} * {b}"
                );
                assert_eq!(
                    u128::from(field.add_mod(a, b)),
                    (u128::from(a) + u128::from(b)) % u128::from(p)
                );
                let e = next() % 64;
                let mut want = 1_u128 % u128::from(p);
                for _ in 0..e {
                    want = want * u128::from(a) % u128::from(p);
                }
                assert_eq!(field.pow_mod(a, e), want as u64, "p = {p}, {a}^{e}");
                if a != 0 {
                    let inv = field.inv(a).expect("non-zero is invertible modulo a prime");
                    assert_eq!(field.mul_mod(a, inv), 1, "p = {p}, a = {a}");
                }
            }
            assert!(field.inv(0).is_none());
        }
        // A modulus below two is not a field; the old representation answered
        // zero rather than panicking, and so does this one.
        let degenerate = Gfp::new(1);
        assert_eq!(degenerate.mul_mod(3, 4), 0);
        assert_eq!(degenerate.pow_mod(3, 2), 0);
        assert!(degenerate.inv(3).is_none());
    }

    /// Batch inversion must equal the element-wise inverses, and the prefix
    /// products must survive a modulus close to the word size.
    #[test]
    fn gfp_batch_inverse_matches_individual_inverses() {
        for p in [7_u64, 1_000_000_007, (1 << 61) - 1] {
            let field = Gfp::new(p);
            // Reduce into `[1, p - 1]`: a zero residue has no inverse, and the
            // point of the test is the prefix-product trick, not the rejection.
            let values: Vec<u64> = (1..=64).map(|v| (v * 7_919 + 3) % (p - 1) + 1).collect();
            let batch = field.batch_inv(&values).expect("all values are non-zero");
            for (i, &v) in values.iter().enumerate() {
                assert_eq!(field.mul_mod(v, batch[i]), 1, "p = {p}, index {i}");
            }
            assert!(field.batch_inv(&[1, 0, 2]).is_none());
        }
    }

    /// The Barrett reducer is a [`Domain`] in its own right: a bare
    /// [`BarrettGfp`] must drive the same elimination as the field wrapper,
    /// which is what lets a caller skip the wrapper entirely.
    #[test]
    fn barrett_domain_eliminates_like_gfp() {
        let p = 1_000_003_u64;
        let mut field = Gfp::new(p);
        let mut reducer = Barrett::new_prime(p);
        let triplets = vec![(0, 0, 2_u64), (0, 1, 1), (1, 0, 1), (1, 1, 3), (2, 2, 5)];
        let via_field = SparseMat::from_triplets(&mut field, 3, 3, triplets.clone()).unwrap();
        let via_reducer = SparseMat::from_triplets(&mut reducer, 3, 3, triplets).unwrap();
        assert_eq!(via_field, via_reducer);
        assert_eq!(
            rank(&mut field, &via_field),
            rank(&mut reducer, &via_reducer)
        );
        // 2x + y = 5, x + 3y = 10  =>  x = 1, y = 3  (mod 1000003)
        let b = vec![5_u64, 10, 0];
        assert_eq!(
            solve(&mut field, &via_field, &b).unwrap(),
            solve(&mut reducer, &via_reducer, &b).unwrap()
        );
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
