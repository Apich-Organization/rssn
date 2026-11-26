//! Rank-revealing QR, randomized SVD, LOBPCG and Anderson acceleration.
//!
//! Rank-revealing and matrix-free numerical linear algebra, and fixed-point
//! acceleration.
//!
//! The dense routines in [`crate::kernels::dense`] factor a whole matrix. The
//! algorithms here answer the questions that come up when that is the wrong
//! tool:
//!
//! * [`qr_pivoted`]: column-pivoted Householder QR (Businger-Golub), the
//!   standard *rank-revealing* factorisation, with a numerical rank and the
//!   basic least-squares solution of a rank-deficient system;
//! * [`interpolative_decomposition`]: a low-rank approximation in terms of
//!   a subset of the original matrix columns;
//! * [`randomized_svd`]: a rank-`k` approximation of an operator that is only
//!   available as products `A x` and `A^T y` (Halko, Martinsson and Tropp),
//!   in `O(k)` products instead of a full decomposition;
//! * [`lobpcg`]: the smallest eigenpairs of a large symmetric operator by
//!   locally optimal block preconditioned conjugate gradients (Knyazev), which
//!   needs products and an optional preconditioner but no factorisation;
//!   [`lobpcg_generalized`] uses a positive-definite mass operator and
//!   optional soft locking for converged approximations;
//! * [`anderson`]: Anderson acceleration of a fixed-point iteration
//!   `x = g(x)`, which turns the linear convergence of a self-consistent
//!   iteration into something close to a quasi-Newton method without a
//!   Jacobian.
//!
//! Every routine is deterministic (random test matrices come from a seeded
//! generator), reports failure through its result rather than by panicking,
//! and treats non-finite input as an error.

#![allow(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_precision_loss,
    clippy::many_single_char_names,
    clippy::module_name_repetitions,
    clippy::single_call_fn,
    clippy::missing_docs_in_private_items,
    clippy::needless_range_loop,
    reason = "numerical kernels index matrices whose extents are checked on entry"
)]

use crate::kernels::dense;
use crate::kernels::dense::Mat;
use crate::kernels::par::for_chunks_mut;

/// Why a routine of this module failed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum LowRankError {
    /// Inconsistent sizes or non-finite input.
    Shape(String),
    /// The iteration did not reach the tolerance; the string says how far it
    /// got.
    NoConvergence(String),
}

impl std::fmt::Display for LowRankError {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        match self {
            | Self::Shape(m) => write!(f, "shape error: {m}"),
            | Self::NoConvergence(m) => write!(f, "no convergence: {m}"),
        }
    }
}

impl std::error::Error for LowRankError {}

/// A seeded generator of standard normal numbers (xorshift and Box-Muller).
#[derive(Clone, Debug)]
pub struct Gaussian {
    state: u64,
    spare: Option<f64>,
}

impl Gaussian {
    /// A generator; equal seeds give equal streams.
    #[must_use]
    pub const fn new(seed: u64) -> Self {
        // The seed goes through the SplitMix64 finalizer so that nearby seeds
        // (`2k` and `2k + 1` collided under the earlier `seed | 1`) give
        // unrelated streams; xorshift needs a nonzero state.
        let mut z = seed.wrapping_add(0x9E37_79B9_7F4A_7C15);
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        Self {
            state: if z == 0 {
                0x9E37_79B9_7F4A_7C15
            } else {
                z
            },
            spare: None,
        }
    }

    fn uniform(&mut self) -> f64 {
        self.state ^= self.state << 13;
        self.state ^= self.state >> 7;
        self.state ^= self.state << 17;
        ((self.state >> 11) as f64 + 0.5) / (1u64 << 53) as f64
    }

    /// The next standard normal number.
    pub fn sample(&mut self) -> f64 {
        if let Some(z) = self.spare.take() {
            return z;
        }
        let (u, v) = (self.uniform(), self.uniform());
        let r = (-2.0 * u.ln()).sqrt();
        let (s, c) = (std::f64::consts::TAU * v).sin_cos();
        self.spare = Some(r * s);
        r * c
    }
}

fn dot(
    a: &[f64],
    b: &[f64],
) -> f64 {
    super::simd::dot_f64(a, b)
}

fn norm(a: &[f64]) -> f64 {
    super::simd::norm2_sq_f64(a).sqrt()
}

/// Orthonormalises `columns` in place by modified Gram-Schmidt with one
/// re-orthogonalisation pass, against `basis` first. Columns that are
/// numerically dependent are removed. Returns the indices kept.
fn orthonormalize(
    basis: &[Vec<f64>],
    columns: &mut Vec<Vec<f64>>,
) -> Vec<usize> {
    let mut kept: Vec<Vec<f64>> = Vec::with_capacity(columns.len());
    let mut indices = Vec::with_capacity(columns.len());
    for (index, mut v) in std::mem::take(columns).into_iter().enumerate() {
        let before = norm(&v);
        if !before.is_finite() || before == 0.0 {
            continue;
        }
        for _ in 0..2 {
            for q in basis.iter().chain(&kept) {
                let c = dot(q, &v);
                super::simd::axpy_f64(-c, q, &mut v);
            }
        }
        let after = norm(&v);
        // What survives the projection must be more than rounding noise.
        if after > 1e-10 * before {
            super::simd::scale_f64(after.recip(), &mut v);
            kept.push(v);
            indices.push(index);
        }
    }
    *columns = kept;
    indices
}

// ---------------------------------------------------------------------------
// Column-pivoted QR
// ---------------------------------------------------------------------------

const PIVOT_STEP_PARALLEL_WORK: usize = 196_608;
const THIN_Q_PARALLEL_WORK: usize = 4_000_000;

fn householder_grain(
    minimum_work: usize,
    work_per_column: usize,
) -> usize {
    minimum_work.div_ceil(work_per_column).max(4)
}

fn apply_householder(
    column: &mut [f64],
    reflector: &[f64],
) {
    let projection = 2.0 * dot(reflector, column);
    super::simd::axpy_f64(-projection, reflector, column);
}

/// A column-pivoted QR factorisation `A P = Q R`.
#[derive(Clone, Debug)]
pub struct PivotedQr {
    /// `Q`, `m x k` with orthonormal columns, `k = min(m, n)`.
    pub q: Mat,
    /// `R`, `k x n`, upper triangular with a non-increasing diagonal in
    /// magnitude.
    pub r: Mat,
    /// The permutation: column `j` of `A P` is column `perm[j]` of `A`.
    pub perm: Vec<usize>,
}

impl PivotedQr {
    /// The numerical rank: the number of diagonal entries of `R` above
    /// `tol * |R[0][0]|`. `None` picks `max(m, n) * eps`.
    #[must_use]
    pub fn rank(
        &self,
        tol: Option<f64>,
    ) -> usize {
        let tol = tol.unwrap_or_else(|| self.q.rows.max(self.r.cols) as f64 * f64::EPSILON);
        numerical_rank(&self.r, self.r.rows.min(self.r.cols), tol)
    }

    /// The basic solution of the least-squares problem `min |A x - b|`: the
    /// solution with zeros in the `n - rank` columns the pivoting put last.
    ///
    /// # Errors
    /// [`LowRankError::Shape`] when `b` has the wrong length.
    pub fn solve(
        &self,
        b: &[f64],
        tol: Option<f64>,
    ) -> Result<Vec<f64>, LowRankError> {
        if b.len() != self.q.rows {
            return Err(LowRankError::Shape(format!(
                "right-hand side has {} entries, the matrix has {} rows",
                b.len(),
                self.q.rows
            )));
        }
        let rank = self.rank(tol);
        // y = Q^T b, then back-substitution on the leading rank x rank block.
        let mut y: Vec<f64> = (0..rank)
            .map(|j| (0..self.q.rows).map(|i| self.q.at(i, j) * b[i]).sum())
            .collect();
        solve_triangular(&mut y, |i, j| self.r.at(i, j));
        let mut x = vec![0.0; self.r.cols];
        for (j, value) in y.into_iter().enumerate() {
            x[self.perm[j]] = value;
        }
        Ok(x)
    }
}

fn numerical_rank(
    r: &Mat,
    limit: usize,
    tolerance: f64,
) -> usize {
    if limit == 0 {
        return 0;
    }
    let top = r.at(0, 0).abs();
    (0..limit)
        .take_while(|&i| r.at(i, i).abs() > tolerance * top && top > 0.0)
        .count()
}

fn solve_triangular(
    y: &mut [f64],
    entry: impl Fn(usize, usize) -> f64,
) {
    for i in (0..y.len()).rev() {
        for j in i + 1..y.len() {
            y[i] -= entry(i, j) * y[j];
        }
        y[i] /= entry(i, i);
    }
}

struct PivotedReduction {
    transposed_r: Mat,
    perm: Vec<usize>,
    reflectors: Vec<Vec<f64>>,
    scale: f64,
}

/// Column-pivoted Householder QR (Businger and Golub).
///
/// At each step the remaining column of largest norm is brought forward, so
/// the diagonal of `R` decreases in magnitude and its decay reveals the
/// numerical rank. Column norms are downdated and recomputed when
/// cancellation makes the downdate unreliable. Each reflector applies to
/// contiguous columns, with large independent column batches running in parallel.
///
/// # Errors
/// [`LowRankError::Shape`] for inconsistent storage or non-finite entries,
/// and [`LowRankError::NoConvergence`] if the resulting factor exceeds the
/// representable floating-point range.
pub fn qr_pivoted(a: &Mat) -> Result<PivotedQr, LowRankError> {
    let (m, n) = (a.rows, a.cols);
    let k = m.min(n);
    let PivotedReduction {
        transposed_r,
        perm,
        reflectors,
        scale,
    } = reduce_columns(a, k, true)?;
    let mut thin = Mat::zeros(k, n);
    for i in 0..k {
        for j in i..n {
            thin.set(i, j, transposed_r.at(j, i) * scale);
        }
    }
    if !thin.is_finite() {
        return Err(LowRankError::NoConvergence(
            "the triangular factor exceeds finite floating-point range".into(),
        ));
    }
    drop(transposed_r);
    // Thin Q: apply the reflectors to the first k columns of the identity.
    let mut transposed_q = Mat::zeros(k, m);
    if k > 0 {
        for_chunks_mut(
            &mut transposed_q.data,
            m,
            householder_grain(THIN_Q_PARALLEL_WORK, m * k),
            |j, column| {
                column[j] = 1.0;
                for (step, reflector) in reflectors.iter().enumerate().rev() {
                    apply_householder(&mut column[step..], reflector);
                }
            },
        );
    }
    let q = if k <= 1 {
        Mat {
            rows: m,
            cols: k,
            data: transposed_q.data,
        }
    } else {
        transposed_q.transpose()
    };
    Ok(PivotedQr { q, r: thin, perm })
}

fn reduce_columns(
    a: &Mat,
    k: usize,
    retain_reflectors: bool,
) -> Result<PivotedReduction, LowRankError> {
    let (m, n) = (a.rows, a.cols);
    let magnitude = a
        .data
        .iter()
        .try_fold(0.0_f64, |max, &x| x.is_finite().then_some(max.max(x.abs())));
    let Some(magnitude) = magnitude.filter(|_| m.checked_mul(n) == Some(a.data.len())) else {
        return Err(LowRankError::Shape(
            "the matrix storage is inconsistent or has non-finite entries".into(),
        ));
    };
    // Scale only extreme inputs, so the fast sum-of-squares norms operate in
    // a safe range without changing the common path's rounding.
    let scale = if magnitude > 0.0 && !(1e-140..=1e140).contains(&magnitude) {
        magnitude
    } else {
        1.0
    };
    if k == 0 {
        return Ok(PivotedReduction {
            transposed_r: Mat::zeros(n, 0),
            perm: if retain_reflectors {
                (0..n).collect()
            } else {
                Vec::new()
            },
            reflectors: Vec::new(),
            scale,
        });
    }
    // Householder application, swaps and norm fallback all traverse columns.
    let mut transposed_r = a.transpose();
    #[allow(
        clippy::float_cmp,
        reason = "only the exact identity scale skips multiplication"
    )]
    if scale != 1.0 {
        for x in &mut transposed_r.data {
            *x /= scale;
        }
    }
    let mut perm: Vec<usize> = (0..n).collect();
    let mut norms: Vec<f64> = transposed_r.data.chunks_exact(m).map(norm).collect();
    let mut reference = norms.clone();
    let mut reflectors = Vec::with_capacity(if retain_reflectors { k } else { 0 });
    let mut reused_reflector = Vec::new();
    for step in 0..k {
        let pivot = (step..n)
            .max_by(|&x, &y| norms[x].total_cmp(&norms[y]))
            .unwrap_or(step);
        if pivot != step {
            let (head, tail) = transposed_r.data.split_at_mut(pivot * m);
            head[step * m..(step + 1) * m].swap_with_slice(&mut tail[..m]);
            perm.swap(step, pivot);
            norms.swap(step, pivot);
            reference.swap(step, pivot);
        }
        // Householder vector for column `step`, rows `step..m`.
        let start = step * m;
        let mut v = std::mem::take(&mut reused_reflector);
        v.extend_from_slice(&transposed_r.data[start + step..start + m]);
        let alpha = norm(&v);
        if alpha > 0.0 {
            let beta = if v[0] >= 0.0 { -alpha } else { alpha };
            v[0] -= beta;
            let vn = norm(&v);
            if vn > 0.0 {
                for x in &mut v {
                    *x /= vn;
                }
                for_chunks_mut(
                    &mut transposed_r.data[start..],
                    m,
                    householder_grain(PIVOT_STEP_PARALLEL_WORK, m - step),
                    |_, column| apply_householder(&mut column[step..], &v),
                );
            }
        } else {
            v.fill(0.0);
        }
        transposed_r.data[start + step + 1..start + m].fill(0.0);
        // Downdate the norms of the remaining columns.
        for j in step + 1..n {
            if norms[j] > 0.0 {
                let ratio = transposed_r.at(j, step) / norms[j];
                let remaining = (1.0 - ratio * ratio).max(0.0);
                norms[j] *= remaining.sqrt();
                if norms[j] < 1e-7 * reference[j] {
                    norms[j] = norm(&transposed_r.data[j * m + step + 1..(j + 1) * m]);
                    reference[j] = norms[j];
                }
            }
        }
        if retain_reflectors {
            reflectors.push(v);
        } else {
            v.clear();
            reused_reflector = v;
        }
    }
    Ok(PivotedReduction {
        transposed_r,
        perm,
        reflectors,
        scale,
    })
}

/// An approximation `A ≈ A[:, columns] * coefficients`.
#[derive(Clone, Debug)]
pub struct ColumnInterpolation {
    /// Indices of the original columns used as the basis, in pivot order.
    pub columns: Vec<usize>,
    /// Dimensionless interpolation weights, one row per selected column.
    /// Columns appear in the same order as in the input matrix.
    pub coefficients: Mat,
}

/// Computes a column interpolative decomposition by truncated pivoted QR.
///
/// The actual rank is at most `max_rank`. The leading triangular block
/// determines the interpolation coefficients; selected columns have identity
/// coefficients. A zero-rank result has a zero-row coefficient matrix.
///
/// # Errors
///
/// Returns an error when the input has inconsistent storage or non-finite
/// entries, the requested rank exceeds either extent, or the coefficients
/// cannot be represented as finite floats.
pub fn interpolative_decomposition(
    a: &Mat,
    max_rank: usize,
) -> Result<ColumnInterpolation, LowRankError> {
    if max_rank > a.rows.min(a.cols) {
        return Err(LowRankError::Shape(
            "interpolative rank exceeds a matrix extent".into(),
        ));
    }
    let PivotedReduction {
        transposed_r,
        mut perm,
        ..
    } = reduce_columns(a, max_rank, false)?;
    let tolerance = a.rows.max(a.cols) as f64 * f64::EPSILON;
    let rank = numerical_rank(&transposed_r, max_rank, tolerance);
    let mut coefficients = Mat::zeros(rank, a.cols);
    let mut rhs = vec![0.0; rank];
    for (pivot, &column) in perm.iter().enumerate() {
        if pivot < rank {
            coefficients.set(pivot, column, 1.0);
            continue;
        }
        for (i, value) in rhs.iter_mut().enumerate() {
            *value = transposed_r.at(pivot, i);
        }
        solve_triangular(&mut rhs, |i, j| transposed_r.at(j, i));
        for (i, &value) in rhs.iter().enumerate() {
            coefficients.set(i, column, value);
        }
    }
    if !coefficients.is_finite() {
        return Err(LowRankError::NoConvergence(
            "interpolation coefficients exceed finite floating-point range".into(),
        ));
    }
    perm.truncate(rank);
    Ok(ColumnInterpolation {
        columns: perm,
        coefficients,
    })
}

// ---------------------------------------------------------------------------
// Randomized SVD
// ---------------------------------------------------------------------------

/// Parameters of [`randomized_svd`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RandomizedOptions {
    /// Extra random directions beyond the requested rank (5 to 10 is the
    /// usual range).
    pub oversample: usize,
    /// Power iterations; each sharpens the spectrum at the price of two more
    /// passes over the operator. Two suffice unless the singular values decay
    /// very slowly.
    pub power_iterations: usize,
    /// Seed of the random test matrix.
    pub seed: u64,
}

impl Default for RandomizedOptions {
    fn default() -> Self {
        Self {
            oversample: 8,
            power_iterations: 2,
            seed: 0x5EED_1234_ABCD_0001,
        }
    }
}

/// A rank-`k` approximate singular value decomposition of an `m x n`
/// operator given by its products: `matvec(x, y)` writes `A x` (`x` of
/// length `n`, `y` of length `m`) and `rmatvec(y, x)` writes `A^T y`.
///
/// The range of `A` is sampled with `k + oversample` random vectors,
/// sharpened by power iterations with re-orthogonalisation between them, and
/// the small projected matrix is decomposed densely. The cost is
/// at most `2 (power_iterations + 1) (k + oversample)` products. Extreme
/// sampled columns are rescaled before orthogonalisation. Singular values
/// retain the scale of the original operator.
///
/// # Errors
/// [`LowRankError::Shape`] when `k` is zero or exceeds `min(m, n)`, or a
/// product returns non-finite numbers; [`LowRankError::NoConvergence`] when
/// the dense decomposition of the projected matrix fails.
pub fn randomized_svd(
    matvec: &dyn Fn(&[f64], &mut [f64]),
    rmatvec: &dyn Fn(&[f64], &mut [f64]),
    m: usize,
    n: usize,
    k: usize,
    options: &RandomizedOptions,
) -> Result<dense::Svd, LowRankError> {
    if k == 0 || k > m.min(n) {
        return Err(LowRankError::Shape(format!(
            "rank {k} requested of a {m} x {n} operator"
        )));
    }
    let l = k.saturating_add(options.oversample).min(m.min(n));
    let mut rng = Gaussian::new(options.seed);
    let write_product = |f: &dyn Fn(&[f64], &mut [f64]), x: &[f64], y: &mut [f64]| {
        f(x, y);
        y.iter()
            .try_fold(0.0_f64, |peak, &v| {
                v.is_finite().then_some(peak.max(v.abs()))
            })
            .ok_or_else(|| LowRankError::Shape("the operator returned a non-finite product".into()))
    };
    let apply =
        |f: &dyn Fn(&[f64], &mut [f64]), x: &[f64], len: usize| -> Result<Vec<f64>, LowRankError> {
            let mut y = vec![0.0; len];
            let peak = write_product(f, x, &mut y)?;
            if peak > 0.0 && !(1e-140..=1e140).contains(&peak) {
                // Rescaling a column preserves the sampled subspace and
                // keeps Gram-Schmidt's squared norms in range.
                for value in &mut y {
                    *value /= peak;
                }
            }
            Ok(y)
        };
    let mut omega = vec![0.0; n];
    let mut q = Vec::with_capacity(l);
    for _ in 0..l {
        omega.fill_with(|| rng.sample());
        q.push(apply(matvec, &omega, m)?);
    }
    orthonormalize(&[], &mut q);
    for _ in 0..options.power_iterations {
        let mut z = q
            .iter()
            .map(|x| apply(rmatvec, x, n))
            .collect::<Result<Vec<_>, _>>()?;
        orthonormalize(&[], &mut z);
        q = z
            .iter()
            .map(|x| apply(matvec, x, m))
            .collect::<Result<Vec<_>, _>>()?;
        orthonormalize(&[], &mut q);
    }
    if q.is_empty() {
        return Err(LowRankError::Shape(
            "the operator returned a numerically zero sampled range".into(),
        ));
    }
    // B = Q^T A, row i being (A^T q_i)^T.
    let mut b = Mat::zeros(q.len(), n);
    for (row, x) in b.data.chunks_exact_mut(n).zip(&q) {
        write_product(rmatvec, x, row)?;
    }
    let small = dense::svd(&b).map_err(|e| LowRankError::NoConvergence(format!("{e:?}")))?;
    let k = k.min(small.s.len());
    let mut u = Mat::zeros(m, k);
    // Lift U = Q * U_small with contiguous output rows. Each element
    // accumulates its contributions from Q in the original order.
    for (qi, coefficients) in q.iter().zip(small.u.data.chunks_exact(small.u.cols)) {
        for (row, &sample) in u.data.chunks_exact_mut(k).zip(qi) {
            for (value, &coefficient) in row.iter_mut().zip(coefficients) {
                *value += coefficient * sample;
            }
        }
    }
    let v = if small.v.cols == k {
        small.v
    } else {
        let mut v = Mat::zeros(n, k);
        for (out, row) in v
            .data
            .chunks_exact_mut(k)
            .zip(small.v.data.chunks_exact(small.v.cols))
        {
            out.copy_from_slice(&row[..k]);
        }
        v
    };
    Ok(dense::Svd {
        u,
        s: small.s.into_iter().take(k).collect(),
        v,
    })
}

// ---------------------------------------------------------------------------
// LOBPCG
// ---------------------------------------------------------------------------

/// Eigenpairs found by [`lobpcg`].
#[derive(Clone, Debug)]
pub struct Eigenpairs {
    /// Eigenvalues, ascending.
    pub values: Vec<f64>,
    /// `vectors[i]` is the unit eigenvector of `values[i]`.
    pub vectors: Vec<Vec<f64>>,
    /// Euclidean residual norms of the eigenvalue equation.
    pub residuals: Vec<f64>,
    /// Iterations used.
    pub iterations: usize,
}

/// An optional preconditioner `r -> T r`, applied as `(input, output)`.
type Preconditioner<'a> = Option<&'a dyn Fn(&[f64], &mut [f64])>;

/// The `k` smallest eigenpairs of a symmetric operator of size `n` by LOBPCG.
///
/// `apply(x, y)` writes `A x`. `precondition`, when given, writes an
/// approximation of `A^-1 r` (it must be symmetric positive definite; the
/// identity is used otherwise). Each iteration performs a Rayleigh-Ritz
/// projection on the span of the current vectors, the preconditioned
/// residuals and the previous search directions, with the basis
/// re-orthonormalised every time, which is what keeps the method stable when
/// the directions become nearly dependent close to convergence.
///
/// # Errors
/// [`LowRankError::Shape`] when `k` is zero or more than `n / 3` would allow
/// (`3 k <= n` is needed for the trial subspace), or a product is not
/// finite, or when the tolerance is not positive and finite.
/// Returns [`LowRankError::NoConvergence`] when `max_iterations` pass with a
/// residual above `tol * max(1, |lambda|)`.
pub fn lobpcg(
    apply: &dyn Fn(&[f64], &mut [f64]),
    precondition: Preconditioner<'_>,
    n: usize,
    k: usize,
    tol: f64,
    max_iterations: usize,
    seed: u64,
) -> Result<Eigenpairs, LowRankError> {
    if k == 0 || k.saturating_mul(3) > n {
        return Err(LowRankError::Shape(format!(
            "{k} eigenpairs requested of an operator of size {n}; LOBPCG needs 3 k <= n"
        )));
    }
    validate_lobpcg_workspace(n, k, tol)?;
    let mul = |x: &[f64]| -> Vec<f64> {
        let mut y = vec![0.0; n];
        apply(x, &mut y);
        y
    };
    let mut rng = Gaussian::new(seed);
    let mut x: Vec<Vec<f64>> = (0..k)
        .map(|_| (0..n).map(|_| rng.sample()).collect())
        .collect();
    orthonormalize(&[], &mut x);
    let mut p: Vec<Vec<f64>> = Vec::new();
    let mut values = vec![0.0; k];
    let mut residuals = vec![f64::INFINITY; k];
    for iteration in 0..=max_iterations {
        // Rayleigh-Ritz on the current block alone gives values and residuals.
        let ax: Vec<Vec<f64>> = x.iter().map(|v| mul(v)).collect();
        if ax.iter().flatten().any(|v| !v.is_finite()) {
            return Err(LowRankError::Shape(
                "the operator returned non-finite products".into(),
            ));
        }
        let mut w: Vec<Vec<f64>> = Vec::with_capacity(k);
        for i in 0..x.len() {
            values[i] = dot(&x[i], &ax[i]);
            let r: Vec<f64> = ax[i]
                .iter()
                .zip(&x[i])
                .map(|(a, v)| a - values[i] * v)
                .collect();
            residuals[i] = norm(&r);
            let mut t = vec![0.0; n];
            match precondition {
                | Some(f) => f(&r, &mut t),
                | None => t.copy_from_slice(&r),
            }
            if t.iter().any(|value| !value.is_finite()) {
                return Err(LowRankError::Shape(
                    "non-finite LOBPCG preconditioner product".into(),
                ));
            }
            w.push(t);
        }
        let converged = values
            .iter()
            .zip(&residuals)
            .all(|(l, r)| *r <= tol * l.abs().max(1.0));
        if converged || iteration == max_iterations {
            let mut order: Vec<usize> = (0..x.len()).collect();
            order.sort_by(|&a, &b| values[a].total_cmp(&values[b]));
            let result = Eigenpairs {
                values: order.iter().map(|&i| values[i]).collect(),
                vectors: order.iter().map(|&i| x[i].clone()).collect(),
                residuals: order.iter().map(|&i| residuals[i]).collect(),
                iterations: iteration,
            };
            return if converged {
                Ok(result)
            } else {
                Err(LowRankError::NoConvergence(format!(
                    "largest residual {:e} after {max_iterations} iterations",
                    result.residuals.iter().copied().fold(0.0, f64::max)
                )))
            };
        }
        // Trial basis [X, W, P], orthonormal.
        orthonormalize(&x, &mut w);
        let mut basis = x.clone();
        basis.extend(w);
        let mut directions = std::mem::take(&mut p);
        orthonormalize(&basis, &mut directions);
        basis.extend(directions);
        let size = basis.len();
        let a_basis: Vec<Vec<f64>> = basis.iter().map(|v| mul(v)).collect();
        let mut h = Mat::zeros(size, size);
        for i in 0..size {
            for j in i..size {
                // Symmetrised, so rounding in the products cannot make the
                // projected matrix unsymmetric.
                let value = f64::midpoint(dot(&basis[i], &a_basis[j]), dot(&basis[j], &a_basis[i]));
                h.set(i, j, value);
                h.set(j, i, value);
            }
        }
        let (_, vectors) = dense::eigen_symmetric(&h)
            .map_err(|e| LowRankError::NoConvergence(format!("projected problem: {e:?}")))?;
        let kept = x.len();
        let mut next: Vec<Vec<f64>> = Vec::with_capacity(kept);
        for column in 0..kept {
            let mut v = vec![0.0; n];
            for (row, b) in basis.iter().enumerate() {
                let c = vectors.at(row, column);
                for (out, value) in v.iter_mut().zip(b) {
                    *out += c * value;
                }
            }
            next.push(v);
        }
        // New search directions: what the step added beyond the old block.
        let mut steps: Vec<Vec<f64>> = next.clone();
        orthonormalize(&x, &mut steps);
        p = steps;
        orthonormalize(&[], &mut next);
        x = next;
    }
    Err(LowRankError::NoConvergence(
        "iteration budget exhausted".into(),
    ))
}

/// Controls for the generalized symmetric-definite LOBPCG solver.
#[derive(Clone, Copy, Debug)]
pub struct LobpcgOptions {
    /// Positive finite tolerance on `||A x - lambda B x||₂ / max(1, |lambda|)`.
    pub tolerance: f64,
    /// Maximum number of trial-space expansions after the initial projection.
    pub max_iterations: usize,
    /// Seed for the starting block.
    pub seed: u64,
    /// Omit converged approximations from residual and search-direction updates.
    ///
    /// All approximations participate in Rayleigh–Ritz. A converged pair can
    /// resume its updates if a later projection increases its residual.
    pub soft_lock: bool,
}

impl Default for LobpcgOptions {
    fn default() -> Self {
        Self {
            tolerance: 1e-8,
            max_iterations: 200,
            seed: 0x5eed,
            soft_lock: true,
        }
    }
}

struct GeneralizedTrial {
    x: Vec<f64>,
    bx: Vec<f64>,
    ax: Vec<f64>,
}

impl GeneralizedTrial {
    fn linear_combination(
        basis: &[Self],
        coefficients: impl Iterator<Item = f64>,
    ) -> Self {
        let n = basis.first().map_or(0, |column| column.x.len());
        let mut result = Self {
            x: vec![0.0; n],
            bx: vec![0.0; n],
            ax: vec![0.0; n],
        };
        for (column, coefficient) in basis.iter().zip(coefficients) {
            axpy(&mut result.x, coefficient, &column.x);
            axpy(&mut result.bx, coefficient, &column.bx);
            axpy(&mut result.ax, coefficient, &column.ax);
        }
        result
    }
}

fn axpy(
    target: &mut [f64],
    coefficient: f64,
    source: &[f64],
) {
    for (out, value) in target.iter_mut().zip(source) {
        *out += coefficient * value;
    }
}

fn generalized_candidate(
    mut x: Vec<f64>,
    apply: &dyn Fn(&[f64], &mut [f64]),
    mass: &dyn Fn(&[f64], &mut [f64]),
    basis: &[GeneralizedTrial],
) -> Result<Option<GeneralizedTrial>, LowRankError> {
    if x.iter().any(|value| !value.is_finite()) {
        return Err(LowRankError::Shape("non-finite LOBPCG trial vector".into()));
    }
    if x.iter().all(|&value| value == 0.0) {
        return Ok(None);
    }
    let mut bx = vec![0.0; x.len()];
    mass(&x, &mut bx);
    let initial_squared = dot(&x, &bx);
    if !initial_squared.is_finite() || initial_squared <= 0.0 {
        return Err(LowRankError::Shape(
            "mass products must be finite and positive definite on trial vectors".into(),
        ));
    }
    for _ in 0..2 {
        for q in basis {
            let coefficient = dot(&q.x, &bx);
            axpy(&mut x, -coefficient, &q.x);
            axpy(&mut bx, -coefficient, &q.bx);
        }
    }
    let squared = dot(&x, &bx);
    if !squared.is_finite() {
        return Err(LowRankError::Shape("non-finite mass inner product".into()));
    }
    if squared <= initial_squared * 1e-20 {
        // A direction in the current subspace contributes no new information.
        return Ok(None);
    }
    let scale = squared.sqrt();
    for value in x.iter_mut().chain(&mut bx) {
        *value /= scale;
    }
    let mut ax = vec![0.0; x.len()];
    apply(&x, &mut ax);
    if x.iter()
        .chain(&bx)
        .chain(&ax)
        .any(|value| !value.is_finite())
    {
        return Err(LowRankError::Shape(
            "non-finite LOBPCG operator product".into(),
        ));
    }
    Ok(Some(GeneralizedTrial { x, bx, ax }))
}

fn generalized_ritz(
    basis: &[GeneralizedTrial],
    k: usize,
) -> Result<(Vec<f64>, Vec<GeneralizedTrial>), LowRankError> {
    let size = basis.len();
    let mut projected = Mat::zeros(size, size);
    for (i, left) in basis.iter().enumerate() {
        for (j, right) in basis.iter().enumerate().skip(i) {
            let entry = f64::midpoint(dot(&left.x, &right.ax), dot(&right.x, &left.ax));
            projected.set(i, j, entry);
            projected.set(j, i, entry);
        }
    }
    let (mut values, coefficients) = dense::eigen_symmetric(&projected).map_err(|error| {
        LowRankError::NoConvergence(format!("projected eigenproblem: {error:?}"))
    })?;
    values.truncate(k);
    let vectors = (0..k)
        .map(|column| {
            GeneralizedTrial::linear_combination(
                basis,
                (0..size).map(|row| coefficients.at(row, column)),
            )
        })
        .collect();
    Ok((values, vectors))
}

fn scaled_residual_norm(residual: &[f64]) -> f64 {
    let scale = residual.iter().map(|value| value.abs()).fold(0.0, f64::max);
    if scale == 0.0 || !scale.is_finite() {
        return scale;
    }
    scale
        * residual
            .iter()
            .map(|value| (value / scale).powi(2))
            .sum::<f64>()
            .sqrt()
}

/// The `k` smallest eigenpairs of a symmetric-definite pencil `A x = lambda B x`.
///
/// Both operators use callbacks of the form `(input, output)`. `apply` writes
/// `A x`, and `mass` writes `B x` for symmetric positive-definite `B`.
/// `precondition` optionally applies a symmetric positive-definite approximation
/// to the inverse of the shifted operator. The solver returns B-orthonormal
/// vectors in ascending eigenvalue order. It forms neither full matrix.
///
/// # Errors
/// Rejects invalid dimensions, a nonpositive or non-finite tolerance,
/// non-finite products, and mass products that are nonpositive on trial
/// vectors. Reports nonconvergence when the work budget is exhausted or
/// the trial subspace cannot advance.
pub fn lobpcg_generalized(
    apply: &dyn Fn(&[f64], &mut [f64]),
    mass: &dyn Fn(&[f64], &mut [f64]),
    precondition: Preconditioner<'_>,
    n: usize,
    k: usize,
    options: LobpcgOptions,
) -> Result<Eigenpairs, LowRankError> {
    validate_lobpcg_workspace(n, k, options.tolerance)?;

    let mut rng = Gaussian::new(options.seed);
    let mut basis = Vec::with_capacity(k);
    for _ in 0..k {
        let mut vector = Vec::new();
        vector
            .try_reserve_exact(n)
            .map_err(|_| LowRankError::Shape("LOBPCG cannot allocate a trial vector".into()))?;
        vector.extend((0..n).map(|_| rng.sample()));
        let q = generalized_candidate(vector, apply, mass, &basis)?.ok_or_else(|| {
            LowRankError::NoConvergence("the initial block is numerically rank deficient".into())
        })?;
        basis.push(q);
    }
    let (mut values, mut vectors) = generalized_ritz(&basis, k)?;
    let mut directions: Vec<Option<Vec<f64>>> = vec![None; k];

    for iteration in 0..=options.max_iterations {
        crate::cancellation::checkpoint()
            .map_err(|error| LowRankError::NoConvergence(error.to_string()))?;
        let (mut residuals, mut norms) = generalized_residuals(&vectors, &values)?;
        let mut converged = convergence_flags(&norms, &values, options.tolerance);
        if converged.iter().all(|&done| done) {
            // Ritz vector products are assembled from cached trial images.
            // Reapply the original operators to certify the reported residuals.
            for column in &mut vectors {
                column.ax.fill(0.0);
                column.bx.fill(0.0);
                apply(&column.x, &mut column.ax);
                mass(&column.x, &mut column.bx);
            }
            (residuals, norms) = generalized_residuals(&vectors, &values)?;
            converged = convergence_flags(&norms, &values, options.tolerance);
            if converged.iter().all(|&done| done) {
                return Ok(Eigenpairs {
                    values,
                    vectors: vectors.into_iter().map(|column| column.x).collect(),
                    residuals: norms,
                    iterations: iteration,
                });
            }
        }
        if iteration == options.max_iterations {
            return Err(LowRankError::NoConvergence(format!(
                "largest residual {:e} after {iteration} iterations",
                norms.iter().copied().fold(0.0, f64::max)
            )));
        }

        // Converged approximations stay in Rayleigh–Ritz, while soft locking
        // keeps their residuals out of the preconditioner and search block.
        let active = |i: usize| !options.soft_lock || !converged[i];
        basis = vectors;
        for (_, rate) in residuals
            .into_iter()
            .enumerate()
            .filter(|(i, _)| active(*i))
        {
            if basis.len() == n {
                break;
            }
            let search = if let Some(precondition) = precondition {
                let mut output = vec![0.0; n];
                precondition(&rate, &mut output);
                output
            } else {
                rate
            };
            if let Some(column) = generalized_candidate(search, apply, mass, &basis)? {
                basis.push(column);
            }
        }
        for (index, slot) in directions.iter_mut().enumerate() {
            if basis.len() == n {
                break;
            }
            if active(index)
                && let Some(search) = slot.take()
                && let Some(column) = generalized_candidate(search, apply, mass, &basis)?
            {
                basis.push(column);
            }
        }
        if basis.len() == k {
            return Err(LowRankError::NoConvergence(
                "no independent LOBPCG search direction".into(),
            ));
        }

        (values, vectors) = generalized_ritz(&basis, k)?;
        for (index, (slot, next)) in directions.iter_mut().zip(&vectors).enumerate() {
            *slot = active(index).then(|| {
                let mut step = next.x.clone();
                for previous in basis.iter().take(k) {
                    let coefficient = dot(&previous.x, &next.bx);
                    axpy(&mut step, -coefficient, &previous.x);
                }
                step
            });
        }
    }
    Err(LowRankError::NoConvergence(
        "iteration budget exhausted".into(),
    ))
}

fn validate_lobpcg_workspace(
    n: usize,
    k: usize,
    tolerance: f64,
) -> Result<(), LowRankError> {
    let error = || LowRankError::Shape("invalid or unrepresentable LOBPCG workspace".into());
    if k == 0 || k > n || tolerance <= 0.0 || !tolerance.is_finite() {
        return Err(error());
    }
    let trial_width = k.saturating_mul(3).min(n);
    let elements = n
        .checked_mul(trial_width)
        .and_then(|basis| {
            trial_width
                .checked_mul(trial_width)
                .and_then(|gram| gram.checked_add(basis))
        })
        .ok_or_else(error)?;
    if elements
        .checked_mul(std::mem::size_of::<f64>())
        .is_none_or(|bytes| isize::try_from(bytes).is_err())
    {
        return Err(error());
    }
    Ok(())
}

fn convergence_flags(
    norms: &[f64],
    values: &[f64],
    tolerance: f64,
) -> Vec<bool> {
    norms
        .iter()
        .zip(values)
        .map(|(&residual, &lambda)| residual <= tolerance * lambda.abs().max(1.0))
        .collect()
}

fn generalized_residuals(
    vectors: &[GeneralizedTrial],
    values: &[f64],
) -> Result<(Vec<Vec<f64>>, Vec<f64>), LowRankError> {
    let mut residuals = Vec::with_capacity(values.len());
    let mut norms = Vec::with_capacity(values.len());
    for (column, &lambda) in vectors.iter().zip(values) {
        let residual: Vec<_> = column
            .ax
            .iter()
            .zip(&column.bx)
            .map(|(&ax, &bx)| ax - lambda * bx)
            .collect();
        if residual.iter().any(|value| !value.is_finite()) {
            return Err(LowRankError::NoConvergence(
                "non-finite generalized eigenpair residual".into(),
            ));
        }
        norms.push(scaled_residual_norm(&residual));
        residuals.push(residual);
    }
    Ok((residuals, norms))
}

// ---------------------------------------------------------------------------
// Anderson acceleration
// ---------------------------------------------------------------------------

/// The outcome of [`anderson`].
#[derive(Clone, Debug, PartialEq)]
pub struct FixedPoint {
    /// The point found.
    pub x: Vec<f64>,
    /// `|g(x) - x|` there.
    pub residual: f64,
    /// Evaluations of `g`.
    pub evaluations: usize,
}

/// Solves `x = g(x)` by Anderson acceleration with a window of `window`
/// previous iterates.
///
/// Each step mixes the last iterates so that the linearised residual is
/// smallest in the least-squares sense; with `window = 0` it is the plain
/// iteration `x <- g(x)`. `damping` in `(0, 1]` scales the step (`1` is the
/// undamped method). A step whose least-squares problem is degenerate, or
/// that produces non-finite numbers, falls back to the plain iterate and
/// restarts the window.
///
/// # Errors
/// [`LowRankError::Shape`] when `g` changes the length of the vector or the
/// start is not finite; [`LowRankError::NoConvergence`] when the residual is
/// still above `tol` after `max_evaluations`.
pub fn anderson(
    g: &dyn Fn(&[f64]) -> Vec<f64>,
    start: &[f64],
    window: usize,
    damping: f64,
    tol: f64,
    max_evaluations: usize,
) -> Result<FixedPoint, LowRankError> {
    let n = start.len();
    if start.iter().any(|v| !v.is_finite()) {
        return Err(LowRankError::Shape("the start is not finite".into()));
    }
    let damping = if damping > 0.0 && damping <= 1.0 {
        damping
    } else {
        1.0
    };
    let mut x = start.to_vec();
    // Differences of residuals and of g-values over the window.
    let mut d_f: Vec<Vec<f64>> = Vec::new();
    let mut d_g: Vec<Vec<f64>> = Vec::new();
    let mut previous: Option<(Vec<f64>, Vec<f64>)> = None;
    for evaluation in 1..=max_evaluations.max(1) {
        let gx = g(&x);
        if gx.len() != n {
            return Err(LowRankError::Shape(format!(
                "g returned {} components for {n}",
                gx.len()
            )));
        }
        let f: Vec<f64> = gx.iter().zip(&x).map(|(a, b)| a - b).collect();
        let residual = norm(&f);
        if !residual.is_finite() {
            return Err(LowRankError::NoConvergence(
                "the iteration left the finite numbers".into(),
            ));
        }
        if residual <= tol {
            return Ok(FixedPoint {
                x,
                residual,
                evaluations: evaluation,
            });
        }
        if let Some((g_old, f_old)) = &previous {
            d_g.push(gx.iter().zip(g_old).map(|(a, b)| a - b).collect());
            d_f.push(f.iter().zip(f_old).map(|(a, b)| a - b).collect());
            if d_f.len() > window {
                d_f.remove(0);
                d_g.remove(0);
            }
        }
        // gamma minimises |f - dF gamma|; x_next = g(x) - dG gamma.
        let mut next = gx.clone();
        if !d_f.is_empty() {
            match least_squares(&d_f, &f) {
                | Some(gamma) => {
                    for (column, c) in d_g.iter().zip(&gamma) {
                        for (out, v) in next.iter_mut().zip(column) {
                            *out -= c * v;
                        }
                    }
                    if damping < 1.0 {
                        // Mix with the damped plain step on the same combination.
                        for ((out, xi), fi) in next.iter_mut().zip(&x).zip(&f) {
                            let plain = xi + damping * fi;
                            *out = damping * *out + (1.0 - damping) * plain;
                        }
                    }
                    if next.iter().any(|v| !v.is_finite()) {
                        next.clone_from(&gx);
                        d_f.clear();
                        d_g.clear();
                    }
                },
                | None => {
                    d_f.clear();
                    d_g.clear();
                },
            }
        } else if damping < 1.0 {
            for ((out, xi), fi) in next.iter_mut().zip(&x).zip(&f) {
                *out = xi + damping * fi;
            }
        }
        previous = Some((gx, f));
        x = next;
    }
    let gx = g(&x);
    let residual = norm(&gx.iter().zip(&x).map(|(a, b)| a - b).collect::<Vec<_>>());
    Err(LowRankError::NoConvergence(format!(
        "residual {residual:e} after {max_evaluations} evaluations"
    )))
}

/// The least-squares solution of `sum c_j columns[j] = rhs` by Gram-Schmidt
/// QR, or `None` when the columns are numerically dependent.
fn least_squares(
    columns: &[Vec<f64>],
    rhs: &[f64],
) -> Option<Vec<f64>> {
    let k = columns.len();
    let mut q: Vec<Vec<f64>> = Vec::with_capacity(k);
    let mut r = vec![vec![0.0; k]; k];
    for (j, column) in columns.iter().enumerate() {
        let mut v = column.clone();
        let before = norm(&v);
        for _ in 0..2 {
            for (i, qi) in q.iter().enumerate() {
                let c = dot(qi, &v);
                r[i][j] += c;
                for (x, y) in v.iter_mut().zip(qi) {
                    *x -= c * y;
                }
            }
        }
        let after = norm(&v);
        if after.partial_cmp(&(1e-12 * before)) != Some(std::cmp::Ordering::Greater) {
            return None;
        }
        r[j][j] = after;
        q.push(v.into_iter().map(|x| x / after).collect());
    }
    let mut c: Vec<f64> = q.iter().map(|qi| dot(qi, rhs)).collect();
    for i in (0..k).rev() {
        for j in i + 1..k {
            c[i] -= r[i][j] * c[j];
        }
        c[i] /= r[i][i];
    }
    c.iter().all(|v| v.is_finite()).then_some(c)
}

#[cfg(test)]
mod tests {
    #![allow(
        clippy::unwrap_used,
        clippy::expect_used,
        clippy::panic,
        clippy::float_cmp
    )]

    use super::*;

    fn mat(rows: &[&[f64]]) -> Mat {
        Mat::from_rows(&rows.iter().map(|r| r.to_vec()).collect::<Vec<_>>()).unwrap()
    }

    /// `A = sum_i s_i u_i v_i^T` with orthonormal factors from a seed.
    fn with_spectrum(
        m: usize,
        n: usize,
        s: &[f64],
        seed: u64,
    ) -> Mat {
        let mut rng = Gaussian::new(seed);
        let mut u: Vec<Vec<f64>> = (0..s.len())
            .map(|_| (0..m).map(|_| rng.sample()).collect())
            .collect();
        let mut v: Vec<Vec<f64>> = (0..s.len())
            .map(|_| (0..n).map(|_| rng.sample()).collect())
            .collect();
        orthonormalize(&[], &mut u);
        orthonormalize(&[], &mut v);
        let mut a = Mat::zeros(m, n);
        for ((ui, vi), si) in u.iter().zip(&v).zip(s) {
            for i in 0..m {
                for j in 0..n {
                    let value = a.at(i, j) + si * ui[i] * vi[j];
                    a.set(i, j, value);
                }
            }
        }
        a
    }

    #[test]
    fn the_gaussian_generator_is_seeded_and_standard() {
        let mut a = Gaussian::new(7);
        let mut b = Gaussian::new(7);
        let xs: Vec<f64> = (0..20_000).map(|_| a.sample()).collect();
        assert!((0..100).all(|i| xs[i] == b.sample()));
        let mean = xs.iter().sum::<f64>() / xs.len() as f64;
        let var = xs.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / xs.len() as f64;
        assert!(
            mean.abs() < 0.03 && (var - 1.0).abs() < 0.05,
            "{mean} {var}"
        );
    }

    #[test]
    fn neighbouring_seeds_give_different_streams() {
        // `seed | 1` made 2k and 2k + 1 the same generator.
        for seed in [0_u64, 2, 42, 1_000_000] {
            let mut a = Gaussian::new(seed);
            let mut b = Gaussian::new(seed + 1);
            let same = (0..32).filter(|_| a.sample() == b.sample()).count();
            assert_eq!(same, 0, "seeds {seed} and {}", seed + 1);
        }
        // Seed zero is a valid, non-degenerate stream.
        let mut z = Gaussian::new(0);
        assert!((0..8).any(|_| z.sample() != 0.0));
    }

    #[test]
    fn pivoted_qr_reconstructs_and_reveals_rank() {
        let a = with_spectrum(9, 6, &[5.0, 2.0, 0.5], 3);
        let f = qr_pivoted(&a).unwrap();
        assert_eq!(f.rank(None), 3);
        // A P = Q R.
        let qr = f.q.matmul(&f.r).unwrap();
        for i in 0..9 {
            for j in 0..6 {
                assert!((qr.at(i, j) - a.at(i, f.perm[j])).abs() < 1e-12);
            }
        }
        // Orthonormal columns and a decreasing diagonal.
        let gram = f.q.transpose().matmul(&f.q).unwrap();
        for i in 0..6 {
            for j in 0..6 {
                let want = f64::from(u8::from(i == j));
                assert!((gram.at(i, j) - want).abs() < 1e-12);
            }
        }
        assert!((1..6).all(|i| f.r.at(i, i).abs() <= f.r.at(i - 1, i - 1).abs() + 1e-12));
        // A consistent rank-deficient system is solved exactly.
        let x_true = [1.0, -2.0, 0.5, 3.0, 0.0, -1.0];
        let b = a.matvec(&x_true);
        let x = f.solve(&b, None).unwrap();
        let back = a.matvec(&x);
        assert!(back.iter().zip(&b).all(|(p, q)| (p - q).abs() < 1e-11));
        assert_eq!(x.iter().filter(|v| **v == 0.0).count(), 3);
        assert!(f.solve(&[1.0], None).is_err());
    }

    #[test]
    fn pivoted_qr_edge_cases() {
        let zero = Mat::zeros(3, 2);
        let f = qr_pivoted(&zero).unwrap();
        assert_eq!(f.rank(None), 0);
        assert_eq!(f.solve(&[1.0, 2.0, 3.0], None).unwrap(), vec![0.0, 0.0]);
        let wide = mat(&[&[1.0, 2.0, 3.0], &[2.0, 4.0, 6.0]]);
        assert_eq!(qr_pivoted(&wide).unwrap().rank(None), 1);
        let empty = Mat::zeros(0, 0);
        assert_eq!(qr_pivoted(&empty).unwrap().rank(None), 0);
        let mut bad = Mat::zeros(2, 2);
        bad.set(0, 0, f64::NAN);
        assert!(qr_pivoted(&bad).is_err());
    }

    #[test]
    fn pivoted_qr_preserves_extreme_uniform_scales() {
        for scale in [1e-250, 1e250] {
            let a = mat(&[&[scale, 2.0 * scale], &[2.0 * scale, -scale]]);
            let factor = qr_pivoted(&a).unwrap();
            assert_eq!(factor.rank(None), 2);
            let product = factor.q.matmul(&factor.r).unwrap();
            for i in 0..a.rows {
                for j in 0..a.cols {
                    let want = a.at(i, factor.perm[j]) / scale;
                    assert!((product.at(i, j) / scale - want).abs() < 1e-14);
                }
            }
        }
        let overflowing = mat(&[&[f64::MAX], &[f64::MAX]]);
        assert!(matches!(
            qr_pivoted(&overflowing),
            Err(LowRankError::NoConvergence(_))
        ));
    }

    #[test]
    fn interpolative_decomposition_uses_original_columns_in_input_order() {
        let a = mat(&[
            &[1.0, 2.0, 3.0, 2.0],
            &[2.0, 4.0, 6.0, 1.0],
            &[0.0, 0.0, 0.0, 3.0],
        ]);
        let result = interpolative_decomposition(&a, 3).unwrap();
        assert_eq!(result.columns, [2, 3]);
        let mut selected = Mat::zeros(a.rows, result.columns.len());
        for (j, &source) in result.columns.iter().enumerate() {
            for i in 0..a.rows {
                selected.set(i, j, a.at(i, source));
            }
        }
        let product = selected.matmul(&result.coefficients).unwrap();
        assert!(
            product
                .data
                .iter()
                .zip(&a.data)
                .all(|(x, y)| (x - y).abs() < 1e-14)
        );
        let rank_one = interpolative_decomposition(&a, 1).unwrap();
        assert_eq!(rank_one.columns, [2]);
        assert_eq!(rank_one.coefficients.at(0, 2), 1.0);
        let zero = interpolative_decomposition(&Mat::zeros(2, 3), 2).unwrap();
        assert_eq!(zero.columns, [] as [usize; 0]);
        assert_eq!((zero.coefficients.rows, zero.coefficients.cols), (0, 3));
        assert!(interpolative_decomposition(&a, 4).is_err());
        assert!(interpolative_decomposition(&mat(&[&[f64::NAN]]), 0).is_err());
    }

    #[test]
    fn randomized_svd_saturates_oversampling_and_rejects_any_invalid_product() {
        let identity = |x: &[f64], y: &mut [f64]| y.copy_from_slice(x);
        let options = RandomizedOptions {
            oversample: usize::MAX,
            ..RandomizedOptions::default()
        };
        let result = randomized_svd(&identity, &identity, 2, 2, 1, &options).unwrap();
        assert!((result.s[0] - 1.0).abs() < 1e-14);
        let first = std::cell::Cell::new(true);
        let invalid = |x: &[f64], y: &mut [f64]| {
            y.copy_from_slice(x);
            if first.replace(false) {
                y[0] = f64::INFINITY;
            }
        };
        assert!(matches!(
            randomized_svd(&invalid, &identity, 2, 2, 1, &options),
            Err(LowRankError::Shape(_))
        ));
    }

    #[test]
    fn randomized_svd_finds_the_leading_singular_triplets() {
        let spectrum: Vec<f64> = (0..30).map(|i| 0.6_f64.powi(i)).collect();
        let a = with_spectrum(120, 80, &spectrum, 11);
        let at = a.transpose();
        let mv = |x: &[f64], y: &mut [f64]| y.copy_from_slice(&a.matvec(x));
        let rmv = |y: &[f64], x: &mut [f64]| x.copy_from_slice(&at.matvec(y));
        let k = 6;
        let svd = randomized_svd(&mv, &rmv, 120, 80, k, &RandomizedOptions::default()).unwrap();
        for i in 0..k {
            assert!(
                (svd.s[i] - spectrum[i]).abs() < 1e-8 * spectrum[i],
                "sigma {i}: {} vs {}",
                svd.s[i],
                spectrum[i]
            );
        }
        // |A - U S V^T|_F is the tail of the spectrum, up to a small factor.
        let mut error = 0.0;
        for i in 0..120 {
            for j in 0..80 {
                let approx: f64 = (0..k)
                    .map(|t| svd.u.at(i, t) * svd.s[t] * svd.v.at(j, t))
                    .sum();
                error += (a.at(i, j) - approx).powi(2);
            }
        }
        let tail: f64 = spectrum[k..].iter().map(|s| s * s).sum::<f64>();
        assert!(
            error.sqrt() < 1.01 * tail.sqrt(),
            "{} vs {}",
            error.sqrt(),
            tail.sqrt()
        );
        assert!(randomized_svd(&mv, &rmv, 120, 80, 0, &RandomizedOptions::default()).is_err());
        assert!(randomized_svd(&mv, &rmv, 120, 80, 81, &RandomizedOptions::default()).is_err());
    }

    /// The 1-D Dirichlet Laplacian `tridiag(-1, 2, -1)`.
    fn laplacian(
        x: &[f64],
        y: &mut [f64],
    ) {
        let n = x.len();
        for i in 0..n {
            let left = if i > 0 { x[i - 1] } else { 0.0 };
            let right = if i + 1 < n { x[i + 1] } else { 0.0 };
            y[i] = 2.0 * x[i] - left - right;
        }
    }

    #[test]
    fn lobpcg_finds_the_smallest_eigenvalues_of_the_laplacian() {
        let n = 200;
        let exact =
            |j: usize| 2.0 - 2.0 * (std::f64::consts::PI * j as f64 / (n as f64 + 1.0)).cos();
        let pairs = lobpcg(&laplacian, None, n, 4, 1e-7, 3000, 5).unwrap();
        for (j, value) in pairs.values.iter().enumerate() {
            assert!(
                (value - exact(j + 1)).abs() < 1e-9,
                "eigenvalue {j}: {value} vs {}",
                exact(j + 1)
            );
        }
        // A good preconditioner (an exact tridiagonal solve) converges in a
        // handful of iterations.
        let thomas = |r: &[f64], z: &mut [f64]| {
            let n = r.len();
            let mut c = vec![0.0; n];
            let mut d = vec![0.0; n];
            c[0] = -0.5;
            d[0] = r[0] / 2.0;
            for i in 1..n {
                let m = 2.0 + c[i - 1];
                c[i] = -1.0 / m;
                d[i] = (r[i] + d[i - 1]) / m;
            }
            z[n - 1] = d[n - 1];
            for i in (0..n - 1).rev() {
                z[i] = d[i] - c[i] * z[i + 1];
            }
        };
        let fast = lobpcg(&laplacian, Some(&thomas), n, 4, 1e-7, 200, 5).unwrap();
        assert!(
            fast.iterations * 5 < pairs.iterations,
            "{} vs {}",
            fast.iterations,
            pairs.iterations
        );
        assert!((fast.values[0] - exact(1)).abs() < 1e-9);
        // Eigenvectors are orthonormal and satisfy the equation.
        for i in 0..4 {
            let mut av = vec![0.0; n];
            laplacian(&fast.vectors[i], &mut av);
            let r: f64 = av
                .iter()
                .zip(&fast.vectors[i])
                .map(|(a, v)| (a - fast.values[i] * v).powi(2))
                .sum::<f64>()
                .sqrt();
            assert!(r < 1e-6);
            assert!((norm(&fast.vectors[i]) - 1.0).abs() < 1e-10);
        }
        assert!(lobpcg(&laplacian, None, 10, 4, 1e-8, 10, 1).is_err());
        assert!(matches!(
            lobpcg(&laplacian, None, n, 4, 1e-12, 2, 5),
            Err(LowRankError::NoConvergence(_))
        ));
    }

    #[test]
    fn anderson_accelerates_a_slow_contraction() {
        // x = cos(x) componentwise with a coupling: linear rate about 0.67.
        let g = |x: &[f64]| -> Vec<f64> {
            let n = x.len();
            (0..n)
                .map(|i| (0.9 * x[i] + 0.1 * x[(i + 1) % n]).cos())
                .collect()
        };
        let start = vec![0.0; 8];
        let plain = anderson(&g, &start, 0, 1.0, 1e-12, 2000).unwrap();
        let fast = anderson(&g, &start, 5, 1.0, 1e-12, 2000).unwrap();
        assert!(
            fast.evaluations * 3 < plain.evaluations,
            "{} vs {}",
            fast.evaluations,
            plain.evaluations
        );
        for (a, b) in plain.x.iter().zip(&fast.x) {
            assert!((a - b).abs() < 1e-10);
            assert!((a - 0.739_085_133_215_160_6).abs() < 1e-10);
        }
        // A linear map with spectral radius just below one: hopeless for the
        // plain iteration within the budget, exact for the accelerated one.
        let h = |x: &[f64]| vec![0.999 * x[0] + 1.0, 0.5 * x[1] - 2.0];
        assert!(matches!(
            anderson(&h, &[0.0, 0.0], 0, 1.0, 1e-10, 200),
            Err(LowRankError::NoConvergence(_))
        ));
        let solved = anderson(&h, &[0.0, 0.0], 3, 1.0, 1e-8, 200).unwrap();
        assert!((solved.x[0] - 1000.0).abs() < 1e-4 && (solved.x[1] + 4.0).abs() < 1e-6);
        // Degenerate inputs.
        assert!(anderson(&g, &[f64::NAN], 2, 1.0, 1e-8, 10).is_err());
        assert!(anderson(&|_| vec![1.0, 2.0], &[0.0], 2, 1.0, 1e-8, 10).is_err());
        assert!(anderson(&|x| vec![x[0] * 2.0 + 1.0], &[1.0], 0, 1.0, 1e-8, 5000).is_err());
    }
}

#[cfg(test)]
mod generalized_lobpcg_tests {
    use super::*;
    use std::cell::Cell;

    fn close(
        left: f64,
        right: f64,
        tolerance: f64,
    ) {
        assert!((left - right).abs() < tolerance, "{left} != {right}");
    }

    #[test]
    fn weighted_lobpcg_solves_a_nondiagonal_pencil() -> Result<(), LowRankError> {
        let n = 18;
        let apply = |x: &[f64], out: &mut [f64]| {
            out[0] = 5.0 * x[0] - x[1];
            out[1] = 5.0 * x[1] - x[0];
            for i in 2..n {
                out[i] = if i == 2 {
                    2.0 * x[i]
                } else {
                    (i as f64 + 2.0) * x[i]
                };
            }
        };
        let mass = |x: &[f64], out: &mut [f64]| {
            out[0] = 3.0 * x[0] + x[1];
            out[1] = 3.0 * x[1] + x[0];
            out[2..].copy_from_slice(&x[2..]);
        };
        let result = lobpcg_generalized(
            &apply,
            &mass,
            None,
            n,
            3,
            LobpcgOptions {
                tolerance: 1e-9,
                max_iterations: 250,
                ..LobpcgOptions::default()
            },
        )?;
        for (&actual, expected) in result.values.iter().zip([1.0, 2.0, 3.0]) {
            close(actual, expected, 1e-10);
        }
        for (i, x) in result.vectors.iter().enumerate() {
            let mut bx = vec![0.0; n];
            let mut ax = vec![0.0; n];
            mass(x, &mut bx);
            apply(x, &mut ax);
            for (j, y) in result.vectors.iter().enumerate() {
                close(dot(&bx, y), if i == j { 1.0 } else { 0.0 }, 1e-9);
            }
            axpy(&mut ax, -result.values[i], &bx);
            close(norm(&ax), result.residuals[i], 1e-13);
            assert!(result.residuals[i] < 3.1e-9);
        }
        Ok(())
    }

    #[test]
    fn soft_locking_can_skip_converged_preconditioner_calls() -> Result<(), LowRankError> {
        let n = 30;
        let diagonal: Vec<_> = (0..n)
            .map(|i| {
                if i == 0 {
                    0.1
                } else {
                    i as f64 + 2.0
                }
            })
            .collect();
        let apply = |x: &[f64], out: &mut [f64]| {
            for ((slot, value), scale) in out.iter_mut().zip(x).zip(&diagonal) {
                *slot = value * scale;
            }
        };
        let mass = |x: &[f64], out: &mut [f64]| out.copy_from_slice(x);
        let calls = Cell::new(0_usize);
        let precondition = |x: &[f64], out: &mut [f64]| {
            calls.set(calls.get().saturating_add(1));
            for ((slot, value), scale) in out.iter_mut().zip(x).zip(&diagonal) {
                *slot = value / scale;
            }
        };
        let mut options = LobpcgOptions {
            tolerance: 1e-7,
            ..LobpcgOptions::default()
        };
        let locked = lobpcg_generalized(&apply, &mass, Some(&precondition), n, 3, options)?;
        let locked_calls = calls.replace(0);
        options.soft_lock = false;
        let free = lobpcg_generalized(&apply, &mass, Some(&precondition), n, 3, options)?;
        for (&lock_value, &free_value) in locked.values.iter().zip(&free.values) {
            close(lock_value, free_value, 1e-10);
        }
        assert!(locked_calls < locked.iterations.saturating_mul(locked.values.len()));
        assert_eq!(
            calls.get(),
            free.iterations.saturating_mul(free.values.len())
        );
        Ok(())
    }

    #[test]
    fn a_full_initial_block_supports_small_indefinite_stiffness() -> Result<(), LowRankError> {
        let apply = |x: &[f64], out: &mut [f64]| {
            out[0] = 2.0 * x[0] - 6.0 * x[1];
            out[1] = 2.0 * x[1] - 6.0 * x[0];
        };
        let mass = |x: &[f64], out: &mut [f64]| {
            out[0] = 3.0 * x[0] + x[1];
            out[1] = 3.0 * x[1] + x[0];
        };
        let options = LobpcgOptions {
            max_iterations: 0,
            ..LobpcgOptions::default()
        };
        let result = lobpcg_generalized(&apply, &mass, None, 2, 2, options)?;
        close(result.values[0], -1.0, 1e-12);
        close(result.values[1], 4.0, 1e-12);
        assert_eq!(result.iterations, 0);
        Ok(())
    }

    #[test]
    fn invalid_mass_and_workspace_are_refused() {
        let identity = |x: &[f64], out: &mut [f64]| out.copy_from_slice(x);
        let indefinite = |x: &[f64], out: &mut [f64]| {
            for (y, v) in out.iter_mut().zip(x) {
                *y = -*v;
            }
        };
        assert!(matches!(
            lobpcg_generalized(&identity, &indefinite, None, 4, 1, LobpcgOptions::default()),
            Err(LowRankError::Shape(_)),
        ));
        assert!(matches!(
            lobpcg_generalized(
                &identity,
                &identity,
                None,
                usize::MAX,
                2,
                LobpcgOptions::default()
            ),
            Err(LowRankError::Shape(_)),
        ));
    }
}
