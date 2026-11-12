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
//! * [`randomized_svd`]: a rank-`k` approximation of an operator that is only
//!   available as products `A x` and `A^T y` (Halko, Martinsson and Tropp),
//!   in `O(k)` products instead of a full decomposition;
//! * [`lobpcg`]: the smallest eigenpairs of a large symmetric operator by
//!   locally optimal block preconditioned conjugate gradients (Knyazev), which
//!   needs products and an optional preconditioner but no factorisation;
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
        Self {
            state: seed | 1,
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
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn norm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
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
                for (x, y) in v.iter_mut().zip(q) {
                    *x -= c * y;
                }
            }
        }
        let after = norm(&v);
        // What survives the projection must be more than rounding noise.
        if after > 1e-10 * before {
            for x in &mut v {
                *x /= after;
            }
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
        let k = self.r.rows.min(self.r.cols);
        let top = if k == 0 {
            0.0
        } else {
            self.r.at(0, 0).abs()
        };
        let tol = tol.unwrap_or(self.q.rows.max(self.r.cols) as f64 * f64::EPSILON);
        (0..k)
            .take_while(|&i| self.r.at(i, i).abs() > tol * top && top > 0.0)
            .count()
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
        for i in (0..rank).rev() {
            for j in i + 1..rank {
                y[i] -= self.r.at(i, j) * y[j];
            }
            y[i] /= self.r.at(i, i);
        }
        let mut x = vec![0.0; self.r.cols];
        for (j, value) in y.into_iter().enumerate() {
            x[self.perm[j]] = value;
        }
        Ok(x)
    }
}

/// Column-pivoted Householder QR (Businger and Golub).
///
/// At each step the remaining column of largest norm is brought forward, so
/// the diagonal of `R` decreases in magnitude and its decay reveals the
/// numerical rank. Column norms are downdated and recomputed when
/// cancellation makes the downdate unreliable.
///
/// # Errors
/// [`LowRankError::Shape`] for inconsistent storage or non-finite entries.
pub fn qr_pivoted(a: &Mat) -> Result<PivotedQr, LowRankError> {
    let (m, n) = (a.rows, a.cols);
    if m.checked_mul(n) != Some(a.data.len()) || !a.is_finite() {
        return Err(LowRankError::Shape(
            "the matrix storage is inconsistent or has non-finite entries".into(),
        ));
    }
    let k = m.min(n);
    let mut r = a.clone();
    let mut perm: Vec<usize> = (0..n).collect();
    let mut norms: Vec<f64> = (0..n)
        .map(|j| (0..m).map(|i| r.at(i, j).powi(2)).sum::<f64>().sqrt())
        .collect();
    let mut reference = norms.clone();
    let mut reflectors: Vec<Vec<f64>> = Vec::with_capacity(k);
    for step in 0..k {
        let pivot = (step..n)
            .max_by(|&x, &y| norms[x].total_cmp(&norms[y]))
            .unwrap_or(step);
        if pivot != step {
            for i in 0..m {
                let (x, y) = (r.at(i, step), r.at(i, pivot));
                r.set(i, step, y);
                r.set(i, pivot, x);
            }
            perm.swap(step, pivot);
            norms.swap(step, pivot);
            reference.swap(step, pivot);
        }
        // Householder vector for column `step`, rows `step..m`.
        let mut v: Vec<f64> = (step..m).map(|i| r.at(i, step)).collect();
        let alpha = norm(&v);
        if alpha > 0.0 {
            let beta = if v[0] >= 0.0 { -alpha } else { alpha };
            v[0] -= beta;
            let vn = norm(&v);
            if vn > 0.0 {
                for x in &mut v {
                    *x /= vn;
                }
                for j in step..n {
                    let s: f64 = 2.0 * (step..m).map(|i| v[i - step] * r.at(i, j)).sum::<f64>();
                    for i in step..m {
                        let value = r.at(i, j) - s * v[i - step];
                        r.set(i, j, value);
                    }
                }
            }
        } else {
            v.fill(0.0);
        }
        for i in step + 1..m {
            r.set(i, step, 0.0);
        }
        // Downdate the norms of the remaining columns.
        for j in step + 1..n {
            if norms[j] > 0.0 {
                let ratio = r.at(step, j) / norms[j];
                let remaining = (1.0 - ratio * ratio).max(0.0);
                norms[j] *= remaining.sqrt();
                if norms[j] < 1e-7 * reference[j] {
                    norms[j] = (step + 1..m)
                        .map(|i| r.at(i, j).powi(2))
                        .sum::<f64>()
                        .sqrt();
                    reference[j] = norms[j];
                }
            }
        }
        reflectors.push(v);
    }
    // Thin Q: apply the reflectors to the first k columns of the identity.
    let mut q = Mat::zeros(m, k);
    for j in 0..k {
        q.set(j, j, 1.0);
    }
    for (step, v) in reflectors.iter().enumerate().rev() {
        for j in 0..k {
            let s: f64 = 2.0 * (step..m).map(|i| v[i - step] * q.at(i, j)).sum::<f64>();
            for i in step..m {
                let value = q.at(i, j) - s * v[i - step];
                q.set(i, j, value);
            }
        }
    }
    let mut thin = Mat::zeros(k, n);
    for i in 0..k {
        for j in i..n {
            thin.set(i, j, r.at(i, j));
        }
    }
    Ok(PivotedQr { q, r: thin, perm })
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
/// `2 (power_iterations + 1) (k + oversample)` products.
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
    let l = (k + options.oversample).min(m.min(n));
    let mut rng = Gaussian::new(options.seed);
    let apply = |f: &dyn Fn(&[f64], &mut [f64]), xs: &[Vec<f64>], len: usize| -> Vec<Vec<f64>> {
        xs.iter()
            .map(|x| {
                let mut y = vec![0.0; len];
                f(x, &mut y);
                y
            })
            .collect()
    };
    let omega: Vec<Vec<f64>> = (0..l)
        .map(|_| (0..n).map(|_| rng.sample()).collect())
        .collect();
    let mut q = apply(matvec, &omega, m);
    orthonormalize(&[], &mut q);
    for _ in 0..options.power_iterations {
        let mut z = apply(rmatvec, &q, n);
        orthonormalize(&[], &mut z);
        q = apply(matvec, &z, m);
        orthonormalize(&[], &mut q);
    }
    if q.is_empty() || q.iter().flatten().any(|v| !v.is_finite()) {
        return Err(LowRankError::Shape(
            "the operator returned zero or non-finite products".into(),
        ));
    }
    // B = Q^T A, row i being (A^T q_i)^T.
    let rows = apply(rmatvec, &q, n);
    let mut b = Mat::zeros(rows.len(), n);
    for (i, row) in rows.iter().enumerate() {
        for (j, v) in row.iter().enumerate() {
            b.set(i, j, *v);
        }
    }
    let small = dense::svd(&b).map_err(|e| LowRankError::NoConvergence(format!("{e:?}")))?;
    let k = k.min(small.s.len());
    let mut u = Mat::zeros(m, k);
    for j in 0..k {
        for (i, qi) in q.iter().enumerate() {
            let c = small.u.at(i, j);
            for row in 0..m {
                let value = u.at(row, j) + c * qi[row];
                u.set(row, j, value);
            }
        }
    }
    let mut v = Mat::zeros(n, k);
    for j in 0..k {
        for i in 0..n {
            v.set(i, j, small.v.at(i, j));
        }
    }
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
    /// Residual norms `|A x - lambda x|`.
    pub residuals: Vec<f64>,
    /// Iterations used.
    pub iterations: usize,
}

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
/// finite; [`LowRankError::NoConvergence`] when `max_iterations` pass with a
/// residual above `tol * max(1, |lambda|)`.
pub fn lobpcg(
    apply: &dyn Fn(&[f64], &mut [f64]),
    precondition: Option<&dyn Fn(&[f64], &mut [f64])>,
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
                let value = 0.5 * (dot(&basis[i], &a_basis[j]) + dot(&basis[j], &a_basis[i]));
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
        if !(after > 1e-12 * before) {
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
