//! Real Schur, QZ (generalized Schur), Sylvester, and Lyapunov decompositions.
//!
//! Real Schur, QZ (generalized Schur), Sylvester, and continuous-time
//! Lyapunov solvers on [`Mat`].

use crate::kernels::dense::DenseError;
use crate::kernels::dense::Mat;
use crate::kernels::dense::{
    self,
};
use faer::dyn_stack::MemBuffer;
use faer::dyn_stack::MemStack;
use faer::linalg::evd::ComputeEigenvectors;
use faer::linalg::gevd::GeneralizedSchurParams;
use faer::linalg::gevd::qz_real;
use num_complex::Complex64;

/// Infinity norm `max_i sum_j |A_{i,j}|`.
#[must_use]
pub fn norm_inf(a: &Mat) -> f64 {
    let mut best = 0.0_f64;
    for i in 0..a.rows {
        let mut s = 0.0;
        for j in 0..a.cols {
            s += a.at(i, j).abs();
        }
        if s > best {
            best = s;
        }
    }
    best
}

/// Result of the real Schur decomposition `A = Q T Q^T`.
#[derive(Clone, Debug)]
pub struct Schur {
    /// Orthogonal matrix `Q`.
    pub q: Mat,
    /// Quasi-upper-triangular matrix `T`.
    pub t: Mat,
}

/// Real Schur decomposition `A = Q T Q^T`, where `Q` is orthogonal and `T` is
/// quasi-upper-triangular (1x1 and 2x2 diagonal blocks).
///
/// # Errors
///
/// [`DenseError::Shape`] for non-square input, [`DenseError::NonFinite`] for
/// non-finite entries.
pub fn schur(a: &Mat) -> Result<Schur, DenseError> {
    let (q, t) = dense::schur(a)?;
    Ok(Schur { q, t })
}

/// A reusable Bartels-Stewart reduction for `A X + X B = C`.
///
/// Keep the plan when the coefficient matrices are fixed and the right-hand side
/// changes. [`Self::for_lyapunov`] prepares the case `B = A^T` from one real Schur
/// decomposition.
#[derive(Debug)]
pub struct SylvesterPlan {
    rows: usize,
    cols: usize,
    reduced: Option<ReducedSylvester>,
}

#[derive(Debug)]
struct ReducedSylvester {
    left: Schur,
    left_transpose: Mat,
    right: Schur,
    right_transpose: Mat,
    row_blocks: Vec<(usize, usize)>,
    col_blocks: Vec<(usize, usize)>,
    large_coefficients: bool,
}

// Allow room for the sum of two coefficients and the growth of a four-row elimination.
const BLOCK_SCALE_LIMIT: f64 = f64::MAX / 16.0;

impl SylvesterPlan {
    pub(crate) fn heap_bytes(&self) -> usize {
        self.reduced.as_ref().map_or(0, |reduced| {
            [
                &reduced.left.q,
                &reduced.left.t,
                &reduced.left_transpose,
                &reduced.right.q,
                &reduced.right.t,
                &reduced.right_transpose,
            ]
            .into_iter()
            .map(|matrix| matrix.data.capacity() * size_of::<f64>())
            .sum::<usize>()
                + (reduced.row_blocks.capacity() + reduced.col_blocks.capacity())
                    * size_of::<(usize, usize)>()
        })
    }

    /// Prepares a Sylvester solver for square matrices `A` and `B`.
    ///
    /// # Errors
    /// Returns an error for invalid matrix storage, non-finite values,
    /// cancellation or a failed Schur reduction.
    pub fn new(
        a: &Mat,
        b: &Mat,
    ) -> Result<Self, DenseError> {
        if !valid_square(a) || !valid_square(b) {
            return Err(DenseError::Shape);
        }
        Self::from_valid_shapes(a, b)
    }

    /// Prepares a reusable solver for `A X + X A^T = Q`.
    ///
    /// # Errors
    /// Returns an error for invalid matrix storage, non-finite values,
    /// cancellation or a failed Schur reduction.
    pub fn for_lyapunov(a: &Mat) -> Result<Self, DenseError> {
        if !valid_square(a) {
            return Err(DenseError::Shape);
        }
        Self::from_valid_lyapunov_shape(a)
    }

    /// Solves the original equation for a new right-hand side with the plan's shape.
    ///
    /// # Errors
    /// Returns an error for inconsistent shape, non-finite input or arithmetic,
    /// cancellation, or a singular Sylvester operator.
    pub fn solve(
        &self,
        c: &Mat,
    ) -> Result<Mat, DenseError> {
        validate_rhs(c, self.rows, self.cols)?;
        self.solve_validated(c)
    }

    fn from_valid_shapes(
        a: &Mat,
        b: &Mat,
    ) -> Result<Self, DenseError> {
        if !a.is_finite() || !b.is_finite() {
            return Err(DenseError::NonFinite);
        }
        let reduced = if a.rows == 0 || b.rows == 0 {
            None
        } else {
            Some(ReducedSylvester::new(
                sylvester_schur(a)?,
                sylvester_schur(b)?,
            ))
        };
        Ok(Self {
            rows: a.rows,
            cols: b.rows,
            reduced,
        })
    }

    fn from_valid_lyapunov_shape(a: &Mat) -> Result<Self, DenseError> {
        if !a.is_finite() {
            return Err(DenseError::NonFinite);
        }
        let n = a.rows;
        let reduced = if n == 0 {
            None
        } else {
            let left = sylvester_schur(a)?;
            // Reversing coordinates makes T^T upper quasi-triangular again.
            // For the reversal J, A^T = (U J) (J T^T J) (U J)^T.
            let mut q = Mat::zeros(n, n);
            let mut t = Mat::zeros(n, n);
            for (i, ri) in (0..n).zip((0..n).rev()) {
                for (j, rj) in (0..n).zip((0..n).rev()) {
                    q.set(i, j, left.q.at(i, rj));
                    t.set(i, j, left.t.at(rj, ri));
                }
            }
            Some(ReducedSylvester::new(left, Schur { q, t }))
        };
        Ok(Self {
            rows: n,
            cols: n,
            reduced,
        })
    }

    fn solve_validated(
        &self,
        c: &Mat,
    ) -> Result<Mat, DenseError> {
        let Some(reduced) = &self.reduced else {
            return Ok(Mat::zeros(self.rows, self.cols));
        };
        reduced.solve(c)
    }
}

impl ReducedSylvester {
    fn new(
        left: Schur,
        right: Schur,
    ) -> Self {
        let left_transpose = left.q.transpose();
        let right_transpose = right.q.transpose();
        let row_blocks = schur_blocks(&left.t);
        let col_blocks = schur_blocks(&right.t);
        let large_coefficients = left
            .t
            .data
            .iter()
            .chain(&right.t.data)
            .any(|v| v.abs() >= BLOCK_SCALE_LIMIT);
        Self {
            left,
            left_transpose,
            right,
            right_transpose,
            row_blocks,
            col_blocks,
            large_coefficients,
        }
    }

    fn solve(
        &self,
        c: &Mat,
    ) -> Result<Mat, DenseError> {
        crate::cancellation::checkpoint()?;
        let f = self.left_transpose.matmul(c)?.matmul(&self.right.q)?;
        let r = &self.left.t;
        let s = &self.right.t;
        let mut y = Mat::zeros(r.rows, s.rows);
        for &(c0, cw) in &self.col_blocks {
            for &(r0, rw) in self.row_blocks.iter().rev() {
                crate::cancellation::checkpoint()?;
                let mut rhs = [0.0; 4];
                for dc in 0..cw {
                    let col = c0.saturating_add(dc);
                    for dr in 0..rw {
                        let row = r0.saturating_add(dr);
                        let mut value = f.at(row, col);
                        for k in r0.saturating_add(rw)..r.rows {
                            value -= r.at(row, k) * y.at(k, col);
                        }
                        for l in 0..c0 {
                            value -= y.at(row, l) * s.at(l, col);
                        }
                        if let Some(slot) = rhs.get_mut(dc.saturating_mul(rw).saturating_add(dr)) {
                            *slot = value;
                        }
                    }
                }
                if rw == 1 && cw == 1 {
                    let left = r.at(r0, r0);
                    let right = s.at(c0, c0);
                    let coefficient = left + right;
                    let [value, ..] = rhs;
                    if !value.is_finite() {
                        return Err(DenseError::NonFinite);
                    }
                    if coefficient == 0.0 {
                        return Err(DenseError::Singular);
                    }
                    let entry = if coefficient.is_finite() {
                        value / coefficient
                    } else {
                        (0.5 * value) / (0.5 * left + 0.5 * right)
                    };
                    y.set(r0, c0, entry);
                } else {
                    let solution =
                        solve_schur_block(r, s, (r0, rw), (c0, cw), rhs, self.large_coefficients)?;
                    for dc in 0..cw {
                        for dr in 0..rw {
                            let entry = dc.saturating_mul(rw).saturating_add(dr);
                            let value = solution.get(entry).copied().ok_or(DenseError::Shape)?;
                            y.set(r0.saturating_add(dr), c0.saturating_add(dc), value);
                        }
                    }
                }
            }
        }
        let result = self.left.q.matmul(&y)?.matmul(&self.right_transpose)?;
        result
            .is_finite()
            .then_some(result)
            .ok_or(DenseError::NonFinite)
    }
}

fn valid_square(a: &Mat) -> bool {
    a.rows == a.cols && a.rows.checked_mul(a.cols) == Some(a.data.len())
}

fn validate_rhs(
    c: &Mat,
    rows: usize,
    cols: usize,
) -> Result<(), DenseError> {
    if c.rows != rows || c.cols != cols || rows.checked_mul(cols) != Some(c.data.len()) {
        return Err(DenseError::Shape);
    }
    if !c.is_finite() {
        return Err(DenseError::NonFinite);
    }
    Ok(())
}

fn sylvester_schur(a: &Mat) -> Result<Schur, DenseError> {
    crate::cancellation::checkpoint()?;
    let peak = a
        .data
        .iter()
        .fold(0.0_f64, |max, value| max.max(value.abs()));
    // The Francis/Householder reduction can form expressions of degree four.
    // Normalize extreme operands independently to protect the Schur iteration.
    let result = if peak > 0.0 && !(1e-75..=1e75).contains(&peak) {
        let scaled = Mat {
            rows: a.rows,
            cols: a.cols,
            data: a.data.iter().map(|v| v / peak).collect(),
        };
        let mut result = schur(&scaled)?;
        for value in &mut result.t.data {
            *value *= peak;
        }
        result
    } else {
        schur(a)?
    };
    (result.q.is_finite() && result.t.is_finite())
        .then_some(result)
        .ok_or(DenseError::NonFinite)
}

fn schur_blocks(t: &Mat) -> Vec<(usize, usize)> {
    let mut blocks = Vec::new();
    let mut i = 0usize;
    while i < t.rows {
        let next = i.saturating_add(1);
        // The real Schur reduction writes exact zeros at deflated boundaries.
        let width = if next < t.rows && t.at(next, i) != 0.0 {
            2
        } else {
            1
        };
        blocks.push((i, width));
        i = i.saturating_add(width);
    }
    blocks
}

#[expect(
    clippy::single_call_fn,
    clippy::indexing_slicing,
    reason = "Schur block widths at most two give at most four equations"
)]
fn solve_schur_block(
    r: &Mat,
    s: &Mat,
    (r0, rw): (usize, usize),
    (c0, cw): (usize, usize),
    mut rhs: [f64; 4],
    large_coefficients: bool,
) -> Result<[f64; 4], DenseError> {
    let peak = if large_coefficients {
        (0..rw)
            .flat_map(|i| (0..rw).map(move |j| r.at(r0.saturating_add(i), r0.saturating_add(j))))
            .chain((0..cw).flat_map(|i| {
                (0..cw).map(move |j| s.at(c0.saturating_add(i), c0.saturating_add(j)))
            }))
            .fold(0.0_f64, |max, value| max.max(value.abs()))
    } else {
        0.0
    };
    let scale = if peak >= BLOCK_SCALE_LIMIT {
        let scale = peak.recip();
        for value in &mut rhs {
            *value *= scale;
        }
        scale
    } else {
        1.0
    };
    let dim = rw.saturating_mul(cw);
    let mut sys = [[0.0; 4]; 4];
    for dc in 0..cw {
        for dr in 0..rw {
            let eq = dc.saturating_mul(rw).saturating_add(dr);
            for dr2 in 0..rw {
                let col = dc.saturating_mul(rw).saturating_add(dr2);
                sys[eq][col] += r.at(r0.saturating_add(dr), r0.saturating_add(dr2)) * scale;
            }
            for dc2 in 0..cw {
                let col = dc2.saturating_mul(rw).saturating_add(dr);
                sys[eq][col] += s.at(c0.saturating_add(dc2), c0.saturating_add(dc)) * scale;
            }
        }
    }
    if !sys.iter().flatten().chain(&rhs).all(|v| v.is_finite()) {
        return Err(DenseError::NonFinite);
    }
    // Apply partial pivoting to the right-hand side as the small system is eliminated.
    for k in 0..dim {
        crate::cancellation::checkpoint()?;
        let below = k.saturating_add(1)..dim;
        let mut pivot = k;
        for row in below.clone() {
            if sys[row][k].abs() > sys[pivot][k].abs() {
                pivot = row;
            }
        }
        if sys[pivot][k] == 0.0 {
            return Err(DenseError::Singular);
        }
        sys.swap(pivot, k);
        rhs.swap(pivot, k);
        for row in below.clone() {
            let multiplier = sys[row][k] / sys[k][k];
            if multiplier != 0.0 {
                for col in below.clone() {
                    sys[row][col] -= multiplier * sys[k][col];
                }
                rhs[row] -= multiplier * rhs[k];
            }
        }
    }
    for row in (0..dim).rev() {
        for col in row.saturating_add(1)..dim {
            rhs[row] -= sys[row][col] * rhs[col];
        }
        rhs[row] /= sys[row][row];
    }
    rhs.iter()
        .all(|v| v.is_finite())
        .then_some(rhs)
        .ok_or(DenseError::NonFinite)
}

/// Solves the continuous-time Sylvester equation `A X + X B = C` by real
/// Schur reduction and block back-substitution.
///
/// # Errors
/// Returns an error for invalid shapes, non-finite input or arithmetic,
/// cancellation, a failed Schur reduction, or a singular Sylvester operator.
pub fn sylvester(
    a: &Mat,
    b: &Mat,
    c: &Mat,
) -> Result<Mat, DenseError> {
    if !valid_square(a) || !valid_square(b) {
        return Err(DenseError::Shape);
    }
    validate_rhs(c, a.rows, b.rows)?;
    SylvesterPlan::from_valid_shapes(a, b)?.solve_validated(c)
}

/// Solves the continuous-time Lyapunov equation `A X + X A^T = Q`.
///
/// # Errors
/// As [`sylvester`].
pub fn lyapunov(
    a: &Mat,
    q: &Mat,
) -> Result<Mat, DenseError> {
    if !valid_square(a) {
        return Err(DenseError::Shape);
    }
    validate_rhs(q, a.rows, a.rows)?;
    SylvesterPlan::from_valid_lyapunov_shape(a)?.solve_validated(q)
}

/// Result of the real QZ (generalized Schur) decomposition of `(A, B)`:
/// `A = Q S Z^T` and `B = Q T Z^T`.
#[derive(Clone, Debug)]
pub struct Qz {
    /// Quasi-upper-triangular `S`.
    pub s: Mat,
    /// Upper-triangular `T`.
    pub t: Mat,
    /// Left orthogonal matrix `Q`.
    pub q: Mat,
    /// Right orthogonal matrix `Z`.
    pub z: Mat,
    /// Generalized eigenvalue numerators `alpha`.
    pub alpha: Vec<Complex64>,
    /// Generalized eigenvalue denominators `beta`.
    pub beta: Vec<f64>,
}

fn givens(
    x: f64,
    y: f64,
) -> (f64, f64) {
    if y == 0.0 {
        (1.0, 0.0)
    } else {
        let r = x.hypot(y);
        (x / r, y / r)
    }
}

fn rot_rows(
    m: &mut Mat,
    i1: usize,
    i2: usize,
    c: f64,
    sg: f64,
    col_start: usize,
) {
    for j in col_start..m.cols {
        let v1 = m.at(i1, j);
        let v2 = m.at(i2, j);
        m.set(i1, j, c * v1 + sg * v2);
        m.set(i2, j, -sg * v1 + c * v2);
    }
}

fn rot_cols(
    m: &mut Mat,
    j1: usize,
    j2: usize,
    c: f64,
    sg: f64,
    row_end: usize,
) {
    for i in 0..row_end {
        let v1 = m.at(i, j1);
        let v2 = m.at(i, j2);
        m.set(i, j1, c * v1 + sg * v2);
        m.set(i, j2, -sg * v1 + c * v2);
    }
}

/// Real QZ (generalized Schur) decomposition of `(A, B)`.
///
/// # Errors
///
/// [`DenseError::Shape`], [`DenseError::NonFinite`], or [`DenseError::NoConvergence`].
pub fn qz(
    a: &Mat,
    b: &Mat,
) -> Result<Qz, DenseError> {
    if a.rows != a.cols || b.rows != b.cols || a.rows != b.rows {
        return Err(DenseError::Shape);
    }
    if !a.is_finite() || !b.is_finite() {
        return Err(DenseError::NonFinite);
    }
    let n = a.rows;
    if n == 0 {
        return Ok(Qz {
            s: Mat::zeros(0, 0),
            t: Mat::zeros(0, 0),
            q: Mat::zeros(0, 0),
            z: Mat::zeros(0, 0),
            alpha: Vec::new(),
            beta: Vec::new(),
        });
    }
    if n == 1 {
        return Ok(Qz {
            s: a.clone(),
            t: b.clone(),
            q: Mat::identity(1),
            z: Mat::identity(1),
            alpha: vec![Complex64::new(a.at(0, 0), 0.0)],
            beta: vec![b.at(0, 0)],
        });
    }

    let (mut q, mut t) = dense::qr(b);
    let mut s = q.transpose().matmul(a)?;
    let mut z = Mat::identity(n);

    // Hessenberg-triangular reduction.
    for j in 0..n.saturating_sub(2) {
        for i in (j.saturating_add(2)..n).rev() {
            let im1 = i.saturating_sub(1);
            let (c, sg) = givens(s.at(im1, j), s.at(i, j));
            rot_rows(&mut s, im1, i, c, sg, j);
            s.set(i, j, 0.0);
            rot_rows(&mut t, im1, i, c, sg, im1);
            rot_cols(&mut q, im1, i, c, sg, n);

            let a1 = t.at(i, i);
            let a2 = t.at(i, im1);
            if a2 != 0.0 {
                let r = a1.hypot(a2);
                let cz = a1 / r;
                let sz = a2 / r;
                rot_cols(&mut t, im1, i, cz, -sz, i.saturating_add(1));
                t.set(i, im1, 0.0);
                rot_cols(&mut s, im1, i, cz, -sz, n);
                rot_cols(&mut z, im1, i, cz, -sz, n);
            }
        }
    }

    // Use faer's multi-shift QZ iteration, which explicitly deflates zero
    // diagonals of B. Treating a singular B as a shifted standard eigenproblem
    // loses infinite eigenvalues and can stall a hand-written single-shift loop.
    let params = <GeneralizedSchurParams as faer::Auto<f64>>::auto();
    let par = faer::Par::Seq;
    let mut fs = faer::Mat::from_fn(n, n, |i, j| s.at(i, j));
    let mut ft = faer::Mat::from_fn(n, n, |i, j| t.at(i, j));
    let mut fq = faer::Mat::from_fn(n, n, |i, j| q.at(i, j));
    let mut fz = faer::Mat::from_fn(n, n, |i, j| z.at(i, j));
    let mut ar = faer::Col::<f64>::zeros(n);
    let mut ai = faer::Col::<f64>::zeros(n);
    let mut be = faer::Col::<f64>::zeros(n);
    let mut scratch = MemBuffer::new(qz_real::hessenberg_to_qz_scratch::<f64>(n, par, params));
    qz_real::hessenberg_to_qz(
        fs.as_mut(),
        ft.as_mut(),
        Some(fq.as_mut()),
        Some(fz.as_mut()),
        ar.as_mut(),
        ai.as_mut(),
        be.as_mut(),
        ComputeEigenvectors::Yes,
        par,
        params,
        MemStack::new(&mut scratch),
    );
    let alpha: Vec<_> = (0..n).map(|i| Complex64::new(ar[i], ai[i])).collect();
    let beta: Vec<_> = (0..n).map(|i| be[i]).collect();
    if alpha.iter().any(|v| !v.re.is_finite() || !v.im.is_finite())
        || beta.iter().any(|v| !v.is_finite())
    {
        return Err(DenseError::NoConvergence);
    }
    for i in 0..n {
        for j in 0..n {
            s.set(i, j, fs[(i, j)]);
            t.set(i, j, ft[(i, j)]);
            q.set(i, j, fq[(i, j)]);
            z.set(i, j, fz[(i, j)]);
        }
    }

    Ok(Qz {
        s,
        t,
        q,
        z,
        alpha,
        beta,
    })
}

/// Generalized eigenvalues `alpha_i / beta_i` of `A v = lambda B v`.
///
/// # Errors
///
/// As [`qz`].
pub fn eigvals_generalized(
    a: &Mat,
    b: &Mat,
) -> Result<Vec<Complex64>, DenseError> {
    let res = qz(a, b)?;
    Ok(res
        .alpha
        .into_iter()
        .zip(res.beta)
        .map(|(al, be)| {
            if be == 0.0 {
                Complex64::new(
                    if al.norm() == 0.0 {
                        f64::NAN
                    } else {
                        f64::INFINITY
                    },
                    0.0,
                )
            } else {
                Complex64::new(al.re / be, al.im / be)
            }
        })
        .collect())
}
