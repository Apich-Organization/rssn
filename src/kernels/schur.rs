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

/// Solves the continuous-time Sylvester equation `A X + X B = C` via the
/// Bartels-Stewart algorithm (real Schur reduction + block back-substitution).
///
/// # Errors
///
/// [`DenseError::Shape`] on incompatible dimensions, [`DenseError::NonFinite`]
/// on non-finite entries, or [`DenseError::Singular`] when `A` and `-B` share
/// an eigenvalue.
pub fn sylvester(
    a: &Mat,
    b: &Mat,
    c: &Mat,
) -> Result<Mat, DenseError> {
    if a.rows != a.cols || b.rows != b.cols || c.rows != a.rows || c.cols != b.cols {
        return Err(DenseError::Shape);
    }
    if !a.is_finite() || !b.is_finite() || !c.is_finite() {
        return Err(DenseError::NonFinite);
    }
    let n = a.rows;
    let m = b.rows;
    if n == 0 || m == 0 {
        return Ok(Mat::zeros(n, m));
    }

    let sa = schur(a)?;
    let sb = schur(b)?;
    let u = sa.q;
    let r = sa.t;
    let v = sb.q;
    let s = sb.t;

    let f = u.transpose().matmul(c)?.matmul(&v)?;

    let mut row_blocks: Vec<(usize, usize)> = Vec::new();
    let mut i = n;
    while i > 0 {
        if i >= 2
            && r.at(i.saturating_sub(1), i.saturating_sub(2)).abs()
                > 1e-12
                    * (r.at(i.saturating_sub(2), i.saturating_sub(2)).abs()
                        + r.at(i.saturating_sub(1), i.saturating_sub(1)).abs()
                        + 1.0)
        {
            let start = i.saturating_sub(2);
            row_blocks.push((start, 2));
            i = start;
        } else {
            let start = i.saturating_sub(1);
            row_blocks.push((start, 1));
            i = start;
        }
    }

    let mut col_blocks: Vec<(usize, usize)> = Vec::new();
    let mut j = 0;
    while j < m {
        if j.saturating_add(1) < m
            && s.at(j.saturating_add(1), j).abs()
                > 1e-12
                    * (s.at(j, j).abs()
                        + s.at(j.saturating_add(1), j.saturating_add(1)).abs()
                        + 1.0)
        {
            col_blocks.push((j, 2));
            j = j.saturating_add(2);
        } else {
            col_blocks.push((j, 1));
            j = j.saturating_add(1);
        }
    }

    let mut y = Mat::zeros(n, m);
    for &(c0, cw) in &col_blocks {
        for &(r0, rw) in &row_blocks {
            let mut rhs = [0.0; 4];
            for dr in 0..rw {
                for dc in 0..cw {
                    let mut val = f.at(r0.saturating_add(dr), c0.saturating_add(dc));
                    for k in r0.saturating_add(rw)..n {
                        val -= r.at(r0.saturating_add(dr), k) * y.at(k, c0.saturating_add(dc));
                    }
                    for l in 0..c0 {
                        val -= y.at(r0.saturating_add(dr), l) * s.at(l, c0.saturating_add(dc));
                    }
                    if let Some(slot) = rhs.get_mut(dr.saturating_mul(cw).saturating_add(dc)) {
                        *slot = val;
                    }
                }
            }
            let dim = rw.saturating_mul(cw);
            let mut sys = Mat::zeros(dim, dim);
            for dc in 0..cw {
                for dr in 0..rw {
                    let eq = dc.saturating_mul(rw).saturating_add(dr);
                    for dr2 in 0..rw {
                        let col = dc.saturating_mul(rw).saturating_add(dr2);
                        sys.set(
                            eq,
                            col,
                            sys.at(eq, col) + r.at(r0.saturating_add(dr), r0.saturating_add(dr2)),
                        );
                    }
                    for dc2 in 0..cw {
                        let col = dc2.saturating_mul(rw).saturating_add(dr);
                        sys.set(
                            eq,
                            col,
                            sys.at(eq, col) + s.at(c0.saturating_add(dc2), c0.saturating_add(dc)),
                        );
                    }
                }
            }
            let mut rhs_vec = vec![0.0; dim];
            for dc in 0..cw {
                for dr in 0..rw {
                    if let Some(slot) = rhs_vec.get_mut(dc.saturating_mul(rw).saturating_add(dr)) {
                        *slot = rhs
                            .get(dr.saturating_mul(cw).saturating_add(dc))
                            .copied()
                            .unwrap_or(0.0);
                    }
                }
            }
            let sol = dense::solve(&sys, &rhs_vec)?;
            for dc in 0..cw {
                for dr in 0..rw {
                    y.set(
                        r0.saturating_add(dr),
                        c0.saturating_add(dc),
                        sol.get(dc.saturating_mul(rw).saturating_add(dr))
                            .copied()
                            .unwrap_or(0.0),
                    );
                }
            }
        }
    }

    u.matmul(&y)?.matmul(&v.transpose())
}

/// Solves the continuous-time Lyapunov equation `A X + X A^T = Q`.
///
/// # Errors
///
/// As [`sylvester`].
pub fn lyapunov(
    a: &Mat,
    q: &Mat,
) -> Result<Mat, DenseError> {
    if a.rows != a.cols || q.rows != a.rows || q.cols != a.cols {
        return Err(DenseError::Shape);
    }
    sylvester(a, &a.transpose(), q)
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
