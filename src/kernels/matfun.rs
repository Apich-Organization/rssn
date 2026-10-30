//! # Matrix functions
//!
//! * [`expm`]: matrix exponential by scaling and squaring with diagonal
//!   Pade approximants of degree 3, 5, 7, 9 or 13 chosen from the 1-norm
//!   (Higham 2005 / Al-Mohy and Higham 2009 thresholds);
//! * [`sqrtm`]: principal square root by the Denman-Beavers iteration with
//!   determinantal scaling;
//! * [`logm`]: principal logarithm by inverse scaling and squaring with a
//!   Gauss-Legendre (Pade) partial-fraction evaluation of `log(I + X)`;
//! * [`phi_matrices`]: the exponential-integrator functions
//!   `phi_k(A)` (`phi_0 = exp`, `phi_{k+1}(z) = (phi_k(z) - 1/k!) / z`)
//!   through one exponential of an augmented matrix;
//! * [`expm_multiply`]: the action `exp(t A) v` for a sparse [`Csr`]
//!   matrix by scaled truncated Taylor series, never forming `exp(t A)`;
//!   [`expm_krylov`] does the same for a matrix-free operator with an
//!   Arnoldi/Krylov projection and adaptive sub-stepping.
//!
//! `sqrtm` and `logm` require that `A` has no eigenvalues on the closed
//! negative real axis (for `logm` also that `A` is nonsingular); they return
//! an error otherwise.

use crate::kernels::dense::DenseError;
use crate::kernels::dense::Mat;
use crate::kernels::dense::{self};
use crate::kernels::krylov::Csr;

fn mm(
    a: &Mat,
    b: &Mat,
) -> Result<Mat, DenseError> {
    a.matmul(b)
}

fn lin(
    a: &Mat,
    ca: f64,
    b: &Mat,
    cb: f64,
) -> Mat {
    let mut r = a.clone();
    for (x, y) in r.data.iter_mut().zip(&b.data) {
        *x = ca * *x + cb * y;
    }
    r
}

fn scaled_identity_add(
    a: &Mat,
    c: f64,
) -> Mat {
    let mut r = a.clone();
    for i in 0..r.rows {
        r.set(i, i, r.at(i, i) + c);
    }
    r
}

/// Matrix 1-norm (maximum absolute column sum).
#[must_use]
pub fn norm1(a: &Mat) -> f64 {
    (0..a.cols)
        .map(|j| (0..a.rows).map(|i| a.at(i, j).abs()).sum::<f64>())
        .fold(0.0, f64::max)
}

const PADE3: [f64; 4] = [120.0, 60.0, 12.0, 1.0];
const PADE5: [f64; 6] = [30240.0, 15120.0, 3360.0, 420.0, 30.0, 1.0];
const PADE7: [f64; 8] = [
    17297280.0, 8648640.0, 1995840.0, 277200.0, 25200.0, 1512.0, 56.0, 1.0,
];
const PADE9: [f64; 10] = [
    17643225600.0,
    8821612800.0,
    2075673600.0,
    302702400.0,
    30270240.0,
    2162160.0,
    110880.0,
    3960.0,
    90.0,
    1.0,
];
const PADE13: [f64; 14] = [
    64764752532480000.0,
    32382376266240000.0,
    7771770303897600.0,
    1187353796428800.0,
    129060195264000.0,
    10559470521600.0,
    670442572800.0,
    33522128640.0,
    1323241920.0,
    40840800.0,
    960960.0,
    16380.0,
    182.0,
    1.0,
];
const THETA: [f64; 5] = [
    1.495585217958292e-2,
    2.539398330063230e-1,
    9.504178996162932e-1,
    2.097847961257068,
    5.371920351148152,
];

/// Matrix exponential `exp(A)` of a square matrix.
///
/// # Errors
/// [`DenseError::Shape`] for non-square input, [`DenseError::Singular`] if a
/// Pade denominator is singular (not expected for finite input).
pub fn expm(a: &Mat) -> Result<Mat, DenseError> {
    if a.rows != a.cols {
        return Err(DenseError::Shape);
    }
    let n = a.rows;
    if n == 0 {
        return Ok(Mat::zeros(0, 0));
    }
    let nrm = norm1(a);
    if !nrm.is_finite() {
        return Err(DenseError::NoConvergence);
    }
    let id = Mat::identity(n);
    let pade_low = |a: &Mat, c: &[f64]| -> Result<(Mat, Mat), DenseError> {
        // U = A * (sum odd), V = sum even
        let a2 = mm(a, a)?;
        let deg = c.len() - 1;
        let mut powers: Vec<Mat> = vec![id.clone(), a2.clone()];
        for k in 2..=deg / 2 {
            let nx = mm(&powers[k - 1], &a2)?;
            powers.push(nx);
        }
        let mut u_in = Mat::zeros(n, n);
        let mut v = Mat::zeros(n, n);
        for k in 0..=deg / 2 {
            if 2 * k < deg {
                u_in = lin(&u_in, 1.0, &powers[k], c[2 * k + 1]);
            }
            v = lin(&v, 1.0, &powers[k], c[2 * k]);
        }
        Ok((mm(a, &u_in)?, v))
    };
    let degrees = [3usize, 5, 7, 9];
    let coefs: [&[f64]; 4] = [&PADE3, &PADE5, &PADE7, &PADE9];
    for (idx, &d) in degrees.iter().enumerate() {
        let _ = d;
        if nrm <= THETA[idx] {
            let (u, v) = pade_low(a, coefs[idx])?;
            return pade_solve(&u, &v);
        }
    }
    // degree 13 with scaling and squaring
    let s = if nrm > THETA[4] {
        (nrm / THETA[4]).log2().ceil().max(0.0) as i32
    } else {
        0
    };
    let a_s = lin(a, 1.0 / 2f64.powi(s), a, 0.0);
    let a2 = mm(&a_s, &a_s)?;
    let a4 = mm(&a2, &a2)?;
    let a6 = mm(&a4, &a2)?;
    let b = &PADE13;
    // U = A [ A6 (b13 A6 + b11 A4 + b9 A2) + b7 A6 + b5 A4 + b3 A2 + b1 I ]
    let t1 = lin(&lin(&a6, b[13], &a4, b[11]), 1.0, &a2, b[9]);
    let mut inner = mm(&a6, &t1)?;
    inner = lin(&inner, 1.0, &a6, b[7]);
    inner = lin(&inner, 1.0, &a4, b[5]);
    inner = lin(&inner, 1.0, &a2, b[3]);
    inner = scaled_identity_add(&inner, b[1]);
    let u = mm(&a_s, &inner)?;
    let t2 = lin(&lin(&a6, b[12], &a4, b[10]), 1.0, &a2, b[8]);
    let mut v = mm(&a6, &t2)?;
    v = lin(&v, 1.0, &a6, b[6]);
    v = lin(&v, 1.0, &a4, b[4]);
    v = lin(&v, 1.0, &a2, b[2]);
    v = scaled_identity_add(&v, b[0]);
    let mut r = pade_solve(&u, &v)?;
    for _ in 0..s {
        r = mm(&r, &r)?;
    }
    Ok(r)
}

fn pade_solve(
    u: &Mat,
    v: &Mat,
) -> Result<Mat, DenseError> {
    // (V - U)^{-1} (V + U)
    let den = lin(v, 1.0, u, -1.0);
    let num = lin(v, 1.0, u, 1.0);
    let lu = dense::lu_factor(&den)?;
    let n = num.rows;
    let mut out = Mat::zeros(n, n);
    for j in 0..n {
        let col: Vec<f64> = (0..n).map(|i| num.at(i, j)).collect();
        let x = lu.solve(&col);
        for i in 0..n {
            out.set(i, j, x[i]);
        }
    }
    Ok(out)
}

/// Principal square root of `A` (Denman-Beavers iteration with
/// determinantal scaling, quadratically convergent).
///
/// # Errors
/// [`DenseError::Shape`], or [`DenseError::NoConvergence`] /
/// [`DenseError::Singular`] when `A` has eigenvalues on the closed negative
/// real axis.
pub fn sqrtm(a: &Mat) -> Result<Mat, DenseError> {
    if a.rows != a.cols {
        return Err(DenseError::Shape);
    }
    let n = a.rows;
    if n == 0 {
        return Ok(Mat::zeros(0, 0));
    }
    let mut y = a.clone();
    let mut z = Mat::identity(n);
    for _ in 0..100 {
        let ly = dense::lu_factor(&y)?;
        let lz = dense::lu_factor(&z)?;
        let (dy, dz) = (ly.det(), lz.det());
        let mu = (dy * dz).abs().powf(-1.0 / (2.0 * n as f64));
        let mu = if mu.is_finite() && mu > 0.0 {
            mu
        } else {
            1.0
        };
        let yi = ly.inverse();
        let zi = lz.inverse();
        let yn = lin(&y, 0.5 * mu, &zi, 0.5 / mu);
        let zn = lin(&z, 0.5 * mu, &yi, 0.5 / mu);
        let diff = norm1(&lin(&yn, 1.0, &y, -1.0));
        y = yn;
        z = zn;
        if !diff.is_finite() {
            return Err(DenseError::NoConvergence);
        }
        if diff <= 1e-15 * norm1(&y).max(1.0) {
            return Ok(y);
        }
    }
    Err(DenseError::NoConvergence)
}

fn gauss_legendre_01(m: usize) -> (Vec<f64>, Vec<f64>) {
    let mut x = vec![0.0; m];
    let mut w = vec![0.0; m];
    for i in 0..m {
        let mut z = (std::f64::consts::PI * (i as f64 + 0.75) / (m as f64 + 0.5)).cos();
        let mut dp = 1.0;
        for _ in 0..100 {
            let (mut p0, mut p1) = (1.0, z);
            for k in 1..m {
                let kf = k as f64;
                let p2 = ((2.0 * kf + 1.0) * z * p1 - kf * p0) / (kf + 1.0);
                p0 = p1;
                p1 = p2;
            }
            dp = m as f64 * (z * p1 - p0) / (z * z - 1.0);
            let dz = p1 / dp;
            z -= dz;
            if dz.abs() < 1e-16 {
                break;
            }
        }
        x[i] = 0.5 * (1.0 - z);
        w[i] = 1.0 / ((1.0 - z * z) * dp * dp);
    }
    (x, w)
}

/// Principal logarithm of a matrix without eigenvalues on the closed
/// negative real axis, by inverse scaling and squaring.
///
/// # Errors
/// [`DenseError::Shape`] or the errors of [`sqrtm`].
pub fn logm(a: &Mat) -> Result<Mat, DenseError> {
    if a.rows != a.cols {
        return Err(DenseError::Shape);
    }
    let n = a.rows;
    if n == 0 {
        return Ok(Mat::zeros(0, 0));
    }
    let id = Mat::identity(n);
    let mut t = a.clone();
    let mut k = 0;
    while norm1(&lin(&t, 1.0, &id, -1.0)) > 0.2 {
        t = sqrtm(&t)?;
        k += 1;
        if k > 60 {
            return Err(DenseError::NoConvergence);
        }
    }
    let x = lin(&t, 1.0, &id, -1.0);
    // log(I + X) = sum_j w_j X (I + x_j X)^{-1}
    let (nodes, weights) = gauss_legendre_01(9);
    let mut r = Mat::zeros(n, n);
    for (xj, wj) in nodes.iter().zip(&weights) {
        let m = lin(&id, 1.0, &x, *xj);
        let lu = dense::lu_factor(&m)?;
        for c in 0..n {
            let col: Vec<f64> = (0..n).map(|i| x.at(i, c)).collect();
            let sol = lu.solve(&col);
            for i in 0..n {
                r.set(i, c, r.at(i, c) + wj * sol[i]);
            }
        }
    }
    // X (I + x X)^{-1} = (I + x X)^{-1} X, columns above are (I + xX)^{-1} X
    Ok(lin(&r, 2f64.powi(k), &r, 0.0))
}

/// The `phi_k(A)` functions, `k = 0..=p`, from one augmented matrix exponential.
///
/// The functions `phi_0(A), ..., phi_p(A)` where
/// `phi_0(z) = e^z` and `phi_k(z) = (phi_{k-1}(z) - 1/(k-1)!) / z`
/// (so `phi_1(z) = (e^z - 1)/z`), from one matrix exponential of the
/// augmented matrix `[[A, I, 0, ..], [0, 0, I, ..], ..]` (Sidje).
///
/// # Errors
/// As [`expm`].
pub fn phi_matrices(
    a: &Mat,
    p: usize,
) -> Result<Vec<Mat>, DenseError> {
    if a.rows != a.cols {
        return Err(DenseError::Shape);
    }
    let n = a.rows;
    let big = n * (p + 1);
    let mut m = Mat::zeros(big, big);
    for i in 0..n {
        for j in 0..n {
            m.set(i, j, a.at(i, j));
        }
    }
    for blk in 0..p {
        for i in 0..n {
            m.set(blk * n + i, (blk + 1) * n + i, 1.0);
        }
    }
    let e = expm(&m)?;
    let mut out = Vec::with_capacity(p + 1);
    for blk in 0..=p {
        let mut f = Mat::zeros(n, n);
        for i in 0..n {
            for j in 0..n {
                f.set(i, j, e.at(i, blk * n + j));
            }
        }
        out.push(f);
    }
    Ok(out)
}

/// Maximum absolute column sum of a sparse matrix.
#[must_use]
pub fn csr_norm1(a: &Csr) -> f64 {
    let mut cols = vec![0.0; a.cols];
    for (k, &c) in a.indices.iter().enumerate() {
        cols[c] += a.values[k].abs();
    }
    cols.into_iter().fold(0.0, f64::max)
}

/// Action `exp(t A) v` for a sparse matrix by a stepped Taylor series.
///
/// The action `exp(t A) v` for a sparse matrix, by splitting `[0, t]` into
/// `s` steps with `|t A|_1 / s <= 1` and summing the Taylor series of each
/// step until the terms are below `1e-16` of the partial sum.
#[must_use]
pub fn expm_multiply(
    a: &Csr,
    v: &[f64],
    t: f64,
) -> Vec<f64> {
    let n = v.len();
    let nrm = csr_norm1(a) * t.abs();
    let s = nrm.ceil().max(1.0) as usize;
    let h = t / s as f64;
    let mut w = v.to_vec();
    let mut term = vec![0.0; n];
    let mut tmp = vec![0.0; n];
    for _ in 0..s {
        term.clone_from(&w);
        let mut sum = w.clone();
        let mut small = 0;
        for j in 1..200 {
            a.matvec(&term, &mut tmp);
            for i in 0..n {
                term[i] = tmp[i] * h / f64::from(j);
                sum[i] += term[i];
            }
            let tn = term.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
            let sn = sum.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
            if tn <= 1e-17 * sn.max(f64::MIN_POSITIVE) {
                small += 1;
                if small >= 2 {
                    break;
                }
            } else {
                small = 0;
            }
        }
        w = sum;
    }
    w
}

/// Action `exp(t A) v` for a matrix-free operator by Krylov projection.
///
/// The action `exp(t A) v` for a matrix-free operator by Krylov projection:
/// `exp(t A) v ~ |v| V_m exp(t H_m) e_1` with an Arnoldi basis of dimension
/// `m`, advancing in sub-steps whose error estimate
/// `|v| h_{m+1,m} |e_m^T exp(tau H_m) e_1|` stays below `tol` per unit time.
///
/// Returns `None` if the sub-step shrinks below `1e-12 |t|`.
pub fn expm_krylov<A: Fn(&[f64], &mut [f64])>(
    apply: A,
    v: &[f64],
    t: f64,
    m: usize,
    tol: f64,
) -> Option<Vec<f64>> {
    let n = v.len();
    let m = m.clamp(2, n.max(2)).min(n.max(1));
    let mut w = v.to_vec();
    let mut done = 0.0;
    let mut tau = t;
    let mut guard = 0;
    while (t - done).abs() > 1e-15 * t.abs() {
        guard += 1;
        if guard > 100_000 {
            return None;
        }
        let rem = t - done;
        if tau.abs() > rem.abs() || tau * rem < 0.0 {
            tau = rem;
        }
        let beta = w.iter().map(|x| x * x).sum::<f64>().sqrt();
        if beta == 0.0 {
            return Some(w);
        }
        let (basis, h, fnorm) = crate::kernels::eigs::arnoldi(&apply, &w, m);
        let mut th = h.clone();
        for x in &mut th.data {
            *x *= tau;
        }
        let e = expm(&th).ok()?;
        let est = beta * fnorm * e.at(m - 1, 0).abs();
        if est > tol * (tau / t).abs().max(1e-300) && tau.abs() > 1e-12 * t.abs() && fnorm > 0.0 {
            tau *= 0.5;
            continue;
        }
        let mut nw = vec![0.0; n];
        for j in 0..m {
            let c = beta * e.at(j, 0);
            for (x, b) in nw.iter_mut().zip(&basis[j]) {
                *x += c * b;
            }
        }
        w = nw;
        done += tau;
        if fnorm == 0.0 {
            // invariant subspace: exact for any step
            continue;
        }
        tau *= 1.5;
    }
    Some(w)
}
