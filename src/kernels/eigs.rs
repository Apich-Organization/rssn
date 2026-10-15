//! # Sparse eigensolvers: Lanczos and Arnoldi with restarts
//!
//! Matrix-free eigenvalue solvers for a few eigenpairs of large operators:
//!
//! * [`eigs_sym`]: thick-restart Lanczos (Krylov-Schur form) with full
//!   reorthogonalisation for symmetric operators;
//! * [`eigs`]: implicitly restarted Arnoldi (IRAM, exact shifts applied
//!   with real single and double shifts so that conjugate pairs stay
//!   together) for general real operators; eigenvalues may be complex and
//!   are reported with the real and imaginary parts of the Ritz vectors;
//! * [`eigs_sym_shift_invert`] / [`eigs_shift_invert`]: spectral
//!   transformation `theta = 1 / (lambda - sigma)`; the caller supplies the
//!   action of `(A - sigma I)^{-1}` (for example through
//!   [`crate::kernels::krylov::Csr::sparse_lu`]) and the eigenvalues closest
//!   to `sigma` are returned;
//! * [`arnoldi`]: the plain Arnoldi factorisation `A V = V H + f e_m^T`.
//!
//! The operator closure writes `A x` into its second argument.

use num_complex::Complex64;

use crate::kernels::dense::{self, Mat};
use crate::kernels::random::Rng;

/// Which part of the spectrum is wanted.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Which {
    /// Eigenvalues of largest modulus.
    LargestMagnitude,
    /// Eigenvalues of smallest modulus (slow without shift-invert).
    SmallestMagnitude,
    /// Largest real part (largest algebraic for symmetric operators).
    LargestReal,
    /// Smallest real part (smallest algebraic for symmetric operators).
    SmallestReal,
    /// Largest absolute imaginary part.
    LargestImag,
}

/// Parameters of the eigensolvers.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EigsOptions {
    /// Krylov subspace size (`0` chooses `max(2k + 1, 20)`, capped at `n`).
    pub ncv: usize,
    /// Relative tolerance on the Ritz residuals.
    pub tol: f64,
    /// Maximal number of restarts.
    pub max_restarts: usize,
    /// Seed of the pseudo-random start vector.
    pub seed: u64,
}

impl Default for EigsOptions {
    fn default() -> Self {
        Self { ncv: 0, tol: 1e-10, max_restarts: 500, seed: 0x5eed }
    }
}

/// Errors of the eigensolvers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EigsError {
    /// Invalid dimensions or `k` out of range.
    Invalid,
    /// The small dense eigenproblem failed.
    Dense,
}

/// Result of an eigensolve.
#[derive(Debug, Clone, PartialEq)]
pub struct EigsResult {
    /// Eigenvalues `(re, im)`, most wanted first.
    pub values: Vec<(f64, f64)>,
    /// Real parts of the eigenvectors (unit 2-norm of the complex vector).
    pub vectors_re: Vec<Vec<f64>>,
    /// Imaginary parts of the eigenvectors (zero for real eigenvalues).
    pub vectors_im: Vec<Vec<f64>>,
    /// Residual estimates `|A x - lambda x|`.
    pub residuals: Vec<f64>,
    /// Whether all requested pairs met the tolerance.
    pub converged: bool,
    /// Number of restarts performed.
    pub restarts: usize,
    /// Number of operator applications.
    pub matvecs: usize,
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn nrm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}

/// Orthogonalises `w` against `v[..cnt]` (classical Gram-Schmidt twice) and
/// returns the coefficients.
fn orthogonalise(w: &mut [f64], v: &[Vec<f64>], cnt: usize) -> Vec<f64> {
    let mut h = vec![0.0; cnt];
    for _ in 0..2 {
        for i in 0..cnt {
            let c = dot(w, &v[i]);
            h[i] += c;
            for (wq, vq) in w.iter_mut().zip(&v[i]) {
                *wq -= c * vq;
            }
        }
    }
    h
}

fn random_unit(n: usize, seed: u64) -> Vec<f64> {
    let mut rng = Rng::new(seed);
    let mut v: Vec<f64> = (0..n).map(|_| rng.normal()).collect();
    let s = nrm(&v);
    for e in &mut v {
        *e /= s;
    }
    v
}

/// Arnoldi factorisation of length `m` started at `v0` (normalised
/// internally): returns `(V, H, beta)` with `V` holding `m + 1` basis
/// vectors (the last one is the normalised residual direction, zero after a
/// lucky breakdown), `H` the `m x m` Hessenberg matrix and `beta` the
/// residual norm, so that `A V_m = V_m H + beta v_{m+1} e_m^T`.
pub fn arnoldi<A: Fn(&[f64], &mut [f64])>(
    apply: A,
    v0: &[f64],
    m: usize,
) -> (Vec<Vec<f64>>, Mat, f64) {
    let n = v0.len();
    let s = nrm(v0).max(f64::MIN_POSITIVE);
    let mut v: Vec<Vec<f64>> = vec![v0.iter().map(|e| e / s).collect()];
    let mut h = Mat::zeros(m, m);
    let mut beta = 0.0;
    let mut w = vec![0.0; n];
    for j in 0..m {
        apply(&v[j], &mut w);
        let c = orthogonalise(&mut w, &v, j + 1);
        for i in 0..=j {
            h.set(i, j, c[i]);
        }
        beta = nrm(&w);
        if j + 1 < m {
            h.set(j + 1, j, beta);
        }
        if beta > 1e-14 * (1.0 + c.iter().map(|x| x.abs()).sum::<f64>()) {
            v.push(w.iter().map(|e| e / beta).collect());
        } else {
            beta = 0.0;
            v.push(vec![0.0; n]);
            if j + 1 < m {
                h.set(j + 1, j, 0.0);
                // continue with an arbitrary orthogonal direction would be needed;
                // stop the factorisation here (invariant subspace found)
                let mut hh = Mat::zeros(m, m);
                for a in 0..m {
                    for b in 0..m {
                        hh.set(a, b, h.at(a, b));
                    }
                }
                while v.len() < m + 1 {
                    v.push(vec![0.0; n]);
                }
                return (v, hh, 0.0);
            }
        }
    }
    (v, h, beta)
}

fn key_of(which: Which, v: (f64, f64)) -> f64 {
    // smaller key = more wanted
    match which {
        Which::LargestMagnitude => -v.0.hypot(v.1),
        Which::SmallestMagnitude => v.0.hypot(v.1),
        Which::LargestReal => -v.0,
        Which::SmallestReal => v.0,
        Which::LargestImag => -v.1.abs(),
    }
}

/// Unit-norm eigenvector of the small matrix `h` for eigenvalue `theta` by
/// complex inverse iteration.
fn small_eigvec(h: &Mat, theta: (f64, f64)) -> Vec<Complex64> {
    let n = h.rows;
    let scale = 1.0 + theta.0.hypot(theta.1);
    let shift = Complex64::new(theta.0 + 1e-10 * scale, theta.1 + 1e-10 * scale * 0.0);
    let mut a = vec![vec![Complex64::new(0.0, 0.0); n]; n];
    for i in 0..n {
        for j in 0..n {
            a[i][j] = Complex64::new(h.at(i, j), 0.0);
        }
        a[i][i] -= shift;
    }
    // LU with partial pivoting
    let mut piv: Vec<usize> = (0..n).collect();
    for k in 0..n {
        let mut p = k;
        for i in k + 1..n {
            if a[i][k].norm() > a[p][k].norm() {
                p = i;
            }
        }
        a.swap(k, p);
        piv.swap(k, p);
        if a[k][k].norm() < 1e-300 {
            a[k][k] = Complex64::new(1e-300, 0.0);
        }
        for i in k + 1..n {
            let f = a[i][k] / a[k][k];
            a[i][k] = f;
            for j in k + 1..n {
                let t = a[k][j];
                a[i][j] -= f * t;
            }
        }
    }
    let mut x: Vec<Complex64> = (0..n).map(|i| Complex64::new(1.0 + 0.1 * i as f64, 0.3)).collect();
    for _ in 0..3 {
        let mut y: Vec<Complex64> = piv.iter().map(|&p| x[p]).collect();
        for i in 0..n {
            for j in 0..i {
                let t = y[j];
                y[i] -= a[i][j] * t;
            }
        }
        for i in (0..n).rev() {
            for j in i + 1..n {
                let t = y[j];
                y[i] -= a[i][j] * t;
            }
            y[i] /= a[i][i];
        }
        let s = y.iter().map(Complex64::norm_sqr).sum::<f64>().sqrt();
        x = y.iter().map(|c| c / s).collect();
    }
    x
}

/// Thick-restart Lanczos for the `k` wanted eigenpairs of a symmetric
/// operator on `R^n`.
///
/// # Errors
/// [`EigsError::Invalid`] for `k == 0` or `k > n`; [`EigsError::Dense`] if
/// the projected eigenproblem fails.
pub fn eigs_sym<A: Fn(&[f64], &mut [f64])>(
    apply: A,
    n: usize,
    k: usize,
    which: Which,
    opts: &EigsOptions,
) -> Result<EigsResult, EigsError> {
    if k == 0 || k > n {
        return Err(EigsError::Invalid);
    }
    let m = (if opts.ncv == 0 { (2 * k + 1).max(20) } else { opts.ncv }).max(k + 1).min(n);
    let mut v: Vec<Vec<f64>> = vec![random_unit(n, opts.seed)];
    // projected matrix T = V^T A V (symmetric), entries filled column by column
    let mut t = Mat::zeros(m, m);
    let mut w = vec![0.0; n];
    let mut matvecs = 0;
    let mut kept = 0usize; // number of locked/thick vectors at the front
    let mut restarts = 0;
    loop {
        // expand from `kept` to m
        let mut beta = 0.0;
        for j in kept..m {
            apply(&v[j], &mut w);
            matvecs += 1;
            let c = orthogonalise(&mut w, &v, j + 1);
            for i in 0..=j {
                t.set(i, j, c[i]);
                t.set(j, i, c[i]);
            }
            beta = nrm(&w);
            if j + 1 < m {
                t.set(j + 1, j, beta);
                t.set(j, j + 1, beta);
            }
            if beta > 1e-13 {
                v.push(w.iter().map(|e| e / beta).collect());
            } else if j + 1 >= m {
                beta = 0.0;
                v.push(vec![0.0; n]);
            } else {
                // invariant subspace: restart with a fresh orthogonal direction
                let mut r = random_unit(n, opts.seed.wrapping_add(j as u64 + 17));
                let _ = orthogonalise(&mut r, &v, j + 1);
                let rn = nrm(&r);
                for e in &mut r {
                    *e /= rn;
                }
                beta = 0.0;
                if j + 1 < m {
                    t.set(j + 1, j, 0.0);
                    t.set(j, j + 1, 0.0);
                }
                v.push(r);
            }
        }
        let (theta, s) = dense::eigen_symmetric(&t).map_err(|_| EigsError::Dense)?;
        let mut idx: Vec<usize> = (0..m).collect();
        idx.sort_by(|&a, &b| {
            key_of(which, (theta[a], 0.0)).total_cmp(&key_of(which, (theta[b], 0.0)))
        });
        let res: Vec<f64> = idx.iter().map(|&i| beta * s.at(m - 1, i).abs()).collect();
        let all_conv = (0..k).all(|q| res[q] <= opts.tol * theta[idx[q]].abs().max(1e-12));
        if all_conv || restarts >= opts.max_restarts || m >= n {
            let mut out = EigsResult {
                values: Vec::new(),
                vectors_re: Vec::new(),
                vectors_im: Vec::new(),
                residuals: Vec::new(),
                converged: all_conv || m >= n,
                restarts,
                matvecs,
            };
            for q in 0..k {
                let col = idx[q];
                let mut x = vec![0.0; n];
                for i in 0..m {
                    let c = s.at(i, col);
                    for (xq, vq) in x.iter_mut().zip(&v[i]) {
                        *xq += c * vq;
                    }
                }
                out.values.push((theta[col], 0.0));
                out.vectors_re.push(x);
                out.vectors_im.push(vec![0.0; n]);
                out.residuals.push(if m >= n { 0.0 } else { res[q] });
            }
            return Ok(out);
        }
        // thick restart: keep the wanted Ritz pairs plus a few extra for faster convergence
        let keep = (k + (m - k) / 2).min(m - 1).max(k);
        let vnext = v[m].clone();
        let mut nv: Vec<Vec<f64>> = Vec::with_capacity(keep + 1);
        for q in 0..keep {
            let col = idx[q];
            let mut x = vec![0.0; n];
            for i in 0..m {
                let c = s.at(i, col);
                for (xq, vq) in x.iter_mut().zip(&v[i]) {
                    *xq += c * vq;
                }
            }
            nv.push(x);
        }
        nv.push(vnext);
        v = nv;
        t = Mat::zeros(m, m);
        for q in 0..keep {
            t.set(q, q, theta[idx[q]]);
            let c = beta * s.at(m - 1, idx[q]);
            t.set(keep, q, c);
            t.set(q, keep, c);
        }
        kept = keep;
        restarts += 1;
    }
}

/// Implicitly restarted Arnoldi for the `k` wanted eigenpairs of a general
/// real operator.
///
/// # Errors
/// [`EigsError::Invalid`] for `k == 0` or `k >= n`; [`EigsError::Dense`] if
/// the small dense eigenproblem fails.
pub fn eigs<A: Fn(&[f64], &mut [f64])>(
    apply: A,
    n: usize,
    k: usize,
    which: Which,
    opts: &EigsOptions,
) -> Result<EigsResult, EigsError> {
    if k == 0 || k > n {
        return Err(EigsError::Invalid);
    }
    let m = (if opts.ncv == 0 { (2 * k + 2).max(24) } else { opts.ncv }).max(k + 2).min(n);
    let mut matvecs = 0;
    let mut restarts = 0;
    let mut w = vec![0.0; n];
    // Arnoldi data: v has m+1 vectors, h is m x m (Hessenberg), f_norm = beta
    let mut v: Vec<Vec<f64>> = vec![random_unit(n, opts.seed)];
    let mut h = Mat::zeros(m, m);
    let mut start = 0usize;
    let mut beta;
    loop {
        // expand the factorisation from `start` to m
        beta = 0.0;
        for j in start..m {
            apply(&v[j], &mut w);
            matvecs += 1;
            let c = orthogonalise(&mut w, &v, j + 1);
            for i in 0..=j {
                h.set(i, j, c[i]);
            }
            beta = nrm(&w);
            let big = 1e-13 * (1.0 + c.iter().map(|x| x.abs()).fold(0.0, f64::max));
            if beta > big {
                v.push(w.iter().map(|e| e / beta).collect());
            } else if j + 1 >= m {
                beta = 0.0;
                v.push(vec![0.0; n]);
            } else {
                beta = 0.0;
                let mut r = random_unit(n, opts.seed.wrapping_add(j as u64 + 31));
                let _ = orthogonalise(&mut r, &v, j + 1);
                let rn = nrm(&r);
                for e in &mut r {
                    *e /= rn;
                }
                v.push(r);
            }
            if j + 1 < m {
                h.set(j + 1, j, beta);
            }
        }
        let ritz = dense::eigenvalues(&h).map_err(|_| EigsError::Dense)?;
        let mut order: Vec<usize> = (0..m).collect();
        order.sort_by(|&a, &b| key_of(which, ritz[a]).total_cmp(&key_of(which, ritz[b])));
        // do not split a conjugate pair
        let mut kk = k;
        let pair = |a: (f64, f64), b: (f64, f64)| {
            a.1 != 0.0
                && (a.0 - b.0).abs() <= 1e-10 * (1.0 + a.0.abs())
                && (a.1 + b.1).abs() <= 1e-10 * (1.0 + a.1.abs())
        };
        if kk < m && pair(ritz[order[kk - 1]], ritz[order[kk]]) {
            kk += 1;
        }
        // residual estimates of the wanted Ritz pairs
        let mut vecs: Vec<Vec<Complex64>> = Vec::new();
        let mut res: Vec<f64> = Vec::new();
        for q in 0..kk.min(m) {
            let y = small_eigvec(&h, ritz[order[q]]);
            res.push(beta * y[m - 1].norm());
            vecs.push(y);
        }
        let conv_count =
            (0..k).filter(|&q| res[q] <= opts.tol * ritz[order[q]].0.hypot(ritz[order[q]].1).max(1e-12)).count();
        let all_conv = conv_count == k;
        if all_conv || restarts >= opts.max_restarts || m >= n {
            let mut out = EigsResult {
                values: Vec::new(),
                vectors_re: Vec::new(),
                vectors_im: Vec::new(),
                residuals: Vec::new(),
                converged: all_conv || m >= n,
                restarts,
                matvecs,
            };
            for q in 0..kk.min(m) {
                let y = &vecs[q];
                let (mut xr, mut xi) = (vec![0.0; n], vec![0.0; n]);
                for i in 0..m {
                    for (e, vq) in v[i].iter().enumerate() {
                        xr[e] += y[i].re * vq;
                        xi[e] += y[i].im * vq;
                    }
                }
                // rotate the phase so that the largest component is real
                let mut big = 0usize;
                let mut bigv = 0.0;
                for e in 0..n {
                    let a = xr[e].hypot(xi[e]);
                    if a > bigv {
                        bigv = a;
                        big = e;
                    }
                }
                if bigv > 0.0 && ritz[order[q]].1 != 0.0 {
                    let ph = Complex64::new(xr[big], xi[big]).conj() / bigv;
                    for e in 0..n {
                        let c = Complex64::new(xr[e], xi[e]) * ph;
                        xr[e] = c.re;
                        xi[e] = c.im;
                    }
                } else {
                    let sgn = if xr[big] < 0.0 { -1.0 } else { 1.0 };
                    for e in 0..n {
                        xr[e] *= sgn;
                        xi[e] *= sgn;
                    }
                    if ritz[order[q]].1 == 0.0 {
                        xi.fill(0.0);
                    }
                }
                out.values.push(ritz[order[q]]);
                out.vectors_re.push(xr);
                out.vectors_im.push(xi);
                out.residuals.push(if m >= n { 0.0 } else { res[q] });
            }
            return Ok(out);
        }
        // implicit restart with the unwanted Ritz values as exact shifts
        let mut real_shifts: Vec<f64> = Vec::new();
        let mut pair_shifts: Vec<(f64, f64)> = Vec::new();
        for &i in &order[kk..] {
            let r = ritz[i];
            if r.1 == 0.0 {
                real_shifts.push(r.0);
            } else if r.1 > 0.0 {
                pair_shifts.push(r);
            }
        }
        let mut shifts: Vec<(f64, f64, bool)> = Vec::new();
        for r in real_shifts {
            shifts.push((r, 0.0, false));
        }
        for p in pair_shifts {
            shifts.push((p.0, p.1, true));
        }
        let mut q_tot = Mat::identity(m);
        let mut hh = h.clone();
        for &(mu_re, mu_im, is_pair) in &shifts {
            let mut mm;
            if is_pair {
                let h2 = hh.matmul(&hh).map_err(|_| EigsError::Dense)?;
                mm = Mat::zeros(m, m);
                for a in 0..m {
                    for b in 0..m {
                        let id = if a == b { 1.0 } else { 0.0 };
                        mm.set(a, b, h2.at(a, b) - 2.0 * mu_re * hh.at(a, b) + (mu_re * mu_re + mu_im * mu_im) * id);
                    }
                }
            } else {
                mm = hh.clone();
                for a in 0..m {
                    mm.set(a, a, mm.at(a, a) - mu_re);
                }
            }
            let (q, _) = dense::qr(&mm);
            let qt = q.transpose();
            hh = qt.matmul(&hh).and_then(|x| x.matmul(&q)).map_err(|_| EigsError::Dense)?;
            q_tot = q_tot.matmul(&q).map_err(|_| EigsError::Dense)?;
            for a in 2..m {
                for b in 0..a - 1 {
                    hh.set(a, b, 0.0);
                }
            }
        }
        // truncate to kk: V_kk = V Q[:, :kk]; new residual direction
        let mut nv: Vec<Vec<f64>> = Vec::with_capacity(kk + 1);
        for c in 0..kk {
            let mut x = vec![0.0; n];
            for i in 0..m {
                let qc = q_tot.at(i, c);
                if qc != 0.0 {
                    for (xe, ve) in x.iter_mut().zip(&v[i]) {
                        *xe += qc * ve;
                    }
                }
            }
            nv.push(x);
        }
        let hk = hh.at(kk, kk - 1);
        let sig = q_tot.at(m - 1, kk - 1);
        let mut f = vec![0.0; n];
        for i in 0..m {
            let qc = q_tot.at(i, kk);
            for (fe, ve) in f.iter_mut().zip(&v[i]) {
                *fe += qc * hk * ve;
            }
        }
        for (fe, ve) in f.iter_mut().zip(&v[m]) {
            *fe += beta * sig * ve;
        }
        let _ = orthogonalise(&mut f, &nv, kk);
        let fb = nrm(&f);
        let vnew = if fb > 1e-14 {
            f.iter().map(|e| e / fb).collect()
        } else {
            let mut r = random_unit(n, opts.seed.wrapping_add(restarts as u64 + 101));
            let _ = orthogonalise(&mut r, &nv, kk);
            let rn = nrm(&r);
            r.iter().map(|e| e / rn).collect()
        };
        nv.push(vnew);
        v = nv;
        let mut hn = Mat::zeros(m, m);
        for a in 0..kk {
            for b in 0..kk {
                hn.set(a, b, hh.at(a, b));
            }
        }
        // the coupling of the new residual direction to column kk-1
        // is carried by beta_k = ||f|| (subdiagonal of the extended matrix)
        h = hn;
        // store beta_k in H[kk, kk-1] for the extension (kk < m by construction)
        h.set(kk, kk - 1, fb);
        start = kk;
        restarts += 1;
    }
}

/// Symmetric shift-invert eigensolver: the eigenvalues of `A` closest to
/// `sigma`, given the action `solve(x, y)` of `(A - sigma I)^{-1}`.
///
/// # Errors
/// As [`eigs_sym`].
pub fn eigs_sym_shift_invert<S: Fn(&[f64], &mut [f64])>(
    solve: S,
    n: usize,
    k: usize,
    sigma: f64,
    opts: &EigsOptions,
) -> Result<EigsResult, EigsError> {
    let mut r = eigs_sym(solve, n, k, Which::LargestMagnitude, opts)?;
    for (val, res) in r.values.iter_mut().zip(&mut r.residuals) {
        let th = val.0;
        *res /= th * th;
        *val = (sigma + 1.0 / th, 0.0);
    }
    Ok(r)
}

/// General shift-invert eigensolver (see [`eigs_sym_shift_invert`]).
///
/// # Errors
/// As [`eigs`].
pub fn eigs_shift_invert<S: Fn(&[f64], &mut [f64])>(
    solve: S,
    n: usize,
    k: usize,
    sigma: f64,
    opts: &EigsOptions,
) -> Result<EigsResult, EigsError> {
    let mut r = eigs(solve, n, k, Which::LargestMagnitude, opts)?;
    for (val, res) in r.values.iter_mut().zip(&mut r.residuals) {
        let th = Complex64::new(val.0, val.1);
        *res /= th.norm_sqr();
        let lam = Complex64::new(sigma, 0.0) + Complex64::new(1.0, 0.0) / th;
        *val = (lam.re, lam.im);
    }
    Ok(r)
}
