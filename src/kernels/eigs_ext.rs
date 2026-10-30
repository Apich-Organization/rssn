//! # Robust sparse eigensolvers: Krylov-Schur, block Lanczos, Jacobi-Davidson
//!
//! Extensions of [`crate::kernels::eigs`] for harder problems:
//!
//! * [`eigs_krylov_schur`]: Krylov-Schur (Stewart) for **non-symmetric**
//!   operators. The Arnoldi factorisation is contracted by reordering the
//!   complex Schur form of the projected matrix, so wanted Ritz pairs
//!   (including complex-conjugate pairs) are kept exactly instead of through
//!   implicit shifts; this is the numerically robust restart for strongly
//!   non-normal matrices, where exact-shift IRAM can stagnate. The small
//!   problem is a complex Schur decomposition ([`complex_schur`]) with
//!   Wilkinson-shifted QR and Givens swaps ([`reorder_schur`]). The symmetric
//!   counterpart is the thick-restart Lanczos method of
//!   [`crate::kernels::eigs::eigs_sym`] (the symmetric Krylov-Schur form);
//! * [`eigs_sym_block`]: block Lanczos with full reorthogonalisation and
//!   thick restart (Ritz vectors plus the residual block). A block as large
//!   as the cluster multiplicity finds every member of tight clusters and
//!   multiple eigenvalues, which a single-vector Lanczos process can only
//!   resolve one rounding error at a time;
//! * [`eigs_jd_sym`] / [`eigs_jd`]: Jacobi-Davidson for **interior** eigenvalues of a
//!   symmetric operator with harmonic Ritz extraction, locking/deflation and
//!   a projected GMRES solution of the correction equation (optionally
//!   preconditioned);
//! * generalised problems `A x = lambda B x`: [`eigs_gen_sym_shift_invert`]
//!   (symmetric `A`, symmetric positive definite `B`; block Lanczos in the
//!   `B` inner product on `(A - sigma B)^{-1} B`) and
//!   [`eigs_gen_shift_invert`] (any `A`, nonsingular `B`; Krylov-Schur).
//!   The caller supplies the action of `(A - sigma B)^{-1}`, for example
//!   through [`crate::kernels::krylov::Csr::sparse_lu`].
//!
//! After `k` eigenvalues are locked the Jacobi-Davidson solver keeps iterating
//! until the best remaining candidate is no closer to the target than the
//! `k`-th locked eigenvalue, so that the *nearest* eigenvalues are returned.
//!
//! [`eigs_jd`] extends Jacobi-Davidson to non-symmetric operators (complex
//! targets, partial Schur deflation, complex harmonic extraction).
//!
//! [`eigs_jd_two_sided`] is the *two-sided* variant: it expands a right and a
//! left search space together, extracts Petrov-Galerkin pairs of the projected
//! pencil, returns left and right eigenvectors with the eigenvalue condition
//! numbers `1/|w^H u|`, uses the (stationary) two-sided Rayleigh quotient, solves
//! two oblique correction equations and deflates converged pairs obliquely
//! (`I - Q (P^H Q)^{-1} P^H`). Its extraction is the standard one (no harmonic
//! variant), so the target should be near the wanted eigenvalues.
//!
//! Limitations: [`eigs_jd_sym`] is for symmetric operators only; `B` must be
//! positive definite in the symmetric generalised solver; the two-sided solver
//! needs the transpose action `A^T` and assumes diagonalisable (non-defective)
//! wanted eigenvalues.

use num_complex::Complex64;

use crate::kernels::dense::Mat;
use crate::kernels::dense::{
    self,
};
use crate::kernels::eigs::EigsError;
use crate::kernels::eigs::EigsOptions;
use crate::kernels::eigs::EigsResult;
use crate::kernels::eigs::Which;
use crate::kernels::eigs::key_of;
use crate::kernels::krylov::gmres;
use crate::kernels::random::Rng;

type C = Complex64;
type Schur = (Vec<Vec<C>>, Vec<Vec<C>>);
type OpFn<'a> = dyn Fn(&[f64], &mut [f64]) + 'a;

fn rdot(
    a: &[f64],
    b: &[f64],
) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn rnrm(a: &[f64]) -> f64 {
    rdot(a, a).sqrt()
}

fn cdot(
    a: &[C],
    b: &[C],
) -> C {
    a.iter().zip(b).map(|(x, y)| x.conj() * y).sum()
}

fn cnrm(a: &[C]) -> f64 {
    a.iter().map(Complex64::norm_sqr).sum::<f64>().sqrt()
}

fn random_vec(
    n: usize,
    rng: &mut Rng,
) -> Vec<f64> {
    (0..n).map(|_| rng.normal()).collect()
}

// ------------------------------------------------------- complex Schur form

/// Complex Givens rotation `[[c, s], [-conj(s), c]]` (`c` real) mapping
/// `(a, b)` to `(r, 0)`.
fn givens(
    a: C,
    b: C,
) -> (f64, C) {
    let (na, nb) = (a.norm(), b.norm());
    if nb == 0.0 {
        return (1.0, C::new(0.0, 0.0));
    }
    let r = na.hypot(nb);
    if na == 0.0 {
        return (0.0, b.conj() / nb);
    }
    (na / r, (a / na) * b.conj() / r)
}

/// `rows k, k+1 <- G rows`, columns `c0..ncols`.
fn rot_rows(
    h: &mut [Vec<C>],
    k: usize,
    c: f64,
    s: C,
    c0: usize,
) {
    let ncols = h[k].len();
    for j in c0..ncols {
        let (x, y) = (h[k][j], h[k + 1][j]);
        h[k][j] = x * c + s * y;
        h[k + 1][j] = -s.conj() * x + y * c;
    }
}

/// `columns k, k+1 <- columns * G^H`, rows `0..r1`.
fn rot_cols(
    h: &mut [Vec<C>],
    k: usize,
    c: f64,
    s: C,
    r1: usize,
) {
    for row in h.iter_mut().take(r1) {
        let (x, y) = (row[k], row[k + 1]);
        row[k] = x * c + y * s.conj();
        row[k + 1] = -x * s + y * c;
    }
}

/// Complex Schur decomposition `A = Q T Q^H` (`T` upper triangular, `Q`
/// unitary) of a dense square matrix by Householder reduction to Hessenberg
/// form and Wilkinson-shifted QR iterations with deflation.
///
/// Returns `None` if the iteration does not converge.
#[must_use]
pub fn complex_schur(a: &[Vec<C>]) -> Option<Schur> {
    let n = a.len();
    let mut h: Vec<Vec<C>> = a.to_vec();
    let mut q: Vec<Vec<C>> = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| {
                    if i == j {
                        C::new(1.0, 0.0)
                    } else {
                        C::new(0.0, 0.0)
                    }
                })
                .collect()
        })
        .collect();
    // Householder reduction to Hessenberg form
    for k in 0..n.saturating_sub(2) {
        let mut x: Vec<C> = (k + 1..n).map(|i| h[i][k]).collect();
        let xn = cnrm(&x);
        if xn == 0.0 {
            continue;
        }
        let phase = if x[0].norm() == 0.0 {
            C::new(1.0, 0.0)
        } else {
            x[0] / x[0].norm()
        };
        x[0] += phase * xn;
        let vn = cnrm(&x);
        if vn == 0.0 {
            continue;
        }
        for e in &mut x {
            *e /= vn;
        }
        // H <- P H P with P = I - 2 v v^H acting on indices k+1..n
        for j in 0..n {
            let s: C = (0..x.len()).map(|i| x[i].conj() * h[k + 1 + i][j]).sum();
            for i in 0..x.len() {
                h[k + 1 + i][j] -= 2.0 * x[i] * s;
            }
        }
        for row in &mut h {
            let s: C = (0..x.len()).map(|i| row[k + 1 + i] * x[i]).sum();
            for i in 0..x.len() {
                row[k + 1 + i] -= 2.0 * s * x[i].conj();
            }
        }
        for row in &mut q {
            let s: C = (0..x.len()).map(|i| row[k + 1 + i] * x[i]).sum();
            for i in 0..x.len() {
                row[k + 1 + i] -= 2.0 * s * x[i].conj();
            }
        }
        for i in k + 2..n {
            h[i][k] = C::new(0.0, 0.0);
        }
    }
    // shifted QR on the Hessenberg matrix
    let eps = f64::EPSILON;
    let mut hi = n.saturating_sub(1);
    let mut iter = 0usize;
    let mut total = 0usize;
    while hi > 0 {
        total += 1;
        if total > 100 * n.max(10) {
            return None;
        }
        let mut l = hi;
        while l > 0 {
            let s = h[l - 1][l - 1].norm() + h[l][l].norm();
            let s = if s == 0.0 { 1.0 } else { s };
            if h[l][l - 1].norm() <= eps * s {
                h[l][l - 1] = C::new(0.0, 0.0);
                break;
            }
            l -= 1;
        }
        if l == hi {
            hi -= 1;
            iter = 0;
            continue;
        }
        iter += 1;
        // Wilkinson shift from the trailing 2x2 block
        let (a11, a12, a21, a22) = (h[hi - 1][hi - 1], h[hi - 1][hi], h[hi][hi - 1], h[hi][hi]);
        let mu = if iter.is_multiple_of(10) {
            a22 + C::new(h[hi][hi - 1].norm(), 0.0)
        } else {
            let tr_half = (a11 + a22) * 0.5;
            let det = a11 * a22 - a12 * a21;
            let disc = (tr_half * tr_half - det).sqrt();
            let (l1, l2) = (tr_half + disc, tr_half - disc);
            if (l1 - a22).norm() < (l2 - a22).norm() {
                l1
            } else {
                l2
            }
        };
        for i in l..=hi {
            h[i][i] -= mu;
        }
        let mut rots: Vec<(f64, C)> = Vec::with_capacity(hi - l);
        for k in l..hi {
            let (c, s) = givens(h[k][k], h[k + 1][k]);
            rot_rows(&mut h, k, c, s, k);
            h[k + 1][k] = C::new(0.0, 0.0);
            rots.push((c, s));
        }
        for (idx, &(c, s)) in rots.iter().enumerate() {
            let k = l + idx;
            rot_cols(&mut h, k, c, s, (k + 2).min(hi + 1));
            rot_cols(&mut q, k, c, s, n);
        }
        for i in l..=hi {
            h[i][i] += mu;
        }
    }
    for (i, row) in h.iter_mut().enumerate() {
        for e in row.iter_mut().take(i) {
            *e = C::new(0.0, 0.0);
        }
    }
    Some((h, q))
}

/// Swaps the adjacent diagonal entries `i, i+1` of a Schur form by a Givens
/// rotation, updating `Q` accordingly.
fn swap_schur(
    t: &mut [Vec<C>],
    q: &mut [Vec<C>],
    i: usize,
) {
    let n = t.len();
    let (a, b, x) = (t[i][i], t[i + 1][i + 1], t[i][i + 1]);
    let (c, s) = givens(x, b - a);
    rot_rows(t, i, c, s, i);
    rot_cols(t, i, c, s, (i + 2).min(n));
    rot_cols(q, i, c, s, n);
    t[i + 1][i] = C::new(0.0, 0.0);
}

/// Reorders a complex Schur form so that its diagonal is sorted by `key`
/// (smallest first) for the first `count` positions; `Q` is updated so that
/// `A = Q T Q^H` still holds.
pub fn reorder_schur(
    t: &mut [Vec<C>],
    q: &mut [Vec<C>],
    count: usize,
    key: &dyn Fn(C) -> f64,
) {
    let n = t.len();
    for p in 0..count.min(n) {
        let mut best = p;
        for j in p + 1..n {
            if key(t[j][j]) < key(t[best][best]) {
                best = j;
            }
        }
        for i in (p..best).rev() {
            swap_schur(t, q, i);
        }
    }
}

// ----------------------------------------------------------- Krylov-Schur

fn apply_complex<A: Fn(&[f64], &mut [f64])>(
    apply: &A,
    v: &[C],
    out: &mut [C],
    re: &mut [f64],
    im: &mut [f64],
) {
    let n = v.len();
    let (mut xr, mut xi) = (vec![0.0; n], vec![0.0; n]);
    for i in 0..n {
        xr[i] = v[i].re;
        xi[i] = v[i].im;
    }
    apply(&xr, re);
    if xi.iter().any(|e| *e != 0.0) {
        apply(&xi, im);
    } else {
        im.fill(0.0);
    }
    for i in 0..n {
        out[i] = C::new(re[i], im[i]);
    }
}

fn orth_complex(
    w: &mut [C],
    v: &[Vec<C>],
    cnt: usize,
) -> Vec<C> {
    let mut h = vec![C::new(0.0, 0.0); cnt];
    for _ in 0..2 {
        for i in 0..cnt {
            let c = cdot(&v[i], w);
            h[i] += c;
            for (wq, vq) in w.iter_mut().zip(&v[i]) {
                *wq -= c * vq;
            }
        }
    }
    h
}

fn random_complex_orth(
    n: usize,
    v: &[Vec<C>],
    cnt: usize,
    seed: u64,
) -> Vec<C> {
    let mut rng = Rng::new(seed);
    let mut r: Vec<C> = (0..n).map(|_| C::new(rng.normal(), 0.0)).collect();
    let _ = orth_complex(&mut r, v, cnt);
    let rn = cnrm(&r);
    for e in &mut r {
        *e /= rn;
    }
    r
}

/// Krylov-Schur eigensolver for a general real operator: the `k` eigenpairs
/// selected by `which`.
///
/// The Arnoldi factorisation `A V = V B + beta v e^T` of dimension `m`
/// (`opts.ncv`, default `max(2k + 1, 20)`) is contracted by computing the
/// Schur form of `B`, moving the wanted eigenvalues to the leading block and
/// truncating; the factorisation is then extended again. Ritz residuals come
/// from the last row of the Schur vectors and are exact for the Arnoldi
/// relation. Conjugate pairs are kept together.
///
/// # Errors
/// [`EigsError::Invalid`] for `k == 0` or `k > n`; [`EigsError::Dense`] if
/// the Schur iteration fails.
pub fn eigs_krylov_schur<A: Fn(&[f64], &mut [f64])>(
    apply: A,
    n: usize,
    k: usize,
    which: Which,
    opts: &EigsOptions,
) -> Result<EigsResult, EigsError> {
    if k == 0 || k > n {
        return Err(EigsError::Invalid);
    }
    let m = (if opts.ncv == 0 {
        (2 * k + 1).max(20)
    } else {
        opts.ncv
    })
    .max(k + 2)
    .min(n);
    let zero = C::new(0.0, 0.0);
    let mut rng = Rng::new(opts.seed);
    let mut v0: Vec<C> = (0..n).map(|_| C::new(rng.normal(), 0.0)).collect();
    let s0 = cnrm(&v0);
    for e in &mut v0 {
        *e /= s0;
    }
    let mut v: Vec<Vec<C>> = vec![v0];
    let mut b: Vec<Vec<C>> = vec![vec![zero; m]; m];
    let mut kept = 0usize;
    let mut restarts = 0usize;
    let mut matvecs = 0usize;
    let (mut wr, mut wi) = (vec![0.0; n], vec![0.0; n]);
    let mut w = vec![zero; n];
    let key = |z: C| key_of(which, (z.re, z.im));
    loop {
        let mut beta = 0.0;
        for j in kept..m {
            apply_complex(&apply, &v[j], &mut w, &mut wr, &mut wi);
            matvecs += 1;
            let h = orth_complex(&mut w, &v, j + 1);
            for (i, hv) in h.iter().enumerate() {
                b[i][j] = *hv;
            }
            beta = cnrm(&w);
            if beta > 1e-13 {
                v.push(w.iter().map(|e| e / beta).collect());
                if j + 1 < m {
                    b[j + 1][j] = C::new(beta, 0.0);
                }
            } else {
                beta = 0.0;
                if j + 1 < m {
                    b[j + 1][j] = zero;
                    let r =
                        random_complex_orth(n, &v, j + 1, opts.seed.wrapping_add(j as u64 + 31));
                    v.push(r);
                } else {
                    v.push(vec![zero; n]);
                }
            }
        }
        let (mut t, mut q) = complex_schur(&b).ok_or(EigsError::Dense)?;
        reorder_schur(&mut t, &mut q, m, &key);
        // residuals of the leading Ritz pairs
        let kk_eval = (k + 1).min(m);
        let mut res = vec![0.0; kk_eval];
        let mut ys: Vec<Vec<C>> = Vec::with_capacity(kk_eval);
        for i in 0..kk_eval {
            let mut y = vec![zero; i + 1];
            y[i] = C::new(1.0, 0.0);
            let lam = t[i][i];
            for r in (0..i).rev() {
                let mut s = zero;
                for c in r + 1..=i {
                    s += t[r][c] * y[c];
                }
                let d = t[r][r] - lam;
                let d = if d.norm() < 1e-14 * lam.norm().max(1.0) {
                    C::new(1e-14 * lam.norm().max(1.0), 0.0)
                } else {
                    d
                };
                y[r] = -s / d;
            }
            let yn = cnrm(&y);
            let cy: C = (0..=i).map(|j| q[m - 1][j] * y[j]).sum();
            res[i] = beta * cy.norm() / yn;
            ys.push(y);
        }
        let conv = |i: usize| res[i] <= opts.tol * t[i][i].norm().max(1e-12);
        let mut kk_out = k;
        // do not split a conjugate pair at the boundary
        if k < m && t[k - 1][k - 1].im.abs() > 1e-12 * t[k - 1][k - 1].norm().max(1.0) {
            let d = (t[k - 1][k - 1] - t[k][k].conj()).norm();
            if d < 1e-8 * t[k][k].norm().max(1.0) {
                kk_out = k + 1;
            }
        }
        let all_conv = (0..k).all(&conv) && (kk_out == k || conv(k));
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
            for i in 0..k.min(kk_eval) {
                let y = &ys[i];
                let yn = cnrm(y);
                let mut x = vec![zero; n];
                for j in 0..=i {
                    // Schur vector j = V q[:, j]
                    let coef = y[j] / yn;
                    for c in 0..m {
                        let qc = q[c][j] * coef;
                        for (xe, ve) in x.iter_mut().zip(&v[c]) {
                            *xe += qc * ve;
                        }
                    }
                }
                let xn = cnrm(&x);
                let big = x
                    .iter()
                    .copied()
                    .max_by(|p, r| p.norm().total_cmp(&r.norm()))
                    .unwrap_or(zero);
                let phase = if big.norm() == 0.0 {
                    C::new(1.0, 0.0)
                } else {
                    big.conj() / big.norm()
                };
                out.values.push((t[i][i].re, t[i][i].im));
                out.vectors_re
                    .push(x.iter().map(|e| (e * phase / xn).re).collect());
                out.vectors_im
                    .push(x.iter().map(|e| (e * phase / xn).im).collect());
                out.residuals.push(if m >= n { 0.0 } else { res[i] });
            }
            return Ok(out);
        }
        // contraction to kk Schur vectors (never splitting a conjugate pair)
        let mut keep = (k + (m - k) / 2).clamp(k + 1, m - 1);
        if keep < m
            && t[keep - 1][keep - 1].im.abs() > 1e-12 * t[keep - 1][keep - 1].norm().max(1.0)
            && (t[keep - 1][keep - 1] - t[keep][keep].conj()).norm()
                < 1e-8 * t[keep][keep].norm().max(1.0)
        {
            keep += 1;
        }
        let keep = keep.min(m - 1);
        let vnext = v[m].clone();
        let mut nv: Vec<Vec<C>> = Vec::with_capacity(keep + 1);
        for i in 0..keep {
            let mut x = vec![zero; n];
            for c in 0..m {
                let qc = q[c][i];
                for (xe, ve) in x.iter_mut().zip(&v[c]) {
                    *xe += qc * ve;
                }
            }
            nv.push(x);
        }
        nv.push(vnext);
        v = nv;
        b = vec![vec![zero; m]; m];
        for i in 0..keep {
            for j in i..keep {
                b[i][j] = t[i][j];
            }
            b[keep][i] = q[m - 1][i] * beta;
        }
        kept = keep;
        restarts += 1;
    }
}

// ------------------------------------------------------------- block Lanczos

/// Basis of the block engine: `B`-orthonormal vectors `v`, `B v` and `op v`.
struct BlockState {
    v: Vec<Vec<f64>>,
    bv: Vec<Vec<f64>>,
    av: Vec<Vec<f64>>,
    matvecs: usize,
    rng: Rng,
}

/// Block Lanczos/Davidson engine in the (optional) `B` inner product.
///
/// `op` must be self-adjoint in that inner product.
struct BlockEngine<'a> {
    n: usize,
    op: &'a dyn Fn(&[f64], &mut [f64]),
    bmul: Option<&'a OpFn<'a>>,
    block: usize,
}

struct BlockOut {
    theta: Vec<f64>,
    vectors: Vec<Vec<f64>>,
    residuals: Vec<f64>,
    converged: bool,
    restarts: usize,
    matvecs: usize,
}

impl BlockEngine<'_> {
    fn b_apply(
        &self,
        x: &[f64],
    ) -> Vec<f64> {
        self.bmul.map_or_else(
            || x.to_vec(),
            |bm| {
                let mut t = vec![0.0; self.n];
                bm(x, &mut t);
                t
            },
        )
    }

    /// Adds the B-orthonormalised candidates to the basis (dependent ones
    /// are replaced by random directions); returns the number added.
    fn extend(
        &self,
        st: &mut BlockState,
        cands: Vec<Vec<f64>>,
        limit: usize,
    ) -> usize {
        let mut added = 0;
        for mut c in cands {
            if st.v.len() >= limit.min(self.n) {
                break;
            }
            for attempt in 0..4 {
                let bc = self.b_apply(&c);
                let before = rdot(&c, &bc).max(0.0).sqrt();
                // two passes of B-Gram-Schmidt
                let mut bw = bc;
                for _ in 0..2 {
                    for i in 0..st.v.len() {
                        let coef = rdot(&st.v[i], &bw);
                        for q in 0..self.n {
                            c[q] -= coef * st.v[i][q];
                            bw[q] -= coef * st.bv[i][q];
                        }
                    }
                }
                let bw = self.b_apply(&c);
                let after = rdot(&c, &bw).max(0.0).sqrt();
                if after > 1e-9 * before.max(1e-300) && after > 1e-14 {
                    for e in &mut c {
                        *e /= after;
                    }
                    let bnew: Vec<f64> = bw.iter().map(|e| e / after).collect();
                    let mut a = vec![0.0; self.n];
                    (self.op)(&c, &mut a);
                    st.matvecs += 1;
                    st.v.push(c);
                    st.bv.push(bnew);
                    st.av.push(a);
                    added += 1;
                    break;
                }
                if attempt == 3 {
                    break;
                }
                c = random_vec(self.n, &mut st.rng);
            }
        }
        added
    }

    /// Replaces the basis by the Ritz combinations `cols` of the first `m`
    /// vectors (columns of `s`).
    fn contract(
        &self,
        st: &mut BlockState,
        s: &Mat,
        cols: &[usize],
    ) {
        let n = self.n;
        let m = st.v.len();
        let (mut nv, mut nb, mut na) = (Vec::new(), Vec::new(), Vec::new());
        for &col in cols {
            let (mut x, mut bx, mut ax) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
            for i in 0..m {
                let c = s.at(i, col);
                for q in 0..n {
                    x[q] += c * st.v[i][q];
                    bx[q] += c * st.bv[i][q];
                    ax[q] += c * st.av[i][q];
                }
            }
            nv.push(x);
            nb.push(bx);
            na.push(ax);
        }
        st.v = nv;
        st.bv = nb;
        st.av = na;
    }

    fn run(
        &self,
        k: usize,
        which: Which,
        opts: &EigsOptions,
    ) -> Result<BlockOut, EigsError> {
        let n = self.n;
        let bs = self.block.max(1);
        let m_max = (if opts.ncv == 0 {
            (2 * k + 2 * bs).max(24).max(4 * bs)
        } else {
            opts.ncv
        })
        .max(k + 2 * bs)
        .min(n);
        let mut st = BlockState {
            v: Vec::new(),
            bv: Vec::new(),
            av: Vec::new(),
            matvecs: 0,
            rng: Rng::new(opts.seed),
        };
        let start: Vec<Vec<f64>> = (0..bs).map(|_| random_vec(n, &mut st.rng)).collect();
        let mut last_block_start = 0usize;
        let mut added = self.extend(&mut st, start, m_max);
        let mut restarts = 0usize;
        loop {
            // Krylov expansion by blocks until the basis is full
            while st.v.len() < m_max && added > 0 {
                let cands: Vec<Vec<f64>> = (last_block_start..st.v.len())
                    .map(|i| st.av[i].clone())
                    .collect();
                last_block_start = st.v.len();
                added = self.extend(&mut st, cands, m_max);
            }
            let m = st.v.len();
            let mut t = Mat::zeros(m, m);
            for i in 0..m {
                for j in 0..m {
                    t.set(i, j, rdot(&st.bv[i], &st.av[j]));
                }
            }
            for i in 0..m {
                for j in 0..i {
                    let s = f64::midpoint(t.at(i, j), t.at(j, i));
                    t.set(i, j, s);
                    t.set(j, i, s);
                }
            }
            let (theta, s) = dense::eigen_symmetric(&t).map_err(|_| EigsError::Dense)?;
            let mut idx: Vec<usize> = (0..m).collect();
            idx.sort_by(|&a, &b| {
                key_of(which, (theta[a], 0.0)).total_cmp(&key_of(which, (theta[b], 0.0)))
            });
            let kk = k.min(m);
            let nres = (kk + bs).min(m);
            let mut xs: Vec<Vec<f64>> = Vec::with_capacity(nres);
            let mut rs: Vec<Vec<f64>> = Vec::with_capacity(nres);
            let mut res = Vec::with_capacity(nres);
            for &col in idx.iter().take(nres) {
                let mut x = vec![0.0; n];
                let mut ax = vec![0.0; n];
                for i in 0..m {
                    let c = s.at(i, col);
                    for q in 0..n {
                        x[q] += c * st.v[i][q];
                        ax[q] += c * st.av[i][q];
                    }
                }
                let r: Vec<f64> = ax.iter().zip(&x).map(|(a, b)| a - theta[col] * b).collect();
                res.push(rnrm(&r));
                xs.push(x);
                rs.push(r);
            }
            let scale = |i: usize| opts.tol * theta[idx[i]].abs().max(1e-12);
            let all_conv = (0..kk).all(|i| res[i] <= scale(i));
            if all_conv || restarts >= opts.max_restarts || m >= n {
                return Ok(BlockOut {
                    theta: idx.iter().take(kk).map(|&i| theta[i]).collect(),
                    vectors: xs.into_iter().take(kk).collect(),
                    residuals: if m >= n {
                        vec![0.0; kk]
                    } else {
                        res.into_iter().take(kk).collect()
                    },
                    converged: all_conv || m >= n,
                    restarts,
                    matvecs: st.matvecs,
                });
            }
            // thick restart: keep Ritz vectors, continue with the residual block
            let keep = (k + (m_max - k) / 2)
                .min(m_max.saturating_sub(bs + 1))
                .max(k)
                .min(m - 1);
            self.contract(&mut st, &s, &idx[..keep]);
            let mut cands: Vec<Vec<f64>> = Vec::new();
            for i in 0..kk {
                if res[i] > scale(i) && cands.len() < bs {
                    cands.push(rs[i].clone());
                }
            }
            for r in rs.iter().skip(kk) {
                if cands.len() < bs {
                    cands.push(r.clone());
                }
            }
            last_block_start = st.v.len();
            added = self.extend(&mut st, cands, m_max);
            if added == 0 {
                last_block_start = st.v.len().saturating_sub(1);
                let c = vec![st.av[last_block_start].clone()];
                added = self.extend(&mut st, c, m_max);
            }
            restarts += 1;
        }
    }
}

fn pack(
    out: BlockOut,
    map: &dyn Fn(f64) -> f64,
    n: usize,
) -> EigsResult {
    EigsResult {
        values: out.theta.iter().map(|&t| (map(t), 0.0)).collect(),
        vectors_re: out.vectors,
        vectors_im: vec![vec![0.0; n]; out.theta.len()],
        residuals: out.residuals,
        converged: out.converged,
        restarts: out.restarts,
        matvecs: out.matvecs,
    }
}

/// Block Lanczos with thick restart for the `k` wanted eigenpairs of a
/// symmetric operator, using blocks of `block` vectors.
///
/// Choose `block` at least as large as the multiplicity or the number of
/// members of the tight clusters of the wanted part of the spectrum. The
/// basis is kept orthonormal (two passes of Gram-Schmidt); the projected
/// matrix is formed explicitly, so residuals are true residuals.
///
/// # Errors
/// [`EigsError::Invalid`] for `k == 0`, `k > n` or `block == 0`;
/// [`EigsError::Dense`] if the projected eigenproblem fails.
pub fn eigs_sym_block<A: Fn(&[f64], &mut [f64])>(
    apply: A,
    n: usize,
    k: usize,
    which: Which,
    block: usize,
    opts: &EigsOptions,
) -> Result<EigsResult, EigsError> {
    if k == 0 || k > n || block == 0 {
        return Err(EigsError::Invalid);
    }
    let eng = BlockEngine {
        n,
        op: &apply,
        bmul: None,
        block,
    };
    let out = eng.run(k, which, opts)?;
    Ok(pack(out, &|t| t, n))
}

// ------------------------------------------------- generalised eigenproblems

/// Options of the generalised shift-invert eigensolvers.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GenEigsOptions {
    /// Shift `sigma`; the eigenvalues closest to it are returned.
    pub sigma: f64,
    /// Block size of the symmetric (block Lanczos) solver.
    pub block: usize,
    /// Krylov and tolerance parameters.
    pub eigs: EigsOptions,
}

fn gen_residuals(
    apply_a: &dyn Fn(&[f64], &mut [f64]),
    apply_b: &dyn Fn(&[f64], &mut [f64]),
    r: &mut EigsResult,
) {
    let n = r.vectors_re.first().map_or(0, Vec::len);
    let (mut ax, mut bx) = (vec![0.0; n], vec![0.0; n]);
    for i in 0..r.values.len() {
        let (lr, li) = r.values[i];
        let (xr, xi) = (&r.vectors_re[i], &r.vectors_im[i]);
        // residual of A x - lambda B x with complex x, lambda
        let mut sq = 0.0;
        let mut parts = [vec![0.0; n], vec![0.0; n]];
        apply_a(xr, &mut ax);
        apply_b(xr, &mut bx);
        for q in 0..n {
            parts[0][q] = ax[q] - lr * bx[q];
            parts[1][q] = -li * bx[q];
        }
        // contribution of the imaginary part of x
        apply_a(xi, &mut ax);
        apply_b(xi, &mut bx);
        for q in 0..n {
            parts[0][q] += li * bx[q];
            parts[1][q] += ax[q] - lr * bx[q];
        }
        for p in &parts {
            sq += rdot(p, p);
        }
        r.residuals[i] = sq.sqrt();
    }
}

/// Symmetric generalised eigenproblem `A x = lambda B x` (`A` symmetric,
/// `B` symmetric positive definite): the `k` eigenvalues closest to
/// `opts.sigma`.
///
/// `solve(x, y)` writes `(A - sigma B)^{-1} x` into `y`. Block Lanczos runs
/// in the `B` inner product on `(A - sigma B)^{-1} B`, which is self-adjoint
/// there; `theta` maps back to `lambda = sigma + 1/theta`. The returned
/// eigenvectors are `B`-orthonormal and `residuals` are the true norms
/// `|A x - lambda B x|`.
///
/// # Errors
/// As [`eigs_sym_block`].
pub fn eigs_gen_sym_shift_invert<A, B, S>(
    apply_a: A,
    apply_b: B,
    solve: S,
    n: usize,
    k: usize,
    opts: &GenEigsOptions,
) -> Result<EigsResult, EigsError>
where
    A: Fn(&[f64], &mut [f64]),
    B: Fn(&[f64], &mut [f64]),
    S: Fn(&[f64], &mut [f64]),
{
    if k == 0 || k > n || opts.block == 0 {
        return Err(EigsError::Invalid);
    }
    let op = |x: &[f64], y: &mut [f64]| {
        let mut bx = vec![0.0; x.len()];
        apply_b(x, &mut bx);
        solve(&bx, y);
    };
    let eng = BlockEngine {
        n,
        op: &op,
        bmul: Some(&apply_b),
        block: opts.block,
    };
    let out = eng.run(k, Which::LargestMagnitude, &opts.eigs)?;
    let sigma = opts.sigma;
    let mut r = pack(out, &|t| sigma + 1.0 / t, n);
    gen_residuals(&apply_a, &apply_b, &mut r);
    Ok(r)
}

/// General generalised eigenproblem `A x = lambda B x`.
///
/// `A` is an arbitrary real matrix and `B` real and nonsingular: the `k`
/// eigenvalues closest to `opts.sigma` by Krylov-Schur on
/// `(A - sigma B)^{-1} B`. Eigenvalues may be complex.
///
/// # Errors
/// As [`eigs_krylov_schur`].
pub fn eigs_gen_shift_invert<A, B, S>(
    apply_a: A,
    apply_b: B,
    solve: S,
    n: usize,
    k: usize,
    opts: &GenEigsOptions,
) -> Result<EigsResult, EigsError>
where
    A: Fn(&[f64], &mut [f64]),
    B: Fn(&[f64], &mut [f64]),
    S: Fn(&[f64], &mut [f64]),
{
    let op = |x: &[f64], y: &mut [f64]| {
        let mut bx = vec![0.0; x.len()];
        apply_b(x, &mut bx);
        solve(&bx, y);
    };
    let mut r = eigs_krylov_schur(op, n, k, Which::LargestMagnitude, &opts.eigs)?;
    for v in &mut r.values {
        let lam = C::new(opts.sigma, 0.0) + C::new(1.0, 0.0) / C::new(v.0, v.1);
        *v = (lam.re, lam.im);
    }
    gen_residuals(&apply_a, &apply_b, &mut r);
    Ok(r)
}

// ----------------------------------------------------------- Jacobi-Davidson

/// Parameters of [`eigs_jd_sym`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct JdOptions {
    /// Residual tolerance `|A x - theta x| <= tol * max(|theta|, 1)`.
    pub tol: f64,
    /// Maximal number of outer iterations.
    pub max_iter: usize,
    /// Dimension the search space is reduced to at a restart.
    pub min_dim: usize,
    /// Maximal dimension of the search space.
    pub max_dim: usize,
    /// GMRES iterations per correction equation.
    pub inner_iters: usize,
    /// Relative tolerance of the correction-equation solve.
    pub inner_tol: f64,
    /// Use harmonic Ritz extraction (recommended for interior eigenvalues).
    pub harmonic: bool,
    /// Seed of the start vector.
    pub seed: u64,
}

impl Default for JdOptions {
    fn default() -> Self {
        Self {
            tol: 1e-9,
            max_iter: 500,
            min_dim: 6,
            max_dim: 20,
            inner_iters: 30,
            inner_tol: 1e-3,
            harmonic: true,
            seed: 0x5eed,
        }
    }
}

/// Gram-Schmidt on the columns of `v`, applying the same combinations to `w`.
fn orthonormalise_pair(
    v: &mut Vec<Vec<f64>>,
    w: &mut Vec<Vec<f64>>,
) {
    let mut keep_v: Vec<Vec<f64>> = Vec::new();
    let mut keep_w: Vec<Vec<f64>> = Vec::new();
    for (mut x, mut y) in v.drain(..).zip(w.drain(..)) {
        for _ in 0..2 {
            for (p, pw) in keep_v.iter().zip(&keep_w) {
                let c = rdot(p, &x);
                for (a, b) in x.iter_mut().zip(p) {
                    *a -= c * b;
                }
                for (a, b) in y.iter_mut().zip(pw) {
                    *a -= c * b;
                }
            }
        }
        let nx = rnrm(&x);
        if nx > 1e-10 {
            for e in &mut x {
                *e /= nx;
            }
            for e in &mut y {
                *e /= nx;
            }
            keep_v.push(x);
            keep_w.push(y);
        }
    }
    *v = keep_v;
    *w = keep_w;
}

/// Jacobi-Davidson for the `k` eigenvalues of a symmetric operator closest
/// to `target` (interior eigenvalues included).
///
/// The search space is expanded with the solution `t` of the projected
/// correction equation `(I - Q Q^T)(A - theta I)(I - Q Q^T) t = -r` (`Q`
/// holds the locked eigenvectors and the current Ritz vector), solved
/// approximately by GMRES; `precond`, if given, writes `K^{-1} r` for an
/// approximation `K` of `A - theta I` and is made consistent with the
/// projector. Harmonic Ritz extraction picks the Ritz vector with the
/// smallest `|A - target|` residual direction, which is reliable for
/// interior targets where ordinary Rayleigh-Ritz values do not converge to
/// the nearest eigenvalues. Converged vectors are locked and deflated.
///
/// # Errors
/// [`EigsError::Invalid`] for `k == 0` or `k > n`; [`EigsError::Dense`] if a
/// projected eigenproblem fails.
pub fn eigs_jd_sym<A>(
    apply: A,
    n: usize,
    k: usize,
    target: f64,
    opts: &JdOptions,
    precond: Option<&OpFn<'_>>,
) -> Result<EigsResult, EigsError>
where
    A: Fn(&[f64], &mut [f64]),
{
    if k == 0 || k > n {
        return Err(EigsError::Invalid);
    }
    let mut rng = Rng::new(opts.seed);
    let max_dim = opts.max_dim.max(opts.min_dim + 2).min(n);
    let min_dim = opts.min_dim.min(max_dim - 1).max(1);
    let shifted = |x: &[f64], y: &mut [f64]| {
        apply(x, y);
        for (a, b) in y.iter_mut().zip(x) {
            *a -= target * b;
        }
    };
    let mut q: Vec<Vec<f64>> = Vec::new(); // locked eigenvectors
    let mut evals: Vec<f64> = Vec::new();
    let mut resid: Vec<f64> = Vec::new();
    let mut v: Vec<Vec<f64>> = Vec::new();
    let mut w: Vec<Vec<f64>> = Vec::new(); // (A - target) v
    let mut matvecs = 0usize;
    let project = |x: &mut [f64], qs: &[Vec<f64>]| {
        for _ in 0..2 {
            for p in qs {
                let c = rdot(p, x);
                for (a, b) in x.iter_mut().zip(p) {
                    *a -= c * b;
                }
            }
        }
    };
    let mut t0 = random_vec(n, &mut rng);
    project(&mut t0, &q);
    let tn = rnrm(&t0);
    for e in &mut t0 {
        *e /= tn;
    }
    let mut wt = vec![0.0; n];
    shifted(&t0, &mut wt);
    matvecs += 1;
    v.push(t0);
    w.push(wt);
    let mut iter = 0usize;
    let mut verify_iters = 0usize;
    while iter < opts.max_iter {
        iter += 1;
        if q.len() >= k {
            verify_iters += 1;
            if verify_iters > 150 {
                break;
            }
        }
        let m = v.len();
        // small eigenproblem
        let mut h2 = Mat::zeros(m, m); // W^T V = V^T (A - target) V
        for i in 0..m {
            for j in 0..m {
                h2.set(i, j, rdot(&v[i], &w[j]));
            }
        }
        let mut cols: Vec<(f64, Vec<f64>)> = Vec::new(); // (selection key, coefficient vector)
        if opts.harmonic {
            let mut h1 = Mat::zeros(m, m); // W^T W
            for i in 0..m {
                for j in 0..m {
                    h1.set(i, j, rdot(&w[i], &w[j]));
                }
            }
            if let Ok(l) = dense::cholesky(&h1) {
                // M = L^{-1} H2 L^{-T}
                let mut mm = Mat::zeros(m, m);
                for j in 0..m {
                    // solve L y = H2[:, j]
                    let mut y: Vec<f64> = (0..m)
                        .map(|i| f64::midpoint(h2.at(i, j), h2.at(j, i)))
                        .collect();
                    for i in 0..m {
                        for c in 0..i {
                            y[i] -= l.at(i, c) * y[c];
                        }
                        y[i] /= l.at(i, i);
                    }
                    for i in 0..m {
                        mm.set(i, j, y[i]);
                    }
                }
                let mut mt = Mat::zeros(m, m);
                for i in 0..m {
                    // solve L z = mm[i, :]^T
                    let mut y: Vec<f64> = (0..m).map(|j| mm.at(i, j)).collect();
                    for r in 0..m {
                        for c in 0..r {
                            y[r] -= l.at(r, c) * y[c];
                        }
                        y[r] /= l.at(r, r);
                    }
                    for j in 0..m {
                        mt.set(i, j, y[j]);
                    }
                }
                for i in 0..m {
                    for j in 0..i {
                        let s = f64::midpoint(mt.at(i, j), mt.at(j, i));
                        mt.set(i, j, s);
                        mt.set(j, i, s);
                    }
                }
                let (nu, z) = dense::eigen_symmetric(&mt).map_err(|_| EigsError::Dense)?;
                for (j, nuj) in nu.iter().enumerate() {
                    // c = L^{-T} z_j
                    let mut c: Vec<f64> = (0..m).map(|i| z.at(i, j)).collect();
                    for r in (0..m).rev() {
                        for s in r + 1..m {
                            c[r] -= l.at(s, r) * c[s];
                        }
                        c[r] /= l.at(r, r);
                    }
                    cols.push((-nuj.abs(), c));
                }
            }
        }
        if cols.is_empty() {
            let mut hs = h2.clone();
            for i in 0..m {
                for j in 0..i {
                    let s = f64::midpoint(hs.at(i, j), hs.at(j, i));
                    hs.set(i, j, s);
                    hs.set(j, i, s);
                }
            }
            let (th, z) = dense::eigen_symmetric(&hs).map_err(|_| EigsError::Dense)?;
            for (j, tj) in th.iter().enumerate() {
                cols.push((tj.abs(), (0..m).map(|i| z.at(i, j)).collect()));
            }
        }
        cols.sort_by(|a, b| a.0.total_cmp(&b.0));
        // best candidate
        let c0 = &cols[0].1;
        let mut u = vec![0.0; n];
        let mut wu = vec![0.0; n];
        for i in 0..m {
            for p in 0..n {
                u[p] += c0[i] * v[i][p];
                wu[p] += c0[i] * w[i][p];
            }
        }
        let un = rnrm(&u);
        for e in &mut u {
            *e /= un;
        }
        for e in &mut wu {
            *e /= un;
        }
        // A u = wu + target u; Rayleigh quotient
        let theta = target + rdot(&u, &wu);
        let r: Vec<f64> = wu
            .iter()
            .zip(&u)
            .map(|(a, b)| a - (theta - target) * b)
            .collect();
        let rn = rnrm(&r);
        if q.len() >= k {
            // verification phase: stop once the best remaining candidate is
            // not closer to the target than the k-th locked eigenvalue
            let mut ds: Vec<f64> = evals.iter().map(|e| (e - target).abs()).collect();
            ds.sort_by(f64::total_cmp);
            if (theta - target).abs() >= ds[k - 1] {
                break;
            }
        }
        if rn <= opts.tol * theta.abs().max(1.0) {
            q.push(u);
            evals.push(theta);
            resid.push(rn);
            // deflate: keep the other Ritz vectors
            let take = min_dim.min(cols.len() - 1);
            let mut nv: Vec<Vec<f64>> = Vec::new();
            let mut nw: Vec<Vec<f64>> = Vec::new();
            for (_, c) in cols.iter().skip(1).take(take) {
                let (mut x, mut y) = (vec![0.0; n], vec![0.0; n]);
                for i in 0..m {
                    for p in 0..n {
                        x[p] += c[i] * v[i][p];
                        y[p] += c[i] * w[i][p];
                    }
                }
                nv.push(x);
                nw.push(y);
            }
            // remove the locked direction from the remaining vectors
            for (x, y) in nv.iter_mut().zip(nw.iter_mut()) {
                let last = &q[q.len() - 1];
                let cdot = rdot(last, x);
                for (a, b) in x.iter_mut().zip(last) {
                    *a -= cdot * b;
                }
                let mut ay = vec![0.0; n];
                shifted(last, &mut ay);
                for (a, b) in y.iter_mut().zip(&ay) {
                    *a -= cdot * b;
                }
            }
            orthonormalise_pair(&mut nv, &mut nw);
            if nv.is_empty() {
                let mut t = random_vec(n, &mut rng);
                project(&mut t, &q);
                let tn = rnrm(&t);
                for e in &mut t {
                    *e /= tn;
                }
                let mut wt = vec![0.0; n];
                shifted(&t, &mut wt);
                matvecs += 1;
                nv.push(t);
                nw.push(wt);
            }
            v = nv;
            w = nw;
            continue;
        }
        // correction equation with projector P = I - Qt Qt^T, Qt = [Q, u]
        let mut qt: Vec<Vec<f64>> = q.clone();
        qt.push(u.clone());
        let theta_c = theta;
        let proj_apply = |x: &[f64], y: &mut [f64]| {
            let mut px = x.to_vec();
            project(&mut px, &qt);
            let mut ay = vec![0.0; n];
            apply(&px, &mut ay);
            for (a, b) in ay.iter_mut().zip(&px) {
                *a -= theta_c * b;
            }
            project(&mut ay, &qt);
            y.copy_from_slice(&ay);
        };
        let (kinv_q, small): (Vec<Vec<f64>>, Option<Mat>) = precond.map_or_else(
            || (Vec::new(), None),
            |pc| {
                let kq: Vec<Vec<f64>> = qt
                    .iter()
                    .map(|x| {
                        let mut y = vec![0.0; n];
                        pc(x, &mut y);
                        y
                    })
                    .collect();
                let nq = qt.len();
                let mut sm = Mat::zeros(nq, nq);
                for i in 0..nq {
                    for j in 0..nq {
                        sm.set(i, j, rdot(&qt[i], &kq[j]));
                    }
                }
                (kq, Some(sm))
            },
        );
        let precond_apply = |x: &[f64], y: &mut [f64]| {
            if let (Some(pc), Some(sm)) = (precond, small.as_ref()) {
                let mut z = vec![0.0; n];
                pc(x, &mut z);
                let rhs: Vec<f64> = qt.iter().map(|p| rdot(p, &z)).collect();
                if let Ok(g) = dense::solve(sm, &rhs) {
                    for (gi, kq) in g.iter().zip(&kinv_q) {
                        for (a, b) in z.iter_mut().zip(kq) {
                            *a -= gi * b;
                        }
                    }
                }
                y.copy_from_slice(&z);
            } else {
                y.copy_from_slice(x);
                project(y, &qt);
            }
        };
        let mut rhs: Vec<f64> = r.iter().map(|e| -e).collect();
        project(&mut rhs, &qt);
        let sol = gmres(
            proj_apply,
            precond_apply,
            &rhs,
            None,
            opts.inner_iters,
            opts.inner_tol,
            opts.inner_iters,
        );
        matvecs += sol.iterations;
        let mut t = sol.x;
        project(&mut t, &qt);
        // orthogonalise against the search space
        for _ in 0..2 {
            for p in &v {
                let c = rdot(p, &t);
                for (a, b) in t.iter_mut().zip(p) {
                    *a -= c * b;
                }
            }
        }
        let mut tn = rnrm(&t);
        if tn.is_nan() || tn <= 1e-12 || !tn.is_finite() {
            t = random_vec(n, &mut rng);
            project(&mut t, &q);
            for p in &v {
                let c = rdot(p, &t);
                for (a, b) in t.iter_mut().zip(p) {
                    *a -= c * b;
                }
            }
            tn = rnrm(&t);
        }
        for e in &mut t {
            *e /= tn;
        }
        let mut wt = vec![0.0; n];
        shifted(&t, &mut wt);
        matvecs += 1;
        v.push(t);
        w.push(wt);
        // implicit restart
        if v.len() > max_dim {
            let mut nv: Vec<Vec<f64>> = Vec::new();
            let mut nw: Vec<Vec<f64>> = Vec::new();
            let mcur = m;
            for (_, c) in cols.iter().take(min_dim) {
                let (mut x, mut y) = (vec![0.0; n], vec![0.0; n]);
                for i in 0..mcur {
                    for p in 0..n {
                        x[p] += c[i] * v[i][p];
                        y[p] += c[i] * w[i][p];
                    }
                }
                nv.push(x);
                nw.push(y);
            }
            // keep the newest direction too
            nv.push(v[v.len() - 1].clone());
            nw.push(w[w.len() - 1].clone());
            orthonormalise_pair(&mut nv, &mut nw);
            v = nv;
            w = nw;
        }
    }
    // order by distance to the target
    let mut order: Vec<usize> = (0..q.len()).collect();
    order.sort_by(|&a, &b| {
        (evals[a] - target)
            .abs()
            .total_cmp(&(evals[b] - target).abs())
    });
    let converged = q.len() >= k;
    order.truncate(k);
    Ok(EigsResult {
        values: order.iter().map(|&i| (evals[i], 0.0)).collect(),
        vectors_re: order.iter().map(|&i| q[i].clone()).collect(),
        vectors_im: order.iter().map(|_| vec![0.0; n]).collect(),
        residuals: order.iter().map(|&i| resid[i]).collect(),
        converged,
        restarts: iter,
        matvecs,
    })
}

// ------------------------------------------- non-symmetric Jacobi-Davidson

fn proj_c(
    x: &mut [C],
    qs: &[Vec<C>],
) {
    for _ in 0..2 {
        for p in qs {
            let c = cdot(p, x);
            for (a, b) in x.iter_mut().zip(p) {
                *a -= c * b;
            }
        }
    }
}

/// Complex Cholesky factor `L` (lower) of a Hermitian positive definite matrix.
fn cholesky_c(h: &[Vec<C>]) -> Option<Vec<Vec<C>>> {
    let m = h.len();
    let mut l = vec![vec![C::new(0.0, 0.0); m]; m];
    for i in 0..m {
        for j in 0..=i {
            let mut s = h[i][j];
            for k in 0..j {
                s -= l[i][k] * l[j][k].conj();
            }
            if i == j {
                if s.re <= 0.0 || !s.re.is_finite() {
                    return None;
                }
                l[i][i] = C::new(s.re.sqrt(), 0.0);
            } else {
                l[i][j] = s / l[j][j];
            }
        }
    }
    Some(l)
}

/// Solves `L y = b` for lower triangular `L`.
fn forward_c(
    l: &[Vec<C>],
    b: &[C],
) -> Vec<C> {
    let m = b.len();
    let mut y = b.to_vec();
    for i in 0..m {
        for k in 0..i {
            let t = l[i][k] * y[k];
            y[i] -= t;
        }
        y[i] /= l[i][i];
    }
    y
}

/// Solves `L^H x = b` for lower triangular `L`.
fn backward_h_c(
    l: &[Vec<C>],
    b: &[C],
) -> Vec<C> {
    let m = b.len();
    let mut x = b.to_vec();
    for i in (0..m).rev() {
        for k in i + 1..m {
            let t = l[k][i].conj() * x[k];
            x[i] -= t;
        }
        x[i] /= l[i][i].conj();
    }
    x
}

/// Unit eigenvector `y` of the leading `(i + 1)` block of an upper triangular
/// `t` for the eigenvalue `t[i][i]`.
fn triangular_eigvec(
    t: &[Vec<C>],
    i: usize,
) -> Vec<C> {
    let mut y = vec![C::new(0.0, 0.0); i + 1];
    y[i] = C::new(1.0, 0.0);
    let lam = t[i][i];
    let floor = 1e-14 * lam.norm().max(1.0);
    for r in (0..i).rev() {
        let mut s = C::new(0.0, 0.0);
        for c in r + 1..=i {
            s += t[r][c] * y[c];
        }
        let mut d = t[r][r] - lam;
        if d.norm() < floor {
            d = C::new(floor, 0.0);
        }
        y[r] = -s / d;
    }
    let yn = cnrm(&y);
    y.iter().map(|e| e / yn).collect()
}

/// Right-preconditioned complex GMRES without restart: approximately solves
/// `A x = b` in at most `iters` iterations (relative tolerance `tol`).
fn gmres_c(
    apply: &dyn Fn(&[C]) -> Vec<C>,
    precond: &dyn Fn(&[C]) -> Vec<C>,
    b: &[C],
    iters: usize,
    tol: f64,
) -> Vec<C> {
    let n = b.len();
    let zero = C::new(0.0, 0.0);
    let beta = cnrm(b);
    if beta == 0.0 {
        return vec![zero; n];
    }
    let mut v: Vec<Vec<C>> = vec![b.iter().map(|e| e / beta).collect()];
    let mut z: Vec<Vec<C>> = Vec::new();
    let mut h = vec![vec![zero; iters]; iters + 1];
    let mut cs = vec![0.0; iters];
    let mut sn = vec![zero; iters];
    let mut g = vec![zero; iters + 1];
    g[0] = C::new(beta, 0.0);
    let mut used = 0;
    for j in 0..iters {
        let zj = precond(&v[j]);
        let mut w = apply(&zj);
        z.push(zj);
        for i in 0..=j {
            let c = cdot(&v[i], &w);
            h[i][j] = c;
            for (a, bb) in w.iter_mut().zip(&v[i]) {
                *a -= c * bb;
            }
        }
        let hn = cnrm(&w);
        h[j + 1][j] = C::new(hn, 0.0);
        for i in 0..j {
            let (x, y) = (h[i][j], h[i + 1][j]);
            h[i][j] = x * cs[i] + sn[i] * y;
            h[i + 1][j] = -sn[i].conj() * x + y * cs[i];
        }
        let (c, s) = givens(h[j][j], h[j + 1][j]);
        cs[j] = c;
        sn[j] = s;
        let (x, y) = (h[j][j], h[j + 1][j]);
        h[j][j] = x * c + s * y;
        h[j + 1][j] = C::new(0.0, 0.0);
        let gj = g[j];
        g[j] = gj * c;
        g[j + 1] = -s.conj() * gj;
        used = j + 1;
        if g[j + 1].norm() <= tol * beta || hn < 1e-300 {
            break;
        }
        v.push(w.iter().map(|e| e / hn).collect());
    }
    let mut y = vec![zero; used];
    for i in (0..used).rev() {
        let mut s = g[i];
        for kk in i + 1..used {
            s -= h[i][kk] * y[kk];
        }
        y[i] = s / h[i][i];
    }
    let mut x = vec![zero; n];
    for (yi, zi) in y.iter().zip(&z) {
        for (a, bb) in x.iter_mut().zip(zi) {
            *a += yi * bb;
        }
    }
    x
}

/// Gram-Schmidt on the columns of `v` applying the same combinations to `w`.
fn orthonormalise_pair_c(
    v: &mut Vec<Vec<C>>,
    w: &mut Vec<Vec<C>>,
) {
    let mut keep_v: Vec<Vec<C>> = Vec::new();
    let mut keep_w: Vec<Vec<C>> = Vec::new();
    for (mut x, mut y) in v.drain(..).zip(w.drain(..)) {
        for _ in 0..2 {
            for (p, pw) in keep_v.iter().zip(&keep_w) {
                let c = cdot(p, &x);
                for (a, b) in x.iter_mut().zip(p) {
                    *a -= c * b;
                }
                for (a, b) in y.iter_mut().zip(pw) {
                    *a -= c * b;
                }
            }
        }
        let nx = cnrm(&x);
        if nx > 1e-10 {
            for e in &mut x {
                *e /= nx;
            }
            for e in &mut y {
                *e /= nx;
            }
            keep_v.push(x);
            keep_w.push(y);
        }
    }
    *v = keep_v;
    *w = keep_w;
}

/// Gaussian elimination with partial pivoting for a small complex system.
fn csolve_small(
    a: &[Vec<C>],
    b: &[C],
) -> Option<Vec<C>> {
    let n = b.len();
    let mut a: Vec<Vec<C>> = a.to_vec();
    let mut b = b.to_vec();
    for kk in 0..n {
        let p = (kk..n).max_by(|&i, &j| a[i][kk].norm().total_cmp(&a[j][kk].norm()))?;
        if a[p][kk].norm() < 1e-300 {
            return None;
        }
        a.swap(kk, p);
        b.swap(kk, p);
        for i in kk + 1..n {
            let f = a[i][kk] / a[kk][kk];
            for c in kk..n {
                let t = a[kk][c];
                a[i][c] -= f * t;
            }
            let t = b[kk];
            b[i] -= f * t;
        }
    }
    let mut x = vec![C::new(0.0, 0.0); n];
    for i in (0..n).rev() {
        let mut s = b[i];
        for c in i + 1..n {
            s -= a[i][c] * x[c];
        }
        x[i] = s / a[i][i];
    }
    Some(x)
}

/// Harmonic (or standard) Ritz candidates of the search space: coefficient
/// vectors `c` (so that `u = V c`) of the `count` best pairs, best first.
fn jd_candidates(
    v: &[Vec<C>],
    w: &[Vec<C>],
    harmonic: bool,
    count: usize,
) -> Result<Vec<Vec<C>>, EigsError> {
    let zero = C::new(0.0, 0.0);
    let m = v.len();
    let ncand = count.min(m);
    let lift = |t: &[Vec<C>], u: &[Vec<C>], i: usize| -> Vec<C> {
        let y = triangular_eigvec(t, i);
        let mut z = vec![zero; m];
        for (j, yj) in y.iter().enumerate() {
            for r in 0..m {
                z[r] += u[r][j] * yj;
            }
        }
        z
    };
    if harmonic {
        let mut h2 = vec![vec![zero; m]; m]; // W^H V
        let mut h1 = vec![vec![zero; m]; m]; // W^H W
        for i in 0..m {
            for j in 0..m {
                h2[i][j] = cdot(&w[i], &v[j]);
                h1[i][j] = cdot(&w[i], &w[j]);
            }
        }
        if let Some(l) = cholesky_c(&h1) {
            // M = L^{-1} H2 L^{-H}
            let x1: Vec<Vec<C>> = (0..m)
                .map(|j| {
                    let col: Vec<C> = (0..m).map(|i| h2[i][j]).collect();
                    forward_c(&l, &col)
                })
                .collect(); // x1[j] = column j of X = L^{-1} H2
            let mut mm = vec![vec![zero; m]; m];
            for i in 0..m {
                let rowx: Vec<C> = (0..m).map(|j| x1[j][i].conj()).collect();
                let y = forward_c(&l, &rowx);
                for j in 0..m {
                    mm[i][j] = y[j].conj();
                }
            }
            if let Some((mut t, mut u)) = complex_schur(&mm) {
                reorder_schur(&mut t, &mut u, ncand, &|z: C| -z.norm());
                return Ok((0..ncand)
                    .map(|i| backward_h_c(&l, &lift(&t, &u, i)))
                    .collect());
            }
        }
    }
    // standard Rayleigh-Ritz: V^H (A - tau) V = V^H W
    let hm: Vec<Vec<C>> = (0..m)
        .map(|i| (0..m).map(|j| cdot(&v[i], &w[j])).collect())
        .collect();
    let (mut t, mut u) = complex_schur(&hm).ok_or(EigsError::Dense)?;
    reorder_schur(&mut t, &mut u, ncand, &|z: C| z.norm());
    Ok((0..ncand).map(|i| lift(&t, &u, i)).collect())
}

/// Jacobi-Davidson for the `k` eigenvalues of a general **non-symmetric**
/// real operator closest to the (complex) `target`, interior eigenvalues
/// included.
///
/// The method works in complex arithmetic with a *partial Schur* deflation:
/// the search space stays orthogonal to the locked Schur vectors `Q`
/// (`A Q = Q R`, `R` upper triangular), harmonic Ritz pairs of the pencil
/// `(W^H V, W^H W)`, `W = (I - Q Q^H)(A - target) V`, are extracted with the
/// complex Schur decomposition [`complex_schur`], and the correction equation
/// `(I - Q~ Q~^H)(A - theta)(I - Q~ Q~^H) t = -r` (`Q~ = [Q, u]`) is solved by
/// complex GMRES, optionally right-preconditioned by `precond`
/// (`K^{-1}` for a real `K ~ A - theta`, applied to the real and imaginary
/// parts). Eigenvalues are the diagonal of `R` and the eigenvectors are
/// `Q y` with `y` the eigenvectors of `R`; `residuals` are true residuals
/// `|A x - lambda x|`. As in [`eigs_jd_sym`], iteration continues after `k`
/// locks until no unlocked candidate is closer to the target.
///
/// # Errors
/// [`EigsError::Invalid`] for `k == 0` or `k > n`; [`EigsError::Dense`] if a
/// projected problem fails.
pub fn eigs_jd<A>(
    apply: A,
    n: usize,
    k: usize,
    target: (f64, f64),
    opts: &JdOptions,
    precond: Option<&OpFn<'_>>,
) -> Result<EigsResult, EigsError>
where
    A: Fn(&[f64], &mut [f64]),
{
    if k == 0 || k > n {
        return Err(EigsError::Invalid);
    }
    let zero = C::new(0.0, 0.0);
    let tau = C::new(target.0, target.1);
    let max_dim = opts.max_dim.max(opts.min_dim + 2).min(n);
    let min_dim = opts.min_dim.min(max_dim - 1).max(1);
    let mut rng = Rng::new(opts.seed);
    let mut matvecs = 0usize;
    let a_plain = |x: &[C], mv: &mut usize| -> Vec<C> {
        let mut out = vec![zero; n];
        let (mut re, mut im) = (vec![0.0; n], vec![0.0; n]);
        apply_complex(&apply, x, &mut out, &mut re, &mut im);
        *mv += 1;
        out
    };
    // deflated shifted operator (I - Q Q^H)(A - tau) x
    let op = |x: &[C], q: &[Vec<C>], mv: &mut usize| -> Vec<C> {
        let mut y = a_plain(x, mv);
        for (a, b) in y.iter_mut().zip(x) {
            *a -= tau * b;
        }
        proj_c(&mut y, q);
        y
    };
    let rand_vec = |rng: &mut Rng, q: &[Vec<C>], v: &[Vec<C>]| -> Vec<C> {
        let mut t: Vec<C> = (0..n).map(|_| C::new(rng.normal(), 0.0)).collect();
        proj_c(&mut t, q);
        proj_c(&mut t, v);
        let tn = cnrm(&t);
        t.iter().map(|e| e / tn).collect()
    };
    let mut q: Vec<Vec<C>> = Vec::new();
    let mut rmat: Vec<Vec<C>> = Vec::new(); // R[i][j]
    let mut v: Vec<Vec<C>> = Vec::new();
    let mut w: Vec<Vec<C>> = Vec::new();
    let t0 = rand_vec(&mut rng, &q, &v);
    w.push(op(&t0, &q, &mut matvecs));
    v.push(t0);
    let mut iter = 0usize;
    let mut verify = 0usize;
    while iter < opts.max_iter {
        iter += 1;
        if q.len() >= k {
            verify += 1;
            if verify > 150 {
                break;
            }
        }
        let m = v.len();
        let cands = jd_candidates(&v, &w, opts.harmonic, min_dim + 1)?;
        let c0 = &cands[0];
        let mut u = vec![zero; n];
        let mut wu = vec![zero; n];
        for i in 0..m {
            for p in 0..n {
                u[p] += c0[i] * v[i][p];
                wu[p] += c0[i] * w[i][p];
            }
        }
        let un = cnrm(&u);
        for e in &mut u {
            *e /= un;
        }
        for e in &mut wu {
            *e /= un;
        }
        let theta = tau + cdot(&u, &wu);
        let r: Vec<C> = wu
            .iter()
            .zip(&u)
            .map(|(a, b)| a - (theta - tau) * b)
            .collect();
        let rn = cnrm(&r);
        if q.len() >= k {
            let mut ds: Vec<f64> = (0..q.len()).map(|i| (rmat[i][i] - tau).norm()).collect();
            ds.sort_by(f64::total_cmp);
            if (theta - tau).norm() >= ds[k - 1] {
                break;
            }
        }
        if rn <= opts.tol * theta.norm().max(1.0) {
            // lock the Schur vector: R gets the column Q^H A u
            let au = a_plain(&u, &mut matvecs);
            let kq = q.len();
            for row in &mut rmat {
                row.push(zero);
            }
            rmat.push(vec![zero; kq + 1]);
            for (i, qi) in q.iter().enumerate() {
                rmat[i][kq] = cdot(qi, &au);
            }
            rmat[kq][kq] = theta;
            q.push(u);
            // restart with the remaining candidates
            let mut nv: Vec<Vec<C>> = Vec::new();
            for c in cands.iter().skip(1).take(min_dim) {
                let mut x = vec![zero; n];
                for i in 0..m {
                    for p in 0..n {
                        x[p] += c[i] * v[i][p];
                    }
                }
                proj_c(&mut x, &q);
                nv.push(x);
            }
            let mut dummy: Vec<Vec<C>> = nv.clone();
            orthonormalise_pair_c(&mut nv, &mut dummy);
            if nv.is_empty() {
                nv.push(rand_vec(&mut rng, &q, &[]));
            }
            let nw: Vec<Vec<C>> = nv.iter().map(|x| op(x, &q, &mut matvecs)).collect();
            v = nv;
            w = nw;
            continue;
        }
        // correction equation with Q~ = [Q, u]
        let mut qt = q.clone();
        qt.push(u.clone());
        let proj_op = |x: &[C]| -> Vec<C> {
            let mut px = x.to_vec();
            proj_c(&mut px, &qt);
            let mut mv = 0usize;
            let mut y = a_plain(&px, &mut mv);
            for (a, b) in y.iter_mut().zip(&px) {
                *a -= theta * b;
            }
            proj_c(&mut y, &qt);
            y
        };
        let real_pc = |pc: &OpFn<'_>, x: &[C]| -> Vec<C> {
            let (xr, xi): (Vec<f64>, Vec<f64>) = x.iter().map(|e| (e.re, e.im)).unzip();
            let (mut yr, mut yi) = (vec![0.0; n], vec![0.0; n]);
            pc(&xr, &mut yr);
            pc(&xi, &mut yi);
            (0..n).map(|i| C::new(yr[i], yi[i])).collect()
        };
        let (kinv_q, small): (Vec<Vec<C>>, Option<Vec<Vec<C>>>) = precond.map_or_else(
            || (Vec::new(), None),
            |pc| {
                let kq: Vec<Vec<C>> = qt.iter().map(|x| real_pc(pc, x)).collect();
                let nq = qt.len();
                let sm: Vec<Vec<C>> = (0..nq)
                    .map(|i| (0..nq).map(|j| cdot(&qt[i], &kq[j])).collect())
                    .collect();
                (kq, Some(sm))
            },
        );
        let precond_apply = |x: &[C]| -> Vec<C> {
            if let (Some(pc), Some(sm)) = (precond, small.as_ref()) {
                let mut z = real_pc(pc, x);
                let rhs: Vec<C> = qt.iter().map(|p| cdot(p, &z)).collect();
                if let Some(g) = csolve_small(sm, &rhs) {
                    for (gi, kqv) in g.iter().zip(&kinv_q) {
                        for (a, b) in z.iter_mut().zip(kqv) {
                            *a -= gi * b;
                        }
                    }
                }
                z
            } else {
                let mut y = x.to_vec();
                proj_c(&mut y, &qt);
                y
            }
        };
        let mut rhs: Vec<C> = r.iter().map(|e| -e).collect();
        proj_c(&mut rhs, &qt);
        let mut t = gmres_c(
            &proj_op,
            &precond_apply,
            &rhs,
            opts.inner_iters,
            opts.inner_tol,
        );
        matvecs += opts.inner_iters;
        proj_c(&mut t, &qt);
        proj_c(&mut t, &v);
        let mut tn = cnrm(&t);
        if tn.is_nan() || tn <= 1e-12 {
            t = rand_vec(&mut rng, &q, &v);
            tn = 1.0;
        }
        for e in &mut t {
            *e /= tn;
        }
        w.push(op(&t, &q, &mut matvecs));
        v.push(t);
        if v.len() > max_dim {
            let mut nv: Vec<Vec<C>> = Vec::new();
            let mut nw: Vec<Vec<C>> = Vec::new();
            for c in cands.iter().take(min_dim) {
                let (mut x, mut y) = (vec![zero; n], vec![zero; n]);
                for i in 0..m {
                    for p in 0..n {
                        x[p] += c[i] * v[i][p];
                        y[p] += c[i] * w[i][p];
                    }
                }
                nv.push(x);
                nw.push(y);
            }
            nv.push(v[v.len() - 1].clone());
            nw.push(w[w.len() - 1].clone());
            orthonormalise_pair_c(&mut nv, &mut nw);
            v = nv;
            w = nw;
        }
    }
    // eigenpairs from the partial Schur form A Q = Q R
    let kq = q.len();
    let mut order: Vec<usize> = (0..kq).collect();
    order.sort_by(|&a, &b| {
        (rmat[a][a] - tau)
            .norm()
            .total_cmp(&(rmat[b][b] - tau).norm())
    });
    order.truncate(k);
    let mut out = EigsResult {
        values: Vec::new(),
        vectors_re: Vec::new(),
        vectors_im: Vec::new(),
        residuals: Vec::new(),
        converged: kq >= k,
        restarts: iter,
        matvecs,
    };
    for &i in &order {
        let y = triangular_eigvec(&rmat, i);
        let mut x = vec![zero; n];
        for (j, yj) in y.iter().enumerate() {
            for p in 0..n {
                x[p] += q[j][p] * yj;
            }
        }
        let xn = cnrm(&x);
        let big = x
            .iter()
            .copied()
            .max_by(|p, r| p.norm().total_cmp(&r.norm()))
            .unwrap_or(zero);
        let phase = if big.norm() == 0.0 {
            C::new(1.0, 0.0)
        } else {
            big.conj() / big.norm()
        };
        let xs: Vec<C> = x.iter().map(|e| e * phase / xn).collect();
        let lam = rmat[i][i];
        let mut mv = 0usize;
        let ax = a_plain(&xs, &mut mv);
        let res: f64 = ax
            .iter()
            .zip(&xs)
            .map(|(a, b)| (a - lam * b).norm_sqr())
            .sum::<f64>()
            .sqrt();
        out.values.push((lam.re, lam.im));
        out.vectors_re.push(xs.iter().map(|e| e.re).collect());
        out.vectors_im.push(xs.iter().map(|e| e.im).collect());
        out.residuals.push(res);
    }
    Ok(out)
}

// ------------------------------------------------- two-sided Jacobi-Davidson

/// Result of [`eigs_jd_two_sided`]: eigenvalues with right and left
/// eigenvectors and the condition numbers read off their angle.
#[derive(Debug, Clone, PartialEq)]
pub struct TwoSidedResult {
    /// Eigenvalues `(re, im)`, nearest to the target first.
    pub values: Vec<(f64, f64)>,
    /// Real parts of the right eigenvectors `A u = lambda u` (unit 2-norm).
    pub right_re: Vec<Vec<f64>>,
    /// Imaginary parts of the right eigenvectors.
    pub right_im: Vec<Vec<f64>>,
    /// Real parts of the left eigenvectors `w^H A = lambda w^H` (unit 2-norm,
    /// phase chosen so that `w^H u > 0`).
    pub left_re: Vec<Vec<f64>>,
    /// Imaginary parts of the left eigenvectors.
    pub left_im: Vec<Vec<f64>>,
    /// True residuals `|A u - lambda u|`.
    pub right_residuals: Vec<f64>,
    /// True residuals `|A^T w - conj(lambda) w|` (for a real operator).
    pub left_residuals: Vec<f64>,
    /// Eigenvalue condition numbers `1 / |w^H u|` of the unit vectors: the
    /// first-order sensitivity of the eigenvalue to a perturbation of `A`.
    pub condition: Vec<f64>,
    /// Whether `k` pairs were locked.
    pub converged: bool,
    /// Outer iterations.
    pub iterations: usize,
    /// Operator applications (of `A` and `A^T` together).
    pub matvecs: usize,
}

/// Oblique projector `x - Q D^{-1} (P^H x)` with `D = P^H Q`.
fn oblique(
    x: &[C],
    q: &[Vec<C>],
    p: &[Vec<C>],
    d: &[Vec<C>],
) -> Vec<C> {
    let mut out = x.to_vec();
    if q.is_empty() {
        return out;
    }
    let rhs: Vec<C> = p.iter().map(|pi| cdot(pi, x)).collect();
    if let Some(c) = csolve_small(d, &rhs) {
        for (ci, qi) in c.iter().zip(q) {
            for (a, b) in out.iter_mut().zip(qi) {
                *a -= ci * b;
            }
        }
    }
    out
}

/// The Gram matrix `P^H Q`.
fn gram_pq(
    p: &[Vec<C>],
    q: &[Vec<C>],
) -> Vec<Vec<C>> {
    p.iter()
        .map(|pi| q.iter().map(|qj| cdot(pi, qj)).collect())
        .collect()
}

/// The conjugate transpose of a small matrix.
fn adjoint_small(a: &[Vec<C>]) -> Vec<Vec<C>> {
    let rows = a.len();
    let cols = a.first().map_or(0, Vec::len);
    (0..cols)
        .map(|j| (0..rows).map(|i| a[i][j].conj()).collect())
        .collect()
}

/// The Petrov-Galerkin data of the current search spaces.
struct PetrovData {
    /// Eigenvalue indices of the small matrix, nearest to zero first.
    order: Vec<usize>,
    /// Schur form and vectors of `G = S^{-1} H`.
    t: Vec<Vec<C>>,
    u: Vec<Vec<C>>,
    /// Matrix whose columns are the right eigenvectors of `G`.
    y_mat: Vec<Vec<C>>,
    /// `S = W^H V`.
    s_mat: Vec<Vec<C>>,
}

impl PetrovData {
    /// Right coefficient vector `y` of the eigenvalue `idx`.
    fn right(
        &self,
        idx: usize,
    ) -> Vec<C> {
        let yy = triangular_eigvec(&self.t, idx);
        let m = self.u.len();
        (0..m)
            .map(|r| (0..yy.len()).map(|j| self.u[r][j] * yy[j]).sum::<C>())
            .collect()
    }

    /// Left coefficient vector `z` of the eigenvalue `idx` (`z^H H = theta z^H S`).
    fn left(
        &self,
        idx: usize,
    ) -> Option<Vec<C>> {
        let m = self.u.len();
        let mut e = vec![C::new(0.0, 0.0); m];
        e[idx] = C::new(1.0, 0.0);
        let zeta = csolve_small(&adjoint_small(&self.y_mat), &e)?;
        csolve_small(&adjoint_small(&self.s_mat), &zeta)
    }
}

/// Two-sided Jacobi-Davidson for the `k` eigenvalues of a general real
/// operator nearest to the (complex) `target`, with **left and right**
/// eigenvectors.
///
/// Two search spaces `V` (right) and `W` (left), each orthonormal, are
/// expanded together. The projected pencil `(W^H (A - tau) V, W^H V)` gives
/// Petrov-Galerkin pairs: the Ritz value `theta`, the right Ritz vector
/// `u = V y` with residual `r = (A - tau - theta) u` orthogonal to `W`, and the
/// left Ritz vector `w = W z` with residual `s = (A - tau - theta)^H w`
/// orthogonal to `V`; the left vectors of the small pencil come from the
/// inverse of its right-eigenvector matrix, so the *two-sided Rayleigh
/// quotient* `theta = w^H A u / w^H u` is stationary (its error is quadratic in
/// the residuals, unlike the one-sided quotient for a non-normal `A`). The two
/// oblique correction equations
/// `(I - u w^H / w^H u)(A - theta)(I - u w^H / w^H u) t = -r` and the adjoint
/// one for `s` are solved by GMRES. Converged pairs are deflated *obliquely*:
/// locked right vectors `Q` and left vectors `P` (scaled so that
/// `p_j^H q_j = 1`, mutually biorthogonal for distinct eigenvalues) define the
/// projector `I - Q (P^H Q)^{-1} P^H` that keeps both spaces away from them.
///
/// `apply` is `x -> A x`, `apply_t` is `x -> A^T x` (the adjoint of the real
/// operator). A pair is accepted when both residuals are below
/// `opts.tol * max(|lambda|, 1)`. The standard (not harmonic) extraction is
/// used, so `target` should lie near the wanted eigenvalues.
///
/// # Errors
/// [`EigsError::Invalid`] for `k == 0` or `k > n`; [`EigsError::Dense`] if a
/// projected problem fails.
#[allow(clippy::too_many_lines, clippy::similar_names)]
pub fn eigs_jd_two_sided<A, AT>(
    apply: A,
    apply_t: AT,
    n: usize,
    k: usize,
    target: (f64, f64),
    opts: &JdOptions,
) -> Result<TwoSidedResult, EigsError>
where
    A: Fn(&[f64], &mut [f64]),
    AT: Fn(&[f64], &mut [f64]),
{
    if k == 0 || k > n {
        return Err(EigsError::Invalid);
    }
    let zero = C::new(0.0, 0.0);
    let tau = C::new(target.0, target.1);
    let max_dim = opts.max_dim.max(opts.min_dim + 2).min(n);
    let min_dim = opts.min_dim.min(max_dim - 1).max(1);
    let mut rng = Rng::new(opts.seed);
    let mut matvecs = 0_usize;
    let run_a = |x: &[C], mv: &mut usize| -> Vec<C> {
        let mut out = vec![zero; n];
        let (mut re, mut im) = (vec![0.0; n], vec![0.0; n]);
        apply_complex(&apply, x, &mut out, &mut re, &mut im);
        *mv += 1;
        out
    };
    let run_t = |x: &[C], mv: &mut usize| -> Vec<C> {
        let mut out = vec![zero; n];
        let (mut re, mut im) = (vec![0.0; n], vec![0.0; n]);
        apply_complex(&apply_t, x, &mut out, &mut re, &mut im);
        *mv += 1;
        out
    };
    // locked right (q) / left (p) vectors, scaled so that p_j^H q_j = 1
    let mut q: Vec<Vec<C>> = Vec::new();
    let mut p: Vec<Vec<C>> = Vec::new();
    let mut lam: Vec<C> = Vec::new();
    let mut unit_pairs: Vec<(Vec<C>, Vec<C>)> = Vec::new();
    let deflate_r =
        |x: &[C], q: &[Vec<C>], p: &[Vec<C>]| -> Vec<C> { oblique(x, q, p, &gram_pq(p, q)) };
    let deflate_l = |x: &[C], q: &[Vec<C>], p: &[Vec<C>]| -> Vec<C> {
        oblique(x, p, q, &adjoint_small(&gram_pq(p, q)))
    };
    // (A - tau) x with the oblique deflation applied to the result
    let op_r = |x: &[C], q: &[Vec<C>], p: &[Vec<C>], mv: &mut usize| -> Vec<C> {
        let mut y = run_a(x, mv);
        for (a, b) in y.iter_mut().zip(x) {
            *a -= tau * b;
        }
        deflate_r(&y, q, p)
    };
    // (A - tau)^H x with the adjoint deflation I - P D^{-H} Q^H
    let op_l = |x: &[C], q: &[Vec<C>], p: &[Vec<C>], mv: &mut usize| -> Vec<C> {
        let mut y = run_t(x, mv);
        for (a, b) in y.iter_mut().zip(x) {
            *a -= tau.conj() * b;
        }
        deflate_l(&y, q, p)
    };
    let random_vec =
        |rng: &mut Rng| -> Vec<C> { (0..n).map(|_| C::new(rng.normal(), 0.0)).collect() };
    let orthonormal_against = |x: &mut Vec<C>, basis: &[Vec<C>]| -> bool {
        for _ in 0..2 {
            for b in basis {
                let c = cdot(b, x);
                for (a, bb) in x.iter_mut().zip(b) {
                    *a -= c * bb;
                }
            }
        }
        let nx = cnrm(x);
        if nx <= 1e-10 || nx.is_nan() {
            return false;
        }
        for e in x.iter_mut() {
            *e /= nx;
        }
        true
    };
    let mut v: Vec<Vec<C>> = Vec::new();
    let mut w: Vec<Vec<C>> = Vec::new();
    let mut av: Vec<Vec<C>> = Vec::new();
    let mut aw: Vec<Vec<C>> = Vec::new();
    {
        let mut x = random_vec(&mut rng);
        let mut y = random_vec(&mut rng);
        let _ = orthonormal_against(&mut x, &[]);
        let _ = orthonormal_against(&mut y, &[]);
        av.push(op_r(&x, &q, &p, &mut matvecs));
        aw.push(op_l(&y, &q, &p, &mut matvecs));
        v.push(x);
        w.push(y);
    }
    // combination `sum_i c_i b_i` of basis vectors
    let combine = |c: &[C], basis: &[Vec<C>]| -> Vec<C> {
        let mut out = vec![zero; n];
        for (ci, b) in c.iter().zip(basis) {
            for (a, bb) in out.iter_mut().zip(b) {
                *a += ci * bb;
            }
        }
        out
    };
    let mut iter = 0_usize;
    let mut verify = 0_usize;
    while iter < opts.max_iter {
        iter += 1;
        if q.len() >= k {
            verify += 1;
            if verify > 150 {
                break;
            }
        }
        let m = v.len();
        let s_mat: Vec<Vec<C>> = (0..m)
            .map(|i| (0..m).map(|j| cdot(&w[i], &v[j])).collect())
            .collect();
        let h_mat: Vec<Vec<C>> = (0..m)
            .map(|i| (0..m).map(|j| cdot(&w[i], &av[j])).collect())
            .collect();
        // G = S^{-1} H, column by column
        let mut g_cols: Vec<Vec<C>> = Vec::with_capacity(m);
        for j in 0..m {
            let col: Vec<C> = (0..m).map(|i| h_mat[i][j]).collect();
            if let Some(c) = csolve_small(&s_mat, &col) {
                g_cols.push(c);
            } else {
                break;
            }
        }
        let petrov = if g_cols.len() == m {
            let g: Vec<Vec<C>> = (0..m)
                .map(|i| (0..m).map(|j| g_cols[j][i]).collect())
                .collect();
            let (t, uu) = complex_schur(&g).ok_or(EigsError::Dense)?;
            let mut order: Vec<usize> = (0..m).collect();
            order.sort_by(|&a, &b| t[a][a].norm().total_cmp(&t[b][b].norm()));
            let mut data = PetrovData {
                order,
                t,
                u: uu,
                y_mat: Vec::new(),
                s_mat,
            };
            let cols: Vec<Vec<C>> = (0..m).map(|i| data.right(i)).collect();
            data.y_mat = (0..m)
                .map(|i| (0..m).map(|j| cols[j][i]).collect())
                .collect();
            Some(data)
        } else {
            None
        };
        // the best Ritz triplet, when the small problem is regular
        let best = petrov.as_ref().and_then(|d| {
            let idx = d.order[0];
            let y = d.right(idx);
            let z = d.left(idx)?;
            let theta_p = d.t[idx][idx];
            let (mut u, mut au) = (combine(&y, &v), combine(&y, &av));
            let (mut wv, mut aw_z) = (combine(&z, &w), combine(&z, &aw));
            let (nu, nw) = (cnrm(&u), cnrm(&wv));
            if nu == 0.0 || nw == 0.0 || nu.is_nan() || nw.is_nan() {
                return None;
            }
            for e in u.iter_mut().chain(au.iter_mut()) {
                *e /= nu;
            }
            for e in wv.iter_mut().chain(aw_z.iter_mut()) {
                *e /= nw;
            }
            let ru: Vec<C> = au.iter().zip(&u).map(|(a, b)| a - theta_p * b).collect();
            let rl: Vec<C> = aw_z
                .iter()
                .zip(&wv)
                .map(|(a, b)| a - theta_p.conj() * b)
                .collect();
            Some((theta_p, u, wv, ru, rl))
        });
        let (Some(petrov), Some((theta_p, u, wv, ru, rl))) = (petrov, best) else {
            // singular small pencil: add random directions
            let mut x = deflate_r(&random_vec(&mut rng), &q, &p);
            let mut y = deflate_l(&random_vec(&mut rng), &q, &p);
            if orthonormal_against(&mut x, &v) && orthonormal_against(&mut y, &w) {
                av.push(op_r(&x, &q, &p, &mut matvecs));
                aw.push(op_l(&y, &q, &p, &mut matvecs));
                v.push(x);
                w.push(y);
            }
            continue;
        };
        let theta = tau + theta_p;
        let (rn, sn) = (cnrm(&ru), cnrm(&rl));
        if q.len() >= k {
            let mut ds: Vec<f64> = lam.iter().map(|l| (l - tau).norm()).collect();
            ds.sort_by(f64::total_cmp);
            if theta_p.norm() >= ds[k - 1] {
                break;
            }
        }
        // restart bases from the Ritz pairs `skip..skip + count` of the current spaces
        let restart = |petrov: &PetrovData,
                       skip: usize,
                       count: usize,
                       q: &[Vec<C>],
                       p: &[Vec<C>],
                       v: &[Vec<C>],
                       w: &[Vec<C>]| {
            let mut nv: Vec<Vec<C>> = Vec::new();
            let mut nw: Vec<Vec<C>> = Vec::new();
            for &idx in petrov.order.iter().skip(skip).take(count) {
                let Some(z) = petrov.left(idx) else { continue };
                let mut x = deflate_r(&combine(&petrov.right(idx), v), q, p);
                let mut yv = deflate_l(&combine(&z, w), q, p);
                if orthonormal_against(&mut x, &nv) && orthonormal_against(&mut yv, &nw) {
                    nv.push(x);
                    nw.push(yv);
                }
            }
            (nv, nw)
        };
        let bound = opts.tol * theta.norm().max(1.0);
        let rho = cdot(&wv, &u);
        if rn <= bound && sn <= bound && rho.norm() > 1e-12 {
            // lock: scale the left vector so that p^H q = 1
            q.push(u.clone());
            p.push(wv.iter().map(|e| e / rho.conj()).collect());
            lam.push(theta);
            unit_pairs.push((u.clone(), wv.clone()));
            let (mut nv, mut nw) = restart(&petrov, 1, min_dim, &q, &p, &v, &w);
            if nv.is_empty() {
                let mut x = deflate_r(&random_vec(&mut rng), &q, &p);
                let mut yv = deflate_l(&random_vec(&mut rng), &q, &p);
                let _ = orthonormal_against(&mut x, &[]);
                let _ = orthonormal_against(&mut yv, &[]);
                nv.push(x);
                nw.push(yv);
            }
            av = nv.iter().map(|x| op_r(x, &q, &p, &mut matvecs)).collect();
            aw = nw.iter().map(|x| op_l(x, &q, &p, &mut matvecs)).collect();
            v = nv;
            w = nw;
            continue;
        }
        // correction equations with the transient pair added to the deflation
        let mut qt = q.clone();
        let mut pt = p.clone();
        qt.push(u.clone());
        pt.push(if rho.norm() > 1e-12 {
            wv.iter().map(|e| e / rho.conj()).collect()
        } else {
            wv.clone()
        });
        let d_t = gram_pq(&pt, &qt);
        let d_th = adjoint_small(&d_t);
        let pi_r = |x: &[C]| oblique(x, &qt, &pt, &d_t);
        let pi_l = |x: &[C]| oblique(x, &pt, &qt, &d_th);
        let right_op = |x: &[C]| -> Vec<C> {
            let px = pi_r(x);
            let mut mv = 0_usize;
            let mut y = run_a(&px, &mut mv);
            for (a, b) in y.iter_mut().zip(&px) {
                *a -= theta * b;
            }
            pi_r(&y)
        };
        let left_op = |x: &[C]| -> Vec<C> {
            let px = pi_l(x);
            let mut mv = 0_usize;
            let mut y = run_t(&px, &mut mv);
            for (a, b) in y.iter_mut().zip(&px) {
                *a -= theta.conj() * b;
            }
            pi_l(&y)
        };
        let ident = |x: &[C]| -> Vec<C> { x.to_vec() };
        let rhs_r = pi_r(&ru.iter().map(|e| -e).collect::<Vec<_>>());
        let rhs_l = pi_l(&rl.iter().map(|e| -e).collect::<Vec<_>>());
        let mut t_r = gmres_c(&right_op, &ident, &rhs_r, opts.inner_iters, opts.inner_tol);
        let mut t_l = gmres_c(&left_op, &ident, &rhs_l, opts.inner_iters, opts.inner_tol);
        matvecs += 2 * opts.inner_iters;
        t_r = deflate_r(&pi_r(&t_r), &q, &p);
        t_l = deflate_l(&pi_l(&t_l), &q, &p);
        let ok_r = orthonormal_against(&mut t_r, &v);
        let ok_l = orthonormal_against(&mut t_l, &w);
        if !(ok_r && ok_l) {
            let mut x = deflate_r(&random_vec(&mut rng), &q, &p);
            let mut y = deflate_l(&random_vec(&mut rng), &q, &p);
            if orthonormal_against(&mut x, &v) && orthonormal_against(&mut y, &w) {
                t_r = x;
                t_l = y;
            } else {
                continue;
            }
        }
        av.push(op_r(&t_r, &q, &p, &mut matvecs));
        aw.push(op_l(&t_l, &q, &p, &mut matvecs));
        v.push(t_r);
        w.push(t_l);
        if v.len() > max_dim {
            // thick restart onto the best Ritz pairs plus the newest directions
            let (mut nv, mut nw) = restart(&petrov, 0, min_dim, &q, &p, &v, &w);
            let (mut x, mut yv) = (v[v.len() - 1].clone(), w[w.len() - 1].clone());
            if orthonormal_against(&mut x, &nv) && orthonormal_against(&mut yv, &nw) {
                nv.push(x);
                nw.push(yv);
            }
            if !nv.is_empty() {
                av = nv.iter().map(|x| op_r(x, &q, &p, &mut matvecs)).collect();
                aw = nw.iter().map(|x| op_l(x, &q, &p, &mut matvecs)).collect();
                v = nv;
                w = nw;
            }
        }
    }
    // output: nearest first, true residuals
    let kq = q.len();
    let mut order: Vec<usize> = (0..kq).collect();
    order.sort_by(|&a, &b| (lam[a] - tau).norm().total_cmp(&(lam[b] - tau).norm()));
    order.truncate(k);
    let mut out = TwoSidedResult {
        values: Vec::new(),
        right_re: Vec::new(),
        right_im: Vec::new(),
        left_re: Vec::new(),
        left_im: Vec::new(),
        right_residuals: Vec::new(),
        left_residuals: Vec::new(),
        condition: Vec::new(),
        converged: kq >= k,
        iterations: iter,
        matvecs,
    };
    for &i in &order {
        let (u, wv) = &unit_pairs[i];
        let big = u
            .iter()
            .copied()
            .max_by(|a, b| a.norm().total_cmp(&b.norm()))
            .unwrap_or(zero);
        let phase = if big.norm() == 0.0 {
            C::new(1.0, 0.0)
        } else {
            big.conj() / big.norm()
        };
        let us: Vec<C> = u.iter().map(|e| e * phase).collect();
        let c = cdot(wv, &us);
        let wphase = if c.norm() == 0.0 {
            C::new(1.0, 0.0)
        } else {
            c / c.norm()
        };
        let ws: Vec<C> = wv.iter().map(|e| e * wphase).collect();
        let l = lam[i];
        let mut mv = 0_usize;
        let au = run_a(&us, &mut mv);
        let aws = run_t(&ws, &mut mv);
        let res_r = au
            .iter()
            .zip(&us)
            .map(|(a, b)| (a - l * b).norm_sqr())
            .sum::<f64>()
            .sqrt();
        let res_l = aws
            .iter()
            .zip(&ws)
            .map(|(a, b)| (a - l.conj() * b).norm_sqr())
            .sum::<f64>()
            .sqrt();
        out.values.push((l.re, l.im));
        out.right_re.push(us.iter().map(|e| e.re).collect());
        out.right_im.push(us.iter().map(|e| e.im).collect());
        out.left_re.push(ws.iter().map(|e| e.re).collect());
        out.left_im.push(ws.iter().map(|e| e.im).collect());
        out.right_residuals.push(res_r);
        out.left_residuals.push(res_l);
        out.condition.push(1.0 / cdot(&ws, &us).norm().max(1e-300));
    }
    Ok(out)
}
