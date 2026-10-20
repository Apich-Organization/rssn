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
//! * [`eigs_jd_sym`]: Jacobi-Davidson for **interior** eigenvalues of a
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
//! Limitations: Jacobi-Davidson is implemented for symmetric operators only
//! (interior eigenvalues of non-symmetric problems are reached by
//! shift-invert with [`eigs_krylov_schur`]); `B` must be positive definite
//! in the symmetric generalised solver.

use num_complex::Complex64;

use crate::kernels::dense::{self, Mat};
use crate::kernels::eigs::{EigsError, EigsOptions, EigsResult, Which, key_of};
use crate::kernels::krylov::gmres;
use crate::kernels::random::Rng;

type C = Complex64;

fn rdot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn rnrm(a: &[f64]) -> f64 {
    rdot(a, a).sqrt()
}

fn cdot(a: &[C], b: &[C]) -> C {
    a.iter().zip(b).map(|(x, y)| x.conj() * y).sum()
}

fn cnrm(a: &[C]) -> f64 {
    a.iter().map(Complex64::norm_sqr).sum::<f64>().sqrt()
}

fn random_vec(n: usize, rng: &mut Rng) -> Vec<f64> {
    (0..n).map(|_| rng.normal()).collect()
}

// ------------------------------------------------------- complex Schur form

/// Complex Givens rotation `[[c, s], [-conj(s), c]]` (`c` real) mapping
/// `(a, b)` to `(r, 0)`.
fn givens(a: C, b: C) -> (f64, C) {
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
fn rot_rows(h: &mut [Vec<C>], k: usize, c: f64, s: C, c0: usize) {
    let ncols = h[k].len();
    for j in c0..ncols {
        let (x, y) = (h[k][j], h[k + 1][j]);
        h[k][j] = x * c + s * y;
        h[k + 1][j] = -s.conj() * x + y * c;
    }
}

/// `columns k, k+1 <- columns * G^H`, rows `0..r1`.
fn rot_cols(h: &mut [Vec<C>], k: usize, c: f64, s: C, r1: usize) {
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
pub fn complex_schur(a: &[Vec<C>]) -> Option<(Vec<Vec<C>>, Vec<Vec<C>>)> {
    let n = a.len();
    let mut h: Vec<Vec<C>> = a.to_vec();
    let mut q: Vec<Vec<C>> =
        (0..n).map(|i| (0..n).map(|j| if i == j { C::new(1.0, 0.0) } else { C::new(0.0, 0.0) }).collect()).collect();
    // Householder reduction to Hessenberg form
    for k in 0..n.saturating_sub(2) {
        let mut x: Vec<C> = (k + 1..n).map(|i| h[i][k]).collect();
        let xn = cnrm(&x);
        if xn == 0.0 {
            continue;
        }
        let phase = if x[0].norm() == 0.0 { C::new(1.0, 0.0) } else { x[0] / x[0].norm() };
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
        for row in h.iter_mut() {
            let s: C = (0..x.len()).map(|i| row[k + 1 + i] * x[i]).sum();
            for i in 0..x.len() {
                row[k + 1 + i] -= 2.0 * s * x[i].conj();
            }
        }
        for row in q.iter_mut() {
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
        let mu = if iter % 10 == 0 {
            a22 + C::new(h[hi][hi - 1].norm(), 0.0)
        } else {
            let tr_half = (a11 + a22) * 0.5;
            let det = a11 * a22 - a12 * a21;
            let disc = (tr_half * tr_half - det).sqrt();
            let (l1, l2) = (tr_half + disc, tr_half - disc);
            if (l1 - a22).norm() < (l2 - a22).norm() { l1 } else { l2 }
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
fn swap_schur(t: &mut [Vec<C>], q: &mut [Vec<C>], i: usize) {
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
pub fn reorder_schur(t: &mut [Vec<C>], q: &mut [Vec<C>], count: usize, key: &dyn Fn(C) -> f64) {
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

fn apply_complex<A: Fn(&[f64], &mut [f64])>(apply: &A, v: &[C], out: &mut [C], re: &mut [f64], im: &mut [f64]) {
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
        im.iter_mut().for_each(|e| *e = 0.0);
    }
    for i in 0..n {
        out[i] = C::new(re[i], im[i]);
    }
}

fn orth_complex(w: &mut [C], v: &[Vec<C>], cnt: usize) -> Vec<C> {
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

fn random_complex_orth(n: usize, v: &[Vec<C>], cnt: usize, seed: u64) -> Vec<C> {
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
    let m = (if opts.ncv == 0 { (2 * k + 1).max(20) } else { opts.ncv }).max(k + 2).min(n);
    let zero = C::new(0.0, 0.0);
    let mut rng = Rng::new(opts.seed);
    let mut v0: Vec<C> = (0..n).map(|_| C::new(rng.normal(), 0.0)).collect();
    let s0 = cnrm(&v0);
    v0.iter_mut().for_each(|e| *e /= s0);
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
                    let r = random_complex_orth(n, &v, j + 1, opts.seed.wrapping_add(j as u64 + 31));
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
                let d = if d.norm() < 1e-14 * lam.norm().max(1.0) { C::new(1e-14 * lam.norm().max(1.0), 0.0) } else { d };
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
                let big = x.iter().copied().max_by(|p, r| p.norm().total_cmp(&r.norm())).unwrap_or(zero);
                let phase = if big.norm() == 0.0 { C::new(1.0, 0.0) } else { big.conj() / big.norm() };
                out.values.push((t[i][i].re, t[i][i].im));
                out.vectors_re.push(x.iter().map(|e| (e * phase / xn).re).collect());
                out.vectors_im.push(x.iter().map(|e| (e * phase / xn).im).collect());
                out.residuals.push(if m >= n { 0.0 } else { res[i] });
            }
            return Ok(out);
        }
        // contraction to kk Schur vectors (never splitting a conjugate pair)
        let mut keep = (k + (m - k) / 2).clamp(k + 1, m - 1);
        if keep < m
            && t[keep - 1][keep - 1].im.abs() > 1e-12 * t[keep - 1][keep - 1].norm().max(1.0)
            && (t[keep - 1][keep - 1] - t[keep][keep].conj()).norm() < 1e-8 * t[keep][keep].norm().max(1.0)
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
    bmul: Option<&'a dyn Fn(&[f64], &mut [f64])>,
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
    fn b_apply(&self, x: &[f64]) -> Vec<f64> {
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
    fn extend(&self, st: &mut BlockState, cands: Vec<Vec<f64>>, limit: usize) -> usize {
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
                    c.iter_mut().for_each(|e| *e /= after);
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
    fn contract(&self, st: &mut BlockState, s: &Mat, cols: &[usize]) {
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

    fn run(&self, k: usize, which: Which, opts: &EigsOptions) -> Result<BlockOut, EigsError> {
        let n = self.n;
        let bs = self.block.max(1);
        let m_max = (if opts.ncv == 0 { (2 * k + 2 * bs).max(24).max(4 * bs) } else { opts.ncv })
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
                let cands: Vec<Vec<f64>> = (last_block_start..st.v.len()).map(|i| st.av[i].clone()).collect();
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
                    let s = 0.5 * (t.at(i, j) + t.at(j, i));
                    t.set(i, j, s);
                    t.set(j, i, s);
                }
            }
            let (theta, s) = dense::eigen_symmetric(&t).map_err(|_| EigsError::Dense)?;
            let mut idx: Vec<usize> = (0..m).collect();
            idx.sort_by(|&a, &b| key_of(which, (theta[a], 0.0)).total_cmp(&key_of(which, (theta[b], 0.0))));
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
                    residuals: if m >= n { vec![0.0; kk] } else { res.into_iter().take(kk).collect() },
                    converged: all_conv || m >= n,
                    restarts,
                    matvecs: st.matvecs,
                });
            }
            // thick restart: keep Ritz vectors, continue with the residual block
            let keep = (k + (m_max - k) / 2).min(m_max.saturating_sub(bs + 1)).max(k).min(m - 1);
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

fn pack(out: BlockOut, map: &dyn Fn(f64) -> f64, n: usize) -> EigsResult {
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
    let eng = BlockEngine { n, op: &apply, bmul: None, block };
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
    let eng = BlockEngine { n, op: &op, bmul: Some(&apply_b), block: opts.block };
    let out = eng.run(k, Which::LargestMagnitude, &opts.eigs)?;
    let sigma = opts.sigma;
    let mut r = pack(out, &|t| sigma + 1.0 / t, n);
    gen_residuals(&apply_a, &apply_b, &mut r);
    Ok(r)
}

/// General generalised eigenproblem `A x = lambda B x` (`A` arbitrary real,
/// `B` real and nonsingular): the `k` eigenvalues closest to `opts.sigma`
/// by Krylov-Schur on `(A - sigma B)^{-1} B`. Eigenvalues may be complex.
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
fn orthonormalise_pair(v: &mut Vec<Vec<f64>>, w: &mut Vec<Vec<f64>>) {
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
            x.iter_mut().for_each(|e| *e /= nx);
            y.iter_mut().for_each(|e| *e /= nx);
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
    precond: Option<&dyn Fn(&[f64], &mut [f64])>,
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
    t0.iter_mut().for_each(|e| *e /= tn);
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
                    let mut y: Vec<f64> = (0..m).map(|i| 0.5 * (h2.at(i, j) + h2.at(j, i))).collect();
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
                        let s = 0.5 * (mt.at(i, j) + mt.at(j, i));
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
                    let s = 0.5 * (hs.at(i, j) + hs.at(j, i));
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
        u.iter_mut().for_each(|e| *e /= un);
        wu.iter_mut().for_each(|e| *e /= un);
        // A u = wu + target u; Rayleigh quotient
        let theta = target + rdot(&u, &wu);
        let r: Vec<f64> = wu.iter().zip(&u).map(|(a, b)| a - (theta - target) * b).collect();
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
                t.iter_mut().for_each(|e| *e /= tn);
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
        let (kinv_q, small): (Vec<Vec<f64>>, Option<Mat>) = if let Some(pc) = precond {
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
        } else {
            (Vec::new(), None)
        };
        let precond_apply = |x: &[f64], y: &mut [f64]| match (precond, small.as_ref()) {
            (Some(pc), Some(sm)) => {
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
            }
            _ => {
                y.copy_from_slice(x);
                project(y, &qt);
            }
        };
        let mut rhs: Vec<f64> = r.iter().map(|e| -e).collect();
        project(&mut rhs, &qt);
        let sol = gmres(proj_apply, precond_apply, &rhs, None, opts.inner_iters, opts.inner_tol, opts.inner_iters);
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
        if !(tn > 1e-12) || !tn.is_finite() {
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
        t.iter_mut().for_each(|e| *e /= tn);
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
    order.sort_by(|&a, &b| (evals[a] - target).abs().total_cmp(&(evals[b] - target).abs()));
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
