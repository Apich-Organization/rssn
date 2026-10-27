//! # Further Krylov solvers
//!
//! Matrix-free iterative solvers complementing [`crate::kernels::krylov`]:
//!
//! * [`minres`]: MINRES for symmetric (possibly indefinite) systems with a
//!   symmetric positive definite preconditioner (Paige and Saunders);
//! * [`bicgstab_l`]: BiCGSTAB(`l`) of Sleijpen and Fokkema, right
//!   preconditioned (`l = 1` is BiCGSTAB, larger `l` copes with complex
//!   spectra);
//! * [`idr_s`]: IDR(`s`) of Sonneveld and van Gijzen with the
//!   biorthogonal residual-smoothing-free variant and right preconditioning;
//! * [`fgmres`]: flexible restarted GMRES whose preconditioner may change
//!   from one iteration to the next (for example an inner iterative solve).
//!
//! `apply(x, y)` writes `A x` into `y`; the preconditioner closures write an
//! approximation of `M^{-1} r` into their second argument.

use crate::kernels::krylov::IterResult;
use crate::kernels::random::Rng;

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn nrm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}

fn axpy(alpha: f64, x: &[f64], y: &mut [f64]) {
    for (yi, xi) in y.iter_mut().zip(x) {
        *yi += alpha * xi;
    }
}

fn residual<A: Fn(&[f64], &mut [f64])>(apply: &A, b: &[f64], x: &[f64]) -> Vec<f64> {
    let mut ax = vec![0.0; b.len()];
    apply(x, &mut ax);
    b.iter().zip(&ax).map(|(bi, ai)| bi - ai).collect()
}

/// MINRES for a symmetric (indefinite allowed) operator with an SPD
/// preconditioner `precond`.
///
/// The reported `residual` is the true relative residual `|b - A x| / |b|`
/// of the returned iterate; the iteration stops when the recurrence
/// estimate (preconditioned residual norm) drops below `tol` relative to
/// its initial value.
pub fn minres<A, M>(
    apply: A,
    precond: M,
    b: &[f64],
    x0: Option<&[f64]>,
    tol: f64,
    max_iter: usize,
) -> IterResult
where
    A: Fn(&[f64], &mut [f64]),
    M: Fn(&[f64], &mut [f64]),
{
    let n = b.len();
    let bn = nrm(b).max(f64::MIN_POSITIVE);
    let mut x = x0.map_or_else(|| vec![0.0; n], <[f64]>::to_vec);
    let r0 = residual(&apply, b, &x);
    let mut y = vec![0.0; n];
    precond(&r0, &mut y);
    let beta1sq = dot(&r0, &y);
    if beta1sq <= 0.0 || nrm(&r0) / bn <= tol {
        let res = nrm(&r0) / bn;
        return IterResult { x, residual: res, iterations: 0, converged: res <= tol };
    }
    let beta1 = beta1sq.sqrt();
    let mut r1 = r0.clone();
    let mut r2 = r0;
    let (mut oldb, mut beta) = (0.0, beta1);
    let (mut dbar, mut epsln, mut phibar) = (0.0, 0.0, beta1);
    let (mut cs, mut sn) = (-1.0, 0.0);
    let mut w = vec![0.0; n];
    let mut w2 = vec![0.0; n];
    let mut v = vec![0.0; n];
    let mut itn = 0;
    while itn < max_iter {
        itn += 1;
        let s = 1.0 / beta;
        for i in 0..n {
            v[i] = s * y[i];
        }
        apply(&v, &mut y);
        if itn >= 2 {
            axpy(-beta / oldb, &r1, &mut y);
        }
        let alfa = dot(&v, &y);
        axpy(-alfa / beta, &r2, &mut y);
        r1 = std::mem::replace(&mut r2, y.clone());
        precond(&r2, &mut y);
        oldb = beta;
        let bsq = dot(&r2, &y);
        beta = bsq.max(0.0).sqrt();
        let oldeps = epsln;
        let delta = cs * dbar + sn * alfa;
        let gbar = sn * dbar - cs * alfa;
        epsln = sn * beta;
        dbar = -cs * beta;
        let gamma = gbar.hypot(beta).max(f64::MIN_POSITIVE);
        cs = gbar / gamma;
        sn = beta / gamma;
        let phi = cs * phibar;
        phibar *= sn;
        let w1 = std::mem::replace(&mut w2, std::mem::take(&mut w));
        w = (0..n).map(|i| (v[i] - oldeps * w1[i] - delta * w2[i]) / gamma).collect();
        axpy(phi, &w, &mut x);
        if phibar / beta1 <= tol || beta == 0.0 {
            break;
        }
    }
    let res = nrm(&residual(&apply, b, &x)) / bn;
    IterResult { x, residual: res, iterations: itn, converged: res <= 10.0 * tol.max(1e-15) }
}

/// `BiCGSTAB`(`l`) with right preconditioning (`l >= 1`).
///
/// `l = 1` is classical `BiCGSTAB`; `l = 2` or more smooths the residual
/// polynomial with a degree-`l` minimal-residual step, which avoids the
/// breakdown of `BiCGSTAB` on matrices with strongly complex spectra. The
/// iteration count is the number of outer cycles; each costs `2 l`
/// matrix-vector products.
pub fn bicgstab_l<A, M>(
    apply: A,
    precond: M,
    l: usize,
    b: &[f64],
    x0: Option<&[f64]>,
    tol: f64,
    max_iter: usize,
) -> IterResult
where
    A: Fn(&[f64], &mut [f64]),
    M: Fn(&[f64], &mut [f64]),
{
    let l = l.max(1);
    let n = b.len();
    let bn = nrm(b).max(f64::MIN_POSITIVE);
    let x0v = x0.map_or_else(|| vec![0.0; n], <[f64]>::to_vec);
    let mut tmp = vec![0.0; n];
    // operator B v = A M^{-1} v; unknown correction y with x = x0 + M^{-1} y
    let op = |v: &[f64], out: &mut [f64], tmp: &mut Vec<f64>| {
        precond(v, tmp);
        apply(tmp, out);
    };
    let r_init = residual(&apply, b, &x0v);
    let mut xc = vec![0.0; n]; // correction in the preconditioned space
    let rt = r_init.clone();
    let mut r: Vec<Vec<f64>> = vec![vec![0.0; n]; l + 1];
    let mut u: Vec<Vec<f64>> = vec![vec![0.0; n]; l + 1];
    r[0] = r_init;
    let (mut rho0, mut alpha, mut omega) = (1.0, 0.0, 1.0);
    let mut it = 0;
    let mut rn = nrm(&r[0]);
    while rn / bn > tol && it < max_iter {
        it += 1;
        rho0 *= -omega;
        for j in 0..l {
            let rho1 = dot(&rt, &r[j]);
            if rho0 == 0.0 {
                return finish(&apply, &precond, b, &x0v, &xc, it);
            }
            let beta = alpha * rho1 / rho0;
            rho0 = rho1;
            for i in 0..=j {
                for q in 0..n {
                    u[i][q] = r[i][q] - beta * u[i][q];
                }
            }
            let uj = u[j].clone();
            let mut out = vec![0.0; n];
            op(&uj, &mut out, &mut tmp);
            u[j + 1] = out;
            let gamma = dot(&u[j + 1], &rt);
            if gamma == 0.0 {
                return finish(&apply, &precond, b, &x0v, &xc, it);
            }
            alpha = rho0 / gamma;
            for i in 0..=j {
                let ui1 = u[i + 1].clone();
                axpy(-alpha, &ui1, &mut r[i]);
            }
            let rj = r[j].clone();
            let mut out = vec![0.0; n];
            op(&rj, &mut out, &mut tmp);
            r[j + 1] = out;
            let u0 = u[0].clone();
            axpy(alpha, &u0, &mut xc);
        }
        // minimal residual part
        let mut tau = vec![vec![0.0; l + 1]; l + 1];
        let mut sigma = vec![0.0; l + 1];
        let mut gp = vec![0.0; l + 1];
        for j in 1..=l {
            for i in 1..j {
                tau[i][j] = dot(&r[j], &r[i]) / sigma[i];
                let ri = r[i].clone();
                axpy(-tau[i][j], &ri, &mut r[j]);
            }
            sigma[j] = dot(&r[j], &r[j]);
            gp[j] = if sigma[j] == 0.0 { 0.0 } else { dot(&r[0], &r[j]) / sigma[j] };
        }
        let mut g = vec![0.0; l + 1];
        g[l] = gp[l];
        omega = g[l];
        for j in (1..l).rev() {
            let mut s = gp[j];
            for i in j + 1..=l {
                s -= tau[j][i] * g[i];
            }
            g[j] = s;
        }
        let mut gpp = vec![0.0; l + 1];
        for j in 1..l {
            let mut s = g[j + 1];
            for i in j + 1..l {
                s += tau[j][i] * g[i + 1];
            }
            gpp[j] = s;
        }
        let r0c = r[0].clone();
        axpy(g[1], &r0c, &mut xc);
        let rl = r[l].clone();
        axpy(-gp[l], &rl, &mut r[0]);
        let ul = u[l].clone();
        axpy(-g[l], &ul, &mut u[0]);
        for j in 1..l {
            let uj = u[j].clone();
            axpy(-g[j], &uj, &mut u[0]);
            let rj = r[j].clone();
            axpy(gpp[j], &rj, &mut xc);
            axpy(-gp[j], &rj, &mut r[0]);
        }
        rn = nrm(&r[0]);
        if !rn.is_finite() {
            break;
        }
    }
    finish(&apply, &precond, b, &x0v, &xc, it)
}

fn finish<A, M>(apply: &A, precond: &M, b: &[f64], x0: &[f64], xc: &[f64], it: usize) -> IterResult
where
    A: Fn(&[f64], &mut [f64]),
    M: Fn(&[f64], &mut [f64]),
{
    let n = b.len();
    let mut corr = vec![0.0; n];
    precond(xc, &mut corr);
    let x: Vec<f64> = x0.iter().zip(&corr).map(|(a, c)| a + c).collect();
    let res = nrm(&residual(apply, b, &x)) / nrm(b).max(f64::MIN_POSITIVE);
    IterResult { x, residual: res, iterations: it, converged: res.is_finite() && res < 1.0 }
}

/// Parameters of [`idr_s_with`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IdrOptions {
    /// Dimension `s` of the shadow space.
    pub s: usize,
    /// Relative residual tolerance.
    pub tol: f64,
    /// Maximal number of matrix-vector products.
    pub max_iter: usize,
    /// Seed of the pseudo-random shadow vectors.
    pub seed: u64,
}

/// IDR(`s`) with right preconditioning.
///
/// The shadow space is `s` fixed pseudo-random orthonormal vectors
/// (`seed`); the method needs only about `n + n/s` matrix-vector products
/// and `s + 4` vectors of storage, and approaches the speed of full GMRES
/// for moderate `s` on nonsymmetric systems. `iterations` counts
/// matrix-vector products.
pub fn idr_s_with<A, M>(
    apply: A,
    precond: M,
    b: &[f64],
    x0: Option<&[f64]>,
    opts: &IdrOptions,
) -> IterResult
where
    A: Fn(&[f64], &mut [f64]),
    M: Fn(&[f64], &mut [f64]),
{
    let IdrOptions { s, tol, max_iter, seed } = *opts;
    let n = b.len();
    let s = s.clamp(1, n.max(1));
    let bn = nrm(b).max(f64::MIN_POSITIVE);
    let mut x = x0.map_or_else(|| vec![0.0; n], <[f64]>::to_vec);
    let mut r = residual(&apply, b, &x);
    let mut mv = 1usize;
    if nrm(&r) / bn <= tol {
        let res = nrm(&r) / bn;
        return IterResult { x, residual: res, iterations: 0, converged: true };
    }
    // orthonormal shadow vectors
    let mut rng = Rng::new(seed);
    let mut p: Vec<Vec<f64>> = Vec::with_capacity(s);
    for _ in 0..s {
        let mut v: Vec<f64> = (0..n).map(|_| rng.normal()).collect();
        for q in &p {
            let c = dot(q, &v);
            axpy(-c, q, &mut v);
        }
        let nv = nrm(&v);
        for e in &mut v {
            *e /= nv;
        }
        p.push(v);
    }
    let mut g: Vec<Vec<f64>> = vec![vec![0.0; n]; s];
    let mut u: Vec<Vec<f64>> = vec![vec![0.0; n]; s];
    let mut m = vec![vec![0.0; s]; s];
    for i in 0..s {
        m[i][i] = 1.0;
    }
    let mut om = 1.0;
    let mut f: Vec<f64> = p.iter().map(|pi| dot(pi, &r)).collect();
    let mut t = vec![0.0; n];
    let mut vv = vec![0.0; n];
    'outer: while nrm(&r) / bn > tol && mv < max_iter {
        for k in 0..s {
            // solve the lower triangular system M[k..,k..] c = f[k..]
            let mut c = vec![0.0; s - k];
            for i in k..s {
                let mut sum = f[i];
                for j in k..i {
                    sum -= m[i][j] * c[j - k];
                }
                if m[i][i] == 0.0 {
                    break 'outer;
                }
                c[i - k] = sum / m[i][i];
            }
            let mut v = r.clone();
            for j in k..s {
                axpy(-c[j - k], &g[j], &mut v);
            }
            precond(&v, &mut vv);
            let mut uk = vec![0.0; n];
            for j in k..s {
                axpy(c[j - k], &u[j], &mut uk);
            }
            axpy(om, &vv, &mut uk);
            apply(&uk, &mut t);
            mv += 1;
            let mut gk = t.clone();
            for j in 0..k {
                let alpha = dot(&p[j], &gk) / m[j][j];
                let gj = g[j].clone();
                axpy(-alpha, &gj, &mut gk);
                let uj = u[j].clone();
                axpy(-alpha, &uj, &mut uk);
            }
            for i in k..s {
                m[i][k] = dot(&p[i], &gk);
            }
            if m[k][k] == 0.0 {
                break 'outer;
            }
            let beta = f[k] / m[k][k];
            axpy(-beta, &gk, &mut r);
            axpy(beta, &uk, &mut x);
            g[k] = gk;
            u[k] = uk;
            for i in k + 1..s {
                f[i] -= beta * m[i][k];
            }
            if nrm(&r) / bn <= tol || mv >= max_iter {
                break 'outer;
            }
        }
        // dimension reduction step
        precond(&r, &mut vv);
        apply(&vv, &mut t);
        mv += 1;
        let (tt, tr) = (dot(&t, &t), dot(&t, &r));
        if tt == 0.0 {
            break;
        }
        om = tr / tt;
        let rho = tr.abs() / (tt.sqrt() * nrm(&r)).max(f64::MIN_POSITIVE);
        if rho < 0.7 && rho > 0.0 {
            om *= 0.7 / rho;
        }
        axpy(-om, &t, &mut r);
        axpy(om, &vv, &mut x);
        for i in 0..s {
            f[i] = dot(&p[i], &r);
        }
    }
    let res = nrm(&residual(&apply, b, &x)) / bn;
    IterResult { x, residual: res, iterations: mv, converged: res <= 10.0 * tol.max(1e-15) }
}

/// [`idr_s_with`] with the options passed as separate arguments.
#[deprecated(note = "use `idr_s_with` with `IdrOptions`")]
#[allow(clippy::too_many_arguments)] // legacy signature kept as a thin wrapper
pub fn idr_s<A, M>(
    apply: A,
    precond: M,
    s: usize,
    b: &[f64],
    x0: Option<&[f64]>,
    tol: f64,
    max_iter: usize,
    seed: u64,
) -> IterResult
where
    A: Fn(&[f64], &mut [f64]),
    M: Fn(&[f64], &mut [f64]),
{
    idr_s_with(apply, precond, b, x0, &IdrOptions { s, tol, max_iter, seed })
}

/// Flexible GMRES(`restart`) with a variable right preconditioner.
///
/// Flexible GMRES(`restart`): right-preconditioned GMRES whose
/// preconditioner `precond` may vary between applications (for instance
/// an inner Krylov solve or a nonlinear smoother). The preconditioned
/// directions are stored, costing one extra vector per iteration.
pub fn fgmres<A, M>(
    apply: A,
    precond: M,
    b: &[f64],
    x0: Option<&[f64]>,
    restart: usize,
    tol: f64,
    max_iter: usize,
) -> IterResult
where
    A: Fn(&[f64], &mut [f64]),
    M: Fn(&[f64], &mut [f64]),
{
    let n = b.len();
    let m = restart.clamp(1, n.max(1));
    let bn = nrm(b).max(f64::MIN_POSITIVE);
    let mut x = x0.map_or_else(|| vec![0.0; n], <[f64]>::to_vec);
    let mut total = 0;
    let mut res;
    loop {
        let r = residual(&apply, b, &x);
        let beta = nrm(&r);
        res = beta / bn;
        if res <= tol || total >= max_iter || beta == 0.0 {
            break;
        }
        let mut v: Vec<Vec<f64>> = vec![r.iter().map(|e| e / beta).collect()];
        let mut z: Vec<Vec<f64>> = Vec::new();
        let mut h = vec![vec![0.0; m]; m + 1];
        let (mut cs, mut sn) = (vec![0.0; m], vec![0.0; m]);
        let mut g = vec![0.0; m + 1];
        g[0] = beta;
        let mut used = 0;
        for j in 0..m {
            let mut zj = vec![0.0; n];
            precond(&v[j], &mut zj);
            let mut w = vec![0.0; n];
            apply(&zj, &mut w);
            z.push(zj);
            total += 1;
            for i in 0..=j {
                h[i][j] = dot(&w, &v[i]);
                let vi = v[i].clone();
                axpy(-h[i][j], &vi, &mut w);
            }
            // reorthogonalise once for robustness
            for i in 0..=j {
                let c = dot(&w, &v[i]);
                h[i][j] += c;
                let vi = v[i].clone();
                axpy(-c, &vi, &mut w);
            }
            h[j + 1][j] = nrm(&w);
            if h[j + 1][j] > 0.0 {
                let hn = h[j + 1][j];
                v.push(w.iter().map(|e| e / hn).collect());
            } else {
                v.push(vec![0.0; n]);
            }
            for i in 0..j {
                let t = cs[i] * h[i][j] + sn[i] * h[i + 1][j];
                h[i + 1][j] = -sn[i] * h[i][j] + cs[i] * h[i + 1][j];
                h[i][j] = t;
            }
            let d = h[j][j].hypot(h[j + 1][j]);
            if d == 0.0 {
                used = j;
                break;
            }
            cs[j] = h[j][j] / d;
            sn[j] = h[j + 1][j] / d;
            h[j][j] = d;
            h[j + 1][j] = 0.0;
            g[j + 1] = -sn[j] * g[j];
            g[j] *= cs[j];
            used = j + 1;
            if g[j + 1].abs() / bn <= tol || total >= max_iter {
                break;
            }
        }
        let mut y = vec![0.0; used];
        for i in (0..used).rev() {
            let mut s = g[i];
            for k in i + 1..used {
                s -= h[i][k] * y[k];
            }
            y[i] = s / h[i][i];
        }
        for i in 0..used {
            axpy(y[i], &z[i], &mut x);
        }
        if used == 0 {
            break;
        }
    }
    let res = nrm(&residual(&apply, b, &x)) / bn;
    IterResult { x, residual: res, iterations: total, converged: res <= tol.max(1e-15) * 10.0 }
}
