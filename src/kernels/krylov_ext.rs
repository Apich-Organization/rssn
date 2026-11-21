//! MINRES, BiCGSTAB(l), IDR(s) and flexible GMRES.
//!
//! # Further Krylov solvers
//!
//! Matrix-free iterative solvers complementing [`crate::kernels::krylov`]:
//!
//! * [`minres`]: MINRES for symmetric (possibly indefinite) systems with a
//!   symmetric positive definite preconditioner (Paige and Saunders);
//! * [`bicgstab_l`]: `BiCGSTAB`(`l`) of Sleijpen and Fokkema, right
//!   preconditioned (`l = 1` is `BiCGSTAB`, larger `l` copes with complex
//!   spectra);
//! * [`idr_s`]: IDR(`s`) of Sonneveld and van Gijzen with the
//!   biorthogonal residual-smoothing-free variant and right preconditioning;
//! * [`fgmres`]: flexible restarted GMRES whose preconditioner may change
//!   from one iteration to the next (for example an inner iterative solve).
//!
//! `apply(x, y)` writes `A x` into `y`; the preconditioner closures write an
//! approximation of `M^{-1} r` into their second argument.

use std::alloc::Layout;

use crate::kernels::krylov::IterResult;
use crate::kernels::krylov::residual_ratio;
use crate::kernels::random::Rng;

fn dot(
    a: &[f64],
    b: &[f64],
) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn nrm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}

fn axpy(
    alpha: f64,
    x: &[f64],
    y: &mut [f64],
) {
    for (yi, xi) in y.iter_mut().zip(x) {
        *yi += alpha * xi;
    }
}

fn residual<A: Fn(&[f64], &mut [f64])>(
    apply: &A,
    b: &[f64],
    x: &[f64],
) -> Vec<f64> {
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
        return IterResult {
            x,
            residual: res,
            iterations: 0,
            converged: res <= tol,
        };
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
        w = (0..n)
            .map(|i| (v[i] - oldeps * w1[i] - delta * w2[i]) / gamma)
            .collect();
        axpy(phi, &w, &mut x);
        if phibar / beta1 <= tol || beta == 0.0 {
            break;
        }
    }
    let res = nrm(&residual(&apply, b, &x)) / bn;
    IterResult {
        x,
        residual: res,
        iterations: itn,
        converged: res <= 10.0 * tol.max(1e-15),
    }
}

/// `BiCGSTAB`(`l`) with right preconditioning (`l >= 1`).
///
/// `l = 1` is classical `BiCGSTAB`; `l = 2` or more smooths the residual
/// polynomial with a degree-`l` minimal-residual step, which avoids the
/// breakdown of `BiCGSTAB` on matrices with strongly complex spectra. The
/// iteration count is the number of outer cycles started. Each complete cycle
/// uses `2 l` matrix-vector products, plus the initial and final verification
/// for the solve.
///
/// Convergence is verified against the original system using the requested
/// nonnegative, finite tolerance. With a zero right-hand side, the reported
/// relative residual is zero for an exact solution and infinite otherwise.
/// Both callbacks must overwrite every output entry with finite values.
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
    let mut x = x0
        .filter(|v| v.len() == n)
        .map_or_else(|| vec![0.0; n], <[f64]>::to_vec);
    if x0.is_some_and(|v| v.len() != n)
        || !tol.is_finite()
        || tol < 0.0
        || b.iter().chain(&x).any(|v| !v.is_finite())
    {
        return bicgstab_l_report(x, f64::INFINITY, 0, tol);
    }

    let mut image = vec![f64::NAN; n];
    if n != 0 {
        apply(&x, &mut image);
    }
    let r_init: Vec<_> = b.iter().zip(&image).map(|(b, a)| b - a).collect();
    let bn = nrm(b).max(f64::MIN_POSITIVE);
    let mut rn = nrm(&r_init);
    let initial = || residual_ratio(b, &image).unwrap_or(f64::INFINITY);
    if max_iter == 0 || !rn.is_finite() || !bn.is_finite() || rn / bn <= tol {
        return bicgstab_l_report(x, initial(), 0, tol);
    }
    let Some(width) = l.checked_add(1) else {
        return bicgstab_l_report(x, initial(), 0, tol);
    };
    let Some(small_matrix_len) = width
        .checked_mul(width)
        .filter(|&len| Layout::array::<f64>(len).is_ok())
    else {
        return bicgstab_l_report(x, initial(), 0, tol);
    };

    let mut tmp = vec![0.0; n];
    // operator B v = A M^{-1} v; unknown correction y with x = x0 + M^{-1} y
    let op = |v: &[f64], out: &mut [f64], tmp: &mut [f64]| {
        if !bicgstab_l_finite(v) {
            return false;
        }
        tmp.fill(f64::NAN);
        precond(v, tmp);
        if !bicgstab_l_finite(tmp) {
            return false;
        }
        // Every output entry participates in a checked dot product or squared
        // norm before that direction is reused.
        out.fill(f64::NAN);
        apply(tmp, out);
        true
    };
    let mut xc = vec![0.0; n]; // correction in the preconditioned space
    let rt = r_init;
    let mut r = vec![vec![0.0; n]; width];
    let mut u = vec![vec![0.0; n]; width];
    r[0].copy_from_slice(&rt);
    let mut tau = vec![0.0; small_matrix_len];
    let mut sigma = vec![0.0; width];
    let mut gp = vec![0.0; width];
    let mut g = vec![0.0; width];
    let mut gpp = vec![0.0; width];
    let (mut rho0, mut alpha, mut omega) = (1.0_f64, 0.0, 1.0);
    let mut it = 0;
    'cycles: while rn / bn > tol && it < max_iter {
        it += 1;
        rho0 *= -omega;
        for j in 0..l {
            let rho1 = dot(&rt, &r[j]);
            if rho0 == 0.0 || !rho0.is_finite() || !rho1.is_finite() {
                break 'cycles;
            }
            let beta = alpha * rho1 / rho0;
            if !beta.is_finite() {
                break 'cycles;
            }
            rho0 = rho1;
            for i in 0..=j {
                for q in 0..n {
                    u[i][q] = r[i][q] - beta * u[i][q];
                }
            }
            let (u_prefix, u_next) = u.split_at_mut(j + 1);
            if !op(&u_prefix[j], &mut u_next[0], &mut tmp) {
                break 'cycles;
            }
            let gamma = dot(&u[j + 1], &rt);
            if gamma == 0.0 || !gamma.is_finite() {
                break 'cycles;
            }
            alpha = rho0 / gamma;
            if !alpha.is_finite() {
                break 'cycles;
            }
            for i in 0..=j {
                axpy(-alpha, &u[i + 1], &mut r[i]);
            }
            let (r_prefix, r_next) = r.split_at_mut(j + 1);
            if !op(&r_prefix[j], &mut r_next[0], &mut tmp) {
                break 'cycles;
            }
            axpy(alpha, &u[0], &mut xc);
        }
        // Minimal residual part. Every used entry of the upper triangle
        // and the coefficient vectors is overwritten in this cycle.
        for j in 1..=l {
            let (prefix, current) = r.split_at_mut(j);
            let rj = &mut current[0];
            for i in 1..j {
                let coefficient = if sigma[i] == 0.0 {
                    0.0
                } else {
                    dot(rj, &prefix[i]) / sigma[i]
                };
                if !coefficient.is_finite() {
                    break 'cycles;
                }
                tau[i * width + j] = coefficient;
                axpy(-coefficient, &prefix[i], rj);
            }
            sigma[j] = dot(rj, rj);
            if !sigma[j].is_finite() || (sigma[j] == 0.0 && rj.iter().any(|&v| v != 0.0)) {
                break 'cycles;
            }
            gp[j] = if sigma[j] == 0.0 {
                0.0
            } else {
                dot(&prefix[0], rj) / sigma[j]
            };
            if !gp[j].is_finite() {
                break 'cycles;
            }
        }
        g[l] = gp[l];
        omega = g[l];
        for j in (1..l).rev() {
            let mut s = gp[j];
            for i in j + 1..=l {
                s -= tau[j * width + i] * g[i];
            }
            g[j] = s;
        }
        for j in 1..l {
            let mut s = g[j + 1];
            for i in j + 1..l {
                s += tau[j * width + i] * g[i + 1];
            }
            gpp[j] = s;
        }
        if g[1..].iter().chain(&gpp[1..l]).any(|v| !v.is_finite()) {
            break;
        }
        let (r0, r_high) = r.split_at_mut(1);
        let (u0, u_high) = u.split_at_mut(1);
        axpy(g[1], &r0[0], &mut xc);
        axpy(-gp[l], &r_high[l - 1], &mut r0[0]);
        axpy(-g[l], &u_high[l - 1], &mut u0[0]);
        for j in 1..l {
            axpy(-g[j], &u_high[j - 1], &mut u0[0]);
            let rj = &r_high[j - 1];
            axpy(gpp[j], rj, &mut xc);
            axpy(-gp[j], rj, &mut r0[0]);
        }
        rn = nrm(&r[0]);
        if !rn.is_finite() {
            break;
        }
    }
    if bicgstab_l_finite(&xc) {
        tmp.fill(f64::NAN);
        precond(&xc, &mut tmp);
        if x.iter().zip(&tmp).all(|(&a, &c)| (a + c).is_finite()) {
            for (value, &correction) in x.iter_mut().zip(&tmp) {
                *value += correction;
            }
            image.fill(f64::NAN);
            apply(&x, &mut image);
        }
    }
    let residual = residual_ratio(b, &image).unwrap_or(f64::INFINITY);
    bicgstab_l_report(x, residual, it, tol)
}

#[allow(clippy::needless_bitwise_bool)] // An eager reduction can use vector instructions.
fn bicgstab_l_finite(values: &[f64]) -> bool {
    values
        .iter()
        .fold(true, |finite, value| finite & value.is_finite())
}

const fn bicgstab_l_report(
    x: Vec<f64>,
    residual: f64,
    it: usize,
    tol: f64,
) -> IterResult {
    IterResult {
        x,
        residual,
        iterations: it,
        converged: residual.is_finite() && residual <= tol,
    }
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
    let IdrOptions {
        s,
        tol,
        max_iter,
        seed,
    } = *opts;
    let n = b.len();
    let s = s.clamp(1, n.max(1));
    let bn = nrm(b).max(f64::MIN_POSITIVE);
    let mut x = x0.map_or_else(|| vec![0.0; n], <[f64]>::to_vec);
    let mut r = residual(&apply, b, &x);
    let mut mv = 1usize;
    if nrm(&r) / bn <= tol {
        let res = nrm(&r) / bn;
        return IterResult {
            x,
            residual: res,
            iterations: 0,
            converged: true,
        };
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
    IterResult {
        x,
        residual: res,
        iterations: mv,
        converged: res <= 10.0 * tol.max(1e-15),
    }
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
    idr_s_with(
        apply,
        precond,
        b,
        x0,
        &IdrOptions {
            s,
            tol,
            max_iter,
            seed,
        },
    )
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
    IterResult {
        x,
        residual: res,
        iterations: total,
        converged: res <= tol.max(1e-15) * 10.0,
    }
}
