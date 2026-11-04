//! # Optimisation kernels
//!
//! * unconstrained: Nelder-Mead, BFGS, L-BFGS, dogleg trust region;
//! * global, derivative free: differential evolution and simulated
//!   annealing (seeded, deterministic);
//! * constrained: sequential quadratic programming (SQP) with a
//!   least-squares restoration phase, which drives `h(x) = 0` and
//!   `g(x) <= 0` to tolerance instead of discouraging violations with a
//!   penalty; an augmented Lagrangian is also available;
//! * linear programming: two-phase tableau simplex with Dantzig or
//!   steepest-edge pricing and Bland's rule as an anti-cycling fallback.
//!
//! # Limitations
//!
//! * Only smooth problems. [`sqp`] differentiates by central differences
//!   and takes a Newton step on the Lagrangian, so a kink (`abs`, `max`,
//!   a piecewise objective) defeats the line search; use [`nelder_mead`],
//!   [`differential_evolution`] or a reformulation for those.
//! * The multipliers are those of the linear-independence constraint
//!   qualification: a redundant or degenerate constraint at the solution
//!   (two rows with proportional gradients) leaves them non-unique, and
//!   only one of the possible multiplier vectors is reported.
//! * Only first-order (KKT) optimality is what [`sqp`] verifies while it
//!   iterates, but the returned [`ConstrainedResult::second_order`] now tests
//!   second-order sufficiency at the reported KKT point: the dense Hessian of
//!   the Lagrangian (central differences on the objective and the multiplier
//!   combination) is projected onto the null space of the active constraint
//!   Jacobian and its smallest eigenvalue is compared with a finite-difference
//!   noise margin ([`SecondOrderStatus`], also available standalone as
//!   [`second_order_check`]). The verdict is [`SecondOrderStatus::Inconclusive`]
//!   when there is no feasible direction (the active gradients span `R^n`), at
//!   a non-KKT outcome, or when the reduced eigenvalue is inside the margin, so
//!   a saddle is flagged rather than silently called optimal, but sufficiency
//!   is never claimed on a degenerate eigenvalue. The Hessian is dense and
//!   finite-difference, not analytic, and the active set uses the reported
//!   multipliers (degenerate/LICQ-violating constraints can make the projection
//!   ill-conditioned and the verdict inconclusive).
//! * Infeasibility is a *local* verdict: [`ConstrainedStatus::Infeasible`]
//!   means the restoration phase reached a stationary point of the summed
//!   squared violation with a positive residual. It is sound for the small
//!   convex problems the routine is aimed at, not a global infeasibility
//!   proof.
//! * No second-order cone, semidefinite or integer constraints, no sparse
//!   linear algebra: the subproblems are dense and cost `O(n^3)` per
//!   iteration. [`augmented_lagrangian`] is retained for problems where
//!   the SQP subproblem is too large, but it stops on a penalty and its
//!   multipliers are only asymptotically the true ones.
//! * [`simplex`] solves linear programs only; a quadratic objective with
//!   linear constraints must go through [`sqp`].

use crate::kernels::dense::Mat;
use crate::kernels::dense::{
    self,
};
use crate::kernels::random::Rng;

/// Result of a minimisation.
#[derive(Debug, Clone, PartialEq)]
pub struct OptResult {
    /// The best point found.
    pub x: Vec<f64>,
    /// Objective value at `x`.
    pub fx: f64,
    /// Iterations (generations for global methods).
    pub iterations: usize,
    /// Number of objective evaluations.
    pub evaluations: usize,
    /// Whether the stopping tolerance was reached.
    pub converged: bool,
}

/// Central-difference gradient of `f` at `x`.
pub fn numerical_gradient<F: Fn(&[f64]) -> f64>(
    f: &F,
    x: &[f64],
) -> Vec<f64> {
    let mut xp = x.to_vec();
    (0..x.len())
        .map(|i| {
            let h = 6e-6 * x[i].abs().max(1.0);
            xp[i] = x[i] + h;
            let a = f(&xp);
            xp[i] = x[i] - h;
            let b = f(&xp);
            xp[i] = x[i];
            (a - b) / (2.0 * h)
        })
        .collect()
}

fn dot(
    a: &[f64],
    b: &[f64],
) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// Nelder-Mead simplex method. `step` is the initial simplex edge length.
pub fn nelder_mead<F: Fn(&[f64]) -> f64>(
    f: F,
    x0: &[f64],
    step: f64,
    tol: f64,
    max_iter: usize,
) -> OptResult {
    let n = x0.len();
    let mut pts: Vec<Vec<f64>> = vec![x0.to_vec()];
    for i in 0..n {
        let mut p = x0.to_vec();
        p[i] += if p[i] == 0.0 {
            step * 0.05_f64.max(step)
        } else {
            step * p[i].abs().max(1.0) * p[i].signum()
        };
        pts.push(p);
    }
    let mut vals: Vec<f64> = pts.iter().map(|p| f(p)).collect();
    let mut evals = n + 1;
    let mut it = 0;
    let mut converged = false;
    while it < max_iter {
        it += 1;
        let mut idx: Vec<usize> = (0..=n).collect();
        idx.sort_by(|&a, &b| vals[a].total_cmp(&vals[b]));
        pts = idx.iter().map(|&i| pts[i].clone()).collect();
        vals = idx.iter().map(|&i| vals[i]).collect();
        let spread = (vals[n] - vals[0]).abs();
        let size = pts[1..]
            .iter()
            .map(|p| {
                p.iter()
                    .zip(&pts[0])
                    .map(|(a, b)| (a - b).abs())
                    .fold(0.0, f64::max)
            })
            .fold(0.0, f64::max);
        if spread <= tol * (1.0 + vals[0].abs()) && size <= tol.sqrt() {
            converged = true;
            break;
        }
        let cen: Vec<f64> = (0..n)
            .map(|j| pts[..n].iter().map(|p| p[j]).sum::<f64>() / n as f64)
            .collect();
        let along =
            |t: f64| -> Vec<f64> { (0..n).map(|j| cen[j] + t * (pts[n][j] - cen[j])).collect() };
        let xr = along(-1.0);
        let fr = f(&xr);
        evals += 1;
        if fr < vals[0] {
            let xe = along(-2.0);
            let fe = f(&xe);
            evals += 1;
            if fe < fr {
                pts[n] = xe;
                vals[n] = fe;
            } else {
                pts[n] = xr;
                vals[n] = fr;
            }
        } else if fr < vals[n - 1] {
            pts[n] = xr;
            vals[n] = fr;
        } else {
            let (xc, fc) = if fr < vals[n] {
                let x = along(-0.5);
                let v = f(&x);
                (x, v)
            } else {
                let x = along(0.5);
                let v = f(&x);
                (x, v)
            };
            evals += 1;
            if fc < vals[n].min(fr) {
                pts[n] = xc;
                vals[n] = fc;
            } else {
                for i in 1..=n {
                    for j in 0..n {
                        pts[i][j] = pts[0][j] + 0.5 * (pts[i][j] - pts[0][j]);
                    }
                    vals[i] = f(&pts[i]);
                    evals += 1;
                }
            }
        }
    }
    let best = (0..=n)
        .min_by(|&a, &b| vals[a].total_cmp(&vals[b]))
        .unwrap_or(0);
    OptResult {
        x: pts[best].clone(),
        fx: vals[best],
        iterations: it,
        evaluations: evals,
        converged,
    }
}

/// Backtracking line search with the Armijo condition and a curvature
/// check by bisection-expansion (weak Wolfe). Returns `(step, f_new)`.
fn line_search<F: Fn(&[f64]) -> f64, G: Fn(&[f64]) -> Vec<f64>>(
    f: &F,
    g: &G,
    x: &[f64],
    fx: f64,
    gx: &[f64],
    d: &[f64],
    evals: &mut usize,
) -> (f64, f64) {
    let slope = dot(gx, d);
    let (mut lo, mut hi) = (0.0, f64::INFINITY);
    let mut t = 1.0;
    let mut best = (0.0, fx);
    for _ in 0..60 {
        let xn: Vec<f64> = x.iter().zip(d).map(|(a, b)| a + t * b).collect();
        let fnw = f(&xn);
        *evals += 1;
        if !fnw.is_finite() || fnw > fx + 1e-4 * t * slope {
            hi = t;
        } else {
            best = (t, fnw);
            let gn = g(&xn);
            if dot(&gn, d) < 0.9 * slope {
                lo = t;
            } else {
                return best;
            }
        }
        t = if hi.is_finite() {
            f64::midpoint(lo, hi)
        } else {
            2.0 * lo.max(t)
        };
    }
    best
}

/// BFGS quasi-Newton method with a weak-Wolfe line search. If `grad` is
/// `None` a central-difference gradient is used.
pub fn bfgs<F: Fn(&[f64]) -> f64, G: Fn(&[f64]) -> Vec<f64>>(
    f: F,
    grad: Option<&G>,
    x0: &[f64],
    tol: f64,
    max_iter: usize,
) -> OptResult {
    let gfun = |x: &[f64]| grad.map_or_else(|| numerical_gradient(&f, x), |g| g(x));
    let n = x0.len();
    let mut x = x0.to_vec();
    let mut fx = f(&x);
    let mut gx = gfun(&x);
    let mut h = Mat::identity(n);
    let mut fresh = true;
    let mut evals = 1;
    let mut it = 0;
    let mut converged = false;
    while it < max_iter {
        if gx.iter().fold(0.0_f64, |m, v| m.max(v.abs())) <= tol {
            converged = true;
            break;
        }
        it += 1;
        let mut d: Vec<f64> = h.matvec(&gx).iter().map(|v| -v).collect();
        if dot(&d, &gx) >= 0.0 {
            h = Mat::identity(n);
            d = gx.iter().map(|v| -v).collect();
        }
        let (t, fnw) = line_search(&f, &gfun, &x, fx, &gx, &d, &mut evals);
        if t == 0.0 {
            if fresh {
                break;
            }
            h = Mat::identity(n);
            fresh = true;
            continue;
        }
        fresh = false;
        let xn: Vec<f64> = x.iter().zip(&d).map(|(a, b)| a + t * b).collect();
        let gn = gfun(&xn);
        let s: Vec<f64> = xn.iter().zip(&x).map(|(a, b)| a - b).collect();
        let y: Vec<f64> = gn.iter().zip(&gx).map(|(a, b)| a - b).collect();
        let sy = dot(&s, &y);
        if sy > 1e-12 {
            if it == 1 {
                let scale = sy / dot(&y, &y);
                h = Mat::identity(n);
                for i in 0..n {
                    h.set(i, i, scale);
                }
            }
            let rho = 1.0 / sy;
            let hy = h.matvec(&y);
            let yhy = dot(&y, &hy);
            for i in 0..n {
                for j in 0..n {
                    let v = h.at(i, j) - rho * (hy[i] * s[j] + s[i] * hy[j])
                        + (rho * rho * yhy + rho) * s[i] * s[j];
                    h.set(i, j, v);
                }
            }
        }
        x = xn;
        fx = fnw;
        gx = gn;
    }
    OptResult {
        x,
        fx,
        iterations: it,
        evaluations: evals,
        converged,
    }
}

/// Limited-memory BFGS with history `m` (two-loop recursion).
pub fn lbfgs<F: Fn(&[f64]) -> f64, G: Fn(&[f64]) -> Vec<f64>>(
    f: F,
    grad: Option<&G>,
    x0: &[f64],
    m: usize,
    tol: f64,
    max_iter: usize,
) -> OptResult {
    let gfun = |x: &[f64]| grad.map_or_else(|| numerical_gradient(&f, x), |g| g(x));
    let mut x = x0.to_vec();
    let mut fx = f(&x);
    let mut gx = gfun(&x);
    let mut hist: Vec<(Vec<f64>, Vec<f64>, f64)> = Vec::new();
    let mut evals = 1;
    let mut it = 0;
    let mut converged = false;
    while it < max_iter {
        if gx.iter().fold(0.0_f64, |a, v| a.max(v.abs())) <= tol {
            converged = true;
            break;
        }
        it += 1;
        let mut q = gx.clone();
        let mut alpha = vec![0.0; hist.len()];
        for (i, (s, y, rho)) in hist.iter().enumerate().rev() {
            alpha[i] = rho * dot(s, &q);
            for j in 0..q.len() {
                q[j] -= alpha[i] * y[j];
            }
        }
        if let Some((s, y, _)) = hist.last() {
            let gamma = dot(s, y) / dot(y, y);
            for v in &mut q {
                *v *= gamma;
            }
        }
        for (i, (s, y, rho)) in hist.iter().enumerate() {
            let beta = rho * dot(y, &q);
            for j in 0..q.len() {
                q[j] += s[j] * (alpha[i] - beta);
            }
        }
        let mut d: Vec<f64> = q.iter().map(|v| -v).collect();
        if dot(&d, &gx) >= 0.0 {
            hist.clear();
            d = gx.iter().map(|v| -v).collect();
        }
        let (t, fnw) = line_search(&f, &gfun, &x, fx, &gx, &d, &mut evals);
        if t == 0.0 {
            if hist.is_empty() {
                break;
            }
            hist.clear();
            continue;
        }
        let xn: Vec<f64> = x.iter().zip(&d).map(|(a, b)| a + t * b).collect();
        let gn = gfun(&xn);
        let s: Vec<f64> = xn.iter().zip(&x).map(|(a, b)| a - b).collect();
        let y: Vec<f64> = gn.iter().zip(&gx).map(|(a, b)| a - b).collect();
        let sy = dot(&s, &y);
        if sy > 1e-12 {
            if hist.len() == m.max(1) {
                hist.remove(0);
            }
            hist.push((s, y, 1.0 / sy));
        }
        x = xn;
        fx = fnw;
        gx = gn;
    }
    OptResult {
        x,
        fx,
        iterations: it,
        evaluations: evals,
        converged,
    }
}

/// Trust-region method with the dogleg step. The Hessian is supplied
/// by `hess`; pass `None` for a finite-difference Hessian from gradients.
pub fn trust_region<F, G, H>(
    f: F,
    grad: G,
    hess: Option<&H>,
    x0: &[f64],
    tol: f64,
    max_iter: usize,
) -> OptResult
where
    F: Fn(&[f64]) -> f64,
    G: Fn(&[f64]) -> Vec<f64>,
    H: Fn(&[f64]) -> Mat,
{
    let n = x0.len();
    let hfun = |x: &[f64]| -> Mat {
        hess.map_or_else(
            || {
                let g0 = grad(x);
                let mut hm = Mat::zeros(n, n);
                let mut xp = x.to_vec();
                for j in 0..n {
                    let h = 1e-6 * x[j].abs().max(1.0);
                    xp[j] = x[j] + h;
                    let g1 = grad(&xp);
                    xp[j] = x[j];
                    for i in 0..n {
                        hm.set(i, j, (g1[i] - g0[i]) / h);
                    }
                }
                for i in 0..n {
                    for j in 0..i {
                        let v = f64::midpoint(hm.at(i, j), hm.at(j, i));
                        hm.set(i, j, v);
                        hm.set(j, i, v);
                    }
                }
                hm
            },
            |h| h(x),
        )
    };
    let mut x = x0.to_vec();
    let mut fx = f(&x);
    let mut radius = 1.0;
    let mut evals = 1;
    let mut it = 0;
    let mut converged = false;
    while it < max_iter {
        let g = grad(&x);
        if g.iter().fold(0.0_f64, |a, v| a.max(v.abs())) <= tol {
            converged = true;
            break;
        }
        it += 1;
        let b = hfun(&x);
        let gn2 = dot(&g, &g);
        let bg = b.matvec(&g);
        let gbg = dot(&g, &bg);
        // Cauchy point
        let tau = if gbg <= 0.0 {
            f64::INFINITY
        } else {
            gn2 / gbg
        };
        let pu: Vec<f64> = g.iter().map(|v| -tau * v).collect();
        let neg: Vec<f64> = g.iter().map(|v| -v).collect();
        let pb = dense::solve(&b, &neg)
            .ok()
            .filter(|p| dot(p, &b.matvec(p)) > 0.0);
        let norm = |v: &[f64]| dot(v, v).sqrt();
        let step: Vec<f64> = match pb {
            | Some(pb) if norm(&pb) <= radius => pb,
            | Some(pb) if tau.is_finite() && norm(&pu) < radius => {
                // intersection of the segment pu -> pb with the boundary
                let d: Vec<f64> = pb.iter().zip(&pu).map(|(a, c)| a - c).collect();
                let (a2, b2, c2) = (
                    dot(&d, &d),
                    2.0 * dot(&pu, &d),
                    dot(&pu, &pu) - radius * radius,
                );
                let s = (-b2 + (b2 * b2 - 4.0 * a2 * c2).max(0.0).sqrt()) / (2.0 * a2);
                pu.iter().zip(&d).map(|(a, c)| a + s * c).collect()
            },
            | _ => {
                let gn = gn2.sqrt();
                g.iter().map(|v| -radius * v / gn).collect()
            },
        };
        let bp = b.matvec(&step);
        let pred = -(dot(&g, &step) + 0.5 * dot(&step, &bp));
        let xn: Vec<f64> = x.iter().zip(&step).map(|(a, c)| a + c).collect();
        let fnw = f(&xn);
        evals += 1;
        let rho = if pred.abs() < 1e-300 {
            0.0
        } else {
            (fx - fnw) / pred
        };
        if rho < 0.25 {
            radius = 0.25 * norm(&step);
        } else if rho > 0.75 && (norm(&step) - radius).abs() < 1e-10 * radius.max(1.0) {
            radius = (2.0 * radius).min(1e10);
        }
        if rho > 1e-4 && fnw.is_finite() {
            x = xn;
            fx = fnw;
        }
        if radius < 1e-14 {
            break;
        }
    }
    OptResult {
        x,
        fx,
        iterations: it,
        evaluations: evals,
        converged,
    }
}

/// Differential evolution (`rand/1/bin`) within box `bounds`, seeded.
pub fn differential_evolution<F: Fn(&[f64]) -> f64>(
    f: F,
    bounds: &[(f64, f64)],
    pop_size: usize,
    generations: usize,
    seed: u64,
) -> OptResult {
    let n = bounds.len();
    let np = pop_size.max(4);
    let mut rng = Rng::new(seed);
    let mut pop: Vec<Vec<f64>> = (0..np)
        .map(|_| {
            bounds
                .iter()
                .map(|&(a, b)| rng.uniform_range(a, b))
                .collect()
        })
        .collect();
    let mut vals: Vec<f64> = pop.iter().map(|p| f(p)).collect();
    let mut evals = np;
    let (fw, cr) = (0.7, 0.9);
    let mut generation = 0;
    while generation < generations {
        generation += 1;
        for i in 0..np {
            let mut pick = || {
                loop {
                    let k = rng.below(np as u64) as usize;
                    if k != i {
                        return k;
                    }
                }
            };
            let (a, b, c) = loop {
                let (a, b, c) = (pick(), pick(), pick());
                if a != b && b != c && a != c {
                    break (a, b, c);
                }
            };
            let jr = rng.below(n as u64) as usize;
            let trial: Vec<f64> = (0..n)
                .map(|j| {
                    if rng.uniform() < cr || j == jr {
                        (pop[a][j] + fw * (pop[b][j] - pop[c][j])).clamp(bounds[j].0, bounds[j].1)
                    } else {
                        pop[i][j]
                    }
                })
                .collect();
            let ft = f(&trial);
            evals += 1;
            if ft <= vals[i] {
                pop[i] = trial;
                vals[i] = ft;
            }
        }
        let mean = vals.iter().sum::<f64>() / np as f64;
        let sd = (vals.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / np as f64).sqrt();
        if sd <= 1e-12 * (1.0 + mean.abs()) {
            break;
        }
    }
    let best = (0..np)
        .min_by(|&a, &b| vals[a].total_cmp(&vals[b]))
        .unwrap_or(0);
    OptResult {
        x: pop[best].clone(),
        fx: vals[best],
        iterations: generation,
        evaluations: evals,
        converged: generation < generations,
    }
}

/// Simulated annealing with geometric cooling and Gaussian proposals of
/// scale `step * T / T0`-adjusted; seeded and deterministic.
pub fn simulated_annealing<F: Fn(&[f64]) -> f64>(
    f: F,
    x0: &[f64],
    step: f64,
    t0: f64,
    cooling: f64,
    iterations: usize,
    seed: u64,
) -> OptResult {
    let mut rng = Rng::new(seed);
    let mut x = x0.to_vec();
    let mut fx = f(&x);
    let mut best = (x.clone(), fx);
    let mut t = t0;
    for _ in 0..iterations {
        let scale = step * (t / t0).sqrt().max(1e-3);
        let cand: Vec<f64> = x.iter().map(|v| v + scale * rng.normal()).collect();
        let fc = f(&cand);
        if fc < fx || rng.uniform() < ((fx - fc) / t).exp() {
            x = cand;
            fx = fc;
            if fx < best.1 {
                best = (x.clone(), fx);
            }
        }
        t *= cooling;
    }
    OptResult {
        x: best.0,
        fx: best.1,
        iterations,
        evaluations: iterations + 1,
        converged: true,
    }
}

/// Augmented Lagrangian for `min f(x)` subject to `h(x) = 0` and
/// `g(x) <= 0`; the inner problems are solved by BFGS.
pub fn augmented_lagrangian<F, H, G>(
    f: F,
    eq: H,
    ineq: G,
    x0: &[f64],
    tol: f64,
    max_outer: usize,
) -> OptResult
where
    F: Fn(&[f64]) -> f64,
    H: Fn(&[f64]) -> Vec<f64>,
    G: Fn(&[f64]) -> Vec<f64>,
{
    let mut x = x0.to_vec();
    let ne = eq(&x).len();
    let ni = ineq(&x).len();
    let mut lam = vec![0.0; ne];
    let mut mu = vec![0.0; ni];
    let mut rho = 10.0;
    let mut evals = 0;
    let mut outer = 0;
    let mut converged = false;
    for _ in 0..max_outer {
        outer += 1;
        let (lc, mc) = (lam.clone(), mu.clone());
        let phi = |z: &[f64]| -> f64 {
            let h = eq(z);
            let g = ineq(z);
            let mut v = f(z);
            for i in 0..ne {
                v += lc[i] * h[i] + 0.5 * rho * h[i] * h[i];
            }
            for i in 0..ni {
                let t = (mc[i] + rho * g[i]).max(0.0);
                v += (t * t - mc[i] * mc[i]) / (2.0 * rho);
            }
            v
        };
        let r = bfgs(phi, None::<&fn(&[f64]) -> Vec<f64>>, &x, 1e-9, 500);
        evals += r.evaluations;
        x = r.x;
        let h = eq(&x);
        let g = ineq(&x);
        for i in 0..ne {
            lam[i] += rho * h[i];
        }
        for i in 0..ni {
            mu[i] = (mu[i] + rho * g[i]).max(0.0);
        }
        let viol = h
            .iter()
            .map(|v| v.abs())
            .chain(g.iter().map(|v| v.max(0.0)))
            .fold(0.0, f64::max);
        if viol <= tol {
            let comp = g
                .iter()
                .zip(&mu)
                .map(|(a, b)| (a * b).abs())
                .fold(0.0, f64::max);
            if comp <= tol.max(1e-6) {
                converged = true;
                break;
            }
        }
        rho = (rho * 4.0).min(1e8);
    }
    let fx = f(&x);
    OptResult {
        x,
        fx,
        iterations: outer,
        evaluations: evals,
        converged,
    }
}

/// Outcome of a constrained minimisation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConstrainedStatus {
    /// A KKT point was reached: the constraint violation is at most the
    /// requested tolerance and the stationarity residual is at most
    /// `max(tolerance, 1e-8)`. This is a *first-order* verdict; see
    /// [`ConstrainedResult::second_order`] for the second-order one.
    Optimal,
    /// The constraints could not be satisfied to tolerance. The
    /// least-squares restoration phase stopped at a point where no
    /// direction lowers the summed squared violation, so no nearby
    /// feasible point exists.
    Infeasible,
    /// The iteration limit was reached before either verdict.
    IterationLimit,
}

/// Second-order sufficiency verdict at the KKT point reported by [`sqp`].
///
/// The check builds the Hessian of the Lagrangian `L = f + lambda . h + mu . g`,
/// projects it onto the null space of the active constraint Jacobian (all
/// equality rows, plus the inequality rows that are active at the point), and
/// asks whether the reduced Hessian is positive definite with a margin.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SecondOrderStatus {
    /// The reduced Hessian is positive definite with a margin: the KKT point is
    /// a strict local minimiser to the accuracy of the finite-difference check.
    Sufficient,
    /// A feasible direction has negative curvature: the KKT point is a saddle
    /// (or a maximum) along the feasible set, not a local minimiser.
    NotSufficient,
    /// Not decided: the outcome was not [`ConstrainedStatus::Optimal`], there
    /// is no feasible direction to test, or the smallest reduced eigenvalue is
    /// inside the finite-difference noise margin. The check never guesses.
    Inconclusive,
}

/// Result of [`sqp`], with feasibility and optimality reported separately.
#[derive(Debug, Clone, PartialEq)]
pub struct ConstrainedResult {
    /// Final point.
    pub x: Vec<f64>,
    /// Objective value at `x`.
    pub fx: f64,
    /// Outcome of the solve.
    pub status: ConstrainedStatus,
    /// Number of SQP iterations taken.
    pub iterations: usize,
    /// Number of objective, constraint and difference evaluations.
    pub evaluations: usize,
    /// Feasibility measure `max(max_i |h_i|, max_i max(0, g_i))` at `x`.
    pub violation: f64,
    /// Optimality measure at `x`: the infinity norm of
    /// `grad f + sum_i lambda_i grad h_i + sum_i mu_i grad g_i`, i.e. of
    /// the gradient of `L = f + lambda . h + mu . g`. This is the KKT
    /// stationarity residual; it is `f64::INFINITY` when no feasible point
    /// was found.
    pub kkt_residual: f64,
    /// Multipliers of `h`, in the order the constraint closure returns.
    /// The sign convention is that of `L = f + lambda . h + mu . g`.
    pub lambda: Vec<f64>,
    /// Multipliers of `g`: non-negative, and zero for every constraint
    /// that is not active at `x` (complementary slackness).
    pub mu: Vec<f64>,
    /// Second-order sufficiency verdict at the reported point. Computed when
    /// the status is [`ConstrainedStatus::Optimal`] (with the reported
    /// multipliers) and [`SecondOrderStatus::Inconclusive`] otherwise; the
    /// evaluation cost of the check is included in
    /// [`ConstrainedResult::evaluations`].
    pub second_order: SecondOrderStatus,
}

/// Solution of the quadratic program solved at each SQP iteration.
#[derive(Debug, Clone, PartialEq)]
pub struct QpSolution {
    /// Step `d` minimising `0.5 d^T H d + c^T d`.
    pub d: Vec<f64>,
    /// Multipliers of the equality rows `A_eq d = b_eq`.
    pub lambda: Vec<f64>,
    /// Multipliers of the inequality rows `A_ineq d <= b_ineq`, zero for
    /// the rows that are not in the working set.
    pub mu: Vec<f64>,
}

/// `v.get(i)`, defaulting to zero, so that no indexing is used.
fn pick(
    v: &[f64],
    i: usize,
) -> f64 {
    v.get(i).copied().unwrap_or(0.0)
}

/// Largest absolute equality violation or positive inequality violation.
fn max_violation(
    h: &[f64],
    g: &[f64],
) -> f64 {
    h.iter()
        .map(|v| v.abs())
        .chain(g.iter().map(|v| v.max(0.0)))
        .fold(0.0, f64::max)
}

/// Central-difference Jacobian of a vector-valued constraint.
fn fd_jacobian<C: Fn(&[f64]) -> Vec<f64>>(
    c: &C,
    x: &[f64],
    rows: usize,
    evals: &mut usize,
) -> Vec<Vec<f64>> {
    let mut jac = vec![vec![0.0; x.len()]; rows];
    let mut xp = x.to_vec();
    for (j, &xj) in x.iter().enumerate() {
        let step = 6e-6 * xj.abs().max(1.0);
        if let Some(slot) = xp.get_mut(j) {
            *slot = xj + step;
        }
        let a = c(&xp);
        if let Some(slot) = xp.get_mut(j) {
            *slot = xj - step;
        }
        let b = c(&xp);
        if let Some(slot) = xp.get_mut(j) {
            *slot = xj;
        }
        *evals += 2;
        for (i, row) in jac.iter_mut().enumerate() {
            if let Some(slot) = row.get_mut(j) {
                *slot = (pick(&a, i) - pick(&b, i)) / (2.0 * step);
            }
        }
    }
    jac
}

/// Solves `min 0.5 d^T H d + c^T d` subject to `A_eq d = b_eq` and
/// `A_ineq d <= b_ineq` by a primal active-set method started from `d = 0`,
/// which must be feasible (`b_ineq >= 0`).
///
/// The KKT system of the working set is solved exactly at every step, so
/// the returned multipliers are those of the Lagrangian
/// `0.5 d^T H d + c^T d + lambda . (A_eq d - b_eq) + mu . (A_ineq d - b_ineq)`.
/// Returns `None` if a KKT system stays singular after regularisation or
/// the working set fails to settle within the iteration budget.
#[must_use]
pub fn qp_active_set(
    hess: &Mat,
    grad: &[f64],
    a_eq: &[Vec<f64>],
    b_eq: &[f64],
    a_ineq: &[Vec<f64>],
    b_ineq: &[f64],
) -> Option<QpSolution> {
    let n = grad.len();
    let ne = a_eq.len();
    let ni = a_ineq.len();
    if hess.rows != n || hess.cols != n || b_eq.len() != ne || b_ineq.len() != ni {
        return None;
    }
    let mut rows: Vec<&[f64]> = Vec::with_capacity(ne + ni);
    rows.extend(a_eq.iter().map(Vec::as_slice));
    rows.extend(a_ineq.iter().map(Vec::as_slice));
    if rows.iter().any(|r| r.len() != n) {
        return None;
    }
    let mut rhs: Vec<f64> = b_eq.to_vec();
    rhs.extend_from_slice(b_ineq);
    let mut ws: Vec<usize> = (0..ne).collect();
    for (j, b) in b_ineq.iter().enumerate() {
        if *b <= 1e-12 {
            ws.push(ne + j);
        }
    }
    let budget = 50 * (n + ne + ni + 1);
    for _ in 0..budget {
        // KKT matrix of the working set: [H A^T; A 0] [d; mult] = [-c; -b].
        let k = n + ws.len();
        let mut mat = Mat::zeros(k, k);
        for i in 0..n {
            for j in 0..n {
                mat.set(i, j, hess.at(i, j));
            }
        }
        for (col, &w) in ws.iter().enumerate() {
            let row = rows.get(w)?;
            for i in 0..n {
                let v = pick(row, i);
                mat.set(i, n + col, v);
                mat.set(n + col, i, v);
            }
        }
        let mut b = vec![0.0; k];
        for i in 0..n {
            if let Some(slot) = b.get_mut(i) {
                *slot = -pick(grad, i);
            }
        }
        for (col, &w) in ws.iter().enumerate() {
            if let Some(slot) = b.get_mut(n + col) {
                *slot = pick(&rhs, w);
            }
        }
        let solved = dense::solve(&mat, &b).or_else(|_| {
            let mut fixed = mat.clone();
            for i in 0..n {
                let v = fixed.at(i, i) + 1e-10 * (1.0 + fixed.at(i, i).abs());
                fixed.set(i, i, v);
            }
            dense::solve(&fixed, &b)
        });
        let sol = solved.ok()?;
        let d: Vec<f64> = sol.iter().take(n).copied().collect();
        let mult: Vec<f64> = sol.iter().skip(n).copied().collect();
        let dnorm = d.iter().fold(0.0_f64, |a, v| a.max(v.abs()));
        // Most negative inequality multiplier in the working set.
        let worst = ws
            .iter()
            .enumerate()
            .filter(|&(_, &w)| w >= ne)
            .map(|(col, _)| col)
            .min_by(|&a, &b| pick(&mult, a).total_cmp(&pick(&mult, b)));
        if dnorm <= 1e-12 {
            match worst {
                | Some(col) if pick(&mult, col) < -1e-10 => {
                    ws.remove(col);
                    continue;
                },
                | _ => {
                    let mut lambda = vec![0.0; ne];
                    let mut mu = vec![0.0; ni];
                    for (col, &w) in ws.iter().enumerate() {
                        if w < ne {
                            if let Some(slot) = lambda.get_mut(w) {
                                *slot = pick(&mult, col);
                            }
                        } else if let Some(slot) = mu.get_mut(w - ne) {
                            *slot = pick(&mult, col);
                        }
                    }
                    return Some(QpSolution { d, lambda, mu });
                },
            }
        }
        // Longest feasible step along `d`, limited by the inactive rows.
        let mut alpha = 1.0_f64;
        let mut blocker: Option<usize> = None;
        for (j, row) in a_ineq.iter().enumerate() {
            if ws.contains(&(ne + j)) {
                continue;
            }
            let slope = row.iter().zip(&d).map(|(a, v)| a * v).sum::<f64>();
            if slope > 1e-12 {
                let limit = pick(b_ineq, j).max(0.0) / slope;
                if limit < alpha - 1e-12 {
                    alpha = limit.max(0.0);
                    blocker = Some(ne + j);
                }
            }
        }
        if let Some(w) = blocker {
            // Add the blocking row and re-solve; the next KKT system yields
            // both the shortened step and its multipliers.
            if !ws.contains(&w) {
                ws.push(w);
            }
            continue;
        }
        if worst.is_some_and(|col| pick(&mult, col) < -1e-10) {
            if let Some(col) = worst {
                ws.remove(col);
            }
            continue;
        }
        let mut lambda = vec![0.0; ne];
        let mut mu = vec![0.0; ni];
        for (col, &w) in ws.iter().enumerate() {
            if w < ne {
                if let Some(slot) = lambda.get_mut(w) {
                    *slot = pick(&mult, col);
                }
            } else if let Some(slot) = mu.get_mut(w - ne) {
                *slot = pick(&mult, col);
            }
        }
        return Some(QpSolution { d, lambda, mu });
    }
    None
}

/// Sequential quadratic programming for `min f(x)` subject to
/// `h(x) = 0` and `g(x) <= 0`.
///
/// Each iteration linearises the constraints, solves the quadratic program
/// of [`qp_active_set`] for a step `d`, and takes a backtracking step on the
/// `l1` merit function `f + rho (|h|_1 + max(0, g)_1)`. A least-squares
/// restoration phase is run whenever the iterate is infeasible, so the
/// constraints are satisfied to `tol` rather than merely penalised. The
/// multipliers returned in [`ConstrainedResult::lambda`] and
/// [`ConstrainedResult::mu`] are the ones of the final quadratic program,
/// i.e. the exact multipliers of the Lagrangian `L = f + lambda . h + mu . g`;
/// in particular `-mu_i` is the derivative of the optimal value with respect
/// to relaxing `g_i(x) <= 0` to `g_i(x) <= eps`.
///
/// `grad` may be `None`, in which case central differences are used for the
/// objective; constraint Jacobians are always central differences.
///
/// Convergence is declared when the violation is at most `tol` and the KKT
/// residual is at most `max(tol, 1e-8)`.
#[must_use]
pub fn sqp<F, G, H, K>(
    f: F,
    grad: Option<&G>,
    eq: H,
    ineq: K,
    x0: &[f64],
    tol: f64,
    max_iter: usize,
) -> ConstrainedResult
where
    F: Fn(&[f64]) -> f64,
    G: Fn(&[f64]) -> Vec<f64>,
    H: Fn(&[f64]) -> Vec<f64>,
    K: Fn(&[f64]) -> Vec<f64>,
{
    let n = x0.len();
    let gfun = |x: &[f64]| grad.map_or_else(|| numerical_gradient(&f, x), |g| g(x));
    let stat_tol = tol.max(1e-8);
    let mut x = x0.to_vec();
    let mut evals = 0;
    let mut hx = eq(&x);
    let mut gx = ineq(&x);
    evals += 2;
    let ne = hx.len();
    let ni = gx.len();
    let mut lambda = vec![0.0; ne];
    let mut mu = vec![0.0; ni];
    let mut hess = Mat::identity(n);
    let mut fx = f(&x);
    evals += 1;
    let mut iterations = 0;
    let mut status = ConstrainedStatus::IterationLimit;
    let mut second_order = SecondOrderStatus::Inconclusive;
    let mut kkt = f64::INFINITY;
    loop {
        let violation = max_violation(&hx, &gx);
        if violation > tol {
            if restore_feasibility(&mut x, &eq, &ineq, tol, 40, &mut evals) {
                hx = eq(&x);
                gx = ineq(&x);
                fx = f(&x);
                evals += 3;
            } else {
                status = ConstrainedStatus::Infeasible;
                kkt = f64::INFINITY;
                break;
            }
        }
        let gf = gfun(&x);
        let jh = fd_jacobian(&eq, &x, ne, &mut evals);
        let jg = fd_jacobian(&ineq, &x, ni, &mut evals);
        let beq: Vec<f64> = hx.iter().map(|v| -v).collect();
        let bineq: Vec<f64> = gx.iter().map(|v| -v).collect();
        let Some(qp) = qp_active_set(&hess, &gf, &jh, &beq, &jg, &bineq) else {
            break;
        };
        kkt = kkt_residual(&gf, &jh, &jg, &qp.lambda, &qp.mu);
        lambda.clone_from(&qp.lambda);
        mu.clone_from(&qp.mu);
        let dnorm = qp.d.iter().fold(0.0_f64, |a, v| a.max(v.abs()));
        let viol = max_violation(&hx, &gx);
        if dnorm <= 1e-10 && viol <= tol && kkt <= stat_tol {
            status = ConstrainedStatus::Optimal;
            let (verdict, check_evals) =
                second_order_status(&f, &eq, &ineq, &x, &lambda, &mu, stat_tol);
            second_order = verdict;
            evals = evals.saturating_add(check_evals);
            break;
        }
        if iterations >= max_iter {
            break;
        }
        iterations += 1;
        // l1 merit function and its linear model along the SQP step.
        let rho = 1.1
            * qp.lambda
                .iter()
                .chain(qp.mu.iter())
                .fold(1.0_f64, |a, v| a.max(v.abs()));
        let phi0 = fx + rho * viol;
        let slope = dot(&gf, &qp.d) - rho * viol;
        let mut alpha = 1.0;
        let mut accepted = false;
        let mut trial = x.clone();
        let mut tfx = fx;
        let mut thx = hx.clone();
        let mut tgx = gx.clone();
        for _ in 0..40 {
            trial = x.iter().zip(&qp.d).map(|(a, b)| a + alpha * b).collect();
            tfx = f(&trial);
            thx = eq(&trial);
            tgx = ineq(&trial);
            evals += 3;
            let phi = tfx + rho * max_violation(&thx, &tgx);
            if phi.is_finite() && phi <= phi0 + 1e-4 * alpha * slope {
                accepted = true;
                break;
            }
            alpha *= 0.5;
        }
        if !accepted {
            break;
        }
        // Damped BFGS update of the Hessian of the Lagrangian.
        let gl_old = lagrangian_gradient(&gf, &jh, &jg, &lambda, &mu);
        let gn = gfun(&trial);
        evals += 1;
        let jhn = fd_jacobian(&eq, &trial, ne, &mut evals);
        let jgn = fd_jacobian(&ineq, &trial, ni, &mut evals);
        let gl_new = lagrangian_gradient(&gn, &jhn, &jgn, &lambda, &mu);
        let s: Vec<f64> = trial.iter().zip(&x).map(|(a, b)| a - b).collect();
        let y: Vec<f64> = gl_new.iter().zip(&gl_old).map(|(a, b)| a - b).collect();
        let sy = dot(&s, &y);
        let bs = hess.matvec(&s);
        let sbs = dot(&s, &bs);
        if sbs > 1e-16 {
            let yy: Vec<f64> = if sy < 0.2 * sbs {
                let theta = 0.8 * sbs / (sbs - sy);
                y.iter()
                    .zip(&bs)
                    .map(|(a, b)| theta * a + (1.0 - theta) * b)
                    .collect()
            } else {
                y
            };
            let syy = dot(&s, &yy);
            if syy > 1e-16 {
                let rho2 = 1.0 / syy;
                let hy = hess.matvec(&yy);
                let yhy = dot(&yy, &hy);
                for i in 0..n {
                    for j in 0..n {
                        let v = hess.at(i, j)
                            - rho2 * (pick(&hy, i) * pick(&s, j) + pick(&s, i) * pick(&hy, j))
                            + (rho2 * rho2 * yhy + rho2) * pick(&s, i) * pick(&s, j);
                        hess.set(i, j, v);
                    }
                }
            }
        }
        x = trial;
        fx = tfx;
        hx = thx;
        gx = tgx;
    }
    ConstrainedResult {
        violation: max_violation(&hx, &gx),
        x,
        fx,
        status,
        iterations,
        evaluations: evals,
        kkt_residual: kkt,
        lambda,
        mu,
        second_order,
    }
}

/// `grad f + J_h^T lambda + J_g^T mu`, the gradient of the Lagrangian.
fn lagrangian_gradient(
    gf: &[f64],
    jh: &[Vec<f64>],
    jg: &[Vec<f64>],
    lambda: &[f64],
    mu: &[f64],
) -> Vec<f64> {
    let mut r = gf.to_vec();
    for (row, l) in jh.iter().zip(lambda) {
        for (a, b) in r.iter_mut().zip(row) {
            *a += l * b;
        }
    }
    for (row, m) in jg.iter().zip(mu) {
        for (a, b) in r.iter_mut().zip(row) {
            *a += m * b;
        }
    }
    r
}

/// Infinity norm of the gradient of the Lagrangian, i.e. the KKT
/// stationarity residual for the multipliers `lambda` (on `h`) and `mu`
/// (on `g`, with `L = f + lambda . h + mu . g`).
#[must_use]
pub fn kkt_residual(
    gf: &[f64],
    jh: &[Vec<f64>],
    jg: &[Vec<f64>],
    lambda: &[f64],
    mu: &[f64],
) -> f64 {
    lagrangian_gradient(gf, jh, jg, lambda, mu)
        .iter()
        .fold(0.0, |a, v| a.max(v.abs()))
}

/// Relative margin inside which the smallest reduced eigenvalue is treated as
/// numerically zero by [`second_order_check`].
const SECOND_ORDER_MARGIN: f64 = 1e-6;

/// Second-order sufficiency check at a KKT point `x`.
///
/// Builds the dense Hessian of the Lagrangian
/// `L = f + sum_i lambda_i h_i + sum_j mu_j g_j` by central differences of `L`
/// itself (so no analytic gradient is required), projects it onto the null
/// space of the active constraint Jacobian and tests the reduced Hessian for
/// positive definiteness. Active inequalities are those with `g_j(x)` at least
/// `-1e-8` or with a strictly positive reported multiplier.
///
/// Returns [`SecondOrderStatus::Sufficient`] when the smallest eigenvalue of the
/// reduced Hessian is larger than a finite-difference margin, and
/// [`SecondOrderStatus::NotSufficient`] when it is more negative than that
/// margin. Returns [`SecondOrderStatus::Inconclusive`] when the active
/// constraint gradients span the whole space (no feasible direction), when the
/// null space or the eigenproblem cannot be formed, or when the eigenvalue is
/// inside the margin.
///
/// The cost is `O(n^2)` evaluations of `f`, `h` and `g` for the Hessian (each
/// call is counted as one evaluation), plus the finite-difference Jacobians of
/// the constraints, and dense `O(n^3)` linear algebra.
#[must_use]
pub fn second_order_check<F, H, K>(
    f: F,
    eq: H,
    ineq: K,
    x: &[f64],
    lambda: &[f64],
    mu: &[f64],
) -> SecondOrderStatus
where
    F: Fn(&[f64]) -> f64,
    H: Fn(&[f64]) -> Vec<f64>,
    K: Fn(&[f64]) -> Vec<f64>,
{
    second_order_status(&f, &eq, &ineq, x, lambda, mu, 1e-8).0
}

/// Implementation of [`second_order_check`] that also reports how many
/// `f`/`h`/`g` evaluations the check used, so [`sqp`] can include them in its
/// count.
fn second_order_status<F, H, K>(
    f: &F,
    eq: &H,
    ineq: &K,
    x: &[f64],
    lambda: &[f64],
    mu: &[f64],
    activation_tol: f64,
) -> (SecondOrderStatus, usize)
where
    F: Fn(&[f64]) -> f64,
    H: Fn(&[f64]) -> Vec<f64>,
    K: Fn(&[f64]) -> Vec<f64>,
{
    let n = x.len();
    if n == 0 {
        return (SecondOrderStatus::Inconclusive, 0);
    }
    let hx = eq(x);
    let gx = ineq(x);
    let ne = hx.len();
    let ni = gx.len();
    let lag_cost = 1 + ne + ni;
    let mut evals = 2;
    let lag = |z: &[f64]| -> f64 {
        let mut value = f(z);
        let h = eq(z);
        let g = ineq(z);
        for (i, l) in lambda.iter().enumerate() {
            value += l * pick(&h, i);
        }
        for (i, m) in mu.iter().enumerate() {
            value += m * pick(&g, i);
        }
        value
    };
    // Hessian of L by central differences. The step is larger than the one
    // used for gradients because a second difference divides by h^2: 1e-4
    // balances truncation (O(h^2)) against round-off (O(eps/h^2)).
    let mut hess = Mat::zeros(n, n);
    let l0 = lag(x);
    evals += lag_cost;
    let mut z = x.to_vec();
    let step = |i: usize| 1e-4 * x[i].abs().max(1.0);
    for i in 0..n {
        let hi = step(i);
        z[i] = x[i] + hi;
        let fp = lag(&z);
        z[i] = x[i] - hi;
        let fm = lag(&z);
        z[i] = x[i];
        hess.set(i, i, (fp - 2.0 * l0 + fm) / (hi * hi));
        for j in 0..i {
            let hj = step(j);
            z[i] = x[i] + hi;
            z[j] = x[j] + hj;
            let pp = lag(&z);
            z[j] = x[j] - hj;
            let pm = lag(&z);
            z[i] = x[i] - hi;
            z[j] = x[j] + hj;
            let mp = lag(&z);
            z[j] = x[j] - hj;
            let mm = lag(&z);
            z[i] = x[i];
            z[j] = x[j];
            let value = (pp - pm - mp + mm) / (4.0 * hi * hj);
            hess.set(i, j, value);
            hess.set(j, i, value);
        }
    }
    // 1 + 2 n^2 Lagrangian evaluations (the diagonal needs 2n, each of the
    // n(n-1)/2 off-diagonal pairs needs 4, and l0 is one).
    evals += lag_cost * (1 + 2 * n * n);
    // Active rows: every equality plus the inequalities active at x.
    let mut rows: Vec<Vec<f64>> = fd_jacobian(eq, x, ne, &mut evals);
    for (j, row) in fd_jacobian(ineq, x, ni, &mut evals).into_iter().enumerate() {
        let g_j = pick(&gx, j);
        let mu_j = pick(mu, j);
        if g_j >= -activation_tol || mu_j > 0.0 {
            rows.push(row);
        }
    }
    // Orthonormal basis of the active row space, two Gram-Schmidt passes.
    let mut basis: Vec<Vec<f64>> = Vec::new();
    for mut row in rows {
        let raw = dot(&row, &row).sqrt();
        for _ in 0..2 {
            for q in &basis {
                let d = dot(&row, q);
                for (a, b) in row.iter_mut().zip(q) {
                    *a -= d * b;
                }
            }
        }
        let norm = dot(&row, &row).sqrt();
        if norm > 1e-10 * raw.max(1.0) {
            for a in &mut row {
                *a /= norm;
            }
            basis.push(row);
        }
    }
    let rank = basis.len();
    let target = n.saturating_sub(rank);
    if target == 0 {
        // The active constraints pin the point: there is no feasible direction.
        return (SecondOrderStatus::Inconclusive, evals);
    }
    // Orthonormal basis of the null space of the active rows.
    let mut null_basis: Vec<Vec<f64>> = Vec::with_capacity(target);
    for k in 0..n {
        let mut v = vec![0.0; n];
        v[k] = 1.0;
        for _ in 0..2 {
            for q in &basis {
                let d = dot(&v, q);
                for (a, b) in v.iter_mut().zip(q) {
                    *a -= d * b;
                }
            }
        }
        let norm = dot(&v, &v).sqrt();
        if norm > 1e-8 {
            for a in &mut v {
                *a /= norm;
            }
            null_basis.push(v);
            if null_basis.len() == target {
                break;
            }
        }
    }
    if null_basis.len() != target {
        return (SecondOrderStatus::Inconclusive, evals);
    }
    // Z^T H Z, symmetrised against round-off before the eigen-decomposition.
    let mut reduced = Mat::zeros(target, target);
    for a in 0..target {
        let hz = hess.matvec(&null_basis[a]);
        for b in 0..target {
            reduced.set(a, b, dot(&hz, &null_basis[b]));
        }
    }
    for a in 0..target {
        for b in 0..a {
            let value = reduced.at(a, b).midpoint(reduced.at(b, a));
            reduced.set(a, b, value);
            reduced.set(b, a, value);
        }
    }
    let Ok((eigenvalues, _)) = dense::eigen_symmetric(&reduced) else {
        return (SecondOrderStatus::Inconclusive, evals);
    };
    let Some(&lambda_min) = eigenvalues.first() else {
        return (SecondOrderStatus::Inconclusive, evals);
    };
    let scale = hess.data.iter().fold(0.0_f64, |acc, v| acc.max(v.abs()));
    let margin = SECOND_ORDER_MARGIN * scale.max(1e-12);
    let verdict = if lambda_min > margin {
        SecondOrderStatus::Sufficient
    } else if lambda_min < -margin {
        SecondOrderStatus::NotSufficient
    } else {
        SecondOrderStatus::Inconclusive
    };
    (verdict, evals)
}

/// Gauss-Newton restoration: moves `x` until the violation is at most
/// `tol`, and reports whether it succeeded. Each step solves
/// `min |A d + c|_2` over the equalities and the violated inequalities and
/// backtracks on the summed squared violation; a step that cannot lower it
/// means that sum is stationary at a positive value.
fn restore_feasibility<H, K>(
    x: &mut Vec<f64>,
    eq: &H,
    ineq: &K,
    tol: f64,
    max_iter: usize,
    evals: &mut usize,
) -> bool
where
    H: Fn(&[f64]) -> Vec<f64>,
    K: Fn(&[f64]) -> Vec<f64>,
{
    let n = x.len();
    let squared = |h: &[f64], g: &[f64]| -> f64 {
        h.iter().map(|v| v * v).sum::<f64>()
            + g.iter().map(|v| v.max(0.0) * v.max(0.0)).sum::<f64>()
    };
    for _ in 0..max_iter {
        let h = eq(x);
        let g = ineq(x);
        *evals += 2;
        if max_violation(&h, &g) <= tol {
            return true;
        }
        let m0 = squared(&h, &g);
        let jh = fd_jacobian(eq, x, h.len(), evals);
        let jg = fd_jacobian(ineq, x, g.len(), evals);
        // Every equality stays in the system so that restoring an
        // inequality does not slide off an equality manifold.
        let mut rows: Vec<&Vec<f64>> = Vec::new();
        let mut res: Vec<f64> = Vec::new();
        for (i, hi) in h.iter().enumerate() {
            if let Some(r) = jh.get(i) {
                rows.push(r);
                res.push(*hi);
            }
        }
        for (i, gi) in g.iter().enumerate() {
            if *gi > tol
                && let Some(r) = jg.get(i)
            {
                rows.push(r);
                res.push(*gi);
            }
        }
        if rows.is_empty() {
            return false;
        }
        // Normal equations of the least-squares step, with a ridge so that
        // an underdetermined or rank-deficient set is still solvable.
        let mut ata = Mat::zeros(n, n);
        let mut atr = vec![0.0; n];
        for (row, r) in rows.iter().zip(&res) {
            for i in 0..n {
                let ri = pick(row, i);
                if let Some(slot) = atr.get_mut(i) {
                    *slot += ri * r;
                }
                for j in 0..n {
                    ata.set(i, j, ata.at(i, j) + ri * pick(row, j));
                }
            }
        }
        let scale = (0..n).map(|i| ata.at(i, i)).fold(0.0_f64, f64::max);
        for i in 0..n {
            ata.set(i, i, ata.at(i, i) + 1e-12 * scale.max(1.0));
        }
        let rhs: Vec<f64> = atr.iter().map(|v| -v).collect();
        let Ok(d) = dense::solve(&ata, &rhs) else {
            return false;
        };
        let mut alpha = 1.0;
        let mut moved = false;
        for _ in 0..30 {
            let xt: Vec<f64> = x.iter().zip(&d).map(|(a, b)| a + alpha * b).collect();
            let mt = squared(&eq(&xt), &ineq(&xt));
            *evals += 2;
            if mt < m0 {
                *x = xt;
                moved = true;
                break;
            }
            alpha *= 0.5;
        }
        if !moved {
            return false;
        }
    }
    max_violation(&eq(x), &ineq(x)) <= tol
}

/// Constraint sense of a linear-programming row.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Relation {
    /// `a . x <= b`
    Le,
    /// `a . x >= b`
    Ge,
    /// `a . x == b`
    Eq,
}

/// Outcome of a linear program.
#[derive(Debug, Clone, PartialEq)]
pub enum LpStatus {
    /// Optimal solution found.
    Optimal {
        /// Optimal point.
        x: Vec<f64>,
        /// Optimal objective value.
        value: f64,
    },
    /// No feasible point exists.
    Infeasible,
    /// The objective is unbounded below.
    Unbounded,
}

/// Entering-variable pricing rule of the simplex method.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Pricing {
    /// Most negative reduced cost (Dantzig's rule).
    Dantzig,
    /// Largest `d_j^2 / (1 + ||B^-1 a_j||^2)`: the steepest edge in the
    /// space of all variables (exact, from the tableau column norms).
    SteepestEdge,
    /// Smallest index with negative reduced cost (Bland's rule; slow but
    /// cycle-free).
    Bland,
}

/// Result of [`simplex_with`]: the status plus pivoting statistics.
#[derive(Debug, Clone, PartialEq)]
pub struct LpResult {
    /// Outcome of the program.
    pub status: LpStatus,
    /// Total number of pivots over both phases.
    pub iterations: usize,
    /// Number of pivots made under the Bland anti-cycling fallback.
    pub bland_pivots: usize,
}

/// Two-phase simplex: minimise `c . x` subject to rows `a_i . x (rel_i) b_i`
/// and `x >= 0`, with steepest-edge pricing and the Bland anti-cycling
/// fallback; see [`simplex_with`].
#[must_use]
pub fn simplex(
    c: &[f64],
    a: &[Vec<f64>],
    rel: &[Relation],
    b: &[f64],
) -> LpStatus {
    simplex_with(c, a, rel, b, Pricing::SteepestEdge).status
}

/// Two-phase tableau simplex with a selectable pricing rule.
///
/// For [`Pricing::Dantzig`] and [`Pricing::SteepestEdge`] the solver
/// monitors degeneracy: after a run of consecutive pivots that do not
/// improve the objective it switches to Bland's rule (smallest-index
/// entering variable, ties in the ratio test broken by smallest basic
/// index) until a strictly improving pivot occurs, which guarantees
/// termination on cycling examples such as Beale's. The leaving-variable
/// tie-break is always Bland's.
#[must_use]
pub fn simplex_with(
    c: &[f64],
    a: &[Vec<f64>],
    rel: &[Relation],
    b: &[f64],
    pricing: Pricing,
) -> LpResult {
    let iters = std::cell::Cell::new(0usize);
    let bland_n = std::cell::Cell::new(0usize);
    let status = simplex_core(c, a, rel, b, pricing, &iters, &bland_n);
    LpResult {
        status,
        iterations: iters.get(),
        bland_pivots: bland_n.get(),
    }
}

fn simplex_core(
    c: &[f64],
    a: &[Vec<f64>],
    rel: &[Relation],
    b: &[f64],
    pricing: Pricing,
    iters: &std::cell::Cell<usize>,
    bland_n: &std::cell::Cell<usize>,
) -> LpStatus {
    let n = c.len();
    let m = a.len();
    // Normalise so that b >= 0.
    let mut rows: Vec<(Vec<f64>, Relation, f64)> = (0..m)
        .map(|i| {
            if b[i] < 0.0 {
                let flipped = match rel[i] {
                    | Relation::Le => Relation::Ge,
                    | Relation::Ge => Relation::Le,
                    | Relation::Eq => Relation::Eq,
                };
                (a[i].iter().map(|v| -v).collect(), flipped, -b[i])
            } else {
                (a[i].clone(), rel[i], b[i])
            }
        })
        .collect();
    let nslack = rows.iter().filter(|r| r.1 != Relation::Eq).count();
    let nart = rows.iter().filter(|r| r.1 != Relation::Le).count();
    let total = n + nslack + nart;
    let mut t = vec![vec![0.0; total + 1]; m];
    let mut basis = vec![0; m];
    let (mut si, mut ai) = (n, n + nslack);
    for (i, (row, r, bv)) in rows.iter_mut().enumerate() {
        t[i][..n].copy_from_slice(row);
        t[i][total] = *bv;
        match r {
            | Relation::Le => {
                t[i][si] = 1.0;
                basis[i] = si;
                si += 1;
            },
            | Relation::Ge => {
                t[i][si] = -1.0;
                si += 1;
                t[i][ai] = 1.0;
                basis[i] = ai;
                ai += 1;
            },
            | Relation::Eq => {
                t[i][ai] = 1.0;
                basis[i] = ai;
                ai += 1;
            },
        }
    }
    let eps = 1e-9;
    // Runs the simplex on cost vector `cost` (length `total`) restricted to columns < `limit`.
    let run = |t: &mut Vec<Vec<f64>>, basis: &mut Vec<usize>, cost: &[f64], limit: usize| -> bool {
        // reduced costs row
        let mut rc: Vec<f64> = (0..=total)
            .map(|j| if j < total { cost[j] } else { 0.0 })
            .collect();
        for i in 0..m {
            let cb = cost[basis[i]];
            if cb != 0.0 {
                for j in 0..=total {
                    rc[j] -= cb * t[i][j];
                }
            }
        }
        let mut degenerate_run = 0usize;
        let stall_limit = 2 * (m + 1);
        let mut last_obj = rc[total];
        for _ in 0..50_000 {
            let use_bland = pricing == Pricing::Bland || degenerate_run >= stall_limit;
            let entering = if use_bland {
                (0..limit).find(|&j| rc[j] < -eps)
            } else if pricing == Pricing::Dantzig {
                let mut best: Option<(usize, f64)> = None;
                for j in 0..limit {
                    if rc[j] < -eps && best.is_none_or(|(_, v)| rc[j] < v) {
                        best = Some((j, rc[j]));
                    }
                }
                best.map(|(j, _)| j)
            } else {
                let mut best: Option<(usize, f64)> = None;
                for j in 0..limit {
                    if rc[j] < -eps {
                        let gamma = 1.0 + (0..m).map(|i| t[i][j] * t[i][j]).sum::<f64>();
                        let score = rc[j] * rc[j] / gamma;
                        if best.is_none_or(|(_, v)| score > v) {
                            best = Some((j, score));
                        }
                    }
                }
                best.map(|(j, _)| j)
            };
            let Some(e) = entering else { return true };
            iters.set(iters.get() + 1);
            if use_bland {
                bland_n.set(bland_n.get() + 1);
            }
            let mut best: Option<(usize, f64)> = None;
            for i in 0..m {
                if t[i][e] > eps {
                    let ratio = t[i][total] / t[i][e];
                    let better = match best {
                        | None => true,
                        | Some((bi, br)) => {
                            ratio < br - 1e-12 || (ratio <= br + 1e-12 && basis[i] < basis[bi])
                        },
                    };
                    if better {
                        best = Some((i, ratio));
                    }
                }
            }
            let Some((r, _)) = best else { return false };
            let p = t[r][e];
            for j in 0..=total {
                t[r][j] /= p;
            }
            for i in 0..m {
                if i != r && t[i][e] != 0.0 {
                    let fct = t[i][e];
                    for j in 0..=total {
                        t[i][j] -= fct * t[r][j];
                    }
                }
            }
            let fct = rc[e];
            for j in 0..=total {
                rc[j] -= fct * t[r][j];
            }
            basis[r] = e;
            // rc[total] tracks minus the objective: reset the stall counter on strict progress
            if (rc[total] - last_obj).abs() > 1e-12 * (1.0 + last_obj.abs()) {
                degenerate_run = 0;
                last_obj = rc[total];
            } else {
                degenerate_run += 1;
            }
        }
        true
    };
    if nart > 0 {
        let mut c1 = vec![0.0; total];
        for v in &mut c1[n + nslack..] {
            *v = 1.0;
        }
        run(&mut t, &mut basis, &c1, total);
        let infeas: f64 = (0..m)
            .filter(|&i| basis[i] >= n + nslack)
            .map(|i| t[i][total])
            .sum();
        if infeas > 1e-7 {
            return LpStatus::Infeasible;
        }
        // Pivot remaining (degenerate) artificials out of the basis.
        for i in 0..m {
            if basis[i] >= n + nslack
                && let Some(e) = (0..n + nslack).find(|&j| t[i][j].abs() > eps)
            {
                let p = t[i][e];
                for j in 0..=total {
                    t[i][j] /= p;
                }
                for k in 0..m {
                    if k != i && t[k][e] != 0.0 {
                        let fct = t[k][e];
                        for j in 0..=total {
                            t[k][j] -= fct * t[i][j];
                        }
                    }
                }
                basis[i] = e;
            }
        }
    }
    let mut c2 = vec![0.0; total];
    c2[..n].copy_from_slice(c);
    if !run(&mut t, &mut basis, &c2, n + nslack) {
        return LpStatus::Unbounded;
    }
    let mut x = vec![0.0; n];
    for i in 0..m {
        if basis[i] < n {
            x[basis[i]] = t[i][total];
        }
    }
    let value = dot(c, &x);
    LpStatus::Optimal { x, value }
}
