//! # Optimisation kernels
//!
//! * unconstrained: Nelder-Mead, BFGS, L-BFGS, dogleg trust region;
//! * global, derivative free: differential evolution and simulated
//!   annealing (seeded, deterministic);
//! * constrained: augmented Lagrangian for equality and inequality
//!   constraints;
//! * linear programming: two-phase simplex method with Bland's rule.
#![allow(
    clippy::manual_midpoint,
    clippy::missing_const_for_fn,
    clippy::struct_field_names,
    clippy::or_fun_call,
    clippy::manual_swap,
    clippy::if_not_else,
    clippy::unnecessary_sort_by,
    clippy::while_float,
    clippy::too_long_first_doc_paragraph,
    clippy::cast_sign_loss,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::needless_pass_by_value,
    clippy::manual_map,
    clippy::unnecessary_map_or,
    clippy::suboptimal_flops,
    clippy::similar_names,
    clippy::unreadable_literal,
    clippy::excessive_precision,
    clippy::needless_range_loop,
    clippy::float_cmp,
    clippy::too_many_lines,
    clippy::cognitive_complexity,
    clippy::option_if_let_else,
    clippy::many_single_char_names,
    clippy::type_complexity,
    clippy::too_many_arguments,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_precision_loss,
    clippy::map_unwrap_or,
    clippy::cast_lossless
)]

use crate::kernels::dense::{self, Mat};
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
pub fn numerical_gradient<F: Fn(&[f64]) -> f64>(f: &F, x: &[f64]) -> Vec<f64> {
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

fn dot(a: &[f64], b: &[f64]) -> f64 {
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
        p[i] += if p[i] == 0.0 { step * 0.05_f64.max(step) } else { step * p[i].abs().max(1.0) * p[i].signum() };
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
        let size = pts[1..].iter().map(|p| {
            p.iter().zip(&pts[0]).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max)
        }).fold(0.0, f64::max);
        if spread <= tol * (1.0 + vals[0].abs()) && size <= tol.sqrt() {
            converged = true;
            break;
        }
        let cen: Vec<f64> = (0..n).map(|j| pts[..n].iter().map(|p| p[j]).sum::<f64>() / n as f64).collect();
        let along = |t: f64| -> Vec<f64> { (0..n).map(|j| cen[j] + t * (pts[n][j] - cen[j])).collect() };
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
    let best = (0..=n).min_by(|&a, &b| vals[a].total_cmp(&vals[b])).unwrap_or(0);
    OptResult { x: pts[best].clone(), fx: vals[best], iterations: it, evaluations: evals, converged }
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
        t = if hi.is_finite() { 0.5 * (lo + hi) } else { 2.0 * lo.max(t) };
    }
    best
}

/// BFGS quasi-Newton method with a weak-Wolfe line search. If `grad` is
/// `None` a central-difference gradient is used.
pub fn bfgs<F: Fn(&[f64]) -> f64, G: Fn(&[f64]) -> Vec<f64>>(
    f: F,
    grad: Option<G>,
    x0: &[f64],
    tol: f64,
    max_iter: usize,
) -> OptResult {
    let gfun = |x: &[f64]| grad.as_ref().map_or_else(|| numerical_gradient(&f, x), |g| g(x));
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
    OptResult { x, fx, iterations: it, evaluations: evals, converged }
}

/// Limited-memory BFGS with history `m` (two-loop recursion).
pub fn lbfgs<F: Fn(&[f64]) -> f64, G: Fn(&[f64]) -> Vec<f64>>(
    f: F,
    grad: Option<G>,
    x0: &[f64],
    m: usize,
    tol: f64,
    max_iter: usize,
) -> OptResult {
    let gfun = |x: &[f64]| grad.as_ref().map_or_else(|| numerical_gradient(&f, x), |g| g(x));
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
    OptResult { x, fx, iterations: it, evaluations: evals, converged }
}

/// Trust-region method with the dogleg step. The Hessian is supplied
/// by `hess`; pass `None` for a finite-difference Hessian from gradients.
pub fn trust_region<F, G, H>(
    f: F,
    grad: G,
    hess: Option<H>,
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
        hess.as_ref().map_or_else(
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
                        let v = 0.5 * (hm.at(i, j) + hm.at(j, i));
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
        let tau = if gbg <= 0.0 { f64::INFINITY } else { gn2 / gbg };
        let pu: Vec<f64> = g.iter().map(|v| -tau * v).collect();
        let neg: Vec<f64> = g.iter().map(|v| -v).collect();
        let pb = dense::solve(&b, &neg).ok().filter(|p| dot(p, &b.matvec(p)) > 0.0);
        let norm = |v: &[f64]| dot(v, v).sqrt();
        let step: Vec<f64> = match pb {
            Some(pb) if norm(&pb) <= radius => pb,
            Some(pb) if tau.is_finite() && norm(&pu) < radius => {
                // intersection of the segment pu -> pb with the boundary
                let d: Vec<f64> = pb.iter().zip(&pu).map(|(a, c)| a - c).collect();
                let (a2, b2, c2) = (dot(&d, &d), 2.0 * dot(&pu, &d), dot(&pu, &pu) - radius * radius);
                let s = (-b2 + (b2 * b2 - 4.0 * a2 * c2).max(0.0).sqrt()) / (2.0 * a2);
                pu.iter().zip(&d).map(|(a, c)| a + s * c).collect()
            }
            _ => {
                let gn = gn2.sqrt();
                g.iter().map(|v| -radius * v / gn).collect()
            }
        };
        let bp = b.matvec(&step);
        let pred = -(dot(&g, &step) + 0.5 * dot(&step, &bp));
        let xn: Vec<f64> = x.iter().zip(&step).map(|(a, c)| a + c).collect();
        let fnw = f(&xn);
        evals += 1;
        let rho = if pred.abs() < 1e-300 { 0.0 } else { (fx - fnw) / pred };
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
    OptResult { x, fx, iterations: it, evaluations: evals, converged }
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
        .map(|_| bounds.iter().map(|&(a, b)| rng.uniform_range(a, b)).collect())
        .collect();
    let mut vals: Vec<f64> = pop.iter().map(|p| f(p)).collect();
    let mut evals = np;
    let (fw, cr) = (0.7, 0.9);
    let mut generation = 0;
    while generation < generations {
        generation += 1;
        for i in 0..np {
            let mut pick = || loop {
                let k = rng.below(np as u64) as usize;
                if k != i {
                    return k;
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
    let best = (0..np).min_by(|&a, &b| vals[a].total_cmp(&vals[b])).unwrap_or(0);
    OptResult { x: pop[best].clone(), fx: vals[best], iterations: generation, evaluations: evals, converged: generation < generations }
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
    OptResult { x: best.0, fx: best.1, iterations, evaluations: iterations + 1, converged: true }
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
        let r = bfgs(phi, None::<fn(&[f64]) -> Vec<f64>>, &x, 1e-9, 500);
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
        let viol = h.iter().map(|v| v.abs()).chain(g.iter().map(|v| v.max(0.0))).fold(0.0, f64::max);
        if viol <= tol {
            let comp = g.iter().zip(&mu).map(|(a, b)| (a * b).abs()).fold(0.0, f64::max);
            if comp <= tol.max(1e-6) {
                converged = true;
                break;
            }
        }
        rho = (rho * 4.0).min(1e8);
    }
    let fx = f(&x);
    OptResult { x, fx, iterations: outer, evaluations: evals, converged }
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

/// Two-phase simplex: minimise `c . x` subject to rows `a_i . x (rel_i) b_i`
/// and `x >= 0`.
#[must_use]
pub fn simplex(c: &[f64], a: &[Vec<f64>], rel: &[Relation], b: &[f64]) -> LpStatus {
    let n = c.len();
    let m = a.len();
    // Normalise so that b >= 0.
    let mut rows: Vec<(Vec<f64>, Relation, f64)> = (0..m)
        .map(|i| {
            if b[i] < 0.0 {
                let flipped = match rel[i] {
                    Relation::Le => Relation::Ge,
                    Relation::Ge => Relation::Le,
                    Relation::Eq => Relation::Eq,
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
            Relation::Le => {
                t[i][si] = 1.0;
                basis[i] = si;
                si += 1;
            }
            Relation::Ge => {
                t[i][si] = -1.0;
                si += 1;
                t[i][ai] = 1.0;
                basis[i] = ai;
                ai += 1;
            }
            Relation::Eq => {
                t[i][ai] = 1.0;
                basis[i] = ai;
                ai += 1;
            }
        }
    }
    let eps = 1e-9;
    // Runs the simplex on cost vector `cost` (length `total`) restricted to columns < `limit`.
    let run = |t: &mut Vec<Vec<f64>>, basis: &mut Vec<usize>, cost: &[f64], limit: usize| -> bool {
        // reduced costs row
        let mut rc: Vec<f64> = (0..=total).map(|j| if j < total { cost[j] } else { 0.0 }).collect();
        for i in 0..m {
            let cb = cost[basis[i]];
            if cb != 0.0 {
                for j in 0..=total {
                    rc[j] -= cb * t[i][j];
                }
            }
        }
        for _ in 0..50_000 {
            let Some(e) = (0..limit).find(|&j| rc[j] < -eps) else { return true };
            let mut best: Option<(usize, f64)> = None;
            for i in 0..m {
                if t[i][e] > eps {
                    let ratio = t[i][total] / t[i][e];
                    let better = match best {
                        None => true,
                        Some((bi, br)) => ratio < br - 1e-12 || (ratio <= br + 1e-12 && basis[i] < basis[bi]),
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
        }
        true
    };
    if nart > 0 {
        let mut c1 = vec![0.0; total];
        for v in &mut c1[n + nslack..] {
            *v = 1.0;
        }
        run(&mut t, &mut basis, &c1, total);
        let infeas: f64 = (0..m).filter(|&i| basis[i] >= n + nslack).map(|i| t[i][total]).sum();
        if infeas > 1e-7 {
            return LpStatus::Infeasible;
        }
        // Pivot remaining (degenerate) artificials out of the basis.
        for i in 0..m {
            if basis[i] >= n + nslack {
                if let Some(e) = (0..n + nslack).find(|&j| t[i][j].abs() > eps) {
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
