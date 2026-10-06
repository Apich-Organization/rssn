//! # Adaptive, stiff, symplectic and boundary-value ODE kernels
//!
//! * [`dopri5`]: Dormand-Prince 5(4) with step-size control, a fourth-order
//!   continuous extension ("dense output") and event detection;
//! * [`rosenbrock23`]: the L-stable Rosenbrock pair of Shampine and
//!   Reichelt for stiff problems (finite-difference Jacobian);
//! * [`leapfrog`] and [`yoshida4`]: symplectic integrators for
//!   `H = p^2/2 + V(q)`;
//! * [`bvp_fd`]: second-order finite-difference BVP solver for
//!   `y'' = f(x, y, y')` with Dirichlet conditions (Newton iteration);
//!   [`shooting`]: general shooting method;
//! * [`dae_index1`]: backward-Euler for semi-explicit index-1 DAEs;
//! * [`dde_rk4`]: constant-delay differential equations by the method of
//!   steps with Hermite history interpolation.
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
use crate::kernels::rootfind;

/// Errors of the ODE kernels.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OdeError {
    /// Step size underflow or too many steps.
    StepLimit,
    /// The state became non-finite.
    NonFinite,
    /// A nonlinear / linear solve failed.
    SolveFailed,
    /// Invalid arguments.
    Invalid,
}

/// Tolerances and limits for adaptive integrators.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OdeOptions {
    /// Relative tolerance.
    pub rtol: f64,
    /// Absolute tolerance.
    pub atol: f64,
    /// Initial step (`None` selects one automatically).
    pub h0: Option<f64>,
    /// Maximal step size.
    pub hmax: f64,
    /// Maximal number of steps.
    pub max_steps: usize,
}

impl Default for OdeOptions {
    fn default() -> Self {
        Self { rtol: 1e-8, atol: 1e-10, h0: None, hmax: f64::INFINITY, max_steps: 1_000_000 }
    }
}

/// An event function `g(t, y)`; a sign change of `g` is an event.
pub struct Event<'a> {
    /// The event function.
    pub g: &'a dyn Fn(f64, &[f64]) -> f64,
    /// Stop the integration at the event.
    pub terminal: bool,
    /// `1` rising only, `-1` falling only, `0` both.
    pub direction: i8,
}

/// A located event.
#[derive(Debug, Clone, PartialEq)]
pub struct EventHit {
    /// Index into the event list.
    pub index: usize,
    /// Event time.
    pub t: f64,
    /// State at the event time.
    pub y: Vec<f64>,
}

#[derive(Debug, Clone)]
struct DenseStep {
    t0: f64,
    h: f64,
    y0: Vec<f64>,
    d: [Vec<f64>; 4],
}

/// Solution of an adaptive integration.
#[derive(Debug, Clone)]
pub struct OdeSolution {
    /// Accepted time points.
    pub t: Vec<f64>,
    /// States at the accepted time points.
    pub y: Vec<Vec<f64>>,
    /// Located events, in time order.
    pub events: Vec<EventHit>,
    /// Number of right-hand-side evaluations.
    pub nfev: usize,
    /// Number of rejected steps.
    pub rejected: usize,
    dense: Vec<DenseStep>,
}

impl OdeSolution {
    /// Evaluates the continuous extension at `t` (available for
    /// [`dopri5`]); `None` outside the integrated range.
    #[must_use]
    pub fn eval(&self, t: f64) -> Option<Vec<f64>> {
        let (&first, &last) = (self.t.first()?, self.t.last()?);
        let (lo, hi) = if first <= last { (first, last) } else { (last, first) };
        if t < lo || t > hi || self.dense.is_empty() {
            return None;
        }
        let idx = self
            .dense
            .partition_point(|s| (s.t0 + s.h - t) * s.h.signum() < 0.0)
            .min(self.dense.len() - 1);
        Some(dense_eval(&self.dense[idx], t))
    }
}

fn dense_eval(s: &DenseStep, t: f64) -> Vec<f64> {
    let x = (t - s.t0) / s.h;
    let (x1, x2, x3, x4) = (x, x * x, x * x * x, x * x * x * x);
    (0..s.y0.len())
        .map(|i| {
            s.y0[i] + s.h * (s.d[0][i] * x1 + s.d[1][i] * x2 + s.d[2][i] * x3 + s.d[3][i] * x4)
        })
        .collect()
}

const DP_A: [[f64; 6]; 6] = [
    [1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0],
    [44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0],
    [19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0, 0.0, 0.0],
    [9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0, 0.0],
    [35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0],
];
const DP_C: [f64; 6] = [0.2, 0.3, 0.8, 8.0 / 9.0, 1.0, 1.0];
const DP_E: [f64; 7] = [
    71.0 / 57600.0,
    0.0,
    -71.0 / 16695.0,
    71.0 / 1920.0,
    -17253.0 / 339200.0,
    22.0 / 525.0,
    -1.0 / 40.0,
];
const DP_P: [[f64; 4]; 7] = [
    [1.0, -8048581381.0 / 2820520608.0, 8663915743.0 / 2820520608.0, -12715105075.0 / 11282082432.0],
    [0.0, 0.0, 0.0, 0.0],
    [0.0, 131558114200.0 / 32700410799.0, -68118460800.0 / 10900136933.0, 87487479700.0 / 32700410799.0],
    [0.0, -1754552775.0 / 470086768.0, 14199869525.0 / 1410260304.0, -10690763975.0 / 1880347072.0],
    [0.0, 127303824393.0 / 49829197408.0, -318862633887.0 / 49829197408.0, 701980252875.0 / 199316789632.0],
    [0.0, -282668133.0 / 205662961.0, 2019193451.0 / 616988883.0, -1453857185.0 / 822651844.0],
    [0.0, 40617522.0 / 29380423.0, -110615467.0 / 29380423.0, 69997945.0 / 29380423.0],
];

fn err_norm(err: &[f64], y0: &[f64], y1: &[f64], rtol: f64, atol: f64) -> f64 {
    let n = err.len().max(1) as f64;
    (err.iter()
        .zip(y0.iter().zip(y1))
        .map(|(e, (a, b))| {
            let sc = atol + rtol * a.abs().max(b.abs());
            (e / sc) * (e / sc)
        })
        .sum::<f64>()
        / n)
        .sqrt()
}

fn initial_step<F: Fn(f64, &[f64], &mut [f64])>(
    f: &F,
    t0: f64,
    y0: &[f64],
    f0: &[f64],
    dir: f64,
    order: f64,
    rtol: f64,
    atol: f64,
) -> f64 {
    let n = y0.len();
    let sc: Vec<f64> = y0.iter().map(|v| atol + rtol * v.abs()).collect();
    let d0 = (y0.iter().zip(&sc).map(|(a, s)| (a / s) * (a / s)).sum::<f64>() / n as f64).sqrt();
    let d1 = (f0.iter().zip(&sc).map(|(a, s)| (a / s) * (a / s)).sum::<f64>() / n as f64).sqrt();
    let h0 = if d0 < 1e-5 || d1 < 1e-5 { 1e-6 } else { 0.01 * d0 / d1 };
    let y1: Vec<f64> = y0.iter().zip(f0).map(|(a, b)| a + dir * h0 * b).collect();
    let mut f1 = vec![0.0; n];
    f(t0 + dir * h0, &y1, &mut f1);
    let d2 = (f1.iter().zip(f0).zip(&sc).map(|((a, b), s)| ((a - b) / s).powi(2)).sum::<f64>()
        / n as f64)
        .sqrt()
        / h0;
    let h1 = if d1.max(d2) <= 1e-15 {
        (h0 * 1e-3).max(1e-6)
    } else {
        (0.01 / d1.max(d2)).powf(1.0 / (order + 1.0))
    };
    (100.0 * h0).min(h1)
}

/// Dormand-Prince 5(4) integration of `y' = f(t, y)` from `t0` to `t1`
/// (either direction) with dense output and event detection.
///
/// `f(t, y, dy)` writes the derivative into `dy`. Events are located by
/// Brent iteration on the continuous extension; a terminal event ends the
/// integration at the event time.
///
/// # Errors
/// [`OdeError::StepLimit`] or [`OdeError::NonFinite`].
pub fn dopri5<F: Fn(f64, &[f64], &mut [f64])>(
    f: F,
    t0: f64,
    t1: f64,
    y0: &[f64],
    opts: &OdeOptions,
    events: &[Event<'_>],
) -> Result<OdeSolution, OdeError> {
    let n = y0.len();
    if n == 0 || t0 == t1 {
        return Err(OdeError::Invalid);
    }
    let dir = if t1 > t0 { 1.0 } else { -1.0 };
    let mut t = t0;
    let mut y = y0.to_vec();
    let mut k = vec![vec![0.0; n]; 7];
    f(t, &y, &mut k[0]);
    let mut nfev = 1;
    let mut h = opts
        .h0
        .unwrap_or_else(|| initial_step(&f, t0, y0, &k[0], dir, 5.0, opts.rtol, opts.atol))
        .abs()
        .min(opts.hmax)
        * dir;
    let mut sol = OdeSolution {
        t: vec![t0],
        y: vec![y0.to_vec()],
        events: Vec::new(),
        nfev: 0,
        rejected: 0,
        dense: Vec::new(),
    };
    let mut gprev: Vec<f64> = events.iter().map(|e| (e.g)(t0, y0)).collect();
    let mut steps = 0;
    let mut tmp = vec![0.0; n];
    while (t1 - t) * dir > 0.0 {
        if steps >= opts.max_steps {
            return Err(OdeError::StepLimit);
        }
        steps += 1;
        if (t + h - t1) * dir > 0.0 {
            h = t1 - t;
        }
        // stages
        for s in 0..6 {
            for i in 0..n {
                let mut acc = 0.0;
                for j in 0..=s {
                    if DP_A[s][j] != 0.0 {
                        acc += DP_A[s][j] * k[j][i];
                    }
                }
                tmp[i] = y[i] + h * acc;
            }
            let (head, tail) = k.split_at_mut(s + 1);
            let _ = head;
            f(t + DP_C[s] * h, &tmp, &mut tail[0]);
            nfev += 1;
        }
        let ynew = tmp.clone();
        // after stage s=5 tmp holds y + h*sum(b k) = y_new (row 6 of A is b)
        let mut err = vec![0.0; n];
        for i in 0..n {
            let mut e = 0.0;
            for j in 0..7 {
                e += DP_E[j] * k[j][i];
            }
            err[i] = h * e;
        }
        if ynew.iter().any(|v| !v.is_finite()) {
            h *= 0.25;
            sol.rejected += 1;
            if h.abs() < 1e-300 {
                return Err(OdeError::NonFinite);
            }
            continue;
        }
        let en = err_norm(&err, &y, &ynew, opts.rtol, opts.atol);
        if en <= 1.0 {
            // dense coefficients
            let d: [Vec<f64>; 4] = std::array::from_fn(|c| {
                (0..n).map(|i| (0..7).map(|j| DP_P[j][c] * k[j][i]).sum()).collect()
            });
            let step = DenseStep { t0: t, h, y0: y.clone(), d };
            let tn = t + h;
            // events
            let mut hits: Vec<EventHit> = Vec::new();
            let mut stop = false;
            for (ei, ev) in events.iter().enumerate() {
                let g1 = (ev.g)(tn, &ynew);
                let g0 = gprev[ei];
                let crossing = (g0 < 0.0 && g1 >= 0.0 && ev.direction >= 0)
                    || (g0 > 0.0 && g1 <= 0.0 && ev.direction <= 0);
                if crossing {
                    let gf = |tt: f64| (ev.g)(tt, &dense_eval(&step, tt));
                    let root = rootfind::brent(gf, t, tn, 1e-14 * (1.0 + tn.abs()), 200).unwrap_or(tn);
                    hits.push(EventHit { index: ei, t: root, y: dense_eval(&step, root) });
                    if ev.terminal {
                        stop = true;
                    }
                }
                gprev[ei] = g1;
            }
            hits.sort_by(|a, b| ((a.t - b.t) * dir).total_cmp(&0.0));
            if stop {
                // keep events up to and including the first terminal one
                let mut kept = Vec::new();
                for h in hits {
                    let term = events[h.index].terminal;
                    kept.push(h);
                    if term {
                        break;
                    }
                }
                let last = kept.last().cloned();
                sol.events.extend(kept);
                sol.dense.push(step);
                if let Some(l) = last {
                    sol.t.push(l.t);
                    sol.y.push(l.y);
                }
                sol.nfev = nfev;
                return Ok(sol);
            }
            sol.events.extend(hits);
            sol.dense.push(step);
            t = tn;
            y = ynew;
            sol.t.push(t);
            sol.y.push(y.clone());
            // FSAL
            let last = k[6].clone();
            k[0] = last;
            let fac = if en == 0.0 { 5.0 } else { (0.9 * en.powf(-0.2)).clamp(0.2, 5.0) };
            h *= fac;
            if h.abs() > opts.hmax {
                h = opts.hmax * dir;
            }
        } else {
            sol.rejected += 1;
            h *= (0.9 * en.powf(-0.2)).clamp(0.1, 0.9);
            if h.abs() < 1e-14 * (1.0 + t.abs()) {
                return Err(OdeError::StepLimit);
            }
        }
    }
    sol.nfev = nfev;
    Ok(sol)
}

/// Rosenbrock 2(3) (the `ode23s` pair) for stiff problems.
///
/// The Jacobian and the explicit time derivative are obtained by finite
/// differences at every step. No dense output is provided.
///
/// # Errors
/// [`OdeError::StepLimit`], [`OdeError::NonFinite`] or [`OdeError::SolveFailed`].
pub fn rosenbrock23<F: Fn(f64, &[f64], &mut [f64])>(
    f: F,
    t0: f64,
    t1: f64,
    y0: &[f64],
    opts: &OdeOptions,
) -> Result<OdeSolution, OdeError> {
    let n = y0.len();
    if n == 0 || t0 == t1 {
        return Err(OdeError::Invalid);
    }
    let dir = if t1 > t0 { 1.0 } else { -1.0 };
    let d = 1.0 / (2.0 + 2.0_f64.sqrt());
    let e32 = 6.0 + 2.0_f64.sqrt();
    let mut t = t0;
    let mut y = y0.to_vec();
    let mut f0 = vec![0.0; n];
    f(t, &y, &mut f0);
    let mut nfev = 1;
    let mut h = opts
        .h0
        .unwrap_or_else(|| initial_step(&f, t0, y0, &f0, dir, 2.0, opts.rtol, opts.atol))
        .abs()
        .min(opts.hmax)
        * dir;
    let mut sol = OdeSolution {
        t: vec![t0],
        y: vec![y0.to_vec()],
        events: Vec::new(),
        nfev: 0,
        rejected: 0,
        dense: Vec::new(),
    };
    let mut steps = 0;
    while (t1 - t) * dir > 0.0 {
        if steps >= opts.max_steps {
            return Err(OdeError::StepLimit);
        }
        steps += 1;
        if (t + h - t1) * dir > 0.0 {
            h = t1 - t;
        }
        // Jacobian J and dF/dt by forward differences.
        let mut jac = Mat::zeros(n, n);
        let mut yp = y.clone();
        let mut fp = vec![0.0; n];
        for j in 0..n {
            let dj = 1.5e-8 * y[j].abs().max(1e-5);
            yp[j] = y[j] + dj;
            f(t, &yp, &mut fp);
            yp[j] = y[j];
            for i in 0..n {
                jac.set(i, j, (fp[i] - f0[i]) / dj);
            }
        }
        let dt = 1.5e-8 * t.abs().max(1e-5);
        f(t + dt, &y, &mut fp);
        let tder: Vec<f64> = (0..n).map(|i| (fp[i] - f0[i]) / dt).collect();
        nfev += n + 1;
        let mut w = Mat::identity(n);
        for i in 0..n {
            for j in 0..n {
                let v = w.at(i, j) - h * d * jac.at(i, j);
                w.set(i, j, v);
            }
        }
        let Ok(lu) = dense::lu_factor(&w) else {
            h *= 0.5;
            sol.rejected += 1;
            if h.abs() < 1e-14 {
                return Err(OdeError::SolveFailed);
            }
            continue;
        };
        let rhs1: Vec<f64> = (0..n).map(|i| f0[i] + h * d * tder[i]).collect();
        let k1 = lu.solve(&rhs1);
        let ymid: Vec<f64> = (0..n).map(|i| y[i] + 0.5 * h * k1[i]).collect();
        let mut f1 = vec![0.0; n];
        f(t + 0.5 * h, &ymid, &mut f1);
        let r2: Vec<f64> = (0..n).map(|i| f1[i] - k1[i]).collect();
        let s2 = lu.solve(&r2);
        let k2: Vec<f64> = (0..n).map(|i| s2[i] + k1[i]).collect();
        let ynew: Vec<f64> = (0..n).map(|i| y[i] + h * k2[i]).collect();
        let mut f2 = vec![0.0; n];
        f(t + h, &ynew, &mut f2);
        nfev += 2;
        let r3: Vec<f64> = (0..n)
            .map(|i| f2[i] - e32 * (k2[i] - f1[i]) - 2.0 * (k1[i] - f0[i]) + h * d * tder[i])
            .collect();
        let k3 = lu.solve(&r3);
        let err: Vec<f64> = (0..n).map(|i| h / 6.0 * (k1[i] - 2.0 * k2[i] + k3[i])).collect();
        if ynew.iter().any(|v| !v.is_finite()) {
            h *= 0.25;
            sol.rejected += 1;
            if h.abs() < 1e-300 {
                return Err(OdeError::NonFinite);
            }
            continue;
        }
        let en = err_norm(&err, &y, &ynew, opts.rtol, opts.atol);
        if en <= 1.0 {
            t += h;
            y = ynew;
            f0 = f2;
            sol.t.push(t);
            sol.y.push(y.clone());
            let fac = if en == 0.0 { 5.0 } else { (0.9 * en.powf(-1.0 / 3.0)).clamp(0.2, 5.0) };
            h *= fac;
            if h.abs() > opts.hmax {
                h = opts.hmax * dir;
            }
        } else {
            sol.rejected += 1;
            h *= (0.9 * en.powf(-1.0 / 3.0)).clamp(0.1, 0.9);
            if h.abs() < 1e-14 * (1.0 + t.abs()) {
                return Err(OdeError::StepLimit);
            }
        }
    }
    sol.nfev = nfev;
    Ok(sol)
}

/// Velocity-Verlet (leapfrog) for `q'' = acc(q)`; returns `(q, p)` after
/// `steps` steps of size `dt` (unit mass, `p = q'`).
pub fn leapfrog<A: Fn(&[f64]) -> Vec<f64>>(
    acc: A,
    q0: &[f64],
    p0: &[f64],
    dt: f64,
    steps: usize,
) -> (Vec<f64>, Vec<f64>) {
    let mut q = q0.to_vec();
    let mut p = p0.to_vec();
    let mut a = acc(&q);
    for _ in 0..steps {
        for i in 0..q.len() {
            p[i] += 0.5 * dt * a[i];
            q[i] += dt * p[i];
        }
        a = acc(&q);
        for i in 0..q.len() {
            p[i] += 0.5 * dt * a[i];
        }
    }
    (q, p)
}

/// Fourth-order symplectic Yoshida (triple-jump) composition of leapfrog.
pub fn yoshida4<A: Fn(&[f64]) -> Vec<f64>>(
    acc: A,
    q0: &[f64],
    p0: &[f64],
    dt: f64,
    steps: usize,
) -> (Vec<f64>, Vec<f64>) {
    let c = 2.0_f64.cbrt();
    let w1 = 1.0 / (2.0 - c);
    let w0 = -c / (2.0 - c);
    let mut q = q0.to_vec();
    let mut p = p0.to_vec();
    for _ in 0..steps {
        for w in [w1, w0, w1] {
            let (qn, pn) = leapfrog(&acc, &q, &p, w * dt, 1);
            q = qn;
            p = pn;
        }
    }
    (q, p)
}

/// Second-order finite-difference solution of the BVP
/// `y'' = f(x, y, y')`, `y(a) = ya`, `y(b) = yb` on `n` interior points,
/// by Newton iteration on the discrete system. Returns `(x, y)` including
/// the end points; `guess` gives the initial profile.
///
/// # Errors
/// [`OdeError::SolveFailed`] if Newton fails.
pub fn bvp_fd<F: Fn(f64, f64, f64) -> f64, G: Fn(f64) -> f64>(
    f: F,
    a: f64,
    b: f64,
    ya: f64,
    yb: f64,
    n: usize,
    guess: G,
) -> Result<(Vec<f64>, Vec<f64>), OdeError> {
    if n < 1 {
        return Err(OdeError::Invalid);
    }
    let h = (b - a) / (n + 1) as f64;
    let x: Vec<f64> = (0..n + 2).map(|i| a + h * i as f64).collect();
    let resid = |u: &[f64]| -> Vec<f64> {
        let full = |i: usize| -> f64 {
            if i == 0 { ya } else if i == n + 1 { yb } else { u[i - 1] }
        };
        (1..=n)
            .map(|i| {
                let (ym, y0, yp) = (full(i - 1), full(i), full(i + 1));
                (yp - 2.0 * y0 + ym) / (h * h) - f(x[i], y0, (yp - ym) / (2.0 * h))
            })
            .collect()
    };
    let u0: Vec<f64> = (1..=n).map(|i| guess(x[i])).collect();
    let sol = rootfind::newton_system(resid, None::<fn(&[f64]) -> Mat>, &u0, 1e-11, 100)
        .map_err(|_| OdeError::SolveFailed)?;
    let mut y = vec![ya];
    y.extend(sol.x);
    y.push(yb);
    Ok((x, y))
}

/// General shooting method: `params -> initial state` and
/// `final state -> residual` define the boundary conditions; the
/// parameters are adjusted by Newton iteration so the residual vanishes.
///
/// # Errors
/// [`OdeError::SolveFailed`] if the parameter iteration or an inner
/// integration fails.
pub fn shooting<F, I, R>(
    f: F,
    t0: f64,
    t1: f64,
    init: I,
    residual: R,
    p0: &[f64],
    opts: &OdeOptions,
) -> Result<Vec<f64>, OdeError>
where
    F: Fn(f64, &[f64], &mut [f64]),
    I: Fn(&[f64]) -> Vec<f64>,
    R: Fn(&[f64]) -> Vec<f64>,
{
    let shoot = |p: &[f64]| -> Vec<f64> {
        let y0 = init(p);
        match dopri5(&f, t0, t1, &y0, opts, &[]) {
            Ok(s) => residual(s.y.last().map_or(&[][..], Vec::as_slice)),
            Err(_) => vec![f64::NAN; p.len()],
        }
    };
    rootfind::newton_system(shoot, None::<fn(&[f64]) -> Mat>, p0, 1e-9, 60)
        .map(|s| s.x)
        .map_err(|_| OdeError::SolveFailed)
}

/// Backward Euler for the semi-explicit index-1 DAE
/// `y' = f(t, y, z)`, `0 = g(t, y, z)` with fixed step `h`.
/// Returns the trajectory of `(y, z)` concatenated.
///
/// # Errors
/// [`OdeError::SolveFailed`] if the nonlinear solve fails.
pub fn dae_index1<F, G>(
    f: F,
    g: G,
    y0: &[f64],
    z0: &[f64],
    t0: f64,
    h: f64,
    steps: usize,
) -> Result<Vec<Vec<f64>>, OdeError>
where
    F: Fn(f64, &[f64], &[f64]) -> Vec<f64>,
    G: Fn(f64, &[f64], &[f64]) -> Vec<f64>,
{
    let ny = y0.len();
    let mut cur: Vec<f64> = y0.iter().chain(z0).copied().collect();
    let mut out = vec![cur.clone()];
    for s in 0..steps {
        let tn = t0 + (s + 1) as f64 * h;
        let prev = cur.clone();
        let sys = |u: &[f64]| -> Vec<f64> {
            let (y, z) = u.split_at(ny);
            let fy = f(tn, y, z);
            let mut r: Vec<f64> = (0..ny).map(|i| y[i] - prev[i] - h * fy[i]).collect();
            r.extend(g(tn, y, z));
            r
        };
        let sol = rootfind::newton_system(sys, None::<fn(&[f64]) -> Mat>, &prev, 1e-12, 50)
            .map_err(|_| OdeError::SolveFailed)?;
        cur = sol.x;
        out.push(cur.clone());
    }
    Ok(out)
}

/// Constant-delay DDE `y'(t) = f(t, y(t), y(t - tau))` with history
/// `hist(t)` for `t <= t0`, integrated by classical RK4 (method of steps)
/// with cubic Hermite interpolation of the stored solution. Requires
/// `tau >= h`. Returns the grid and states.
///
/// # Errors
/// [`OdeError::Invalid`] if `tau < h` or `h <= 0`.
pub fn dde_rk4<F, H>(
    f: F,
    hist: H,
    tau: f64,
    t0: f64,
    t1: f64,
    h: f64,
) -> Result<(Vec<f64>, Vec<Vec<f64>>), OdeError>
where
    F: Fn(f64, &[f64], &[f64]) -> Vec<f64>,
    H: Fn(f64) -> Vec<f64>,
{
    if h <= 0.0 || tau < h {
        return Err(OdeError::Invalid);
    }
    let mut ts = vec![t0];
    let mut ys = vec![hist(t0)];
    let mut fs = vec![f(t0, &ys[0], &hist(t0 - tau))];
    let n = ys[0].len();
    let delayed = |t: f64, ts: &[f64], ys: &[Vec<f64>], fs: &[Vec<f64>]| -> Vec<f64> {
        if t <= t0 {
            return hist(t);
        }
        let idx = ts.partition_point(|&x| x <= t).saturating_sub(1).min(ts.len() - 2);
        let (ta, tb) = (ts[idx], ts[idx + 1]);
        let s = (t - ta) / (tb - ta);
        let (h00, h10, h01, h11) = (
            (1.0 + 2.0 * s) * (1.0 - s) * (1.0 - s),
            s * (1.0 - s) * (1.0 - s),
            s * s * (3.0 - 2.0 * s),
            s * s * (s - 1.0),
        );
        (0..n)
            .map(|i| {
                h00 * ys[idx][i] + h10 * (tb - ta) * fs[idx][i] + h01 * ys[idx + 1][i]
                    + h11 * (tb - ta) * fs[idx + 1][i]
            })
            .collect()
    };
    let steps = ((t1 - t0) / h).ceil() as usize;
    for s in 0..steps {
        let t = t0 + s as f64 * h;
        let y = ys[s].clone();
        // delayed values at t, t+h/2, t+h; only history available up to t is
        // needed when tau >= h, except when the delayed point equals t (tau == h):
        // handled by using the Hermite data of the interval (t-tau, t).
        let dl = |tt: f64, ts: &[f64], ys: &[Vec<f64>], fs: &[Vec<f64>]| {
            let target = tt - tau;
            if ts.len() < 2 { hist(target) } else { delayed(target, ts, ys, fs) }
        };
        let k1 = f(t, &y, &dl(t, &ts, &ys, &fs));
        let y2: Vec<f64> = (0..n).map(|i| y[i] + 0.5 * h * k1[i]).collect();
        let k2 = f(t + 0.5 * h, &y2, &dl(t + 0.5 * h, &ts, &ys, &fs));
        let y3: Vec<f64> = (0..n).map(|i| y[i] + 0.5 * h * k2[i]).collect();
        let k3 = f(t + 0.5 * h, &y3, &dl(t + 0.5 * h, &ts, &ys, &fs));
        let y4: Vec<f64> = (0..n).map(|i| y[i] + h * k3[i]).collect();
        let k4 = f(t + h, &y4, &dl(t + h, &ts, &ys, &fs));
        let yn: Vec<f64> =
            (0..n).map(|i| y[i] + h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i])).collect();
        let tn = t + h;
        ts.push(tn);
        let fn_ = f(tn, &yn, &dl(tn, &ts, &ys, &fs));
        ys.push(yn);
        fs.push(fn_);
    }
    Ok((ts, ys))
}
