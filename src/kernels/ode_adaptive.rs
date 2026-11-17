//! # Adaptive, stiff, symplectic and boundary-value ODE kernels
//!
//! * [`dopri5`]: Dormand-Prince 5(4) with step-size control, a fourth-order
//!   continuous extension ("dense output") and event detection;
//! * [`rosenbrock23`] / [`rosenbrock23_events`]: the L-stable Rosenbrock pair
//!   of Shampine and Reichelt for stiff problems (finite-difference Jacobian,
//!   cubic Hermite dense output, event detection);
//! * [`leapfrog`] and [`yoshida4`]: symplectic integrators for
//!   `H = p^2/2 + V(q)`;
//! * [`bvp_fd`]: second-order finite-difference BVP solver for
//!   `y'' = f(x, y, y')` with Dirichlet conditions (Newton iteration);
//!   [`shooting`]: general shooting method;
//! * [`dae_index1`]: backward-Euler for semi-explicit index-1 DAEs;
//! * [`dde_rk4`]: constant-delay differential equations by the method of
//!   steps with Hermite history interpolation.
//!
//! # Event location
//!
//! An [`Event`] is a scalar function `g(t, y)`; a sign change is an event. It is
//! located by Brent iteration on the integrator's *own* continuous extension, so
//! the accuracy of an event time is the accuracy of the extension: fourth order
//! for [`dopri5`], the collocation polynomial for Radau IIA, the
//! backward-difference interpolant for BDF, and the cubic Hermite interpolant for
//! [`rosenbrock23`] (third order — its event times are correspondingly coarser
//! than its step tolerance).
//!
//! Every integrator with a continuous extension therefore has an event-aware
//! entry point: [`dopri5`], [`rosenbrock23_events`],
//! [`radau5_events`](crate::kernels::ode_stiff::radau5_events) and
//! [`bdf_events`](crate::kernels::ode_stiff::bdf_events). A `terminal` event ends
//! the integration at the event point, which is recorded as the last sample; a
//! non-terminal one is recorded and the integration continues.
//!
//! The fixed-step solvers in [`crate::kernels::ode`] have no extension, so they
//! locate events by step bisection instead — see
//! [`solve_fixed_events`](crate::kernels::ode::solve_fixed_events).

use super::ode_contract::InitialStepOptions;
use super::ode_contract::initial_events;
use super::ode_contract::initial_step;
use super::ode_contract::validate_problem;
use crate::kernels::dense::Mat;
use crate::kernels::dense::{
    self,
};
use crate::kernels::rootfind;

/// Errors of the ODE kernels.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OdeError {
    /// The embedding host cancelled the computation.
    Cancelled,
    /// Step size underflow or too many steps.
    StepLimit,
    /// The state became non-finite.
    NonFinite,
    /// A nonlinear / linear solve failed.
    SolveFailed,
    /// Invalid arguments.
    Invalid,
}

impl std::fmt::Display for OdeError {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        f.write_str(match self {
            | Self::Cancelled => "computation cancelled",
            | Self::StepLimit => "ODE step size underflow or step budget exhausted",
            | Self::NonFinite => "ODE state or rate became nonfinite",
            | Self::SolveFailed => "ODE nonlinear or linear solve failed",
            | Self::Invalid => "invalid ODE arguments",
        })
    }
}
impl std::error::Error for OdeError {}
impl From<crate::cancellation::Cancelled> for OdeError {
    fn from(_: crate::cancellation::Cancelled) -> Self {
        Self::Cancelled
    }
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
    /// Minimum step magnitude, except a final clipped interval; zero uses only representability.
    pub hmin: f64,
    /// Maximal step size.
    pub hmax: f64,
    /// Maximal number of attempted steps, including rejected steps.
    pub max_steps: usize,
}

impl Default for OdeOptions {
    fn default() -> Self {
        Self {
            rtol: 1e-8,
            atol: 1e-10,
            h0: None,
            hmin: 0.0,
            hmax: f64::INFINITY,
            max_steps: 1_000_000,
        }
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

/// Representation of the continuous extension over one accepted step.
#[derive(Debug, Clone)]
pub(crate) enum DenseKind {
    /// `y0 + h * sum_k d[k] x^(k+1)` with `x = (t - t0) / h`.
    Poly(Vec<Vec<f64>>),
    /// Lagrange interpolation through the nodes `x_j` (in units of the step)
    /// with the given values (the collocation polynomial of Radau IIA).
    Lagrange {
        nodes: Vec<f64>,
        vals: Vec<Vec<f64>>,
    },
    /// Newton backward-difference interpolation `sum_k D_k prod_{j<k} (s+j)/(j+1)`
    /// with `s = x - 1` (the BDF interpolant); `D_0 = y(t0 + h)`.
    Backward(Vec<Vec<f64>>),
    /// Cubic Hermite interpolation from the end value and both end derivatives.
    Hermite {
        y1: Vec<f64>,
        f0: Vec<f64>,
        f1: Vec<f64>,
    },
}

#[derive(Debug, Clone)]
pub(crate) struct DenseStep {
    pub(crate) t0: f64,
    pub(crate) h: f64,
    pub(crate) y0: Vec<f64>,
    pub(crate) kind: DenseKind,
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
    /// Order of the method used for each accepted step (variable-order
    /// integrators only; empty otherwise).
    pub orders: Vec<usize>,
    pub(crate) dense: Vec<DenseStep>,
}

impl OdeSolution {
    /// Creates a solution holding only the initial point.
    pub(crate) fn start(
        t0: f64,
        y0: &[f64],
    ) -> Self {
        Self {
            t: vec![t0],
            y: vec![y0.to_vec()],
            events: Vec::new(),
            nfev: 0,
            rejected: 0,
            orders: Vec::new(),
            dense: Vec::new(),
        }
    }

    /// This solution as a structured result, naming each state component.
    ///
    /// The integrator works with anonymous columns; a model knows what they
    /// are. Passing the names here is what turns `y[3][1]` into
    /// `velocity`, and it is the only place the two views meet.
    ///
    /// A name list shorter than the state is padded with `y{k}`, and one
    /// longer is truncated: a mismatch is a caller's slip, not a reason to
    /// lose the solution.
    #[must_use]
    pub fn view(
        &self,
        components: &[&str],
    ) -> crate::kernel::result::OdeSolutionView {
        let width = self.y.first().map_or(0, Vec::len);
        let flatten_order = (0..width)
            .map(|k| {
                components
                    .get(k)
                    .map_or_else(|| format!("y{k}"), |name| (*name).to_owned())
            })
            .collect();
        let events = self
            .events
            .iter()
            .map(|hit| {
                crate::kernel::result::EventRecord {
                    t: hit.t,
                    component: components.get(hit.index).map(|name| (*name).to_owned()),
                    y: hit.y.clone(),
                }
            })
            .collect();
        crate::kernel::result::OdeSolutionView {
            t: self.t.clone(),
            y: self.y.clone(),
            flatten_order,
            events,
            nfev: self.nfev,
            rejected: self.rejected,
            dimensions: vec![None; width],
        }
    }

    /// Evaluates the continuous extension at `t` (available for
    /// [`dopri5`], [`rosenbrock23`], the Radau and the BDF integrators); `None` outside the integrated range.
    #[must_use]
    pub fn eval(
        &self,
        t: f64,
    ) -> Option<Vec<f64>> {
        let (&first, &last) = (self.t.first()?, self.t.last()?);
        let (lo, hi) = if first <= last {
            (first, last)
        } else {
            (last, first)
        };
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

pub(crate) fn dense_eval(
    s: &DenseStep,
    t: f64,
) -> Vec<f64> {
    let x = (t - s.t0) / s.h;
    let n = s.y0.len();
    match &s.kind {
        | DenseKind::Poly(d) => {
            let (x1, x2, x3, x4) = (x, x * x, x * x * x, x * x * x * x);
            (0..n)
                .map(|i| {
                    s.y0[i] + s.h * (d[0][i] * x1 + d[1][i] * x2 + d[2][i] * x3 + d[3][i] * x4)
                })
                .collect()
        },
        | DenseKind::Lagrange { nodes, vals } => {
            let mut out = vec![0.0; n];
            for (j, vj) in vals.iter().enumerate() {
                let mut l = 1.0;
                for (k, xk) in nodes.iter().enumerate() {
                    if k != j {
                        l *= (x - xk) / (nodes[j] - xk);
                    }
                }
                for i in 0..n {
                    out[i] += l * vj[i];
                }
            }
            out
        },
        | DenseKind::Hermite { y1, f0, f1 } => {
            let (x2, x3) = (x * x, x * x * x);
            let (h00, h10, h01, h11) = (
                2.0 * x3 - 3.0 * x2 + 1.0,
                x3 - 2.0 * x2 + x,
                -2.0 * x3 + 3.0 * x2,
                x3 - x2,
            );
            (0..n)
                .map(|i| h00 * s.y0[i] + h10 * s.h * f0[i] + h01 * y1[i] + h11 * s.h * f1[i])
                .collect()
        },
        | DenseKind::Backward(dd) => {
            let sv = x - 1.0;
            let mut out = dd[0].clone();
            let mut p = 1.0;
            for (k, dk) in dd.iter().enumerate().skip(1) {
                p *= (sv + (k - 1) as f64) / k as f64;
                for i in 0..n {
                    out[i] += p * dk[i];
                }
            }
            out
        },
    }
}

/// The accepted step `[t_old, tn]` (integration direction `dir`) ending at `ynew`.
#[derive(Clone, Copy)]
pub(crate) struct StepWindow<'a> {
    pub(crate) t_old: f64,
    pub(crate) tn: f64,
    pub(crate) ynew: &'a [f64],
    pub(crate) dir: f64,
}

/// Event location tolerance: relative to the step, but never finer than the
/// floating-point spacing at its endpoints.
pub(crate) fn event_time_tolerance(
    t0: f64,
    t1: f64,
) -> f64 {
    let magnitude = t0.abs().max(t1.abs());
    let next = magnitude.next_up();
    let spacing = if next.is_finite() {
        next - magnitude
    } else {
        magnitude - magnitude.next_down()
    };
    ((t1 - t0).abs() * 1e-12)
        .max(spacing)
        .max(f64::from_bits(1))
}

/// Locates the events crossed during an accepted step and updates the
/// solution. Returns `true` when a terminal event ended the integration
/// (the step and the event point are then already recorded).
pub(crate) fn process_events(
    events: &[Event<'_>],
    gprev: &mut [f64],
    sol: &mut OdeSolution,
    step: DenseStep,
    window: &StepWindow<'_>,
) -> Result<bool, OdeError> {
    let StepWindow { t_old, tn, ynew, dir } = *window;
    let mut hits: Vec<EventHit> = Vec::new();
    let mut stop = false;
    for (ei, ev) in events.iter().enumerate() {
        let g1 = (ev.g)(tn, ynew);
        let g0 = gprev[ei];
        if !g0.is_finite() || !g1.is_finite() {
            return Err(OdeError::NonFinite);
        }
        let crossing = (g0 < 0.0 && g1 >= 0.0 && ev.direction >= 0)
            || (g0 > 0.0 && g1 <= 0.0 && ev.direction <= 0);
        if crossing {
            let gf = |tt: f64| (ev.g)(tt, &dense_eval(&step, tt));
            let root = rootfind::brent(gf, t_old, tn, event_time_tolerance(t_old, tn), 200)
                .map_err(|error| {
                    match error {
                        | rootfind::RootError::Invalid => OdeError::NonFinite,
                        | _ => OdeError::SolveFailed,
                    }
                })?;
            hits.push(EventHit {
                index: ei,
                t: root,
                y: dense_eval(&step, root),
            });
            if ev.terminal {
                stop = true;
            }
        }
        gprev[ei] = g1;
    }
    hits.sort_by(|a, b| ((a.t - b.t) * dir).total_cmp(&0.0));
    if stop {
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
        return Ok(true);
    }
    sol.events.extend(hits);
    sol.dense.push(step);
    Ok(false)
}

const DP_A: [[f64; 6]; 6] = [
    [1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0],
    [44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0],
    [
        19372.0 / 6561.0,
        -25360.0 / 2187.0,
        64448.0 / 6561.0,
        -212.0 / 729.0,
        0.0,
        0.0,
    ],
    [
        9017.0 / 3168.0,
        -355.0 / 33.0,
        46732.0 / 5247.0,
        49.0 / 176.0,
        -5103.0 / 18656.0,
        0.0,
    ],
    [
        35.0 / 384.0,
        0.0,
        500.0 / 1113.0,
        125.0 / 192.0,
        -2187.0 / 6784.0,
        11.0 / 84.0,
    ],
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
    [
        1.0,
        -8048581381.0 / 2820520608.0,
        8663915743.0 / 2820520608.0,
        -12715105075.0 / 11282082432.0,
    ],
    [0.0, 0.0, 0.0, 0.0],
    [
        0.0,
        131558114200.0 / 32700410799.0,
        -68118460800.0 / 10900136933.0,
        87487479700.0 / 32700410799.0,
    ],
    [
        0.0,
        -1754552775.0 / 470086768.0,
        14199869525.0 / 1410260304.0,
        -10690763975.0 / 1880347072.0,
    ],
    [
        0.0,
        127303824393.0 / 49829197408.0,
        -318862633887.0 / 49829197408.0,
        701980252875.0 / 199316789632.0,
    ],
    [
        0.0,
        -282668133.0 / 205662961.0,
        2019193451.0 / 616988883.0,
        -1453857185.0 / 822651844.0,
    ],
    [
        0.0,
        40617522.0 / 29380423.0,
        -110615467.0 / 29380423.0,
        69997945.0 / 29380423.0,
    ],
];

pub(crate) fn err_norm(
    err: &[f64],
    y0: &[f64],
    y1: &[f64],
    rtol: f64,
    atol: f64,
) -> f64 {
    super::ode_contract::rms(
        err.iter().zip(y0.iter().zip(y1)).map(|(&error, (&a, &b))| {
            super::ode_contract::scaled_ratio(error, atol + rtol * a.abs().max(b.abs()))
        }),
        err.len(),
    )
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
// Equality detects an unrepresentable step/interval endpoint, not approximate convergence.
#[allow(clippy::float_cmp)]
pub fn dopri5<F: Fn(f64, &[f64], &mut [f64])>(
    f: F,
    t0: f64,
    t1: f64,
    y0: &[f64],
    opts: &OdeOptions,
    events: &[Event<'_>],
) -> Result<OdeSolution, OdeError> {
    dopri5_controlled(f, t0, t1, y0, opts, events, None)
}

// Exact comparisons detect unrepresentable steps and zero local error.
#[allow(clippy::float_cmp)]
pub(crate) fn dopri5_controlled<F: Fn(f64, &[f64], &mut [f64])>(
    f: F,
    t0: f64,
    t1: f64,
    y0: &[f64],
    opts: &OdeOptions,
    events: &[Event<'_>],
    stop: Option<&dyn Fn() -> bool>,
) -> Result<OdeSolution, OdeError> {
    validate_problem(t0, t1, y0, opts, events)?;
    super::ode_contract::abort_if_requested(stop)?;
    let f = super::ode_contract::checked_rhs_with_stop(f, || stop.is_some_and(|stop| stop()));
    let n = y0.len();
    let dir = if t1 > t0 { 1.0 } else { -1.0 };
    let mut t = t0;
    let mut y = y0.to_vec();
    let mut k = vec![vec![0.0; n]; 7];
    f(t, &y, &mut k[0]);
    if k[0].iter().any(|value| !value.is_finite()) {
        return Err(OdeError::NonFinite);
    }
    let mut nfev = 1;
    let mut h = opts
        .h0
        .unwrap_or_else(|| {
            let estimate = initial_step(
                &f,
                t0,
                y0,
                &k[0],
                InitialStepOptions {
                    rtol: opts.rtol,
                    atol: opts.atol,
                    bound: t1,
                    max_step: opts.hmax,
                    order: 5.0,
                },
            );
            nfev += estimate.evaluations;
            estimate.step
        })
        .abs()
        .min(opts.hmax)
        * dir;
    let mut sol = OdeSolution {
        t: vec![t0],
        y: vec![y0.to_vec()],
        events: Vec::new(),
        nfev: 0,
        rejected: 0,
        orders: Vec::new(),
        dense: Vec::new(),
    };
    let mut gprev = initial_events(events, t0, y0)?;
    let mut steps = 0;
    let mut tmp = vec![0.0; n];
    while (t1 - t) * dir > 0.0 {
        super::ode_contract::abort_if_requested(stop)?;
        if steps >= opts.max_steps {
            return Err(OdeError::StepLimit);
        }
        steps += 1;
        if (t + h - t1) * dir > 0.0 {
            h = t1 - t;
        }
        if !h.is_finite() || t + h == t || (h.abs() < opts.hmin && (t1 - t).abs() > opts.hmin) {
            return Err(OdeError::StepLimit);
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
            super::ode_contract::abort_if_requested(stop)?;
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
            if !h.is_finite() || crate::kernels::exact_eq(t + h, t) {
                return Err(OdeError::NonFinite);
            }
            continue;
        }
        let en = err_norm(&err, &y, &ynew, opts.rtol, opts.atol);
        if en.is_nan() {
            return Err(OdeError::NonFinite);
        }
        if en <= 1.0 {
            // dense coefficients
            let d: [Vec<f64>; 4] = std::array::from_fn(|c| {
                (0..n)
                    .map(|i| (0..7).map(|j| DP_P[j][c] * k[j][i]).sum())
                    .collect()
            });
            let step = DenseStep {
                t0: t,
                h,
                y0: y.clone(),
                kind: DenseKind::Poly(d.to_vec()),
            };
            let tn = t + h;
            if process_events(
                events,
                &mut gprev,
                &mut sol,
                step,
                &StepWindow {
                    t_old: t,
                    tn,
                    ynew: &ynew,
                    dir,
                },
            )? {
                sol.nfev = nfev;
                return Ok(sol);
            }
            t = tn;
            y = ynew;
            sol.t.push(t);
            sol.y.push(y.clone());
            // FSAL
            let last = k[6].clone();
            k[0] = last;
            let fac = if en == 0.0 {
                5.0
            } else {
                (0.9 * en.powf(-0.2)).clamp(0.2, 5.0)
            };
            h *= fac;
            if h.abs() > opts.hmax {
                h = opts.hmax * dir;
            }
        } else {
            sol.rejected += 1;
            h *= (0.9 * en.powf(-0.2)).clamp(0.1, 0.9);
            if !h.is_finite() || crate::kernels::exact_eq(t + h, t) {
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
/// differences at every step. Dense output is the cubic Hermite
/// interpolant of the step ([`OdeSolution::eval`]); events are not supported.
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
    rosenbrock23_events(f, t0, t1, y0, opts, &[])
}

/// The same integrator, locating events on its Hermite extension.
///
/// The extension is third order, so an event time is located to the accuracy of
/// the step rather than to the step size; the bracket is refined by Brent
/// iteration to `1e-14`.
///
/// # Errors
/// As [`rosenbrock23`].
// Equality detects an unrepresentable step/interval endpoint, not approximate convergence.
#[allow(clippy::float_cmp)]
pub fn rosenbrock23_events<F: Fn(f64, &[f64], &mut [f64])>(
    f: F,
    t0: f64,
    t1: f64,
    y0: &[f64],
    opts: &OdeOptions,
    events: &[Event<'_>],
) -> Result<OdeSolution, OdeError> {
    rosenbrock23_events_controlled(f, t0, t1, y0, opts, events, None)
}

// Exact comparisons detect unrepresentable steps and zero local error.
#[allow(clippy::float_cmp)]
pub(crate) fn rosenbrock23_events_controlled<F: Fn(f64, &[f64], &mut [f64])>(
    f: F,
    t0: f64,
    t1: f64,
    y0: &[f64],
    opts: &OdeOptions,
    events: &[Event<'_>],
    stop: Option<&dyn Fn() -> bool>,
) -> Result<OdeSolution, OdeError> {
    validate_problem(t0, t1, y0, opts, events)?;
    super::ode_contract::abort_if_requested(stop)?;
    let f = super::ode_contract::checked_rhs_with_stop(f, || stop.is_some_and(|stop| stop()));
    let n = y0.len();
    let dir = if t1 > t0 { 1.0 } else { -1.0 };
    let d = 1.0 / (2.0 + 2.0_f64.sqrt());
    let e32 = 6.0 + 2.0_f64.sqrt();
    let mut t = t0;
    let mut y = y0.to_vec();
    let mut f0 = vec![0.0; n];
    f(t, &y, &mut f0);
    if f0.iter().any(|value| !value.is_finite()) {
        return Err(OdeError::NonFinite);
    }
    let mut nfev = 1;
    let mut h = opts
        .h0
        .unwrap_or_else(|| {
            let estimate = initial_step(
                &f,
                t0,
                y0,
                &f0,
                InitialStepOptions {
                    rtol: opts.rtol,
                    atol: opts.atol,
                    bound: t1,
                    max_step: opts.hmax,
                    order: 2.0,
                },
            );
            nfev += estimate.evaluations;
            estimate.step
        })
        .abs()
        .min(opts.hmax)
        * dir;
    let mut sol = OdeSolution {
        t: vec![t0],
        y: vec![y0.to_vec()],
        events: Vec::new(),
        nfev: 0,
        rejected: 0,
        orders: Vec::new(),
        dense: Vec::new(),
    };
    // The event functions at the start of the interval; a crossing is a sign
    // change from one accepted step's end to the next.
    let mut gprev = initial_events(events, t0, y0)?;
    let mut steps = 0;
    while (t1 - t) * dir > 0.0 {
        super::ode_contract::abort_if_requested(stop)?;
        if steps >= opts.max_steps {
            return Err(OdeError::StepLimit);
        }
        steps += 1;
        if (t + h - t1) * dir > 0.0 {
            h = t1 - t;
        }
        if !h.is_finite() || t + h == t || (h.abs() < opts.hmin && (t1 - t).abs() > opts.hmin) {
            return Err(OdeError::StepLimit);
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
        let probe = super::ode_contract::time_probe(t, t0, t1);
        let dt = probe - t;
        f(probe, &y, &mut fp);
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
            if !h.is_finite() || crate::kernels::exact_eq(t + h, t) {
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
        let err: Vec<f64> = (0..n)
            .map(|i| h / 6.0 * (k1[i] - 2.0 * k2[i] + k3[i]))
            .collect();
        if ynew.iter().any(|v| !v.is_finite()) {
            h *= 0.25;
            sol.rejected += 1;
            if !h.is_finite() || crate::kernels::exact_eq(t + h, t) {
                return Err(OdeError::NonFinite);
            }
            continue;
        }
        let en = err_norm(&err, &y, &ynew, opts.rtol, opts.atol);
        if en.is_nan() {
            return Err(OdeError::NonFinite);
        }
        if en <= 1.0 {
            let step = DenseStep {
                t0: t,
                h,
                y0: y.clone(),
                kind: DenseKind::Hermite {
                    y1: ynew.clone(),
                    f0: f0.clone(),
                    f1: f2.clone(),
                },
            };
            // `process_events` records the step and the events; a terminal
            // event ends the integration where it happened.
            let window = StepWindow {
                t_old: t,
                tn: t + h,
                ynew: &ynew,
                dir,
            };
            if process_events(events, &mut gprev, &mut sol, step, &window)? {
                sol.nfev = nfev;
                return Ok(sol);
            }
            t += h;
            y = ynew;
            f0 = f2;
            sol.t.push(t);
            sol.y.push(y.clone());
            let fac = if en == 0.0 {
                5.0
            } else {
                (0.9 * en.powf(-1.0 / 3.0)).clamp(0.2, 5.0)
            };
            h *= fac;
            if h.abs() > opts.hmax {
                h = opts.hmax * dir;
            }
        } else {
            sol.rejected += 1;
            h *= (0.9 * en.powf(-1.0 / 3.0)).clamp(0.1, 0.9);
            if !h.is_finite() || crate::kernels::exact_eq(t + h, t) {
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

/// Finite-difference BVP solver for `y'' = f(x, y, y')` with Dirichlet ends.
///
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
            if i == 0 {
                ya
            } else if i == n + 1 {
                yb
            } else {
                u[i - 1]
            }
        };
        (1..=n)
            .map(|i| {
                let (ym, y0, yp) = (full(i - 1), full(i), full(i + 1));
                (yp - 2.0 * y0 + ym) / (h * h) - f(x[i], y0, (yp - ym) / (2.0 * h))
            })
            .collect()
    };
    let u0: Vec<f64> = (1..=n).map(|i| guess(x[i])).collect();
    let sol = rootfind::newton_system(resid, None::<&fn(&[f64]) -> Mat>, &u0, 1e-11, 100)
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
            | Ok(s) => residual(s.y.last().map_or(&[][..], Vec::as_slice)),
            | Err(_) => vec![f64::NAN; p.len()],
        }
    };
    rootfind::newton_system(shoot, None::<&fn(&[f64]) -> Mat>, p0, 1e-9, 60)
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
        crate::cancellation::checkpoint()?;
        let tn = t0 + (s + 1) as f64 * h;
        let prev = cur.clone();
        let sys = |u: &[f64]| -> Vec<f64> {
            let (y, z) = u.split_at(ny);
            let fy = f(tn, y, z);
            let mut r: Vec<f64> = (0..ny).map(|i| y[i] - prev[i] - h * fy[i]).collect();
            r.extend(g(tn, y, z));
            r
        };
        let sol = rootfind::newton_system(sys, None::<&fn(&[f64]) -> Mat>, &prev, 1e-12, 50)
            .map_err(|_| OdeError::SolveFailed)?;
        cur = sol.x;
        out.push(cur.clone());
    }
    Ok(out)
}

/// Constant-delay DDE integrator (method of steps with RK4).
///
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
    if !h.is_finite()
        || !tau.is_finite()
        || !t0.is_finite()
        || !t1.is_finite()
        || t1 < t0
        || h <= 0.0
        || tau < h
    {
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
        // The endpoint time is appended before its rate is evaluated. Restrict
        // interpolation to samples whose state and derivative already exist;
        // tau == h lands exactly on the last complete sample.
        let len = ts.len().min(ys.len()).min(fs.len());
        if len < 2 || t >= ts[len - 1] {
            return ys[len - 1].clone();
        }
        let idx = ts[..len]
            .partition_point(|&x| x <= t)
            .saturating_sub(1)
            .min(len - 2);
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
                h00 * ys[idx][i]
                    + h10 * (tb - ta) * fs[idx][i]
                    + h01 * ys[idx + 1][i]
                    + h11 * (tb - ta) * fs[idx + 1][i]
            })
            .collect()
    };
    let steps = ((t1 - t0) / h).ceil() as usize;
    for s in 0..steps {
        crate::cancellation::checkpoint()?;
        let t = t0 + s as f64 * h;
        let h = h.min(t1 - t);
        let y = ys[s].clone();
        // delayed values at t, t+h/2, t+h; only history available up to t is
        // needed when tau >= h, except when the delayed point equals t (tau == h):
        // handled by using the Hermite data of the interval (t-tau, t).
        let dl = |tt: f64, ts: &[f64], ys: &[Vec<f64>], fs: &[Vec<f64>]| {
            let target = tt - tau;
            if ts.len() < 2 {
                hist(target)
            } else {
                delayed(target, ts, ys, fs)
            }
        };
        let k1 = f(t, &y, &dl(t, &ts, &ys, &fs));
        let y2: Vec<f64> = (0..n).map(|i| y[i] + 0.5 * h * k1[i]).collect();
        let k2 = f(t + 0.5 * h, &y2, &dl(t + 0.5 * h, &ts, &ys, &fs));
        let y3: Vec<f64> = (0..n).map(|i| y[i] + 0.5 * h * k2[i]).collect();
        let k3 = f(t + 0.5 * h, &y3, &dl(t + 0.5 * h, &ts, &ys, &fs));
        let y4: Vec<f64> = (0..n).map(|i| y[i] + h * k3[i]).collect();
        let k4 = f(t + h, &y4, &dl(t + h, &ts, &ys, &fs));
        let yn: Vec<f64> = (0..n)
            .map(|i| y[i] + h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]))
            .collect();
        let tn = t + h;
        ts.push(tn);
        let fn_ = f(tn, &yn, &dl(tn, &ts, &ys, &fs));
        ys.push(yn);
        fs.push(fn_);
    }
    Ok((ts, ys))
}

/// Fixed-step Dormand-Prince 5(4) integrator with dense output and event
/// location.
///
/// Integrates `y' = f(t, y)` from `t0` to `t1` in `steps` equal steps with no
/// step-size rejection. Populates the fourth-order continuous extension
/// ([`OdeSolution::eval`]) and locates events by Brent iteration.
///
/// # Errors
/// [`OdeError::Invalid`] if `steps == 0` or inputs are invalid,
/// [`OdeError::NonFinite`] on non-finite values, or [`OdeError::StepLimit`] on
/// underflow.
#[allow(clippy::float_cmp)]
pub fn dopri5_fixed<F: Fn(f64, &[f64], &mut [f64])>(
    f: F,
    t0: f64,
    t1: f64,
    y0: &[f64],
    steps: usize,
    events: &[Event<'_>],
) -> Result<OdeSolution, OdeError> {
    if steps == 0 {
        return Err(OdeError::Invalid);
    }
    let default_opts = OdeOptions::default();
    validate_problem(t0, t1, y0, &default_opts, events)?;
    if crate::kernels::exact_eq(t0, t1) {
        return Ok(OdeSolution::start(t0, y0));
    }
    let f = super::ode_contract::checked_rhs(f);
    let n = y0.len();
    let dir = if t1 > t0 { 1.0 } else { -1.0 };
    let mut t = t0;
    let mut y = y0.to_vec();
    let mut k = vec![vec![0.0; n]; 7];
    let Some(k0) = k.first_mut() else {
        return Err(OdeError::Invalid);
    };
    f(t, &y, k0);
    if k0.iter().any(|value| !value.is_finite()) {
        return Err(OdeError::NonFinite);
    }
    let mut nfev = 1usize;
    let mut sol = OdeSolution {
        t: vec![t0],
        y: vec![y0.to_vec()],
        events: Vec::new(),
        nfev: 0,
        rejected: 0,
        orders: Vec::new(),
        dense: Vec::new(),
    };
    let mut gprev = initial_events(events, t0, y0)?;
    let mut tmp = vec![0.0; n];
    for step_idx in 1..=steps {
        let (tn, h) = super::ode_contract::fixed_step(t0, t1, step_idx, steps, t)?;
        for s in 0..6 {
            let row = DP_A.get(s).ok_or(OdeError::Invalid)?;
            for (i, tmp_i) in tmp.iter_mut().enumerate() {
                let mut acc = 0.0;
                for j in 0..=s {
                    let a_sj = row.get(j).copied().unwrap_or(0.0);
                    if a_sj != 0.0 {
                        let k_ji = k.get(j).and_then(|kj| kj.get(i)).copied().unwrap_or(0.0);
                        acc += a_sj * k_ji;
                    }
                }
                *tmp_i = y.get(i).copied().unwrap_or(0.0) + h * acc;
            }
            let c_s = DP_C.get(s).copied().unwrap_or(0.0);
            let (_, tail) = k.split_at_mut(s.saturating_add(1));
            let Some(next_k) = tail.first_mut() else {
                return Err(OdeError::Invalid);
            };
            f(t + c_s * h, &tmp, next_k);
            nfev = nfev.saturating_add(1);
        }
        let ynew = tmp.clone();
        if ynew.iter().any(|v| !v.is_finite()) {
            return Err(OdeError::NonFinite);
        }
        let d: [Vec<f64>; 4] = std::array::from_fn(|c| {
            (0..n)
                .map(|i| {
                    (0..7)
                        .map(|j| {
                            let p_jc = DP_P
                                .get(j)
                                .and_then(|row| row.get(c))
                                .copied()
                                .unwrap_or(0.0);
                            let k_ji = k.get(j).and_then(|kj| kj.get(i)).copied().unwrap_or(0.0);
                            p_jc * k_ji
                        })
                        .sum()
                })
                .collect()
        });
        let step_dense = DenseStep {
            t0: t,
            h,
            y0: y.clone(),
            kind: DenseKind::Poly(d.to_vec()),
        };
        if process_events(
            events,
            &mut gprev,
            &mut sol,
            step_dense,
            &StepWindow {
                t_old: t,
                tn,
                ynew: &ynew,
                dir,
            },
        )? {
            sol.nfev = nfev;
            return Ok(sol);
        }
        sol.t.push(tn);
        sol.y.push(ynew.clone());
        t = tn;
        y = ynew;
        k.swap(0, 6);
    }
    sol.nfev = nfev;
    Ok(sol)
}

/// Fixed-step Rosenbrock 2(3) integrator without event detection.
///
/// # Errors
/// [`OdeError::Invalid`], [`OdeError::NonFinite`], [`OdeError::SolveFailed`]
/// or [`OdeError::StepLimit`].
pub fn rosenbrock23_fixed<F: Fn(f64, &[f64], &mut [f64])>(
    f: F,
    t0: f64,
    t1: f64,
    y0: &[f64],
    steps: usize,
) -> Result<OdeSolution, OdeError> {
    rosenbrock23_fixed_events(f, t0, t1, y0, steps, &[])
}

/// Fixed-step Rosenbrock 2(3) integrator with cubic Hermite dense output and
/// event location.
///
/// # Errors
/// As [`rosenbrock23_fixed`].
#[allow(clippy::float_cmp)]
pub fn rosenbrock23_fixed_events<F: Fn(f64, &[f64], &mut [f64])>(
    f: F,
    t0: f64,
    t1: f64,
    y0: &[f64],
    steps: usize,
    events: &[Event<'_>],
) -> Result<OdeSolution, OdeError> {
    if steps == 0 {
        return Err(OdeError::Invalid);
    }
    let default_opts = OdeOptions::default();
    validate_problem(t0, t1, y0, &default_opts, events)?;
    if crate::kernels::exact_eq(t0, t1) {
        return Ok(OdeSolution::start(t0, y0));
    }
    let f = super::ode_contract::checked_rhs(f);
    let n = y0.len();
    let dir = if t1 > t0 { 1.0 } else { -1.0 };
    let d = 1.0 / (2.0 + 2.0_f64.sqrt());
    let mut t = t0;
    let mut y = y0.to_vec();
    let mut f0 = vec![0.0; n];
    f(t, &y, &mut f0);
    if f0.iter().any(|value| !value.is_finite()) {
        return Err(OdeError::NonFinite);
    }
    let mut nfev = 1usize;
    let mut sol = OdeSolution {
        t: vec![t0],
        y: vec![y0.to_vec()],
        events: Vec::new(),
        nfev: 0,
        rejected: 0,
        orders: Vec::new(),
        dense: Vec::new(),
    };
    let mut gprev = initial_events(events, t0, y0)?;
    for step_idx in 1..=steps {
        let (tn, h) = super::ode_contract::fixed_step(t0, t1, step_idx, steps, t)?;
        let mut jac = Mat::zeros(n, n);
        let mut yp = y.clone();
        let mut fp = vec![0.0; n];
        for j in 0..n {
            let yj = y.get(j).copied().unwrap_or(0.0);
            let dj = 1.5e-8 * yj.abs().max(1e-5);
            if let Some(slot) = yp.get_mut(j) {
                *slot = yj + dj;
            }
            f(t, &yp, &mut fp);
            if let Some(slot) = yp.get_mut(j) {
                *slot = yj;
            }
            for i in 0..n {
                let df = fp.get(i).copied().unwrap_or(0.0) - f0.get(i).copied().unwrap_or(0.0);
                jac.set(i, j, df / dj);
            }
        }
        let probe = super::ode_contract::time_probe(t, t0, t1);
        let dt = probe - t;
        f(probe, &y, &mut fp);
        let tder: Vec<f64> = fp
            .iter()
            .zip(&f0)
            .map(|(fpi, f0i)| (*fpi - *f0i) / dt)
            .collect();
        nfev = nfev.saturating_add(n.saturating_add(1));
        let mut w = Mat::identity(n);
        for i in 0..n {
            for j in 0..n {
                let v = w.at(i, j) - h * d * jac.at(i, j);
                w.set(i, j, v);
            }
        }
        let lu = dense::lu_factor(&w).map_err(|_| OdeError::SolveFailed)?;
        let rhs1: Vec<f64> = f0
            .iter()
            .zip(&tder)
            .map(|(f0i, tdi)| *f0i + h * d * *tdi)
            .collect();
        let k1 = lu.solve(&rhs1);
        let ymid: Vec<f64> = y
            .iter()
            .zip(&k1)
            .map(|(yi, k1i)| *yi + 0.5 * h * *k1i)
            .collect();
        let mut f1 = vec![0.0; n];
        f(t + 0.5 * h, &ymid, &mut f1);
        let rhs2: Vec<f64> = f1
            .iter()
            .zip(&k1)
            .zip(&tder)
            .map(|((f1i, k1i), tdi)| *f1i - *k1i + h * d * *tdi)
            .collect();
        let k2 = lu
            .solve(&rhs2)
            .iter()
            .zip(&k1)
            .map(|(a, b)| *a + *b)
            .collect::<Vec<_>>();
        let ynew: Vec<f64> = y
            .iter()
            .zip(&k1)
            .zip(&k2)
            .map(|((yi, k1i), k2i)| *yi + h * (*k1i + *k2i) / 2.0)
            .collect();
        if ynew.iter().any(|v| !v.is_finite()) {
            return Err(OdeError::NonFinite);
        }
        let mut f2 = vec![0.0; n];
        f(tn, &ynew, &mut f2);
        nfev = nfev.saturating_add(2);
        if f2.iter().any(|value| !value.is_finite()) {
            return Err(OdeError::NonFinite);
        }
        let step_dense = DenseStep {
            t0: t,
            h,
            y0: y.clone(),
            kind: DenseKind::Hermite {
                y1: ynew.clone(),
                f0: f0.clone(),
                f1: f2.clone(),
            },
        };
        if process_events(
            events,
            &mut gprev,
            &mut sol,
            step_dense,
            &StepWindow {
                t_old: t,
                tn,
                ynew: &ynew,
                dir,
            },
        )? {
            sol.nfev = nfev;
            return Ok(sol);
        }
        sol.t.push(tn);
        sol.y.push(ynew.clone());
        t = tn;
        y = ynew;
        f0 = f2;
    }
    sol.nfev = nfev;
    Ok(sol)
}

/// Adaptive velocity-Verlet (leapfrog) via Richardson step-doubling error control.
///
/// Integrates `q'' = acc(q)` over `[t0, t1]`, returning an [`OdeSolution`] whose
/// state vector at each node is `[q..., p...]`.
///
/// # Errors
/// [`OdeError::Invalid`] for mismatched or empty states, [`OdeError::NonFinite`]
/// on non-finite values, or [`OdeError::StepLimit`] if `max_steps` or `hmin` is
/// exceeded.
pub fn leapfrog_adaptive<A: Fn(&[f64]) -> Vec<f64>>(
    acc: A,
    q0: &[f64],
    p0: &[f64],
    t0: f64,
    t1: f64,
    opts: &OdeOptions,
) -> Result<OdeSolution, OdeError> {
    symplectic_doubling_adaptive(
        |q, p, dt| leapfrog(&acc, q, p, dt, 1),
        2.0,
        q0,
        p0,
        t0,
        t1,
        opts,
    )
}

/// Adaptive fourth-order Yoshida composition via Richardson step-doubling error
/// control.
///
/// Integrates `q'' = acc(q)` over `[t0, t1]`, returning an [`OdeSolution`] whose
/// state vector at each node is `[q..., p...]`.
///
/// # Errors
/// As [`leapfrog_adaptive`].
pub fn yoshida4_adaptive<A: Fn(&[f64]) -> Vec<f64>>(
    acc: A,
    q0: &[f64],
    p0: &[f64],
    t0: f64,
    t1: f64,
    opts: &OdeOptions,
) -> Result<OdeSolution, OdeError> {
    symplectic_doubling_adaptive(
        |q, p, dt| yoshida4(&acc, q, p, dt, 1),
        4.0,
        q0,
        p0,
        t0,
        t1,
        opts,
    )
}

fn symplectic_doubling_adaptive<S>(
    step_fn: S,
    order: f64,
    q0: &[f64],
    p0: &[f64],
    t0: f64,
    t1: f64,
    opts: &OdeOptions,
) -> Result<OdeSolution, OdeError>
where
    S: Fn(&[f64], &[f64], f64) -> (Vec<f64>, Vec<f64>),
{
    if q0.is_empty()
        || q0.len() != p0.len()
        || !t0.is_finite()
        || !t1.is_finite()
        || !opts.rtol.is_finite()
        || opts.rtol < 0.0
        || !opts.atol.is_finite()
        || opts.atol <= 0.0
        || opts.max_steps == 0
        || q0.iter().chain(p0).any(|v| !v.is_finite())
    {
        return Err(OdeError::Invalid);
    }
    let y_init: Vec<f64> = q0.iter().chain(p0).copied().collect();
    if crate::kernels::exact_eq(t0, t1) {
        return Ok(OdeSolution::start(t0, &y_init));
    }
    let dir = if t1 > t0 { 1.0 } else { -1.0 };
    let span = (t1 - t0).abs();
    let mut h = opts
        .h0
        .unwrap_or_else(|| (span * 0.01).max(1e-6))
        .abs()
        .clamp(opts.hmin.max(1e-14), opts.hmax)
        * dir;
    let mut t = t0;
    let mut q = q0.to_vec();
    let mut p = p0.to_vec();
    let mut sol = OdeSolution::start(t0, &y_init);
    let denom_err = 2.0_f64.powf(order) - 1.0;
    let mut attempts = 0usize;
    while (t1 - t) * dir > 0.0 {
        crate::cancellation::checkpoint()?;
        if attempts >= opts.max_steps {
            return Err(OdeError::StepLimit);
        }
        attempts = attempts.saturating_add(1);
        if (t + h - t1) * dir > 0.0 {
            h = t1 - t;
        }
        if !h.is_finite()
            || crate::kernels::exact_eq(t + h, t)
            || (h.abs() < opts.hmin && (t1 - t).abs() > opts.hmin)
        {
            return Err(OdeError::StepLimit);
        }
        let (qf, pf) = step_fn(&q, &p, h);
        let (qm, pm) = step_fn(&q, &p, 0.5 * h);
        let (qh, ph) = step_fn(&qm, &pm, 0.5 * h);
        sol.nfev = sol.nfev.saturating_add(3);
        if qh
            .iter()
            .chain(&ph)
            .chain(&qf)
            .chain(&pf)
            .any(|v| !v.is_finite())
        {
            h *= 0.5;
            sol.rejected = sol.rejected.saturating_add(1);
            continue;
        }
        let y_cur: Vec<f64> = q.iter().chain(&p).copied().collect();
        let y_half: Vec<f64> = qh.iter().chain(&ph).copied().collect();
        let y_full: Vec<f64> = qf.iter().chain(&pf).copied().collect();
        let err_vec: Vec<f64> = y_half
            .iter()
            .zip(&y_full)
            .map(|(yh, yf)| (*yh - *yf).abs() / denom_err)
            .collect();
        let en = err_norm(&err_vec, &y_cur, &y_half, opts.rtol, opts.atol);
        if en.is_nan() {
            return Err(OdeError::NonFinite);
        }
        if en <= 1.0 {
            t += h;
            q = qh;
            p = ph;
            sol.t.push(t);
            sol.y.push(y_half);
        } else {
            sol.rejected = sol.rejected.saturating_add(1);
        }
        let factor = if en == 0.0 {
            5.0
        } else {
            (0.9 * en.powf(-1.0 / (order + 1.0))).clamp(0.2, 5.0)
        };
        let next_h = (h.abs() * factor).min(opts.hmax);
        h = next_h * dir;
    }
    Ok(sol)
}

/// Unified velocity-Verlet integrator supporting both fixed and adaptive step
/// sizes via [`OdeStepMode`](super::ode::OdeStepMode).
///
/// # Errors
/// As [`leapfrog_adaptive`].
pub fn leapfrog_with<A: Fn(&[f64]) -> Vec<f64>>(
    acc: A,
    q0: &[f64],
    p0: &[f64],
    t0: f64,
    t1: f64,
    step_mode: super::ode::OdeStepMode,
) -> Result<OdeSolution, OdeError> {
    match step_mode {
        | super::ode::OdeStepMode::Fixed { dt } => {
            symplectic_fixed(|q, p, h| leapfrog(&acc, q, p, h, 1), q0, p0, t0, t1, dt)
        },
        | super::ode::OdeStepMode::Adaptive { .. } => {
            let opts = step_mode.to_ode_options()?;
            leapfrog_adaptive(acc, q0, p0, t0, t1, &opts)
        },
    }
}

/// Unified fourth-order Yoshida integrator supporting both fixed and adaptive
/// step sizes via [`OdeStepMode`](super::ode::OdeStepMode).
///
/// # Errors
/// As [`yoshida4_adaptive`].
pub fn yoshida4_with<A: Fn(&[f64]) -> Vec<f64>>(
    acc: A,
    q0: &[f64],
    p0: &[f64],
    t0: f64,
    t1: f64,
    step_mode: super::ode::OdeStepMode,
) -> Result<OdeSolution, OdeError> {
    match step_mode {
        | super::ode::OdeStepMode::Fixed { dt } => {
            symplectic_fixed(|q, p, h| yoshida4(&acc, q, p, h, 1), q0, p0, t0, t1, dt)
        },
        | super::ode::OdeStepMode::Adaptive { .. } => {
            let opts = step_mode.to_ode_options()?;
            yoshida4_adaptive(acc, q0, p0, t0, t1, &opts)
        },
    }
}

fn symplectic_fixed<S>(
    step_fn: S,
    q0: &[f64],
    p0: &[f64],
    t0: f64,
    t1: f64,
    dt: f64,
) -> Result<OdeSolution, OdeError>
where
    S: Fn(&[f64], &[f64], f64) -> (Vec<f64>, Vec<f64>),
{
    if q0.is_empty()
        || q0.len() != p0.len()
        || !t0.is_finite()
        || !t1.is_finite()
        || !dt.is_finite()
        || dt <= 0.0
        || q0.iter().chain(p0).any(|v| !v.is_finite())
    {
        return Err(OdeError::Invalid);
    }
    let y_init: Vec<f64> = q0.iter().chain(p0).copied().collect();
    if crate::kernels::exact_eq(t0, t1) {
        return Ok(OdeSolution::start(t0, &y_init));
    }
    let steps = ((t1 - t0).abs() / dt).ceil().max(1.0) as usize;
    let mut t = t0;
    let mut q = q0.to_vec();
    let mut p = p0.to_vec();
    let mut sol = OdeSolution::start(t0, &y_init);
    for idx in 1..=steps {
        crate::cancellation::checkpoint()?;
        let (tn, h) = super::ode_contract::fixed_step(t0, t1, idx, steps, t)?;
        let (qn, pn) = step_fn(&q, &p, h);
        if qn.iter().chain(&pn).any(|v| !v.is_finite()) {
            return Err(OdeError::NonFinite);
        }
        q = qn;
        p = pn;
        t = tn;
        sol.t.push(t);
        sol.y.push(q.iter().chain(&p).copied().collect());
        sol.nfev = sol.nfev.saturating_add(1);
    }
    Ok(sol)
}

/// Adaptive semi-explicit index-1 DAE solver (`y' = f(t, y, z)`,
/// `0 = g(t, y, z)`) using backward Euler with Richardson step-doubling.
///
/// # Errors
/// [`OdeError::Invalid`], [`OdeError::SolveFailed`], [`OdeError::NonFinite`],
/// or [`OdeError::StepLimit`].
#[allow(clippy::float_cmp)]
pub fn dae_index1_adaptive<F, G>(
    f: F,
    g: G,
    y0: &[f64],
    z0: &[f64],
    t0: f64,
    t1: f64,
    opts: &OdeOptions,
) -> Result<OdeSolution, OdeError>
where
    F: Fn(f64, &[f64], &[f64]) -> Vec<f64>,
    G: Fn(f64, &[f64], &[f64]) -> Vec<f64>,
{
    if y0.is_empty()
        || !t0.is_finite()
        || !t1.is_finite()
        || !opts.rtol.is_finite()
        || opts.rtol < 0.0
        || !opts.atol.is_finite()
        || opts.atol <= 0.0
        || opts.max_steps == 0
        || y0.iter().chain(z0).any(|v| !v.is_finite())
    {
        return Err(OdeError::Invalid);
    }
    let ny = y0.len();
    let init: Vec<f64> = y0.iter().chain(z0).copied().collect();
    if crate::kernels::exact_eq(t0, t1) {
        return Ok(OdeSolution::start(t0, &init));
    }
    let dir = if t1 > t0 { 1.0 } else { -1.0 };
    let mut h = opts
        .h0
        .unwrap_or_else(|| ((t1 - t0).abs() * 0.01).max(1e-6))
        .abs()
        .clamp(opts.hmin.max(1e-14), opts.hmax)
        * dir;
    let mut t = t0;
    let mut cur = init;
    let mut sol = OdeSolution::start(t0, &cur);
    let mut attempts = 0usize;
    while (t1 - t) * dir > 0.0 {
        crate::cancellation::checkpoint()?;
        if attempts >= opts.max_steps {
            return Err(OdeError::StepLimit);
        }
        attempts = attempts.saturating_add(1);
        if (t + h - t1) * dir > 0.0 {
            h = t1 - t;
        }
        if !h.is_finite()
            || crate::kernels::exact_eq(t + h, t)
            || (h.abs() < opts.hmin && (t1 - t).abs() > opts.hmin)
        {
            return Err(OdeError::StepLimit);
        }
        let (y_cur, z_cur) = cur.split_at(ny);
        let Ok(full_traj) = dae_index1(&f, &g, y_cur, z_cur, t, h, 1) else {
            h *= 0.5;
            sol.rejected = sol.rejected.saturating_add(1);
            continue;
        };
        let Ok(half_traj) = dae_index1(&f, &g, y_cur, z_cur, t, 0.5 * h, 2) else {
            h *= 0.5;
            sol.rejected = sol.rejected.saturating_add(1);
            continue;
        };
        let Some(full) = full_traj.last() else {
            return Err(OdeError::SolveFailed);
        };
        let Some(half) = half_traj.last() else {
            return Err(OdeError::SolveFailed);
        };
        sol.nfev = sol.nfev.saturating_add(3);
        let err_vec: Vec<f64> = half
            .iter()
            .zip(full)
            .map(|(a, b)| (*a - *b).abs())
            .collect();
        let en = err_norm(&err_vec, &cur, half, opts.rtol, opts.atol);
        if en.is_nan() {
            return Err(OdeError::NonFinite);
        }
        if en <= 1.0 {
            t += h;
            cur = half.clone();
            sol.t.push(t);
            sol.y.push(cur.clone());
        } else {
            sol.rejected = sol.rejected.saturating_add(1);
        }
        let factor = if en == 0.0 {
            5.0
        } else {
            (0.9 * en.powf(-0.5)).clamp(0.2, 5.0)
        };
        h = (h.abs() * factor).min(opts.hmax) * dir;
    }
    Ok(sol)
}

/// Unified semi-explicit index-1 DAE solver supporting both fixed and adaptive
/// step sizes via [`OdeStepMode`](super::ode::OdeStepMode).
///
/// # Errors
/// As [`dae_index1_adaptive`].
pub fn dae_index1_with<F, G>(
    f: F,
    g: G,
    y0: &[f64],
    z0: &[f64],
    t0: f64,
    t1: f64,
    step_mode: super::ode::OdeStepMode,
) -> Result<OdeSolution, OdeError>
where
    F: Fn(f64, &[f64], &[f64]) -> Vec<f64>,
    G: Fn(f64, &[f64], &[f64]) -> Vec<f64>,
{
    match step_mode {
        | super::ode::OdeStepMode::Fixed { dt } => {
            if !dt.is_finite() || dt <= 0.0 || !t0.is_finite() || !t1.is_finite() {
                return Err(OdeError::Invalid);
            }
            if crate::kernels::exact_eq(t0, t1) {
                let init: Vec<f64> = y0.iter().chain(z0).copied().collect();
                return Ok(OdeSolution::start(t0, &init));
            }
            let steps = ((t1 - t0).abs() / dt).ceil().max(1.0) as usize;
            let h = (t1 - t0) / steps as f64;
            let states = dae_index1(f, g, y0, z0, t0, h, steps)?;
            let ts: Vec<f64> = (0..=steps)
                .map(|i| {
                    if i == steps {
                        t1
                    } else {
                        t0 + i as f64 * h
                    }
                })
                .collect();
            Ok(OdeSolution {
                t: ts,
                y: states,
                events: Vec::new(),
                nfev: steps,
                rejected: 0,
                orders: Vec::new(),
                dense: Vec::new(),
            })
        },
        | super::ode::OdeStepMode::Adaptive { .. } => {
            let opts = step_mode.to_ode_options()?;
            dae_index1_adaptive(f, g, y0, z0, t0, t1, &opts)
        },
    }
}

/// Adaptive constant-delay DDE integrator using RK4 with Richardson
/// step-doubling error control and cubic Hermite history interpolation.
///
/// # Errors
/// [`OdeError::Invalid`], [`OdeError::NonFinite`], or [`OdeError::StepLimit`].
#[allow(clippy::float_cmp)]
pub fn dde_rk4_adaptive<F, H>(
    f: F,
    hist: H,
    tau: f64,
    t0: f64,
    t1: f64,
    opts: &OdeOptions,
) -> Result<OdeSolution, OdeError>
where
    F: Fn(f64, &[f64], &[f64]) -> Vec<f64>,
    H: Fn(f64) -> Vec<f64>,
{
    if !tau.is_finite()
        || tau <= 0.0
        || !t0.is_finite()
        || !t1.is_finite()
        || t1 < t0
        || !opts.rtol.is_finite()
        || opts.rtol < 0.0
        || !opts.atol.is_finite()
        || opts.atol <= 0.0
        || opts.max_steps == 0
    {
        return Err(OdeError::Invalid);
    }
    let y0 = hist(t0);
    if y0.is_empty() || y0.iter().any(|v| !v.is_finite()) {
        return Err(OdeError::Invalid);
    }
    if crate::kernels::exact_eq(t0, t1) {
        return Ok(OdeSolution::start(t0, &y0));
    }
    let n = y0.len();
    let f0 = f(t0, &y0, &hist(t0 - tau));
    let mut ts = vec![t0];
    let mut ys = vec![y0.clone()];
    let mut fs = vec![f0];
    let mut h = opts
        .h0
        .unwrap_or_else(|| ((t1 - t0) * 0.02).min(tau * 0.5).max(1e-6))
        .abs()
        .clamp(opts.hmin.max(1e-14), opts.hmax.min(tau));
    let mut t = t0;
    let mut y = y0;
    let mut rejected = 0usize;
    let mut nfev = 1usize;
    let mut attempts = 0usize;
    let eval_delay = |target: f64, ts: &[f64], ys: &[Vec<f64>], fs: &[Vec<f64>]| -> Vec<f64> {
        if target <= t0 {
            return hist(target);
        }
        let len = ts.len().min(ys.len()).min(fs.len());
        let Some(last_t) = len.checked_sub(1).and_then(|i| ts.get(i)).copied() else {
            return hist(target);
        };
        if len < 2 || target >= last_t {
            return len
                .checked_sub(1)
                .and_then(|i| ys.get(i))
                .cloned()
                .unwrap_or_else(|| hist(target));
        }
        let idx = ts
            .get(..len)
            .map_or(0, |slice| slice.partition_point(|&x| x <= target))
            .saturating_sub(1)
            .min(len.saturating_sub(2));
        let ta = ts.get(idx).copied().unwrap_or(t0);
        let tb = ts.get(idx.saturating_add(1)).copied().unwrap_or(last_t);
        let dt_seg = tb - ta;
        let s = (target - ta) / dt_seg;
        let (h00, h10, h01, h11) = (
            (1.0 + 2.0 * s) * (1.0 - s) * (1.0 - s),
            s * (1.0 - s) * (1.0 - s),
            s * s * (3.0 - 2.0 * s),
            s * s * (s - 1.0),
        );
        (0..n)
            .map(|i| {
                let ya = ys.get(idx).and_then(|v| v.get(i)).copied().unwrap_or(0.0);
                let fa = fs.get(idx).and_then(|v| v.get(i)).copied().unwrap_or(0.0);
                let yb = ys
                    .get(idx.saturating_add(1))
                    .and_then(|v| v.get(i))
                    .copied()
                    .unwrap_or(0.0);
                let fb = fs
                    .get(idx.saturating_add(1))
                    .and_then(|v| v.get(i))
                    .copied()
                    .unwrap_or(0.0);
                h00 * ya + h10 * dt_seg * fa + h01 * yb + h11 * dt_seg * fb
            })
            .collect()
    };
    let rk4_substep =
        |t_cur: f64, y_cur: &[f64], dt: f64, ts: &[f64], ys: &[Vec<f64>], fs: &[Vec<f64>]| {
            let k1 = f(t_cur, y_cur, &eval_delay(t_cur - tau, ts, ys, fs));
            let y2: Vec<f64> = y_cur
                .iter()
                .zip(&k1)
                .map(|(yi, ki)| *yi + 0.5 * dt * *ki)
                .collect();
            let k2 = f(
                t_cur + 0.5 * dt,
                &y2,
                &eval_delay(t_cur + 0.5 * dt - tau, ts, ys, fs),
            );
            let y3: Vec<f64> = y_cur
                .iter()
                .zip(&k2)
                .map(|(yi, ki)| *yi + 0.5 * dt * *ki)
                .collect();
            let k3 = f(
                t_cur + 0.5 * dt,
                &y3,
                &eval_delay(t_cur + 0.5 * dt - tau, ts, ys, fs),
            );
            let y4: Vec<f64> = y_cur
                .iter()
                .zip(&k3)
                .map(|(yi, ki)| *yi + dt * *ki)
                .collect();
            let k4 = f(t_cur + dt, &y4, &eval_delay(t_cur + dt - tau, ts, ys, fs));
            (0..n)
                .map(|i| {
                    let yi = y_cur.get(i).copied().unwrap_or(0.0);
                    let k1i = k1.get(i).copied().unwrap_or(0.0);
                    let k2i = k2.get(i).copied().unwrap_or(0.0);
                    let k3i = k3.get(i).copied().unwrap_or(0.0);
                    let k4i = k4.get(i).copied().unwrap_or(0.0);
                    yi + dt / 6.0 * (k1i + 2.0 * k2i + 2.0 * k3i + k4i)
                })
                .collect::<Vec<f64>>()
        };
    while t < t1 {
        crate::cancellation::checkpoint()?;
        if attempts >= opts.max_steps {
            return Err(OdeError::StepLimit);
        }
        attempts = attempts.saturating_add(1);
        h = h.min(tau).min(t1 - t);
        if !h.is_finite()
            || crate::kernels::exact_eq(t + h, t)
            || (h < opts.hmin && (t1 - t) > opts.hmin)
        {
            return Err(OdeError::StepLimit);
        }
        let y_full = rk4_substep(t, &y, h, &ts, &ys, &fs);
        let y_mid = rk4_substep(t, &y, 0.5 * h, &ts, &ys, &fs);
        let y_half = rk4_substep(t + 0.5 * h, &y_mid, 0.5 * h, &ts, &ys, &fs);
        nfev = nfev.saturating_add(12);
        if y_half.iter().chain(&y_full).any(|v| !v.is_finite()) {
            h *= 0.5;
            rejected = rejected.saturating_add(1);
            continue;
        }
        let err_vec: Vec<f64> = y_half
            .iter()
            .zip(&y_full)
            .map(|(yh, yf)| (*yh - *yf).abs() / 15.0)
            .collect();
        let en = err_norm(&err_vec, &y, &y_half, opts.rtol, opts.atol);
        if en.is_nan() {
            return Err(OdeError::NonFinite);
        }
        if en <= 1.0 {
            t += h;
            y = y_half;
            let fn_ = f(t, &y, &eval_delay(t - tau, &ts, &ys, &fs));
            nfev = nfev.saturating_add(1);
            ts.push(t);
            ys.push(y.clone());
            fs.push(fn_);
        } else {
            rejected = rejected.saturating_add(1);
        }
        let factor = if en == 0.0 {
            5.0
        } else {
            (0.9 * en.powf(-0.2)).clamp(0.2, 5.0)
        };
        h = (h * factor).min(opts.hmax).min(tau);
    }
    Ok(OdeSolution {
        t: ts,
        y: ys,
        events: Vec::new(),
        nfev,
        rejected,
        orders: Vec::new(),
        dense: Vec::new(),
    })
}

/// Unified constant-delay DDE integrator supporting both fixed and adaptive
/// step sizes via [`OdeStepMode`](super::ode::OdeStepMode).
///
/// # Errors
/// As [`dde_rk4_adaptive`] and [`dde_rk4`].
pub fn dde_rk4_with<F, H>(
    f: F,
    hist: H,
    tau: f64,
    t0: f64,
    t1: f64,
    step_mode: super::ode::OdeStepMode,
) -> Result<OdeSolution, OdeError>
where
    F: Fn(f64, &[f64], &[f64]) -> Vec<f64>,
    H: Fn(f64) -> Vec<f64>,
{
    match step_mode {
        | super::ode::OdeStepMode::Fixed { dt } => {
            let (ts, ys) = dde_rk4(f, hist, tau, t0, t1, dt)?;
            let nfev = ts.len().saturating_sub(1).saturating_mul(5);
            Ok(OdeSolution {
                t: ts,
                y: ys,
                events: Vec::new(),
                nfev,
                rejected: 0,
                orders: Vec::new(),
                dense: Vec::new(),
            })
        },
        | super::ode::OdeStepMode::Adaptive { .. } => {
            let opts = step_mode.to_ode_options()?;
            dde_rk4_adaptive(f, hist, tau, t0, t1, &opts)
        },
    }
}
