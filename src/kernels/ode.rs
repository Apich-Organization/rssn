//! # Initial value problems
//!
//! Integrators for first-order systems `y' = f(t, y)`. The right-hand side
//! is a closure writing the derivative into a caller-provided buffer, so the
//! steppers allocate nothing per step and work the same whether `f` is a
//! hand-written function, an interpreted term or JIT-compiled code.
//!
//! Fixed-step methods return the state at every step. The adaptive
//! Dormand–Prince 5(4) integrator chooses its own steps to meet a tolerance
//! and returns the accepted points.

use crate::kernels::ode_adaptive::Event;
use crate::kernels::ode_adaptive::event_time_tolerance;

use serde::Deserialize;
use serde::Serialize;

/// Fixed-step methods.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
pub enum OdeSolverMethod {
    /// Forward Euler (order 1).
    Euler,
    /// Heun's method (order 2).
    Heun,
    /// Classical Runge–Kutta (order 4).
    RungeKutta4,
}

impl OdeSolverMethod {
    pub(crate) const fn stages(self) -> usize {
        match self {
            | Self::Euler => 1,
            | Self::Heun => 2,
            | Self::RungeKutta4 => 4,
        }
    }
}

/// A discrete trajectory: times and the state at each time.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Trajectory {
    /// Sample times in integration order.
    pub t: Vec<f64>,
    /// `y[k]` is the state at `t[k]`.
    pub y: Vec<Vec<f64>>,
}

impl Trajectory {
    /// The final state, or an empty slice for an empty trajectory.
    #[must_use]
    pub fn last(&self) -> &[f64] {
        self.y.last().map_or(&[], Vec::as_slice)
    }
}

/// Integrates `y' = f(t, y)` from `t_span.0` to `t_span.1` in `steps` equal
/// steps.
///
/// # Errors
/// Rejects zero steps, non-finite inputs, unrepresentable steps, or a
/// trajectory that cannot be reserved. Returns an error when the state stops
/// being finite.
pub fn solve_fixed(
    f: impl Fn(f64, &[f64], &mut [f64]),
    y0: &[f64],
    t_span: (f64, f64),
    steps: usize,
    method: OdeSolverMethod,
) -> Result<Trajectory, String> {
    crate::cancellation::checkpoint().map_err(|error| error.to_string())?;
    let (h, mut out) = fixed_grid(y0, t_span, steps)?;
    out.t.push(t_span.0);
    out.y.push(y0.to_vec());
    integrate_fixed(&f, y0, t_span, steps, method, h, |time, state| {
        out.t.push(time);
        out.y.push(state.to_vec());
    })?;
    Ok(out)
}

/// Integrates on the same fixed grid as [`solve_fixed`] and returns its final
/// state using scratch space proportional to the state dimension.
///
/// # Errors
/// Rejects zero steps, non-finite inputs, unrepresentable steps, non-finite
/// states, or cancellation.
pub fn solve_fixed_final(
    f: impl Fn(f64, &[f64], &mut [f64]),
    y0: &[f64],
    t_span: (f64, f64),
    steps: usize,
    method: OdeSolverMethod,
) -> Result<Vec<f64>, String> {
    crate::cancellation::checkpoint().map_err(|error| error.to_string())?;
    let h = fixed_step_size(y0, t_span, steps)?;
    integrate_fixed(&f, y0, t_span, steps, method, h, |_, _| {})
}

#[allow(clippy::float_cmp, reason = "detects a step rounded to the same time")]
fn integrate_fixed<F: Fn(f64, &[f64], &mut [f64])>(
    f: &F,
    y0: &[f64],
    t_span: (f64, f64),
    steps: usize,
    method: OdeSolverMethod,
    h: f64,
    mut record: impl FnMut(f64, &[f64]),
) -> Result<Vec<f64>, String> {
    let mut y = y0.to_vec();
    let (mut k, mut stage) = fixed_scratch(y.len(), method);
    for step in 0..steps {
        crate::cancellation::checkpoint().map_err(|error| error.to_string())?;
        let t = t_span.0 + h * step as f64;
        if t_span.0 != t_span.1 && (h == 0.0 || t + h == t) {
            return Err("The step is too small to advance time.".to_string());
        }
        advance(method, f, t, &mut y, h, &mut k, &mut stage);
        crate::cancellation::checkpoint().map_err(|error| error.to_string())?;
        if y.iter().any(|v| !v.is_finite()) {
            return Err("Overflow or invalid value encountered during ODE solving.".to_string());
        }
        record(t_span.0 + h * (step + 1) as f64, &y);
    }
    Ok(y)
}

#[allow(clippy::float_cmp, reason = "detects a step rounded to the same time")]
fn fixed_step_size(
    y0: &[f64],
    t_span: (f64, f64),
    steps: usize,
) -> Result<f64, String> {
    if steps == 0 {
        return Err("The number of steps must be positive.".to_owned());
    }
    let duration = t_span.1 - t_span.0;
    if !duration.is_finite() || y0.iter().any(|value| !value.is_finite()) {
        return Err("The time interval length and initial state must be finite.".to_owned());
    }
    let h = duration / steps as f64;
    if duration != 0.0 && t_span.0 + h == t_span.0 {
        return Err("The step is too small to advance time.".to_owned());
    }
    Ok(h)
}

fn fixed_grid(
    y0: &[f64],
    t_span: (f64, f64),
    steps: usize,
) -> Result<(f64, Trajectory), String> {
    let h = fixed_step_size(y0, t_span, steps)?;
    let capacity = steps
        .checked_add(1)
        .ok_or_else(|| "The number of trajectory samples is too large.".to_owned())?;
    let mut t = Vec::new();
    let mut y = Vec::new();
    t.try_reserve_exact(capacity)
        .and_then(|()| y.try_reserve_exact(capacity))
        .map_err(|error| format!("Cannot reserve the ODE trajectory: {error}"))?;
    Ok((h, Trajectory { t, y }))
}

fn fixed_scratch(
    n: usize,
    method: OdeSolverMethod,
) -> ([Vec<f64>; 4], Vec<f64>) {
    let stages = method.stages();
    (
        std::array::from_fn(|i| {
            if i < stages {
                vec![0.0; n]
            } else {
                Vec::new()
            }
        }),
        if stages > 1 {
            vec![0.0; n]
        } else {
            Vec::new()
        },
    )
}

/// One step of `method` from `(t, y)` with step `h`, in place.
///
/// The scratch buffers are passed in so that a caller taking many steps does
/// not allocate per step; the event locator below relies on being able to take
/// a *partial* step from the same state, which is what makes it exact for a
/// one-step method.
fn advance<F: Fn(f64, &[f64], &mut [f64])>(
    method: OdeSolverMethod,
    f: &F,
    t: f64,
    y: &mut [f64],
    h: f64,
    k: &mut [Vec<f64>; 4],
    stage: &mut [f64],
) {
    let [k0, k1, k2, k3] = k;
    match method {
        | OdeSolverMethod::Euler => {
            f(t, y, k0);
            for (value, &rate) in y.iter_mut().zip(k0.iter()) {
                *value += h * rate;
            }
        },
        | OdeSolverMethod::Heun => {
            f(t, y, k0);
            for ((stage, &value), &rate) in stage.iter_mut().zip(y.iter()).zip(k0.iter()) {
                *stage = value + h * rate;
            }
            f(t + h, stage, k1);
            for ((value, &first), &second) in y.iter_mut().zip(k0.iter()).zip(k1.iter()) {
                *value += 0.5 * h * (first + second);
            }
        },
        | OdeSolverMethod::RungeKutta4 => {
            f(t, y, k0);
            for ((stage, &value), &rate) in stage.iter_mut().zip(y.iter()).zip(k0.iter()) {
                *stage = value + 0.5 * h * rate;
            }
            f(t + 0.5 * h, stage, k1);
            for ((stage, &value), &rate) in stage.iter_mut().zip(y.iter()).zip(k1.iter()) {
                *stage = value + 0.5 * h * rate;
            }
            f(t + 0.5 * h, stage, k2);
            for ((stage, &value), &rate) in stage.iter_mut().zip(y.iter()).zip(k2.iter()) {
                *stage = value + h * rate;
            }
            f(t + h, stage, k3);
            for ((((value, &a), &b), &c), &d) in y
                .iter_mut()
                .zip(k0.iter())
                .zip(k1.iter())
                .zip(k2.iter())
                .zip(k3.iter())
            {
                *value += h / 6.0 * (a + 2.0 * b + 2.0 * c + d);
            }
        },
    }
}

/// One event located by [`solve_fixed_events`].
#[derive(Clone, Debug, PartialEq)]
pub struct FixedEvent {
    /// Index into the event list.
    pub index: usize,
    /// The time it happened.
    pub t: f64,
    /// The state at that time, as `method` computes it.
    pub y: Vec<f64>,
    /// Whether the integration stopped here.
    pub terminal: bool,
}

/// Integrates in `steps` equal steps and locates events.
///
/// # How an event is located without a dense output
///
/// A fixed-step method has no interpolant, so there is nothing to Brent-iterate
/// on. It is a *one-step* method, however, which means the state at any
/// intermediate time is just another application of the same step from the same
/// start: the locator bisects the step length, re-integrating the partial step
/// each time, until the bracket is below `1e-13`. That is exact for the method —
/// the located state is a point of the method's own trajectory, not of an
/// interpolant — and it costs about fifty evaluations per event, which is the
/// right trade for a method that is cheap by construction.
///
/// A `terminal` event stops the integration; the trajectory then ends at the
/// event point.
///
/// # Errors
/// As [`solve_fixed`].
// Equality detects an unrepresentable step/interval endpoint, not approximate convergence.
#[allow(clippy::float_cmp)]
pub fn solve_fixed_events<F: Fn(f64, &[f64], &mut [f64])>(
    f: F,
    y0: &[f64],
    t_span: (f64, f64),
    steps: usize,
    method: OdeSolverMethod,
    events: &[Event<'_>],
) -> Result<(Trajectory, Vec<FixedEvent>), String> {
    crate::cancellation::checkpoint().map_err(|error| error.to_string())?;
    let (h, mut out) = fixed_grid(y0, t_span, steps)?;
    let n = y0.len();
    let mut y = y0.to_vec();
    out.t.push(t_span.0);
    out.y.push(y.clone());
    let mut hits: Vec<FixedEvent> = Vec::new();
    let (mut k, mut stage) = fixed_scratch(n, method);
    let (mut probe_k, mut probe_stage) = fixed_scratch(n, method);
    let mut gprev: Vec<f64> = events.iter().map(|event| (event.g)(t_span.0, y0)).collect();
    let mut next = y.clone();
    let mut step_hits: Vec<FixedEvent> = Vec::new();
    for step in 0..steps {
        crate::cancellation::checkpoint().map_err(|error| error.to_string())?;
        let t = t_span.0 + h * step as f64;
        if t_span.0 != t_span.1 && (h == 0.0 || t + h == t) {
            return Err("The step is too small to advance time.".to_string());
        }
        next.copy_from_slice(&y);
        advance(method, &f, t, &mut next, h, &mut k, &mut stage);
        crate::cancellation::checkpoint().map_err(|error| error.to_string())?;
        if next.iter().any(|v| !v.is_finite()) {
            return Err("Overflow or invalid value encountered during ODE solving.".to_string());
        }
        // Locate every crossing inside this step, in time order.
        for (index, event) in events.iter().enumerate() {
            let g1 = (event.g)(t + h, &next);
            let g0 = gprev.get(index).copied().unwrap_or(g1);
            let crossing = (g0 < 0.0 && g1 >= 0.0 && event.direction >= 0)
                || (g0 > 0.0 && g1 <= 0.0 && event.direction <= 0);
            if crossing {
                let (fraction, state) = locate_in_step(
                    &f,
                    method,
                    t,
                    &y,
                    h,
                    g0,
                    event,
                    &mut probe_k,
                    &mut probe_stage,
                );
                step_hits.push(FixedEvent {
                    index,
                    t: t + fraction * h,
                    y: state,
                    terminal: event.terminal,
                });
            }
            if let Some(slot) = gprev.get_mut(index) {
                *slot = g1;
            }
        }
        step_hits.sort_by(|a, b| {
            let (ta, tb) = if h >= 0.0 {
                (a.t, b.t)
            } else {
                (-a.t, -b.t)
            };
            ta.total_cmp(&tb)
        });
        if let Some(position) = step_hits.iter().position(|hit| hit.terminal) {
            if let Some(last) = step_hits.get(position) {
                out.t.push(last.t);
                out.y.push(last.y.clone());
            }
            hits.extend(step_hits.into_iter().take(position + 1));
            return Ok((out, hits));
        }
        hits.append(&mut step_hits);
        std::mem::swap(&mut y, &mut next);
        out.t.push(t_span.0 + h * (step + 1) as f64);
        out.y.push(y.clone());
    }
    Ok((out, hits))
}

/// Bisects the step `[t, t + h]` for the event's root.
///
/// Returns the fraction of the step at which `g` changes sign, and the state
/// there as `method` computes it.
#[allow(clippy::too_many_arguments)]
// Equality detects an unrepresentable step/interval endpoint, not approximate convergence.
#[allow(clippy::float_cmp)]
fn locate_in_step<F: Fn(f64, &[f64], &mut [f64])>(
    f: &F,
    method: OdeSolverMethod,
    t: f64,
    y: &[f64],
    h: f64,
    g0: f64,
    event: &Event<'_>,
    k: &mut [Vec<f64>; 4],
    stage: &mut [f64],
) -> (f64, Vec<f64>) {
    let (mut lo, mut hi) = (0.0_f64, 1.0_f64);
    let tolerance = event_time_tolerance(t, t + h);
    let mut state = y.to_vec();
    for _ in 0..80 {
        let mid = f64::midpoint(lo, hi);
        if mid == lo || mid == hi {
            break;
        }
        state.copy_from_slice(y);
        advance(method, f, t, &mut state, mid * h, k, stage);
        let gm = (event.g)(t + mid * h, &state);
        if gm == 0.0 {
            return (mid, state);
        }
        if (gm < 0.0) == (g0 < 0.0) {
            lo = mid;
        } else {
            hi = mid;
        }
        if (hi - lo) * h.abs() <= tolerance {
            break;
        }
    }
    // Report the state at the bracket's far end, which is on the other side.
    state.copy_from_slice(y);
    advance(method, f, t, &mut state, hi * h, k, stage);
    (hi, state)
}

// Dormand–Prince 5(4) coefficients.
const C: [f64; 7] = [0.0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0];
const A: [[f64; 6]; 7] = [
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0],
    [44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0],
    [
        19_372.0 / 6_561.0,
        -25_360.0 / 2_187.0,
        64_448.0 / 6_561.0,
        -212.0 / 729.0,
        0.0,
        0.0,
    ],
    [
        9_017.0 / 3_168.0,
        -355.0 / 33.0,
        46_732.0 / 5_247.0,
        49.0 / 176.0,
        -5_103.0 / 18_656.0,
        0.0,
    ],
    [
        35.0 / 384.0,
        0.0,
        500.0 / 1_113.0,
        125.0 / 192.0,
        -2_187.0 / 6_784.0,
        11.0 / 84.0,
    ],
];
/// Difference between the 5th- and 4th-order weights.
const E: [f64; 7] = [
    71.0 / 57_600.0,
    0.0,
    -71.0 / 16_695.0,
    71.0 / 1_920.0,
    -17_253.0 / 339_200.0,
    22.0 / 525.0,
    -1.0 / 40.0,
];

/// Integrates `y' = f(t, y)` over `t_span` with the adaptive Dormand–Prince
/// 5(4) method, keeping the local error below `atol + rtol * |y|`.
///
/// # Errors
/// Returns an error when the step size underflows, `max_steps` accepted or
/// rejected steps are exhausted, or the state stops being finite.
#[allow(clippy::float_cmp)] // exact comparison against a sentinel / integer-valued input is intended
pub fn solve_adaptive(
    f: impl Fn(f64, &[f64], &mut [f64]),
    y0: &[f64],
    t_span: (f64, f64),
    rtol: f64,
    atol: f64,
    max_steps: usize,
) -> Result<Trajectory, String> {
    crate::cancellation::checkpoint().map_err(|error| error.to_string())?;
    let n = y0.len();
    let (t0, t1) = t_span;
    let direction = if t1 >= t0 { 1.0 } else { -1.0 };
    let mut t = t0;
    let mut y = y0.to_vec();
    let mut out = Trajectory {
        t: vec![t],
        y: vec![y.clone()],
    };
    if t0 == t1 {
        return Ok(out);
    }
    let mut k = vec![vec![0.0; n]; 7];
    let mut stage = vec![0.0; n];
    let mut next = vec![0.0; n];
    let mut h = direction * ((t1 - t0).abs() * 1e-2).max(1e-8);
    f(t, &y, &mut k[0]);
    for _ in 0..max_steps {
        crate::cancellation::checkpoint().map_err(|error| error.to_string())?;
        if (t1 - t) * direction <= 0.0 {
            return Ok(out);
        }
        if (t + h - t1) * direction > 0.0 {
            h = t1 - t;
        }
        for s in 1..7 {
            for i in 0..n {
                let mut acc = 0.0;
                for (j, a) in A[s].iter().enumerate().take(s) {
                    acc += a * k[j][i];
                }
                stage[i] = y[i] + h * acc;
            }
            if s == 6 {
                // The last stage is evaluated at the new point (FSAL).
                next.copy_from_slice(&stage);
            }
            let (done, rest) = k.split_at_mut(s);
            let _ = done;
            f(t + C[s] * h, &stage, &mut rest[0]);
        }
        let mut err = 0.0_f64;
        for i in 0..n {
            let mut e = 0.0;
            for (j, w) in E.iter().enumerate() {
                e += w * k[j][i];
            }
            let scale = atol + rtol * y[i].abs().max(next[i].abs());
            err = err.max((h * e).abs() / scale);
        }
        if !err.is_finite() || next.iter().any(|v| !v.is_finite()) {
            return Err("Overflow or invalid value encountered during ODE solving.".to_string());
        }
        if err <= 1.0 {
            t += h;
            y.copy_from_slice(&next);
            out.t.push(t);
            out.y.push(y.clone());
            k.swap(0, 6);
        }
        let factor = if err == 0.0 {
            5.0
        } else {
            (0.9 * err.powf(-0.2)).clamp(0.2, 5.0)
        };
        h *= factor;
        if h.abs() < 1e-14 * t.abs().max(1.0) {
            return Err("Step size underflow: the problem may be stiff or singular.".to_string());
        }
    }
    if (t1 - t) * direction <= 0.0 {
        Ok(out)
    } else {
        Err("Maximum number of steps exceeded.".to_string())
    }
}

/// Step-size control mode for unified ODE, geometric, symplectic, and PDE
/// solvers.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub enum OdeStepMode {
    /// Fixed step size `dt` (the final step may be adjusted to land on `t1`).
    Fixed {
        /// Positive step size magnitude.
        dt: f64,
    },
    /// Adaptive step-size control with local error tolerances.
    Adaptive {
        /// Relative error tolerance.
        rtol: f64,
        /// Absolute error tolerance.
        atol: f64,
        /// Optional initial step size (`None` selects automatically).
        dt_init: Option<f64>,
        /// Optional minimum step size (`None` uses representability / `0.0`).
        dt_min: Option<f64>,
        /// Optional maximum step size (`None` is unbounded).
        dt_max: Option<f64>,
        /// Maximum number of attempted steps (accepted plus rejected).
        max_steps: usize,
    },
}

impl Default for OdeStepMode {
    fn default() -> Self {
        Self::Adaptive {
            rtol: 1e-8,
            atol: 1e-10,
            dt_init: None,
            dt_min: None,
            dt_max: None,
            max_steps: 1_000_000,
        }
    }
}

/// The most steps a fixed-step integration may be asked for (the solution
/// stores every one of them).
pub(crate) const MAX_FIXED_STEPS: usize = 100_000_000;

impl OdeStepMode {
    /// Constructs a fixed-step mode with step size `dt`.
    #[must_use]
    pub const fn fixed(dt: f64) -> Self {
        Self::Fixed { dt }
    }

    /// Constructs an adaptive-step mode with the given tolerances and default
    /// step bounds.
    #[must_use]
    pub const fn adaptive(
        rtol: f64,
        atol: f64,
    ) -> Self {
        Self::Adaptive {
            rtol,
            atol,
            dt_init: None,
            dt_min: None,
            dt_max: None,
            max_steps: 1_000_000,
        }
    }

    /// Converts an [`OdeStepMode::Adaptive`] into
    /// [`OdeOptions`](crate::kernels::ode_adaptive::OdeOptions).
    ///
    /// # Errors
    /// Returns [`OdeError::Invalid`](crate::kernels::ode_adaptive::OdeError::Invalid)
    /// when called on [`OdeStepMode::Fixed`] or with invalid tolerances/bounds.
    pub fn to_ode_options(
        &self
    ) -> Result<crate::kernels::ode_adaptive::OdeOptions, crate::kernels::ode_adaptive::OdeError>
    {
        match *self {
            | Self::Fixed { .. } => Err(crate::kernels::ode_adaptive::OdeError::Invalid),
            | Self::Adaptive {
                rtol,
                atol,
                dt_init,
                dt_min,
                dt_max,
                max_steps,
            } => {
                let hmin = dt_min.unwrap_or(0.0);
                let hmax = dt_max.unwrap_or(f64::INFINITY);
                if !rtol.is_finite()
                    || rtol < 0.0
                    || !atol.is_finite()
                    || atol <= 0.0
                    || !hmin.is_finite()
                    || hmin < 0.0
                    || hmax.is_nan()
                    || hmax < hmin
                    || max_steps == 0
                    || dt_init.is_some_and(|v| !v.is_finite() || v <= 0.0)
                {
                    return Err(crate::kernels::ode_adaptive::OdeError::Invalid);
                }
                Ok(crate::kernels::ode_adaptive::OdeOptions {
                    rtol,
                    atol,
                    h0: dt_init,
                    hmin,
                    hmax,
                    max_steps,
                })
            },
        }
    }

    /// Computes the number of fixed steps for `[t0, t1]` when `self` is
    /// [`OdeStepMode::Fixed`].
    ///
    /// # Errors
    /// Returns [`OdeError::Invalid`](crate::kernels::ode_adaptive::OdeError::Invalid)
    /// for non-fixed modes, non-positive/non-finite `dt`, or a step count above
    /// [`MAX_FIXED_STEPS`].
    pub(crate) fn fixed_steps(
        &self,
        t0: f64,
        t1: f64,
    ) -> Result<usize, crate::kernels::ode_adaptive::OdeError> {
        match *self {
            | Self::Fixed { dt } => {
                if !dt.is_finite() || dt <= 0.0 || !t0.is_finite() || !t1.is_finite() {
                    return Err(crate::kernels::ode_adaptive::OdeError::Invalid);
                }
                let steps = ((t1 - t0).abs() / dt).ceil().max(1.0);
                // A tiny `dt` over a long span would otherwise allocate until
                // memory runs out (`inf` and huge counts saturate the cast).
                if steps > MAX_FIXED_STEPS as f64 {
                    return Err(crate::kernels::ode_adaptive::OdeError::Invalid);
                }
                Ok(steps as usize)
            },
            | Self::Adaptive { .. } => Err(crate::kernels::ode_adaptive::OdeError::Invalid),
        }
    }
}

/// Unified selection of initial-value ODE integrators.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
pub enum OdeMethod {
    /// Forward Euler (order 1; adaptive via Richardson step-doubling).
    Euler,
    /// Heun's method (order 2; adaptive via Richardson step-doubling).
    Heun,
    /// Classical Runge–Kutta (order 4; adaptive via Richardson step-doubling).
    RungeKutta4,
    /// Dormand–Prince 5(4) explicit Runge–Kutta pair.
    Dopri5,
    /// Shampine–Reichelt L-stable Rosenbrock 2(3) pair.
    Rosenbrock23,
    /// Alexander L-stable SDIRK3(2) with finite-difference Jacobian.
    Sdirk3,
    /// 3-stage, order-5 Radau IIA implicit Runge–Kutta.
    Radau5,
    /// Variable-order (1..5) Backward Differentiation Formulas.
    Bdf,
    /// Exponential Rosenbrock 3(2) method.
    Exprb32,
}

fn explicit_fixed_with_events<F: Fn(f64, &[f64], &mut [f64])>(
    method: OdeSolverMethod,
    f: F,
    t0: f64,
    t1: f64,
    y0: &[f64],
    steps: usize,
    events: &[Event<'_>],
) -> Result<crate::kernels::ode_adaptive::OdeSolution, crate::kernels::ode_adaptive::OdeError> {
    use crate::kernels::ode_adaptive::DenseKind;
    use crate::kernels::ode_adaptive::DenseStep;
    use crate::kernels::ode_adaptive::OdeError;
    use crate::kernels::ode_adaptive::OdeOptions;
    use crate::kernels::ode_adaptive::OdeSolution;
    use crate::kernels::ode_adaptive::StepWindow;
    use crate::kernels::ode_adaptive::process_events;
    use crate::kernels::ode_contract::fixed_step;
    use crate::kernels::ode_contract::initial_events;
    use crate::kernels::ode_contract::validate_problem;

    if steps == 0 {
        return Err(OdeError::Invalid);
    }
    validate_problem(t0, t1, y0, &OdeOptions::default(), events)?;
    if crate::kernels::exact_eq(t0, t1) {
        return Ok(OdeSolution::start(t0, y0));
    }
    let f = crate::kernels::ode_contract::checked_rhs(f);
    let n = y0.len();
    let dir = if t1 > t0 { 1.0 } else { -1.0 };
    let mut t = t0;
    let mut y = y0.to_vec();
    let (mut k, mut stage) = fixed_scratch(n, method);
    let mut f0 = vec![0.0; n];
    let mut f1 = vec![0.0; n];
    f(t, &y, &mut f0);
    if f0.iter().any(|v| !v.is_finite()) {
        return Err(OdeError::NonFinite);
    }
    let evals_per_step = match method {
        | OdeSolverMethod::Euler => 1usize,
        | OdeSolverMethod::Heun => 2usize,
        | OdeSolverMethod::RungeKutta4 => 4usize,
    };
    let mut nfev = 1usize;
    let mut sol = OdeSolution::start(t0, y0);
    let mut gprev = initial_events(events, t0, y0)?;
    for idx in 1..=steps {
        crate::cancellation::checkpoint()?;
        let (tn, h) = fixed_step(t0, t1, idx, steps, t)?;
        let mut ynew = y.clone();
        advance(method, &f, t, &mut ynew, h, &mut k, &mut stage);
        nfev = nfev.saturating_add(evals_per_step);
        if ynew.iter().any(|v| !v.is_finite()) {
            return Err(OdeError::NonFinite);
        }
        f(tn, &ynew, &mut f1);
        nfev = nfev.saturating_add(1);
        if f1.iter().any(|v| !v.is_finite()) {
            return Err(OdeError::NonFinite);
        }
        let step_dense = DenseStep {
            t0: t,
            h,
            y0: y.clone(),
            kind: DenseKind::Hermite {
                y1: ynew.clone(),
                f0: f0.clone(),
                f1: f1.clone(),
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
        f0.copy_from_slice(&f1);
    }
    sol.nfev = nfev;
    Ok(sol)
}

#[allow(clippy::float_cmp)]
fn explicit_doubling_adaptive<F: Fn(f64, &[f64], &mut [f64])>(
    method: OdeSolverMethod,
    f: F,
    t0: f64,
    t1: f64,
    y0: &[f64],
    opts: &crate::kernels::ode_adaptive::OdeOptions,
    events: &[Event<'_>],
) -> Result<crate::kernels::ode_adaptive::OdeSolution, crate::kernels::ode_adaptive::OdeError> {
    use crate::kernels::ode_adaptive::DenseKind;
    use crate::kernels::ode_adaptive::DenseStep;
    use crate::kernels::ode_adaptive::OdeError;
    use crate::kernels::ode_adaptive::OdeSolution;
    use crate::kernels::ode_adaptive::StepWindow;
    use crate::kernels::ode_adaptive::err_norm;
    use crate::kernels::ode_adaptive::process_events;
    use crate::kernels::ode_contract::InitialStepOptions;
    use crate::kernels::ode_contract::initial_events;
    use crate::kernels::ode_contract::initial_step;
    use crate::kernels::ode_contract::validate_problem;

    validate_problem(t0, t1, y0, opts, events)?;
    if crate::kernels::exact_eq(t0, t1) {
        return Ok(OdeSolution::start(t0, y0));
    }
    let f = crate::kernels::ode_contract::checked_rhs(f);
    let n = y0.len();
    let dir = if t1 > t0 { 1.0 } else { -1.0 };
    let (order, evals_per_substep) = match method {
        | OdeSolverMethod::Euler => (1.0_f64, 1usize),
        | OdeSolverMethod::Heun => (2.0_f64, 2usize),
        | OdeSolverMethod::RungeKutta4 => (4.0_f64, 4usize),
    };
    let denom = 2.0_f64.powf(order) - 1.0;
    let mut t = t0;
    let mut y = y0.to_vec();
    let mut f0 = vec![0.0; n];
    f(t, &y, &mut f0);
    if f0.iter().any(|v| !v.is_finite()) {
        return Err(OdeError::NonFinite);
    }
    let mut nfev = 1usize;
    let mut h = opts
        .h0
        .unwrap_or_else(|| {
            let est = initial_step(
                &f,
                t0,
                y0,
                &f0,
                InitialStepOptions {
                    rtol: opts.rtol,
                    atol: opts.atol,
                    bound: t1,
                    max_step: opts.hmax,
                    order,
                },
            );
            nfev = nfev.saturating_add(est.evaluations);
            est.step
        })
        .abs()
        .min(opts.hmax)
        * dir;
    let mut sol = OdeSolution::start(t0, y0);
    let mut gprev = initial_events(events, t0, y0)?;
    let (mut k, mut stage) = fixed_scratch(n, method);
    let mut f1 = vec![0.0; n];
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
        let mut y_full = y.clone();
        advance(method, &f, t, &mut y_full, h, &mut k, &mut stage);
        let mut y_half = y.clone();
        advance(method, &f, t, &mut y_half, 0.5 * h, &mut k, &mut stage);
        advance(
            method,
            &f,
            t + 0.5 * h,
            &mut y_half,
            0.5 * h,
            &mut k,
            &mut stage,
        );
        nfev = nfev.saturating_add(evals_per_substep.saturating_mul(3));
        if y_half.iter().chain(&y_full).any(|v| !v.is_finite()) {
            h *= 0.25;
            sol.rejected = sol.rejected.saturating_add(1);
            continue;
        }
        let err_vec: Vec<f64> = y_half
            .iter()
            .zip(&y_full)
            .map(|(yh, yf)| (*yh - *yf).abs() / denom)
            .collect();
        let en = err_norm(&err_vec, &y, &y_half, opts.rtol, opts.atol);
        if en.is_nan() {
            return Err(OdeError::NonFinite);
        }
        if en <= 1.0 {
            let tn = t + h;
            f(tn, &y_half, &mut f1);
            nfev = nfev.saturating_add(1);
            if f1.iter().any(|v| !v.is_finite()) {
                return Err(OdeError::NonFinite);
            }
            let step_dense = DenseStep {
                t0: t,
                h,
                y0: y.clone(),
                kind: DenseKind::Hermite {
                    y1: y_half.clone(),
                    f0: f0.clone(),
                    f1: f1.clone(),
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
                    ynew: &y_half,
                    dir,
                },
            )? {
                sol.nfev = nfev;
                return Ok(sol);
            }
            sol.t.push(tn);
            sol.y.push(y_half.clone());
            t = tn;
            y = y_half;
            f0.copy_from_slice(&f1);
        } else {
            sol.rejected = sol.rejected.saturating_add(1);
        }
        let factor = if en == 0.0 {
            5.0
        } else {
            (0.9 * en.powf(-1.0 / (order + 1.0))).clamp(0.2, 5.0)
        };
        h = (h.abs() * factor).min(opts.hmax) * dir;
    }
    sol.nfev = nfev;
    Ok(sol)
}

/// Unified initial-value ODE solver supporting every [`OdeMethod`] in both
/// fixed and adaptive [`OdeStepMode`]s without events.
///
/// # Errors
/// Returns [`OdeError`](crate::kernels::ode_adaptive::OdeError) when inputs are
/// invalid, the step budget or step floor is reached, the state becomes
/// non-finite, or an implicit solve fails.
pub fn solve_ode<F: Fn(f64, &[f64], &mut [f64])>(
    f: F,
    y0: &[f64],
    t_span: (f64, f64),
    method: OdeMethod,
    step_mode: OdeStepMode,
) -> Result<crate::kernels::ode_adaptive::OdeSolution, crate::kernels::ode_adaptive::OdeError> {
    solve_ode_events(f, y0, t_span, method, step_mode, &[])
}

/// Unified initial-value ODE solver supporting every [`OdeMethod`] in both
/// fixed and adaptive [`OdeStepMode`]s with direction-aware / terminal [`Event`]
/// detection.
///
/// # Errors
/// As [`solve_ode`], or [`OdeError::Invalid`](crate::kernels::ode_adaptive::OdeError::Invalid)
/// if an event callback returns a non-finite value.
pub fn solve_ode_events<F: Fn(f64, &[f64], &mut [f64])>(
    f: F,
    y0: &[f64],
    t_span: (f64, f64),
    method: OdeMethod,
    step_mode: OdeStepMode,
    events: &[Event<'_>],
) -> Result<crate::kernels::ode_adaptive::OdeSolution, crate::kernels::ode_adaptive::OdeError> {
    let (t0, t1) = t_span;
    match (method, step_mode) {
        | (OdeMethod::Euler, OdeStepMode::Fixed { .. }) => {
            let steps = step_mode.fixed_steps(t0, t1)?;
            explicit_fixed_with_events(OdeSolverMethod::Euler, f, t0, t1, y0, steps, events)
        },
        | (OdeMethod::Euler, OdeStepMode::Adaptive { .. }) => {
            let opts = step_mode.to_ode_options()?;
            explicit_doubling_adaptive(OdeSolverMethod::Euler, f, t0, t1, y0, &opts, events)
        },
        | (OdeMethod::Heun, OdeStepMode::Fixed { .. }) => {
            let steps = step_mode.fixed_steps(t0, t1)?;
            explicit_fixed_with_events(OdeSolverMethod::Heun, f, t0, t1, y0, steps, events)
        },
        | (OdeMethod::Heun, OdeStepMode::Adaptive { .. }) => {
            let opts = step_mode.to_ode_options()?;
            explicit_doubling_adaptive(OdeSolverMethod::Heun, f, t0, t1, y0, &opts, events)
        },
        | (OdeMethod::RungeKutta4, OdeStepMode::Fixed { .. }) => {
            let steps = step_mode.fixed_steps(t0, t1)?;
            explicit_fixed_with_events(OdeSolverMethod::RungeKutta4, f, t0, t1, y0, steps, events)
        },
        | (OdeMethod::RungeKutta4, OdeStepMode::Adaptive { .. }) => {
            let opts = step_mode.to_ode_options()?;
            explicit_doubling_adaptive(OdeSolverMethod::RungeKutta4, f, t0, t1, y0, &opts, events)
        },
        | (OdeMethod::Dopri5, OdeStepMode::Fixed { .. }) => {
            let steps = step_mode.fixed_steps(t0, t1)?;
            crate::kernels::ode_adaptive::dopri5_fixed(f, t0, t1, y0, steps, events)
        },
        | (OdeMethod::Dopri5, OdeStepMode::Adaptive { .. }) => {
            let opts = step_mode.to_ode_options()?;
            crate::kernels::ode_adaptive::dopri5(f, t0, t1, y0, &opts, events)
        },
        | (OdeMethod::Rosenbrock23, OdeStepMode::Fixed { .. }) => {
            let steps = step_mode.fixed_steps(t0, t1)?;
            crate::kernels::ode_adaptive::rosenbrock23_fixed_events(f, t0, t1, y0, steps, events)
        },
        | (OdeMethod::Rosenbrock23, OdeStepMode::Adaptive { .. }) => {
            let opts = step_mode.to_ode_options()?;
            crate::kernels::ode_adaptive::rosenbrock23_events(f, t0, t1, y0, &opts, events)
        },
        | (OdeMethod::Sdirk3, OdeStepMode::Fixed { .. }) => {
            let steps = step_mode.fixed_steps(t0, t1)?;
            crate::kernels::ode_sdirk::sdirk3_fd_fixed_events(f, t0, t1, y0, steps, events)
        },
        | (OdeMethod::Sdirk3, OdeStepMode::Adaptive { .. }) => {
            let opts = step_mode.to_ode_options()?;
            let sdirk_opts = crate::kernels::ode_sdirk::SdirkOptions {
                rtol: opts.rtol,
                atol: opts.atol,
                h0: opts.h0,
                hmax: opts.hmax,
                max_steps: opts.max_steps,
                ..crate::kernels::ode_sdirk::SdirkOptions::default()
            };
            crate::kernels::ode_sdirk::sdirk3_fd_events(f, t0, t1, y0, &sdirk_opts, events)
        },
        | (OdeMethod::Radau5, OdeStepMode::Fixed { .. }) => {
            let steps = step_mode.fixed_steps(t0, t1)?;
            crate::kernels::ode_stiff::radau5_fixed_events(f, t0, t1, y0, steps, events)
        },
        | (OdeMethod::Radau5, OdeStepMode::Adaptive { .. }) => {
            let opts = step_mode.to_ode_options()?;
            crate::kernels::ode_stiff::radau5_events(f, t0, t1, y0, &opts, events)
        },
        | (OdeMethod::Bdf, OdeStepMode::Fixed { .. }) => {
            let steps = step_mode.fixed_steps(t0, t1)?;
            crate::kernels::ode_stiff::bdf_fixed_events(f, t0, t1, y0, steps, events)
        },
        | (OdeMethod::Bdf, OdeStepMode::Adaptive { .. }) => {
            let opts = step_mode.to_ode_options()?;
            crate::kernels::ode_stiff::bdf_events(f, t0, t1, y0, &opts, 5, events)
        },
        | (OdeMethod::Exprb32, _) => {
            let _ = events;
            crate::kernels::ode_ext::exprb32_with(f, None, t0, t1, y0, step_mode)
        },
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn fixed_steps_refuses_an_unbounded_step_count() {
        use crate::kernels::ode_adaptive::OdeError;
        let tiny = OdeStepMode::Fixed { dt: 1e-300 };
        assert_eq!(tiny.fixed_steps(0.0, 1.0), Err(OdeError::Invalid));
        let ok = OdeStepMode::Fixed { dt: 0.25 };
        assert_eq!(ok.fixed_steps(0.0, 1.0), Ok(4));
    }

    use super::*;

    fn decay(
        _t: f64,
        y: &[f64],
        dy: &mut [f64],
    ) {
        dy[0] = -2.0 * y[0];
    }

    #[test]
    fn fixed_step_orders_of_convergence() {
        let exact = (-2.0_f64).exp();
        let error = |method, steps| {
            let tr = solve_fixed(decay, &[1.0], (0.0, 1.0), steps, method)
                .unwrap_or_else(|e| panic!("{e}"));
            (tr.last()[0] - exact).abs()
        };
        for (method, order) in [
            (OdeSolverMethod::Euler, 1.0),
            (OdeSolverMethod::Heun, 2.0),
            (OdeSolverMethod::RungeKutta4, 4.0),
        ] {
            let ratio = error(method, 50) / error(method, 100);
            let observed = ratio.log2();
            assert!(
                (observed - order).abs() < 0.2,
                "{method:?}: observed order {observed}"
            );
        }
        assert!(solve_fixed(decay, &[1.0], (0.0, 1.0), 0, OdeSolverMethod::Euler).is_err());
    }

    #[test]
    fn adaptive_meets_tolerance() {
        // Harmonic oscillator over ten periods.
        let f = |_t: f64, y: &[f64], dy: &mut [f64]| {
            dy[0] = y[1];
            dy[1] = -y[0];
        };
        let t_end = 20.0 * std::f64::consts::PI;
        let tr = solve_adaptive(f, &[1.0, 0.0], (0.0, t_end), 1e-10, 1e-12, 100_000)
            .unwrap_or_else(|e| panic!("{e}"));
        assert!((tr.last()[0] - 1.0).abs() < 1e-7, "x = {}", tr.last()[0]);
        assert!(tr.last()[1].abs() < 1e-7);
        assert!((tr.t.last().copied().unwrap_or(0.0) - t_end).abs() < 1e-12);
        assert!(tr.t.windows(2).all(|w| w[1] > w[0]));
        // Far fewer steps than a fixed-step method would need.
        assert!(tr.t.len() < 5_000, "{} steps", tr.t.len());
    }

    #[test]
    fn adaptive_backwards_and_degenerate_spans() {
        let tr = solve_adaptive(decay, &[1.0], (1.0, 0.0), 1e-9, 1e-12, 10_000)
            .unwrap_or_else(|e| panic!("{e}"));
        assert!((tr.last()[0] - 2.0_f64.exp()).abs() < 1e-7);
        let tr = solve_adaptive(decay, &[1.0], (0.5, 0.5), 1e-9, 1e-12, 10)
            .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(tr.t, vec![0.5]);
    }

    #[test]
    fn blow_up_is_reported() {
        // y' = y^2, y(0) = 1 blows up at t = 1.
        let f = |_t: f64, y: &[f64], dy: &mut [f64]| dy[0] = y[0] * y[0];
        assert!(solve_adaptive(f, &[1.0], (0.0, 2.0), 1e-8, 1e-10, 100_000).is_err());
    }
}
