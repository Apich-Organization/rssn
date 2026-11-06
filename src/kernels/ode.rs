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

/// A discrete trajectory: times and the state at each time.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Trajectory {
    /// Sample times, strictly monotone.
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
/// Returns an error when `steps` is zero or the state stops being finite.
// Equality detects an unrepresentable step/interval endpoint, not approximate convergence.
#[allow(clippy::float_cmp)]
pub fn solve_fixed(
    f: impl Fn(f64, &[f64], &mut [f64]),
    y0: &[f64],
    t_span: (f64, f64),
    steps: usize,
    method: OdeSolverMethod,
) -> Result<Trajectory, String> {
    if steps == 0 {
        return Err("The number of steps must be positive.".to_string());
    }
    let n = y0.len();
    let h = (t_span.1 - t_span.0) / steps as f64;
    let mut y = y0.to_vec();
    let mut out = Trajectory {
        t: Vec::with_capacity(steps + 1),
        y: Vec::with_capacity(steps + 1),
    };
    out.t.push(t_span.0);
    out.y.push(y.clone());
    let mut k = vec![vec![0.0; n]; 4];
    let mut stage = vec![0.0; n];
    for step in 0..steps {
        let t = t_span.0 + h * step as f64;
        if t_span.0 != t_span.1 && (h == 0.0 || t + h == t) {
            return Err("The step is too small to advance time.".to_string());
        }
        advance(method, &f, t, &mut y, h, &mut k, &mut stage);
        if y.iter().any(|v| !v.is_finite()) {
            return Err("Overflow or invalid value encountered during ODE solving.".to_string());
        }
        out.t.push(t_span.0 + h * (step + 1) as f64);
        out.y.push(y.clone());
    }
    Ok(out)
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
    k: &mut [Vec<f64>],
    stage: &mut [f64],
) {
    let n = y.len();
    match method {
        | OdeSolverMethod::Euler => {
            f(t, y, &mut k[0]);
            for i in 0..n {
                y[i] += h * k[0][i];
            }
        },
        | OdeSolverMethod::Heun => {
            f(t, y, &mut k[0]);
            for i in 0..n {
                stage[i] = y[i] + h * k[0][i];
            }
            f(t + h, stage, &mut k[1]);
            for i in 0..n {
                y[i] += 0.5 * h * (k[0][i] + k[1][i]);
            }
        },
        | OdeSolverMethod::RungeKutta4 => {
            f(t, y, &mut k[0]);
            for i in 0..n {
                stage[i] = y[i] + 0.5 * h * k[0][i];
            }
            f(t + 0.5 * h, stage, &mut k[1]);
            for i in 0..n {
                stage[i] = y[i] + 0.5 * h * k[1][i];
            }
            f(t + 0.5 * h, stage, &mut k[2]);
            for i in 0..n {
                stage[i] = y[i] + h * k[2][i];
            }
            f(t + h, stage, &mut k[3]);
            for i in 0..n {
                y[i] += h / 6.0 * (k[0][i] + 2.0 * k[1][i] + 2.0 * k[2][i] + k[3][i]);
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
    if steps == 0 {
        return Err("The number of steps must be positive.".to_string());
    }
    let n = y0.len();
    let h = (t_span.1 - t_span.0) / steps as f64;
    let mut y = y0.to_vec();
    let mut out = Trajectory {
        t: Vec::with_capacity(steps + 1),
        y: Vec::with_capacity(steps + 1),
    };
    out.t.push(t_span.0);
    out.y.push(y.clone());
    let mut hits: Vec<FixedEvent> = Vec::new();
    let mut k = vec![vec![0.0; n]; 4];
    let mut stage = vec![0.0; n];
    let mut probe_k = vec![vec![0.0; n]; 4];
    let mut probe_stage = vec![0.0; n];
    let mut gprev: Vec<f64> = events.iter().map(|event| (event.g)(t_span.0, y0)).collect();
    for step in 0..steps {
        let t = t_span.0 + h * step as f64;
        if t_span.0 != t_span.1 && (h == 0.0 || t + h == t) {
            return Err("The step is too small to advance time.".to_string());
        }
        let mut next = y.clone();
        advance(method, &f, t, &mut next, h, &mut k, &mut stage);
        if next.iter().any(|v| !v.is_finite()) {
            return Err("Overflow or invalid value encountered during ODE solving.".to_string());
        }
        // Locate every crossing inside this step, in time order.
        let mut step_hits: Vec<FixedEvent> = Vec::new();
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
            let kept: Vec<FixedEvent> = step_hits.into_iter().take(position + 1).collect();
            if let Some(last) = kept.last() {
                out.t.push(last.t);
                out.y.push(last.y.clone());
            }
            hits.extend(kept);
            return Ok((out, hits));
        }
        hits.extend(step_hits);
        y = next;
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
    k: &mut [Vec<f64>],
    stage: &mut [f64],
) -> (f64, Vec<f64>) {
    let (mut lo, mut hi) = (0.0_f64, 1.0_f64);
    let tolerance = event_time_tolerance(t, t + h);
    for _ in 0..80 {
        let mid = f64::midpoint(lo, hi);
        if mid == lo || mid == hi { break; }
        let mut state = y.to_vec();
        advance(method, f, t, &mut state, mid * h, k, stage);
        let gm = (event.g)(t + mid * h, &state);
        if gm == 0.0 { return (mid, state); }
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
    let mut final_state = y.to_vec();
    advance(method, f, t, &mut final_state, hi * h, k, stage);
    (hi, final_state)
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

#[cfg(test)]
mod tests {
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
