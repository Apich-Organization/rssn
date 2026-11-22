//! PID control: first-order-plus-dead-time identification, tuning rules and
//! closed-loop simulation.
//!
//! Tuning starts from the plant step response. The time to reach 28.3 % and
//! 63.2 % of the final value fixes a first-order-plus-dead-time model
//! `K exp(-theta s) / (tau s + 1)` by the two-point method; a classical rule
//! then maps `(K, tau, theta)` to `(kp, ki, kd)`.

use super::linalg;
use super::transfer::TransferFunction;
use super::transfer::step_response;

/// PID parameterisation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PidForm {
    /// `kp e + ki integral(e) + kd de/dt`.
    Parallel,
    /// `kp (e + integral(e)/ti + td de/dt)`.
    Standard,
    /// `kp (1 + 1/(ti s))(1 + td s)`, the interacting form.
    Series,
}

/// Tuning rule applied to the identified first-order-plus-dead-time model.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TuningRule {
    /// Ziegler-Nichols reaction-curve rule.
    ZieglerNichols,
    /// Cohen-Coon rule.
    CohenCoon,
    /// Internal model control with a closed-loop time constant.
    Imc,
    /// AMIGO rule.
    Amigo,
    /// Chien-Hrones-Reswick rule (setpoint, 0 % overshoot).
    ChienHronesReswick,
}

/// A discrete PID controller in parallel form.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Pid {
    /// Proportional gain.
    pub kp: f64,
    /// Integral gain.
    pub ki: f64,
    /// Derivative gain.
    pub kd: f64,
    /// Derivative filter time constant; `None` disables filtering.
    pub tf: Option<f64>,
}

/// First-order-plus-dead-time model identified from a step response.
struct Fopdt {
    gain: f64,
    tau: f64,
    theta: f64,
}

/// Interpolated time at which `values` first crosses `level`.
fn crossing_time(
    t: &[f64],
    values: &[f64],
    level: f64,
) -> Option<f64> {
    for i in 1..values.len() {
        if values[i - 1] < level && values[i] >= level {
            let span = values[i] - values[i - 1];
            if span.abs() <= f64::MIN_POSITIVE {
                return Some(t[i]);
            }
            let frac = (level - values[i - 1]) / span;
            return Some(t[i - 1] + frac * (t[i] - t[i - 1]));
        }
    }
    None
}

/// Identifies a first-order-plus-dead-time model from the step response.
fn identify_fopdt(plant: &TransferFunction) -> Result<Fopdt, String> {
    let poles = plant.poles();
    let mut fastest_decay = f64::INFINITY;
    for pole in &poles {
        if pole.re < 0.0 {
            fastest_decay = fastest_decay.min(-pole.re);
        }
    }
    let horizon = if fastest_decay.is_finite() && fastest_decay > 1e-6 {
        (20.0 / fastest_decay).min(500.0)
    } else {
        100.0
    };
    let samples = 2001;
    let t: Vec<f64> = (0..samples)
        .map(|i| horizon * i as f64 / (samples - 1) as f64)
        .collect();
    let y = step_response(plant, &t)?;
    let tail: f64 = y.iter().rev().take(5).sum::<f64>() / 5.0;
    if !tail.is_finite() || tail.abs() < 1e-12 {
        return Err("the plant has no finite steady-state gain".to_string());
    }
    let normalised: Vec<f64> = y.iter().map(|v| v / tail).collect();
    let t28 = crossing_time(&t, &normalised, 0.283)
        .ok_or_else(|| "the plant does not reach 28.3% of its final value".to_string())?;
    let t63 = crossing_time(&t, &normalised, 0.632)
        .ok_or_else(|| "the plant does not reach 63.2% of its final value".to_string())?;
    let tau = 1.5 * (t63 - t28);
    if !(tau > 0.0) {
        return Err("the identified time constant is not positive".to_string());
    }
    let theta = (t63 - tau).max(1e-3 * tau);
    Ok(Fopdt {
        gain: tail,
        tau,
        theta,
    })
}

/// Tunes a PID controller for `plant` with the selected rule.
///
/// `lambda` is the closed-loop time constant of the IMC rule; it is ignored by
/// the other rules and defaults to `max(0.1 tau, theta)`.
///
/// # Errors
/// Returns an error when the plant step response does not identify a finite
/// first-order-plus-dead-time model or when `lambda` is not positive.
pub fn tune(
    plant: &TransferFunction,
    rule: TuningRule,
    lambda: Option<f64>,
) -> Result<Pid, String> {
    let model = identify_fopdt(plant)?;
    let k = model.gain;
    let tau = model.tau;
    let theta = model.theta;
    if lambda.is_some_and(|l| !(l.is_finite() && l > 0.0)) {
        return Err("lambda must be positive".to_string());
    }
    let (kp, ti, td) = match rule {
        | TuningRule::ZieglerNichols => (1.2 * tau / (k * theta), 2.0 * theta, 0.5 * theta),
        | TuningRule::CohenCoon => {
            let ratio = theta / tau;
            (
                (tau / (k * theta)) * (4.0 / 3.0 + ratio / 4.0),
                theta * (32.0 + 6.0 * ratio) / (13.0 + 8.0 * ratio),
                theta * 4.0 / (11.0 + 2.0 * ratio),
            )
        },
        | TuningRule::Imc => {
            let lam = lambda.unwrap_or_else(|| theta.max(0.1 * tau));
            (
                (tau + 0.5 * theta) / (k * (lam + 0.5 * theta)),
                tau + 0.5 * theta,
                tau * theta / (2.0 * tau + theta),
            )
        },
        | TuningRule::Amigo => {
            let ratio = tau / theta;
            (
                (0.2 + 0.45 * ratio) / k,
                theta * (0.4 * theta + 0.8 * tau) / (theta + 0.1 * tau),
                0.5 * theta * tau / (0.3 * theta + tau),
            )
        },
        | TuningRule::ChienHronesReswick => (0.6 * tau / (k * theta), tau, 0.5 * theta),
    };
    let ki = if ti.abs() > f64::MIN_POSITIVE {
        kp / ti
    } else {
        0.0
    };
    let kd = kp * td;
    Ok(Pid {
        kp,
        ki,
        kd,
        tf: Some(0.1 * td),
    })
}

/// Closed-loop step response of `plant` under `pid` on the time grid `t`.
///
/// The derivative is filtered with a first-order filter of time constant
/// [`Pid::tf`] when it is set. The plant is integrated with fourth-order
/// Runge-Kutta on each interval.
///
/// # Errors
/// Returns an error when the plant cannot be converted to state space or the
/// time grid is not strictly increasing.
pub fn pid_response(
    pid: &Pid,
    plant: &TransferFunction,
    t: &[f64],
    setpoint: f64,
) -> Result<Vec<f64>, String> {
    let ss = plant.to_state_space()?;
    let n = ss.n_states();
    let mut x = vec![0.0; n];
    let mut out = vec![0.0; t.len()];
    if t.is_empty() {
        return Ok(out);
    }
    let mut y_prev = 0.0_f64;
    let mut integral = 0.0_f64;
    let mut filtered = 0.0_f64;
    let mut error_prev = setpoint;
    for k in 1..t.len() {
        let h = t[k] - t[k - 1];
        if !(h > 0.0) {
            return Err("time vector must be strictly increasing".to_string());
        }
        let error = setpoint - y_prev;
        integral += error * h;
        let raw = (error - error_prev) / h;
        filtered = match pid.tf {
            | Some(tf) if tf > 0.0 => {
                let alpha = tf / (tf + h);
                alpha * filtered + (1.0 - alpha) * raw
            },
            | _ => raw,
        };
        let u = pid.kp * error + pid.ki * integral + pid.kd * filtered;
        x = linalg::rk4_linear(&ss.a, &ss.b, &x, &[u], h);
        let cx = linalg::matvec(&ss.c, &x);
        let du = ss.d.first().and_then(|r| r.first()).copied().unwrap_or(0.0);
        y_prev = cx.first().copied().unwrap_or(0.0) + du * u;
        out[k] = y_prev;
        error_prev = error;
    }
    Ok(out)
}

/// Simulates a positional PID over an error sequence with optional output
/// limits and conditional-integration anti-windup.
///
/// Returns the control signal at each sample.
pub fn simulate_pid(
    pid: &Pid,
    error: &[f64],
    dt: f64,
    u_limits: Option<(f64, f64)>,
) -> Vec<f64> {
    let mut out = Vec::with_capacity(error.len());
    if error.is_empty() || !(dt > 0.0) {
        return out;
    }
    let mut integral = 0.0_f64;
    let mut previous = error[0];
    for (k, &e) in error.iter().enumerate() {
        let derivative = if k == 0 {
            0.0
        } else {
            (e - previous) / dt
        };
        let mut u = pid.kp * e + pid.ki * integral + pid.kd * derivative;
        let mut saturated = false;
        if let Some((lo, hi)) = u_limits {
            if u > hi {
                u = hi;
                saturated = e > 0.0;
            } else if u < lo {
                u = lo;
                saturated = e < 0.0;
            }
        }
        out.push(u);
        if !saturated {
            integral += e * dt;
        }
        previous = e;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernels::control::transfer::TransferFunction;
    use crate::kernels::control::transfer::feedback;
    use crate::kernels::control::transfer::series;

    fn plant_third_order() -> TransferFunction {
        TransferFunction::new(vec![1.0], vec![1.0, 3.0, 3.0, 1.0]).unwrap()
    }

    fn pid_transfer(pid: &Pid) -> TransferFunction {
        let tf = pid.tf.unwrap_or(0.0);
        // C(s) = kp + ki/s + kd s/(1 + tf s).
        let num = vec![pid.kp * tf + pid.kd, pid.kp + pid.ki * tf, pid.ki];
        let den = vec![tf, 1.0, 0.0];
        TransferFunction::new(num, den).unwrap()
    }

    #[test]
    fn ziegler_nichols_returns_positive_gains() {
        let plant = plant_third_order();
        let pid = tune(&plant, TuningRule::ZieglerNichols, None).unwrap();
        assert!(pid.kp > 0.0);
        assert!(pid.ki > 0.0);
        assert!(pid.kd > 0.0);
        assert!(pid.tf.is_some_and(|v| v > 0.0));
    }

    #[test]
    fn ziegler_nichols_closed_loop_is_stable() {
        let plant = plant_third_order();
        let pid = tune(&plant, TuningRule::ZieglerNichols, None).unwrap();
        let controller = pid_transfer(&pid);
        let loop_tf = series(&controller, &plant).unwrap();
        let one = TransferFunction::new(vec![1.0], vec![1.0]).unwrap();
        let closed = feedback(&loop_tf, &one).unwrap();
        assert!(closed.is_stable());
    }

    #[test]
    fn all_rules_produce_finite_positive_proportional_gain() {
        let plant = plant_third_order();
        for rule in [
            TuningRule::ZieglerNichols,
            TuningRule::CohenCoon,
            TuningRule::Imc,
            TuningRule::Amigo,
            TuningRule::ChienHronesReswick,
        ] {
            let pid = tune(&plant, rule, None).unwrap();
            assert!(pid.kp.is_finite() && pid.kp > 0.0, "rule {rule:?}");
        }
    }

    #[test]
    fn simulate_pid_clamps_output_and_does_not_wind_up() {
        let pid = Pid {
            kp: 1.0,
            ki: 1.0,
            kd: 0.0,
            tf: None,
        };
        let error = vec![1.0; 200];
        let u = simulate_pid(&pid, &error, 0.1, Some((-2.0, 2.0)));
        assert!(u.iter().all(|v| *v >= -2.0 - 1e-12 && *v <= 2.0 + 1e-12));
        // Conditional integration keeps the saturated output bounded at 2.
        assert!((u.last().copied().unwrap() - 2.0).abs() < 1e-9);
    }

    #[test]
    fn pid_response_tracks_a_constant_setpoint() {
        let plant = TransferFunction::new(vec![1.0], vec![1.0, 1.0]).unwrap();
        let pid = Pid {
            kp: 2.0,
            ki: 2.0,
            kd: 0.1,
            tf: Some(0.01),
        };
        let t: Vec<f64> = (0..501).map(|i| i as f64 * 0.02).collect();
        let y = pid_response(&pid, &plant, &t, 1.0).unwrap();
        assert!((y.last().copied().unwrap() - 1.0).abs() < 0.05);
    }
}
