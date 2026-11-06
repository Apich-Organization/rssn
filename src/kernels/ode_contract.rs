//! Shared contracts and scale-safe startup for adaptive ODE solvers.
//!
//! No callback is evaluated while malformed problem data is being validated.
use super::ode_adaptive::Event;
use super::ode_adaptive::OdeError;
use super::ode_adaptive::OdeOptions;

/// Checks common adaptive inputs before a RHS, event, mass or Jacobian callback.
pub(super) fn validate_problem(
    t0: f64,
    t1: f64,
    state: &[f64],
    options: &OdeOptions,
    events: &[Event<'_>],
) -> Result<(), OdeError> {
    if !t0.is_finite()
        || !t1.is_finite()
        || super::exact_eq(t0, t1)
        || state.is_empty()
        || state.iter().any(|value| !value.is_finite())
        || !options.rtol.is_finite()
        || !options.atol.is_finite()
        || options.rtol < 0.0
        || options.atol < 0.0
        || (super::exact_eq(options.rtol, 0.0) && super::exact_eq(options.atol, 0.0))
        || options
            .h0
            .is_some_and(|step| !step.is_finite() || super::exact_eq(step, 0.0))
        || options.hmax.is_nan()
        || options.hmax <= 0.0
        || options.max_steps == 0
        || events
            .iter()
            .any(|event| !matches!(event.direction, -1..=1))
    {
        return Err(OdeError::Invalid);
    }
    Ok(())
}

/// A RHS must write every output component at every invocation. Reusing zero
/// or stale entries would otherwise turn an incomplete callback into a valid
/// differential equation. Auto-traits follow F, including Sync for parallel FD.
pub(super) fn checked_rhs<F: Fn(f64, &[f64], &mut [f64])>(
    function: F
) -> impl Fn(f64, &[f64], &mut [f64]) {
    checked_rhs_with_stop(function, || false)
}

/// The optional host cancellation check inherits auto-traits from its closure,
/// preserving Sync for the existing parallel finite-difference entry point.
pub(super) fn checked_rhs_with_stop<F, C>(
    function: F,
    stop: C,
) -> impl Fn(f64, &[f64], &mut [f64])
where
    F: Fn(f64, &[f64], &mut [f64]),
    C: Fn() -> bool,
{
    move |time: f64, state: &[f64], derivative: &mut [f64]| {
        derivative.fill(f64::NAN);
        if !stop() && time.is_finite() && state.iter().all(|value| value.is_finite()) {
            function(time, state, derivative);
        }
    }
}

/// A zero absolute scale with a zero error means no error; a nonzero error
/// at that scale cannot satisfy a relative-only tolerance.
#[allow(clippy::float_cmp)] // Exact zero defines the zero-tolerance case.
pub(super) fn scaled_ratio(
    value: f64,
    scale: f64,
) -> f64 {
    if scale == 0.0 {
        if value == 0.0 {
            0.0
        } else {
            f64::INFINITY
        }
    } else {
        value / scale
    }
}

/// RMS without squaring huge/small unscaled components or overflowing a norm
/// that only becomes representable after dividing by sqrt(count).
pub(super) fn rms(
    values: impl IntoIterator<Item = f64>,
    count: usize,
) -> f64 {
    let (mut largest, mut sum) = (0.0_f64, 0.0_f64);
    for value in values {
        let magnitude = value.abs();
        if !magnitude.is_finite() {
            return magnitude;
        }
        if magnitude > largest {
            sum = 1.0 + sum * (largest / magnitude).powi(2);
            largest = magnitude;
        } else if magnitude > 0.0 {
            sum += (magnitude / largest).powi(2);
        }
    }
    largest * (sum / count.max(1) as f64).sqrt()
}

#[derive(Clone, Copy)]
pub(super) struct InitialStepOptions {
    pub(super) rtol: f64,
    pub(super) atol: f64,
    pub(super) bound: f64,
    pub(super) max_step: f64,
    pub(super) order: f64,
}

/// The startup probe is a real RHS evaluation, included in public counters.
#[derive(Clone, Copy)]
pub(super) struct StepEstimate {
    pub(super) step: f64,
    pub(super) evaluations: usize,
}

/// Startup probes remain inside the requested interval and maximal step.
/// The estimate is a heuristic: the integrator still validates every step.
pub(super) fn initial_step<F: Fn(f64, &[f64], &mut [f64])>(
    function: &F,
    time: f64,
    state: &[f64],
    derivative: &[f64],
    options: InitialStepOptions,
) -> StepEstimate {
    let InitialStepOptions {
        rtol,
        atol,
        bound,
        max_step,
        order,
    } = options;
    let direction = if bound > time { 1.0 } else { -1.0 };
    let scales: Vec<f64> = state
        .iter()
        .map(|value| atol + rtol * value.abs())
        .collect();
    let d0 = rms(
        state
            .iter()
            .zip(&scales)
            .map(|(&value, &scale)| scaled_ratio(value, scale)),
        state.len(),
    );
    let d1 = rms(
        derivative
            .iter()
            .zip(&scales)
            .map(|(&value, &scale)| scaled_ratio(value, scale)),
        state.len(),
    );
    let estimate = if d0 < 1e-5 || d1 < 1e-5 || !d0.is_finite() || !d1.is_finite() {
        1e-6
    } else {
        0.01 * d0 / d1
    };
    // A step must be resolvable at the current time: far from the origin the
    // spacing of representable times can exceed the heuristic estimate, and a
    // step that cannot advance `time` would be reported as a step-limit failure
    // before the first step is taken.
    let resolution = 8.0 * (time.abs().next_up() - time.abs());
    let step = estimate
        .max(resolution)
        .min(max_step)
        .min((bound - time).abs());
    if !step.is_finite() || step <= 0.0 {
        return StepEstimate { step, evaluations: 0 };
    }
    let trial_time = if direction > 0.0 {
        (time + step).min(bound)
    } else {
        (time - step).max(bound)
    };
    if super::exact_eq(trial_time, time) {
        return StepEstimate { step, evaluations: 0 };
    }
    let step = (trial_time - time).abs();
    let trial: Vec<f64> = state
        .iter()
        .zip(derivative)
        .map(|(&value, &slope)| value + direction * step * slope)
        .collect();
    let mut trial_derivative = vec![f64::NAN; state.len()];
    function(trial_time, &trial, &mut trial_derivative);
    let d2 = rms(
        trial_derivative
            .iter()
            .zip(derivative)
            .zip(&scales)
            .map(|((&a, &b), &scale)| scaled_ratio(a - b, scale)),
        state.len(),
    ) / step;
    let curvature = if d2.is_finite() {
        d1.max(d2)
    } else {
        f64::INFINITY
    };
    let adjusted = if !curvature.is_finite() {
        let smaller = step * 1e-3;
        if smaller > 0.0 {
            smaller
        } else {
            step
        }
    } else if curvature <= 1e-15 {
        (step * 1e-3).max(1e-6)
    } else {
        (0.01 / curvature).powf(1.0 / (order + 1.0))
    };
    StepEstimate {
        step: (100.0 * step)
            .min(adjusted)
            .min(max_step)
            .min((bound - time).abs()),
        evaluations: 1,
    }
}

/// A finite-difference time probe stays within the model's declared interval.
/// At the upper endpoint a backward difference replaces the forward one.
pub(super) fn time_probe(
    time: f64,
    t0: f64,
    t1: f64,
) -> f64 {
    let lo = t0.min(t1);
    let hi = t0.max(t1);
    let delta = 1.5e-8 * time.abs().max(1e-3);
    let forward = (time + delta).min(hi);
    if forward > time {
        forward
    } else {
        (time - delta).max(lo)
    }
}

/// The event values at the initial sample must be valid as well as later ones.
pub(super) fn initial_events(
    events: &[Event<'_>],
    time: f64,
    state: &[f64],
) -> Result<Vec<f64>, OdeError> {
    events
        .iter()
        .map(|event| {
            let value = (event.g)(time, state);
            if value.is_finite() {
                Ok(value)
            } else {
                Err(OdeError::NonFinite)
            }
        })
        .collect()
}

/// A fixed grid is anchored to the original interval, including the exact
/// endpoint. Reject a substep that has no representable time advance.
pub(super) fn fixed_step(
    t0: f64,
    t1: f64,
    index: usize,
    count: usize,
    current: f64,
) -> Result<(f64, f64), OdeError> {
    let next = if index == count {
        t1
    } else {
        let fraction = index as f64 / count as f64;
        let width = t1 - t0;
        if width.is_finite() {
            t0 + fraction * width
        } else {
            (1.0 - fraction) * t0 + fraction * t1
        }
    };
    let step = next - current;
    if !next.is_finite()
        || !step.is_finite()
        || super::exact_eq(next, current)
        || (step > 0.0) != (t1 > t0)
    {
        return Err(OdeError::StepLimit);
    }
    Ok((next, step))
}

/// Internal host adapters keep the original error in their own slot. This
/// escape stops adaptive retries after that failure, without panic unwinding
/// or changing the public infallible-RHS API.
pub(super) fn abort_if_requested(stop: Option<&dyn Fn() -> bool>) -> Result<(), OdeError> {
    if stop.is_some_and(|stop| stop()) {
        Err(OdeError::NonFinite)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rms_preserves_representable_huge_and_tiny_norms() {
        for magnitude in [1e300, 1e-300] {
            let norm = rms([magnitude, magnitude], 2);
            assert!((norm / magnitude - 1.0).abs() < 1e-14);
        }
        let norm = rms([1e308, 1e308, 1e308, 1e308], 4);
        assert!(norm.is_finite());
        assert!((norm / 1e308 - 1.0).abs() < 1e-14);
    }

    #[test]
    fn relative_only_zero_scales_have_defined_error_ratios() {
        assert!(super::super::exact_eq(scaled_ratio(0.0, 0.0), 0.0));
        assert!(scaled_ratio(1.0, 0.0).is_infinite());
        assert!(scaled_ratio(f64::NAN, 1.0).is_nan());
    }
}
