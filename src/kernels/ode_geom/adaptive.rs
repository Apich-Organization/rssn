//! Error-controlled stepping for the requested geometric method, without method substitution.
use super::EventSpec;
use super::Events;
use super::FixedGeomSolution;
use super::GeometricMethod;
use super::Hit;
use super::OdeGeomError;
use super::Rate;
use super::evaluate;
use super::gauss_legendre_step;
use super::partitioned_ruth4_step;
use super::partitioned_verlet_step;
use std::cell::Cell;

/// Tolerances and hard work limits for adaptive geometric integration.
#[derive(Clone, Copy, Debug)]
pub struct AdaptiveGeomOptions {
    /// Relative componentwise local-error tolerance.
    pub rtol: f64,
    /// Positive absolute componentwise local-error tolerance.
    pub atol: f64,
    /// Initial trial step; zero selects one hundredth of the span.
    pub h0: f64,
    /// Minimum trial step, except a final clipped interval.
    pub h_min: f64,
    /// Maximum trial step (positive infinity is permitted).
    pub h_max: f64,
    /// Maximum attempted steps, including rejected steps.
    pub max_steps: usize,
    /// Maximum actual RHS callback invocations, including rejected trials.
    pub max_evaluations: usize,
}
impl Default for AdaptiveGeomOptions {
    fn default() -> Self {
        Self {
            rtol: 1e-8,
            atol: 1e-10,
            h0: 0.0,
            h_min: 1e-14,
            h_max: f64::INFINITY,
            max_steps: 100_000,
            max_evaluations: 10_000_000,
        }
    }
}
impl AdaptiveGeomOptions {
    pub(super) fn validate(self) -> Result<(), OdeGeomError> {
        if !self.rtol.is_finite()
            || self.rtol < 0.0
            || !self.atol.is_finite()
            || self.atol <= 0.0
            || !self.h0.is_finite()
            || self.h0 < 0.0
            || !self.h_min.is_finite()
            || self.h_min <= 0.0
            || self.h_max.is_nan()
            || self.h_max < self.h_min
            || self.max_steps == 0
            || self.max_evaluations == 0
        {
            return Err(OdeGeomError::Domain(
                "invalid adaptive tolerances, step bounds, or work limits".into(),
            ));
        }
        Ok(())
    }
}

pub(super) fn geometric_step(
    method: GeometricMethod,
    rhs: &mut Rate<'_>,
    t: f64,
    y: &[f64],
    h: f64,
) -> Result<Vec<f64>, OdeGeomError> {
    match method {
        | GeometricMethod::Rk4 => {
            let n = y.len();
            let k1 = evaluate(rhs, t, y, n)?;
            let point = |k: &[f64], fraction: f64| {
                y.iter()
                    .zip(k)
                    .map(|(a, b)| a + fraction * h * b)
                    .collect::<Vec<_>>()
            };
            let k2 = evaluate(rhs, t + 0.5 * h, &point(&k1, 0.5), n)?;
            let k3 = evaluate(rhs, t + 0.5 * h, &point(&k2, 0.5), n)?;
            let k4 = evaluate(rhs, t + h, &point(&k3, 1.0), n)?;
            let out: Vec<f64> = (0..n)
                .map(|i| y[i] + h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]))
                .collect();
            if out.iter().any(|x| !x.is_finite()) {
                return Err(OdeGeomError::NotConverged);
            }
            Ok(out)
        },
        | GeometricMethod::StormerVerlet => partitioned_verlet_step(rhs, t, y, h).map(|x| x.0),
        | GeometricMethod::Ruth4 => partitioned_ruth4_step(rhs, t, y, h).map(|x| x.0),
        | GeometricMethod::ImplicitMidpoint => gauss_legendre_step(rhs, t, y, h, 1).map(|x| x.0),
        | GeometricMethod::GaussLegendre4 => gauss_legendre_step(rhs, t, y, h, 2).map(|x| x.0),
        | GeometricMethod::GaussLegendre6 => gauss_legendre_step(rhs, t, y, h, 3).map(|x| x.0),
    }
}
const fn order(method: GeometricMethod) -> i32 {
    match method {
        | GeometricMethod::StormerVerlet | GeometricMethod::ImplicitMidpoint => 2,
        | GeometricMethod::GaussLegendre6 => 6,
        | _ => 4,
    }
}
fn hermite(
    y0: &[f64],
    y1: &[f64],
    f0: &[f64],
    f1: &[f64],
    h: f64,
    theta: f64,
) -> Vec<f64> {
    let t2 = theta * theta;
    let t3 = t2 * theta;
    (0..y0.len())
        .map(|i| {
            (2.0 * t3 - 3.0 * t2 + 1.0) * y0[i]
                + (t3 - 2.0 * t2 + theta) * h * f0[i]
                + (-2.0 * t3 + 3.0 * t2) * y1[i]
                + (t3 - t2) * h * f1[i]
        })
        .collect()
}
fn guard(
    events: &mut Events<'_>,
    index: usize,
    t: f64,
    y: &[f64],
    invalid: &Cell<bool>,
) -> Result<f64, OdeGeomError> {
    let value = events(index, t, y);
    if !value.is_finite() {
        invalid.set(true);
        return Err(OdeGeomError::Domain(
            "event callback returned a nonfinite value".into(),
        ));
    }
    Ok(value)
}
#[allow(clippy::too_many_arguments)] // One accepted interval and its callback contracts.
fn append_segment(
    solution: &mut FixedGeomSolution,
    previous: &mut [f64],
    t1: f64,
    y1: Vec<f64>,
    f1: Vec<f64>,
    specs: &[EventSpec],
    events: &mut Events<'_>,
    invalid: &Cell<bool>,
) -> Result<bool, OdeGeomError> {
    let t0 = *solution.t.last().ok_or(OdeGeomError::NotConverged)?;
    let y0 = solution.y.last().ok_or(OdeGeomError::NotConverged)?;
    let f0 = solution.rates.last().ok_or(OdeGeomError::NotConverged)?;
    let h = t1 - t0;
    let mut hits = Vec::new();
    for (index, spec) in specs.iter().enumerate() {
        let g0 = previous[index];
        let g1 = guard(events, index, t1, &y1, invalid)?;
        let rising = g0 < 0.0 && g1 >= 0.0;
        let falling = g0 > 0.0 && g1 <= 0.0;
        if (spec.direction >= 0 && rising) || (spec.direction <= 0 && falling) {
            let (mut lo, mut hi) = (0.0, 1.0);
            let mut gl = g0;
            if super::super::exact_eq(g1, 0.0) {
                lo = 1.0;
            } else {
                for _ in 0..60 {
                    let mid = f64::midpoint(lo, hi);
                    let ym = hermite(y0, &y1, f0, &f1, h, mid);
                    let gm = guard(events, index, t0 + mid * h, &ym, invalid)?;
                    if super::super::exact_eq(gm, 0.0) {
                        lo = mid;
                        hi = mid;
                        break;
                    }
                    if gm.is_sign_positive() == gl.is_sign_positive() {
                        lo = mid;
                        gl = gm;
                    } else {
                        hi = mid;
                    }
                }
            }
            let theta = f64::midpoint(lo, hi);
            let t = t0 + theta * h;
            hits.push(Hit {
                index,
                t,
                y: hermite(y0, &y1, f0, &f1, h, theta),
            });
        }
        previous[index] = g1;
    }
    hits.sort_by(|a, b| a.t.total_cmp(&b.t));
    for hit in hits {
        let terminal = specs[hit.index].terminal;
        if terminal {
            // Differentiate the same continuous extension so terminal dense
            // output retains its endpoint derivatives without another RHS call.
            let theta = (hit.t - t0) / h;
            let t2 = theta * theta;
            let derivative = (0..y0.len())
                .map(|i| {
                    ((6.0 * t2 - 6.0 * theta) * y0[i] + (-6.0 * t2 + 6.0 * theta) * y1[i]) / h
                        + (3.0 * t2 - 4.0 * theta + 1.0) * f0[i]
                        + (3.0 * t2 - 2.0 * theta) * f1[i]
                })
                .collect();
            solution.t.push(hit.t);
            solution.y.push(hit.y.clone());
            solution.rates.push(derivative);
            solution.hits.push(hit);
            return Ok(true);
        }
        solution.hits.push(hit);
    }
    solution.t.push(t1);
    solution.y.push(y1);
    solution.rates.push(f1);
    Ok(false)
}

/// Integrate the selected RK4 or geometric method using step doubling.
///
/// The accepted state is two half steps of the requested method, without
/// extrapolation or method substitution. Both the endpoint estimate and the
/// midpoint defect of its cubic Hermite extension control acceptance. Adaptive
/// state-dependent step selection does not retain fixed-step symplecticity.
///
/// # Errors
/// Invalid inputs or callback outputs, exhausted attempts/evaluations, a
/// tolerance requiring a smaller step than `h_min`, or failed implicit stages.
#[allow(clippy::too_many_arguments)] // Stable problem/options/events boundary.
pub fn solve_adaptive_geometric(
    rhs: &mut Rate<'_>,
    y0: &[f64],
    span: (f64, f64),
    method: GeometricMethod,
    options: AdaptiveGeomOptions,
    specs: &[EventSpec],
    events: &mut Events<'_>,
    invalid_event: &Cell<bool>,
) -> Result<FixedGeomSolution, OdeGeomError> {
    options.validate()?;
    let (start, end) = span;
    if y0.is_empty() {
        return Err(OdeGeomError::Empty);
    }
    if !start.is_finite()
        || !end.is_finite()
        || end < start
        || !(end - start).is_finite()
        || y0.iter().any(|x| !x.is_finite())
        || specs.iter().any(|s| !matches!(s.direction, -1..=1))
    {
        return Err(OdeGeomError::Domain("invalid adaptive ODE problem".into()));
    }
    if matches!(
        method,
        GeometricMethod::StormerVerlet | GeometricMethod::Ruth4
    ) && !y0.len().is_multiple_of(2)
    {
        return Err(OdeGeomError::ShapeMismatch(
            "partitioned geometric state must have even length".into(),
        ));
    }
    let used = Cell::new(0usize);
    let exhausted = Cell::new(false);
    let mut counted = |t: f64, y: &[f64]| {
        if exhausted.get() || used.get() >= options.max_evaluations {
            exhausted.set(true);
            return vec![f64::NAN; y.len()];
        }
        used.set(used.get() + 1);
        rhs(t, y)
    };
    let r0 = evaluate(&mut counted, start, y0, y0.len())?;
    let mut solution = FixedGeomSolution {
        t: vec![start],
        y: vec![y0.to_vec()],
        rates: vec![r0],
        hits: Vec::new(),
        nfev: 0,
        rejected: 0,
    };
    let mut previous = Vec::with_capacity(specs.len());
    for (index, spec) in specs.iter().enumerate() {
        let g = guard(events, index, start, y0, invalid_event)?;
        previous.push(g);
        if super::super::exact_eq(g, 0.0) {
            solution.hits.push(Hit {
                index,
                t: start,
                y: y0.to_vec(),
            });
            if spec.terminal {
                solution.nfev = used.get();
                return Ok(solution);
            }
        }
    }
    let mut t = start;
    let mut y = y0.to_vec();
    let mut h = if options.h0 > 0.0 {
        options.h0
    } else {
        (end - start) / 100.0
    }
    .min(options.h_max)
    .max(options.h_min);
    let exponent = 1.0 / f64::from((order(method) + 1).min(4));
    let divisor = 2.0_f64.powi(order(method)) - 1.0;
    for _ in 0..options.max_steps {
        if t >= end {
            solution.nfev = used.get();
            return Ok(solution);
        }
        let remaining = end - t;
        h = h.min(remaining).min(options.h_max);
        if crate::kernels::exact_eq(t + h, t)
            || !h.is_finite()
            || (h < options.h_min && remaining > options.h_min)
        {
            return Err(OdeGeomError::NotConverged);
        }
        let next = if h >= remaining { end } else { t + h };
        let mid = f64::midpoint(t, next);
        let half = mid - t;
        let trial = (|| {
            let coarse = geometric_step(method, &mut counted, t, &y, next - t)?;
            let middle = geometric_step(method, &mut counted, t, &y, half)?;
            let fine = geometric_step(method, &mut counted, mid, &middle, next - mid)?;
            let fm = evaluate(&mut counted, mid, &middle, y.len())?;
            let ff = evaluate(&mut counted, next, &fine, y.len())?;
            Ok::<_, OdeGeomError>((coarse, middle, fine, fm, ff))
        })();
        if exhausted.get() {
            return Err(OdeGeomError::NotConverged);
        }
        let (coarse, middle, fine, fm, ff) = match trial {
            | Ok(values) => values,
            | Err(OdeGeomError::NotConverged | OdeGeomError::Singular) => {
                solution.rejected += 1;
                h *= 0.5;
                continue;
            },
            | Err(error) => return Err(error),
        };
        let f0 = solution.rates.last().ok_or(OdeGeomError::NotConverged)?;
        let dense_mid = hermite(&y, &fine, f0, &ff, next - t, 0.5);
        let mut error = 0.0_f64;
        for i in 0..y.len() {
            let scale =
                options.atol + options.rtol * y[i].abs().max(fine[i].abs()).max(middle[i].abs());
            error = error
                .max(((fine[i] - coarse[i]) / divisor / scale).abs())
                .max(((middle[i] - dense_mid[i]) / scale).abs());
        }
        if !error.is_finite() {
            return Err(OdeGeomError::NotConverged);
        }
        if error <= 1.0 {
            if append_segment(
                &mut solution,
                &mut previous,
                mid,
                middle,
                fm,
                specs,
                events,
                invalid_event,
            )? || append_segment(
                &mut solution,
                &mut previous,
                next,
                fine.clone(),
                ff,
                specs,
                events,
                invalid_event,
            )? {
                solution.nfev = used.get();
                return Ok(solution);
            }
            t = next;
            y = fine;
            h *= if super::super::exact_eq(error, 0.0) {
                3.0
            } else {
                (0.9 * error.powf(-exponent)).clamp(0.2, 3.0)
            };
        } else {
            solution.rejected += 1;
            h *= (0.9 * error.powf(-exponent)).clamp(0.1, 0.8);
        }
    }
    if t >= end {
        solution.nfev = used.get();
        Ok(solution)
    } else {
        Err(OdeGeomError::NotConverged)
    }
}

/// Unified geometric ODE solver supporting both fixed and adaptive step sizes
/// via [`OdeStepMode`](crate::kernels::ode::OdeStepMode).
///
/// # Errors
/// As [`solve_fixed_geometric`](super::solve_fixed_geometric) and
/// [`solve_adaptive_geometric`].
#[allow(clippy::too_many_arguments)]
pub fn solve_geometric(
    rhs: &mut Rate<'_>,
    y0: &[f64],
    span: (f64, f64),
    method: GeometricMethod,
    step_mode: crate::kernels::ode::OdeStepMode,
    specs: &[EventSpec],
    events: &mut Events<'_>,
    invalid_event: &Cell<bool>,
) -> Result<FixedGeomSolution, OdeGeomError> {
    match step_mode {
        | crate::kernels::ode::OdeStepMode::Fixed { dt } => {
            if !dt.is_finite() || dt <= 0.0 || !span.0.is_finite() || !span.1.is_finite() {
                return Err(OdeGeomError::Domain("invalid fixed step size dt".into()));
            }
            let steps = ((span.1 - span.0).abs() / dt).ceil().max(1.0) as usize;
            super::solve_fixed_geometric(
                rhs,
                y0,
                span,
                steps,
                method,
                None,
                specs,
                events,
                invalid_event,
            )
        },
        | crate::kernels::ode::OdeStepMode::Adaptive {
            rtol,
            atol,
            dt_init,
            dt_min,
            dt_max,
            max_steps,
        } => {
            let options = AdaptiveGeomOptions {
                rtol,
                atol,
                h0: dt_init.unwrap_or(0.0),
                h_min: dt_min.unwrap_or(1e-14),
                h_max: dt_max.unwrap_or(f64::INFINITY),
                max_steps,
                max_evaluations: max_steps.saturating_mul(100),
            };
            solve_adaptive_geometric(rhs, y0, span, method, options, specs, events, invalid_event)
        },
    }
}
