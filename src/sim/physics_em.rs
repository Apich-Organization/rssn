//! Euler-type integrators for mechanical systems.
//!
//! # Euler-type time integrators
//!
//! Forward, midpoint, Heun, semi-implicit (symplectic) and backward Euler for
//! ODE systems and mechanical systems. Not a PDE module: spatial
//! discretisations and boundary conditions live in the solver modules
//! ([`crate::sim::structured`], [`crate::sim::physics_fem`], ...).

use rayon::prelude::*;
use serde::Deserialize;
use serde::Serialize;

use crate::kernels::par::for_chunks_mut;
use crate::kernels::reuse::f64_scratch;
use crate::sim::physics_rkm::OdeSystem;

/// Configuration for the Euler solver.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct EulerSolverConfig {
    /// The time step.
    pub dt: f64,
}

/// Time integration formula used by [`solve_explicit_euler`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ExplicitEulerMethod {
    /// First-order forward Euler with one right-hand-side evaluation per step.
    Forward,
    /// Second-order explicit midpoint with two evaluations per step.
    Midpoint,
    /// Second-order Heun (explicit trapezoid) with two evaluations per step.
    Heun,
}

/// Solves an ODE with a fixed Euler formula, ending exactly at the requested time.
///
/// The sign of `dt` must agree with the direction of travel. Times are located
/// on the nominal grid from the start, without accumulating rounding drift.
/// State updates use the represented interval between successive grid points.
/// The last step ends at `t_span.1`, shortening the interval if necessary.
/// A zero-length span returns the initial state without evaluating the RHS.
///
/// # Errors
///
/// Returns an error for invalid times, a zero or oppositely signed step,
/// mismatched dimensions, non-finite states or right-hand sides, unrepresentable
/// progress, or a history whose size cannot be represented or reserved.
pub fn solve_explicit_euler<S: OdeSystem>(
    system: &S,
    y0: &[f64],
    t_span: (f64, f64),
    dt: f64,
    method: ExplicitEulerMethod,
) -> Result<Vec<(f64, Vec<f64>)>, String> {
    explicit_history::<S, true>(system, y0, t_span, dt, method).map_err(str::to_owned)
}

/// Solves an ODE with first-order forward Euler and a complete step at each sample.
///
/// For an endpoint between grid points the final sample uses a full step. Use
/// [`solve_explicit_euler`] for a checked solve ending at the requested time.
pub fn solve_forward_euler<S: OdeSystem>(
    system: &S,
    y0: &[f64],
    t_span: (f64, f64),
    dt: f64,
) -> Vec<(f64, Vec<f64>)> {
    legacy_explicit_history(system, y0, t_span, dt, ExplicitEulerMethod::Forward)
}

/// Solves an ODE with second-order explicit midpoint and complete steps.
///
/// Use [`solve_explicit_euler`] to terminate exactly at the requested time.
pub fn solve_midpoint_euler<S: OdeSystem>(
    system: &S,
    y0: &[f64],
    t_span: (f64, f64),
    dt: f64,
) -> Vec<(f64, Vec<f64>)> {
    legacy_explicit_history(system, y0, t_span, dt, ExplicitEulerMethod::Midpoint)
}

/// Solves an ODE with second-order Heun and complete steps.
///
/// Use [`solve_explicit_euler`] to terminate exactly at the requested time.
pub fn solve_heun_euler<S: OdeSystem>(
    system: &S,
    y0: &[f64],
    t_span: (f64, f64),
    dt: f64,
) -> Vec<(f64, Vec<f64>)> {
    legacy_explicit_history(system, y0, t_span, dt, ExplicitEulerMethod::Heun)
}

fn legacy_explicit_history<S: OdeSystem>(
    system: &S,
    y0: &[f64],
    t_span: (f64, f64),
    dt: f64,
    method: ExplicitEulerMethod,
) -> Vec<(f64, Vec<f64>)> {
    explicit_history::<S, false>(system, y0, t_span, dt, method)
        .unwrap_or_else(|_| vec![(t_span.0, y0.to_vec())])
}

struct EulerSchedule {
    start: f64,
    end: f64,
    dt: f64,
    steps: usize,
}

// The raw solver and native accounting both use this upper bound.
pub(crate) fn explicit_euler_step_bound(
    span: (f64, f64),
    dt: f64,
) -> Result<usize, &'static str> {
    euler_step_bound::<true>(span, dt)
}

fn euler_step_bound<const CHECKED: bool>(
    (start, end): (f64, f64),
    dt: f64,
) -> Result<usize, &'static str> {
    if !start.is_finite() || !end.is_finite() || !dt.is_finite() || dt == 0.0 {
        return Err("Euler times and step must be finite and the step nonzero");
    }
    if (end > start && dt < 0.0) || (end < start && dt > 0.0) {
        return Err("Euler step must advance toward the endpoint");
    }
    let span = end - start;
    let mut count = if span.is_finite() {
        (span / dt).ceil()
    } else {
        // Difference can overflow for finite endpoints of opposite signs.
        (end / dt - start / dt).ceil()
    };
    if CHECKED && span != 0.0 {
        count = count.max(1.0);
    }
    // A grid index above 2^53 loses integer precision when converted to f64.
    if !count.is_finite() || count > 9_007_199_254_740_992.0 || count >= usize::MAX as f64 {
        return Err("Euler step count cannot be represented");
    }
    Ok(count as usize)
}

impl EulerSchedule {
    fn new<const CHECKED: bool>(
        span: (f64, f64),
        dt: f64,
    ) -> Result<Self, &'static str> {
        let steps = if CHECKED {
            explicit_euler_step_bound(span, dt)?
        } else {
            euler_step_bound::<false>(span, dt)?
        };
        Ok(Self {
            start: span.0,
            end: span.1,
            dt,
            steps,
        })
    }

    fn next<const CHECKED: bool>(
        &self,
        index: usize,
        t: f64,
    ) -> Result<(f64, f64), &'static str> {
        let mut next = if CHECKED {
            self.dt.mul_add(index as f64, self.start)
        } else {
            t + self.dt
        };
        let forward = self.dt > 0.0;
        if CHECKED
            && (index == self.steps
                || (forward && next >= self.end)
                || (!forward && next <= self.end))
        {
            next = self.end;
        }
        let h = if CHECKED { next - t } else { self.dt };
        if !next.is_finite() || !h.is_finite() || (forward && next <= t) || (!forward && next >= t)
        {
            return Err("Euler step does not advance the represented time");
        }
        Ok((next, h))
    }
}

const EULER_BLOCK: usize = 8192;
const FORWARD_EULER_GRAIN: usize = 4;
// Two-stage formulas enter the pool twice per step. At the smaller grain,
// those two barriers cost more than the coordinate work in both methods.
const TWO_STAGE_EULER_GRAIN: usize = 8;

// A quarter of the finite range leaves room for both Heun rates and the state update.
const EULER_SAFE_VALUE: f64 = f64::MAX / 4.0;

fn euler_values_within(
    values: &[f64],
    limit: f64,
) -> bool {
    // A full, branchless scan lets the finite/range check vectorize.
    values
        .iter()
        .fold(true, |valid, value| valid & (value.abs() <= limit))
}

fn finite_euler_state<const CHECKED: bool>(values: &[f64]) -> Result<(), &'static str> {
    if CHECKED && !euler_values_within(values, f64::MAX) {
        return Err("Euler produced a non-finite state or derivative");
    }
    Ok(())
}

fn euler_rate_bounded<const CHECKED: bool>(
    rate: &[f64],
    limit: f64,
) -> Result<bool, &'static str> {
    if !CHECKED || euler_values_within(rate, limit) {
        return Ok(true);
    }
    finite_euler_state::<CHECKED>(rate)?;
    Ok(false)
}

fn euler_add(
    y: &mut [f64],
    rate: &[f64],
    h: f64,
    exceptional: bool,
    grain: usize,
) {
    for_chunks_mut(y, EULER_BLOCK, grain, |block, out| {
        let start = block * EULER_BLOCK;
        let end = start + out.len();
        let pairs = out.iter_mut().zip(&rate[start..end]);
        if exceptional {
            for (value, &rhs) in pairs {
                let tentative = *value + h * rhs;
                *value = if tentative.is_finite() {
                    tentative
                } else {
                    h.mul_add(rhs, *value)
                };
            }
        } else {
            for (value, &rhs) in pairs {
                *value += h * rhs;
            }
        }
    });
}

fn euler_heun_add(
    y: &mut [f64],
    k1: &[f64],
    k2: &[f64],
    h: f64,
    exceptional: bool,
) {
    let half = 0.5 * h;
    for_chunks_mut(y, EULER_BLOCK, TWO_STAGE_EULER_GRAIN, |block, out| {
        let start = block * EULER_BLOCK;
        let range = start..start + out.len();
        let triples = out.iter_mut().zip(&k1[range.clone()]).zip(&k2[range]);
        if exceptional {
            for ((value, &a), &b) in triples {
                let tentative = *value + half * (a + b);
                *value = if half.is_normal() && tentative.is_finite() {
                    tentative
                } else {
                    h.mul_add(f64::midpoint(a, b), *value)
                };
            }
        } else {
            for ((value, &a), &b) in triples {
                *value += half * (a + b);
            }
        }
    });
}

fn euler_step<S: OdeSystem, const CHECKED: bool>(
    system: &S,
    method: ExplicitEulerMethod,
    y: &mut [f64],
    workspace: &mut [f64],
    (t, next, h): (f64, f64, f64),
) -> Result<(), &'static str> {
    let dim = y.len();
    let (k1, rest) = workspace.split_at_mut(dim);
    let state_bounded = !CHECKED || euler_values_within(y, EULER_SAFE_VALUE);
    let rate_limit = EULER_SAFE_VALUE / h.abs().max(1.0);
    system.eval(t, y, k1);
    let first_bounded = euler_rate_bounded::<CHECKED>(k1, rate_limit)?;
    let exceptional_first = !state_bounded || !first_bounded;

    let exceptional_result = match method {
        | ExplicitEulerMethod::Forward => {
            euler_add(y, k1, h, exceptional_first, FORWARD_EULER_GRAIN);
            exceptional_first
        },
        | ExplicitEulerMethod::Midpoint | ExplicitEulerMethod::Heun => {
            let (stage, k2) = rest.split_at_mut(dim);
            let half = 0.5 * h;
            let tiny_half = CHECKED && !half.is_normal();
            let heun = matches!(method, ExplicitEulerMethod::Heun);
            let stage_h = if heun { h } else { half };
            for_chunks_mut(stage, EULER_BLOCK, TWO_STAGE_EULER_GRAIN, |block, out| {
                let start = block * EULER_BLOCK;
                let range = start..start + out.len();
                let triples = out.iter_mut().zip(&y[range.clone()]).zip(&k1[range]);
                if tiny_half && !heun {
                    for ((value, &yi), &a) in triples {
                        *value = h.mul_add(0.5 * a, yi);
                    }
                } else {
                    for ((value, &yi), &a) in triples {
                        *value = stage_h.mul_add(a, yi);
                    }
                }
            });
            if exceptional_first {
                finite_euler_state::<CHECKED>(stage)?;
            }
            let stage_time = if heun {
                next
            } else if CHECKED {
                f64::midpoint(t, next)
            } else {
                0.5f64.mul_add(h, t)
            };
            system.eval(stage_time, stage, k2);
            let second_bounded = euler_rate_bounded::<CHECKED>(k2, rate_limit)?;

            let exceptional = !state_bounded || !second_bounded;
            if heun {
                let exceptional = exceptional || !first_bounded || tiny_half;
                euler_heun_add(y, k1, k2, h, exceptional);
                exceptional
            } else {
                euler_add(y, k2, h, exceptional, TWO_STAGE_EULER_GRAIN);
                exceptional
            }
        },
    };
    if exceptional_result {
        finite_euler_state::<CHECKED>(y)?;
    }
    Ok(())
}

fn explicit_history<S: OdeSystem, const CHECKED: bool>(
    system: &S,
    y0: &[f64],
    t_span: (f64, f64),
    dt: f64,
    method: ExplicitEulerMethod,
) -> Result<Vec<(f64, Vec<f64>)>, &'static str> {
    use std::mem::size_of;

    let schedule = EulerSchedule::new::<CHECKED>(t_span, dt)?;
    let dim = y0.len();
    if dim != system.dim() {
        return Err("Euler state length does not match the system dimension");
    }
    finite_euler_state::<CHECKED>(y0)?;
    if schedule.steps != 0 {
        schedule.next::<CHECKED>(1, schedule.start)?;
    }
    let rows = schedule.steps + 1;
    let record_bytes = dim
        .checked_mul(size_of::<f64>())
        .and_then(|bytes| bytes.checked_add(size_of::<(f64, Vec<f64>)>()))
        .and_then(|row| row.checked_mul(rows));
    if record_bytes.is_none_or(|size| size > isize::MAX as usize) {
        return Err("Euler history exceeds addressable storage");
    }

    let mut history = Vec::new();
    history
        .try_reserve_exact(rows)
        .map_err(|_| "Euler history reservation failed")?;
    history.push((schedule.start, y0.to_vec()));
    if schedule.steps == 0 {
        return Ok(history);
    }

    let buffers = if matches!(method, ExplicitEulerMethod::Forward) {
        1
    } else {
        3
    };
    let len = dim
        .checked_mul(buffers)
        .ok_or("Euler workspace size overflow")?;
    if len > isize::MAX as usize / size_of::<f64>() {
        return Err("Euler workspace exceeds addressable storage");
    }
    let mut workspace = f64_scratch(len);
    workspace.resize(len, 0.0);
    let mut t = schedule.start;
    let mut y = y0.to_vec();
    for i in 1..=schedule.steps {
        let (next, h) = match schedule.next::<CHECKED>(i, t) {
            | Ok(step) => step,
            | Err(_) if !CHECKED => break,
            | Err(error) => return Err(error),
        };
        euler_step::<S, CHECKED>(system, method, &mut y, &mut workspace, (t, next, h))?;
        history.push((next, y.clone()));
        t = next;
        if CHECKED && (schedule.end - t) == 0.0 {
            break;
        }
    }
    Ok(history)
}

/// Defines a special kind of ODE system for mechanics: dv/dt = a(x), dx/dt = v.
/// The state vector is partitioned into `[positions..., velocities...]`.
pub trait MechanicalSystem {
    /// The number of spatial dimensions (e.g., 1, 2, 3).
    fn spatial_dim(&self) -> usize;

    /// Evaluates the acceleration `a(x)` given the positions `x`.
    fn eval_acceleration(
        &self,
        x: &[f64],
        a: &mut [f64],
    );
}

/// Solves a second-order ODE system (typically mechanical systems) using the semi-implicit Euler method.
///
/// This method is often used in mechanics for its energy conservation properties (it is symplectic).
/// It is particularly well-suited for systems where the state vector can be partitioned into
/// positions and velocities, and accelerations depend only on positions.
/// It assumes the state vector `y` is ordered as `[x0, x1, ..., v0, v1, ...]`.
///
/// # Arguments
/// * `system` - The mechanical system to solve, implementing the `MechanicalSystem` trait.
/// * `y0` - The initial state `[initial_positions..., initial_velocities...]`.
/// * `dt` - The fixed time step.
///
/// # Returns
/// A `Vec` of tuples `(time, state_vector)` representing the solution path.
///
/// # Errors
///
/// This function will return an error if the length of the initial state vector `y0`
/// is not twice the spatial dimension of the system.
pub fn solve_semi_implicit_euler<S: MechanicalSystem>(
    system: &S,
    y0: &[f64],
    t_span: (f64, f64),
    dt: f64,
) -> Result<Vec<(f64, Vec<f64>)>, String> {
    let s_dim = system.spatial_dim();

    if y0.len() != 2 * s_dim {
        return Err("State vector \
                    length must be \
                    twice the spatial \
                    dimension."
            .to_string());
    }

    let (t_start, t_end) = t_span;

    let steps: usize = (((t_end - t_start) / dt).ceil() as i64)
        .try_into()
        .unwrap_or(0);

    let mut t = t_start;

    let mut y = y0.to_vec();

    let mut history = Vec::with_capacity(steps + 1);

    history.push((t, y.clone()));

    let mut a = vec![0.0; s_dim];

    for _ in 0..steps {
        let (x, v) = y.split_at_mut(s_dim);

        system.eval_acceleration(x, &mut a);

        v.par_iter_mut().zip(&a).for_each(|(vi, &ai)| {
            *vi += dt * ai;
        });

        x.par_iter_mut().zip(&*v).for_each(|(xi, &vi)| {
            *xi += dt * vi;
        });

        t += dt;

        history.push((t, y.clone()));
    }

    Ok(history)
}

use crate::sim::physics_rkm::DampedOscillatorSystem;

/// Solves the damped harmonic oscillator with the less stable forward Euler method.
///
/// This scenario demonstrates the numerical energy gain over time that can occur
/// with the forward Euler method when solving oscillatory systems, highlighting
/// its limitations for such problems.
///
/// # Returns
/// A `Vec` of tuples `(time, state_vector)` representing the solution path.
#[must_use]
pub fn simulate_oscillator_forward_euler_scenario() -> Vec<(f64, Vec<f64>)> {
    let system = DampedOscillatorSystem {
        omega: 2.0 * std::f64::consts::PI,
        zeta: 0.0,
    };

    let y0 = &[1.0, 0.0];

    let t_span = (0.0, 10.0);

    let dt = 0.01;

    solve_forward_euler(&system, y0, t_span, dt)
}

/// A simple 2D orbital system (e.g., planet around a star).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct OrbitalSystem {
    /// The gravitational constant.
    pub gravitational_constant: f64,
    /// The mass of the star.
    pub star_mass: f64,
}

impl MechanicalSystem for OrbitalSystem {
    fn spatial_dim(&self) -> usize {
        2
    }

    fn eval_acceleration(
        &self,
        x: &[f64],
        a: &mut [f64],
    ) {
        let (px, py) = (x[0], x[1]);

        let dist_sq = py.mul_add(py, px.powi(2));

        let dist_cubed = dist_sq.sqrt().powi(3);

        let force_magnitude = -self.gravitational_constant * self.star_mass / dist_cubed;

        a[0] = force_magnitude * px;

        a[1] = force_magnitude * py;
    }
}

/// Solves a 2D orbital mechanics problem using the energy-conserving semi-implicit Euler method.
///
/// This scenario simulates the orbit of a celestial body around a much more massive star,
/// demonstrating the semi-implicit Euler method's ability to preserve energy and produce
/// stable, realistic orbits over long simulation times.
///
/// # Returns
/// A `Vec` of tuples `(time, state_vector)` representing the solution path.
///
/// # Errors
///
/// This function will return an error if the underlying `solve_semi_implicit_euler`
/// function encounters an error, e.g., if the initial state vector dimensions are incorrect.
pub fn simulate_gravity_semi_implicit_euler_scenario() -> Result<Vec<(f64, Vec<f64>)>, String> {
    let system = OrbitalSystem {
        gravitational_constant: 1.0,
        star_mass: 1000.0,
    };

    let y0 = &[10.0, 0.0, 0.0, 30.0];

    let t_span = (0.0, 2.5);

    let dt = 0.001;

    solve_semi_implicit_euler(&system, y0, t_span, dt)
}

use crate::kernels::matrix::Matrix;

/// Defines a linear ODE system: dy/dt = A * y.
pub trait LinearOdeSystem {
    /// The dimension of the system.
    fn dim(&self) -> usize;

    /// Provides the system matrix A.
    fn get_matrix(&self) -> Matrix<f64>;
}

/// Solves a linear ODE system y' = Ay using the backward (implicit) Euler method.
/// This method is A-stable, making it excellent for stiff equations.
/// The update rule is y_{n+1} = (I - dt*A)^-1 * `y_n`.
///
/// # Returns
/// A `Result` containing the solution path, or an error string if matrix inversion fails.
/// Solves a linear ODE system `y' = Ay` using the backward (implicit) Euler method.
///
/// This method is A-stable, making it excellent for stiff equations where explicit methods
/// would require prohibitively small time steps. The update rule is `y_{n+1} = (I - dt*A)^-1 * y_n`.
///
/// # Arguments
/// * `system` - The linear ODE system to solve, implementing the `LinearOdeSystem` trait.
/// * `y0` - The initial state vector.
/// * `t_span` - The time interval `(t_start, t_end)`.
///
/// # Returns
/// A `Result` containing the solution path as `Vec<(f64, Vec<f64>)>`.
///
/// # Errors
///
/// This function will return an error if the matrix `(I - dt*A)` is not invertible,
/// which can occur if the system is singular for the given time step `dt`.
pub fn solve_backward_euler_linear<S: LinearOdeSystem>(
    system: &S,
    y0: &[f64],
    t_span: (f64, f64),
    dt: f64,
) -> Result<Vec<(f64, Vec<f64>)>, String> {
    let (t_start, t_end) = t_span;

    let steps: usize = (((t_end - t_start) / dt).ceil() as i64)
        .try_into()
        .unwrap_or(0);

    let mut t = t_start;

    let mut y = y0.to_vec();

    let mut history = Vec::with_capacity(steps + 1);

    history.push((t, y.clone()));

    let dim = system.dim();

    let a = system.get_matrix();

    let identity = Matrix::identity(dim);

    let m = identity - (a * dt);

    let m_inv = m.inverse().ok_or(
        "Matrix (I - dt*A) is not \
         invertible.",
    )?;

    for _ in 0..steps {
        let y_matrix = Matrix::new(dim, 1, y.clone());

        let y_new_matrix = m_inv.clone() * y_matrix;

        y.clone_from(&y_new_matrix.get_cols()[0]);

        t += dt;

        history.push((t, y.clone()));
    }

    Ok(history)
}

/// A stiff system where one component decays much faster than the other.
pub struct StiffDecaySystem;

impl LinearOdeSystem for StiffDecaySystem {
    fn dim(&self) -> usize {
        2
    }

    fn get_matrix(&self) -> Matrix<f64> {
        Matrix::new(2, 2, vec![-20.0, 0.0, 0.0, -0.5])
    }
}

/// Solves a stiff ODE system, demonstrating the stability of the backward Euler method.
///
/// This scenario highlights how the backward Euler method can handle stiff equations
/// with relatively large time steps, whereas a forward method would become unstable.
///
/// # Returns
/// A `Result` containing the solution path as `Vec<(f64, Vec<f64>)>`.
///
/// # Errors
///
/// This function will return an error if the underlying `solve_backward_euler_linear`
/// function encounters an error, e.g., if matrix inversion fails.
pub fn simulate_stiff_decay_scenario() -> Result<Vec<(f64, Vec<f64>)>, String> {
    let system = StiffDecaySystem;

    let y0 = &[1.0, 1.0];

    let t_span = (0.0, 5.0);

    let dt = 0.2;

    solve_backward_euler_linear(&system, y0, t_span, dt)
}
