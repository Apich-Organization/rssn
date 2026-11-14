//! Tolerance-controlled time evolution of semidiscrete PDEs.
//!
//! Spatial truncation error is independent of this temporal error estimate.
//! Every accepted state is produced by two half steps of the selected method;
//! the controller never substitutes another integrator or extrapolates states.
use super::Boundary3D;
use super::BoundaryCondition;
use super::FdmError;
use super::Grid3D;
use super::TimeStepping;
use super::enforce_dirichlet_3d;
use super::fixed_value_3d;
use super::homogeneous_3d;
use super::laplacian_3d;
use super::solve_line_implicit;
use super::upwind_delta_3d;

/// Adaptive time controls; `atol` has the dimension of a field value.
#[derive(Clone, Copy, Debug)]
pub struct PdeAdaptiveOptions {
    /// Initial time.
    pub t0: f64,
    /// Requested final time, at or after `t0`.
    pub t1: f64,
    /// Relative local temporal error tolerance.
    pub rtol: f64,
    /// Absolute local temporal error tolerance.
    pub atol: f64,
    /// Initial trial step; zero selects one tenth of the interval.
    pub h0: f64,
    /// Smallest ordinary trial step; the final remainder may be smaller.
    pub h_min: f64,
    /// Largest trial step.
    pub h_max: f64,
    /// Maximum accepted steps.
    pub max_steps: usize,
    /// Maximum calls to the selected step kernel, including event localization.
    pub max_step_calls: usize,
    /// Maximum rejected trials, including failed implicit solves.
    pub max_rejections: usize,
}
impl Default for PdeAdaptiveOptions {
    fn default() -> Self {
        Self {
            t0: 0.0,
            t1: 1.0,
            rtol: 1e-6,
            atol: 1e-9,
            h0: 0.0,
            h_min: 1e-12,
            h_max: 1.0,
            max_steps: 100_000,
            max_step_calls: 300_000,
            max_rejections: 1000,
        }
    }
}

/// Checked scalar terminal-event callback for PDE states.
pub type PdeEventFunction<'a> = dyn FnMut(f64, &[f64]) -> Result<f64, FdmError> + 'a;

/// A terminal event evaluated on accepted PDE states.
pub struct PdeTerminalEvent<'a> {
    /// Positive for rising crossings, negative for falling, zero for either.
    pub direction: i8,
    /// Finite scalar event function; its zero terminates integration.
    pub evaluate: &'a mut PdeEventFunction<'a>,
}

/// Accepted nodes and actual work of a tolerance-controlled PDE evolution.
#[derive(Clone, Debug)]
pub struct PdeAdaptiveSolution {
    /// Accepted times, including the exact final endpoint or terminal event.
    pub times: Vec<f64>,
    /// Flat states at the corresponding times.
    pub states: Vec<Vec<f64>>,
    /// Accepted steps (event-truncated steps count once).
    pub accepted: usize,
    /// Rejected trials; rejected states never appear in `states`.
    pub rejected: usize,
    /// Calls of the selected step kernel, including event localization.
    pub step_calls: usize,
    /// Final accepted normalized error estimate.
    pub error_ratio: f64,
    /// Terminal event time when integration stopped at a crossing.
    pub event_time: Option<f64>,
}

fn validate_options(options: &PdeAdaptiveOptions) -> Result<(), FdmError> {
    if !options.t0.is_finite()
        || !options.t1.is_finite()
        || options.t1 < options.t0
        || !(options.t1 - options.t0).is_finite()
        || !options.rtol.is_finite()
        || options.rtol < 0.0
        || !options.atol.is_finite()
        || options.atol <= 0.0
        || !options.h0.is_finite()
        || options.h0 < 0.0
        || !options.h_min.is_finite()
        || options.h_min <= 0.0
        || !options.h_max.is_finite()
        || options.h_max < options.h_min
        || options.max_step_calls == 0
        || options.max_steps == 0
        || options.max_rejections == 0
    {
        return Err(FdmError::Domain(
            "invalid adaptive PDE interval, tolerances, step bounds, or budgets".into(),
        ));
    }
    Ok(())
}
fn checked_state(
    value: Vec<f64>,
    expected: usize,
) -> Result<Vec<f64>, FdmError> {
    if value.len() != expected {
        return Err(FdmError::ShapeMismatch(
            "PDE step changed state length".into(),
        ));
    }
    if value.iter().any(|v| !v.is_finite()) {
        return Err(FdmError::NotConverged);
    }
    Ok(value)
}
fn take_step<S>(
    step: &mut S,
    t: f64,
    y: &[f64],
    h: f64,
    calls: &mut usize,
    limit: usize,
) -> Result<Vec<f64>, FdmError>
where
    S: FnMut(f64, &[f64], f64) -> Result<Vec<f64>, FdmError>,
{
    if *calls >= limit {
        return Err(FdmError::Domain(
            "adaptive PDE step-call budget exhausted".into(),
        ));
    }
    *calls += 1;
    checked_state(step(t, y, h)?, y.len())
}
fn error_ratio(
    coarse: &[f64],
    fine: &[f64],
    old: &[f64],
    atols: &[f64],
    rtol: f64,
    divisor: f64,
) -> f64 {
    coarse
        .iter()
        .zip(fine)
        .zip(old)
        .zip(atols)
        .map(|(((coarse, fine), old), atol)| {
            let value =
                ((fine - coarse) / divisor).abs() / (atol + rtol * old.abs().max(fine.abs()));
            if value.is_finite() {
                value
            } else {
                f64::INFINITY
            }
        })
        .fold(0.0, f64::max)
}
fn checked_event(
    event: &mut PdeTerminalEvent<'_>,
    t: f64,
    y: &[f64],
) -> Result<f64, FdmError> {
    let value = (event.evaluate)(t, y)?;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(FdmError::Domain(
            "PDE event returned a non-finite value".into(),
        ))
    }
}

/// Adapt any restartable PDE step kernel with a step-doubling error estimate.
///
/// `order` is the selected method's temporal order. `absolute_tolerances` has
/// one positive entry per state coefficient. `stability` returns the maximum
/// stable step for the current state, or positive infinity for an implicit
/// scheme. Temporal control cannot certify spatial discretization error.
///
/// # Errors
/// Invalid controls or state shape, exhausted budgets, underflow, a failed step
/// kernel, or a non-finite event. Tolerance violations are never force-accepted.
#[allow(clippy::too_many_arguments)] // Explicit independent controller contracts.
pub fn adaptive_evolve<S, C>(
    initial: &[f64],
    order: u32,
    absolute_tolerances: &[f64],
    options: &PdeAdaptiveOptions,
    mut step: S,
    mut stability: C,
    mut event: Option<&mut PdeTerminalEvent<'_>>,
) -> Result<PdeAdaptiveSolution, FdmError>
where
    S: FnMut(f64, &[f64], f64) -> Result<Vec<f64>, FdmError>,
    C: FnMut(f64, &[f64]) -> Result<f64, FdmError>,
{
    validate_options(options)?;
    if initial.is_empty()
        || initial.len() != absolute_tolerances.len()
        || order == 0
        || order > 8
        || initial.iter().any(|v| !v.is_finite())
        || absolute_tolerances
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
    {
        return Err(FdmError::Domain(
            "invalid adaptive PDE initial state or component tolerances".into(),
        ));
    }
    let mut solution = PdeAdaptiveSolution {
        times: vec![options.t0],
        states: vec![initial.to_vec()],
        accepted: 0,
        rejected: 0,
        step_calls: 0,
        error_ratio: 0.0,
        event_time: None,
    };
    let mut t = options.t0;
    let mut y = initial.to_vec();
    let mut previous_event = match event.as_deref_mut() {
        | Some(e) => Some(checked_event(e, t, &y)?),
        | None => None,
    };
    if previous_event.is_some_and(|v| crate::kernels::exact_eq(v, 0.0)) {
        solution.event_time = Some(t);
        return Ok(solution);
    }
    let mut h = if options.h0 > 0.0 {
        options.h0
    } else {
        (options.t1 - t) / 10.0
    };
    h = h.clamp(options.h_min, options.h_max);
    let divisor = 2.0_f64.powi(order as i32) - 1.0;
    while t < options.t1 {
        if solution.accepted >= options.max_steps {
            return Err(FdmError::NotConverged);
        }
        let cap = stability(t, &y)?;
        if cap.is_nan() || cap <= 0.0 {
            return Err(FdmError::Domain(
                "PDE stability bound must be positive".into(),
            ));
        }
        let remaining = options.t1 - t;
        let trial = h.min(options.h_max).min(cap).min(remaining);
        if (trial < options.h_min && remaining > options.h_min) || t + trial <= t {
            return Err(FdmError::NotConverged);
        }
        // Count attempted kernels even when an implicit solve fails.
        let coarse = take_step(
            &mut step,
            t,
            &y,
            trial,
            &mut solution.step_calls,
            options.max_step_calls,
        );
        let candidate = match coarse {
            | Ok(coarse) => {
                let half = take_step(
                    &mut step,
                    t,
                    &y,
                    0.5 * trial,
                    &mut solution.step_calls,
                    options.max_step_calls,
                );
                match half {
                    | Ok(half) => {
                        take_step(
                            &mut step,
                            t + 0.5 * trial,
                            &half,
                            0.5 * trial,
                            &mut solution.step_calls,
                            options.max_step_calls,
                        )
                        .map(|fine| (coarse, fine))
                    },
                    | Err(error) => Err(error),
                }
            },
            | Err(error) => Err(error),
        };
        let (coarse, fine) = match candidate {
            | Ok(pair) => pair,
            | Err(FdmError::NotConverged) => {
                solution.rejected += 1;
                if solution.rejected > options.max_rejections || trial <= options.h_min {
                    return Err(FdmError::NotConverged);
                }
                h = (0.25 * trial).max(options.h_min);
                continue;
            },
            | Err(error) => return Err(error),
        };
        let ratio = error_ratio(
            &coarse,
            &fine,
            &y,
            absolute_tolerances,
            options.rtol,
            divisor,
        );
        let factor = if ratio <= f64::MIN_POSITIVE {
            4.0
        } else {
            (0.9 * ratio.powf(-1.0 / f64::from(order + 1))).clamp(0.1, 4.0)
        };
        if !ratio.is_finite() || ratio > 1.0 {
            solution.rejected += 1;
            if solution.rejected > options.max_rejections || trial <= options.h_min {
                return Err(FdmError::NotConverged);
            }
            h = (trial * factor.min(0.8)).max(options.h_min);
            continue;
        }
        let endpoint = if crate::kernels::exact_eq(trial, remaining) {
            options.t1
        } else {
            t + trial
        };
        if let Some(e) = event.as_deref_mut() {
            let next = checked_event(e, endpoint, &fine)?;
            let previous = previous_event.unwrap_or(next);
            let crossing = (previous < 0.0 && next >= 0.0 && e.direction >= 0)
                || (previous > 0.0 && next <= 0.0 && e.direction <= 0);
            if crossing {
                let (mut lo, mut hi) = (0.0, trial);
                let mut left = previous;
                let mut found = fine.clone();
                let mut found_ratio = ratio;
                for _ in 0..64 {
                    let mid = f64::midpoint(lo, hi);
                    if t + mid <= t || hi - lo <= 8.0 * f64::EPSILON * endpoint.abs().max(1.0) {
                        break;
                    }
                    let coarse_probe = take_step(
                        &mut step,
                        t,
                        &y,
                        mid,
                        &mut solution.step_calls,
                        options.max_step_calls,
                    )?;
                    let half = take_step(
                        &mut step,
                        t,
                        &y,
                        0.5 * mid,
                        &mut solution.step_calls,
                        options.max_step_calls,
                    )?;
                    let state = take_step(
                        &mut step,
                        t + 0.5 * mid,
                        &half,
                        0.5 * mid,
                        &mut solution.step_calls,
                        options.max_step_calls,
                    )?;
                    let probe_ratio = error_ratio(
                        &coarse_probe,
                        &state,
                        &y,
                        absolute_tolerances,
                        options.rtol,
                        divisor,
                    );
                    if !probe_ratio.is_finite() || probe_ratio > 1.0 {
                        return Err(FdmError::NotConverged);
                    }
                    let value = checked_event(e, t + mid, &state)?;
                    if crate::kernels::exact_eq(value, 0.0) {
                        hi = mid;
                        found = state;
                        found_ratio = probe_ratio;
                        break;
                    }
                    if (left < 0.0) == (value < 0.0) {
                        lo = mid;
                        left = value;
                    } else {
                        hi = mid;
                        found = state;
                        found_ratio = probe_ratio;
                    }
                }
                solution.accepted += 1;
                solution.error_ratio = found_ratio;
                solution.event_time = Some(t + hi);
                solution.times.push(t + hi);
                solution.states.push(found);
                return Ok(solution);
            }
            previous_event = Some(next);
        }
        t = endpoint;
        y = fine;
        solution.accepted += 1;
        solution.error_ratio = ratio;
        solution.times.push(t);
        solution.states.push(y.clone());
        h = (trial * factor).min(options.h_max);
    }
    Ok(solution)
}

fn geometry(
    extents: &[usize],
    spacing: f64,
    bc: &Boundary3D,
) -> Result<(Grid3D, Boundary3D), FdmError> {
    if !(1..=3).contains(&extents.len())
        || extents.iter().any(|n| *n < 3)
        || !spacing.is_finite()
        || spacing <= 0.0
    {
        return Err(FdmError::Domain("adaptive PDE grids require 1-3 axes, at least three points each, and positive finite spacing".into()));
    }
    let mut dims = [1; 3];
    dims[..extents.len()].copy_from_slice(extents);
    dims.iter()
        .try_fold(1usize, |n, d| n.checked_mul(*d))
        .ok_or_else(|| FdmError::Domain("PDE grid size overflows".into()))?;
    let mut conditions = bc.clone();
    if extents.len() < 3 {
        conditions.z_min = BoundaryCondition::Periodic;
        conditions.z_max = BoundaryCondition::Periodic;
    }
    if extents.len() < 2 {
        conditions.y_min = BoundaryCondition::Periodic;
        conditions.y_max = BoundaryCondition::Periodic;
    }
    conditions.validate()?;
    let grid = Grid3D {
        nx: dims[0],
        ny: dims[1],
        nz: dims[2],
        h: spacing,
    };
    for (face, count) in [
        (&conditions.x_min, dims[1] * dims[2]),
        (&conditions.x_max, dims[1] * dims[2]),
        (&conditions.y_min, dims[0] * dims[2]),
        (&conditions.y_max, dims[0] * dims[2]),
        (&conditions.z_min, dims[0] * dims[1]),
        (&conditions.z_max, dims[0] * dims[1]),
    ] {
        let finite = match face {
            | BoundaryCondition::Dirichlet(v) | BoundaryCondition::Neumann(v) => v.is_finite(),
            | BoundaryCondition::DirichletProfile(v) => v.iter().all(|v| v.is_finite()),
            | BoundaryCondition::Robin { alpha, beta, gamma } => {
                alpha.is_finite() && beta.is_finite() && gamma.is_finite()
            },
            | BoundaryCondition::Periodic => true,
        };
        if !finite {
            return Err(FdmError::Domain(
                "adaptive PDE boundary data must be finite".into(),
            ));
        }
        if let BoundaryCondition::DirichletProfile(values) = face
            && values.len() != count
        {
            return Err(FdmError::ShapeMismatch(
                "Dirichlet profile length differs from its boundary face".into(),
            ));
        }
    }
    Ok((grid, conditions))
}
fn initial_field(
    initial: &[f64],
    grid: Grid3D,
    bc: &Boundary3D,
) -> Result<Vec<f64>, FdmError> {
    if initial.len() != grid.len() {
        return Err(FdmError::ShapeMismatch(
            "initial field does not fit PDE grid".into(),
        ));
    }
    let mut field = initial.to_vec();
    for z in 0..grid.nz {
        for y in 0..grid.ny {
            for x in 0..grid.nx {
                if let Some(value) = fixed_value_3d(grid, bc, x, y, z) {
                    let actual = initial[grid.index(x, y, z)];
                    if (actual - value).abs() > 64.0 * f64::EPSILON * actual.abs().max(value.abs())
                    {
                        return Err(FdmError::Domain(
                            "initial field conflicts with static Dirichlet boundary data".into(),
                        ));
                    }
                }
            }
        }
    }
    enforce_dirichlet_3d(&mut field, grid, bc);
    // Contradictory edge/corner Dirichlet data is a malformed problem.
    for z in 0..grid.nz {
        for y in 0..grid.ny {
            for x in 0..grid.nx {
                let faces = [
                    (x == 0, &bc.x_min, z * grid.ny + y),
                    (x + 1 == grid.nx, &bc.x_max, z * grid.ny + y),
                    (y == 0, &bc.y_min, z * grid.nx + x),
                    (y + 1 == grid.ny, &bc.y_max, z * grid.nx + x),
                    (z == 0, &bc.z_min, y * grid.nx + x),
                    (z + 1 == grid.nz, &bc.z_max, y * grid.nx + x),
                ];
                let mut prescribed = None;
                for (active, face, index) in faces {
                    if active {
                        let value = match face {
                            | BoundaryCondition::Dirichlet(_)
                            | BoundaryCondition::DirichletProfile(_) => {
                                Some(face.dirichlet_at(index))
                            },
                            | BoundaryCondition::Robin { alpha, beta, gamma }
                                if beta.abs() < 1e-15 =>
                            {
                                Some(gamma / alpha)
                            },
                            | _ => None,
                        };
                        if let Some(value) = value {
                            if let Some(previous) = prescribed
                                && !crate::kernels::exact_eq(previous, value)
                            {
                                return Err(FdmError::Domain(
                                    "conflicting Dirichlet values meet at a grid edge or corner"
                                        .into(),
                                ));
                            }
                            prescribed = Some(value);
                        }
                    }
                }
            }
        }
    }
    Ok(field)
}

pub(super) fn heat_step(
    grid: Grid3D,
    bc: &Boundary3D,
    alpha: f64,
    method: TimeStepping,
    state: &[f64],
    dt: f64,
) -> Result<Vec<f64>, FdmError> {
    let lap = laplacian_3d(state, grid, bc, None);
    let spacing_squared = grid.h * grid.h;
    let r = alpha * dt / spacing_squared;
    let mut next = match method {
        | TimeStepping::Explicit => {
            state
                .iter()
                .zip(lap)
                .map(|(u, l)| u + alpha * dt * l)
                .collect()
        },
        | TimeStepping::BackwardEuler | TimeStepping::CrankNicolson => {
            let theta = if method == TimeStepping::BackwardEuler {
                1.0
            } else {
                0.5
            };
            let homogeneous = homogeneous_3d(bc);
            let boundary = laplacian_3d(&vec![0.0; state.len()], grid, bc, None);
            let rhs: Vec<_> = (0..state.len())
                .map(|i| state[i] + alpha * dt * ((1.0 - theta) * lap[i] + theta * boundary[i]))
                .collect();
            let apply = |x: &[f64], out: &mut [f64]| {
                let l = laplacian_3d(x, grid, &homogeneous, None);
                for i in 0..x.len() {
                    out[i] = x[i] - theta * alpha * dt * l[i];
                }
            };
            let result = crate::kernels::krylov::gmres(
                apply,
                |x, y| y.copy_from_slice(x),
                &rhs,
                Some(state),
                40.min(state.len()),
                1e-12,
                20 * state.len(),
            );
            if !result.converged {
                return Err(FdmError::NotConverged);
            }
            result.x
        },
        | TimeStepping::Adi => {
            let mut stage: Vec<_> = state
                .iter()
                .zip(&lap)
                .map(|(u, l)| u + alpha * dt * l)
                .collect();
            let dims = [grid.nx, grid.ny, grid.nz];
            for axis in 0..3 {
                if dims[axis] == 1 {
                    continue;
                }
                let old = laplacian_3d(state, grid, bc, Some(axis));
                for i in 0..stage.len() {
                    stage[i] -= 0.5 * alpha * dt * old[i];
                }
                let (left, right) = match axis {
                    | 0 => (&bc.x_min, &bc.x_max),
                    | 1 => (&bc.y_min, &bc.y_max),
                    | _ => (&bc.z_min, &bc.z_max),
                };
                for a in 0..dims[(axis + 1) % 3] {
                    for b in 0..dims[(axis + 2) % 3] {
                        let mut point = [0; 3];
                        point[(axis + 1) % 3] = a;
                        point[(axis + 2) % 3] = b;
                        let line: Vec<_> = (0..dims[axis])
                            .map(|j| {
                                point[axis] = j;
                                stage[grid.index(point[0], point[1], point[2])]
                            })
                            .collect();
                        let face = match axis {
                            | 0 => point[2] * grid.ny + point[1],
                            | 1 => point[2] * grid.nx + point[0],
                            | _ => point[1] * grid.nx + point[0],
                        };
                        let solved =
                            solve_line_implicit(&line, 0.5 * r, grid.h, left, right, face)?;
                        for (j, value) in solved.into_iter().enumerate() {
                            point[axis] = j;
                            stage[grid.index(point[0], point[1], point[2])] = value;
                        }
                    }
                }
                enforce_dirichlet_3d(&mut stage, grid, bc);
            }
            stage
        },
    };
    enforce_dirichlet_3d(&mut next, grid, bc);
    Ok(next)
}

/// Adaptive heat evolution in one, two, or three spatial dimensions.
///
/// BE and CN solve their full implicit operators; ADI uses Douglas splitting.
/// Explicit trials are additionally bounded by the diffusion CFL condition.
///
/// # Errors
/// Malformed grids/BCs, invalid controls, inconsistent state, or failed convergence.
#[allow(clippy::too_many_arguments)]
pub fn heat_adaptive(
    extents: &[usize],
    spacing: f64,
    alpha: f64,
    initial: &[f64],
    bc: &Boundary3D,
    method: TimeStepping,
    options: &PdeAdaptiveOptions,
    event: Option<&mut PdeTerminalEvent<'_>>,
) -> Result<PdeAdaptiveSolution, FdmError> {
    let (grid, bc) = geometry(extents, spacing, bc)?;
    let initial = initial_field(initial, grid, &bc)?;
    if !alpha.is_finite() || alpha < 0.0 {
        return Err(FdmError::Domain(
            "heat diffusivity must be finite and nonnegative".into(),
        ));
    }
    let cap = if method == TimeStepping::Explicit && alpha > 0.0 {
        0.49 * spacing * spacing / (extents.len() as f64 * alpha)
    } else {
        f64::INFINITY
    };
    let order = if matches!(method, TimeStepping::Explicit | TimeStepping::BackwardEuler) {
        1
    } else {
        2
    };
    adaptive_evolve(
        &initial,
        order,
        &vec![options.atol; initial.len()],
        options,
        |_, u, h| heat_step(grid, &bc, alpha, method, u, h),
        |_, _| Ok(cap),
        event,
    )
}

/// Adaptive wave evolution by variable-step velocity Verlet.
///
/// Each state concatenates displacement and velocity. `velocity_atol` has
/// displacement/time dimensions, independently of `options.atol`.
///
/// # Errors
/// Malformed grid/BCs, invalid state/tolerances, or exhausted error/CFL budgets.
#[allow(clippy::too_many_arguments)]
pub fn wave_adaptive(
    extents: &[usize],
    spacing: f64,
    c: f64,
    initial: &[f64],
    velocity: &[f64],
    velocity_atol: f64,
    bc: &Boundary3D,
    options: &PdeAdaptiveOptions,
    event: Option<&mut PdeTerminalEvent<'_>>,
) -> Result<PdeAdaptiveSolution, FdmError> {
    let (grid, bc) = geometry(extents, spacing, bc)?;
    let initial = initial_field(initial, grid, &bc)?;
    let n = initial.len();
    if velocity.len() != n || !c.is_finite() || !velocity_atol.is_finite() || velocity_atol <= 0.0 {
        return Err(FdmError::Domain(
            "invalid wave speed, velocity, or velocity tolerance".into(),
        ));
    }
    let mut state = initial;
    let mut v = velocity.to_vec();
    for z in 0..grid.nz {
        for y in 0..grid.ny {
            for x in 0..grid.nx {
                if fixed_value_3d(grid, &bc, x, y, z).is_some() {
                    if !crate::kernels::exact_eq(v[grid.index(x, y, z)], 0.0) {
                        return Err(FdmError::Domain("static Dirichlet wave boundary requires zero initial boundary velocity".into()));
                    }
                    v[grid.index(x, y, z)] = 0.0;
                }
            }
        }
    }
    state.extend(v);
    let mut atol = vec![options.atol; n];
    atol.extend(vec![velocity_atol; n]);
    let cap = if c.abs() > 0.0 {
        0.99 * spacing / (c.abs() * (extents.len() as f64).sqrt())
    } else {
        f64::INFINITY
    };
    adaptive_evolve(
        &state,
        2,
        &atol,
        options,
        |_, state, h| {
            let lap = laplacian_3d(&state[..n], grid, &bc, None);
            let mut next = state.to_vec();
            for i in 0..n {
                next[n + i] += 0.5 * h * c * c * lap[i];
                next[i] += h * next[n + i];
            }
            enforce_dirichlet_3d(&mut next[..n], grid, &bc);
            let lap = laplacian_3d(&next[..n], grid, &bc, None);
            for i in 0..n {
                next[n + i] += 0.5 * h * c * c * lap[i];
            }
            Ok(next)
        },
        |_, _| Ok(cap),
        event,
    )
}

/// Adaptive first-order upwind advection with a multidimensional CFL bound.
///
/// # Errors
/// Invalid grid, face data, velocities, controls, or exhausted temporal budgets.
#[allow(clippy::too_many_arguments)]
pub fn advection_adaptive(
    extents: &[usize],
    spacing: f64,
    velocity: &[f64],
    initial: &[f64],
    bc: &Boundary3D,
    options: &PdeAdaptiveOptions,
    event: Option<&mut PdeTerminalEvent<'_>>,
) -> Result<PdeAdaptiveSolution, FdmError> {
    let (grid, bc) = geometry(extents, spacing, bc)?;
    let initial = initial_field(initial, grid, &bc)?;
    if velocity.len() != extents.len() || velocity.iter().any(|v| !v.is_finite()) {
        return Err(FdmError::Domain(
            "one finite advection speed is required per axis".into(),
        ));
    }
    let speed = velocity.iter().map(|v| v.abs()).sum::<f64>();
    let cap = if speed > 0.0 {
        0.99 * spacing / speed
    } else {
        f64::INFINITY
    };
    adaptive_evolve(
        &initial,
        1,
        &vec![options.atol; initial.len()],
        options,
        |_, state, h| {
            let mut next = state.to_vec();
            for z in 0..grid.nz {
                for y in 0..grid.ny {
                    for x in 0..grid.nx {
                        if fixed_value_3d(grid, &bc, x, y, z).is_some() {
                            continue;
                        }
                        let i = grid.index(x, y, z);
                        for (axis, v) in velocity.iter().enumerate() {
                            let courant = v * h / spacing;
                            next[i] -= courant
                                * upwind_delta_3d(state, grid, &bc, [x, y, z], axis, courant);
                        }
                    }
                }
            }
            enforce_dirichlet_3d(&mut next, grid, &bc);
            Ok(next)
        },
        |_, _| Ok(cap),
        event,
    )
}
