//! Advanced spatial discretizations share the verified temporal controller.
use super::FdmError;
use super::FvmScheme;
use super::solve_fvm_conservation_1d;
use super::solve_spectral_etdrk4;
use crate::sim::fdm_nd::adaptive::PdeAdaptiveOptions;
use crate::sim::fdm_nd::adaptive::PdeAdaptiveSolution;
use crate::sim::fdm_nd::adaptive::PdeTerminalEvent;
use crate::sim::fdm_nd::adaptive::adaptive_evolve;
use crate::sim::fem_simplex::FemProblem;
use crate::sim::fem_simplex::FemSpace;
use crate::sim::structured::LinearSolver;

/// Adaptive FEM heat evolution from existing nodal degrees of freedom.
///
/// Theta one is backward Euler; theta one half is Crank-Nicolson. This controls
/// temporal error on a fixed mesh; spatial mesh refinement is a separate solve.
///
/// # Errors
/// Invalid nodal data, theta, controls, boundary data or linear-system failure.
pub fn fem_heat_adaptive(
    space: &FemSpace,
    problem: &FemProblem<'_>,
    initial: &[f64],
    theta: f64,
    options: &PdeAdaptiveOptions,
    event: Option<&mut PdeTerminalEvent<'_>>,
) -> Result<PdeAdaptiveSolution, FdmError> {
    if initial.len() != space.dofs.len() || !(0.5..=1.0).contains(&theta) {
        return Err(FdmError::Domain(
            "invalid FEM adaptive state or theta".into(),
        ));
    }
    let order = if crate::kernels::exact_eq(theta, 0.5) {
        2
    } else {
        1
    };
    adaptive_evolve(
        initial,
        order,
        &vec![options.atol; initial.len()],
        options,
        |t, u, h| {
            crate::sim::fem_simplex::heat_step(space, problem, u, t, h, theta, LinearSolver::Auto)
                .map_err(FdmError::Domain)
        },
        |_, _| Ok(f64::INFINITY),
        event,
    )
}

/// Adaptive scalar finite-volume conservation law on fixed uniform cells.
///
/// The selected flux/reconstruction is retained on every full and half step.
/// Classical Lax-Friedrichs has timestep-dependent numerical viscosity; strict
/// temporal tolerances may therefore fail at a fixed spatial grid, explicitly.
///
/// # Errors
/// Invalid domain/state, non-finite flux/speed, CFL underflow or error-budget exhaustion.
#[allow(clippy::too_many_arguments)]
pub fn fvm_adaptive<F, G>(
    domain_size: f64,
    initial: &[f64],
    periodic: bool,
    scheme: FvmScheme,
    flux: F,
    wave_speed: G,
    options: &PdeAdaptiveOptions,
    event: Option<&mut PdeTerminalEvent<'_>>,
) -> Result<PdeAdaptiveSolution, FdmError>
where
    F: Fn(f64) -> f64,
    G: Fn(f64) -> f64,
{
    if initial.len() < 3 || !domain_size.is_finite() || domain_size <= 0.0 {
        return Err(FdmError::Domain("invalid adaptive FVM domain".into()));
    }
    let dx = domain_size / initial.len() as f64;
    let order = if matches!(scheme, FvmScheme::MusclMinmod | FvmScheme::MusclVanLeer) {
        2
    } else {
        1
    };
    adaptive_evolve(
        initial,
        order,
        &vec![options.atol; initial.len()],
        options,
        |_, state, h| {
            let (_, next) = solve_fvm_conservation_1d(
                state.len(),
                domain_size,
                h,
                1,
                periodic,
                scheme,
                |x| state[((x / dx).floor() as usize).min(state.len() - 1)],
                &flux,
                &wave_speed,
            )?;
            Ok(next)
        },
        |_, state| {
            let mut speed = 0.0_f64;
            for u in state {
                let v = wave_speed(*u).abs();
                if !v.is_finite() {
                    return Err(FdmError::Domain("FVM wave speed must be finite".into()));
                }
                speed = speed.max(v);
            }
            Ok(if speed > 0.0 {
                0.9 * dx / speed
            } else {
                f64::INFINITY
            })
        },
        event,
    )
}

/// Adaptive periodic scalar DG advection with the existing SSPRK3 time method.
///
/// # Errors
/// Invalid mesh/nodal data, failed CFL checks or exhausted temporal budgets.
pub fn dg_advection_adaptive(
    mesh: &crate::sim::dg::DgMesh,
    velocity: f64,
    initial: &[f64],
    options: &PdeAdaptiveOptions,
    event: Option<&mut PdeTerminalEvent<'_>>,
) -> Result<PdeAdaptiveSolution, FdmError> {
    dg_advection_nd_adaptive(mesh, [velocity, 0.0], initial, options, event)
}

/// Adaptive Fourier ETDRK4 reaction-diffusion evolution in one to three dimensions.
///
/// Both linear diffusion/reaction and nonlinear callback physical times are
/// retained on each restarted ETDRK4 step. The temporal order remains four.
///
/// # Errors
/// Invalid periodic grid, shape, coefficients, callback shape/domain or temporal controls.
#[allow(clippy::too_many_arguments)]
pub fn spectral_etdrk4_adaptive<N>(
    extents: &[usize],
    lengths: &[f64],
    alpha: f64,
    sigma: f64,
    initial: &[f64],
    nonlinear: N,
    options: &PdeAdaptiveOptions,
    event: Option<&mut PdeTerminalEvent<'_>>,
) -> Result<PdeAdaptiveSolution, FdmError>
where
    N: Fn(f64, &[f64]) -> Vec<f64>,
{
    if lengths.iter().any(|v| !v.is_finite() || *v <= 0.0)
        || extents.iter().try_fold(1usize, |n, d| n.checked_mul(*d)) != Some(initial.len())
    {
        return Err(FdmError::Domain(
            "invalid spectral lengths or extent product".into(),
        ));
    }
    let grid =
        crate::sim::spectral::SpectralGrid::new(extents, lengths).map_err(FdmError::Domain)?;
    if initial.len() != grid.size() || !alpha.is_finite() || alpha < 0.0 || !sigma.is_finite() {
        return Err(FdmError::Domain("invalid adaptive spectral data".into()));
    }
    adaptive_evolve(
        initial,
        4,
        &vec![options.atol; initial.len()],
        options,
        |t, state, h| {
            let callback = &nonlinear;
            let expected = state.len();
            solve_spectral_etdrk4(
                extents,
                lengths,
                alpha,
                sigma,
                h,
                1,
                state,
                move |local, u| {
                    let value = callback(t + local, u);
                    if value.len() == expected {
                        value
                    } else {
                        vec![f64::NAN; expected]
                    }
                },
            )
        },
        |_, _| Ok(f64::INFINITY),
        event,
    )
}

/// Adaptive periodic one- or two-dimensional nodal DG advection, using SSPRK3.
///
/// # Errors
/// Invalid mesh/velocity/state, CFL bound, temporal controls or exhausted budgets.
pub fn dg_advection_nd_adaptive(
    mesh: &crate::sim::dg::DgMesh,
    velocity: [f64; 2],
    initial: &[f64],
    options: &PdeAdaptiveOptions,
    event: Option<&mut PdeTerminalEvent<'_>>,
) -> Result<PdeAdaptiveSolution, FdmError> {
    if !(1..=2).contains(&mesh.dim)
        || mesh.p == 0
        || mesh.ne[..mesh.dim].contains(&0)
        || velocity.iter().any(|v| !v.is_finite())
        || (0..mesh.dim)
            .any(|a| !mesh.lo[a].is_finite() || !mesh.hi[a].is_finite() || mesh.hi[a] <= mesh.lo[a])
    {
        return Err(FdmError::Domain(
            "invalid adaptive DG mesh or velocity".into(),
        ));
    }
    let rate = (0..mesh.dim)
        .map(|a| velocity[a].abs() * mesh.ne[a] as f64 / (mesh.hi[a] - mesh.lo[a]))
        .sum::<f64>()
        * (2.0 * mesh.p as f64 + 1.0);
    let cap = if rate > 0.0 {
        0.9 / rate
    } else {
        f64::INFINITY
    };
    adaptive_evolve(
        initial,
        3,
        &vec![options.atol; initial.len()],
        options,
        |t, state, h| {
            crate::sim::dg::advection_step_nd(mesh, velocity, state, t, h).map_err(FdmError::Domain)
        },
        |_, _| Ok(cap),
        event,
    )
}

/// Adaptive FEM heat with a serial vectorized nodal source callback.
///
/// `source(t, coordinates)` returns one physical field/time value per node.
/// Its consistent mass load is added to the problem's existing source load.
///
/// # Errors
/// Invalid source/state, theta/controls, callback errors or failed linear solve.
#[allow(clippy::too_many_arguments)]
pub fn fem_heat_source_adaptive<S>(
    space: &FemSpace,
    problem: &FemProblem<'_>,
    initial: &[f64],
    theta: f64,
    source: S,
    options: &PdeAdaptiveOptions,
    event: Option<&mut PdeTerminalEvent<'_>>,
) -> Result<PdeAdaptiveSolution, FdmError>
where
    S: Fn(f64, &[[f64; 3]]) -> Result<Vec<f64>, FdmError>,
{
    if initial.len() != space.dofs.len() || !(0.5..=1.0).contains(&theta) {
        return Err(FdmError::Domain(
            "invalid adaptive FEM source problem".into(),
        ));
    }
    let order = if crate::kernels::exact_eq(theta, 0.5) {
        2
    } else {
        1
    };
    adaptive_evolve(
        initial,
        order,
        &vec![options.atol; initial.len()],
        options,
        |t, state, h| {
            let before = source(t, &space.dofs)?;
            if before.len() != state.len() || before.iter().any(|v| !v.is_finite()) {
                return Err(FdmError::Domain(
                    "FEM source callback returned invalid nodal data".into(),
                ));
            }
            let after = source(t + h, &space.dofs)?;
            crate::sim::fem_simplex::heat_step_nodal_source(
                space,
                problem,
                state,
                t,
                h,
                theta,
                LinearSolver::Auto,
                &before,
                &after,
            )
            .map_err(FdmError::Domain)
        },
        |_, _| Ok(f64::INFINITY),
        event,
    )
}
