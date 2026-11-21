//! # Geodesic Relativity Simulation
//!
//! This module provides tools for simulating geodesic motion in the spacetime of a Schwarzschild black hole.
//! It calculates the paths of particles (orbits) and light rays by numerically integrating the geodesic equations
//! derived from the Schwarzschild metric.
//!
//! # Overview
//!
//! The simulation assumes a spherically symmetric, non-rotating massive body (a Schwarzschild black hole).
//! The dynamics are governed by the relativistic geodesic equations, which describe how massive particles and
//! photons move under the influence of the black hole's curvature.
//!
//! Key features include:
//! - **Effective Potential Calculation**: Helper methods to calculate the effective potential for analysis of orbits.
//! - **ODE System Definition**: `SchwarzschildSystem` defines the equations of motion for use with numerical ODE solvers.
//! - **Integration**: Uses implicit midpoint for equatorial timelike trajectories and retains adaptive Dormand–Prince 5(4) for Schwarzschild.
//! - **Scenario Simulation**: Preconfigured timelike orbits, including plunges. The historical CSV name `photon_orbit` refers to a timelike example.
//!
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/black_hole_orbits_rust.png)
//!
//! ## Coverage
//!
//! | Capability | Support |
//! |---|---|
//! | Spacetime | `run_geodesic_simulation`: Schwarzschild; [`run_kerr_newman_geodesic`]: Kerr-Newman (Kerr, Reissner-Nordstroem, Schwarzschild as special cases); equatorial plane |
//! | Integration | both delegate to the symplectic implicit-midpoint engine of [`crate::sim::relativity`] with mass-shell projection (energy and angular momentum are conserved exactly); the legacy Dormand-Prince integrator is the fallback when the initial data lie where the engine cannot start (inside the horizon) |
//! | Events | integration stops slightly outside the horizon (`r <= 1.001 r_+`) or when the particle escapes (`r >= 1e4`) |
//!
//! The initial state `[r, dr/dtau, phi, dphi/dtau]` is converted to the
//! conserved energy and angular momentum of the chosen spacetime by the
//! timelike normalisation `g_{mu nu} u^mu u^nu = -1`. Photon orbits, ISCO and
//! general (non-equatorial) geodesics are in [`crate::sim::relativity`].

use std::fs::File;
use std::io::Write;
use std::num::NonZeroUsize;
use std::path::Path;

use rayon::prelude::*;
use serde::Deserialize;
use serde::Serialize;

use crate::sim::physics_rkm::DormandPrince54;
use crate::sim::physics_rkm::OdeSystem;
use crate::sim::relativity::Geodesic;
use crate::sim::relativity::KerrNewman;
use crate::sim::relativity::Termination;

// Scoped native accounting bounds both implicit steps and legacy adaptive RHS calls.
// The Arc is explicitly inherited by parallel file-scenario workers.
std::thread_local! {static NATIVE_GEODESIC_BUDGET:std::cell::RefCell<Option<std::sync::Arc<std::sync::atomic::AtomicUsize>>>=const{std::cell::RefCell::new(None)};}
struct GeodesicBudgetExhausted;
pub(crate) fn shared_native_budget<T>(
    budget: Option<std::sync::Arc<std::sync::atomic::AtomicUsize>>,
    run: impl FnOnce() -> T,
) -> Result<T, String> {
    struct Restore(Option<std::sync::Arc<std::sync::atomic::AtomicUsize>>);
    impl Drop for Restore {
        fn drop(&mut self) {
            NATIVE_GEODESIC_BUDGET.with(|b| {
                b.replace(self.0.take());
            });
        }
    }
    let prior = NATIVE_GEODESIC_BUDGET.with(|b| b.replace(budget));
    let _restore = Restore(prior);
    match std::panic::catch_unwind(std::panic::AssertUnwindSafe(run)) {
        | Ok(v) => Ok(v),
        | Err(e) if e.is::<GeodesicBudgetExhausted>() => {
            Err("native geodesic step/RHS budget exhausted".into())
        },
        | Err(e) => std::panic::resume_unwind(e),
    }
}
pub(crate) fn with_native_budget<T>(
    budget: usize,
    run: impl FnOnce() -> T,
) -> Result<T, String> {
    shared_native_budget(
        Some(std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(
            budget,
        ))),
        run,
    )
}
pub(crate) fn inherited_native_budget() -> Option<std::sync::Arc<std::sync::atomic::AtomicUsize>> {
    NATIVE_GEODESIC_BUDGET.with(|b| b.borrow().clone())
}
pub(crate) fn native_charge() {
    NATIVE_GEODESIC_BUDGET.with(|b| {
        if let Some(budget) = b.borrow().as_ref() {
            let mut old = budget.load(std::sync::atomic::Ordering::Relaxed);
            loop {
                let Some(next) = old.checked_sub(1) else {
                    std::panic::resume_unwind(Box::new(GeodesicBudgetExhausted));
                };
                match budget.compare_exchange_weak(
                    old,
                    next,
                    std::sync::atomic::Ordering::Relaxed,
                    std::sync::atomic::Ordering::Relaxed,
                ) {
                    | Ok(_) => break,
                    | Err(value) => old = value,
                }
            }
        }
    });
}

/// Parameters for the geodesic simulation.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct GeodesicParameters {
    /// The mass of the black hole.
    pub black_hole_mass: f64,
    /// Initial state: `[r, dr/dτ, φ, dφ/dτ]`
    pub initial_state: [f64; 4],
    /// Total proper time for the simulation.
    pub proper_time_end: f64,
    /// Initial time step for the adaptive solver.
    pub initial_dt: f64,
}

impl GeodesicParameters {
    /// Calculates the effective potential for a Schwarzschild black hole.
    /// `V_eff(r)` = -M/r + L^2/(2r^2) - ML^2/r^3
    #[must_use]
    pub fn effective_potential(
        &self,
        r: f64,
        l: f64,
    ) -> f64 {
        let m = self.black_hole_mass;

        -m / r + l * l / (2.0 * r * r) - m * l * l / r.powi(3)
    }
}

/// Represents the Schwarzschild geodesic equations as a system of first-order ODEs.
pub struct SchwarzschildSystem {
    mass: f64,
}

impl OdeSystem for SchwarzschildSystem {
    fn dim(&self) -> usize {
        4
    }

    fn eval(
        &self,
        _t: f64,
        y: &[f64],
        dy: &mut [f64],
    ) {
        native_charge();
        let (r, r_dot, _phi, phi_dot) = (y[0], y[1], y[2], y[3]);

        let l = r * r * phi_dot;

        let r_ddot = -self.mass / r.powi(2) + l.powi(2) / r.powi(3)
            - 3.0 * self.mass * l.powi(2) / r.powi(4);

        let phi_ddot = -2.0 * r_dot * phi_dot / r;

        dy[0] = r_dot;

        dy[1] = r_ddot;

        dy[2] = phi_dot;

        dy[3] = phi_ddot;
    }
}

/// The sign of the conserved Killing energy of a timelike equatorial state.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EquatorialEnergyBranch {
    /// Nonnegative energy. Continues the unique future-directed exterior root.
    Nonnegative,
    /// The second future-directed root in the ergoregion, with nonpositive energy.
    Nonpositive,
}

/// Prepares a timelike equatorial Kerr–Newman geodesic in geometric units.
///
/// The proper-time state is `[r, dr/dτ, φ, dφ/dτ]`. Normalization selects
/// a root with `dt/dτ > 0` in Boyer–Lindquist coordinates and converts it
/// to the conserved energy per unit mass, angular momentum and radial
/// momentum used by the Hamiltonian integrator. The radius must be
/// outside the outer horizon when a horizon exists.
///
/// The nonnegative branch includes the stationary-limit surface. Within the
/// ergoregion, a compatible spatial velocity can give two future-directed
/// roots, one for each sign of the Killing energy. The branches coincide
/// when the discriminant vanishes. The nonpositive branch is available only
/// where the stationary Killing vector is spacelike.
///
/// # Errors
///
/// Returns an error for nonfinite inputs, a radius on or inside the outer
/// horizon, or a state without a representable future-directed root on
/// the selected branch.
pub fn equatorial_timelike_geodesic(
    metric: KerrNewman,
    state: [f64; 4],
    branch: EquatorialEnergyBranch,
) -> Result<Geodesic, String> {
    let [r, rd, phi, phid] = state;
    let (m, a, q) = (metric.m, metric.a, metric.q);
    if ![m, a, q, r, rd, phi, phid].into_iter().all(f64::is_finite) {
        return Err("geodesic metric and initial state must be finite".into());
    }
    let delta = r * r - 2.0 * m * r + a * a + q * q;
    if r <= 0.0 || !delta.is_finite() || delta <= 0.0 || metric.horizon().is_some_and(|rh| r <= rh)
    {
        return Err("initial radius is not outside the horizon".into());
    }

    let w = (2.0 * m * r - q * q) / (r * r);
    let gtt = -(1.0 - w);
    let gtp = -a * w;
    let gpp = r * r + a * a + a * a * w;
    let grr = r * r / delta;
    let b = gtp * phid;
    let c = gpp * phid * phid + grr * rd * rd + 1.0;
    if !(gpp > 0.0 && c.is_finite() && gtt.is_finite() && b.is_finite()) {
        return Err(
            "initial metric or velocity coefficients are outside the numerical range".into(),
        );
    }
    if (gtt >= 0.0 && b >= 0.0)
        || (matches!(branch, EquatorialEnergyBranch::Nonpositive) && gtt <= 0.0)
    {
        return Err("no future-directed timelike velocity on this energy branch".into());
    }

    let disc = if gtt >= 0.0 {
        // The equatorial block satisfies g_tphi² - g_tt * g_phiphi = Δ.
        // Using that identity avoids subtracting two large rotation terms.
        let azimuth = delta.sqrt() * phid;
        azimuth.mul_add(azimuth, -gtt * (grr * rd * rd + 1.0))
    } else {
        b.powi(2) - gtt * c
    };
    if !disc.is_finite() || disc < 0.0 {
        return Err("no future-directed timelike velocity with this spatial velocity".into());
    }
    let root = disc.sqrt().abs();
    let (ut, energy) = match branch {
        | EquatorialEnergyBranch::Nonpositive => ((root - b) / gtt, -root),
        | EquatorialEnergyBranch::Nonnegative if gtt >= 0.0 || b < -0.5 * root => {
            // Rationalizing the positive-energy root gives the finite
            // linear limit at g_tt = 0 and avoids cancellation there.
            (c / (root - b), root)
        },
        | EquatorialEnergyBranch::Nonnegative => {
            let ut = (-b - root) / gtt;
            (ut, -(gtt * ut + b))
        },
    };
    let lz = gtp * ut + gpp * phid;
    let pr = grr * rd;
    if !(ut > 0.0 && ut.is_finite() && energy.is_finite() && lz.is_finite() && pr.is_finite()) {
        return Err("future-directed timelike velocity is outside the numerical range".into());
    }
    Ok(Geodesic {
        metric,
        t: 0.0,
        r,
        theta: std::f64::consts::FRAC_PI_2,
        phi,
        pr,
        ptheta: 0.0,
        energy,
        lz,
        lambda: 0.0,
        mu2: 1.0,
    })
}

/// Equatorial timelike Kerr-Newman geodesic on the symplectic engine.
///
/// The spacetime has (mass `params.black_hole_mass`, spin
/// parameter `spin = a`, charge `charge = Q`, geometric units) integrated for the proper time
/// `params.proper_time_end` with the symplectic implicit-midpoint engine of
/// [`crate::sim::relativity`] (step `params.initial_dt`, mass-shell projection). This uses
/// the [`EquatorialEnergyBranch::Nonnegative`] branch. The initial state is
/// `[r, dr/dtau, phi, dphi/dtau]` as in [`GeodesicParameters`]. Returns the Cartesian path
/// `(r cos phi, r sin phi)`; the run stops at the horizon (`r <= 1.001 r_+`) or on escape.
///
/// # Errors
/// Fails when the initial data are not a valid timelike state outside the horizon or the
/// implicit step does not converge.
pub fn run_kerr_newman_geodesic(
    params: &GeodesicParameters,
    spin: f64,
    charge: f64,
) -> Result<Vec<(f64, f64)>, String> {
    run_kerr_newman_geodesic_sampled(
        params,
        spin,
        charge,
        EquatorialEnergyBranch::Nonnegative,
        NonZeroUsize::MIN,
    )
}

/// Runs an equatorial timelike Kerr–Newman geodesic and retains fewer positions.
///
/// Returns the initial Cartesian position, every `retain_every` integration
/// steps, and the final reached position, including an event between samples.
/// A stride of one gives the full trajectory. The stride changes output
/// storage only. The same initial step, proper-time target, horizon and
/// escape tests, and model step limit apply at every stride. The selected
/// energy branch has the meaning of [`equatorial_timelike_geodesic`].
///
/// # Errors
///
/// Returns an error for an invalid initial state, a nonfinite or negative
/// proper-time target, a nonfinite or nonpositive initial step, or a failed
/// implicit integration step.
pub fn run_kerr_newman_geodesic_sampled(
    params: &GeodesicParameters,
    spin: f64,
    charge: f64,
    branch: EquatorialEnergyBranch,
    retain_every: NonZeroUsize,
) -> Result<Vec<(f64, f64)>, String> {
    if !valid_time_span(params) {
        return Err(
            "proper time must be finite and nonnegative and step must be finite and positive"
                .into(),
        );
    }
    let metric = KerrNewman {
        m: params.black_hole_mass,
        a: spin,
        q: charge,
    };
    let mut g = equatorial_timelike_geodesic(metric, params.initial_state, branch)?;
    let r_stop = metric.horizon().map_or(0.05 * metric.m, |rh| 1.001 * rh);
    if g.r <= r_stop {
        return Err("initial radius is inside the stopping radius".into());
    }
    let mut path = vec![(g.r * g.phi.cos(), g.r * g.phi.sin())];
    let h0 = params.initial_dt.max(1e-6);
    let mut tau = 0.0;
    let mut guard = 0usize;
    let stride = retain_every.get();
    let mut until_sample = stride;
    while tau < params.proper_time_end && guard < 2_000_000 {
        native_charge();
        guard += 1;
        let h = h0.min(params.proper_time_end - tau);
        let before = g.lambda;
        let term = g.integrate(h, 1, true, r_stop, 1e4)?;
        let elapsed = g.lambda - before;
        if !(elapsed > 0.0 && elapsed.is_finite() && g.r.is_finite() && g.phi.is_finite()) {
            return Err("geodesic step did not produce a finite, advancing state".into());
        }
        tau += elapsed;
        until_sample -= 1;
        if until_sample == 0 {
            path.push((g.r * g.phi.cos(), g.r * g.phi.sin()));
            until_sample = stride;
        }
        if term != Termination::Steps {
            break;
        }
    }
    if until_sample != stride {
        path.push((g.r * g.phi.cos(), g.r * g.phi.sin()));
    }
    Ok(path)
}

fn valid_time_span(params: &GeodesicParameters) -> bool {
    params.proper_time_end.is_finite()
        && params.proper_time_end >= 0.0
        && params.initial_dt.is_finite()
        && params.initial_dt > 0.0
}

/// Runs a geodesic simulation around a Schwarzschild black hole.
///
/// Delegates to [`run_kerr_newman_geodesic`] (`a = Q = 0`): the symplectic implicit-midpoint engine
/// of [`crate::sim::relativity`]. If the engine cannot start from the given state it falls back to
/// the adaptive Dormand-Prince 5(4) integration of the reduced Schwarzschild equations.
///
/// # Arguments
/// * `params` - A reference to `GeodesicParameters` containing the black hole mass,
///   initial state of the particle, total proper time, and initial time step.
///
/// # Returns
/// A `Vec` of `(f64, f64)` tuples, where each tuple is an `(x, y)` coordinate
/// in Cartesian space, representing the simulated orbit.
#[must_use]
pub fn run_geodesic_simulation(params: &GeodesicParameters) -> Vec<(f64, f64)> {
    if let Ok(path) = run_kerr_newman_geodesic(params, 0.0, 0.0) {
        return path;
    }
    run_geodesic_simulation_rk45(params)
}

/// Adaptive Dormand-Prince 5(4) integration of the reduced Schwarzschild equations (the legacy
/// path, kept as an independent reference and as the fallback of [`run_geodesic_simulation`]).
///
/// An invalid proper-time target or initial step returns only the initial Cartesian position.
#[must_use]
pub fn run_geodesic_simulation_rk45(params: &GeodesicParameters) -> Vec<(f64, f64)> {
    if !valid_time_span(params) {
        let [r, _, phi, _] = params.initial_state;
        return vec![(r * phi.cos(), r * phi.sin())];
    }
    let system = SchwarzschildSystem {
        mass: params.black_hole_mass,
    };

    let solver = DormandPrince54::default();

    let t_span = (0.0, params.proper_time_end);

    let tolerance = (1e-7, 1e-7);

    let history = solver.solve(
        &system,
        &params.initial_state,
        t_span,
        params.initial_dt,
        tolerance,
    );

    history
        .iter()
        .map(|(_t, state)| {
            let r = state[0];

            let phi = state[2];

            (r * phi.cos(), r * phi.sin())
        })
        .collect()
}

/// An example scenario that simulates several types of orbits around a black hole.
///
/// This function sets up and runs simulations for:
/// - A stable, precessing orbit.
/// - A plunging orbit (where the particle falls into the black hole).
/// - A fast timelike approach, saved under the historical name `photon_orbit`.
///
/// The results are saved to `.csv` files for external visualization.
///
/// # Errors
///
/// This function will return an error if it fails to create the output directory.
/// The files `orbit_<name>.csv` are written into `output_dir`, which is created
/// if missing; a file that cannot be created is skipped.
pub fn simulate_black_hole_orbits_scenario(output_dir: &Path) -> std::io::Result<()> {
    std::fs::create_dir_all(output_dir)?;

    println!(
        "Running Black Hole orbit \
         simulation..."
    );

    let black_hole_mass = 1.0;

    let stable_orbit_params = GeodesicParameters {
        black_hole_mass,
        initial_state: [10.0, 0.0, 0.0, 0.035],
        proper_time_end: 1500.0,
        initial_dt: 0.1,
    };

    let plunging_orbit_params = GeodesicParameters {
        black_hole_mass,
        initial_state: [10.0, 0.0, 0.0, 0.02],
        proper_time_end: 500.0,
        initial_dt: 0.1,
    };

    let photon_orbit_params = GeodesicParameters {
        black_hole_mass,
        initial_state: [10.0, -1.0, 0.0, 0.03],
        proper_time_end: 50.0,
        initial_dt: 0.01,
    };

    let orbits = vec![
        ("stable_orbit", stable_orbit_params),
        ("plunging_orbit", plunging_orbit_params),
        ("photon_orbit", photon_orbit_params),
    ];

    let budget = NATIVE_GEODESIC_BUDGET.with(|b| b.borrow().clone());
    let results: Vec<Result<(), String>> = orbits
        .into_par_iter()
        .map(|(name, params)| {
            shared_native_budget(budget.clone(), || {
                println!("Simulating {name}...");

                let path = run_geodesic_simulation(&params);

                let filename = output_dir.join(format!("orbit_{name}.csv"));

                if let Ok(mut file) = File::create(&filename) {
                    let _ = writeln!(file, "x,y");

                    for (x, y) in path {
                        let _ = writeln!(file, "{x},{y}");
                    }

                    println!(
                        "Saved path to \
                     {}",
                        filename.display()
                    );
                }
            })
        })
        .collect();
    for result in results {
        result.map_err(std::io::Error::other)?;
    }

    Ok(())
}
