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
//! - **Adaptive Integration**: Uses a high-order adaptive Runge-Kutta method (Dormand-Prince 5(4)) for precision.
//! - **Scenario Simulation**: Pre-configured scenarios for stable orbits, plunging trajectories, and photon deflection.
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
use std::path::Path;

use rayon::prelude::*;
use serde::Deserialize;
use serde::Serialize;

use crate::sim::physics_rkm::DormandPrince54;
use crate::sim::physics_rkm::OdeSystem;
use crate::sim::relativity::{Geodesic, KerrNewman, Termination};

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

/// Converts the proper-time state `[r, dr/dtau, phi, dphi/dtau]` of a timelike equatorial geodesic of
/// the Kerr-Newman spacetime `metric` into the conserved energy / angular momentum and the radial
/// momentum of the Hamiltonian engine.
fn equatorial_geodesic(
    metric: KerrNewman,
    state: [f64; 4],
) -> Result<Geodesic, String> {
    let [r, rd, phi, phid] = state;
    let (m, a, q) = (metric.m, metric.a, metric.q);
    let delta = r * r - 2.0 * m * r + a * a + q * q;
    if delta <= 0.0 || r <= 0.0 {
        return Err("initial radius is not outside the horizon".into());
    }
    let w = (2.0 * m * r - q * q) / (r * r);
    let gtt = -(1.0 - w);
    let gtp = -a * w;
    let gpp = r * r + a * a + a * a * w;
    let grr = r * r / delta;
    let b = gtp * phid;
    let c = gpp * phid * phid + grr * rd * rd + 1.0;
    let disc = b.powi(2) - gtt * c;
    if gtt >= 0.0 || disc < 0.0 {
        return Err("no future-directed timelike velocity with this spatial velocity".into());
    }
    let ut = (-b - disc.sqrt()) / gtt;
    let energy = -(gtt * ut + gtp * phid);
    let lz = gtp * ut + gpp * phid;
    Ok(Geodesic {
        metric,
        t: 0.0,
        r,
        theta: std::f64::consts::FRAC_PI_2,
        phi,
        pr: grr * rd,
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
/// [`crate::sim::relativity`] (step `params.initial_dt`, mass-shell projection). The state is
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
    let metric = KerrNewman { m: params.black_hole_mass, a: spin, q: charge };
    let mut g = equatorial_geodesic(metric, params.initial_state)?;
    let r_stop = metric.horizon().map_or(0.05 * metric.m, |rh| 1.001 * rh);
    if g.r <= r_stop {
        return Err("initial radius is inside the stopping radius".into());
    }
    let mut path = vec![(g.r * g.phi.cos(), g.r * g.phi.sin())];
    let h0 = params.initial_dt.max(1e-6);
    let mut tau = 0.0;
    let mut guard = 0usize;
    while tau < params.proper_time_end && guard < 2_000_000 {
        guard += 1;
        let h = h0.min(params.proper_time_end - tau);
        let before = g.lambda;
        let term = g.integrate(h, 1, true, r_stop, 1e4)?;
        tau += g.lambda - before;
        path.push((g.r * g.phi.cos(), g.r * g.phi.sin()));
        if term != Termination::Steps {
            break;
        }
    }
    Ok(path)
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
#[must_use]
pub fn run_geodesic_simulation_rk45(params: &GeodesicParameters) -> Vec<(f64, f64)> {
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
/// - A photon orbit (light bending).
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

    orbits.into_par_iter().for_each(|(name, params)| {
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
    });

    Ok(())
}
