
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/black_hole_orbits_rust.png)

use std::fs::File;
use std::io::Write;

use rayon::prelude::*;
use serde::Deserialize;
use serde::Serialize;

use crate::physics::physics_rkm::DormandPrince54;
use crate::physics::physics_rkm::OdeSystem;

/// Parameters for the geodesic simulation.
#[derive(
    Clone, Debug, Serialize, Deserialize,
)]

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

        -m / r + l * l / (2.0 * r * r)
            - m * l * l / r.powi(3)
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

        let (r, r_dot, _phi, phi_dot) = (
            y[0], y[1], y[2], y[3],
        );

        let l = r * r * phi_dot;

        let r_ddot = -self.mass
            / r.powi(2)
            + l.powi(2) / r.powi(3)
            - 3.0
                * self.mass
                * l.powi(2)
                / r.powi(4);

        let phi_ddot =
            -2.0 * r_dot * phi_dot / r;

        dy[0] = r_dot;

        dy[1] = r_ddot;

        dy[2] = phi_dot;

        dy[3] = phi_ddot;
    }
}

/// Runs a geodesic simulation around a Schwarzschild black hole.
///
/// This function uses an adaptive Runge-Kutta solver (Dormand-Prince 5(4)) to integrate
/// the Schwarzschild geodesic equations. The output is a series of `(x, y)` coordinates
/// representing the path of a particle in the black hole's spacetime.
///
/// # Arguments
/// * `params` - A reference to `GeodesicParameters` containing the black hole mass,
///   initial state of the particle, total proper time, and initial time step.
///
/// # Returns
/// A `Vec` of `(f64, f64)` tuples, where each tuple is an `(x, y)` coordinate
/// in Cartesian space, representing the simulated orbit.

#[must_use]

pub fn run_geodesic_simulation(
    params: &GeodesicParameters
) -> Vec<(f64, f64)> {

    let system = SchwarzschildSystem {
        mass: params.black_hole_mass,
    };

    let solver =
        DormandPrince54::default();

    let t_span = (
        0.0,
        params.proper_time_end,
    );

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

            (
                r * phi.cos(),
                r * phi.sin(),
            )
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
/// This function will return an error if it fails to create or write to the output CSV files.

pub fn simulate_black_hole_orbits_scenario()
-> std::io::Result<()> {

    println!(
        "Running Black Hole orbit \
         simulation..."
    );

    let black_hole_mass = 1.0;

    let stable_orbit_params =
        GeodesicParameters {
            black_hole_mass,
            initial_state: [
                10.0, 0.0, 0.0, 0.035,
            ],
            proper_time_end: 1500.0,
            initial_dt: 0.1,
        };

    let plunging_orbit_params =
        GeodesicParameters {
            black_hole_mass,
            initial_state: [
                10.0, 0.0, 0.0, 0.02,
            ],
            proper_time_end: 500.0,
            initial_dt: 0.1,
        };

    let photon_orbit_params =
        GeodesicParameters {
            black_hole_mass,
            initial_state: [
                10.0, -1.0, 0.0, 0.03,
            ],
            proper_time_end: 50.0,
            initial_dt: 0.01,
        };

    let orbits = vec![
        (
            "stable_orbit",
            stable_orbit_params,
        ),
        (
            "plunging_orbit",
            plunging_orbit_params,
        ),
        (
            "photon_orbit",
            photon_orbit_params,
        ),
    ];

    orbits
        .into_par_iter()
        .for_each(|(name, params)| {

            println!(
                "Simulating {name}..."
            );

            let path =
                run_geodesic_simulation(
                    &params,
                );

            let filename = format!(
                "orbit_{name}.csv"
            );

            if let Ok(mut file) =
                File::create(&filename)
            {

                let _ = writeln!(
                    file,
                    "x,y"
                );

                for (x, y) in path {

                    let _ = writeln!(
                        file,
                        "{x},{y}"
                    );
                }

                println!(
                    "Saved path to \
                     {filename}"
                );
            }
        });

    Ok(())
}
