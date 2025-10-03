// src/physics/physics_sim/geodesic_relativity.rs
// Simulator for calculating geodesic paths around a Schwarzschild black hole.

use crate::physics::physics_rkm::{DormandPrince54, OdeSystem};
use std::fs::File;
use std::io::Write;

/// Parameters for the geodesic simulation.
pub struct GeodesicParameters {
    pub black_hole_mass: f64,
    /// Initial state: `[r, dr/dτ, φ, dφ/dτ]`
    pub initial_state: [f64; 4],
    /// Total proper time for the simulation.
    pub proper_time_end: f64,
    /// Initial time step for the adaptive solver.
    pub initial_dt: f64,
}

/// Represents the Schwarzschild geodesic equations as a system of first-order ODEs.
pub struct SchwarzschildSystem {
    mass: f64,
}

impl OdeSystem for SchwarzschildSystem {
    fn dim(&self) -> usize {
        4
    }

    fn eval(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
        let (r, r_dot, _phi, phi_dot) = (y[0], y[1], y[2], y[3]);

        // From the geodesic equations in the orbital plane:
        // d²r/dτ² = -M/r² + L²/r³ - 3ML²/r⁴
        // d²φ/dτ² = -2(dr/dτ)(dφ/dτ)/r
        // Where L is the specific angular momentum, L = r²(dφ/dτ). It's a conserved quantity.
        let l = r * r * phi_dot;

        let r_ddot = -self.mass / r.powi(2) + l.powi(2) / r.powi(3)
            - 3.0 * self.mass * l.powi(2) / r.powi(4);
        let phi_ddot = -2.0 * r_dot * phi_dot / r;

        // Convert to a system of first-order ODEs:
        // y[0] = r, y[1] = dr/dτ, y[2] = φ, y[3] = dφ/dτ
        dy[0] = r_dot; // dr/dτ
        dy[1] = r_ddot; // d²r/dτ²
        dy[2] = phi_dot; // dφ/dτ
        dy[3] = phi_ddot; // d²φ/dτ²
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
pub fn run_geodesic_simulation(params: &GeodesicParameters) -> Vec<(f64, f64)> {
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

    // Convert polar (r, φ) coordinates from the simulation result to Cartesian (x, y) for plotting.
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
pub fn simulate_black_hole_orbits_scenario() {
    println!("Running Black Hole orbit simulation...");

    let black_hole_mass = 1.0;

    // Case 1: A stable, precessing orbit
    let stable_orbit_params = GeodesicParameters {
        black_hole_mass,
        initial_state: [10.0, 0.0, 0.0, 0.035], // r, dr/dτ, φ, dφ/dτ
        proper_time_end: 1500.0,
        initial_dt: 0.1,
    };

    // Case 2: A plunging orbit
    let plunging_orbit_params = GeodesicParameters {
        black_hole_mass,
        initial_state: [10.0, 0.0, 0.0, 0.02],
        proper_time_end: 500.0,
        initial_dt: 0.1,
    };

    // Case 3: A photon orbit (light bending)
    // For photons, we use a different equation, but can approximate with a high L/E ratio.
    let photon_orbit_params = GeodesicParameters {
        black_hole_mass,
        initial_state: [10.0, -1.0, 0.0, 0.03], // Start with inward velocity
        proper_time_end: 50.0,
        initial_dt: 0.01,
    };

    let orbits = vec![
        ("stable_orbit", stable_orbit_params),
        ("plunging_orbit", plunging_orbit_params),
        ("photon_orbit", photon_orbit_params),
    ];

    for (name, params) in orbits {
        println!("Simulating {}...", name);
        let path = run_geodesic_simulation(&params);

        let filename = format!("orbit_{}.csv", name);
        let mut file = File::create(&filename).unwrap();
        writeln!(file, "x,y").unwrap();
        path.iter().for_each(|(x, y)| {
            writeln!(file, "{},{}", x, y).unwrap();
        });
        println!("Saved path to {}", filename);
    }
}
