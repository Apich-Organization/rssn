// src/physics/physics_sim/fdtd_electrodynamics.rs
// 2D FDTD simulation for Maxwell's equations.

use crate::output::io::write_npy_file;
use ndarray::Array2;

/// Parameters for the FDTD simulation.
pub struct FdtdParameters {
    pub width: usize,
    pub height: usize,
    pub time_steps: usize,
    pub source_pos: (usize, usize),
    pub source_freq: f64,
}

/// Runs a 2D FDTD simulation for the Transverse Magnetic (TM) mode.
/// This solves for Ez, Hx, and Hy fields.
///
/// # Arguments
/// * `params` - The simulation parameters.
///
/// # Returns
/// A `Vec` containing snapshots of the Ez field at specified intervals.
pub fn run_fdtd_simulation(params: &FdtdParameters) -> Vec<Array2<f64>> {
    let (nx, ny) = (params.width, params.height);
    let mut ez = Array2::<f64>::zeros((nx, ny));
    let mut hx = Array2::<f64>::zeros((nx, ny));
    let mut hy = Array2::<f64>::zeros((nx, ny));

    // Simple PML (Perfectly Matched Layer) for absorbing boundary conditions
    let pml_thickness = 10;
    let _ez_x_pml = Array2::<f64>::zeros((pml_thickness, ny));
    let _ez_y_pml = Array2::<f64>::zeros((nx, pml_thickness));

    let mut snapshots = Vec::new();

    for t in 0..params.time_steps {
        // Update H fields
        for i in 0..nx - 1 {
            for j in 0..ny - 1 {
                hx[[i, j]] -= 0.5 * (ez[[i, j + 1]] - ez[[i, j]]);
                hy[[i, j]] += 0.5 * (ez[[i + 1, j]] - ez[[i, j]]);
            }
        }

        // Update E field
        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                ez[[i, j]] += 0.5 * (hy[[i, j]] - hy[[i - 1, j]] - (hx[[i, j]] - hx[[i, j - 1]]));
            }
        }

        // Add a simple sinusoidal source (a "Gaussian pulse")
        let pulse = (-((t as f64 - 30.0) * (t as f64 - 30.0)) / 100.0).exp()
            * (2.0 * std::f64::consts::PI * params.source_freq * (t as f64)).sin();
        ez[[params.source_pos.0, params.source_pos.1]] += pulse;

        // Apply simple absorbing boundary conditions by damping the edges
        for i in 0..nx {
            for j in 0..ny {
                if i < pml_thickness
                    || i >= nx - pml_thickness
                    || j < pml_thickness
                    || j >= ny - pml_thickness
                {
                    ez[[i, j]] *= 0.95; // Dampen fields at the boundary
                }
            }
        }

        // Save a snapshot at intervals
        if t % 5 == 0 {
            snapshots.push(ez.clone());
        }
    }

    snapshots
}

/// An example scenario that runs an FDTD simulation and saves the final state.
pub fn simulate_wave_propagation_scenario() {
    println!("Running FDTD wave propagation simulation...");

    let params = FdtdParameters {
        width: 200,
        height: 200,
        time_steps: 400,
        source_pos: (100, 100),
        source_freq: 0.1,
    };

    let snapshots = run_fdtd_simulation(&params);

    if let Some(final_state) = snapshots.last() {
        let filename = "fdtd_wave_final_state.npy";
        println!("Saving final state to {}", filename);
        write_npy_file(filename, final_state);
    } else {
        println!("Simulation produced no snapshots.");
    }
}
