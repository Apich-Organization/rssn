// src/physics/physics_sim/ising_statistical.rs
// 2D Ising model simulation using the Metropolis Monte Carlo method.

use crate::output::io::write_npy_file;
use ndarray::Array2;
use rand::{thread_rng, Rng};
use std::fs::File;
use std::io::Write;

/// Parameters for the Ising model simulation.
pub struct IsingParameters {
    pub width: usize,
    pub height: usize,
    pub temperature: f64,
    pub mc_steps: usize, // Monte Carlo steps
}

/// Runs a 2D Ising model simulation for a given temperature.
///
/// # Arguments
/// * `params` - The simulation parameters.
///
/// # Returns
/// A tuple containing the final grid state (`Vec<i8>`) and the final magnetization.
pub fn run_ising_simulation(params: &IsingParameters) -> (Vec<i8>, f64) {
    let mut rng = thread_rng();
    let mut grid: Vec<i8> = (0..params.width * params.height)
        .map(|_| if rng.gen::<bool>() { 1 } else { -1 })
        .collect();

    let n_spins = (params.width * params.height) as f64;

    // Main Monte Carlo loop
    for _ in 0..params.mc_steps {
        for _ in 0..grid.len() {
            // One step = N trial flips
            // Pick a random spin
            let i = rng.gen_range(0..params.height);
            let j = rng.gen_range(0..params.width);
            let idx = i * params.width + j;

            // Calculate energy change if this spin is flipped
            // Energy E = -J * sum(s_i * s_j), with J=1
            // We use periodic boundary conditions.
            let top = grid[((i.wrapping_sub(1)) % params.height) * params.width + j];
            let bottom = grid[((i + 1) % params.height) * params.width + j];
            let left = grid[i * params.width + (j.wrapping_sub(1)) % params.width];
            let right = grid[i * params.width + (j + 1) % params.width];

            let sum_neighbors = (top + bottom + left + right) as f64;
            let delta_e = 2.0 * grid[idx] as f64 * sum_neighbors;

            // Metropolis acceptance criterion
            if delta_e < 0.0 || rng.gen::<f64>() < (-delta_e / params.temperature).exp() {
                grid[idx] *= -1; // Flip the spin
            }
        }
    }

    let magnetization: f64 = grid.iter().map(|&s| s as f64).sum::<f64>() / n_spins;
    (grid, magnetization.abs())
}

/// An example scenario that simulates the Ising model across a range of temperatures
/// to observe the phase transition.
pub fn simulate_ising_phase_transition_scenario() {
    println!("Running Ising model phase transition simulation...");

    let temperatures = (0..=40).map(|i| 0.1 + i as f64 * 0.1);
    let mut results = String::from("temperature,magnetization\n");

    for (i, temp) in temperatures.enumerate() {
        println!("Simulating at T = {:.2}", temp);
        let params = IsingParameters {
            width: 50,
            height: 50,
            temperature: temp,
            mc_steps: 2000, // Equilibration steps + measurement steps
        };

        let (grid, mag) = run_ising_simulation(&params);
        results.push_str(&format!("{},{}\n", temp, mag));

        // Save grid state for a low and a high temperature
        if i == 5 {
            // T = 1.5 (ordered)
            let arr: Array2<f64> = Array2::from_shape_vec(
                (params.height, params.width),
                grid.iter().map(|&s| s as f64).collect(),
            )
            .unwrap();
            write_npy_file("ising_low_temp_state.npy", &arr);
            println!("Saved low temperature state to ising_low_temp_state.npy");
        }
        if i == 35 {
            // T = 3.6 (disordered)
            let arr: Array2<f64> = Array2::from_shape_vec(
                (params.height, params.width),
                grid.iter().map(|&s| s as f64).collect(),
            )
            .unwrap();
            write_npy_file("ising_high_temp_state.npy", &arr);
            println!("Saved high temperature state to ising_high_temp_state.npy");
        }
    }

    // Save magnetization data
    let mut file = File::create("ising_magnetization_vs_temp.csv").unwrap();
    file.write_all(results.as_bytes()).unwrap();
    println!("Saved magnetization data to ising_magnetization_vs_temp.csv");
}
