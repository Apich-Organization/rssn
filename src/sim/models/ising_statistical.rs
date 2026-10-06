//! #  Ising Model Simulation
//!
//! This module implements a Monte Carlo simulation of the 2D Ising model, a mathematical model
//! of ferromagnetism in statistical mechanics. It uses the Metropolis-Hastings algorithm to explore
//! the configuration space of spins on a grid.
//!
//! # Overview
//!
//! The simulation models a lattice of magnetic spins that can be in one of two states (+1 or -1).
//! The system evolves towards thermal equilibrium by randomly flipping spins based on the energy cost
//! and the temperature (Metropolis criterion).
//!
//! Key features include:
//! - **Monte Carlo Method**: Uses the Metropolis algorithm to sample states according to the Boltzmann distribution.
//! - **Parallel Update**: Implements a checkerboard (Red-Black) update scheme to allow parallel updates of non-interacting spins.
//! - **Phase Transition**: Can simulate the system across a range of temperatures to observe the transition from ordered (ferromagnetic) to disordered (paramagnetic) phases.
//! - **Magnetization Tracking**: Calculates the average magnetization of the system to characterize the phase.
//!
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/ising_phase_transition.png)

use std::fmt::Write as OtherWrite;
use std::fs::File;
use std::io::Write;
use std::path::Path;

use ndarray::Array2;
use rand_v10::prelude::*;
use rand_v10::rngs::StdRng;
use rayon::prelude::*;
use serde::Deserialize;
use serde::Serialize;

use crate::io::write_npy_file;

/// Parameters for the Ising model simulation.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct IsingParameters {
    /// The width of the simulation grid.
    pub width: usize,
    /// The height of the simulation grid.
    pub height: usize,
    /// The temperature of the system.
    pub temperature: f64,
    /// The number of Monte Carlo steps to perform.
    pub mc_steps: usize,
}

/// Seed used by [`run_ising_simulation`], so that plain runs are reproducible.
pub const DEFAULT_ISING_SEED: u64 = 0x1517_0000_5EED;

/// Derives an independent RNG for one (step, colour, row) of the checkerboard sweep.
fn row_rng(
    seed: u64,
    step: usize,
    colour: u64,
    row: usize,
) -> StdRng {
    let mixed = seed
        ^ ((step as u64) * 2 + colour).wrapping_mul(0x9E37_79B9_7F4A_7C15)
        ^ (row as u64 + 1).wrapping_mul(0xBF58_476D_1CE4_E5B9);

    StdRng::seed_from_u64(mixed)
}

/// Runs an Ising model simulation with the fixed seed [`DEFAULT_ISING_SEED`].
///
/// Equivalent to [`run_ising_simulation_seeded`]; repeated calls with the same
/// parameters return identical results.
#[must_use]
pub fn run_ising_simulation(params: &IsingParameters) -> (Vec<i8>, f64) {
    run_ising_simulation_seeded(params, DEFAULT_ISING_SEED)
}

/// Runs an Ising model simulation whose random numbers all derive from `seed`.
///
/// The initial configuration and every row update use their own generator
/// derived from `seed`, so the result does not depend on thread scheduling.
#[must_use]
pub fn run_ising_simulation_seeded(
    params: &IsingParameters,
    seed: u64,
) -> (Vec<i8>, f64) {
    let mut local_rng = StdRng::seed_from_u64(seed);

    let mut grid: Vec<i8> = (0..params.width * params.height)
        .map(|_| {
            if local_rng.random::<bool>() {
                1
            } else {
                -1
            }
        })
        .collect();

    let n_spins = (params.width * params.height) as f64;

    let b = 1.0 / params.temperature;

    for step in 0..params.mc_steps {
        // Checkerboard update for parallelism
        let grid_ptr = grid.as_mut_ptr() as usize;

        let width = params.width;

        let height = params.height;

        // Red points
        (0..height).into_par_iter().for_each(|i| {
            let mut local_rng = row_rng(seed, step, 0, i);

            for j in 0..width {
                if (i + j) % 2 == 0 {
                    let idx = i * width + j;

                    unsafe {
                        let g = grid_ptr as *mut i8;

                        let top = *g.add(((i + height - 1) % height) * width + j);

                        let bottom = *g.add(((i + 1) % height) * width + j);

                        let left = *g.add(i * width + (j + width - 1) % width);

                        let right = *g.add(i * width + (j + 1) % width);

                        let sum_neighbors = f64::from(top + bottom + left + right);

                        let delta_e = 2.0 * f64::from(*g.add(idx)) * sum_neighbors;

                        if delta_e < 0.0 || local_rng.random::<f64>() < (-delta_e * b).exp() {
                            *g.add(idx) *= -1;
                        }
                    }
                }
            }
        });

        // Black points
        (0..height).into_par_iter().for_each(|i| {
            let mut local_rng = row_rng(seed, step, 1, i);

            for j in 0..width {
                if (i + j) % 2 != 0 {
                    let idx = i * width + j;

                    unsafe {
                        let g = grid_ptr as *mut i8;

                        let top = *g.add(((i + height - 1) % height) * width + j);

                        let bottom = *g.add(((i + 1) % height) * width + j);

                        let left = *g.add(i * width + (j + width - 1) % width);

                        let right = *g.add(i * width + (j + 1) % width);

                        let sum_neighbors = f64::from(top + bottom + left + right);

                        let delta_e = 2.0 * f64::from(*g.add(idx)) * sum_neighbors;

                        if delta_e < 0.0 || local_rng.random::<f64>() < (-delta_e * b).exp() {
                            *g.add(idx) *= -1;
                        }
                    }
                }
            }
        });
    }

    let magnetization: f64 = grid.par_iter().map(|&s| f64::from(s)).sum::<f64>() / n_spins;

    (grid, magnetization.abs())
}

/// An example scenario that simulates the Ising model across a range of temperatures
/// to observe the phase transition.
///
/// # Errors
///
/// This function will return an error if it fails to reshape the `Array2<f64>` for NPY
/// output or if it fails to create or write to the output files (CSV or NPY).
/// The files `ising_low_temp_state.npy`, `ising_high_temp_state.npy` and
/// `ising_magnetization_vs_temp.csv` are written into `output_dir`, which is
/// created if missing.
pub fn simulate_ising_phase_transition_scenario(output_dir: &Path) -> Result<(), String> {
    std::fs::create_dir_all(output_dir).map_err(|e| e.to_string())?;

    println!(
        "Running Ising model phase \
         transition simulation..."
    );

    let temperatures: Vec<f64> = (0..=40).map(|i| f64::from(i).mul_add(0.1, 0.1)).collect();

    let scenario_results: Vec<(f64, f64, Vec<i8>)> = temperatures
        .par_iter()
        .map(|&temp| {
            let params = IsingParameters {
                width: 50,
                height: 50,
                temperature: temp,
                mc_steps: 2000,
            };

            let (grid, mag) = run_ising_simulation(&params);

            (temp, mag, grid)
        })
        .collect();

    let mut results = String::from("temperature,magnetization\n");

    for (i, (temp, mag, grid)) in scenario_results.iter().enumerate() {
        writeln!(results, "{temp},{mag}").expect("String transition failed.");

        if i == 5 {
            let arr: Array2<f64> =
                Array2::from_shape_vec((50, 50), grid.iter().map(|&s| f64::from(s)).collect())
                    .map_err(|e| e.to_string())?;

            write_npy_file(output_dir.join("ising_low_temp_state.npy"), &arr)?;
        }

        if i == 35 {
            let arr: Array2<f64> =
                Array2::from_shape_vec((50, 50), grid.iter().map(|&s| f64::from(s)).collect())
                    .map_err(|e| e.to_string())?;

            write_npy_file(output_dir.join("ising_high_temp_state.npy"), &arr)?;
        }
    }

    let mut file = File::create(output_dir.join("ising_magnetization_vs_temp.csv"))
        .map_err(|e| e.to_string())?;

    file.write_all(results.as_bytes())
        .map_err(|e| e.to_string())?;

    println!(
        "Saved results to CSV and NPY \
         files."
    );

    Ok(())
}

/// Thermodynamic averages of a Monte Carlo run (per spin, `J = k_B = 1`).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct IsingObservables {
    /// `⟨E⟩ / N`.
    pub energy: f64,
    /// `⟨|M|⟩ / N`.
    pub magnetization: f64,
    /// Specific heat `(⟨E²⟩ - ⟨E⟩²) / (N T²)`.
    pub specific_heat: f64,
    /// Susceptibility `(⟨M²⟩ - ⟨|M|⟩²) / (N T)`.
    pub susceptibility: f64,
    /// Binder cumulant `1 - ⟨M⁴⟩ / (3 ⟨M²⟩²)`.
    pub binder: f64,
    /// Mean Wolff cluster size (fraction of the lattice).
    pub mean_cluster: f64,
}

fn ising_energy(
    grid: &[i8],
    width: usize,
    height: usize,
) -> f64 {
    let mut e = 0_i64;
    for i in 0..height {
        for j in 0..width {
            let s = i64::from(grid[i * width + j]);
            let right = i64::from(grid[i * width + (j + 1) % width]);
            let down = i64::from(grid[((i + 1) % height) * width + j]);
            e -= s * (right + down);
        }
    }
    e as f64
}

/// The Wolff single-cluster algorithm:
///
/// a cluster grown from a random seed spin, each aligned neighbour added
/// with probability `1 - exp(-2/T)`, is flipped as a whole. Rejection-free
/// and free of critical slowing down, so it samples near `T_c` far better
/// than local Metropolis updates. `sweeps` cluster updates are measured
/// after `thermalisation` discarded ones.
#[must_use]
pub fn run_wolff_simulation(
    width: usize,
    height: usize,
    temperature: f64,
    thermalisation: usize,
    sweeps: usize,
    seed: u64,
) -> IsingObservables {
    let n = width * height;
    let mut rng = StdRng::seed_from_u64(seed);
    let mut grid = vec![1_i8; n];
    let add_probability = 1.0 - (-2.0 / temperature).exp();
    let mut stack = Vec::with_capacity(n);
    let mut in_cluster = vec![false; n];
    let (mut e_sum, mut e2_sum, mut m_sum, mut m2_sum, mut m4_sum, mut cluster_sum) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    for step in 0..thermalisation + sweeps {
        let start = rng.random_range(0..n);
        let spin = grid[start];
        stack.clear();
        stack.push(start);
        in_cluster[start] = true;
        let mut members = Vec::new();
        while let Some(site) = stack.pop() {
            members.push(site);
            let (i, j) = (site / width, site % width);
            let neighbours = [
                ((i + height - 1) % height) * width + j,
                ((i + 1) % height) * width + j,
                i * width + (j + width - 1) % width,
                i * width + (j + 1) % width,
            ];
            for nb in neighbours {
                if !in_cluster[nb] && grid[nb] == spin && rng.random::<f64>() < add_probability {
                    in_cluster[nb] = true;
                    stack.push(nb);
                }
            }
        }
        for &site in &members {
            grid[site] = -grid[site];
            in_cluster[site] = false;
        }
        if step >= thermalisation {
            let e = ising_energy(&grid, width, height);
            let m = grid.iter().map(|&s| f64::from(s)).sum::<f64>().abs();
            e_sum += e;
            e2_sum += e * e;
            m_sum += m;
            m2_sum += m * m;
            m4_sum += m.powi(4);
            cluster_sum += members.len() as f64;
        }
    }
    let samples = sweeps.max(1) as f64;
    let nf = n as f64;
    let (e, e2, m, m2, m4) = (e_sum / samples, e2_sum / samples, m_sum / samples, m2_sum / samples, m4_sum / samples);
    IsingObservables {
        energy: e / nf,
        magnetization: m / nf,
        specific_heat: (e2 - e * e) / (nf * temperature * temperature),
        susceptibility: (m2 - m * m) / (nf * temperature),
        binder: 1.0 - m4 / (3.0 * m2 * m2),
        mean_cluster: cluster_sum / samples / nf,
    }
}

#[cfg(test)]
mod wolff_tests {
    use super::*;

    #[test]
    fn ordered_and_disordered_phases() {
        // T_c = 2/ln(1 + √2) ≈ 2.269.
        let cold = run_wolff_simulation(24, 24, 1.5, 200, 1000, 7);
        let hot = run_wolff_simulation(24, 24, 4.0, 200, 3000, 7);
        assert!(cold.magnetization > 0.95 && cold.binder > 0.6, "{cold:?}");
        assert!(hot.magnetization < 0.3 && hot.binder < 0.4, "{hot:?}");
        // Onsager's energy per spin at T = 1.5 is ≈ -1.947.
        assert!((cold.energy + 1.947).abs() < 0.02, "{cold:?}");
        // The specific heat peaks near T_c.
        let near = run_wolff_simulation(24, 24, 2.27, 200, 3000, 7);
        assert!(near.specific_heat > cold.specific_heat && near.specific_heat > hot.specific_heat, "{near:?}");
    }
}
