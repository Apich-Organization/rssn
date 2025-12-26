use crate::output::io::write_npy_file;
use ndarray::Array2;
use rand::{thread_rng, Rng, SeedableRng};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::fmt::Write as OtherWrite;
use std::fs::File;
use std::io::Write;

/// Parameters for the Ising model simulation.
#[derive(Clone, Debug, Serialize, Deserialize)]

pub struct IsingParameters {
    pub width: usize,
    pub height: usize,
    pub temperature: f64,
    pub mc_steps: usize,
}

pub fn run_ising_simulation(params: &IsingParameters) -> (Vec<i8>, f64) {

    let mut rng = thread_rng();

    let mut grid: Vec<i8> = (0..params.width * params.height)
        .map(|_| {
            if rng.gen::<bool>() {

                1
            } else {

                -1
            }
        })
        .collect();

    let n_spins = (params.width * params.height) as f64;

    let b = 1.0 / params.temperature;

    for _ in 0..params.mc_steps {

        // Checkerboard update for parallelism
        let grid_ptr = grid.as_mut_ptr() as usize;

        let width = params.width;

        let height = params.height;

        // Red points
        (0..height)
            .into_par_iter()
            .for_each(|i| {

                let mut local_rng = thread_rng();

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

                            if delta_e < 0.0 || local_rng.gen::<f64>() < (-delta_e * b).exp() {

                                *g.add(idx) *= -1;
                            }
                        }
                    }
                }
            });

        // Black points
        (0..height)
            .into_par_iter()
            .for_each(|i| {

                let mut local_rng = thread_rng();

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

                            if delta_e < 0.0 || local_rng.gen::<f64>() < (-delta_e * b).exp() {

                                *g.add(idx) *= -1;
                            }
                        }
                    }
                }
            });
    }

    let magnetization: f64 = grid
        .par_iter()
        .map(|&s| f64::from(s))
        .sum::<f64>()
        / n_spins;

    (
        grid,
        magnetization.abs(),
    )
}

/// An example scenario that simulates the Ising model across a range of temperatures
/// to observe the phase transition.

pub fn simulate_ising_phase_transition_scenario() -> Result<(), String> {

    println!("Running Ising model phase transition simulation...");

    let temperatures: Vec<f64> = (0..=40)
        .map(|i| 0.1 + f64::from(i) * 0.1)
        .collect();

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

    for (i, (temp, mag, grid)) in scenario_results
        .iter()
        .enumerate()
    {

        writeln!(
            results,
            "{},{}",
            temp, mag
        )
        .expect("String transition failed.");

        if i == 5 {

            let arr: Array2<f64> = Array2::from_shape_vec(
                (50, 50),
                grid.iter()
                    .map(|&s| f64::from(s))
                    .collect(),
            )
            .map_err(|e| e.to_string())?;

            write_npy_file(
                "ising_low_temp_state.npy",
                &arr,
            )?;
        }

        if i == 35 {

            let arr: Array2<f64> = Array2::from_shape_vec(
                (50, 50),
                grid.iter()
                    .map(|&s| f64::from(s))
                    .collect(),
            )
            .map_err(|e| e.to_string())?;

            write_npy_file(
                "ising_high_temp_state.npy",
                &arr,
            )?;
        }
    }

    let mut file = File::create("ising_magnetization_vs_temp.csv").map_err(|e| e.to_string())?;

    file.write_all(results.as_bytes())
        .map_err(|e| e.to_string())?;

    println!("Saved results to CSV and NPY files.");

    Ok(())
}
