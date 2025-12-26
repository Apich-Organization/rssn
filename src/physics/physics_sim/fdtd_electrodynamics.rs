use crate::output::io::write_npy_file;
use ndarray::{Array1, Array2};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// Parameters for the FDTD simulation.
#[derive(Clone, Debug, Serialize, Deserialize)]
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
    let pml_thickness = 10;
    let mut snapshots = Vec::new();

    for t in 0..params.time_steps {
        // Update Hx and Hy (Magnetic field)
        let ez_ptr = ez.as_ptr() as usize;
        let hx_ptr = hx.as_mut_ptr() as usize;
        let hy_ptr = hy.as_mut_ptr() as usize;

        // Parallel update for Hx
        (0..nx).into_par_iter().for_each(|i| {
            if i < nx - 1 {
                for j in 0..ny - 1 {
                    unsafe {
                        let ez = ez_ptr as *const f64;
                        let hx = hx_ptr as *mut f64;
                        let val_j1 = *ez.add(i * ny + (j + 1));
                        let val_j = *ez.add(i * ny + j);
                        *hx.add(i * ny + j) -= 0.5 * (val_j1 - val_j);
                    }
                }
            }
        });

        // Parallel update for Hy
        (0..nx).into_par_iter().for_each(|i| {
            if i < nx - 1 {
                for j in 0..ny - 1 {
                    unsafe {
                        let ez = ez_ptr as *const f64;
                        let hy = hy_ptr as *mut f64;
                        let val_i1 = *ez.add((i + 1) * ny + j);
                        let val_j = *ez.add(i * ny + j);
                        *hy.add(i * ny + j) += 0.5 * (val_i1 - val_j);
                    }
                }
            }
        });

        // Update Ez (Electric field)
        let ez_mut_ptr = ez.as_mut_ptr() as usize;
        let hx_const_ptr = hx.as_ptr() as usize;
        let hy_const_ptr = hy.as_ptr() as usize;
        (1..nx - 1).into_par_iter().for_each(|i| {
            for j in 1..ny - 1 {
                unsafe {
                    let ez = ez_mut_ptr as *mut f64;
                    let hx = hx_const_ptr as *const f64;
                    let hy = hy_const_ptr as *const f64;
                    let hy_val = *hy.add(i * ny + j);
                    let hy_prev = *hy.add((i - 1) * ny + j);
                    let hx_val = *hx.add(i * ny + j);
                    let hx_prev = *hx.add(i * ny + (j - 1));
                    *ez.add(i * ny + j) += 0.5 * (hy_val - hy_prev - (hx_val - hx_prev));
                }
            }
        });

        // Add soft source
        let pulse = (-((t as f64 - 30.0).powi(2)) / 100.0).exp()
            * (2.0 * std::f64::consts::PI * params.source_freq * (t as f64)).sin();
        ez[[params.source_pos.0, params.source_pos.1]] += pulse;

        // Apply PML/Boundary Damping in parallel
        (0..nx).into_par_iter().for_each(|i| {
            for j in 0..ny {
                if i < pml_thickness
                    || i >= nx - pml_thickness
                    || j < pml_thickness
                    || j >= ny - pml_thickness
                {
                    unsafe {
                        let ez = ez_mut_ptr as *mut f64;
                        *ez.add(i * ny + j) *= 0.95;
                    }
                }
            }
        });

        if t % 5 == 0 {
            snapshots.push(ez.clone());
        }
    }
    snapshots
}
/// An example scenario that runs an FDTD simulation and saves the final state.
pub fn simulate_and_save_final_state(
    grid_size: usize,
    time_steps: usize,
    filename: &str,
) -> Result<(), String> {
    let mut ez = Array1::<f64>::zeros(grid_size);
    let mut hy = Array1::<f64>::zeros(grid_size - 1);
    for _ in 0..time_steps {
        for i in 0..grid_size - 1 {
            hy[i] += 0.5 * (ez[i + 1] - ez[i]);
        }
        for i in 1..grid_size - 1 {
            ez[i] += 0.5 * (hy[i] - hy[i - 1]);
        }
    }
    let final_state_2d = ez
        .into_shape_with_order((grid_size, 1))
        .map_err(|e| e.to_string())?;
    write_npy_file(filename, &final_state_2d)?;
    Ok(())
}
