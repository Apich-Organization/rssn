//! # FDTD electrodynamics (2D TM)
//!
//! | Equations | Dimension | Boundaries | Sources |
//! |---|---|---|---|
//! | Maxwell TM (`Ez, Hx, Hy`) | 2D | grid edge | point source |
//!
//! [`crate::sim::yee_fdtd`] is the general Yee solver in 1D/2D/3D with PEC,
//! periodic and convolutional PML faces, dielectrics, current sources and
//! TF/SF plane waves.

use ndarray::Array1;
use ndarray::Array2;
use rayon::prelude::*;
use serde::Deserialize;
use serde::Serialize;

use crate::io::write_npy_file;

/// Parameters for the FDTD simulation.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FdtdParameters {
    /// The width of the simulation grid.
    pub width: usize,
    /// The height of the simulation grid.
    pub height: usize,
    /// The number of time steps to simulate.
    pub time_steps: usize,
    /// The position of the source.
    pub source_pos: (usize, usize),
    /// The frequency of the source.
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
///
/// # Panics
/// Panics if `source_pos` lies outside the grid.
#[must_use]
pub fn run_fdtd_simulation(params: &FdtdParameters) -> Vec<Array2<f64>> {
    let (nx, ny) = (params.width, params.height);

    assert!(
        params.source_pos.0 < nx && params.source_pos.1 < ny,
        "FDTD source position lies outside the grid"
    );

    // Row-major flat buffers: row `i` is the chunk `[i * ny, (i + 1) * ny)`.
    let mut ez = vec![0.0_f64; nx * ny];

    let mut hx = vec![0.0_f64; nx * ny];

    let mut hy = vec![0.0_f64; nx * ny];

    let pml_thickness = 10;

    let mut snapshots = Vec::new();

    let row_len = ny.max(1);

    for t in 0..params.time_steps {
        // Update Hx and Hy (Magnetic field): row `i` only reads Ez and writes its own H row.
        hx.par_chunks_mut(row_len).zip(ez.par_chunks(row_len)).take(nx.saturating_sub(1)).for_each(
            |(hx_row, ez_row)| {
                for j in 0..ny - 1 {
                    hx_row[j] -= 0.5 * (ez_row[j + 1] - ez_row[j]);
                }
            },
        );

        let ez_next_rows = ez.par_chunks(row_len).skip(1);

        hy.par_chunks_mut(row_len)
            .zip(ez.par_chunks(row_len))
            .zip(ez_next_rows)
            .for_each(|((hy_row, ez_row), ez_next)| {
                for j in 0..ny - 1 {
                    hy_row[j] += 0.5 * (ez_next[j] - ez_row[j]);
                }
            });

        // Update Ez (Electric field) from the freshly updated H rows (interior only).
        ez.par_chunks_mut(row_len)
            .enumerate()
            .skip(1)
            .take(nx.saturating_sub(2))
            .for_each(|(i, ez_row)| {
                let hy_row = &hy[i * ny..(i + 1) * ny];

                let hy_prev_row = &hy[(i - 1) * ny..i * ny];

                let hx_row = &hx[i * ny..(i + 1) * ny];

                for j in 1..ny - 1 {
                    ez_row[j] += 0.5 * (hy_row[j] - hy_prev_row[j] - (hx_row[j] - hx_row[j - 1]));
                }
            });

        // Add soft source
        let pulse = (-((t as f64 - 30.0).powi(2)) / 100.0).exp()
            * (2.0 * std::f64::consts::PI * params.source_freq * (t as f64)).sin();

        ez[params.source_pos.0 * ny + params.source_pos.1] += pulse;

        // Apply PML/Boundary Damping in parallel
        ez.par_chunks_mut(row_len).enumerate().for_each(|(i, ez_row)| {
            for (j, v) in ez_row.iter_mut().enumerate() {
                if i < pml_thickness
                    || i >= nx - pml_thickness
                    || j < pml_thickness
                    || j >= ny - pml_thickness
                {
                    *v *= 0.95;
                }
            }
        });

        if t % 5 == 0 {
            snapshots.push(
                Array2::from_shape_vec((nx, ny), ez.clone())
                    .unwrap_or_else(|_| Array2::zeros((nx, ny))),
            );
        }
    }

    snapshots
}

/// An example scenario that runs an FDTD simulation and saves the final state.
///
/// # Errors
///
/// This function will return an error if the `ez` array cannot be reshaped
/// or if the final state cannot be written to the file.
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

/// General Yee solver (1D–3D, PEC/periodic/PML faces, TF/SF): see [`crate::sim::yee_fdtd`].
pub use crate::sim::yee_fdtd::{Yee, YeeConfig};
