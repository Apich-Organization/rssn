//! # FDTD electrodynamics (2D TM)
//!
//! | Equations | Dimension | Boundaries | Sources |
//! |---|---|---|---|
//! | Maxwell TM (`Ez, Hx, Hy`) | 2D | ten-cell damping belt | soft point source |
//!
//! [`crate::sim::yee_fdtd`] is the general Yee solver in 1D/2D/3D with PEC,
//! periodic and convolutional PML faces, dielectrics, current sources and
//! TF/SF plane waves.

use ndarray::Array1;
use ndarray::Array2;
use serde::Deserialize;
use serde::Serialize;

use crate::io::write_npy_file;
use crate::kernels::par::for_chunks_mut;

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
    /// The source frequency in cycles per time step. Integer shifts alias.
    pub source_freq: f64,
}

const DAMPING_WIDTH: usize = 10;
const DAMPING_FACTOR: f64 = 0.95;
const PARALLEL_CELLS: usize = 8_192;

/// Runs the two-dimensional transverse-magnetic model from zero fields.
///
/// The returned arrays use shape `(width, height)` and record the electric
/// field after time steps 1, 6, 11, and so on. The soft point source is added
/// before the ten-cell boundary belt damps the electric field.
///
/// # Panics
/// Panics if the source lies outside the grid or its dimensions overflow
/// addressable field storage.
#[must_use]
pub fn run_fdtd_simulation(params: &FdtdParameters) -> Vec<Array2<f64>> {
    let size = fdtd_grid_size(params).unwrap_or_else(|error| panic!("{error}"));
    let mut snapshots = Vec::new();
    if params.time_steps != 0 {
        drop(propagate_fdtd(params, size, |step, ez| {
            if step % 5 == 0 {
                let shape = (params.width, params.height);
                snapshots.push(
                    Array2::from_shape_vec(shape, ez.to_vec())
                        .unwrap_or_else(|_| Array2::zeros(shape)),
                );
            }
        }));
    }
    snapshots
}

/// Returns the electric field after the last requested time step.
///
/// The model starts from zero fields and records only the final array in
/// `(width, height)` order. Zero steps return the unexcited field. The
/// half-Courant update, point source and boundary damping are shared with
/// [`run_fdtd_simulation`].
///
/// # Errors
/// Returns an error for an out-of-grid source, oversized dimensions or a
/// non-finite frequency.
pub fn run_fdtd_final_field(params: &FdtdParameters) -> Result<Array2<f64>, String> {
    let size = fdtd_grid_size(params)?;
    if !params.source_freq.is_finite() {
        return Err("FDTD source frequency must be finite".into());
    }
    let field = propagate_fdtd(params, size, |_, _| {});
    Array2::from_shape_vec((params.width, params.height), field).map_err(|error| error.to_string())
}

fn fdtd_grid_size(params: &FdtdParameters) -> Result<usize, &'static str> {
    if params.source_pos.0 >= params.width || params.source_pos.1 >= params.height {
        return Err("FDTD source position lies outside the grid");
    }
    params
        .width
        .checked_mul(params.height)
        .filter(|&size| size <= isize::MAX as usize / size_of::<[f64; 2]>())
        .ok_or("FDTD grid dimensions overflow addressable storage")
}

fn propagate_fdtd(
    params: &FdtdParameters,
    size: usize,
    mut record: impl FnMut(usize, &[f64]),
) -> Vec<f64> {
    let (nx, ny) = (params.width, params.height);
    let (source_row, source_column) = params.source_pos;
    let mut ez = vec![0.0_f64; size];
    if params.time_steps == 0 {
        return ez;
    }

    // The source is sampled at integer steps, so integer shifts in cycles
    // per step alias exactly. Reduce before constructing the phase.
    let angular_frequency = std::f64::consts::TAU * params.source_freq.fract();
    let pulse = |step| {
        let time = step as f64;
        (-((time - 30.0).powi(2)) / 100.0).exp() * (angular_frequency * time).sin()
    };
    if nx < 3 || ny < 3 {
        // There are no electric interior cells and every site is damped.
        let source = source_row * ny + source_column;
        for step in 0..params.time_steps {
            ez[source] += pulse(step);
            ez[source] *= DAMPING_FACTOR;
            record(step, &ez);
        }
        return ez;
    }

    // Each magnetic cell stores Hx and Hy, both computed from the same
    // two electric rows before any electric cell is advanced.
    let mut magnetic = vec![[0.0_f64; 2]; size];
    let rows_per_parallel_job = PARALLEL_CELLS.div_ceil(ny).max(2);
    let last_magnetic = (nx - 1) * ny;
    let undamped_rows = DAMPING_WIDTH..nx.saturating_sub(DAMPING_WIDTH);
    let damp = |field: &mut [f64]| field.iter_mut().for_each(|v| *v *= DAMPING_FACTOR);

    for step in 0..params.time_steps {
        for_chunks_mut(
            &mut magnetic[..last_magnetic],
            ny,
            rows_per_parallel_job,
            |i, row| {
                let previous = &ez[i * ny..(i + 1) * ny];
                let next = &ez[(i + 1) * ny..(i + 2) * ny];
                for (j, field) in row.iter_mut().take(ny - 1).enumerate() {
                    field[0] -= 0.5 * (previous[j + 1] - previous[j]);
                    field[1] += 0.5 * (next[j] - previous[j]);
                }
            },
        );

        let excitation = pulse(step);
        for_chunks_mut(&mut ez, ny, rows_per_parallel_job, |i, row| {
            if i != 0 && i < nx - 1 {
                let current = &magnetic[i * ny..(i + 1) * ny];
                let previous = &magnetic[(i - 1) * ny..i * ny];
                for j in 1..ny - 1 {
                    row[j] += 0.5
                        * (current[j][1] - previous[j][1] - (current[j][0] - current[j - 1][0]));
                }
            }
            if i == source_row {
                row[source_column] += excitation;
            }
            if !undamped_rows.contains(&i) || ny <= 2 * DAMPING_WIDTH {
                damp(row);
            } else {
                damp(&mut row[..DAMPING_WIDTH]);
                damp(&mut row[ny - DAMPING_WIDTH..]);
            }
        });
        record(step, &ez);
    }
    ez
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
pub use crate::sim::yee_fdtd::Yee;
/// General Yee solver (1D–3D, PEC/periodic/PML faces, TF/SF): see [`crate::sim::yee_fdtd`].
pub use crate::sim::yee_fdtd::YeeConfig;
