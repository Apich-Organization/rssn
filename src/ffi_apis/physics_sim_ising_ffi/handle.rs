//! Handle-based FFI API for physics sim Ising statistical functions.

use crate::numerical::matrix::Matrix;
use crate::physics::physics_sim::ising_statistical::IsingParameters;
use crate::physics::physics_sim::ising_statistical::{
    self,
};

/// Result handle for the Ising model simulation containing spin grid and magnetization.
///
/// This C-compatible struct encapsulates the output of an Ising model Monte Carlo
/// simulation, providing both the final spin configuration and the computed magnetization.
#[repr(C)]

pub struct IsingResultHandle {
    /// Pointer to a Matrix containing the final spin configuration as f64 values (±1.0).
    pub grid: *mut Matrix<f64>,
    /// Average magnetization M = ⟨∑ᵢsᵢ⟩/N, ranging from -1 (all spins down) to +1 (all spins up).
    pub magnetization: f64,
}

/// Runs a 2D Ising model simulation and returns the final grid as a Matrix handle and the magnetization.
#[no_mangle]

pub extern "C" fn rssn_physics_sim_ising_run(
    width: usize,
    height: usize,
    temperature: f64,
    mc_steps: usize,
) -> IsingResultHandle {

    let params = IsingParameters {
        width,
        height,
        temperature,
        mc_steps,
    };

    let (grid, mag) = ising_statistical::run_ising_simulation(&params);

    let grid_f64: Vec<f64> = grid
        .into_iter()
        .map(|s| f64::from(s))
        .collect();

    let matrix = Matrix::new(
        height,
        width,
        grid_f64,
    );

    IsingResultHandle {
        grid: Box::into_raw(Box::new(
            matrix,
        )),
        magnetization: mag,
    }
}

/// Frees the Ising result handle.
#[no_mangle]

pub unsafe extern "C" fn rssn_physics_sim_ising_free_result(
    handle: IsingResultHandle
) {

    if !handle
        .grid
        .is_null()
    {

        let _ =
            Box::from_raw(handle.grid);
    }
}
