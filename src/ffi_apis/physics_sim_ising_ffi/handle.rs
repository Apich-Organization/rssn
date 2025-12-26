//! Handle-based FFI API for physics sim Ising statistical functions.

use crate::numerical::matrix::Matrix;
use crate::physics::physics_sim::ising_statistical::{
    self,
    IsingParameters,
};

#[repr(C)]

pub struct IsingResultHandle {
    pub grid: *mut Matrix<f64>,
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

    let (grid, mag) =
        ising_statistical::run_ising_simulation(&params);

    let grid_f64: Vec<f64> = grid
        .into_iter()
        .map(|s| s as f64)
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
