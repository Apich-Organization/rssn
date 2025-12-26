//! Handle-based FFI API for physics sim GPE superfluidity functions.

use crate::numerical::matrix::Matrix;
use crate::physics::physics_sim::gpe_superfluidity::GpeParameters;
use crate::physics::physics_sim::gpe_superfluidity::{
    self,
};

/// Runs the GPE ground state finder and returns the result as a Matrix handle (Nx x Ny).
#[no_mangle]

pub extern "C" fn rssn_physics_sim_gpe_run_ground_state_finder(
    nx: usize,
    ny: usize,
    lx: f64,
    ly: f64,
    d_tau: f64,
    time_steps: usize,
    g: f64,
    trap_strength: f64,
) -> *mut Matrix<f64> {

    let params = GpeParameters {
        nx,
        ny,
        lx,
        ly,
        d_tau,
        time_steps,
        g,
        trap_strength,
    };

    match gpe_superfluidity::run_gpe_ground_state_finder(&params) {
        | Ok(res) => {

            let rows = res.nrows();

            let cols = res.ncols();

            Box::into_raw(Box::new(
                Matrix::new(
                    rows,
                    cols,
                    res.into_raw_vec(),
                ),
            ))
        },
        | Err(_) => std::ptr::null_mut(),
    }
}
