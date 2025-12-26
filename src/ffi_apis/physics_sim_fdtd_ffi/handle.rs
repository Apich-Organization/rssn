//! Handle-based FFI API for physics sim FDTD electrodynamics functions.

use crate::numerical::matrix::Matrix;
use crate::physics::physics_sim::fdtd_electrodynamics::{self, FdtdParameters};

/// Runs a 2D FDTD simulation and returns the final Ez field as a Matrix handle (WxH).
#[no_mangle]

pub extern "C" fn rssn_physics_sim_fdtd_run_2d(
    width: usize,
    height: usize,
    time_steps: usize,
    source_x: usize,
    source_y: usize,
    source_freq: f64,
) -> *mut Matrix<f64> {

    let params = FdtdParameters {
        width,
        height,
        time_steps,
        source_pos: (source_x, source_y),
        source_freq,
    };

    let snapshots = fdtd_electrodynamics::run_fdtd_simulation(&params);

    if let Some(final_ez) = snapshots.last() {

        let rows = final_ez.nrows();

        let cols = final_ez.ncols();

        let data = final_ez
            .clone()
            .into_raw_vec();

        Box::into_raw(Box::new(Matrix::new(rows, cols, data)))
    } else {

        std::ptr::null_mut()
    }
}
