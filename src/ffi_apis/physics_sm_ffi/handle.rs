//! Handle-based FFI API for physics SM (Spectral Methods) functions.

use crate::numerical::matrix::Matrix;
use crate::physics::physics_sm;

/// Simulates the 1D advection-diffusion scenario and returns the results as a Matrix handle (1xN).
#[unsafe(no_mangle)]

pub extern "C" fn rssn_physics_sm_simulate_1d_advection(
) -> *mut Matrix<f64> {

    let results = physics_sm::simulate_1d_advection_diffusion_scenario();

    let n = results.len();

    Box::into_raw(Box::new(
        Matrix::new(1, n, results),
    ))
}

/// Simulates the 2D advection-diffusion scenario and returns the results as a Matrix handle (`WxH`).
#[unsafe(no_mangle)]

pub extern "C" fn rssn_physics_sm_simulate_2d_advection(
) -> *mut Matrix<f64> {

    let results = physics_sm::simulate_2d_advection_diffusion_scenario();

    // Assuming NxN for simplicity in the scenario, let's just make it a square matrix if it fits,
    // otherwise a single row. The scenario uses 64x64.
    let n = results.len();

    let dim =
        (n as f64).sqrt() as usize;

    if dim * dim == n {

        Box::into_raw(Box::new(
            Matrix::new(
                dim,
                dim,
                results,
            ),
        ))
    } else {

        Box::into_raw(Box::new(
            Matrix::new(1, n, results),
        ))
    }
}
