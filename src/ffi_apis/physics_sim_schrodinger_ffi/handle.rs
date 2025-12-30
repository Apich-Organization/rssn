//! Handle-based FFI API for physics sim Schrodinger quantum functions.

use num_complex::Complex;

use crate::numerical::matrix::Matrix;
use crate::physics::physics_sim::schrodinger_quantum::SchrodingerParameters;
use crate::physics::physics_sim::schrodinger_quantum::{
    self,
};

/// Runs a Schrodinger simulation and returns the final probability density as a Matrix handle (`NxxNy`).
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_physics_sim_schrodinger_run_2d(
    nx: usize,
    ny: usize,
    lx: f64,
    ly: f64,
    dt: f64,
    time_steps: usize,
    hbar: f64,
    mass: f64,
    potential_ptr: *const f64,
    initial_psi_re_ptr: *const f64,
    initial_psi_im_ptr: *const f64,
) -> *mut Matrix<f64> { unsafe {

    if potential_ptr.is_null()
        || initial_psi_re_ptr.is_null()
        || initial_psi_im_ptr.is_null()
    {

        return std::ptr::null_mut();
    }

    let n = nx * ny;

    let potential =
        std::slice::from_raw_parts(
            potential_ptr,
            n,
        )
        .to_vec();

    let re = std::slice::from_raw_parts(
        initial_psi_re_ptr,
        n,
    );

    let im = std::slice::from_raw_parts(
        initial_psi_im_ptr,
        n,
    );

    let mut initial_psi: Vec<
        Complex<f64>,
    > = re
        .iter()
        .zip(im.iter())
        .map(|(&r, &i)| {

            Complex::new(r, i)
        })
        .collect();

    let params =
        SchrodingerParameters {
            nx,
            ny,
            lx,
            ly,
            dt,
            time_steps,
            hbar,
            mass,
            potential,
        };

    match schrodinger_quantum::run_schrodinger_simulation(
        &params,
        &mut initial_psi,
    ) {
        | Ok(snapshots) => {
            if let Some(final_state) = snapshots.last() {

                let rows = final_state.nrows();

                let cols = final_state.ncols();

                let (data, _offset) = final_state
                    .clone()
                    .into_raw_vec_and_offset();

                Box::into_raw(Box::new(
                    Matrix::new(rows, cols, data),
                ))
            } else {

                std::ptr::null_mut()
            }
        },
        | Err(_) => std::ptr::null_mut(),
    }
}}
