//! Handle-based FFI API for physics sim Navier-Stokes functions.

use crate::numerical::matrix::Matrix;
use crate::physics::physics_sim::navier_stokes_fluid;

/// Result handles for the Navier-Stokes simulation containing velocity and pressure fields.
///
/// This C-compatible struct encapsulates the output of an incompressible Navier-Stokes
/// solver, providing matrix handles to the computed velocity components and pressure.
#[repr(C)]

pub struct NavierStokesResultHandles {
    /// Pointer to a Matrix containing the horizontal velocity field u(x,y) in m/s.
    pub u: *mut Matrix<f64>,
    /// Pointer to a Matrix containing the vertical velocity field v(x,y) in m/s.
    pub v: *mut Matrix<f64>,
    /// Pointer to a Matrix containing the pressure field p(x,y) in Pa.
    pub p: *mut Matrix<f64>,
}

/// Runs the lid-driven cavity simulation and returns handles to the U, V, and P matrices.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_physics_sim_navier_stokes_run_lid_driven_cavity(
    nx: usize,
    ny: usize,
    re: f64,
    dt: f64,
    n_iter: usize,
    lid_velocity: f64,
) -> NavierStokesResultHandles {

    let params = navier_stokes_fluid::NavierStokesParameters {
        nx,
        ny,
        re,
        dt,
        n_iter,
        lid_velocity,
    };

    match navier_stokes_fluid::run_lid_driven_cavity(&params) {
        | Ok((u, v, p)) => {

            let u_mat = Matrix::new(
                ny,
                nx,
                u.into_raw_vec_and_offset().0,
            );

            let v_mat = Matrix::new(
                ny,
                nx,
                v.into_raw_vec_and_offset().0,
            );

            let p_mat = Matrix::new(
                ny,
                nx,
                p.into_raw_vec_and_offset().0,
            );

            NavierStokesResultHandles {
                u : Box::into_raw(Box::new(u_mat)),
                v : Box::into_raw(Box::new(v_mat)),
                p : Box::into_raw(Box::new(p_mat)),
            }
        },
        | Err(_) => {
            NavierStokesResultHandles {
                u : std::ptr::null_mut(),
                v : std::ptr::null_mut(),
                p : std::ptr::null_mut(),
            }
        },
    }
}

/// Frees the result handles.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_physics_sim_navier_stokes_free_results(
    handles: NavierStokesResultHandles
) { unsafe {

    if !handles.u.is_null() {

        let _ =
            Box::from_raw(handles.u);
    }

    if !handles.v.is_null() {

        let _ =
            Box::from_raw(handles.v);
    }

    if !handles.p.is_null() {

        let _ =
            Box::from_raw(handles.p);
    }
}}
