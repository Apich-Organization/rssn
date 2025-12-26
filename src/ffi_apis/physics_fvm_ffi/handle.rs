//! Handle-based FFI API for physics FVM functions.

use crate::physics::physics_fvm::{
    self,
    Cell,
    Mesh,
};
use std::ptr;

/// Creates a new Mesh handle.
#[no_mangle]

pub extern "C" fn rssn_physics_fvm_mesh_new(
    num_cells: usize,
    domain_size: f64,
) -> *mut Mesh {

    Box::into_raw(Box::new(Mesh::new(
        num_cells,
        domain_size,
        |_| 0.0,
    )))
}

/// Frees a Mesh handle.
#[no_mangle]

pub unsafe extern "C" fn rssn_physics_fvm_mesh_free(mesh: *mut Mesh) {

    if !mesh.is_null() {

        let _ = Box::from_raw(mesh);
    }
}

/// Returns a pointer to the mesh data.
#[no_mangle]

pub unsafe extern "C" fn rssn_physics_fvm_mesh_data(mesh: *mut Mesh) -> *mut f64 {

    if mesh.is_null() {

        return ptr::null_mut();
    }

    (*mesh)
        .cells
        .as_mut_ptr() as *mut f64
}

/// Simulates 1D advection and returns the final values in a new buffer.
#[no_mangle]

pub extern "C" fn rssn_physics_fvm_simulate_advection_1d() -> *mut f64 {

    let result = physics_fvm::simulate_1d_advection_scenario();

    let mut result = result.into_boxed_slice();

    let ptr = result.as_mut_ptr();

    std::mem::forget(result);

    ptr
}
