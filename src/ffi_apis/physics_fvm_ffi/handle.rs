//! Handle-based FFI API for physics FVM functions.

use std::ptr;

use crate::physics::physics_fvm::Mesh;
use crate::physics::physics_fvm::{
    self,
};

/// Creates a new Mesh handle.
#[unsafe(no_mangle)]

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
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_physics_fvm_mesh_free(
    mesh: *mut Mesh
) { unsafe {

    if !mesh.is_null() {

        let _ = Box::from_raw(mesh);
    }
}}

/// Returns a pointer to the mesh data.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_physics_fvm_mesh_data(
    mesh: *mut Mesh
) -> *mut f64 { unsafe {

    if mesh.is_null() {

        return ptr::null_mut();
    }

    (*mesh)
        .cells
        .as_mut_ptr()
        .cast::<f64>()
}}

/// Simulates 1D advection and returns the final values in a new buffer.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_physics_fvm_simulate_advection_1d(
) -> *mut f64 {

    let result = physics_fvm::simulate_1d_advection_scenario();

    let mut result =
        result.into_boxed_slice();

    let ptr = result.as_mut_ptr();

    std::mem::forget(result);

    ptr
}
