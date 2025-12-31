//! Handle-based FFI API for physics FDM functions.

use std::ptr;

use crate::physics::physics_fdm::Dimensions;
use crate::physics::physics_fdm::FdmGrid;
use crate::physics::physics_fdm::{
    self,
};

/// Creates a new `FdmGrid` handle with the given dimensions.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_physics_fdm_grid_new(
    d1: usize,
    d2: usize,
    d3: usize,
) -> *mut FdmGrid<f64> {

    let dims = if d3 > 0 {

        Dimensions::D3(d1, d2, d3)
    } else if d2 > 0 {

        Dimensions::D2(d1, d2)
    } else {

        Dimensions::D1(d1)
    };

    Box::into_raw(Box::new(
        FdmGrid::new(dims),
    ))
}

/// Frees a `FdmGrid` handle.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_physics_fdm_grid_free(
    grid: *mut FdmGrid<f64>
) {

    unsafe {

        if !grid.is_null() {

            let _ = Box::from_raw(grid);
        }
    }
}

/// Returns the size of the grid data.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_physics_fdm_grid_len(
    grid: *mut FdmGrid<f64>
) -> usize {

    unsafe {

        if grid.is_null() {

            return 0;
        }

        (*grid).len()
    }
}

/// Returns a pointer to the grid data.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_physics_fdm_grid_data(
    grid: *mut FdmGrid<f64>
) -> *mut f64 {

    unsafe {

        if grid.is_null() {

            return ptr::null_mut();
        }

        (*grid)
            .as_mut_slice()
            .as_mut_ptr()
    }
}

/// Simulates 2D heat conduction and returns a new `FdmGrid` handle.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_physics_fdm_simulate_heat_2d()
-> *mut FdmGrid<f64> {

    Box::into_raw(Box::new(
        physics_fdm::simulate_2d_heat_conduction_scenario(),
    ))
}

/// Simulates 2D wave propagation and returns a new `FdmGrid` handle.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_physics_fdm_simulate_wave_2d()
-> *mut FdmGrid<f64> {

    Box::into_raw(Box::new(
        physics_fdm::simulate_2d_wave_propagation_scenario(),
    ))
}
