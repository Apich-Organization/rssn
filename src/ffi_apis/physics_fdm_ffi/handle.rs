//! Handle-based FFI API for physics FDM functions.

use crate::physics::physics_fdm::{
    self,
    Dimensions,
    FdmGrid,
};
use std::ptr;

/// Creates a new FdmGrid handle with the given dimensions.
#[no_mangle]

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

/// Frees a FdmGrid handle.
#[no_mangle]

pub unsafe extern "C" fn rssn_physics_fdm_grid_free(grid: *mut FdmGrid<f64>) {

    if !grid.is_null() {

        let _ = Box::from_raw(grid);
    }
}

/// Returns the size of the grid data.
#[no_mangle]

pub unsafe extern "C" fn rssn_physics_fdm_grid_len(grid: *mut FdmGrid<f64>) -> usize {

    if grid.is_null() {

        return 0;
    }

    (*grid).len()
}

/// Returns a pointer to the grid data.
#[no_mangle]

pub unsafe extern "C" fn rssn_physics_fdm_grid_data(grid: *mut FdmGrid<f64>) -> *mut f64 {

    if grid.is_null() {

        return ptr::null_mut();
    }

    (*grid)
        .as_mut_slice()
        .as_mut_ptr()
}

/// Simulates 2D heat conduction and returns a new FdmGrid handle.
#[no_mangle]

pub extern "C" fn rssn_physics_fdm_simulate_heat_2d() -> *mut FdmGrid<f64> {

    Box::into_raw(Box::new(
        physics_fdm::simulate_2d_heat_conduction_scenario(),
    ))
}

/// Simulates 2D wave propagation and returns a new FdmGrid handle.
#[no_mangle]

pub extern "C" fn rssn_physics_fdm_simulate_wave_2d() -> *mut FdmGrid<f64> {

    Box::into_raw(Box::new(
        physics_fdm::simulate_2d_wave_propagation_scenario(),
    ))
}
