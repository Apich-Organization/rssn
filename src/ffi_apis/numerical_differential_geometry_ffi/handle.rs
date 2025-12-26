//! Handle-based FFI API for numerical differential geometry.

use std::ptr;

use crate::ffi_apis::ffi_api::update_last_error;
use crate::numerical::differential_geometry;
use crate::numerical::matrix::Matrix;
use crate::symbolic::coordinates::CoordinateSystem;

/// Computes the metric tensor at a given point.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_dg_metric_tensor(
    system: CoordinateSystem,
    point: *const f64,
    n_vars: usize,
) -> *mut Matrix<f64> {

    if point.is_null() {

        return ptr::null_mut();
    }

    let point_slice =
        std::slice::from_raw_parts(
            point,
            n_vars,
        );

    match differential_geometry::metric_tensor_at_point(
        system,
        point_slice,
    ) {
        | Ok(g) => {

            let rows = g.len();

            let cols =
                if rows > 0 { g[0].len() } else { 0 };

            let flattened: Vec<f64> = g
                .into_iter()
                .flatten()
                .collect();

            Box::into_raw(Box::new(
                Matrix::new(
                    rows, cols, flattened,
                ),
            ))
        },
        | Err(e) => {

            update_last_error(e);

            ptr::null_mut()
        },
    }
}

/// Computes the Christoffel symbols at a given point.
/// Returns a flattened vector of size dim^3.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_dg_christoffel_symbols(
    system: CoordinateSystem,
    point: *const f64,
    n_vars: usize,
) -> *mut Vec<f64> {

    if point.is_null() {

        return ptr::null_mut();
    }

    let point_slice =
        std::slice::from_raw_parts(
            point,
            n_vars,
        );

    match differential_geometry::christoffel_symbols(
        system,
        point_slice,
    ) {
        | Ok(c) => {

            let flattened: Vec<f64> = c
                .into_iter()
                .flatten()
                .flatten()
                .collect();

            Box::into_raw(Box::new(flattened))
        },
        | Err(e) => {

            update_last_error(e);

            ptr::null_mut()
        },
    }
}

/// Computes the Ricci tensor at a given point.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_dg_ricci_tensor(
    system: CoordinateSystem,
    point: *const f64,
    n_vars: usize,
) -> *mut Matrix<f64> {

    if point.is_null() {

        return ptr::null_mut();
    }

    let point_slice =
        std::slice::from_raw_parts(
            point,
            n_vars,
        );

    match differential_geometry::ricci_tensor(
        system,
        point_slice,
    ) {
        | Ok(r) => {

            let rows = r.len();

            let cols =
                if rows > 0 { r[0].len() } else { 0 };

            let flattened: Vec<f64> = r
                .into_iter()
                .flatten()
                .collect();

            Box::into_raw(Box::new(
                Matrix::new(
                    rows, cols, flattened,
                ),
            ))
        },
        | Err(e) => {

            update_last_error(e);

            ptr::null_mut()
        },
    }
}

/// Computes the Ricci scalar at a given point.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_dg_ricci_scalar(
    system: CoordinateSystem,
    point: *const f64,
    n_vars: usize,
    result: *mut f64,
) -> i32 {

    if point.is_null()
        || result.is_null()
    {

        return -1;
    }

    let point_slice =
        std::slice::from_raw_parts(
            point,
            n_vars,
        );

    match differential_geometry::ricci_scalar(
        system,
        point_slice,
    ) {
        | Ok(r) => {

            *result = r;

            0
        },
        | Err(e) => {

            update_last_error(e);

            -1
        },
    }
}
