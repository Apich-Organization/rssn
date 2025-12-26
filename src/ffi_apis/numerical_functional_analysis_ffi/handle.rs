//! Handle-based FFI API for numerical functional analysis.

use std::ptr;

use crate::ffi_apis::ffi_api::update_last_error;
use crate::numerical::functional_analysis;

/// Helper to convert raw pointers to (x, y) tuples.

unsafe fn to_points(
    x: *const f64,
    y: *const f64,
    len: usize,
) -> Vec<(f64, f64)> {

    let x_slice =
        std::slice::from_raw_parts(
            x, len,
        );

    let y_slice =
        std::slice::from_raw_parts(
            y, len,
        );

    x_slice
        .iter()
        .zip(y_slice.iter())
        .map(|(&xi, &yi)| (xi, yi))
        .collect()
}

/// Calculates the L1 norm.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_fa_l1_norm(
    x: *const f64,
    y: *const f64,
    len: usize,
) -> f64 {

    if x.is_null() || y.is_null() {

        return 0.0;
    }

    let points = to_points(x, y, len);

    functional_analysis::l1_norm(
        &points,
    )
}

/// Calculates the L2 norm.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_fa_l2_norm(
    x: *const f64,
    y: *const f64,
    len: usize,
) -> f64 {

    if x.is_null() || y.is_null() {

        return 0.0;
    }

    let points = to_points(x, y, len);

    functional_analysis::l2_norm(
        &points,
    )
}

/// Calculates the L-infinity norm.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_fa_infinity_norm(
    x: *const f64,
    y: *const f64,
    len: usize,
) -> f64 {

    if x.is_null() || y.is_null() {

        return 0.0;
    }

    let points = to_points(x, y, len);

    functional_analysis::infinity_norm(
        &points,
    )
}

/// Calculates the inner product.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_fa_inner_product(
    x1: *const f64,
    y1: *const f64,
    len1: usize,
    x2: *const f64,
    y2: *const f64,
    len2: usize,
    result: *mut f64,
) -> i32 {

    if x1.is_null()
        || y1.is_null()
        || x2.is_null()
        || y2.is_null()
        || result.is_null()
    {

        return -1;
    }

    let p1 = to_points(x1, y1, len1);

    let p2 = to_points(x2, y2, len2);

    match functional_analysis::inner_product(&p1, &p2) {
        | Ok(val) => {

            *result = val;

            0
        },
        | Err(e) => {

            update_last_error(e);

            -1
        },
    }
}
