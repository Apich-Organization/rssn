//! Handle-based FFI API for numerical interpolation.

use std::ptr;
use std::sync::Arc;

use crate::ffi_apis::ffi_api::update_last_error;
use crate::numerical::interpolate;
use crate::numerical::polynomial::Polynomial;

/// Opaque type for cubic spline closure.

pub type CubicSplineHandle =
    Arc<dyn Fn(f64) -> f64>;

/// Computes Lagrange interpolation and returns a Polynomial pointer.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_lagrange_interpolation(
    x_coords: *const f64,
    y_coords: *const f64,
    len: usize,
) -> *mut Polynomial { unsafe {

    if x_coords.is_null()
        || y_coords.is_null()
    {

        update_last_error("Null pointer passed to rssn_num_lagrange_interpolation".to_string());

        return ptr::null_mut();
    }

    let x_slice =  {

        std::slice::from_raw_parts(
            x_coords,
            len,
        )
    };

    let y_slice =  {

        std::slice::from_raw_parts(
            y_coords,
            len,
        )
    };

    let points: Vec<(f64, f64)> =
        x_slice
            .iter()
            .zip(y_slice.iter())
            .map(|(&x, &y)| (x, y))
            .collect();

    match interpolate::lagrange_interpolation(&points) {
        | Ok(poly) => Box::into_raw(Box::new(poly)),
        | Err(e) => {

            update_last_error(e);

            ptr::null_mut()
        },
    }
}}

/// Creates a cubic spline interpolator handle.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_cubic_spline_interpolation(
    x_coords: *const f64,
    y_coords: *const f64,
    len: usize,
) -> *mut CubicSplineHandle { unsafe {

    if x_coords.is_null()
        || y_coords.is_null()
    {

        update_last_error("Null pointer passed to rssn_num_cubic_spline_interpolation".to_string());

        return ptr::null_mut();
    }

    let x_slice =  {

        std::slice::from_raw_parts(
            x_coords,
            len,
        )
    };

    let y_slice =  {

        std::slice::from_raw_parts(
            y_coords,
            len,
        )
    };

    let points: Vec<(f64, f64)> =
        x_slice
            .iter()
            .zip(y_slice.iter())
            .map(|(&x, &y)| (x, y))
            .collect();

    match interpolate::cubic_spline_interpolation(&points) {
        | Ok(spline) => Box::into_raw(Box::new(spline)),
        | Err(e) => {

            update_last_error(e);

            ptr::null_mut()
        },
    }
}}

/// Evaluates a cubic spline at a given x coordinate.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_cubic_spline_evaluate(
    handle: *const CubicSplineHandle,
    x: f64,
) -> f64 { unsafe {

    if handle.is_null() {

        update_last_error("Null pointer passed to rssn_num_cubic_spline_evaluate".to_string());

        return f64::NAN;
    }

    let spline =  {

        &*handle
    };

    spline(x)
}}

/// Frees a cubic spline handle.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_cubic_spline_free(
    handle: *mut CubicSplineHandle
) {

    if !handle.is_null() {

        unsafe {

            let _ =
                Box::from_raw(handle);
        }
    }
}

/// Evaluates a Bézier curve at parameter t.
/// `control_points` is a flattened array of size `n_points` * dim.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_bezier_curve(
    control_points: *const f64,
    n_points: usize,
    dim: usize,
    t: f64,
    out_point: *mut f64,
) -> i32 {

    if control_points.is_null()
        || out_point.is_null()
    {

        return -1;
    }

    let data = unsafe {

        std::slice::from_raw_parts(
            control_points,
            n_points * dim,
        )
    };

    let mut cp_vecs =
        Vec::with_capacity(n_points);

    for i in 0 .. n_points {

        cp_vecs.push(
            data[i * dim
                .. (i + 1) * dim]
                .to_vec(),
        );
    }

    let result =
        interpolate::bezier_curve(
            &cp_vecs,
            t,
        );

    unsafe {

        std::ptr::copy_nonoverlapping(
            result.as_ptr(),
            out_point,
            dim,
        );
    }

    0
}

/// Evaluates a B-spline curve at parameter t.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_b_spline(
    control_points: *const f64,
    n_points: usize,
    dim: usize,
    degree: usize,
    knots: *const f64,
    n_knots: usize,
    t: f64,
    out_point: *mut f64,
) -> i32 {

    if control_points.is_null()
        || knots.is_null()
        || out_point.is_null()
    {

        return -1;
    }

    let data = unsafe {

        std::slice::from_raw_parts(
            control_points,
            n_points * dim,
        )
    };

    let mut cp_vecs =
        Vec::with_capacity(n_points);

    for i in 0 .. n_points {

        cp_vecs.push(
            data[i * dim
                .. (i + 1) * dim]
                .to_vec(),
        );
    }

    let knot_slice = unsafe {

        std::slice::from_raw_parts(
            knots,
            n_knots,
        )
    };

    match interpolate::b_spline(
        &cp_vecs,
        degree,
        knot_slice,
        t,
    ) {
        | Some(result) => {

            unsafe {

                std::ptr::copy_nonoverlapping(
                    result.as_ptr(),
                    out_point,
                    dim,
                );
            }

            0
        },
        | None => -1,
    }
}
