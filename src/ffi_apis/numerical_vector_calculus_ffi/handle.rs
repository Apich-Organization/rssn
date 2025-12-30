//! Handle-based FFI API for numerical vector calculus.

use std::ffi::CStr;
use std::os::raw::c_char;
use std::ptr;

use crate::ffi_apis::ffi_api::update_last_error;
use crate::numerical::vector_calculus;
use crate::symbolic::core::Expr;

/// Computes the numerical divergence of a vector field at a point.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_vector_calculus_divergence(
    funcs: *const *const Expr,
    n_funcs: usize,
    vars: *const *const c_char,
    point: *const f64,
    n_vars: usize,
    result: *mut f64,
) -> i32 { unsafe {

    if funcs.is_null()
        || vars.is_null()
        || point.is_null()
        || result.is_null()
    {

        return -1;
    }

    let mut funcs_list =
        Vec::with_capacity(n_funcs);

    for i in 0 .. n_funcs {

        funcs_list.push(
            (*(*funcs.add(i))).clone(),
        );
    }

    let mut vars_list =
        Vec::with_capacity(n_vars);

    for i in 0 .. n_vars {

        let v_ptr = *vars.add(i);

        match CStr::from_ptr(v_ptr)
            .to_str()
        {
            | Ok(s) => {

                vars_list.push(s);
            },
            | Err(_) => {

                update_last_error(format!(
                    "Invalid UTF-8 for variable at index {i}"
                ));

                return -1;
            },
        }
    }

    let point_slice =
        std::slice::from_raw_parts(
            point,
            n_vars,
        );

    match vector_calculus::divergence_expr(
        &funcs_list,
        &vars_list,
        point_slice,
    ) {
        | Ok(v) => {

            *result = v;

            0
        },
        | Err(e) => {

            update_last_error(e);

            -1
        },
    }
}}

/// Computes the numerical curl of a 3D vector field at a point.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_vector_calculus_curl(
    funcs: *const *const Expr,
    vars: *const *const c_char,
    point: *const f64,
) -> *mut Vec<f64> { unsafe {

    if funcs.is_null()
        || vars.is_null()
        || point.is_null()
    {

        return ptr::null_mut();
    }

    let mut funcs_list =
        Vec::with_capacity(3);

    for i in 0 .. 3 {

        funcs_list.push(
            (*(*funcs.add(i))).clone(),
        );
    }

    let mut vars_list =
        Vec::with_capacity(3);

    for i in 0 .. 3 {

        let v_ptr = *vars.add(i);

        match CStr::from_ptr(v_ptr)
            .to_str()
        {
            | Ok(s) => {

                vars_list.push(s);
            },
            | Err(_) => {

                update_last_error(format!(
                    "Invalid UTF-8 for variable at index {i}"
                ));

                return ptr::null_mut();
            },
        }
    }

    let point_slice =
        std::slice::from_raw_parts(
            point, 3,
        );

    match vector_calculus::curl_expr(
        &funcs_list,
        &vars_list,
        point_slice,
    ) {
        | Ok(v) => {
            Box::into_raw(Box::new(v))
        },
        | Err(e) => {

            update_last_error(e);

            ptr::null_mut()
        },
    }
}}

/// Computes the numerical Laplacian of a scalar field at a point.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_vector_calculus_laplacian(
    f: *const Expr,
    vars: *const *const c_char,
    point: *const f64,
    n_vars: usize,
    result: *mut f64,
) -> i32 { unsafe {

    if f.is_null()
        || vars.is_null()
        || point.is_null()
        || result.is_null()
    {

        return -1;
    }

    let mut vars_list =
        Vec::with_capacity(n_vars);

    for i in 0 .. n_vars {

        let v_ptr = *vars.add(i);

        match CStr::from_ptr(v_ptr)
            .to_str()
        {
            | Ok(s) => {

                vars_list.push(s);
            },
            | Err(_) => {

                update_last_error(format!(
                    "Invalid UTF-8 for variable at index {i}"
                ));

                return -1;
            },
        }
    }

    let point_slice =
        std::slice::from_raw_parts(
            point,
            n_vars,
        );

    match vector_calculus::laplacian(
        &*f,
        &vars_list,
        point_slice,
    ) {
        | Ok(v) => {

            *result = v;

            0
        },
        | Err(e) => {

            update_last_error(e);

            -1
        },
    }
}}

/// Computes the numerical directional derivative of a function at a point.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_vector_calculus_directional_derivative(
    f: *const Expr,
    vars: *const *const c_char,
    point: *const f64,
    direction: *const f64,
    n_vars: usize,
    result: *mut f64,
) -> i32 { unsafe {

    if f.is_null()
        || vars.is_null()
        || point.is_null()
        || direction.is_null()
        || result.is_null()
    {

        return -1;
    }

    let mut vars_list =
        Vec::with_capacity(n_vars);

    for i in 0 .. n_vars {

        let v_ptr = *vars.add(i);

        match CStr::from_ptr(v_ptr)
            .to_str()
        {
            | Ok(s) => {

                vars_list.push(s);
            },
            | Err(_) => {

                update_last_error(format!(
                    "Invalid UTF-8 for variable at index {i}"
                ));

                return -1;
            },
        }
    }

    let point_slice =
        std::slice::from_raw_parts(
            point,
            n_vars,
        );

    let direction_slice =
        std::slice::from_raw_parts(
            direction,
            n_vars,
        );

    match vector_calculus::directional_derivative(
        &*f,
        &vars_list,
        point_slice,
        direction_slice,
    ) {
        | Ok(v) => {

            *result = v;

            0
        },
        | Err(e) => {

            update_last_error(e);

            -1
        },
    }
}}
