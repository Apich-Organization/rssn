//! Handle-based FFI API for numerical complex analysis.

use std::collections::HashMap;
use std::ffi::CStr;
use std::os::raw::c_char;

use num_complex::Complex;

use crate::ffi_apis::ffi_api::update_last_error;
use crate::numerical::complex_analysis;
use crate::symbolic::core::Expr;

/// Evaluates a symbolic expression to a complex number.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_complex_eval(
    expr_ptr: *const Expr,
    var_names: *const *const c_char,
    var_re: *const f64,
    var_im: *const f64,
    n_vars: usize,
    res_re: *mut f64,
    res_im: *mut f64,
) -> i32 { unsafe {

    if expr_ptr.is_null()
        || res_re.is_null()
        || res_im.is_null()
    {

        update_last_error(
            "Null pointer passed to \
             rssn_num_complex_eval"
                .to_string(),
        );

        return -1;
    }

    let expr = &*expr_ptr;

    let mut vars = HashMap::new();

    for i in 0 .. n_vars {

        let name = CStr::from_ptr(
            *var_names.add(i),
        )
        .to_string_lossy()
        .into_owned();

        let val = Complex::new(
            *var_re.add(i),
            *var_im.add(i),
        );

        vars.insert(name, val);
    }

    match complex_analysis::eval_complex_expr(expr, &vars) {
        | Ok(res) => {

            *res_re = res.re;

            *res_im = res.im;

            0
        },
        | Err(e) => {

            update_last_error(e);

            -1
        },
    }
}}

/// Computes a contour integral of a symbolic expression.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_complex_contour_integral(
    expr_ptr: *const Expr,
    var_ptr: *const c_char,
    path_re: *const f64,
    path_im: *const f64,
    path_len: usize,
    res_re: *mut f64,
    res_im: *mut f64,
) -> i32 { unsafe {

    if expr_ptr.is_null()
        || var_ptr.is_null()
        || path_re.is_null()
        || path_im.is_null()
        || res_re.is_null()
        || res_im.is_null()
    {

        update_last_error("Null pointer passed to rssn_num_complex_contour_integral".to_string());

        return -1;
    }

    let expr = &*expr_ptr;

    let var = CStr::from_ptr(var_ptr)
        .to_string_lossy();

    let path: Vec<Complex<f64>> = (0
        .. path_len)
        .map(|i| {

            Complex::new(
                *path_re.add(i),
                *path_im.add(i),
            )
        })
        .collect();

    match complex_analysis::contour_integral_expr(expr, &var, &path) {
        | Ok(res) => {

            *res_re = res.re;

            *res_im = res.im;

            0
        },
        | Err(e) => {

            update_last_error(e);

            -1
        },
    }
}}

/// Computes the residue of a symbolic expression.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_complex_residue(
    expr_ptr: *const Expr,
    var_ptr: *const c_char,
    z0_re: f64,
    z0_im: f64,
    radius: f64,
    n_points: usize,
    res_re: *mut f64,
    res_im: *mut f64,
) -> i32 { unsafe {

    if expr_ptr.is_null()
        || var_ptr.is_null()
        || res_re.is_null()
        || res_im.is_null()
    {

        update_last_error(
            "Null pointer passed to \
             rssn_num_complex_residue"
                .to_string(),
        );

        return -1;
    }

    let expr = &*expr_ptr;

    let var = CStr::from_ptr(var_ptr)
        .to_string_lossy();

    let z0 = Complex::new(z0_re, z0_im);

    match complex_analysis::residue_expr(
        expr,
        &var,
        z0,
        radius,
        n_points,
    ) {
        | Ok(res) => {

            *res_re = res.re;

            *res_im = res.im;

            0
        },
        | Err(e) => {

            update_last_error(e);

            -1
        },
    }
}}
