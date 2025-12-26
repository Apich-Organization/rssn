//! Handle-based FFI API for numerical calculus operations.

use std::ffi::CStr;
use std::os::raw::c_char;
use std::ptr;

use crate::ffi_apis::ffi_api::update_last_error;
use crate::numerical::calculus;
use crate::numerical::matrix::Matrix;
use crate::symbolic::core::Expr;

/// Computes the numerical partial derivative of a function with respect to a variable at a point.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_calculus_partial_derivative(
    f: *const Expr,
    var: *const c_char,
    x: f64,
    result: *mut f64,
) -> i32 {

    if f.is_null()
        || var.is_null()
        || result.is_null()
    {

        return -1;
    }

    let f_expr = &*f;

    let var_str =
        match CStr::from_ptr(var)
            .to_str()
        {
            | Ok(s) => s,
            | Err(_) => {

                update_last_error(
                    "Invalid UTF-8 \
                     string for \
                     variable name"
                        .to_string(),
                );

                return -1;
            },
        };

    match calculus::partial_derivative(
        f_expr,
        var_str,
        x,
    ) {
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

/// Computes the numerical gradient of a function at a point.
/// Returns a pointer to a Vec<f64> containing the gradient.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_calculus_gradient(
    f: *const Expr,
    vars: *const *const c_char,
    point: *const f64,
    n_vars: usize,
) -> *mut Vec<f64> {

    if f.is_null()
        || vars.is_null()
        || point.is_null()
    {

        return ptr::null_mut();
    }

    let f_expr = &*f;

    let mut vars_list =
        Vec::with_capacity(n_vars);

    for i in 0..n_vars {

        let v_ptr = *vars.add(i);

        if v_ptr.is_null() {

            update_last_error(format!(
                "Variable at index {} \
                 is null",
                i
            ));

            return ptr::null_mut();
        }

        match CStr::from_ptr(v_ptr)
            .to_str()
        {
            | Ok(s) => {
                vars_list.push(s)
            },
            | Err(_) => {

                update_last_error(format!(
                    "Invalid UTF-8 for variable at index \
                     {}",
                    i
                ));

                return ptr::null_mut();
            },
        }
    }

    let point_slice =
        std::slice::from_raw_parts(
            point,
            n_vars,
        );

    match calculus::gradient(
        f_expr,
        &vars_list,
        point_slice,
    ) {
        | Ok(grad) => {
            Box::into_raw(Box::new(
                grad,
            ))
        },
        | Err(e) => {

            update_last_error(e);

            ptr::null_mut()
        },
    }
}

/// Computes the numerical Jacobian matrix of a vector-valued function at a point.
/// Returns a pointer to a Matrix<f64>.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_calculus_jacobian(
    funcs: *const *const Expr,
    n_funcs: usize,
    vars: *const *const c_char,
    point: *const f64,
    n_vars: usize,
) -> *mut Matrix<f64> {

    if funcs.is_null()
        || vars.is_null()
        || point.is_null()
    {

        return ptr::null_mut();
    }

    let mut funcs_list =
        Vec::with_capacity(n_funcs);

    for i in 0..n_funcs {

        funcs_list.push(
            (*(*funcs.add(i))).clone(),
        );
    }

    let mut vars_list =
        Vec::with_capacity(n_vars);

    for i in 0..n_vars {

        let v_ptr = *vars.add(i);

        match CStr::from_ptr(v_ptr)
            .to_str()
        {
            | Ok(s) => {
                vars_list.push(s)
            },
            | Err(_) => {

                update_last_error(format!(
                    "Invalid UTF-8 for variable at index \
                     {}",
                    i
                ));

                return ptr::null_mut();
            },
        }
    }

    let point_slice =
        std::slice::from_raw_parts(
            point,
            n_vars,
        );

    match calculus::jacobian(
        &funcs_list,
        &vars_list,
        point_slice,
    ) {
        | Ok(jac_vecs) => {

            let rows = jac_vecs.len();

            let cols = if rows > 0 {

                jac_vecs[0].len()
            } else {

                0
            };

            let mut flattened =
                Vec::with_capacity(
                    rows * cols,
                );

            for row in jac_vecs {

                flattened.extend(row);
            }

            Box::into_raw(Box::new(
                Matrix::new(
                    rows,
                    cols,
                    flattened,
                ),
            ))
        },
        | Err(e) => {

            update_last_error(e);

            ptr::null_mut()
        },
    }
}

/// Computes the numerical Hessian matrix of a scalar function at a point.
/// Returns a pointer to a Matrix<f64>.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_calculus_hessian(
    f: *const Expr,
    vars: *const *const c_char,
    point: *const f64,
    n_vars: usize,
) -> *mut Matrix<f64> {

    if f.is_null()
        || vars.is_null()
        || point.is_null()
    {

        return ptr::null_mut();
    }

    let f_expr = &*f;

    let mut vars_list =
        Vec::with_capacity(n_vars);

    for i in 0..n_vars {

        let v_ptr = *vars.add(i);

        match CStr::from_ptr(v_ptr)
            .to_str()
        {
            | Ok(s) => {
                vars_list.push(s)
            },
            | Err(_) => {

                update_last_error(format!(
                    "Invalid UTF-8 for variable at index \
                     {}",
                    i
                ));

                return ptr::null_mut();
            },
        }
    }

    let point_slice =
        std::slice::from_raw_parts(
            point,
            n_vars,
        );

    match calculus::hessian(
        f_expr,
        &vars_list,
        point_slice,
    ) {
        | Ok(hess_vecs) => {

            let rows = hess_vecs.len();

            let cols = if rows > 0 {

                hess_vecs[0].len()
            } else {

                0
            };

            let mut flattened =
                Vec::with_capacity(
                    rows * cols,
                );

            for row in hess_vecs {

                flattened.extend(row);
            }

            Box::into_raw(Box::new(
                Matrix::new(
                    rows,
                    cols,
                    flattened,
                ),
            ))
        },
        | Err(e) => {

            update_last_error(e);

            ptr::null_mut()
        },
    }
}
