use std::os::raw::c_char;
use std::os::raw::c_int;
use std::sync::Arc;

use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::stats_regression;

unsafe fn collect_pairs(
    x_data: *const *const Expr,
    y_data: *const *const Expr,
    len: usize,
) -> Vec<(Expr, Expr)> {

    let mut data =
        Vec::with_capacity(len);

    for i in 0 .. len {

        let x_ptr = *x_data.add(i);

        let y_ptr = *y_data.add(i);

        if !x_ptr.is_null()
            && !y_ptr.is_null()
        {

            data.push((
                (*x_ptr).clone(),
                (*y_ptr).clone(),
            ));
        }
    }

    data
}

/// Performs a simple linear regression.

///

/// Takes raw pointers to arrays of `Expr` (x and y data) and the length of the data.

/// Returns a raw pointer to an `Expr` (vector) containing the intercept and slope coefficients.

#[no_mangle]

pub unsafe extern "C" fn rssn_simple_linear_regression(
    x_data: *const *const Expr,
    y_data: *const *const Expr,
    len: usize,
) -> *mut Expr {

    if x_data.is_null()
        || y_data.is_null()
    {

        return std::ptr::null_mut();
    }

    let data = collect_pairs(
        x_data,
        y_data,
        len,
    );

    let (b0, b1) = stats_regression::simple_linear_regression_symbolic(&data);

    Box::into_raw(Box::new(
        Expr::Vector(vec![b0, b1]),
    ))
}

/// Performs a polynomial regression.

///

/// Takes raw pointers to arrays of `Expr` (x and y data), the length of the data,

/// and the degree of the polynomial.

/// Returns a raw pointer to an `Expr` (vector) containing the coefficients of the polynomial.

#[no_mangle]

pub unsafe extern "C" fn rssn_polynomial_regression(
    x_data: *const *const Expr,
    y_data: *const *const Expr,
    len: usize,
    degree: usize,
) -> *mut Expr {

    if x_data.is_null()
        || y_data.is_null()
    {

        return std::ptr::null_mut();
    }

    let data = collect_pairs(
        x_data,
        y_data,
        len,
    );

    match stats_regression::polynomial_regression_symbolic(&data, degree) {
        | Ok(coeffs) => {
            Box::into_raw(Box::new(
                Expr::Vector(coeffs),
            ))
        },
        | Err(_) => std::ptr::null_mut(),
    }
}

/// Performs a nonlinear regression.

///

/// Takes raw pointers to arrays of `Expr` (x and y data), the length of the data,

/// a raw pointer to an `Expr` (model), raw pointers to arrays of C-style strings (variables and parameters),

/// and their respective lengths.

/// Returns a raw pointer to an `Expr` representing the solutions (optimized parameter values).

#[no_mangle]

pub unsafe extern "C" fn rssn_nonlinear_regression(
    x_data: *const *const Expr,
    y_data: *const *const Expr,
    len: usize,
    model: *const Expr,
    vars: *const *const c_char,
    vars_len: usize,
    params: *const *const c_char,
    params_len: usize,
) -> *mut Expr {

    if x_data.is_null()
        || y_data.is_null()
        || model.is_null()
    {

        return std::ptr::null_mut();
    }

    let data = collect_pairs(
        x_data,
        y_data,
        len,
    );

    let model_expr = &*model;

    // Collect strings
    let mut vars_vec =
        Vec::with_capacity(vars_len);

    for i in 0 .. vars_len {

        if let Some(s) =
            c_str_to_str(*vars.add(i))
        {

            vars_vec.push(s);
        }
    }

    let mut params_vec =
        Vec::with_capacity(params_len);

    for i in 0 .. params_len {

        if let Some(s) =
            c_str_to_str(*params.add(i))
        {

            params_vec.push(s);
        }
    }

    match stats_regression::nonlinear_regression_symbolic(
        &data,
        model_expr,
        &vars_vec,
        &params_vec,
    ) {
        | Some(solutions) => {

            let eqs = solutions
                .into_iter()
                .map(|(p, v)| {

                    Expr::Eq(
                        Arc::new(p),
                        Arc::new(v),
                    )
                })
                .collect();

            Box::into_raw(Box::new(
                Expr::Solutions(eqs),
            ))
        },
        | None => std::ptr::null_mut(),
    }
}
