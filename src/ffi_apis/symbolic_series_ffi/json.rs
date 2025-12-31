use std::ffi::c_char;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_json_string;
use crate::symbolic::core::Expr;
use crate::symbolic::series::analytic_continuation;
use crate::symbolic::series::analyze_convergence;
use crate::symbolic::series::asymptotic_expansion;
use crate::symbolic::series::fourier_series;
use crate::symbolic::series::laurent_series;
use crate::symbolic::series::product;
use crate::symbolic::series::summation;
use crate::symbolic::series::taylor_series;

/// Computes the Taylor series expansion of an expression.

///

/// Takes JSON strings representing `Expr` (expression), `String` (variable),

/// `Expr` (center), and `usize` (order).

/// Returns a JSON string representing the `Expr` of the Taylor series.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_taylor_series(
    expr_json: *const c_char,
    var_json: *const c_char,
    center_json: *const c_char,
    order_json: *const c_char,
) -> *mut c_char {

    let expr: Option<Expr> =
        from_json_string(expr_json);

    let var: Option<String> =
        from_json_string(var_json);

    let center: Option<Expr> =
        from_json_string(center_json);

    let order: Option<usize> =
        from_json_string(order_json);

    match (
        expr,
        var,
        center,
        order,
    ) {
        | (
            Some(e),
            Some(v),
            Some(c),
            Some(o),
        ) => {

            let result = taylor_series(
                &e, &v, &c, o,
            );

            to_json_string(&result)
        },
        | _ => std::ptr::null_mut(),
    }
}

/// Computes the Laurent series expansion of an expression.

///

/// Takes JSON strings representing `Expr` (expression), `String` (variable),

/// `Expr` (center), and `usize` (order).

/// Returns a JSON string representing the `Expr` of the Laurent series.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_laurent_series(
    expr_json: *const c_char,
    var_json: *const c_char,
    center_json: *const c_char,
    order_json: *const c_char,
) -> *mut c_char {

    let expr: Option<Expr> =
        from_json_string(expr_json);

    let var: Option<String> =
        from_json_string(var_json);

    let center: Option<Expr> =
        from_json_string(center_json);

    let order: Option<usize> =
        from_json_string(order_json);

    match (
        expr,
        var,
        center,
        order,
    ) {
        | (
            Some(e),
            Some(v),
            Some(c),
            Some(o),
        ) => {

            let result = laurent_series(
                &e, &v, &c, o,
            );

            to_json_string(&result)
        },
        | _ => std::ptr::null_mut(),
    }
}

/// Computes the Fourier series expansion of an expression.

///

/// Takes JSON strings representing `Expr` (expression), `String` (variable),

/// `Expr` (period), and `usize` (order).

/// Returns a JSON string representing the `Expr` of the Fourier series.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_fourier_series(
    expr_json: *const c_char,
    var_json: *const c_char,
    period_json: *const c_char,
    order_json: *const c_char,
) -> *mut c_char {

    let expr: Option<Expr> =
        from_json_string(expr_json);

    let var: Option<String> =
        from_json_string(var_json);

    let period: Option<Expr> =
        from_json_string(period_json);

    let order: Option<usize> =
        from_json_string(order_json);

    match (
        expr,
        var,
        period,
        order,
    ) {
        | (
            Some(e),
            Some(v),
            Some(p),
            Some(o),
        ) => {

            let result = fourier_series(
                &e, &v, &p, o,
            );

            to_json_string(&result)
        },
        | _ => std::ptr::null_mut(),
    }
}

/// Computes the summation of an expression.

///

/// Takes JSON strings representing `Expr` (expression), `String` (variable),

/// `Expr` (lower bound), and `Expr` (upper bound).

/// Returns a JSON string representing the `Expr` of the summation.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_summation(
    expr_json: *const c_char,
    var_json: *const c_char,
    lower_json: *const c_char,
    upper_json: *const c_char,
) -> *mut c_char {

    let expr: Option<Expr> =
        from_json_string(expr_json);

    let var: Option<String> =
        from_json_string(var_json);

    let lower: Option<Expr> =
        from_json_string(lower_json);

    let upper: Option<Expr> =
        from_json_string(upper_json);

    match (
        expr, var, lower, upper,
    ) {
        | (
            Some(e),
            Some(v),
            Some(l),
            Some(u),
        ) => {

            let result = summation(
                &e, &v, &l, &u,
            );

            to_json_string(&result)
        },
        | _ => std::ptr::null_mut(),
    }
}

/// Computes the product of an expression.

///

/// Takes JSON strings representing `Expr` (expression), `String` (variable),

/// `Expr` (lower bound), and `Expr` (upper bound).

/// Returns a JSON string representing the `Expr` of the product.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_product(
    expr_json: *const c_char,
    var_json: *const c_char,
    lower_json: *const c_char,
    upper_json: *const c_char,
) -> *mut c_char {

    let expr: Option<Expr> =
        from_json_string(expr_json);

    let var: Option<String> =
        from_json_string(var_json);

    let lower: Option<Expr> =
        from_json_string(lower_json);

    let upper: Option<Expr> =
        from_json_string(upper_json);

    match (
        expr, var, lower, upper,
    ) {
        | (
            Some(e),
            Some(v),
            Some(l),
            Some(u),
        ) => {

            let result =
                product(&e, &v, &l, &u);

            to_json_string(&result)
        },
        | _ => std::ptr::null_mut(),
    }
}

/// Analyzes the convergence of a series.

///

/// Takes JSON strings representing `Expr` (series expression) and `String` (variable).

/// Returns a JSON string representing the `Expr` of the convergence analysis result.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_series_json_analyze_convergence(
    series_json: *const c_char,
    var_json: *const c_char,
) -> *mut c_char {

    let series: Option<Expr> =
        from_json_string(series_json);

    let var: Option<String> =
        from_json_string(var_json);

    match (series, var) {
        | (Some(s), Some(v)) => {

            let result =
                analyze_convergence(
                    &s, &v,
                );

            to_json_string(&result)
        },
        | _ => std::ptr::null_mut(),
    }
}

/// Computes the asymptotic expansion of an expression.

///

/// Takes JSON strings representing `Expr` (expression), `String` (variable),

/// `Expr` (point), and `usize` (order).

/// Returns a JSON string representing the `Expr` of the asymptotic expansion.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_asymptotic_expansion(
    expr_json: *const c_char,
    var_json: *const c_char,
    point_json: *const c_char,
    order_json: *const c_char,
) -> *mut c_char {

    let expr: Option<Expr> =
        from_json_string(expr_json);

    let var: Option<String> =
        from_json_string(var_json);

    let point: Option<Expr> =
        from_json_string(point_json);

    let order: Option<usize> =
        from_json_string(order_json);

    match (
        expr, var, point, order,
    ) {
        | (
            Some(e),
            Some(v),
            Some(p),
            Some(o),
        ) => {

            let result =
                asymptotic_expansion(
                    &e, &v, &p, o,
                );

            to_json_string(&result)
        },
        | _ => std::ptr::null_mut(),
    }
}

/// Computes the analytic continuation of a series.

///

/// Takes JSON strings representing `Expr` (expression), `String` (variable),

/// `Expr` (original center), `Expr` (new center), and `usize` (order).

/// Returns a JSON string representing the `Expr` of the analytic continuation.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_analytic_continuation(
    expr_json: *const c_char,
    var_json: *const c_char,
    orig_center_json: *const c_char,
    new_center_json: *const c_char,
    order_json: *const c_char,
) -> *mut c_char {

    let expr: Option<Expr> =
        from_json_string(expr_json);

    let var: Option<String> =
        from_json_string(var_json);

    let orig_center: Option<Expr> =
        from_json_string(
            orig_center_json,
        );

    let new_center: Option<Expr> =
        from_json_string(
            new_center_json,
        );

    let order: Option<usize> =
        from_json_string(order_json);

    match (
        expr,
        var,
        orig_center,
        new_center,
        order,
    ) {
        | (
            Some(e),
            Some(v),
            Some(oc),
            Some(nc),
            Some(o),
        ) => {

            let result =
                analytic_continuation(
                    &e, &v, &oc, &nc, o,
                );

            to_json_string(&result)
        },
        | _ => std::ptr::null_mut(),
    }
}
