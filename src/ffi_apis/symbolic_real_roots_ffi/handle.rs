use std::ffi::CStr;
use std::os::raw::c_char;
use std::os::raw::c_double;

use crate::symbolic::core::Expr;
use crate::symbolic::polynomial::expr_to_sparse_poly;
use crate::symbolic::polynomial::sparse_poly_to_expr;
use crate::symbolic::real_roots::count_real_roots_in_interval;
use crate::symbolic::real_roots::isolate_real_roots;
use crate::symbolic::real_roots::sturm_sequence;

/// Generates the Sturm sequence for a given polynomial (Handle)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_sturm_sequence_handle(
    expr_ptr: *const Expr,
    var_ptr: *const c_char,
) -> *mut Vec<Expr> {

    if expr_ptr.is_null()
        || var_ptr.is_null()
    {

        return std::ptr::null_mut();
    }

    unsafe {

        let expr = &*expr_ptr;

        let var_cstr =
            CStr::from_ptr(var_ptr);

        let var_str = match var_cstr.to_str() {
            | Ok(s) => s,
            | Err(_) => return std::ptr::null_mut(),
        };

        // Convert Expr to SparsePolynomial
        // Note: expr_to_sparse_poly requires a list of variables if we want to be precise,
        // but here we are treating it as a univariate polynomial in `var`.
        // The implementation in integration.rs used `expr_to_sparse_poly(expr, &[var])`.
        let poly = expr_to_sparse_poly(
            expr,
            &[var_str],
        );

        let seq = sturm_sequence(
            &poly,
            var_str,
        );

        // Convert back to Exprs
        let expr_seq: Vec<Expr> = seq
            .into_iter()
            .map(|p| {

                sparse_poly_to_expr(&p)
            })
            .collect();

        Box::into_raw(Box::new(
            expr_seq,
        ))
    }
}

/// Counts the number of distinct real roots in an interval (Handle)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_count_real_roots_in_interval_handle(
    expr_ptr: *const Expr,
    var_ptr: *const c_char,
    a: c_double,
    b: c_double,
) -> i64 {

    if expr_ptr.is_null()
        || var_ptr.is_null()
    {

        return -1;
    }

    unsafe {

        let expr = &*expr_ptr;

        let var_cstr =
            CStr::from_ptr(var_ptr);

        let var_str =
            match var_cstr.to_str() {
                | Ok(s) => s,
                | Err(_) => return -1,
            };

        let poly = expr_to_sparse_poly(
            expr,
            &[var_str],
        );

        match count_real_roots_in_interval(&poly, var_str, a, b) {
            | Ok(count) => count as i64,
            | Err(_) => -1,
        }
    }
}

/// Isolates real roots in an interval (Handle)
/// Returns a pointer to a Vec<(f64, f64)>
#[unsafe(no_mangle)]

pub extern "C" fn rssn_isolate_real_roots_handle(
    expr_ptr: *const Expr,
    var_ptr: *const c_char,
    precision: c_double,
) -> *mut Vec<(f64, f64)> {

    if expr_ptr.is_null()
        || var_ptr.is_null()
    {

        return std::ptr::null_mut();
    }

    unsafe {

        let expr = &*expr_ptr;

        let var_cstr =
            CStr::from_ptr(var_ptr);

        let var_str = match var_cstr.to_str() {
            | Ok(s) => s,
            | Err(_) => return std::ptr::null_mut(),
        };

        let poly = expr_to_sparse_poly(
            expr,
            &[var_str],
        );

        match isolate_real_roots(
            &poly,
            var_str,
            precision,
        ) {
            | Ok(roots) => {
                Box::into_raw(Box::new(
                    roots,
                ))
            },
            | Err(_) => {
                std::ptr::null_mut()
            },
        }
    }
}

/// Frees a `Vec<Expr>` handle
#[unsafe(no_mangle)]

pub extern "C" fn rssn_free_expr_vec_handle(
    ptr: *mut Vec<Expr>
) {

    if !ptr.is_null() {

        unsafe {

            let _ = Box::from_raw(ptr);
        }
    }
}

/// Frees a Vec<(f64, f64)> handle
#[unsafe(no_mangle)]

pub extern "C" fn rssn_free_interval_vec_handle(
    ptr: *mut Vec<(f64, f64)>
) {

    if !ptr.is_null() {

        unsafe {

            let _ = Box::from_raw(ptr);
        }
    }
}
