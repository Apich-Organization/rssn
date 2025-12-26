use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::polynomial::{
    expr_to_sparse_poly,
    sparse_poly_to_expr,
};
use crate::symbolic::real_roots::*;
use std::ffi::CStr;
use std::os::raw::c_char;

/// Generates the Sturm sequence for a given polynomial (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_sturm_sequence(
    expr_buf: BincodeBuffer,
    var_ptr: *const c_char,
) -> BincodeBuffer {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    if let Some(e) = expr {

        unsafe {

            if var_ptr.is_null() {

                return BincodeBuffer::empty();
            }

            let var_cstr =
                CStr::from_ptr(var_ptr);

            let var_str = match var_cstr.to_str() {
                | Ok(s) => s,
                | Err(_) => return BincodeBuffer::empty(),
            };

            let poly =
                expr_to_sparse_poly(
                    &e,
                    &[var_str],
                );

            let seq = sturm_sequence(
                &poly, var_str,
            );

            let expr_seq: Vec<Expr> = seq
                .into_iter()
                .map(|p| sparse_poly_to_expr(&p))
                .collect();

            to_bincode_buffer(&expr_seq)
        }
    } else {

        BincodeBuffer::empty()
    }
}

/// Counts the number of distinct real roots in an interval (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_count_real_roots_in_interval(
    expr_buf: BincodeBuffer,
    var_ptr: *const c_char,
    a: f64,
    b: f64,
) -> i64 {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    if let Some(e) = expr {

        unsafe {

            if var_ptr.is_null() {

                return -1;
            }

            let var_cstr =
                CStr::from_ptr(var_ptr);

            let var_str = match var_cstr
                .to_str()
            {
                | Ok(s) => s,
                | Err(_) => return -1,
            };

            let poly =
                expr_to_sparse_poly(
                    &e,
                    &[var_str],
                );

            match count_real_roots_in_interval(
                &poly, var_str, a, b,
            ) {
                | Ok(count) => count as i64,
                | Err(_) => -1,
            }
        }
    } else {

        -1
    }
}

/// Isolates real roots in an interval (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_isolate_real_roots(
    expr_buf: BincodeBuffer,
    var_ptr: *const c_char,
    precision: f64,
) -> BincodeBuffer {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    if let Some(e) = expr {

        unsafe {

            if var_ptr.is_null() {

                return BincodeBuffer::empty();
            }

            let var_cstr =
                CStr::from_ptr(var_ptr);

            let var_str = match var_cstr.to_str() {
                | Ok(s) => s,
                | Err(_) => return BincodeBuffer::empty(),
            };

            let poly =
                expr_to_sparse_poly(
                    &e,
                    &[var_str],
                );

            match isolate_real_roots(
                &poly, var_str,
                precision,
            ) {
                | Ok(roots) => {
                    to_bincode_buffer(
                        &roots,
                    )
                },
                | Err(_) => {
                    BincodeBuffer::empty(
                    )
                },
            }
        }
    } else {

        BincodeBuffer::empty()
    }
}
