use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::optimize::*;
use std::collections::HashMap;
use std::ffi::CStr;
use std::os::raw::c_char;

/// Finds extrema of a function (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_find_extrema(
    expr_buf: BincodeBuffer,
    vars_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr: Option<Expr> = from_bincode_buffer(&expr_buf);

    let vars: Option<Vec<String>> = from_bincode_buffer(&vars_buf);

    if let (Some(e), Some(v)) = (expr, vars) {

        let vars_refs: Vec<&str> = v
            .iter()
            .map(|s| s.as_str())
            .collect();

        match find_extrema(&e, &vars_refs) {
            Ok(points) => to_bincode_buffer(&points),
            Err(_) => BincodeBuffer::empty(),
        }
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes Hessian matrix (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_hessian_matrix(
    expr_buf: BincodeBuffer,
    vars_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr: Option<Expr> = from_bincode_buffer(&expr_buf);

    let vars: Option<Vec<String>> = from_bincode_buffer(&vars_buf);

    if let (Some(e), Some(v)) = (expr, vars) {

        let vars_refs: Vec<&str> = v
            .iter()
            .map(|s| s.as_str())
            .collect();

        let hessian = hessian_matrix(&e, &vars_refs);

        to_bincode_buffer(&hessian)
    } else {

        BincodeBuffer::empty()
    }
}

/// Finds constrained extrema (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_find_constrained_extrema(
    expr_buf: BincodeBuffer,
    constraints_buf: BincodeBuffer,
    vars_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr: Option<Expr> = from_bincode_buffer(&expr_buf);

    let constraints: Option<Vec<Expr>> = from_bincode_buffer(&constraints_buf);

    let vars: Option<Vec<String>> = from_bincode_buffer(&vars_buf);

    if let (Some(e), Some(c), Some(v)) = (expr, constraints, vars) {

        let vars_refs: Vec<&str> = v
            .iter()
            .map(|s| s.as_str())
            .collect();

        match find_constrained_extrema(&e, &c, &vars_refs) {
            Ok(solutions) => to_bincode_buffer(&solutions),
            Err(_) => BincodeBuffer::empty(),
        }
    } else {

        BincodeBuffer::empty()
    }
}
