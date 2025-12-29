use std::ffi::c_char;

use crate::ffi_apis::common::{from_json_string, to_json_string};
use crate::symbolic::tensor::Tensor;

/// Performs tensor addition.

///

/// Takes two JSON strings representing `Tensor` objects as input,

/// and returns a JSON string representing their sum.

#[no_mangle]

pub extern "C" fn rssn_json_tensor_add(
    t1_json: *const c_char,
    t2_json: *const c_char,
) -> *mut c_char {

    let t1: Option<Tensor> =
        from_json_string(t1_json);

    let t2: Option<Tensor> =
        from_json_string(t2_json);

    if let (
        Some(tensor1),
        Some(tensor2),
    ) = (t1, t2)
    {

        match tensor1.add(&tensor2) {
            | Ok(result) => {
                to_json_string(&result)
            },
            | Err(_) => {
                std::ptr::null_mut()
            },
        }
    } else {

        std::ptr::null_mut()
    }
}

/// Performs scalar multiplication on a tensor.

///

/// Takes a JSON string representing a `Tensor` and a JSON string representing an `Expr` (scalar).

/// Returns a JSON string representing the resulting `Tensor`.

#[no_mangle]

pub extern "C" fn rssn_json_tensor_scalar_mul(
    t_json: *const c_char,
    scalar_json: *const c_char,
) -> *mut c_char {

    let t: Option<Tensor> =
        from_json_string(t_json);

    let scalar: Option<
        crate::symbolic::core::Expr,
    > = from_json_string(scalar_json);

    if let (Some(tensor), Some(s)) =
        (t, scalar)
    {

        match tensor.scalar_mul(&s) {
            | Ok(result) => {
                to_json_string(&result)
            },
            | Err(_) => {
                std::ptr::null_mut()
            },
        }
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the outer product of two tensors.

///

/// Takes two JSON strings representing `Tensor` objects as input,

/// and returns a JSON string representing their outer product.

#[no_mangle]

pub extern "C" fn rssn_json_tensor_outer_product(
    t1_json: *const c_char,
    t2_json: *const c_char,
) -> *mut c_char {

    let t1: Option<Tensor> =
        from_json_string(t1_json);

    let t2: Option<Tensor> =
        from_json_string(t2_json);

    if let (
        Some(tensor1),
        Some(tensor2),
    ) = (t1, t2)
    {

        match tensor1
            .outer_product(&tensor2)
        {
            | Ok(result) => {
                to_json_string(&result)
            },
            | Err(_) => {
                std::ptr::null_mut()
            },
        }
    } else {

        std::ptr::null_mut()
    }
}

/// Performs tensor contraction.

///

/// Takes a JSON string representing a `Tensor`, and two `usize` values representing the axes to contract.

/// Returns a JSON string representing the contracted `Tensor`.

#[no_mangle]

pub extern "C" fn rssn_json_tensor_contract(
    t_json: *const c_char,
    axis1: usize,
    axis2: usize,
) -> *mut c_char {

    let t: Option<Tensor> =
        from_json_string(t_json);

    if let Some(tensor) = t {

        match tensor
            .contract(axis1, axis2)
        {
            | Ok(result) => {
                to_json_string(&result)
            },
            | Err(_) => {
                std::ptr::null_mut()
            },
        }
    } else {

        std::ptr::null_mut()
    }
}
