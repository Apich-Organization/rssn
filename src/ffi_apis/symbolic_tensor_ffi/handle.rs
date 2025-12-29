use crate::symbolic::tensor::Tensor;

/// Performs tensor addition.

///

/// Takes raw pointers to two `Tensor` objects as input,

/// and returns a raw pointer to a new `Tensor` representing their sum.

#[no_mangle]

pub extern "C" fn rssn_tensor_add_handle(
    t1: *const Tensor,
    t2: *const Tensor,
) -> *mut Tensor {

    let t1_ref = unsafe {

        &*t1
    };

    let t2_ref = unsafe {

        &*t2
    };

    match t1_ref.add(t2_ref) {
        | Ok(result) => {
            Box::into_raw(Box::new(
                result,
            ))
        },
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
}

/// Performs scalar multiplication on a tensor.

///

/// Takes a raw pointer to a `Tensor` and a raw pointer to an `Expr` (scalar).

/// Returns a raw pointer to a new `Tensor` representing the result.

#[no_mangle]

pub extern "C" fn rssn_tensor_scalar_mul_handle(
    t: *const Tensor,
    scalar : *const crate::symbolic::core::Expr,
) -> *mut Tensor {

    let t_ref = unsafe {

        &*t
    };

    let scalar_ref = unsafe {

        &*scalar
    };

    match t_ref.scalar_mul(scalar_ref) {
        | Ok(result) => {
            Box::into_raw(Box::new(
                result,
            ))
        },
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
}

/// Computes the outer product of two tensors.

///

/// Takes raw pointers to two `Tensor` objects as input,

/// and returns a raw pointer to a new `Tensor` representing their outer product.

#[no_mangle]

pub extern "C" fn rssn_tensor_outer_product_handle(
    t1: *const Tensor,
    t2: *const Tensor,
) -> *mut Tensor {

    let t1_ref = unsafe {

        &*t1
    };

    let t2_ref = unsafe {

        &*t2
    };

    match t1_ref.outer_product(t2_ref) {
        | Ok(result) => {
            Box::into_raw(Box::new(
                result,
            ))
        },
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
}

/// Performs tensor contraction.

///

/// Takes a raw pointer to a `Tensor`, and two `usize` values representing the axes to contract.

/// Returns a raw pointer to a new `Tensor` representing the contracted tensor.

#[no_mangle]

pub extern "C" fn rssn_tensor_contract_handle(
    t: *const Tensor,
    axis1: usize,
    axis2: usize,
) -> *mut Tensor {

    let t_ref = unsafe {

        &*t
    };

    match t_ref.contract(axis1, axis2) {
        | Ok(result) => {
            Box::into_raw(Box::new(
                result,
            ))
        },
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
}
