//! Handle-based FFI API for numerical tensor operations.

use crate::ffi_apis::ffi_api::update_last_error;
use crate::numerical::tensor::{
    self,
    TensorData,
};
use ndarray::{
    ArrayD,
    IxDyn,
};
use std::ptr;

/// Creates a new tensor from shape and data.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_tensor_create(
    shape: *const usize,
    ndim: usize,
    data: *const f64,
    data_len: usize,
) -> *mut ArrayD<f64> {

    if shape.is_null() || data.is_null() {

        update_last_error("Null pointer passed to rssn_num_tensor_create".to_string());

        return ptr::null_mut();
    }

    let s = unsafe {

        std::slice::from_raw_parts(shape, ndim)
    };

    let d = unsafe {

        std::slice::from_raw_parts(data, data_len)
    };

    match ArrayD::from_shape_vec(IxDyn(s), d.to_vec()) {
        Ok(arr) => Box::into_raw(Box::new(arr)),
        Err(e) => {

            update_last_error(e.to_string());

            ptr::null_mut()
        }
    }
}

/// Frees a tensor object.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_tensor_free(tensor: *mut ArrayD<f64>) {

    if !tensor.is_null() {

        unsafe {

            let _ = Box::from_raw(tensor);
        }
    }
}

/// Returns the number of dimensions.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_tensor_get_ndim(tensor: *const ArrayD<f64>) -> usize {

    if tensor.is_null() {

        return 0;
    }

    unsafe {

        (*tensor).ndim()
    }
}

/// Returns the shape of the tensor.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_tensor_get_shape(
    tensor: *const ArrayD<f64>,
    out_shape: *mut usize,
) -> i32 {

    if tensor.is_null() || out_shape.is_null() {

        return -1;
    }

    let t = unsafe {

        &*tensor
    };

    let shape = t.shape();

    unsafe {

        ptr::copy_nonoverlapping(
            shape.as_ptr(),
            out_shape,
            shape.len(),
        );
    }

    0
}

/// Tensor contraction (tensordot).
#[no_mangle]

pub unsafe extern "C" fn rssn_num_tensor_tensordot(
    a: *const ArrayD<f64>,
    b: *const ArrayD<f64>,
    axes_a: *const usize,
    axes_a_len: usize,
    axes_b: *const usize,
    axes_b_len: usize,
) -> *mut ArrayD<f64> {

    if a.is_null() || b.is_null() || axes_a.is_null() || axes_b.is_null() {

        return ptr::null_mut();
    }

    let ta = unsafe {

        &*a
    };

    let tb = unsafe {

        &*b
    };

    let aa = unsafe {

        std::slice::from_raw_parts(axes_a, axes_a_len)
    };

    let ab = unsafe {

        std::slice::from_raw_parts(axes_b, axes_b_len)
    };

    match tensor::tensordot(ta, tb, aa, ab) {
        Ok(res) => Box::into_raw(Box::new(res)),
        Err(e) => {

            update_last_error(e);

            ptr::null_mut()
        }
    }
}

/// Outer product of two tensors.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_tensor_outer_product(
    a: *const ArrayD<f64>,
    b: *const ArrayD<f64>,
) -> *mut ArrayD<f64> {

    if a.is_null() || b.is_null() {

        return ptr::null_mut();
    }

    let ta = unsafe {

        &*a
    };

    let tb = unsafe {

        &*b
    };

    match tensor::outer_product(ta, tb) {
        Ok(res) => Box::into_raw(Box::new(res)),
        Err(e) => {

            update_last_error(e);

            ptr::null_mut()
        }
    }
}

/// Frobenius norm of a tensor.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_tensor_norm(tensor: *const ArrayD<f64>) -> f64 {

    if tensor.is_null() {

        return 0.0;
    }

    let t = unsafe {

        &*tensor
    };

    tensor::norm(t)
}
