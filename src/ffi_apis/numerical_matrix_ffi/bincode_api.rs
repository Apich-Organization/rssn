//! Bincode-based FFI API for numerical matrix operations.

use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::matrix::Matrix;
use serde::{Deserialize, Serialize};

#[derive(Deserialize)]
struct MatrixOpRequest {
    m1: Matrix<f64>,
    m2: Option<Matrix<f64>>,
}

fn decode<T: for<'de> Deserialize<'de>>(data: *const u8, len: usize) -> Option<T> {
    if data.is_null() {
        return None;
    }
    let slice = unsafe { std::slice::from_raw_parts(data, len) };
    bincode_next::serde::decode_from_slice(slice, bincode_next::config::standard())
        .ok()
        .map(|(v, _)| v)
}

fn encode<T: Serialize>(val: &T) -> BincodeBuffer {
    match bincode_next::serde::encode_to_vec(val, bincode_next::config::standard()) {
        Ok(bytes) => BincodeBuffer::from_vec(bytes),
        Err(_) => BincodeBuffer::empty(),
    }
}

/// Matrix addition via Bincode.
#[no_mangle]
pub unsafe extern "C" fn rssn_num_matrix_add_bincode(data: *const u8, len: usize) -> BincodeBuffer {
    let req: MatrixOpRequest = match decode(data, len) {
        Some(r) => r,
        None => {
            return encode(&FfiResult::<Matrix<f64>, String> {
                ok: None,
                err: Some("Bincode decode error".to_string()),
            })
        }
    };

    let m2 = match req.m2 {
        Some(m) => m,
        None => {
            return encode(&FfiResult::<Matrix<f64>, String> {
                ok: None,
                err: Some("m2 required".to_string()),
            })
        }
    };

    if req.m1.rows() != m2.rows() || req.m1.cols() != m2.cols() {
        return encode(&FfiResult::<Matrix<f64>, String> {
            ok: None,
            err: Some("Dimension mismatch".to_string()),
        });
    }

    let result = req.m1 + m2;
    encode(&FfiResult::<Matrix<f64>, String> {
        ok: Some(result),
        err: None,
    })
}

/// Matrix multiplication via Bincode.
#[no_mangle]
pub unsafe extern "C" fn rssn_num_matrix_mul_bincode(data: *const u8, len: usize) -> BincodeBuffer {
    let req: MatrixOpRequest = match decode(data, len) {
        Some(r) => r,
        None => {
            return encode(&FfiResult::<Matrix<f64>, String> {
                ok: None,
                err: Some("Bincode decode error".to_string()),
            })
        }
    };

    let m2 = match req.m2 {
        Some(m) => m,
        None => {
            return encode(&FfiResult::<Matrix<f64>, String> {
                ok: None,
                err: Some("m2 required".to_string()),
            })
        }
    };

    if req.m1.cols() != m2.rows() {
        return encode(&FfiResult::<Matrix<f64>, String> {
            ok: None,
            err: Some("Dimension mismatch".to_string()),
        });
    }

    let result = req.m1 * m2;
    encode(&FfiResult::<Matrix<f64>, String> {
        ok: Some(result),
        err: None,
    })
}
