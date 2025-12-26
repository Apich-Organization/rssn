//! Bincode-based FFI API for numerical sparse matrix operations.

use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::sparse::{self, SparseMatrixData};
use serde::{Deserialize, Serialize};

#[derive(Deserialize)]

struct SpMvRequest {
    matrix: SparseMatrixData,
    vector: Vec<f64>,
}

fn decode<T: for<'de> Deserialize<'de>>(
    data: *const u8,
    len: usize,
) -> Option<T> {

    if data.is_null() {

        return None;
    }

    let slice = unsafe {

        std::slice::from_raw_parts(data, len)
    };

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

/// Sparse matrix-vector multiplication via Bincode.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_sparse_spmv_bincode(
    data: *const u8,
    len: usize,
) -> BincodeBuffer {

    let req: SpMvRequest = match decode(data, len) {
        Some(r) => r,
        None => {
            return encode(&FfiResult::<Vec<f64>, String> {
                ok: None,
                err: Some("Bincode decode error".to_string()),
            })
        }
    };

    let mat = req
        .matrix
        .to_csmat();

    match sparse::sp_mat_vec_mul(&mat, &req.vector) {
        Ok(res) => encode(&FfiResult::<Vec<f64>, String> {
            ok: Some(res),
            err: None,
        }),
        Err(e) => encode(&FfiResult::<Vec<f64>, String> {
            ok: None,
            err: Some(e),
        }),
    }
}
