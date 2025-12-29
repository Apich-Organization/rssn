//! Bincode-based FFI API for numerical sparse matrix operations.

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::sparse::SparseMatrixData;
use crate::numerical::sparse::{
    self,
};

#[derive(Deserialize)]

struct SpMvRequest {
    matrix: SparseMatrixData,
    vector: Vec<f64>,
}

fn decode<
    T: for<'de> Deserialize<'de>,
>(
    data: *const u8,
    len: usize,
) -> Option<T> {

    if data.is_null() {

        return None;
    }

    let slice = unsafe {

        std::slice::from_raw_parts(
            data, len,
        )
    };

    bincode_next::serde::decode_from_slice(
        slice,
        bincode_next::config::standard(),
    )
    .ok()
    .map(|(v, _)| v)
}

fn encode<T: Serialize>(
    val: &T
) -> BincodeBuffer {

    match bincode_next::serde::encode_to_vec(
        val,
        bincode_next::config::standard(),
    ) {
        | Ok(bytes) => BincodeBuffer::from_vec(bytes),
        | Err(_) => BincodeBuffer::empty(),
    }
}

/// Sparse matrix-vector multiplication via Bincode.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_sparse_spmv_bincode(
    data: *const u8,
    len: usize,
) -> BincodeBuffer {

    let req: SpMvRequest = match decode(
        data, len,
    ) {
        | Some(r) => r,
        | None => return encode(
            &FfiResult::<
                Vec<f64>,
                String,
            > {
                ok: None,
                err: Some(
                    "Bincode decode \
                     error"
                        .to_string(),
                ),
            },
        ),
    };

    let mat = req
        .matrix
        .to_csmat();

    match sparse::sp_mat_vec_mul(
        &mat,
        &req.vector,
    ) {
        | Ok(res) => {
            encode(&FfiResult::<
                Vec<f64>,
                String,
            > {
                ok: Some(res),
                err: None,
            })
        },
        | Err(e) => {
            encode(&FfiResult::<
                Vec<f64>,
                String,
            > {
                ok: None,
                err: Some(e),
            })
        },
    }
}
