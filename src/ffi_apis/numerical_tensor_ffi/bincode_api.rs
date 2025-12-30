//! Bincode-based FFI API for numerical tensor operations.

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::tensor::TensorData;
use crate::numerical::tensor::{
    self,
};

#[derive(Deserialize)]

struct TensordotRequest {
    a: TensorData,
    b: TensorData,
    axes_a: Vec<usize>,
    axes_b: Vec<usize>,
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

/// Tensor contraction via Bincode.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_tensor_tensordot_bincode(
    data: *const u8,
    len: usize,
) -> BincodeBuffer {

    let req : TensordotRequest = match decode(data, len) {
        | Some(r) => r,
        | None => {
            return encode(
                &FfiResult::<TensorData, String> {
                    ok : None,
                    err : Some("Bincode decode error".to_string()),
                },
            )
        },
    };

    let a = match req.a.to_arrayd() {
        | Ok(arr) => arr,
        | Err(e) => {
            return encode(
                &FfiResult::<
                    TensorData,
                    String,
                > {
                    ok: None,
                    err: Some(e),
                },
            )
        },
    };

    let b = match req.b.to_arrayd() {
        | Ok(arr) => arr,
        | Err(e) => {
            return encode(
                &FfiResult::<
                    TensorData,
                    String,
                > {
                    ok: None,
                    err: Some(e),
                },
            )
        },
    };

    match tensor::tensordot(
        &a,
        &b,
        &req.axes_a,
        &req.axes_b,
    ) {
        | Ok(res) => {
            encode(&FfiResult::<
                TensorData,
                String,
            > {
                ok: Some(
                    TensorData::from(
                        &res,
                    ),
                ),
                err: None,
            })
        },
        | Err(e) => {
            encode(&FfiResult::<
                TensorData,
                String,
            > {
                ok: None,
                err: Some(e),
            })
        },
    }
}
