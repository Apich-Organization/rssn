//! Bincode-based FFI API for numerical polynomial operations.

use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::polynomial::Polynomial;
use serde::{Deserialize, Serialize};

#[derive(Deserialize)]

struct PolyBinaryOpRequest {
    a: Polynomial,
    b: Polynomial,
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

    bincode_next::serde::decode_from_slice(
        slice,
        bincode_next::config::standard(),
    )
    .ok()
    .map(|(v, _)| v)
}

fn encode<T: Serialize>(val: &T) -> BincodeBuffer {

    match bincode_next::serde::encode_to_vec(
        val,
        bincode_next::config::standard(),
    ) {
        Ok(bytes) => BincodeBuffer::from_vec(bytes),
        Err(_) => BincodeBuffer::empty(),
    }
}

/// Adds two polynomials via Bincode.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_poly_add_bincode(
    data: *const u8,
    len: usize,
) -> BincodeBuffer {

    let req: PolyBinaryOpRequest = match decode(data, len) {
        Some(r) => r,
        None => {
            return encode(
                &FfiResult::<Polynomial, String> {
                    ok: None,
                    err: Some("Bincode decode error".to_string()),
                },
            )
        }
    };

    let res = req.a + req.b;

    encode(
        &FfiResult::<Polynomial, String> {
            ok: Some(res),
            err: None,
        },
    )
}
