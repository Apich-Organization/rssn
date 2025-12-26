//! Bincode-based FFI API for numerical finite field arithmetic.

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::finite_field::PrimeFieldElement;
use crate::numerical::finite_field::{
    self,
};

#[derive(Deserialize)]

struct PfeBinaryOpRequest {
    a: PrimeFieldElement,
    b: PrimeFieldElement,
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

/// GF(p) addition via Bincode.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_ff_pfe_add_bincode(
    data: *const u8,
    len: usize,
) -> BincodeBuffer {

    let req : PfeBinaryOpRequest = match decode(data, len) {
        | Some(r) => r,
        | None => {
            return encode(
                &FfiResult::<PrimeFieldElement, String> {
                    ok : None,
                    err : Some("Bincode decode error".to_string()),
                },
            )
        },
    };

    let res = req.a + req.b;

    encode(&FfiResult::<
        PrimeFieldElement,
        String,
    > {
        ok: Some(res),
        err: None,
    })
}

/// GF(p) multiplication via Bincode.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_ff_pfe_mul_bincode(
    data: *const u8,
    len: usize,
) -> BincodeBuffer {

    let req : PfeBinaryOpRequest = match decode(data, len) {
        | Some(r) => r,
        | None => {
            return encode(
                &FfiResult::<PrimeFieldElement, String> {
                    ok : None,
                    err : Some("Bincode decode error".to_string()),
                },
            )
        },
    };

    let res = req.a * req.b;

    encode(&FfiResult::<
        PrimeFieldElement,
        String,
    > {
        ok: Some(res),
        err: None,
    })
}
