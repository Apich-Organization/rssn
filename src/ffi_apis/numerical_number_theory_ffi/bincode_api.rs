//! Bincode-based FFI API for numerical number theory operations.

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::number_theory as nt;

#[derive(Deserialize)]

struct FactorizeRequest {
    n : u64,
}

fn decode<
    T : for<'de> Deserialize<'de>,
>(
    data : *const u8,
    len : usize,
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

fn encode<T : Serialize>(
    val : &T
) -> BincodeBuffer {

    match bincode_next::serde::encode_to_vec(
        val,
        bincode_next::config::standard(),
    ) {
        | Ok(bytes) => BincodeBuffer::from_vec(bytes),
        | Err(_) => BincodeBuffer::empty(),
    }
}

/// Factorizes a number via Bincode.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_nt_factorize_bincode(
    data : *const u8,
    len : usize,
) -> BincodeBuffer {

    let req : FactorizeRequest = match decode(data, len) {
        | Some(r) => r,
        | None => {
            return encode(
                &FfiResult::<Vec<u64>, String> {
                    ok : None,
                    err : Some("Bincode decode error".to_string()),
                },
            )
        },
    };

    let factors = nt::factorize(req.n);

    encode(&FfiResult::<
        Vec<u64>,
        String,
    > {
        ok : Some(factors),
        err : None,
    })
}
