//! Bincode-based FFI API for numerical real root finding.

use crate::ffi_apis::common::BincodeBuffer;
use crate::numerical::polynomial::Polynomial;
use crate::numerical::real_roots;
use serde::{Deserialize, Serialize};

#[derive(Deserialize)]
struct FindRootsInput {
    coeffs: Vec<f64>,
    tolerance: f64,
}

#[derive(Serialize)]
struct FfiResult<T> {
    ok: Option<T>,
    err: Option<String>,
}

fn decode<T: for<'de> Deserialize<'de>>(buffer: BincodeBuffer) -> Option<T> {
    let slice = unsafe { buffer.as_slice() };
    bincode_next::serde::decode_from_slice(slice, bincode_next::config::standard())
        .ok()
        .map(|(v, _)| v)
}

fn encode<T: Serialize>(val: T) -> BincodeBuffer {
    match bincode_next::serde::encode_to_vec(&val, bincode_next::config::standard()) {
        Ok(bytes) => BincodeBuffer::from_vec(bytes),
        Err(_) => BincodeBuffer::empty(),
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_real_roots_find_roots_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {
    let input: FindRootsInput = match decode(buffer) {
        Some(v) => v,
        None => {
            return encode(FfiResult::<Vec<f64>> {
                ok: None,
                err: Some("Bincode decode error".to_string()),
            })
        }
    };

    let poly = Polynomial::new(input.coeffs);
    match real_roots::find_roots(&poly, input.tolerance) {
        Ok(roots) => encode(FfiResult {
            ok: Some(roots),
            err: None::<String>,
        }),
        Err(e) => encode(FfiResult::<Vec<f64>> {
            ok: None,
            err: Some(e),
        }),
    }
}
