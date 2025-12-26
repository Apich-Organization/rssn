//! Bincode-based FFI API for numerical transforms.

use crate::ffi_apis::common::{from_bincode_buffer, to_bincode_buffer, BincodeBuffer};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::transforms;
use num_complex::Complex;
use serde::{Deserialize, Serialize};

#[derive(Deserialize)]

struct TransformInput {
    data: Vec<Complex<f64>>,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_fft_bincode(buffer: BincodeBuffer) -> BincodeBuffer {

    let mut input: TransformInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<Vec<Complex<f64>>, String> {
                ok: None,
                err: Some("Invalid Bincode input".to_string()),
            })
        }
    };

    transforms::fft(&mut input.data);

    let ffi_res = FfiResult {
        ok: Some(input.data),
        err: None::<String>,
    };

    to_bincode_buffer(&ffi_res)
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_ifft_bincode(buffer: BincodeBuffer) -> BincodeBuffer {

    let mut input: TransformInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<Vec<Complex<f64>>, String> {
                ok: None,
                err: Some("Invalid Bincode input".to_string()),
            })
        }
    };

    transforms::ifft(&mut input.data);

    let ffi_res = FfiResult {
        ok: Some(input.data),
        err: None::<String>,
    };

    to_bincode_buffer(&ffi_res)
}
