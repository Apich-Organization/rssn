//! JSON-based FFI API for numerical transforms.

use std::os::raw::c_char;

use num_complex::Complex;
use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::transforms;

#[derive(Deserialize)]

struct TransformInput {
    data : Vec<Complex<f64>>,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_fft_json(
    input_json : *const c_char
) -> *mut c_char {

    let mut input: TransformInput =
        match from_json_string(input_json) {
            | Some(i) => i,
            | None => {
                return to_c_string(
                    serde_json::to_string(&FfiResult::<
                        Vec<Complex<f64>>,
                        String,
                    > {
                        ok: None,
                        err: Some(
                            "Invalid JSON input"
                                .to_string(),
                        ),
                    })
                    .unwrap(),
                )
            },
        };

    transforms::fft(&mut input.data);

    let ffi_res = FfiResult {
        ok : Some(input.data),
        err : None::<String>,
    };

    to_c_string(
        serde_json::to_string(&ffi_res)
            .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_ifft_json(
    input_json : *const c_char
) -> *mut c_char {

    let mut input: TransformInput =
        match from_json_string(input_json) {
            | Some(i) => i,
            | None => {
                return to_c_string(
                    serde_json::to_string(&FfiResult::<
                        Vec<Complex<f64>>,
                        String,
                    > {
                        ok: None,
                        err: Some(
                            "Invalid JSON input"
                                .to_string(),
                        ),
                    })
                    .unwrap(),
                )
            },
        };

    transforms::ifft(&mut input.data);

    let ffi_res = FfiResult {
        ok : Some(input.data),
        err : None::<String>,
    };

    to_c_string(
        serde_json::to_string(&ffi_res)
            .unwrap(),
    )
}
