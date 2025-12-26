//! JSON-based FFI API for numerical signal processing.

use std::os::raw::c_char;

use rustfft::num_complex::Complex;
use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::signal;

#[derive(Deserialize)]

struct FftInput {
    data : Vec<Complex<f64>>,
}

#[derive(Deserialize)]

struct ConvolveInput {
    a : Vec<f64>,
    v : Vec<f64>,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_signal_fft_json(
    input_json : *const c_char
) -> *mut c_char {

    let mut input: FftInput =
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

    let result =
        signal::fft(&mut input.data);

    let ffi_res = FfiResult {
        ok : Some(result),
        err : None::<String>,
    };

    to_c_string(
        serde_json::to_string(&ffi_res)
            .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_signal_convolve_json(
    input_json : *const c_char
) -> *mut c_char {

    let input: ConvolveInput =
        match from_json_string(input_json) {
            | Some(i) => i,
            | None => {
                return to_c_string(
                    serde_json::to_string(&FfiResult::<
                        Vec<f64>,
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

    let result = signal::convolve(
        &input.a,
        &input.v,
    );

    let ffi_res = FfiResult {
        ok : Some(result),
        err : None::<String>,
    };

    to_c_string(
        serde_json::to_string(&ffi_res)
            .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_signal_cross_correlation_json(
    input_json : *const c_char
) -> *mut c_char {

    let input: ConvolveInput =
        match from_json_string(input_json) {
            | Some(i) => i,
            | None => {
                return to_c_string(
                    serde_json::to_string(&FfiResult::<
                        Vec<f64>,
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

    let result =
        signal::cross_correlation(
            &input.a,
            &input.v,
        );

    let ffi_res = FfiResult {
        ok : Some(result),
        err : None::<String>,
    };

    to_c_string(
        serde_json::to_string(&ffi_res)
            .unwrap(),
    )
}
