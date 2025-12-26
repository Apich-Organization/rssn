//! JSON-based FFI API for numerical functional analysis.

use std::os::raw::c_char;

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::functional_analysis;

#[derive(Deserialize)]

struct PointsInput {
    points : Vec<(f64, f64)>,
}

#[derive(Deserialize)]

struct InnerProductInput {
    f : Vec<(f64, f64)>,
    g : Vec<(f64, f64)>,
}

#[derive(Deserialize)]

struct GramSchmidtInput {
    basis : Vec<Vec<(f64, f64)>>,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_fa_l2_norm_json(input_json : *const c_char) -> *mut c_char {

    let input : PointsInput = match from_json_string(input_json) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON input".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let res = functional_analysis::l2_norm(&input.points);

    let ffi_res = FfiResult {
        ok : Some(res),
        err : None::<String>,
    };

    to_c_string(serde_json::to_string(&ffi_res).unwrap())
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_fa_inner_product_json(input_json : *const c_char) -> *mut c_char {

    let input : InnerProductInput = match from_json_string(input_json) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON input".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    match functional_analysis::inner_product(&input.f, &input.g) {
        | Ok(res) => {
            to_c_string(
                serde_json::to_string(&FfiResult {
                    ok : Some(res),
                    err : None::<String>,
                })
                .unwrap(),
            )
        },
        | Err(e) => {
            to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some(e),
                    },
                )
                .unwrap(),
            )
        },
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_fa_gram_schmidt_json(input_json : *const c_char) -> *mut c_char {

    let input : GramSchmidtInput = match from_json_string(input_json) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<Vec<(f64, f64)>>, String> {
                        ok : None,
                        err : Some("Invalid JSON input".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    match functional_analysis::gram_schmidt(&input.basis) {
        | Ok(res) => {
            to_c_string(
                serde_json::to_string(&FfiResult {
                    ok : Some(res),
                    err : None::<String>,
                })
                .unwrap(),
            )
        },
        | Err(e) => {
            to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<Vec<(f64, f64)>>, String> {
                        ok : None,
                        err : Some(e),
                    },
                )
                .unwrap(),
            )
        },
    }
}
