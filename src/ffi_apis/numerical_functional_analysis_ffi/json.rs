//! JSON-based FFI API for numerical functional analysis.

use std::os::raw::c_char;

use serde::Deserialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::functional_analysis;

#[derive(Deserialize)]

struct PointsInput {
    points: Vec<(f64, f64)>,
}

#[derive(Deserialize)]

struct InnerProductInput {
    f: Vec<(f64, f64)>,
    g: Vec<(f64, f64)>,
}

#[derive(Deserialize)]

struct GramSchmidtInput {
    basis: Vec<Vec<(f64, f64)>>,
}

/// Computes the L2 norm of a function (represented by a series of points) using JSON for serialization.

#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_num_fa_l2_norm_json(
    input_json: *const c_char
) -> *mut c_char {

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

    let res =
        functional_analysis::l2_norm(
            &input.points,
        );

    let ffi_res = FfiResult {
        ok: Some(res),
        err: None::<String>,
    };

    to_c_string(
        serde_json::to_string(&ffi_res)
            .unwrap(),
    )
}

/// Computes the inner product of two functions (represented by series of points) using JSON for serialization.

#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_num_fa_inner_product_json(
    input_json: *const c_char
) -> *mut c_char {

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

/// Applies the Gram-Schmidt orthonormalization process to a set of basis functions using JSON for serialization.

#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_num_fa_gram_schmidt_json(
    input_json: *const c_char
) -> *mut c_char {

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
