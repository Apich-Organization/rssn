//! JSON-based FFI API for numerical combinatorics.

use std::os::raw::c_char;

use serde::Deserialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::combinatorics;

#[derive(Deserialize)]

struct NInput {
    n: u64,
}

#[derive(Deserialize)]

struct NKInput {
    n: u64,
    k: u64,
}

#[derive(Deserialize)]

struct RecurrenceInput {
    coeffs: Vec<f64>,
    initial_conditions: Vec<f64>,
    target_n: usize,
}

/// Computes the factorial of a number using JSON for serialization.

#[no_mangle]

pub unsafe extern "C" fn rssn_num_comb_factorial_json(
    input_json: *const c_char
) -> *mut c_char {

    let input : NInput = match from_json_string(input_json) {
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

    let res = combinatorics::factorial(
        input.n,
    );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(res),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Computes the number of permutations using JSON for serialization.

#[no_mangle]

pub unsafe extern "C" fn rssn_num_comb_permutations_json(
    input_json: *const c_char
) -> *mut c_char {

    let input : NKInput = match from_json_string(input_json) {
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
        combinatorics::permutations(
            input.n,
            input.k,
        );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(res),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Computes the number of combinations using JSON for serialization.

#[no_mangle]

pub unsafe extern "C" fn rssn_num_comb_combinations_json(
    input_json: *const c_char
) -> *mut c_char {

    let input : NKInput = match from_json_string(input_json) {
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
        combinatorics::combinations(
            input.n,
            input.k,
        );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(res),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Solves a linear recurrence relation using JSON for serialization.

#[no_mangle]

pub unsafe extern "C" fn rssn_num_comb_solve_recurrence_json(
    input_json: *const c_char
) -> *mut c_char {

    let input : RecurrenceInput = match from_json_string(input_json) {
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

    match combinatorics::solve_recurrence_numerical(
        &input.coeffs,
        &input.initial_conditions,
        input.target_n,
    ) {
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

#[derive(Deserialize)]

struct XNInput {
    x: f64,
    n: u64,
}

/// Computes the Stirling number of the second kind using JSON for serialization.

#[no_mangle]

pub unsafe extern "C" fn rssn_num_comb_stirling_second_json(
    input_json: *const c_char
) -> *mut c_char {

    let input : NKInput = match from_json_string(input_json) {
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
        combinatorics::stirling_second(
            input.n,
            input.k,
        );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(res),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Computes the Bell number using JSON for serialization.

#[no_mangle]

pub unsafe extern "C" fn rssn_num_comb_bell_json(
    input_json: *const c_char
) -> *mut c_char {

    let input : NInput = match from_json_string(input_json) {
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
        combinatorics::bell(input.n);

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(res),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Computes the Catalan number using JSON for serialization.

#[no_mangle]

pub unsafe extern "C" fn rssn_num_comb_catalan_json(
    input_json: *const c_char
) -> *mut c_char {

    let input : NInput = match from_json_string(input_json) {
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
        combinatorics::catalan(input.n);

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(res),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Computes the rising factorial using JSON for serialization.

#[no_mangle]

pub unsafe extern "C" fn rssn_num_comb_rising_factorial_json(
    input_json: *const c_char
) -> *mut c_char {

    let input : XNInput = match from_json_string(input_json) {
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
        combinatorics::rising_factorial(
            input.x,
            input.n,
        );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(res),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Computes the falling factorial using JSON for serialization.

#[no_mangle]

pub unsafe extern "C" fn rssn_num_comb_falling_factorial_json(
    input_json: *const c_char
) -> *mut c_char {

    let input : XNInput = match from_json_string(input_json) {
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

    let res = combinatorics::falling_factorial(input.x, input.n);

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(res),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}
