//! Bincode-based FFI API for numerical combinatorics.

use crate::ffi_apis::common::{from_bincode_buffer, to_bincode_buffer, BincodeBuffer};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::combinatorics;
use serde::Deserialize;

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

#[no_mangle]

pub unsafe extern "C" fn rssn_num_comb_factorial_bincode(buffer: BincodeBuffer) -> BincodeBuffer {

    let input: NInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<f64, String> {
                ok: None,
                err: Some("Invalid Bincode input".to_string()),
            })
        }
    };

    let res = combinatorics::factorial(input.n);

    to_bincode_buffer(&FfiResult {
        ok: Some(res),
        err: None::<String>,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_comb_permutations_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {

    let input: NKInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<f64, String> {
                ok: None,
                err: Some("Invalid Bincode input".to_string()),
            })
        }
    };

    let res = combinatorics::permutations(input.n, input.k);

    to_bincode_buffer(&FfiResult {
        ok: Some(res),
        err: None::<String>,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_comb_combinations_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {

    let input: NKInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<f64, String> {
                ok: None,
                err: Some("Invalid Bincode input".to_string()),
            })
        }
    };

    let res = combinatorics::combinations(input.n, input.k);

    to_bincode_buffer(&FfiResult {
        ok: Some(res),
        err: None::<String>,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_comb_solve_recurrence_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {

    let input: RecurrenceInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<f64, String> {
                ok: None,
                err: Some("Invalid Bincode input".to_string()),
            })
        }
    };

    match combinatorics::solve_recurrence_numerical(
        &input.coeffs,
        &input.initial_conditions,
        input.target_n,
    ) {
        Ok(res) => to_bincode_buffer(&FfiResult {
            ok: Some(res),
            err: None::<String>,
        }),
        Err(e) => to_bincode_buffer(&FfiResult::<f64, String> {
            ok: None,
            err: Some(e),
        }),
    }
}

#[derive(Deserialize)]

struct XNInput {
    x: f64,
    n: u64,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_comb_stirling_second_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {

    let input: NKInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<f64, String> {
                ok: None,
                err: Some("Invalid Bincode input".to_string()),
            })
        }
    };

    let res = combinatorics::stirling_second(input.n, input.k);

    to_bincode_buffer(&FfiResult {
        ok: Some(res),
        err: None::<String>,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_comb_bell_bincode(buffer: BincodeBuffer) -> BincodeBuffer {

    let input: NInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<f64, String> {
                ok: None,
                err: Some("Invalid Bincode input".to_string()),
            })
        }
    };

    let res = combinatorics::bell(input.n);

    to_bincode_buffer(&FfiResult {
        ok: Some(res),
        err: None::<String>,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_comb_catalan_bincode(buffer: BincodeBuffer) -> BincodeBuffer {

    let input: NInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<f64, String> {
                ok: None,
                err: Some("Invalid Bincode input".to_string()),
            })
        }
    };

    let res = combinatorics::catalan(input.n);

    to_bincode_buffer(&FfiResult {
        ok: Some(res),
        err: None::<String>,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_comb_rising_factorial_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {

    let input: XNInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<f64, String> {
                ok: None,
                err: Some("Invalid Bincode input".to_string()),
            })
        }
    };

    let res = combinatorics::rising_factorial(input.x, input.n);

    to_bincode_buffer(&FfiResult {
        ok: Some(res),
        err: None::<String>,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_comb_falling_factorial_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {

    let input: XNInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<f64, String> {
                ok: None,
                err: Some("Invalid Bincode input".to_string()),
            })
        }
    };

    let res = combinatorics::falling_factorial(input.x, input.n);

    to_bincode_buffer(&FfiResult {
        ok: Some(res),
        err: None::<String>,
    })
}
