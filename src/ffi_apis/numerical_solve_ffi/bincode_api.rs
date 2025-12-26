//! Bincode-based FFI API for numerical equation solvers.

use crate::ffi_apis::common::BincodeBuffer;
use crate::numerical::matrix::Matrix;
use crate::numerical::solve::{self, LinearSolution};
use serde::{Deserialize, Serialize};

#[derive(Deserialize)]

struct SolveLinearInput {
    matrix: Matrix<f64>,
    vector: Vec<f64>,
}

#[derive(Serialize)]

struct FfiResult<T> {
    ok: Option<T>,
    err: Option<String>,
}

fn decode<T: for<'de> Deserialize<'de>>(buffer: BincodeBuffer) -> Option<T> {

    let slice = unsafe {

        buffer.as_slice()
    };

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

/// Bincode FFI for solving linear systems.
#[no_mangle]

pub unsafe extern "C" fn rssn_solve_linear_system_bincode(buffer: BincodeBuffer) -> BincodeBuffer {

    let input: SolveLinearInput = match decode(buffer) {
        Some(v) => v,
        None => {
            return encode(FfiResult::<LinearSolution> {
                ok: None,
                err: Some("Bincode decode error".to_string()),
            })
        }
    };

    match solve::solve_linear_system(&input.matrix, &input.vector) {
        Ok(sol) => encode(FfiResult {
            ok: Some(sol),
            err: None::<String>,
        }),
        Err(e) => encode(FfiResult::<LinearSolution> {
            ok: None,
            err: Some(e),
        }),
    }
}
