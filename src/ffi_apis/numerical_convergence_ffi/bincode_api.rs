//! Bincode-based FFI API for numerical convergence operations.

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::BincodeBuffer;
use crate::numerical::convergence;

#[derive(Deserialize)]

struct SeqInput {
    sequence: Vec<f64>,
}

#[derive(Serialize)]

struct FfiResult<T> {
    ok: Option<T>,
    err: Option<String>,
}

fn decode<
    T: for<'de> Deserialize<'de>,
>(
    buffer: BincodeBuffer
) -> Option<T> {

    let slice = unsafe {

        buffer.as_slice()
    };

    bincode_next::serde::decode_from_slice(
        slice,
        bincode_next::config::standard(),
    )
    .ok()
    .map(|(v, _)| v)
}

fn encode<T: Serialize>(
    val: T
) -> BincodeBuffer {

    match bincode_next::serde::encode_to_vec(
        &val,
        bincode_next::config::standard(),
    ) {
        | Ok(bytes) => BincodeBuffer::from_vec(bytes),
        | Err(_) => BincodeBuffer::empty(),
    }
}

/// Bincode FFI for Aitken acceleration.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_convergence_aitken_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: SeqInput = match decode(
        buffer,
    ) {
        | Some(v) => v,
        | None => {
            return encode(FfiResult::<
                Vec<f64>,
            > {
                ok: None,
                err: Some(
                    "Bincode decode \
                     error"
                        .to_string(),
                ),
            })
        },
    };

    encode(FfiResult {
        ok : Some(convergence::aitken_acceleration(&input.sequence)),
        err : None::<String>,
    })
}

/// Bincode FFI for Richardson extrapolation.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_convergence_richardson_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: SeqInput = match decode(
        buffer,
    ) {
        | Some(v) => v,
        | None => {
            return encode(FfiResult::<
                Vec<f64>,
            > {
                ok: None,
                err: Some(
                    "Bincode decode \
                     error"
                        .to_string(),
                ),
            })
        },
    };

    encode(FfiResult {
        ok : Some(convergence::richardson_extrapolation(&input.sequence)),
        err : None::<String>,
    })
}

/// Bincode FFI for Wynn's epsilon algorithm.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_convergence_wynn_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: SeqInput = match decode(
        buffer,
    ) {
        | Some(v) => v,
        | None => {
            return encode(FfiResult::<
                Vec<f64>,
            > {
                ok: None,
                err: Some(
                    "Bincode decode \
                     error"
                        .to_string(),
                ),
            })
        },
    };

    encode(FfiResult {
        ok: Some(
            convergence::wynn_epsilon(
                &input.sequence,
            ),
        ),
        err: None::<String>,
    })
}
