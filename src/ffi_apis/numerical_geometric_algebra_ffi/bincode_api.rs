//! Bincode-based FFI API for numerical geometric algebra operations.

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::BincodeBuffer;
use crate::numerical::geometric_algebra::Multivector3D;

#[derive(Deserialize)]

struct GaInput {
    mv: Multivector3D,
}

#[derive(Deserialize)]

struct TwoGaInput {
    mv1: Multivector3D,
    mv2: Multivector3D,
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

/// Bincode FFI for `ga_add`.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_ga_add_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: TwoGaInput = match decode(
        buffer,
    ) {
        | Some(v) => v,
        | None => {
            return encode(FfiResult::<
                Multivector3D,
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
        ok: Some(input.mv1 + input.mv2),
        err: None::<String>,
    })
}

/// Bincode FFI for `ga_sub`.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_ga_sub_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: TwoGaInput = match decode(
        buffer,
    ) {
        | Some(v) => v,
        | None => {
            return encode(FfiResult::<
                Multivector3D,
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
        ok: Some(input.mv1 - input.mv2),
        err: None::<String>,
    })
}

/// Bincode FFI for `ga_mul`.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_ga_mul_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: TwoGaInput = match decode(
        buffer,
    ) {
        | Some(v) => v,
        | None => {
            return encode(FfiResult::<
                Multivector3D,
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
        ok: Some(input.mv1 * input.mv2),
        err: None::<String>,
    })
}

/// Bincode FFI for `ga_wedge`.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_ga_wedge_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: TwoGaInput = match decode(
        buffer,
    ) {
        | Some(v) => v,
        | None => {
            return encode(FfiResult::<
                Multivector3D,
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
            input
                .mv1
                .wedge(input.mv2),
        ),
        err: None::<String>,
    })
}

/// Bincode FFI for `ga_dot`.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_ga_dot_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: TwoGaInput = match decode(
        buffer,
    ) {
        | Some(v) => v,
        | None => {
            return encode(FfiResult::<
                Multivector3D,
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
            input
                .mv1
                .dot(input.mv2),
        ),
        err: None::<String>,
    })
}

/// Bincode FFI for `ga_reverse`.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_ga_reverse_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: GaInput = match decode(
        buffer,
    ) {
        | Some(v) => v,
        | None => {
            return encode(FfiResult::<
                Multivector3D,
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
        ok: Some(input.mv.reverse()),
        err: None::<String>,
    })
}

/// Bincode FFI for `ga_norm`.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_ga_norm_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: GaInput = match decode(
        buffer,
    ) {
        | Some(v) => v,
        | None => {
            return encode(FfiResult::<
                f64,
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
        ok: Some(input.mv.norm()),
        err: None::<String>,
    })
}

/// Bincode FFI for `ga_inv`.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_ga_inv_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: GaInput = match decode(
        buffer,
    ) {
        | Some(v) => v,
        | None => {
            return encode(FfiResult::<
                Multivector3D,
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

    let res =
        match input.mv.inv() {
            | Some(v) => {
                FfiResult {
                    ok: Some(v),
                    err: None,
                }
            },
            | None => FfiResult {
                ok: None,
                err: Some(
                    "Multivector is \
                     not invertible"
                        .to_string(),
                ),
            },
        };

    encode(res)
}
