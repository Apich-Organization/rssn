//! Bincode-based FFI API for numerical vector operations.

use crate::ffi_apis::common::BincodeBuffer;
use crate::numerical::vector;
use serde::{Deserialize, Serialize};

#[derive(Deserialize)]

struct TwoVecInput {
    v1: Vec<f64>,
    v2: Vec<f64>,
}

#[derive(Deserialize)]

struct VecScalarInput {
    v: Vec<f64>,
    s: f64,
}

#[derive(Deserialize)]

struct VecNormInput {
    v: Vec<f64>,
    p: f64,
}

#[derive(Deserialize)]

struct VecEpsilonInput {
    v1: Vec<f64>,
    v2: Vec<f64>,
    epsilon: f64,
}

#[derive(Deserialize)]

struct LerpInput {
    v1: Vec<f64>,
    v2: Vec<f64>,
    t: f64,
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

    bincode_next::serde::decode_from_slice(
        slice,
        bincode_next::config::standard(),
    )
    .ok()
    .map(|(v, _)| v)
}

fn encode<T: Serialize>(val: T) -> BincodeBuffer {

    match bincode_next::serde::encode_to_vec(
        &val,
        bincode_next::config::standard(),
    ) {
        Ok(bytes) => BincodeBuffer::from_vec(bytes),
        Err(_) => BincodeBuffer::empty(),
    }
}

/// Bincode FFI for vec_add.
#[no_mangle]

pub unsafe extern "C" fn rssn_vec_add_bincode(buffer: BincodeBuffer) -> BincodeBuffer {

    let input: TwoVecInput = match decode(buffer) {
        Some(v) => v,
        None => {
            return encode(
                FfiResult::<Vec<f64>> {
                    ok: None,
                    err: Some("Bincode decode error".to_string()),
                },
            )
        }
    };

    let res = match vector::vec_add(&input.v1, &input.v2) {
        Ok(v) => {
            FfiResult {
                ok: Some(v),
                err: None,
            }
        }
        Err(e) => {
            FfiResult {
                ok: None,
                err: Some(e),
            }
        }
    };

    encode(res)
}

/// Bincode FFI for vec_sub.
#[no_mangle]

pub unsafe extern "C" fn rssn_vec_sub_bincode(buffer: BincodeBuffer) -> BincodeBuffer {

    let input: TwoVecInput = match decode(buffer) {
        Some(v) => v,
        None => {
            return encode(
                FfiResult::<Vec<f64>> {
                    ok: None,
                    err: Some("Bincode decode error".to_string()),
                },
            )
        }
    };

    let res = match vector::vec_sub(&input.v1, &input.v2) {
        Ok(v) => {
            FfiResult {
                ok: Some(v),
                err: None,
            }
        }
        Err(e) => {
            FfiResult {
                ok: None,
                err: Some(e),
            }
        }
    };

    encode(res)
}

/// Bincode FFI for scalar_mul.
#[no_mangle]

pub unsafe extern "C" fn rssn_vec_scalar_mul_bincode(buffer: BincodeBuffer) -> BincodeBuffer {

    let input: VecScalarInput = match decode(buffer) {
        Some(v) => v,
        None => {
            return encode(
                FfiResult::<Vec<f64>> {
                    ok: None,
                    err: Some("Bincode decode error".to_string()),
                },
            )
        }
    };

    let v = vector::scalar_mul(&input.v, input.s);

    encode(FfiResult {
        ok: Some(v),
        err: None,
    })
}

/// Bincode FFI for dot_product.
#[no_mangle]

pub unsafe extern "C" fn rssn_vec_dot_product_bincode(buffer: BincodeBuffer) -> BincodeBuffer {

    let input: TwoVecInput = match decode(buffer) {
        Some(v) => v,
        None => {
            return encode(FfiResult::<f64> {
                ok: None,
                err: Some("Bincode decode error".to_string()),
            })
        }
    };

    let res = match vector::dot_product(&input.v1, &input.v2) {
        Ok(v) => {
            FfiResult {
                ok: Some(v),
                err: None,
            }
        }
        Err(e) => {
            FfiResult {
                ok: None,
                err: Some(e),
            }
        }
    };

    encode(res)
}

/// Bincode FFI for norm.
#[no_mangle]

pub unsafe extern "C" fn rssn_vec_norm_bincode(buffer: BincodeBuffer) -> BincodeBuffer {

    let input: Vec<f64> = match decode(buffer) {
        Some(v) => v,
        None => {
            return encode(FfiResult::<f64> {
                ok: None,
                err: Some("Bincode decode error".to_string()),
            })
        }
    };

    let v = vector::norm(&input);

    encode(FfiResult {
        ok: Some(v),
        err: None,
    })
}

/// Bincode FFI for lp_norm.
#[no_mangle]

pub unsafe extern "C" fn rssn_vec_lp_norm_bincode(buffer: BincodeBuffer) -> BincodeBuffer {

    let input: VecNormInput = match decode(buffer) {
        Some(v) => v,
        None => {
            return encode(FfiResult::<f64> {
                ok: None,
                err: Some("Bincode decode error".to_string()),
            })
        }
    };

    let v = vector::lp_norm(&input.v, input.p);

    encode(FfiResult {
        ok: Some(v),
        err: None,
    })
}

/// Bincode FFI for normalize.
#[no_mangle]

pub unsafe extern "C" fn rssn_vec_normalize_bincode(buffer: BincodeBuffer) -> BincodeBuffer {

    let input: Vec<f64> = match decode(buffer) {
        Some(v) => v,
        None => {
            return encode(
                FfiResult::<Vec<f64>> {
                    ok: None,
                    err: Some("Bincode decode error".to_string()),
                },
            )
        }
    };

    let res = match vector::normalize(&input) {
        Ok(v) => {
            FfiResult {
                ok: Some(v),
                err: None,
            }
        }
        Err(e) => {
            FfiResult {
                ok: None,
                err: Some(e),
            }
        }
    };

    encode(res)
}

/// Bincode FFI for cross_product.
#[no_mangle]

pub unsafe extern "C" fn rssn_vec_cross_product_bincode(buffer: BincodeBuffer) -> BincodeBuffer {

    let input: TwoVecInput = match decode(buffer) {
        Some(v) => v,
        None => {
            return encode(
                FfiResult::<Vec<f64>> {
                    ok: None,
                    err: Some("Bincode decode error".to_string()),
                },
            )
        }
    };

    let res = match vector::cross_product(&input.v1, &input.v2) {
        Ok(v) => {
            FfiResult {
                ok: Some(v),
                err: None,
            }
        }
        Err(e) => {
            FfiResult {
                ok: None,
                err: Some(e),
            }
        }
    };

    encode(res)
}

/// Bincode FFI for distance.
#[no_mangle]

pub unsafe extern "C" fn rssn_vec_distance_bincode(buffer: BincodeBuffer) -> BincodeBuffer {

    let input: TwoVecInput = match decode(buffer) {
        Some(v) => v,
        None => {
            return encode(FfiResult::<f64> {
                ok: None,
                err: Some("Bincode decode error".to_string()),
            })
        }
    };

    let res = match vector::distance(&input.v1, &input.v2) {
        Ok(v) => {
            FfiResult {
                ok: Some(v),
                err: None,
            }
        }
        Err(e) => {
            FfiResult {
                ok: None,
                err: Some(e),
            }
        }
    };

    encode(res)
}

/// Bincode FFI for angle.
#[no_mangle]

pub unsafe extern "C" fn rssn_vec_angle_bincode(buffer: BincodeBuffer) -> BincodeBuffer {

    let input: TwoVecInput = match decode(buffer) {
        Some(v) => v,
        None => {
            return encode(FfiResult::<f64> {
                ok: None,
                err: Some("Bincode decode error".to_string()),
            })
        }
    };

    let res = match vector::angle(&input.v1, &input.v2) {
        Ok(v) => {
            FfiResult {
                ok: Some(v),
                err: None,
            }
        }
        Err(e) => {
            FfiResult {
                ok: None,
                err: Some(e),
            }
        }
    };

    encode(res)
}

/// Bincode FFI for project.
#[no_mangle]

pub unsafe extern "C" fn rssn_vec_project_bincode(buffer: BincodeBuffer) -> BincodeBuffer {

    let input: TwoVecInput = match decode(buffer) {
        Some(v) => v,
        None => {
            return encode(
                FfiResult::<Vec<f64>> {
                    ok: None,
                    err: Some("Bincode decode error".to_string()),
                },
            )
        }
    };

    let res = match vector::project(&input.v1, &input.v2) {
        Ok(v) => {
            FfiResult {
                ok: Some(v),
                err: None,
            }
        }
        Err(e) => {
            FfiResult {
                ok: None,
                err: Some(e),
            }
        }
    };

    encode(res)
}

/// Bincode FFI for reflect.
#[no_mangle]

pub unsafe extern "C" fn rssn_vec_reflect_bincode(buffer: BincodeBuffer) -> BincodeBuffer {

    let input: TwoVecInput = match decode(buffer) {
        Some(v) => v,
        None => {
            return encode(
                FfiResult::<Vec<f64>> {
                    ok: None,
                    err: Some("Bincode decode error".to_string()),
                },
            )
        }
    };

    let res = match vector::reflect(&input.v1, &input.v2) {
        Ok(v) => {
            FfiResult {
                ok: Some(v),
                err: None,
            }
        }
        Err(e) => {
            FfiResult {
                ok: None,
                err: Some(e),
            }
        }
    };

    encode(res)
}

/// Bincode FFI for lerp.
#[no_mangle]

pub unsafe extern "C" fn rssn_vec_lerp_bincode(buffer: BincodeBuffer) -> BincodeBuffer {

    let input: LerpInput = match decode(buffer) {
        Some(v) => v,
        None => {
            return encode(
                FfiResult::<Vec<f64>> {
                    ok: None,
                    err: Some("Bincode decode error".to_string()),
                },
            )
        }
    };

    let res = match vector::lerp(
        &input.v1, &input.v2, input.t,
    ) {
        Ok(v) => {
            FfiResult {
                ok: Some(v),
                err: None,
            }
        }
        Err(e) => {
            FfiResult {
                ok: None,
                err: Some(e),
            }
        }
    };

    encode(res)
}

/// Bincode FFI for is_orthogonal.
#[no_mangle]

pub unsafe extern "C" fn rssn_vec_is_orthogonal_bincode(buffer: BincodeBuffer) -> BincodeBuffer {

    let input: VecEpsilonInput = match decode(buffer) {
        Some(v) => v,
        None => {
            return encode(FfiResult::<bool> {
                ok: None,
                err: Some("Bincode decode error".to_string()),
            })
        }
    };

    let res = match vector::is_orthogonal(
        &input.v1,
        &input.v2,
        input.epsilon,
    ) {
        Ok(v) => {
            FfiResult {
                ok: Some(v),
                err: None,
            }
        }
        Err(e) => {
            FfiResult {
                ok: None,
                err: Some(e),
            }
        }
    };

    encode(res)
}

/// Bincode FFI for is_parallel.
#[no_mangle]

pub unsafe extern "C" fn rssn_vec_is_parallel_bincode(buffer: BincodeBuffer) -> BincodeBuffer {

    let input: VecEpsilonInput = match decode(buffer) {
        Some(v) => v,
        None => {
            return encode(FfiResult::<bool> {
                ok: None,
                err: Some("Bincode decode error".to_string()),
            })
        }
    };

    let res = match vector::is_parallel(
        &input.v1,
        &input.v2,
        input.epsilon,
    ) {
        Ok(v) => {
            FfiResult {
                ok: Some(v),
                err: None,
            }
        }
        Err(e) => {
            FfiResult {
                ok: None,
                err: Some(e),
            }
        }
    };

    encode(res)
}

/// Bincode FFI for cosine_similarity.
#[no_mangle]

pub unsafe extern "C" fn rssn_vec_cosine_similarity_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: TwoVecInput = match decode(buffer) {
        Some(v) => v,
        None => {
            return encode(FfiResult::<f64> {
                ok: None,
                err: Some("Bincode decode error".to_string()),
            })
        }
    };

    let res = match vector::cosine_similarity(&input.v1, &input.v2) {
        Ok(v) => {
            FfiResult {
                ok: Some(v),
                err: None,
            }
        }
        Err(e) => {
            FfiResult {
                ok: None,
                err: Some(e),
            }
        }
    };

    encode(res)
}
