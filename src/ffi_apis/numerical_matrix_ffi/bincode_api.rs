//! Bincode-based FFI API for numerical matrix operations.

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::matrix::Backend;
use crate::numerical::matrix::FaerDecompositionResult;
use crate::numerical::matrix::Matrix;

#[derive(Deserialize)]

struct MatrixOpRequest {
    m1: Matrix<f64>,
    m2: Option<Matrix<f64>>,
}

#[derive(Deserialize)]

struct MatrixBackendRequest {
    matrix: Matrix<f64>,
    backend_id: i32, /* 0: Native, 1: Faer */
}

#[derive(Deserialize)]

struct MatrixDecompositionRequest {
    matrix: Matrix<f64>,
    kind: crate::numerical::matrix::FaerDecompositionType,
}

fn decode<
    T: for<'de> Deserialize<'de>,
>(
    data: *const u8,
    len: usize,
) -> Option<T> {

    if data.is_null() {

        return None;
    }

    let slice = unsafe {

        std::slice::from_raw_parts(
            data, len,
        )
    };

    bincode_next::serde::decode_from_slice(slice, bincode_next::config::standard())
        .ok()
        .map(|(v, _)| v)
}

fn encode<T: Serialize>(
    val: &T
) -> BincodeBuffer {

    match bincode_next::serde::encode_to_vec(val, bincode_next::config::standard()) {
        Ok(bytes) => BincodeBuffer::from_vec(bytes),
        Err(_) => BincodeBuffer::empty(),
    }
}

/// Matrix addition via Bincode.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_matrix_add_bincode(
    data: *const u8,
    len: usize,
) -> BincodeBuffer {

    let req: MatrixOpRequest = match decode(data, len) {
        Some(r) => r,
        None => {
            return encode(&FfiResult::<Matrix<f64>, String> {
                ok: None,
                err: Some("Bincode decode error".to_string()),
            })
        }
    };

    let m2 =
        match req.m2 {
            | Some(m) => m,
            | None => return encode(
                &FfiResult::<
                    Matrix<f64>,
                    String,
                > {
                    ok: None,
                    err: Some(
                        "m2 required"
                            .to_string(
                            ),
                    ),
                },
            ),
        };

    if req.m1.rows() != m2.rows()
        || req.m1.cols() != m2.cols()
    {

        return encode(&FfiResult::<
            Matrix<f64>,
            String,
        > {
            ok: None,
            err: Some(
                "Dimension mismatch"
                    .to_string(),
            ),
        });
    }

    let result = req.m1 + m2;

    encode(&FfiResult::<
        Matrix<f64>,
        String,
    > {
        ok: Some(result),
        err: None,
    })
}

/// Matrix multiplication via Bincode.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_matrix_mul_bincode(
    data: *const u8,
    len: usize,
) -> BincodeBuffer {

    let req: MatrixOpRequest = match decode(data, len) {
        Some(r) => r,
        None => {
            return encode(&FfiResult::<Matrix<f64>, String> {
                ok: None,
                err: Some("Bincode decode error".to_string()),
            })
        }
    };

    let m2 =
        match req.m2 {
            | Some(m) => m,
            | None => return encode(
                &FfiResult::<
                    Matrix<f64>,
                    String,
                > {
                    ok: None,
                    err: Some(
                        "m2 required"
                            .to_string(
                            ),
                    ),
                },
            ),
        };

    if req.m1.cols() != m2.rows() {

        return encode(&FfiResult::<
            Matrix<f64>,
            String,
        > {
            ok: None,
            err: Some(
                "Dimension mismatch"
                    .to_string(),
            ),
        });
    }

    let result = req.m1 * m2;

    encode(&FfiResult::<
        Matrix<f64>,
        String,
    > {
        ok: Some(result),
        err: None,
    })
}

/// Sets backend for a matrix via Bincode.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_matrix_set_backend_bincode(
    data: *const u8,
    len: usize,
) -> BincodeBuffer {

    let req: MatrixBackendRequest = match decode(data, len) {
        Some(r) => r,
        None => return encode(&FfiResult::<Matrix<f64>, String> { ok: None, err: Some("Bincode decode error".to_string()) }),
    };

    let backend = match req.backend_id {
        | 1 => Backend::Faer,
        | _ => Backend::Native,
    };

    let m = req
        .matrix
        .with_backend(backend);

    encode(&FfiResult::<
        Matrix<f64>,
        String,
    > {
        ok: Some(m),
        err: None,
    })
}

/// Decomposes a matrix via Bincode.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_matrix_decompose_bincode(
    data: *const u8,
    len: usize,
) -> BincodeBuffer {

    let req: MatrixDecompositionRequest = match decode(data, len) {
        Some(r) => r,
        None => return encode(&FfiResult::<FaerDecompositionResult<f64>, String> { ok: None, err: Some("Bincode decode error".to_string()) }),
    };

    match req
        .matrix
        .decompose(req.kind)
    {
        | Some(result) => {
            encode(&FfiResult::<
                FaerDecompositionResult<
                    f64,
                >,
                String,
            > {
                ok: Some(result),
                err: None,
            })
        },
        | None => {
            encode(&FfiResult::<
                FaerDecompositionResult<
                    f64,
                >,
                String,
            > {
                ok: None,
                err: Some(
                    "Decomposition \
                     failed or backend \
                     not supported"
                        .to_string(),
                ),
            })
        },
    }
}
