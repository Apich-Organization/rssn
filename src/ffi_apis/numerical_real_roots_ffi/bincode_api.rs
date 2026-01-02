//! Bincode-based FFI API for numerical real root finding.

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::BincodeBuffer;
use crate::numerical::polynomial::Polynomial;
use crate::numerical::real_roots;

#[derive(Deserialize)]

struct FindRootsInput {
    coeffs: Vec<f64>,
    tolerance: f64,
}

#[derive(Serialize)]

struct FfiResult<T> {
    ok: Option<T>,
    err: Option<String>,
}

pub(crate) fn decode<
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

/// Finds all real roots of a polynomial using numerical methods and bincode serialization.
///
/// Uses root-finding algorithms to locate all real zeros of the polynomial.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `FindRootsInput` with:
///   - `coeffs`: Polynomial coefficients [a₀, a₁, ..., aₙ] for a₀ + a₁x + ... + aₙxⁿ
///   - `tolerance`: Convergence tolerance for root finding
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<Vec<f64>, String>` with either:
/// - `ok`: Array of real roots found
/// - `err`: Error message if root finding failed
///
/// # Safety
///
/// This function is unsafe because it receives raw pointers through FFI.
/// The caller must ensure the input buffer contains valid bincode data.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_real_roots_find_roots_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : FindRootsInput = match decode(buffer) {
        | Some(v) => v,
        | None => {
            return encode(
                FfiResult::<Vec<f64>> {
                    ok : None,
                    err : Some("Bincode decode error".to_string()),
                },
            )
        },
    };

    let poly =
        Polynomial::new(input.coeffs);

    match real_roots::find_roots(
        &poly,
        input.tolerance,
    ) {
        | Ok(roots) => {
            encode(FfiResult {
                ok: Some(roots),
                err: None::<String>,
            })
        },
        | Err(e) => {
            encode(FfiResult::<
                Vec<f64>,
            > {
                ok: None,
                err: Some(e),
            })
        },
    }
}
