//! Bincode-based FFI API for numerical transforms.

use num_complex::Complex;
use serde::Deserialize;

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::transforms;

#[derive(Deserialize)]

struct TransformInput {
    data: Vec<Complex<f64>>,
}

/// Computes the Fast Fourier Transform (FFT) in-place using bincode serialization.
///
/// The FFT converts a sequence from the time/space domain to the frequency domain,
/// computing X(k) = Σx(n)e^(-2πikn/N) for k = 0, ..., N-1.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `TransformInput` with:
///   - `data`: Vector of complex numbers to transform
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<Vec<Complex<f64>>, String>` with either:
/// - `ok`: FFT-transformed data in frequency domain
/// - `err`: Error message if computation failed
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_fft_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let mut input : TransformInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<Vec<Complex<f64>>, String> {
                    ok : None,
                    err : Some("Invalid Bincode input".to_string()),
                },
            )
        },
    };

    transforms::fft(&mut input.data);

    let ffi_res = FfiResult {
        ok: Some(input.data),
        err: None::<String>,
    };

    to_bincode_buffer(&ffi_res)
}

/// Computes the Inverse Fast Fourier Transform (IFFT) in-place using bincode serialization.
///
/// The IFFT converts a sequence from the frequency domain back to the time/space domain,
/// computing x(n) = (1/N) ΣX(k)e^(2πikn/N) for n = 0, ..., N-1.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `TransformInput` with:
///   - `data`: Vector of complex numbers in frequency domain
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<Vec<Complex<f64>>, String>` with either:
/// - `ok`: IFFT-transformed data in time/space domain
/// - `err`: Error message if computation failed
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_ifft_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let mut input : TransformInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<Vec<Complex<f64>>, String> {
                    ok : None,
                    err : Some("Invalid Bincode input".to_string()),
                },
            )
        },
    };

    transforms::ifft(&mut input.data);

    let ffi_res = FfiResult {
        ok: Some(input.data),
        err: None::<String>,
    };

    to_bincode_buffer(&ffi_res)
}
