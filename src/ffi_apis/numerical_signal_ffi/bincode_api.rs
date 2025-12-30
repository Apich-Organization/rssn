//! Bincode-based FFI API for numerical signal processing.

use rustfft::num_complex::Complex;
use serde::Deserialize;

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::signal;

#[derive(Deserialize)]

struct FftInput {
    data: Vec<Complex<f64>>,
}

#[derive(Deserialize)]

struct ConvolveInput {
    a: Vec<f64>,
    v: Vec<f64>,
}

/// Computes the Fast Fourier Transform (FFT) of complex data using bincode serialization.
///
/// The FFT converts a signal from time/space domain to frequency domain using the
/// Cooley-Tukey algorithm in O(N log N) time.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `FftInput` with:
///   - `data`: Array of complex numbers to transform
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<Vec<Complex<f64>>, String>` with
/// the frequency domain representation.
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

pub unsafe extern "C" fn rssn_num_signal_fft_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let mut input : FftInput = match from_bincode_buffer(&buffer) {
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

    let result =
        signal::fft(&mut input.data);

    let ffi_res = FfiResult {
        ok: Some(result),
        err: None::<String>,
    };

    to_bincode_buffer(&ffi_res)
}

/// Computes the discrete convolution of two signals using bincode serialization.
///
/// The convolution is defined as (a * v)[n] = Σ a[k]v[n-k], representing the
/// combined effect of two systems or filtering operation.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `ConvolveInput` with:
///   - `a`: First signal array
///   - `v`: Second signal array (kernel)
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<Vec<f64>, String>` with
/// the convolution result.
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

pub unsafe extern "C" fn rssn_num_signal_convolve_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : ConvolveInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<Vec<f64>, String> {
                    ok : None,
                    err : Some("Invalid Bincode input".to_string()),
                },
            )
        },
    };

    let result = signal::convolve(
        &input.a,
        &input.v,
    );

    let ffi_res = FfiResult {
        ok: Some(result),
        err: None::<String>,
    };

    to_bincode_buffer(&ffi_res)
}

/// Computes the cross-correlation of two signals using bincode serialization.
///
/// Cross-correlation measures similarity between signals as a function of lag:
/// (a ⋆ v)[n] = Σ a[k]v[n+k].
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `ConvolveInput` with:
///   - `a`: First signal array
///   - `v`: Second signal array
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<Vec<f64>, String>` with
/// the cross-correlation result.
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

pub unsafe extern "C" fn rssn_num_signal_cross_correlation_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : ConvolveInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<Vec<f64>, String> {
                    ok : None,
                    err : Some("Invalid Bincode input".to_string()),
                },
            )
        },
    };

    let result =
        signal::cross_correlation(
            &input.a,
            &input.v,
        );

    let ffi_res = FfiResult {
        ok: Some(result),
        err: None::<String>,
    };

    to_bincode_buffer(&ffi_res)
}
