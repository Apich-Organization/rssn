//! JSON-based FFI API for numerical special functions.

use std::os::raw::c_char;

use serde::Deserialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::special;

#[derive(Deserialize)]

struct SingleInput {
    x: f64,
}

#[derive(Deserialize)]

struct TwoInput {
    a: f64,
    b: f64,
}

#[derive(Deserialize)]

struct ThreeInput {
    x: f64,
    a: f64,
    b: f64,
}

#[derive(Deserialize)]

struct PolyInput {
    n: u32,
    x: f64,
}

#[derive(Deserialize)]

struct IntInput {
    n: u64,
}

#[derive(Deserialize)]

struct BinomialInput {
    n: u64,
    k: u64,
}

// Gamma functions
/// Computes the Gamma function Γ(x) via JSON serialization.
///
/// The Gamma function extends the factorial to real and complex numbers: Γ(n) = (n-1)!.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `x`: Argument of the Gamma function (must be positive)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the value Γ(x).
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_special_gamma_json(
    input: *const c_char
) -> *mut c_char {

    let input : SingleInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(special::gamma_numerical(input.x)),
            err : None::<String>,
        })
        .unwrap(),
    )
}

/// Computes the natural logarithm of the Gamma function ln(Γ(x)) via JSON serialization.
///
/// More numerically stable than ln(gamma(x)) for large x.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `x`: Argument of the log-gamma function
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the value ln(Γ(x)).
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_special_ln_gamma_json(
    input: *const c_char
) -> *mut c_char {

    let input : SingleInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(special::ln_gamma_numerical(input.x)),
            err : None::<String>,
        })
        .unwrap(),
    )
}

/// Computes the Digamma function ψ(x) = d/dx[ln(Γ(x))] via JSON serialization.
///
/// The logarithmic derivative of the Gamma function.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `x`: Argument of the Digamma function
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the value ψ(x).
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_special_digamma_json(
    input: *const c_char
) -> *mut c_char {

    let input : SingleInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(special::digamma_numerical(input.x)),
            err : None::<String>,
        })
        .unwrap(),
    )
}

// Beta functions
/// Computes the Beta function B(a, b) = Γ(a)Γ(b) / Γ(a+b) via JSON serialization.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `a`: First shape parameter (must be positive)
///   - `b`: Second shape parameter (must be positive)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the value B(a, b).
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_special_beta_json(
    input: *const c_char
) -> *mut c_char {

    let input : TwoInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(special::beta_numerical(input.a, input.b)),
            err : None::<String>,
        })
        .unwrap(),
    )
}

/// Computes the regularized incomplete Beta function I(x; a, b) via JSON serialization.
///
/// The cumulative distribution function of the Beta distribution.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `x`: Upper limit of integration (0 ≤ x ≤ 1)
///   - `a`: First shape parameter (must be positive)
///   - `b`: Second shape parameter (must be positive)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the value I(x; a, b) ∈ [0, 1].
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_special_regularized_beta_json(
    input: *const c_char
) -> *mut c_char {

    let input : ThreeInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(
                special::regularized_beta(
                    input.x,
                    input.a,
                    input.b,
                ),
            ),
            err : None::<String>,
        })
        .unwrap(),
    )
}

// Error functions
/// Computes the error function erf(x) via JSON serialization.
///
/// Defined as erf(x) = (2/√π) ∫₀^x e^(-t²) dt.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `x`: Argument of the error function
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the value erf(x) ∈ [-1, 1].
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_special_erf_json(
    input: *const c_char
) -> *mut c_char {

    let input : SingleInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(special::erf_numerical(input.x)),
            err : None::<String>,
        })
        .unwrap(),
    )
}

/// Computes the complementary error function erfc(x) = 1 - erf(x) via JSON serialization.
///
/// More accurate than 1 - erf(x) for large positive x.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `x`: Argument of the complementary error function
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the value erfc(x) ∈ [0, 2].
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_special_erfc_json(
    input: *const c_char
) -> *mut c_char {

    let input : SingleInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(special::erfc_numerical(input.x)),
            err : None::<String>,
        })
        .unwrap(),
    )
}

// Bessel functions
/// Computes the Bessel function of the first kind of order zero J₀(x) via JSON serialization.
///
/// Solution to Bessel's differential equation for ν = 0.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `x`: Argument of the Bessel function
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the value J₀(x).
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_special_bessel_j0_json(
    input: *const c_char
) -> *mut c_char {

    let input : SingleInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(
                    special::bessel_j0(
                        input.x,
                    ),
                ),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Computes the Bessel function of the first kind of order one J₁(x) via JSON serialization.
///
/// Solution to Bessel's differential equation for ν = 1.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `x`: Argument of the Bessel function
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the value J₁(x).
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_special_bessel_j1_json(
    input: *const c_char
) -> *mut c_char {

    let input : SingleInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(
                    special::bessel_j1(
                        input.x,
                    ),
                ),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

// Orthogonal polynomials
/// Computes the Legendre polynomial Pₙ(x) via JSON serialization.
///
/// Orthogonal on [-1, 1], used in multipole expansions and spherical harmonics.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `n`: Polynomial degree (non-negative integer)
///   - `x`: Evaluation point (typically in [-1, 1])
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the value Pₙ(x).
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_special_legendre_p_json(
    input: *const c_char
) -> *mut c_char {

    let input : PolyInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(
                    special::legendre_p(
                        input.n,
                        input.x,
                    ),
                ),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Computes the Chebyshev polynomial of the first kind Tₙ(x) via JSON serialization.
///
/// Satisfies Tₙ(cos θ) = cos(nθ), minimizes polynomial interpolation error.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `n`: Polynomial degree (non-negative integer)
///   - `x`: Evaluation point (typically in [-1, 1])
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the value Tₙ(x).
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_special_chebyshev_t_json(
    input: *const c_char
) -> *mut c_char {

    let input : PolyInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(special::chebyshev_t(input.n, input.x)),
            err : None::<String>,
        })
        .unwrap(),
    )
}

/// Computes the Hermite polynomial Hₙ(x) via JSON serialization.
///
/// Orthogonal with weight e^(-x²), appears in quantum harmonic oscillator wavefunctions.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `n`: Polynomial degree (non-negative integer)
///   - `x`: Evaluation point
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the value Hₙ(x).
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_special_hermite_h_json(
    input: *const c_char
) -> *mut c_char {

    let input : PolyInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(
                    special::hermite_h(
                        input.n,
                        input.x,
                    ),
                ),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

// Other special functions
/// Computes the factorial n! via JSON serialization.
///
/// For large n, computed using the Gamma function.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `n`: Non-negative integer
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the value n! (as f64).
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_special_factorial_json(
    input: *const c_char
) -> *mut c_char {

    let input : IntInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(
                    special::factorial(
                        input.n,
                    ),
                ),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Computes the binomial coefficient C(n, k) via JSON serialization.
///
/// The number of ways to choose k items from n items.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `n`: Total number of items (non-negative integer)
///   - `k`: Number of items to choose (0 ≤ k ≤ n)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the value C(n, k) (as f64).
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_special_binomial_json(
    input: *const c_char
) -> *mut c_char {

    let input : BinomialInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(
                    special::binomial(
                        input.n,
                        input.k,
                    ),
                ),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Computes the sigmoid function σ(x) = 1 / (1 + e^(-x)) via JSON serialization.
///
/// Common activation function in neural networks, maps ℝ to (0, 1).
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `x`: Input value
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the value σ(x) ∈ (0, 1).
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_special_sigmoid_json(
    input: *const c_char
) -> *mut c_char {

    let input : SingleInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(
                    special::sigmoid(
                        input.x,
                    ),
                ),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Computes the normalized sinc function sinc(x) = sin(x) / x via JSON serialization.
///
/// With sinc(0) = 1 by continuity. Appears in signal processing and Fourier analysis.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `x`: Input value (radians)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the value sinc(x).
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_special_sinc_json(
    input: *const c_char
) -> *mut c_char {

    let input : SingleInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(
                    special::sinc(
                        input.x,
                    ),
                ),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}
