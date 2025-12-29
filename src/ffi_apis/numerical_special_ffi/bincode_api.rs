//! Bincode-based FFI API for numerical special functions.

use serde::Deserialize;

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
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
/// Computes the Gamma function Γ(x) using bincode serialization.
///
/// The Gamma function is defined as Γ(x) = ∫₀^∞ t^(x-1) e^(-t) dt for x > 0,
/// generalizing the factorial function to real and complex numbers: Γ(n) = (n-1)!.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `SingleInput` with:
///   - `x`: Argument of the Gamma function (must be positive)
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<f64, String>` with either:
/// - `ok`: The value Γ(x)
/// - `err`: Error message if input is invalid
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_gamma_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : SingleInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    to_bincode_buffer(&FfiResult {
        ok: Some(
            special::gamma_numerical(
                input.x,
            ),
        ),
        err: None::<String>,
    })
}

/// Computes the natural logarithm of the Gamma function ln(Γ(x)) using bincode serialization.
///
/// This function is more numerically stable than computing ln(Γ(x)) = ln(gamma(x)),
/// especially for large values of x where Γ(x) would overflow.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `SingleInput` with:
///   - `x`: Argument of the log-gamma function
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<f64, String>` with either:
/// - `ok`: The value ln(Γ(x))
/// - `err`: Error message if input is invalid
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_ln_gamma_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : SingleInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    to_bincode_buffer(&FfiResult {
        ok: Some(
            special::ln_gamma_numerical(
                input.x,
            ),
        ),
        err: None::<String>,
    })
}

/// Computes the Digamma function ψ(x) using bincode serialization.
///
/// The Digamma function is the logarithmic derivative of the Gamma function:
/// ψ(x) = d/dx[ln(Γ(x))] = Γ'(x) / Γ(x).
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `SingleInput` with:
///   - `x`: Argument of the Digamma function
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<f64, String>` with either:
/// - `ok`: The value ψ(x)
/// - `err`: Error message if input is invalid
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_digamma_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : SingleInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    to_bincode_buffer(&FfiResult {
        ok: Some(
            special::digamma_numerical(
                input.x,
            ),
        ),
        err: None::<String>,
    })
}

// Beta functions
/// Computes the Beta function B(a, b) using bincode serialization.
///
/// The Beta function is defined as B(a, b) = ∫₀^1 t^(a-1) (1-t)^(b-1) dt,
/// which can be expressed in terms of Gamma functions: B(a, b) = Γ(a)Γ(b) / Γ(a+b).
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `TwoInput` with:
///   - `a`: First shape parameter (must be positive)
///   - `b`: Second shape parameter (must be positive)
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<f64, String>` with either:
/// - `ok`: The value B(a, b)
/// - `err`: Error message if input is invalid
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_beta_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : TwoInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    to_bincode_buffer(&FfiResult {
        ok: Some(
            special::beta_numerical(
                input.a,
                input.b,
            ),
        ),
        err: None::<String>,
    })
}

/// Computes the regularized incomplete Beta function I(x; a, b) using bincode serialization.
///
/// The regularized incomplete Beta function is defined as:
/// I(x; a, b) = B(x; a, b) / B(a, b), where B(x; a, b) = ∫₀^x t^(a-1) (1-t)^(b-1) dt.
/// This function is the cumulative distribution function of the Beta distribution.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `ThreeInput` with:
///   - `x`: Upper limit of integration (0 ≤ x ≤ 1)
///   - `a`: First shape parameter (must be positive)
///   - `b`: Second shape parameter (must be positive)
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<f64, String>` with either:
/// - `ok`: The value I(x; a, b)
/// - `err`: Error message if input is invalid
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_regularized_beta_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : ThreeInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    to_bincode_buffer(&FfiResult {
        ok: Some(
            special::regularized_beta(
                input.x,
                input.a,
                input.b,
            ),
        ),
        err: None::<String>,
    })
}

// Error functions
/// Computes the error function erf(x) using bincode serialization.
///
/// The error function is defined as erf(x) = (2/√π) ∫₀^x e^(-t²) dt.
/// It represents the probability that a random variable from a standard normal distribution
/// falls within [-x, x].
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `SingleInput` with:
///   - `x`: Argument of the error function
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<f64, String>` with either:
/// - `ok`: The value erf(x) (range: [-1, 1])
/// - `err`: Error message if input is invalid
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_erf_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : SingleInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    to_bincode_buffer(&FfiResult {
        ok: Some(
            special::erf_numerical(
                input.x,
            ),
        ),
        err: None::<String>,
    })
}

/// Computes the complementary error function erfc(x) using bincode serialization.
///
/// The complementary error function is defined as erfc(x) = 1 - erf(x) = (2/√π) ∫ₓ^∞ e^(-t²) dt.
/// It is more accurate than computing 1 - erf(x) for large positive x.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `SingleInput` with:
///   - `x`: Argument of the complementary error function
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<f64, String>` with either:
/// - `ok`: The value erfc(x) (range: [0, 2])
/// - `err`: Error message if input is invalid
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_erfc_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : SingleInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    to_bincode_buffer(&FfiResult {
        ok: Some(
            special::erfc_numerical(
                input.x,
            ),
        ),
        err: None::<String>,
    })
}

// Bessel functions
/// Computes the Bessel function of the first kind of order zero J₀(x) using bincode serialization.
///
/// The Bessel function J₀(x) is a solution to Bessel's differential equation for ν = 0:
/// x²y'' + xy' + x²y = 0. It appears in problems with cylindrical symmetry.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `SingleInput` with:
///   - `x`: Argument of the Bessel function
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<f64, String>` with either:
/// - `ok`: The value J₀(x)
/// - `err`: Error message if input is invalid
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_bessel_j0_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : SingleInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    to_bincode_buffer(&FfiResult {
        ok: Some(special::bessel_j0(
            input.x,
        )),
        err: None::<String>,
    })
}

/// Computes the Bessel function of the first kind of order one J₁(x) using bincode serialization.
///
/// The Bessel function J₁(x) is a solution to Bessel's differential equation for ν = 1:
/// x²y'' + xy' + (x² - 1)y = 0. It appears in wave propagation and vibration problems.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `SingleInput` with:
///   - `x`: Argument of the Bessel function
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<f64, String>` with either:
/// - `ok`: The value J₁(x)
/// - `err`: Error message if input is invalid
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_bessel_j1_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : SingleInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    to_bincode_buffer(&FfiResult {
        ok: Some(special::bessel_j1(
            input.x,
        )),
        err: None::<String>,
    })
}

// Orthogonal polynomials
/// Computes the Legendre polynomial Pₙ(x) using bincode serialization.
///
/// The Legendre polynomials are solutions to Legendre's differential equation:
/// (1 - x²)y'' - 2xy' + n(n+1)y = 0. They are orthogonal on [-1, 1] and appear
/// in multipole expansions and spherical harmonics.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `PolyInput` with:
///   - `n`: Polynomial degree (non-negative integer)
///   - `x`: Evaluation point (typically in [-1, 1])
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<f64, String>` with either:
/// - `ok`: The value Pₙ(x)
/// - `err`: Error message if input is invalid
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_legendre_p_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : PolyInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    to_bincode_buffer(&FfiResult {
        ok: Some(special::legendre_p(
            input.n,
            input.x,
        )),
        err: None::<String>,
    })
}

/// Computes the Chebyshev polynomial of the first kind Tₙ(x) using bincode serialization.
///
/// The Chebyshev polynomials satisfy Tₙ(cos θ) = cos(nθ) and are orthogonal on [-1, 1]
/// with weight function 1/√(1-x²). They minimize polynomial interpolation error.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `PolyInput` with:
///   - `n`: Polynomial degree (non-negative integer)
///   - `x`: Evaluation point (typically in [-1, 1])
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<f64, String>` with either:
/// - `ok`: The value Tₙ(x)
/// - `err`: Error message if input is invalid
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_chebyshev_t_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : PolyInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    to_bincode_buffer(&FfiResult {
        ok: Some(
            special::chebyshev_t(
                input.n,
                input.x,
            ),
        ),
        err: None::<String>,
    })
}

/// Computes the Hermite polynomial Hₙ(x) using bincode serialization.
///
/// The Hermite polynomials (physicists' version) are orthogonal on (-∞, ∞) with weight
/// function e^(-x²). They appear in quantum harmonic oscillator wavefunctions and
/// probability theory.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `PolyInput` with:
///   - `n`: Polynomial degree (non-negative integer)
///   - `x`: Evaluation point
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<f64, String>` with either:
/// - `ok`: The value Hₙ(x)
/// - `err`: Error message if input is invalid
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_hermite_h_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : PolyInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    to_bincode_buffer(&FfiResult {
        ok: Some(special::hermite_h(
            input.n,
            input.x,
        )),
        err: None::<String>,
    })
}

// Other special functions
/// Computes the factorial n! using bincode serialization.
///
/// The factorial is defined as n! = n × (n-1) × ... × 2 × 1 for n ≥ 1, with 0! = 1.
/// For large n, the result is computed using the Gamma function to avoid overflow.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `IntInput` with:
///   - `n`: Non-negative integer
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<f64, String>` with either:
/// - `ok`: The value n! (as f64)
/// - `err`: Error message if input is invalid
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_factorial_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : IntInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    to_bincode_buffer(&FfiResult {
        ok: Some(special::factorial(
            input.n,
        )),
        err: None::<String>,
    })
}

/// Computes the binomial coefficient C(n, k) = n! / (k!(n-k)!) using bincode serialization.
///
/// The binomial coefficient represents the number of ways to choose k items from n items
/// without regard to order. Also appears as coefficients in the binomial expansion.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `BinomialInput` with:
///   - `n`: Total number of items (non-negative integer)
///   - `k`: Number of items to choose (0 ≤ k ≤ n)
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<f64, String>` with either:
/// - `ok`: The value C(n, k) (as f64)
/// - `err`: Error message if input is invalid
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_binomial_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : BinomialInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    to_bincode_buffer(&FfiResult {
        ok: Some(special::binomial(
            input.n,
            input.k,
        )),
        err: None::<String>,
    })
}

/// Computes the sigmoid function σ(x) = 1 / (1 + e^(-x)) using bincode serialization.
///
/// The sigmoid is a smooth, S-shaped activation function commonly used in neural networks
/// and logistic regression. It maps any real value to the range (0, 1).
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `SingleInput` with:
///   - `x`: Input value
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<f64, String>` with either:
/// - `ok`: The value σ(x) ∈ (0, 1)
/// - `err`: Error message if input is invalid
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_sigmoid_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : SingleInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    to_bincode_buffer(&FfiResult {
        ok: Some(special::sigmoid(
            input.x,
        )),
        err: None::<String>,
    })
}

/// Computes the normalized sinc function sinc(x) = sin(x) / x using bincode serialization.
///
/// The sinc function is defined as sinc(x) = sin(x)/x for x ≠ 0, and sinc(0) = 1 by continuity.
/// It appears in signal processing (Fourier analysis) and optics (diffraction patterns).
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `SingleInput` with:
///   - `x`: Input value (radians)
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<f64, String>` with either:
/// - `ok`: The value sinc(x)
/// - `err`: Error message if input is invalid
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_sinc_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : SingleInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    to_bincode_buffer(&FfiResult {
        ok: Some(special::sinc(
            input.x,
        )),
        err: None::<String>,
    })
}
