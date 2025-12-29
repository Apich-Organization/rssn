//! Handle-based FFI API for numerical special functions.

use crate::numerical::special;

// Gamma functions
/// Computes the Gamma function Γ(x) via handle-based FFI.
///
/// The Gamma function extends the factorial to real and complex numbers: Γ(n) = (n-1)!.
///
/// # Arguments
///
/// * `x` - Argument of the Gamma function (must be positive)
///
/// # Returns
///
/// The value Γ(x).
#[no_mangle]

pub extern "C" fn rssn_num_special_gamma(
    x: f64
) -> f64 {

    special::gamma_numerical(x)
}

/// Computes the natural logarithm of the Gamma function ln(Γ(x)) via handle-based FFI.
///
/// More numerically stable than ln(gamma(x)) for large x.
///
/// # Arguments
///
/// * `x` - Argument of the log-gamma function
///
/// # Returns
///
/// The value ln(Γ(x)).
#[no_mangle]

pub extern "C" fn rssn_num_special_ln_gamma(
    x: f64
) -> f64 {

    special::ln_gamma_numerical(x)
}

/// Computes the Digamma function ψ(x) = d/dx[ln(Γ(x))] via handle-based FFI.
///
/// The logarithmic derivative of the Gamma function.
///
/// # Arguments
///
/// * `x` - Argument of the Digamma function
///
/// # Returns
///
/// The value ψ(x).
#[no_mangle]

pub extern "C" fn rssn_num_special_digamma(
    x: f64
) -> f64 {

    special::digamma_numerical(x)
}

/// Computes the lower incomplete Gamma function γ(s, x) via handle-based FFI.
///
/// Defined as γ(s, x) = ∫₀^x t^(s-1) e^(-t) dt.
///
/// # Arguments
///
/// * `s` - Shape parameter
/// * `x` - Upper limit of integration
///
/// # Returns
///
/// The value γ(s, x).
#[no_mangle]

pub extern "C" fn rssn_num_special_lower_incomplete_gamma(
    s: f64,
    x: f64,
) -> f64 {

    special::lower_incomplete_gamma(
        s, x,
    )
}

/// Computes the upper incomplete Gamma function Γ(s, x) via handle-based FFI.
///
/// Defined as Γ(s, x) = ∫ₓ^∞ t^(s-1) e^(-t) dt.
///
/// # Arguments
///
/// * `s` - Shape parameter
/// * `x` - Lower limit of integration
///
/// # Returns
///
/// The value Γ(s, x).
#[no_mangle]

pub extern "C" fn rssn_num_special_upper_incomplete_gamma(
    s: f64,
    x: f64,
) -> f64 {

    special::upper_incomplete_gamma(
        s, x,
    )
}

// Beta functions
/// Computes the Beta function B(a, b) = Γ(a)Γ(b) / Γ(a+b) via handle-based FFI.
///
/// # Arguments
///
/// * `a` - First shape parameter (must be positive)
/// * `b` - Second shape parameter (must be positive)
///
/// # Returns
///
/// The value B(a, b).
#[no_mangle]

pub extern "C" fn rssn_num_special_beta(
    a: f64,
    b: f64,
) -> f64 {

    special::beta_numerical(a, b)
}

/// Computes the natural logarithm of the Beta function ln(B(a, b)) via handle-based FFI.
///
/// More numerically stable than ln(beta(a, b)) for large parameters.
///
/// # Arguments
///
/// * `a` - First shape parameter (must be positive)
/// * `b` - Second shape parameter (must be positive)
///
/// # Returns
///
/// The value ln(B(a, b)).
#[no_mangle]

pub extern "C" fn rssn_num_special_ln_beta(
    a: f64,
    b: f64,
) -> f64 {

    special::ln_beta_numerical(a, b)
}

/// Computes the regularized incomplete Beta function I(x; a, b) via handle-based FFI.
///
/// The cumulative distribution function of the Beta distribution.
///
/// # Arguments
///
/// * `x` - Upper limit of integration (0 ≤ x ≤ 1)
/// * `a` - First shape parameter (must be positive)
/// * `b` - Second shape parameter (must be positive)
///
/// # Returns
///
/// The value I(x; a, b) ∈ [0, 1].
#[no_mangle]

pub extern "C" fn rssn_num_special_regularized_beta(
    x: f64,
    a: f64,
    b: f64,
) -> f64 {

    special::regularized_beta(x, a, b)
}

// Error functions
/// Computes the error function erf(x) via handle-based FFI.
///
/// Defined as erf(x) = (2/√π) ∫₀^x e^(-t²) dt.
///
/// # Arguments
///
/// * `x` - Argument of the error function
///
/// # Returns
///
/// The value erf(x) ∈ [-1, 1].
#[no_mangle]

pub extern "C" fn rssn_num_special_erf(
    x: f64
) -> f64 {

    special::erf_numerical(x)
}

/// Computes the complementary error function erfc(x) = 1 - erf(x) via handle-based FFI.
///
/// More accurate than 1 - erf(x) for large positive x.
///
/// # Arguments
///
/// * `x` - Argument of the complementary error function
///
/// # Returns
///
/// The value erfc(x) ∈ [0, 2].
#[no_mangle]

pub extern "C" fn rssn_num_special_erfc(
    x: f64
) -> f64 {

    special::erfc_numerical(x)
}

/// Computes the inverse error function erf⁻¹(x) via handle-based FFI.
///
/// Finds y such that erf(y) = x.
///
/// # Arguments
///
/// * `x` - Value in the range (-1, 1)
///
/// # Returns
///
/// The value y such that erf(y) = x.
#[no_mangle]

pub extern "C" fn rssn_num_special_inverse_erf(
    x: f64
) -> f64 {

    special::inverse_erf_numerical(x)
}

// Bessel functions
/// Computes the Bessel function of the first kind of order zero J₀(x) via handle-based FFI.
///
/// Solution to Bessel's differential equation for ν = 0.
///
/// # Arguments
///
/// * `x` - Argument of the Bessel function
///
/// # Returns
///
/// The value J₀(x).
#[no_mangle]

pub extern "C" fn rssn_num_special_bessel_j0(
    x: f64
) -> f64 {

    special::bessel_j0(x)
}

/// Computes the Bessel function of the first kind of order one J₁(x) via handle-based FFI.
///
/// Solution to Bessel's differential equation for ν = 1.
///
/// # Arguments
///
/// * `x` - Argument of the Bessel function
///
/// # Returns
///
/// The value J₁(x).
#[no_mangle]

pub extern "C" fn rssn_num_special_bessel_j1(
    x: f64
) -> f64 {

    special::bessel_j1(x)
}

/// Computes the Bessel function of the second kind of order zero Y₀(x) via handle-based FFI.
///
/// Also known as the Neumann function or Weber function.
///
/// # Arguments
///
/// * `x` - Argument of the Bessel function (must be positive)
///
/// # Returns
///
/// The value Y₀(x).
#[no_mangle]

pub extern "C" fn rssn_num_special_bessel_y0(
    x: f64
) -> f64 {

    special::bessel_y0(x)
}

/// Computes the Bessel function of the second kind of order one Y₁(x) via handle-based FFI.
///
/// Also known as the Neumann function or Weber function.
///
/// # Arguments
///
/// * `x` - Argument of the Bessel function (must be positive)
///
/// # Returns
///
/// The value Y₁(x).
#[no_mangle]

pub extern "C" fn rssn_num_special_bessel_y1(
    x: f64
) -> f64 {

    special::bessel_y1(x)
}

/// Computes the modified Bessel function of the first kind of order zero I₀(x) via handle-based FFI.
///
/// Solution to the modified Bessel equation for ν = 0.
///
/// # Arguments
///
/// * `x` - Argument of the modified Bessel function
///
/// # Returns
///
/// The value I₀(x).
#[no_mangle]

pub extern "C" fn rssn_num_special_bessel_i0(
    x: f64
) -> f64 {

    special::bessel_i0(x)
}

/// Computes the modified Bessel function of the first kind of order one I₁(x) via handle-based FFI.
///
/// Solution to the modified Bessel equation for ν = 1.
///
/// # Arguments
///
/// * `x` - Argument of the modified Bessel function
///
/// # Returns
///
/// The value I₁(x).
#[no_mangle]

pub extern "C" fn rssn_num_special_bessel_i1(
    x: f64
) -> f64 {

    special::bessel_i1(x)
}

// Orthogonal polynomials
/// Computes the Legendre polynomial Pₙ(x) via handle-based FFI.
///
/// Orthogonal on [-1, 1], used in multipole expansions and spherical harmonics.
///
/// # Arguments
///
/// * `n` - Polynomial degree (non-negative integer)
/// * `x` - Evaluation point (typically in [-1, 1])
///
/// # Returns
///
/// The value Pₙ(x).
#[no_mangle]

pub extern "C" fn rssn_num_special_legendre_p(
    n: u32,
    x: f64,
) -> f64 {

    special::legendre_p(n, x)
}

/// Computes the Chebyshev polynomial of the first kind Tₙ(x) via handle-based FFI.
///
/// Satisfies Tₙ(cos θ) = cos(nθ), minimizes polynomial interpolation error.
///
/// # Arguments
///
/// * `n` - Polynomial degree (non-negative integer)
/// * `x` - Evaluation point (typically in [-1, 1])
///
/// # Returns
///
/// The value Tₙ(x).
#[no_mangle]

pub extern "C" fn rssn_num_special_chebyshev_t(
    n: u32,
    x: f64,
) -> f64 {

    special::chebyshev_t(n, x)
}

/// Computes the Chebyshev polynomial of the second kind Uₙ(x) via handle-based FFI.
///
/// Satisfies Uₙ(cos θ) = sin((n+1)θ) / sin(θ).
///
/// # Arguments
///
/// * `n` - Polynomial degree (non-negative integer)
/// * `x` - Evaluation point (typically in [-1, 1])
///
/// # Returns
///
/// The value Uₙ(x).
#[no_mangle]

pub extern "C" fn rssn_num_special_chebyshev_u(
    n: u32,
    x: f64,
) -> f64 {

    special::chebyshev_u(n, x)
}

/// Computes the Hermite polynomial Hₙ(x) via handle-based FFI.
///
/// Orthogonal with weight e^(-x²), appears in quantum harmonic oscillator wavefunctions.
///
/// # Arguments
///
/// * `n` - Polynomial degree (non-negative integer)
/// * `x` - Evaluation point
///
/// # Returns
///
/// The value Hₙ(x).
#[no_mangle]

pub extern "C" fn rssn_num_special_hermite_h(
    n: u32,
    x: f64,
) -> f64 {

    special::hermite_h(n, x)
}

/// Computes the Laguerre polynomial Lₙ(x) via handle-based FFI.
///
/// Orthogonal with weight e^(-x) on [0, ∞), appears in quantum mechanics (hydrogen atom).
///
/// # Arguments
///
/// * `n` - Polynomial degree (non-negative integer)
/// * `x` - Evaluation point (typically ≥ 0)
///
/// # Returns
///
/// The value Lₙ(x).
#[no_mangle]

pub extern "C" fn rssn_num_special_laguerre_l(
    n: u32,
    x: f64,
) -> f64 {

    special::laguerre_l(n, x)
}

// Other special functions
/// Computes the factorial n! via handle-based FFI.
///
/// For large n, computed using the Gamma function.
///
/// # Arguments
///
/// * `n` - Non-negative integer
///
/// # Returns
///
/// The value n! (as f64).
#[no_mangle]

pub extern "C" fn rssn_num_special_factorial(
    n: u64
) -> f64 {

    special::factorial(n)
}

/// Computes the double factorial n!! via handle-based FFI.
///
/// Defined as n!! = n × (n-2) × (n-4) × ... (down to 1 or 2).
///
/// # Arguments
///
/// * `n` - Non-negative integer
///
/// # Returns
///
/// The value n!! (as f64).
#[no_mangle]

pub extern "C" fn rssn_num_special_double_factorial(
    n: u64
) -> f64 {

    special::double_factorial(n)
}

/// Computes the binomial coefficient C(n, k) via handle-based FFI.
///
/// The number of ways to choose k items from n items.
///
/// # Arguments
///
/// * `n` - Total number of items (non-negative integer)
/// * `k` - Number of items to choose (0 ≤ k ≤ n)
///
/// # Returns
///
/// The value C(n, k) (as f64).
#[no_mangle]

pub extern "C" fn rssn_num_special_binomial(
    n: u64,
    k: u64,
) -> f64 {

    special::binomial(n, k)
}

/// Computes the Riemann zeta function ζ(s) via handle-based FFI.
///
/// Defined as ζ(s) = ∑ₙ₌₁^∞ 1/n^s for Re(s) > 1.
///
/// # Arguments
///
/// * `s` - Argument of the zeta function
///
/// # Returns
///
/// The value ζ(s).
#[no_mangle]

pub extern "C" fn rssn_num_special_zeta(
    s: f64
) -> f64 {

    special::riemann_zeta(s)
}

/// Computes the normalized sinc function sinc(x) = sin(x) / x via handle-based FFI.
///
/// With sinc(0) = 1 by continuity. Appears in signal processing and Fourier analysis.
///
/// # Arguments
///
/// * `x` - Input value (radians)
///
/// # Returns
///
/// The value sinc(x).
#[no_mangle]

pub extern "C" fn rssn_num_special_sinc(
    x: f64
) -> f64 {

    special::sinc(x)
}

/// Computes the sigmoid function σ(x) = 1 / (1 + e^(-x)) via handle-based FFI.
///
/// Common activation function in neural networks, maps ℝ to (0, 1).
///
/// # Arguments
///
/// * `x` - Input value
///
/// # Returns
///
/// The value σ(x) ∈ (0, 1).
#[no_mangle]

pub extern "C" fn rssn_num_special_sigmoid(
    x: f64
) -> f64 {

    special::sigmoid(x)
}

/// Computes the softplus function softplus(x) = ln(1 + e^x) via handle-based FFI.
///
/// A smooth approximation to the `ReLU` activation function.
///
/// # Arguments
///
/// * `x` - Input value
///
/// # Returns
///
/// The value softplus(x) > 0.
#[no_mangle]

pub extern "C" fn rssn_num_special_softplus(
    x: f64
) -> f64 {

    special::softplus(x)
}

/// Computes the logit function logit(p) = ln(p / (1 - p)) via handle-based FFI.
///
/// The inverse of the sigmoid function, maps (0, 1) to ℝ.
///
/// # Arguments
///
/// * `p` - Probability value in the range (0, 1)
///
/// # Returns
///
/// The value logit(p).
#[no_mangle]

pub extern "C" fn rssn_num_special_logit(
    p: f64
) -> f64 {

    special::logit(p)
}
