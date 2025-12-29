//! FFI API for numerical special functions.
//!
//! This module provides three FFI interface styles:
//! - `handle`: Direct C-compatible API for direct memory management
//! - `json`: JSON string-based API for language-agnostic integration
//! - `bincode_api`: Binary serialization API for high-performance applications
//!
//! ## Supported Functions
//!
//! ### Gamma and Beta Functions
//! - `gamma_numerical` - Gamma function Γ(x)
//! - `ln_gamma_numerical` - Natural log of gamma ln(Γ(x))
//! - `digamma_numerical` - Digamma function ψ(x)
//! - `beta_numerical` - Beta function B(a, b)
//! - `ln_beta_numerical` - Natural log of beta ln(B(a, b))
//! - `regularized_incomplete_beta` - Regularized incomplete beta Iₓ(a, b)
//! - `regularized_gamma_p` - Lower incomplete gamma P(a, x)
//! - `regularized_gamma_q` - Upper incomplete gamma Q(a, x)
//!
//! ### Error Functions
//! - `erf_numerical` - Error function erf(x)
//! - `erfc_numerical` - Complementary error function erfc(x)
//! - `inverse_erf` - Inverse error function erf⁻¹(x)
//! - `inverse_erfc` - Inverse complementary error function erfc⁻¹(x)
//!
//! ### Combinatorial Functions
//! - `factorial` - Factorial n!
//! - `double_factorial` - Double factorial n!!
//! - `binomial` - Binomial coefficient C(n, k)
//! - `rising_factorial` - Pochhammer symbol (x)ₙ
//! - `falling_factorial` - Falling factorial (x)₍ₙ₎
//! - `ln_factorial` - Natural log of factorial ln(n!)
//!
//! ### Bessel Functions
//! - `bessel_j0`, `bessel_j1` - Bessel functions of the first kind
//! - `bessel_y0`, `bessel_y1` - Bessel functions of the second kind
//! - `bessel_i0`, `bessel_i1` - Modified Bessel functions of the first kind
//! - `bessel_k0`, `bessel_k1` - Modified Bessel functions of the second kind
//!
//! ### Other Functions
//! - `sinc` - Normalized sinc function sin(πx)/(πx)
//! - `zeta` - Riemann zeta function ζ(s)
/// bincode-based FFI bindings for numerical special functions using serialized `Expr` values.
pub mod bincode_api;
/// Handle-based FFI bindings for numerical special functions using opaque `Expr` handles.
pub mod handle;
/// JSON-based FFI bindings for numerical special functions using serialized `Expr` values.
pub mod json;
