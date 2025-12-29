//! # FFI API for Symbolic Special Functions
//!
//! This module provides Foreign Function Interface (FFI) bindings for the symbolic
//! special functions module. Three API styles are available:
//!
//! ## Supported Functions
//!
//! ### Gamma and Related Functions
//! - `gamma(z)` - Gamma function Γ(z)
//! - `beta(a, b)` - Beta function B(a, b)
//! - `digamma(z)` - Digamma function ψ(z)
//! - `polygamma(n, z)` - Polygamma function ψ⁽ⁿ⁾(z)
//! - `ln_gamma(z)` - Log-gamma function ln(Γ(z))
//!
//! ### Error Functions
//! - `erf(z)` - Error function
//! - `erfc(z)` - Complementary error function
//! - `erfi(z)` - Imaginary error function
//!
//! ### Number Theory
//! - `zeta(s)` - Riemann zeta function ζ(s)
//!
//! ### Bessel Functions
//! - `bessel_j(n, x)` - Bessel function of first kind Jₙ(x)
//! - `bessel_y(n, x)` - Bessel function of second kind Yₙ(x)
//! - `bessel_i(n, x)` - Modified Bessel function of first kind Iₙ(x)
//! - `bessel_k(n, x)` - Modified Bessel function of second kind Kₙ(x)
//!
//! ### Orthogonal Polynomials
//! - `legendre_p(n, x)` - Legendre polynomial Pₙ(x)
//! - `laguerre_l(n, x)` - Laguerre polynomial Lₙ(x)
//! - `generalized_laguerre(n, α, x)` - Generalized Laguerre Lₙᵅ(x)
//! - `hermite_h(n, x)` - Hermite polynomial Hₙ(x)
//! - `chebyshev_t(n, x)` - Chebyshev polynomial of first kind Tₙ(x)
//! - `chebyshev_u(n, x)` - Chebyshev polynomial of second kind Uₙ(x)
//!
//! ### Differential Equations
//! - `bessel_differential_equation(y, x, n)` - x²y'' + xy' + (x² - n²)y = 0
//! - `legendre_differential_equation(y, x, n)` - (1-x²)y'' - 2xy' + n(n+1)y = 0
//! - `laguerre_differential_equation(y, x, n)` - xy'' + (1-x)y' + ny = 0
//! - `hermite_differential_equation(y, x, n)` - y'' - 2xy' + 2ny = 0
//! - `chebyshev_differential_equation(y, x, n)` - (1-x²)y'' - xy' + n²y = 0
//!
//! ### Rodrigues Formulas
//! - `legendre_rodrigues_formula(n, x)`
//! - `hermite_rodrigues_formula(n, x)`
//!
//! ## API Styles
//!
//! - **Handle API** (`handle.rs`): Raw pointer-based API for maximum performance
//! - **JSON API** (`json.rs`): JSON string-based API for language interop
//! - **Bincode API** (`bincode_api.rs`): Binary serialization for efficiency

//! - **Bincode API** (`bincode_api.rs`): Binary serialization for efficiency

/// Bincode API for symbolic special functions.
pub mod bincode_api;
/// Handle API for symbolic special functions.
pub mod handle;
/// JSON API for symbolic special functions.
pub mod json;
