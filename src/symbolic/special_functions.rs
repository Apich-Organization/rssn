//! # Symbolic Special Functions
//!
//! This module provides symbolic representations and "smart" constructors for various
//! special functions commonly encountered in mathematics, physics, and engineering.
//! These constructors include built-in simplification logic for common arguments and identities.
//!
//! ## Supported Functions
//!
//! ### Gamma and Related Functions
//! - `gamma(z)` - Gamma function Γ(z)
//! - `beta(a, b)` - Beta function B(a, b)
//! - `digamma(z)` - Digamma function ψ(z) = Γ'(z)/Γ(z)
//! - `polygamma(n, z)` - Polygamma function ψ⁽ⁿ⁾(z)
//! - `ln_gamma(z)` - Log-gamma function ln(Γ(z))
//!
//! ### Error Functions
//! - `erf(z)` - Error function
//! - `erfc(z)` - Complementary error function 1 - erf(z)
//! - `erfi(z)` - Imaginary error function -i·erf(iz)
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
//! - `hermite_h(n, x)` - Hermite polynomial Hₙ(x)
//! - `chebyshev_t(n, x)` - Chebyshev polynomial of first kind Tₙ(x)
//! - `chebyshev_u(n, x)` - Chebyshev polynomial of second kind Uₙ(x)
//!
//! ### Differential Equations
//! - `bessel_differential_equation` - x²y'' + xy' + (x² - n²)y = 0
//! - `legendre_differential_equation` - (1-x²)y'' - 2xy' + n(n+1)y = 0
//! - `laguerre_differential_equation` - xy'' + (1-x)y' + ny = 0
//! - `hermite_differential_equation` - y'' - 2xy' + 2ny = 0
//!
//! ### Rodrigues Formulas
//! - `legendre_rodrigues_formula` - Pₙ(x) = (1/2ⁿn!)·dⁿ/dxⁿ[(x²-1)ⁿ]
//! - `hermite_rodrigues_formula` - Hₙ(x) = (-1)ⁿeˣ²·dⁿ/dxⁿ[e⁻ˣ²]
//!
//! ## Examples
//!
//! ### Gamma Function
//! ```
//! 
//! use rssn::symbolic::core::Expr;
//! use rssn::symbolic::special_functions::gamma;
//!
//! // Γ(5) = 4! = 24
//! let g = gamma(Expr::Constant(5.0));
//!
//! assert_eq!(
//!     g,
//!     Expr::Constant(24.0)
//! );
//!
//! // Γ(0.5) = √π
//! let g_half = gamma(Expr::Constant(0.5));
//! // Returns Expr::Sqrt(Expr::Pi)
//! ```
//!
//! ### Orthogonal Polynomials
//! ```
//! 
//! use rssn::symbolic::core::Expr;
//! use rssn::symbolic::special_functions::hermite_h;
//! use rssn::symbolic::special_functions::legendre_p;
//!
//! // P₀(x) = 1
//! let p0 = legendre_p(
//!     Expr::Constant(0.0),
//!     Expr::Variable("x".to_string()),
//! );
//!
//! assert_eq!(
//!     p0,
//!     Expr::Constant(1.0)
//! );
//!
//! // H₀(x) = 1, H₁(x) = 2x
//! let h0 = hermite_h(
//!     Expr::Constant(0.0),
//!     Expr::Variable("x".to_string()),
//! );
//!
//! assert_eq!(
//!     h0,
//!     Expr::Constant(1.0)
//! );
//! ```

use std::sync::Arc;

use crate::symbolic::calculus::differentiate;
use crate::symbolic::calculus::factorial;
use crate::symbolic::core::Expr;
use crate::symbolic::simplify::is_zero;
use crate::symbolic::simplify_dag::simplify;

// ============================================================================
// Gamma and Related Functions
// ============================================================================

/// Symbolic representation and smart constructor for the Gamma function, `Γ(z)`.
///
/// The Gamma function is an extension of the factorial function to complex numbers.
/// This constructor includes simplification rules for integer and half-integer arguments,
/// as well as the recurrence relation `Γ(z+1) = zΓ(z)`.
///
/// # Arguments
/// * `arg` - The argument `z` of the Gamma function.
///
/// # Returns
/// An `Expr` representing `Γ(z)`.
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::core::Expr;
/// use rssn::symbolic::special_functions::gamma;
///
/// // Γ(5) = 4! = 24
/// let g = gamma(Expr::Constant(5.0));
///
/// assert_eq!(
///     g,
///     Expr::Constant(24.0)
/// );
/// ```
#[must_use]

pub fn gamma(arg : Expr) -> Expr {

    let s_arg = simplify(&arg);

    if let Some(n) = s_arg.to_f64() {

        if n > 0.0 && n.fract() == 0.0 {

            return Expr::Constant(factorial(
                (n - 1.0) as usize,
            ));
        }

        if (n - 0.5).abs() < 1e-9 {

            return Expr::new_sqrt(Expr::Pi);
        }
    }

    if s_arg == Expr::Constant(1.0) {

        return Expr::Constant(1.0);
    }

    if let Expr::Add(a, b) = &s_arg {

        if **b == Expr::Constant(1.0) {

            return simplify(&Expr::new_mul(
                a.clone(),
                gamma(a.as_ref().clone()),
            ));
        }

        if **a == Expr::Constant(1.0) {

            return simplify(&Expr::new_mul(
                b.clone(),
                gamma(b.as_ref().clone()),
            ));
        }
    }

    Expr::new_gamma(s_arg)
}

/// Symbolic representation for the log-gamma function, `ln(Γ(z))`.
///
/// The log-gamma function is useful for avoiding overflow when computing
/// the gamma function for large arguments.
///
/// # Arguments
/// * `arg` - The argument `z`.
///
/// # Returns
/// An `Expr` representing `ln(Γ(z))`.
#[must_use]

pub fn ln_gamma(arg : Expr) -> Expr {

    let g = gamma(arg);

    // If gamma simplified to a constant, compute the log
    if let Expr::Constant(v) = &g {

        if *v > 0.0 {

            return Expr::Constant(v.ln());
        }
    }

    Expr::new_log(g)
}

/// Symbolic representation and smart constructor for the Beta function, `B(x, y)`.
///
/// The Beta function is closely related to the Gamma function by the identity:
/// `B(x, y) = Γ(x)Γ(y) / Γ(x+y)`.
///
/// # Arguments
/// * `a` - The first argument `x`.
/// * `b` - The second argument `y`.
///
/// # Returns
/// An `Expr` representing `B(x, y)`.
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::core::Expr;
/// use rssn::symbolic::special_functions::beta;
///
/// // B(1, 1) = 1
/// let b = beta(
///     Expr::Constant(1.0),
///     Expr::Constant(1.0),
/// );
///
/// assert_eq!(
///     b,
///     Expr::Constant(1.0)
/// );
/// ```
#[must_use]

pub fn beta(
    a : Expr,
    b : Expr,
) -> Expr {

    let gamma_a = gamma(a.clone());

    let gamma_b = gamma(b.clone());

    let gamma_a_plus_b = gamma(simplify(
        &Expr::new_add(a, b),
    ));

    simplify(&Expr::new_div(
        Expr::new_mul(gamma_a, gamma_b),
        gamma_a_plus_b,
    ))
}

/// Symbolic representation and smart constructor for the Digamma function, `ψ(z)`.
///
/// The Digamma function is the logarithmic derivative of the Gamma function: `ψ(z) = Γ'(z)/Γ(z)`.
/// This constructor includes simplification rules for `ψ(1)` (related to Euler-Mascheroni constant)
/// and the recurrence relation `ψ(z+1) = ψ(z) + 1/z`.
///
/// # Arguments
/// * `arg` - The argument `z` of the Digamma function.
///
/// # Returns
/// An `Expr` representing `ψ(z)`.
#[must_use]

pub fn digamma(arg : Expr) -> Expr {

    let s_arg = simplify(&arg);

    if let Some(n) = s_arg.to_f64() {

        if (n - 1.0).abs() < 1e-9 {

            return Expr::Variable("-gamma".to_string());
        }
    }

    if let Expr::Add(a, b) = &s_arg {

        if **b == Expr::Constant(1.0) {

            return simplify(&Expr::new_add(
                digamma(a.as_ref().clone()),
                Expr::new_div(
                    Expr::Constant(1.0),
                    a.clone(),
                ),
            ));
        }

        if **a == Expr::Constant(1.0) {

            return simplify(&Expr::new_add(
                digamma(b.as_ref().clone()),
                Expr::new_div(
                    Expr::Constant(1.0),
                    b.clone(),
                ),
            ));
        }
    }

    Expr::new_digamma(s_arg)
}

/// Symbolic representation for the Polygamma function, `ψ⁽ⁿ⁾(z)`.
///
/// The polygamma function is the n-th derivative of the digamma function.
/// `ψ⁽⁰⁾(z) = ψ(z)` (digamma)
/// `ψ⁽¹⁾(z)` is called the trigamma function
/// `ψ⁽²⁾(z)` is called the tetragamma function, etc.
///
/// # Arguments
/// * `n` - The order of the derivative (non-negative integer).
/// * `z` - The argument of the function.
///
/// # Returns
/// An `Expr` representing `ψ⁽ⁿ⁾(z)`.
#[must_use]

pub fn polygamma(
    n : Expr,
    z : Expr,
) -> Expr {

    let s_n = simplify(&n);

    let s_z = simplify(&z);

    // ψ⁽⁰⁾(z) = ψ(z)
    if let Some(order) = s_n.to_f64() {

        if order.abs() < 1e-9 {

            return digamma(s_z);
        }
    }

    // Use BinaryList for polygamma as it's not in the core constructors
    Expr::BinaryList(
        "polygamma".to_string(),
        Arc::new(s_n),
        Arc::new(s_z),
    )
}

// ============================================================================
// Error Functions
// ============================================================================

/// Symbolic representation and smart constructor for the Error Function, `erf(z)`.
///
/// The error function is a special function of sigmoid shape that arises in probability,
/// statistics, and partial differential equations. This constructor includes simplification
/// rules for `erf(0)` and `erf(-z)`.
///
/// # Arguments
/// * `arg` - The argument `z` of the error function.
///
/// # Returns
/// An `Expr` representing `erf(z)`.
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::core::Expr;
/// use rssn::symbolic::special_functions::erf;
///
/// // erf(0) = 0
/// let e = erf(Expr::Constant(0.0));
///
/// assert_eq!(
///     e,
///     Expr::Constant(0.0)
/// );
/// ```
#[must_use]

pub fn erf(arg : Expr) -> Expr {

    let s_arg = simplify(&arg);

    if is_zero(&s_arg) {

        return Expr::Constant(0.0);
    }

    if let Expr::Neg(inner) = s_arg {

        return Expr::new_neg(erf((*inner).clone()));
    }

    Expr::new_erf(s_arg)
}

/// Symbolic representation and smart constructor for the Complementary Error Function, `erfc(z)`.
///
/// Defined as `erfc(z) = 1 - erf(z)`.
///
/// # Arguments
/// * `arg` - The argument `z`.
///
/// # Returns
/// An `Expr` representing `erfc(z)`.
#[must_use]

pub fn erfc(arg : Expr) -> Expr {

    simplify(&Expr::new_sub(
        Expr::Constant(1.0),
        erf(arg),
    ))
}

/// Symbolic representation and smart constructor for the Imaginary Error Function, `erfi(z)`.
///
/// Defined as `erfi(z) = -i * erf(iz)`.
///
/// # Arguments
/// * `arg` - The argument `z`.
///
/// # Returns
/// An `Expr` representing `erfi(z)`.
#[must_use]

pub fn erfi(arg : Expr) -> Expr {

    let i = Expr::new_complex(
        Expr::Constant(0.0),
        Expr::Constant(1.0),
    );

    simplify(&Expr::new_mul(
        Expr::new_neg(i.clone()),
        erf(Expr::new_mul(
            i, arg,
        )),
    ))
}

// ============================================================================
// Riemann Zeta Function
// ============================================================================

/// Symbolic representation and smart constructor for the Riemann Zeta function, `ζ(s)`.
///
/// The Riemann Zeta function is a central object in number theory, with connections
/// to the distribution of prime numbers. This constructor includes simplification
/// rules for specific integer arguments (e.g., `ζ(0)`, `ζ(1)`, `ζ(2)`) and trivial zeros.
///
/// # Arguments
/// * `arg` - The argument `s` of the Zeta function.
///
/// # Returns
/// An `Expr` representing `ζ(s)`.
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::core::Expr;
/// use rssn::symbolic::special_functions::zeta;
///
/// // ζ(0) = -1/2
/// let z0 = zeta(Expr::Constant(0.0));
///
/// assert_eq!(
///     z0,
///     Expr::Constant(-0.5)
/// );
/// ```
#[must_use]

pub fn zeta(arg : Expr) -> Expr {

    let s_arg = simplify(&arg);

    if let Some(n) = s_arg.to_f64() {

        if n.fract() == 0.0 {

            let n_int = n as i32;

            if n_int == 0 {

                return Expr::Constant(-0.5);
            }

            if n_int == 1 {

                return Expr::Infinity;
            }

            if n_int == 2 {

                return simplify(&Expr::new_div(
                    Expr::new_pow(
                        Expr::Pi,
                        Expr::Constant(2.0),
                    ),
                    Expr::Constant(6.0),
                ));
            }

            if n_int == 4 {

                return simplify(&Expr::new_div(
                    Expr::new_pow(
                        Expr::Pi,
                        Expr::Constant(4.0),
                    ),
                    Expr::Constant(90.0),
                ));
            }

            if n_int < 0 && n_int % 2 == 0 {

                return Expr::Constant(0.0);
            }
        }
    }

    Expr::new_zeta(s_arg)
}

// ============================================================================
// Bessel Functions
// ============================================================================

/// Symbolic representation and smart constructor for the Bessel function of the first kind, `J_n(x)`.
///
/// Bessel functions are solutions to Bessel's differential equation and arise in many
/// problems of wave propagation and potential theory. This constructor includes simplification
/// rules for `J_n(0)` and `J_{-n}(x)`.
///
/// # Arguments
/// * `order` - The order `n` of the Bessel function.
/// * `arg` - The argument `x` of the Bessel function.
///
/// # Returns
/// An `Expr` representing `J_n(x)`.
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::core::Expr;
/// use rssn::symbolic::special_functions::bessel_j;
///
/// // J_0(0) = 1
/// let j0 = bessel_j(
///     Expr::Constant(0.0),
///     Expr::Constant(0.0),
/// );
///
/// assert_eq!(
///     j0,
///     Expr::Constant(1.0)
/// );
/// ```
#[must_use]

pub fn bessel_j(
    order : Expr,
    arg : Expr,
) -> Expr {

    let s_order = simplify(&order);

    let s_arg = simplify(&arg);

    if is_zero(&s_arg) {

        if let Some(n) = s_order.to_f64() {

            if n == 0.0 {

                return Expr::Constant(1.0);
            }

            if n > 0.0 {

                return Expr::Constant(0.0);
            }
        }
    }

    if let Expr::Neg(inner_order) = &s_order {

        if let Some(n) = inner_order.to_f64() {

            if n.fract() == 0.0 {

                let factor = Expr::new_pow(
                    Expr::Constant(-1.0),
                    Expr::Constant(n),
                );

                return simplify(&Expr::new_mul(
                    factor,
                    bessel_j(
                        inner_order
                            .as_ref()
                            .clone(),
                        s_arg,
                    ),
                ));
            }
        }
    }

    Expr::new_bessel_j(s_order, s_arg)
}

/// Symbolic representation and smart constructor for the Bessel function of the second kind, `Y_n(x)`.
///
/// Bessel functions of the second kind are also solutions to Bessel's differential equation,
/// linearly independent from `J_n(x)`. They typically have a singularity at the origin.
///
/// # Arguments
/// * `order` - The order `n` of the Bessel function.
/// * `arg` - The argument `x` of the Bessel function.
///
/// # Returns
/// An `Expr` representing `Y_n(x)`.
#[must_use]

pub fn bessel_y(
    order : Expr,
    arg : Expr,
) -> Expr {

    let s_order = simplify(&order);

    let s_arg = simplify(&arg);

    if is_zero(&s_arg) {

        return Expr::NegativeInfinity;
    }

    Expr::new_bessel_y(s_order, s_arg)
}

/// Symbolic representation for the Modified Bessel function of the first kind, `I_n(x)`.
///
/// Modified Bessel functions satisfy the modified Bessel equation and arise in
/// problems with cylindrical symmetry involving non-oscillatory behavior.
///
/// # Arguments
/// * `order` - The order `n` of the Bessel function.
/// * `arg` - The argument `x` of the Bessel function.
///
/// # Returns
/// An `Expr` representing `I_n(x)`.
#[must_use]

pub fn bessel_i(
    order : Expr,
    arg : Expr,
) -> Expr {

    let s_order = simplify(&order);

    let s_arg = simplify(&arg);

    if is_zero(&s_arg) {

        if let Some(n) = s_order.to_f64() {

            if n == 0.0 {

                return Expr::Constant(1.0);
            }

            if n > 0.0 {

                return Expr::Constant(0.0);
            }
        }
    }

    // Use BinaryList for bessel_i as it's not in the core constructors
    Expr::BinaryList(
        "bessel_i".to_string(),
        Arc::new(s_order),
        Arc::new(s_arg),
    )
}

/// Symbolic representation for the Modified Bessel function of the second kind, `K_n(x)`.
///
/// This function is also known as the Macdonald function or Basset function.
///
/// # Arguments
/// * `order` - The order `n` of the Bessel function.
/// * `arg` - The argument `x` of the Bessel function.
///
/// # Returns
/// An `Expr` representing `K_n(x)`.
#[must_use]

pub fn bessel_k(
    order : Expr,
    arg : Expr,
) -> Expr {

    let s_order = simplify(&order);

    let s_arg = simplify(&arg);

    if is_zero(&s_arg) {

        return Expr::Infinity;
    }

    // Use BinaryList for bessel_k as it's not in the core constructors
    Expr::BinaryList(
        "bessel_k".to_string(),
        Arc::new(s_order),
        Arc::new(s_arg),
    )
}

// ============================================================================
// Orthogonal Polynomials
// ============================================================================

/// Symbolic representation and smart constructor for the Legendre Polynomials, `P_n(x)`.
///
/// Legendre polynomials are solutions to Legendre's differential equation and form a complete
/// orthogonal set over the interval `[-1, 1]`. They are widely used in physics and engineering.
/// This constructor includes simplification rules for `P_0(x)`, `P_1(x)`, and the recurrence relation.
///
/// # Arguments
/// * `degree` - The degree `n` of the Legendre polynomial.
/// * `arg` - The argument `x` of the polynomial.
///
/// # Returns
/// An `Expr` representing `P_n(x)`.
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::core::Expr;
/// use rssn::symbolic::special_functions::legendre_p;
///
/// // P_0(x) = 1
/// let p0 = legendre_p(
///     Expr::Constant(0.0),
///     Expr::Variable("x".to_string()),
/// );
///
/// assert_eq!(
///     p0,
///     Expr::Constant(1.0)
/// );
/// ```
#[must_use]

pub fn legendre_p(
    degree : Expr,
    arg : Expr,
) -> Expr {

    let s_degree = simplify(&degree);

    if let Some(n) = s_degree.to_f64() {

        let n_int = n as i32;

        if n >= 0.0 && n.fract() == 0.0 {

            if n_int == 0 {

                return Expr::Constant(1.0);
            }

            if n_int == 1 {

                return arg;
            }

            // Use recurrence for small n only to avoid stack overflow
            if n_int <= 10 {

                let p_n = legendre_p(
                    Expr::Constant(n - 1.0),
                    arg.clone(),
                );

                let p_n_minus_1 = legendre_p(
                    Expr::Constant(n - 2.0),
                    arg.clone(),
                );

                let term1 = Expr::new_mul(
                    Expr::Constant(2.0f64.mul_add(n, -1.0)),
                    Expr::new_mul(arg, p_n),
                );

                let term2 = Expr::new_mul(
                    Expr::Constant(n - 1.0),
                    p_n_minus_1,
                );

                return simplify(&Expr::new_div(
                    Expr::new_sub(term1, term2),
                    Expr::Constant(n),
                ));
            }
        }
    }

    Expr::new_legendre_p(s_degree, arg)
}

/// Symbolic representation and smart constructor for the Laguerre Polynomials, `L_n(x)`.
///
/// Laguerre polynomials are solutions to Laguerre's differential equation and are important
/// in quantum mechanics (e.g., hydrogen atom wave functions) and other areas.
/// This constructor includes simplification rules for `L_0(x)`, `L_1(x)`, and the recurrence relation.
///
/// # Arguments
/// * `degree` - The degree `n` of the Laguerre polynomial.
/// * `arg` - The argument `x` of the polynomial.
///
/// # Returns
/// An `Expr` representing `L_n(x)`.
#[must_use]

pub fn laguerre_l(
    degree : Expr,
    arg : Expr,
) -> Expr {

    let s_degree = simplify(&degree);

    if let Some(n) = s_degree.to_f64() {

        let n_int = n as i32;

        if n >= 0.0 && n.fract() == 0.0 {

            if n_int == 0 {

                return Expr::Constant(1.0);
            }

            if n_int == 1 {

                return simplify(&Expr::new_sub(
                    Expr::Constant(1.0),
                    arg,
                ));
            }

            // Use recurrence for small n only
            if n_int <= 10 {

                let l_n = laguerre_l(
                    Expr::Constant(n - 1.0),
                    arg.clone(),
                );

                let l_n_minus_1 = laguerre_l(
                    Expr::Constant(n - 2.0),
                    arg.clone(),
                );

                let term1_factor = simplify(&Expr::new_sub(
                    Expr::Constant(2.0f64.mul_add(n, -1.0)),
                    arg,
                ));

                let term1 = Expr::new_mul(term1_factor, l_n);

                let term2 = Expr::new_mul(
                    Expr::Constant(n - 1.0),
                    l_n_minus_1,
                );

                return simplify(&Expr::new_div(
                    Expr::new_sub(term1, term2),
                    Expr::Constant(n),
                ));
            }
        }
    }

    Expr::new_laguerre_l(s_degree, arg)
}

/// Symbolic representation for the Generalized Laguerre Polynomials, `L_n^α(x)`.
///
/// Generalized (associated) Laguerre polynomials appear in the radial part of
/// the hydrogen atom wave function.
///
/// # Arguments
/// * `n` - The degree.
/// * `alpha` - The generalization parameter.
/// * `x` - The argument.
///
/// # Returns
/// An `Expr` representing `L_n^α(x)`.
#[must_use]

pub fn generalized_laguerre(
    n : Expr,
    alpha : Expr,
    x : Expr,
) -> Expr {

    let s_n = simplify(&n);

    let s_alpha = simplify(&alpha);

    let s_x = simplify(&x);

    // L_n^0(x) = L_n(x)
    if let Some(a) = s_alpha.to_f64() {

        if a.abs() < 1e-9 {

            return laguerre_l(s_n, s_x);
        }
    }

    // L_0^α(x) = 1
    if let Some(n_val) = s_n.to_f64() {

        if n_val.abs() < 1e-9 {

            return Expr::Constant(1.0);
        }
    }

    // Use NaryList for generalized_laguerre
    Expr::NaryList(
        "generalized_laguerre".to_string(),
        vec![s_n, s_alpha, s_x],
    )
}

/// Symbolic representation and smart constructor for the Hermite Polynomials, `H_n(x)`.
///
/// Hermite polynomials are solutions to Hermite's differential equation and are crucial
/// in quantum mechanics (e.g., harmonic oscillator) and probability theory.
/// This constructor includes simplification rules for `H_0(x)`, `H_1(x)`, and the recurrence relation.
///
/// # Arguments
/// * `degree` - The degree `n` of the Hermite polynomial.
/// * `arg` - The argument `x` of the polynomial.
///
/// # Returns
/// An `Expr` representing `H_n(x)`.
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::core::Expr;
/// use rssn::symbolic::special_functions::hermite_h;
///
/// // H_0(x) = 1
/// let h0 = hermite_h(
///     Expr::Constant(0.0),
///     Expr::Variable("x".to_string()),
/// );
///
/// assert_eq!(
///     h0,
///     Expr::Constant(1.0)
/// );
/// ```
#[must_use]

pub fn hermite_h(
    degree : Expr,
    arg : Expr,
) -> Expr {

    let s_degree = simplify(&degree);

    if let Some(n) = s_degree.to_f64() {

        let n_int = n as i32;

        if n >= 0.0 && n.fract() == 0.0 {

            if n_int == 0 {

                return Expr::Constant(1.0);
            }

            if n_int == 1 {

                return simplify(&Expr::new_mul(
                    Expr::Constant(2.0),
                    arg,
                ));
            }

            // Use recurrence for small n only
            if n_int <= 10 {

                let h_n = hermite_h(
                    Expr::Constant(n - 1.0),
                    arg.clone(),
                );

                let h_n_minus_1 = hermite_h(
                    Expr::Constant(n - 2.0),
                    arg.clone(),
                );

                let term1 = Expr::new_mul(
                    Expr::Constant(2.0),
                    Expr::new_mul(arg, h_n),
                );

                let term2 = Expr::new_mul(
                    Expr::Constant(2.0 * (n - 1.0)),
                    h_n_minus_1,
                );

                return simplify(&Expr::new_sub(
                    term1, term2,
                ));
            }
        }
    }

    Expr::new_hermite_h(s_degree, arg)
}

/// Symbolic representation for the Chebyshev Polynomial of the first kind, `T_n(x)`.
///
/// Chebyshev polynomials are defined by `T_n(cos(θ)) = cos(nθ)` and are important
/// in approximation theory and numerical analysis.
///
/// # Arguments
/// * `n` - The degree.
/// * `x` - The argument.
///
/// # Returns
/// An `Expr` representing `T_n(x)`.
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::core::Expr;
/// use rssn::symbolic::special_functions::chebyshev_t;
///
/// // T_0(x) = 1, T_1(x) = x
/// let t0 = chebyshev_t(
///     Expr::Constant(0.0),
///     Expr::Variable("x".to_string()),
/// );
///
/// assert_eq!(
///     t0,
///     Expr::Constant(1.0)
/// );
/// ```
#[must_use]

pub fn chebyshev_t(
    n : Expr,
    x : Expr,
) -> Expr {

    let s_n = simplify(&n);

    let s_x = simplify(&x);

    if let Some(n_val) = s_n.to_f64() {

        let n_int = n_val as i32;

        if n_val >= 0.0 && n_val.fract() == 0.0 {

            if n_int == 0 {

                return Expr::Constant(1.0);
            }

            if n_int == 1 {

                return s_x;
            }

            // Recurrence: T_n(x) = 2x*T_{n-1}(x) - T_{n-2}(x)
            if n_int <= 10 {

                let t_n1 = chebyshev_t(
                    Expr::Constant(n_val - 1.0),
                    s_x.clone(),
                );

                let t_n2 = chebyshev_t(
                    Expr::Constant(n_val - 2.0),
                    s_x.clone(),
                );

                return simplify(&Expr::new_sub(
                    Expr::new_mul(
                        Expr::Constant(2.0),
                        Expr::new_mul(s_x, t_n1),
                    ),
                    t_n2,
                ));
            }
        }
    }

    // Use BinaryList for chebyshev_t
    Expr::BinaryList(
        "chebyshev_t".to_string(),
        Arc::new(s_n),
        Arc::new(s_x),
    )
}

/// Symbolic representation for the Chebyshev Polynomial of the second kind, `U_n(x)`.
///
/// Chebyshev polynomials of the second kind satisfy `U_n(cos(θ))·sin(θ) = sin((n+1)θ)`.
///
/// # Arguments
/// * `n` - The degree.
/// * `x` - The argument.
///
/// # Returns
/// An `Expr` representing `U_n(x)`.
#[must_use]

pub fn chebyshev_u(
    n : Expr,
    x : Expr,
) -> Expr {

    let s_n = simplify(&n);

    let s_x = simplify(&x);

    if let Some(n_val) = s_n.to_f64() {

        let n_int = n_val as i32;

        if n_val >= 0.0 && n_val.fract() == 0.0 {

            if n_int == 0 {

                return Expr::Constant(1.0);
            }

            if n_int == 1 {

                return simplify(&Expr::new_mul(
                    Expr::Constant(2.0),
                    s_x,
                ));
            }

            // Recurrence: U_n(x) = 2x*U_{n-1}(x) - U_{n-2}(x)
            if n_int <= 10 {

                let u_n1 = chebyshev_u(
                    Expr::Constant(n_val - 1.0),
                    s_x.clone(),
                );

                let u_n2 = chebyshev_u(
                    Expr::Constant(n_val - 2.0),
                    s_x.clone(),
                );

                return simplify(&Expr::new_sub(
                    Expr::new_mul(
                        Expr::Constant(2.0),
                        Expr::new_mul(s_x, u_n1),
                    ),
                    u_n2,
                ));
            }
        }
    }

    // Use BinaryList for chebyshev_u
    Expr::BinaryList(
        "chebyshev_u".to_string(),
        Arc::new(s_n),
        Arc::new(s_x),
    )
}

// ============================================================================
// Differential Equations
// ============================================================================

/// Represents Bessel's differential equation: `x²y'' + xy' + (x² - n²)y = 0`.
///
/// This second-order linear ordinary differential equation is fundamental in many areas
/// of physics and engineering, particularly in problems involving cylindrical symmetry.
/// Its solutions are known as Bessel functions.
///
/// # Arguments
/// * `y` - The unknown function `y(x)`.
/// * `x` - The independent variable `x`.
/// * `n` - The order `n` of the Bessel equation.
///
/// # Returns
/// An `Expr::Eq` representing Bessel's differential equation.
#[must_use]

pub fn bessel_differential_equation(
    y : &Expr,
    x : &Expr,
    n : &Expr,
) -> Expr {

    let y_prime = differentiate(y, "x");

    let y_double_prime = differentiate(&y_prime, "x");

    let term1 = Expr::new_mul(
        Expr::new_pow(
            x.clone(),
            Expr::Constant(2.0),
        ),
        y_double_prime,
    );

    let term2 = Expr::new_mul(x.clone(), y_prime);

    let term3 = Expr::new_mul(
        Expr::new_sub(
            Expr::new_pow(
                x.clone(),
                Expr::Constant(2.0),
            ),
            Expr::new_pow(
                n.clone(),
                Expr::Constant(2.0),
            ),
        ),
        y.clone(),
    );

    Expr::Eq(
        Arc::new(Expr::new_add(
            term1,
            Expr::new_add(term2, term3),
        )),
        Arc::new(Expr::Constant(0.0)),
    )
}

/// Represents Legendre's differential equation: `(1-x²)y'' - 2xy' + n(n+1)y = 0`.
///
/// This second-order linear ordinary differential equation arises in the solution of
/// Laplace's equation in spherical coordinates. Its solutions are Legendre polynomials.
///
/// # Arguments
/// * `y` - The unknown function `y(x)`.
/// * `x` - The independent variable `x`.
/// * `n` - The degree `n` of the Legendre equation.
///
/// # Returns
/// An `Expr::Eq` representing Legendre's differential equation.
#[must_use]

pub fn legendre_differential_equation(
    y : &Expr,
    x : &Expr,
    n : &Expr,
) -> Expr {

    let y_prime = differentiate(y, "x");

    let y_double_prime = differentiate(&y_prime, "x");

    let term1 = Expr::new_mul(
        Expr::new_sub(
            Expr::Constant(1.0),
            Expr::new_pow(
                x.clone(),
                Expr::Constant(2.0),
            ),
        ),
        y_double_prime,
    );

    let term2 = Expr::new_mul(
        Expr::Constant(-2.0),
        Expr::new_mul(x.clone(), y_prime),
    );

    let term3 = Expr::new_mul(
        Expr::new_mul(
            n.clone(),
            Expr::new_add(
                n.clone(),
                Expr::Constant(1.0),
            ),
        ),
        y.clone(),
    );

    Expr::Eq(
        Arc::new(Expr::new_add(
            term1,
            Expr::new_add(term2, term3),
        )),
        Arc::new(Expr::Constant(0.0)),
    )
}

/// Represents Rodrigues' Formula for Legendre Polynomials: `P_n(x) = (1/(2^n n!)) * d^n/dx^n [(x^2 - 1)^n]`.
///
/// This formula provides a compact way to define Legendre polynomials and is useful
/// for deriving their properties.
///
/// # Arguments
/// * `n` - The degree `n` of the Legendre polynomial.
/// * `x` - The independent variable `x`.
///
/// # Returns
/// An `Expr::Eq` representing Rodrigues' formula.
#[must_use]

pub fn legendre_rodrigues_formula(
    n : &Expr,
    x : &Expr,
) -> Expr {

    let n_f64 = if let Expr::Constant(val) = n {

        *val
    } else {

        return Expr::new_legendre_p(n.clone(), x.clone());
    };

    let n_factorial = Expr::Constant(factorial(
        n_f64 as usize,
    ));

    Expr::Eq(
        Arc::new(legendre_p(
            n.clone(),
            x.clone(),
        )),
        Arc::new(Expr::new_mul(
            Expr::new_div(
                Expr::Constant(1.0),
                Expr::new_mul(
                    Expr::new_pow(
                        Expr::Constant(2.0),
                        n.clone(),
                    ),
                    n_factorial,
                ),
            ),
            Expr::DerivativeN(
                Arc::new(Expr::new_pow(
                    Expr::new_sub(
                        Expr::new_pow(
                            x.clone(),
                            Expr::Constant(2.0),
                        ),
                        Expr::Constant(1.0),
                    ),
                    n.clone(),
                )),
                "x".to_string(),
                Arc::new(n.clone()),
            ),
        )),
    )
}

/// Represents Laguerre's differential equation: `xy'' + (1-x)y' + ny = 0`.
///
/// This second-order linear ordinary differential equation arises in the quantum mechanical
/// treatment of the hydrogen atom. Its solutions are Laguerre polynomials.
///
/// # Arguments
/// * `y` - The unknown function `y(x)`.
/// * `x` - The independent variable `x`.
/// * `n` - The degree `n` of the Laguerre equation.
///
/// # Returns
/// An `Expr::Eq` representing Laguerre's differential equation.
#[must_use]

pub fn laguerre_differential_equation(
    y : &Expr,
    x : &Expr,
    n : &Expr,
) -> Expr {

    let y_prime = differentiate(y, "x");

    let y_double_prime = differentiate(&y_prime, "x");

    let term1 = Expr::new_mul(
        x.clone(),
        y_double_prime,
    );

    let term2 = Expr::new_mul(
        Expr::new_sub(
            Expr::Constant(1.0),
            x.clone(),
        ),
        y_prime,
    );

    let term3 = Expr::new_mul(n.clone(), y.clone());

    Expr::Eq(
        Arc::new(Expr::new_add(
            term1,
            Expr::new_add(term2, term3),
        )),
        Arc::new(Expr::Constant(0.0)),
    )
}

/// Represents Hermite's differential equation: `y'' - 2xy' + 2ny = 0`.
///
/// This second-order linear ordinary differential equation is important in quantum mechanics
/// (e.g., for the quantum harmonic oscillator) and probability theory. Its solutions are Hermite polynomials.
///
/// # Arguments
/// * `y` - The unknown function `y(x)`.
/// * `x` - The independent variable `x`.
/// * `n` - The degree `n` of the Hermite equation.
///
/// # Returns
/// An `Expr::Eq` representing Hermite's differential equation.
#[must_use]

pub fn hermite_differential_equation(
    y : &Expr,
    x : &Expr,
    n : &Expr,
) -> Expr {

    let y_prime = differentiate(y, "x");

    let y_double_prime = differentiate(&y_prime, "x");

    let term1 = y_double_prime;

    let term2 = Expr::new_mul(
        Expr::Constant(-2.0),
        Expr::new_mul(x.clone(), y_prime),
    );

    let term3 = Expr::new_mul(
        Expr::new_mul(
            Expr::Constant(2.0),
            n.clone(),
        ),
        y.clone(),
    );

    Expr::Eq(
        Arc::new(Expr::new_add(
            term1,
            Expr::new_add(term2, term3),
        )),
        Arc::new(Expr::Constant(0.0)),
    )
}

/// Represents Rodrigues' Formula for Hermite Polynomials: `H_n(x) = (-1)^n e^(x^2) d^n/dx^n [e^(-x^2)]`.
///
/// This formula provides a compact way to define Hermite polynomials and is useful
/// for deriving their properties.
///
/// # Arguments
/// * `n` - The degree `n` of the Hermite polynomial.
/// * `x` - The independent variable `x`.
///
/// # Returns
/// An `Expr::Eq` representing Rodrigues' formula.
#[must_use]

pub fn hermite_rodrigues_formula(
    n : &Expr,
    x : &Expr,
) -> Expr {

    Expr::Eq(
        Arc::new(hermite_h(
            n.clone(),
            x.clone(),
        )),
        Arc::new(Expr::new_mul(
            Expr::new_pow(
                Expr::Constant(-1.0),
                n.clone(),
            ),
            Expr::new_mul(
                Expr::new_exp(Expr::new_pow(
                    x.clone(),
                    Expr::Constant(2.0),
                )),
                Expr::DerivativeN(
                    Arc::new(Expr::new_exp(
                        Expr::new_neg(Expr::new_pow(
                            x.clone(),
                            Expr::Constant(2.0),
                        )),
                    )),
                    "x".to_string(),
                    Arc::new(n.clone()),
                ),
            ),
        )),
    )
}

/// Represents Chebyshev's differential equation: `(1-x²)y'' - xy' + n²y = 0`.
///
/// This differential equation has Chebyshev polynomials as its solutions.
///
/// # Arguments
/// * `y` - The unknown function.
/// * `x` - The independent variable.
/// * `n` - The degree.
///
/// # Returns
/// An `Expr::Eq` representing Chebyshev's differential equation.
#[must_use]

pub fn chebyshev_differential_equation(
    y : &Expr,
    x : &Expr,
    n : &Expr,
) -> Expr {

    let y_prime = differentiate(y, "x");

    let y_double_prime = differentiate(&y_prime, "x");

    let term1 = Expr::new_mul(
        Expr::new_sub(
            Expr::Constant(1.0),
            Expr::new_pow(
                x.clone(),
                Expr::Constant(2.0),
            ),
        ),
        y_double_prime,
    );

    let term2 = Expr::new_neg(Expr::new_mul(
        x.clone(),
        y_prime,
    ));

    let term3 = Expr::new_mul(
        Expr::new_pow(
            n.clone(),
            Expr::Constant(2.0),
        ),
        y.clone(),
    );

    Expr::Eq(
        Arc::new(Expr::new_add(
            term1,
            Expr::new_add(term2, term3),
        )),
        Arc::new(Expr::Constant(0.0)),
    )
}
