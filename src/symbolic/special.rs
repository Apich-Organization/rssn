//! # Numerical Special Functions
//!
//! This module provides numerical implementations of various special functions
//! commonly encountered in mathematics, physics, and engineering.
//!
//! ## Supported Functions
//!
//! ### Gamma and Beta Functions
//! - `gamma_numerical` - Gamma function Γ(x)
//! - `ln_gamma_numerical` - Natural log of gamma function ln(Γ(x))
//! - `digamma_numerical` - Digamma function ψ(x) = d/dx ln(Γ(x))
//! - `beta_numerical` - Beta function B(a, b)
//! - `ln_beta_numerical` - Natural log of beta function ln(B(a, b))
//! - `regularized_incomplete_beta` - Regularized incomplete beta Iₓ(a, b)
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
//! - `rising_factorial` - Rising factorial (Pochhammer symbol)
//! - `falling_factorial` - Falling factorial
//!
//! ### Bessel Functions
//! - `bessel_j0` - Bessel function of first kind J₀(x)
//! - `bessel_j1` - Bessel function of first kind J₁(x)
//! - `bessel_y0` - Bessel function of second kind Y₀(x)
//! - `bessel_y1` - Bessel function of second kind Y₁(x)
//! - `bessel_i0` - Modified Bessel function of first kind I₀(x)
//! - `bessel_i1` - Modified Bessel function of first kind I₁(x)
//! - `bessel_k0` - Modified Bessel function of second kind K₀(x)
//! - `bessel_k1` - Modified Bessel function of second kind K₁(x)
//!
//! ### Other Special Functions
//! - `sinc` - Normalized sinc function sin(πx)/(πx)
//! - `zeta` - Riemann zeta function ζ(s)
//!
//! ## Examples
//!
//! ### Gamma Function
//! ```
//! 
//! use rssn::symbolic::special::factorial;
//! use rssn::symbolic::special::gamma_numerical;
//!
//! // Γ(5) = 4! = 24
//! assert!((gamma_numerical(5.0) - 24.0).abs() < 1e-10);
//!
//! // factorial(4) = 24
//! assert_eq!(factorial(4), 24);
//! ```
//!
//! ### Error Functions
//! ```
//! 
//! use rssn::symbolic::special::erf_numerical;
//! use rssn::symbolic::special::erfc_numerical;
//!
//! // erf(0) = 0, erfc(0) = 1
//! assert!((erf_numerical(0.0)).abs() < 1e-10);
//!
//! assert!((erfc_numerical(0.0) - 1.0).abs() < 1e-10);
//!
//! // erf(x) + erfc(x) = 1
//! let x = 1.5;
//!
//! assert!(
//!     (erf_numerical(x) + erfc_numerical(x) - 1.0).abs()
//!         < 1e-10
//! );
//! ```
//!
//! ### Bessel Functions
//! ```
//! 
//! use rssn::symbolic::special::bessel_j0;
//! use rssn::symbolic::special::bessel_j1;
//!
//! // J₀(0) = 1, J₁(0) = 0
//! assert!((bessel_j0(0.0) - 1.0).abs() < 1e-6);
//!
//! assert!((bessel_j1(0.0)).abs() < 1e-6);
//! ```

use std::f64::consts::FRAC_2_PI;
use std::f64::consts::FRAC_PI_4;

use statrs::function::beta::beta;
use statrs::function::beta::ln_beta;
use statrs::function::erf::erf;
use statrs::function::erf::erfc;
use statrs::function::gamma::digamma;
use statrs::function::gamma::gamma;
use statrs::function::gamma::ln_gamma;

// ============================================================================
// Gamma and Related Functions
// ============================================================================

/// Computes the gamma function, `Γ(x)`.
///
/// The gamma function is an extension of the factorial function to real and complex numbers.
/// For positive integers `n`, `Γ(n) = (n-1)!`.
///
/// # Arguments
/// * `x` - The input value.
///
/// # Returns
/// The numerical value of `Γ(x)`.
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::special::gamma_numerical;
///
/// assert!((gamma_numerical(5.0) - 24.0).abs() < 1e-10);
///
/// assert!(
///     (gamma_numerical(0.5)
///         - std::f64::consts::PI.sqrt())
///     .abs()
///         < 1e-10
/// );
/// ```
#[must_use]

pub fn gamma_numerical(x: f64) -> f64 {

    gamma(x)
}

/// Computes the natural logarithm of the gamma function, `ln(Γ(x))`.
///
/// This function is often used to avoid overflow/underflow issues when `Γ(x)` itself
/// is very large or very small.
///
/// # Arguments
/// * `x` - The input value.
///
/// # Returns
/// The numerical value of `ln(Γ(x))`.
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::special::ln_gamma_numerical;
///
/// assert!(
///     (ln_gamma_numerical(5.0) - 24.0_f64.ln()).abs()
///         < 1e-10
/// );
/// ```
#[must_use]

pub fn ln_gamma_numerical(
    x: f64
) -> f64 {

    ln_gamma(x)
}

/// Computes the digamma function, `ψ(x) = d/dx ln(Γ(x))`.
///
/// The digamma function is the logarithmic derivative of the gamma function.
/// It appears in various formulas involving the gamma function.
///
/// # Arguments
/// * `x` - The input value.
///
/// # Returns
/// The numerical value of `ψ(x)`.
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::special::digamma_numerical;
///
/// // ψ(1) = -γ (Euler-Mascheroni constant)
/// let euler_mascheroni = 0.5772156649015329;
///
/// assert!(
///     (digamma_numerical(1.0) + euler_mascheroni).abs()
///         < 1e-10
/// );
/// ```
#[must_use]

pub fn digamma_numerical(
    x: f64
) -> f64 {

    digamma(x)
}

/// Computes the beta function, `B(a, b)`.
///
/// The beta function is closely related to the gamma function:
/// `B(a, b) = Γ(a)Γ(b) / Γ(a+b)`.
/// It appears in probability theory and mathematical physics.
///
/// # Arguments
/// * `a` - The first input value.
/// * `b` - The second input value.
///
/// # Returns
/// The numerical value of `B(a, b)`.
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::special::beta_numerical;
///
/// assert!((beta_numerical(1.0, 1.0) - 1.0).abs() < 1e-10);
///
/// assert!(
///     (beta_numerical(2.0, 2.0) - 1.0 / 6.0).abs()
///         < 1e-10
/// );
/// ```
#[must_use]

pub fn beta_numerical(
    a: f64,
    b: f64,
) -> f64 {

    beta(a, b)
}

/// Computes the natural logarithm of the beta function, `ln(B(a, b))`.
///
/// # Arguments
/// * `a` - The first input value.
/// * `b` - The second input value.
///
/// # Returns
/// The numerical value of `ln(B(a, b))`.
#[must_use]

pub fn ln_beta_numerical(
    a: f64,
    b: f64,
) -> f64 {

    ln_beta(a, b)
}

/// Computes the regularized incomplete beta function, `Iₓ(a, b)`.
///
/// The regularized incomplete beta function is defined as:
/// `Iₓ(a, b) = B(x; a, b) / B(a, b)`
/// where B(x; a, b) is the incomplete beta function.
///
/// This function is widely used in statistics, particularly for the
/// cumulative distribution function of the beta distribution and
/// related distributions.
///
/// # Arguments
/// * `a` - First shape parameter (> 0)
/// * `b` - Second shape parameter (> 0)
/// * `x` - Upper limit of integration (0 ≤ x ≤ 1)
///
/// # Returns
/// The value of the regularized incomplete beta function.
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::special::regularized_incomplete_beta;
///
/// // I₀(a, b) = 0 for any a, b > 0
/// assert!((regularized_incomplete_beta(2.0, 3.0, 0.0)).abs() < 1e-10);
///
/// // I₁(a, b) = 1 for any a, b > 0
/// assert!((regularized_incomplete_beta(2.0, 3.0, 1.0) - 1.0).abs() < 1e-10);
/// ```
#[must_use]

pub fn regularized_incomplete_beta(
    a: f64,
    b: f64,
    x: f64,
) -> f64 {

    statrs::function::beta::beta_reg(
        a, b, x,
    )
}

// ============================================================================
// Error Functions
// ============================================================================

/// Computes the error function, `erf(x)`.
///
/// The error function is a special function of sigmoid shape that arises in probability,
/// statistics, and partial differential equations.
///
/// `erf(x) = (2/√π) ∫₀ˣ e^(-t²) dt`
///
/// # Arguments
/// * `x` - The input value.
///
/// # Returns
/// The numerical value of `erf(x)`.
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::special::erf_numerical;
///
/// assert!((erf_numerical(0.0)).abs() < 1e-10);
///
/// assert!(
///     (erf_numerical(1.0) - 0.8427007929497148).abs()
///         < 1e-10
/// );
/// ```
#[must_use]

pub fn erf_numerical(x: f64) -> f64 {

    erf(x)
}

/// Computes the complementary error function, `erfc(x) = 1 - erf(x)`.
///
/// # Arguments
/// * `x` - The input value.
///
/// # Returns
/// The numerical value of `erfc(x)`.
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::special::erfc_numerical;
///
/// assert!((erfc_numerical(0.0) - 1.0).abs() < 1e-10);
/// ```
#[must_use]

pub fn erfc_numerical(x: f64) -> f64 {

    erfc(x)
}

/// Computes the inverse error function, `erf⁻¹(x)`.
///
/// Returns y such that erf(y) = x.
/// Uses Newton-Raphson iteration for computation.
///
/// # Arguments
/// * `x` - The input value (-1 < x < 1).
///
/// # Returns
/// The numerical value of `erf⁻¹(x)`.
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::special::erf_numerical;
/// use rssn::symbolic::special::inverse_erf;
///
/// let x = 0.5;
///
/// let y = inverse_erf(x);
///
/// assert!((erf_numerical(y) - x).abs() < 1e-6);
/// ```
#[must_use]

pub fn inverse_erf(x: f64) -> f64 {

    // Use rational approximation followed by Newton-Raphson refinement
    if x.abs() >= 1.0 {

        if (x - 1.0).abs()
            < f64::EPSILON
        {

            return f64::INFINITY;
        }

        if (x + 1.0).abs()
            < f64::EPSILON
        {

            return f64::NEG_INFINITY;
        }

        return f64::NAN;
    }

    if x == 0.0 {

        return 0.0;
    }

    // Initial approximation using Winitzki's formula
    let a = 0.147;

    let ln_term = (x
        .mul_add(-x, 1.0)
        .ln())
    .abs();

    let term1 = 2.0
        / (std::f64::consts::PI * a)
        + ln_term / 2.0;

    let sign = if x < 0.0 {

        -1.0
    } else {

        1.0
    };

    let approx = sign
        * (term1.sqrt()
            - (2.0
                / (std::f64::consts::PI
                    * a)
                + ln_term / 2.0)
                .sqrt())
        .sqrt();

    // Refine with Newton-Raphson: y_{n+1} = y_n - (erf(y_n) - x) / erf'(y_n)
    // erf'(x) = 2/sqrt(pi) * exp(-x^2)
    let two_over_sqrt_pi = 2.0
        / std::f64::consts::PI.sqrt();

    let mut y = approx;

    for _ in 0 .. 5 {

        let erf_y = erf(y);

        let derivative =
            two_over_sqrt_pi
                * (-y * y).exp();

        if derivative.abs() < 1e-15 {

            break;
        }

        y -= (erf_y - x) / derivative;
    }

    y
}

/// Computes the inverse complementary error function, `erfc⁻¹(x)`.
///
/// Returns y such that erfc(y) = x.
///
/// # Arguments
/// * `x` - The input value (0 < x < 2).
///
/// # Returns
/// The numerical value of `erfc⁻¹(x)`.
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::special::erfc_numerical;
/// use rssn::symbolic::special::inverse_erfc;
///
/// let x = 0.5;
///
/// let y = inverse_erfc(x);
///
/// assert!((erfc_numerical(y) - x).abs() < 1e-6);
/// ```
#[must_use]

pub fn inverse_erfc(x: f64) -> f64 {

    // erfc^{-1}(x) = erf^{-1}(1 - x)
    inverse_erf(1.0 - x)
}

// ============================================================================
// Combinatorial Functions
// ============================================================================

/// Computes the factorial, `n!`.
///
/// For non-negative integers: n! = n × (n-1) × ... × 2 × 1
/// By convention, 0! = 1.
///
/// # Arguments
/// * `n` - A non-negative integer.
///
/// # Returns
/// The factorial of n as u64. May overflow for n > 20.
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::special::factorial;
///
/// assert_eq!(factorial(0), 1);
///
/// assert_eq!(factorial(5), 120);
///
/// assert_eq!(
///     factorial(10),
///     3628800
/// );
/// ```
#[must_use]

pub fn factorial(n: u64) -> u64 {

    (1 ..= n).product()
}

/// Computes the double factorial, `n!!`.
///
/// For odd n: n!! = n × (n-2) × ... × 3 × 1
/// For even n: n!! = n × (n-2) × ... × 4 × 2
/// By convention, 0!! = 1 and (-1)!! = 1.
///
/// # Arguments
/// * `n` - A non-negative integer.
///
/// # Returns
/// The double factorial of n.
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::special::double_factorial;
///
/// assert_eq!(
///     double_factorial(0),
///     1
/// );
///
/// assert_eq!(
///     double_factorial(5),
///     15
/// ); // 5 * 3 * 1
/// assert_eq!(
///     double_factorial(6),
///     48
/// ); // 6 * 4 * 2
/// ```
#[must_use]

pub const fn double_factorial(
    n: u64
) -> u64 {

    if n <= 1 {

        return 1;
    }

    let mut result = 1u64;

    let mut k = n;

    while k > 1 {

        result *= k;

        k -= 2;
    }

    result
}

/// Computes the binomial coefficient, C(n, k) = n! / (k! × (n-k)!).
///
/// Also known as "n choose k", this counts the number of ways to choose
/// k items from n items without regard to order.
///
/// # Arguments
/// * `n` - Total number of items.
/// * `k` - Number of items to choose.
///
/// # Returns
/// The binomial coefficient as u64. Returns 0 if k > n.
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::special::binomial;
///
/// assert_eq!(binomial(5, 0), 1);
///
/// assert_eq!(binomial(5, 2), 10);
///
/// assert_eq!(binomial(5, 5), 1);
///
/// assert_eq!(binomial(10, 3), 120);
/// ```
#[must_use]

pub fn binomial(
    n: u64,
    k: u64,
) -> u64 {

    if k > n {

        return 0;
    }

    if k == 0 || k == n {

        return 1;
    }

    // Use the smaller of k and n-k for efficiency
    let k = if k > n - k {

        n - k
    } else {

        k
    };

    let mut result = 1u64;

    for i in 0 .. k {

        result =
            result * (n - i) / (i + 1);
    }

    result
}

/// Computes the rising factorial (Pochhammer symbol), `(x)ₙ`.
///
/// The rising factorial is defined as:
/// `(x)ₙ = x × (x+1) × (x+2) × ... × (x+n-1)`
///
/// # Arguments
/// * `x` - Base value.
/// * `n` - Number of terms.
///
/// # Returns
/// The rising factorial as f64.
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::special::rising_factorial;
///
/// assert!(
///     (rising_factorial(3.0, 4) - 360.0).abs() < 1e-10
/// ); // 3 * 4 * 5 * 6
/// assert!(
///     (rising_factorial(1.0, 5) - 120.0).abs() < 1e-10
/// ); // 1 * 2 * 3 * 4 * 5 = 5!
/// ```
#[must_use]

pub fn rising_factorial(
    x: f64,
    n: u32,
) -> f64 {

    if n == 0 {

        return 1.0;
    }

    (0 .. n).fold(1.0, |acc, i| {

        acc * (x + f64::from(i))
    })
}

/// Computes the falling factorial, `(x)₍ₙ₎`.
///
/// The falling factorial is defined as:
/// `(x)₍ₙ₎ = x × (x-1) × (x-2) × ... × (x-n+1)`
///
/// # Arguments
/// * `x` - Base value.
/// * `n` - Number of terms.
///
/// # Returns
/// The falling factorial as f64.
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::special::falling_factorial;
///
/// assert!(
///     (falling_factorial(5.0, 3) - 60.0).abs() < 1e-10
/// ); // 5 * 4 * 3
/// assert!(
///     (falling_factorial(10.0, 4) - 5040.0).abs() < 1e-10
/// ); // 10 * 9 * 8 * 7
/// ```
#[must_use]

pub fn falling_factorial(
    x: f64,
    n: u32,
) -> f64 {

    if n == 0 {

        return 1.0;
    }

    (0 .. n).fold(1.0, |acc, i| {

        acc * (x - f64::from(i))
    })
}

// ============================================================================
// Bessel Functions
// ============================================================================

/// Computes the Bessel function of the first kind, J₀(x).
///
/// Uses polynomial approximation for small x and asymptotic expansion for large x.
///
/// # Arguments
/// * `x` - The input value.
///
/// # Returns
/// The numerical value of J₀(x).
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::special::bessel_j0;
///
/// assert!((bessel_j0(0.0) - 1.0).abs() < 1e-6);
/// ```
#[must_use]

pub fn bessel_j0(x: f64) -> f64 {

    let ax = x.abs();

    if ax < 8.0 {
        let y = x * x;
        
        // Numerator coefficients
        let ans1 = y.mul_add(
            y.mul_add(
                y.mul_add(
                    y.mul_add(y.mul_add(-184.9052456, 77392.33017), -11214424.18),
                    651619640.7,
                ),
                -13362590354.0,
            ),
            57568490574.0,
        );

        // Denominator coefficients
        let ans2 = y.mul_add(
            y.mul_add(
                y.mul_add(
                    y.mul_add(y.mul_add(1.0, 267.8532712), 59272.64853),
                    9494680.718,
                ),
                1029532985.0,
            ),
            57568490411.0,
        );

        return ans1 / ans2;
    } else {
        // Asymptotic form
        let z = 8.0 / ax;

        let y = z * z;

        let xx = ax - FRAC_PI_4; // ax - pi/4
        let ans1 = 1.0
            + y * (-0.109_862_862_7e-2
                + y * (0.273_451_040_7e-4 + y * (-0.207_337_063_9e-5 + y * 0.209_388_721_1e-6)));

        let ans2 = -0.015_624_999_95
            + y * (0.143_048_876_5e-3
                + y * (-0.691_114_765_1e-5 + y * (0.762_109_516_1e-6 - y * 0.934_945_152e-7)));

        (FRAC_2_PI / ax).sqrt()
            * xx.cos().mul_add(
                ans1,
                -(z * xx.sin() * ans2),
            )
    }
}

/// Computes the Bessel function of the first kind, J₁(x).
///
/// Uses polynomial approximation for small x and asymptotic expansion for large x.
///
/// # Arguments
/// * `x` - The input value.
///
/// # Returns
/// The numerical value of J₁(x).
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::special::bessel_j1;
///
/// assert!((bessel_j1(0.0)).abs() < 1e-10);
/// ```
#[must_use]

pub fn bessel_j1(x: f64) -> f64 {

    let ax = x.abs();

    if ax < 8.0 {

        let y = x * x;

        let ans1 = x
            * y.mul_add(
                y.mul_add(
                    y.mul_add(
                        y.mul_add(
                            y.mul_add(-30.160_366_06, 15_704.482_60),
                            -2_972_611.439,
                        ),
                        242_396_853.1,
                    ),
                    -7_895_059_235.0,
                ),
                72_362_614_232.0,
            );

        let ans2 = y.mul_add(
            y.mul_add(
                y.mul_add(
                    y.mul_add(
                        y.mul_add(1.0, 376.999_139_7),
                        99_447.433_94,
                    ),
                    18_583_304.74,
                ),
                2_300_535_178.0,
            ),
            144_725_228_442.0,
        );

        ans1 / ans2
    } else {

        let z = 8.0 / ax;

        let y = z * z;

        let xx = ax - 2.356_194_491; // ax - 3*pi/4
        let ans1 = 1.0
            + y * (0.183_105e-2
                + y * (-0.351_639_649_6e-4 + y * (0.245_752_017_4e-5 - y * 0.240_337_019e-6)));

        let ans2 = 0.046_874_999_95
            + y * (-0.200_269_087_3e-3
                + y * (0.844_919_909_6e-5 + y * (-0.882_289_87e-6 + y * 0.105_787_412e-6)));

        let result = (FRAC_2_PI / ax)
            .sqrt()
            * xx.cos().mul_add(
                ans1,
                -(z * xx.sin() * ans2),
            );

        if x < 0.0 {

            -result
        } else {

            result
        }
    }
}

/// Computes the Bessel function of the second kind, Y₀(x).
///
/// Also known as the Neumann function N₀(x).
///
/// # Arguments
/// * `x` - The input value (x > 0).
///
/// # Returns
/// The numerical value of Y₀(x).
#[must_use]

pub fn bessel_y0(x: f64) -> f64 {

    if x < 0.0 {

        return f64::NAN;
    }

    if x < 8.0 {

        let y = x * x;

        let ans1 = y.mul_add(
            y.mul_add(
                y.mul_add(
                    y.mul_add(
                        y.mul_add(
                            228.462_273_3,
                            -86_327.927_57,
                        ),
                        10_879_881.29,
                    ),
                    -512_359_803.6,
                ),
                7_062_834_065.0,
            ),
            -2_957_821_389.0,
        );

        let ans2 = y.mul_add(
            y.mul_add(
                y.mul_add(
                    y.mul_add(
                        y.mul_add(1.0, 226.103_024_4),
                        47_447.264_70,
                    ),
                    7_189_466.438,
                ),
                745_249_964.8,
            ),
            40_076_544_269.0,
        );

        (FRAC_2_PI * bessel_j0(x))
            .mul_add(
                x.ln(),
                ans1 / ans2,
            )
    } else {

        let z = 8.0 / x;

        let y = z * z;

        let xx = x - FRAC_PI_4;

        let ans1 = 1.0
            + y * (-0.001_098_628_627e-2
                + y * (0.273_451_040_7e-4 + y * (-0.207_337_063_9e-5 + y * 0.209_388_721_1e-6)));

        let ans2 = -0.015_624_999_95
            + y * (0.143_048_876_5e-3
                + y * (-0.691_114_765_1e-5 + y * (0.762_109_516_1e-6 - y * 0.934_945_152e-7)));

        (FRAC_2_PI / x).sqrt()
            * xx.sin().mul_add(
                ans1,
                z * xx.cos() * ans2,
            )
    }
}

/// Computes the Bessel function of the second kind, Y₁(x).
///
/// Also known as the Neumann function N₁(x).
///
/// # Arguments
/// * `x` - The input value (x > 0).
///
/// # Returns
/// The numerical value of Y₁(x).
#[must_use]

pub fn bessel_y1(x: f64) -> f64 {

    if x < 0.0 {

        return f64::NAN;
    }

    if x < 8.0 {

        let y = x * x;

        let ans1 = x
            * y.mul_add(
                y.mul_add(
                    y.mul_add(
                        y.mul_add(
                            y.mul_add(8_511.937_935, -4_237_922.726),
                            734_926_455.1,
                        ),
                        -51_534_381_390.0,
                    ),
                    1_275_274_390_000.0,
                ),
                -4_900_604_943_000.0,
            );

        let ans2 = y.mul_add(
            y.mul_add(
                y.mul_add(
                    y.mul_add(
                        y.mul_add(
                            y.mul_add(1.0, 354.963_288_5),
                            102_042.605,
                        ),
                        22_459_040.02,
                    ),
                    3_733_650_367.0,
                ),
                4_244_419_664_000.0,
            ),
            24_909_857_500_000.0,
        );

        FRAC_2_PI.mul_add(
            bessel_j1(x).mul_add(
                x.ln(),
                -(1.0 / x),
            ),
            ans1 / ans2,
        )
    } else {

        let z = 8.0 / x;

        let y = z * z;

        let xx = x - 2.356_194_491;

        let ans1 = 1.0
            + y * (0.183_105e-2
                + y * (-0.351_639_649_6e-4 + y * (0.245_752_017_4e-5 - y * 0.240_337_019e-6)));

        let ans2 = 0.046_874_999_95
            + y * (-0.200_269_087_3e-3
                + y * (0.844_919_909_6e-5 + y * (-0.882_289_87e-6 + y * 0.105_787_412e-6)));

        (FRAC_2_PI / x).sqrt()
            * xx.sin().mul_add(
                ans1,
                z * xx.cos() * ans2,
            )
    }
}

/// Computes the modified Bessel function of the first kind, I₀(x).
///
/// # Arguments
/// * `x` - The input value.
///
/// # Returns
/// The numerical value of I₀(x).
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::special::bessel_i0;
///
/// assert!((bessel_i0(0.0) - 1.0).abs() < 1e-10);
/// ```
#[must_use]

pub fn bessel_i0(x: f64) -> f64 {

    let ax = x.abs();

    if ax < 3.75 {

        let y = (x / 3.75).powi(2);

        1.0 + y
            * (3.515_622_9
                + y * (3.089_942_4
                    + y * (1.206_749_2 + y * (0.265_973_2 + y * (0.036_076_8 + y * 0.004_581_3)))))
    } else {

        let y = 3.75 / ax;

        let ans = 0.398_942_28
            + y * (0.013_285_92
                + y * (0.002_253_19
                    + y * (-0.001_575_65
                        + y * (0.009_162_81
                            + y * (-0.020_577_06
                                + y * (0.026_355_37 + y * (-0.016_476_33 + y * 0.003_923_77)))))));

        ans * ax.exp() / ax.sqrt()
    }
}

/// Computes the modified Bessel function of the first kind, I₁(x).
///
/// # Arguments
/// * `x` - The input value.
///
/// # Returns
/// The numerical value of I₁(x).
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::special::bessel_i1;
///
/// assert!((bessel_i1(0.0)).abs() < 1e-10);
/// ```
#[must_use]

pub fn bessel_i1(x: f64) -> f64 {

    let ax = x.abs();

    if ax < 3.75 {

        let y = (x / 3.75).powi(2);

        let ans = ax
            * (0.5
                + y * (0.878_905_94
                    + y * (0.514_988_69
                        + y * (0.150_849_34
                            + y * (0.026_587_33 + y * (0.003_015_32 + y * 0.000_324_11))))));

        if x < 0.0 {

            -ans
        } else {

            ans
        }
    } else {

        let y = 3.75 / ax;

        let ans = 0.398_942_28
            + y * (-0.039_880_24
                + y * (-0.003_620_18
                    + y * (0.001_638_01
                        + y * (-0.010_315_55
                            + y * (0.022_829_67
                                + y * (-0.028_953_12 + y * (0.017_876_54 - y * 0.004_200_59)))))));

        let result =
            ans * ax.exp() / ax.sqrt();

        if x < 0.0 {

            -result
        } else {

            result
        }
    }
}

/// Computes the modified Bessel function of the second kind, K₀(x).
///
/// # Arguments
/// * `x` - The input value (x > 0).
///
/// # Returns
/// The numerical value of K₀(x).
#[must_use]

pub fn bessel_k0(x: f64) -> f64 {

    if x <= 0.0 {

        return f64::NAN;
    }

    if x <= 2.0 {

        let y = x * x / 4.0;

        (-x.ln() + 0.577_215_664_901_532_9).mul_add(bessel_i0(x), -0.577_215_664_901_532_9 + y * (0.422_784_2
                    + y * (0.230_697_56
                        + y * (0.034_885_9 + y * (0.002_626_98 + y * (0.000_107_5 + y * 0.000_007_4))))))
    } else {

        let y = 2.0 / x;

        (std::f64::consts::FRAC_PI_2 / x).sqrt()
            * (-x).exp()
            * (1.253_314_14
                + y * (-0.078_323_58
                    + y * (0.021_895_68
                        + y * (-0.010_624_46
                            + y * (0.005_878_72 + y * (-0.002_515_40 + y * 0.000_532_08))))))
    }
}

/// Computes the modified Bessel function of the second kind, K₁(x).
///
/// # Arguments
/// * `x` - The input value (x > 0).
///
/// # Returns
/// The numerical value of K₁(x).
#[must_use]

pub fn bessel_k1(x: f64) -> f64 {

    if x <= 0.0 {

        return f64::NAN;
    }

    if x <= 2.0 {

        let y = x * x / 4.0;

        (x.ln() - 0.577_215_664_901_532_9).mul_add(bessel_i1(x), (1.0 / x) * (1.0
                    + y * (0.154_431_44
                        + y * (-0.672_785_79
                            + y * (-0.181_568_97
                                + y * (-0.019_194_02 + y * (-0.001_104_04 - y * 0.000_046_86)))))))
    } else {

        let y = 2.0 / x;

        (std::f64::consts::FRAC_PI_2 / x).sqrt()
            * (-x).exp()
            * (1.253_314_14
                + y * (0.234_986_19
                    + y * (-0.036_556_20
                        + y * (0.015_042_68
                            + y * (-0.007_803_53 + y * (0.003_256_14 - y * 0.000_682_45))))))
    }
}

// ============================================================================
// Other Special Functions
// ============================================================================

/// Computes the normalized sinc function, `sinc(x) = sin(πx) / (πx)`.
///
/// By convention, sinc(0) = 1.
///
/// # Arguments
/// * `x` - The input value.
///
/// # Returns
/// The numerical value of sinc(x).
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::special::sinc;
///
/// assert!((sinc(0.0) - 1.0).abs() < 1e-10);
///
/// assert!((sinc(1.0)).abs() < 1e-10); // sin(π)/π = 0
/// ```
#[must_use]

pub fn sinc(x: f64) -> f64 {

    if x.abs() < 1e-15 {

        1.0
    } else {

        let pi_x =
            std::f64::consts::PI * x;

        pi_x.sin() / pi_x
    }
}

/// Computes the Riemann zeta function, `ζ(s)`.
///
/// The Riemann zeta function is defined as:
/// `ζ(s) = Σₙ₌₁^∞ 1/nˢ` for Re(s) > 1
///
/// This implementation uses numerical approximation.
///
/// # Arguments
/// * `s` - The input value (s > 1 for convergence).
///
/// # Returns
/// The numerical value of ζ(s).
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::special::zeta;
///
/// // ζ(2) = π²/6 (approximate numerical computation)
/// let expected = std::f64::consts::PI.powi(2) / 6.0;
///
/// assert!((zeta(2.0) - expected).abs() < 0.001); // ~0.1% relative error
/// ```
#[must_use]

pub fn zeta(s: f64) -> f64 {

    if s <= 1.0 {

        return f64::NAN; // Not defined for s <= 1 in this simple implementation
    }

    // Use the Euler-Maclaurin formula for better convergence
    let n_terms = 100;

    let mut sum = 0.0;

    // Direct summation for first n_terms
    for n in 1 ..= n_terms {

        sum +=
            1.0 / f64::from(n).powf(s);
    }

    // Euler-Maclaurin correction
    let n = f64::from(n_terms);

    sum += n.powf(1.0 - s) / (s - 1.0);

    sum += 0.5 * n.powf(-s);

    sum += s / 12.0 * n.powf(-s - 1.0);

    sum
}

/// Computes the logarithm of the factorial, `ln(n!)`.
///
/// Uses the log-gamma function for computation, avoiding overflow for large n.
///
/// # Arguments
/// * `n` - A non-negative integer.
///
/// # Returns
/// The natural logarithm of n!.
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::special::ln_factorial;
///
/// assert!(
///     (ln_factorial(5) - 120.0_f64.ln()).abs() < 1e-10
/// );
/// ```
#[must_use]

pub fn ln_factorial(n: u64) -> f64 {

    if n <= 1 {

        return 0.0;
    }

    ln_gamma(n as f64 + 1.0)
}

/// Computes the regularized lower incomplete gamma function, P(a, x).
///
/// P(a, x) = γ(a, x) / Γ(a) where γ(a, x) is the lower incomplete gamma.
///
/// # Arguments
/// * `a` - Shape parameter (a > 0).
/// * `x` - Upper limit of integration (x >= 0).
///
/// # Returns
/// The value of the regularized lower incomplete gamma.
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::special::regularized_gamma_p;
///
/// // P(1, x) = 1 - e^(-x) for the exponential distribution
/// let x = 2.0;
///
/// assert!(
///     (regularized_gamma_p(1.0, x) - (1.0 - (-x).exp()))
///         .abs()
///         < 1e-10
/// );
/// ```
#[must_use]

pub fn regularized_gamma_p(
    a: f64,
    x: f64,
) -> f64 {

    statrs::function::gamma::gamma_lr(
        a, x,
    )
}

/// Computes the regularized upper incomplete gamma function, Q(a, x).
///
/// Q(a, x) = Γ(a, x) / Γ(a) = 1 - P(a, x)
///
/// # Arguments
/// * `a` - Shape parameter (a > 0).
/// * `x` - Lower limit of integration (x >= 0).
///
/// # Returns
/// The value of the regularized upper incomplete gamma.
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::special::regularized_gamma_q;
///
/// // Q(a, x) + P(a, x) = 1
/// let a = 2.0;
///
/// let x = 1.5;
///
/// use rssn::symbolic::special::regularized_gamma_p;
///
/// assert!(
///     (regularized_gamma_p(a, x)
///         + regularized_gamma_q(a, x)
///         - 1.0)
///         .abs()
///         < 1e-10
/// );
/// ```
#[must_use]

pub fn regularized_gamma_q(
    a: f64,
    x: f64,
) -> f64 {

    statrs::function::gamma::gamma_ur(
        a, x,
    )
}
