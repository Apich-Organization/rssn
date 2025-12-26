//! # Numerical Integration (Quadrature)
//!
//! This module provides various numerical integration (quadrature) methods for approximating
//! definite integrals of functions. It includes implementations of the Trapezoidal rule,
//! Simpson's rule, Adaptive Quadrature (Adaptive Simpson's), Romberg Integration, and
//! Gauss-Legendre Quadrature.
//!
//! # Methods
//!
//! * **Trapezoidal Rule**: Simple and robust, linear approximation.
//! * **Simpson's Rule**: Quadratic approximation, generally more accurate than Trapezoidal for smooth functions.
//! * **Adaptive Quadrature**: Recursively subdivides intervals to achieve a specified error tolerance.
//! * **Romberg Integration**: Uses Richardson extrapolation on Trapezoidal approximation to achieve high order accuracy.
//! * **Gauss-Legendre Quadrature**: Uses optimal sample points (roots of Legendre polynomials) for high accuracy with fewer function evaluations, ideal for smooth functions.

use std::collections::HashMap;

use serde::Deserialize;
use serde::Serialize;

use crate::numerical::elementary::eval_expr;
use crate::symbolic::core::Expr;

/// Enum to select the numerical integration method.
#[derive(
    Debug,
    Clone,
    Copy,
    PartialEq,
    Eq,
    Serialize,
    Deserialize,
)]

pub enum QuadratureMethod {
    /// Trapezoidal rule.
    Trapezoidal,
    /// Simpson's 1/3 rule.
    Simpson,
    /// Adaptive Simpson's quadrature.
    Adaptive,
    /// Romberg integration.
    Romberg,
    /// Gauss-Legendre quadrature (n=5).
    GaussLegendre,
}

/// # Pure Numerical Trapezoidal Rule
///
/// Performs numerical integration using the trapezoidal rule for a given function.
///
/// ## Arguments
/// * `f` - The function to integrate, as a closure.
/// * `range` - A tuple `(a, b)` representing the integration interval.
/// * `n_steps` - The number of steps to use.
///
/// ## Returns
/// The approximate value of the definite integral.
///
/// ## Example
/// ```
/// 
/// use rssn::numerical::integrate::trapezoidal_rule;
///
/// let f = |x: f64| x * x;
///
/// let res = trapezoidal_rule(f, (0.0, 1.0), 1000);
///
/// assert!((res - 1.0 / 3.0).abs() < 1e-4);
/// ```

pub fn trapezoidal_rule<F>(
    f : F,
    range : (f64, f64),
    n_steps : usize,
) -> f64
where
    F : Fn(f64) -> f64,
{

    let (a, b) = range;

    if n_steps == 0 {

        return 0.0;
    }

    // Correctly handle a == b (integral is 0)
    // and a > b (integral is negative of b to a)
    // The logic below works for a > b as h will be negative.
    if (a - b).abs() < f64::EPSILON {

        return 0.0;
    }

    let h = (b - a) / (n_steps as f64);

    let mut sum = 0.5 * (f(a) + f(b));

    for i in 1 .. n_steps {

        let x = a + (i as f64) * h;

        sum += f(x);
    }

    h * sum
}

/// # Pure Numerical Simpson's Rule
///
/// Performs numerical integration using Simpson's 1/3 rule.
///
/// ## Arguments
/// * `f` - The function to integrate.
/// * `range` - Interval `(a, b)`.
/// * `n_steps` - Number of steps (must be usually even, but we handle it).
///
/// ## Returns
/// `Result<f64, String>`
///
/// ## Example
/// ```
/// 
/// use rssn::numerical::integrate::simpson_rule;
///
/// let f = |x: f64| x * x;
///
/// let res = simpson_rule(f, (0.0, 1.0), 10).unwrap();
///
/// assert!((res - 1.0 / 3.0).abs() < 1e-10);
/// ```

pub fn simpson_rule<F>(
    f : F,
    range : (f64, f64),
    n_steps : usize,
) -> Result<f64, String>
where
    F : Fn(f64) -> f64,
{

    let (a, b) = range;

    if n_steps == 0 {

        return Ok(0.0);
    }

    if (a - b).abs() < f64::EPSILON {

        return Ok(0.0);
    }

    // Simpson's rule requires even number of intervals for the strict global formula.
    // If odd, we can warn or adjust. For now, enforce even.
    let steps = if n_steps % 2 != 0 {

        n_steps + 1
    } else {

        n_steps
    };

    let h = (b - a) / (steps as f64);

    let mut sum = f(a) + f(b);

    for i in 1 .. steps {

        let x = a + (i as f64) * h;

        let weight = if i % 2 == 0 {

            2.0
        } else {

            4.0
        };

        sum += weight * f(x);
    }

    Ok((h / 3.0) * sum)
}

/// # Adaptive Quadrature (Adaptive Simpson's Method)
///
/// Recursively refines the interval until the estimated error satisfies the tolerance.
///
/// ## Arguments
/// * `f` - Function to integrate.
/// * `range` - Interval `(a, b)`.
/// * `tolerance` - Error tolerance (e.g., 1e-6).
///
/// ## Example
/// ```
/// 
/// use rssn::numerical::integrate::adaptive_quadrature;
///
/// let f = |x: f64| x.sin();
///
/// let res = adaptive_quadrature(
///     f,
///     (
///         0.0,
///         std::f64::consts::PI,
///     ),
///     1e-6,
/// );
///
/// assert!((res - 2.0).abs() < 1e-6);
/// ```

pub fn adaptive_quadrature<F>(
    f : F,
    range : (f64, f64),
    tolerance : f64,
) -> f64
where
    F : Fn(f64) -> f64,
{

    // Inner recursive function
    fn adaptive_recursive<F>(
        f : &F,
        a : f64,
        b : f64,
        eps : f64,
        whole_simpson : f64,
        limit : usize,
    ) -> f64
    where
        F : Fn(f64) -> f64,
    {

        if limit == 0 {

            // Recursion limit reached, return current best guess
            return whole_simpson;
        }

        let mid = (a + b) / 2.0;

        let sub_mid_left =
            (a + mid) / 2.0;

        let sub_mid_right =
            (mid + b) / 2.0;

        let fa = f(a);

        let fb = f(b);

        let fm = f(mid);

        let fml = f(sub_mid_left);

        let fmr = f(sub_mid_right);

        // Simp(a, b) = (b-a)/6 * (f(a) + 4f(m) + f(b))
        let left_simpson = (mid - a)
            / 6.0
            * (fa + 4.0 * fml + fm);

        let right_simpson = (b - mid)
            / 6.0
            * (fm + 4.0 * fmr + fb);

        let sum_halves = left_simpson
            + right_simpson;

        // Error estimate (1/15 rule)
        let error = (sum_halves
            - whole_simpson)
            .abs()
            / 15.0;

        if error <= eps {

            // Richardson extrapolation: S + (S - S_whole)/15
            sum_halves
                + (sum_halves
                    - whole_simpson)
                    / 15.0
        } else {

            adaptive_recursive(
                f,
                a,
                mid,
                eps / 2.0,
                left_simpson,
                limit - 1,
            ) + adaptive_recursive(
                f,
                mid,
                b,
                eps / 2.0,
                right_simpson,
                limit - 1,
            )
        }
    }

    let (a, b) = range;

    if (a - b).abs() < f64::EPSILON {

        return 0.0;
    }

    // Initial Simpson estimate
    let mid = (a + b) / 2.0;

    let fm = f(mid);

    let initial_simpson = (b - a) / 6.0
        * (f(a) + 4.0 * fm + f(b));

    adaptive_recursive(
        &f,
        a,
        b,
        tolerance,
        initial_simpson,
        100,
    ) // limit depth to avoid stack overflow
}

/// # Romberg Integration
///
/// Uses Richardson extrapolation on the Trapezoidal rule to improve accuracy.
///
/// ## Arguments
/// * `f` - Function to integrate.
/// * `range` - Interval.
/// * `max_steps` - Order of extrapolation (e.g., 5-10).
///
/// ## Example
/// ```
/// 
/// use rssn::numerical::integrate::romberg_integration;
///
/// let f = |x: f64| x.exp();
///
/// let res = romberg_integration(f, (0.0, 1.0), 6);
///
/// assert!(
///     (res - (std::f64::consts::E - 1.0)).abs() < 1e-10
/// );
/// ```

pub fn romberg_integration<F>(
    f : F,
    range : (f64, f64),
    max_steps : usize,
) -> f64
where
    F : Fn(f64) -> f64,
{

    let (a, b) = range;

    if max_steps == 0 {

        return 0.0;
    }

    if (a - b).abs() < f64::EPSILON {

        return 0.0;
    }

    let mut r =
        vec![
            vec![0.0; max_steps];
            max_steps
        ];

    // R[0][0]
    let h = b - a;

    r[0][0] = 0.5 * h * (f(a) + f(b));

    for i in 1 .. max_steps {

        // Calculate R[i][0] using Trapezoidal rule with 2^i segments
        // But we can update from R[i-1][0] efficiently
        let steps_prev = 1 << (i - 1);

        let h_i = h / (1 << i) as f64;

        let mut sum = 0.0;

        for k in 1 ..= steps_prev {

            let x = a
                + (2 * k - 1) as f64
                    * h_i;

            sum += f(x);
        }

        r[i][0] = 0.5 * r[i - 1][0]
            + h_i * sum;

        // Richardson extrapolation
        for j in 1 ..= i {

            let k =
                4.0_f64.powi(j as i32);

            r[i][j] = (k * r[i][j - 1]
                - r[i - 1][j - 1])
                / (k - 1.0);
        }
    }

    r[max_steps - 1][max_steps - 1]
}

/// # Gauss-Legendre Quadrature
///
/// Uses standard weights and nodes for n = 5 (hardcoded for now as it's efficient for general purpose).
/// Can be extended to arbitrary n in the future.
///
/// ## Example
/// ```
/// 
/// use rssn::numerical::integrate::gauss_legendre_quadrature;
///
/// let f = |x: f64| x.powi(3);
///
/// let res = gauss_legendre_quadrature(f, (0.0, 1.0));
///
/// assert!((res - 0.25).abs() < 1e-10);
/// ```

pub fn gauss_legendre_quadrature<F>(
    f : F,
    range : (f64, f64),
) -> f64
where
    F : Fn(f64) -> f64,
{

    let (a, b) = range;

    if (a - b).abs() < f64::EPSILON {

        return 0.0;
    }

    let mid = (a + b) / 2.0;

    let half_len = (b - a) / 2.0;

    // Nodes and weights for n=5 (from standard tables)
    // x_i are for interval [-1, 1]
    let nodes = [
        0.0,
        0.5384693101056831,
        -0.5384693101056831,
        0.906179845938664,
        -0.906179845938664,
    ];

    let weights = [
        0.5688888888888889,
        0.4786286704993665,
        0.4786286704993665,
        0.2369268850561891,
        0.2369268850561891,
    ];

    let mut sum = 0.0;

    for i in 0 .. 5 {

        // Transform x from [-1, 1] to [a, b]
        let x =
            mid + half_len * nodes[i];

        sum += weights[i] * f(x);
    }

    half_len * sum
}

/// Performs numerical integration (quadrature) of a function `f(x)` over an interval `[a, b]`.
///
/// # Arguments
/// * `f` - The expression to integrate.
/// * `var` - The variable of integration.
/// * `range` - A tuple `(a, b)` representing the integration interval.
/// * `n_steps` - The number of steps to use for the integration (for non-adaptive methods).
/// * `method` - The quadrature method to use.
///
/// # Returns
/// A `Result` containing the numerical value of the integral, or an error string.

pub fn quadrature(
    f : &Expr,
    var : &str,
    range : (f64, f64),
    n_steps : usize,
    method : &QuadratureMethod,
) -> Result<f64, String> {

    let func = |x : f64| -> f64 {

        let mut vars = HashMap::new();

        vars.insert(var.to_string(), x);

        eval_expr(f, &vars)
            .unwrap_or(f64::NAN)
    };

    let result = match method {
        | QuadratureMethod::Trapezoidal => {
            trapezoidal_rule(func, range, n_steps)
        },
        | QuadratureMethod::Simpson => {
            simpson_rule(func, range, n_steps)?
        },
        | QuadratureMethod::Adaptive => {
            adaptive_quadrature(func, range, 1e-6)
        },
        | QuadratureMethod::Romberg => {
            romberg_integration(func, range, 6)
        },
        | QuadratureMethod::GaussLegendre => {
            gauss_legendre_quadrature(func, range)
        },
    };

    if result.is_nan() {

        return Err("Integration \
                    resulted in \
                    NaN, likely due \
                    to evaluation \
                    error."
            .to_string());
    }

    Ok(result)
}
