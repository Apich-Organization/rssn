//! # Numerical Integration (Quadrature)
//!
//! This module provides numerical integration (quadrature) methods for approximating
//! definite integrals of functions. It includes implementations of the Trapezoidal
//! rule, Simpson's rule, and an adaptive quadrature method.

use crate::numerical::elementary::eval_expr;
use crate::symbolic::core::Expr;
use std::collections::HashMap;

/// Enum to select the numerical integration method.
pub enum QuadratureMethod {
    Trapezoidal,
    Simpson,
    Adaptive,
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
pub fn trapezoidal_rule<F>(f: F, range: (f64, f64), n_steps: usize) -> f64
where
    F: Fn(f64) -> f64,
{
    let (a, b) = range;
    if a >= b {
        return 0.0;
    }
    let h = (b - a) / (n_steps as f64);
    let mut sum = 0.5 * (f(a) + f(b));
    for i in 1..n_steps {
        let x = a + (i as f64) * h;
        sum += f(x);
    }
    h * sum
}

/// # Pure Numerical Simpson's Rule
///
/// Performs numerical integration using Simpson's rule for a given function.
///
/// ## Arguments
/// * `f` - The function to integrate, as a closure.
/// * `range` - A tuple `(a, b)` representing the integration interval.
/// * `n_steps` - The number of steps to use (must be even).
///
/// ## Returns
/// The approximate value of the definite integral.
pub fn simpson_rule<F>(f: F, range: (f64, f64), n_steps: usize) -> Result<f64, String>
where
    F: Fn(f64) -> f64,
{
    let (a, b) = range;
    if a >= b {
        return Ok(0.0);
    }
    if n_steps % 2 != 0 {
        return Err("Simpson's rule requires an even number of steps.".to_string());
    }
    let h = (b - a) / (n_steps as f64);
    let mut sum = f(a) + f(b);
    for i in 1..n_steps {
        let x = a + (i as f64) * h;
        let factor = if i % 2 == 0 { 2.0 } else { 4.0 };
        sum += factor * f(x);
    }
    Ok((h / 3.0) * sum)
}

/// # Adaptive Quadrature
///
/// Performs numerical integration using an adaptive quadrature method based on Simpson's rule.
/// It recursively refines the integration subintervals to achieve a specified tolerance.
///
/// ## Arguments
/// * `f` - The function to integrate.
/// * `range` - The integration interval `(a, b)`.
/// * `tolerance` - The desired accuracy.
///
/// ## Returns
/// The approximate value of the definite integral.
pub fn adaptive_quadrature<F>(f: F, range: (f64, f64), tolerance: f64) -> f64
where
    F: Fn(f64) -> f64,
{
    fn adaptive_quadrature_recursive<F>(
        f: &F,
        range: (f64, f64),
        tolerance: f64,
        whole_integral: f64,
    ) -> f64
    where
        F: Fn(f64) -> f64,
    {
        let (a, b) = range;
        let mid = (a + b) / 2.0;
        let left_half = simpson_rule(f, (a, mid), 2).expect("Cannot apply rules to the left_half.");
        let right_half =
            simpson_rule(f, (mid, b), 2).expect("Cannot apply rules to the right_half.");
        let combined = left_half + right_half;

        if (combined - whole_integral).abs() <= 15.0 * tolerance {
            combined + (combined - whole_integral) / 15.0
        } else {
            let new_tolerance = tolerance / 2.0;
            adaptive_quadrature_recursive(f, (a, mid), new_tolerance, left_half)
                + adaptive_quadrature_recursive(f, (mid, b), new_tolerance, right_half)
        }
    }

    let (a, b) = range;
    let whole_integral =
        simpson_rule(&f, (a, b), 2).expect("Cannot apply rules to the whole_integral.");
    adaptive_quadrature_recursive(&f, range, tolerance, whole_integral)
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
    f: &Expr,
    var: &str,
    range: (f64, f64),
    n_steps: usize,
    method: &QuadratureMethod,
) -> Result<f64, String> {
    let func = |x: f64| -> f64 {
        let mut vars = HashMap::new();
        vars.insert(var.to_string(), x);
        eval_expr(f, &vars).unwrap_or(0.0)
    };

    match method {
        QuadratureMethod::Trapezoidal => Ok(trapezoidal_rule(func, range, n_steps)),
        QuadratureMethod::Simpson => simpson_rule(func, range, n_steps),
        QuadratureMethod::Adaptive => Ok(adaptive_quadrature(func, range, 1e-6)),
    }
}
