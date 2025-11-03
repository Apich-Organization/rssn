//! # Numerical Elementary Operations
//!
//! This module provides numerical evaluation of symbolic expressions.
//! It includes a core function `eval_expr` that recursively evaluates an `Expr`
//! to an `f64` value, handling basic arithmetic, trigonometric, and exponential functions.
use crate::symbolic::core::Expr;
use num_traits::ToPrimitive;
use std::collections::HashMap;
/// Evaluates a symbolic expression to a numerical `f64` value.
///
/// This function recursively traverses the expression tree and computes the numerical value.
/// It handles basic arithmetic, trigonometric, and exponential functions.
///
/// # Arguments
/// * `expr` - The expression to evaluate.
/// * `vars` - A `HashMap` containing the numerical `f64` values for the variables in the expression.
///
/// # Returns
/// A `Result` containing the numerical value if the evaluation is successful, otherwise an error string.
pub fn eval_expr<S: ::std::hash::BuildHasher>(
    expr: &Expr,
    vars: &HashMap<String, f64, S>,
) -> Result<f64, String> {
    let value = match expr {
        Expr::Dag(node) => {
            let converted_expr = node
                .to_expr()
                .map_err(|e| format!("Failed to convert DAG node to Expr: {}", e))?;
            eval_expr(&converted_expr, vars)?
        }
        Expr::Constant(c) => *c,
        Expr::BigInt(i) => i
            .to_f64()
            .ok_or_else(|| "BigInt conversion to f64 failed".to_string())?,
        Expr::Variable(v) => vars
            .get(v)
            .copied()
            .ok_or_else(|| format!("Variable '{}' not found", v))?,
        Expr::Add(a, b) => {
            let val_a = eval_expr(a, vars)?;
            let val_b = eval_expr(b, vars)?;
            let result = val_a + val_b;
            if result.is_infinite() {
                return Err("Addition resulted in infinite value".to_string());
            }
            if result.is_nan() {
                return Err("Addition resulted in NaN".to_string());
            }
            result
        }
        Expr::Sub(a, b) => {
            let val_a = eval_expr(a, vars)?;
            let val_b = eval_expr(b, vars)?;
            let result = val_a - val_b;
            if result.is_infinite() {
                return Err("Subtraction resulted in infinite value".to_string());
            }
            if result.is_nan() {
                return Err("Subtraction resulted in NaN".to_string());
            }
            result
        }
        Expr::Mul(a, b) => {
            let val_a = eval_expr(a, vars)?;
            let val_b = eval_expr(b, vars)?;
            let result = val_a * val_b;
            if result.is_infinite()
                && !(val_a.is_infinite() || val_b.is_infinite())
                && !(val_a == 0.0 || val_b == 0.0)
            {
                // Only report infinity if it's not expected from infinite operands or 0*inf cases
                return Err("Multiplication resulted in infinite value".to_string());
            }
            if result.is_nan() {
                return Err("Multiplication resulted in NaN".to_string());
            }
            result
        }
        Expr::Div(a, b) => {
            let val_a = eval_expr(a, vars)?;
            let val_b = eval_expr(b, vars)?;
            if val_b == 0.0 {
                if val_a == 0.0 {
                    return Err("Division by zero: 0/0 is undefined".to_string());
                } else {
                    return Err("Division by zero".to_string());
                }
            }
            let result = val_a / val_b;
            if result.is_infinite() {
                return Err("Division resulted in infinite value".to_string());
            }
            if result.is_nan() {
                return Err("Division resulted in NaN".to_string());
            }
            result
        }
        Expr::Power(base, exp) => {
            let base_val = eval_expr(base, vars)?;
            let exp_val = eval_expr(exp, vars)?;
            // Check for common problematic cases
            if base_val == 0.0 && exp_val < 0.0 {
                return Err("Power operation: 0 raised to negative power is undefined".to_string());
            }
            if base_val < 0.0 && exp_val.fract() != 0.0 {
                // Non-integer power of negative number
                return Err(
                    "Power operation: negative base raised to non-integer exponent is not real"
                        .to_string(),
                );
            }
            let result = base_val.powf(exp_val);
            if result.is_infinite() {
                return Err("Power operation resulted in infinite value".to_string());
            }
            if result.is_nan() {
                return Err("Power operation resulted in NaN".to_string());
            }
            result
        }
        Expr::Neg(a) => {
            let val = eval_expr(a, vars)?;
            let result = -val;
            if result.is_infinite() {
                return Err("Negation resulted in infinite value".to_string());
            }
            result
        }
        Expr::Sqrt(a) => {
            let val = eval_expr(a, vars)?;
            if val < 0.0 {
                return Err("Square root of negative number".to_string());
            }
            let result = val.sqrt();
            if result.is_infinite() || result.is_nan() {
                return Err("Square root resulted in invalid value".to_string());
            }
            result
        }
        Expr::Abs(a) => {
            let val = eval_expr(a, vars)?;
            let result = val.abs();
            if result.is_infinite() || result.is_nan() {
                return Err("Absolute value resulted in invalid value".to_string());
            }
            result
        }
        Expr::Sin(a) => {
            let val = eval_expr(a, vars)?;
            let result = val.sin();
            if result.is_infinite() || result.is_nan() {
                return Err("Sine operation resulted in invalid value".to_string());
            }
            result
        }
        Expr::Cos(a) => {
            let val = eval_expr(a, vars)?;
            let result = val.cos();
            if result.is_infinite() || result.is_nan() {
                return Err("Cosine operation resulted in invalid value".to_string());
            }
            result
        }
        Expr::Tan(a) => {
            let val = eval_expr(a, vars)?;
            let result = val.tan(); // May produce infinity for values near (2k+1)*π/2
            if result.is_nan() {
                return Err("Tangent operation resulted in invalid value".to_string());
            }
            result
        }
        Expr::Log(a) => {
            let val = eval_expr(a, vars)?;
            if val <= 0.0 {
                return Err("Logarithm of non-positive number".to_string());
            }
            let result = val.ln();
            if result.is_infinite() || result.is_nan() {
                return Err("Logarithm resulted in invalid value".to_string());
            }
            result
        }
        Expr::Exp(a) => {
            let val = eval_expr(a, vars)?;
            let result = val.exp();
            if result.is_infinite() && val.is_finite() {
                // Only report if it's not expected from an infinite input
                return Err("Exponential resulted in infinite value".to_string());
            }
            if result.is_nan() {
                return Err("Exponential resulted in NaN".to_string());
            }
            result
        }
        Expr::Pi => std::f64::consts::PI,
        Expr::E => std::f64::consts::E,
        _ => {
            return Err(format!(
                "Numerical evaluation for expression {:?} is not implemented",
                expr
            ))
        }
    };

    // Final check for validity of the result
    if value.is_nan() {
        return Err("Evaluation resulted in NaN".to_string());
    }

    Ok(value)
}
/// # Pure Numerical Elementary Functions
///
/// This module provides pure numerical implementations of elementary mathematical functions.
/// These functions operate directly on `f64` values and are used for numerical calculations
/// without any symbolic manipulation.
pub mod pure {
    /// Computes the sine of a number.
    pub fn sin(x: f64) -> f64 {
        x.sin()
    }
    /// Computes the cosine of a number.
    pub fn cos(x: f64) -> f64 {
        x.cos()
    }
    /// Computes the tangent of a number.
    pub fn tan(x: f64) -> f64 {
        x.tan()
    }
    /// Computes the absolute value of a number.
    pub fn abs(x: f64) -> f64 {
        x.abs()
    }
    /// Computes the natural logarithm of a number.
    pub fn ln(x: f64) -> f64 {
        x.ln()
    }
    /// Computes the exponential of a number.
    pub fn exp(x: f64) -> f64 {
        x.exp()
    }
    /// Computes the square root of a number.
    pub fn sqrt(x: f64) -> f64 {
        x.sqrt()
    }
    /// Raises a number to a power.
    pub fn pow(base: f64, exp: f64) -> f64 {
        base.powf(exp)
    }
}
