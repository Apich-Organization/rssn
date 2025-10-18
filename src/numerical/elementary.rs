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
    match expr {
		Expr::Dag(node) => {
			return eval_expr(&node.to_expr().unwrap(), vars);
		}
        Expr::Constant(c) => Ok(*c),
        Expr::BigInt(i) => Ok(i
            .to_f64()
            .ok_or_else(|| "BigInt conversion to f64 failed".to_string())?),
        Expr::Variable(v) => vars
            .get(v)
            .copied()
            .ok_or_else(|| format!("Variable '{}' not found", v)),
        Expr::Add(a, b) => Ok(eval_expr(a, vars)? + eval_expr(b, vars)?),
        Expr::Sub(a, b) => Ok(eval_expr(a, vars)? - eval_expr(b, vars)?),
        Expr::Mul(a, b) => Ok(eval_expr(a, vars)? * eval_expr(b, vars)?),
        Expr::Div(a, b) => Ok(eval_expr(a, vars)? / eval_expr(b, vars)?),
        Expr::Power(b, e) => Ok(eval_expr(b, vars)?.powf(eval_expr(e, vars)?)),
        Expr::Neg(a) => Ok(-eval_expr(a, vars)?),
        Expr::Sqrt(a) => Ok(eval_expr(a, vars)?.sqrt()),
        Expr::Abs(a) => Ok(eval_expr(a, vars)?.abs()),
        Expr::Sin(a) => Ok(eval_expr(a, vars)?.sin()),
        Expr::Cos(a) => Ok(eval_expr(a, vars)?.cos()),
        Expr::Tan(a) => Ok(eval_expr(a, vars)?.tan()),
        Expr::Log(a) => Ok(eval_expr(a, vars)?.ln()),
        Expr::Exp(a) => Ok(eval_expr(a, vars)?.exp()),
        Expr::Pi => Ok(std::f64::consts::PI),
        Expr::E => Ok(std::f64::consts::E),
        _ => Err(format!(
            "Numerical evaluation for expression {:?} is not implemented",
            expr
        )),
    }
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
