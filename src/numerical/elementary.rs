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
pub fn eval_expr(expr: &Expr, vars: &HashMap<String, f64>) -> Result<f64, String> {
    match expr {
        Expr::Constant(c) => Ok(*c),
        Expr::BigInt(i) => Ok(i
            .to_f64()
            .ok_or_else(|| "BigInt conversion to f64 failed".to_string())?),
        Expr::Variable(v) => vars
            .get(v)
            .cloned()
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
