//! # Numerical Complex Analysis
//!
//! This module provides numerical tools for complex analysis.
//! It includes functions for evaluating symbolic expressions to complex numbers,
//! which is fundamental for numerical computations involving complex functions.

use crate::symbolic::core::Expr;
use num_complex::Complex;
use num_traits::ToPrimitive;
use std::collections::HashMap;

/// Evaluates a symbolic expression to a numerical `Complex<f64>` value.
///
/// This function recursively traverses the expression tree and computes the complex numerical value.
/// It handles basic arithmetic, trigonometric, exponential, and logarithmic functions for complex numbers.
///
/// # Arguments
/// * `expr` - The expression to evaluate.
/// * `vars` - A `HashMap` containing the numerical `Complex<f64>` values for the variables in the expression.
///
/// # Returns
/// A `Result` containing the complex numerical value if the evaluation is successful, otherwise an error string.
pub fn eval_complex_expr(
    expr: &Expr,
    vars: &HashMap<String, Complex<f64>>,
) -> Result<Complex<f64>, String> {
    match expr {
        Expr::Constant(c) => Ok(Complex::new(*c, 0.0)),
        Expr::BigInt(i) => Ok(Complex::new(
            i.to_f64().ok_or("f64 conversion failed")?,
            0.0,
        )),
        Expr::Variable(v) => vars
            .get(v)
            .cloned()
            .ok_or_else(|| format!("Variable '{}' not found", v)),
        Expr::Complex(re, im) => {
            let re_val = eval_complex_expr(re, vars)?.re;
            let im_val = eval_complex_expr(im, vars)?.re; // eval_complex_expr returns Complex, so we take its real part
            Ok(Complex::new(re_val, im_val))
        }
        Expr::Add(a, b) => Ok(eval_complex_expr(a, vars)? + eval_complex_expr(b, vars)?),
        Expr::Sub(a, b) => Ok(eval_complex_expr(a, vars)? - eval_complex_expr(b, vars)?),
        Expr::Mul(a, b) => Ok(eval_complex_expr(a, vars)? * eval_complex_expr(b, vars)?),
        Expr::Div(a, b) => Ok(eval_complex_expr(a, vars)? / eval_complex_expr(b, vars)?),
        Expr::Power(b, e) => Ok(eval_complex_expr(b, vars)?.powc(eval_complex_expr(e, vars)?)),
        Expr::Neg(a) => Ok(-eval_complex_expr(a, vars)?),
        Expr::Sqrt(a) => Ok(eval_complex_expr(a, vars)?.sqrt()),
        Expr::Abs(a) => Ok(Complex::new(eval_complex_expr(a, vars)?.norm(), 0.0)),
        Expr::Sin(a) => Ok(eval_complex_expr(a, vars)?.sin()),
        Expr::Cos(a) => Ok(eval_complex_expr(a, vars)?.cos()),
        Expr::Tan(a) => Ok(eval_complex_expr(a, vars)?.tan()),
        Expr::Log(a) => Ok(eval_complex_expr(a, vars)?.ln()),
        Expr::Exp(a) => Ok(eval_complex_expr(a, vars)?.exp()),
        Expr::Pi => Ok(Complex::new(std::f64::consts::PI, 0.0)),
        Expr::E => Ok(Complex::new(std::f64::consts::E, 0.0)),
        _ => Err(format!(
            "Numerical complex evaluation for expression {:?} is not implemented",
            expr
        )),
    }
}
