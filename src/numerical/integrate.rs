//! # Numerical Integration (Quadrature)
//!
//! This module provides numerical integration (quadrature) methods for approximating
//! definite integrals of functions. It includes implementations of the Trapezoidal
//! rule and Simpson's rule, which are widely used for their accuracy and efficiency.

use crate::symbolic::core::Expr;
//use crate::numerical::elementary::eval_expr as other_eval_expr;
use num_traits::ToPrimitive;
use std::collections::HashMap;

/// Enum to select the numerical integration method.
pub enum QuadratureMethod {
    Trapezoidal,
    Simpson,
}

/// Performs numerical integration (quadrature) of a function `f(x)` over an interval `[a, b]`.
///
/// # Arguments
/// * `f` - The expression to integrate.
/// * `var` - The variable of integration.
/// * `range` - A tuple `(a, b)` representing the integration interval.
/// * `n_steps` - The number of steps to use for the integration.
/// * `method` - The quadrature method to use.
///
/// # Returns
/// A `Result` containing the numerical value of the integral, or an error string.
pub fn quadrature(
    f: &Expr,
    var: &str,
    range: (f64, f64),
    n_steps: usize,
    method: QuadratureMethod,
) -> Result<f64, String> {
    let (a, b) = range;
    if a >= b {
        return Ok(0.0);
    }

    let h = (b - a) / (n_steps as f64);
    let mut vars = HashMap::new();

    match method {
        QuadratureMethod::Trapezoidal => {
            let mut sum = 0.5 * (eval_at(f, var, a, &mut vars)? + eval_at(f, var, b, &mut vars)?);
            for i in 1..n_steps {
                let x = a + (i as f64) * h;
                sum += eval_at(f, var, x, &mut vars)?;
            }
            Ok(h * sum)
        }
        QuadratureMethod::Simpson => {
            if n_steps % 2 != 0 {
                return Err("Simpson's rule requires an even number of steps.".to_string());
            }
            let mut sum = eval_at(f, var, a, &mut vars)? + eval_at(f, var, b, &mut vars)?;
            for i in 1..n_steps {
                let x = a + (i as f64) * h;
                let factor = if i % 2 == 0 { 2.0 } else { 4.0 };
                sum += factor * eval_at(f, var, x, &mut vars)?;
            }
            Ok((h / 3.0) * sum)
        }
    }
}

/// Helper to evaluate an expression at a single point.
pub(crate) fn eval_at(
    expr: &Expr,
    var: &str,
    val: f64,
    vars: &mut HashMap<String, f64>,
) -> Result<f64, String> {
    vars.insert(var.to_string(), val);
    eval_expr(expr, vars)
}

/// Recursive helper to evaluate an expression to a numerical f64 value.
pub(crate) fn eval_expr(expr: &Expr, vars: &HashMap<String, f64>) -> Result<f64, String> {
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
