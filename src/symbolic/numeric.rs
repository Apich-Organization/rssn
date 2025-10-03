//! # Numerical Evaluation
//!
//! This module provides functions for the numerical evaluation of symbolic expressions.
//! It attempts to simplify expressions to a single floating-point number where possible,
//! handling well-known constants and function values. It includes support for basic
//! arithmetic, trigonometric, hyperbolic, and special functions.

use crate::symbolic::core::Expr;
use num_traits::ToPrimitive;
use std::f64::consts;

const F64_EPSILON: f64 = 1e-9;

/// Evaluates a symbolic expression to a numerical `f64` value, if possible.
///
/// This function recursively traverses the expression tree and computes the numerical value.
/// It handles:
/// - Constants and numeric types (`Constant`, `BigInt`, `Rational`).
/// - Symbolic constants (`Pi`, `E`).
/// - Basic arithmetic operations (`Add`, `Sub`, `Mul`, `Div`, `Power`).
/// - A wide range of elementary and special functions.
///
/// For trigonometric functions, it specifically checks for common angles (e.g., π/6, π/4, π/2)
/// to provide exact results instead of floating-point approximations.
///
/// # Arguments
/// * `expr` - The expression to evaluate.
///
/// # Returns
/// An `Option<f64>` containing the numerical value if the evaluation is successful,
/// otherwise `None`.
pub fn evaluate_numerical(expr: &Expr) -> Option<f64> {
    match expr {
        Expr::Constant(c) => Some(*c),
        Expr::BigInt(i) => i.to_f64(),
        Expr::Rational(r) => r.to_f64(),
        Expr::Pi => Some(consts::PI),
        Expr::E => Some(consts::E),
        Expr::Add(a, b) => Some(evaluate_numerical(a)? + evaluate_numerical(b)?),
        Expr::Sub(a, b) => Some(evaluate_numerical(a)? - evaluate_numerical(b)?),
        Expr::Mul(a, b) => Some(evaluate_numerical(a)? * evaluate_numerical(b)?),
        Expr::Div(a, b) => Some(evaluate_numerical(a)? / evaluate_numerical(b)?),
        Expr::Power(b, e) => Some(evaluate_numerical(b)?.powf(evaluate_numerical(e)?)),
        Expr::Sqrt(a) => Some(evaluate_numerical(a)?.sqrt()),
        Expr::Log(a) => Some(evaluate_numerical(a)?.ln()),
        Expr::Exp(a) => Some(evaluate_numerical(a)?.exp()),
        Expr::Abs(a) => Some(evaluate_numerical(a)?.abs()),

        // Trigonometric Functions
        Expr::Sin(a) => {
            let val = evaluate_numerical(a)?;
            if (val.abs() - consts::PI).abs() < F64_EPSILON {
                Some(0.0)
            } else if (val - consts::FRAC_PI_6).abs() < F64_EPSILON {
                Some(0.5)
            } else if (val - consts::FRAC_PI_4).abs() < F64_EPSILON {
                Some(consts::FRAC_1_SQRT_2)
            } else if (val - consts::FRAC_PI_3).abs() < F64_EPSILON {
                Some(3.0f64.sqrt() / 2.0)
            } else if (val - consts::FRAC_PI_2).abs() < F64_EPSILON {
                Some(1.0)
            } else {
                Some(val.sin())
            }
        }
        Expr::Cos(a) => {
            let val = evaluate_numerical(a)?;
            if (val.abs() - consts::FRAC_PI_2).abs() < F64_EPSILON {
                Some(0.0)
            } else if (val - consts::FRAC_PI_6).abs() < F64_EPSILON {
                Some(3.0f64.sqrt() / 2.0)
            } else if (val - consts::FRAC_PI_4).abs() < F64_EPSILON {
                Some(consts::FRAC_1_SQRT_2)
            } else if (val - consts::FRAC_PI_3).abs() < F64_EPSILON {
                Some(0.5)
            } else if (val.abs() - consts::PI).abs() < F64_EPSILON {
                Some(-1.0)
            } else {
                Some(val.cos())
            }
        }
        Expr::Tan(a) => {
            let val = evaluate_numerical(a)?;
            if val.abs() < F64_EPSILON {
                Some(0.0)
            } else if (val - consts::FRAC_PI_4).abs() < F64_EPSILON {
                Some(1.0)
            } else if (val.abs() - consts::FRAC_PI_2).abs() < F64_EPSILON {
                Some(f64::INFINITY)
            } else {
                Some(val.tan())
            }
        }
        Expr::ArcSin(a) => Some(evaluate_numerical(a)?.asin()),
        Expr::ArcCos(a) => Some(evaluate_numerical(a)?.acos()),
        Expr::ArcTan(a) => Some(evaluate_numerical(a)?.atan()),

        // Hyperbolic Functions
        Expr::Sinh(a) => Some(evaluate_numerical(a)?.sinh()),
        Expr::Cosh(a) => Some(evaluate_numerical(a)?.cosh()),
        Expr::Tanh(a) => Some(evaluate_numerical(a)?.tanh()),
        Expr::ArcSinh(a) => Some(evaluate_numerical(a)?.asinh()),
        Expr::ArcCosh(a) => Some(evaluate_numerical(a)?.acosh()),
        Expr::ArcTanh(a) => Some(evaluate_numerical(a)?.atanh()),

        // Factorial
        Expr::Factorial(a) => {
            let n = evaluate_numerical(a)?;
            if n.fract() == 0.0 && n >= 0.0 {
                let mut result = 1.0;
                for i in 2..=n as u64 {
                    result *= i as f64;
                }
                Some(result)
            } else {
                // For non-integers, this would be the Gamma function Γ(n+1)
                None
            }
        }

        // Other special functions - delegate to external crates if available, otherwise None
        Expr::Gamma(a) => {
            // In a real implementation, you would link to a gamma function library like `libm` or `statrs`
            let _val = evaluate_numerical(a)?;
            // Placeholder for gamma function
            None
        }

        Expr::Floor(a) => Some(evaluate_numerical(a)?.floor()),
        Expr::Neg(a) => Some(-evaluate_numerical(a)?),

        _ => None, // Return None for expressions that cannot be evaluated to a single number
    }
}
