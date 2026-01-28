//! # Numerical Elementary Operations
//!
//! This module provides numerical evaluation and implementations of elementary
//! mathematical functions. It serves as the bridge between symbolic expressions
//! and their numerical values.
//!
//! ## Core Functionalities
//!
//! - **Expression Evaluation**: Recursively computes the numerical value of any
//!   symbolic `Expr` given a set of variable assignments.
//! - **Domain Validation**: Rigorously checks for mathematical domain violations
//!   (e.g., division by zero, logarithm of non-positive numbers).
//! - **Numerical Stability**: Handles edge cases like infinite values and NaN
//!   to ensure robust calculations.
//! - **Pure Numerical Functions**: High-performance wrappers around standard
//!   library functions for direct use on `f64`.
//!
//! ## Mathematical Formulas
//!
//! The evaluation follows standard mathematical definitions:
//!
//! - **Addition**: $f(x, y) = x + y$
//! - **Multiplication**: $f(x, y) = x \cdot y$
//! - **Power**: $f(x, y) = x^y$
//! - **Exponential**: $f(x) = e^x$
//! - **Logarithm**: $f(x) = \ln(x)$
//! - **Trigonometry**: $\sin(x), \cos(x), \tan(x)$, etc.

use std::collections::HashMap;

use num_traits::ToPrimitive;

use crate::symbolic::core::Expr;

/// Evaluates a symbolic expression to a numerical `f64` value.
///
/// This function recursively traverses the expression tree and computes the numerical value.
/// It comprehensively supports arithmetic, trigonometric, hyperbolic, and special constants.
///
/// # Arguments
///
/// * `expr` - The symbolic expression to evaluate.
/// * `vars` - A map containing numerical values for the variables in the expression.
///
/// # Returns
///
/// * `Ok(f64)` - The computed numerical result.
///
/// # Errors
///
/// Returns an error if:
/// - A variable in the expression is not provided in the `vars` map.
/// - A domain violation occurs (e.g., division by zero, log of non-positive number).
/// - An unsupported expression variant is encountered.
///
/// # Examples
///
/// ```rust
/// 
/// use std::collections::HashMap;
///
/// use rssn::numerical::elementary::eval_expr;
/// use rssn::symbolic::core::Expr;
///
/// let x = Expr::new_variable("x");
///
/// let expr = Expr::new_add(
///     x,
///     Expr::new_constant(2.0),
/// );
///
/// let mut vars = HashMap::new();
///
/// vars.insert("x".to_string(), 3.0);
///
/// let result = eval_expr(&expr, &vars).unwrap();
///
/// assert_eq!(result, 5.0);
/// ```

pub fn eval_expr<
    S: ::std::hash::BuildHasher,
>(
    expr: &Expr,
    vars: &HashMap<String, f64, S>,
) -> Result<f64, String> {

    match expr {
        | Expr::Dag(node) => {

            let converted_expr = node
                .to_expr()
                .map_err(|e| format!("Invalid DAG node: {e}"))?;

            eval_expr(
                &converted_expr,
                vars,
            )
        },
        | Expr::new_constant(c) => Ok(*c),
        | Expr::new_bigint(i) => {
            i.to_f64()
                .ok_or_else(|| "BigInt overflow during evaluation".to_string())
        },
        | Expr::new_rational(r) => {
            r.to_f64()
                .ok_or_else(|| "Rational overflow during evaluation".to_string())
        },
        | Expr::Variable(v) => {
            vars.get(v)
                .copied()
                .ok_or_else(|| format!("Unknown variable: '{v}'"))
        },

        // Arithmetic
        | Expr::Add(a, b) => Ok(eval_expr(a, vars)? + eval_expr(b, vars)?),
        | Expr::AddList(list) => {

            let mut sum = 0.0;

            for e in list {

                sum += eval_expr(e, vars)?;
            }

            Ok(sum)
        },
        | Expr::Sub(a, b) => Ok(eval_expr(a, vars)? - eval_expr(b, vars)?),
        | Expr::Mul(a, b) => Ok(eval_expr(a, vars)? * eval_expr(b, vars)?),
        | Expr::MulList(list) => {

            let mut prod = 1.0;

            for e in list {

                prod *= eval_expr(e, vars)?;
            }

            Ok(prod)
        },
        | Expr::Div(a, b) => {

            let den = eval_expr(b, vars)?;

            if den == 0.0 {

                return Err("Division by zero".to_string());
            }

            Ok(eval_expr(a, vars)? / den)
        },
        | Expr::Neg(a) => Ok(-eval_expr(a, vars)?),
        | Expr::Power(b, e) => {

            let base = eval_expr(b, vars)?;

            let exp = eval_expr(e, vars)?;

            if base == 0.0 && exp < 0.0 {

                return Err("Undefined operation: 0^negative power".to_string());
            }

            if base < 0.0 && exp.fract() != 0.0 {

                return Err(
                    "Complex result: negative base raised to non-integer power".to_string(),
                );
            }

            Ok(base.powf(exp))
        },
        | Expr::Abs(a) => Ok(eval_expr(a, vars)?.abs()),
        | Expr::Sqrt(a) => {

            let val = eval_expr(a, vars)?;

            if val < 0.0 {

                return Err("Square root of negative number".to_string());
            }

            Ok(val.sqrt())
        },

        // Trigonometric
        | Expr::Sin(a) => Ok(eval_expr(a, vars)?.sin()),
        | Expr::Cos(a) => Ok(eval_expr(a, vars)?.cos()),
        | Expr::Tan(a) => Ok(eval_expr(a, vars)?.tan()),
        | Expr::Sec(a) => Ok(1.0 / eval_expr(a, vars)?.cos()),
        | Expr::Csc(a) => Ok(1.0 / eval_expr(a, vars)?.sin()),
        | Expr::Cot(a) => Ok(1.0 / eval_expr(a, vars)?.tan()),

        // Inverse Trigonometric
        | Expr::ArcSin(a) => {

            let val = eval_expr(a, vars)?;

            if !(-1.0 ..= 1.0).contains(&val) {

                return Err("Inverse sine argument out of domain [-1, 1]".to_string());
            }

            Ok(val.asin())
        },
        | Expr::ArcCos(a) => {

            let val = eval_expr(a, vars)?;

            if !(-1.0 ..= 1.0).contains(&val) {

                return Err("Inverse cosine argument out of domain [-1, 1]".to_string());
            }

            Ok(val.acos())
        },
        | Expr::ArcTan(a) => Ok(eval_expr(a, vars)?.atan()),
        | Expr::Atan2(y, x) => Ok(eval_expr(y, vars)?.atan2(eval_expr(x, vars)?)),
        | Expr::ArcSec(a) => {

            let val = eval_expr(a, vars)?;

            if val.abs() < 1.0 {

                return Err(
                    "Inverse secant argument out of domain (-inf, -1] U [1, inf)".to_string(),
                );
            }

            Ok((1.0 / val).acos())
        },
        | Expr::ArcCsc(a) => {

            let val = eval_expr(a, vars)?;

            if val.abs() < 1.0 {

                return Err(
                    "Inverse cosecant argument out of domain (-inf, -1] U [1, inf)".to_string(),
                );
            }

            Ok((1.0 / val).asin())
        },
        | Expr::ArcCot(a) => Ok((1.0 / eval_expr(a, vars)?).atan()),

        // Hyperbolic
        | Expr::Sinh(a) => Ok(eval_expr(a, vars)?.sinh()),
        | Expr::Cosh(a) => Ok(eval_expr(a, vars)?.cosh()),
        | Expr::Tanh(a) => Ok(eval_expr(a, vars)?.tanh()),
        | Expr::Sech(a) => Ok(1.0 / eval_expr(a, vars)?.cosh()),
        | Expr::Csch(a) => Ok(1.0 / eval_expr(a, vars)?.sinh()),
        | Expr::Coth(a) => Ok(1.0 / eval_expr(a, vars)?.tanh()),

        // Inverse Hyperbolic
        | Expr::ArcSinh(a) => Ok(eval_expr(a, vars)?.asinh()),
        | Expr::ArcCosh(a) => {

            let val = eval_expr(a, vars)?;

            if val < 1.0 {

                return Err("Inverse hyperbolic cosine argument < 1".to_string());
            }

            Ok(val.acosh())
        },
        | Expr::ArcTanh(a) => {

            let val = eval_expr(a, vars)?;

            if val <= -1.0 || val >= 1.0 {

                return Err(
                    "Inverse hyperbolic tangent argument out of domain (-1, 1)".to_string(),
                );
            }

            Ok(val.atanh())
        },

        // Exponential and Logarithmic
        | Expr::Exp(a) => Ok(eval_expr(a, vars)?.exp()),
        | Expr::Log(a) => {

            let val = eval_expr(a, vars)?;

            if val <= 0.0 {

                return Err("Logarithm of non-positive number".to_string());
            }

            Ok(val.ln())
        },
        | Expr::LogBase(b, a) => {

            let base = eval_expr(b, vars)?;

            let val = eval_expr(a, vars)?;

            if base <= 0.0 || (base - 1.0).abs() < f64::EPSILON {

                return Err("Invalid logarithm base".to_string());
            }

            if val <= 0.0 {

                return Err("Logarithm of non-positive number".to_string());
            }

            Ok(val.log(base))
        },

        // Constants
        | Expr::Pi => Ok(std::f64::consts::PI),
        | Expr::E => Ok(std::f64::consts::E),
        | Expr::Infinity => Ok(f64::INFINITY),
        | Expr::NegativeInfinity => Ok(f64::NEG_INFINITY),

        // Rounding
        | Expr::Floor(a) => Ok(eval_expr(a, vars)?.floor()),

        // Comparison/Selection
        | Expr::Max(a, b) => Ok(eval_expr(a, vars)?.max(eval_expr(b, vars)?)),

        // Fallback or Unimplemented
        | _ => {
            Err(format!(
                "Numerical evaluation of {expr:?} is not supported"
            ))
        },
    }
}

/// Evaluates an expression with a single variable `x`.
///
/// # Errors
///
/// Returns an error if the expression evaluation fails.

pub fn eval_expr_single(
    expr: &Expr,
    x_name: &str,
    x_val: f64,
) -> Result<f64, String> {

    let mut vars = HashMap::new();

    vars.insert(
        x_name.to_string(),
        x_val,
    );

    eval_expr(expr, &vars)
}

/// # Pure Numerical Elementary Functions
///
/// Direct `f64` implementations of elementary functions for high-performance use.

pub mod pure {

    /// Sine function.
    #[must_use]

    pub fn sin(x: f64) -> f64 {

        x.sin()
    }

    /// Cosine function.
    #[must_use]

    pub fn cos(x: f64) -> f64 {

        x.cos()
    }

    /// Tangent function.
    #[must_use]

    pub fn tan(x: f64) -> f64 {

        x.tan()
    }

    /// Inverse sine.
    #[must_use]

    pub fn asin(x: f64) -> f64 {

        x.asin()
    }

    /// Inverse cosine.
    #[must_use]

    pub fn acos(x: f64) -> f64 {

        x.acos()
    }

    /// Inverse tangent.
    #[must_use]

    pub fn atan(x: f64) -> f64 {

        x.atan()
    }

    /// Two-argument inverse tangent.
    #[must_use]

    pub fn atan2(
        y: f64,
        x: f64,
    ) -> f64 {

        y.atan2(x)
    }

    /// Hyperbolic sine.
    #[must_use]

    pub fn sinh(x: f64) -> f64 {

        x.sinh()
    }

    /// Hyperbolic cosine.
    #[must_use]

    pub fn cosh(x: f64) -> f64 {

        x.cosh()
    }

    /// Hyperbolic tangent.
    #[must_use]

    pub fn tanh(x: f64) -> f64 {

        x.tanh()
    }

    /// Inverse hyperbolic sine.
    #[must_use]

    pub fn asinh(x: f64) -> f64 {

        x.asinh()
    }

    /// Inverse hyperbolic cosine.
    #[must_use]

    pub fn acosh(x: f64) -> f64 {

        x.acosh()
    }

    /// Inverse hyperbolic tangent.
    #[must_use]

    pub fn atanh(x: f64) -> f64 {

        x.atanh()
    }

    /// Absolute value.
    #[must_use]

    pub const fn abs(x: f64) -> f64 {

        x.abs()
    }

    /// Square root.
    #[must_use]

    pub fn sqrt(x: f64) -> f64 {

        x.sqrt()
    }

    /// Natural logarithm.
    #[must_use]

    pub fn ln(x: f64) -> f64 {

        x.ln()
    }

    /// Logarithm with base.
    #[must_use]

    pub fn log(
        x: f64,
        base: f64,
    ) -> f64 {

        x.log(base)
    }

    /// Exponential.
    #[must_use]

    pub fn exp(x: f64) -> f64 {

        x.exp()
    }

    /// Power.
    #[must_use]

    pub fn pow(
        base: f64,
        exp: f64,
    ) -> f64 {

        base.powf(exp)
    }

    /// Floor rounding.
    #[must_use]

    pub fn floor(x: f64) -> f64 {

        x.floor()
    }

    /// Ceil rounding.
    #[must_use]

    pub fn ceil(x: f64) -> f64 {

        x.ceil()
    }

    /// Round to nearest integer.
    #[must_use]

    pub fn round(x: f64) -> f64 {

        x.round()
    }

    /// Signum function.
    #[must_use]

    pub const fn signum(x: f64) -> f64 {

        x.signum()
    }
}
