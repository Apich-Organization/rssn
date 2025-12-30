//! # Numerical Evaluation
//!
//! This module provides functions for the numerical evaluation of symbolic expressions.
//! It attempts to simplify expressions to a single floating-point number where possible,
//! handling well-known constants and function values. It includes support for basic
//! arithmetic, trigonometric, hyperbolic, and special functions.

use std::f64::consts;

use num_complex::Complex64;
use num_traits::ToPrimitive;

use crate::symbolic::core::Expr;

const F64_EPSILON: f64 = 1e-9;

/// Evaluates a symbolic expression to a numerical `f64` value, if possible.
///
/// This function recursively traverses the expression tree and computes the numerical value.
/// It uses complex number evaluation internally to handle cases where real results
/// are obtained through complex intermediate steps (e.g., Cardano's formula for cubic roots).
///
/// # Arguments
/// * `expr` - The expression to evaluate.
///
/// # Returns
/// An `Option<f64>` containing the numerical value if the evaluation is successful
/// and the result is real (imaginary part near zero), otherwise `None`.
#[must_use]

pub fn evaluate_numerical(
    expr: &Expr
) -> Option<f64> {

    let complex_val =
        evaluate_complex(expr)?;

    if complex_val.im.abs()
        < F64_EPSILON
    {

        Some(complex_val.re)
    } else {

        None
    }
}

/// Evaluates a symbolic expression to a complex numerical `Complex64` value, if possible.
///
/// # Arguments
/// * `expr` - The expression to evaluate.
///
/// # Returns
/// An `Option<Complex64>` containing the complex numerical value if successful.
///
/// # Panics
///
/// Panics if a `Dag` node cannot be converted to an `Expr`, which indicates an
/// internal inconsistency in the expression representation. This should ideally
/// not happen in a well-formed expression DAG.
#[must_use]
#[allow(clippy::suboptimal_flops)]

pub fn evaluate_complex(
    expr: &Expr
) -> Option<Complex64> {

    match expr {
        | Expr::Dag(node) => {
            evaluate_complex(
                &node
                    .to_expr()
                    .expect(
                        "Eva Complex",
                    ),
            )
        },
        | Expr::Constant(c) => {
            Some(Complex64::new(
                *c, 0.0,
            ))
        },
        | Expr::BigInt(i) => {
            Some(Complex64::new(
                i.to_f64()?,
                0.0,
            ))
        },
        | Expr::Rational(r) => {
            Some(Complex64::new(
                r.to_f64()?,
                0.0,
            ))
        },
        | Expr::Pi => {
            Some(Complex64::new(
                consts::PI,
                0.0,
            ))
        },
        | Expr::E => {
            Some(Complex64::new(
                consts::E,
                0.0,
            ))
        },
        | Expr::Add(a, b) => {
            Some(
                evaluate_complex(a)?
                    + evaluate_complex(
                        b,
                    )?,
            )
        },
        | Expr::Sub(a, b) => {
            Some(
                evaluate_complex(a)?
                    - evaluate_complex(
                        b,
                    )?,
            )
        },
        | Expr::Mul(a, b) => {
            Some(
                evaluate_complex(a)?
                    * evaluate_complex(
                        b,
                    )?,
            )
        },
        | Expr::Div(a, b) => {
            Some(
                evaluate_complex(a)?
                    / evaluate_complex(
                        b,
                    )?,
            )
        },
        | Expr::Power(b, e) => {

            let base =
                evaluate_complex(b)?;

            let exp =
                evaluate_complex(e)?;

            if exp.im == 0.0 {

                Some(base.powf(exp.re))
            } else {

                Some(base.powc(exp))
            }
        },
        | Expr::Sqrt(a) => {
            Some(
                evaluate_complex(a)?
                    .sqrt(),
            )
        },
        | Expr::Log(a) => {
            Some(
                evaluate_complex(a)?
                    .ln(),
            )
        },
        | Expr::Exp(a) => {
            Some(
                evaluate_complex(a)?
                    .exp(),
            )
        },
        | Expr::Abs(a) => {
            Some(Complex64::new(
                evaluate_complex(a)?
                    .norm(),
                0.0,
            ))
        },
        | Expr::Sin(a) => {
            Some(
                evaluate_complex(a)?
                    .sin(),
            )
        },
        | Expr::Cos(a) => {
            Some(
                evaluate_complex(a)?
                    .cos(),
            )
        },
        | Expr::Tan(a) => {
            Some(
                evaluate_complex(a)?
                    .tan(),
            )
        },
        | Expr::ArcSin(a) => {
            Some(
                evaluate_complex(a)?
                    .asin(),
            )
        },
        | Expr::ArcCos(a) => {
            Some(
                evaluate_complex(a)?
                    .acos(),
            )
        },
        | Expr::ArcTan(a) => {
            Some(
                evaluate_complex(a)?
                    .atan(),
            )
        },
        | Expr::Sinh(a) => {
            Some(
                evaluate_complex(a)?
                    .sinh(),
            )
        },
        | Expr::Cosh(a) => {
            Some(
                evaluate_complex(a)?
                    .cosh(),
            )
        },
        | Expr::Tanh(a) => {
            Some(
                evaluate_complex(a)?
                    .tanh(),
            )
        },
        | Expr::ArcSinh(a) => {
            Some(
                evaluate_complex(a)?
                    .asinh(),
            )
        },
        | Expr::ArcCosh(a) => {
            Some(
                evaluate_complex(a)?
                    .acosh(),
            )
        },
        | Expr::ArcTanh(a) => {
            Some(
                evaluate_complex(a)?
                    .atanh(),
            )
        },
        | Expr::Complex(re, im) => {
            Some(Complex64::new(
                evaluate_complex(re)?
                    .re,
                evaluate_complex(im)?
                    .re,
            ))
        },
        | Expr::AddList(list) => {

            let mut sum =
                Complex64::new(
                    0.0, 0.0,
                );

            for item in list {

                sum +=
                    evaluate_complex(
                        item,
                    )?;
            }

            Some(sum)
        },
        | Expr::MulList(list) => {

            let mut prod =
                Complex64::new(
                    1.0, 0.0,
                );

            for item in list {

                prod *=
                    evaluate_complex(
                        item,
                    )?;
            }

            Some(prod)
        },
        | Expr::LogBase(base, arg) => {

            let b =
                evaluate_complex(base)?;

            let a =
                evaluate_complex(arg)?;

            Some(a.ln() / b.ln())
        },
        | Expr::Factorial(a) => {

            let val =
                evaluate_complex(a)?;

            if val.im == 0.0
                && val.re.fract() == 0.0
                && val.re >= 0.0
            {

                let mut result = 1.0;

                for i in
                    2 ..= (val.re as i64).try_into().unwrap_or(0)
                {

                    result *= f64::from(i);
                }

                Some(Complex64::new(
                    result,
                    0.0,
                ))
            } else {

                None
            }
        },
        | Expr::Floor(a) => {

            let val =
                evaluate_complex(a)?;

            if val.im == 0.0 {

                Some(Complex64::new(
                    val.re.floor(),
                    0.0,
                ))
            } else {

                None
            }
        },
        | Expr::Neg(a) => {
            Some(-evaluate_complex(
                a,
            )?)
        },
        | _ => None,
    }
}
