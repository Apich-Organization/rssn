//! # Numerical Complex Analysis
//!
//! This module provides numerical tools for complex analysis.
//! It includes functions for evaluating symbolic expressions to complex numbers,
//! which is fundamental for numerical computations involving complex functions.
use crate::symbolic::core::Expr;
use num_complex::Complex;
use num_traits::{ToPrimitive, Zero};
use std::collections::HashMap;
/// # Numerical Contour Integration (Simpson's Rule)
///
/// This function computes the contour integral of a complex function `f` along a given path.
/// It uses Simpson's rule for improved accuracy over the trapezoidal rule.
///
/// ## Arguments
/// * `f` - The complex function to integrate, represented as a closure.
/// * `path` - A slice of `Complex<f64>` points defining the contour.
///
/// ## Returns
/// The numerical result of the contour integral.
pub fn contour_integral<F>(f: F, path: &[Complex<f64>]) -> Complex<f64>
where
    F: Fn(Complex<f64>) -> Complex<f64>,
{
    let mut integral = Complex::zero();
    for i in 0..path.len() - 1 {
        let z1 = path[i];
        let z2 = path[i + 1];
        let mid = (z1 + z2) / 2.0;
        let dz = z2 - z1;
        integral += (f(z1) + 4.0 * f(mid) + f(z2)) / 6.0 * dz;
    }
    integral
}
/// # Residue Calculation
///
/// This function calculates the residue of a complex function `f` at a point `z0`.
/// It uses a small circular contour around `z0` to compute the residue via Cauchy's Residue Theorem.
///
/// ## Arguments
/// * `f` - The complex function, represented as a closure.
/// * `z0` - The point at which to calculate the residue.
/// * `radius` - The radius of the circular contour to use for the calculation.
/// * `n_points` - The number of points to use for the contour integration.
///
/// ## Returns
/// The residue of the function at `z0`.
pub fn residue<F>(f: F, z0: Complex<f64>, radius: f64, n_points: usize) -> Complex<f64>
where
    F: Fn(Complex<f64>) -> Complex<f64>,
{
    let mut path = Vec::with_capacity(n_points + 1);
    for i in 0..=n_points {
        let angle = 2.0 * std::f64::consts::PI * (i as f64) / (n_points as f64);
        let point = z0 + radius * Complex::new(angle.cos(), angle.sin());
        path.push(point);
    }
    let integral = contour_integral(|z| f(z), &path);
    integral / (2.0 * std::f64::consts::PI * Complex::new(0.0, 1.0))
}
/// # Cauchy's Argument Principle
///
/// Calculates the number of zeros minus the number of poles (N - P) of a function `f`
/// inside a given contour. It does so by computing the winding number of f(z) around the origin.
///
/// ## Arguments
/// * `f` - The complex function.
/// * `contour` - A closed path in the complex plane.
///
/// ## Returns
/// The difference between the number of zeros and poles, which should be an integer.
pub fn count_zeros_poles<F>(f: F, contour: &[Complex<f64>]) -> Complex<f64>
where
    F: Fn(Complex<f64>) -> Complex<f64>,
{
    let integral = contour_integral(|z| complex_derivative(&f, z) / f(z), contour);
    integral / (2.0 * std::f64::consts::PI * Complex::new(0.0, 1.0))
}
/// # Numerical Differentiation
///
/// Computes the derivative of a complex function `f` at a point `z` using the central difference formula.
///
/// ## Arguments
/// * `f` - The complex function.
/// * `z` - The point at which to compute the derivative.
///
/// ## Returns
/// The numerical derivative of `f` at `z`.
pub fn complex_derivative<F>(f: &F, z: Complex<f64>) -> Complex<f64>
where
    F: Fn(Complex<f64>) -> Complex<f64>,
{
    let h = 1e-6;
    let h_complex = Complex::new(h, h);
    (f(z + h_complex) - f(z - h_complex)) / (2.0 * h_complex)
}
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
pub fn eval_complex_expr<S: ::std::hash::BuildHasher>(
    expr: &Expr,
    vars: &HashMap<String, Complex<f64>, S>,
) -> Result<Complex<f64>, String> {
    match expr {
		Expr::Dag(node) => {
			return eval_complex_expr(&node.to_expr().unwrap(), vars);
		}
        Expr::Constant(c) => Ok(Complex::new(*c, 0.0)),
        Expr::BigInt(i) => Ok(Complex::new(
            i.to_f64().ok_or("f64 conversion failed")?,
            0.0,
        )),
        Expr::Variable(v) => vars
            .get(v)
            .copied()
            .ok_or_else(|| format!("Variable '{}' not found", v)),
        Expr::Complex(re, im) => {
            let re_val = eval_complex_expr(re, vars)?.re;
            let im_val = eval_complex_expr(im, vars)?.re;
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
