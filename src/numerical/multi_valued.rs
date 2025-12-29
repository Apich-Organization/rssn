//! # Numerical Multi-valued Functions in Complex Analysis
//!
//! This module provides numerical methods for handling multi-valued functions
//! in complex analysis, particularly focusing on finding roots of complex functions
//! using Newton's method.

use std::collections::HashMap;

use num_complex::Complex;

use crate::numerical::complex_analysis::eval_complex_expr;
use crate::symbolic::core::Expr;

/// Finds a root of a complex function `f(z) = 0` using Newton's method.
///
/// Newton's method is an iterative root-finding algorithm. For complex functions,
/// it uses the formula `z_{n+1} = z_n - f(z_n) / f'(z_n)`.
///
/// # Arguments
/// * `f` - The complex function as a symbolic expression.
/// * `f_prime` - The derivative of the function, `f'`.
/// * `start_point` - An initial guess for the root in the complex plane.
/// * `tolerance` - The desired precision of the root.
/// * `max_iter` - The maximum number of iterations.
///
/// # Returns
/// An `Option` containing the complex root if found, otherwise `None`.
#[must_use]

pub fn newton_method_complex(
    f: &Expr,
    f_prime: &Expr,
    start_point: Complex<f64>,
    tolerance: f64,
    max_iter: usize,
) -> Option<Complex<f64>> {

    let mut z = start_point;

    let mut vars = HashMap::new();

    for _ in 0 .. max_iter {

        vars.insert("z".to_string(), z);

        let f_val =
            match eval_complex_expr(
                f, &vars,
            ) {
                | Ok(val) => val,
                | Err(_) => {
                    return None
                },
            };

        let f_prime_val =
            match eval_complex_expr(
                f_prime,
                &vars,
            ) {
                | Ok(val) => val,
                | Err(_) => {
                    return None
                },
            };

        if f_prime_val.norm_sqr()
            < 1e-12
        {

            return None;
        }

        let delta = f_val / f_prime_val;

        z -= delta;

        if delta.norm() < tolerance {

            return Some(z);
        }
    }

    None
}

/// Computes the k-th branch of the complex logarithm.

#[must_use]

pub fn complex_log_k(
    z: Complex<f64>,
    k: i32,
) -> Complex<f64> {

    let ln_r = z.norm().ln();

    let theta = (2.0
        * std::f64::consts::PI)
        .mul_add(
            f64::from(k),
            z.arg(),
        );

    Complex::new(ln_r, theta)
}

/// Computes the k-th branch of the complex square root.

#[must_use]

pub fn complex_sqrt_k(
    z: Complex<f64>,
    k: i32,
) -> Complex<f64> {

    let r_sqrt = z.norm().sqrt();

    let theta = (2.0
        * std::f64::consts::PI)
        .mul_add(
            f64::from(k),
            z.arg(),
        )
        / 2.0;

    Complex::new(
        r_sqrt * theta.cos(),
        r_sqrt * theta.sin(),
    )
}

/// Computes the k-th branch of the complex power z^w.

#[must_use]

pub fn complex_pow_k(
    z: Complex<f64>,
    w: Complex<f64>,
    k: i32,
) -> Complex<f64> {

    let log_z_k = complex_log_k(z, k);

    (w * log_z_k).exp()
}

/// Computes the k-th branch of the complex n-th root.

#[must_use]

pub fn complex_nth_root_k(
    z: Complex<f64>,
    n: u32,
    k: i32,
) -> Complex<f64> {

    let r_root = z
        .norm()
        .powf(1.0 / f64::from(n));

    let theta = (2.0
        * std::f64::consts::PI)
        .mul_add(
            f64::from(k),
            z.arg(),
        )
        / f64::from(n);

    Complex::new(
        r_root * theta.cos(),
        r_root * theta.sin(),
    )
}

/// Computes the k-th branch of the complex arcsine.

#[must_use]

pub fn complex_arcsin_k(
    z: Complex<f64>,
    k: i32,
) -> Complex<f64> {

    let pi = std::f64::consts::PI;

    let principal = z.asin();

    let k_pi = Complex::new(
        f64::from(k) * pi,
        0.0,
    );

    if k % 2 == 0 {

        k_pi + principal
    } else {

        k_pi - principal
    }
}

/// Computes the k-th branch of the complex arccosine.
/// s is +1 or -1.

#[must_use]

pub fn complex_arccos_k(
    z: Complex<f64>,
    k: i32,
    s: i32,
) -> Complex<f64> {

    let pi = std::f64::consts::PI;

    let principal = z.acos();

    let two_k_pi = Complex::new(
        2.0 * f64::from(k) * pi,
        0.0,
    );

    if s >= 0 {

        two_k_pi + principal
    } else {

        two_k_pi - principal
    }
}

/// Computes the k-th branch of the complex arctangent.

#[must_use]

pub fn complex_arctan_k(
    z: Complex<f64>,
    k: i32,
) -> Complex<f64> {

    let pi = std::f64::consts::PI;

    let principal = z.atan();

    let k_pi = Complex::new(
        f64::from(k) * pi,
        0.0,
    );

    k_pi + principal
}
