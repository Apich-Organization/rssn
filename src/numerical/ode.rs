//! # Numerical Ordinary Differential Equation (ODE) Solvers
//!
//! This module provides numerical methods for solving ordinary differential equations (ODEs).
//! It includes an implementation of the fourth-order Runge-Kutta (RK4) method for solving
//! systems of first-order ODEs, which is a widely used and accurate technique.

use crate::numerical::elementary::eval_expr;
use crate::symbolic::core::Expr;
use std::collections::HashMap;

/// Solves a system of first-order ODEs `Y' = F(x, Y)` using the fourth-order Runge-Kutta method.
///
/// The RK4 method is a popular and accurate numerical technique for approximating solutions
/// to ordinary differential equations. It involves calculating four intermediate slopes
/// to estimate the next point in the solution.
///
/// # Arguments
/// * `funcs` - A slice of expressions for the functions `F_i(x, y_0, y_1, ...)`.
/// * `y0` - The initial values `[y0_0, y0_1, ...]`.
/// * `x_range` - The interval `(x0, x_end)` over which to solve.
/// * `num_steps` - The number of steps to take.
///
/// # Returns
/// A `Result` containing a vector of states, where each state is a vector of `y` values
/// at a given `x`, or an error string if evaluation fails.
pub fn solve_ode_system_rk4(
    funcs: &[Expr],
    y0: &[f64],
    x_range: (f64, f64),
    num_steps: usize,
) -> Result<Vec<Vec<f64>>, String> {
    let (x0, x_end) = x_range;
    let h = (x_end - x0) / (num_steps as f64);
    let mut x = x0;
    let mut y_vec = y0.to_vec();
    let mut results = vec![y_vec.clone()];
    let mut vars = HashMap::new();

    for _ in 0..num_steps {
        let k1 = eval_f(funcs, x, &y_vec, &mut vars)?;
        let k2 = eval_f(
            funcs,
            x + h / 2.0,
            &add_vec(&y_vec, &scale_vec(&k1, h / 2.0)),
            &mut vars,
        )?;
        let k3 = eval_f(
            funcs,
            x + h / 2.0,
            &add_vec(&y_vec, &scale_vec(&k2, h / 2.0)),
            &mut vars,
        )?;
        let k4 = eval_f(
            funcs,
            x + h,
            &add_vec(&y_vec, &scale_vec(&k3, h)),
            &mut vars,
        )?;

        let weighted_sum = add_vec(&scale_vec(&k2, 2.0), &scale_vec(&k3, 2.0));
        let weighted_sum = add_vec(&weighted_sum, &k1);
        let weighted_sum = add_vec(&weighted_sum, &k4);

        y_vec = add_vec(&y_vec, &scale_vec(&weighted_sum, h / 6.0));
        x += h;
        results.push(y_vec.clone());
    }

    Ok(results)
}

pub(crate) fn eval_f(
    funcs: &[Expr],
    x: f64,
    y_vec: &[f64],
    vars: &mut HashMap<String, f64>,
) -> Result<Vec<f64>, String> {
    vars.insert("x".to_string(), x);
    for (i, y_val) in y_vec.iter().enumerate() {
        vars.insert(format!("y{}", i), *y_val);
    }
    let mut results = Vec::new();
    for f in funcs {
        results.push(eval_expr(f, vars)?);
    }
    Ok(results)
}

pub(crate) fn add_vec(v1: &[f64], v2: &[f64]) -> Vec<f64> {
    v1.iter().zip(v2.iter()).map(|(a, b)| a + b).collect()
}

pub(crate) fn scale_vec(v: &[f64], s: f64) -> Vec<f64> {
    v.iter().map(|a| a * s).collect()
}
