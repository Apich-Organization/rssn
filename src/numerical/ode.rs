use std::collections::HashMap;

use serde::Deserialize;
use serde::Serialize;

use crate::numerical::elementary::eval_expr;
use crate::symbolic::core::Expr;

/// Methods for solving ordinary differential equations.
#[derive(
    Debug,
    Clone,
    Copy,
    Serialize,
    Deserialize,
    PartialEq,
    Eq,
)]

pub enum OdeSolverMethod {
    /// Forward Euler method (1st order).
    Euler,
    /// Heun's method (2nd order Runge-Kutta).
    Heun,
    /// Classic RK4 method (4th order Runge-Kutta).
    RungeKutta4,
}

/// Solves a system of first-order ODEs `Y' = F(x, Y)`.
///
/// This is a general wrapper for different numerical ODE integration methods.
///
/// # Arguments
/// * `funcs` - Expressions for the derivatives `dy_i/dx = F_i(x, y_0, y_1, ...)`.
/// * `y0` - Initial values `[y0_0, y0_1, ...]`.
/// * `x_range` - The interval `(x0, x_end)`.
/// * `num_steps` - Total number of integration steps.
/// * `method` - The solver method to use.
///
/// # Example
/// ```rust
/// 
/// use rssn::numerical::ode::solve_ode_system;
/// use rssn::numerical::ode::OdeSolverMethod;
/// use rssn::symbolic::core::Expr;
///
/// let y0 = Expr::new_variable("y0"); // dy/dx = y
/// let y_init = vec![1.0];
///
/// let res = solve_ode_system(
///     &[y0],
///     &y_init,
///     (0.0, 1.0),
///     10,
///     OdeSolverMethod::RungeKutta4,
/// )
/// .unwrap();
///
/// // At x=1, y should be approx e^1 = 2.718...
/// assert!((res.last().unwrap()[0] - 2.718).abs() < 0.01);
/// ```

pub fn solve_ode_system(
    funcs: &[Expr],
    y0: &[f64],
    x_range: (f64, f64),
    num_steps: usize,
    method: OdeSolverMethod,
) -> Result<Vec<Vec<f64>>, String> {

    match method {
        | OdeSolverMethod::Euler => {
            solve_ode_euler(
                funcs,
                y0,
                x_range,
                num_steps,
            )
        },
        | OdeSolverMethod::Heun => {
            solve_ode_heun(
                funcs,
                y0,
                x_range,
                num_steps,
            )
        },
        | OdeSolverMethod::RungeKutta4 => {
            solve_ode_system_rk4(
                funcs,
                y0,
                x_range,
                num_steps,
            )
        },
    }
}

/// Solves an ODE system using the Euler method (first-order).

pub fn solve_ode_euler(
    funcs: &[Expr],
    y0: &[f64],
    x_range: (f64, f64),
    num_steps: usize,
) -> Result<Vec<Vec<f64>>, String> {

    let (x0, x_end) = x_range;

    let h = (x_end - x0)
        / (num_steps as f64);

    let mut x = x0;

    let mut y = y0.to_vec();

    let mut results = vec![y.clone()];

    let mut vars = HashMap::new();

    for _ in 0 .. num_steps {

        let dy = eval_f(
            funcs,
            x,
            &y,
            &mut vars,
        )?;

        y = add_vec(
            &y,
            &scale_vec(&dy, h),
        );

        x += h;

        results.push(y.clone());
    }

    Ok(results)
}

/// Solves an ODE system using Heun's method (second-order).

pub fn solve_ode_heun(
    funcs: &[Expr],
    y0: &[f64],
    x_range: (f64, f64),
    num_steps: usize,
) -> Result<Vec<Vec<f64>>, String> {

    let (x0, x_end) = x_range;

    let h = (x_end - x0)
        / (num_steps as f64);

    let mut x = x0;

    let mut y = y0.to_vec();

    let mut results = vec![y.clone()];

    let mut vars = HashMap::new();

    for _ in 0 .. num_steps {

        let k1 = eval_f(
            funcs,
            x,
            &y,
            &mut vars,
        )?;

        let y_pred = add_vec(
            &y,
            &scale_vec(&k1, h),
        );

        let k2 = eval_f(
            funcs,
            x + h,
            &y_pred,
            &mut vars,
        )?;

        let dy = scale_vec(
            &add_vec(&k1, &k2),
            0.5,
        );

        y = add_vec(
            &y,
            &scale_vec(&dy, h),
        );

        x += h;

        results.push(y.clone());
    }

    Ok(results)
}

/// Solves an ODE system using the fourth-order Runge-Kutta method.

pub fn solve_ode_system_rk4(
    funcs: &[Expr],
    y0: &[f64],
    x_range: (f64, f64),
    num_steps: usize,
) -> Result<Vec<Vec<f64>>, String> {

    let (x0, x_end) = x_range;

    let h = (x_end - x0)
        / (num_steps as f64);

    let mut x = x0;

    let mut y_vec = y0.to_vec();

    let mut results =
        vec![y_vec.clone()];

    let mut vars = HashMap::new();

    for _ in 0 .. num_steps {

        let k1 = eval_f(
            funcs,
            x,
            &y_vec,
            &mut vars,
        )?;

        let k2 = eval_f(
            funcs,
            x + h / 2.0,
            &add_vec(
                &y_vec,
                &scale_vec(
                    &k1,
                    h / 2.0,
                ),
            ),
            &mut vars,
        )?;

        let k3 = eval_f(
            funcs,
            x + h / 2.0,
            &add_vec(
                &y_vec,
                &scale_vec(
                    &k2,
                    h / 2.0,
                ),
            ),
            &mut vars,
        )?;

        let k4 = eval_f(
            funcs,
            x + h,
            &add_vec(
                &y_vec,
                &scale_vec(&k3, h),
            ),
            &mut vars,
        )?;

        let weighted_sum = add_vec(
            &add_vec(
                &k1,
                &scale_vec(&k2, 2.0),
            ),
            &add_vec(
                &scale_vec(&k3, 2.0),
                &k4,
            ),
        );

        y_vec = add_vec(
            &y_vec,
            &scale_vec(
                &weighted_sum,
                h / 6.0,
            ),
        );

        x += h;

        if y_vec
            .iter()
            .any(|&val| {

                !val.is_finite()
            })
        {

            return Err("Overflow or \
                        invalid value \
                        encountered \
                        during ODE \
                        solving."
                .to_string());
        }

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

    for (i, y_val) in y_vec
        .iter()
        .enumerate()
    {

        vars.insert(
            format!("y{i}"),
            *y_val,
        );
    }

    let mut results = Vec::new();

    for f in funcs {

        results
            .push(eval_expr(f, vars)?);
    }

    Ok(results)
}

pub(crate) fn add_vec(
    v1: &[f64],
    v2: &[f64],
) -> Vec<f64> {

    v1.iter()
        .zip(v2.iter())
        .map(|(a, b)| a + b)
        .collect()
}

pub(crate) fn scale_vec(
    v: &[f64],
    s: f64,
) -> Vec<f64> {

    v.iter()
        .map(|a| a * s)
        .collect()
}
