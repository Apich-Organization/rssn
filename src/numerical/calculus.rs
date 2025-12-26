//! # Numerical Calculus
//!
//! This module provides numerical calculus operations, primarily focusing on
//! finite difference methods for approximating derivatives. It includes functions
//! for computing the numerical gradient of multivariate functions.

use std::collections::HashMap;

use crate::numerical::elementary::eval_expr;
use crate::symbolic::core::Expr;

/// Computes the numerical partial derivative of a function with respect to a variable at a point.
///
/// # Example
/// ```rust
/// 
/// use rssn::numerical::calculus::partial_derivative;
/// use rssn::symbolic::core::Expr;
///
/// let x = Expr::new_variable("x");
///
/// let f = Expr::new_pow(
///     x,
///     Expr::new_constant(2.0),
/// ); // f(x) = x^2
/// let val = partial_derivative(&f, "x", 2.0).unwrap();
///
/// // Expect derivative at x=2 to be approx 4.0
/// assert!((val - 4.0).abs() < 1e-5);
/// ```

pub fn partial_derivative(
    f : &Expr,
    var : &str,
    x : f64,
) -> Result<f64, String> {

    let h = 1e-6;

    let mut vars_plus = HashMap::new();

    vars_plus.insert(
        var.to_string(),
        x + h,
    );

    let f_plus =
        eval_expr(f, &vars_plus)?;

    let mut vars_minus = HashMap::new();

    vars_minus.insert(
        var.to_string(),
        x - h,
    );

    let f_minus =
        eval_expr(f, &vars_minus)?;

    Ok((f_plus - f_minus) / (2.0 * h))
}

/// Computes the numerical gradient of a multivariate function at a given point.
///
/// The gradient is approximated using finite differences. For each variable, the partial
/// derivative is estimated by evaluating the function at `x + h` and `x - h`.
///
/// # Example
/// ```rust
/// 
/// use rssn::numerical::calculus::gradient;
/// use rssn::symbolic::core::Expr;
///
/// let x = Expr::new_variable("x");
///
/// let y = Expr::new_variable("y");
///
/// let f = Expr::new_add(
///     Expr::new_pow(
///         x,
///         Expr::new_constant(2.0),
///     ),
///     y,
/// ); // f(x,y) = x^2 + y
/// let grad = gradient(
///     &f,
///     &["x", "y"],
///     &[2.0, 3.0],
/// )
/// .unwrap();
///
/// // grad = [2x, 1] at (2,3) = [4, 1]
/// assert!((grad[0] - 4.0).abs() < 1e-5);
///
/// assert!((grad[1] - 1.0).abs() < 1e-5);
/// ```

pub fn gradient(
    f : &Expr,
    vars : &[&str],
    point : &[f64],
) -> Result<Vec<f64>, String> {

    if vars.len() != point.len() {

        return Err(
            "Number of variables must \
             match number of point \
             dimensions"
                .to_string(),
        );
    }

    let mut grad =
        Vec::with_capacity(vars.len());

    let h = 1e-6;

    for i in 0 .. vars.len() {

        let mut point_plus_h =
            point.to_vec();

        point_plus_h[i] += h;

        let mut point_minus_h =
            point.to_vec();

        point_minus_h[i] -= h;

        let f_plus_h = eval_at_point(
            f,
            vars,
            &point_plus_h,
        )?;

        let f_minus_h = eval_at_point(
            f,
            vars,
            &point_minus_h,
        )?;

        let partial_deriv = (f_plus_h
            - f_minus_h)
            / (2.0 * h);

        grad.push(partial_deriv);
    }

    Ok(grad)
}

/// Computes the numerical Jacobian matrix of a vector-valued function at a given point.
///
/// # Example
/// ```rust
/// 
/// use rssn::numerical::calculus::jacobian;
/// use rssn::symbolic::core::Expr;
///
/// let x = Expr::new_variable("x");
///
/// let y = Expr::new_variable("y");
///
/// // f1 = x*y, f2 = x + y
/// let f1 = Expr::new_mul(x.clone(), y.clone());
///
/// let f2 = Expr::new_add(x, y);
///
/// let jac = jacobian(
///     &[f1, f2],
///     &["x", "y"],
///     &[1.0, 2.0],
/// )
/// .unwrap();
///
/// // J = [[y, x], [1, 1]] at (1,2) = [[2, 1], [1, 1]]
/// assert!((jac[0][0] - 2.0).abs() < 1e-5);
/// ```

pub fn jacobian(
    funcs : &[Expr],
    vars : &[&str],
    point : &[f64],
) -> Result<Vec<Vec<f64>>, String> {

    let mut jac =
        Vec::with_capacity(funcs.len());

    for f in funcs {

        jac.push(gradient(
            f, vars, point,
        )?);
    }

    Ok(jac)
}

/// Computes the numerical Hessian matrix of a scalar function at a given point.
///
/// # Example
/// ```rust
/// 
/// use rssn::numerical::calculus::hessian;
/// use rssn::symbolic::core::Expr;
///
/// let x = Expr::new_variable("x");
///
/// let y = Expr::new_variable("y");
///
/// // f = x^2 + x*y + y^2
/// let f = Expr::new_add(
///     Expr::new_add(
///         Expr::new_pow(
///             x.clone(),
///             Expr::new_constant(2.0),
///         ),
///         Expr::new_mul(x, y.clone()),
///     ),
///     Expr::new_pow(
///         y,
///         Expr::new_constant(2.0),
///     ),
/// );
///
/// let hess = hessian(
///     &f,
///     &["x", "y"],
///     &[1.0, 1.0],
/// )
/// .unwrap();
///
/// // H = [[2, 1], [1, 2]]
/// assert!((hess[0][0] - 2.0).abs() < 1e-5);
///
/// assert!((hess[0][1] - 1.0).abs() < 1e-5);
/// ```

pub fn hessian(
    f : &Expr,
    vars : &[&str],
    point : &[f64],
) -> Result<Vec<Vec<f64>>, String> {

    if vars.len() != point.len() {

        return Err(
            "Number of variables must \
             match number of point \
             dimensions"
                .to_string(),
        );
    }

    let n = vars.len();

    let h = 1e-4; // Slightly larger h for second derivatives
    let mut hess =
        vec![vec![0.0; n]; n];

    for i in 0 .. n {

        for j in 0 ..= i {

            let val = if i == j {

                // fxx approx (f(x+h) - 2f(x) + f(x-h)) / h^2
                let mut p_plus =
                    point.to_vec();

                p_plus[i] += h;

                let mut p_minus =
                    point.to_vec();

                p_minus[i] -= h;

                let f_plus =
                    eval_at_point(
                        f,
                        vars,
                        &p_plus,
                    )?;

                let f_center =
                    eval_at_point(
                        f, vars, point,
                    )?;

                let f_minus =
                    eval_at_point(
                        f,
                        vars,
                        &p_minus,
                    )?;

                (2.0f64.mul_add(
                    -f_center,
                    f_plus,
                ) + f_minus)
                    / (h * h)
            } else {

                // fxy approx (f(x+h, y+h) - f(x+h, y-h) - f(x-h, y+h) + f(x-h, y-h)) / (4h^2)
                let mut p_pp =
                    point.to_vec();

                p_pp[i] += h;

                p_pp[j] += h;

                let mut p_pm =
                    point.to_vec();

                p_pm[i] += h;

                p_pm[j] -= h;

                let mut p_mp =
                    point.to_vec();

                p_mp[i] -= h;

                p_mp[j] += h;

                let mut p_mm =
                    point.to_vec();

                p_mm[i] -= h;

                p_mm[j] -= h;

                let f_pp =
                    eval_at_point(
                        f, vars, &p_pp,
                    )?;

                let f_pm =
                    eval_at_point(
                        f, vars, &p_pm,
                    )?;

                let f_mp =
                    eval_at_point(
                        f, vars, &p_mp,
                    )?;

                let f_mm =
                    eval_at_point(
                        f, vars, &p_mm,
                    )?;

                (f_pp - f_pm - f_mp
                    + f_mm)
                    / (4.0 * h * h)
            };

            hess[i][j] = val;

            hess[j][i] = val;
        }
    }

    Ok(hess)
}

/// Helper to evaluate a multivariate expression at a point.
///
/// This function substitutes the numerical values from `point` into the `vars`
/// of the `expr` and then numerically evaluates the resulting expression.
///
/// # Arguments
/// * `expr` - The expression to evaluate.
/// * `vars` - The variables of the expression.
/// * `point` - The numerical values for the variables.
///
/// # Returns
/// A `Result` containing the numerical value of the expression, or an error string.

pub(crate) fn eval_at_point(
    expr : &Expr,
    vars : &[&str],
    point : &[f64],
) -> Result<f64, String> {

    let mut vars_map = HashMap::new();

    for (i, &var) in vars
        .iter()
        .enumerate()
    {

        vars_map.insert(
            var.to_string(),
            point[i],
        );
    }

    eval_expr(expr, &vars_map)
}
