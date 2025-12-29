//! # Numerical Calculus of Variations
//!
//! This module provides numerical tools for the calculus of variations.
//! It focuses on evaluating the action of a functional for a given path,
//! which is a fundamental step in solving problems like finding geodesics
//! or optimizing paths in physics and engineering.

use crate::numerical::integrate::quadrature;
use crate::numerical::integrate::QuadratureMethod;
use crate::symbolic::calculus::differentiate;
use crate::symbolic::calculus::substitute;
use crate::symbolic::core::Expr;

/// Evaluates the action of a functional for a given path.
///
/// The functional is `S[y] = integral from a to b of L(t, y, y_dot) dt`.
/// This function numerically computes this integral by first substituting the `path`
/// and its derivative into the `lagrangian`, and then performing numerical quadrature.
///
/// # Arguments
/// * `lagrangian` - The Lagrangian expression `L`. It should be an expression in terms of
///   `t_var`, `path_var`, and `path_dot_var`.
/// * `path` - The actual path `y(t)` as an expression.
/// * `t_var` - The name of the independent variable (e.g., "t").
/// * `path_var` - The name of the path variable used in the Lagrangian (e.g., "y").
/// * `path_dot_var` - The name of the path's derivative used in the Lagrangian (e.g., "`y_dot`").
/// * `t_range` - The interval of integration `(a, b)`.
///
/// # Returns
/// A `Result` containing the numerical value of the action integral.
///
/// # Example
/// ```rust
/// 
/// use rssn::numerical::calculus_of_variations::evaluate_action;
/// use rssn::symbolic::core::Expr;
///
/// // L = 1/2 * (y_dot)^2 (Free particle action)
/// let t = Expr::new_variable("t");
///
/// let y = Expr::new_variable("y");
///
/// let y_dot = Expr::new_variable("y_dot");
///
/// let lagrangian = Expr::new_mul(
///     Expr::new_constant(0.5),
///     Expr::new_pow(
///         y_dot.clone(),
///         Expr::new_constant(2.0),
///     ),
/// );
///
/// // Path: y(t) = t
/// let path = t.clone();
///
/// let action = evaluate_action(
///     &lagrangian,
///     &path,
///     "t",
///     "y",
///     "y_dot",
///     (0.0, 1.0),
/// )
/// .unwrap();
///
/// // Integral of 0.5 from 0 to 1 is 0.5
/// assert!((action - 0.5).abs() < 1e-5);
/// ```
///
/// # Errors
///
/// Returns an error if the numerical quadrature fails.

pub fn evaluate_action(
    lagrangian: &Expr,
    path: &Expr,
    t_var: &str,
    path_var: &str,
    path_dot_var: &str,
    t_range: (f64, f64),
) -> Result<f64, String> {

    let path_dot =
        differentiate(path, t_var);

    let integrand_with_y = substitute(
        lagrangian,
        path_var,
        path,
    );

    let integrand = substitute(
        &integrand_with_y,
        path_dot_var,
        &path_dot,
    );

    quadrature(
        &integrand,
        t_var,
        t_range,
        1000,
        &QuadratureMethod::Simpson,
    )
}

/// Computes the Euler-Lagrange expression for a given Lagrangian.
///
/// The Euler-Lagrange equation is: `d/dt(dL/dy_dot) - dL/dy = 0`.
/// This function returns the symbolic expression `d/dt(dL/dy_dot) - dL/dy`.
///
/// # Arguments
/// * `lagrangian` - The Lagrangian expression `L(t, y, y_dot)`.
/// * `t_var` - The independent variable (e.g., "t").
/// * `path_var` - The path variable (e.g., "y").
/// * `path_dot_var` - The derivative of the path variable (e.g., "`y_dot`").
///
/// # Example
/// ```rust
/// 
/// use rssn::numerical::calculus_of_variations::euler_lagrange;
/// use rssn::symbolic::core::Expr;
///
/// // L = 1/2 * m * y_dot^2 - m * g * y
/// // EL: m * y_ddot + m * g = 0
/// ```

#[must_use]

pub fn euler_lagrange(
    lagrangian: &Expr,
    t_var: &str,
    path_var: &str,
    path_dot_var: &str,
) -> Expr {

    let dl_dy = differentiate(
        lagrangian,
        path_var,
    );

    let dl_dy_dot = differentiate(
        lagrangian,
        path_dot_var,
    );

    // We need to take d/dt (dl_dy_dot).
    // Since dl_dy_dot depends on t, y(t), and y_dot(t),
    // d/dt (dl_dy_dot) = d/dt (dl_dy_dot) + d/dy(dl_dy_dot)*y_dot + d/dy_dot(dl_dy_dot)*y_ddot
    // However, our symbolic differentiator only handles explicit t dependence.
    // For now, we'll return a simplified version or just the partials if full EL ODE generation is complex.
    // But we can actually perform the full chain rule if we define y_ddot as another variable.

    let y_dot_sym = Expr::new_variable(
        path_dot_var,
    );

    let y_ddot_sym = Expr::new_variable(
        &format!("{path_dot_var}_dot"),
    );

    let d_dt_explicit = differentiate(
        &dl_dy_dot,
        t_var,
    );

    let d_dy = differentiate(
        &dl_dy_dot,
        path_var,
    );

    let d_dy_dot = differentiate(
        &dl_dy_dot,
        path_dot_var,
    );

    // d/dt (dL/dy_dot) = dL_explicit_t/dy_dot + dL/dy/dy_dot * y_dot + dL/dy_dot/dy_dot * y_ddot
    let d_dt_total = Expr::new_add(
        Expr::new_add(
            d_dt_explicit,
            Expr::new_mul(
                d_dy,
                y_dot_sym,
            ),
        ),
        Expr::new_mul(
            d_dy_dot,
            y_ddot_sym,
        ),
    );

    Expr::new_sub(d_dt_total, dl_dy)
}
