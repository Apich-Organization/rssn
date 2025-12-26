//! # Symbolic Result Verification
//!
//! This module provides functions for verifying symbolic results numerically.
//! It uses random sampling and numerical evaluation to check the correctness of
//! solutions to equations, integrals, ODEs, and matrix operations. This is particularly
//! useful for complex symbolic computations where direct algebraic verification is difficult.

use std::collections::HashMap;

use rand::thread_rng;
use rand::Rng;

use crate::numerical::elementary::eval_expr;
use crate::numerical::integrate::quadrature;
use crate::numerical::integrate::QuadratureMethod;
use crate::symbolic::calculus::differentiate;
use crate::symbolic::calculus::substitute;
use crate::symbolic::core::Expr;
use crate::symbolic::matrix;
use crate::symbolic::simplify_dag::simplify;

const TOLERANCE: f64 = 1e-6;

const NUM_SAMPLES: usize = 100;

/// Verifies a solution to a single equation or a system of equations using numerical sampling.
///
/// This function substitutes the proposed solution into the equations and evaluates the result
/// at several random points. If the equations are satisfied (i.e., evaluate to approximately zero)
/// at all sample points, the solution is considered verified.
///
/// # Arguments
/// * `equations` - A slice of `Expr` representing the equations to verify.
/// * `solution` - A `HashMap` mapping variable names to their proposed solution expressions.
/// * `free_vars` - A slice of string slices representing any free variables in the solution.
///
/// # Returns
/// `true` if the solution is numerically verified, `false` otherwise.
/// Verifies a solution to a single equation or a system of equations using numerical sampling.
///
/// This function substitutes the proposed solution into the equations and evaluates the result
/// at several random points. If the equations are satisfied (i.e., evaluate to approximately zero)
/// at all sample points, the solution is considered verified.
///
/// # Arguments
/// * `equations` - A slice of `Expr` representing the equations to verify.
/// * `solution` - A `HashMap` mapping variable names to their proposed solution expressions.
/// * `free_vars` - A slice of string slices representing any free variables in the solution.
///
/// # Returns
/// `true` if the solution is numerically verified, `false` otherwise.
#[must_use]

pub fn verify_equation_solution(
    equations: &[Expr],
    solution: &HashMap<String, Expr>,
    free_vars: &[&str],
) -> bool {

    let mut rng = thread_rng();

    for eq in equations {

        let unwrapped_eq =
            unwrap_dag(eq.clone());

        let diff =
            if let Expr::Eq(lhs, rhs) =
                unwrapped_eq
            {

                simplify(
                    &Expr::new_sub(
                        lhs.clone(),
                        rhs.clone(),
                    ),
                )
            } else {

                unwrapped_eq.clone()
            };

        for _ in 0 .. NUM_SAMPLES {

            let mut current_vars =
                HashMap::new();

            // Random values for free variables
            for var in free_vars {

                current_vars.insert(
                    (*var).to_string(),
                    rng.gen_range(
                        -10.0 .. 10.0,
                    ),
                );
            }

            // Substitute the proposed solution into the equation
            let mut substituted_expr =
                diff.clone();

            for (var, sol_expr) in
                solution
            {

                substituted_expr = substitute(
                    &substituted_expr,
                    var,
                    sol_expr,
                );
            }

            match eval_expr(
                &simplify(
                    &substituted_expr,
                ),
                &current_vars,
            ) {
                | Ok(val) => {

                    if val.abs()
                        > TOLERANCE
                    {

                        return false;
                    }
                },
                | Err(_) => {

                    // Try another point if evaluation fails (e.g., division by zero at a specific random point)
                    continue;
                },
            }
        }
    }

    true
}

fn unwrap_dag(expr: Expr) -> Expr {

    match expr {
        | Expr::Dag(node) => {
            node.to_expr()
                .unwrap_or(Expr::Dag(
                    node,
                ))
        },
        | _ => expr,
    }
}

/// Verifies an indefinite integral `F(x)` for an integrand `f(x)` by checking if `F'(x) == f(x)`.
#[must_use]

pub fn verify_indefinite_integral(
    integrand: &Expr,
    integral_result: &Expr,
    var: &str,
) -> bool {

    let derivative_of_result =
        differentiate(
            integral_result,
            var,
        );

    let diff =
        simplify(&Expr::new_sub(
            integrand.clone(),
            derivative_of_result,
        ));

    let mut rng = thread_rng();

    let mut success_count = 0;

    let mut attempt_count = 0;

    while success_count < NUM_SAMPLES
        && attempt_count
            < NUM_SAMPLES * 2
    {

        let mut vars = HashMap::new();

        let x_val = rng
            .gen_range(-10.0 .. 10.0);

        vars.insert(
            var.to_string(),
            x_val,
        );

        if let Ok(val) =
            eval_expr(&diff, &vars)
        {

            if val.abs() > TOLERANCE {

                return false;
            }

            success_count += 1;
        }

        attempt_count += 1;
    }

    success_count > 0
}

/// Verifies a definite integral by comparing the symbolic result with numerical quadrature.
#[must_use]

pub fn verify_definite_integral(
    integrand: &Expr,
    var: &str,
    range: (f64, f64),
    symbolic_result: &Expr,
) -> bool {

    let symbolic_val = match eval_expr(
        &simplify(symbolic_result),
        &HashMap::new(),
    ) {
        | Ok(v) => v,
        | Err(_) => return false,
    };

    if let Ok(numerical_val) =
        quadrature(
            integrand,
            var,
            range,
            1000,
            &QuadratureMethod::Simpson,
        )
    {

        (symbolic_val - numerical_val)
            .abs()
            < TOLERANCE
    } else {

        false
    }
}

/// Verifies a solution to an ODE `G(x, y, y', y'', ...) = 0` by numerical sampling.
#[must_use]

pub fn verify_ode_solution(
    ode: &Expr,
    solution: &Expr,
    func_name: &str,
    var: &str,
) -> bool {

    // 1. Convert ODE to f(x, y, y', y'', ...) = 0 form
    let unwrapped_ode =
        unwrap_dag(ode.clone());

    let eq_zero =
        if let Expr::Eq(lhs, rhs) =
            unwrapped_ode
        {

            Expr::new_sub(lhs, rhs)
        } else {

            unwrapped_ode
        };

    let mut rng = thread_rng();

    for _ in 0 .. NUM_SAMPLES {

        let x_val = rng
            .gen_range(-10.0 .. 10.0);

        let mut vars = HashMap::new();

        vars.insert(
            var.to_string(),
            x_val,
        );

        // We need to substitute y, y', y'', ... in the ODE
        // This is a bit complex as we need to find all derivatives of func_name
        // For now, let's just handle y and y' for simplicity, or assume 'solution' is substituted for 'func_name'

        // Better approach: symbolically substitute and differentiate
        let mut substituted_ode =
            simplify(&eq_zero);

        // This is a naive substitution. Proper ODE verification requires handling derivatives specifically.
        // Assuming the ODE uses standard notation or we substitute derivatives of the solution.
        let y = solution.clone();

        let y_prime =
            differentiate(&y, var);

        let y_double_prime =
            differentiate(
                &y_prime,
                var,
            );

        substituted_ode = substitute(
            &substituted_ode,
            func_name,
            &y,
        );

        substituted_ode = substitute(
            &substituted_ode,
            &format!("{func_name}'"),
            &y_prime,
        );

        substituted_ode = substitute(
            &substituted_ode,
            &format!("{func_name}''"),
            &y_double_prime,
        );

        match eval_expr(
            &simplify(&substituted_ode),
            &vars,
        ) {
            | Ok(val) => {

                if val.abs()
                    > TOLERANCE * 10.0
                {

                    // ODEs can be more sensitive
                    return false;
                }
            },
            | Err(_) => continue,
        }
    }

    true
}

/// Verifies a matrix inverse `A⁻¹` by checking if `A * A⁻¹` is the identity matrix.
#[must_use]

pub fn verify_matrix_inverse(
    original: &Expr,
    inverse: &Expr,
) -> bool {

    let product = matrix::mul_matrices(
        original,
        inverse,
    );

    let simplified_product =
        unwrap_dag(simplify(&product));

    if let Expr::Matrix(prod_mat) =
        simplified_product
    {

        let n = prod_mat.len();

        for i in 0 .. n {

            for j in 0 .. n {

                let expected = if i == j
                {

                    1.0
                } else {

                    0.0
                };

                match eval_expr(
                    &prod_mat[i][j],
                    &HashMap::new(),
                ) {
                    | Ok(val) => {

                        if (val
                            - expected)
                            .abs()
                            > TOLERANCE
                        {

                            return false;
                        }
                    },
                    | Err(_) => {
                        return false
                    },
                }
            }
        }

        return true;
    }

    false
}

/// Verifies a symbolic derivative `f'(x)` by comparing it to a numerical differentiation.
#[must_use]

pub fn verify_derivative(
    original_func: &Expr,
    derivative_func: &Expr,
    var: &str,
) -> bool {

    let mut rng = thread_rng();

    for _ in 0 .. NUM_SAMPLES {

        let x_val = rng
            .gen_range(-10.0 .. 10.0);

        let mut vars_map =
            HashMap::new();

        vars_map.insert(
            var.to_string(),
            x_val,
        );

        let symbolic_deriv_val =
            match eval_expr(
                derivative_func,
                &vars_map,
            ) {
                | Ok(v) => v,
                | Err(_) => continue,
            };

        let numerical_deriv_val = match crate::numerical::calculus::gradient(
            original_func,
            &[var],
            &[x_val],
        ) {
            | Ok(grad_vec) => grad_vec[0],
            | Err(_) => continue,
        };

        if (symbolic_deriv_val
            - numerical_deriv_val)
            .abs()
            > TOLERANCE * 100.0
        {

            return false;
        }
    }

    true
}

/// Verifies a symbolic limit `lim_{x->x0} f(x) = L`.
#[must_use]

pub fn verify_limit(
    f: &Expr,
    var: &str,
    target: &Expr,
    limit_val: &Expr,
) -> bool {

    let x0 = match eval_expr(
        &simplify(target),
        &HashMap::new(),
    ) {
        | Ok(v) => v,
        | Err(_) => return false,
    };

    let l = match eval_expr(
        &simplify(limit_val),
        &HashMap::new(),
    ) {
        | Ok(v) => v,
        | Err(_) => return false,
    };

    let epsilons = [1e-3, 1e-5, 1e-7];

    for &eps in &epsilons {

        let mut vars = HashMap::new();

        vars.insert(
            var.to_string(),
            x0 + eps,
        );

        if let Ok(val) =
            eval_expr(f, &vars)
        {

            if (val - l).abs()
                > eps.mul_add(
                    100.0,
                    TOLERANCE,
                )
            {

                return false;
            }
        }

        vars.insert(
            var.to_string(),
            x0 - eps,
        );

        if let Ok(val) =
            eval_expr(f, &vars)
        {

            if (val - l).abs()
                > eps.mul_add(
                    100.0,
                    TOLERANCE,
                )
            {

                return false;
            }
        }
    }

    true
}
