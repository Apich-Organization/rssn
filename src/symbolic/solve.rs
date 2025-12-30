//! # Symbolic Equation Solving
//!
//! This module provides a powerful set of tools for solving equations and systems of equations.
//! It includes dispatchers that can handle polynomial, transcendental, linear, and multivariate
//! polynomial systems by selecting the appropriate algorithm, such as substitution, Gaussian
//! elimination, or Grobner bases.
//!
//! ## Examples
//!
//! ### Solving a Linear Equation
//! ```
//! 
//! use rssn::symbolic::core::Expr;
//! use rssn::symbolic::solve::solve;
//!
//! // Solve 2x + 4 = 0 for x
//! let eq = Expr::new_add(
//!     Expr::new_mul(
//!         Expr::new_constant(2.0),
//!         Expr::new_variable("x"),
//!     ),
//!     Expr::new_constant(4.0),
//! );
//!
//! let solutions = solve(&eq, "x");
//! // solutions contains [-2.0]
//! ```
//!
//! ### Solving a System of Linear Equations
//! ```
//! 
//! use rssn::symbolic::core::Expr;
//! use rssn::symbolic::solve::solve_linear_system;
//!
//! // System:
//! // x + y = 3
//! // x - y = 1
//! let eq1 = Expr::new_sub(
//!     Expr::new_add(
//!         Expr::new_variable("x"),
//!         Expr::new_variable("y"),
//!     ),
//!     Expr::new_constant(3.0),
//! );
//!
//! let eq2 = Expr::new_sub(
//!     Expr::new_sub(
//!         Expr::new_variable("x"),
//!         Expr::new_variable("y"),
//!     ),
//!     Expr::new_constant(1.0),
//! );
//!
//! let system = Expr::System(vec![eq1, eq2]);
//!
//! let solutions = solve_linear_system(
//!     &system,
//!     &[
//!         "x".to_string(),
//!         "y".to_string(),
//!     ],
//! )
//! .unwrap();
//! // solutions = [2.0, 1.0]
//! ```

use std::collections::HashMap;
use std::sync::Arc;

use num_traits::ToPrimitive;

use crate::symbolic::calculus::substitute;
use crate::symbolic::core::Expr;
use crate::symbolic::core::SparsePolynomial;
use crate::symbolic::grobner::buchberger;
use crate::symbolic::grobner::MonomialOrder;
use crate::symbolic::matrix::create_empty_matrix;
use crate::symbolic::matrix::get_matrix_dims;
use crate::symbolic::matrix::null_space;
use crate::symbolic::matrix::rref;
use crate::symbolic::polynomial::expr_to_sparse_poly;
use crate::symbolic::polynomial::sparse_poly_to_expr;
use crate::symbolic::simplify::collect_and_order_terms;
use crate::symbolic::simplify::is_zero;
use crate::symbolic::simplify_dag::simplify;

/// Solves a single equation for a given variable.
///
/// This function acts as a dispatcher, attempting to solve the equation by first
/// simplifying it and then trying different specialized solvers:
/// - **Polynomial Solver**: For algebraic equations up to quartic degree.
/// - **Transcendental Solver**: For equations involving trigonometric or exponential functions.
///
/// If the equation is not explicitly an `Expr::Eq`, it is treated as `expr = 0`.
///
/// # Arguments
/// * `expr` - The equation to solve (e.g., `Expr::Eq(lhs, rhs)` or `lhs - rhs`).
/// * `var` - The variable to solve for.
///
/// # Returns
/// A `Vec<Expr>` containing the symbolic solutions. If no explicit solution is found,
/// it may return an unevaluated `Expr::Solve` expression.
#[must_use]

pub fn solve(
    expr: &Expr,
    var: &str,
) -> Vec<Expr> {

    let equation =
        if let Expr::Eq(left, right) =
            expr
        {

            simplify(&Expr::new_sub(
                left.clone(),
                right.clone(),
            ))
        } else {

            expr.clone()
        };

    if let Some(solutions) =
        solve_polynomial(&equation, var)
    {

        return solutions;
    }

    if let Some(solutions) =
        solve_transcendental(
            &equation,
            var,
        )
    {

        return solutions;
    }

    vec![Expr::Solve(
        Arc::new(equation),
        var.to_string(),
    )]
}

/// Solves a system of multivariate equations.
///
/// This function acts as a dispatcher, attempting to solve the system using different strategies:
/// - **Substitution**: Iteratively solves for variables and substitutes them into other equations.
/// - **Grobner Bases**: For polynomial systems, computes a Grobner basis to simplify the system.
///
/// # Arguments
/// * `equations` - A slice of `Expr` representing the equations in the system.
/// * `vars` - A slice of string slices representing the variables to solve for.
///
/// # Returns
/// An `Option<Vec<(String, Expr)>>` containing a vector of `(variable_name, solution_expression)`
/// pairs if a solution is found, or `None` if the system cannot be solved by the implemented methods.
#[must_use]

pub fn solve_system(
    equations: &[Expr],
    vars: &[&str],
) -> Option<Vec<(Expr, Expr)>> {

    if let Some(solutions) =
        solve_system_by_substitution(
            equations,
            vars,
        )
    {

        return Some(solutions);
    }

    if let Some(solutions) =
        solve_system_with_grobner(
            equations,
            vars,
        )
    {

        return Some(solutions);
    }

    None
}

/// Solves a system of multivariate equations using iterative substitution and elimination.
///
/// This function attempts to solve for one variable at a time and substitute its solution
/// into the remaining equations. It is particularly effective for systems where variables
/// can be easily isolated.
///
/// # Arguments
/// * `equations` - A slice of `Expr` representing the equations in the system.
/// * `vars` - A slice of string slices representing the variables to solve for.
///
/// # Returns
/// An `Option<Vec<(Expr, Expr)>>` containing a vector of `(variable_expression, solution_expression)`
/// pairs if a partial or complete solution is found, or `None` if the system cannot be solved.
#[must_use]

pub fn solve_system_parcial(
    equations: &[Expr],
    vars: &[&str],
) -> Option<Vec<(Expr, Expr)>> {

    let mut remaining_eqs: Vec<Expr> =
        equations.to_vec();

    let mut solutions: HashMap<
        String,
        Expr,
    > = HashMap::new();

    let mut progress = true;

    while progress
        && !remaining_eqs.is_empty()
    {

        progress = false;

        let mut solved_eq_index : Option<usize> = None;

        for (i, eq) in remaining_eqs
            .iter()
            .enumerate()
        {

            let mut current_eq =
                eq.clone();

            for (
                solved_var,
                solution_expr,
            ) in &solutions
            {

                current_eq = substitute(
                    &current_eq,
                    solved_var,
                    solution_expr,
                );
            }

            let remaining_vars: Vec<
                &str,
            > = vars
                .iter()
                .filter(|v| {

                    !solutions
                        .contains_key(
                            **v,
                        )
                })
                .copied()
                .collect();

            if remaining_vars.len() == 1
            {

                let var_to_solve =
                    remaining_vars[0];

                let mut new_solutions =
                    solve(
                        &current_eq,
                        var_to_solve,
                    );

                if !new_solutions
                    .is_empty()
                {

                    let solution =
                        new_solutions
                            .remove(0);

                    solutions.insert(
                        var_to_solve
                            .to_string(
                            ),
                        solution,
                    );

                    solved_eq_index =
                        Some(i);

                    progress = true;

                    break;
                }
            }
        }

        if let Some(index) =
            solved_eq_index
        {

            remaining_eqs.remove(index);
        }
    }

    let mut final_solutions =
        HashMap::new();

    for var_name in vars
        .iter()
        .map(|s| (*s).to_string())
    {

        if let Some(mut solution) =
            solutions
                .get(&var_name)
                .cloned()
        {

            let mut changed = true;

            while changed {

                changed = false;

                for (
                    solved_var,
                    sol_expr,
                ) in &solutions
                {

                    if solved_var
                        != &var_name
                    {

                        let new_solution = substitute(
                            &solution,
                            solved_var,
                            sol_expr,
                        );

                        if new_solution
                            != solution
                        {

                            solution = new_solution;

                            changed =
                                true;
                        }
                    }
                }
            }

            final_solutions.insert(
                var_name,
                simplify(&solution),
            );
        }
    }

    if final_solutions.len()
        == vars.len()
    {

        Some(
            vars.iter()
                .map(|&v| {

                    (
                        Expr::Variable(v.to_string()),
                        match final_solutions.get(v) {
                            | Some(s) => s.clone(),
                            | _none => unreachable!(),
                        },
                    )
                })
                .collect(),
        )
    } else {

        None
    }
}

/// Solves a system of linear equations `Ax = b` for any `M x N` matrix `A`.
///
/// This function constructs an augmented matrix `[A | b]`, computes its Reduced Row Echelon Form (RREF),
/// and then analyzes the RREF to determine the nature of the solution:
/// - **Unique Solution**: Returns a column vector `x`.
/// - **Infinite Solutions**: Returns a parametric solution (particular solution + null space basis).
/// - **No Solution**: Returns `Expr::NoSolution`.
///
/// # Arguments
/// * `a` - An `Expr::Matrix` representing the coefficient matrix `A`.
/// * `b` - An `Expr::Matrix` representing the constant vector `b` (must be a column vector).
///
/// # Returns
/// A `Result` containing an `Expr` representing the solution (matrix, system, or no solution).
///
/// # Errors
///
/// This function will return an error if:
/// - `A` or `b` are not valid matrices.
/// - The row dimensions of `A` and `b` are incompatible.
/// - `b` is not a column vector.
/// - The `rref` computation fails.
/// - The `null_space` computation fails.

pub fn solve_linear_system_mat(
    a: &Expr,
    b: &Expr,
) -> Result<Expr, String> {

    let (a_rows, a_cols) =
        get_matrix_dims(a).ok_or_else(
            || {

                "A is not a valid \
                 matrix"
                    .to_string()
            },
        )?;

    let (b_rows, b_cols) =
        get_matrix_dims(b).ok_or_else(
            || {

                "b is not a valid \
                 matrix"
                    .to_string()
            },
        )?;

    if a_rows != b_rows {

        return Err("Matrix A and \
                    vector b have \
                    incompatible \
                    row dimensions"
            .to_string());
    }

    if b_cols != 1 {

        return Err("b must be a \
                    column vector"
            .to_string());
    }

    let Expr::Matrix(a_mat) = a else {

        unreachable!()
    };

    let Expr::Matrix(b_mat) = b else {

        unreachable!()
    };

    let mut augmented_mat =
        a_mat.clone();

    for i in 0 .. a_rows {

        augmented_mat[i]
            .push(b_mat[i][0].clone());
    }

    let rref_expr = rref(
        &Expr::Matrix(augmented_mat),
    )?;

    let Expr::Matrix(rref_mat) =
        rref_expr
    else {

        unreachable!()
    };

    for (i, _row) in rref_mat
        .iter()
        .take(a_rows)
        .enumerate()
    {

        let is_lhs_zero = rref_mat[i]
            [0 .. a_cols]
            .iter()
            .all(is_zero);

        if is_lhs_zero
            && !is_zero(
                &rref_mat[i][a_cols],
            )
        {

            return Ok(
                Expr::NoSolution,
            );
        }
    }

    let mut pivot_cols = Vec::new();

    let mut lead = 0;

    for row in rref_mat.iter().take(a_rows) {
        if lead >= a_cols {
            break;
        }

        let mut i = lead;

        while i < a_cols && is_zero(&row[i]) {
            i += 1;
        }

        if i < a_cols {
            pivot_cols.push(i);
            lead = i + 1;
        }
    }

    if (0..a_cols).all(|c| pivot_cols.contains(&c)) {

        let mut solution =
            create_empty_matrix(
                a_cols,
                1,
            );

        for (i, &p_col) in pivot_cols
            .iter()
            .enumerate()
        {

            solution[p_col][0] =
                rref_mat[i][a_cols]
                    .clone();
        }

        Ok(Expr::Matrix(
            solution,
        ))
    } else {

        let particular_solution = {

            let mut sol =
                create_empty_matrix(
                    a_cols,
                    1,
                );

            for (i, &p_col) in
                pivot_cols
                    .iter()
                    .enumerate()
            {

                sol[p_col][0] =
                    rref_mat[i][a_cols]
                        .clone();
            }

            sol
        };

        let null_space_basis =
            null_space(a)?;

        Ok(Expr::System(vec![
            Expr::Matrix(
                particular_solution,
            ),
            null_space_basis,
        ]))
    }
}

/// Solves a system of linear equations symbolically using Gaussian elimination.
///
/// This function takes a system of equations and a list of variables, and attempts
/// to find symbolic solutions for each variable. It leverages the `solve_system`
/// dispatcher internally.
///
/// # Arguments
/// * `system` - An `Expr::System` containing `Expr::Eq` expressions.
/// * `vars` - A slice of strings representing the variables to solve for.
///
/// # Returns
/// A `Result` containing a vector of `Expr` representing the solutions for `vars`.
///
/// # Errors
///
/// This function will return an error if:
/// - The input `system` is not an `Expr::System`.
/// - The underlying `solve_system` dispatcher fails to find a solution.

pub fn solve_linear_system(
    system: &Expr,
    vars: &[String],
) -> Result<Vec<Expr>, String> {

    if let Expr::System(eqs) = system {

        let vars_str : Vec<&str> = vars
            .iter()
            .map(std::string::String::as_str)
            .collect();

        match solve_system(
            eqs,
            &vars_str,
        ) {
            | Some(solutions) => {

                let mut sol_map : HashMap<Expr, Expr> = solutions
                    .into_iter()
                    .collect();

                let ordered_solutions : Vec<Expr> = vars
                    .iter()
                    .map(|var| {

                        sol_map
                            .remove(&Expr::Variable(
                                var.clone(),
                            ))
                            .unwrap_or(Expr::Variable(
                                "NotFound".to_string(),
                            ))
                    })
                    .collect();

                Ok(ordered_solutions)
            },
            | _none => {
                Err("System could not \
                     be solved."
                    .to_string())
            },
        }
    } else {

        Err(
            "Input must be a system \
             of equations."
                .to_string(),
        )
    }
}

/// Solves a system of linear equations symbolically using Gaussian elimination.
///
/// This function constructs an augmented matrix from the system of equations and
/// performs Gaussian elimination to transform it into row echelon form, from which
/// the solutions can be directly read.
///
/// # Arguments
/// * `system` - An `Expr::System` containing `Expr::Eq` expressions.
/// * `vars` - A slice of strings representing the variables to solve for.
///
/// # Returns
/// A `Result` containing a vector of `Expr` representing the solutions for `vars`.
///
/// # Errors
///
/// This function will return an error if:
/// - The input is not an `Expr::System`.
/// - The number of equations does not match the number of variables.
/// - An element in the system is not a valid `Expr::Eq`.
/// - The system matrix is singular or underdetermined, preventing a unique solution.

pub fn solve_linear_system_gauss(
    system: &Expr,
    vars: &[String],
) -> Result<Vec<Expr>, String> {

    if let Expr::System(eqs) = system {

        let n = vars.len();

        if eqs.len() != n {

            return Err(format!(
                "Number of equations \
                 ({}) does not match \
                 number of variables \
                 ({})",
                eqs.len(),
                n
            ));
        }

        let mut matrix_a = vec![
            vec![
                    Expr::Constant(0.0);
                    n
                ];
            n
        ];

        let mut vector_b = vec![
                Expr::Constant(0.0);
                n
            ];

        for (i, eq) in eqs
            .iter()
            .enumerate()
        {

            let (lhs, rhs) = match eq {
                | Expr::Eq(l, r) => {
                    (l, r)
                },
                | _ => {
                    return Err(
                        format!(
                        "Item {i} is \
                         not a valid \
                         equation"
                    ),
                    )
                },
            };

            vector_b[i] =
                rhs.as_ref().clone();

            if let Some(coeffs) = extract_polynomial_coeffs(lhs, "") {

                for (_term_str, _coeff) in coeffs
                    .iter()
                    .zip(vars.iter())
                {}
            }

            let (_, terms) =
                collect_and_order_terms(
                    lhs,
                );

            for (term, coeff) in terms {

                if let Some(j) = vars
                    .iter()
                    .position(|v| {

                        v == &term
                            .to_string()
                    })
                {

                    matrix_a[i][j] =
                        coeff;
                } else if !is_zero(
                    &coeff,
                ) && term
                    .to_string()
                    != "1"
                {

                    vector_b[i] = simplify(&Expr::new_sub(
                        vector_b[i].clone(),
                        Expr::new_mul(coeff, term),
                    ));
                } else if term
                    .to_string()
                    == "1"
                {

                    vector_b[i] = simplify(&Expr::new_sub(
                        vector_b[i].clone(),
                        coeff,
                    ));
                }
            }
        }

        for i in 0 .. n {

            let mut max_row = i;

            for (k, _item) in matrix_a
                .iter()
                .enumerate()
                .take(n)
                .skip(i + 1)
            {

                if !is_zero(
                    &matrix_a[k][i],
                ) {

                    max_row = k;

                    break;
                }
            }

            matrix_a.swap(i, max_row);

            vector_b.swap(i, max_row);

            let pivot =
                matrix_a[i][i].clone();

            if is_zero(&pivot) {

                return Err("Matrix is singular or underdetermined".to_string());
            }

            for item in matrix_a[i].iter_mut().take(n).skip(i) {
                *item = simplify(
                    &Expr::new_div(
                        item.clone(),
                        pivot.clone(),
                    ),
                );
            }

            vector_b[i] = simplify(
                &Expr::new_div(
                    vector_b[i].clone(),
                    pivot.clone(),
                ),
            );

            for k in 0 .. n {

                if i != k {

                    let factor =
                        matrix_a[k][i]
                            .clone();

                    let (row_i, row_k) = if i < k {
                        let (start, end) = matrix_a.split_at_mut(k);
                        (&start[i], &mut end[0])
                    } else {
                        let (start, end) = matrix_a.split_at_mut(i);
                        (&end[0], &mut start[k])
                    };

                    for (item_k, item_i) in row_k.iter_mut().zip(row_i.iter()).take(n).skip(i) {
                        let term = simplify(&Expr::new_mul(
                            factor.clone(),
                            item_i.clone(),
                        ));

                        *item_k = simplify(&Expr::new_sub(
                            item_k.clone(),
                            term,
                        ));
                    }

                    let term_b = simplify(&Expr::new_mul(
                        factor.clone(),
                        vector_b[i].clone(),
                    ));

                    vector_b[k] = simplify(&Expr::new_sub(
                        vector_b[k].clone(),
                        term_b,
                    ));
                }
            }
        }

        Ok(vector_b)
    } else {

        Err(
            "Input expression is not \
             a system of equations"
                .to_string(),
        )
    }
}

pub(crate) fn solve_system_by_substitution(
    equations: &[Expr],
    vars: &[&str],
) -> Option<Vec<(Expr, Expr)>> {

    let mut remaining_eqs: Vec<Expr> =
        equations.to_vec();

    let mut solutions: HashMap<
        Expr,
        Expr,
    > = HashMap::new();

    let mut progress = true;

    while progress
        && !remaining_eqs.is_empty()
    {

        progress = false;

        let mut solved_eq_index : Option<usize> = None;

        for (i, eq) in remaining_eqs
            .iter()
            .enumerate()
        {

            let mut current_eq =
                eq.clone();

            for (
                solved_var,
                solution_expr,
            ) in &solutions
            {

                current_eq = substitute(
                    &current_eq,
                    &solved_var
                        .to_string(),
                    solution_expr,
                );
            }

            let remaining_vars : Vec<&str> = vars
                .iter()
                .filter(|v| {

                    !solutions.contains_key(&Expr::Variable(
                        (**v).to_string(),
                    ))
                })
                .copied()
                .collect();

            if remaining_vars.len() == 1
            {

                let var_to_solve =
                    remaining_vars[0];

                let mut new_solutions =
                    solve(
                        &current_eq,
                        var_to_solve,
                    );

                if !new_solutions
                    .is_empty()
                {

                    let solution =
                        new_solutions
                            .remove(0);

                    solutions.insert(
                        Expr::Variable(var_to_solve.to_string()),
                        solution,
                    );

                    solved_eq_index =
                        Some(i);

                    progress = true;

                    break;
                }
            }
        }

        if let Some(index) =
            solved_eq_index
        {

            remaining_eqs.remove(index);
        }
    }

    if solutions.len() != vars.len() {

        return None;
    }

    let mut final_solutions =
        HashMap::new();

    for &var_name_str in vars {

        let var_expr = Expr::Variable(
            var_name_str.to_string(),
        );

        if let Some(mut solution) =
            solutions
                .get(&var_expr)
                .cloned()
        {

            for (
                solved_var,
                sol_expr,
            ) in &solutions
            {

                if solved_var
                    != &var_expr
                {

                    solution = substitute(
                        &solution,
                        &solved_var.to_string(),
                        sol_expr,
                    );
                }
            }

            final_solutions.insert(
                var_expr,
                simplify(&solution),
            );
        }
    }

    Some(
        final_solutions
            .into_iter()
            .collect(),
    )
}

pub(crate) fn solve_system_with_grobner(
    equations: &[Expr],
    vars: &[&str],
) -> Option<Vec<(Expr, Expr)>> {

    let basis: Vec<SparsePolynomial> =
        equations
            .iter()
            .map(|eq| {

                expr_to_sparse_poly(
                    eq, vars,
                )
            })
            .collect();

    let grobner_basis = match buchberger(
        &basis,
        MonomialOrder::Lexicographical,
    ) {
        | Ok(basis) => basis,
        | Err(_) => return None,
    };

    let mut solutions: HashMap<
        Expr,
        Expr,
    > = HashMap::new();

    for poly in grobner_basis
        .iter()
        .rev()
    {

        let mut current_eq =
            sparse_poly_to_expr(poly);

        for (var, val) in &solutions {

            current_eq = substitute(
                &current_eq,
                &var.to_string(),
                val,
            );
        }

        let remaining_vars: Vec<&str> =
            vars.iter()
                .filter(|v| {

                    contains_var(
                        &current_eq,
                        v,
                    )
                })
                .copied()
                .collect();

        if remaining_vars.len() == 1 {

            let roots = solve(
                &current_eq,
                remaining_vars[0],
            );

            if roots.is_empty() {

                return None;
            }

            solutions.insert(
                Expr::Variable(
                    remaining_vars[0]
                        .to_string(),
                ),
                roots[0].clone(),
            );
        } else if !remaining_vars
            .is_empty()
            && !is_zero(&current_eq)
        {

            return None;
        }
    }

    if solutions.len() == vars.len() {

        Some(
            solutions
                .into_iter()
                .collect(),
        )
    } else {

        None
    }
}

pub(crate) fn solve_polynomial(
    expr: &Expr,
    var: &str,
) -> Option<Vec<Expr>> {

    // Handle Expr::Eq by converting to lhs - rhs
    let normalized_expr =
        if let Expr::Eq(left, right) =
            expr
        {

            Expr::new_sub(
                left.clone(),
                right.clone(),
            )
        } else {

            expr.clone()
        };

    let poly = expr_to_sparse_poly(
        &normalized_expr,
        &[var],
    );

    let expanded_expr =
        sparse_poly_to_expr(&poly);

    // eprintln!("solve_polynomial: expr={:?}, var={}, expanded={:?}", expr, var, expanded_expr);
    let coeffs =
        extract_polynomial_coeffs(
            &expanded_expr,
            var,
        )?;

    // eprintln!("solve_polynomial: coeffs={:?}", coeffs);
    let degree = coeffs.len() - 1;

    match degree {
        | 0 => Some(vec![]),
        | 1 => {
            Some(solve_linear(
                &coeffs,
            ))
        },
        | 2 => {
            Some(solve_quadratic(
                &coeffs,
            ))
        },
        | 3 => {
            Some(solve_cubic(&coeffs))
        },
        | 4 => {
            Some(solve_quartic(
                &coeffs,
            ))
        },
        | _ => {

            let poly_expr =
                expr.clone();

            let mut roots = Vec::new();

            for i in 0 .. degree {

                roots.push(
                    Expr::RootOf {
                        poly: Arc::new(
                            poly_expr
                                .clone(
                                ),
                        ),
                        index: i as u32,
                    },
                );
            }

            Some(roots)
        },
    }
}

pub(crate) fn solve_linear(
    coeffs: &[Expr]
) -> Vec<Expr> {

    let a = &coeffs[0];

    let b = &coeffs[1];

    vec![simplify(
        &Expr::Neg(Arc::new(
            Expr::Div(
                Arc::new(b.clone()),
                Arc::new(a.clone()),
            ),
        )),
    )]
}

pub(crate) fn solve_quadratic(
    coeffs: &[Expr]
) -> Vec<Expr> {

    let a = &coeffs[0];

    let b = &coeffs[1];

    let c = &coeffs[2];

    let discriminant =
        simplify(&Expr::new_sub(
            Expr::new_pow(
                b.clone(),
                Expr::Constant(2.0),
            ),
            Expr::new_mul(
                Expr::Constant(4.0),
                Expr::new_mul(
                    a.clone(),
                    c.clone(),
                ),
            ),
        ));

    let sqrt_d = simplify(
        &Expr::new_sqrt(discriminant),
    );

    let two_a =
        simplify(&Expr::new_mul(
            Expr::Constant(2.0),
            a.clone(),
        ));

    vec![
        simplify(&Expr::Div(
            Arc::new(Expr::Add(
                Arc::new(Expr::Neg(
                    Arc::new(b.clone()),
                )),
                Arc::new(
                    sqrt_d.clone(),
                ),
            )),
            Arc::new(two_a.clone()),
        )),
        simplify(&Expr::Div(
            Arc::new(Expr::Sub(
                Arc::new(Expr::Neg(
                    Arc::new(b.clone()),
                )),
                Arc::new(sqrt_d),
            )),
            Arc::new(two_a),
        )),
    ]
}

pub(crate) fn solve_cubic(
    coeffs: &[Expr]
) -> Vec<Expr> {

    let a = &coeffs[0];

    let b = &simplify(&Expr::new_div(
        coeffs[1].clone(),
        a.clone(),
    ));

    let c = &simplify(&Expr::new_div(
        coeffs[2].clone(),
        a.clone(),
    ));

    let d = &simplify(&Expr::new_div(
        coeffs[3].clone(),
        a.clone(),
    ));

    // Reduced cubic: y^3 + py + q = 0, where x = y - b/3
    let p = simplify(&Expr::new_sub(
        c.clone(),
        Expr::new_div(
            Expr::new_pow(
                b.clone(),
                Expr::Constant(2.0),
            ),
            Expr::Constant(3.0),
        ),
    ));

    // q = d - bc/3 + 2b^3/27
    let q = simplify(&Expr::new_add(
        d.clone(),
        Expr::new_add(
            Expr::new_neg(
                Expr::new_div(
                    Expr::new_mul(
                        b.clone(),
                        c.clone(),
                    ),
                    Expr::Constant(3.0),
                ),
            ),
            Expr::new_mul(
                Expr::Constant(
                    2.0 / 27.0,
                ),
                Expr::new_pow(
                    b.clone(),
                    Expr::Constant(3.0),
                ),
            ),
        ),
    ));

    let inner_sqrt =
        simplify(&Expr::new_add(
            Expr::new_pow(
                Expr::new_div(
                    q.clone(),
                    Expr::Constant(2.0),
                ),
                Expr::Constant(2.0),
            ),
            Expr::new_pow(
                Expr::new_div(
                    p,
                    Expr::Constant(3.0),
                ),
                Expr::Constant(3.0),
            ),
        ));

    let sqrt_inner =
        Expr::new_sqrt(inner_sqrt);

    let u_term =
        simplify(&Expr::new_add(
            Expr::new_neg(
                Expr::new_div(
                    q.clone(),
                    Expr::Constant(2.0),
                ),
            ),
            sqrt_inner.clone(),
        ));

    let v_term =
        simplify(&Expr::new_sub(
            Expr::new_neg(
                Expr::new_div(
                    q,
                    Expr::Constant(2.0),
                ),
            ),
            sqrt_inner,
        ));

    let u = simplify(&Expr::new_pow(
        u_term,
        Expr::Constant(1.0 / 3.0),
    ));

    let v = simplify(&Expr::new_pow(
        v_term,
        Expr::Constant(1.0 / 3.0),
    ));

    let omega = Expr::new_complex(
        Expr::Constant(-0.5),
        Expr::new_div(
            Expr::new_sqrt(
                Expr::Constant(3.0),
            ),
            Expr::Constant(2.0),
        ),
    );

    let omega2 = Expr::new_complex(
        Expr::Constant(-0.5),
        Expr::new_neg(Expr::new_div(
            Expr::new_sqrt(
                Expr::Constant(3.0),
            ),
            Expr::Constant(2.0),
        )),
    );

    let sub_term =
        simplify(&Expr::new_div(
            b.clone(),
            Expr::Constant(3.0),
        ));

    let root1 =
        simplify(&Expr::new_sub(
            Expr::new_add(
                u.clone(),
                v.clone(),
            ),
            sub_term.clone(),
        ));

    let root2 =
        simplify(&Expr::new_sub(
            Expr::new_add(
                Expr::new_mul(
                    omega.clone(),
                    u.clone(),
                ),
                Expr::new_mul(
                    omega2.clone(),
                    v.clone(),
                ),
            ),
            sub_term.clone(),
        ));

    let root3 =
        simplify(&Expr::new_sub(
            Expr::new_add(
                Expr::new_mul(
                    omega2,
                    u,
                ),
                Expr::new_mul(omega, v),
            ),
            sub_term,
        ));

    vec![root1, root2, root3]
}

pub(crate) fn solve_quartic(
    coeffs: &[Expr]
) -> Vec<Expr> {

    if coeffs.len() < 5 {

        return vec![];
    }

    let a = &coeffs[0];

    let b = &simplify(&Expr::new_div(
        coeffs[1].clone(),
        a.clone(),
    ));

    let c = &simplify(&Expr::new_div(
        coeffs[2].clone(),
        a.clone(),
    ));

    let d = &simplify(&Expr::new_div(
        coeffs[3].clone(),
        a.clone(),
    ));

    let e = &simplify(&Expr::new_div(
        coeffs[4].clone(),
        a.clone(),
    ));

    // Reduced quartic: y^4 + py^2 + qy + r = 0, where x = y - b/4
    let b2 = Expr::new_pow(
        b.clone(),
        Expr::Constant(2.0),
    );

    let b3 = Expr::new_pow(
        b.clone(),
        Expr::Constant(3.0),
    );

    let b4 = Expr::new_pow(
        b.clone(),
        Expr::Constant(4.0),
    );

    let p = simplify(&Expr::new_sub(
        c.clone(),
        Expr::new_mul(
            Expr::Constant(3.0 / 8.0),
            b2.clone(),
        ),
    ));

    let q = simplify(&Expr::new_add(
        d.clone(),
        Expr::new_add(
            Expr::new_neg(
                Expr::new_mul(
                    Expr::Constant(0.5),
                    Expr::new_mul(
                        b.clone(),
                        c.clone(),
                    ),
                ),
            ),
            Expr::new_mul(
                Expr::Constant(0.125),
                b3,
            ),
        ),
    ));

    let r = simplify(&Expr::new_add(
        e.clone(),
        Expr::new_add(
            Expr::new_neg(
                Expr::new_mul(
                    Expr::Constant(
                        0.25,
                    ),
                    Expr::new_mul(
                        b.clone(),
                        d.clone(),
                    ),
                ),
            ),
            Expr::new_add(
                Expr::new_mul(
                    Expr::Constant(
                        1.0 / 16.0,
                    ),
                    Expr::new_mul(
                        b2,
                        c.clone(),
                    ),
                ),
                Expr::new_neg(
                    Expr::new_mul(
                        Expr::Constant(
                            3.0 / 256.0,
                        ),
                        b4,
                    ),
                ),
            ),
        ),
    ));

    if is_zero(&q) {

        // Biquadratic case: y^4 + py^2 + r = 0
        let discriminant =
            simplify(&Expr::new_sub(
                Expr::new_pow(
                    p.clone(),
                    Expr::Constant(2.0),
                ),
                Expr::new_mul(
                    Expr::Constant(4.0),
                    r,
                ),
            ));

        let sqrt_disc = Expr::new_sqrt(
            discriminant,
        );

        let two = Expr::Constant(2.0);

        let y2_1 =
            simplify(&Expr::new_div(
                Expr::new_add(
                    Expr::new_neg(
                        p.clone(),
                    ),
                    sqrt_disc.clone(),
                ),
                two.clone(),
            ));

        let y2_2 =
            simplify(&Expr::new_div(
                Expr::new_sub(
                    Expr::new_neg(p),
                    sqrt_disc,
                ),
                two,
            ));

        let y1 = Expr::new_sqrt(
            y2_1.clone(),
        );

        let y2 = Expr::new_neg(
            Expr::new_sqrt(y2_1),
        );

        let y3 = Expr::new_sqrt(
            y2_2.clone(),
        );

        let y4 = Expr::new_neg(
            Expr::new_sqrt(y2_2),
        );

        let b_over_4 =
            simplify(&Expr::new_div(
                b.clone(),
                Expr::Constant(4.0),
            ));

        return vec![
            simplify(&Expr::new_sub(
                y1,
                b_over_4.clone(),
            )),
            simplify(&Expr::new_sub(
                y2,
                b_over_4.clone(),
            )),
            simplify(&Expr::new_sub(
                y3,
                b_over_4.clone(),
            )),
            simplify(&Expr::new_sub(
                y4,
                b_over_4,
            )),
        ];
    }

    // Resolvent cubic: 8m^3 + 8pm^2 + (2p^2 - 8r)m - q^2 = 0
    let cubic_coeffs = [
        Expr::Constant(8.0),
        Expr::new_mul(
            Expr::Constant(8.0),
            p.clone(),
        ),
        Expr::new_sub(
            Expr::new_mul(
                Expr::Constant(2.0),
                Expr::new_pow(
                    p.clone(),
                    Expr::Constant(2.0),
                ),
            ),
            Expr::new_mul(
                Expr::Constant(8.0),
                r,
            ),
        ),
        Expr::new_neg(Expr::new_pow(
            q.clone(),
            Expr::Constant(2.0),
        )),
    ];

    let m_roots =
        solve_cubic(&cubic_coeffs);

    let m = m_roots[0].clone(); // Just pick one root

    let sqrt_2m =
        Expr::new_sqrt(Expr::new_mul(
            Expr::Constant(2.0),
            m.clone(),
        ));

    let q_over_2sqrt_2m = Expr::new_div(
        q,
        Expr::new_mul(
            Expr::Constant(2.0),
            sqrt_2m.clone(),
        ),
    );

    // Quadratic 1: y^2 - sqrt(2m)y + (p/2 + m + q/(2sqrt(2m))) = 0
    let quad1_coeffs = [
        Expr::Constant(1.0),
        Expr::new_neg(sqrt_2m.clone()),
        Expr::new_add(
            Expr::new_mul(
                Expr::Constant(0.5),
                p.clone(),
            ),
            Expr::new_add(
                m.clone(),
                q_over_2sqrt_2m.clone(),
            ),
        ),
    ];

    // Quadratic 2: y^2 + sqrt(2m)y + (p/2 + m - q/(2sqrt(2m))) = 0
    let quad2_coeffs = [
        Expr::Constant(1.0),
        sqrt_2m,
        Expr::new_sub(
            Expr::new_add(
                Expr::new_mul(
                    Expr::Constant(0.5),
                    p,
                ),
                m,
            ),
            q_over_2sqrt_2m,
        ),
    ];

    let y_roots1 =
        solve_quadratic(&quad1_coeffs);

    let y_roots2 =
        solve_quadratic(&quad2_coeffs);

    let b_over_4 =
        simplify(&Expr::new_div(
            b.clone(),
            Expr::Constant(4.0),
        ));

    let mut solutions = Vec::new();

    for y in y_roots1
        .into_iter()
        .chain(y_roots2.into_iter())
    {

        solutions.push(simplify(
            &Expr::new_sub(
                y,
                b_over_4.clone(),
            ),
        ));
    }

    solutions
}

pub(crate) fn solve_transcendental(
    expr: &Expr,
    var: &str,
) -> Option<Vec<Expr>> {

    if let Expr::Sub(lhs, rhs) = expr {

        return solve_transcendental_pattern(lhs, rhs, var);
    }

    if let Expr::Add(lhs, rhs) = expr {

        return solve_transcendental_pattern(
            lhs,
            &Expr::new_neg(rhs.clone()),
            var,
        );
    }

    None
}

pub(crate) fn solve_transcendental_pattern(
    lhs: &Expr,
    rhs: &Expr,
    var: &str,
) -> Option<Vec<Expr>> {

    let n =
        Expr::Variable("k".to_string());

    let pi = Expr::Pi;

    let (func_part, const_part) =
        if contains_var(lhs, var)
            && !contains_var(rhs, var)
        {

            (lhs, rhs)
        } else if !contains_var(
            lhs, var,
        ) && contains_var(
            rhs, var,
        ) {

            (rhs, lhs)
        } else {

            return None;
        };

    match func_part {
        | Expr::Sin(arg) => {

            let inner_solutions = solve(
                &Expr::Eq(
                    arg.clone(),
                    Arc::new(Expr::new_add(
                        Expr::new_mul(n.clone(), pi),
                        Expr::new_mul(
                            Expr::new_pow(
                                Expr::Constant(-1.0),
                                n,
                            ),
                            Expr::new_arcsin(const_part.clone()),
                        ),
                    )),
                ),
                var,
            );

            Some(inner_solutions)
        },
        | Expr::Cos(arg) => {

            let sol1 = solve(
                &Expr::Eq(
                    arg.clone(),
                    Arc::new(Expr::new_add(
                        Expr::new_mul(
                            Expr::Constant(2.0),
                            Expr::new_mul(
                                n.clone(),
                                pi.clone(),
                            ),
                        ),
                        Expr::new_arccos(const_part.clone()),
                    )),
                ),
                var,
            );

            let sol2 = solve(
                &Expr::Eq(
                    arg.clone(),
                    Arc::new(Expr::new_sub(
                        Expr::new_mul(
                            Expr::Constant(2.0),
                            Expr::new_mul(n, pi),
                        ),
                        Expr::new_arccos(const_part.clone()),
                    )),
                ),
                var,
            );

            Some([sol1, sol2].concat())
        },
        | Expr::Tan(arg) => {

            let inner_solutions = solve(
                &Expr::Eq(
                    arg.clone(),
                    Arc::new(Expr::new_add(
                        Expr::new_mul(n, pi),
                        Expr::new_arctan(const_part.clone()),
                    )),
                ),
                var,
            );

            Some(inner_solutions)
        },
        | Expr::Exp(arg) => {

            let i = Expr::new_complex(
                Expr::Constant(0.0),
                Expr::Constant(1.0),
            );

            let log_sol = Expr::new_add(
                Expr::new_log(
                    const_part.clone(),
                ),
                Expr::new_mul(
                    Expr::new_mul(
                        Expr::Constant(
                            2.0,
                        ),
                        Expr::new_mul(
                            pi, i,
                        ),
                    ),
                    n,
                ),
            );

            Some(solve(
                &Expr::Eq(
                    arg.clone(),
                    Arc::new(log_sol),
                ),
                var,
            ))
        },
        | _ => None,
    }
}

pub(crate) fn contains_var(
    expr: &Expr,
    var: &str,
) -> bool {

    let mut found = false;

    expr.pre_order_walk(&mut |e| {
        if let Expr::Variable(v) = e {

            if v == var {

                found = true;
            }
        }
    });

    found
}

#[must_use]

/// Extracts polynomial coefficients from an expression with respect to a variable.
///
/// This function attempts to represent the given expression as a polynomial in `var`
/// and returns the coefficients in descending order of degree (e.g., [`a_n`, ..., `a_1`, `a_0`]).
///
/// # Arguments
/// * `expr` - The symbolic expression to analyze.
/// * `var` - The variable with respect to which coefficients are extracted.
///
/// # Returns
/// An `Option<Vec<Expr>>` containing the coefficients if the expression is a polynomial,
/// or `None` if it contains transcendental terms or other non-polynomial structures in `var`.

pub fn extract_polynomial_coeffs(
    expr: &Expr,
    var: &str,
) -> Option<Vec<Expr>> {

    let mut coeffs_map = HashMap::new();

    collect_coeffs(
        expr,
        var,
        &mut coeffs_map,
        &Expr::Constant(1.0),
    )?;

    if coeffs_map.is_empty() {

        if contains_var(expr, var) {

            return None;
        } else {

            let mut map =
                HashMap::new();

            map.insert(0, expr.clone());

            coeffs_map = map;
        }
    }

    let max_degree = *coeffs_map
        .keys()
        .max()
        .unwrap_or(&0);

    let mut coeffs =
        vec![
            Expr::Constant(0.0);
            max_degree as usize + 1
        ];

    for (degree, coeff) in coeffs_map {

        coeffs[degree as usize] =
            simplify(&coeff);
    }

    coeffs.reverse();

    Some(coeffs)
}

pub(crate) fn collect_coeffs(
    expr: &Expr,
    var: &str,
    coeffs: &mut HashMap<u32, Expr>,
    factor: &Expr,
) -> Option<()> {

    match expr {
        | Expr::Dag(node) => {
            collect_coeffs(
                &node
                    .to_expr()
                    .expect(
                        "Dag Coeffs",
                    ),
                var,
                coeffs,
                factor,
            )
        },
        | Expr::Variable(v)
            if v == var =>
        {

            let entry = coeffs
                .entry(1)
                .or_insert_with(|| {

                    Expr::Constant(0.0)
                });

            *entry = simplify(
                &Expr::new_add(
                    entry.clone(),
                    factor.clone(),
                ),
            );

            Some(())
        },
        | Expr::Power(b, e) => {

            if let (
                Expr::Variable(v),
                Expr::Constant(p),
            ) = (&**b, &**e)
            {

                if v == var {

                    let degree =
                        p.to_u32()?;

                    let entry = coeffs
                        .entry(degree)
                        .or_insert_with(|| Expr::Constant(0.0));

                    *entry = simplify(
                        &Expr::new_add(
                            entry
                                .clone(
                                ),
                            factor
                                .clone(
                                ),
                        ),
                    );

                    return Some(());
                }
            }

            let entry = coeffs
                .entry(0)
                .or_insert_with(|| {

                    Expr::Constant(0.0)
                });

            *entry = simplify(
                &Expr::new_add(
                    entry.clone(),
                    Expr::new_mul(
                        expr.clone(),
                        factor.clone(),
                    ),
                ),
            );

            Some(())
        },
        | Expr::Add(a, b) => {

            collect_coeffs(
                a,
                var,
                coeffs,
                factor,
            )?;

            collect_coeffs(
                b,
                var,
                coeffs,
                factor,
            )
        },
        | Expr::Sub(a, b) => {

            collect_coeffs(
                a,
                var,
                coeffs,
                factor,
            )?;

            collect_coeffs(
                b,
                var,
                coeffs,
                &simplify(
                    &Expr::new_neg(
                        factor.clone(),
                    ),
                ),
            )
        },
        | Expr::Mul(a, b) => {
            if !contains_var(a, var) {

                collect_coeffs(
                    b,
                    var,
                    coeffs,
                    &simplify(
                        &Expr::new_mul(
                            factor
                                .clone(
                                ),
                            a.clone(),
                        ),
                    ),
                )
            } else if !contains_var(
                b, var,
            ) {

                collect_coeffs(
                    a,
                    var,
                    coeffs,
                    &simplify(
                        &Expr::new_mul(
                            factor
                                .clone(
                                ),
                            b.clone(),
                        ),
                    ),
                )
            } else {

                None
            }
        },
        | Expr::Neg(e) => {
            collect_coeffs(
                e,
                var,
                coeffs,
                &simplify(
                    &Expr::new_neg(
                        factor.clone(),
                    ),
                ),
            )
        },
        | _ if !contains_var(
            expr, var,
        ) =>
        {

            let entry = coeffs
                .entry(0)
                .or_insert_with(|| {

                    Expr::Constant(0.0)
                });

            *entry = simplify(
                &Expr::new_add(
                    entry.clone(),
                    Expr::new_mul(
                        expr.clone(),
                        factor.clone(),
                    ),
                ),
            );

            Some(())
        },
        | _ => None,
    }
}
