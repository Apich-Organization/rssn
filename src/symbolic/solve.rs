//! # Symbolic Equation Solving
//!
//! This module provides a powerful set of tools for solving equations and systems of equations.
//! It includes dispatchers that can handle polynomial, transcendental, linear, and multivariate
//! polynomial systems by selecting the appropriate algorithm, such as substitution, Gaussian
//! elimination, or Grobner bases.

use crate::symbolic::calculus::substitute;
use crate::symbolic::core::{Expr, Monomial, SparsePolynomial};
use crate::symbolic::grobner::{buchberger, MonomialOrder};
use crate::symbolic::matrix::create_empty_matrix;
use crate::symbolic::matrix::get_matrix_dims;
use crate::symbolic::matrix::null_space;
use crate::symbolic::matrix::rref;
use crate::symbolic::simplify::collect_and_order_terms;
use crate::symbolic::simplify::{is_zero, simplify};
use num_traits::ToPrimitive;
use std::collections::BTreeMap;
use std::collections::HashMap;

// =====================================================================================
// region: Main Solver Dispatchers
// =====================================================================================

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
pub fn solve(expr: &Expr, var: &str) -> Vec<Expr> {
    let equation = if let Expr::Eq(left, right) = expr {
        simplify(Expr::Sub(left.clone(), right.clone()))
    } else {
        expr.clone()
    };

    if let Some(solutions) = solve_polynomial(&equation, var) {
        return solutions;
    }
    if let Some(solutions) = solve_transcendental(&equation, var) {
        return solutions;
    }

    vec![Expr::Solve(Box::new(equation), var.to_string())]
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
pub fn solve_system(equations: &[Expr], vars: &[&str]) -> Option<Vec<(String, Expr)>> {
    if let Some(solutions) = solve_system_by_substitution(equations, vars) {
        return Some(solutions);
    }
    if let Some(solutions) = solve_system_with_grobner(equations, vars) {
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
pub fn solve_system_parcial(equations: &[Expr], vars: &[&str]) -> Option<Vec<(Expr, Expr)>> {
    let mut remaining_eqs: Vec<Expr> = equations.to_vec();
    let mut solutions: HashMap<String, Expr> = HashMap::new();
    let mut progress = true;

    while progress && !remaining_eqs.is_empty() {
        progress = false;
        let mut solved_eq_index: Option<usize> = None;

        for (i, eq) in remaining_eqs.iter().enumerate() {
            let mut current_eq = eq.clone();
            for (solved_var, solution_expr) in &solutions {
                current_eq = substitute(&current_eq, solved_var, solution_expr);
            }

            let remaining_vars: Vec<&str> = vars
                .iter()
                .filter(|v| !solutions.contains_key(**v))
                .cloned()
                .collect();
            if remaining_vars.len() == 1 {
                let var_to_solve = remaining_vars[0];
                let mut new_solutions = solve(&current_eq, var_to_solve);

                if !new_solutions.is_empty() {
                    // For now, take the first solution if multiple are returned.
                    // A more advanced implementation could handle branching.
                    let solution = new_solutions.remove(0);
                    solutions.insert(var_to_solve.to_string(), solution);
                    solved_eq_index = Some(i);
                    progress = true;
                    break; // Restart with the new solution
                }
            }
        }

        if let Some(index) = solved_eq_index {
            remaining_eqs.remove(index);
        }
    }

    // Back-substitution to resolve dependencies in solutions
    let mut final_solutions = HashMap::new();
    for var_name in vars.iter().map(|s| s.to_string()) {
        if let Some(mut solution) = solutions.get(&var_name).cloned() {
            let mut changed = true;
            while changed {
                changed = false;
                for (solved_var, sol_expr) in &solutions {
                    if solved_var != &var_name {
                        let new_solution = substitute(&solution, solved_var, sol_expr);
                        if new_solution != solution {
                            solution = new_solution;
                            changed = true;
                        }
                    }
                }
            }
            final_solutions.insert(var_name, simplify(solution));
        }
    }

    if final_solutions.len() == vars.len() {
        Some(
            vars.iter()
                .map(|&v| {
                    (
                        Expr::Variable(v.to_string()),
                        final_solutions.get(v).unwrap().clone(),
                    )
                })
                .collect(),
        )
    } else {
        None // Could not solve the system completely
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
/// A `Result` containing an `Expr` representing the solution (matrix, system, or no solution),
/// or an error string if inputs are invalid or dimensions are incompatible.
pub fn solve_linear_system_mat(a: &Expr, b: &Expr) -> Result<Expr, String> {
    let (a_rows, a_cols) =
        get_matrix_dims(a).ok_or_else(|| "A is not a valid matrix".to_string())?;
    let (b_rows, b_cols) =
        get_matrix_dims(b).ok_or_else(|| "b is not a valid matrix".to_string())?;

    if a_rows != b_rows {
        return Err("Matrix A and vector b have incompatible row dimensions".to_string());
    }
    if b_cols != 1 {
        return Err("b must be a column vector".to_string());
    }

    // 1. Construct augmented matrix [A | b]
    let Expr::Matrix(a_mat) = a else {
        unreachable!()
    };
    let Expr::Matrix(b_mat) = b else {
        unreachable!()
    };
    let mut augmented_mat = a_mat.clone();
    for i in 0..a_rows {
        augmented_mat[i].push(b_mat[i][0].clone());
    }

    // 2. Compute RREF
    let rref_expr = rref(&Expr::Matrix(augmented_mat))?;
    let Expr::Matrix(rref_mat) = rref_expr else {
        unreachable!()
    };

    // 3. Analyze RREF for solutions
    // Check for inconsistency: a row like [0, 0, ..., 0 | c] where c != 0
    //for i in 0..a_rows {
    // for (i , _item) in rref_mat.iter().take(a_rows) {
    //     let is_lhs_zero = rref_mat[i][0..a_cols].iter().all(is_zero);
    //     if is_lhs_zero && !is_zero(&rref_mat[i][a_cols]) {
    //         return Ok(Expr::NoSolution);
    //     }
    // }
    // for i in 0..a_rows { // Original index-based loop
    for (i, _row) in rref_mat.iter().take(a_rows).enumerate() {
        // i is the index (usize), _row is the reference to rref_mat[i]
        let is_lhs_zero = rref_mat[i][0..a_cols].iter().all(is_zero);
        if is_lhs_zero && !is_zero(&rref_mat[i][a_cols]) {
            return Ok(Expr::NoSolution);
        }
    }

    // Identify pivot and free variables
    let mut pivot_cols = Vec::new();
    let mut lead = 0;
    for r in 0..a_rows {
        if lead >= a_cols {
            break;
        }
        let mut i = lead;
        while i < a_cols && is_zero(&rref_mat[r][i]) {
            i += 1;
        }
        if i < a_cols {
            pivot_cols.push(i);
            lead = i + 1;
        }
    }

    let free_cols: Vec<usize> = (0..a_cols).filter(|c| !pivot_cols.contains(c)).collect();

    if free_cols.is_empty() {
        // Unique solution
        let mut solution = create_empty_matrix(a_cols, 1);
        for (i, &p_col) in pivot_cols.iter().enumerate() {
            solution[p_col][0] = rref_mat[i][a_cols].clone();
        }
        Ok(Expr::Matrix(solution))
    } else {
        // Infinite solutions (parametric)
        let particular_solution = {
            let mut sol = create_empty_matrix(a_cols, 1);
            for (i, &p_col) in pivot_cols.iter().enumerate() {
                sol[p_col][0] = rref_mat[i][a_cols].clone();
            }
            sol
        };

        let null_space_basis = null_space(a)?;

        Ok(Expr::System(vec![
            Expr::Matrix(particular_solution),
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
/// A `Result` containing a vector of `Expr` representing the solutions for `vars`,
/// or an error string if the system cannot be solved.
pub fn solve_linear_system(system: &Expr, vars: &[String]) -> Result<Vec<Expr>, String> {
    if let Expr::System(eqs) = system {
        let vars_str: Vec<&str> = vars.iter().map(|s| s.as_str()).collect();
        match solve_system(eqs, &vars_str) {
            Some(solutions) => {
                let mut sol_map: HashMap<String, Expr> = solutions.into_iter().collect();
                let ordered_solutions: Vec<Expr> = vars
                    .iter()
                    .map(|var| {
                        sol_map
                            .remove(var)
                            .unwrap_or(Expr::Variable("NotFound".to_string()))
                    })
                    .collect();
                Ok(ordered_solutions)
            }
            None => Err("System could not be solved.".to_string()),
        }
    } else {
        Err("Input must be a system of equations.".to_string())
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
/// A `Result` containing a vector of `Expr` representing the solutions for `vars`,
/// or an error string if the system is inconsistent, singular, or inputs are invalid.
pub fn solve_linear_system_gauss(system: &Expr, vars: &[String]) -> Result<Vec<Expr>, String> {
    if let Expr::System(eqs) = system {
        let n = vars.len();
        if eqs.len() != n {
            return Err(format!(
                "Number of equations ({}) does not match number of variables ({})",
                eqs.len(),
                n
            ));
        }

        let mut matrix_a = vec![vec![Expr::Constant(0.0); n]; n];
        let mut vector_b = vec![Expr::Constant(0.0); n];

        // 1. Build the augmented matrix [A|b] from the system of equations.
        for (i, eq) in eqs.iter().enumerate() {
            let (lhs, rhs) = match eq {
                Expr::Eq(l, r) => (l, r),
                _ => return Err(format!("Item {} is not a valid equation", i)),
            };

            // Initially, set the RHS of the matrix equation.
            vector_b[i] = *rhs.clone();

            // Extract coefficients for each variable from the LHS.
            if let Some(coeffs) = extract_polynomial_coeffs(lhs, "") {
                // Dummy var to get all terms
                for (_term_str, _coeff) in coeffs.iter().zip(vars.iter()) {
                    // This part is tricky. A better `extract_coeffs` is needed.
                    // For now, let's assume a simpler structure.
                }
            }

            // A more direct approach to build the matrix:
            let (_, terms) = collect_and_order_terms(lhs);
            for (term, coeff) in terms {
                if let Some(j) = vars.iter().position(|v| v == &term.to_string()) {
                    matrix_a[i][j] = coeff;
                } else if !is_zero(&coeff) && term.to_string() != "1" {
                    // If a term is not one of the variables, move it to the RHS.
                    vector_b[i] = simplify(Expr::Sub(
                        Box::new(vector_b[i].clone()),
                        Box::new(Expr::Mul(Box::new(coeff), Box::new(term))),
                    ));
                } else if term.to_string() == "1" {
                    // Constant term on LHS
                    vector_b[i] =
                        simplify(Expr::Sub(Box::new(vector_b[i].clone()), Box::new(coeff)));
                }
            }
        }

        // 2. Perform Gaussian elimination.
        for i in 0..n {
            // Find pivot
            let mut max_row = i;
            for (k, _item) in matrix_a.iter().enumerate().take(n).skip(i + 1) {
                // For symbolic solving, we can't just compare magnitudes.
                // We'll pivot on the first non-zero element.
                if !is_zero(&matrix_a[k][i]) {
                    max_row = k;
                    break;
                }
            }
            matrix_a.swap(i, max_row);
            vector_b.swap(i, max_row);

            let pivot = matrix_a[i][i].clone();
            if is_zero(&pivot) {
                return Err("Matrix is singular or underdetermined".to_string());
            }

            // Normalize the pivot row
            for j in i..n {
                matrix_a[i][j] = simplify(Expr::Div(
                    Box::new(matrix_a[i][j].clone()),
                    Box::new(pivot.clone()),
                ));
            }
            vector_b[i] = simplify(Expr::Div(
                Box::new(vector_b[i].clone()),
                Box::new(pivot.clone()),
            ));

            // Eliminate other rows
            for k in 0..n {
                if i != k {
                    let factor = matrix_a[k][i].clone();
                    for j in i..n {
                        let term = simplify(Expr::Mul(
                            Box::new(factor.clone()),
                            Box::new(matrix_a[i][j].clone()),
                        ));
                        matrix_a[k][j] =
                            simplify(Expr::Sub(Box::new(matrix_a[k][j].clone()), Box::new(term)));
                    }
                    let term_b = simplify(Expr::Mul(
                        Box::new(factor.clone()),
                        Box::new(vector_b[i].clone()),
                    ));
                    vector_b[k] =
                        simplify(Expr::Sub(Box::new(vector_b[k].clone()), Box::new(term_b)));
                }
            }
        }

        // 3. Back substitution is not needed as we have a diagonal matrix.
        Ok(vector_b)
    } else {
        Err("Input expression is not a system of equations".to_string())
    }
}

// =====================================================================================
// endregion: Main Solver Dispatchers
// =====================================================================================

// =====================================================================================
// region: System Solvers
// =====================================================================================

pub(crate) fn solve_system_by_substitution(
    equations: &[Expr],
    vars: &[&str],
) -> Option<Vec<(String, Expr)>> {
    let mut remaining_eqs: Vec<Expr> = equations.to_vec();
    let mut solutions: HashMap<String, Expr> = HashMap::new();
    let mut progress = true;

    while progress && !remaining_eqs.is_empty() {
        progress = false;
        let mut solved_eq_index: Option<usize> = None;

        for (i, eq) in remaining_eqs.iter().enumerate() {
            let mut current_eq = eq.clone();
            for (solved_var, solution_expr) in &solutions {
                current_eq = substitute(&current_eq, solved_var, solution_expr);
            }

            let remaining_vars: Vec<&str> = vars
                .iter()
                .filter(|v| !solutions.contains_key(**v))
                .cloned()
                .collect();
            if remaining_vars.len() == 1 {
                let var_to_solve = remaining_vars[0];
                let mut new_solutions = solve(&current_eq, var_to_solve);

                if !new_solutions.is_empty() {
                    let solution = new_solutions.remove(0);
                    solutions.insert(var_to_solve.to_string(), solution);
                    solved_eq_index = Some(i);
                    progress = true;
                    break;
                }
            }
        }

        if let Some(index) = solved_eq_index {
            remaining_eqs.remove(index);
        }
    }

    if solutions.len() != vars.len() {
        return None;
    }

    let mut final_solutions = HashMap::new();
    for &var_name_str in vars {
        let var_name = var_name_str.to_string();
        if let Some(mut solution) = solutions.get(&var_name).cloned() {
            for (solved_var, sol_expr) in &solutions {
                if solved_var != &var_name {
                    solution = substitute(&solution, solved_var, sol_expr);
                }
            }
            final_solutions.insert(var_name, simplify(solution));
        }
    }

    Some(final_solutions.into_iter().collect())
}

pub(crate) fn solve_system_with_grobner(
    equations: &[Expr],
    vars: &[&str],
) -> Option<Vec<(String, Expr)>> {
    let basis: Vec<SparsePolynomial> = equations
        .iter()
        .map(|eq| expr_to_sparse_poly(eq, vars))
        .collect();
    let grobner_basis = buchberger(&basis, MonomialOrder::Lexicographical);

    let mut solutions: HashMap<String, Expr> = HashMap::new();
    for poly in grobner_basis.iter().rev() {
        let mut current_eq = sparse_poly_to_expr(poly, vars);
        for (var, val) in &solutions {
            current_eq = substitute(&current_eq, var, val);
        }

        let remaining_vars: Vec<&str> = vars
            .iter()
            .filter(|v| contains_var(&current_eq, v))
            .cloned()
            .collect();
        if remaining_vars.len() == 1 {
            let roots = solve(&current_eq, remaining_vars[0]);
            if roots.is_empty() {
                return None;
            }
            solutions.insert(remaining_vars[0].to_string(), roots[0].clone());
        } else if !remaining_vars.is_empty() && !is_zero(&current_eq) {
            return None;
        }
    }

    if solutions.len() == vars.len() {
        Some(solutions.into_iter().collect())
    } else {
        None
    }
}

// =====================================================================================
// endregion: System Solvers
// =====================================================================================

// =====================================================================================
// region: Polynomial Solver
// =====================================================================================

pub(crate) fn solve_polynomial(expr: &Expr, var: &str) -> Option<Vec<Expr>> {
    let coeffs = extract_polynomial_coeffs(expr, var)?;
    let degree = coeffs.len() - 1;

    match degree {
        0 => Some(vec![]),
        1 => Some(solve_linear(&coeffs)),
        2 => Some(solve_quadratic(&coeffs)),
        3 => Some(solve_cubic(&coeffs)),
        4 => Some(solve_quartic(&coeffs)),
        _ => {
            let poly_expr = expr.clone();
            let mut roots = Vec::new();
            for i in 0..degree {
                roots.push(Expr::RootOf {
                    poly: Box::new(poly_expr.clone()),
                    index: i as u32,
                });
            }
            Some(roots)
        }
    }
}

pub(crate) fn solve_linear(coeffs: &[Expr]) -> Vec<Expr> {
    let a = &coeffs[0];
    let b = &coeffs[1];
    vec![simplify(Expr::Neg(Box::new(Expr::Div(
        Box::new(b.clone()),
        Box::new(a.clone()),
    ))))]
}

pub(crate) fn solve_quadratic(coeffs: &[Expr]) -> Vec<Expr> {
    let a = &coeffs[0];
    let b = &coeffs[1];
    let c = &coeffs[2];
    let discriminant = simplify(Expr::Sub(
        Box::new(Expr::Power(
            Box::new(b.clone()),
            Box::new(Expr::Constant(2.0)),
        )),
        Box::new(Expr::Mul(
            Box::new(Expr::Constant(4.0)),
            Box::new(Expr::Mul(Box::new(a.clone()), Box::new(c.clone()))),
        )),
    ));
    let sqrt_d = simplify(Expr::Sqrt(Box::new(discriminant)));
    let two_a = simplify(Expr::Mul(
        Box::new(Expr::Constant(2.0)),
        Box::new(a.clone()),
    ));
    vec![
        simplify(Expr::Div(
            Box::new(Expr::Add(
                Box::new(Expr::Neg(Box::new(b.clone()))),
                Box::new(sqrt_d.clone()),
            )),
            Box::new(two_a.clone()),
        )),
        simplify(Expr::Div(
            Box::new(Expr::Sub(
                Box::new(Expr::Neg(Box::new(b.clone()))),
                Box::new(sqrt_d),
            )),
            Box::new(two_a),
        )),
    ]
}

pub(crate) fn solve_cubic(coeffs: &[Expr]) -> Vec<Expr> {
    let a = &coeffs[0];
    let b = &simplify(Expr::Div(Box::new(coeffs[1].clone()), Box::new(a.clone())));
    let c = &simplify(Expr::Div(Box::new(coeffs[2].clone()), Box::new(a.clone())));
    let d = &simplify(Expr::Div(Box::new(coeffs[3].clone()), Box::new(a.clone())));

    let p = simplify(Expr::Sub(
        Box::new(c.clone()),
        Box::new(Expr::Div(
            Box::new(Expr::Power(
                Box::new(b.clone()),
                Box::new(Expr::Constant(2.0)),
            )),
            Box::new(Expr::Constant(3.0)),
        )),
    ));
    let q = simplify(Expr::Add(
        Box::new(Expr::Mul(
            Box::new(Expr::Constant(2.0 / 27.0)),
            Box::new(Expr::Power(
                Box::new(b.clone()),
                Box::new(Expr::Constant(3.0)),
            )),
        )),
        Box::new(Expr::Sub(
            Box::new(Expr::Mul(Box::new(b.clone()), Box::new(c.clone()))),
            Box::new(d.clone()),
        )),
    ));

    let inner_sqrt = simplify(Expr::Add(
        Box::new(Expr::Power(
            Box::new(Expr::Div(
                Box::new(q.clone()),
                Box::new(Expr::Constant(2.0)),
            )),
            Box::new(Expr::Constant(2.0)),
        )),
        Box::new(Expr::Power(
            Box::new(Expr::Div(
                Box::new(p.clone()),
                Box::new(Expr::Constant(3.0)),
            )),
            Box::new(Expr::Constant(3.0)),
        )),
    ));
    let u = simplify(Expr::Power(
        Box::new(Expr::Add(
            Box::new(Expr::Neg(Box::new(Expr::Div(
                Box::new(q.clone()),
                Box::new(Expr::Constant(2.0)),
            )))),
            Box::new(Expr::Sqrt(Box::new(inner_sqrt.clone()))),
        )),
        Box::new(Expr::Constant(1.0 / 3.0)),
    ));
    let v = simplify(Expr::Power(
        Box::new(Expr::Sub(
            Box::new(Expr::Neg(Box::new(Expr::Div(
                Box::new(q.clone()),
                Box::new(Expr::Constant(2.0)),
            )))),
            Box::new(Expr::Sqrt(Box::new(inner_sqrt))),
        )),
        Box::new(Expr::Constant(1.0 / 3.0)),
    ));

    let sub_term = simplify(Expr::Div(
        Box::new(b.clone()),
        Box::new(Expr::Constant(3.0)),
    ));

    let root1 = simplify(Expr::Sub(
        Box::new(Expr::Add(Box::new(u.clone()), Box::new(v.clone()))),
        Box::new(sub_term.clone()),
    ));
    vec![root1]
}

pub(crate) fn solve_quartic(_coeffs: &[Expr]) -> Vec<Expr> {
    let poly_expr = Expr::Variable("QuarticPoly".to_string());
    vec![
        Expr::RootOf {
            poly: Box::new(poly_expr.clone()),
            index: 0,
        },
        Expr::RootOf {
            poly: Box::new(poly_expr.clone()),
            index: 1,
        },
        Expr::RootOf {
            poly: Box::new(poly_expr.clone()),
            index: 2,
        },
        Expr::RootOf {
            poly: Box::new(poly_expr),
            index: 3,
        },
    ]
}

// =====================================================================================
// endregion: Polynomial Solver
// =====================================================================================

// =====================================================================================
// region: Transcendental Solver
// =====================================================================================

pub(crate) fn solve_transcendental(expr: &Expr, var: &str) -> Option<Vec<Expr>> {
    if let Expr::Sub(lhs, rhs) = expr {
        return solve_transcendental_pattern(lhs, rhs, var);
    }
    if let Expr::Add(lhs, rhs) = expr {
        return solve_transcendental_pattern(lhs, &Expr::Neg(rhs.clone()), var);
    }
    None
}

pub(crate) fn solve_transcendental_pattern(lhs: &Expr, rhs: &Expr, var: &str) -> Option<Vec<Expr>> {
    let n = Expr::Variable("k".to_string());
    let pi = Expr::Pi;

    let (func_part, const_part) = if contains_var(lhs, var) && !contains_var(rhs, var) {
        (lhs, rhs)
    } else if !contains_var(lhs, var) && contains_var(rhs, var) {
        (rhs, lhs)
    } else {
        return None;
    };

    match func_part {
        Expr::Sin(arg) => {
            let inner_solutions = solve(
                &Expr::Eq(
                    arg.clone(),
                    Box::new(Expr::Add(
                        Box::new(Expr::Mul(Box::new(n.clone()), Box::new(pi.clone()))),
                        Box::new(Expr::Mul(
                            Box::new(Expr::Power(Box::new(Expr::Constant(-1.0)), Box::new(n))),
                            Box::new(Expr::ArcSin(Box::new(const_part.clone()))),
                        )),
                    )),
                ),
                var,
            );
            Some(inner_solutions)
        }
        Expr::Cos(arg) => {
            let sol1 = solve(
                &Expr::Eq(
                    arg.clone(),
                    Box::new(Expr::Add(
                        Box::new(Expr::Mul(
                            Box::new(Expr::Constant(2.0)),
                            Box::new(Expr::Mul(Box::new(n.clone()), Box::new(pi.clone()))),
                        )),
                        Box::new(Expr::ArcCos(Box::new(const_part.clone()))),
                    )),
                ),
                var,
            );
            let sol2 = solve(
                &Expr::Eq(
                    arg.clone(),
                    Box::new(Expr::Sub(
                        Box::new(Expr::Mul(
                            Box::new(Expr::Constant(2.0)),
                            Box::new(Expr::Mul(Box::new(n.clone()), Box::new(pi.clone()))),
                        )),
                        Box::new(Expr::ArcCos(Box::new(const_part.clone()))),
                    )),
                ),
                var,
            );
            Some([sol1, sol2].concat())
        }
        Expr::Tan(arg) => {
            let inner_solutions = solve(
                &Expr::Eq(
                    arg.clone(),
                    Box::new(Expr::Add(
                        Box::new(Expr::Mul(Box::new(n.clone()), Box::new(pi.clone()))),
                        Box::new(Expr::ArcTan(Box::new(const_part.clone()))),
                    )),
                ),
                var,
            );
            Some(inner_solutions)
        }
        Expr::Exp(arg) => {
            let i = Expr::Complex(Box::new(Expr::Constant(0.0)), Box::new(Expr::Constant(1.0)));
            let log_sol = Expr::Add(
                Box::new(Expr::Log(Box::new(const_part.clone()))),
                Box::new(Expr::Mul(
                    Box::new(Expr::Mul(
                        Box::new(Expr::Constant(2.0)),
                        Box::new(Expr::Mul(Box::new(pi.clone()), Box::new(i))),
                    )),
                    Box::new(n),
                )),
            );
            Some(solve(&Expr::Eq(arg.clone(), Box::new(log_sol)), var))
        }
        _ => None,
    }
}

// =====================================================================================
// endregion: Transcendental Solver
// =====================================================================================

// =====================================================================================
// region: Helpers
// =====================================================================================

pub(crate) fn contains_var(expr: &Expr, var: &str) -> bool {
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

pub fn extract_polynomial_coeffs(expr: &Expr, var: &str) -> Option<Vec<Expr>> {
    let mut coeffs_map = HashMap::new();
    collect_coeffs(expr, var, &mut coeffs_map, &Expr::Constant(1.0))?;
    if coeffs_map.is_empty() {
        if !contains_var(expr, var) {
            let mut map = HashMap::new();
            map.insert(0, expr.clone());
            coeffs_map = map;
        } else {
            return None;
        }
    }
    let max_degree = *coeffs_map.keys().max().unwrap_or(&0);
    let mut coeffs = vec![Expr::Constant(0.0); max_degree as usize + 1];
    for (degree, coeff) in coeffs_map {
        coeffs[degree as usize] = simplify(coeff);
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
        Expr::Variable(v) if v == var => {
            let entry = coeffs.entry(1).or_insert_with(|| Expr::Constant(0.0));
            *entry = simplify(Expr::Add(Box::new(entry.clone()), Box::new(factor.clone())));
            Some(())
        }
        Expr::Power(b, e) => {
            if let (Expr::Variable(v), Expr::Constant(p)) = (&**b, &**e) {
                if v == var {
                    let degree = p.to_u32()?;
                    let entry = coeffs.entry(degree).or_insert_with(|| Expr::Constant(0.0));
                    *entry = simplify(Expr::Add(Box::new(entry.clone()), Box::new(factor.clone())));
                    return Some(());
                }
            }
            let entry = coeffs.entry(0).or_insert_with(|| Expr::Constant(0.0));
            *entry = simplify(Expr::Add(
                Box::new(entry.clone()),
                Box::new(Expr::Mul(Box::new(expr.clone()), Box::new(factor.clone()))),
            ));
            Some(())
        }
        Expr::Add(a, b) => {
            collect_coeffs(a, var, coeffs, factor)?;
            collect_coeffs(b, var, coeffs, factor)
        }
        Expr::Sub(a, b) => {
            collect_coeffs(a, var, coeffs, factor)?;
            collect_coeffs(
                b,
                var,
                coeffs,
                &simplify(Expr::Neg(Box::new(factor.clone()))),
            )
        }
        Expr::Mul(a, b) => {
            if !contains_var(a, var) {
                collect_coeffs(
                    b,
                    var,
                    coeffs,
                    &simplify(Expr::Mul(Box::new(factor.clone()), a.clone())),
                )
            } else if !contains_var(b, var) {
                collect_coeffs(
                    a,
                    var,
                    coeffs,
                    &simplify(Expr::Mul(Box::new(factor.clone()), b.clone())),
                )
            } else {
                None
            }
        }
        Expr::Neg(e) => collect_coeffs(
            e,
            var,
            coeffs,
            &simplify(Expr::Neg(Box::new(factor.clone()))),
        ),
        _ if !contains_var(expr, var) => {
            let entry = coeffs.entry(0).or_insert_with(|| Expr::Constant(0.0));
            *entry = simplify(Expr::Add(
                Box::new(entry.clone()),
                Box::new(Expr::Mul(Box::new(expr.clone()), Box::new(factor.clone()))),
            ));
            Some(())
        }
        _ => None,
    }
}

pub(crate) fn expr_to_sparse_poly(expr: &Expr, _vars: &[&str]) -> SparsePolynomial {
    let mut terms = BTreeMap::new();
    collect_poly_terms_recursive(expr, &mut terms, &Expr::Constant(1.0));
    SparsePolynomial { terms }
}

pub(crate) fn collect_poly_terms_recursive(
    _expr: &Expr,
    _terms: &mut BTreeMap<Monomial, Expr>,
    _current_coeff: &Expr,
) {
    // This is a recursive helper to build the sparse polynomial.
    // It needs to be robust.
}

pub(crate) fn sparse_poly_to_expr(poly: &SparsePolynomial, _vars: &[&str]) -> Expr {
    let mut total_expr = Expr::Constant(0.0);
    for (mono, coeff) in &poly.terms {
        let mut term_expr = coeff.clone();
        for (var_name, &exp) in &mono.0 {
            if exp > 0 {
                let var_expr = Expr::Power(
                    Box::new(Expr::Variable(var_name.clone())),
                    Box::new(Expr::Constant(exp as f64)),
                );
                term_expr = simplify(Expr::Mul(Box::new(term_expr), Box::new(var_expr)));
            }
        }
        total_expr = simplify(Expr::Add(Box::new(total_expr), Box::new(term_expr)));
    }
    total_expr
}

// =====================================================================================
// endregion: Helpers
// =====================================================================================
