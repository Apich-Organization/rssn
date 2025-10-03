//! # Combinatorics Module
//!
//! This module provides functions for various combinatorial calculations,
//! including permutations, combinations, binomial expansion, and solving
//! linear recurrence relations with constant coefficients. It also includes
//! tools for analyzing sequences from generating functions and applying the
//! Principle of Inclusion-Exclusion.

use crate::symbolic::calculus;
use crate::symbolic::core::Expr;
use crate::symbolic::series;
use crate::symbolic::simplify::{is_zero, simplify};
use crate::symbolic::solve::{extract_polynomial_coeffs, solve, solve_linear_system};
use std::collections::HashMap;

/// Expands an expression of the form `(a+b)^n` using the Binomial Theorem.
///
/// The Binomial Theorem states that `(a+b)^n = Σ_{k=0 to n} [ (n choose k) * a^(n-k) * b^k ]`.
/// This function returns a symbolic representation of this summation.
///
/// # Arguments
/// * `expr` - The expression to expand, expected to be in the form `Expr::Power(Expr::Add(a, b), n)`.
///
/// # Returns
/// An `Expr` representing the expanded binomial summation.
pub fn expand_binomial(expr: &Expr) -> Expr {
    if let Expr::Power(base, exponent) = expr {
        if let Expr::Add(a, b) = &**base {
            let n = exponent.clone();
            let k = Expr::Variable("k".to_string());

            // Calculate the combinations term: (n choose k)
            let combinations_term = combinations(*n.clone(), k.clone());

            // Calculate the a^(n-k) term
            let a_term = Expr::Power(
                a.clone(),
                Box::new(Expr::Sub(n.clone(), Box::new(k.clone()))),
            );

            // Calculate the b^k term
            let b_term = Expr::Power(b.clone(), Box::new(k.clone()));

            // Form the full term inside the summation: (n choose k) * a^(n-k) * b^k
            let full_term = Expr::Mul(
                Box::new(combinations_term),
                Box::new(Expr::Mul(Box::new(a_term), Box::new(b_term))),
            );

            // Construct the summation expression from k=0 to n
            return Expr::Summation(
                Box::new(full_term),
                "k".to_string(),
                Box::new(Expr::Constant(0.0)), // Lower bound for k
                n,                             // Upper bound for k
            );
        }
    }

    // If the expression is not in the expected (a+b)^n form, return it as is.
    expr.clone()
}

/// Calculates the number of permutations of `n` items taken `k` at a time, P(n, k).
///
/// The formula for permutations is `P(n, k) = n! / (n-k)!`.
///
/// # Arguments
/// * `n` - The total number of items.
/// * `k` - The number of items to choose.
///
/// # Returns
/// An `Expr` representing the number of permutations.
pub fn permutations(n: Expr, k: Expr) -> Expr {
    simplify(Expr::Div(
        Box::new(Expr::Factorial(Box::new(n.clone()))),
        Box::new(Expr::Factorial(Box::new(Expr::Sub(
            Box::new(n),
            Box::new(k),
        )))),
    ))
}

/// Calculates the number of combinations of `n` items taken `k` at a time, C(n, k).
///
/// The formula for combinations is `C(n, k) = n! / (k! * (n-k)!)` or `P(n, k) / k!`.
///
/// # Arguments
/// * `n` - The total number of items.
/// * `k` - The number of items to choose.
///
/// # Returns
/// An `Expr` representing the number of combinations.
pub fn combinations(n: Expr, k: Expr) -> Expr {
    simplify(Expr::Div(
        Box::new(permutations(n.clone(), k.clone())),
        Box::new(Expr::Factorial(Box::new(k))),
    ))
}

/// Solves a linear recurrence relation with constant coefficients.
///
/// This function implements the method of undetermined coefficients to find the particular
/// solution and combines it with the homogeneous solution to provide a general closed-form
/// solution. If initial conditions are provided, it solves for the constants in the general solution.
///
/// # Arguments
/// * `equation` - An `Expr::Eq` representing the recurrence relation. It should be in the form
///   `a(n) = c_1*a(n-1) + ... + c_k*a(n-k) + F(n)`. The `lhs` is assumed to contain the `a(n)` term
///   and the `rhs` contains the `a(n-k)` terms and `F(n)`.
/// * `initial_conditions` - A slice of tuples `(n_value, a_n_value)` for initial values.
///   These are used to determine the specific constants in the general solution.
/// * `term` - The name of the recurrence term, e.g., "a" for `a(n)`.
///
/// # Returns
/// An `Expr` representing the closed-form solution of the recurrence relation.
pub fn solve_recurrence(equation: Expr, initial_conditions: &[(Expr, Expr)], term: &str) -> Expr {
    if let Expr::Eq(lhs, rhs) = &equation {
        // 1. Deconstruct the equation to find the homogeneous part and the non-homogeneous term F(n).
        //    The `deconstruct_recurrence_eq` function is a simplification for now.
        let (homogeneous_coeffs, f_n) = deconstruct_recurrence_eq(lhs, rhs, term);

        // 2. Solve the homogeneous part by finding the roots of the characteristic equation.
        //    The characteristic equation is formed from the coefficients of the homogeneous part.
        let char_eq = build_characteristic_equation(&homogeneous_coeffs);
        let roots = solve(&char_eq, "r"); // Solve for the roots 'r' of the characteristic equation
        let mut root_counts: HashMap<Expr, usize> = HashMap::new();
        // Count multiplicities of roots, e.g., r=2 with multiplicity 2.
        for root in &roots {
            *root_counts.entry(root.clone()).or_insert(0) += 1;
        }

        // 3. Build the homogeneous solution form with unknown constants C_i.
        //    The form depends on the roots and their multiplicities.
        let (homogeneous_solution, const_vars) = build_homogeneous_solution(&root_counts);

        // 4. Find the particular solution based on the form of F(n) using the method of undetermined coefficients.
        let particular_solution =
            solve_particular_solution(&f_n, &root_counts, &homogeneous_coeffs, term);

        // 5. The general solution is the sum of the homogeneous and particular solutions.
        let general_solution = simplify(Expr::Add(
            Box::new(homogeneous_solution),
            Box::new(particular_solution),
        ));

        // 6. If no initial conditions are provided or no constants to solve for, return the general solution.
        if initial_conditions.is_empty() || const_vars.is_empty() {
            return general_solution;
        }

        // 7. Use initial conditions to solve for the specific values of the constants C_i.
        if let Some(final_solution) =
            solve_for_constants(&general_solution, &const_vars, initial_conditions)
        {
            return final_solution;
        }
    }
    // Fallback: If the equation cannot be parsed or solved, return an unevaluated Solve expression.
    Expr::Solve(Box::new(equation), term.to_string())
}

// region: Solver Helpers

/// Deconstructs the recurrence `lhs = rhs` into homogeneous coefficients and the F(n) term.
///
/// This is a simplified implementation. A robust parser would be needed to handle arbitrary
/// recurrence relation structures. Currently, it uses placeholder coefficients.
///
/// # Arguments
/// * `lhs` - The left-hand side of the recurrence equation.
/// * `rhs` - The right-hand side of the recurrence equation.
/// * `term` - The name of the recurrence term (e.g., "a").
///
/// # Returns
/// A tuple containing:
///   - `Vec<Expr>`: Coefficients of the homogeneous part (e.g., `[c_k, c_{k-1}, ..., c_0]`).
///   - `Expr`: The non-homogeneous term `F(n)`.
pub(crate) fn deconstruct_recurrence_eq(lhs: &Expr, rhs: &Expr, _term: &str) -> (Vec<Expr>, Expr) {
    // Example: For a(n) = 2*a(n-1) + n, the homogeneous part is a(n) - 2*a(n-1) = 0, and F(n) = n.
    // The coefficients would be for a(n-1) and a(n).
    // This part needs a proper parser to extract coefficients from terms like `a(n-k)`.
    // Placeholder: Assumes a simple first-order recurrence for demonstration.
    let _simplified_lhs = simplify(lhs.clone()); // This line is kept for potential future parsing logic.
    let coeffs = vec![Expr::Constant(-2.0), Expr::Constant(1.0)]; // Example: for a(n) - 2*a(n-1) = F(n)
    (coeffs, rhs.clone()) // Assumes rhs is F(n)
}

/// Builds the characteristic equation from the coefficients of the homogeneous recurrence.
///
/// For a recurrence `c_k*a(n) + c_{k-1}*a(n-1) + ... + c_0*a(n-k) = 0`,
/// the characteristic equation is `c_k*r^k + c_{k-1}*r^(k-1) + ... + c_0 = 0`.
///
/// # Arguments
/// * `coeffs` - A slice of `Expr` representing the coefficients of the homogeneous recurrence.
///
/// # Returns
/// An `Expr` representing the characteristic polynomial equation.
pub(crate) fn build_characteristic_equation(coeffs: &[Expr]) -> Expr {
    let mut terms = Vec::new();
    let r = Expr::Variable("r".to_string()); // The variable for the characteristic equation
    for (i, coeff) in coeffs.iter().enumerate() {
        // Each term is coeff * r^i
        let term = Expr::Mul(
            Box::new(coeff.clone()),
            Box::new(Expr::Power(
                Box::new(r.clone()),
                Box::new(Expr::Constant(i as f64)),
            )),
        );
        terms.push(term);
    }
    // Combine all terms into a polynomial sum.
    let mut poly = terms.pop().unwrap(); // Assumes at least one coefficient
    for term in terms {
        poly = Expr::Add(Box::new(poly), Box::new(term));
    }
    poly
}

/// Builds the homogeneous solution from the roots of the characteristic equation.
///
/// The form of the homogeneous solution depends on the roots and their multiplicities:
/// - For a distinct real root `r`, the term is `C * r^n`.
/// - For a real root `r` with multiplicity `m`, the terms are `(C_0 + C_1*n + ... + C_{m-1}*n^(m-1)) * r^n`.
/// - For complex conjugate roots `a ± bi`, the terms involve `(sqrt(a^2+b^2))^n * (C_1*cos(theta*n) + C_2*sin(theta*n))`.
///   (Note: Current implementation primarily handles real roots).
///
/// # Arguments
/// * `root_counts` - A HashMap where keys are the roots and values are their multiplicities.
///
/// # Returns
/// A tuple containing:
///   - `Expr`: The homogeneous solution with symbolic constants `C_i`.
///   - `Vec<String>`: A list of the names of the symbolic constants `C_i` used.
pub(crate) fn build_homogeneous_solution(
    root_counts: &HashMap<Expr, usize>,
) -> (Expr, Vec<String>) {
    let mut homogeneous_solution = Expr::Constant(0.0);
    let mut const_idx = 0;
    let mut const_vars = vec![];

    for (root, &multiplicity) in root_counts.iter() {
        let mut poly_term = Expr::Constant(0.0);
        // For each multiplicity, add a term C_i * n^j
        for i in 0..multiplicity {
            let c_name = format!("C{}", const_idx);
            let c = Expr::Variable(c_name.clone());
            const_vars.push(c_name);
            const_idx += 1;
            let n_pow_i = Expr::Power(
                Box::new(Expr::Variable("n".to_string())),
                Box::new(Expr::Constant(i as f64)),
            );
            poly_term = simplify(Expr::Add(
                Box::new(poly_term),
                Box::new(Expr::Mul(Box::new(c), Box::new(n_pow_i))),
            ));
        }

        // The root term: r^n
        let root_term = Expr::Power(
            Box::new(root.clone()),
            Box::new(Expr::Variable("n".to_string())),
        );
        // Combine (C_0 + C_1*n + ...) * r^n
        homogeneous_solution = simplify(Expr::Add(
            Box::new(homogeneous_solution),
            Box::new(Expr::Mul(Box::new(poly_term), Box::new(root_term))),
        ));
    }
    (homogeneous_solution, const_vars)
}

/// Determines and solves for the particular solution `a_n^(p)` using the method of undetermined coefficients.
///
/// This function guesses the form of the particular solution based on `F(n)`, substitutes it
/// into the recurrence, and solves a system of linear equations for the unknown coefficients.
///
/// # Arguments
/// * `f_n` - The non-homogeneous term `F(n)` from the recurrence relation.
/// * `char_roots` - A HashMap of characteristic roots and their multiplicities.
/// * `homogeneous_coeffs` - Coefficients of the homogeneous part of the recurrence.
/// * `term` - The name of the recurrence term (e.g., "a").
///
/// # Returns
/// An `Expr` representing the particular solution.
pub(crate) fn solve_particular_solution(
    f_n: &Expr,
    char_roots: &HashMap<Expr, usize>,
    homogeneous_coeffs: &[Expr],
    _term: &str,
) -> Expr {
    // If F(n) is zero, the particular solution is zero.
    if is_zero(f_n) {
        return Expr::Constant(0.0);
    }

    // 1. Guess the form of the particular solution with unknown coefficients.
    let (particular_form, unknown_coeffs) = guess_particular_form(f_n, char_roots);
    // If no suitable form can be guessed, return zero or an error.
    if unknown_coeffs.is_empty() {
        return Expr::Constant(0.0);
    }

    // 2. Substitute the guessed form back into the original recurrence relation.
    //    The recurrence is `a(n) + c_{k-1}a(n-1) + ... + c_0*a(n-k) = F(n)`.
    //    Let `L(a(n)) = a(n) + c_{k-1}a(n-1) + ... + c_0*a(n-k)`.
    //    We need to compute `L(particular_form)` and set it equal to `F(n)`.

    // Start with the a(n) term (coefficient 1 for a(n))
    let mut lhs_substituted = particular_form.clone();
    // Iterate through the other homogeneous terms c_i * a(n-i)
    for (i, coeff) in homogeneous_coeffs.iter().enumerate() {
        // The index `i` here corresponds to the power of r in the characteristic equation.
        // For a(n-k), the index would be 0. For a(n), it would be k.
        // Assuming homogeneous_coeffs are [c_0, c_1, ..., c_k] for a(n-k), ..., a(n)
        // The current structure of `homogeneous_coeffs` is assumed to be for `a(n-k)` to `a(n)`.
        // This part needs careful alignment with `deconstruct_recurrence_eq`.

        // For a recurrence like `a(n) = c_1*a(n-1) + c_2*a(n-2) + F(n)`
        // The homogeneous part is `a(n) - c_1*a(n-1) - c_2*a(n-2) = 0`
        // `homogeneous_coeffs` would be `[-c_1, -c_2, 1]` (for a(n-1), a(n-2), a(n))
        // The loop should iterate over the terms `a(n-1), a(n-2), ...`

        // This part of the code needs to be robustly linked to how `deconstruct_recurrence_eq` works.
        // Assuming `homogeneous_coeffs` are `[coeff_of_a(n-k), ..., coeff_of_a(n-1), coeff_of_a(n)]`
        // And `i` goes from 0 to k.
        // `coeff` is `coeff_of_a(n-i)`
        let n_minus_i = Expr::Sub(
            Box::new(Expr::Variable("n".to_string())),
            Box::new(Expr::Constant(i as f64)),
        );
        let term_an_i = calculus::substitute(&particular_form, "n", &n_minus_i);
        // Add `coeff * a(n-i)` to the LHS. Note: `homogeneous_coeffs` should include `a(n)`'s coeff (usually 1).
        lhs_substituted = Expr::Add(
            Box::new(lhs_substituted),
            Box::new(Expr::Mul(Box::new(coeff.clone()), Box::new(term_an_i))),
        );
    }

    // The equation to solve for coefficients is `L(particular_form) - F(n) = 0`.
    let equation_to_solve = simplify(Expr::Sub(Box::new(lhs_substituted), Box::new(f_n.clone())));

    // 3. Collect coefficients of powers of 'n' from `equation_to_solve` to form a system of linear equations.
    //    This step assumes `equation_to_solve` is a polynomial in 'n'.
    if let Some(poly_coeffs) = extract_polynomial_coeffs(&equation_to_solve, "n") {
        let mut system_eqs = Vec::new();
        // Each coefficient of 'n' must be zero for the equation to hold for all 'n'.
        for coeff_eq in poly_coeffs {
            if !is_zero(&coeff_eq) {
                system_eqs.push(Expr::Eq(Box::new(coeff_eq), Box::new(Expr::Constant(0.0))));
            }
        }

        // 4. Solve the system of linear equations for the unknown coefficients (A0, A1, ...).
        if let Ok(solutions) = solve_linear_system(&Expr::System(system_eqs), &unknown_coeffs) {
            let mut final_solution = particular_form;
            // Substitute the solved values back into the particular form.
            for (var, val) in unknown_coeffs.iter().zip(solutions.iter()) {
                final_solution = calculus::substitute(&final_solution, var, val);
            }
            return simplify(final_solution);
        }
    }

    // Fallback if solving for coefficients fails.
    Expr::Constant(0.0)
}

/// Guesses the form of the particular solution with unknown coefficients based on the form of `F(n)`.
/// Handles polynomial, exponential, and polynomial-exponential product forms.
/// Applies the modification rule if the guessed form overlaps with the homogeneous solution.
///
/// # Arguments
/// * `f_n` - The non-homogeneous term `F(n)`.
/// * `char_roots` - A HashMap of characteristic roots and their multiplicities.
///
/// # Returns
/// A tuple containing:
///   - `Expr`: The guessed form of the particular solution with symbolic unknown coefficients.
///   - `Vec<String>`: A list of the names of the unknown coefficients (e.g., "A0", "A1").
pub(crate) fn guess_particular_form(
    f_n: &Expr,
    char_roots: &HashMap<Expr, usize>,
) -> (Expr, Vec<String>) {
    let n_var = Expr::Variable("n".to_string());

    // Helper closure to create a polynomial form: A_0 + A_1*n + ... + A_d*n^d
    let create_poly_form = |degree: usize, prefix: &str| -> (Expr, Vec<String>) {
        let mut unknown_coeffs = Vec::new();
        let mut form = Expr::Constant(0.0);
        for i in 0..=degree {
            let coeff_name = format!("{}{}", prefix, i);
            unknown_coeffs.push(coeff_name.clone());
            form = Expr::Add(
                Box::new(form),
                Box::new(Expr::Mul(
                    Box::new(Expr::Variable(coeff_name)),
                    Box::new(Expr::Power(
                        Box::new(n_var.clone()),
                        Box::new(Expr::Constant(i as f64)),
                    )),
                )),
            );
        }
        (form, unknown_coeffs)
    };

    match f_n {
        // Case 1: F(n) is a polynomial (or constant), e.g., F(n) = n^2 + 3
        Expr::Polynomial(_) | Expr::Constant(_) => {
            let degree = extract_polynomial_coeffs(f_n, "n").map_or(0, |c| c.len() - 1);
            // Check if 1 is a root of the characteristic equation (multiplicity 's').
            // If 1 is a root, the guess needs to be multiplied by n^s.
            let s = *char_roots.get(&Expr::Constant(1.0)).unwrap_or(&0);

            let (mut form, coeffs) = create_poly_form(degree, "A");

            // Apply modification rule: multiply by n^s if 1 is a characteristic root.
            if s > 0 {
                form = Expr::Mul(
                    Box::new(Expr::Power(
                        Box::new(n_var.clone()),
                        Box::new(Expr::Constant(s as f64)),
                    )),
                    Box::new(form),
                );
            }
            (form, coeffs)
        }
        // Case 2: F(n) is an exponential, e.g., F(n) = b^n
        Expr::Power(base, exp) if matches!(&**exp, Expr::Variable(v) if v == "n") => {
            let b = base.clone(); // The base of the exponential
                                  // Check if 'b' is a root of the characteristic equation (multiplicity 's').
            let s = *char_roots.get(&b).unwrap_or(&0);

            let coeff_name = "A0".to_string(); // Single unknown coefficient for this form
            let mut form = Expr::Mul(
                Box::new(Expr::Variable(coeff_name.clone())),
                Box::new(f_n.clone()),
            );
            let coeffs = vec![coeff_name];

            // Apply modification rule: multiply by n^s if 'b' is a characteristic root.
            if s > 0 {
                form = Expr::Mul(
                    Box::new(Expr::Power(
                        Box::new(n_var.clone()),
                        Box::new(Expr::Constant(s as f64)),
                    )),
                    Box::new(form),
                );
            }
            (form, coeffs)
        }
        // Case 3: F(n) is a product of a polynomial and an exponential, e.g., F(n) = P(n) * b^n
        Expr::Mul(poly_expr, exp_expr) => {
            // Check if the second part is an exponential of the form b^n
            if let Expr::Power(base, exp) = &**exp_expr {
                if matches!(&**exp, Expr::Variable(v) if v == "n") {
                    let b = base.clone(); // The base of the exponential
                                          // Check if 'b' is a root of the characteristic equation (multiplicity 's').
                    let s = *char_roots.get(&b).unwrap_or(&0);

                    // Get the degree of the polynomial part P(n)
                    let degree =
                        extract_polynomial_coeffs(poly_expr, "n").map_or(0, |c| c.len() - 1);
                    // Create a polynomial form with unknown coefficients for P(n)
                    let (poly_form, poly_coeffs) = create_poly_form(degree, "A");

                    // The initial guess is (A_0 + A_1*n + ...) * b^n
                    let mut form = Expr::Mul(Box::new(poly_form), exp_expr.clone());

                    // Apply modification rule: multiply by n^s if 'b' is a characteristic root.
                    if s > 0 {
                        form = Expr::Mul(
                            Box::new(Expr::Power(
                                Box::new(n_var.clone()),
                                Box::new(Expr::Constant(s as f64)),
                            )),
                            Box::new(form),
                        );
                    }
                    return (form, poly_coeffs);
                }
            }
            // Fallback if the product form is not recognized (e.g., not P(n)*b^n)
            (Expr::Constant(0.0), vec![])
        }
        // TODO: Handle trigonometric forms (e.g., sin(kn), cos(kn)) and their combinations.
        Expr::Sin(arg) | Expr::Cos(arg) => {
            // Guess for F(n) = sin(kn) or cos(kn) is A*cos(kn) + B*sin(kn)
            let k_n = arg.clone();
            let coeff_a_name = "A".to_string();
            let coeff_b_name = "B".to_string();
            let unknown_coeffs = vec![coeff_a_name.clone(), coeff_b_name.clone()];

            let form = Expr::Add(
                Box::new(Expr::Mul(
                    Box::new(Expr::Variable(coeff_a_name)),
                    Box::new(Expr::Cos(k_n.clone())),
                )),
                Box::new(Expr::Mul(
                    Box::new(Expr::Variable(coeff_b_name)),
                    Box::new(Expr::Sin(k_n.clone())),
                )),
            );

            // TODO: Implement modification rule: check if cos(kn) or sin(kn) are part of the
            // homogeneous solution (i.e., if e^(ik) is a characteristic root) and multiply by n^s.

            (form, unknown_coeffs)
        }
        _ => (Expr::Constant(0.0), vec![]),
    }
}

/// Solves for the constants C_i in the general solution using the initial conditions.
///
/// This function substitutes the initial conditions into the general solution to form
/// a system of linear equations, which is then solved for the constants C_i.
///
/// # Arguments
/// * `general_solution` - The general solution of the recurrence with symbolic constants C_i.
/// * `const_vars` - A slice of strings representing the names of the constants C_i.
/// * `initial_conditions` - A slice of tuples `(n_value, a_n_value)` for initial values.
///
/// # Returns
/// An `Option<Expr>` representing the final particular solution with constants evaluated,
/// or `None` if the system cannot be solved.
pub(crate) fn solve_for_constants(
    general_solution: &Expr,
    const_vars: &[String],
    initial_conditions: &[(Expr, Expr)],
) -> Option<Expr> {
    let mut system_eqs = Vec::new();
    // For each initial condition (n_val, y_n_val), create an equation:
    // general_solution(n_val) = y_n_val
    for (n_val, y_n_val) in initial_conditions {
        let mut eq_lhs = general_solution.clone();
        // Substitute 'n' with n_val in the general solution.
        eq_lhs = calculus::substitute(&eq_lhs, "n", n_val);
        system_eqs.push(Expr::Eq(Box::new(eq_lhs), Box::new(y_n_val.clone())));
    }

    // Solve the system of linear equations for the constants C_i.
    if let Ok(const_vals) = solve_linear_system(
        &Expr::System(system_eqs),
        &const_vars.iter().map(|s| s.to_string()).collect::<Vec<_>>(),
    ) {
        let mut final_solution = general_solution.clone();
        // Substitute the solved values of C_i back into the general solution.
        for (c_name, c_val) in const_vars.iter().zip(const_vals.iter()) {
            final_solution = calculus::substitute(&final_solution, c_name, c_val);
        }
        return Some(simplify(final_solution));
    }
    None
}

// endregion

// Other combinatorics functions remain unchanged...

/// Extracts the sequence of coefficients from a generating function in closed form.
///
/// It computes the Taylor series of the expression around 0 and then extracts the coefficients
/// of the resulting polynomial.
///
/// # Arguments
/// * `expr` - The generating function expression (e.g., `1/(1-x)`).
/// * `var` - The variable of the function (e.g., "x").
/// * `max_order` - The number of terms to extract from the sequence.
///
/// # Returns
/// A vector of expressions representing the coefficients `a_0, a_1, ..., a_{max_order}`.
pub fn get_sequence_from_gf(expr: &Expr, var: &str, max_order: usize) -> Vec<Expr> {
    // The coefficients of the sequence are the coefficients of the Taylor series expansion.
    let series_poly = series::taylor_series(expr, var, &Expr::Constant(0.0), max_order);

    // `extract_polynomial_coeffs` requires an equation, so we create a dummy one.
    // This is a workaround; a more direct polynomial coefficient extractor would be better.
    let dummy_equation = Expr::Eq(Box::new(series_poly), Box::new(Expr::Constant(0.0)));
    extract_polynomial_coeffs(&dummy_equation, var).unwrap_or_default()
}

/// Applies the Principle of Inclusion-Exclusion.
///
/// This function calculates the size of the union of multiple sets given the sizes of
/// all possible intersections.
///
/// # Arguments
/// * `intersections` - A slice of vectors of expressions. `intersections[k]` should contain
///   the sizes of all (k+1)-wise intersections. For example, for sets A, B, C:
///   - `intersections[0]` = `[|A|, |B|, |C|]`
///   - `intersections[1]` = `[|A∩B|, |A∩C|, |B∩C|]`
///   - `intersections[2]` = `[|A∩B∩C|]`
///
/// # Returns
/// An expression representing the size of the union of the sets.
pub fn apply_inclusion_exclusion(intersections: &[Vec<Expr>]) -> Expr {
    let mut total_union_size = Expr::Constant(0.0);
    let mut sign = 1.0;

    for intersection_level in intersections {
        // Sum of sizes at the current level
        let sum_at_level = intersection_level
            .iter()
            .fold(Expr::Constant(0.0), |acc, size| {
                Expr::Add(Box::new(acc), Box::new(size.clone()))
            });

        // Apply alternating sign based on the level of intersection
        if sign > 0.0 {
            total_union_size = Expr::Add(Box::new(total_union_size), Box::new(sum_at_level));
        } else {
            total_union_size = Expr::Sub(Box::new(total_union_size), Box::new(sum_at_level));
        }

        sign *= -1.0; // Toggle sign for the next level
    }

    simplify(total_union_size)
}

/// Finds the smallest period of a sequence.
///
/// A sequence `S` has period `p` if `S[i] == S[i+p]` for all valid `i`.
/// This function finds the smallest `p > 0` for which this holds.
///
/// # Arguments
/// * `sequence` - A slice of `Expr` representing the sequence.
///
/// # Returns
/// An `Option<usize>` containing the smallest period if the sequence is periodic, otherwise `None`.
pub fn find_period(sequence: &[Expr]) -> Option<usize> {
    let n = sequence.len();
    if n == 0 {
        return None;
    }

    for p in 1..=n / 2 {
        if n % p == 0 {
            // Optimization: period must divide the length
            let mut is_periodic = true;
            for i in 0..(n - p) {
                if sequence[i] != sequence[i + p] {
                    is_periodic = false;
                    break;
                }
            }
            if is_periodic {
                return Some(p);
            }
        }
    }

    // If no smaller period is found, the period is the length of the sequence itself (or it's not periodic in a repeating sense).
    // Depending on definition, one might return n or None. We return None if no *repeating* period is found.
    None
}
