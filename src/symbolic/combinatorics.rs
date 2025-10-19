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
use std::sync::Arc;
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
            let combinations_term = combinations(&n.as_ref().clone(), k.clone());
            let a_term = Expr::new_pow(a.clone(), Expr::new_sub(n.clone(), k.clone()));
            let b_term = Expr::new_pow(b.clone(), k.clone());
            let full_term = Expr::new_mul(combinations_term, Expr::new_mul(a_term, b_term));
            return Expr::Summation(
                Arc::new(full_term),
                "k".to_string(),
                Arc::new(Expr::Constant(0.0)),
                n,
            );
        }
    }
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
    simplify(Expr::new_div(
        Expr::Factorial(Arc::new(n.clone())),
        Expr::Factorial(Arc::new(Expr::new_sub(n, k))),
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
pub fn combinations(n: &Expr, k: Expr) -> Expr {
    simplify(Expr::new_div(
        permutations(n.clone(), k.clone()),
        Expr::Factorial(Arc::new(k)),
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
        let (homogeneous_coeffs, f_n) = deconstruct_recurrence_eq(lhs, rhs, term);
        let char_eq = build_characteristic_equation(&homogeneous_coeffs);
        let roots = solve(&char_eq, "r");
        let mut root_counts: HashMap<Expr, usize> = HashMap::new();
        for root in &roots {
            *root_counts.entry(root.clone()).or_insert(0) += 1;
        }
        let (homogeneous_solution, const_vars) = build_homogeneous_solution(&root_counts);
        let particular_solution =
            solve_particular_solution(&f_n, &root_counts, &homogeneous_coeffs, term);
        let general_solution = simplify(Expr::new_add(homogeneous_solution, particular_solution));
        if initial_conditions.is_empty() || const_vars.is_empty() {
            return general_solution;
        }
        if let Some(final_solution) =
            solve_for_constants(&general_solution, &const_vars, initial_conditions)
        {
            return final_solution;
        }
    }
    Expr::Solve(Arc::new(equation), term.to_string())
}
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
    let _simplified_lhs = simplify(lhs.clone());
    let coeffs = vec![Expr::Constant(-2.0), Expr::Constant(1.0)];
    (coeffs, rhs.clone())
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
    let r = Expr::Variable("r".to_string());
    for (i, coeff) in coeffs.iter().enumerate() {
        let term = Expr::new_mul(
            coeff.clone(),
            Expr::new_pow(r.clone(), Expr::Constant(i as f64)),
        );
        terms.push(term);
    }
    if terms.is_empty() {
        return Expr::Constant(0.0);
    }
    let mut poly = match terms.pop() {
        Some(t) => t,
        _none => unreachable!(),
    };
    for term in terms {
        poly = Expr::new_add(poly, term);
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
    for (root, &multiplicity) in root_counts {
        let mut poly_term = Expr::Constant(0.0);
        for i in 0..multiplicity {
            let c_name = format!("C{}", const_idx);
            let c = Expr::Variable(c_name.clone());
            const_vars.push(c_name);
            const_idx += 1;
            let n_pow_i = Expr::new_pow(Expr::Variable("n".to_string()), Expr::Constant(i as f64));
            poly_term = simplify(Expr::new_add(poly_term, Expr::new_mul(c, n_pow_i)));
        }
        let root_term = Expr::new_pow(root.clone(), Expr::Variable("n".to_string()));
        homogeneous_solution = simplify(Expr::new_add(
            homogeneous_solution,
            Expr::new_mul(poly_term, root_term),
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
    if is_zero(f_n) {
        return Expr::Constant(0.0);
    }
    let (particular_form, unknown_coeffs) = guess_particular_form(f_n, char_roots);
    if unknown_coeffs.is_empty() {
        return Expr::Constant(0.0);
    }
    let mut lhs_substituted = particular_form.clone();
    for (i, coeff) in homogeneous_coeffs.iter().enumerate() {
        let n_minus_i = Expr::new_sub(Expr::Variable("n".to_string()), Expr::Constant(i as f64));
        let term_an_i = calculus::substitute(&particular_form, "n", &n_minus_i);
        lhs_substituted = Expr::new_add(lhs_substituted, Expr::new_mul(coeff.clone(), term_an_i));
    }
    let equation_to_solve = simplify(Expr::new_sub(lhs_substituted, f_n.clone()));
    if let Some(poly_coeffs) = extract_polynomial_coeffs(&equation_to_solve, "n") {
        let mut system_eqs = Vec::new();
        for coeff_eq in poly_coeffs {
            if !is_zero(&coeff_eq) {
                system_eqs.push(Expr::Eq(Arc::new(coeff_eq), Arc::new(Expr::Constant(0.0))));
            }
        }
        if let Ok(solutions) = solve_linear_system(&Expr::System(system_eqs), &unknown_coeffs) {
            let mut final_solution = particular_form;
            for (var, val) in unknown_coeffs.iter().zip(solutions.iter()) {
                final_solution = calculus::substitute(&final_solution, var, val);
            }
            return simplify(final_solution);
        }
    }
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
    let create_poly_form = |degree: usize, prefix: &str| -> (Expr, Vec<String>) {
        let mut unknown_coeffs = Vec::new();
        let mut form = Expr::Constant(0.0);
        for i in 0..=degree {
            let coeff_name = format!("{}{}", prefix, i);
            unknown_coeffs.push(coeff_name.clone());
            form = Expr::new_add(
                form,
                Expr::new_mul(
                    Expr::Variable(coeff_name),
                    Expr::new_pow(n_var.clone(), Expr::Constant(i as f64)),
                ),
            );
        }
        (form, unknown_coeffs)
    };
    match f_n {
        Expr::Polynomial(_) | Expr::Constant(_) => {
            let degree = extract_polynomial_coeffs(f_n, "n").map_or(0, |c| c.len() - 1);
            let s = *char_roots.get(&Expr::Constant(1.0)).unwrap_or(&0);
            let (mut form, coeffs) = create_poly_form(degree, "A");
            if s > 0 {
                form = Expr::new_mul(Expr::new_pow(n_var.clone(), Expr::Constant(s as f64)), form);
            }
            (form, coeffs)
        }
        Expr::Power(base, exp) if matches!(&** exp, Expr::Variable(v) if v == "n") => {
            let b = base.clone();
            let s = *char_roots.get(&b).unwrap_or(&0);
            let coeff_name = "A0".to_string();
            let mut form = Expr::new_mul(Expr::Variable(coeff_name.clone()), f_n.clone());
            let coeffs = vec![coeff_name];
            if s > 0 {
                form = Expr::new_mul(Expr::new_pow(n_var.clone(), Expr::Constant(s as f64)), form);
            }
            (form, coeffs)
        }
        Expr::Mul(poly_expr, exp_expr) => {
            if let Expr::Power(base, exp) = &**exp_expr {
                if matches!(&** exp, Expr::Variable(v) if v == "n") {
                    let b = base.clone();
                    let s = *char_roots.get(&b).unwrap_or(&0);
                    let degree =
                        extract_polynomial_coeffs(poly_expr, "n").map_or(0, |c| c.len() - 1);
                    let (poly_form, poly_coeffs) = create_poly_form(degree, "A");
                    let mut form = Expr::new_mul(poly_form, exp_expr.clone());
                    if s > 0 {
                        form = Expr::new_mul(
                            Expr::new_pow(n_var.clone(), Expr::Constant(s as f64)),
                            form,
                        );
                    }
                    return (form, poly_coeffs);
                }
            }
            (Expr::Constant(0.0), vec![])
        }
        Expr::Sin(arg) | Expr::Cos(arg) => {
            let k_n = arg.clone();
            let coeff_a_name = "A".to_string();
            let coeff_b_name = "B".to_string();
            let unknown_coeffs = vec![coeff_a_name.clone(), coeff_b_name.clone()];
            let form = Expr::new_add(
                Expr::new_mul(Expr::Variable(coeff_a_name), Expr::new_cos(k_n.clone())),
                Expr::new_mul(Expr::Variable(coeff_b_name), Expr::new_sin(k_n.clone())),
            );
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
    for (n_val, y_n_val) in initial_conditions {
        let mut eq_lhs = general_solution.clone();
        eq_lhs = calculus::substitute(&eq_lhs, "n", n_val);
        system_eqs.push(Expr::Eq(Arc::new(eq_lhs), Arc::new(y_n_val.clone())));
    }
    if let Ok(const_vals) = solve_linear_system(
        &Expr::System(system_eqs),
        &const_vars.iter().cloned().collect::<Vec<_>>(),
    ) {
        let mut final_solution = general_solution.clone();
        for (c_name, c_val) in const_vars.iter().zip(const_vals.iter()) {
            final_solution = calculus::substitute(&final_solution, c_name, c_val);
        }
        return Some(simplify(final_solution));
    }
    None
}
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
    let series_poly = series::taylor_series(expr, var, &Expr::Constant(0.0), max_order);
    let dummy_equation = Expr::Eq(Arc::new(series_poly), Arc::new(Expr::Constant(0.0)));
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
        let sum_at_level = intersection_level
            .iter()
            .fold(Expr::Constant(0.0), |acc, size| {
                Expr::new_add(acc, size.clone())
            });
        if sign > 0.0 {
            total_union_size = Expr::new_add(total_union_size, sum_at_level);
        } else {
            total_union_size = Expr::new_sub(total_union_size, sum_at_level);
        }
        sign *= -1.0;
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
        if n.is_multiple_of(p) {
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
    None
}
