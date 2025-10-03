//! # Symbolic Polynomial Manipulation
//!
//! This module provides a suite of functions for the symbolic manipulation of polynomials.
//! It supports operations such as addition, multiplication, differentiation, and long division.
//! It also includes tools for analyzing polynomial properties like degree and leading coefficients,
//! and for converting between symbolic expressions and coefficient-based representations.

use crate::symbolic::core::{Expr, Monomial, SparsePolynomial};
use crate::symbolic::grobner::subtract_poly;
use crate::symbolic::real_roots::eval_expr;
use crate::symbolic::simplify::as_f64;
use crate::symbolic::simplify::{is_zero, simplify};
use num_bigint::BigInt;
use num_traits::{One, ToPrimitive, Zero};
use std::collections::BTreeMap;
use std::collections::HashMap;
use std::ops::Add;
use std::ops::Mul;
use std::ops::Neg;
use std::ops::Sub;

/// Adds two sparse polynomials.
///
/// It iterates through the terms of the second polynomial and adds its coefficients
/// to the corresponding terms of the first polynomial. If a term does not exist
/// in the first polynomial, it is inserted.
///
/// # Arguments
/// * `p1` - The first sparse polynomial.
/// * `p2` - The second sparse polynomial.
///
/// # Returns
/// A new `SparsePolynomial` representing the sum.
pub fn add_poly(p1: &SparsePolynomial, p2: &SparsePolynomial) -> SparsePolynomial {
    let mut result_terms = p1.terms.clone();

    for (monomial, coeff2) in &p2.terms {
        let coeff1 = result_terms
            .entry(monomial.clone())
            .or_insert(Expr::Constant(0.0));
        *coeff1 = Expr::Add(Box::new(coeff1.clone()), Box::new(coeff2.clone()));
    }

    SparsePolynomial {
        terms: result_terms,
    }
}

/// Multiplies two sparse polynomials.
///
/// This function computes the product by iterating through all pairs of terms
/// from the two input polynomials. For each pair `(m1, c1)` and `(m2, c2)`:
/// - The new coefficient is `c1 * c2`.
/// - The new monomial is formed by adding the exponents of the variables from `m1` and `m2`.
///
/// # Arguments
/// * `p1` - The first sparse polynomial.
/// * `p2` - The second sparse polynomial.
///
/// # Returns
/// A new `SparsePolynomial` representing the product.
pub fn mul_poly(p1: &SparsePolynomial, p2: &SparsePolynomial) -> SparsePolynomial {
    let mut result_terms: BTreeMap<Monomial, Expr> = BTreeMap::new();

    for (m1, c1) in &p1.terms {
        for (m2, c2) in &p2.terms {
            // Multiply coefficients
            let new_coeff = Expr::Mul(Box::new(c1.clone()), Box::new(c2.clone()));

            // Add exponents for the new monomial
            let mut new_mono_map = m1.0.clone();
            for (var, exp2) in &m2.0 {
                let exp1 = new_mono_map.entry(var.clone()).or_insert(0);
                *exp1 += exp2;
            }
            let new_mono = Monomial(new_mono_map);

            // Add the new term to the result
            let existing_coeff = result_terms.entry(new_mono).or_insert(Expr::Constant(0.0));
            *existing_coeff = Expr::Add(Box::new(existing_coeff.clone()), Box::new(new_coeff));
        }
    }

    SparsePolynomial {
        terms: result_terms,
    }
}

/// Differentiates a sparse polynomial with respect to a given variable.
///
/// It applies the power rule to each term in the polynomial. For a term `c * x^n`,
/// the derivative is `(c * n) * x^(n-1)`.
/// Terms not containing the variable are eliminated, as their derivative is zero.
///
/// # Arguments
/// * `p` - The sparse polynomial to differentiate.
/// * `var` - The name of the variable to differentiate with respect to.
///
/// # Returns
/// A new `SparsePolynomial` representing the derivative.
pub fn differentiate_poly(p: &SparsePolynomial, var: &str) -> SparsePolynomial {
    let mut result_terms: BTreeMap<Monomial, Expr> = BTreeMap::new();

    for (monomial, coeff) in &p.terms {
        if let Some(&exp) = monomial.0.get(var) {
            if exp > 0 {
                // New coefficient is old_coeff * exponent
                let new_coeff = Expr::Mul(
                    Box::new(coeff.clone()),
                    Box::new(Expr::Constant(exp as f64)),
                );

                // New monomial has the exponent of 'var' decreased by 1
                let mut new_mono_map = monomial.0.clone();
                if exp == 1 {
                    new_mono_map.remove(var);
                } else {
                    *new_mono_map.get_mut(var).unwrap() -= 1;
                }
                let new_mono = Monomial(new_mono_map);

                result_terms.insert(new_mono, new_coeff);
            }
        }
        // If the monomial doesn't contain the variable, its derivative is zero, so we do nothing.
    }

    SparsePolynomial {
        terms: result_terms,
    }
}

/// Checks if an expression tree contains a specific variable.
///
/// This function performs a pre-order traversal of the expression tree and returns
/// `true` as soon as it finds a `Expr::Variable` node with the specified name.
///
/// # Arguments
/// * `expr` - The expression to search within.
/// * `var` - The name of the variable to look for.
///
/// # Returns
/// `true` if the variable is found, `false` otherwise.
pub fn contains_var(expr: &Expr, var: &str) -> bool {
    let mut found = false;
    expr.pre_order_walk(&mut |e| {
        if let Expr::Variable(name) = e {
            if name == var {
                found = true;
            }
        }
    });
    found
}

/// Checks if a given expression is a polynomial with respect to a specific variable.
///
/// A an expression is considered a polynomial in `var` if it is composed of
/// sums, products, and non-negative integer powers of `var`. Division is only
/// permitted if the denominator is a constant expression (i.e., does not contain `var`).
/// Transcendental functions (sin, cos, log, etc.) of `var` are not permitted.
///
/// # Arguments
/// * `expr` - The expression to check.
/// * `var` - The variable to check for polynomial properties against.
///
/// # Returns
/// `true` if the expression is a polynomial in `var`, `false` otherwise.
pub fn is_polynomial(expr: &Expr, var: &str) -> bool {
    match expr {
        // Constants are polynomials of degree 0.
        Expr::Constant(_) | Expr::BigInt(_) | Expr::Rational(_) => true,
        // A variable is a polynomial. If it's `var`, it's degree 1. If not, it's a constant (degree 0).
        Expr::Variable(_) => true,

        // Arithmetic operations are polynomial if their operands are.
        Expr::Add(a, b) | Expr::Sub(a, b) | Expr::Mul(a, b) => {
            is_polynomial(a, var) && is_polynomial(b, var)
        }

        Expr::Div(a, b) => {
            // Division is only allowed if the denominator is a constant expression (does not contain `var`).
            is_polynomial(a, var) && !contains_var(b, var)
        }

        Expr::Power(base, exp) => {
            // For x^n, n must be a non-negative integer.
            if let Expr::BigInt(n) = &**exp {
                if n >= &BigInt::zero() {
                    return is_polynomial(base, var);
                }
            }
            // Any other kind of power is not a polynomial (e.g., x^2.5, x^y).
            false
        }
        Expr::Neg(a) => is_polynomial(a, var),

        // Transcendental functions of the variable are not polynomials.
        Expr::Sin(_)
        | Expr::Cos(_)
        | Expr::Tan(_)
        | Expr::Log(_)
        | Expr::Exp(_)
        | Expr::Sec(_)
        | Expr::Csc(_)
        | Expr::Cot(_) => !contains_var(expr, var),

        // Default to false for unhandled complex cases.
        _ => false,
    }
}

/// Calculates the degree of a polynomial expression with respect to a variable.
///
/// This function determines the highest power of `var` in the expression by
/// recursively analyzing the symbolic tree. It handles addition, subtraction,
/// multiplication, division, and powers.
///
/// - For `+` or `-`, the degree is the maximum of the operands' degrees.
/// - For `*`, the degree is the sum of the operands' degrees.
/// - For `/`, the degree is the difference of the operands' degrees.
///
/// # Arguments
/// * `expr` - The polynomial expression.
/// * `var` - The variable of interest.
///
/// # Returns
/// An `i64` representing the degree of the polynomial. Returns `-1` if the expression
/// is not a simple polynomial in the specified variable.
pub fn polynomial_degree(expr: &Expr, var: &str) -> i64 {
    let s_expr = simplify(expr.clone());
    match s_expr {
        Expr::Add(a, b) | Expr::Sub(a, b) => {
            std::cmp::max(polynomial_degree(&a, var), polynomial_degree(&b, var))
        }
        Expr::Mul(a, b) => polynomial_degree(&a, var) + polynomial_degree(&b, var),
        Expr::Div(a, b) => polynomial_degree(&a, var) - polynomial_degree(&b, var),
        Expr::Power(ref base, ref exp) => {
            if let (Expr::Variable(v), Expr::BigInt(n)) = (base.as_ref(), exp.as_ref()) {
                if v == var {
                    return n.to_i64().unwrap_or(0);
                }
            }
            if !contains_var(&s_expr, var) {
                0
            } else {
                -1
            } // Not a simple polynomial
        }
        Expr::Variable(name) if name == var => 1,
        _ => 0, // Constants, other variables
    }
}

/// Finds the leading coefficient of a polynomial expression with respect to a variable.
///
/// The leading coefficient is the coefficient of the term with the highest degree.
/// This function works by recursively traversing the symbolic tree and determining
/// the leading coefficient based on the operation.
///
/// # Arguments
/// * `expr` - The polynomial expression.
/// * `var` - The variable of interest.
///
/// # Returns
/// An `Expr` representing the leading coefficient.
pub fn leading_coefficient(expr: &Expr, var: &str) -> Expr {
    let s_expr = simplify(expr.clone());
    match s_expr {
        Expr::Add(a, b) => {
            let deg_a = polynomial_degree(&a, var);
            let deg_b = polynomial_degree(&b, var);
            if deg_a > deg_b {
                leading_coefficient(&a, var)
            } else if deg_b > deg_a {
                leading_coefficient(&b, var)
            } else {
                simplify(Expr::Add(
                    Box::new(leading_coefficient(&a, var)),
                    Box::new(leading_coefficient(&b, var)),
                ))
            }
        }
        Expr::Sub(a, b) => {
            let deg_a = polynomial_degree(&a, var);
            let deg_b = polynomial_degree(&b, var);
            if deg_a > deg_b {
                leading_coefficient(&a, var)
            } else if deg_b > deg_a {
                simplify(Expr::Neg(Box::new(leading_coefficient(&b, var))))
            } else {
                simplify(Expr::Sub(
                    Box::new(leading_coefficient(&a, var)),
                    Box::new(leading_coefficient(&b, var)),
                ))
            }
        }
        Expr::Mul(a, b) => simplify(Expr::Mul(
            Box::new(leading_coefficient(&a, var)),
            Box::new(leading_coefficient(&b, var)),
        )),
        Expr::Div(a, b) => simplify(Expr::Div(
            Box::new(leading_coefficient(&a, var)),
            Box::new(leading_coefficient(&b, var)),
        )),
        Expr::Power(base, exp) => {
            if let (Expr::Variable(v), Expr::BigInt(_)) = (&*base, &*exp) {
                if v == var {
                    return Expr::BigInt(BigInt::one());
                }
            }
            simplify(Expr::Power(Box::new(leading_coefficient(&base, var)), exp))
        }
        Expr::Variable(name) if name == var => Expr::BigInt(BigInt::one()),
        _ => s_expr, // It's its own leading coefficient
    }
}

/// Performs polynomial long division on two expressions by symbolic manipulation.
///
/// This function implements the classic long division algorithm for polynomials.
/// It repeatedly subtracts multiples of the divisor from the remainder until the
/// remainder's degree is less than the divisor's degree.
///
/// # Arguments
/// * `n` - The numerator expression (the dividend).
/// * `d` - The denominator expression (the divisor).
/// * `var` - The variable of the polynomials.
///
/// # Returns
/// A tuple `(quotient, remainder)` where both are `Expr`.
pub fn polynomial_long_division(n: &Expr, d: &Expr, var: &str) -> (Expr, Expr) {
    pub(crate) fn is_zero_local(expr: &Expr) -> bool {
        match expr {
            Expr::Constant(c) => *c == 0.0,
            Expr::BigInt(i) => i.is_zero(),
            Expr::Rational(r) => r.is_zero(),
            _ => false,
        }
    }

    let mut q = Expr::BigInt(BigInt::zero());
    let mut r = n.clone();
    let d_deg = polynomial_degree(d, var);

    if d_deg < 0 {
        return (Expr::BigInt(BigInt::zero()), r);
    }

    let mut r_deg = polynomial_degree(&r, var);
    let mut iterations = 0;

    while r_deg >= d_deg && !is_zero_local(&r) {
        let lead_r = leading_coefficient(&r, var);
        let lead_d = leading_coefficient(d, var);

        let t_deg = r_deg - d_deg;
        let t_coeff = simplify(Expr::Div(Box::new(lead_r), Box::new(lead_d)));

        let t = if t_deg == 0 {
            t_coeff
        } else {
            simplify(Expr::Mul(
                Box::new(t_coeff),
                Box::new(Expr::Power(
                    Box::new(Expr::Variable(var.to_string())),
                    Box::new(Expr::BigInt(BigInt::from(t_deg))),
                )),
            ))
        };

        q = simplify(Expr::Add(Box::new(q.clone()), Box::new(t.clone())));
        let t_times_d = simplify(Expr::Mul(Box::new(t), Box::new(d.clone())));
        r = simplify(Expr::Sub(Box::new(r), Box::new(t_times_d)));
        let new_r_deg = polynomial_degree(&r, var);

        if new_r_deg >= r_deg {
            iterations += 1;
            if iterations > 5 {
                break;
            }
        } else {
            iterations = 0;
        }
        r_deg = new_r_deg;
    }

    (q, r)
}

// --- Alternate Implementations using Coefficient Vectors ---

/// Recursively collects coefficients of a polynomial expression into a map of degree -> coefficient.
pub(crate) fn collect_coeffs_recursive(expr: &Expr, var: &str) -> BTreeMap<u32, Expr> {
    match &simplify(expr.clone()) {
        Expr::Add(a, b) => {
            let mut map_a = collect_coeffs_recursive(a, var);
            let map_b = collect_coeffs_recursive(b, var);
            for (deg, coeff_b) in map_b {
                let coeff_a = map_a
                    .entry(deg)
                    .or_insert_with(|| Expr::BigInt(BigInt::zero()));
                *coeff_a = simplify(Expr::Add(Box::new(coeff_a.clone()), Box::new(coeff_b)));
            }
            map_a
        }
        Expr::Sub(a, b) => {
            let mut map_a = collect_coeffs_recursive(a, var);
            let map_b = collect_coeffs_recursive(b, var);
            for (deg, coeff_b) in map_b {
                let coeff_a = map_a
                    .entry(deg)
                    .or_insert_with(|| Expr::BigInt(BigInt::zero()));
                *coeff_a = simplify(Expr::Sub(Box::new(coeff_a.clone()), Box::new(coeff_b)));
            }
            map_a
        }
        Expr::Mul(a, b) => {
            let map_a = collect_coeffs_recursive(a, var);
            let map_b = collect_coeffs_recursive(b, var);
            let mut result_map = BTreeMap::new();
            for (deg_a, coeff_a) in &map_a {
                for (deg_b, coeff_b) in &map_b {
                    let new_deg = deg_a + deg_b;
                    let new_coeff_term = simplify(Expr::Mul(
                        Box::new(coeff_a.clone()),
                        Box::new(coeff_b.clone()),
                    ));
                    let entry = result_map
                        .entry(new_deg)
                        .or_insert_with(|| Expr::BigInt(BigInt::zero()));
                    *entry = simplify(Expr::Add(Box::new(entry.clone()), Box::new(new_coeff_term)));
                }
            }
            result_map
        }
        Expr::Power(base, exp) => {
            if let (Expr::Variable(v), Expr::BigInt(n)) = (base.as_ref(), exp.as_ref()) {
                if v == var {
                    let mut map = BTreeMap::new();
                    map.insert(n.to_u32().unwrap_or(0), Expr::BigInt(BigInt::one()));
                    return map;
                }
            }
            if !contains_var(base, var) {
                let mut map = BTreeMap::new();
                map.insert(0, expr.clone());
                return map;
            }
            BTreeMap::new()
        }
        Expr::Variable(v) if v == var => {
            let mut map = BTreeMap::new();
            map.insert(1, Expr::BigInt(BigInt::one()));
            map
        }
        Expr::Neg(a) => {
            let map_a = collect_coeffs_recursive(a, var);
            let mut result_map = BTreeMap::new();
            for (deg, coeff) in map_a {
                result_map.insert(deg, simplify(Expr::Neg(Box::new(coeff))));
            }
            result_map
        }
        e if !contains_var(e, var) => {
            let mut map = BTreeMap::new();
            map.insert(0, e.clone());
            map
        }
        _ => BTreeMap::new(),
    }
}

/// Converts a polynomial expression into a dense vector of its coefficients.
///
/// The coefficients are ordered from the highest degree term to the constant term.
/// This function first collects coefficients into a map from degree to expression,
/// then constructs a dense vector, filling in zero for any missing terms.
///
/// # Arguments
/// * `expr` - The polynomial expression.
/// * `var` - The variable of the polynomial.
///
/// # Returns
/// A `Vec<Expr>` containing the coefficients. Returns an empty vector if the
/// expression is not a valid polynomial.
pub fn to_polynomial_coeffs_vec(expr: &Expr, var: &str) -> Vec<Expr> {
    let map = collect_coeffs_recursive(expr, var);
    if map.is_empty() {
        if !contains_var(expr, var) {
            return vec![expr.clone()];
        }
        return vec![];
    }
    let max_deg = map.keys().max().cloned().unwrap_or(0);
    let mut result = vec![Expr::BigInt(BigInt::zero()); max_deg as usize + 1];
    for (deg, coeff) in map {
        result[deg as usize] = coeff;
    }
    result
}

/// Converts a dense vector of coefficients back into a polynomial expression.
///
/// The coefficient vector is assumed to be ordered from the constant term `c0`
/// to the highest degree term `cn`. The function constructs the expression
/// `c0 + c1*x + c2*x^2 + ... + cn*x^n`.
///
/// # Arguments
/// * `coeffs` - A slice of `Expr` representing the coefficients `[c0, c1, ...]`.
/// * `var` - The variable name for the polynomial.
///
/// # Returns
/// An `Expr` representing the constructed polynomial.
pub fn from_coeffs_to_expr(coeffs: &[Expr], var: &str) -> Expr {
    let mut expr = Expr::BigInt(BigInt::zero());
    for (i, coeff) in coeffs.iter().enumerate() {
        if !is_zero(&simplify(coeff.clone())) {
            let power = if i == 0 {
                Expr::BigInt(BigInt::one())
            } else {
                Expr::Power(
                    Box::new(Expr::Variable(var.to_string())),
                    Box::new(Expr::BigInt(BigInt::from(i))),
                )
            };
            let term = if i == 0 {
                coeff.clone()
            } else if let Expr::BigInt(b) = coeff {
                if b.is_one() {
                    power
                } else {
                    Expr::Mul(Box::new(coeff.clone()), Box::new(power))
                }
            } else {
                Expr::Mul(Box::new(coeff.clone()), Box::new(power))
            };
            expr = simplify(Expr::Add(Box::new(expr), Box::new(term)));
        }
    }
    expr
}

/// Performs polynomial long division using coefficient vectors.
///
/// This function provides an alternative to the symbolic `polynomial_long_division`.
/// It first converts the numerator and denominator expressions into dense coefficient
/// vectors and then performs the division algorithm on the vectors.
///
/// # Arguments
/// * `n` - The numerator expression.
/// * `d` - The denominator expression.
/// * `var` - The variable of the polynomials.
///
/// # Returns
/// A tuple `(quotient, remainder)` as `Expr`.
///
/// # Panics
/// Panics if the denominator is the zero polynomial.
pub fn polynomial_long_division_coeffs(n: &Expr, d: &Expr, var: &str) -> (Expr, Expr) {
    let mut num_coeffs = to_polynomial_coeffs_vec(n, var);
    let mut den_coeffs = to_polynomial_coeffs_vec(d, var);

    while den_coeffs
        .last()
        .is_some_and(|c| is_zero(&simplify(c.clone())))
    {
        den_coeffs.pop();
    }

    if den_coeffs.is_empty() {
        panic!("Polynomial division by zero");
    }

    let den_deg = den_coeffs.len() - 1;
    let mut num_deg = num_coeffs.len() - 1;

    if num_deg < den_deg {
        return (Expr::BigInt(BigInt::zero()), n.clone());
    }

    let lead_den = den_coeffs.last().unwrap().clone();
    let mut quot_coeffs = vec![Expr::BigInt(BigInt::zero()); num_deg - den_deg + 1];

    while num_deg >= den_deg {
        let lead_num = num_coeffs[num_deg].clone();
        let coeff = simplify(Expr::Div(Box::new(lead_num), Box::new(lead_den.clone())));

        let deg_diff = num_deg - den_deg;
        if deg_diff < quot_coeffs.len() {
            quot_coeffs[deg_diff] = coeff.clone();
        }

        for (i, _item) in den_coeffs.iter().enumerate().take(den_deg + 1) {
            if let Some(num_coeff) = num_coeffs.get_mut(deg_diff + i) {
                let term_to_sub = simplify(Expr::Mul(
                    Box::new(coeff.clone()),
                    Box::new(den_coeffs[i].clone()),
                ));
                *num_coeff = simplify(Expr::Sub(
                    Box::new(num_coeff.clone()),
                    Box::new(term_to_sub),
                ));
            }
        }

        while num_coeffs
            .last()
            .is_some_and(|c| is_zero(&simplify(c.clone())))
        {
            num_coeffs.pop();
        }

        if num_coeffs.is_empty() {
            num_deg = 0;
            let _help = num_deg;
            break;
        } else {
            num_deg = num_coeffs.len() - 1;
        }
    }

    let quotient = from_coeffs_to_expr(&quot_coeffs, var);
    let remainder = from_coeffs_to_expr(&num_coeffs, var);

    (quotient, remainder)
}

/// Converts a multivariate expression into a `SparsePolynomial` representation.
///
/// This function is designed to handle expressions with multiple variables, as specified
/// in the `vars` slice. It recursively processes the expression tree to identify terms
/// and their corresponding multivariate monomials.
///
/// # Arguments
/// * `expr` - The symbolic expression to convert.
/// * `vars` - A slice of variable names to be treated as parts of the polynomial's monomials.
///
/// # Returns
/// A `SparsePolynomial` representing the multivariate expression.
pub fn expr_to_sparse_poly(expr: &Expr, vars: &[&str]) -> SparsePolynomial {
    let mut terms = BTreeMap::new();
    collect_terms_recursive(expr, vars, &mut terms);
    SparsePolynomial { terms }
}

pub(crate) fn collect_terms_recursive(
    expr: &Expr,
    vars: &[&str],
    terms: &mut BTreeMap<Monomial, Expr>,
) {
    match &simplify(expr.clone()) {
        Expr::Add(a, b) => {
            collect_terms_recursive(a, vars, terms);
            collect_terms_recursive(b, vars, terms);
        }
        Expr::Sub(a, b) => {
            collect_terms_recursive(a, vars, terms);
            let mut neg_terms = BTreeMap::new();
            collect_terms_recursive(b, vars, &mut neg_terms);
            for (mono, coeff) in neg_terms {
                let entry = terms.entry(mono).or_insert_with(|| Expr::Constant(0.0));
                *entry = simplify(Expr::Sub(Box::new(entry.clone()), Box::new(coeff)));
            }
        }
        Expr::Mul(a, b) => {
            let mut p1_terms = BTreeMap::new();
            collect_terms_recursive(a, vars, &mut p1_terms);
            let mut p2_terms = BTreeMap::new();
            collect_terms_recursive(b, vars, &mut p2_terms);
            let p1 = SparsePolynomial { terms: p1_terms };
            let p2 = SparsePolynomial { terms: p2_terms };
            let product = p1 * p2;
            for (mono, coeff) in product.terms {
                let entry = terms.entry(mono).or_insert_with(|| Expr::Constant(0.0));
                *entry = simplify(Expr::Add(Box::new(entry.clone()), Box::new(coeff)));
            }
        }
        Expr::Power(base, exp) => {
            if let Some(e) = as_f64(exp) {
                if e.fract() == 0.0 && e >= 0.0 {
                    let mut p_base_terms = BTreeMap::new();
                    collect_terms_recursive(base, vars, &mut p_base_terms);
                    let p_base = SparsePolynomial {
                        terms: p_base_terms,
                    };
                    let mut result = SparsePolynomial {
                        terms: BTreeMap::from([(Monomial(BTreeMap::new()), Expr::Constant(1.0))]),
                    };
                    for _ in 0..(e as u32) {
                        result = result * p_base.clone();
                    }
                    for (mono, coeff) in result.terms {
                        let entry = terms.entry(mono).or_insert_with(|| Expr::Constant(0.0));
                        *entry = simplify(Expr::Add(Box::new(entry.clone()), Box::new(coeff)));
                    }
                    return;
                }
            }
            // Fallback for non-integer powers
            add_term(expr, &Expr::Constant(1.0), terms, vars);
        }
        Expr::Neg(a) => {
            let mut neg_terms = BTreeMap::new();
            collect_terms_recursive(a, vars, &mut neg_terms);
            for (mono, coeff) in neg_terms {
                let entry = terms.entry(mono).or_insert_with(|| Expr::Constant(0.0));
                *entry = simplify(Expr::Sub(Box::new(entry.clone()), Box::new(coeff)));
            }
        }
        _ => {
            add_term(expr, &Expr::Constant(1.0), terms, vars);
        }
    }
}

pub(crate) fn add_term(
    expr: &Expr,
    factor: &Expr,
    terms: &mut BTreeMap<Monomial, Expr>,
    vars: &[&str],
) {
    let mut is_poly_in_vars = false;
    for var in vars {
        if contains_var(expr, var) {
            is_poly_in_vars = true;
            break;
        }
    }

    if !is_poly_in_vars {
        // Treat as part of the coefficient
        let entry = terms
            .entry(Monomial(BTreeMap::new()))
            .or_insert(Expr::Constant(0.0));
        *entry = simplify(Expr::Add(
            Box::new(entry.clone()),
            Box::new(Expr::Mul(Box::new(factor.clone()), Box::new(expr.clone()))),
        ));
        return;
    }

    if let Expr::Variable(v) = expr {
        if vars.contains(&v.as_str()) {
            let mut mono_map = BTreeMap::new();
            mono_map.insert(v.clone(), 1);
            let entry = terms
                .entry(Monomial(mono_map))
                .or_insert(Expr::Constant(0.0));
            *entry = simplify(Expr::Add(Box::new(entry.clone()), Box::new(factor.clone())));
            return;
        }
    }
    // Fallback for complex terms
    let entry = terms
        .entry(Monomial(BTreeMap::new()))
        .or_insert(Expr::Constant(0.0));
    *entry = simplify(Expr::Add(
        Box::new(entry.clone()),
        Box::new(Expr::Mul(Box::new(factor.clone()), Box::new(expr.clone()))),
    ));
}

impl Neg for SparsePolynomial {
    type Output = Self;
    fn neg(self) -> Self {
        let mut new_terms = BTreeMap::new();
        for (mono, coeff) in self.terms {
            new_terms.insert(mono, simplify(Expr::Neg(Box::new(coeff))));
        }
        SparsePolynomial { terms: new_terms }
    }
}

impl SparsePolynomial {
    pub fn eval(&self, vars: &HashMap<String, f64>) -> f64 {
        self.terms
            .iter()
            .map(|(mono, coeff)| {
                let coeff_val = eval_expr(coeff, vars);
                let mono_val = mono.0.iter().fold(1.0, |acc, (var, exp)| {
                    let val = vars.get(var).cloned().unwrap_or(0.0);
                    acc * val.powi(*exp as i32)
                });
                coeff_val * mono_val
            })
            .sum()
    }
}

/// Multiplies a sparse polynomial by a scalar expression.
///
/// This function iterates through each term of the polynomial and multiplies its
/// coefficient by the given scalar expression. The monomials of the polynomial are unchanged.
///
/// # Arguments
/// * `poly` - The sparse polynomial.
/// * `scalar` - The scalar expression to multiply by.
///
/// # Returns
/// A new `SparsePolynomial` which is the result of the scalar multiplication.
pub fn poly_mul_scalar_expr(poly: &SparsePolynomial, scalar: &Expr) -> SparsePolynomial {
    let mut new_terms = BTreeMap::new();
    for (mono, coeff) in &poly.terms {
        new_terms.insert(
            mono.clone(),
            simplify(Expr::Mul(Box::new(coeff.clone()), Box::new(scalar.clone()))),
        );
    }
    SparsePolynomial { terms: new_terms }
}

/// Computes the greatest common divisor (GCD) of two sparse, single-variable polynomials.
///
/// This function uses the Euclidean algorithm, adapted for polynomials. It repeatedly
/// replaces the larger polynomial with the remainder of the division of the two polynomials
/// until the remainder is zero. The last non-zero remainder is the GCD.
///
/// # Arguments
/// * `a` - The first polynomial.
/// * `b` - The second polynomial.
/// * `var` - The variable of the polynomials.
///
/// # Returns
/// A new `SparsePolynomial` representing the greatest common divisor.
pub fn gcd(a: SparsePolynomial, b: SparsePolynomial, var: &str) -> SparsePolynomial {
    if b.terms.is_empty() {
        a
    } else {
        gcd(b.clone(), a.long_division(b, var).1, var)
    }
}

pub(crate) fn is_divisible(m1: &Monomial, m2: &Monomial) -> bool {
    m2.0.iter()
        .all(|(var, exp2)| m1.0.get(var).is_some_and(|exp1| exp1 >= exp2))
}

pub(crate) fn subtract_monomials(m1: &Monomial, m2: &Monomial) -> Monomial {
    let mut result = m1.0.clone();
    for (var, exp2) in &m2.0 {
        let exp1 = result.entry(var.clone()).or_insert(0);
        *exp1 -= exp2;
    }
    Monomial(result.into_iter().filter(|(_, exp)| *exp > 0).collect())
}

impl SparsePolynomial {
    pub fn degree(&self, var: &str) -> isize {
        self.terms
            .keys()
            .map(|m| m.0.get(var).cloned().unwrap_or(0) as isize)
            .max()
            .unwrap_or(-1)
    }

    pub fn leading_term(&self, var: &str) -> Option<(Monomial, Expr)> {
        self.terms
            .iter()
            .max_by_key(|(m, _)| m.0.get(var).cloned().unwrap_or(0))
            .map(|(m, c)| (m.clone(), c.clone()))
    }

    pub fn long_division(self, divisor: Self, var: &str) -> (Self, Self) {
        if divisor.terms.is_empty() {
            return (
                SparsePolynomial {
                    terms: BTreeMap::new(),
                },
                self,
            );
        }
        let mut quotient = SparsePolynomial {
            terms: BTreeMap::new(),
        };
        let mut remainder = self;
        let divisor_deg = divisor.degree(var);

        while remainder.degree(var) >= divisor_deg {
            let (lm_d, lc_d) = match divisor.leading_term(var) {
                Some(term) => term,
                None => break, // Divisor is zero
            };
            let (lm_r, lc_r) = match remainder.leading_term(var) {
                Some(term) => term,
                None => break, // Remainder is zero
            };

            if !is_divisible(&lm_r, &lm_d) {
                break;
            }

            let t_coeff = simplify(Expr::Div(Box::new(lc_r), Box::new(lc_d.clone())));
            let t_mono = subtract_monomials(&lm_r, &lm_d);
            let mut t = SparsePolynomial {
                terms: BTreeMap::new(),
            };
            t.terms.insert(t_mono, t_coeff);

            quotient = add_poly(&quotient, &t);
            let sub_term = mul_poly(&t, &divisor);
            remainder = subtract_poly(&remainder, &sub_term);
        }
        (quotient, remainder)
    }

    pub fn get_coeffs_as_vec(&self, var: &str) -> Vec<Expr> {
        let deg = self.degree(var);
        if deg < 0 {
            return vec![];
        }
        let mut coeffs = vec![Expr::Constant(0.0); (deg + 1) as usize];
        for (mono, coeff) in &self.terms {
            let d = mono.0.get(var).cloned().unwrap_or(0) as usize;
            if d < coeffs.len() {
                coeffs[d] = coeff.clone();
            }
        }
        coeffs.reverse();
        coeffs
    }

    pub fn get_coeff_for_power(&self, var: &str, power: usize) -> Option<Expr> {
        let mut mono_map = BTreeMap::new();
        if power > 0 {
            mono_map.insert(var.to_string(), power as u32);
        }
        self.terms.get(&Monomial(mono_map)).cloned()
    }
}

impl Add for SparsePolynomial {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        add_poly(&self, &rhs)
    }
}

impl Sub for SparsePolynomial {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        let neg_rhs = mul_poly(&rhs, &poly_from_coeffs(&[Expr::Constant(-1.0)], ""));
        add_poly(&self, &neg_rhs)
    }
}

impl Mul for SparsePolynomial {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        mul_poly(&self, &rhs)
    }
}

/*
impl SparsePolynomial {
    pub fn degree(&self, var: &str) -> isize {
        self.terms.keys()
            .map(|m| m.0.get(var).cloned().unwrap_or(0) as isize)
            .max().unwrap_or(-1)
    }

    pub fn get_coeffs_as_vec(&self, var: &str) -> Vec<Expr> {
        let deg = self.degree(var);
        if deg < 0 { return vec![]; }
        let mut coeffs = vec![Expr::Constant(0.0); (deg + 1) as usize];
        for (mono, coeff) in &self.terms {
            let d = mono.0.get(var).cloned().unwrap_or(0) as usize;
            if d < coeffs.len() { coeffs[d] = coeff.clone(); }
        }
        coeffs.reverse();
        coeffs
    }

    pub fn get_coeff_for_power(&self, var: &str, power: usize) -> Option<Expr> {
        let mut mono_map = BTreeMap::new();
        if power > 0 {
            mono_map.insert(var.to_string(), power as u32);
        }
        self.terms.get(&Monomial(mono_map)).cloned()
    }
}
*/

/// Creates a `SparsePolynomial` from a dense vector of coefficients.
///
/// The coefficients are assumed to be ordered from the highest degree term to the constant term
/// (e.g., `[c_n, c_{n-1}, ..., c_0]`).
///
/// # Arguments
/// * `coeffs` - A slice of `Expr` representing the coefficients.
/// * `var` - The variable name for the polynomial.
///
/// # Returns
/// A `SparsePolynomial` created from the coefficients.
pub fn poly_from_coeffs(coeffs: &[Expr], var: &str) -> SparsePolynomial {
    let mut terms = BTreeMap::new();
    let n = coeffs.len() - 1;
    for (i, coeff) in coeffs.iter().enumerate() {
        if !is_zero(&simplify(coeff.clone())) {
            let mut mono_map = BTreeMap::new();
            let power = (n - i) as u32;
            if power > 0 {
                mono_map.insert(var.to_string(), power);
            }
            terms.insert(Monomial(mono_map), coeff.clone());
        }
    }
    SparsePolynomial { terms }
}
/// Converts a sparse polynomial back into a symbolic expression.
///
/// This function iterates through the terms of the sparse polynomial and constructs
/// a symbolic `Expr` tree by summing up all `coefficient * monomial` terms.
///
/// # Arguments
/// * `poly` - The sparse polynomial to convert.
///
/// # Returns
/// An `Expr` representing the polynomial.
pub fn sparse_poly_to_expr(poly: &SparsePolynomial) -> Expr {
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

/*
impl Add for SparsePolynomial {
    type Output = Self;
    pub(crate) fn add(self, rhs: Self) -> Self {
        let mut result_terms = self.terms.clone();
        for (mono, coeff) in rhs.terms {
            let entry = result_terms.entry(mono).or_insert_with(|| Expr::Constant(0.0));
            *entry = simplify(Expr::Add(Box::new(entry.clone()), Box::new(coeff)));
        }
        SparsePolynomial { terms: result_terms }
    }
}

impl Sub for SparsePolynomial {
    type Output = Self;
    pub(crate) fn sub(self, rhs: Self) -> Self {
        let mut result_terms = self.terms.clone();
        for (mono, coeff) in rhs.terms {
            let entry = result_terms.entry(mono).or_insert_with(|| Expr::Constant(0.0));
            *entry = simplify(Expr::Sub(Box::new(entry.clone()), Box::new(coeff)));
        }
        SparsePolynomial { terms: result_terms }
    }
}

impl Mul for SparsePolynomial {
    type Output = Self;
    pub(crate) fn mul(self, rhs: Self) -> Self {
        let mut result_terms = BTreeMap::new();
        for (m1, c1) in &self.terms {
            for (m2, c2) in &rhs.terms {
                let new_coeff = simplify(Expr::Mul(Box::new(c1.clone()), Box::new(c2.clone())));
                let mut new_mono_map = m1.0.clone();
                for (var, exp2) in &m2.0 {
                    *new_mono_map.entry(var.clone()).or_insert(0) += exp2;
                }
                let new_mono = Monomial(new_mono_map);
                let entry = result_terms.entry(new_mono).or_insert_with(|| Expr::Constant(0.0));
                *entry = simplify(Expr::Add(Box::new(entry.clone()), Box::new(new_coeff)));
            }
        }
        SparsePolynomial { terms: result_terms }
    }
}
*/

/*
impl SparsePolynomial {
    pub fn degree(&self, var: &str) -> isize {
        self.terms.keys()
            .map(|m| m.0.get(var).cloned().unwrap_or(0) as isize)
            .max().unwrap_or(-1)
    }

    pub fn get_coeffs_as_vec(&self, var: &str) -> Vec<Expr> {
        let deg = self.degree(var);
        if deg < 0 { return vec![]; }
        let mut coeffs = vec![Expr::Constant(0.0); (deg + 1) as usize];
        for (mono, coeff) in &self.terms {
            let d = mono.0.get(var).cloned().unwrap_or(0) as usize;
            if d < coeffs.len() { coeffs[d] = coeff.clone(); }
        }
        coeffs.reverse();
        coeffs
    }

    pub fn get_coeff_for_power(&self, var: &str, power: usize) -> Option<Expr> {
        let mut mono_map = BTreeMap::new();
        if power > 0 {
            mono_map.insert(var.to_string(), power as u32);
        }
        self.terms.get(&Monomial(mono_map)).cloned()
    }
}
*/

/*
/// Creates a SparsePolynomial from a dense vector of coefficients (c_n, c_{n-1}, ..., c0).
pub fn poly_from_coeffs(coeffs: &[Expr], var: &str) -> SparsePolynomial {
    let mut terms = BTreeMap::new();
    let n = coeffs.len() - 1;
    for (i, coeff) in coeffs.iter().enumerate() {
        if !is_zero(&simplify(coeff.clone())) {
            let mut mono_map = BTreeMap::new();
            let power = (n - i) as u32;
            if power > 0 {
                mono_map.insert(var.to_string(), power);
            }
            terms.insert(Monomial(mono_map), coeff.clone());
        }
    }
    SparsePolynomial { terms }
}

/// Converts a sparse polynomial back into a symbolic expression.
pub fn sparse_poly_to_expr(poly: &SparsePolynomial) -> Expr {
    let mut total_expr = Expr::Constant(0.0);
    for (mono, coeff) in &poly.terms {
        let mut term_expr = coeff.clone();
        for (var_name, &exp) in &mono.0 {
            if exp > 0 {
                let var_expr = Expr::Power(Box::new(Expr::Variable(var_name.clone())), Box::new(Expr::Constant(exp as f64)));
                term_expr = simplify(Expr::Mul(Box::new(term_expr), Box::new(var_expr)));
            }
        }
        total_expr = simplify(Expr::Add(Box::new(total_expr), Box::new(term_expr)));
    }
    total_expr
}

/*
impl SparsePolynomial {
    pub fn degree(&self, var: &str) -> isize {
        self.terms.keys()
            .map(|m| m.0.get(var).cloned().unwrap_or(0) as isize)
            .max().unwrap_or(-1)
    }

    pub fn get_coeffs_as_vec(&self, var: &str) -> Vec<Expr> {
        let deg = self.degree(var);
        if deg < 0 { return vec![Expr::Constant(0.0)]; }
        let mut coeffs = vec![Expr::Constant(0.0); (deg + 1) as usize];
        for (mono, coeff) in &self.terms {
            let d = mono.0.get(var).cloned().unwrap_or(0) as usize;
            if d < coeffs.len() { coeffs[d] = coeff.clone(); }
        }
        coeffs.reverse();
        coeffs
    }

    pub fn get_coeff_for_power(&self, var: &str, power: usize) -> Option<Expr> {
        let mut mono_map = BTreeMap::new();
        if power > 0 {
            mono_map.insert(var.to_string(), power as u32);
        }
        self.terms.get(&Monomial(mono_map)).cloned()
    }
}
*/

impl Sub for SparsePolynomial {
    type Output = Self;
    pub(crate) fn sub(self, rhs: Self) -> Self {
        let neg_rhs = mul_poly(&rhs, &poly_from_coeffs(&[Expr::Constant(-1.0)], ""));
        add_poly(&self, &neg_rhs)
    }
}

/*
/// Creates a SparsePolynomial from a dense vector of coefficients (c_n, c_{n-1}, ..., c0).
pub fn poly_from_coeffs(coeffs: &[Expr], var: &str) -> SparsePolynomial {
    let mut terms = BTreeMap::new();
    let n = coeffs.len() - 1;
    for (i, coeff) in coeffs.iter().enumerate() {
        if !is_zero(&simplify(coeff.clone())) {
            let mut mono_map = BTreeMap::new();
            let power = (n - i) as u32;
            if power > 0 {
                mono_map.insert(var.to_string(), power);
            }
            terms.insert(Monomial(mono_map), coeff.clone());
        }
    }
    SparsePolynomial { terms }
}

/// Converts a sparse polynomial back into a symbolic expression.
pub fn sparse_poly_to_expr(poly: &SparsePolynomial) -> Expr {
    let mut total_expr = Expr::Constant(0.0);
    for (mono, coeff) in &poly.terms {
        let mut term_expr = coeff.clone();
        for (var_name, &exp) in &mono.0 {
            if exp > 0 {
                let var_expr = Expr::Power(Box::new(Expr::Variable(var_name.clone())), Box::new(Expr::Constant(exp as f64)));
                term_expr = simplify(Expr::Mul(Box::new(term_expr), Box::new(var_expr)));
            }
        }
        total_expr = simplify(Expr::Add(Box::new(total_expr), Box::new(term_expr)));
    }
    total_expr
}

impl SparsePolynomial {
    pub fn degree(&self, var: &str) -> isize {
        self.terms.keys()
            .map(|m| m.0.get(var).cloned().unwrap_or(0) as isize)
            .max().unwrap_or(-1)
    }

    pub fn get_coeffs_as_vec(&self, var: &str) -> Vec<Expr> {
        let deg = self.degree(var);
        if deg < 0 { return vec![Expr::Constant(0.0)]; }
        let mut coeffs = vec![Expr::Constant(0.0); (deg + 1) as usize];
        for (mono, coeff) in &self.terms {
            let d = mono.0.get(var).cloned().unwrap_or(0) as usize;
            coeffs[d] = coeff.clone();
        }
        coeffs.reverse();
        coeffs
    }

    pub fn get_coeff_for_power(&self, var: &str, power: usize) -> Option<Expr> {
        let mut mono_map = BTreeMap::new();
        if power > 0 {
            mono_map.insert(var.to_string(), power as u32);
        }
        self.terms.get(&Monomial(mono_map)).cloned()
    }
}

impl Sub for SparsePolynomial {
    type Output = Self;
    pub(crate) fn sub(self, rhs: Self) -> Self {
        let neg_rhs = mul_poly(&rhs, &poly_from_coeffs(&[Expr::Constant(-1.0)], ""));
        add_poly(&self, &neg_rhs)
    }
}

/// Creates a SparsePolynomial from a dense vector of coefficients (c0, c1, c2, ...).
pub fn poly_from_coeffs(coeffs: &[Expr], var: &str) -> SparsePolynomial {
    let mut terms = BTreeMap::new();
    for (i, coeff) in coeffs.iter().enumerate() {
        if !is_zero(&simplify(coeff.clone())) {
            let mut mono_map = BTreeMap::new();
            if i > 0 {
                mono_map.insert(var.to_string(), i as u32);
            }
            terms.insert(Monomial(mono_map), coeff.clone());
        }
    }
    SparsePolynomial { terms }
}

impl SparsePolynomial {
    /// Returns the coefficient of a specific power of a variable.
    pub fn get_coeff_for_power(&self, power: usize) -> Option<Expr> {
        // This is a simplified version for single-variable polynomials.
        let mono = Monomial(BTreeMap::from([(self.get_main_var()?, power as u32)]));
        self.terms.get(&mono).cloned()
    }

    /// Gets the main variable of a single-variable polynomial.
    pub(crate) fn get_main_var(&self) -> Option<String> {
        self.terms.keys().next()?.0.keys().next().cloned()
    }

    /// Converts a sparse polynomial to a dense vector of coefficients.
    pub fn get_coeffs_as_vec(&self, max_deg: usize) -> Vec<Expr> {
        let mut vec = vec![Expr::Constant(0.0); max_deg + 1];
        let var = self.get_main_var().unwrap_or_default();
        for (mono, coeff) in &self.terms {
            let deg = mono.0.get(&var).cloned().unwrap_or(0) as usize;
            if deg <= max_deg {
                vec[deg] = coeff.clone();
            }
        }
        vec
    }
}
*/
*/
