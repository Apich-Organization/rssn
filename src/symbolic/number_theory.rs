//! # Number Theory Module
//!
//! This module provides a collection of tools for symbolic number theory.
//! It includes solvers for various types of Diophantine equations (linear, Pell's, Pythagorean),
//! functions for primality testing, continued fraction expansion, and the Chinese Remainder Theorem.

use crate::symbolic::core::{Expr, Monomial, SparsePolynomial};
use crate::symbolic::simplify::{is_one, simplify};
use num_bigint::{BigInt, ToBigInt as OtherToBigInt};
use num_traits::{One, ToPrimitive, Zero};
use std::collections::{BTreeMap, HashMap};

// =====================================================================================
// region: Helpers
// =====================================================================================

trait ToBigInt {
    fn to_bigint(&self) -> Option<BigInt>;
}

impl ToBigInt for Expr {
    fn to_bigint(&self) -> Option<BigInt> {
        match self {
            Expr::BigInt(i) => Some(i.clone()),
            Expr::Constant(f) => f.to_bigint(),
            Expr::Rational(r) => r.to_integer().into(),
            _ => None,
        }
    }
}

/// Converts a symbolic expression into a `SparsePolynomial` representation.
///
/// This function recursively traverses the expression tree and collects terms,
/// grouping them by their monomial basis. It is a foundational step for many
/// polynomial-specific algorithms.
///
/// **Note:** This is a simplified implementation. It handles basic arithmetic
/// (`+`, `-`, `*`, `^`) but may not correctly expand all complex nested expressions.
/// For instance, `(x+1)*(x+2)` is handled, but more complex scenarios might be
/// treated as non-polynomial terms.
///
/// # Arguments
/// * `expr` - The symbolic expression to convert.
///
/// # Returns
/// A `SparsePolynomial` representing the input expression.
pub fn expr_to_sparse_poly(expr: &Expr) -> SparsePolynomial {
    let mut terms = BTreeMap::new();
    collect_poly_terms_recursive(expr, &mut terms, &Expr::BigInt(BigInt::one()));
    SparsePolynomial { terms }
}

/// Recursive helper to collect terms for the sparse polynomial conversion.
pub(crate) fn collect_poly_terms_recursive(
    expr: &Expr,
    terms: &mut BTreeMap<Monomial, Expr>,
    current_coeff: &Expr,
) {
    match &simplify(expr.clone()) {
        Expr::Add(a, b) => {
            collect_poly_terms_recursive(a, terms, current_coeff);
            collect_poly_terms_recursive(b, terms, current_coeff);
        }
        Expr::Sub(a, b) => {
            collect_poly_terms_recursive(a, terms, current_coeff);
            let neg_coeff = simplify(Expr::Neg(Box::new(current_coeff.clone())));
            collect_poly_terms_recursive(b, terms, &neg_coeff);
        }
        Expr::Mul(a, b) => {
            // This is a simplification. A full implementation would require expanding
            // products of sums, which is complex. We assume one factor is the coefficient.
            if let Some(val) = a.to_bigint() {
                let next_coeff = simplify(Expr::Mul(
                    Box::new(current_coeff.clone()),
                    Box::new(Expr::BigInt(val)),
                ));
                collect_poly_terms_recursive(b, terms, &next_coeff);
            } else if let Some(val) = b.to_bigint() {
                let next_coeff = simplify(Expr::Mul(
                    Box::new(current_coeff.clone()),
                    Box::new(Expr::BigInt(val)),
                ));
                collect_poly_terms_recursive(a, terms, &next_coeff);
            } else {
                // Fallback for complex multiplications, treat as non-polynomial term for now
                let mono = Monomial(BTreeMap::new()); // Constant term
                let entry = terms
                    .entry(mono)
                    .or_insert_with(|| Expr::BigInt(BigInt::zero()));
                *entry = simplify(Expr::Add(Box::new(entry.clone()), Box::new(expr.clone())));
            }
        }
        Expr::Power(base, exp) => {
            if let (Expr::Variable(v), Some(e)) = (&**base, exp.to_bigint()) {
                if let Some(e_u32) = e.to_u32() {
                    let mut mono_map = BTreeMap::new();
                    mono_map.insert(v.clone(), e_u32);
                    let mono = Monomial(mono_map);
                    let entry = terms
                        .entry(mono)
                        .or_insert_with(|| Expr::BigInt(BigInt::zero()));
                    *entry = simplify(Expr::Add(
                        Box::new(entry.clone()),
                        Box::new(current_coeff.clone()),
                    ));
                }
            }
        }
        Expr::Variable(v) => {
            let mut mono_map = BTreeMap::new();
            mono_map.insert(v.clone(), 1);
            let mono = Monomial(mono_map);
            let entry = terms
                .entry(mono)
                .or_insert_with(|| Expr::BigInt(BigInt::zero()));
            *entry = simplify(Expr::Add(
                Box::new(entry.clone()),
                Box::new(current_coeff.clone()),
            ));
        }
        Expr::Neg(a) => {
            let neg_coeff = simplify(Expr::Neg(Box::new(current_coeff.clone())));
            collect_poly_terms_recursive(a, terms, &neg_coeff);
        }
        // Default: treat as a constant term
        e => {
            let mono = Monomial(BTreeMap::new()); // Represents the constant term
            let entry = terms
                .entry(mono)
                .or_insert_with(|| Expr::BigInt(BigInt::zero()));
            let term = simplify(Expr::Mul(
                Box::new(current_coeff.clone()),
                Box::new(e.clone()),
            ));
            *entry = simplify(Expr::Add(Box::new(entry.clone()), Box::new(term)));
        }
    }
}

// =====================================================================================
// endregion: Helpers
// =====================================================================================

// =====================================================================================
// region: Core Solvers
// =====================================================================================

/// Solves the linear Diophantine equation `ax + by = c`.
pub(crate) fn solve_linear_diophantine(
    coeffs: &HashMap<String, Expr>,
    c: &Expr,
    vars: &[&str],
) -> Result<Vec<Expr>, String> {
    let a = coeffs.get(vars[0]).ok_or("Var not found")?.clone();
    let b = coeffs.get(vars[1]).ok_or("Var not found")?.clone();

    let (a_int, b_int, c_int) = (
        a.to_bigint().ok_or("Coefficient 'a' must be an integer.")?,
        b.to_bigint().ok_or("Coefficient 'b' must be an integer.")?,
        c.to_bigint().ok_or("Constant 'c' must be an integer.")?,
    );

    let (g, x0, y0) = extended_gcd_inner(a_int.clone(), b_int.clone());

    if &c_int % &g != BigInt::zero() {
        return Err("No integer solutions exist (gcd does not divide constant).".to_string());
    }

    let factor = &c_int / &g;
    let x0_sol = &x0 * &factor;
    let y0_sol = &y0 * &factor;

    let t = Expr::Variable("t".to_string());
    let b_div_g = &b_int / &g;
    let a_div_g = &a_int / &g;

    let x_sol_expr = Expr::Add(
        Box::new(Expr::BigInt(x0_sol)),
        Box::new(Expr::Mul(
            Box::new(Expr::BigInt(b_div_g)),
            Box::new(t.clone()),
        )),
    );

    let y_sol_expr = Expr::Sub(
        Box::new(Expr::BigInt(y0_sol)),
        Box::new(Expr::Mul(Box::new(Expr::BigInt(a_div_g)), Box::new(t))),
    );

    Ok(vec![x_sol_expr, y_sol_expr])
}

/// Solves Pell's equation x^2 - n*y^2 = 1.
pub(crate) fn solve_pell(n: &Expr) -> Result<(Expr, Expr), String> {
    let n_val = n
        .to_bigint()
        .ok_or("n in Pell's equation must be an integer.")?;
    let n_f64 = n_val.to_f64().ok_or("n is too large to process.")?;

    if n_f64.sqrt().fract() == 0.0 {
        return Err("n cannot be a perfect square.".to_string());
    }

    let (a0, period) = sqrt_continued_fraction(n).ok_or("Failed to compute continued fraction.")?;

    let m = period.len();
    let k = if m % 2 == 0 { m - 1 } else { 2 * m - 1 };

    let (h, p) = get_convergent(a0, &period, k);

    Ok((Expr::BigInt(h), Expr::BigInt(p)))
}

/// Generates the parametric solution for the Pythagorean equation x^2 + y^2 = z^2.
pub(crate) fn solve_pythagorean(
    coeffs: &HashMap<String, Expr>,
    vars: &[&str],
) -> Result<Vec<Expr>, String> {
    let x_var = vars[0];
    let y_var = vars[1];
    let z_var = vars[2];

    // Find which variable has the negative coefficient (the hypotenuse)
    let (c1, c2, c3) = (
        coeffs.get(x_var).unwrap(),
        coeffs.get(y_var).unwrap(),
        coeffs.get(z_var).unwrap(),
    );

    let (x, y, z) = if is_neg_one(c3) {
        (x_var, y_var, z_var)
    } else if is_neg_one(c2) {
        (x_var, z_var, y_var)
    } else if is_neg_one(c1) {
        (y_var, z_var, x_var)
    } else {
        return Err("Equation is not in a recognized Pythagorean form.".to_string());
    };

    let m = Expr::Variable("m".to_string());
    let n = Expr::Variable("n".to_string());
    let k = Expr::Variable("k".to_string());

    let m_sq = Expr::Power(Box::new(m.clone()), Box::new(Expr::BigInt(BigInt::from(2))));
    let n_sq = Expr::Power(Box::new(n.clone()), Box::new(Expr::BigInt(BigInt::from(2))));

    let x_sol = Expr::Mul(
        Box::new(k.clone()),
        Box::new(Expr::Sub(Box::new(m_sq.clone()), Box::new(n_sq.clone()))),
    );
    let y_sol = Expr::Mul(
        Box::new(k.clone()),
        Box::new(Expr::Mul(
            Box::new(Expr::BigInt(BigInt::from(2))),
            Box::new(Expr::Mul(Box::new(m.clone()), Box::new(n.clone()))),
        )),
    );
    let z_sol = Expr::Mul(
        Box::new(k),
        Box::new(Expr::Add(Box::new(m_sq), Box::new(n_sq))),
    );

    // Return solutions in the order of the original variables
    let mut solutions = vec![Expr::Constant(0.0); 3];
    solutions[vars.iter().position(|&v| v == x).unwrap()] = x_sol;
    solutions[vars.iter().position(|&v| v == y).unwrap()] = y_sol;
    solutions[vars.iter().position(|&v| v == z).unwrap()] = z_sol;

    Ok(solutions)
}

// =====================================================================================
// endregion: Core Solvers
// =====================================================================================

// =====================================================================================
// region: Public API
// =====================================================================================

/// Acts as a dispatcher to solve various types of Diophantine equations.
///
/// This function analyzes the structure of the input equation to determine its type
/// and then calls the appropriate specialized solver.
///
/// It currently supports:
/// - **Linear Diophantine Equations**: `ax + by = c` for two variables.
/// - **Pythagorean Equations**: `x^2 + y^2 = z^2` (and its variations).
/// - **Pell's Equations**: `x^2 - ny^2 = 1`.
///
/// The function first simplifies the equation `LHS = RHS` into the form `f(vars) = 0`
/// and converts it to a sparse polynomial to analyze its properties (degree, coefficients).
///
/// # Arguments
/// * `equation` - An `Expr::Eq` representing the Diophantine equation.
/// * `vars` - A slice of strings representing the variables to solve for.
///
/// # Returns
/// * `Ok(Vec<Expr>)` containing the parametric solutions for the variables if a solver is found.
/// * `Err(String)` if the equation type is not recognized or not supported.
pub fn solve_diophantine(equation: &Expr, vars: &[&str]) -> Result<Vec<Expr>, String> {
    let (lhs, rhs) = match equation {
        Expr::Eq(l, r) => (l, r),
        _ => return Err("Input must be an equation.".to_string()),
    };

    let simplified_expr = simplify(Expr::Sub(lhs.clone(), rhs.clone()));
    let poly = expr_to_sparse_poly(&simplified_expr);

    let mut var_coeffs = HashMap::new();
    let mut constant_term = Expr::BigInt(BigInt::zero());
    let mut degrees = HashMap::new();
    let mut is_linear = true;

    for (mono, coeff) in &poly.terms {
        if mono.0.is_empty() {
            constant_term = simplify(Expr::Add(Box::new(constant_term), Box::new(coeff.clone())));
            continue;
        }
        let total_degree = mono.0.values().sum::<u32>();
        if total_degree > 1 {
            is_linear = false;
        }

        if mono.0.len() == 1 {
            let (var, &deg) = mono.0.iter().next().unwrap();
            degrees.entry(var.clone()).or_insert(Vec::new()).push(deg);
            if deg == 1 {
                let entry = var_coeffs
                    .entry(var.clone())
                    .or_insert_with(|| Expr::BigInt(BigInt::zero()));
                *entry = simplify(Expr::Add(Box::new(entry.clone()), Box::new(coeff.clone())));
            }
        }
    }

    // --- Linear Solver: ax + by = c ---
    if is_linear && vars.len() == 2 {
        let c = simplify(Expr::Neg(Box::new(constant_term)));
        return solve_linear_diophantine(&var_coeffs, &c, vars);
    }

    // --- Pythagorean Solver: x^2 + y^2 = z^2 ---
    if vars.len() == 3
        && poly.terms.len() == 3
        && poly
            .terms
            .values()
            .all(|c| is_one(&simplify(Expr::Abs(Box::new(c.clone())))))
    {
        let mut neg_count = 0;
        for c in poly.terms.values() {
            if is_neg_one(c) {
                neg_count += 1;
            }
        }

        if neg_count == 1 {
            let mut coeffs_map = HashMap::new();
            for (mono, coeff) in poly.terms {
                if let Some((var, &2)) = mono.0.iter().next() {
                    coeffs_map.insert(var.clone(), coeff);
                }
            }
            return solve_pythagorean(&coeffs_map, vars);
        }
    }

    // --- Pell's Solver: x^2 - ny^2 = 1 ---
    if vars.len() == 2 && poly.terms.len() == 3 {
        // Simplified check for Pell's equation
        // A robust check would analyze the sparse polynomial terms in more detail
        if let Ok(solutions) = solve_pell_from_poly(&poly, vars) {
            return Ok(solutions);
        }
    }

    Err("Equation type not recognized or not supported.".to_string())
}

/// Attempts to solve Pell's equation given a `SparsePolynomial` representation.
///
/// This function inspects the terms of a sparse polynomial to see if it matches
/// the form `x^2 - ny^2 - 1 = 0`. It extracts the coefficient `n` and calls the
/// core Pell's equation solver.
///
/// # Arguments
/// * `poly` - The `SparsePolynomial` representation of the equation.
/// * `vars` - A slice containing the two variable names (e.g., `["x", "y"]`).
///
/// # Returns
/// * `Ok(Vec<Expr>)` with the fundamental solution `(x, y)` if successful.
/// * `Err(String)` if the polynomial does not match the recognized form of Pell's equation.
pub fn solve_pell_from_poly(poly: &SparsePolynomial, vars: &[&str]) -> Result<Vec<Expr>, String> {
    let mut n_coeff: Option<Expr> = None;
    let mut const_term: Option<Expr> = None;

    for (mono, coeff) in &poly.terms {
        if mono.0.len() == 1 && mono.0.get(vars[0]) == Some(&2) && is_one(coeff) {
            // x^2 term
        } else if mono.0.len() == 1 && mono.0.get(vars[1]) == Some(&2) {
            // y^2 term
            n_coeff = Some(simplify(Expr::Neg(Box::new(coeff.clone()))));
        } else if mono.0.is_empty() {
            const_term = Some(coeff.clone());
        }
    }

    if let (Some(n), Some(k)) = (n_coeff, const_term) {
        if is_neg_one(&k) {
            // Equation was x^2 - ny^2 - 1 = 0
            let (x_sol, y_sol) = solve_pell(&n)?;
            return Ok(vec![x_sol, y_sol]);
        }
    }
    Err("Not a recognized Pell's equation form".to_string())
}

/// Checks if an expression is numerically equal to -1.
///
/// Handles `Expr::Constant` and `Expr::BigInt` variants.
pub fn is_neg_one(expr: &Expr) -> bool {
    matches!(expr, Expr::Constant(val) if *val == -1.0)
        || matches!(expr, Expr::BigInt(val) if *val == BigInt::from(-1))
}

/// Checks if an expression is numerically equal to 2.
///
/// Handles `Expr::Constant` and `Expr::BigInt` variants.
pub fn is_two(expr: &Expr) -> bool {
    matches!(expr, Expr::Constant(val) if *val == 2.0)
        || matches!(expr, Expr::BigInt(val) if *val == BigInt::from(2))
}

/// The internal, recursive implementation of the Extended Euclidean Algorithm for `BigInt`.
///
/// This function computes the greatest common divisor (GCD) of `a` and `b`
/// and also finds the Bézout coefficients `x` and `y` such that `a*x + b*y = gcd(a, b)`.
///
/// # Arguments
/// * `a` - The first integer.
/// * `b` - The second integer.
///
/// # Returns
/// A tuple `(g, x, y)` where `g` is the GCD, and `x`, `y` are the coefficients.
pub fn extended_gcd_inner(a: BigInt, b: BigInt) -> (BigInt, BigInt, BigInt) {
    if b.is_zero() {
        return (a, BigInt::one(), BigInt::zero());
    }
    let (g, x, y) = extended_gcd_inner(b.clone(), &a % &b);
    (g, y.clone(), x - (a / b) * y)
}

/// Solves a system of congruences using the Chinese Remainder Theorem.
///
/// Given a set of congruences of the form `x ≡ a_i (mod n_i)`,
/// this function finds a single solution `x` that satisfies all of them simultaneously,
/// provided that the moduli `n_i` are pairwise coprime.
///
/// # Arguments
/// * `congruences` - A slice of tuples, where each tuple `(a, n)` represents
///   the congruence `x ≡ a (mod n)`.
///
/// # Returns
/// * `Some(Expr)` containing the smallest non-negative solution `x` if the moduli are pairwise coprime.
/// * `None` if the moduli are not pairwise coprime, in which case a unique solution is not guaranteed by this method.
pub fn chinese_remainder(congruences: &[(Expr, Expr)]) -> Option<Expr> {
    let mut n_total = Expr::BigInt(BigInt::one());
    for (_, n) in congruences {
        n_total = simplify(Expr::Mul(Box::new(n_total), Box::new(n.clone())));
    }

    let mut result = Expr::BigInt(BigInt::zero());
    for (a, n) in congruences {
        let n_i = simplify(Expr::Div(Box::new(n_total.clone()), Box::new(n.clone())));
        let (g, m_i, _) = extended_gcd(&n_i, n);
        if !is_one(&g) {
            return None;
        }
        result = simplify(Expr::Add(
            Box::new(result),
            Box::new(simplify(Expr::Mul(
                Box::new(a.clone()),
                Box::new(simplify(Expr::Mul(Box::new(n_i), Box::new(m_i)))),
            ))),
        ));
    }

    Some(simplify(Expr::Mod(Box::new(result), Box::new(n_total))))
}

/// Performs a primality test on a given expression.
///
/// If the expression can be evaluated to a `BigInt`, this function uses a deterministic
/// trial division method to check for primality. It is efficient for reasonably sized numbers.
///
/// If the expression is symbolic (e.g., contains variables), it returns a symbolic
/// representation `Expr::IsPrime(Box(expr))`, which represents the primality test
/// as an unevaluated function.
///
/// # Arguments
/// * `n` - The expression to test for primality.
///
/// # Returns
/// * `Expr::Boolean(true)` if `n` is a prime number.
/// * `Expr::Boolean(false)` if `n` is a composite number.
/// * `Expr::IsPrime(...)` if `n` is a symbolic expression.
pub fn is_prime(n: &Expr) -> Expr {
    if let Some(n_bigint) = n.to_bigint() {
        if n_bigint <= BigInt::one() {
            return Expr::Boolean(false);
        }
        if n_bigint <= BigInt::from(3) {
            return Expr::Boolean(true);
        }
        if &n_bigint % 2 == BigInt::zero() || &n_bigint % 3 == BigInt::zero() {
            return Expr::Boolean(false);
        }
        let mut i = BigInt::from(5);
        while &i * &i <= n_bigint {
            if &n_bigint % &i == BigInt::zero() || &n_bigint % (&i + 2) == BigInt::zero() {
                return Expr::Boolean(false);
            }
            i += 6;
        }
        Expr::Boolean(true)
    } else {
        Expr::IsPrime(Box::new(n.clone()))
    }
}

/// Computes the continued fraction expansion of the square root of an integer.
///
/// For a non-square integer `n`, the continued fraction of `sqrt(n)` is periodic:
/// `sqrt(n) = [a0; a1, a2, ..., an, 2*a0, a1, ...]`
/// This function finds the initial term `a0` and the periodic part `[a1, ..., an]`.
/// This is a key step in solving Pell's equation.
///
/// # Arguments
/// * `n_expr` - An expression representing the integer `n`.
///
/// # Returns
/// * `Some((a0, period))` where `a0` is the integer part of the square root and
///   `period` is a `Vec<i64>` containing the repeating block of the continued fraction.
/// * `None` if `n` is a perfect square or cannot be converted to an integer.
pub fn sqrt_continued_fraction(n_expr: &Expr) -> Option<(i64, Vec<i64>)> {
    let n = n_expr.to_bigint()?.to_f64()? as i64;
    let sqrt_n_floor = (n as f64).sqrt().floor() as i64;

    if sqrt_n_floor * sqrt_n_floor == n {
        return None;
    }

    let mut m = 0;
    let mut d = 1;
    let mut a = sqrt_n_floor;

    let mut history = HashMap::new();
    let mut periodic_part = Vec::new();

    loop {
        let state = (m, d, a);
        if history.contains_key(&state) {
            break;
        }
        history.insert(state, periodic_part.len());

        m = d * a - m;
        d = (n - m * m) / d;
        if d == 0 {
            return None;
        }
        a = (sqrt_n_floor + m) / d;
        periodic_part.push(a);
    }

    Some((sqrt_n_floor, periodic_part))
}

/// Calculates the k-th convergent of a continued fraction.
///
/// Given the initial term `a0` and the periodic part of a continued fraction,
/// this function computes the numerator `h_k` and denominator `p_k` of the k-th
/// convergent fraction `h_k / p_k`.
///
/// The convergents are calculated using the recurrence relations:
/// `h_i = a_i * h_{i-1} + h_{i-2}`
/// `p_i = a_i * p_{i-1} + p_{i-2}`
///
/// # Arguments
/// * `a0` - The initial, non-periodic term of the continued fraction.
/// * `period` - A slice representing the periodic part of the fraction.
/// * `k` - The index of the convergent to compute.
///
/// # Returns
/// A tuple `(h, p)` containing the numerator and denominator of the k-th convergent as `BigInt`s.
pub fn get_convergent(a0: i64, period: &[i64], k: usize) -> (BigInt, BigInt) {
    let mut h_minus_2 = BigInt::from(0);
    let mut h_minus_1 = BigInt::from(1);
    let mut p_minus_2 = BigInt::from(1);
    let mut p_minus_1 = BigInt::from(0);

    for i in 0..=k {
        let a_i = if i == 0 {
            BigInt::from(a0)
        } else {
            BigInt::from(period[(i - 1) % period.len()])
        };
        let h_i = &a_i * &h_minus_1 + &h_minus_2;
        let p_i = &a_i * &p_minus_1 + &p_minus_2;
        h_minus_2 = h_minus_1;
        h_minus_1 = h_i;
        p_minus_2 = p_minus_1;
        p_minus_1 = p_i;
    }
    (h_minus_1, p_minus_1)
}

/// Public wrapper for the Extended Euclidean Algorithm.
///
/// This function attempts to convert the input expressions `a` and `b` to `BigInt`s.
/// - If successful, it calls the concrete `extended_gcd_inner` implementation and returns the result as `Expr::BigInt`.
/// - If the inputs are symbolic, it returns a tuple of symbolic variables `(g, x, y)`
///   representing the result of the computation.
///
/// # Arguments
/// * `a` - The first expression.
/// * `b` - The second expression.
///
/// # Returns
/// A tuple `(g, x, y)` of expressions representing the GCD and Bézout coefficients.
pub fn extended_gcd(a: &Expr, b: &Expr) -> (Expr, Expr, Expr) {
    if let (Some(a_int), Some(b_int)) = (a.to_bigint(), b.to_bigint()) {
        let (g, x, y) = extended_gcd_inner(a_int, b_int);
        return (Expr::BigInt(g), Expr::BigInt(x), Expr::BigInt(y));
    }
    (
        Expr::Variable("g".to_string()),
        Expr::Variable("x".to_string()),
        Expr::Variable("y".to_string()),
    )
}

// =====================================================================================
// endregion: Public API
// =====================================================================================
