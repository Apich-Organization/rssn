//! # Symbolic Polynomial Manipulation
//!
//! This module provides comprehensive tools for symbolic manipulation of polynomials,
//! supporting both univariate and multivariate polynomials in various representations.
//!
//! ## Overview
//!
//! The module offers two main polynomial representations:
//!
//! 1. **Expression-based**: Polynomials as symbolic [`Expr`] trees
//! 2. **Sparse representation**: Polynomials as [`SparsePolynomial`] with explicit monomial-coefficient pairs
//!
//! ## Key Features
//!
//! ### Basic Operations
//! - **Addition/Subtraction**: [`add_poly`], [`subtract_poly`](crate::symbolic::grobner::subtract_poly)
//! - **Multiplication**: [`mul_poly`]
//! - **Differentiation**: [`differentiate_poly`]
//! - **Division**: [`polynomial_long_division`], [`polynomial_long_division_coeffs`]
//!
//! ### Analysis
//! - **Degree computation**: [`polynomial_degree`]
//! - **Leading coefficient**: [`leading_coefficient`]
//! - **Polynomial detection**: [`is_polynomial`]
//! - **Variable detection**: [`contains_var`]
//!
//! ### Representation Conversion
//! - **Expression to sparse**: [`expr_to_sparse_poly`]
//! - **Sparse to expression**: [`sparse_poly_to_expr`]
//! - **Coefficient vectors**: [`to_polynomial_coeffs_vec`], [`from_coeffs_to_expr`]
//!
//! ### Advanced Operations
//! - **GCD computation**: [`gcd`]
//! - **Scalar multiplication**: [`poly_mul_scalar_expr`]
//! - **Evaluation**: [`SparsePolynomial::eval`]
//!
//! ## Representations
//!
//! ### Expression-Based Polynomials
//!
//! Polynomials can be represented as standard [`Expr`] trees:
//!
//! ```rust
//! 
//! use rssn::symbolic::core::Expr;
//!
//! // x^2 + 2x + 1
//! let poly = Expr::new_add(
//!     Expr::new_add(
//!         Expr::new_pow(
//!             Expr::new_variable("x"),
//!             Expr::new_constant(2.0),
//!         ),
//!         Expr::new_mul(
//!             Expr::new_constant(2.0),
//!             Expr::new_variable("x"),
//!         ),
//!     ),
//!     Expr::new_constant(1.0),
//! );
//! ```
//!
//! ### Sparse Polynomial Representation
//!
//! For multivariate polynomials, the sparse representation is more efficient:
//!
//! ```rust
//! 
//! use rssn::symbolic::core::{
//!     Expr,
//!     Monomial,
//!     SparsePolynomial,
//! };
//! use rssn::symbolic::polynomial::expr_to_sparse_poly;
//!
//! let expr = Expr::new_add(
//!     Expr::new_mul(
//!         Expr::new_variable("x"),
//!         Expr::new_variable("y"),
//!     ),
//!     Expr::new_constant(1.0),
//! );
//!
//! let sparse = expr_to_sparse_poly(&expr, &["x", "y"]);
//! // Represents: x*y + 1
//! ```
//!
//! ## Examples
//!
//! ### Polynomial Long Division
//!
//! ```rust
//! 
//! use rssn::symbolic::core::Expr;
//! use rssn::symbolic::polynomial::polynomial_long_division;
//!
//! // Divide x^2 + 3x + 2 by x + 1
//! let dividend = Expr::new_add(
//!     Expr::new_add(
//!         Expr::new_pow(
//!             Expr::new_variable("x"),
//!             Expr::new_constant(2.0),
//!         ),
//!         Expr::new_mul(
//!             Expr::new_constant(3.0),
//!             Expr::new_variable("x"),
//!         ),
//!     ),
//!     Expr::new_constant(2.0),
//! );
//!
//! let divisor = Expr::new_add(
//!     Expr::new_variable("x"),
//!     Expr::new_constant(1.0),
//! );
//!
//! let (quotient, remainder) = polynomial_long_division(
//!     &dividend, &divisor, "x",
//! );
//! // quotient = x + 2, remainder = 0
//! ```
//!
//! ### Differentiation
//!
//! ```rust
//! 
//! use rssn::symbolic::core::{
//!     Expr,
//!     SparsePolynomial,
//! };
//! use rssn::symbolic::polynomial::{
//!     differentiate_poly,
//!     expr_to_sparse_poly,
//!     sparse_poly_to_expr,
//! };
//!
//! let expr = Expr::new_pow(
//!     Expr::new_variable("x"),
//!     Expr::new_constant(3.0),
//! );
//!
//! let poly = expr_to_sparse_poly(
//!     &expr,
//!     &["x"],
//! );
//!
//! let derivative =
//!     differentiate_poly(&poly, "x");
//! // Result: 3x^2
//! ```
//!
//! ### GCD Computation
//!
//! ```rust
//! 
//! use rssn::symbolic::core::SparsePolynomial;
//! use rssn::symbolic::polynomial::gcd;
//!
//! // Find GCD of two polynomials
//! // let gcd_poly = gcd(poly1, poly2, "x");
//! ```
//!
//! ## Performance Considerations
//!
//! - **Sparse representation**: More efficient for multivariate polynomials with few terms
//! - **Expression-based**: Better for symbolic manipulation and simplification
//! - **Coefficient vectors**: Fastest for dense univariate polynomials
//!
//! ## See Also
//!
//! - [`grobner`](crate::symbolic::grobner) - Gröbner basis computation
//! - [`poly_factorization`](crate::symbolic::poly_factorization) - Polynomial factorization
//! - [`real_roots`](crate::symbolic::real_roots) - Finding real roots of polynomials
//! - [`core`](crate::symbolic::core) - Core expression types

use crate::symbolic::core::{
    Expr,
    Monomial,
    SparsePolynomial,
};
use crate::symbolic::grobner::subtract_poly;
use crate::symbolic::real_roots::eval_expr;
use crate::symbolic::simplify::as_f64;
use crate::symbolic::simplify::is_zero;
use crate::symbolic::simplify_dag::simplify;
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{
    One,
    ToPrimitive,
    Zero,
};
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
#[must_use]

pub fn add_poly(
    p1: &SparsePolynomial,
    p2: &SparsePolynomial,
) -> SparsePolynomial {

    let mut result_terms =
        p1.terms.clone();

    for (monomial, coeff2) in &p2.terms
    {

        let coeff1 = result_terms
            .entry(monomial.clone())
            .or_insert(Expr::Constant(
                0.0,
            ));

        *coeff1 = Expr::new_add(
            coeff1.clone(),
            coeff2.clone(),
        );
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
#[must_use]

pub fn mul_poly(
    p1: &SparsePolynomial,
    p2: &SparsePolynomial,
) -> SparsePolynomial {

    let mut result_terms: BTreeMap<
        Monomial,
        Expr,
    > = BTreeMap::new();

    for (m1, c1) in &p1.terms {

        for (m2, c2) in &p2.terms {

            let new_coeff =
                Expr::new_mul(
                    c1.clone(),
                    c2.clone(),
                );

            let mut new_mono_map =
                m1.0.clone();

            for (var, exp2) in &m2.0 {

                let exp1 = new_mono_map
                    .entry(var.clone())
                    .or_insert(0);

                *exp1 += exp2;
            }

            let new_mono =
                Monomial(new_mono_map);

            let existing_coeff =
                result_terms
                    .entry(new_mono)
                    .or_insert(
                        Expr::Constant(
                            0.0,
                        ),
                    );

            *existing_coeff =
                Expr::new_add(
                    existing_coeff
                        .clone(),
                    new_coeff,
                );
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
#[must_use]

pub fn differentiate_poly(
    p: &SparsePolynomial,
    var: &str,
) -> SparsePolynomial {

    let mut result_terms: BTreeMap<
        Monomial,
        Expr,
    > = BTreeMap::new();

    for (monomial, coeff) in &p.terms {

        if let Some(&exp) =
            monomial.0.get(var)
        {

            if exp > 0 {

                let new_coeff =
                    Expr::new_mul(
                        coeff.clone(),
                        Expr::Constant(
                            f64::from(
                                exp,
                            ),
                        ),
                    );

                let mut new_mono_map =
                    monomial.0.clone();

                if exp == 1 {

                    new_mono_map
                        .remove(var);
                } else if let Some(e) =
                    new_mono_map
                        .get_mut(var)
                {

                    *e -= 1;
                }

                let new_mono = Monomial(
                    new_mono_map,
                );

                result_terms.insert(
                    new_mono, new_coeff,
                );
            }
        }
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
#[must_use]

pub fn contains_var(
    expr: &Expr,
    var: &str,
) -> bool {

    let mut found = false;

    expr.pre_order_walk(&mut |e| {
        if let Expr::Variable(name) = e
        {

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
#[must_use]

pub fn is_polynomial(
    expr: &Expr,
    var: &str,
) -> bool {

    match expr {
        | Expr::Dag(node) => {
            is_polynomial(
                &node
                    .to_expr()
                    .expect(
                        "Is Polynomial",
                    ),
                var,
            )
        },
        | Expr::Constant(_)
        | Expr::BigInt(_)
        | Expr::Rational(_) => true,
        | Expr::Variable(_) => true,
        | Expr::Add(a, b)
        | Expr::Sub(a, b)
        | Expr::Mul(a, b) => {
            is_polynomial(a, var)
                && is_polynomial(b, var)
        },
        | Expr::Div(a, b) => {
            is_polynomial(a, var)
                && !contains_var(b, var)
        },
        | Expr::Power(base, exp) => {

            // Check if exponent is a non-negative integer
            let is_valid_exp = match &**exp {
                | Expr::BigInt(n) => n >= &BigInt::zero(),
                | Expr::Constant(c) => *c >= 0.0 && c.fract() == 0.0,
                | Expr::Rational(r) => r >= &BigRational::zero() && r.is_integer(),
                | _ => false,
            };

            if is_valid_exp {

                is_polynomial(base, var)
            } else {

                // Negative or non-integer exponent
                !contains_var(base, var)
            }
        },
        | Expr::Neg(a) => {
            is_polynomial(a, var)
        },
        // N-ary list variants
        | Expr::AddList(terms)
        | Expr::MulList(terms) => {
            terms
                .iter()
                .all(|t| {
                    is_polynomial(
                        t, var,
                    )
                })
        },
        // Generic list variants - check if they don't contain the variable
        | Expr::UnaryList(_, a) => {
            is_polynomial(a, var)
        },
        | Expr::BinaryList(_, a, b) => {
            is_polynomial(a, var)
                && is_polynomial(b, var)
        },
        | Expr::NaryList(_, args) => {
            args.iter()
                .all(|arg| {
                    is_polynomial(
                        arg, var,
                    )
                })
        },
        | Expr::Sin(_)
        | Expr::Cos(_)
        | Expr::Tan(_)
        | Expr::Log(_)
        | Expr::Exp(_)
        | Expr::Sec(_)
        | Expr::Csc(_)
        | Expr::Cot(_) => {
            !contains_var(expr, var)
        },
        | _ => false,
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
#[must_use]

pub fn polynomial_degree(
    expr: &Expr,
    var: &str,
) -> i64 {

    let s_expr =
        simplify(&expr.clone());

    // Convert to AST if it's a DAG to properly match on variants
    let s_expr = if s_expr.is_dag() {

        s_expr
            .to_ast()
            .unwrap_or(s_expr)
    } else {

        s_expr
    };

    match s_expr {
        | Expr::Add(a, b)
        | Expr::Sub(a, b) => {
            std::cmp::max(
                polynomial_degree(
                    &a, var,
                ),
                polynomial_degree(
                    &b, var,
                ),
            )
        },
        | Expr::Mul(a, b) => {
            polynomial_degree(&a, var)
                + polynomial_degree(
                    &b, var,
                )
        },
        | Expr::Div(a, b) => {
            polynomial_degree(&a, var)
                - polynomial_degree(
                    &b, var,
                )
        },
        | Expr::Power(
            ref base,
            ref exp,
        ) => {

            // Check if base is the variable and exponent is a non-negative integer
            if let Expr::Variable(v) =
                base.as_ref()
            {

                if v == var {

                    // Extract the exponent value
                    return match exp.as_ref() {
                        | Expr::BigInt(n) => {
                            n.to_i64()
                                .unwrap_or(0)
                        },
                        | Expr::Constant(c) if c.fract() == 0.0 && *c >= 0.0 => *c as i64,
                        | Expr::Rational(r) if r.is_integer() && r >= &BigRational::zero() => {
                            r.to_integer()
                                .to_i64()
                                .unwrap_or(0)
                        },
                        | _ => 0,
                    };
                }
            }

            if contains_var(
                &s_expr, var,
            ) {

                -1
            } else {

                0
            }
        },
        // N-ary list variants
        | Expr::AddList(terms) => {

            // For AddList, degree is max of all terms
            terms
                .iter()
                .map(|t| {
                    polynomial_degree(
                        t, var,
                    )
                })
                .max()
                .unwrap_or(0)
        },
        | Expr::MulList(factors) => {

            // For MulList, degree is sum of all factors
            factors
                .iter()
                .map(|f| {
                    polynomial_degree(
                        f, var,
                    )
                })
                .sum()
        },
        | Expr::Variable(name)
            if name == var =>
        {
            1
        },
        | _ => 0,
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
#[must_use]

pub fn leading_coefficient(
    expr: &Expr,
    var: &str,
) -> Expr {

    let s_expr =
        simplify(&expr.clone());

    // Convert to AST if it's a DAG to properly match on variants
    let s_expr = if s_expr.is_dag() {

        s_expr
            .to_ast()
            .unwrap_or(s_expr)
    } else {

        s_expr
    };

    match s_expr {
        | Expr::Add(a, b) => {

            let deg_a =
                polynomial_degree(
                    &a, var,
                );

            let deg_b =
                polynomial_degree(
                    &b, var,
                );

            if deg_a > deg_b {

                leading_coefficient(
                    &a, var,
                )
            } else if deg_b > deg_a {

                leading_coefficient(
                    &b, var,
                )
            } else {

                simplify(&Expr::new_add(
                    leading_coefficient(&a, var),
                    leading_coefficient(&b, var),
                ))
            }
        },
        | Expr::Sub(a, b) => {

            let deg_a =
                polynomial_degree(
                    &a, var,
                );

            let deg_b =
                polynomial_degree(
                    &b, var,
                );

            if deg_a > deg_b {

                leading_coefficient(
                    &a, var,
                )
            } else if deg_b > deg_a {

                simplify(&Expr::new_neg(
                    leading_coefficient(&b, var),
                ))
            } else {

                simplify(&Expr::new_sub(
                    leading_coefficient(&a, var),
                    leading_coefficient(&b, var),
                ))
            }
        },
        | Expr::Mul(a, b) => {
            simplify(&Expr::new_mul(
                leading_coefficient(
                    &a, var,
                ),
                leading_coefficient(
                    &b, var,
                ),
            ))
        },
        | Expr::Div(a, b) => {
            simplify(&Expr::new_div(
                leading_coefficient(
                    &a, var,
                ),
                leading_coefficient(
                    &b, var,
                ),
            ))
        },
        | Expr::Power(base, exp) => {

            if let (
                Expr::Variable(v),
                Expr::BigInt(_),
            ) = (&*base, &*exp)
            {

                if v == var {

                    return Expr::BigInt(BigInt::one());
                }
            }

            simplify(&Expr::new_pow(
                leading_coefficient(
                    &base, var,
                ),
                exp,
            ))
        },
        | Expr::Variable(name)
            if name == var =>
        {
            Expr::BigInt(BigInt::one())
        },
        | _ => s_expr,
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
#[must_use]

pub fn polynomial_long_division(
    n: &Expr,
    d: &Expr,
    var: &str,
) -> (Expr, Expr) {

    pub(crate) fn is_zero_local(
        expr: &Expr
    ) -> bool {

        match expr {
            | Expr::Dag(node) => {
                is_zero_local(
                    &node
                        .to_expr()
                        .expect(
                            "Is Zero",
                        ),
                )
            },
            | Expr::Constant(c) => {
                *c == 0.0
            },
            | Expr::BigInt(i) => {
                i.is_zero()
            },
            | Expr::Rational(r) => {
                r.is_zero()
            },
            | _ => false,
        }
    }

    let mut q =
        Expr::BigInt(BigInt::zero());

    let mut r = n.clone();

    let d_deg =
        polynomial_degree(d, var);

    if d_deg < 0 {

        return (
            Expr::BigInt(BigInt::zero()),
            r,
        );
    }

    let mut r_deg =
        polynomial_degree(&r, var);

    let mut iterations = 0;

    let mut total_iterations = 0;

    const MAX_TOTAL_ITERATIONS: usize =
        100;

    while r_deg >= d_deg
        && !is_zero_local(&r)
        && total_iterations
            < MAX_TOTAL_ITERATIONS
    {

        let lead_r =
            leading_coefficient(
                &r, var,
            );

        let lead_d =
            leading_coefficient(d, var);

        let t_deg = r_deg - d_deg;

        let t_coeff =
            simplify(&Expr::new_div(
                lead_r, lead_d,
            ));

        let t = if t_deg == 0 {

            t_coeff
        } else {

            simplify(&Expr::new_mul(
                t_coeff,
                Expr::new_pow(
                    Expr::Variable(
                        var.to_string(),
                    ),
                    Expr::BigInt(
                        BigInt::from(
                            t_deg,
                        ),
                    ),
                ),
            ))
        };

        q = simplify(&Expr::new_add(
            q.clone(),
            t.clone(),
        ));

        let t_times_d =
            simplify(&Expr::new_mul(
                t,
                d.clone(),
            ));

        r = simplify(&Expr::new_sub(
            r, t_times_d,
        ));

        let new_r_deg =
            polynomial_degree(&r, var);

        if new_r_deg >= r_deg {

            iterations += 1;

            if iterations > 10 {

                break;
            }
        } else {

            iterations = 0;
        }

        r_deg = new_r_deg;

        total_iterations += 1;
    }

    (q, r)
}

/// Recursively collects coefficients of a polynomial expression into a map of degree -> coefficient.

pub(crate) fn collect_coeffs_recursive(
    expr: &Expr,
    var: &str,
) -> BTreeMap<u32, Expr> {

    let simplified =
        simplify(&expr.clone());

    // Convert to AST if it's a DAG to properly match on variants
    let s_expr = if simplified.is_dag()
    {

        simplified
            .to_ast()
            .unwrap_or(simplified)
    } else {

        simplified
    };

    match &s_expr {
        | Expr::Add(a, b) => {

            let mut map_a = collect_coeffs_recursive(a, var);

            let map_b = collect_coeffs_recursive(b, var);

            for (deg, coeff_b) in map_b
            {

                let coeff_a = map_a
                    .entry(deg)
                    .or_insert_with(|| Expr::BigInt(BigInt::zero()));

                *coeff_a = simplify(
                    &Expr::new_add(
                        coeff_a.clone(),
                        coeff_b,
                    ),
                );
            }

            map_a
        },
        | Expr::Sub(a, b) => {

            let mut map_a = collect_coeffs_recursive(a, var);

            let map_b = collect_coeffs_recursive(b, var);

            for (deg, coeff_b) in map_b
            {

                let coeff_a = map_a
                    .entry(deg)
                    .or_insert_with(|| Expr::BigInt(BigInt::zero()));

                *coeff_a = simplify(
                    &Expr::new_sub(
                        coeff_a.clone(),
                        coeff_b,
                    ),
                );
            }

            map_a
        },
        | Expr::Mul(a, b) => {

            let map_a = collect_coeffs_recursive(a, var);

            let map_b = collect_coeffs_recursive(b, var);

            let mut result_map =
                BTreeMap::new();

            for (deg_a, coeff_a) in
                &map_a
            {

                for (deg_b, coeff_b) in
                    &map_b
                {

                    let new_deg =
                        deg_a + deg_b;

                    let new_coeff_term = simplify(&Expr::new_mul(
                        coeff_a.clone(),
                        coeff_b.clone(),
                    ));

                    let entry = result_map
                        .entry(new_deg)
                        .or_insert_with(|| Expr::BigInt(BigInt::zero()));

                    *entry = simplify(&Expr::new_add(
                        entry.clone(),
                        new_coeff_term,
                    ));
                }
            }

            result_map
        },
        | Expr::Power(base, exp) => {

            if let (
                Expr::Variable(v),
                Expr::BigInt(n),
            ) = (
                base.as_ref(),
                exp.as_ref(),
            ) {

                if v == var {

                    let mut map =
                        BTreeMap::new();

                    map.insert(
                        n.to_u32()
                            .unwrap_or(
                                0,
                            ),
                        Expr::BigInt(
                            BigInt::one(
                            ),
                        ),
                    );

                    return map;
                }
            }

            if !contains_var(base, var)
            {

                let mut map =
                    BTreeMap::new();

                map.insert(
                    0,
                    expr.clone(),
                );

                return map;
            }

            BTreeMap::new()
        },
        | Expr::Variable(v)
            if v == var =>
        {

            let mut map =
                BTreeMap::new();

            map.insert(
                1,
                Expr::BigInt(
                    BigInt::one(),
                ),
            );

            map
        },
        | Expr::Neg(a) => {

            let map_a = collect_coeffs_recursive(a, var);

            let mut result_map =
                BTreeMap::new();

            for (deg, coeff) in map_a {

                result_map.insert(
                    deg,
                    simplify(
                        &Expr::new_neg(
                            coeff,
                        ),
                    ),
                );
            }

            result_map
        },
        | e if !contains_var(e, var) =>
        {

            let mut map =
                BTreeMap::new();

            map.insert(0, e.clone());

            map
        },
        | _ => BTreeMap::new(),
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
#[must_use]

pub fn to_polynomial_coeffs_vec(
    expr: &Expr,
    var: &str,
) -> Vec<Expr> {

    let map = collect_coeffs_recursive(
        expr, var,
    );

    if map.is_empty() {

        if !contains_var(expr, var) {

            return vec![expr.clone()];
        }

        return vec![];
    }

    let max_deg = map
        .keys()
        .max()
        .copied()
        .unwrap_or(0);

    let mut result = vec![
        Expr::BigInt(
            BigInt::zero()
        );
        max_deg
            as usize
            + 1
    ];

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
#[must_use]

pub fn from_coeffs_to_expr(
    coeffs: &[Expr],
    var: &str,
) -> Expr {

    let mut expr =
        Expr::BigInt(BigInt::zero());

    for (i, coeff) in coeffs
        .iter()
        .enumerate()
    {

        if !is_zero(&simplify(
            &coeff.clone(),
        )) {

            let power = if i == 0 {

                Expr::BigInt(
                    BigInt::one(),
                )
            } else {

                Expr::new_pow(
                    Expr::Variable(
                        var.to_string(),
                    ),
                    Expr::BigInt(
                        BigInt::from(i),
                    ),
                )
            };

            let term = if i == 0 {

                coeff.clone()
            } else if let Expr::BigInt(
                b,
            ) = coeff
            {

                if b.is_one() {

                    power
                } else {

                    Expr::new_mul(
                        coeff.clone(),
                        power,
                    )
                }
            } else {

                Expr::new_mul(
                    coeff.clone(),
                    power,
                )
            };

            expr = simplify(
                &Expr::new_add(
                    expr, term,
                ),
            );
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

pub fn polynomial_long_division_coeffs(
    n: &Expr,
    d: &Expr,
    var: &str,
) -> Result<(Expr, Expr), String> {

    let mut num_coeffs =
        to_polynomial_coeffs_vec(
            n, var,
        );

    let mut den_coeffs =
        to_polynomial_coeffs_vec(
            d, var,
        );

    while den_coeffs
        .last()
        .is_some_and(|c| {

            is_zero(&simplify(
                &c.clone(),
            ))
        })
    {

        den_coeffs.pop();
    }

    if den_coeffs.is_empty() {

        return Err("Polynomial \
                    division by zero"
            .to_string());
    }

    let den_deg = den_coeffs.len() - 1;

    let mut num_deg =
        num_coeffs.len() - 1;

    if num_deg < den_deg {

        return Ok((
            Expr::BigInt(BigInt::zero()),
            n.clone(),
        ));
    }

    let lead_den =
        match den_coeffs.last() {
            | Some(c) => c.clone(),
            | None => unreachable!(),
        };

    let mut quot_coeffs =
        vec![
            Expr::BigInt(BigInt::zero());
            num_deg - den_deg + 1
        ];

    while num_deg >= den_deg {

        let lead_num =
            num_coeffs[num_deg].clone();

        let coeff =
            simplify(&Expr::new_div(
                lead_num,
                lead_den.clone(),
            ));

        let deg_diff =
            num_deg - den_deg;

        if deg_diff < quot_coeffs.len()
        {

            quot_coeffs[deg_diff] =
                coeff.clone();
        }

        for (i, _item) in den_coeffs
            .iter()
            .enumerate()
            .take(den_deg + 1)
        {

            if let Some(num_coeff) =
                num_coeffs.get_mut(
                    deg_diff + i,
                )
            {

                let term_to_sub =
                    simplify(
                        &Expr::new_mul(
                            coeff
                                .clone(
                                ),
                            den_coeffs
                                [i]
                                .clone(
                                ),
                        ),
                    );

                *num_coeff = simplify(
                    &Expr::new_sub(
                        num_coeff
                            .clone(),
                        term_to_sub,
                    ),
                );
            }
        }

        while num_coeffs
            .last()
            .is_some_and(|c| {

                is_zero(&simplify(
                    &c.clone(),
                ))
            })
        {

            num_coeffs.pop();
        }

        if num_coeffs.is_empty() {

            num_deg = 0;

            let _help = num_deg;

            break;
        }

        num_deg = num_coeffs.len() - 1;
    }

    let quotient = from_coeffs_to_expr(
        &quot_coeffs,
        var,
    );

    let remainder = from_coeffs_to_expr(
        &num_coeffs,
        var,
    );

    Ok((quotient, remainder))
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
#[must_use]

pub fn expr_to_sparse_poly(
    expr: &Expr,
    vars: &[&str],
) -> SparsePolynomial {

    let mut terms = BTreeMap::new();

    collect_terms_recursive(
        expr, vars, &mut terms,
    );

    SparsePolynomial { terms }
}

pub(crate) fn collect_terms_recursive(
    expr: &Expr,
    vars: &[&str],
    terms: &mut BTreeMap<
        Monomial,
        Expr,
    >,
) {

    let simplified =
        simplify(&expr.clone());

    // Convert to AST if it's a DAG to properly match on variants
    let s_expr = if simplified.is_dag()
    {

        simplified
            .to_ast()
            .unwrap_or(simplified)
    } else {

        simplified
    };

    match &s_expr {
        | Expr::Add(a, b) => {

            collect_terms_recursive(
                a, vars, terms,
            );

            collect_terms_recursive(
                b, vars, terms,
            );
        },
        | Expr::Sub(a, b) => {

            collect_terms_recursive(
                a, vars, terms,
            );

            let mut neg_terms =
                BTreeMap::new();

            collect_terms_recursive(
                b,
                vars,
                &mut neg_terms,
            );

            for (mono, coeff) in
                neg_terms
            {

                let entry = terms
                    .entry(mono)
                    .or_insert_with(|| Expr::Constant(0.0));

                *entry = simplify(
                    &Expr::new_sub(
                        entry.clone(),
                        coeff,
                    ),
                );
            }
        },
        | Expr::Mul(a, b) => {

            let mut p1_terms =
                BTreeMap::new();

            collect_terms_recursive(
                a,
                vars,
                &mut p1_terms,
            );

            let mut p2_terms =
                BTreeMap::new();

            collect_terms_recursive(
                b,
                vars,
                &mut p2_terms,
            );

            let p1 = SparsePolynomial {
                terms: p1_terms,
            };

            let p2 = SparsePolynomial {
                terms: p2_terms,
            };

            let product = p1 * p2;

            for (mono, coeff) in
                product.terms
            {

                let entry = terms
                    .entry(mono)
                    .or_insert_with(|| Expr::Constant(0.0));

                *entry = simplify(
                    &Expr::new_add(
                        entry.clone(),
                        coeff,
                    ),
                );
            }
        },
        | Expr::Power(base, exp) => {

            if let Some(e) = as_f64(exp)
            {

                if e.fract() == 0.0
                    && e >= 0.0
                {

                    let mut
                    p_base_terms =
                        BTreeMap::new();

                    collect_terms_recursive(
                        base,
                        vars,
                        &mut p_base_terms,
                    );

                    let p_base = SparsePolynomial {
                        terms: p_base_terms,
                    };

                    let mut result = SparsePolynomial {
                        terms: BTreeMap::from([(
                            Monomial(BTreeMap::new()),
                            Expr::Constant(1.0),
                        )]),
                    };

                    for _ in
                        0..(e as u32)
                    {

                        result = result
                            * p_base
                                .clone(
                                );
                    }

                    for (mono, coeff) in
                        result.terms
                    {

                        let entry = terms
                            .entry(mono)
                            .or_insert_with(|| Expr::Constant(0.0));

                        *entry = simplify(&Expr::new_add(
                            entry.clone(),
                            coeff,
                        ));
                    }

                    return;
                }
            }

            add_term(
                expr,
                &Expr::Constant(1.0),
                terms,
                vars,
            );
        },
        | Expr::Neg(a) => {

            let mut neg_terms =
                BTreeMap::new();

            collect_terms_recursive(
                a,
                vars,
                &mut neg_terms,
            );

            for (mono, coeff) in
                neg_terms
            {

                let entry = terms
                    .entry(mono)
                    .or_insert_with(|| Expr::Constant(0.0));

                *entry = simplify(
                    &Expr::new_sub(
                        entry.clone(),
                        coeff,
                    ),
                );
            }
        },
        | _ => {

            add_term(
                expr,
                &Expr::Constant(1.0),
                terms,
                vars,
            );
        },
    }
}

pub(crate) fn add_term(
    expr: &Expr,
    factor: &Expr,
    terms: &mut BTreeMap<
        Monomial,
        Expr,
    >,
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

        let entry = terms
            .entry(Monomial(
                BTreeMap::new(),
            ))
            .or_insert(Expr::Constant(
                0.0,
            ));

        *entry =
            simplify(&Expr::new_add(
                entry.clone(),
                Expr::new_mul(
                    factor.clone(),
                    expr.clone(),
                ),
            ));

        return;
    }

    if let Expr::Variable(v) = expr {

        if vars.contains(&v.as_str()) {

            let mut mono_map =
                BTreeMap::new();

            mono_map
                .insert(v.clone(), 1);

            let entry = terms
                .entry(Monomial(
                    mono_map,
                ))
                .or_insert(
                    Expr::Constant(0.0),
                );

            *entry = simplify(
                &Expr::new_add(
                    entry.clone(),
                    factor.clone(),
                ),
            );

            return;
        }
    }

    let entry = terms
        .entry(Monomial(
            BTreeMap::new(),
        ))
        .or_insert(Expr::Constant(0.0));

    *entry = simplify(&Expr::new_add(
        entry.clone(),
        Expr::new_mul(
            factor.clone(),
            expr.clone(),
        ),
    ));
}

impl Neg for SparsePolynomial {
    type Output = Self;

    fn neg(self) -> Self {

        let mut new_terms =
            BTreeMap::new();

        for (mono, coeff) in self.terms
        {

            new_terms.insert(
                mono,
                simplify(
                    &Expr::new_neg(
                        coeff,
                    ),
                ),
            );
        }

        Self { terms: new_terms }
    }
}

impl SparsePolynomial {
    #[must_use]

    pub fn eval(
        &self,
        vars: &HashMap<String, f64>,
    ) -> f64 {

        self.terms
            .iter()
            .map(|(mono, coeff)| {

                let coeff_val = eval_expr(coeff, vars);

                let mono_val = mono.0.iter().fold(
                    1.0,
                    |acc, (var, exp)| {

                        let val = vars
                            .get(var)
                            .copied()
                            .unwrap_or(0.0);

                        acc * val.powi(*exp as i32)
                    },
                );

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
#[must_use]

pub fn poly_mul_scalar_expr(
    poly: &SparsePolynomial,
    scalar: &Expr,
) -> SparsePolynomial {

    let mut new_terms = BTreeMap::new();

    for (mono, coeff) in &poly.terms {

        new_terms.insert(
            mono.clone(),
            simplify(&Expr::new_mul(
                coeff.clone(),
                scalar.clone(),
            )),
        );
    }

    SparsePolynomial {
        terms: new_terms,
    }
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
#[must_use]

pub fn gcd(
    mut a: SparsePolynomial,
    mut b: SparsePolynomial,
    var: &str,
) -> SparsePolynomial {

    const MAX_ITERATIONS: usize = 100;

    let mut iterations = 0;

    while !b.terms.is_empty()
        && iterations < MAX_ITERATIONS
    {

        // Check if b is effectively zero (all coefficients are zero)
        let b_is_zero = b
            .terms
            .values()
            .all(|coeff| {

                matches!(coeff, Expr::Constant(c) if c.abs() < 1e-10)
                    || matches!(coeff, Expr::BigInt(n) if n.is_zero())
            });

        if b_is_zero {

            break;
        }

        let (_, remainder) = a
            .long_division(
                b.clone(),
                var,
            );

        // Check if remainder is effectively zero
        if remainder
            .terms
            .is_empty()
        {

            return b;
        }

        // Euclidean algorithm: gcd(a, b) = gcd(b, a mod b)
        a = b;

        b = remainder;

        iterations += 1;
    }

    a
}

pub(crate) fn is_divisible(
    m1: &Monomial,
    m2: &Monomial,
) -> bool {

    m2.0.iter()
        .all(|(var, exp2)| {

            m1.0.get(var)
                .is_some_and(|exp1| {
                    exp1 >= exp2
                })
        })
}

pub(crate) fn subtract_monomials(
    m1: &Monomial,
    m2: &Monomial,
) -> Monomial {

    let mut result = m1.0.clone();

    for (var, exp2) in &m2.0 {

        let exp1 = result
            .entry(var.clone())
            .or_insert(0);

        *exp1 -= exp2;
    }

    Monomial(
        result
            .into_iter()
            .filter(|(_, exp)| *exp > 0)
            .collect(),
    )
}

impl SparsePolynomial {
    #[must_use]

    pub fn degree(
        &self,
        var: &str,
    ) -> isize {

        self.terms
            .keys()
            .map(|m| {

                m.0.get(var)
                    .copied()
                    .unwrap_or(0)
                    as isize
            })
            .max()
            .unwrap_or(-1)
    }

    #[must_use]

    pub fn leading_term(
        &self,
        var: &str,
    ) -> Option<(Monomial, Expr)> {

        self.terms
            .iter()
            .max_by_key(|(m, _)| {

                m.0.get(var)
                    .copied()
                    .unwrap_or(0)
            })
            .map(|(m, c)| {
                (m.clone(), c.clone())
            })
    }

    #[must_use]

    pub fn long_division(
        self,
        divisor: Self,
        var: &str,
    ) -> (Self, Self) {

        if divisor
            .terms
            .is_empty()
        {

            return (
                Self {
                    terms:
                        BTreeMap::new(),
                },
                self,
            );
        }

        let mut quotient = Self {
            terms: BTreeMap::new(),
        };

        let mut remainder = self;

        let divisor_deg =
            divisor.degree(var);

        let mut iterations = 0;

        const MAX_ITERATIONS: usize =
            1000;

        while remainder.degree(var)
            >= divisor_deg
            && iterations
                < MAX_ITERATIONS
        {

            let (lm_d, lc_d) =
                match divisor
                    .leading_term(var)
                {
                    | Some(term) => {
                        term
                    },
                    | None => break,
                };

            let (lm_r, lc_r) =
                match remainder
                    .leading_term(var)
                {
                    | Some(term) => {
                        term
                    },
                    | None => break,
                };

            if !is_divisible(
                &lm_r, &lm_d,
            ) {

                break;
            }

            let t_coeff = simplify(
                &Expr::new_div(
                    lc_r,
                    lc_d.clone(),
                ),
            );

            let t_mono =
                subtract_monomials(
                    &lm_r, &lm_d,
                );

            let mut t = Self {
                terms: BTreeMap::new(),
            };

            t.terms.insert(
                t_mono, t_coeff,
            );

            quotient =
                add_poly(&quotient, &t);

            let sub_term =
                mul_poly(&t, &divisor);

            remainder = subtract_poly(
                &remainder, &sub_term,
            );

            iterations += 1;
        }

        (quotient, remainder)
    }

    #[must_use]

    pub fn get_coeffs_as_vec(
        &self,
        var: &str,
    ) -> Vec<Expr> {

        let deg = self.degree(var);

        if deg < 0 {

            return vec![];
        }

        let mut coeffs =
            vec![
                Expr::Constant(0.0);
                (deg + 1) as usize
            ];

        for (mono, coeff) in &self.terms
        {

            let d = mono
                .0
                .get(var)
                .copied()
                .unwrap_or(0)
                as usize;

            if d < coeffs.len() {

                let mut other_vars =
                    mono.0.clone();

                other_vars.remove(var);

                let term_coeff =
                    if other_vars
                        .is_empty()
                    {

                        coeff.clone()
                    } else {

                        let mut other_terms = BTreeMap::new();

                        other_terms.insert(
                        Monomial(other_vars),
                        coeff.clone(),
                    );

                        sparse_poly_to_expr(&Self { terms: other_terms })
                    };

                coeffs[d] = simplify(
                    &Expr::new_add(
                        coeffs[d]
                            .clone(),
                        term_coeff,
                    ),
                );
            }
        }

        coeffs.reverse();

        coeffs
    }

    #[must_use]

    pub fn get_coeff_for_power(
        &self,
        var: &str,
        power: usize,
    ) -> Option<Expr> {

        let mut mono_map =
            BTreeMap::new();

        if power > 0 {

            mono_map.insert(
                var.to_string(),
                power as u32,
            );
        }

        self.terms
            .get(&Monomial(mono_map))
            .cloned()
    }

    /// Removes terms with zero coefficients from the polynomial.

    pub fn prune_zeros(&mut self) {

        self.terms.retain(
            |_, coeff| {

                !is_zero(&simplify(
                    &coeff.clone(),
                ))
            },
        );
    }
}

impl Add for SparsePolynomial {
    type Output = Self;

    fn add(
        self,
        rhs: Self,
    ) -> Self {

        add_poly(&self, &rhs)
    }
}

impl Sub for SparsePolynomial {
    type Output = Self;

    fn sub(
        self,
        rhs: Self,
    ) -> Self {

        let neg_rhs = mul_poly(
            &rhs,
            &poly_from_coeffs(
                &[Expr::Constant(-1.0)],
                "",
            ),
        );

        add_poly(&self, &neg_rhs)
    }
}

impl Mul for SparsePolynomial {
    type Output = Self;

    fn mul(
        self,
        rhs: Self,
    ) -> Self {

        mul_poly(&self, &rhs)
    }
}

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
#[must_use]

pub fn poly_from_coeffs(
    coeffs: &[Expr],
    var: &str,
) -> SparsePolynomial {

    let mut terms = BTreeMap::new();

    let n = coeffs.len() - 1;

    for (i, coeff) in coeffs
        .iter()
        .enumerate()
    {

        if !is_zero(&simplify(
            &coeff.clone(),
        )) {

            let mut mono_map =
                BTreeMap::new();

            let power = (n - i) as u32;

            if power > 0 {

                mono_map.insert(
                    var.to_string(),
                    power,
                );
            }

            terms.insert(
                Monomial(mono_map),
                coeff.clone(),
            );
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
#[must_use]

pub fn sparse_poly_to_expr(
    poly: &SparsePolynomial
) -> Expr {

    let mut total_expr =
        Expr::Constant(0.0);

    for (mono, coeff) in &poly.terms {

        let mut term_expr =
            coeff.clone();

        for (var_name, &exp) in &mono.0
        {

            if exp > 0 {

                let var_expr =
                    Expr::new_pow(
                        Expr::Variable(
                            var_name
                                .clone(
                                ),
                        ),
                        Expr::Constant(
                            f64::from(
                                exp,
                            ),
                        ),
                    );

                term_expr = simplify(
                    &Expr::new_mul(
                        term_expr,
                        var_expr,
                    ),
                );
            }
        }

        total_expr =
            simplify(&Expr::new_add(
                total_expr, term_expr,
            ));
    }

    total_expr
}
