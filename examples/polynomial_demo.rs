//! Polynomial Manipulation Example
//!
//! This example demonstrates various polynomial operations including:
//! - Creating polynomials in different representations
//! - Polynomial arithmetic (addition, multiplication, division)
//! - Differentiation
//! - GCD computation
//! - Representation conversions

use rssn::symbolic::core::{Expr, Monomial, SparsePolynomial};
use rssn::symbolic::polynomial::*;
use std::collections::BTreeMap;

fn main() {

    println!("=== Polynomial Manipulation Examples ===\n");

    // Example 1: Creating polynomials
    println!("1. Creating Polynomials");

    println!("   ---------------------");

    // Create x^2 + 2x + 1 using expression form
    let poly1 = Expr::new_add(
        Expr::new_add(
            Expr::new_pow(Expr::new_variable("x"), Expr::new_constant(2.0)),
            Expr::new_mul(Expr::new_constant(2.0), Expr::new_variable("x")),
        ),
        Expr::new_constant(1.0),
    );

    println!("   Polynomial 1: x^2 + 2x + 1");

    println!("   Degree: {}", polynomial_degree(&poly1, "x"));

    println!("   Is polynomial in x: {}\n", is_polynomial(&poly1, "x"));

    // Example 2: Polynomial long division
    println!("2. Polynomial Long Division");

    println!("   -------------------------");

    // Divide x^2 + 3x + 2 by x + 1
    let dividend = Expr::new_add(
        Expr::new_add(
            Expr::new_pow(Expr::new_variable("x"), Expr::new_constant(2.0)),
            Expr::new_mul(Expr::new_constant(3.0), Expr::new_variable("x")),
        ),
        Expr::new_constant(2.0),
    );

    let divisor = Expr::new_add(Expr::new_variable("x"), Expr::new_constant(1.0));

    println!("   Dividend: x^2 + 3x + 2");

    println!("   Divisor: x + 1");

    let (quotient, remainder) = polynomial_long_division(&dividend, &divisor, "x");

    println!("   Quotient: {}", quotient);

    println!("   Remainder: {}\n", remainder);

    // Example 3: Differentiation
    println!("3. Polynomial Differentiation");

    println!("   ---------------------------");

    // d/dx(x^3 + 2x^2 + x)
    let poly_to_diff = Expr::new_add(
        Expr::new_add(
            Expr::new_pow(Expr::new_variable("x"), Expr::new_constant(3.0)),
            Expr::new_mul(
                Expr::new_constant(2.0),
                Expr::new_pow(Expr::new_variable("x"), Expr::new_constant(2.0)),
            ),
        ),
        Expr::new_variable("x"),
    );

    println!("   Original: x^3 + 2x^2 + x");

    // Convert to sparse form for differentiation
    let sparse = expr_to_sparse_poly(&poly_to_diff, &["x"]);

    let derivative = differentiate_poly(&sparse, "x");

    let derivative_expr = sparse_poly_to_expr(&derivative);

    println!("   Derivative: {}\n", derivative_expr);

    // Example 4: Sparse polynomial operations
    println!("4. Sparse Polynomial Operations");

    println!("   -----------------------------");

    // Create two sparse polynomials: x + 1 and x + 2
    let mut terms1 = BTreeMap::new();

    let mut mono_x1 = BTreeMap::new();

    mono_x1.insert("x".to_string(), 1);

    terms1.insert(Monomial(mono_x1), Expr::Constant(1.0));

    terms1.insert(Monomial(BTreeMap::new()), Expr::Constant(1.0));

    let mut terms2 = BTreeMap::new();

    let mut mono_x2 = BTreeMap::new();

    mono_x2.insert("x".to_string(), 1);

    terms2.insert(Monomial(mono_x2), Expr::Constant(1.0));

    terms2.insert(Monomial(BTreeMap::new()), Expr::Constant(2.0));

    let p1 = SparsePolynomial { terms: terms1 };

    let p2 = SparsePolynomial { terms: terms2 };

    println!("   p1: x + 1");

    println!("   p2: x + 2");

    let sum = add_poly(&p1, &p2);

    let product = mul_poly(&p1, &p2);

    println!("   p1 + p2: {}", sparse_poly_to_expr(&sum));

    println!("   p1 * p2: {}\n", sparse_poly_to_expr(&product));

    // Example 5: Multivariate polynomials
    println!("5. Multivariate Polynomials");

    println!("   -------------------------");

    // Create x*y + 2x + 3y + 4
    let multivar = Expr::new_add(
        Expr::new_add(
            Expr::new_add(
                Expr::new_mul(Expr::new_variable("x"), Expr::new_variable("y")),
                Expr::new_mul(Expr::new_constant(2.0), Expr::new_variable("x")),
            ),
            Expr::new_mul(Expr::new_constant(3.0), Expr::new_variable("y")),
        ),
        Expr::new_constant(4.0),
    );

    println!("   Polynomial: x*y + 2x + 3y + 4");

    let sparse_multi = expr_to_sparse_poly(&multivar, &["x", "y"]);

    println!("   Number of terms: {}", sparse_multi.terms.len());

    println!("   Degree in x: {}", sparse_multi.degree("x"));

    println!("   Degree in y: {}\n", sparse_multi.degree("y"));

    // Example 6: Coefficient extraction
    println!("6. Coefficient Extraction");

    println!("   -----------------------");

    let poly_coeffs = Expr::new_add(
        Expr::new_add(
            Expr::new_mul(
                Expr::new_constant(3.0),
                Expr::new_pow(Expr::new_variable("x"), Expr::new_constant(2.0)),
            ),
            Expr::new_mul(Expr::new_constant(2.0), Expr::new_variable("x")),
        ),
        Expr::new_constant(1.0),
    );

    println!("   Polynomial: 3x^2 + 2x + 1");

    let coeffs = to_polynomial_coeffs_vec(&poly_coeffs, "x");

    println!("   Coefficients [c0, c1, c2]: {:?}", coeffs);

    // Reconstruct from coefficients
    let reconstructed = from_coeffs_to_expr(&coeffs, "x");

    println!("   Reconstructed: {}\n", reconstructed);

    // Example 7: GCD computation
    println!("7. GCD Computation");

    println!("   ----------------");

    // GCD of x^2 - 1 and x^2 - 2x + 1
    let gcd_poly1 = Expr::new_sub(
        Expr::new_pow(Expr::new_variable("x"), Expr::new_constant(2.0)),
        Expr::new_constant(1.0),
    );

    let gcd_poly2 = Expr::new_add(
        Expr::new_sub(
            Expr::new_pow(Expr::new_variable("x"), Expr::new_constant(2.0)),
            Expr::new_mul(Expr::new_constant(2.0), Expr::new_variable("x")),
        ),
        Expr::new_constant(1.0),
    );

    println!("   p1: x^2 - 1");

    println!("   p2: x^2 - 2x + 1");

    let sparse1 = expr_to_sparse_poly(&gcd_poly1, &["x"]);

    let sparse2 = expr_to_sparse_poly(&gcd_poly2, &["x"]);

    let gcd_result = gcd(sparse1, sparse2, "x");

    println!("   GCD: {}\n", sparse_poly_to_expr(&gcd_result));

    println!("=== Examples Complete ===");
}
