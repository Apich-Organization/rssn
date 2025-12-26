use num_bigint::BigInt;
use num_traits::{One, Zero};
use rssn::symbolic::finite_field::{FiniteFieldPolynomial, PrimeField, PrimeFieldElement};
use rssn::symbolic::poly_factorization::*;
use std::sync::Arc;

fn create_test_field(modulus: i64) -> Arc<PrimeField> {

    PrimeField::new(BigInt::from(modulus))
}

fn create_poly(
    coeffs: Vec<i64>,
    field: Arc<PrimeField>,
) -> FiniteFieldPolynomial {

    let elements: Vec<PrimeFieldElement> = coeffs
        .into_iter()
        .map(|c| PrimeFieldElement::new(BigInt::from(c), field.clone()))
        .collect();

    FiniteFieldPolynomial::new(elements, field)
}

#[test]

fn test_poly_derivative_gf() {

    // Test derivative of x^2 + 2x + 1 over GF(5)
    // Derivative should be 2x + 2
    let field = create_test_field(5);

    let poly = create_poly(vec![1, 2, 1], field.clone()); // x^2 + 2x + 1

    let derivative = poly_derivative_gf(&poly);

    // Derivative: 2x + 2
    assert_eq!(derivative.degree(), 1);

    assert_eq!(
        derivative
            .coeffs
            .len(),
        2
    );

    assert_eq!(derivative.coeffs[0].value, BigInt::from(2)); // coefficient of x
    assert_eq!(derivative.coeffs[1].value, BigInt::from(2)); // constant term
}

#[test]

fn test_poly_gcd_gf() {

    // Test GCD of (x^2 - 1) and (x - 1) over GF(5)
    // x^2 - 1 = (x-1)(x+1), so GCD should be (x-1)
    let field = create_test_field(5);

    let poly1 = create_poly(vec![1, 0, -1], field.clone()); // x^2 - 1
    let poly2 = create_poly(vec![1, -1], field.clone()); // x - 1

    let gcd = poly_gcd_gf(poly1, poly2).unwrap();

    // GCD should be x - 1 (or a scalar multiple)
    assert_eq!(gcd.degree(), 1);
}

#[test]

fn test_square_free_factorization_simple() {

    // Test with x^2 over GF(3)
    // x^2 = x * x, so it's not square-free
    let field = create_test_field(3);

    let poly = create_poly(vec![1, 0, 0], field.clone()); // x^2

    let result = square_free_factorization_gf(poly);

    assert!(result.is_ok());

    let factors = result.unwrap();

    // Should have at least one factor
    assert!(!factors.is_empty());
}

#[test]

fn test_factor_gf_linear() {

    // Test factoring a linear polynomial x - 1 over GF(5)
    // Should return itself as it's already irreducible
    let field = create_test_field(5);

    let poly = create_poly(vec![1, -1], field.clone()); // x - 1

    let factors = factor_gf(&poly).unwrap();

    assert_eq!(factors.len(), 1);

    assert_eq!(factors[0].degree(), 1);
}

#[test]

fn test_factor_gf_quadratic() {

    // Test factoring x^2 - 1 over GF(5)
    // x^2 - 1 = (x-1)(x+1) = (x-1)(x+4) in GF(5)
    let field = create_test_field(5);

    let poly = create_poly(vec![1, 0, -1], field.clone()); // x^2 - 1

    let factors = factor_gf(&poly).unwrap();

    // Should factor into two linear factors
    assert!(factors.len() >= 2 || factors[0].degree() == 2);
}

#[test]

fn test_berlekamp_factorization_simple() {

    // Test Berlekamp on x^2 - 1 over GF(3)
    let field = create_test_field(3);

    let poly = create_poly(vec![1, 0, -1], field.clone()); // x^2 - 1

    let factors = berlekamp_factorization(&poly).unwrap();

    // Should produce factors
    assert!(!factors.is_empty());
}

#[test]

fn test_distinct_degree_factorization() {

    // Test DDF on x^3 - x over GF(2)
    // x^3 - x = x(x^2 - 1) = x(x-1)(x+1) in GF(2)
    let field = create_test_field(2);

    let poly = create_poly(vec![1, 0, -1, 0], field.clone()); // x^3 - x

    let result = distinct_degree_factorization(&poly);

    assert!(result.is_ok());

    let ddf_factors = result.unwrap();

    // Should have factors grouped by degree
    assert!(!ddf_factors.is_empty());
}

#[test]

fn test_cantor_zassenhaus() {

    // Test Cantor-Zassenhaus on a polynomial over a larger field
    let field = create_test_field(101); // Larger prime
    let poly = create_poly(vec![1, 0, -1], field.clone()); // x^2 - 1

    let factors = cantor_zassenhaus(&poly).unwrap();

    // Should produce factors
    assert!(!factors.is_empty());
}

#[test]

fn test_poly_pow_mod() {

    // Test x^2 mod (x^2 + 1) over GF(5)
    let field = create_test_field(5);

    let base = create_poly(vec![1, 0], field.clone()); // x
    let modulus = create_poly(vec![1, 0, 1], field.clone()); // x^2 + 1
    let exp = BigInt::from(2);

    let result = poly_pow_mod(base, &exp, &modulus).unwrap();

    // x^2 mod (x^2 + 1) = -1 = 4 in GF(5)
    assert_eq!(result.degree(), 0);

    assert_eq!(result.coeffs[0].value, BigInt::from(4));
}

#[test]

fn test_poly_mul_scalar() {

    // Test multiplying polynomial by scalar
    let field = create_test_field(7);

    let poly = create_poly(vec![1, 2, 3], field.clone()); // x^2 + 2x + 3
    let scalar = BigInt::from(2);

    let result = poly_mul_scalar(&poly, &scalar);

    // Should be 2x^2 + 4x + 6
    assert_eq!(result.coeffs[0].value, BigInt::from(2));

    assert_eq!(result.coeffs[1].value, BigInt::from(4));

    assert_eq!(result.coeffs[2].value, BigInt::from(6));
}

#[test]

fn test_poly_extended_gcd() {

    // Test extended GCD: find s, t such that a*s + b*t = gcd(a, b)
    let field = create_test_field(5);

    let a = create_poly(vec![1, 0, -1], field.clone()); // x^2 - 1
    let b = create_poly(vec![1, -1], field.clone()); // x - 1

    let result = poly_extended_gcd(a.clone(), b.clone());

    assert!(result.is_ok());

    let (gcd, s, t) = result.unwrap();

    // Verify: a*s + b*t = gcd
    let left_side = a * s + b * t;

    // The result should equal the GCD (up to scalar multiple)
    assert_eq!(left_side.degree(), gcd.degree());
}

#[test]

fn test_factor_gf_irreducible() {

    // Test factoring an irreducible polynomial
    // x^2 + x + 1 is irreducible over GF(2)
    let field = create_test_field(2);

    let poly = create_poly(vec![1, 1, 1], field.clone()); // x^2 + x + 1

    let factors = factor_gf(&poly).unwrap();

    // Should return itself as it's irreducible
    assert_eq!(factors.len(), 1);

    assert_eq!(factors[0].degree(), 2);
}

#[test]

fn test_zero_polynomial_derivative() {

    // Test derivative of zero polynomial
    let field = create_test_field(5);

    let poly = create_poly(vec![], field.clone());

    let derivative = poly_derivative_gf(&poly);

    assert_eq!(
        derivative
            .coeffs
            .len(),
        0
    );
}

#[test]

fn test_constant_polynomial_derivative() {

    // Test derivative of constant polynomial
    let field = create_test_field(5);

    let poly = create_poly(vec![3], field.clone()); // constant 3

    let derivative = poly_derivative_gf(&poly);

    // Derivative of constant is 0
    assert_eq!(
        derivative
            .coeffs
            .len(),
        0
    );
}
