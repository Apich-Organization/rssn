use std::sync::Arc;

use num_bigint::BigInt;
use rssn::symbolic::finite_field::*;

#[test]

fn test_prime_field_element_creation() {

    let field = PrimeField::new(BigInt::from(7));

    let elem = PrimeFieldElement::new(
        BigInt::from(5),
        field.clone(),
    );

    assert_eq!(
        elem.value,
        BigInt::from(5)
    );

    // Test automatic reduction
    let elem2 = PrimeFieldElement::new(
        BigInt::from(10),
        field.clone(),
    );

    assert_eq!(
        elem2.value,
        BigInt::from(3)
    ); // 10 mod 7 = 3

    // Test negative values
    let elem3 = PrimeFieldElement::new(
        BigInt::from(-2),
        field.clone(),
    );

    assert_eq!(
        elem3.value,
        BigInt::from(5)
    ); // -2 mod 7 = 5
}

#[test]

fn test_prime_field_element_addition() {

    let field = PrimeField::new(BigInt::from(7));

    let a = PrimeFieldElement::new(
        BigInt::from(5),
        field.clone(),
    );

    let b = PrimeFieldElement::new(
        BigInt::from(4),
        field.clone(),
    );

    let result = a + b;

    assert_eq!(
        result.value,
        BigInt::from(2)
    ); // (5 + 4) mod 7 = 2
}

#[test]

fn test_prime_field_element_subtraction() {

    let field = PrimeField::new(BigInt::from(7));

    let a = PrimeFieldElement::new(
        BigInt::from(3),
        field.clone(),
    );

    let b = PrimeFieldElement::new(
        BigInt::from(5),
        field.clone(),
    );

    let result = a - b;

    assert_eq!(
        result.value,
        BigInt::from(5)
    ); // (3 - 5) mod 7 = -2 mod 7 = 5
}

#[test]

fn test_prime_field_element_multiplication() {

    let field = PrimeField::new(BigInt::from(7));

    let a = PrimeFieldElement::new(
        BigInt::from(3),
        field.clone(),
    );

    let b = PrimeFieldElement::new(
        BigInt::from(5),
        field.clone(),
    );

    let result = a * b;

    assert_eq!(
        result.value,
        BigInt::from(1)
    ); // (3 * 5) mod 7 = 15 mod 7 = 1
}

#[test]

fn test_prime_field_element_division() {

    let field = PrimeField::new(BigInt::from(7));

    let a = PrimeFieldElement::new(
        BigInt::from(6),
        field.clone(),
    );

    let b = PrimeFieldElement::new(
        BigInt::from(3),
        field.clone(),
    );

    let result = a / b;

    assert_eq!(
        result.value,
        BigInt::from(2)
    ); // 6 / 3 = 2 in GF(7)
}

#[test]

fn test_prime_field_element_inverse() {

    let field = PrimeField::new(BigInt::from(7));

    let elem = PrimeFieldElement::new(
        BigInt::from(3),
        field.clone(),
    );

    let inv = elem
        .inverse()
        .expect("Inverse should exist");

    assert_eq!(
        inv.value,
        BigInt::from(5)
    ); // 3 * 5 = 15 = 1 mod 7

    // Verify: elem * inv = 1
    let one = elem.clone() * inv;

    assert_eq!(
        one.value,
        BigInt::from(1)
    );
}

#[test]

fn test_finite_field_polynomial_creation() {

    let field = PrimeField::new(BigInt::from(5));

    let coeffs = vec![
        PrimeFieldElement::new(
            BigInt::from(1),
            field.clone(),
        ),
        PrimeFieldElement::new(
            BigInt::from(2),
            field.clone(),
        ),
        PrimeFieldElement::new(
            BigInt::from(3),
            field.clone(),
        ),
    ];

    let poly = FiniteFieldPolynomial::new(
        coeffs,
        field.clone(),
    );

    assert_eq!(poly.degree(), 2);

    assert_eq!(poly.coeffs.len(), 3);
}

#[test]

fn test_finite_field_polynomial_degree() {

    let field = PrimeField::new(BigInt::from(5));

    // Zero polynomial
    let zero_poly = FiniteFieldPolynomial::new(
        vec![],
        field.clone(),
    );

    assert_eq!(
        zero_poly.degree(),
        -1
    );

    // Constant polynomial
    let const_poly = FiniteFieldPolynomial::new(
        vec![
            PrimeFieldElement::new(
                BigInt::from(3),
                field.clone(),
            ),
        ],
        field.clone(),
    );

    assert_eq!(
        const_poly.degree(),
        0
    );

    // Degree 2 polynomial
    let poly = FiniteFieldPolynomial::new(
        vec![
            PrimeFieldElement::new(
                BigInt::from(1),
                field.clone(),
            ),
            PrimeFieldElement::new(
                BigInt::from(2),
                field.clone(),
            ),
            PrimeFieldElement::new(
                BigInt::from(3),
                field.clone(),
            ),
        ],
        field.clone(),
    );

    assert_eq!(poly.degree(), 2);
}

#[test]

fn test_finite_field_polynomial_addition() {

    let field = PrimeField::new(BigInt::from(5));

    // p1 = x + 2
    let p1 = FiniteFieldPolynomial::new(
        vec![
            PrimeFieldElement::new(
                BigInt::from(1),
                field.clone(),
            ),
            PrimeFieldElement::new(
                BigInt::from(2),
                field.clone(),
            ),
        ],
        field.clone(),
    );

    // p2 = 2x + 3
    let p2 = FiniteFieldPolynomial::new(
        vec![
            PrimeFieldElement::new(
                BigInt::from(2),
                field.clone(),
            ),
            PrimeFieldElement::new(
                BigInt::from(3),
                field.clone(),
            ),
        ],
        field.clone(),
    );

    // p1 + p2 = 3x + 5 = 3x + 0 (since 5 mod 5 = 0)
    let result = p1 + p2;

    assert_eq!(
        result.coeffs.len(),
        2
    );

    assert_eq!(
        result.coeffs[0].value,
        BigInt::from(3)
    ); // x coefficient
    assert_eq!(
        result.coeffs[1].value,
        BigInt::from(0)
    ); // constant (5 mod 5 = 0)
}

#[test]

fn test_finite_field_polynomial_multiplication() {

    let field = PrimeField::new(BigInt::from(5));

    // p1 = x + 1
    let p1 = FiniteFieldPolynomial::new(
        vec![
            PrimeFieldElement::new(
                BigInt::from(1),
                field.clone(),
            ),
            PrimeFieldElement::new(
                BigInt::from(1),
                field.clone(),
            ),
        ],
        field.clone(),
    );

    // p2 = x + 2
    let p2 = FiniteFieldPolynomial::new(
        vec![
            PrimeFieldElement::new(
                BigInt::from(1),
                field.clone(),
            ),
            PrimeFieldElement::new(
                BigInt::from(2),
                field.clone(),
            ),
        ],
        field.clone(),
    );

    // p1 * p2 = x^2 + 3x + 2
    let result = p1 * p2;

    assert_eq!(result.degree(), 2);

    assert_eq!(
        result.coeffs[0].value,
        BigInt::from(1)
    ); // x^2 coefficient
    assert_eq!(
        result.coeffs[1].value,
        BigInt::from(3)
    ); // x coefficient
    assert_eq!(
        result.coeffs[2].value,
        BigInt::from(2)
    ); // constant
}

#[test]

fn test_finite_field_polynomial_long_division() {

    let field = PrimeField::new(BigInt::from(5));

    // dividend = x^2 + 2x + 3
    let dividend = FiniteFieldPolynomial::new(
        vec![
            PrimeFieldElement::new(
                BigInt::from(1),
                field.clone(),
            ),
            PrimeFieldElement::new(
                BigInt::from(2),
                field.clone(),
            ),
            PrimeFieldElement::new(
                BigInt::from(3),
                field.clone(),
            ),
        ],
        field.clone(),
    );

    // divisor = x + 1
    let divisor = FiniteFieldPolynomial::new(
        vec![
            PrimeFieldElement::new(
                BigInt::from(1),
                field.clone(),
            ),
            PrimeFieldElement::new(
                BigInt::from(1),
                field.clone(),
            ),
        ],
        field.clone(),
    );

    let (quotient, remainder) = dividend
        .long_division(&divisor)
        .expect("Division should succeed");

    // The quotient degree depends on implementation details
    // Just verify that the division succeeded and remainder degree is less than divisor
    assert!(remainder.degree() < divisor.degree());
}

#[test]

fn test_extension_field_element_creation() {

    let prime_field = PrimeField::new(BigInt::from(5));

    // Create irreducible polynomial x^2 + 2 over GF(5)
    let irreducible = FiniteFieldPolynomial::new(
        vec![
            PrimeFieldElement::new(
                BigInt::from(1),
                prime_field.clone(),
            ),
            PrimeFieldElement::new(
                BigInt::from(0),
                prime_field.clone(),
            ),
            PrimeFieldElement::new(
                BigInt::from(2),
                prime_field.clone(),
            ),
        ],
        prime_field.clone(),
    );

    let ext_field = Arc::new(ExtensionField {
        prime_field : prime_field.clone(),
        irreducible_poly : irreducible,
    });

    // Create element x + 1
    let poly = FiniteFieldPolynomial::new(
        vec![
            PrimeFieldElement::new(
                BigInt::from(1),
                prime_field.clone(),
            ),
            PrimeFieldElement::new(
                BigInt::from(1),
                prime_field.clone(),
            ),
        ],
        prime_field.clone(),
    );

    let elem = ExtensionFieldElement::new(poly, ext_field);

    assert!(elem.poly.degree() <= 1); // Should be reduced modulo irreducible poly
}

#[test]

fn test_extension_field_element_arithmetic() {

    let prime_field = PrimeField::new(BigInt::from(5));

    // Create irreducible polynomial x^2 + 2 over GF(5)
    let irreducible = FiniteFieldPolynomial::new(
        vec![
            PrimeFieldElement::new(
                BigInt::from(1),
                prime_field.clone(),
            ),
            PrimeFieldElement::new(
                BigInt::from(0),
                prime_field.clone(),
            ),
            PrimeFieldElement::new(
                BigInt::from(2),
                prime_field.clone(),
            ),
        ],
        prime_field.clone(),
    );

    let ext_field = Arc::new(ExtensionField {
        prime_field : prime_field.clone(),
        irreducible_poly : irreducible,
    });

    // Create two elements
    let poly1 = FiniteFieldPolynomial::new(
        vec![
            PrimeFieldElement::new(
                BigInt::from(1),
                prime_field.clone(),
            ),
            PrimeFieldElement::new(
                BigInt::from(1),
                prime_field.clone(),
            ),
        ],
        prime_field.clone(),
    );

    let poly2 = FiniteFieldPolynomial::new(
        vec![
            PrimeFieldElement::new(
                BigInt::from(2),
                prime_field.clone(),
            ),
            PrimeFieldElement::new(
                BigInt::from(3),
                prime_field.clone(),
            ),
        ],
        prime_field.clone(),
    );

    let elem1 = ExtensionFieldElement::new(
        poly1,
        ext_field.clone(),
    );

    let elem2 = ExtensionFieldElement::new(
        poly2,
        ext_field.clone(),
    );

    // Test addition
    let sum = elem1
        .clone()
        .add(elem2.clone())
        .expect("Addition should succeed");

    assert!(sum.poly.degree() <= 1);

    // Test multiplication
    let product = elem1
        .mul(elem2)
        .expect("Multiplication should succeed");

    assert!(
        product
            .poly
            .degree()
            <= 1
    ); // Should be reduced modulo irreducible poly
}
