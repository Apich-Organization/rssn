use num_bigint::BigInt;
use rssn::symbolic::finite_field::*;
use std::sync::Arc;

fn main() {

    let field = PrimeField::new(BigInt::from(5));

    // p1 = x + 2
    let p1 = FiniteFieldPolynomial::new(
        vec![
            PrimeFieldElement::new(BigInt::from(1), field.clone()),
            PrimeFieldElement::new(BigInt::from(2), field.clone()),
        ],
        field.clone(),
    );

    // p2 = 2x + 3
    let p2 = FiniteFieldPolynomial::new(
        vec![
            PrimeFieldElement::new(BigInt::from(2), field.clone()),
            PrimeFieldElement::new(BigInt::from(3), field.clone()),
        ],
        field.clone(),
    );

    // p1 + p2 = 3x + 5 = 3x (since 5 mod 5 = 0)
    let result = p1 + p2;

    println!(
        "Addition result: degree={}, coeffs={:?}",
        result.degree(),
        result.coeffs.iter().map(|c| &c.value).collect::<Vec<_>>()
    );

    // Test division
    let field2 = PrimeField::new(BigInt::from(5));

    // dividend = x^2 + 2x + 3
    let dividend = FiniteFieldPolynomial::new(
        vec![
            PrimeFieldElement::new(BigInt::from(1), field2.clone()),
            PrimeFieldElement::new(BigInt::from(2), field2.clone()),
            PrimeFieldElement::new(BigInt::from(3), field2.clone()),
        ],
        field2.clone(),
    );

    // divisor = x + 1
    let divisor = FiniteFieldPolynomial::new(
        vec![
            PrimeFieldElement::new(BigInt::from(1), field2.clone()),
            PrimeFieldElement::new(BigInt::from(1), field2.clone()),
        ],
        field2.clone(),
    );

    let (quotient, remainder) = dividend
        .long_division(&divisor)
        .expect("Division should succeed");

    println!(
        "Division: quotient degree={}, remainder degree={}",
        quotient.degree(),
        remainder.degree()
    );

    println!(
        "Quotient coeffs: {:?}",
        quotient.coeffs.iter().map(|c| &c.value).collect::<Vec<_>>()
    );

    println!(
        "Remainder coeffs: {:?}",
        remainder
            .coeffs
            .iter()
            .map(|c| &c.value)
            .collect::<Vec<_>>()
    );
}
