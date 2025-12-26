//! Test suite for error_correction_helper module (finite field operations).

use num_bigint::BigInt;
use rssn::symbolic::core::Expr;
use rssn::symbolic::error_correction_helper::*;

#[test]

fn test_gf256_add() {

    // Addition in GF(2^8) is XOR
    assert_eq!(gf256_add(0, 0), 0);

    assert_eq!(
        gf256_add(0x12, 0x34),
        0x12 ^ 0x34
    );

    assert_eq!(
        gf256_add(0xFF, 0xFF),
        0
    );
}

#[test]

fn test_gf256_mul() {

    // Multiplication in GF(2^8)
    assert_eq!(gf256_mul(0, 5), 0);

    assert_eq!(gf256_mul(5, 0), 0);

    assert_eq!(gf256_mul(1, 5), 5);

    assert_eq!(gf256_mul(5, 1), 5);

    // Check associativity: (a*b)*c == a*(b*c)
    let a = 0x53u8;

    let b = 0xCAu8;

    let c = 0x02u8;

    assert_eq!(
        gf256_mul(gf256_mul(a, b), c),
        gf256_mul(a, gf256_mul(b, c))
    );
}

#[test]

fn test_gf256_inv() {

    // 0 has no inverse
    assert!(gf256_inv(0).is_err());

    // a * inv(a) == 1
    for a in 1u8..=255u8 {

        let inv_a =
            gf256_inv(a).unwrap();

        assert_eq!(
            gf256_mul(a, inv_a),
            1,
            "Inverse of {} failed",
            a
        );
    }
}

#[test]

fn test_gf256_div() {

    // Division by zero fails
    assert!(gf256_div(5, 0).is_err());

    // a / 1 == a
    assert_eq!(
        gf256_div(5, 1).unwrap(),
        5
    );

    // a / a == 1
    for a in 1u8..=255u8 {

        assert_eq!(
            gf256_div(a, a).unwrap(),
            1
        );
    }
}

#[test]

fn test_gf256_log() {

    // Log of 0 is undefined
    assert!(gf256_log(0).is_err());

    // log(1) = 0 (since alpha^0 = 1)
    assert_eq!(
        gf256_log(1).unwrap(),
        0
    );

    // exp(log(a)) = a for all non-zero a
    for a in 1u8..=255u8 {

        let log_a =
            gf256_log(a).unwrap();

        assert_eq!(gf256_exp(log_a), a);
    }
}

#[test]

fn test_gf256_pow() {

    // a^0 = 1 for any a
    assert_eq!(gf256_pow(5, 0), 1);

    assert_eq!(gf256_pow(0, 0), 1);

    // a^1 = a
    assert_eq!(gf256_pow(5, 1), 5);

    // 0^n = 0 for n > 0
    assert_eq!(gf256_pow(0, 5), 0);

    // a^2 = a * a
    let a = 0x53u8;

    assert_eq!(
        gf256_pow(a, 2),
        gf256_mul(a, a)
    );

    // a^3 = a * a * a
    assert_eq!(
        gf256_pow(a, 3),
        gf256_mul(gf256_mul(a, a), a)
    );
}

#[test]

fn test_poly_eval_gf256() {

    // p(x) = 1 + x + x^2 evaluated at x = 2
    // In GF(256): coefficients are [1, 1, 1] (highest degree first: x^2 + x + 1)
    // Horner: ((1 * 2) ^ 1) * 2 ^ 1 = (2 ^ 1) * 2 ^ 1 = 3 * 2 ^ 1 = 6 ^ 1 = 7
    let poly = vec![1u8, 1, 1]; // x^2 + x + 1
    let result =
        poly_eval_gf256(&poly, 2);

    // Manual check: 1*4 ^ 1*2 ^ 1 = 4 ^ 2 ^ 1 = 7
    assert_eq!(result, 7);
}

#[test]

fn test_poly_add_gf256() {

    // (x + 1) + (x^2 + 1) = x^2 + x
    let p1 = vec![1u8, 1]; // x + 1
    let p2 = vec![1u8, 0, 1]; // x^2 + 1
    let result =
        poly_add_gf256(&p1, &p2);

    // Result: x^2 + x + (1 ^ 1) = x^2 + x = [1, 1, 0]
    assert_eq!(
        result,
        vec![1u8, 1, 0]
    );
}

#[test]

fn test_poly_mul_gf256() {

    // (x + 1) * (x + 1) = x^2 + 2x + 1 = x^2 + 1 (since 2 = 0 in GF(2))
    let p1 = vec![1u8, 1]; // x + 1
    let p2 = vec![1u8, 1]; // x + 1
    let result =
        poly_mul_gf256(&p1, &p2);

    // x^2 + 2x + 1 => [1, 0, 1] in GF(2^8)
    assert_eq!(
        result,
        vec![1u8, 0, 1]
    );
}

#[test]

fn test_poly_scale_gf256() {

    let poly = vec![1u8, 2, 3];

    let scaled =
        poly_scale_gf256(&poly, 0);

    assert_eq!(
        scaled,
        vec![0, 0, 0]
    );

    let scaled1 =
        poly_scale_gf256(&poly, 1);

    assert_eq!(scaled1, poly);
}

#[test]

fn test_poly_derivative_gf256() {

    // Derivative of x^2 + x + 1 = 2x + 1 = 0 + 1 = 1 (in char 2)
    // Actually in GF(2^8), d/dx[a*x^n] = n*a*x^(n-1), and n mod 2 matters
    let poly = vec![1u8, 1, 1]; // x^2 + x + 1
    let deriv =
        poly_derivative_gf256(&poly);

    // x^2 has power 2 (even, vanishes), x has power 1 (odd, survives)
    assert!(!deriv.is_empty());
}

#[test]

fn test_poly_gcd_gf256() {

    // GCD of p and p should be p (monic)
    let p = vec![1u8, 1]; // x + 1
    let gcd = poly_gcd_gf256(&p, &p);

    // Should be monic version of p
    assert!(!gcd.is_empty());

    assert_eq!(gcd[0], 1); // Monic

    // GCD of a polynomial and 1 is 1
    let one = vec![1u8];

    let gcd2 = poly_gcd_gf256(&p, &one);

    assert_eq!(gcd2, vec![1u8]);
}

#[test]

fn test_finite_field_element_arithmetic(
) {

    let field = FiniteField::new(7); // GF(7)
    let a = FieldElement::new(
        BigInt::from(3),
        field.clone(),
    );

    let b = FieldElement::new(
        BigInt::from(5),
        field.clone(),
    );

    // 3 + 5 = 8 = 1 (mod 7)
    let sum = (a.clone() + b.clone())
        .unwrap();

    assert_eq!(
        sum.value,
        BigInt::from(1)
    );

    // 3 * 5 = 15 = 1 (mod 7)
    let prod = (a.clone() * b.clone())
        .unwrap();

    assert_eq!(
        prod.value,
        BigInt::from(1)
    );

    // 3 - 5 = -2 = 5 (mod 7)
    let diff = (a.clone() - b.clone())
        .unwrap();

    assert_eq!(
        diff.value,
        BigInt::from(5)
    );
}

#[test]

fn test_field_element_is_zero_is_one() {

    let field = FiniteField::new(7);

    let zero = FieldElement::new(
        BigInt::from(0),
        field.clone(),
    );

    let one = FieldElement::new(
        BigInt::from(1),
        field.clone(),
    );

    let three = FieldElement::new(
        BigInt::from(3),
        field.clone(),
    );

    assert!(zero.is_zero());

    assert!(!zero.is_one());

    assert!(!one.is_zero());

    assert!(one.is_one());

    assert!(!three.is_zero());

    assert!(!three.is_one());
}

#[test]

fn test_field_element_pow() {

    let field = FiniteField::new(7);

    let a = FieldElement::new(
        BigInt::from(3),
        field.clone(),
    );

    // 3^0 = 1
    let pow0 = a.pow(0);

    assert!(pow0.is_one());

    // 3^1 = 3
    let pow1 = a.pow(1);

    assert_eq!(
        pow1.value,
        BigInt::from(3)
    );

    // 3^2 = 9 = 2 (mod 7)
    let pow2 = a.pow(2);

    assert_eq!(
        pow2.value,
        BigInt::from(2)
    );

    // 3^6 = 729 = 1 (mod 7) by Fermat's little theorem
    let pow6 = a.pow(6);

    assert!(pow6.is_one());
}

#[test]

fn test_field_element_inverse() {

    let field = FiniteField::new(7);

    let a = FieldElement::new(
        BigInt::from(3),
        field.clone(),
    );

    let inv_a = a.inverse().expect(
        "3 should be invertible in \
         GF(7)",
    );

    // 3 * 5 = 15 = 1 (mod 7), so inv(3) = 5
    assert_eq!(
        inv_a.value,
        BigInt::from(5)
    );
}

#[test]

fn test_poly_operations_gf() {

    let field = FiniteField::new(7);

    // p1 = x + 1 = [1, 1]
    let p1 = Expr::Polynomial(vec![
        Expr::BigInt(BigInt::from(1)),
        Expr::BigInt(BigInt::from(1)),
    ]);

    // p2 = x + 2 = [1, 2]
    let p2 = Expr::Polynomial(vec![
        Expr::BigInt(BigInt::from(1)),
        Expr::BigInt(BigInt::from(2)),
    ]);

    // Add: (x+1) + (x+2) = 2x + 3
    let sum =
        poly_add_gf(&p1, &p2, &field)
            .unwrap();

    if let Expr::Polynomial(coeffs) =
        sum
    {

        assert_eq!(coeffs.len(), 2);

        // coeffs[0] = 2, coeffs[1] = 3
        assert_eq!(
            coeffs[0],
            Expr::BigInt(BigInt::from(
                2
            ))
        );

        assert_eq!(
            coeffs[1],
            Expr::BigInt(BigInt::from(
                3
            ))
        );
    } else {

        panic!("Expected Polynomial");
    }
}
