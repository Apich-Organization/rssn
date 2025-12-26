use rssn::symbolic::core::{Expr, Monomial, SparsePolynomial};
use rssn::symbolic::grobner::{buchberger, poly_division_multivariate, MonomialOrder};
use std::collections::BTreeMap;

fn create_monomial(vars: &[(&str, u32)]) -> Monomial {

    let mut map = BTreeMap::new();

    for (var, exp) in vars {

        map.insert(var.to_string(), *exp);
    }

    Monomial(map)
}

fn create_sparse_poly(terms: &[(Vec<(&str, u32)>, f64)]) -> SparsePolynomial {

    let mut poly_terms = BTreeMap::new();

    for (vars, coeff) in terms {

        let mono = create_monomial(vars);

        poly_terms.insert(mono, Expr::new_constant(*coeff));
    }

    SparsePolynomial { terms: poly_terms }
}

#[test]

fn test_buchberger_simple() {

    // Test with a simple ideal: <x^2 - 1, xy - 1>
    let poly1 = create_sparse_poly(&[
        (vec![("x", 2)], 1.0),
        (vec![], -1.0),
    ]);

    let poly2 = create_sparse_poly(&[
        (vec![("x", 1), ("y", 1)], 1.0),
        (vec![], -1.0),
    ]);

    let basis = vec![poly1, poly2];

    let result = buchberger(&basis, MonomialOrder::Lexicographical);

    assert!(result.is_ok());

    let grobner = result.unwrap();

    assert!(!grobner.is_empty());
}

#[test]

fn test_poly_division_simple() {

    // Divide x^2 by x
    let dividend = create_sparse_poly(&[(vec![("x", 2)], 1.0)]);

    let divisor = create_sparse_poly(&[(vec![("x", 1)], 1.0)]);

    let result = poly_division_multivariate(&dividend, &[divisor], MonomialOrder::Lexicographical);

    assert!(result.is_ok());

    let (quotients, remainder) = result.unwrap();

    assert_eq!(quotients.len(), 1);

    assert!(remainder.terms.is_empty());
}

#[test]

fn test_poly_division_with_remainder() {

    // Divide x^2 + 1 by x
    let dividend = create_sparse_poly(&[
        (vec![("x", 2)], 1.0),
        (vec![], 1.0),
    ]);

    let divisor = create_sparse_poly(&[(vec![("x", 1)], 1.0)]);

    let result = poly_division_multivariate(&dividend, &[divisor], MonomialOrder::Lexicographical);

    assert!(result.is_ok());

    let (quotients, remainder) = result.unwrap();

    assert_eq!(quotients.len(), 1);

    assert!(!remainder.terms.is_empty());
}

#[test]

fn test_buchberger_empty() {

    let basis: Vec<SparsePolynomial> = vec![];

    let result = buchberger(&basis, MonomialOrder::Lexicographical);

    assert!(result.is_ok());

    assert!(result.unwrap().is_empty());
}

#[test]

fn test_monomial_order() {

    // Test that different monomial orders work
    let poly = create_sparse_poly(&[
        (vec![("x", 2), ("y", 1)], 1.0),
        (vec![("x", 1)], -1.0),
    ]);

    let basis = vec![poly];

    let lex = buchberger(&basis, MonomialOrder::Lexicographical);

    let grlex = buchberger(&basis, MonomialOrder::GradedLexicographical);

    let grevlex = buchberger(&basis, MonomialOrder::GradedReverseLexicographical);

    assert!(lex.is_ok());

    assert!(grlex.is_ok());

    assert!(grevlex.is_ok());
}
