//! Grobner Basis Computation Example
//!
//! This example demonstrates how to compute the Grobner basis for a system of
//! polynomial equations using Buchberger's algorithm.
//!
//! We will compute the Grobner basis for the intersection of two circles:
//! Equation 1: x^2 + y^2 - 1 = 0
//! Equation 2: (x-1)^2 + y^2 - 1 = 0  => x^2 - 2x + y^2 = 0

use std::collections::BTreeMap;
use std::ops::Add;
use std::ops::Mul;
use std::ops::Sub;

use rssn::symbolic::core::Expr;
use rssn::symbolic::core::Monomial;
use rssn::symbolic::core::SparsePolynomial;
use rssn::symbolic::grobner::MonomialOrder;
use rssn::symbolic::grobner::buchberger;
use rssn::symbolic::polynomial::sparse_poly_to_expr;

// Helper function to create a SparsePolynomial from a single variable.
fn var_poly(
    name: &str
) -> SparsePolynomial {

    let mut terms = BTreeMap::new();

    let mut mono = BTreeMap::new();

    mono.insert(name.to_string(), 1);

    terms.insert(
        Monomial(mono),
        Expr::Constant(1.0),
    );

    SparsePolynomial {
        terms,
    }
}

// Helper function to create a SparsePolynomial from a constant.
fn const_poly(
    value: f64
) -> SparsePolynomial {

    let mut terms = BTreeMap::new();

    terms.insert(
        Monomial(BTreeMap::new()),
        Expr::Constant(value),
    );

    SparsePolynomial {
        terms,
    }
}

fn main() {

    println!(
        "=== Grobner Basis \
         Computation Example ===\n"
    );

    // Define variables and constants as SparsePolynomials
    let x = var_poly("x");

    let y = var_poly("y");

    let one = const_poly(1.0);

    let two = const_poly(2.0);

    // Equation 1: x^2 + y^2 - 1 = 0
    let f1_sparse = (x.clone()
        * x.clone())
    .add(y.clone() * y.clone())
    .sub(one.clone());

    println!(
        "Polynomial 1 (f1): {}",
        sparse_poly_to_expr(&f1_sparse)
    );

    // Equation 2: x^2 - 2x + y^2 = 0
    let f2_sparse = (x.clone()
        * x.clone())
    .sub(two.clone() * x.clone())
    .add(y.clone() * y.clone());

    println!(
        "Polynomial 2 (f2): {}\n",
        sparse_poly_to_expr(&f2_sparse)
    );

    // Create the basis for Buchberger's algorithm
    let basis =
        vec![f1_sparse, f2_sparse];

    println!(
        "Computing Grobner basis \
         using Lexicographical \
         order...\n"
    );

    // Compute the Grobner basis
    match buchberger(
        &basis,
        MonomialOrder::Lexicographical,
    ) {
        | Ok(grobner_basis) => {

            println!("Grobner Basis:");

            for (i, poly) in
                grobner_basis
                    .iter()
                    .enumerate()
            {

                println!(
                    "  g{}: {}",
                    i + 1,
                    sparse_poly_to_expr(
                        poly
                    )
                );
            }
        },
        | Err(e) => {

            eprintln!(
                "Error computing \
                 Grobner basis: {}",
                e
            );
        },
    }

    println!(
        "\n=== Example Complete ==="
    );
}
