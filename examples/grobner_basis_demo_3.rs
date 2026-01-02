//! Grobner Basis Computation Example
//!
//! This example demonstrates how to compute the Grobner basis for a system of
//! polynomial equations using Buchberger's algorithm.
//!
//! We will compute the Grobner basis for the intersection of two circles:
//! Equation 1: x^2 + y^2 - 1 = 0
//! Equation 2: (x-1)^2 + y^2 - 1 = 0  => x^2 - 2x + y^2 = 0

use std::collections::BTreeMap;

use rssn::symbolic::core::Expr;
use rssn::symbolic::core::Monomial;
use rssn::symbolic::core::SparsePolynomial;
use rssn::symbolic::grobner::MonomialOrder;
use rssn::symbolic::grobner::buchberger;
use rssn::symbolic::polynomial::expr_to_sparse_poly;
use rssn::symbolic::polynomial::sparse_poly_to_expr;

fn main() {

    println!(
        "=== Grobner Basis \
         Computation Example ===\n"
    );

    // Define the variables
    let x = Expr::new_variable("x");

    let y = Expr::new_variable("y");

    let one = Expr::new_constant(1.0);

    let two = Expr::new_constant(2.0);

    // Equation 1: x^2 + y^2 - 1 = 0
    let f1_expr = Expr::new_sub(
        Expr::new_add(
            Expr::new_pow(
                x.clone(),
                two.clone(),
            ),
            Expr::new_pow(
                y.clone(),
                two.clone(),
            ),
        ),
        one.clone(),
    );

    println!(
        "Polynomial 1 (f1): {}",
        f1_expr
    );

    // Equation 2: (x-1)^2 + y^2 - 1 = 0  => x^2 - 2x + y^2 = 0
    let f2_expr = Expr::new_add(
        Expr::new_sub(
            Expr::new_pow(
                x.clone(),
                two.clone(),
            ),
            Expr::new_mul(
                two.clone(),
                x.clone(),
            ),
        ),
        Expr::new_pow(
            y.clone(),
            two.clone(),
        ),
    );

    println!(
        "Polynomial 2 (f2): {}\n",
        f2_expr
    );

    // Convert expressions to SparsePolynomial
    let f1_sparse = expr_to_sparse_poly(
        &f1_expr,
        &["x", "y"],
    );

    let f2_sparse = expr_to_sparse_poly(
        &f2_expr,
        &["x", "y"],
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
