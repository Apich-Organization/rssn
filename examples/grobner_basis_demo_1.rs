//! Grobner Basis Computation Example
//!
//! This example demonstrates how to compute the Grobner basis for a system of
//! polynomial equations using Buchberger's algorithm.
//!
//! We will compute the Grobner basis for the intersection of two circles:
//! Equation 1: x^2 + y^2 - 1 = 0
//! Equation 2: (x-1)^2 + y^2 - 1 = 0  => x^2 - 2x + y^2 = 0

use rssn::input::parser::parse_expr;
use rssn::symbolic::core::Expr;
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

    // Define the polynomial equations as strings
    let input1 = "x^2 + y^2 - 1";

    let input2 = "x^2 - 2*x + y^2";

    // Parse the strings into expressions
    let f1_expr =
        match parse_expr(input1) {
            | Ok(("", expr)) => expr,
            | Ok((rem, _)) => {

                panic!(
                    "Unparsed input: \
                     '{}'",
                    rem
                )
            },
            | Err(e) => {

                panic!(
                    "Failed to parse \
                     expression '{}': \
                     {:?}",
                    input1, e
                )
            },
        };

    let f2_expr =
        match parse_expr(input2) {
            | Ok(("", expr)) => expr,
            | Ok((rem, _)) => {

                panic!(
                    "Unparsed input: \
                     '{}'",
                    rem
                )
            },
            | Err(e) => {

                panic!(
                    "Failed to parse \
                     expression '{}': \
                     {:?}",
                    input2, e
                )
            },
        };

    println!(
        "Polynomial 1 (f1): {}",
        f1_expr
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
