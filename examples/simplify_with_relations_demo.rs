//! Simplification with Side-Relations Example
//!
//! This example demonstrates how to simplify an expression given a side-relation
//! using the `simplify_with_relations` function.
//!
//! We will simplify the expression `2*x^2` given the side-relation `x^2 + y^2 - 1 = 0`.

use rssn::input::parser::parse_expr;
use rssn::symbolic::cas_foundations::simplify_with_relations;
use rssn::symbolic::core::Expr;
use rssn::symbolic::grobner::MonomialOrder;

fn main() {

    println!(
        "=== Simplification with \
         Side-Relations Example ===\n"
    );

    // Define the expression to simplify as a string
    let input_expr_str = "2*x^2";

    // Define the side-relation as a string (LHS of equation set to 0)
    let relation_str = "x^2 + y^2 - 1";

    // let aa = parse_expr(relation_str);
    // let bb = parse_expr("$#dsazZ");
    // println!("aa: {:?}", &aa);
    // println!("aa: {}", &aa.unwrap().1);
    // println!("bb: {:?}", &bb);
    // println!("bb: {}", &bb.unwrap().1);

    // Parse the strings into expressions
    let expr_to_simplify =
        match parse_expr(input_expr_str)
        {
            | Ok(("", expr)) => expr,
            | Ok((rem, _)) => {
                panic!(
                    "Unparsed input \
                     for expression: \
                     '{}'",
                    rem
                )
            },
            | Err(e) => {
                panic!(
                    "Failed to parse \
                     expression '{}': \
                     {:?}",
                    input_expr_str, e
                )
            },
        };

    let relation_expr = match parse_expr(
        relation_str,
    ) {
        | Ok(("", expr)) => expr,
        | Ok((rem, _)) => {
            panic!(
                "Unparsed input for \
                 relation: '{}'",
                rem
            )
        },
        | Err(e) => {
            panic!(
                "Failed to parse \
                 relation '{}': {:?}",
                relation_str, e
            )
        },
    };

    println!(
        "Original Expression: {}",
        expr_to_simplify
    );

    println!(
        "Side-Relation: {} = 0\n",
        relation_expr
    );

    // Define the variables and monomial order
    let vars = &["x", "y"];

    let order =
        MonomialOrder::Lexicographical;

    // Call the simplification function
    let simplified_expr =
        simplify_with_relations(
            &expr_to_simplify,
            &[relation_expr], /* Relations are passed as a slice of Expr */
            vars,
            order,
        );

    println!(
        "Simplified Expression: {}",
        simplified_expr
    );

    println!(
        "\n=== Example Complete ==="
    );
}
