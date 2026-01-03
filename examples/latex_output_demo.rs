//! LaTeX Output Generation Example
//!
//! This example demonstrates how to render a complex symbolic expression as a LaTeX string
//! for inclusion in scientific papers and documents.

use rssn::input::parser::parse_expr;
use rssn::output::latex::to_latex;
use rssn::symbolic::core::Expr;
use std::sync::Arc;

fn main() {
    println!("=== LaTeX Output Generation Example ===\n");

    // 1. Create a complex symbolic expression.
    // Let's create the definite integral of x^2 * sin(x) from 0 to pi.
    // Expression: ∫[0,π] x²sin(x) dx

    let integrand_str = "x^2 * sin(x)";
    let integrand = match parse_expr(integrand_str) {
        Ok(("", expr)) => expr,
        Ok((rem, _)) => panic!("Unparsed input: '{}'", rem),
        Err(e) => panic!("Failed to parse '{}': {:?}", integrand_str, e),
    };

    let integral_expr = Expr::Integral {
        integrand: Arc::new(integrand),
        var: Arc::new(Expr::Variable("x".to_string())),
        lower_bound: Arc::new(Expr::Constant(0.0)),
        upper_bound: Arc::new(Expr::Pi),
    };

    println!("Symbolic Expression: \n  {}\n", integral_expr);

    // 2. Convert the expression to a LaTeX string.
    let latex_string = to_latex(&integral_expr);

    // 3. Print the LaTeX string.
    println!("Generated LaTeX String:\n  {}\n", latex_string);

    println!("To use this in a LaTeX document, you would typically place it inside a math environment, like:");
    println!("  \\begin{{equation}}");
    println!("    {}", latex_string);
    println!("  \\end{{equation}}");

    println!("\n=== Example Complete ===");
}
