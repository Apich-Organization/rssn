use rssn::symbolic::calculus::differentiate;
use rssn::symbolic::core::Expr;

fn main() {

    let x = Expr::new_variable("x");

    // Define the expression: sin(x) using the new method syntax
    let expr = x.sin();

    // Differentiate with respect to 'x'
    let derivative =
        differentiate(&expr, "x");

    // The result will be cos(x)
    println!(
        "The derivative of {} is: {}",
        expr, derivative
    );
}
