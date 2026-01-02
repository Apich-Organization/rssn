use rssn::symbolic::calculus::differentiate;
use rssn::symbolic::core::Expr;

fn main() {

    // Define a variable 'x'
    let x = Expr::new_variable("x");

    // Define the expression: sin(x)
    let expr = Expr::new_sin(x);

    // Differentiate with respect to 'x'
    let derivative =
        differentiate(&expr, "x");

    // The result will be cos(x)
    println!(
        "The derivative of {} is: {}",
        expr, derivative
    );
}
