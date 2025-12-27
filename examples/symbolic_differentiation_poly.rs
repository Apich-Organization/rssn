use rssn::symbolic::calculus::differentiate;
use rssn::symbolic::core::{Expr, ToConstant};

fn main() {
    let x = Expr::new_variable("x");

    // Define the expression: x^3 + 2*x
    // Using operator overloading and helper methods
    // We can now use .constant() directly on f64 values for cleaner syntax
    let expr = x.pow(3.0.constant()) + x.clone() * 2.0;

    // Differentiate with respect to 'x'
    let derivative = differentiate(&expr, "x");

    // The result should strictly be 3.0 * x^2 + 2.0 (after simplification)
    println!("Expression: {}", expr);
    println!("Derivative: {}", derivative);
}
