use rssn::symbolic::core::Expr;

fn main() {
    println!("=== Symbolic Operator Overloading Demo ===");

    // 1. Basic Arithmetic
    let x = Expr::new_variable("x");
    let y = Expr::new_variable("y");
    let c2 = Expr::new_constant(2.0);

    let sum = x.clone() + y.clone();
    let diff = x.clone() - y.clone();
    let prod = x.clone() * y.clone();
    let quot = x.clone() / y.clone();
    let neg = -x.clone();

    println!("x + y = {}", sum);
    println!("x - y = {}", diff);
    println!("x * y = {}", prod);
    println!("x / y = {}", quot);
    println!("-x = {}", neg);

    // 2. Mixing with f64
    let expr_plus_float = x.clone() + 5.0;
    let float_plus_expr = 10.0 + x.clone();
    let expr_times_float = y.clone() * 3.14;

    println!("x + 5.0 = {}", expr_plus_float);
    println!("10.0 + x = {}", float_plus_expr);
    println!("y * 3.14 = {}", expr_times_float);

    // 3. Mathematical Methods
    let sin_x = x.sin();
    let cos_x = x.cos();
    let tan_x = x.tan();
    let exp_x = x.exp();
    let ln_x = x.ln();
    let sqrt_x = x.sqrt();
    let abs_x = x.abs();

    println!("sin(x) = {}", sin_x);
    println!("cos(x) = {}", cos_x);
    println!("tan(x) = {}", tan_x);
    println!("exp(x) = {}", exp_x);
    println!("ln(x) = {}", ln_x);
    println!("sqrt(x) = {}", sqrt_x);
    println!("|x| = {}", abs_x);

    // 4. Power and Inverse Trig
    let x_squared = x.pow(Expr::new_constant(2.0));
    let asin_x = x.asin();
    let acos_x = x.acos();
    let atan_x = x.atan();

    println!("x^2 = {}", x_squared);
    println!("asin(x) = {}", asin_x);
    println!("acos(x) = {}", acos_x);
    println!("atan(x) = {}", atan_x);

    // 5. Hyperbolic Functions
    let sinh_x = x.sinh();
    let cosh_x = x.cosh();
    let tanh_x = x.tanh();

    println!("sinh(x) = {}", sinh_x);
    println!("cosh(x) = {}", cosh_x);
    println!("tanh(x) = {}", tanh_x);

    // 6. Complex Expressions
    // sin(x)^2 + cos(x)^2
    let trig_identity = x.sin().pow(Expr::new_constant(2.0)) + x.cos().pow(Expr::new_constant(2.0));
    println!("sin(x)^2 + cos(x)^2 = {}", trig_identity);

    // Quadratic Formula numerator: -b + sqrt(b^2 - 4ac)
    let a = Expr::new_variable("a");
    let b = Expr::new_variable("b");
    let c = Expr::new_variable("c");

    let discriminant = b.pow(Expr::new_constant(2.0)) - Expr::new_constant(4.0) * a.clone() * c.clone();
    let numerator = -b.clone() + discriminant.sqrt();
    println!("Quadratic numerator: {}", numerator);

    println!("\n=== Demo Complete ===");
}
