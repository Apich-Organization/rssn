use rssn::symbolic::calculus::*;
use rssn::symbolic::core::Expr;

fn main() {

    let x = Expr::new_variable("x");

    // Differentiation
    let f = Expr::new_sin(x.clone());

    println!("f(x) = {}", f);

    let df = differentiate(&f, "x");

    println!("df/dx = {}", df);

    // Integration
    let g = Expr::new_pow(
        x.clone(),
        Expr::Constant(2.0),
    );

    println!("g(x) = {}", g);

    let int_g = integrate(&g, "x", None, None);

    println!(
        "int(g) dx = {}",
        int_g
    );

    // Analytic check
    let z = Expr::new_variable("z");

    let h = Expr::new_exp(z.clone());

    println!("h(z) = {}", h);

    let is_analytic = check_analytic(&h, "z");

    println!(
        "Is h(z) analytic? {}",
        is_analytic
    );
}
