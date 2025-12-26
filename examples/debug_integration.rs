use rssn::symbolic::calculus::definite_integrate;
use rssn::symbolic::core::Expr;
use rssn::symbolic::simplify_dag::simplify;

fn main() {

    let x = Expr::Variable("x".to_string());

    let lower = Expr::Constant(-1.0);

    let upper = Expr::Constant(1.0);

    // Integrate x from -1 to 1
    let integral = definite_integrate(
        &x,
        "x",
        &lower,
        &upper,
    );

    println!(
        "Integral of x from -1 to 1: {:?}",
        integral
    );

    let simplified = simplify(&integral);

    println!(
        "Simplified: {:?}",
        simplified
    );
}
