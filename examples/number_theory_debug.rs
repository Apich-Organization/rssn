use std::sync::Arc;

use rssn::symbolic::core::Expr;
use rssn::symbolic::number_theory::solve_diophantine;

fn main() {

    // Debug Pythagorean
    println!("--- Pythagorean ---");

    let x = Expr::new_variable("x");

    let y = Expr::new_variable("y");

    let z = Expr::new_variable("z");

    let eq = Expr::Eq(
        Arc::new(Expr::new_add(
            Expr::new_pow(
                x.clone(),
                Expr::new_constant(2.0),
            ),
            Expr::new_pow(
                y.clone(),
                Expr::new_constant(2.0),
            ),
        )),
        Arc::new(Expr::new_pow(
            z.clone(),
            Expr::new_constant(2.0),
        )),
    );

    let vars = vec!["x", "y", "z"];

    let result = solve_diophantine(&eq, &vars);

    println!(
        "Result Pythagorean: {:?}",
        result
    );
}
