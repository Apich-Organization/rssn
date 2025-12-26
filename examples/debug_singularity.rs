use rssn::symbolic::core::Expr;
use rssn::symbolic::simplify_dag::simplify;

fn main() {

    let z = Expr::Variable("z".to_string());

    let func = Expr::new_div(
        Expr::Constant(1.0),
        z.clone(),
    );

    println!(
        "Function: {:?}",
        func
    );

    let singularity = Expr::Constant(0.0);

    let factor = Expr::new_sub(
        z.clone(),
        singularity.clone(),
    );

    let factor_simplified = simplify(&factor);

    println!(
        "Factor (z - 0): {:?}",
        factor
    );

    println!(
        "Factor simplified: {:?}",
        factor_simplified
    );

    // Check what the denominator actually is
    if let Expr::Div(num, den) = &func {

        println!(
            "Numerator: {:?}",
            num
        );

        println!(
            "Denominator: {:?}",
            den
        );

        println!(
            "Den == z? {}",
            den.as_ref() == &z
        );

        println!(
            "Den == factor_simplified? {}",
            den.as_ref() == &factor_simplified
        );
    }
}
