use rssn::Expr;
use rssn::definite_integrate;
use rssn::output::typst::to_typst;
use rssn::parse_expr;
use rssn::prelude::numeric_evaluate_numerical;

fn main() {
    let expr = parse_expr("-x^2 + 1").unwrap().1;

    let result = definite_integrate(
        &expr,
        "x",
        &Expr::new_constant(-1.),
        &Expr::new_constant(1.),
    );

    println!(
        "{} = {:?}",
        to_typst(&result),
        numeric_evaluate_numerical(&result)
    );
}
