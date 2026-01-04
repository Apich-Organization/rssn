use rssn::symbolic::cas_foundations::simplify_with_relations;
use rssn::symbolic::core::Expr;
use rssn::symbolic::grobner::MonomialOrder;

fn main() {

    // Define variables x and y
    let x = Expr::new_variable("x");

    let y = Expr::new_variable("y");

    // Expression to simplify: 2*x^2
    let two =
        Expr::new_bigint(2.into());

    let x_sq = Expr::new_pow(
        x.clone(),
        Expr::new_bigint(2.into()),
    );

    let expr_to_simplify =
        Expr::new_mul(two, x_sq);

    // Define the side-relation: x^2 + y^2 - 1 = 0
    let y_sq = Expr::new_pow(
        y.clone(),
        Expr::new_bigint(2.into()),
    );

    let one =
        Expr::new_bigint(1.into());

    let relation = Expr::new_sub(
        Expr::new_add(
            x.clone(),
            y.clone(),
        ),
        one,
    );

    // Simplify the expression with respect to the relation
    let simplified_expr = simplify_with_relations(
    &expr_to_simplify,
    &[relation],
    &["x", "y"],
    MonomialOrder::Lexicographical,
);

    // The result will be 2 - 4y + 2y^2
    // Note: The exact output format and canonical form may vary.
    println!(
        "Original expression: {}",
        expr_to_simplify
    );

    println!(
        "Simplified expression: {}",
        simplified_expr
    );
}
