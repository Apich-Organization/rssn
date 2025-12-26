use std::collections::HashMap;

use rssn::numerical::elementary::eval_expr;
use rssn::symbolic::core::Expr;
use rssn::symbolic::matrix;
use rssn::symbolic::simplify_dag::simplify;

fn main() {

    let a = Expr::Matrix(vec![
        vec![
            Expr::Constant(2.0),
            Expr::Constant(0.0),
        ],
        vec![
            Expr::Constant(0.0),
            Expr::Constant(2.0),
        ],
    ]);

    let inv = Expr::Matrix(vec![
        vec![
            Expr::Constant(0.5),
            Expr::Constant(0.0),
        ],
        vec![
            Expr::Constant(0.0),
            Expr::Constant(0.5),
        ],
    ]);

    let product = matrix::mul_matrices(&a, &inv);

    println!(
        "Product: {:?}",
        product
    );

    let simplified = simplify(&product);

    println!(
        "Simplified: {:?}",
        simplified
    );

    if let Expr::Matrix(prod_mat) = simplified {

        for i in 0 .. prod_mat.len() {

            for j in 0 .. prod_mat[i].len() {

                println!(
                    "Element [{}][{}]: {:?}",
                    i,
                    j,
                    eval_expr(
                        &prod_mat[i][j],
                        &HashMap::new()
                    )
                );
            }
        }
    }
}
