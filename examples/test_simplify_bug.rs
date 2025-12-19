use rssn::symbolic::core::Expr;
use rssn::symbolic::simplify_dag::simplify;

fn main() {
    let a = Expr::new_variable("a");
    let a2 = Expr::new_pow(a.clone(), Expr::Constant(2.0));
    let a3 = Expr::new_pow(a.clone(), Expr::Constant(3.0));
    let div = Expr::new_div(a2, a3);
    println!("Raw: {}", div);
    println!("Simplified: {}", simplify(&div));
}
