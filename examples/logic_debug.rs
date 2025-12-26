use rssn::symbolic::core::Expr;
use rssn::symbolic::logic::{is_satisfiable, simplify_logic};
use std::sync::Arc;

fn main() {

    // Test 1: And flattening
    let a = Expr::Predicate {
        name: "A".to_string(),
        args: vec![],
    };

    let b = Expr::Predicate {
        name: "B".to_string(),
        args: vec![],
    };

    let c = Expr::Predicate {
        name: "C".to_string(),
        args: vec![],
    };

    let expr = Expr::And(vec![
        Expr::And(vec![a.clone(), b.clone()]),
        c.clone(),
    ]);

    let result = simplify_logic(&expr);

    println!("And flattening result: {:?}", result);

    // Test 2: Complex SAT
    let expr2 = Expr::And(vec![
        Expr::Or(vec![a.clone(), b.clone()]),
        Expr::Or(vec![
            Expr::Not(Arc::new(a.clone())),
            c.clone(),
        ]),
        Expr::Or(vec![
            Expr::Not(Arc::new(b.clone())),
            Expr::Not(Arc::new(c.clone())),
        ]),
    ]);

    println!("\nComplex SAT expression: {:?}", expr2);

    let result2 = is_satisfiable(&expr2);

    println!("Complex SAT result: {:?}", result2);

    // Let me try a simpler SAT problem
    let simple = Expr::And(vec![a.clone(), b.clone()]);

    println!("\nSimple SAT (A And B): {:?}", is_satisfiable(&simple));
}
