use rssn::symbolic::simplify_dag::simplify;
use rssn::symbolic::core::Expr;
use num_traits::Zero;

#[test]
fn test_simplify_dag_consistency() {
    let x = Expr::new_variable("x");
    let y = Expr::new_variable("y");
    
    // (x + y) - x -> y
    let expr1 = Expr::new_sub(Expr::new_add(x.clone(), y.clone()), x.clone());
    let simplified1 = simplify(&expr1);
    assert_eq!(simplified1, y);
    
    // (2*x) - (2*x) -> 0
    let two_x = Expr::new_mul(Expr::Constant(2.0), x.clone());
    let expr2 = Expr::new_sub(two_x.clone(), two_x.clone());
    let simplified2 = simplify(&expr2);
    // Check if it simplifies to zero
    println!("Simplified (2x - 2x): {:?}", simplified2);
    
    let is_zero = |expr: &Expr| -> bool {
        match expr {
            Expr::Constant(c) => c.abs() < 1e-9,
            Expr::BigInt(n) => n.is_zero(),
            Expr::Dag(node) => {
                match &node.op {
                    rssn::symbolic::core::DagOp::Constant(c) => c.abs() < 1e-9,
                    rssn::symbolic::core::DagOp::BigInt(n) => n.is_zero(),
                    _ => false,
                }
            }
            _ => false,
        }
    };
    
    assert!(is_zero(&simplified2), "Expected 0, got {:?}", simplified2);
}
