use rssn::symbolic::core::Expr;
use rssn::symbolic::radicals::simplify_radicals;
use std::sync::Arc;

fn is_one(expr: &Expr) -> bool {
    match expr {
        Expr::Constant(c) => (c - 1.0).abs() < f64::EPSILON,
        Expr::Sqrt(inner) => is_one(inner),
        _ => false,
    }
}

fn is_sqrt_k(expr: &Expr, k: f64) -> bool {
    match expr {
        Expr::Sqrt(inner) => match inner.as_ref() {
            Expr::Constant(c) => (c - k).abs() < f64::EPSILON,
            _ => false,
        },
        _ => false,
    }
}

fn is_neg_sqrt_k(expr: &Expr, k: f64) -> bool {
    match expr {
        Expr::Neg(inner) => is_sqrt_k(inner, k),
        Expr::Mul(a, b) => {
            // Check for -1 * sqrt(k)
            if let Expr::Constant(c) = a.as_ref() {
                if (c + 1.0).abs() < f64::EPSILON {
                    return is_sqrt_k(b, k);
                }
            }
            if let Expr::Constant(c) = b.as_ref() {
                if (c + 1.0).abs() < f64::EPSILON {
                    return is_sqrt_k(a, k);
                }
            }
            false
        }
        _ => false,
    }
}

fn resolve(expr: Expr) -> Expr {
    if let Expr::Dag(node) = expr {
        node.to_expr().unwrap_or_else(|_| Expr::Dag(node))
    } else {
        expr
    }
}

#[test]
fn test_denest_sqrt_add() {
    // sqrt(3 + 2*sqrt(2)) = 1 + sqrt(2)
    let inner = Expr::new_add(
        Expr::new_constant(3.0),
        Expr::new_mul(
            Expr::new_constant(2.0),
            Expr::new_sqrt(Expr::new_constant(2.0))
        )
    );
    let expr = Expr::new_sqrt(inner);
    
    let simplified = resolve(simplify_radicals(&expr));
    println!("Simplified: {:?}", simplified);
    
    if let Expr::Add(a, b) = simplified {
        let has_one = is_one(&a) || is_one(&b);
        let has_sqrt2 = is_sqrt_k(&a, 2.0) || is_sqrt_k(&b, 2.0);
        assert!(has_one && has_sqrt2, "Expected 1 + sqrt(2)");
    } else {
        panic!("Expected Add expression");
    }
}

#[test]
fn test_denest_sqrt_sub() {
    // sqrt(5 - 2*sqrt(6)) = sqrt(3) - sqrt(2)
    let inner = Expr::new_sub(
        Expr::new_constant(5.0),
        Expr::new_mul(
            Expr::new_constant(2.0),
            Expr::new_sqrt(Expr::new_constant(6.0))
        )
    );
    let expr = Expr::new_sqrt(inner);
    
    let simplified = resolve(simplify_radicals(&expr));
    println!("Simplified: {:?}", simplified);
    
    // Can be Sub(sqrt(3), sqrt(2)) or Add(sqrt(3), -sqrt(2))
    match simplified {
        Expr::Sub(a, b) => {
            let is_sqrt3 = is_sqrt_k(&a, 3.0);
            let is_sqrt2 = is_sqrt_k(&b, 2.0);
            assert!(is_sqrt3 && is_sqrt2, "Expected sqrt(3) - sqrt(2)");
        },
        Expr::Add(a, b) => {
            let has_sqrt3 = is_sqrt_k(&a, 3.0) || is_sqrt_k(&b, 3.0);
            let has_neg_sqrt2 = is_neg_sqrt_k(&a, 2.0) || is_neg_sqrt_k(&b, 2.0);
            assert!(has_sqrt3 && has_neg_sqrt2, "Expected sqrt(3) - sqrt(2)");
        },
        _ => panic!("Expected Sub or Add expression"),
    }
}

#[test]
fn test_denest_sqrt_recursive() {
    // sqrt(sqrt(3 + 2*sqrt(2))) = sqrt(1 + sqrt(2))
    let inner_inner = Expr::new_add(
        Expr::new_constant(3.0),
        Expr::new_mul(
            Expr::new_constant(2.0),
            Expr::new_sqrt(Expr::new_constant(2.0))
        )
    );
    let inner = Expr::new_sqrt(inner_inner);
    let expr = Expr::new_sqrt(inner);
    
    let simplified = resolve(simplify_radicals(&expr));
    println!("Simplified recursive: {:?}", simplified);
    
    if let Expr::Sqrt(s) = simplified {
        let inner_s = resolve(s.as_ref().clone());
        if let Expr::Add(a, b) = inner_s {
             let has_one = is_one(&a) || is_one(&b);
             let has_sqrt2 = is_sqrt_k(&a, 2.0) || is_sqrt_k(&b, 2.0);
            assert!(has_one && has_sqrt2, "Expected sqrt(1 + sqrt(2))");
        } else {
            panic!("Expected Add inside Sqrt");
        }
    } else {
        panic!("Expected Sqrt expression");
    }
}

#[test]
fn test_denest_sqrt_simple() {
    // sqrt(3 + sqrt(8)) = 1 + sqrt(2)
    let inner = Expr::new_add(
        Expr::new_constant(3.0),
        Expr::new_sqrt(Expr::new_constant(8.0))
    );
    let expr = Expr::new_sqrt(inner);
    
    let simplified = resolve(simplify_radicals(&expr));
    println!("Simplified simple: {:?}", simplified);
    
    if let Expr::Add(a, b) = simplified {
        let has_one = is_one(&a) || is_one(&b);
        let has_sqrt2 = is_sqrt_k(&a, 2.0) || is_sqrt_k(&b, 2.0);
        assert!(has_one && has_sqrt2, "Expected 1 + sqrt(2)");
    } else {
        panic!("Expected Add expression");
    }
}
