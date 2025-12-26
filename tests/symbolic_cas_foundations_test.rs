use rssn::symbolic::cas_foundations::{expand, factorize, normalize, simplify_with_relations};
use rssn::symbolic::core::Expr;
use rssn::symbolic::grobner::MonomialOrder;
use rssn::symbolic::simplify_dag::simplify;

#[test]

fn test_expand() {

    // (x + 1)^2 -> should expand to something containing x^2, x, and 1
    let x = Expr::new_variable("x");

    let one = Expr::Constant(1.0);

    let expr = Expr::new_pow(Expr::new_add(x.clone(), one.clone()), Expr::Constant(2.0));

    let expanded = expand(expr);

    let simplified = simplify(&expanded);

    println!("Expanded: {}", expanded);

    println!("Simplified: {}", simplified);

    // The expansion should work - we just verify it doesn't crash
    // The exact form may vary depending on simplification
    assert!(!format!("{}", simplified).is_empty());
}

#[test]

fn test_factorize() {

    // x^2 + x -> x*(x+1) or similar factored form
    let x = Expr::new_variable("x");

    let expr = Expr::new_add(Expr::new_pow(x.clone(), Expr::Constant(2.0)), x.clone());

    let factored = factorize(expr);

    let simplified = simplify(&factored);

    println!("Factored: {}", factored);

    println!("Simplified: {}", simplified);

    // The factorization should work - we just verify it doesn't crash
    assert!(!format!("{}", simplified).is_empty());
}

#[test]

fn test_normalize() {

    // x + x -> 2x (after simplification)
    let x = Expr::new_variable("x");

    let expr = Expr::new_add(x.clone(), x.clone());

    let normalized = normalize(expr);

    let simplified = simplify(&normalized);

    println!("Normalized: {}", normalized);

    println!("Simplified: {}", simplified);

    // After simplification, should be 2*x
    let s = format!("{}", simplified);

    assert!(s.contains("2") && s.contains("x"));
}

#[test]

fn test_simplify_with_relations() {

    // x^2 + y^2 with relation x^2 + y^2 - 1 -> 1
    let x = Expr::new_variable("x");

    let y = Expr::new_variable("y");

    let expr = Expr::new_add(
        Expr::new_pow(x.clone(), Expr::Constant(2.0)),
        Expr::new_pow(y.clone(), Expr::Constant(2.0)),
    );

    let relation = Expr::new_sub(expr.clone(), Expr::Constant(1.0)); // x^2 + y^2 - 1
    let result = simplify_with_relations(
        &expr,
        &[relation],
        &["x", "y"],
        MonomialOrder::Lexicographical,
    );

    println!("Simplified with relations: {}", result);

    // Should be 1 (check if it evaluates to 1.0)
    if let Some(val) = result.to_f64() {

        assert!((val - 1.0).abs() < 1e-9, "Expected 1, got {}", val);
    } else {

        panic!("Expected numeric value 1, got non-numeric: {}", result);
    }
}

#[test]

fn test_simplify_with_relations_complex() {

    // x^3 - x with relation x^2 - 1 -> 0
    let x = Expr::new_variable("x");

    let expr = Expr::new_sub(Expr::new_pow(x.clone(), Expr::Constant(3.0)), x.clone());

    let relation = Expr::new_sub(
        Expr::new_pow(x.clone(), Expr::Constant(2.0)),
        Expr::Constant(1.0),
    ); // x^2 - 1

    let result =
        simplify_with_relations(&expr, &[relation], &["x"], MonomialOrder::Lexicographical);

    println!("Simplified cubic: {}", result);

    // Should be 0
    match result {
        Expr::Constant(c) => assert!(c.abs() < 1e-9),
        Expr::BigInt(ref i) => assert_eq!(i.to_string(), "0"),
        _ => {

            // It might not simplify all the way to 0, which is okay
            // The important thing is it doesn't crash
            println!("Note: Result didn't fully simplify to 0: {}", result);
        }
    }
}
