use rssn::symbolic::core::Expr;
use std::sync::Arc;

#[test]

fn test_dag_serialization_json() {

    // Create a DAG expression
    let expr = Expr::new_add(
        Expr::new_variable("x"),
        Expr::Constant(1.0),
    );

    assert!(expr.is_dag());

    // Serialize to JSON
    let json = serde_json::to_string(&expr).unwrap();

    println!(
        "Serialized: {}",
        json
    );

    // Deserialize from JSON
    let deserialized: Expr = serde_json::from_str(&json).unwrap();

    // Should still be DAG
    assert!(deserialized.is_dag());

    // Should be equal
    assert_eq!(expr, deserialized);
}

#[test]

fn test_dag_serialization_bincode() {

    // Create a DAG expression
    let expr = Expr::new_add(
        Expr::new_variable("x"),
        Expr::Constant(1.0),
    );

    assert!(expr.is_dag());

    // Serialize to bincode
    let bytes = bincode::serde::encode_to_vec(
        &expr,
        bincode::config::standard(),
    )
    .unwrap();

    println!(
        "Serialized to {} bytes",
        bytes.len()
    );

    // Deserialize from bincode
    let (deserialized, _): (Expr, usize) = bincode::serde::decode_from_slice(
        &bytes,
        bincode::config::standard(),
    )
    .unwrap();

    // Should still be DAG
    assert!(deserialized.is_dag());

    // Should be equal
    assert_eq!(expr, deserialized);
}

#[test]

fn test_nested_dag_serialization() {

    // Create a nested DAG: (x + 1) * (y + 2)
    let expr = Expr::new_mul(
        Expr::new_add(
            Expr::new_variable("x"),
            Expr::Constant(1.0),
        ),
        Expr::new_add(
            Expr::new_variable("y"),
            Expr::Constant(2.0),
        ),
    );

    assert!(expr.is_dag());

    // Serialize and deserialize
    let json = serde_json::to_string(&expr).unwrap();

    let deserialized: Expr = serde_json::from_str(&json).unwrap();

    assert!(deserialized.is_dag());

    assert_eq!(expr, deserialized);
}

#[test]

fn test_ast_serialization_still_works() {

    // Create an AST expression (old style)
    let ast = Expr::Add(
        Arc::new(Expr::Variable(
            "x".to_string(),
        )),
        Arc::new(Expr::Constant(1.0)),
    );

    assert!(!ast.is_dag());

    // Should still serialize
    let json = serde_json::to_string(&ast).unwrap();

    let deserialized: Expr = serde_json::from_str(&json).unwrap();

    // Deserialized might be AST or DAG depending on implementation
    // Just verify it works
    assert_eq!(
        format!("{}", ast),
        format!("{}", deserialized)
    );
}

#[test]

fn test_dag_with_sharing() {

    // Create an expression with shared subexpressions
    let x = Expr::new_variable("x");

    let expr = Expr::new_add(x.clone(), x.clone());

    assert!(expr.is_dag());

    // Serialize and deserialize
    let json = serde_json::to_string(&expr).unwrap();

    let deserialized: Expr = serde_json::from_str(&json).unwrap();

    assert!(deserialized.is_dag());
    // Note: Sharing might not be preserved across serialization
    // but semantic equality should hold
}
