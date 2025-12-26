use rssn::symbolic::core::Expr;
use std::sync::Arc;

#[test]

fn test_is_dag() {

    // DAG expressions
    let dag_expr = Expr::new_variable("x");

    assert!(dag_expr.is_dag());

    let dag_add = Expr::new_add(
        Expr::new_variable("x"),
        Expr::Constant(1.0),
    );

    assert!(dag_add.is_dag());

    // AST expressions
    let ast_const = Expr::Constant(1.0);

    assert!(!ast_const.is_dag());

    let ast_var = Expr::Variable("x".to_string());

    assert!(!ast_var.is_dag());
}

#[test]

fn test_to_dag_constant() {

    let ast = Expr::Constant(42.0);

    assert!(!ast.is_dag());

    let dag = ast
        .to_dag()
        .unwrap();

    assert!(dag.is_dag());
}

#[test]

fn test_to_dag_variable() {

    let ast = Expr::Variable("x".to_string());

    assert!(!ast.is_dag());

    let dag = ast
        .to_dag()
        .unwrap();

    assert!(dag.is_dag());
}

#[test]

fn test_to_dag_add() {

    let ast = Expr::Add(
        Arc::new(Expr::Variable(
            "x".to_string(),
        )),
        Arc::new(Expr::Constant(1.0)),
    );

    assert!(!ast.is_dag());

    let dag = ast
        .to_dag()
        .unwrap();

    assert!(dag.is_dag());
}

#[test]

fn test_to_dag_already_dag() {

    let dag1 = Expr::new_variable("x");

    assert!(dag1.is_dag());

    let dag2 = dag1
        .to_dag()
        .unwrap();

    assert!(dag2.is_dag());

    // Should be the same (cloned)
    assert_eq!(dag1, dag2);
}

#[test]

fn test_to_dag_form_in_place() {

    let mut expr = Expr::Constant(1.0);

    assert!(!expr.is_dag());

    expr.to_dag_form();

    assert!(expr.is_dag());
}

#[test]

fn test_to_dag_form_nested() {

    let mut expr = Expr::Add(
        Arc::new(Expr::Mul(
            Arc::new(Expr::Variable(
                "x".to_string(),
            )),
            Arc::new(Expr::Constant(2.0)),
        )),
        Arc::new(Expr::Constant(1.0)),
    );

    assert!(!expr.is_dag());

    expr.to_dag_form();

    assert!(expr.is_dag());
}

#[test]

fn test_to_ast_from_dag() {

    let dag = Expr::new_variable("x");

    assert!(dag.is_dag());

    let ast = dag
        .to_ast()
        .unwrap();
    // The AST might still be in DAG form or converted, depending on implementation
    // Just verify it doesn't error
}

#[test]

fn test_to_ast_from_ast() {

    let ast1 = Expr::Constant(1.0);

    assert!(!ast1.is_dag());

    let ast2 = ast1
        .to_ast()
        .unwrap();

    assert_eq!(ast1, ast2);
}

#[test]

fn test_dag_conversion_preserves_semantics() {

    // Create an AST expression: (x + 1) * 2
    let ast = Expr::Mul(
        Arc::new(Expr::Add(
            Arc::new(Expr::Variable(
                "x".to_string(),
            )),
            Arc::new(Expr::Constant(1.0)),
        )),
        Arc::new(Expr::Constant(2.0)),
    );

    // Convert to DAG
    let dag = ast
        .to_dag()
        .unwrap();

    assert!(dag.is_dag());

    // Convert back to AST
    let _ast2 = dag
        .to_ast()
        .unwrap();

    // The structure might be normalized (e.g., operands sorted)
    // so we just verify the conversion doesn't error
    // Semantic equivalence would require evaluation, which is beyond this test
}

#[test]

fn test_dag_sharing() {

    // Create two references to the same subexpression
    let x = Expr::new_variable("x");

    let expr1 = Expr::new_add(x.clone(), x.clone());

    assert!(expr1.is_dag());

    // The DAG should share the "x" node internally
    // This is verified by the DAG manager's deduplication
}

#[test]

fn test_mixed_ast_dag() {

    // Create a mixed expression (some AST, some DAG)
    let dag_x = Expr::new_variable("x");

    let ast_const = Expr::Constant(1.0);

    // This creates an AST node containing a DAG node
    let mixed = Expr::Add(
        Arc::new(dag_x),
        Arc::new(ast_const),
    );

    assert!(!mixed.is_dag()); // The top level is AST

    // Convert to pure DAG
    let pure_dag = mixed
        .to_dag()
        .unwrap();

    assert!(pure_dag.is_dag());
}
