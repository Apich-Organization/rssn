use std::sync::Arc;
use num_bigint::BigInt;
use ordered_float::OrderedFloat;

use super::*;
use crate::symbolic::core::Expr;

#[test]
fn test_union_find_basic() {
    let mut uf = union_find::UnionFind::new();
    let s0 = uf.make_set();
    let s1 = uf.make_set();
    let s2 = uf.make_set();

    assert_eq!(uf.find(s0), s0);
    assert_eq!(uf.find(s1), s1);
    assert_ne!(uf.find(s0), uf.find(s1));

    let (winner, changed) = uf.union(s0, s1);
    assert!(changed);
    assert_eq!(uf.find(s0), winner);
    assert_eq!(uf.find(s1), winner);

    let (_, changed_again) = uf.union(s0, s1);
    assert!(!changed_again);

    uf.union(s1, s2);
    assert_eq!(uf.find(s0), uf.find(s2));
}

#[test]
fn test_egraph_hashconsing() {
    let mut egraph = EGraph::new();
    let c1 = egraph.add_node(ENode::Constant(OrderedFloat(42.0)));
    let c2 = egraph.add_node(ENode::Constant(OrderedFloat(42.0)));
    assert_eq!(c1, c2, "Hashconsing should return identical Id for identical nodes");

    let x1 = egraph.add_node(ENode::Variable("x".to_string()));
    let x2 = egraph.add_node(ENode::Variable("x".to_string()));
    assert_eq!(x1, x2);

    let add1 = egraph.add_node(ENode::Add(x1, c1));
    let add2 = egraph.add_node(ENode::Add(x2, c2));
    assert_eq!(add1, add2, "Structural subexpression sharing must hold");
}

#[test]
fn test_egraph_congruence_closure() {
    let mut egraph = EGraph::new();
    let a = egraph.add_node(ENode::Variable("a".to_string()));
    let b = egraph.add_node(ENode::Variable("b".to_string()));

    let f_a = egraph.add_node(ENode::Sin(a));
    let f_b = egraph.add_node(ENode::Sin(b));
    assert_ne!(egraph.find(f_a), egraph.find(f_b));

    // Union a and b
    egraph.union(a, b);
    egraph.rebuild();

    // Congruence closure must deduce that sin(a) == sin(b)
    assert_eq!(egraph.find(f_a), egraph.find(f_b));
}

#[test]
fn test_constant_folding() {
    let pipeline = EGraphPipeline::builder()
        .with_algebra()
        .build();

    // 2 + 3 -> 5
    let expr = Expr::Add(
        Arc::new(Expr::BigInt(BigInt::from(2))),
        Arc::new(Expr::BigInt(BigInt::from(3))),
    );
    let simplified = pipeline.simplify(&expr);
    assert_eq!(simplified, Expr::BigInt(BigInt::from(5)));

    // 10 * 4 -> 40
    let expr2 = Expr::Mul(
        Arc::new(Expr::BigInt(BigInt::from(10))),
        Arc::new(Expr::BigInt(BigInt::from(4))),
    );
    let simplified2 = pipeline.simplify(&expr2);
    assert_eq!(simplified2, Expr::BigInt(BigInt::from(40)));
}

#[test]
fn test_identity_reductions() {
    let pipeline = EGraphPipeline::builder()
        .with_algebra()
        .build();

    let x = Arc::new(Expr::Variable("x".to_string()));
    let zero = Arc::new(Expr::BigInt(BigInt::from(0)));
    let one = Arc::new(Expr::BigInt(BigInt::from(1)));

    // x + 0 -> x
    let add_zero = Expr::Add(x.clone(), zero.clone());
    assert_eq!(pipeline.simplify(&add_zero), *x);

    // x * 1 -> x
    let mul_one = Expr::Mul(x.clone(), one.clone());
    assert_eq!(pipeline.simplify(&mul_one), *x);

    // x * 0 -> 0
    let mul_zero = Expr::Mul(x.clone(), zero.clone());
    assert_eq!(pipeline.simplify(&mul_zero), *zero);

    // x - x -> 0
    let sub_self = Expr::Sub(x.clone(), x.clone());
    assert_eq!(pipeline.simplify(&sub_self), *zero);

    // x / x -> 1
    let div_self = Expr::Div(x.clone(), x.clone());
    assert_eq!(pipeline.simplify(&div_self), *one);
}

#[test]
fn test_differentiation_de_cocooning() {
    let pipeline = EGraphPipeline::builder()
        .with_algebra()
        .with_calculus()
        .build();

    let x = Arc::new(Expr::Variable("x".to_string()));
    let y = Arc::new(Expr::Variable("y".to_string()));
    let c5 = Arc::new(Expr::BigInt(BigInt::from(5)));

    // d/dx(5) = 0
    let d_const = Expr::Derivative(c5, "x".to_string());
    assert_eq!(pipeline.simplify(&d_const), Expr::BigInt(BigInt::from(0)));

    // d/dx(x) = 1
    let d_x = Expr::Derivative(x.clone(), "x".to_string());
    assert_eq!(pipeline.simplify(&d_x), Expr::BigInt(BigInt::from(1)));

    // d/dx(y) remains d/dx(y) to support dependent functions in ODEs/PDEs/Euler-Lagrange
    let d_y = Expr::Derivative(y, "x".to_string());
    assert_eq!(pipeline.simplify(&d_y), d_y);

    // d/dx(sin(x)) = cos(x)
    let sin_x = Arc::new(Expr::Sin(x.clone()));
    let d_sin = Expr::Derivative(sin_x, "x".to_string());
    let simplified_d_sin = pipeline.simplify(&d_sin);
    assert_eq!(simplified_d_sin, Expr::Cos(x.clone()));

    // d/dx(x^2) = 2*x
    let two = Arc::new(Expr::BigInt(BigInt::from(2)));
    let x_sq = Arc::new(Expr::Power(x.clone(), two.clone()));
    let d_x_sq = Expr::Derivative(x_sq, "x".to_string());
    let simplified_x_sq = pipeline.simplify(&d_x_sq);
    assert_eq!(
        simplified_x_sq,
        Expr::Mul(two, x.clone())
    );
}

#[test]
fn test_trig_pythagorean_identity() {
    let pipeline = EGraphPipeline::builder()
        .with_algebra()
        .with_trigonometry()
        .build();

    let x = Arc::new(Expr::Variable("x".to_string()));
    let two = Arc::new(Expr::BigInt(BigInt::from(2)));

    // sin^2(x) + cos^2(x) -> 1
    let sin_x = Arc::new(Expr::Sin(x.clone()));
    let cos_x = Arc::new(Expr::Cos(x.clone()));
    let sin_sq = Arc::new(Expr::Power(sin_x, two.clone()));
    let cos_sq = Arc::new(Expr::Power(cos_x, two));
    let pythagoras = Expr::Add(sin_sq, cos_sq);

    assert_eq!(
        pipeline.simplify(&pythagoras),
        Expr::BigInt(BigInt::from(1))
    );
}

#[test]
fn test_facade_diff_function() {
    let x = Expr::Variable("x".to_string());
    let expr = Expr::Mul(
        Arc::new(Expr::BigInt(BigInt::from(3))),
        Arc::new(Expr::Sin(Arc::new(x))),
    );

    // d/dx(3 * sin(x)) = 3 * cos(x)
    let res = diff(&expr, "x");
    assert_eq!(
        res,
        Expr::Mul(
            Arc::new(Expr::BigInt(BigInt::from(3))),
            Arc::new(Expr::Cos(Arc::new(Expr::Variable("x".to_string()))))
        )
    );
}

#[test]
fn test_solve_oracle_integration() {
    let pipeline = EGraphPipeline::builder()
        .with_algebra()
        .with_equation_solver()
        .build();

    // Solve(2*x - 4, "x")
    let eq = Expr::Sub(
        Arc::new(Expr::Mul(
            Arc::new(Expr::Constant(2.0)),
            Arc::new(Expr::Variable("x".to_string())),
        )),
        Arc::new(Expr::Constant(4.0)),
    );
    let solve_expr = Expr::Solve(Arc::new(eq), "x".to_string());
    let res = pipeline.simplify(&solve_expr);

    // Should resolve to 2.0 or Solutions([2.0])
    assert!(
        matches!(res, Expr::Constant(c) if (c - 2.0).abs() < 1e-6)
            || matches!(res, Expr::Solutions(ref sols) if !sols.is_empty())
    );
}

#[test]
fn test_chained_typestate_compilation() {
    // Valid typestate compilation: Unconfigured -> HasRules -> Ready
    let custom_pipeline = EGraphPipeline::builder()
        .with_algebra()
        .with_calculus()
        .with_heuristic_budget(
            HeuristicBudget::new()
                .max_depth(5)
                .max_iterations(10)
                .prioritize_high_weight_ops(true),
        )
        .build();

    let x = Arc::new(Expr::Variable("x".to_string()));
    let sum = Expr::Add(x.clone(), x.clone());
    let res = custom_pipeline.simplify(&sum);
    assert!(res == *x || matches!(res, Expr::Mul(_, _) | Expr::Add(_, _)));
}
