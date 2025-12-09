use num_traits::ToPrimitive;
use rssn::symbolic::core::{DagOp, Expr};
use rssn::symbolic::stats::{correlation, covariance, mean, std_dev, variance};

fn evaluate_expr(expr: &Expr) -> Option<f64> {
    match expr {
        Expr::Constant(v) => Some(*v),
        Expr::BigInt(v) => v.to_f64(),
        Expr::Rational(v) => v.to_f64(),
        Expr::Sqrt(a) => evaluate_expr(a).map(|v| v.sqrt()),
        Expr::Add(a, b) => Some(evaluate_expr(a)? + evaluate_expr(b)?),
        Expr::Sub(a, b) => Some(evaluate_expr(a)? - evaluate_expr(b)?),
        Expr::Mul(a, b) => Some(evaluate_expr(a)? * evaluate_expr(b)?),
        Expr::Div(a, b) => Some(evaluate_expr(a)? / evaluate_expr(b)?),
        Expr::Dag(node) => evaluate_dag(node),
        _ => None,
    }
}

fn evaluate_dag(node: &rssn::symbolic::core::DagNode) -> Option<f64> {
    match &node.op {
        DagOp::Constant(v) => Some(v.into_inner()),
        DagOp::BigInt(v) => v.to_f64(),
        DagOp::Rational(v) => v.to_f64(),
        DagOp::Sqrt => evaluate_dag(&node.children[0]).map(|v| v.sqrt()),
        DagOp::Add => {
            let mut sum = 0.0;
            for c in &node.children {
                sum += evaluate_dag(c)?;
            }
            Some(sum)
        }
        DagOp::Mul => {
            let mut prod = 1.0;
            for c in &node.children {
                prod *= evaluate_dag(c)?;
            }
            Some(prod)
        }
        DagOp::Sub => {
            // Binary
            if node.children.len() == 2 {
                Some(evaluate_dag(&node.children[0])? - evaluate_dag(&node.children[1])?)
            } else {
                None
            }
        }
        DagOp::Div => {
            // Binary
            if node.children.len() == 2 {
                Some(evaluate_dag(&node.children[0])? / evaluate_dag(&node.children[1])?)
            } else {
                None
            }
        }
        DagOp::Power => {
            // Binary
            if node.children.len() == 2 {
                Some(evaluate_dag(&node.children[0])?.powf(evaluate_dag(&node.children[1])?))
            } else {
                None
            }
        }
        _ => None,
    }
}

fn assert_approx_eq(expr: &Expr, expected: f64) {
    if let Some(val) = evaluate_expr(expr) {
        assert!(
            (val - expected).abs() < 1e-9,
            "Expected {}, got {} (from {:?})",
            expected,
            val,
            expr
        );
    } else {
        panic!(
            "Expected numeric result {}, got non-numeric {:?}",
            expected, expr
        );
    }
}

#[test]
fn test_mean() {
    let data = vec![
        Expr::Constant(1.0),
        Expr::Constant(2.0),
        Expr::Constant(3.0),
    ];
    let m = mean(&data);
    assert_approx_eq(&m, 2.0);
}

#[test]
fn test_variance() {
    let data = vec![
        Expr::Constant(1.0),
        Expr::Constant(2.0),
        Expr::Constant(3.0),
    ];
    let v = variance(&data);
    assert_approx_eq(&v, 2.0 / 3.0);
}

#[test]
fn test_std_dev() {
    let data = vec![
        Expr::Constant(1.0),
        Expr::Constant(2.0),
        Expr::Constant(3.0),
    ];
    let s = std_dev(&data);
    assert_approx_eq(&s, (2.0f64 / 3.0).sqrt());
}

#[test]
fn test_covariance_correlation() {
    let data1 = vec![
        Expr::Constant(1.0),
        Expr::Constant(2.0),
        Expr::Constant(3.0),
    ];
    let data2 = vec![
        Expr::Constant(1.0),
        Expr::Constant(2.0),
        Expr::Constant(3.0),
    ];

    // Perfect correlation
    let cov = covariance(&data1, &data2);
    assert_approx_eq(&cov, 2.0 / 3.0);

    let corr = correlation(&data1, &data2);
    assert_approx_eq(&corr, 1.0);
}

#[test]
fn test_symbolic_mean() {
    let x = Expr::new_variable("x");
    let y = Expr::new_variable("y");
    let data = vec![x.clone(), y.clone()];

    let m = mean(&data);
    // Just ensure it returns something complex-ish or symbolic
    println!("Symbolic mean: {}", m);
}
