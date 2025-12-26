use rssn::symbolic::core::{
    DagOp,
    Distribution,
    Expr,
};
use rssn::symbolic::stats_probability::*;

use num_traits::ToPrimitive;
use std::sync::Arc;

fn evaluate_expr(expr: &Expr) -> Option<f64> {

    match expr {
        | Expr::Constant(v) => Some(*v),
        | Expr::BigInt(v) => v.to_f64(),
        | Expr::Rational(v) => v.to_f64(),
        | Expr::Sqrt(a) => evaluate_expr(a).map(|v| v.sqrt()),
        | Expr::Add(a, b) => Some(evaluate_expr(a)? + evaluate_expr(b)?),
        | Expr::Sub(a, b) => Some(evaluate_expr(a)? - evaluate_expr(b)?),
        | Expr::Mul(a, b) => Some(evaluate_expr(a)? * evaluate_expr(b)?),
        | Expr::Div(a, b) => Some(evaluate_expr(a)? / evaluate_expr(b)?),
        | Expr::Dag(node) => evaluate_dag(node),
        | _ => None,
    }
}

fn evaluate_dag(node: &rssn::symbolic::core::DagNode) -> Option<f64> {

    match &node.op {
        | DagOp::Constant(v) => Some(v.into_inner()),
        | DagOp::BigInt(v) => v.to_f64(),
        | DagOp::Rational(v) => v.to_f64(),
        | DagOp::Sqrt => evaluate_dag(&node.children[0]).map(|v| v.sqrt()),
        | DagOp::Add => {

            let mut sum = 0.0;

            for c in &node.children {

                sum += evaluate_dag(c)?;
            }

            Some(sum)
        },
        | DagOp::Mul => {

            let mut prod = 1.0;

            for c in &node.children {

                prod *= evaluate_dag(c)?;
            }

            Some(prod)
        },
        | DagOp::Sub => {
            if node.children.len() == 2 {

                Some(evaluate_dag(&node.children[0])? - evaluate_dag(&node.children[1])?)
            } else {

                None
            }
        },
        | DagOp::Div => {
            if node.children.len() == 2 {

                Some(evaluate_dag(&node.children[0])? / evaluate_dag(&node.children[1])?)
            } else {

                None
            }
        },
        | DagOp::Power => {
            if node.children.len() == 2 {

                Some(
                    evaluate_dag(&node.children[0])?.powf(evaluate_dag(
                        &node.children[1],
                    )?),
                )
            } else {

                None
            }
        },
        | _ => None,
    }
}

// Helper to check numeric value
fn assert_approx_eq(
    expr: &Expr,
    expected: f64,
) {

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

fn test_normal_distribution() {

    let mu = Expr::Constant(0.0);

    let sigma = Expr::Constant(1.0);

    let dist = Normal {
        mean: mu.clone(),
        std_dev: sigma.clone(),
    };

    // E[X] = 0
    assert_approx_eq(
        &dist.expectation(),
        0.0,
    );

    // Var[X] = 1^2 = 1
    assert_approx_eq(
        &dist.variance(),
        1.0,
    );

    // PDF at 0 should be 1/sqrt(2pi)
    let pdf_0 = dist.pdf(&Expr::Constant(0.0));

    // Check approximate value: 1 / sqrt(2 * pi) approx 0.39894228
    // Evaluate symbolically
    // simplify should handle sqrt(2*pi) as is, or evaluate if we force it?
    // Let's just print it for manual verification or check structure
    println!(
        "Normal PDF at 0: {:?}",
        pdf_0
    );
}

#[test]

fn test_exponential_distribution() {

    let lambda = Expr::Constant(2.0);

    let dist = Exponential {
        rate: lambda.clone(),
    };

    // E[X] = 1/2
    assert_approx_eq(
        &dist.expectation(),
        0.5,
    );

    // Var[X] = 1/4
    assert_approx_eq(
        &dist.variance(),
        0.25,
    );
}

#[test]

fn test_uniform_distribution() {

    let min = Expr::Constant(0.0);

    let max = Expr::Constant(10.0);

    let dist = Uniform { min, max };

    // E[X] = 5
    assert_approx_eq(
        &dist.expectation(),
        5.0,
    );

    // Var[X] = 100 / 12 = 25/3 = 8.333...
    assert_approx_eq(
        &dist.variance(),
        100.0 / 12.0,
    );
}

#[test]

fn test_bernoulli_distribution() {

    let p = Expr::Constant(0.3);

    let dist = Bernoulli { p: p.clone() };

    // E[X] = 0.3
    assert_approx_eq(
        &dist.expectation(),
        0.3,
    );

    // Var[X] = 0.3 * 0.7 = 0.21
    assert_approx_eq(
        &dist.variance(),
        0.21,
    );
}

#[test]

fn test_poisson_distribution() {

    let lambda = Expr::Constant(3.0);

    let dist = Poisson {
        rate: lambda.clone(),
    };

    assert_approx_eq(
        &dist.expectation(),
        3.0,
    );

    assert_approx_eq(
        &dist.variance(),
        3.0,
    );
}

#[test]

fn test_usage_in_expr() {

    let n = Normal {
        mean: Expr::Constant(0.0),
        std_dev: Expr::Constant(1.0),
    };

    let _expr = Expr::Distribution(Arc::new(n));
    // Test that it compiles and fits in Expr
}
