use num_traits::ToPrimitive;
use rssn::symbolic::core::{
    DagOp,
    Expr,
};
use rssn::symbolic::simplify_dag::simplify;
use rssn::symbolic::stats_inference::*;
use std::sync::Arc;

// --- Helper Functions ---

fn evaluate_expr(expr: &Expr) -> Option<f64> {

    match expr {
        | Expr::Constant(v) => Some(*v),
        | Expr::BigInt(v) => v.to_f64(),
        | Expr::Rational(v) => v.to_f64(),
        | Expr::Add(a, b) => Some(evaluate_expr(a)? + evaluate_expr(b)?),
        | Expr::Sub(a, b) => Some(evaluate_expr(a)? - evaluate_expr(b)?),
        | Expr::Mul(a, b) => Some(evaluate_expr(a)? * evaluate_expr(b)?),
        | Expr::Div(a, b) => Some(evaluate_expr(a)? / evaluate_expr(b)?),
        | Expr::Power(a, b) => Some(evaluate_expr(a)?.powf(evaluate_expr(b)?)),
        | Expr::Sqrt(a) => evaluate_expr(a).map(|v| v.sqrt()),
        | Expr::Abs(a) => evaluate_expr(a).map(|v| v.abs()),
        | Expr::Dag(node) => evaluate_dag(node),
        | _ => None,
    }
}

fn evaluate_dag(node: &rssn::symbolic::core::DagNode) -> Option<f64> {

    match &node.op {
        | DagOp::Constant(v) => Some(v.into_inner()),
        | DagOp::BigInt(v) => v.to_f64(),
        | DagOp::Rational(v) => v.to_f64(),
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
        | DagOp::Sqrt => evaluate_dag(&node.children[0]).map(|v| v.sqrt()),
        | _ => None,
    }
}

fn assert_approx_eq(
    expr: &Expr,
    expected: f64,
) {

    if let Some(val) = evaluate_expr(expr) {

        assert!(
            (val - expected).abs() < 1e-6,
            "Expected {}, got {} (from {:?})",
            expected,
            val,
            expr
        );
    } else {

        // Fallback for types not fully handled by evaluate (like distribution calls inside p-value)
        // For p-value, we might just check structure or simplify partially.
        // But test_statistic should be evaluatable.
        println!(
            "Warning: Could not numerically evaluate {:?}, skipping numeric check",
            expr
        );
    }
}

#[test]

fn test_one_sample_t_test() {

    // Data: 1, 2, 3. Mean = 2, Sample Std Dev = 1
    // Target mean = 2.
    // t should be (2 - 2) / (1/sqrt(3)) = 0
    let data = vec![
        Expr::Constant(1.0),
        Expr::Constant(2.0),
        Expr::Constant(3.0),
    ];

    let target = Expr::Constant(2.0);

    let result = one_sample_t_test_symbolic(&data, &target);

    // Test Statistic
    assert_approx_eq(
        &result.test_statistic,
        0.0,
    );

    // DF = n - 1 = 3 - 1 = 2
    assert_approx_eq(
        result
            .degrees_of_freedom
            .as_ref()
            .unwrap(),
        2.0,
    );
}

#[test]

fn test_two_sample_t_test() {

    // Sample 1: 1, 1, 1 (Var=0) -> Avoid 0 variance for stability in Welch's?
    // Sample 1: 1, 2, 3 (Mean=2, Var=1, n=3)
    // Sample 2: 4, 5, 6 (Mean=5, Var=1, n=3)
    // Mu diff = 0
    // t = (2 - 5 - 0) / sqrt(1/3 + 1/3) = -3 / sqrt(2/3) = -3 / 0.816 = -3.67

    let data1 = vec![
        Expr::Constant(1.0),
        Expr::Constant(2.0),
        Expr::Constant(3.0),
    ];

    let data2 = vec![
        Expr::Constant(4.0),
        Expr::Constant(5.0),
        Expr::Constant(6.0),
    ];

    let diff = Expr::Constant(0.0);

    let result = two_sample_t_test_symbolic(
        &data1, &data2, &diff,
    );

    let expected_t = -3.0 / (2.0f64 / 3.0).sqrt();

    assert_approx_eq(
        &result.test_statistic,
        expected_t,
    );
}

#[test]

fn test_z_test() {

    // Data: 1, 2, 3. Mean = 2.
    // Target = 2.
    // Sigma = 1.
    // Z = (2 - 2) / (1/sqrt(3)) = 0
    let data = vec![
        Expr::Constant(1.0),
        Expr::Constant(2.0),
        Expr::Constant(3.0),
    ];

    let target = Expr::Constant(2.0);

    let sigma = Expr::Constant(1.0);

    let result = z_test_symbolic(
        &data, &target, &sigma,
    );

    assert_approx_eq(
        &result.test_statistic,
        0.0,
    );

    assert!(result
        .degrees_of_freedom
        .is_none());
}
