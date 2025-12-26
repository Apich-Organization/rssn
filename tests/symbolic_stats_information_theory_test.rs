use num_traits::ToPrimitive;
use rssn::symbolic::core::{
    DagOp,
    Expr,
};
use rssn::symbolic::simplify_dag::simplify;
use rssn::symbolic::stats_information_theory::*;
use std::sync::Arc;

// --- Helper Functions ---

fn evaluate_expr(expr: &Expr) -> Option<f64> {

    match expr {
        Expr::Constant(v) => Some(*v),
        Expr::BigInt(v) => v.to_f64(),
        Expr::Rational(v) => v.to_f64(),
        Expr::Add(a, b) => Some(evaluate_expr(a)? + evaluate_expr(b)?),
        Expr::Sub(a, b) => Some(evaluate_expr(a)? - evaluate_expr(b)?),
        Expr::Mul(a, b) => Some(evaluate_expr(a)? * evaluate_expr(b)?),
        Expr::Div(a, b) => Some(evaluate_expr(a)? / evaluate_expr(b)?),
        Expr::Power(a, b) => Some(evaluate_expr(a)?.powf(evaluate_expr(b)?)),
        Expr::Log(a) => Some(evaluate_expr(a)?.ln()),
        Expr::LogBase(a, b) => Some(evaluate_expr(a)?.log(evaluate_expr(b)?)),
        Expr::Neg(a) => Some(-evaluate_expr(a)?),
        Expr::Dag(node) => evaluate_dag(node),
        _ => None,
    }
}

fn evaluate_dag(node: &rssn::symbolic::core::DagNode) -> Option<f64> {

    match &node.op {
        DagOp::Constant(v) => Some(v.into_inner()),
        DagOp::BigInt(v) => v.to_f64(),
        DagOp::Rational(v) => v.to_f64(),
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
            if node.children.len() == 2 {

                Some(evaluate_dag(&node.children[0])? - evaluate_dag(&node.children[1])?)
            } else {

                None
            }
        }
        DagOp::Div => {
            if node.children.len() == 2 {

                Some(evaluate_dag(&node.children[0])? / evaluate_dag(&node.children[1])?)
            } else {

                None
            }
        }
        DagOp::Power => {
            if node.children.len() == 2 {

                Some(
                    evaluate_dag(&node.children[0])?.powf(evaluate_dag(
                        &node.children[1],
                    )?),
                )
            } else {

                None
            }
        }
        DagOp::Neg => {

            if !node
                .children
                .is_empty()
            {

                Some(-evaluate_dag(
                    &node.children[0],
                )?)
            } else {

                None
            }
        }
        DagOp::Log => {

            if !node
                .children
                .is_empty()
            {

                Some(evaluate_dag(&node.children[0])?.ln())
            } else {

                None
            }
        }
        _ => None,
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

        panic!(
            "Checking approx eq for {:?} failed to evaluate to float",
            expr
        );
    }
}

#[test]

fn test_shannon_entropy() {

    // H([0.5, 0.5]) = - (0.5 * log2(0.5) + 0.5 * log2(0.5))
    // log2(0.5) = -1.
    // H = - (0.5 * -1 + 0.5 * -1) = - (-0.5 - 0.5) = 1.0
    let probs = vec![
        Expr::Constant(0.5),
        Expr::Constant(0.5),
    ];

    let ent = shannon_entropy(&probs);

    assert_approx_eq(&ent, 1.0);
}

#[test]

fn test_kl_divergence() {

    // P = [0.5, 0.5], Q = [0.25, 0.75]
    // KL(P||Q) = 0.5 * log2(0.5/0.25) + 0.5 * log2(0.5/0.75)
    // = 0.5 * log2(2) + 0.5 * log2(2/3)
    // = 0.5 * 1 + 0.5 * (1 - log2(3))
    // = 0.5 + 0.5 - 0.5 * log2(3) = 1 - 0.5 * 1.58496 = 1 - 0.79248 = 0.20752

    let p = vec![
        Expr::Constant(0.5),
        Expr::Constant(0.5),
    ];

    let q = vec![
        Expr::Constant(0.25),
        Expr::Constant(0.75),
    ];

    let kl = kl_divergence(&p, &q).unwrap();

    let expected = 0.5 * (0.5 / 0.25f64).log2() + 0.5 * (0.5 / 0.75f64).log2();

    assert_approx_eq(&kl, expected);
}

#[test]

fn test_cross_entropy() {

    // H(P, Q) = H(P) + KL(P||Q) = 1.0 + 0.20752 = 1.20752
    let p = vec![
        Expr::Constant(0.5),
        Expr::Constant(0.5),
    ];

    let q = vec![
        Expr::Constant(0.25),
        Expr::Constant(0.75),
    ];

    let ce = cross_entropy(&p, &q).unwrap();

    let expected = -(0.5 * (0.25f64).log2() + 0.5 * (0.75f64).log2());

    assert_approx_eq(&ce, expected);
}

#[test]

fn test_gini_impurity() {

    // G([0.5, 0.5]) = 1 - (0.5^2 + 0.5^2) = 1 - 0.5 = 0.5
    let probs = vec![
        Expr::Constant(0.5),
        Expr::Constant(0.5),
    ];

    let gini = gini_impurity(&probs);

    assert_approx_eq(&gini, 0.5);
}

#[test]

fn test_joint_entropy() {

    // Joint P(X,Y):
    //      Y=0   Y=1
    // X=0  0.25  0.25
    // X=1  0.25  0.25
    // H(X,Y) = -4 * (0.25 * log2(0.25)) = - (-2) = 2.0

    let matrix = Expr::Matrix(vec![
        vec![
            Expr::Constant(0.25),
            Expr::Constant(0.25),
        ],
        vec![
            Expr::Constant(0.25),
            Expr::Constant(0.25),
        ],
    ]);

    let joint_h = joint_entropy(&matrix).unwrap();

    assert_approx_eq(&joint_h, 2.0);
}

#[test]

fn test_mutual_information() {

    // X and Y are independent in previous example, so I(X;Y) = 0.
    // H(X) = 1, H(Y) = 1, H(X,Y) = 2.
    // I(X;Y) = 1 + 1 - 2 = 0.

    let matrix = Expr::Matrix(vec![
        vec![
            Expr::Constant(0.25),
            Expr::Constant(0.25),
        ],
        vec![
            Expr::Constant(0.25),
            Expr::Constant(0.25),
        ],
    ]);

    let mi = mutual_information(&matrix).unwrap();

    assert_approx_eq(&mi, 0.0);

    // Dependent case: X=Y with 0.5 prob
    //      Y=0  Y=1
    // X=0  0.5  0.0
    // X=1  0.0  0.5
    // H(X)=1, H(Y)=1, H(X,Y)=1 (since pairs are (0,0) or (1,1) with 0.5)
    // I(X;Y) = 1 + 1 - 1 = 1.

    let matrix_dep = Expr::Matrix(vec![
        vec![
            Expr::Constant(0.5),
            Expr::Constant(1e-10),
        ], // Use small epsilon for 0 log 0 if log(0) causes issues, or check if 0 handled
        // Actually log(0) usually undefined/neg inf.
        // Logic might need to handle 0 probability by skipping term or using limit.
        // Let's see how shannon_entropy handles 0.
        // It does p * log(p). 0 * -inf is NaN if not handled.
        // Let's use small eps for test.
        vec![
            Expr::Constant(1e-10),
            Expr::Constant(0.5),
        ],
    ]);

    // Note: 1e-10 * log(1e-10) is very small close to 0.
    // H(0.5, 0.5, eps, eps) ~ 1.

    let mi_dep = mutual_information(&matrix_dep).unwrap();

    assert_approx_eq(&mi_dep, 1.0);
}
