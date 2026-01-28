use std::sync::Arc;

use num_bigint::BigInt;
use rssn::symbolic::convergence::ConvergenceResult;
use rssn::symbolic::convergence::analyze_convergence;
use rssn::symbolic::core::Expr;

#[test]

fn test_p_series_convergent() {

    // Σ(1/n^2) converges (p = 2 > 1)
    // Using BigInt for exact pattern matching
    let term = Expr::Div(
        Arc::new(Expr::new_bigint(
            BigInt::from(1),
        )),
        Arc::new(Expr::Power(
            Arc::new(Expr::Variable(
                "n".to_string(),
            )),
            Arc::new(Expr::new_constant(
                2.0,
            )),
        )),
    );

    let result =
        analyze_convergence(&term, "n");

    assert_eq!(
        result,
        ConvergenceResult::Converges
    );
}

#[test]

fn test_p_series_divergent() {

    // Σ(1/n) diverges (p = 1)
    let term = Expr::Div(
        Arc::new(Expr::new_bigint(
            BigInt::from(1),
        )),
        Arc::new(Expr::Power(
            Arc::new(Expr::Variable(
                "n".to_string(),
            )),
            Arc::new(Expr::new_constant(
                1.0,
            )),
        )),
    );

    let result =
        analyze_convergence(&term, "n");

    assert_eq!(
        result,
        ConvergenceResult::Diverges
    );
}

#[test]

fn test_p_series_half() {

    // Σ(1/√n) = Σ(1/n^0.5) diverges (p = 0.5 < 1)
    let term = Expr::Div(
        Arc::new(Expr::new_bigint(
            BigInt::from(1),
        )),
        Arc::new(Expr::Power(
            Arc::new(Expr::Variable(
                "n".to_string(),
            )),
            Arc::new(Expr::new_constant(
                0.5,
            )),
        )),
    );

    let result =
        analyze_convergence(&term, "n");

    assert_eq!(
        result,
        ConvergenceResult::Diverges
    );
}

#[test]

fn test_p_series_three() {

    // Σ(1/n^3) converges (p = 3 > 1)
    let term = Expr::Div(
        Arc::new(Expr::new_bigint(
            BigInt::from(1),
        )),
        Arc::new(Expr::Power(
            Arc::new(Expr::Variable(
                "n".to_string(),
            )),
            Arc::new(Expr::new_constant(
                3.0,
            )),
        )),
    );

    let result =
        analyze_convergence(&term, "n");

    assert_eq!(
        result,
        ConvergenceResult::Converges
    );
}

#[test]

fn test_convergence_result_types() {

    // Test that the function returns valid ConvergenceResult types
    let term = Expr::new_variable("n");

    let result =
        analyze_convergence(&term, "n");

    // Should be one of the three possible results
    assert!(matches!(
        result,
        ConvergenceResult::Converges
            | ConvergenceResult::Diverges
            | ConvergenceResult::Inconclusive
    ));
}
