use std::sync::Arc;

use rssn::symbolic::core::Expr;
use rssn::symbolic::logic::is_satisfiable;
use rssn::symbolic::logic::simplify_logic;
use rssn::symbolic::logic::to_cnf;
use rssn::symbolic::logic::to_dnf;

#[test]

fn test_simplify_double_negation_boolean(
) {

    // Not(Not(True)) -> True
    let expr = Expr::Not(Arc::new(
        Expr::Not(Arc::new(
            Expr::Boolean(true),
        )),
    ));

    let result = simplify_logic(&expr);

    assert_eq!(
        result,
        Expr::Boolean(true)
    );
}

#[test]

fn test_simplify_and_with_false() {

    // A And False -> False
    let a = Expr::Predicate {
        name : "A".to_string(),
        args : vec![],
    };

    let expr = Expr::And(vec![
        a,
        Expr::Boolean(false),
    ]);

    let result = simplify_logic(&expr);

    assert_eq!(
        result,
        Expr::Boolean(false)
    );
}

#[test]

fn test_simplify_and_with_true() {

    // A And True -> A
    let a = Expr::Predicate {
        name : "A".to_string(),
        args : vec![],
    };

    let expr = Expr::And(vec![
        a.clone(),
        Expr::Boolean(true),
    ]);

    let result = simplify_logic(&expr);

    assert_eq!(result, a);
}

#[test]

fn test_simplify_or_with_true() {

    // A Or True -> True
    let a = Expr::Predicate {
        name : "A".to_string(),
        args : vec![],
    };

    let expr = Expr::Or(vec![
        a,
        Expr::Boolean(true),
    ]);

    let result = simplify_logic(&expr);

    assert_eq!(
        result,
        Expr::Boolean(true)
    );
}

#[test]

fn test_simplify_or_with_false() {

    // A Or False -> A
    let a = Expr::Predicate {
        name : "A".to_string(),
        args : vec![],
    };

    let expr = Expr::Or(vec![
        a.clone(),
        Expr::Boolean(false),
    ]);

    let result = simplify_logic(&expr);

    assert_eq!(result, a);
}

#[test]

fn test_simplify_tautology() {

    // A Or Not(A) -> True
    let a = Expr::Predicate {
        name : "A".to_string(),
        args : vec![],
    };

    let expr = Expr::Or(vec![
        a.clone(),
        Expr::Not(Arc::new(a)),
    ]);

    let result = simplify_logic(&expr);

    assert_eq!(
        result,
        Expr::Boolean(true)
    );
}

#[test]

fn test_to_cnf_simple() {

    // (A Or B) And (C Or D) is already in CNF
    let a = Expr::Predicate {
        name : "A".to_string(),
        args : vec![],
    };

    let b = Expr::Predicate {
        name : "B".to_string(),
        args : vec![],
    };

    let c = Expr::Predicate {
        name : "C".to_string(),
        args : vec![],
    };

    let d = Expr::Predicate {
        name : "D".to_string(),
        args : vec![],
    };

    let expr = Expr::And(vec![
        Expr::Or(vec![
            a.clone(),
            b.clone(),
        ]),
        Expr::Or(vec![
            c.clone(),
            d.clone(),
        ]),
    ]);

    let result = to_cnf(&expr);

    // Should remain in CNF form or simplify
    // Just verify it doesn't panic
    let _ = result;
}

#[test]

fn test_to_cnf_distribution() {

    // A Or (B And C) -> (A Or B) And (A Or C)
    let a = Expr::Predicate {
        name : "A".to_string(),
        args : vec![],
    };

    let b = Expr::Predicate {
        name : "B".to_string(),
        args : vec![],
    };

    let c = Expr::Predicate {
        name : "C".to_string(),
        args : vec![],
    };

    let expr = Expr::Or(vec![
        a.clone(),
        Expr::And(vec![
            b.clone(),
            c.clone(),
        ]),
    ]);

    let result = to_cnf(&expr);

    // Verify it produces a result without panicking
    let _ = result;
}

#[test]

fn test_to_dnf_simple() {

    // (A And B) Or (C And D) is already in DNF
    let a = Expr::Predicate {
        name : "A".to_string(),
        args : vec![],
    };

    let b = Expr::Predicate {
        name : "B".to_string(),
        args : vec![],
    };

    let c = Expr::Predicate {
        name : "C".to_string(),
        args : vec![],
    };

    let d = Expr::Predicate {
        name : "D".to_string(),
        args : vec![],
    };

    let expr = Expr::Or(vec![
        Expr::And(vec![
            a.clone(),
            b.clone(),
        ]),
        Expr::And(vec![
            c.clone(),
            d.clone(),
        ]),
    ]);

    let result = to_dnf(&expr);

    // Verify it produces a result without panicking
    let _ = result;
}

#[test]

fn test_is_satisfiable_true() {

    // A Or Not(A) is always satisfiable (tautology)
    let a = Expr::Predicate {
        name : "A".to_string(),
        args : vec![],
    };

    let expr = Expr::Or(vec![
        a.clone(),
        Expr::Not(Arc::new(a)),
    ]);

    let result = is_satisfiable(&expr);

    assert_eq!(result, Some(true));
}

#[test]

fn test_is_satisfiable_false() {

    // A And Not(A) is unsatisfiable (contradiction)
    let a = Expr::Predicate {
        name : "A".to_string(),
        args : vec![],
    };

    let expr = Expr::And(vec![
        a.clone(),
        Expr::Not(Arc::new(a)),
    ]);

    let result = is_satisfiable(&expr);

    assert_eq!(result, Some(false));
}

#[test]

fn test_is_satisfiable_simple_sat() {

    // A Or B is satisfiable
    let a = Expr::Predicate {
        name : "A".to_string(),
        args : vec![],
    };

    let b = Expr::Predicate {
        name : "B".to_string(),
        args : vec![],
    };

    let expr = Expr::Or(vec![a, b]);

    let result = is_satisfiable(&expr);

    assert_eq!(result, Some(true));
}

#[test]

fn test_is_satisfiable_with_quantifier()
{

    // ForAll(x, P(x)) should return None (undecidable)
    let p = Expr::Predicate {
        name : "P".to_string(),
        args : vec![Expr::Variable(
            "x".to_string(),
        )],
    };

    let expr = Expr::ForAll(
        "x".to_string(),
        Arc::new(p),
    );

    let result = is_satisfiable(&expr);

    assert_eq!(result, None);
}

#[test]

fn test_simplify_implies() {

    // A => B should be converted and simplified
    let a = Expr::Predicate {
        name : "A".to_string(),
        args : vec![],
    };

    let b = Expr::Predicate {
        name : "B".to_string(),
        args : vec![],
    };

    let expr = Expr::Implies(
        Arc::new(a.clone()),
        Arc::new(b.clone()),
    );

    let result = simplify_logic(&expr);

    // Verify it produces a result without panicking
    let _ = result;
}

#[test]

fn test_simplify_equivalent() {

    // A <=> B should be converted and simplified
    let a = Expr::Predicate {
        name : "A".to_string(),
        args : vec![],
    };

    let b = Expr::Predicate {
        name : "B".to_string(),
        args : vec![],
    };

    let expr = Expr::Equivalent(
        Arc::new(a.clone()),
        Arc::new(b.clone()),
    );

    let result = simplify_logic(&expr);

    // Verify it produces a result without panicking
    let _ = result;
}

#[test]

fn test_simplify_xor() {

    // A Xor B should be converted and simplified
    let a = Expr::Predicate {
        name : "A".to_string(),
        args : vec![],
    };

    let b = Expr::Predicate {
        name : "B".to_string(),
        args : vec![],
    };

    let expr = Expr::Xor(
        Arc::new(a.clone()),
        Arc::new(b.clone()),
    );

    let result = simplify_logic(&expr);

    // Verify it produces a result without panicking
    let _ = result;
}

#[test]

fn test_de_morgan_forall() {

    // Not(ForAll(x, P(x))) -> Exists(x, Not(P(x)))
    let p = Expr::Predicate {
        name : "P".to_string(),
        args : vec![Expr::Variable(
            "x".to_string(),
        )],
    };

    let expr = Expr::Not(Arc::new(
        Expr::ForAll(
            "x".to_string(),
            Arc::new(p),
        ),
    ));

    let result = simplify_logic(&expr);

    assert!(matches!(
        result,
        Expr::Exists(_, _)
    ));
}

#[test]

fn test_de_morgan_exists() {

    // Not(Exists(x, P(x))) -> ForAll(x, Not(P(x)))
    let p = Expr::Predicate {
        name : "P".to_string(),
        args : vec![Expr::Variable(
            "x".to_string(),
        )],
    };

    let expr = Expr::Not(Arc::new(
        Expr::Exists(
            "x".to_string(),
            Arc::new(p),
        ),
    ));

    let result = simplify_logic(&expr);

    assert!(matches!(
        result,
        Expr::ForAll(_, _)
    ));
}

#[test]

fn test_quantifier_reduction() {

    // ForAll(x, True) -> True (x is not free in True)
    let expr = Expr::ForAll(
        "x".to_string(),
        Arc::new(Expr::Boolean(true)),
    );

    let result = simplify_logic(&expr);

    assert_eq!(
        result,
        Expr::Boolean(true)
    );
}

#[test]

fn test_complex_sat_problem() {

    // (A Or B) And (Not(A) Or C) And (Not(B) Or Not(C))
    // This is satisfiable with A=true, C=true, B=false
    let a = Expr::Predicate {
        name : "A".to_string(),
        args : vec![],
    };

    let b = Expr::Predicate {
        name : "B".to_string(),
        args : vec![],
    };

    let c = Expr::Predicate {
        name : "C".to_string(),
        args : vec![],
    };

    let expr = Expr::And(vec![
        Expr::Or(vec![
            a.clone(),
            b.clone(),
        ]),
        Expr::Or(vec![
            Expr::Not(Arc::new(
                a.clone(),
            )),
            c.clone(),
        ]),
        Expr::Or(vec![
            Expr::Not(Arc::new(b)),
            Expr::Not(Arc::new(c)),
        ]),
    ]);

    let result = is_satisfiable(&expr);

    // This should be satisfiable
    assert_eq!(result, Some(true));
}

#[test]

fn test_boolean_simplification() {

    // Not(True) -> False
    let expr = Expr::Not(Arc::new(
        Expr::Boolean(true),
    ));

    let result = simplify_logic(&expr);

    assert_eq!(
        result,
        Expr::Boolean(false)
    );

    // Not(False) -> True
    let expr = Expr::Not(Arc::new(
        Expr::Boolean(false),
    ));

    let result = simplify_logic(&expr);

    assert_eq!(
        result,
        Expr::Boolean(true)
    );
}

#[test]

fn test_and_flattening() {

    // And(And(A, B), C) should flatten to And(A, B, C) or similar
    let a = Expr::Predicate {
        name : "A".to_string(),
        args : vec![],
    };

    let b = Expr::Predicate {
        name : "B".to_string(),
        args : vec![],
    };

    let c = Expr::Predicate {
        name : "C".to_string(),
        args : vec![],
    };

    let expr = Expr::And(vec![
        Expr::And(vec![
            a.clone(),
            b.clone(),
        ]),
        c.clone(),
    ]);

    let result = simplify_logic(&expr);

    // Should flatten nested Ands
    if let Expr::And(terms) = result {

        assert!(terms.len() >= 2); // At least some terms
    } else {

        panic!(
            "Expected And expression, \
             got {:?}",
            result
        );
    }
}
