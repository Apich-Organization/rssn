use std::sync::Arc;

use num_traits::ToPrimitive;
use rssn::symbolic::core::Expr;
use rssn::symbolic::numeric::evaluate_numerical;
use rssn::symbolic::solve::solve;
use rssn::symbolic::solve::solve_linear_system;
use rssn::symbolic::solve::solve_system;

fn assert_is_value(
    expr: &Expr,
    expected: f64,
) {

    let val = evaluate_numerical(expr)
        .unwrap_or_else(|| {

            panic!(
                "Failed to evaluate \
                 expression: {:?}",
                expr
            )
        });

    assert!(
        (val - expected).abs() < 1e-10,
        "Expected {}, got {}",
        expected,
        val
    );
}

#[test]

fn test_solve_linear() {

    // 2x + 4 = 0  =>  x = -2
    let eq = Expr::new_add(
        Expr::new_mul(
            Expr::new_constant(2.0),
            Expr::new_variable("x"),
        ),
        Expr::new_constant(4.0),
    );

    let solutions = solve(&eq, "x");

    assert_eq!(solutions.len(), 1);

    assert_is_value(
        &solutions[0],
        -2.0,
    );
}

#[test]

fn test_solve_quadratic() {

    // x^2 - 4 = 0  =>  x = 2, x = -2
    let eq = Expr::new_sub(
        Expr::new_pow(
            Expr::new_variable("x"),
            Expr::new_constant(2.0),
        ),
        Expr::new_constant(4.0),
    );

    let solutions = solve(&eq, "x");

    assert_eq!(solutions.len(), 2);

    // Solutions might be in any order
    let has_2 = solutions
        .iter()
        .any(|s| {

            evaluate_numerical(s)
                .map(|v| {

                    (v - 2.0).abs()
                        < 1e-10
                })
                .unwrap_or(false)
        });

    let has_neg_2 = solutions
        .iter()
        .any(|s| {

            evaluate_numerical(s)
                .map(|v| {

                    (v - -2.0).abs()
                        < 1e-10
                })
                .unwrap_or(false)
        });

    assert!(
        has_2,
        "Should contain 2.0, got {:?}",
        solutions
    );

    assert!(
        has_neg_2,
        "Should contain -2.0, got {:?}",
        solutions
    );
}

#[test]

fn test_solve_linear_system() {

    // x + y = 3
    // x - y = 1
    // => x = 2, y = 1
    let eq1 = Expr::new_sub(
        Expr::new_add(
            Expr::new_variable("x"),
            Expr::new_variable("y"),
        ),
        Expr::new_constant(3.0),
    );

    let eq2 = Expr::new_sub(
        Expr::new_sub(
            Expr::new_variable("x"),
            Expr::new_variable("y"),
        ),
        Expr::new_constant(1.0),
    );

    let system =
        Expr::System(vec![eq1, eq2]);

    let vars = vec![
        "x".to_string(),
        "y".to_string(),
    ];

    let solutions =
        solve_linear_system(
            &system,
            &vars,
        )
        .unwrap();

    assert_eq!(solutions.len(), 2);

    // solutions are ordered by vars: [x, y]
    assert_is_value(&solutions[0], 2.0);

    assert_is_value(&solutions[1], 1.0);
}

#[test]

fn test_solve_system_substitution() {

    // Simple system solvable by substitution
    // y = x + 1
    // x = 2
    // => x = 2, y = 3
    let eq1 = Expr::new_sub(
        Expr::new_variable("y"),
        Expr::new_add(
            Expr::new_variable("x"),
            Expr::new_constant(1.0),
        ),
    );

    let eq2 = Expr::new_sub(
        Expr::new_variable("x"),
        Expr::new_constant(2.0),
    );

    // solve_system tries substitution first
    let solutions = solve_system(
        &[eq1, eq2],
        &["x", "y"],
    )
    .unwrap();

    // Returns Vec<(Expr, Expr)> pairs of (Variable, Value)
    let x_sol = solutions
        .iter()
        .find(|(v, _)| matches!(v, Expr::Variable(name) if name == "x"));

    let y_sol = solutions
        .iter()
        .find(|(v, _)| matches!(v, Expr::Variable(name) if name == "y"));

    assert!(x_sol.is_some());

    assert!(y_sol.is_some());

    assert_is_value(
        &x_sol.unwrap().1,
        2.0,
    );

    assert_is_value(
        &y_sol.unwrap().1,
        3.0,
    );
}

#[test]

fn test_solve_cubic() {

    // x^3 - 6x^2 + 11x - 6 = 0 => x = 1, 2, 3
    let x = Expr::new_variable("x");

    let eq = Expr::new_add(
        Expr::new_sub(
            Expr::new_add(
                Expr::new_pow(
                    x.clone(),
                    Expr::new_constant(3.0),
                ),
                Expr::new_mul(
                    Expr::new_constant(-6.0),
                    Expr::new_pow(
                        x.clone(),
                        Expr::new_constant(2.0),
                    ),
                ),
            ),
            Expr::new_mul(
                Expr::new_constant(-11.0),
                x.clone(),
            ),
        ),
        Expr::new_constant(-6.0),
    );

    // Wait, x^3 - 6x^2 + 11x - 6 = 0
    // + 11x, not -11x.
    let eq = Expr::new_add(
        Expr::new_add(
            Expr::new_add(
                Expr::new_pow(
                    x.clone(),
                    Expr::new_constant(3.0),
                ),
                Expr::new_mul(
                    Expr::new_constant(-6.0),
                    Expr::new_pow(
                        x.clone(),
                        Expr::new_constant(2.0),
                    ),
                ),
            ),
            Expr::new_mul(
                Expr::new_constant(11.0),
                x.clone(),
            ),
        ),
        Expr::new_constant(-6.0),
    );

    let solutions = solve(&eq, "x");

    assert_eq!(solutions.len(), 3);

    eprintln!(
        "Solutions: {:?}",
        solutions
    );

    // Solutions might be in any order
    let mut vals: Vec<f64> = solutions
        .iter()
        .map(|s| {

            evaluate_numerical(s)
                .expect(
                "Failed to evaluate",
            )
        })
        .collect();

    vals.sort_by(|a, b| {

        a.partial_cmp(b)
            .unwrap()
    });

    assert!(
        (vals[0] - 1.0).abs() < 1e-7
    );

    assert!(
        (vals[1] - 2.0).abs() < 1e-7
    );

    assert!(
        (vals[2] - 3.0).abs() < 1e-7
    );
}

#[test]

fn test_solve_quartic_biquadratic() {

    // x^4 - 5x^2 + 4 = 0 => x = 1, -1, 2, -2
    let x = Expr::new_variable("x");

    let eq = Expr::new_add(
        Expr::new_add(
            Expr::new_pow(
                x.clone(),
                Expr::new_constant(4.0),
            ),
            Expr::new_mul(
                Expr::new_constant(
                    -5.0,
                ),
                Expr::new_pow(
                    x,
                    Expr::new_constant(
                        2.0,
                    ),
                ),
            ),
        ),
        Expr::new_constant(4.0),
    );

    let solutions = solve(&eq, "x");

    assert_eq!(solutions.len(), 4);

    let mut vals: Vec<f64> = solutions
        .iter()
        .map(|s| {

            evaluate_numerical(s)
                .expect(
                "Failed to evaluate",
            )
        })
        .collect();

    vals.sort_by(|a, b| {

        a.partial_cmp(b)
            .unwrap()
    });

    assert!(
        (vals[0] - -2.0).abs() < 1e-7
    );

    assert!(
        (vals[1] - -1.0).abs() < 1e-7
    );

    assert!(
        (vals[2] - 1.0).abs() < 1e-7
    );

    assert!(
        (vals[3] - 2.0).abs() < 1e-7
    );
}

#[test]

fn test_solve_quartic_general() {

    // x^4 - 10x^3 + 35x^2 - 50x + 24 = 0 => x = 1, 2, 3, 4
    let x = Expr::new_variable("x");

    let eq = Expr::new_add(
        Expr::new_add(
            Expr::new_add(
                Expr::new_add(
                    Expr::new_pow(
                        x.clone(),
                        Expr::new_constant(
                            4.0,
                        ),
                    ),
                    Expr::new_mul(
                        Expr::new_constant(
                            -10.0,
                        ),
                        Expr::new_pow(
                            x.clone(),
                            Expr::new_constant(
                                3.0,
                            ),
                        ),
                    ),
                ),
                Expr::new_mul(
                    Expr::new_constant(35.0),
                    Expr::new_pow(
                        x.clone(),
                        Expr::new_constant(
                            2.0,
                        ),
                    ),
                ),
            ),
            Expr::new_mul(
                Expr::new_constant(-50.0),
                x,
            ),
        ),
        Expr::new_constant(24.0),
    );

    let solutions = solve(&eq, "x");

    assert_eq!(solutions.len(), 4);

    let mut vals: Vec<f64> = solutions
        .iter()
        .map(|s| {

            evaluate_numerical(s)
                .expect(
                "Failed to evaluate",
            )
        })
        .collect();

    vals.sort_by(|a, b| {

        a.partial_cmp(b)
            .unwrap()
    });

    assert!(
        (vals[0] - 1.0).abs() < 1e-7
    );

    assert!(
        (vals[1] - 2.0).abs() < 1e-7
    );

    assert!(
        (vals[2] - 3.0).abs() < 1e-7
    );

    assert!(
        (vals[3] - 4.0).abs() < 1e-7
    );
}
