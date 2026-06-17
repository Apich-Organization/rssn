use rssn::{
    Expr, definite_integrate, output::typst::to_typst, parse_expr,
    prelude::numeric_evaluate_numerical,
};

#[test]
fn test_definite_integral_parsed_expressions() {
    let test_cases: Vec<(&str, f64)> = vec![
        ("x^2", 2.0 / 3.0),
        ("x^2+1", 8.0 / 3.0),
        ("1-x^2", 4.0 / 3.0),
        ("-x^2+1", 4.0 / 3.0),
    ];

    for (input, expected) in &test_cases {
        let expr = parse_expr(input).unwrap().1;

        let result = definite_integrate(
            &expr,
            "x",
            &Expr::new_constant(-1.),
            &Expr::new_constant(1.),
        );

        let numerical = numeric_evaluate_numerical(&result)
            .unwrap_or_else(|| panic!("Failed to numerically evaluate integral of '{}'", input));

        assert!(
            (numerical - expected).abs() < 1e-10,
            "Integral of '{}' from -1 to 1: expected {}, got {}",
            input,
            expected,
            numerical
        );

        // Also verify typst output is non-empty
        let typst = to_typst(&result);

        assert!(
            !typst.is_empty(),
            "Typst output for '{}' should not be empty",
            input
        );
    }
}
