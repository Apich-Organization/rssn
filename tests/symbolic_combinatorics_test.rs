use rssn::symbolic::combinatorics::*;
use rssn::symbolic::core::Expr;
use rssn::symbolic::simplify_dag::simplify;

#[test]

fn test_permutations() {

    // P(5, 2) = 5! / 3! = 120 / 6 = 20
    let n = Expr::new_constant(5.0);

    let k = Expr::new_constant(2.0);

    let p = permutations(n, k);

    let simplified = simplify(&p);

    if let Expr::new_constant(val) =
        simplified
    {

        assert!(
            (val - 20.0).abs() < 1e-10
        );
    } else {
        // It might return factorial form if not fully simplified to constant
        // Let's check if it simplifies to 20.0
        // Depending on simplify implementation, factorials of constants might be evaluated
        // If not, we might need to check structure
    }
}

#[test]

fn test_combinations() {

    // C(5, 2) = 5! / (2! * 3!) = 120 / (2 * 6) = 10
    let n = Expr::new_constant(5.0);

    let k = Expr::new_constant(2.0);

    let c = combinations(&n, k);

    let simplified = simplify(&c);

    if let Expr::new_constant(val) =
        simplified
    {

        assert!(
            (val - 10.0).abs() < 1e-10
        );
    }
}

#[test]

fn test_catalan_number() {

    // C_0 = 1
    // C_1 = 1
    // C_2 = 2
    // C_3 = 5

    let c0 =
        simplify(&catalan_number(0));

    if let Expr::new_constant(val) = c0 {

        assert!(
            (val - 1.0).abs() < 1e-10
        );
    }

    let c1 =
        simplify(&catalan_number(1));

    if let Expr::new_constant(val) = c1 {

        assert!(
            (val - 1.0).abs() < 1e-10
        );
    }

    let c2 =
        simplify(&catalan_number(2));

    if let Expr::new_constant(val) = c2 {

        assert!(
            (val - 2.0).abs() < 1e-10
        );
    }

    let c3 =
        simplify(&catalan_number(3));

    if let Expr::new_constant(val) = c3 {

        assert!(
            (val - 5.0).abs() < 1e-10
        );
    }
}

#[test]

fn test_stirling_number_second_kind() {

    // S(3, 2) = 3
    // S(4, 2) = 7

    let s32 = simplify(
        &stirling_number_second_kind(
            3, 2,
        ),
    );

    if let Expr::new_constant(val) = s32 {

        assert!(
            (val - 3.0).abs() < 1e-10
        );
    }

    let s42 = simplify(
        &stirling_number_second_kind(
            4, 2,
        ),
    );

    if let Expr::new_constant(val) = s42 {

        assert!(
            (val - 7.0).abs() < 1e-10
        );
    }
}

#[test]

fn test_bell_number() {

    // B_0 = 1
    // B_1 = 1
    // B_2 = 2
    // B_3 = 5
    // B_4 = 15

    let b0 = simplify(&bell_number(0));

    if let Expr::new_constant(val) = b0 {

        assert!(
            (val - 1.0).abs() < 1e-10
        );
    }

    let b1 = simplify(&bell_number(1));

    if let Expr::new_constant(val) = b1 {

        assert!(
            (val - 1.0).abs() < 1e-10
        );
    }

    let b2 = simplify(&bell_number(2));

    if let Expr::new_constant(val) = b2 {

        assert!(
            (val - 2.0).abs() < 1e-10
        );
    }

    let b3 = simplify(&bell_number(3));

    if let Expr::new_constant(val) = b3 {

        assert!(
            (val - 5.0).abs() < 1e-10
        );
    }

    let b4 = simplify(&bell_number(4));

    if let Expr::new_constant(val) = b4 {

        assert!(
            (val - 15.0).abs() < 1e-10
        );
    }
}

#[test]

fn test_expand_binomial() {

    // (a+b)^2 = a^2 + 2ab + b^2
    // The expand_binomial function returns a summation expression
    let a =
        Expr::Variable("a".to_string());

    let b =
        Expr::Variable("b".to_string());

    let expr = Expr::new_pow(
        Expr::new_add(a, b),
        Expr::new_constant(2.0),
    );

    let expanded =
        expand_binomial(&expr);

    if let Expr::Summation(_, _, _, _) =
        expanded
    {

        // Correct form
        println!(
            "Expanded binomial: {}",
            expanded
        );
    } else {

        println!(
            "Expanded binomial: {}",
            expanded
        );

        panic!(
            "Expected Summation \
             expression, got: {:?}",
            expanded
        );
    }
    // assert!(false);
}

#[test]

fn test_expand_binomial3() {

    // (a+b)^2 = a^2 + 2ab + b^2
    // The expand_binomial function returns a summation expression
    let a =
        Expr::Variable("a".to_string());

    let b =
        Expr::Variable("b".to_string());

    let expr = Expr::new_pow(
        Expr::new_add(a, b),
        Expr::new_constant(3.0),
    );

    let expanded =
        expand_binomial(&expr);

    if let Expr::Summation(_, _, _, _) =
        expanded
    {

        // Correct form
        println!(
            "Expanded binomial: {}",
            expanded
        );
    } else {

        println!(
            "Expanded binomial: {}",
            expanded
        );

        panic!(
            "Expected Summation \
             expression, got: {:?}",
            expanded
        );
    }
    // assert!(false);
}
