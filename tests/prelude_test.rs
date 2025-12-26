#[allow(unused_imports)]
use rssn::prelude::*;

#[test]

fn test_prelude_imports() {

    // Verify that we can use some common types and functions from prelude
    let _expr = Expr::new_variable("x");

    let _mat = DMatrix::<f64>::identity(3, 3);

    let _vec = DVector::<f64>::zeros(3);

    // Verify constant
    let _date = get_build_date();

    // Verify numerical
    // let _pi = numerical::numerical_pi(); // Wait, is numerical_pi exported?
    // Let's check what's exported in numerical.
    // numerical::numerical_sin(1.0);
}

#[test]

fn test_numerical_prelude() {

    let val = numerical::numerical_sin(0.0);

    assert_eq!(val, 0.0);
}
