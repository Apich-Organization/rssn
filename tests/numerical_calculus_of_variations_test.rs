use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;
use rssn::numerical::calculus_of_variations::*;
use rssn::symbolic::core::Expr;

#[test]

fn test_evaluate_action_free_particle()
{

    // L = 0.5 * y_dot^2
    let t = Expr::new_variable("t");

    let y = Expr::new_variable("y");

    let y_dot =
        Expr::new_variable("y_dot");

    let lagrangian = Expr::new_mul(
        Expr::new_constant(0.5),
        Expr::new_pow(
            y_dot.clone(),
            Expr::new_constant(2.0),
        ),
    );

    // Path: y(t) = 2t -> y_dot(t) = 2
    // S = integral from 0 to 1 of 0.5 * 2^2 dt = 2.0
    let path = Expr::new_mul(
        Expr::new_constant(2.0),
        t.clone(),
    );

    let action = evaluate_action(
        &lagrangian,
        &path,
        "t",
        "y",
        "y_dot",
        (0.0, 1.0),
    )
    .unwrap();

    assert_approx_eq!(
        action,
        2.0,
        1e-5
    );
}

#[test]

fn test_evaluate_action_harmonic_oscillator(
) {

    // L = 0.5 * (y_dot^2 - y^2)
    let t = Expr::new_variable("t");

    let y = Expr::new_variable("y");

    let y_dot =
        Expr::new_variable("y_dot");

    let lagrangian = Expr::new_mul(
        Expr::new_constant(0.5),
        Expr::new_sub(
            Expr::new_pow(
                y_dot.clone(),
                Expr::new_constant(2.0),
            ),
            Expr::new_pow(
                y.clone(),
                Expr::new_constant(2.0),
            ),
        ),
    );

    // Path: y(t) = sin(t) -> y_dot(t) = cos(t)
    // S = integral from 0 to pi of 0.5 * (cos^2(t) - sin^2(t)) dt = 0.5 * integral cos(2t) dt = 0
    let path = Expr::new_sin(t.clone());

    let action = evaluate_action(
        &lagrangian,
        &path,
        "t",
        "y",
        "y_dot",
        (
            0.0,
            std::f64::consts::PI,
        ),
    )
    .unwrap();

    assert_approx_eq!(
        action,
        0.0,
        1e-5
    );
}

#[test]

fn test_euler_lagrange_free_particle() {

    // L = 0.5 * y_dot^2
    let lagrangian = Expr::new_mul(
        Expr::new_constant(0.5),
        Expr::new_pow(
            Expr::new_variable("y_dot"),
            Expr::new_constant(2.0),
        ),
    );

    let el = euler_lagrange(
        &lagrangian,
        "t",
        "y",
        "y_dot",
    );

    // d/dt (y_dot) - 0 = y_ddot
    // Our implementation returns detailed expression:
    // d/dt(dL/dy_dot)_explicit + d/dy(dL/dy_dot)*y_dot + d/dy_dot(dL/dy_dot)*y_ddot - dL/dy
    // For L = 0.5 * y_dot^2:
    // dL/dy_dot = y_dot
    // d/dt(explicit) = 0
    // d/dy(y_dot) = 0
    // d/dy_dot(y_dot) = 1
    // EL = 1 * y_dot_dot - 0 = y_dot_dot
    let s = format!("{:?}", el);

    assert!(s.contains("y_dot_dot"));
}

proptest! {
    #[test]
    fn proptest_action_linearity(c in 0.1..100.0f64) {
        // L = c * y_dot
        let t = Expr::new_variable("t");
        let y_dot = Expr::new_variable("y_dot");
        let lagrangian = Expr::new_mul(Expr::new_constant(c), y_dot);

        // Path y(t) = t^2 -> y_dot = 2t
        // S = integral 0 to 1 of c * 2t dt = c * [t^2] = c
        let path = Expr::new_pow(t.clone(), Expr::new_constant(2.0));
        let action = evaluate_action(&lagrangian, &path, "t", "y", "y_dot", (0.0, 1.0)).unwrap();
        assert_approx_eq!(action, c, 1e-5);
    }
}
