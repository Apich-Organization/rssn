use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;
use rssn::numerical::ode::*;
use rssn::symbolic::core::Expr;

#[test]

fn test_solve_ode_rk4_exponential() {

    let y0 = Expr::new_variable("y0"); // dy/dx = y
    let y_init = vec![1.0];

    let res = solve_ode_system_rk4(
        &[y0],
        &y_init,
        (0.0, 1.0),
        100,
    )
    .unwrap();

    // y(1) = e^1 approx 2.71828
    assert_approx_eq!(
        res.last().unwrap()[0],
        2.71828,
        1e-5
    );
}

#[test]

fn test_solve_ode_euler_exponential() {

    let y0 = Expr::new_variable("y0");

    let y_init = vec![1.0];

    let res = solve_ode_euler(
        &[y0],
        &y_init,
        (0.0, 1.0),
        1000,
    )
    .unwrap();

    // Euler is less accurate, needs more steps
    assert_approx_eq!(
        res.last().unwrap()[0],
        2.71828,
        1e-2
    );
}

#[test]

fn test_solve_ode_heun_exponential() {

    let y0 = Expr::new_variable("y0");

    let y_init = vec![1.0];

    let res = solve_ode_heun(
        &[y0],
        &y_init,
        (0.0, 1.0),
        100,
    )
    .unwrap();

    assert_approx_eq!(
        res.last().unwrap()[0],
        2.71828,
        1e-3
    );
}

#[test]

fn test_solve_ode_system_oscillator() {

    // dy0/dx = y1
    // dy1/dx = -y0
    // Initial: y0=1, y1=0 -> y0(x)=cos(x), y1(x)=-sin(x)
    let f1 = Expr::new_variable("y1");

    let f2 = Expr::new_mul(
        Expr::new_constant(-1.0),
        Expr::new_variable("y0"),
    );

    let y_init = vec![1.0, 0.0];

    let res = solve_ode_system(
        &[f1, f2],
        &y_init,
        (
            0.0,
            std::f64::consts::PI,
        ),
        100,
        OdeSolverMethod::RungeKutta4,
    )
    .unwrap();

    // At x=pi, y0 = cos(pi) = -1, y1 = -sin(pi) = 0
    assert_approx_eq!(
        res.last().unwrap()[0],
        -1.0,
        1e-5
    );

    assert_approx_eq!(
        res.last().unwrap()[1],
        0.0,
        1e-5
    );
}

proptest! {
    #[test]
    fn proptest_ode_linear(a in -2.0..2.0f64, x_end in 0.1..2.0f64) {
        let f = Expr::new_constant(a); // dy/dx = a -> y(x) = ax + y0
        let y_init = vec![0.0];
        let res = solve_ode_system(&[f], &y_init, (0.0, x_end), 50, OdeSolverMethod::RungeKutta4).unwrap();
        assert_approx_eq!(res.last().unwrap()[0], a * x_end, 1e-9);
    }
}
