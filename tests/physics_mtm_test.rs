//! Unit and property-based tests for the physics MTM (Multigrid) module.

use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;
use rssn::physics::physics_mtm::*;

#[test]

fn test_solve_poisson_1d_multigrid_basic()
 {

    let result = simulate_1d_poisson_multigrid_scenario().unwrap();

    let n = result.len();

    let dx = 1.0 / (n - 1) as f64;

    // Check u(x) ≈ x(1-x)
    for i in 0 .. n {

        let x = i as f64 * dx;

        let expected = x * (1.0 - x);

        println!(
            "i-{} result[i]-{} \
             expected-{}",
            i, result[i], expected
        );

        assert_approx_eq!(
            result[i],
            expected,
            0.05
        );
    }
}

#[test]

fn test_solve_poisson_2d_multigrid_basic()
 {

    let result = simulate_2d_poisson_multigrid_scenario().unwrap();

    let n = (result.len() as f64).sqrt()
        as usize;

    let h = 1.0 / (n - 1) as f64;

    // Check u(x,y) ≈ sin(pi*x)*sin(pi*y)
    for i in 0 .. n {

        for j in 0 .. n {

            let x = i as f64 * h;

            let y = j as f64 * h;

            let expected = (std::f64::consts::PI * x).sin() * (std::f64::consts::PI * y).sin();

            assert_approx_eq!(
                result[i * n + j],
                expected,
                0.1
            );
        }
    }
}

proptest! {
    #[test]
    fn prop_poisson_1d_linearity(scale in 0.1..10.0f64) {
        const K: usize = 5;
        const N_INTERIOR: usize = 2_usize.pow(K as u32) - 1;
        let f1 = vec![1.0; N_INTERIOR];
        let f2 = f1.iter().map(|&x| x * scale).collect::<Vec<_>>();

        let u1 = solve_poisson_1d_multigrid(N_INTERIOR, &f1, 5).unwrap();
        let u2 = solve_poisson_1d_multigrid(N_INTERIOR, &f2, 5).unwrap();

        for i in 0..u1.len() {
            prop_assert!((u1[i] * scale - u2[i]).abs() < 1e-10);
        }
    }
}
