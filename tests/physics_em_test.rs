//! Unit and property-based tests for the physics EM (Euler Methods) module.

use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;
use rssn::physics::physics_em::*;
use rssn::physics::physics_rkm::DampedOscillatorSystem;

#[test]

fn test_solve_forward_euler_decay() {

    // y' = -y, y(0) = 1 => y(t) = e^-t
    struct DecaySystem;

    impl rssn::physics::physics_rkm::OdeSystem for DecaySystem {
        fn dim(&self) -> usize {

            1
        }

        fn eval(
            &self,
            _t : f64,
            y : &[f64],
            dy : &mut [f64],
        ) {

            dy[0] = -y[0];
        }
    }

    let sys = DecaySystem;

    let res = solve_forward_euler(
        &sys,
        &[1.0],
        (0.0, 1.0),
        0.1,
    );

    let final_y = res
        .last()
        .unwrap()
        .1[0];

    // Exact: e^-1 approx 0.3678
    // Forward Euler approach: (1 - 0.1)^10 = 0.9^10 approx 0.3486
    assert_approx_eq!(
        final_y,
        0.3486784401,
        1e-10
    );
}

#[test]

fn test_solve_midpoint_euler_decay() {

    struct DecaySystem;

    impl rssn::physics::physics_rkm::OdeSystem for DecaySystem {
        fn dim(&self) -> usize {

            1
        }

        fn eval(
            &self,
            _t : f64,
            y : &[f64],
            dy : &mut [f64],
        ) {

            dy[0] = -y[0];
        }
    }

    let sys = DecaySystem;

    let res = solve_midpoint_euler(
        &sys,
        &[1.0],
        (0.0, 1.0),
        0.1,
    );

    let final_y = res
        .last()
        .unwrap()
        .1[0];

    // Midpoint method is second order, should be closer to 0.3678 than forward Euler.
    assert!(
        (final_y - 0.3678).abs()
            < (0.3486f64 - 0.3678f64)
                .abs()
    );
}

#[test]

fn test_solve_heun_euler_decay() {

    struct DecaySystem;

    impl rssn::physics::physics_rkm::OdeSystem for DecaySystem {
        fn dim(&self) -> usize {

            1
        }

        fn eval(
            &self,
            _t : f64,
            y : &[f64],
            dy : &mut [f64],
        ) {

            dy[0] = -y[0];
        }
    }

    let sys = DecaySystem;

    let res = solve_heun_euler(
        &sys,
        &[1.0],
        (0.0, 1.0),
        0.1,
    );

    let final_y = res
        .last()
        .unwrap()
        .1[0];

    // Heun's method is also second order.
    assert!(
        (final_y - 0.3678).abs()
            < (0.3486f64 - 0.3678f64)
                .abs()
    );
}

proptest! {
    #[test]
    fn prop_euler_stability_oscillator(dt in 0.001..0.05f64) {
        let system = DampedOscillatorSystem {
            omega: 1.0,
            zeta: 0.1,
        };
        let res = solve_forward_euler(&system, &[1.0, 0.0], (0.0, 1.0), dt);
        for (_, y) in res {
            prop_assert!(y[0].is_finite());
            prop_assert!(y[1].is_finite());
        }
    }
}
