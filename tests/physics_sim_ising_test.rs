//! Unit and property-based tests for the physics sim Ising statistical module.

use proptest::prelude::*;
use rssn::physics::physics_sim::ising_statistical::*;

#[test]

fn test_ising_simulation_low_temp_smoke() {

    let params = IsingParameters {
        width : 10,
        height : 10,
        temperature : 0.1,
        mc_steps : 100,
    };

    let (grid, mag) = run_ising_simulation(&params);

    println!("Grid: {:?}", grid);

    println!(
        "Magnetization: {}",
        mag
    );

    println!(
        "grid Len: {}",
        grid.len()
    );

    assert_eq!(grid.len(), 100);

    // At very low temp, magnetization should be high (close to 1.0)
    assert!(mag > 0.5);
}

#[test]

fn test_ising_simulation_high_temp_smoke() {

    let params = IsingParameters {
        width : 10,
        height : 10,
        temperature : 10.0,
        mc_steps : 100,
    };

    let (_grid, mag) = run_ising_simulation(&params);

    // At high temp, magnetization should be low (close to 0.0)
    assert!(mag < 0.5);
}

#[test]

fn test_ising_phase_transition_scenario_run() {

    // This is a slow scenario, but we can run a reduced version or just smoke test it
    let res = simulate_ising_phase_transition_scenario();

    assert!(res.is_ok());
}

proptest! {
    #[test]
    fn prop_ising_stability(t in 0.1..5.0f64, steps in 1..10usize) {
        let params = IsingParameters {
            width: 8,
            height: 8,
            temperature: t,
            mc_steps: steps,
        };

        let (grid, mag) = run_ising_simulation(&params);
        prop_assert_eq!(grid.len(), 64);
        prop_assert!(mag >= 0.0 && mag <= 1.0);
        for &s in grid.iter() {
            prop_assert!(s == 1 || s == -1);
        }
    }
}
