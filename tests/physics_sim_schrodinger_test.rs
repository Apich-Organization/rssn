//! Unit and property-based tests for the physics sim Schrodinger quantum module.

use num_complex::Complex;
use proptest::prelude::*;
use rssn::physics::physics_sim::schrodinger_quantum::*;

#[test]

fn test_schrodinger_simulation_box_smoke(
) {

    let nx = 32;

    let ny = 32;

    let params =
        SchrodingerParameters {
            nx,
            ny,
            lx: 10.0,
            ly: 10.0,
            dt: 0.1,
            time_steps: 20,
            hbar: 1.0,
            mass: 1.0,
            potential: vec![
                0.0;
                nx * ny
            ],
        };

    let mut initial_psi =
        vec![
            Complex::new(1.0, 0.0);
            nx * ny
        ];

    let res =
        run_schrodinger_simulation(
            &params,
            &mut initial_psi,
        )
        .unwrap();

    // Check if the wave function evolves (not just zeros)
    assert!(res.len() > 0);

    let final_density =
        res.last().unwrap();

    assert_eq!(
        final_density.nrows(),
        ny
    );

    assert_eq!(
        final_density.ncols(),
        nx
    );

    let sum_density: f64 =
        final_density
            .iter()
            .sum();

    assert!(sum_density > 0.0);
}

#[test]

fn test_double_slit_scenario_run() {

    let res =
        simulate_double_slit_scenario();

    // This scenario might fail if "schrodinger_double_slit.npy" cannot be written,
    // but the logic should hold.
    assert!(res.is_ok());
}

proptest! {
    #[test]
    fn prop_schrodinger_stability(steps in 1..10usize, dt in 0.01..0.1f64) {
        let nx = 16;
        let ny = 16;
        let params = SchrodingerParameters {
            nx,
            ny,
            lx: 5.0,
            ly: 5.0,
            dt,
            time_steps: steps,
            hbar: 1.0,
            mass: 1.0,
            potential: vec![0.0; nx * ny],
        };

        let mut initial_psi = vec![Complex::new(0.5, 0.5); nx * ny];
        let res = run_schrodinger_simulation(&params, &mut initial_psi).unwrap();

        if let Some(density) = res.last() {
            for &val in density.iter() {
                prop_assert!(val.is_finite());
                prop_assert!(val >= 0.0);
            }
        }
    }
}
