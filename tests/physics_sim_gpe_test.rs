//! Unit and property-based tests for the physics sim GPE superfluidity module.

use rssn::physics::physics_sim::gpe_superfluidity::*;
use proptest::prelude::*;

#[test]
fn test_simulate_bose_einstein_vortex_scenario() {
    simulate_bose_einstein_vortex_scenario();
}


#[test]
fn test_gpe_ground_state_smoke() {
    let nx = 32;
    let ny = 32;
    let params = GpeParameters {
        nx,
        ny,
        lx: 10.0,
        ly: 10.0,
        d_tau: 0.1,
        time_steps: 20,
        g: 100.0,
        trap_strength: 1.0,
    };
    
    let res = run_gpe_ground_state_finder(&params).unwrap();
    assert_eq!(res.nrows(), ny);
    assert_eq!(res.ncols(), nx);
    
    // Check if the density is normalized
    let sum_density: f64 = res.iter().sum();
    let dx = params.lx / nx as f64;
    let dy = params.ly / ny as f64;
    // The normalization in the code is: norm_factor = (n_total / (norm * dx * dy)).sqrt();
    // And n_total is nx * ny.
    // So sum(density) * dx * dy should be n_total.
    let total_prob = sum_density * dx * dy;
    assert!((total_prob - (nx * ny) as f64).abs() < 1e-5);
}

#[test]
fn test_bec_vortex_scenario_run() {
    // Reduced iterations for faster test
    let params = GpeParameters {
        nx: 64,
        ny: 64,
        lx: 10.0,
        ly: 10.0,
        d_tau: 0.05,
        time_steps: 50,
        g: 100.0,
        trap_strength: 1.0,
    };
    let res = run_gpe_ground_state_finder(&params);
    assert!(res.is_ok());
}

proptest! {
    #[test]
    fn prop_gpe_stability(steps in 1..5usize, d_tau in 0.01..0.05f64) {
        let nx = 16;
        let ny = 16;
        let params = GpeParameters {
            nx,
            ny,
            lx: 5.0,
            ly: 5.0,
            d_tau,
            time_steps: steps,
            g: 50.0,
            trap_strength: 1.0,
        };
        
        let res = run_gpe_ground_state_finder(&params).unwrap();
        for &val in res.iter() {
            prop_assert!(val.is_finite());
            prop_assert!(val >= 0.0);
        }
    }
}
