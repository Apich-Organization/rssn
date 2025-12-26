//! Unit and property-based tests for the physics SM (Spectral Methods) module.

use rssn::physics::physics_sm::*;
use proptest::prelude::*;
use assert_approx_eq::assert_approx_eq;

#[test]
fn test_solve_advection_diffusion_1d_basic() {
    let result = simulate_1d_advection_diffusion_scenario();
    // Check if the Gaussian pulse has shifted or diffused
    assert_eq!(result.len(), 128);
    // Mass conservation (roughly)
    let sum: f64 = result.iter().sum();
    assert!(sum > 0.0);
}

#[test]
fn test_solve_advection_diffusion_2d_basic() {
    let result = simulate_2d_advection_diffusion_scenario();
    assert_eq!(result.len(), 64 * 64);
    let sum: f64 = result.iter().sum();
    assert!(sum > 0.0);
}

proptest! {
    #[test]
    fn prop_spectral_linearity_1d(scale in 0.1..10.0f64) {
        const N: usize = 64;
        let mut ic1 = vec![0.0; N];
        ic1[N/2] = 1.0;
        let ic2: Vec<_> = ic1.iter().map(|&x| x * scale).collect();
        
        let res1 = solve_advection_diffusion_1d(&ic1, 1.0, 1.0, 0.01, 0.01, 10);
        let res2 = solve_advection_diffusion_1d(&ic2, 1.0, 1.0, 0.01, 0.01, 10);
        
        for i in 0..N {
            prop_assert!((res1[i] * scale - res2[i]).abs() < 1e-10);
        }
    }
}
