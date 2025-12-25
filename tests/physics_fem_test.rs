//! Tests for physics FEM module.

use rssn::physics::physics_fem::*;
use proptest::prelude::*;

#[test]
fn test_poisson_1d_basic() {
    let result = solve_poisson_1d(10, 1.0, |_| 2.0).unwrap();
    assert_eq!(result.len(), 11);
    assert_eq!(result[0], 0.0);
    assert_eq!(result[10], 0.0);
    assert!(result[5] > 0.0); // Maximum at center for -u'' = 2
}

#[test]
fn test_poisson_2d_basic() {
    let result = solve_poisson_2d(5, 5, |_, _| 2.0).unwrap();
    assert_eq!(result.len(), 36);
}

proptest! {
    #[test]
    fn prop_poisson_1d_symmetry(n in 10..50usize) {
        let result = solve_poisson_1d(n, 1.0, |_| 1.0).unwrap();
        let mid = n / 2;
        // Solution should be symmetric for symmetric force
        for i in 0..mid {
            prop_assert!((result[i] - result[n-i]).abs() < 1e-10);
        }
    }
}
