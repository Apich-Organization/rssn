//! Tests for physics FVM module.

use proptest::prelude::*;
use rssn::physics::physics_fvm::*;

#[test]

fn test_advection_1d_stability() {

    let mut mesh = Mesh::new(100, 1.0, |x| {
        if x > 0.4 && x < 0.6 {

            1.0
        } else {

            0.0
        }
    });

    let velocity = 1.0;

    let dx = 1.0 / 100.0;

    let dt = 0.5 * dx / velocity;

    let result = solve_advection_1d(
        &mut mesh,
        velocity,
        dt,
        10,
        || (0.0, 0.0),
    );

    for &val in &result {

        assert!(val.is_finite());

        assert!(val >= -1e-10); // Should be non-negative
        assert!(val <= 1.0 + 1e-10); // Should not exceed initial max
    }
}

#[test]

fn test_burgers_1d_shock() {

    let mut mesh = Mesh::new(100, 1.0, |x| {
        if x < 0.5 {

            1.0
        } else {

            0.0
        }
    });

    let dx = 1.0 / 100.0;

    let dt = 0.001;

    let result = solve_burgers_1d(&mut mesh, dt, 50);

    // Shock should have moved to the right and stayed sharp-ish but slightly smoothed by LF
    assert!(result[60] > 0.0);
}

#[test]

fn test_swe_1d_dam_break() {

    let n = 100;

    let mut h = vec![1.0; n];

    for i in 50 .. n {

        h[i] = 0.5;
    } // Dam break setup
    let hu = vec![0.0; n];

    let dx = 1.0 / n as f64;

    let dt = 0.001;

    let result = solve_shallow_water_1d(
        h, hu, dx, dt, 50, 9.81,
    );

    assert!(result[45].h < 1.0); // Rarefaction wave
    assert!(result[55].h > 0.5); // Bore (shock) wave
}

proptest! {
    #[test]
    fn prop_advection_conserved(
        v in 0.1..2.0f64,
        steps in 1usize..20usize
    ) {
        let n = 50;
        let mut mesh = Mesh::new(n, 1.0, |_| 1.0);
        let dx = 1.0 / n as f64;
        let dt = 0.1 * dx / v;
        let initial_sum: f64 = mesh.cells.iter().map(|c| c.value).sum();
        let result = solve_advection_1d(&mut mesh, v, dt, steps, || (1.0, 1.0));
        let final_sum: f64 = result.iter().sum();
        // Since we have constant inflow/outflow of 1.0, sum should stay roughly constant
        prop_assert!((final_sum - initial_sum).abs() < 1e-10);
    }
}
