//! Unit and property-based tests for the physics sim geodesic relativity module.

use proptest::prelude::*;
use rssn::physics::physics_sim::geodesic_relativity::*;

#[test]

fn test_geodesic_simulation_smoke() {

    let params = GeodesicParameters {
        black_hole_mass: 1.0,
        initial_state: [
            10.0, 0.0, 0.0, 0.035,
        ],
        proper_time_end: 100.0,
        initial_dt: 0.1,
    };

    let path = run_geodesic_simulation(
        &params,
    );

    assert!(path.len() > 0);

    // Initial position should be (10, 0) in Cartesian
    let (x0, y0) = path[0];

    assert!((x0 - 10.0).abs() < 1e-10);

    assert!(y0.abs() < 1e-10);
}

#[test]

fn test_effective_potential_scaling() {

    let params = GeodesicParameters {
        black_hole_mass: 1.0,
        initial_state: [
            10.0, 0.0, 0.0, 0.035,
        ],
        proper_time_end: 10.0,
        initial_dt: 0.1,
    };

    // Test that effective potential is correct at some point
    // V_eff(10, 0.35) = -1/10 + (10^2 * 0.35^2)/(2*10^2) - 1*(100*0.1225)/1000
    // L = r^2 * phi_dot = 100 * 0.035 = 3.5
    // V_eff = -1/10 + 3.5^2 / 200 - 3.5^2 / 1000
    // V_eff = -0.1 + 12.25/200 - 12.25/1000 = -0.1 + 0.06125 - 0.01225 = -0.051
    let v = params
        .effective_potential(10.0, 3.5);

    assert!(
        (v - (-0.051)).abs() < 1e-10
    );
}

#[test]

fn test_black_hole_orbits_scenario_run()
{

    let res = simulate_black_hole_orbits_scenario();

    assert!(res.is_ok());
}

proptest! {
    #[test]
    fn prop_geodesic_stability(mass in 0.5..2.0f64, r0 in 6.0..20.0f64) {
        let params = GeodesicParameters {
            black_hole_mass: mass,
            initial_state: [r0, 0.0, 0.0, 0.01],
            proper_time_end: 50.0,
            initial_dt: 0.1,
        };

        let path = run_geodesic_simulation(&params);
        for (x, y) in path {
            prop_assert!(x.is_finite());
            prop_assert!(y.is_finite());
        }
    }
}
