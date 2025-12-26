//! Unit and property-based tests for the physics sim linear elasticity module.

use proptest::prelude::*;
use rssn::physics::physics_sim::linear_elasticity::*;

#[test]

fn test_stiffness_matrix_symmetry() {

    let e = 1e7;

    let nu = 0.3;

    let p1 = (0.0, 0.0);

    let p2 = (1.0, 0.0);

    let p3 = (1.0, 1.0);

    let p4 = (0.0, 1.0);

    let k = element_stiffness_matrix(
        p1, p2, p3, p4, e, nu,
    );

    for i in 0..8 {

        for j in 0..8 {

            assert!((k[[i, j]] - k[[j, i]]).abs() < 1e-10);
        }
    }
}

#[test]

fn test_simulate_cantilever_beam_scenario() {

    println!("Running 2D Cantilever Beam simulation...");

    simulate_cantilever_beam_scenario().unwrap();
}

#[test]

fn test_run_elasticity_simulation_basic() {

    let nodes = vec![
        (0.0, 0.0),
        (1.0, 0.0),
        (1.0, 1.0),
        (0.0, 1.0),
    ];

    let elements = vec![[0, 1, 2, 3]];

    let fixed_nodes = vec![0, 3];

    let loads = vec![
        (1, 1000.0, 0.0),
        (2, 1000.0, 0.0),
    ];

    let params = ElasticityParameters {
        nodes,
        elements,
        youngs_modulus: 1e7,
        poissons_ratio: 0.3,
        fixed_nodes,
        loads,
    };

    let res = run_elasticity_simulation(&params).unwrap();

    assert_eq!(res.len(), 8);

    // Fixed nodes should have zero displacement
    assert_eq!(res[0], 0.0);

    assert_eq!(res[1], 0.0);

    assert_eq!(res[6], 0.0);

    assert_eq!(res[7], 0.0);

    // Pulling to the right, x-displacements of nodes 1 and 2 should be positive
    assert!(res[2] > 0.0);

    assert!(res[4] > 0.0);
}

proptest! {
    #[test]
    fn prop_stiffness_scaling(scale in 0.1..10.0f64) {
        let e1 = 1e7;
        let e2 = e1 * scale;
        let nu = 0.3;
        let p1 = (0.0, 0.0);
        let p2 = (1.0, 0.0);
        let p3 = (1.0, 1.0);
        let p4 = (0.0, 1.0);
        let k1 = element_stiffness_matrix(p1, p2, p3, p4, e1, nu);
        let k2 = element_stiffness_matrix(p1, p2, p3, p4, e2, nu);

        for i in 0..8 {
            for j in 0..8 {
                prop_assert!((k1[[i, j]] * scale - k2[[i, j]]).abs() < 1e-5);
            }
        }
    }
}
