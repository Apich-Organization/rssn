//! Rigid MD rotor preflight and composition regressions.

use rssn::sim::md_engine::ForceField;
use rssn::sim::md_engine::MdSystem;
use rssn::sim::md_engine::SimBox;
use rssn::sim::md_rigid::RigidMdSimulation;
use rssn::sim::symplectic::RigidBody;
use rssn::sim::symplectic::composition_weights;
use rssn::sim::symplectic::quat_to_matrix;

fn free_rotor(dt: f64) -> RigidMdSimulation {
    let positions = vec![
        [3.1, 5.2, 7.3],
        [3.8, 5.1, 7.25],
        [3.0, 5.8, 7.45],
        [3.32, 5.58, 7.86],
    ];
    let mut system = MdSystem::new(positions, SimBox::periodic([24.0; 3]));
    system.mass = vec![1.0, 1.25, 0.75, 0.5];
    system.vel = vec![
        [0.3, -0.2, 0.1],
        [0.1, 0.5, -0.2],
        [0.2, 0.3, 0.56],
        [-0.2, 0.1, 0.46],
    ];
    RigidMdSimulation::new(
        system,
        ForceField::new(Vec::new(), 2.0),
        0.3,
        dt,
        &[vec![0, 1, 2, 3]],
    )
    .unwrap()
}

#[derive(Debug, PartialEq)]
struct Snapshot {
    time: u64,
    orientation: [f64; 4],
    angular_momentum: [f64; 3],
    center: [f64; 3],
    center_velocity: [f64; 3],
    wrapped: Vec<[f64; 3]>,
    unwrapped: Vec<[f64; 3]>,
    velocity: Vec<[f64; 3]>,
    forces: Vec<[f64; 3]>,
}

fn snapshot(simulation: &RigidMdSimulation) -> Snapshot {
    let molecule = &simulation.mols[0];
    Snapshot {
        time: simulation.sys.time.to_bits(),
        orientation: molecule.rb.q,
        angular_momentum: molecule.rb.l_body,
        center: molecule.com,
        center_velocity: molecule.vcom,
        wrapped: simulation.sys.pos.clone(),
        unwrapped: simulation.sys.unwrapped.clone(),
        velocity: simulation.sys.vel.clone(),
        forces: simulation.forces.f.clone(),
    }
}

#[test]
fn invalid_composition_does_not_change_the_body_or_sites() {
    let mut simulation = free_rotor(0.125);
    simulation.order = 3;
    let before = snapshot(&simulation);
    for exact in [false, true] {
        simulation.exact_rotor = exact;
        assert!(simulation.step().unwrap_err().contains("order"));
        assert_eq!(snapshot(&simulation), before);
        assert!(simulation.run(3).unwrap_err().contains("order"));
        assert_eq!(snapshot(&simulation), before);
    }
    simulation.run(0).unwrap();
    assert_eq!(snapshot(&simulation), before);
}

#[test]
fn invalid_time_and_substep_overflow_are_detected_before_motion() {
    let mut simulation = free_rotor(0.125);
    simulation.dt = f64::NAN;
    let before = snapshot(&simulation);
    assert!(simulation.step().unwrap_err().contains("must be finite"));
    assert_eq!(snapshot(&simulation), before);

    simulation.order = 4;
    simulation.dt = f64::MAX;
    assert!(
        simulation
            .run(2)
            .unwrap_err()
            .contains("rotational substep")
    );
    assert_eq!(snapshot(&simulation), before);

    simulation.order = 2;
    simulation.sys.time = f64::MAX;
    let before_overflow = snapshot(&simulation);
    assert!(simulation.step().unwrap_err().contains("step time"));
    assert_eq!(snapshot(&simulation), before_overflow);
}

#[test]
fn all_rotor_modes_follow_their_reference_composition_on_signed_steps() {
    for (order, top, exact) in [
        (2, false, false),
        (6, true, false),
        (8, false, true),
        (10, true, false),
    ] {
        let mut simulation = free_rotor(0.0125);
        let initial = &simulation.mols[0];
        let mut rotor = RigidBody {
            q: initial.rb.q,
            l_body: initial.rb.l_body,
            inertia: initial.rb.inertia,
        };
        let mut center = initial.com;
        let velocity = initial.vcom;
        let reference_weights = composition_weights(order).unwrap();
        let mut time = 0.0_f64;
        simulation.order = order;
        simulation.top_split = top;
        simulation.exact_rotor = exact;
        for (dt, steps) in [(0.0125, 7), (-0.0175, 5)] {
            simulation.dt = dt;
            simulation.run(steps).unwrap();
            for _ in 0..steps {
                for (coordinate, velocity) in center.iter_mut().zip(velocity) {
                    *coordinate += dt * velocity;
                }
                for &coefficient in &reference_weights {
                    let interval = coefficient * dt;
                    if exact {
                        rotor.free_flow_exact(interval);
                    } else if top {
                        rotor.free_step_top(interval);
                    } else {
                        rotor.free_step(interval);
                    }
                }
                time += dt;
            }
            let actual = &simulation.mols[0];
            assert_eq!(simulation.sys.time.to_bits(), time.to_bits());
            assert_eq!(actual.com, center);
            assert_eq!(actual.vcom, velocity);
            assert_eq!(actual.rb.q, rotor.q);
            assert_eq!(actual.rb.l_body, rotor.l_body);
            let matrix = quat_to_matrix(&rotor.q);
            for (&site, body) in actual.sites.iter().zip(&actual.body) {
                let expected: [f64; 3] = std::array::from_fn(|axis| {
                    center[axis]
                        + matrix[axis]
                            .iter()
                            .zip(body)
                            .map(|(a, b)| a * b)
                            .sum::<f64>()
                });
                for (unwrapped, expected) in simulation.sys.unwrapped[site].iter().zip(expected) {
                    assert!((unwrapped - expected).abs() < 2e-14, "site {site}");
                }
            }
        }
    }
}
