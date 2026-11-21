//! Regression tests for Hermite time-step selection and scheduling.

use rssn::sim::hermite::HermiteConfig;
use rssn::sim::hermite::HermiteOrder;
use rssn::sim::hermite::HermiteSystem;
use rssn::sim::hermite::StepPolicy;

fn system(
    g: f64,
    order: HermiteOrder,
    policy: StepPolicy,
    t0: f64,
) -> Result<HermiteSystem, String> {
    HermiteSystem::new(
        vec![1.0, 0.0],
        vec![[0.0; 3], [1.0, 0.0, 0.0]],
        vec![[0.0; 3], [0.0, 1.0, 0.0]],
        &HermiteConfig {
            g,
            eps: 0.0,
            order,
            policy,
            t0,
        },
    )
}

fn block(
    eta: f64,
    dt_min: f64,
    dt_max: f64,
) -> StepPolicy {
    StepPolicy::Block { eta, dt_max, dt_min }
}

fn binary_steps() -> impl Iterator<Item = f64> {
    let min_step = f64::MIN_POSITIVE * (1.0 / 4_503_599_627_370_496.0);
    std::iter::successors(Some(min_step), |&step| {
        let next = 2.0 * step;
        next.is_finite().then_some(next)
    })
}

#[test]
fn block_upper_bound_uses_the_exact_binary_lattice() {
    let order = HermiteOrder::Fourth;
    let powers: Vec<_> = binary_steps().collect();
    for (index, &power) in powers.iter().enumerate() {
        let upper = if index == 0 {
            power
        } else {
            power.next_down()
        };
        let expected = if index == 0 {
            power
        } else {
            powers[index - 1]
        };
        let at_edge = system(0.0, order, block(0.1, powers[0], upper), 0.0).unwrap();
        assert_eq!(at_edge.dt, [expected; 2], "upper bound {upper:e}");

        let at_power = system(0.0, order, block(0.1, powers[0], power), 0.0).unwrap();
        assert_eq!(at_power.dt, [power; 2], "upper bound {power:e}");
    }
}

#[test]
fn doubling_requires_a_permitted_power_at_the_block_time() {
    let max = 0.125_f64.next_down();
    let mut simulation = system(1.0, HermiteOrder::Fourth, block(1e8, 1e-6, max), 0.0625).unwrap();
    assert_eq!(simulation.dt, [0.0625; 2]);
    simulation.step().unwrap();
    assert_eq!(simulation.time, 0.125);
    assert_eq!(simulation.dt, [0.0625; 2]);

    // With the next power now allowed, an unaligned endpoint still prevents doubling.
    simulation.policy = block(1e8, 1e-6, 0.125);
    simulation.step().unwrap();
    assert_eq!(simulation.time, 0.1875);
    assert_eq!(simulation.dt, [0.0625; 2]);
    simulation.step().unwrap();
    assert_eq!(simulation.time, 0.25);
    assert_eq!(simulation.dt, [0.125; 2]);
}

#[test]
fn first_step_and_following_steps_respect_nonbinary_lower_bound() {
    let requested_min = 1e-6;
    let requested_max = 0.125_f64.next_down();
    let permitted: Vec<_> = binary_steps()
        .filter(|&h| (requested_min..=requested_max).contains(&h))
        .collect();
    let order = HermiteOrder::Sixth;
    let mut simulation =
        system(1.0, order, block(1e-12, requested_min, requested_max), 0.0).unwrap();
    assert_eq!(simulation.dt, [*permitted.last().unwrap(), permitted[0]]);

    let end = 0.0001_f64;
    let n = simulation.mass.len() as u64;
    let attempted = (end / requested_min).ceil() as u64 + 1;
    let mut blocks = 0;
    while simulation.next_time() <= end {
        simulation.step().unwrap();
        assert!(simulation.dt.iter().all(|h| permitted.contains(h)));
        blocks += 1;
        assert!(blocks <= attempted);
    }
    assert!(blocks > 1);
    // At most N active targets per nominal attempted block. The native charge allows
    // N^2 work per target (including the optional 1PN correction).
    assert!(simulation.updates * n * n <= attempted * n * n * n);
}

#[test]
fn no_block_is_shortened_to_reach_an_off_grid_horizon() {
    let mut simulation = system(0.0, HermiteOrder::Fourth, block(0.1, 0.08, 0.14), -0.25).unwrap();
    assert_eq!(simulation.dt, [0.125; 2]);
    let requested_end = 0.03;
    simulation.advance_to(requested_end).unwrap();
    assert_eq!(simulation.time, 0.0);
    assert_eq!(simulation.next_time(), 0.125);
    assert_eq!(simulation.updates, 4);
    let (positions, velocities) = simulation.state_at(requested_end).unwrap();
    let analytic = requested_end + 0.25;
    assert!((positions[1][1] - analytic).abs() <= f64::EPSILON);
    assert_eq!(velocities, [[0.0; 3], [0.0, 1.0, 0.0]]);
}

#[test]
fn infeasible_bounds_and_nonfinite_configurations_are_rejected() {
    let order = HermiteOrder::Fourth;
    for policy in [
        block(0.1, 0.063, 0.124),
        block(0.1, f64::MAX, f64::MAX),
        block(0.1, f64::from_bits(3), f64::from_bits(3)),
    ] {
        assert!(
            system(0.0, order, policy, 0.0)
                .unwrap_err()
                .contains("power of two")
        );
    }
    for policy in [
        block(f64::NAN, 0.01, 1.0),
        block(f64::INFINITY, 0.01, 1.0),
        block(0.1, 0.0, 1.0),
        block(0.1, 0.01, f64::INFINITY),
        block(0.1, f64::NAN, 1.0),
        StepPolicy::Shared(f64::NAN),
        StepPolicy::Shared(f64::INFINITY),
    ] {
        assert!(system(0.0, order, policy, 0.0).is_err());
    }
    for initial_time in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(system(0.0, order, StepPolicy::Shared(1.0), initial_time).is_err());
    }
    let mut simulation = system(0.0, order, StepPolicy::Shared(1.0), 0.0).unwrap();
    for end_time in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(simulation.advance_to(end_time).is_err());
        assert_eq!(simulation.time, 0.0);
        assert_eq!(simulation.updates, 0);
    }
}

#[test]
fn an_unrepresentable_next_block_leaves_the_system_unchanged() {
    let stalled_time = 9_007_199_254_740_992.0;
    for (initial, step) in [(stalled_time, 1.0), (f64::MAX, f64::MAX)] {
        let mut simulation =
            system(0.0, HermiteOrder::Fourth, StepPolicy::Shared(step), initial).unwrap();
        let original = simulation.clone();
        assert!(simulation.step().unwrap_err().contains("later finite time"));
        if initial == stalled_time {
            assert!(simulation.advance_to(initial).is_err());
        }
        assert_eq!(simulation.time, original.time);
        assert_eq!(simulation.t_last, original.t_last);
        assert_eq!(simulation.dt, original.dt);
        assert_eq!(simulation.pos, original.pos);
        assert_eq!(simulation.vel, original.vel);
        assert_eq!(simulation.acc, original.acc);
        assert_eq!(simulation.jerk, original.jerk);
        assert_eq!(simulation.snap, original.snap);
        assert_eq!(simulation.crackle, original.crackle);
        assert_eq!(simulation.interactions, 0);
        assert_eq!(simulation.updates, 0);
    }
}
