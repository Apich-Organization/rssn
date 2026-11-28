//! Pair-field checks for the unsoftened chain coordinates.

use rssn::sim::chain::ChainSystem;

fn check_relative(
    actual: f64,
    expected: f64,
) {
    assert!(actual.is_finite(), "non-finite output for {expected:e}");
    if expected == 0.0 {
        assert_eq!(actual.abs().to_bits(), 0);
    } else {
        assert!(
            ((actual - expected) / expected).abs() <= 2e-14,
            "{actual:e} vs {expected:e}"
        );
    }
}

fn check_pair(
    g: f64,
    masses: [f64; 2],
    length: f64,
    direction: [f64; 3],
) {
    let positions = [[0.0; 3], direction.map(|v| v * length)];
    let chain = ChainSystem::new(g, &masses, &positions, &[[0.0; 3]; 2]).unwrap();

    // Evaluate in units of the supplied length, independent of the
    // implementation's raw distance and inverse-cube calculations.
    let radius: f64 = direction.iter().map(|v| v * v).sum::<f64>().sqrt();
    let energy = -(g / length) / radius * masses[0] * masses[1];
    check_relative(chain.potential(), energy);
    check_relative(chain.b, -energy);
    let scaled_g = ((g / length) / length) / (radius * radius * radius);
    let a = chain.accelerations();
    for (axis, &direction) in direction.iter().enumerate() {
        let field = scaled_g * direction;
        check_relative(a[0][axis], field * masses[1]);
        check_relative(a[1][axis], -field * masses[0]);
    }
}

#[test]
fn pair_accelerations_and_binding_survive_extreme_separations_and_mass_products() {
    check_pair(1e300, [2.0, 3.0], 1e200, [1.0, 2.0, -1.0]);
    check_pair(1e-300, [2.0, 3.0], 1e-200, [1.0, -2.0, 1.0]);
    check_pair(1e300, [2e10, 3e10], 1e160, [1.0, 0.0, 0.0]);
    // Euclidean length exceeds f64::MAX, even though the components and
    // both physical outputs can be represented.
    check_pair(1e308, [1.0, 1.0], 1e308, [1.4, 1.5, 0.0]);
}

#[test]
fn a_small_force_component_can_be_recovered_after_numerator_underflow() {
    check_pair(1e-60, [1e-60, 2e-60], 1e-100, [1.0, 1e-150, 0.0]);
}

#[test]
fn reported_original_order_preserves_pair_directions_and_mass_weighted_balance() {
    let g = 1.75;
    let masses = [1.0, 0.75, 2.5, 1.4];
    let positions = [
        [2.0, 1.0, 8.0],
        [0.0, 0.0, 0.0],
        [-4.0, -1.0, 2.0],
        [0.25, 0.1, -0.1],
    ];
    let c = ChainSystem::new(g, &masses, &positions, &[[0.0; 3]; 4]).unwrap();
    assert_ne!(c.order, vec![0, 1, 2, 3]);
    let a = c.accelerations();
    let mut balance = [0.0_f64; 3];
    for i in 0..positions.len() {
        let mut reference = [0.0; 3];
        for j in 0..positions.len() {
            if i == j {
                continue;
            }
            let r = std::array::from_fn::<_, 3, _>(|d| positions[j][d] - positions[i][d]);
            let radius = r.iter().map(|v| v * v).sum::<f64>().sqrt();
            let factor = g * masses[j] / radius.powi(3);
            for d in 0..3 {
                reference[d] += factor * r[d];
            }
        }
        for d in 0..3 {
            check_relative(a[i][d], reference[d]);
            balance[d] += masses[i] * a[i][d];
        }
    }
    assert!(balance.iter().all(|v| v.abs() < 1e-13));
}

#[test]
fn closest_pair_is_selected_for_large_finite_distances() {
    let positions = [[0.0; 3], [3e160, 0.0, 0.0], [3.1e160, 0.0, 0.0]];
    let c = ChainSystem::new(1e300, &[1.0; 3], &positions, &[[0.0; 3]; 3]).unwrap();
    let slot1 = c.order.iter().position(|&i| i == 1).unwrap();
    let slot2 = c.order.iter().position(|&i| i == 2).unwrap();
    assert_eq!(slot1.abs_diff(slot2), 1);
    assert!(c.potential().is_finite() && c.potential() < 0.0);
}

#[test]
fn coincident_initial_states_error_and_separated_zero_gravity_has_zero_force() {
    let pos = [[0.0; 3], [1.0, 2.0, 3.0]];
    for g in [0.0, 1.0] {
        let error = ChainSystem::new(g, &[2.0, 3.0], &[[1.0; 3]; 2], &[[0.0; 3]; 2]).unwrap_err();
        assert!(error.contains("coincident"), "{error}");
    }
    let c = ChainSystem::new(0.0, &[2.0, 3.0], &pos, &[[0.0; 3]; 2]).unwrap();
    assert_eq!(c.potential().to_bits(), (-0.0_f64).to_bits());
    assert_eq!(c.accelerations(), vec![[0.0; 3]; 2]);

    let mut singular = ChainSystem::new(1.0, &[2.0, 3.0], &pos, &[[0.0; 3]; 2]).unwrap();
    singular.xc[0] = [0.0; 3];
    assert_eq!(singular.potential().to_bits(), f64::NEG_INFINITY.to_bits());
    assert!(
        singular
            .accelerations()
            .iter()
            .flatten()
            .all(|v| v.is_nan())
    );
}

#[test]
fn nonfinite_cartesian_data_and_unrepresentable_links_are_errors() {
    let velocities = [[0.0; 3]; 2];
    let pos = [[0.0; 3], [1.0, 2.0, 3.0]];
    assert!(ChainSystem::new(1.0, &[f64::NAN, 1.0], &pos, &velocities).is_err());
    assert!(ChainSystem::new(f64::INFINITY, &[1.0; 2], &pos, &velocities).is_err());
    assert!(
        ChainSystem::new(
            1.0,
            &[1.0; 2],
            &[[-f64::MAX, 0.0, 0.0], [f64::MAX, 0.0, 0.0]],
            &velocities,
        )
        .is_err()
    );
}
