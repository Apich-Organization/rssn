//! Scale-aware local and global Euler–Bernoulli beam stiffness.

use std::f64::consts::FRAC_PI_2;

use rssn::sim::physics_fea::BeamElement2D;
use rssn::sim::physics_fea::solve_static_structural;

fn close(
    got: f64,
    expected: f64,
) {
    assert!(
        (got - expected).abs() <= expected.abs() * 3e-14 + f64::MIN_POSITIVE,
        "{got:e} vs {expected:e}"
    );
}

#[test]
fn axial_and_bending_stiffness_survive_both_ends_of_the_scale_range() {
    for (length, young, area, inertia, axial) in [
        (1e110, 1e200, 1e-80, 1e130, 1e10),
        (1e-110, 1e-200, 1e80, 1e-130, 1e-10),
    ] {
        let beam = BeamElement2D::new(length, young, area, inertia, 0.0);
        let local = beam.local_stiffness_matrix();
        // EI/L³ equals one at both scales, even though the separate terms do not fit.
        assert!(local.data().iter().all(|value| value.is_finite()));
        close(*local.get(0, 0), axial);
        close(*local.get(1, 1), 12.0);
        close(*local.get(1, 2), 6.0 * length);
        close(*local.get(2, 2), 4.0 * length * length);
        close(*local.get(2, 5), 2.0 * length * length);
        let global = beam.global_stiffness_entries();
        for (a, b) in global.iter().flatten().zip(local.data()) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
    }

    for scale in [1e-200, 1e200] {
        // EA falls outside range on its own, and I = 0 suppresses bending.
        let axial = BeamElement2D::new(scale, scale, scale, 0.0, 0.0).local_stiffness_matrix();
        close(*axial.get(0, 0), scale);
        assert_eq!(axial.get(1, 1).to_bits(), 0);
        assert!(axial.data().iter().all(|value| value.is_finite()));
    }
}

#[test]
fn global_stiffness_has_the_analytic_energy_and_dof_order() {
    let (length, young, area, inertia) = (2.0, 12.0, 0.75, 0.25);
    let displacement = [0.27, -0.63, 0.17, 0.91, 0.55, -0.42];
    for angle in [-FRAC_PI_2, -0.54, 0.0, 0.63, FRAC_PI_2] {
        let beam = BeamElement2D::new(length, young, area, inertia, angle);
        let a = beam.global_stiffness_entries();
        let matrix = beam.global_stiffness_matrix();
        for (i, row) in a.iter().enumerate() {
            for (j, value) in row.iter().enumerate() {
                assert_eq!(value.to_bits(), matrix.get(i, j).to_bits());
                assert!((value - a[j][i]).abs() < 1e-13);
            }
        }

        let (s, c) = angle.sin_cos();
        let dx = displacement[3] - displacement[0];
        let dy = displacement[4] - displacement[1];
        let parallel = c * dx + s * dy;
        let transverse = -s * dx + c * dy;
        let bending = transverse - (length / 2.0) * (displacement[2] + displacement[5]);
        let twist = displacement[5] - displacement[2];
        let exact = (young * area / length) * parallel * parallel
            + (12.0 * young * inertia / length.powi(3)) * bending * bending
            + (young * inertia / length) * twist * twist;
        let energy = a
            .iter()
            .zip(displacement)
            .map(|(row, u)| {
                u * row
                    .iter()
                    .zip(displacement)
                    .map(|(k, v)| k * v)
                    .sum::<f64>()
            })
            .sum::<f64>();
        assert!(energy >= 0.0);
        close(energy, exact);
        if angle.abs() > 1.5 {
            // At ±90°, the global y DOF is axial and the x–rotation term flips sign.
            close(a[1][1], young * area / length);
            close(
                a[0][2],
                -angle.signum() * 6.0 * young * inertia / length.powi(2),
            );
        }
        let rigid_rotation = [0.0, 0.0, 1.0, -s * length, c * length, 1.0];
        assert!(a.iter().all(|row| {
            row.iter()
                .zip(rigid_rotation)
                .map(|(k, u)| k * u)
                .sum::<f64>()
                .abs()
                < 2e-14
        }));
    }
}

#[test]
fn tilted_cantilever_tip_displacement_matches_separate_axial_and_bending_solutions() {
    let (length, young, area, inertia) = (2.5_f64, 3000.0, 0.2, 0.013);
    let angle = 0.63_f64;
    let (s, c) = angle.sin_cos();
    let (axial_force, transverse_force) = (3.1, 5.7);
    let beam = BeamElement2D::new(length, young, area, inertia, angle);
    let force = vec![
        0.0,
        0.0,
        0.0,
        axial_force * c - transverse_force * s,
        axial_force * s + transverse_force * c,
        0.0,
    ];
    let u = solve_static_structural(
        beam.global_stiffness_matrix(),
        force,
        &[(0, 0.0), (1, 0.0), (2, 0.0)],
    )
    .unwrap();
    close(c * u[3] + s * u[4], axial_force * length / (young * area));
    close(
        -s * u[3] + c * u[4],
        transverse_force * length.powi(3) / (3.0 * young * inertia),
    );
    close(
        u[5],
        transverse_force * length.powi(2) / (2.0 * young * inertia),
    );
}

#[test]
fn signed_zero_and_negative_length_keep_the_legacy_ieee_convention() {
    let beam = BeamElement2D::new(-2.0, 4.0, 3.0, 0.25, -0.0);
    let negative = beam.local_stiffness_matrix();
    for ((row, col), expected) in [
        ((0, 0), -6.0),
        ((1, 1), -1.5),
        ((1, 2), 1.5),
        ((2, 2), -2.0),
    ] {
        close(*negative.get(row, col), expected);
    }
    for (a, b) in beam
        .global_stiffness_entries()
        .iter()
        .flatten()
        .zip(negative.data())
    {
        assert_eq!(a.to_bits(), b.to_bits());
    }
    for length in [0.0_f64, -0.0] {
        let beam = BeamElement2D::new(length, 4.0, 3.0, 0.25, 0.0);
        let local = beam.local_stiffness_matrix();
        assert_eq!(
            local.get(0, 0).to_bits(),
            f64::INFINITY.copysign(length).to_bits()
        );
        assert_eq!(
            local.get(1, 1).to_bits(),
            f64::INFINITY.copysign(length).to_bits()
        );
        assert_eq!(local.get(1, 2).to_bits(), f64::INFINITY.to_bits());
        assert!(
            beam.global_stiffness_entries()
                .iter()
                .flatten()
                .all(|value| value.is_nan())
        );
    }
}
