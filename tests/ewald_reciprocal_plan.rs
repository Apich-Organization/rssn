//! Periodicity, scaling and ordered reduction of the direct Ewald plan.

use std::f64::consts::PI;

use rssn::sim::ewald::EwaldReciprocalPlan;
use rssn::sim::ewald::Recip;
use rssn::sim::ewald::diag_cell;
use rssn::sim::ewald::ewald_reciprocal_cell;

const CELL: [[f64; 3]; 3] = [[4.0, 0.0, 0.0], [1.0, 3.0, 0.0], [0.5, 0.75, 5.0]];
const POSITION: [[f64; 3]; 4] = [
    [0.125, 0.375, 0.75],
    [1.25, 0.625, 1.125],
    [2.125, 1.0, 2.75],
    [3.25, 0.5, 3.75],
];
const CHARGE: [f64; 4] = [0.75, -0.5, 0.375, -1.0];

fn close(
    reference: &Recip,
    actual: &Recip,
    energy_factor: f64,
    force_factor: f64,
    tolerance: f64,
) {
    let check = |x: f64, y: f64| {
        assert!((x - y).abs() <= tolerance * (1.0 + x.abs()), "{x} vs {y}");
    };
    check(reference.energy, actual.energy * energy_factor);
    for (a, b) in reference
        .virial
        .iter()
        .flatten()
        .zip(actual.virial.iter().flatten())
    {
        check(*a, *b * energy_factor);
    }
    for (a, b) in reference
        .forces
        .iter()
        .flatten()
        .zip(actual.forces.iter().flatten())
    {
        check(*a, *b * force_factor);
    }
}

#[test]
fn plans_and_legacy_entry_points_match_under_large_lattice_shifts() {
    for h in [CELL, diag_cell([4.0, 3.0, 5.0])] {
        let plan = EwaldReciprocalPlan::new(&h, 1.1, 8).unwrap();
        let reference = plan.evaluate(&POSITION, &CHARGE, 0.75).unwrap();
        let shifted = std::array::from_fn::<_, 4, _>(|i| {
            let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
            let turns = [
                sign * (1_u64 << 34) as f64,
                -sign * (1_u64 << 35) as f64,
                31.0,
            ];
            std::array::from_fn(|a| {
                POSITION[i][a] + turns[0] * h[0][a] + turns[1] * h[1][a] + turns[2] * h[2][a]
            })
        });
        let result = plan.evaluate(&shifted, &CHARGE, 0.75).unwrap();
        close(&reference, &result, 1.0, 1.0, 1.5e-13);
        let wrapped = ewald_reciprocal_cell(&shifted, &CHARGE, &h, 1.1, 8, 0.75).unwrap();
        close(&reference, &wrapped, 1.0, 1.0, 1.5e-13);
        for a in 0..3 {
            let net: f64 = result.forces.iter().map(|force| force[a]).sum();
            assert!(net.abs() < 2e-14);
            for b in 0..3 {
                assert_eq!(result.virial[a][b].to_bits(), result.virial[b][a].to_bits());
            }
        }
    }
}

#[test]
fn normalized_geometry_and_charges_preserve_physical_homogeneity() {
    let reference = EwaldReciprocalPlan::new(&CELL, 1.0, 7)
        .unwrap()
        .evaluate(&POSITION, &CHARGE, 0.8)
        .unwrap();
    for power in [-400, 400] {
        let length = 2.0_f64.powi(power);
        let h = CELL.map(|row| row.map(|value| value * length));
        let pos = POSITION.map(|row| row.map(|value| value * length));
        let result = EwaldReciprocalPlan::new(&h, 1.0 / length, 7)
            .unwrap()
            .evaluate(&pos, &CHARGE, 0.8)
            .unwrap();
        close(&reference, &result, length, length * length, 2e-14);
    }
    let charge_scale = 2.0_f64.powi(550);
    let kc_scale = 2.0_f64.powi(-900);
    let result = EwaldReciprocalPlan::new(&CELL, 1.0, 7)
        .unwrap()
        .evaluate(&POSITION, &CHARGE.map(|q| q * charge_scale), 0.8 * kc_scale)
        .unwrap();
    close(
        &reference,
        &result,
        2.0_f64.powi(-200),
        2.0_f64.powi(-200),
        2e-14,
    );
}

#[test]
fn plan_tracks_reciprocal_energy_force_and_strain() {
    let plan = EwaldReciprocalPlan::new(&CELL, 1.25, 8).unwrap();
    let q = [1.0, -1.0, 0.5, -0.5];
    let original = plan.evaluate(&POSITION, &q, 1.0).unwrap();
    let eps = 2.0e-6;
    for (particle, axis) in [(0, 1), (1, 0), (3, 2)] {
        let mut pos = POSITION;
        pos[particle][axis] += eps;
        let plus = plan.evaluate(&pos, &q, 1.0).unwrap().energy;
        pos[particle][axis] -= 2.0 * eps;
        let minus = plan.evaluate(&pos, &q, 1.0).unwrap().energy;
        let force = -(plus - minus) / (2.0 * eps);
        assert!((original.forces[particle][axis] - force).abs() < 5e-9);
    }
    // An x stretch of both cell and particle coordinates leaves the self term
    // invariant. Neutral charges also remove changes in the background energy.
    let strain = |factor| {
        let stretch = |[x, y, z]: [f64; 3]| [factor * x, y, z];
        EwaldReciprocalPlan::new(&CELL.map(stretch), 1.25, 8)
            .unwrap()
            .evaluate(&POSITION.map(stretch), &q, 1.0)
            .unwrap()
            .energy
    };
    let derivative = -(strain(1.0 + eps) - strain(1.0 - eps)) / (2.0 * eps);
    assert!((original.virial[0][0] - derivative).abs() < 5e-9);
}

#[test]
fn bounded_block_reduction_keeps_bits_for_distinct_worker_pools() {
    let pos: Vec<[f64; 3]> = (0..96)
        .map(|i| {
            [
                (i * 19 % 103) as f64 / 32.0,
                (i * 17 % 97) as f64 / 32.0,
                (i * 23 % 137) as f64 / 32.0,
            ]
        })
        .collect();
    let q: Vec<f64> = (0..pos.len())
        .map(|i| (i as f64 % 5.0) / 4.0 - 0.5)
        .collect();
    let plan = EwaldReciprocalPlan::new(&CELL, 1.2, 10).unwrap();
    let run = |n| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build()
            .unwrap()
            .install(|| plan.evaluate(&pos, &q, 1.0).unwrap())
    };
    let a = run(1);
    for n in [2, 4] {
        let b = run(n);
        assert_eq!(a.energy.to_bits(), b.energy.to_bits());
        for (x, y) in a.forces.iter().flatten().zip(b.forces.iter().flatten()) {
            assert_eq!(x.to_bits(), y.to_bits());
        }
        for (x, y) in a.virial.iter().flatten().zip(b.virial.iter().flatten()) {
            assert_eq!(x.to_bits(), y.to_bits());
        }
    }
}

#[test]
fn zero_modes_and_checked_errors_preserve_self_and_background_semantics() {
    let len = [2.0, 3.0, 4.0];
    let h = diag_cell(len);
    let plan = EwaldReciprocalPlan::new(&h, 1.5, 0).unwrap();
    let result = plan.evaluate(&[[0.0; 3]], &[2.0], -0.5).unwrap();
    let expected = 3.0 / PI.sqrt() + PI / (24.0 * 1.5 * 1.5);
    assert!((result.energy - expected).abs() < 2e-15);
    assert_eq!(result.forces, vec![[0.0; 3]]);
    assert_eq!(result.virial, [[0.0; 3]; 3]);
    assert!(plan.evaluate(&POSITION, &CHARGE[..3], 1.0).is_err());
    assert!(plan.evaluate(&[[f64::NAN, 0.0, 0.0]], &[1.0], 1.0).is_err());
    assert!(plan.evaluate(&[[0.0; 3]], &[f64::INFINITY], 1.0).is_err());
    assert!(EwaldReciprocalPlan::new(&h, 0.0, 2).is_err());
    assert!(EwaldReciprocalPlan::new(&h, 1.0, i32::MAX).is_err());
    assert!(EwaldReciprocalPlan::new(&h, 1.0, -1).is_err());
    assert!(EwaldReciprocalPlan::new(&diag_cell([1.0, 1.0, 0.0]), 1.0, 1).is_err());
}

#[test]
fn exactly_neutral_charges_have_no_background_at_tiny_screening() {
    let alpha = 2.0_f64.powi(-500);
    let plan = EwaldReciprocalPlan::new(&CELL, alpha, 6).unwrap();
    let result = plan
        .evaluate(&POSITION, &[3.0, -1.0, -1.0, -1.0], 1.0)
        .unwrap();
    let expected = -12.0 * alpha / PI.sqrt();
    assert!((result.energy / expected - 1.0).abs() < 2e-15);
    assert_eq!(result.forces, vec![[0.0; 3]; 4]);
    assert_eq!(result.virial, [[0.0; 3]; 3]);
}

#[test]
fn orthorhombic_phase_wrap_does_not_form_an_overflowing_image_index() {
    let length = 2.0_f64.powi(-400);
    let plan = EwaldReciprocalPlan::new(&diag_cell([length; 3]), 1.0 / length, 5).unwrap();
    let q = [1.0, -1.0];
    let mut pos = [[0.0; 3], [length / 4.0, length / 8.0, 0.0]];
    let expected = plan.evaluate(&pos, &q, 1.0).unwrap();
    pos[0][0] = f64::MAX;
    let wrapped = plan.evaluate(&pos, &q, 1.0).unwrap();
    assert_eq!(expected.energy.to_bits(), wrapped.energy.to_bits());
    assert_eq!(expected.forces, wrapped.forces);
    assert_eq!(expected.virial, wrapped.virial);
}
