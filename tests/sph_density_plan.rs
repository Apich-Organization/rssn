//! Fixed-geometry SPH density and its periodic and reflecting image terms.

use rssn::sim::sph::KernelKind;
use rssn::sim::sph::Sph;
use rssn::sim::sph::SphBoundary;
use rssn::sim::sph::SphDensityPlan;

fn close(
    expected: f64,
    actual: f64,
) {
    assert!(
        (actual / expected - 1.0).abs() < 2.0e-13,
        "{expected:e} versus {actual:e}"
    );
}

fn direct_mixed(
    pos: &[[f64; 2]],
    masses: &[f64],
    h: &[f64],
    kernel: KernelKind,
) -> Vec<f64> {
    pos.iter()
        .enumerate()
        .map(|(i, &[x, y])| {
            let x = x.rem_euclid(1.0);
            let mut sum = 0.0;
            for (&[px, py], &m) in pos.iter().zip(masses) {
                for sx in [-1.0, 0.0, 1.0] {
                    for image_y in [-py, py, 2.0 - py] {
                        let dx = x - (px.rem_euclid(1.0) + sx);
                        let dy = y - image_y;
                        sum += m * kernel.w(2, dx.hypot(dy), h[i]);
                    }
                }
            }
            sum
        })
        .collect()
}

#[test]
fn reusable_density_matches_independent_periodic_and_mirror_images() {
    let pos = [
        [0.0625, 0.03125],
        [0.9375, 0.0625],
        [0.5, 0.375],
        [0.125, 0.9375],
        [0.875, 0.96875],
    ];
    let bc = [SphBoundary::Periodic, SphBoundary::Reflecting];
    let mut shifted = pos;
    shifted[0][0] += (1_u64 << 40) as f64;
    shifted[1][0] -= (1_u64 << 39) as f64;
    let a = SphDensityPlan::new(&pos, [0.0; 2], [1.0; 2], bc, 0.45).unwrap();
    let b = SphDensityPlan::new(&shifted, [0.0; 2], [1.0; 2], bc, 0.45).unwrap();
    let mass = [0.4, 0.2, 0.1, 0.0, 0.8];
    for kernel in [KernelKind::CubicSpline, KernelKind::WendlandC2] {
        for h in [
            [0.15, 0.22, 0.32, 0.13, 0.18],
            [0.43, 0.35, 0.27, 0.4, 0.33],
        ] {
            let expected = direct_mixed(&shifted, &mass, &h, kernel);
            for actual in [a.evaluate(&mass, &h, kernel), b.evaluate(&mass, &h, kernel)] {
                for (expected, actual) in expected.iter().zip(actual.unwrap()) {
                    close(*expected, actual);
                }
            }
        }
    }
}

#[test]
fn kernel_and_smoothing_scaling_survive_unrepresentable_intermediates() {
    let base = [
        [0.11, 0.13, 0.26],
        [0.29, 0.17, 0.45],
        [0.3, 0.5, 0.6],
        [0.8, 0.9, 0.8],
    ];
    let m = [0.5, 0.4, 0.6, 0.3];
    let h = [0.3, 0.25, 0.28, 0.31];
    let bc = [SphBoundary::Open; 3];
    let plan = SphDensityPlan::new(&base, [0.0; 3], [1.0; 3], bc, 0.35).unwrap();
    let base_sph = Sph::new(
        base.to_vec(),
        m.to_vec(),
        vec![1.0; 4],
        [0.0; 3],
        [1.0; 3],
        bc,
    );
    for sign in [-1, 1] {
        let length = 2.0_f64.powi(sign * 550);
        let mass = 2.0_f64.powi(sign * 800);
        let density = 2.0_f64.powi(-sign * 850);
        let pos = base.map(|p| p.map(|x| x * length));
        let scaled_m = m.map(|value| value * mass);
        let scaled_h = h.map(|value| value * length);
        let scaled_plan =
            SphDensityPlan::new(&pos, [0.0; 3], [length; 3], bc, 0.35 * length).unwrap();
        for kernel in [KernelKind::CubicSpline, KernelKind::WendlandC2] {
            let want = plan.evaluate(&m, &h, kernel).unwrap();
            let actual = scaled_plan.evaluate(&scaled_m, &scaled_h, kernel).unwrap();
            for (want, actual) in want.into_iter().zip(actual) {
                close(want, actual / density);
            }
        }
        let scaled_sph = Sph::new(
            pos.to_vec(),
            scaled_m.to_vec(),
            vec![1.0; 4],
            [0.0; 3],
            [length; 3],
            bc,
        );
        for i in 0..4 {
            close(base_sph.h[i], scaled_sph.h[i] / length);
            close(base_sph.rho[i], scaled_sph.rho[i] / density);
        }
    }
}

#[test]
fn far_open_particle_keeps_its_self_term_when_cell_indices_exceed_i64() {
    let pos = [[0.0], [0.375], [2.0_f64.powi(100)]];
    let h = [0.25; 3];
    let m = [0.8, 0.6, 0.4];
    let plan = SphDensityPlan::new(&pos, [0.0], [1.0], [SphBoundary::Open], h[0]).unwrap();
    let kind = KernelKind::CubicSpline;
    let rho = plan.evaluate(&m, &h, kind).unwrap();
    close(
        m[0] * kind.w(1, 0.0, h[0]) + m[1] * kind.w(1, 0.375, h[0]),
        rho[0],
    );
    close(
        m[1] * kind.w(1, 0.0, h[1]) + m[0] * kind.w(1, 0.375, h[1]),
        rho[1],
    );
    close(m[2] * kind.w(1, 0.0, h[2]), rho[2]);
}

#[test]
fn a_periodic_remainder_rounded_to_the_upper_face_maps_to_the_origin() {
    let solve = |x| {
        SphDensityPlan::new(&[[x], [0.25]], [0.0], [1.0], [SphBoundary::Periodic], 0.2)
            .unwrap()
            .evaluate(&[1.0, 0.75], &[0.2; 2], KernelKind::CubicSpline)
            .unwrap()
    };
    let expected = solve(0.0);
    assert_eq!(solve(-f64::from_bits(1)), expected);
}

#[test]
fn dense_neighbor_fallback_is_deterministic_across_worker_counts() {
    let n = 740;
    let points = vec![[0.0; 2]; n];
    let h = vec![0.25; n];
    let mass: Vec<f64> = (0..n).map(|i| (i % 9 + 1) as f64 / 8.0).collect();
    let kind = KernelKind::WendlandC2;
    let plan =
        SphDensityPlan::new(&points, [0.0; 2], [1.0; 2], [SphBoundary::Open; 2], 0.25).unwrap();
    let run = |workers| {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(workers)
            .build()
            .unwrap();
        pool.install(|| plan.evaluate(&mass, &h, kind).unwrap())
    };
    let scalar = run(1);
    let parallel = run(2);
    let exact = mass.iter().sum::<f64>() * kind.w(2, 0.0, 0.25);
    for (a, b) in scalar.into_iter().zip(parallel) {
        assert_eq!(a.to_bits(), b.to_bits());
        close(exact, a);
    }
}

#[test]
fn checked_plan_handles_empty_fields_and_rejects_invalid_data() {
    let plan = SphDensityPlan::new(&[[0.5]], [0.0], [1.0], [SphBoundary::Periodic], 0.25).unwrap();
    let kind = KernelKind::CubicSpline;
    assert!(plan.evaluate(&[-1.0], &[0.25], kind).is_err());
    assert!(plan.evaluate(&[1.0], &[0.5], kind).is_err());
    assert!(plan.evaluate(&[1.0], &[0.0], kind).is_err());
    assert!(plan.evaluate(&[], &[0.25], kind).is_err());
    assert!(SphDensityPlan::new(&[[0.5]], [1.0], [0.0], [SphBoundary::Reflecting], 0.25).is_err());
    assert!(SphDensityPlan::new(&[[f64::NAN]], [0.0], [1.0], [SphBoundary::Open], 0.25).is_err());
    let empty = SphDensityPlan::new(&[], [0.0], [1.0], [SphBoundary::Open], 0.25).unwrap();
    assert!(empty.evaluate(&[], &[], kind).unwrap().is_empty());
}
