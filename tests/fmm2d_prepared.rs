//! Reusable open-boundary two-dimensional FMM geometry.

use rssn::sim::fmm::FmmOptions;
use rssn::sim::fmm::PreparedFmm2d;
use rssn::sim::fmm::fmm2d;

fn direct(
    points: &[[f64; 2]],
    values: &[f64],
    eps: f64,
) -> Vec<[f64; 2]> {
    points
        .iter()
        .enumerate()
        .map(|(i, target)| {
            let mut value = [0.0; 2];
            for (j, (source, &strength)) in points.iter().zip(values).enumerate() {
                if i != j {
                    let dx = source[0] - target[0];
                    let dy = source[1] - target[1];
                    let r = dx.hypot(dy).hypot(eps);
                    value[0] += (dx / r) * (strength / r);
                    value[1] += (dy / r) * (strength / r);
                }
            }
            value
        })
        .collect()
}

fn clusters() -> Vec<[f64; 2]> {
    let mut points = Vec::new();
    for (a, b) in [(-0.5, -0.5), (-0.5, 0.5), (0.5, -0.5), (0.5, 0.5)] {
        for i in 0..16 {
            let t = i as f64 * 0.387;
            points.push([a + 0.025 * t.cos(), b + 0.025 * t.sin()]);
        }
    }
    points.extend([[-1.0, -1.0], [1.0, 1.0]]);
    points
}

#[test]
fn near_field_retains_finite_components_across_the_float_range() {
    let src = [1.0, -2.0];
    for scale in [1e-200, 1.0, 1e200] {
        let points = [[0.0, 0.0], [scale, 0.0]];
        for softening in [0.0, 2.0] {
            let opts = FmmOptions {
                eps: softening * scale,
                ..FmmOptions::default()
            };
            let result = fmm2d(&points, &src, &opts).unwrap();
            let normalizer = 1.0 + softening * softening;
            assert!((result[0][0] * scale + 2.0 / normalizer).abs() < 2e-15);
            assert!((result[1][0] * scale + 1.0 / normalizer).abs() < 2e-15);
            assert!(result.iter().all(|[_, y]| *y == 0.0));
        }
    }
    let points = [[-1e308, 0.0], [1e308, 0.0]];
    let result = fmm2d(&points, &src, &FmmOptions::default()).unwrap();
    assert!((result[0][0] * 1e308 + 1.0).abs() < 2e-15);
    assert!((result[1][0] * 1e308 + 0.5).abs() < 2e-15);

    let points = [[0.0, 0.0], [1e-200, 0.0], [1e200, 0.0]];
    let opts = FmmOptions {
        leaf_size: 3,
        ..FmmOptions::default()
    };
    let result = fmm2d(&points, &[1.0, 2.0, 0.0], &opts).unwrap();
    assert!((result[0][0] * 1e-200 - 2.0).abs() < 2e-15);
    assert!((result[1][0] * 1e-200 + 1.0).abs() < 2e-15);
    assert!((result[2][0] * 1e200 + 3.0).abs() < 2e-15);
}

#[test]
fn prepared_expansions_match_direct_sums_when_the_length_unit_changes() {
    let points = clusters();
    let first: Vec<_> = (0..points.len()).map(|i| (i as f64 * 0.22).cos()).collect();
    let second: Vec<_> = (0..points.len()).map(|i| (i as f64 * 0.67).sin()).collect();
    for (scale, eps) in [(1e-150, 0.0), (1.0, 0.0), (1e150, 0.0), (1e150, 1e-4)] {
        let scaled: Vec<_> = points
            .iter()
            .map(|p| [p[0] * scale, p[1] * scale])
            .collect();
        let opts = FmmOptions {
            order: 10,
            leaf_size: 8,
            theta: 0.5,
            eps: eps * scale,
        };
        let prepared = PreparedFmm2d::new(&scaled, &opts).unwrap();
        assert_eq!(prepared.particle_count(), scaled.len());
        for values in [&first, &second] {
            let reference = direct(&points, values, eps);
            let warm = prepared.evaluate(values).unwrap();
            let cold = fmm2d(&scaled, values, &opts).unwrap();
            for ((actual, cold), exact) in warm
                .iter()
                .flatten()
                .zip(cold.iter().flatten())
                .zip(reference.iter().flatten())
            {
                assert_eq!(actual.to_bits(), cold.to_bits());
                let normalized = actual * scale;
                assert!(
                    (normalized - exact).abs() < 2e-4 * (1.0 + exact.abs()),
                    "{scale:e}: {normalized:e} vs {exact:e}"
                );
            }
        }
    }
}

#[test]
fn coincident_points_require_positive_softening_and_inputs_are_checked() {
    let points = [[-0.0, 0.0], [0.0, -0.0]];
    let opts = FmmOptions::default();
    assert!(
        fmm2d(&points, &[1.0, -2.0], &opts)
            .unwrap_err()
            .contains("distinct")
    );
    let plan = PreparedFmm2d::new(&points, &FmmOptions { eps: 1e-250, ..opts }).unwrap();
    assert!(
        plan.evaluate(&[1.0, -2.0])
            .unwrap()
            .iter()
            .flatten()
            .all(|x| *x == 0.0)
    );
    assert!(plan.evaluate(&[1.0]).is_err());
    assert!(plan.evaluate(&[1.0, f64::NAN]).is_err());
    assert!(PreparedFmm2d::new(&[[f64::INFINITY, 0.0]], &opts).is_err());
    assert!(PreparedFmm2d::new(&points, &FmmOptions { eps: -1.0, ..opts }).is_err());
    assert!(
        PreparedFmm2d::new(
            &points,
            &FmmOptions {
                order: usize::MAX,
                ..opts
            }
        )
        .is_err()
    );
}

#[test]
fn softened_quadratic_work_has_a_bounded_interaction_cache() {
    let points: Vec<_> = (0..768)
        .map(|i| {
            let t = i as f64 * 0.019;
            [t.sin(), (t * 1.37).cos()]
        })
        .collect();
    let opts = FmmOptions {
        eps: 8.0,
        leaf_size: 1,
        ..FmmOptions::default()
    };
    let plan = PreparedFmm2d::new(&points, &opts).unwrap();
    assert!(plan.retained_bytes() < 700 * points.len());
    let src: Vec<_> = (0..points.len()).map(|i| (i as f64 * 0.27).cos()).collect();
    let exact = direct(&points, &src, opts.eps);
    let actual = plan.evaluate(&src).unwrap();
    for (a, b) in actual.iter().flatten().zip(exact.iter().flatten()) {
        assert!((a - b).abs() < 1e-13);
    }
}
