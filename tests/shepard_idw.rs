//! Scale-stable Shepard inverse-distance interpolation.

use rssn::kernels::interp_ext::InterpError;
use rssn::kernels::interp_ext::shepard_idw;

fn close(
    actual: f64,
    expected: f64,
    tolerance: f64,
) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} vs {expected}"
    );
}

#[test]
fn distinct_small_distances_and_large_distances_give_the_same_weights() {
    for scale in [1e-200, 1.0, 1e200] {
        let left = [0.0, 0.0];
        let right = [2.0 * scale, 0.0];
        let points: [&[f64]; 2] = [&left, &right];
        close(
            shepard_idw(&points, &[0.0, 10.0], &[0.5 * scale, 0.0], 2.0).unwrap(),
            1.0,
            3e-15,
        );
        close(
            shepard_idw(&points, &[0.0, 10.0], &right, 2.0).unwrap(),
            10.0,
            0.0,
        );
    }
    close(
        shepard_idw(&[&[1e308], &[0.0]], &[10.0, 0.0], &[-1e308], 2.0).unwrap(),
        2.0,
        1e-14,
    );
}

#[test]
fn weights_and_ordinates_remain_finite_at_scale_extremes() {
    let points: [&[f64]; 2] = [&[-0.1], &[0.1]];
    let large = f64::MAX;
    close(
        shepard_idw(&points, &[large, large], &[0.0], 3.0).unwrap(),
        large,
        0.0,
    );
    close(
        shepard_idw(&points, &[large, -large], &[0.0], 3.0).unwrap(),
        0.0,
        0.0,
    );

    let tiny = f64::from_bits(1);
    let far = 1e308_f64;
    let power = 1e-6;
    let weight = (power * (tiny.ln() - far.ln())).exp();
    let expected = (2.0 + 4.0 * weight) / (1.0 + weight);
    close(
        shepard_idw(&[&[tiny], &[far]], &[2.0, 4.0], &[0.0], power).unwrap(),
        expected,
        3e-15,
    );
    close(
        shepard_idw(&[&[1.0], &[2.0]], &[2.0, 4.0], &[0.0], 1e308).unwrap(),
        2.0,
        0.0,
    );
}

#[test]
fn exact_matches_keep_first_value_after_validating_all_samples() {
    let point = [2e-250, -3e-250];
    close(
        shepard_idw(&[&point, &point], &[2.0, 4.0], &point, 2.0).unwrap(),
        2.0,
        0.0,
    );
    assert!(matches!(
        shepard_idw(&[&point, &[1.0]], &[2.0, 3.0], &point, 2.0),
        Err(InterpError::ShapeMismatch(_))
    ));
    for power in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            shepard_idw(&[&point], &[2.0], &point, power),
            Err(InterpError::Domain(_))
        ));
    }
    assert!(matches!(
        shepard_idw(&[&point], &[f64::NAN], &point, 2.0),
        Err(InterpError::Domain(_))
    ));
    assert!(matches!(
        shepard_idw(&[&[f64::INFINITY, 0.0]], &[2.0], &point, 2.0),
        Err(InterpError::Domain(_))
    ));
    assert!(matches!(
        shepard_idw(&[&point], &[2.0], &[0.0, f64::NAN], 2.0),
        Err(InterpError::Domain(_))
    ));
}

#[test]
fn block_reduction_matches_an_independent_weighted_residual() {
    for (dimension, count) in [(1, 5), (7, 257), (32, 5200)] {
        let data: Vec<Vec<f64>> = (0..count)
            .map(|i| {
                (0..dimension)
                    .map(|axis| f64::from((i * (axis + 17) + axis * 21) % 8191) / 8191.0)
                    .collect()
            })
            .collect();
        let points: Vec<_> = data.iter().map(Vec::as_slice).collect();
        let values: Vec<_> = (0..count).map(|i| (f64::from(i) * 0.71).sin()).collect();
        let query: Vec<_> = (0..dimension)
            .map(|i| f64::from(i) * 0.041 + 0.128)
            .collect();
        let power = 3.25;
        let (mut numerator, mut denominator) = (0.0, 0.0);
        for (point, &value) in points.iter().zip(&values) {
            let distance = point
                .iter()
                .zip(&query)
                .fold(0.0_f64, |norm, (&a, &b)| norm.hypot(a - b));
            let weight = 1.0 / distance.powf(power);
            numerator += weight * value;
            denominator += weight;
        }
        let evaluate = |threads| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(|| shepard_idw(&points, &values, &query, power).unwrap())
        };
        let single = evaluate(1);
        close(single, numerator / denominator, 2e-13);
        assert_eq!(single.to_bits(), evaluate(4).to_bits());
        if count > 4096 {
            let tiny: Vec<Vec<_>> = data
                .iter()
                .map(|point| point.iter().map(|coordinate| coordinate * 1e-200).collect())
                .collect();
            let tiny_points: Vec<_> = tiny.iter().map(Vec::as_slice).collect();
            let tiny_query: Vec<_> = query.iter().map(|coordinate| coordinate * 1e-200).collect();
            let evaluate = |threads| {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .unwrap()
                    .install(|| shepard_idw(&tiny_points, &values, &tiny_query, power).unwrap())
            };
            let scaled = evaluate(1);
            close(scaled, single, 2e-13);
            assert_eq!(scaled.to_bits(), evaluate(4).to_bits());
        }
    }
}
