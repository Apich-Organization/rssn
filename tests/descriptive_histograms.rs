//! Descriptive histogram edge membership and bounded parallel counting.

use rssn::kernels::descriptive::histogram;
use rssn::kernels::descriptive::histogram_edges;
use rssn::kernels::numeric_common::NumericError;

fn reference_counts(
    values: &[f64],
    edges: &[f64],
) -> Vec<i64> {
    let mut counts = vec![0; edges.len() - 1];
    for &value in values {
        let bin = edges[1..edges.len() - 1].partition_point(|&edge| edge <= value);
        counts[bin] += 1;
    }
    counts
}

#[test]
fn custom_edges_retain_interior_boundaries_and_both_outer_endpoints() {
    let edges = [-10.0, -1.5, 0.0, 2.0, 12.0];
    let values = [12.0, -10.0, -1.5, -0.0, 0.0, 2.0, 1.0, -9.0];
    assert_eq!(histogram_edges(&values, &edges).unwrap(), [2, 1, 3, 2]);
    assert_eq!(histogram_edges(&values, &[-10.0, 12.0]).unwrap(), [8]);
    assert_eq!(
        histogram_edges(&[-f64::MAX, 0.0, f64::MAX], &[-f64::MAX, 0.0, f64::MAX]).unwrap(),
        [1, 2]
    );
}

#[test]
fn custom_edges_reject_bad_ranges_and_nonfinite_inputs() {
    assert_eq!(histogram_edges(&[], &[0.0, 1.0]), Err(NumericError::Empty));
    for edges in [
        &[][..],
        &[0.0],
        &[1.0, 0.0],
        &[0.0, -0.0],
        &[-1.0, f64::NAN, 1.0],
        &[f64::NEG_INFINITY, 1.0],
        &[0.0, f64::INFINITY],
    ] {
        assert!(matches!(
            histogram_edges(&[0.0], edges),
            Err(NumericError::Domain(_))
        ));
    }
    for value in [-f64::EPSILON, 1.0 + f64::EPSILON, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            histogram_edges(&[0.5, value], &[0.0, 1.0]),
            Err(NumericError::Domain(_))
        ));
    }
}

#[test]
fn uniform_edges_handle_constant_extreme_and_unresolvable_widths() {
    assert_eq!(
        histogram(&[2.0; 3], 3).unwrap(),
        (vec![3, 0, 0], vec![1.5, 2.5, 3.5, 4.5])
    );
    assert_eq!(histogram(&[], 1), Err(NumericError::Empty));
    for bins in [0, usize::MAX] {
        assert!(matches!(
            histogram(&[0.0, 1.0], bins),
            Err(NumericError::Domain(_))
        ));
    }
    assert!(matches!(
        histogram(&[0.0, f64::NAN], 1),
        Err(NumericError::Domain(_))
    ));
    for values in [
        [-f64::MAX, 0.0, f64::MAX],
        [1.0, 1.0_f64.next_up(), 1.0_f64.next_up()],
        [0.0, f64::from_bits(1), f64::from_bits(2)],
    ] {
        let (counts, edges) = histogram(&values, 7).unwrap();
        assert!(edges.iter().all(|edge| edge.is_finite()));
        assert!(edges.is_sorted());
        assert_eq!(edges.first(), values.first());
        assert_eq!(edges.last(), values.last());
        assert_eq!(counts, reference_counts(&values, &edges));
    }
}

#[test]
fn nested_rayon_pools_and_wide_histograms_match_sequential_counts() {
    let data: Vec<_> = (0_u32..(262_144 + 19))
        .map(|index| f64::from(index.wrapping_mul(577) % 1009) - 504.0)
        .collect();
    let explicit = [-504.0, -17.0, -0.5, 1.0, 100.0, 504.0];
    let expected_explicit = reference_counts(&data, &explicit);
    let sequential = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap();
    let workers = rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build()
        .unwrap();
    for bins in [1, 127, 40_001] {
        let expected = sequential.install(|| histogram(&data, bins).unwrap());
        let actual = workers.install(|| histogram(&data, bins).unwrap());
        assert_eq!(actual.0, reference_counts(&data, &actual.1));
        assert_eq!(actual, expected);
    }
    for pool in [&sequential, &workers] {
        assert_eq!(
            pool.install(|| histogram_edges(&data, &explicit).unwrap()),
            expected_explicit
        );
    }
}
