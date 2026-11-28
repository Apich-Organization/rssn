//! Euclidean neighbor search across ordinary and extreme coordinate scales.

use rssn::kernels::geometry::GeometryError;
use rssn::kernels::geometry::KdTree;

#[test]
fn kd_tree_does_not_lose_small_distances_when_squaring() {
    let s = f64::from_bits(1);
    let tree = KdTree::build(vec![
        vec![s, s],
        vec![2.0 * s, 0.0],
        vec![s, 0.0],
        vec![0.0, 0.0],
    ])
    .unwrap();
    let query = [0.0, 0.0];
    assert_eq!(tree.nearest(&query).unwrap(), Some((3, 0.0)));
    assert_eq!(
        tree.k_nearest(&query, usize::MAX).unwrap(),
        [(3, 0.0), (2, s), (0, s), (1, 2.0 * s)],
    );
    assert_eq!(tree.radius_search(&query, 0.0).unwrap(), [(3, 0.0)]);
    // sqrt(2) * s rounds to s on output. It still lies outside radius s.
    assert_eq!(tree.radius_search(&query, s).unwrap(), [(3, 0.0), (2, s)]);
}

#[test]
fn kd_tree_compares_large_distances_before_rounding_the_output() {
    let tree = KdTree::build(vec![
        vec![2e200, 0.0],
        vec![1e200, 1e200],
        vec![-3e200, 0.0],
        vec![0.0, 1e200],
    ])
    .unwrap();
    let query = [0.0, 0.0];
    let all = tree.k_nearest(&query, 4).unwrap();
    assert_eq!(
        all.iter().map(|entry| entry.0).collect::<Vec<_>>(),
        [3, 1, 0, 2]
    );
    assert!(all.iter().all(|entry| entry.1.is_finite()));
    assert_eq!(tree.nearest(&query).unwrap(), Some(all[0]));
    assert_eq!(tree.radius_search(&query, 1.5e200).unwrap(), all[..2]);

    let max = f64::MAX;
    let overflowing =
        KdTree::build(vec![vec![-max], vec![-0.75 * max], vec![-0.25 * max]]).unwrap();
    let unbounded = [(2, f64::INFINITY), (1, f64::INFINITY), (0, f64::INFINITY)];
    assert_eq!(overflowing.nearest(&[max]).unwrap(), Some(unbounded[0]));
    assert_eq!(overflowing.k_nearest(&[max], 3).unwrap(), unbounded);
    assert!(overflowing.radius_search(&[max], max).unwrap().is_empty());
    assert_eq!(
        overflowing.radius_search(&[max], f64::INFINITY).unwrap(),
        unbounded
    );
}

#[test]
fn kd_search_and_split_pruning_agree_across_binary_scales() {
    let scales = [
        f64::from_bits(1),
        f64::MIN_POSITIVE,
        f64::from_bits((1023 - 600) << 52),
        f64::from_bits((1023 - 512) << 52),
        f64::from_bits((1023 - 450) << 52),
        1.0,
        f64::from_bits((1023 + 450) << 52),
        f64::from_bits((1023 + 512) << 52),
        f64::from_bits((1023 + 1019) << 52),
    ];
    let indices = |entries: Vec<(usize, f64)>| {
        entries
            .into_iter()
            .map(|(index, _)| index)
            .collect::<Vec<_>>()
    };
    for scale in scales {
        let points = [[-4.0, 0.0], [3.0, 4.0], [1.0, -2.0], [5.0, 4.0], [0.0, 0.0]];
        let points = points
            .into_iter()
            .map(|point| vec![point[0] * scale, point[1] * scale]);
        let tree = KdTree::build(points.collect()).unwrap();
        let query = [2.0 * scale, -2.0 * scale];
        assert_eq!(tree.nearest(&query).unwrap(), Some((2, scale)));
        assert_eq!(indices(tree.k_nearest(&query, 5).unwrap()), [2, 4, 1, 0, 3]);
        assert_eq!(
            indices(tree.radius_search(&query, 3.0 * scale).unwrap()),
            [2, 4]
        );
    }
}

#[test]
fn kd_tree_nearest_batch_preserves_order_ties_and_validation() {
    let tree = KdTree::build(vec![
        vec![-1.0, 0.0],
        vec![1.0, 0.0],
        vec![2.0, -1.0],
        vec![1.0, 0.0],
    ])
    .unwrap();
    let inputs: Vec<_> = (0..512)
        .map(|i| vec![f64::from(i % 5) - 2.0, f64::from(i % 9) - 4.0])
        .collect();
    let expected: Vec<_> = inputs
        .iter()
        .map(|query| tree.nearest(query).unwrap())
        .collect();
    assert_eq!(tree.nearest_batch(&inputs).unwrap(), expected);
    for (query, nearest) in inputs.iter().zip(expected) {
        assert_eq!(tree.k_nearest(query, 1).unwrap().first().copied(), nearest);
    }
    let mixed = [vec![1.0, 0.0], vec![1e200, -1e200]];
    assert_eq!(
        tree.nearest_batch(&mixed).unwrap(),
        mixed
            .iter()
            .map(|query| tree.nearest(query).unwrap())
            .collect::<Vec<_>>(),
    );
    assert!(tree.nearest_batch(&[]).unwrap().is_empty());
    assert!(matches!(
        tree.nearest_batch(&[vec![0.0, 0.0], vec![1.0], vec![f64::NAN; 2]]),
        Err(GeometryError::ShapeMismatch(_)),
    ));
    assert!(matches!(
        tree.k_nearest(&[f64::NAN], 0),
        Err(GeometryError::Domain(_))
    ));
    assert!(tree.radius_search(&[0.0, 0.0], -1.0).unwrap().is_empty());
}
