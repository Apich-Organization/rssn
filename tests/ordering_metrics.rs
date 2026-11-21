//! Edge cases for sparse pattern traversal and ordering diagnostics.

use rssn::kernels::numeric_common::NumericError;
use rssn::kernels::ordering::OrderingMetrics;
use rssn::kernels::ordering::amd_order;
use rssn::kernels::ordering::ordering_metrics;

#[test]
fn amd_preserves_isolate_order_and_ignores_repeated_entries() {
    assert!(amd_order(&[]).is_empty());
    assert_eq!(
        amd_order(&vec![Vec::new(); 40]),
        (0..40).rev().collect::<Vec<_>>(),
    );

    let mut edges = vec![Vec::new(); 32];
    for range in [0..18, 20..24] {
        for v in range.start + 1..range.end {
            edges[v].push(v - 1);
            edges[v - 1].push(v);
        }
    }
    let simple = amd_order(&edges);
    for (v, neighbors) in edges.iter_mut().enumerate() {
        let unique = neighbors.clone();
        neighbors.push(v);
        for _ in 0..64 {
            neighbors.extend_from_slice(&unique);
        }
    }
    assert_eq!(amd_order(&edges), simple);
    assert_eq!(simple[..10], [31, 30, 29, 28, 27, 26, 25, 24, 19, 18]);
}

#[test]
fn bandwidth_and_profile_use_new_positions_and_count_edges_once() {
    let adj = [
        vec![0, 1, 1],
        vec![2, 0],
        vec![3, 3],
        vec![3],
        vec![3],
        vec![],
    ];
    let natural = [0, 1, 2, 3, 4, 5];
    let shuffled = [0, 2, 4, 1, 3, 5];
    assert_eq!(
        ordering_metrics(&adj, &natural).unwrap(),
        OrderingMetrics {
            half_bandwidth: 1,
            profile: 4,
        },
    );
    assert_eq!(
        ordering_metrics(&adj, &shuffled).unwrap(),
        OrderingMetrics {
            half_bandwidth: 3,
            profile: 6,
        },
    );
}

#[test]
fn ordering_metrics_validate_shape_and_indices_before_computing_widths() {
    let adj = [Vec::new(), Vec::new()];
    assert_eq!(
        ordering_metrics(&[], &[]).unwrap(),
        OrderingMetrics::default(),
    );
    assert_eq!(
        ordering_metrics(&adj, &[1, 0]).unwrap(),
        OrderingMetrics::default(),
    );
    assert!(matches!(
        ordering_metrics(&adj, &[1]),
        Err(NumericError::ShapeMismatch(_)),
    ));
    for permutation in [[0, 0], [0, 2]] {
        assert!(matches!(
            ordering_metrics(&adj, &permutation),
            Err(NumericError::Domain(_)),
        ));
    }
    assert!(matches!(
        ordering_metrics(&[vec![usize::MAX]], &[0]),
        Err(NumericError::Domain(_)),
    ));
}

#[test]
fn sparse_profile_does_not_truncate_to_a_narrow_integer() {
    let half = 65_536usize;
    let mut adj = vec![Vec::new(); 2 * half];
    for (i, row) in adj.iter_mut().take(half).enumerate() {
        row.push(i + half);
    }
    let natural: Vec<_> = (0..adj.len()).collect();
    let result = ordering_metrics(&adj, &natural);
    if let Some(profile) = half.checked_mul(half) {
        assert_eq!(
            result.unwrap(),
            OrderingMetrics {
                half_bandwidth: half,
                profile,
            },
        );
    } else {
        assert!(matches!(result, Err(NumericError::Domain(_))));
    }
}
