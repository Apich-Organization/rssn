use rssn::numerical::graph::Graph;
use rssn::numerical::topology::*;

#[test]
fn test_connected_components() {
    let mut g = Graph::new(4);
    g.add_edge(0, 1, 1.0);
    g.add_edge(2, 3, 1.0);
    let comps = find_connected_components(&g);
    assert_eq!(comps.len(), 2);
}

#[test]
fn test_vietoris_rips() {
    let p1 = [0.0, 0.0];
    let p2 = [0.5, 0.0];
    let p3 = [0.0, 0.5];
    let points = vec![&p1 as &[f64], &p2 as &[f64], &p3 as &[f64]];

    // With epsilon = 0.8, all points are connected to each other
    // Triangle should be formed.
    let simplices = vietoris_rips_complex(&points, 0.8, 2);
    // 0-simplices: [0], [1], [2]
    // 1-simplices: [0,1], [0,2], [1,2]
    // 2-simplices: [0,1,2]
    // Total: 7
    assert_eq!(simplices.len(), 7);
}

#[test]
fn test_betti_numbers() {
    let p1 = [0.0, 0.0];
    let p2 = [1.0, 0.0];
    let p3 = [0.0, 1.0];
    let p4 = [1.0, 1.0];
    let points = vec![&p1 as &[f64], &p2 as &[f64], &p3 as &[f64], &p4 as &[f64]];

    // epsilon = 1.1: Square shape with points connected but no middle connections.
    // Except it's a square, so it might have a hole if we don't triangulate.
    // Vietoris-Rips: if all points in a set are within epsilon, it's a simplex.
    // For a unit square:
    // d(0,1)=1, d(0,2)=1, d(1,3)=1, d(2,3)=1, d(0,3)=sqrt(2)=1.414, d(1,2)=sqrt(2)=1.414

    // epsilon = 1.1:
    // 0-simplices: 0,1,2,3
    // 1-simplices: (0,1), (0,2), (1,3), (2,3)
    // 2-simplices: none
    // B0 = 1 (one component)
    // B1 = 1 (one hole)
    let betti = betti_numbers_at_radius(&points, 1.1, 1);
    assert_eq!(betti[0], 1);
    assert_eq!(betti[1], 1);
}

#[test]
fn test_persistence() {
    let points = vec![
        vec![0.0, 0.0],
        vec![1.0, 0.0],
        vec![0.0, 1.0],
        vec![1.0, 1.0],
    ];
    let diagrams = compute_persistence(&points, 1.5, 15, 1);
    // Dimension 1 hole should be born around 1.0 and die around 1.414
    let d1 = &diagrams[1];
    assert!(
        d1.intervals.len() >= 1,
        "Expected at least one interval in dimension 1, found {:?}",
        d1.intervals
    );
    let found = d1.intervals.iter().any(|interval| {
        interval.birth >= 0.0
            && interval.birth <= 1.2
            && interval.death >= 1.3
            && interval.death <= 1.6
    });
    assert!(
        found,
        "No interval matched birth ~1.0, death ~1.4. Intervals: {:?}",
        d1.intervals
    );
}

#[cfg(test)]
mod proptests {
    use super::*;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn prop_euclidean_distance_symmetric(x1 in -100.0..100.0, y1 in -100.0..100.0, x2 in -100.0..100.0, y2 in -100.0..100.0) {
            let p1 = vec![x1, y1];
            let p2 = vec![x2, y2];
            let d12 = euclidean_distance(&p1, &p2);
            let d21 = euclidean_distance(&p2, &p1);
            prop_assert!((d12 - d21).abs() < 1e-9);
        }

        #[test]
        fn prop_euclidean_distance_triangle_inequality(
            x1 in -100.0..100.0, y1 in -100.0..100.0,
            x2 in -100.0..100.0, y2 in -100.0..100.0,
            x3 in -100.0..100.0, y3 in -100.0..100.0
        ) {
            let p1 = vec![x1, y1];
            let p2 = vec![x2, y2];
            let p3 = vec![x3, y3];
            let d12 = euclidean_distance(&p1, &p2);
            let d23 = euclidean_distance(&p2, &p3);
            let d13 = euclidean_distance(&p1, &p3);
            prop_assert!(d13 <= d12 + d23 + 1e-9);
        }

        #[test]
        fn prop_connected_components_sum(n in 1usize..20, edges in proptest::collection::vec((0usize..20, 0usize..20), 0..30)) {
            let mut g = Graph::new(n);
            for (u, v) in edges {
                if u < n && v < n {
                    g.add_edge(u, v, 1.0);
                }
            }
            let comps = find_connected_components(&g);
            let total_nodes: usize = comps.iter().map(|c| c.len()).sum();
            prop_assert_eq!(total_nodes, n);

            // Each node exactly once
            let mut nodes: Vec<usize> = comps.into_iter().flatten().collect();
            nodes.sort_unstable();
            for i in 0..n {
                prop_assert_eq!(nodes[i], i);
            }
        }
    }
}
