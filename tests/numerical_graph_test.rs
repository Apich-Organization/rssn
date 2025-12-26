use rssn::numerical::graph::{bfs, dijkstra, floyd_warshall, page_rank, Graph};
use std::f64::INFINITY;

#[test]

fn test_graph_creation_and_bfs() {

    let mut graph = Graph::new(5);

    // 0 -> 1 -> 2
    // |    ^
    // v    |
    // 3 -> 4
    graph.add_edge(0, 1, 1.0);

    graph.add_edge(1, 2, 1.0);

    graph.add_edge(0, 3, 1.0);

    graph.add_edge(3, 4, 1.0);

    graph.add_edge(4, 1, 1.0);

    let dist = bfs(&graph, 0);

    assert_eq!(dist[0], 0);

    assert_eq!(dist[1], 1);

    assert_eq!(dist[2], 2);

    assert_eq!(dist[3], 1);

    assert_eq!(dist[4], 2);
}

#[test]

fn test_dijkstra() {

    let mut graph = Graph::new(3);

    graph.add_edge(0, 1, 1.0);

    graph.add_edge(1, 2, 2.0);

    graph.add_edge(0, 2, 4.0); // 1+2 = 3 < 4

    let (dist, prev) = dijkstra(&graph, 0);

    assert_eq!(dist[0], 0.0);

    assert_eq!(dist[1], 1.0);

    assert_eq!(dist[2], 3.0);

    assert_eq!(prev[0], None);

    assert_eq!(prev[1], Some(0));

    assert_eq!(prev[2], Some(1));
}

#[test]

fn test_page_rank() {

    let mut graph = Graph::new(3);

    // 0 -> 1 -> 2 -> 0 circle
    graph.add_edge(0, 1, 1.0);

    graph.add_edge(1, 2, 1.0);

    graph.add_edge(2, 0, 1.0);

    let scores = page_rank(&graph, 0.85, 1e-6, 100);

    // Should be equal due to symmetry
    assert!((scores[0] - scores[1]).abs() < 1e-4);

    assert!((scores[1] - scores[2]).abs() < 1e-4);

    assert!((scores[0] + scores[1] + scores[2] - 1.0).abs() < 1e-6);

    // Star graph: 1->0, 2->0
    let mut graph2 = Graph::new(3);

    graph2.add_edge(1, 0, 1.0);

    graph2.add_edge(2, 0, 1.0);

    let scores2 = page_rank(&graph2, 0.85, 1e-6, 100);

    // Node 0 should have highest score
    assert!(scores2[0] > scores2[1]);

    assert!(scores2[0] > scores2[2]);
}

#[test]

fn test_floyd_warshall() {

    let mut graph = Graph::new(4);

    // 0 -> 1 (1)
    // 1 -> 2 (2)
    // 2 -> 3 (3)
    // 3 -> 0 (10)
    graph.add_edge(0, 1, 1.0);

    graph.add_edge(1, 2, 2.0);

    graph.add_edge(2, 3, 3.0);

    graph.add_edge(3, 0, 10.0);

    let dist = floyd_warshall(&graph);

    let n = 4;

    // 0 to 3: 0->1->2->3 = 1+2+3 = 6
    assert_eq!(dist[0 * n + 3], 6.0);

    // 3 to 0: 3->0 = 10
    assert_eq!(dist[3 * n + 0], 10.0);

    // 0 to 2: 0->1->2 = 3
    assert_eq!(dist[0 * n + 2], 3.0);

    // Diagonal
    assert_eq!(dist[0 * n + 0], 0.0);

    // Unconnected: 1 -> 0 ? 1->2->3->0 = 2+3+10 = 15
    assert_eq!(dist[1 * n + 0], 15.0);
}

#[test]

fn test_connected_components() {

    let mut graph = Graph::new(5);

    // 0-1-2, 3-4
    graph.add_edge(0, 1, 1.0);

    graph.add_edge(1, 0, 1.0);

    graph.add_edge(1, 2, 1.0);

    graph.add_edge(2, 1, 1.0);

    graph.add_edge(3, 4, 1.0);

    graph.add_edge(4, 3, 1.0);

    let comp = rssn::numerical::graph::connected_components(&graph);

    assert_eq!(comp[0], comp[1]);

    assert_eq!(comp[1], comp[2]);

    assert_eq!(comp[3], comp[4]);

    assert_ne!(comp[0], comp[3]);
}

#[test]

fn test_minimum_spanning_tree() {

    let mut graph = Graph::new(4);

    // 0-1 (1), 1-2 (2), 2-3 (3), 0-3 (10)
    // MST should be 0-1, 1-2, 2-3. Total weight = 1+2+3 = 6.
    graph.add_edge(0, 1, 1.0);

    graph.add_edge(1, 0, 1.0);

    graph.add_edge(1, 2, 2.0);

    graph.add_edge(2, 1, 2.0);

    graph.add_edge(2, 3, 3.0);

    graph.add_edge(3, 2, 3.0);

    graph.add_edge(0, 3, 10.0);

    graph.add_edge(3, 0, 10.0);

    let mst = rssn::numerical::graph::minimum_spanning_tree(&graph);

    // Check total weight
    let mut total_weight = 0.0;

    let mut edges_count = 0;

    for u in 0..4 {

        for &(v, w) in mst.adj(u) {

            if u < v {

                // count each edge once
                total_weight += w;

                edges_count += 1;
            }
        }
    }

    assert_eq!(total_weight, 6.0);

    assert_eq!(edges_count, 3);
}

// Property-based tests using proptest
proptest::proptest! {
    /// PageRank scores should always sum to 1.0
    #[test]
    fn prop_page_rank_sums_to_one(n in 2usize..10) {
        let mut graph = Graph::new(n);
        // Create a simple cycle graph
        for i in 0..n {
            graph.add_edge(i, (i + 1) % n, 1.0);
        }
        let scores = page_rank(&graph, 0.85, 1e-9, 200);
        let sum: f64 = scores.iter().sum();
        proptest::prop_assert!((sum - 1.0).abs() < 1e-6, "PageRank sum was {} instead of 1.0", sum);
    }

    /// BFS distance to start node should always be 0
    #[test]
    fn prop_bfs_start_distance_zero(n in 1usize..20, start in 0usize..20) {
        let n = n.max(1);
        let start = start % n;
        let graph = Graph::new(n);
        let dist = bfs(&graph, start);
        proptest::prop_assert_eq!(dist[start], 0);
    }

    /// Dijkstra distance to start node should always be 0
    #[test]
    fn prop_dijkstra_start_distance_zero(n in 1usize..20, start in 0usize..20) {
        let n = n.max(1);
        let start = start % n;
        let graph = Graph::new(n);
        let (dist, _) = dijkstra(&graph, start);
        proptest::prop_assert_eq!(dist[start], 0.0);
    }

    /// Floyd-Warshall diagonal should always be 0
    #[test]
    fn prop_floyd_warshall_diagonal_zero(n in 1usize..10) {
        let graph = Graph::new(n);
        let dist = floyd_warshall(&graph);
        for i in 0..n {
            proptest::prop_assert_eq!(dist[i * n + i], 0.0);
        }
    }

    /// Connected components: each node should have a valid component ID
    #[test]
    fn prop_connected_components_valid(n in 1usize..20) {
        let graph = Graph::new(n);
        let comp = rssn::numerical::graph::connected_components(&graph);
        proptest::prop_assert_eq!(comp.len(), n);
        for i in 0..n {
            // Each isolated node gets its own component
            proptest::prop_assert!(comp[i] < n);
        }
    }

    /// MST of a connected graph with n nodes should have n-1 edges (when counted once)
    #[test]
    fn prop_mst_edge_count(n in 2usize..8) {
        let mut graph = Graph::new(n);
        // Create a fully connected undirected graph
        for i in 0..n {
            for j in (i+1)..n {
                let weight = ((i + j) % 10 + 1) as f64;
                graph.add_edge(i, j, weight);
                graph.add_edge(j, i, weight);
            }
        }
        let mst = rssn::numerical::graph::minimum_spanning_tree(&graph);

        let mut edges_count = 0;
        for u in 0..n {
            for &(v, _) in mst.adj(u) {
                if u < v {
                    edges_count += 1;
                }
            }
        }
        proptest::prop_assert_eq!(edges_count, n - 1, "MST should have n-1 edges");
    }
}
