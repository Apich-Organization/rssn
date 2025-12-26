use std::collections::HashSet;

use rssn::symbolic::core::Expr;
use rssn::symbolic::graph::Graph;
use rssn::symbolic::graph_algorithms::*;

#[test]

fn test_graph_basics() {

    let mut g = Graph::new(false);

    g.add_edge(
        &"A",
        &"B",
        Expr::Constant(1.0),
    );

    g.add_edge(
        &"B",
        &"C",
        Expr::Constant(2.0),
    );

    assert_eq!(g.node_count(), 3);

    assert_eq!(
        g.get_edges().len(),
        2
    );

    assert_eq!(
        g.out_degree(
            g.get_node_id(&"B")
                .unwrap()
        ),
        2
    );
}

#[test]

fn test_bfs_dfs() {

    let mut g = Graph::new(false);

    g.add_edge(
        &0,
        &1,
        Expr::Constant(1.0),
    );

    g.add_edge(
        &0,
        &2,
        Expr::Constant(1.0),
    );

    g.add_edge(
        &1,
        &3,
        Expr::Constant(1.0),
    );

    g.add_edge(
        &2,
        &4,
        Expr::Constant(1.0),
    );

    let bfs_order = bfs(&g, 0);

    assert_eq!(bfs_order.len(), 5);

    // BFS from 0 should visit 0 first
    assert_eq!(bfs_order[0], 0);

    let dfs_order = dfs(&g, 0);

    assert_eq!(dfs_order.len(), 5);

    assert_eq!(dfs_order[0], 0);
}

#[test]

fn test_bipartite_matching() {

    let mut g = Graph::new(false);

    // Partition 0: 0, 1, 2
    // Partition 1: 3, 4, 5
    g.add_node(0);

    g.add_node(1);

    g.add_node(2);

    g.add_node(3);

    g.add_node(4);

    g.add_node(5);

    g.add_edge(
        &0,
        &3,
        Expr::Constant(1.0),
    );

    g.add_edge(
        &0,
        &4,
        Expr::Constant(1.0),
    );

    g.add_edge(
        &1,
        &3,
        Expr::Constant(1.0),
    );

    g.add_edge(
        &2,
        &5,
        Expr::Constant(1.0),
    );

    let partition =
        vec![0, 0, 0, 1, 1, 1];

    let matching =
        bipartite_maximum_matching(
            &g,
            &partition,
        );

    // Max matching size should be 3: (0-4), (1-3), (2-5) or similar
    assert_eq!(matching.len(), 3);
}

#[test]

fn test_max_flow() {

    let mut g = Graph::new(true);

    // s=0, t=3
    g.add_edge(
        &0,
        &1,
        Expr::Constant(10.0),
    );

    g.add_edge(
        &0,
        &2,
        Expr::Constant(10.0),
    );

    g.add_edge(
        &1,
        &2,
        Expr::Constant(2.0),
    );

    g.add_edge(
        &1,
        &3,
        Expr::Constant(4.0),
    );

    g.add_edge(
        &2,
        &3,
        Expr::Constant(8.0),
    );

    let flow =
        edmonds_karp_max_flow(&g, 0, 3);

    assert!(
        (flow - 12.0).abs() < 1e-10
    );
}

#[test]

fn test_serialization() {

    let mut g = Graph::new(true);

    g.add_edge(
        &"A".to_string(),
        &"B".to_string(),
        Expr::Constant(1.0),
    );

    let serialized =
        serde_json::to_string(&g)
            .unwrap();

    let deserialized : Graph<String> =
        serde_json::from_str(
            &serialized,
        )
        .unwrap();

    assert_eq!(
        deserialized.node_count(),
        2
    );

    assert!(deserialized.is_directed());
}
