use rssn::symbolic::core::Expr;
use rssn::symbolic::graph::Graph;
use rssn::symbolic::graph_algorithms::*;

#[test]

fn test_dfs_traversal() {

    let mut g = Graph::new(false);

    g.add_edge(
        &"A",
        &"B",
        Expr::Constant(1.0),
    );

    g.add_edge(
        &"B",
        &"C",
        Expr::Constant(1.0),
    );

    g.add_edge(
        &"C",
        &"D",
        Expr::Constant(1.0),
    );

    g.add_edge(
        &"A",
        &"D",
        Expr::Constant(1.0),
    );

    let result = dfs(&g, 0);

    assert_eq!(result.len(), 4);

    assert_eq!(result[0], 0); // Starts at node 0
}

#[test]

fn test_bfs_traversal() {

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

    let result = bfs(&g, 0);

    assert_eq!(result.len(), 5);

    assert_eq!(result[0], 0);

    // BFS should visit level by level
    assert!(result[1] == 1 || result[1] == 2);
}

#[test]

fn test_connected_components() {

    let mut g = Graph::new(false);

    // Component 1: 0-1-2
    g.add_edge(
        &0,
        &1,
        Expr::Constant(1.0),
    );

    g.add_edge(
        &1,
        &2,
        Expr::Constant(1.0),
    );

    // Component 2: 3-4
    g.add_edge(
        &3,
        &4,
        Expr::Constant(1.0),
    );

    // Isolated node: 5
    g.add_node(5);

    let components = connected_components(&g);

    assert_eq!(components.len(), 3);
}

#[test]

fn test_is_connected() {

    let mut g1 = Graph::new(false);

    g1.add_edge(
        &0,
        &1,
        Expr::Constant(1.0),
    );

    g1.add_edge(
        &1,
        &2,
        Expr::Constant(1.0),
    );

    assert!(is_connected(&g1));

    let mut g2 = Graph::new(false);

    g2.add_edge(
        &0,
        &1,
        Expr::Constant(1.0),
    );

    g2.add_node(2); // Isolated node
    assert!(!is_connected(&g2));
}

#[test]

fn test_strongly_connected_components() {

    let mut g = Graph::new(true);

    // SCC 1: 0 <-> 1
    g.add_edge(
        &0,
        &1,
        Expr::Constant(1.0),
    );

    g.add_edge(
        &1,
        &0,
        Expr::Constant(1.0),
    );

    // SCC 2: 2 <-> 3
    g.add_edge(
        &2,
        &3,
        Expr::Constant(1.0),
    );

    g.add_edge(
        &3,
        &2,
        Expr::Constant(1.0),
    );

    // Connection between SCCs
    g.add_edge(
        &1,
        &2,
        Expr::Constant(1.0),
    );

    let sccs = strongly_connected_components(&g);

    assert!(sccs.len() >= 2);
}

#[test]

fn test_has_cycle_undirected() {

    let mut g1 = Graph::new(false);

    g1.add_edge(
        &0,
        &1,
        Expr::Constant(1.0),
    );

    g1.add_edge(
        &1,
        &2,
        Expr::Constant(1.0),
    );

    assert!(!has_cycle(&g1)); // Tree, no cycle

    let mut g2 = Graph::new(false);

    g2.add_edge(
        &0,
        &1,
        Expr::Constant(1.0),
    );

    g2.add_edge(
        &1,
        &2,
        Expr::Constant(1.0),
    );

    g2.add_edge(
        &2,
        &0,
        Expr::Constant(1.0),
    );

    assert!(has_cycle(&g2)); // Triangle, has cycle
}

#[test]

fn test_has_cycle_directed() {

    let mut g1 = Graph::new(true);

    g1.add_edge(
        &0,
        &1,
        Expr::Constant(1.0),
    );

    g1.add_edge(
        &1,
        &2,
        Expr::Constant(1.0),
    );

    assert!(!has_cycle(&g1)); // DAG, no cycle

    let mut g2 = Graph::new(true);

    g2.add_edge(
        &0,
        &1,
        Expr::Constant(1.0),
    );

    g2.add_edge(
        &1,
        &2,
        Expr::Constant(1.0),
    );

    g2.add_edge(
        &2,
        &0,
        Expr::Constant(1.0),
    );

    assert!(has_cycle(&g2)); // Cycle
}

#[test]

fn test_bridges_and_articulation_points() {

    let mut g = Graph::new(false);

    // 0-1-2 where 1 is an articulation point and 0-1, 1-2 are bridges
    g.add_edge(
        &0,
        &1,
        Expr::Constant(1.0),
    );

    g.add_edge(
        &1,
        &2,
        Expr::Constant(1.0),
    );

    let (bridges, aps) = find_bridges_and_articulation_points(&g);

    assert_eq!(bridges.len(), 2);

    assert_eq!(aps.len(), 1);

    assert!(aps.contains(&1));
}

#[test]

fn test_kruskal_mst() {

    let mut g = Graph::new(false);

    g.add_edge(
        &0,
        &1,
        Expr::Constant(1.0),
    );

    g.add_edge(
        &0,
        &2,
        Expr::Constant(4.0),
    );

    g.add_edge(
        &1,
        &2,
        Expr::Constant(2.0),
    );

    g.add_edge(
        &1,
        &3,
        Expr::Constant(5.0),
    );

    g.add_edge(
        &2,
        &3,
        Expr::Constant(3.0),
    );

    let mst = kruskal_mst(&g);

    assert_eq!(mst.len(), 3); // MST has n-1 edges for n nodes

    // Check total weight is minimal (1 + 2 + 3 = 6)
    let total_weight: f64 = mst
        .iter()
        .filter_map(|(_, _, w)| {

            match w {
                Expr::Constant(v) => Some(*v),
                _ => None,
            }
        })
        .sum();

    assert!((total_weight - 6.0).abs() < 1e-10);
}

#[test]

fn test_edmonds_karp_max_flow() {

    let mut g = Graph::new(true);

    // Classic max flow example
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
        &1,
        &4,
        Expr::Constant(8.0),
    );

    g.add_edge(
        &2,
        &4,
        Expr::Constant(9.0),
    );

    g.add_edge(
        &3,
        &5,
        Expr::Constant(10.0),
    );

    g.add_edge(
        &4,
        &3,
        Expr::Constant(6.0),
    );

    g.add_edge(
        &4,
        &5,
        Expr::Constant(10.0),
    );

    let flow = edmonds_karp_max_flow(&g, 0, 5);

    assert!(flow > 0.0);
}

#[test]

fn test_dinic_max_flow() {

    let mut g = Graph::new(true);

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
        &3,
        Expr::Constant(4.0),
    );

    g.add_edge(
        &2,
        &3,
        Expr::Constant(8.0),
    );

    let flow = dinic_max_flow(&g, 0, 3);

    assert!((flow - 12.0).abs() < 1e-10);
}

#[test]

fn test_bellman_ford_symbolic() {

    let mut g = Graph::new(true);

    g.add_edge(
        &0,
        &1,
        Expr::Constant(5.0),
    );

    g.add_edge(
        &0,
        &2,
        Expr::Constant(3.0),
    );

    g.add_edge(
        &1,
        &3,
        Expr::Constant(2.0),
    );

    g.add_edge(
        &2,
        &1,
        Expr::Constant(1.0),
    );

    g.add_edge(
        &2,
        &3,
        Expr::Constant(6.0),
    );

    let result = bellman_ford(&g, 0);

    assert!(result.is_ok());

    let (distances, _) = result.unwrap();

    assert_eq!(distances.len(), 4);

    // Distance to node 0 should be 0
    match &distances[&0] {
        Expr::Constant(v) => assert_eq!(*v, 0.0),
        _ => panic!("Expected constant"),
    }
}

#[test]

fn test_bellman_ford_negative_cycle() {

    let mut g = Graph::new(true);

    g.add_edge(
        &0,
        &1,
        Expr::Constant(1.0),
    );

    g.add_edge(
        &1,
        &2,
        Expr::Constant(1.0),
    );

    g.add_edge(
        &2,
        &0,
        Expr::Constant(-5.0),
    ); // Creates negative cycle

    let result = bellman_ford(&g, 0);

    assert!(result.is_err());
}

#[test]

fn test_is_bipartite() {

    let mut g1 = Graph::new(false);

    // Bipartite graph: 0-1, 0-2, 1-3, 2-3
    g1.add_edge(
        &0,
        &1,
        Expr::Constant(1.0),
    );

    g1.add_edge(
        &0,
        &2,
        Expr::Constant(1.0),
    );

    g1.add_edge(
        &1,
        &3,
        Expr::Constant(1.0),
    );

    g1.add_edge(
        &2,
        &3,
        Expr::Constant(1.0),
    );

    let partition = is_bipartite(&g1);

    assert!(partition.is_some());

    let mut g2 = Graph::new(false);

    // Not bipartite: triangle
    g2.add_edge(
        &0,
        &1,
        Expr::Constant(1.0),
    );

    g2.add_edge(
        &1,
        &2,
        Expr::Constant(1.0),
    );

    g2.add_edge(
        &2,
        &0,
        Expr::Constant(1.0),
    );

    assert!(is_bipartite(&g2).is_none());
}

#[test]

fn test_bipartite_maximum_matching() {

    let mut g = Graph::new(false);

    // Left: 0, 1, 2; Right: 3, 4, 5
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
        &1,
        &5,
        Expr::Constant(1.0),
    );

    g.add_edge(
        &2,
        &5,
        Expr::Constant(1.0),
    );

    let partition = vec![0, 0, 0, 1, 1, 1];

    let matching = bipartite_maximum_matching(&g, &partition);

    assert!(matching.len() >= 2);
}

#[test]

fn test_topological_sort() {

    let mut g = Graph::new(true);

    // DAG: 0 -> 1 -> 3, 0 -> 2 -> 3
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
        &3,
        Expr::Constant(1.0),
    );

    let sorted = topological_sort(&g);

    assert!(sorted.is_some());

    let order = sorted.unwrap();

    assert_eq!(order.len(), 4);

    // 0 should come before 1, 2, 3
    let pos_0 = order
        .iter()
        .position(|&x| x == 0)
        .unwrap();

    let pos_3 = order
        .iter()
        .position(|&x| x == 3)
        .unwrap();

    assert!(pos_0 < pos_3);
}

#[test]

fn test_topological_sort_with_cycle() {

    let mut g = Graph::new(true);

    // Graph with cycle
    g.add_edge(
        &0,
        &1,
        Expr::Constant(1.0),
    );

    g.add_edge(
        &1,
        &2,
        Expr::Constant(1.0),
    );

    g.add_edge(
        &2,
        &0,
        Expr::Constant(1.0),
    );

    let sorted = topological_sort(&g);

    assert!(sorted.is_none()); // Should fail due to cycle
}

#[test]

fn test_symbolic_weights() {

    let mut g = Graph::new(false);

    // Use symbolic expressions as weights
    let x = Expr::Variable("x".to_string());

    let y = Expr::Variable("y".to_string());

    g.add_edge(&0, &1, x.clone());

    g.add_edge(&1, &2, y.clone());

    // Verify edges preserve symbolic structure
    let edges = g.get_edges();

    assert_eq!(edges.len(), 2);

    match &edges[0].2 {
        Expr::Variable(name) => assert!(name == "x" || name == "y"),
        _ => panic!("Expected variable"),
    }
}

#[test]

fn test_min_cost_max_flow() {

    let mut g = Graph::new(true);

    // Edge weights are (capacity, cost) tuples
    g.add_edge(
        &0,
        &1,
        Expr::Tuple(vec![
            Expr::Constant(10.0),
            Expr::Constant(2.0),
        ]),
    );

    g.add_edge(
        &0,
        &2,
        Expr::Tuple(vec![
            Expr::Constant(10.0),
            Expr::Constant(4.0),
        ]),
    );

    g.add_edge(
        &1,
        &3,
        Expr::Tuple(vec![
            Expr::Constant(5.0),
            Expr::Constant(1.0),
        ]),
    );

    g.add_edge(
        &2,
        &3,
        Expr::Tuple(vec![
            Expr::Constant(5.0),
            Expr::Constant(3.0),
        ]),
    );

    let (flow, cost) = min_cost_max_flow(&g, 0, 3);

    assert!(flow > 0.0);

    assert!(cost > 0.0);
}
