use rssn::symbolic::core::Expr;
use rssn::symbolic::graph::Graph;
use rssn::symbolic::graph_operations::*;

#[test]
fn test_induced_subgraph() {
    let mut g = Graph::new(false);
    g.add_edge(&"A", &"B", Expr::Constant(1.0));
    g.add_edge(&"B", &"C", Expr::Constant(1.0));
    g.add_edge(&"C", &"A", Expr::Constant(1.0));
    g.add_edge(&"C", &"D", Expr::Constant(1.0));

    let sub = induced_subgraph(&g, &["A", "B", "C"]);
    assert_eq!(sub.node_count(), 3);
    assert_eq!(sub.get_edges().len(), 3); // A-B, B-C, C-A
    assert!(sub.get_node_id(&"D").is_none());
}

#[test]
fn test_union() {
    let mut g1 = Graph::new(false);
    g1.add_edge(&"A", &"B", Expr::Constant(1.0));

    let mut g2 = Graph::new(false);
    g2.add_edge(&"B", &"C", Expr::Constant(1.0));

    let u = union(&g1, &g2);
    assert_eq!(u.node_count(), 3); // A, B, C
    assert_eq!(u.get_edges().len(), 2); // A-B, B-C
}

#[test]
fn test_intersection() {
    let mut g1 = Graph::new(false);
    g1.add_edge(&"A", &"B", Expr::Constant(1.0));
    g1.add_edge(&"B", &"C", Expr::Constant(1.0));

    let mut g2 = Graph::new(false);
    g2.add_edge(&"B", &"C", Expr::Constant(1.0));
    g2.add_edge(&"C", &"D", Expr::Constant(1.0));

    let i = intersection(&g1, &g2);
    assert_eq!(i.node_count(), 2); // B, C (intersection of nodes)
    assert_eq!(i.get_edges().len(), 1); // B-C (intersection of edges)
}

#[test]
fn test_cartesian_product() {
    let mut g1 = Graph::new(false);
    g1.add_edge(&"A", &"B", Expr::Constant(1.0)); // Path graph P2

    let mut g2 = Graph::new(false);
    g2.add_edge(&"1", &"2", Expr::Constant(1.0)); // Path graph P2

    // P2 x P2 should be C4 (Cycle graph with 4 nodes)
    let prod = cartesian_product(&g1, &g2);
    assert_eq!(prod.node_count(), 4);
    assert_eq!(prod.get_edges().len(), 4);
}

#[test]
fn test_complement() {
    let mut g = Graph::new(false);
    // Triangle graph: A-B, B-C, C-A
    g.add_edge(&"A", &"B", Expr::Constant(1.0));
    g.add_edge(&"B", &"C", Expr::Constant(1.0));
    g.add_edge(&"C", &"A", Expr::Constant(1.0));

    // Complement of K3 is empty graph (with 3 nodes)
    let comp = complement(&g);
    assert_eq!(comp.node_count(), 3);
    assert_eq!(comp.get_edges().len(), 0);

    // Path graph P3: A-B-C
    let mut g2 = Graph::new(false);
    g2.add_edge(&"A", &"B", Expr::Constant(1.0));
    g2.add_edge(&"B", &"C", Expr::Constant(1.0));

    // Complement of P3 should have edge A-C
    let comp2 = complement(&g2);
    assert_eq!(comp2.get_edges().len(), 1);
    // Check if edge is A-C
    // Since we don't know exact indices, we check if there is an edge.
}

#[test]
fn test_disjoint_union() {
    let mut g1 = Graph::new(false);
    g1.add_edge(&"A", &"B", Expr::Constant(1.0));

    let mut g2 = Graph::new(false);
    g2.add_edge(&"A", &"B", Expr::Constant(1.0));

    let du = disjoint_union(&g1, &g2);
    assert_eq!(du.node_count(), 4); // A1, B1, A2, B2
    assert_eq!(du.get_edges().len(), 2); // A1-B1, A2-B2
}

#[test]
fn test_join() {
    let mut g1 = Graph::new(false);
    g1.add_node("A");

    let mut g2 = Graph::new(false);
    g2.add_node("B");

    let j = join(&g1, &g2);
    assert_eq!(j.node_count(), 2);
    assert_eq!(j.get_edges().len(), 1); // A-B
}
