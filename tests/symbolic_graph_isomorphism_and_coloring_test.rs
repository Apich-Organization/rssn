use rssn::symbolic::core::Expr;
use rssn::symbolic::graph::Graph;
use rssn::symbolic::graph_isomorphism_and_coloring::*;

#[test]

fn test_isomorphism_heuristic_isomorphic()
 {

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

    g1.add_edge(
        &2,
        &0,
        Expr::Constant(1.0),
    );

    let mut g2 = Graph::new(false);

    g2.add_edge(
        &"A",
        &"B",
        Expr::Constant(1.0),
    );

    g2.add_edge(
        &"B",
        &"C",
        Expr::Constant(1.0),
    );

    g2.add_edge(
        &"C",
        &"A",
        Expr::Constant(1.0),
    );

    assert!(
        are_isomorphic_heuristic(
            &g1, &g2
        )
    );
}

#[test]

fn test_isomorphism_heuristic_non_isomorphic()
 {

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

    // Path graph 0-1-2

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

    // Cycle graph 0-1-2-0

    assert!(
        !are_isomorphic_heuristic(
            &g1, &g2
        )
    );
}

#[test]

fn test_greedy_coloring_triangle() {

    let mut g = Graph::new(false);

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

    let colors = greedy_coloring(&g);

    assert_eq!(colors.len(), 3);

    // Verify valid coloring
    for u in 0 .. g.node_count() {

        for (v, _) in g.neighbors(u) {

            assert_ne!(
                colors[&u],
                colors[&*v]
            );
        }
    }
}

#[test]

fn test_greedy_coloring_bipartite() {

    let mut g = Graph::new(false);

    // Square: 0-1-2-3-0
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
        &3,
        Expr::Constant(1.0),
    );

    g.add_edge(
        &3,
        &0,
        Expr::Constant(1.0),
    );

    let colors = greedy_coloring(&g);

    // Greedy might not be optimal (2), but should be valid.

    // Verify valid coloring
    for u in 0 .. g.node_count() {

        for (v, _) in g.neighbors(u) {

            assert_ne!(
                colors[&u],
                colors[&*v]
            );
        }
    }
}

#[test]

fn test_chromatic_number_exact_k4() {

    let mut g = Graph::new(false);

    // K4
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
        &0,
        &3,
        Expr::Constant(1.0),
    );

    g.add_edge(
        &1,
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

    let chi =
        chromatic_number_exact(&g);

    assert_eq!(chi, 4);
}

#[test]

fn test_chromatic_number_exact_petersen()
 {

    // Petersen graph has chromatic number 3
    let mut g = Graph::new(false);

    // Outer cycle
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
        &3,
        Expr::Constant(1.0),
    );

    g.add_edge(
        &3,
        &4,
        Expr::Constant(1.0),
    );

    g.add_edge(
        &4,
        &0,
        Expr::Constant(1.0),
    );

    // Inner star
    g.add_edge(
        &5,
        &7,
        Expr::Constant(1.0),
    );

    g.add_edge(
        &7,
        &9,
        Expr::Constant(1.0),
    );

    g.add_edge(
        &9,
        &6,
        Expr::Constant(1.0),
    );

    g.add_edge(
        &6,
        &8,
        Expr::Constant(1.0),
    );

    g.add_edge(
        &8,
        &5,
        Expr::Constant(1.0),
    );

    // Spokes
    g.add_edge(
        &0,
        &5,
        Expr::Constant(1.0),
    );

    g.add_edge(
        &1,
        &6,
        Expr::Constant(1.0),
    );

    g.add_edge(
        &2,
        &7,
        Expr::Constant(1.0),
    );

    g.add_edge(
        &3,
        &8,
        Expr::Constant(1.0),
    );

    g.add_edge(
        &4,
        &9,
        Expr::Constant(1.0),
    );

    let chi =
        chromatic_number_exact(&g);

    assert_eq!(chi, 3);
}
