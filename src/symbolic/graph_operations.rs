use crate::symbolic::core::Expr;
use crate::symbolic::graph::Graph;
use std::collections::{HashMap, HashSet};
use std::fmt::Debug;
use std::hash::Hash;

/// Creates an induced subgraph from a given set of node labels.
pub fn induced_subgraph<V: Eq + Hash + Clone + Debug>(
    graph: &Graph<V>,
    node_labels: &[V],
) -> Graph<V> {
    let mut sub = Graph::new(graph.is_directed);
    let node_set: HashSet<_> = node_labels.iter().cloned().collect();

    for label in node_labels {
        sub.add_node(label.clone());
    }

    for label in node_labels {
        if let Some(u) = graph.get_node_id(label) {
            if let Some(neighbors) = graph.adj.get(u) {
                for &(v_id, ref weight) in neighbors {
                    let v_label = &graph.nodes[v_id];
                    if node_set.contains(v_label) {
                        sub.add_edge(label, v_label, weight.clone());
                    }
                }
            }
        }
    }
    sub
}

/// Computes the union of two graphs.
pub fn union<V: Eq + Hash + Clone + Debug>(g1: &Graph<V>, g2: &Graph<V>) -> Graph<V> {
    let mut new_graph = g1.clone();
    for (_id, label) in g2.nodes.iter().enumerate() {
        new_graph.add_node(label.clone());
    }
    for (u, v, weight) in g2.get_edges() {
        new_graph.add_edge(&g2.nodes[u], &g2.nodes[v], weight);
    }
    new_graph
}

/// Computes the intersection of two graphs.
pub fn intersection<V: Eq + Hash + Clone + Debug>(g1: &Graph<V>, g2: &Graph<V>) -> Graph<V> {
    let mut new_graph = Graph::new(g1.is_directed && g2.is_directed);
    let g1_nodes: HashSet<_> = g1.nodes.iter().collect();
    let g2_nodes: HashSet<_> = g2.nodes.iter().collect();

    for &node_label in g1_nodes.intersection(&g2_nodes) {
        new_graph.add_node((*node_label).clone());
    }

    for (u, v, weight) in g1.get_edges() {
        let u_label = &g1.nodes[u];
        let v_label = &g1.nodes[v];
        if new_graph.get_node_id(u_label).is_some() && new_graph.get_node_id(v_label).is_some() {
            // Check if this edge also exists in g2
            if let Some(u2_id) = g2.get_node_id(u_label) {
                if let Some(neighbors) = g2.adj.get(u2_id) {
                    if neighbors
                        .iter()
                        .any(|&(v2_id, ref w2)| g2.nodes[v2_id] == *v_label && *w2 == weight)
                    {
                        new_graph.add_edge(u_label, v_label, weight);
                    }
                }
            }
        }
    }
    new_graph
}

/// Computes the Cartesian product of two graphs.
pub fn cartesian_product<V: Eq + Hash + Clone + Debug>(
    g1: &Graph<V>,
    g2: &Graph<V>,
) -> Graph<Expr> {
    let mut new_graph = Graph::new(g1.is_directed || g2.is_directed);
    let mut node_map = HashMap::new();

    // Create new nodes as tuples
    for (u_id, u_label) in g1.nodes.iter().enumerate() {
        for (v_id, v_label) in g2.nodes.iter().enumerate() {
            let new_label = Expr::Tuple(vec![
                Expr::Variable(format!("{:?}", u_label)),
                Expr::Variable(format!("{:?}", v_label)),
            ]);
            node_map.insert((u_id, v_id), new_label.clone());
            new_graph.add_node(new_label);
        }
    }

    // Create edges
    for (u1, v1) in node_map.keys() {
        for (u2, v2) in node_map.keys() {
            let n1 = &node_map[&(*u1, *v1)];
            let n2 = &node_map[&(*u2, *v2)];

            // Case 1: u1 = u2 and (v1, v2) is an edge in g2
            if u1 == u2 && g2.adj[*v1].iter().any(|(n, _)| *n == *v2) {
                new_graph.add_edge(n1, n2, Expr::Constant(1.0));
            }
            // Case 2: v1 = v2 and (u1, u2) is an edge in g1
            if v1 == v2 && g1.adj[*u1].iter().any(|(n, _)| *n == *u2) {
                new_graph.add_edge(n1, n2, Expr::Constant(1.0));
            }
        }
    }
    new_graph
}

/// Computes the Tensor product of two graphs.
pub fn tensor_product<V: Eq + Hash + Clone + Debug>(g1: &Graph<V>, g2: &Graph<V>) -> Graph<Expr> {
    let mut new_graph = Graph::new(g1.is_directed || g2.is_directed);
    let mut node_map = HashMap::new();

    for (u_id, u_label) in g1.nodes.iter().enumerate() {
        for (v_id, v_label) in g2.nodes.iter().enumerate() {
            let new_label = Expr::Tuple(vec![
                Expr::Variable(format!("{:?}", u_label)),
                Expr::Variable(format!("{:?}", v_label)),
            ]);
            node_map.insert((u_id, v_id), new_label.clone());
            new_graph.add_node(new_label);
        }
    }

    for (u1, v1) in node_map.keys() {
        for (u2, v2) in node_map.keys() {
            if g1.adj[*u1].iter().any(|(n, _)| *n == *u2)
                && g2.adj[*v1].iter().any(|(n, _)| *n == *v2)
            {
                let n1 = &node_map[&(*u1, *v1)];
                let n2 = &node_map[&(*u2, *v2)];
                new_graph.add_edge(n1, n2, Expr::Constant(1.0));
            }
        }
    }
    new_graph
}
