use std::collections::HashMap;
use std::collections::HashSet;
use std::fmt::Debug;
use std::hash::Hash;

use crate::symbolic::core::Expr;
use crate::symbolic::graph::Graph;

/// Trait to convert a value to a symbolic expression.
/// This allows preserving the structure of node labels when creating product graphs.

pub trait ToExpr {
    fn to_expr(&self) -> Expr;
}

impl ToExpr for String {
    fn to_expr(&self) -> Expr {

        Expr::Variable(self.clone())
    }
}

impl ToExpr for &str {
    fn to_expr(&self) -> Expr {

        Expr::Variable((*self).to_string())
    }
}

impl ToExpr for Expr {
    fn to_expr(&self) -> Expr {

        self.clone()
    }
}

impl ToExpr for usize {
    fn to_expr(&self) -> Expr {

        Expr::Constant(*self as f64)
    }
}

impl ToExpr for i32 {
    fn to_expr(&self) -> Expr {

        Expr::Constant(f64::from(*self))
    }
}

/// Creates an induced subgraph from a given set of node labels.
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::core::Expr;
/// use rssn::symbolic::graph::Graph;
/// use rssn::symbolic::graph_operations::induced_subgraph;
///
/// let mut g = Graph::new(false);
///
/// g.add_edge(
///     &"A",
///     &"B",
///     Expr::Constant(1.0),
/// );
///
/// g.add_edge(
///     &"B",
///     &"C",
///     Expr::Constant(1.0),
/// );
///
/// let sub = induced_subgraph(&g, &["A", "B"]);
///
/// assert_eq!(sub.node_count(), 2);
/// ```

pub fn induced_subgraph<V : Eq + Hash + Clone + Debug>(
    graph : &Graph<V>,
    node_labels : &[V],
) -> Graph<V> {

    let mut sub = Graph::new(graph.is_directed);

    let node_set : HashSet<_> = node_labels
        .iter()
        .cloned()
        .collect();

    for label in node_labels {

        sub.add_node(label.clone());
    }

    for label in node_labels {

        if let Some(u) = graph.get_node_id(label) {

            if let Some(neighbors) = graph.adj.get(u) {

                for &(v_id, ref weight) in neighbors {

                    // For undirected graphs, we only add if u <= v_id to avoid duplicates,
                    // because add_edge adds both directions.
                    if !graph.is_directed && u > v_id {

                        continue;
                    }

                    let v_label = &graph.nodes[v_id];

                    if node_set.contains(v_label) {

                        sub.add_edge(
                            label,
                            v_label,
                            weight.clone(),
                        );
                    }
                }
            }
        }
    }

    sub
}

/// Computes the union of two graphs.
/// This assumes the graphs share the same vertex universe (labels match = same vertex).
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::core::Expr;
/// use rssn::symbolic::graph::Graph;
/// use rssn::symbolic::graph_operations::union;
///
/// let mut g1 = Graph::new(false);
///
/// g1.add_edge(
///     &"A",
///     &"B",
///     Expr::Constant(1.0),
/// );
///
/// let mut g2 = Graph::new(false);
///
/// g2.add_edge(
///     &"B",
///     &"C",
///     Expr::Constant(1.0),
/// );
///
/// let u = union(&g1, &g2);
///
/// assert_eq!(u.node_count(), 3);
/// ```
#[must_use]

pub fn union<V : Eq + Hash + Clone + Debug>(
    g1 : &Graph<V>,
    g2 : &Graph<V>,
) -> Graph<V> {

    let mut new_graph = g1.clone();

    for label in &g2.nodes {

        new_graph.add_node(label.clone());
    }

    for (u, v, weight) in g2.get_edges() {

        new_graph.add_edge(
            &g2.nodes[u],
            &g2.nodes[v],
            weight,
        );
    }

    new_graph
}

/// Computes the intersection of two graphs.
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::core::Expr;
/// use rssn::symbolic::graph::Graph;
/// use rssn::symbolic::graph_operations::intersection;
///
/// let mut g1 = Graph::new(false);
///
/// g1.add_edge(
///     &"A",
///     &"B",
///     Expr::Constant(1.0),
/// );
///
/// let mut g2 = Graph::new(false);
///
/// g2.add_edge(
///     &"A",
///     &"B",
///     Expr::Constant(1.0),
/// );
///
/// g2.add_edge(
///     &"B",
///     &"C",
///     Expr::Constant(1.0),
/// );
///
/// let i = intersection(&g1, &g2);
///
/// assert_eq!(
///     i.get_edges().len(),
///     1
/// );
/// ```
#[must_use]

pub fn intersection<V : Eq + Hash + Clone + Debug>(
    g1 : &Graph<V>,
    g2 : &Graph<V>,
) -> Graph<V> {

    let mut new_graph = Graph::new(g1.is_directed && g2.is_directed);

    let g1_nodes : HashSet<_> = g1
        .nodes
        .iter()
        .collect();

    let g2_nodes : HashSet<_> = g2
        .nodes
        .iter()
        .collect();

    for &node_label in g1_nodes.intersection(&g2_nodes) {

        new_graph.add_node((*node_label).clone());
    }

    // Iterate edges of g1
    for (u, v, weight) in g1.get_edges() {

        let u_label = &g1.nodes[u];

        let v_label = &g1.nodes[v];

        // Check if both nodes exist in g2
        if let (Some(u2), Some(v2)) = (
            g2.get_node_id(u_label),
            g2.get_node_id(v_label),
        ) {

            let has_edge = if g2.is_directed {

                g2.adj[u2]
                    .iter()
                    .any(|(n, w)| *n == v2 && *w == weight)
            } else {

                g2.adj[u2]
                    .iter()
                    .any(|(n, w)| *n == v2 && *w == weight)
            };

            if has_edge {

                new_graph.add_edge(
                    u_label,
                    v_label,
                    weight,
                );
            }
        }
    }

    new_graph
}

/// Computes the Cartesian product of two graphs.
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::core::Expr;
/// use rssn::symbolic::graph::Graph;
/// use rssn::symbolic::graph_operations::cartesian_product;
///
/// let mut g1 = Graph::new(false);
///
/// g1.add_edge(
///     &"A",
///     &"B",
///     Expr::Constant(1.0),
/// );
///
/// let mut g2 = Graph::new(false);
///
/// g2.add_edge(
///     &"1",
///     &"2",
///     Expr::Constant(1.0),
/// );
///
/// let prod = cartesian_product(&g1, &g2);
///
/// assert_eq!(prod.node_count(), 4);
/// ```
#[must_use]

pub fn cartesian_product<V : Eq + Hash + Clone + Debug + ToExpr>(
    g1 : &Graph<V>,
    g2 : &Graph<V>,
) -> Graph<Expr> {

    let mut new_graph = Graph::new(g1.is_directed || g2.is_directed);

    let mut node_map = HashMap::new();

    // Create nodes
    for (u_id, u_label) in g1
        .nodes
        .iter()
        .enumerate()
    {

        for (v_id, v_label) in g2
            .nodes
            .iter()
            .enumerate()
        {

            let new_label = Expr::Tuple(vec![
                u_label.to_expr(),
                v_label.to_expr(),
            ]);

            node_map.insert(
                (u_id, v_id),
                new_label.clone(),
            );

            new_graph.add_node(new_label);
        }
    }

    // Edges from G1: for each v in G2, copy G1
    for (v_id, _) in g2
        .nodes
        .iter()
        .enumerate()
    {

        for (u1, u2, weight) in g1.get_edges() {

            let n1 = &node_map[&(u1, v_id)];

            let n2 = &node_map[&(u2, v_id)];

            new_graph.add_edge(n1, n2, weight);
        }
    }

    // Edges from G2: for each u in G1, copy G2
    for (u_id, _) in g1
        .nodes
        .iter()
        .enumerate()
    {

        for (v1, v2, weight) in g2.get_edges() {

            let n1 = &node_map[&(u_id, v1)];

            let n2 = &node_map[&(u_id, v2)];

            new_graph.add_edge(n1, n2, weight);
        }
    }

    new_graph
}

/// Computes the Tensor product of two graphs.
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::core::Expr;
/// use rssn::symbolic::graph::Graph;
/// use rssn::symbolic::graph_operations::tensor_product;
///
/// let mut g1 = Graph::new(false);
///
/// g1.add_edge(
///     &"A",
///     &"B",
///     Expr::Constant(1.0),
/// );
///
/// let mut g2 = Graph::new(false);
///
/// g2.add_edge(
///     &"1",
///     &"2",
///     Expr::Constant(1.0),
/// );
///
/// let prod = tensor_product(&g1, &g2);
///
/// assert_eq!(prod.node_count(), 4);
/// ```
#[must_use]

pub fn tensor_product<V : Eq + Hash + Clone + Debug + ToExpr>(
    g1 : &Graph<V>,
    g2 : &Graph<V>,
) -> Graph<Expr> {

    let mut new_graph = Graph::new(g1.is_directed || g2.is_directed);

    let mut node_map = HashMap::new();

    for (u_id, u_label) in g1
        .nodes
        .iter()
        .enumerate()
    {

        for (v_id, v_label) in g2
            .nodes
            .iter()
            .enumerate()
        {

            let new_label = Expr::Tuple(vec![
                u_label.to_expr(),
                v_label.to_expr(),
            ]);

            node_map.insert(
                (u_id, v_id),
                new_label.clone(),
            );

            new_graph.add_node(new_label);
        }
    }

    // Edge (u1, v1) ~ (u2, v2) if u1~u2 in G1 AND v1~v2 in G2
    for (u1, u2, w1) in g1.get_edges() {

        for (v1, v2, w2) in g2.get_edges() {

            let n1 = &node_map[&(u1, v1)];

            let n2 = &node_map[&(u2, v2)];

            let weight = Expr::Mul(
                std::sync::Arc::new(w1.clone()),
                std::sync::Arc::new(w2.clone()),
            );

            new_graph.add_edge(
                n1,
                n2,
                weight.clone(),
            );

            if !g1.is_directed && !g2.is_directed {

                let n3 = &node_map[&(u1, v2)];

                let n4 = &node_map[&(u2, v1)];

                new_graph.add_edge(
                    n3,
                    n4,
                    weight.clone(),
                );
            } else if !g1.is_directed {

                let n3 = &node_map[&(u2, v1)];

                let n4 = &node_map[&(u1, v2)];

                new_graph.add_edge(
                    n3,
                    n4,
                    weight.clone(),
                );
            } else if !g2.is_directed {

                let n3 = &node_map[&(u1, v2)];

                let n4 = &node_map[&(u2, v1)];

                new_graph.add_edge(
                    n3,
                    n4,
                    weight.clone(),
                );
            }
        }
    }

    new_graph
}

/// Computes the complement of a graph.
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::core::Expr;
/// use rssn::symbolic::graph::Graph;
/// use rssn::symbolic::graph_operations::complement;
///
/// let mut g = Graph::new(false);
///
/// g.add_edge(
///     &"A",
///     &"B",
///     Expr::Constant(1.0),
/// );
///
/// g.add_edge(
///     &"B",
///     &"C",
///     Expr::Constant(1.0),
/// );
///
/// g.add_edge(
///     &"C",
///     &"A",
///     Expr::Constant(1.0),
/// );
///
/// let comp = complement(&g);
///
/// assert_eq!(
///     comp.get_edges()
///         .len(),
///     0
/// );
/// ```
#[must_use]

pub fn complement<V : Eq + Hash + Clone + Debug>(graph : &Graph<V>) -> Graph<V> {

    let mut new_graph = Graph::new(graph.is_directed);

    for label in &graph.nodes {

        new_graph.add_node(label.clone());
    }

    let n = graph.nodes.len();

    for i in 0 .. n {

        for j in 0 .. n {

            if i == j {

                continue;
            }

            if !graph.is_directed && i > j {

                continue;
            }

            let has_edge = graph.adj[i]
                .iter()
                .any(|(neighbor, _)| *neighbor == j);

            if !has_edge {

                new_graph.add_edge(
                    &graph.nodes[i],
                    &graph.nodes[j],
                    Expr::Constant(1.0),
                );
            }
        }
    }

    new_graph
}

/// Computes the disjoint union of two graphs.
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::core::Expr;
/// use rssn::symbolic::graph::Graph;
/// use rssn::symbolic::graph_operations::disjoint_union;
///
/// let mut g1 = Graph::new(false);
///
/// g1.add_node("A");
///
/// let mut g2 = Graph::new(false);
///
/// g2.add_node("A");
///
/// let du = disjoint_union(&g1, &g2);
///
/// assert_eq!(du.node_count(), 2);
/// ```
#[must_use]

pub fn disjoint_union<V : Eq + Hash + Clone + Debug + ToExpr>(
    g1 : &Graph<V>,
    g2 : &Graph<V>,
) -> Graph<Expr> {

    let mut new_graph = Graph::new(g1.is_directed && g2.is_directed);

    for label in &g1.nodes {

        let new_label = Expr::Tuple(vec![
            Expr::Constant(0.0),
            label.to_expr(),
        ]);

        new_graph.add_node(new_label);
    }

    for label in &g2.nodes {

        let new_label = Expr::Tuple(vec![
            Expr::Constant(1.0),
            label.to_expr(),
        ]);

        new_graph.add_node(new_label);
    }

    for (u, v, weight) in g1.get_edges() {

        let u_label = Expr::Tuple(vec![
            Expr::Constant(0.0),
            g1.nodes[u].to_expr(),
        ]);

        let v_label = Expr::Tuple(vec![
            Expr::Constant(0.0),
            g1.nodes[v].to_expr(),
        ]);

        new_graph.add_edge(
            &u_label,
            &v_label,
            weight,
        );
    }

    for (u, v, weight) in g2.get_edges() {

        let u_label = Expr::Tuple(vec![
            Expr::Constant(1.0),
            g2.nodes[u].to_expr(),
        ]);

        let v_label = Expr::Tuple(vec![
            Expr::Constant(1.0),
            g2.nodes[v].to_expr(),
        ]);

        new_graph.add_edge(
            &u_label,
            &v_label,
            weight,
        );
    }

    new_graph
}

/// Computes the join of two graphs (Zykov sum).
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::core::Expr;
/// use rssn::symbolic::graph::Graph;
/// use rssn::symbolic::graph_operations::join;
///
/// let mut g1 = Graph::new(false);
///
/// g1.add_node("A");
///
/// let mut g2 = Graph::new(false);
///
/// g2.add_node("B");
///
/// let j = join(&g1, &g2);
///
/// assert_eq!(
///     j.get_edges().len(),
///     1
/// );
/// ```
#[must_use]

pub fn join<V : Eq + Hash + Clone + Debug + ToExpr>(
    g1 : &Graph<V>,
    g2 : &Graph<V>,
) -> Graph<Expr> {

    let mut new_graph = disjoint_union(g1, g2);

    let g1_labels : Vec<Expr> = g1
        .nodes
        .iter()
        .map(|l| {

            Expr::Tuple(vec![
                Expr::Constant(0.0),
                l.to_expr(),
            ])
        })
        .collect();

    let g2_labels : Vec<Expr> = g2
        .nodes
        .iter()
        .map(|l| {

            Expr::Tuple(vec![
                Expr::Constant(1.0),
                l.to_expr(),
            ])
        })
        .collect();

    for u_label in &g1_labels {

        for v_label in &g2_labels {

            new_graph.add_edge(
                u_label,
                v_label,
                Expr::Constant(1.0),
            );

            if new_graph.is_directed {

                new_graph.add_edge(
                    v_label,
                    u_label,
                    Expr::Constant(1.0),
                );
            }
        }
    }

    new_graph
}
