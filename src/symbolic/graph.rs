//! # Symbolic Graph Data Structure
//!
//! This module provides a generic symbolic graph data structure that can represent
//! both directed and undirected graphs, as well as hypergraphs. Nodes can be labeled
//! with any type `V`, and edges can have symbolic weights (`Expr`). It includes
//! functionalities for adding nodes and edges, retrieving neighbors, calculating degrees,
//! and converting the graph to various matrix representations (adjacency, incidence, Laplacian).

use crate::symbolic::core::Expr;
use std::collections::HashMap;
use std::fmt::Debug;
use std::hash::Hash;

/// Represents a generic symbolic graph.
/// V is the type for vertex labels (e.g., String, Expr).
#[derive(Debug, Clone)]
pub struct Graph<V>
where
    V: Eq + Hash + Clone + Debug,
{
    pub(crate) nodes: Vec<V>,
    pub(crate) node_map: HashMap<V, usize>,
    pub(crate) adj: Vec<Vec<(usize, Expr)>>, // (neighbor_index, edge_weight)
    pub(crate) rev_adj: Vec<Vec<(usize, Expr)>>, // For directed graphs
    pub(crate) hyperedges: Vec<(std::collections::HashSet<usize>, Expr)>,
    pub(crate) is_directed: bool,
}

impl<V> Graph<V>
where
    V: Eq + Hash + Clone + Debug,
{
    /// Creates a new graph.
    ///
    /// # Arguments
    /// * `is_directed` - A boolean indicating whether the graph is directed (`true`) or undirected (`false`).
    ///
    /// # Returns
    /// A new `Graph` instance.
    pub fn new(is_directed: bool) -> Self {
        Graph {
            nodes: Vec::new(),
            node_map: HashMap::new(),
            adj: Vec::new(),
            rev_adj: Vec::new(),
            hyperedges: Vec::new(),
            is_directed,
        }
    }

    /// Adds a node with a given label to the graph.
    ///
    /// If a node with the same label already exists, its existing ID is returned.
    ///
    /// # Arguments
    /// * `label` - The label of the node (type `V`).
    ///
    /// # Returns
    /// The internal `usize` ID of the node.
    pub fn add_node(&mut self, label: V) -> usize {
        if let Some(&id) = self.node_map.get(&label) {
            return id;
        }
        let id = self.nodes.len();
        self.nodes.push(label.clone());
        self.node_map.insert(label, id);
        self.adj.push(Vec::new());
        self.rev_adj.push(Vec::new());
        id
    }

    /// Adds an edge between two nodes.
    ///
    /// If the graph is undirected, an edge is added in both directions.
    /// If the nodes do not exist, they are added to the graph.
    ///
    /// # Arguments
    /// * `from_label` - The label of the source node.
    /// * `to_label` - The label of the destination node.
    /// * `weight` - The symbolic weight of the edge.
    pub fn add_edge(&mut self, from_label: &V, to_label: &V, weight: Expr) {
        let from_id = self.add_node(from_label.clone());
        let to_id = self.add_node(to_label.clone());

        self.adj[from_id].push((to_id, weight.clone()));
        self.rev_adj[to_id].push((from_id, weight.clone()));

        if !self.is_directed {
            self.adj[to_id].push((from_id, weight.clone()));
            self.rev_adj[from_id].push((to_id, weight));
        }
    }

    /// Gets the internal ID of a node given its label.
    ///
    /// # Arguments
    /// * `label` - The label of the node.
    ///
    /// # Returns
    /// An `Option<usize>` containing the node's ID if found, `None` otherwise.
    pub fn get_node_id(&self, label: &V) -> Option<usize> {
        self.node_map.get(label).cloned()
    }

    /// Gets the neighbors of a node.
    ///
    /// # Arguments
    /// * `node_id` - The internal ID of the node.
    ///
    /// # Returns
    /// An iterator over `(neighbor_id, edge_weight)` tuples.
    pub fn neighbors(&self, node_id: usize) -> impl Iterator<Item = &(usize, Expr)> {
        self.adj.get(node_id).into_iter().flatten()
    }

    /// Gets the out-degree of a node.
    ///
    /// The out-degree is the number of edges originating from the node.
    /// For undirected graphs, this is equivalent to the degree.
    ///
    /// # Arguments
    /// * `node_id` - The internal ID of the node.
    ///
    /// # Returns
    /// The out-degree as a `usize`.
    pub fn out_degree(&self, node_id: usize) -> usize {
        self.adj.get(node_id).map_or(0, |v| v.len())
    }

    /// Gets the in-degree of a node.
    ///
    /// The in-degree is the number of edges terminating at the node.
    /// For undirected graphs, this is equivalent to the degree.
    ///
    /// # Arguments
    /// * `node_id` - The internal ID of the node.
    ///
    /// # Returns
    /// The in-degree as a `usize`.
    pub fn in_degree(&self, node_id: usize) -> usize {
        self.rev_adj.get(node_id).map_or(0, |v| v.len())
    }

    /// Returns a list of all edges in the graph.
    ///
    /// For undirected graphs, each edge is listed only once (e.g., `(u, v)` but not `(v, u)`).
    ///
    /// # Returns
    /// A `Vec<(usize, usize, Expr)>` where each tuple is `(from_node_id, to_node_id, edge_weight)`.
    pub fn get_edges(&self) -> Vec<(usize, usize, Expr)> {
        let mut edges = Vec::new();
        for (u, neighbors) in self.adj.iter().enumerate() {
            for &(v, ref weight) in neighbors {
                if !self.is_directed && u > v {
                    continue;
                } // Avoid duplicates in undirected graphs
                edges.push((u, v, weight.clone()));
            }
        }
        edges
    }

    /// Adds a hyperedge that connects a set of vertices.
    ///
    /// A hyperedge is an edge that can connect any number of vertices.
    ///
    /// # Arguments
    /// * `labels` - A slice of node labels (`V`) that the hyperedge connects.
    /// * `weight` - The symbolic weight of the hyperedge.
    pub fn add_hyperedge(&mut self, labels: &[V], weight: Expr) {
        let ids: std::collections::HashSet<usize> = labels
            .iter()
            .map(|label| self.add_node(label.clone()))
            .collect();
        self.hyperedges.push((ids, weight));
    }

    /// Returns the adjacency matrix of the graph.
    ///
    /// The adjacency matrix `A` is a square matrix where `A[i][j]` represents
    /// the weight of the edge from node `i` to node `j`. If no edge exists, the value is 0.
    ///
    /// # Returns
    /// An `Expr::Matrix` representing the adjacency matrix.
    pub fn to_adjacency_matrix(&self) -> Expr {
        let n = self.nodes.len();
        let mut matrix = vec![vec![Expr::Constant(0.0); n]; n];
        for u in 0..n {
            if let Some(neighbors) = self.adj.get(u) {
                for &(v, ref weight) in neighbors {
                    matrix[u][v] = weight.clone();
                }
            }
        }
        Expr::Matrix(matrix)
    }

    /// Returns the incidence matrix of the graph.
    ///
    /// The incidence matrix `B` is a matrix where rows correspond to nodes and columns
    /// correspond to edges. `B[i][j]` is 1 if node `i` is incident to edge `j`.
    /// For directed graphs, it can be -1 for the source and 1 for the destination.
    ///
    /// # Returns
    /// An `Expr::Matrix` representing the incidence matrix.
    pub fn to_incidence_matrix(&self) -> Expr {
        let n = self.nodes.len();
        let edges = self.get_edges();
        let m = edges.len();
        let mut matrix = vec![vec![Expr::Constant(0.0); m]; n];

        for (j, &(u, v, _)) in edges.iter().enumerate() {
            if self.is_directed {
                matrix[u][j] = Expr::Constant(-1.0);
                matrix[v][j] = Expr::Constant(1.0);
            } else {
                matrix[u][j] = Expr::Constant(1.0);
                matrix[v][j] = Expr::Constant(1.0);
            }
        }
        Expr::Matrix(matrix)
    }

    /// Returns the Laplacian matrix of the graph (`L = D - A`).
    ///
    /// The Laplacian matrix `L` is defined as the difference between the degree matrix `D`
    /// (a diagonal matrix with node degrees on the diagonal) and the adjacency matrix `A`.
    /// It is a fundamental matrix in spectral graph theory.
    ///
    /// # Returns
    /// An `Expr::Matrix` representing the Laplacian matrix.
    pub fn to_laplacian_matrix(&self) -> Expr {
        let n = self.nodes.len();
        let adj_matrix_expr = self.to_adjacency_matrix();
        let _adj_matrix = if let Expr::Matrix(m) = &adj_matrix_expr {
            m
        } else {
            return Expr::Variable("Error".to_string());
        };

        let mut deg_matrix = vec![vec![Expr::Constant(0.0); n]; n];
        for i in 0..n {
            let degree = self.out_degree(i); // For undirected, in_degree == out_degree
            deg_matrix[i][i] = Expr::Constant(degree as f64);
        }

        crate::symbolic::matrix::sub_matrices(&Expr::Matrix(deg_matrix), &adj_matrix_expr)
    }
}

/*
}
or(0, |v| v.len())
    }
}
*/
