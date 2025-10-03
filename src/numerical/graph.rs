//! # Numerical Graph Algorithms
//!
//! This module provides graph data structures and algorithms tailored for numerical
//! applications. It includes a weighted graph representation and an implementation
//! of Dijkstra's algorithm for finding shortest paths.

use std::cmp::Ordering;
use std::collections::BinaryHeap;

/// Represents a graph with weighted edges for numerical algorithms.
/// The graph is represented by an adjacency list.
pub struct Graph {
    adj: Vec<Vec<(usize, f64)>>, // (neighbor_index, edge_weight)
}

#[derive(Copy, Clone, PartialEq)]
pub struct State {
    cost: f64,
    position: usize,
}

// Manual implementation of Ord for State to make it a min-heap.
impl Ord for State {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .cost
            .partial_cmp(&self.cost)
            .unwrap_or(Ordering::Equal)
    }
}

impl PartialOrd for State {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Eq for State {}

impl Graph {
    /// Creates a new graph with a specified number of nodes.
    ///
    /// The graph is initialized with no edges.
    ///
    /// # Arguments
    /// * `num_nodes` - The total number of nodes in the graph.
    ///
    /// # Returns
    /// A new `Graph` instance.
    pub fn new(num_nodes: usize) -> Self {
        Graph {
            adj: vec![vec![]; num_nodes],
        }
    }

    /// Adds a directed edge with a weight between two nodes.
    ///
    /// # Arguments
    /// * `u` - The index of the source node.
    /// * `v` - The index of the destination node.
    /// * `weight` - The weight of the edge.
    pub fn add_edge(&mut self, u: usize, v: usize, weight: f64) {
        self.adj[u].push((v, weight));
    }

    /// Returns the total number of nodes in the graph.
    pub fn num_nodes(&self) -> usize {
        self.adj.len()
    }

    /// Returns an immutable slice of the neighbors and edge weights for a given node.
    ///
    /// # Arguments
    /// * `u` - The index of the node.
    ///
    /// # Returns
    /// A slice of `(usize, f64)` tuples, where each tuple is `(neighbor_index, edge_weight)`.
    pub fn adj(&self, u: usize) -> &[(usize, f64)] {
        &self.adj[u]
    }
}

/// Finds the shortest paths from a single source node to all other nodes
/// using Dijkstra's algorithm.
///
/// Dijkstra's algorithm is a greedy algorithm that solves the single-source
/// shortest path problem for a graph with non-negative edge weights.
///
/// # Arguments
/// * `graph` - The graph to search.
/// * `start_node` - The index of the starting node.
///
/// # Returns
/// A tuple containing:
///   - A vector of distances from the start node to each node.
///   - A vector of predecessors to reconstruct the shortest paths.
pub fn dijkstra(graph: &Graph, start_node: usize) -> (Vec<f64>, Vec<Option<usize>>) {
    let num_nodes = graph.adj.len();
    let mut dist: Vec<f64> = vec![f64::INFINITY; num_nodes];
    let mut prev: Vec<Option<usize>> = vec![None; num_nodes];
    let mut heap = BinaryHeap::new();

    dist[start_node] = 0.0;
    heap.push(State {
        cost: 0.0,
        position: start_node,
    });

    while let Some(State { cost, position }) = heap.pop() {
        if cost > dist[position] {
            continue;
        }

        for &(neighbor, weight) in &graph.adj[position] {
            if dist[position] + weight < dist[neighbor] {
                dist[neighbor] = dist[position] + weight;
                prev[neighbor] = Some(position);
                heap.push(State {
                    cost: dist[neighbor],
                    position: neighbor,
                });
            }
        }
    }
    (dist, prev)
}
