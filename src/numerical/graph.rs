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
    adj: Vec<Vec<(usize, f64)>>,
}

/// Represents a state in Dijkstra's algorithm.
#[derive(Copy, Clone, PartialEq)]

pub struct State {
    cost: f64,
    position: usize,
}

impl Ord for State {
    fn cmp(
        &self,
        other: &Self,
    ) -> Ordering {

        other
            .cost
            .partial_cmp(&self.cost)
            .unwrap_or(Ordering::Equal)
    }
}

impl PartialOrd for State {
    fn partial_cmp(
        &self,
        other: &Self,
    ) -> Option<Ordering> {

        Some(self.cmp(other))
    }
}

impl Eq for State {
}

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
    #[must_use]

    pub fn new(
        num_nodes: usize
    ) -> Self {

        Self {
            adj: vec![
                vec![];
                num_nodes
            ],
        }
    }

    /// Adds a directed edge with a weight between two nodes.
    ///
    /// # Arguments
    /// * `u` - The index of the source node.
    /// * `v` - The index of the destination node.
    /// * `weight` - The weight of the edge.

    pub fn add_edge(
        &mut self,
        u: usize,
        v: usize,
        weight: f64,
    ) {

        self.adj[u].push((v, weight));
    }

    /// Returns the total number of nodes in the graph.
    #[must_use]

    pub const fn num_nodes(
        &self
    ) -> usize {

        self.adj.len()
    }

    /// Returns an immutable slice of the neighbors and edge weights for a given node.
    ///
    /// # Arguments
    /// * `u` - The index of the node.
    ///
    /// # Returns
    /// A slice of `(usize, f64)` tuples, where each tuple is `(neighbor_index, edge_weight)`.
    #[must_use]

    pub fn adj(
        &self,
        u: usize,
    ) -> &[(usize, f64)] {

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
#[must_use]

pub fn dijkstra(
    graph: &Graph,
    start_node: usize,
) -> (
    Vec<f64>,
    Vec<Option<usize>>,
) {

    let num_nodes = graph.adj.len();

    let mut dist: Vec<f64> =
        vec![f64::INFINITY; num_nodes];

    let mut prev: Vec<Option<usize>> =
        vec![None; num_nodes];

    let mut heap = BinaryHeap::new();

    dist[start_node] = 0.0;

    heap.push(State {
        cost: 0.0,
        position: start_node,
    });

    while let Some(State {
        cost,
        position,
    }) = heap.pop()
    {

        if cost > dist[position] {

            continue;
        }

        for &(neighbor, weight) in
            &graph.adj[position]
        {

            if dist[position] + weight
                < dist[neighbor]
            {

                dist[neighbor] = dist
                    [position]
                    + weight;

                prev[neighbor] =
                    Some(position);

                heap.push(State {
                    cost: dist
                        [neighbor],
                    position: neighbor,
                });
            }
        }
    }

    (dist, prev)
}

/// Performs a Breadth-First Search (BFS) starting from a node.
/// Returns the unweighted shortest path distances from the start node to all other reachable nodes.
/// Returns `usize::MAX` for unreachable nodes.

#[must_use]

pub fn bfs(
    graph: &Graph,
    start_node: usize,
) -> Vec<usize> {

    let num_nodes = graph.num_nodes();

    let mut dist =
        vec![usize::MAX; num_nodes];

    let mut queue =
        std::collections::VecDeque::new(
        );

    dist[start_node] = 0;

    queue.push_back(start_node);

    while let Some(u) =
        queue.pop_front()
    {

        for &(v, _) in graph.adj(u) {

            if dist[v] == usize::MAX {

                dist[v] = dist[u] + 1;

                queue.push_back(v);
            }
        }
    }

    dist
}

/// Computes the `PageRank` of the nodes in the graph.
///
/// # Arguments
/// * `graph` - The graph.
/// * `damping_factor` - The probability of continuing following links (usually 0.85).
/// * `tolerance` - The convergence tolerance.
/// * `max_iter` - The maximum number of iterations.
///
/// # Returns
/// A vector of scores summing to 1.

#[must_use]

pub fn page_rank(
    graph: &Graph,
    damping_factor: f64,
    tolerance: f64,
    max_iter: usize,
) -> Vec<f64> {

    let num_nodes = graph.num_nodes();

    if num_nodes == 0 {

        return vec![];
    }

    let initial_score =
        1.0 / num_nodes as f64;

    let mut scores =
        vec![initial_score; num_nodes];

    let mut new_scores =
        vec![0.0; num_nodes];

    // Calculate out-degree for each node
    let  _out_degree =
        vec![0; num_nodes];
    let out_degree: Vec<usize> = (0..num_nodes).map(|u| graph.adj(u).len()).collect();

    for _ in 0 .. max_iter {

        let mut total_sink_score = 0.0;

        for u in 0 .. num_nodes {

            if out_degree[u] == 0 {

                total_sink_score +=
                    scores[u];
            }
        }

        let base_score = (1.0
            - damping_factor)
            / num_nodes as f64;

        let sink_share = damping_factor
            * total_sink_score
            / num_nodes as f64;

        new_scores.iter_mut().for_each(|score| *score = base_score + sink_share);

        for u in 0 .. num_nodes {

            if out_degree[u] > 0 {

                let share =
                    damping_factor
                        * scores[u]
                        / out_degree[u]
                            as f64;

                for &(v, _) in
                    graph.adj(u)
                {

                    new_scores[v] +=
                        share;
                }
            }
        }

        // Check convergence
        let mut diff = 0.0;

        for i in 0 .. num_nodes {

            diff += (new_scores[i]
                - scores[i])
                .abs();
        }

        scores.copy_from_slice(
            &new_scores,
        );

        if diff < tolerance {

            break;
        }
    }

    scores
}

/// Solves the All-Pairs Shortest Path problem using the Floyd-Warshall algorithm.
///
/// # Returns
/// A flattened vector of size `n * n` representing the distance matrix.
/// `dist[i * n + j]` is the distance from i to j.

#[must_use]

pub fn floyd_warshall(
    graph: &Graph
) -> Vec<f64> {

    let n = graph.num_nodes();

    let mut dist =
        vec![f64::INFINITY; n * n];

    // Initialize distances
    for i in 0 .. n {

        dist[i * n + i] = 0.0;

        for &(j, w) in graph.adj(i) {

            dist[i * n + j] =
                dist[i * n + j].min(w);
        }
    }

    for k in 0 .. n {

        for i in 0 .. n {

            for j in 0 .. n {

                let d_ik =
                    dist[i * n + k];

                let d_kj =
                    dist[k * n + j];

                if d_ik + d_kj
                    < dist[i * n + j]
                {

                    dist[i * n + j] =
                        d_ik + d_kj;
                }
            }
        }
    }

    dist
}

/// Finds the connected components of the graph.
/// Returns a vector where each element corresponds to a node and contains its component ID.

#[must_use]

pub fn connected_components(
    graph: &Graph
) -> Vec<usize> {

    let num_nodes = graph.num_nodes();

    let mut component =
        vec![usize::MAX; num_nodes];

    let mut current_component = 0;

    for i in 0 .. num_nodes {

        if component[i] == usize::MAX {

            let mut queue = std::collections::VecDeque::new();

            queue.push_back(i);

            component[i] =
                current_component;

            while let Some(u) =
                queue.pop_front()
            {

                for &(v, _) in
                    graph.adj(u)
                {

                    if component[v]
                        == usize::MAX
                    {

                        component[v] = current_component;

                        queue
                            .push_back(
                                v,
                            );
                    }
                }
            }

            current_component += 1;
        }
    }

    component
}

/// Computes the Minimum Spanning Tree (MST) using Prim's algorithm.
/// Returns a Graph representing the MST.
/// Assumes graph is connected (or computes MST forest).

#[must_use]

pub fn minimum_spanning_tree(
    graph: &Graph
) -> Graph {

    let num_nodes = graph.num_nodes();

    let mut mst = Graph::new(num_nodes);

    if num_nodes == 0 {

        return mst;
    }

    let mut visited =
        vec![false; num_nodes];

    let mut min_edge =
        vec![f64::INFINITY; num_nodes];

    let mut parent =
        vec![None; num_nodes];

    let mut heap = BinaryHeap::new();

    // Start from node 0 (or iterate if disconnected)
    for start_node in 0 .. num_nodes {

        if visited[start_node] {

            continue;
        }

        min_edge[start_node] = 0.0;

        heap.push(State {
            cost: 0.0,
            position: start_node,
        });

        while let Some(State {
            cost,
            position: u,
        }) = heap.pop()
        {

            if visited[u] {

                continue;
            }

            visited[u] = true;

            if let Some(p) = parent[u] {

                // For undirected logic on directed graph, we'd need to be careful.
                // Usually MST is defined for undirected graphs.
                // Here our Graph is directed.
                // If it represents undirected, edges are doubled.
                // We add edge p->u and u->p to MST.
                mst.add_edge(
                    p, u, cost,
                );

                mst.add_edge(
                    u, p, cost,
                );
            }

            for &(v, weight) in
                graph.adj(u)
            {

                if !visited[v]
                    && weight
                        < min_edge[v]
                {

                    min_edge[v] =
                        weight;

                    parent[v] = Some(u);

                    heap.push(State {
                        cost: weight,
                        position: v,
                    });
                }
            }
        }
    }

    mst
}
