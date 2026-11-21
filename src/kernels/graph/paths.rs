use std::collections::VecDeque;

use super::Graph;
use super::GraphStructureError;
use super::topological_sort;

/// Finds unweighted distances from the nearest of several source vertices.
///
/// The result is indexed by vertex. Unreachable vertices have `usize::MAX`
/// distance, including every vertex when `sources` is empty. Repeated
/// sources are visited once and weights are ignored. Work is linear in
/// the graph size and the length of `sources`.
///
/// # Errors
///
/// Reports an invalid source or a reached edge with an invalid endpoint.
#[allow(
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    reason = "the queue contains only valid vertices and a reachable simple path uses fewer edges than vertices"
)]
pub fn multi_source_bfs(
    graph: &Graph,
    sources: &[usize],
) -> Result<Vec<usize>, GraphStructureError> {
    let mut distances = vec![usize::MAX; graph.num_nodes()];
    let mut queue = VecDeque::new();
    for &source in sources {
        let distance = distances
            .get_mut(source)
            .ok_or(GraphStructureError::InvalidVertex)?;
        if *distance == usize::MAX {
            *distance = 0;
            queue.push_back(source);
        }
    }

    while let Some(vertex) = queue.pop_front() {
        let next_distance = distances[vertex] + 1;
        for &(target, _) in graph.adj(vertex) {
            let distance = distances
                .get_mut(target)
                .ok_or(GraphStructureError::InvalidVertex)?;
            if *distance == usize::MAX {
                *distance = next_distance;
                queue.push_back(target);
            }
        }
    }
    Ok(distances)
}

/// Distances and predecessors from one source in vertex order.
#[derive(Debug, Clone, PartialEq)]
pub struct ShortestPaths {
    /// Shortest path lengths, with positive infinity for unreachable vertices.
    pub distances: Vec<f64>,
    /// Predecessor vertices. The source and unreachable vertices have none.
    pub previous: Vec<Option<usize>>,
}

impl ShortestPaths {
    fn new(
        n: usize,
        start: usize,
    ) -> Result<Self, GraphStructureError> {
        if start >= n {
            return Err(GraphStructureError::InvalidVertex);
        }
        let mut distances = vec![f64::INFINITY; n];
        if let Some(distance) = distances.get_mut(start) {
            *distance = 0.0;
        }
        Ok(Self {
            distances,
            previous: vec![None; n],
        })
    }

    #[allow(
        clippy::arithmetic_side_effects,
        clippy::indexing_slicing,
        reason = "edge endpoints were validated before relaxation and the sum is checked for overflow"
    )]
    fn relax(
        &mut self,
        source: usize,
        target: usize,
        weight: f64,
    ) -> Result<bool, GraphStructureError> {
        let candidate = self.distances[source] + weight;
        if !candidate.is_finite() {
            return Err(GraphStructureError::NumericOverflow);
        }
        if candidate < self.distances[target] {
            self.distances[target] = candidate;
            self.previous[target] = Some(source);
            return Ok(true);
        }
        Ok(false)
    }
}

fn validate_weights(graph: &Graph) -> Result<(), GraphStructureError> {
    let n = graph.num_nodes();
    for &(target, weight) in graph.adj.iter().flatten() {
        if target >= n {
            return Err(GraphStructureError::InvalidVertex);
        }
        if !weight.is_finite() {
            return Err(GraphStructureError::NonFiniteWeight);
        }
    }
    Ok(())
}

/// Solves single-source shortest paths with Bellman–Ford.
///
/// Negative weights are allowed. The relaxation sweep exits early once
/// distances stabilize. A negative cycle anywhere reachable from `start`
/// is an error. Unreachable cycles leave the result unaffected.
///
/// # Errors
///
/// Returns an error for an invalid vertex, non-finite weight, arithmetic
/// overflow, or reachable negative cycle.
#[allow(
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    reason = "vertices index equally sized graph and distance arrays and the source is in range"
)]
pub fn bellman_ford(
    graph: &Graph,
    start: usize,
) -> Result<ShortestPaths, GraphStructureError> {
    let n = graph.num_nodes();
    let mut paths = ShortestPaths::new(n, start)?;
    validate_weights(graph)?;
    for pass in 0..n {
        let mut changed = false;
        for (source, edges) in graph.adj.iter().enumerate() {
            if paths.distances[source].is_finite() {
                for &(target, weight) in edges {
                    changed |= paths.relax(source, target, weight)?;
                }
            }
        }
        if !changed {
            return Ok(paths);
        }
        if pass == n - 1 {
            return Err(GraphStructureError::NegativeCycle);
        }
    }
    Ok(paths)
}

/// Solves single-source shortest paths on a directed acyclic graph.
///
/// Relaxes edges once in topological order. Work and space are linear in
/// the graph size, and finite negative weights are supported.
///
/// # Errors
///
/// Returns an error for an invalid vertex, non-finite weight, arithmetic
/// overflow, or a directed cycle anywhere in the graph.
#[allow(
    clippy::indexing_slicing,
    reason = "the order and endpoints were validated before relaxation"
)]
pub fn dag_shortest_paths(
    graph: &Graph,
    start: usize,
) -> Result<ShortestPaths, GraphStructureError> {
    let mut paths = ShortestPaths::new(graph.num_nodes(), start)?;
    let order = topological_sort(graph)?;
    validate_weights(graph)?;
    for source in order {
        if paths.distances[source].is_finite() {
            for &(target, weight) in graph.adj(source) {
                paths.relax(source, target, weight)?;
            }
        }
    }
    Ok(paths)
}
