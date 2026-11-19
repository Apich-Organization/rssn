use std::collections::VecDeque;
use std::fmt;

use super::Graph;

/// Structural failure in a directed graph operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GraphStructureError {
    /// An edge points outside the vertex range.
    InvalidVertex,
    /// A directed cycle prevents a topological order.
    DirectedCycle,
}

impl fmt::Display for GraphStructureError {
    fn fmt(
        &self,
        formatter: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        formatter.write_str(match self {
            | Self::InvalidVertex => "edge points outside the graph",
            | Self::DirectedCycle => "graph contains a directed cycle",
        })
    }
}

impl std::error::Error for GraphStructureError {}

/// Returns a topological ordering of the vertices using Kahn's algorithm.
///
/// An edge from `u` to `v` puts `u` before `v`. Parallel edges are counted
/// separately, all vertices are returned, and edge weights are ignored.
/// The traversal takes linear time in the number of vertices and edges.
///
/// # Errors
///
/// Reports an invalid edge endpoint or a directed cycle, including a self-loop.
pub fn topological_sort(graph: &Graph) -> Result<Vec<usize>, GraphStructureError> {
    topological_by(&graph.adj, |edge| edge.0)
}

/// Kahn traversal shared by numerical and symbolic graph views.
#[allow(
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    reason = "endpoints are validated before traversal and degree counters are bounded by the adjacency lists"
)]
pub fn topological_by<Edge>(
    adjacency: &[Vec<Edge>],
    target: impl Fn(&Edge) -> usize,
) -> Result<Vec<usize>, GraphStructureError> {
    let n = adjacency.len();
    let mut indegrees = vec![0_usize; n];
    for edge in adjacency.iter().flatten() {
        *indegrees
            .get_mut(target(edge))
            .ok_or(GraphStructureError::InvalidVertex)? += 1;
    }
    let mut ready: VecDeque<_> = indegrees
        .iter()
        .enumerate()
        .filter_map(|(vertex, &degree)| (degree == 0).then_some(vertex))
        .collect();
    let mut order = Vec::with_capacity(n);
    while let Some(vertex) = ready.pop_front() {
        order.push(vertex);
        for edge in &adjacency[vertex] {
            let target = target(edge);
            indegrees[target] -= 1;
            if indegrees[target] == 0 {
                ready.push_back(target);
            }
        }
    }

    if order.len() == n {
        Ok(order)
    } else {
        Err(GraphStructureError::DirectedCycle)
    }
}

/// Labels vertices by strongly connected component in linear time.
///
/// Vertices with the same label can reach each other along directed edges.
/// Labels are dense and topological: every edge crossing components points
/// from a smaller label to a larger one. Edge weights are ignored. The
/// traversal uses explicit stacks so deep paths do not grow the call stack.
///
/// # Errors
///
/// Reports an invalid edge endpoint.
pub fn strongly_connected_components(graph: &Graph) -> Result<Vec<usize>, GraphStructureError> {
    strongly_connected_by(&graph.adj, |edge| edge.0).ok_or(GraphStructureError::InvalidVertex)
}

struct Traversal {
    sequence: usize,
    discovered: Vec<Option<usize>>,
    low: Vec<usize>,
    active: Vec<usize>,
    labels: Vec<usize>,
    components: usize,
}

impl Traversal {
    #[allow(
        clippy::arithmetic_side_effects,
        clippy::indexing_slicing,
        reason = "the vertex is in range and discovery counts at most one visit per vertex"
    )]
    fn enter(
        &mut self,
        vertex: usize,
    ) {
        self.discovered[vertex] = Some(self.sequence);
        self.low[vertex] = self.sequence;
        self.sequence += 1;
        self.active.push(vertex);
    }
}

/// Iterative Tarjan traversal shared by numerical and symbolic graph views.
#[allow(
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    reason = "vertices come from checked edges or the range of equally sized arrays, and each receives one label"
)]
pub fn strongly_connected_by<Edge>(
    adjacency: &[Vec<Edge>],
    target: impl Fn(&Edge) -> usize,
) -> Option<Vec<usize>> {
    let n = adjacency.len();
    let pending = usize::MAX;
    let mut state = Traversal {
        sequence: 0,
        discovered: vec![None; n],
        low: vec![0; n],
        active: Vec::new(),
        labels: vec![pending; n],
        components: 0,
    };
    let mut search = Vec::new();

    for (root, edges) in adjacency.iter().enumerate() {
        if state.discovered[root].is_some() {
            continue;
        }
        state.enter(root);
        search.push((root, edges.iter()));
        while let Some((vertex, remaining)) = search.last_mut() {
            let u = *vertex;
            if let Some(edge) = remaining.next() {
                let v = target(edge);
                match *state.discovered.get(v)? {
                    | None => {
                        state.enter(v);
                        search.push((v, adjacency[v].iter()));
                    },
                    | Some(index) if state.labels[v] == pending => {
                        state.low[u] = state.low[u].min(index);
                    },
                    | Some(_) => {},
                }
                continue;
            }

            if state.discovered[u] == Some(state.low[u]) {
                while let Some(vertex) = state.active.pop() {
                    state.labels[vertex] = state.components;
                    if vertex == u {
                        break;
                    }
                }
                state.components += 1;
            }
            search.pop();
            if let Some((parent, _)) = search.last() {
                state.low[*parent] = state.low[*parent].min(state.low[u]);
            }
        }
    }

    for label in &mut state.labels {
        *label = state.components - 1 - *label;
    }
    Some(state.labels)
}
