use std::collections::VecDeque;

use super::Graph;
use super::GraphStructureError;

/// A maximum flow value and the partition of a corresponding minimum cut.
#[derive(Debug, Clone, PartialEq)]
pub struct MaxFlow {
    /// Total flow delivered to the sink in the units of the edge capacities.
    pub value: f64,
    /// One boolean per vertex, true for vertices on the source side.
    pub source_side: Vec<bool>,
}

#[derive(Clone, Copy)]
struct Arc {
    target: usize,
    reverse: usize,
    remaining: f64,
}

struct Residual {
    adjacency: Vec<Vec<Arc>>,
    level: Vec<usize>,
    next: Vec<usize>,
    queue: VecDeque<usize>,
    path: Vec<(usize, usize)>,
}

#[allow(
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    reason = "residual edges have paired reverse indices and level counters are bounded by the vertex count"
)]
impl Residual {
    fn new(graph: &Graph) -> Result<Self, GraphStructureError> {
        let n = graph.num_nodes();
        let mut adjacency = vec![Vec::new(); n];
        for (source, edges) in graph.adj.iter().enumerate() {
            for &(target, capacity) in edges {
                if target >= n {
                    return Err(GraphStructureError::InvalidVertex);
                }
                if !capacity.is_finite() {
                    return Err(GraphStructureError::NonFiniteWeight);
                }
                if capacity < 0.0 {
                    return Err(GraphStructureError::NegativeCapacity);
                }
                if capacity == 0.0 || source == target {
                    continue;
                }
                let forward = adjacency[source].len();
                let reverse = adjacency[target].len();
                adjacency[source].push(Arc {
                    target,
                    reverse,
                    remaining: capacity,
                });
                adjacency[target].push(Arc {
                    target: source,
                    reverse: forward,
                    remaining: 0.0,
                });
            }
        }
        Ok(Self {
            adjacency,
            level: vec![usize::MAX; n],
            next: vec![0; n],
            queue: VecDeque::new(),
            path: Vec::new(),
        })
    }

    fn build_level_graph(
        &mut self,
        source: usize,
    ) {
        self.level.fill(usize::MAX);
        self.level[source] = 0;
        self.queue.push_back(source);
        while let Some(vertex) = self.queue.pop_front() {
            let next_level = self.level[vertex] + 1;
            for edge in &self.adjacency[vertex] {
                if edge.remaining > 0.0 && self.level[edge.target] == usize::MAX {
                    self.level[edge.target] = next_level;
                    self.queue.push_back(edge.target);
                }
            }
        }
    }

    fn blocking_flow(
        &mut self,
        source: usize,
        sink: usize,
        value: &mut f64,
    ) -> Result<(), GraphStructureError> {
        self.next.fill(0);
        let mut vertex = source;
        loop {
            if vertex == sink {
                let amount = self
                    .path
                    .iter()
                    .map(|&(u, i)| self.adjacency[u][i].remaining)
                    .fold(f64::INFINITY, f64::min);
                *value += amount;
                if !value.is_finite() {
                    return Err(GraphStructureError::NumericOverflow);
                }
                for (u, i) in self.path.drain(..) {
                    let edge = &mut self.adjacency[u][i];
                    edge.remaining -= amount;
                    let (target, reverse) = (edge.target, edge.reverse);
                    let remaining = &mut self.adjacency[target][reverse].remaining;
                    *remaining += amount;
                    if !remaining.is_finite() {
                        return Err(GraphStructureError::NumericOverflow);
                    }
                }
                vertex = source;
                continue;
            }

            let next = &mut self.next[vertex];
            let next_level = self.level[vertex] + 1;
            while let Some(edge) = self.adjacency[vertex].get(*next) {
                if edge.remaining > 0.0 && self.level[edge.target] == next_level {
                    break;
                }
                *next += 1;
            }
            if let Some(edge) = self.adjacency[vertex].get(*next) {
                self.path.push((vertex, *next));
                vertex = edge.target;
            } else if let Some((parent, edge_index)) = self.path.pop() {
                self.next[parent] = edge_index + 1;
                vertex = parent;
            } else {
                return Ok(());
            }
        }
    }
}

/// Computes a maximum flow and a minimum source–sink cut using Dinic's
/// algorithm with an iterative blocking-flow traversal.
///
/// Arc weights are capacities. Parallel and antiparallel arcs are independent.
/// Loops and zero-capacity arcs make no contribution. The algorithm uses
/// `O(V + E)` auxiliary space and has an `O(V² E)` worst-case time bound.
///
/// # Errors
///
/// Returns an error for invalid endpoints, equal source and sink, non-finite
/// or negative capacities, and a flow value outside finite floating-point range.
#[allow(
    clippy::indexing_slicing,
    reason = "the sink index was checked against the length of the level array"
)]
pub fn maximum_flow(
    graph: &Graph,
    source: usize,
    sink: usize,
) -> Result<MaxFlow, GraphStructureError> {
    if source >= graph.num_nodes() || sink >= graph.num_nodes() {
        return Err(GraphStructureError::InvalidVertex);
    }
    if source == sink {
        return Err(GraphStructureError::IdenticalTerminals);
    }
    let mut residual = Residual::new(graph)?;
    let mut value = 0.0;
    loop {
        residual.build_level_graph(source);
        if residual.level[sink] == usize::MAX {
            return Ok(MaxFlow {
                value,
                source_side: residual
                    .level
                    .into_iter()
                    .map(|level| level != usize::MAX)
                    .collect(),
            });
        }
        residual.blocking_flow(source, sink, &mut value)?;
    }
}
