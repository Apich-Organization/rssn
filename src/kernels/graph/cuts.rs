use super::GraphStructureError;

/// Bridges and articulation vertices of an undirected multigraph.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UndirectedCuts {
    /// Bridge endpoints in depth-first completion order, with the parent first.
    pub bridges: Vec<[usize; 2]>,
    /// Vertices whose removal increases the component count, in vertex order.
    pub articulation_points: Vec<usize>,
}

/// Finds bridges and articulation points in an undirected multigraph.
///
/// Each endpoint pair describes one undirected edge. Parallel edges are
/// independent and prevent one another from being bridges. Self-loops have
/// no effect. A depth-first traversal considers roots in vertex order and
/// incident edges in input order. Runtime and storage are both linear in
/// the number of vertices and edges.
///
/// # Errors
///
/// Reports any edge endpoint outside `0..num_vertices`.
pub fn undirected_cuts(
    num_vertices: usize,
    edges: &[[usize; 2]],
) -> Result<UndirectedCuts, GraphStructureError> {
    let mut adjacency = vec![Vec::new(); num_vertices];
    for (edge_id, &[u, v]) in edges.iter().enumerate() {
        adjacency
            .get_mut(u)
            .ok_or(GraphStructureError::InvalidVertex)?
            .push((v, edge_id));
        adjacency
            .get_mut(v)
            .ok_or(GraphStructureError::InvalidVertex)?
            .push((u, edge_id));
    }
    Ok(undirected_cuts_by(&adjacency, |&edge| edge))
}

struct Frame<'a, Edge> {
    vertex: usize,
    remaining: std::slice::Iter<'a, Edge>,
    parent_edge: Option<usize>,
    children: usize,
}

/// Low-link traversal shared by checked endpoint and symbolic adjacency inputs.
#[allow(
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    reason = "both callers validate endpoints, and tree-edge counters are bounded by the vertex count"
)]
pub fn undirected_cuts_by<Edge>(
    adjacency: &[Vec<Edge>],
    target_and_id: impl Fn(&Edge) -> (usize, usize),
) -> UndirectedCuts {
    let unvisited = usize::MAX;
    let mut discovered = vec![unvisited; adjacency.len()];
    let mut low = vec![0; adjacency.len()];
    let mut articulation = vec![false; adjacency.len()];
    let mut time = 0;
    let mut bridges = Vec::new();
    let mut search = Vec::new();

    for (root, edges) in adjacency.iter().enumerate() {
        if discovered[root] != unvisited {
            continue;
        }
        discovered[root] = time;
        low[root] = time;
        time += 1;
        search.push(Frame {
            vertex: root,
            remaining: edges.iter(),
            parent_edge: None,
            children: 0,
        });
        while let Some(frame) = search.last_mut() {
            let u = frame.vertex;
            if let Some(edge) = frame.remaining.next() {
                let (v, id) = target_and_id(edge);
                if frame.parent_edge == Some(id) {
                    continue;
                }
                if discovered[v] == unvisited {
                    frame.children += 1;
                    discovered[v] = time;
                    low[v] = time;
                    time += 1;
                    search.push(Frame {
                        vertex: v,
                        remaining: adjacency[v].iter(),
                        parent_edge: Some(id),
                        children: 0,
                    });
                } else {
                    low[u] = low[u].min(discovered[v]);
                }
                continue;
            }

            if frame.parent_edge.is_none() {
                articulation[u] = frame.children > 1;
            }
            search.pop();
            if let Some(parent) = search.last() {
                let p = parent.vertex;
                low[p] = low[p].min(low[u]);
                if low[u] > discovered[p] {
                    bridges.push([p, u]);
                }
                if parent.parent_edge.is_some() && low[u] >= discovered[p] {
                    articulation[p] = true;
                }
            }
        }
    }

    UndirectedCuts {
        bridges,
        articulation_points: articulation
            .into_iter()
            .enumerate()
            .filter_map(|(vertex, is_cut)| is_cut.then_some(vertex))
            .collect(),
    }
}
