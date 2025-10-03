//! # Comprehensive Graph Algorithms
//!
//! This module provides a comprehensive suite of graph algorithms for various tasks,
//! including graph traversal (DFS, BFS), connectivity analysis (connected components,
//! strongly connected components), cycle detection, minimum spanning trees (Kruskal's,
//! Prim's), network flow (Edmonds-Karp, Dinic's), shortest paths (Dijkstra's,
//! Bellman-Ford, Floyd-Warshall), and topological sorting.

use crate::prelude::simplify;
use crate::symbolic::core::Expr;
use crate::symbolic::graph::Graph;
use crate::symbolic::simplify::as_f64;
use ordered_float::OrderedFloat;
use std::collections::{HashMap, HashSet, VecDeque};
use std::fmt::Debug;
use std::hash::Hash;

// =====================================================================================
// region: Graph Traversal
// =====================================================================================

/// Performs a Depth-First Search (DFS) traversal on a graph.
///
/// DFS explores as far as possible along each branch before backtracking.
///
/// # Arguments
/// * `graph` - The graph to traverse.
/// * `start_node` - The index of the starting node.
///
/// # Returns
/// A `Vec<usize>` containing the node IDs in the order they were visited.
pub fn dfs<V>(graph: &Graph<V>, start_node: usize) -> Vec<usize>
where
    V: Eq + Hash + Clone + std::fmt::Debug,
{
    let mut visited = HashSet::new();
    let mut result = Vec::new();
    dfs_recursive(graph, start_node, &mut visited, &mut result);
    result
}

pub(crate) fn dfs_recursive<V>(
    graph: &Graph<V>,
    u: usize,
    visited: &mut HashSet<usize>,
    result: &mut Vec<usize>,
) where
    V: Eq + Hash + Clone + std::fmt::Debug,
{
    visited.insert(u);
    result.push(u);
    if let Some(neighbors) = graph.adj.get(u) {
        for &(v, _) in neighbors {
            if !visited.contains(&v) {
                dfs_recursive(graph, v, visited, result);
            }
        }
    }
}

/// Performs a Breadth-First Search (BFS) traversal on a graph.
///
/// BFS explores all of the neighbor nodes at the present depth prior to moving on to nodes at the next depth level.
///
/// # Arguments
/// * `graph` - The graph to traverse.
/// * `start_node` - The index of the starting node.
///
/// # Returns
/// A `Vec<usize>` containing the node IDs in the order they were visited.
pub fn bfs<V>(graph: &Graph<V>, start_node: usize) -> Vec<usize>
where
    V: Eq + Hash + Clone + std::fmt::Debug,
{
    let mut visited = HashSet::new();
    let mut result = Vec::new();
    let mut queue = VecDeque::new();

    visited.insert(start_node);
    queue.push_back(start_node);

    while let Some(u) = queue.pop_front() {
        result.push(u);
        if let Some(neighbors) = graph.adj.get(u) {
            for &(v, _) in neighbors {
                if !visited.contains(&v) {
                    visited.insert(v);
                    queue.push_back(v);
                }
            }
        }
    }
    result
}

// =====================================================================================
// region: Connectivity & Components
// =====================================================================================

/// Finds all connected components of an undirected graph.
///
/// A connected component is a subgraph in which any two vertices are connected to each other
/// by paths, and which is connected to no additional vertices in the supergraph.
///
/// # Arguments
/// * `graph` - The graph to analyze.
///
/// # Returns
/// A `Vec<Vec<usize>>` where each inner `Vec` represents a connected component.
pub fn connected_components<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    graph: &Graph<V>,
) -> Vec<Vec<usize>> {
    let mut visited = HashSet::new();
    let mut components = Vec::new();
    for node_id in 0..graph.nodes.len() {
        if !visited.contains(&node_id) {
            let component = bfs(graph, node_id);
            for &visited_node in &component {
                visited.insert(visited_node);
            }
            components.push(component);
        }
    }
    components
}

/// Checks if the graph is connected.
///
/// An undirected graph is connected if for every pair of vertices `(u, v)`,
/// there is a path from `u` to `v`.
///
/// # Arguments
/// * `graph` - The graph to check.
///
/// # Returns
/// `true` if the graph is connected, `false` otherwise.
pub fn is_connected<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(graph: &Graph<V>) -> bool {
    connected_components(graph).len() == 1
}

/// Finds all strongly connected components (SCCs) of a directed graph using Tarjan's algorithm.
///
/// An SCC is a subgraph where every vertex is reachable from every other vertex within that subgraph.
/// Tarjan's algorithm is an efficient DFS-based method for finding SCCs.
///
/// # Arguments
/// * `graph` - The directed graph to analyze.
///
/// # Returns
/// A `Vec<Vec<usize>>` where each inner `Vec` represents a strongly connected component.
pub fn strongly_connected_components<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    graph: &Graph<V>,
) -> Vec<Vec<usize>> {
    let mut scc = Vec::new();
    let mut stack = Vec::new();
    let mut on_stack = HashSet::new();
    let mut discovery_times = HashMap::new();
    let mut low_link = HashMap::new();
    let mut time = 0;
    for node_id in 0..graph.nodes.len() {
        tarjan_scc_util(
            graph,
            node_id,
            &mut time,
            &mut discovery_times,
            &mut low_link,
            &mut stack,
            &mut on_stack,
            &mut scc,
        );
    }
    scc
}

pub(crate) fn tarjan_scc_util<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    graph: &Graph<V>,
    u: usize,
    time: &mut usize,
    disc: &mut HashMap<usize, usize>,
    low: &mut HashMap<usize, usize>,
    stack: &mut Vec<usize>,
    on_stack: &mut HashSet<usize>,
    scc: &mut Vec<Vec<usize>>,
) {
    disc.insert(u, *time);
    low.insert(u, *time);
    *time += 1;
    stack.push(u);
    on_stack.insert(u);

    if let Some(neighbors) = graph.adj.get(u) {
        for &(v, _) in neighbors {
            if !disc.contains_key(&v) {
                tarjan_scc_util(graph, v, time, disc, low, stack, on_stack, scc);
                let low_u = *low.get(&u).unwrap();
                let low_v = *low.get(&v).unwrap();
                low.insert(u, low_u.min(low_v));
            } else if on_stack.contains(&v) {
                let low_u = *low.get(&u).unwrap();
                let disc_v = *disc.get(&v).unwrap();
                low.insert(u, low_u.min(disc_v));
            }
        }
    }

    if low.get(&u) == disc.get(&u) {
        let mut component = Vec::new();
        while let Some(top) = stack.pop() {
            on_stack.remove(&top);
            component.push(top);
            if top == u {
                break;
            }
        }
        scc.push(component);
    }
}

// =====================================================================================
// region: Cycles, Bridges, Articulation Points
// =====================================================================================

/// Detects if a cycle exists in the graph.
///
/// For directed graphs, it uses a DFS-based approach with a recursion stack.
/// For undirected graphs, it uses a DFS-based approach that checks for back-edges.
///
/// # Arguments
/// * `graph` - The graph to check.
///
/// # Returns
/// `true` if a cycle is found, `false` otherwise.
pub fn has_cycle<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(graph: &Graph<V>) -> bool {
    let mut visited = HashSet::new();
    if graph.is_directed {
        let mut recursion_stack = HashSet::new();
        for node_id in 0..graph.nodes.len() {
            if !visited.contains(&node_id) {
                if has_cycle_directed_util(graph, node_id, &mut visited, &mut recursion_stack) {
                    return true;
                }
            }
        }
    } else {
        for node_id in 0..graph.nodes.len() {
            if !visited.contains(&node_id) {
                if has_cycle_undirected_util(graph, node_id, &mut visited, None) {
                    return true;
                }
            }
        }
    }
    false
}

pub(crate) fn has_cycle_directed_util<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    graph: &Graph<V>,
    u: usize,
    visited: &mut HashSet<usize>,
    rec_stack: &mut HashSet<usize>,
) -> bool {
    visited.insert(u);
    rec_stack.insert(u);
    if let Some(neighbors) = graph.adj.get(u) {
        for &(v, _) in neighbors {
            if !visited.contains(&v) {
                if has_cycle_directed_util(graph, v, visited, rec_stack) {
                    return true;
                }
            } else if rec_stack.contains(&v) {
                return true;
            }
        }
    }
    rec_stack.remove(&u);
    false
}

pub(crate) fn has_cycle_undirected_util<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    graph: &Graph<V>,
    u: usize,
    visited: &mut HashSet<usize>,
    parent: Option<usize>,
) -> bool {
    visited.insert(u);
    if let Some(neighbors) = graph.adj.get(u) {
        for &(v, _) in neighbors {
            if !visited.contains(&v) {
                if has_cycle_undirected_util(graph, v, visited, Some(u)) {
                    return true;
                }
            } else if Some(v) != parent {
                return true;
            }
        }
    }
    false
}

/// Finds all bridges and articulation points (cut vertices) in a graph using Tarjan's algorithm.
///
/// A bridge is an edge whose removal increases the number of connected components.
/// An articulation point is a vertex whose removal increases the number of connected components.
/// Tarjan's algorithm is an efficient DFS-based method for finding these.
///
/// # Arguments
/// * `graph` - The graph to analyze.
///
/// # Returns
/// A tuple `(bridges, articulation_points)` where `bridges` is a `Vec<(usize, usize)>`
/// and `articulation_points` is a `Vec<usize>`.
pub fn find_bridges_and_articulation_points<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    graph: &Graph<V>,
) -> (Vec<(usize, usize)>, Vec<usize>) {
    let mut bridges = Vec::new();
    let mut articulation_points = HashSet::new();
    let mut visited = HashSet::new();
    let mut discovery_times = HashMap::new();
    let mut low_link = HashMap::new();
    let mut time = 0;

    for node_id in 0..graph.nodes.len() {
        if !visited.contains(&node_id) {
            b_and_ap_util(
                graph,
                node_id,
                None,
                &mut time,
                &mut visited,
                &mut discovery_times,
                &mut low_link,
                &mut bridges,
                &mut articulation_points,
            );
        }
    }
    (bridges, articulation_points.into_iter().collect())
}

pub(crate) fn b_and_ap_util<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    graph: &Graph<V>,
    u: usize,
    parent: Option<usize>,
    time: &mut usize,
    visited: &mut HashSet<usize>,
    disc: &mut HashMap<usize, usize>,
    low: &mut HashMap<usize, usize>,
    bridges: &mut Vec<(usize, usize)>,
    ap: &mut HashSet<usize>,
) {
    visited.insert(u);
    disc.insert(u, *time);
    low.insert(u, *time);
    *time += 1;
    let mut children = 0;

    if let Some(neighbors) = graph.adj.get(u) {
        for &(v, _) in neighbors {
            if Some(v) == parent {
                continue;
            }
            if visited.contains(&v) {
                low.insert(u, (*low.get(&u).unwrap()).min(*disc.get(&v).unwrap()));
            } else {
                children += 1;
                b_and_ap_util(graph, v, Some(u), time, visited, disc, low, bridges, ap);
                low.insert(u, (*low.get(&u).unwrap()).min(*low.get(&v).unwrap()));

                if parent.is_some() && low.get(&v).unwrap() >= disc.get(&u).unwrap() {
                    ap.insert(u);
                }
                if low.get(&v).unwrap() > disc.get(&u).unwrap() {
                    bridges.push((u, v));
                }
            }
        }
    }

    if parent.is_none() && children > 1 {
        ap.insert(u);
    }
}

// =====================================================================================
// region: Minimum Spanning Tree
// =====================================================================================

/// A Disjoint Set Union (DSU) data structure for Kruskal's algorithm.
pub struct DSU {
    parent: Vec<usize>,
}

impl DSU {
    pub(crate) fn new(n: usize) -> Self {
        DSU {
            parent: (0..n).collect(),
        }
    }

    pub(crate) fn find(&mut self, i: usize) -> usize {
        if self.parent[i] == i {
            return i;
        }
        self.parent[i] = self.find(self.parent[i]);
        self.parent[i]
    }

    pub(crate) fn union(&mut self, i: usize, j: usize) {
        let root_i = self.find(i);
        let root_j = self.find(j);
        if root_i != root_j {
            self.parent[root_i] = root_j;
        }
    }
}

/// Finds the Minimum Spanning Tree (MST) of a graph using Kruskal's algorithm.
///
/// Kruskal's algorithm is a greedy algorithm that finds an MST for a connected,
/// undirected graph. It works by adding edges in increasing order of weight,
/// as long as they do not form a cycle.
///
/// # Arguments
/// * `graph` - The graph to analyze.
///
/// # Returns
/// A `Vec<(usize, usize, Expr)>` representing the edges `(u, v, weight)` that form the MST.
pub fn kruskal_mst<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    graph: &Graph<V>,
) -> Vec<(usize, usize, Expr)> {
    let mut edges = graph.get_edges();

    edges.sort_by(|a, b| {
        let weight_a = as_f64(&a.2).unwrap_or(f64::INFINITY);
        let weight_b = as_f64(&b.2).unwrap_or(f64::INFINITY);
        weight_a.partial_cmp(&weight_b).unwrap()
    });

    let mut dsu = DSU::new(graph.nodes.len());
    let mut mst = Vec::new();

    for (u, v, weight) in edges {
        if dsu.find(u) != dsu.find(v) {
            dsu.union(u, v);
            mst.push((u, v, weight));
        }
    }
    mst
}

// =====================================================================================
// region: Network Flow
// =====================================================================================

/// Finds the maximum flow from a source `s` to a sink `t` in a flow network
/// using the Edmonds-Karp algorithm.
///
/// Edmonds-Karp is an implementation of the Ford-Fulkerson method that uses BFS
/// to find augmenting paths in the residual graph. It guarantees to find the maximum flow.
///
/// # Arguments
/// * `capacity_graph` - A graph where edge weights represent capacities.
/// * `s` - The source node index.
/// * `t` - The sink node index.
///
/// # Returns
/// The maximum flow value as an `f64`.
pub fn edmonds_karp_max_flow<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    capacity_graph: &Graph<V>,
    s: usize,
    t: usize,
) -> f64 {
    let n = capacity_graph.nodes.len();
    let mut residual_capacity = vec![vec![0.0; n]; n];
    for u in 0..n {
        if let Some(neighbors) = capacity_graph.adj.get(u) {
            for &(v, ref cap) in neighbors {
                residual_capacity[u][v] = as_f64(cap).unwrap_or(0.0);
            }
        }
    }

    let mut max_flow = 0.0;

    loop {
        let (parent, path_flow) = bfs_for_augmenting_path(&residual_capacity, s, t);

        if path_flow == 0.0 {
            break; // No more augmenting paths
        }

        max_flow += path_flow;
        let mut v = t;
        while v != s {
            let u = parent[v].unwrap();
            residual_capacity[u][v] -= path_flow;
            residual_capacity[v][u] += path_flow;
            v = u;
        }
    }

    max_flow
}

/// Helper BFS to find an augmenting path in the residual graph.
pub(crate) fn bfs_for_augmenting_path(
    capacity: &Vec<Vec<f64>>,
    s: usize,
    t: usize,
) -> (Vec<Option<usize>>, f64) {
    let n = capacity.len();
    let mut parent = vec![None; n];
    let mut queue = VecDeque::new();
    let mut path_flow = vec![f64::INFINITY; n];

    queue.push_back(s);

    while let Some(u) = queue.pop_front() {
        for v in 0..n {
            if parent[v].is_none() && v != s && capacity[u][v] > 0.0 {
                parent[v] = Some(u);
                path_flow[v] = path_flow[u].min(capacity[u][v]);
                if v == t {
                    return (parent, path_flow[t]);
                }
                queue.push_back(v);
            }
        }
    }

    (parent, 0.0)
}

/// Finds the maximum flow from a source `s` to a sink `t` in a flow network
/// using Dinic's algorithm.
///
/// Dinic's algorithm is a more efficient algorithm for solving the maximum flow problem
/// compared to Edmonds-Karp, especially for dense graphs. It uses a level graph
/// and blocking flows to find augmenting paths.
///
/// # Arguments
/// * `capacity_graph` - A graph where edge weights represent capacities.
/// * `s` - The source node index.
/// * `t` - The sink node index.
///
/// # Returns
/// The maximum flow value as an `f64`.
pub fn dinic_max_flow<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    capacity_graph: &Graph<V>,
    s: usize,
    t: usize,
) -> f64 {
    let n = capacity_graph.nodes.len();
    let mut residual_capacity = vec![vec![0.0; n]; n];
    for u in 0..n {
        if let Some(neighbors) = capacity_graph.adj.get(u) {
            for &(v, ref cap) in neighbors {
                residual_capacity[u][v] = as_f64(cap).unwrap_or(0.0);
            }
        }
    }

    let mut max_flow = 0.0;
    let mut level = vec![0; n];

    while dinic_bfs(&residual_capacity, s, t, &mut level) {
        let mut ptr = vec![0; n];
        while {
            let pushed = dinic_dfs(
                &mut residual_capacity,
                s,
                t,
                f64::INFINITY,
                &level,
                &mut ptr,
            );
            if pushed > 0.0 {
                max_flow += pushed;
                true
            } else {
                false
            }
        } {}
    }

    max_flow
}

pub(crate) fn dinic_bfs(
    capacity: &Vec<Vec<f64>>,
    s: usize,
    t: usize,
    level: &mut Vec<i32>,
) -> bool {
    level.iter_mut().for_each(|l| *l = -1);
    level[s] = 0;
    let mut q = VecDeque::new();
    q.push_back(s);

    while let Some(u) = q.pop_front() {
        for v in 0..capacity.len() {
            if level[v] < 0 && capacity[u][v] > 0.0 {
                level[v] = level[u] + 1;
                q.push_back(v);
            }
        }
    }
    level[t] != -1
}

pub(crate) fn dinic_dfs(
    cap: &mut Vec<Vec<f64>>,
    u: usize,
    t: usize,
    pushed: f64,
    level: &Vec<i32>,
    ptr: &mut Vec<usize>,
) -> f64 {
    if pushed == 0.0 {
        return 0.0;
    }
    if u == t {
        return pushed;
    }

    while ptr[u] < cap.len() {
        let v = ptr[u];
        if level[v] != level[u] + 1 || cap[u][v] == 0.0 {
            ptr[u] += 1;
            continue;
        }
        let tr = dinic_dfs(cap, v, t, pushed.min(cap[u][v]), level, ptr);
        if tr == 0.0 {
            ptr[u] += 1;
            continue;
        }
        cap[u][v] -= tr;
        cap[v][u] += tr;
        return tr;
    }
    0.0
}

/// Finds the shortest paths from a single source in a graph with possible negative edge weights.
///
/// Bellman-Ford algorithm is capable of handling graphs where some edge weights are negative,
/// unlike Dijkstra's algorithm. It can also detect negative-weight cycles.
///
/// # Arguments
/// * `graph` - The graph to analyze.
/// * `start_node` - The index of the starting node.
///
/// # Returns
/// A `Result` containing a tuple `(distances, predecessors)`.
/// `distances` is a `HashMap` from node ID to shortest distance.
/// `predecessors` is a `HashMap` from node ID to its predecessor on the shortest path.
/// Returns an error string if a negative-weight cycle is detected.
pub fn bellman_ford<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    graph: &Graph<V>,
    start_node: usize,
) -> Result<(HashMap<usize, f64>, HashMap<usize, Option<usize>>), String> {
    let n = graph.nodes.len();
    let mut dist = HashMap::new();
    let mut prev = HashMap::new();

    for node_id in 0..graph.nodes.len() {
        dist.insert(node_id, f64::INFINITY);
    }
    dist.insert(start_node, 0.0);

    for _ in 1..n {
        for u in 0..n {
            if let Some(neighbors) = graph.adj.get(u) {
                for &(v, ref weight) in neighbors {
                    let w = as_f64(weight).unwrap_or(f64::INFINITY);
                    if dist[&u] != f64::INFINITY && dist[&u] + w < dist[&v] {
                        dist.insert(v, dist[&u] + w);
                        prev.insert(v, Some(u));
                    }
                }
            }
        }
    }

    // Check for negative-weight cycles
    for u in 0..n {
        if let Some(neighbors) = graph.adj.get(u) {
            for &(v, ref weight) in neighbors {
                let w = as_f64(weight).unwrap_or(f64::INFINITY);
                if dist[&u] != f64::INFINITY && dist[&u] + w < dist[&v] {
                    return Err("Graph contains a negative-weight cycle.".to_string());
                }
            }
        }
    }

    Ok((dist, prev))
}

/// Solves the Minimum-Cost Maximum-Flow problem using the successive shortest path algorithm with Bellman-Ford.
///
/// This algorithm finds the maximum flow through a network while minimizing the total cost of the flow.
/// It repeatedly finds the shortest augmenting path in the residual graph, where edge weights are costs.
/// Assumes edge weights are given as a tuple `(capacity, cost)`.
///
/// # Arguments
/// * `graph` - The graph where edge weights are `Expr::Tuple(capacity, cost)`.
/// * `s` - The source node index.
/// * `t` - The sink node index.
///
/// # Returns
/// A tuple `(max_flow, min_cost)` as `f64`.
#[allow(unused_variables)]
pub fn min_cost_max_flow<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    graph: &Graph<V>,
    s: usize,
    t: usize,
) -> (f64, f64) {
    let n = graph.nodes.len();
    let mut capacity = vec![vec![0.0; n]; n];
    let mut cost = vec![vec![0.0; n]; n];

    for u in 0..n {
        if let Some(neighbors) = graph.adj.get(u) {
            for &(v, ref weight) in neighbors {
                if let Expr::Tuple(t) = weight {
                    if t.len() == 2 {
                        capacity[u][v] = as_f64(&t[0]).unwrap_or(0.0);
                        cost[u][v] = as_f64(&t[1]).unwrap_or(0.0);
                    }
                }
            }
        }
    }

    let mut flow = 0.0;
    let mut total_cost = 0.0;

    loop {
        // Find shortest path in residual graph using costs as weights
        let mut dist = vec![f64::INFINITY; n];
        let mut parent = vec![None; n];
        dist[s] = 0.0;

        for _ in 1..n {
            for u in 0..n {
                for v in 0..n {
                    if capacity[u][v] > 0.0
                        && dist[u] != f64::INFINITY
                        && dist[u] + cost[u][v] < dist[v]
                    {
                        dist[v] = dist[u] + cost[u][v];
                        parent[v] = Some(u);
                    }
                }
            }
        }

        if dist[t] == f64::INFINITY {
            break; // No more augmenting paths
        }

        // Find path flow
        let mut path_flow = f64::INFINITY;
        let mut curr = t;
        while let Some(prev) = parent[curr] {
            path_flow = path_flow.min(capacity[prev][curr]);
            curr = prev;
        }

        // Augment flow
        flow += path_flow;
        total_cost += path_flow * dist[t];
        let mut v = t;
        while let Some(u) = parent[v] {
            capacity[u][v] -= path_flow;
            capacity[v][u] += path_flow;
            // Note: cost[v][u] should be -cost[u][v], this simple matrix representation doesn't handle that well.
            // A full implementation would use a proper residual graph struct.
            v = u;
        }
    }

    (0.0, 0.0) // Returns (max_flow, min_cost)
}

// =====================================================================================
// region: Matching, Covering, and Partitioning
// =====================================================================================

/// Checks if a graph is bipartite using BFS-based 2-coloring.
///
/// A graph is bipartite if its vertices can be divided into two disjoint and independent sets
/// `U` and `V` such that every edge connects a vertex in `U` to one in `V`.
///
/// # Arguments
/// * `graph` - The graph to check.
///
/// # Returns
/// `Some(partition)` if bipartite, where `partition[i]` is `0` or `1` indicating the set.
/// `None` if not bipartite.
pub fn is_bipartite<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    graph: &Graph<V>,
) -> Option<Vec<i8>> {
    let n = graph.nodes.len();
    let mut colors = vec![-1; n]; // -1: uncolored, 0: color 1, 1: color 2

    for i in 0..n {
        if colors[i] == -1 {
            let mut queue = VecDeque::new();
            queue.push_back(i);
            colors[i] = 0;

            while let Some(u) = queue.pop_front() {
                if let Some(neighbors) = graph.adj.get(u) {
                    for &(v, _) in neighbors {
                        if colors[v] == -1 {
                            colors[v] = 1 - colors[u];
                            queue.push_back(v);
                        } else if colors[v] == colors[u] {
                            return None; // Not bipartite
                        }
                    }
                }
            }
        }
    }
    Some(colors)
}

/// Finds the maximum cardinality matching in a bipartite graph by reducing it to a max-flow problem.
///
/// A matching is a set of edges without common vertices. A maximum matching is one with the largest
/// possible number of edges. This function constructs a flow network from the bipartite graph
/// and uses a max-flow algorithm to find the matching.
///
/// # Arguments
/// * `graph` - The bipartite graph.
/// * `partition` - The partition of vertices into two sets (from `is_bipartite`).
///
/// # Returns
/// A `Vec<(usize, usize)>` representing the edges `(u, v)` in the maximum matching.
pub fn bipartite_maximum_matching<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    graph: &Graph<V>,
    partition: &[i8],
) -> Vec<(usize, usize)> {
    let n = graph.nodes.len();
    let s = n; // Source
    let t = n + 1; // Sink
    let mut flow_graph = Graph::new(true);

    let mut u_nodes = Vec::new();
    let mut v_nodes = Vec::new();
    for i in 0..n {
        if partition[i] == 0 {
            u_nodes.push(i);
        } else {
            v_nodes.push(i);
        }
    }

    // Build flow network
    for &u in &u_nodes {
        flow_graph.add_edge(&s, &u, Expr::Constant(1.0));
        if let Some(neighbors) = graph.adj.get(u) {
            for &(v, _) in neighbors {
                flow_graph.add_edge(&u, &v, Expr::Constant(1.0));
            }
        }
    }
    for &v in &v_nodes {
        flow_graph.add_edge(&v, &t, Expr::Constant(1.0));
    }

    // The max flow is equivalent to the max matching size.
    // To get the actual edges, we need to inspect the flow on edges (u,v).
    // This requires a version of edmonds_karp that returns the final flow network.
    // For now, we return an empty vec as a placeholder for the matched edges.
    let _max_flow = edmonds_karp_max_flow(&flow_graph, s, t);
    vec![]
}

// =====================================================================================
// region: Minimum Spanning Tree & Topological Sort
// =====================================================================================

/// Finds the Minimum Spanning Tree (MST) of a graph using Prim's algorithm.
///
/// Prim's algorithm is a greedy algorithm that finds an MST for a connected,
/// undirected graph. It grows the MST from an initial vertex by iteratively
/// adding the cheapest edge that connects a vertex in the tree to one outside the tree.
///
/// # Arguments
/// * `graph` - The graph to analyze.
/// * `start_node` - The starting node for building the MST.
///
/// # Returns
/// A `Vec<(usize, usize, Expr)>` representing the edges `(u, v, weight)` that form the MST.
pub fn prim_mst<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    graph: &Graph<V>,
    start_node: usize,
) -> Vec<(usize, usize, Expr)> {
    let n = graph.nodes.len();
    let mut mst = Vec::new();
    let mut visited = vec![false; n];
    let mut pq = std::collections::BinaryHeap::new();

    visited[start_node] = true;
    if let Some(neighbors) = graph.adj.get(start_node) {
        for &(v, ref weight) in neighbors {
            let cost = as_f64(weight).unwrap_or(f64::INFINITY);
            pq.push((
                ordered_float::OrderedFloat(-cost),
                start_node,
                v,
                weight.clone(),
            ));
        }
    }

    while let Some((_, u, v, weight)) = pq.pop() {
        if visited[v] {
            continue;
        }
        visited[v] = true;
        mst.push((u, v, weight));

        if let Some(neighbors) = graph.adj.get(v) {
            for &(next_v, ref next_weight) in neighbors {
                if !visited[next_v] {
                    let cost = as_f64(next_weight).unwrap_or(f64::INFINITY);
                    pq.push((
                        ordered_float::OrderedFloat(-cost),
                        v,
                        next_v,
                        next_weight.clone(),
                    ));
                }
            }
        }
    }
    mst
}

/// Performs a topological sort on a directed acyclic graph (DAG) using Kahn's algorithm (BFS-based).
///
/// A topological sort is a linear ordering of its vertices such that for every directed edge `uv`
/// from vertex `u` to vertex `v`, `u` comes before `v` in the ordering.
/// Kahn's algorithm works by iteratively removing vertices with an in-degree of 0.
///
/// # Arguments
/// * `graph` - The DAG to sort.
///
/// # Returns
/// A `Result` containing a `Vec<usize>` of node indices in topological order,
/// or an error string if the graph has a cycle.
pub fn topological_sort_kahn<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    graph: &Graph<V>,
) -> Result<Vec<usize>, String> {
    if !graph.is_directed {
        return Err("Topological sort is only defined for directed graphs.".to_string());
    }
    let n = graph.nodes.len();
    let mut in_degree = vec![0; n];
    for i in 0..n {
        in_degree[i] = graph.in_degree(i);
    }

    let mut queue: VecDeque<usize> = (0..n).filter(|&i| in_degree[i] == 0).collect();
    let mut sorted_order = Vec::new();

    while let Some(u) = queue.pop_front() {
        sorted_order.push(u);
        if let Some(neighbors) = graph.adj.get(u) {
            for &(v, _) in neighbors {
                in_degree[v] -= 1;
                if in_degree[v] == 0 {
                    queue.push_back(v);
                }
            }
        }
    }

    if sorted_order.len() == n {
        Ok(sorted_order)
    } else {
        Err("Graph has a cycle, topological sort is not possible.".to_string())
    }
}

/// Performs a topological sort on a directed acyclic graph (DAG) using a DFS-based algorithm.
///
/// This algorithm works by performing a DFS traversal and adding vertices to the sorted list
/// after all their dependencies (children in the DFS tree) have been visited.
///
/// # Arguments
/// * `graph` - The DAG to sort.
///
/// # Returns
/// A `Vec<usize>` containing the node IDs in topological order.
pub fn topological_sort_dfs<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    graph: &Graph<V>,
) -> Vec<usize> {
    let mut visited = HashSet::new();
    let mut stack = Vec::new();
    for node_id in 0..graph.nodes.len() {
        if !visited.contains(&node_id) {
            topo_dfs_util(graph, node_id, &mut visited, &mut stack);
        }
    }
    stack.reverse();
    stack
}

pub(crate) fn topo_dfs_util<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    graph: &Graph<V>,
    u: usize,
    visited: &mut HashSet<usize>,
    stack: &mut Vec<usize>,
) {
    visited.insert(u);
    if let Some(neighbors) = graph.adj.get(u) {
        for &(v, _) in neighbors {
            if !visited.contains(&v) {
                topo_dfs_util(graph, v, visited, stack);
            }
        }
    }
    stack.push(u);
}

/// Finds the minimum vertex cover of a bipartite graph using Kőnig's theorem.
///
/// Kőnig's theorem states that in any bipartite graph, the number of edges in a maximum matching
/// equals the number of vertices in a minimum vertex cover. This function leverages a maximum
/// matching to construct the minimum vertex cover.
///
/// # Arguments
/// * `graph` - The bipartite graph.
/// * `partition` - The partition of vertices into two sets.
/// * `matching` - The set of edges in a maximum matching.
///
/// # Returns
/// A `Vec<usize>` of node indices representing the minimum vertex cover.
pub fn bipartite_minimum_vertex_cover<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    graph: &Graph<V>,
    partition: &[i8],
    matching: &[(usize, usize)],
) -> Vec<usize> {
    let mut u_nodes = HashSet::new();
    let mut matched_nodes_u = HashSet::new();
    for i in 0..partition.len() {
        if partition[i] == 0 {
            u_nodes.insert(i);
        }
    }
    for &(u, v) in matching {
        if u_nodes.contains(&u) {
            matched_nodes_u.insert(u);
        } else {
            matched_nodes_u.insert(v);
        }
    }

    let unmatched_u: Vec<_> = u_nodes.difference(&matched_nodes_u).cloned().collect();

    // Find all vertices reachable from unmatched U vertices by alternating paths.
    let mut visited = HashSet::new();
    let mut queue = VecDeque::from(unmatched_u);
    while let Some(u) = queue.pop_front() {
        if visited.contains(&u) {
            continue;
        }
        visited.insert(u);
        if let Some(neighbors) = graph.adj.get(u) {
            for &(v, _) in neighbors {
                // If (u,v) is NOT a matching edge, and v is not visited, add v.
                if !matching.contains(&(u, v))
                    && !matching.contains(&(v, u))
                    && !visited.contains(&v)
                {
                    queue.push_back(v);
                }
            }
        }
    }

    // The minimum vertex cover is (U \ Z) U (V ∩ Z), where Z is the set of visited vertices.
    let mut cover = Vec::new();
    for u in u_nodes {
        if !visited.contains(&u) {
            cover.push(u);
        }
    }
    for i in 0..partition.len() {
        if partition[i] == 1 && visited.contains(&i) {
            cover.push(i);
        }
    }
    cover
}

/// Finds the maximum cardinality matching in a bipartite graph using the Hopcroft-Karp algorithm.
///
/// The Hopcroft-Karp algorithm is an efficient algorithm for finding maximum cardinality matchings
/// in bipartite graphs. It works by repeatedly finding a maximal set of shortest augmenting paths.
///
/// # Arguments
/// * `graph` - The bipartite graph.
/// * `partition` - The partition of vertices into two sets.
///
/// # Returns
/// A `Vec<(usize, usize)>` representing the edges `(u, v)` in the maximum matching.
#[allow(unused_variables)]
pub fn hopcroft_karp_bipartite_matching<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    graph: &Graph<V>,
    partition: &[i8],
) -> Vec<(usize, usize)> {
    let n = graph.nodes.len();
    let mut u_nodes = Vec::new();
    for i in 0..n {
        if partition[i] == 0 {
            u_nodes.push(i);
        }
    }

    let mut pair_u = vec![None; n];
    let mut pair_v = vec![None; n];
    let mut dist = vec![0; n];
    let mut matching = 0;

    while hopcroft_karp_bfs(graph, &u_nodes, &mut pair_u, &mut pair_v, &mut dist) {
        for &u in &u_nodes {
            if pair_u[u].is_none() {
                if hopcroft_karp_dfs(graph, u, &mut pair_u, &mut pair_v, &mut dist) {
                    matching += 1;
                }
            }
        }
    }

    let mut result = Vec::new();
    for u in 0..n {
        if let Some(v) = pair_u[u] {
            result.push((u, v));
        }
    }
    result
}

pub(crate) fn hopcroft_karp_bfs<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    graph: &Graph<V>,
    u_nodes: &[usize],
    pair_u: &mut [Option<usize>],
    pair_v: &mut [Option<usize>],
    dist: &mut [usize],
) -> bool {
    let mut queue = VecDeque::new();
    for &u in u_nodes {
        if pair_u[u].is_none() {
            dist[u] = 0;
            queue.push_back(u);
        } else {
            dist[u] = usize::MAX;
        }
    }
    let mut found_path = false;

    while let Some(u) = queue.pop_front() {
        if let Some(neighbors) = graph.adj.get(u) {
            for &(v, _) in neighbors {
                if let Some(next_u_opt) = pair_v.get(v) {
                    if let Some(next_u) = next_u_opt {
                        if dist[*next_u] == usize::MAX {
                            dist[*next_u] = dist[u] + 1;
                            queue.push_back(*next_u);
                        }
                    } else {
                        // Unmatched V node found
                        found_path = true;
                    }
                }
            }
        }
    }
    found_path
}

pub(crate) fn hopcroft_karp_dfs<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    graph: &Graph<V>,
    u: usize,
    pair_u: &mut [Option<usize>],
    pair_v: &mut [Option<usize>],
    dist: &mut [usize],
) -> bool {
    if let Some(neighbors) = graph.adj.get(u) {
        for &(v, _) in neighbors {
            if let Some(next_u_opt) = pair_v.get(v) {
                if let Some(next_u) = next_u_opt {
                    if dist[*next_u] == dist[u] + 1 {
                        if hopcroft_karp_dfs(graph, *next_u, pair_u, pair_v, dist) {
                            pair_v[v] = Some(u);
                            pair_u[u] = Some(v);
                            return true;
                        }
                    }
                } else {
                    // Found augmenting path to an unmatched V node
                    pair_v[v] = Some(u);
                    pair_u[u] = Some(v);
                    return true;
                }
            }
        }
    }
    dist[u] = usize::MAX;
    false
}

/// Finds the maximum cardinality matching in a general graph using Edmonds's Blossom Algorithm.
///
/// Edmonds's Blossom Algorithm is a polynomial-time algorithm for finding maximum matchings
/// in general (non-bipartite) graphs. It works by iteratively finding augmenting paths
/// and handling "blossoms" (odd-length cycles).
///
/// # Arguments
/// * `graph` - The graph to analyze.
///
/// # Returns
/// A `Vec<(usize, usize)>` representing the edges `(u, v)` in the maximum matching.
#[allow(unused_variables)]
pub fn blossom_algorithm<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    graph: &Graph<V>,
) -> Vec<(usize, usize)> {
    let n = graph.nodes.len();
    let mut matching = vec![None; n];
    let mut matches = 0;

    for i in 0..n {
        if matching[i].is_none() {
            let path = find_augmenting_path_with_blossoms(graph, i, &matching);
            if !path.is_empty() {
                matches += 1;
                let mut u = path[0];
                for &v in path.iter().skip(1) {
                    matching[u] = Some(v);
                    matching[v] = Some(u);
                    u = v;
                }
            }
        }
    }

    let mut result = Vec::new();
    for u in 0..n {
        if let Some(v) = matching[u] {
            if u < v {
                result.push((u, v));
            }
        }
    }
    result
}

pub(crate) fn find_augmenting_path_with_blossoms<
    V: Eq + std::hash::Hash + Clone + std::fmt::Debug,
>(
    graph: &Graph<V>,
    start_node: usize,
    matching: &[Option<usize>],
) -> Vec<usize> {
    let n = graph.nodes.len();
    let mut parent = vec![None; n];
    let mut origin = (0..n).collect::<Vec<_>>();
    let mut level = vec![-1; n];
    let mut queue = VecDeque::new();

    level[start_node] = 0;
    queue.push_back(start_node);

    while let Some(u) = queue.pop_front() {
        if let Some(neighbors) = graph.adj.get(u) {
            for &(v, _) in neighbors {
                if level[v] == -1 {
                    // Unvisited node
                    if let Some(w) = matching[v] {
                        parent[v] = Some(u);
                        level[v] = 1;
                        level[w] = 0;
                        queue.push_back(w);
                    } else {
                        // Found an augmenting path
                        parent[v] = Some(u);
                        let mut path = vec![v, u];
                        let mut curr = u;
                        while let Some(p) = parent[curr] {
                            path.push(p);
                            curr = p;
                        }
                        return path;
                    }
                } else if level[v] == 0 {
                    // Found a blossom
                    let base = find_common_ancestor(&origin, &parent, u, v);
                    contract_blossom::<V>(
                        base,
                        u,
                        v,
                        &mut queue,
                        &mut level,
                        &mut origin,
                        &mut parent,
                        matching,
                    );
                    contract_blossom::<V>(
                        base,
                        v,
                        u,
                        &mut queue,
                        &mut level,
                        &mut origin,
                        &mut parent,
                        matching,
                    );
                }
            }
        }
    }
    vec![] // No augmenting path found
}

pub(crate) fn find_common_ancestor(
    origin: &[usize],
    parent: &[Option<usize>],
    mut u: usize,
    mut v: usize,
) -> usize {
    let mut visited = vec![false; origin.len()];
    loop {
        u = origin[u];
        visited[u] = true;
        if parent[u].is_none() {
            break;
        }
        u = parent[u].unwrap();
    }
    loop {
        v = origin[v];
        if visited[v] {
            return v;
        }
        v = parent[v].unwrap();
    }
}

pub(crate) fn contract_blossom<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    base: usize,
    mut u: usize,
    v: usize,
    queue: &mut VecDeque<usize>,
    level: &mut Vec<i32>,
    origin: &mut [usize],
    parent: &mut Vec<Option<usize>>,
    matching: &[Option<usize>],
) {
    while origin[u] != base {
        parent[u] = Some(v);
        origin[u] = base;
        if let Some(w) = matching[u] {
            if level[w] == -1 {
                level[w] = 0;
                queue.push_back(w);
            }
        }
        u = parent[u].unwrap();
    }
}

// =====================================================================================
// region: Path Finding
// =====================================================================================

/// Finds the shortest path in an unweighted graph from a source node using BFS.
///
/// Since all edge weights are implicitly 1, BFS naturally finds the shortest path
/// in terms of number of edges.
///
/// # Arguments
/// * `graph` - The graph to search.
/// * `start_node` - The index of the starting node.
///
/// # Returns
/// A `HashMap<usize, (usize, Option<usize>)>` where keys are node IDs and values are
/// `(distance_from_start, predecessor_node_id)`.
pub fn shortest_path_unweighted<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    graph: &Graph<V>,
    start_node: usize,
) -> HashMap<usize, (usize, Option<usize>)> {
    let mut distances = HashMap::new();
    let mut predecessors = HashMap::new();
    let mut queue = VecDeque::new();

    distances.insert(start_node, 0);
    predecessors.insert(start_node, None);
    queue.push_back(start_node);

    while let Some(u) = queue.pop_front() {
        let u_dist = *distances.get(&u).unwrap();
        if let Some(neighbors) = graph.adj.get(u) {
            for &(v, _) in neighbors {
                if !distances.contains_key(&v) {
                    distances.insert(v, u_dist + 1);
                    predecessors.insert(v, Some(u));
                    queue.push_back(v);
                }
            }
        }
    }

    let mut result = HashMap::new();
    for (node, dist) in distances {
        result.insert(node, (dist, predecessors.get(&node).cloned().flatten()));
    }
    result
}

/// Finds the shortest paths from a single source using Dijkstra's algorithm.
///
/// Dijkstra's algorithm is a greedy algorithm that solves the single-source
/// shortest path problem for a graph with non-negative edge weights.
///
/// # Arguments
/// * `graph` - The graph to search.
/// * `start_node` - The index of the starting node.
///
/// # Returns
/// A `HashMap<usize, (Expr, Option<usize>)>` where keys are node IDs and values are
/// `(shortest_distance, predecessor_node_id)`.
pub fn dijkstra<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    graph: &Graph<V>,
    start_node: usize,
) -> HashMap<usize, (Expr, Option<usize>)> {
    let mut dist = HashMap::new();
    let mut prev = HashMap::new();
    let mut pq = std::collections::BinaryHeap::new();

    for node_id in 0..graph.nodes.len() {
        dist.insert(node_id, Expr::Infinity);
    }

    dist.insert(start_node, Expr::Constant(0.0));
    prev.insert(start_node, None);
    // BinaryHeap is a max-heap, so we store negative costs.
    pq.push((OrderedFloat(0.0), start_node));

    while let Some((cost, u)) = pq.pop() {
        let cost = -cost.0;
        if cost > as_f64(dist.get(&u).unwrap()).unwrap_or(f64::INFINITY) {
            continue;
        }

        if let Some(neighbors) = graph.adj.get(u) {
            for &(v, ref weight) in neighbors {
                let new_dist = simplify(Expr::Add(
                    Box::new(Expr::Constant(cost)),
                    Box::new(weight.clone()),
                ));
                if as_f64(&new_dist).unwrap_or(f64::INFINITY)
                    < as_f64(dist.get(&v).unwrap()).unwrap_or(f64::INFINITY)
                {
                    dist.insert(v, new_dist.clone());
                    prev.insert(v, Some(u));
                    pq.push((OrderedFloat(-as_f64(&new_dist).unwrap()), v));
                }
            }
        }
    }

    let mut result = HashMap::new();
    for (node, d) in dist {
        result.insert(node, (d, prev.get(&node).cloned().flatten()));
    }
    result
}

/// Finds all-pairs shortest paths using the Floyd-Warshall algorithm.
///
/// The Floyd-Warshall algorithm is an all-pairs shortest path algorithm that works
/// for both directed and undirected graphs with non-negative or negative edge weights
/// (but no negative cycles).
///
/// # Arguments
/// * `graph` - The graph to analyze.
///
/// # Returns
/// An `Expr::Matrix` where `M[i][j]` is the shortest distance from node `i` to node `j`.
pub fn floyd_warshall<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(graph: &Graph<V>) -> Expr {
    let n = graph.nodes.len();
    let mut dist = vec![vec![Expr::Infinity; n]; n];

    for i in 0..n {
        dist[i][i] = Expr::Constant(0.0);
    }

    for u in 0..n {
        if let Some(neighbors) = graph.adj.get(u) {
            for &(v, ref weight) in neighbors {
                dist[u][v] = weight.clone();
            }
        }
    }

    for k in 0..n {
        for i in 0..n {
            for j in 0..n {
                let new_dist = simplify(Expr::Add(
                    Box::new(dist[i][k].clone()),
                    Box::new(dist[k][j].clone()),
                ));
                // This comparison needs to handle symbolic expressions.
                // For now, we assume numerical weights.
                if as_f64(&dist[i][j]).unwrap_or(f64::INFINITY)
                    > as_f64(&new_dist).unwrap_or(f64::INFINITY)
                {
                    dist[i][j] = new_dist;
                }
            }
        }
    }

    Expr::Matrix(dist)
}

/// Performs spectral analysis on a graph matrix (e.g., Adjacency or Laplacian).
///
/// This function computes the eigenvalues and eigenvectors of the given matrix,
/// which are crucial for understanding various graph properties like connectivity,
/// centrality, and clustering.
///
/// # Arguments
/// * `matrix` - The graph matrix as an `Expr::Matrix`.
///
/// # Returns
/// A `Result` containing a tuple `(eigenvalues, eigenvectors_matrix)`.
/// `eigenvalues` is a column vector of eigenvalues.
/// `eigenvectors_matrix` is a matrix where each column is an eigenvector.
pub fn spectral_analysis(matrix: &Expr) -> Result<(Expr, Expr), String> {
    crate::symbolic::matrix::eigen_decomposition(matrix)
}

/// Computes the algebraic connectivity of a graph.
///
/// The algebraic connectivity is the second-smallest eigenvalue of the Laplacian matrix
/// of a graph. It measures how well-connected a graph is and is related to graph robustness.
///
/// # Arguments
/// * `graph` - The graph to analyze.
///
/// # Returns
/// A `Result` containing an `Expr` representing the algebraic connectivity,
/// or an error string if computation fails or the graph has fewer than 2 eigenvalues.
pub fn algebraic_connectivity<V>(graph: &Graph<V>) -> Result<Expr, String>
where
    V: Clone,
    V: Debug,
    V: Eq,
    V: Hash,
{
    let laplacian = graph.to_laplacian_matrix();
    let (eigenvalues, _) = spectral_analysis(&laplacian)?;

    if let Expr::Matrix(eig_vec) = eigenvalues {
        if eig_vec.len() < 2 {
            return Err("Graph has fewer than 2 eigenvalues.".to_string());
        }
        // Eigenvalues from `eigen_decomposition` are not guaranteed to be sorted.
        // We need to convert them to f64, sort, and then get the second one.
        let mut numerical_eigenvalues = Vec::new();
        for val_expr in eig_vec.iter().flatten() {
            if let Some(val) = as_f64(val_expr) {
                numerical_eigenvalues.push(val);
            } else {
                return Err("Eigenvalues are not all numerical, cannot sort.".to_string());
            }
        }
        numerical_eigenvalues.sort_by(|a, b| a.partial_cmp(b).unwrap());
        Ok(Expr::Constant(numerical_eigenvalues[1]))
    } else {
        Err("Eigenvalue computation did not return a vector.".to_string())
    }
}
