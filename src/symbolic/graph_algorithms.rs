//! # Comprehensive Graph Algorithms
//!
//! This module provides a comprehensive suite of graph algorithms for various tasks,
//! including graph traversal (DFS, BFS), connectivity analysis (connected components,
//! strongly connected components), cycle detection, minimum spanning trees (Kruskal's,
//! Prim's), network flow (Edmonds-Karp, Dinic's), shortest paths (Dijkstra's,
//! Bellman-Ford, Floyd-Warshall), and topological sorting.
//!
//! **Note on Symbolic Computation**: This module maintains symbolic expressions throughout
//! computations where possible. For algorithms that require numeric comparisons (e.g., shortest
//! path, max flow), we provide a `try_numeric_value` helper that attempts to extract numeric
//! values from symbolic expressions when needed for comparison, while preserving the symbolic
//! structure in results.

use crate::symbolic::core::Expr;
use crate::symbolic::graph::Graph;
use crate::symbolic::simplify_dag::simplify;
use num_traits::ToPrimitive;
use ordered_float::OrderedFloat;
use std::collections::{HashMap, HashSet, VecDeque};
use std::fmt::Debug;
use std::hash::Hash;
use std::sync::Arc;

/// Helper function to extract a numeric value from a symbolic expression for comparison purposes.
/// This is used internally by algorithms that need to compare weights.
/// Returns None if the expression cannot be evaluated to a number.

fn try_numeric_value(expr: &Expr) -> Option<f64> {

    match expr {
        Expr::Constant(val) => Some(*val),
        Expr::BigInt(val) => val.to_f64(),
        Expr::Rational(val) => val.to_f64(),
        Expr::Add(lhs, rhs) => {

            let l = try_numeric_value(lhs)?;

            let r = try_numeric_value(rhs)?;

            Some(l + r)
        }
        Expr::Sub(lhs, rhs) => {

            let l = try_numeric_value(lhs)?;

            let r = try_numeric_value(rhs)?;

            Some(l - r)
        }
        Expr::Mul(lhs, rhs) => {

            let l = try_numeric_value(lhs)?;

            let r = try_numeric_value(rhs)?;

            Some(l * r)
        }
        Expr::Div(lhs, rhs) => {

            let l = try_numeric_value(lhs)?;

            let r = try_numeric_value(rhs)?;

            if r == 0.0 {

                None
            } else {

                Some(l / r)
            }
        }
        Expr::Neg(inner) => {

            let v = try_numeric_value(inner)?;

            Some(-v)
        }
        Expr::AddList(list) => {

            let mut sum = 0.0;

            for e in list {

                sum += try_numeric_value(e)?;
            }

            Some(sum)
        }
        Expr::MulList(list) => {

            let mut prod = 1.0;

            for e in list {

                prod *= try_numeric_value(e)?;
            }

            Some(prod)
        }
        Expr::Dag(node) => {
            if let Ok(inner) = node.to_expr() {

                try_numeric_value(&inner)
            } else {

                None
            }
        }
        // Try to simplify and extract
        _ => {

            let simplified = simplify(&expr.clone());

            match simplified {
                Expr::Constant(val) => Some(val),
                Expr::BigInt(val) => val.to_f64(),
                Expr::Rational(val) => val.to_f64(),
                _ => None,
            }
        }
    }
}

/// Symbolically adds two expressions.

fn symbolic_add(a: &Expr, b: &Expr) -> Expr {

    Expr::Add(Arc::new(a.clone()), Arc::new(b.clone()))
}

/// Symbolically compares two expressions for ordering.
/// Returns Some(Ordering) if a numeric comparison can be made, None otherwise.
#[allow(dead_code)]

fn symbolic_compare(a: &Expr, b: &Expr) -> Option<std::cmp::Ordering> {

    let a_val = try_numeric_value(a)?;

    let b_val = try_numeric_value(b)?;

    a_val.partial_cmp(&b_val)
}

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
#[must_use]

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
#[must_use]

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
#[must_use]

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
#[must_use]

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
#[must_use]

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

                if let (Some(&low_u), Some(&low_v)) = (low.get(&u), low.get(&v)) {

                    low.insert(u, low_u.min(low_v));
                }
            } else if on_stack.contains(&v) {

                if let (Some(&low_u), Some(&disc_v)) = (low.get(&u), disc.get(&v)) {

                    low.insert(u, low_u.min(disc_v));
                }
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
#[must_use]

pub fn has_cycle<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(graph: &Graph<V>) -> bool {

    let mut visited = HashSet::new();

    if graph.is_directed {

        let mut recursion_stack = HashSet::new();

        for node_id in 0..graph.nodes.len() {

            if !visited.contains(&node_id)
                && has_cycle_directed_util(graph, node_id, &mut visited, &mut recursion_stack)
            {

                return true;
            }
        }
    } else {

        for node_id in 0..graph.nodes.len() {

            if !visited.contains(&node_id)
                && has_cycle_undirected_util(graph, node_id, &mut visited, None)
            {

                return true;
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
#[must_use]

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

    (
        bridges,
        articulation_points
            .into_iter()
            .collect(),
    )
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

                if let (Some(&low_u), Some(&disc_v)) = (low.get(&u), disc.get(&v)) {

                    low.insert(u, low_u.min(disc_v));
                }
            } else {

                children += 1;

                b_and_ap_util(graph, v, Some(u), time, visited, disc, low, bridges, ap);

                if let (Some(&low_u), Some(&low_v)) = (low.get(&u), low.get(&v)) {

                    low.insert(u, low_u.min(low_v));
                }

                if parent.is_some() {

                    if let (Some(&low_v), Some(&disc_u)) = (low.get(&v), disc.get(&u)) {

                        if low_v >= disc_u {

                            ap.insert(u);
                        }
                    }
                }

                if let (Some(&low_v), Some(&disc_u)) = (low.get(&v), disc.get(&u)) {

                    if low_v > disc_u {

                        bridges.push((u, v));
                    }
                }
            }
        }
    }

    if parent.is_none() && children > 1 {

        ap.insert(u);
    }
}

/// A Disjoint Set Union (DSU) data structure for Kruskal's algorithm.

pub struct DSU {
    parent: Vec<usize>,
}

impl DSU {
    pub(crate) fn new(n: usize) -> Self {

        Self {
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
#[must_use]

pub fn kruskal_mst<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    graph: &Graph<V>,
) -> Vec<(usize, usize, Expr)> {

    let mut edges = graph.get_edges();

    edges.sort_by(|a, b| {

        let weight_a = try_numeric_value(&a.2).unwrap_or(f64::INFINITY);

        let weight_b = try_numeric_value(&b.2).unwrap_or(f64::INFINITY);

        weight_a
            .partial_cmp(&weight_b)
            .unwrap_or(std::cmp::Ordering::Equal)
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
#[must_use]

pub fn edmonds_karp_max_flow<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    capacity_graph: &Graph<V>,
    s: usize,
    t: usize,
) -> f64 {

    let n = capacity_graph
        .nodes
        .len();

    let mut residual_capacity = vec![vec![0.0; n]; n];

    for u in 0..n {

        if let Some(neighbors) = capacity_graph
            .adj
            .get(u)
        {

            for &(v, ref cap) in neighbors {

                residual_capacity[u][v] = try_numeric_value(cap).unwrap_or(0.0);
            }
        }
    }

    let mut max_flow = 0.0;

    loop {

        let (parent, path_flow) = bfs_for_augmenting_path(&residual_capacity, s, t);

        if path_flow == 0.0 {

            break;
        }

        max_flow += path_flow;

        let mut v = t;

        while v != s {

            if let Some(u) = parent[v] {

                residual_capacity[u][v] -= path_flow;

                residual_capacity[v][u] += path_flow;

                v = u;
            } else {

                break;
            }
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
#[must_use]

pub fn dinic_max_flow<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    capacity_graph: &Graph<V>,
    s: usize,
    t: usize,
) -> f64 {

    let n = capacity_graph
        .nodes
        .len();

    let mut residual_capacity = vec![vec![0.0; n]; n];

    for u in 0..n {

        if let Some(neighbors) = capacity_graph
            .adj
            .get(u)
        {

            for &(v, ref cap) in neighbors {

                residual_capacity[u][v] = try_numeric_value(cap).unwrap_or(0.0);
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

    for l in level.iter_mut() {

        *l = -1;
    }

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
/// **Note**: This function maintains symbolic expressions for distances where possible.
/// Numeric comparisons are used internally, but the result preserves symbolic structure.
///
/// # Arguments
/// * `graph` - The graph to analyze.
/// * `start_node` - The index of the starting node.
///
/// # Returns
/// A `Result` containing a tuple `(distances, predecessors)`.
/// `distances` is a `HashMap` from node ID to shortest distance (as `Expr`).
/// `predecessors` is a `HashMap` from node ID to its predecessor on the shortest path.
/// Returns an error string if a negative-weight cycle is detected.

pub fn bellman_ford<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    graph: &Graph<V>,
    start_node: usize,
) -> Result<(HashMap<usize, Expr>, HashMap<usize, Option<usize>>), String> {

    let n = graph.nodes.len();

    let mut dist: HashMap<usize, Expr> = HashMap::new();

    let mut prev = HashMap::new();

    let infinity = Expr::Infinity;

    for node_id in 0..graph.nodes.len() {

        dist.insert(node_id, infinity.clone());
    }

    dist.insert(start_node, Expr::Constant(0.0));

    for _ in 1..n {

        for u in 0..n {

            if let Some(neighbors) = graph.adj.get(u) {

                for &(v, ref weight) in neighbors {

                    let dist_u = &dist[&u];

                    let dist_v = &dist[&v];

                    // Check if dist[u] + weight < dist[v]
                    let new_dist = symbolic_add(dist_u, weight);

                    // For comparison, we need numeric values
                    let dist_u_val = try_numeric_value(dist_u).unwrap_or(f64::INFINITY);

                    let weight_val = try_numeric_value(weight).unwrap_or(f64::INFINITY);

                    let dist_v_val = try_numeric_value(dist_v).unwrap_or(f64::INFINITY);

                    if dist_u_val != f64::INFINITY && dist_u_val + weight_val < dist_v_val {

                        dist.insert(v, new_dist);

                        prev.insert(v, Some(u));
                    }
                }
            }
        }
    }

    // Check for negative cycles
    for u in 0..n {

        if let Some(neighbors) = graph.adj.get(u) {

            for &(v, ref weight) in neighbors {

                let dist_u_val = try_numeric_value(&dist[&u]).unwrap_or(f64::INFINITY);

                let weight_val = try_numeric_value(weight).unwrap_or(f64::INFINITY);

                let dist_v_val = try_numeric_value(&dist[&v]).unwrap_or(f64::INFINITY);

                if dist_u_val != f64::INFINITY && dist_u_val + weight_val < dist_v_val {

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
#[must_use]

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

                        capacity[u][v] = try_numeric_value(&t[0]).unwrap_or(0.0);

                        cost[u][v] = try_numeric_value(&t[1]).unwrap_or(0.0);
                    }
                }
            }
        }
    }

    let mut flow = 0.0;

    let mut total_cost = 0.0;

    loop {

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

            break;
        }

        let mut path_flow = f64::INFINITY;

        let mut curr = t;

        while let Some(prev) = parent[curr] {

            path_flow = path_flow.min(capacity[prev][curr]);

            curr = prev;
        }

        flow += path_flow;

        total_cost += path_flow * dist[t];

        let mut v = t;

        while let Some(u) = parent[v] {

            capacity[u][v] -= path_flow;

            capacity[v][u] += path_flow;

            v = u;
        }
    }

    (flow, total_cost)
}

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
#[must_use]

pub fn is_bipartite<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    graph: &Graph<V>,
) -> Option<Vec<i8>> {

    let n = graph.nodes.len();

    let mut colors = vec![-1; n];

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

                            return None;
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
#[must_use]

pub fn bipartite_maximum_matching<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    graph: &Graph<V>,
    partition: &[i8],
) -> Vec<(usize, usize)> {

    let n = graph.nodes.len();

    let mut match_r = vec![None; n]; // Stores the matching for nodes in partition 1
    let mut matching_edges = Vec::new();

    // Identify nodes in partition 0
    let u_nodes: Vec<usize> = (0..n)
        .filter(|&i| partition[i] == 0)
        .collect();

    for &u in &u_nodes {

        let mut visited = vec![false; n];

        bpm_dfs(graph, u, &mut visited, &mut match_r);
    }

    for (v, u_opt) in match_r
        .iter()
        .enumerate()
    {

        if let Some(u) = u_opt {

            if partition[v] == 1 {

                matching_edges.push((*u, v));
            }
        }
    }

    matching_edges
}

fn bpm_dfs<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    graph: &Graph<V>,
    u: usize,
    visited: &mut Vec<bool>,
    match_r: &mut Vec<Option<usize>>,
) -> bool {

    if let Some(neighbors) = graph.adj.get(u) {

        for &(v, _) in neighbors {

            // We only care about edges to unvisited nodes in the other partition
            // (Assumes graph edges are directed or undirected correctly, but for undirected we check visited)
            if !visited[v] {

                visited[v] = true;

                // If v is not matched or its match can find another match
                if match_r[v].is_none() || bpm_dfs(graph, match_r[v].unwrap(), visited, match_r) {

                    match_r[v] = Some(u);

                    return true;
                }
            }
        }
    }

    false
}

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
#[must_use]

pub fn prim_mst<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    graph: &Graph<V>,
    start_node: usize,
) -> Vec<(usize, usize, Expr)> {

    let n = graph.nodes.len();

    let mut mst = Vec::new();

    let mut visited = vec![false; n];

    let mut pq = std::collections::BinaryHeap::new();

    visited[start_node] = true;

    if let Some(neighbors) = graph
        .adj
        .get(start_node)
    {

        for &(v, ref weight) in neighbors {

            let cost = try_numeric_value(weight).unwrap_or(f64::INFINITY);

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

                    let cost = try_numeric_value(next_weight).unwrap_or(f64::INFINITY);

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

    let mut queue: VecDeque<usize> = (0..n)
        .filter(|&i| in_degree[i] == 0)
        .collect();

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
#[must_use]

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

/// Performs a topological sort on a directed acyclic graph (DAG).
///
/// This is a convenience wrapper that checks for cycles and returns None if found.
/// Uses Kahn's algorithm internally.
///
/// # Arguments
/// * `graph` - The directed graph to sort.
///
/// # Returns
/// `Some(Vec<usize>)` containing the topologically sorted node IDs if the graph is a DAG.
/// `None` if the graph contains a cycle.
#[must_use]

pub fn topological_sort<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    graph: &Graph<V>,
) -> Option<Vec<usize>> {

    topological_sort_kahn(graph).ok()
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
#[must_use]

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

    let unmatched_u: Vec<_> = u_nodes
        .difference(&matched_nodes_u)
        .copied()
        .collect();

    let mut visited = HashSet::new();

    let mut queue = VecDeque::from(unmatched_u);

    while let Some(u) = queue.pop_front() {

        if visited.contains(&u) {

            continue;
        }

        visited.insert(u);

        if let Some(neighbors) = graph.adj.get(u) {

            for &(v, _) in neighbors {

                if !matching.contains(&(u, v))
                    && !matching.contains(&(v, u))
                    && !visited.contains(&v)
                {

                    queue.push_back(v);
                }
            }
        }
    }

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
#[must_use]

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

            if pair_u[u].is_none()
                && hopcroft_karp_dfs(graph, u, &mut pair_u, &mut pair_v, &mut dist)
            {

                matching += 1;
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

                    if dist[*next_u] == dist[u] + 1
                        && hopcroft_karp_dfs(graph, *next_u, pair_u, pair_v, dist)
                    {

                        pair_v[v] = Some(u);

                        pair_u[u] = Some(v);

                        return true;
                    }
                } else {

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
) -> Result<Vec<(usize, usize)>, String> {

    let n = graph.nodes.len();

    let mut matching = vec![None; n];

    let mut matches = 0;

    for i in 0..n {

        if matching[i].is_none() {

            let path = find_augmenting_path_with_blossoms(graph, i, &matching)?;

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

    Ok(result)
}

pub(crate) fn find_augmenting_path_with_blossoms<
    V: Eq + std::hash::Hash + Clone + std::fmt::Debug,
>(
    graph: &Graph<V>,
    start_node: usize,
    matching: &[Option<usize>],
) -> Result<Vec<usize>, String> {

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

                    if let Some(w) = matching[v] {

                        parent[v] = Some(u);

                        level[v] = 1;

                        level[w] = 0;

                        queue.push_back(w);
                    } else {

                        parent[v] = Some(u);

                        let mut path = vec![v, u];

                        let mut curr = u;

                        while let Some(p) = parent[curr] {

                            path.push(p);

                            curr = p;
                        }

                        return Ok(path);
                    }
                } else if level[v] == 0 {

                    let base = find_common_ancestor(&origin, &parent, u, v)?;

                    contract_blossom(
                        base,
                        u,
                        v,
                        &mut queue,
                        &mut level,
                        &mut origin,
                        &mut parent,
                        matching,
                    );

                    contract_blossom(
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

    Ok(vec![])
}

pub(crate) fn find_common_ancestor(
    origin: &[usize],
    parent: &[Option<usize>],
    mut u: usize,
    mut v: usize,
) -> Result<usize, String> {

    let mut visited = vec![false; origin.len()];

    loop {

        u = origin[u];

        visited[u] = true;

        if let Some(p) = parent[u] {

            u = p;
        } else {

            break;
        }
    }

    loop {

        v = origin[v];

        if visited[v] {

            return Ok(v);
        }

        if let Some(p) = parent[v] {

            v = p;
        } else {

            return Err("Could not find a common ancestor in blossom algorithm.".to_string());
        }
    }
}

pub(crate) fn contract_blossom(
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

        if let Some(p) = parent[u] {

            u = p;
        } else {

            break;
        }
    }
}

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
#[must_use]

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

        let u_dist = if let Some(d) = distances.get(&u) {

            *d
        } else {

            continue;
        };

        if let Some(neighbors) = graph.adj.get(u) {

            for &(v, _) in neighbors {

                if let std::collections::hash_map::Entry::Vacant(e) = distances.entry(v) {

                    e.insert(u_dist + 1);

                    predecessors.insert(v, Some(u));

                    queue.push_back(v);
                }
            }
        }
    }

    let mut result = HashMap::new();

    for (node, dist) in distances {

        result.insert(
            node,
            (
                dist,
                predecessors
                    .get(&node)
                    .copied()
                    .flatten(),
            ),
        );
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

    pq.push((OrderedFloat(0.0), start_node));

    while let Some((cost, u)) = pq.pop() {

        let cost = -cost.0;

        if cost
            > dist
                .get(&u)
                .and_then(try_numeric_value)
                .unwrap_or(f64::INFINITY)
        {

            continue;
        }

        if let Some(neighbors) = graph.adj.get(u) {

            for &(v, ref weight) in neighbors {

                let new_dist = simplify(&Expr::new_add(Expr::Constant(cost), weight.clone()));

                if try_numeric_value(&new_dist).unwrap_or(f64::INFINITY)
                    < dist
                        .get(&v)
                        .and_then(try_numeric_value)
                        .unwrap_or(f64::INFINITY)
                {

                    dist.insert(v, new_dist.clone());

                    prev.insert(v, Some(u));

                    if let Some(cost) = try_numeric_value(&new_dist) {

                        pq.push((OrderedFloat(-cost), v));
                    }
                }
            }
        }
    }

    let mut result = HashMap::new();

    for (node, d) in dist {

        result.insert(
            node,
            (
                d,
                prev.get(&node)
                    .copied()
                    .flatten(),
            ),
        );
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
#[must_use]

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

                let new_dist = simplify(&Expr::new_add(dist[i][k].clone(), dist[k][j].clone()));

                if try_numeric_value(&dist[i][j]).unwrap_or(f64::INFINITY)
                    > try_numeric_value(&new_dist).unwrap_or(f64::INFINITY)
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

        let mut numerical_eigenvalues = Vec::new();

        for val_expr in eig_vec
            .iter()
            .flatten()
        {

            if let Some(val) = try_numeric_value(val_expr) {

                numerical_eigenvalues.push(val);
            } else {

                return Err("Eigenvalues are not all numerical, cannot sort.".to_string());
            }
        }

        numerical_eigenvalues.sort_by(|a, b| {
            a.partial_cmp(b)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        Ok(Expr::Constant(numerical_eigenvalues[1]))
    } else {

        Err("Eigenvalue computation did not return a vector.".to_string())
    }
}
