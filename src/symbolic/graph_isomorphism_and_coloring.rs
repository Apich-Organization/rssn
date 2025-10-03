//! # Graph Isomorphism and Coloring Algorithms
//!
//! This module provides algorithms for two fundamental problems in graph theory:
//! graph isomorphism testing and graph coloring. It includes a heuristic approach
//! for isomorphism using the Weisfeiler-Lehman test and greedy algorithms for
//! vertex coloring, as well as an exact (NP-hard) chromatic number solver.

use crate::prelude::Expr;
use crate::symbolic::graph::Graph;
use std::collections::HashMap;

// =====================================================================================
// region: Graph Isomorphism
// =====================================================================================

/// Checks if two graphs are potentially isomorphic using the Weisfeiler-Lehman test (Color Refinement).
///
/// This is a powerful heuristic. It returns `true` if the color histograms of the two graphs match
/// after a sufficient number of refinement iterations, and `false` otherwise.
///
/// **Note**: Passing this test does not guarantee isomorphism for all graph classes,
/// as there exist non-isomorphic graphs that cannot be distinguished by the WL test.
///
/// # Arguments
/// * `g1` - The first graph.
/// * `g2` - The second graph.
///
/// # Returns
/// `true` if the graphs are indistinguishable by the WL test, `false` otherwise.
pub fn are_isomorphic_heuristic<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    g1: &Graph<V>,
    g2: &Graph<V>,
) -> bool {
    if g1.nodes.len() != g2.nodes.len() || g1.get_edges().len() != g2.get_edges().len() {
        return false;
    }

    let h1 = wl_test(g1);
    let h2 = wl_test(g2);

    h1 == h2
}

// Performs the 1-dimensional Weisfeiler-Lehman test.
pub(crate) fn wl_test<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    graph: &Graph<V>,
) -> HashMap<String, usize> {
    let n = graph.nodes.len();
    let mut colors: Vec<String> = (0..n).map(|i| graph.in_degree(i).to_string()).collect();

    for _ in 0..n {
        // Iterate a number of times related to the number of vertices
        let mut next_colors = Vec::with_capacity(n);
        for i in 0..n {
            let mut neighbor_colors = Vec::new();
            if let Some(neighbors) = graph.adj.get(i) {
                for &(v, _) in neighbors {
                    neighbor_colors.push(colors[v].clone());
                }
            }
            neighbor_colors.sort();
            let new_color_signature = format!("{}-{}", colors[i], neighbor_colors.join(","));
            // In a real implementation, we would hash this signature to get a new color label.
            // For simplicity, we use the signature string itself.
            next_colors.push(new_color_signature);
        }
        colors = next_colors;
    }

    let mut histogram = HashMap::new();
    for color in colors {
        *histogram.entry(color).or_insert(0) += 1;
    }
    histogram
}

// =====================================================================================
// region: Graph Coloring
// =====================================================================================

/// Finds a valid vertex coloring using a greedy heuristic (Welsh-Powell algorithm).
///
/// This algorithm sorts the vertices by degree in descending order and then assigns
/// the smallest available color to each vertex. It is a heuristic and does not
/// guarantee an optimal coloring (i.e., finding the chromatic number).
///
/// # Arguments
/// * `graph` - The graph to color.
///
/// # Returns
/// A `HashMap<usize, usize>` where keys are node IDs and values are their assigned colors.
pub fn greedy_coloring<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    graph: &Graph<V>,
) -> HashMap<usize, usize> {
    let mut nodes: Vec<usize> = (0..graph.nodes.len()).collect();
    // Sort nodes by degree in descending order
    nodes.sort_by(|a, b| graph.out_degree(*b).cmp(&graph.out_degree(*a)));

    let mut colors = HashMap::new();
    let mut color_counter = 0;

    for &node_id in &nodes {
        if !colors.contains_key(&node_id) {
            colors.insert(node_id, color_counter);
            for &other_node_id in &nodes {
                if !colors.contains_key(&other_node_id) {
                    // Check if other_node_id is not adjacent to any node already in the current color class
                    let is_safe = graph.adj.get(other_node_id).map_or(true, |neighbors: &_| {
                        neighbors
                            .iter()
                            .all(|(n, _): &(usize, Expr)| colors.get(n) != Some(&color_counter))
                    });
                    if is_safe {
                        colors.insert(other_node_id, color_counter);
                    }
                }
            }
            color_counter += 1;
        }
    }
    colors
}

/// Finds the chromatic number of a graph using exhaustive backtracking.
///
/// The chromatic number `χ(G)` is the minimum number of colors needed to color the vertices
/// of a graph such that no two adjacent vertices share the same color.
/// This is an NP-hard problem, so this exact algorithm will be very slow for large graphs.
///
/// # Arguments
/// * `graph` - The graph to analyze.
///
/// # Returns
/// The chromatic number as a `usize`.
pub fn chromatic_number_exact<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    graph: &Graph<V>,
) -> usize {
    let n = graph.nodes.len();
    if n == 0 {
        return 0;
    }
    for k in 1..=n {
        let mut colors = vec![0; n];
        if can_color_with_k(graph, k, &mut colors, 0) {
            return k;
        }
    }
    n
}

pub(crate) fn can_color_with_k<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    graph: &Graph<V>,
    k: usize,
    colors: &mut [usize],
    node_idx: usize,
) -> bool {
    if node_idx == graph.nodes.len() {
        return true; // All nodes have been colored
    }

    for c in 1..=k {
        if is_safe_to_color(graph, node_idx, c, colors) {
            colors[node_idx] = c;
            if can_color_with_k(graph, k, colors, node_idx + 1) {
                return true;
            }
            colors[node_idx] = 0; // Backtrack
        }
    }
    false
}

pub(crate) fn is_safe_to_color<V: Eq + std::hash::Hash + Clone + std::fmt::Debug>(
    graph: &Graph<V>,
    u: usize,
    color: usize,
    colors: &[usize],
) -> bool {
    if let Some(neighbors) = graph.adj.get(u) {
        for &(v, _) in neighbors {
            if colors[v] == color {
                return false;
            }
        }
    }
    true
}
