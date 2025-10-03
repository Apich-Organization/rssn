//! # Numerical Computational Topology
//!
//! This module provides numerical tools for computational topology.
//! It includes algorithms for finding connected components in graphs and constructing
//! Vietoris-Rips simplicial complexes from point clouds.

use crate::numerical::graph::Graph;
use std::collections::VecDeque;

/// Finds the connected components of a graph using Breadth-First Search (BFS).
///
/// # Arguments
/// * `graph` - The graph to analyze.
///
/// # Returns
/// A vector of vectors, where each inner vector contains the node indices of a connected component.
pub fn find_connected_components(graph: &Graph) -> Vec<Vec<usize>> {
    let num_nodes = graph.num_nodes();
    let mut visited = vec![false; num_nodes];
    let mut components = Vec::new();

    for i in 0..num_nodes {
        if !visited[i] {
            let mut component = Vec::new();
            let mut queue = VecDeque::new();

            visited[i] = true;
            queue.push_back(i);

            while let Some(u) = queue.pop_front() {
                component.push(u);
                for &(v, _) in graph.adj(u) {
                    if !visited[v] {
                        visited[v] = true;
                        queue.push_back(v);
                    }
                }
            }
            components.push(component);
        }
    }
    components
}

/// Represents a simplex in a simplicial complex (e.g., \[0\], \[0,1\], \[0,1,2\]).
pub type Simplex = Vec<usize>;

/// Constructs a Vietoris-Rips simplicial complex from a set of points for a given radius.
///
/// A Vietoris-Rips complex is a type of simplicial complex constructed from a set of points.
/// A simplex is included if all its vertices are within `epsilon` distance of each other.
///
/// # Arguments
/// * `points` - A slice of points, where each point is a slice of `f64`.
/// * `epsilon` - The distance threshold (radius).
/// * `max_dim` - The maximum dimension of simplices to compute (e.g., 2 for triangles).
///
/// # Returns
/// A vector of all simplices in the complex up to the specified dimension.
pub fn vietoris_rips_complex(points: &[&[f64]], epsilon: f64, max_dim: usize) -> Vec<Simplex> {
    let n_points = points.len();
    let mut simplices = Vec::new();

    // 0-simplices (vertices)
    for i in 0..n_points {
        simplices.push(vec![i]);
    }
    if max_dim == 0 {
        return simplices;
    }

    // 1-simplices (edges)
    let mut edges = Vec::new();
    for i in 0..n_points {
        for j in (i + 1)..n_points {
            let dist = euclidean_distance(points[i], points[j]);
            if dist < epsilon {
                edges.push(vec![i, j]);
            }
        }
    }
    simplices.extend(edges.clone());
    if max_dim == 1 {
        return simplices;
    }

    // 2-simplices (triangles) and higher, built from lower-dimensional simplices
    let mut current_simplices = edges;
    for _dim in 2..=max_dim {
        let mut next_simplices = Vec::new();
        for simplex in &current_simplices {
            for i in (simplex.last().unwrap() + 1)..n_points {
                // Check if the new point `i` is connected to all vertices in the current simplex
                let is_connected_to_all = simplex
                    .iter()
                    .all(|&v| euclidean_distance(points[v], points[i]) < epsilon);

                if is_connected_to_all {
                    let mut new_simplex = simplex.clone();
                    new_simplex.push(i);
                    next_simplices.push(new_simplex);
                }
            }
        }
        if next_simplices.is_empty() {
            break;
        }
        simplices.extend(next_simplices.clone());
        current_simplices = next_simplices;
    }

    simplices
}

pub(crate) fn euclidean_distance(p1: &[f64], p2: &[f64]) -> f64 {
    p1.iter()
        .zip(p2.iter())
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>()
        .sqrt()
}
