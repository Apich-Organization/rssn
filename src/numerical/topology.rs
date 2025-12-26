//! # Numerical Computational Topology
//!
//! This module provides numerical tools for computational topology.
//! It includes algorithms for finding connected components in graphs,
//! constructing Vietoris-Rips simplicial complexes, and computing persistent homology.

use crate::numerical::graph::Graph;
use crate::symbolic::topology::{
    ChainComplex,
    SimplicialComplex,
};
use serde::{
    Deserialize,
    Serialize,
};
use std::collections::VecDeque;

/// Finds the connected components of a graph using Breadth-First Search (BFS).
///
/// # Arguments
/// * `graph` - The graph to analyze.
///
/// # Returns
/// A vector of vectors, where each inner vector contains the node indices of a connected component.
#[must_use]

pub fn find_connected_components(
    graph: &Graph
) -> Vec<Vec<usize>> {

    let num_nodes = graph.num_nodes();

    let mut visited =
        vec![false; num_nodes];

    let mut components = Vec::new();

    for i in 0..num_nodes {

        if !visited[i] {

            let mut component =
                Vec::new();

            let mut queue =
                VecDeque::new();

            visited[i] = true;

            queue.push_back(i);

            while let Some(u) =
                queue.pop_front()
            {

                component.push(u);

                for &(v, _) in
                    graph.adj(u)
                {

                    if !visited[v] {

                        visited[v] =
                            true;

                        queue
                            .push_back(
                                v,
                            );
                    }
                }
            }

            components.push(component);
        }
    }

    components
}

/// Represents a simplex in a simplicial complex.

pub type Simplex = Vec<usize>;

/// Represents a persistence interval (birth, death).
#[derive(
    Debug,
    Clone,
    Copy,
    Serialize,
    Deserialize,
    PartialEq,
)]

pub struct PersistenceInterval {
    pub birth: f64,
    pub death: f64,
}

/// Represents a persistence diagram for a specific dimension.
#[derive(
    Debug,
    Clone,
    Serialize,
    Deserialize,
    PartialEq,
)]

pub struct PersistenceDiagram {
    pub dimension: usize,
    pub intervals:
        Vec<PersistenceInterval>,
}

/// Constructs a Vietoris-Rips simplicial complex from a set of points for a given radius.
///
/// # Arguments
/// * `points` - A slice of points, where each point is a slice of `f64`.
/// * `epsilon` - The distance threshold (radius).
/// * `max_dim` - The maximum dimension of simplices to compute.
#[must_use]

pub fn vietoris_rips_complex(
    points: &[&[f64]],
    epsilon: f64,
    max_dim: usize,
) -> Vec<Simplex> {

    let n_points = points.len();

    let mut simplices = Vec::new();

    for i in 0..n_points {

        simplices.push(vec![i]);
    }

    if max_dim == 0 {

        return simplices;
    }

    let mut current_simplices =
        Vec::new();

    for i in 0..n_points {

        for j in (i + 1)..n_points {

            if euclidean_distance(
                points[i], points[j],
            ) <= epsilon
            {

                current_simplices
                    .push(vec![i, j]);
            }
        }
    }

    simplices.extend(
        current_simplices.clone(),
    );

    for _dim in 2..=max_dim {

        let mut next_simplices =
            Vec::new();

        for simplex in
            &current_simplices
        {

            let last_v = simplex
                [simplex.len() - 1];

            for i in
                (last_v + 1)..n_points
            {

                let mut all_ok = true;

                for &v in simplex {

                    if euclidean_distance(points[v], points[i]) > epsilon {

                        all_ok = false;

                        break;
                    }
                }

                if all_ok {

                    let mut new_simplex =
                        simplex.clone();

                    new_simplex.push(i);

                    next_simplices
                        .push(
                            new_simplex,
                        );
                }
            }
        }

        if next_simplices.is_empty() {

            break;
        }

        simplices.extend(
            next_simplices.clone(),
        );

        current_simplices =
            next_simplices;
    }

    simplices
}

/// Computes the Betti numbers for a point cloud at a given radius.
///
/// # Arguments
/// * `points` - A slice of points.
/// * `epsilon` - The distance threshold.
/// * `max_dim` - Maximum dimension to compute Betti numbers for.

pub fn betti_numbers_at_radius(
    points: &[&[f64]],
    epsilon: f64,
    max_dim: usize,
) -> Vec<usize> {

    let simplices =
        vietoris_rips_complex(
            points, epsilon, max_dim,
        );

    let mut complex =
        SimplicialComplex::new();

    for s in simplices {

        complex.add_simplex(&s);
    }

    let chain_complex =
        ChainComplex::new(complex);

    let mut betti = Vec::new();

    for k in 0..=max_dim {

        betti.push(
            chain_complex
                .compute_homology_betti_number(k)
                .unwrap_or(0),
        );
    }

    betti
}

/// Computes the persistent homology (persistence diagram) for a point cloud.
/// This is a simplified version using a fixed number of steps.

pub fn compute_persistence(
    points: &[Vec<f64>],
    max_epsilon: f64,
    steps: usize,
    max_dim: usize,
) -> Vec<PersistenceDiagram> {

    let mut diagrams = vec![
        PersistenceDiagram {
            dimension: 0,
            intervals: Vec::new()
        };
        max_dim + 1
    ];

    for d in 0..=max_dim {

        diagrams[d].dimension = d;
    }

    // This is a naive implementation: we track Betti number changes.
    // Real persistent homology requires tracking individual cycles, but for FFI start
    // we provide a way to see topological features appearing and disappearing.

    let mut prev_betti =
        vec![0; max_dim + 1];

    let mut open_intervals: Vec<
        Vec<f64>,
    > = vec![Vec::new(); max_dim + 1];

    for step in 0..=steps {

        let epsilon = max_epsilon
            * (step as f64
                / steps as f64);

        let current_betti =
            betti_numbers_at_radius(
                &points
                    .iter()
                    .map(|v| {
                        v.as_slice()
                    })
                    .collect::<Vec<_>>(
                    ),
                epsilon,
                max_dim,
            );

        for d in 0..=max_dim {

            let cb = current_betti[d];

            let pb = prev_betti[d];

            if cb > pb {

                // New features born
                for _ in 0..(cb - pb) {

                    open_intervals[d]
                        .push(epsilon);
                }
            } else if cb < pb {

                // Features died
                for _ in 0..(pb - cb) {

                    if let Some(birth) =
                        open_intervals
                            [d]
                            .pop()
                    {

                        diagrams[d]
                            .intervals
                            .push(
                                PersistenceInterval {
                                    birth,
                                    death: epsilon,
                                },
                            );
                    }
                }
            }
        }

        prev_betti = current_betti;
    }

    // Close remaining intervals
    for d in 0..=max_dim {

        while let Some(birth) =
            open_intervals[d].pop()
        {

            diagrams[d]
                .intervals
                .push(
                    PersistenceInterval {
                        birth,
                        death: max_epsilon,
                    },
                );
        }
    }

    diagrams
}

/// Computes the Euclidean distance between two points.

pub fn euclidean_distance(
    p1: &[f64],
    p2: &[f64],
) -> f64 {

    p1.iter()
        .zip(p2.iter())
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>()
        .sqrt()
}
