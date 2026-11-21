//! # Finite Element Method (FEM) Solvers
//!
//! This module implements Finite Element Analysis (FEM) tools for solving partial differential equations,
//! specifically focusing on the Poisson equation in one, two, and three dimensions.
//!
//! # Overview
//!
//! The Finite Element Method is a powerful numerical technique for solving problems of engineering and mathematical physics.
//! This module provides a set of solvers that discretize the domain into simpler elements (lines, quads, hexes) to approximate
//! the solution of the Poisson equation with Dirichlet boundary conditions.
//!
//! Key components include:
//! - **Dimensional Solvers**: Specialized solvers for 1D (`solve_poisson_1d`), 2D (`solve_poisson_2d`), and 3D (`solve_poisson_3d`) problems.
//! - **Parallel Assembly**: Independent elements compute their load contributions on the active CPU pool for large grids. Contributions are added in grid order.
//! - **Numerical Integration**: Uses Gaussian quadrature for quadrilateral and hexahedral loads.
//! - **Linear Solves**: The one-dimensional path uses sparse conjugate gradients.
//!   Uniform quadrilateral and hexahedral meshes share a prepared sine-basis solver.
//! - **Simulation Scenarios**: Ready-to-run examples demonstrating the application of the solvers to standard test problems.
//!
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/fem_bridge_3d_disp.png)
//!
//! # Coverage audit
//!
//! | Equation | 1D | 2D | 3D | Elements | Boundary conditions |
//! |---|---|---|---|---|---|
//! | Poisson | yes | yes (quads) | yes (hexahedra) | linear / bilinear / trilinear | zero Dirichlet |
//!
//! Unstructured simplex meshes (intervals, triangles, tetrahedra) with P1/P2
//! elements, marker-based Dirichlet/Neumann/Robin conditions, time-dependent data
//! and the θ-method heat equation are in [`crate::sim::fem_simplex`]; elasticity is
//! in [`crate::sim::elasticity_fem`].

use ndarray::Array1;
use rayon::prelude::*;
use serde::Deserialize;
use serde::Serialize;

use crate::kernels::sparse::csr_from_triplets;
use crate::kernels::sparse::solve_conjugate_gradient;

type FemTriplet = (usize, usize, f64);

type ElementData1D = (Vec<FemTriplet>, [f64; 2]);

mod cartesian;

pub use cartesian::PreparedCartesianPoisson;

/// A struct for Gaussian quadrature.
#[allow(dead_code)]
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct GaussQuadrature {
    points: Vec<f64>,
    weights: Vec<f64>,
}

impl GaussQuadrature {
    pub(crate) fn new() -> Self {
        let points = vec![-1.0 / 3.0_f64.sqrt(), 1.0 / 3.0_f64.sqrt()];

        let weights = vec![1.0, 1.0];

        Self { points, weights }
    }
}

/// Solves the 1D Poisson equation: -d^2u/dx^2 = f(x)
/// with Dirichlet boundary conditions u(0) = 0, u(L) = 0.
///
/// # Arguments
/// * `n_elements` - Number of linear elements.
/// * `domain_length` - The length L of the domain.
/// * `force_fn` - The forcing function f(x).
///
/// # Returns
/// A vector representing the solution u at each node.
///
/// # Errors
///
/// This function will return an error if `n_elements` is 0, or if the underlying
/// linear system solver fails to converge or encounters numerical instability.
pub fn solve_poisson_1d<F>(
    n_elements: usize,
    domain_length: f64,
    force_fn: F,
) -> Result<Vec<f64>, String>
where
    F: Fn(f64) -> f64 + Send + Sync,
{
    let n_nodes = n_elements + 1;

    let h = domain_length / n_elements as f64;

    let force_fn = &force_fn;

    // Parallel element assembly
    let element_data: Vec<ElementData1D> = (0..n_elements)
        .into_par_iter()
        .map(move |i| {
            let x1 = i as f64 * h;

            let x2 = (i + 1) as f64 * h;

            let k_local = [[1.0 / h, -1.0 / h], [-1.0 / h, 1.0 / h]];

            let f_local = [h / 2.0 * force_fn(x1), h / 2.0 * force_fn(x2)];

            let nodes = [i, i + 1];

            let mut local_triplets = Vec::with_capacity(4);

            for r in 0..2 {
                for c in 0..2 {
                    local_triplets.push((nodes[r], nodes[c], k_local[r][c]));
                }
            }

            (local_triplets, f_local)
        })
        .collect();

    let mut triplets = Vec::with_capacity(n_elements * 4);

    let mut f = vec![0.0; n_nodes];

    for (i, (local_triplets, f_vals)) in element_data.into_iter().enumerate() {
        triplets.extend(local_triplets);

        f[i] += f_vals[0];

        f[i + 1] += f_vals[1];
    }

    let last_node = n_nodes - 1;

    triplets.retain(|(r, _, _)| *r != 0 && *r != last_node);

    triplets.push((0, 0, 1.0));

    triplets.push((last_node, last_node, 1.0));

    f[0] = 0.0;

    f[last_node] = 0.0;

    let k_sparse = csr_from_triplets(n_nodes, n_nodes, &triplets);

    let f_array = Array1::from(f);

    let u_array = solve_conjugate_gradient(&k_sparse, &f_array, None, 1000, 1e-9)?;

    Ok(u_array.to_vec())
}

/// Example scenario for the 1D FEM Poisson solver.
///
/// # Errors
///
/// This function will return an error if the underlying `solve_poisson_1d` function
/// encounters an error.
pub fn simulate_1d_poisson_scenario() -> Result<Vec<f64>, String> {
    const N_ELEMENTS: usize = 50;

    const L: f64 = 1.0;

    let force = |_x: f64| 2.0;

    solve_poisson_1d(N_ELEMENTS, L, force)
}

/// Solves the 2D Poisson equation on a unit square with zero Dirichlet boundaries.
///
/// Uses the same tensor-product Q1 stiffness assembled from bilinear or
/// trilinear elements, diagonalized by the discrete sine basis.
///
/// # Errors
/// Reports invalid extents, non-finite source values, cancellation or an
/// unrepresentable result.
pub fn solve_poisson_2d<F>(
    n_elements_x: usize,
    n_elements_y: usize,
    force_fn: F,
) -> Result<Vec<f64>, String>
where
    F: Fn(f64, f64) -> f64 + Send + Sync,
{
    let plan = PreparedCartesianPoisson::new(&[n_elements_x, n_elements_y])?;
    plan.solve_legacy(|[x, y, _]| force_fn(x, y))
}

/// Example scenario for the 2D FEM Poisson solver.
///
/// # Errors
///
/// This function will return an error if the underlying `solve_poisson_2d` function
/// encounters an error.
#[allow(clippy::unnecessary_cast)]
pub fn simulate_2d_poisson_scenario() -> Result<Vec<f64>, String> {
    const N_ELEMENTS: usize = 20;

    let force = |x, y| {
        2.0 * std::f64::consts::PI.powi(2)
            * (std::f64::consts::PI * (x as f64)).sin()
            * (std::f64::consts::PI * (y as f64)).sin()
    };

    solve_poisson_2d(N_ELEMENTS, N_ELEMENTS, force)
}

/// Solves the 3D Poisson equation on a unit cube with zero Dirichlet boundaries.
///
/// Uses the same tensor-product Q1 stiffness assembled from bilinear or
/// trilinear elements, diagonalized by the discrete sine basis.
///
/// # Errors
/// Reports invalid extents, non-finite source values, cancellation or an
/// unrepresentable result.
pub fn solve_poisson_3d<F>(
    n_elements: usize,
    force_fn: F,
) -> Result<Vec<f64>, String>
where
    F: Fn(f64, f64, f64) -> f64 + Send + Sync,
{
    let plan = PreparedCartesianPoisson::new(&[n_elements; 3])?;
    plan.solve_legacy(|[x, y, z]| force_fn(x, y, z))
}

/// Example scenario for the 3D FEM Poisson solver.
///
/// # Errors
///
/// This function will return an error if the underlying `solve_poisson_3d` function
/// encounters an error.
#[allow(clippy::unnecessary_cast)]
pub fn simulate_3d_poisson_scenario() -> Result<Vec<f64>, String> {
    const N_ELEMENTS: usize = 10;

    let force = |x, y, z| {
        3.0 * std::f64::consts::PI.powi(2)
            * (std::f64::consts::PI * (x as f64)).sin()
            * (std::f64::consts::PI * (y as f64)).sin()
            * (std::f64::consts::PI * (z as f64)).sin()
    };

    solve_poisson_3d(N_ELEMENTS, force)
}

/// Poisson/reaction–diffusion `-Δu + c u = f` on interval, triangle or
/// tetrahedron meshes with P1/P2 elements and marker boundary conditions.
///
/// # Errors
/// For unsupported boundary conditions or a singular system.
pub fn solve_poisson_simplex(
    space: &crate::sim::fem_simplex::FemSpace,
    c: f64,
    source: &(dyn Fn(f64, [f64; 3]) -> f64 + Sync),
    bcs: &crate::sim::boundary::MarkerBcs<crate::sim::boundary::BoundaryCondition>,
) -> Result<Vec<f64>, String> {
    let prob = crate::sim::fem_simplex::FemProblem {
        k: 1.0,
        c,
        f: source,
        bcs,
    };
    crate::sim::fem_simplex::solve_static(
        space,
        &prob,
        0.0,
        crate::sim::structured::LinearSolver::Auto,
    )
}
