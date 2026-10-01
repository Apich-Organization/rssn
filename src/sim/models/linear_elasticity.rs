use std::fs::File;
use std::io::Write;
use std::path::Path;

use ndarray::Array1;
use ndarray::Array2;
use ndarray::array;
use rayon::prelude::*;
use serde::Deserialize;
use serde::Serialize;
use sprs_rssn::CsMat;

use crate::kernels::sparse::csr_from_triplets;
use crate::kernels::sparse::solve_conjugate_gradient;

/// Defines the node points of the mesh.
pub type Nodes = Vec<(f64, f64)>;

/// Defines the elements by indexing into the nodes vector.
pub type Elements = Vec<[usize; 4]>;

/// Parameters for the linear elasticity simulation.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ElasticityParameters {
    /// The nodes of the mesh.
    pub nodes: Nodes,
    /// The elements of the mesh.
    pub elements: Elements,
    /// The Young's modulus of the material.
    pub youngs_modulus: f64,
    /// The Poisson's ratio of the material.
    pub poissons_ratio: f64,
    /// The indices of the nodes that are fixed.
    pub fixed_nodes: Vec<usize>,
    /// The loads applied to the nodes.
    pub loads: Vec<(usize, f64, f64)>,
}

/// Calculates the element stiffness matrix for a 2D quadrilateral element (plane stress).
///
/// This is the standard bilinear isoparametric Q4 element of unit thickness:
/// `K = ∫ Bᵀ C B det(J) dξ dη` evaluated with 2x2 Gauss quadrature
/// (Zienkiewicz & Taylor, "The Finite Element Method", Vol. 1, ch. 6;
/// Cook et al., "Concepts and Applications of Finite Element Analysis", ch. 6).
/// The corners `p1..p4` must be given counter-clockwise (bottom-left,
/// bottom-right, top-right, top-left for an axis-aligned element) and the
/// degrees of freedom are ordered `[u1, v1, u2, v2, u3, v3, u4, v4]`.
/// Multiply by the thickness for a thicker plate.
#[must_use]
pub fn element_stiffness_matrix(
    p1: (f64, f64),
    p2: (f64, f64),
    p3: (f64, f64),
    p4: (f64, f64),
    e: f64,
    nu: f64,
) -> Array2<f64> {
    let pts = [p1, p2, p3, p4];

    // Natural coordinates of the corners.
    let xi_n = [-1.0, 1.0, 1.0, -1.0];

    let eta_n = [-1.0, -1.0, 1.0, 1.0];

    let c_mat = (e / nu.mul_add(-nu, 1.0))
        * array![[1.0, nu, 0.0], [nu, 1.0, 0.0], [0.0, 0.0, (1.0 - nu) / 2.0]];

    let g = 1.0 / 3.0_f64.sqrt();

    let mut k = Array2::<f64>::zeros((8, 8));

    for &xi in &[-g, g] {
        for &eta in &[-g, g] {
            // Shape-function derivatives in natural coordinates.
            let mut dn_dxi = [0.0; 4];

            let mut dn_deta = [0.0; 4];

            for n in 0..4 {
                dn_dxi[n] = 0.25 * xi_n[n] * (1.0 + eta_n[n] * eta);

                dn_deta[n] = 0.25 * eta_n[n] * (1.0 + xi_n[n] * xi);
            }

            // Jacobian J = [[dx/dxi, dy/dxi], [dx/deta, dy/deta]].
            let mut j = [[0.0; 2]; 2];

            for n in 0..4 {
                j[0][0] += dn_dxi[n] * pts[n].0;

                j[0][1] += dn_dxi[n] * pts[n].1;

                j[1][0] += dn_deta[n] * pts[n].0;

                j[1][1] += dn_deta[n] * pts[n].1;
            }

            let det = j[0][0].mul_add(j[1][1], -(j[0][1] * j[1][0]));

            let inv = [
                [j[1][1] / det, -j[0][1] / det],
                [-j[1][0] / det, j[0][0] / det],
            ];

            // Strain-displacement matrix B (3 x 8): [ex, ey, gxy].
            let mut b_mat = Array2::<f64>::zeros((3, 8));

            for n in 0..4 {
                let dn_dx = inv[0][0].mul_add(dn_dxi[n], inv[0][1] * dn_deta[n]);

                let dn_dy = inv[1][0].mul_add(dn_dxi[n], inv[1][1] * dn_deta[n]);

                b_mat[[0, 2 * n]] = dn_dx;

                b_mat[[1, 2 * n + 1]] = dn_dy;

                b_mat[[2, 2 * n]] = dn_dy;

                b_mat[[2, 2 * n + 1]] = dn_dx;
            }

            // Gauss weights are 1 for the 2-point rule.
            k = k + b_mat.t().dot(&c_mat.dot(&b_mat)) * det.abs();
        }
    }

    k
}

/// Runs a 2D linear elasticity simulation using the Finite Element Method.
///
/// This function assembles the global stiffness matrix and force vector based on
/// the provided nodes, elements, material properties, boundary conditions, and loads.
/// It then solves the resulting linear system to find the nodal displacements.
///
/// # Arguments
/// * `params` - An `ElasticityParameters` struct containing all simulation inputs.
///
/// # Returns
/// A `Result` containing a `Vec<f64>` of nodal displacements (u, v for each node).
///
/// # Errors
///
/// This function will return an error if the Conjugate Gradient solver fails to converge
/// or if the global stiffness matrix is ill-conditioned.
pub fn run_elasticity_simulation(params: &ElasticityParameters) -> Result<Vec<f64>, String> {
    let n_nodes = params.nodes.len();

    let n_dofs = n_nodes * 2;

    // Parallel element stiffness matrix assembly
    let triplets: Vec<(usize, usize, f64)> = params
        .elements
        .par_iter()
        .flat_map(|element| {
            let p1 = params.nodes[element[0]];

            let p2 = params.nodes[element[1]];

            let p3 = params.nodes[element[2]];

            let p4 = params.nodes[element[3]];

            let k_element = element_stiffness_matrix(
                p1,
                p2,
                p3,
                p4,
                params.youngs_modulus,
                params.poissons_ratio,
            );

            let dof_indices = [
                element[0] * 2,
                element[0] * 2 + 1,
                element[1] * 2,
                element[1] * 2 + 1,
                element[2] * 2,
                element[2] * 2 + 1,
                element[3] * 2,
                element[3] * 2 + 1,
            ];

            let mut element_triplets = Vec::with_capacity(64);

            for r in 0..8 {
                for c in 0..8 {
                    element_triplets.push((dof_indices[r], dof_indices[c], k_element[[r, c]]));
                }
            }

            element_triplets
        })
        .collect();

    let mut f_global = Array1::<f64>::zeros(n_dofs);

    for &(node_idx, fx, fy) in &params.loads {
        f_global[node_idx * 2] += fx;

        f_global[node_idx * 2 + 1] += fy;
    }

    // Handle boundary conditions: zero out rows and columns of fixed DOFs
    let fixed_dofs: std::collections::HashSet<usize> = params
        .fixed_nodes
        .iter()
        .flat_map(|&node_idx| vec![node_idx * 2, node_idx * 2 + 1])
        .collect();

    let mut filtered_triplets: Vec<(usize, usize, f64)> = triplets
        .into_par_iter()
        .filter(|(r, c, _)| !fixed_dofs.contains(r) && !fixed_dofs.contains(c))
        .collect();

    for &node_idx in &params.fixed_nodes {
        let dof1 = node_idx * 2;

        let dof2 = node_idx * 2 + 1;

        filtered_triplets.push((dof1, dof1, 1.0));

        filtered_triplets.push((dof2, dof2, 1.0));

        f_global[dof1] = 0.0;

        f_global[dof2] = 0.0;
    }

    let k_global: CsMat<f64> = csr_from_triplets(n_dofs, n_dofs, &filtered_triplets);

    let displacements = solve_conjugate_gradient(&k_global, &f_global, None, 5000, 1e-9)?;

    Ok(displacements.to_vec())
}

/// An example scenario for a cantilever beam under a point load.
///
/// This function sets up a mesh for a 2D cantilever beam, defines fixed boundary
/// conditions at one end and applies a point load at the free end. It then runs
/// the elasticity simulation and saves the original and deformed node positions
/// to CSV files for visualization.
///
/// # Errors
///
/// This function will return an error if the underlying `run_elasticity_simulation`
/// fails or if it cannot create or write to the output CSV files.
/// `beam_original.csv` and `beam_deformed.csv` are written into `output_dir`,
/// which is created if it does not exist.
pub fn simulate_cantilever_beam_scenario(output_dir: &Path) -> Result<(), String> {
    println!(
        "Running 2D Cantilever Beam \
         simulation..."
    );

    let beam_length = 10.0;

    let beam_height = 2.0;

    let nx = 20;

    let ny = 4;

    let mut nodes: Nodes = Vec::new();

    for j in 0..=ny {
        for i in 0..=nx {
            nodes.push((
                i as f64 * beam_length / nx as f64,
                j as f64 * beam_height / ny as f64,
            ));
        }
    }

    let mut elements: Elements = Vec::new();

    for j in 0..ny {
        for i in 0..nx {
            let n1 = j * (nx + 1) + i;

            let n2 = j * (nx + 1) + i + 1;

            let n3 = (j + 1) * (nx + 1) + i + 1;

            let n4 = (j + 1) * (nx + 1) + i;

            elements.push([n1, n2, n3, n4]);
        }
    }

    let fixed_nodes: Vec<usize> = (0..=ny).map(|j| j * (nx + 1)).collect();

    let loads = vec![((ny / 2) * (nx + 1) + nx, 0.0, -1e3)];

    let params = ElasticityParameters {
        nodes: nodes.clone(),
        elements,
        youngs_modulus: 1e7,
        poissons_ratio: 0.3,
        fixed_nodes,
        loads,
    };

    let d = run_elasticity_simulation(&params)?;

    println!(
        "Simulation finished. Saving \
         results..."
    );

    let mut new_nodes = nodes.clone();

    for i in 0..nodes.len() {
        new_nodes[i].0 += d[i * 2];

        new_nodes[i].1 += d[i * 2 + 1];
    }

    std::fs::create_dir_all(output_dir).map_err(|e| e.to_string())?;

    let mut orig_file =
        File::create(output_dir.join("beam_original.csv")).map_err(|e| e.to_string())?;

    let mut def_file =
        File::create(output_dir.join("beam_deformed.csv")).map_err(|e| e.to_string())?;

    writeln!(orig_file, "x,y").map_err(|e| e.to_string())?;

    writeln!(def_file, "x,y").map_err(|e| e.to_string())?;

    for n in &nodes {
        writeln!(orig_file, "{},{}", n.0, n.1).map_err(|e| e.to_string())?;
    }

    for n in &new_nodes {
        writeln!(def_file, "{},{}", n.0, n.1).map_err(|e| e.to_string())?;
    }

    println!(
        "Original and deformed node \
         positions saved to .csv \
         files."
    );

    Ok(())
}
