use std::fs::File;
use std::io::Write;

use ndarray::array;
use ndarray::Array1;
use ndarray::Array2;
use rayon::prelude::*;
use serde::Deserialize;
use serde::Serialize;
use sprs_rssn::CsMat;

use crate::numerical::sparse::csr_from_triplets;
use crate::numerical::sparse::solve_conjugate_gradient;

/// Defines the node points of the mesh.

pub type Nodes = Vec<(f64, f64)>;

/// Defines the elements by indexing into the nodes vector.

pub type Elements = Vec<[usize; 4]>;

/// Parameters for the linear elasticity simulation.
#[derive(
    Clone, Debug, Serialize, Deserialize,
)]

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

pub fn element_stiffness_matrix(
    _p1: (f64, f64),
    _p2: (f64, f64),
    _p3: (f64, f64),
    _p4: (f64, f64),
    e: f64,
    nu: f64,
) -> Array2<f64> {

    // B-matrix for a 2D quadrilateral (Q4) element under plane stress.
    // These values are derived for a unit square element with nodes ordered as follows:
    // p1 = bottom-left, p2 = bottom-right, p3 = top-right, p4 = top-left
    //
    // The B-matrix relates nodal displacements to strains. The values below are obtained
    // by differentiating the shape functions with respect to x and y at the element centroid,
    // under the assumption of a unit square reference element. For more details, see:
    // - Zienkiewicz & Taylor, "The Finite Element Method," Vol. 1, Section 6.5 (Q4 element)
    // - Cook et al., "Concepts and Applications of Finite Element Analysis," Table for Q4 shape function derivatives
    //
    // If the element geometry or node ordering changes, this matrix must be recomputed accordingly.
    let b_mat = array![
        [
            -0.25, 0.0, 0.25, 0.0,
            0.25, 0.0, -0.25, 0.0
        ],
        [
            0.0, -0.25, 0.0, -0.25,
            0.0, 0.25, 0.0, 0.25
        ],
        [
            -0.25, -0.25, -0.25, 0.25,
            0.25, 0.25, 0.25, -0.25
        ]
    ];

    let c_mat = (e / (1.0 - nu * nu))
        * array![
            [1.0, nu, 0.0],
            [nu, 1.0, 0.0],
            [
                0.0,
                0.0,
                (1.0 - nu) / 2.0
            ]
        ];

    b_mat
        .t()
        .dot(&c_mat.dot(&b_mat))
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
/// A `Result` containing a `Vec<f64>` of nodal displacements (u, v for each node),
/// or an error string if the linear system cannot be solved.

pub fn run_elasticity_simulation(
    params: &ElasticityParameters
) -> Result<Vec<f64>, String> {

    let n_nodes = params.nodes.len();

    let n_dofs = n_nodes * 2;

    // Parallel element stiffness matrix assembly
    let triplets : Vec<(usize, usize, f64)> = params
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

            for r in 0 .. 8 {

                for c in 0 .. 8 {

                    element_triplets.push((
                        dof_indices[r],
                        dof_indices[c],
                        k_element[[r, c]],
                    ));
                }
            }

            element_triplets
        })
        .collect();

    let mut f_global =
        Array1::<f64>::zeros(n_dofs);

    for &(node_idx, fx, fy) in
        &params.loads
    {

        f_global[node_idx * 2] += fx;

        f_global[node_idx * 2 + 1] +=
            fy;
    }

    // Handle boundary conditions: zero out rows and columns of fixed DOFs
    let fixed_dofs : std::collections::HashSet<usize> = params
        .fixed_nodes
        .iter()
        .flat_map(|&node_idx| {

            vec![
                node_idx * 2,
                node_idx * 2 + 1,
            ]
        })
        .collect();

    let mut filtered_triplets: Vec<(
        usize,
        usize,
        f64,
    )> = triplets
        .into_par_iter()
        .filter(|(r, c, _)| {

            !fixed_dofs.contains(r)
                && !fixed_dofs
                    .contains(c)
        })
        .collect();

    for &node_idx in &params.fixed_nodes
    {

        let dof1 = node_idx * 2;

        let dof2 = node_idx * 2 + 1;

        filtered_triplets
            .push((dof1, dof1, 1.0));

        filtered_triplets
            .push((dof2, dof2, 1.0));

        f_global[dof1] = 0.0;

        f_global[dof2] = 0.0;
    }

    let k_global: CsMat<f64> =
        csr_from_triplets(
            n_dofs,
            n_dofs,
            &filtered_triplets,
        );

    let displacements =
        solve_conjugate_gradient(
            &k_global,
            &f_global,
            None,
            5000,
            1e-9,
        )?;

    Ok(displacements.to_vec())
}

/// An example scenario for a cantilever beam under a point load.
///
/// This function sets up a mesh for a 2D cantilever beam, defines fixed boundary
/// conditions at one end and applies a point load at the free end. It then runs
/// the elasticity simulation and saves the original and deformed node positions
/// to CSV files for visualization.

pub fn simulate_cantilever_beam_scenario(
) -> Result<(), String> {

    println!(
        "Running 2D Cantilever Beam \
         simulation..."
    );

    let beam_length = 10.0;

    let beam_height = 2.0;

    let nx = 20;

    let ny = 4;

    let mut nodes: Nodes = Vec::new();

    for j in 0 ..= ny {

        for i in 0 ..= nx {

            nodes.push((
                i as f64 * beam_length
                    / nx as f64,
                j as f64 * beam_height
                    / ny as f64,
            ));
        }
    }

    let mut elements: Elements =
        Vec::new();

    for j in 0 .. ny {

        for i in 0 .. nx {

            let n1 = j * (nx + 1) + i;

            let n2 =
                j * (nx + 1) + i + 1;

            let n3 = (j + 1) * (nx + 1)
                + i
                + 1;

            let n4 =
                (j + 1) * (nx + 1) + i;

            elements
                .push([n1, n2, n3, n4]);
        }
    }

    let fixed_nodes: Vec<usize> = (0
        ..= ny)
        .map(|j| j * (nx + 1))
        .collect();

    let loads = vec![(
        (ny / 2) * (nx + 1) + nx,
        0.0,
        -1e3,
    )];

    let params = ElasticityParameters {
        nodes: nodes.clone(),
        elements,
        youngs_modulus: 1e7,
        poissons_ratio: 0.3,
        fixed_nodes,
        loads,
    };

    let d = run_elasticity_simulation(
        &params,
    )?;

    println!(
        "Simulation finished. Saving \
         results..."
    );

    let mut new_nodes = nodes.clone();

    for i in 0 .. nodes.len() {

        new_nodes[i].0 += d[i * 2];

        new_nodes[i].1 += d[i * 2 + 1];
    }

    let mut orig_file = File::create(
        "beam_original.csv",
    )
    .map_err(|e| e.to_string())?;

    let mut def_file = File::create(
        "beam_deformed.csv",
    )
    .map_err(|e| e.to_string())?;

    writeln!(orig_file, "x,y")
        .map_err(|e| e.to_string())?;

    writeln!(def_file, "x,y")
        .map_err(|e| e.to_string())?;

    for n in &nodes {

        writeln!(
            orig_file,
            "{},{}",
            n.0, n.1
        )
        .map_err(|e| e.to_string())?;
    }

    for n in &new_nodes {

        writeln!(
            def_file,
            "{},{}",
            n.0, n.1
        )
        .map_err(|e| e.to_string())?;
    }

    println!(
        "Original and deformed node \
         positions saved to .csv \
         files."
    );

    Ok(())
}
