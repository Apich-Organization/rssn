// src/physics/physics_sim/linear_elasticity.rs
// 2D Linear Elasticity solver using the Finite Element Method.

use crate::numerical::sparse::{csr_from_triplets, solve_conjugate_gradient};
use ndarray::{array, Array1, Array2};
use sprs::CsMat;
use std::fs::File;
use std::io::Write;

/// Defines the node points of the mesh.
pub type Nodes = Vec<(f64, f64)>;
/// Defines the elements by indexing into the nodes vector.
pub type Elements = Vec<[usize; 4]>; // Quadrilateral elements

/// Parameters for the linear elasticity simulation.
pub struct ElasticityParameters {
    pub nodes: Nodes,
    pub elements: Elements,
    pub youngs_modulus: f64,
    pub poissons_ratio: f64,
    pub fixed_nodes: Vec<usize>,
    pub loads: Vec<(usize, f64, f64)>, // (node_index, fx, fy)
}

/// Calculates the element stiffness matrix for a 2D quadrilateral element (plane stress).
pub(crate) fn element_stiffness_matrix(
    _p1: (f64, f64),
    _p2: (f64, f64),
    _p3: (f64, f64),
    _p4: (f64, f64),
    e: f64,
    nu: f64,
) -> Array2<f64> {
    // This is a simplified formulation using an assumed shape function derivative (B matrix)
    // A full implementation would use isoparametric mapping and numerical quadrature.
    let b_mat = array![
        [-0.25, 0.0, 0.25, 0.0, 0.25, 0.0, -0.25, 0.0],
        [0.0, -0.25, 0.0, -0.25, 0.0, 0.25, 0.0, 0.25],
        [-0.25, -0.25, -0.25, 0.25, 0.25, 0.25, 0.25, -0.25]
    ];

    let c_mat = (e / (1.0 - nu * nu))
        * array![[1.0, nu, 0.0], [nu, 1.0, 0.0], [0.0, 0.0, (1.0 - nu) / 2.0]];

    b_mat.t().dot(&c_mat.dot(&b_mat))
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
pub fn run_elasticity_simulation(params: &ElasticityParameters) -> Result<Vec<f64>, String> {
    let n_nodes = params.nodes.len();
    let n_dofs = n_nodes * 2; // 2 degrees of freedom (u, v) per node

    let mut triplets = Vec::new();
    let mut f_global = Array1::<f64>::zeros(n_dofs);

    // Assemble global stiffness matrix and force vector
    for element in &params.elements {
        let p1 = params.nodes[element[0]];
        let p2 = params.nodes[element[1]];
        let p3 = params.nodes[element[2]];
        let p4 = params.nodes[element[3]];

        let k_element =
            element_stiffness_matrix(p1, p2, p3, p4, params.youngs_modulus, params.poissons_ratio);

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

        for r in 0..8 {
            for c in 0..8 {
                triplets.push((dof_indices[r], dof_indices[c], k_element[[r, c]]));
            }
        }
    }

    // Apply loads
    for &(node_idx, fx, fy) in &params.loads {
        f_global[node_idx * 2] += fx;
        f_global[node_idx * 2 + 1] += fy;
    }

    // Apply boundary conditions (fixed nodes)
    for &node_idx in &params.fixed_nodes {
        let dof1 = node_idx * 2;
        let dof2 = node_idx * 2 + 1;
        // Remove rows/cols corresponding to fixed DOFs and modify RHS
        // A simpler way (used here) is to zero out the row/col and put 1 on the diagonal.
        triplets.retain(|(r, c, _)| *r != dof1 && *c != dof1 && *r != dof2 && *c != dof2);
        f_global[dof1] = 0.0;
        f_global[dof2] = 0.0;
    }
    for &node_idx in &params.fixed_nodes {
        triplets.push((node_idx * 2, node_idx * 2, 1.0));
        triplets.push((node_idx * 2 + 1, node_idx * 2 + 1, 1.0));
    }

    // Create sparse matrix and solve
    let k_global: CsMat<f64> = csr_from_triplets(n_dofs, n_dofs, &triplets);
    let displacements = solve_conjugate_gradient(&k_global, &f_global, None, 5000, 1e-9)?;

    Ok(displacements.to_vec())
}

/// An example scenario for a cantilever beam under a point load.
///
/// This function sets up a mesh for a 2D cantilever beam, defines fixed boundary
/// conditions at one end and applies a point load at the free end. It then runs
/// the elasticity simulation and saves the original and deformed node positions
/// to CSV files for visualization.
pub fn simulate_cantilever_beam_scenario() {
    println!("Running 2D Cantilever Beam simulation...");

    // 1. Create the mesh (nodes and elements)
    let beam_length = 10.0;
    let beam_height = 2.0;
    let nx = 20; // elements in x
    let ny = 4; // elements in y

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

    // 2. Define boundary conditions and loads
    let fixed_nodes: Vec<usize> = (0..=ny).map(|j| j * (nx + 1)).collect();
    let loads = vec![((ny / 2) * (nx + 1) + nx, 0.0, -1e3)]; // Point load at the middle of the free end

    let params = ElasticityParameters {
        nodes: nodes.clone(),
        elements,
        youngs_modulus: 1e7,
        poissons_ratio: 0.3,
        fixed_nodes,
        loads,
    };

    // 3. Run simulation
    match run_elasticity_simulation(&params) {
        Ok(d) => {
            println!("Simulation finished. Saving results...");
            let mut new_nodes = nodes.clone();
            for i in 0..nodes.len() {
                new_nodes[i].0 += d[i * 2];
                new_nodes[i].1 += d[i * 2 + 1];
            }

            // Save as simple CSV for easy plotting
            let mut orig_file = File::create("beam_original.csv").unwrap();
            let mut def_file = File::create("beam_deformed.csv").unwrap();
            writeln!(orig_file, "x,y").unwrap();
            writeln!(def_file, "x,y").unwrap();
            nodes.iter().for_each(|n| {
                writeln!(orig_file, "{},{}", n.0, n.1).unwrap();
            });
            new_nodes.iter().for_each(|n| {
                writeln!(def_file, "{},{}", n.0, n.1).unwrap();
            });

            println!("Original and deformed node positions saved to .csv files.");
        }
        Err(e) => eprintln!("An error occurred: {}", e),
    }
}
