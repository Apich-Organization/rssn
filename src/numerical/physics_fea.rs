//! # Numerical Finite Element Analysis (FEA)
//!
//! This module provides numerical methods for Finite Element Analysis (FEA).
//! It includes basic implementations for solving static structural problems,
//! such as assembling global stiffness matrices and solving for displacements.

use crate::numerical::matrix::Matrix;

/// Represents a 1D linear finite element.
pub struct LinearElement1D {
    pub length: f64,
    pub youngs_modulus: f64,
    pub area: f64,
}

impl LinearElement1D {
    /// Computes the local stiffness matrix for a 1D linear element.
    ///
    /// The local stiffness matrix relates the forces and displacements at the nodes
    /// of a single element. For a 1D linear element, it is given by:
    /// `[EA/L  -EA/L]`
    /// `[-EA/L  EA/L]`
    ///
    /// # Returns
    /// A `Matrix<f64>` representing the local stiffness matrix.
    pub fn local_stiffness_matrix(&self) -> Matrix<f64> {
        let k = self.youngs_modulus * self.area / self.length;
        Matrix::new(2, 2, vec![k, -k, -k, k])
    }
}

/// Assembles the global stiffness matrix for a 1D structure composed of linear elements.
///
/// This function takes a list of local stiffness matrices and their connectivity
/// to construct the global stiffness matrix for the entire structure.
///
/// # Arguments
/// * `num_nodes` - The total number of nodes in the structure.
/// * `elements` - A vector of tuples `(element_matrix, node1_idx, node2_idx)`,
///   where `element_matrix` is the local stiffness matrix and `node1_idx`, `node2_idx`
///   are the global indices of the nodes connected by this element.
///
/// # Returns
/// A `Matrix<f64>` representing the global stiffness matrix.
pub fn assemble_global_stiffness_matrix(
    num_nodes: usize,
    elements: &[(Matrix<f64>, usize, usize)],
) -> Matrix<f64> {
    let mut global_k = Matrix::zeros(num_nodes, num_nodes);

    for (local_k, n1, n2) in elements {
        // Add contributions from local stiffness matrix to global matrix
        *global_k.get_mut(*n1, *n1) += local_k.get(0, 0);
        *global_k.get_mut(*n1, *n2) += local_k.get(0, 1);
        *global_k.get_mut(*n2, *n1) += local_k.get(1, 0);
        *global_k.get_mut(*n2, *n2) += local_k.get(1, 1);
    }
    global_k
}

/// Solves a static structural problem for displacements.
///
/// This function takes the global stiffness matrix, applied forces, and boundary conditions
/// to solve for the nodal displacements. It modifies the global stiffness matrix and force
/// vector to incorporate Dirichlet boundary conditions.
///
/// # Arguments
/// * `global_k` - The global stiffness matrix.
/// * `forces` - The global force vector.
/// * `fixed_dofs` - A vector of tuples `(node_idx, prescribed_displacement)` for fixed degrees of freedom.
///
/// # Returns
/// A `Result` containing a `Vec<f64>` of nodal displacements, or an error string if the system is singular.
pub fn solve_static_structural(
    mut global_k: Matrix<f64>,
    mut forces: Vec<f64>,
    fixed_dofs: &[(usize, f64)],
) -> Result<Vec<f64>, String> {
    let n = global_k.rows();
    if forces.len() != n {
        return Err("Force vector dimension mismatch.".to_string());
    }

    // Apply boundary conditions (Dirichlet)
    for &(node_idx, prescribed_disp) in fixed_dofs {
        // Modify force vector
        for i in 0..n {
            if i != node_idx {
                forces[i] -= global_k.get(i, node_idx) * prescribed_disp;
            }
        }
        // Modify stiffness matrix
        for i in 0..n {
            *global_k.get_mut(node_idx, i) = 0.0;
            *global_k.get_mut(i, node_idx) = 0.0;
        }
        *global_k.get_mut(node_idx, node_idx) = 1.0;
        forces[node_idx] = prescribed_disp;
    }

    // Solve the modified system K * U = F
    let solution = crate::numerical::solve::solve_linear_system(&global_k, &forces)?;

    if let crate::numerical::solve::LinearSolution::Unique(u) = solution {
        Ok(u)
    } else {
        Err(
            "System is singular or has infinite solutions after applying boundary conditions."
                .to_string(),
        )
    }
}
