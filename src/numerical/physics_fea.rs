//! # Numerical Finite Element Analysis (FEA)
//!
//! This module provides numerical methods for Finite Element Analysis (FEA),
//! a powerful computational technique for solving problems in structural mechanics,
//! heat transfer, fluid dynamics, and electromagnetics.
//!
//! ## Element Types
//!
//! ### 1D Elements
//! - **Linear Bar/Truss**: Axial deformation only
//! - **Euler-Bernoulli Beam**: Bending and axial deformation
//!
//! ### 2D Elements
//! - **CST (Constant Strain Triangle)**: 3-node triangular element
//! - **LST (Linear Strain Triangle)**: 6-node triangular element
//! - **Quadrilateral**: 4-node bilinear element
//!
//! ## Analysis Types
//! - Static structural analysis
//! - Steady-state thermal analysis
//! - Natural frequency (modal) analysis
//!
//! ## Features
//! - Global stiffness matrix assembly
//! - Boundary condition application
//! - Nodal displacement and stress computation
//! - Element stress recovery
//!
//! ## Example
//!
//! ```rust
//! 
//! use rssn::numerical::physics_fea::*;
//!
//! // Create a simple 1D bar
//! let element = LinearElement1D {
//!     length : 1.0,
//!     youngs_modulus : 200e9, // Steel
//!     area : 0.001,
//! };
//!
//! let k = element.local_stiffness_matrix();
//! ```

use serde::Deserialize;
use serde::Serialize;

use crate::numerical::matrix::Matrix;

// ============================================================================
// Material Properties
// ============================================================================

/// Material properties for structural analysis.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]

pub struct Material {
    /// Young's modulus (Pa)
    pub youngs_modulus : f64,
    /// Poisson's ratio (dimensionless)
    pub poissons_ratio : f64,
    /// Density (kg/m³)
    pub density : f64,
    /// Thermal conductivity (W/(m·K))
    pub thermal_conductivity : f64,
    /// Coefficient of thermal expansion (1/K)
    pub thermal_expansion : f64,
    /// Yield strength (Pa)
    pub yield_strength : f64,
}

impl Material {
    /// Creates a new material with specified properties.
    #[must_use]

    pub const fn new(
        youngs_modulus : f64,
        poissons_ratio : f64,
        density : f64,
        thermal_conductivity : f64,
        thermal_expansion : f64,
        yield_strength : f64,
    ) -> Self {

        Self {
            youngs_modulus,
            poissons_ratio,
            density,
            thermal_conductivity,
            thermal_expansion,
            yield_strength,
        }
    }

    /// Creates steel material with typical properties.
    #[must_use]

    pub const fn steel() -> Self {

        Self {
            youngs_modulus : 200e9,
            poissons_ratio : 0.3,
            density : 7850.0,
            thermal_conductivity : 50.0,
            thermal_expansion : 12e-6,
            yield_strength : 250e6,
        }
    }

    /// Creates aluminum material with typical properties.
    #[must_use]

    pub const fn aluminum() -> Self {

        Self {
            youngs_modulus : 70e9,
            poissons_ratio : 0.33,
            density : 2700.0,
            thermal_conductivity : 205.0,
            thermal_expansion : 23e-6,
            yield_strength : 270e6,
        }
    }

    /// Creates copper material with typical properties.
    #[must_use]

    pub const fn copper() -> Self {

        Self {
            youngs_modulus : 117e9,
            poissons_ratio : 0.34,
            density : 8960.0,
            thermal_conductivity : 401.0,
            thermal_expansion : 17e-6,
            yield_strength : 70e6,
        }
    }

    /// Shear modulus G = E / (2(1+ν))
    #[must_use]

    pub fn shear_modulus(&self) -> f64 {

        self.youngs_modulus / (2.0 * (1.0 + self.poissons_ratio))
    }

    /// Bulk modulus K = E / (3(1-2ν))
    #[must_use]

    pub fn bulk_modulus(&self) -> f64 {

        self.youngs_modulus / (3.0 * 2.0f64.mul_add(-self.poissons_ratio, 1.0))
    }
}

// ============================================================================
// Node Definition
// ============================================================================

/// A node in the finite element mesh.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]

pub struct Node2D {
    pub id : usize,
    pub x : f64,
    pub y : f64,
}

impl Node2D {
    /// Creates a new 2D node.
    #[must_use]

    pub const fn new(
        id : usize,
        x : f64,
        y : f64,
    ) -> Self {

        Self {
            id,
            x,
            y,
        }
    }

    /// Distance to another node.
    #[must_use]

    pub fn distance_to(
        &self,
        other : &Self,
    ) -> f64 {

        let dx = other.x - self.x;

        let dy = other.y - self.y;

        dx.hypot(dy)
    }
}

/// A node in 3D space.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]

pub struct Node3D {
    pub id : usize,
    pub x : f64,
    pub y : f64,
    pub z : f64,
}

impl Node3D {
    /// Creates a new 3D node.
    #[must_use]

    pub const fn new(
        id : usize,
        x : f64,
        y : f64,
        z : f64,
    ) -> Self {

        Self {
            id,
            x,
            y,
            z,
        }
    }
}

// ============================================================================
// 1D Elements
// ============================================================================

/// Represents a 1D linear finite element.

pub struct LinearElement1D {
    pub length : f64,
    pub youngs_modulus : f64,
    pub area : f64,
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
    #[must_use]

    pub fn local_stiffness_matrix(&self) -> Matrix<f64> {

        let k = self.youngs_modulus * self.area / self.length;

        Matrix::new(
            2,
            2,
            vec![k, -k, -k, k],
        )
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
#[must_use]

pub fn assemble_global_stiffness_matrix(
    num_nodes : usize,
    elements : &[(
        Matrix<f64>,
        usize,
        usize,
    )],
) -> Matrix<f64> {

    let mut global_k = Matrix::zeros(num_nodes, num_nodes);

    for (local_k, n1, n2) in elements {

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
    mut global_k : Matrix<f64>,
    mut forces : Vec<f64>,
    fixed_dofs : &[(usize, f64)],
) -> Result<Vec<f64>, String> {

    let n = global_k.rows();

    if forces.len() != n {

        return Err("Force vector dimension mismatch.".to_string());
    }

    for &(node_idx, prescribed_disp) in fixed_dofs {

        for (i, var) in forces
            .iter_mut()
            .enumerate()
            .take(n)
        {

            if i != node_idx {

                *var -= global_k.get(i, node_idx) * prescribed_disp;
            }
        }

        for i in 0 .. n {

            *global_k.get_mut(node_idx, i) = 0.0;

            *global_k.get_mut(i, node_idx) = 0.0;
        }

        *global_k.get_mut(node_idx, node_idx) = 1.0;

        forces[node_idx] = prescribed_disp;
    }

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

// ============================================================================
// 2D Triangular Element (CST - Constant Strain Triangle)
// ============================================================================

/// Constant Strain Triangle (CST) element for 2D plane stress/strain analysis.
#[derive(Debug, Clone)]

pub struct TriangleElement2D {
    /// Node indices (3 nodes)
    pub nodes : [usize; 3],
    /// Node coordinates
    pub coords : [(f64, f64); 3],
    /// Thickness
    pub thickness : f64,
    /// Material
    pub material : Material,
    /// Plane stress (true) or plane strain (false)
    pub plane_stress : bool,
}

impl TriangleElement2D {
    /// Creates a new triangular element.
    #[must_use]

    pub const fn new(
        nodes : [usize; 3],
        coords : [(f64, f64); 3],
        thickness : f64,
        material : Material,
        plane_stress : bool,
    ) -> Self {

        Self {
            nodes,
            coords,
            thickness,
            material,
            plane_stress,
        }
    }

    /// Calculates the area of the triangle.
    #[must_use]

    pub fn area(&self) -> f64 {

        let (x1, y1) = self.coords[0];

        let (x2, y2) = self.coords[1];

        let (x3, y3) = self.coords[2];

        0.5 * (x2 - x1).mul_add(y3 - y1, -((x3 - x1) * (y2 - y1))).abs()
    }

    /// Computes the constitutive (D) matrix for plane stress or plane strain.
    #[must_use]

    pub fn constitutive_matrix(&self) -> Matrix<f64> {

        let e = self
            .material
            .youngs_modulus;

        let nu = self
            .material
            .poissons_ratio;

        if self.plane_stress {

            let factor = e / nu.mul_add(-nu, 1.0);

            Matrix::new(
                3,
                3,
                vec![
                    factor,
                    factor * nu,
                    0.0,
                    factor * nu,
                    factor,
                    0.0,
                    0.0,
                    0.0,
                    factor * (1.0 - nu) / 2.0,
                ],
            )
        } else {

            // Plane strain
            let factor = e / ((1.0 + nu) * 2.0f64.mul_add(-nu, 1.0));

            Matrix::new(
                3,
                3,
                vec![
                    factor * (1.0 - nu),
                    factor * nu,
                    0.0,
                    factor * nu,
                    factor * (1.0 - nu),
                    0.0,
                    0.0,
                    0.0,
                    factor * 2.0f64.mul_add(-nu, 1.0) / 2.0,
                ],
            )
        }
    }

    /// Computes the strain-displacement (B) matrix.
    #[must_use]

    pub fn b_matrix(&self) -> Matrix<f64> {

        let (x1, y1) = self.coords[0];

        let (x2, y2) = self.coords[1];

        let (x3, y3) = self.coords[2];

        let a = self.area();

        let two_a = 2.0 * a;

        // Beta and gamma coefficients
        let b1 = (y2 - y3) / two_a;

        let b2 = (y3 - y1) / two_a;

        let b3 = (y1 - y2) / two_a;

        let g1 = (x3 - x2) / two_a;

        let g2 = (x1 - x3) / two_a;

        let g3 = (x2 - x1) / two_a;

        // B matrix (3 x 6)
        Matrix::new(
            3,
            6,
            vec![
                b1, 0.0, b2, 0.0, b3, 0.0, // Row 1: ∂u/∂x
                0.0, g1, 0.0, g2, 0.0, g3, // Row 2: ∂v/∂y
                g1, b1, g2, b2, g3, b3, // Row 3: ∂u/∂y + ∂v/∂x
            ],
        )
    }

    /// Computes the local stiffness matrix (6x6).
    #[must_use]

    pub fn local_stiffness_matrix(&self) -> Matrix<f64> {

        let b = self.b_matrix();

        let d = self.constitutive_matrix();

        let a = self.area();

        let t = self.thickness;

        // K = t * A * B^T * D * B
        let bt = b.transpose();

        let db = d * b;

        let btdb = bt * db;

        btdb * (t * a)
    }

    /// Computes element stresses given nodal displacements.
    #[must_use]

    pub fn compute_stress(
        &self,
        displacements : &[f64],
    ) -> Vec<f64> {

        assert_eq!(
            displacements.len(),
            6,
            "Need 6 displacement values for CST element"
        );

        let b = self.b_matrix();

        let d = self.constitutive_matrix();

        // Strain = B * u
        let mut strain = [0.0; 3];

        for i in 0 .. 3 {

            for j in 0 .. 6 {

                strain[i] += *b.get(i, j) * displacements[j];
            }
        }

        // Stress = D * strain
        let mut stress = vec![0.0; 3];

        for i in 0 .. 3 {

            for j in 0 .. 3 {

                stress[i] += *d.get(i, j) * strain[j];
            }
        }

        stress
    }

    /// Computes von Mises stress from stress components [σx, σy, τxy].
    #[must_use]

    pub fn von_mises_stress(stress : &[f64]) -> f64 {

        assert_eq!(stress.len(), 3);

        let sx = stress[0];

        let sy = stress[1];

        let txy = stress[2];

        ((sx.mul_add(sx, -(sx * sy)) + sy * sy) + 3.0 * txy * txy).sqrt()
    }
}

// ============================================================================
// Euler-Bernoulli Beam Element
// ============================================================================

/// 2D Euler-Bernoulli beam element with axial and bending DOFs.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]

pub struct BeamElement2D {
    /// Element length
    pub length : f64,
    /// Young's modulus
    pub youngs_modulus : f64,
    /// Cross-sectional area
    pub area : f64,
    /// Second moment of area (moment of inertia)
    pub moment_of_inertia : f64,
    /// Angle of orientation (radians from global x-axis)
    pub angle : f64,
}

impl BeamElement2D {
    /// Creates a new beam element.
    #[must_use]

    pub const fn new(
        length : f64,
        youngs_modulus : f64,
        area : f64,
        moment_of_inertia : f64,
        angle : f64,
    ) -> Self {

        Self {
            length,
            youngs_modulus,
            area,
            moment_of_inertia,
            angle,
        }
    }

    /// Computes the local stiffness matrix in local coordinates (6x6).
    /// DOFs: [u1, v1, θ1, u2, v2, θ2]
    #[must_use]

    pub fn local_stiffness_matrix(&self) -> Matrix<f64> {

        let e = self.youngs_modulus;

        let a = self.area;

        let i = self.moment_of_inertia;

        let l = self.length;

        let l2 = l * l;

        let l3 = l * l * l;

        let ea_l = e * a / l;

        let ei_l3 = e * i / l3;

        let ei_l2 = e * i / l2;

        let ei_l = e * i / l;

        // Local stiffness matrix
        Matrix::new(
            6,
            6,
            vec![
                ea_l,
                0.0,
                0.0,
                -ea_l,
                0.0,
                0.0,
                0.0,
                12.0 * ei_l3,
                6.0 * ei_l2,
                0.0,
                -12.0 * ei_l3,
                6.0 * ei_l2,
                0.0,
                6.0 * ei_l2,
                4.0 * ei_l,
                0.0,
                -6.0 * ei_l2,
                2.0 * ei_l,
                -ea_l,
                0.0,
                0.0,
                ea_l,
                0.0,
                0.0,
                0.0,
                -12.0 * ei_l3,
                -6.0 * ei_l2,
                0.0,
                12.0 * ei_l3,
                -6.0 * ei_l2,
                0.0,
                6.0 * ei_l2,
                2.0 * ei_l,
                0.0,
                -6.0 * ei_l2,
                4.0 * ei_l,
            ],
        )
    }

    /// Computes the transformation matrix from local to global coordinates.
    #[must_use]

    pub fn transformation_matrix(&self) -> Matrix<f64> {

        let c = self.angle.cos();

        let s = self.angle.sin();

        Matrix::new(
            6,
            6,
            vec![
                c, s, 0.0, 0.0, 0.0, 0.0, -s, c, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, c, s, 0.0, 0.0, 0.0, 0.0, -s, c, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
            ],
        )
    }

    /// Computes the global stiffness matrix.
    #[must_use]

    pub fn global_stiffness_matrix(&self) -> Matrix<f64> {

        let k_local = self.local_stiffness_matrix();

        let t = self.transformation_matrix();

        let tt = t.transpose();

        let kt = k_local * t;

        tt * kt
    }

    /// Computes the consistent mass matrix.
    #[must_use]

    pub fn mass_matrix(
        &self,
        density : f64,
    ) -> Matrix<f64> {

        let a = self.area;

        let l = self.length;

        let m = density * a * l;

        // Consistent mass matrix (simplified, diagonal terms for lumped mass)
        let factor = m / 420.0;

        Matrix::new(
            6,
            6,
            vec![
                140.0 * factor,
                0.0,
                0.0,
                70.0 * factor,
                0.0,
                0.0,
                0.0,
                156.0 * factor,
                22.0 * l * factor,
                0.0,
                54.0 * factor,
                -13.0 * l * factor,
                0.0,
                22.0 * l * factor,
                4.0 * l * l * factor,
                0.0,
                13.0 * l * factor,
                -3.0 * l * l * factor,
                70.0 * factor,
                0.0,
                0.0,
                140.0 * factor,
                0.0,
                0.0,
                0.0,
                54.0 * factor,
                13.0 * l * factor,
                0.0,
                156.0 * factor,
                -22.0 * l * factor,
                0.0,
                -13.0 * l * factor,
                -3.0 * l * l * factor,
                0.0,
                -22.0 * l * factor,
                4.0 * l * l * factor,
            ],
        )
    }
}

// ============================================================================
// Thermal Element
// ============================================================================

/// 1D thermal element for heat conduction analysis.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]

pub struct ThermalElement1D {
    /// Element length
    pub length : f64,
    /// Thermal conductivity
    pub conductivity : f64,
    /// Cross-sectional area
    pub area : f64,
}

impl ThermalElement1D {
    /// Creates a new 1D thermal element.
    #[must_use]

    pub const fn new(
        length : f64,
        conductivity : f64,
        area : f64,
    ) -> Self {

        Self {
            length,
            conductivity,
            area,
        }
    }

    /// Computes the conductivity (stiffness) matrix.
    #[must_use]

    pub fn conductivity_matrix(&self) -> Matrix<f64> {

        let k = self.conductivity * self.area / self.length;

        Matrix::new(
            2,
            2,
            vec![k, -k, -k, k],
        )
    }
}

/// 2D triangular thermal element.
#[derive(Debug, Clone)]

pub struct ThermalTriangle2D {
    /// Node coordinates
    pub coords : [(f64, f64); 3],
    /// Thickness
    pub thickness : f64,
    /// Thermal conductivity
    pub conductivity : f64,
}

impl ThermalTriangle2D {
    /// Creates a new 2D triangular thermal element.
    #[must_use]

    pub const fn new(
        coords : [(f64, f64); 3],
        thickness : f64,
        conductivity : f64,
    ) -> Self {

        Self {
            coords,
            thickness,
            conductivity,
        }
    }

    /// Calculates the area of the triangle.
    #[must_use]

    pub fn area(&self) -> f64 {

        let (x1, y1) = self.coords[0];

        let (x2, y2) = self.coords[1];

        let (x3, y3) = self.coords[2];

        0.5 * (x2 - x1).mul_add(y3 - y1, -((x3 - x1) * (y2 - y1))).abs()
    }

    /// Computes the conductivity matrix (3x3).
    #[must_use]

    pub fn conductivity_matrix(&self) -> Matrix<f64> {

        let (x1, y1) = self.coords[0];

        let (x2, y2) = self.coords[1];

        let (x3, y3) = self.coords[2];

        let a = self.area();

        let k = self.conductivity;

        let t = self.thickness;

        // B coefficients (temperature gradient)
        let b1 = y2 - y3;

        let b2 = y3 - y1;

        let b3 = y1 - y2;

        let c1 = x3 - x2;

        let c2 = x1 - x3;

        let c3 = x2 - x1;

        let factor = k * t / (4.0 * a);

        Matrix::new(
            3,
            3,
            vec![
                factor * b1.mul_add(b1, c1 * c1),
                factor * b1.mul_add(b2, c1 * c2),
                factor * b1.mul_add(b3, c1 * c3),
                factor * b2.mul_add(b1, c2 * c1),
                factor * b2.mul_add(b2, c2 * c2),
                factor * b2.mul_add(b3, c2 * c3),
                factor * b3.mul_add(b1, c3 * c1),
                factor * b3.mul_add(b2, c3 * c2),
                factor * b3.mul_add(b3, c3 * c3),
            ],
        )
    }
}

// ============================================================================
// Utility Functions
// ============================================================================

/// Applies boundary conditions to a global stiffness matrix and force vector.
/// Uses the penalty method for prescribed DOFs.

pub fn apply_boundary_conditions_penalty(
    global_k : &mut Matrix<f64>,
    forces : &mut [f64],
    fixed_dofs : &[(usize, f64)],
    penalty : f64,
) {

    for &(dof, value) in fixed_dofs {

        *global_k.get_mut(dof, dof) += penalty;

        forces[dof] += penalty * value;
    }
}

/// Assembles a 2D global stiffness matrix from triangular elements.
#[must_use]

pub fn assemble_2d_stiffness_matrix(
    num_dofs : usize,
    elements : &[(
        Matrix<f64>,
        [usize; 6],
    )],
) -> Matrix<f64> {

    let mut global_k = Matrix::zeros(num_dofs, num_dofs);

    for (local_k, dof_map) in elements {

        for i in 0 .. 6 {

            for j in 0 .. 6 {

                let gi = dof_map[i];

                let gj = dof_map[j];

                *global_k.get_mut(gi, gj) += local_k.get(i, j);
            }
        }
    }

    global_k
}

/// Calculates element strains from nodal displacements.
#[must_use]

pub fn compute_element_strain(
    b_matrix : &Matrix<f64>,
    displacements : &[f64],
) -> Vec<f64> {

    let mut strain = vec![0.0; b_matrix.rows()];

    for i in 0 .. b_matrix.rows() {

        for j in 0 .. b_matrix.cols() {

            strain[i] += *b_matrix.get(i, j) * displacements[j];
        }
    }

    strain
}

/// Converts stress to principal stresses.
/// Returns (sigma1, sigma2, angle) where sigma1 >= sigma2.
#[must_use]

pub fn principal_stresses(stress : &[f64]) -> (f64, f64, f64) {

    assert_eq!(
        stress.len(),
        3,
        "Need [σx, σy, τxy]"
    );

    let sx = stress[0];

    let sy = stress[1];

    let txy = stress[2];

    let avg = f64::midpoint(sx, sy);

    let diff = (sx - sy) / 2.0;

    let r = diff.hypot(txy);

    let sigma1 = avg + r;

    let sigma2 = avg - r;

    let angle = 0.5 * txy.atan2(diff);

    (
        sigma1,
        sigma2,
        angle,
    )
}

/// Computes the maximum shear stress from principal stresses.
#[must_use]

pub fn max_shear_stress(
    sigma1 : f64,
    sigma2 : f64,
) -> f64 {

    (sigma1 - sigma2).abs() / 2.0
}

/// Safety factor based on von Mises criterion.
#[must_use]

pub fn safety_factor_von_mises(
    stress : &[f64],
    yield_strength : f64,
) -> f64 {

    let vm = TriangleElement2D::von_mises_stress(stress);

    if vm > 0.0 {

        yield_strength / vm
    } else {

        f64::INFINITY
    }
}

/// Creates a simple rectangular mesh of triangular elements.
/// Returns (nodes, elements) where elements are node indices.
#[must_use]

pub fn create_rectangular_mesh(
    width : f64,
    height : f64,
    nx : usize,
    ny : usize,
) -> (
    Vec<Node2D>,
    Vec<[usize; 3]>,
) {

    let mut nodes = Vec::with_capacity((nx + 1) * (ny + 1));

    let dx = width / nx as f64;

    let dy = height / ny as f64;

    // Create nodes
    for j in 0 ..= ny {

        for i in 0 ..= nx {

            let id = j * (nx + 1) + i;

            nodes.push(Node2D::new(
                id,
                i as f64 * dx,
                j as f64 * dy,
            ));
        }
    }

    // Create triangular elements (2 per rectangle)
    let mut elements = Vec::with_capacity(2 * nx * ny);

    for j in 0 .. ny {

        for i in 0 .. nx {

            let n1 = j * (nx + 1) + i;

            let n2 = n1 + 1;

            let n3 = n1 + (nx + 1);

            let n4 = n3 + 1;

            // Lower-left triangle
            elements.push([n1, n2, n3]);

            // Upper-right triangle
            elements.push([n2, n4, n3]);
        }
    }

    (nodes, elements)
}

/// Refines a triangular mesh by subdividing each triangle into 4 smaller triangles.
#[must_use]

pub fn refine_mesh(
    nodes : &[Node2D],
    elements : &[[usize; 3]],
) -> (
    Vec<Node2D>,
    Vec<[usize; 3]>,
) {

    let mut new_nodes = nodes.to_vec();

    let mut new_elements = Vec::with_capacity(elements.len() * 4);

    let mut edge_midpoints : std::collections::HashMap<(usize, usize), usize> =
        std::collections::HashMap::new();

    for &[n1, n2, n3] in elements {

        // Get or create midpoint nodes
        let mut get_midpoint = |a : usize, b : usize| -> usize {

            let key = if a < b {

                (a, b)
            } else {

                (b, a)
            };

            if let Some(&idx) = edge_midpoints.get(&key) {

                idx
            } else {

                let mid_x = f64::midpoint(new_nodes[a].x, new_nodes[b].x);

                let mid_y = f64::midpoint(new_nodes[a].y, new_nodes[b].y);

                let idx = new_nodes.len();

                new_nodes.push(Node2D::new(
                    idx, mid_x, mid_y,
                ));

                edge_midpoints.insert(key, idx);

                idx
            }
        };

        let m12 = get_midpoint(n1, n2);

        let m23 = get_midpoint(n2, n3);

        let m31 = get_midpoint(n3, n1);

        // Create 4 new triangles
        new_elements.push([n1, m12, m31]);

        new_elements.push([m12, n2, m23]);

        new_elements.push([m31, m23, n3]);

        new_elements.push([m12, m23, m31]);
    }

    (
        new_nodes,
        new_elements,
    )
}
