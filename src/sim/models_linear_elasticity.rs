//! # Linear elasticity (2D quadrilaterals)
//!
//! | Case | Dimension | Elements | Boundaries |
//! |---|---|---|---|
//! | plane stress | 2D | bilinear quadrilaterals | clamped nodes, nodal loads |
//!
//! Plane strain, 3D tetrahedra, P2 elements and mixed displacement/traction
//! markers are in [`crate::sim::elasticity_fem`].

use std::fs::File;
use std::io::Write;
use std::path::Path;

use ndarray::Array1;
use ndarray::Array2;
use serde::Deserialize;
use serde::Serialize;

use crate::kernels::sparse::SparseMat;
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
    let block = stiffness([p1, p2, p3, p4], e, nu, false).unwrap_or([[f64::NAN; 8]; 8]);
    Array2::from_shape_fn((8, 8), |(i, j)| block[i][j])
}

const GAUSS: f64 = 0.577_350_269_189_625_8;
const XI: [f64; 4] = [-1.0, 1.0, 1.0, -1.0];
const ETA: [f64; 4] = [-1.0, -1.0, 1.0, 1.0];

fn local_points(points: [(f64, f64); 4]) -> Result<[[f64; 2]; 4], &'static str> {
    if points
        .iter()
        .any(|&(x, y)| !x.is_finite() || !y.is_finite())
    {
        return Err("Q4 vertices must be finite");
    }
    let origin = [points[0].0, points[0].1];
    let mut local = points.map(|(x, y)| [x - origin[0], y - origin[1]]);
    if local.iter().flatten().any(|x| !x.is_finite()) {
        let mid = |coordinate: fn((f64, f64)) -> f64| {
            let mut range = [f64::INFINITY, f64::NEG_INFINITY];
            for point in points {
                let x = coordinate(point);
                range = [range[0].min(x), range[1].max(x)];
            }
            f64::midpoint(range[0], range[1])
        };
        let center = [mid(|p| p.0), mid(|p| p.1)];
        local = points.map(|(x, y)| [x - center[0], y - center[1]]);
    }
    let scale = local
        .iter()
        .flatten()
        .fold(0.0_f64, |m, value| m.max(value.abs()));
    if scale == 0.0 {
        return Err("Q4 vertices coincide");
    }
    // A common in-plane scale cancels between the gradients and the area.
    if !(1e-135..=1e135).contains(&scale) {
        for value in local.iter_mut().flatten() {
            *value /= scale;
        }
    }
    Ok(local)
}

fn material(
    young: f64,
    poisson: f64,
) -> Result<(), &'static str> {
    if !young.is_finite() || young <= 0.0 || !(-1.0 < poisson && poisson < 0.5) {
        return Err("Q4 material requires finite E > 0 and -1 < nu < 0.5");
    }
    Ok(())
}

fn stiffness(
    points: [(f64, f64); 4],
    young: f64,
    poisson: f64,
    checked: bool,
) -> Result<[[f64; 8]; 8], &'static str> {
    let points = local_points(points)?;
    let normal = poisson.mul_add(-poisson, 1.0).recip();
    let transverse = poisson * normal;
    let shear = 0.5 * (1.0 - poisson) * normal;
    let mut k = [[0.0; 8]; 8];

    for xi in [-GAUSS, GAUSS] {
        for eta in [-GAUSS, GAUSS] {
            let dxi = std::array::from_fn::<_, 4, _>(|i| 0.25 * XI[i] * (1.0 + ETA[i] * eta));
            let deta = std::array::from_fn::<_, 4, _>(|i| 0.25 * ETA[i] * (1.0 + XI[i] * xi));
            let mut j = [[0.0; 2]; 2];
            for (i, p) in points.iter().enumerate() {
                j[0][0] += dxi[i] * p[0];
                j[0][1] += dxi[i] * p[1];
                j[1][0] += deta[i] * p[0];
                j[1][1] += deta[i] * p[1];
            }
            let det = j[0][0].mul_add(j[1][1], -(j[0][1] * j[1][0]));
            if !det.is_finite() || det == 0.0 || (checked && det < 0.0) {
                return Err(
                    "Q4 mapping must have a finite nonzero Jacobian with positive orientation",
                );
            }
            let area_weight = det.abs().sqrt();
            let mut strain = [[0.0; 8]; 3];
            for i in 0..4 {
                // Combine inverse-Jacobian and quadrature factors before multiplying by the material.
                let dx = j[1][1].mul_add(dxi[i], -j[0][1] * deta[i]) / area_weight;
                let dy = (-j[1][0]).mul_add(dxi[i], j[0][0] * deta[i]) / area_weight;
                strain[0][2 * i] = dx;
                strain[1][2 * i + 1] = dy;
                strain[2][2 * i] = dy;
                strain[2][2 * i + 1] = dx;
            }
            let stress: [[f64; 8]; 3] = [
                std::array::from_fn(|i| normal * strain[0][i] + transverse * strain[1][i]),
                std::array::from_fn(|i| transverse * strain[0][i] + normal * strain[1][i]),
                std::array::from_fn(|i| shear * strain[2][i]),
            ];
            for (row, entries) in k.iter_mut().enumerate() {
                for (col, slot) in entries.iter_mut().enumerate().skip(row) {
                    *slot += strain[0][row] * stress[0][col]
                        + strain[1][row] * stress[1][col]
                        + strain[2][row] * stress[2][col];
                }
            }
        }
    }
    for row in 0..8 {
        for col in row..8 {
            let entry = k[row][col] * young;
            if checked && !entry.is_finite() {
                return Err("Q4 stiffness is outside the floating-point range");
            }
            k[row][col] = entry;
            k[col][row] = entry;
        }
    }
    Ok(k)
}

pub(crate) fn checked_stiffness(
    points: [(f64, f64); 4],
    young: f64,
    poisson: f64,
) -> Result<[[f64; 8]; 8], &'static str> {
    material(young, poisson)?;
    stiffness(points, young, poisson, true)
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
    let n_dofs = n_nodes.checked_mul(2).ok_or("Q4 degree count overflow")?;
    material(params.youngs_modulus, params.poissons_ratio)?;
    if params.elements.iter().flatten().any(|&i| i >= n_nodes)
        || params.fixed_nodes.iter().any(|&i| i >= n_nodes)
        || params
            .loads
            .iter()
            .any(|&(i, fx, fy)| i >= n_nodes || !fx.is_finite() || !fy.is_finite())
    {
        return Err(
            "Q4 connectivity, fixed nodes and loads must refer to finite valid data".into(),
        );
    }
    let mut fixed = vec![false; n_nodes];
    for &node in &params.fixed_nodes {
        fixed[node] = true;
    }
    let mut f_global = Array1::<f64>::zeros(n_dofs);
    for &(node, fx, fy) in &params.loads {
        if !fixed[node] {
            f_global[node * 2] += fx;
            f_global[node * 2 + 1] += fy;
        }
    }
    if f_global.iter().any(|value| !value.is_finite()) {
        return Err("Q4 loads sum outside the floating-point range".into());
    }

    let capacity = params
        .elements
        .len()
        .checked_mul(64)
        .and_then(|n| n.checked_add(params.fixed_nodes.len().saturating_mul(2)))
        .ok_or("Q4 stiffness entry count overflow")?;
    let mut triplets = Vec::new();
    triplets
        .try_reserve_exact(capacity)
        .map_err(|e| format!("Cannot reserve Q4 assembly: {e}"))?;
    let blocks = crate::kernels::par::map_slice(&params.elements, 32, |element| {
        stiffness(
            element.map(|i| params.nodes[i]),
            params.youngs_modulus,
            params.poissons_ratio,
            true,
        )
    });
    for (element, block) in params.elements.iter().zip(blocks) {
        let block = block?;
        for (row, values) in block.iter().enumerate() {
            let node = element[row / 2];
            if fixed[node] {
                continue;
            }
            for (col, &value) in values.iter().enumerate() {
                let other = element[col / 2];
                if !fixed[other] {
                    triplets.push((node * 2 + row % 2, other * 2 + col % 2, value));
                }
            }
        }
    }
    for (node, &pinned) in fixed.iter().enumerate() {
        if pinned {
            let dof = node * 2;
            triplets.extend([(dof, dof, 1.0), (dof + 1, dof + 1, 1.0)]);
        }
    }

    let k_global: SparseMat = csr_from_triplets(n_dofs, n_dofs, &triplets);

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

/// Simplex-mesh elasticity (plane stress/strain, 3D, P1/P2): see [`crate::sim::elasticity_fem`].
pub use crate::sim::elasticity_fem::ElasticityProblem;
/// Simplex-mesh elasticity (plane stress/strain, 3D, P1/P2): see [`crate::sim::elasticity_fem`].
pub use crate::sim::elasticity_fem::Material as SimplexMaterial;
/// Simplex-mesh elasticity (plane stress/strain, 3D, P1/P2): see [`crate::sim::elasticity_fem`].
pub use crate::sim::elasticity_fem::solve as solve_simplex_elasticity;
