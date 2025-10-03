// src/physics/physics_fem.rs

use crate::numerical::sparse::{csr_from_triplets, solve_conjugate_gradient};
use ndarray::Array1;

// --- Gaussian Quadrature ---
// Helper for numerical integration over elements.
#[allow(dead_code)]
pub struct GaussQuadrature {
    points: Vec<f64>,
    weights: Vec<f64>,
}

impl GaussQuadrature {
    // 2-point Gaussian quadrature is sufficient for linear elements.
    pub(crate) fn new() -> Self {
        let points = vec![-1.0 / 3.0_f64.sqrt(), 1.0 / 3.0_f64.sqrt()];
        let weights = vec![1.0, 1.0];
        GaussQuadrature { points, weights }
    }
}

// --- 1D FEM Implementation ---

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
pub fn solve_poisson_1d<F>(
    n_elements: usize,
    domain_length: f64,
    force_fn: F,
) -> Result<Vec<f64>, String>
where
    F: Fn(f64) -> f64,
{
    let n_nodes = n_elements + 1;
    let h = domain_length / n_elements as f64; // Element length

    let mut triplets = Vec::new();
    let mut f = vec![0.0; n_nodes];

    // Assemble stiffness matrix K and load vector F
    for i in 0..n_elements {
        let x1 = i as f64 * h;
        let x2 = (i + 1) as f64 * h;

        // Local stiffness matrix for a linear element
        let k_local = (1.0 / h) * ndarray::arr2(&[[1.0, -1.0], [-1.0, 1.0]]);

        // Local load vector (using trapezoidal rule for simplicity)
        let f_local = (h / 2.0) * ndarray::arr1(&[force_fn(x1), force_fn(x2)]);

        // Add to global matrix and vector
        let nodes = [i, i + 1];
        for r in 0..2 {
            f[nodes[r]] += f_local[r];
            for c in 0..2 {
                triplets.push((nodes[r], nodes[c], k_local[[r, c]]));
            }
        }
    }

    // Apply Dirichlet boundary conditions u(0)=0, u(L)=0
    // For node 0
    triplets.retain(|(r, _, _)| *r != 0);
    triplets.push((0, 0, 1.0));
    f[0] = 0.0;

    // For last node
    let last_node = n_nodes - 1;
    triplets.retain(|(r, _, _)| *r != last_node);
    triplets.push((last_node, last_node, 1.0));
    f[last_node] = 0.0;

    // Create sparse matrix and solve
    let k_sparse = csr_from_triplets(n_nodes, n_nodes, &triplets);
    let f_array = Array1::from(f);

    let u_array = solve_conjugate_gradient(&k_sparse, &f_array, None, 1000, 1e-9)?;
    Ok(u_array.to_vec())
}

/// Example scenario for the 1D FEM Poisson solver.
pub fn simulate_1d_poisson_scenario() -> Result<Vec<f64>, String> {
    const N_ELEMENTS: usize = 50;
    const L: f64 = 1.0;
    // f(x) = 2
    let force = |_x: f64| 2.0;
    // Exact solution is u(x) = x * (1-x)
    solve_poisson_1d(N_ELEMENTS, L, force)
}

// --- 2D FEM Implementation ---

/// Solves the 2D Poisson equation on a unit square with zero Dirichlet boundaries.
pub fn solve_poisson_2d<F>(
    n_elements_x: usize,
    n_elements_y: usize,
    force_fn: F,
) -> Result<Vec<f64>, String>
where
    F: Fn(f64, f64) -> f64,
{
    let (nx, ny) = (n_elements_x, n_elements_y);
    let (n_nodes_x, n_nodes_y) = (nx + 1, ny + 1);
    let n_nodes = n_nodes_x * n_nodes_y;
    let (hx, hy) = (1.0 / nx as f64, 1.0 / ny as f64);

    let mut triplets = Vec::new();
    let mut f = vec![0.0; n_nodes];
    let gauss = GaussQuadrature::new();

    // Loop over elements
    for j in 0..ny {
        for i in 0..nx {
            let mut k_local = ndarray::Array2::<f64>::zeros((4, 4));
            let mut f_local = ndarray::Array1::<f64>::zeros(4);

            // Gaussian quadrature
            for gp_y in &gauss.points {
                for gp_x in &gauss.points {
                    // Shape functions and their derivatives in natural coordinates (-1, 1)
                    let n = [
                        0.25 * (1.0 - gp_x) * (1.0 - gp_y),
                        0.25 * (1.0 + gp_x) * (1.0 - gp_y),
                        0.25 * (1.0 + gp_x) * (1.0 + gp_y),
                        0.25 * (1.0 - gp_x) * (1.0 + gp_y),
                    ];
                    let d_n_dxi = [
                        -0.25 * (1.0 - gp_y),
                        0.25 * (1.0 - gp_y),
                        0.25 * (1.0 + gp_y),
                        -0.25 * (1.0 + gp_y),
                    ];
                    let d_n_deta = [
                        -0.25 * (1.0 - gp_x),
                        -0.25 * (1.0 + gp_x),
                        0.25 * (1.0 + gp_x),
                        0.25 * (1.0 - gp_x),
                    ];

                    // Jacobian of the transformation
                    let j_inv = [[2.0 / hx, 0.0], [0.0, 2.0 / hy]];
                    let det_j = (hx * hy) / 4.0;

                    // Derivatives in physical coordinates
                    let mut d_n_dx = [0.0; 4];
                    let mut d_n_dy = [0.0; 4];
                    for k in 0..4 {
                        d_n_dx[k] = d_n_dxi[k] * j_inv[0][0];
                        d_n_dy[k] = d_n_deta[k] * j_inv[1][1];
                    }

                    // B matrix (gradient of shape functions)
                    let b = ndarray::arr2(&[d_n_dx, d_n_dy]);

                    // Add contribution to local stiffness matrix
                    k_local += &(&b.t().dot(&b) * det_j);

                    // Map gauss point to physical coords and add to load vector
                    let x = (i as f64 + (1.0 + gp_x) / 2.0) * hx;
                    let y = (j as f64 + (1.0 + gp_y) / 2.0) * hy;
                    for k in 0..4 {
                        f_local[k] += n[k] * force_fn(x, y) * det_j;
                    }
                }
            }

            // Map local to global
            let nodes = [
                j * n_nodes_x + i,
                (j * n_nodes_x) + i + 1,
                (j + 1) * n_nodes_x + i + 1,
                (j + 1) * n_nodes_x + i,
            ];
            for r in 0..4 {
                f[nodes[r]] += f_local[r];
                for c in 0..4 {
                    triplets.push((nodes[r], nodes[c], k_local[[r, c]]));
                }
            }
        }
    }

    // Apply boundary conditions
    let mut boundary_nodes = std::collections::HashSet::new();
    for j in 0..n_nodes_y {
        for i in 0..n_nodes_x {
            if i == 0 || i == n_nodes_x - 1 || j == 0 || j == n_nodes_y - 1 {
                boundary_nodes.insert(j * n_nodes_x + i);
            }
        }
    }

    triplets.retain(|(r, c, _)| !boundary_nodes.contains(r) && !boundary_nodes.contains(c));
    for node_idx in &boundary_nodes {
        triplets.push((*node_idx, *node_idx, 1.0));
        f[*node_idx] = 0.0;
    }

    // Solve
    let k_sparse = csr_from_triplets(n_nodes, n_nodes, &triplets);
    let f_array = Array1::from(f);
    let u_array = solve_conjugate_gradient(&k_sparse, &f_array, None, 2000, 1e-9)?;
    Ok(u_array.to_vec())
}

/// Example scenario for the 2D FEM Poisson solver.
pub fn simulate_2d_poisson_scenario() -> Result<Vec<f64>, String> {
    const N_ELEMENTS: usize = 20;
    // f(x,y) = 2*pi^2 * sin(pi*x) * sin(pi*y)
    let force = |x, y| {
        2.0 * std::f64::consts::PI.powi(2)
            * (std::f64::consts::PI * (x as f64)).sin()
            * (std::f64::consts::PI * (y as f64)).sin()
    };
    // Exact solution is u(x,y) = sin(pi*x) * sin(pi*y)
    solve_poisson_2d(N_ELEMENTS, N_ELEMENTS, force)
}

// --- 3D FEM Implementation ---
// Note: 3D FEM is significantly more complex. This is a simplified version on a structured grid.

/// Solves the 3D Poisson equation on a unit cube with zero Dirichlet boundaries.
pub fn solve_poisson_3d<F>(n_elements: usize, force_fn: F) -> Result<Vec<f64>, String>
where
    F: Fn(f64, f64, f64) -> f64,
{
    let (nx, ny, nz) = (n_elements, n_elements, n_elements);
    let (n_nodes_x, n_nodes_y, n_nodes_z) = (nx + 1, ny + 1, nz + 1);
    let n_nodes = n_nodes_x * n_nodes_y * n_nodes_z;
    let (hx, hy, hz) = (1.0 / nx as f64, 1.0 / ny as f64, 1.0 / nz as f64);

    let mut triplets = Vec::new();
    let mut f = vec![0.0; n_nodes];
    let gauss = GaussQuadrature::new();

    // Loop over elements
    for k_el in 0..nz {
        for j_el in 0..ny {
            for i_el in 0..nx {
                let mut k_local = ndarray::Array2::<f64>::zeros((8, 8));

                // Gaussian quadrature
                for gp_z in &gauss.points {
                    for gp_y in &gauss.points {
                        for gp_x in &gauss.points {
                            // Shape functions and derivatives for a hexahedral element
                            let mut n = [0.0; 8];
                            let mut d_n_dxi = [0.0; 8];
                            let mut d_n_deta = [0.0; 8];
                            let mut d_n_dzeta = [0.0; 8];
                            let xi = [-1.0, 1.0];
                            for l in 0..8 {
                                let i = l & 1;
                                let j = (l >> 1) & 1;
                                let m = (l >> 2) & 1;
                                n[l] = 0.125
                                    * (1.0 + xi[i] * gp_x)
                                    * (1.0 + xi[j] * gp_y)
                                    * (1.0 + xi[m] * gp_z);
                                d_n_dxi[l] =
                                    0.125 * xi[i] * (1.0 + xi[j] * gp_y) * (1.0 + xi[m] * gp_z);
                                d_n_deta[l] =
                                    0.125 * (1.0 + xi[i] * gp_x) * xi[j] * (1.0 + xi[m] * gp_z);
                                d_n_dzeta[l] =
                                    0.125 * (1.0 + xi[i] * gp_x) * (1.0 + xi[j] * gp_y) * xi[m];
                            }

                            let det_j = (hx * hy * hz) / 8.0;
                            let d_n_dx: Vec<f64> = d_n_dxi.iter().map(|&d| d * 2.0 / hx).collect();
                            let d_n_dy: Vec<f64> = d_n_deta.iter().map(|&d| d * 2.0 / hy).collect();
                            let d_n_dz: Vec<f64> =
                                d_n_dzeta.iter().map(|&d| d * 2.0 / hz).collect();

                            // Assemble local stiffness matrix
                            for r in 0..8 {
                                for c in 0..8 {
                                    let val = d_n_dx[r] * d_n_dx[c]
                                        + d_n_dy[r] * d_n_dy[c]
                                        + d_n_dz[r] * d_n_dz[c];
                                    k_local[[r, c]] += val * det_j;
                                }
                            }
                        }
                    }
                }

                // Map local to global
                let nodes = [
                    (k_el * n_nodes_y + j_el) * n_nodes_x + i_el,
                    (k_el * n_nodes_y + j_el) * n_nodes_x + i_el + 1,
                    (k_el * n_nodes_y + j_el + 1) * n_nodes_x + i_el + 1,
                    (k_el * n_nodes_y + j_el + 1) * n_nodes_x + i_el,
                    ((k_el + 1) * n_nodes_y + j_el) * n_nodes_x + i_el,
                    ((k_el + 1) * n_nodes_y + j_el) * n_nodes_x + i_el + 1,
                    ((k_el + 1) * n_nodes_y + j_el + 1) * n_nodes_x + i_el + 1,
                    ((k_el + 1) * n_nodes_y + j_el + 1) * n_nodes_x + i_el,
                ];
                for r in 0..8 {
                    for c in 0..8 {
                        triplets.push((nodes[r], nodes[c], k_local[[r, c]]));
                    }
                }
            }
        }
    }

    let mut boundary_nodes = std::collections::HashSet::new();

    // Simplified load vector (lumped) and boundary conditions
    for k in 0..n_nodes_z {
        for j in 0..n_nodes_y {
            for i in 0..n_nodes_x {
                let idx = (k * n_nodes_y + j) * n_nodes_x + i;
                let is_boundary = i == 0 || j == 0 || k == 0 || i == nx || j == ny || k == nz;
                if is_boundary {
                    boundary_nodes.insert(idx);
                } else {
                    let x = i as f64 * hx;
                    let y = j as f64 * hy;
                    let z = k as f64 * hz;
                    f[idx] = force_fn(x, y, z);
                }
            }
        }
    }

    // Apply boundary conditions
    triplets.retain(|(r, c, _)| !boundary_nodes.contains(r) && !boundary_nodes.contains(c));
    for node_idx in &boundary_nodes {
        triplets.push((*node_idx, *node_idx, 1.0));
        f[*node_idx] = 0.0;
    }

    // Solve
    let k_sparse = csr_from_triplets(n_nodes, n_nodes, &triplets);
    let f_array = Array1::from(f);
    let u_array = solve_conjugate_gradient(&k_sparse, &f_array, None, 3000, 1e-9)?;
    Ok(u_array.to_vec())
}

/// Example scenario for the 3D FEM Poisson solver.
pub fn simulate_3d_poisson_scenario() -> Result<Vec<f64>, String> {
    const N_ELEMENTS: usize = 10; // Keep it small, 3D is expensive
    let force = |x, y, z| {
        3.0 * std::f64::consts::PI.powi(2)
            * (std::f64::consts::PI * (x as f64)).sin()
            * (std::f64::consts::PI * (y as f64)).sin()
            * (std::f64::consts::PI * (z as f64)).sin()
    };
    // Exact solution is u(x,y,z) = sin(pi*x)*sin(pi*y)*sin(pi*z)
    solve_poisson_3d(N_ELEMENTS, force)
}
