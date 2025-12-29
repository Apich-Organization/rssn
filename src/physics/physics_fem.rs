use ndarray::Array1;
use rayon::prelude::*;
use serde::Deserialize;
use serde::Serialize;

use crate::numerical::sparse::csr_from_triplets;
use crate::numerical::sparse::solve_conjugate_gradient;

#[allow(dead_code)]
#[derive(
    Clone, Debug, Serialize, Deserialize,
)]
/// A struct for Gaussian quadrature.

pub struct GaussQuadrature {
    points: Vec<f64>,
    weights: Vec<f64>,
}

impl GaussQuadrature {
    pub(crate) fn new() -> Self {

        let points = vec![
            -1.0 / 3.0_f64.sqrt(),
            1.0 / 3.0_f64.sqrt(),
        ];

        let weights = vec![1.0, 1.0];

        Self {
            points,
            weights,
        }
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

    let h = domain_length
        / n_elements as f64;

    let force_fn = &force_fn;

    // Parallel element assembly
    let element_data: Vec<(
        Vec<(usize, usize, f64)>,
        [f64; 2],
    )> = (0 .. n_elements)
        .into_par_iter()
        .map(move |i| {

            let x1 = i as f64 * h;

            let x2 = (i + 1) as f64 * h;

            let k_local = [
                [1.0 / h, -1.0 / h],
                [-1.0 / h, 1.0 / h],
            ];

            let f_local = [
                h / 2.0 * force_fn(x1),
                h / 2.0 * force_fn(x2),
            ];

            let nodes = [i, i + 1];

            let mut local_triplets =
                Vec::with_capacity(4);

            for r in 0 .. 2 {

                for c in 0 .. 2 {

                    local_triplets
                        .push((
                            nodes[r],
                            nodes[c],
                            k_local[r]
                                [c],
                        ));
                }
            }

            (
                local_triplets,
                f_local,
            )
        })
        .collect();

    let mut triplets =
        Vec::with_capacity(
            n_elements * 4,
        );

    let mut f = vec![0.0; n_nodes];

    for (i, (local_triplets, f_vals)) in
        element_data
            .into_iter()
            .enumerate()
    {

        triplets.extend(local_triplets);

        f[i] += f_vals[0];

        f[i + 1] += f_vals[1];
    }

    let last_node = n_nodes - 1;

    triplets.retain(|(r, _, _)| {

        *r != 0 && *r != last_node
    });

    triplets.push((0, 0, 1.0));

    triplets.push((
        last_node,
        last_node,
        1.0,
    ));

    f[0] = 0.0;

    f[last_node] = 0.0;

    let k_sparse = csr_from_triplets(
        n_nodes,
        n_nodes,
        &triplets,
    );

    let f_array = Array1::from(f);

    let u_array =
        solve_conjugate_gradient(
            &k_sparse,
            &f_array,
            None,
            1000,
            1e-9,
        )?;

    Ok(u_array.to_vec())
}

/// Example scenario for the 1D FEM Poisson solver.
///
/// # Errors
///
/// This function will return an error if the underlying `solve_poisson_1d` function
/// encounters an error.

pub fn simulate_1d_poisson_scenario(
) -> Result<Vec<f64>, String> {

    const N_ELEMENTS: usize = 50;

    const L: f64 = 1.0;

    let force = |_x: f64| 2.0;

    solve_poisson_1d(
        N_ELEMENTS,
        L,
        force,
    )
}

/// Solves the 2D Poisson equation on a unit square with zero Dirichlet boundaries.
///
/// # Errors
///
/// This function will return an error if the Conjugate Gradient solver fails to converge
/// or if the linear system is ill-conditioned.

pub fn solve_poisson_2d<F>(
    n_elements_x: usize,
    n_elements_y: usize,
    force_fn: F,
) -> Result<Vec<f64>, String>
where
    F: Fn(f64, f64) -> f64
        + Send
        + Sync,
{

    let (nx, ny) = (
        n_elements_x,
        n_elements_y,
    );

    let (n_nodes_x, n_nodes_y) =
        (nx + 1, ny + 1);

    let n_nodes = n_nodes_x * n_nodes_y;

    let (hx, hy) = (
        1.0 / nx as f64,
        1.0 / ny as f64,
    );

    let force_fn = &force_fn;

    let element_data : Vec<(
        Vec<(usize, usize, f64)>,
        [f64; 4],
        [usize; 4],
    )> = (0 .. ny)
        .into_par_iter()
        .flat_map(move |j| {

            (0 .. nx)
                .into_par_iter()
                .map(move |i| {

                    let gauss = GaussQuadrature::new();

                    let mut k_local = ndarray::Array2::<f64>::zeros((4, 4));

                    let mut f_local = [0.0; 4];

                    for gp_y in &gauss.points {

                        for gp_x in &gauss.points {

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

                            let det_j = (hx * hy) / 4.0;

                            let d_n_dx : Vec<f64> = d_n_dxi
                                .iter()
                                .map(|&d| d * 2.0 / hx)
                                .collect();

                            let d_n_dy : Vec<f64> = d_n_deta
                                .iter()
                                .map(|&d| d * 2.0 / hy)
                                .collect();

                            for r in 0 .. 4 {

                                for c in 0 .. 4 {

                                    k_local[[r, c]] +=
                                        d_n_dx[r].mul_add(d_n_dx[c], d_n_dy[r] * d_n_dy[c]) * det_j;
                                }
                            }

                            let x = (i as f64 + (1.0 + gp_x) / 2.0) * hx;

                            let y = (j as f64 + (1.0 + gp_y) / 2.0) * hy;

                            for k in 0 .. 4 {

                                f_local[k] += n[k] * force_fn(x, y) * det_j;
                            }
                        }
                    }

                    let nodes = [
                        j * n_nodes_x + i,
                        j * n_nodes_x + i + 1,
                        (j + 1) * n_nodes_x + i + 1,
                        (j + 1) * n_nodes_x + i,
                    ];

                    let mut local_triplets = Vec::with_capacity(16);

                    for r in 0 .. 4 {

                        for c in 0 .. 4 {

                            local_triplets.push((
                                nodes[r],
                                nodes[c],
                                k_local[[r, c]],
                            ));
                        }
                    }

                    (
                        local_triplets,
                        f_local,
                        nodes,
                    )
                })
        })
        .collect();

    let mut triplets =
        Vec::with_capacity(
            nx * ny * 16,
        );

    let mut f = vec![0.0; n_nodes];

    for (
        local_triplets,
        f_vals,
        nodes,
    ) in element_data
    {

        triplets.extend(local_triplets);

        for k in 0 .. 4 {

            f[nodes[k]] += f_vals[k];
        }
    }

    let mut boundary_nodes =
        std::collections::HashSet::new(
        );

    for j in 0 .. n_nodes_y {

        for i in 0 .. n_nodes_x {

            if i == 0
                || i == n_nodes_x - 1
                || j == 0
                || j == n_nodes_y - 1
            {

                boundary_nodes.insert(
                    j * n_nodes_x + i,
                );
            }
        }
    }

    triplets.retain(|(r, c, _)| {

        !boundary_nodes.contains(r)
            && !boundary_nodes
                .contains(c)
    });

    for node_idx in &boundary_nodes {

        triplets.push((
            *node_idx,
            *node_idx,
            1.0,
        ));

        f[*node_idx] = 0.0;
    }

    let k_sparse = csr_from_triplets(
        n_nodes,
        n_nodes,
        &triplets,
    );

    let f_array = Array1::from(f);

    let u_array =
        solve_conjugate_gradient(
            &k_sparse,
            &f_array,
            None,
            2000,
            1e-9,
        )?;

    Ok(u_array.to_vec())
}

/// Example scenario for the 2D FEM Poisson solver.
///
/// # Errors
///
/// This function will return an error if the underlying `solve_poisson_2d` function
/// encounters an error.
#[allow(clippy::unnecessary_cast)]

pub fn simulate_2d_poisson_scenario(
) -> Result<Vec<f64>, String> {

    const N_ELEMENTS: usize = 20;

    let force = |x, y| {

        2.0 * std::f64::consts::PI
            .powi(2)
            * (std::f64::consts::PI
                * (x as f64))
                .sin()
            * (std::f64::consts::PI
                * (y as f64))
                .sin()
    };

    solve_poisson_2d(
        N_ELEMENTS,
        N_ELEMENTS,
        force,
    )
}

/// Solves the 3D Poisson equation on a unit cube with zero Dirichlet boundaries.
///
/// # Errors
///
/// This function will return an error if the Conjugate Gradient solver fails to converge
/// or if the linear system is ill-conditioned.

pub fn solve_poisson_3d<F>(
    n_elements: usize,
    force_fn: F,
) -> Result<Vec<f64>, String>
where
    F: Fn(f64, f64, f64) -> f64
        + Send
        + Sync,
{

    let (nx, ny, nz) = (
        n_elements,
        n_elements,
        n_elements,
    );

    let (
        n_nodes_x,
        n_nodes_y,
        n_nodes_z,
    ) = (
        nx + 1,
        ny + 1,
        nz + 1,
    );

    let n_nodes = n_nodes_x
        * n_nodes_y
        * n_nodes_z;

    let (hx, hy, hz) = (
        1.0 / nx as f64,
        1.0 / ny as f64,
        1.0 / nz as f64,
    );

    let force_fn = &force_fn;

    let element_data : Vec<(
        Vec<(usize, usize, f64)>,
        [f64; 8],
        [usize; 8],
    )> = (0 .. nz)
        .into_par_iter()
        .flat_map(move |k_el| {

            (0 .. ny)
                .into_par_iter()
                .flat_map(move |j_el| {

                    (0 .. nx)
                        .into_par_iter()
                        .map(move |i_el| {

                            let mut k_local = ndarray::Array2::<f64>::zeros((8, 8));

                            let mut f_local = [0.0; 8];

                            let gauss = GaussQuadrature::new();

                            for gp_z in &gauss.points {

                                for gp_y in &gauss.points {

                                    for gp_x in &gauss.points {

                                        let mut n = [0.0; 8];

                                        let mut d_n_dxi = [0.0; 8];

                                        let mut d_n_deta = [0.0; 8];

                                        let mut d_n_dzeta = [0.0; 8];

                                        let xi = [-1.0, 1.0];

                                        for l in 0 .. 8 {

                                            let i = l & 1;

                                            let j = (l >> 1) & 1;

                                            let m = (l >> 2) & 1;

                                            n[l] = 0.125
                                                * (1.0 + xi[i] * gp_x)
                                                * (1.0 + xi[j] * gp_y)
                                                * (1.0 + xi[m] * gp_z);

                                            d_n_dxi[l] = 0.125
                                                * xi[i]
                                                * (1.0 + xi[j] * gp_y)
                                                * (1.0 + xi[m] * gp_z);

                                            d_n_deta[l] = 0.125
                                                * (1.0 + xi[i] * gp_x)
                                                * xi[j]
                                                * (1.0 + xi[m] * gp_z);

                                            d_n_dzeta[l] = 0.125
                                                * (1.0 + xi[i] * gp_x)
                                                * (1.0 + xi[j] * gp_y)
                                                * xi[m];
                                        }

                                        let det_j = (hx * hy * hz) / 8.0;

                                        let d_n_dx : Vec<f64> = d_n_dxi
                                            .iter()
                                            .map(|&d| d * 2.0 / hx)
                                            .collect();

                                        let d_n_dy : Vec<f64> = d_n_deta
                                            .iter()
                                            .map(|&d| d * 2.0 / hy)
                                            .collect();

                                        let d_n_dz : Vec<f64> = d_n_dzeta
                                            .iter()
                                            .map(|&d| d * 2.0 / hz)
                                            .collect();

                                        for r in 0 .. 8 {

                                            for c in 0 .. 8 {

                                                k_local[[r, c]] += d_n_dz[r].mul_add(d_n_dz[c], d_n_dx[r].mul_add(d_n_dx[c], d_n_dy[r] * d_n_dy[c]))
                                                    * det_j;
                                            }
                                        }

                                        let x = (i_el as f64 + (1.0 + gp_x) / 2.0) * hx;

                                        let y = (j_el as f64 + (1.0 + gp_y) / 2.0) * hy;

                                        let z = (k_el as f64 + (1.0 + gp_z) / 2.0) * hz;

                                        for l in 0 .. 8 {

                                            f_local[l] += n[l] * force_fn(x, y, z) * det_j;
                                        }
                                    }
                                }
                            }

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

                            let mut local_triplets = Vec::with_capacity(64);

                            for r in 0 .. 8 {

                                for c in 0 .. 8 {

                                    local_triplets.push((
                                        nodes[r],
                                        nodes[c],
                                        k_local[[r, c]],
                                    ));
                                }
                            }

                            (
                                local_triplets,
                                f_local,
                                nodes,
                            )
                        })
                })
        })
        .collect();

    let mut triplets =
        Vec::with_capacity(
            nx * ny * nz * 64,
        );

    let mut f = vec![0.0; n_nodes];

    for (
        local_triplets,
        f_vals,
        nodes,
    ) in element_data
    {

        triplets.extend(local_triplets);

        for l in 0 .. 8 {

            f[nodes[l]] += f_vals[l];
        }
    }

    let mut boundary_nodes =
        std::collections::HashSet::new(
        );

    for k in 0 .. n_nodes_z {

        for j in 0 .. n_nodes_y {

            for i in 0 .. n_nodes_x {

                let idx =
                    (k * n_nodes_y + j)
                        * n_nodes_x
                        + i;

                let is_boundary = i
                    == 0
                    || j == 0
                    || k == 0
                    || i == nx
                    || j == ny
                    || k == nz;

                if is_boundary {

                    boundary_nodes
                        .insert(idx);
                } else {

                    let x =
                        i as f64 * hx;

                    let y =
                        j as f64 * hy;

                    let z =
                        k as f64 * hz;

                    f[idx] = force_fn(
                        x, y, z,
                    );
                }
            }
        }
    }

    triplets.retain(|(r, c, _)| {

        !boundary_nodes.contains(r)
            && !boundary_nodes
                .contains(c)
    });

    for node_idx in &boundary_nodes {

        triplets.push((
            *node_idx,
            *node_idx,
            1.0,
        ));

        f[*node_idx] = 0.0;
    }

    let k_sparse = csr_from_triplets(
        n_nodes,
        n_nodes,
        &triplets,
    );

    let f_array = Array1::from(f);

    let u_array =
        solve_conjugate_gradient(
            &k_sparse,
            &f_array,
            None,
            3000,
            1e-9,
        )?;

    Ok(u_array.to_vec())
}

/// Example scenario for the 3D FEM Poisson solver.
///
/// # Errors
///
/// This function will return an error if the underlying `solve_poisson_3d` function
/// encounters an error.
#[allow(clippy::unnecessary_cast)]

pub fn simulate_3d_poisson_scenario(
) -> Result<Vec<f64>, String> {

    const N_ELEMENTS: usize = 10;

    let force = |x, y, z| {

        3.0 * std::f64::consts::PI
            .powi(2)
            * (std::f64::consts::PI
                * (x as f64))
                .sin()
            * (std::f64::consts::PI
                * (y as f64))
                .sin()
            * (std::f64::consts::PI
                * (z as f64))
                .sin()
    };

    solve_poisson_3d(N_ELEMENTS, force)
}
