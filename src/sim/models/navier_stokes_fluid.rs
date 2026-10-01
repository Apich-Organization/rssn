//! # Navier-Stokes Fluid Simulation
//!
//! This module provides Computational Fluid Dynamics (CFD) solvers for incompressible, viscous fluid flow
//! governed by the Navier-Stokes equations. It focuses on 2D simulations using standard numerical techniques
//! for pressure-velocity coupling.
//!
//! # Overview
//!
//! The simulation engine primarily uses the projection method (Chorin, 1968) to decouple the velocity
//! and pressure fields. This involves an intermediate velocity prediction followed by a pressure correction step
//! to enforce the incompressibility constraint (continuity equation).
//!
//! Key features include:
//! - **Stable Solvers**: Implements the projection method for stable time-stepping.
//! - **Multigrid Acceleration**: Utilizes a geometric multigrid solver (V-cycles) for the pressure Poisson equation, critical for performance on fine grids.
//! - **Parallel Computation**: Leverages `rayon` for parallelized grid operations (advection, diffusion updates).
//! - **Diverse Scenarios**:
//!     - **Channel Flow**: Simulates flow past obstacles with configurable Reynolds numbers.
//!     - **Lid-Driven Cavity**: A classic CFD benchmark problem for testing internal flow and vortex formation.
//!
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/karman_velocity_mag.png)

use std::path::Path;

use ndarray::Array2;
use rayon::prelude::*;
use serde::Deserialize;
use serde::Serialize;

use crate::io::write_npy_file;
use crate::sim::physics_mtm::solve_poisson_2d_multigrid;

/// Parameters for the Navier-Stokes simulation.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct NavierStokesParameters {
    /// Number of grid points in the x-direction.
    pub nx: usize,
    /// Number of grid points in the y-direction.
    pub ny: usize,
    /// Reynolds number.
    pub re: f64,
    /// Time step size.
    pub dt: f64,
    /// Number of simulation iterations.
    pub n_iter: usize,
    /// Velocity of the lid (driving force).
    pub lid_velocity: f64,
}

/// Type of `NavierStokesOutput`.
pub type NavierStokesOutput = Result<(Array2<f64>, Array2<f64>, Array2<f64>), String>;

/// Solves the 2D Navier-Stokes equations for channel flow with an obstacle.
/// Uses the projection method (Chorin 1968).
///
/// # Arguments
/// * `nx` - Grid width.
/// * `ny` - Grid height.
/// * `re` - Reynolds number.
/// * `dt` - Time step.
/// * `n_iter` - Number of iterations.
/// * `obstacle_mask` - A boolean mask where true indicates an obstacle (u=0).
///
/// # Returns
/// Tuple of (u, v, p) arrays.
///
/// The solver works on a square `nx` x `nx` grid (`ny` must equal `nx`) whose
/// size is `2^k + 1`, as required by the multigrid pressure solver.
///
/// # Errors
/// Returns an error if the grid is smaller than 3 points, not square, or if
/// the obstacle mask does not have shape `(ny, nx)`; also if the multigrid
/// solver rejects the grid size.
///
/// # Panics
/// Panics if an intermediate array is not contiguous, which cannot happen for
/// arrays this function allocates itself.
pub fn run_channel_flow(
    nx: usize,
    ny: usize,
    re: f64,
    dt: f64,
    n_iter: usize,
    obstacle_mask: &Array2<bool>,
) -> NavierStokesOutput {
    if nx < 3 || ny < 3 {
        return Err("Grid must have at least 3 points in each direction.".to_string());
    }

    if nx != ny {
        return Err(format!(
            "Channel flow requires a square grid, got {nx} x {ny}."
        ));
    }

    if obstacle_mask.dim() != (ny, nx) {
        return Err(format!(
            "Obstacle mask has shape {:?} but the grid is ({ny}, {nx}).",
            obstacle_mask.dim()
        ));
    }

    let n = nx;

    let h = 1.0 / (n as f64 - 1.0);

    let nu = 1.0 / re;

    let mut u = Array2::<f64>::zeros((n, n));

    let mut v = Array2::<f64>::zeros((n, n));

    let mut p = Array2::<f64>::zeros((n, n));

    // Helper for indexing
    let idx = |i: usize, j: usize| j * n + i;

    // Multigrid size
    let mg_size = n;

    for _iter in 0..n_iter {
        let u_old = u.clone();

        let v_old = v.clone();

        // 1. Advection-Diffusion (Explicit) -> Intermediate Velocity (u*, v*)
        // u_star = u_old + dt * ( - (u grad) u + nu laplacian u )
        let mut u_star = u.clone();

        let mut v_star = v.clone();

        u_star
            .as_slice_mut()
            .expect("Contiguous array")
            .par_iter_mut()
            .zip(
                v_star
                    .as_slice_mut()
                    .expect("Contiguous array")
                    .par_iter_mut(),
            )
            .enumerate()
            .for_each(|(id, (u_val, v_val))| {
                let i = id % n;
                let j = id / n;

                // Skip boundaries and obstacles for now
                if i == 0 || i == n - 1 || j == 0 || j == n - 1 || obstacle_mask[[j, i]] {
                    return;
                }

                // Upwind Advection
                let u_curr = u_old[[j, i]];
                let v_curr = v_old[[j, i]];

                let du_dx = if u_curr > 0.0 {
                    u_curr - u_old[[j, i - 1]]
                } else {
                    u_old[[j, i + 1]] - u_curr
                } / h;

                let du_dy = if v_curr > 0.0 {
                    u_curr - u_old[[j - 1, i]]
                } else {
                    u_old[[j + 1, i]] - u_curr
                } / h;

                let dv_dx = if u_curr > 0.0 {
                    v_curr - v_old[[j, i - 1]]
                } else {
                    v_old[[j, i + 1]] - v_curr
                } / h;

                let dv_dy = if v_curr > 0.0 {
                    v_curr - v_old[[j - 1, i]]
                } else {
                    v_old[[j + 1, i]] - v_curr
                } / h;

                // Diffusion (Central Difference)
                let lap_u = (2.0f64.mul_add(
                    -u_curr,
                    2.0f64.mul_add(-u_curr, u_old[[j, i + 1]])
                        + u_old[[j, i - 1]]
                        + u_old[[j + 1, i]],
                ) + u_old[[j - 1, i]])
                    / (h * h);

                let lap_v = (2.0f64.mul_add(
                    -v_curr,
                    2.0f64.mul_add(-v_curr, v_old[[j, i + 1]])
                        + v_old[[j, i - 1]]
                        + v_old[[j + 1, i]],
                ) + v_old[[j - 1, i]])
                    / (h * h);

                *u_val = dt.mul_add(-u_curr.mul_add(du_dx, v_curr * du_dy) + nu * lap_u, u_curr);
                *v_val = dt.mul_add(-u_curr.mul_add(dv_dx, v_curr * dv_dy) + nu * lap_v, v_curr);
            });

        // Apply BCs to u_star, v_star
        // Inflow (Left)
        for j in 0..n {
            u_star[[j, 0]] = 1.0;

            v_star[[j, 0]] = 0.0;
        }

        // Outflow (Right) - Zero Gradient
        for j in 0..n {
            u_star[[j, n - 1]] = u_star[[j, n - 2]];

            v_star[[j, n - 1]] = v_star[[j, n - 2]];
        }

        // Walls (Top/Bottom) - We do Slip for channel here
        for i in 0..n {
            u_star[[0, i]] = u_star[[1, i]]; // Slip
            v_star[[0, i]] = 0.0;

            u_star[[n - 1, i]] = u_star[[n - 2, i]]; // Slip
            v_star[[n - 1, i]] = 0.0;
        }

        // Obstacle - No Slip
        for j in 0..n {
            for i in 0..n {
                if obstacle_mask[[j, i]] {
                    u_star[[j, i]] = 0.0;

                    v_star[[j, i]] = 0.0;
                }
            }
        }

        // 2. Pressure Correction (Poisson Step)
        // laplacian p = div(u*) / dt
        let mut rhs_vec = vec![0.0; n * n];

        // Calculate Div U*
        for j in 1..n - 1 {
            for i in 1..n - 1 {
                if obstacle_mask[[j, i]] {
                    continue;
                }

                let div = (u_star[[j, i + 1]] - u_star[[j, i - 1]]) / (2.0 * h)
                    + (v_star[[j + 1, i]] - v_star[[j - 1, i]]) / (2.0 * h);

                rhs_vec[idx(i, j)] = div / dt;
            }
        }

        // Solve Poisson
        // Note: Our MG solver expects "f" where "-laplacian u = f".
        // Our equation is "laplacian p = R". So we pass "-R" as f.
        for val in &mut rhs_vec {
            *val = -*val;
        }

        let p_flat = solve_poisson_2d_multigrid(mg_size, &rhs_vec, 2)?; // 2 V-cycles

        // Copy back to P array
        for j in 0..n {
            for i in 0..n {
                p[[j, i]] = p_flat[idx(i, j)];
            }
        }

        // 3. Velocity Correction
        // u_new = u* - dt * grad p
        for j in 1..n - 1 {
            for i in 1..n - 1 {
                if obstacle_mask[[j, i]] {
                    u[[j, i]] = 0.0;

                    v[[j, i]] = 0.0;

                    continue;
                }

                let dp_dx = (p[[j, i + 1]] - p[[j, i - 1]]) / (2.0 * h);

                let dp_dy = (p[[j + 1, i]] - p[[j - 1, i]]) / (2.0 * h);

                u[[j, i]] = dt.mul_add(-dp_dx, u_star[[j, i]]);

                v[[j, i]] = dt.mul_add(-dp_dy, v_star[[j, i]]);
            }
        }

        // Re-apply BCs for u, v
        // Inflow (Left)
        for j in 0..n {
            u[[j, 0]] = 1.0;

            v[[j, 0]] = 0.0;
        }

        // Outflow (Right)
        for j in 0..n {
            u[[j, n - 1]] = u[[j, n - 2]];

            v[[j, n - 1]] = v[[j, n - 2]];
        }

        // Walls (Top/Bottom)
        for i in 0..n {
            u[[0, i]] = u[[1, i]];

            v[[0, i]] = 0.0;

            u[[n - 1, i]] = u[[n - 2, i]];

            v[[n - 1, i]] = 0.0;
        }
    }

    Ok((u, v, p))
}

/// Main solver for the 2D lid-driven cavity problem.
///
/// Chorin's projection method on a staggered (MAC) grid: `p[j][i]` lives at
/// cell centres, `u[j][i]` on the left face of cell `(j, i)` and `v[j][i]` on
/// its bottom face. Each step
/// 1. advances the intermediate velocity explicitly with first-order upwind
///    advection and central diffusion with `nu = 1 / re`,
/// 2. solves `lap(phi) = div(u*) / dt` with the multigrid solver (the same
///    solver and zero-Dirichlet closure as [`run_channel_flow`]), and
/// 3. corrects `u = u* - dt grad(phi)` and sets `p = phi`.
///
/// The top row of `u` is the moving lid (`lid_velocity`), the bottom row and
/// the side faces are at rest. The explicit scheme needs
/// `dt <= h^2 re / 4` (diffusion) and `dt <= h / |u|` (advection) to stay stable.
///
/// # Errors
///
/// Returns an error if the grid is smaller than 3 points, if the cells are
/// not square (`nx != ny`), or if the multigrid solver fails.
pub fn run_lid_driven_cavity(params: &NavierStokesParameters) -> NavierStokesOutput {
    let (nx, ny, re, dt) = (params.nx, params.ny, params.re, params.dt);

    if nx < 3 || ny < 3 {
        return Err("Grid must have at least 3 points in each direction.".to_string());
    }

    if nx != ny {
        return Err(format!(
            "The cavity solver requires square cells, got {nx} x {ny}."
        ));
    }

    let hx = 1.0 / (nx - 1) as f64;

    let hy = 1.0 / (ny - 1) as f64;

    let nu = 1.0 / re;

    let mut u = Array2::<f64>::zeros((ny, nx + 1));

    let mut v = Array2::<f64>::zeros((ny + 1, nx));

    let mut p = Array2::<f64>::zeros((ny, nx));

    // Boundary conditions: lid velocity
    for j in 0..=nx {
        u[[ny - 1, j]] = params.lid_velocity;
    }

    let mg_size_k = (((nx.max(ny) - 1) as f64).log2().ceil() as i64)
        .try_into()
        .unwrap_or(0);

    let mg_size = 2_usize.pow(mg_size_k) + 1;

    // The multigrid solver assumes spacing 1 / (mg_size - 1); rescale the
    // right-hand side so that it solves the Poisson problem for spacing hx.
    let h_mg = 1.0 / (mg_size - 1) as f64;

    let rhs_scale = (hx / h_mg).powi(2);

    for _ in 0..params.n_iter {
        // 1. Intermediate velocity (explicit advection + diffusion).
        let mut u_star = u.clone();

        let mut v_star = v.clone();

        for j in 1..ny - 1 {
            for i in 1..nx {
                let uc = u[[j, i]];

                // v averaged to the u-face.
                let vc = 0.25 * (v[[j, i - 1]] + v[[j, i]] + v[[j + 1, i - 1]] + v[[j + 1, i]]);

                let du_dx = if uc > 0.0 {
                    uc - u[[j, i - 1]]
                } else {
                    u[[j, i + 1]] - uc
                } / hx;

                let du_dy = if vc > 0.0 {
                    uc - u[[j - 1, i]]
                } else {
                    u[[j + 1, i]] - uc
                } / hy;

                let lap = (u[[j, i + 1]] - 2.0 * uc + u[[j, i - 1]]) / (hx * hx)
                    + (u[[j + 1, i]] - 2.0 * uc + u[[j - 1, i]]) / (hy * hy);

                u_star[[j, i]] = uc + dt * (-(uc * du_dx + vc * du_dy) + nu * lap);
            }
        }

        for j in 1..ny {
            for i in 1..nx - 1 {
                let vc = v[[j, i]];

                // u averaged to the v-face.
                let uc = 0.25 * (u[[j - 1, i]] + u[[j - 1, i + 1]] + u[[j, i]] + u[[j, i + 1]]);

                let dv_dx = if uc > 0.0 {
                    vc - v[[j, i - 1]]
                } else {
                    v[[j, i + 1]] - vc
                } / hx;

                let dv_dy = if vc > 0.0 {
                    vc - v[[j - 1, i]]
                } else {
                    v[[j + 1, i]] - vc
                } / hy;

                let lap = (v[[j, i + 1]] - 2.0 * vc + v[[j, i - 1]]) / (hx * hx)
                    + (v[[j + 1, i]] - 2.0 * vc + v[[j - 1, i]]) / (hy * hy);

                v_star[[j, i]] = vc + dt * (-(uc * dv_dx + vc * dv_dy) + nu * lap);
            }
        }

        // 2. Pressure Poisson: lap(phi) = div(u*) / dt, i.e. -lap(phi) = -div/dt.
        let mut rhs_padded = vec![0.0; mg_size * mg_size];

        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let div_u_star = (u_star[[j, i + 1]] - u_star[[j, i]]) / hx
                    + (v_star[[j + 1, i]] - v_star[[j, i]]) / hy;

                rhs_padded[j * mg_size + i] = -rhs_scale * div_u_star / dt;
            }
        }

        let phi_vec = solve_poisson_2d_multigrid(mg_size, &rhs_padded, 10)?;

        let phi = Array2::from_shape_vec((mg_size, mg_size), phi_vec).map_err(|e| e.to_string())?;

        // 3. Projection: u = u* - dt grad(phi); the pressure is phi itself.
        for j in 0..ny {
            for i in 0..nx {
                p[[j, i]] = phi[[j, i]];
            }
        }

        u = u_star;

        v = v_star;

        for j in 1..ny - 1 {
            for i in 1..nx {
                u[[j, i]] -= dt / hx * (phi[[j, i]] - phi[[j, i - 1]]);
            }
        }

        for j in 1..ny {
            for i in 1..nx - 1 {
                v[[j, i]] -= dt / hy * (phi[[j, i]] - phi[[j - 1, i]]);
            }
        }
    }

    // Interpolate the staggered velocities to cell centres.
    let mut u_centered = Array2::<f64>::zeros((ny, nx));

    let mut v_centered = Array2::<f64>::zeros((ny, nx));

    for j in 0..ny {
        for i in 0..nx {
            u_centered[[j, i]] = f64::midpoint(u[[j, i]], u[[j, i + 1]]);

            v_centered[[j, i]] = f64::midpoint(v[[j, i]], v[[j + 1, i]]);
        }
    }

    Ok((u_centered, v_centered, p))
}

/// An example scenario for the lid-driven cavity simulation.
///
/// The fields are written as `cavity_u_velocity.npy`, `cavity_v_velocity.npy`
/// and `cavity_pressure.npy` inside `output_dir`, which is created if missing.
pub fn simulate_lid_driven_cavity_scenario(output_dir: &Path) {
    const K: usize = 6;

    const N: usize = 2_usize.pow(K as u32) + 1;

    println!(
        "Running 2D Lid-Driven Cavity \
         simulation..."
    );

    let params = NavierStokesParameters {
        nx: N,
        ny: N,
        re: 100.0,
        dt: 0.01,
        n_iter: 200,
        lid_velocity: 1.0,
    };

    match run_lid_driven_cavity(&params) {
        | Ok((u, v, p)) => {
            println!(
                "Simulation finished. \
                 Saving results..."
            );

            let save_result = (|| -> Result<(), String> {
                std::fs::create_dir_all(output_dir).map_err(|e| e.to_string())?;

                write_npy_file(output_dir.join("cavity_u_velocity.npy"), &u)?;

                write_npy_file(output_dir.join("cavity_v_velocity.npy"), &v)?;

                write_npy_file(output_dir.join("cavity_pressure.npy"), &p)?;

                Ok(())
            })();

            if let Err(e) = save_result {
                eprintln!(
                    "Failed to save \
                     results: {e}"
                );
            } else {
                println!(
                    "Results saved to \
                     .npy files."
                );
            }
        },
        | Err(e) => {
            eprintln!(
                "An error occurred \
                 during simulation: \
                 {e}"
            );
        },
    }
}
