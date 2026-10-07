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
use crate::kernels::krylov::cg;
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

/// Inflow velocity profile on the left boundary of a channel.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub enum InflowProfile {
    /// Uniform axial velocity `u = u0`.
    Uniform(f64),
    /// Poiseuille profile `u = 4 u_max y (1 - y)` on the unit-height channel.
    Parabolic(f64),
}

impl InflowProfile {
    /// Axial velocity at height `y` in `[0, 1]`.
    #[must_use]
    pub fn velocity(
        &self,
        y: f64,
    ) -> f64 {
        match *self {
            | Self::Uniform(u0) => u0,
            | Self::Parabolic(u_max) => 4.0 * u_max * y * (1.0 - y),
        }
    }
}

/// Condition on the top and bottom channel walls.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum WallCondition {
    /// Free slip: `v = 0` and `du/dy = 0`.
    Slip,
    /// No slip: `u = v = 0`.
    NoSlip,
}

/// Configuration of [`run_channel_flow_with`].
#[derive(Clone, Debug)]
pub struct ChannelFlowConfig<'a> {
    /// Grid points per direction (the domain is the unit square).
    pub n: usize,
    /// Reynolds number.
    pub re: f64,
    /// Time step.
    pub dt: f64,
    /// Number of time steps.
    pub n_iter: usize,
    /// Mask of solid grid points (`true` = no-slip obstacle), shape `(n, n)`.
    pub obstacle_mask: &'a Array2<bool>,
    /// Inflow profile on the left boundary.
    pub inflow: InflowProfile,
    /// Condition on the top and bottom walls.
    pub walls: WallCondition,
}

/// Solves the 2D Navier-Stokes equations for channel flow with an obstacle.
///
/// Free-slip walls and a uniform unit inflow; see [`run_channel_flow_with`]
/// for no-slip walls and a Poiseuille inflow, and for the boundary
/// treatment.
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
/// The solver works on a square `nx` x `nx` grid (`ny` must equal `nx`).
///
/// # Errors
/// Returns an error if the grid is smaller than 3 points, not square, or if
/// the obstacle mask does not have shape `(ny, nx)`; also if the pressure
/// solver fails.
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

    run_channel_flow_with(&ChannelFlowConfig {
        n: nx,
        re,
        dt,
        n_iter,
        obstacle_mask,
        inflow: InflowProfile::Uniform(1.0),
        walls: WallCondition::Slip,
    })
}

/// Applies the channel boundary conditions to a velocity field: prescribed
/// inflow on the left, zero-gradient outflow on the right, the wall
/// condition on top and bottom (walls win at the corners) and no-slip
/// (`u = v = 0`) on every obstacle point.
fn apply_channel_bcs(
    u: &mut Array2<f64>,
    v: &mut Array2<f64>,
    cfg: &ChannelFlowConfig<'_>,
) {
    let n = cfg.n;

    let h = 1.0 / (n as f64 - 1.0);

    for j in 0..n {
        u[[j, 0]] = cfg.inflow.velocity(j as f64 * h);

        v[[j, 0]] = 0.0;
    }

    for j in 0..n {
        u[[j, n - 1]] = u[[j, n - 2]];

        v[[j, n - 1]] = v[[j, n - 2]];
    }

    for i in 0..n {
        match cfg.walls {
            | WallCondition::Slip => {
                u[[0, i]] = u[[1, i]];

                u[[n - 1, i]] = u[[n - 2, i]];
            },
            | WallCondition::NoSlip => {
                u[[0, i]] = 0.0;

                u[[n - 1, i]] = 0.0;
            },
        }

        v[[0, i]] = 0.0;

        v[[n - 1, i]] = 0.0;
    }

    for ((j, i), &solid) in cfg.obstacle_mask.indexed_iter() {
        if solid {
            u[[j, i]] = 0.0;

            v[[j, i]] = 0.0;
        }
    }
}

/// Solves `-lap(p) = f` for the pressure on the interior nodes of the unit
/// square with conjugate gradients (Jacobi preconditioned).
///
/// Neumann (`dp/dn = 0`) conditions hold on the inflow, the walls and every
/// obstacle face, i.e. a neighbour that is outside the fluid is dropped from
/// the stencil; the outflow column `i = n - 1` is a Dirichlet condition
/// `p = 0`. `p` is the warm start and receives the solution; boundary
/// nodes are filled by copying the adjacent interior value (Neumann) or
/// zero (outflow).
fn solve_channel_pressure(
    p: &mut Array2<f64>,
    rhs: &Array2<f64>,
    solid: &Array2<bool>,
    h: f64,
) -> Result<(), String> {
    let n = p.nrows();

    let inv_h2 = 1.0 / (h * h);

    // Interior fluid unknowns, numbered row-major.
    let mut id = Array2::<usize>::from_elem((n, n), usize::MAX);

    let mut nodes = Vec::new();

    for j in 1..n - 1 {
        for i in 1..n - 1 {
            if !solid[[j, i]] {
                id[[j, i]] = nodes.len();

                nodes.push((j, i));
            }
        }
    }

    let m = nodes.len();

    if m == 0 {
        p.fill(0.0);

        return Ok(());
    }

    // Number of non-Neumann neighbours and the interior-fluid neighbours.
    let neighbours = |j: usize, i: usize| -> (f64, [Option<usize>; 4]) {
        let mut count = 0.0;

        let mut list = [None; 4];

        for (slot, (dj, di)) in [(0isize, -1isize), (0, 1), (-1, 0), (1, 0)].into_iter().enumerate() {
            let jj = j as isize + dj;

            let ii = i as isize + di;

            if ii == n as isize - 1 {
                // Dirichlet outflow column (p = 0): contributes only to the diagonal.
                count += 1.0;
            } else if (1..n as isize - 1).contains(&jj) && (1..n as isize - 1).contains(&ii) {
                let k = id[[jj as usize, ii as usize]];

                if k != usize::MAX {
                    count += 1.0;

                    list[slot] = Some(k);
                }
            }
        }

        (count, list)
    };

    let stencils: Vec<(f64, [Option<usize>; 4])> = nodes.iter().map(|&(j, i)| neighbours(j, i)).collect();

    let apply = |x: &[f64], y: &mut [f64]| {
        for (k, (count, list)) in stencils.iter().enumerate() {
            let mut s = count * x[k];

            for nb in list.iter().flatten() {
                s -= x[*nb];
            }

            y[k] = s * inv_h2;
        }
    };

    let diag: Vec<f64> = stencils.iter().map(|(c, _)| (c * inv_h2).max(f64::MIN_POSITIVE)).collect();

    let precond = |r: &[f64], z: &mut [f64]| {
        for k in 0..r.len() {
            z[k] = r[k] / diag[k];
        }
    };

    let b: Vec<f64> = nodes.iter().map(|&(j, i)| rhs[[j, i]]).collect();

    let x0: Vec<f64> = nodes.iter().map(|&(j, i)| p[[j, i]]).collect();

    let res = cg(apply, precond, &b, Some(&x0), 1e-10, 20 * m + 200);

    if !res.residual.is_finite() {
        return Err("Pressure solve produced a non-finite residual.".to_string());
    }

    for (k, &(j, i)) in nodes.iter().enumerate() {
        p[[j, i]] = res.x[k];
    }

    // Boundary and solid values: Neumann copy, Dirichlet outflow.
    for j in 0..n {
        p[[j, n - 1]] = 0.0;
    }

    for j in 0..n {
        p[[j, 0]] = p[[j.clamp(1, n - 2), 1]];
    }

    for i in 0..n - 1 {
        p[[0, i]] = p[[1, i]];

        p[[n - 1, i]] = p[[n - 2, i]];
    }

    for ((j, i), &s) in solid.indexed_iter() {
        if s {
            p[[j, i]] = 0.0;
        }
    }

    Ok(())
}

/// Solves the 2D incompressible Navier-Stokes equations in the unit square
/// with an optional no-slip obstacle, by Chorin's projection method.
///
/// Each step
/// 1. advances an intermediate velocity explicitly (first-order upwind
///    advection, central diffusion) at every fluid point; solid points and
///    boundary points are then set by the boundary conditions,
/// 2. imposes the boundary conditions: the inflow profile on the left,
///    zero-gradient outflow on the right, the [`WallCondition`] on top and
///    bottom, and *no slip* (`u = v = 0`) on every point of the obstacle
///    mask (so fluid points next to the obstacle see zero velocity in
///    their diffusion and advection stencils),
/// 3. solves the pressure Poisson equation `lap(p) = div(u*) / dt` over the
///    fluid points only, with `dp/dn = 0` on the inflow, the walls and all
///    obstacle faces and `p = 0` at the outflow, by preconditioned
///    conjugate gradients, and
/// 4. corrects `u = u* - dt grad(p)` at fluid points, where a gradient
///    across an obstacle face uses the Neumann (mirrored) pressure.
///
/// The explicit scheme needs `dt <= h^2 re / 4` and `dt <= h / |u|`.
///
/// # Errors
/// Returns an error if the grid has fewer than 3 points, the mask shape is
/// not `(n, n)`, or the pressure solve produces a non-finite residual.
pub fn run_channel_flow_with(cfg: &ChannelFlowConfig<'_>) -> NavierStokesOutput {
    let n = cfg.n;

    if n < 3 {
        return Err("Grid must have at least 3 points in each direction.".to_string());
    }

    if cfg.obstacle_mask.dim() != (n, n) {
        return Err(format!(
            "Obstacle mask has shape {:?} but the grid is ({n}, {n}).",
            cfg.obstacle_mask.dim()
        ));
    }

    let (re, dt) = (cfg.re, cfg.dt);

    let obstacle_mask = cfg.obstacle_mask;

    let h = 1.0 / (n as f64 - 1.0);

    let nu = 1.0 / re;

    let mut u = Array2::<f64>::zeros((n, n));

    let mut v = Array2::<f64>::zeros((n, n));

    let mut p = Array2::<f64>::zeros((n, n));

    // Start from the boundary data so the first step already sees the inflow.
    apply_channel_bcs(&mut u, &mut v, cfg);

    let mut rhs = Array2::<f64>::zeros((n, n));

    for _iter in 0..cfg.n_iter {
        let u_old = u.clone();

        let v_old = v.clone();

        // 1. Advection-diffusion at fluid points -> intermediate velocity.
        let mut u_star = u.clone();

        let mut v_star = v.clone();

        u_star
            .as_slice_mut()
            .ok_or_else(|| "velocity array is not contiguous".to_string())?
            .par_iter_mut()
            .zip(
                v_star
                    .as_slice_mut()
                    .ok_or_else(|| "velocity array is not contiguous".to_string())?
                    .par_iter_mut(),
            )
            .enumerate()
            .for_each(|(id, (u_val, v_val))| {
                let i = id % n;
                let j = id / n;

                // Boundary and obstacle points are set by the boundary
                // conditions below, not by the momentum equation.
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

                // Diffusion (Central Difference); solid neighbours hold zero velocity.
                let lap_u = (u_old[[j, i + 1]] + u_old[[j, i - 1]] + u_old[[j + 1, i]] + u_old[[j - 1, i]]
                    - 4.0 * u_curr)
                    / (h * h);

                let lap_v = (v_old[[j, i + 1]] + v_old[[j, i - 1]] + v_old[[j + 1, i]] + v_old[[j - 1, i]]
                    - 4.0 * v_curr)
                    / (h * h);

                *u_val = dt.mul_add(-u_curr.mul_add(du_dx, v_curr * du_dy) + nu * lap_u, u_curr);
                *v_val = dt.mul_add(-u_curr.mul_add(dv_dx, v_curr * dv_dy) + nu * lap_v, v_curr);
            });

        // 2. Boundary conditions on the intermediate velocity.
        apply_channel_bcs(&mut u_star, &mut v_star, cfg);

        // 3. Pressure Poisson over fluid points: -lap(p) = -div(u*) / dt.
        rhs.fill(0.0);

        for j in 1..n - 1 {
            for i in 1..n - 1 {
                if obstacle_mask[[j, i]] {
                    continue;
                }

                let div = (u_star[[j, i + 1]] - u_star[[j, i - 1]]) / (2.0 * h)
                    + (v_star[[j + 1, i]] - v_star[[j - 1, i]]) / (2.0 * h);

                rhs[[j, i]] = -div / dt;
            }
        }

        solve_channel_pressure(&mut p, &rhs, obstacle_mask, h)?;

        // 4. Velocity correction with Neumann-mirrored pressure at solid faces.
        let pressure_at = |j: usize, i: usize, jc: usize, ic: usize| -> f64 {
            if obstacle_mask[[j, i]] { p[[jc, ic]] } else { p[[j, i]] }
        };

        for j in 1..n - 1 {
            for i in 1..n - 1 {
                if obstacle_mask[[j, i]] {
                    u[[j, i]] = 0.0;

                    v[[j, i]] = 0.0;

                    continue;
                }

                let dp_dx = (pressure_at(j, i + 1, j, i) - pressure_at(j, i - 1, j, i)) / (2.0 * h);

                let dp_dy = (pressure_at(j + 1, i, j, i) - pressure_at(j - 1, i, j, i)) / (2.0 * h);

                u[[j, i]] = dt.mul_add(-dp_dx, u_star[[j, i]]);

                v[[j, i]] = dt.mul_add(-dp_dy, v_star[[j, i]]);
            }
        }

        apply_channel_bcs(&mut u, &mut v, cfg);
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
