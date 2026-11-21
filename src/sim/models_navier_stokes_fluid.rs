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
//!
//! # Coverage audit
//!
//! | Case | Dimension | Boundaries | Notes |
//! |---|---|---|---|
//! | channel flow with obstacles | 2D | inflow profile, outflow, walls | first-order projection |
//! | lid-driven cavity | 2D | moving lid, no-slip | |
//!
//! [`crate::sim::ns_mac`] is the second-order MAC solver for 2D and 3D with
//! no-slip, slip, inflow, outflow and periodic faces, body forces and obstacles.

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

const PRESSURE_TOLERANCE: f64 = 1e-10;

struct PressureRow {
    cell: usize,
    neighbors: [usize; 4],
    degree: u8,
    diagonal: u8,
}

impl PressureRow {
    fn adjacent(&self) -> &[usize] {
        &self.neighbors[..usize::from(self.degree)]
    }
}

/// Prepared pressure graph and per-step vectors for one fixed obstacle mask.
struct ChannelPressure {
    n: usize,
    h: f64,
    inv_h2: f64,
    rows: Vec<PressureRow>,
    diagonal: Vec<f64>,
    closed: Vec<usize>,
    closed_ends: Vec<usize>,
    rhs: Vec<f64>,
    initial: Vec<f64>,
}

impl ChannelPressure {
    fn new(
        solid: &Array2<bool>,
        h: f64,
    ) -> Self {
        let n = solid.nrows();
        let mut id = Array2::from_elem((n, n), usize::MAX);
        let mut rows = Vec::new();
        for j in 1..n - 1 {
            for i in 1..n - 1 {
                if !solid[[j, i]] {
                    id[[j, i]] = rows.len();
                    rows.push(PressureRow {
                        cell: j * n + i,
                        neighbors: [0; 4],
                        degree: 0,
                        diagonal: 0,
                    });
                }
            }
        }
        for row in &mut rows {
            let (j, i) = (row.cell / n, row.cell % n);
            for (jj, ii) in [(j, i - 1), (j, i + 1), (j - 1, i), (j + 1, i)] {
                let k = id[[jj, ii]];
                if k != usize::MAX {
                    row.neighbors[usize::from(row.degree)] = k;
                    row.degree += 1;
                }
            }
            // A solid point on the outflow is a wall face, not a pressure anchor.
            row.diagonal = row.degree + u8::from(i == n - 2 && !solid[[j, n - 1]]);
        }

        let width = n - 2;
        let connected_rectangle =
            rows.len() == width * width && (1..n - 1).any(|j| !solid[[j, n - 1]]);
        let (closed, closed_ends) = if connected_rectangle {
            (Vec::new(), Vec::new())
        } else {
            Self::pin_neumann_components(&mut rows)
        };
        let inv_h2 = 1.0 / (h * h);
        let diagonal = rows
            .iter()
            .map(|row| f64::from(row.diagonal) * inv_h2)
            .collect();
        let count = rows.len();
        Self {
            n,
            h,
            inv_h2,
            rows,
            diagonal,
            closed,
            closed_ends,
            rhs: vec![0.0; count],
            initial: vec![0.0; count],
        }
    }

    fn pin_neumann_components(rows: &mut [PressureRow]) -> (Vec<usize>, Vec<usize>) {
        let mut visited = vec![false; rows.len()];
        let mut pins = vec![false; rows.len()];
        let mut component = Vec::new();
        let mut closed = Vec::new();
        let mut closed_ends = Vec::new();
        for seed in 0..rows.len() {
            if visited[seed] {
                continue;
            }
            visited[seed] = true;
            component.clear();
            component.push(seed);
            let mut offset = 0;
            let mut has_outflow = false;
            while let Some(&k) = component.get(offset) {
                let row = &rows[k];
                has_outflow |= row.diagonal > row.degree;
                for &neighbor in row.adjacent() {
                    if !visited[neighbor] {
                        visited[neighbor] = true;
                        component.push(neighbor);
                    }
                }
                offset += 1;
            }
            if !has_outflow {
                pins[seed] = true;
                closed.extend_from_slice(&component);
                closed_ends.push(closed.len());
            }
        }

        if !closed.is_empty() {
            for (i, row) in rows.iter_mut().enumerate() {
                if pins[i] {
                    row.degree = 0;
                    row.diagonal = 1;
                } else {
                    // A pinned neighbor is zero Dirichlet: keep its diagonal
                    // contribution and drop its off-diagonal coefficient.
                    let mut degree = 0;
                    for j in 0..usize::from(row.degree) {
                        let k = row.neighbors[j];
                        if !pins[k] {
                            row.neighbors[usize::from(degree)] = k;
                            degree += 1;
                        }
                    }
                    row.degree = degree;
                }
            }
        }
        (closed, closed_ends)
    }

    fn solve(
        &mut self,
        pressure: &mut Array2<f64>,
        u: &[f64],
        v: &[f64],
        solid: &Array2<bool>,
        dt: f64,
    ) -> Result<(), String> {
        let n = self.n;
        let h = 2.0 * self.h;
        let p = pressure
            .as_slice_mut()
            .ok_or("channel pressure storage is not contiguous")?;
        for (k, row) in self.rows.iter().enumerate() {
            let i = row.cell;
            let div = (u[i + 1] - u[i - 1]) / h + (v[i + n] - v[i - n]) / h;
            let force = -div / dt;
            if !force.is_finite() {
                return Err("channel pressure source is non-finite".into());
            }
            self.rhs[k] = force;
            self.initial[k] = p[i];
        }
        self.check_closed_flux()?;
        if self.rhs.iter().all(|&value| value == 0.0) {
            pressure.fill(0.0);
            return Ok(());
        }

        let apply = |x: &[f64], y: &mut [f64]| {
            for (k, row) in self.rows.iter().enumerate() {
                let mut value = f64::from(row.diagonal) * x[k];
                for &neighbor in row.adjacent() {
                    value -= x[neighbor];
                }
                y[k] = value * self.inv_h2;
            }
        };
        let precond = |r: &[f64], z: &mut [f64]| {
            for ((z, r), d) in z.iter_mut().zip(r).zip(&self.diagonal) {
                *z = r / d;
            }
        };
        let cap = self.rows.len().saturating_mul(20).saturating_add(200);
        let result = cg(
            apply,
            precond,
            &self.rhs,
            Some(&self.initial),
            PRESSURE_TOLERANCE,
            cap,
        );
        if !result.converged {
            return Err(format!(
                "channel pressure solve did not converge after {} iterations (relative residual {})",
                result.iterations, result.residual
            ));
        }
        for (row, &value) in self.rows.iter().zip(&result.x) {
            p[row.cell] = value;
        }
        Self::fill_boundaries(pressure, solid);
        Ok(())
    }

    fn check_closed_flux(&mut self) -> Result<(), String> {
        let mut start = 0;
        for &end in &self.closed_ends {
            let component = &self.closed[start..end];
            let scale = component
                .iter()
                .map(|&k| self.rhs[k].abs())
                .fold(0.0, f64::max);
            if scale > 0.0 {
                let (sum, squared) = component.iter().fold((0.0, 0.0_f64), |(sum, squared), &k| {
                    let term = self.rhs[k] / scale;
                    (sum + term, squared + term * term)
                });
                if sum.abs() > PRESSURE_TOLERANCE * squared.sqrt() {
                    return Err("closed channel region has incompatible net flux".into());
                }
            }
            let gauge = component[0];
            self.rhs[gauge] = 0.0;
            self.initial[gauge] = 0.0;
            start = end;
        }
        Ok(())
    }

    fn fill_boundaries(
        pressure: &mut Array2<f64>,
        solid: &Array2<bool>,
    ) {
        let n = pressure.nrows();
        for j in 0..n {
            pressure[[j, n - 1]] = 0.0;
            pressure[[j, 0]] = pressure[[j.clamp(1, n - 2), 1]];
        }
        for i in 0..n - 1 {
            pressure[[0, i]] = pressure[[1, i]];
            pressure[[n - 1, i]] = pressure[[n - 2, i]];
        }
        for ((j, i), &obstacle) in solid.indexed_iter() {
            if obstacle {
                pressure[[j, i]] = 0.0;
            }
        }
    }
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
/// Returns an error for invalid grid size, mask shape, time step, Reynolds
/// number or inflow, non-finite fields, incompatible flux in a closed
/// fluid region, or an unconverged pressure solve.
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

    if !re.is_finite() || re <= 0.0 || !nu.is_finite() {
        return Err("channel Reynolds number must give finite positive viscosity".into());
    }
    if !dt.is_finite() || dt < 0.0 || (cfg.n_iter != 0 && dt == 0.0) {
        return Err("channel time step must be finite and positive when advancing".into());
    }
    let inflow_speed = match cfg.inflow {
        | InflowProfile::Uniform(speed) | InflowProfile::Parabolic(speed) => speed,
    };
    if !inflow_speed.is_finite() {
        return Err("channel inflow speed must be finite".into());
    }

    let mut u = Array2::<f64>::zeros((n, n));

    let mut v = Array2::<f64>::zeros((n, n));

    let mut p = Array2::<f64>::zeros((n, n));

    // Start from the boundary data so the first step already sees the inflow.
    apply_channel_bcs(&mut u, &mut v, cfg);

    if u.iter().any(|x| !x.is_finite()) {
        return Err("channel inflow profile must be finite on the grid".into());
    }
    if cfg.n_iter == 0 {
        return Ok((u, v, p));
    }

    let mut pressure = ChannelPressure::new(obstacle_mask, h);
    let mut u_star = Array2::<f64>::zeros((n, n));
    let mut v_star = Array2::<f64>::zeros((n, n));
    for _iter in 0..cfg.n_iter {
        // 1. Advection-diffusion at fluid points -> intermediate velocity.
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
                let u_curr = u[[j, i]];
                let v_curr = v[[j, i]];

                let du_dx = if u_curr > 0.0 {
                    u_curr - u[[j, i - 1]]
                } else {
                    u[[j, i + 1]] - u_curr
                } / h;

                let du_dy = if v_curr > 0.0 {
                    u_curr - u[[j - 1, i]]
                } else {
                    u[[j + 1, i]] - u_curr
                } / h;

                let dv_dx = if u_curr > 0.0 {
                    v_curr - v[[j, i - 1]]
                } else {
                    v[[j, i + 1]] - v_curr
                } / h;

                let dv_dy = if v_curr > 0.0 {
                    v_curr - v[[j - 1, i]]
                } else {
                    v[[j + 1, i]] - v_curr
                } / h;

                // Diffusion (Central Difference); solid neighbours hold zero velocity.
                let lap_u = (u[[j, i + 1]] + u[[j, i - 1]] + u[[j + 1, i]] + u[[j - 1, i]]
                    - 4.0 * u_curr)
                    / (h * h);

                let lap_v = (v[[j, i + 1]] + v[[j, i - 1]] + v[[j + 1, i]] + v[[j - 1, i]]
                    - 4.0 * v_curr)
                    / (h * h);

                *u_val = dt.mul_add(-u_curr.mul_add(du_dx, v_curr * du_dy) + nu * lap_u, u_curr);
                *v_val = dt.mul_add(-u_curr.mul_add(dv_dx, v_curr * dv_dy) + nu * lap_v, v_curr);
            });

        // 2. Boundary conditions on the intermediate velocity.
        apply_channel_bcs(&mut u_star, &mut v_star, cfg);

        // 3. Pressure Poisson over fluid points: -lap(p) = -div(u*) / dt.
        pressure.solve(
            &mut p,
            u_star
                .as_slice()
                .ok_or("channel velocity storage is not contiguous")?,
            v_star
                .as_slice()
                .ok_or("channel velocity storage is not contiguous")?,
            obstacle_mask,
            dt,
        )?;

        // 4. Velocity correction with Neumann-mirrored pressure at solid faces.
        let pressure_at = |j: usize, i: usize, jc: usize, ic: usize| -> f64 {
            if obstacle_mask[[j, i]] {
                p[[jc, ic]]
            } else {
                p[[j, i]]
            }
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

                let next_u = dt.mul_add(-dp_dx, u_star[[j, i]]);
                let next_v = dt.mul_add(-dp_dy, v_star[[j, i]]);
                if !next_u.is_finite() || !next_v.is_finite() {
                    return Err("channel velocity became non-finite".into());
                }
                u[[j, i]] = next_u;
                v[[j, i]] = next_v;
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

/// Second-order 2D/3D MAC solver with obstacles and immersed bodies: see [`crate::sim::ns_mac`].
pub use crate::sim::ns_mac::ImmersedBody;
/// Second-order 2D/3D MAC solver with obstacles and immersed bodies: see [`crate::sim::ns_mac`].
pub use crate::sim::ns_mac::NsConfig;
/// Second-order 2D/3D MAC solver with obstacles and immersed bodies: see [`crate::sim::ns_mac`].
pub use crate::sim::ns_mac::NsSolver;

#[cfg(test)]
mod pressure_tests {
    use super::ChannelPressure;
    use ndarray::Array2;

    #[test]
    fn a_closed_two_cell_region_solves_compatible_nonzero_pressure() {
        let n = 5;
        let (j, i) = (2, 2);
        let h = 0.25;
        let dt = 0.01;
        let solid = Array2::from_shape_fn((n, n), |(y, x)| y != j || (x != i && x != i + 1));
        let mut u = vec![0.0; n * n];
        u[j * n + i] = 1.0;
        u[j * n + i + 1] = 1.0;
        let v = vec![0.0; n * n];
        let mut p = Array2::<f64>::zeros((n, n));
        ChannelPressure::new(&solid, h)
            .solve(&mut p, &u, &v, &solid, dt)
            .unwrap_or_else(|error| panic!("{error}"));
        assert_eq!(p[[j, i]].to_bits(), 0.0_f64.to_bits());
        assert!((p[[j, i + 1]] - h / (2.0 * dt)).abs() < 1e-13);
    }
}
