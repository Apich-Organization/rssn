//! # Finite Volume Method (FVM) Solvers
//!
//! This module implements Finite Volume Method (FVM) solvers for hyperbolic partial differential equations (PDEs),
//! focusing on conservation laws. FVM is particularly well-suited for problems involving discontinuities (like shocks)
//! and for ensuring that physical quantities (mass, momentum, energy) are strictly conserved locally and globally.
//!
//! # Overview
//!
//! The module provides a framework for simulating fluid dynamics and wave propagation problems on structured grids
//! in 1D, 2D, and 3D. It discretizes the integral form of conservation laws, making it robust for calculating
//! flows with high gradients.
//!
//! Key components include:
//! - **Mesh Structures**: `Mesh` (1D), `Mesh2D`, and `Mesh3D` define the simulation domains composed of finite volumes (cells).
//! - **Conservation Solvers**:
//!     - **Linear Advection**: Solvers for 1D, 2D, and 3D advection equations using upwind schemes.
//!     - **Non-Linear Equations**: `solve_burgers_1d` for Burgers' equation and `solve_shallow_water_1d` for shallow water equations.
//! - **Numerical Fluxes**: Implementation of standard flux functions like `upwind_flux` and `lax_friedrichs_flux`.
//! - **Limiters**: Includes `minmod` and `van_leer` limiters (foundational for constructing high-resolution schemes like MUSCL).
//!
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/fvm_advection_frame_00.png)
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/fvm_advection_frame_03.png)
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/fvm_advection_frame_05.png)
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/fvm_advection_frame_07.png)
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/fvm_advection_frame_09.png)
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/fvm_advection_frame_12.png)
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/fvm_advection_frame_14.png)
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/fvm_advection_frame_16.png)
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/fvm_advection_frame_18.png)
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/fvm_advection_frame_19.png)
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/fvm_advection_frame_20.png)
//!
//! # Coverage audit
//!
//! | Equation | 1D | 2D | 3D | Scheme | Boundary conditions |
//! |---|---|---|---|---|---|
//! | linear advection | yes | yes | yes | first-order upwind | callback-selected fixed value |
//! | Burgers | yes | no | no | Lax–Friedrichs | fixed |
//! | shallow water | yes | yes | yes (barotropic, `γ = 2`) | Lax–Friedrichs | fixed |
//!
//! High-order finite volumes (MUSCL with four limiters, WENO5; HLLC, Roe and
//! Rusanov solvers; Euler and shallow water in 1D–3D; transmissive,
//! reflective, periodic and fixed-state faces) are in [`crate::sim::euler_nd`];
//! nodal discontinuous Galerkin (1D/2D) is in [`crate::sim::dg`].

use rayon::prelude::*;
use serde::Deserialize;
use serde::Serialize;

/// Represents a single cell or control volume in the mesh.
#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct Cell {
    /// The average value of the conserved quantity (e.g., density, concentration) in this cell.
    pub value: f64,
}

/// Represents a 1D simulation domain, composed of a series of cells.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Mesh {
    /// A vector of cells that make up the mesh.
    pub cells: Vec<Cell>,
    /// The size of each cell. Assumed to be uniform.
    pub dx: f64,
}

impl Mesh {
    /// Creates a new 1D mesh.
    ///
    /// # Arguments
    /// * `num_cells` - The number of cells in the mesh.
    /// * `domain_size` - The total length of the simulation domain.
    /// * `initial_conditions` - A function to set the initial value for each cell.
    pub fn new<F>(
        num_cells: usize,
        domain_size: f64,
        initial_conditions: F,
    ) -> Self
    where
        F: Fn(f64) -> f64,
    {
        let dx = domain_size / num_cells as f64;

        let cells = (0..num_cells)
            .map(|i| {
                let cell_center_x = (i as f64 + 0.5) * dx;

                Cell {
                    value: initial_conditions(cell_center_x),
                }
            })
            .collect();

        Self { cells, dx }
    }

    /// Returns the number of cells in the mesh.
    #[inline]
    #[must_use]
    pub const fn num_cells(&self) -> usize {
        self.cells.len()
    }
}

/// Calculates the numerical flux between two cells using the first-order upwind scheme.
#[inline]
pub(crate) fn upwind_flux(
    u_left: f64,
    u_right: f64,
    velocity: f64,
) -> f64 {
    if velocity > 0.0 {
        velocity * u_left
    } else {
        velocity * u_right
    }
}

/// Lax-Friedrichs numerical flux for a general conservation law `u_t` + f(u)_x = 0.
#[inline]
pub fn lax_friedrichs_flux<F>(
    u_left: f64,
    u_right: f64,
    dt: f64,
    dx: f64,
    flux_fn: F,
) -> f64
where
    F: Fn(f64) -> f64,
{
    0.5f64.mul_add(
        flux_fn(u_left) + flux_fn(u_right),
        -(0.5 * (dx / dt) * (u_right - u_left)),
    )
}

/// Minmod limiter for MUSCL reconstruction.
#[inline]
#[must_use]
pub fn minmod(
    a: f64,
    b: f64,
) -> f64 {
    if a * b <= 0.0 {
        0.0
    } else if a.abs() < b.abs() {
        a
    } else {
        b
    }
}

/// Van Leer limiter for MUSCL reconstruction.
#[inline]
#[must_use]
pub fn van_leer(
    a: f64,
    b: f64,
) -> f64 {
    if a * b <= 0.0 {
        0.0
    } else {
        2.0 * a * b / (a + b)
    }
}

/// Solves the 1D linear advection equation (`u_t + a * u_x = 0`) using the Finite Volume Method.
pub fn solve_advection_1d<F>(
    mesh: &mut Mesh,
    velocity: f64,
    dt: f64,
    steps: usize,
    boundary_conditions: F,
) -> Vec<f64>
where
    F: Fn() -> (f64, f64) + Sync,
{
    let num_cells = mesh.num_cells();

    let dx = mesh.dx;

    let mut current_values: Vec<f64> = mesh.cells.iter().map(|c| c.value).collect();

    let mut next_values = vec![0.0; num_cells];

    for _ in 0..steps {
        next_values
            .par_iter_mut()
            .enumerate()
            .for_each(|(i, next_val)| {
                let (bc_left, bc_right) = boundary_conditions();

                let u_i = current_values[i];

                let u_left = if i > 0 {
                    current_values[i - 1]
                } else {
                    bc_left
                };

                let u_right = if i < num_cells - 1 {
                    current_values[i + 1]
                } else {
                    bc_right
                };

                let flux_left = upwind_flux(u_left, u_i, velocity);

                let flux_right = upwind_flux(u_i, u_right, velocity);

                *next_val = (dt / dx).mul_add(-(flux_right - flux_left), u_i);
            });

        current_values.copy_from_slice(&next_values);
    }

    current_values
}

/// Example scenario: Simulates the advection of a square wave (top-hat profile).
#[must_use]
pub fn simulate_1d_advection_scenario() -> Vec<f64> {
    const NUM_CELLS: usize = 200;

    const DOMAIN_SIZE: f64 = 1.0;

    const VELOCITY: f64 = 1.0;

    const CFL: f64 = 0.5;

    let dx = DOMAIN_SIZE / NUM_CELLS as f64;

    let dt = CFL * dx / VELOCITY.abs();

    let total_time = 0.5;

    let steps: usize = ((total_time / dt).ceil() as i64).try_into().unwrap_or(0);

    let mut mesh = Mesh::new(NUM_CELLS, DOMAIN_SIZE, |x| {
        if x > 0.2 && x < 0.4 {
            1.0
        } else {
            0.0
        }
    });

    let boundary_conditions = || (0.0, 0.0);

    solve_advection_1d(&mut mesh, VELOCITY, dt, steps, boundary_conditions)
}

/// Solves the 1D Burgers' equation `u_t + (u^2/2)_x = 0` using FVM and Lax-Friedrichs flux.
pub fn solve_burgers_1d(
    mesh: &mut Mesh,
    dt: f64,
    steps: usize,
) -> Vec<f64> {
    let num_cells = mesh.num_cells();

    let dx = mesh.dx;

    let mut current_values: Vec<f64> = mesh.cells.iter().map(|c| c.value).collect();

    let mut next_values = vec![0.0; num_cells];

    let flux_fn = |u: f64| 0.5 * u * u;

    for _ in 0..steps {
        next_values
            .par_iter_mut()
            .enumerate()
            .for_each(|(i, next_val)| {
                let u_i = current_values[i];

                let u_left = if i > 0 {
                    current_values[i - 1]
                } else {
                    current_values[0]
                };

                let u_right = if i < num_cells - 1 {
                    current_values[i + 1]
                } else {
                    current_values[num_cells - 1]
                };

                let f_left = lax_friedrichs_flux(u_left, u_i, dt, dx, flux_fn);

                let f_right = lax_friedrichs_flux(u_i, u_right, dt, dx, flux_fn);

                *next_val = (dt / dx).mul_add(-(f_right - f_left), u_i);
            });

        current_values.copy_from_slice(&next_values);
    }

    current_values
}

/// State for 1D Shallow Water Equations: (h, hu)
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct SweState {
    /// The height of the water.
    pub h: f64,
    /// The momentum of the water.
    pub hu: f64,
}

/// Solves 1D Shallow Water Equations using FVM and Lax-Friedrichs flux.
#[must_use]
pub fn solve_shallow_water_1d(
    initial_h: Vec<f64>,
    initial_hu: Vec<f64>,
    dx: f64,
    dt: f64,
    steps: usize,
    g: f64,
) -> Vec<SweState> {
    let n = initial_h.len();

    let mut current: Vec<SweState> = initial_h
        .into_iter()
        .zip(initial_hu)
        .map(|(h, hu)| SweState { h, hu })
        .collect();

    let mut next = current.clone();

    let flux_fn = |s: SweState| {
        let u = if s.h > 1e-6 {
            s.hu / s.h
        } else {
            0.0
        };

        (s.hu, s.hu.mul_add(u, 0.5 * g * s.h * s.h))
    };

    for _ in 0..steps {
        next.par_iter_mut().enumerate().for_each(|(i, next_s)| {
            let s_i = current[i];

            let s_l = if i > 0 {
                current[i - 1]
            } else {
                current[0]
            };

            let s_r = if i < n - 1 {
                current[i + 1]
            } else {
                current[n - 1]
            };

            let (f_l_h, f_l_hu) = {
                let fl = flux_fn(s_l);

                let fi = flux_fn(s_i);

                (
                    0.5f64.mul_add(fl.0 + fi.0, -(0.5 * (dx / dt) * (s_i.h - s_l.h))),
                    0.5f64.mul_add(fl.1 + fi.1, -(0.5 * (dx / dt) * (s_i.hu - s_l.hu))),
                )
            };

            let (f_r_h, f_r_hu) = {
                let fi = flux_fn(s_i);

                let fr = flux_fn(s_r);

                (
                    0.5f64.mul_add(fi.0 + fr.0, -(0.5 * (dx / dt) * (s_r.h - s_i.h))),
                    0.5f64.mul_add(fi.1 + fr.1, -(0.5 * (dx / dt) * (s_r.hu - s_i.hu))),
                )
            };

            next_s.h = (dt / dx).mul_add(-(f_r_h - f_l_h), s_i.h);

            next_s.hu = (dt / dx).mul_add(-(f_r_hu - f_l_hu), s_i.hu);
        });

        current.copy_from_slice(&next);
    }

    current
}

/// Represents a 2D simulation domain as a grid of cells.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Mesh2D {
    /// A flattened vector of cells representing the 2D grid.
    pub cells: Vec<Cell>,
    /// The number of cells in the x-direction.
    pub width: usize,
    /// The number of cells in the y-direction.
    pub height: usize,
    /// The size of each cell in the x-direction.
    pub dx: f64,
    /// The size of each cell in the y-direction.
    pub dy: f64,
}

impl Mesh2D {
    /// Creates a new 2D mesh.
    pub fn new<F>(
        width: usize,
        height: usize,
        domain_size: (f64, f64),
        initial_conditions: F,
    ) -> Self
    where
        F: Fn(f64, f64) -> f64,
    {
        let dx = domain_size.0 / width as f64;

        let dy = domain_size.1 / height as f64;

        let mut cells = Vec::with_capacity(width * height);

        for j in 0..height {
            for i in 0..width {
                let center_x = (i as f64 + 0.5) * dx;

                let center_y = (j as f64 + 0.5) * dy;

                cells.push(Cell {
                    value: initial_conditions(center_x, center_y),
                });
            }
        }

        Self {
            cells,
            width,
            height,
            dx,
            dy,
        }
    }
}

/// Solves the 2D linear advection equation using the Finite Volume Method.
pub fn solve_advection_2d<F>(
    mesh: &mut Mesh2D,
    velocity: (f64, f64),
    dt: f64,
    steps: usize,
    boundary_conditions: F,
) -> Vec<f64>
where
    F: Fn(usize, usize, usize, usize) -> bool + Sync,
{
    let (width, height) = (mesh.width, mesh.height);

    let (dx, dy) = (mesh.dx, mesh.dy);

    let mut current_values: Vec<f64> = mesh.cells.iter().map(|c| c.value).collect();

    let mut next_values = vec![0.0; width * height];

    for _ in 0..steps {
        next_values
            .par_iter_mut()
            .enumerate()
            .for_each(|(idx, next_val)| {
                let i = idx % width;

                let j = idx / width;

                if boundary_conditions(i, j, width, height) {
                    *next_val = 0.0;

                    return;
                }

                let u_ij = current_values[idx];

                let u_left = current_values[idx - 1];

                let u_right = current_values[idx + 1];

                let u_down = current_values[idx - width];

                let u_up = current_values[idx + width];

                let flux_west = upwind_flux(u_left, u_ij, velocity.0);

                let flux_east = upwind_flux(u_ij, u_right, velocity.0);

                let flux_south = upwind_flux(u_down, u_ij, velocity.1);

                let flux_north = upwind_flux(u_ij, u_up, velocity.1);

                *next_val = (dt / dy).mul_add(
                    -(flux_north - flux_south),
                    (dt / dx).mul_add(-(flux_east - flux_west), u_ij),
                );
            });

        current_values.copy_from_slice(&next_values);
    }

    current_values
}

/// Example scenario: Simulates the advection of a 2D Gaussian blob.
#[must_use]
pub fn simulate_2d_advection_scenario() -> Vec<f64> {
    const WIDTH: usize = 100;

    const HEIGHT: usize = 100;

    const DOMAIN_SIZE: (f64, f64) = (1.0, 1.0);

    const CFL: f64 = 0.4;

    let velocity = (0.5, 0.3);

    let dx = DOMAIN_SIZE.0 / WIDTH as f64;

    let dy = DOMAIN_SIZE.1 / HEIGHT as f64;

    let dt = CFL * (dx.min(dy)) / (f64::abs(velocity.0) + f64::abs(velocity.1)).max(1e-6_f64);

    let total_time = 0.6;

    let steps: usize = ((total_time / dt).ceil() as i64).try_into().unwrap_or(0);

    let mut mesh = Mesh2D::new(WIDTH, HEIGHT, DOMAIN_SIZE, |x, y| {
        let (cx, cy) = (0.25, 0.5);

        let sigma_sq = 0.005;

        let dist_sq = (y - cy).mul_add(y - cy, (x - cx).powi(2));

        (-dist_sq / (2.0 * sigma_sq)).exp()
    });

    let boundary_conditions = |i: usize, j: usize, width: usize, height: usize| -> bool {
        i == 0 || i == width - 1 || j == 0 || j == height - 1
    };

    solve_advection_2d(&mut mesh, velocity, dt, steps, boundary_conditions)
}

/// Represents a 3D simulation domain as a grid of cells.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Mesh3D {
    /// A flattened vector of cells representing the 3D grid.
    pub cells: Vec<Cell>,
    /// The number of cells in the x-direction.
    pub width: usize,
    /// The number of cells in the y-direction.
    pub height: usize,
    /// The number of cells in the z-direction.
    pub depth: usize,
    /// The size of each cell in the x-direction.
    pub dx: f64,
    /// The size of each cell in the y-direction.
    pub dy: f64,
    /// The size of each cell in the z-direction.
    pub dz: f64,
}

impl Mesh3D {
    /// Creates a new 3D mesh.
    pub fn new<F>(
        width: usize,
        height: usize,
        depth: usize,
        domain_size: (f64, f64, f64),
        initial_conditions: F,
    ) -> Self
    where
        F: Fn(f64, f64, f64) -> f64,
    {
        let dx = domain_size.0 / width as f64;

        let dy = domain_size.1 / height as f64;

        let dz = domain_size.2 / depth as f64;

        let mut cells = Vec::with_capacity(width * height * depth);

        for k in 0..depth {
            for j in 0..height {
                for i in 0..width {
                    let center_x = (i as f64 + 0.5) * dx;

                    let center_y = (j as f64 + 0.5) * dy;

                    let center_z = (k as f64 + 0.5) * dz;

                    cells.push(Cell {
                        value: initial_conditions(center_x, center_y, center_z),
                    });
                }
            }
        }

        Self {
            cells,
            width,
            height,
            depth,
            dx,
            dy,
            dz,
        }
    }
}

/// Solves the 3D linear advection equation using FVM.
pub fn solve_advection_3d<F>(
    mesh: &mut Mesh3D,
    velocity: (f64, f64, f64),
    dt: f64,
    steps: usize,
    boundary_conditions: F,
) -> Vec<f64>
where
    F: Fn(usize, usize, usize, usize, usize, usize) -> bool + Sync,
{
    let (width, height, depth) = (mesh.width, mesh.height, mesh.depth);

    let (dx, dy, dz) = (mesh.dx, mesh.dy, mesh.dz);

    let mut current_values: Vec<f64> = mesh.cells.iter().map(|c| c.value).collect();

    let mut next_values = vec![0.0; width * height * depth];

    let plane_size = width * height;

    for _ in 0..steps {
        next_values
            .par_iter_mut()
            .enumerate()
            .for_each(|(idx, next_val)| {
                let i = idx % width;

                let j = (idx / width) % height;

                let k = idx / plane_size;

                if boundary_conditions(i, j, k, width, height, depth) {
                    *next_val = 0.0;

                    return;
                }

                let u_ijk = current_values[idx];

                let u_left = current_values[idx - 1];

                let u_right = current_values[idx + 1];

                let u_down = current_values[idx - width];

                let u_up = current_values[idx + width];

                let u_back = current_values[idx - plane_size];

                let u_front = current_values[idx + plane_size];

                let flux_west = upwind_flux(u_left, u_ijk, velocity.0);

                let flux_east = upwind_flux(u_ijk, u_right, velocity.0);

                let flux_south = upwind_flux(u_down, u_ijk, velocity.1);

                let flux_north = upwind_flux(u_ijk, u_up, velocity.1);

                let flux_back = upwind_flux(u_back, u_ijk, velocity.2);

                let flux_front = upwind_flux(u_ijk, u_front, velocity.2);

                *next_val = (dt / dz).mul_add(
                    -(flux_front - flux_back),
                    (dt / dy).mul_add(
                        -(flux_north - flux_south),
                        (dt / dx).mul_add(-(flux_east - flux_west), u_ijk),
                    ),
                );
            });

        current_values.copy_from_slice(&next_values);
    }

    current_values
}

/// Example scenario: Simulates the advection of a 3D Gaussian blob.
#[must_use]
pub fn simulate_3d_advection_scenario() -> Vec<f64> {
    const WIDTH: usize = 30;

    const HEIGHT: usize = 30;

    const DEPTH: usize = 30;

    const DOMAIN_SIZE: (f64, f64, f64) = (1.0, 1.0, 1.0);

    const CFL: f64 = 0.3;

    let velocity = (0.5, 0.3, 0.1);

    let dx = DOMAIN_SIZE.0 / WIDTH as f64;

    let dy = DOMAIN_SIZE.1 / HEIGHT as f64;

    let dz = DOMAIN_SIZE.2 / DEPTH as f64;

    let min_dim = dx.min(dy).min(dz);

    let vel_mag =
        (f64::abs(velocity.0) + f64::abs(velocity.1) + f64::abs(velocity.2)).max(1e-6_f64);

    let dt = CFL * min_dim / vel_mag;

    let total_time = 0.7;

    let steps: usize = ((total_time / dt).ceil() as i64).try_into().unwrap_or(0);

    let mut mesh = Mesh3D::new(WIDTH, HEIGHT, DEPTH, DOMAIN_SIZE, |x, y, z| {
        let (cx, cy, cz) = (0.3, 0.5, 0.5);

        let sigma_sq = 0.01;

        let dist_sq = (z - cz).mul_add(z - cz, (y - cy).mul_add(y - cy, (x - cx).powi(2)));

        (-dist_sq / (2.0 * sigma_sq)).exp()
    });

    let boundary_conditions =
        |i, j, k, w, h, d| i == 0 || i == w - 1 || j == 0 || j == h - 1 || k == 0 || k == d - 1;

    solve_advection_3d(&mut mesh, velocity, dt, steps, boundary_conditions)
}

/// High-order finite volumes (MUSCL / WENO5, HLLC / Roe, Euler and shallow
/// water in 1D–3D) of [`crate::sim::euler_nd`], the successor of the
/// first-order solvers of this module.
pub use crate::sim::euler_nd::FvConfig;
/// High-order finite volumes (MUSCL / WENO5, HLLC / Roe, Euler and shallow
/// water in 1D–3D) of [`crate::sim::euler_nd`], the successor of the
/// first-order solvers of this module.
pub use crate::sim::euler_nd::FvGrid;
/// High-order finite volumes (MUSCL / WENO5, HLLC / Roe, Euler and shallow
/// water in 1D–3D) of [`crate::sim::euler_nd`], the successor of the
/// first-order solvers of this module.
pub use crate::sim::euler_nd::GasBc;
/// High-order finite volumes (MUSCL / WENO5, HLLC / Roe, Euler and shallow
/// water in 1D–3D) of [`crate::sim::euler_nd`], the successor of the
/// first-order solvers of this module.
pub use crate::sim::euler_nd::Law;
/// High-order finite volumes (MUSCL / WENO5, HLLC / Roe, Euler and shallow
/// water in 1D–3D) of [`crate::sim::euler_nd`], the successor of the
/// first-order solvers of this module.
pub use crate::sim::euler_nd::Reconstruction;
/// High-order finite volumes (MUSCL / WENO5, HLLC / Roe, Euler and shallow
/// water in 1D–3D) of [`crate::sim::euler_nd`], the successor of the
/// first-order solvers of this module.
pub use crate::sim::euler_nd::RiemannSolver;
/// High-order finite volumes (MUSCL / WENO5, HLLC / Roe, Euler and shallow
/// water in 1D–3D) of [`crate::sim::euler_nd`], the successor of the
/// first-order solvers of this module.
pub use crate::sim::euler_nd::SlopeLimiter;
/// High-order finite volumes (MUSCL / WENO5, HLLC / Roe, Euler and shallow
/// water in 1D–3D) of [`crate::sim::euler_nd`], the successor of the
/// first-order solvers of this module.
pub use crate::sim::euler_nd::solve as solve_conservation_law_nd;

// ============================================================================
// Space-dimension coverage: Burgers 1D/2D/3D and shallow water 2D/3D
// ============================================================================
//
// First-order Lax-Friedrichs finite volumes on cell-centred grids with
// zero-gradient (outflow) boundaries. High-order MUSCL/WENO solvers for Euler
// and shallow water in 1D-3D are in [`crate::sim::euler_nd`] (re-exported at the end of this module).

/// Burgers' equation `u_t + Σ_a (u²/2)_{x_a} = 0` on a 1D, 2D or 3D grid.
///
/// Cells (`dims` cells per axis, `h` cell sizes, flat, `x` fastest) with the
/// Lax-Friedrichs flux (numerical viscosity shared between the axes, as the
/// unsplit multi-dimensional scheme requires) and zero-gradient boundaries. With
/// `dims.len() == 1` it reproduces [`solve_burgers_1d`]; the scheme is stable for
/// `max|u| dt Σ_a 1/h_a <= 1`.
///
/// # Errors
/// For 0 or more than 3 axes, a length mismatch or non-positive `dt`/`h`.
pub fn solve_burgers_nd(
    u0: &[f64],
    dims: &[usize],
    h: &[f64],
    dt: f64,
    steps: usize,
) -> Result<Vec<f64>, String> {
    if dims.is_empty() || dims.len() > 3 || dims.len() != h.len() || dims.contains(&0) {
        return Err("1 to 3 axes with one spacing and a positive size each".to_owned());
    }

    if u0.len() != dims.iter().product::<usize>()
        || dt.is_nan()
        || dt <= 0.0
        || h.iter().any(|&v| v.is_nan() || v <= 0.0)
    {
        return Err("array length / step sizes inconsistent".to_owned());
    }

    let dim = dims.len();

    let mut n = [1usize; 3];

    n[..dim].copy_from_slice(dims);

    let stride = [1, n[0], n[0] * n[1]];

    let flux = |u: f64| 0.5 * u * u;

    let mut cur = u0.to_vec();

    let mut next = vec![0.0; cur.len()];

    for _ in 0..steps {
        next.par_iter_mut().enumerate().for_each(|(p, out)| {
            let u = cur[p];

            let q = [p % n[0], (p / stride[1]) % n[1], p / stride[2]];

            let mut value = u;

            for a in 0..dim {
                let left = if q[a] > 0 {
                    cur[p - stride[a]]
                } else {
                    u
                };

                let right = if q[a] + 1 < n[a] {
                    cur[p + stride[a]]
                } else {
                    u
                };

                // the numerical viscosity is shared between the axes (dim * dt), the
                // unsplit multi-dimensional Lax-Friedrichs condition
                let dt_lf = dt * dim as f64;

                let f_left = lax_friedrichs_flux(left, u, dt_lf, h[a], flux);

                let f_right = lax_friedrichs_flux(u, right, dt_lf, h[a], flux);

                value -= dt / h[a] * (f_right - f_left);
            }

            *out = value;
        });

        std::mem::swap(&mut cur, &mut next);
    }

    Ok(cur)
}

/// 2D Burgers' equation; 2D counterpart of [`solve_burgers_1d`]. See
/// [`solve_burgers_nd`].
///
/// # Errors
/// As [`solve_burgers_nd`].
pub fn solve_burgers_2d(
    u0: &[f64],
    dims: [usize; 2],
    h: [f64; 2],
    dt: f64,
    steps: usize,
) -> Result<Vec<f64>, String> {
    solve_burgers_nd(u0, &dims, &h, dt, steps)
}

/// 3D Burgers' equation; see [`solve_burgers_nd`].
///
/// # Errors
/// As [`solve_burgers_nd`].
pub fn solve_burgers_3d(
    u0: &[f64],
    dims: [usize; 3],
    h: [f64; 3],
    dt: f64,
    steps: usize,
) -> Result<Vec<f64>, String> {
    solve_burgers_nd(u0, &dims, &h, dt, steps)
}

/// State of the 2D shallow water equations: depth and the two momenta.
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct SweState2D {
    /// The height of the water.
    pub h: f64,
    /// Momentum in x.
    pub hu: f64,
    /// Momentum in y.
    pub hv: f64,
}

/// 2D shallow water equations `h_t + (hu)_x + (hv)_y = 0`, ...
///
/// On a grid of `dims = [nx, ny]` cells (flat, `x` fastest) with the Lax-Friedrichs flux and
/// zero-gradient boundaries. 2D counterpart of [`solve_shallow_water_1d`]; the
/// numerical viscosity is shared between the two axes (unsplit scheme), so it
/// smears slightly more than the 1D routine on y-independent data. Stable for
/// `(|u| + sqrt(g h)) dt (1/dx + 1/dy) <= 1` (roughly).
///
/// # Errors
/// For a length mismatch or non-positive step sizes.
#[allow(clippy::too_many_arguments)]
pub fn solve_shallow_water_2d(
    initial_h: &[f64],
    initial_hu: &[f64],
    initial_hv: &[f64],
    dims: [usize; 2],
    spacing: [f64; 2],
    dt: f64,
    steps: usize,
    g: f64,
) -> Result<Vec<SweState2D>, String> {
    let [nx, ny] = dims;

    if nx == 0
        || ny == 0
        || initial_h.len() != nx * ny
        || initial_hu.len() != nx * ny
        || initial_hv.len() != nx * ny
        || dt.is_nan()
        || dt <= 0.0
        || spacing.iter().any(|&v| v.is_nan() || v <= 0.0)
    {
        return Err("array length / step sizes inconsistent".to_owned());
    }

    let mut current: Vec<SweState2D> = (0..nx * ny)
        .map(|p| {
            SweState2D {
                h: initial_h[p],
                hu: initial_hu[p],
                hv: initial_hv[p],
            }
        })
        .collect();

    let mut next = current.clone();

    let flux = |s: SweState2D, axis: usize| -> [f64; 3] {
        let (u, v) = if s.h > 1e-6 {
            (s.hu / s.h, s.hv / s.h)
        } else {
            (0.0, 0.0)
        };

        let pressure = 0.5 * g * s.h * s.h;

        if axis == 0 {
            [s.hu, s.hu.mul_add(u, pressure), s.hu * v]
        } else {
            [s.hv, s.hv * u, s.hv.mul_add(v, pressure)]
        }
    };

    let components = |s: SweState2D| [s.h, s.hu, s.hv];

    for _ in 0..steps {
        next.par_iter_mut().enumerate().for_each(|(p, out)| {
            let (i, j) = (p % nx, p / nx);

            let c = current[p];

            let mut value = components(c);

            for axis in 0..2 {
                let (stride, index, count) = if axis == 0 {
                    (1, i, nx)
                } else {
                    (nx, j, ny)
                };

                let left = if index > 0 {
                    current[p - stride]
                } else {
                    c
                };

                let right = if index + 1 < count {
                    current[p + stride]
                } else {
                    c
                };

                let h_a = spacing[axis];

                let numerical = |a: SweState2D, b: SweState2D| -> [f64; 3] {
                    let (fa, fb) = (flux(a, axis), flux(b, axis));

                    let (ca, cb) = (components(a), components(b));

                    [0, 1, 2]
                        .map(|k| f64::midpoint(fa[k], fb[k]) - 0.25 * (h_a / dt) * (cb[k] - ca[k]))
                };

                let (f_left, f_right) = (numerical(left, c), numerical(c, right));

                for k in 0..3 {
                    value[k] -= dt / h_a * (f_right[k] - f_left[k]);
                }
            }

            *out = SweState2D {
                h: value[0],
                hu: value[1],
                hv: value[2],
            };
        });

        std::mem::swap(&mut current, &mut next);
    }

    Ok(current)
}

/// State of the 3D barotropic shallow-water system: depth and the three momenta.
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct SweState3D {
    /// The height of the water.
    pub h: f64,
    /// Momentum in x.
    pub hu: f64,
    /// Momentum in y.
    pub hv: f64,
    /// Momentum in z.
    pub hw: f64,
}

/// 3D barotropic shallow-water system `h_t + ∇·(h u) = 0`,
/// `(h u)_t + ∇·(h u ⊗ u) + ∇(g h²/2) = 0`.
///
/// Euler with `γ = 2`, `p = g h²/2` and no
/// energy equation: the 3D analogue of the planar shallow-water equations.
///
/// On a grid of `dims = [nx, ny, nz]` cells (flat, `x` fastest) with the
/// Lax-Friedrichs flux and zero-gradient boundaries. 3D counterpart of
/// [`solve_shallow_water_2d`]; the numerical viscosity is shared between the
/// three axes (unsplit scheme). Stable for
/// `(|u| + sqrt(g h)) dt (1/dx + 1/dy + 1/dz) <= 1` (roughly).
///
/// # Errors
/// For a length mismatch or non-positive step sizes.
#[allow(clippy::too_many_arguments)]
pub fn solve_shallow_water_3d(
    initial_h: &[f64],
    initial_hu: &[f64],
    initial_hv: &[f64],
    initial_hw: &[f64],
    dims: [usize; 3],
    spacing: [f64; 3],
    dt: f64,
    steps: usize,
    g: f64,
) -> Result<Vec<SweState3D>, String> {
    let [nx, ny, nz] = dims;
    let total = nx * ny * nz;

    if nx == 0
        || ny == 0
        || nz == 0
        || initial_h.len() != total
        || initial_hu.len() != total
        || initial_hv.len() != total
        || initial_hw.len() != total
        || dt.is_nan()
        || dt <= 0.0
        || spacing.iter().any(|&v| v.is_nan() || v <= 0.0)
    {
        return Err("array length / step sizes inconsistent".to_owned());
    }

    let mut current: Vec<SweState3D> = (0..total)
        .map(|p| {
            SweState3D {
                h: initial_h[p],
                hu: initial_hu[p],
                hv: initial_hv[p],
                hw: initial_hw[p],
            }
        })
        .collect();

    let mut next = current.clone();

    let components = |s: SweState3D| [s.h, s.hu, s.hv, s.hw];

    let flux = |s: SweState3D, axis: usize| -> [f64; 4] {
        let m = components(s);

        let pressure = 0.5 * g * s.h * s.h;

        let mut f = [0.0; 4];

        f[0] = m[1 + axis];

        for k in 1..4 {
            let velocity = if s.h > 1e-6 {
                m[1 + axis] / s.h
            } else {
                0.0
            };

            f[k] = m[k] * velocity;
        }

        f[1 + axis] += pressure;

        f
    };

    for _ in 0..steps {
        next.par_iter_mut().enumerate().for_each(|(p, out)| {
            let q = [p % nx, (p / nx) % ny, p / (nx * ny)];

            let stride = [1, nx, nx * ny];

            let c = current[p];

            let mut value = components(c);

            for axis in 0..3 {
                let count = dims[axis];

                let left = if q[axis] > 0 {
                    current[p - stride[axis]]
                } else {
                    c
                };

                let right = if q[axis] + 1 < count {
                    current[p + stride[axis]]
                } else {
                    c
                };

                let h_a = spacing[axis];

                // 0.5 * h / (3 dt): the viscosity shared between the three axes
                let numerical = |a: SweState3D, b: SweState3D| -> [f64; 4] {
                    let (fa, fb) = (flux(a, axis), flux(b, axis));

                    let (ca, cb) = (components(a), components(b));

                    [0, 1, 2, 3]
                        .map(|k| f64::midpoint(fa[k], fb[k]) - (h_a / (6.0 * dt)) * (cb[k] - ca[k]))
                };

                let (f_left, f_right) = (numerical(left, c), numerical(c, right));

                for k in 0..4 {
                    value[k] -= dt / h_a * (f_right[k] - f_left[k]);
                }
            }

            *out = SweState3D {
                h: value[0],
                hu: value[1],
                hv: value[2],
                hw: value[3],
            };
        });

        std::mem::swap(&mut current, &mut next);
    }

    Ok(current)
}
