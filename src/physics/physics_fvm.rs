//! # Finite Volume Method (FVM) Module
//!
//! This module provides tools for solving partial differential equations, particularly
//! those that can be expressed in a conservative form, like advection equations.
//! It is well-suited for problems in fluid dynamics where the conservation of
//! quantities like mass, momentum, and energy is crucial.

use rayon::prelude::*;
use serde::{Deserialize, Serialize};

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
    pub fn new<F>(num_cells: usize, domain_size: f64, initial_conditions: F) -> Self
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
        Mesh { cells, dx }
    }

    /// Returns the number of cells in the mesh.
    #[inline]
    pub fn num_cells(&self) -> usize {
        self.cells.len()
    }
}

/// Calculates the numerical flux between two cells using the first-order upwind scheme.
#[inline]
pub(crate) fn upwind_flux(u_left: f64, u_right: f64, velocity: f64) -> f64 {
    if velocity > 0.0 {
        velocity * u_left
    } else {
        velocity * u_right
    }
}

/// Lax-Friedrichs numerical flux for a general conservation law u_t + f(u)_x = 0.
#[inline]
pub fn lax_friedrichs_flux<F>(u_left: f64, u_right: f64, dt: f64, dx: f64, flux_fn: F) -> f64 
where F: Fn(f64) -> f64 
{
    0.5 * (flux_fn(u_left) + flux_fn(u_right)) - 0.5 * (dx / dt) * (u_right - u_left)
}

/// Minmod limiter for MUSCL reconstruction.
#[inline]
pub fn minmod(a: f64, b: f64) -> f64 {
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
pub fn van_leer(a: f64, b: f64) -> f64 {
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
                *next_val = u_i - (dt / dx) * (flux_right - flux_left);
            });
        current_values.copy_from_slice(&next_values);
    }
    current_values
}

/// Example scenario: Simulates the advection of a square wave (top-hat profile).
pub fn simulate_1d_advection_scenario() -> Vec<f64> {
    const NUM_CELLS: usize = 200;
    const DOMAIN_SIZE: f64 = 1.0;
    const VELOCITY: f64 = 1.0;
    const CFL: f64 = 0.5;
    let dx = DOMAIN_SIZE / NUM_CELLS as f64;
    let dt = CFL * dx / VELOCITY.abs();
    let total_time = 0.5;
    let steps = (total_time / dt).ceil() as usize;
    let mut mesh = Mesh::new(NUM_CELLS, DOMAIN_SIZE, |x| {
        if x > 0.2 && x < 0.4 { 1.0 } else { 0.0 }
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
        next_values.par_iter_mut().enumerate().for_each(|(i, next_val)| {
            let u_i = current_values[i];
            let u_left = if i > 0 { current_values[i - 1] } else { current_values[0] };
            let u_right = if i < num_cells - 1 { current_values[i + 1] } else { current_values[num_cells - 1] };
            
            let f_left = lax_friedrichs_flux(u_left, u_i, dt, dx, flux_fn);
            let f_right = lax_friedrichs_flux(u_i, u_right, dt, dx, flux_fn);
            
            *next_val = u_i - (dt / dx) * (f_right - f_left);
        });
        current_values.copy_from_slice(&next_values);
    }
    current_values
}

/// State for 1D Shallow Water Equations: (h, hu)
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct SweState {
    pub h: f64,
    pub hu: f64,
}

/// Solves 1D Shallow Water Equations using FVM and Lax-Friedrichs flux.
pub fn solve_shallow_water_1d(
    initial_h: Vec<f64>,
    initial_hu: Vec<f64>,
    dx: f64,
    dt: f64,
    steps: usize,
    g: f64,
) -> Vec<SweState> {
    let n = initial_h.len();
    let mut current: Vec<SweState> = initial_h.into_iter().zip(initial_hu).map(|(h, hu)| SweState { h, hu }).collect();
    let mut next = current.clone();

    let flux_fn = |s: SweState| {
        let u = if s.h > 1e-6 { s.hu / s.h } else { 0.0 };
        (s.hu, s.hu * u + 0.5 * g * s.h * s.h)
    };

    for _ in 0..steps {
        next.par_iter_mut().enumerate().for_each(|(i, next_s)| {
            let s_i = current[i];
            let s_l = if i > 0 { current[i - 1] } else { current[0] };
            let s_r = if i < n - 1 { current[i + 1] } else { current[n - 1] };
            
            let (f_l_h, f_l_hu) = {
                let fl = flux_fn(s_l);
                let fi = flux_fn(s_i);
                (
                    0.5 * (fl.0 + fi.0) - 0.5 * (dx / dt) * (s_i.h - s_l.h),
                    0.5 * (fl.1 + fi.1) - 0.5 * (dx / dt) * (s_i.hu - s_l.hu)
                )
            };
            
            let (f_r_h, f_r_hu) = {
                let fi = flux_fn(s_i);
                let fr = flux_fn(s_r);
                (
                    0.5 * (fi.0 + fr.0) - 0.5 * (dx / dt) * (s_r.h - s_i.h),
                    0.5 * (fi.1 + fr.1) - 0.5 * (dx / dt) * (s_r.hu - s_i.hu)
                )
            };
            
            next_s.h = s_i.h - (dt / dx) * (f_r_h - f_l_h);
            next_s.hu = s_i.hu - (dt / dx) * (f_r_hu - f_l_hu);
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
        Mesh2D {
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
                *next_val = u_ij
                    - (dt / dx) * (flux_east - flux_west)
                    - (dt / dy) * (flux_north - flux_south);
            });
        current_values.copy_from_slice(&next_values);
    }
    current_values
}

/// Example scenario: Simulates the advection of a 2D Gaussian blob.
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
    let steps = (total_time / dt).ceil() as usize;
    let mut mesh = Mesh2D::new(WIDTH, HEIGHT, DOMAIN_SIZE, |x, y| {
        let (cx, cy) = (0.25, 0.5);
        let sigma_sq = 0.005;
        let dist_sq = (x - cx).powi(2) + (y - cy).powi(2);
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
        Mesh3D {
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
                *next_val = u_ijk
                    - (dt / dx) * (flux_east - flux_west)
                    - (dt / dy) * (flux_north - flux_south)
                    - (dt / dz) * (flux_front - flux_back);
            });
        current_values.copy_from_slice(&next_values);
    }
    current_values
}

/// Example scenario: Simulates the advection of a 3D Gaussian blob.
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
    let steps = (total_time / dt).ceil() as usize;
    let mut mesh = Mesh3D::new(WIDTH, HEIGHT, DEPTH, DOMAIN_SIZE, |x, y, z| {
        let (cx, cy, cz) = (0.3, 0.5, 0.5);
        let sigma_sq = 0.01;
        let dist_sq = (x - cx).powi(2) + (y - cy).powi(2) + (z - cz).powi(2);
        (-dist_sq / (2.0 * sigma_sq)).exp()
    });
    let boundary_conditions =
        |i, j, k, w, h, d| i == 0 || i == w - 1 || j == 0 || j == h - 1 || k == 0 || k == d - 1;
    solve_advection_3d(&mut mesh, velocity, dt, steps, boundary_conditions)
}
