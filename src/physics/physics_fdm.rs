//! # Finite Difference Method (FDM) Module
//!
//! This module provides tools for solving partial differential equations (PDEs)
//! using the finite difference method. It includes a generic grid structure
//! and a solver for the 2D heat equation as an example.

use std::ops::Index;
use std::ops::IndexMut;

use rayon::prelude::*;
use serde::Deserialize;
use serde::Serialize;

/// Represents the dimensions of the simulation grid.
/// cbindgen:ignore
#[derive(
    Clone,
    Debug,
    Serialize,
    Deserialize,
    PartialEq,
    Eq,
)]

pub enum Dimensions {
    /// 1-dimensional grid with a given size.
    D1(usize),
    /// 2-dimensional grid with width and height.
    D2(usize, usize),
    /// 3-dimensional grid with width, height, and depth.
    D3(usize, usize, usize),
}

/// A generic grid structure for finite difference method simulations.
/// It can represent a 1D, 2D, or 3D grid.
/// cbindgen:ignore
#[derive(
    Clone, Debug, Serialize, Deserialize,
)]

pub struct FdmGrid<T> {
    pub(crate) data: Vec<T>,
    pub(crate) dims: Dimensions,
}

impl<
        T: Clone + Default + Send + Sync,
    > FdmGrid<T>
{
    /// Creates a new grid from existing data and dimensions.

    #[must_use] 
    pub const fn from_data(
        data: Vec<T>,
        dims: Dimensions,
    ) -> Self {

        Self {
            data,
            dims,
        }
    }

    /// Creates a new grid with the given dimensions, initialized with a default value.

    #[must_use] 
    pub fn new(
        dims: Dimensions
    ) -> Self {

        let size = match dims {
            | Dimensions::D1(x) => x,
            | Dimensions::D2(x, y) => {
                x * y
            },
            | Dimensions::D3(
                x,
                y,
                z,
            ) => x * y * z,
        };

        Self {
            data: vec![
                T::default();
                size
            ],
            dims,
        }
    }

    /// Creates a new grid with the given dimensions, initialized with a specified value.

    pub fn with_value(
        dims: Dimensions,
        value: T,
    ) -> Self {

        let size = match dims {
            | Dimensions::D1(x) => x,
            | Dimensions::D2(x, y) => {
                x * y
            },
            | Dimensions::D3(
                x,
                y,
                z,
            ) => x * y * z,
        };

        Self {
            data: vec![value; size],
            dims,
        }
    }

    /// Returns the dimensions of the grid.
    #[inline]

    #[must_use] 
    pub const fn dimensions(
        &self
    ) -> &Dimensions {

        &self.dims
    }

    /// Returns a slice to the underlying data.
    #[inline]

    #[must_use] 
    pub fn as_slice(&self) -> &[T] {

        &self.data
    }

    /// Returns a mutable slice to the underlying data.
    #[inline]

    pub fn as_mut_slice(
        &mut self
    ) -> &mut [T] {

        &mut self.data
    }

    /// Size of the grid.

    #[must_use] 
    pub const fn len(&self) -> usize {

        self.data.len()
    }

    /// Is the grid empty.

    #[must_use] 
    pub const fn is_empty(&self) -> bool {

        self.data.is_empty()
    }
}

impl<T> Index<usize> for FdmGrid<T> {
    type Output = T;

    #[inline]

    fn index(
        &self,
        index: usize,
    ) -> &Self::Output {

        &self.data[index]
    }
}

impl<T> IndexMut<usize> for FdmGrid<T> {
    #[inline]

    fn index_mut(
        &mut self,
        index: usize,
    ) -> &mut Self::Output {

        &mut self.data[index]
    }
}

impl<T> Index<(usize, usize)>
    for FdmGrid<T>
{
    type Output = T;

    #[inline]

    fn index(
        &self,
        (x, y): (usize, usize),
    ) -> &Self::Output {

        if let Dimensions::D2(
            width,
            _,
        ) = self.dims
        {

            &self.data[y * width + x]
        } else {

            panic!(
                "Attempted to use 2D \
                 indexing on a non-2D \
                 grid."
            );
        }
    }
}

impl<T> IndexMut<(usize, usize)>
    for FdmGrid<T>
{
    #[inline]

    fn index_mut(
        &mut self,
        (x, y): (usize, usize),
    ) -> &mut Self::Output {

        if let Dimensions::D2(
            width,
            _,
        ) = self.dims
        {

            &mut self.data
                [y * width + x]
        } else {

            panic!(
                "Attempted to use 2D \
                 indexing on a non-2D \
                 grid."
            );
        }
    }
}

impl<T> Index<(usize, usize, usize)>
    for FdmGrid<T>
{
    type Output = T;

    #[inline]

    fn index(
        &self,
        (x, y, z): (
            usize,
            usize,
            usize,
        ),
    ) -> &Self::Output {

        if let Dimensions::D3(
            width,
            height,
            _,
        ) = self.dims
        {

            &self.data[z
                * width
                * height
                + y * width
                + x]
        } else {

            panic!(
                "Attempted to use 3D \
                 indexing on a non-3D \
                 grid."
            );
        }
    }
}

impl<T> IndexMut<(usize, usize, usize)>
    for FdmGrid<T>
{
    #[inline]

    fn index_mut(
        &mut self,
        (x, y, z): (
            usize,
            usize,
            usize,
        ),
    ) -> &mut Self::Output {

        if let Dimensions::D3(
            width,
            height,
            _,
        ) = self.dims
        {

            &mut self.data[z
                * width
                * height
                + y * width
                + x]
        } else {

            panic!(
                "Attempted to use 3D \
                 indexing on a non-3D \
                 grid."
            );
        }
    }
}

// ============================================================================
// Solvers
// ============================================================================

/// Solves a 2D heat equation `u_t = alpha * ∇²u` using the finite difference method.

pub fn solve_heat_equation_2d<F>(
    width: usize,
    height: usize,
    alpha: f64,
    dx: f64,
    dy: f64,
    dt: f64,
    steps: usize,
    initial_conditions: F,
) -> FdmGrid<f64>
where
    F: Fn(usize, usize) -> f64 + Sync,
{

    let dims =
        Dimensions::D2(width, height);

    let mut grid = FdmGrid::new(dims);

    let mut next_grid = grid.clone();

    grid.as_mut_slice()
        .par_iter_mut()
        .enumerate()
        .for_each(|(i, val)| {

            let x = i % width;

            let y = i / width;

            *val = initial_conditions(
                x, y,
            );
        });

    let r_x = alpha * dt / dx.powi(2);

    let r_y = alpha * dt / dy.powi(2);

    for _ in 0 .. steps {

        next_grid
            .as_mut_slice()
            .par_iter_mut()
            .enumerate()
            .for_each(
                |(i, next_val)| {

                    let x = i % width;

                    let y = i / width;

                    if x == 0
                        || x == width
                            - 1
                        || y == 0
                        || y == height
                            - 1
                    {

                        *next_val =
                            grid[i]; // Boundary remains constant
                        return;
                    }

                    let lap_x = 2.0f64.mul_add(-grid[(
                                x, y,
                            )], grid
                        [(x + 1, y)])
                        + grid[(
                            x - 1,
                            y,
                        )];

                    let lap_y = 2.0f64.mul_add(-grid[(
                                x, y,
                            )], grid
                        [(x, y + 1)])
                        + grid[(
                            x,
                            y - 1,
                        )];

                    *next_val = grid[i]
                        + r_x * lap_x
                        + r_y * lap_y;
                },
            );

        std::mem::swap(
            &mut grid,
            &mut next_grid,
        );
    }

    grid
}

/// Solves 2D Wave equation `u_tt = c^2 * ∇²u`.

pub fn solve_wave_equation_2d<F>(
    width: usize,
    height: usize,
    c: f64,
    dx: f64,
    dy: f64,
    dt: f64,
    steps: usize,
    initial_u: F,
) -> FdmGrid<f64>
where
    F: Fn(usize, usize) -> f64 + Sync,
{

    let dims =
        Dimensions::D2(width, height);

    let mut u_prev =
        FdmGrid::new(dims.clone());

    let mut u_curr =
        FdmGrid::new(dims.clone());

    let mut u_next =
        FdmGrid::new(dims);

    u_curr
        .as_mut_slice()
        .par_iter_mut()
        .enumerate()
        .for_each(|(i, val)| {

            let x = i % width;

            let y = i / width;

            *val = initial_u(x, y);
        });

    // First step (assuming u_t = 0 at t=0)
    let s_x = (c * dt / dx).powi(2);

    let s_y = (c * dt / dy).powi(2);

    u_prev
        .data
        .copy_from_slice(&u_curr.data);

    for _ in 0 .. steps {

        u_next
            .as_mut_slice()
            .par_iter_mut()
            .enumerate()
            .for_each(
                |(i, next_val)| {

                    let x = i % width;

                    let y = i / width;

                    if x == 0
                        || x == width
                            - 1
                        || y == 0
                        || y == height
                            - 1
                    {

                        *next_val = 0.0;

                        return;
                    }

                    let lap_x = 2.0f64.mul_add(-u_curr[(
                                x, y,
                            )], u_curr
                        [(x + 1, y)])
                        + u_curr[(
                            x - 1,
                            y,
                        )];

                    let lap_y = 2.0f64.mul_add(-u_curr[(
                                x, y,
                            )], u_curr
                        [(x, y + 1)])
                        + u_curr[(
                            x,
                            y - 1,
                        )];

                    *next_val = 2.0f64.mul_add(u_curr[i], -u_prev[i])
                        + s_x * lap_x
                        + s_y * lap_y;
                },
            );

        std::mem::swap(
            &mut u_prev,
            &mut u_curr,
        );

        std::mem::swap(
            &mut u_curr,
            &mut u_next,
        );
    }

    u_curr
}

/// Solves Poisson equation `∇²u = f` using Successive Over-Relaxation (SOR).
///
/// # Panics
///
/// Panics if the `source` grid's dimensions do not match the `width` and `height` provided.
/// Panics if internal indexing operations go out of bounds due to incorrect `FdmGrid` usage,
/// though typical usage within this function should prevent this.

#[must_use] 
pub fn solve_poisson_2d(
    width: usize,
    height: usize,
    source: &FdmGrid<f64>,
    dx: f64,
    dy: f64,
    omega: f64,
    max_iter: usize,
    tolerance: f64,
) -> FdmGrid<f64> {

    let dims =
        Dimensions::D2(width, height);

    let mut u: FdmGrid<f64> =
        FdmGrid::new(dims);

    let dx2 = dx * dx;

    let dy2 = dy * dy;

    let beta = dx2 / dy2;

    let factor =
        1.0 / (2.0 * (1.0 + beta));

    for _ in 0 .. max_iter {

        let mut max_diff: f64 = 0.0;

        for pass in 0 .. 2 {

            let u_ptr = u
                .as_mut_slice()
                .as_mut_ptr()
                as usize;

            let diffs : Vec<f64> = u
                .as_mut_slice()
                .par_iter_mut()
                .enumerate()
                .filter_map(|(i, val)| {

                    let x = i % width;

                    let y = i / width;

                    if x == 0 || x == width - 1 || y == 0 || y == height - 1 {

                        return None;
                    }

                    if (x + y) % 2 != pass {

                        return None;
                    }

                    // Safety: In Red-Black SOR, the neighbors we read are always from the
                    // opposite color set, which is not being modified in this pass.
                    let get_u = |nx : usize, ny : usize| unsafe {

                        let ptr = u_ptr as *const f64;

                        *ptr.add(ny * width + nx)
                    };

                    let old_val : f64 = *val;

                    let target = factor
                        * (beta.mul_add(get_u(x, y + 1) + get_u(x, y - 1), get_u(x + 1, y) + get_u(x - 1, y))
                            - dx2 * source[i]);

                    *val = omega.mul_add(target - old_val, old_val);

                    let diff : f64 = (*val - old_val).abs();

                    Some(diff)
                })
                .collect();

            if let Some(d) = diffs
                .into_iter()
                .max_by(|a, b| {

                    a.partial_cmp(b)
                        .unwrap()
                })
            {

                if d > max_diff {

                    max_diff = d;
                }
            }
        }

        if max_diff < tolerance {

            break;
        }
    }

    u
}

/// Solves 1D Burgers' equation `u_t + u*u_x = nu*u_xx`.

#[must_use] 
pub fn solve_burgers_1d(
    initial_u: &[f64],
    dx: f64,
    nu: f64,
    dt: f64,
    steps: usize,
) -> Vec<f64> {

    let mut u = initial_u.to_vec();

    let n = u.len();

    if n < 2 {

        return u;
    }

    let mut u_next = vec![0.0; n];

    for _ in 0 .. steps {

        u_next
            .par_iter_mut()
            .enumerate()
            .for_each(|(i, next_val)| {

                if i == 0 || i == n - 1 {

                    *next_val = u[i];

                    return;
                }

                // Central difference for convection and diffusion
                let convection = u[i] * (u[i + 1] - u[i - 1]) / (2.0 * dx);

                let diffusion = nu * (2.0f64.mul_add(-u[i], u[i + 1]) + u[i - 1]) / dx.powi(2);

                *next_val = dt.mul_add(diffusion - convection, u[i]);
            });

        u.copy_from_slice(&u_next);
    }

    u
}

/// Solves a 1D advection-diffusion equation `u_t + c*u_x = d*u_xx` using FDM.

#[must_use] 
pub fn solve_advection_diffusion_1d(
    initial_cond: &[f64],
    dx: f64,
    c: f64,
    d: f64,
    dt: f64,
    steps: usize,
) -> Vec<f64> {

    let mut u = initial_cond.to_vec();

    let n = u.len();

    if n < 2 {

        return u;
    }

    let mut u_next = vec![0.0; n];

    for _ in 0 .. steps {

        u_next
            .par_iter_mut()
            .enumerate()
            .for_each(|(i, next_val)| {

                if i == 0 || i == n - 1 {

                    *next_val = u[i];

                    return;
                }

                // Upwind for advection
                let advection = if c >= 0.0 {

                    -c * (u[i] - u[i - 1]) / dx
                } else {

                    -c * (u[i + 1] - u[i]) / dx
                };

                let diffusion = d * (2.0f64.mul_add(-u[i], u[i + 1]) + u[i - 1]) / dx.powi(2);

                *next_val = dt.mul_add(advection + diffusion, u[i]);
            });

        u.copy_from_slice(&u_next);
    }

    u
}

// ============================================================================
// Scenarios
// ============================================================================

/// Example scenario: Simulates heat conduction on a 100x100 plate.

#[must_use] 
pub fn simulate_2d_heat_conduction_scenario(
) -> FdmGrid<f64> {

    const WIDTH: usize = 100;

    const HEIGHT: usize = 100;

    const ALPHA: f64 = 0.02;

    const DX: f64 = 1.0;

    const DY: f64 = 1.0;

    const DT: f64 = 0.1;

    const STEPS: usize = 1000;

    solve_heat_equation_2d(
        WIDTH,
        HEIGHT,
        ALPHA,
        DX,
        DY,
        DT,
        STEPS,
        |x, y| {

            let dx_cen = x as f64
                - (WIDTH / 2) as f64;

            let dy_cen = y as f64
                - (HEIGHT / 2) as f64;

            if dy_cen.mul_add(dy_cen, dx_cen.powi(2))
                < 25.0
            {

                100.0
            } else {

                0.0
            }
        },
    )
}

/// Example scenario: Simulates wave propagation on a 120x120 grid.

#[must_use] 
pub fn simulate_2d_wave_propagation_scenario(
) -> FdmGrid<f64> {

    const WIDTH: usize = 120;

    const HEIGHT: usize = 120;

    const C: f64 = 1.0;

    const DX: f64 = 1.0;

    const DY: f64 = 1.0;

    const DT: f64 = 0.5;

    const STEPS: usize = 200;

    solve_wave_equation_2d(
        WIDTH,
        HEIGHT,
        C,
        DX,
        DY,
        DT,
        STEPS,
        |x, y| {

            let dx_cen = x as f64
                - (WIDTH / 2) as f64;

            let dy_cen = y as f64
                - (HEIGHT / 2) as f64;

            let dist2 = dy_cen.mul_add(dy_cen, dx_cen.powi(2));

            (-dist2 / 50.0).exp() // Gaussian pulse
        },
    )
}
