//! # Finite Difference Method (FDM) Module
//!
//! This module provides tools for solving partial differential equations (PDEs)
//! using the finite difference method. It includes a generic grid structure
//! and a solver for the 2D heat equation as an example.

use rayon::prelude::*;
use std::ops::{Index, IndexMut};

/// Represents the dimensions of the simulation grid.
#[derive(Clone)]
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
#[derive(Clone)]
pub struct Grid<T> {
    data: Vec<T>,
    dims: Dimensions,
}

impl<T: Clone + Default + Send + Sync> Grid<T> {
    /// Creates a new grid with the given dimensions, initialized with a default value.
    pub fn new(dims: Dimensions) -> Self {
        let size = match dims {
            Dimensions::D1(x) => x,
            Dimensions::D2(x, y) => x * y,
            Dimensions::D3(x, y, z) => x * y * z,
        };
        Grid {
            data: vec![T::default(); size],
            dims,
        }
    }

    /// Returns the dimensions of the grid.
    #[inline]
    pub fn dimensions(&self) -> &Dimensions {
        &self.dims
    }

    /// Returns a slice to the underlying data.
    #[inline]
    pub fn as_slice(&self) -> &[T] {
        &self.data
    }

    /// Returns a mutable slice to the underlying data.
    #[inline]
    pub fn as_mut_slice(&mut self) -> &mut [T] {
        &mut self.data
    }
}

// Indexing for 1D grids
impl<T> Index<usize> for Grid<T> {
    type Output = T;
    #[inline]
    fn index(&self, index: usize) -> &Self::Output {
        &self.data[index]
    }
}

impl<T> IndexMut<usize> for Grid<T> {
    #[inline]
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.data[index]
    }
}

// Indexing for 2D grids
impl<T> Index<(usize, usize)> for Grid<T> {
    type Output = T;
    #[inline]
    fn index(&self, (x, y): (usize, usize)) -> &Self::Output {
        if let Dimensions::D2(width, _) = self.dims {
            &self.data[y * width + x]
        } else {
            panic!("Attempted to use 2D indexing on a non-2D grid.");
        }
    }
}

impl<T> IndexMut<(usize, usize)> for Grid<T> {
    #[inline]
    fn index_mut(&mut self, (x, y): (usize, usize)) -> &mut Self::Output {
        if let Dimensions::D2(width, _) = self.dims {
            &mut self.data[y * width + x]
        } else {
            panic!("Attempted to use 2D indexing on a non-2D grid.");
        }
    }
}

// Indexing for 3D grids
impl<T> Index<(usize, usize, usize)> for Grid<T> {
    type Output = T;
    #[inline]
    fn index(&self, (x, y, z): (usize, usize, usize)) -> &Self::Output {
        if let Dimensions::D3(width, height, _) = self.dims {
            &self.data[z * width * height + y * width + x]
        } else {
            panic!("Attempted to use 3D indexing on a non-3D grid.");
        }
    }
}

impl<T> IndexMut<(usize, usize, usize)> for Grid<T> {
    #[inline]
    fn index_mut(&mut self, (x, y, z): (usize, usize, usize)) -> &mut Self::Output {
        if let Dimensions::D3(width, height, _) = self.dims {
            &mut self.data[z * width * height + y * width + x]
        } else {
            panic!("Attempted to use 3D indexing on a non-3D grid.");
        }
    }
}

/// Solves a 2D heat equation `u_t = alpha * ∇²u` using the finite difference method.
///
/// This function simulates heat conduction on a 2D plate with fixed zero-temperature boundaries.
/// It uses a simple forward-Euler in time, central-difference in space scheme.
///
/// # Arguments
/// * `width` - The width of the grid.
/// * `height` - The height of the grid.
/// * `alpha` - The thermal diffusivity constant.
/// * `dt` - The time step. For stability, `dt` must be small, typically `dt <= dx² / (4 * alpha)`.
/// * `steps` - The number of time steps to simulate.
/// * `initial_conditions` - A function `(x, y) -> temperature` that sets the initial temperature for each grid cell.
///
/// # Returns
/// A `Grid<f64>` containing the final temperature distribution after all time steps.
pub fn solve_heat_equation_2d<F>(
    width: usize,
    height: usize,
    alpha: f64,
    dt: f64,
    steps: usize,
    initial_conditions: F,
) -> Grid<f64>
where
    F: Fn(usize, usize) -> f64 + Sync,
{
    let dims = Dimensions::D2(width, height);
    let mut grid = Grid::new(dims);
    let mut next_grid = grid.clone();

    // Apply initial conditions in parallel
    grid.as_mut_slice()
        .par_iter_mut()
        .enumerate()
        .for_each(|(i, val)| {
            let x = i % width;
            let y = i / width;
            *val = initial_conditions(x, y);
        });

    // Main simulation loop
    for _ in 0..steps {
        next_grid
            .as_mut_slice()
            .par_iter_mut()
            .enumerate()
            .for_each(|(i, next_val)| {
                let x = i % width;
                let y = i / width;

                // Apply fixed boundary conditions (temperature = 0 at edges)
                if x == 0 || x == width - 1 || y == 0 || y == height - 1 {
                    *next_val = 0.0;
                    return;
                }

                // Central difference scheme for the Laplacian
                let laplacian =
                    grid[(x + 1, y)] + grid[(x - 1, y)] + grid[(x, y + 1)] + grid[(x, y - 1)]
                        - 4.0 * grid[(x, y)];

                // Forward Euler method for time evolution
                // This is a simple but potentially unstable method. A smaller dt is better.
                // dx and dy are assumed to be 1 for simplicity.
                *next_val = grid[(x, y)] + alpha * dt * laplacian;
            });

        std::mem::swap(&mut grid, &mut next_grid);
    }

    grid
}

/// Example scenario: Simulates heat conduction on a 100x100 plate
/// with a central heat source. This function demonstrates a concrete application
/// of the FDM solver.
///
/// # Returns
/// A `Grid<f64>` containing the final temperature distribution.
pub fn simulate_2d_heat_conduction_scenario() -> Grid<f64> {
    const WIDTH: usize = 100;
    const HEIGHT: usize = 100;
    const ALPHA: f64 = 0.02; // Thermal diffusivity
    const DT: f64 = 0.1; // Time step, chosen for stability
    const STEPS: usize = 1000; // Number of simulation steps

    solve_heat_equation_2d(WIDTH, HEIGHT, ALPHA, DT, STEPS, |x, y| {
        // Initial condition: a hot spot in the center
        let dx = x as f64 - (WIDTH / 2) as f64;
        let dy = y as f64 - (HEIGHT / 2) as f64;
        if dx.powi(2) + dy.powi(2) < 25.0 {
            // Radius 5
            100.0
        } else {
            0.0
        }
    })
}
