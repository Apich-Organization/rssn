//! # Finite Difference Method (FDM) Module
//!
//! This module provides tools for solving partial differential equations (PDEs)
//! using the finite difference method. It includes a generic grid structure
//! and a solver for the 2D heat equation as an example.
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/ripple_3d_slice_frame_02.png)
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/wave_surface_frame_00.png)
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/wave_surface_frame_01.png)
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/wave_surface_frame_02.png)
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/wave_surface_frame_03.png)
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/wave_surface_frame_04.png)
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/wave_surface_frame_05.png)
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/wave_surface_frame_06.png)
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/wave_surface_frame_07.png)
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/wave_surface_frame_08.png)
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/wave_surface_frame_09.png)
//!
//! # Coverage audit
//!
//! | Equation | 1D | 2D | 3D | Boundary conditions | Initial data |
//! |---|---|---|---|---|---|
//! | heat, explicit | no | yes | no | boundary values frozen (Dirichlet) | grid |
//! | wave, explicit | no | yes | yes | homogeneous Dirichlet | displacement, optional sampled initial velocity |
//! | Poisson | no | yes | no | Dirichlet | source grid |
//! | Burgers, advection–diffusion | yes | no | no | boundary values frozen | grid |
//!
//! The dimension-generic layer lives in [`crate::sim::structured`]: Laplacian
//! operators on 1D/2D/3D vertex-centred grids with Dirichlet, Neumann, Robin,
//! periodic and symmetry faces (per face, time-dependent data), Poisson/Helmholtz
//! with direct or Krylov solvers, θ-method heat with step-doubling error control,
//! leapfrog waves with Mur absorbing faces and Newton–Krylov nonlinear diffusion.
//! Boundary descriptions are shared through [`crate::sim::boundary`].

use std::ops::Index;
use std::ops::IndexMut;

use rayon::prelude::*;
use serde::Deserialize;
use serde::Serialize;

/// Represents the dimensions of the simulation grid.
/// cbindgen:ignore
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
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
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FdmGrid<T> {
    pub(crate) data: Vec<T>,
    pub(crate) dims: Dimensions,
}

impl<T: Clone + Default + Send + Sync> FdmGrid<T> {
    /// Creates a new grid from existing data and dimensions.
    #[must_use]
    pub const fn from_data(
        data: Vec<T>,
        dims: Dimensions,
    ) -> Self {
        Self { data, dims }
    }

    /// Creates a new grid with the given dimensions, initialized with a default value.
    #[must_use]
    pub fn new(dims: Dimensions) -> Self {
        let size = match dims {
            | Dimensions::D1(x) => x,
            | Dimensions::D2(x, y) => x * y,
            | Dimensions::D3(x, y, z) => x * y * z,
        };

        Self {
            data: vec![T::default(); size],
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
            | Dimensions::D2(x, y) => x * y,
            | Dimensions::D3(x, y, z) => x * y * z,
        };

        Self {
            data: vec![value; size],
            dims,
        }
    }

    /// Returns the dimensions of the grid.
    #[inline]
    #[must_use]
    pub const fn dimensions(&self) -> &Dimensions {
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
    pub fn as_mut_slice(&mut self) -> &mut [T] {
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

impl<T> Index<(usize, usize)> for FdmGrid<T> {
    type Output = T;

    #[inline]
    fn index(
        &self,
        (x, y): (usize, usize),
    ) -> &Self::Output {
        if let Dimensions::D2(width, _) = self.dims {
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

impl<T> IndexMut<(usize, usize)> for FdmGrid<T> {
    #[inline]
    fn index_mut(
        &mut self,
        (x, y): (usize, usize),
    ) -> &mut Self::Output {
        if let Dimensions::D2(width, _) = self.dims {
            &mut self.data[y * width + x]
        } else {
            panic!(
                "Attempted to use 2D \
                 indexing on a non-2D \
                 grid."
            );
        }
    }
}

impl<T> Index<(usize, usize, usize)> for FdmGrid<T> {
    type Output = T;

    #[inline]
    fn index(
        &self,
        (x, y, z): (usize, usize, usize),
    ) -> &Self::Output {
        if let Dimensions::D3(width, height, _) = self.dims {
            &self.data[z * width * height + y * width + x]
        } else {
            panic!(
                "Attempted to use 3D \
                 indexing on a non-3D \
                 grid."
            );
        }
    }
}

impl<T> IndexMut<(usize, usize, usize)> for FdmGrid<T> {
    #[inline]
    fn index_mut(
        &mut self,
        (x, y, z): (usize, usize, usize),
    ) -> &mut Self::Output {
        if let Dimensions::D3(width, height, _) = self.dims {
            &mut self.data[z * width * height + y * width + x]
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
// Solver Configuration
// ============================================================================

/// Configuration for 2D FDM solvers (heat, wave equations).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FdmSolverConfig2D {
    /// Grid width (number of points in x direction).
    pub width: usize,
    /// Grid height (number of points in y direction).
    pub height: usize,
    /// Grid spacing in x direction.
    pub dx: f64,
    /// Grid spacing in y direction.
    pub dy: f64,
    /// Time step.
    pub dt: f64,
    /// Number of time steps.
    pub steps: usize,
}

/// Configuration for 3D FDM solvers.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FdmSolverConfig3D {
    /// Grid width (number of points in x direction).
    pub width: usize,
    /// Grid height (number of points in y direction).
    pub height: usize,
    /// Grid depth (number of points in z direction).
    pub depth: usize,
    /// Grid spacing in x direction.
    pub dx: f64,
    /// Grid spacing in y direction.
    pub dy: f64,
    /// Grid spacing in z direction.
    pub dz: f64,
    /// Time step.
    pub dt: f64,
    /// Number of time steps.
    pub steps: usize,
}

/// Configuration for 2D Poisson solver.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PoissonSolverConfig2D {
    /// Grid width (number of points in x direction).
    pub width: usize,
    /// Grid height (number of points in y direction).
    pub height: usize,
    /// Grid spacing in x direction.
    pub dx: f64,
    /// Grid spacing in y direction.
    pub dy: f64,
    /// Relaxation factor for SOR (typically 1.0 < omega < 2.0).
    pub omega: f64,
    /// Maximum number of iterations.
    pub max_iter: usize,
    /// Convergence tolerance.
    pub tolerance: f64,
}

// ============================================================================
// Solvers
// ============================================================================

/// Solves a 2D heat equation `u_t = alpha * ∇²u` using the finite difference method.
pub fn solve_heat_equation_2d<F>(
    config: &FdmSolverConfig2D,
    alpha: f64,
    initial_conditions: F,
) -> FdmGrid<f64>
where
    F: Fn(usize, usize) -> f64 + Sync,
{
    let width = config.width;

    let height = config.height;

    let dx = config.dx;

    let dy = config.dy;

    let dt = config.dt;

    let steps = config.steps;

    let dims = Dimensions::D2(width, height);

    let mut grid = FdmGrid::new(dims);

    let mut next_grid = grid.clone();

    grid.as_mut_slice()
        .par_iter_mut()
        .enumerate()
        .for_each(|(i, val)| {
            let x = i % width;

            let y = i / width;

            *val = initial_conditions(x, y);
        });

    let r_x = alpha * dt / dx.powi(2);

    let r_y = alpha * dt / dy.powi(2);

    for _ in 0..steps {
        next_grid
            .as_mut_slice()
            .par_iter_mut()
            .enumerate()
            .for_each(|(i, next_val)| {
                let x = i % width;

                let y = i / width;

                if x == 0 || x == width - 1 || y == 0 || y == height - 1 {
                    *next_val = grid[i]; // Boundary remains constant
                    return;
                }

                let lap_x = 2.0f64.mul_add(-grid[(x, y)], grid[(x + 1, y)]) + grid[(x - 1, y)];

                let lap_y = 2.0f64.mul_add(-grid[(x, y)], grid[(x, y + 1)]) + grid[(x, y - 1)];

                *next_val = grid[i] + r_x * lap_x + r_y * lap_y;
            });

        std::mem::swap(&mut grid, &mut next_grid);
    }

    grid
}

#[derive(Clone, Copy)]
struct WaveBox<const VOLUME: bool> {
    n: [usize; 3],
    courant2: [f64; 3],
}

impl<const VOLUME: bool> WaveBox<VOLUME> {
    const fn has_interior(self) -> bool {
        self.n[0] > 2 && self.n[1] > 2 && (!VOLUME || self.n[2] > 2)
    }

    fn has_spatial_step(self) -> bool {
        self.courant2.iter().any(|&coefficient| coefficient != 0.0)
    }

    fn sample(
        self,
        field: impl Fn(usize, usize, usize) -> f64 + Sync,
    ) -> Vec<f64> {
        let [nx, ny, nz] = self.n;
        let mut values = vec![0.0; nx * ny * nz];
        if values.is_empty() {
            return values;
        }
        let rows_per_chunk = (8192 / nx).max(1);
        crate::kernels::par::for_chunks_mut(
            &mut values,
            rows_per_chunk * nx,
            4,
            |chunk, values| {
                for (local_row, row) in values.chunks_exact_mut(nx).enumerate() {
                    let index = chunk * rows_per_chunk + local_row;
                    let (y, z) = (index % ny, index / ny);
                    for (x, value) in row.iter_mut().enumerate() {
                        *value = field(x, y, z);
                    }
                }
            },
        );
        values
    }

    fn project(
        self,
        values: &mut [f64],
    ) {
        if !self.has_interior() {
            values.fill(0.0);
            return;
        }
        let [nx, ny, nz] = self.n;
        for (z, plane) in values.chunks_exact_mut(nx * ny).enumerate() {
            if VOLUME && (z == 0 || z + 1 == nz) {
                plane.fill(0.0);
                continue;
            }
            for (y, row) in plane.chunks_exact_mut(nx).enumerate() {
                if y == 0 || y + 1 == ny {
                    row.fill(0.0);
                } else {
                    row[0] = 0.0;
                    row[nx - 1] = 0.0;
                }
            }
        }
    }

    fn update_interior(
        self,
        destination: &mut [f64],
        value: impl Fn(usize, f64) -> f64 + Sync,
    ) {
        if !self.has_interior() {
            return;
        }
        let [nx, ny, nz] = self.n;
        let rows_per_chunk = (8192 / nx).max(1);
        crate::kernels::par::for_chunks_mut(
            &mut *destination,
            rows_per_chunk * nx,
            4,
            |chunk, values| {
                let rows = values.chunks_exact_mut(nx);
                for (row_index, row) in (chunk * rows_per_chunk..).zip(rows) {
                    let (y, z) = (row_index % ny, row_index / ny);
                    if y > 0 && y + 1 < ny && (!VOLUME || (z > 0 && z + 1 < nz)) {
                        let offset = row_index * nx;
                        for x in 1..nx - 1 {
                            row[x] = value(offset + x, row[x]);
                        }
                    }
                }
            },
        );
    }

    fn terms(
        self,
        u: &[f64],
        p: usize,
    ) -> [f64; 3] {
        let center = u[p];
        let laplacian = |stride| 2.0_f64.mul_add(-center, u[p + stride]) + u[p - stride];
        [
            self.courant2[0] * laplacian(1),
            self.courant2[1] * laplacian(self.n[0]),
            if VOLUME {
                self.courant2[2] * laplacian(self.n[0] * self.n[1])
            } else {
                0.0
            },
        ]
    }
}

struct DirichletWave<const VOLUME: bool> {
    grid: WaveBox<VOLUME>,
    current: Vec<f64>,
    previous: crate::kernels::reuse::Scratch<'static, f64>,
}

impl<const VOLUME: bool> DirichletWave<VOLUME> {
    fn new(
        grid: WaveBox<VOLUME>,
        current: Vec<f64>,
        velocity: Option<&[f64]>,
        dt: f64,
    ) -> Self {
        let mut previous = crate::kernels::reuse::f64_scratch(current.len());
        previous.resize(current.len(), 0.0);
        // Center the time derivative at t=0 with a half-acceleration
        // correction in the fictitious previous step.
        grid.update_interior(&mut previous, |p, _| {
            let [x, y, z] = grid.terms(&current, p);
            let half = f64::midpoint(x, y);
            let half = if VOLUME {
                half + 0.5 * z
            } else {
                half
            };
            let position = current[p] + half;
            velocity.map_or(position, |v| position - dt * v[p])
        });
        Self {
            grid,
            current,
            previous,
        }
    }

    fn advance(&mut self) {
        let current = &self.current;
        let grid = self.grid;
        // Each update needs only its own previous value, so the older
        // buffer is also the output. Boundary entries stay zero in both.
        grid.update_interior(&mut self.previous, |p, previous| {
            let [x, y, z] = grid.terms(current, p);
            let next = 2.0_f64.mul_add(current[p], -previous) + x + y;
            if VOLUME { next + z } else { next }
        });
        std::mem::swap(&mut self.current, &mut *self.previous);
    }
}

/// Solves the two-dimensional wave equation with zero initial velocity.
///
/// Positive step counts impose homogeneous Dirichlet values from the
/// initial time. With zero steps the function returns the raw initial samples.
/// Use time steps satisfying `(c dt / dx)^2 + (c dt / dy)^2 <= 1`.
pub fn solve_wave_equation_2d<F>(
    config: &FdmSolverConfig2D,
    c: f64,
    initial_u: F,
) -> FdmGrid<f64>
where
    F: Fn(usize, usize) -> f64 + Sync,
{
    let dims = Dimensions::D2(config.width, config.height);
    let grid = WaveBox::<false> {
        n: [config.width, config.height, 1],
        courant2: [
            (c * config.dt / config.dx).powi(2),
            (c * config.dt / config.dy).powi(2),
            0.0,
        ],
    };
    let initial = grid.sample(|x, y, _| initial_u(x, y));
    let data = evolve_dirichlet(grid, initial, config.dt, config.steps);
    FdmGrid::from_data(data, dims)
}

/// Solves the three-dimensional wave equation with zero initial velocity.
///
/// Uses a second-order start and homogeneous Dirichlet faces for positive
/// step counts. With zero steps the function returns the raw initial samples.
/// Use time steps satisfying `Σ_a (c dt / h_a)^2 <= 1`.
pub fn solve_wave_equation_3d<F>(
    config: &FdmSolverConfig3D,
    c: f64,
    initial_u: F,
) -> FdmGrid<f64>
where
    F: Fn(usize, usize, usize) -> f64 + Sync,
{
    let dims = Dimensions::D3(config.width, config.height, config.depth);
    let grid = WaveBox::<true> {
        n: [config.width, config.height, config.depth],
        courant2: [
            (c * config.dt / config.dx).powi(2),
            (c * config.dt / config.dy).powi(2),
            (c * config.dt / config.dz).powi(2),
        ],
    };
    let initial = grid.sample(initial_u);
    let data = evolve_dirichlet(grid, initial, config.dt, config.steps);
    FdmGrid::from_data(data, dims)
}

fn evolve_dirichlet<const VOLUME: bool>(
    grid: WaveBox<VOLUME>,
    mut data: Vec<f64>,
    dt: f64,
    steps: usize,
) -> Vec<f64> {
    if steps == 0 {
        return data;
    }
    grid.project(&mut data);
    if !grid.has_interior() || !grid.has_spatial_step() {
        return data;
    }
    let mut solver = DirichletWave::new(grid, data, None, dt);
    for _ in 0..steps {
        solver.advance();
    }
    solver.current
}

/// One retained displacement field from the homogeneous Dirichlet wave solver.
#[derive(Clone, Debug)]
pub struct DirichletWaveSample {
    /// Number of completed steps since the initial time.
    pub step: usize,
    /// Displacement values at the grid nodes, x fastest.
    pub displacement: FdmGrid<f64>,
}

/// Evolves a two- or three-dimensional Dirichlet wave and retains step samples.
///
/// Initial displacement and velocity share a grid. Faces are zero from the
/// initial sample onward. Retains step zero, every multiple of `save_every`,
/// and the last step. Time advances by the signed step `dt`.
///
/// # Errors
/// Invalid shapes, spacing, speed, time step or cadence, a violated leapfrog
/// CFL bound, non-finite data or states, allocation limits, or cancellation.
pub fn solve_wave_dirichlet_sampled(
    u0: &FdmGrid<f64>,
    v0: &FdmGrid<f64>,
    spacing: &[f64],
    c: f64,
    dt: f64,
    steps: usize,
    save_every: usize,
) -> Result<Vec<DirichletWaveSample>, String> {
    let (n, dim) = match *u0.dimensions() {
        | Dimensions::D2(x, y) => ([x, y, 1], 2),
        | Dimensions::D3(x, y, z) => ([x, y, z], 3),
        | Dimensions::D1(_) => {
            return Err("sampled Dirichlet waves require two or three axes".into());
        },
    };
    let point_count = n
        .iter()
        .try_fold(1usize, |size, extent| size.checked_mul(*extent));
    if n.contains(&0)
        || point_count != Some(u0.len())
        || v0.dimensions() != u0.dimensions()
        || v0.len() != u0.len()
        || spacing.len() != dim
        || spacing.iter().any(|&dx| !dx.is_finite() || dx <= 0.0)
        || !c.is_finite()
        || c < 0.0
        || !dt.is_finite()
        || u0
            .as_slice()
            .iter()
            .chain(v0.as_slice())
            .any(|value| !value.is_finite())
    {
        return Err("invalid sampled Dirichlet wave data or mesh".into());
    }
    let stride =
        std::num::NonZeroUsize::new(save_every).ok_or("wave sample cadence must be positive")?;
    let count = crate::kernels::numeric_common::trajectory_sample_count(steps, stride)
        .ok_or("too many sampled wave states")?;
    let mut samples = Vec::new();
    samples
        .try_reserve_exact(count)
        .map_err(|_| "too many sampled wave states")?;
    let mut courant2 = [0.0; 3];
    if c > 0.0 && dt != 0.0 {
        let travel = c * dt.abs();
        for (ratio2, &dx) in courant2.iter_mut().zip(spacing) {
            let ratio = if travel.is_normal() {
                travel / dx
            } else {
                crate::kernels::numeric_common::scale_product(c, &[dt.abs()], &[dx])
                    .map_err(|error| error.to_string())?
            };
            *ratio2 = ratio * ratio;
        }
    }
    if courant2.iter().sum::<f64>() > 1.0 {
        return Err("sampled Dirichlet wave CFL bound exceeded".into());
    }
    let mut retain = |step, values: &[f64]| {
        samples.push(DirichletWaveSample {
            step,
            displacement: FdmGrid::from_data(values.to_vec(), u0.dimensions().clone()),
        });
    };
    if dim == 2 {
        sample_dirichlet(
            WaveBox::<false> { n, courant2 },
            u0,
            v0,
            dt,
            steps,
            stride,
            &mut retain,
        )?;
    } else {
        sample_dirichlet(
            WaveBox::<true> { n, courant2 },
            u0,
            v0,
            dt,
            steps,
            stride,
            &mut retain,
        )?;
    }
    Ok(samples)
}

fn sample_dirichlet<const VOLUME: bool>(
    grid: WaveBox<VOLUME>,
    u0: &FdmGrid<f64>,
    v0: &FdmGrid<f64>,
    dt: f64,
    steps: usize,
    stride: std::num::NonZeroUsize,
    retain: &mut impl FnMut(usize, &[f64]),
) -> Result<(), String> {
    crate::cancellation::checkpoint().map_err(|error| error.to_string())?;
    let mut initial = u0.as_slice().to_vec();
    grid.project(&mut initial);
    retain(0, &initial);
    if steps == 0 {
        return Ok(());
    }
    let moving = dt != 0.0 && grid.has_interior();
    if !moving || !grid.has_spatial_step() {
        let mut completed = 0;
        while completed < steps {
            crate::cancellation::checkpoint().map_err(|error| error.to_string())?;
            completed = completed.saturating_add(stride.get()).min(steps);
            if moving {
                let time = completed as f64 * dt;
                grid.update_interior(&mut initial, |p, _| {
                    let rate = v0[p];
                    if time.is_finite() {
                        time.mul_add(rate, u0[p])
                    } else if rate == 0.0 {
                        u0[p]
                    } else {
                        let displacement = crate::kernels::numeric_common::scale_product(
                            rate,
                            &[completed as f64, dt.abs()],
                            &[],
                        )
                        .unwrap_or(f64::NAN);
                        u0[p] + f64::copysign(1.0, dt) * displacement
                    }
                });
                if initial.iter().any(|value| !value.is_finite()) {
                    return Err("sampled Dirichlet wave state is not representable".into());
                }
            }
            retain(completed, &initial);
        }
        return crate::cancellation::checkpoint().map_err(|error| error.to_string());
    }
    let mut solver = DirichletWave::new(grid, initial, Some(v0.as_slice()), dt);
    let mut completed = 0;
    while completed < steps {
        let next = completed.saturating_add(stride.get()).min(steps);
        for _ in completed..next {
            crate::cancellation::checkpoint().map_err(|error| error.to_string())?;
            solver.advance();
        }
        if solver.current.iter().any(|value| !value.is_finite()) {
            return Err("sampled Dirichlet wave state is not representable".into());
        }
        completed = next;
        retain(completed, &solver.current);
    }
    crate::cancellation::checkpoint().map_err(|error| error.to_string())
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
    config: &PoissonSolverConfig2D,
    source: &FdmGrid<f64>,
) -> FdmGrid<f64> {
    let width = config.width;

    let height = config.height;

    let dx = config.dx;

    let dy = config.dy;

    let omega = config.omega;

    let max_iter = config.max_iter;

    let tolerance = config.tolerance;

    let dims = Dimensions::D2(width, height);

    let mut u: FdmGrid<f64> = FdmGrid::new(dims);

    let dx2 = dx * dx;

    let dy2 = dy * dy;

    let beta = dx2 / dy2;

    let factor = 1.0 / (2.0 * (1.0 + beta));

    for _ in 0..max_iter {
        let mut max_diff: f64 = 0.0;

        for pass in 0..2 {
            // Red-black SOR: neighbours read here belong to the opposite colour, which is
            // untouched in this pass, so computing all updates from a shared view and then
            // applying them is exactly the in-place sweep.
            let u_view = u.as_slice();

            let updates: Vec<(usize, f64, f64)> = u_view
                .par_iter()
                .enumerate()
                .filter_map(|(i, &old_val)| {
                    let x = i % width;

                    let y = i / width;

                    if x == 0 || x == width - 1 || y == 0 || y == height - 1 {
                        return None;
                    }

                    if (x + y) % 2 != pass {
                        return None;
                    }

                    let get_u = |nx: usize, ny: usize| u_view[ny * width + nx];

                    let target = factor
                        * dx2.mul_add(
                            -source[i],
                            beta.mul_add(
                                get_u(x, y + 1) + get_u(x, y - 1),
                                get_u(x + 1, y) + get_u(x - 1, y),
                            ),
                        );

                    let new_val = omega.mul_add(target - old_val, old_val);

                    Some((i, new_val, (new_val - old_val).abs()))
                })
                .collect();

            let u_mut = u.as_mut_slice();

            for &(i, new_val, _) in &updates {
                u_mut[i] = new_val;
            }

            if let Some(d) = updates
                .iter()
                .map(|&(_, _, d)| d)
                .max_by(|a, b| a.partial_cmp(b).unwrap())
                && d > max_diff
            {
                max_diff = d;
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

    for _ in 0..steps {
        u_next.par_iter_mut().enumerate().for_each(|(i, next_val)| {
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

    for _ in 0..steps {
        u_next.par_iter_mut().enumerate().for_each(|(i, next_val)| {
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
pub fn simulate_2d_heat_conduction_scenario() -> FdmGrid<f64> {
    const WIDTH: usize = 100;

    const HEIGHT: usize = 100;

    const ALPHA: f64 = 0.02;

    let config = FdmSolverConfig2D {
        width: WIDTH,
        height: HEIGHT,
        dx: 1.0,
        dy: 1.0,
        dt: 0.1,
        steps: 1000,
    };

    solve_heat_equation_2d(&config, ALPHA, |x, y| {
        let dx_cen = x as f64 - (WIDTH / 2) as f64;

        let dy_cen = y as f64 - (HEIGHT / 2) as f64;

        if dy_cen.mul_add(dy_cen, dx_cen.powi(2)) < 25.0 {
            100.0
        } else {
            0.0
        }
    })
}

/// Example scenario: Simulates wave propagation on a 120x120 grid.
#[must_use]
pub fn simulate_2d_wave_propagation_scenario() -> FdmGrid<f64> {
    const WIDTH: usize = 120;

    const HEIGHT: usize = 120;

    const C: f64 = 1.0;

    let config = FdmSolverConfig2D {
        width: WIDTH,
        height: HEIGHT,
        dx: 1.0,
        dy: 1.0,
        dt: 0.5,
        steps: 200,
    };

    solve_wave_equation_2d(&config, C, |x, y| {
        let dx_cen = x as f64 - (WIDTH / 2) as f64;

        let dy_cen = y as f64 - (HEIGHT / 2) as f64;

        let dist2 = dy_cen.mul_add(dy_cen, dx_cen.powi(2));

        (-dist2 / 50.0).exp() // Gaussian pulse
    })
}

fn dims_of(d: &Dimensions) -> Vec<usize> {
    match *d {
        | Dimensions::D1(x) => vec![x],
        | Dimensions::D2(x, y) => vec![x, y],
        | Dimensions::D3(x, y, z) => vec![x, y, z],
    }
}

fn structured_grid(
    dims: &[usize],
    spacing: &[f64],
    bcs: &crate::sim::boundary::FaceBcs,
) -> Result<crate::sim::structured::Grid, String> {
    if dims.len() != spacing.len() {
        return Err("one spacing per axis is required".into());
    }
    let hi: Vec<f64> = (0..dims.len())
        .map(|a| {
            spacing[a]
                * if bcs.axis_periodic(a).unwrap_or(false) {
                    dims[a]
                } else {
                    dims[a] - 1
                } as f64
        })
        .collect();
    crate::sim::structured::Grid::new(dims, &vec![0.0; dims.len()], &hi, bcs)
}

fn node_value(
    grid: &crate::sim::structured::Grid,
    data: &[f64],
    spacing: &[f64],
    x: [f64; 3],
) -> f64 {
    let mut idx = [0usize; 3];
    for a in 0..grid.dim {
        idx[a] = ((x[a] / spacing[a]).round() as usize).min(grid.n[a] - 1);
    }
    data[grid.index(idx)]
}

/// Heat equation `u_t = α Δu` on the grid of `initial` (1D, 2D or 3D).
///
/// Face conditions `bcs` may be Dirichlet, Neumann, Robin, periodic or time
/// dependent; the Crank–Nicolson method of [`crate::sim::structured`] is used.
///
/// # Errors
/// For inconsistent sizes or boundary descriptions.
pub fn solve_heat_equation_bc(
    initial: &FdmGrid<f64>,
    spacing: &[f64],
    alpha: f64,
    dt: f64,
    steps: usize,
    bcs: &crate::sim::boundary::FaceBcs,
) -> Result<FdmGrid<f64>, String> {
    use crate::sim::structured::HeatProblem;
    use crate::sim::structured::LinearSolver;
    use crate::sim::structured::heat_fixed;
    let dims = dims_of(initial.dimensions());
    let grid = structured_grid(&dims, spacing, bcs)?;
    let init = |x: [f64; 3]| node_value(&grid, initial.as_slice(), spacing, x);
    let prob = HeatProblem {
        grid: &grid,
        bcs,
        alpha,
        init: &init,
        source: None,
        solver: LinearSolver::Auto,
    };
    let u = heat_fixed(&prob, dt * steps as f64, steps, 0.5)?;
    Ok(FdmGrid::from_data(u, initial.dimensions().clone()))
}

/// Wave equation `u_tt = c² Δu` on a 1D, 2D or 3D grid with face conditions
/// `bcs` (including first-order absorbing faces) by leapfrog.
///
/// # Errors
/// If the CFL condition is violated or the description is inconsistent.
pub fn solve_wave_equation_bc(
    u0: &FdmGrid<f64>,
    v0: &FdmGrid<f64>,
    spacing: &[f64],
    c: f64,
    dt: f64,
    steps: usize,
    bcs: &crate::sim::boundary::FaceBcs,
) -> Result<FdmGrid<f64>, String> {
    use crate::sim::structured::WaveProblem;
    use crate::sim::structured::wave_leapfrog;
    let dims = dims_of(u0.dimensions());
    let grid = structured_grid(&dims, spacing, bcs)?;
    let (f0, f1) = (
        |x: [f64; 3]| node_value(&grid, u0.as_slice(), spacing, x),
        |x: [f64; 3]| node_value(&grid, v0.as_slice(), spacing, x),
    );
    let prob = WaveProblem {
        grid: &grid,
        bcs,
        c,
        u0: &f0,
        v0: &f1,
        source: None,
    };
    let u = wave_leapfrog(&prob, dt * steps as f64, steps)?;
    Ok(FdmGrid::from_data(u, u0.dimensions().clone()))
}

/// Poisson/Helmholtz equation `-Δu + σ u = f` on a 1D, 2D or 3D grid with
/// face conditions `bcs`, solved by geometric multigrid.
///
/// # Errors
/// For inconsistent sizes or boundary descriptions.
pub fn solve_poisson_bc(
    source: &FdmGrid<f64>,
    spacing: &[f64],
    sigma: f64,
    bcs: &crate::sim::boundary::FaceBcs,
) -> Result<FdmGrid<f64>, String> {
    let dims = dims_of(source.dimensions());
    let r = crate::sim::physics_mtm::solve_poisson_multigrid_nd(
        &dims,
        spacing,
        source.as_slice(),
        sigma,
        bcs,
        crate::sim::multigrid::MgOptions::default(),
    )?;
    Ok(FdmGrid::from_data(r.u, source.dimensions().clone()))
}

// ============================================================================
// Space-dimension coverage: Poisson 1D/3D, Burgers 2D/3D, advection-diffusion 2D/3D
// ============================================================================
//
// Kernels: [`crate::sim::dimensional`]. Arrays are flat with `x` fastest and the
// boundary values are frozen (Poisson: zero Dirichlet).

/// Configuration of the dimension-generic SOR Poisson solver
/// ([`solve_poisson_nd`]); the 1D/2D/3D special cases are
/// [`solve_poisson_1d`], [`solve_poisson_2d`] and [`solve_poisson_3d`].
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PoissonSolverConfigNd {
    /// Nodes per axis (1 to 3 entries).
    pub dims: Vec<usize>,
    /// Grid spacing per axis.
    pub spacing: Vec<f64>,
    /// SOR relaxation factor.
    pub omega: f64,
    /// Maximum number of sweeps.
    pub max_iter: usize,
    /// Convergence tolerance on the largest update.
    pub tolerance: f64,
}

/// Solves `Δu = f` with zero Dirichlet boundaries on a 1D, 2D or 3D grid by SOR.
///
/// # Errors
/// For inconsistent shapes or lengths.
pub fn solve_poisson_nd(
    config: &PoissonSolverConfigNd,
    source: &[f64],
) -> Result<Vec<f64>, String> {
    let zero = vec![0.0; source.len()];

    crate::sim::dimensional::solve_poisson_iterative_nd(
        source,
        &zero,
        &config.dims,
        &config.spacing,
        crate::sim::dimensional::IterativeMethod::Sor(config.omega),
        config.max_iter,
        config.tolerance,
    )
}

/// 1D Poisson equation `u'' = f`, zero Dirichlet ends, SOR. 1D counterpart of
/// [`solve_poisson_2d`].
///
/// # Errors
/// As [`solve_poisson_nd`].
pub fn solve_poisson_1d(
    source: &[f64],
    dx: f64,
    omega: f64,
    max_iter: usize,
    tolerance: f64,
) -> Result<Vec<f64>, String> {
    solve_poisson_nd(
        &PoissonSolverConfigNd {
            dims: vec![source.len()],
            spacing: vec![dx],
            omega,
            max_iter,
            tolerance,
        },
        source,
    )
}

/// 3D Poisson equation `Δu = f`, zero Dirichlet boundaries, SOR. 3D counterpart
/// of [`solve_poisson_2d`].
///
/// # Errors
/// As [`solve_poisson_nd`].
pub fn solve_poisson_3d(
    source: &[f64],
    dims: [usize; 3],
    h: [f64; 3],
    omega: f64,
    max_iter: usize,
    tolerance: f64,
) -> Result<Vec<f64>, String> {
    solve_poisson_nd(
        &PoissonSolverConfigNd {
            dims: dims.to_vec(),
            spacing: h.to_vec(),
            omega,
            max_iter,
            tolerance,
        },
        source,
    )
}

/// 2D Burgers' equation `u_t + u (u_x + u_y) = nu Δu` with central convection
/// differences and frozen boundaries; returns the state after `steps` steps.
/// 2D counterpart of [`solve_burgers_1d`].
///
/// # Errors
/// For inconsistent shapes or a non-positive `dt`.
pub fn solve_burgers_2d(
    initial_u: &[f64],
    dims: [usize; 2],
    h: [f64; 2],
    nu: f64,
    dt: f64,
    steps: usize,
) -> Result<Vec<f64>, String> {
    crate::sim::dimensional::advance_transport_nd(
        initial_u,
        &dims,
        &h,
        &crate::sim::dimensional::Transport::BurgersCentral { nu },
        dt,
        steps,
    )
}

/// 3D Burgers' equation; see [`solve_burgers_2d`].
///
/// # Errors
/// As [`solve_burgers_2d`].
pub fn solve_burgers_3d(
    initial_u: &[f64],
    dims: [usize; 3],
    h: [f64; 3],
    nu: f64,
    dt: f64,
    steps: usize,
) -> Result<Vec<f64>, String> {
    crate::sim::dimensional::advance_transport_nd(
        initial_u,
        &dims,
        &h,
        &crate::sim::dimensional::Transport::BurgersCentral { nu },
        dt,
        steps,
    )
}

/// 2D advection-diffusion `u_t + c·∇u = d Δu` (upwind + central, frozen
/// boundaries); returns the final state. 2D counterpart of
/// [`solve_advection_diffusion_1d`].
///
/// # Errors
/// As [`solve_burgers_2d`].
pub fn solve_advection_diffusion_2d(
    initial_cond: &[f64],
    dims: [usize; 2],
    h: [f64; 2],
    c: [f64; 2],
    d: f64,
    dt: f64,
    steps: usize,
) -> Result<Vec<f64>, String> {
    crate::sim::dimensional::advance_transport_nd(
        initial_cond,
        &dims,
        &h,
        &crate::sim::dimensional::Transport::AdvectionDiffusion {
            velocity: c.to_vec(),
            alpha: d,
        },
        dt,
        steps,
    )
}

/// 3D advection-diffusion; see [`solve_advection_diffusion_2d`].
///
/// # Errors
/// As [`solve_burgers_2d`].
pub fn solve_advection_diffusion_3d(
    initial_cond: &[f64],
    dims: [usize; 3],
    h: [f64; 3],
    c: [f64; 3],
    d: f64,
    dt: f64,
    steps: usize,
) -> Result<Vec<f64>, String> {
    crate::sim::dimensional::advance_transport_nd(
        initial_cond,
        &dims,
        &h,
        &crate::sim::dimensional::Transport::AdvectionDiffusion {
            velocity: c.to_vec(),
            alpha: d,
        },
        dt,
        steps,
    )
}
