//! Multigrid methods.
//!
//! # Multigrid (Poisson)
//!
//! | Equation | 1D | 2D | 3D | Boundary conditions | Cycles |
//! |---|---|---|---|---|---|
//! | Poisson | yes | yes | no | zero Dirichlet | V |
//!
//! [`crate::sim::multigrid`] provides geometric multigrid for
//! `-Δu + σu = f` in 1D/2D/3D with Dirichlet, Neumann, Robin and periodic faces,
//! V, W and full multigrid cycles.

use rayon::prelude::*;
use serde::Deserialize;
use serde::Serialize;

/// A 1D grid for the multigrid solver.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Grid {
    /// The solution vector.
    pub u: Vec<f64>,
    /// The right-hand side vector.
    pub f: Vec<f64>,
    /// The grid spacing.
    pub h: f64,
}

impl Grid {
    pub(crate) fn new(
        size: usize,
        h: f64,
    ) -> Self {
        Self {
            u: vec![0.0; size],
            f: vec![0.0; size],
            h,
        }
    }

    pub(crate) const fn size(&self) -> usize {
        self.u.len()
    }
}

/// Computes the residual r = f - Au for the 1D Poisson problem.
pub(crate) fn calculate_residual(grid: &Grid) -> Vec<f64> {
    let n = grid.size();

    let h_sq_inv = 1.0 / (grid.h * grid.h);

    let mut residual = vec![0.0; n];

    for (i, vars) in residual.iter_mut().enumerate().take(n - 1).skip(1) {
        let a_u = (2.0f64.mul_add(grid.u[i], -grid.u[i - 1]) - grid.u[i + 1]) * h_sq_inv;

        *vars = grid.f[i] - a_u;
    }

    residual
}

/// Smoother: Applies a few iterations of the weighted Jacobi method.
pub(crate) fn smooth(
    grid: &mut Grid,
    num_sweeps: usize,
) {
    let n = grid.size();

    let h_sq = grid.h * grid.h;

    let omega = 2.0 / 3.0;

    for _ in 0..num_sweeps {
        let u_old = grid.u.clone();

        for i in 1..n - 1 {
            let prev = u_old[i - 1];

            let next = u_old[i + 1];

            let f_i = grid.f[i];

            let new_u_i = 0.5 * h_sq.mul_add(f_i, prev + next);

            grid.u[i] = (1.0_f64 - omega).mul_add(u_old[i], omega * new_u_i);
        }
    }
}

/// Restriction: Transfers a fine-grid residual to a coarse grid using full weighting.
pub(crate) fn restrict(fine_residual: &[f64]) -> Vec<f64> {
    let fine_n = fine_residual.len();

    let coarse_n = (fine_n / 2) + 1;

    let mut coarse_f = vec![0.0; coarse_n];

    for (i, vars) in coarse_f.iter_mut().enumerate().take(coarse_n - 1).skip(1) {
        let j = 2 * i;

        *vars = 0.25f64.mul_add(
            fine_residual[j + 1],
            0.25f64.mul_add(fine_residual[j - 1], 0.5 * fine_residual[j]),
        );
    }

    coarse_f
}

/// Prolongation: Interpolates a coarse-grid correction to a fine grid.
pub(crate) fn prolongate(coarse_correction: &[f64]) -> Vec<f64> {
    let coarse_n = coarse_correction.len();

    let fine_n = 2 * (coarse_n - 1) + 1;

    let mut fine_correction = vec![0.0; fine_n];

    for i in 0..coarse_n {
        fine_correction[2 * i] = coarse_correction[i];
    }

    for i in 0..coarse_n - 1 {
        fine_correction[2 * i + 1] = f64::midpoint(coarse_correction[i], coarse_correction[i + 1]);
    }

    fine_correction
}

/// Performs a single multigrid V-cycle.
pub(crate) fn v_cycle(
    grid: &mut Grid,
    level: usize,
    max_levels: usize,
) {
    let pre_sweeps = 2;

    let post_sweeps = 2;

    smooth(grid, pre_sweeps);

    if level < max_levels - 1 {
        let residual = calculate_residual(grid);

        let coarse_f = restrict(&residual);

        let coarse_n = coarse_f.len();

        let mut coarse_grid = Grid::new(coarse_n, grid.h * 2.0);

        coarse_grid.f = coarse_f;

        v_cycle(&mut coarse_grid, level + 1, max_levels);

        let correction = prolongate(&coarse_grid.u);

        let n = grid.size();

        for (u_i, &corr_i) in grid.u.iter_mut().zip(&correction).take(n) {
            *u_i += corr_i;
        }
    }

    smooth(grid, post_sweeps);
}

/// Solves the 1D Poisson equation `-u_xx = f` using the multigrid method.
///
/// # Arguments
/// * `n` - Number of interior points on the finest grid. Must be `2^k - 1`.
/// * `f` - The right-hand side function values on the finest grid.
/// * `num_cycles` - The number of V-cycles to perform.
///
/// # Returns
/// The solution vector `u`.
///
/// # Errors
///
/// This function will return an error if the grid size `n` is not of the form `2^k - 1`.
pub fn solve_poisson_1d_multigrid(
    n: usize,
    f: &[f64],
    num_cycles: usize,
) -> Result<Vec<f64>, String> {
    let num_levels: usize = ((n as f64 + 1.0).log2() as i64).try_into().unwrap_or(0);

    if (2_usize.pow(num_levels as u32) - 1) != n {
        return Err("Grid size `n` \
                    must be of the \
                    form 2^k - 1."
            .to_string());
    }

    let mut finest_grid = Grid::new(n + 2, 1.0 / (n + 1) as f64);

    finest_grid.f[1..=n].copy_from_slice(f);

    for _ in 0..num_cycles {
        v_cycle(&mut finest_grid, 0, num_levels);
    }

    Ok(finest_grid.u)
}

/// Solves a 1D Poisson problem with a known analytical solution.
/// `-u_xx = 2` on `[0, 1]` with `u(0)=u(1)=0`. Exact solution is `u(x) = x(1-x)`.
///
/// # Errors
///
/// This function will return an error if the underlying `solve_poisson_1d_multigrid` function
/// encounters an error, e.g., if the grid size is invalid.
pub fn simulate_1d_poisson_multigrid_scenario() -> Result<Vec<f64>, String> {
    const K: usize = 7;

    const N_INTERIOR: usize = 2_usize.pow(K as u32) - 1;

    let f = vec![2.0; N_INTERIOR];

    let num_v_cycles = 10;

    solve_poisson_1d_multigrid(N_INTERIOR, &f, num_v_cycles)
}

/// A 2D grid for the multigrid solver.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Grid2D {
    /// The solution vector.
    pub u: Vec<f64>,
    /// The right-hand side vector.
    pub f: Vec<f64>,
    /// The number of grid points in each dimension.
    pub n: usize,
    /// The grid spacing.
    pub h: f64,
}

impl Grid2D {
    fn new(
        n: usize,
        h: f64,
    ) -> Self {
        Self {
            u: vec![0.0; n * n],
            f: vec![0.0; n * n],
            n,
            h,
        }
    }
}

// The two independent row passes can split large grids in the active pool.
const PARALLEL_POISSON_SIDE: usize = 257;
const POISSON_ROW_BATCH: usize = 16;

fn smooth_2d(grid: &mut Grid2D) {
    let n = grid.n;
    let h_sq = grid.h * grid.h;

    // A pre- or post-smoothing pass has two red-black sweeps. Nodes of one
    // colour read only nodes of the opposite colour, so each update is in place.
    for _ in 0..2 {
        for colour in 0..2 {
            for i in 1..n - 1 {
                let first = 1 + ((i + 1 + colour) % 2);
                for j in (first..n - 1).step_by(2) {
                    let k = i * n + j;
                    let neighbors = grid.u[k - n] + grid.u[k + n] + grid.u[k - 1] + grid.u[k + 1];
                    grid.u[k] = 0.25 * h_sq.mul_add(grid.f[k], neighbors);
                }
            }
        }
    }
}

fn residual_2d(
    grid: &Grid2D,
    residual: &mut [f64],
    parallel: bool,
) {
    let n = grid.n;
    let h_sq_inv = 1.0 / (grid.h * grid.h);
    let interior = &mut residual[n..n * (n - 1)];
    let residual_row = |(i, row): (usize, &mut [f64])| {
        let first = (i + 1) * n;
        for (j, val) in row.iter_mut().enumerate().take(n - 1).skip(1) {
            let k = first + j;
            let a_u = (4.0f64.mul_add(grid.u[k], -grid.u[k - n])
                - grid.u[k + n]
                - grid.u[k - 1]
                - grid.u[k + 1])
                * h_sq_inv;
            *val = grid.f[k] - a_u;
        }
    };

    if parallel {
        interior
            .par_chunks_mut(n)
            .with_min_len(POISSON_ROW_BATCH)
            .enumerate()
            .for_each(residual_row);
    } else {
        interior.chunks_mut(n).enumerate().for_each(residual_row);
    }
}

fn restrict_2d(
    residual: &[f64],
    fine_n: usize,
    coarse: &mut Grid2D,
) {
    for i in 1..coarse.n - 1 {
        for j in 1..coarse.n - 1 {
            let k = 2 * (i * fine_n + j);
            coarse.f[i * coarse.n + j] = (residual[k].mul_add(
                4.0,
                (residual[k - fine_n] + residual[k + fine_n] + residual[k - 1] + residual[k + 1])
                    * 2.0,
            ) + (residual[k - fine_n - 1]
                + residual[k + fine_n - 1]
                + residual[k - fine_n + 1]
                + residual[k + fine_n + 1]))
                / 16.0;
        }
    }
}

// Preserve the interpolation grouping and rounding of the fixed V-cycle.
#[allow(clippy::manual_midpoint)]
fn correct_2d(
    fine: &mut Grid2D,
    coarse: &Grid2D,
    parallel: bool,
) {
    let n = fine.n;
    let w = coarse.n;
    let interior = &mut fine.u[n..n * (n - 1)];
    let correct_row = |(i, row): (usize, &mut [f64])| {
        let i = i + 1;
        let base = (i / 2) * w;
        for (j, value) in row.iter_mut().enumerate().take(n - 1).skip(1) {
            let k = base + j / 2;
            *value += match (i % 2, j % 2) {
                | (0, 0) => coarse.u[k],
                | (1, 0) => 0.5 * (coarse.u[k] + coarse.u[k + w]),
                | (0, 1) => 0.5 * (coarse.u[k] + coarse.u[k + 1]),
                | _ => {
                    0.25 * (coarse.u[k] + coarse.u[k + w] + coarse.u[k + 1] + coarse.u[k + w + 1])
                },
            };
        }
    };
    if parallel {
        interior
            .par_chunks_mut(n)
            .with_min_len(POISSON_ROW_BATCH)
            .enumerate()
            .for_each(correct_row);
    } else {
        interior.chunks_mut(n).enumerate().for_each(correct_row);
    }
}

fn v_cycle_2d(
    grid: &mut Grid2D,
    lower_levels: &mut [Grid2D],
    residual: &mut [f64],
    can_split: bool,
) {
    let parallel = can_split && grid.n >= PARALLEL_POISSON_SIDE;
    smooth_2d(grid);

    if let Some((coarse, rest)) = lower_levels.split_first_mut() {
        residual_2d(grid, residual, parallel);
        restrict_2d(residual, grid.n, coarse);
        coarse.u.fill(0.0);
        v_cycle_2d(coarse, rest, residual, can_split);
        correct_2d(grid, coarse, parallel);
    }

    smooth_2d(grid);
}

/// Solves `-∇²u = f` on the unit square by running `num_cycles` multigrid V-cycles.
///
/// Input and output are both flattened with x varying fastest. All boundary
/// values are fixed to zero, and source entries on those boundaries are ignored.
/// With zero cycles, or with a boundary-only 2×2 grid, the result is all zeros.
///
/// # Errors
///
/// Returns an error if `n` is not of the form `2^k + 1` for a nonnegative integer
/// `k`, or if `f` does not have exactly `n * n` entries.
pub fn solve_poisson_2d_multigrid(
    n: usize,
    f: &[f64],
    num_cycles: usize,
) -> Result<Vec<f64>, String> {
    if n < 2 || !(n - 1).is_power_of_two() {
        return Err("Grid size `n` must be of the form 2^k + 1.".into());
    }
    if n.checked_mul(n) != Some(f.len()) {
        return Err("`f` must have n^2 entries.".into());
    }
    if n == 2 || num_cycles == 0 {
        return Ok(vec![0.0; f.len()]);
    }

    let mut grid = Grid2D::new(n, 1.0 / (n - 1) as f64);
    grid.f.copy_from_slice(f);

    let levels = (n - 1).trailing_zeros() as usize;
    let mut lower_levels: Vec<Grid2D> = (1..levels)
        .map(|level| {
            let side = ((n - 1) >> level) + 1;
            Grid2D::new(side, 1.0 / (side - 1) as f64)
        })
        .collect();
    let mut residual = if lower_levels.is_empty() {
        Vec::new()
    } else {
        vec![0.0; f.len()]
    };
    let can_split = crate::execution::should_split(n, PARALLEL_POISSON_SIDE);

    for _ in 0..num_cycles {
        v_cycle_2d(&mut grid, &mut lower_levels, &mut residual, can_split);
    }

    Ok(grid.u)
}

/// Solves a 2D Poisson problem with a known analytical solution.
///
/// # Errors
///
/// This function will return an error if the underlying `solve_poisson_2d_multigrid` function
/// encounters an error, e.g., if the grid size is invalid.
pub fn simulate_2d_poisson_multigrid_scenario() -> Result<Vec<f64>, String> {
    const K: usize = 5;

    const N: usize = 2_usize.pow(K as u32) + 1;

    let h = 1.0 / (N - 1) as f64;

    let mut f = vec![0.0; N * N];

    for i in 0..N {
        for j in 0..N {
            let x = i as f64 * h;

            let y = j as f64 * h;

            f[i * N + j] = 2.0
                * std::f64::consts::PI.powi(2)
                * (std::f64::consts::PI * x).sin()
                * (std::f64::consts::PI * y).sin();
        }
    }

    solve_poisson_2d_multigrid(N, &f, 10)
}

/// Geometric multigrid for `-Δu + σu = f` on a 1D, 2D or 3D grid.
///
/// Face conditions may be Dirichlet, Neumann, Robin or periodic; V, W or full
/// multigrid cycles. Data are stored with `x` fastest.
///
/// # Errors
/// For inconsistent sizes or boundary descriptions.
pub fn solve_poisson_multigrid_nd(
    dims: &[usize],
    spacing: &[f64],
    f: &[f64],
    sigma: f64,
    bcs: &crate::sim::boundary::FaceBcs,
    options: crate::sim::multigrid::MgOptions,
) -> Result<crate::sim::multigrid::MgResult, String> {
    if dims.len() != spacing.len() || f.len() != dims.iter().product::<usize>() {
        return Err("inconsistent grid description".into());
    }
    let hi: Vec<f64> = (0..dims.len())
        .map(|a| {
            let periodic = bcs.axis_periodic(a).unwrap_or(false);
            spacing[a]
                * if periodic {
                    dims[a]
                } else {
                    dims[a] - 1
                } as f64
        })
        .collect();
    let grid = crate::sim::structured::Grid::new(dims, &vec![0.0; dims.len()], &hi, bcs)?;
    crate::sim::multigrid::MultigridSolver::new(&grid, bcs, sigma, options)?.solve(f, 0.0)
}

/// Solves the 3D Poisson equation `-Δu = f` on the unit cube.
///
/// Zero Dirichlet boundary, `n` nodes per axis including the boundary, `x` fastest, `n = 2^k + 1`)
/// with `num_cycles` multigrid V-cycles. 3D counterpart of
/// [`solve_poisson_2d_multigrid`]; the general 1D/2D/3D entry point with any
/// boundary conditions is [`solve_poisson_multigrid_nd`].
///
/// # Errors
///
/// This function will return an error if `n` is not of the form `2^k + 1`, if
/// `f.len() != n^3` or if the multigrid hierarchy cannot be built.
pub fn solve_poisson_3d_multigrid(
    n: usize,
    f: &[f64],
    num_cycles: usize,
) -> Result<Vec<f64>, String> {
    if n < 3 || !(n - 1).is_power_of_two() {
        return Err("Grid size `n` must be of the form 2^k + 1.".to_string());
    }

    if f.len() != n * n * n {
        return Err("`f` must have n^3 entries.".to_string());
    }

    let h = 1.0 / (n - 1) as f64;

    let options = crate::sim::multigrid::MgOptions {
        max_cycles: num_cycles,
        tol: 0.0,
        ..crate::sim::multigrid::MgOptions::default()
    };

    solve_poisson_multigrid_nd(
        &[n, n, n],
        &[h, h, h],
        f,
        0.0,
        &crate::sim::boundary::FaceBcs::dirichlet_zero(),
        options,
    )
    .map(|r| r.u)
}
