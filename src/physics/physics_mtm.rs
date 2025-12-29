use rayon::prelude::*;
use serde::Deserialize;
use serde::Serialize;

#[derive(
    Clone, Debug, Serialize, Deserialize,
)]
/// A 1D grid for the multigrid solver.

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

    pub(crate) const fn size(
        &self
    ) -> usize {

        self.u.len()
    }
}

/// Computes the residual r = f - Au for the 1D Poisson problem.

pub(crate) fn calculate_residual(
    grid: &Grid
) -> Vec<f64> {

    let n = grid.size();

    let h_sq_inv =
        1.0 / (grid.h * grid.h);

    let mut residual = vec![0.0; n];

    for (i, vars) in residual
        .iter_mut()
        .enumerate()
        .take(n - 1)
        .skip(1)
    {

        let a_u = (2.0f64.mul_add(
            grid.u[i],
            -grid.u[i - 1],
        ) - grid.u[i + 1])
            * h_sq_inv;

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

    for _ in 0 .. num_sweeps {

        let u_old = grid.u.clone();

        for i in 1 .. n - 1 {

            let prev = u_old[i - 1];

            let next = u_old[i + 1];

            let f_i = grid.f[i];

            let new_u_i = 0.5
                * h_sq.mul_add(
                    f_i,
                    prev + next,
                );

            grid.u[i] = (1.0_f64
                - omega)
                .mul_add(
                    u_old[i],
                    omega * new_u_i,
                );
        }
    }
}

/// Restriction: Transfers a fine-grid residual to a coarse grid using full weighting.

pub(crate) fn restrict(
    fine_residual: &[f64]
) -> Vec<f64> {

    let fine_n = fine_residual.len();

    let coarse_n = (fine_n / 2) + 1;

    let mut coarse_f =
        vec![0.0; coarse_n];

    for (i, vars) in coarse_f
        .iter_mut()
        .enumerate()
        .take(coarse_n - 1)
        .skip(1)
    {

        let j = 2 * i;

        *vars = 0.25f64.mul_add(
            fine_residual[j + 1],
            0.25f64.mul_add(
                fine_residual[j - 1],
                0.5 * fine_residual[j],
            ),
        );
    }

    coarse_f
}

/// Prolongation: Interpolates a coarse-grid correction to a fine grid.

pub(crate) fn prolongate(
    coarse_correction: &[f64]
) -> Vec<f64> {

    let coarse_n =
        coarse_correction.len();

    let fine_n = 2 * (coarse_n - 1) + 1;

    let mut fine_correction =
        vec![0.0; fine_n];

    for i in 0 .. coarse_n {

        fine_correction[2 * i] =
            coarse_correction[i];
    }

    for i in 0 .. coarse_n - 1 {

        fine_correction[2 * i + 1] = 0.5
            * (coarse_correction[i]
                + coarse_correction
                    [i + 1]);
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

        let residual =
            calculate_residual(grid);

        let coarse_f =
            restrict(&residual);

        let coarse_n = coarse_f.len();

        let mut coarse_grid = Grid::new(
            coarse_n,
            grid.h * 2.0,
        );

        coarse_grid.f = coarse_f;

        v_cycle(
            &mut coarse_grid,
            level + 1,
            max_levels,
        );

        let correction =
            prolongate(&coarse_grid.u);

        for (i, _vars) in correction
            .iter()
            .enumerate()
            .take(grid.size())
        {

            grid.u[i] += correction[i];
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

    let num_levels = (n as f64 + 1.0)
        .log2()
        as usize;

    if (2_usize.pow(num_levels as u32)
        - 1)
        != n
    {

        return Err("Grid size `n` \
                    must be of the \
                    form 2^k - 1."
            .to_string());
    }

    let mut finest_grid = Grid::new(
        n + 2,
        1.0 / (n + 1) as f64,
    );

    finest_grid.f[1 ..= n]
        .copy_from_slice(f);

    for _ in 0 .. num_cycles {

        v_cycle(
            &mut finest_grid,
            0,
            num_levels,
        );
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

pub fn simulate_1d_poisson_multigrid_scenario(
) -> Result<Vec<f64>, String> {

    const K: usize = 7;

    const N_INTERIOR: usize =
        2_usize.pow(K as u32) - 1;

    let f = vec![2.0; N_INTERIOR];

    let num_v_cycles = 10;

    solve_poisson_1d_multigrid(
        N_INTERIOR,
        &f,
        num_v_cycles,
    )
}

#[derive(
    Clone, Debug, Serialize, Deserialize,
)]
/// A 2D grid for the multigrid solver.

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
    pub(crate) fn new(
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

    pub(crate) const fn idx(
        &self,
        i: usize,
        j: usize,
    ) -> usize {

        i * self.n + j
    }
}

/// 2D Smoother: Red-Black Gauss-Seidel.

pub(crate) fn smooth_2d(
    grid: &mut Grid2D,
    num_sweeps: usize,
) {

    let n = grid.n;

    let h_sq = grid.h * grid.h;

    for _ in 0 .. num_sweeps {

        // Red points
        let u_ptr =
            grid.u.as_ptr() as usize;

        (1 .. n - 1)
            .into_par_iter()
            .for_each(|i| {
                for j in 1 .. n - 1 {

                    if (i + j) % 2 == 0 {

                        unsafe {

                            let u = u_ptr as *mut f64;

                            let u_neighbors = *u.add((i - 1) * n + j)
                                + *u.add((i + 1) * n + j)
                                + *u.add(i * n + (j - 1))
                                + *u.add(i * n + (j + 1));

                            *u.add(i * n + j) = 0.25 * h_sq.mul_add(grid.f[i * n + j], u_neighbors);
                        }
                    }
                }
            });

        // Black points
        (1 .. n - 1)
            .into_par_iter()
            .for_each(|i| {
                for j in 1 .. n - 1 {

                    if (i + j) % 2 != 0 {

                        unsafe {

                            let u = u_ptr as *mut f64;

                            let u_neighbors = *u.add((i - 1) * n + j)
                                + *u.add((i + 1) * n + j)
                                + *u.add(i * n + (j - 1))
                                + *u.add(i * n + (j + 1));

                            *u.add(i * n + j) = 0.25 * h_sq.mul_add(grid.f[i * n + j], u_neighbors);
                        }
                    }
                }
            });
    }
}

/// 2D Residual Calculation.

pub(crate) fn calculate_residual_2d(
    grid: &Grid2D
) -> Grid2D {

    let n = grid.n;

    let h_sq_inv =
        1.0 / (grid.h * grid.h);

    let mut residual_grid =
        Grid2D::new(n, grid.h);

    residual_grid.f[n .. n * (n - 1)] // Interior rows
        .par_chunks_mut(n)
        .enumerate()
        .for_each(|(i_offset, row)| {

            let i = i_offset + 1;

            for j in 1 .. n - 1 {

                let a_u = (4.0f64
                    .mul_add(
                        grid.u
                            [i * n + j],
                        -grid.u[(i
                            - 1)
                            * n
                            + j],
                    )
                    - grid.u[(i + 1)
                        * n
                        + j]
                    - grid.u[i * n
                        + (j - 1)]
                    - grid.u[i * n
                        + (j + 1)])
                    * h_sq_inv;

                row[j] = grid.f
                    [i * n + j]
                    - a_u;
            }
        });

    residual_grid
}

/// 2D Restriction: Full-weighting.

pub(crate) fn restrict_2d(
    fine_grid: &Grid2D
) -> Grid2D {

    let fine_n = fine_grid.n;

    let coarse_n = (fine_n - 1) / 2 + 1;

    let mut coarse_grid = Grid2D::new(
        coarse_n,
        fine_grid.h * 2.0,
    );

    for i in 1 .. coarse_n - 1 {

        for j in 1 .. coarse_n - 1 {

            let fi = 2 * i;

            let fj = 2 * j;

            let val = (fine_grid.f
                [fine_grid
                    .idx(fi, fj)]
            .mul_add(
                4.0,
                (fine_grid.f
                    [fine_grid.idx(
                        fi - 1,
                        fj,
                    )]
                    + fine_grid.f
                        [fine_grid
                            .idx(
                                fi + 1,
                                fj,
                            )]
                    + fine_grid.f
                        [fine_grid
                            .idx(
                                fi,
                                fj - 1,
                            )]
                    + fine_grid.f
                        [fine_grid
                            .idx(
                                fi,
                                fj + 1,
                            )])
                    * 2.0,
            ) + (fine_grid
                .f[fine_grid
                .idx(fi - 1, fj - 1)]
                + fine_grid.f
                    [fine_grid.idx(
                        fi + 1,
                        fj - 1,
                    )]
                + fine_grid.f
                    [fine_grid.idx(
                        fi - 1,
                        fj + 1,
                    )]
                + fine_grid.f
                    [fine_grid.idx(
                        fi + 1,
                        fj + 1,
                    )]))
                / 16.0;

            let helper =
                coarse_grid.idx(i, j);

            coarse_grid.f[helper] = val;
        }
    }

    coarse_grid
}

/// 2D Prolongation: Bilinear interpolation.

pub(crate) fn prolongate_2d(
    coarse_grid: &Grid2D
) -> Grid2D {

    let coarse_n = coarse_grid.n;

    let fine_n = 2 * (coarse_n - 1) + 1;

    let mut fine_grid = Grid2D::new(
        fine_n,
        coarse_grid.h / 2.0,
    );

    for i in 0 .. coarse_n {

        for j in 0 .. coarse_n {

            let helper = fine_grid
                .idx(2 * i, 2 * j);

            fine_grid.u[helper] =
                coarse_grid.u
                    [coarse_grid
                        .idx(i, j)];
        }
    }

    for i in 0 .. fine_n {

        for j in 0 .. fine_n {

            if i % 2 == 1 && j % 2 == 0
            {

                let helper_a = fine_grid
                    .u[fine_grid
                    .idx(i - 1, j)]
                    + fine_grid.u
                        [fine_grid
                            .idx(
                                i + 1,
                                j,
                            )];

                let a_helper =
                    fine_grid.idx(i, j);

                fine_grid.u[a_helper] =
                    0.5 * (helper_a);
            } else if i % 2 == 0
                && j % 2 == 1
            {

                let helper_b = fine_grid
                    .u[fine_grid
                    .idx(i, j - 1)]
                    + fine_grid.u
                        [fine_grid
                            .idx(
                                i,
                                j + 1,
                            )];

                let b_helper =
                    fine_grid.idx(i, j);

                fine_grid.u[b_helper] =
                    0.5 * (helper_b);
            } else if i % 2 == 1
                && j % 2 == 1
            {

                let helper_c = fine_grid
                    .u[fine_grid
                    .idx(i - 1, j - 1)]
                    + fine_grid.u
                        [fine_grid.idx(
                            i + 1,
                            j - 1,
                        )]
                    + fine_grid.u
                        [fine_grid.idx(
                            i - 1,
                            j + 1,
                        )]
                    + fine_grid.u
                        [fine_grid
                            .idx(
                                i + 1,
                                j + 1,
                            )];

                let c_helper =
                    fine_grid.idx(i, j);

                fine_grid.u[c_helper] =
                    0.25 * (helper_c);
            }
        }
    }

    fine_grid
}

/// 2D V-Cycle.

pub(crate) fn v_cycle_2d(
    grid: &mut Grid2D,
    level: usize,
    max_levels: usize,
) {

    smooth_2d(grid, 2);

    if level < max_levels - 1 {

        let residual_grid =
            calculate_residual_2d(grid);

        let mut coarse_grid =
            restrict_2d(&residual_grid);

        v_cycle_2d(
            &mut coarse_grid,
            level + 1,
            max_levels,
        );

        let correction_grid =
            prolongate_2d(&coarse_grid);

        for i in 0 .. grid.u.len() {

            grid.u[i] +=
                correction_grid.u[i];
        }
    }

    smooth_2d(grid, 2);
}

/// Solves the 2D Poisson equation `-∇²u = f` using the multigrid method.
///
/// # Errors
///
/// This function will return an error if the grid size `n` is not of the form `2^k + 1`.

pub fn solve_poisson_2d_multigrid(
    n: usize,
    f: &[f64],
    num_cycles: usize,
) -> Result<Vec<f64>, String> {

    let num_levels = (n as f64 - 1.0)
        .log2()
        as usize;

    if (2_usize.pow(num_levels as u32)
        + 1)
        != n
    {

        return Err("Grid size `n` \
                    must be of the \
                    form 2^k + 1."
            .to_string());
    }

    let mut finest_grid = Grid2D::new(
        n,
        1.0 / (n - 1) as f64,
    );

    finest_grid
        .f
        .copy_from_slice(f);

    for _ in 0 .. num_cycles {

        v_cycle_2d(
            &mut finest_grid,
            0,
            num_levels,
        );
    }

    Ok(finest_grid.u)
}

/// Solves a 2D Poisson problem with a known analytical solution.
///
/// # Errors
///
/// This function will return an error if the underlying `solve_poisson_2d_multigrid` function
/// encounters an error, e.g., if the grid size is invalid.

pub fn simulate_2d_poisson_multigrid_scenario(
) -> Result<Vec<f64>, String> {

    const K: usize = 5;

    const N: usize =
        2_usize.pow(K as u32) + 1;

    let h = 1.0 / (N - 1) as f64;

    let mut f = vec![0.0; N * N];

    for i in 0 .. N {

        for j in 0 .. N {

            let x = i as f64 * h;

            let y = j as f64 * h;

            f[i * N + j] = 2.0
                * std::f64::consts::PI
                    .powi(2)
                * (std::f64::consts::PI
                    * x)
                    .sin()
                * (std::f64::consts::PI
                    * y)
                    .sin();
        }
    }

    solve_poisson_2d_multigrid(
        N, &f, 10,
    )
}
