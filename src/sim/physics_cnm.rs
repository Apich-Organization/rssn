//! Crank–Nicolson time stepping.
//!
//! # Crank–Nicolson solvers
//!
//! | Equation | 1D | 2D | 3D | Boundary conditions |
//! |---|---|---|---|---|
//! | Schrödinger | yes | no | no | vanishing wavefunction |
//! | heat | yes | yes (ADI) | yes (ADI) | Dirichlet |
//!
//! The θ-method on 1D/2D/3D grids with Robin, Neumann, periodic faces,
//! time-dependent data and adaptive step control is in
//! [`crate::sim::structured`].

use num_complex::Complex;
use serde::Deserialize;
use serde::Serialize;

/// Solves a system of linear equations Ax = d where A is a tridiagonal matrix.
/// `a`: sub-diagonal (n-1 elements), `b`: main diagonal (n elements), `c`: super-diagonal (n-1 elements).
pub(crate) fn solve_tridiagonal_system(
    a: &[f64],
    b: &[f64],
    c: &[f64],
    d: &mut [f64],
) -> Vec<f64> {
    let n = b.len();

    let mut c_prime = vec![0.0; n];

    let mut x = vec![0.0; n];

    c_prime[0] = c[0] / b[0];

    d[0] /= b[0];

    for i in 1..n {
        let m = 1.0 / a[i - 1].mul_add(-c_prime[i - 1], b[i]);

        c_prime[i] = if i < n - 1 { c[i] * m } else { 0.0 };

        d[i] = a[i - 1].mul_add(-d[i - 1], d[i]) * m;
    }

    x[n - 1] = d[n - 1];

    for i in (0..n - 1).rev() {
        x[i] = c_prime[i].mul_add(-x[i + 1], d[i]);
    }

    x
}

/// Complex version of tridiagonal solver for Schrödinger equation.
pub(crate) fn solve_tridiagonal_system_complex(
    a: &[Complex<f64>],
    b: &[Complex<f64>],
    c: &[Complex<f64>],
    d: &mut [Complex<f64>],
) -> Vec<Complex<f64>> {
    let n = b.len();

    let mut c_prime = vec![Complex::new(0.0, 0.0); n];

    let mut x = vec![Complex::new(0.0, 0.0); n];

    c_prime[0] = c[0] / b[0];

    d[0] /= b[0];

    for i in 1..n {
        let m = Complex::new(1.0, 0.0) / (b[i] - a[i - 1] * c_prime[i - 1]);

        c_prime[i] = if i < n - 1 {
            c[i] * m
        } else {
            Complex::new(0.0, 0.0)
        };

        d[i] = (d[i] - a[i - 1] * d[i - 1]) * m;
    }

    x[n - 1] = d[n - 1];

    for i in (0..n - 1).rev() {
        x[i] = d[i] - c_prime[i] * x[i + 1];
    }

    x
}

/// Solves the 1D Schrödinger equation i * `psi_t` = -0.5 * `psi_xx` + V * psi
/// using the Crank-Nicolson method.
///
/// # Arguments
/// * `psi_initial` - The initial wave function (complex).
/// * `v` - The potential energy distribution (real).
/// * `dx` - The spatial step size.
/// * `dt` - The time step.
/// * `steps` - The number of time steps.
#[must_use]
pub fn solve_schrodinger_1d_cn(
    psi_initial: &[Complex<f64>],
    v: &[f64],
    dx: f64,
    dt: f64,
    steps: usize,
) -> Vec<Complex<f64>> {
    let n = psi_initial.len();

    let mut psi = psi_initial.to_vec();

    let r = Complex::new(0.0, dt / (4.0 * dx * dx));

    let mut a = vec![-r; n - 1];

    let mut b = vec![Complex::new(0.0, 0.0); n];

    let mut c = vec![-r; n - 1];

    // Dirichlet boundary conditions (wave function vanishes at boundaries)
    b[0] = Complex::new(1.0, 0.0);

    c[0] = Complex::new(0.0, 0.0);

    b[n - 1] = Complex::new(1.0, 0.0);

    a[n - 2] = Complex::new(0.0, 0.0);

    let mut d = vec![Complex::new(0.0, 0.0); n];

    for _ in 0..steps {
        for i in 1..n - 1 {
            b[i] = Complex::new(1.0, 0.5 * dt * v[i]) + Complex::new(2.0, 0.0) * r;

            d[i] = r * psi[i - 1]
                + (Complex::new(1.0, -0.5 * dt * v[i]) - Complex::new(2.0, 0.0) * r) * psi[i]
                + r * psi[i + 1];
        }

        d[0] = Complex::new(0.0, 0.0);

        d[n - 1] = Complex::new(0.0, 0.0);

        psi = solve_tridiagonal_system_complex(&a, &b, &c, &mut d);
    }

    psi
}

/// Solves the 1D heat equation `u_t` = D * `u_xx` using the Crank-Nicolson method.
///
/// # Arguments
/// * `initial_condition` - The initial temperature distribution.
/// * `dx` - The spatial step size.
/// * `dt` - The time step.
/// * `d_coeff` - The diffusion coefficient D.
/// * `steps` - The number of time steps to simulate.
///
/// # Returns
/// The final temperature distribution.
#[must_use]
pub fn solve_heat_equation_1d_cn(
    initial_condition: &[f64],
    dx: f64,
    dt: f64,
    d_coeff: f64,
    steps: usize,
) -> Vec<f64> {
    let n = initial_condition.len();
    if n >= 3 && dx > 0.0 && dt > 0.0 && steps > 0 {
        // delegate to the dimension-generic θ-method (θ = 1/2, zero Dirichlet ends)
        let bcs = crate::sim::boundary::FaceBcs::dirichlet_zero();
        if let Ok(grid) =
            crate::sim::structured::Grid::new(&[n], &[0.0], &[dx * (n - 1) as f64], &bcs)
        {
            let init = |x: [f64; 3]| initial_condition[((x[0] / dx).round() as usize).min(n - 1)];
            let prob = crate::sim::structured::HeatProblem {
                grid: &grid,
                bcs: &bcs,
                alpha: d_coeff,
                init: &init,
                source: None,
                solver: crate::sim::structured::LinearSolver::Direct,
            };
            if let Ok(u) = crate::sim::structured::heat_fixed(&prob, dt * steps as f64, steps, 0.5)
            {
                return u;
            }
        }
    }
    solve_heat_equation_1d_cn_tridiagonal(initial_condition, dx, dt, d_coeff, steps)
}

/// Direct tridiagonal Crank–Nicolson reference (zero Dirichlet ends).
#[must_use]
pub fn solve_heat_equation_1d_cn_tridiagonal(
    initial_condition: &[f64],
    dx: f64,
    dt: f64,
    d_coeff: f64,
    steps: usize,
) -> Vec<f64> {
    let n = initial_condition.len();

    let mut u = initial_condition.to_vec();

    let alpha = d_coeff * dt / (2.0 * dx * dx);

    let mut a = vec![-alpha; n - 1];

    let mut b = vec![2.0f64.mul_add(alpha, 1.0); n];

    let mut c = vec![-alpha; n - 1];

    b[0] = 1.0;

    c[0] = 0.0;

    b[n - 1] = 1.0;

    a[n - 2] = 0.0;

    let mut d = vec![0.0; n];

    for _ in 0..steps {
        for i in 1..n - 1 {
            d[i] = alpha * u[i - 1] + 2.0f64.mul_add(-alpha, 1.0) * u[i] + alpha * u[i + 1];
        }

        d[0] = 0.0;

        d[n - 1] = 0.0;

        u = solve_tridiagonal_system(&a, &b, &c, &mut d);
    }

    u
}

/// Scenario for the 1D Crank-Nicolson solver.
#[must_use]
pub fn simulate_1d_heat_conduction_cn_scenario() -> Vec<f64> {
    const N: usize = 100;

    const L: f64 = 1.0;

    let dx = L / (N - 1) as f64;

    let dt = 0.001;

    let d_coeff = 0.01;

    let mut u0 = vec![0.0; N];

    for (i, vars) in u0.iter_mut().enumerate().take(N) {
        *vars = (std::f64::consts::PI * i as f64 * dx).sin();
    }

    solve_heat_equation_1d_cn(&u0, dx, dt, d_coeff, 50)
}

/// Configuration for the 2D heat equation solver.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct HeatEquationSolverConfig {
    /// The number of grid points in the x direction.
    pub nx: usize,
    /// The number of grid points in the y direction.
    pub ny: usize,
    /// The grid spacing in the x direction.
    pub dx: f64,
    /// The grid spacing in the y direction.
    pub dy: f64,
    /// The time step.
    pub dt: f64,
    /// The diffusion coefficient.
    pub d_coeff: f64,
    /// The number of time steps to simulate.
    pub steps: usize,
}

/// Solves the 2D heat equation `u_t` = D * (`u_xx` + `u_yy`) using the ADI method.
///
/// Uses homogeneous Dirichlet boundaries from the first half-step. With
/// zero steps, returns the initial values verbatim. Grid data is x-fastest.
///
/// # Panics
///
/// Panics for a positive-step nonempty grid whose extents overflow or
/// exceed the supplied data length.
#[must_use]
pub fn solve_heat_equation_2d_cn_adi(
    initial_condition: &[f64],
    config: &HeatEquationSolverConfig,
) -> Vec<f64> {
    if config.steps == 0 || config.nx == 0 || config.ny == 0 {
        return initial_condition.to_vec();
    }
    let mut u = initial_condition.to_vec();
    let size = config
        .nx
        .checked_mul(config.ny)
        .unwrap_or_else(|| panic!("Crank–Nicolson ADI grid extents overflow"));
    let grid = &mut u[..size];
    let [nx, ny] = [config.nx, config.ny];
    for (j, row) in grid.chunks_exact_mut(nx).enumerate() {
        row[0] = 0.0;
        row[nx - 1] = 0.0;
        if j == 0 || j + 1 == ny {
            row.fill(0.0);
        }
    }
    heat_adi(grid, config, adi_coefficients(config), || Ok(()))
        .unwrap_or_else(|()| unreachable!("unchecked Crank–Nicolson ADI cannot be cancelled"));

    u
}

/// Evolves a two-dimensional heat field with fixed Dirichlet boundaries.
///
/// The x-fastest initial field supplies the four boundary edges. They remain
/// fixed during both ADI half-steps. The grid has `nx * ny` nodes, including
/// the boundaries. Spacings must be positive, and diffusivity and the time
/// step must be nonnegative. A grid without interior nodes is stationary.
///
/// # Errors
///
/// Returns an error for shape mismatches, non-finite data or results, invalid
/// physical parameters, unrepresentable coefficients, or cancellation.
pub fn solve_heat_equation_2d_cn_adi_fixed_boundary(
    initial_condition: &[f64],
    config: &HeatEquationSolverConfig,
) -> Result<Vec<f64>, String> {
    if config.nx == 0
        || config.ny == 0
        || config.nx.checked_mul(config.ny) != Some(initial_condition.len())
        || initial_condition.iter().any(|value| !value.is_finite())
        || [config.dx, config.dy]
            .iter()
            .any(|&spacing| !spacing.is_finite() || spacing <= 0.0)
        || [config.d_coeff, config.dt]
            .iter()
            .any(|&parameter| !parameter.is_finite() || parameter < 0.0)
    {
        return Err("invalid Crank–Nicolson ADI field or controls".into());
    }
    let checkpoint = || crate::cancellation::checkpoint().map_err(|error| error.to_string());
    checkpoint()?;
    if config.steps == 0 || config.nx <= 2 || config.ny <= 2 {
        return Ok(initial_condition.to_vec());
    }
    let alpha = adi_coefficients(config);
    if alpha
        .iter()
        .any(|factor| !factor.mul_add(2.0, 1.0).is_finite())
    {
        return Err("Crank–Nicolson ADI coefficient is not representable".into());
    }
    let mut u = initial_condition.to_vec();
    heat_adi(&mut u, config, alpha, checkpoint)?;
    if u.iter().any(|value| !value.is_finite()) {
        return Err("Crank–Nicolson ADI state is not representable".into());
    }
    Ok(u)
}

fn adi_coefficients(config: &HeatEquationSolverConfig) -> [f64; 2] {
    let product = config.d_coeff * config.dt;
    [config.dx, config.dy].map(|spacing| {
        let denominator = 2.0 * spacing * spacing;
        let direct = product / denominator;
        if (product.is_normal() && denominator.is_normal())
            || !config.d_coeff.is_finite()
            || !config.dt.is_finite()
            || !spacing.is_finite()
            || spacing <= 0.0
        {
            direct
        } else if config.d_coeff == 0.0 || config.dt == 0.0 {
            0.0
        } else {
            crate::kernels::numeric_common::scale_product(
                config.d_coeff.copysign(product),
                &[config.dt.abs()],
                &[2.0, spacing, spacing],
            )
            .unwrap_or(direct)
        }
    })
}

struct AdiFactors {
    lower: f64,
    // Reciprocal pivot and modified upper diagonal at each interior node.
    interior: Vec<[f64; 2]>,
}

impl AdiFactors {
    fn new(
        nodes: usize,
        alpha: f64,
    ) -> Self {
        let lower = -alpha;
        let diagonal = 2.0_f64.mul_add(alpha, 1.0);
        let mut interior = Vec::with_capacity(nodes - 2);
        let mut upper = 0.0;
        for _ in 1..nodes - 1 {
            let inverse = 1.0 / lower.mul_add(-upper, diagonal);
            upper = lower * inverse;
            interior.push([inverse, upper]);
        }
        Self { lower, interior }
    }

    fn solve(
        &self,
        line: &mut [f64],
    ) {
        let last = line.len() - 1;
        // The identity rows for the fixed endpoints seed the two substitutions.
        let mut previous = line[0];
        let mut next = line[last];
        let interior = &mut line[1..last];
        for (value, &[inverse, _]) in interior.iter_mut().zip(&self.interior) {
            previous = self.lower.mul_add(-previous, *value) * inverse;
            *value = previous;
        }
        for (value, &[_, upper]) in interior.iter_mut().zip(&self.interior).rev() {
            next = upper.mul_add(-next, *value);
            *value = next;
        }
    }
}

fn heat_adi<E>(
    u: &mut [f64],
    config: &HeatEquationSolverConfig,
    alpha: [f64; 2],
    checkpoint: impl Fn() -> Result<(), E>,
) -> Result<(), E> {
    let [nx, ny] = [config.nx, config.ny];
    let [alpha_x, alpha_y] = alpha;
    if nx <= 2 || ny <= 2 || (alpha_x == 0.0 && alpha_y == 0.0) {
        return Ok(());
    }
    let x_factors = AdiFactors::new(nx, alpha_x);
    let y_factors = AdiFactors::new(ny, alpha_y);
    let mut half = crate::kernels::reuse::f64_scratch(u.len());
    let mut transposed = crate::kernels::reuse::f64_scratch(u.len());
    half.resize(u.len(), 0.0);
    transposed.resize(u.len(), 0.0);
    for _ in 0..config.steps {
        checkpoint()?;
        adi_sweep(u, &mut half, nx, alpha_y, &x_factors);
        transpose_heat(&half, &mut transposed, nx, ny);
        adi_sweep(&transposed, &mut half, ny, alpha_x, &y_factors);
        transpose_heat(&half, u, ny, nx);
    }
    checkpoint()
}

fn adi_sweep(
    field: &[f64],
    result: &mut [f64],
    line_len: usize,
    alpha: f64,
    implicit: &AdiFactors,
) {
    let lines = result.len() / line_len;
    let rows_per_chunk = (4096 / line_len).max(1);
    let center = 2.0_f64.mul_add(-alpha, 1.0);
    crate::kernels::par::for_chunks_mut(result, rows_per_chunk * line_len, 4, |chunk, rows| {
        for (local, row) in rows.chunks_exact_mut(line_len).enumerate() {
            let j = chunk * rows_per_chunk + local;
            let start = j * line_len;
            let current = &field[start..start + line_len];
            if j == 0 || j + 1 == lines {
                row.copy_from_slice(current);
                continue;
            }
            let lower = &field[start - line_len..start];
            let upper = &field[start + line_len..start + 2 * line_len];
            for i in 1..line_len - 1 {
                row[i] = alpha * lower[i] + center * current[i] + alpha * upper[i];
            }
            row[0] = current[0];
            row[line_len - 1] = current[line_len - 1];
            implicit.solve(row);
        }
    });
}

fn transpose_heat(
    source: &[f64],
    target: &mut [f64],
    nx: usize,
    ny: usize,
) {
    const TILE: usize = 32;
    let columns_per_chunk = nx.min(TILE);
    crate::kernels::par::for_chunks_mut(target, columns_per_chunk * ny, 4, |block, columns| {
        let x0 = block * columns_per_chunk;
        for y0 in (0..ny).step_by(TILE) {
            let y_end = (y0 + TILE).min(ny);
            for (column, output) in columns.chunks_exact_mut(ny).enumerate() {
                let x = x0 + column;
                for y in y0..y_end {
                    output[y] = source[y * nx + x];
                }
            }
        }
    });
}

/// Scenario for the 2D Crank-Nicolson ADI solver.
#[must_use]
pub fn simulate_2d_heat_conduction_cn_adi_scenario() -> Vec<f64> {
    const NX: usize = 50;

    const NY: usize = 50;

    let dx = 1.0 / (NX - 1) as f64;

    let dy = 1.0 / (NY - 1) as f64;

    let dt = 0.01;

    let d_coeff = 0.05;

    let config = HeatEquationSolverConfig {
        nx: NX,
        ny: NY,
        dx,
        dy,
        dt,
        d_coeff,
        steps: 50,
    };

    let mut u0 = vec![0.0; NX * NY];

    for j in 0..NY {
        for i in 0..NX {
            let x = i as f64 * dx;

            let y = j as f64 * dy;

            if (y - 0.5).mul_add(y - 0.5, (x - 0.5).powi(2)) < 0.05 {
                u0[j * NX + i] = 100.0;
            }
        }
    }

    solve_heat_equation_2d_cn_adi(&u0, &config)
}

/// Crank–Nicolson (θ = 1/2) heat equation on a 1D, 2D or 3D grid.
///
/// Any face conditions of [`crate::sim::boundary`] are allowed (Dirichlet,
/// Neumann, Robin, periodic, time-dependent data); data are stored with `x` fastest.
///
/// # Errors
/// For inconsistent sizes or boundary descriptions.
pub fn solve_heat_equation_nd_cn(
    dims: &[usize],
    spacing: &[f64],
    initial_condition: &[f64],
    d_coeff: f64,
    dt: f64,
    steps: usize,
    bcs: &crate::sim::boundary::FaceBcs,
) -> Result<Vec<f64>, String> {
    use crate::sim::structured::Grid;
    use crate::sim::structured::HeatProblem;
    use crate::sim::structured::LinearSolver;
    use crate::sim::structured::heat_fixed;
    if dims.len() != spacing.len()
        || initial_condition.len() != dims.iter().product::<usize>()
        || steps == 0
    {
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
    let grid = Grid::new(dims, &vec![0.0; dims.len()], &hi, bcs)?;
    let init = |x: [f64; 3]| {
        let mut idx = [0usize; 3];
        for a in 0..dims.len() {
            idx[a] = ((x[a] / spacing[a]).round() as usize).min(dims[a] - 1);
        }
        initial_condition[grid.index(idx)]
    };
    let prob = HeatProblem {
        grid: &grid,
        bcs,
        alpha: d_coeff,
        init: &init,
        source: None,
        solver: LinearSolver::Auto,
    };
    heat_fixed(&prob, dt * steps as f64, steps, 0.5)
}

// ============================================================================
// Space-dimension coverage: 3D heat (ADI) and 2D/3D Schrodinger
// ============================================================================

/// Dimension-generic factorised Crank-Nicolson heat solver for 1D, 2D and 3D
/// grids (zero Dirichlet).
pub use crate::sim::dimensional::solve_heat_adi_nd;
/// Dimension-generic Cayley-split Crank-Nicolson Schrodinger solver for 1D, 2D
/// and 3D grids (zero Dirichlet).
pub use crate::sim::dimensional::solve_schrodinger_cn_nd;

/// Configuration for the 3D Crank-Nicolson ADI heat solver.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct HeatEquationSolverConfig3D {
    /// Number of grid points in x.
    pub nx: usize,
    /// Number of grid points in y.
    pub ny: usize,
    /// Number of grid points in z.
    pub nz: usize,
    /// Grid spacing in x.
    pub dx: f64,
    /// Grid spacing in y.
    pub dy: f64,
    /// Grid spacing in z.
    pub dz: f64,
    /// Time step.
    pub dt: f64,
    /// Diffusion coefficient.
    pub d_coeff: f64,
    /// Number of time steps.
    pub steps: usize,
}

/// 3D heat equation `u_t = D (u_xx + u_yy + u_zz)` by Crank-Nicolson.
///
/// Uses approximate factorisation (ADI): three tridiagonal sweeps per step,
/// unconditionally stable, second order, zero Dirichlet boundaries. `x` is the
/// fastest index. 3D counterpart of [`solve_heat_equation_2d_cn_adi`].
///
/// # Errors
/// For a wrong array length or non-positive step sizes.
pub fn solve_heat_equation_3d_cn_adi(
    initial_condition: &[f64],
    config: &HeatEquationSolverConfig3D,
) -> Result<Vec<f64>, String> {
    solve_heat_adi_nd(
        initial_condition,
        &[config.nx, config.ny, config.nz],
        &[config.dx, config.dy, config.dz],
        config.d_coeff,
        config.dt,
        config.steps,
    )
}

/// 2D Schrodinger equation `i psi_t = -1/2 Laplacian psi + V psi`.
///
/// Zero Dirichlet, by the Cayley/Crank-Nicolson split of
/// [`solve_schrodinger_cn_nd`]; `x` fastest. 2D counterpart of [`solve_schrodinger_1d_cn`]; conserves the norm.
///
/// # Errors
/// For inconsistent lengths or a non-positive step.
pub fn solve_schrodinger_2d_cn(
    psi_initial: &[Complex<f64>],
    v: &[f64],
    dims: [usize; 2],
    h: [f64; 2],
    dt: f64,
    steps: usize,
) -> Result<Vec<Complex<f64>>, String> {
    solve_schrodinger_cn_nd(psi_initial, v, &dims, &h, dt, steps)
}

/// 3D Schrodinger equation by the Cayley/Crank-Nicolson split; see
/// [`solve_schrodinger_2d_cn`].
///
/// # Errors
/// As [`solve_schrodinger_2d_cn`].
pub fn solve_schrodinger_3d_cn(
    psi_initial: &[Complex<f64>],
    v: &[f64],
    dims: [usize; 3],
    h: [f64; 3],
    dt: f64,
    steps: usize,
) -> Result<Vec<Complex<f64>>, String> {
    solve_schrodinger_cn_nd(psi_initial, v, &dims, &h, dt, steps)
}
