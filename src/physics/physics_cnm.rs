use num_complex::Complex;
use rayon::prelude::*;
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

    for i in 1 .. n {

        let m = 1.0
            / a[i - 1].mul_add(
                -c_prime[i - 1],
                b[i],
            );

        c_prime[i] = if i < n - 1 {

            c[i] * m
        } else {

            0.0
        };

        d[i] = a[i - 1]
            .mul_add(-d[i - 1], d[i])
            * m;
    }

    x[n - 1] = d[n - 1];

    for i in (0 .. n - 1).rev() {

        x[i] = c_prime[i]
            .mul_add(-x[i + 1], d[i]);
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

    let mut c_prime =
        vec![Complex::new(0.0, 0.0); n];

    let mut x =
        vec![Complex::new(0.0, 0.0); n];

    c_prime[0] = c[0] / b[0];

    d[0] /= b[0];

    for i in 1 .. n {

        let m = Complex::new(1.0, 0.0)
            / (b[i]
                - a[i - 1]
                    * c_prime[i - 1]);

        c_prime[i] = if i < n - 1 {

            c[i] * m
        } else {

            Complex::new(0.0, 0.0)
        };

        d[i] = (d[i]
            - a[i - 1] * d[i - 1])
            * m;
    }

    x[n - 1] = d[n - 1];

    for i in (0 .. n - 1).rev() {

        x[i] = d[i]
            - c_prime[i] * x[i + 1];
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

    let r = Complex::new(
        0.0,
        dt / (4.0 * dx * dx),
    );

    let mut a = vec![-r; n - 1];

    let mut b =
        vec![Complex::new(0.0, 0.0); n];

    let mut c = vec![-r; n - 1];

    // Dirichlet boundary conditions (wave function vanishes at boundaries)
    b[0] = Complex::new(1.0, 0.0);

    c[0] = Complex::new(0.0, 0.0);

    b[n - 1] = Complex::new(1.0, 0.0);

    a[n - 2] = Complex::new(0.0, 0.0);

    let mut d =
        vec![Complex::new(0.0, 0.0); n];

    for _ in 0 .. steps {

        for i in 1 .. n - 1 {

            b[i] = Complex::new(
                1.0,
                0.5 * dt * v[i],
            ) + Complex::new(
                2.0, 0.0,
            ) * r;

            d[i] = r * psi[i - 1]
                + (Complex::new(
                    1.0,
                    -0.5 * dt * v[i],
                ) - Complex::new(
                    2.0, 0.0,
                ) * r)
                    * psi[i]
                + r * psi[i + 1];
        }

        d[0] = Complex::new(0.0, 0.0);

        d[n - 1] =
            Complex::new(0.0, 0.0);

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

    let mut u =
        initial_condition.to_vec();

    let alpha =
        d_coeff * dt / (2.0 * dx * dx);

    let mut a = vec![-alpha; n - 1];

    let mut b =
        vec![
            2.0f64.mul_add(alpha, 1.0);
            n
        ];

    let mut c = vec![-alpha; n - 1];

    b[0] = 1.0;

    c[0] = 0.0;

    b[n - 1] = 1.0;

    a[n - 2] = 0.0;

    let mut d = vec![0.0; n];

    for _ in 0 .. steps {

        for i in 1 .. n - 1 {

            d[i] = alpha * u[i - 1]
                + 2.0f64.mul_add(
                    -alpha,
                    1.0,
                ) * u[i]
                + alpha * u[i + 1];
        }

        d[0] = 0.0;

        d[n - 1] = 0.0;

        u = solve_tridiagonal_system(
            &a,
            &b,
            &c,
            &mut d,
        );
    }

    u
}

/// Scenario for the 1D Crank-Nicolson solver.

#[must_use]

pub fn simulate_1d_heat_conduction_cn_scenario()
-> Vec<f64> {

    const N: usize = 100;

    const L: f64 = 1.0;

    let dx = L / (N - 1) as f64;

    let dt = 0.001;

    let d_coeff = 0.01;

    let mut u0 = vec![0.0; N];

    for (i, vars) in u0
        .iter_mut()
        .enumerate()
        .take(N)
    {

        *vars = (std::f64::consts::PI
            * i as f64
            * dx)
            .sin();
    }

    solve_heat_equation_1d_cn(
        &u0,
        dx,
        dt,
        d_coeff,
        50,
    )
}

#[derive(
    Clone, Debug, Serialize, Deserialize,
)]
/// Configuration for the 2D heat equation solver.

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
/// ADI splits the problem into two half-steps, each solving one dimension implicitly.

#[must_use]

pub fn solve_heat_equation_2d_cn_adi(
    initial_condition: &[f64],
    config: &HeatEquationSolverConfig,
) -> Vec<f64> {

    let mut u =
        initial_condition.to_vec();

    let alpha_x = config.d_coeff
        * config.dt
        / (2.0 * config.dx * config.dx);

    let alpha_y = config.d_coeff
        * config.dt
        / (2.0 * config.dy * config.dy);

    let mut ax =
        vec![-alpha_x; config.nx - 1];

    let mut bx =
        vec![
            2.0f64
                .mul_add(alpha_x, 1.0);
            config.nx
        ];

    let mut cx =
        vec![-alpha_x; config.nx - 1];

    bx[0] = 1.0;

    cx[0] = 0.0;

    bx[config.nx - 1] = 1.0;

    ax[config.nx - 2] = 0.0;

    let mut ay =
        vec![-alpha_y; config.ny - 1];

    let mut by =
        vec![
            2.0f64
                .mul_add(alpha_y, 1.0);
            config.ny
        ];

    let mut cy =
        vec![-alpha_y; config.ny - 1];

    by[0] = 1.0;

    cy[0] = 0.0;

    by[config.ny - 1] = 1.0;

    ay[config.ny - 2] = 0.0;

    let mut u_half =
        vec![
            0.0;
            config.nx * config.ny
        ];

    for _ in 0 .. config.steps {

        // Step 1: Solve implicitly in x, explicitly in y
        u_half
            .par_chunks_mut(config.nx)
            .enumerate()
            .for_each(|(j, row_half)| {

                if j == 0 || j == config.ny - 1 {

                    for i in 0 .. config.nx {

                        row_half[i] = u[j * config.nx + i];
                    }

                    return;
                }

                let mut d = vec![0.0; config.nx];

                for i in 1 .. config.nx - 1 {

                    let u_ij = u[j * config.nx + i];

                    let u_ijm1 = u[(j - 1) * config.nx + i];

                    let u_ijp1 = u[(j + 1) * config.nx + i];

                    d[i] = alpha_y * u_ijm1 + 2.0f64.mul_add(-alpha_y, 1.0) * u_ij + alpha_y * u_ijp1;
                }

                let row_sol = solve_tridiagonal_system(
                    &ax,
                    &bx,
                    &cx,
                    &mut d,
                );

                row_half[..config.nx].copy_from_slice(&row_sol[..config.nx]);
            });

        // Step 2: Solve implicitly in y, explicitly in x
        // Transpose u_half into a temporary buffer for efficient memory access
        let mut u_half_t = vec![
                0.0;
                config.nx * config.ny
            ];

        for j in 0 .. config.ny {

            for i in 0 .. config.nx {

                u_half_t[i * config
                    .ny
                    + j] = u_half
                    [j * config.nx + i];
            }
        }

        let mut u_next_t = vec![
                0.0;
                config.nx * config.ny
            ];

        u_next_t
            .par_chunks_mut(config.ny)
            .enumerate()
            .for_each(|(i, row_next_t)| {

                if i == 0 || i == config.nx - 1 {

                    for j in 0 .. config.ny {

                        row_next_t[j] = u_half_t[i * config.ny + j];
                    }

                    return;
                }

                let mut d_transposed = vec![0.0; config.ny];

                for j in 1 .. config.ny - 1 {

                    let u_im1j = u_half_t[(i - 1) * config.ny + j];

                    let u_ij = u_half_t[i * config.ny + j];

                    let u_ip1j = u_half_t[(i + 1) * config.ny + j];

                    d_transposed[j] =
                        alpha_x * u_im1j + 2.0f64.mul_add(-alpha_x, 1.0) * u_ij + alpha_x * u_ip1j;
                }

                let col_sol = solve_tridiagonal_system(
                    &ay,
                    &by,
                    &cy,
                    &mut d_transposed,
                );

                row_next_t[..config.ny].copy_from_slice(&col_sol[..config.ny]);
            });

        // Transpose back
        for i in 0 .. config.nx {

            for j in 0 .. config.ny {

                u[j * config.nx + i] =
                    u_next_t[i
                        * config.ny
                        + j];
            }
        }
    }

    u
}

/// Scenario for the 2D Crank-Nicolson ADI solver.

#[must_use]

pub fn simulate_2d_heat_conduction_cn_adi_scenario()
-> Vec<f64> {

    const NX: usize = 50;

    const NY: usize = 50;

    let dx = 1.0 / (NX - 1) as f64;

    let dy = 1.0 / (NY - 1) as f64;

    let dt = 0.01;

    let d_coeff = 0.05;

    let config =
        HeatEquationSolverConfig {
            nx: NX,
            ny: NY,
            dx,
            dy,
            dt,
            d_coeff,
            steps: 50,
        };

    let mut u0 = vec![0.0; NX * NY];

    for j in 0 .. NY {

        for i in 0 .. NX {

            let x = i as f64 * dx;

            let y = j as f64 * dy;

            if (y - 0.5).mul_add(
                y - 0.5,
                (x - 0.5).powi(2),
            ) < 0.05
            {

                u0[j * NX + i] = 100.0;
            }
        }
    }

    solve_heat_equation_2d_cn_adi(
        &u0,
        &config,
    )
}
