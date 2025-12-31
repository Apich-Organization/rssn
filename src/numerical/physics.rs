//! # Numerical Physics
//!
//! This module provides numerical methods for solving physics problems including:
//!
//! ## Classical Mechanics
//! - Particle motion simulation under force fields
//! - N-body gravitational dynamics
//! - Projectile motion with drag
//! - Harmonic oscillator systems
//!
//! ## Quantum Mechanics
//! - 1D/2D/3D Schrödinger equation solvers (finite difference)
//! - Quantum harmonic oscillator energy levels
//!
//! ## Statistical Mechanics
//! - Ising model simulation
//! - Maxwell-Boltzmann distribution
//! - Partition function calculations
//!
//! ## Thermodynamics
//! - Heat equation (Crank-Nicolson method)
//! - Ideal gas law calculations
//!
//! ## Waves
//! - 1D wave equation solver
//! - Standing wave analysis
//!
//! ## Electromagnetism
//! - Coulomb's law
//! - Biot-Savart law (magnetic field)
//! - Electric field from charge distributions
//!
//! ## Physical Constants (SI units)
//! All constants are provided in standard SI units.

use std::collections::HashMap;
use std::sync::Arc;

use rand::Rng;
use rand::thread_rng;
use serde::Deserialize;
use serde::Serialize;

use crate::numerical::elementary::eval_expr;
use crate::numerical::matrix::Matrix;
use crate::numerical::ode::solve_ode_system_rk4;
use crate::symbolic::core::Expr;

// ============================================================================
// Physical Constants (SI units)
// ============================================================================

/// Speed of light in vacuum (m/s)

pub const SPEED_OF_LIGHT: f64 =
    299_792_458.0;

/// Planck's constant (J·s)

pub const PLANCK_CONSTANT: f64 =
    6.626_070_15e-34;

/// Reduced Planck's constant ħ = h/(2π) (J·s)

pub const HBAR: f64 = 1.054_571_817e-34;

/// Elementary charge (C)

pub const ELEMENTARY_CHARGE: f64 =
    1.602_176_634e-19;

/// Electron mass (kg)

pub const ELECTRON_MASS: f64 =
    9.109_383_56e-31;

/// Proton mass (kg)

pub const PROTON_MASS: f64 =
    1.672_621_898e-27;

/// Neutron mass (kg)

pub const NEUTRON_MASS: f64 =
    1.674_927_351e-27;

/// Gravitational constant (m³/(kg·s²))

pub const GRAVITATIONAL_CONSTANT: f64 =
    6.674_30e-11;

/// Avogadro's number (mol⁻¹)

pub const AVOGADRO_NUMBER: f64 =
    6.022_140_76e23;

/// Boltzmann constant (J/K)

pub const BOLTZMANN_CONSTANT: f64 =
    1.380_649e-23;

/// Gas constant R = NA × kB (J/(mol·K))

pub const GAS_CONSTANT: f64 =
    8.314_462_618;

/// Stefan-Boltzmann constant (W/(m²·K⁴))

pub const STEFAN_BOLTZMANN: f64 =
    5.670_374_419e-8;

/// Vacuum permittivity ε₀ (F/m)

pub const VACUUM_PERMITTIVITY: f64 =
    8.854_187_817e-12;

/// Vacuum permeability μ₀ (H/m)

pub const VACUUM_PERMEABILITY: f64 =
    1.256_637_061e-6;

/// Coulomb constant k = 1/(4πε₀) (N·m²/C²)

pub const COULOMB_CONSTANT: f64 =
    8.987_551_787e9;

/// Standard Earth gravity (m/s²)

pub const STANDARD_GRAVITY: f64 =
    9.806_65;

/// Atomic mass unit (kg)

pub const ATOMIC_MASS_UNIT: f64 =
    1.660_539_067e-27;

/// Bohr radius (m)

pub const BOHR_RADIUS: f64 =
    5.291_772_109e-11;

/// Fine structure constant

pub const FINE_STRUCTURE_CONSTANT: f64 =
    7.297_352_566e-3;

// ============================================================================
// Classical Mechanics Types
// ============================================================================

/// A particle with mass, position, and velocity in 3D space.
#[derive(
    Debug,
    Clone,
    Copy,
    PartialEq,
    Serialize,
    Deserialize,
)]

pub struct Particle3D {
    /// Mass of the particle.
    pub mass: f64,
    /// X position.
    pub x: f64,
    /// Y position.
    pub y: f64,
    /// Z position.
    pub z: f64,
    /// X velocity.
    pub vx: f64,
    /// Y velocity.
    pub vy: f64,
    /// Z velocity.
    pub vz: f64,
}

impl Particle3D {
    /// Creates a new particle.
    #[must_use]

    pub const fn new(
        mass: f64,
        x: f64,
        y: f64,
        z: f64,
        vx: f64,
        vy: f64,
        vz: f64,
    ) -> Self {

        Self {
            mass,
            x,
            y,
            z,
            vx,
            vy,
            vz,
        }
    }

    /// Kinetic energy of the particle.
    #[must_use]

    pub fn kinetic_energy(
        &self
    ) -> f64 {

        0.5 * self.mass
            * self.vz.mul_add(
                self.vz,
                self.vx.mul_add(
                    self.vx,
                    self.vy * self.vy,
                ),
            )
    }

    /// Momentum magnitude.
    #[must_use]

    pub fn momentum(&self) -> f64 {

        self.mass
            * self
                .vz
                .mul_add(
                    self.vz,
                    self.vx.mul_add(
                        self.vx,
                        self.vy
                            * self.vy,
                    ),
                )
                .sqrt()
    }
}

// ============================================================================
// Classical Mechanics Functions
// ============================================================================

/// Simulates the 3D motion of a particle under a force field.
///
/// This function integrates Newton's second law (`F=ma`) for a single particle
/// in a given force field using the fourth-order Runge-Kutta method.
///
/// # Arguments
/// * `force_exprs` - A tuple of `Expr` representing the force components `(Fx, Fy, Fz)`.
/// * `mass` - The mass of the particle.
/// * `initial_pos` - The initial position `(x, y, z)`.
/// * `initial_vel` - The initial velocity `(vx, vy, vz)`.
/// * `t_range` - The time interval `(t_start, t_end)` for the simulation.
/// * `num_steps` - The number of time steps.
///
/// # Returns
/// A `Result` containing a `Vec<Vec<f64>>` representing the trajectory over time,
/// or an error string if evaluation fails.
///
/// # Errors
/// Returns an error if symbolic expression evaluation or ODE solving fails.

pub fn simulate_particle_motion(
    force_exprs: (&Expr, &Expr, &Expr),
    mass: f64,
    initial_pos: (f64, f64, f64),
    initial_vel: (f64, f64, f64),
    t_range: (f64, f64),
    num_steps: usize,
) -> Result<Vec<Vec<f64>>, String> {

    let (fx_expr, fy_expr, fz_expr) =
        force_exprs;

    let m_expr = Expr::Constant(mass);

    let ode_funcs: Vec<Expr> = vec![
        Expr::Variable(
            "y3".to_string(),
        ),
        Expr::Variable(
            "y4".to_string(),
        ),
        Expr::Variable(
            "y5".to_string(),
        ),
        Expr::Div(
            Arc::new(fx_expr.clone()),
            Arc::new(m_expr.clone()),
        ),
        Expr::Div(
            Arc::new(fy_expr.clone()),
            Arc::new(m_expr.clone()),
        ),
        Expr::Div(
            Arc::new(fz_expr.clone()),
            Arc::new(m_expr),
        ),
    ];

    let y0 = vec![
        initial_pos.0,
        initial_pos.1,
        initial_pos.2,
        initial_vel.0,
        initial_vel.1,
        initial_vel.2,
    ];

    solve_ode_system_rk4(
        &ode_funcs,
        &y0,
        t_range,
        num_steps,
    )
}

/// Simulates a 2D Ising model using the Metropolis-Hastings algorithm.
///
/// The Ising model is a mathematical model of ferromagnetism in statistical mechanics.
/// This function simulates the evolution of a 2D lattice of spins using the Metropolis-Hastings
/// algorithm to reach thermal equilibrium at a given temperature.
///
/// # Arguments
/// * `size` - The dimension of the square lattice (e.g., `size x size`).
/// * `temperature` - The temperature of the system.
/// * `steps` - The number of Monte Carlo steps to perform.
///
/// # Returns
/// A `Vec<Vec<i8>>` representing the final spin configuration of the lattice.
#[allow(clippy::needless_range_loop)]
#[must_use]

pub fn simulate_ising_model(
    size: usize,
    temperature: f64,
    steps: usize,
) -> Vec<Vec<i8>> {

    let mut rng = thread_rng();

    let mut lattice =
        vec![vec![0i8; size]; size];

    for i in 0 .. steps {

        for j in 0 .. size {

            lattice[i][j] = if rng
                .r#gen::<bool>()
            {

                1
            } else {

                -1
            };
        }
    }

    for _ in 0 .. steps {

        let i =
            rng.gen_range(0 .. size);

        let j =
            rng.gen_range(0 .. size);

        let top = lattice
            [(i + size - 1) % size][j];

        let bottom =
            lattice[(i + 1) % size][j];

        let left = lattice[i]
            [(j + size - 1) % size];

        let right =
            lattice[i][(j + 1) % size];

        let neighbor_sum =
            top + bottom + left + right;

        let delta_e = 2.0
            * f64::from(lattice[i][j])
            * f64::from(neighbor_sum);

        if delta_e < 0.0
            || rng.r#gen::<f64>()
                < (-delta_e
                    / temperature)
                    .exp()
        {

            lattice[i][j] *= -1;
        }
    }

    lattice
}

/// Solves the 1D time-independent Schrödinger equation `Hψ = Eψ` using the finite difference method.
///
/// This function discretizes the 1D Schrödinger equation on a grid, transforming it into
/// an eigenvalue problem for the Hamiltonian matrix. The eigenvalues correspond to energy
/// levels, and the eigenvectors correspond to the wave functions.
/// Assumes `hbar = 1` and `m = 1`.
///
/// # Arguments
/// * `potential_expr` - The symbolic expression for the potential `V(x)`.
/// * `var` - The independent variable (e.g., "x").
/// * `range` - The spatial range `(x_min, x_max)`.
/// * `num_points` - The number of grid points.
///
/// # Returns
/// A `Result` containing a tuple `(eigenvalues, eigenvectors_matrix)`,
/// or an error string if the number of points is too small or evaluation fails.
///
/// # Errors
/// Returns an error if `num_points` is less than 3, if symbolic evaluation fails, or if the decomposition fails to converge.

pub fn solve_1d_schrodinger(
    potential_expr: &Expr,
    var: &str,
    range: (f64, f64),
    num_points: usize,
) -> Result<
    (
        Vec<f64>,
        Matrix<f64>,
    ),
    String,
> {

    let (x_min, x_max) = range;

    if num_points < 3 {

        return Err("num_points must \
                    be >= 3"
            .to_string());
    }

    let dx = (x_max - x_min)
        / (num_points as f64 - 1.0);

    let points: Vec<f64> = (0
        .. num_points)
        .map(|i| {

            (i as f64)
                .mul_add(dx, x_min)
        })
        .collect();

    let mut potential_values =
        Vec::with_capacity(num_points);

    for &x in &points {

        let mut vars = HashMap::new();

        vars.insert(var.to_string(), x);

        potential_values.push(
            eval_expr(
                potential_expr,
                &vars,
            )?,
        );
    }

    let n = num_points;

    let mut h_data = vec![0.0; n * n];

    let factor = -0.5 / (dx * dx);

    for i in 0 .. n {

        h_data[i * n + i] = (-2.0f64)
            .mul_add(
                factor,
                potential_values[i],
            );

        if i > 0 {

            h_data[i * n + i - 1] =
                factor;
        }

        if i + 1 < n {

            h_data[i * n + i + 1] =
                factor;
        }
    }

    let hamiltonian =
        Matrix::new(n, n, h_data);

    hamiltonian
        .jacobi_eigen_decomposition(
            2000, 1e-12,
        )
}

/// Solves the 2D time-independent Schrödinger equation `Hψ = Eψ` on a rectangular grid.
///
/// This function discretizes the 2D Schrödinger equation using a 5-point Laplacian
/// finite-difference stencil, converting it into a large eigenvalue problem.
/// Dirichlet boundary conditions (`ψ=0`) are enforced by setting a large potential at boundary points.
/// Warning: For large grids, this dense matrix approach can be very slow and memory-intensive.
///
/// # Arguments
/// * `potential_expr` - The symbolic expression for the potential `V(x,y)`.
/// * `var_x`, `var_y` - The independent variables (e.g., "x", "y").
/// * `ranges` - The spatial ranges `(x_min, x_max, y_min, y_max)`.
/// * `grid` - The number of grid points `(nx, ny)` in each direction (including boundaries).
///
/// # Returns
/// A `Result` containing a tuple `(eigenvalues, eigenvectors_matrix)`,
/// or an error string if grid dimensions are too small or evaluation fails.
///
/// # Errors
/// Returns an error if grid dimensions are less than 3, if symbolic evaluation fails, or if the decomposition fails to converge.

pub fn solve_2d_schrodinger(
    potential_expr: &Expr,
    var_x: &str,
    var_y: &str,
    ranges: (f64, f64, f64, f64),
    grid: (usize, usize),
) -> Result<
    (
        Vec<f64>,
        Matrix<f64>,
    ),
    String,
> {

    let (x_min, x_max, y_min, y_max) =
        ranges;

    let (nx, ny) = grid;

    if nx < 3 || ny < 3 {

        return Err("grid dimensions \
                    must be at \
                    least 3 in each \
                    direction"
            .to_string());
    }

    let dx = (x_max - x_min)
        / (nx as f64 - 1.0);

    let dy = (y_max - y_min)
        / (ny as f64 - 1.0);

    let n = nx * ny;

    let mut potential =
        vec![0.0_f64; n];

    for ix in 0 .. nx {

        for iy in 0 .. ny {

            let x = (ix as f64)
                .mul_add(dx, x_min);

            let y = (iy as f64)
                .mul_add(dy, y_min);

            let mut vars =
                HashMap::new();

            vars.insert(
                var_x.to_string(),
                x,
            );

            vars.insert(
                var_y.to_string(),
                y,
            );

            potential[ix * ny + iy] =
                eval_expr(
                    potential_expr,
                    &vars,
                )?;
        }
    }

    let mut h_data =
        vec![0.0_f64; n * n];

    let fx = -0.5 / (dx * dx);

    let fy = -0.5 / (dy * dy);

    let large = 1e12;

    for ix in 0 .. nx {

        for iy in 0 .. ny {

            let idx = ix * ny + iy;

            if ix == 0
                || ix == nx - 1
                || iy == 0
                || iy == ny - 1
            {

                h_data[idx * n + idx] =
                    large
                        + potential[idx];

                continue;
            }

            h_data[idx * n + idx] =
                (-2.0f64).mul_add(
                    fx + fy,
                    potential[idx],
                );

            let idx_left =
                (ix - 1) * ny + iy;

            let idx_right =
                (ix + 1) * ny + iy;

            h_data
                [idx * n + idx_left] =
                fx;

            h_data
                [idx * n + idx_right] =
                fx;

            let idx_down =
                ix * ny + (iy - 1);

            let idx_up =
                ix * ny + (iy + 1);

            h_data
                [idx * n + idx_down] =
                fy;

            h_data[idx * n + idx_up] =
                fy;
        }
    }

    let hamiltonian =
        Matrix::new(n, n, h_data);

    hamiltonian
        .jacobi_eigen_decomposition(
            5000, 1e-10,
        )
}

/// Solves the 3D time-independent Schrödinger equation `Hψ = Eψ` on a rectangular grid.
///
/// This function discretizes the 3D Schrödinger equation using a 7-point Laplacian
/// finite-difference stencil, converting it into a large eigenvalue problem.
/// Dirichlet boundary conditions (`ψ=0`) are enforced by setting a large potential at boundary points.
/// Warning: For large grids, this dense matrix approach can be extremely slow and memory-intensive.
/// Consider smaller grids or sparse/iterative methods for higher performance.
///
/// # Arguments
/// * `potential_expr` - The symbolic expression for the potential `V(x,y,z)`.
/// * `var_x`, `var_y`, `var_z` - The independent variables (e.g., "x", "y", "z").
/// * `ranges` - The spatial ranges `(x_min, x_max, y_min, y_max, z_min, z_max)`.
/// * `grid` - The number of grid points `(nx, ny, nz)` in each direction (including boundaries).
///
/// # Returns
/// A `Result` containing a tuple `(eigenvalues, eigenvectors_matrix)`,
/// or an error string if grid dimensions are too small, too large, or evaluation fails.
///
/// # Errors
/// Returns an error if grid dimensions are less than 3, if the total number of points exceeds 25,000, if symbolic evaluation fails, or if the decomposition fails to converge.

pub fn solve_3d_schrodinger(
    potential_expr: &Expr,
    var_x: &str,
    var_y: &str,
    var_z: &str,
    ranges: (
        f64,
        f64,
        f64,
        f64,
        f64,
        f64,
    ),
    grid: (usize, usize, usize),
) -> Result<
    (
        Vec<f64>,
        Matrix<f64>,
    ),
    String,
> {

    let (
        x_min,
        x_max,
        y_min,
        y_max,
        z_min,
        z_max,
    ) = ranges;

    let (nx, ny, nz) = grid;

    if nx < 3 || ny < 3 || nz < 3 {

        return Err("grid dimensions \
                    must be at \
                    least 3 in each \
                    direction"
            .to_string());
    }

    let dx = (x_max - x_min)
        / (nx as f64 - 1.0);

    let dy = (y_max - y_min)
        / (ny as f64 - 1.0);

    let dz = (z_max - z_min)
        / (nz as f64 - 1.0);

    let n = nx * ny * nz;

    if n > 25000 {

        return Err(format!(
            "Grid too large (nx*ny*nz \
             = {n}). Dense 3D solver \
             will be extremely slow \
             and memory-heavy. \
             Consider smaller grid or \
             sparse/iterative methods."
        ));
    }

    let mut potential =
        vec![0.0_f64; n];

    for ix in 0 .. nx {

        for iy in 0 .. ny {

            for iz in 0 .. nz {

                let x = (ix as f64)
                    .mul_add(dx, x_min);

                let y = (iy as f64)
                    .mul_add(dy, y_min);

                let z = (iz as f64)
                    .mul_add(dz, z_min);

                let mut vars =
                    HashMap::new();

                vars.insert(
                    var_x.to_string(),
                    x,
                );

                vars.insert(
                    var_y.to_string(),
                    y,
                );

                vars.insert(
                    var_z.to_string(),
                    z,
                );

                potential[(ix * ny
                    + iy)
                    * nz
                    + iz] = eval_expr(
                    potential_expr,
                    &vars,
                )?;
            }
        }
    }

    let mut h_data =
        vec![0.0_f64; n * n];

    let fx = -0.5 / (dx * dx);

    let fy = -0.5 / (dy * dy);

    let fz = -0.5 / (dz * dz);

    let large = 1e12;

    for ix in 0 .. nx {

        for iy in 0 .. ny {

            for iz in 0 .. nz {

                let idx =
                    (ix * ny + iy) * nz
                        + iz;

                if ix == 0
                    || ix == nx - 1
                    || iy == 0
                    || iy == ny - 1
                    || iz == 0
                    || iz == nz - 1
                {

                    h_data[idx * n
                        + idx] = large
                        + potential
                            [idx];

                    continue;
                }

                h_data[idx * n + idx] =
                    (-2.0f64).mul_add(
                        fx + fy + fz,
                        potential[idx],
                    );

                let idx_xm = ((ix - 1)
                    * ny
                    + iy)
                    * nz
                    + iz;

                let idx_xp = ((ix + 1)
                    * ny
                    + iy)
                    * nz
                    + iz;

                let idx_ym = (ix * ny
                    + (iy - 1))
                    * nz
                    + iz;

                let idx_yp = (ix * ny
                    + (iy + 1))
                    * nz
                    + iz;

                let idx_zm =
                    (ix * ny + iy) * nz
                        + (iz - 1);

                let idx_zp =
                    (ix * ny + iy) * nz
                        + (iz + 1);

                h_data[idx * n
                    + idx_xm] = fx;

                h_data[idx * n
                    + idx_xp] = fx;

                h_data[idx * n
                    + idx_ym] = fy;

                h_data[idx * n
                    + idx_yp] = fy;

                h_data[idx * n
                    + idx_zm] = fz;

                h_data[idx * n
                    + idx_zp] = fz;
            }
        }
    }

    let hamiltonian =
        Matrix::new(n, n, h_data);

    hamiltonian
        .jacobi_eigen_decomposition(
            5000, 1e-10,
        )
}

/// Solves the 1D Heat Equation `u_t = alpha * u_xx` using the Crank-Nicolson method.
///
/// The Crank-Nicolson method is a finite difference method for numerically solving
/// the heat equation and similar partial differential equations. It is implicit,
/// unconditionally stable, and second-order accurate in both space and time.
/// Assumes Dirichlet boundary conditions `u(a,t)=u(b,t)=0`.
///
/// # Arguments
/// * `init_func` - A closure representing the initial condition `u(x,0)`.
/// * `alpha` - The thermal diffusivity constant.
/// * `range` - The spatial range `(a, b)`.
/// * `nx` - The number of spatial grid points.
/// * `t_range` - The time range `(t_start, t_end)`.
/// * `nt` - The number of time steps.
///
/// # Returns
/// A `Result` containing a `Vec<Vec<f64>>` where each inner `Vec` is the solution `u`
/// at a given time step, or an error string if input dimensions are invalid.
///
/// # Errors
/// Returns an error if `nx` < 3 or `nt` < 1.
#[allow(clippy::suspicious_operation_groupings)]

pub fn solve_heat_equation_1d_crank_nicolson(
    init_func: &dyn Fn(f64) -> f64,
    alpha: f64,
    range: (f64, f64),
    nx: usize,
    t_range: (f64, f64),
    nt: usize,
) -> Result<Vec<Vec<f64>>, String> {

    let (a, b) = range;

    if nx < 3 || nt < 1 {

        return Err("nx must be >=3 \
                    and nt >= 1"
            .to_string());
    }

    let dx =
        (b - a) / (nx as f64 - 1.0);

    let dt = (t_range.1 - t_range.0)
        / (nt as f64);

    let r = alpha * dt / (dx * dx);

    let interior = nx - 2;

    let a_diag =
        vec![1.0 + r; interior];

    let a_lower =
        vec![-r / 2.0; interior - 1];

    let a_upper =
        vec![-r / 2.0; interior - 1];

    let b_diag =
        vec![1.0 - r; interior];

    let b_lower =
        vec![r / 2.0; interior - 1];

    let b_upper =
        vec![r / 2.0; interior - 1];

    let mut u0 = vec![0.0_f64; nx];

    for (i, var) in u0
        .iter_mut()
        .enumerate()
        .take(nx)
    {

        let x =
            (i as f64).mul_add(dx, a);

        *var = init_func(x);
    }

    u0[0] = 0.0;

    u0[nx - 1] = 0.0;

    let mut results: Vec<Vec<f64>> =
        Vec::with_capacity(nt + 1);

    results.push(u0.clone());

    let solve_tridiag =
        |a_l: Vec<f64>,
         mut a_d: Vec<f64>,
         a_u: Vec<f64>,
         mut d: Vec<f64>|
         -> Vec<f64> {

            let n = a_d.len();

            for i in 1 .. n {

                let m = a_l[i - 1]
                    / a_d[i - 1];

                a_d[i] -=
                    m * a_u[i - 1];

                d[i] -= m * d[i - 1];
            }

            let mut x =
                vec![0.0_f64; n];

            x[n - 1] =
                d[n - 1] / a_d[n - 1];

            for i in (0 .. n - 1).rev()
            {

                x[i] = a_u[i].mul_add(
                    -x[i + 1],
                    d[i],
                ) / a_d[i];
            }

            x
        };

    for _step in 0 .. nt {

        let mut d =
            vec![0.0_f64; interior];

        for i in 0 .. interior {

            let global_i = i + 1;

            let mut val = b_diag[i]
                * u0[global_i];

            if i > 0 {

                val += b_lower[i - 1]
                    * u0[global_i - 1];
            }

            if i + 1 < interior {

                val += b_upper[i]
                    * u0[global_i + 1];
            }

            d[i] = val;
        }

        let u_interior = solve_tridiag(
            a_lower.clone(),
            a_diag.clone(),
            a_upper.clone(),
            d,
        );

        let mut u_new =
            vec![0.0_f64; nx];

        u_new[0] = 0.0;

        u_new[1 ..= interior]
            .copy_from_slice(
                &u_interior
                    [.. interior],
            );

        u_new[nx - 1] = 0.0;

        results.push(u_new.clone());

        u0 = u_new;
    }

    Ok(results)
}

/// Solves the 1D Wave Equation `u_tt = c^2 * u_xx` using an explicit finite-difference method.
///
/// This function implements a basic explicit finite-difference scheme for the 1D wave equation.
/// It requires initial displacement and velocity profiles and applies Dirichlet boundary conditions.
/// Stability is governed by the CFL condition (`c*dt/dx <= 1`).
///
/// # Arguments
/// * `initial_u` - Initial displacement vector (length `num_points`).
/// * `initial_ut` - Initial velocity vector (length `num_points`).
/// * `c` - Wave speed.
/// * `dx` - Spatial step size (uniform).
/// * `dt` - Time step size (must satisfy CFL condition).
/// * `num_steps` - Number of time steps to march.
///
/// # Returns
/// A `Result` containing a `Vec<Vec<f64>>` where each inner `Vec` is the solution `u`
/// at a given time step, or an error string if input dimensions mismatch or CFL is violated.
///
/// # Errors
/// Returns an error if initial vectors have different lengths, if `num_points` < 3, or if the CFL stability condition is violated.

pub fn solve_wave_equation_1d(
    initial_u: &[f64],
    initial_ut: &[f64],
    c: f64,
    dx: f64,
    dt: f64,
    num_steps: usize,
) -> Result<Vec<Vec<f64>>, String> {

    let n = initial_u.len();

    if initial_ut.len() != n {

        return Err("initial_ut \
                    length mismatch"
            .to_string());
    }

    if n < 3 {

        return Err("need at least 3 \
                    grid points"
            .to_string());
    }

    let cfl = c * dt / dx;

    if cfl.abs() > 1.0 {

        return Err(format!(
            "CFL violation: c*dt/dx = \
             {cfl} > 1. Reduce dt or \
             increase dx."
        ));
    }

    let mut u_prev =
        initial_u.to_owned();

    let mut u_curr = vec![0.0; n];

    for i in 1 .. (n - 1) {

        let u_xx = (2.0f64.mul_add(
            -u_prev[i],
            u_prev[i - 1],
        ) + u_prev[i + 1])
            / (dx * dx);

        u_curr[i] =
            (0.5 * (c * c) * (dt * dt))
                .mul_add(
                    u_xx,
                    dt.mul_add(
                        initial_ut[i],
                        u_prev[i],
                    ),
                );
    }

    u_curr[0] = 0.0;

    u_curr[n - 1] = 0.0;

    let mut snapshots: Vec<Vec<f64>> =
        Vec::with_capacity(
            num_steps + 1,
        );

    snapshots.push(u_prev.clone());

    snapshots.push(u_curr.clone());

    for _step in 2 ..= num_steps {

        let mut u_next = vec![0.0; n];

        for i in 1 .. (n - 1) {

            let u_xx =
                (2.0f64.mul_add(
                    -u_curr[i],
                    u_curr[i - 1],
                ) + u_curr[i + 1])
                    / (dx * dx);

            u_next[i] = ((c * c)
                * (dt * dt))
                .mul_add(
                    u_xx,
                    2.0f64.mul_add(
                        u_curr[i],
                        -u_prev[i],
                    ),
                );
        }

        u_next[0] = 0.0;

        u_next[n - 1] = 0.0;

        snapshots.push(u_next.clone());

        u_prev = u_curr;

        u_curr = u_next;
    }

    Ok(snapshots)
}

// ============================================================================
// Additional Classical Mechanics
// ============================================================================

/// Parameters for projectile motion simulation.
#[derive(
    Debug,
    Clone,
    Copy,
    Serialize,
    Deserialize,
)]

pub struct ProjectileParams {
    /// Initial velocity (m/s)
    pub v0: f64,
    /// Launch angle (radians)
    pub angle: f64,
    /// Mass of the projectile (kg)
    pub mass: f64,
    /// Drag coefficient (dimensionless)
    pub drag_coeff: f64,
    /// Cross-sectional area (m²)
    pub area: f64,
    /// Air density (kg/m³)
    pub air_density: f64,
    /// Time step (s)
    pub dt: f64,
    /// Maximum simulation time (s)
    pub max_time: f64,
}

/// Simulates projectile motion with air resistance (drag).
///
/// Uses Euler integration to solve equations of motion:
/// F_drag = -0.5 * ρ * A * Cd * v² * v_hat
///
/// # Arguments
/// * `params` - Simulation parameters
///
/// # Returns
/// Vector of (time, x, y, vx, vy) tuples
#[must_use]

pub fn projectile_motion_with_drag(
    params: ProjectileParams
) -> Vec<(
    f64,
    f64,
    f64,
    f64,
    f64,
)> {

    let ProjectileParams {
        v0,
        angle,
        mass,
        drag_coeff,
        area,
        air_density,
        dt,
        max_time,
    } = params;

    let mut results = Vec::new();

    let mut t = 0.0;

    let mut x = 0.0;

    let mut y = 0.0;

    let mut vx = v0 * angle.cos();

    let mut vy = v0 * angle.sin();

    let k = 0.5
        * drag_coeff
        * air_density
        * area
        / mass;

    while t <= max_time && y >= 0.0 {

        results.push((t, x, y, vx, vy));

        let v = vx.hypot(vy);

        if v > 1e-10 {

            let ax = -k * v * vx;

            let ay = (k * v).mul_add(
                -vy,
                -STANDARD_GRAVITY,
            );

            vx += ax * dt;

            vy += ay * dt;
        } else {

            vy -= STANDARD_GRAVITY * dt;
        }

        x += vx * dt;

        y += vy * dt;

        t += dt;
    }

    results
}

/// Simple harmonic oscillator solution.
///
/// Returns position x(t) = A * cos(ωt + φ)
///
/// # Arguments
/// * `amplitude` - Amplitude A
/// * `omega` - Angular frequency ω (rad/s)
/// * `phase` - Initial phase φ (radians)
/// * `time` - Time t
#[must_use]

pub fn simple_harmonic_oscillator(
    amplitude: f64,
    omega: f64,
    phase: f64,
    time: f64,
) -> f64 {

    amplitude
        * omega
            .mul_add(time, phase)
            .cos()
}

/// Damped harmonic oscillator solution.
///
/// Returns position for underdamped case: x(t) = A * e^(-γt) * cos(ω't + φ)
/// where ω' = √(ω₀² - γ²)
///
/// # Arguments
/// * `amplitude` - Initial amplitude
/// * `omega0` - Natural angular frequency ω₀
/// * `gamma` - Damping coefficient γ
/// * `phase` - Initial phase
/// * `time` - Time
#[must_use]

pub fn damped_harmonic_oscillator(
    amplitude: f64,
    omega0: f64,
    gamma: f64,
    phase: f64,
    time: f64,
) -> f64 {

    let omega_sq = omega0.mul_add(
        omega0,
        -(gamma * gamma),
    );

    if omega_sq > 0.0 {

        // Underdamped
        let omega_prime =
            omega_sq.sqrt();

        amplitude
            * (-gamma * time).exp()
            * omega_prime
                .mul_add(time, phase)
                .cos()
    } else if omega_sq < 0.0 {

        // Overdamped
        let beta = (-omega_sq).sqrt();

        amplitude
            * (-gamma * time).exp()
            * ((-beta * time).exp()
                + (beta * time).exp())
            / 2.0
    } else {

        // Critically damped
        amplitude
            * gamma.mul_add(time, 1.0)
            * (-gamma * time).exp()
    }
}

/// Simulates N-body gravitational dynamics using the leapfrog integrator.
///
/// # Arguments
/// * `particles` - Vector of particles with mass, position, and velocity
/// * `dt` - Time step
/// * `num_steps` - Number of steps to simulate
/// * `g` - Gravitational constant (default: use `GRAVITATIONAL_CONSTANT`)
///
/// # Returns
/// Vector of particle state snapshots at each time step
#[must_use]

pub fn simulate_n_body(
    mut particles: Vec<Particle3D>,
    dt: f64,
    num_steps: usize,
    g: f64,
) -> Vec<Vec<Particle3D>> {

    let n = particles.len();

    let mut snapshots =
        Vec::with_capacity(num_steps);

    snapshots.push(particles.clone());

    for _ in 0 .. num_steps {

        // Compute accelerations
        let mut ax = vec![0.0; n];

        let mut ay = vec![0.0; n];

        let mut az = vec![0.0; n];

        for i in 0 .. n {

            for j in 0 .. n {

                if i != j {

                    let dx = particles
                        [j]
                        .x
                        - particles[i]
                            .x;

                    let dy = particles
                        [j]
                        .y
                        - particles[i]
                            .y;

                    let dz = particles
                        [j]
                        .z
                        - particles[i]
                            .z;

                    let r_sq =
                        dz.mul_add(
                            dz,
                            dx.mul_add(
                                dx,
                                dy * dy,
                            ),
                        ) + 1e-10; // Softening
                    let r = r_sq.sqrt();

                    let f = g
                        * particles[j]
                            .mass
                        / (r_sq * r);

                    ax[i] += f * dx;

                    ay[i] += f * dy;

                    az[i] += f * dz;
                }
            }
        }

        // Leapfrog integration
        for i in 0 .. n {

            particles[i].vx +=
                ax[i] * dt;

            particles[i].vy +=
                ay[i] * dt;

            particles[i].vz +=
                az[i] * dt;

            particles[i].x +=
                particles[i].vx * dt;

            particles[i].y +=
                particles[i].vy * dt;

            particles[i].z +=
                particles[i].vz * dt;
        }

        snapshots
            .push(particles.clone());
    }

    snapshots
}

/// Calculates the total gravitational potential energy of a system of particles.
#[must_use]

pub fn gravitational_potential_energy(
    particles: &[Particle3D],
    g: f64,
) -> f64 {

    let n = particles.len();

    let mut energy = 0.0;

    for i in 0 .. n {

        for j in (i + 1) .. n {

            let dx = particles[j].x
                - particles[i].x;

            let dy = particles[j].y
                - particles[i].y;

            let dz = particles[j].z
                - particles[i].z;

            let r = dz
                .mul_add(
                    dz,
                    dx.mul_add(
                        dx,
                        dy * dy,
                    ),
                )
                .sqrt();

            if r > 1e-10 {

                energy -= g
                    * particles[i].mass
                    * particles[j].mass
                    / r;
            }
        }
    }

    energy
}

/// Calculates the total kinetic energy of a system of particles.
#[must_use]

pub fn total_kinetic_energy(
    particles: &[Particle3D]
) -> f64 {

    particles
        .iter()
        .map(Particle3D::kinetic_energy)
        .sum()
}

// ============================================================================
// Electromagnetism
// ============================================================================

/// Calculates the electric force between two point charges using Coulomb's law.
///
/// F = k * q1 * q2 / r²
///
/// # Arguments
/// * `q1`, `q2` - Charges (C)
/// * `r` - Distance between charges (m)
///
/// # Returns force magnitude (N), positive for repulsion
#[must_use]

pub fn coulomb_force(
    q1: f64,
    q2: f64,
    r: f64,
) -> f64 {

    if r.abs() < 1e-15 {

        return f64::INFINITY;
    }

    COULOMB_CONSTANT * q1 * q2 / (r * r)
}

/// Calculates the electric field magnitude from a point charge.
///
/// E = k * q / r²
///
/// # Arguments
/// * `q` - Charge (C)
/// * `r` - Distance from charge (m)
#[must_use]

pub fn electric_field_point_charge(
    q: f64,
    r: f64,
) -> f64 {

    if r.abs() < 1e-15 {

        return f64::INFINITY;
    }

    COULOMB_CONSTANT * q.abs() / (r * r)
}

/// Calculates the electric potential from a point charge.
///
/// V = k * q / r
///
/// # Arguments
/// * `q` - Charge (C)
/// * `r` - Distance from charge (m)
#[must_use]

pub fn electric_potential_point_charge(
    q: f64,
    r: f64,
) -> f64 {

    if r.abs() < 1e-15 {

        return f64::INFINITY
            * q.signum();
    }

    COULOMB_CONSTANT * q / r
}

/// Calculates the magnetic field magnitude at distance r from an infinite wire.
///
/// B = μ₀ * I / (2π * r)
///
/// # Arguments
/// * `current` - Current in wire (A)
/// * `r` - Perpendicular distance from wire (m)
#[must_use]

pub fn magnetic_field_infinite_wire(
    current: f64,
    r: f64,
) -> f64 {

    if r.abs() < 1e-15 {

        return f64::INFINITY;
    }

    VACUUM_PERMEABILITY * current.abs()
        / (2.0
            * std::f64::consts::PI
            * r)
}

/// Calculates the Lorentz force on a charged particle.
///
/// F = q(E + v × B) - returns force magnitude assuming v ⟂ B
///
/// # Arguments
/// * `charge` - Particle charge (C)
/// * `velocity` - Particle speed (m/s)
/// * `e_field` - Electric field magnitude (V/m)
/// * `b_field` - Magnetic field magnitude (T)
#[must_use]

pub fn lorentz_force(
    charge: f64,
    velocity: f64,
    e_field: f64,
    b_field: f64,
) -> f64 {

    charge.abs()
        * velocity
            .mul_add(b_field, e_field)
}

/// Calculates the cyclotron radius (Larmor radius) for a charged particle.
///
/// r = mv / (|q|B)
///
/// # Arguments
/// * `mass` - Particle mass (kg)
/// * `velocity` - Particle speed (m/s)
/// * `charge` - Particle charge (C)
/// * `b_field` - Magnetic field magnitude (T)
#[must_use]

pub fn cyclotron_radius(
    mass: f64,
    velocity: f64,
    charge: f64,
    b_field: f64,
) -> f64 {

    if charge.abs() < 1e-30
        || b_field.abs() < 1e-30
    {

        return f64::INFINITY;
    }

    mass * velocity
        / (charge.abs() * b_field)
}

// ============================================================================
// Thermodynamics
// ============================================================================

/// Ideal gas law: PV = nRT, solving for pressure.
///
/// # Arguments
/// * `n` - Amount of substance (mol)
/// * `t` - Temperature (K)
/// * `v` - Volume (m³)
#[must_use]

pub fn ideal_gas_pressure(
    n: f64,
    t: f64,
    v: f64,
) -> f64 {

    if v.abs() < 1e-30 {

        return f64::INFINITY;
    }

    n * GAS_CONSTANT * t / v
}

/// Ideal gas law: PV = nRT, solving for volume.
#[must_use]

pub fn ideal_gas_volume(
    n: f64,
    t: f64,
    p: f64,
) -> f64 {

    if p.abs() < 1e-30 {

        return f64::INFINITY;
    }

    n * GAS_CONSTANT * t / p
}

/// Ideal gas law: PV = nRT, solving for temperature.
#[must_use]

pub fn ideal_gas_temperature(
    p: f64,
    v: f64,
    n: f64,
) -> f64 {

    if n.abs() < 1e-30 {

        return f64::INFINITY;
    }

    p * v / (n * GAS_CONSTANT)
}

/// Maxwell-Boltzmann speed distribution probability density.
///
/// f(v) = 4π * (m/(2πkT))^(3/2) * v² * exp(-mv²/(2kT))
///
/// # Arguments
/// * `v` - Speed (m/s)
/// * `mass` - Particle mass (kg)
/// * `temperature` - Temperature (K)
#[must_use]

pub fn maxwell_boltzmann_speed_distribution(
    v: f64,
    mass: f64,
    temperature: f64,
) -> f64 {

    let kt = BOLTZMANN_CONSTANT
        * temperature;

    if kt < 1e-30 {

        return 0.0;
    }

    let a = mass
        / (2.0
            * std::f64::consts::PI
            * kt);

    4.0 * std::f64::consts::PI
        * a.powf(1.5)
        * v
        * v
        * (-mass * v * v / (2.0 * kt))
            .exp()
}

/// Mean speed from Maxwell-Boltzmann distribution.
///
/// ⟨v⟩ = √(8kT/(πm))
#[must_use]

pub fn maxwell_boltzmann_mean_speed(
    mass: f64,
    temperature: f64,
) -> f64 {

    (8.0 * BOLTZMANN_CONSTANT
        * temperature
        / (std::f64::consts::PI * mass))
        .sqrt()
}

/// RMS speed from Maxwell-Boltzmann distribution.
///
/// `v_rms` = √(3kT/m)
#[must_use]

pub fn maxwell_boltzmann_rms_speed(
    mass: f64,
    temperature: f64,
) -> f64 {

    (3.0 * BOLTZMANN_CONSTANT
        * temperature
        / mass)
        .sqrt()
}

/// Stefan-Boltzmann law: total power radiated by a blackbody.
///
/// P = σ * A * T⁴
///
/// # Arguments
/// * `area` - Surface area (m²)
/// * `temperature` - Temperature (K)
#[must_use]

pub fn blackbody_power(
    area: f64,
    temperature: f64,
) -> f64 {

    STEFAN_BOLTZMANN
        * area
        * temperature.powi(4)
}

/// Wien's displacement law: peak wavelength of blackbody radiation.
///
/// `λ_max` = b / T, where b ≈ 2.898 × 10⁻³ m·K
#[must_use]

pub fn wien_displacement_wavelength(
    temperature: f64
) -> f64 {

    const WIEN_CONSTANT: f64 =
        2.897_771_955e-3;

    if temperature < 1e-10 {

        return f64::INFINITY;
    }

    WIEN_CONSTANT / temperature
}

// ============================================================================
// Special Relativity
// ============================================================================

/// Lorentz factor γ = 1 / √(1 - v²/c²)
///
/// # Arguments
/// * `velocity` - Velocity (m/s)
#[must_use]

pub fn lorentz_factor(
    velocity: f64
) -> f64 {

    let beta =
        velocity / SPEED_OF_LIGHT;

    if beta.abs() >= 1.0 {

        return f64::INFINITY;
    }

    1.0 / beta
        .mul_add(-beta, 1.0)
        .sqrt()
}

/// Time dilation: Δt = γ * Δt₀
///
/// # Arguments
/// * `proper_time` - Proper time interval (s)
/// * `velocity` - Relative velocity (m/s)
#[must_use]

pub fn time_dilation(
    proper_time: f64,
    velocity: f64,
) -> f64 {

    proper_time
        * lorentz_factor(velocity)
}

/// Length contraction: L = L₀ / γ
///
/// # Arguments
/// * `proper_length` - Proper length (m)
/// * `velocity` - Relative velocity (m/s)
#[must_use]

pub fn length_contraction(
    proper_length: f64,
    velocity: f64,
) -> f64 {

    proper_length
        / lorentz_factor(velocity)
}

/// Relativistic momentum: p = γmv
///
/// # Arguments
/// * `mass` - Rest mass (kg)
/// * `velocity` - Velocity (m/s)
#[must_use]

pub fn relativistic_momentum(
    mass: f64,
    velocity: f64,
) -> f64 {

    lorentz_factor(velocity)
        * mass
        * velocity
}

/// Relativistic kinetic energy: KE = (γ - 1)mc²
///
/// # Arguments
/// * `mass` - Rest mass (kg)
/// * `velocity` - Velocity (m/s)
#[must_use]

pub fn relativistic_kinetic_energy(
    mass: f64,
    velocity: f64,
) -> f64 {

    (lorentz_factor(velocity) - 1.0)
        * mass
        * SPEED_OF_LIGHT
        * SPEED_OF_LIGHT
}

/// Total relativistic energy: E = γmc²
///
/// # Arguments
/// * `mass` - Rest mass (kg)
/// * `velocity` - Velocity (m/s)
#[must_use]

pub fn relativistic_total_energy(
    mass: f64,
    velocity: f64,
) -> f64 {

    lorentz_factor(velocity)
        * mass
        * SPEED_OF_LIGHT
        * SPEED_OF_LIGHT
}

/// Mass-energy equivalence: E = mc²
#[must_use]

pub fn mass_energy(mass: f64) -> f64 {

    mass * SPEED_OF_LIGHT
        * SPEED_OF_LIGHT
}

/// Relativistic velocity addition: u = (v + w) / (1 + vw/c²)
#[must_use]
#[allow(clippy::suspicious_operation_groupings)]

pub fn relativistic_velocity_addition(
    v: f64,
    w: f64,
) -> f64 {

    (v + w)
        / (1.0
            + v * w
                / (SPEED_OF_LIGHT
                    * SPEED_OF_LIGHT))
}

// ============================================================================
// Quantum Mechanics Helpers
// ============================================================================

/// Quantum harmonic oscillator energy levels.
///
/// `E_n` = ħω(n + 1/2)
///
/// # Arguments
/// * `n` - Quantum number (0, 1, 2, ...)
/// * `omega` - Angular frequency (rad/s)
#[must_use]

pub fn quantum_harmonic_oscillator_energy(
    n: u64,
    omega: f64,
) -> f64 {

    HBAR * omega * (n as f64 + 0.5)
}

/// Hydrogen atom energy levels (Bohr model).
///
/// `E_n` = -13.6 eV / n²
///
/// # Arguments
/// * `n` - Principal quantum number (1, 2, 3, ...)
///
/// # Returns energy in Joules
#[must_use]

pub fn hydrogen_energy_level(
    n: u64
) -> f64 {

    const RYDBERG_ENERGY_J: f64 =
        2.179_872_361e-18; // 13.6 eV in Joules

    if n == 0 {

        return f64::NEG_INFINITY;
    }

    -RYDBERG_ENERGY_J
        / (n as f64 * n as f64)
}

/// De Broglie wavelength: λ = h/p
///
/// # Arguments
/// * `momentum` - Particle momentum (kg·m/s)
#[must_use]

pub fn de_broglie_wavelength(
    momentum: f64
) -> f64 {

    if momentum.abs() < 1e-40 {

        return f64::INFINITY;
    }

    PLANCK_CONSTANT / momentum
}

/// Heisenberg uncertainty principle minimum: Δx·Δp ≥ ħ/2
///
/// Returns minimum position uncertainty given momentum uncertainty.
#[must_use]

pub fn heisenberg_position_uncertainty(
    momentum_uncertainty: f64
) -> f64 {

    if momentum_uncertainty.abs()
        < 1e-40
    {

        return 0.0;
    }

    HBAR / (2.0 * momentum_uncertainty)
}

/// Photon energy: E = hf = hc/λ
///
/// # Arguments
/// * `wavelength` - Wavelength (m)
#[must_use]

pub fn photon_energy(
    wavelength: f64
) -> f64 {

    if wavelength.abs() < 1e-20 {

        return f64::INFINITY;
    }

    PLANCK_CONSTANT * SPEED_OF_LIGHT
        / wavelength
}

/// Photon wavelength from energy.
#[must_use]

pub fn photon_wavelength(
    energy: f64
) -> f64 {

    if energy.abs() < 1e-40 {

        return f64::INFINITY;
    }

    PLANCK_CONSTANT * SPEED_OF_LIGHT
        / energy
}

/// Compton wavelength of a particle: `λ_C` = h/(mc)
#[must_use]

pub fn compton_wavelength(
    mass: f64
) -> f64 {

    if mass.abs() < 1e-40 {

        return f64::INFINITY;
    }

    PLANCK_CONSTANT
        / (mass * SPEED_OF_LIGHT)
}
