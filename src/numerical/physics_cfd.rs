//! # Numerical Computational Fluid Dynamics (CFD)
//!
//! This module provides numerical methods for Computational Fluid Dynamics (CFD),
//! the branch of fluid mechanics that uses numerical analysis to solve and analyze
//! problems involving fluid flows.
//!
//! ## Equations Solved
//!
//! ### Transport Equations
//! - **Advection equation**: `∂u/∂t + c·∇u = 0`
//! - **Diffusion equation**: `∂u/∂t = α∇²u`
//! - **Advection-diffusion equation**: `∂u/∂t + c·∇u = α∇²u`
//!
//! ### Incompressible Flow
//! - **Poisson equation for pressure**: `∇²p = f`
//! - **Stokes equations**: Creeping flow (Re << 1)
//! - **Navier-Stokes equations**: Full viscous flow
//!
//! ## Numerical Methods
//!
//! - **Finite Difference Method (FDM)**: Central, upwind, Lax-Wendroff schemes
//! - **Iterative Solvers**: Jacobi, Gauss-Seidel, SOR
//! - **Time Integration**: Explicit Euler, Runge-Kutta
//!
//! ## Features
//!
//! - Velocity field computation
//! - Pressure field solving
//! - Vorticity calculation
//! - Stream function computation
//! - Reynolds number estimation
//! - CFL condition checking
//!
//! ## Example
//!
//! ```rust
//! 
//! use rssn::numerical::physics_cfd::*;
//!
//! // Solve 1D advection
//! let u0 = vec![
//!     0.0, 0.0, 1.0, 1.0, 0.0, 0.0,
//! ];
//!
//! let results = solve_advection_1d(
//!     &u0, 1.0, 0.1, 0.01, 10,
//! );
//! ```

use serde::Deserialize;
use serde::Serialize;

use crate::numerical::matrix::Matrix;

// ============================================================================
// Fluid Properties
// ============================================================================

/// Common fluid properties at standard conditions.
#[derive(
    Debug,
    Clone,
    Copy,
    PartialEq,
    Serialize,
    Deserialize,
)]

pub struct FluidProperties {
    /// Density (kg/m³)
    pub density : f64,
    /// Dynamic viscosity (Pa·s)
    pub dynamic_viscosity : f64,
    /// Thermal conductivity (W/(m·K))
    pub thermal_conductivity : f64,
    /// Specific heat capacity (J/(kg·K))
    pub specific_heat : f64,
}

impl FluidProperties {
    /// Creates new fluid properties.
    #[must_use]

    pub const fn new(
        density : f64,
        dynamic_viscosity : f64,
        thermal_conductivity : f64,
        specific_heat : f64,
    ) -> Self {

        Self {
            density,
            dynamic_viscosity,
            thermal_conductivity,
            specific_heat,
        }
    }

    /// Air at 20°C and 1 atm.
    #[must_use]

    pub const fn air() -> Self {

        Self {
            density : 1.204,
            dynamic_viscosity:
                1.825e-5,
            thermal_conductivity:
                0.0257,
            specific_heat : 1005.0,
        }
    }

    /// Water at 20°C.
    #[must_use]

    pub const fn water() -> Self {

        Self {
            density : 998.2,
            dynamic_viscosity:
                1.002e-3,
            thermal_conductivity:
                0.598,
            specific_heat : 4182.0,
        }
    }

    /// Kinematic viscosity ν = μ/ρ (m²/s)
    #[must_use]

    pub fn kinematic_viscosity(
        &self
    ) -> f64 {

        self.dynamic_viscosity
            / self.density
    }

    /// Thermal diffusivity α = k/(ρ·Cp) (m²/s)
    #[must_use]

    pub fn thermal_diffusivity(
        &self
    ) -> f64 {

        self.thermal_conductivity
            / (self.density
                * self.specific_heat)
    }

    /// Prandtl number Pr = ν/α = μ·Cp/k
    #[must_use]

    pub fn prandtl_number(
        &self
    ) -> f64 {

        self.dynamic_viscosity
            * self.specific_heat
            / self.thermal_conductivity
    }
}

// ============================================================================
// Dimensionless Numbers
// ============================================================================

/// Calculates the Reynolds number: Re = ρVL/μ = VL/ν
///
/// # Arguments
/// * `velocity` - Characteristic velocity (m/s)
/// * `length` - Characteristic length (m)
/// * `kinematic_viscosity` - Kinematic viscosity ν (m²/s)
#[must_use]

pub fn reynolds_number(
    velocity : f64,
    length : f64,
    kinematic_viscosity : f64,
) -> f64 {

    velocity * length
        / kinematic_viscosity
}

/// Calculates the Mach number: Ma = V/c
///
/// # Arguments
/// * `velocity` - Flow velocity (m/s)
/// * `speed_of_sound` - Speed of sound in the medium (m/s)
#[must_use]

pub fn mach_number(
    velocity : f64,
    speed_of_sound : f64,
) -> f64 {

    velocity / speed_of_sound
}

/// Calculates the Froude number: Fr = V/√(gL)
///
/// # Arguments
/// * `velocity` - Flow velocity (m/s)
/// * `length` - Characteristic length (m)
/// * `gravity` - Gravitational acceleration (m/s²)
#[must_use]

pub fn froude_number(
    velocity : f64,
    length : f64,
    gravity : f64,
) -> f64 {

    velocity / (gravity * length).sqrt()
}

/// Calculates the CFL (Courant-Friedrichs-Lewy) number.
///
/// For stability, typically CFL ≤ 1 for explicit methods.
///
/// # Arguments
/// * `velocity` - Flow velocity (m/s)
/// * `dt` - Time step (s)
/// * `dx` - Grid spacing (m)
#[must_use]

pub fn cfl_number(
    velocity : f64,
    dt : f64,
    dx : f64,
) -> f64 {

    velocity.abs() * dt / dx
}

/// Checks if the CFL condition is satisfied for stability.
#[must_use]

pub fn check_cfl_stability(
    velocity : f64,
    dt : f64,
    dx : f64,
    max_cfl : f64,
) -> bool {

    cfl_number(velocity, dt, dx)
        <= max_cfl
}

/// Calculates the diffusion number for stability analysis.
/// For stability, typically r ≤ 0.5 for 1D explicit diffusion.
#[must_use]

pub fn diffusion_number(
    alpha : f64,
    dt : f64,
    dx : f64,
) -> f64 {

    alpha * dt / (dx * dx)
}

// ============================================================================
// 1D Solvers
// ============================================================================

/// Solves the 1D advection equation `du/dt + c * du/dx = 0` using an explicit finite difference scheme.
///
/// This function implements a simple first-order upwind scheme for stability.
///
/// # Arguments
/// * `u0` - Initial condition (vector of `u` values at `t=0`).
/// * `c` - Advection speed.
/// * `dx` - Spatial step size.
/// * `dt` - Time step size.
/// * `num_steps` - Number of time steps to simulate.
///
/// # Returns
/// A `Vec<Vec<f64>>` where each inner `Vec` is the solution `u` at a given time step.
#[must_use]

pub fn solve_advection_1d(
    u0 : &[f64],
    c : f64,
    dx : f64,
    dt : f64,
    num_steps : usize,
) -> Vec<Vec<f64>> {

    let n = u0.len();

    let mut u = u0.to_vec();

    let mut results =
        Vec::with_capacity(
            num_steps + 1,
        );

    results.push(u.clone());

    let nu = c * dt / dx;

    for _ in 0 .. num_steps {

        let mut u_next = vec![0.0; n];

        for i in 1 .. (n - 1) {

            if c > 0.0 {

                u_next[i] = nu.mul_add(
                    -(u[i] - u[i - 1]),
                    u[i],
                );
            } else {

                u_next[i] = nu.mul_add(
                    -(u[i + 1] - u[i]),
                    u[i],
                );
            }
        }

        u_next[0] = u_next[n - 2];

        u_next[n - 1] = u_next[1];

        u = u_next;

        results.push(u.clone());
    }

    results
}

/// Solves the 1D diffusion equation `du/dt = alpha * d2u/dx2` using an explicit finite difference scheme.
///
/// # Arguments
/// * `u0` - Initial condition.
/// * `alpha` - Diffusion coefficient.
/// * `dx` - Spatial step size.
/// * `dt` - Time step size.
/// * `num_steps` - Number of time steps.
///
/// # Returns
/// A `Vec<Vec<f64>>` where each inner `Vec` is the solution `u` at a given time step.
#[must_use]

pub fn solve_diffusion_1d(
    u0 : &[f64],
    alpha : f64,
    dx : f64,
    dt : f64,
    num_steps : usize,
) -> Vec<Vec<f64>> {

    let n = u0.len();

    let mut u = u0.to_vec();

    let mut results =
        Vec::with_capacity(
            num_steps + 1,
        );

    results.push(u.clone());

    let r = alpha * dt / (dx * dx);

    for _ in 0 .. num_steps {

        let mut u_next = vec![0.0; n];

        u_next[0] = u[0];

        u_next[n - 1] = u[n - 1];

        for i in 1 .. (n - 1) {

            u_next[i] = r.mul_add(
                2.0f64.mul_add(
                    -u[i],
                    u[i - 1],
                ) + u[i + 1],
                u[i],
            );
        }

        u = u_next;

        results.push(u.clone());
    }

    results
}

/// Solves the 2D Poisson equation `∇²u = f` using Jacobi iteration.
///
/// This function implements the Jacobi iterative method to solve the Poisson equation
/// on a 2D grid with Dirichlet boundary conditions (implicitly handled by the iteration
/// not updating boundary points). It is suitable for steady-state problems.
///
/// # Arguments
/// * `f` - Source term (2D grid) as a `Matrix<f64>`.
/// * `u0` - Initial guess for `u` (2D grid) as a `Matrix<f64>`.
/// * `dx`, `dy` - Grid spacing in x and y directions.
/// * `max_iter` - Maximum number of iterations.
/// * `tolerance` - Convergence tolerance for the maximum difference between successive iterations.
///
/// # Returns
/// A `Matrix<f64>` representing the solution `u`.
#[must_use]

pub fn solve_poisson_2d_jacobi(
    f : &Matrix<f64>,
    u0 : &Matrix<f64>,
    dx : f64,
    dy : f64,
    max_iter : usize,
    tolerance : f64,
) -> Matrix<f64> {

    let nx = u0.rows();

    let ny = u0.cols();

    let mut u = u0.clone();

    let mut u_new = u0.clone();

    let dx2 = dx * dx;

    let dy2 = dy * dy;

    for _iter in 0 .. max_iter {

        let mut max_diff = 0.0;

        for i in 1 .. (nx - 1) {

            for j in 1 .. (ny - 1) {

                let val = 0.5
                    * (dy2.mul_add(
                        u.get(i + 1, j)
                            + u.get(
                                i - 1,
                                j,
                            ),
                        dx2 * (u.get(
                            i,
                            j + 1,
                        ) + u
                            .get(
                                i,
                                j - 1,
                            )),
                    ) - (dx2
                        * dy2
                        * f.get(i, j)))
                    / (dx2 + dy2);

                let diff = (val
                    - u.get(i, j))
                .abs();

                if diff > max_diff {

                    max_diff = diff;
                }

                *u_new.get_mut(i, j) =
                    val;
            }
        }

        u = u_new.clone();

        if max_diff < tolerance {

            break;
        }
    }

    u
}

/// Solves the 2D Poisson equation using Gauss-Seidel iteration.
/// Generally converges faster than Jacobi.
#[must_use]

pub fn solve_poisson_2d_gauss_seidel(
    f : &Matrix<f64>,
    u0 : &Matrix<f64>,
    dx : f64,
    dy : f64,
    max_iter : usize,
    tolerance : f64,
) -> Matrix<f64> {

    let nx = u0.rows();

    let ny = u0.cols();

    let mut u = u0.clone();

    let dx2 = dx * dx;

    let dy2 = dy * dy;

    let factor = 2.0 * (dx2 + dy2);

    for _iter in 0 .. max_iter {

        let mut max_diff = 0.0;

        for i in 1 .. (nx - 1) {

            for j in 1 .. (ny - 1) {

                let old_val =
                    *u.get(i, j);

                let new_val =
                    (dy2.mul_add(
                        *u.get(
                            i + 1,
                            j,
                        ) + *u.get(
                            i - 1,
                            j,
                        ),
                        dx2 * (*u.get(
                            i,
                            j + 1,
                        ) + *u
                            .get(
                                i,
                                j - 1,
                            )),
                    ) - dx2
                        * dy2
                        * *f.get(i, j))
                        / factor;

                *u.get_mut(i, j) =
                    new_val;

                let diff = (new_val
                    - old_val)
                    .abs();

                if diff > max_diff {

                    max_diff = diff;
                }
            }
        }

        if max_diff < tolerance {

            break;
        }
    }

    u
}

/// Solves the 2D Poisson equation using Successive Over-Relaxation (SOR).
///
/// # Arguments
/// * `omega` - Relaxation parameter (1 < ω < 2 for over-relaxation)
#[must_use]

pub fn solve_poisson_2d_sor(
    f : &Matrix<f64>,
    u0 : &Matrix<f64>,
    dx : f64,
    dy : f64,
    omega : f64,
    max_iter : usize,
    tolerance : f64,
) -> Matrix<f64> {

    let nx = u0.rows();

    let ny = u0.cols();

    let mut u = u0.clone();

    let dx2 = dx * dx;

    let dy2 = dy * dy;

    let factor = 2.0 * (dx2 + dy2);

    for _iter in 0 .. max_iter {

        let mut max_diff = 0.0;

        for i in 1 .. (nx - 1) {

            for j in 1 .. (ny - 1) {

                let old_val =
                    *u.get(i, j);

                let gs_val =
                    (dy2.mul_add(
                        *u.get(
                            i + 1,
                            j,
                        ) + *u.get(
                            i - 1,
                            j,
                        ),
                        dx2 * (*u.get(
                            i,
                            j + 1,
                        ) + *u
                            .get(
                                i,
                                j - 1,
                            )),
                    ) - dx2
                        * dy2
                        * *f.get(i, j))
                        / factor;

                let new_val = omega
                    .mul_add(
                        gs_val
                            - old_val,
                        old_val,
                    );

                *u.get_mut(i, j) =
                    new_val;

                let diff = (new_val
                    - old_val)
                    .abs();

                if diff > max_diff {

                    max_diff = diff;
                }
            }
        }

        if max_diff < tolerance {

            break;
        }
    }

    u
}

// ============================================================================
// Advection-Diffusion
// ============================================================================

/// Solves the 1D advection-diffusion equation: `∂u/∂t + c·∂u/∂x = α·∂²u/∂x²`
#[must_use]

pub fn solve_advection_diffusion_1d(
    u0 : &[f64],
    c : f64,
    alpha : f64,
    dx : f64,
    dt : f64,
    num_steps : usize,
) -> Vec<Vec<f64>> {

    let n = u0.len();

    let mut u = u0.to_vec();

    let mut results =
        Vec::with_capacity(
            num_steps + 1,
        );

    results.push(u.clone());

    let nu = c * dt / dx;

    let r = alpha * dt / (dx * dx);

    for _ in 0 .. num_steps {

        let mut u_next = vec![0.0; n];

        u_next[0] = u[0];

        u_next[n - 1] = u[n - 1];

        for i in 1 .. (n - 1) {

            // Upwind for advection + central for diffusion
            let advection = if c > 0.0 {

                -nu * (u[i] - u[i - 1])
            } else {

                -nu * (u[i + 1] - u[i])
            };

            let diffusion = r
                * (2.0f64.mul_add(
                    -u[i],
                    u[i + 1],
                ) + u[i - 1]);

            u_next[i] = u[i]
                + advection
                + diffusion;
        }

        u = u_next;

        results.push(u.clone());
    }

    results
}

// ============================================================================
// Burgers' Equation
// ============================================================================

/// Solves the 1D viscous Burgers' equation: `∂u/∂t + u·∂u/∂x = ν·∂²u/∂x²`
///
/// This nonlinear equation models shock formation and viscous dissipation.
#[must_use]

pub fn solve_burgers_1d(
    u0 : &[f64],
    nu : f64,
    dx : f64,
    dt : f64,
    num_steps : usize,
) -> Vec<Vec<f64>> {

    let n = u0.len();

    let mut u = u0.to_vec();

    let mut results =
        Vec::with_capacity(
            num_steps + 1,
        );

    results.push(u.clone());

    let r = nu * dt / (dx * dx);

    for _ in 0 .. num_steps {

        let mut u_next = vec![0.0; n];

        u_next[0] = u[0];

        u_next[n - 1] = u[n - 1];

        for i in 1 .. (n - 1) {

            // Nonlinear advection (conservative form with upwinding)
            let advection = if u[i]
                > 0.0
            {

                u[i] * (u[i] - u[i - 1])
                    / dx
            } else {

                u[i] * (u[i + 1] - u[i])
                    / dx
            };

            // Diffusion
            let diffusion = r
                * (2.0f64.mul_add(
                    -u[i],
                    u[i + 1],
                ) + u[i - 1]);

            u_next[i] = dt.mul_add(
                -advection,
                u[i],
            ) + diffusion;
        }

        u = u_next;

        results.push(u.clone());
    }

    results
}

// ============================================================================
// Vorticity and Stream Function
// ============================================================================

/// Computes the vorticity ω = ∂v/∂x - ∂u/∂y from velocity field.
///
/// # Arguments
/// * `u` - x-component of velocity (Matrix)
/// * `v` - y-component of velocity (Matrix)
/// * `dx`, `dy` - Grid spacing
#[must_use]

pub fn compute_vorticity(
    u : &Matrix<f64>,
    v : &Matrix<f64>,
    dx : f64,
    dy : f64,
) -> Matrix<f64> {

    let nx = u.rows();

    let ny = u.cols();

    let mut omega =
        Matrix::zeros(nx, ny);

    for i in 1 .. (nx - 1) {

        for j in 1 .. (ny - 1) {

            let dv_dx = (*v
                .get(i + 1, j)
                - *v.get(i - 1, j))
                / (2.0 * dx);

            let du_dy = (*u
                .get(i, j + 1)
                - *u.get(i, j - 1))
                / (2.0 * dy);

            *omega.get_mut(i, j) =
                dv_dx - du_dy;
        }
    }

    omega
}

/// Computes the stream function from vorticity using Poisson equation.
///
/// ∇²ψ = -ω
#[must_use]

pub fn compute_stream_function(
    omega : &Matrix<f64>,
    dx : f64,
    dy : f64,
    max_iter : usize,
    tolerance : f64,
) -> Matrix<f64> {

    let nx = omega.rows();

    let ny = omega.cols();

    // Negate vorticity for the Poisson equation
    let mut neg_omega =
        Matrix::zeros(nx, ny);

    for i in 0 .. nx {

        for j in 0 .. ny {

            *neg_omega.get_mut(i, j) =
                -*omega.get(i, j);
        }
    }

    let psi0 = Matrix::zeros(nx, ny);

    solve_poisson_2d_gauss_seidel(
        &neg_omega,
        &psi0,
        dx,
        dy,
        max_iter,
        tolerance,
    )
}

/// Computes velocity field from stream function.
///
/// u = ∂ψ/∂y, v = -∂ψ/∂x
#[must_use]

pub fn velocity_from_stream_function(
    psi : &Matrix<f64>,
    dx : f64,
    dy : f64,
) -> (
    Matrix<f64>,
    Matrix<f64>,
) {

    let nx = psi.rows();

    let ny = psi.cols();

    let mut u = Matrix::zeros(nx, ny);

    let mut v = Matrix::zeros(nx, ny);

    for i in 1 .. (nx - 1) {

        for j in 1 .. (ny - 1) {

            *u.get_mut(i, j) = (*psi
                .get(i, j + 1)
                - *psi.get(i, j - 1))
                / (2.0 * dy);

            *v.get_mut(i, j) = -(*psi
                .get(i + 1, j)
                - *psi.get(i - 1, j))
                / (2.0 * dx);
        }
    }

    (u, v)
}

// ============================================================================
// Pressure and Velocity Correction
// ============================================================================

/// Computes the divergence of a 2D velocity field.
///
/// div(V) = ∂u/∂x + ∂v/∂y
#[must_use]

pub fn compute_divergence(
    u : &Matrix<f64>,
    v : &Matrix<f64>,
    dx : f64,
    dy : f64,
) -> Matrix<f64> {

    let nx = u.rows();

    let ny = u.cols();

    let mut div = Matrix::zeros(nx, ny);

    for i in 1 .. (nx - 1) {

        for j in 1 .. (ny - 1) {

            let du_dx = (*u
                .get(i + 1, j)
                - *u.get(i - 1, j))
                / (2.0 * dx);

            let dv_dy = (*v
                .get(i, j + 1)
                - *v.get(i, j - 1))
                / (2.0 * dy);

            *div.get_mut(i, j) =
                du_dx + dv_dy;
        }
    }

    div
}

/// Computes the gradient of a scalar field.
///
/// Returns (∂p/∂x, ∂p/∂y)
#[must_use]

pub fn compute_gradient(
    p : &Matrix<f64>,
    dx : f64,
    dy : f64,
) -> (
    Matrix<f64>,
    Matrix<f64>,
) {

    let nx = p.rows();

    let ny = p.cols();

    let mut dp_dx =
        Matrix::zeros(nx, ny);

    let mut dp_dy =
        Matrix::zeros(nx, ny);

    for i in 1 .. (nx - 1) {

        for j in 1 .. (ny - 1) {

            *dp_dx.get_mut(i, j) = (*p
                .get(i + 1, j)
                - *p.get(i - 1, j))
                / (2.0 * dx);

            *dp_dy.get_mut(i, j) = (*p
                .get(i, j + 1)
                - *p.get(i, j - 1))
                / (2.0 * dy);
        }
    }

    (dp_dx, dp_dy)
}

/// Computes the Laplacian of a scalar field.
///
/// ∇²f = ∂²f/∂x² + ∂²f/∂y²
#[must_use]

pub fn compute_laplacian(
    f : &Matrix<f64>,
    dx : f64,
    dy : f64,
) -> Matrix<f64> {

    let nx = f.rows();

    let ny = f.cols();

    let mut lap = Matrix::zeros(nx, ny);

    let dx2 = dx * dx;

    let dy2 = dy * dy;

    for i in 1 .. (nx - 1) {

        for j in 1 .. (ny - 1) {

            let d2f_dx2 =
                (2.0f64.mul_add(
                    -*f.get(i, j),
                    *f.get(i + 1, j),
                ) + *f.get(i - 1, j))
                    / dx2;

            let d2f_dy2 =
                (2.0f64.mul_add(
                    -*f.get(i, j),
                    *f.get(i, j + 1),
                ) + *f.get(i, j - 1))
                    / dy2;

            *lap.get_mut(i, j) =
                d2f_dx2 + d2f_dy2;
        }
    }

    lap
}

// ============================================================================
// Lid-Driven Cavity Flow
// ============================================================================

/// Simplified lid-driven cavity flow solver using stream function-vorticity formulation.
///
/// # Arguments
/// * `nx`, `ny` - Grid dimensions
/// * `re` - Reynolds number
/// * `lid_velocity` - Top lid velocity
/// * `num_steps` - Number of time steps
/// * `dt` - Time step
///
/// # Returns
/// (stream function, vorticity)
#[must_use]

pub fn lid_driven_cavity_simple(
    nx : usize,
    ny : usize,
    re : f64,
    lid_velocity : f64,
    num_steps : usize,
    dt : f64,
) -> (
    Matrix<f64>,
    Matrix<f64>,
) {

    let dx = 1.0 / (nx - 1) as f64;

    let dy = 1.0 / (ny - 1) as f64;

    let nu = lid_velocity / re;

    let mut psi = Matrix::zeros(nx, ny);

    let mut omega =
        Matrix::zeros(nx, ny);

    // Set boundary condition for vorticity at lid
    for i in 0 .. nx {

        *omega.get_mut(i, ny - 1) = -2.0
            * *psi.get(i, ny - 2)
            / (dy * dy)
            - 2.0 * lid_velocity / dy;
    }

    for _ in 0 .. num_steps {

        // Solve for stream function
        psi = compute_stream_function(
            &omega,
            dx,
            dy,
            100,
            1e-6,
        );

        // Update vorticity
        let mut omega_new =
            omega.clone();

        for i in 1 .. (nx - 1) {

            for j in 1 .. (ny - 1) {

                // Advection (using stream function)
                let u = (*psi
                    .get(i, j + 1)
                    - *psi
                        .get(i, j - 1))
                    / (2.0 * dy);

                let v = -(*psi
                    .get(i + 1, j)
                    - *psi
                        .get(i - 1, j))
                    / (2.0 * dx);

                let domega_dx = (*omega
                    .get(i + 1, j)
                    - *omega
                        .get(i - 1, j))
                    / (2.0 * dx);

                let domega_dy = (*omega
                    .get(i, j + 1)
                    - *omega
                        .get(i, j - 1))
                    / (2.0 * dy);

                // Diffusion
                let d2omega_dx2 =
                    (2.0f64.mul_add(
                        -*omega
                            .get(i, j),
                        *omega.get(
                            i + 1,
                            j,
                        ),
                    ) + *omega
                        .get(i - 1, j))
                        / (dx * dx);

                let d2omega_dy2 =
                    (2.0f64.mul_add(
                        -*omega
                            .get(i, j),
                        *omega.get(
                            i,
                            j + 1,
                        ),
                    ) + *omega
                        .get(i, j - 1))
                        / (dy * dy);

                *omega_new.get_mut(i, j) = *omega.get(i, j)
                    + dt * nu.mul_add(d2omega_dx2 + d2omega_dy2, -(u * domega_dx + v * domega_dy));
            }
        }

        omega = omega_new;

        // Update boundary conditions
        for i in 0 .. nx {

            *omega.get_mut(i, ny - 1) =
                -2.0 * *psi
                    .get(i, ny - 2)
                    / (dy * dy)
                    - 2.0 * lid_velocity
                        / dy;

            *omega.get_mut(i, 0) = -2.0
                * *psi.get(i, 1)
                / (dy * dy);
        }

        for j in 0 .. ny {

            *omega.get_mut(0, j) = -2.0
                * *psi.get(1, j)
                / (dx * dx);

            *omega.get_mut(nx - 1, j) =
                -2.0 * *psi
                    .get(nx - 2, j)
                    / (dx * dx);
        }
    }

    (psi, omega)
}

// ============================================================================
// Utility Functions
// ============================================================================

/// Applies Dirichlet boundary conditions to a 2D field.

pub fn apply_dirichlet_bc(
    field : &mut Matrix<f64>,
    boundary_value : f64,
) {

    let nx = field.rows();

    let ny = field.cols();

    // Top and bottom
    for i in 0 .. nx {

        *field.get_mut(i, 0) =
            boundary_value;

        *field.get_mut(i, ny - 1) =
            boundary_value;
    }

    // Left and right
    for j in 0 .. ny {

        *field.get_mut(0, j) =
            boundary_value;

        *field.get_mut(nx - 1, j) =
            boundary_value;
    }
}

/// Applies Neumann boundary conditions (zero gradient) to a 2D field.

pub fn apply_neumann_bc(
    field : &mut Matrix<f64>
) {

    let nx = field.rows();

    let ny = field.cols();

    // Copy from interior to boundaries
    for i in 0 .. nx {

        *field.get_mut(i, 0) =
            *field.get(i, 1);

        *field.get_mut(i, ny - 1) =
            *field.get(i, ny - 2);
    }

    for j in 0 .. ny {

        *field.get_mut(0, j) =
            *field.get(1, j);

        *field.get_mut(nx - 1, j) =
            *field.get(nx - 2, j);
    }
}

/// Computes the maximum velocity magnitude from a 2D velocity field.
#[must_use]

pub fn max_velocity_magnitude(
    u : &Matrix<f64>,
    v : &Matrix<f64>,
) -> f64 {

    let nx = u.rows();

    let ny = u.cols();

    let mut max_vel = 0.0;

    for i in 0 .. nx {

        for j in 0 .. ny {

            let vel = (*u.get(i, j))
                .mul_add(
                    *u.get(i, j),
                    *v.get(i, j)
                        * *v.get(i, j),
                )
                .sqrt();

            if vel > max_vel {

                max_vel = vel;
            }
        }
    }

    max_vel
}

/// Computes the L2 norm of a field for convergence checking.
#[must_use]

pub fn l2_norm(
    field : &Matrix<f64>
) -> f64 {

    let nx = field.rows();

    let ny = field.cols();

    let mut sum = 0.0;

    for i in 0 .. nx {

        for j in 0 .. ny {

            sum += *field.get(i, j)
                * *field.get(i, j);
        }
    }

    (sum / (nx * ny) as f64).sqrt()
}

/// Computes the maximum absolute value in a field.
#[must_use]

pub fn max_abs(
    field : &Matrix<f64>
) -> f64 {

    let nx = field.rows();

    let ny = field.cols();

    let mut max_val = 0.0;

    for i in 0 .. nx {

        for j in 0 .. ny {

            let abs_val = field
                .get(i, j)
                .abs();

            if abs_val > max_val {

                max_val = abs_val;
            }
        }
    }

    max_val
}
