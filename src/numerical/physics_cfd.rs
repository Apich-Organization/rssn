//! # Numerical Computational Fluid Dynamics (CFD)
//!
//! This module provides numerical methods for Computational Fluid Dynamics (CFD).
//! It includes solvers for fundamental fluid dynamics equations like the Navier-Stokes
//! equations, typically using finite difference or finite volume methods.

use crate::numerical::matrix::Matrix;

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
pub fn solve_advection_1d(u0: &[f64], c: f64, dx: f64, dt: f64, num_steps: usize) -> Vec<Vec<f64>> {
    let n = u0.len();
    let mut u = u0.to_vec();
    let mut results = Vec::with_capacity(num_steps + 1);
    results.push(u.clone());

    let nu = c * dt / dx; // Courant number

    for _ in 0..num_steps {
        let mut u_next = vec![0.0; n];
        for i in 1..(n - 1) {
            // Upwind scheme for c > 0
            if c > 0.0 {
                u_next[i] = u[i] - nu * (u[i] - u[i - 1]);
            } else {
                // Downwind scheme for c < 0
                u_next[i] = u[i] - nu * (u[i + 1] - u[i]);
            }
        }
        // Boundary conditions (e.g., periodic or fixed)
        u_next[0] = u_next[n - 2]; // Simple periodic boundary
        u_next[n - 1] = u_next[1]; // Simple periodic boundary

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
pub fn solve_diffusion_1d(
    u0: &[f64],
    alpha: f64,
    dx: f64,
    dt: f64,
    num_steps: usize,
) -> Vec<Vec<f64>> {
    let n = u0.len();
    let mut u = u0.to_vec();
    let mut results = Vec::with_capacity(num_steps + 1);
    results.push(u.clone());

    let r = alpha * dt / (dx * dx); // Diffusion number

    for _ in 0..num_steps {
        let mut u_next = vec![0.0; n];
        // Boundary conditions (e.g., fixed at ends)
        u_next[0] = u[0];
        u_next[n - 1] = u[n - 1];

        for i in 1..(n - 1) {
            u_next[i] = u[i] + r * (u[i - 1] - 2.0 * u[i] + u[i + 1]);
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
pub fn solve_poisson_2d_jacobi(
    f: &Matrix<f64>,
    u0: &Matrix<f64>,
    dx: f64,
    dy: f64,
    max_iter: usize,
    tolerance: f64,
) -> Matrix<f64> {
    let nx = u0.rows();
    let ny = u0.cols();
    let mut u = u0.clone();
    let mut u_new = u0.clone();

    let dx2 = dx * dx;
    let dy2 = dy * dy;

    for _iter in 0..max_iter {
        let mut max_diff = 0.0;
        for i in 1..(nx - 1) {
            for j in 1..(ny - 1) {
                let val = 0.5
                    * ((dy2 * (u.get(i + 1, j) + u.get(i - 1, j)))
                        + (dx2 * (u.get(i, j + 1) + u.get(i, j - 1)))
                        - (dx2 * dy2 * f.get(i, j)))
                    / (dx2 + dy2);
                let diff = (val - u.get(i, j)).abs();
                if diff > max_diff {
                    max_diff = diff;
                }
                *u_new.get_mut(i, j) = val;
            }
        }
        u = u_new.clone();
        if max_diff < tolerance {
            break;
        }
    }
    u
}
