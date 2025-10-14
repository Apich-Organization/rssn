use std::sync::Arc;

use crate::numerical::elementary::eval_expr;
use crate::numerical::matrix::Matrix;
use crate::numerical::ode::solve_ode_system_rk4;
use crate::symbolic::core::Expr;
use rand::{thread_rng, Rng};
use std::collections::HashMap;

// =====================================================================================
// region: Classical Mechanics & Electromagnetism
// =====================================================================================

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
pub fn simulate_particle_motion(
    force_exprs: (&Expr, &Expr, &Expr),
    mass: f64,
    initial_pos: (f64, f64, f64),
    initial_vel: (f64, f64, f64),
    t_range: (f64, f64),
    num_steps: usize,
) -> Result<Vec<Vec<f64>>, String> {
    let (fx_expr, fy_expr, fz_expr) = force_exprs;
    let m_expr = Expr::Constant(mass);

    let ode_funcs: Vec<Expr> = vec![
        Expr::Variable("y3".to_string()), // dx/dt = vx
        Expr::Variable("y4".to_string()), // dy/dt = vy
        Expr::Variable("y5".to_string()), // dz/dt = vz
        Expr::Div(Arc::new(fx_expr.clone()), Arc::new(m_expr.clone())),
        Expr::Div(Arc::new(fy_expr.clone()), Arc::new(m_expr.clone())),
        Expr::Div(Arc::new(fz_expr.clone()), Arc::new(m_expr.clone())),
    ];

    let y0 = vec![
        initial_pos.0,
        initial_pos.1,
        initial_pos.2,
        initial_vel.0,
        initial_vel.1,
        initial_vel.2,
    ];

    solve_ode_system_rk4(&ode_funcs, &y0, t_range, num_steps)
}

// =====================================================================================
// region: Statistical Mechanics
// =====================================================================================

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
pub fn simulate_ising_model(size: usize, temperature: f64, steps: usize) -> Vec<Vec<i8>> {
    let mut rng = thread_rng();
    let mut lattice = vec![vec![0i8; size]; size];
    for i in 0..steps {
        for j in 0..size {
            lattice[i][j] = if rng.gen::<bool>() { 1 } else { -1 };
        }
    }

    for _ in 0..steps {
        let i = rng.gen_range(0..size);
        let j = rng.gen_range(0..size);

        // Sum of neighbor spins
        let top = lattice[(i + size - 1) % size][j];
        let bottom = lattice[(i + 1) % size][j];
        let left = lattice[i][(j + size - 1) % size];
        let right = lattice[i][(j + 1) % size];
        let neighbor_sum = top + bottom + left + right;

        // Change in energy if spin is flipped
        let delta_e = 2.0 * f64::from(lattice[i][j]) * f64::from(neighbor_sum);

        // Metropolis-Hastings acceptance criterion
        if delta_e < 0.0 || rng.gen::<f64>() < (-delta_e / temperature).exp() {
            lattice[i][j] *= -1; // Flip the spin
        }
    }

    lattice
}

// =====================================================================================
// region: Partial Differential Equations (Implementations)
// =====================================================================================

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
pub fn solve_1d_schrodinger(
    potential_expr: &Expr,
    var: &str,
    range: (f64, f64),
    num_points: usize,
) -> Result<(Vec<f64>, Matrix<f64>), String> {
    let (x_min, x_max) = range;
    if num_points < 3 {
        return Err("num_points must be >= 3".to_string());
    }
    let dx = (x_max - x_min) / (num_points as f64 - 1.0);
    let points: Vec<f64> = (0..num_points).map(|i| x_min + i as f64 * dx).collect();

    // Evaluate potential
    let mut potential_values = Vec::with_capacity(num_points);
    for &x in &points {
        let mut vars = HashMap::new();
        vars.insert(var.to_string(), x);
        potential_values.push(eval_expr(potential_expr, &vars)?);
    }

    // Build Hamiltonian dense matrix
    let n = num_points;
    let mut h_data = vec![0.0; n * n];
    let factor = -0.5 / (dx * dx);

    for i in 0..n {
        h_data[i * n + i] = -2.0 * factor + potential_values[i];
        if i > 0 {
            h_data[i * n + i - 1] = factor;
        }
        if i + 1 < n {
            h_data[i * n + i + 1] = factor;
        }
    }

    let hamiltonian = Matrix::new(n, n, h_data);
    hamiltonian.jacobi_eigen_decomposition(2000, 1e-12)
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
pub fn solve_2d_schrodinger(
    potential_expr: &Expr,
    var_x: &str,
    var_y: &str,
    ranges: (f64, f64, f64, f64),
    grid: (usize, usize),
) -> Result<(Vec<f64>, Matrix<f64>), String> {
    let (x_min, x_max, y_min, y_max) = ranges;
    let (nx, ny) = grid;
    if nx < 3 || ny < 3 {
        return Err("grid dimensions must be at least 3 in each direction".to_string());
    }

    let dx = (x_max - x_min) / (nx as f64 - 1.0);
    let dy = (y_max - y_min) / (ny as f64 - 1.0);

    // Precompute potential on grid (row-major: i -> x index, j -> y index)
    let n = nx * ny;
    let mut potential = vec![0.0_f64; n];
    for ix in 0..nx {
        for iy in 0..ny {
            let x = x_min + ix as f64 * dx;
            let y = y_min + iy as f64 * dy;
            let mut vars = HashMap::new();
            vars.insert(var_x.to_string(), x);
            vars.insert(var_y.to_string(), y);
            potential[ix * ny + iy] = eval_expr(potential_expr, &vars)?;
        }
    }

    // Build Hamiltonian dense matrix H (n x n)
    let mut h_data = vec![0.0_f64; n * n];
    // finite-difference factors
    let fx = -0.5 / (dx * dx);
    let fy = -0.5 / (dy * dy);

    // large potential to enforce Dirichlet boundary ~ psi=0 at boundaries
    let large = 1e12;

    for ix in 0..nx {
        for iy in 0..ny {
            let idx = ix * ny + iy;
            // If boundary point, force Dirichlet approximately
            if ix == 0 || ix == nx - 1 || iy == 0 || iy == ny - 1 {
                h_data[idx * n + idx] = large + potential[idx];
                continue;
            }

            // diagonal contribution: -2*(fx + fy) + V = (-2*fx -2*fy) + V
            // recall fx,fy are negative numbers (fx = -0.5/dx^2), so -2*fx = 1/dx^2 etc.
            h_data[idx * n + idx] = -2.0 * (fx + fy) + potential[idx];

            // neighbors in x direction
            let idx_left = (ix - 1) * ny + iy;
            let idx_right = (ix + 1) * ny + iy;
            h_data[idx * n + idx_left] = fx;
            h_data[idx * n + idx_right] = fx;

            // neighbors in y direction
            let idx_down = ix * ny + (iy - 1);
            let idx_up = ix * ny + (iy + 1);
            h_data[idx * n + idx_down] = fy;
            h_data[idx * n + idx_up] = fy;
        }
    }

    let hamiltonian = Matrix::new(n, n, h_data);
    // eigen decomposition
    hamiltonian.jacobi_eigen_decomposition(5000, 1e-10)
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
pub fn solve_3d_schrodinger(
    potential_expr: &Expr,
    var_x: &str,
    var_y: &str,
    var_z: &str,
    ranges: (f64, f64, f64, f64, f64, f64),
    grid: (usize, usize, usize),
) -> Result<(Vec<f64>, Matrix<f64>), String> {
    let (x_min, x_max, y_min, y_max, z_min, z_max) = ranges;
    let (nx, ny, nz) = grid;
    if nx < 3 || ny < 3 || nz < 3 {
        return Err("grid dimensions must be at least 3 in each direction".to_string());
    }

    let dx = (x_max - x_min) / (nx as f64 - 1.0);
    let dy = (y_max - y_min) / (ny as f64 - 1.0);
    let dz = (z_max - z_min) / (nz as f64 - 1.0);

    let n = nx * ny * nz;

    if n > 25000 {
        return Err(format!(
            "Grid too large (nx*ny*nz = {}). Dense 3D solver will be extremely slow and memory-heavy. Consider smaller grid or sparse/iterative methods.",
            n
        ));
    }

    let mut potential = vec![0.0_f64; n];
    for ix in 0..nx {
        for iy in 0..ny {
            for iz in 0..nz {
                let x = x_min + ix as f64 * dx;
                let y = y_min + iy as f64 * dy;
                let z = z_min + iz as f64 * dz;
                let mut vars = HashMap::new();
                vars.insert(var_x.to_string(), x);
                vars.insert(var_y.to_string(), y);
                vars.insert(var_z.to_string(), z);
                potential[(ix * ny + iy) * nz + iz] = eval_expr(potential_expr, &vars)?;
            }
        }
    }

    let mut h_data = vec![0.0_f64; n * n];
    let fx = -0.5 / (dx * dx);
    let fy = -0.5 / (dy * dy);
    let fz = -0.5 / (dz * dz);
    let large = 1e12;

    for ix in 0..nx {
        for iy in 0..ny {
            for iz in 0..nz {
                let idx = (ix * ny + iy) * nz + iz;
                // boundary
                if ix == 0 || ix == nx - 1 || iy == 0 || iy == ny - 1 || iz == 0 || iz == nz - 1 {
                    h_data[idx * n + idx] = large + potential[idx];
                    continue;
                }

                h_data[idx * n + idx] = -2.0 * (fx + fy + fz) + potential[idx];

                // six neighbors
                let idx_xm = ((ix - 1) * ny + iy) * nz + iz;
                let idx_xp = ((ix + 1) * ny + iy) * nz + iz;
                let idx_ym = (ix * ny + (iy - 1)) * nz + iz;
                let idx_yp = (ix * ny + (iy + 1)) * nz + iz;
                let idx_zm = (ix * ny + iy) * nz + (iz - 1);
                let idx_zp = (ix * ny + iy) * nz + (iz + 1);

                h_data[idx * n + idx_xm] = fx;
                h_data[idx * n + idx_xp] = fx;
                h_data[idx * n + idx_ym] = fy;
                h_data[idx * n + idx_yp] = fy;
                h_data[idx * n + idx_zm] = fz;
                h_data[idx * n + idx_zp] = fz;
            }
        }
    }

    let hamiltonian = Matrix::new(n, n, h_data);
    hamiltonian.jacobi_eigen_decomposition(5000, 1e-10)
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
        return Err("nx must be >=3 and nt >= 1".to_string());
    }
    let dx = (b - a) / (nx as f64 - 1.0);
    let dt = (t_range.1 - t_range.0) / (nt as f64);
    let r = alpha * dt / (dx * dx);

    // Build tridiagonal system matrices A (left) and B (right) for Crank-Nicolson:
    // A u^{n+1} = B u^{n} ; A = I + r/2 * L ; B = I - r/2 * L ; L is 2nd-difference operator
    // We'll solve with simple Gauss elimination for small nx.
    let interior = nx - 2; // number of unknowns (excluding Dirichlet boundary)
                           // Prebuild A (tridiagonal)
    let a_diag = vec![1.0 + r; interior];
    let a_lower = vec![-r / 2.0; interior - 1];
    let a_upper = vec![-r / 2.0; interior - 1];

    // B applied to u^n (interior only)
    let b_diag = vec![1.0 - r; interior];
    let b_lower = vec![r / 2.0; interior - 1];
    let b_upper = vec![r / 2.0; interior - 1];

    // initial condition
    let mut u0 = vec![0.0_f64; nx];
    for (i, var) in u0.iter_mut().enumerate().take(nx) {
        let x = a + i as f64 * dx;
        *var = init_func(x);
    }
    // apply Dirichlet boundaries u[0]=u[nx-1]=0
    u0[0] = 0.0;
    u0[nx - 1] = 0.0;

    let mut results: Vec<Vec<f64>> = Vec::with_capacity(nt + 1);
    results.push(u0.clone());

    // Helper: solve tridiagonal A x = d (Thomas algorithm)
    let solve_tridiag =
        |a_l: Vec<f64>, mut a_d: Vec<f64>, a_u: Vec<f64>, mut d: Vec<f64>| -> Vec<f64> {
            let n = a_d.len();
            // forward
            for i in 1..n {
                let m = a_l[i - 1] / a_d[i - 1];
                a_d[i] -= m * a_u[i - 1];
                d[i] -= m * d[i - 1];
            }
            let mut x = vec![0.0_f64; n];
            x[n - 1] = d[n - 1] / a_d[n - 1];
            for i in (0..n - 1).rev() {
                x[i] = (d[i] - a_u[i] * x[i + 1]) / a_d[i];
            }
            x
        };

    // time stepping
    for _step in 0..nt {
        // compute right-hand side d = B * u^n (interior)
        let mut d = vec![0.0_f64; interior];
        for i in 0..interior {
            // global index = i+1
            let global_i = i + 1;
            let mut val = b_diag[i] * u0[global_i];
            if i > 0 {
                val += b_lower[i - 1] * u0[global_i - 1];
            }
            if i + 1 < interior {
                val += b_upper[i] * u0[global_i + 1];
            }
            // boundaries are zero so no extra terms
            d[i] = val;
        }

        // solve A u^{n+1}_{interior} = d
        let u_interior = solve_tridiag(a_lower.clone(), a_diag.clone(), a_upper.clone(), d);

        // assemble full solution with Dirichlet boundaries
        let mut u_new = vec![0.0_f64; nx];
        u_new[0] = 0.0;
        /*for i in 0..interior {
            u_new[i + 1] = u_interior[i];
        }*/
		u_new[1..=interior].copy_from_slice(&u_interior[..interior]);
        u_new[nx - 1] = 0.0;

        results.push(u_new.clone());
        u0 = u_new;
    }

    Ok(results)
}

// ==========================
// Extended PDE / QM Solvers
// ==========================

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
        return Err("initial_ut length mismatch".to_string());
    }
    if n < 3 {
        return Err("need at least 3 grid points".to_string());
    }
    // CFL check (simple)
    let cfl = c * dt / dx;
    if cfl.abs() > 1.0 {
        return Err(format!(
            "CFL violation: c*dt/dx = {} > 1. Reduce dt or increase dx.",
            cfl
        ));
    }

    // u^{n-1}, u^{n}
    let mut u_prev = initial_u.to_owned();
    // compute u at first time step using Taylor expansion:
    // u^{1} = u^{0} + dt * ut^{0} + 0.5 * (c dt)^2 * u_xx^{0}
    let mut u_curr = vec![0.0; n];
    // interior second derivative
    for i in 1..(n - 1) {
        let u_xx = (u_prev[i - 1] - 2.0 * u_prev[i] + u_prev[i + 1]) / (dx * dx);
        u_curr[i] = u_prev[i] + dt * initial_ut[i] + 0.5 * (c * c) * (dt * dt) * u_xx;
    }
    // Dirichlet BCs: keep boundary zero (or same as initial)
    u_curr[0] = 0.0;
    u_curr[n - 1] = 0.0;

    let mut snapshots: Vec<Vec<f64>> = Vec::with_capacity(num_steps + 1);
    snapshots.push(u_prev.clone()); // t=0
    snapshots.push(u_curr.clone()); // t=dt

    // Time stepping (leap-frog / central time difference)
    for _step in 2..=num_steps {
        let mut u_next = vec![0.0; n];
        for i in 1..(n - 1) {
            let u_xx = (u_curr[i - 1] - 2.0 * u_curr[i] + u_curr[i + 1]) / (dx * dx);
            u_next[i] = 2.0 * u_curr[i] - u_prev[i] + (c * c) * (dt * dt) * u_xx;
        }
        u_next[0] = 0.0;
        u_next[n - 1] = 0.0;

        snapshots.push(u_next.clone());
        u_prev = u_curr;
        u_curr = u_next;
    }

    Ok(snapshots)
}
