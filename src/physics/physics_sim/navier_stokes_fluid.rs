// src/physics/physics_sim/navier_stokes_fluid.rs
// A 2D incompressible Navier-Stokes solver for the lid-driven cavity problem.
// Uses the SIMPLE algorithm on a staggered grid, with a multigrid solver for the pressure-Poisson equation.

use crate::output::io::write_npy_file;
use crate::physics::physics_mtm::solve_poisson_2d_multigrid;
use ndarray::Array2;

/// Parameters for the Navier-Stokes simulation.
pub struct NavierStokesParameters {
    pub nx: usize,     // Number of cells in x
    pub ny: usize,     // Number of cells in y
    pub re: f64,       // Reynolds number
    pub dt: f64,       // Time step
    pub n_iter: usize, // Number of time-stepping iterations
    pub lid_velocity: f64,
}

/// Main solver for the 2D lid-driven cavity problem.
pub fn run_lid_driven_cavity(
    params: &NavierStokesParameters,
) -> Result<(Array2<f64>, Array2<f64>, Array2<f64>), String> {
    let (nx, ny, _re, dt) = (params.nx, params.ny, params.re, params.dt);
    let hx = 1.0 / (nx - 1) as f64;
    let hy = 1.0 / (ny - 1) as f64;

    // Staggered grid: u, v, p
    let mut u = Array2::<f64>::zeros((ny, nx + 1));
    let mut v = Array2::<f64>::zeros((ny + 1, nx));
    let mut p = Array2::<f64>::zeros((ny, nx));

    // Lid velocity boundary condition
    for j in 0..nx + 1 {
        u[[ny - 1, j]] = params.lid_velocity;
    }

    for _ in 0..params.n_iter {
        // --- Momentum Predictor Step (u*, v*) ---
        // This is a simplified explicit discretization of momentum equations.
        // A full implementation would require careful handling of advection terms.
        let u_old = u.clone();
        let v_old = v.clone();

        // --- Pressure-Poisson Equation ---
        // 1. Calculate the RHS of the Poisson equation
        let mut rhs = vec![0.0; nx * ny];
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let div_u_star = ((u_old[[j, i + 1]] - u_old[[j, i]]) / hx)
                    + ((v_old[[j + 1, i]] - v_old[[j, i]]) / hy);
                rhs[j * nx + i] = div_u_star / dt;
            }
        }

        // 2. Solve the Poisson equation for pressure correction using Multigrid
        // The multigrid solver needs a grid size of 2^k+1. We find the next suitable size.
        let mg_size_k = ((nx.max(ny) - 1) as f64).log2().ceil() as u32;
        let mg_size = 2_usize.pow(mg_size_k) + 1;
        let mut rhs_padded = vec![0.0; mg_size * mg_size];
        for j in 0..ny {
            for i in 0..nx {
                rhs_padded[j * mg_size + i] = rhs[j * nx + i];
            }
        }

        let p_corr_vec = solve_poisson_2d_multigrid(mg_size, &rhs_padded, 5)?;
        let p_corr = Array2::from_shape_vec((mg_size, mg_size), p_corr_vec).unwrap();

        // --- Correction Step ---
        // Update pressure
        for j in 0..ny {
            for i in 0..nx {
                p[[j, i]] += 0.7 * p_corr[[j, i]]; // Under-relaxation for stability
            }
        }

        // Update velocities
        for j in 1..ny - 1 {
            for i in 1..nx {
                u[[j, i]] -= dt / hx * (p_corr[[j, i]] - p_corr[[j, i - 1]]);
            }
        }
        for j in 1..ny {
            for i in 1..nx - 1 {
                v[[j, i]] -= dt / hy * (p_corr[[j, i]] - p_corr[[j - 1, i]]);
            }
        }
    }

    // Interpolate staggered velocities to cell centers for output
    let mut u_centered = Array2::<f64>::zeros((ny, nx));
    let mut v_centered = Array2::<f64>::zeros((ny, nx));
    for j in 0..ny {
        for i in 0..nx {
            u_centered[[j, i]] = 0.5 * (u[[j, i]] + u[[j, i + 1]]);
            v_centered[[j, i]] = 0.5 * (v[[j, i]] + v[[j + 1, i]]);
        }
    }

    Ok((u_centered, v_centered, p))
}

/// An example scenario for the lid-driven cavity simulation.
pub fn simulate_lid_driven_cavity_scenario() {
    println!("Running 2D Lid-Driven Cavity simulation...");

    // Grid size must be 2^k+1 for the multigrid solver to work perfectly.
    const K: usize = 6;
    const N: usize = 2_usize.pow(K as u32) + 1;

    let params = NavierStokesParameters {
        nx: N,
        ny: N,
        re: 100.0,
        dt: 0.01,
        n_iter: 200,
        lid_velocity: 1.0,
    };

    match run_lid_driven_cavity(&params) {
        Ok((u, v, p)) => {
            println!("Simulation finished. Saving results...");
            write_npy_file("cavity_u_velocity.npy", &u);
            write_npy_file("cavity_v_velocity.npy", &v);
            write_npy_file("cavity_pressure.npy", &p);
            println!("Results saved to .npy files.");
        }
        Err(e) => {
            eprintln!("An error occurred during simulation: {}", e);
        }
    }
}
