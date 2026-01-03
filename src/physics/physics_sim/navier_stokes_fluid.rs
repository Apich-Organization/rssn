use ndarray::Array2;
use rayon::prelude::*;
use serde::Deserialize;
use serde::Serialize;

use crate::output::io::write_npy_file;
use crate::physics::physics_mtm::solve_poisson_2d_multigrid;

/// Parameters for the Navier-Stokes simulation.
#[derive(
    Clone, Debug, Serialize, Deserialize,
)]

pub struct NavierStokesParameters {
    /// Number of grid points in the x-direction.
    pub nx: usize,
    /// Number of grid points in the y-direction.
    pub ny: usize,
    /// Reynolds number.
    pub re: f64,
    /// Time step size.
    pub dt: f64,
    /// Number of simulation iterations.
    pub n_iter: usize,
    /// Velocity of the lid (driving force).
    pub lid_velocity: f64,
}

/// Type of `NavierStokesOutput`.

pub type NavierStokesOutput = Result<
    (
        Array2<f64>,
        Array2<f64>,
        Array2<f64>,
    ),
    String,
>;

/// Solves the 2D Navier-Stokes equations for channel flow with an obstacle.
/// Uses the projection method (Chorin 1968).
///
/// # Arguments
/// * `nx` - Grid width.
/// * `ny` - Grid height.
/// * `re` - Reynolds number.
/// * `dt` - Time step.
/// * `n_iter` - Number of iterations.
/// * `obstacle_mask` - A boolean mask where true indicates an obstacle (u=0).
///
/// # Returns
/// Tuple of (u, v, p) arrays.

pub fn run_channel_flow(
    nx: usize,
    _ny: usize,
    re: f64,
    dt: f64,
    n_iter: usize,
    obstacle_mask: &Array2<bool>,
) -> NavierStokesOutput {

    let n = nx;

    let h = 1.0 / (n as f64 - 1.0);

    let nu = 1.0 / re;

    let mut u =
        Array2::<f64>::zeros((n, n));

    let mut v =
        Array2::<f64>::zeros((n, n));

    let mut p =
        Array2::<f64>::zeros((n, n));

    // Helper for indexing
    let idx =
        |i: usize, j: usize| j * n + i;

    // Multigrid size
    let mg_size = n;

    for _iter in 0 .. n_iter {

        let u_old = u.clone();

        let v_old = v.clone();

        // 1. Advection-Diffusion (Explicit) -> Intermediate Velocity (u*, v*)
        // u_star = u_old + dt * ( - (u grad) u + nu laplacian u )
        let mut u_star = u.clone();

        let mut v_star = v.clone();

        u_star
            .as_slice_mut()
            .expect("Contiguous array")
            .par_iter_mut()
            .zip(v_star.as_slice_mut().expect("Contiguous array").par_iter_mut())
            .enumerate()
            .for_each(|(id, (u_val, v_val))| {
                let i = id % n;
                let j = id / n;

                // Skip boundaries and obstacles for now
                if i == 0 || i == n - 1 || j == 0 || j == n - 1 || obstacle_mask[[j, i]] {
                    return;
                }

                // Upwind Advection
                let u_curr = u_old[[j, i]];
                let v_curr = v_old[[j, i]];

                let du_dx = if u_curr > 0.0 {
                    u_curr - u_old[[j, i - 1]]
                } else {
                    u_old[[j, i + 1]] - u_curr
                } / h;

                let du_dy = if v_curr > 0.0 {
                    u_curr - u_old[[j - 1, i]]
                } else {
                    u_old[[j + 1, i]] - u_curr
                } / h;

                let dv_dx = if u_curr > 0.0 {
                    v_curr - v_old[[j, i - 1]]
                } else {
                    v_old[[j, i + 1]] - v_curr
                } / h;

                let dv_dy = if v_curr > 0.0 {
                    v_curr - v_old[[j - 1, i]]
                } else {
                    v_old[[j + 1, i]] - v_curr
                } / h;

                // Diffusion (Central Difference)
                let lap_u = (u_old[[j, i + 1]] - 2.0 * u_curr + u_old[[j, i - 1]]
                    + u_old[[j + 1, i]] - 2.0 * u_curr + u_old[[j - 1, i]]) / (h * h);

                let lap_v = (v_old[[j, i + 1]] - 2.0 * v_curr + v_old[[j, i - 1]]
                    + v_old[[j + 1, i]] - 2.0 * v_curr + v_old[[j - 1, i]]) / (h * h);

                *u_val = u_curr + dt * (-(u_curr * du_dx + v_curr * du_dy) + nu * lap_u);
                *v_val = v_curr + dt * (-(u_curr * dv_dx + v_curr * dv_dy) + nu * lap_v);
            });

        // Apply BCs to u_star, v_star
        // Inflow (Left)
        for j in 0 .. n {

            u_star[[j, 0]] = 1.0;

            v_star[[j, 0]] = 0.0;
        }

        // Outflow (Right) - Zero Gradient
        for j in 0 .. n {

            u_star[[j, n - 1]] =
                u_star[[j, n - 2]];

            v_star[[j, n - 1]] =
                v_star[[j, n - 2]];
        }

        // Walls (Top/Bottom) - We do Slip for channel here
        for i in 0 .. n {

            u_star[[0, i]] =
                u_star[[1, i]]; // Slip
            v_star[[0, i]] = 0.0;

            u_star[[n - 1, i]] =
                u_star[[n - 2, i]]; // Slip
            v_star[[n - 1, i]] = 0.0;
        }

        // Obstacle - No Slip
        for j in 0 .. n {

            for i in 0 .. n {

                if obstacle_mask[[j, i]]
                {

                    u_star[[j, i]] =
                        0.0;

                    v_star[[j, i]] =
                        0.0;
                }
            }
        }

        // 2. Pressure Correction (Poisson Step)
        // laplacian p = div(u*) / dt
        let mut rhs_vec =
            vec![0.0; n * n];

        // Calculate Div U*
        for j in 1 .. n - 1 {

            for i in 1 .. n - 1 {

                if obstacle_mask[[j, i]]
                {

                    continue;
                }

                let div = (u_star
                    [[j, i + 1]]
                    - u_star
                        [[j, i - 1]])
                    / (2.0 * h)
                    + (v_star
                        [[j + 1, i]]
                        - v_star[[
                            j - 1,
                            i,
                        ]])
                        / (2.0 * h);

                rhs_vec[idx(i, j)] =
                    div / dt;
            }
        }

        // Solve Poisson
        // Note: Our MG solver expects "f" where "-laplacian u = f".
        // Our equation is "laplacian p = R". So we pass "-R" as f.
        for val in rhs_vec.iter_mut() {

            *val = -*val;
        }

        let p_flat =
            solve_poisson_2d_multigrid(
                mg_size,
                &rhs_vec,
                2,
            )?; // 2 V-cycles

        // Copy back to P array
        for j in 0 .. n {

            for i in 0 .. n {

                p[[j, i]] =
                    p_flat[idx(i, j)];
            }
        }

        // 3. Velocity Correction
        // u_new = u* - dt * grad p
        for j in 1 .. n - 1 {

            for i in 1 .. n - 1 {

                if obstacle_mask[[j, i]]
                {

                    u[[j, i]] = 0.0;

                    v[[j, i]] = 0.0;

                    continue;
                }

                let dp_dx = (p
                    [[j, i + 1]]
                    - p[[j, i - 1]])
                    / (2.0 * h);

                let dp_dy = (p
                    [[j + 1, i]]
                    - p[[j - 1, i]])
                    / (2.0 * h);

                u[[j, i]] = u_star
                    [[j, i]]
                    - dt * dp_dx;

                v[[j, i]] = v_star
                    [[j, i]]
                    - dt * dp_dy;
            }
        }

        // Re-apply BCs for u, v
        // Inflow (Left)
        for j in 0 .. n {

            u[[j, 0]] = 1.0;

            v[[j, 0]] = 0.0;
        }

        // Outflow (Right)
        for j in 0 .. n {

            u[[j, n - 1]] =
                u[[j, n - 2]];

            v[[j, n - 1]] =
                v[[j, n - 2]];
        }

        // Walls (Top/Bottom)
        for i in 0 .. n {

            u[[0, i]] = u[[1, i]];

            v[[0, i]] = 0.0;

            u[[n - 1, i]] =
                u[[n - 2, i]];

            v[[n - 1, i]] = 0.0;
        }
    }

    Ok((u, v, p))
}

/// Main solver for the 2D lid-driven cavity problem.
///
/// # Errors
///
/// This function will return an error if the underlying Poisson solver fails,
/// or if there are issues reshaping the pressure correction array.

pub fn run_lid_driven_cavity(
    params: &NavierStokesParameters
) -> NavierStokesOutput {

    let (nx, ny, _re, dt) = (
        params.nx,
        params.ny,
        params.re,
        params.dt,
    );

    let hx = 1.0 / (nx - 1) as f64;

    let hy = 1.0 / (ny - 1) as f64;

    let mut u = Array2::<f64>::zeros((
        ny,
        nx + 1,
    ));

    let mut v = Array2::<f64>::zeros((
        ny + 1,
        nx,
    ));

    let mut p =
        Array2::<f64>::zeros((ny, nx));

    // Boundary conditions: lid velocity
    for j in 0 ..= nx {

        u[[ny - 1, j]] =
            params.lid_velocity;
    }

    let mg_size_k =
        (((nx.max(ny) - 1) as f64)
            .log2()
            .ceil() as i64)
            .try_into()
            .unwrap_or(0);

    let mg_size =
        2_usize.pow(mg_size_k) + 1;

    for _ in 0 .. params.n_iter {

        let u_old = u.clone();

        let v_old = v.clone();

        // Calculate RHS in parallel
        let mut rhs_padded = vec![
                0.0;
                mg_size * mg_size
            ];

        let rhs_ptr = rhs_padded
            .as_mut_ptr()
            as usize;

        (1 .. ny - 1)
            .into_par_iter()
            .for_each(|j| {
                for i in 1 .. nx - 1 {

                    let div_u_star = ((u_old[[j, i + 1]] - u_old[[j, i]]) / hx)
                        + ((v_old[[j + 1, i]] - v_old[[j, i]]) / hy);

                    unsafe {

                        *(rhs_ptr as *mut f64).add(j * mg_size + i) = div_u_star / dt;
                    }
                }
            });

        // Solve Poisson for pressure correction
        let p_corr_vec =
            solve_poisson_2d_multigrid(
                mg_size,
                &rhs_padded,
                10,
            )?; // More V-cycles for accuracy
        let p_corr =
            Array2::from_shape_vec(
                (mg_size, mg_size),
                p_corr_vec,
            )
            .map_err(
                |e| e.to_string(),
            )?;

        // Update pressure and velocities in parallel
        let p_ptr =
            p.as_mut_ptr() as usize;

        let u_ptr =
            u.as_mut_ptr() as usize;

        let v_ptr =
            v.as_mut_ptr() as usize;

        // Update P
        (0 .. ny)
            .into_par_iter()
            .for_each(|j| {
                for i in 0 .. nx {

                    unsafe {

                        *(p_ptr as *mut f64).add(j * nx + i) += 0.7 * p_corr[[j, i]];
                    }
                }
            });

        // Update U
        (1 .. ny - 1)
            .into_par_iter()
            .for_each(|j| {
                for i in 1 .. nx {

                    unsafe {

                        *(u_ptr as *mut f64).add(j * (nx + 1) + i) -=
                            dt / hx * (p_corr[[j, i]] - p_corr[[j, i - 1]]);
                    }
                }
            });

        // Update V
        (1 .. ny)
            .into_par_iter()
            .for_each(|j| {
                for i in 1 .. nx - 1 {

                    unsafe {

                        *(v_ptr as *mut f64).add(j * nx + i) -=
                            dt / hy * (p_corr[[j, i]] - p_corr[[j - 1, i]]);
                    }
                }
            });
    }

    // ... centering ...
    let mut u_centered =
        Array2::<f64>::zeros((ny, nx));

    let mut v_centered =
        Array2::<f64>::zeros((ny, nx));

    let uc_ptr = u_centered.as_mut_ptr()
        as usize;

    let vc_ptr = v_centered.as_mut_ptr()
        as usize;

    (0 .. ny)
        .into_par_iter()
        .for_each(|j| {
            for i in 0 .. nx {

                unsafe {

                    *(uc_ptr
                        as *mut f64)
                        .add(
                            j * nx + i,
                        ) = 0.5
                        * (u[[j, i]]
                            + u[[
                                j,
                                i + 1,
                            ]]);

                    *(vc_ptr
                        as *mut f64)
                        .add(
                            j * nx + i,
                        ) = 0.5
                        * (v[[j, i]]
                            + v[[
                                j + 1,
                                i,
                            ]]);
                }
            }
        });

    Ok((
        u_centered,
        v_centered,
        p,
    ))
}

/// An example scenario for the lid-driven cavity simulation.

pub fn simulate_lid_driven_cavity_scenario()
 {

    const K: usize = 6;

    const N: usize =
        2_usize.pow(K as u32) + 1;

    println!(
        "Running 2D Lid-Driven Cavity \
         simulation..."
    );

    let params =
        NavierStokesParameters {
            nx: N,
            ny: N,
            re: 100.0,
            dt: 0.01,
            n_iter: 200,
            lid_velocity: 1.0,
        };

    match run_lid_driven_cavity(&params)
    {
        | Ok((u, v, p)) => {

            println!(
                "Simulation finished. \
                 Saving results..."
            );

            let save_result = (|| -> Result<(), String> {

                write_npy_file(
                    "cavity_u_velocity.npy",
                    &u,
                )?;

                write_npy_file(
                    "cavity_v_velocity.npy",
                    &v,
                )?;

                write_npy_file(
                    "cavity_pressure.npy",
                    &p,
                )?;

                Ok(())
            })();

            if let Err(e) = save_result
            {

                eprintln!(
                    "Failed to save \
                     results: {e}"
                );
            } else {

                println!(
                    "Results saved to \
                     .npy files."
                );
            }
        },
        | Err(e) => {

            eprintln!(
                "An error occurred \
                 during simulation: \
                 {e}"
            );
        },
    }
}
