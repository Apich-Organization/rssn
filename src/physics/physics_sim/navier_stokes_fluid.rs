use crate::output::io::write_npy_file;
use crate::physics::physics_mtm::solve_poisson_2d_multigrid;
use ndarray::Array2;
use rayon::prelude::*;
use serde::{
    Deserialize,
    Serialize,
};

/// Parameters for the Navier-Stokes simulation.
#[derive(
    Clone, Debug, Serialize, Deserialize,
)]

pub struct NavierStokesParameters {
    pub nx: usize,
    pub ny: usize,
    pub re: f64,
    pub dt: f64,
    pub n_iter: usize,
    pub lid_velocity: f64,
}

/// Type of NavierStokesOutput.

pub type NavierStokesOutput = Result<
    (
        Array2<f64>,
        Array2<f64>,
        Array2<f64>,
    ),
    String,
>;

/// Main solver for the 2D lid-driven cavity problem.

pub fn run_lid_driven_cavity(
    params: &NavierStokesParameters
) -> NavierStokesOutput {

    let (nx, ny, _re, dt) = (
        params.nx, params.ny,
        params.re, params.dt,
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
    for j in 0..=nx {

        u[[ny - 1, j]] =
            params.lid_velocity;
    }

    let mg_size_k =
        ((nx.max(ny) - 1) as f64)
            .log2()
            .ceil() as u32;

    let mg_size =
        2_usize.pow(mg_size_k) + 1;

    for _ in 0..params.n_iter {

        let u_old = u.clone();

        let v_old = v.clone();

        // Calculate RHS in parallel
        let mut rhs_padded =
            vec![
                0.0;
                mg_size * mg_size
            ];

        let rhs_ptr = rhs_padded
            .as_mut_ptr()
            as usize;

        (1..ny - 1)
            .into_par_iter()
            .for_each(|j| {
                for i in 1..nx - 1 {

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
        (0..ny)
            .into_par_iter()
            .for_each(|j| {
                for i in 0..nx {

                    unsafe {

                        *(p_ptr as *mut f64).add(j * nx + i) += 0.7 * p_corr[[j, i]];
                    }
                }
            });

        // Update U
        (1..ny - 1)
            .into_par_iter()
            .for_each(|j| {
                for i in 1..nx {

                    unsafe {

                        *(u_ptr as *mut f64).add(j * (nx + 1) + i) -=
                            dt / hx * (p_corr[[j, i]] - p_corr[[j, i - 1]]);
                    }
                }
            });

        // Update V
        (1..ny)
            .into_par_iter()
            .for_each(|j| {
                for i in 1..nx - 1 {

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

    (0..ny)
        .into_par_iter()
        .for_each(|j| {
            for i in 0..nx {

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
        u_centered, v_centered, p,
    ))
}

/// An example scenario for the lid-driven cavity simulation.

pub fn simulate_lid_driven_cavity_scenario(
) {

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
                     results: {}",
                    e
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
                 during simulation: {}",
                e
            );
        },
    }
}
