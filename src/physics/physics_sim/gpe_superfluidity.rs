use ndarray::Array2;
use num_complex::Complex;
use rayon::prelude::*;
use serde::Deserialize;
use serde::Serialize;

use crate::output::io::write_npy_file;
use crate::physics::physics_sm::create_k_grid;
use crate::physics::physics_sm::fft2d;
use crate::physics::physics_sm::ifft2d;

/// Parameters for the GPE simulation.
#[derive(
    Clone, Debug, Serialize, Deserialize,
)]

pub struct GpeParameters {
    /// The number of grid points in the x-direction.
    pub nx: usize,
    /// The number of grid points in the y-direction.
    pub ny: usize,
    /// The length of the simulation domain in the x-direction.
    pub lx: f64,
    /// The length of the simulation domain in the y-direction.
    pub ly: f64,
    /// The imaginary time step.
    pub d_tau: f64,
    /// The number of time steps to simulate.
    pub time_steps: usize,
    /// The interaction strength.
    pub g: f64,
    /// The strength of the trapping potential.
    pub trap_strength: f64,
}

/// Finds the ground state of a 2D BEC by evolving the GPE in imaginary time.
///
/// # Returns
/// The final wave function `psi` representing the ground state.

pub fn run_gpe_ground_state_finder(
    params: &GpeParameters
) -> Result<Array2<f64>, String> {

    let dx =
        params.lx / params.nx as f64;

    let dy =
        params.ly / params.ny as f64;

    let n_total =
        (params.nx * params.ny) as f64;

    let mut potential = vec![
            0.0;
            params.nx * params.ny
        ];

    let nx = params.nx;

    let ny = params.ny;

    let trap_strength =
        params.trap_strength;

    potential
        .par_chunks_mut(nx)
        .enumerate()
        .for_each(|(j, row)| {

            let y = (j as f64
                - ny as f64 / 2.0)
                * dy;

            let y_sq = y.powi(2);

            for (i, val) in row
                .iter_mut()
                .enumerate()
            {

                let x = (i as f64
                    - nx as f64 / 2.0)
                    * dx;

                *val = 0.5
                    * trap_strength
                    * (x.powi(2)
                        + y_sq);
            }
        });

    let kx =
        create_k_grid(params.nx, dx);

    let ky =
        create_k_grid(params.ny, dy);

    let mut kinetic_operator = vec![
            0.0;
            params.nx * params.ny
        ];

    let d_tau = params.d_tau;

    kinetic_operator
        .par_chunks_mut(nx)
        .enumerate()
        .for_each(|(j, row)| {

            let ky_sq = ky[j].powi(2);

            for (i, val) in row
                .iter_mut()
                .enumerate()
            {

                let k_sq = kx[i]
                    .powi(2)
                    + ky_sq;

                *val = (-0.5
                    * k_sq
                    * d_tau)
                    .exp();
            }
        });

    let mut psi: Vec<Complex<f64>> =
        potential
            .iter()
            .map(|&v| {

                Complex::new(
                    (-v * 0.1).exp(),
                    0.0,
                )
            })
            .collect();

    for _ in 0 .. params.time_steps {

        psi.par_iter_mut()
            .enumerate()
            .for_each(|(idx, p)| {

                let v_eff = potential
                    [idx]
                    + params.g
                        * p.norm_sqr();

                *p *= (-v_eff
                    * params.d_tau
                    / 2.0)
                    .exp();
            });

        fft2d(
            &mut psi,
            params.nx,
            params.ny,
        );

        psi.par_iter_mut()
            .zip(&kinetic_operator)
            .for_each(|(p, k_op)| {

                *p *= k_op
            });

        ifft2d(
            &mut psi,
            params.nx,
            params.ny,
        );

        psi.par_iter_mut()
            .enumerate()
            .for_each(|(idx, p)| {

                let v_eff = potential
                    [idx]
                    + params.g
                        * p.norm_sqr();

                *p *= (-v_eff
                    * params.d_tau
                    / 2.0)
                    .exp();
            });

        let norm: f64 = psi
            .par_iter()
            .map(|p| p.norm_sqr())
            .sum();

        let norm_factor = (n_total
            / (norm * dx * dy))
            .sqrt();

        psi.par_iter_mut()
            .for_each(|p| {

                *p *= norm_factor
            });
    }

    let probability_density: Vec<f64> =
        psi.iter()
            .map(|p| p.norm_sqr())
            .collect();

    Array2::from_shape_vec(
        (params.ny, params.nx),
        probability_density,
    )
    .map_err(|e| e.to_string())
}

/// An example scenario that finds the ground state of a BEC, which may contain a vortex.

pub fn simulate_bose_einstein_vortex_scenario(
) -> Result<(), String> {

    println!(
        "Running GPE simulation to \
         find BEC ground state..."
    );

    let params = GpeParameters {
        nx: 128,
        ny: 128,
        lx: 20.0,
        ly: 20.0,
        d_tau: 0.01,
        time_steps: 500,
        g: 500.0,
        trap_strength: 1.0,
    };

    let final_density =
        run_gpe_ground_state_finder(
            &params,
        )?;

    let filename =
        "gpe_vortex_state.npy";

    println!(
        "Simulation finished. Saving \
         final density to {}",
        filename
    );

    write_npy_file(
        filename,
        &final_density,
    )?;

    Ok(())
}
