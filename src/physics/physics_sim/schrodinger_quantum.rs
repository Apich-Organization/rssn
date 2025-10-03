// src/physics/physics_sim/schrodinger_quantum.rs
// 2D Time-Dependent Schrodinger Equation solver using the Split-Step Fourier Method.

use crate::output::io::write_npy_file;
use crate::physics::physics_sm::{create_k_grid, fft2d, ifft2d};
use ndarray::Array2;
use num_complex::Complex;
use rayon::prelude::*;

/// Parameters for the Schrodinger simulation.
pub struct SchrodingerParameters {
    pub nx: usize, // Number of grid points in x
    pub ny: usize, // Number of grid points in y
    pub lx: f64,   // Domain length in x
    pub ly: f64,   // Domain length in y
    pub dt: f64,   // Time step
    pub time_steps: usize,
    pub hbar: f64,           // Reduced Planck's constant
    pub mass: f64,           // Particle mass
    pub potential: Vec<f64>, // Potential V(x, y) flattened
}

/// Runs a 2D Schrodinger simulation using the Split-Step Fourier method.
///
/// # Arguments
/// * `params` - The simulation parameters.
/// * `initial_psi` - The initial wave function (complex values, flattened).
///
/// # Returns
/// A `Vec` containing snapshots of the probability density `|psi|^2`.
pub fn run_schrodinger_simulation(
    params: &SchrodingerParameters,
    initial_psi: &mut Vec<Complex<f64>>,
) -> Vec<Array2<f64>> {
    let dx = params.lx / params.nx as f64;
    let dy = params.ly / params.ny as f64;

    // Create momentum-space grids
    let kx = create_k_grid(params.nx, dx);
    let ky = create_k_grid(params.ny, dy);

    // Pre-calculate the kinetic energy evolution operator
    let mut kinetic_operator = vec![Complex::default(); params.nx * params.ny];
    for j in 0..params.ny {
        for i in 0..params.nx {
            let k_sq = kx[i].powi(2) + ky[j].powi(2);
            let phase = -params.hbar * k_sq / (2.0 * params.mass) * params.dt;
            kinetic_operator[j * params.nx + i] = Complex::from_polar(1.0, phase);
        }
    }

    // Pre-calculate the potential energy evolution operator
    let potential_operator: Vec<_> = params
        .potential
        .par_iter()
        .map(|&v| {
            let phase = -v * params.dt / (2.0 * params.hbar);
            Complex::from_polar(1.0, phase)
        })
        .collect();

    let mut snapshots = Vec::new();
    let mut psi = initial_psi.clone();

    for t_step in 0..params.time_steps {
        // First half-step in potential space
        psi.par_iter_mut()
            .zip(&potential_operator)
            .for_each(|(p, v_op)| *p *= v_op);

        // Step in momentum space
        fft2d(&mut psi, params.nx, params.ny);
        psi.par_iter_mut()
            .zip(&kinetic_operator)
            .for_each(|(p, k_op)| *p *= k_op);
        ifft2d(&mut psi, params.nx, params.ny);

        // Second half-step in potential space
        psi.par_iter_mut()
            .zip(&potential_operator)
            .for_each(|(p, v_op)| *p *= v_op);

        // Save snapshot
        if t_step % 10 == 0 {
            let probability_density: Vec<f64> = psi.par_iter().map(|p| p.norm_sqr()).collect();
            snapshots
                .push(Array2::from_shape_vec((params.ny, params.nx), probability_density).unwrap());
        }
    }

    snapshots
}

/// An example scenario simulating a wave packet hitting a double slit.
pub fn simulate_double_slit_scenario() {
    println!("Running 2D Schrodinger simulation for a double slit...");

    const NX: usize = 256;
    const NY: usize = 256;

    // Setup potential V(x,y) for a double slit
    let mut potential = vec![0.0; NX * NY];
    let slit_center_y = NY / 2;
    let slit_width = 10;
    let slit_spacing = 40;
    let barrier_x = NX / 5;
    for j in 0..NY {
        let is_slit1 = j > slit_center_y - slit_spacing / 2 - slit_width
            && j < slit_center_y - slit_spacing / 2;
        let is_slit2 = j > slit_center_y + slit_spacing / 2
            && j < slit_center_y + slit_spacing / 2 + slit_width;
        if !is_slit1 && !is_slit2 {
            potential[j * NX + barrier_x] = 1e5; // High potential barrier
        }
    }

    let params = SchrodingerParameters {
        nx: NX,
        ny: NY,
        lx: NX as f64,
        ly: NY as f64,
        dt: 0.1,
        time_steps: 300,
        hbar: 1.0,
        mass: 1.0,
        potential,
    };

    // Setup initial wave function: a Gaussian packet moving towards the slits
    let mut initial_psi = vec![Complex::default(); NX * NY];
    let initial_pos = (NX as f64 / 10.0, NY as f64 / 2.0);
    let initial_momentum = (5.0, 0.0);
    let packet_width_sq = 100.0;
    for j in 0..NY {
        for i in 0..NX {
            let dx = i as f64 - initial_pos.0;
            let dy = j as f64 - initial_pos.1;
            let norm_sq = dx * dx + dy * dy;
            let phase = initial_momentum.0 * dx + initial_momentum.1 * dy;
            let envelope = (-norm_sq / (2.0 * packet_width_sq)).exp();
            initial_psi[j * NX + i] = Complex::from_polar(envelope, phase);
        }
    }

    let snapshots = run_schrodinger_simulation(&params, &mut initial_psi);

    if let Some(final_state) = snapshots.last() {
        let filename = "schrodinger_double_slit.npy";
        println!("Saving final probability density to {}", filename);
        write_npy_file(filename, final_state);
    } else {
        println!("Simulation produced no snapshots.");
    }
}
