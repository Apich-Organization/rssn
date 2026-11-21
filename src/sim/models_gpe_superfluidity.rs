//! # Gross–Pitaevskii ground states (2D)
//!
//! | Problem | Dimension | Method |
//! |---|---|---|
//! | imaginary-time ground state | 2D | split-step Fourier, periodic box |
//!
//! [`crate::sim::split_step`] covers 1D–3D, real-time evolution (Strang and
//! fourth-order composition) and the linear Schrödinger equation.

use ndarray::Array2;
use num_complex::Complex;
use rayon::prelude::*;
use serde::Deserialize;
use serde::Serialize;

use std::path::Path;

use crate::io::write_npy_file;
use crate::kernels::fft::Norm;
use crate::kernels::fft::transform_axes_in_place;
use crate::kernels::numeric_common::scale_product;
use crate::kernels::par::for_chunks_mut;
use crate::kernels::par::for_each_mut;
use crate::kernels::par::sum_fixed;
use crate::kernels::reuse::Scratch;
use crate::kernels::reuse::fft_scratch_for;
use crate::sim::physics_sm::create_k_grid;

const POINTWISE_GRAIN: usize = 16_384;

/// Parameters for the GPE simulation.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct GpeParameters {
    /// The number of grid points in the x-direction.
    pub nx: usize,
    /// The number of grid points in the y-direction.
    pub ny: usize,
    /// The length of the simulation domain in the x-direction.
    pub lx: f64,
    /// The length of the simulation domain in the y-direction.
    pub ly: f64,
    /// The nonnegative imaginary time increment in units with `hbar = m = 1`.
    pub d_tau: f64,
    /// The number of time steps to simulate.
    pub time_steps: usize,
    /// The contact coupling, with units of energy times area.
    pub g: f64,
    /// The harmonic trap coefficient, with units of energy per squared length.
    pub trap_strength: f64,
}

/// Evolves the 2D GPE for a fixed number of imaginary-time steps on a periodic box.
///
/// The ground-state probability density `|psi|^2` has rows along `y` and
/// columns along `x`, with physical normalization `sum(|psi|^2) dx dy = 1`.
/// Zero steps or a zero time increment return the normalized initial profile.
///
/// # Errors
/// Returns an error for an empty or oversized grid, unrepresentable spacing,
/// non-positive box lengths, negative imaginary time, non-finite coefficients,
/// or a propagated state whose normalized density cannot be represented.
pub fn run_gpe_ground_state_finder(params: &GpeParameters) -> Result<Array2<f64>, String> {
    let nx = params.nx;
    let ny = params.ny;
    if nx == 0 || ny == 0 {
        return Err("GPE grid extents must be positive".into());
    }
    let size = nx
        .checked_mul(ny)
        .filter(|&size| size <= isize::MAX as usize / size_of::<Complex<f64>>())
        .ok_or("GPE grid size overflows addressable storage")?;
    let dx = params.lx / nx as f64;
    let dy = params.ly / ny as f64;
    let d_tau = params.d_tau;
    if dx <= 0.0
        || dy <= 0.0
        || !dx.is_finite()
        || !dy.is_finite()
        || !d_tau.is_finite()
        || d_tau < 0.0
        || !params.g.is_finite()
        || !params.trap_strength.is_finite()
    {
        return Err("invalid GPE spacing, imaginary time or coefficient".into());
    }

    let mut potential = vec![0.0; size];
    if params.trap_strength != 0.0 {
        let strength = params.trap_strength;
        let half_trap = 0.5 * strength;
        potential
            .par_chunks_mut(nx)
            .enumerate()
            .try_for_each(|(j, row)| {
                let y = (j as f64 - ny as f64 / 2.0) * dy;
                let y_sq = y.powi(2);
                for (i, val) in row.iter_mut().enumerate() {
                    let x = (i as f64 - nx as f64 / 2.0) * dx;
                    let radius_sq = x.mul_add(x, y_sq);
                    *val = if radius_sq.is_normal() && half_trap.is_normal() {
                        half_trap * radius_sq
                    } else {
                        trapped_axis_energy(x, strength)? + trapped_axis_energy(y, strength)?
                    };
                    if !val.is_finite() {
                        return Err("GPE trapping potential is not representable");
                    }
                }
                Ok(())
            })?;
    }

    let propagate = params.time_steps != 0 && d_tau > 0.0;
    // The initial amplitude enters the nonlinear propagation. At zero time
    // only the normalized profile is needed, so its largest weight can be one.
    let initial_offset = if propagate {
        0.0
    } else {
        potential.iter().copied().fold(f64::INFINITY, f64::min)
    };
    let mut psi: Vec<_> = potential
        .iter()
        .map(|&v| Complex::new((-(v - initial_offset) * 0.1).exp(), 0.0))
        .collect();
    if propagate {
        let kx = create_k_grid(nx, dx);
        let ky = create_k_grid(ny, dy);
        let mut kinetic_operator = vec![0.0; size];
        kinetic_operator
            .par_chunks_mut(nx)
            .enumerate()
            .for_each(|(j, row)| {
                let ky_sq = ky[j].powi(2);
                for (i, val) in row.iter_mut().enumerate() {
                    let k_sq = kx[i].mul_add(kx[i], ky_sq);
                    let damping = if k_sq.is_normal() {
                        0.5 * k_sq * d_tau
                    } else {
                        gpe_axis_damping(i, nx, params.lx, d_tau)
                            + gpe_axis_damping(j, ny, params.ly, d_tau)
                    };
                    *val = (-damping).exp();
                }
            });

        let mut fft = GpeFft::new([ny, nx]);
        for _ in 0..params.time_steps {
            gpe_potential_step(&mut psi, &potential, params.g, d_tau);
            fft.transform(&mut psi, false);
            for_each_mut(&mut psi, POINTWISE_GRAIN, |i, p| *p *= kinetic_operator[i]);
            fft.transform(&mut psi, true);
            gpe_potential_step(&mut psi, &potential, params.g, d_tau);
            gpe_normalize(&mut psi, dx, dy)?;
        }
    } else {
        gpe_normalize(&mut psi, dx, dy)?;
    }

    let probability_density: Vec<f64> = psi.iter().map(Complex::norm_sqr).collect();
    if probability_density.iter().any(|p| !p.is_finite())
        || probability_density.iter().all(|&p| p == 0.0)
    {
        return Err("GPE normalized density is not representable".into());
    }
    Array2::from_shape_vec((ny, nx), probability_density).map_err(|e| e.to_string())
}

fn trapped_axis_energy(
    position: f64,
    strength: f64,
) -> Result<f64, &'static str> {
    if position == 0.0 {
        return Ok(0.0);
    }
    let position = position.abs();
    scale_product(strength, &[0.5, position, position], &[])
        .map_err(|_| "GPE trapping potential is not representable")
}

fn gpe_axis_damping(
    index: usize,
    extent: usize,
    length: f64,
    step: f64,
) -> f64 {
    let mode = std::f64::consts::TAU * index.min(extent - index) as f64;
    // A positive damping exponent too large for f64 gives a zero multiplier.
    scale_product(0.5 * mode * mode, &[step], &[length, length]).unwrap_or(f64::INFINITY)
}

fn gpe_potential_step(
    psi: &mut [Complex<f64>],
    potential: &[f64],
    coupling: f64,
    step: f64,
) {
    for_each_mut(psi, POINTWISE_GRAIN, |i, p| {
        let v = potential[i];
        let v_eff = if coupling == 0.0 {
            v
        } else {
            coupling.mul_add(p.norm_sqr(), v)
        };
        *p *= (-v_eff * step / 2.0).exp();
    });
}

fn gpe_normalize(
    psi: &mut [Complex<f64>],
    dx: f64,
    dy: f64,
) -> Result<(), String> {
    // A fixed association order makes the normalization independent of workers.
    let mass = |psi: &[Complex<f64>]| sum_fixed(psi.len(), 1_024, 32_768, |i| psi[i].norm_sqr());
    let mut norm = mass(psi);
    if !norm.is_normal() {
        let peak = psi
            .par_iter()
            .map(|p| {
                if p.re.is_finite() && p.im.is_finite() {
                    p.re.abs().max(p.im.abs())
                } else {
                    f64::INFINITY
                }
            })
            .reduce(|| 0.0, f64::max);
        if peak == 0.0 || !peak.is_finite() {
            return Err("GPE imaginary-time state has zero or non-finite norm".into());
        }
        for_each_mut(psi, POINTWISE_GRAIN, |_, p| *p /= peak);
        norm = mass(psi);
    }

    let row_norm = norm * dx;
    let cell_norm = row_norm * dy;
    let reciprocal = cell_norm.recip();
    let factor = if row_norm.is_normal() && cell_norm.is_normal() && reciprocal.is_normal() {
        reciprocal.sqrt()
    } else {
        scale_product(1.0, &[], &[norm.sqrt(), dx.sqrt(), dy.sqrt()])
            .map_err(|_| "GPE normalization factor is not representable")?
    };
    if factor == 0.0 {
        return Err("GPE normalization factor is not representable".into());
    }
    for_each_mut(psi, POINTWISE_GRAIN, |_, p| *p *= factor);
    Ok(())
}

struct GpeFft {
    shape: [usize; 2],
    batch_rows: bool,
    transposed: Option<Scratch<'static, Complex<f64>>>,
}

impl GpeFft {
    fn new(shape: [usize; 2]) -> Self {
        let [ny, nx] = shape;
        let size = ny * nx;
        // A mixed-radix rectangle can be dominated by its more expensive
        // axis. The shared wide FFT scheduler covers large grids.
        let mixed_axes = nx.is_power_of_two() != ny.is_power_of_two();
        let batch_axes = mixed_axes && nx > 1 && ny > 1 && (8_192..=131_072).contains(&size);
        let transposed = (batch_axes && !ny.is_power_of_two()).then(|| {
            let mut buffer = fft_scratch_for(size);
            buffer.resize(size, Complex::ZERO);
            buffer
        });
        Self {
            shape,
            batch_rows: batch_axes && !nx.is_power_of_two(),
            transposed,
        }
    }

    fn transform(
        &mut self,
        data: &mut [Complex<f64>],
        inverse: bool,
    ) {
        if !self.batch_rows && self.transposed.is_none() {
            transform_axes_in_place(data, &self.shape, &[1, 0], inverse, Norm::Backward);
            return;
        }

        // Use an unscaled transform for each axis and apply the same
        // full-grid inverse normalization as the shared two-axis path.
        let norm = if inverse {
            Norm::Forward
        } else {
            Norm::Backward
        };
        let [ny, nx] = self.shape;
        if self.batch_rows {
            gpe_transform_rows(data, nx, inverse, norm);
        } else {
            transform_axes_in_place(data, &self.shape, &[1], inverse, norm);
        }
        if let Some(work) = self.transposed.as_mut() {
            gpe_transpose(data, work, nx, ny);
            gpe_transform_rows(work, ny, inverse, norm);
            gpe_transpose(work, data, ny, nx);
        } else {
            transform_axes_in_place(data, &self.shape, &[0], inverse, norm);
        }
        if inverse {
            let scale = (data.len() as f64).recip();
            for_each_mut(data, POINTWISE_GRAIN, |_, p| *p *= scale);
        }
    }
}

fn gpe_transform_rows(
    data: &mut [Complex<f64>],
    width: usize,
    inverse: bool,
    norm: Norm,
) {
    let batch = width * (2_048 / width).max(1);
    for_chunks_mut(data, batch, 2, |_, rows| {
        transform_axes_in_place(rows, &[rows.len() / width, width], &[1], inverse, norm);
    });
}

fn gpe_transpose(
    data: &[Complex<f64>],
    output: &mut [Complex<f64>],
    width: usize,
    height: usize,
) {
    const COLUMNS: usize = 4;
    for_chunks_mut(output, COLUMNS * height, 4, |batch, columns| {
        for (offset, column) in columns.chunks_exact_mut(height).enumerate() {
            let x = batch * COLUMNS + offset;
            for (dest, &source) in column.iter_mut().zip(data[x..].iter().step_by(width)) {
                *dest = source;
            }
        }
    });
}

/// An example scenario that finds the ground state of a BEC, which may contain a vortex.
///
/// # Errors
///
/// This function will return an error if the GPE ground state finder fails
/// or if the final density cannot be written to the NPY file.
/// The density is saved as `gpe_vortex_state.npy` inside `output_dir`, which is
/// created if missing.
pub fn simulate_bose_einstein_vortex_scenario(output_dir: &Path) -> Result<(), String> {
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

    let final_density = run_gpe_ground_state_finder(&params)?;

    std::fs::create_dir_all(output_dir).map_err(|e| e.to_string())?;

    let filename = output_dir.join("gpe_vortex_state.npy");

    println!(
        "Simulation finished. Saving \
         final density to {}",
        filename.display()
    );

    write_npy_file(&filename, &final_density)?;

    Ok(())
}

/// GPE ground states and dynamics in 1D–3D: see [`crate::sim::split_step`].
pub use crate::sim::split_step::GroundState as GpeGroundState;
/// GPE ground states and dynamics in 1D–3D: see [`crate::sim::split_step`].
pub use crate::sim::split_step::SplitStep as GpeSolver;
