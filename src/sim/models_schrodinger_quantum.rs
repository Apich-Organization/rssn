//! # Schrodinger Equation Simulation
//!
//! This module provides a solver for the time-dependent Schrödinger equation, a fundamental equation
//! in quantum mechanics that describes how the quantum state of a physical system changes over time.
//!
//! # Overview
//!
//! The simulation solves the 2D Schrödinger equation on a grid using the Split-Step Fourier Method.
//! This method is particularly efficient because it handles the kinetic energy operator in the
//! frequency domain (using FFTs) and the potential energy operator in the spatial domain.
//!
//! Key features include:
//! - **Spectral Solver**: Utilizes Fast Fourier Transforms (FFT) for high accuracy and stability.
//! - **Arbitrary Potentials**: Supports user-defined potential energy landscapes (e.g., barriers, wells).
//! - **Wave Packet Dynamics**: Simulates the evolution of Gaussian wave packets.
//! - **Scenario Simulation**: Includes a pre-configured double-slit experiment scenario.
//!
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/quantum_tunneling_frame_001.png)
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/quantum_tunneling_frame_005.png)
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/quantum_tunneling_frame_007.png)
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/quantum_tunneling_frame_009.png)
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/quantum_tunneling_frame_010.png)
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/quantum_tunneling_frame_020.png)
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/quantum_tunneling_frame_025.png)
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/quantum_tunneling_frame_032.png)
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/quantum_tunneling_frame_033.png)
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/quantum_tunneling_frame_039.png)
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/quantum_tunneling_3d.png)
//!
//! # Coverage audit
//!
//! | Equation | Dimension | Boundaries | Method |
//! |---|---|---|---|
//! | time-dependent Schrödinger | 2D | periodic | split-step Fourier |
//!
//! [`crate::sim::split_step`] covers 1D–3D, the GPE nonlinearity, fourth-order
//! composition and imaginary-time ground states.

use ndarray::Array2;
use num_complex::Complex;
use serde::Deserialize;
use serde::Serialize;

use std::f64::consts::TAU;
use std::num::NonZeroUsize;
use std::path::Path;

use crate::io::write_npy_file;
use crate::kernels::fft::Norm;
use crate::kernels::fft::transform_axes_in_place;
use crate::kernels::numeric_common::scale_product;
use crate::kernels::par::for_chunks_mut;
use crate::kernels::par::for_each_mut;
use crate::kernels::reuse::fft_scratch_for;

/// Parameters for the Schrodinger simulation.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SchrodingerParameters {
    /// Number of grid points in the x-direction.
    pub nx: usize,
    /// Number of grid points in the y-direction.
    pub ny: usize,
    /// Domain length in the x-direction.
    pub lx: f64,
    /// Domain length in the y-direction.
    pub ly: f64,
    /// Time step size.
    pub dt: f64,
    /// Total number of simulation time steps.
    pub time_steps: usize,
    /// Reduced Planck constant (hbar).
    pub hbar: f64,
    /// Mass of the particle.
    pub mass: f64,
    /// Potential energy distribution (flattened).
    pub potential: Vec<f64>,
}

const PARALLEL_POINTS: usize = 32_768;

fn grid_size(params: &SchrodingerParameters) -> Result<usize, String> {
    if params.nx == 0
        || params.ny == 0
        || [params.lx, params.ly, params.hbar, params.mass]
            .iter()
            .any(|value| !value.is_finite() || *value <= 0.0)
        || !params.dt.is_finite()
    {
        return Err("invalid Schrodinger grid, lengths, mass, Planck constant or time step".into());
    }
    let size = params
        .nx
        .checked_mul(params.ny)
        .filter(|&n| n <= isize::MAX as usize / size_of::<Complex<f64>>())
        .ok_or("Schrodinger grid size overflows addressable storage")?;
    if params.potential.len() != size || params.potential.iter().any(|v| !v.is_finite()) {
        return Err("Schrodinger needs one finite potential sample per grid point".into());
    }
    Ok(size)
}

fn wave_peak(
    wave: &[Complex<f64>],
    size: usize,
) -> Result<f64, String> {
    if wave.len() != size {
        return Err("Schrodinger wave function does not match the grid".into());
    }
    let mut peak = 0.0_f64;
    for value in wave {
        if !value.re.is_finite() || !value.im.is_finite() {
            return Err("Schrodinger wave function must be finite".into());
        }
        peak = peak.max(value.re.abs()).max(value.im.abs());
    }
    Ok(peak)
}

fn reserve<T>(count: usize) -> Result<Vec<T>, String> {
    let mut buffer = Vec::new();
    buffer
        .try_reserve_exact(count)
        .map_err(|_| "Schrodinger storage allocation failed")?;
    Ok(buffer)
}

/// Prepared real-time propagation in a fixed real potential on a periodic 2D grid.
///
/// The model retains the physical mass and reduced Planck constant. The
/// potential and wave function are flattened with the x coordinate varying
/// fastest. Each evolution reuses the same Strang splitting phases.
#[derive(Debug)]
pub struct SchrodingerPlan2D {
    shape: [usize; 2],
    size: usize,
    packed_kinetic: bool,
    phases: Option<SchrodingerPhases>,
}

#[derive(Debug)]
struct SchrodingerPhases {
    kinetic: Vec<Complex<f64>>,
    potential: PotentialPhases,
}

#[derive(Debug)]
enum PotentialPhases {
    Spatial {
        half: Vec<Complex<f64>>,
        whole: Vec<Complex<f64>>,
    },
    Uniform {
        axes: [Vec<f64>; 2],
        kinetic_bound: f64,
        half_angle: f64,
        half: Complex<f64>,
        whole: Complex<f64>,
    },
}

impl PotentialPhases {
    fn apply(
        &self,
        wave: &mut [Complex<f64>],
        whole_step: bool,
    ) {
        match self {
            | Self::Spatial { half, whole } => {
                apply_phase(wave, if whole_step { whole } else { half });
            },
            | Self::Uniform {
                half_angle,
                half,
                whole,
                ..
            } => {
                if *half_angle != 0.0 {
                    let phase = if whole_step { *whole } else { *half };
                    for_each_mut(wave, PARALLEL_POINTS, |_, value| *value *= phase);
                }
            },
        }
    }
}

fn axis_phases(
    angles: &[f64],
    advance: f64,
) -> Result<Vec<Complex<f64>>, String> {
    let mut phases = reserve(angles.len())?;
    phases.extend(
        angles
            .iter()
            .map(|&angle| Complex::from_polar(1.0, angle * advance)),
    );
    Ok(phases)
}

impl SchrodingerPlan2D {
    pub(crate) fn heap_bytes(&self) -> usize {
        self.phases.as_ref().map_or(0, |phases| {
            let potential = match &phases.potential {
                | PotentialPhases::Spatial { half, whole } => {
                    (half.capacity() + whole.capacity()) * size_of::<Complex<f64>>()
                },
                | PotentialPhases::Uniform { axes, .. } => {
                    axes.iter().map(Vec::capacity).sum::<usize>() * size_of::<f64>()
                },
            };
            phases.kinetic.capacity() * size_of::<Complex<f64>>() + potential
        })
    }

    /// Prepares the signed time increment, grid and physical potential.
    ///
    /// The evolution call supplies the step count. The `time_steps` field
    /// controls the legacy snapshot-returning entry point.
    ///
    /// # Errors
    /// Returns an error for invalid geometry or coefficients, an
    /// unrepresentable phase angle, or failure to reserve phase storage.
    pub fn new(params: &SchrodingerParameters) -> Result<Self, String> {
        let size = grid_size(params)?;
        let shape = [params.ny, params.nx];
        if params.dt == 0.0 {
            return Ok(Self {
                shape,
                size,
                packed_kinetic: false,
                phases: None,
            });
        }

        let half_step = -0.5 * params.dt;
        let potential_factor = half_step / params.hbar;
        let direct_potential = half_step.is_normal() && potential_factor.is_normal();
        let angle = |value: f64| {
            if direct_potential {
                let phase = value * potential_factor;
                if phase.is_finite() {
                    return Ok(phase);
                }
            }
            scale_product(value, &[0.5, params.dt.abs()], &[params.hbar])
                .map(|phase| -params.dt.signum() * phase)
                .map_err(|_| "Schrodinger potential phase is not representable".to_string())
        };
        let v0 = params.potential[0];
        let uniform_angle = params
            .potential
            .iter()
            .all(|&v| crate::kernels::exact_eq(v, v0))
            .then(|| angle(v0))
            .transpose()?;
        let axis_angles = [
            Self::kinetic_axis(params.nx, params.lx, params)?,
            Self::kinetic_axis(params.ny, params.ly, params)?,
        ];
        let kx = axis_phases(&axis_angles[0], 1.0)?;
        let ky = axis_phases(&axis_angles[1], 1.0)?;
        let packed_kinetic = (8_192..=131_072).contains(&size)
            && params.nx > 1
            && params.ny > 1
            && !(params.nx.is_power_of_two() && params.ny.is_power_of_two());
        let mut kinetic = reserve(size)?;
        // The kinetic propagation separates into two commuting axes. Form
        // the two phases before multiplying to preserve a small axis phase
        // even when adding it to the other angle would round it away.
        let (inner, outer) = if packed_kinetic {
            (ky, kx)
        } else {
            (kx, ky)
        };
        for outer in outer {
            kinetic.extend(inner.iter().map(|&inner| inner * outer));
        }
        let potential = if let Some(half_angle) = uniform_angle {
            let kinetic_bound = axis_angles
                .iter()
                .flatten()
                .map(|value| value.abs())
                .fold(0.0, f64::max);
            let half = Complex::from_polar(1.0, half_angle);
            PotentialPhases::Uniform {
                axes: axis_angles,
                kinetic_bound,
                half_angle,
                half,
                whole: half * half,
            }
        } else {
            let mut half = reserve(size)?;
            for &value in &params.potential {
                half.push(Complex::from_polar(1.0, angle(value)?));
            }
            let mut whole = reserve(size)?;
            whole.extend(half.iter().map(|&phase| phase * phase));
            PotentialPhases::Spatial { half, whole }
        };
        Ok(Self {
            shape,
            size,
            packed_kinetic,
            phases: Some(SchrodingerPhases { kinetic, potential }),
        })
    }

    fn kinetic_axis(
        extent: usize,
        length: f64,
        params: &SchrodingerParameters,
    ) -> Result<Vec<f64>, String> {
        let mut result = reserve(extent)?;
        for i in 0..extent {
            let mode = TAU * i.min(extent - i) as f64;
            let angle = scale_product(
                0.5 * mode * mode,
                &[params.hbar, params.dt.abs()],
                &[params.mass, length, length],
            )
            .map_err(|_| "Schrodinger kinetic phase is not representable")?;
            result.push(-params.dt.signum() * angle);
        }
        Ok(result)
    }

    fn transform(
        &self,
        wave: &mut [Complex<f64>],
        inverse: bool,
    ) {
        transform_axes_in_place(wave, &self.shape, &[1, 0], inverse, Norm::Backward);
    }

    fn kinetic_step(
        &self,
        wave: &mut [Complex<f64>],
        kinetic: &[Complex<f64>],
        packed: Option<&mut [Complex<f64>]>,
    ) {
        if let Some(packed) = packed {
            let [ny, nx] = self.shape;
            // Keep the packed column layout across the kinetic multiplication
            // and the inverse column transform.
            transform_rows(wave, nx, false);
            transpose_wave(wave, packed, nx, ny);
            transform_rows(packed, ny, false);
            apply_phase(packed, kinetic);
            transform_rows(packed, ny, true);
            transpose_wave(packed, wave, ny, nx);
            transform_rows(wave, nx, true);
            let scale = 1.0 / self.size as f64;
            for_each_mut(wave, PARALLEL_POINTS, |_, value| *value *= scale);
        } else {
            self.transform(wave, false);
            apply_phase(wave, kinetic);
            self.transform(wave, true);
        }
    }

    /// Returns the final wave and optionally streams earlier wave samples.
    ///
    /// With a stride `s`, the observer sees completed steps `1, 1 + s, ...`
    /// within the requested interval. A zero step count emits no samples.
    /// Pass `None` for final-state-only evolution. A uniform potential can
    /// then be advanced directly in one pair of Fourier transforms.
    ///
    /// # Errors
    /// Invalid initial samples or shape, an unrepresentable evolved state,
    /// failed allocations, cancellation or an error from the observer.
    pub fn evolve(
        &self,
        initial: &[Complex<f64>],
        steps: usize,
        snapshot_stride: Option<NonZeroUsize>,
        mut observe: impl FnMut(usize, &[Complex<f64>]) -> Result<(), String>,
    ) -> Result<Vec<Complex<f64>>, String> {
        let peak = wave_peak(initial, self.size)?;
        let mut wave = reserve(self.size)?;
        wave.extend_from_slice(initial);
        let Some(phases) = self.phases.as_ref().filter(|_| steps != 0 && peak != 0.0) else {
            if let Some(stride) = snapshot_stride {
                for step in (0..steps).step_by(stride.get()) {
                    crate::cancellation::checkpoint().map_err(|error| error.to_string())?;
                    observe(step + 1, initial)?;
                }
            }
            return Ok(wave);
        };

        // A power of two keeps the Fourier intermediates in range while the
        // callback and return value retain the physical sample scale.
        let scale = (!(1e-140..=1e140).contains(&peak))
            .then(|| f64::from_bits(peak.to_bits() & 0x7ff0_0000_0000_0000).max(f64::MIN_POSITIVE));
        if let Some(scale) = scale {
            for value in &mut wave {
                *value /= scale;
            }
        }
        if steps > 1
            && snapshot_stride.is_none()
            && let PotentialPhases::Uniform {
                axes,
                kinetic_bound,
                half_angle,
                ..
            } = &phases.potential
            && steps as u128 <= (1_u128 << 53)
            && (kinetic_bound * steps as f64).is_finite()
            && (half_angle * (2.0 * steps as f64)).is_finite()
        {
            crate::cancellation::checkpoint().map_err(|error| error.to_string())?;
            let advance = steps as f64;
            let phase = Complex::from_polar(1.0, half_angle * (2.0 * advance));
            let kx = axis_phases(&axes[0], advance)?;
            let mut ky = axis_phases(&axes[1], advance)?;
            for value in &mut ky {
                *value *= phase;
            }
            let nx = self.shape[1];
            self.transform(&mut wave, false);
            for_each_mut(&mut wave, PARALLEL_POINTS, |i, value| {
                *value *= kx[i % nx] * ky[i / nx];
            });
            self.transform(&mut wave, true);
        } else {
            let mut packed = self.packed_kinetic.then(|| {
                let mut scratch = fft_scratch_for(self.size);
                scratch.resize(self.size, Complex::ZERO);
                scratch
            });
            let mut restored = Vec::new();
            phases.potential.apply(&mut wave, false);
            for index in 0..steps {
                crate::cancellation::checkpoint().map_err(|error| error.to_string())?;
                self.kinetic_step(
                    &mut wave,
                    &phases.kinetic,
                    packed.as_deref_mut().map(Vec::as_mut_slice),
                );
                let step = index + 1;
                let sampled =
                    snapshot_stride.is_some_and(|stride| index.is_multiple_of(stride.get()));
                if sampled || step == steps {
                    phases.potential.apply(&mut wave, false);
                    if sampled {
                        if let Some(scale) = scale {
                            restored.clear();
                            restored
                                .try_reserve(self.size)
                                .map_err(|_| "Schrodinger observer storage allocation failed")?;
                            restored.extend(wave.iter().map(|&value| value * scale));
                            wave_peak(&restored, self.size)?;
                            observe(step, &restored)?;
                        } else {
                            observe(step, &wave)?;
                        }
                    }
                    if step != steps {
                        phases.potential.apply(&mut wave, false);
                    }
                } else {
                    phases.potential.apply(&mut wave, true);
                }
            }
        }
        if let Some(scale) = scale {
            for value in &mut wave {
                *value *= scale;
            }
        }
        wave_peak(&wave, self.size)?;
        Ok(wave)
    }
}

fn transform_rows(
    data: &mut [Complex<f64>],
    width: usize,
    inverse: bool,
) {
    let batch = width * (2_048 / width).max(1);
    let norm = if inverse {
        Norm::Forward
    } else {
        Norm::Backward
    };
    for_chunks_mut(data, batch, 2, |_, rows| {
        transform_axes_in_place(rows, &[rows.len() / width, width], &[1], inverse, norm);
    });
}

fn transpose_wave(
    input: &[Complex<f64>],
    output: &mut [Complex<f64>],
    width: usize,
    height: usize,
) {
    const COLUMNS: usize = 8;
    for_chunks_mut(output, COLUMNS * height, 4, |batch, columns| {
        let first = batch * COLUMNS;
        for (row, samples) in input.chunks_exact(width).enumerate() {
            for (target, &value) in columns.chunks_exact_mut(height).zip(&samples[first..]) {
                target[row] = value;
            }
        }
    });
}

fn apply_phase(
    wave: &mut [Complex<f64>],
    phase: &[Complex<f64>],
) {
    for_each_mut(wave, PARALLEL_POINTS, |i, value| *value *= phase[i]);
}

/// Evolves the periodic 2D wave by Strang splitting and returns density frames.
///
/// Frames are sampled after steps 1, 11, 21, and so on. Each frame has `ny`
/// rows and `nx` columns. On success the supplied wave becomes the final state.
///
/// # Errors
/// Invalid geometry or samples, an unrepresentable phase, wave or density,
/// failed allocations, or cancellation.
pub fn run_schrodinger_simulation(
    params: &SchrodingerParameters,
    initial_psi: &mut [Complex<f64>],
) -> Result<Vec<Array2<f64>>, String> {
    if params.time_steps == 0 {
        wave_peak(initial_psi, grid_size(params)?)?;
        return Ok(Vec::new());
    }
    let plan = SchrodingerPlan2D::new(params)?;
    let mut snapshots = reserve(params.time_steps.div_ceil(10))?;
    let wave = plan.evolve(
        initial_psi,
        params.time_steps,
        NonZeroUsize::new(10),
        |_, wave| {
            let mut density = reserve(wave.len())?;
            for value in wave {
                let probability = value.norm_sqr();
                if !probability.is_finite() {
                    return Err("Schrodinger probability density is not representable".into());
                }
                density.push(probability);
            }
            let frame = Array2::from_shape_vec((params.ny, params.nx), density)
                .map_err(|error| error.to_string())?;
            snapshots.push(frame);
            Ok(())
        },
    )?;
    initial_psi.copy_from_slice(&wave);
    Ok(snapshots)
}

/// An example scenario simulating a wave packet hitting a double slit.
///
/// # Errors
///
/// This function will return an error if the underlying `run_schrodinger_simulation` fails
/// or if the final state cannot be written to the NPY file.
/// The last sampled density is saved as `schrodinger_double_slit.npy` inside
/// `output_dir`, which is created if missing.
pub fn simulate_double_slit_scenario(output_dir: &Path) -> Result<(), String> {
    const NX: usize = 256;

    const NY: usize = 256;

    println!(
        "Running 2D Schrodinger \
         simulation for a double \
         slit..."
    );

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
            potential[j * NX + barrier_x] = 1e5;
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

    let mut initial_psi = vec![Complex::default(); NX * NY];

    let initial_pos = (NX as f64 / 10.0, NY as f64 / 2.0);

    let initial_momentum = (5.0_f64, 0.0_f64);

    let packet_width_sq = 100.0;

    for j in 0..NY {
        for i in 0..NX {
            let dx = i as f64 - initial_pos.0;

            let dy = j as f64 - initial_pos.1;

            let norm_sq = dx * dx + dy * dy;

            let phase = initial_momentum.0.mul_add(dx, initial_momentum.1 * dy);

            let envelope = (-norm_sq / (2.0 * packet_width_sq)).exp();

            initial_psi[j * NX + i] = Complex::from_polar(envelope, phase);
        }
    }

    let snapshots = run_schrodinger_simulation(&params, &mut initial_psi)?;

    if let Some(final_state) = snapshots.last() {
        std::fs::create_dir_all(output_dir).map_err(|e| e.to_string())?;

        let filename = output_dir.join("schrodinger_double_slit.npy");

        println!(
            "Saving final probability \
             density to {}",
            filename.display()
        );

        write_npy_file(&filename, final_state)?;
    } else {
        println!(
            "Simulation produced no \
             snapshots."
        );
    }

    Ok(())
}

/// Split-step Fourier solver in 1D–3D with real- and imaginary-time propagation: see [`crate::sim::split_step`].
pub use crate::sim::split_step::GroundState;
/// Split-step Fourier solver in 1D–3D with real- and imaginary-time propagation: see [`crate::sim::split_step`].
pub use crate::sim::split_step::SplitOrder;
/// Split-step Fourier solver in 1D–3D with real- and imaginary-time propagation: see [`crate::sim::split_step`].
pub use crate::sim::split_step::SplitStep;
