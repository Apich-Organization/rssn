//! #  Ising Model Simulation
//!
//! This module implements a Monte Carlo simulation of the 2D Ising model, a mathematical model
//! of ferromagnetism in statistical mechanics. It uses the Metropolis-Hastings algorithm to explore
//! the configuration space of spins on a grid.
//!
//! # Overview
//!
//! The simulation models a lattice of magnetic spins that can be in one of two states (+1 or -1).
//! The system evolves towards thermal equilibrium by randomly flipping spins based on the energy cost
//! and the temperature (Metropolis criterion).
//!
//! Key features include:
//! - **Monte Carlo Method**: Uses the Metropolis algorithm to sample states according to the Boltzmann distribution.
//! - **Parallel Update**: Implements a checkerboard (Red-Black) update scheme to allow parallel updates of non-interacting spins.
//! - **Phase Transition**: Can simulate the system across a range of temperatures to observe the transition from ordered (ferromagnetic) to disordered (paramagnetic) phases.
//! - **Magnetization Tracking**: Calculates the average magnetization of the system to characterize the phase.
//!
//! ![refer to this image](https://raw.githubusercontent.com/Apich-Organization/rssn/refs/heads/dev/doc/ising_phase_transition.png)
//!
//! ## Coverage
//!
//! | Capability | Support |
//! |---|---|
//! | Dimension | `run_ising_simulation*`: 2-D parallel checkerboard Metropolis; [`run_lattice_simulation`] and `run_wolff_simulation` delegate to the dimension-generic engine of [`crate::sim::spin_models`] (any `D`) |
//! | Models | Ising, q-state Potts, XY, Heisenberg ([`SpinModel`]) through [`run_lattice_simulation`] |
//! | Algorithms | Metropolis (checkerboard, 2-D Ising), Metropolis / heat bath / Wolff / Swendsen-Wang for every model through [`SpinAlgorithm`]; tempering, Wang-Landau, reweighting in [`crate::sim::spin_models`] |
//! | Boundaries | periodic (checkerboard with sequential odd edges or single-site sweeps for axes of length one or two), periodic or open (`run_lattice_simulation`) |
//! | Initial conditions | random (checkerboard), hot or cold (`run_lattice_simulation`) |
//! | Exact references | Onsager `Tc`, energy and magnetisation in [`crate::sim::spin_models`] |

use std::fmt::Write as OtherWrite;
use std::fs::File;
use std::io::Write;
use std::num::NonZeroUsize;
use std::path::Path;

use ndarray::Array2;
use rand_v10::prelude::*;
use rand_v10::rngs::StdRng;
use rayon::prelude::*;
use serde::Deserialize;
use serde::Serialize;

use crate::io::write_npy_file;
use crate::sim::spin_models::OnModel;
use crate::sim::spin_models::Potts;
use crate::sim::spin_models::Update;
use crate::sim::spin_models::binder_cumulant;
use crate::sim::spin_models::specific_heat;
use crate::sim::spin_models::susceptibility;

/// Parameters for the Ising model simulation.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct IsingParameters {
    /// The width of the simulation grid.
    pub width: usize,
    /// The height of the simulation grid.
    pub height: usize,
    /// The temperature of the system.
    pub temperature: f64,
    /// The number of Monte Carlo steps to perform.
    pub mc_steps: usize,
}

/// Seed used by [`run_ising_simulation`], so that plain runs are reproducible.
pub const DEFAULT_ISING_SEED: u64 = 0x1517_0000_5EED;

/// Derives an independent RNG for one (step, colour, row) of the checkerboard sweep.
fn row_rng(
    seed: u64,
    step: usize,
    colour: u64,
    row: usize,
) -> StdRng {
    let mixed = seed
        ^ (step as u64)
            .wrapping_mul(2)
            .wrapping_add(colour)
            .wrapping_mul(0x9E37_79B9_7F4A_7C15)
        ^ (row as u64)
            .wrapping_add(1)
            .wrapping_mul(0xBF58_476D_1CE4_E5B9);

    StdRng::seed_from_u64(mixed)
}

/// Runs an Ising model simulation with the fixed seed [`DEFAULT_ISING_SEED`].
///
/// Equivalent to [`run_ising_simulation_seeded`]; repeated calls with the same
/// parameters return identical results.
#[must_use]
pub fn run_ising_simulation(params: &IsingParameters) -> (Vec<i8>, f64) {
    run_ising_simulation_seeded(params, DEFAULT_ISING_SEED)
}

/// Runs an Ising model simulation whose random numbers all derive from `seed`.
///
/// The initial configuration and every row update use their own generator
/// derived from `seed`, so the result does not depend on thread scheduling.
/// Grids with an axis of length one or two use ordered single-site sweeps.
/// At positive temperature on a multi-site grid, each visit in such a
/// sweep proposes a flip with probability one half.
#[must_use]
pub fn run_ising_simulation_seeded(
    params: &IsingParameters,
    seed: u64,
) -> (Vec<i8>, f64) {
    let Ok(mut model) = Checkerboard::new(params, seed) else {
        return (Vec::new(), f64::NAN);
    };
    for step in 0..params.mc_steps {
        model.sweep(step);
    }
    let sum: i64 = model.spins.iter().map(|&s| i64::from(s)).sum();
    let magnetization = (sum as f64 / model.spins.len() as f64).abs();
    (model.spins, magnetization)
}

#[derive(Clone, Copy)]
struct MetropolisWeights {
    uphill: [f64; 2],
}

impl MetropolisWeights {
    fn new(temperature: f64) -> Self {
        let beta = temperature.abs().recip();
        Self {
            uphill: [(-4.0 * beta).exp(), (-8.0 * beta).exp()],
        }
    }

    fn update(
        self,
        spin: i8,
        neighbors: i8,
        rng: &mut StdRng,
    ) -> i8 {
        let delta = 2 * spin * neighbors;
        if delta < 0 {
            return -spin;
        }
        let accept = match delta {
            | 0 => 1.0,
            | 4 => self.uphill[0],
            | _ => self.uphill[1],
        };
        // Nonnegative energy changes always draw from the row's stream,
        // including a zero-cost move at zero temperature.
        if rng.random::<f64>() < accept {
            -spin
        } else {
            spin
        }
    }
}

struct Checkerboard {
    spins: Vec<i8>,
    next: Vec<i8>,
    width: usize,
    height: usize,
    inner_width: usize,
    inner_height: usize,
    weights: MetropolisWeights,
    seed: u64,
    parallel: bool,
    lazy_narrow: bool,
}

const fn previous(
    index: usize,
    size: usize,
) -> usize {
    if index == 0 {
        size - 1
    } else {
        index - 1
    }
}

const fn following(
    index: usize,
    size: usize,
) -> usize {
    if index + 1 == size {
        0
    } else {
        index + 1
    }
}

impl Checkerboard {
    fn new(
        params: &IsingParameters,
        seed: u64,
    ) -> Result<Self, &'static str> {
        const MIN_PARALLEL_SITES: usize = 16_384;

        let width = params.width;
        let height = params.height;
        let sites = width
            .checked_mul(height)
            .filter(|&n| n > 0)
            .ok_or("invalid Ising grid shape")?;
        if params.temperature.is_nan() || params.temperature < 0.0 {
            return Err("Ising temperature must be nonnegative");
        }

        let mut spins = Vec::new();
        spins
            .try_reserve_exact(sites)
            .map_err(|_| "cannot allocate Ising spins")?;
        let mut local_rng = StdRng::seed_from_u64(seed);
        spins.extend((0..sites).map(|_| {
            if local_rng.random::<bool>() {
                1
            } else {
                -1
            }
        }));

        let inner_width = width - width % 2;
        let inner_height = height - height % 2;
        let mut next = Vec::new();
        let parallel = if params.mc_steps > 0 && width > 2 && height > 2 {
            let slots = inner_height * (inner_width / 2);
            next.try_reserve_exact(slots)
                .map_err(|_| "cannot allocate Ising sweep buffer")?;
            next.resize(slots, 0);
            crate::execution::should_split(inner_width * inner_height, MIN_PARALLEL_SITES)
        } else {
            false
        };
        Ok(Self {
            spins,
            next,
            width,
            height,
            inner_width,
            inner_height,
            weights: MetropolisWeights::new(params.temperature),
            seed,
            parallel,
            lazy_narrow: sites > 1 && params.temperature > 0.0,
        })
    }

    fn sweep(
        &mut self,
        step: usize,
    ) {
        if self.width <= 2 || self.height <= 2 {
            // Repeated neighbors on thin tori can leave closed classes
            // under unconditional systematic proposals. Skipping some
            // proposals keeps the single-site Gibbs balance and lets the
            // sweep reach configurations that would otherwise be skipped.
            let mut rng = row_rng(self.seed, step, 2, 0);
            for row in 0..self.height {
                for column in 0..self.width {
                    if !self.lazy_narrow || rng.random::<bool>() {
                        self.update_site(row, column, &mut rng);
                    }
                }
            }
            return;
        }
        if !self.next.is_empty() {
            self.sweep_colour(step, 0);
            self.sweep_colour(step, 1);
        }

        // The two colours leave any extra odd row or column untouched.
        // Sweep those sites in their original flat order.
        if self.inner_width < self.width || self.inner_height < self.height {
            let mut rng = row_rng(self.seed, step, 2, 0);
            if self.inner_width < self.width {
                for row in 0..self.inner_height {
                    self.update_site(row, self.inner_width, &mut rng);
                }
            }
            if self.inner_height < self.height {
                for column in 0..self.width {
                    self.update_site(self.inner_height, column, &mut rng);
                }
            }
        }
    }

    fn sweep_colour(
        &mut self,
        step: usize,
        colour: usize,
    ) {
        let width = self.width;
        let height = self.height;
        let seed = self.seed;
        let weights = self.weights;
        let grid = &self.spins;
        let columns = self.inner_width / 2;

        let update_row = |row, output: &mut [i8]| {
            let mut rng = row_rng(seed, step, colour as u64, row);
            let first = row * width;
            let up = previous(row, height) * width;
            let down = following(row, height) * width;
            let current = &grid[first..first + width];
            let top = &grid[up..up + width];
            let bottom = &grid[down..down + width];
            let start = (row % 2) ^ colour;
            for (index, new) in output.iter_mut().enumerate() {
                let column = start + 2 * index;
                let neighbors = top[column]
                    + bottom[column]
                    + current[previous(column, width)]
                    + current[following(column, width)];
                *new = weights.update(current[column], neighbors, &mut rng);
            }
        };

        // Every row reads the same pre-update configuration. Reusing the
        // half-size buffer avoids allocating one flip list per row and colour.
        if self.parallel {
            crate::kernels::par::for_chunks_mut(&mut self.next, columns, 2, update_row);
        } else {
            for (row, output) in self.next.chunks_exact_mut(columns).enumerate() {
                update_row(row, output);
            }
        }
        for (row, (spins, updated)) in self
            .spins
            .chunks_exact_mut(width)
            .zip(self.next.chunks_exact(columns))
            .enumerate()
        {
            let start = (row % 2) ^ colour;
            for (spin, &new) in spins[start..self.inner_width]
                .iter_mut()
                .step_by(2)
                .zip(updated)
            {
                *spin = new;
            }
        }
    }

    fn update_site(
        &mut self,
        row: usize,
        column: usize,
        rng: &mut StdRng,
    ) {
        let width = self.width;
        let height = self.height;
        let start = row * width;
        let idx = start + column;
        let v = &mut self.spins;
        // Self-bonds on a unit axis have constant energy. The two bonds on
        // an extent-two axis remain distinct, although they share endpoints.
        let vertical = if height == 1 {
            0
        } else {
            v[previous(row, height) * width + column] + v[following(row, height) * width + column]
        };
        let horizontal = if width == 1 {
            0
        } else {
            v[start + previous(column, width)] + v[start + following(column, width)]
        };
        v[idx] = self.weights.update(v[idx], vertical + horizontal, rng);
    }

    fn energy_and_magnetization(&self) -> (f64, f64) {
        let mut horizontal = 0_i64;
        let mut vertical = 0_i64;
        let mut magnetization = 0_i64;
        let width = self.width;
        for (row, spins) in self.spins.chunks_exact(width).enumerate() {
            let down = following(row, self.height) * width;
            let below = &self.spins[down..down + width];
            let right = spins.iter().skip(1).chain(&spins[..1]);
            for ((&spin, &right), &bottom) in spins.iter().zip(right).zip(below) {
                horizontal += i64::from(spin * right);
                vertical += i64::from(spin * bottom);
                magnetization += i64::from(spin);
            }
        }
        let n = self.spins.len() as f64;
        (
            -(horizontal as f64 + vertical as f64) / n,
            (magnetization as f64 / n).abs(),
        )
    }
}

#[derive(Default)]
struct Moment {
    mean: f64,
    centered_square: f64,
}

impl Moment {
    fn push(
        &mut self,
        value: f64,
        count: f64,
    ) {
        let shift = value - self.mean;
        self.mean += shift / count;
        self.centered_square += shift * (value - self.mean);
    }
}

/// Measures a seeded periodic 2D checkerboard Ising run after thermalization.
///
/// Discards `thermalisation` sweeps, then runs `params.mc_steps` more.
/// Measures after every `sample_every` sweeps in the latter interval and
/// returns the final spins alongside the per-site averages. Incomplete
/// sampling intervals still advance the spins. Measurements use constant
/// extra storage and do not consume random draws.
///
/// The bond energy counts the right and lower bond of every spin. Unit
/// axes therefore contribute a constant self-bond energy, and extent-two
/// axes have two bonds between the same pair. Susceptibility uses the
/// variance of absolute magnetization. Binder is NaN when all measured
/// magnetizations vanish.
///
/// # Errors
///
/// Returns an error for an invalid lattice, an overflowing sweep total,
/// no complete measurement interval, or a nonpositive temperature or
/// nonfinite inverse temperature.
pub fn measure_ising_simulation_seeded(
    params: &IsingParameters,
    seed: u64,
    thermalisation: usize,
    sample_every: NonZeroUsize,
) -> Result<(Vec<i8>, LatticeObservables), String> {
    let beta = params.temperature.recip();
    if params.temperature <= 0.0 || !beta.is_finite() {
        return Err("Ising measurements require positive temperature with finite inverse".into());
    }
    let stride = sample_every.get();
    if params.mc_steps < stride {
        return Err("Ising measurements require a complete sampling interval".into());
    }
    let end = thermalisation
        .checked_add(params.mc_steps)
        .ok_or("Ising sweep count overflow")?;
    let mut model = Checkerboard::new(params, seed)?;
    for step in 0..thermalisation {
        model.sweep(step);
    }

    let (mut energy, mut magnetization) = (Moment::default(), Moment::default());
    let (mut m2, mut m4) = (0.0, 0.0);
    let mut samples = 0.0;
    let mut until_sample = stride;
    for step in thermalisation..end {
        model.sweep(step);
        until_sample -= 1;
        if until_sample == 0 {
            until_sample = stride;
            let (e, m) = model.energy_and_magnetization();
            samples += 1.0;
            energy.push(e, samples);
            magnetization.push(m, samples);
            let square = m * m;
            m2 += square;
            m4 += square * square;
        }
    }
    let n = model.spins.len() as f64;
    m2 /= samples;
    m4 /= samples;
    let result = LatticeObservables {
        energy: energy.mean,
        magnetization: magnetization.mean,
        specific_heat: (energy.centered_square / samples) * beta * beta * n,
        susceptibility: (magnetization.centered_square / samples) * beta * n,
        binder: 1.0 - m4 / (3.0 * m2 * m2),
    };
    Ok((model.spins, result))
}

/// An example scenario that simulates the Ising model across a range of temperatures
/// to observe the phase transition.
///
/// # Errors
///
/// This function will return an error if it fails to reshape the `Array2<f64>` for NPY
/// output or if it fails to create or write to the output files (CSV or NPY).
/// The files `ising_low_temp_state.npy`, `ising_high_temp_state.npy` and
/// `ising_magnetization_vs_temp.csv` are written into `output_dir`, which is
/// created if missing.
pub fn simulate_ising_phase_transition_scenario(output_dir: &Path) -> Result<(), String> {
    std::fs::create_dir_all(output_dir).map_err(|e| e.to_string())?;

    println!(
        "Running Ising model phase \
         transition simulation..."
    );

    let temperatures: Vec<f64> = (0..=40).map(|i| f64::from(i).mul_add(0.1, 0.1)).collect();

    let scenario_results: Vec<(f64, f64, Vec<i8>)> = temperatures
        .par_iter()
        .map(|&temp| {
            let params = IsingParameters {
                width: 50,
                height: 50,
                temperature: temp,
                mc_steps: 2000,
            };

            let (grid, mag) = run_ising_simulation(&params);

            (temp, mag, grid)
        })
        .collect();

    let mut results = String::from("temperature,magnetization\n");

    for (i, (temp, mag, grid)) in scenario_results.iter().enumerate() {
        writeln!(results, "{temp},{mag}").expect("String transition failed.");

        if i == 5 {
            let arr: Array2<f64> =
                Array2::from_shape_vec((50, 50), grid.iter().map(|&s| f64::from(s)).collect())
                    .map_err(|e| e.to_string())?;

            write_npy_file(output_dir.join("ising_low_temp_state.npy"), &arr)?;
        }

        if i == 35 {
            let arr: Array2<f64> =
                Array2::from_shape_vec((50, 50), grid.iter().map(|&s| f64::from(s)).collect())
                    .map_err(|e| e.to_string())?;

            write_npy_file(output_dir.join("ising_high_temp_state.npy"), &arr)?;
        }
    }

    let mut file = File::create(output_dir.join("ising_magnetization_vs_temp.csv"))
        .map_err(|e| e.to_string())?;

    file.write_all(results.as_bytes())
        .map_err(|e| e.to_string())?;

    println!(
        "Saved results to CSV and NPY \
         files."
    );

    Ok(())
}

/// Thermodynamic averages of a Monte Carlo run (per spin, `J = k_B = 1`).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct IsingObservables {
    /// `⟨E⟩ / N`.
    pub energy: f64,
    /// `⟨|M|⟩ / N`.
    pub magnetization: f64,
    /// Specific heat `(⟨E²⟩ - ⟨E⟩²) / (N T²)`.
    pub specific_heat: f64,
    /// Susceptibility `(⟨M²⟩ - ⟨|M|⟩²) / (N T)`.
    pub susceptibility: f64,
    /// Binder cumulant `1 - ⟨M⁴⟩ / (3 ⟨M²⟩²)`.
    pub binder: f64,
    /// Mean Wolff cluster size (fraction of the lattice).
    pub mean_cluster: f64,
}

/// The Wolff single-cluster algorithm:
///
/// a cluster grown from a random seed spin, each aligned neighbour added
/// with probability `1 - exp(-2/T)`, is flipped as a whole. Rejection-free
/// and free of critical slowing down, so it samples near `T_c` far better
/// than local Metropolis updates. `sweeps` cluster updates are measured
/// after `thermalisation` discarded ones.
///
/// Runs on the dimension-generic engine [`crate::sim::spin_models::OnModel`] (cold start,
/// periodic `width x height` lattice); the observables are `NaN` for lattices the engine
/// rejects (periodic sides below 3).
#[must_use]
pub fn run_wolff_simulation(
    width: usize,
    height: usize,
    temperature: f64,
    thermalisation: usize,
    sweeps: usize,
    seed: u64,
) -> IsingObservables {
    let nan = IsingObservables {
        energy: f64::NAN,
        magnetization: f64::NAN,
        specific_heat: f64::NAN,
        susceptibility: f64::NAN,
        binder: f64::NAN,
        mean_cluster: f64::NAN,
    };
    let Ok(mut model) = OnModel::<2, 1>::new([width, height], true, 1.0 / temperature, seed) else {
        return nan;
    };
    model.cold_start();
    let nf = (width * height) as f64;
    for _ in 0..thermalisation {
        if model.wolff_step().is_err() {
            return nan;
        }
    }
    let (mut e, mut m, mut cl) = (Vec::new(), Vec::new(), 0.0);
    for _ in 0..sweeps {
        match model.wolff_step() {
            | Ok(size) => cl += size as f64,
            | Err(_) => return nan,
        }
        e.push(model.energy());
        m.push(model.magnetisation()[0].abs());
    }
    let samples = sweeps.max(1) as f64;
    let mean = |v: &[f64]| v.iter().sum::<f64>() / samples;
    let mps: Vec<f64> = m.iter().map(|x| x / nf).collect();
    let m1 = mean(&m);
    let m2 = m.iter().map(|x| x * x).sum::<f64>() / samples;
    IsingObservables {
        energy: mean(&e) / nf,
        magnetization: m1 / nf,
        specific_heat: specific_heat(&e, 1.0 / temperature, width * height),
        susceptibility: (m2 - m1 * m1) / (nf * temperature),
        binder: binder_cumulant(&mps),
        mean_cluster: cl / samples / nf,
    }
}

/// Spin-model family of [`run_lattice_simulation`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SpinModel {
    /// Ising (`N = 1`).
    Ising,
    /// q-state Potts.
    Potts(usize),
    /// Planar rotor (O(2)).
    Xy,
    /// Classical Heisenberg (O(3)).
    Heisenberg,
}

/// Monte Carlo update of [`run_lattice_simulation`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SpinAlgorithm {
    /// Single-spin Metropolis with the given proposal width (ignored for Ising / Potts).
    Metropolis(f64),
    /// Heat bath (exact conditional sampling; all four models).
    HeatBath,
    /// Wolff single cluster.
    Wolff,
    /// Swendsen-Wang.
    SwendsenWang,
}

/// Per-site thermodynamic averages returned by [`run_lattice_simulation`] (`J = k_B = 1`).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LatticeObservables {
    /// `<E> / N`.
    pub energy: f64,
    /// `<|m|>` (magnetisation per site; Potts order parameter).
    pub magnetization: f64,
    /// Specific heat per site.
    pub specific_heat: f64,
    /// Susceptibility per site.
    pub susceptibility: f64,
    /// Binder cumulant.
    pub binder: f64,
}

fn average_series(
    e: &[f64],
    m: &[f64],
    beta: f64,
    sites: usize,
) -> LatticeObservables {
    let n = e.len().max(1) as f64;
    LatticeObservables {
        energy: e.iter().sum::<f64>() / n / sites as f64,
        magnetization: m.iter().map(|x| x.abs()).sum::<f64>() / n,
        specific_heat: specific_heat(e, beta, sites),
        susceptibility: susceptibility(m, beta, sites),
        binder: binder_cumulant(m),
    }
}

fn run_on<const D: usize, const N: usize>(
    dims: [usize; D],
    periodic: bool,
    algo: SpinAlgorithm,
    temperature: f64,
    (thermalisation, sweeps): (usize, usize),
    seed: u64,
) -> Result<LatticeObservables, String> {
    let beta = 1.0 / temperature;
    let mut model = OnModel::<D, N>::new(dims, periodic, beta, seed)?;
    let update = match algo {
        | SpinAlgorithm::Metropolis(d) => Update::Metropolis(d),
        | SpinAlgorithm::HeatBath => Update::HeatBath,
        | SpinAlgorithm::Wolff => Update::Wolff,
        | SpinAlgorithm::SwendsenWang => Update::SwendsenWang,
    };
    for _ in 0..thermalisation {
        model.update(update)?;
    }
    let sites = model.lattice.sites();
    let (mut e, mut m) = (Vec::new(), Vec::new());
    for _ in 0..sweeps {
        model.update(update)?;
        e.push(model.energy());
        m.push(model.abs_magnetisation_per_site());
    }
    Ok(average_series(&e, &m, beta, sites))
}

/// Monte Carlo of a spin model on a hypercubic lattice of any dimension.
///
/// Ising, Potts, XY or Heisenberg on a lattice of **any**
/// dimension `D` (periodic or open), delegating to [`crate::sim::spin_models`]. Returns the
/// per-site averages over `sweeps` measurements after `thermalisation` updates.
///
/// # Errors
/// Fails for invalid lattices, `temperature <= 0` or unsupported model/algorithm combinations
/// (for example cluster updates with a field).
pub fn run_lattice_simulation<const D: usize>(
    dims: [usize; D],
    periodic: bool,
    model: SpinModel,
    algo: SpinAlgorithm,
    temperature: f64,
    steps: (usize, usize),
    seed: u64,
) -> Result<LatticeObservables, String> {
    if temperature <= 0.0 {
        return Err("temperature must be positive".into());
    }
    match model {
        | SpinModel::Ising => run_on::<D, 1>(dims, periodic, algo, temperature, steps, seed),
        | SpinModel::Xy => run_on::<D, 2>(dims, periodic, algo, temperature, steps, seed),
        | SpinModel::Heisenberg => run_on::<D, 3>(dims, periodic, algo, temperature, steps, seed),
        | SpinModel::Potts(q) => {
            let beta = 1.0 / temperature;
            let mut p = Potts::<D>::new(dims, periodic, q, beta, seed)?;
            let update = match algo {
                | SpinAlgorithm::Metropolis(d) => Update::Metropolis(d),
                | SpinAlgorithm::HeatBath => Update::HeatBath,
                | SpinAlgorithm::Wolff => Update::Wolff,
                | SpinAlgorithm::SwendsenWang => Update::SwendsenWang,
            };
            for _ in 0..steps.0 {
                p.update(update);
            }
            let sites = p.spins.len();
            let (mut e, mut m) = (Vec::new(), Vec::new());
            for _ in 0..steps.1 {
                p.update(update);
                e.push(p.energy());
                m.push(p.order_parameter());
            }
            Ok(average_series(&e, &m, beta, sites))
        },
    }
}

#[cfg(test)]
mod wolff_tests {
    use super::*;

    #[test]
    fn ordered_and_disordered_phases() {
        // T_c = 2/ln(1 + √2) ≈ 2.269.
        let cold = run_wolff_simulation(24, 24, 1.5, 200, 1000, 7);
        let hot = run_wolff_simulation(24, 24, 4.0, 200, 3000, 7);
        assert!(cold.magnetization > 0.95 && cold.binder > 0.6, "{cold:?}");
        assert!(hot.magnetization < 0.3 && hot.binder < 0.4, "{hot:?}");
        // Onsager's energy per spin at T = 1.5 is ≈ -1.947.
        assert!((cold.energy + 1.947).abs() < 0.02, "{cold:?}");
        // The specific heat peaks near T_c.
        let near = run_wolff_simulation(24, 24, 2.27, 200, 3000, 7);
        assert!(
            near.specific_heat > cold.specific_heat && near.specific_heat > hot.specific_heat,
            "{near:?}"
        );
    }
}
