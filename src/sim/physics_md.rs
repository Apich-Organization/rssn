//! # Numerical Molecular Dynamics (MD)
//!
//! This module provides numerical methods for Molecular Dynamics (MD) simulations,
//! a computational method for analyzing the physical movements of atoms and molecules.
//!
//! ## Features
//!
//! ### Interatomic Potentials
//! - **Lennard-Jones**: Standard 12-6 potential for noble gases
//! - **Morse potential**: Better for covalent bonds
//! - **Harmonic potential**: Simple spring model
//! - **Coulomb potential**: Electrostatic interactions
//!
//! ### Integration Algorithms
//! - **Velocity Verlet**: Standard symplectic integrator
//! - **Leapfrog**: Alternative formulation
//! - **Euler**: Simple first-order (reference only)
//!
//! ### Thermodynamics
//! - Temperature calculation from kinetic energy
//! - Pressure calculation from virial
//! - Thermostats (velocity rescaling, Berendsen)
//!
//! ### Analysis
//! - Radial distribution function
//! - Mean square displacement
//! - Kinetic and potential energy
//!
//! ## Example
//!
//! ```rust
//! use rssn::sim::physics_md::*;
//!
//! let p1 = Particle::new(0, 1.0, vec![0.0, 0.0, 0.0], vec![0.0, 0.0, 0.0]);
//!
//! let p2 = Particle::new(1, 1.0, vec![1.5, 0.0, 0.0], vec![0.0, 0.0, 0.0]);
//!
//! let (potential, force) = lennard_jones_interaction(&p1, &p2, 1.0, 1.0).unwrap();
//! ```

use serde::Deserialize;
use serde::Serialize;

use crate::kernels::vector::norm;
use crate::kernels::vector::scalar_mul;
use crate::kernels::vector::vec_add;
use crate::kernels::vector::vec_sub;

// ============================================================================
// Physical Constants (Reduced Units)
// ============================================================================

/// Boltzmann constant in SI units (J/K)
pub const BOLTZMANN_CONSTANT_SI: f64 = 1.380_649e-23;

/// Avogadro's number (1/mol)
pub const AVOGADRO_NUMBER: f64 = 6.022_140_76e23;

/// Reduced unit for temperature (using argon as reference)
/// 1 reduced temperature = `ε/k_B` ≈ 120 K for argon
pub const TEMPERATURE_UNIT_ARGON: f64 = 119.8;

/// Reduced unit for length (using argon as reference)
/// 1 reduced length = σ ≈ 3.4 Å for argon
pub const LENGTH_UNIT_ARGON: f64 = 3.4e-10;

/// Reduced unit for energy (using argon as reference)
/// 1 reduced energy = ε ≈ 1.65e-21 J for argon
pub const ENERGY_UNIT_ARGON: f64 = 1.65e-21;

// ============================================================================
// Particle Definition
// ============================================================================

/// Represents a particle in a molecular dynamics simulation.
/// cbindgen:ignore
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Particle {
    /// Unique identifier
    pub id: usize,
    /// Mass in reduced units
    pub mass: f64,
    /// Position vector
    pub position: Vec<f64>,
    /// Velocity vector
    pub velocity: Vec<f64>,
    /// Force vector (computed during simulation)
    pub force: Vec<f64>,
    /// Charge (for Coulomb interactions)
    pub charge: f64,
    /// Particle type (for multi-species simulations)
    pub particle_type: usize,
}

impl Particle {
    /// Creates a new particle with default charge and type.
    #[must_use]
    pub fn new(
        id: usize,
        mass: f64,
        position: Vec<f64>,
        velocity: Vec<f64>,
    ) -> Self {
        let dim = position.len();

        Self {
            id,
            mass,
            position,
            velocity,
            force: vec![0.0; dim],
            charge: 0.0,
            particle_type: 0,
        }
    }

    /// Creates a new particle with charge.
    #[must_use]
    pub fn with_charge(
        id: usize,
        mass: f64,
        position: Vec<f64>,
        velocity: Vec<f64>,
        charge: f64,
    ) -> Self {
        let dim = position.len();

        Self {
            id,
            mass,
            position,
            velocity,
            force: vec![0.0; dim],
            charge,
            particle_type: 0,
        }
    }

    /// Returns the kinetic energy of the particle: KE = 0.5 * m * v²
    #[must_use]
    pub fn kinetic_energy(&self) -> f64 {
        let v2: f64 = self.velocity.iter().map(|v| v * v).sum();

        0.5 * self.mass * v2
    }

    /// Returns the momentum of the particle: p = m * v
    #[must_use]
    pub fn momentum(&self) -> Vec<f64> {
        scalar_mul(&self.velocity, self.mass)
    }

    /// Returns the speed (magnitude of velocity)
    #[must_use]
    pub fn speed(&self) -> f64 {
        norm(&self.velocity)
    }

    /// Distance to another particle
    ///
    /// # Errors
    /// Returns an error if the particles have different dimensions.
    pub fn distance_to(
        &self,
        other: &Self,
    ) -> Result<f64, String> {
        let r_vec = vec_sub(&self.position, &other.position)?;

        Ok(norm(&r_vec))
    }
}

// ============================================================================
// Interatomic Potentials
// ============================================================================

/// Computes the Lennard-Jones potential and force between two particles.
///
/// The Lennard-Jones potential is a simple mathematical model that describes the
/// interaction between a pair of neutral atoms or molecules. It has a repulsive
/// term at short distances and an attractive term at long distances.
///
/// # Arguments
/// * `p1`, `p2` - The two particles.
/// * `epsilon` - The depth of the potential well.
/// * `sigma` - The finite distance at which the inter-particle potential is zero.
///
/// # Returns
/// A tuple `(potential, force_on_p1)` where `potential` is the scalar potential energy
/// and `force_on_p1` is the force vector acting on `p1` due to `p2`.
///
/// # Errors
/// Returns an error if the particles have different dimensions.
pub fn lennard_jones_interaction(
    p1: &Particle,
    p2: &Particle,
    epsilon: f64,
    sigma: f64,
) -> Result<(f64, Vec<f64>), String> {
    let r_vec = vec_sub(&p1.position, &p2.position)?;

    let r = norm(&r_vec);

    if r < 1e-9 {
        return Ok((f64::INFINITY, vec![0.0; r_vec.len()]));
    }

    let sigma_over_r = sigma / r;

    let sigma_over_r6 = sigma_over_r.powi(6);

    let sigma_over_r12 = sigma_over_r6.powi(2);

    let potential = 4.0 * epsilon * (sigma_over_r12 - sigma_over_r6);

    let force_magnitude = 24.0 * epsilon * 2.0f64.mul_add(sigma_over_r12, -sigma_over_r6) / r;

    let force_on_p1 = scalar_mul(&r_vec, force_magnitude / r);

    Ok((potential, force_on_p1))
}

/// Integrates the equations of motion for a system of particles using the Velocity Verlet algorithm.
///
/// The Velocity Verlet algorithm is a popular numerical integration scheme for molecular dynamics
/// simulations. It is time-reversible and preserves phase space volume, making it suitable
/// for long-term simulations.
///
/// # Arguments
/// * `particles` - A mutable vector of `Particle`s.
/// * `dt` - The time step.
/// * `num_steps` - The number of simulation steps.
/// * `force_calculator` - A closure that computes the total force on each particle.
///   It takes `&mut Vec<Particle>` and returns `Result<(), String>`.
///
/// # Returns
/// A `Vec<Vec<Particle>>` representing the trajectory of the particles over time.
///
/// # Errors
/// Returns an error if the force calculator fails or if vector operations fail due to dimension mismatch.
pub fn integrate_velocity_verlet<F>(
    particles: &mut Vec<Particle>,
    dt: f64,
    num_steps: usize,
    mut force_calculator: F,
) -> Result<Vec<Vec<Particle>>, String>
where
    F: FnMut(&mut Vec<Particle>) -> Result<(), String>,
{
    let mut trajectory = Vec::with_capacity(num_steps + 1);

    trajectory.push(particles.clone());

    force_calculator(particles)?;

    for _step in 0..num_steps {
        for p in particles.iter_mut() {
            let acc = scalar_mul(&p.force, 1.0 / p.mass);

            p.velocity = vec_add(&p.velocity, &scalar_mul(&acc, 0.5 * dt))?;

            p.position = vec_add(&p.position, &scalar_mul(&p.velocity, dt))?;
        }

        force_calculator(particles)?;

        for p in particles.iter_mut() {
            let acc = scalar_mul(&p.force, 1.0 / p.mass);

            p.velocity = vec_add(&p.velocity, &scalar_mul(&acc, 0.5 * dt))?;
        }

        trajectory.push(particles.clone());
    }

    Ok(trajectory)
}

// ============================================================================
// Additional Potentials
// ============================================================================

/// Morse potential for diatomic molecules.
///
/// V(r) = De * (1 - exp(-a(r - re)))²
///
/// # Arguments
/// * `p1`, `p2` - The two particles
/// * `de` - Dissociation energy
/// * `a` - Controls the width of the potential well
/// * `re` - Equilibrium bond distance
///
/// # Errors
/// Returns an error if the particles have different dimensions.
pub fn morse_interaction(
    p1: &Particle,
    p2: &Particle,
    de: f64,
    a: f64,
    re: f64,
) -> Result<(f64, Vec<f64>), String> {
    let r_vec = vec_sub(&p1.position, &p2.position)?;

    let r = norm(&r_vec);

    if r < 1e-9 {
        return Ok((f64::INFINITY, vec![0.0; r_vec.len()]));
    }

    let exp_term = (-a * (r - re)).exp();

    let one_minus_exp = 1.0 - exp_term;

    let potential = de * one_minus_exp * one_minus_exp;

    // V = De (1 - e)^2 with e = exp(-a (r - re)), so dV/dr = 2 De a (1 - e) e
    // and the force along r_vec (from p2 to p1) is -dV/dr.
    let force_magnitude = -2.0 * de * a * one_minus_exp * exp_term;

    let force_on_p1 = scalar_mul(&r_vec, force_magnitude / r);

    Ok((potential, force_on_p1))
}

/// Harmonic (spring) potential.
///
/// V(r) = 0.5 * k * (r - r0)²
///
/// # Arguments
/// * `p1`, `p2` - The two particles
/// * `k` - Spring constant
/// * `r0` - Equilibrium distance
///
/// # Errors
/// Returns an error if the particles have different dimensions.
pub fn harmonic_interaction(
    p1: &Particle,
    p2: &Particle,
    k: f64,
    r0: f64,
) -> Result<(f64, Vec<f64>), String> {
    let r_vec = vec_sub(&p1.position, &p2.position)?;

    let r = norm(&r_vec);

    if r < 1e-9 {
        return Ok((0.0, vec![0.0; r_vec.len()]));
    }

    let dr = r - r0;

    let potential = 0.5 * k * dr * dr;

    let force_magnitude = -k * dr;

    let force_on_p1 = scalar_mul(&r_vec, force_magnitude / r);

    Ok((potential, force_on_p1))
}

/// Coulomb (electrostatic) potential.
///
/// V(r) = `k_e` * q1 * q2 / r
///
/// # Arguments
/// * `p1`, `p2` - The two particles (with charge fields)
/// * `k_coulomb` - Coulomb constant (in appropriate units)
///
/// # Errors
/// Returns an error if the particles have different dimensions.
pub fn coulomb_interaction(
    p1: &Particle,
    p2: &Particle,
    k_coulomb: f64,
) -> Result<(f64, Vec<f64>), String> {
    let r_vec = vec_sub(&p1.position, &p2.position)?;

    let r = norm(&r_vec);

    if r < 1e-9 {
        return Ok((f64::INFINITY, vec![0.0; r_vec.len()]));
    }

    let potential = k_coulomb * p1.charge * p2.charge / r;

    // Force = -dV/dr * r_hat = k * q1 * q2 / r² * r_hat
    let force_magnitude = k_coulomb * p1.charge * p2.charge / (r * r);

    let force_on_p1 = scalar_mul(&r_vec, force_magnitude / r);

    Ok((potential, force_on_p1))
}

/// Soft-sphere potential for avoiding particle overlap.
///
/// V(r) = ε * (σ/r)^n for r < σ, 0 otherwise
///
/// # Errors
/// Returns an error if the particles have different dimensions.
pub fn soft_sphere_interaction(
    p1: &Particle,
    p2: &Particle,
    epsilon: f64,
    sigma: f64,
    n: i32,
) -> Result<(f64, Vec<f64>), String> {
    let r_vec = vec_sub(&p1.position, &p2.position)?;

    let r = norm(&r_vec);

    if r >= sigma {
        return Ok((0.0, vec![0.0; r_vec.len()]));
    }

    if r < 1e-9 {
        return Ok((f64::INFINITY, vec![0.0; r_vec.len()]));
    }

    let sigma_over_r = sigma / r;

    let potential = epsilon * sigma_over_r.powi(n);

    let force_magnitude = f64::from(n) * epsilon * sigma_over_r.powi(n) / r;

    let force_on_p1 = scalar_mul(&r_vec, force_magnitude / r);

    Ok((potential, force_on_p1))
}

// ============================================================================
// System Properties
// ============================================================================

/// Calculates the total kinetic energy of the system.
#[must_use]
pub fn total_kinetic_energy(particles: &[Particle]) -> f64 {
    particles.iter().map(Particle::kinetic_energy).sum()
}

/// Calculates the total momentum of the system.
///
/// # Errors
/// Returns an error if the particle list is empty.
pub fn total_momentum(particles: &[Particle]) -> Result<Vec<f64>, String> {
    if particles.is_empty() {
        return Err("Empty particle \
                    list"
            .to_string());
    }

    let dim = particles[0].position.len();

    let mut total = vec![0.0; dim];

    for p in particles {
        let mom = p.momentum();

        for (i, m) in mom.iter().enumerate() {
            total[i] += m;
        }
    }

    Ok(total)
}

/// Calculates the center of mass of the system.
///
/// # Errors
/// Returns an error if the particle list is empty.
pub fn center_of_mass(particles: &[Particle]) -> Result<Vec<f64>, String> {
    if particles.is_empty() {
        return Err("Empty particle \
                    list"
            .to_string());
    }

    let dim = particles[0].position.len();

    let mut com = vec![0.0; dim];

    let mut total_mass = 0.0;

    for p in particles {
        total_mass += p.mass;

        for (i, pos) in p.position.iter().enumerate() {
            com[i] += p.mass * pos;
        }
    }

    for c in &mut com {
        *c /= total_mass;
    }

    Ok(com)
}

/// Calculates the temperature from kinetic energy.
///
/// T = 2 * KE / (dim * N * `k_B`)
/// In reduced units with `k_B` = 1: T = 2 * KE / (dim * N)
#[must_use]
pub fn temperature(particles: &[Particle]) -> f64 {
    if particles.is_empty() {
        return 0.0;
    }

    let dim = particles[0].position.len();

    let ke = total_kinetic_energy(particles);

    let n = particles.len();

    // In reduced units (k_B = 1)
    2.0 * ke / (dim * n) as f64
}

/// Calculates instantaneous pressure using the virial theorem.
///
/// P = (N * `k_B` * T + virial) / V
#[must_use]
pub fn pressure(
    particles: &[Particle],
    volume: f64,
    virial: f64,
) -> f64 {
    if particles.is_empty() || volume <= 0.0 {
        return 0.0;
    }

    let n = particles.len() as f64;

    let t = temperature(particles);

    // In reduced units (k_B = 1)
    n.mul_add(t, virial) / volume
}

/// Removes center of mass velocity from the system.
///
/// # Errors
/// Returns an error if the particle list is empty.
pub fn remove_com_velocity(particles: &mut [Particle]) -> Result<(), String> {
    if particles.is_empty() {
        return Ok(());
    }

    let dim = particles[0].position.len();

    let mut total_momentum = vec![0.0; dim];

    let mut total_mass = 0.0;

    for p in particles.iter() {
        total_mass += p.mass;

        for (i, v) in p.velocity.iter().enumerate() {
            total_momentum[i] += p.mass * v;
        }
    }

    let com_velocity: Vec<f64> = total_momentum.iter().map(|m| m / total_mass).collect();

    for p in particles.iter_mut() {
        for (i, v) in p.velocity.iter_mut().enumerate() {
            *v -= com_velocity[i];
        }
    }

    Ok(())
}

// ============================================================================
// Thermostats
// ============================================================================

/// Velocity rescaling thermostat.
///
/// Rescales velocities to achieve target temperature.
pub fn velocity_rescale(
    particles: &mut [Particle],
    target_temp: f64,
) {
    let current_temp = temperature(particles);

    if current_temp <= 0.0 {
        return;
    }

    let scale = (target_temp / current_temp).sqrt();

    for p in particles.iter_mut() {
        for v in &mut p.velocity {
            *v *= scale;
        }
    }
}

/// Berendsen thermostat.
///
/// Gently couples the system to a heat bath.
///
/// # Arguments
/// * `particles` - System particles
/// * `target_temp` - Target temperature
/// * `tau` - Coupling time constant
/// * `dt` - Time step
pub fn berendsen_thermostat(
    particles: &mut [Particle],
    target_temp: f64,
    tau: f64,
    dt: f64,
) {
    let current_temp = temperature(particles);

    if current_temp <= 0.0 {
        return;
    }

    let scale = (dt / tau)
        .mul_add(target_temp / current_temp - 1.0, 1.0)
        .sqrt();

    for p in particles.iter_mut() {
        for v in &mut p.velocity {
            *v *= scale;
        }
    }
}

// ============================================================================
// Periodic Boundary Conditions
// ============================================================================

/// Applies periodic boundary conditions to a position.
#[must_use]
pub fn apply_pbc(
    position: &[f64],
    box_size: &[f64],
) -> Vec<f64> {
    position
        .iter()
        .zip(box_size.iter())
        .map(|(&x, &l)| {
            let mut wrapped = x % l;

            if wrapped < 0.0 {
                wrapped += l;
            }

            wrapped
        })
        .collect()
}

/// Applies minimum image convention for distance calculation.
#[must_use]
pub fn minimum_image_distance(
    r: &[f64],
    box_size: &[f64],
) -> Vec<f64> {
    r.iter()
        .zip(box_size.iter())
        .map(|(&dx, &l)| {
            let mut d = dx % l;

            if d > l / 2.0 {
                d -= l;
            } else if d < -l / 2.0 {
                d += l;
            }

            d
        })
        .collect()
}

// ============================================================================
// Analysis Functions
// ============================================================================

/// Calculates the radial distribution function g(r).
///
/// # Arguments
/// * `particles` - System particles
/// * `box_size` - Box dimensions
/// * `num_bins` - Number of histogram bins
/// * `r_max` - Maximum distance to consider
///
/// # Returns
/// (`r_values`, `g_r`) where `r_values` are bin centers and `g_r` are g(r) values
#[must_use]
pub fn radial_distribution_function(
    particles: &[Particle],
    box_size: &[f64],
    num_bins: usize,
    r_max: f64,
) -> (Vec<f64>, Vec<f64>) {
    let n = particles.len();

    if n < 2 {
        return (vec![], vec![]);
    }

    let dr = r_max / num_bins as f64;

    let mut histogram = vec![0usize; num_bins];

    // Count pairs in each bin
    for i in 0..n {
        for j in (i + 1)..n {
            if let Ok(r_vec) = vec_sub(&particles[i].position, &particles[j].position) {
                let r_mic = minimum_image_distance(&r_vec, box_size);

                let r = norm(&r_mic);

                if r < r_max {
                    // Safe to cast as r and dr are positive
                    let bin: usize = ((r / dr) as i64).try_into().unwrap_or(0);

                    if bin < num_bins {
                        histogram[bin] += 2; // Count both i-j and j-i
                    }
                }
            }
        }
    }

    // Normalize by ideal gas distribution
    let volume: f64 = box_size.iter().product();

    let rho = n as f64 / volume;

    let pi = std::f64::consts::PI;

    let r_values: Vec<f64> = (0..num_bins).map(|i| (i as f64 + 0.5) * dr).collect();

    let g_r: Vec<f64> = histogram
        .iter()
        .enumerate()
        .map(|(i, &count)| {
            let _r = (i as f64 + 0.5) * dr;

            let shell_volume =
                (4.0 / 3.0) * pi * (((i + 1) as f64 * dr).powi(3) - (i as f64 * dr).powi(3));

            let ideal_count = rho * shell_volume * n as f64;

            if ideal_count > 0.0 {
                count as f64 / ideal_count
            } else {
                0.0
            }
        })
        .collect();

    (r_values, g_r)
}

/// Calculates mean square displacement.
///
/// MSD(t) = <|r(t) - r(0)|²>
#[must_use]
pub fn mean_square_displacement(
    initial: &[Particle],
    current: &[Particle],
) -> f64 {
    if initial.len() != current.len() || initial.is_empty() {
        return 0.0;
    }

    let mut msd = 0.0;

    for (p0, p) in initial.iter().zip(current.iter()) {
        if let Ok(dr) = vec_sub(&p.position, &p0.position) {
            let dr2: f64 = dr.iter().map(|x| x * x).sum();

            msd += dr2;
        }
    }

    msd / initial.len() as f64
}

/// Initializes particle velocities from Maxwell-Boltzmann distribution.
pub fn initialize_velocities_maxwell_boltzmann(
    particles: &mut [Particle],
    target_temp: f64,
    rng_seed: u64,
) {
    // Simple pseudo-random generator (LCG)
    let mut rng_state = rng_seed;

    let _next_random = || {
        rng_state = rng_state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);

        (rng_state >> 33) as f64 / (1u64 << 31) as f64
    };

    // Box-Muller transform for Gaussian random numbers
    let gaussian = |rng: &mut dyn FnMut() -> f64| -> f64 {
        let u1 = rng();

        let u2 = rng();

        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    };

    for p in particles.iter_mut() {
        let sigma = (target_temp / p.mass).sqrt();

        for v in &mut p.velocity {
            // Create a mutable closure to use with gaussian
            let mut rng_fn = || {
                rng_state = rng_state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);

                (rng_state >> 33) as f64 / (1u64 << 31) as f64
            };

            *v = sigma * gaussian(&mut rng_fn);
        }
    }

    // Remove center of mass velocity
    let _ = remove_com_velocity(particles);

    // Rescale to exact target temperature
    velocity_rescale(particles, target_temp);
}

/// Creates a simple cubic lattice of particles.
#[must_use]
pub fn create_cubic_lattice(
    n_per_side: usize,
    lattice_constant: f64,
    mass: f64,
) -> Vec<Particle> {
    let mut particles = Vec::with_capacity(n_per_side * n_per_side * n_per_side);

    let mut id = 0;

    for i in 0..n_per_side {
        for j in 0..n_per_side {
            for k in 0..n_per_side {
                let position = vec![
                    i as f64 * lattice_constant,
                    j as f64 * lattice_constant,
                    k as f64 * lattice_constant,
                ];

                let velocity = vec![0.0, 0.0, 0.0];

                particles.push(Particle::new(id, mass, position, velocity));

                id += 1;
            }
        }
    }

    particles
}

/// Creates an FCC (face-centered cubic) lattice of particles.
#[must_use]
pub fn create_fcc_lattice(
    n_cells: usize,
    lattice_constant: f64,
    mass: f64,
) -> Vec<Particle> {
    let mut particles = Vec::with_capacity(4 * n_cells * n_cells * n_cells);

    let mut id = 0;

    // FCC basis positions (in units of lattice constant)
    let basis = [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5],
    ];

    for i in 0..n_cells {
        for j in 0..n_cells {
            for k in 0..n_cells {
                for b in &basis {
                    let position = vec![
                        (i as f64 + b[0]) * lattice_constant,
                        (j as f64 + b[1]) * lattice_constant,
                        (k as f64 + b[2]) * lattice_constant,
                    ];

                    let velocity = vec![0.0, 0.0, 0.0];

                    particles.push(Particle::new(id, mass, position, velocity));

                    id += 1;
                }
            }
        }
    }

    particles
}

// ============================================================================
// Stochastic and extended-system thermostats, neighbour lists
// ============================================================================

/// A deterministic Gaussian sampler (`SplitMix64` + Box–Muller), so that
/// Langevin runs are reproducible from a seed.
#[derive(Clone, Debug)]
pub struct GaussianStream {
    state: u64,
    spare: Option<f64>,
}

impl GaussianStream {
    /// A stream seeded with `seed`.
    #[must_use]
    pub const fn new(seed: u64) -> Self {
        Self { state: seed, spare: None }
    }

    fn uniform(&mut self) -> f64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        ((z >> 11) as f64 + 0.5) / (1_u64 << 53) as f64
    }

    /// A standard normal deviate.
    pub fn next_normal(&mut self) -> f64 {
        if let Some(s) = self.spare.take() {
            return s;
        }
        let (u1, u2) = (self.uniform(), self.uniform());
        let r = (-2.0 * u1.ln()).sqrt();
        let theta = 2.0 * std::f64::consts::PI * u2;
        self.spare = Some(r * theta.sin());
        r * theta.cos()
    }
}

/// Langevin dynamics with the BAOAB splitting (Leimkuhler–Matthews):
///
/// half kick (B), half drift (A), exact Ornstein–Uhlenbeck velocity update
/// with friction `gamma` at temperature `kt` (O), half drift, half kick.
/// Samples the canonical ensemble with configurational error `O(dt²)`.
///
/// # Errors
/// Propagates errors of the force calculation.
pub fn integrate_langevin_baoab<F>(
    particles: &mut [Particle],
    dt: f64,
    num_steps: usize,
    gamma: f64,
    kt: f64,
    noise: &mut GaussianStream,
    mut force_calculator: F,
) -> Result<(), String>
where
    F: FnMut(&mut [Particle]) -> Result<(), String>,
{
    let c1 = (-gamma * dt).exp();
    force_calculator(particles)?;
    for _ in 0..num_steps {
        for p in particles.iter_mut() {
            for k in 0..p.velocity.len() {
                p.velocity[k] += 0.5 * dt * p.force[k] / p.mass;
                p.position[k] += 0.5 * dt * p.velocity[k];
            }
            let c2 = ((1.0 - c1 * c1) * kt / p.mass).sqrt();
            for v in &mut p.velocity {
                *v = c1 * *v + c2 * noise.next_normal();
            }
            for k in 0..p.position.len() {
                p.position[k] += 0.5 * dt * p.velocity[k];
            }
        }
        force_calculator(particles)?;
        for p in particles.iter_mut() {
            for k in 0..p.velocity.len() {
                p.velocity[k] += 0.5 * dt * p.force[k] / p.mass;
            }
        }
    }
    Ok(())
}

/// The Nosé–Hoover thermostat:
///
/// the friction `ξ` obeys `Q dξ/dt = Σ m v² - N_f kT`, a deterministic
/// extended system whose invariant `H + Q ξ²/2 + N_f kT s` (with `ds/dt =
/// ξ`) is conserved. Velocity-Verlet with the friction applied as an exact
/// exponential scaling in half steps. Returns the final `(ξ, s)`.
///
/// # Errors
/// Propagates errors of the force calculation.
#[allow(clippy::too_many_arguments)]
pub fn integrate_nose_hoover<F>(
    particles: &mut [Particle],
    dt: f64,
    num_steps: usize,
    kt: f64,
    q_mass: f64,
    xi0: f64,
    mut force_calculator: F,
) -> Result<(f64, f64), String>
where
    F: FnMut(&mut [Particle]) -> Result<(), String>,
{
    let degrees: f64 = particles.iter().map(|p| p.velocity.len() as f64).sum();
    let kinetic2 = |ps: &[Particle]| -> f64 { ps.iter().map(|p| p.mass * p.velocity.iter().map(|v| v * v).sum::<f64>()).sum() };
    let (mut xi, mut s) = (xi0, 0.0);
    force_calculator(particles)?;
    for _ in 0..num_steps {
        // Thermostat half step.
        xi += 0.25 * dt * (kinetic2(particles) - degrees * kt) / q_mass;
        let scale = (-0.5 * dt * xi).exp();
        for p in particles.iter_mut() {
            p.velocity.iter_mut().for_each(|v| *v *= scale);
        }
        s += 0.5 * dt * xi;
        xi += 0.25 * dt * (kinetic2(particles) - degrees * kt) / q_mass;
        // Velocity Verlet.
        for p in particles.iter_mut() {
            for k in 0..p.velocity.len() {
                p.velocity[k] += 0.5 * dt * p.force[k] / p.mass;
                p.position[k] += dt * p.velocity[k];
            }
        }
        force_calculator(particles)?;
        for p in particles.iter_mut() {
            for k in 0..p.velocity.len() {
                p.velocity[k] += 0.5 * dt * p.force[k] / p.mass;
            }
        }
        // Thermostat half step.
        xi += 0.25 * dt * (kinetic2(particles) - degrees * kt) / q_mass;
        let scale = (-0.5 * dt * xi).exp();
        for p in particles.iter_mut() {
            p.velocity.iter_mut().for_each(|v| *v *= scale);
        }
        s += 0.5 * dt * xi;
        xi += 0.25 * dt * (kinetic2(particles) - degrees * kt) / q_mass;
    }
    Ok((xi, s))
}

/// Pairs `(i, j)`, `i < j`, closer than `cutoff` under periodic
/// boundaries in a box of side `box_length`.
///
/// They are found with a linked-cell list in `O(N)` time (cells of side at
/// least `cutoff`; a plain double loop when the box holds fewer than three
/// cells per side).
#[must_use]
#[allow(clippy::cast_sign_loss, clippy::cast_possible_truncation, clippy::cast_possible_wrap, clippy::cast_precision_loss)]
pub fn neighbour_pairs(
    particles: &[Particle],
    box_length: f64,
    cutoff: f64,
) -> Vec<(usize, usize)> {
    let dim = particles.first().map_or(0, |p| p.position.len());
    let distance2 = |a: &[f64], b: &[f64]| -> f64 {
        a.iter()
            .zip(b)
            .map(|(x, y)| {
                let mut d = x - y;
                d -= box_length * (d / box_length).round();
                d * d
            })
            .sum()
    };
    let cells_per_side = (box_length / cutoff).floor() as usize;
    let mut pairs = Vec::new();
    if cells_per_side < 3 || dim == 0 || dim > 3 {
        for i in 0..particles.len() {
            for j in i + 1..particles.len() {
                if distance2(&particles[i].position, &particles[j].position) < cutoff * cutoff {
                    pairs.push((i, j));
                }
            }
        }
        return pairs;
    }
    let side = box_length / cells_per_side as f64;
    let cell_of = |pos: &[f64]| -> Vec<usize> {
        pos.iter()
            .map(|&x| {
                let wrapped = x.rem_euclid(box_length);
                ((wrapped / side) as usize).min(cells_per_side - 1)
            })
            .collect()
    };
    let flat = |c: &[usize]| c.iter().fold(0, |acc, &k| acc * cells_per_side + k);
    let mut cells: Vec<Vec<usize>> = vec![Vec::new(); cells_per_side.pow(dim as u32)];
    for (i, p) in particles.iter().enumerate() {
        cells[flat(&cell_of(&p.position))].push(i);
    }
    // Neighbouring cell offsets in {-1, 0, 1}^dim.
    let offsets: Vec<Vec<isize>> = (0..3_usize.pow(dim as u32))
        .map(|mut code| {
            (0..dim)
                .map(|_| {
                    let o = (code % 3) as isize - 1;
                    code /= 3;
                    o
                })
                .collect()
        })
        .collect();
    for (i, p) in particles.iter().enumerate() {
        let home = cell_of(&p.position);
        for offset in &offsets {
            let neighbour: Vec<usize> =
                home.iter().zip(offset).map(|(&c, &o)| (c as isize + o).rem_euclid(cells_per_side as isize) as usize).collect();
            for &j in &cells[flat(&neighbour)] {
                if j > i && distance2(&p.position, &particles[j].position) < cutoff * cutoff {
                    pairs.push((i, j));
                }
            }
        }
    }
    pairs.sort_unstable();
    pairs.dedup();
    pairs
}

#[cfg(test)]
mod thermostat_tests {
    use super::*;

    #[allow(clippy::unnecessary_wraps)] // the force-callback signature
    fn harmonic(ps: &mut [Particle]) -> Result<(), String> {
        for p in ps.iter_mut() {
            p.force = p.position.iter().map(|x| -x).collect();
        }
        Ok(())
    }

    #[test]
    fn langevin_samples_the_canonical_distribution() {
        // Harmonic oscillators: <v²> = <x²> = kT (m = k = 1).
        let mut ps: Vec<Particle> = (0..200).map(|i| Particle::new(i, 1.0, vec![0.0], vec![0.0])).collect();
        let mut noise = GaussianStream::new(42);
        integrate_langevin_baoab(&mut ps, 0.05, 2000, 1.0, 0.7, &mut noise, harmonic).unwrap();
        let (mut v2, mut x2) = (0.0, 0.0);
        let mut samples = 0.0;
        for _ in 0..200 {
            integrate_langevin_baoab(&mut ps, 0.05, 10, 1.0, 0.7, &mut noise, harmonic).unwrap();
            for p in &ps {
                v2 += p.velocity[0] * p.velocity[0];
                x2 += p.position[0] * p.position[0];
                samples += 1.0;
            }
        }
        assert!((v2 / samples - 0.7).abs() < 0.03 && (x2 / samples - 0.7).abs() < 0.03, "{} {}", v2 / samples, x2 / samples);
    }

    #[test]
    fn nose_hoover_conserves_its_extended_energy() {
        let mut ps: Vec<Particle> =
            (0..5).map(|i| Particle::new(i, 1.0, vec![0.1 * i as f64, -0.2], vec![0.3, 0.1 * i as f64 - 0.2])).collect();
        let (kt, q) = (0.5, 2.0);
        let energy = |ps: &[Particle]| -> f64 {
            ps.iter().map(|p| 0.5 * p.mass * p.velocity.iter().map(|v| v * v).sum::<f64>() + 0.5 * p.position.iter().map(|x| x * x).sum::<f64>()).sum()
        };
        let e0 = energy(&ps);
        let (xi, s) = integrate_nose_hoover(&mut ps, 0.01, 5000, kt, q, 0.0, harmonic).unwrap();
        let degrees = 10.0;
        let invariant = energy(&ps) + 0.5 * q * xi * xi + degrees * kt * s;
        assert!((invariant - e0).abs() < 1e-3 * (1.0 + e0), "{invariant} vs {e0}");
    }

    #[test]
    fn cell_lists_find_every_close_pair() {
        let mut noise = GaussianStream::new(3);
        let box_length = 10.0;
        let ps: Vec<Particle> = (0..300)
            .map(|i| {
                let pos: Vec<f64> = (0..3).map(|_| (noise.next_normal() * 3.0).rem_euclid(box_length)).collect();
                Particle::new(i, 1.0, pos, vec![0.0; 3])
            })
            .collect();
        let fast = neighbour_pairs(&ps, box_length, 2.5);
        let slow = neighbour_pairs(&ps, box_length, 6.0); // brute force path
        let brute: Vec<(usize, usize)> = slow
            .into_iter()
            .filter(|&(i, j)| {
                let d2: f64 = ps[i]
                    .position
                    .iter()
                    .zip(&ps[j].position)
                    .map(|(a, b)| {
                        let mut d = a - b;
                        d -= box_length * (d / box_length).round();
                        d * d
                    })
                    .sum();
                d2 < 2.5 * 2.5
            })
            .collect();
        assert_eq!(fast, brute);
        assert!(!fast.is_empty());
    }
}
