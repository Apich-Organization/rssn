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
//! 
//! use rssn::numerical::physics_md::*;
//!
//! let p1 = Particle::new(
//!     0,
//!     1.0,
//!     vec![0.0, 0.0, 0.0],
//!     vec![0.0, 0.0, 0.0],
//! );
//!
//! let p2 = Particle::new(
//!     1,
//!     1.0,
//!     vec![1.5, 0.0, 0.0],
//!     vec![0.0, 0.0, 0.0],
//! );
//!
//! let (potential, force) =
//!     lennard_jones_interaction(
//!         &p1, &p2, 1.0, 1.0,
//!     )
//!     .unwrap();
//! ```

use serde::Deserialize;
use serde::Serialize;

use crate::numerical::vector::norm;
use crate::numerical::vector::scalar_mul;
use crate::numerical::vector::vec_add;
use crate::numerical::vector::vec_sub;

// ============================================================================
// Physical Constants (Reduced Units)
// ============================================================================

/// Boltzmann constant in SI units (J/K)

pub const BOLTZMANN_CONSTANT_SI: f64 =
    1.380_649e-23;

/// Avogadro's number (1/mol)

pub const AVOGADRO_NUMBER: f64 =
    6.022_140_76e23;

/// Reduced unit for temperature (using argon as reference)
/// 1 reduced temperature = `ε/k_B` ≈ 120 K for argon

pub const TEMPERATURE_UNIT_ARGON: f64 =
    119.8;

/// Reduced unit for length (using argon as reference)
/// 1 reduced length = σ ≈ 3.4 Å for argon

pub const LENGTH_UNIT_ARGON: f64 =
    3.4e-10;

/// Reduced unit for energy (using argon as reference)
/// 1 reduced energy = ε ≈ 1.65e-21 J for argon

pub const ENERGY_UNIT_ARGON: f64 =
    1.65e-21;

// ============================================================================
// Particle Definition
// ============================================================================

/// Represents a particle in a molecular dynamics simulation.
/// cbindgen:ignore
#[derive(
    Clone, Debug, Serialize, Deserialize,
)]

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

    pub fn kinetic_energy(
        &self
    ) -> f64 {

        let v2: f64 = self
            .velocity
            .iter()
            .map(|v| v * v)
            .sum();

        0.5 * self.mass * v2
    }

    /// Returns the momentum of the particle: p = m * v
    #[must_use]

    pub fn momentum(&self) -> Vec<f64> {

        scalar_mul(
            &self.velocity,
            self.mass,
        )
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

        let r_vec = vec_sub(
            &self.position,
            &other.position,
        )?;

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

    let r_vec = vec_sub(
        &p1.position,
        &p2.position,
    )?;

    let r = norm(&r_vec);

    if r < 1e-9 {

        return Ok((
            f64::INFINITY,
            vec![0.0; r_vec.len()],
        ));
    }

    let sigma_over_r = sigma / r;

    let sigma_over_r6 =
        sigma_over_r.powi(6);

    let sigma_over_r12 =
        sigma_over_r6.powi(2);

    let potential = 4.0
        * epsilon
        * (sigma_over_r12
            - sigma_over_r6);

    let force_magnitude = 24.0
        * epsilon
        * 2.0f64.mul_add(
            sigma_over_r12,
            -sigma_over_r6,
        )
        / r;

    let force_on_p1 = scalar_mul(
        &r_vec,
        force_magnitude / r,
    );

    Ok((
        potential,
        force_on_p1,
    ))
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
    F: FnMut(
        &mut Vec<Particle>,
    ) -> Result<(), String>,
{

    let mut trajectory =
        Vec::with_capacity(
            num_steps + 1,
        );

    trajectory.push(particles.clone());

    force_calculator(particles)?;

    for _step in 0 .. num_steps {

        for p in particles.iter_mut() {

            let acc = scalar_mul(
                &p.force,
                1.0 / p.mass,
            );

            p.velocity = vec_add(
                &p.velocity,
                &scalar_mul(
                    &acc,
                    0.5 * dt,
                ),
            )?;

            p.position = vec_add(
                &p.position,
                &scalar_mul(
                    &p.velocity,
                    dt,
                ),
            )?;
        }

        force_calculator(particles)?;

        for p in particles.iter_mut() {

            let acc = scalar_mul(
                &p.force,
                1.0 / p.mass,
            );

            p.velocity = vec_add(
                &p.velocity,
                &scalar_mul(
                    &acc,
                    0.5 * dt,
                ),
            )?;
        }

        trajectory
            .push(particles.clone());
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

    let r_vec = vec_sub(
        &p1.position,
        &p2.position,
    )?;

    let r = norm(&r_vec);

    if r < 1e-9 {

        return Ok((
            f64::INFINITY,
            vec![0.0; r_vec.len()],
        ));
    }

    let exp_term =
        (-a * (r - re)).exp();

    let one_minus_exp = 1.0 - exp_term;

    let potential = de
        * one_minus_exp
        * one_minus_exp;

    // Force = -dV/dr = 2 * De * a * (1 - exp) * exp
    let force_magnitude = 2.0
        * de
        * a
        * one_minus_exp
        * exp_term;

    let force_on_p1 = scalar_mul(
        &r_vec,
        force_magnitude / r,
    );

    Ok((
        potential,
        force_on_p1,
    ))
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

    let r_vec = vec_sub(
        &p1.position,
        &p2.position,
    )?;

    let r = norm(&r_vec);

    if r < 1e-9 {

        return Ok((
            0.0,
            vec![0.0; r_vec.len()],
        ));
    }

    let dr = r - r0;

    let potential = 0.5 * k * dr * dr;

    let force_magnitude = -k * dr;

    let force_on_p1 = scalar_mul(
        &r_vec,
        force_magnitude / r,
    );

    Ok((
        potential,
        force_on_p1,
    ))
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

    let r_vec = vec_sub(
        &p1.position,
        &p2.position,
    )?;

    let r = norm(&r_vec);

    if r < 1e-9 {

        return Ok((
            f64::INFINITY,
            vec![0.0; r_vec.len()],
        ));
    }

    let potential = k_coulomb
        * p1.charge
        * p2.charge
        / r;

    // Force = -dV/dr * r_hat = k * q1 * q2 / r² * r_hat
    let force_magnitude = k_coulomb
        * p1.charge
        * p2.charge
        / (r * r);

    let force_on_p1 = scalar_mul(
        &r_vec,
        force_magnitude / r,
    );

    Ok((
        potential,
        force_on_p1,
    ))
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

    let r_vec = vec_sub(
        &p1.position,
        &p2.position,
    )?;

    let r = norm(&r_vec);

    if r >= sigma {

        return Ok((
            0.0,
            vec![0.0; r_vec.len()],
        ));
    }

    if r < 1e-9 {

        return Ok((
            f64::INFINITY,
            vec![0.0; r_vec.len()],
        ));
    }

    let sigma_over_r = sigma / r;

    let potential =
        epsilon * sigma_over_r.powi(n);

    let force_magnitude = f64::from(n)
        * epsilon
        * sigma_over_r.powi(n)
        / r;

    let force_on_p1 = scalar_mul(
        &r_vec,
        force_magnitude / r,
    );

    Ok((
        potential,
        force_on_p1,
    ))
}

// ============================================================================
// System Properties
// ============================================================================

/// Calculates the total kinetic energy of the system.
#[must_use]

pub fn total_kinetic_energy(
    particles: &[Particle]
) -> f64 {

    particles
        .iter()
        .map(Particle::kinetic_energy)
        .sum()
}

/// Calculates the total momentum of the system.
///
/// # Errors
/// Returns an error if the particle list is empty.

pub fn total_momentum(
    particles: &[Particle]
) -> Result<Vec<f64>, String> {

    if particles.is_empty() {

        return Err("Empty particle \
                    list"
            .to_string());
    }

    let dim = particles[0]
        .position
        .len();

    let mut total = vec![0.0; dim];

    for p in particles {

        let mom = p.momentum();

        for (i, m) in mom
            .iter()
            .enumerate()
        {

            total[i] += m;
        }
    }

    Ok(total)
}

/// Calculates the center of mass of the system.
///
/// # Errors
/// Returns an error if the particle list is empty.

pub fn center_of_mass(
    particles: &[Particle]
) -> Result<Vec<f64>, String> {

    if particles.is_empty() {

        return Err("Empty particle \
                    list"
            .to_string());
    }

    let dim = particles[0]
        .position
        .len();

    let mut com = vec![0.0; dim];

    let mut total_mass = 0.0;

    for p in particles {

        total_mass += p.mass;

        for (i, pos) in p
            .position
            .iter()
            .enumerate()
        {

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

pub fn temperature(
    particles: &[Particle]
) -> f64 {

    if particles.is_empty() {

        return 0.0;
    }

    let dim = particles[0]
        .position
        .len();

    let ke =
        total_kinetic_energy(particles);

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

    if particles.is_empty()
        || volume <= 0.0
    {

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

pub fn remove_com_velocity(
    particles: &mut [Particle]
) -> Result<(), String> {

    if particles.is_empty() {

        return Ok(());
    }

    let dim = particles[0]
        .position
        .len();

    let mut total_momentum =
        vec![0.0; dim];

    let mut total_mass = 0.0;

    for p in particles.iter() {

        total_mass += p.mass;

        for (i, v) in p
            .velocity
            .iter()
            .enumerate()
        {

            total_momentum[i] +=
                p.mass * v;
        }
    }

    let com_velocity: Vec<f64> =
        total_momentum
            .iter()
            .map(|m| m / total_mass)
            .collect();

    for p in particles.iter_mut() {

        for (i, v) in p
            .velocity
            .iter_mut()
            .enumerate()
        {

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

    let current_temp =
        temperature(particles);

    if current_temp <= 0.0 {

        return;
    }

    let scale = (target_temp
        / current_temp)
        .sqrt();

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

    let current_temp =
        temperature(particles);

    if current_temp <= 0.0 {

        return;
    }

    let scale = (dt / tau)
        .mul_add(
            target_temp / current_temp
                - 1.0,
            1.0,
        )
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

    let mut histogram =
        vec![0usize; num_bins];

    // Count pairs in each bin
    for i in 0 .. n {

        for j in (i + 1) .. n {

            if let Ok(r_vec) = vec_sub(
                &particles[i].position,
                &particles[j].position,
            ) {

                let r_mic = minimum_image_distance(&r_vec, box_size);

                let r = norm(&r_mic);

                if r < r_max {
                    // Safe to cast as r and dr are positive
                    let bin = ((r / dr) as isize) as usize;

                    if bin < num_bins {

                        histogram
                            [bin] += 2; // Count both i-j and j-i
                    }
                }
            }
        }
    }

    // Normalize by ideal gas distribution
    let volume: f64 = box_size
        .iter()
        .product();

    let rho = n as f64 / volume;

    let pi = std::f64::consts::PI;

    let r_values: Vec<f64> = (0
        .. num_bins)
        .map(|i| (i as f64 + 0.5) * dr)
        .collect();

    let g_r: Vec<f64> = histogram
        .iter()
        .enumerate()
        .map(|(i, &count)| {

            let _r =
                (i as f64 + 0.5) * dr;

            let shell_volume = (4.0
                / 3.0)
                * pi
                * (((i + 1) as f64
                    * dr)
                    .powi(3)
                    - (i as f64 * dr)
                        .powi(3));

            let ideal_count = rho
                * shell_volume
                * n as f64;

            if ideal_count > 0.0 {

                count as f64
                    / ideal_count
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

    if initial.len() != current.len()
        || initial.is_empty()
    {

        return 0.0;
    }

    let mut msd = 0.0;

    for (p0, p) in initial
        .iter()
        .zip(current.iter())
    {

        if let Ok(dr) = vec_sub(
            &p.position,
            &p0.position,
        ) {

            let dr2: f64 = dr
                .iter()
                .map(|x| x * x)
                .sum();

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
            .wrapping_mul(
                6_364_136_223_846_793_005,
            )
            .wrapping_add(
                1_442_695_040_888_963_407,
            );

        (rng_state >> 33) as f64
            / (1u64 << 31) as f64
    };

    // Box-Muller transform for Gaussian random numbers
    let gaussian =
        |rng: &mut dyn FnMut()
            -> f64|
         -> f64 {

            let u1 = rng();

            let u2 = rng();

            (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
        };

    for p in particles.iter_mut() {

        let sigma = (target_temp
            / p.mass)
            .sqrt();

        for v in &mut p.velocity {

            // Create a mutable closure to use with gaussian
            let mut rng_fn = || {

                rng_state = rng_state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);

                (rng_state >> 33) as f64
                    / (1u64 << 31)
                        as f64
            };

            *v = sigma
                * gaussian(&mut rng_fn);
        }
    }

    // Remove center of mass velocity
    let _ =
        remove_com_velocity(particles);

    // Rescale to exact target temperature
    velocity_rescale(
        particles,
        target_temp,
    );
}

/// Creates a simple cubic lattice of particles.
#[must_use]

pub fn create_cubic_lattice(
    n_per_side: usize,
    lattice_constant: f64,
    mass: f64,
) -> Vec<Particle> {

    let mut particles =
        Vec::with_capacity(
            n_per_side
                * n_per_side
                * n_per_side,
        );

    let mut id = 0;

    for i in 0 .. n_per_side {

        for j in 0 .. n_per_side {

            for k in 0 .. n_per_side {

                let position = vec![
                    i as f64 * lattice_constant,
                    j as f64 * lattice_constant,
                    k as f64 * lattice_constant,
                ];

                let velocity =
                    vec![0.0, 0.0, 0.0];

                particles.push(
                    Particle::new(
                        id,
                        mass,
                        position,
                        velocity,
                    ),
                );

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

    let mut particles =
        Vec::with_capacity(
            4 * n_cells
                * n_cells
                * n_cells,
        );

    let mut id = 0;

    // FCC basis positions (in units of lattice constant)
    let basis = [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5],
    ];

    for i in 0 .. n_cells {

        for j in 0 .. n_cells {

            for k in 0 .. n_cells {

                for b in &basis {

                    let position = vec![
                        (i as f64 + b[0]) * lattice_constant,
                        (j as f64 + b[1]) * lattice_constant,
                        (k as f64 + b[2]) * lattice_constant,
                    ];

                    let velocity = vec![
                        0.0, 0.0, 0.0,
                    ];

                    particles.push(
                        Particle::new(
                            id,
                            mass,
                            position,
                            velocity,
                        ),
                    );

                    id += 1;
                }
            }
        }
    }

    particles
}
