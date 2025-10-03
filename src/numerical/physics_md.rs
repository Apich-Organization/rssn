//! # Numerical Molecular Dynamics (MD)
//!
//! This module provides numerical methods for Molecular Dynamics (MD) simulations.
//! It includes implementations of common integration algorithms like Velocity Verlet
//! for simulating the motion of particles under interatomic forces.

use crate::numerical::vector::{norm, scalar_mul, vec_add, vec_sub};

/// Represents a particle in a molecular dynamics simulation.
#[derive(Clone)]
pub struct Particle {
    pub id: usize,
    pub mass: f64,
    pub position: Vec<f64>,
    pub velocity: Vec<f64>,
    pub force: Vec<f64>,
}

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
pub fn lennard_jones_interaction(
    p1: &Particle,
    p2: &Particle,
    epsilon: f64,
    sigma: f64,
) -> (f64, Vec<f64>) {
    let r_vec = vec_sub(&p1.position, &p2.position).unwrap();
    let r = norm(&r_vec);

    if r < 1e-9 {
        // Avoid division by zero for overlapping particles
        return (f64::INFINITY, vec![0.0; r_vec.len()]);
    }

    let sigma_over_r = sigma / r;
    let sigma_over_r6 = sigma_over_r.powi(6);
    let sigma_over_r12 = sigma_over_r6.powi(2);

    let potential = 4.0 * epsilon * (sigma_over_r12 - sigma_over_r6);

    let force_magnitude = 24.0 * epsilon * (2.0 * sigma_over_r12 - sigma_over_r6) / r;
    let force_on_p1 = scalar_mul(&r_vec, force_magnitude / r);

    (potential, force_on_p1)
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

    // Initial forces
    force_calculator(particles)?;

    for _step in 0..num_steps {
        for p in particles.iter_mut() {
            // Update velocity (half step)
            let acc = scalar_mul(&p.force, 1.0 / p.mass);
            p.velocity = vec_add(&p.velocity, &scalar_mul(&acc, 0.5 * dt)).unwrap();

            // Update position
            p.position = vec_add(&p.position, &scalar_mul(&p.velocity, dt)).unwrap();
        }

        // Compute new forces at new positions
        force_calculator(particles)?;

        for p in particles.iter_mut() {
            // Update velocity (full step)
            let acc = scalar_mul(&p.force, 1.0 / p.mass);
            p.velocity = vec_add(&p.velocity, &scalar_mul(&acc, 0.5 * dt)).unwrap();
        }
        trajectory.push(particles.clone());
    }

    Ok(trajectory)
}
