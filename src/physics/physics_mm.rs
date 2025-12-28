use std::ops::Add;
use std::ops::Div;
use std::ops::Mul;
use std::ops::Sub;

use rayon::prelude::*;
use serde::Deserialize;
use serde::Serialize;

#[derive(
    Debug,
    Clone,
    Copy,
    Default,
    Serialize,
    Deserialize,
)]
/// A 2D vector.
pub struct Vector2D {
    /// The x component of the vector.
    pub x: f64,
    /// The y component of the vector.
    pub y: f64,
}

impl Vector2D {
    /// Creates a new 2D vector.
    pub fn new(
        x: f64,
        y: f64,
    ) -> Self {

        Self {
            x,
            y,
        }
    }

    pub(crate) fn norm_sq(
        &self
    ) -> f64 {

        self.x * self.x
            + self.y * self.y
    }
}

impl Add for Vector2D {
    type Output = Self;

    fn add(
        self,
        rhs: Self,
    ) -> Self {

        Self {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
        }
    }
}

impl Sub for Vector2D {
    type Output = Self;

    fn sub(
        self,
        rhs: Self,
    ) -> Self {

        Self {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
        }
    }
}

impl Mul<f64> for Vector2D {
    type Output = Self;

    fn mul(
        self,
        rhs: f64,
    ) -> Self {

        Self {
            x: self.x * rhs,
            y: self.y * rhs,
        }
    }
}

impl Div<f64> for Vector2D {
    type Output = Self;

    fn div(
        self,
        rhs: f64,
    ) -> Self {

        Self {
            x: self.x / rhs,
            y: self.y / rhs,
        }
    }
}

/// cbindgen:ignore
#[derive(
    Debug, Clone, Serialize, Deserialize,
)]
/// A particle in the SPH simulation.
pub struct Particle {
    /// The position of the particle.
    pub pos: Vector2D,
    /// The velocity of the particle.
    pub vel: Vector2D,
    /// The force acting on the particle.
    pub force: Vector2D,
    /// The density of the particle.
    pub density: f64,
    /// The pressure of the particle.
    pub pressure: f64,
    /// The mass of the particle.
    pub mass: f64,
}

#[derive(
    Debug, Clone, Serialize, Deserialize,
)]
/// The Poly6 kernel for SPH simulations.
pub struct Poly6Kernel {
    /// The squared smoothing radius.
    pub h_sq: f64,
    /// The normalization factor.
    pub factor: f64,
}

impl Poly6Kernel {
    /// Creates a new Poly6 kernel.
    pub fn new(h: f64) -> Self {

        Self {
            h_sq : h * h,
            factor : 315.0 / (64.0 * std::f64::consts::PI * h.powi(9)),
        }
    }

    pub(crate) fn value(
        &self,
        r_sq: f64,
    ) -> f64 {

        if r_sq >= self.h_sq {

            return 0.0;
        }

        let diff = self.h_sq - r_sq;

        self.factor * diff * diff * diff
    }
}

#[derive(
    Debug, Clone, Serialize, Deserialize,
)]
/// The spiky kernel for SPH simulations.
pub struct SpikyKernel {
    /// The smoothing radius.
    pub h: f64,
    /// The normalization factor.
    pub factor: f64,
}

impl SpikyKernel {
    /// Creates a new spiky kernel.
    pub fn new(h: f64) -> Self {

        Self {
            h,
            factor: -45.0
                / (std::f64::consts::PI
                    * h.powi(6)),
        }
    }

    pub(crate) fn gradient(
        &self,
        r_vec: Vector2D,
        r_norm: f64,
    ) -> Vector2D {

        if r_norm >= self.h
            || r_norm == 0.0
        {

            return Vector2D::default();
        }

        let diff = self.h - r_norm;

        r_vec
            * (self.factor
                * diff
                * diff
                / r_norm)
    }
}

#[derive(
    Debug, Clone, Serialize, Deserialize,
)]
/// The SPH system.
pub struct SPHSystem {
    /// The particles in the system.
    pub particles: Vec<Particle>,
    /// The Poly6 kernel.
    pub poly6: Poly6Kernel,
    /// The spiky kernel.
    pub spiky: SpikyKernel,
    /// The gravity vector.
    pub gravity: Vector2D,
    /// The viscosity of the fluid.
    pub viscosity: f64,
    /// The gas constant.
    pub gas_const: f64,
    /// The rest density of the fluid.
    pub rest_density: f64,
    /// The bounds of the simulation.
    pub bounds: Vector2D,
}

impl SPHSystem {
    /// Computes the density and pressure of each particle.
    pub fn compute_density_pressure(
        &mut self
    ) {

        let poly6 = &self.poly6;

        let gas_const = self.gas_const;

        let rest_density =
            self.rest_density;

        // Use a raw pointer to provide an immutable view of all particles to each thread
        // to avoid conflicts with the mutable iterator.
        let particles_ptr = self
            .particles
            .as_ptr()
            as usize;

        let n = self.particles.len();

        self.particles
            .par_iter_mut()
            .for_each(|p_i| {

                let particles_ref = unsafe {

                    std::slice::from_raw_parts(
                        particles_ptr as *const Particle,
                        n,
                    )
                };

                let mut density = 0.0;

                for p_j in particles_ref {

                    let r_vec = p_i.pos - p_j.pos;

                    density += p_j.mass * poly6.value(r_vec.norm_sq());
                }

                p_i.density = density;

                p_i.pressure = gas_const * (density - rest_density).max(0.0);
            });
    }

    /// Computes the forces acting on each particle.
    pub fn compute_forces(&mut self) {

        let spiky = &self.spiky;

        let poly6 = &self.poly6;

        let viscosity = self.viscosity;

        let gravity = self.gravity;

        let particles_ptr = self
            .particles
            .as_ptr()
            as usize;

        let n = self.particles.len();

        self.particles
            .par_iter_mut()
            .enumerate()
            .for_each(|(i, p_i)| {

                let particles_ref = unsafe {

                    std::slice::from_raw_parts(
                        particles_ptr as *const Particle,
                        n,
                    )
                };

                let mut force = Vector2D::default();

                for (j, p_j) in particles_ref
                    .iter()
                    .enumerate()
                {

                    if i == j {

                        continue;
                    }

                    let r_vec = p_i.pos - p_j.pos;

                    let r_norm = (r_vec.norm_sq()).sqrt();

                    if r_norm < spiky.h {

                        let avg_pressure = (p_i.pressure + p_j.pressure) / 2.0;

                        force =
                            force - spiky.gradient(r_vec, r_norm) * (avg_pressure / p_j.density);

                        let vel_diff = p_j.vel - p_i.vel;

                        force = force + vel_diff * (viscosity * poly6.value(r_vec.norm_sq()));
                    }
                }

                p_i.force = force + gravity * p_i.density;
            });
    }

    /// Integrates the equations of motion for each particle.
    pub fn integrate(
        &mut self,
        dt: f64,
    ) {

        let bounds = self.bounds;

        self.particles
            .par_iter_mut()
            .for_each(|p| {

                p.vel = p.vel
                    + p.force
                        * (dt / p
                            .density);

                p.pos =
                    p.pos + p.vel * dt;

                if p.pos.x < 0.0 {

                    p.vel.x *= -0.5;

                    p.pos.x = 0.0;
                }

                if p.pos.x > bounds.x {

                    p.vel.x *= -0.5;

                    p.pos.x = bounds.x;
                }

                if p.pos.y < 0.0 {

                    p.vel.y *= -0.5;

                    p.pos.y = 0.0;
                }

                if p.pos.y > bounds.y {

                    p.vel.y *= -0.5;

                    p.pos.y = bounds.y;
                }
            });
    }

    /// Updates the SPH system by one time step.
    pub fn update(
        &mut self,
        dt: f64,
    ) {

        self.compute_density_pressure();

        self.compute_forces();

        self.integrate(dt);
    }
}

/// Simulates a 2D dam break scenario using Smoothed-Particle Hydrodynamics (SPH).
///
/// This function initializes a block of fluid particles (the "dam") and simulates
/// its collapse and flow under gravity. It demonstrates the SPH method's ability
/// to model fluid dynamics without a fixed mesh.
///
/// # Returns
/// A `Vec` of tuples `(x, y)` representing the final positions of the particles.

pub fn simulate_dam_break_2d_scenario(
) -> Vec<(f64, f64)> {

    let h = 0.1;

    let mut system = SPHSystem {
        particles: Vec::new(),
        poly6: Poly6Kernel::new(h),
        spiky: SpikyKernel::new(h),
        gravity: Vector2D::new(
            0.0, -9.8,
        ),
        viscosity: 0.01,
        gas_const: 2000.0,
        rest_density: 1000.0,
        bounds: Vector2D::new(4.0, 4.0),
    };

    let particle_mass = 1.0;

    for y in (0 .. 20)
        .map(|v| f64::from(v) * h * 0.8)
    {

        for x in (0 .. 10).map(|v| {

            f64::from(v) * h * 0.8
        }) {

            system
                .particles
                .push(Particle {
                pos: Vector2D::new(
                    x,
                    y + 0.1,
                ),
                vel: Vector2D::default(
                ),
                force:
                    Vector2D::default(),
                density: 0.0,
                pressure: 0.0,
                mass: particle_mass,
            });
        }
    }

    let dt = 0.005;

    for _ in 0 .. 200 {

        system.update(dt);
    }

    system
        .particles
        .iter()
        .map(|p| (p.pos.x, p.pos.y))
        .collect()
}
