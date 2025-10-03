// src/physics/physics_mm.rs
// A module for Meshfree Methods, implementing a Smoothed-Particle Hydrodynamics (SPH) solver.

use std::ops::{Add, Mul, Sub};

// --- Core Structs ---

#[derive(Debug, Clone, Copy, Default)]
pub struct Vector2D {
    pub x: f64,
    pub y: f64,
}

impl Vector2D {
    pub(crate) fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }
    pub(crate) fn norm_sq(&self) -> f64 {
        self.x * self.x + self.y * self.y
    }
}

impl Add for Vector2D {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
        }
    }
}
impl Sub for Vector2D {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
        }
    }
}
impl Mul<f64> for Vector2D {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Self {
            x: self.x * rhs,
            y: self.y * rhs,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Particle {
    pub pos: Vector2D,   // Position
    pub vel: Vector2D,   // Velocity
    pub force: Vector2D, // Force accumulator
    pub density: f64,
    pub pressure: f64,
    pub mass: f64,
}

// --- SPH Kernels ---

pub struct Poly6Kernel {
    h_sq: f64,
    factor: f64,
}

impl Poly6Kernel {
    pub(crate) fn new(h: f64) -> Self {
        Self {
            h_sq: h * h,
            factor: 315.0 / (64.0 * std::f64::consts::PI * h.powi(9)),
        }
    }

    pub(crate) fn value(&self, r_sq: f64) -> f64 {
        if r_sq >= self.h_sq {
            return 0.0;
        }
        let diff = self.h_sq - r_sq;
        self.factor * diff * diff * diff
    }
}

pub struct SpikyKernel {
    h: f64,
    factor: f64,
}

impl SpikyKernel {
    pub(crate) fn new(h: f64) -> Self {
        Self {
            h,
            factor: -45.0 / (std::f64::consts::PI * h.powi(6)),
        }
    }

    pub(crate) fn gradient(&self, r_vec: Vector2D, r_norm: f64) -> Vector2D {
        if r_norm >= self.h || r_norm == 0.0 {
            return Vector2D::default();
        }
        let diff = self.h - r_norm;
        r_vec * (self.factor * diff * diff / r_norm)
    }
}

// --- SPH System ---

pub struct SPHSystem {
    particles: Vec<Particle>,
    poly6: Poly6Kernel,
    spiky: SpikyKernel,
    gravity: Vector2D,
    viscosity: f64,
    gas_const: f64,
    rest_density: f64,
    bounds: Vector2D,
}

impl SPHSystem {
    pub(crate) fn compute_density_pressure(&mut self) {
        for i in 0..self.particles.len() {
            let mut density = 0.0;
            for j in 0..self.particles.len() {
                let r_vec = self.particles[i].pos - self.particles[j].pos;
                density += self.particles[j].mass * self.poly6.value(r_vec.norm_sq());
            }
            self.particles[i].density = density;
            self.particles[i].pressure = self.gas_const * (density - self.rest_density).max(0.0);
        }
    }

    pub(crate) fn compute_forces(&mut self) {
        for i in 0..self.particles.len() {
            let mut force = Vector2D::default();
            for j in 0..self.particles.len() {
                if i == j {
                    continue;
                }
                let r_vec = self.particles[i].pos - self.particles[j].pos;
                let r_norm = (r_vec.norm_sq()).sqrt();

                // Pressure force
                let avg_pressure = (self.particles[i].pressure + self.particles[j].pressure) / 2.0;
                force = force
                    - self.spiky.gradient(r_vec, r_norm)
                        * (avg_pressure / self.particles[j].density);

                // Viscosity force
                let vel_diff = self.particles[j].vel - self.particles[i].vel;
                force = force + vel_diff * (self.viscosity * self.poly6.value(r_vec.norm_sq()));
            }
            self.particles[i].force = force + self.gravity * self.particles[i].density;
        }
    }

    pub(crate) fn integrate(&mut self, dt: f64) {
        for p in &mut self.particles {
            // Leapfrog integration (or Verlet-style)
            p.vel = p.vel + p.force * (dt / p.density);
            p.pos = p.pos + p.vel * dt;

            // Boundary conditions (collision with walls)
            if p.pos.x < 0.0 {
                p.vel.x *= -0.5;
                p.pos.x = 0.0;
            }
            if p.pos.x > self.bounds.x {
                p.vel.x *= -0.5;
                p.pos.x = self.bounds.x;
            }
            if p.pos.y < 0.0 {
                p.vel.y *= -0.5;
                p.pos.y = 0.0;
            }
            if p.pos.y > self.bounds.y {
                p.vel.y *= -0.5;
                p.pos.y = self.bounds.y;
            }
        }
    }

    pub fn update(&mut self, dt: f64) {
        self.compute_density_pressure();
        self.compute_forces();
        self.integrate(dt);
    }
}

// --- Example Scenario ---

/// Simulates a 2D dam break scenario using Smoothed-Particle Hydrodynamics (SPH).
///
/// This function initializes a block of fluid particles (the "dam") and simulates
/// its collapse and flow under gravity. It demonstrates the SPH method's ability
/// to model fluid dynamics without a fixed mesh.
///
/// # Returns
/// A `Vec` of tuples `(x, y)` representing the final positions of the particles.
pub fn simulate_dam_break_2d_scenario() -> Vec<(f64, f64)> {
    let h = 0.1; // Smoothing radius
    let mut system = SPHSystem {
        particles: Vec::new(),
        poly6: Poly6Kernel::new(h),
        spiky: SpikyKernel::new(h),
        gravity: Vector2D::new(0.0, -9.8),
        viscosity: 0.01,
        gas_const: 2000.0,
        rest_density: 1000.0,
        bounds: Vector2D::new(4.0, 4.0),
    };

    // Initialize particles in a block (the "dam")
    let particle_mass = 1.0;
    for y in (0..20).map(|v| v as f64 * h * 0.8) {
        for x in (0..10).map(|v| v as f64 * h * 0.8) {
            system.particles.push(Particle {
                pos: Vector2D::new(x, y + 0.1),
                vel: Vector2D::default(),
                force: Vector2D::default(),
                density: 0.0,
                pressure: 0.0,
                mass: particle_mass,
            });
        }
    }

    // Simulation loop
    let dt = 0.005;
    for _ in 0..200 {
        // Run for 200 steps
        system.update(dt);
    }

    // Return final positions
    system
        .particles
        .iter()
        .map(|p| (p.pos.x, p.pos.y))
        .collect()
}
