use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::ops::{Add, Div, Mul, Sub};

#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize)]
pub struct Vector2D {
    pub x: f64,
    pub y: f64,
}
impl Vector2D {
    pub fn new(x: f64, y: f64) -> Self {
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
impl Div<f64> for Vector2D {
    type Output = Self;
    fn div(self, rhs: f64) -> Self {
        Self {
            x: self.x / rhs,
            y: self.y / rhs,
        }
    }
}
/// cbindgen:ignore
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Particle {
    pub pos: Vector2D,
    pub vel: Vector2D,
    pub force: Vector2D,
    pub density: f64,
    pub pressure: f64,
    pub mass: f64,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Poly6Kernel {
    pub h_sq: f64,
    pub factor: f64,
}
impl Poly6Kernel {
    pub fn new(h: f64) -> Self {
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
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpikyKernel {
    pub h: f64,
    pub factor: f64,
}
impl SpikyKernel {
    pub fn new(h: f64) -> Self {
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
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SPHSystem {
    pub particles: Vec<Particle>,
    pub poly6: Poly6Kernel,
    pub spiky: SpikyKernel,
    pub gravity: Vector2D,
    pub viscosity: f64,
    pub gas_const: f64,
    pub rest_density: f64,
    pub bounds: Vector2D,
}
impl SPHSystem {
    pub fn compute_density_pressure(&mut self) {
        let poly6 = &self.poly6;
        let gas_const = self.gas_const;
        let rest_density = self.rest_density;
        
        // Use a raw pointer to provide an immutable view of all particles to each thread
        // to avoid conflicts with the mutable iterator.
        let particles_ptr = self.particles.as_ptr() as usize;
        let n = self.particles.len();

        self.particles.par_iter_mut().for_each(|p_i| {
            let particles_ref = unsafe { std::slice::from_raw_parts(particles_ptr as *const Particle, n) };
            let mut density = 0.0;
            for p_j in particles_ref {
                let r_vec = p_i.pos - p_j.pos;
                density += p_j.mass * poly6.value(r_vec.norm_sq());
            }
            p_i.density = density;
            p_i.pressure = gas_const * (density - rest_density).max(0.0);
        });
    }

    pub fn compute_forces(&mut self) {
        let spiky = &self.spiky;
        let poly6 = &self.poly6;
        let viscosity = self.viscosity;
        let gravity = self.gravity;
        
        let particles_ptr = self.particles.as_ptr() as usize;
        let n = self.particles.len();

        self.particles.par_iter_mut().enumerate().for_each(|(i, p_i)| {
            let particles_ref = unsafe { std::slice::from_raw_parts(particles_ptr as *const Particle, n) };
            let mut force = Vector2D::default();
            for (j, p_j) in particles_ref.iter().enumerate() {
                if i == j {
                    continue;
                }
                let r_vec = p_i.pos - p_j.pos;
                let r_norm = (r_vec.norm_sq()).sqrt();
                if r_norm < spiky.h {
                    let avg_pressure = (p_i.pressure + p_j.pressure) / 2.0;
                    force = force
                        - spiky.gradient(r_vec, r_norm)
                            * (avg_pressure / p_j.density);
                    let vel_diff = p_j.vel - p_i.vel;
                    force = force + vel_diff * (viscosity * poly6.value(r_vec.norm_sq()));
                }
            }
            p_i.force = force + gravity * p_i.density;
        });
    }
    pub fn integrate(&mut self, dt: f64) {
        let bounds = self.bounds;
        self.particles.par_iter_mut().for_each(|p| {
            p.vel = p.vel + p.force * (dt / p.density);
            p.pos = p.pos + p.vel * dt;
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
    pub fn update(&mut self, dt: f64) {
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
pub fn simulate_dam_break_2d_scenario() -> Vec<(f64, f64)> {
    let h = 0.1;
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
    let particle_mass = 1.0;
    for y in (0..20).map(|v| f64::from(v) * h * 0.8) {
        for x in (0..10).map(|v| f64::from(v) * h * 0.8) {
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
    let dt = 0.005;
    for _ in 0..200 {
        system.update(dt);
    }
    system
        .particles
        .iter()
        .map(|p| (p.pos.x, p.pos.y))
        .collect()
}
