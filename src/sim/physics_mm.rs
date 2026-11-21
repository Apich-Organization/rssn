//! # Meshless methods (smoothed-particle hydrodynamics)
//!
//! | Method | Dimension | Equations | Boundaries |
//! |---|---|---|---|
//! | SPH (poly6 / spiky kernels) | 2D | weakly compressible fluid, dam break | reflecting box |
//!
//! No 1D/3D SPH, no boundary particles, no viscosity models beyond the
//! artificial term.
//!
//! Large, sparse particle sets use compact neighbor rows in particle-index
//! order. The combined density and force update builds those rows once.
//! Small or dense sets and coordinates outside the safe cell range use pair sums.

use std::ops::Add;
use std::ops::Div;
use std::ops::Mul;
use std::ops::Sub;

use rayon::prelude::*;
use serde::Deserialize;
use serde::Serialize;

use crate::kernels::par;
use crate::kernels::reuse::Scratch;
use crate::kernels::reuse::ScratchPool;
use crate::sim::sph::Grid;

/// A 2D vector.
#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize)]
pub struct Vector2D {
    /// The x component of the vector.
    pub x: f64,
    /// The y component of the vector.
    pub y: f64,
}

impl Vector2D {
    /// Creates a new 2D vector.
    #[must_use]
    pub const fn new(
        x: f64,
        y: f64,
    ) -> Self {
        Self { x, y }
    }

    pub(crate) fn norm_sq(&self) -> f64 {
        self.x.mul_add(self.x, self.y * self.y)
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
/// A particle in the SPH simulation.
#[derive(Debug, Clone, Serialize, Deserialize)]
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

/// The Poly6 kernel for SPH simulations.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Poly6Kernel {
    /// The squared smoothing radius.
    pub h_sq: f64,
    /// The normalization factor.
    pub factor: f64,
}

impl Poly6Kernel {
    /// Creates a new Poly6 kernel.
    #[must_use]
    pub fn new(h: f64) -> Self {
        Self {
            h_sq: h * h,
            factor: 315.0 / (64.0 * std::f64::consts::PI * h.powi(9)),
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

/// The spiky kernel for SPH simulations.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpikyKernel {
    /// The smoothing radius.
    pub h: f64,
    /// The normalization factor.
    pub factor: f64,
}

impl SpikyKernel {
    /// Creates a new spiky kernel.
    #[must_use]
    pub fn new(h: f64) -> Self {
        Self {
            h,
            factor: -45.0 / (std::f64::consts::PI * h.powi(6)),
        }
    }

    pub(crate) fn gradient(
        &self,
        r_vec: Vector2D,
        r_norm: f64,
    ) -> Vector2D {
        if r_norm >= self.h || r_norm == 0.0 {
            return Vector2D::default();
        }

        let diff = self.h - r_norm;

        r_vec * (self.factor * diff * diff / r_norm)
    }
}

/// The SPH system.
#[derive(Debug, Clone, Serialize)]
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

impl<'de> Deserialize<'de> for SPHSystem {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        #[derive(Deserialize)]
        struct SPHSystemData {
            particles: Vec<Particle>,
            poly6: Poly6Kernel,
            spiky: SpikyKernel,
            gravity: Vector2D,
            viscosity: f64,
            gas_const: f64,
            rest_density: f64,
            bounds: Vector2D,
        }

        let data = SPHSystemData::deserialize(deserializer)?;

        Ok(Self {
            particles: data.particles,
            poly6: data.poly6,
            spiky: data.spiky,
            gravity: data.gravity,
            viscosity: data.viscosity,
            gas_const: data.gas_const,
            rest_density: data.rest_density,
            bounds: data.bounds,
        })
    }
}

struct Neighbors {
    offsets: Scratch<'static, usize>,
    indices: Scratch<'static, usize>,
}

impl Neighbors {
    fn new(
        particles: &[Particle],
        density_radius_sq: Option<f64>,
        force_radius: Option<f64>,
    ) -> Option<Self> {
        const MAX_RECORDS: usize = 512 * 1024;
        static WORK: ScratchPool<usize> = ScratchPool::new(4, 4 * 1024 * 1024);

        let n = particles.len();
        let min_particles = if force_radius.is_some() {
            1024
        } else {
            1280
        };
        // An underflowed squared distance can fall inside the force cutoff
        // even when its coordinates are too far apart to share a cell.
        if n < min_particles
            || (density_radius_sq.is_some() && particles.iter().any(|p| !p.mass.is_finite()))
            || density_radius_sq.is_some_and(|r| !r.is_normal() || r <= 0.0)
            || force_radius.is_some_and(|r| !(r * r).is_normal() || r <= 0.0)
        {
            return None;
        }
        let density_radius_sq = density_radius_sq.unwrap_or(0.0);
        let force_radius = force_radius.unwrap_or(0.0);
        let reach = density_radius_sq.sqrt().max(force_radius);
        let points: Vec<[f64; 2]> = particles.iter().map(|p| [p.pos.x, p.pos.y]).collect();
        let grid = Grid::for_reach(&points, reach)?;
        let mut offsets = WORK.checkout_for(n + 1);
        let mut indices = WORK.checkout_for(n.saturating_mul(8).min(MAX_RECORDS));
        offsets.push(0);
        for (i, point) in points.iter().enumerate() {
            let start = indices.len();
            let mut visited = 0;
            grid.for_near(point, |j| {
                visited += 1;
                if indices.len() < MAX_RECORDS {
                    let square = (particles[i].pos - particles[j].pos).norm_sq();
                    if square < density_radius_sq || square.sqrt() < force_radius {
                        indices.push(j);
                    }
                }
            });
            if visited > n / 2 || indices.len() == MAX_RECORDS {
                return None;
            }
            indices[start..].sort_unstable();
            offsets.push(indices.len());
        }
        Some(Self { offsets, indices })
    }

    fn row(
        &self,
        i: usize,
    ) -> &[usize] {
        &self.indices[self.offsets[i]..self.offsets[i + 1]]
    }
}

impl SPHSystem {
    /// Computes the density and pressure of each particle.
    pub fn compute_density_pressure(&mut self) {
        let neighbors = Neighbors::new(&self.particles, Some(self.poly6.h_sq), None);
        self.density_pressure_with(neighbors.as_ref());
    }

    fn density_pressure_with(
        &mut self,
        neighbors: Option<&Neighbors>,
    ) {
        let poly6 = &self.poly6;
        let gas_const = self.gas_const;
        let rest_density = self.rest_density;
        let results: Vec<_> = neighbors.map_or_else(
            || self.direct_density_pressure(),
            |rows| {
                par::map_indexed(self.particles.len(), 128, |i| {
                    let points = rows.row(i).iter().map(|&j| &self.particles[j]);
                    density_at(&self.particles[i], points, poly6, gas_const, rest_density)
                })
            },
        );
        for (p, (density, pressure)) in self.particles.iter_mut().zip(results) {
            p.density = density;

            p.pressure = pressure;
        }
    }

    // Keep the contiguous all-pairs loop separate from indirect neighbor lookups.
    #[inline(never)]
    fn direct_density_pressure(&self) -> Vec<(f64, f64)> {
        let poly6 = &self.poly6;
        let gas_const = self.gas_const;
        let rest_density = self.rest_density;
        self.particles
            .par_iter()
            .map(|p_i| {
                let mut density = 0.0;
                for p_j in &self.particles {
                    let r_vec = p_i.pos - p_j.pos;
                    density += p_j.mass * poly6.value(r_vec.norm_sq());
                }
                (density, gas_const * (density - rest_density).max(0.0))
            })
            .collect()
    }

    /// Computes the forces acting on each particle.
    pub fn compute_forces(&mut self) {
        let neighbors = Neighbors::new(&self.particles, None, Some(self.spiky.h));
        self.forces_with(neighbors.as_ref());
    }

    fn forces_with(
        &mut self,
        neighbors: Option<&Neighbors>,
    ) {
        let n = self.particles.len();
        let forces = neighbors.map_or_else(
            || {
                par::map_indexed(n, 48, |i| {
                    self.force_at(i, self.particles.iter().enumerate())
                })
            },
            |rows| {
                par::map_indexed(n, 128, |i| {
                    let particles = rows.row(i).iter().map(|&j| (j, &self.particles[j]));
                    self.force_at(i, particles)
                })
            },
        );

        for (p, force) in self.particles.iter_mut().zip(forces) {
            p.force = force;
        }
    }

    fn force_at<'a>(
        &self,
        i: usize,
        particles: impl Iterator<Item = (usize, &'a Particle)>,
    ) -> Vector2D {
        let p_i = &self.particles[i];
        let mut force = Vector2D::default();
        for (j, p_j) in particles {
            if i == j {
                continue;
            }
            let r_vec = p_i.pos - p_j.pos;
            let square = r_vec.norm_sq();
            let r_norm = square.sqrt();
            if r_norm < self.spiky.h {
                let avg_pressure = f64::midpoint(p_i.pressure, p_j.pressure);
                force = force - self.spiky.gradient(r_vec, r_norm) * (avg_pressure / p_j.density);
                let vel_diff = p_j.vel - p_i.vel;
                force = force + vel_diff * (self.viscosity * self.poly6.value(square));
            }
        }
        force + self.gravity * p_i.density
    }

    /// Updates density, pressure and force for the current positions and velocities.
    /// Reuses one neighborhood search for both interaction sweeps.
    pub fn compute_density_and_forces(&mut self) {
        let neighbors = Neighbors::new(&self.particles, Some(self.poly6.h_sq), Some(self.spiky.h));
        self.density_pressure_with(neighbors.as_ref());
        self.forces_with(neighbors.as_ref());
    }

    /// Integrates the equations of motion for each particle.
    pub fn integrate(
        &mut self,
        dt: f64,
    ) {
        let bounds = self.bounds;

        par::for_each_mut(&mut self.particles, 128, |_, p| {
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

    /// Updates the SPH system by one time step.
    pub fn update(
        &mut self,
        dt: f64,
    ) {
        self.compute_density_and_forces();
        self.integrate(dt);
    }
}

#[inline]
fn density_at<'a>(
    p_i: &Particle,
    particles: impl Iterator<Item = &'a Particle>,
    poly6: &Poly6Kernel,
    gas_const: f64,
    rest_density: f64,
) -> (f64, f64) {
    let mut density = 0.0;
    for p_j in particles {
        let r_vec = p_i.pos - p_j.pos;
        density += p_j.mass * poly6.value(r_vec.norm_sq());
    }
    (density, gas_const * (density - rest_density).max(0.0))
}

/// Simulates a 2D dam break scenario using Smoothed-Particle Hydrodynamics (SPH).
///
/// This function initializes a block of fluid particles (the "dam") and simulates
/// its collapse and flow under gravity. It demonstrates the SPH method's ability
/// to model fluid dynamics without a fixed mesh.
///
/// # Returns
/// A `Vec` of tuples `(x, y)` representing the final positions of the particles.
#[must_use]
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
