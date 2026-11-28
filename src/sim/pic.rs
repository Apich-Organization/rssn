//! Electrostatic particle-in-cell with `CIC` and four-point cubic assignment.
//!
//! # Electrostatic particle-in-cell on a periodic one-dimensional mesh
//!
//! Particles carry a charge and a mass; each step deposits the charge on a
//! uniform periodic mesh, solves the Poisson equation in Fourier space and
//! gathers the mesh field back to the particles (a leapfrog `PIC` cycle).
//!
//! The charge-assignment and field-gathering shape functions are selectable at
//! configuration time:
//!
//! | Shape | Alias | Order | Cells per particle |
//! |---|---|---|---|
//! | [`Shape::Linear`] | `CIC` (cloud-in-cell) | 2 | 2 |
//! | [`Shape::Cubic`] | four-point cubic (Lagrange) | 4 | 4 |
//!
//! Deposition and gathering use the *same* shape, so the total force on a
//! periodic neutral system vanishes to roundoff and the single-particle
//! self-force is zero. The linear mesh field is the spectral Poisson solution
//! with the `CIC` structure factor `sinc^2(k dx / 2)` divided out. The cubic
//! shape is the four-point cubic interpolation kernel, which reproduces
//! cubics, so it needs no deconvolution and its smooth-field error is fourth
//! order. The positive B-spline clouds (`TSC`, `PQS`) are deliberately not
//! used: their second moment is non-zero, so they stay second order for a
//! smooth field and only lower the aliasing constant.
//!
//! The linear path is the reference `CIC` implementation and is dispatched
//! separately from the cubic stencil, so selecting [`Shape::Linear`] reproduces
//! the old code bit for bit.
//!
//! Limitations: one spatial dimension, a uniform periodic mesh and
//! electrostatics only (no magnetic field, no relativistic pusher, no
//! two- or three-dimensional tensor-product shapes). The cubic kernel has
//! negative weights, so the deposited density can be locally negative even for
//! positive charges (it conserves the total charge and the first moment but is
//! not a positive cloud). `rho` and `J` are the assignment-smoothed fields, so
//! the continuity residual is limited by the time difference and by the
//! assignment aliasing rather than being enforced to roundoff, and there is no
//! charge-conserving (Esirkepov) current deposition. There are no collisions,
//! no external fields and no boundary interaction beyond the periodic wrap.
//! The cubic shape touches 4 cells per particle instead of 2, so deposition and
//! gathering cost about twice the linear path.

use std::f64::consts::PI;

use crate::kernels::fftconv::Complex;
use crate::kernels::fftconv::fft;
use crate::kernels::numeric_common::scale_product;

const PARTICLE_GRAIN: usize = 16_384;

/// Charge-assignment and field-gathering shape function.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Shape {
    /// Cloud-in-cell (`CIC`): the linear shape, second order, two cells.
    Linear,
    /// Four-point cubic (Lagrange) shape, fourth order, four cells.
    Cubic,
}

impl Shape {
    /// Number of mesh cells a particle contributes to (the stencil width).
    #[must_use]
    pub const fn cells(self) -> usize {
        match self {
            | Self::Linear => 2,
            | Self::Cubic => 4,
        }
    }

    /// Order of the shape's interpolation error (the field error at particles).
    #[must_use]
    pub const fn order(self) -> usize {
        match self {
            | Self::Linear => 2,
            | Self::Cubic => 4,
        }
    }
}

/// A periodic uniform mesh.
#[derive(Clone, Copy, Debug)]
pub struct Grid {
    cells: usize,
    length: f64,
    spacing: f64,
}

impl Grid {
    /// Creates a periodic mesh with `cells` cells on a box of length `length`.
    ///
    /// # Errors
    /// If `cells` is not a power of two of at least 8, `length` is not
    /// finite and positive, or the spacing rounds to zero.
    pub fn new(
        cells: usize,
        length: f64,
    ) -> Result<Self, String> {
        if cells < 8 || !cells.is_power_of_two() {
            return Err("PIC needs a power-of-two cell count of at least 8".into());
        }
        if !length.is_finite() || length <= 0.0 {
            return Err("PIC needs a finite positive box length".into());
        }
        let spacing = length / cells as f64;
        if spacing == 0.0 {
            return Err("PIC grid spacing rounds to zero".into());
        }
        Ok(Self {
            cells,
            length,
            spacing,
        })
    }

    /// Number of cells.
    #[must_use]
    pub const fn cells(&self) -> usize {
        self.cells
    }

    /// Box length.
    #[must_use]
    pub const fn length(&self) -> f64 {
        self.length
    }

    /// Cell size `length / cells`.
    #[must_use]
    pub const fn spacing(&self) -> f64 {
        self.spacing
    }

    fn phase(
        &self,
        x: f64,
    ) -> (usize, f64) {
        let wrapped = if (0.0..self.length).contains(&x) {
            x
        } else {
            x.rem_euclid(self.length)
        };
        let u = if self.spacing.is_normal() {
            wrapped / self.spacing
        } else {
            // A rounded subnormal spacing can move a point to a different cell.
            (wrapped / self.length) * self.cells as f64
        };
        let floor = u.floor();
        (floor as usize & (self.cells - 1), u - floor)
    }

    fn linear(
        &self,
        x: f64,
    ) -> Stencil<2> {
        let (start, f) = self.phase(x);
        Stencil {
            start,
            weights: [1.0 - f, f],
        }
    }

    fn cubic(
        &self,
        x: f64,
    ) -> Stencil<4> {
        let (cell, f) = self.phase(x);
        let a = f + 1.0;
        let c = f - 1.0;
        let d = f - 2.0;
        Stencil {
            start: cell.wrapping_sub(1) & (self.cells - 1),
            weights: [
                -f * c * d / 6.0,
                a * c * d / 2.0,
                -a * f * d / 2.0,
                a * f * c / 6.0,
            ],
        }
    }

    fn charge_density<const N: usize>(
        &self,
        stencils: impl Iterator<Item = Stencil<N>>,
        q: &[f64],
    ) -> Vec<f64> {
        let mut rho = vec![0.0; self.cells];
        let mask = self.cells - 1;
        if self.spacing.is_normal() {
            for (stencil, &charge) in stencils.zip(q) {
                stencil.scatter(&mut rho, mask, |w| charge * w);
            }
            for density in &mut rho {
                *density /= self.spacing;
            }
        } else {
            // On a subnormal mesh, forming q*w can erase a finite density.
            for (stencil, &charge) in stencils.zip(q) {
                stencil.scatter(&mut rho, mask, |weight| {
                    let direct = charge * weight / self.spacing;
                    if charge == 0.0 || weight == 0.0 {
                        direct
                    } else {
                        scale_product(
                            charge.signum() * weight.signum(),
                            &[charge.abs(), weight.abs(), self.cells as f64],
                            &[self.length],
                        )
                        .unwrap_or(direct)
                    }
                });
            }
        }
        rho
    }

    /// Deposits `q` at `pos` and returns the charge density per unit length.
    ///
    /// Positions are wrapped periodically. The two slices are zipped, so a
    /// trailing element of the longer slice is ignored.
    #[must_use]
    pub fn deposit(
        &self,
        pos: &[f64],
        q: &[f64],
        shape: Shape,
    ) -> Vec<f64> {
        match shape {
            | Shape::Linear => self.charge_density(pos.iter().map(|&x| self.linear(x)), q),
            | Shape::Cubic => self.charge_density(pos.iter().map(|&x| self.cubic(x)), q),
        }
    }

    /// Gathers a mesh quantity to the particles with the same shape weights.
    ///
    /// # Errors
    /// If `field` does not have one entry per cell.
    pub fn gather(
        &self,
        field: &[f64],
        pos: &[f64],
        shape: Shape,
    ) -> Result<Vec<f64>, String> {
        if field.len() != self.cells {
            return Err("PIC gather needs one field value per cell".into());
        }
        let mask = self.cells - 1;
        Ok(match shape {
            | Shape::Linear => {
                crate::kernels::par::map_slice(pos, PARTICLE_GRAIN, |&x| {
                    self.linear(x).sample(field, mask)
                })
            },
            | Shape::Cubic => {
                crate::kernels::par::map_slice(pos, PARTICLE_GRAIN, |&x| {
                    self.cubic(x).sample(field, mask)
                })
            },
        })
    }

    /// Prepares the periodic weights for repeated transfers at fixed positions.
    ///
    /// # Errors
    ///
    /// A position is non-finite.
    pub fn prepare_assignment(
        &self,
        pos: &[f64],
        shape: Shape,
    ) -> Result<AssignmentPlan, String> {
        if pos.iter().any(|x| !x.is_finite()) {
            return Err("PIC assignment positions must be finite".into());
        }
        let stencils = match shape {
            | Shape::Linear => {
                AssignmentStencils::Linear(crate::kernels::par::map_slice(
                    pos,
                    PARTICLE_GRAIN,
                    |&x| self.linear(x),
                ))
            },
            | Shape::Cubic => {
                AssignmentStencils::Cubic(crate::kernels::par::map_slice(
                    pos,
                    PARTICLE_GRAIN,
                    |&x| self.cubic(x),
                ))
            },
        };
        Ok(AssignmentPlan {
            grid: *self,
            stencils,
            particles: pos.len(),
        })
    }

    /// Solves `-d^2 phi / dx^2 = rho` and returns `phi` and `E = -d phi / dx`.
    ///
    /// The zero mode is dropped, i.e. a uniform neutralising background is
    /// assumed. The `CIC` structure factor is divided out of the linear mesh
    /// density; the cubic kernel reproduces cubics, so it needs no correction.
    ///
    /// # Errors
    /// If `rho` does not have one entry per cell.
    pub fn solve(
        &self,
        rho: &[f64],
        shape: Shape,
    ) -> Result<Field, String> {
        if rho.len() != self.cells {
            return Err("PIC Poisson solve needs one density per cell".into());
        }
        let n = self.cells;
        let samples: Vec<Complex> = rho.iter().map(|&r| (r, 0.0)).collect();
        let rk = fft(&samples, false);
        let mut phik = vec![(0.0, 0.0); n];
        let mut ek = vec![(0.0, 0.0); n];
        for (m, (p, e)) in phik.iter_mut().zip(ek.iter_mut()).enumerate() {
            if m == 0 {
                continue;
            }
            let mf = if m <= n / 2 {
                m as f64
            } else {
                m as f64 - n as f64
            };
            let k = 2.0 * PI * mf / self.length;
            let x = 0.5 * k * self.spacing;
            let sinc = if x.abs() < 1.0e-12 {
                1.0
            } else {
                x.sin() / x
            };
            let corr = match shape {
                | Shape::Linear => sinc.powi(2),
                | Shape::Cubic => 1.0,
            };
            let rkm = rk.get(m).copied().unwrap_or((0.0, 0.0));
            // The inverse FFT carries the 1/n, so `a` is n times the Fourier
            // coefficient rho_k: phi_k = rho_k / k^2 and E_k = -i rho_k / k.
            let a = (rkm.0 / corr, rkm.1 / corr);
            *p = (a.0 / (k * k), a.1 / (k * k));
            *e = (a.1 / k, -a.0 / k);
        }
        let phi = fft(&phik, true).iter().map(|c| c.0).collect();
        let e = fft(&ek, true).iter().map(|c| c.0).collect();
        Ok(Field { phi, e })
    }
}

/// Mesh potential and electric field.
#[derive(Clone, Debug)]
pub struct Field {
    /// Potential at the cell centres.
    pub phi: Vec<f64>,
    /// Electric field at the cell centres (`-d phi / dx`).
    pub e: Vec<f64>,
}

#[derive(Clone, Copy, Debug)]
struct Stencil<const N: usize> {
    start: usize,
    weights: [f64; N],
}

impl<const N: usize> Stencil<N> {
    fn scatter(
        self,
        rho: &mut [f64],
        mask: usize,
        scale: impl Fn(f64) -> f64,
    ) {
        for (k, weight) in self.weights.into_iter().enumerate() {
            rho[self.start.wrapping_add(k) & mask] += scale(weight);
        }
    }

    fn sample(
        self,
        field: &[f64],
        mask: usize,
    ) -> f64 {
        if N == 2 {
            // Keep the original CIC association, including signed zero.
            return self.weights[0] * field[self.start]
                + self.weights[1] * field[self.start.wrapping_add(1) & mask];
        }
        self.weights
            .into_iter()
            .enumerate()
            .fold(0.0, |acc, (k, w)| {
                acc + w * field[self.start.wrapping_add(k) & mask]
            })
    }
}

#[derive(Clone, Debug)]
enum AssignmentStencils {
    Linear(Vec<Stencil<2>>),
    Cubic(Vec<Stencil<4>>),
}

/// Prepared periodic particle-to-mesh and mesh-to-particle transfers.
///
/// Constructed by [`Grid::prepare_assignment`], a plan retains the stencil of
/// each particle in the original order. Deposition and gathering share the
/// same weights for every field evaluation.
#[derive(Clone, Debug)]
pub struct AssignmentPlan {
    grid: Grid,
    stencils: AssignmentStencils,
    particles: usize,
}

impl AssignmentPlan {
    pub(crate) const fn heap_bytes(&self) -> usize {
        match &self.stencils {
            | AssignmentStencils::Linear(s) => s.capacity() * size_of::<Stencil<2>>(),
            | AssignmentStencils::Cubic(s) => s.capacity() * size_of::<Stencil<4>>(),
        }
    }

    /// Deposits one charge per prepared particle as density per unit length.
    ///
    /// # Errors
    ///
    /// The charge array length differs from the prepared particle count.
    pub fn deposit(
        &self,
        q: &[f64],
    ) -> Result<Vec<f64>, String> {
        if q.len() != self.particles {
            return Err("PIC assignment needs one charge per prepared particle".into());
        }
        Ok(match &self.stencils {
            | AssignmentStencils::Linear(s) => self.grid.charge_density(s.iter().copied(), q),
            | AssignmentStencils::Cubic(s) => self.grid.charge_density(s.iter().copied(), q),
        })
    }

    /// Evaluates the prepared stencils on a mesh field.
    ///
    /// # Errors
    ///
    /// The field does not have one entry per grid cell.
    pub fn gather(
        &self,
        field: &[f64],
    ) -> Result<Vec<f64>, String> {
        if field.len() != self.grid.cells {
            return Err("PIC gather needs one field value per cell".into());
        }
        let mask = self.grid.cells - 1;
        Ok(match &self.stencils {
            | AssignmentStencils::Linear(s) => {
                crate::kernels::par::map_slice(s, PARTICLE_GRAIN, |&stencil| {
                    stencil.sample(field, mask)
                })
            },
            | AssignmentStencils::Cubic(s) => {
                crate::kernels::par::map_slice(s, PARTICLE_GRAIN, |&stencil| {
                    stencil.sample(field, mask)
                })
            },
        })
    }
}

/// Configuration of a [`Pic`] run.
#[derive(Clone, Copy, Debug)]
pub struct PicConfig {
    /// Number of mesh cells (a power of two, at least 8).
    pub cells: usize,
    /// Box length (finite and positive).
    pub length: f64,
    /// Shape used for both deposition and gathering.
    pub shape: Shape,
}

impl PicConfig {
    /// Validates a configuration.
    ///
    /// # Errors
    /// If the mesh or the box length is invalid.
    pub fn new(
        cells: usize,
        length: f64,
        shape: Shape,
    ) -> Result<Self, String> {
        Grid::new(cells, length)?;
        Ok(Self { cells, length, shape })
    }
}

/// A one-dimensional electrostatic `PIC` system with periodic boundaries.
#[derive(Clone, Debug)]
pub struct Pic {
    grid: Grid,
    shape: Shape,
    pos: Vec<f64>,
    vel: Vec<f64>,
    charge: Vec<f64>,
    mass: Vec<f64>,
    field: Vec<f64>,
}

impl Pic {
    /// Creates a system from equal-length particle arrays.
    ///
    /// # Errors
    /// If the configuration is invalid, the arrays differ in length, there are
    /// no particles, or a mass is not finite and positive.
    pub fn new(
        cfg: &PicConfig,
        pos: Vec<f64>,
        vel: Vec<f64>,
        charge: Vec<f64>,
        mass: Vec<f64>,
    ) -> Result<Self, String> {
        let grid = Grid::new(cfg.cells, cfg.length)?;
        let n = pos.len();
        if n == 0 {
            return Err("PIC needs at least one particle".into());
        }
        if vel.len() != n || charge.len() != n || mass.len() != n {
            return Err("PIC particle arrays must have equal length".into());
        }
        if mass.iter().any(|&m| !m.is_finite() || m <= 0.0) {
            return Err("PIC masses must be finite and positive".into());
        }
        let mut pos = pos;
        for x in &mut pos {
            *x = x.rem_euclid(cfg.length);
        }
        Ok(Self {
            grid,
            shape: cfg.shape,
            pos,
            vel,
            charge,
            mass,
            field: vec![0.0; cfg.cells],
        })
    }

    /// The mesh.
    #[must_use]
    pub const fn grid(&self) -> &Grid {
        &self.grid
    }

    /// The shape function.
    #[must_use]
    pub const fn shape(&self) -> Shape {
        self.shape
    }

    /// Particle positions.
    #[must_use]
    pub fn positions(&self) -> &[f64] {
        &self.pos
    }

    /// Particle velocities.
    #[must_use]
    pub fn velocities(&self) -> &[f64] {
        &self.vel
    }

    /// The electric field at the cell centres from the last half step.
    #[must_use]
    pub fn field(&self) -> &[f64] {
        &self.field
    }

    /// Total charge.
    #[must_use]
    pub fn total_charge(&self) -> f64 {
        self.charge.iter().sum()
    }

    /// Total momentum `sum m v`.
    #[must_use]
    pub fn total_momentum(&self) -> f64 {
        self.mass
            .iter()
            .zip(self.vel.iter())
            .map(|(&m, &v)| m * v)
            .sum()
    }

    /// One leapfrog step: half kick, drift, field solve, half kick.
    ///
    /// # Errors
    /// If `dt` is not finite and positive, or the mesh solve fails.
    pub fn step(
        &mut self,
        dt: f64,
    ) -> Result<(), String> {
        if !dt.is_finite() || dt <= 0.0 {
            return Err("PIC needs a finite positive time step".into());
        }
        self.accelerate(0.5 * dt)?;
        let length = self.grid.length();
        for (x, &v) in self.pos.iter_mut().zip(self.vel.iter()) {
            *x = (*x + v * dt).rem_euclid(length);
        }
        self.accelerate(0.5 * dt)
    }

    /// Half a kick with the field at the current positions.
    fn accelerate(
        &mut self,
        dt: f64,
    ) -> Result<(), String> {
        let rho = self.grid.deposit(&self.pos, &self.charge, self.shape);
        let solved = self.grid.solve(&rho, self.shape)?;
        let e = self.grid.gather(&solved.e, &self.pos, self.shape)?;
        for ((v, &ei), (&qi, &mi)) in self
            .vel
            .iter_mut()
            .zip(e.iter())
            .zip(self.charge.iter().zip(self.mass.iter()))
        {
            *v += qi / mi * ei * dt;
        }
        self.field = solved.e;
        Ok(())
    }
}
