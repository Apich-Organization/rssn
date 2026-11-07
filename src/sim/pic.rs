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
    /// If `cells` is not a power of two of at least 8, or `length` is not
    /// finite and positive.
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
        Ok(Self {
            cells,
            length,
            spacing: length / cells as f64,
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
        let mut rho = vec![0.0; self.cells];
        match shape {
            // The reference CIC path: kept separate from the generic stencils
            // so that selecting `Linear` is bit-identical to the old scheme.
            | Shape::Linear => {
                for (&x, &qi) in pos.iter().zip(q.iter()) {
                    let u = x.rem_euclid(self.length) / self.spacing;
                    let fl = u.floor();
                    let f = u - fl;
                    let i0 = (fl as i64).rem_euclid(self.cells as i64) as usize;
                    let i1 = i0.wrapping_add(1).rem_euclid(self.cells);
                    let w0 = 1.0 - f;
                    if let Some(c) = rho.get_mut(i0) {
                        *c += qi * w0;
                    }
                    if let Some(c) = rho.get_mut(i1) {
                        *c += qi * f;
                    }
                }
            },
            | Shape::Cubic => {
                for (&x, &qi) in pos.iter().zip(q.iter()) {
                    let u = x.rem_euclid(self.length) / self.spacing;
                    let (base, w, count) = stencil(u, self.cells);
                    for k in 0..count {
                        let idx = base.wrapping_add(k).rem_euclid(self.cells);
                        if let (Some(c), Some(&wk)) = (rho.get_mut(idx), w.get(k)) {
                            *c += qi * wk;
                        }
                    }
                }
            },
        }
        for c in &mut rho {
            *c /= self.spacing;
        }
        rho
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
        let mut out = Vec::with_capacity(pos.len());
        match shape {
            | Shape::Linear => {
                for &x in pos {
                    let u = x.rem_euclid(self.length) / self.spacing;
                    let fl = u.floor();
                    let f = u - fl;
                    let i0 = (fl as i64).rem_euclid(self.cells as i64) as usize;
                    let i1 = i0.wrapping_add(1).rem_euclid(self.cells);
                    let w0 = 1.0 - f;
                    let v0 = field.get(i0).copied().unwrap_or(0.0);
                    let v1 = field.get(i1).copied().unwrap_or(0.0);
                    out.push(w0 * v0 + f * v1);
                }
            },
            | Shape::Cubic => {
                for &x in pos {
                    let u = x.rem_euclid(self.length) / self.spacing;
                    let (base, w, count) = stencil(u, self.cells);
                    let mut acc = 0.0;
                    for k in 0..count {
                        let idx = base.wrapping_add(k).rem_euclid(self.cells);
                        if let (Some(&vk), Some(&wk)) = (field.get(idx), w.get(k)) {
                            acc += wk * vk;
                        }
                    }
                    out.push(acc);
                }
            },
        }
        Ok(out)
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

/// The first cell and weights of the four-point cubic stencil.
///
/// The weight at index `k` belongs to the cell `start + k` (modulo the mesh);
/// the stencil starts one cell below the floor cell. The weights are the
/// Lagrange basis of the four nodes `{-1, 0, 1, 2}` in cell units, evaluated at
/// the fractional offset `f`, so they reproduce cubics and the smooth-field
/// error is fourth order. `Linear` never reaches this helper; its dedicated
/// path lives in [`Grid::deposit`] and [`Grid::gather`].
fn stencil(
    u: f64,
    cells: usize,
) -> (usize, [f64; 4], usize) {
    let fl = u.floor();
    let f = u - fl;
    let start = (fl as i64)
        .rem_euclid(cells as i64)
        .wrapping_sub(1)
        .rem_euclid(cells as i64) as usize;
    let a = f + 1.0;
    let c = f - 1.0;
    let d = f - 2.0;
    let w0 = -f * c * d / 6.0;
    let w1 = a * c * d / 2.0;
    let w2 = -a * f * d / 2.0;
    let w3 = a * f * c / 6.0;
    (start, [w0, w1, w2, w3], 4)
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
