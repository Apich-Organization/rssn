//! Smoothed-particle hydrodynamics in 1-3 dimensions.
//!
//! # Smoothed-particle hydrodynamics in 1, 2 and 3 dimensions
//!
//! Compressible, inviscid SPH for an ideal gas with density summation,
//! variable smoothing lengths `h = eta (m / rho)^(1/D)` and Monaghan
//! artificial viscosity.
//!
//! | Feature | Support |
//! |---|---|
//! | Dimensions | any `D` (const generic); 1-3 tested |
//! | Kernels | cubic B-spline, Wendland C2 (support radius `2h`, normalised in every `D`) |
//! | Density | summation `rho_i = sum m_j W(r_ij, h_i)` with fixed-point `h` iteration |
//! | Momentum / energy | symmetric pressure gradient with symmetrised kernel gradient, `du/dt` consistent with `dv/dt` |
//! | Artificial viscosity | Monaghan `alpha`, `beta`, optional **Balsara switch** ([`Sph::balsara`]) that suppresses the viscosity in shear flow |
//! | Boundaries per axis ([`SphBoundary`]) | open, periodic (ghost images), reflecting walls (mirror ghost particles); fixed boundary particles via [`Sph::fixed`] |
//! | Integrator | Heun (RK2) with CFL helper |
//! | Neighbour search | uniform cell grid, `O(N)` |
//!
//! The variational **grad-h** equations of motion (Springel & Hernquist
//! 2002) are available through [`Sph::grad_h`]; they include the factor
//! `Ω_i = 1 + h_i/(D ρ_i) Σ_j m_j ∂W_ij/∂h_i` and conserve energy and
//! momentum with variable smoothing lengths. Dissipative particle dynamics
//! lives in [`crate::sim::dpd`].
//!
//! Self-gravity: set [`Sph::gravity`] to a [`SphGravity`]; the accelerations of the Barnes-Hut
//! tree code [`crate::sim::nbody::gravity_accelerations`] (Plummer-softened, parallel) are added
//! to the hydrodynamic ones. [`Sph::gravitational_energy`] gives the matching 3D pair energy.
//! The gravitational softening is a fixed `eps`, independent of the smoothing lengths, so the
//! force is not the variationally consistent adaptive-softening one (Price & Monaghan 2007):
//! energy is conserved to the tree's accuracy and the time integration, not to round-off.
//!
//! Limitations: gravity is the free-space (non-periodic) tree force, with no ghost images for
//! periodic or reflecting boundaries.

use std::collections::HashMap;
use std::f64::consts::PI;

use crate::kernels::numeric_common::scale_product;

/// Smoothing kernel family (support radius `2h`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum KernelKind {
    /// Cubic B-spline.
    CubicSpline,
    /// Wendland C2.
    WendlandC2,
}

impl KernelKind {
    /// Kernel value `W(r, h)` in `d` dimensions.
    #[must_use]
    pub fn w(
        self,
        d: usize,
        r: f64,
        h: f64,
    ) -> f64 {
        match self {
            | Self::CubicSpline => {
                let q = r / h;
                let sigma = match d {
                    | 1 => 2.0 / (3.0 * h),
                    | 2 => 10.0 / (7.0 * PI * h * h),
                    | _ => 1.0 / (PI * h * h * h),
                };
                sigma
                    * if q < 1.0 {
                        1.0 - 1.5 * q * q + 0.75 * q * q * q
                    } else if q < 2.0 {
                        0.25 * (2.0 - q).powi(3)
                    } else {
                        0.0
                    }
            },
            | Self::WendlandC2 => {
                let hh = 2.0 * h;
                let q = r / hh;
                if q >= 1.0 {
                    return 0.0;
                }
                match d {
                    | 1 => 5.0 / (4.0 * hh) * (1.0 - q).powi(3) * (1.0 + 3.0 * q),
                    | 2 => 7.0 / (PI * hh * hh) * (1.0 - q).powi(4) * (4.0 * q + 1.0),
                    | _ => 21.0 / (2.0 * PI * hh * hh * hh) * (1.0 - q).powi(4) * (4.0 * q + 1.0),
                }
            },
        }
    }

    /// Radial derivative `dW/dr`.
    #[must_use]
    pub fn dw(
        self,
        d: usize,
        r: f64,
        h: f64,
    ) -> f64 {
        match self {
            | Self::CubicSpline => {
                let q = r / h;
                let sigma = match d {
                    | 1 => 2.0 / (3.0 * h),
                    | 2 => 10.0 / (7.0 * PI * h * h),
                    | _ => 1.0 / (PI * h * h * h),
                };
                sigma
                    * if q < 1.0 {
                        -3.0 * q + 2.25 * q * q
                    } else if q < 2.0 {
                        -0.75 * (2.0 - q).powi(2)
                    } else {
                        0.0
                    }
                    / h
            },
            | Self::WendlandC2 => {
                let hh = 2.0 * h;
                let q = r / hh;
                if q >= 1.0 {
                    return 0.0;
                }
                match d {
                    | 1 => 5.0 / (4.0 * hh) * (-12.0 * q * (1.0 - q).powi(2)) / hh,
                    | 2 => 7.0 / (PI * hh * hh) * (-20.0 * q * (1.0 - q).powi(3)) / hh,
                    | _ => 21.0 / (2.0 * PI * hh * hh * hh) * (-20.0 * q * (1.0 - q).powi(3)) / hh,
                }
            },
        }
    }
}

/// Boundary treatment of one axis.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SphBoundary {
    /// No boundary.
    Open,
    /// Periodic images.
    Periodic,
    /// Mirror ghost particles (free-slip reflecting wall).
    Reflecting,
}

/// Self-gravity parameters of [`Sph`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SphGravity {
    /// Gravitational constant `G`.
    pub g: f64,
    /// Plummer softening length.
    pub eps: f64,
    /// Barnes-Hut opening angle.
    pub theta: f64,
}

/// SPH particle system in `D` dimensions.
#[derive(Clone, Debug)]
pub struct Sph<const D: usize> {
    /// Positions.
    pub pos: Vec<[f64; D]>,
    /// Velocities.
    pub vel: Vec<[f64; D]>,
    /// Particle masses.
    pub mass: Vec<f64>,
    /// Specific internal energies.
    pub u: Vec<f64>,
    /// Densities (updated by [`Sph::update_density`]).
    pub rho: Vec<f64>,
    /// Smoothing lengths.
    pub h: Vec<f64>,
    /// Frozen boundary particles (do not move, still contribute).
    pub fixed: Vec<bool>,
    /// Adiabatic index.
    pub gamma: f64,
    /// Smoothing length factor `eta`.
    pub eta: f64,
    /// Linear artificial viscosity coefficient.
    pub alpha: f64,
    /// Quadratic artificial viscosity coefficient.
    pub beta: f64,
    /// Kernel.
    pub kernel: KernelKind,
    /// Domain lower corner.
    pub lo: [f64; D],
    /// Domain upper corner.
    pub hi: [f64; D],
    /// Boundary per axis.
    pub bc: [SphBoundary; D],
    /// Use the grad-h (Springel–Hernquist) corrected, variational equations of motion.
    pub grad_h: bool,
    /// Balsara (1995) shear switch `f = |∇·v| / (|∇·v| + |∇×v| + 10⁻⁴ c / h)` on the artificial viscosity.
    pub balsara: bool,
    /// Self-gravity (`None`: pure hydrodynamics).
    pub gravity: Option<SphGravity>,
}

#[derive(Clone, Copy)]
struct Ghost<const D: usize> {
    src: usize,
    pos: [f64; D],
    flip: [f64; D],
}

fn make_ghosts<const D: usize>(
    pos: &[[f64; D]],
    lo: [f64; D],
    hi: [f64; D],
    bc: [SphBoundary; D],
    reach: f64,
) -> Vec<Ghost<D>> {
    let mut all: Vec<_> = pos
        .iter()
        .enumerate()
        .map(|(src, &pos)| {
            Ghost {
                src,
                pos,
                flip: [1.0; D],
            }
        })
        .collect();
    for a in 0..D {
        if bc[a] == SphBoundary::Open {
            continue;
        }
        let original = all.len();
        let len = hi[a] - lo[a];
        for i in 0..original {
            let g = all[i];
            if g.pos[a] - lo[a] < reach {
                let mut near = g;
                if bc[a] == SphBoundary::Periodic {
                    near.pos[a] += len;
                } else {
                    near.pos[a] = 2.0 * lo[a] - g.pos[a];
                    near.flip[a] = -g.flip[a];
                }
                all.push(near);
            }
            if hi[a] - g.pos[a] < reach {
                let mut far = g;
                if bc[a] == SphBoundary::Periodic {
                    far.pos[a] -= len;
                } else {
                    far.pos[a] = 2.0 * hi[a] - g.pos[a];
                    far.flip[a] = -g.flip[a];
                }
                all.push(far);
            }
        }
    }
    all
}

/// Cell grid for geometric candidate lookups at a bounded distance.
pub(crate) struct Grid<const D: usize> {
    cell: f64,
    map: HashMap<[i64; D], Vec<usize>>,
}

impl<const D: usize> Grid<D> {
    /// Indexes near points when the coordinates can be binned with room for rounding.
    pub(crate) fn for_reach(
        pts: &[[f64; D]],
        reach: f64,
    ) -> Option<Self> {
        let width = reach * 1.125;
        (D <= 3
            && width.is_finite()
            && pts
                .iter()
                .flatten()
                .all(|&d| (d / width).abs() < (1_u64 << 47) as f64))
        .then(|| Self::new(pts, width))
    }

    fn key(
        &self,
        p: &[f64; D],
    ) -> [i64; D] {
        let mut k = [0i64; D];
        for a in 0..D {
            k[a] = (p[a] / self.cell).floor() as i64;
        }
        k
    }

    fn new(
        pts: &[[f64; D]],
        cell: f64,
    ) -> Self {
        let mut g = Self {
            cell,
            map: HashMap::new(),
        };
        for (i, p) in pts.iter().enumerate() {
            let k = g.key(p);
            g.map.entry(k).or_default().push(i);
        }
        g
    }

    /// Visits every point from the query cell and its adjacent cells.
    pub(crate) fn for_near(
        &self,
        p: &[f64; D],
        mut f: impl FnMut(usize),
    ) {
        let k = self.key(p);
        let count = 3usize.pow(D as u32);
        for mut r in 0..count {
            let mut kk = k;
            for a in 0..D {
                kk[a] += (r % 3) as i64 - 1;
                r /= 3;
            }
            if let Some(v) = self.map.get(&kk) {
                for &j in v {
                    f(j);
                }
            }
        }
    }
}

struct DensityCells<const D: usize> {
    ghosts: Vec<Ghost<D>>,
    origin: [f64; D],
    grid: Option<Grid<D>>,
}

impl<const D: usize> DensityCells<D> {
    fn visit(
        &self,
        index: usize,
        reach: f64,
        mut f: impl FnMut(usize, f64),
    ) {
        let pos = self.ghosts[index].pos;
        let reach_squared = reach * reach;
        let visit = |j: usize| {
            let ghost = &self.ghosts[j];
            let delta: [f64; D] = std::array::from_fn(|a| pos[a] - ghost.pos[a]);
            let r2: f64 = delta.iter().map(|d| d * d).sum();
            if reach_squared.is_normal() && r2 > reach_squared {
                return;
            }
            let r = if r2.is_normal() || delta.iter().all(|&d| d == 0.0) {
                r2.sqrt()
            } else {
                delta.into_iter().fold(0.0, f64::hypot)
            };
            if r < reach {
                f(ghost.src, r);
            }
        };
        if let Some(grid) = &self.grid {
            let coords = std::array::from_fn(|a| pos[a] - self.origin[a]);
            grid.for_near(&coords, visit);
        } else {
            (0..self.ghosts.len()).for_each(visit);
        }
    }
}

/// Reusable SPH density summation for fixed particle positions.
///
/// Periodic positions are wrapped into the primary cell. The sum includes
/// the particle itself and one image per side of each periodic or reflecting
/// axis, including their combinations at corners. On a periodic axis a
/// support radius up to the period includes every contributing image.
pub struct SphDensityPlan<const D: usize> {
    n: usize,
    max_h: f64,
    offsets: Vec<usize>,
    neighbors: Vec<(usize, f64)>,
    remaining: Option<DensityCells<D>>,
}

impl<const D: usize> SphDensityPlan<D> {
    pub(crate) fn cached_heap_bytes(&self) -> Option<usize> {
        self.remaining.is_none().then(|| {
            self.offsets.capacity() * size_of::<usize>()
                + self.neighbors.capacity() * size_of::<(usize, f64)>()
        })
    }

    /// Prepares a density plan with smoothing lengths up to `max_h`.
    ///
    /// # Errors
    ///
    /// Unsupported dimension, non-finite positions, invalid smoothing
    /// length or invalid periodic or reflecting extents.
    pub fn new(
        pos: &[[f64; D]],
        lo: [f64; D],
        hi: [f64; D],
        bc: [SphBoundary; D],
        max_h: f64,
    ) -> Result<Self, String> {
        if !(1..=3).contains(&D) {
            return Err("SPH density requires one to three dimensions".into());
        }
        Self::build(pos, lo, hi, bc, max_h)
    }

    fn build(
        pos: &[[f64; D]],
        mut lo: [f64; D],
        mut hi: [f64; D],
        bc: [SphBoundary; D],
        max_h: f64,
    ) -> Result<Self, String> {
        const MAX_RECORDS: usize = 524_288;
        let reach = 2.0 * max_h;
        if reach <= 0.0 || !reach.is_finite() || pos.iter().flatten().any(|x| !x.is_finite()) {
            return Err("invalid SPH density coordinates or smoothing limit".into());
        }
        let mut positions = pos.to_vec();
        for a in 0..D {
            if bc[a] == SphBoundary::Open {
                continue;
            }
            let length = hi[a] - lo[a];
            if !length.is_finite() || length <= 0.0 {
                return Err("invalid SPH density boundary extent".into());
            }
            for point in &mut positions {
                let x = point[a] - lo[a];
                point[a] = if bc[a] == SphBoundary::Periodic {
                    let wrapped = if (0.0..length).contains(&x) {
                        x
                    } else if x.is_finite() {
                        x.rem_euclid(length)
                    } else {
                        (point[a].rem_euclid(length) - lo[a].rem_euclid(length)).rem_euclid(length)
                    };
                    if wrapped >= length {
                        0.0
                    } else {
                        wrapped
                    }
                } else {
                    x
                };
            }
            lo[a] = 0.0;
            hi[a] = length;
        }
        let ghosts = make_ghosts(&positions, lo, hi, bc, reach);
        let origin = positions.first().copied().unwrap_or([0.0; D]);
        let relative: Vec<[f64; D]> = ghosts
            .iter()
            .map(|g| std::array::from_fn(|a| g.pos[a] - origin[a]))
            .collect();
        let grid = Grid::for_reach(&relative, reach);
        let cells = DensityCells { ghosts, origin, grid };
        let n = pos.len();
        let mut offsets = Vec::with_capacity(n + 1);
        offsets.push(0);
        let mut neighbors = Vec::with_capacity(n.saturating_mul(16).min(MAX_RECORDS));
        for i in 0..n {
            let start = neighbors.len();
            let mut full = false;
            cells.visit(i, reach, |j, distance| {
                if neighbors.len() < MAX_RECORDS {
                    neighbors.push((j, distance));
                } else {
                    full = true;
                }
            });
            if full {
                neighbors.truncate(start);
                break;
            }
            offsets.push(neighbors.len());
        }
        let remaining = (offsets.len() <= n).then_some(cells);
        Ok(Self {
            n,
            max_h,
            offsets,
            neighbors,
            remaining,
        })
    }

    /// Evaluates the mass density for a new mass, smoothing and kernel field.
    ///
    /// # Errors
    ///
    /// Wrong field lengths, negative or non-finite masses, invalid smoothing
    /// lengths or an unrepresentable density.
    pub fn evaluate(
        &self,
        mass: &[f64],
        smoothing: &[f64],
        kernel: KernelKind,
    ) -> Result<Vec<f64>, String> {
        if mass.len() != self.n
            || smoothing.len() != self.n
            || mass.iter().any(|&m| !m.is_finite() || m < 0.0)
            || smoothing
                .iter()
                .any(|&h| !h.is_finite() || h <= 0.0 || h > self.max_h)
        {
            return Err("invalid SPH density masses or smoothing lengths".into());
        }
        let mut rho = vec![0.0; self.n];
        self.write_density(mass, smoothing, kernel, &mut rho);
        if rho.iter().any(|x| !x.is_finite()) {
            return Err("SPH density is not representable".into());
        }
        Ok(rho)
    }

    fn write_density(
        &self,
        mass: &[f64],
        smoothing: &[f64],
        kernel: KernelKind,
        rho: &mut [f64],
    ) {
        crate::kernels::par::for_each_mut(rho, 512, |i, density| {
            let h = smoothing[i];
            let reach = 2.0 * h;
            let mut sum = 0.0;
            let mut add = |j: usize, distance: f64| {
                let m = mass[j];
                if m == 0.0 || distance >= reach {
                    return;
                }
                let weight = kernel.w(D, distance, h);
                sum += if weight.is_normal() {
                    m * weight
                } else {
                    scale_product(
                        kernel.w(D, distance / h, 1.0) * m.signum(),
                        &[m.abs()],
                        &[h; 3][..D.clamp(1, 3)],
                    )
                    .unwrap_or(f64::INFINITY)
                };
            };
            if let Some(&end) = self.offsets.get(i + 1) {
                for &(j, distance) in &self.neighbors[self.offsets[i]..end] {
                    add(j, distance);
                }
            } else if let Some(cells) = &self.remaining {
                cells.visit(i, reach, add);
            }
            *density = sum;
        });
    }
}

fn sph_root(
    value: f64,
    dimension: usize,
) -> f64 {
    match dimension {
        | 1 => value,
        | 2 => value.sqrt(),
        | 3 => value.cbrt(),
        | _ => value.powf(1.0 / dimension as f64),
    }
}

impl<const D: usize> Sph<D> {
    /// Creates a gas of equal-mass particles at rest with given positions and specific energy.
    #[must_use]
    pub fn new(
        pos: Vec<[f64; D]>,
        mass: Vec<f64>,
        u: Vec<f64>,
        lo: [f64; D],
        hi: [f64; D],
        bc: [SphBoundary; D],
    ) -> Self {
        let n = pos.len();
        let mut s = Self {
            vel: vec![[0.0; D]; n],
            rho: vec![1.0; n],
            h: vec![1.0; n],
            fixed: vec![false; n],
            pos,
            mass,
            u,
            gamma: 1.4,
            eta: 1.2,
            alpha: 1.0,
            beta: 2.0,
            kernel: KernelKind::CubicSpline,
            lo,
            hi,
            bc,
            grad_h: false,
            balsara: false,
            gravity: None,
        };
        // Estimate h from the domain volume per particle.
        let lengths: [f64; D] = std::array::from_fn(|a| {
            let length = hi[a] - lo[a];
            if length.is_finite() && length > 0.0 {
                length
            } else {
                1e-12
            }
        });
        let volume: f64 = lengths.iter().product();
        let per_particle = volume / n as f64;
        let direct = s.eta * per_particle.powf(1.0 / D as f64);
        let h0 = if n == 0 || per_particle.is_normal() || !(1..=3).contains(&D) {
            direct
        } else {
            let roots = lengths.map(|length| sph_root(length, D));
            scale_product(s.eta, &roots, &[sph_root(n as f64, D)]).unwrap_or(direct)
        };
        s.h.fill(h0);
        s.update_density();
        s
    }

    fn ghosts(
        &self,
        reach: f64,
    ) -> Vec<Ghost<D>> {
        make_ghosts(&self.pos, self.lo, self.hi, self.bc, reach)
    }

    /// Recomputes `rho` and `h` by summation with a fixed-point iteration on `h`.
    pub fn update_density(&mut self) {
        self.rho.resize(self.pos.len(), 0.0);
        if self.pos.is_empty() {
            return;
        }
        let mut plan: Option<SphDensityPlan<D>> = None;
        let mut converged = false;
        for iteration in 0..=4 {
            let max_h = self.h.iter().copied().fold(0.0, f64::max);
            if plan.as_ref().is_none_or(|plan| plan.max_h < max_h) {
                let horizon = (1.125 * max_h).min(0.5 * f64::MAX);
                plan = SphDensityPlan::build(&self.pos, self.lo, self.hi, self.bc, horizon).ok();
            }
            let Some(plan) = &plan else {
                self.rho.fill(f64::NAN);
                return;
            };
            plan.write_density(&self.mass, &self.h, self.kernel, &mut self.rho);
            if converged || iteration == 4 {
                break;
            }
            let mut change = 0.0_f64;
            for ((h, &m), &rho) in self.h.iter_mut().zip(&self.mass).zip(&self.rho) {
                let ratio = m / rho;
                let mut next = self.eta * ratio.powf(1.0 / D as f64);
                if !ratio.is_normal() && m > 0.0 && rho > 0.0 && self.eta > 0.0 {
                    next = scale_product(self.eta, &[sph_root(m, D)], &[sph_root(rho, D)])
                        .unwrap_or(next);
                }
                change = change.max(((next - *h) / *h).abs());
                *h = next;
            }
            converged = change < 1e-6;
        }
    }

    /// Pressure of particle `i`.
    #[must_use]
    pub fn pressure(
        &self,
        i: usize,
    ) -> f64 {
        (self.gamma - 1.0) * self.rho[i] * self.u[i]
    }

    /// Sound speed of particle `i`.
    #[must_use]
    pub fn sound_speed(
        &self,
        i: usize,
    ) -> f64 {
        (self.gamma * (self.gamma - 1.0) * self.u[i])
            .max(0.0)
            .sqrt()
    }

    /// Accelerations and energy rates for the current `(pos, vel, u)` (densities are refreshed).
    pub fn derivatives(&mut self) -> (Vec<[f64; D]>, Vec<f64>) {
        self.update_density();
        let n = self.pos.len();
        let hmax = self.h.iter().copied().fold(0.0, f64::max);
        let gh = self.ghosts(2.0 * hmax);
        let pts: Vec<[f64; D]> = gh.iter().map(|g| g.pos).collect();
        let grid = Grid::new(&pts, 2.0 * hmax);
        let mut acc = vec![[0.0; D]; n];
        let mut du = vec![0.0; n];
        let p: Vec<f64> = (0..n).map(|i| self.pressure(i)).collect();
        let c: Vec<f64> = (0..n).map(|i| self.sound_speed(i)).collect();
        // grad-h factors Omega_i = 1 + h_i / (D rho_i) * sum_j m_j dW_ij/dh_i
        let mut omega = vec![1.0; n];
        if self.grad_h {
            for i in 0..n {
                let mut sum = 0.0;
                grid.for_near(&self.pos[i], |jj| {
                    let j = gh[jj].src;
                    let mut r2 = 0.0;
                    for a in 0..D {
                        r2 += (self.pos[i][a] - gh[jj].pos[a]).powi(2);
                    }
                    let r = r2.sqrt();
                    let hi = self.h[i];
                    if r < 2.0 * hi {
                        let dwdh = -(D as f64 * self.kernel.w(D, r, hi)
                            + r * self.kernel.dw(D, r, hi))
                            / hi;
                        sum += self.mass[j] * dwdh;
                    }
                });
                omega[i] = 1.0 + self.h[i] / (D as f64 * self.rho[i]) * sum;
            }
        }
        // Balsara switch from the SPH estimates of the divergence and curl of the velocity
        let mut fbal = vec![1.0; n];
        if self.balsara {
            for i in 0..n {
                let (mut div, mut curl) = (0.0, [0.0; 3]);
                grid.for_near(&self.pos[i], |jj| {
                    let g = &gh[jj];
                    let j = g.src;
                    if jj == i {
                        return;
                    }
                    let mut rij = [0.0; D];
                    let mut vij = [0.0; D];
                    let mut r2 = 0.0;
                    for a in 0..D {
                        rij[a] = self.pos[i][a] - g.pos[a];
                        vij[a] = self.vel[i][a] - g.flip[a] * self.vel[j][a];
                        r2 += rij[a] * rij[a];
                    }
                    let r = r2.sqrt();
                    let hi = self.h[i];
                    if r >= 2.0 * hi || r == 0.0 {
                        return;
                    }
                    let dw = self.kernel.dw(D, r, hi);
                    let mut gradw = [0.0; 3];
                    let mut v3 = [0.0; 3];
                    for a in 0..D {
                        gradw[a] = dw * rij[a] / r;
                        v3[a] = vij[a];
                    }
                    div -= self.mass[j] * (0..D).map(|a| v3[a] * gradw[a]).sum::<f64>();
                    curl[0] += self.mass[j] * (v3[1] * gradw[2] - v3[2] * gradw[1]);
                    curl[1] += self.mass[j] * (v3[2] * gradw[0] - v3[0] * gradw[2]);
                    curl[2] += self.mass[j] * (v3[0] * gradw[1] - v3[1] * gradw[0]);
                });
                let rho = self.rho[i];
                let divv = (div / rho).abs();
                let curlv =
                    (curl[0] * curl[0] + curl[1] * curl[1] + curl[2] * curl[2]).sqrt() / rho;
                fbal[i] = divv / (divv + curlv + 1e-4 * c[i] / self.h[i]);
            }
        }
        for i in 0..n {
            grid.for_near(&self.pos[i], |jj| {
                let g = &gh[jj];
                let j = g.src;
                if jj == i {
                    return;
                }
                let mut rij = [0.0; D];
                let mut vij = [0.0; D];
                let mut r2 = 0.0;
                for a in 0..D {
                    rij[a] = self.pos[i][a] - g.pos[a];
                    vij[a] = self.vel[i][a] - g.flip[a] * self.vel[j][a];
                    r2 += rij[a] * rij[a];
                }
                let r = r2.sqrt();
                let (hi, hj) = (self.h[i], self.h[j]);
                if r >= 2.0 * hi.max(hj) || r == 0.0 {
                    return;
                }
                let (dwi, dwj) = (self.kernel.dw(D, r, hi), self.kernel.dw(D, r, hj));
                let gw = f64::midpoint(dwi, dwj);
                let vr: f64 = (0..D).map(|a| vij[a] * rij[a]).sum();
                let pi = if vr < 0.0 {
                    let hb = f64::midpoint(hi, hj);
                    let mu = hb * vr / (r2 + 0.01 * hb * hb);
                    let cb = f64::midpoint(c[i], c[j]);
                    let rb = f64::midpoint(self.rho[i], self.rho[j]);
                    (-self.alpha * cb * mu + self.beta * mu * mu) / rb
                        * f64::midpoint(fbal[i], fbal[j])
                } else {
                    0.0
                };
                if self.grad_h {
                    let ti = p[i] / (omega[i] * self.rho[i] * self.rho[i]);
                    let tj = p[j] / (omega[j] * self.rho[j] * self.rho[j]);
                    let mut vdot_i = 0.0;
                    let mut vdot_avg = 0.0;
                    for a in 0..D {
                        let n_a = rij[a] / r;
                        acc[i][a] -= self.mass[j] * (ti * dwi + tj * dwj + pi * gw) * n_a;
                        vdot_i += vij[a] * dwi * n_a;
                        vdot_avg += vij[a] * gw * n_a;
                    }
                    du[i] += self.mass[j] * ti * vdot_i + 0.5 * self.mass[j] * pi * vdot_avg;
                    return;
                }
                let fac =
                    p[i] / (self.rho[i] * self.rho[i]) + p[j] / (self.rho[j] * self.rho[j]) + pi;
                let mut vdotg = 0.0;
                for a in 0..D {
                    let gradw = gw * rij[a] / r;
                    acc[i][a] -= self.mass[j] * fac * gradw;
                    vdotg += vij[a] * gradw;
                }
                du[i] += 0.5 * self.mass[j] * fac * vdotg;
            });
        }
        if let Some(gr) = self.gravity
            && let Ok(ag) = crate::sim::nbody::gravity_accelerations(
                &self.pos, &self.mass, gr.g, gr.eps, gr.theta,
            )
        {
            for (a, g) in acc.iter_mut().zip(&ag) {
                for k in 0..D {
                    a[k] += g[k];
                }
            }
        }
        for i in 0..n {
            if self.fixed[i] {
                acc[i] = [0.0; D];
                du[i] = 0.0;
            }
        }
        (acc, du)
    }

    /// Advances one Heun (RK2) step of size `dt`.
    pub fn step(
        &mut self,
        dt: f64,
    ) {
        let n = self.pos.len();
        let (x0, v0, u0) = (self.pos.clone(), self.vel.clone(), self.u.clone());
        let (a1, du1) = self.derivatives();
        for i in 0..n {
            if self.fixed[i] {
                continue;
            }
            for a in 0..D {
                self.pos[i][a] = x0[i][a] + dt * v0[i][a];
                self.vel[i][a] = v0[i][a] + dt * a1[i][a];
            }
            self.u[i] = (u0[i] + dt * du1[i]).max(1e-12);
        }
        self.reflect_and_wrap();
        let (a2, du2) = self.derivatives();
        for i in 0..n {
            if self.fixed[i] {
                continue;
            }
            for a in 0..D {
                self.pos[i][a] = x0[i][a] + 0.5 * dt * (v0[i][a] + self.vel[i][a]);
                self.vel[i][a] = v0[i][a] + 0.5 * dt * (a1[i][a] + a2[i][a]);
            }
            self.u[i] = (u0[i] + 0.5 * dt * (du1[i] + du2[i])).max(1e-12);
        }
        self.reflect_and_wrap();
        self.update_density();
    }

    fn reflect_and_wrap(&mut self) {
        for i in 0..self.pos.len() {
            for a in 0..D {
                let len = self.hi[a] - self.lo[a];
                match self.bc[a] {
                    | SphBoundary::Open => {},
                    | SphBoundary::Periodic => {
                        let w = ((self.pos[i][a] - self.lo[a]) / len).floor();
                        self.pos[i][a] -= w * len;
                    },
                    | SphBoundary::Reflecting => {
                        if self.pos[i][a] < self.lo[a] {
                            self.pos[i][a] = 2.0 * self.lo[a] - self.pos[i][a];
                            self.vel[i][a] = -self.vel[i][a];
                        } else if self.pos[i][a] > self.hi[a] {
                            self.pos[i][a] = 2.0 * self.hi[a] - self.pos[i][a];
                            self.vel[i][a] = -self.vel[i][a];
                        }
                    },
                }
            }
        }
    }

    /// CFL-limited time step `cfl * min_i h_i / (c_i + |v_i| + viscous term)`.
    #[must_use]
    pub fn suggest_dt(
        &self,
        cfl: f64,
    ) -> f64 {
        let mut dt = f64::INFINITY;
        for i in 0..self.pos.len() {
            let v: f64 = self.vel[i].iter().map(|x| x * x).sum::<f64>().sqrt();
            let s = self.sound_speed(i) * (1.0 + 1.2 * self.alpha) + v;
            dt = dt.min(cfl * self.h[i] / s.max(1e-12));
        }
        dt
    }

    /// Total kinetic plus thermal energy.
    #[must_use]
    pub fn total_energy(&self) -> f64 {
        (0..self.pos.len())
            .map(|i| {
                self.mass[i] * (0.5 * self.vel[i].iter().map(|x| x * x).sum::<f64>() + self.u[i])
            })
            .sum()
    }

    /// Plummer-softened pair energy `-G sum_{i<j} m_i m_j / sqrt(r^2 + eps^2)` of the self-gravity
    /// (`O(N^2)`; three dimensions only; `None` without gravity or for `D != 3`).
    #[must_use]
    pub fn gravitational_energy(&self) -> Option<f64> {
        let gr = self.gravity?;
        if D != 3 {
            return None;
        }
        let n = self.pos.len();
        let mut e = 0.0;
        for i in 0..n {
            for j in i + 1..n {
                let r2: f64 = (0..D)
                    .map(|a| (self.pos[i][a] - self.pos[j][a]).powi(2))
                    .sum();
                e -= gr.g * self.mass[i] * self.mass[j] / (r2 + gr.eps * gr.eps).sqrt();
            }
        }
        Some(e)
    }

    /// Total momentum.
    #[must_use]
    pub fn momentum(&self) -> [f64; D] {
        let mut p = [0.0; D];
        for i in 0..self.pos.len() {
            for a in 0..D {
                p[a] += self.mass[i] * self.vel[i][a];
            }
        }
        p
    }
}
