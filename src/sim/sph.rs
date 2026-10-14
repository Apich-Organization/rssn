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
//! | Artificial viscosity | Monaghan `alpha`, `beta` with Balsara-free form |
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
//! Limitations: no Balsara switch, no self-gravity.
#![allow(
    clippy::arithmetic_side_effects,
    clippy::assigning_clones,
    clippy::cast_lossless,
    clippy::cast_possible_truncation,
    clippy::chunks_exact_to_as_chunks,
    clippy::cognitive_complexity,
    clippy::derive_partial_eq_without_eq,
    clippy::doc_markdown,
    clippy::excessive_precision,
    clippy::float_cmp,
    clippy::if_not_else,
    clippy::imprecise_flops,
    clippy::indexing_slicing,
    clippy::int_plus_one,
    clippy::items_after_statements,
    clippy::manual_let_else,
    clippy::manual_memcpy,
    clippy::manual_midpoint,
    clippy::many_single_char_names,
    clippy::missing_const_for_fn,
    clippy::module_name_repetitions,
    clippy::needless_pass_by_value,
    clippy::needless_range_loop,
    clippy::option_if_let_else,
    clippy::similar_names,
    clippy::single_call_fn,
    clippy::single_match_else,
    clippy::suboptimal_flops,
    clippy::suspicious_operation_groupings,
    clippy::too_long_first_doc_paragraph,
    clippy::too_many_arguments,
    clippy::too_many_lines,
    clippy::type_complexity,
    clippy::unreadable_literal,
    clippy::useless_let_if_seq,
    clippy::while_float
)]

use std::collections::HashMap;
use std::f64::consts::PI;

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
            Self::CubicSpline => {
                let q = r / h;
                let sigma = match d {
                    1 => 2.0 / (3.0 * h),
                    2 => 10.0 / (7.0 * PI * h * h),
                    _ => 1.0 / (PI * h * h * h),
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
            Self::WendlandC2 => {
                let hh = 2.0 * h;
                let q = r / hh;
                if q >= 1.0 {
                    return 0.0;
                }
                match d {
                    1 => 5.0 / (4.0 * hh) * (1.0 - q).powi(3) * (1.0 + 3.0 * q),
                    2 => 7.0 / (PI * hh * hh) * (1.0 - q).powi(4) * (4.0 * q + 1.0),
                    _ => 21.0 / (2.0 * PI * hh * hh * hh) * (1.0 - q).powi(4) * (4.0 * q + 1.0),
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
            Self::CubicSpline => {
                let q = r / h;
                let sigma = match d {
                    1 => 2.0 / (3.0 * h),
                    2 => 10.0 / (7.0 * PI * h * h),
                    _ => 1.0 / (PI * h * h * h),
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
            Self::WendlandC2 => {
                let hh = 2.0 * h;
                let q = r / hh;
                if q >= 1.0 {
                    return 0.0;
                }
                match d {
                    1 => 5.0 / (4.0 * hh) * (-12.0 * q * (1.0 - q).powi(2)) / hh,
                    2 => 7.0 / (PI * hh * hh) * (-20.0 * q * (1.0 - q).powi(3)) / hh,
                    _ => 21.0 / (2.0 * PI * hh * hh * hh) * (-20.0 * q * (1.0 - q).powi(3)) / hh,
                }
            },
        }
    }
}

/// Boundary treatment of one axis.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SphBoundary {
    /// No boundary.
    Open,
    /// Periodic images.
    Periodic,
    /// Mirror ghost particles (free-slip reflecting wall).
    Reflecting,
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
}

#[derive(Clone, Copy)]
struct Ghost<const D: usize> {
    src: usize,
    pos: [f64; D],
    flip: [f64; D],
}

struct Grid<const D: usize> {
    cell: f64,
    map: HashMap<[i64; D], Vec<usize>>,
}

impl<const D: usize> Grid<D> {
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
        let mut g = Self { cell, map: HashMap::new() };
        for (i, p) in pts.iter().enumerate() {
            let k = g.key(p);
            g.map.entry(k).or_default().push(i);
        }
        g
    }

    fn for_near(
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
        };
        // initial h guess from the bounding volume
        let vol: f64 = (0..D).map(|a| (hi[a] - lo[a]).max(1e-12)).product();
        let h0 = s.eta * (vol / n as f64).powf(1.0 / D as f64);
        s.h = vec![h0; n];
        s.update_density();
        s
    }

    fn ghosts(
        &self,
        reach: f64,
    ) -> Vec<Ghost<D>> {
        let mut all: Vec<Ghost<D>> = (0..self.pos.len()).map(|i| Ghost { src: i, pos: self.pos[i], flip: [1.0; D] }).collect();
        for a in 0..D {
            if self.bc[a] == SphBoundary::Open {
                continue;
            }
            let len = self.hi[a] - self.lo[a];
            let mut extra = Vec::new();
            for g in &all {
                let x = g.pos[a];
                if x - self.lo[a] < reach {
                    let mut n = *g;
                    match self.bc[a] {
                        SphBoundary::Periodic => n.pos[a] = x + len,
                        _ => {
                            n.pos[a] = 2.0 * self.lo[a] - x;
                            n.flip[a] = -g.flip[a];
                        },
                    }
                    extra.push(n);
                }
                if self.hi[a] - x < reach {
                    let mut n = *g;
                    match self.bc[a] {
                        SphBoundary::Periodic => n.pos[a] = x - len,
                        _ => {
                            n.pos[a] = 2.0 * self.hi[a] - x;
                            n.flip[a] = -g.flip[a];
                        },
                    }
                    extra.push(n);
                }
            }
            all.extend(extra);
        }
        all
    }

    /// Recomputes `rho` and `h` by summation with a fixed-point iteration on `h`.
    pub fn update_density(&mut self) {
        let n = self.pos.len();
        for _ in 0..4 {
            let hmax = self.h.iter().copied().fold(0.0, f64::max);
            let gh = self.ghosts(2.0 * hmax);
            let pts: Vec<[f64; D]> = gh.iter().map(|g| g.pos).collect();
            let grid = Grid::new(&pts, 2.0 * hmax);
            let mut new_rho = vec![0.0; n];
            for i in 0..n {
                let hi = self.h[i];
                let mut s = 0.0;
                grid.for_near(&self.pos[i], |j| {
                    let mut r2 = 0.0;
                    for a in 0..D {
                        r2 += (self.pos[i][a] - pts[j][a]).powi(2);
                    }
                    let r = r2.sqrt();
                    if r < 2.0 * hi {
                        s += self.mass[gh[j].src] * self.kernel.w(D, r, hi);
                    }
                });
                new_rho[i] = s;
            }
            let mut change = 0.0_f64;
            for i in 0..n {
                let hn = self.eta * (self.mass[i] / new_rho[i]).powf(1.0 / D as f64);
                change = change.max(((hn - self.h[i]) / self.h[i]).abs());
                self.h[i] = hn;
            }
            self.rho = new_rho;
            if change < 1e-6 {
                break;
            }
        }
        // final density with converged h
        let hmax = self.h.iter().copied().fold(0.0, f64::max);
        let gh = self.ghosts(2.0 * hmax);
        let pts: Vec<[f64; D]> = gh.iter().map(|g| g.pos).collect();
        let grid = Grid::new(&pts, 2.0 * hmax);
        for i in 0..n {
            let hi = self.h[i];
            let mut s = 0.0;
            grid.for_near(&self.pos[i], |j| {
                let mut r2 = 0.0;
                for a in 0..D {
                    r2 += (self.pos[i][a] - pts[j][a]).powi(2);
                }
                let r = r2.sqrt();
                if r < 2.0 * hi {
                    s += self.mass[gh[j].src] * self.kernel.w(D, r, hi);
                }
            });
            self.rho[i] = s;
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
        (self.gamma * (self.gamma - 1.0) * self.u[i]).max(0.0).sqrt()
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
                        let dwdh = -(D as f64 * self.kernel.w(D, r, hi) + r * self.kernel.dw(D, r, hi)) / hi;
                        sum += self.mass[j] * dwdh;
                    }
                });
                omega[i] = 1.0 + self.h[i] / (D as f64 * self.rho[i]) * sum;
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
                let gw = 0.5 * (dwi + dwj);
                let vr: f64 = (0..D).map(|a| vij[a] * rij[a]).sum();
                let mut pi = 0.0;
                if vr < 0.0 {
                    let hb = 0.5 * (hi + hj);
                    let mu = hb * vr / (r2 + 0.01 * hb * hb);
                    let cb = 0.5 * (c[i] + c[j]);
                    let rb = 0.5 * (self.rho[i] + self.rho[j]);
                    pi = (-self.alpha * cb * mu + self.beta * mu * mu) / rb;
                }
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
                let fac = p[i] / (self.rho[i] * self.rho[i]) + p[j] / (self.rho[j] * self.rho[j]) + pi;
                let mut vdotg = 0.0;
                for a in 0..D {
                    let gradw = gw * rij[a] / r;
                    acc[i][a] -= self.mass[j] * fac * gradw;
                    vdotg += vij[a] * gradw;
                }
                du[i] += 0.5 * self.mass[j] * fac * vdotg;
            });
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
                    SphBoundary::Open => {},
                    SphBoundary::Periodic => {
                        let w = ((self.pos[i][a] - self.lo[a]) / len).floor();
                        self.pos[i][a] -= w * len;
                    },
                    SphBoundary::Reflecting => {
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
        (0..self.pos.len()).map(|i| self.mass[i] * (0.5 * self.vel[i].iter().map(|x| x * x).sum::<f64>() + self.u[i])).sum()
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
