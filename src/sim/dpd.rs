//! # Dissipative particle dynamics (DPD) in 2D and 3D
//!
//! Groot–Warren DPD fluid in a periodic box or a walled slab: soft conservative repulsion
//! `F^C = a (1 - r/r_c) r̂`, pairwise dissipative friction
//! `F^D = -γ w²(r) (r̂·v_ij) r̂` and random force
//! `F^R = σ w(r) ξ_ij Δt^{-1/2} r̂` with `w = 1 - r/r_c`, `σ² = 2 γ k_B T`
//! (fluctuation–dissipation), integrated by the modified velocity-Verlet
//! scheme (`λ = 0.65`). Every pair force is antisymmetric, so momentum is
//! conserved to round-off, and the thermostat is local and Galilean
//! invariant.
//!
//! | Feature | Support |
//! |---|---|
//! | Dimensions | 2 and 3 |
//! | Boundary | periodic box (minimum image) or one walled slab, cell-list neighbour search |
//! | Diagnostics | kinetic temperature, momentum, virial pressure |
//! | Equation of state | `p = ρ k_B T + α a ρ²`; measured `α = 0.0918 ± 0.0003` at `a = 25`, `ρ = 3` (3D, `r_c = 1`, `Δt` from 0.01 to 0.04, boxes 6 and 8) |
//!
//! **Pressure check.** Groot and Warren quote the fit `α ≈ 0.101` (mean field
//! `π/30 = 0.1047`) for `ρ ≳ 2`. The code was checked for a virial, cut-off or
//! weight-function error by computing the same pressure in an independent
//! way: a Metropolis Monte Carlo simulation of the pure conservative soft
//! potential `U = a r_c (1 - r/r_c)² / 2` (343 particles, `ρ = 3`, the
//! equilibrium measure of the DPD thermostat) gives `α = 0.0918`, identical
//! to the dynamics (`p = 23.65`). The virial (`Σ a w r / (d V)`), the cell
//! list, the weight `w = 1 - r/r_c` and the temperature (`T = 0.9999`) are
//! therefore correct; the correlation reduction from the mean-field value
//! `0.1047` at `a = 25` is larger than the fit value of the paper
//! (`p = 25.7`), so that single number should not be used as a reference
//! for `a = 25`, `ρ = 3` (the fit is accurate for larger `a ρ`).
//!
//! **Walls.** [`Dpd::set_walls`] replaces the periodic boundary along one axis by a slab between
//! two planes at `0` and `L`, with an effective soft repulsion `a_wall (1 - d/r_c)` ([`DpdWall`])
//! and a reflection rule for particles that cross ([`WallKind`]): specular (free slip),
//! bounce-back (the whole velocity reverses) or [`WallKind::NoSlip`] (specular plus a Langevin
//! tangential collision, so a drift decays and the wall stays at `k_B T`). The wall layer has no
//! structure of its own (no frozen particles, no density matching), so the fluid density next to a
//! wall deviates from the bulk within about `r_c`. The wall force enters the virial, so
//! [`Dpd::pressure`] stays consistent in a walled system, and the collision is applied to the
//! pre-collision half of the velocity-Verlet step, so a specular wall does not pump energy in.
//!
//! Limitations: single species and no bonds (polymers). The wall is a boundary rule and not a
//! material: there is no frozen, structured or thermostatted layer of wall particles and no wall
//! density profile, so the density next to a wall is whatever the soft repulsion produces.
//! Exactly one slab is supported — one walled axis with two planes — so walls on two axes at once,
//! more than two planes, and moving walls are not available, and the wall carries no velocity, so
//! [`WallKind::NoSlip`] can damp a drift but cannot drive a shear flow. [`WallKind::Specular`] and
//! [`WallKind::BounceBack`] exchange no energy with the fluid beyond the elastic collision, and
//! only [`WallKind::NoSlip`] thermostats at the wall. [`Dpd::pressure`] returns the trace over `d`
//! and does not resolve the normal and tangential components separately, which differ in a slab.

use crate::kernels::random::Rng;
use std::collections::HashMap;

/// DPD parameters.
#[derive(Clone, Copy, Debug)]
pub struct DpdParams {
    /// Dimension (2 or 3).
    pub dim: usize,
    /// Periodic box edge lengths.
    pub box_len: [f64; 3],
    /// Conservative repulsion `a`.
    pub a: f64,
    /// Friction coefficient `γ`.
    pub gamma: f64,
    /// Temperature `k_B T`.
    pub kt: f64,
    /// Cut-off radius.
    pub rc: f64,
    /// Particle mass.
    pub mass: f64,
}

/// How a particle that reaches a wall is reflected.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum WallKind {
    /// Specular reflection: the normal velocity changes sign and the tangential velocity is
    /// untouched (free slip).
    Specular,
    /// Bounce-back: the whole velocity changes sign (no slip at the wall).
    BounceBack,
    /// Specular reflection plus a Langevin tangential collision: the normal velocity changes
    /// sign and each tangential component becomes `(1 - friction) v_t + ξ` with
    /// `⟨ξ²⟩ = friction (2 - friction) k_B T / m`, so a uniform drift decays while a Maxwell
    /// distribution at `k_B T` is left invariant. `friction = 0` is [`WallKind::Specular`] and
    /// `friction = 1` replaces the tangential velocity by a fresh thermal draw at the wall.
    NoSlip {
        /// Fraction of the tangential velocity removed at each collision, in `[0, 1]`.
        friction: f64,
    },
}

/// A pair of parallel walls at `0` and `box_len[axis]` normal to `axis` (see [`Dpd::set_walls`]).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DpdWall {
    /// Axis the walls are normal to (the box is not periodic along it).
    pub axis: usize,
    /// Effective wall repulsion: a particle at distance `d < r_c` from a wall feels the normal
    /// force `a_wall (1 - d / r_c)` (the soft DPD potential integrated over a frozen wall
    /// layer, up to the wall density); `0` gives a purely reflecting wall.
    pub a_wall: f64,
    /// Reflection rule for particles that cross a wall.
    pub kind: WallKind,
}

/// DPD particle system.
pub struct Dpd {
    /// Parameters.
    pub params: DpdParams,
    /// Positions.
    pub pos: Vec<[f64; 3]>,
    /// Velocities.
    pub vel: Vec<[f64; 3]>,
    force: Vec<[f64; 3]>,
    rng: Rng,
    /// Virial of the conservative force from the last force evaluation.
    virial: f64,
    /// Elapsed time.
    pub t: f64,
    wall: Option<DpdWall>,
    /// Particles that crossed a wall (an odd number of times) in the last position update.
    crossed: Vec<bool>,
}

impl Dpd {
    /// Random uniform positions and Maxwell velocities with zero total momentum.
    ///
    /// # Errors
    /// For a bad dimension, a box smaller than three cut-offs, or fewer than two particles.
    pub fn new(
        n: usize,
        params: DpdParams,
        seed: u64,
    ) -> Result<Self, String> {
        if !(2..=3).contains(&params.dim) || n < 2 {
            return Err("DPD needs dimension 2 or 3 and at least two particles".into());
        }
        if (0..params.dim).any(|a| params.box_len[a] < 3.0 * params.rc) {
            return Err("box must be at least three cut-off radii wide".into());
        }
        let mut rng = Rng::new(seed);
        let mut pos = vec![[0.0; 3]; n];
        let mut vel = vec![[0.0; 3]; n];
        let sd = (params.kt / params.mass).sqrt();
        for i in 0..n {
            for a in 0..params.dim {
                pos[i][a] = rng.uniform() * params.box_len[a];
                vel[i][a] = rng.normal() * sd;
            }
        }
        let mut s = Self {
            params,
            pos,
            vel,
            force: vec![[0.0; 3]; n],
            rng,
            virial: 0.0,
            t: 0.0,
            wall: None,
            crossed: vec![false; n],
        };
        let mom = s.momentum();
        for v in &mut s.vel {
            for a in 0..params.dim {
                v[a] -= mom[a] / (n as f64 * params.mass);
            }
        }
        s.compute_forces(1.0, true);
        Ok(s)
    }

    /// Replaces the periodic boundary along `wall.axis` by a slab bounded by two walls at `0`
    /// and `L`. Particles outside are folded back and reflected; the wall force acts on later
    /// steps. Pair forces along that axis use the plain (not minimum-image) distance, and
    /// [`Dpd::pressure`] includes the virial of the wall force. Only one axis can be walled.
    ///
    /// # Errors
    /// For an axis outside the dimension, a negative `a_wall`, or a
    /// [`WallKind::NoSlip`] friction outside `[0, 1]`.
    pub fn set_walls(
        &mut self,
        wall: DpdWall,
    ) -> Result<(), String> {
        if wall.axis >= self.params.dim || wall.a_wall < 0.0 {
            return Err("invalid wall".into());
        }
        if let WallKind::NoSlip { friction } = wall.kind
            && !(0.0..=1.0).contains(&friction)
        {
            return Err("wall friction must lie in [0, 1]".into());
        }
        self.wall = Some(wall);
        self.reflect();
        self.compute_forces(1.0, true);
        Ok(())
    }

    /// Mean velocity component `axis` over all particles.
    #[must_use]
    pub fn mean_velocity(
        &self,
        axis: usize,
    ) -> f64 {
        self.vel.iter().map(|v| v[axis]).sum::<f64>() / self.len() as f64
    }

    /// Number of particles.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.pos.len()
    }

    /// Whether there are no particles.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.pos.is_empty()
    }

    /// Number density.
    #[must_use]
    pub fn density(&self) -> f64 {
        let vol: f64 = (0..self.params.dim)
            .map(|a| self.params.box_len[a])
            .product();
        self.len() as f64 / vol
    }

    /// Total momentum.
    #[must_use]
    pub fn momentum(&self) -> [f64; 3] {
        let mut m = [0.0; 3];
        for v in &self.vel {
            for a in 0..3 {
                m[a] += self.params.mass * v[a];
            }
        }
        m
    }

    /// Kinetic temperature `m Σ v² / (d (N - 1))`.
    #[must_use]
    pub fn temperature(&self) -> f64 {
        let ke: f64 = self
            .vel
            .iter()
            .map(|v| v.iter().map(|x| x * x).sum::<f64>())
            .sum();
        self.params.mass * ke / (self.params.dim as f64 * (self.len() - 1) as f64)
    }

    /// Pressure `ρ k_B T_kin + (virial + wall virial) / (d V)` from the last force evaluation.
    ///
    /// The virial sums `Σ a w r` over interacting pairs plus `Σ f_wall d` over the particles
    /// within `r_c` of a wall, `d` being the distance to that wall. The result is the trace of
    /// the pressure tensor over `d`; a slab is not isotropic, so the normal and tangential
    /// pressures differ even though this single number does not say so.
    #[must_use]
    pub fn pressure(&self) -> f64 {
        let vol: f64 = (0..self.params.dim)
            .map(|a| self.params.box_len[a])
            .product();
        self.density() * self.temperature() + self.virial / (self.params.dim as f64 * vol)
    }

    fn wrap(&mut self) {
        let wall_axis = self.wall.map(|w| w.axis);
        for p in &mut self.pos {
            for a in 0..self.params.dim {
                if Some(a) != wall_axis {
                    p[a] = p[a].rem_euclid(self.params.box_len[a]);
                }
            }
        }
        self.reflect();
    }

    /// Applies the reflection rule of the wall to one velocity: the normal component changes
    /// sign and bounce-back also reverses the tangential ones. The tangential friction of a
    /// [`WallKind::NoSlip`] wall is not applied here (see [`Self::friction_velocity`]).
    fn reflect_velocity(
        &self,
        v: &mut [f64; 3],
    ) {
        let Some(w) = self.wall else { return };
        v[w.axis] = -v[w.axis];
        if w.kind == WallKind::BounceBack {
            for a in (0..self.params.dim).filter(|&a| a != w.axis) {
                v[a] = -v[a];
            }
        }
    }

    /// Applies the tangential friction of a [`WallKind::NoSlip`] wall to the final velocity of a
    /// step, once per crossing. The friction is a Langevin collision, `v_t -> (1 - f) v_t + ξ`
    /// with `⟨ξ²⟩ = f (2 - f) k_B T / m`, so a drift decays while a Maxwell distribution at `k_B T`
    /// is left invariant; a purely dissipative rule would cool the fluid.
    fn friction_velocity(
        &mut self,
        v: &mut [f64; 3],
    ) {
        let Some(w) = self.wall else { return };
        let WallKind::NoSlip { friction } = w.kind else { return };
        let keep = 1.0 - friction;
        let sd = (friction * (2.0 - friction) * self.params.kt / self.params.mass).sqrt();
        for a in (0..self.params.dim).filter(|&a| a != w.axis) {
            v[a] = keep * v[a] + sd * self.rng.normal();
        }
    }

    /// Folds particles that left the wall slab back inside and reflects their velocity; records
    /// which particles crossed an odd number of times (the final velocity of the step is
    /// reflected for those as well).
    fn reflect(&mut self) {
        let Some(w) = self.wall else { return };
        let len = self.params.box_len[w.axis];
        self.crossed = vec![false; self.pos.len()];
        for i in 0..self.pos.len() {
            // every crossing folds the position and reflects the velocity once
            while self.pos[i][w.axis] < 0.0 || self.pos[i][w.axis] > len {
                let z = self.pos[i][w.axis];
                self.pos[i][w.axis] = if z < 0.0 { -z } else { 2.0 * len - z };
                let mut v = self.vel[i];
                self.reflect_velocity(&mut v);
                self.vel[i] = v;
                self.crossed[i] = !self.crossed[i];
            }
        }
    }

    fn compute_forces(
        &mut self,
        dt: f64,
        rebuild: bool,
    ) {
        let _ = rebuild;
        let p = self.params;
        let n = self.len();
        let ncell: Vec<i64> = (0..p.dim)
            .map(|a| ((p.box_len[a] / p.rc).floor() as i64).max(1))
            .collect();
        let mut cells: HashMap<[i64; 3], Vec<usize>> = HashMap::new();
        let wall_axis = self.wall.map(|w| w.axis);
        let key = |x: &[f64; 3]| -> [i64; 3] {
            let mut k = [0i64; 3];
            for a in 0..p.dim {
                let c = ((x[a] / p.box_len[a]) * ncell[a] as f64).floor() as i64;
                k[a] = if Some(a) == wall_axis {
                    c.clamp(0, ncell[a] - 1)
                } else {
                    c.rem_euclid(ncell[a])
                };
            }
            k
        };
        for i in 0..n {
            cells.entry(key(&self.pos[i])).or_default().push(i);
        }
        let sigma = (2.0 * p.gamma * p.kt).sqrt();
        let mut f = vec![[0.0; 3]; n];
        let mut virial = 0.0;
        let ox: Vec<i64> = vec![-1, 0, 1];
        let oy: Vec<i64> = if p.dim >= 2 {
            vec![-1, 0, 1]
        } else {
            vec![0]
        };
        let oz: Vec<i64> = if p.dim >= 3 {
            vec![-1, 0, 1]
        } else {
            vec![0]
        };
        for i in 0..n {
            let ki = key(&self.pos[i]);
            let mut seen: Vec<[i64; 3]> = Vec::new();
            for &dx in &ox {
                for &dy in &oy {
                    for &dz in &oz {
                        let off = [dx, dy, dz];
                        let mut kk = [0i64; 3];
                        let mut outside = false;
                        for a in 0..3 {
                            kk[a] = if a < p.dim {
                                if Some(a) == wall_axis {
                                    let c = ki[a] + off[a];
                                    outside |= c < 0 || c >= ncell[a];
                                    c.clamp(0, ncell[a] - 1)
                                } else {
                                    (ki[a] + off[a]).rem_euclid(ncell[a])
                                }
                            } else {
                                0
                            };
                        }
                        if outside || seen.contains(&kk) {
                            continue;
                        }
                        seen.push(kk);
                        let Some(list) = cells.get(&kk) else { continue };
                        for &j in list {
                            if j <= i {
                                continue;
                            }
                            let mut d = [0.0; 3];
                            let mut r2 = 0.0;
                            for a in 0..p.dim {
                                let mut x = self.pos[i][a] - self.pos[j][a];
                                if Some(a) != wall_axis {
                                    x -= p.box_len[a] * (x / p.box_len[a]).round();
                                }
                                d[a] = x;
                                r2 += x * x;
                            }
                            if r2 >= p.rc * p.rc || r2 == 0.0 {
                                continue;
                            }
                            let r = r2.sqrt();
                            let w = 1.0 - r / p.rc;
                            let mut vr = 0.0;
                            for a in 0..p.dim {
                                vr += (self.vel[i][a] - self.vel[j][a]) * d[a] / r;
                            }
                            let xi = self.rng.normal();
                            let mag = p.a * w - p.gamma * w * w * vr + sigma * w * xi / dt.sqrt();
                            virial += p.a * w * r;
                            for a in 0..p.dim {
                                let fa = mag * d[a] / r;
                                f[i][a] += fa;
                                f[j][a] -= fa;
                            }
                        }
                    }
                }
            }
        }
        if let Some(w) = self.wall
            && w.a_wall > 0.0
        {
            let len = p.box_len[w.axis];
            for (i, fi) in f.iter_mut().enumerate() {
                let z = self.pos[i][w.axis];
                if z < p.rc {
                    let fw = w.a_wall * (1.0 - z / p.rc);
                    fi[w.axis] += fw;
                    // the wall force is external, so its virial is the force times the distance
                    // from the wall plane it acts through, not a pair separation
                    virial += fw * z;
                }
                let d = len - z;
                if d < p.rc {
                    let fw = w.a_wall * (1.0 - d / p.rc);
                    fi[w.axis] -= fw;
                    virial += fw * d;
                }
            }
        }
        self.force = f;
        self.virial = virial;
    }

    /// One modified velocity-Verlet step (Groot–Warren, `λ = 0.65`).
    pub fn step(
        &mut self,
        dt: f64,
    ) {
        let lambda = 0.65;
        let m = self.params.mass;
        let dim = self.params.dim;
        let f_old = self.force.clone();
        let v_old = self.vel.clone();
        for i in 0..self.len() {
            for a in 0..dim {
                self.pos[i][a] += dt * self.vel[i][a] + 0.5 * dt * dt * f_old[i][a] / m;
                self.vel[i][a] += lambda * dt * f_old[i][a] / m;
            }
        }
        self.wrap();
        self.compute_forces(dt, true);
        for i in 0..self.len() {
            if self.wall.is_some() && self.crossed[i] {
                // The particle is now on the reflected branch of its path, so the part of the
                // velocity that predates the collision is mirrored (the half kick from `f_old`
                // included) and only the post-collision impulse from `f_new` is added on top.
                // Mirroring `f_new` as well would apply the wall's own impulse in the wrong
                // direction and pump energy in on every crossing.
                let mut v = [0.0; 3];
                for a in 0..dim {
                    v[a] = v_old[i][a] + 0.5 * dt * f_old[i][a] / m;
                }
                self.reflect_velocity(&mut v);
                for a in 0..dim {
                    v[a] += 0.5 * dt * self.force[i][a] / m;
                }
                self.friction_velocity(&mut v);
                self.vel[i] = v;
            } else {
                for a in 0..dim {
                    self.vel[i][a] = v_old[i][a] + 0.5 * dt * (f_old[i][a] + self.force[i][a]) / m;
                }
            }
        }
        self.t += dt;
    }

    /// Runs `steps` steps.
    pub fn run(
        &mut self,
        dt: f64,
        steps: usize,
    ) {
        for _ in 0..steps {
            self.step(dt);
        }
    }
}
