//! # Dissipative particle dynamics (DPD) in 2D and 3D
//!
//! Groot–Warren DPD fluid in a periodic box: soft conservative repulsion
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
//! | Boundary | periodic box (minimum image), cell-list neighbour search |
//! | Diagnostics | kinetic temperature, momentum, virial pressure |
//! | Equation of state | `p = ρ k_B T + α a ρ²`; measured `α ≈ 0.092` at `a = 25`, `ρ = 3` (3D, `r_c = 1`, independent of `Δt`), against `0.101` quoted by Groot and Warren |
//!
//! Limitations: single species, no bonds (polymers) and no walls.
#![allow(
    clippy::manual_memcpy,
    clippy::too_long_first_doc_paragraph,
    clippy::cast_lossless,
    clippy::too_many_lines,
    clippy::too_many_arguments,
    clippy::needless_range_loop,
    clippy::missing_const_for_fn,
    clippy::many_single_char_names,
    clippy::similar_names,
    clippy::cognitive_complexity,
    clippy::option_if_let_else
)]

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
}

impl Dpd {
    /// Random uniform positions and Maxwell velocities with zero total momentum.
    ///
    /// # Errors
    /// For a bad dimension, a box smaller than three cut-offs, or fewer than two particles.
    pub fn new(n: usize, params: DpdParams, seed: u64) -> Result<Self, String> {
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
        let mut s = Self { params, pos, vel, force: vec![[0.0; 3]; n], rng, virial: 0.0, t: 0.0 };
        let mom = s.momentum();
        for v in &mut s.vel {
            for a in 0..params.dim {
                v[a] -= mom[a] / (n as f64 * params.mass);
            }
        }
        s.compute_forces(1.0, true);
        Ok(s)
    }

    /// Number of particles.
    #[must_use]
    pub fn len(&self) -> usize {
        self.pos.len()
    }

    /// Whether there are no particles.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.pos.is_empty()
    }

    /// Number density.
    #[must_use]
    pub fn density(&self) -> f64 {
        let vol: f64 = (0..self.params.dim).map(|a| self.params.box_len[a]).product();
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
        let ke: f64 = self.vel.iter().map(|v| v.iter().map(|x| x * x).sum::<f64>()).sum();
        self.params.mass * ke / (self.params.dim as f64 * (self.len() - 1) as f64)
    }

    /// Pressure `ρ k_B T_kin + virial / (d V)` from the last force evaluation.
    #[must_use]
    pub fn pressure(&self) -> f64 {
        let vol: f64 = (0..self.params.dim).map(|a| self.params.box_len[a]).product();
        self.density() * self.temperature() + self.virial / (self.params.dim as f64 * vol)
    }

    fn wrap(&mut self) {
        for p in &mut self.pos {
            for a in 0..self.params.dim {
                p[a] = p[a].rem_euclid(self.params.box_len[a]);
            }
        }
    }

    fn compute_forces(&mut self, dt: f64, rebuild: bool) {
        let _ = rebuild;
        let p = self.params;
        let n = self.len();
        let ncell: Vec<i64> = (0..p.dim).map(|a| ((p.box_len[a] / p.rc).floor() as i64).max(1)).collect();
        let mut cells: HashMap<[i64; 3], Vec<usize>> = HashMap::new();
        let key = |x: &[f64; 3]| -> [i64; 3] {
            let mut k = [0i64; 3];
            for a in 0..p.dim {
                k[a] = (((x[a] / p.box_len[a]) * ncell[a] as f64).floor() as i64).rem_euclid(ncell[a]);
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
        let oy: Vec<i64> = if p.dim >= 2 { vec![-1, 0, 1] } else { vec![0] };
        let oz: Vec<i64> = if p.dim >= 3 { vec![-1, 0, 1] } else { vec![0] };
        for i in 0..n {
            let ki = key(&self.pos[i]);
            let mut seen: Vec<[i64; 3]> = Vec::new();
            for &dx in &ox {
                for &dy in &oy {
                    for &dz in &oz {
                        let off = [dx, dy, dz];
                        let mut kk = [0i64; 3];
                        for a in 0..3 {
                            kk[a] = if a < p.dim { (ki[a] + off[a]).rem_euclid(ncell[a]) } else { 0 };
                        }
                        if seen.contains(&kk) {
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
                                x -= p.box_len[a] * (x / p.box_len[a]).round();
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
        self.force = f;
        self.virial = virial;
    }

    /// One modified velocity-Verlet step (Groot–Warren, `λ = 0.65`).
    pub fn step(&mut self, dt: f64) {
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
            for a in 0..dim {
                self.vel[i][a] = v_old[i][a] + 0.5 * dt * (f_old[i][a] + self.force[i][a]) / m;
            }
        }
        self.t += dt;
    }

    /// Runs `steps` steps.
    pub fn run(&mut self, dt: f64, steps: usize) {
        for _ in 0..steps {
            self.step(dt);
        }
    }
}
