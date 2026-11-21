//! Chain and algorithmic regularisation for few-body gravity.
//!
//! # Chain regularisation and algorithmic regularisation for few-body gravity
//!
//! Gravitational `N`-body dynamics (3-D, `G` explicit) integrated in the fictitious time
//! `s` with `dt = ds / (T + B)` (`B = -E_0`, the logarithmic-Hamiltonian time
//! transformation of Mikkola and Tanikawa), which makes the step size follow the
//! shortest dynamical time of the system. Combined with
//!
//! * the **Mikkola-Aarseth chain**: positions and velocities are stored as `N - 1`
//!   vectors between neighbours of a chain that links every body to its nearest
//!   companions (re-selected as the configuration evolves, exactly, by summing chain
//!   vectors), so tight binaries keep full relative precision even far from the origin;
//! * **Gragg-Bulirsch-Stoer extrapolation** of the (symmetric, reversible) leapfrog in
//!   `s`, giving an arbitrary-order, error-controlled scheme, with adaptive `ds`.
//!
//! The Cartesian counterpart [`ttl_gbs_step`] applies the same extrapolation to the
//! generic [`crate::sim::symplectic::TtlState`] leapfrog, and
//! [`crate::sim::symplectic::ttl_step_composed`] provides composition schemes of order 4, 6, 8
//! for the time-transformed leapfrog.

use crate::kernels::numeric_common::scale_product;
use crate::sim::symplectic::TtlState;
use crate::sim::symplectic::ttl_step;

/// Gragg-Bulirsch-Stoer extrapolation of a symmetric one-step method.
///
/// `run(n)` returns the solution after one macro step computed with `n` equal substeps.
/// The error expansion of a symmetric method is even in the step, so the Richardson
/// table in `h^2` is used with the sequence `n = 2, 4, ..., 2 kmax`. Returns the
/// extrapolated state and a relative error estimate (difference of the last two diagonal entries).
fn gbs<F: FnMut(usize) -> Vec<f64>>(
    kmax: usize,
    mut run: F,
) -> (Vec<f64>, f64) {
    let kmax = kmax.max(1);
    let seq: Vec<usize> = (1..=kmax).map(|k| 2 * k).collect();
    let mut t: Vec<Vec<Vec<f64>>> = Vec::new();
    for k in 0..kmax {
        let mut row = vec![run(seq[k])];
        for j in 1..=k {
            let ratio = (seq[k] as f64 / seq[k - j] as f64).powi(2) - 1.0;
            let prev = &row[j - 1];
            let up = &t[k - 1][j - 1];
            let next: Vec<f64> = (0..prev.len())
                .map(|i| prev[i] + (prev[i] - up[i]) / ratio)
                .collect();
            row.push(next);
        }
        t.push(row);
    }
    let last = &t[kmax - 1];
    let best = last[kmax - 1].clone();
    let err = if kmax >= 2 {
        let other = &last[kmax - 2];
        let scale = best.iter().fold(0.0_f64, |m, x| m.max(x.abs())).max(1e-300);
        best.iter()
            .zip(other)
            .fold(0.0_f64, |m, (a, b)| m.max((a - b).abs()))
            / scale
    } else {
        0.0
    };
    (best, err)
}

/// One Gragg-Bulirsch-Stoer extrapolated step of the time-transformed leapfrog.
///
/// One Gragg-Bulirsch-Stoer extrapolated step of the time-transformed leapfrog
/// ([`ttl_step`]) of size `ds` in the fictitious time, using `kmax` columns.
/// Returns the relative error estimate. The scheme has order `2 kmax` in `ds`.
pub fn ttl_gbs_step<F>(
    s: &mut TtlState,
    mass: &[f64],
    b: f64,
    ds: f64,
    kmax: usize,
    potential_grad: &mut F,
) -> f64
where
    F: FnMut(&[f64], &mut [f64]) -> f64,
{
    let n = s.q.len();
    let start = s.clone();
    let (y, err) = gbs(kmax, |m| {
        let mut st = start.clone();
        for _ in 0..m {
            ttl_step(&mut st, mass, b, ds / m as f64, potential_grad);
        }
        let mut y = st.q.clone();
        y.extend_from_slice(&st.p);
        y.push(st.t);
        y
    });
    s.q.copy_from_slice(&y[..n]);
    s.p.copy_from_slice(&y[n..2 * n]);
    s.t = y[2 * n];
    err
}

type V3 = [f64; 3];

fn add3(
    a: V3,
    b: V3,
) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn sub3(
    a: V3,
    b: V3,
) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn squared_norm3(a: V3) -> f64 {
    a[0] * a[0] + a[1] * a[1] + a[2] * a[2]
}

fn norm3(a: V3) -> f64 {
    let squared = squared_norm3(a);
    if squared.is_normal() {
        squared.sqrt()
    } else {
        a[0].hypot(a[1]).hypot(a[2])
    }
}

fn scaled_radius(
    r: V3,
    g: f64,
    mass_a: f64,
    mass_b: f64,
) -> Option<(f64, f64)> {
    if !g.is_finite()
        || !mass_a.is_finite()
        || !mass_b.is_finite()
        || mass_a <= 0.0
        || mass_b <= 0.0
        || r.iter().any(|x| !x.is_finite())
    {
        return None;
    }
    let scale = r.iter().map(|x| x.abs()).fold(0.0, f64::max);
    (scale > 0.0).then(|| (scale, squared_norm3(r.map(|x| x / scale)).sqrt()))
}

fn scaled_interaction(
    value: f64,
    numerator: &[f64],
    denominator: &[f64],
) -> f64 {
    // Pair geometry and reference factors have already been validated. The
    // remaining failure is an interaction outside the f64 output range.
    scale_product(value, numerator, denominator).unwrap_or_else(|_| f64::INFINITY.copysign(value))
}

/// Gravitational `N`-body system in the Mikkola-Aarseth chain representation.
///
/// The centre of mass is separated (`com_pos`, `com_vel`); the relative motion is stored in
/// the chain vectors `xc[k] = x_{order[k+1]} - x_{order[k]}` and the velocity chain `vc`.
#[derive(Clone, Debug)]
pub struct ChainSystem {
    /// Gravitational constant.
    pub g: f64,
    /// Masses by body index.
    pub mass: Vec<f64>,
    /// Chain order: `order[k]` is the body at chain position `k`.
    pub order: Vec<usize>,
    /// Chain position vectors.
    pub xc: Vec<V3>,
    /// Chain velocity vectors.
    pub vc: Vec<V3>,
    /// Centre-of-mass position at `t = 0`.
    pub com_pos: V3,
    /// Centre-of-mass velocity.
    pub com_vel: V3,
    /// Physical time.
    pub t: f64,
    /// Binding parameter `B = -E_0` of the time transformation `dt = ds / (T + B)`.
    pub b: f64,
}

impl ChainSystem {
    /// Builds the chain system from Cartesian data.
    ///
    /// # Errors
    /// Fails for fewer than two bodies, length mismatch, non-positive
    /// or non-finite masses, non-finite inputs, or an unrepresentable initial
    /// chain or energy. The unsoftened interaction is singular at coincidence.
    pub fn new(
        g: f64,
        mass: &[f64],
        pos: &[V3],
        vel: &[V3],
    ) -> Result<Self, String> {
        let n = mass.len();
        if n < 2 || pos.len() != n || vel.len() != n {
            return Err("need at least two bodies with matching position/velocity data".into());
        }
        if mass.iter().any(|&m| !m.is_finite() || m <= 0.0) {
            return Err("masses must be finite and positive".into());
        }
        if !g.is_finite() || pos.iter().chain(vel).flatten().any(|v| !v.is_finite()) {
            return Err("chain gravity, positions and velocities must be finite".into());
        }
        let mt: f64 = mass.iter().sum();
        if !mt.is_finite() {
            return Err("total chain mass is not representable".into());
        }
        let mut com_pos = [0.0; 3];
        let mut com_vel = [0.0; 3];
        for i in 0..n {
            for a in 0..3 {
                com_pos[a] += mass[i] * pos[i][a] / mt;
                com_vel[a] += mass[i] * vel[i][a] / mt;
            }
        }
        if com_pos.iter().chain(&com_vel).any(|v| !v.is_finite()) {
            return Err("chain centre of mass is not representable".into());
        }
        let rel: Vec<V3> = pos.iter().map(|p| sub3(*p, com_pos)).collect();
        let order = greedy_chain(&rel);
        let xc: Vec<V3> = (0..n - 1)
            .map(|k| sub3(pos[order[k + 1]], pos[order[k]]))
            .collect();
        let vc: Vec<V3> = (0..n - 1)
            .map(|k| sub3(vel[order[k + 1]], vel[order[k]]))
            .collect();
        if xc.iter().chain(&vc).flatten().any(|v| !v.is_finite()) {
            return Err("relative chain state is not representable".into());
        }
        let mut s = Self {
            g,
            mass: mass.to_vec(),
            order,
            xc,
            vc,
            com_pos,
            com_vel,
            t: 0.0,
            b: 0.0,
        };
        s.b = -s.energy();
        if !s.b.is_finite() {
            return Err(
                "initial chain energy is not finite (coincident bodies are singular)".into(),
            );
        }
        Ok(s)
    }

    fn masses_in_chain(&self) -> Vec<f64> {
        self.order.iter().map(|&i| self.mass[i]).collect()
    }

    /// Velocities per chain position (centre-of-mass frame).
    fn chain_velocities(
        &self,
        vc: &[V3],
    ) -> Vec<V3> {
        let n = self.order.len();
        let mc = self.masses_in_chain();
        let mt: f64 = mc.iter().sum();
        let mut w = vec![[0.0; 3]; n];
        for k in 1..n {
            w[k] = add3(w[k - 1], vc[k - 1]);
        }
        let mut mean = [0.0; 3];
        for k in 0..n {
            for a in 0..3 {
                mean[a] += mc[k] * w[k][a] / mt;
            }
        }
        for wk in &mut w {
            *wk = sub3(*wk, mean);
        }
        w
    }

    /// Kinetic energy.
    #[must_use]
    pub fn kinetic(&self) -> f64 {
        let mc = self.masses_in_chain();
        let w = self.chain_velocities(&self.vc);
        (0..mc.len())
            .map(|k| 0.5 * mc[k] * w[k].iter().map(|x| x * x).sum::<f64>())
            .sum()
    }

    /// Accelerations per body in the original input order.
    ///
    /// The returned vectors are Cartesian accelerations. In states mutated
    /// after construction, coincidences retain the gravitational singularity.
    #[must_use]
    pub fn accelerations(&self) -> Vec<V3> {
        let mut acceleration = vec![[0.0; 3]; self.order.len()];
        self.pair_interactions::<true, true>(&mut acceleration);
        acceleration
    }

    /// Traverses pairs in chain order and optionally writes their accelerations.
    /// The accumulated potential is the positive binding contribution `U = -V`.
    fn pair_interactions<const FORCE: bool, const BODY_ORDER: bool>(
        &self,
        acceleration: &mut [V3],
    ) -> f64 {
        if FORCE {
            acceleration.fill([0.0; 3]);
        }
        let ordinary_mass = (self.g == 0.0 || (1e-70..=1e70).contains(&self.g.abs()))
            && self.mass.iter().all(|m| (1e-70..=1e70).contains(m));
        // Nonzero links in this range have dyadic spacing above 1e-96.
        // Summing links cannot create finer nonzero spacing. At the other
        // end an increment bounded by 1e75 stops changing a sum long
        // before its cube could overflow.
        let ordinary_links = self.xc.iter().flatten().all(|x| {
            let component = x.abs();
            component == 0.0 || (1e-80..=1e75).contains(&component)
        });
        if ordinary_mass && ordinary_links {
            self.sum_pairs::<FORCE, BODY_ORDER, false>(acceleration)
        } else {
            self.sum_pairs::<FORCE, BODY_ORDER, true>(acceleration)
        }
    }

    fn sum_pairs<const FORCE: bool, const BODY_ORDER: bool, const SCALED: bool>(
        &self,
        acceleration: &mut [V3],
    ) -> f64 {
        let g = self.g;
        let mut binding = 0.0;
        for (a, &body_a) in self.order.iter().enumerate() {
            let mass_a = self.mass[body_a];
            let force_a = if BODY_ORDER { body_a } else { a };
            let gm_a = g * mass_a;
            let mut r = [0.0; 3];
            for (offset, &body_b) in self.order[a + 1..].iter().enumerate() {
                let b = a + 1 + offset;
                r = add3(r, self.xc[b - 1]);
                let mass_b = self.mass[body_b];
                let force_b = if BODY_ORDER { body_b } else { b };
                let radius = if SCALED {
                    scaled_radius(r, g, mass_a, mass_b)
                } else {
                    None
                };
                if let Some((scale, unit)) = radius {
                    binding += scaled_interaction(g, &[mass_a, mass_b], &[scale, unit]);
                    if FORCE {
                        let factors = [scale, scale, scale, unit, unit, unit];
                        for (i, &ri) in r.iter().enumerate() {
                            if ri != 0.0 {
                                let signed_g = if ri.is_sign_negative() {
                                    -g
                                } else {
                                    g
                                };
                                acceleration[force_a][i] +=
                                    scaled_interaction(signed_g, &[mass_b, ri.abs()], &factors);
                                acceleration[force_b][i] -=
                                    scaled_interaction(signed_g, &[mass_a, ri.abs()], &factors);
                            }
                        }
                    }
                } else {
                    // Ordinary arithmetic also preserves the singularity
                    // in a manually supplied coincident state.
                    let distance = squared_norm3(r).sqrt();
                    binding += gm_a * mass_b / distance;
                    if FORCE {
                        let cube = distance * distance * distance;
                        let gm_b = g * mass_b;
                        for (i, &ri) in r.iter().enumerate() {
                            acceleration[force_a][i] += gm_b * ri / cube;
                            acceleration[force_b][i] -= gm_a * ri / cube;
                        }
                    }
                }
            }
        }
        binding
    }

    /// Potential energy (negative for attractive gravity).
    #[must_use]
    pub fn potential(&self) -> f64 {
        -self.pair_interactions::<false, false>(&mut [])
    }

    /// Total energy.
    #[must_use]
    pub fn energy(&self) -> f64 {
        self.kinetic() + self.potential()
    }

    /// Cartesian positions and velocities (original body order, absolute frame) at the current time.
    #[must_use]
    pub fn cartesian(&self) -> (Vec<V3>, Vec<V3>) {
        let n = self.order.len();
        let mc = self.masses_in_chain();
        let mt: f64 = mc.iter().sum();
        let mut y = vec![[0.0; 3]; n];
        for k in 1..n {
            y[k] = add3(y[k - 1], self.xc[k - 1]);
        }
        let mut mean = [0.0; 3];
        for k in 0..n {
            for a in 0..3 {
                mean[a] += mc[k] * y[k][a] / mt;
            }
        }
        let w = self.chain_velocities(&self.vc);
        let mut pos = vec![[0.0; 3]; n];
        let mut vel = vec![[0.0; 3]; n];
        for k in 0..n {
            for a in 0..3 {
                pos[self.order[k]][a] =
                    y[k][a] - mean[a] + self.com_pos[a] + self.com_vel[a] * self.t;
                vel[self.order[k]][a] = w[k][a] + self.com_vel[a];
            }
        }
        (pos, vel)
    }

    /// Angular momentum about the centre of mass.
    #[must_use]
    pub fn angular_momentum(&self) -> V3 {
        let n = self.order.len();
        let mc = self.masses_in_chain();
        let mt: f64 = mc.iter().sum();
        let mut y = vec![[0.0; 3]; n];
        for k in 1..n {
            y[k] = add3(y[k - 1], self.xc[k - 1]);
        }
        let mut mean = [0.0; 3];
        for k in 0..n {
            for a in 0..3 {
                mean[a] += mc[k] * y[k][a] / mt;
            }
        }
        let w = self.chain_velocities(&self.vc);
        let mut l = [0.0; 3];
        for k in 0..n {
            let r = sub3(y[k], mean);
            l[0] += mc[k] * (r[1] * w[k][2] - r[2] * w[k][1]);
            l[1] += mc[k] * (r[2] * w[k][0] - r[0] * w[k][2]);
            l[2] += mc[k] * (r[0] * w[k][1] - r[1] * w[k][0]);
        }
        l
    }

    /// Re-selects the chain from the current configuration; the new chain vectors are
    /// exact sums of the old ones, so no precision is lost.
    pub fn rechain(&mut self) {
        let n = self.order.len();
        let mut y = vec![[0.0; 3]; n];
        for k in 1..n {
            y[k] = add3(y[k - 1], self.xc[k - 1]);
        }
        let mut by_body = vec![[0.0; 3]; n];
        let mut pos_of = vec![0usize; n];
        for k in 0..n {
            by_body[self.order[k]] = y[k];
            pos_of[self.order[k]] = k;
        }
        let new_order = greedy_chain(&by_body);
        if new_order == self.order {
            return;
        }
        let range_sum = |v: &[V3], from: usize, to: usize| -> V3 {
            let (lo, hi, sign) = if from <= to {
                (from, to, 1.0)
            } else {
                (to, from, -1.0)
            };
            let mut s = [0.0; 3];
            for item in &v[lo..hi] {
                s = add3(s, *item);
            }
            [sign * s[0], sign * s[1], sign * s[2]]
        };
        let mut nx = Vec::with_capacity(n - 1);
        let mut nv = Vec::with_capacity(n - 1);
        for k in 0..n - 1 {
            let (a, b) = (pos_of[new_order[k]], pos_of[new_order[k + 1]]);
            nx.push(range_sum(&self.xc, a, b));
            nv.push(range_sum(&self.vc, a, b));
        }
        self.order = new_order;
        self.xc = nx;
        self.vc = nv;
    }

    fn half_drift(
        &mut self,
        h: f64,
    ) {
        let dt = h / (self.kinetic() + self.b);
        for k in 0..self.xc.len() {
            for a in 0..3 {
                self.xc[k][a] += dt * self.vc[k][a];
            }
        }
        self.t += dt;
    }

    /// One second-order time-transformed leapfrog step of size `h` in `s` (chain coordinates).
    pub fn leapfrog(
        &mut self,
        h: f64,
    ) {
        let mut acceleration = vec![[0.0; 3]; self.order.len()];
        self.leapfrog_with_acceleration(h, &mut acceleration);
    }

    fn leapfrog_with_acceleration(
        &mut self,
        h: f64,
        acc: &mut [V3],
    ) {
        self.half_drift(0.5 * h);
        let u = self.pair_interactions::<true, false>(acc);
        let dt = h / u;
        for k in 0..self.vc.len() {
            for a in 0..3 {
                self.vc[k][a] += dt * (acc[k + 1][a] - acc[k][a]);
            }
        }
        self.half_drift(0.5 * h);
    }

    fn pack(&self) -> Vec<f64> {
        let mut y = Vec::with_capacity(6 * self.xc.len() + 1);
        for x in &self.xc {
            y.extend_from_slice(x);
        }
        for v in &self.vc {
            y.extend_from_slice(v);
        }
        y.push(self.t);
        y
    }

    fn unpack(
        &mut self,
        y: &[f64],
    ) {
        let m = self.xc.len();
        for k in 0..m {
            for a in 0..3 {
                self.xc[k][a] = y[3 * k + a];
                self.vc[k][a] = y[3 * m + 3 * k + a];
            }
        }
        self.t = y[6 * m];
    }

    /// One Gragg-Bulirsch-Stoer extrapolated step of size `ds` with `kmax` columns (order
    /// `2 kmax`) at fixed chain. Returns the relative error estimate.
    pub fn gbs_step(
        &mut self,
        ds: f64,
        kmax: usize,
    ) -> f64 {
        let start = self.clone();
        let mut acceleration = vec![[0.0; 3]; self.order.len()];
        let (y, err) = gbs(kmax, |m| {
            let mut s = start.clone();
            for _ in 0..m {
                s.leapfrog_with_acceleration(ds / m as f64, &mut acceleration);
            }
            s.pack()
        });
        self.unpack(&y);
        err
    }

    /// Advances to the physical time `t_end` with adaptive fictitious steps (initial `ds0`)
    /// and local error tolerance `tol`; the chain is re-selected before every step and the
    /// last step is solved for so that the final time is hit to round-off. Returns the
    /// number of accepted steps.
    ///
    /// # Errors
    /// Fails if the step size collapses or the final-time solve does not converge.
    pub fn evolve_until(
        &mut self,
        t_end: f64,
        ds0: f64,
        kmax: usize,
        tol: f64,
    ) -> Result<usize, String> {
        let mut ds = ds0;
        let mut steps = 0usize;
        let mut guard = 0usize;
        while self.t < t_end {
            guard += 1;
            if guard > 5_000_000 {
                return Err("too many steps".into());
            }
            let backup = self.clone();
            self.rechain();
            let err = self.gbs_step(ds, kmax);
            if err.is_nan() || err > tol || !self.t.is_finite() {
                *self = backup;
                ds *= 0.5;
                if ds < 1e-14 {
                    return Err("fictitious step size collapsed".into());
                }
                continue;
            }
            if self.t > t_end {
                *self = backup;
                self.rechain();
                self.land(t_end, ds, kmax)?;
                steps += 1;
                break;
            }
            steps += 1;
            if err < 0.05 * tol {
                ds *= 1.5;
            }
        }
        Ok(steps)
    }

    fn land(
        &mut self,
        t_end: f64,
        ds_max: f64,
        kmax: usize,
    ) -> Result<(), String> {
        let base = self.clone();
        let mut d = ((t_end - base.t) * (base.kinetic() + base.b)).min(ds_max);
        for _ in 0..60 {
            let mut s = base.clone();
            s.gbs_step(d, kmax);
            let miss = t_end - s.t;
            if miss.abs() <= 1e-13 * t_end.abs().max(1.0) {
                *self = s;
                return Ok(());
            }
            d += miss * (s.kinetic() + s.b);
            if d <= 0.0 {
                d = 1e-12;
            }
        }
        Err("final-time solve did not converge".into())
    }
}

/// Greedy chain: start from the closest pair and repeatedly attach the nearest remaining
/// body to either end.
fn greedy_chain(pos: &[V3]) -> Vec<usize> {
    let n = pos.len();
    let dist = |a: usize, b: usize| norm3(sub3(pos[a], pos[b]));
    let (mut bi, mut bj, mut bd) = (0, 1, f64::INFINITY);
    for i in 0..n {
        for j in i + 1..n {
            let d = dist(i, j);
            if d < bd {
                (bi, bj, bd) = (i, j, d);
            }
        }
    }
    let mut chain = std::collections::VecDeque::from([bi, bj]);
    let mut used = vec![false; n];
    used[bi] = true;
    used[bj] = true;
    while chain.len() < n {
        let (head, tail) = (chain[0], chain[chain.len() - 1]);
        let mut best = (usize::MAX, false, f64::INFINITY);
        for c in 0..n {
            if used[c] {
                continue;
            }
            let (dh, dt) = (dist(c, head), dist(c, tail));
            // Even a separation outside f64 range has a deterministic
            // candidate. This avoids indexing with the absent-body sentinel.
            if best.0 == usize::MAX || dh < best.2 {
                best = (c, true, dh);
            }
            if dt < best.2 {
                best = (c, false, dt);
            }
        }
        used[best.0] = true;
        if best.1 {
            chain.push_front(best.0);
        } else {
            chain.push_back(best.0);
        }
    }
    chain.into_iter().collect()
}
