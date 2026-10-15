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



use crate::sim::symplectic::{ttl_step, TtlState};

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
            let next: Vec<f64> = (0..prev.len()).map(|i| prev[i] + (prev[i] - up[i]) / ratio).collect();
            row.push(next);
        }
        t.push(row);
    }
    let last = &t[kmax - 1];
    let best = last[kmax - 1].clone();
    let err = if kmax >= 2 {
        let other = &last[kmax - 2];
        let scale = best.iter().fold(0.0_f64, |m, x| m.max(x.abs())).max(1e-300);
        best.iter().zip(other).fold(0.0_f64, |m, (a, b)| m.max((a - b).abs())) / scale
    } else {
        0.0
    };
    (best, err)
}

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

fn norm3(a: V3) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
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
    /// Fails for fewer than two bodies, length mismatch or non-positive masses.
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
        if mass.iter().any(|&m| m <= 0.0) {
            return Err("masses must be positive".into());
        }
        let mt: f64 = mass.iter().sum();
        let mut com_pos = [0.0; 3];
        let mut com_vel = [0.0; 3];
        for i in 0..n {
            for a in 0..3 {
                com_pos[a] += mass[i] * pos[i][a] / mt;
                com_vel[a] += mass[i] * vel[i][a] / mt;
            }
        }
        let rel: Vec<V3> = pos.iter().map(|p| sub3(*p, com_pos)).collect();
        let order = greedy_chain(&rel);
        let xc: Vec<V3> = (0..n - 1).map(|k| sub3(pos[order[k + 1]], pos[order[k]])).collect();
        let vc: Vec<V3> = (0..n - 1).map(|k| sub3(vel[order[k + 1]], vel[order[k]])).collect();
        let mut s = Self { g, mass: mass.to_vec(), order, xc, vc, com_pos, com_vel, t: 0.0, b: 0.0 };
        s.b = -s.energy();
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
        (0..mc.len()).map(|k| 0.5 * mc[k] * w[k].iter().map(|x| x * x).sum::<f64>()).sum()
    }

    /// Accelerations per chain position and the positive potential `U = -V`.
    fn chain_accelerations(
        &self,
        xc: &[V3],
    ) -> (Vec<V3>, f64) {
        let n = self.order.len();
        let mc = self.masses_in_chain();
        let mut acc = vec![[0.0; 3]; n];
        let mut u = 0.0;
        for a in 0..n {
            let mut r = [0.0; 3];
            for b in a + 1..n {
                r = add3(r, xc[b - 1]); // chain sum of vectors between a and b
                let d = norm3(r);
                let d3 = d * d * d;
                u += self.g * mc[a] * mc[b] / d;
                for k in 0..3 {
                    acc[a][k] += self.g * mc[b] * r[k] / d3;
                    acc[b][k] -= self.g * mc[a] * r[k] / d3;
                }
            }
        }
        (acc, u)
    }

    /// Potential energy (negative).
    #[must_use]
    pub fn potential(&self) -> f64 {
        -self.chain_accelerations(&self.xc).1
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
                pos[self.order[k]][a] = y[k][a] - mean[a] + self.com_pos[a] + self.com_vel[a] * self.t;
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
            let (lo, hi, sign) = if from <= to { (from, to, 1.0) } else { (to, from, -1.0) };
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
        self.half_drift(0.5 * h);
        let (acc, u) = self.chain_accelerations(&self.xc);
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
        let (y, err) = gbs(kmax, |m| {
            let mut s = start.clone();
            for _ in 0..m {
                s.leapfrog(ds / m as f64);
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
            if dh < best.2 {
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
