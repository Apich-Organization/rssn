//! # Stochastic series expansion with directed-loop updates for spin-1/2 XXZ models
//!
//! Quantum Monte Carlo for `H = sum_<ij> [J_z S^z_i S^z_j + (J_xy/2)(S^+_i S^-_j + h.c.)]` with
//! antiferromagnetic transverse coupling `J_xy > 0` on any **bipartite** graph (chains, square
//! lattices, ladders; the sign of the off-diagonal terms is removed by the sublattice rotation).
//! The partition function is expanded in a power series of `H = N_b C - sum_b (H_{1b} + H_{2b})`
//! with `H_{1b} = C - J_z S^z_i S^z_j` and `H_{2b}` the (rotated) spin flip, and sampled by
//!
//! * a diagonal update (insertion / removal of diagonal bond operators with the usual
//!   `N_b beta W / (M - n)` acceptance ratios) and
//! * **directed-loop** (worm) updates: the loop head flips the spin of the leg it enters and
//!   leaves the vertex through a leg drawn by the heat-bath rule `P(out) ~ W(vertex after the
//!   flips)`, which satisfies the directed-loop equations (bounces are allowed). Loops are
//!   repeated until about `2 n` vertex legs have been visited.
//!
//! The operator-string cutoff `M` grows automatically during thermalisation.
//!
//! | Quantity | Estimator |
//! |---|---|
//! | Energy | `E = N_b C - <n>/beta` |
//! | Specific heat | `C_v = <n^2> - <n>^2 - <n>` |
//! | Magnetisation, staggered structure | from the propagated spin state at imaginary time 0 |
//!
//! Limitations: spin-1/2 only, antiferromagnetic `J_xy`, bipartite graphs, `J_z` of either sign;
//! no winding-number (stiffness) estimators and no measurement at general imaginary times.

use crate::kernels::random::Rng;

/// Graph helper: nearest-neighbour bonds of a chain with `l` sites.
#[must_use]
pub fn chain_bonds(
    l: usize,
    periodic: bool,
) -> Vec<(usize, usize)> {
    let mut b: Vec<(usize, usize)> = (0..l.saturating_sub(1)).map(|i| (i, i + 1)).collect();
    if periodic && l > 2 {
        b.push((l - 1, 0));
    }
    b
}

/// Graph helper: nearest-neighbour bonds of an `lx x ly` square lattice (site `x + lx y`).
#[must_use]
pub fn square_bonds(
    lx: usize,
    ly: usize,
    periodic: bool,
) -> Vec<(usize, usize)> {
    let mut b = Vec::new();
    for y in 0..ly {
        for x in 0..lx {
            let i = x + lx * y;
            if x + 1 < lx {
                b.push((i, i + 1));
            } else if periodic && lx > 2 {
                b.push((i, lx * y));
            }
            if y + 1 < ly {
                b.push((i, i + lx));
            } else if periodic && ly > 2 {
                b.push((i, x));
            }
        }
    }
    b
}

const NONE: usize = usize::MAX;

/// SSE simulation of the spin-1/2 XXZ model.
#[derive(Clone, Debug)]
pub struct SseXxz {
    n_sites: usize,
    bonds: Vec<(usize, usize)>,
    /// Ising anisotropy coupling `J_z` (the transverse coupling is `J_xy = 1`).
    pub jz: f64,
    /// Inverse temperature.
    pub beta: f64,
    shift: f64,
    spins: Vec<i8>,
    ops: Vec<u32>,
    n_ops: usize,
    rng: Rng,
    link: Vec<usize>,
    legspin: Vec<i8>,
    first: Vec<usize>,
    last: Vec<usize>,
    colour: Vec<i8>,
}

/// Measurements of one sweep.
#[derive(Clone, Copy, Debug)]
pub struct SseSample {
    /// Number of operators `n`.
    pub n_ops: usize,
    /// Total `S^z` of the propagated state at imaginary time 0.
    pub magnetisation: f64,
    /// Staggered magnetisation (sublattice sign from the bond graph colouring).
    pub staggered: f64,
}

impl SseXxz {
    /// Creates a simulation on the given bonds.
    ///
    /// # Errors
    /// Fails for an empty bond list, out-of-range sites, a non-positive `beta`, or a graph that is
    /// not bipartite.
    pub fn new(
        n_sites: usize,
        bonds: Vec<(usize, usize)>,
        jz: f64,
        beta: f64,
        seed: u64,
    ) -> Result<Self, String> {
        if bonds.is_empty()
            || beta <= 0.0
            || bonds
                .iter()
                .any(|&(i, j)| i >= n_sites || j >= n_sites || i == j)
        {
            return Err("need a non-empty valid bond list and beta > 0".into());
        }
        // bipartite check
        let mut colour = vec![-1i8; n_sites];
        for start in 0..n_sites {
            if colour[start] >= 0 {
                continue;
            }
            colour[start] = 0;
            let mut stack = vec![start];
            while let Some(u) = stack.pop() {
                for &(a, b) in &bonds {
                    let other = if a == u {
                        b
                    } else if b == u {
                        a
                    } else {
                        continue;
                    };
                    if colour[other] < 0 {
                        colour[other] = 1 - colour[u];
                        stack.push(other);
                    } else if colour[other] == colour[u] {
                        return Err("the lattice must be bipartite".into());
                    }
                }
            }
        }
        let mut rng = Rng::new(seed);
        let spins = (0..n_sites)
            .map(|_| if rng.uniform() < 0.5 { 1 } else { -1 })
            .collect();
        let shift = 0.25 * jz.abs() + 0.05;
        let cap = 16;
        Ok(Self {
            n_sites,
            bonds,
            jz,
            beta,
            shift,
            spins,
            ops: vec![0; cap],
            n_ops: 0,
            rng,
            link: Vec::new(),
            legspin: Vec::new(),
            first: vec![NONE; n_sites],
            last: vec![NONE; n_sites],
            colour,
        })
    }

    /// Number of bonds.
    #[must_use]
    pub const fn n_bonds(&self) -> usize {
        self.bonds.len()
    }

    /// Diagonal-operator constant `C`.
    #[must_use]
    pub const fn constant(&self) -> f64 {
        self.shift
    }

    /// Weight of a vertex with leg spins `(in_i, in_j, out_i, out_j)` (each `+-1`).
    fn weight(
        &self,
        s: [i8; 4],
    ) -> f64 {
        let [a, b, c, d] = s;
        if a + b != c + d {
            return 0.0;
        }
        if a == c && b == d {
            if a == b {
                self.shift - 0.25 * self.jz
            } else {
                self.shift + 0.25 * self.jz
            }
        } else if a == d && b == c && a != b {
            0.5
        } else {
            0.0
        }
    }

    fn diag_update(&mut self) {
        let nb = self.bonds.len() as f64;
        let m = self.ops.len();
        for p in 0..m {
            let op = self.ops[p];
            if op == 0 {
                let b = self.rng.below(self.bonds.len() as u64) as usize;
                let (i, j) = self.bonds[b];
                let w = self.weight([self.spins[i], self.spins[j], self.spins[i], self.spins[j]]);
                let prob = nb * self.beta * w / (m - self.n_ops) as f64;
                if prob >= 1.0 || self.rng.uniform() < prob {
                    self.ops[p] = (b as u32 + 1) * 2;
                    self.n_ops += 1;
                }
            } else if op.is_multiple_of(2) {
                let b = (op / 2 - 1) as usize;
                let (i, j) = self.bonds[b];
                let w = self.weight([self.spins[i], self.spins[j], self.spins[i], self.spins[j]]);
                let prob = (m - self.n_ops + 1) as f64 / (nb * self.beta * w);
                if prob >= 1.0 || self.rng.uniform() < prob {
                    self.ops[p] = 0;
                    self.n_ops -= 1;
                }
            } else {
                let b = (op / 2 - 1) as usize;
                let (i, j) = self.bonds[b];
                self.spins[i] = -self.spins[i];
                self.spins[j] = -self.spins[j];
            }
        }
    }

    fn build_links(&mut self) {
        let m = self.ops.len();
        self.link.clear();
        self.link.resize(4 * m, NONE);
        self.legspin.clear();
        self.legspin.resize(4 * m, 0);
        self.first.fill(NONE);
        self.last.fill(NONE);
        let mut sigma = self.spins.clone();
        for p in 0..m {
            let op = self.ops[p];
            if op == 0 {
                continue;
            }
            let b = (op / 2 - 1) as usize;
            let (i, j) = self.bonds[b];
            let v = 4 * p;
            self.legspin[v] = sigma[i];
            self.legspin[v + 1] = sigma[j];
            if op % 2 == 1 {
                sigma[i] = -sigma[i];
                sigma[j] = -sigma[j];
            }
            self.legspin[v + 2] = sigma[i];
            self.legspin[v + 3] = sigma[j];
            for (site, leg_in, leg_out) in [(i, v, v + 2), (j, v + 1, v + 3)] {
                if self.last[site] == NONE {
                    self.first[site] = leg_in;
                } else {
                    let l = self.last[site];
                    self.link[l] = leg_in;
                    self.link[leg_in] = l;
                }
                self.last[site] = leg_out;
            }
        }
        for site in 0..self.n_sites {
            if self.first[site] != NONE {
                let (f, l) = (self.first[site], self.last[site]);
                self.link[l] = f;
                self.link[f] = l;
            }
        }
    }

    /// One loop-update pass; returns the number of vertex legs visited.
    fn loop_update(&mut self) -> usize {
        self.build_links();
        let m = self.ops.len();
        if self.n_ops == 0 {
            return 0;
        }
        let target = 2 * self.n_ops.max(8);
        let mut visited = 0;
        while visited < target {
            // random leg of a random existing operator
            let v0 = loop {
                let v = self.rng.below(4 * m as u64) as usize;
                if self.ops[v / 4] != 0 {
                    break v;
                }
            };
            let mut v = v0;
            loop {
                let p = v / 4;
                let l_in = v % 4;
                let mut s = [
                    self.legspin[4 * p],
                    self.legspin[4 * p + 1],
                    self.legspin[4 * p + 2],
                    self.legspin[4 * p + 3],
                ];
                s[l_in] = -s[l_in];
                let mut w = [0.0; 4];
                for (l_out, wv) in w.iter_mut().enumerate() {
                    let mut t = s;
                    t[l_out] = -t[l_out];
                    *wv = self.weight(t);
                }
                let total: f64 = w.iter().sum();
                let mut r = self.rng.uniform() * total;
                let mut l_out = 3;
                for (k, wv) in w.iter().enumerate() {
                    if r < *wv {
                        l_out = k;
                        break;
                    }
                    r -= wv;
                }
                s[l_out] = -s[l_out];
                for k in 0..4 {
                    self.legspin[4 * p + k] = s[k];
                }
                visited += 1;
                v = self.link[4 * p + l_out];
                if v == v0 {
                    break;
                }
            }
        }
        // write back operator types and the reference state
        for p in 0..m {
            if self.ops[p] == 0 {
                continue;
            }
            let b = self.ops[p] / 2 - 1;
            let diagonal = self.legspin[4 * p] == self.legspin[4 * p + 2];
            self.ops[p] = b * 2 + 2 + u32::from(!diagonal);
        }
        for site in 0..self.n_sites {
            if self.first[site] == NONE {
                if self.rng.uniform() < 0.5 {
                    self.spins[site] = -self.spins[site];
                }
            } else {
                self.spins[site] = self.legspin[self.first[site]];
            }
        }
        visited
    }

    /// One Monte Carlo sweep (diagonal update, directed loops, cutoff adjustment).
    pub fn sweep(&mut self) {
        self.diag_update();
        self.loop_update();
        let want = self.n_ops + self.n_ops / 3 + 8;
        if want > self.ops.len() {
            self.ops.resize(want, 0);
        }
    }

    /// Measurements of the current configuration.
    #[must_use]
    pub fn sample(&self) -> SseSample {
        let mag: f64 = self.spins.iter().map(|&s| 0.5 * f64::from(s)).sum();
        let stag: f64 = self
            .spins
            .iter()
            .zip(&self.colour)
            .map(|(&s, &c)| 0.5 * f64::from(s) * if c == 0 { 1.0 } else { -1.0 })
            .sum();
        SseSample {
            n_ops: self.n_ops,
            magnetisation: mag,
            staggered: stag,
        }
    }

    /// Runs `therm` thermalisation sweeps and `sweeps` measured sweeps and returns the averages
    /// `(energy, specific_heat, <m^2>, <m_stag^2>)` with the energy statistical error estimate.
    #[must_use]
    pub fn run(
        &mut self,
        therm: usize,
        sweeps: usize,
    ) -> SseResult {
        for _ in 0..therm {
            self.sweep();
        }
        let nb = self.bonds.len() as f64;
        let (mut sn, mut sn2, mut sm2, mut ss2) = (0.0, 0.0, 0.0, 0.0);
        let mut series = Vec::with_capacity(sweeps);
        for _ in 0..sweeps {
            self.sweep();
            let s = self.sample();
            let n = s.n_ops as f64;
            sn += n;
            sn2 += n * n;
            sm2 += s.magnetisation * s.magnetisation;
            ss2 += s.staggered * s.staggered;
            series.push(n);
        }
        let cnt = sweeps.max(1) as f64;
        let nbar = sn / cnt;
        let (n2_mean, nbar_sq) = (sn2 / cnt, nbar * nbar);
        let energy = nb * self.shift - nbar / self.beta;
        // blocked error of <n>
        let blocks = 20.min(series.len().max(1));
        let bl = series.len() / blocks.max(1);
        let means: Vec<f64> = (0..blocks)
            .map(|k| series[k * bl..(k + 1) * bl].iter().sum::<f64>() / bl.max(1) as f64)
            .collect();
        let mm = means.iter().sum::<f64>() / blocks as f64;
        let var = means.iter().map(|x| (x - mm).powi(2)).sum::<f64>() / (blocks.max(2) - 1) as f64;
        SseResult {
            energy,
            energy_error: (var / blocks as f64).sqrt() / self.beta,
            specific_heat: n2_mean - nbar_sq - nbar,
            magnetisation_sq: sm2 / cnt,
            staggered_sq: ss2 / cnt,
            mean_ops: nbar,
        }
    }
}

/// Averages of an SSE run.
#[derive(Clone, Copy, Debug)]
pub struct SseResult {
    /// Energy `<H>`.
    pub energy: f64,
    /// Statistical error of the energy (blocking).
    pub energy_error: f64,
    /// Specific heat (total, not per site).
    pub specific_heat: f64,
    /// `<(sum S^z)^2>`.
    pub magnetisation_sq: f64,
    /// `<(sum (-1)^i S^z_i)^2>`.
    pub staggered_sq: f64,
    /// Mean expansion order.
    pub mean_ops: f64,
}
