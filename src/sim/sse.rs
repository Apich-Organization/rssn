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
//!   flips)`, which satisfies the directed-loop equations (bounces are allowed).
//!   A nonempty sweep starts `max(8, N_sites)` worms. This count depends only on
//!   the lattice size. Sites untouched by the string are sampled as free spins,
//!   including every site when the expansion order is zero.
//!
//! The operator-string cutoff `M` grows as the expansion order increases, including during thermalisation.
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

// Native language calls bound random-walk work and operator storage together.
// The original direct Rust API remains unchanged unless a scoped budget is installed.
std::thread_local! {static NATIVE_SSE_BUDGET:std::cell::Cell<Option<(usize,usize)>>=const{std::cell::Cell::new(None)};}
struct SseBudgetExhausted;
pub(crate) fn with_native_budget<T>(
    work: usize,
    operators: usize,
    run: impl FnOnce() -> T,
) -> Result<T, String> {
    struct Restore(Option<(usize, usize)>);
    impl Drop for Restore {
        fn drop(&mut self) {
            NATIVE_SSE_BUDGET.set(self.0);
        }
    }
    let old = NATIVE_SSE_BUDGET.replace(Some((work, operators)));
    let _restore = Restore(old);
    match std::panic::catch_unwind(std::panic::AssertUnwindSafe(run)) {
        | Ok(v) => Ok(v),
        | Err(e) if e.is::<SseBudgetExhausted>() => {
            Err("SSE native operation/storage budget exhausted".into())
        },
        | Err(e) => std::panic::resume_unwind(e),
    }
}
fn budget_tick(
    operations: usize,
    storage: usize,
) {
    NATIVE_SSE_BUDGET.with(|slot| {
        if let Some((left, maximum)) = slot.get() {
            if let Some(next) = left.checked_sub(operations)
                && storage <= maximum
            {
                slot.set(Some((next, maximum)));
            } else {
                std::panic::resume_unwind(Box::new(SseBudgetExhausted));
            }
        }
    });
}

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

#[derive(Clone, Copy, Debug)]
struct LoopExits {
    weights: [f64; 4],
    total: f64,
}

#[derive(Clone, Debug)]
struct VertexWeights {
    jz_bits: u64,
    // The first weight is for opposite spins, the second for aligned spins.
    diagonal: [f64; 2],
    exits: [LoopExits; 64],
}

impl VertexWeights {
    fn new(
        jz: f64,
        shift: f64,
    ) -> Self {
        let diagonal = [shift + 0.25 * jz, shift - 0.25 * jz];
        let weights: [f64; 16] = std::array::from_fn(|mask| {
            let input = mask & 3;
            let output = mask >> 2;
            let aligned = input == 0 || input == 3;
            if input == output {
                diagonal[usize::from(aligned)]
            } else if !aligned && input ^ output == 3 {
                0.5
            } else {
                0.0
            }
        });
        let exits = std::array::from_fn(|index| {
            let open = (index / 4) ^ (1 << (index % 4));
            let weights = std::array::from_fn(|exit| weights[open ^ (1 << exit)]);
            LoopExits {
                total: weights.iter().sum(),
                weights,
            }
        });
        Self {
            jz_bits: jz.to_bits(),
            diagonal,
            exits,
        }
    }
}

/// SSE simulation of the spin-1/2 XXZ model.
#[derive(Clone, Debug)]
pub struct SseXxz {
    n_sites: usize,
    bonds: Vec<(usize, usize)>,
    /// Ising anisotropy coupling `J_z` (the transverse coupling is `J_xy = 1`).
    ///
    /// If changed between sweeps, both diagonal weights must remain nonnegative
    /// with the value returned by [`Self::constant`].
    pub jz: f64,
    /// Positive, finite inverse temperature.
    pub beta: f64,
    shift: f64,
    spins: Vec<i8>,
    ops: Vec<u32>,
    n_ops: usize,
    rng: Rng,
    link: Vec<usize>,
    vertices: Vec<u8>,
    propagated: Vec<i8>,
    weights: VertexWeights,
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
    /// Fails for an empty or unencodable bond list, invalid sites or
    /// couplings, a non-positive inverse temperature, or a graph that is
    /// not bipartite.
    pub fn new(
        n_sites: usize,
        bonds: Vec<(usize, usize)>,
        jz: f64,
        beta: f64,
        seed: u64,
    ) -> Result<Self, String> {
        if bonds.is_empty()
            || !jz.is_finite()
            || !beta.is_finite()
            || beta <= 0.0
            || bonds
                .iter()
                .any(|&(i, j)| i >= n_sites || j >= n_sites || i == j)
        {
            return Err("need valid bonds, finite coupling and finite beta > 0".into());
        }
        if bonds.len() > u32::MAX as usize / 2 {
            return Err("bond count exceeds the SSE operator encoding".into());
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
            vertices: Vec::new(),
            propagated: Vec::new(),
            weights: VertexWeights::new(jz, shift),
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

    fn diag_update(&mut self) {
        let nb = self.bonds.len() as f64;
        let m = self.ops.len();
        budget_tick(m, m);
        for p in 0..m {
            let op = self.ops[p];
            if op == 0 {
                let b = self.rng.below(self.bonds.len() as u64) as usize;
                let (i, j) = self.bonds[b];
                let w = self.weights.diagonal[usize::from(self.spins[i] == self.spins[j])];
                let prob = nb * self.beta * w / (m - self.n_ops) as f64;
                if prob >= 1.0 || self.rng.uniform() < prob {
                    self.ops[p] = (b as u32 + 1) * 2;
                    self.n_ops += 1;
                }
            } else if op.is_multiple_of(2) {
                let b = (op / 2 - 1) as usize;
                let (i, j) = self.bonds[b];
                let w = self.weights.diagonal[usize::from(self.spins[i] == self.spins[j])];
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
        budget_tick(m, m);
        // Each active leg is linked during this pass. Identity slots are
        // never read as vertices, so the rest of the allocated array can stay.
        self.link.resize(4 * m, NONE);
        self.vertices.resize(m, 0);
        self.first.fill(NONE);
        self.last.fill(NONE);
        self.propagated.clone_from(&self.spins);
        for p in 0..m {
            let op = self.ops[p];
            if op == 0 {
                continue;
            }
            let b = (op / 2 - 1) as usize;
            let (i, j) = self.bonds[b];
            let v = 4 * p;
            let input = u8::from(self.propagated[i] > 0) | (u8::from(self.propagated[j] > 0) << 1);
            let output = if op % 2 == 1 {
                self.propagated[i] = -self.propagated[i];
                self.propagated[j] = -self.propagated[j];
                input ^ 3
            } else {
                input
            };
            self.vertices[p] = input | (output << 2);
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

    fn loop_update(&mut self) {
        self.build_links();
        let m = self.ops.len();
        budget_tick(m, m);
        if self.n_ops != 0 {
            // A predetermined number of reversible worm transitions preserves
            // the distribution. A stop based on their lengths can bias it.
            for _ in 0..self.n_sites.max(8) {
                // random leg of a random existing operator
                let v0 = loop {
                    budget_tick(1, m);
                    let v = self.rng.below(4 * m as u64) as usize;
                    if self.ops[v / 4] != 0 {
                        break v;
                    }
                };
                let mut v = v0;
                loop {
                    budget_tick(1, m);
                    let p = v / 4;
                    let l_in = v % 4;
                    let mask = self.vertices[p];
                    let exits = &self.weights.exits[usize::from(mask) * 4 + l_in];
                    let mut r = self.rng.uniform() * exits.total;
                    let mut l_out = 3;
                    for (k, wv) in exits.weights.iter().enumerate() {
                        if r < *wv {
                            l_out = k;
                            break;
                        }
                        r -= wv;
                    }
                    self.vertices[p] = mask ^ (1 << l_in) ^ (1 << l_out);
                    v = self.link[4 * p + l_out];
                    if v == v0 {
                        break;
                    }
                }
            }
            // Write back operator types once all worms have closed.
            for p in 0..m {
                if self.ops[p] == 0 {
                    continue;
                }
                let b = self.ops[p] / 2 - 1;
                let vertex = self.vertices[p];
                let diagonal = vertex & 1 == (vertex >> 2) & 1;
                self.ops[p] = b * 2 + 2 + u32::from(!diagonal);
            }
        }

        budget_tick(self.n_sites, m);
        for site in 0..self.n_sites {
            if self.first[site] == NONE {
                if self.rng.uniform() < 0.5 {
                    self.spins[site] = -self.spins[site];
                }
            } else {
                let first = self.first[site];
                let bit = (self.vertices[first / 4] >> (first % 4)) & 1;
                self.spins[site] = if bit == 0 { -1 } else { 1 };
            }
        }
    }

    /// One Monte Carlo sweep (diagonal update, directed loops, cutoff adjustment).
    ///
    /// Fixed worm counts and the zero-order spin refresh reproduce the same
    /// trajectory for the same seed and inputs under this schedule. They
    /// change the sequences produced by earlier visit-capped releases.
    pub fn sweep(&mut self) {
        if self.weights.jz_bits != self.jz.to_bits() {
            self.weights = VertexWeights::new(self.jz, self.shift);
        }
        self.diag_update();
        self.loop_update();
        let want = self.n_ops + self.n_ops / 3 + 8;
        if want > self.ops.len() {
            budget_tick(0, want);
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
        // The error estimate uses only complete equal-size blocks; the
        // trailing remainder still contributes to the global estimators.
        let blocks = sweeps.clamp(1, 20);
        let block_len = sweeps / blocks;
        let covered = blocks * block_len;
        let mut block_sums = [0.0; 20];
        for index in 0..sweeps {
            self.sweep();
            let s = self.sample();
            let n = s.n_ops as f64;
            sn += n;
            sn2 += n * n;
            sm2 += s.magnetisation * s.magnetisation;
            ss2 += s.staggered * s.staggered;
            if index < covered {
                block_sums[index / block_len] += n;
            }
        }
        let cnt = sweeps.max(1) as f64;
        let nbar = sn / cnt;
        let (n2_mean, nbar_sq) = (sn2 / cnt, nbar * nbar);
        let energy = nb * self.shift - nbar / self.beta;
        let means = &mut block_sums[..blocks];
        for sum in means.iter_mut() {
            *sum /= block_len.max(1) as f64;
        }
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
