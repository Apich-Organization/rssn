//! # Lattice spin models: Ising, Potts, XY, Heisenberg in 1-3 (any) dimensions
//!
//! Monte Carlo for `H = -J sum_<ij> s_i . s_j - h sum_i s_i[0]` (O(N) models:
//! N = 1 Ising, 2 XY, 3 Heisenberg) and the q-state Potts model
//! `H = -J sum_<ij> delta(s_i, s_j)` on hypercubic lattices of any
//! dimension `D` (const generic), with periodic or open boundaries.
//!
//! | Algorithm | Ising (O(1)) | XY / Heisenberg (O(N), N>=2) | Potts |
//! |---|---|---|---|
//! | Metropolis | yes | yes (random-walk proposal) | yes |
//! | Heat bath | yes | no (needs von Mises / Gibbs sampling) | yes |
//! | Wolff cluster (embedded Ising for O(N)), `h = 0`, `J > 0` | yes | yes | yes |
//! | Swendsen-Wang (embedded) | yes | yes | yes |
//! | Over-relaxation (micro-canonical) | - | yes | - |
//!
//! Further tools: [`ParallelTempering`] (replica exchange), [`wang_landau_ising`]
//! (density of states), single-histogram reweighting ([`reweight`],
//! [`canonical_from_dos`]), finite-size-scaling helpers
//! ([`binder_crossing`], [`fit_exponent`], [`collapse_variable`]), blocking
//! error estimates, and exact references: 1-D chain ([`ising1d_energy_per_site`]),
//! 2-D Onsager ([`onsager_tc`], [`onsager_energy_per_site`],
//! [`onsager_magnetisation`]) and brute-force enumeration
//! ([`enumerate_ising_dos`]).
//!
//! Boundary / initial conditions: periodic and open boundaries (lattice
//! sides >= 3 when periodic), hot (random) and cold (aligned) starts, user
//! supplied spin arrays via the public `spins` field.
//!
//! Limitations: only nearest-neighbour hypercubic lattices; the field `h`
//! couples to component 0 and is incompatible with cluster updates.
#![allow(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::needless_range_loop,
    clippy::many_single_char_names,
    clippy::too_many_arguments,
    clippy::too_many_lines,
    clippy::missing_const_for_fn,
    clippy::float_cmp,
    clippy::single_call_fn,
    clippy::module_name_repetitions,
    clippy::cognitive_complexity,
    clippy::type_complexity
)]

use crate::kernels::random::Rng;

/// Hypercubic lattice with nearest-neighbour table.
#[derive(Clone, Debug)]
pub struct Lattice<const D: usize> {
    /// Side lengths.
    pub dims: [usize; D],
    /// Whether all boundaries are periodic.
    pub periodic: bool,
    nbr: Vec<usize>,
    n: usize,
}

impl<const D: usize> Lattice<D> {
    /// Builds a lattice; periodic lattices need every side `>= 3`.
    ///
    /// # Errors
    /// Returns a message for empty or too small sides.
    pub fn new(
        dims: [usize; D],
        periodic: bool,
    ) -> Result<Self, String> {
        if D == 0 {
            return Err("dimension must be positive".into());
        }
        let min = if periodic { 3 } else { 1 };
        if dims.iter().any(|&l| l < min) {
            return Err(format!("every side must be >= {min}"));
        }
        let n: usize = dims.iter().product();
        let mut stride = [1usize; D];
        for k in 1..D {
            stride[k] = stride[k - 1] * dims[k - 1];
        }
        let mut nbr = vec![usize::MAX; n * 2 * D];
        for i in 0..n {
            for k in 0..D {
                let x = (i / stride[k]) % dims[k];
                let plus = if x + 1 < dims[k] {
                    Some(i + stride[k])
                } else if periodic {
                    Some(i - x * stride[k])
                } else {
                    None
                };
                let minus = if x > 0 {
                    Some(i - stride[k])
                } else if periodic {
                    Some(i + (dims[k] - 1) * stride[k])
                } else {
                    None
                };
                if let Some(p) = plus {
                    nbr[i * 2 * D + 2 * k] = p;
                }
                if let Some(m) = minus {
                    nbr[i * 2 * D + 2 * k + 1] = m;
                }
            }
        }
        Ok(Self { dims, periodic, nbr, n })
    }

    /// Number of sites.
    #[must_use]
    pub fn sites(&self) -> usize {
        self.n
    }

    /// Neighbour slice of a site (`usize::MAX` marks a missing open-boundary neighbour).
    #[must_use]
    pub fn neighbours(
        &self,
        i: usize,
    ) -> &[usize] {
        &self.nbr[i * 2 * D..(i + 1) * 2 * D]
    }

    /// Number of bonds.
    #[must_use]
    pub fn bonds(&self) -> usize {
        self.nbr.iter().filter(|&&j| j != usize::MAX).count() / 2
    }

    fn for_each_bond(
        &self,
        mut f: impl FnMut(usize, usize),
    ) {
        for i in 0..self.n {
            for k in 0..D {
                let j = self.nbr[i * 2 * D + 2 * k];
                if j != usize::MAX {
                    f(i, j);
                }
            }
        }
    }
}

fn dot<const N: usize>(
    a: &[f64; N],
    b: &[f64; N],
) -> f64 {
    let mut s = 0.0;
    for k in 0..N {
        s += a[k] * b[k];
    }
    s
}

fn random_unit<const N: usize>(rng: &mut Rng) -> [f64; N] {
    loop {
        let mut v = [0.0; N];
        for c in &mut v {
            *c = rng.normal();
        }
        let n = dot(&v, &v).sqrt();
        if n > 1e-12 {
            for c in &mut v {
                *c /= n;
            }
            return v;
        }
    }
}

/// Update algorithm selector for the generic drivers.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Update {
    /// Single-spin Metropolis sweep with the given proposal width (ignored for Ising/Potts).
    Metropolis(f64),
    /// Heat-bath sweep.
    HeatBath,
    /// One Wolff cluster.
    Wolff,
    /// One Swendsen-Wang sweep.
    SwendsenWang,
}

struct UnionFind {
    parent: Vec<usize>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        Self { parent: (0..n).collect() }
    }

    fn find(
        &mut self,
        mut x: usize,
    ) -> usize {
        while self.parent[x] != x {
            self.parent[x] = self.parent[self.parent[x]];
            x = self.parent[x];
        }
        x
    }

    fn union(
        &mut self,
        a: usize,
        b: usize,
    ) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra != rb {
            self.parent[ra] = rb;
        }
    }
}

/// O(N) classical spin model on a `D`-dimensional lattice.
#[derive(Clone, Debug)]
pub struct OnModel<const D: usize, const N: usize> {
    /// Lattice geometry.
    pub lattice: Lattice<D>,
    /// Unit spin vectors (for `N = 1`: `+-1`).
    pub spins: Vec<[f64; N]>,
    /// Inverse temperature.
    pub beta: f64,
    /// Exchange coupling.
    pub j: f64,
    /// Field along component 0.
    pub field: f64,
    /// Random number generator.
    pub rng: Rng,
}

impl<const D: usize, const N: usize> OnModel<D, N> {
    /// Creates a model with a hot (random) start.
    ///
    /// # Errors
    /// Propagates lattice construction errors.
    pub fn new(
        dims: [usize; D],
        periodic: bool,
        beta: f64,
        seed: u64,
    ) -> Result<Self, String> {
        let lattice = Lattice::new(dims, periodic)?;
        let mut rng = Rng::new(seed);
        let spins = (0..lattice.sites())
            .map(|_| {
                if N == 1 {
                    let mut s = [0.0; N];
                    s[0] = if rng.uniform() < 0.5 { 1.0 } else { -1.0 };
                    s
                } else {
                    random_unit::<N>(&mut rng)
                }
            })
            .collect();
        Ok(Self { lattice, spins, beta, j: 1.0, field: 0.0, rng })
    }

    /// Aligns all spins along component 0 (cold start).
    pub fn cold_start(&mut self) {
        for s in &mut self.spins {
            *s = [0.0; N];
            s[0] = 1.0;
        }
    }

    /// Total energy.
    #[must_use]
    pub fn energy(&self) -> f64 {
        let mut e = 0.0;
        self.lattice.for_each_bond(|i, j| e -= self.j * dot(&self.spins[i], &self.spins[j]));
        for s in &self.spins {
            e -= self.field * s[0];
        }
        e
    }

    /// Total magnetisation vector.
    #[must_use]
    pub fn magnetisation(&self) -> [f64; N] {
        let mut m = [0.0; N];
        for s in &self.spins {
            for k in 0..N {
                m[k] += s[k];
            }
        }
        m
    }

    /// Magnetisation per site (norm).
    #[must_use]
    pub fn abs_magnetisation_per_site(&self) -> f64 {
        let m = self.magnetisation();
        dot(&m, &m).sqrt() / self.lattice.sites() as f64
    }

    fn local_field(
        &self,
        i: usize,
    ) -> [f64; N] {
        let mut f = [0.0; N];
        f[0] = self.field;
        for &j in self.lattice.neighbours(i) {
            if j != usize::MAX {
                for k in 0..N {
                    f[k] += self.j * self.spins[j][k];
                }
            }
        }
        f
    }

    /// One Metropolis sweep (`N` randomly chosen sites); returns the acceptance rate.
    pub fn metropolis_sweep(
        &mut self,
        delta: f64,
    ) -> f64 {
        let n = self.lattice.sites();
        let mut acc = 0usize;
        for _ in 0..n {
            // random site selection: sequential sweeps make domain walls move ballistically in 1-D
            let i = self.rng.below(n as u64) as usize;
            let f = self.local_field(i);
            let old = self.spins[i];
            let mut new = old;
            if N == 1 {
                new[0] = -old[0];
            } else {
                for k in 0..N {
                    new[k] = old[k] + delta * self.rng.normal();
                }
                let nn = dot(&new, &new).sqrt();
                for c in &mut new {
                    *c /= nn;
                }
            }
            let mut de = 0.0;
            for k in 0..N {
                de -= (new[k] - old[k]) * f[k];
            }
            if de <= 0.0 || self.rng.uniform() < (-self.beta * de).exp() {
                self.spins[i] = new;
                acc += 1;
            }
        }
        acc as f64 / n as f64
    }

    /// One heat-bath sweep (Ising only).
    ///
    /// # Errors
    /// Returns an error for `N != 1`.
    pub fn heat_bath_sweep(&mut self) -> Result<(), String> {
        if N != 1 {
            return Err("heat bath is implemented for N = 1 (Ising) only".into());
        }
        for i in 0..self.lattice.sites() {
            let f = self.local_field(i)[0];
            let p_up = 1.0 / (1.0 + (-2.0 * self.beta * f).exp());
            self.spins[i][0] = if self.rng.uniform() < p_up { 1.0 } else { -1.0 };
        }
        Ok(())
    }

    /// One sweep of micro-canonical over-relaxation (`N >= 2`): reflects every spin about its local field.
    pub fn overrelax_sweep(&mut self) {
        for i in 0..self.lattice.sites() {
            let f = self.local_field(i);
            let ff = dot(&f, &f);
            if ff > 1e-24 {
                let s = self.spins[i];
                let c = 2.0 * dot(&s, &f) / ff;
                for k in 0..N {
                    self.spins[i][k] = c * f[k] - s[k];
                }
            }
        }
    }

    fn cluster_ready(&self) -> Result<(), String> {
        if self.field != 0.0 {
            return Err("cluster updates require zero field".into());
        }
        if self.j <= 0.0 {
            return Err("cluster updates require a ferromagnetic coupling".into());
        }
        Ok(())
    }

    /// Grows and flips one Wolff cluster; returns its size.
    ///
    /// # Errors
    /// Fails for non-zero field or `J <= 0`.
    pub fn wolff_step(&mut self) -> Result<usize, String> {
        self.cluster_ready()?;
        let r = random_unit::<N>(&mut self.rng);
        let n = self.lattice.sites();
        let seed = self.rng.below(n as u64) as usize;
        let mut in_cluster = vec![false; n];
        let mut stack = vec![seed];
        in_cluster[seed] = true;
        let mut size = 0;
        while let Some(i) = stack.pop() {
            let pi = dot(&self.spins[i], &r);
            // flip i
            for k in 0..N {
                self.spins[i][k] -= 2.0 * pi * r[k];
            }
            size += 1;
            for t in 0..2 * D {
                let j = self.lattice.nbr[i * 2 * D + t];
                if j == usize::MAX || in_cluster[j] {
                    continue;
                }
                let pj = dot(&self.spins[j], &r);
                // after flipping i its projection is -pi; bond uses the pre-flip product
                let prod = pi * pj;
                if prod > 0.0 && self.rng.uniform() < 1.0 - (-2.0 * self.beta * self.j * prod).exp() {
                    in_cluster[j] = true;
                    stack.push(j);
                }
            }
        }
        Ok(size)
    }

    /// One Swendsen-Wang sweep with a random embedding direction; returns the number of clusters.
    ///
    /// # Errors
    /// Fails for non-zero field or `J <= 0`.
    pub fn swendsen_wang_step(&mut self) -> Result<usize, String> {
        self.cluster_ready()?;
        let r = random_unit::<N>(&mut self.rng);
        let n = self.lattice.sites();
        let proj: Vec<f64> = self.spins.iter().map(|s| dot(s, &r)).collect();
        let mut uf = UnionFind::new(n);
        let (beta, jj) = (self.beta, self.j);
        let mut links = Vec::new();
        self.lattice.for_each_bond(|i, j| links.push((i, j)));
        for (i, j) in links {
            let prod = proj[i] * proj[j];
            if prod > 0.0 && self.rng.uniform() < 1.0 - (-2.0 * beta * jj * prod).exp() {
                uf.union(i, j);
            }
        }
        let mut flip = vec![None; n];
        let mut clusters = 0;
        for i in 0..n {
            let root = uf.find(i);
            let f = match flip[root] {
                Some(f) => f,
                None => {
                    let f = self.rng.uniform() < 0.5;
                    flip[root] = Some(f);
                    clusters += 1;
                    f
                },
            };
            if f {
                for k in 0..N {
                    self.spins[i][k] -= 2.0 * proj[i] * r[k];
                }
            }
        }
        Ok(clusters)
    }

    /// Applies one update of the chosen kind (one sweep, or one Wolff cluster).
    ///
    /// # Errors
    /// Propagates unsupported-combination errors.
    pub fn update(
        &mut self,
        u: Update,
    ) -> Result<(), String> {
        match u {
            Update::Metropolis(d) => {
                self.metropolis_sweep(d);
            },
            Update::HeatBath => self.heat_bath_sweep()?,
            Update::Wolff => {
                self.wolff_step()?;
            },
            Update::SwendsenWang => {
                self.swendsen_wang_step()?;
            },
        }
        Ok(())
    }
}

/// q-state Potts model.
#[derive(Clone, Debug)]
pub struct Potts<const D: usize> {
    /// Lattice geometry.
    pub lattice: Lattice<D>,
    /// Number of states.
    pub q: usize,
    /// State of every site.
    pub spins: Vec<usize>,
    /// Inverse temperature.
    pub beta: f64,
    /// Coupling.
    pub j: f64,
    /// Random number generator.
    pub rng: Rng,
}

impl<const D: usize> Potts<D> {
    /// Creates a Potts model with a random start.
    ///
    /// # Errors
    /// Propagates lattice errors; requires `q >= 2`.
    pub fn new(
        dims: [usize; D],
        periodic: bool,
        q: usize,
        beta: f64,
        seed: u64,
    ) -> Result<Self, String> {
        if q < 2 {
            return Err("q must be at least 2".into());
        }
        let lattice = Lattice::new(dims, periodic)?;
        let mut rng = Rng::new(seed);
        let spins = (0..lattice.sites()).map(|_| rng.below(q as u64) as usize).collect();
        Ok(Self { lattice, q, spins, beta, j: 1.0, rng })
    }

    /// Total energy `-J * (number of equal bonds)`.
    #[must_use]
    pub fn energy(&self) -> f64 {
        let mut c = 0usize;
        self.lattice.for_each_bond(|i, j| {
            if self.spins[i] == self.spins[j] {
                c += 1;
            }
        });
        -self.j * c as f64
    }

    /// Order parameter `(q max_s n_s / N - 1)/(q - 1)`.
    #[must_use]
    pub fn order_parameter(&self) -> f64 {
        let mut cnt = vec![0usize; self.q];
        for &s in &self.spins {
            cnt[s] += 1;
        }
        let mx = cnt.iter().copied().max().unwrap_or(0) as f64;
        (self.q as f64 * mx / self.spins.len() as f64 - 1.0) / (self.q as f64 - 1.0)
    }

    fn equal_neighbours(
        &self,
        i: usize,
        state: usize,
    ) -> usize {
        self.lattice.neighbours(i).iter().filter(|&&j| j != usize::MAX && self.spins[j] == state).count()
    }

    /// Metropolis sweep; returns the acceptance rate.
    pub fn metropolis_sweep(&mut self) -> f64 {
        let n = self.spins.len();
        let mut acc = 0;
        for _ in 0..n {
            let i = self.rng.below(n as u64) as usize;
            let old = self.spins[i];
            let mut new = self.rng.below(self.q as u64 - 1) as usize;
            if new >= old {
                new += 1;
            }
            let de = -self.j * (self.equal_neighbours(i, new) as f64 - self.equal_neighbours(i, old) as f64);
            if de <= 0.0 || self.rng.uniform() < (-self.beta * de).exp() {
                self.spins[i] = new;
                acc += 1;
            }
        }
        acc as f64 / n as f64
    }

    /// Heat-bath sweep.
    pub fn heat_bath_sweep(&mut self) {
        let mut w = vec![0.0; self.q];
        for i in 0..self.spins.len() {
            let mut tot = 0.0;
            for s in 0..self.q {
                w[s] = (self.beta * self.j * self.equal_neighbours(i, s) as f64).exp();
                tot += w[s];
            }
            let mut u = self.rng.uniform() * tot;
            let mut pick = self.q - 1;
            for s in 0..self.q {
                if u < w[s] {
                    pick = s;
                    break;
                }
                u -= w[s];
            }
            self.spins[i] = pick;
        }
    }

    /// One Wolff cluster (`J > 0`); returns its size.
    pub fn wolff_step(&mut self) -> usize {
        let n = self.spins.len();
        let seed = self.rng.below(n as u64) as usize;
        let old = self.spins[seed];
        let mut new = self.rng.below(self.q as u64 - 1) as usize;
        if new >= old {
            new += 1;
        }
        let p = 1.0 - (-self.beta * self.j).exp();
        let mut stack = vec![seed];
        self.spins[seed] = new;
        let mut size = 1;
        while let Some(i) = stack.pop() {
            for t in 0..2 * D {
                let j = self.lattice.nbr[i * 2 * D + t];
                if j != usize::MAX && self.spins[j] == old && self.rng.uniform() < p {
                    self.spins[j] = new;
                    stack.push(j);
                    size += 1;
                }
            }
        }
        size
    }

    /// One Swendsen-Wang sweep (`J > 0`); returns the number of clusters.
    pub fn swendsen_wang_step(&mut self) -> usize {
        let n = self.spins.len();
        let p = 1.0 - (-self.beta * self.j).exp();
        let mut uf = UnionFind::new(n);
        let mut links = Vec::new();
        self.lattice.for_each_bond(|i, j| links.push((i, j)));
        for (i, j) in links {
            if self.spins[i] == self.spins[j] && self.rng.uniform() < p {
                uf.union(i, j);
            }
        }
        let mut state = vec![usize::MAX; n];
        let mut clusters = 0;
        for i in 0..n {
            let r = uf.find(i);
            if state[r] == usize::MAX {
                state[r] = self.rng.below(self.q as u64) as usize;
                clusters += 1;
            }
            self.spins[i] = state[r];
        }
        clusters
    }

    /// Applies one update of the chosen kind.
    pub fn update(
        &mut self,
        u: Update,
    ) {
        match u {
            Update::Metropolis(_) => {
                self.metropolis_sweep();
            },
            Update::HeatBath => self.heat_bath_sweep(),
            Update::Wolff => {
                self.wolff_step();
            },
            Update::SwendsenWang => {
                self.swendsen_wang_step();
            },
        }
    }
}

// ---------------------------------------------------------------------------
// Parallel tempering
// ---------------------------------------------------------------------------

/// Replica exchange over a ladder of inverse temperatures for O(N) models.
#[derive(Clone, Debug)]
pub struct ParallelTempering<const D: usize, const N: usize> {
    /// One replica per temperature (index = temperature slot, `betas[k]` ascending or any order).
    pub replicas: Vec<OnModel<D, N>>,
    /// Inverse temperatures.
    pub betas: Vec<f64>,
    /// Number of attempted/accepted swaps per neighbouring pair.
    pub swap_stats: Vec<(usize, usize)>,
    rng: Rng,
}

impl<const D: usize, const N: usize> ParallelTempering<D, N> {
    /// Creates replicas at the given inverse temperatures.
    ///
    /// # Errors
    /// Propagates lattice errors; requires at least two temperatures.
    pub fn new(
        dims: [usize; D],
        periodic: bool,
        betas: &[f64],
        seed: u64,
    ) -> Result<Self, String> {
        if betas.len() < 2 {
            return Err("need at least two temperatures".into());
        }
        let mut replicas = Vec::new();
        for (k, &b) in betas.iter().enumerate() {
            replicas.push(OnModel::new(dims, periodic, b, seed.wrapping_add(1000 + k as u64))?);
        }
        Ok(Self { replicas, betas: betas.to_vec(), swap_stats: vec![(0, 0); betas.len() - 1], rng: Rng::new(seed ^ 0xABCD) })
    }

    /// Performs `rounds` rounds of (`update` on every replica, then an exchange sweep).
    ///
    /// # Errors
    /// Propagates unsupported update errors.
    pub fn run(
        &mut self,
        rounds: usize,
        update: Update,
        mut measure: impl FnMut(usize, &[OnModel<D, N>]),
    ) -> Result<(), String> {
        for round in 0..rounds {
            for r in &mut self.replicas {
                r.update(update)?;
            }
            let start = round % 2;
            let mut k = start;
            while k + 1 < self.replicas.len() {
                let (e1, e2) = (self.replicas[k].energy(), self.replicas[k + 1].energy());
                let d = (self.betas[k] - self.betas[k + 1]) * (e1 - e2);
                self.swap_stats[k].0 += 1;
                if d >= 0.0 || self.rng.uniform() < d.exp() {
                    // exchange configurations, keep betas attached to slots
                    let (a, b) = self.replicas.split_at_mut(k + 1);
                    std::mem::swap(&mut a[k].spins, &mut b[0].spins);
                    self.swap_stats[k].1 += 1;
                }
                k += 2;
            }
            measure(round, &self.replicas);
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Density of states, reweighting
// ---------------------------------------------------------------------------

/// Brute-force density of states of the Ising model `E = -sum s_i s_j` by enumerating all
/// `2^N` states (`N <= 26`). Returns `(energy, count)` pairs sorted by energy.
///
/// # Errors
/// Fails for lattices larger than 26 sites or invalid geometry.
pub fn enumerate_ising_dos<const D: usize>(
    dims: [usize; D],
    periodic: bool,
) -> Result<Vec<(i64, f64)>, String> {
    let lat = Lattice::new(dims, periodic)?;
    let n = lat.sites();
    if n > 26 {
        return Err("enumeration limited to 26 sites".into());
    }
    let mut bonds = Vec::new();
    lat.for_each_bond(|i, j| bonds.push((i, j)));
    let mut map = std::collections::BTreeMap::new();
    for state in 0u64..(1u64 << n) {
        let mut e = 0i64;
        for &(i, j) in &bonds {
            let si = ((state >> i) & 1) as i64 * 2 - 1;
            let sj = ((state >> j) & 1) as i64 * 2 - 1;
            e -= si * sj;
        }
        *map.entry(e).or_insert(0.0) += 1.0;
    }
    Ok(map.into_iter().collect())
}

/// Wang-Landau estimate of the Ising density of states, normalised to `sum g = 2^N`.
///
/// The move set is single spin flips; the flatness criterion is `min H > flatness * mean H`
/// over the energies visited so far; `ln f` starts at 1 and is halved until `ln_f_final`.
///
/// # Errors
/// Propagates lattice errors.
pub fn wang_landau_ising<const D: usize>(
    dims: [usize; D],
    periodic: bool,
    ln_f_final: f64,
    flatness: f64,
    seed: u64,
) -> Result<Vec<(i64, f64)>, String> {
    let lat = Lattice::new(dims, periodic)?;
    let n = lat.sites();
    let nb = lat.bonds() as i64;
    // energy index (E + nb) / 2
    let levels = (nb + 1) as usize;
    let mut rng = Rng::new(seed);
    let mut spins: Vec<i64> = (0..n).map(|_| if rng.uniform() < 0.5 { 1 } else { -1 }).collect();
    let energy_of = |s: &[i64]| -> i64 {
        let mut e = 0;
        lat.for_each_bond(|i, j| e -= s[i] * s[j]);
        e
    };
    let mut e = energy_of(&spins);
    let mut ln_g = vec![0.0f64; levels];
    let mut visited = vec![false; levels];
    let mut ln_f = 1.0;
    let idx = |e: i64| ((e + nb) / 2) as usize;
    while ln_f > ln_f_final {
        let mut hist = vec![0usize; levels];
        let mut checks = 0usize;
        loop {
            for _ in 0..n * 100 {
                let i = rng.below(n as u64) as usize;
                let mut s = 0;
                for &j in lat.neighbours(i) {
                    if j != usize::MAX {
                        s += spins[j];
                    }
                }
                let e_new = e + 2 * spins[i] * s;
                let (a, b) = (idx(e), idx(e_new));
                if ln_g[a] - ln_g[b] >= 0.0 || rng.uniform() < (ln_g[a] - ln_g[b]).exp() {
                    spins[i] = -spins[i];
                    e = e_new;
                }
                let k = idx(e);
                ln_g[k] += ln_f;
                hist[k] += 1;
                visited[k] = true;
            }
            checks += 1;
            let vis: Vec<usize> = (0..levels).filter(|&k| visited[k]).map(|k| hist[k]).collect();
            let mean = vis.iter().sum::<usize>() as f64 / vis.len() as f64;
            let mn = vis.iter().copied().min().unwrap_or(0) as f64;
            if mn > flatness * mean && checks > 2 {
                break;
            }
            if checks > 200_000 {
                break;
            }
        }
        ln_f /= 2.0;
    }
    let mx = (0..levels).filter(|&k| visited[k]).map(|k| ln_g[k]).fold(f64::NEG_INFINITY, f64::max);
    let total: f64 = (0..levels).filter(|&k| visited[k]).map(|k| (ln_g[k] - mx).exp()).sum();
    let shift = (n as f64) * std::f64::consts::LN_2 - (mx + total.ln());
    Ok((0..levels).filter(|&k| visited[k]).map(|k| (2 * k as i64 - nb, ln_g[k] + shift)).collect())
}

/// Canonical quantities from a density of states `(E, ln g)`: returns `(ln Z, <E>, C)` (`C` is the heat capacity `beta^2 var(E)`).
#[must_use]
pub fn canonical_from_dos(
    dos: &[(f64, f64)],
    beta: f64,
) -> (f64, f64, f64) {
    let m = dos.iter().map(|&(e, lg)| lg - beta * e).fold(f64::NEG_INFINITY, f64::max);
    let (mut z, mut e1, mut e2) = (0.0, 0.0, 0.0);
    for &(e, lg) in dos {
        let w = (lg - beta * e - m).exp();
        z += w;
        e1 += w * e;
        e2 += w * e * e;
    }
    let me = e1 / z;
    (m + z.ln(), me, beta * beta * (e2 / z - me * me))
}

/// Single-histogram reweighting: weighted mean of `observable` at `beta1` from samples drawn at `beta0`
/// with energies `energies` (total energy per sample).
#[must_use]
pub fn reweight(
    energies: &[f64],
    observable: &[f64],
    beta0: f64,
    beta1: f64,
) -> f64 {
    let d = beta1 - beta0;
    let m = energies.iter().map(|&e| -d * e).fold(f64::NEG_INFINITY, f64::max);
    let (mut num, mut den) = (0.0, 0.0);
    for (&e, &o) in energies.iter().zip(observable) {
        let w = (-d * e - m).exp();
        num += w * o;
        den += w;
    }
    num / den
}

// ---------------------------------------------------------------------------
// Statistics and finite-size scaling
// ---------------------------------------------------------------------------

/// Mean and blocking (binning) error estimate of a correlated series.
#[must_use]
pub fn mean_and_error(x: &[f64]) -> (f64, f64) {
    let n = x.len();
    let mean = x.iter().sum::<f64>() / n as f64;
    let mut cur: Vec<f64> = x.to_vec();
    let mut best = 0.0_f64;
    while cur.len() >= 16 {
        let m = cur.len() as f64;
        let mu = cur.iter().sum::<f64>() / m;
        let var = cur.iter().map(|v| (v - mu).powi(2)).sum::<f64>() / (m - 1.0);
        best = best.max((var / m).sqrt());
        cur = cur.chunks_exact(2).map(|c| 0.5 * (c[0] + c[1])).collect();
    }
    (mean, best)
}

/// Binder cumulant `U4 = 1 - <m^4> / (3 <m^2>^2)` from magnetisation-per-site samples.
#[must_use]
pub fn binder_cumulant(m: &[f64]) -> f64 {
    let n = m.len() as f64;
    let m2 = m.iter().map(|v| v * v).sum::<f64>() / n;
    let m4 = m.iter().map(|v| v.powi(4)).sum::<f64>() / n;
    1.0 - m4 / (3.0 * m2 * m2)
}

/// Susceptibility per site `beta N (<m^2> - <|m|>^2)` from magnetisation-per-site samples.
#[must_use]
pub fn susceptibility(
    m: &[f64],
    beta: f64,
    sites: usize,
) -> f64 {
    let n = m.len() as f64;
    let m2 = m.iter().map(|v| v * v).sum::<f64>() / n;
    let m1 = m.iter().map(|v| v.abs()).sum::<f64>() / n;
    beta * sites as f64 * (m2 - m1 * m1)
}

/// Specific heat per site `beta^2 (<E^2> - <E>^2) / N` from total-energy samples.
#[must_use]
pub fn specific_heat(
    e: &[f64],
    beta: f64,
    sites: usize,
) -> f64 {
    let n = e.len() as f64;
    let mu = e.iter().sum::<f64>() / n;
    let var = e.iter().map(|v| (v - mu).powi(2)).sum::<f64>() / n;
    beta * beta * var / sites as f64
}

/// First crossing of two curves sampled at the same temperatures (linear interpolation).
#[must_use]
pub fn binder_crossing(
    t: &[f64],
    u_small: &[f64],
    u_large: &[f64],
) -> Option<f64> {
    for k in 0..t.len().saturating_sub(1) {
        let d0 = u_small[k] - u_large[k];
        let d1 = u_small[k + 1] - u_large[k + 1];
        if d0 == 0.0 {
            return Some(t[k]);
        }
        if d0 * d1 < 0.0 {
            return Some(t[k] + (t[k + 1] - t[k]) * d0 / (d0 - d1));
        }
    }
    None
}

/// Least-squares exponent `a` in `y ~ L^a` from `(L, y)` pairs.
#[must_use]
pub fn fit_exponent(
    sizes: &[f64],
    values: &[f64],
) -> f64 {
    let n = sizes.len() as f64;
    let xs: Vec<f64> = sizes.iter().map(|v| v.ln()).collect();
    let ys: Vec<f64> = values.iter().map(|v| v.ln()).collect();
    let (mx, my) = (xs.iter().sum::<f64>() / n, ys.iter().sum::<f64>() / n);
    let sxy: f64 = xs.iter().zip(&ys).map(|(x, y)| (x - mx) * (y - my)).sum();
    let sxx: f64 = xs.iter().map(|x| (x - mx).powi(2)).sum();
    sxy / sxx
}

/// Scaling variable `(T - Tc) L^(1/nu)` for data collapse.
#[must_use]
pub fn collapse_variable(
    t: f64,
    tc: f64,
    nu: f64,
    l: f64,
) -> f64 {
    (t - tc) * l.powf(1.0 / nu)
}

// ---------------------------------------------------------------------------
// Exact references
// ---------------------------------------------------------------------------

/// Critical temperature of the 2-D square-lattice Ising model (`k_B = 1`).
#[must_use]
pub fn onsager_tc(j: f64) -> f64 {
    2.0 * j / (1.0 + std::f64::consts::SQRT_2).ln()
}

fn ellip_k(k: f64) -> f64 {
    let (mut a, mut b) = (1.0, (1.0 - k * k).sqrt());
    for _ in 0..40 {
        let (an, bn) = (0.5 * (a + b), (a * b).sqrt());
        a = an;
        b = bn;
    }
    std::f64::consts::FRAC_PI_2 / a
}

/// Onsager internal energy per site of the 2-D square-lattice Ising model.
#[must_use]
pub fn onsager_energy_per_site(
    t: f64,
    j: f64,
) -> f64 {
    let k = j / t;
    let th = (2.0 * k).tanh();
    let modulus = 2.0 * (2.0 * k).sinh() / (2.0 * k).cosh().powi(2);
    -j / th * (1.0 + 2.0 / std::f64::consts::PI * (2.0 * th * th - 1.0) * ellip_k(modulus.min(1.0 - 1e-15)))
}

/// Onsager-Yang spontaneous magnetisation of the 2-D Ising model (`0` above `Tc`).
#[must_use]
pub fn onsager_magnetisation(
    t: f64,
    j: f64,
) -> f64 {
    if t >= onsager_tc(j) {
        0.0
    } else {
        (1.0 - (2.0 * j / t).sinh().powi(-4)).powf(0.125)
    }
}

/// Exact energy per site of the periodic 1-D Ising chain of `n` sites.
#[must_use]
pub fn ising1d_energy_per_site(
    beta: f64,
    j: f64,
    n: usize,
) -> f64 {
    let k = beta * j;
    let (c, s) = (2.0 * k.cosh(), 2.0 * k.sinh());
    let nn = n as f64;
    let z = c.powf(nn) + s.powf(nn);
    let dz = nn * c.powf(nn - 1.0) * s + nn * s.powf(nn - 1.0) * c;
    -j * dz / z / nn
}
