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
//! | Heat bath | yes | yes (exact: von Mises for XY, inverse-CDF `exp(a cos theta)` for Heisenberg, Wood's von Mises-Fisher algorithm for `N >= 4`; field included) | yes |
//! | Wolff cluster (embedded Ising for O(N)), `J > 0`, field through a ghost spin | yes | yes | yes |
//! | Swendsen-Wang (embedded; field through a ghost spin) | yes | yes | yes |
//! | Over-relaxation (micro-canonical) | - | yes | - |
//!
//! Further tools: [`ParallelTempering`] (replica exchange for O(N) models),
//! [`PottsTempering`] (replica exchange for Potts), [`replica_exchange_wang_landau`] (multi-window, multi-walker REWL for Ising and O(N)), [`wang_landau_ising`] (density of
//! states of the Ising model) and [`wang_landau_on`] (density of states of continuous O(N)
//! spins on a binned energy axis), single-histogram reweighting ([`reweight`],
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
//! couples to component 0 (cluster updates include it through a ghost spin); the quantum spin-1/2 XXZ
//! model (SSE with directed loops) lives in [`crate::sim::sse`]; single-walker Wang-Landau for continuous spins
//! ([`wang_landau_on`]) is restricted to a user-supplied energy window (moves leaving it are
//! rejected) and uses single-spin moves ([`replica_exchange_wang_landau`] adds windows and
//! walkers).

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
    pub const fn sites(&self) -> usize {
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

/// Draws a heat-bath direction for an O(N) spin.
///
/// The unit vector has density `P(s) ~ exp(beta s . f)` on the sphere `S^(N-1)`: the heat-bath conditional distribution of an O(N) spin in the local field `f`
/// (von Mises-Fisher distribution; `N = 1` logistic, `N = 2` Best-Fisher, `N = 3` exact inverse
/// CDF, `N >= 4` Wood's rejection algorithm).
///
/// # Errors
/// Fails for `N = 0`.
pub fn heat_bath_direction<const N: usize>(
    rng: &mut Rng,
    f: &[f64; N],
    beta: f64,
) -> Result<[f64; N], String> {
    let mut out = [0.0; N];
    let fm = dot(f, f).sqrt();
    let a = beta * fm;
    match N {
        1 => {
            let p_up = 1.0 / (1.0 + (-2.0 * beta * f[0]).exp());
            out[0] = if rng.uniform() < p_up { 1.0 } else { -1.0 };
        },
        2 => {
            let alpha = f[1].atan2(f[0]);
            let theta = if a < 1e-10 {
                (2.0 * rng.uniform() - 1.0) * std::f64::consts::PI
            } else {
                // Best-Fisher rejection sampler for the von Mises distribution
                let tau = 1.0 + (1.0 + 4.0 * a * a).sqrt();
                let rho = (tau - (2.0 * tau).sqrt()) / (2.0 * a);
                let r = (1.0 + rho * rho) / (2.0 * rho);
                loop {
                    let (u1, u2, u3) = (rng.uniform(), rng.uniform(), rng.uniform());
                    let z = (std::f64::consts::PI * u1).cos();
                    let fz = (1.0 + r * z) / (r + z);
                    let c = a * (r - fz);
                    if c * (2.0 - c) - u2 > 0.0 || (c / u2).ln() + 1.0 - c >= 0.0 {
                        break if u3 > 0.5 { fz.clamp(-1.0, 1.0).acos() } else { -fz.clamp(-1.0, 1.0).acos() };
                    }
                }
            };
            out[0] = (alpha + theta).cos();
            out[1] = (alpha + theta).sin();
        },
        3 => {
            let u = rng.uniform();
            let x = if a < 1e-10 { 2.0 * u - 1.0 } else { (1.0 + (u + (1.0 - u) * (-2.0 * a).exp()).ln() / a).clamp(-1.0, 1.0) };
            let sx = (1.0 - x * x).max(0.0).sqrt();
            let psi = 2.0 * std::f64::consts::PI * rng.uniform();
            let axis = if fm > 0.0 { [f[0] / fm, f[1] / fm, f[2] / fm] } else { [0.0, 0.0, 1.0] };
            let helper = if axis[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
            let mut e1 = [
                axis[1] * helper[2] - axis[2] * helper[1],
                axis[2] * helper[0] - axis[0] * helper[2],
                axis[0] * helper[1] - axis[1] * helper[0],
            ];
            let n1 = (e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]).sqrt();
            for c in &mut e1 {
                *c /= n1;
            }
            let e2 = [
                axis[1] * e1[2] - axis[2] * e1[1],
                axis[2] * e1[0] - axis[0] * e1[2],
                axis[0] * e1[1] - axis[1] * e1[0],
            ];
            for k in 0..3 {
                out[k] = x * axis[k] + sx * (psi.cos() * e1[k] + psi.sin() * e2[k]);
            }
        },
        0 => return Err("heat bath needs N >= 1".into()),
        _ => {
            // von Mises-Fisher distribution on S^(N-1) with concentration `a` (Wood, 1994):
            // sample the cosine `w` to the mean direction by rejection from a Beta envelope,
            // then a uniform direction orthogonal to it.
            let nm1 = (N - 1) as f64;
            let w = if a < 1e-10 {
                2.0 * rng.beta(0.5 * nm1, 0.5 * nm1) - 1.0
            } else {
                let b = nm1 / (2.0 * a + (4.0 * a * a + nm1 * nm1).sqrt());
                let x0 = (1.0 - b) / (1.0 + b);
                let c = a * x0 + nm1 * (1.0 - x0 * x0).ln();
                loop {
                    let z = rng.beta(0.5 * nm1, 0.5 * nm1);
                    let u = rng.uniform();
                    let wc = (1.0 - (1.0 + b) * z) / (1.0 - (1.0 - b) * z);
                    if a * wc + nm1 * (1.0 - x0 * wc).ln() - c >= u.ln() {
                        break wc;
                    }
                }
            };
            let mut mu = [0.0; N];
            if fm > 0.0 {
                for k in 0..N {
                    mu[k] = f[k] / fm;
                }
            } else {
                mu[0] = 1.0;
            }
            // uniform unit vector orthogonal to mu
            let perp = loop {
                let mut g = [0.0; N];
                for c in &mut g {
                    *c = rng.normal();
                }
                let proj = dot(&g, &mu);
                for k in 0..N {
                    g[k] -= proj * mu[k];
                }
                let n = dot(&g, &g).sqrt();
                if n > 1e-12 {
                    for c in &mut g {
                        *c /= n;
                    }
                    break g;
                }
            };
            let sw = (1.0 - w * w).max(0.0).sqrt();
            for k in 0..N {
                out[k] = w * mu[k] + sw * perp[k];
            }
        },
    }
    Ok(out)
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

    /// Local field `h e_0 + J sum_nn s_j` acting on site `i`.
    #[must_use]
    pub fn local_field_public(
        &self,
        i: usize,
    ) -> [f64; N] {
        self.local_field(i)
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

    /// One heat-bath sweep: every site is redrawn from its exact conditional distribution
    /// `P(s) ~ exp(beta s . f)` (`f` = local field). Ising: logistic; XY: von Mises
    /// (Best-Fisher rejection); Heisenberg: exact inverse CDF of `exp(a cos theta)`.
    ///
    /// # Errors
    /// Returns an error for `N = 0`.
    pub fn heat_bath_sweep(&mut self) -> Result<(), String> {
        for i in 0..self.lattice.sites() {
            let f = self.local_field(i);
            self.spins[i] = heat_bath_direction::<N>(&mut self.rng, &f, self.beta)?;
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
        if self.j <= 0.0 {
            return Err("cluster updates require a ferromagnetic coupling".into());
        }
        Ok(())
    }

    /// Grows and flips one Wolff cluster; returns its size.
    ///
    /// A non-zero field `h` is handled by the **ghost-spin** construction: the field couples every
    /// site to an extra ghost spin `sign(h) e_0` with strength `|h|`, the cluster may absorb the
    /// ghost, and when it does the cluster's complement is flipped instead (the global symmetry of the
    /// extended system), which keeps detailed balance.
    ///
    /// # Errors
    /// Fails for `J <= 0`.
    pub fn wolff_step(&mut self) -> Result<usize, String> {
        self.cluster_ready()?;
        let r = random_unit::<N>(&mut self.rng);
        let n = self.lattice.sites();
        let seed = self.rng.below(n as u64) as usize;
        let mut in_cluster = vec![false; n];
        let mut stack = vec![seed];
        in_cluster[seed] = true;
        let mut size = 0;
        let hh = self.field.abs();
        let g0 = r[0] * self.field.signum();
        let mut ghost_in = false;
        while let Some(i) = stack.pop() {
            let pi = dot(&self.spins[i], &r);
            // ghost bond: the cluster absorbs the ghost, which then bonds to every other site
            if !ghost_in && hh > 0.0 && pi * g0 > 0.0 && self.rng.uniform() < 1.0 - (-2.0 * self.beta * hh * pi * g0).exp() {
                ghost_in = true;
                for j in 0..n {
                    if in_cluster[j] || j == i {
                        continue;
                    }
                    let pj = dot(&self.spins[j], &r);
                    if pj * g0 > 0.0 && self.rng.uniform() < 1.0 - (-2.0 * self.beta * hh * pj * g0).exp() {
                        in_cluster[j] = true;
                        stack.push(j);
                    }
                }
            }
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
        if ghost_in {
            // the cluster contained the ghost: reflecting everything undoes the cluster flip and
            // reflects its complement
            for sp in &mut self.spins {
                let p = dot(sp, &r);
                for k in 0..N {
                    sp[k] -= 2.0 * p * r[k];
                }
            }
        }
        Ok(size)
    }

    /// One Swendsen-Wang sweep with a random embedding direction; returns the number of clusters.
    ///
    /// A non-zero field is included through the ghost spin `sign(h) e_0`: the cluster attached to
    /// the ghost is never flipped.
    ///
    /// # Errors
    /// Fails for `J <= 0`.
    pub fn swendsen_wang_step(&mut self) -> Result<usize, String> {
        self.cluster_ready()?;
        let r = random_unit::<N>(&mut self.rng);
        let n = self.lattice.sites();
        let proj: Vec<f64> = self.spins.iter().map(|s| dot(s, &r)).collect();
        let mut uf = UnionFind::new(n + 1);
        let (beta, jj) = (self.beta, self.j);
        let hh = self.field.abs();
        let g0 = r[0] * self.field.signum();
        if hh > 0.0 {
            for i in 0..n {
                let prod = proj[i] * g0;
                if prod > 0.0 && self.rng.uniform() < 1.0 - (-2.0 * beta * hh * prod).exp() {
                    uf.union(i, n);
                }
            }
        }
        let mut links = Vec::new();
        self.lattice.for_each_bond(|i, j| links.push((i, j)));
        for (i, j) in links {
            let prod = proj[i] * proj[j];
            if prod > 0.0 && self.rng.uniform() < 1.0 - (-2.0 * beta * jj * prod).exp() {
                uf.union(i, j);
            }
        }
        let mut flip = vec![None; n + 1];
        let ghost_root = uf.find(n);
        flip[ghost_root] = Some(false);
        let mut clusters = 0;
        for i in 0..n {
            let root = uf.find(i);
            let f = if let Some(f) = flip[root] { f } else {
                let f = self.rng.uniform() < 0.5;
                flip[root] = Some(f);
                clusters += 1;
                f
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
        f64::from(acc) / n as f64
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

/// Replica exchange for the q-state Potts model.
#[derive(Clone, Debug)]
pub struct PottsTempering<const D: usize> {
    /// One replica per temperature slot.
    pub replicas: Vec<Potts<D>>,
    /// Inverse temperatures.
    pub betas: Vec<f64>,
    /// Attempted / accepted swaps per neighbouring pair.
    pub swap_stats: Vec<(usize, usize)>,
    rng: Rng,
}

impl<const D: usize> PottsTempering<D> {
    /// Creates replicas at the given inverse temperatures.
    ///
    /// # Errors
    /// Propagates lattice errors; requires at least two temperatures and `q >= 2`.
    pub fn new(
        dims: [usize; D],
        periodic: bool,
        q: usize,
        betas: &[f64],
        seed: u64,
    ) -> Result<Self, String> {
        if betas.len() < 2 {
            return Err("need at least two temperatures".into());
        }
        let mut replicas = Vec::new();
        for (k, &b) in betas.iter().enumerate() {
            replicas.push(Potts::new(dims, periodic, q, b, seed.wrapping_add(2000 + k as u64))?);
        }
        Ok(Self { replicas, betas: betas.to_vec(), swap_stats: vec![(0, 0); betas.len() - 1], rng: Rng::new(seed ^ 0x5EED) })
    }

    /// `rounds` rounds of (`update` on every replica, then alternating even/odd exchanges
    /// accepted with probability `min(1, exp((b_k - b_k+1)(E_k - E_k+1)))`).
    pub fn run(
        &mut self,
        rounds: usize,
        update: Update,
        mut measure: impl FnMut(usize, &[Potts<D>]),
    ) {
        for round in 0..rounds {
            for r in &mut self.replicas {
                r.update(update);
            }
            let mut k = round % 2;
            while k + 1 < self.replicas.len() {
                let (e1, e2) = (self.replicas[k].energy(), self.replicas[k + 1].energy());
                let d = (self.betas[k] - self.betas[k + 1]) * (e1 - e2);
                self.swap_stats[k].0 += 1;
                if d >= 0.0 || self.rng.uniform() < d.exp() {
                    let (a, b) = self.replicas.split_at_mut(k + 1);
                    std::mem::swap(&mut a[k].spins, &mut b[0].spins);
                    self.swap_stats[k].1 += 1;
                }
                k += 2;
            }
            measure(round, &self.replicas);
        }
    }
}

/// Wang-Landau density of states of an O(N) model with **continuous** spins on a binned
/// energy axis `[e_min, e_max)` with `nbins` bins (moves leaving the window are rejected).
///
/// Single-spin moves of width `delta`; `ln f` starts at 1 and is halved whenever the histogram
/// of visited bins satisfies `min H > flatness * mean H` (checked every `check` attempted
/// moves) until `ln_f_final`. Returns `(bin centre, ln g)` for the visited bins, normalised
/// to `sum g = 1`.
///
/// # Errors
/// Fails for invalid bins or lattice errors.
pub fn wang_landau_on<const D: usize, const N: usize>(
    dims: [usize; D],
    periodic: bool,
    (e_min, e_max, nbins): (f64, f64, usize),
    delta: f64,
    ln_f_final: f64,
    flatness: f64,
    seed: u64,
) -> Result<Vec<(f64, f64)>, String> {
    if nbins < 2 || e_max <= e_min {
        return Err("invalid energy window".into());
    }
    let mut m = OnModel::<D, N>::new(dims, periodic, 1.0, seed)?;
    let n = m.lattice.sites();
    let width = (e_max - e_min) / nbins as f64;
    let bin = |e: f64| -> Option<usize> {
        if e < e_min || e >= e_max {
            None
        } else {
            Some((((e - e_min) / width) as usize).min(nbins - 1))
        }
    };
    let mut e = m.energy();
    let mut cur = bin(e).ok_or_else(|| "initial energy outside the window; widen it or start cold".to_string())?;
    let mut ln_g = vec![0.0f64; nbins];
    let mut visited = vec![false; nbins];
    let mut ln_f = 1.0;
    let check = 20_000 * n;
    let mut guard = 0usize;
    while ln_f > ln_f_final {
        let mut hist = vec![0usize; nbins];
        loop {
            for _ in 0..check {
                let i = m.rng.below(n as u64) as usize;
                let f = m.local_field_public(i);
                let old = m.spins[i];
                let mut new = old;
                if N == 1 {
                    new[0] = -old[0];
                } else {
                    for k in 0..N {
                        new[k] = old[k] + delta * m.rng.normal();
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
                if let Some(nb) = bin(e + de) {
                    let d = ln_g[cur] - ln_g[nb];
                    if d >= 0.0 || m.rng.uniform() < d.exp() {
                        m.spins[i] = new;
                        e += de;
                        cur = nb;
                    }
                }
                ln_g[cur] += ln_f;
                hist[cur] += 1;
                visited[cur] = true;
            }
            // resynchronise the running energy against round-off
            e = m.energy();
            cur = bin(e).unwrap_or(cur);
            guard += 1;
            let vis: Vec<usize> = (0..nbins).filter(|&k| visited[k]).map(|k| hist[k]).collect();
            let mean = vis.iter().sum::<usize>() as f64 / vis.len() as f64;
            let mn = vis.iter().copied().min().unwrap_or(0) as f64;
            if mn > flatness * mean || guard > 2000 {
                break;
            }
        }
        ln_f /= 2.0;
    }
    let mx = (0..nbins).filter(|&k| visited[k]).map(|k| ln_g[k]).fold(f64::NEG_INFINITY, f64::max);
    let total: f64 = (0..nbins).filter(|&k| visited[k]).map(|k| (ln_g[k] - mx).exp()).sum();
    let shift = -(mx + total.ln());
    Ok((0..nbins).filter(|&k| visited[k]).map(|k| (e_min + (k as f64 + 0.5) * width, ln_g[k] + shift)).collect())
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
    let idx = |e: i64| i64::midpoint(e, nb) as usize;
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
    let (mean_e2, me_sq, beta_sq) = (e2 / z, me * me, beta * beta);
    let var = mean_e2 - me_sq;
    (m + z.ln(), me, beta_sq * var)
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
        cur = cur.as_chunks::<2>().0.iter().map(|c| f64::midpoint(c[0], c[1])).collect();
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
    2.0 * j / std::f64::consts::SQRT_2.ln_1p()
}

fn ellip_k(k: f64) -> f64 {
    let (mut a, mut b) = (1.0, (1.0 - k * k).sqrt());
    for _ in 0..40 {
        let (an, bn) = (f64::midpoint(a, b), (a * b).sqrt());
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

// ---------------------------------------------------------------------------
// Replica-exchange Wang-Landau with multiple walkers per window
// ---------------------------------------------------------------------------

/// Parameters of [`replica_exchange_wang_landau`].
#[derive(Clone, Copy, Debug)]
pub struct RewlOptions {
    /// Energy window `(e_min, e_max)` covered by all energy windows together.
    pub range: (f64, f64),
    /// Number of energy bins over `range`.
    pub nbins: usize,
    /// Number of overlapping energy windows.
    pub windows: usize,
    /// Fraction of a window shared with its neighbour, in `[0, 1)`.
    pub overlap: f64,
    /// Walkers per window; all walkers of a window share one `ln g` and histogram.
    pub walkers: usize,
    /// Single-spin proposal width for continuous spins (ignored for Ising).
    pub delta: f64,
    /// Final modification factor `ln f`.
    pub ln_f_final: f64,
    /// Flatness criterion `min H > flatness mean H` per window.
    pub flatness: f64,
    /// Sweeps of every walker between two exchange attempts and flatness checks.
    pub sweeps_per_exchange: usize,
    /// Random seed.
    pub seed: u64,
}

/// One walker of the replica-exchange Wang-Landau simulation.
struct WlWalker<const D: usize, const N: usize> {
    model: OnModel<D, N>,
    energy: f64,
    bin: usize,
}

/// Greedy steering of a configuration into the energy interval `[lo, hi)`.
fn steer_into_window<const D: usize, const N: usize>(
    m: &mut OnModel<D, N>,
    lo: f64,
    hi: f64,
    delta: f64,
) -> Result<f64, String> {
    let n = m.lattice.sites();
    let target = f64::midpoint(lo, hi);
    let mut e = m.energy();
    for _ in 0..4000 * n.max(16) {
        if e >= lo && e < hi {
            return Ok(e);
        }
        let i = m.rng.below(n as u64) as usize;
        let f = m.local_field_public(i);
        let old = m.spins[i];
        let mut new = old;
        if N == 1 {
            new[0] = -old[0];
        } else {
            for k in 0..N {
                new[k] = old[k] + delta * m.rng.normal();
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
        if (e + de - target).abs() <= (e - target).abs() || m.rng.uniform() < 0.01 {
            m.spins[i] = new;
            e += de;
        }
    }
    let e = m.energy();
    if e >= lo && e < hi {
        Ok(e)
    } else {
        Err("could not steer a walker into its energy window".into())
    }
}

/// Replica-exchange Wang-Landau with several walkers per energy window.
///
/// After Vogel, Li, Wuest and Landau (2013): the energy range is split into overlapping windows, each window runs its own
/// Wang-Landau iteration with `walkers` walkers sharing `ln g`, neighbouring windows exchange
/// configurations with probability `min(1, g_i(E_a) g_j(E_b) / (g_i(E_b) g_j(E_a)))`, and the
/// per-window densities of states are joined at the overlap bin where their logarithmic
/// derivatives agree best. Works for the Ising model (`N = 1`) and continuous O(N) spins.
///
/// Returns `(bin centre, ln g)` for the visited bins, normalised to `sum g = 1`.
///
/// # Errors
/// Fails for invalid options, lattice errors, or if a walker cannot be placed in its window.
pub fn replica_exchange_wang_landau<const D: usize, const N: usize>(
    dims: [usize; D],
    periodic: bool,
    opt: &RewlOptions,
) -> Result<Vec<(f64, f64)>, String> {
    let (e_min, e_max) = opt.range;
    let nbins = opt.nbins;
    let nwin = opt.windows.max(1);
    if nbins < 2 || e_max <= e_min || !(0.0..1.0).contains(&opt.overlap) || opt.walkers == 0 {
        return Err("invalid replica-exchange Wang-Landau options".into());
    }
    let width = (e_max - e_min) / nbins as f64;
    // window layout in bins: [lo_w, hi_w)
    let span = 1.0 + (nwin - 1) as f64 * (1.0 - opt.overlap);
    let nw = ((nbins as f64 / span).ceil() as usize).clamp(2, nbins);
    let stride = if nwin > 1 { (((nbins - nw) as f64) / (nwin - 1) as f64).ceil().max(1.0) as usize } else { 0 };
    let ranges: Vec<(usize, usize)> = (0..nwin)
        .map(|w| {
            let lo = (w * stride).min(nbins - nw);
            let hi = if w + 1 == nwin { nbins } else { (lo + nw).min(nbins) };
            (lo, hi)
        })
        .collect();
    for w in 1..nwin {
        if ranges[w].0 >= ranges[w - 1].1 {
            return Err("energy windows do not overlap; increase `overlap` or `nbins`".into());
        }
    }
    let bin_of = |e: f64| -> Option<usize> {
        if e < e_min || e >= e_max {
            None
        } else {
            Some((((e - e_min) / width) as usize).min(nbins - 1))
        }
    };
    let mut walkers: Vec<Vec<WlWalker<D, N>>> = Vec::new();
    for (w, &(lo, hi)) in ranges.iter().enumerate() {
        let mut ws = Vec::new();
        for k in 0..opt.walkers {
            let seed = opt.seed.wrapping_add(7919 * (w * opt.walkers + k) as u64 + 1);
            let mut model = OnModel::<D, N>::new(dims, periodic, 1.0, seed)?;
            let elo = e_min + lo as f64 * width;
            let ehi = e_min + hi as f64 * width;
            let energy = steer_into_window(&mut model, elo, ehi, opt.delta)?;
            let bin = bin_of(energy).ok_or("walker left the energy range")?;
            ws.push(WlWalker { model, energy, bin });
        }
        walkers.push(ws);
    }
    let n = walkers[0][0].model.lattice.sites();
    let mut ln_g = vec![vec![0.0f64; nbins]; nwin];
    let mut visited = vec![vec![false; nbins]; nwin];
    let mut hist = vec![vec![0usize; nbins]; nwin];
    let mut ln_f = 1.0;
    let mut xrng = Rng::new(opt.seed ^ 0x9E37_79B9);
    let steps = opt.sweeps_per_exchange.max(1) * n;
    let mut round = 0usize;
    let mut guard = 0usize;
    while ln_f > opt.ln_f_final {
        for h in &mut hist {
            h.fill(0);
        }
        loop {
            for (w, &(lo, hi)) in ranges.iter().enumerate() {
                for wk in &mut walkers[w] {
                    for _ in 0..steps {
                        let i = wk.model.rng.below(n as u64) as usize;
                        let f = wk.model.local_field_public(i);
                        let old = wk.model.spins[i];
                        let mut new = old;
                        if N == 1 {
                            new[0] = -old[0];
                        } else {
                            for k in 0..N {
                                new[k] = old[k] + opt.delta * wk.model.rng.normal();
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
                        if let Some(nb) = bin_of(wk.energy + de)
                            && nb >= lo
                            && nb < hi
                        {
                            let d = ln_g[w][wk.bin] - ln_g[w][nb];
                            if d >= 0.0 || wk.model.rng.uniform() < d.exp() {
                                wk.model.spins[i] = new;
                                wk.energy += de;
                                wk.bin = nb;
                            }
                        }
                        ln_g[w][wk.bin] += ln_f;
                        hist[w][wk.bin] += 1;
                        visited[w][wk.bin] = true;
                    }
                    wk.energy = wk.model.energy();
                    wk.bin = bin_of(wk.energy).unwrap_or(wk.bin).clamp(lo, hi - 1);
                }
            }
            // replica exchange between neighbouring windows (alternating parity)
            let mut w = round % 2;
            while w + 1 < nwin {
                let a = xrng.below(opt.walkers as u64) as usize;
                let b = xrng.below(opt.walkers as u64) as usize;
                let (ba, bb) = (walkers[w][a].bin, walkers[w + 1][b].bin);
                let (lo1, hi1) = ranges[w + 1];
                let (lo0, hi0) = ranges[w];
                if ba >= lo1 && ba < hi1 && bb >= lo0 && bb < hi0 {
                    let d = ln_g[w][ba] - ln_g[w][bb] + ln_g[w + 1][bb] - ln_g[w + 1][ba];
                    if d >= 0.0 || xrng.uniform() < d.exp() {
                        let (left, right) = walkers.split_at_mut(w + 1);
                        std::mem::swap(&mut left[w][a].model.spins, &mut right[0][b].model.spins);
                        std::mem::swap(&mut left[w][a].energy, &mut right[0][b].energy);
                        std::mem::swap(&mut left[w][a].bin, &mut right[0][b].bin);
                    }
                }
                w += 2;
            }
            round += 1;
            guard += 1;
            let flat = ranges.iter().enumerate().all(|(w, &(lo, hi))| {
                let vals: Vec<usize> = (lo..hi).filter(|&k| visited[w][k]).map(|k| hist[w][k]).collect();
                if vals.is_empty() {
                    return false;
                }
                let mean = vals.iter().sum::<usize>() as f64 / vals.len() as f64;
                let mn = vals.iter().copied().min().unwrap_or(0) as f64;
                mn > opt.flatness * mean
            });
            if flat || guard > 200_000 {
                break;
            }
        }
        ln_f /= 2.0;
    }
    // join the windows: shift each window onto the merged left part by the mean offset over the
    // interior overlap bins (window edge bins carry the Wang-Landau boundary artefact), and
    // switch from the left to the right window in the middle of the overlap
    let mut merged = vec![f64::NAN; nbins];
    let mut known = vec![false; nbins];
    let (lo0, hi0) = ranges[0];
    for k in lo0..hi0 {
        if visited[0][k] {
            merged[k] = ln_g[0][k];
            known[k] = true;
        }
    }
    for w in 1..nwin {
        let (lo, hi) = ranges[w];
        let prev_hi = ranges[w - 1].1;
        let mut common: Vec<usize> = (lo..prev_hi).filter(|&k| visited[w][k] && known[k]).collect();
        if common.is_empty() {
            return Err("adjacent windows share no visited bin".into());
        }
        let all = common.clone();
        if common.len() >= 4 {
            common.remove(0);
            common.pop();
        }
        let shift = common.iter().map(|&k| merged[k] - ln_g[w][k]).sum::<f64>() / common.len() as f64;
        let mid = all[all.len() / 2];
        for k in mid..hi {
            if visited[w][k] {
                merged[k] = ln_g[w][k] + shift;
                known[k] = true;
            }
        }
    }
    let mx = (0..nbins).filter(|&k| known[k]).map(|k| merged[k]).fold(f64::NEG_INFINITY, f64::max);
    let total: f64 = (0..nbins).filter(|&k| known[k]).map(|k| (merged[k] - mx).exp()).sum();
    let shift = -(mx + total.ln());
    Ok((0..nbins).filter(|&k| known[k]).map(|k| (e_min + (k as f64 + 0.5) * width, merged[k] + shift)).collect())
}
