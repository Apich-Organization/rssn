//! # Hermite N-body integrators with individual block time steps and post-Newtonian gravity
//!
//! Direct-summation gravity in 3-D (`G` explicit, optional Plummer softening) integrated with
//! the Hermite predictor-corrector family used by collisional stellar dynamics codes:
//!
//! | Feature | Notes |
//! |---|---|
//! | 4th-order Hermite (acceleration and jerk, Makino-Aarseth) | corrector with snap and crackle from the two end points, optional `P(EC)^n` iteration |
//! | 6th-order Hermite (Nitadori-Makino: acceleration, jerk and snap, analytic pairwise snap) | predictor to crackle order, 6th-order corrector |
//! | Individual **block time steps** | power-of-two hierarchy, Aarseth criterion `dt = sqrt(eta (|a||s| + |j|^2)/(|j||c| + |s|^2))`, only commensurable doubling, only active particles are updated |
//! | Shared time step | for convergence studies ([`StepPolicy::Shared`]) |
//! | 1PN Einstein-Infeld-Hoffmann (EIH) corrections ([`HermiteSystem::pn_c`]) | acceleration and **jerk** of the PN terms by forward-mode (dual-number) differentiation; 4th-order Hermite only; [`eih_energy`] is the conserved 1PN energy |
//! | Hierarchical hybrid ([`Binary`]) | tight binaries advanced by Kustaanheimo-Stiefel regularisation inside the block-step loop, see [`HermiteSystem::add_binary`] |
//!
//! Limitations: direct `O(N^2)` summation (no neighbour scheme), no GPU/vectorised kernels; PN
//! terms are 1PN only (no spin, no radiation reaction) and need `eps = 0`; the binary hybrid
//! supports the 4th-order scheme and treats other binaries as point masses in the
//! perturbation of a binary.

use std::ops::{Add, Div, Mul, Neg, Sub};

type V3 = [f64; 3];

/// Acceleration, jerk, snap, Newtonian acceleration and Newtonian jerk of one particle.
type ForceSet = (V3, V3, V3, V3, V3);

/// Order of the Hermite scheme.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HermiteOrder {
    /// Fourth order (acceleration and jerk).
    Fourth,
    /// Sixth order (acceleration, jerk and snap).
    Sixth,
}

/// Time-step policy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum StepPolicy {
    /// One shared constant step for all particles.
    Shared(f64),
    /// Individual block steps, powers of two of `dt_max` (a power of two itself),
    /// accuracy parameter `eta`, smallest allowed step `dt_min`.
    Block {
        /// Accuracy parameter of the Aarseth criterion.
        eta: f64,
        /// Largest step (rounded down to a power of two).
        dt_max: f64,
        /// Smallest step.
        dt_min: f64,
    },
}

/// Configuration of a [`HermiteSystem`].
#[derive(Clone, Copy, Debug)]
pub struct HermiteConfig {
    /// Gravitational constant.
    pub g: f64,
    /// Plummer softening length.
    pub eps: f64,
    /// Scheme order.
    pub order: HermiteOrder,
    /// Time-step policy.
    pub policy: StepPolicy,
    /// Initial time.
    pub t0: f64,
}

/// A tight binary advanced by KS regularisation inside the Hermite loop.
#[derive(Clone, Debug)]
pub struct Binary {
    /// Indices of the two components in the *component* arrays (not in the particle arrays).
    pub members: [usize; 2],
    /// Index of the centre-of-mass particle in the particle arrays.
    pub cm: usize,
    /// Index of the (inactive) partner slot that was removed from the dynamics.
    pub partner: usize,
    /// Relative position `x_2 - x_1` at the last time of the centre of mass.
    pub rel: V3,
    /// Relative velocity `v_2 - v_1`.
    pub rel_vel: V3,
    /// KS accuracy parameter (fraction of the harmonic period per RK4 step).
    pub ks_eta: f64,
}

/// Hermite N-body system.
#[derive(Clone, Debug)]
pub struct HermiteSystem {
    /// Scheme order.
    pub order: HermiteOrder,
    /// Time-step policy.
    pub policy: StepPolicy,
    /// Gravitational constant.
    pub g: f64,
    /// Plummer softening length.
    pub eps: f64,
    /// Speed of light for 1PN corrections (`None` = Newtonian).
    pub pn_c: Option<f64>,
    /// Number of `P(EC)` iterations of the 4th-order corrector (1 or 2).
    pub iterations: usize,
    /// Particle masses (centre-of-mass particles carry the total binary mass).
    pub mass: Vec<f64>,
    /// Positions at the particles' own last update times.
    pub pos: Vec<V3>,
    /// Velocities at the particles' own last update times.
    pub vel: Vec<V3>,
    /// Accelerations at the last update times.
    pub acc: Vec<V3>,
    /// Jerks at the last update times.
    pub jerk: Vec<V3>,
    /// Snaps (second derivative of the acceleration).
    pub snap: Vec<V3>,
    /// Crackles (third derivative of the acceleration).
    pub crackle: Vec<V3>,
    /// Last update time of every particle.
    pub t_last: Vec<f64>,
    /// Current individual time step of every particle.
    pub dt: Vec<f64>,
    /// Current system time (time of the latest block).
    pub time: f64,
    /// Number of particle-particle force evaluations so far.
    pub interactions: u64,
    /// Number of particle updates so far.
    pub updates: u64,
    /// Tight binaries (hierarchical hybrid).
    pub binaries: Vec<Binary>,
    /// Masses of the binary components, indexed by [`Binary::members`].
    pub comp_mass: Vec<f64>,
    accn: Vec<V3>,
    jerkn: Vec<V3>,
}

fn sub(
    a: V3,
    b: V3,
) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn add(
    a: V3,
    b: V3,
) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn scale(
    a: V3,
    s: f64,
) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn dot(
    a: V3,
    b: V3,
) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn norm(a: V3) -> f64 {
    dot(a, a).sqrt()
}

// ---------------------------------------------------------------------------
// Forward-mode dual numbers for the PN jerk
// ---------------------------------------------------------------------------

/// Scalar field used by the generic EIH acceleration (plain reals or dual numbers).
trait Scalar:
    Copy + Add<Output = Self> + Sub<Output = Self> + Mul<Output = Self> + Div<Output = Self> + Neg<Output = Self>
{
    fn from_f64(x: f64) -> Self;
    fn sqrt(self) -> Self;
}

impl Scalar for f64 {
    fn from_f64(x: f64) -> Self {
        x
    }

    fn sqrt(self) -> Self {
        self.sqrt()
    }
}

/// First-order dual number `v + d eps`.
#[derive(Clone, Copy, Debug)]
struct Dual {
    v: f64,
    d: f64,
}

impl Add for Dual {
    type Output = Self;

    fn add(
        self,
        o: Self,
    ) -> Self {
        Self { v: self.v + o.v, d: self.d + o.d }
    }
}

impl Sub for Dual {
    type Output = Self;

    fn sub(
        self,
        o: Self,
    ) -> Self {
        Self { v: self.v - o.v, d: self.d - o.d }
    }
}

impl Mul for Dual {
    type Output = Self;

    fn mul(
        self,
        o: Self,
    ) -> Self {
        Self { v: self.v * o.v, d: self.v * o.d + self.d * o.v }
    }
}

impl Div for Dual {
    type Output = Self;

    fn div(
        self,
        o: Self,
    ) -> Self {
        Self { v: self.v / o.v, d: (self.d * o.v - self.v * o.d) / (o.v * o.v) }
    }
}

impl Neg for Dual {
    type Output = Self;

    fn neg(self) -> Self {
        Self { v: -self.v, d: -self.d }
    }
}

impl Scalar for Dual {
    fn from_f64(x: f64) -> Self {
        Self { v: x, d: 0.0 }
    }

    fn sqrt(self) -> Self {
        let s = self.v.sqrt();
        Self { v: s, d: 0.5 * self.d / s }
    }
}

fn dot_s<T: Scalar>(
    a: &[T; 3],
    b: &[T; 3],
) -> T {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// 1PN Einstein-Infeld-Hoffmann acceleration (Newtonian part included) of particle `i` (harmonic
/// coordinates, Will 1993 / Weinberg): `a_i^PN`, with the Newtonian accelerations `an` of the
/// sources entering through the `a_j` terms.
fn eih_correction<T: Scalar>(
    i: usize,
    x: &[[T; 3]],
    v: &[[T; 3]],
    an: &[[T; 3]],
    mass: &[f64],
    g: f64,
    c: f64,
) -> [T; 3] {
    let n = x.len();
    let c2inv = T::from_f64(1.0 / (c * c));
    let dist = |a: usize, b: usize| {
        let d = [x[a][0] - x[b][0], x[a][1] - x[b][1], x[a][2] - x[b][2]];
        (d, dot_s(&d, &d).sqrt())
    };
    let potential_at = |a: usize| {
        let mut s = T::from_f64(0.0);
        for k in 0..n {
            if k != a {
                let (_, r) = dist(a, k);
                s = s + T::from_f64(g * mass[k]) / r;
            }
        }
        s
    };
    let s_i = potential_at(i);
    let vi2 = dot_s(&v[i], &v[i]);
    let mut acc = [T::from_f64(0.0); 3];
    for j in 0..n {
        if j == i {
            continue;
        }
        let (rij, r) = dist(i, j);
        let gm = T::from_f64(g * mass[j]);
        let nij = [rij[0] / r, rij[1] / r, rij[2] / r];
        let s_j = potential_at(j);
        let vj2 = dot_s(&v[j], &v[j]);
        let nvj = dot_s(&nij, &v[j]);
        let rij_an = dot_s(&rij, &an[j]);
        // (r_j - r_i) . a_j = - rij . a_j
        let bracket = T::from_f64(1.0)
            - c2inv
                * (T::from_f64(4.0) * s_i + s_j - vi2 - T::from_f64(2.0) * vj2
                    + T::from_f64(4.0) * dot_s(&v[i], &v[j])
                    + T::from_f64(1.5) * nvj * nvj
                    + T::from_f64(0.5) * rij_an);
        let r2 = r * r;
        let radial = -(gm / r2) * bracket;
        let four_vi_minus_3vj = [
            T::from_f64(4.0) * v[i][0] - T::from_f64(3.0) * v[j][0],
            T::from_f64(4.0) * v[i][1] - T::from_f64(3.0) * v[j][1],
            T::from_f64(4.0) * v[i][2] - T::from_f64(3.0) * v[j][2],
        ];
        let proj = dot_s(&nij, &four_vi_minus_3vj);
        let dv = [v[i][0] - v[j][0], v[i][1] - v[j][1], v[i][2] - v[j][2]];
        for k in 0..3 {
            acc[k] = acc[k]
                + radial * nij[k]
                + c2inv * (gm / r2) * proj * dv[k]
                + T::from_f64(3.5) * c2inv * (gm / r) * an[j][k];
        }
    }
    acc
}

/// Full 1PN (EIH) accelerations (Newtonian part included) of all particles, given the Newtonian
/// accelerations `acc_newton` of the sources.
#[must_use]
pub fn eih_accelerations(
    mass: &[f64],
    pos: &[V3],
    vel: &[V3],
    acc_newton: &[V3],
    g: f64,
    c: f64,
) -> Vec<V3> {
    (0..mass.len()).map(|i| eih_correction::<f64>(i, pos, vel, acc_newton, mass, g, c)).collect()
}

/// Conserved 1PN energy of the EIH N-body Lagrangian (Landau-Lifshitz section 106).
#[must_use]
pub fn eih_energy(
    mass: &[f64],
    pos: &[V3],
    vel: &[V3],
    g: f64,
    c: f64,
) -> f64 {
    let n = mass.len();
    let c2 = c * c;
    let mut e = 0.0;
    for a in 0..n {
        let v2 = dot(vel[a], vel[a]);
        e += 0.5 * mass[a] * v2 + 3.0 / (8.0 * c2) * mass[a] * v2 * v2;
        for b in 0..n {
            if a == b {
                continue;
            }
            let d = sub(pos[a], pos[b]);
            let r = norm(d);
            let nn = scale(d, 1.0 / r);
            let (va2, vb2) = (v2, dot(vel[b], vel[b]));
            let vel_term = 1.5 * (va2 + vb2) - 3.5 * dot(vel[a], vel[b]) - 0.5 * dot(nn, vel[a]) * dot(nn, vel[b]);
            e -= 0.5 * g * mass[a] * mass[b] / r * (1.0 - vel_term / c2);
            for cc in 0..n {
                if cc != a {
                    let rc = norm(sub(pos[a], pos[cc]));
                    e += 0.5 * g * g * mass[a] * mass[b] * mass[cc] / (r * rc * c2);
                }
            }
        }
    }
    e
}

// ---------------------------------------------------------------------------
// Pairwise Newtonian kernels
// ---------------------------------------------------------------------------

/// Acceleration and jerk of `i` due to a source: `(a, j)`.
fn pair_aj(
    dx: V3,
    dv: V3,
    gm: f64,
    eps2: f64,
) -> (V3, V3) {
    let r2 = dot(dx, dx) + eps2;
    let rinv = 1.0 / r2.sqrt();
    let r3inv = rinv / r2;
    let rv = dot(dx, dv) / r2;
    let f = gm * r3inv;
    (scale(dx, f), scale(sub(dv, scale(dx, 3.0 * rv)), f))
}

/// Snap of `i` due to a source, given the relative acceleration `da = a_source - a_i`.
fn pair_snap(
    dx: V3,
    dv: V3,
    da: V3,
    gm: f64,
    eps2: f64,
) -> V3 {
    let r2 = dot(dx, dx) + eps2;
    let rinv = 1.0 / r2.sqrt();
    let f = gm * rinv / r2;
    let alpha = dot(dx, dv) / r2;
    let b = (dot(dv, dv) + dot(dx, da)) / r2;
    let mut s = da;
    for k in 0..3 {
        s[k] += -6.0 * alpha * dv[k] + (-3.0 * b + 15.0 * alpha * alpha) * dx[k];
    }
    scale(s, f)
}

/// Largest power of two not larger than `x` (`x > 0`).
fn pow2_floor(x: f64) -> f64 {
    2.0_f64.powf(x.log2().floor())
}

/// Predicted state of a particle at time `t`.
#[derive(Clone, Copy, Debug)]
struct Predicted {
    x: V3,
    v: V3,
    a: V3,
    j: V3,
    s: V3,
}

impl HermiteSystem {
    /// Creates a system (all particles at time `t0`) and evaluates the initial derivatives.
    ///
    /// # Errors
    /// Fails for mismatched array lengths, fewer than two particles, a non-positive step
    /// parameter, or PN combined with softening or the sixth-order scheme.
    pub fn new(
        mass: Vec<f64>,
        pos: Vec<V3>,
        vel: Vec<V3>,
        cfg: &HermiteConfig,
    ) -> Result<Self, String> {
        let HermiteConfig { g, eps, order, policy, t0 } = *cfg;
        let n = mass.len();
        if n < 2 || pos.len() != n || vel.len() != n {
            return Err("need at least two particles with matching arrays".into());
        }
        match policy {
            StepPolicy::Shared(h) if h <= 0.0 => return Err("shared step must be positive".into()),
            StepPolicy::Block { eta, dt_max, dt_min } if eta <= 0.0 || dt_max <= 0.0 || dt_min <= 0.0 || dt_min > dt_max => {
                return Err("invalid block step parameters".into());
            },
            _ => {},
        }
        let z = vec![[0.0; 3]; n];
        let mut s = Self {
            order,
            policy,
            g,
            eps,
            pn_c: None,
            iterations: 1,
            mass,
            pos,
            vel,
            acc: z.clone(),
            jerk: z.clone(),
            snap: z.clone(),
            crackle: z.clone(),
            t_last: vec![t0; n],
            dt: vec![0.0; n],
            time: t0,
            interactions: 0,
            updates: 0,
            binaries: Vec::new(),
            comp_mass: Vec::new(),
            accn: z.clone(),
            jerkn: z,
        };
        s.initialise()?;
        Ok(s)
    }

    /// Enables the 1PN EIH corrections with speed of light `c` and re-initialises.
    ///
    /// # Errors
    /// Fails with the sixth-order scheme, with softening, or with binaries.
    pub fn enable_post_newtonian(
        &mut self,
        c: f64,
    ) -> Result<(), String> {
        if self.order != HermiteOrder::Fourth || self.eps != 0.0 || !self.binaries.is_empty() || c <= 0.0 {
            return Err("post-Newtonian terms need the 4th-order scheme, eps = 0, no binaries and c > 0".into());
        }
        self.pn_c = Some(c);
        self.initialise()
    }

    fn initialise(&mut self) -> Result<(), String> {
        let n = self.mass.len();
        let t = self.time;
        let pred: Vec<Predicted> = (0..n)
            .map(|i| Predicted { x: self.pos[i], v: self.vel[i], a: [0.0; 3], j: [0.0; 3], s: [0.0; 3] })
            .collect();
        let all: Vec<usize> = (0..n).filter(|&i| !self.is_dead(i)).collect();
        let out = self.forces(&pred, &all, t)?;
        for (k, &i) in all.iter().enumerate() {
            self.acc[i] = out[k].0;
            self.jerk[i] = out[k].1;
            self.snap[i] = out[k].2;
            self.crackle[i] = [0.0; 3];
            self.accn[i] = out[k].3;
            self.jerkn[i] = out[k].4;
            self.t_last[i] = t;
        }
        let live: Vec<usize> = (0..n).filter(|&i| !self.is_dead(i)).collect();
        for i in live {
            self.dt[i] = match self.policy {
                StepPolicy::Shared(h) => h,
                StepPolicy::Block { eta, dt_max, dt_min } => {
                    let (a, j) = (norm(self.acc[i]), norm(self.jerk[i]));
                    let guess = if j > 0.0 { 0.01 * eta.sqrt() * a / j } else { dt_max };
                    Self::quantise(guess.min(dt_max), dt_max, dt_min)
                },
            };
        }
        Ok(())
    }

    fn quantise(
        proposed: f64,
        dt_max: f64,
        dt_min: f64,
    ) -> f64 {
        let mut dt = pow2_floor(dt_max);
        while dt > proposed && dt > dt_min {
            dt *= 0.5;
        }
        dt
    }

    /// Prediction of particle `i` to time `t` (Taylor polynomial of its stored derivatives).
    fn predict(
        &self,
        i: usize,
        t: f64,
    ) -> Predicted {
        let h = t - self.t_last[i];
        let (x0, v0, a0, j0, s0, c0) = (self.pos[i], self.vel[i], self.acc[i], self.jerk[i], self.snap[i], self.crackle[i]);
        let sixth = self.order == HermiteOrder::Sixth;
        let mut p = Predicted { x: [0.0; 3], v: [0.0; 3], a: [0.0; 3], j: [0.0; 3], s: [0.0; 3] };
        for k in 0..3 {
            p.x[k] = x0[k] + h * (v0[k] + h * (a0[k] / 2.0 + h * (j0[k] / 6.0)));
            p.v[k] = v0[k] + h * (a0[k] + h * (j0[k] / 2.0));
            p.a[k] = a0[k] + h * (j0[k] + h * (s0[k] / 2.0 + h * c0[k] / 6.0));
            p.j[k] = j0[k] + h * (s0[k] + h * c0[k] / 2.0);
            p.s[k] = s0[k] + h * c0[k];
            if sixth {
                p.x[k] += h * h * h * h * (s0[k] / 24.0 + h * c0[k] / 120.0);
                p.v[k] += h * h * h * (s0[k] / 6.0 + h * c0[k] / 24.0);
            }
        }
        p
    }

    /// Source list at time `t` from predicted states: (position, velocity, mass) of every point
    /// source; binary centre-of-mass particles are replaced by their two components.
    fn sources(
        &self,
        pred: &[Predicted],
        t: f64,
    ) -> Result<Vec<(V3, V3, f64, usize)>, String> {
        let mut out = Vec::with_capacity(pred.len() + self.binaries.len());
        let mut is_cm = vec![None; pred.len()];
        for (b, bin) in self.binaries.iter().enumerate() {
            is_cm[bin.cm] = Some(b);
        }
        for (i, p) in pred.iter().enumerate() {
            if self.is_dead(i) {
                continue;
            }
            if let Some(b) = is_cm[i] {
                let (x1, v1, x2, v2) = self.binary_components(b, p, t)?;
                let bin = &self.binaries[b];
                out.push((x1, v1, self.comp_mass[bin.members[0]], i));
                out.push((x2, v2, self.comp_mass[bin.members[1]], i));
            } else {
                out.push((p.x, p.v, self.mass[i], i));
            }
        }
        Ok(out)
    }

    /// Accelerations, jerks, snaps (sixth order only), Newtonian acceleration and jerk of the
    /// `active` particles at time `t` from the predicted states.
    fn forces(
        &self,
        pred: &[Predicted],
        active: &[usize],
        t: f64,
    ) -> Result<Vec<ForceSet>, String> {
        let eps2 = self.eps * self.eps;
        let sources = self.sources(pred, t)?;
        let mut is_cm = vec![None; pred.len()];
        for (b, bin) in self.binaries.iter().enumerate() {
            is_cm[bin.cm] = Some(b);
        }
        // target points of every active particle: (x, v, weight)
        let targets = |i: usize| -> Result<Vec<(V3, V3, f64)>, String> {
            if let Some(b) = is_cm[i] {
                let (x1, v1, x2, v2) = self.binary_components(b, &pred[i], t)?;
                let bin = &self.binaries[b];
                let (m1, m2) = (self.comp_mass[bin.members[0]], self.comp_mass[bin.members[1]]);
                Ok(vec![(x1, v1, m1 / (m1 + m2)), (x2, v2, m2 / (m1 + m2))])
            } else {
                Ok(vec![(pred[i].x, pred[i].v, 1.0)])
            }
        };
        let mut res = Vec::with_capacity(active.len());
        let mut an_all: Vec<V3> = pred.iter().map(|p| p.a).collect();
        let mut jn_all: Vec<V3> = pred.iter().map(|p| p.j).collect();
        for &i in active {
            let mut a = [0.0; 3];
            let mut j = [0.0; 3];
            for (xt, vt, w) in targets(i)? {
                for &(xs, vs, ms, owner) in &sources {
                    if owner == i {
                        continue;
                    }
                    let (da, dj) = pair_aj(sub(xs, xt), sub(vs, vt), self.g * ms, eps2);
                    a = add(a, scale(da, w));
                    j = add(j, scale(dj, w));
                }
            }
            an_all[i] = a;
            jn_all[i] = j;
        }
        let sixth = self.order == HermiteOrder::Sixth;
        let mut out = Vec::with_capacity(active.len());
        for &i in active {
            let mut snap = [0.0; 3];
            if sixth {
                let t0 = &targets(i)?[0];
                for &(xs, vs, ms, owner) in &sources {
                    if owner == i {
                        continue;
                    }
                    let da = sub(an_all[owner], an_all[i]);
                    snap = add(snap, pair_snap(sub(xs, t0.0), sub(vs, t0.1), da, self.g * ms, eps2));
                }
            }
            res.push((an_all[i], jn_all[i], snap));
        }
        // post-Newtonian corrections on top
        if let Some(c) = self.pn_c {
            let x: Vec<[Dual; 3]> = pred.iter().map(|p| std::array::from_fn(|k| Dual { v: p.x[k], d: p.v[k] })).collect();
            let v: Vec<[Dual; 3]> = pred
                .iter()
                .enumerate()
                .map(|(i, p)| std::array::from_fn(|k| Dual { v: p.v[k], d: an_all[i][k] }))
                .collect();
            let an: Vec<[Dual; 3]> = (0..pred.len()).map(|i| std::array::from_fn(|k| Dual { v: an_all[i][k], d: jn_all[i][k] })).collect();
            for (idx, &i) in active.iter().enumerate() {
                let pn = eih_correction::<Dual>(i, &x, &v, &an, &self.mass, self.g, c);
                let (a_n, j_n, s_n) = res[idx];
                let a_pn: V3 = std::array::from_fn(|k| pn[k].v);
                let j_pn: V3 = std::array::from_fn(|k| pn[k].d);
                out.push((a_pn, j_pn, s_n, a_n, j_n));
            }
        } else {
            for (idx, _) in active.iter().enumerate() {
                let (a_n, j_n, s_n) = res[idx];
                out.push((a_n, j_n, s_n, a_n, j_n));
            }
        }
        Ok(out)
    }

    /// Components `(x1, v1, x2, v2)` of binary `b` at time `t` given the predicted state of its
    /// centre of mass: the internal KS state is advanced from the last update of the centre of
    /// mass without committing it.
    fn binary_components(
        &self,
        b: usize,
        cm: &Predicted,
        t: f64,
    ) -> Result<(V3, V3, V3, V3), String> {
        let bin = &self.binaries[b];
        let (m1, m2) = (self.comp_mass[bin.members[0]], self.comp_mass[bin.members[1]]);
        let (rel, rel_vel) = self.advance_binary(b, t)?;
        let mt = m1 + m2;
        let x1 = sub(cm.x, scale(rel, m2 / mt));
        let x2 = add(cm.x, scale(rel, m1 / mt));
        let v1 = sub(cm.v, scale(rel_vel, m2 / mt));
        let v2 = add(cm.v, scale(rel_vel, m1 / mt));
        Ok((x1, v1, x2, v2))
    }

    /// Relative state of binary `b` at time `t` by KS propagation under the tidal perturbation of
    /// the other particles (predicted along their Taylor polynomials).
    fn advance_binary(
        &self,
        b: usize,
        t: f64,
    ) -> Result<(V3, V3), String> {
        let bin = &self.binaries[b];
        let t0 = self.t_last[bin.cm];
        let dt = t - t0;
        if dt == 0.0 {
            return Ok((bin.rel, bin.rel_vel));
        }
        let (m1, m2) = (self.comp_mass[bin.members[0]], self.comp_mass[bin.members[1]]);
        let mt = m1 + m2;
        let mu = self.g * mt;
        let cm_idx = bin.cm;
        let others: Vec<usize> = (0..self.mass.len()).filter(|&k| k != cm_idx && !self.is_dead(k)).collect();
        let pert = |rel: &V3, _v: &V3, tau: f64| -> V3 {
            let tt = t0 + tau;
            let hc = tt - self.t_last[cm_idx];
            let cm_x: V3 = std::array::from_fn(|k| {
                self.pos[cm_idx][k] + hc * (self.vel[cm_idx][k] + hc * (self.acc[cm_idx][k] / 2.0 + hc * self.jerk[cm_idx][k] / 6.0))
            });
            let x1 = sub(cm_x, scale(*rel, m2 / mt));
            let x2 = add(cm_x, scale(*rel, m1 / mt));
            let mut p = [0.0; 3];
            for &k in &others {
                let h = tt - self.t_last[k];
                let xk: V3 = std::array::from_fn(|c| self.pos[k][c] + h * (self.vel[k][c] + h * (self.acc[k][c] / 2.0 + h * self.jerk[k][c] / 6.0)));
                let d1 = sub(xk, x1);
                let d2 = sub(xk, x2);
                let (r1, r2) = (norm(d1), norm(d2));
                for c in 0..3 {
                    p[c] += self.g * self.mass[k] * (d2[c] / (r2 * r2 * r2) - d1[c] / (r1 * r1 * r1));
                }
            }
            p
        };
        let (x, v) = crate::sim::symplectic::ks_advance(&bin.rel, &bin.rel_vel, mu, &pert, dt, bin.ks_eta)?;
        Ok((x, v))
    }

    /// Registers particles `a` and `b` (indices into the *current* particle arrays) as a
    /// KS-regularised binary: they are replaced by one centre-of-mass particle (stored at index
    /// `a`; index `b` becomes a massless passive placeholder that is excluded from all sums).
    ///
    /// # Errors
    /// Fails with the sixth-order scheme, PN, equal indices or out-of-range indices.
    pub fn add_binary(
        &mut self,
        a: usize,
        b: usize,
        ks_eta: f64,
    ) -> Result<(), String> {
        let n = self.mass.len();
        if a == b || a >= n || b >= n || self.order != HermiteOrder::Fourth || self.pn_c.is_some() {
            return Err("binaries need two distinct particles, the 4th-order scheme and no PN".into());
        }
        if (self.t_last[a] - self.t_last[b]).abs() > 1e-12 * (1.0 + self.t_last[a].abs()) {
            return Err("binary components must be synchronised".into());
        }
        let (m1, m2) = (self.mass[a], self.mass[b]);
        let mt = m1 + m2;
        let cm_x = scale(add(scale(self.pos[a], m1), scale(self.pos[b], m2)), 1.0 / mt);
        let cm_v = scale(add(scale(self.vel[a], m1), scale(self.vel[b], m2)), 1.0 / mt);
        let rel = sub(self.pos[b], self.pos[a]);
        let rel_vel = sub(self.vel[b], self.vel[a]);
        let ia = self.comp_mass.len();
        self.comp_mass.push(m1);
        self.comp_mass.push(m2);
        self.binaries.push(Binary { members: [ia, ia + 1], cm: a, partner: b, rel, rel_vel, ks_eta });
        // `b` is removed from the dynamics: zero mass at the centre of mass, never active
        self.mass[a] = mt;
        self.mass[b] = 0.0;
        self.pos[a] = cm_x;
        self.vel[a] = cm_v;
        self.pos[b] = cm_x;
        self.vel[b] = cm_v;
        self.t_last[b] = f64::INFINITY;
        self.initialise()
    }

    fn is_dead(&self, i: usize) -> bool {
        self.t_last[i].is_infinite()
    }

    /// Positions and velocities of all particles synchronised to time `t` (binary components are
    /// reported at the position of their regularised internal state; the placeholder particle of a
    /// binary reports its second component).
    ///
    /// # Errors
    /// Propagates binary propagation errors.
    pub fn state_at(
        &self,
        t: f64,
    ) -> Result<(Vec<V3>, Vec<V3>), String> {
        let n = self.mass.len();
        let mut x = vec![[0.0; 3]; n];
        let mut v = vec![[0.0; 3]; n];
        let mut slot = vec![None; n];
        for (b, bin) in self.binaries.iter().enumerate() {
            slot[bin.cm] = Some(b);
        }
        for i in 0..n {
            if self.is_dead(i) {
                continue;
            }
            let p = self.predict(i, t);
            if let Some(b) = slot[i] {
                let (x1, v1, x2, v2) = self.binary_components(b, &p, t)?;
                x[i] = x1;
                v[i] = v1;
                let k = self.binaries[b].partner;
                x[k] = x2;
                v[k] = v2;
            } else {
                x[i] = p.x;
                v[i] = p.v;
            }
        }
        Ok((x, v))
    }

    /// Total energy at the current system time (Newtonian, or the 1PN energy with PN enabled).
    ///
    /// # Errors
    /// Propagates binary propagation errors.
    pub fn total_energy(&self) -> Result<f64, String> {
        let t = self.time;
        let pred: Vec<Predicted> = (0..self.mass.len()).map(|i| if self.is_dead(i) { self.predict_dead(i) } else { self.predict(i, t) }).collect();
        let src = self.sources(&pred, t)?;
        let masses: Vec<f64> = src.iter().map(|s| s.2).collect();
        let pos: Vec<V3> = src.iter().map(|s| s.0).collect();
        let vel: Vec<V3> = src.iter().map(|s| s.1).collect();
        if let Some(c) = self.pn_c {
            return Ok(eih_energy(&masses, &pos, &vel, self.g, c));
        }
        let eps2 = self.eps * self.eps;
        let mut e = 0.0;
        for a in 0..masses.len() {
            e += 0.5 * masses[a] * dot(vel[a], vel[a]);
            for b in 0..a {
                let d = sub(pos[a], pos[b]);
                e -= self.g * masses[a] * masses[b] / (dot(d, d) + eps2).sqrt();
            }
        }
        Ok(e)
    }

    fn predict_dead(
        &self,
        i: usize,
    ) -> Predicted {
        Predicted { x: self.pos[i], v: self.vel[i], a: [0.0; 3], j: [0.0; 3], s: [0.0; 3] }
    }

    /// Time at which the next block ends.
    #[must_use]
    pub fn next_time(&self) -> f64 {
        (0..self.mass.len()).filter(|&i| !self.is_dead(i)).map(|i| self.t_last[i] + self.dt[i]).fold(f64::INFINITY, f64::min)
    }

    /// Advances the system by one block: the particles that are due at the earliest end time.
    ///
    /// # Errors
    /// Propagates binary propagation errors.
    pub fn step(&mut self) -> Result<(), String> {
        let n = self.mass.len();
        let t_new = self.next_time();
        let active: Vec<usize> = (0..n).filter(|&i| !self.is_dead(i) && (self.t_last[i] + self.dt[i] - t_new).abs() <= 1e-12 * self.dt[i]).collect();
        let mut pred: Vec<Predicted> = (0..n).map(|i| if self.is_dead(i) { self.predict_dead(i) } else { self.predict(i, t_new) }).collect();
        let base: Vec<Predicted> = pred.clone();
        let sixth = self.order == HermiteOrder::Sixth;
        let iters = if sixth { 1 } else { self.iterations.max(1) };
        let mut results = Vec::new();
        for _ in 0..iters {
            let out = self.forces(&pred, &active, t_new)?;
            results.clear();
            for (k, &i) in active.iter().enumerate() {
                let h = t_new - self.t_last[i];
                let (a1, j1, s1, _, _) = out[k];
                let (a0, j0, s0) = (self.acc[i], self.jerk[i], self.snap[i]);
                let mut x = [0.0; 3];
                let mut v = [0.0; 3];
                let mut snap1 = [0.0; 3];
                let mut crackle1 = [0.0; 3];
                for c in 0..3 {
                    if sixth {
                        v[c] = self.vel[i][c] + h / 2.0 * (a0[c] + a1[c]) - h * h / 10.0 * (j1[c] - j0[c]) + h * h * h / 120.0 * (s0[c] + s1[c]);
                        x[c] = self.pos[i][c] + h / 2.0 * (self.vel[i][c] + v[c]) - h * h / 10.0 * (a1[c] - a0[c]) + h * h * h / 120.0 * (j0[c] + j1[c]);
                        let r1 = (a1[c] - a0[c] - j0[c] * h - s0[c] * h * h / 2.0) / (h * h * h);
                        let r2 = (j1[c] - j0[c] - s0[c] * h) / (h * h);
                        let r3 = (s1[c] - s0[c]) / h;
                        let acoef = (r3 + 20.0 * r1 - 8.0 * r2) / 2.0;
                        snap1[c] = s1[c];
                        crackle1[c] = 18.0 * acoef - 120.0 * r1 + 36.0 * r2;
                    } else {
                        let snap0 = (-6.0 * (a0[c] - a1[c]) - h * (4.0 * j0[c] + 2.0 * j1[c])) / (h * h);
                        let crackle = (12.0 * (a0[c] - a1[c]) + 6.0 * h * (j0[c] + j1[c])) / (h * h * h);
                        x[c] = base[i].x[c] + h.powi(4) / 24.0 * snap0 + h.powi(5) / 120.0 * crackle;
                        v[c] = base[i].v[c] + h.powi(3) / 6.0 * snap0 + h.powi(4) / 24.0 * crackle;
                        snap1[c] = snap0 + h * crackle;
                        crackle1[c] = crackle;
                    }
                }
                results.push((i, x, v, a1, j1, snap1, crackle1, out[k].3, out[k].4));
            }
            // corrected states feed the next iteration
            for r in &results {
                pred[r.0].x = r.1;
                pred[r.0].v = r.2;
            }
        }
        let mut bin_updates = Vec::new();
        for r in &results {
            if let Some(b) = self.binaries.iter().position(|bin| bin.cm == r.0) {
                bin_updates.push((b, self.advance_binary(b, t_new)?));
            }
        }
        for (b, (rel, rel_vel)) in bin_updates {
            self.binaries[b].rel = rel;
            self.binaries[b].rel_vel = rel_vel;
        }
        for r in &results {
            let (i, x, v, a1, j1, s1, c1, an, jn) = *r;
            self.pos[i] = x;
            self.vel[i] = v;
            self.acc[i] = a1;
            self.jerk[i] = j1;
            self.snap[i] = s1;
            self.crackle[i] = c1;
            self.accn[i] = an;
            self.jerkn[i] = jn;
            self.t_last[i] = t_new;
            self.updates += 1;
            self.dt[i] = self.next_dt(i, t_new);
        }
        self.interactions += (active.len() * n) as u64;
        self.time = t_new;
        Ok(())
    }

    fn next_dt(
        &self,
        i: usize,
        t_new: f64,
    ) -> f64 {
        match self.policy {
            StepPolicy::Shared(h) => h,
            StepPolicy::Block { eta, dt_max, dt_min } => {
                let (a, j, s, c) = (norm(self.acc[i]), norm(self.jerk[i]), norm(self.snap[i]), norm(self.crackle[i]));
                let ss = s * s;
                let denom = j * c + ss;
                let jj = j * j;
                let aj = a * s + jj;
                let proposed = if denom > 0.0 { (eta * aj / denom).sqrt() } else { dt_max };
                let old = self.dt[i];
                let mut dt = old;
                if proposed < old {
                    while dt > proposed && dt > dt_min {
                        dt *= 0.5;
                    }
                } else if proposed >= 2.0 * old && 2.0 * old <= pow2_floor(dt_max) && (t_new / (2.0 * old)).fract() == 0.0 {
                    dt = 2.0 * old;
                }
                dt.max(dt_min)
            },
        }
    }

    /// Runs blocks until the next block would end after `t_end`.
    ///
    /// # Errors
    /// Propagates [`step`](Self::step) errors.
    pub fn advance_to(
        &mut self,
        t_end: f64,
    ) -> Result<(), String> {
        while self.next_time() <= t_end {
            self.step()?;
        }
        Ok(())
    }
}
