//! # Composition symplectic integrators, time transformation, regularisation, rigid bodies
//!
//! Building blocks for separable Hamiltonians `H = |p|^2/2 + V(q)` in any
//! dimension (state `q`, `p` are flat slices of coordinates).
//!
//! | Capability | Dimensions | Notes |
//! |---|---|---|
//! | Composition schemes [`Scheme`]: Verlet (2), Forest-Ruth (4), Yoshida (4, 6, 8), Blanes-Moan (4), triple jump of any even order | any | stage tables, tested by convergence order |
//! | Time-transformed leapfrog [`ttl_step`] (Sundman-type, extended phase space, exactly reversible) | any | attractive potentials `V < 0` |
//! | Kustaanheimo-Stiefel regularisation [`ks_propagate`] | 3 (via a 4-D oscillator) | perturbed Kepler problem, regular through collision |
//! | Rigid body (quaternion, NO_SQUISH / Dullweber-Leimkuhler-McLachlan splitting) [`RigidBody`] | 3 | symplectic, reversible, exact `|q| = 1` and `|L|` conservation in free rotation |
//!
//! Limitations: the TTL scheme uses a fixed binding energy `B` (the energy
//! error is `O(h^2)`, not zero); chain regularisation is not implemented
//! (KS is); the rigid-body torque kick is a second-order reversible
//! splitting rather than an implicit RATTLE solve.
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

/// A symplectic composition scheme for separable Hamiltonians.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Scheme {
    /// Velocity Verlet (order 2).
    Verlet,
    /// Forest-Ruth (order 4, 3 force evaluations).
    ForestRuth,
    /// Yoshida triple jump (order 4).
    Yoshida4,
    /// Yoshida solution A (order 6, 7 stages).
    Yoshida6,
    /// Yoshida solution A (order 8, 15 stages).
    Yoshida8,
    /// Blanes-Moan SRKNb6 (order 4, 6 force evaluations, small error constant).
    BlanesMoan4,
}

impl Scheme {
    /// Order of accuracy of the scheme.
    #[must_use]
    pub const fn order(self) -> usize {
        match self {
            Self::Verlet => 2,
            Self::ForestRuth | Self::Yoshida4 | Self::BlanesMoan4 => 4,
            Self::Yoshida6 => 6,
            Self::Yoshida8 => 8,
        }
    }

    /// Stage table: pairs `(drift, kick)` applied in sequence (drift first).
    #[must_use]
    pub fn stages(self) -> Vec<(f64, f64)> {
        match self {
            Self::Verlet => vec![(0.5, 1.0), (0.5, 0.0)],
            Self::ForestRuth => {
                let t = 1.0 / (2.0 - 2.0_f64.cbrt());
                vec![(t / 2.0, t), ((1.0 - t) / 2.0, 1.0 - 2.0 * t), ((1.0 - t) / 2.0, t), (t / 2.0, 0.0)]
            },
            Self::Yoshida4 => triple_jump(&Self::Verlet.stages(), 2),
            Self::Yoshida6 => {
                let w = [-1.17767998417887, 0.235573213359357, 0.784513610477560];
                yoshida_weights(&w)
            },
            Self::Yoshida8 => {
                let w = [
                    -1.61582374150097,
                    -2.44699182370524,
                    -0.00716989419708120,
                    2.44002732616735,
                    0.157739928123617,
                    1.82020630970714,
                    1.04242620869991,
                ];
                yoshida_weights(&w)
            },
            Self::BlanesMoan4 => {
                let a1 = 0.0792036964311957;
                let a2 = 0.353172906049774;
                let a3 = -0.0420650803577195;
                let a4 = 1.0 - 2.0 * (a1 + a2 + a3);
                let b1 = 0.209515106613362;
                let b2 = -0.143851773179818;
                let b3 = 0.5 - (b1 + b2);
                vec![(a1, b1), (a2, b2), (a3, b3), (a4, b3), (a3, b2), (a2, b1), (a1, 0.0)]
            },
        }
    }
}

/// Step-size multipliers `c_i` (summing to 1) of the symmetric composition
/// `S(c_1 h) S(c_2 h) ...` that lifts a symmetric second-order step `S` to the even `order`
/// (Yoshida sequences for 6 and 8, Suzuki triple jumps beyond).
///
/// # Errors
/// Fails for odd or zero `order`.
pub fn composition_weights(order: usize) -> Result<Vec<f64>, String> {
    if order == 0 || order % 2 != 0 {
        return Err("composition order must be even and positive".into());
    }
    let seq = |w: &[f64]| -> Vec<f64> {
        let w0 = 1.0 - 2.0 * w.iter().sum::<f64>();
        let mut s: Vec<f64> = w.iter().rev().copied().collect();
        s.push(w0);
        s.extend(w.iter().copied());
        s
    };
    let triple = |base: &[f64], p: usize| -> Vec<f64> {
        let g = 1.0 / (2.0 - 2.0_f64.powf(1.0 / (p as f64 + 1.0)));
        let mid = 1.0 - 2.0 * g;
        let mut out = Vec::new();
        for c in [g, mid, g] {
            out.extend(base.iter().map(|w| w * c));
        }
        out
    };
    Ok(match order {
        2 => vec![1.0],
        4 => triple(&[1.0], 2),
        6 => seq(&[-1.17767998417887, 0.235573213359357, 0.784513610477560]),
        8 => seq(&[-1.61582374150097, -2.44699182370524, -0.00716989419708120, 2.44002732616735, 0.157739928123617, 1.82020630970714, 1.04242620869991]),
        _ => {
            let mut w = composition_weights(8)?;
            let mut p = 8;
            while p < order {
                w = triple(&w, p);
                p += 2;
            }
            w
        },
    })
}

/// Time-transformed leapfrog composed to the even `order` (4, 6, 8, ...): the TTL step is a
/// symmetric second-order integrator of the extended Hamiltonian `ln(T + B) - ln(-V)`, so
/// composition with [`composition_weights`] removes the `O(h^2)` energy error of the plain scheme.
///
/// # Errors
/// Fails for an invalid `order`.
pub fn ttl_step_composed<F>(
    s: &mut TtlState,
    mass: &[f64],
    b: f64,
    h: f64,
    order: usize,
    potential_grad: &mut F,
) -> Result<(), String>
where
    F: FnMut(&[f64], &mut [f64]) -> f64,
{
    for c in composition_weights(order)? {
        ttl_step(s, mass, b, c * h, potential_grad);
    }
    Ok(())
}

/// Symmetric composition `S(w_k h) ... S(w_1 h) S(w_0 h) S(w_1 h) ... S(w_k h)` of Verlet.
fn yoshida_weights(w: &[f64]) -> Vec<(f64, f64)> {
    let w0 = 1.0 - 2.0 * w.iter().sum::<f64>();
    let mut seq: Vec<f64> = w.iter().rev().copied().collect();
    seq.push(w0);
    seq.extend(w.iter().copied());
    let base = Scheme::Verlet.stages();
    let mut out = Vec::new();
    for c in seq {
        out.extend(base.iter().map(|&(a, b)| (a * c, b * c)));
    }
    merge_stages(out)
}

fn merge_stages(stages: Vec<(f64, f64)>) -> Vec<(f64, f64)> {
    let mut out: Vec<(f64, f64)> = Vec::with_capacity(stages.len());
    let mut pending_drift = 0.0;
    for (a, b) in stages {
        pending_drift += a;
        if b != 0.0 {
            out.push((pending_drift, b));
            pending_drift = 0.0;
        }
    }
    if pending_drift != 0.0 {
        out.push((pending_drift, 0.0));
    }
    out
}

/// Suzuki triple jump: lifts a symmetric scheme of order `p` to order `p + 2`.
#[must_use]
pub fn triple_jump(
    base: &[(f64, f64)],
    p: usize,
) -> Vec<(f64, f64)> {
    let g = 1.0 / (2.0 - 2.0_f64.powf(1.0 / (p as f64 + 1.0)));
    let mid = 1.0 - 2.0 * g;
    let mut out = Vec::new();
    for c in [g, mid, g] {
        out.extend(base.iter().map(|&(a, b)| (a * c, b * c)));
    }
    merge_stages(out)
}

/// One step of a composition scheme; `acc(q, a)` writes the acceleration.
pub fn step<F>(
    stages: &[(f64, f64)],
    q: &mut [f64],
    p: &mut [f64],
    h: f64,
    acc: &mut F,
) where
    F: FnMut(&[f64], &mut [f64]),
{
    let mut a = vec![0.0; q.len()];
    for &(d, k) in stages {
        for i in 0..q.len() {
            q[i] += d * h * p[i];
        }
        if k != 0.0 {
            acc(q, &mut a);
            for i in 0..q.len() {
                p[i] += k * h * a[i];
            }
        }
    }
}

/// Integrates `n` steps of size `h` with `scheme`.
pub fn integrate<F>(
    scheme: Scheme,
    q: &mut [f64],
    p: &mut [f64],
    h: f64,
    n: usize,
    acc: &mut F,
) where
    F: FnMut(&[f64], &mut [f64]),
{
    let st = scheme.stages();
    for _ in 0..n {
        step(&st, q, p, h, acc);
    }
}

// ---------------------------------------------------------------------------
// Time-transformed leapfrog (Mikkola-Tanikawa)
// ---------------------------------------------------------------------------

/// State of the time-transformed leapfrog (`t` is the physical time).
#[derive(Clone, Debug)]
pub struct TtlState {
    /// Positions.
    pub q: Vec<f64>,
    /// Momenta.
    pub p: Vec<f64>,
    /// Physical time.
    pub t: f64,
}

/// Time-transformed leapfrog for `H = T(p) + V(q)` with `V < 0`.
///
/// Uses the extended Hamiltonian `K = ln(T + B) - ln(-V)` with `B = -E0`
/// (`E0` the energy to preserve); the fictitious step `h` corresponds to the
/// physical step `h / (T + B)`, which is `h r / (GM)` for Kepler (a Sundman
/// transformation `dt = r ds`). The scheme is symmetric and symplectic in
/// the extended phase space. `mass` holds one mass per coordinate;
/// `potential_grad(q, g)` writes `grad V` and returns `V`.
pub fn ttl_step<F>(
    s: &mut TtlState,
    mass: &[f64],
    b: f64,
    h: f64,
    potential_grad: &mut F,
) where
    F: FnMut(&[f64], &mut [f64]) -> f64,
{
    let n = s.q.len();
    let half_drift = |s: &mut TtlState| {
        let kin: f64 = (0..n).map(|i| s.p[i] * s.p[i] / (2.0 * mass[i])).sum();
        let w = 1.0 / (kin + b);
        for i in 0..n {
            s.q[i] += 0.5 * h * w * s.p[i] / mass[i];
        }
        s.t += 0.5 * h * w;
    };
    half_drift(s);
    let mut g = vec![0.0; n];
    let v = potential_grad(&s.q, &mut g);
    for i in 0..n {
        s.p[i] += h * g[i] / v;
    }
    half_drift(s);
}

// ---------------------------------------------------------------------------
// Kustaanheimo-Stiefel regularisation
// ---------------------------------------------------------------------------

/// The KS matrix `L(u)` (4x4, rows give `x1, x2, x3, 0`).
#[must_use]
pub fn ks_matrix(u: &[f64; 4]) -> [[f64; 4]; 4] {
    [
        [u[0], -u[1], -u[2], u[3]],
        [u[1], u[0], -u[3], -u[2]],
        [u[2], u[3], u[0], u[1]],
        [u[3], -u[2], u[1], -u[0]],
    ]
}

/// KS map of a regular state: returns `(u, u')` for position `x` and velocity `v`
/// (`u'` is the derivative with respect to the fictitious time `ds = dt / r`).
#[must_use]
pub fn ks_from_cartesian(
    x: &[f64; 3],
    v: &[f64; 3],
) -> ([f64; 4], [f64; 4]) {
    let r = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
    let mut u = [0.0; 4];
    if x[0] >= 0.0 {
        u[0] = (0.5 * (r + x[0])).sqrt();
        u[3] = 0.0;
        u[1] = (x[1] * u[0] + x[2] * u[3]) / (r + x[0]);
        u[2] = (x[2] * u[0] - x[1] * u[3]) / (r + x[0]);
    } else {
        u[1] = (0.5 * (r - x[0])).sqrt();
        u[2] = 0.0;
        u[0] = (x[1] * u[1] + x[2] * u[2]) / (r - x[0]);
        u[3] = (x[2] * u[1] - x[1] * u[2]) / (r - x[0]);
    }
    // x' = 2 L(u) u'  =>  u' = L(u)^T [r v, 0] / (2 r)  since L L^T = r I.
    let l = ks_matrix(&u);
    let rv = [r * v[0], r * v[1], r * v[2], 0.0];
    let mut up = [0.0; 4];
    for i in 0..4 {
        for k in 0..4 {
            up[i] += l[k][i] * rv[k];
        }
        up[i] /= 2.0 * r;
    }
    (u, up)
}

/// Cartesian position and velocity from a KS state.
#[must_use]
pub fn ks_to_cartesian(
    u: &[f64; 4],
    up: &[f64; 4],
) -> ([f64; 3], [f64; 3]) {
    let l = ks_matrix(u);
    let r = u.iter().map(|c| c * c).sum::<f64>();
    let mut x = [0.0; 3];
    let mut v = [0.0; 3];
    for i in 0..3 {
        for k in 0..4 {
            x[i] += l[i][k] * u[k];
            v[i] += 2.0 * l[i][k] * up[k] / r;
        }
    }
    (x, v)
}

/// Result of a KS propagation: final regular state, physical time and energy.
#[derive(Clone, Debug)]
pub struct KsState {
    /// KS coordinates.
    pub u: [f64; 4],
    /// `du/ds`.
    pub up: [f64; 4],
    /// Physical time.
    pub t: f64,
    /// Specific energy `v^2/2 - 1/r`.
    pub energy: f64,
}

fn ks_rhs<P: Fn(&[f64; 3], &[f64; 3], f64) -> [f64; 3]>(
    y: &[f64; 10],
    pert: &P,
) -> [f64; 10] {
    let u = [y[0], y[1], y[2], y[3]];
    let up = [y[4], y[5], y[6], y[7]];
    let e = y[8];
    let r = u.iter().map(|c| c * c).sum::<f64>();
    let (x, v) = ks_to_cartesian(&u, &up);
    let pf = pert(&x, &v, y[9]);
    let l = ks_matrix(&u);
    let mut ltp = [0.0; 4];
    for i in 0..4 {
        for k in 0..3 {
            ltp[i] += l[k][i] * pf[k];
        }
    }
    // dE/ds = 2 (L u')_{1..3} . P
    let mut de = 0.0;
    for i in 0..3 {
        let mut lu = 0.0;
        for k in 0..4 {
            lu += l[i][k] * up[k];
        }
        de += 2.0 * lu * pf[i];
    }
    let mut d = [0.0; 10];
    for i in 0..4 {
        d[i] = up[i];
        d[4 + i] = 0.5 * e * u[i] + 0.5 * r * ltp[i];
    }
    d[8] = de;
    d[9] = r;
    d
}

/// Propagates a (perturbed) Kepler problem with gravitational parameter 1 in KS
/// variables with RK4 in the fictitious time `s`, for `n` steps of size `ds`.
/// `pert(x, v, t)` is the perturbing acceleration. Without perturbation the
/// equations are a harmonic oscillator of frequency `sqrt(-E/2)`, so
/// near-collision orbits are integrated without singularity.
#[must_use]
pub fn ks_propagate<P: Fn(&[f64; 3], &[f64; 3], f64) -> [f64; 3]>(
    x0: &[f64; 3],
    v0: &[f64; 3],
    pert: &P,
    ds: f64,
    n: usize,
) -> KsState {
    let (u, up) = ks_from_cartesian(x0, v0);
    let r0 = (x0[0] * x0[0] + x0[1] * x0[1] + x0[2] * x0[2]).sqrt();
    let e0 = 0.5 * (v0[0] * v0[0] + v0[1] * v0[1] + v0[2] * v0[2]) - 1.0 / r0;
    let mut y = [u[0], u[1], u[2], u[3], up[0], up[1], up[2], up[3], e0, 0.0];
    for _ in 0..n {
        let k1 = ks_rhs(&y, pert);
        let mut t = y;
        for i in 0..10 {
            t[i] = y[i] + 0.5 * ds * k1[i];
        }
        let k2 = ks_rhs(&t, pert);
        for i in 0..10 {
            t[i] = y[i] + 0.5 * ds * k2[i];
        }
        let k3 = ks_rhs(&t, pert);
        for i in 0..10 {
            t[i] = y[i] + ds * k3[i];
        }
        let k4 = ks_rhs(&t, pert);
        for i in 0..10 {
            y[i] += ds / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
    }
    KsState { u: [y[0], y[1], y[2], y[3]], up: [y[4], y[5], y[6], y[7]], t: y[9], energy: y[8] }
}

// ---------------------------------------------------------------------------
// Rigid body
// ---------------------------------------------------------------------------

/// Quaternion product `a * b` (scalar first).
#[must_use]
pub fn quat_mul(
    a: &[f64; 4],
    b: &[f64; 4],
) -> [f64; 4] {
    [
        a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
        a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
        a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1],
        a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0],
    ]
}

/// Rotation matrix (body to space) of a unit quaternion.
#[must_use]
pub fn quat_to_matrix(q: &[f64; 4]) -> [[f64; 3]; 3] {
    let [w, x, y, z] = *q;
    [
        [w * w + x * x - y * y - z * z, 2.0 * (x * y - w * z), 2.0 * (x * z + w * y)],
        [2.0 * (x * y + w * z), w * w - x * x + y * y - z * z, 2.0 * (y * z - w * x)],
        [2.0 * (x * z - w * y), 2.0 * (y * z + w * x), w * w - x * x - y * y + z * z],
    ]
}

/// A rigid body with principal moments of inertia, orientation and body-frame angular momentum.
#[derive(Clone, Debug)]
pub struct RigidBody {
    /// Principal moments of inertia.
    pub inertia: [f64; 3],
    /// Orientation quaternion (body to space), kept at unit norm.
    pub q: [f64; 4],
    /// Angular momentum in the body frame.
    pub l_body: [f64; 3],
}

impl RigidBody {
    /// Rotational kinetic energy.
    #[must_use]
    pub fn kinetic_energy(&self) -> f64 {
        (0..3).map(|k| self.l_body[k] * self.l_body[k] / (2.0 * self.inertia[k])).sum()
    }

    /// Angular momentum in the space frame.
    #[must_use]
    pub fn l_space(&self) -> [f64; 3] {
        let r = quat_to_matrix(&self.q);
        let mut l = [0.0; 3];
        for i in 0..3 {
            for k in 0..3 {
                l[i] += r[i][k] * self.l_body[k];
            }
        }
        l
    }

    /// Rotation about body axis `k` by `h * L_k / I_k` (exact flow of the split Hamiltonian).
    fn rotate_axis(
        &mut self,
        k: usize,
        h: f64,
    ) {
        let zeta = h * self.l_body[k] / self.inertia[k];
        let (s, c) = (0.5 * zeta).sin_cos();
        let mut e = [c, 0.0, 0.0, 0.0];
        e[k + 1] = s;
        self.q = quat_mul(&self.q, &e);
        let n = self.q.iter().map(|v| v * v).sum::<f64>().sqrt();
        for v in &mut self.q {
            *v /= n;
        }
        let (sz, cz) = zeta.sin_cos();
        let (i, j) = ((k + 1) % 3, (k + 2) % 3);
        let (li, lj) = (self.l_body[i], self.l_body[j]);
        self.l_body[i] = cz * li + sz * lj;
        self.l_body[j] = -sz * li + cz * lj;
    }

    /// One symmetric (NO_SQUISH / DLM) step of the free rotor, second order.
    pub fn free_step(
        &mut self,
        h: f64,
    ) {
        self.rotate_axis(2, 0.5 * h);
        self.rotate_axis(1, 0.5 * h);
        self.rotate_axis(0, h);
        self.rotate_axis(1, 0.5 * h);
        self.rotate_axis(2, 0.5 * h);
    }

    /// One step with a space-frame torque `torque(q)`: half kick, free rotation, half kick.
    pub fn step<T: Fn(&[f64; 4]) -> [f64; 3]>(
        &mut self,
        h: f64,
        torque: &T,
    ) {
        self.kick(0.5 * h, torque);
        self.free_step(h);
        self.kick(0.5 * h, torque);
    }

    fn kick<T: Fn(&[f64; 4]) -> [f64; 3]>(
        &mut self,
        h: f64,
        torque: &T,
    ) {
        let ts = torque(&self.q);
        let r = quat_to_matrix(&self.q);
        for k in 0..3 {
            let mut tb = 0.0;
            for i in 0..3 {
                tb += r[i][k] * ts[i];
            }
            self.l_body[k] += h * tb;
        }
    }
}
