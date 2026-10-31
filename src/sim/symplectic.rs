//! # Composition symplectic integrators, time transformation, regularisation, rigid bodies
//!
//! Building blocks for separable Hamiltonians `H = |p|^2/2 + V(q)` in any
//! dimension (state `q`, `p` are flat slices of coordinates).
//!
//! | Capability | Dimensions | Notes |
//! |---|---|---|
//! | Composition schemes [`Scheme`]: Verlet (2), Forest-Ruth (4), Yoshida (4, 6, 8), Blanes-Moan (4), triple jump of any even order | any | stage tables, tested by convergence order |
//! | Time-transformed leapfrog [`ttl_step`] (Sundman-type, extended phase space, exactly reversible) and its composition to any even order [`ttl_step_composed`] / [`composition_weights`] | any | attractive potentials `V < 0`; order 4, 6, 8 remove the `O(h^2)` energy error (tested) |
//! | Chain regularisation (Mikkola-Aarseth chain), algorithmic regularisation with Gragg-Bulirsch-Stoer extrapolation | 3 | [`crate::sim::chain`] (`N` bodies, adaptive `ds`) |
//! | Kustaanheimo-Stiefel regularisation [`ks_propagate`] | 3 (via a 4-D oscillator) | perturbed Kepler problem, regular through collision |
//! | Rigid body (quaternion, NO_SQUISH / Dullweber-Leimkuhler-McLachlan splitting) [`RigidBody`], exact symmetric-top flows ([`RigidBody::free_step_top`], exact for symmetric bodies), composition to any even order ([`RigidBody::free_step_composed`]), space-torque kicks | 3 | symplectic, reversible, exact `\|q\| = 1` and `\|L\|` conservation in free rotation; rigid molecules in MD: [`crate::sim::md_rigid`] |
//!
//! | Exact asymmetric top [`RigidBody::free_flow_exact`] (Jacobi `sn, cn, dn` for the angular momentum, elliptic integral of the third kind [`crate::sim::elliptic::elliptic_pi`] for the precession angle), Strang torque step on the exact flow [`RigidBody::step_exact`] | 3 | exact for any time (checked against RK4 with tiny steps), group property, symmetric / steady cases handled |
//! | Implicit RATTLE-type torque step [`RigidBody::rattle_step`] (Moser-Veselov discrete rotor, Newton solve of `F J - J F^T = h Pi`) | 3 | symplectic, reversible, second order, exact space-`L` conservation |
//!
//! Limitations: the TTL scheme uses a fixed binding energy `B` (the energy of the exact flow
//! is preserved, but the plain scheme has `O(h^2)` energy error; use the composed or GBS
//! versions); the explicit torque kick is the Strang splitting of symplectic exact flows (reversible,
//! second order, composable to higher order; the implicit RATTLE alternative is second order only);
//! the exact asymmetric flow uses the symmetric-top flow when two moments agree to 1e-12; non-separable or perturbed problems with velocity-dependent forces are not
//! covered by the time-transformed schemes.

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
    /// Blanes-Moan `SRKNb6` (order 4, 6 force evaluations, small error constant).
    BlanesMoan4,
}

impl Scheme {
    /// Order of accuracy of the scheme.
    #[must_use]
    pub const fn order(self) -> usize {
        match self {
            | Self::Verlet => 2,
            | Self::ForestRuth | Self::Yoshida4 | Self::BlanesMoan4 => 4,
            | Self::Yoshida6 => 6,
            | Self::Yoshida8 => 8,
        }
    }

    /// Stage table: pairs `(drift, kick)` applied in sequence (drift first).
    #[must_use]
    pub fn stages(self) -> Vec<(f64, f64)> {
        match self {
            | Self::Verlet => vec![(0.5, 1.0), (0.5, 0.0)],
            | Self::ForestRuth => {
                let t = 1.0 / (2.0 - 2.0_f64.cbrt());
                vec![
                    (t / 2.0, t),
                    ((1.0 - t) / 2.0, 1.0 - 2.0 * t),
                    ((1.0 - t) / 2.0, t),
                    (t / 2.0, 0.0),
                ]
            },
            | Self::Yoshida4 => triple_jump(&Self::Verlet.stages(), 2),
            | Self::Yoshida6 => {
                let w = [-1.17767998417887, 0.235573213359357, 0.784513610477560];
                yoshida_weights(&w)
            },
            | Self::Yoshida8 => {
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
            | Self::BlanesMoan4 => {
                let a1 = 0.0792036964311957;
                let a2 = 0.353172906049774;
                let a3 = -0.0420650803577195;
                let a4 = 1.0 - 2.0 * (a1 + a2 + a3);
                let b1 = 0.209515106613362;
                let b2 = -0.143851773179818;
                let b3 = 0.5 - (b1 + b2);
                vec![
                    (a1, b1),
                    (a2, b2),
                    (a3, b3),
                    (a4, b3),
                    (a3, b2),
                    (a2, b1),
                    (a1, 0.0),
                ]
            },
        }
    }
}

/// Step-size multipliers of a symmetric composition of even `order`.
///
/// Step-size multipliers `c_i` (summing to 1) of the symmetric composition
/// `S(c_1 h) S(c_2 h) ...` that lifts a symmetric second-order step `S` to the even `order`
/// (Yoshida sequences for 6 and 8, Suzuki triple jumps beyond).
///
/// # Errors
/// Fails for odd or zero `order`.
pub fn composition_weights(order: usize) -> Result<Vec<f64>, String> {
    if order == 0 || !order.is_multiple_of(2) {
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
        | 2 => vec![1.0],
        | 4 => triple(&[1.0], 2),
        | 6 => seq(&[-1.17767998417887, 0.235573213359357, 0.784513610477560]),
        | 8 => {
            seq(&[
                -1.61582374150097,
                -2.44699182370524,
                -0.00716989419708120,
                2.44002732616735,
                0.157739928123617,
                1.82020630970714,
                1.04242620869991,
            ])
        },
        | _ => {
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

/// Time-transformed leapfrog composed to an even `order`.
///
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
        u[0] = f64::midpoint(r, x[0]).sqrt();
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

/// Propagates a (perturbed) Kepler problem in KS variables with RK4.
///
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
    KsState {
        u: [y[0], y[1], y[2], y[3]],
        up: [y[4], y[5], y[6], y[7]],
        t: y[9],
        energy: y[8],
    }
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
        [
            w * w + x * x - y * y - z * z,
            2.0 * (x * y - w * z),
            2.0 * (x * z + w * y),
        ],
        [
            2.0 * (x * y + w * z),
            w * w - x * x + y * y - z * z,
            2.0 * (y * z - w * x),
        ],
        [
            2.0 * (x * z - w * y),
            2.0 * (y * z + w * x),
            w * w - x * x - y * y + z * z,
        ],
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
        (0..3)
            .map(|k| self.l_body[k] * self.l_body[k] / (2.0 * self.inertia[k]))
            .sum()
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

    /// Rotation about body axis `k` by the angle `zeta` (exact flow of `H = c L_k^2 / 2`
    /// with `zeta = h c L_k`).
    fn rotate_by(
        &mut self,
        k: usize,
        zeta: f64,
    ) {
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

    /// Rotation about body axis `k` by `h * L_k / I_k` (exact flow of the split Hamiltonian).
    /// An infinite moment of inertia marks a degenerate axis of a linear body (no rotation).
    fn rotate_axis(
        &mut self,
        k: usize,
        h: f64,
    ) {
        let zeta = h * self.l_body[k] / self.inertia[k];
        self.rotate_by(k, zeta);
    }

    /// Exact flow for time `h` of the **symmetric top** obtained by replacing the two inverse
    /// moments about the axes other than `ax` with their mean: `L` is constant in space, so
    /// `q(h) = q_L(h) * q(0) * q_ax(h)` with `q_L` a rotation about `L` by `|L| h / I_perp` and
    /// `q_ax` a rotation about the symmetry axis by `h L_ax (1/I_ax - 1/I_perp)`.
    pub fn symmetric_top_flow(
        &mut self,
        ax: usize,
        h: f64,
    ) {
        let inv = [
            1.0 / self.inertia[0],
            1.0 / self.inertia[1],
            1.0 / self.inertia[2],
        ];
        let (a, b) = ((ax + 1) % 3, (ax + 2) % 3);
        let inv_perp = f64::midpoint(inv[a], inv[b]);
        let ls = self.l_space();
        let lmag = (ls[0] * ls[0] + ls[1] * ls[1] + ls[2] * ls[2]).sqrt();
        if lmag > 0.0 {
            let th = 0.5 * lmag * h * inv_perp;
            let (sn, cs) = th.sin_cos();
            let ql = [cs, sn * ls[0] / lmag, sn * ls[1] / lmag, sn * ls[2] / lmag];
            self.q = quat_mul(&ql, &self.q);
        }
        let ang = 0.5 * h * self.l_body[ax] * (inv[ax] - inv_perp);
        let (sn, cs) = ang.sin_cos();
        let mut e = [cs, 0.0, 0.0, 0.0];
        e[ax + 1] = sn;
        self.q = quat_mul(&self.q, &e);
        let n = self.q.iter().map(|v| v * v).sum::<f64>().sqrt();
        for v in &mut self.q {
            *v /= n;
        }
        let r = quat_to_matrix(&self.q);
        for k in 0..3 {
            self.l_body[k] = (0..3).map(|i| r[i][k] * ls[i]).sum();
        }
    }

    /// Axis whose two partners have the closest inverse moments (best symmetric-top approximation).
    fn symmetry_axis(&self) -> usize {
        let inv = [
            1.0 / self.inertia[0],
            1.0 / self.inertia[1],
            1.0 / self.inertia[2],
        ];
        (0..3)
            .min_by(|&x, &y| {
                let dx = (inv[(x + 1) % 3] - inv[(x + 2) % 3]).abs();
                let dy = (inv[(y + 1) % 3] - inv[(y + 2) % 3]).abs();
                dx.total_cmp(&dy)
            })
            .unwrap_or(0)
    }

    /// Second-order symmetric step of the free rotor built on exact **symmetric-top flows**:
    /// `exp(h/2 R_a) exp(h/2 R_b) exp(h S) exp(h/2 R_b) exp(h/2 R_a)` where `S` is the symmetric
    /// top and `R_a`, `R_b` the single-axis remainders. Exact when the body is itself a symmetric
    /// top and much more accurate than [`free_step`](Self::free_step) for nearly symmetric bodies.
    pub fn free_step_top(
        &mut self,
        h: f64,
    ) {
        let ax = self.symmetry_axis();
        let inv = [
            1.0 / self.inertia[0],
            1.0 / self.inertia[1],
            1.0 / self.inertia[2],
        ];
        let (a, b) = ((ax + 1) % 3, (ax + 2) % 3);
        let inv_perp = f64::midpoint(inv[a], inv[b]);
        for k in [a, b] {
            let z = 0.5 * h * (inv[k] - inv_perp) * self.l_body[k];
            self.rotate_by(k, z);
        }
        self.symmetric_top_flow(ax, h);
        for k in [b, a] {
            let z = 0.5 * h * (inv[k] - inv_perp) * self.l_body[k];
            self.rotate_by(k, z);
        }
    }

    /// Free rotation for time `h` with the chosen splitting composed to the even `order`
    /// (2, 4, 6, 8, ...) by [`composition_weights`].
    ///
    /// # Errors
    /// Fails for an invalid `order`.
    pub fn free_step_composed(
        &mut self,
        h: f64,
        top: bool,
        order: usize,
    ) -> Result<(), String> {
        for c in composition_weights(order)? {
            if top {
                self.free_step_top(c * h);
            } else {
                self.free_step(c * h);
            }
        }
        Ok(())
    }

    /// Adds `h * torque` (space-frame torque) to the angular momentum.
    pub fn kick_space(
        &mut self,
        h: f64,
        torque: &[f64; 3],
    ) {
        let r = quat_to_matrix(&self.q);
        for k in 0..3 {
            let mut tb = 0.0;
            for i in 0..3 {
                tb += r[i][k] * torque[i];
            }
            self.l_body[k] += h * tb;
        }
    }

    /// One symmetric (`NO_SQUISH` / DLM) step of the free rotor, second order.
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

// ---------------------------------------------------------------------------
// Exact free asymmetric top (Jacobi elliptic functions)
// ---------------------------------------------------------------------------

type M3 = [[f64; 3]; 3];

fn mat3_mul(
    a: &M3,
    b: &M3,
) -> M3 {
    let mut c = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = (0..3).map(|k| a[i][k] * b[k][j]).sum();
        }
    }
    c
}

fn mat3_t(a: &M3) -> M3 {
    let mut c = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = a[j][i];
        }
    }
    c
}

/// Rotation matrix about the axis `k` (0, 1, 2) by the angle `a`.
fn rot_axis(
    k: usize,
    a: f64,
) -> M3 {
    let (s, c) = a.sin_cos();
    let (i, j) = ((k + 1) % 3, (k + 2) % 3);
    let mut m = [[0.0; 3]; 3];
    m[k][k] = 1.0;
    m[i][i] = c;
    m[j][j] = c;
    m[j][i] = s;
    m[i][j] = -s;
    m
}

/// Unit quaternion (scalar first) of a proper rotation matrix.
#[must_use]
pub fn matrix_to_quat(r: &[[f64; 3]; 3]) -> [f64; 4] {
    let tr = r[0][0] + r[1][1] + r[2][2];
    let q = if tr > 0.0 {
        let s = 2.0 * (tr + 1.0).sqrt();
        [
            0.25 * s,
            (r[2][1] - r[1][2]) / s,
            (r[0][2] - r[2][0]) / s,
            (r[1][0] - r[0][1]) / s,
        ]
    } else if r[0][0] > r[1][1] && r[0][0] > r[2][2] {
        let s = 2.0 * (1.0 + r[0][0] - r[1][1] - r[2][2]).sqrt();
        [
            (r[2][1] - r[1][2]) / s,
            0.25 * s,
            (r[0][1] + r[1][0]) / s,
            (r[0][2] + r[2][0]) / s,
        ]
    } else if r[1][1] > r[2][2] {
        let s = 2.0 * (1.0 + r[1][1] - r[0][0] - r[2][2]).sqrt();
        [
            (r[0][2] - r[2][0]) / s,
            (r[0][1] + r[1][0]) / s,
            0.25 * s,
            (r[1][2] + r[2][1]) / s,
        ]
    } else {
        let s = 2.0 * (1.0 + r[2][2] - r[0][0] - r[1][1]).sqrt();
        [
            (r[1][0] - r[0][1]) / s,
            (r[0][2] + r[2][0]) / s,
            (r[1][2] + r[2][1]) / s,
            0.25 * s,
        ]
    };
    let n = q.iter().map(|x| x * x).sum::<f64>().sqrt();
    [q[0] / n, q[1] / n, q[2] / n, q[3] / n]
}

/// Constants of the Jacobi-function solution of the Euler equations.
struct EulerSolution {
    /// Axis of the `dn` component.
    d: usize,
    /// Axis of the `cn` component.
    p: usize,
    /// Axis of the `sn` component.
    s: usize,
    /// Amplitude of the `dn` component.
    ad: f64,
    /// Amplitude of the `cn` component.
    ap: f64,
    /// Amplitude of the `sn` component.
    asn: f64,
    /// Sign of the `dn` component.
    sign_d: f64,
    /// Parameter `m = k^2`.
    m: f64,
    /// Signed angular frequency of the Jacobi argument.
    omega: f64,
    /// Argument at `t = 0`.
    tau0: f64,
    /// Amplitude at `t = 0`.
    phi0: f64,
}

impl RigidBody {
    /// Solves the Euler equations of the free asymmetric top in closed form. Returns `None`
    /// for (near-)symmetric bodies and steady rotations, which are handled separately.
    fn euler_solution(&self) -> Option<EulerSolution> {
        let l = self.l_body;
        let w = [
            1.0 / self.inertia[0],
            1.0 / self.inertia[1],
            1.0 / self.inertia[2],
        ];
        let l2 = l.iter().map(|x| x * x).sum::<f64>();
        let two_t = (0..3).map(|k| w[k] * l[k] * l[k]).sum::<f64>();
        let mut ord = [0usize, 1, 2];
        ord.sort_by(|&x, &y| self.inertia[x].total_cmp(&self.inertia[y]));
        let [a, b, c] = ord;
        let scale = w[a].abs().max(w[c].abs());
        if (w[a] - w[b]).abs() < 1e-12 * scale || (w[b] - w[c]).abs() < 1e-12 * scale {
            return None;
        }
        let d = if l2 >= two_t * self.inertia[b] {
            c
        } else {
            a
        };
        let others: [usize; 2] = if d == c { [a, b] } else { [c, b] };
        let mut best: Option<EulerSolution> = None;
        for (p, s) in [(others[0], others[1]), (others[1], others[0])] {
            let ad2 = (two_t - w[p] * l2) / (w[d] - w[p]);
            let ap2 = (w[d] * l2 - two_t) / (w[d] - w[p]);
            let as2 = (two_t - w[d] * l2) / (w[s] - w[d]);
            if ad2 <= 0.0 || ap2 < 0.0 || as2 < 0.0 {
                continue;
            }
            let m = (as2 - ap2) / ad2;
            if !(-1e-12..=1.0 + 1e-12).contains(&m) {
                continue;
            }
            best = Some(EulerSolution {
                d,
                p,
                s,
                ad: ad2.sqrt(),
                ap: ap2.sqrt(),
                asn: as2.sqrt(),
                sign_d: if l[d] < 0.0 { -1.0 } else { 1.0 },
                m: m.clamp(0.0, 1.0),
                omega: 0.0,
                tau0: 0.0,
                phi0: 0.0,
            });
            break;
        }
        let mut sol = best?;
        if sol.ap < 1e-14 * l2.sqrt() || sol.asn < 1e-14 * l2.sqrt() {
            return None;
        }
        let m = sol.m;
        let sn0 = l[sol.s] / sol.asn;
        let cn0 = l[sol.p] / sol.ap;
        sol.phi0 = sn0.atan2(cn0);
        sol.tau0 = crate::kernels::special_ext::elliptic_f(sol.phi0, m);
        let dn0 = (1.0 - m * sn0 * sn0).max(0.0).sqrt();
        let om = [l[0] * w[0], l[1] * w[1], l[2] * w[2]];
        let ldot = [
            l[1] * om[2] - l[2] * om[1],
            l[2] * om[0] - l[0] * om[2],
            l[0] * om[1] - l[1] * om[0],
        ];
        let mag = sol.ad * ((w[sol.d] - w[sol.p]).abs() * (w[sol.s] - w[sol.d]).abs()).sqrt();
        let from_s = sol.asn * cn0.abs() * dn0;
        let from_p = sol.ap * sn0.abs() * dn0;
        let sign = if from_s >= from_p {
            ldot[sol.s] * cn0.signum()
        } else {
            -ldot[sol.p] * sn0.signum()
        };
        sol.omega = if sign < 0.0 { -mag } else { mag };
        Some(sol)
    }

    /// **Exact** free rotation for time `t` of an asymmetric rigid body, in closed form.
    ///
    /// With `d` the principal axis about which the polhode circulates, the body-frame angular
    /// momentum is `L_d = +-A_d dn`, `L_p = A_p cn`, `L_s = A_s sn` (Jacobi functions of
    /// `tau = tau_0 + omega t`, parameter `m`; Landau and Lifshitz, section 37). The orientation
    /// follows from the Euler angles with the (space-fixed) angular momentum as polar axis:
    /// `theta` and `psi` are algebraic in `L`, and the precession angle is
    /// `phi(t) = |L| t / I_d + |L| (2T - L^2 / I_d) int_0^t dt' / (L^2 - L_d^2)`, an elliptic
    /// integral of the **third kind** ([`crate::sim::elliptic::elliptic_pi`]). Symmetric tops and
    /// steady rotations use the exact symmetric-top flow. Valid for any `t` (including negative).
    pub fn free_flow_exact(
        &mut self,
        t: f64,
    ) {
        let l = self.l_body;
        let l2 = l.iter().map(|x| x * x).sum::<f64>();
        if l2 == 0.0 || t == 0.0 {
            return;
        }
        let lmag = l2.sqrt();
        for k in 0..3 {
            let (i, j) = ((k + 1) % 3, (k + 2) % 3);
            if l[i].hypot(l[j]) < 1e-14 * lmag {
                self.rotate_by(k, t * l[k] / self.inertia[k]);
                return;
            }
        }
        let w = [
            1.0 / self.inertia[0],
            1.0 / self.inertia[1],
            1.0 / self.inertia[2],
        ];
        let Some(sol) = self.euler_solution() else {
            let ax = self.symmetry_axis();
            self.symmetric_top_flow(ax, t);
            return;
        };
        let two_t = (0..3).map(|k| w[k] * l[k] * l[k]).sum::<f64>();
        let m = sol.m;
        let tau = sol.omega.mul_add(t, sol.tau0);
        let (sn, cn, dn) = crate::kernels::special_ext::jacobi_elliptic(tau, m);
        let phi_mod = sn.atan2(cn);
        let kk = crate::kernels::special_ext::elliptic_k(m);
        let n_wrap =
            ((tau - crate::kernels::special_ext::elliptic_f(phi_mod, m)) / (4.0 * kk)).round();
        let phi_t = phi_mod + 2.0 * std::f64::consts::PI * n_wrap;
        let mut ln = [0.0; 3];
        ln[sol.p] = sol.ap * cn;
        ln[sol.s] = sol.asn * sn;
        ln[sol.d] = sol.sign_d * sol.ad * dn;
        let nu = sol.ad * sol.ad * m / (sol.ap * sol.ap);
        let j_int = (crate::sim::elliptic::elliptic_pi(-nu, phi_t, m)
            - crate::sim::elliptic::elliptic_pi(-nu, sol.phi0, m))
            / (sol.omega * sol.ap * sol.ap);
        let dphi = lmag * w[sol.d] * t + lmag * (two_t - w[sol.d] * l2) * j_int;
        let d = sol.d;
        let (a1, a2) = ((d + 1) % 3, (d + 2) % 3);
        let angles = |lv: &[f64; 3]| ((lv[d] / lmag).clamp(-1.0, 1.0).acos(), lv[a1].atan2(lv[a2]));
        let (th0, psi0) = angles(&l);
        let (th1, psi1) = angles(&ln);
        let mut pm = [[0.0; 3]; 3];
        pm[a1][0] = 1.0;
        pm[a2][1] = 1.0;
        pm[d][2] = 1.0;
        let r0 = quat_to_matrix(&self.q);
        // R(t) = R0 P Rz(-psi0) Rx(-th0) Rz(dphi) Rx(th1) Rz(psi1) P^T
        let mut r = mat3_mul(&r0, &pm);
        for (axis, ang) in [(2, -psi0), (0, -th0), (2, dphi), (0, th1), (2, psi1)] {
            r = mat3_mul(&r, &rot_axis(axis, ang));
        }
        r = mat3_mul(&r, &mat3_t(&pm));
        let mut q = matrix_to_quat(&r);
        if q.iter().zip(&self.q).map(|(x, y)| x * y).sum::<f64>() < 0.0 {
            for v in &mut q {
                *v = -*v;
            }
        }
        self.q = q;
        self.l_body = ln;
    }

    /// Strang step with space-frame torque built on the **exact** free flow
    /// ([`free_flow_exact`](Self::free_flow_exact)): half kick, exact free rotation, half kick.
    /// Symplectic, time reversible and second order, with no splitting error in the free part.
    pub fn step_exact<T: Fn(&[f64; 4]) -> [f64; 3]>(
        &mut self,
        h: f64,
        torque: &T,
    ) {
        self.kick(0.5 * h, torque);
        self.free_flow_exact(h);
        self.kick(0.5 * h, torque);
    }
}

impl RigidBody {
    /// Rotation `F` (as a vector-parametrised proper rotation) solving the Moser-Veselov
    /// equation `F J - J F^T = h Pi^` for the discrete angular velocity.
    fn moser_veselov_rotation(
        &self,
        h: f64,
        pi: &[f64; 3],
    ) -> Result<[f64; 3], String> {
        let i = self.inertia;
        let jd = [
            0.5 * (i[1] + i[2] - i[0]),
            0.5 * (i[2] + i[0] - i[1]),
            0.5 * (i[0] + i[1] - i[2]),
        ];
        let resid = |v: &[f64; 3]| {
            let f = rotation_from_vector(v);
            // S = F J - J F^T is skew; its axial vector is (S21, S02, S10)
            let mut s = [[0.0; 3]; 3];
            for a in 0..3 {
                for b in 0..3 {
                    s[a][b] = f[a][b] * jd[b] - jd[a] * f[b][a];
                }
            }
            [
                s[2][1] - h * pi[0],
                s[0][2] - h * pi[1],
                s[1][0] - h * pi[2],
            ]
        };
        let mut v = [h * pi[0] / i[0], h * pi[1] / i[1], h * pi[2] / i[2]];
        let scale = (h * (pi[0].abs() + pi[1].abs() + pi[2].abs())).max(1e-300);
        for _ in 0..60 {
            let r = resid(&v);
            let rn = r[0].abs() + r[1].abs() + r[2].abs();
            if rn <= 1e-15 * scale {
                return Ok(v);
            }
            // finite-difference Jacobian
            let mut jac = [[0.0; 3]; 3];
            for c in 0..3 {
                let eps = 1e-7 * (1.0 + v[c].abs());
                let mut vp = v;
                vp[c] += eps;
                let mut vm = v;
                vm[c] -= eps;
                let (rp, rm) = (resid(&vp), resid(&vm));
                for a in 0..3 {
                    jac[a][c] = (rp[a] - rm[a]) / (2.0 * eps);
                }
            }
            let dv = solve3(&jac, &r).ok_or("singular Moser-Veselov Jacobian")?;
            for c in 0..3 {
                v[c] -= dv[c];
            }
        }
        let r = resid(&v);
        if r[0].abs() + r[1].abs() + r[2].abs() <= 1e-11 * scale {
            Ok(v)
        } else {
            Err("Moser-Veselov iteration did not converge".to_string())
        }
    }

    /// One **implicit RATTLE-type** step with a space-frame torque: the discrete rigid-body
    /// (Moser-Veselov) free step, in which the rotation increment `F` solves the implicit
    /// constraint equation `F J - J F^T = h Pi^` (`J = tr(I)/2 - I` in the principal frame,
    /// Newton iteration), sandwiched between half kicks of the torque. The map is symplectic,
    /// time reversible and second order; it conserves the space-frame angular momentum exactly
    /// in the absence of torque and keeps `R^T R = 1` to round-off.
    ///
    /// # Errors
    /// Fails if the Newton iteration for the rotation increment does not converge (very large
    /// steps) or a principal moment is not positive.
    pub fn rattle_step<T: Fn(&[f64; 4]) -> [f64; 3]>(
        &mut self,
        h: f64,
        torque: &T,
    ) -> Result<(), String> {
        if self.inertia.iter().any(|&x| x <= 0.0) {
            return Err("rattle_step needs positive principal moments".to_string());
        }
        self.kick(0.5 * h, torque);
        let pi = self.l_body;
        let v = self.moser_veselov_rotation(h, &pi)?;
        let f = rotation_from_vector(&v);
        let angle = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        let (s, c) = (0.5 * angle).sin_cos();
        let k = if angle > 0.0 {
            s / angle
        } else {
            0.5
        };
        let e = [c, k * v[0], k * v[1], k * v[2]];
        let mut q = quat_mul(&self.q, &e);
        let n = q.iter().map(|x| x * x).sum::<f64>().sqrt();
        for x in &mut q {
            *x /= n;
        }
        self.q = q;
        // Pi_{n+1} = F^T Pi_n
        let mut lb = [0.0; 3];
        for a in 0..3 {
            lb[a] = (0..3).map(|b| f[b][a] * pi[b]).sum();
        }
        self.l_body = lb;
        self.kick(0.5 * h, torque);
        Ok(())
    }
}

/// Rotation matrix `exp(v^)` (Rodrigues).
fn rotation_from_vector(v: &[f64; 3]) -> M3 {
    let th = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    let (a, b) = if th < 1e-6 {
        (1.0 - th * th / 6.0, 0.5 - th * th / 24.0)
    } else {
        (th.sin() / th, (1.0 - th.cos()) / (th * th))
    };
    let k = [[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]];
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let kk: f64 = (0..3).map(|m| k[i][m] * k[m][j]).sum();
            r[i][j] = f64::from(u8::from(i == j)) + a * k[i][j] + b * kk;
        }
    }
    r
}

/// Solves the 3x3 system `a x = b` by Cramer's rule.
fn solve3(
    a: &M3,
    b: &[f64; 3],
) -> Option<[f64; 3]> {
    let det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    if det.abs() < 1e-300 {
        return None;
    }
    let mut x = [0.0; 3];
    for c in 0..3 {
        let mut m = *a;
        for r in 0..3 {
            m[r][c] = b[r];
        }
        let d = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
        x[c] = d / det;
    }
    Some(x)
}

/// Advances a perturbed two-body problem in KS variables for a given physical time.
///
/// The gravitational parameter is `mu`, the physical time `dt > 0` (RK4 in the fictitious time, step `ds = eta / omega` with the
/// harmonic frequency `omega` of the regularised oscillator; the last steps are solved for so
/// that the physical time lands exactly on `dt`). `pert(x, v, tau)` is the perturbing relative
/// acceleration at elapsed time `tau`. Returns the relative position and velocity.
///
/// # Errors
/// Fails for a non-positive `dt`, `mu` or `eta`, or a collision (`r = 0`).
pub fn ks_advance<P: Fn(&[f64; 3], &[f64; 3], f64) -> [f64; 3]>(
    x0: &[f64; 3],
    v0: &[f64; 3],
    mu: f64,
    pert: &P,
    dt: f64,
    eta: f64,
) -> Result<([f64; 3], [f64; 3]), String> {
    if dt <= 0.0 || mu <= 0.0 || eta <= 0.0 {
        return Err("ks_advance needs positive dt, mu and eta".into());
    }
    let r0 = (x0[0] * x0[0] + x0[1] * x0[1] + x0[2] * x0[2]).sqrt();
    if r0 == 0.0 {
        return Err("ks_advance: zero separation".into());
    }
    let (u, up) = ks_from_cartesian(x0, v0);
    let e0 = 0.5 * (v0[0] * v0[0] + v0[1] * v0[1] + v0[2] * v0[2]) - mu / r0;
    let mut y = [u[0], u[1], u[2], u[3], up[0], up[1], up[2], up[3], e0, 0.0];
    let rk4 = |y: &[f64; 10], ds: f64| {
        let k1 = ks_rhs(y, pert);
        let mut t = *y;
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
        let mut out = *y;
        for i in 0..10 {
            out[i] += ds / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
        out
    };
    for _ in 0..10_000_000 {
        let remaining = dt - y[9];
        if !remaining.is_finite() {
            return Err("ks_advance: non-finite state".into());
        }
        if remaining.abs() <= 1e-14 * dt {
            break;
        }
        let r = y[0] * y[0] + y[1] * y[1] + y[2] * y[2] + y[3] * y[3];
        if r <= 0.0 {
            return Err("ks_advance: collision".into());
        }
        let e = y[8];
        let omega = if e < -1e-12 * mu / r {
            (-0.5 * e).sqrt()
        } else {
            (0.25 * mu / r).sqrt()
        };
        let ds_nominal = eta / omega;
        let ds = if remaining.abs() < r * ds_nominal * 1.05 {
            remaining / r
        } else {
            ds_nominal
        };
        y = rk4(&y, ds);
    }
    let (x, v) = ks_to_cartesian(&[y[0], y[1], y[2], y[3]], &[y[4], y[5], y[6], y[7]]);
    Ok((x, v))
}
