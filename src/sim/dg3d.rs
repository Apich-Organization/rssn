//! # Nodal discontinuous Galerkin on hexahedra (3D)
//!
//! Collocation (DGSEM) discontinuous Galerkin on a tensor-product mesh of
//! hexahedral elements: Legendre–Gauss–Lobatto nodes of degree `p >= 1` per
//! axis (`(p + 1)³` nodes per element), strong form with surface lifting,
//! SSP-RK3 time stepping and the same [`DgLaw`] interface as the 1D/2D solver
//! of [`dg`](crate::sim::dg), so every law with a flux in three directions
//! works.
//!
//! * laws: [`LinearAdvection3d`], [`Euler3d`] (`ρ, ρu, ρv, ρw, E`);
//! * interface fluxes: Rusanov (default) or **entropy stable**
//!   (flux-differencing volume terms with the Chandrashekar entropy-conservative
//!   two-point flux and Lax–Friedrichs dissipation, [`Dg3Options::entropy_stable`]);
//! * boundaries: periodic, outflow, slip wall, fixed (inflow) state;
//! * shock limiters: characteristic, hierarchical and WENO
//!   ([`Dg3Options::limiter`]), applied after every SSP-RK3 stage;
//! * Zhang–Shu positivity-preserving scaling of density and pressure.
//!
//! Every limiter preserves the element averages, so a periodic computation
//! keeps total mass, momentum and energy to round-off. The three limiters are
//! the 3D tensor-product analogues of the 1D/2D limiters of
//! [`dg`](crate::sim::dg):
//!
//! | [`Dg3Limiter`] | Method |
//! |---|---|
//! | [`Dg3Limiter::Characteristic`] | TVB-tested characteristic moment limiter: troubled elements are detected from the face traces and the neighbour means, then the Legendre modes are bounded from degree `p` down in the characteristic variables of the normal Jacobian ([`DgLaw::to_characteristic`], Euler eigenvectors for [`Euler3d`]) |
//! | [`Dg3Limiter::Hierarchical`] | the same moment limiter, componentwise, on the elements flagged by the componentwise TVB test |
//! | [`Dg3Limiter::Weno`] | Zhong–Shu WENO reconstruction of troubled elements along every axis (element plus left/right neighbours, componentwise) |
//!
//! Smooth problems converge with order `p + 1` (advection) and `p`–`p + 1`
//! (Euler vortex). Limitations: elements are axis-aligned boxes of one
//! uniform size (no curved or unstructured hexahedra, no hanging nodes); the
//! moment limiters use ghost elements built from the boundary condition at
//! non-periodic faces (mirror image for walls, the exterior state for fixed
//! states, even mirror image for outflow (the flux's exterior state is the interior trace); a direction is skipped only on an
//! axis with a single element), the characteristic limiter needs the law's eigenvectors for the full benefit
//! (it falls back to the componentwise moment limiter otherwise), and the
//! WENO limiter is componentwise (its directional reconstructions are
//! averaged, not a full multidimensional WENO reconstruction).
//!
//! Parallelism: the element-local right-hand side and the troubled-element test run concurrently
//! over the elements.

use crate::sim::dg::DgLaw;
use crate::sim::dg::euler_ec_flux;
use crate::sim::dg::euler_prim;
use crate::sim::dg::euler_wall_state;
use crate::sim::dg::invert;
use crate::sim::dg::legendre;
use crate::sim::dg::lgl;
use crate::sim::dg::lgl_diff_matrix;
use std::sync::Arc;

/// Linear advection `u_t + a·∇u = 0` in 3D.
#[derive(Clone, Copy, Debug)]
pub struct LinearAdvection3d {
    /// Advection velocity.
    pub a: [f64; 3],
}

impl DgLaw for LinearAdvection3d {
    fn nvar(&self) -> usize {
        1
    }

    fn flux(
        &self,
        q: &[f64],
        dir: usize,
        out: &mut [f64],
    ) {
        out[0] = self.a[dir] * q[0];
    }

    fn max_speed(
        &self,
        _q: &[f64],
        dir: usize,
    ) -> f64 {
        self.a[dir].abs()
    }
}

/// Three-dimensional Euler equations `(ρ, ρu, ρv, ρw, E)` of an ideal gas.
#[derive(Clone, Copy, Debug)]
pub struct Euler3d {
    /// Adiabatic exponent.
    pub gamma: f64,
}

impl Euler3d {
    /// Pressure of a conserved state.
    #[must_use]
    pub fn pressure(
        &self,
        q: &[f64],
    ) -> f64 {
        (self.gamma - 1.0) * (q[4] - 0.5 * (q[1] * q[1] + q[2] * q[2] + q[3] * q[3]) / q[0])
    }
}

impl DgLaw for Euler3d {
    fn nvar(&self) -> usize {
        5
    }

    fn flux(
        &self,
        q: &[f64],
        dir: usize,
        out: &mut [f64],
    ) {
        let p = self.pressure(q);
        let vel = q[1 + dir] / q[0];
        out[0] = q[1 + dir];
        for k in 0..3 {
            out[1 + k] = q[1 + k] * vel;
        }
        out[1 + dir] += p;
        out[4] = (q[4] + p) * vel;
    }

    fn max_speed(
        &self,
        q: &[f64],
        dir: usize,
    ) -> f64 {
        let p = self.pressure(q).max(0.0);
        (q[1 + dir] / q[0]).abs() + (self.gamma * p / q[0]).sqrt()
    }

    fn rho_p(
        &self,
        q: &[f64],
    ) -> Option<(f64, f64)> {
        Some((q[0], self.pressure(q)))
    }

    fn has_characteristics(&self) -> bool {
        true
    }

    fn to_characteristic(
        &self,
        q_ref: &[f64],
        dir: usize,
        v: &[f64],
        out: &mut [f64],
    ) {
        euler3_eigen_transform(self.gamma, dir, q_ref, v, out, true);
    }

    fn from_characteristic(
        &self,
        q_ref: &[f64],
        dir: usize,
        w: &[f64],
        out: &mut [f64],
    ) {
        euler3_eigen_transform(self.gamma, dir, q_ref, w, out, false);
    }

    fn two_point_flux(
        &self,
        ql: &[f64],
        qr: &[f64],
        dir: usize,
        out: &mut [f64],
    ) {
        euler_ec_flux(self.gamma, 3, dir, ql, qr, out);
    }

    fn wall_state(
        &self,
        q: &[f64],
        dir: usize,
        out: &mut [f64],
    ) {
        euler_wall_state(3, dir, q, out);
    }
}

/// Applies the left (`to_char`) or right eigenvector matrix of the 3D Euler
/// flux Jacobian in direction `dir`, evaluated at the state `q`, to `v`.
///
/// The wave order is `(u_n - c, u_n (shear 1), u_n (shear 2), u_n (entropy),
/// u_n + c)`, i.e. the 2D transform of [`dg`](crate::sim::dg) extended by the
/// second transverse velocity.
fn euler3_eigen_transform(
    gamma: f64,
    dir: usize,
    q: &[f64],
    v: &[f64],
    out: &mut [f64],
    to_char: bool,
) {
    let (rho, vel, p) = euler_prim(gamma, 3, q);
    let c = (gamma * p.max(1e-300) / rho).sqrt();
    let un = vel[dir];
    let (t1, t2) = match dir {
        | 0 => (1, 2),
        | 1 => (0, 2),
        | _ => (0, 1),
    };
    let (ut1, ut2) = (vel[t1], vel[t2]);
    let q2 = un * un + ut1 * ut1 + ut2 * ut2;
    let h = (q[4] + p) / rho;
    let b1 = (gamma - 1.0) / (c * c);
    let b2 = 0.5 * b1 * q2;
    if to_char {
        // conserved components in sweep order: (ρ, ρu_n, ρu_t1, ρu_t2, E)
        let w = [v[0], v[1 + dir], v[1 + t1], v[1 + t2], v[4]];
        out[0] = f64::midpoint(b2, un / c) * w[0]
            - f64::midpoint(b1 * un, 1.0 / c) * w[1]
            - 0.5 * b1 * ut1 * w[2]
            - 0.5 * b1 * ut2 * w[3]
            + 0.5 * b1 * w[4];
        out[1] = -ut1 * w[0] + w[2];
        out[2] = -ut2 * w[0] + w[3];
        out[3] = (1.0 - b2) * w[0] + b1 * un * w[1] + b1 * ut1 * w[2] + b1 * ut2 * w[3] - b1 * w[4];
        out[4] = 0.5 * (b2 - un / c) * w[0]
            - 0.5 * (b1 * un - 1.0 / c) * w[1]
            - 0.5 * b1 * ut1 * w[2]
            - 0.5 * b1 * ut2 * w[3]
            + 0.5 * b1 * w[4];
    } else {
        // `v` is already in wave order (acoustic-, shear 1, shear 2, entropy, acoustic+)
        let d = v[0] + v[3] + v[4];
        out[0] = d;
        out[1 + dir] = (un - c) * v[0] + un * v[3] + (un + c) * v[4];
        out[1 + t1] = ut1 * d + v[1];
        out[1 + t2] = ut2 * d + v[2];
        out[4] =
            (h - un * c) * v[0] + ut1 * v[1] + ut2 * v[2] + 0.5 * q2 * v[3] + (h + un * c) * v[4];
    }
}

/// Specific thermodynamic entropy `s = ln(p / ρ^γ)` of an Euler state with `nvel` velocity components.
#[must_use]
pub fn euler_entropy(
    gamma: f64,
    nvel: usize,
    q: &[f64],
) -> f64 {
    let (rho, _, p) = euler_prim(gamma, nvel, q);
    (p / rho.powf(gamma)).ln()
}

/// Boundary treatment of one face of the box.
#[derive(Clone)]
pub enum Dg3Bc {
    /// Periodic (set on both faces of the axis).
    Periodic,
    /// Zero-gradient outflow.
    Outflow,
    /// Slip wall ([`DgLaw::wall_state`]).
    Wall,
    /// Fixed exterior state `q(t, x)`.
    Fixed(Arc<dyn Fn(f64, [f64; 3]) -> Vec<f64> + Send + Sync>),
}

/// Uniform box mesh of hexahedra.
#[derive(Clone, Copy, Debug)]
pub struct Dg3Mesh {
    /// Elements per axis.
    pub ne: [usize; 3],
    /// Lower corner.
    pub lo: [f64; 3],
    /// Upper corner.
    pub hi: [f64; 3],
    /// Polynomial degree.
    pub p: usize,
}

impl Dg3Mesh {
    fn h(
        &self,
        a: usize,
    ) -> f64 {
        (self.hi[a] - self.lo[a]) / self.ne[a] as f64
    }

    /// Number of elements.
    #[must_use]
    pub const fn elements(&self) -> usize {
        self.ne[0] * self.ne[1] * self.ne[2]
    }

    /// Nodes per element.
    #[must_use]
    pub const fn nodes_per_element(&self) -> usize {
        (self.p + 1) * (self.p + 1) * (self.p + 1)
    }
}

/// Shock-capturing limiter applied after every SSP-RK3 stage.
///
/// All three limiters preserve the element averages and therefore the total
/// mass, momentum and energy of a periodic computation to round-off. They run
/// before the Zhang–Shu positivity scaling of [`Dg3Options::positivity`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Dg3Limiter {
    /// No limiting (default, preserves the plain DGSEM scheme).
    #[default]
    None,
    /// Characteristic moment limiter: elements flagged by the TVB test have
    /// their Legendre modes bounded from degree `p` down to degree `1` in the
    /// characteristic variables of the flux Jacobian ([`DgLaw::to_characteristic`];
    /// Euler eigenvectors for [`Euler3d`], componentwise for laws without them).
    Characteristic,
    /// Krivodonova's hierarchical moment limiter on the Legendre coefficients,
    /// componentwise, along all three directions, applied to the elements
    /// flagged by the componentwise TVB troubled-element test.
    Hierarchical,
    /// Zhong–Shu WENO reconstruction of the troubled elements along each axis,
    /// componentwise, using the element and its left/right neighbours.
    Weno,
}

/// Options of [`solve_dg3d`].
#[derive(Clone, Copy, Debug)]
pub struct Dg3Options {
    /// Entropy-stable flux differencing (needs [`DgLaw::two_point_flux`]).
    pub entropy_stable: bool,
    /// Positivity-preserving scaling of density and pressure.
    pub positivity: bool,
    /// Shock-capturing limiter (default [`Dg3Limiter::None`], i.e. the plain
    /// DGSEM scheme).
    pub limiter: Dg3Limiter,
    /// TVB constant `M` of the troubled-element test shared by all three
    /// limiters (default `50`): jumps smaller than `M h²` are treated as
    /// smooth and left untouched, which keeps smooth extrema intact.
    pub limiter_tvd: f64,
    /// Relative regularisation `ε` of the WENO oscillation indicator
    /// (default `1e-6`, i.e. `ε = 1e-6 * mean²`).
    pub limiter_weno_eps: f64,
}

impl Default for Dg3Options {
    fn default() -> Self {
        Self {
            entropy_stable: false,
            positivity: false,
            limiter: Dg3Limiter::None,
            limiter_tvd: 50.0,
            limiter_weno_eps: 1e-6,
        }
    }
}

/// Solution snapshot.
#[derive(Clone, Debug)]
pub struct Dg3Solution {
    /// Mesh.
    pub mesh: Dg3Mesh,
    /// Reference nodes.
    pub xi: Vec<f64>,
    /// Reference weights.
    pub w: Vec<f64>,
    /// Values, layout `[element][node][variable]` with node index `i + n1 (j + n1 k)`.
    pub q: Vec<f64>,
    /// Variables.
    pub nvar: usize,
    /// Final time.
    pub t: f64,
    /// Number of steps.
    pub steps: usize,
}

impl Dg3Solution {
    const fn node_ijk(
        &self,
        node: usize,
    ) -> [usize; 3] {
        let n1 = self.mesh.p + 1;
        [node % n1, (node / n1) % n1, node / (n1 * n1)]
    }

    const fn elem_ijk(
        &self,
        e: usize,
    ) -> [usize; 3] {
        let ne = self.mesh.ne;
        [e % ne[0], (e / ne[0]) % ne[1], e / (ne[0] * ne[1])]
    }

    /// Physical coordinates of a node.
    #[must_use]
    pub fn coord(
        &self,
        e: usize,
        node: usize,
    ) -> [f64; 3] {
        let m = &self.mesh;
        let (ij, ex) = (self.node_ijk(node), self.elem_ijk(e));
        std::array::from_fn(|a| {
            m.lo[a] + (ex[a] as f64 + f64::midpoint(1.0, self.xi[ij[a]])) * m.h(a)
        })
    }

    /// Value of a variable at a node.
    #[must_use]
    pub fn value(
        &self,
        e: usize,
        node: usize,
        var: usize,
    ) -> f64 {
        self.q[(e * self.mesh.nodes_per_element() + node) * self.nvar + var]
    }

    fn quad_weight(
        &self,
        node: usize,
    ) -> f64 {
        let ij = self.node_ijk(node);
        let jac: f64 = (0..3).map(|a| 0.5 * self.mesh.h(a)).product();
        jac * (0..3).map(|a| self.w[ij[a]]).product::<f64>()
    }

    /// Quadrature `L2` norm of `q_var - exact`.
    pub fn l2_error(
        &self,
        var: usize,
        exact: impl Fn([f64; 3]) -> f64,
    ) -> f64 {
        let m = &self.mesh;
        let mut s = 0.0;
        for e in 0..m.elements() {
            for nd in 0..m.nodes_per_element() {
                let d = self.value(e, nd, var) - exact(self.coord(e, nd));
                s += self.quad_weight(nd) * d * d;
            }
        }
        s.sqrt()
    }

    /// Integral of a function of the state over the domain.
    pub fn integral(
        &self,
        f: impl Fn(&[f64]) -> f64,
    ) -> f64 {
        let m = &self.mesh;
        let nv = self.nvar;
        let mut s = 0.0;
        for e in 0..m.elements() {
            for nd in 0..m.nodes_per_element() {
                let o = (e * m.nodes_per_element() + nd) * nv;
                s += self.quad_weight(nd) * f(&self.q[o..o + nv]);
            }
        }
        s
    }
}

/// Elements from which the DG right-hand side is evaluated in parallel.
const RHS_PAR_ELEMENTS: usize = 8;

struct Ctx<'a> {
    mesh: Dg3Mesh,
    law: &'a dyn DgLaw,
    bcs: &'a [Dg3Bc; 6],
    d: Vec<Vec<f64>>,
    w: Vec<f64>,
    xi: Vec<f64>,
    nv: usize,
    np: usize,
    n1: usize,
    es: bool,
    limiter: Dg3Limiter,
    tvd: f64,
    weno_eps: f64,
    positivity: bool,
    vmat: Vec<Vec<f64>>,
    vinv: Vec<Vec<f64>>,
}

impl Ctx<'_> {
    const fn neighbour(
        &self,
        e: usize,
        axis: usize,
        max: bool,
    ) -> Option<usize> {
        let ne = self.mesh.ne;
        let mut c = [e % ne[0], (e / ne[0]) % ne[1], e / (ne[0] * ne[1])];
        let periodic = matches!(self.bcs[2 * axis], Dg3Bc::Periodic);
        if max {
            if c[axis] + 1 < ne[axis] {
                c[axis] += 1;
            } else if periodic {
                c[axis] = 0;
            } else {
                return None;
            }
        } else if c[axis] > 0 {
            c[axis] -= 1;
        } else if periodic {
            c[axis] = ne[axis] - 1;
        } else {
            return None;
        }
        Some(c[0] + ne[0] * (c[1] + ne[1] * c[2]))
    }

    fn rhs(
        &self,
        q: &[f64],
        t: f64,
    ) -> Vec<f64> {
        let m = &self.mesh;
        let (nv, np, n1) = (self.nv, self.np, self.n1);
        let mut r = vec![0.0; q.len()];
        let at = |e: usize, nd: usize| -> &[f64] {
            let o = (e * np + nd) * nv;
            &q[o..o + nv]
        };
        let node = |i: usize, j: usize, k: usize| i + n1 * (j + n1 * k);
        // Each element writes only its own block of `r` and reads `q`, so the
        // elements are independent: they are processed in parallel (the
        // result is identical to the serial loop, element by element).
        crate::kernels::par::for_chunks_mut(&mut r, np * nv, RHS_PAR_ELEMENTS, |e, re| {
            let mut fbuf = vec![0.0; nv];
            let mut flux = vec![0.0; np * nv];
            for a in 0..3 {
                let sc = 2.0 / m.h(a);
                if !self.es {
                    for nd in 0..np {
                        self.law.flux(at(e, nd), a, &mut fbuf);
                        flux[nd * nv..(nd + 1) * nv].copy_from_slice(&fbuf);
                    }
                }
                // index of the line node `s` for transverse indices (t1, t2)
                let idx = |s: usize, t1: usize, t2: usize| {
                    match a {
                        | 0 => node(s, t1, t2),
                        | 1 => node(t1, s, t2),
                        | _ => node(t1, t2, s),
                    }
                };
                for t2 in 0..n1 {
                    for t1 in 0..n1 {
                        for i in 0..n1 {
                            let nd = idx(i, t1, t2);
                            for k in 0..n1 {
                                let nk = idx(k, t1, t2);
                                if self.es {
                                    self.law.two_point_flux(at(e, nd), at(e, nk), a, &mut fbuf);
                                    for c in 0..nv {
                                        re[nd * nv + c] -= 2.0 * sc * self.d[i][k] * fbuf[c];
                                    }
                                } else {
                                    for c in 0..nv {
                                        re[nd * nv + c] -= sc * self.d[i][k] * flux[nk * nv + c];
                                    }
                                }
                            }
                        }
                        for max in [false, true] {
                            let edge = if max { n1 - 1 } else { 0 };
                            let nd = idx(edge, t1, t2);
                            let qi = at(e, nd);
                            let mut qo = vec![0.0; nv];
                            match self.neighbour(e, a, max) {
                                | Some(en) => {
                                    qo.copy_from_slice(at(en, idx(n1 - 1 - edge, t1, t2)));
                                },
                                | None => {
                                    match &self.bcs[2 * a + usize::from(max)] {
                                        | Dg3Bc::Fixed(f) => {
                                            let ex = [
                                                e % m.ne[0],
                                                (e / m.ne[0]) % m.ne[1],
                                                e / (m.ne[0] * m.ne[1]),
                                            ];
                                            let ij = match a {
                                                | 0 => [edge, t1, t2],
                                                | 1 => [t1, edge, t2],
                                                | _ => [t1, t2, edge],
                                            };
                                            let x: [f64; 3] = std::array::from_fn(|b| {
                                                m.lo[b]
                                                    + (ex[b] as f64
                                                        + f64::midpoint(1.0, self.xi[ij[b]]))
                                                        * m.h(b)
                                            });
                                            qo.copy_from_slice(&f(t, x));
                                        },
                                        | Dg3Bc::Wall => self.law.wall_state(qi, a, &mut qo),
                                        | _ => qo.copy_from_slice(qi),
                                    }
                                },
                            }
                            let mut fstar = vec![0.0; nv];
                            match (max, self.es) {
                                | (true, false) => self.law.numerical_flux(qi, &qo, a, &mut fstar),
                                | (false, false) => self.law.numerical_flux(&qo, qi, a, &mut fstar),
                                | (true, true) => {
                                    self.law.entropy_stable_flux(qi, &qo, a, &mut fstar);
                                },
                                | (false, true) => {
                                    self.law.entropy_stable_flux(&qo, qi, a, &mut fstar);
                                },
                            }
                            self.law.flux(qi, a, &mut fbuf);
                            let sign = if max { -1.0 } else { 1.0 };
                            for c in 0..nv {
                                re[nd * nv + c] += sign * sc * (fstar[c] - fbuf[c]) / self.w[edge];
                            }
                        }
                    }
                }
            }
        });
        r
    }

    fn max_rate(
        &self,
        q: &[f64],
    ) -> f64 {
        q.chunks(self.nv)
            .map(|c| {
                (0..3)
                    .map(|a| self.law.max_speed(c, a) / self.mesh.h(a))
                    .sum::<f64>()
            })
            .fold(0.0, f64::max)
    }

    fn node_weight(
        &self,
        nd: usize,
    ) -> f64 {
        let n1 = self.n1;
        let ij = [nd % n1, (nd / n1) % n1, nd / (n1 * n1)];
        (0..3).map(|a| self.w[ij[a]] / 2.0).product()
    }

    /// Zhang–Shu scaling of density and pressure towards the element mean.
    fn positivity(
        &self,
        q: &mut [f64],
    ) {
        let (nv, np) = (self.nv, self.np);
        let eps = 1e-10;
        for e in 0..self.mesh.elements() {
            let mean: Vec<f64> = (0..nv)
                .map(|v| {
                    (0..np)
                        .map(|nd| self.node_weight(nd) * q[(e * np + nd) * nv + v])
                        .sum()
                })
                .collect();
            let Some((rm, pm)) = self.law.rho_p(&mean) else {
                return;
            };
            if rm <= eps || pm <= eps {
                for nd in 0..np {
                    for v in 0..nv {
                        q[(e * np + nd) * nv + v] = mean[v];
                    }
                }
                continue;
            }
            let rho_min = (0..np)
                .map(|nd| q[(e * np + nd) * nv])
                .fold(f64::INFINITY, f64::min);
            let th1 = if rho_min < eps {
                ((rm - eps) / (rm - rho_min)).clamp(0.0, 1.0)
            } else {
                1.0
            };
            for nd in 0..np {
                for v in 0..nv {
                    let i = (e * np + nd) * nv + v;
                    q[i] = mean[v] + th1 * (q[i] - mean[v]);
                }
            }
            let mut th2: f64 = 1.0;
            for nd in 0..np {
                let base = (e * np + nd) * nv;
                let node = q[base..base + nv].to_vec();
                let ok = |s: f64| {
                    let st: Vec<f64> = (0..nv).map(|v| mean[v] + s * (node[v] - mean[v])).collect();
                    self.law
                        .rho_p(&st)
                        .is_some_and(|(r, p)| r >= eps && p >= eps)
                };
                if !ok(1.0) {
                    let (mut lo, mut hi) = (0.0, 1.0);
                    for _ in 0..40 {
                        let mid = f64::midpoint(lo, hi);
                        if ok(mid) {
                            lo = mid;
                        } else {
                            hi = mid;
                        }
                    }
                    th2 = th2.min(lo);
                }
            }
            if th2 < 1.0 {
                for nd in 0..np {
                    for v in 0..nv {
                        let i = (e * np + nd) * nv + v;
                        q[i] = mean[v] + th2 * (q[i] - mean[v]);
                    }
                }
            }
        }
    }
}

/// Smallest-magnitude argument when all three arguments share a sign, zero otherwise.
fn minmod3(
    a: f64,
    b: f64,
    c: f64,
) -> f64 {
    if a > 0.0 && b > 0.0 && c > 0.0 {
        a.min(b).min(c)
    } else if a < 0.0 && b < 0.0 && c < 0.0 {
        a.max(b).max(c)
    } else {
        0.0
    }
}

/// Value at reference coordinate `x` of the interpolant through the LGL nodes
/// `xi` with nodal values `f` (the 1D WENO stencil shift of the 1D/2D solver).
fn lagrange_eval(
    xi: &[f64],
    f: &[f64],
    x: f64,
) -> f64 {
    let n = xi.len();
    let mut s = 0.0;
    for j in 0..n {
        let mut l = 1.0;
        for k in 0..n {
            if k != j {
                l *= (x - xi[k]) / (xi[j] - xi[k]);
            }
        }
        s += f[j] * l;
    }
    s
}

impl Ctx<'_> {
    /// Node index on a line of an element: `s` along `axis`, `(t1, t2)` transverse.
    const fn axis_node(
        &self,
        axis: usize,
        s: usize,
        t1: usize,
        t2: usize,
    ) -> usize {
        let n1 = self.n1;
        match axis {
            | 0 => s + n1 * (t1 + n1 * t2),
            | 1 => t1 + n1 * (s + n1 * t2),
            | _ => t1 + n1 * (t2 + n1 * s),
        }
    }

    fn element_mean(
        &self,
        q: &[f64],
        e: usize,
        var: usize,
    ) -> f64 {
        (0..self.np)
            .map(|nd| self.node_weight(nd) * q[(e * self.np + nd) * self.nv + var])
            .sum()
    }

    /// Element averages, layout `[element][variable]`.
    fn means(
        &self,
        q: &[f64],
    ) -> Vec<f64> {
        let (nv, ne) = (self.nv, self.mesh.elements());
        let mut out = vec![0.0; ne * nv];
        for e in 0..ne {
            for v in 0..nv {
                out[e * nv + v] = self.element_mean(q, e, v);
            }
        }
        out
    }

    /// Transverse-weighted average of a variable over one face of an element.
    fn face_average(
        &self,
        q: &[f64],
        e: usize,
        axis: usize,
        upper: bool,
        var: usize,
    ) -> f64 {
        let s = if upper { self.n1 - 1 } else { 0 };
        let mut acc = 0.0;
        for t2 in 0..self.n1 {
            for t1 in 0..self.n1 {
                let nd = self.axis_node(axis, s, t1, t2);
                acc += 0.25 * self.w[t1] * self.w[t2] * q[(e * self.np + nd) * self.nv + var];
            }
        }
        acc
    }

    /// Applies the 1D matrix `mat` along every line of an element (tensor-product
    /// nodal-to-modal transform, layout `[node][variable]` with `node` only).
    fn transform_element(
        &self,
        vals: &[f64],
        mat: &[Vec<f64>],
    ) -> Vec<f64> {
        let n1 = self.n1;
        let mut cur = vals.to_vec();
        let mut nxt = vec![0.0; self.np];
        for axis in 0..3 {
            for t2 in 0..n1 {
                for t1 in 0..n1 {
                    for i in 0..n1 {
                        let nd = self.axis_node(axis, i, t1, t2);
                        nxt[nd] = (0..n1)
                            .map(|j| mat[i][j] * cur[self.axis_node(axis, j, t1, t2)])
                            .sum();
                    }
                }
            }
            std::mem::swap(&mut cur, &mut nxt);
        }
        cur
    }

    /// Legendre coefficients of one element (per-element layout `[node]`).
    fn to_modal(
        &self,
        vals: &[f64],
    ) -> Vec<f64> {
        self.transform_element(vals, &self.vinv)
    }

    /// Nodal values from Legendre coefficients (per-element layout `[node]`).
    #[allow(clippy::wrong_self_convention)] // inverse of `to_modal`, follows the 1D/2D solver
    fn from_modal(
        &self,
        c: &[f64],
    ) -> Vec<f64> {
        self.transform_element(c, &self.vmat)
    }

    /// TVB troubled-element indicator used by the three limiters.
    ///
    /// Following the 1D/2D characteristic limiter, the element is compared with
    /// its neighbours through the face traces `u_L`, `u_R` and the element
    /// averages, componentwise or in characteristic variables when
    /// `characteristic` is set and the law provides them. An element is
    /// troubled when one of its face slopes differs from the minmod of the
    /// forward/backward average slopes by more than the TVB threshold `M h²`.
    fn troubled_elements(
        &self,
        q: &[f64],
        mt: f64,
        characteristic: bool,
        t: f64,
    ) -> Vec<bool> {
        let (nv, ne) = (self.nv, self.mesh.elements());
        let m = &self.mesh;
        let means = self.means(q);
        let has_char = characteristic && self.law.has_characteristics();
        // the test of one element only reads `q` and `means`
        crate::kernels::par::map_indexed(ne, RHS_PAR_ELEMENTS * 4, |e| {
            let mut trouble_e = false;
            let mut tmp = vec![0.0; nv];
            let mut up = vec![0.0; nv];
            let mut um = vec![0.0; nv];
            let mut dp = vec![0.0; nv];
            let mut dm = vec![0.0; nv];
            let mean = &means[e * nv..(e + 1) * nv];
            for axis in 0..3 {
                for v in 0..nv {
                    up[v] = self.face_average(q, e, axis, true, v) - mean[v];
                    um[v] = mean[v] - self.face_average(q, e, axis, false, v);
                    dp[v] = 0.0;
                    dm[v] = 0.0;
                }
                match (
                    self.neighbour(e, axis, true),
                    self.neighbour(e, axis, false),
                ) {
                    | (Some(p), Some(mi)) => {
                        for v in 0..nv {
                            dp[v] = means[p * nv + v] - mean[v];
                            dm[v] = mean[v] - means[mi * nv + v];
                        }
                    },
                    | (Some(p), None) => {
                        // the missing neighbour is the ghost element of the boundary condition
                        let ghost = self.ghost_modal(e, axis, false, q, t);
                        for v in 0..nv {
                            dp[v] = means[p * nv + v] - mean[v];
                            dm[v] = ghost.as_ref().map_or(dp[v], |g| mean[v] - g[v * self.np]);
                        }
                    },
                    | (None, Some(mi)) => {
                        let ghost = self.ghost_modal(e, axis, true, q, t);
                        for v in 0..nv {
                            dm[v] = mean[v] - means[mi * nv + v];
                            dp[v] = ghost.as_ref().map_or(dm[v], |g| g[v * self.np] - mean[v]);
                        }
                    },
                    | (None, None) => {},
                }
                if has_char {
                    self.law.to_characteristic(mean, axis, &up, &mut tmp);
                    up.copy_from_slice(&tmp);
                    self.law.to_characteristic(mean, axis, &um, &mut tmp);
                    um.copy_from_slice(&tmp);
                    self.law.to_characteristic(mean, axis, &dp, &mut tmp);
                    dp.copy_from_slice(&tmp);
                    self.law.to_characteristic(mean, axis, &dm, &mut tmp);
                    dm.copy_from_slice(&tmp);
                }
                let thr = mt * m.h(axis) * m.h(axis);
                let scale = 1e-9
                    * (mean.iter().map(|x| x.abs()).sum::<f64>()
                        + up.iter().map(|x| x.abs()).sum::<f64>()
                        + um.iter().map(|x| x.abs()).sum::<f64>()
                        + 1e-300);
                for k in 0..nv {
                    let tvb = |x: f64| {
                        if x.abs() <= thr {
                            x
                        } else {
                            minmod3(x, dp[k], dm[k])
                        }
                    };
                    if (tvb(up[k]) - up[k]).abs() > scale || (tvb(um[k]) - um[k]).abs() > scale {
                        trouble_e = true;
                    }
                }
                if trouble_e {
                    break;
                }
            }
            trouble_e
        })
    }

    /// Characteristic / hierarchical moment limiter.
    ///
    /// The Legendre modes of each troubled element are bounded from the top
    /// total degree down with Krivodonova's minmod of the mode and the
    /// neighbour differences of the next-lower mode; with `characteristic` the
    /// coefficient vectors are transformed with the eigenvectors of the normal
    /// Jacobian at the element mean first. Modes are limited along each axis
    /// using the left and right neighbour of that axis (periodic wrap-around
    /// for periodic faces, the ghost element of `Ctx::ghost_modal` at a
    /// non-periodic boundary). The element
    /// average (mode zero) is never touched, so the scheme stays conservative.
    fn limit_moment(
        &self,
        q: &mut [f64],
        trouble: &[bool],
        characteristic: bool,
        t: f64,
    ) {
        if !trouble.iter().any(|&t| t) {
            return;
        }
        let (nv, np, n1) = (self.nv, self.np, self.n1);
        let m = &self.mesh;
        let p = n1 - 1;
        let ne = m.elements();
        let has_char = characteristic && self.law.has_characteristics();
        let stride = [1, n1, n1 * n1];
        // Modal coefficients of every element: a troubled element is limited
        // against the lower modes of its (possibly untroubled) neighbours.
        let mut modal = vec![0.0; ne * nv * np];
        for e in 0..ne {
            for v in 0..nv {
                let vals: Vec<f64> = (0..np).map(|nd| q[(e * np + nd) * nv + v]).collect();
                let c = self.to_modal(&vals);
                modal[(e * nv + v) * np..(e * nv + v + 1) * np].copy_from_slice(&c);
            }
        }
        let old = modal.clone();
        let means = self.means(q);
        // Ghost elements behind non-periodic faces of the troubled elements.
        let mut ghosts: std::collections::HashMap<(usize, usize, bool), Vec<f64>> =
            std::collections::HashMap::new();
        for e in (0..ne).filter(|&e| trouble[e]) {
            for axis in 0..3 {
                for max in [false, true] {
                    if self.neighbour(e, axis, max).is_none()
                        && let Some(g) = self.ghost_modal(e, axis, max, q, t)
                    {
                        ghosts.insert((e, axis, max), g);
                    }
                }
            }
        }
        let mut raw_deg = vec![0.0; nv];
        let mut raw_low = vec![0.0; nv];
        let mut raw_plus = vec![0.0; nv];
        let mut raw_minus = vec![0.0; nv];
        let mut w_deg = vec![0.0; nv];
        let mut w_low = vec![0.0; nv];
        let mut w_plus = vec![0.0; nv];
        let mut w_minus = vec![0.0; nv];
        let mut new_deg = vec![0.0; nv];
        let mut cons = vec![0.0; nv];
        for e in 0..ne {
            if !trouble[e] {
                continue;
            }
            let mean = &means[e * nv..(e + 1) * nv];
            let mut changed = false;
            for total in (1..=3 * p).rev() {
                for k in 0..=p {
                    for l in 0..=p {
                        for mm in 0..=p {
                            if k + l + mm != total {
                                continue;
                            }
                            let index = [k, l, mm];
                            for axis in 0..3 {
                                let order = index[axis];
                                if order == 0 {
                                    continue;
                                }
                                let has = |max: bool| {
                                    self.neighbour(e, axis, max).is_some()
                                        || ghosts.contains_key(&(e, axis, max))
                                };
                                if !(has(true) && has(false)) {
                                    continue;
                                }
                                let idx = k + n1 * (l + n1 * mm);
                                let lidx = idx - stride[axis];
                                let fetch = |max: bool, v: usize| -> f64 {
                                    match self.neighbour(e, axis, max) {
                                        | Some(en) => old[(en * nv + v) * np + lidx],
                                        | None => {
                                            ghosts
                                                .get(&(e, axis, max))
                                                .map_or(0.0, |g| g[v * np + lidx])
                                        },
                                    }
                                };
                                for v in 0..nv {
                                    raw_deg[v] = old[(e * nv + v) * np + idx];
                                    raw_low[v] = old[(e * nv + v) * np + lidx];
                                    raw_plus[v] = fetch(true, v);
                                    raw_minus[v] = fetch(false, v);
                                }
                                if has_char {
                                    self.law.to_characteristic(mean, axis, &raw_deg, &mut w_deg);
                                    self.law.to_characteristic(mean, axis, &raw_low, &mut w_low);
                                    self.law
                                        .to_characteristic(mean, axis, &raw_plus, &mut w_plus);
                                    self.law.to_characteristic(
                                        mean,
                                        axis,
                                        &raw_minus,
                                        &mut w_minus,
                                    );
                                } else {
                                    w_deg.copy_from_slice(&raw_deg);
                                    w_low.copy_from_slice(&raw_low);
                                    w_plus.copy_from_slice(&raw_plus);
                                    w_minus.copy_from_slice(&raw_minus);
                                }
                                let f = (2 * order - 1) as f64;
                                let mut any = false;
                                for j in 0..nv {
                                    let lim = minmod3(
                                        f * w_deg[j],
                                        w_plus[j] - w_low[j],
                                        w_low[j] - w_minus[j],
                                    ) / f;
                                    if (lim - w_deg[j]).abs()
                                        > 1e-12 * (w_deg[j].abs() + 1e-300) + 1e-14
                                    {
                                        new_deg[j] = lim;
                                        any = true;
                                    } else {
                                        new_deg[j] = w_deg[j];
                                    }
                                }
                                if !any {
                                    continue;
                                }
                                if has_char {
                                    self.law
                                        .from_characteristic(mean, axis, &new_deg, &mut cons);
                                } else {
                                    cons.copy_from_slice(&new_deg);
                                }
                                for v in 0..nv {
                                    modal[(e * nv + v) * np + idx] = cons[v];
                                }
                                changed = true;
                            }
                        }
                    }
                }
            }
            if changed {
                for v in 0..nv {
                    let vals = self.from_modal(&modal[(e * nv + v) * np..(e * nv + v + 1) * np]);
                    for nd in 0..np {
                        q[(e * np + nd) * nv + v] = vals[nd];
                    }
                }
            }
        }
    }

    /// Modal coefficients (`nv * np`, variable-major) of the ghost element behind the
    /// non-periodic face `(axis, max)` of element `e`: the mirror image with the law's wall state
    /// for walls, the exterior state evaluated at the nodes of the virtual element for fixed
    /// states, and the even mirror image for outflow.
    fn ghost_modal(
        &self,
        e: usize,
        axis: usize,
        max: bool,
        q: &[f64],
        t: f64,
    ) -> Option<Vec<f64>> {
        let (nv, np, n1) = (self.nv, self.np, self.n1);
        let m = &self.mesh;
        let own = |nd: usize| &q[(e * np + nd) * nv..(e * np + nd + 1) * nv];
        let ijk = |nd: usize| [nd % n1, (nd / n1) % n1, nd / (n1 * n1)];
        let mirror = |nd: usize| {
            let mut c = ijk(nd);
            c[axis] = n1 - 1 - c[axis];
            c[0] + n1 * (c[1] + n1 * c[2])
        };
        let mut vals = vec![0.0; np * nv];
        match &self.bcs[2 * axis + usize::from(max)] {
            | Dg3Bc::Periodic => return None,
            | Dg3Bc::Wall => {
                for nd in 0..np {
                    self.law
                        .wall_state(own(mirror(nd)), axis, &mut vals[nd * nv..(nd + 1) * nv]);
                }
            },
            | Dg3Bc::Fixed(f) => {
                let ex = [
                    e % m.ne[0],
                    (e / m.ne[0]) % m.ne[1],
                    e / (m.ne[0] * m.ne[1]),
                ];
                for nd in 0..np {
                    let c = ijk(nd);
                    let x: [f64; 3] = std::array::from_fn(|b| {
                        let cell = if b == axis {
                            if max {
                                m.ne[axis] as f64
                            } else {
                                -1.0
                            }
                        } else {
                            ex[b] as f64
                        };
                        m.lo[b] + (cell + f64::midpoint(1.0, self.xi[c[b]])) * m.h(b)
                    });
                    vals[nd * nv..(nd + 1) * nv].copy_from_slice(&f(t, x));
                }
            },
            | Dg3Bc::Outflow => {
                // The outflow flux uses the interior trace as exterior state, so the ghost
                // element is the even mirror image of the element.
                for nd in 0..np {
                    vals[nd * nv..(nd + 1) * nv].copy_from_slice(own(mirror(nd)));
                }
            },
        }
        let mut out = vec![0.0; nv * np];
        for v in 0..nv {
            let nodal: Vec<f64> = (0..np).map(|nd| vals[nd * nv + v]).collect();
            out[v * np..(v + 1) * np].copy_from_slice(&self.to_modal(&nodal));
        }
        Some(out)
    }

    /// WENO combination of `ncand` nodal candidate vectors of `n1` values,
    /// stored at `cands[k * stride ..]`, using the 1D Zhong–Shu oscillation
    /// indicator. `out` receives the reconstructed values; `shifts`, `d1`,
    /// `d2` and `omegas` are caller-provided scratch buffers.
    #[allow(clippy::too_many_arguments)] // small numerical kernel kept allocation-free
    fn weno_combine(
        &self,
        cands: &[f64],
        stride: usize,
        ncand: usize,
        weights: &[f64],
        mean_i: f64,
        out: &mut [f64],
        shifts: &mut [f64],
        d1: &mut [f64],
        d2: &mut [f64],
        omegas: &mut [f64],
    ) {
        let n1 = self.n1;
        let eps = self.weno_eps * (mean_i * mean_i + 1e-12);
        for k in 0..ncand {
            let vals = &cands[k * stride..k * stride + n1];
            let mc: f64 = (0..n1).map(|s| 0.5 * self.w[s] * vals[s]).sum();
            shifts[k] = mean_i - mc;
            for s in 0..n1 {
                d1[s] = vals[s] + shifts[k];
            }
            let mut beta = 0.0;
            for l in 1..=self.mesh.p {
                for i in 0..n1 {
                    d2[i] = (0..n1).map(|j| self.d[i][j] * d1[j]).sum();
                }
                let scale = f64::from(1u32 << (2 * l - 1));
                beta += scale * (0..n1).map(|j| self.w[j] * d2[j] * d2[j]).sum::<f64>();
                d1.copy_from_slice(d2);
            }
            omegas[k] = weights[k] / ((eps + beta) * (eps + beta));
        }
        let tot: f64 = (0..ncand).map(|k| omegas[k]).sum();
        for s in 0..n1 {
            out[s] = (0..ncand)
                .map(|k| omegas[k] / tot * (cands[k * stride + s] + shifts[k]))
                .sum();
        }
    }

    /// Zhong–Shu WENO limiter: a troubled element is reconstructed along every
    /// axis from its own polynomial and the shifted polynomials of its
    /// left/right neighbours, and the directional reconstructions are averaged
    /// (each preserves the element mean). The reconstruction is componentwise,
    /// as in the 1D/2D solver.
    fn limit_weno(
        &self,
        q: &mut [f64],
        time: f64,
    ) {
        let (nv, np, n1) = (self.nv, self.np, self.n1);
        let m = &self.mesh;
        let ne = m.elements();
        let src = q.to_vec();
        let trouble = self.troubled_elements(&src, self.tvd, false, time);
        if !trouble.iter().any(|&t| t) {
            return;
        }
        let mut cand = vec![0.0; 3 * nv * n1];
        let mut out_comp = vec![0.0; nv * n1];
        let mut acc = vec![0.0; np * nv];
        let mut weights = [0.0; 3];
        let mut shifts = [0.0; 3];
        let mut omegas = [0.0; 3];
        let mut d1 = vec![0.0; n1];
        let mut d2 = vec![0.0; n1];
        let mut nb_line = vec![0.0; n1];
        for e in 0..ne {
            if !trouble[e] {
                continue;
            }
            acc.fill(0.0);
            let mut axes_done: usize = 0;
            for axis in 0..3 {
                let mut axis_used = false;
                for t2 in 0..n1 {
                    for t1 in 0..n1 {
                        let mut ncand = 1;
                        weights[0] = 0.998;
                        for v in 0..nv {
                            for s in 0..n1 {
                                cand[v * n1 + s] =
                                    src[(e * np + self.axis_node(axis, s, t1, t2)) * nv + v];
                            }
                        }
                        for (upper, shift) in [(false, -2.0), (true, 2.0)] {
                            let Some(en) = self.neighbour(e, axis, upper) else {
                                continue;
                            };
                            for v in 0..nv {
                                for j in 0..n1 {
                                    nb_line[j] =
                                        src[(en * np + self.axis_node(axis, j, t1, t2)) * nv + v];
                                }
                                for s in 0..n1 {
                                    cand[(ncand * nv + v) * n1 + s] =
                                        lagrange_eval(&self.xi, &nb_line, self.xi[s] - shift);
                                }
                            }
                            weights[ncand] = 0.001;
                            ncand += 1;
                        }
                        if ncand == 1 {
                            continue;
                        }
                        axis_used = true;
                        if ncand == 2 {
                            weights[0] = 0.999;
                        }
                        for v in 0..nv {
                            let mean_v: f64 =
                                (0..n1).map(|s| 0.5 * self.w[s] * cand[v * n1 + s]).sum();
                            self.weno_combine(
                                &cand[v * n1..],
                                nv * n1,
                                ncand,
                                &weights,
                                mean_v,
                                &mut out_comp[v * n1..],
                                &mut shifts,
                                &mut d1,
                                &mut d2,
                                &mut omegas,
                            );
                        }
                        for s in 0..n1 {
                            let nd = self.axis_node(axis, s, t1, t2);
                            for v in 0..nv {
                                acc[nd * nv + v] += out_comp[v * n1 + s];
                            }
                        }
                    }
                }
                if axis_used {
                    axes_done += 1;
                }
            }
            if axes_done > 0 {
                for k in 0..np * nv {
                    q[e * np * nv + k] = acc[k] / axes_done as f64;
                }
            }
        }
    }

    /// Applies the configured shock limiter and then the Zhang–Shu positivity scaling.
    fn limit_all(
        &self,
        q: &mut [f64],
        t: f64,
    ) {
        match self.limiter {
            | Dg3Limiter::None => {},
            | Dg3Limiter::Characteristic => {
                let trouble = self.troubled_elements(q, self.tvd, true, t);
                self.limit_moment(q, &trouble, true, t);
            },
            | Dg3Limiter::Hierarchical => {
                let trouble = self.troubled_elements(q, self.tvd, false, t);
                self.limit_moment(q, &trouble, false, t);
            },
            | Dg3Limiter::Weno => self.limit_weno(q, t),
        }
        if self.positivity {
            self.positivity(q);
        }
    }
}

/// Integrates `q_t + ∇·F(q) = 0` on a hexahedral mesh to `t_end` with SSP-RK3.
///
/// `bcs` lists `[XMin, XMax, YMin, YMax, ZMin, ZMax]`.
///
/// # Errors
/// For an invalid mesh, inconsistent periodic faces, wrong initial data size
/// or a non-finite wave speed.
pub fn solve_dg3d(
    mesh: &Dg3Mesh,
    law: &dyn DgLaw,
    bcs: &[Dg3Bc; 6],
    init: &dyn Fn([f64; 3]) -> Vec<f64>,
    t_end: f64,
    cfl: f64,
    opt: &Dg3Options,
) -> Result<Dg3Solution, String> {
    if mesh.p == 0 || mesh.ne.contains(&0) {
        return Err("invalid DG mesh".into());
    }
    for a in 0..3 {
        if matches!(bcs[2 * a], Dg3Bc::Periodic) != matches!(bcs[2 * a + 1], Dg3Bc::Periodic) {
            return Err("periodic must be set on both faces of an axis".into());
        }
    }
    let (xi, w) = lgl(mesh.p);
    let d = lgl_diff_matrix(&xi);
    let nv = law.nvar();
    let n1 = mesh.p + 1;
    let np = mesh.nodes_per_element();
    let vmat: Vec<Vec<f64>> = xi.iter().map(|&x| legendre(mesh.p, x)).collect();
    let vinv = invert(vmat.clone());
    let ctx = Ctx {
        mesh: *mesh,
        law,
        bcs,
        d,
        w: w.clone(),
        xi: xi.clone(),
        nv,
        np,
        n1,
        es: opt.entropy_stable,
        limiter: opt.limiter,
        tvd: opt.limiter_tvd,
        weno_eps: opt.limiter_weno_eps,
        positivity: opt.positivity,
        vmat,
        vinv,
    };
    let mut sol = Dg3Solution {
        mesh: *mesh,
        xi,
        w,
        q: vec![0.0; mesh.elements() * np * nv],
        nvar: nv,
        t: 0.0,
        steps: 0,
    };
    for e in 0..mesh.elements() {
        for nd in 0..np {
            let v = init(sol.coord(e, nd));
            if v.len() != nv {
                return Err("initial data has the wrong number of components".into());
            }
            sol.q[(e * np + nd) * nv..(e * np + nd + 1) * nv].copy_from_slice(&v);
        }
    }
    let limiting = opt.limiter != Dg3Limiter::None || opt.positivity;
    let mut q = sol.q.clone();
    let (mut t, mut steps) = (0.0, 0);
    while t < t_end * (1.0 - 1e-14) {
        let rate = ctx.max_rate(&q);
        if !rate.is_finite() || rate <= 0.0 {
            return Err("non-finite wave speed".into());
        }
        let dt = (cfl / ((2 * mesh.p + 1) as f64 * rate)).min(t_end - t);
        let l0 = ctx.rhs(&q, t);
        let mut q1: Vec<f64> = q.iter().zip(&l0).map(|(a, b)| a + dt * b).collect();
        if limiting {
            ctx.limit_all(&mut q1, t + dt);
        }
        let l1 = ctx.rhs(&q1, t + dt);
        let mut q2: Vec<f64> = (0..q.len())
            .map(|i| 0.75 * q[i] + 0.25 * (q1[i] + dt * l1[i]))
            .collect();
        if limiting {
            ctx.limit_all(&mut q2, t + 0.5 * dt);
        }
        let l2 = ctx.rhs(&q2, t + 0.5 * dt);
        q = (0..q.len())
            .map(|i| q[i] / 3.0 + 2.0 / 3.0 * (q2[i] + dt * l2[i]))
            .collect();
        if limiting {
            ctx.limit_all(&mut q, t + dt);
        }
        t += dt;
        steps += 1;
        if steps > 5_000_000 {
            return Err("too many steps".into());
        }
    }
    sol.q = q;
    sol.t = t;
    sol.steps = steps;
    Ok(sol)
}
