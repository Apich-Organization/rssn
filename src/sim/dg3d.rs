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
//! * Zhang–Shu positivity-preserving scaling of density and pressure.
//!
//! Smooth problems converge with order `p + 1` (advection) and `p`–`p + 1`
//! (Euler vortex). Limitations: elements are axis-aligned boxes of one
//! uniform size (no curved or unstructured hexahedra, no hanging nodes), the
//! characteristic / hierarchical / WENO limiters of the 1D/2D solver are not
//! available in 3D.

use crate::sim::dg::{DgLaw, euler_ec_flux, euler_prim, euler_wall_state, lgl, lgl_diff_matrix};
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
    fn flux(&self, q: &[f64], dir: usize, out: &mut [f64]) {
        out[0] = self.a[dir] * q[0];
    }
    fn max_speed(&self, _q: &[f64], dir: usize) -> f64 {
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
    pub fn pressure(&self, q: &[f64]) -> f64 {
        (self.gamma - 1.0) * (q[4] - 0.5 * (q[1] * q[1] + q[2] * q[2] + q[3] * q[3]) / q[0])
    }
}

impl DgLaw for Euler3d {
    fn nvar(&self) -> usize {
        5
    }
    fn flux(&self, q: &[f64], dir: usize, out: &mut [f64]) {
        let p = self.pressure(q);
        let vel = q[1 + dir] / q[0];
        out[0] = q[1 + dir];
        for k in 0..3 {
            out[1 + k] = q[1 + k] * vel;
        }
        out[1 + dir] += p;
        out[4] = (q[4] + p) * vel;
    }
    fn max_speed(&self, q: &[f64], dir: usize) -> f64 {
        let p = self.pressure(q).max(0.0);
        (q[1 + dir] / q[0]).abs() + (self.gamma * p / q[0]).sqrt()
    }
    fn rho_p(&self, q: &[f64]) -> Option<(f64, f64)> {
        Some((q[0], self.pressure(q)))
    }
    fn two_point_flux(&self, ql: &[f64], qr: &[f64], dir: usize, out: &mut [f64]) {
        euler_ec_flux(self.gamma, 3, dir, ql, qr, out);
    }
    fn wall_state(&self, q: &[f64], dir: usize, out: &mut [f64]) {
        euler_wall_state(3, dir, q, out);
    }
}

/// Specific thermodynamic entropy `s = ln(p / ρ^γ)` of an Euler state with `nvel` velocity components.
#[must_use]
pub fn euler_entropy(gamma: f64, nvel: usize, q: &[f64]) -> f64 {
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
    fn h(&self, a: usize) -> f64 {
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

/// Options of [`solve_dg3d`].
#[derive(Clone, Copy, Debug, Default)]
pub struct Dg3Options {
    /// Entropy-stable flux differencing (needs [`DgLaw::two_point_flux`]).
    pub entropy_stable: bool,
    /// Positivity-preserving scaling of density and pressure.
    pub positivity: bool,
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
    const fn node_ijk(&self, node: usize) -> [usize; 3] {
        let n1 = self.mesh.p + 1;
        [node % n1, (node / n1) % n1, node / (n1 * n1)]
    }

    const fn elem_ijk(&self, e: usize) -> [usize; 3] {
        let ne = self.mesh.ne;
        [e % ne[0], (e / ne[0]) % ne[1], e / (ne[0] * ne[1])]
    }

    /// Physical coordinates of a node.
    #[must_use]
    pub fn coord(&self, e: usize, node: usize) -> [f64; 3] {
        let m = &self.mesh;
        let (ij, ex) = (self.node_ijk(node), self.elem_ijk(e));
        std::array::from_fn(|a| m.lo[a] + (ex[a] as f64 + f64::midpoint(1.0, self.xi[ij[a]])) * m.h(a))
    }

    /// Value of a variable at a node.
    #[must_use]
    pub fn value(&self, e: usize, node: usize, var: usize) -> f64 {
        self.q[(e * self.mesh.nodes_per_element() + node) * self.nvar + var]
    }

    fn quad_weight(&self, node: usize) -> f64 {
        let ij = self.node_ijk(node);
        let jac: f64 = (0..3).map(|a| 0.5 * self.mesh.h(a)).product();
        jac * (0..3).map(|a| self.w[ij[a]]).product::<f64>()
    }

    /// Quadrature `L2` norm of `q_var - exact`.
    pub fn l2_error(&self, var: usize, exact: impl Fn([f64; 3]) -> f64) -> f64 {
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
    pub fn integral(&self, f: impl Fn(&[f64]) -> f64) -> f64 {
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
}

impl Ctx<'_> {
    const fn neighbour(&self, e: usize, axis: usize, max: bool) -> Option<usize> {
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

    fn rhs(&self, q: &[f64], t: f64) -> Vec<f64> {
        let m = &self.mesh;
        let (nv, np, n1) = (self.nv, self.np, self.n1);
        let mut r = vec![0.0; q.len()];
        let mut fbuf = vec![0.0; nv];
        let mut flux = vec![0.0; np * nv];
        let at = |e: usize, nd: usize| -> &[f64] {
            let o = (e * np + nd) * nv;
            &q[o..o + nv]
        };
        let node = |i: usize, j: usize, k: usize| i + n1 * (j + n1 * k);
        for e in 0..m.elements() {
            for a in 0..3 {
                let sc = 2.0 / m.h(a);
                if !self.es {
                    for nd in 0..np {
                        self.law.flux(at(e, nd), a, &mut fbuf);
                        flux[nd * nv..(nd + 1) * nv].copy_from_slice(&fbuf);
                    }
                }
                // index of the line node `s` for transverse indices (t1, t2)
                let idx = |s: usize, t1: usize, t2: usize| match a {
                    0 => node(s, t1, t2),
                    1 => node(t1, s, t2),
                    _ => node(t1, t2, s),
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
                                        r[(e * np + nd) * nv + c] -= 2.0 * sc * self.d[i][k] * fbuf[c];
                                    }
                                } else {
                                    for c in 0..nv {
                                        r[(e * np + nd) * nv + c] -= sc * self.d[i][k] * flux[nk * nv + c];
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
                                Some(en) => qo.copy_from_slice(at(en, idx(n1 - 1 - edge, t1, t2))),
                                None => match &self.bcs[2 * a + usize::from(max)] {
                                    Dg3Bc::Fixed(f) => {
                                        let ex = [e % m.ne[0], (e / m.ne[0]) % m.ne[1], e / (m.ne[0] * m.ne[1])];
                                        let ij = match a {
                                            0 => [edge, t1, t2],
                                            1 => [t1, edge, t2],
                                            _ => [t1, t2, edge],
                                        };
                                        let x: [f64; 3] =
                                            std::array::from_fn(|b| m.lo[b] + (ex[b] as f64 + f64::midpoint(1.0, self.xi[ij[b]])) * m.h(b));
                                        qo.copy_from_slice(&f(t, x));
                                    },
                                    Dg3Bc::Wall => self.law.wall_state(qi, a, &mut qo),
                                    _ => qo.copy_from_slice(qi),
                                },
                            }
                            let mut fstar = vec![0.0; nv];
                            match (max, self.es) {
                                (true, false) => self.law.numerical_flux(qi, &qo, a, &mut fstar),
                                (false, false) => self.law.numerical_flux(&qo, qi, a, &mut fstar),
                                (true, true) => self.law.entropy_stable_flux(qi, &qo, a, &mut fstar),
                                (false, true) => self.law.entropy_stable_flux(&qo, qi, a, &mut fstar),
                            }
                            self.law.flux(qi, a, &mut fbuf);
                            let sign = if max { -1.0 } else { 1.0 };
                            for c in 0..nv {
                                r[(e * np + nd) * nv + c] += sign * sc * (fstar[c] - fbuf[c]) / self.w[edge];
                            }
                        }
                    }
                }
            }
        }
        r
    }

    fn max_rate(&self, q: &[f64]) -> f64 {
        q.chunks(self.nv)
            .map(|c| (0..3).map(|a| self.law.max_speed(c, a) / self.mesh.h(a)).sum::<f64>())
            .fold(0.0, f64::max)
    }

    fn node_weight(&self, nd: usize) -> f64 {
        let n1 = self.n1;
        let ij = [nd % n1, (nd / n1) % n1, nd / (n1 * n1)];
        (0..3).map(|a| self.w[ij[a]] / 2.0).product()
    }

    /// Zhang–Shu scaling of density and pressure towards the element mean.
    fn positivity(&self, q: &mut [f64]) {
        let (nv, np) = (self.nv, self.np);
        let eps = 1e-10;
        for e in 0..self.mesh.elements() {
            let mean: Vec<f64> = (0..nv).map(|v| (0..np).map(|nd| self.node_weight(nd) * q[(e * np + nd) * nv + v]).sum()).collect();
            let Some((rm, pm)) = self.law.rho_p(&mean) else { return };
            if rm <= eps || pm <= eps {
                for nd in 0..np {
                    for v in 0..nv {
                        q[(e * np + nd) * nv + v] = mean[v];
                    }
                }
                continue;
            }
            let rho_min = (0..np).map(|nd| q[(e * np + nd) * nv]).fold(f64::INFINITY, f64::min);
            let th1 = if rho_min < eps { ((rm - eps) / (rm - rho_min)).clamp(0.0, 1.0) } else { 1.0 };
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
                    self.law.rho_p(&st).is_some_and(|(r, p)| r >= eps && p >= eps)
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
    let ctx = Ctx { mesh: *mesh, law, bcs, d, w: w.clone(), xi: xi.clone(), nv, np, n1, es: opt.entropy_stable };
    let mut sol = Dg3Solution { mesh: *mesh, xi, w, q: vec![0.0; mesh.elements() * np * nv], nvar: nv, t: 0.0, steps: 0 };
    for e in 0..mesh.elements() {
        for nd in 0..np {
            let v = init(sol.coord(e, nd));
            if v.len() != nv {
                return Err("initial data has the wrong number of components".into());
            }
            sol.q[(e * np + nd) * nv..(e * np + nd + 1) * nv].copy_from_slice(&v);
        }
    }
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
        if opt.positivity {
            ctx.positivity(&mut q1);
        }
        let l1 = ctx.rhs(&q1, t + dt);
        let mut q2: Vec<f64> = (0..q.len()).map(|i| 0.75 * q[i] + 0.25 * (q1[i] + dt * l1[i])).collect();
        if opt.positivity {
            ctx.positivity(&mut q2);
        }
        let l2 = ctx.rhs(&q2, t + 0.5 * dt);
        q = (0..q.len()).map(|i| q[i] / 3.0 + 2.0 / 3.0 * (q2[i] + dt * l2[i])).collect();
        if opt.positivity {
            ctx.positivity(&mut q);
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
