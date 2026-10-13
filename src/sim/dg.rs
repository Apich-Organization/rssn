//! # Nodal discontinuous Galerkin for conservation laws (1D and 2D)
//!
//! Collocation (DGSEM) discontinuous Galerkin on Cartesian element meshes:
//! Legendre–Gauss–Lobatto nodes of any polynomial degree `p >= 1`, strong
//! form with surface lifting, local Lax–Friedrichs (Rusanov) interface
//! fluxes, SSP-RK3 time stepping and an optional exponential modal filter
//! for shock capturing.
//!
//! | Law | 1D | 2D | Boundary conditions |
//! |---|---|---|---|
//! | linear advection ([`LinearAdvection`]) | yes | yes | periodic, outflow, fixed inflow state |
//! | Burgers ([`Burgers`]) | yes | yes | same |
//! | Euler ([`Euler1d`]) | yes | n/a | same |
//!
//! Custom systems implement [`DgLaw`]. Meshes are uniform tensor-product
//! quadrilaterals (intervals in 1D); no hanging nodes, no p-adaptivity.
//! Smooth problems converge with order `p + 1`; shocks need the filter and
//! are first-order accurate locally.
#![allow(
    clippy::manual_memcpy,
    clippy::too_long_first_doc_paragraph,
    clippy::cast_lossless,
    clippy::single_match_else,
    clippy::manual_contains,
    clippy::too_many_lines,
    clippy::too_many_arguments,
    clippy::needless_range_loop,
    clippy::missing_const_for_fn,
    clippy::float_cmp,
    clippy::type_complexity,
    clippy::many_single_char_names,
    clippy::manual_midpoint,
    clippy::while_float
)]

use std::sync::Arc;

/// A system of conservation laws `q_t + ∇·F(q) = 0`.
pub trait DgLaw: Send + Sync {
    /// Number of conserved variables.
    fn nvar(&self) -> usize;
    /// Flux `F_dir(q)` written to `out`.
    fn flux(&self, q: &[f64], dir: usize, out: &mut [f64]);
    /// Largest signal speed in direction `dir`.
    fn max_speed(&self, q: &[f64], dir: usize) -> f64;
    /// Density and pressure of a state for positivity limiting (`None` for laws without such notion).
    fn rho_p(&self, _q: &[f64]) -> Option<(f64, f64)> {
        None
    }
    /// Numerical flux between a left and right state (default: Rusanov).
    fn numerical_flux(&self, ql: &[f64], qr: &[f64], dir: usize, out: &mut [f64]) {
        let n = self.nvar();
        let mut fl = vec![0.0; n];
        let mut fr = vec![0.0; n];
        self.flux(ql, dir, &mut fl);
        self.flux(qr, dir, &mut fr);
        let s = self.max_speed(ql, dir).max(self.max_speed(qr, dir));
        for c in 0..n {
            out[c] = 0.5 * (fl[c] + fr[c]) - 0.5 * s * (qr[c] - ql[c]);
        }
    }
}

/// Linear advection `u_t + a·∇u = 0`.
#[derive(Clone, Copy, Debug)]
pub struct LinearAdvection {
    /// Advection velocity (second entry used in 2D).
    pub a: [f64; 2],
}

impl DgLaw for LinearAdvection {
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

/// Burgers equation `u_t + ∇·(u²/2 (1, 1)) = 0`.
#[derive(Clone, Copy, Debug)]
pub struct Burgers;

impl DgLaw for Burgers {
    fn nvar(&self) -> usize {
        1
    }
    fn flux(&self, q: &[f64], _dir: usize, out: &mut [f64]) {
        out[0] = 0.5 * q[0] * q[0];
    }
    fn max_speed(&self, q: &[f64], _dir: usize) -> f64 {
        q[0].abs()
    }
}

/// One-dimensional Euler equations `(ρ, ρu, E)`.
#[derive(Clone, Copy, Debug)]
pub struct Euler1d {
    /// Adiabatic exponent.
    pub gamma: f64,
}

impl DgLaw for Euler1d {
    fn nvar(&self) -> usize {
        3
    }
    fn flux(&self, q: &[f64], _dir: usize, out: &mut [f64]) {
        let u = q[1] / q[0];
        let p = (self.gamma - 1.0) * (q[2] - 0.5 * q[1] * u);
        out[0] = q[1];
        out[1] = q[1] * u + p;
        out[2] = (q[2] + p) * u;
    }
    fn max_speed(&self, q: &[f64], _dir: usize) -> f64 {
        let u = q[1] / q[0];
        let p = (self.gamma - 1.0) * (q[2] - 0.5 * q[1] * u);
        u.abs() + (self.gamma * p.max(0.0) / q[0]).sqrt()
    }
    fn rho_p(&self, q: &[f64]) -> Option<(f64, f64)> {
        Some((q[0], self.pressure(q)))
    }
}

impl Euler1d {
    fn pressure(self, q: &[f64]) -> f64 {
        (self.gamma - 1.0) * (q[2] - 0.5 * q[1] * q[1] / q[0])
    }
}

/// Two-dimensional Euler equations `(ρ, ρu, ρv, E)`.
#[derive(Clone, Copy, Debug)]
pub struct Euler2d {
    /// Adiabatic exponent.
    pub gamma: f64,
}

impl Euler2d {
    fn pressure(self, q: &[f64]) -> f64 {
        (self.gamma - 1.0) * (q[3] - 0.5 * (q[1] * q[1] + q[2] * q[2]) / q[0])
    }
}

impl DgLaw for Euler2d {
    fn nvar(&self) -> usize {
        4
    }
    fn flux(&self, q: &[f64], dir: usize, out: &mut [f64]) {
        let p = self.pressure(q);
        let vel = q[1 + dir] / q[0];
        out[0] = q[1 + dir];
        out[1] = q[1] * vel;
        out[2] = q[2] * vel;
        out[1 + dir] += p;
        out[3] = (q[3] + p) * vel;
    }
    fn max_speed(&self, q: &[f64], dir: usize) -> f64 {
        let p = self.pressure(q).max(0.0);
        (q[1 + dir] / q[0]).abs() + (self.gamma * p / q[0]).sqrt()
    }
    fn rho_p(&self, q: &[f64]) -> Option<(f64, f64)> {
        Some((q[0], self.pressure(q)))
    }
}

/// Boundary treatment of one mesh face.
#[derive(Clone)]
pub enum DgBc {
    /// Periodic (set on both faces of the axis).
    Periodic,
    /// Zero-gradient outflow.
    Outflow,
    /// Fixed exterior state `q(t, x)` (inflow).
    Fixed(Arc<dyn Fn(f64, [f64; 2]) -> Vec<f64> + Send + Sync>),
}

/// Uniform Cartesian element mesh.
#[derive(Clone, Copy, Debug)]
pub struct DgMesh {
    /// Dimension (1 or 2).
    pub dim: usize,
    /// Elements per axis.
    pub ne: [usize; 2],
    /// Lower corner.
    pub lo: [f64; 2],
    /// Upper corner.
    pub hi: [f64; 2],
    /// Polynomial degree.
    pub p: usize,
}

impl DgMesh {
    /// One-dimensional mesh of `ne` elements on `[lo, hi]`.
    #[must_use]
    pub fn line(ne: usize, lo: f64, hi: f64, p: usize) -> Self {
        Self { dim: 1, ne: [ne, 1], lo: [lo, 0.0], hi: [hi, 1.0], p }
    }

    /// Two-dimensional mesh.
    #[must_use]
    pub fn rect(ne: [usize; 2], lo: [f64; 2], hi: [f64; 2], p: usize) -> Self {
        Self { dim: 2, ne, lo, hi, p }
    }

    fn h(&self, a: usize) -> f64 {
        (self.hi[a] - self.lo[a]) / self.ne[a] as f64
    }

    /// Nodes per element.
    #[must_use]
    pub fn nodes_per_element(&self) -> usize {
        (self.p + 1).pow(self.dim as u32)
    }

    /// Number of elements.
    #[must_use]
    pub fn elements(&self) -> usize {
        self.ne[0] * if self.dim == 2 { self.ne[1] } else { 1 }
    }
}

/// Legendre polynomials `P_0..=P_n` at `x`.
fn legendre(n: usize, x: f64) -> Vec<f64> {
    let mut v = vec![1.0, x];
    for k in 1..n {
        let kf = k as f64;
        v.push(((2.0 * kf + 1.0) * x * v[k] - kf * v[k - 1]) / (kf + 1.0));
    }
    v.truncate(n + 1);
    v
}

/// Legendre–Gauss–Lobatto nodes and weights with `p + 1` points.
#[must_use]
pub fn lgl(p: usize) -> (Vec<f64>, Vec<f64>) {
    let n = p;
    let mut x: Vec<f64> = (0..=n).map(|i| -(std::f64::consts::PI * i as f64 / n as f64).cos()).collect();
    for i in 1..n {
        for _ in 0..100 {
            let l = legendre(n, x[i]);
            let dx = (x[i] * l[n] - l[n - 1]) / ((n as f64 + 1.0) * l[n]);
            x[i] -= dx;
            if dx.abs() < 1e-15 {
                break;
            }
        }
    }
    x[0] = -1.0;
    x[n] = 1.0;
    let w: Vec<f64> = x
        .iter()
        .map(|&xi| {
            let l = legendre(n, xi);
            2.0 / (n as f64 * (n as f64 + 1.0) * l[n] * l[n])
        })
        .collect();
    (x, w)
}

/// Differentiation matrix on the LGL nodes.
#[must_use]
pub fn lgl_diff_matrix(x: &[f64]) -> Vec<Vec<f64>> {
    let n = x.len() - 1;
    let l: Vec<f64> = x.iter().map(|&xi| legendre(n, xi)[n]).collect();
    let mut d = vec![vec![0.0; n + 1]; n + 1];
    for i in 0..=n {
        for j in 0..=n {
            d[i][j] = if i != j {
                l[i] / (l[j] * (x[i] - x[j]))
            } else if i == 0 {
                -(n as f64) * (n as f64 + 1.0) / 4.0
            } else if i == n {
                (n as f64) * (n as f64 + 1.0) / 4.0
            } else {
                0.0
            };
        }
    }
    d
}

/// Exponential modal filter `σ_k = exp(-α (k/p)^order)`.
#[derive(Clone, Copy, Debug)]
pub struct ModalFilter {
    /// Strength (`36` is the machine-precision damping of the top mode).
    pub alpha: f64,
    /// Even order of the filter.
    pub order: u32,
}

fn invert(mut a: Vec<Vec<f64>>) -> Vec<Vec<f64>> {
    let n = a.len();
    let mut inv: Vec<Vec<f64>> = (0..n).map(|i| (0..n).map(|j| f64::from(u8::from(i == j))).collect()).collect();
    for c in 0..n {
        let piv = (c..n).max_by(|&i, &j| a[i][c].abs().total_cmp(&a[j][c].abs())).unwrap_or(c);
        a.swap(c, piv);
        inv.swap(c, piv);
        let d = a[c][c];
        for j in 0..n {
            a[c][j] /= d;
            inv[c][j] /= d;
        }
        for r in 0..n {
            if r != c {
                let f = a[r][c];
                for j in 0..n {
                    a[r][j] -= f * a[c][j];
                    inv[r][j] -= f * inv[c][j];
                }
            }
        }
    }
    inv
}

/// Solution snapshot.
#[derive(Clone, Debug)]
pub struct DgSolution {
    /// Mesh.
    pub mesh: DgMesh,
    /// Reference nodes on `[-1, 1]`.
    pub xi: Vec<f64>,
    /// Quadrature weights on `[-1, 1]`.
    pub w: Vec<f64>,
    /// Values, layout `[element][node][variable]`.
    pub q: Vec<f64>,
    /// Variables.
    pub nvar: usize,
    /// Time reached.
    pub t: f64,
    /// Steps taken.
    pub steps: usize,
}

impl DgSolution {
    fn node_ij(&self, node: usize) -> [usize; 2] {
        let n1 = self.mesh.p + 1;
        [node % n1, node / n1]
    }

    /// Physical coordinates of a node.
    #[must_use]
    pub fn coord(&self, e: usize, node: usize) -> [f64; 2] {
        let m = &self.mesh;
        let ij = self.node_ij(node);
        let ex = [e % m.ne[0], e / m.ne[0]];
        let mut x = [0.0; 2];
        for a in 0..m.dim {
            x[a] = m.lo[a] + (ex[a] as f64 + 0.5 * (1.0 + self.xi[ij[a]])) * m.h(a);
        }
        x
    }

    /// Value of a variable.
    #[must_use]
    pub fn value(&self, e: usize, node: usize, var: usize) -> f64 {
        self.q[(e * self.mesh.nodes_per_element() + node) * self.nvar + var]
    }

    /// Quadrature `L2` norm of `q_var - exact`.
    pub fn l2_error(&self, var: usize, exact: impl Fn([f64; 2]) -> f64) -> f64 {
        let m = &self.mesh;
        let jac: f64 = (0..m.dim).map(|a| 0.5 * m.h(a)).product();
        let mut s = 0.0;
        for e in 0..m.elements() {
            for nd in 0..m.nodes_per_element() {
                let ij = self.node_ij(nd);
                let w: f64 = (0..m.dim).map(|a| self.w[ij[a]]).product();
                let d = self.value(e, nd, var) - exact(self.coord(e, nd));
                s += w * jac * d * d;
            }
        }
        s.sqrt()
    }

    /// Integral of a variable over the domain.
    #[must_use]
    pub fn integral(&self, var: usize) -> f64 {
        let m = &self.mesh;
        let jac: f64 = (0..m.dim).map(|a| 0.5 * m.h(a)).product();
        let mut s = 0.0;
        for e in 0..m.elements() {
            for nd in 0..m.nodes_per_element() {
                let ij = self.node_ij(nd);
                let w: f64 = (0..m.dim).map(|a| self.w[ij[a]]).product();
                s += w * jac * self.value(e, nd, var);
            }
        }
        s
    }
}

struct Ctx<'a> {
    mesh: DgMesh,
    law: &'a dyn DgLaw,
    bcs: &'a [DgBc; 4],
    d: Vec<Vec<f64>>,
    w: Vec<f64>,
    nv: usize,
    np: usize,
    n1: usize,
}

impl Ctx<'_> {
    fn node(&self, i: usize, j: usize) -> usize {
        i + self.n1 * j
    }

    fn rhs(&self, q: &[f64], t: f64, xi: &[f64]) -> Vec<f64> {
        let m = &self.mesh;
        let (nv, np, n1) = (self.nv, self.np, self.n1);
        let ne = m.elements();
        let mut r = vec![0.0; q.len()];
        let mut fbuf = vec![0.0; nv];
        let mut flux = vec![0.0; np * nv];
        let at = |e: usize, nd: usize| -> &[f64] {
            let o = (e * np + nd) * nv;
            &q[o..o + nv]
        };
        let nbr = |e: usize, axis: usize, max: bool| -> Option<usize> {
            let ex = [e % m.ne[0], e / m.ne[0]];
            let mut c = ex;
            let periodic = matches!(self.bcs[2 * axis], DgBc::Periodic);
            if max {
                if c[axis] + 1 < m.ne[axis] {
                    c[axis] += 1;
                } else if periodic {
                    c[axis] = 0;
                } else {
                    return None;
                }
            } else if c[axis] > 0 {
                c[axis] -= 1;
            } else if periodic {
                c[axis] = m.ne[axis] - 1;
            } else {
                return None;
            }
            Some(c[0] + m.ne[0] * c[1])
        };
        for e in 0..ne {
            for a in 0..m.dim {
                let sc = 2.0 / m.h(a);
                for nd in 0..np {
                    self.law.flux(at(e, nd), a, &mut fbuf);
                    flux[nd * nv..(nd + 1) * nv].copy_from_slice(&fbuf);
                }
                // volume term
                let tr = if m.dim == 2 { n1 } else { 1 };
                for t_idx in 0..tr {
                    for i in 0..n1 {
                        let nd = if a == 0 { self.node(i, t_idx) } else { self.node(t_idx, i) };
                        for k in 0..n1 {
                            let nk = if a == 0 { self.node(k, t_idx) } else { self.node(t_idx, k) };
                            for c in 0..nv {
                                r[(e * np + nd) * nv + c] -= sc * self.d[i][k] * flux[nk * nv + c];
                            }
                        }
                    }
                    // surface terms
                    for max in [false, true] {
                        let edge = if max { self.n1 - 1 } else { 0 };
                        let nd = if a == 0 { self.node(edge, t_idx) } else { self.node(t_idx, edge) };
                        let qi = at(e, nd);
                        let nb = nbr(e, a, max);
                        let mut qo = vec![0.0; nv];
                        match nb {
                            Some(en) => {
                                let ond = if a == 0 { self.node(self.n1 - 1 - edge, t_idx) } else { self.node(t_idx, self.n1 - 1 - edge) };
                                qo.copy_from_slice(at(en, ond));
                            },
                            None => match &self.bcs[2 * a + usize::from(max)] {
                                DgBc::Fixed(f) => {
                                    let x = {
                                        let ij = if a == 0 { [edge, t_idx] } else { [t_idx, edge] };
                                        let ex = [e % m.ne[0], e / m.ne[0]];
                                        let mut x = [0.0; 2];
                                        for b in 0..m.dim {
                                            x[b] = m.lo[b] + (ex[b] as f64 + 0.5 * (1.0 + xi[ij[b]])) * m.h(b);
                                        }
                                        x
                                    };
                                    qo.copy_from_slice(&f(t, x));
                                },
                                _ => qo.copy_from_slice(qi),
                            },
                        }
                        let mut fstar = vec![0.0; nv];
                        if max {
                            self.law.numerical_flux(qi, &qo, a, &mut fstar);
                        } else {
                            self.law.numerical_flux(&qo, qi, a, &mut fstar);
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
        r
    }

    fn max_rate(&self, q: &[f64]) -> f64 {
        let mut best: f64 = 0.0;
        for chunk in q.chunks(self.nv) {
            let s: f64 = (0..self.mesh.dim).map(|a| self.law.max_speed(chunk, a) / self.mesh.h(a)).sum();
            best = best.max(s);
        }
        best
    }
}

/// Slope and positivity limiters applied after every Runge–Kutta stage.
#[derive(Clone, Copy, Debug)]
pub struct DgLimiter {
    /// TVB constant `M` of the minmod slope limiter (`None` disables it);
    /// `0` is the strict TVD minmod limiter, larger values keep smooth extrema.
    pub tvb_m: Option<f64>,
    /// Zhang–Shu positivity-preserving scaling of density and pressure
    /// (needs a law implementing [`DgLaw::rho_p`]).
    pub positivity: bool,
}

fn minmod3(a: f64, b: f64, c: f64) -> f64 {
    if a > 0.0 && b > 0.0 && c > 0.0 {
        a.min(b).min(c)
    } else if a < 0.0 && b < 0.0 && c < 0.0 {
        a.max(b).max(c)
    } else {
        0.0
    }
}

impl Ctx<'_> {
    fn neighbour(&self, e: usize, axis: usize, max: bool) -> Option<usize> {
        let m = &self.mesh;
        let mut c = [e % m.ne[0], e / m.ne[0]];
        let periodic = matches!(self.bcs[2 * axis], DgBc::Periodic);
        if max {
            if c[axis] + 1 < m.ne[axis] {
                c[axis] += 1;
            } else if periodic {
                c[axis] = 0;
            } else {
                return None;
            }
        } else if c[axis] > 0 {
            c[axis] -= 1;
        } else if periodic {
            c[axis] = m.ne[axis] - 1;
        } else {
            return None;
        }
        Some(c[0] + m.ne[0] * c[1])
    }

    fn node_weight(&self, nd: usize) -> f64 {
        let n1 = self.n1;
        let ij = [nd % n1, nd / n1];
        (0..self.mesh.dim).map(|a| self.w[ij[a]] / 2.0).product()
    }

    fn limit(&self, q: &mut [f64], lim: &DgLimiter, xi: &[f64]) {
        let (nv, np, n1) = (self.nv, self.np, self.n1);
        let m = &self.mesh;
        let ne = m.elements();
        let mean_of = |q: &[f64], e: usize, var: usize| -> f64 {
            (0..np).map(|nd| self.node_weight(nd) * q[(e * np + nd) * nv + var]).sum()
        };
        if let Some(mt) = lim.tvb_m {
            let means: Vec<f64> = (0..ne * nv).map(|k| mean_of(q, k / nv, k % nv)).collect();
            let old = q.to_vec();
            for e in 0..ne {
                for var in 0..nv {
                    let mean = means[e * nv + var];
                    let mut slopes = [0.0; 2];
                    let mut troubled = false;
                    for a in 0..m.dim {
                        let (mut ul, mut ur) = (0.0, 0.0);
                        let tr = if m.dim == 2 { n1 } else { 1 };
                        for t in 0..tr {
                            let wt = if m.dim == 2 { self.w[t] / 2.0 } else { 1.0 };
                            let (nl, nr) = if a == 0 { (self.node(0, t), self.node(n1 - 1, t)) } else { (self.node(t, 0), self.node(t, n1 - 1)) };
                            ul += wt * old[(e * np + nl) * nv + var];
                            ur += wt * old[(e * np + nr) * nv + var];
                        }
                        let mp = self.neighbour(e, a, true).map(|k| means[k * nv + var]);
                        let mm = self.neighbour(e, a, false).map(|k| means[k * nv + var]);
                        let (dp, dm) = match (mp, mm) {
                            (Some(p), Some(mi)) => (p - mean, mean - mi),
                            (Some(p), None) => (p - mean, p - mean),
                            (None, Some(mi)) => (mean - mi, mean - mi),
                            (None, None) => (0.0, 0.0),
                        };
                        let thr = mt * m.h(a) * m.h(a);
                        let tvb = |x: f64| if x.abs() <= thr { x } else { minmod3(x, dp, dm) };
                        let (up, um) = (ur - mean, mean - ul);
                        let (sp, sm) = (tvb(up), tvb(um));
                        let scale = 1e-9 * (mean.abs() + up.abs() + um.abs() + 1e-300);
                        if (sp - up).abs() > scale || (sm - um).abs() > scale {
                            troubled = true;
                            slopes[a] = 0.5 * (sp + sm);
                        } else {
                            slopes[a] = 0.5 * (up + um);
                        }
                    }
                    if troubled {
                        for nd in 0..np {
                            let ij = [nd % n1, nd / n1];
                            let mut v = mean;
                            for a in 0..m.dim {
                                v += slopes[a] * xi[ij[a]];
                            }
                            q[(e * np + nd) * nv + var] = v;
                        }
                    }
                }
            }
        }
        if lim.positivity {
            let eps = 1e-10;
            for e in 0..ne {
                let mean: Vec<f64> = (0..nv).map(|v| mean_of(q, e, v)).collect();
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
                    let node = &q[base..base + nv];
                    let ok = |s: f64| {
                        let st: Vec<f64> = (0..nv).map(|v| mean[v] + s * (node[v] - mean[v])).collect();
                        self.law.rho_p(&st).is_some_and(|(r, p)| r >= eps && p >= eps)
                    };
                    if !ok(1.0) {
                        let (mut lo, mut hi) = (0.0, 1.0);
                        for _ in 0..40 {
                            let mid = 0.5 * (lo + hi);
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
}

/// Integrates `q_t + ∇·F(q) = 0` with nodal DG to `t_end`.
///
/// `bcs` lists `[XMin, XMax, YMin, YMax]`; 1D meshes use the first two.
///
/// # Errors
/// If the mesh is invalid, periodic faces are inconsistent or the initial
/// data have the wrong number of components.
pub fn solve_dg(
    mesh: &DgMesh,
    law: &dyn DgLaw,
    bcs: &[DgBc; 4],
    init: &dyn Fn([f64; 2]) -> Vec<f64>,
    t_end: f64,
    cfl: f64,
    filter: Option<ModalFilter>,
) -> Result<DgSolution, String> {
    solve_dg_limited(mesh, law, bcs, init, t_end, cfl, filter, None)
}

/// As [`solve_dg`] with an optional slope / positivity limiter applied after
/// every Runge–Kutta stage.
///
/// # Errors
/// See [`solve_dg`].
pub fn solve_dg_limited(
    mesh: &DgMesh,
    law: &dyn DgLaw,
    bcs: &[DgBc; 4],
    init: &dyn Fn([f64; 2]) -> Vec<f64>,
    t_end: f64,
    cfl: f64,
    filter: Option<ModalFilter>,
    limiter: Option<DgLimiter>,
) -> Result<DgSolution, String> {
    if !(1..=2).contains(&mesh.dim) || mesh.p == 0 || mesh.ne.iter().any(|&n| n == 0) {
        return Err("invalid DG mesh".into());
    }
    for a in 0..mesh.dim {
        if matches!(bcs[2 * a], DgBc::Periodic) != matches!(bcs[2 * a + 1], DgBc::Periodic) {
            return Err("periodic must be set on both faces of an axis".into());
        }
    }
    let (xi, w) = lgl(mesh.p);
    let d = lgl_diff_matrix(&xi);
    let nv = law.nvar();
    let n1 = mesh.p + 1;
    let np = mesh.nodes_per_element();
    let ctx = Ctx { mesh: *mesh, law, bcs, d, w: w.clone(), nv, np, n1 };
    let mut sol = DgSolution { mesh: *mesh, xi: xi.clone(), w, q: vec![0.0; mesh.elements() * np * nv], nvar: nv, t: 0.0, steps: 0 };
    for e in 0..mesh.elements() {
        for nd in 0..np {
            let v = init(sol.coord(e, nd));
            if v.len() != nv {
                return Err("initial data has the wrong number of components".into());
            }
            for c in 0..nv {
                sol.q[(e * np + nd) * nv + c] = v[c];
            }
        }
    }
    // modal filter matrix
    let fmat = filter.map(|f| {
        let v: Vec<Vec<f64>> = xi.iter().map(|&x| legendre(mesh.p, x)).collect();
        let vi = invert(v.clone());
        let sig: Vec<f64> = (0..n1).map(|k| (-f.alpha * (k as f64 / mesh.p as f64).powi(f.order as i32)).exp()).collect();
        let mut m = vec![vec![0.0; n1]; n1];
        for i in 0..n1 {
            for j in 0..n1 {
                m[i][j] = (0..n1).map(|k| v[i][k] * sig[k] * vi[k][j]).sum();
            }
        }
        m
    });
    let apply_filter = |q: &mut Vec<f64>| {
        let Some(fm) = &fmat else { return };
        for a in 0..mesh.dim {
            let src = q.clone();
            let tr = if mesh.dim == 2 { n1 } else { 1 };
            for e in 0..mesh.elements() {
                for t_idx in 0..tr {
                    for i in 0..n1 {
                        let nd = if a == 0 { i + n1 * t_idx } else { t_idx + n1 * i };
                        for c in 0..nv {
                            let mut s = 0.0;
                            for k in 0..n1 {
                                let nk = if a == 0 { k + n1 * t_idx } else { t_idx + n1 * k };
                                s += fm[i][k] * src[(e * np + nk) * nv + c];
                            }
                            q[(e * np + nd) * nv + c] = s;
                        }
                    }
                }
            }
        }
    };
    let mut t = 0.0;
    let mut q = sol.q.clone();
    let mut steps = 0;
    while t < t_end * (1.0 - 1e-14) {
        let rate = ctx.max_rate(&q);
        if !rate.is_finite() || rate <= 0.0 {
            return Err("non-finite wave speed".into());
        }
        let dt = (cfl / ((2 * mesh.p + 1) as f64 * rate)).min(t_end - t);
        let l0 = ctx.rhs(&q, t, &xi);
        let mut q1: Vec<f64> = q.iter().zip(&l0).map(|(a, b)| a + dt * b).collect();
        if let Some(lm) = &limiter {
            ctx.limit(&mut q1, lm, &xi);
        }
        let l1 = ctx.rhs(&q1, t + dt, &xi);
        let mut q2: Vec<f64> = (0..q.len()).map(|i| 0.75 * q[i] + 0.25 * (q1[i] + dt * l1[i])).collect();
        if let Some(lm) = &limiter {
            ctx.limit(&mut q2, lm, &xi);
        }
        let l2 = ctx.rhs(&q2, t + 0.5 * dt, &xi);
        q = (0..q.len()).map(|i| q[i] / 3.0 + 2.0 / 3.0 * (q2[i] + dt * l2[i])).collect();
        if let Some(lm) = &limiter {
            ctx.limit(&mut q, lm, &xi);
        }
        if fmat.is_some() {
            apply_filter(&mut q);
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
