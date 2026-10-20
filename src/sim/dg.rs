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
//! | Euler ([`Euler1d`], [`Euler2d`]) | yes | yes | same |
//!
//! Shocks are handled by [`solve_dg_limited`]: a TVB minmod slope limiter on
//! troubled elements and the Zhang–Shu positivity-preserving scaling of density
//! and pressure, applied after every Runge–Kutta stage (an exponential modal
//! filter is also available). Custom systems implement [`DgLaw`]. Meshes are uniform tensor-product
//! quadrilaterals (intervals in 1D); no hanging nodes, no p-adaptivity.
//! Smooth problems converge with order `p + 1`; shocks need the filter and
//! are first-order accurate locally.

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
            out[c] = f64::midpoint(fl[c], fr[c]) - 0.5 * s * (qr[c] - ql[c]);
        }
    }
    /// Two-point volume flux `F#(q_l, q_r)` for flux differencing (consistent
    /// and symmetric; the default arithmetic mean reproduces the strong form,
    /// laws with an entropy-conservative flux override it).
    fn two_point_flux(&self, ql: &[f64], qr: &[f64], dir: usize, out: &mut [f64]) {
        let n = self.nvar();
        let mut fl = vec![0.0; n];
        let mut fr = vec![0.0; n];
        self.flux(ql, dir, &mut fl);
        self.flux(qr, dir, &mut fr);
        for c in 0..n {
            out[c] = f64::midpoint(fl[c], fr[c]);
        }
    }
    /// Entropy-stable interface flux: the two-point flux with a local
    /// Lax–Friedrichs dissipation on the conserved variables.
    fn entropy_stable_flux(&self, ql: &[f64], qr: &[f64], dir: usize, out: &mut [f64]) {
        self.two_point_flux(ql, qr, dir, out);
        let s = self.max_speed(ql, dir).max(self.max_speed(qr, dir));
        for c in 0..self.nvar() {
            out[c] -= 0.5 * s * (qr[c] - ql[c]);
        }
    }
    /// Whether [`DgLaw::to_characteristic`] is implemented.
    fn has_characteristics(&self) -> bool {
        false
    }
    /// Characteristic variables `L(q_ref) v` of the flux Jacobian in direction `dir` (identity by default).
    fn to_characteristic(&self, _q_ref: &[f64], _dir: usize, v: &[f64], out: &mut [f64]) {
        out[..v.len()].copy_from_slice(v);
    }
    /// Inverse of [`DgLaw::to_characteristic`]: `R(q_ref) w`.
    fn from_characteristic(&self, _q_ref: &[f64], _dir: usize, w: &[f64], out: &mut [f64]) {
        out[..w.len()].copy_from_slice(w);
    }
    /// Exterior state of a slip wall normal to `dir` (copy by default).
    fn wall_state(&self, q: &[f64], _dir: usize, out: &mut [f64]) {
        out[..q.len()].copy_from_slice(q);
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
    fn two_point_flux(&self, ql: &[f64], qr: &[f64], dir: usize, out: &mut [f64]) {
        euler_ec_flux(self.gamma, 1, dir, ql, qr, out);
    }
    fn has_characteristics(&self) -> bool {
        true
    }
    fn to_characteristic(&self, q_ref: &[f64], dir: usize, v: &[f64], out: &mut [f64]) {
        euler_eigen_transform(self.gamma, 1, dir, q_ref, v, out, true);
    }
    fn from_characteristic(&self, q_ref: &[f64], dir: usize, w: &[f64], out: &mut [f64]) {
        euler_eigen_transform(self.gamma, 1, dir, q_ref, w, out, false);
    }
    fn wall_state(&self, q: &[f64], dir: usize, out: &mut [f64]) {
        euler_wall_state(1, dir, q, out);
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
        (self.gamma - 1.0) * (q[3] - f64::midpoint(q[1] * q[1], q[2] * q[2]) / q[0])
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
    fn two_point_flux(&self, ql: &[f64], qr: &[f64], dir: usize, out: &mut [f64]) {
        euler_ec_flux(self.gamma, 2, dir, ql, qr, out);
    }
    fn has_characteristics(&self) -> bool {
        true
    }
    fn to_characteristic(&self, q_ref: &[f64], dir: usize, v: &[f64], out: &mut [f64]) {
        euler_eigen_transform(self.gamma, 2, dir, q_ref, v, out, true);
    }
    fn from_characteristic(&self, q_ref: &[f64], dir: usize, w: &[f64], out: &mut [f64]) {
        euler_eigen_transform(self.gamma, 2, dir, q_ref, w, out, false);
    }
    fn wall_state(&self, q: &[f64], dir: usize, out: &mut [f64]) {
        euler_wall_state(2, dir, q, out);
    }
    fn max_speed(&self, q: &[f64], dir: usize) -> f64 {
        let p = self.pressure(q).max(0.0);
        (q[1 + dir] / q[0]).abs() + (self.gamma * p / q[0]).sqrt()
    }
    fn rho_p(&self, q: &[f64]) -> Option<(f64, f64)> {
        Some((q[0], self.pressure(q)))
    }
}


/// Logarithmic mean `(a - b) / (ln a - ln b)` with a series for close arguments.
#[must_use]
pub fn log_mean(a: f64, b: f64) -> f64 {
    let xi = (a - b) / (a + b);
    let f = xi * xi;
    if f < 1e-3 {
        (a + b) / (2.0 * (1.0 + f / 3.0 + f * f / 5.0 + f * f * f / 7.0))
    } else {
        (a - b) / (a.ln() - b.ln())
    }
}

/// Primitive variables `(ρ, velocities, p)` of an Euler state with `nvel` velocity components.
pub(crate) fn euler_prim(gamma: f64, nvel: usize, q: &[f64]) -> (f64, [f64; 3], f64) {
    let rho = q[0];
    let mut vel = [0.0; 3];
    let mut ke = 0.0;
    for k in 0..nvel {
        vel[k] = q[1 + k] / rho;
        ke += 0.5 * rho * vel[k] * vel[k];
    }
    (rho, vel, (gamma - 1.0) * (q[1 + nvel] - ke))
}

/// Entropy-conservative flux of Chandrashekar (2013) for the Euler equations
/// with `nvel` (1 or 2) velocity components in direction `dir`.
pub(crate) fn euler_ec_flux(gamma: f64, nvel: usize, dir: usize, ql: &[f64], qr: &[f64], out: &mut [f64]) {
    let (rl, vl, pl) = euler_prim(gamma, nvel, ql);
    let (rr, vr, pr) = euler_prim(gamma, nvel, qr);
    let (bl, br) = (rl / (2.0 * pl), rr / (2.0 * pr));
    let rho_ln = log_mean(rl, rr);
    let beta_ln = log_mean(bl, br);
    let rho_avg = 0.5 * (rl + rr);
    let beta_avg = 0.5 * (bl + br);
    let p_hat = rho_avg / (2.0 * beta_avg);
    let mut ubar = [0.0; 3];
    let mut v2 = 0.0;
    for k in 0..nvel {
        ubar[k] = 0.5 * (vl[k] + vr[k]);
        v2 += 0.5 * (vl[k] * vl[k] + vr[k] * vr[k]);
    }
    let f_rho = rho_ln * ubar[dir];
    out[0] = f_rho;
    let mut dot = 0.0;
    for k in 0..nvel {
        let m = f_rho * ubar[k] + if k == dir { p_hat } else { 0.0 };
        out[1 + k] = m;
        dot += ubar[k] * m;
    }
    out[1 + nvel] = f_rho * (1.0 / (2.0 * (gamma - 1.0) * beta_ln) - 0.5 * v2) + dot;
}

/// Applies the left (`to_char`) or right eigenvector matrix of the Euler flux
/// Jacobian in direction `dir`, evaluated at the state `q`, to `v`.
pub(crate) fn euler_eigen_transform(gamma: f64, nvel: usize, dir: usize, q: &[f64], v: &[f64], out: &mut [f64], to_char: bool) {
    let (rho, vel, p) = euler_prim(gamma, nvel, q);
    let c = (gamma * p.max(1e-300) / rho).sqrt();
    let un = vel[dir];
    let ut = if nvel == 2 { vel[1 - dir] } else { 0.0 };
    let q2 = un * un + ut * ut;
    let h = (q[1 + nvel] + p) / rho;
    let b1 = (gamma - 1.0) / (c * c);
    let b2 = 0.5 * b1 * q2;
    if nvel == 1 {
        if to_char {
            let l = [
                [0.5 * (b2 + un / c), -0.5 * (b1 * un + 1.0 / c), 0.5 * b1],
                [1.0 - b2, b1 * un, -b1],
                [0.5 * (b2 - un / c), -0.5 * (b1 * un - 1.0 / c), 0.5 * b1],
            ];
            for i in 0..3 {
                out[i] = l[i][0] * v[0] + l[i][1] * v[1] + l[i][2] * v[2];
            }
        } else {
            let r = [[1.0, 1.0, 1.0], [un - c, un, un + c], [h - un * c, 0.5 * un * un, h + un * c]];
            for i in 0..3 {
                out[i] = r[i][0] * v[0] + r[i][1] * v[1] + r[i][2] * v[2];
            }
        }
        return;
    }
    // 2D: rotate so that the sweep direction is the first velocity component
    let t = 1 - dir;
    if to_char {
        let w = [v[0], v[1 + dir], v[1 + t], v[3]];
        let l = [
            [0.5 * (b2 + un / c), -0.5 * (b1 * un + 1.0 / c), -0.5 * b1 * ut, 0.5 * b1],
            [-ut, 0.0, 1.0, 0.0],
            [1.0 - b2, b1 * un, b1 * ut, -b1],
            [0.5 * (b2 - un / c), -0.5 * (b1 * un - 1.0 / c), -0.5 * b1 * ut, 0.5 * b1],
        ];
        for i in 0..4 {
            out[i] = (0..4).map(|k| l[i][k] * w[k]).sum();
        }
    } else {
        let r = [
            [1.0, 0.0, 1.0, 1.0],
            [un - c, 0.0, un, un + c],
            [ut, 1.0, ut, ut],
            [h - un * c, ut, 0.5 * q2, h + un * c],
        ];
        let mut w = [0.0; 4];
        for i in 0..4 {
            w[i] = (0..4).map(|k| r[i][k] * v[k]).sum();
        }
        out[0] = w[0];
        out[1 + dir] = w[1];
        out[1 + t] = w[2];
        out[3] = w[3];
    }
}

/// Ghost state of a slip wall normal to `dir`: the normal momentum is reversed.
pub(crate) fn euler_wall_state(nvel: usize, dir: usize, q: &[f64], out: &mut [f64]) {
    out[..nvel + 2].copy_from_slice(&q[..nvel + 2]);
    out[1 + dir] = -q[1 + dir];
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
    /// Slip wall: the exterior state is [`DgLaw::wall_state`] of the interior one.
    Wall,
    /// Exterior state `q_ext(t, x, q_interior)` computed from the interior
    /// trace (mixed wall / inflow faces, time-dependent states).
    Custom(Arc<dyn Fn(f64, [f64; 2], &[f64]) -> Vec<f64> + Send + Sync>),
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
    pub const fn line(ne: usize, lo: f64, hi: f64, p: usize) -> Self {
        Self { dim: 1, ne: [ne, 1], lo: [lo, 0.0], hi: [hi, 1.0], p }
    }

    /// Two-dimensional mesh.
    #[must_use]
    pub const fn rect(ne: [usize; 2], lo: [f64; 2], hi: [f64; 2], p: usize) -> Self {
        Self { dim: 2, ne, lo, hi, p }
    }

    fn h(&self, a: usize) -> f64 {
        (self.hi[a] - self.lo[a]) / self.ne[a] as f64
    }

    /// Nodes per element.
    #[must_use]
    pub const fn nodes_per_element(&self) -> usize {
        (self.p + 1).pow(self.dim as u32)
    }

    /// Number of elements.
    #[must_use]
    pub const fn elements(&self) -> usize {
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
    const fn node_ij(&self, node: usize) -> [usize; 2] {
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
            x[a] = m.lo[a] + (ex[a] as f64 + f64::midpoint(1.0, self.xi[ij[a]])) * m.h(a);
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
    es: bool,
    vmat: Vec<Vec<f64>>,
    vinv: Vec<Vec<f64>>,
}

impl Ctx<'_> {
    const fn node(&self, i: usize, j: usize) -> usize {
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
                            None => {
                                let bc = &self.bcs[2 * a + usize::from(max)];
                                let xpos = || {
                                    let ij = if a == 0 { [edge, t_idx] } else { [t_idx, edge] };
                                    let ex = [e % m.ne[0], e / m.ne[0]];
                                    let mut x = [0.0; 2];
                                    for b in 0..m.dim {
                                        x[b] = m.lo[b] + (ex[b] as f64 + f64::midpoint(1.0, xi[ij[b]])) * m.h(b);
                                    }
                                    x
                                };
                                match bc {
                                    DgBc::Fixed(f) => qo.copy_from_slice(&f(t, xpos())),
                                    DgBc::Custom(f) => qo.copy_from_slice(&f(t, xpos(), qi)),
                                    DgBc::Wall => self.law.wall_state(qi, a, &mut qo),
                                    _ => qo.copy_from_slice(qi),
                                }
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
    const fn neighbour(&self, e: usize, axis: usize, max: bool) -> Option<usize> {
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
                            slopes[a] = f64::midpoint(sp, sm);
                        } else {
                            slopes[a] = f64::midpoint(up, um);
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
#[allow(clippy::too_many_arguments)] // full solver entry point; `solve_dg` is the shorter wrapper
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
    let opt = DgOptions {
        filter,
        shock: limiter.and_then(|l| l.tvb_m).map_or(ShockLimiter::None, |m| ShockLimiter::Tvb { m }),
        positivity: limiter.is_some_and(|l| l.positivity),
        entropy_stable: false,
    };
    solve_dg_with(mesh, law, bcs, init, t_end, cfl, &opt)
}

/// Shock-capturing limiter applied after every Runge–Kutta stage.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ShockLimiter {
    /// No limiting.
    None,
    /// TVB minmod slope limiter on each conserved variable (see [`DgLimiter`]).
    Tvb {
        /// TVB constant `M`.
        m: f64,
    },
    /// TVB minmod limiter applied to the **characteristic variables** built
    /// from the Euler eigenvectors at the element mean ([`DgLaw::has_characteristics`]).
    TvbCharacteristic {
        /// TVB constant `M`.
        m: f64,
    },
    /// Krivodonova's hierarchical moment limiter (modal coefficients are
    /// bounded from the highest degree down), 1D and 2D.
    Hierarchical,
    /// Zhong–Shu WENO reconstruction of troubled elements (1D). Troubled
    /// elements are detected with the TVB test of constant `m`.
    Weno {
        /// TVB constant `M` of the troubled-element indicator.
        m: f64,
    },
}

/// Options of [`solve_dg_with`].
#[derive(Clone, Copy, Debug)]
pub struct DgOptions {
    /// Optional exponential modal filter.
    pub filter: Option<ModalFilter>,
    /// Shock-capturing limiter.
    pub shock: ShockLimiter,
    /// Zhang–Shu positivity-preserving scaling of density and pressure.
    pub positivity: bool,
    /// Entropy-stable DGSEM: flux-differencing volume terms with the
    /// two-point flux of the law ([`DgLaw::two_point_flux`]) and interface
    /// fluxes `F_ec - λ/2 [[q]]` ([`DgLaw::entropy_stable_flux`]).
    pub entropy_stable: bool,
}

impl Default for DgOptions {
    fn default() -> Self {
        Self { filter: None, shock: ShockLimiter::None, positivity: false, entropy_stable: false }
    }
}

impl Ctx<'_> {
    fn element_mean(&self, q: &[f64], e: usize, var: usize) -> f64 {
        (0..self.np).map(|nd| self.node_weight(nd) * q[(e * self.np + nd) * self.nv + var]).sum()
    }

    /// TVB limiter on characteristic variables.
    fn limit_tvb_characteristic(&self, q: &mut [f64], mt: f64, xi: &[f64]) {
        let (nv, np, n1) = (self.nv, self.np, self.n1);
        let m = &self.mesh;
        let ne = m.elements();
        let means: Vec<f64> = (0..ne * nv).map(|k| self.element_mean(q, k / nv, k % nv)).collect();
        let old = q.to_vec();
        let mut buf = vec![0.0; nv];
        for e in 0..ne {
            let mean: Vec<f64> = means[e * nv..(e + 1) * nv].to_vec();
            let mut slopes = vec![[0.0; 2]; nv];
            let mut troubled = false;
            for a in 0..m.dim {
                let (mut ul, mut ur) = (vec![0.0; nv], vec![0.0; nv]);
                let tr = if m.dim == 2 { n1 } else { 1 };
                for t in 0..tr {
                    let wt = if m.dim == 2 { self.w[t] / 2.0 } else { 1.0 };
                    let (nl, nr) = if a == 0 { (self.node(0, t), self.node(n1 - 1, t)) } else { (self.node(t, 0), self.node(t, n1 - 1)) };
                    for v in 0..nv {
                        ul[v] += wt * old[(e * np + nl) * nv + v];
                        ur[v] += wt * old[(e * np + nr) * nv + v];
                    }
                }
                let mp = self.neighbour(e, a, true);
                let mm = self.neighbour(e, a, false);
                let mut up = vec![0.0; nv];
                let mut um = vec![0.0; nv];
                let mut dp = vec![0.0; nv];
                let mut dm = vec![0.0; nv];
                for v in 0..nv {
                    up[v] = ur[v] - mean[v];
                    um[v] = mean[v] - ul[v];
                    let (p, mi) = (mp.map(|k| means[k * nv + v]), mm.map(|k| means[k * nv + v]));
                    let (a1, a2) = match (p, mi) {
                        (Some(p), Some(mi)) => (p - mean[v], mean[v] - mi),
                        (Some(p), None) => (p - mean[v], p - mean[v]),
                        (None, Some(mi)) => (mean[v] - mi, mean[v] - mi),
                        (None, None) => (0.0, 0.0),
                    };
                    dp[v] = a1;
                    dm[v] = a2;
                }
                let tf = |vec: &[f64], buf: &mut [f64]| self.law.to_characteristic(&mean, a, vec, buf);
                let (mut cup, mut cum, mut cdp, mut cdm) = (vec![0.0; nv], vec![0.0; nv], vec![0.0; nv], vec![0.0; nv]);
                tf(&up, &mut buf);
                cup.copy_from_slice(&buf);
                tf(&um, &mut buf);
                cum.copy_from_slice(&buf);
                tf(&dp, &mut buf);
                cdp.copy_from_slice(&buf);
                tf(&dm, &mut buf);
                cdm.copy_from_slice(&buf);
                let thr = mt * m.h(a) * m.h(a);
                let mut cs = vec![0.0; nv];
                for k in 0..nv {
                    let tvb = |x: f64| if x.abs() <= thr { x } else { minmod3(x, cdp[k], cdm[k]) };
                    let (sp, sm) = (tvb(cup[k]), tvb(cum[k]));
                    let scale = 1e-9 * (mean.iter().map(|x| x.abs()).sum::<f64>() + cup[k].abs() + cum[k].abs() + 1e-300);
                    if (sp - cup[k]).abs() > scale || (sm - cum[k]).abs() > scale {
                        troubled = true;
                        cs[k] = f64::midpoint(sp, sm);
                    } else {
                        cs[k] = f64::midpoint(cup[k], cum[k]);
                    }
                }
                self.law.from_characteristic(&mean, a, &cs, &mut buf);
                for v in 0..nv {
                    slopes[v][a] = buf[v];
                }
            }
            if troubled {
                for nd in 0..np {
                    let ij = [nd % n1, nd / n1];
                    for v in 0..nv {
                        let mut val = mean[v];
                        for a in 0..m.dim {
                            val += slopes[v][a] * xi[ij[a]];
                        }
                        q[(e * np + nd) * nv + v] = val;
                    }
                }
            }
        }
    }

    /// Nodal-to-modal (Legendre) coefficients of one element and variable, tensor layout `k + n1 l`.
    fn to_modal(&self, vals: &[f64]) -> Vec<f64> {
        let n1 = self.n1;
        let vi = &self.vinv;
        if self.mesh.dim == 1 {
            return (0..n1).map(|k| (0..n1).map(|j| vi[k][j] * vals[j]).sum()).collect();
        }
        let mut tmp = vec![0.0; self.np];
        for l in 0..n1 {
            for k in 0..n1 {
                tmp[k + n1 * l] = (0..n1).map(|i| vi[k][i] * vals[i + n1 * l]).sum();
            }
        }
        let mut out = vec![0.0; self.np];
        for l in 0..n1 {
            for k in 0..n1 {
                out[k + n1 * l] = (0..n1).map(|j| vi[l][j] * tmp[k + n1 * j]).sum();
            }
        }
        out
    }

    fn from_modal(&self, c: &[f64]) -> Vec<f64> {
        let n1 = self.n1;
        let v = &self.vmat;
        if self.mesh.dim == 1 {
            return (0..n1).map(|i| (0..n1).map(|k| v[i][k] * c[k]).sum()).collect();
        }
        let mut tmp = vec![0.0; self.np];
        for l in 0..n1 {
            for i in 0..n1 {
                tmp[i + n1 * l] = (0..n1).map(|k| v[i][k] * c[k + n1 * l]).sum();
            }
        }
        let mut out = vec![0.0; self.np];
        for j in 0..n1 {
            for i in 0..n1 {
                out[i + n1 * j] = (0..n1).map(|l| v[j][l] * tmp[i + n1 * l]).sum();
            }
        }
        out
    }

    /// Krivodonova's hierarchical moment limiter.
    fn limit_hierarchical(&self, q: &mut [f64]) {
        let (nv, np, n1) = (self.nv, self.np, self.n1);
        let m = &self.mesh;
        let ne = m.elements();
        let p = n1 - 1;
        let mut modal = vec![0.0; ne * nv * np];
        for e in 0..ne {
            for v in 0..nv {
                let vals: Vec<f64> = (0..np).map(|nd| q[(e * np + nd) * nv + v]).collect();
                let c = self.to_modal(&vals);
                modal[(e * nv + v) * np..(e * nv + v + 1) * np].copy_from_slice(&c);
            }
        }
        let old = modal.clone();
        let mm3 = |a: f64, b: f64, c: f64| minmod3(a, b, c);
        for e in 0..ne {
            for v in 0..nv {
                let base = (e * nv + v) * np;
                let mut changed_element = false;
                let top = m.dim * p;
                for deg in (1..=top).rev() {
                    let mut any = false;
                    for l in 0..if m.dim == 2 { n1 } else { 1 } {
                        for k in 0..n1 {
                            if k + l != deg {
                                continue;
                            }
                            let idx = k + n1 * l;
                            let c = old[base + idx];
                            let mut cands: Vec<f64> = Vec::new();
                            for a in 0..m.dim {
                                let order = if a == 0 { k } else { l };
                                if order == 0 {
                                    continue;
                                }
                                let lower = if a == 0 { (k - 1) + n1 * l } else { k + n1 * (l - 1) };
                                let own = old[base + lower];
                                let (np_, nm_) = (self.neighbour(e, a, true), self.neighbour(e, a, false));
                                let (Some(ep), Some(em)) = (np_, nm_) else { continue };
                                let plus = old[(ep * nv + v) * np + lower];
                                let minus = old[(em * nv + v) * np + lower];
                                let f = (2 * order - 1) as f64;
                                cands.push(mm3(f * c, plus - own, own - minus) / f);
                            }
                            if cands.is_empty() {
                                continue;
                            }
                            let lim = if cands.len() == 1 {
                                cands[0]
                            } else if cands.iter().all(|&x| x > 0.0) {
                                cands.iter().fold(f64::INFINITY, |a, &b| a.min(b))
                            } else if cands.iter().all(|&x| x < 0.0) {
                                cands.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b))
                            } else {
                                0.0
                            };
                            if (lim - c).abs() > 1e-12 * (c.abs() + 1e-300) + 1e-14 {
                                modal[base + idx] = lim;
                                any = true;
                            }
                        }
                    }
                    if any {
                        changed_element = true;
                    } else {
                        break;
                    }
                }
                if changed_element {
                    let vals = self.from_modal(&modal[base..base + np]);
                    for nd in 0..np {
                        q[(e * np + nd) * nv + v] = vals[nd];
                    }
                }
            }
        }
    }

    /// Value at reference coordinate `x` of the interpolant through the LGL nodes `xi` with nodal values `f`.
    fn lagrange_eval(xi: &[f64], f: &[f64], x: f64) -> f64 {
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

    /// Zhong–Shu WENO limiter for troubled elements (1D).
    fn limit_weno(&self, q: &mut [f64], mt: f64, xi: &[f64]) {
        let (nv, np) = (self.nv, self.np);
        let ne = self.mesh.elements();
        let mut trial = q.to_vec();
        let tvb = DgLimiter { tvb_m: Some(mt), positivity: false };
        self.limit(&mut trial, &tvb, xi);
        let old = q.to_vec();
        for e in 0..ne {
            let troubled = (0..np * nv).any(|k| (trial[e * np * nv + k] - old[e * np * nv + k]).abs() > 1e-9 * (old[e * np * nv + k].abs() + 1e-12));
            if !troubled {
                continue;
            }
            for v in 0..nv {
                let own: Vec<f64> = (0..np).map(|nd| old[(e * np + nd) * nv + v]).collect();
                let mean_i: f64 = (0..np).map(|nd| 0.5 * self.w[nd] * own[nd]).sum();
                let mut cands: Vec<(Vec<f64>, f64)> = vec![(own.clone(), 0.998)];
                for (max, shift) in [(false, -2.0), (true, 2.0)] {
                    if let Some(en) = self.neighbour(e, 0, max) {
                        let nb: Vec<f64> = (0..np).map(|nd| old[(en * np + nd) * nv + v]).collect();
                        let vals: Vec<f64> = (0..np).map(|j| Self::lagrange_eval(xi, &nb, xi[j] - shift)).collect();
                        cands.push((vals, 0.001));
                    }
                }
                if cands.len() == 2 {
                    cands[0].1 = 0.999;
                }
                if cands.len() == 1 {
                    continue;
                }
                let eps = 1e-6 * (mean_i * mean_i + 1e-12);
                let mut omegas = Vec::new();
                for (vals, g) in &mut cands {
                    let mc: f64 = (0..np).map(|nd| 0.5 * self.w[nd] * vals[nd]).sum();
                    for x in vals.iter_mut() {
                        *x += mean_i - mc;
                    }
                    // oscillation indicator
                    let mut beta = 0.0;
                    let mut d = vals.clone();
                    for l in 1..=self.mesh.p {
                        d = (0..np).map(|i| (0..np).map(|k| self.d[i][k] * d[k]).sum::<f64>()).collect();
                        let scale = f64::from(1u32 << (2 * l - 1));
                        beta += scale * (0..np).map(|j| self.w[j] * d[j] * d[j]).sum::<f64>();
                    }
                    omegas.push(*g / ((eps + beta) * (eps + beta)));
                }
                let tot: f64 = omegas.iter().sum();
                for nd in 0..np {
                    q[(e * np + nd) * nv + v] = cands.iter().zip(&omegas).map(|((vals, _), w)| w / tot * vals[nd]).sum();
                }
            }
        }
    }

    /// Applies the shock limiter and the positivity scaling.
    fn limit_all(&self, q: &mut [f64], opt: &DgOptions, xi: &[f64]) {
        match opt.shock {
            ShockLimiter::None => {},
            ShockLimiter::Tvb { m } => self.limit(q, &DgLimiter { tvb_m: Some(m), positivity: false }, xi),
            ShockLimiter::TvbCharacteristic { m } => {
                if self.law.has_characteristics() {
                    self.limit_tvb_characteristic(q, m, xi);
                } else {
                    self.limit(q, &DgLimiter { tvb_m: Some(m), positivity: false }, xi);
                }
            },
            ShockLimiter::Hierarchical => self.limit_hierarchical(q),
            ShockLimiter::Weno { m } => {
                if self.mesh.dim == 1 {
                    self.limit_weno(q, m, xi);
                } else {
                    self.limit(q, &DgLimiter { tvb_m: Some(m), positivity: false }, xi);
                }
            },
        }
        if opt.positivity {
            self.limit(q, &DgLimiter { tvb_m: None, positivity: true }, xi);
        }
    }
}

/// Integrates `q_t + ∇·F(q) = 0` with nodal DG and the full option set
/// (characteristic / hierarchical / WENO limiters, entropy-stable fluxes).
///
/// # Errors
/// If the mesh is invalid, periodic faces are inconsistent or the initial
/// data have the wrong number of components.
pub fn solve_dg_with(
    mesh: &DgMesh,
    law: &dyn DgLaw,
    bcs: &[DgBc; 4],
    init: &dyn Fn([f64; 2]) -> Vec<f64>,
    t_end: f64,
    cfl: f64,
    opt: &DgOptions,
) -> Result<DgSolution, String> {
    if !(1..=2).contains(&mesh.dim) || mesh.p == 0 || mesh.ne.contains(&0) {
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
    let vmat: Vec<Vec<f64>> = xi.iter().map(|&x| legendre(mesh.p, x)).collect();
    let vinv = invert(vmat.clone());
    let ctx = Ctx { mesh: *mesh, law, bcs, d, w: w.clone(), nv, np, n1, es: opt.entropy_stable, vmat: vmat.clone(), vinv: vinv.clone() };
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
    let fmat = opt.filter.map(|f| {
        let sig: Vec<f64> = (0..n1).map(|k| (-f.alpha * (k as f64 / mesh.p as f64).powi(f.order as i32)).exp()).collect();
        let mut m = vec![vec![0.0; n1]; n1];
        for i in 0..n1 {
            for j in 0..n1 {
                m[i][j] = (0..n1).map(|k| vmat[i][k] * sig[k] * vinv[k][j]).sum();
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
    let limiting = opt.shock != ShockLimiter::None || opt.positivity;
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
        if limiting {
            ctx.limit_all(&mut q1, opt, &xi);
        }
        let l1 = ctx.rhs(&q1, t + dt, &xi);
        let mut q2: Vec<f64> = (0..q.len()).map(|i| 0.75 * q[i] + 0.25 * (q1[i] + dt * l1[i])).collect();
        if limiting {
            ctx.limit_all(&mut q2, opt, &xi);
        }
        let l2 = ctx.rhs(&q2, t + 0.5 * dt, &xi);
        q = (0..q.len()).map(|i| q[i] / 3.0 + 2.0 / 3.0 * (q2[i] + dt * l2[i])).collect();
        if limiting {
            ctx.limit_all(&mut q, opt, &xi);
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
