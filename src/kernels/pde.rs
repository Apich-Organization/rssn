//! # Finite-difference PDE solvers
//!
//! Second-order finite-difference solvers on uniform grids:
//!
//! * [`heat_1d`] / [`heat_2d`]: `u_t = alpha * Laplace(u)` with the
//!   unconditionally stable Crank-Nicolson scheme (the sparse system is
//!   factored once with [`Csr::sparse_lu`]);
//! * [`wave_1d`] / [`wave_2d`]: `u_tt = c^2 * Laplace(u)` with the explicit
//!   leapfrog scheme and a CFL stability check;
//! * [`poisson_2d`] / [`laplace_2d`]: `-Laplace(u) = f` with the 5-point
//!   stencil solved by sparse LU or conjugate gradients;
//! * [`advection_1d`]: `u_t + a u_x = 0` on a periodic grid with first-order
//!   upwind or second-order Lax-Wendroff;
//! * [`method_of_lines_1d`]: semi-discretises `u_t = F(t, x, u, u_x, u_xx)`
//!   and integrates the resulting ODE system with the stiff solvers of
//!   [`crate::kernels::ode_stiff`].
//!
//! Boundary conditions are described by [`Boundary`]: Dirichlet (`u = g`) or
//! Neumann (`du/dn = g` with `n` the outward normal, imposed with a ghost
//! node so that the scheme stays second-order). Boundary data are functions
//! of the position `(x, y)` on the boundary (`y = 0` in one dimension).
//!
//! [`pde_solver`] is a single entry point dispatching on a [`PdeProblem`].
#![allow(
    clippy::manual_midpoint,
    clippy::missing_const_for_fn,
    clippy::struct_field_names,
    clippy::too_many_lines,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::cast_possible_truncation,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::needless_range_loop,
    clippy::too_many_arguments,
    clippy::many_single_char_names,
    clippy::type_complexity,
    clippy::suboptimal_flops,
    clippy::float_cmp,
    clippy::similar_names,
    clippy::while_float,
    clippy::needless_pass_by_value,
    clippy::neg_cmp_op_on_partial_ord,
    clippy::too_long_first_doc_paragraph
)]

use crate::kernels::krylov::{Csr, SparseOrdering, cg};
use crate::kernels::ode_adaptive::{OdeError, OdeOptions};
use crate::kernels::ode_stiff::{bdf, radau5};
use std::fmt;
use std::sync::Arc;

/// Boundary data as a function of the position `(x, y)` on the boundary.
pub type BcFn = dyn Fn(f64, f64) -> f64 + Send + Sync;

/// A boundary condition on one side of the domain.
#[derive(Clone)]
pub enum Boundary {
    /// `u = g(x, y)`.
    Dirichlet(Arc<BcFn>),
    /// `du/dn = g(x, y)` with `n` the outward unit normal.
    Neumann(Arc<BcFn>),
}

impl fmt::Debug for Boundary {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Dirichlet(_) => f.write_str("Boundary::Dirichlet(..)"),
            Self::Neumann(_) => f.write_str("Boundary::Neumann(..)"),
        }
    }
}

impl Boundary {
    /// Dirichlet condition with a position-dependent value.
    pub fn dirichlet(g: impl Fn(f64, f64) -> f64 + Send + Sync + 'static) -> Self {
        Self::Dirichlet(Arc::new(g))
    }

    /// Neumann condition with a position-dependent outward flux.
    pub fn neumann(g: impl Fn(f64, f64) -> f64 + Send + Sync + 'static) -> Self {
        Self::Neumann(Arc::new(g))
    }

    /// Constant Dirichlet value.
    #[must_use]
    pub fn dirichlet_const(v: f64) -> Self {
        Self::dirichlet(move |_, _| v)
    }

    /// Constant Neumann flux (`0.0` for an insulated / reflecting wall).
    #[must_use]
    pub fn neumann_const(v: f64) -> Self {
        Self::neumann(move |_, _| v)
    }

    fn is_dirichlet(&self) -> bool {
        matches!(self, Self::Dirichlet(_))
    }

    fn value(&self, x: f64, y: f64) -> f64 {
        match self {
            Self::Dirichlet(g) | Self::Neumann(g) => g(x, y),
        }
    }
}

/// Errors reported by the PDE solvers.
#[derive(Debug, Clone, PartialEq)]
pub enum PdeError {
    /// Invalid grid, step count, coefficient or callback result.
    InvalidInput(&'static str),
    /// The explicit scheme violates its stability (CFL) limit.
    Unstable {
        /// Courant number of the request.
        courant: f64,
        /// Largest stable Courant number.
        limit: f64,
    },
    /// The requested combination is not supported (e.g. singular pure-Neumann
    /// Poisson problem, or CG with Neumann sides).
    Unsupported(&'static str),
    /// The linear solve failed or did not converge.
    SolveFailed,
    /// The stiff ODE integrator failed.
    Ode(OdeError),
    /// The solution became non-finite.
    NonFinite,
}

impl fmt::Display for PdeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidInput(m) => write!(f, "invalid input: {m}"),
            Self::Unstable { courant, limit } => {
                write!(f, "unstable: Courant number {courant} exceeds the limit {limit}")
            }
            Self::Unsupported(m) => write!(f, "unsupported: {m}"),
            Self::SolveFailed => f.write_str("linear solve failed"),
            Self::Ode(e) => write!(f, "ODE integration failed: {e:?}"),
            Self::NonFinite => f.write_str("solution became non-finite"),
        }
    }
}

impl std::error::Error for PdeError {}

impl From<OdeError> for PdeError {
    fn from(e: OdeError) -> Self {
        Self::Ode(e)
    }
}

/// Uniform 1D grid on `[x0, x1]` with `n` nodes (including both ends;
/// for [`advection_1d`] the grid is periodic and `x1` is excluded).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Grid1D {
    /// Left end.
    pub x0: f64,
    /// Right end.
    pub x1: f64,
    /// Number of nodes.
    pub n: usize,
}

impl Grid1D {
    /// Grid spacing for a bounded grid.
    #[must_use]
    pub fn dx(&self) -> f64 {
        (self.x1 - self.x0) / (self.n.max(2) - 1) as f64
    }

    /// Node coordinates of a bounded grid.
    #[must_use]
    pub fn nodes(&self) -> Vec<f64> {
        let h = self.dx();
        (0..self.n).map(|i| self.x0 + i as f64 * h).collect()
    }

    fn validate(&self, min_n: usize) -> Result<(), PdeError> {
        if self.n < min_n || !(self.x1 > self.x0) || !self.x0.is_finite() || !self.x1.is_finite() {
            return Err(PdeError::InvalidInput("grid needs enough nodes and x1 > x0"));
        }
        Ok(())
    }
}

/// Uniform 2D grid on `[x0, x1] x [y0, y1]` with `nx * ny` nodes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Grid2D {
    /// Left end.
    pub x0: f64,
    /// Right end.
    pub x1: f64,
    /// Bottom end.
    pub y0: f64,
    /// Top end.
    pub y1: f64,
    /// Nodes in `x`.
    pub nx: usize,
    /// Nodes in `y`.
    pub ny: usize,
}

impl Grid2D {
    /// Unit square with `n x n` nodes.
    #[must_use]
    pub fn unit_square(n: usize) -> Self {
        Self { x0: 0.0, x1: 1.0, y0: 0.0, y1: 1.0, nx: n, ny: n }
    }

    /// Spacing in `x`.
    #[must_use]
    pub fn dx(&self) -> f64 {
        (self.x1 - self.x0) / (self.nx.max(2) - 1) as f64
    }

    /// Spacing in `y`.
    #[must_use]
    pub fn dy(&self) -> f64 {
        (self.y1 - self.y0) / (self.ny.max(2) - 1) as f64
    }

    fn validate(&self) -> Result<(), PdeError> {
        if self.nx < 3 || self.ny < 3 || !(self.x1 > self.x0) || !(self.y1 > self.y0) {
            return Err(PdeError::InvalidInput("2D grid needs >= 3 nodes per axis and positive extents"));
        }
        Ok(())
    }
}

/// A scalar field on a 1D grid.
#[derive(Debug, Clone, PartialEq)]
pub struct Field1D {
    /// Node coordinates.
    pub x: Vec<f64>,
    /// Values at the nodes.
    pub u: Vec<f64>,
    /// Time of the snapshot.
    pub t: f64,
}

/// A scalar field on a 2D grid, stored row-major (`u[j * nx + i]`).
#[derive(Debug, Clone, PartialEq)]
pub struct Field2D {
    /// Nodes in `x`.
    pub nx: usize,
    /// Nodes in `y`.
    pub ny: usize,
    /// Left end.
    pub x0: f64,
    /// Bottom end.
    pub y0: f64,
    /// Spacing in `x`.
    pub dx: f64,
    /// Spacing in `y`.
    pub dy: f64,
    /// Values at the nodes.
    pub u: Vec<f64>,
    /// Time of the snapshot (`0` for stationary problems).
    pub t: f64,
}

impl Field2D {
    /// Value at node `(i, j)`.
    #[must_use]
    pub fn at(&self, i: usize, j: usize) -> f64 {
        self.u[j * self.nx + i]
    }

    /// Coordinates of node `(i, j)`.
    #[must_use]
    pub fn point(&self, i: usize, j: usize) -> (f64, f64) {
        (self.x0 + i as f64 * self.dx, self.y0 + j as f64 * self.dy)
    }

    /// Maximum absolute difference to `exact(x, y)` over all nodes.
    pub fn max_error(&self, exact: impl Fn(f64, f64) -> f64) -> f64 {
        let mut m = 0.0_f64;
        for j in 0..self.ny {
            for i in 0..self.nx {
                let (x, y) = self.point(i, j);
                m = m.max((self.at(i, j) - exact(x, y)).abs());
            }
        }
        m
    }
}

impl Field1D {
    /// Maximum absolute difference to `exact(x)` over all nodes.
    pub fn max_error(&self, exact: impl Fn(f64) -> f64) -> f64 {
        self.x.iter().zip(&self.u).map(|(&x, &u)| (u - exact(x)).abs()).fold(0.0, f64::max)
    }
}

// ---------------------------------------------------------------------------
// Discrete Laplacian with boundary conditions
// ---------------------------------------------------------------------------

/// Discrete Laplacian `L u + b` with Dirichlet nodes marked (their rows are
/// empty: the value is imposed separately).
struct Laplacian {
    a: Csr,
    b: Vec<f64>,
    dirichlet: Vec<Option<f64>>,
}

/// Coefficients `(prev, centre, next, rhs)` of the second difference along one
/// axis at node `i` of `n`, using a ghost node at Neumann ends.
fn axis_stencil(i: usize, n: usize, h: f64, lo: &Boundary, hi: &Boundary, glo: f64, ghi: f64) -> (f64, f64, f64, f64) {
    let h2 = h * h;
    let (mut prev, mut next, mut rhs) = (1.0 / h2, 1.0 / h2, 0.0);
    if i == 0 && !lo.is_dirichlet() {
        // ghost u_{-1} = u_1 + 2 h g (outward derivative -u_x = g)
        (prev, next, rhs) = (0.0, 2.0 / h2, 2.0 * glo / h);
    }
    if i == n - 1 && !hi.is_dirichlet() {
        (prev, next, rhs) = (2.0 / h2, 0.0, 2.0 * ghi / h);
    }
    (prev, -2.0 / h2, next, rhs)
}

/// Builds the Laplacian of a grid. `bcs` is `[left, right, bottom, top]`; the
/// bottom/top entries are ignored when `two_d` is false.
fn build_laplacian(
    nx: usize,
    ny: usize,
    x0: f64,
    y0: f64,
    hx: f64,
    hy: f64,
    two_d: bool,
    bcs: &[Boundary; 4],
) -> Laplacian {
    let n = nx * ny;
    let mut trip = Vec::with_capacity(5 * n);
    let mut b = vec![0.0; n];
    let mut dirichlet = vec![None; n];
    for j in 0..ny {
        for i in 0..nx {
            let k = j * nx + i;
            let x = x0 + i as f64 * hx;
            let y = y0 + j as f64 * hy;
            let mut dval = None;
            if i == 0 && bcs[0].is_dirichlet() {
                dval = Some(bcs[0].value(x, y));
            } else if i == nx - 1 && bcs[1].is_dirichlet() {
                dval = Some(bcs[1].value(x, y));
            } else if two_d && j == 0 && bcs[2].is_dirichlet() {
                dval = Some(bcs[2].value(x, y));
            } else if two_d && j == ny - 1 && bcs[3].is_dirichlet() {
                dval = Some(bcs[3].value(x, y));
            }
            if dval.is_some() {
                dirichlet[k] = dval;
                continue;
            }
            let glo = if i == 0 { bcs[0].value(x, y) } else { 0.0 };
            let ghi = if i == nx - 1 { bcs[1].value(x, y) } else { 0.0 };
            let (p, c, nn, r) = axis_stencil(i, nx, hx, &bcs[0], &bcs[1], glo, ghi);
            let mut centre = c;
            if i > 0 {
                trip.push((k, k - 1, p));
            }
            if i + 1 < nx {
                trip.push((k, k + 1, nn));
            }
            b[k] += r;
            if two_d {
                let glo = if j == 0 { bcs[2].value(x, y) } else { 0.0 };
                let ghi = if j == ny - 1 { bcs[3].value(x, y) } else { 0.0 };
                let (p, c, nn, r) = axis_stencil(j, ny, hy, &bcs[2], &bcs[3], glo, ghi);
                centre += c;
                if j > 0 {
                    trip.push((k, k - nx, p));
                }
                if j + 1 < ny {
                    trip.push((k, k + nx, nn));
                }
                b[k] += r;
            }
            trip.push((k, k, centre));
        }
    }
    Laplacian { a: Csr::from_triplets(n, n, &trip), b, dirichlet }
}

fn check_steps(t_end: f64, steps: usize) -> Result<f64, PdeError> {
    if steps == 0 || !(t_end > 0.0) || !t_end.is_finite() {
        return Err(PdeError::InvalidInput("need t_end > 0 and at least one time step"));
    }
    Ok(t_end / steps as f64)
}

fn check_finite(u: &[f64]) -> Result<(), PdeError> {
    if u.iter().all(|v| v.is_finite()) { Ok(()) } else { Err(PdeError::NonFinite) }
}

fn initial_1d(grid: &Grid1D, u0: &dyn Fn(f64) -> f64, lap: &Laplacian) -> Vec<f64> {
    let h = grid.dx();
    (0..grid.n).map(|i| lap.dirichlet[i].unwrap_or_else(|| u0(grid.x0 + i as f64 * h))).collect()
}

fn initial_2d(g: &Grid2D, u0: &dyn Fn(f64, f64) -> f64, lap: &Laplacian) -> Vec<f64> {
    let (hx, hy) = (g.dx(), g.dy());
    let mut u = vec![0.0; g.nx * g.ny];
    for j in 0..g.ny {
        for i in 0..g.nx {
            let k = j * g.nx + i;
            u[k] = lap.dirichlet[k].unwrap_or_else(|| u0(g.x0 + i as f64 * hx, g.y0 + j as f64 * hy));
        }
    }
    u
}

fn bcs_1d(bc: &[Boundary; 2]) -> [Boundary; 4] {
    [bc[0].clone(), bc[1].clone(), bc[0].clone(), bc[1].clone()]
}

fn field_1d(g: &Grid1D, u: Vec<f64>, t: f64) -> Field1D {
    Field1D { x: g.nodes(), u, t }
}

fn field_2d(g: &Grid2D, u: Vec<f64>, t: f64) -> Field2D {
    Field2D { nx: g.nx, ny: g.ny, x0: g.x0, y0: g.y0, dx: g.dx(), dy: g.dy(), u, t }
}

// ---------------------------------------------------------------------------
// Heat equation (Crank-Nicolson)
// ---------------------------------------------------------------------------

fn crank_nicolson(lap: &Laplacian, alpha: f64, dt: f64, steps: usize, mut u: Vec<f64>) -> Result<Vec<f64>, PdeError> {
    let n = u.len();
    let mut lhs = Vec::with_capacity(lap.a.values.len() + n);
    let mut rhs_op = Vec::with_capacity(lap.a.values.len() + n);
    for i in 0..n {
        lhs.push((i, i, 1.0));
        rhs_op.push((i, i, 1.0));
        for q in lap.a.indptr[i]..lap.a.indptr[i + 1] {
            let v = 0.5 * dt * alpha * lap.a.values[q];
            lhs.push((i, lap.a.indices[q], -v));
            rhs_op.push((i, lap.a.indices[q], v));
        }
    }
    let m1 = Csr::from_triplets(n, n, &lhs);
    let m2 = Csr::from_triplets(n, n, &rhs_op);
    let lu = m1.sparse_lu(SparseOrdering::Rcm).map_err(|_| PdeError::SolveFailed)?;
    let mut rhs = vec![0.0; n];
    for _ in 0..steps {
        m2.matvec(&u, &mut rhs);
        for i in 0..n {
            rhs[i] += dt * alpha * lap.b[i];
        }
        u = lu.solve(&rhs);
    }
    check_finite(&u)?;
    Ok(u)
}

/// Solves `u_t = alpha u_xx` on `grid` with Crank-Nicolson time stepping
/// (`steps` equal steps up to `t_end`); `bc = [left, right]`.
///
/// The scheme is second-order in space and time and unconditionally stable.
///
/// # Errors
/// [`PdeError::InvalidInput`] for a bad grid, step count or `alpha <= 0`;
/// [`PdeError::SolveFailed`] or [`PdeError::NonFinite`] otherwise.
pub fn heat_1d(
    alpha: f64,
    grid: &Grid1D,
    bc: &[Boundary; 2],
    u0: &dyn Fn(f64) -> f64,
    t_end: f64,
    steps: usize,
) -> Result<Field1D, PdeError> {
    grid.validate(3)?;
    if !(alpha > 0.0) {
        return Err(PdeError::InvalidInput("diffusivity must be positive"));
    }
    let dt = check_steps(t_end, steps)?;
    let lap = build_laplacian(grid.n, 1, grid.x0, 0.0, grid.dx(), 1.0, false, &bcs_1d(bc));
    let u = crank_nicolson(&lap, alpha, dt, steps, initial_1d(grid, u0, &lap))?;
    Ok(field_1d(grid, u, t_end))
}

/// Solves `u_t = alpha (u_xx + u_yy)` with Crank-Nicolson;
/// `bc = [left, right, bottom, top]`.
///
/// # Errors
/// As [`heat_1d`].
pub fn heat_2d(
    alpha: f64,
    grid: &Grid2D,
    bc: &[Boundary; 4],
    u0: &dyn Fn(f64, f64) -> f64,
    t_end: f64,
    steps: usize,
) -> Result<Field2D, PdeError> {
    grid.validate()?;
    if !(alpha > 0.0) {
        return Err(PdeError::InvalidInput("diffusivity must be positive"));
    }
    let dt = check_steps(t_end, steps)?;
    let lap = build_laplacian(grid.nx, grid.ny, grid.x0, grid.y0, grid.dx(), grid.dy(), true, bc);
    let u = crank_nicolson(&lap, alpha, dt, steps, initial_2d(grid, u0, &lap))?;
    Ok(field_2d(grid, u, t_end))
}

// ---------------------------------------------------------------------------
// Wave equation (leapfrog)
// ---------------------------------------------------------------------------

fn leapfrog(lap: &Laplacian, c: f64, dt: f64, steps: usize, u0: Vec<f64>, v0: &[f64]) -> Result<Vec<f64>, PdeError> {
    let n = u0.len();
    let c2dt2 = c * c * dt * dt;
    let mut lu = vec![0.0; n];
    lap.a.matvec(&u0, &mut lu);
    // second-order start: u1 = u0 + dt v0 + dt^2/2 c^2 (L u0 + b)
    let mut prev = u0.clone();
    let mut cur = vec![0.0; n];
    for i in 0..n {
        cur[i] = lap.dirichlet[i].unwrap_or(u0[i] + dt * v0[i] + 0.5 * c2dt2 * (lu[i] + lap.b[i]));
    }
    for _ in 1..steps {
        lap.a.matvec(&cur, &mut lu);
        let mut next = vec![0.0; n];
        for i in 0..n {
            next[i] = lap.dirichlet[i].unwrap_or(2.0 * cur[i] - prev[i] + c2dt2 * (lu[i] + lap.b[i]));
        }
        prev = cur;
        cur = next;
    }
    check_finite(&cur)?;
    Ok(cur)
}

fn cfl_check(courant: f64, limit: f64) -> Result<(), PdeError> {
    if courant > limit * (1.0 + 1e-12) { Err(PdeError::Unstable { courant, limit }) } else { Ok(()) }
}

/// Solves `u_tt = c^2 u_xx` with the explicit leapfrog scheme; `v0` is the
/// initial velocity `u_t(x, 0)`.
///
/// The scheme is second-order and stable for the Courant number
/// `c dt / dx <= 1`, which is checked up front.
///
/// # Errors
/// [`PdeError::Unstable`] when the CFL condition fails, or the errors of
/// [`heat_1d`].
pub fn wave_1d(
    c: f64,
    grid: &Grid1D,
    bc: &[Boundary; 2],
    u0: &dyn Fn(f64) -> f64,
    v0: &dyn Fn(f64) -> f64,
    t_end: f64,
    steps: usize,
) -> Result<Field1D, PdeError> {
    grid.validate(3)?;
    if !(c > 0.0) {
        return Err(PdeError::InvalidInput("wave speed must be positive"));
    }
    let dt = check_steps(t_end, steps)?;
    cfl_check(c * dt / grid.dx(), 1.0)?;
    let lap = build_laplacian(grid.n, 1, grid.x0, 0.0, grid.dx(), 1.0, false, &bcs_1d(bc));
    let h = grid.dx();
    let vel: Vec<f64> = (0..grid.n).map(|i| v0(grid.x0 + i as f64 * h)).collect();
    let u = leapfrog(&lap, c, dt, steps, initial_1d(grid, u0, &lap), &vel)?;
    Ok(field_1d(grid, u, t_end))
}

/// Solves `u_tt = c^2 (u_xx + u_yy)` with leapfrog; stable for
/// `c dt sqrt(1/dx^2 + 1/dy^2) <= 1`.
///
/// # Errors
/// As [`wave_1d`].
pub fn wave_2d(
    c: f64,
    grid: &Grid2D,
    bc: &[Boundary; 4],
    u0: &dyn Fn(f64, f64) -> f64,
    v0: &dyn Fn(f64, f64) -> f64,
    t_end: f64,
    steps: usize,
) -> Result<Field2D, PdeError> {
    grid.validate()?;
    if !(c > 0.0) {
        return Err(PdeError::InvalidInput("wave speed must be positive"));
    }
    let dt = check_steps(t_end, steps)?;
    let (hx, hy) = (grid.dx(), grid.dy());
    cfl_check(c * dt * (1.0 / (hx * hx) + 1.0 / (hy * hy)).sqrt(), 1.0)?;
    let lap = build_laplacian(grid.nx, grid.ny, grid.x0, grid.y0, hx, hy, true, bc);
    let mut vel = vec![0.0; grid.nx * grid.ny];
    for j in 0..grid.ny {
        for i in 0..grid.nx {
            vel[j * grid.nx + i] = v0(grid.x0 + i as f64 * hx, grid.y0 + j as f64 * hy);
        }
    }
    let u = leapfrog(&lap, c, dt, steps, initial_2d(grid, u0, &lap), &vel)?;
    Ok(field_2d(grid, u, t_end))
}

// ---------------------------------------------------------------------------
// Poisson / Laplace
// ---------------------------------------------------------------------------

/// Linear solver for [`poisson_2d`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PoissonSolver {
    /// Sparse LU with reverse Cuthill-McKee ordering (any boundary mix).
    SparseLu,
    /// Conjugate gradients on the symmetric interior system with Jacobi
    /// preconditioning (Dirichlet sides only), to the given relative
    /// tolerance.
    Cg {
        /// Relative residual tolerance.
        tol: f64,
    },
}

/// Solves `-(u_xx + u_yy) = f` with the 5-point stencil;
/// `bc = [left, right, bottom, top]`.
///
/// A problem with Neumann conditions on all four sides is singular and is
/// rejected with [`PdeError::Unsupported`].
///
/// # Errors
/// [`PdeError::InvalidInput`], [`PdeError::Unsupported`] (pure Neumann, or CG
/// with a Neumann side) or [`PdeError::SolveFailed`] (singular matrix or CG
/// not converged).
pub fn poisson_2d(
    grid: &Grid2D,
    bc: &[Boundary; 4],
    f: &dyn Fn(f64, f64) -> f64,
    solver: PoissonSolver,
) -> Result<Field2D, PdeError> {
    grid.validate()?;
    let n_dir = bc.iter().filter(|b| b.is_dirichlet()).count();
    if n_dir == 0 {
        return Err(PdeError::Unsupported("pure Neumann Poisson problem is singular"));
    }
    let (nx, ny) = (grid.nx, grid.ny);
    let (hx, hy) = (grid.dx(), grid.dy());
    let lap = build_laplacian(nx, ny, grid.x0, grid.y0, hx, hy, true, bc);
    let n = nx * ny;
    // rhs of -L u = f + b for free nodes, g for Dirichlet nodes
    let rhs: Vec<f64> = (0..n)
        .map(|k| {
            lap.dirichlet[k].unwrap_or_else(|| {
                let (i, j) = (k % nx, k / nx);
                f(grid.x0 + i as f64 * hx, grid.y0 + j as f64 * hy) + lap.b[k]
            })
        })
        .collect();
    let u = match solver {
        PoissonSolver::SparseLu => {
            let mut trip = Vec::with_capacity(lap.a.values.len() + n);
            for i in 0..n {
                if lap.dirichlet[i].is_some() {
                    trip.push((i, i, 1.0));
                } else {
                    for q in lap.a.indptr[i]..lap.a.indptr[i + 1] {
                        trip.push((i, lap.a.indices[q], -lap.a.values[q]));
                    }
                }
            }
            let m = Csr::from_triplets(n, n, &trip);
            let lu = m.sparse_lu(SparseOrdering::Rcm).map_err(|_| PdeError::SolveFailed)?;
            lu.solve(&rhs)
        }
        PoissonSolver::Cg { tol } => {
            if n_dir != 4 {
                return Err(PdeError::Unsupported("CG needs Dirichlet conditions on all sides"));
            }
            solve_poisson_cg(&lap, &rhs, tol)?
        }
    };
    check_finite(&u)?;
    Ok(field_2d(grid, u, 0.0))
}

/// Solves `Laplace(u) = 0` (see [`poisson_2d`]).
///
/// # Errors
/// As [`poisson_2d`].
pub fn laplace_2d(grid: &Grid2D, bc: &[Boundary; 4], solver: PoissonSolver) -> Result<Field2D, PdeError> {
    poisson_2d(grid, bc, &|_, _| 0.0, solver)
}

fn solve_poisson_cg(lap: &Laplacian, rhs: &[f64], tol: f64) -> Result<Vec<f64>, PdeError> {
    let n = rhs.len();
    let free: Vec<usize> = (0..n).filter(|&k| lap.dirichlet[k].is_none()).collect();
    let mut idx = vec![usize::MAX; n];
    for (m, &k) in free.iter().enumerate() {
        idx[k] = m;
    }
    let mut trip = Vec::new();
    let mut b = vec![0.0; free.len()];
    for (m, &k) in free.iter().enumerate() {
        b[m] = rhs[k];
        for q in lap.a.indptr[k]..lap.a.indptr[k + 1] {
            let c = lap.a.indices[q];
            let v = -lap.a.values[q];
            match lap.dirichlet[c] {
                Some(g) => b[m] -= v * g,
                None => trip.push((m, idx[c], v)),
            }
        }
    }
    let a = Csr::from_triplets(free.len(), free.len(), &trip);
    let dinv = a.jacobi_inverse();
    let res = cg(
        |x, y| a.matvec(x, y),
        |r, z| {
            for i in 0..r.len() {
                z[i] = dinv[i] * r[i];
            }
        },
        &b,
        None,
        tol,
        10 * free.len() + 100,
    );
    if !res.converged {
        return Err(PdeError::SolveFailed);
    }
    let mut u = vec![0.0; n];
    for k in 0..n {
        u[k] = match lap.dirichlet[k] {
            Some(g) => g,
            None => res.x[idx[k]],
        };
    }
    Ok(u)
}

// ---------------------------------------------------------------------------
// Advection
// ---------------------------------------------------------------------------

/// Finite-difference scheme for [`advection_1d`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdvectionScheme {
    /// First-order upwind (monotone, dissipative).
    Upwind,
    /// Second-order Lax-Wendroff (dispersive near discontinuities).
    LaxWendroff,
}

/// Solves `u_t + a u_x = 0` on the periodic domain `[x0, x1)` with `n` nodes
/// (`x1` is identified with `x0`). Stable for `|a| dt / dx <= 1`, checked
/// up front.
///
/// # Errors
/// [`PdeError::Unstable`], [`PdeError::InvalidInput`] or
/// [`PdeError::NonFinite`].
pub fn advection_1d(
    a: f64,
    grid: &Grid1D,
    scheme: AdvectionScheme,
    u0: &dyn Fn(f64) -> f64,
    t_end: f64,
    steps: usize,
) -> Result<Field1D, PdeError> {
    grid.validate(3)?;
    let dt = check_steps(t_end, steps)?;
    let n = grid.n;
    let h = (grid.x1 - grid.x0) / n as f64;
    let nu = a * dt / h;
    cfl_check(nu.abs(), 1.0)?;
    let mut u: Vec<f64> = (0..n).map(|i| u0(grid.x0 + i as f64 * h)).collect();
    let mut w = vec![0.0; n];
    for _ in 0..steps {
        for i in 0..n {
            let um = u[(i + n - 1) % n];
            let up = u[(i + 1) % n];
            w[i] = match scheme {
                AdvectionScheme::Upwind => {
                    if nu >= 0.0 { u[i] - nu * (u[i] - um) } else { u[i] - nu * (up - u[i]) }
                }
                AdvectionScheme::LaxWendroff => u[i] - 0.5 * nu * (up - um) + 0.5 * nu * nu * (up - 2.0 * u[i] + um),
            };
        }
        std::mem::swap(&mut u, &mut w);
    }
    check_finite(&u)?;
    let x = (0..n).map(|i| grid.x0 + i as f64 * h).collect();
    Ok(Field1D { x, u, t: t_end })
}

// ---------------------------------------------------------------------------
// Method of lines
// ---------------------------------------------------------------------------

/// Stiff integrator used by [`method_of_lines_1d`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StiffMethod {
    /// Radau IIA(5).
    Radau5,
    /// Variable-order BDF with the given maximum order (1 to 5).
    Bdf(usize),
}

/// Result of [`method_of_lines_1d`].
#[derive(Debug, Clone, PartialEq)]
pub struct MolSolution {
    /// Final field.
    pub field: Field1D,
    /// Number of accepted ODE steps.
    pub steps: usize,
    /// Number of right-hand-side evaluations.
    pub nfev: usize,
}

/// Semi-discretises `u_t = F(t, x, u, u_x, u_xx)` on `grid` with second-order
/// central differences (ghost nodes for Neumann ends, Dirichlet nodes held
/// fixed) and integrates the ODE system to `t_end` with a stiff solver.
///
/// # Errors
/// [`PdeError::InvalidInput`] for a bad grid or `t_end`, [`PdeError::Ode`] if
/// the integrator fails.
pub fn method_of_lines_1d(
    f: &dyn Fn(f64, f64, f64, f64, f64) -> f64,
    grid: &Grid1D,
    bc: &[Boundary; 2],
    u0: &dyn Fn(f64) -> f64,
    t_end: f64,
    method: StiffMethod,
    opts: &OdeOptions,
) -> Result<MolSolution, PdeError> {
    grid.validate(3)?;
    if !(t_end > 0.0) || !t_end.is_finite() {
        return Err(PdeError::InvalidInput("t_end must be positive"));
    }
    let n = grid.n;
    let h = grid.dx();
    let x = grid.nodes();
    let dl = match &bc[0] {
        Boundary::Dirichlet(g) => Some(g(x[0], 0.0)),
        Boundary::Neumann(_) => None,
    };
    let dr = match &bc[1] {
        Boundary::Dirichlet(g) => Some(g(x[n - 1], 0.0)),
        Boundary::Neumann(_) => None,
    };
    let gl = if dl.is_none() { bc[0].value(x[0], 0.0) } else { 0.0 };
    let gr = if dr.is_none() { bc[1].value(x[n - 1], 0.0) } else { 0.0 };
    let first = usize::from(dl.is_some());
    let last = if dr.is_some() { n - 1 } else { n };
    let m = last - first;
    let expand = |y: &[f64]| -> Vec<f64> {
        let mut u = vec![0.0; n];
        if let Some(v) = dl {
            u[0] = v;
        }
        if let Some(v) = dr {
            u[n - 1] = v;
        }
        u[first..last].copy_from_slice(y);
        u
    };
    let rhs = |t: f64, y: &[f64], dy: &mut [f64]| {
        let u = expand(y);
        for i in first..last {
            let um = if i == 0 { u[1] + 2.0 * h * gl } else { u[i - 1] };
            let up = if i == n - 1 { u[n - 2] + 2.0 * h * gr } else { u[i + 1] };
            let ux = (up - um) / (2.0 * h);
            let uxx = (up - 2.0 * u[i] + um) / (h * h);
            dy[i - first] = f(t, x[i], u[i], ux, uxx);
        }
    };
    let y0: Vec<f64> = (first..last).map(|i| u0(x[i])).collect();
    let sol = match method {
        StiffMethod::Radau5 => radau5(rhs, 0.0, t_end, &y0, opts)?,
        StiffMethod::Bdf(order) => bdf(rhs, 0.0, t_end, &y0, opts, order.clamp(1, 5))?,
    };
    let yl = sol.y.last().ok_or(PdeError::SolveFailed)?;
    if yl.len() != m {
        return Err(PdeError::SolveFailed);
    }
    let u = expand(yl);
    check_finite(&u)?;
    Ok(MolSolution { field: Field1D { x, u, t: t_end }, steps: sol.t.len().saturating_sub(1), nfev: sol.nfev })
}

// ---------------------------------------------------------------------------
// Dispatcher
// ---------------------------------------------------------------------------

/// A PDE problem for [`pde_solver`].
pub enum PdeProblem<'a> {
    /// 1D heat equation, see [`heat_1d`].
    Heat1D {
        /// Diffusivity.
        alpha: f64,
        /// Grid.
        grid: Grid1D,
        /// `[left, right]` boundaries.
        bc: [Boundary; 2],
        /// Initial condition.
        u0: &'a dyn Fn(f64) -> f64,
        /// Final time.
        t_end: f64,
        /// Number of time steps.
        steps: usize,
    },
    /// 2D heat equation, see [`heat_2d`].
    Heat2D {
        /// Diffusivity.
        alpha: f64,
        /// Grid.
        grid: Grid2D,
        /// `[left, right, bottom, top]` boundaries.
        bc: [Boundary; 4],
        /// Initial condition.
        u0: &'a dyn Fn(f64, f64) -> f64,
        /// Final time.
        t_end: f64,
        /// Number of time steps.
        steps: usize,
    },
    /// 1D wave equation, see [`wave_1d`].
    Wave1D {
        /// Wave speed.
        c: f64,
        /// Grid.
        grid: Grid1D,
        /// `[left, right]` boundaries.
        bc: [Boundary; 2],
        /// Initial displacement.
        u0: &'a dyn Fn(f64) -> f64,
        /// Initial velocity.
        v0: &'a dyn Fn(f64) -> f64,
        /// Final time.
        t_end: f64,
        /// Number of time steps.
        steps: usize,
    },
    /// 2D wave equation, see [`wave_2d`].
    Wave2D {
        /// Wave speed.
        c: f64,
        /// Grid.
        grid: Grid2D,
        /// `[left, right, bottom, top]` boundaries.
        bc: [Boundary; 4],
        /// Initial displacement.
        u0: &'a dyn Fn(f64, f64) -> f64,
        /// Initial velocity.
        v0: &'a dyn Fn(f64, f64) -> f64,
        /// Final time.
        t_end: f64,
        /// Number of time steps.
        steps: usize,
    },
    /// 2D Poisson equation `-Laplace(u) = f`, see [`poisson_2d`].
    Poisson2D {
        /// Grid.
        grid: Grid2D,
        /// `[left, right, bottom, top]` boundaries.
        bc: [Boundary; 4],
        /// Source term.
        f: &'a dyn Fn(f64, f64) -> f64,
        /// Linear solver.
        solver: PoissonSolver,
    },
    /// 1D advection, see [`advection_1d`].
    Advection1D {
        /// Advection speed.
        a: f64,
        /// Periodic grid.
        grid: Grid1D,
        /// Scheme.
        scheme: AdvectionScheme,
        /// Initial condition.
        u0: &'a dyn Fn(f64) -> f64,
        /// Final time.
        t_end: f64,
        /// Number of time steps.
        steps: usize,
    },
}

/// Solution returned by [`pde_solver`].
#[derive(Debug, Clone, PartialEq)]
pub enum PdeSolution {
    /// Field on a 1D grid.
    D1(Field1D),
    /// Field on a 2D grid.
    D2(Field2D),
}

/// Solves a [`PdeProblem`] with the matching finite-difference scheme.
///
/// # Errors
/// Any [`PdeError`] of the underlying solver.
pub fn pde_solver(problem: &PdeProblem<'_>) -> Result<PdeSolution, PdeError> {
    match problem {
        PdeProblem::Heat1D { alpha, grid, bc, u0, t_end, steps } => {
            heat_1d(*alpha, grid, bc, *u0, *t_end, *steps).map(PdeSolution::D1)
        }
        PdeProblem::Heat2D { alpha, grid, bc, u0, t_end, steps } => {
            heat_2d(*alpha, grid, bc, *u0, *t_end, *steps).map(PdeSolution::D2)
        }
        PdeProblem::Wave1D { c, grid, bc, u0, v0, t_end, steps } => {
            wave_1d(*c, grid, bc, *u0, *v0, *t_end, *steps).map(PdeSolution::D1)
        }
        PdeProblem::Wave2D { c, grid, bc, u0, v0, t_end, steps } => {
            wave_2d(*c, grid, bc, *u0, *v0, *t_end, *steps).map(PdeSolution::D2)
        }
        PdeProblem::Poisson2D { grid, bc, f, solver } => poisson_2d(grid, bc, *f, *solver).map(PdeSolution::D2),
        PdeProblem::Advection1D { a, grid, scheme, u0, t_end, steps } => {
            advection_1d(*a, grid, *scheme, *u0, *t_end, *steps).map(PdeSolution::D1)
        }
    }
}
