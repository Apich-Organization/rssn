//! # Structured-grid finite differences in 1D, 2D and 3D
//!
//! A dimension-generic vertex-centred finite-difference layer shared by the
//! FDM, multigrid, Crank–Nicolson and wave solvers. Boundary conditions come
//! from [`crate::sim::boundary`]; Neumann and Robin data are imposed with a
//! ghost node so that every scheme stays second order.
//!
//! | Capability | 1D | 2D | 3D | Boundary conditions |
//! |---|---|---|---|---|
//! | Laplacian operator ([`LaplacianOp`]) | yes | yes | yes | Dirichlet, Neumann, Robin, periodic, symmetry; per face, time dependent |
//! | Poisson / Helmholtz `-Δu + σu = f` ([`solve_helmholtz`]) | yes | yes | yes | all of the above; LU, ILU-BiCGStab or GMRES |
//! | Heat `u_t = αΔu + s` θ-method ([`heat_fixed`]) | yes | yes | yes | all of the above |
//! | Heat with step-doubling error control ([`heat_adaptive`]) | yes | yes | yes | all of the above |
//! | Wave `u_tt = c²Δu + s` leapfrog ([`wave_leapfrog`]) | yes | yes | yes | Dirichlet, Neumann, Robin, periodic, first-order Mur absorbing |
//! | Nonlinear diffusion `u_t = ∇·(D(u)∇u) + s` ([`nonlinear_diffusion`]) | yes | yes | yes | Dirichlet, Neumann, Robin, periodic; Newton–Krylov (matrix-free GMRES) implicit stepping, optional adaptive steps |
//!
//! Grids are vertex centred: `n` nodes per axis including both boundary
//! nodes (spacing `(hi - lo)/(n - 1)`); a periodic axis has `n` distinct
//! nodes with spacing `(hi - lo)/n`.
//!
//! Limitations: constant-coefficient Laplacian only (variable diffusion is
//! provided by [`nonlinear_diffusion`]); the scalar wave
//! equation with a PML is in [`crate::sim::acoustics`]; the Mur condition
//! here is first order.

use crate::kernels::krylov::{Csr, Ilu0, SparseLu, SparseOrdering, bicgstab, gmres, identity_precond};
use crate::sim::boundary::{BcValue, BoundaryCondition, Face, FaceBcs};
use std::collections::HashMap;
use std::sync::Arc;

/// Vertex-centred structured grid with 1 to 3 axes.
#[derive(Clone, Debug)]
pub struct Grid {
    /// Number of active axes (1, 2 or 3).
    pub dim: usize,
    /// Nodes per axis (1 for inactive axes).
    pub n: [usize; 3],
    /// Lower corner.
    pub lo: [f64; 3],
    /// Upper corner.
    pub hi: [f64; 3],
    /// Periodic flags per axis.
    pub periodic: [bool; 3],
}

impl Grid {
    /// Builds a grid; periodicity is read from `bcs`.
    ///
    /// # Errors
    /// For a bad dimension, too few nodes, `lo >= hi`, or one-sided periodicity.
    pub fn new(n: &[usize], lo: &[f64], hi: &[f64], bcs: &FaceBcs) -> Result<Self, String> {
        let dim = n.len();
        if !(1..=3).contains(&dim) || lo.len() != dim || hi.len() != dim {
            return Err("grid needs 1 to 3 axes with matching bounds".into());
        }
        let mut g = Self { dim, n: [1; 3], lo: [0.0; 3], hi: [1.0; 3], periodic: [false; 3] };
        for a in 0..dim {
            let per = bcs.axis_periodic(a)?;
            if n[a] < if per { 3 } else { 2 } || lo[a] >= hi[a] {
                return Err(format!("axis {a}: invalid node count or bounds"));
            }
            g.n[a] = n[a];
            g.lo[a] = lo[a];
            g.hi[a] = hi[a];
            g.periodic[a] = per;
        }
        Ok(g)
    }

    /// Unit hypercube `[0,1]^dim` with `n` nodes per axis.
    ///
    /// # Errors
    /// See [`Grid::new`].
    pub fn unit(dim: usize, n: usize, bcs: &FaceBcs) -> Result<Self, String> {
        Self::new(&vec![n; dim], &vec![0.0; dim], &vec![1.0; dim], bcs)
    }

    /// Spacing along `axis`.
    #[must_use]
    pub fn h(&self, axis: usize) -> f64 {
        let span = self.hi[axis] - self.lo[axis];
        if self.periodic[axis] { span / self.n[axis] as f64 } else { span / (self.n[axis] - 1) as f64 }
    }

    /// Smallest spacing over the active axes.
    #[must_use]
    pub fn h_min(&self) -> f64 {
        (0..self.dim).map(|a| self.h(a)).fold(f64::INFINITY, f64::min)
    }

    /// Number of nodes.
    #[must_use]
    pub fn len(&self) -> usize {
        self.n.iter().product()
    }

    /// Whether the grid has no nodes (never true for a valid grid).
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Flat index of node `(i, j, k)` (x fastest).
    #[must_use]
    pub const fn index(&self, i: [usize; 3]) -> usize {
        i[0] + self.n[0] * (i[1] + self.n[1] * i[2])
    }

    /// Node multi-index of a flat index.
    #[must_use]
    pub const fn ijk(&self, p: usize) -> [usize; 3] {
        let i = p % self.n[0];
        let r = p / self.n[0];
        [i, r % self.n[1], r / self.n[1]]
    }

    /// Coordinates of node `p` (inactive axes are `0`).
    #[must_use]
    pub fn coord(&self, p: usize) -> [f64; 3] {
        let ijk = self.ijk(p);
        let mut x = [0.0; 3];
        for a in 0..self.dim {
            x[a] = self.lo[a] + ijk[a] as f64 * self.h(a);
        }
        x
    }

    /// Samples `f` at every node.
    pub fn sample(&self, f: impl Fn([f64; 3]) -> f64) -> Vec<f64> {
        (0..self.len()).map(|p| f(self.coord(p))).collect()
    }

    /// Largest absolute difference between `u` and `exact` at the nodes.
    pub fn max_error(&self, u: &[f64], exact: impl Fn([f64; 3]) -> f64) -> f64 {
        u.iter().enumerate().map(|(p, v)| (v - exact(self.coord(p))).abs()).fold(0.0, f64::max)
    }

    /// Discrete L2 norm of the error (`sqrt(h^d Σ e²)`).
    pub fn l2_error(&self, u: &[f64], exact: impl Fn([f64; 3]) -> f64) -> f64 {
        let vol: f64 = (0..self.dim).map(|a| self.h(a)).product();
        (u.iter().enumerate().map(|(p, v)| (v - exact(self.coord(p))).powi(2)).sum::<f64>() * vol).sqrt()
    }
}

/// Dirichlet node data.
#[derive(Clone, Debug)]
struct DirNode {
    value: BcValue,
    x: [f64; 3],
}

/// The discrete Laplacian `L` with boundary conditions folded in.
///
/// For a free node `p`, `(L u)_p + forcing_p(t)` approximates `Δu`; Dirichlet
/// nodes carry `u = g(t, x)` and have no row.
#[derive(Clone, Debug)]
pub struct LaplacianOp {
    /// Grid the operator lives on.
    pub grid: Grid,
    /// Matrix of `L` without Dirichlet rows (those are empty).
    pub l: Csr,
    dirichlet: Vec<Option<DirNode>>,
    forcing: Vec<(usize, f64, BcValue, [f64; 3])>,
    absorbing: Vec<(usize, usize, bool)>,
}

impl LaplacianOp {
    /// Assembles the operator for the given boundary conditions.
    ///
    /// # Errors
    /// If a Robin condition is degenerate (`alpha = beta = 0`).
    pub fn new(grid: &Grid, bcs: &FaceBcs) -> Result<Self, String> {
        let nn = grid.len();
        let mut dirichlet: Vec<Option<DirNode>> = vec![None; nn];
        // Pass 1: Dirichlet nodes.
        for p in 0..nn {
            let ijk = grid.ijk(p);
            for a in 0..grid.dim {
                if grid.periodic[a] {
                    continue;
                }
                for max in [false, true] {
                    let on = if max { ijk[a] == grid.n[a] - 1 } else { ijk[a] == 0 };
                    if !on || dirichlet[p].is_some() {
                        continue;
                    }
                    let face = Face::of_axis(a, max);
                    match bcs.get(face) {
                        BoundaryCondition::Dirichlet(g) => {
                            dirichlet[p] = Some(DirNode { value: g.clone(), x: grid.coord(p) });
                        },
                        BoundaryCondition::Robin { alpha, beta, g } if *beta == 0.0 => {
                            if *alpha == 0.0 {
                                return Err("degenerate Robin condition".into());
                            }
                            let (al, g) = (*alpha, g.clone());
                            dirichlet[p] = Some(DirNode {
                                value: BcValue::spacetime(move |t, x| g.eval(t, x) / al),
                                x: grid.coord(p),
                            });
                        },
                        _ => {},
                    }
                }
            }
        }
        let mut trip: Vec<(usize, usize, f64)> = Vec::new();
        let mut forcing = Vec::new();
        let mut absorbing = Vec::new();
        for p in 0..nn {
            if dirichlet[p].is_some() {
                continue;
            }
            let ijk = grid.ijk(p);
            let x = grid.coord(p);
            for a in 0..grid.dim {
                let h = grid.h(a);
                let h2 = h * h;
                let n = grid.n[a];
                let i = ijk[a];
                let at = |ii: usize| {
                    let mut q = ijk;
                    q[a] = ii;
                    grid.index(q)
                };
                if grid.periodic[a] {
                    let (im, ip) = ((i + n - 1) % n, (i + 1) % n);
                    trip.push((p, at(im), 1.0 / h2));
                    trip.push((p, at(ip), 1.0 / h2));
                    trip.push((p, p, -2.0 / h2));
                    continue;
                }
                if i > 0 && i < n - 1 {
                    trip.push((p, at(i - 1), 1.0 / h2));
                    trip.push((p, at(i + 1), 1.0 / h2));
                    trip.push((p, p, -2.0 / h2));
                    continue;
                }
                let max = i == n - 1;
                let inner = if max { at(n - 2) } else { at(1) };
                let bc = bcs.get(Face::of_axis(a, max));
                // ghost node: u_gh = u_in + 2h dudn
                trip.push((p, inner, 2.0 / h2));
                trip.push((p, p, -2.0 / h2));
                match bc {
                    BoundaryCondition::Neumann(g) => forcing.push((p, 2.0 / h, g.clone(), x)),
                    BoundaryCondition::Symmetry => {},
                    BoundaryCondition::Absorbing => absorbing.push((p, a, max)),
                    BoundaryCondition::Robin { alpha, beta, g } => {
                        trip.push((p, p, -2.0 * alpha / (beta * h)));
                        forcing.push((p, 2.0 / (beta * h), g.clone(), x));
                    },
                    BoundaryCondition::Dirichlet(_) | BoundaryCondition::Periodic => {
                        return Err("inconsistent boundary description".into());
                    },
                }
            }
        }
        Ok(Self { grid: grid.clone(), l: Csr::from_triplets(nn, nn, &trip), dirichlet, forcing, absorbing })
    }

    /// Whether any face carries an absorbing (wave-only) condition.
    #[must_use]
    pub const fn has_absorbing(&self) -> bool {
        !self.absorbing.is_empty()
    }

    /// Whether node `p` is a Dirichlet node.
    #[must_use]
    pub fn is_dirichlet(&self, p: usize) -> bool {
        self.dirichlet.get(p).is_some_and(Option::is_some)
    }

    /// Number of Dirichlet nodes.
    #[must_use]
    pub fn dirichlet_count(&self) -> usize {
        self.dirichlet.iter().filter(|d| d.is_some()).count()
    }

    /// `out = L u` at free nodes, `0` at Dirichlet nodes.
    pub fn apply(&self, u: &[f64], out: &mut [f64]) {
        self.l.matvec(u, out);
    }

    /// Boundary-datum contribution to `Δu` for Neumann and Robin faces at time `t`.
    #[must_use]
    pub fn forcing(&self, t: f64) -> Vec<f64> {
        let mut f = vec![0.0; self.grid.len()];
        for (p, c, g, x) in &self.forcing {
            f[*p] += c * g.eval(t, *x);
        }
        f
    }

    /// Overwrites Dirichlet nodes of `u` with `g(t, x)`.
    pub fn set_dirichlet(&self, u: &mut [f64], t: f64) {
        for (p, d) in self.dirichlet.iter().enumerate() {
            if let Some(d) = d {
                u[p] = d.value.eval(t, d.x);
            }
        }
    }

    /// Matrix `shift I + scale L` with identity rows at Dirichlet nodes.
    #[must_use]
    pub fn system(&self, shift: f64, scale: f64) -> Csr {
        let nn = self.grid.len();
        let mut trip = Vec::with_capacity(self.l.values.len() + nn);
        for p in 0..nn {
            if self.is_dirichlet(p) {
                trip.push((p, p, 1.0));
                continue;
            }
            trip.push((p, p, shift));
            for k in self.l.indptr[p]..self.l.indptr[p + 1] {
                trip.push((p, self.l.indices[k], scale * self.l.values[k]));
            }
        }
        Csr::from_triplets(nn, nn, &trip)
    }
}

/// Choice of linear solver for the sparse systems of this module.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum LinearSolver {
    /// Sparse LU for small systems, ILU(0)-`BiCGStab` for large ones.
    Auto,
    /// Sparse LU (minimum-degree ordering).
    Direct,
    /// ILU(0)-preconditioned `BiCGStab` with the given relative tolerance.
    BiCgStab(f64),
    /// ILU(0)-preconditioned restarted GMRES(30) with the given tolerance.
    Gmres(f64),
}

/// A factored or preconditioned linear system.
pub struct LinearSystem {
    a: Csr,
    lu: Option<SparseLu>,
    ilu: Option<Ilu0>,
    gmres: bool,
    tol: f64,
}

impl LinearSystem {
    /// Prepares `a` for repeated solves.
    ///
    /// # Errors
    /// If the matrix is singular (direct solver).
    pub fn new(a: Csr, solver: LinearSolver) -> Result<Self, String> {
        let n = a.rows;
        let direct = match solver {
            LinearSolver::Direct => true,
            LinearSolver::Auto => n <= 30_000,
            _ => false,
        };
        if direct {
            let lu = a.sparse_lu(SparseOrdering::MinDegree).map_err(|e| format!("sparse LU failed: {e:?}"))?;
            return Ok(Self { a, lu: Some(lu), ilu: None, gmres: false, tol: 0.0 });
        }
        let (gmres, tol) = match solver {
            LinearSolver::Gmres(t) => (true, t),
            LinearSolver::BiCgStab(t) => (false, t),
            _ => (false, 1e-12),
        };
        let ilu = a.ilu0();
        Ok(Self { a, lu: None, ilu, gmres, tol })
    }

    /// Solves `A x = b` (`x0` is the Krylov initial guess).
    ///
    /// # Errors
    /// If an iterative solve does not reach its tolerance.
    pub fn solve(&self, b: &[f64], x0: Option<&[f64]>) -> Result<Vec<f64>, String> {
        if let Some(lu) = &self.lu {
            return Ok(lu.solve(b));
        }
        let apply = |x: &[f64], y: &mut [f64]| self.a.matvec(x, y);
        let pre = |r: &[f64], z: &mut [f64]| match &self.ilu {
            Some(m) => m.apply(r, z),
            None => identity_precond(r, z),
        };
        let res = if self.gmres {
            gmres(apply, pre, b, x0, 30, self.tol, 20 * b.len().max(100))
        } else {
            bicgstab(apply, pre, b, x0, self.tol, 20 * b.len().max(100))
        };
        if res.converged { Ok(res.x) } else { Err(format!("iterative solve stagnated at residual {}", res.residual)) }
    }
}

/// Solves `-Δu + σ u = f` on the grid with the given boundary conditions
/// (`σ` may be negative for the indefinite Helmholtz equation).
///
/// # Errors
/// On a singular problem (e.g. pure Neumann with `σ = 0`), a solver failure or
/// an absorbing boundary.
pub fn solve_helmholtz(
    grid: &Grid,
    bcs: &FaceBcs,
    sigma: f64,
    f: impl Fn([f64; 3]) -> f64,
    solver: LinearSolver,
) -> Result<Vec<f64>, String> {
    let op = LaplacianOp::new(grid, bcs)?;
    if op.has_absorbing() {
        return Err("absorbing boundaries are only valid for waves".into());
    }
    let sys = LinearSystem::new(op.system(sigma, -1.0), solver)?;
    let force = op.forcing(0.0);
    let mut b: Vec<f64> = (0..grid.len()).map(|p| f(grid.coord(p)) + force[p]).collect();
    op.set_dirichlet(&mut b, 0.0);
    sys.solve(&b, None)
}

/// Heat / diffusion problem `u_t = α Δu + s(t, x)`.
pub struct HeatProblem<'a> {
    /// Grid.
    pub grid: &'a Grid,
    /// Boundary conditions.
    pub bcs: &'a FaceBcs,
    /// Diffusivity.
    pub alpha: f64,
    /// Initial condition.
    pub init: &'a dyn Fn([f64; 3]) -> f64,
    /// Optional volumetric source.
    pub source: Option<&'a dyn Fn(f64, [f64; 3]) -> f64>,
    /// Linear solver.
    pub solver: LinearSolver,
}

fn source_vec(grid: &Grid, s: Option<&dyn Fn(f64, [f64; 3]) -> f64>, t: f64) -> Vec<f64> {
    (0..grid.len()).map(|p| s.map_or(0.0, |f| f(t, grid.coord(p)))).collect()
}

struct HeatStepper<'a> {
    prob: &'a HeatProblem<'a>,
    op: LaplacianOp,
    theta: f64,
    cache: HashMap<u64, Arc<LinearSystem>>,
}

impl<'a> HeatStepper<'a> {
    fn new(prob: &'a HeatProblem<'a>, theta: f64) -> Result<Self, String> {
        let op = LaplacianOp::new(prob.grid, prob.bcs)?;
        if op.has_absorbing() {
            return Err("absorbing boundaries are only valid for waves".into());
        }
        Ok(Self { prob, op, theta, cache: HashMap::new() })
    }

    fn step(&mut self, u: &[f64], t: f64, dt: f64) -> Result<Vec<f64>, String> {
        let nn = u.len();
        let a = self.prob.alpha;
        let th = self.theta;
        let key = dt.to_bits();
        let sys = if let Some(s) = self.cache.get(&key) {
            Arc::clone(s)
        } else {
            if self.cache.len() > 8 {
                self.cache.clear();
            }
            let s = Arc::new(LinearSystem::new(self.op.system(1.0, -dt * th * a), self.prob.solver)?);
            self.cache.insert(key, Arc::clone(&s));
            s
        };
        let mut lu = vec![0.0; nn];
        self.op.apply(u, &mut lu);
        let f0 = self.op.forcing(t);
        let f1 = self.op.forcing(t + dt);
        let s0 = source_vec(self.prob.grid, self.prob.source, t);
        let s1 = source_vec(self.prob.grid, self.prob.source, t + dt);
        let mut b = vec![0.0; nn];
        for p in 0..nn {
            b[p] = u[p]
                + dt * (1.0 - th) * (a * (lu[p] + f0[p]) + s0[p])
                + dt * th * (a * f1[p] + s1[p]);
        }
        self.op.set_dirichlet(&mut b, t + dt);
        sys.solve(&b, Some(u))
    }

    fn initial(&self) -> Vec<f64> {
        let mut u = self.prob.grid.sample(|x| (self.prob.init)(x));
        self.op.set_dirichlet(&mut u, 0.0);
        u
    }
}

/// Integrates the heat equation with the θ-method (`θ = 1` backward Euler,
/// `θ = 1/2` Crank–Nicolson) using `steps` equal steps.
///
/// # Errors
/// On solver failure or an invalid `theta` / absorbing boundary.
pub fn heat_fixed(prob: &HeatProblem<'_>, t_end: f64, steps: usize, theta: f64) -> Result<Vec<f64>, String> {
    if !(0.5..=1.0).contains(&theta) || steps == 0 {
        return Err("theta must be in [1/2, 1] and steps > 0".into());
    }
    let mut st = HeatStepper::new(prob, theta)?;
    let mut u = st.initial();
    let dt = t_end / steps as f64;
    for k in 0..steps {
        u = st.step(&u, k as f64 * dt, dt)?;
    }
    Ok(u)
}

/// Controls for adaptive step-doubling integration.
#[derive(Clone, Copy, Debug)]
pub struct AdaptiveOptions {
    /// Error tolerance per step relative to `1 + max|u|`.
    pub tol: f64,
    /// First step.
    pub dt0: f64,
    /// Smallest allowed step.
    pub dt_min: f64,
    /// Largest allowed step.
    pub dt_max: f64,
    /// Cap on attempted steps.
    pub max_steps: usize,
}

impl Default for AdaptiveOptions {
    fn default() -> Self {
        Self { tol: 1e-6, dt0: 1e-3, dt_min: 1e-12, dt_max: f64::INFINITY, max_steps: 100_000 }
    }
}

/// Statistics of an adaptive run.
#[derive(Clone, Copy, Debug, Default)]
pub struct AdaptiveStats {
    /// Accepted steps.
    pub accepted: usize,
    /// Rejected steps.
    pub rejected: usize,
    /// Smallest accepted step.
    pub dt_min_used: f64,
    /// Largest accepted step.
    pub dt_max_used: f64,
}

/// Generic step-doubling driver: `step(u, t, dt)` is one step of a method of
/// the given `order`; the difference between one step of `dt` and two of
/// `dt/2` estimates the local error.
fn adaptive_drive(
    mut u: Vec<f64>,
    t_end: f64,
    order: i32,
    opts: &AdaptiveOptions,
    mut step: impl FnMut(&[f64], f64, f64) -> Result<Vec<f64>, String>,
) -> Result<(Vec<f64>, AdaptiveStats), String> {
    let mut t = 0.0;
    let mut dt = opts.dt0.min(t_end);
    let mut stats = AdaptiveStats { dt_min_used: f64::INFINITY, ..AdaptiveStats::default() };
    let fac = f64::from(2_i32.pow(order.unsigned_abs()) - 1);
    let mut attempts = 0;
    while t < t_end * (1.0 - 1e-14) {
        attempts += 1;
        if attempts > opts.max_steps {
            return Err("adaptive integration exceeded max_steps".into());
        }
        let dt_try = dt.min(t_end - t).min(opts.dt_max);
        let big = step(&u, t, dt_try)?;
        let mid = step(&u, t, 0.5 * dt_try)?;
        let small = step(&mid, t + 0.5 * dt_try, 0.5 * dt_try)?;
        let scale = 1.0 + small.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
        let err = big.iter().zip(&small).fold(0.0_f64, |m, (a, b)| m.max((a - b).abs())) / (fac * scale);
        if err <= opts.tol || dt_try <= opts.dt_min {
            t += dt_try;
            u = small;
            stats.accepted += 1;
            stats.dt_min_used = stats.dt_min_used.min(dt_try);
            stats.dt_max_used = stats.dt_max_used.max(dt_try);
        } else {
            stats.rejected += 1;
        }
        let ratio = if err > 0.0 { 0.9 * (opts.tol / err).powf(1.0 / f64::from(order + 1)) } else { 5.0 };
        dt = (dt_try * ratio.clamp(0.2, 5.0)).max(opts.dt_min);
    }
    if stats.accepted == 0 {
        stats.dt_min_used = 0.0;
    }
    Ok((u, stats))
}

/// Crank–Nicolson heat integration with step-doubling error control.
///
/// # Errors
/// On solver failure or when `max_steps` is exceeded.
pub fn heat_adaptive(
    prob: &HeatProblem<'_>,
    t_end: f64,
    opts: &AdaptiveOptions,
) -> Result<(Vec<f64>, AdaptiveStats), String> {
    let mut st = HeatStepper::new(prob, 0.5)?;
    let u0 = st.initial();
    adaptive_drive(u0, t_end, 2, opts, |u, t, dt| st.step(u, t, dt))
}

/// Wave problem `u_tt = c² Δu + s(t, x)`.
pub struct WaveProblem<'a> {
    /// Grid.
    pub grid: &'a Grid,
    /// Boundary conditions (Dirichlet, Neumann, Robin, periodic, absorbing).
    pub bcs: &'a FaceBcs,
    /// Wave speed.
    pub c: f64,
    /// Initial displacement.
    pub u0: &'a dyn Fn([f64; 3]) -> f64,
    /// Initial velocity.
    pub v0: &'a dyn Fn([f64; 3]) -> f64,
    /// Optional source.
    pub source: Option<&'a dyn Fn(f64, [f64; 3]) -> f64>,
}

/// Leapfrog integration of the wave equation to `t_end` with `steps` steps.
///
/// # Errors
/// If the CFL condition `c dt <= h_min / sqrt(dim)` is violated.
pub fn wave_leapfrog(prob: &WaveProblem<'_>, t_end: f64, steps: usize) -> Result<Vec<f64>, String> {
    let op = LaplacianOp::new(prob.grid, prob.bcs)?;
    let g = prob.grid;
    let nn = g.len();
    let dt = t_end / steps as f64;
    let c2 = prob.c * prob.c;
    if prob.c * dt > g.h_min() / (g.dim as f64).sqrt() * (1.0 + 1e-12) {
        return Err("CFL condition violated".into());
    }
    let mut prev = g.sample(|x| (prob.u0)(x));
    op.set_dirichlet(&mut prev, 0.0);
    let v0 = g.sample(|x| (prob.v0)(x));
    let mut lu = vec![0.0; nn];
    op.apply(&prev, &mut lu);
    let f0 = op.forcing(0.0);
    let s0 = source_vec(g, prob.source, 0.0);
    let mut cur: Vec<f64> = (0..nn)
        .map(|p| prev[p] + dt * v0[p] + 0.5 * dt * dt * (c2 * (lu[p] + f0[p]) + s0[p]))
        .collect();
    mur_update(&op, prob.c, dt, &prev, &mut cur);
    op.set_dirichlet(&mut cur, dt);
    for k in 1..steps {
        let t = k as f64 * dt;
        op.apply(&cur, &mut lu);
        let f = op.forcing(t);
        let s = source_vec(g, prob.source, t);
        let mut next: Vec<f64> =
            (0..nn).map(|p| 2.0 * cur[p] - prev[p] + dt * dt * (c2 * (lu[p] + f[p]) + s[p])).collect();
        mur_update(&op, prob.c, dt, &cur, &mut next);
        op.set_dirichlet(&mut next, t + dt);
        prev = cur;
        cur = next;
    }
    Ok(cur)
}

/// First-order Mur radiation update of the absorbing boundary nodes.
fn mur_update(op: &LaplacianOp, c: f64, dt: f64, old: &[f64], new: &mut [f64]) {
    let g = &op.grid;
    for &(p, a, max) in &op.absorbing {
        let mut q = g.ijk(p);
        q[a] = if max { q[a] - 1 } else { q[a] + 1 };
        let qi = g.index(q);
        let r = (c * dt - g.h(a)) / (c * dt + g.h(a));
        new[p] = old[qi] + r * (new[qi] - old[p]);
    }
}

/// Nonlinear diffusion problem `u_t = ∇·(D(u)∇u) + s(t, x, u)`.
pub struct NonlinearDiffusion<'a> {
    /// Grid.
    pub grid: &'a Grid,
    /// Boundary conditions (Dirichlet, Neumann, Robin, periodic).
    pub bcs: &'a FaceBcs,
    /// Diffusivity `D(u)`.
    pub d: &'a dyn Fn(f64) -> f64,
    /// Optional source `s(t, x, u)`.
    pub source: Option<&'a dyn Fn(f64, [f64; 3], f64) -> f64>,
    /// Initial condition.
    pub init: &'a dyn Fn([f64; 3]) -> f64,
}

/// Newton–Krylov diagnostics.
#[derive(Clone, Copy, Debug, Default)]
pub struct NewtonStats {
    /// Newton iterations over the whole run.
    pub newton_iters: usize,
    /// Matrix-free Jacobian-vector products (GMRES iterations).
    pub krylov_iters: usize,
}

struct NlStepper<'a> {
    p: &'a NonlinearDiffusion<'a>,
    op: LaplacianOp,
    bc: Vec<[Option<(f64, f64, BcValue)>; 2]>, // per axis: (alpha, beta, g) at min / max (Neumann/Robin)
    stats: NewtonStats,
}

impl<'a> NlStepper<'a> {
    fn new(p: &'a NonlinearDiffusion<'a>) -> Result<Self, String> {
        let op = LaplacianOp::new(p.grid, p.bcs)?;
        if op.has_absorbing() {
            return Err("absorbing boundaries are only valid for waves".into());
        }
        let mut bc = Vec::new();
        for a in 0..p.grid.dim {
            let mut pair: [Option<(f64, f64, BcValue)>; 2] = [None, None];
            for (m, slot) in pair.iter_mut().enumerate() {
                *slot = match p.bcs.get(Face::of_axis(a, m == 1)) {
                    BoundaryCondition::Neumann(g) => Some((0.0, 1.0, g.clone())),
                    BoundaryCondition::Symmetry => Some((0.0, 1.0, BcValue::constant(0.0))),
                    BoundaryCondition::Robin { alpha, beta, g } if *beta != 0.0 => Some((*alpha, *beta, g.clone())),
                    _ => None,
                };
            }
            bc.push(pair);
        }
        Ok(Self { p, op, bc, stats: NewtonStats::default() })
    }

    /// Spatial right-hand side `R(u, t)` at free nodes (zero at Dirichlet nodes).
    fn rhs(&self, u: &[f64], t: f64) -> Vec<f64> {
        let g = self.p.grid;
        let nn = g.len();
        let dv: Vec<f64> = u.iter().map(|&v| (self.p.d)(v)).collect();
        let mut r = vec![0.0; nn];
        for p in 0..nn {
            if self.op.is_dirichlet(p) {
                continue;
            }
            let ijk = g.ijk(p);
            let x = g.coord(p);
            let mut acc = 0.0;
            for a in 0..g.dim {
                let h = g.h(a);
                let n = g.n[a];
                let i = ijk[a];
                let at = |ii: usize| {
                    let mut q = ijk;
                    q[a] = ii;
                    g.index(q)
                };
                if g.periodic[a] || (i > 0 && i < n - 1) {
                    let (qm, qp) = (at((i + n - 1) % n), at((i + 1) % n));
                    let dp = f64::midpoint(dv[p], dv[qp]);
                    let dm = f64::midpoint(dv[p], dv[qm]);
                    acc += (dp * (u[qp] - u[p]) - dm * (u[p] - u[qm])) / (h * h);
                } else {
                    let max = i == n - 1;
                    let qi = if max { at(n - 2) } else { at(1) };
                    let din = f64::midpoint(dv[p], dv[qi]);
                    let dudn = match &self.bc[a][usize::from(max)] {
                        Some((al, be, gv)) => (gv.eval(t, x) - al * u[p]) / be,
                        None => 0.0,
                    };
                    acc += (din * (u[qi] - u[p]) / h + dv[p] * dudn) / (0.5 * h);
                }
            }
            if let Some(s) = self.p.source {
                acc += s(t, x, u[p]);
            }
            r[p] = acc;
        }
        r
    }

    /// One θ-method step solved by Newton–Krylov.
    fn step(&mut self, u: &[f64], t: f64, dt: f64, theta: f64) -> Result<Vec<f64>, String> {
        let nn = u.len();
        let r0 = self.rhs(u, t);
        let mut gd = vec![0.0; nn];
        self.op.set_dirichlet(&mut gd, t + dt);
        let resid = |w: &[f64], this: &Self| -> Vec<f64> {
            let r1 = this.rhs(w, t + dt);
            (0..nn)
                .map(|p| {
                    if this.op.is_dirichlet(p) {
                        w[p] - gd[p]
                    } else {
                        w[p] - u[p] - dt * (theta * r1[p] + (1.0 - theta) * r0[p])
                    }
                })
                .collect()
        };
        let mut w = u.to_vec();
        self.op.set_dirichlet(&mut w, t + dt);
        let mut f = resid(&w, self);
        let norm = |v: &[f64]| v.iter().map(|x| x * x).sum::<f64>().sqrt();
        let f0n = norm(&f).max(1e-300);
        for _ in 0..25 {
            if norm(&f) <= 1e-10 * (1.0 + norm(&w)) || norm(&f) <= 1e-12 * f0n {
                break;
            }
            self.stats.newton_iters += 1;
            // frozen-coefficient preconditioner
            let dmean = w.iter().map(|&v| (self.p.d)(v)).sum::<f64>() / nn as f64;
            let sys = LinearSystem::new(self.op.system(1.0, -theta * dt * dmean), LinearSolver::Direct)?;
            let wn = norm(&w);
            let this: &Self = self;
            let jv = |v: &[f64], out: &mut [f64]| {
                let vn = norm(v).max(1e-300);
                let eps = 1e-7 * (1.0 + wn) / vn;
                let wp: Vec<f64> = w.iter().zip(v).map(|(a, b)| a + eps * b).collect();
                let fp = resid(&wp, this);
                for p in 0..nn {
                    out[p] = (fp[p] - f[p]) / eps;
                }
            };
            let pre = |r: &[f64], z: &mut [f64]| {
                let s = sys.solve(r, None).unwrap_or_else(|_| r.to_vec());
                z.copy_from_slice(&s);
            };
            let rhs_vec: Vec<f64> = f.iter().map(|x| -x).collect();
            let k = gmres(jv, pre, &rhs_vec, None, 40, 1e-8, 400);
            self.stats.krylov_iters += k.iterations;
            for p in 0..nn {
                w[p] += k.x[p];
            }
            f = resid(&w, self);
        }
        if norm(&f) > 1e-6 * (1.0 + norm(&w)) {
            return Err("Newton iteration did not converge".into());
        }
        Ok(w)
    }

    fn initial(&self) -> Vec<f64> {
        let mut u = self.p.grid.sample(|x| (self.p.init)(x));
        self.op.set_dirichlet(&mut u, 0.0);
        u
    }
}

/// Fixed-step implicit integration (θ-method) of nonlinear diffusion by
/// Newton–Krylov.
///
/// # Errors
/// If Newton fails to converge or `theta` is outside `[1/2, 1]`.
pub fn nonlinear_diffusion(
    prob: &NonlinearDiffusion<'_>,
    t_end: f64,
    steps: usize,
    theta: f64,
) -> Result<(Vec<f64>, NewtonStats), String> {
    if !(0.5..=1.0).contains(&theta) || steps == 0 {
        return Err("theta must be in [1/2, 1] and steps > 0".into());
    }
    let mut st = NlStepper::new(prob)?;
    let mut u = st.initial();
    let dt = t_end / steps as f64;
    for k in 0..steps {
        u = st.step(&u, k as f64 * dt, dt, theta)?;
    }
    Ok((u, st.stats))
}

/// Adaptive (step-doubling, Crank–Nicolson) Newton–Krylov integration of
/// nonlinear diffusion.
///
/// # Errors
/// If Newton fails to converge or `max_steps` is exceeded.
pub fn nonlinear_diffusion_adaptive(
    prob: &NonlinearDiffusion<'_>,
    t_end: f64,
    opts: &AdaptiveOptions,
) -> Result<(Vec<f64>, AdaptiveStats, NewtonStats), String> {
    let mut st = NlStepper::new(prob)?;
    let u0 = st.initial();
    let cell = std::cell::RefCell::new(&mut st);
    let (u, stats) = adaptive_drive(u0, t_end, 2, opts, |u, t, dt| cell.borrow_mut().step(u, t, dt, 0.5))?;
    let ns = cell.borrow().stats;
    Ok((u, stats, ns))
}
