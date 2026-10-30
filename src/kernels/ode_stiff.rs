//! # Stiff ODE integrators: Radau IIA(5/9/13) and variable-order BDF
//!
//! * [`radau5`] / [`radau5_jac`]: the three-stage, order-5, L-stable Radau IIA
//!   collocation method with simplified Newton iteration (one LU
//!   factorisation of the stage matrix reused across Newton iterations and
//!   steps), an embedded error estimator and adaptive step control;
//! * [`radau_with`] (options in [`RadauOptions`]): the general driver with 3,
//!   5 and 7 stages (orders 5, 9, 13) selected adaptively from the Newton
//!   convergence rate, dense output (the collocation polynomial),
//!   direction-aware/terminal event detection, and mass matrices
//!   `M y' = f` (dense or sparse) for differential-algebraic systems of
//!   index 1, 2 and 3 (see below); [`radau_fixed`] / [`radau5_fixed`] are
//!   constant-step variants used to verify the convergence orders;
//! * [`bdf`] / [`bdf_jac`] / [`bdf_with`]: variable-order (1 to 5) backward
//!   differentiation formulas (fixed-leading-coefficient form, NDF
//!   modification, adaptive order/step, Jacobian/LU reuse); [`bdf_with`]
//!   (options in [`BdfOptions`]) adds dense output (backward-difference
//!   interpolant), event detection, structured Jacobians and a mass matrix
//!   (`M y' = f`, index-1 differential-algebraic systems, see below).
//!
//! ## Structured Jacobians
//!
//! [`JacStructure`] declares the Jacobian as dense, banded (`kl`, `ku`) or
//! sparse (a pattern). Banded and sparse Jacobians are obtained either from a
//! user callback ([`RadauOptions::sparse_jac`], [`BdfOptions::sparse_jac`])
//! or by finite differences with one right-hand-side evaluation per colour of
//! the column-intersection graph ([`crate::kernels::jacobian`]); the stage
//! system `I (x) M - h A (x) J` is factored by a banded LU (stages
//! interleaved per variable so the band stays `s (kl + 1) + ...` wide) or by
//! the sparse LU of [`crate::kernels::krylov`] (AMD ordering), so the cost no
//! longer scales as `(s n)^3`.
//!
//! ## Differential-algebraic systems
//!
//! With a (possibly singular) mass matrix [`radau_with`] solves
//! `M y' = f(t, y)`. For index-2 and index-3 problems in Hessenberg form the
//! variables must be ordered by index (index-1 and differential variables
//! first, then the index-2 block, then the index-3 block) and the block sizes
//! given in [`RadauOptions::index`]; exactly as in RADAU5 the error and Newton
//! norms of the index-2 components are weighted with `h` and those of the
//! index-3 components with `h^2`, which removes the `O(h^-1)`/`O(h^-2)`
//! amplification of the error estimator and restores step-size control.
//! Initial values must be consistent.
//!
//! [`bdf_with`] accepts the same mass matrices ([`BdfOptions::mass`],
//! [`BdfOptions::sparse_mass`]) for index-1 systems: the corrector equation
//! `M (d + psi) = (h / alpha) f(t, y_pred + d)` is solved by the simplified
//! Newton iteration with the matrix `M - (h / alpha) J`, the error estimate
//! covers the algebraic components as well (as in `ode15s`), and the initial
//! derivative needed to start the backward-difference array is
//! `M^-1 f(t0, y0)` when `M` is invertible and `f_q / M_qq` for the
//! differential variables with a diagonal row otherwise (zero for the
//! algebraic ones); the first step is `1e-5 (t1 - t0)` unless
//! [`OdeOptions::h0`] is given. The initial values must be consistent.
//!
//! Limitations: BDF handles index 1 only (use Radau for higher-index
//! Hessenberg systems, whose error estimate needs the `h` and `h^2`
//! weights); the mass matrix must be constant; Radau's
//! algebraic-variable error estimate is the plain collocation estimate for
//! index 1 (no extra scaling needed there).
//!
//! The classical index-3 mechanical benchmark with position-level constraints,
//! Andrews' squeezer (Haug's seven-body planar mechanism), is validated in
//! `tests/kernels/andrews_squeezer.rs` against the reference solution of the
//! Test Set for IVP Solvers (`andrews.f`, computed with PSIDE at
//! `rtol = atol = 1e-14`, `t = 0.03`), using the Hessenberg blocks
//! `n2 = 7`, `n3 = 13`; the test module records the mechanism data, the
//! reference values and the formulation.

use crate::kernels::dense::Mat;
use crate::kernels::dense::{
    self,
};
use crate::kernels::jacobian::Factor;
use crate::kernels::jacobian::Jac;
use crate::kernels::jacobian::JacStructure;
use crate::kernels::jacobian::JacobianProvider;
use crate::kernels::jacobian::Mass;
use crate::kernels::jacobian::SolverKind;
use crate::kernels::jacobian::SparseJacobianFn;
use crate::kernels::jacobian::block_triplets;
use crate::kernels::jacobian::fd_jacobian_dense;
use crate::kernels::krylov::Csr;
use crate::kernels::ode_adaptive::DenseKind;
use crate::kernels::ode_adaptive::DenseStep;
use crate::kernels::ode_adaptive::Event;
use crate::kernels::ode_adaptive::OdeError;
use crate::kernels::ode_adaptive::OdeOptions;
use crate::kernels::ode_adaptive::OdeSolution;
use crate::kernels::ode_adaptive::StepWindow;
use crate::kernels::ode_adaptive::Tol;
use crate::kernels::ode_adaptive::initial_step;
use crate::kernels::ode_adaptive::process_events;

/// User-supplied Jacobian `J(t, y)` of the right-hand side.
pub type JacobianFn<'a> = dyn Fn(f64, &[f64]) -> Mat + 'a;

/// Integration interval `[t0, t1]` and initial state `y0` of an initial-value problem.
#[derive(Clone, Copy)]
struct IvpSpan<'a> {
    t0: f64,
    t1: f64,
    y0: &'a [f64],
}

/// Fixed inputs of one simplified-Newton solve for the Radau stage increments.
#[derive(Clone, Copy)]
struct NewtonCtx<'a> {
    tab: &'a RadauTab,
    mass: Option<Mass<'a>>,
    lu: &'a Factor,
    scale: &'a [f64],
    tol: f64,
    maxit: usize,
}

fn rms(
    v: &[f64],
    scale: &[f64],
) -> f64 {
    let n = scale.len().max(1);
    let s: f64 = v
        .iter()
        .enumerate()
        .map(|(i, x)| (x / scale[i % n]).powi(2))
        .sum();
    (s / v.len().max(1) as f64).sqrt()
}

// ------------------------------------------------------------ Radau IIA

/// Legendre polynomial `P_k(x)` by the three-term recurrence.
fn legendre(
    k: usize,
    x: f64,
) -> f64 {
    let (mut p0, mut p1) = (1.0, x);
    if k == 0 {
        return 1.0;
    }
    for m in 1..k {
        let mf = m as f64;
        let p2 = ((2.0 * mf + 1.0) * x * p1 - mf * p0) / (mf + 1.0);
        p0 = p1;
        p1 = p2;
    }
    p1
}

/// Tableau of the `s`-stage Radau IIA collocation method (order `2s - 1`)
/// together with the data of the embedded error estimator.
#[derive(Debug, Clone)]
struct RadauTab {
    s: usize,
    c: Vec<f64>,
    a: Vec<Vec<f64>>,
    gamma0: f64,
    /// Weights of the stage values `Z_k` in the error estimator.
    w: Vec<f64>,
}

impl RadauTab {
    fn new(s: usize) -> Option<Self> {
        // Nodes: the roots of P_{s-1} - P_s other than x = 1, mapped to [0, 1].
        let g = |x: f64| legendre(s - 1, x) - legendre(s, x);
        let grid = 8000usize;
        let mut c = Vec::with_capacity(s);
        let mut xa = -1.0;
        let mut ga = g(xa);
        for i in 1..=grid {
            let xb = -1.0 + 2.0 * (1.0 - 1e-7) * i as f64 / grid as f64;
            let gb = g(xb);
            if ga == 0.0 {
                c.push(f64::midpoint(xa, 1.0));
            } else if ga * gb < 0.0 {
                let (mut lo, mut hi, mut glo) = (xa, xb, ga);
                for _ in 0..200 {
                    let mid = f64::midpoint(lo, hi);
                    let gm = g(mid);
                    if gm * glo <= 0.0 {
                        hi = mid;
                    } else {
                        lo = mid;
                        glo = gm;
                    }
                }
                c.push(f64::midpoint(f64::midpoint(lo, hi), 1.0));
            }
            xa = xb;
            ga = gb;
        }
        c.push(1.0);
        if c.len() != s {
            return None;
        }
        // A_ij = int_0^{c_i} l_j(tau) d tau, with l_j the Lagrange basis.
        let mut a = vec![vec![0.0; s]; s];
        for j in 0..s {
            let mut poly = vec![1.0];
            for k in 0..s {
                if k == j {
                    continue;
                }
                let den = c[j] - c[k];
                let mut np = vec![0.0; poly.len() + 1];
                for (m, pm) in poly.iter().enumerate() {
                    np[m + 1] += pm / den;
                    np[m] -= pm * c[k] / den;
                }
                poly = np;
            }
            for i in 0..s {
                a[i][j] = poly
                    .iter()
                    .enumerate()
                    .map(|(m, pm)| pm * c[i].powi(m as i32 + 1) / (m as f64 + 1.0))
                    .sum();
            }
        }
        let am = Mat::from_rows(&a).ok()?;
        let ainv = dense::lu_factor(&am).ok()?.inverse();
        let eig = dense::eigenvalues(&ainv).ok()?;
        let mu = eig.iter().min_by(|x, y| x.1.abs().total_cmp(&y.1.abs()))?.0;
        if mu <= 0.0 {
            return None;
        }
        let gamma0 = 1.0 / mu;
        // Embedded weights: gamma0 f(y0) + sum bh_j f_j integrates degree < s exactly.
        let mut v = Mat::zeros(s, s);
        let mut rhs = vec![0.0; s];
        for k in 0..s {
            for j in 0..s {
                v.set(k, j, c[j].powi(k as i32));
            }
            rhs[k] = if k == 0 {
                1.0 - gamma0
            } else {
                1.0 / (k as f64 + 1.0)
            };
        }
        let bh = dense::solve(&v, &rhs).ok()?;
        let w: Vec<f64> = (0..s)
            .map(|k| (0..s).map(|j| (bh[j] - a[s - 1][j]) * ainv.at(j, k)).sum())
            .collect();
        Some(Self { s, c, a, gamma0, w })
    }
}

/// Linear-algebra set-up of the stage systems (dimension, back end, band).
#[derive(Clone, Copy)]
struct SolverSetup {
    n: usize,
    kind: SolverKind,
    band: (usize, usize),
}

/// Factors the full `s n x s n` stage matrix `I (x) M - h (A (x) J)` and the
/// `n x n` matrix `M - h gamma0 J` of the error estimator; `mass = None`
/// stands for the identity.
fn radau_factors(
    tab: &RadauTab,
    mass: Option<Mass<'_>>,
    jac: &Jac,
    h: f64,
    setup: SolverSetup,
) -> Option<(Factor, Factor)> {
    let n = setup.n;
    let s = tab.s;
    let diag = |bi: usize, bj: usize| if bi == bj { 1.0 } else { 0.0 };
    let full = block_triplets(n, s, &diag, &|bi, bj| -h * tab.a[bi][bj], mass, jac);
    let real = block_triplets(n, 1, &|_, _| 1.0, &|_, _| -h * tab.gamma0, mass, jac);
    let band = setup.band;
    let l1 = Factor::from_triplets(n, s, &full, setup.kind, band)?;
    let l2 = Factor::from_triplets(n, 1, &real, setup.kind, band)?;
    Some((l1, l2))
}

/// Simplified Newton iteration for the stage increments `Z`
/// (`M Z_i = h sum_j a_ij f(t + c_j h, y + Z_j)`), started at `z`.
/// Returns the increments, the iteration count and the observed rate.
fn radau_newton<F: Fn(f64, &[f64], &mut [f64])>(
    f: &F,
    t: f64,
    y: &[f64],
    h: f64,
    ctx: &NewtonCtx<'_>,
    mut z: Vec<f64>,
    nfev: &mut usize,
) -> Option<(Vec<f64>, usize, f64)> {
    let NewtonCtx {
        tab,
        mass,
        lu,
        scale,
        tol,
        maxit,
    } = *ctx;
    let n = y.len();
    let s = tab.s;
    let mut fv = vec![0.0; s * n];
    let mut yj = vec![0.0; n];
    let mut dn_old: Option<f64> = None;
    let mut theta: f64 = 0.0;
    for k in 0..maxit {
        for j in 0..s {
            for q in 0..n {
                yj[q] = y[q] + z[j * n + q];
            }
            f(t + tab.c[j] * h, &yj, &mut fv[j * n..(j + 1) * n]);
        }
        *nfev += s;
        if fv.iter().any(|v| !v.is_finite()) {
            return None;
        }
        let mut rhs = vec![0.0; s * n];
        for i in 0..s {
            let mz = mass.map(|mm| mm.matvec(&z[i * n..(i + 1) * n]));
            for q in 0..n {
                let acc: f64 = (0..s).map(|j| tab.a[i][j] * fv[j * n + q]).sum();
                let lhs = mz.as_ref().map_or(z[i * n + q], |v| v[q]);
                rhs[i * n + q] = -(lhs - h * acc);
            }
        }
        let dz = lu.solve(&rhs);
        let dn = rms(&dz, scale);
        let rate = dn_old.map(|o| dn / o);
        if let Some(r) = rate {
            theta = theta.max(r);
            if r >= 1.0 || r.powi((maxit - k) as i32) / (1.0 - r) * dn > tol {
                return None;
            }
        }
        for (zi, d) in z.iter_mut().zip(&dz) {
            *zi += d;
        }
        if dn == 0.0 || rate.is_some_and(|r| r / (1.0 - r) * dn < tol) {
            return Some((z, k + 1, theta));
        }
        dn_old = Some(dn);
    }
    None
}

/// Configuration of the deprecated [`radau`]: number of stages (3, 5 or 7,
/// i.e. method orders 5, 9 and 13) and an optional mass matrix. Prefer
/// [`RadauOptions`].
#[derive(Debug, Clone, Copy)]
pub struct RadauConfig<'a> {
    /// Mass matrix `M` of `M y' = f(t, y)` (`None` is the identity).
    /// `M` may be singular (index-1 differential-algebraic systems with
    /// consistent initial values).
    pub mass: Option<&'a Mat>,
    /// Lowest number of stages used (3, 5 or 7).
    pub min_stages: usize,
    /// Highest number of stages used (3, 5 or 7).
    pub max_stages: usize,
}

impl Default for RadauConfig<'_> {
    fn default() -> Self {
        Self {
            mass: None,
            min_stages: 3,
            max_stages: 7,
        }
    }
}

/// Block sizes of the higher-index variables of a Hessenberg DAE.
///
/// The DAE is `M y' = f(t, y)`. The variables are ordered by index: first the index-1
/// (differential and index-1 algebraic) variables, then `n2` index-2 variables
/// and finally `n3` index-3 variables (the last `n3` components of `y`).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct DaeIndex {
    /// Number of index-2 variables (error weighted with `h`).
    pub n2: usize,
    /// Number of index-3 variables (error weighted with `h^2`).
    pub n3: usize,
}

/// Options of [`radau_with`]: Jacobian source and structure, mass matrix,
/// DAE index blocks, stage range and events.
#[derive(Clone, Copy)]
pub struct RadauOptions<'a> {
    /// Dense Jacobian callback `df/dy` (`None`: finite differences exploiting
    /// [`RadauOptions::structure`]).
    pub jac: Option<&'a JacobianFn<'a>>,
    /// Sparse Jacobian callback; takes precedence over `jac` and implies a
    /// sparse LU when `structure` is [`JacStructure::Dense`].
    pub sparse_jac: Option<&'a SparseJacobianFn<'a>>,
    /// Declared Jacobian structure (dense, banded or sparse pattern).
    pub structure: JacStructure<'a>,
    /// Dense mass matrix `M` of `M y' = f(t, y)` (`None` is the identity).
    pub mass: Option<&'a Mat>,
    /// Sparse mass matrix (used when `mass` is `None`).
    pub sparse_mass: Option<&'a Csr>,
    /// Sizes of the index-2 and index-3 variable blocks.
    pub index: DaeIndex,
    /// Lowest number of stages used (3, 5 or 7).
    pub min_stages: usize,
    /// Highest number of stages used (3, 5 or 7).
    pub max_stages: usize,
    /// Event functions.
    pub events: &'a [Event<'a>],
}

impl Default for RadauOptions<'_> {
    fn default() -> Self {
        Self {
            jac: None,
            sparse_jac: None,
            structure: JacStructure::Dense,
            mass: None,
            sparse_mass: None,
            index: DaeIndex::default(),
            min_stages: 3,
            max_stages: 7,
            events: &[],
        }
    }
}

/// Radau IIA(5) integration of `y' = f(t, y)` with a finite-difference
/// Jacobian; see [`radau5_jac`].
///
/// # Errors
/// [`OdeError::StepLimit`], [`OdeError::NonFinite`], [`OdeError::SolveFailed`]
/// or [`OdeError::Invalid`].
pub fn radau5<F: Fn(f64, &[f64], &mut [f64])>(
    f: F,
    t0: f64,
    t1: f64,
    y0: &[f64],
    opts: &OdeOptions,
) -> Result<OdeSolution, OdeError> {
    let ro = RadauOptions {
        min_stages: 3,
        max_stages: 3,
        ..RadauOptions::default()
    };
    radau_core(&f, IvpSpan { t0, t1, y0 }, opts, &ro)
}

/// Radau IIA(5) integration with a user-supplied Jacobian `jac(t, y)`.
///
/// The method is L-stable, stiffly accurate and of order 5. The stage
/// system is solved by simplified Newton iteration with a single LU
/// factorisation reused across iterations and across steps whose step size
/// is (nearly) unchanged; the Jacobian is refreshed only when the Newton
/// iteration fails to converge. The error estimate is the embedded
/// third-order formula of Hairer and Wanner.
///
/// # Errors
/// [`OdeError::StepLimit`], [`OdeError::NonFinite`], [`OdeError::SolveFailed`]
/// or [`OdeError::Invalid`].
pub fn radau5_jac<F, J>(
    f: F,
    jac: J,
    t0: f64,
    t1: f64,
    y0: &[f64],
    opts: &OdeOptions,
) -> Result<OdeSolution, OdeError>
where
    F: Fn(f64, &[f64], &mut [f64]),
    J: Fn(f64, &[f64]) -> Mat,
{
    let ro = RadauOptions {
        jac: Some(&jac),
        min_stages: 3,
        max_stages: 3,
        ..RadauOptions::default()
    };
    radau_core(&f, IvpSpan { t0, t1, y0 }, opts, &ro)
}

/// Variable-order Radau IIA integration (orders 5, 9 and 13) of
/// `M y' = f(t, y)` with dense output (the collocation polynomial), event
/// detection, structured Jacobians and index-1, 2 and 3 DAEs.
///
/// The number of stages is raised while the simplified Newton iteration
/// converges very fast (smooth, tightly toleranced problems profit from
/// the high orders) and lowered when it converges slowly; the accepted
/// orders are reported in [`OdeSolution::orders`]. Events are handled exactly
/// as in [`crate::kernels::ode_adaptive::dopri5`]: direction-aware, located
/// on the collocation polynomial, terminal events end the run.
///
/// With a singular mass matrix the system is a DAE; the initial values must
/// be consistent, and unless `opts.h0` is given the first step is
/// `1e-5 |t1 - t0|`. For index 2 or 3 give the variable blocks in
/// [`RadauOptions::index`] (see the module documentation).
///
/// # Errors
/// [`OdeError::StepLimit`], [`OdeError::NonFinite`], [`OdeError::SolveFailed`]
/// or [`OdeError::Invalid`] (also for unsupported stage counts, mismatched
/// mass/pattern/index sizes).
pub fn radau_with<F: Fn(f64, &[f64], &mut [f64])>(
    f: F,
    t0: f64,
    t1: f64,
    y0: &[f64],
    opts: &OdeOptions,
    ro: &RadauOptions<'_>,
) -> Result<OdeSolution, OdeError> {
    radau_core(&f, IvpSpan { t0, t1, y0 }, opts, ro)
}

/// Variable-order Radau IIA integration; superseded by [`radau_with`].
///
/// # Errors
/// As [`radau_with`].
#[deprecated(note = "use `radau_with` with `RadauOptions`")]
#[allow(clippy::too_many_arguments)] // legacy signature kept as a thin wrapper
pub fn radau<F: Fn(f64, &[f64], &mut [f64])>(
    f: F,
    jac: Option<&JacobianFn<'_>>,
    t0: f64,
    t1: f64,
    y0: &[f64],
    opts: &OdeOptions,
    events: &[Event<'_>],
    cfg: &RadauConfig<'_>,
) -> Result<OdeSolution, OdeError> {
    let ro = RadauOptions {
        jac,
        mass: cfg.mass,
        min_stages: cfg.min_stages,
        max_stages: cfg.max_stages,
        events,
        ..RadauOptions::default()
    };
    radau_core(&f, IvpSpan { t0, t1, y0 }, opts, &ro)
}

fn radau_core<F: Fn(f64, &[f64], &mut [f64])>(
    f: &F,
    span: IvpSpan<'_>,
    opts: &OdeOptions,
    ro: &RadauOptions<'_>,
) -> Result<OdeSolution, OdeError> {
    let IvpSpan { t0, t1, y0 } = span;
    let n = y0.len();
    if n == 0 || crate::kernels::exact_eq(t0, t1) {
        return Err(OdeError::Invalid);
    }
    let valid = |s: usize| matches!(s, 3 | 5 | 7);
    if !valid(ro.min_stages) || !valid(ro.max_stages) || ro.min_stages > ro.max_stages {
        return Err(OdeError::Invalid);
    }
    if ro.index.n2 + ro.index.n3 > n {
        return Err(OdeError::Invalid);
    }
    let mass: Option<Mass<'_>> = ro
        .mass
        .map(Mass::Dense)
        .or_else(|| ro.sparse_mass.map(Mass::Sparse));
    if let Some(m) = mass
        && (m.rows() != n || m.cols() != n)
    {
        return Err(OdeError::Invalid);
    }
    let provider =
        JacobianProvider::new(n, ro.structure, ro.jac, ro.sparse_jac).ok_or(OdeError::Invalid)?;
    let mut band = provider.band;
    if let Some(m) = mass {
        let (ml, mu) = m.band();
        band = (band.0.max(ml), band.1.max(mu));
    }
    let setup = SolverSetup {
        n,
        kind: provider.kind,
        band,
    };
    let events = ro.events;
    let mut tabs: Vec<RadauTab> = Vec::new();
    for s in [3usize, 5, 7] {
        tabs.push(RadauTab::new(s).ok_or(OdeError::SolveFailed)?);
    }
    let tab_index = |s: usize| (s - 3) / 2;
    let maxit_of = |s: usize| {
        match s {
            | 3 => 7,
            | 5 => 10,
            | _ => 12,
        }
    };
    let dir = if t1 > t0 { 1.0 } else { -1.0 };
    let rtol = opts.rtol.max(100.0 * f64::EPSILON);
    let atol = opts.atol;
    let newton_tol = (10.0 * f64::EPSILON / rtol).max(0.03_f64.min(rtol.sqrt()));
    let mut stages = ro.min_stages;
    let mut t = t0;
    let mut y = y0.to_vec();
    let mut f0 = vec![0.0; n];
    f(t, &y, &mut f0);
    let mut nfev = 1;
    let mut h_abs = match (opts.h0, mass) {
        | (Some(h), _) => h.abs(),
        | (None, None) => initial_step(f, t0, y0, &f0, dir, 5.0, Tol { rtol, atol }),
        | (None, Some(_)) => 1e-5 * (t1 - t0).abs(),
    }
    .min(opts.hmax);
    let mut jac = provider.eval(f, t, &y, &f0, &mut nfev);
    let mut current_jac = true;
    let mut lus: Option<(Factor, Factor)> = None;
    let mut sol = OdeSolution::start(t0, y0);
    let mut gprev: Vec<f64> = events.iter().map(|e| (e.g)(t0, y0)).collect();
    let mut prev_dense: Option<crate::kernels::ode_adaptive::DenseStep> = None;
    let mut since_change = 0usize;
    let mut steps = 0;
    let mut rejected_last = false;
    let (n2, n3) = (ro.index.n2, ro.index.n3);
    // weight of component q in the error/Newton norms: 1, h or h^2 by index
    let index_weight = |q: usize, hh: f64| -> f64 {
        if q + n3 >= n {
            hh * hh
        } else if q + n3 + n2 >= n {
            hh
        } else {
            1.0
        }
    };
    while (t1 - t) * dir > 0.0 {
        if steps >= opts.max_steps {
            return Err(OdeError::StepLimit);
        }
        steps += 1;
        let min_step = 10.0 * f64::EPSILON * t.abs().max(1e-300);
        h_abs = h_abs.min(opts.hmax);
        let mut accepted = false;
        while !accepted {
            if h_abs < min_step || h_abs < 1e-300 {
                return Err(OdeError::StepLimit);
            }
            let tab = &tabs[tab_index(stages)];
            let s = tab.s;
            let mut h = h_abs * dir;
            let mut t_new = t + h;
            if (t_new - t1) * dir > 0.0 {
                t_new = t1;
                h = t_new - t;
                h_abs = h.abs();
                lus = None;
            }
            if lus.is_none() {
                if let Some(pair) = radau_factors(tab, mass, &jac, h, setup) {
                    lus = Some(pair);
                } else {
                    h_abs *= 0.5;
                    sol.rejected += 1;
                    continue;
                }
            }
            let Some((lu_full, lu_real)) = lus.as_ref() else {
                return Err(OdeError::SolveFailed);
            };
            let scale0: Vec<f64> = (0..n)
                .map(|q| (atol + rtol * y[q].abs()) / index_weight(q, h_abs))
                .collect();
            let default_guess = || -> Vec<f64> {
                let mut z = vec![0.0; s * n];
                if mass.is_none() {
                    for i in 0..s {
                        for k in 0..n {
                            z[i * n + k] = tab.c[i] * h * f0[k];
                        }
                    }
                }
                z
            };
            let extrap_guess = || -> Option<Vec<f64>> {
                let pd = prev_dense.as_ref()?;
                let mut z = vec![0.0; s * n];
                for i in 0..s {
                    let v = crate::kernels::ode_adaptive::dense_eval(pd, t + tab.c[i] * h);
                    for k in 0..n {
                        z[i * n + k] = v[k] - y[k];
                    }
                }
                z.iter().all(|v| v.is_finite()).then_some(z)
            };
            let maxit = maxit_of(s);
            let ctx = NewtonCtx {
                tab,
                mass,
                lu: lu_full,
                scale: &scale0,
                tol: newton_tol,
                maxit,
            };
            let mut newton =
                extrap_guess().and_then(|z| radau_newton(f, t, &y, h, &ctx, z, &mut nfev));
            if newton.is_none() {
                newton = radau_newton(f, t, &y, h, &ctx, default_guess(), &mut nfev);
            }
            let Some((z, iters, theta)) = newton else {
                if current_jac {
                    h_abs *= 0.5;
                } else {
                    jac = provider.eval(f, t, &y, &f0, &mut nfev);
                    current_jac = true;
                }
                lus = None;
                sol.rejected += 1;
                rejected_last = true;
                continue;
            };
            let y_new: Vec<f64> = (0..n).map(|q| y[q] + z[(s - 1) * n + q]).collect();
            if y_new.iter().any(|v| !v.is_finite()) {
                return Err(OdeError::NonFinite);
            }
            // error estimate: (M - h g0 J)^{-1} (h g0 f0 + M sum_k w_k Z_k)
            let mut wz = vec![0.0; n];
            for k in 0..s {
                for q in 0..n {
                    wz[q] += tab.w[k] * z[k * n + q];
                }
            }
            let mwz = mass.map_or_else(|| wz.clone(), |mm| mm.matvec(&wz));
            let rhs: Vec<f64> = (0..n).map(|q| h * tab.gamma0 * f0[q] + mwz[q]).collect();
            let mut err = lu_real.solve(&rhs);
            let scale: Vec<f64> = (0..n)
                .map(|q| (atol + rtol * y[q].abs().max(y_new[q].abs())) / index_weight(q, h_abs))
                .collect();
            let mut en = rms(&err, &scale);
            if en > 1.0 && mass.is_none() && (rejected_last || sol.t.len() == 1) {
                let yp: Vec<f64> = (0..n).map(|q| y[q] + err[q]).collect();
                let mut fp = vec![0.0; n];
                f(t, &yp, &mut fp);
                nfev += 1;
                let rhs2: Vec<f64> = (0..n).map(|q| h * tab.gamma0 * fp[q] + wz[q]).collect();
                err = lu_real.solve(&rhs2);
                en = rms(&err, &scale);
            }
            let safety = 0.9 * (2.0 * maxit as f64 + 1.0) / (2.0 * maxit as f64 + iters as f64);
            let expo = 1.0 / (s as f64 + 1.0);
            let factor = if en == 0.0 {
                10.0
            } else {
                (safety * en.powf(-expo)).clamp(0.2, 10.0)
            };
            if en > 1.0 {
                h_abs *= factor.min(0.9);
                lus = None;
                sol.rejected += 1;
                rejected_last = true;
                continue;
            }
            // accept: collocation polynomial through (0, y) and (c_j, y + Z_j)
            let mut nodes = vec![0.0];
            nodes.extend_from_slice(&tab.c);
            let mut vals = vec![y.clone()];
            for j in 0..s {
                vals.push((0..n).map(|q| y[q] + z[j * n + q]).collect());
            }
            let step = crate::kernels::ode_adaptive::DenseStep {
                t0: t,
                h,
                y0: y.clone(),
                kind: crate::kernels::ode_adaptive::DenseKind::Lagrange { nodes, vals },
            };
            prev_dense = Some(step.clone());
            let t_old = t;
            let window = StepWindow {
                t_old,
                tn: t_new,
                ynew: &y_new,
                dir,
            };
            if process_events(events, &mut gprev, &mut sol, step, &window) {
                sol.orders.push(2 * s - 1);
                sol.nfev = nfev;
                return Ok(sol);
            }
            accepted = true;
            t = t_new;
            y = y_new;
            sol.t.push(t);
            sol.y.push(y.clone());
            sol.orders.push(2 * s - 1);
            f(t, &y, &mut f0);
            nfev += 1;
            current_jac = false;
            rejected_last = false;
            since_change += 1;
            let h_changed = !(1.0..1.2).contains(&factor);
            if h_changed {
                h_abs *= factor;
                lus = None;
            }
            // order selection by the Newton convergence rate
            if since_change >= 5 {
                if theta < 2e-3 && stages < ro.max_stages {
                    stages += 2;
                    since_change = 0;
                    lus = None;
                } else if theta > 0.2 && stages > ro.min_stages {
                    stages -= 2;
                    since_change = 0;
                    lus = None;
                }
            }
        }
    }
    sol.nfev = nfev;
    Ok(sol)
}

/// Constant-step Radau IIA(5): `n_steps` equal steps from `t0` to `t1`.
///
/// The Jacobian is a finite-difference one refreshed every step and the
/// Newton iteration is fully converged. Returns the final state; used to
/// measure the order of convergence (5).
///
/// # Errors
/// [`OdeError::Invalid`], [`OdeError::SolveFailed`] or [`OdeError::NonFinite`].
pub fn radau5_fixed<F: Fn(f64, &[f64], &mut [f64])>(
    f: F,
    t0: f64,
    t1: f64,
    y0: &[f64],
    n_steps: usize,
) -> Result<Vec<f64>, OdeError> {
    radau_fixed(f, t0, t1, y0, n_steps, 3)
}

/// Constant-step Radau IIA with `stages` (3, 5 or 7) stages, i.e. order
/// `2 stages - 1` (5, 9, 13).
///
/// # Errors
/// [`OdeError::Invalid`], [`OdeError::SolveFailed`] or [`OdeError::NonFinite`].
pub fn radau_fixed<F: Fn(f64, &[f64], &mut [f64])>(
    f: F,
    t0: f64,
    t1: f64,
    y0: &[f64],
    n_steps: usize,
    stages: usize,
) -> Result<Vec<f64>, OdeError> {
    let n = y0.len();
    if n == 0 || n_steps == 0 || crate::kernels::exact_eq(t0, t1) || !matches!(stages, 3 | 5 | 7) {
        return Err(OdeError::Invalid);
    }
    let tab = RadauTab::new(stages).ok_or(OdeError::SolveFailed)?;
    let h = (t1 - t0) / n_steps as f64;
    let mut y = y0.to_vec();
    let mut f0 = vec![0.0; n];
    let mut nfev = 0;
    let scale = vec![1.0; n];
    for s in 0..n_steps {
        let t = t0 + s as f64 * h;
        f(t, &y, &mut f0);
        let jac = Jac::Dense(fd_jacobian_dense(&f, t, &y, &f0));
        let setup = SolverSetup {
            n,
            kind: SolverKind::Dense,
            band: (0, 0),
        };
        let (lu, _) = radau_factors(&tab, None, &jac, h, setup).ok_or(OdeError::SolveFailed)?;
        let mut z0 = vec![0.0; stages * n];
        for i in 0..stages {
            for k in 0..n {
                z0[i * n + k] = tab.c[i] * h * f0[k];
            }
        }
        let (z, _, _) = radau_newton(
            &f,
            t,
            &y,
            h,
            &NewtonCtx {
                tab: &tab,
                mass: None,
                lu: &lu,
                scale: &scale,
                tol: 1e-13,
                maxit: 60,
            },
            z0,
            &mut nfev,
        )
        .ok_or(OdeError::SolveFailed)?;
        for q in 0..n {
            y[q] += z[(stages - 1) * n + q];
        }
        if y.iter().any(|v| !v.is_finite()) {
            return Err(OdeError::NonFinite);
        }
    }
    Ok(y)
}

// ------------------------------------------------------------------ BDF

const BDF_MAX_ORDER: usize = 5;
const MAXIT: usize = 4;
const BDF_KAPPA: [f64; 6] = [0.0, -0.1850, -1.0 / 9.0, -0.0823, -0.0415, 0.0];

fn bdf_consts() -> ([f64; 6], [f64; 6], [f64; 7]) {
    let mut gamma = [0.0; 6];
    for k in 1..=BDF_MAX_ORDER {
        gamma[k] = gamma[k - 1] + 1.0 / k as f64;
    }
    let mut alpha = [0.0; 6];
    let mut errc = [0.0; 7];
    for k in 0..6 {
        alpha[k] = (1.0 - BDF_KAPPA[k]) * gamma[k];
    }
    for k in 0..7 {
        let kap = if k < 6 { BDF_KAPPA[k] } else { 0.0 };
        let gam = if k < 6 { gamma[k] } else { 0.0 };
        errc[k] = kap * gam + 1.0 / (k as f64 + 1.0);
    }
    (gamma, alpha, errc)
}

fn bdf_r(
    order: usize,
    factor: f64,
) -> Vec<Vec<f64>> {
    let m = order + 1;
    let mut r = vec![vec![0.0; m]; m];
    for j in 0..m {
        r[0][j] = 1.0;
    }
    for i in 1..m {
        for j in 1..m {
            r[i][j] = (i as f64 - 1.0 - factor * j as f64) / i as f64;
        }
    }
    for i in 1..m {
        for j in 0..m {
            r[i][j] *= r[i - 1][j];
        }
    }
    r
}

fn bdf_change_d(
    d: &mut [Vec<f64>],
    order: usize,
    factor: f64,
) {
    let m = order + 1;
    let r = bdf_r(order, factor);
    let u = bdf_r(order, 1.0);
    let mut ru = vec![vec![0.0; m]; m];
    for i in 0..m {
        for j in 0..m {
            ru[i][j] = (0..m).map(|k| r[i][k] * u[k][j]).sum();
        }
    }
    let n = d[0].len();
    let old: Vec<Vec<f64>> = d[..m].to_vec();
    for i in 0..m {
        for q in 0..n {
            d[i][q] = (0..m).map(|k| ru[k][i] * old[k][q]).sum();
        }
    }
}

/// Variable-order BDF (orders 1 to 5, NDF coefficients) with a
/// finite-difference Jacobian; see [`bdf_jac`].
///
/// # Errors
/// [`OdeError::StepLimit`], [`OdeError::NonFinite`], [`OdeError::SolveFailed`]
/// or [`OdeError::Invalid`].
pub fn bdf<F: Fn(f64, &[f64], &mut [f64])>(
    f: F,
    t0: f64,
    t1: f64,
    y0: &[f64],
    opts: &OdeOptions,
    max_order: usize,
) -> Result<OdeSolution, OdeError> {
    let bo = BdfOptions {
        max_order,
        ..BdfOptions::default()
    };
    bdf_core(&f, IvpSpan { t0, t1, y0 }, opts, &bo)
}

/// Variable-order BDF with a user-supplied Jacobian `jac(t, y)`.
///
/// The solution is advanced in the fixed-leading-coefficient form on a
/// backward-difference array; the order (capped at `max_order`, `1..=5`)
/// and the step size are chosen after every `order + 1` equal steps by
/// comparing the error estimates of the neighbouring orders. The
/// Jacobian and its LU factorisation are reused while the Newton
/// iteration converges.
///
/// # Errors
/// [`OdeError::StepLimit`], [`OdeError::NonFinite`], [`OdeError::SolveFailed`]
/// or [`OdeError::Invalid`].
pub fn bdf_jac<F, J>(
    f: F,
    jac: J,
    t0: f64,
    t1: f64,
    y0: &[f64],
    opts: &OdeOptions,
    max_order: usize,
) -> Result<OdeSolution, OdeError>
where
    F: Fn(f64, &[f64], &mut [f64]),
    J: Fn(f64, &[f64]) -> Mat,
{
    let bo = BdfOptions {
        jac: Some(&jac),
        max_order,
        ..BdfOptions::default()
    };
    bdf_core(&f, IvpSpan { t0, t1, y0 }, opts, &bo)
}

/// Options of [`bdf_with`]: Jacobian source and structure, maximal order and
/// events.
#[derive(Clone, Copy)]
pub struct BdfOptions<'a> {
    /// Dense Jacobian callback (`None`: finite differences exploiting
    /// [`BdfOptions::structure`]).
    pub jac: Option<&'a JacobianFn<'a>>,
    /// Sparse Jacobian callback; takes precedence over `jac`.
    pub sparse_jac: Option<&'a SparseJacobianFn<'a>>,
    /// Declared Jacobian structure (dense, banded or sparse pattern).
    pub structure: JacStructure<'a>,
    /// Highest order used (`1..=5`).
    pub max_order: usize,
    /// Event functions.
    pub events: &'a [Event<'a>],
    /// Dense mass matrix `M` of `M y' = f(t, y)` (`None` is the identity).
    /// `M` may be singular: index-1 differential-algebraic systems with
    /// consistent initial values.
    pub mass: Option<&'a Mat>,
    /// Sparse mass matrix (used when `mass` is `None`).
    pub sparse_mass: Option<&'a Csr>,
}

impl Default for BdfOptions<'_> {
    fn default() -> Self {
        Self {
            jac: None,
            sparse_jac: None,
            structure: JacStructure::Dense,
            max_order: BDF_MAX_ORDER,
            events: &[],
            mass: None,
            sparse_mass: None,
        }
    }
}

/// The derivative `y'(t0)` that starts the backward-difference array of a
/// system `M y' = f0`: `M^-1 f0` for an invertible dense `M`; otherwise the
/// differential variables whose row of `M` is a single diagonal entry get
/// `f_q / M_qq` and the algebraic ones zero.
fn initial_derivative(
    mass: Mass<'_>,
    f0: &[f64],
) -> Vec<f64> {
    let n = f0.len();
    if let Mass::Dense(m) = mass
        && let Ok(lu) = dense::lu_factor(m)
    {
        return lu.solve(f0);
    }
    let mut diagonal = vec![0.0; n];
    let mut off_diagonal = vec![false; n];
    mass.for_each(|p, q, v| {
        if v != 0.0 {
            if p == q {
                diagonal[p] = v;
            } else {
                off_diagonal[p] = true;
            }
        }
    });
    (0..n)
        .map(|q| {
            if diagonal[q] != 0.0 && !off_diagonal[q] {
                f0[q] / diagonal[q]
            } else {
                0.0
            }
        })
        .collect()
}

/// Variable-order BDF with dense output, event detection and structured
/// Jacobians.
///
/// Identical to [`bdf_jac`] (the Jacobian is optional) but additionally
/// locates `bo.events` on the backward-difference interpolating polynomial of
/// the last step (direction-aware; a terminal event ends the run), uses a
/// banded or sparse LU when [`BdfOptions::structure`] says so, and
/// [`OdeSolution::eval`] interpolates anywhere in the integrated range.
///
/// # Errors
/// [`OdeError::StepLimit`], [`OdeError::NonFinite`], [`OdeError::SolveFailed`]
/// or [`OdeError::Invalid`].
pub fn bdf_with<F: Fn(f64, &[f64], &mut [f64])>(
    f: F,
    t0: f64,
    t1: f64,
    y0: &[f64],
    opts: &OdeOptions,
    bo: &BdfOptions<'_>,
) -> Result<OdeSolution, OdeError> {
    bdf_core(&f, IvpSpan { t0, t1, y0 }, opts, bo)
}

/// Variable-order BDF with dense output and event detection; superseded by
/// [`bdf_with`].
///
/// # Errors
/// As [`bdf_with`].
#[deprecated(note = "use `bdf_with` with `BdfOptions`")]
#[allow(clippy::too_many_arguments)] // legacy signature kept as a thin wrapper
pub fn bdf_ex<F: Fn(f64, &[f64], &mut [f64])>(
    f: F,
    jac: Option<&JacobianFn<'_>>,
    t0: f64,
    t1: f64,
    y0: &[f64],
    opts: &OdeOptions,
    max_order: usize,
    events: &[Event<'_>],
) -> Result<OdeSolution, OdeError> {
    let bo = BdfOptions {
        jac,
        max_order,
        events,
        ..BdfOptions::default()
    };
    bdf_core(&f, IvpSpan { t0, t1, y0 }, opts, &bo)
}

fn bdf_core<F: Fn(f64, &[f64], &mut [f64])>(
    f: &F,
    span: IvpSpan<'_>,
    opts: &OdeOptions,
    bo: &BdfOptions<'_>,
) -> Result<OdeSolution, OdeError> {
    let IvpSpan { t0, t1, y0 } = span;
    let max_order = bo.max_order;
    let events = bo.events;
    let n = y0.len();
    if n == 0 || crate::kernels::exact_eq(t0, t1) || max_order == 0 {
        return Err(OdeError::Invalid);
    }
    let max_order = max_order.min(BDF_MAX_ORDER);
    let mass: Option<Mass<'_>> = bo
        .mass
        .map(Mass::Dense)
        .or_else(|| bo.sparse_mass.map(Mass::Sparse));
    if let Some(m) = mass
        && (m.rows() != n || m.cols() != n)
    {
        return Err(OdeError::Invalid);
    }
    let provider =
        JacobianProvider::new(n, bo.structure, bo.jac, bo.sparse_jac).ok_or(OdeError::Invalid)?;
    let mut band = provider.band;
    if let Some(m) = mass {
        let (ml, mu) = m.band();
        band = (band.0.max(ml), band.1.max(mu));
    }
    let dir = if t1 > t0 { 1.0 } else { -1.0 };
    let (gamma, alpha, errc) = bdf_consts();
    let rtol = opts.rtol.max(100.0 * f64::EPSILON);
    let atol = opts.atol;
    let newton_tol = (10.0 * f64::EPSILON / rtol).max(0.03_f64.min(rtol.sqrt()));
    let mut t = t0;
    let mut y = y0.to_vec();
    let mut f0 = vec![0.0; n];
    f(t, &y, &mut f0);
    let mut nfev = 1;
    let mut h_abs = match (opts.h0, mass) {
        | (Some(h), _) => h,
        | (None, None) => initial_step(f, t0, y0, &f0, dir, 1.0, Tol { rtol, atol }),
        | (None, Some(_)) => 1e-5 * (t1 - t0).abs(),
    }
    .abs()
    .min(opts.hmax);
    let ydot0 = mass.map_or_else(|| f0.clone(), |m| initial_derivative(m, &f0));
    let mut d = vec![vec![0.0; n]; BDF_MAX_ORDER + 3];
    d[0].clone_from(&y);
    for q in 0..n {
        d[1][q] = ydot0[q] * h_abs * dir;
    }
    let mut jac = provider.eval(f, t, &y, &f0, &mut nfev);
    let mut order = 1usize;
    let mut n_equal = 0usize;
    let mut lu: Option<Factor> = None;
    let mut sol = OdeSolution::start(t0, y0);
    let mut gprev: Vec<f64> = events.iter().map(|e| (e.g)(t0, y0)).collect();
    let mut steps = 0;
    while (t1 - t) * dir > 0.0 {
        if steps >= opts.max_steps {
            return Err(OdeError::StepLimit);
        }
        steps += 1;
        let min_step = 10.0 * f64::EPSILON * t.abs().max(1e-300);
        if h_abs > opts.hmax {
            bdf_change_d(&mut d, order, opts.hmax / h_abs);
            n_equal = 0;
            h_abs = opts.hmax;
        } else if h_abs < min_step {
            bdf_change_d(&mut d, order, min_step / h_abs);
            n_equal = 0;
            h_abs = min_step;
        }
        let mut current_jac = false;
        let mut accepted = false;
        let mut d_inc = vec![0.0; n];
        let mut y_new = vec![0.0; n];
        let mut t_new = t;
        let mut err_norm_v = 0.0;
        let mut safety = 0.9;
        let mut scale = vec![1.0; n];
        while !accepted {
            if h_abs < min_step {
                return Err(OdeError::StepLimit);
            }
            let mut h = h_abs * dir;
            t_new = t + h;
            if (t_new - t1) * dir > 0.0 {
                t_new = t1;
                bdf_change_d(&mut d, order, (t_new - t).abs() / h_abs);
                n_equal = 0;
                lu = None;
            }
            h = t_new - t;
            h_abs = h.abs();
            let y_pred: Vec<f64> = (0..n).map(|q| (0..=order).map(|k| d[k][q]).sum()).collect();
            let sc: Vec<f64> = y_pred.iter().map(|v| atol + rtol * v.abs()).collect();
            let psi: Vec<f64> = (0..n)
                .map(|q| (1..=order).map(|k| d[k][q] * gamma[k]).sum::<f64>() / alpha[order])
                .collect();
            let cc = h / alpha[order];
            let mut converged = false;
            let mut n_iter = 0;
            while !converged {
                if lu.is_none() {
                    let trip = block_triplets(n, 1, &|_, _| 1.0, &|_, _| -cc, mass, &jac);
                    match Factor::from_triplets(n, 1, &trip, provider.kind, band) {
                        | Some(l) => lu = Some(l),
                        | None => break,
                    }
                }
                let Some(l) = lu.as_ref() else { break };
                // Newton iteration
                let mut yy = y_pred.clone();
                let mut dd = vec![0.0; n];
                let mut dn_old: Option<f64> = None;
                let mut fv = vec![0.0; n];
                converged = false;
                for k in 0..MAXIT {
                    f(t_new, &yy, &mut fv);
                    nfev += 1;
                    n_iter = k + 1;
                    if fv.iter().any(|v| !v.is_finite()) {
                        break;
                    }
                    // `c f - M (psi + d)`.
                    let shifted: Vec<f64> = (0..n).map(|q| psi[q] + dd[q]).collect();
                    let m_shifted = mass.map_or_else(|| shifted.clone(), |m| m.matvec(&shifted));
                    let rhs: Vec<f64> = (0..n).map(|q| cc * fv[q] - m_shifted[q]).collect();
                    let dy = l.solve(&rhs);
                    let dn = rms(&dy, &sc);
                    let rate = dn_old.map(|o| dn / o);
                    if let Some(r) = rate
                        && (r >= 1.0 || r.powi((MAXIT - k) as i32) / (1.0 - r) * dn > newton_tol)
                    {
                        break;
                    }
                    for q in 0..n {
                        yy[q] += dy[q];
                        dd[q] += dy[q];
                    }
                    if dn == 0.0 || rate.is_some_and(|r| r / (1.0 - r) * dn < newton_tol) {
                        converged = true;
                        break;
                    }
                    dn_old = Some(dn);
                }
                if converged {
                    y_new = yy;
                    d_inc = dd;
                } else {
                    if current_jac || lu.is_none() {
                        break;
                    }
                    let mut fp = vec![0.0; n];
                    f(t_new, &y_pred, &mut fp);
                    nfev += 1;
                    jac = provider.eval(f, t_new, &y_pred, &fp, &mut nfev);
                    lu = None;
                    current_jac = true;
                }
            }
            if !converged {
                h_abs *= 0.5;
                bdf_change_d(&mut d, order, 0.5);
                n_equal = 0;
                lu = None;
                sol.rejected += 1;
                continue;
            }
            safety = 0.9 * (2.0 * MAXIT as f64 + 1.0) / (2.0 * MAXIT as f64 + n_iter as f64);
            scale = y_new.iter().map(|v| atol + rtol * v.abs()).collect();
            let err: Vec<f64> = d_inc.iter().map(|v| errc[order] * v).collect();
            err_norm_v = rms(&err, &scale);
            if err_norm_v > 1.0 {
                let factor = (safety * err_norm_v.powf(-1.0 / (order as f64 + 1.0))).max(0.2);
                h_abs *= factor;
                bdf_change_d(&mut d, order, factor);
                n_equal = 0;
                sol.rejected += 1;
            } else {
                accepted = true;
            }
        }
        n_equal += 1;
        let t_old = t;
        let y_old = y.clone();
        let h_step = t_new - t_old;
        t = t_new;
        y = y_new;
        // jac was refreshed during this step only if Newton failed; keep it
        for q in 0..n {
            d[order + 2][q] = d_inc[q] - d[order + 1][q];
            d[order + 1][q] = d_inc[q];
        }
        for i in (0..=order).rev() {
            for q in 0..n {
                let v = d[i + 1][q];
                d[i][q] += v;
            }
        }
        let step = DenseStep {
            t0: t_old,
            h: h_step,
            y0: y_old,
            kind: DenseKind::Backward(d[..=order].to_vec()),
        };
        if process_events(
            events,
            &mut gprev,
            &mut sol,
            step,
            &StepWindow {
                t_old,
                tn: t,
                ynew: &y,
                dir,
            },
        ) {
            sol.nfev = nfev;
            return Ok(sol);
        }
        sol.t.push(t);
        sol.y.push(y.clone());
        sol.orders.push(order);
        if n_equal < order + 1 {
            continue;
        }
        let err_m = if order > 1 {
            let e: Vec<f64> = d[order].iter().map(|v| errc[order - 1] * v).collect();
            rms(&e, &scale)
        } else {
            f64::INFINITY
        };
        let err_p = if order < max_order {
            let e: Vec<f64> = d[order + 2].iter().map(|v| errc[order + 1] * v).collect();
            rms(&e, &scale)
        } else {
            f64::INFINITY
        };
        let norms = [err_m, err_norm_v, err_p];
        let mut factors = [0.0; 3];
        for (k, nv) in norms.iter().enumerate() {
            factors[k] = nv.powf(-1.0 / (order as f64 + k as f64));
        }
        let mut best = 0;
        for k in 1..3 {
            if factors[k] > factors[best] {
                best = k;
            }
        }
        order = (order + best).saturating_sub(1).clamp(1, max_order);
        let factor = (safety * factors[best]).min(10.0);
        h_abs *= factor;
        bdf_change_d(&mut d, order, factor);
        n_equal = 0;
        lu = None;
    }
    sol.nfev = nfev;
    Ok(sol)
}
