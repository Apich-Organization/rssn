//! Named, JSON-parameterised scenarios for the dimension-generic continuum
//! and particle engines.
//!
//! Each scenario builds a standard validation problem of one engine, runs
//! it, and returns the error against the exact solution (where one exists)
//! together with a one-dimensional slice of the field. Every parameter has a
//! default, so `{}` is a valid request.
//!
//! | Scenario | Engine | Problem |
//! |---|---|---|
//! | `poisson_multigrid` | [`multigrid`](crate::sim::multigrid) | `-Δu + σu = f`, V/W/FMG, 1D–3D |
//! | `heat_nd` | [`structured`](crate::sim::structured) | decaying mode, Dirichlet/Neumann/periodic, fixed or adaptive steps |
//! | `wave_nd` | [`structured`](crate::sim::structured) | standing wave, leapfrog |
//! | `euler_fv` | [`euler_nd`](crate::sim::euler_nd) | planar Sod tube in 1D–3D, MUSCL/WENO5, HLLC/Roe |
//! | `dg_conservation` | [`dg`](crate::sim::dg) | advection, Burgers, Euler (Sod) |
//! | `acoustics_pml` | [`acoustics`](crate::sim::acoustics) | pulse leaving through a PML |
//! | `yee_fdtd` | [`yee_fdtd`](crate::sim::yee_fdtd) | point source with PML |
//! | `navier_stokes_mac` | [`ns_mac`](crate::sim::ns_mac) | cavity, Taylor–Green, Poiseuille |
//! | `split_step` | [`split_step`](crate::sim::split_step) | harmonic-trap Schrödinger/GPE |
//! | `fem_simplex` | [`fem_simplex`](crate::sim::fem_simplex) | P1/P2 Poisson on intervals, triangles, tetrahedra |
//! | `elasticity_fem` | [`elasticity_fem`](crate::sim::elasticity_fem) | cantilever beam 2D/3D |
//! | `bem_laplace` | [`bem_aca`](crate::sim::bem_aca) | sphere, interior/exterior, ACA or P1 |
//! | `sph` | [`sph`](crate::sim::sph) | 1D Sod tube, optional grad-h |
//! | `dpd` | [`dpd`](crate::sim::dpd) | Groot–Warren fluid |
//! | `md_lj` | [`md_engine`](crate::sim::md_engine) | Lennard-Jones fluid, NVE or Langevin |
//! | `spin_model` | [`spin_models`](crate::sim::spin_models) | Ising/XY/Heisenberg lattice |
//! | `nbody_gravity` | [`nbody`](crate::sim::nbody) | Barnes–Hut accelerations |
//! | `kerr_geodesic` | [`relativity`](crate::sim::relativity) | circular orbit in Kerr–Newman |

use serde::Deserialize;
use serde::de::DeserializeOwned;
use serde_json::Value;
use serde_json::json;
use std::f64::consts::PI;

use crate::sim::boundary::{BcValue, BoundaryCondition as Bc, Face, FaceBcs, FluidBoundary, FluidFaceBcs, MarkerBcs};

/// The scenario names handled here.
pub const ENGINE_SCENARIOS: &[&str] = &[
    "poisson_multigrid",
    "heat_nd",
    "wave_nd",
    "euler_fv",
    "dg_conservation",
    "acoustics_pml",
    "yee_fdtd",
    "navier_stokes_mac",
    "split_step",
    "fem_simplex",
    "elasticity_fem",
    "bem_laplace",
    "sph",
    "dpd",
    "md_lj",
    "spin_model",
    "nbody_gravity",
    "kerr_geodesic",
];

fn parse<T: DeserializeOwned>(params: &str) -> Result<T, String> {
    let text = if params.trim().is_empty() { "{}" } else { params };
    serde_json::from_str(text).map_err(|e| format!("invalid parameters: {e}"))
}

fn check_dim(dim: usize, lo: usize, hi: usize) -> Result<(), String> {
    if (lo..=hi).contains(&dim) { Ok(()) } else { Err(format!("dim must be between {lo} and {hi}")) }
}

fn bc_set(name: &str) -> Result<FaceBcs, String> {
    Ok(match name {
        "dirichlet" => FaceBcs::dirichlet_zero(),
        "neumann" => FaceBcs::uniform(&Bc::neumann(0.0)),
        "periodic" => FaceBcs::uniform(&Bc::Periodic),
        other => return Err(format!("unknown bc `{other}` (dirichlet, neumann, periodic)")),
    })
}

/// Runs an engine scenario; `None` when `name` is not one of [`ENGINE_SCENARIOS`].
pub(super) fn run(name: &str, params: &str) -> Option<Result<Value, String>> {
    Some(match name {
        "poisson_multigrid" => poisson_multigrid(params),
        "heat_nd" => heat_nd(params),
        "wave_nd" => wave_nd(params),
        "euler_fv" => euler_fv(params),
        "dg_conservation" => dg_conservation(params),
        "acoustics_pml" => acoustics_pml(params),
        "yee_fdtd" => yee_fdtd(params),
        "navier_stokes_mac" => navier_stokes_mac(params),
        "split_step" => split_step(params),
        "fem_simplex" => fem_simplex(params),
        "elasticity_fem" => elasticity_fem(params),
        "bem_laplace" => bem_laplace(params),
        "sph" => sph(params),
        "dpd" => dpd(params),
        "md_lj" => md_lj(params),
        "spin_model" => spin_model(params),
        "nbody_gravity" => nbody_gravity(params),
        "kerr_geodesic" => kerr_geodesic(params),
        _ => return None,
    })
}

// ------------------------------------------------------------------ continuum

#[derive(Deserialize)]
#[serde(default)]
struct PoissonReq {
    dim: usize,
    n: usize,
    sigma: f64,
    cycle: String,
    tol: f64,
}

impl Default for PoissonReq {
    fn default() -> Self {
        Self { dim: 2, n: 33, sigma: 0.0, cycle: "v".into(), tol: 1e-10 }
    }
}

fn poisson_multigrid(params: &str) -> Result<Value, String> {
    use crate::sim::multigrid::{MgCycle, MgOptions, MultigridSolver};
    use crate::sim::structured::Grid;
    let r: PoissonReq = parse(params)?;
    check_dim(r.dim, 1, 3)?;
    let cycle = match r.cycle.as_str() {
        "v" => MgCycle::V,
        "w" => MgCycle::W,
        "fmg" => MgCycle::Fmg,
        other => return Err(format!("unknown cycle `{other}` (v, w, fmg)")),
    };
    let bcs = FaceBcs::dirichlet_zero();
    let g = Grid::unit(r.dim, r.n, &bcs)?;
    let d = r.dim;
    let exact = move |x: [f64; 3]| (0..d).map(|a| (PI * x[a]).sin()).product::<f64>();
    let f = g.sample(|x| (d as f64 * PI * PI + r.sigma) * exact(x));
    let mg = MultigridSolver::new(&g, &bcs, r.sigma, MgOptions { cycle, tol: r.tol, ..MgOptions::default() })?;
    let res = mg.solve(&f, 0.0)?;
    Ok(json!({
        "levels": mg.num_levels(),
        "residuals": res.residuals,
        "converged": res.converged,
        "convergence_factor": res.convergence_factor(),
        "max_error": g.max_error(&res.u, exact),
    }))
}

#[derive(Deserialize)]
#[serde(default)]
struct HeatReq {
    dim: usize,
    n: usize,
    alpha: f64,
    t_end: f64,
    steps: usize,
    theta: f64,
    bc: String,
    adaptive_tol: Option<f64>,
}

impl Default for HeatReq {
    fn default() -> Self {
        Self { dim: 2, n: 33, alpha: 1.0, t_end: 0.05, steps: 50, theta: 0.5, bc: "dirichlet".into(), adaptive_tol: None }
    }
}

/// Boxed scalar field on `[0, 1]^3`.
type ScalarField = Box<dyn Fn([f64; 3]) -> f64>;

fn mode_fn(bc: &str, dim: usize) -> (ScalarField, f64) {
    match bc {
        "neumann" => (Box::new(move |x| (0..dim).map(|a| (PI * x[a]).cos()).product()), dim as f64 * PI * PI),
        "periodic" => (Box::new(move |x| (0..dim).map(|a| (2.0 * PI * x[a]).sin()).product()), 4.0 * dim as f64 * PI * PI),
        _ => (Box::new(move |x| (0..dim).map(|a| (PI * x[a]).sin()).product()), dim as f64 * PI * PI),
    }
}

fn heat_nd(params: &str) -> Result<Value, String> {
    use crate::sim::structured::{AdaptiveOptions, Grid, HeatProblem, LinearSolver, heat_adaptive, heat_fixed};
    let r: HeatReq = parse(params)?;
    check_dim(r.dim, 1, 3)?;
    let bcs = bc_set(&r.bc)?;
    let g = Grid::unit(r.dim, r.n, &bcs)?;
    let (mode, rate) = mode_fn(&r.bc, r.dim);
    let prob = HeatProblem { grid: &g, bcs: &bcs, alpha: r.alpha, init: &*mode, source: None, solver: LinearSolver::Auto };
    let (u, extra) = if let Some(tol) = r.adaptive_tol {
        let (u, st) = heat_adaptive(&prob, r.t_end, &AdaptiveOptions { tol, dt0: r.t_end / 100.0, ..AdaptiveOptions::default() })?;
        (u, json!({ "accepted": st.accepted, "rejected": st.rejected }))
    } else {
        (heat_fixed(&prob, r.t_end, r.steps, r.theta)?, json!({ "steps": r.steps }))
    };
    let decay = (-r.alpha * rate * r.t_end).exp();
    let err = g.max_error(&u, |x| decay * mode(x));
    let mid = g.n[0] / 2;
    let line: Vec<f64> = (0..g.n[0]).map(|i| u[g.index([i, if r.dim > 1 { mid } else { 0 }, if r.dim > 2 { mid } else { 0 }])]).collect();
    Ok(json!({ "max_error": err, "line": line, "stats": extra }))
}

#[derive(Deserialize)]
#[serde(default)]
struct WaveReq {
    dim: usize,
    n: usize,
    c: f64,
    t_end: f64,
}

impl Default for WaveReq {
    fn default() -> Self {
        Self { dim: 2, n: 41, c: 1.0, t_end: 0.4 }
    }
}

fn wave_nd(params: &str) -> Result<Value, String> {
    use crate::sim::structured::{Grid, WaveProblem, wave_leapfrog};
    let r: WaveReq = parse(params)?;
    check_dim(r.dim, 1, 3)?;
    let bcs = FaceBcs::dirichlet_zero();
    let g = Grid::unit(r.dim, r.n, &bcs)?;
    let d = r.dim;
    let mode = move |x: [f64; 3]| (0..d).map(|a| (PI * x[a]).sin()).product::<f64>();
    let zero = |_: [f64; 3]| 0.0;
    let dt_max = 0.9 * g.h_min() / (r.c * (d as f64).sqrt());
    let steps = (r.t_end / dt_max).ceil() as usize;
    let prob = WaveProblem { grid: &g, bcs: &bcs, c: r.c, u0: &mode, v0: &zero, source: None };
    let u = wave_leapfrog(&prob, r.t_end, steps)?;
    let w = r.c * PI * (d as f64).sqrt();
    Ok(json!({ "steps": steps, "max_error": g.max_error(&u, |x| mode(x) * (w * r.t_end).cos()) }))
}

#[derive(Deserialize)]
#[serde(default)]
struct EulerReq {
    dim: usize,
    n: usize,
    t_end: f64,
    recon: String,
    riemann: String,
    gamma: f64,
}

impl Default for EulerReq {
    fn default() -> Self {
        Self { dim: 1, n: 100, t_end: 0.2, recon: "weno5".into(), riemann: "hllc".into(), gamma: 1.4 }
    }
}

fn euler_fv(params: &str) -> Result<Value, String> {
    use crate::sim::euler_nd::{FvConfig, FvGrid, Reconstruction, RiemannSolver, SlopeLimiter, solve};
    use crate::sim::gas_dynamics::{Primitive, exact_riemann};
    let r: EulerReq = parse(params)?;
    check_dim(r.dim, 1, 3)?;
    let mut ns = vec![4usize; r.dim];
    ns[0] = r.n;
    let grid = FvGrid::new(&ns, &vec![0.0; r.dim], &vec![1.0; r.dim])?;
    let mut cfg = FvConfig::euler(r.gamma);
    cfg.recon = match r.recon.as_str() {
        "weno5" => Reconstruction::Weno5,
        "first_order" => Reconstruction::FirstOrder,
        "muscl_minmod" => Reconstruction::Muscl(SlopeLimiter::Minmod),
        "muscl_vanleer" => Reconstruction::Muscl(SlopeLimiter::VanLeer),
        "muscl_mc" => Reconstruction::Muscl(SlopeLimiter::Mc),
        "muscl_superbee" => Reconstruction::Muscl(SlopeLimiter::Superbee),
        other => return Err(format!("unknown recon `{other}`")),
    };
    cfg.riemann = match r.riemann.as_str() {
        "hllc" => RiemannSolver::Hllc,
        "roe" => RiemannSolver::Roe,
        "rusanov" => RiemannSolver::Rusanov,
        other => return Err(format!("unknown riemann solver `{other}`")),
    };
    let init = |x: [f64; 3]| if x[0] < 0.5 { [1.0, 0.0, 0.0, 0.0, 1.0] } else { [0.125, 0.0, 0.0, 0.0, 0.1] };
    let sol = solve(&grid, &cfg, &init, r.t_end)?;
    let (l, rr) = (Primitive { rho: 1.0, u: 0.0, p: 1.0 }, Primitive { rho: 0.125, u: 0.0, p: 0.1 });
    let exact = |x: f64| exact_riemann(l, rr, r.gamma, (x - 0.5) / r.t_end).rho;
    let x: Vec<f64> = (0..r.n).map(|i| grid.centre(grid.index([i, 0, 0]))[0]).collect();
    let rho: Vec<f64> = (0..r.n).map(|i| sol.primitive(grid.index([i, 0, 0]))[0]).collect();
    let l1 = sol.l1_error(0, |p| exact(p[0]));
    Ok(json!({ "x": x, "density": rho, "exact": x.iter().map(|&v| exact(v)).collect::<Vec<_>>(), "l1_error": l1, "steps": sol.steps }))
}

#[derive(Deserialize)]
#[serde(default)]
struct DgReq {
    law: String,
    ne: usize,
    p: usize,
    t_end: f64,
    limiter: bool,
}

impl Default for DgReq {
    fn default() -> Self {
        Self { law: "advection".into(), ne: 20, p: 3, t_end: 0.5, limiter: false }
    }
}

fn dg_conservation(params: &str) -> Result<Value, String> {
    use crate::sim::dg::{Burgers, DgBc, DgLaw, DgLimiter, DgMesh, Euler1d, LinearAdvection, solve_dg_limited};
    use crate::sim::gas_dynamics::{Primitive, exact_riemann};
    let r: DgReq = parse(params)?;
    let mesh = DgMesh::line(r.ne, 0.0, 1.0, r.p);
    let per = [DgBc::Periodic, DgBc::Periodic, DgBc::Periodic, DgBc::Periodic];
    let out = [DgBc::Outflow, DgBc::Outflow, DgBc::Outflow, DgBc::Outflow];
    let lim = r.limiter.then_some(DgLimiter { tvb_m: Some(20.0), positivity: true });
    let (law, bcs): (Box<dyn DgLaw>, &[DgBc; 4]) = match r.law.as_str() {
        "advection" => (Box::new(LinearAdvection { a: [1.0, 0.0] }), &per),
        "burgers" => (Box::new(Burgers), &per),
        "euler" => (Box::new(Euler1d { gamma: 1.4 }), &out),
        other => return Err(format!("unknown law `{other}` (advection, burgers, euler)")),
    };
    let init = |x: [f64; 2]| -> Vec<f64> {
        match r.law.as_str() {
            "advection" => vec![(2.0 * PI * x[0]).sin()],
            "burgers" => vec![0.5 + 0.3 * (2.0 * PI * x[0]).sin()],
            _ => {
                let (rho, p) = if x[0] < 0.5 { (1.0, 1.0) } else { (0.125, 0.1) };
                vec![rho, 0.0, p / 0.4]
            },
        }
    };
    let s = solve_dg_limited(&mesh, law.as_ref(), bcs, &init, r.t_end, 0.3, None, lim)?;
    let mut x = Vec::new();
    let mut u = Vec::new();
    for e in 0..mesh.elements() {
        for nd in 0..mesh.nodes_per_element() {
            x.push(s.coord(e, nd)[0]);
            u.push(s.value(e, nd, 0));
        }
    }
    let error = match r.law.as_str() {
        "advection" => Some(s.l2_error(0, |p| (2.0 * PI * (p[0] - r.t_end)).sin())),
        "euler" => {
            let (l, rr) = (Primitive { rho: 1.0, u: 0.0, p: 1.0 }, Primitive { rho: 0.125, u: 0.0, p: 0.1 });
            Some(s.l2_error(0, |p| exact_riemann(l, rr, 1.4, (p[0] - 0.5) / r.t_end).rho))
        },
        _ => None,
    };
    Ok(json!({ "x": x, "u": u, "l2_error": error, "integral": s.integral(0), "steps": s.steps }))
}

#[derive(Deserialize)]
#[serde(default)]
struct AcousticsReq {
    dim: usize,
    n: usize,
    pml_cells: usize,
    t_end: f64,
    sigma: f64,
    absorbing: bool,
}

impl Default for AcousticsReq {
    fn default() -> Self {
        Self { dim: 2, n: 80, pml_cells: 10, t_end: 1.4, sigma: 0.04, absorbing: true }
    }
}

fn acoustics_pml(params: &str) -> Result<Value, String> {
    use crate::sim::acoustics::{Acoustics, AcousticsConfig};
    let r: AcousticsReq = parse(params)?;
    check_dim(r.dim, 1, 3)?;
    let h = 1.0 / r.n as f64;
    let mut cfg = AcousticsConfig::new(&vec![r.n; r.dim], &vec![h; r.dim], 1.0)?;
    if r.absorbing {
        cfg.bcs = FaceBcs::uniform(&Bc::Absorbing);
        cfg.pml_cells = r.pml_cells;
    }
    let mut a = Acoustics::new(cfg)?;
    let (d, s) = (r.dim, r.sigma);
    a.set_pressure(move |x| (-(0..d).map(|k| (x[k] - 0.5).powi(2)).sum::<f64>() / (2.0 * s * s)).exp());
    let e0 = a.energy();
    a.run((r.t_end / a.dt).round() as usize);
    let nn = a.n();
    let mid = nn[0] / 2;
    let line: Vec<f64> = (0..nn[0]).map(|i| a.p[a.idx(i, if d > 1 { mid } else { 0 }, if d > 2 { mid } else { 0 })]).collect();
    Ok(json!({ "energy_ratio": a.energy() / e0, "line": line, "time": a.t }))
}

#[derive(Deserialize)]
#[serde(default)]
struct YeeReq {
    dim: usize,
    n: usize,
    steps: usize,
    pml: bool,
}

impl Default for YeeReq {
    fn default() -> Self {
        Self { dim: 2, n: 60, steps: 120, pml: true }
    }
}

fn yee_fdtd(params: &str) -> Result<Value, String> {
    use crate::sim::yee_fdtd::{Comp, CurrentSource, Yee, YeeConfig};
    let r: YeeReq = parse(params)?;
    check_dim(r.dim, 1, 3)?;
    let h = 1.0 / r.n as f64;
    let mut cfg = YeeConfig::new(&vec![r.n; r.dim], &vec![h; r.dim])?;
    if r.pml {
        cfg.bcs = FaceBcs::uniform(&Bc::Absorbing);
        cfg.pml_cells = (r.n / 8).max(4);
    }
    let mut y = Yee::new(cfg)?;
    let c = r.n / 2;
    y.add_source(CurrentSource {
        cell: [c, if r.dim > 1 { c } else { 0 }, if r.dim > 2 { c } else { 0 }],
        comp: Comp::Ez,
        waveform: Box::new(|t| {
            let z = (t - 0.45) / 0.12;
            -2.0 * z * (-z * z).exp()
        }),
    });
    let mut peak: f64 = 0.0;
    for _ in 0..r.steps {
        y.step();
        peak = peak.max(y.energy());
    }
    let mid = r.n / 2;
    let line: Vec<f64> = (0..r.n).map(|i| y.field(Comp::Ez)[y.idx(i, if r.dim > 1 { mid } else { 0 }, if r.dim > 2 { mid } else { 0 })]).collect();
    Ok(json!({ "dt": y.dt, "energy": y.energy(), "peak_energy": peak, "ez_line": line }))
}

#[derive(Deserialize)]
#[serde(default)]
struct NsReq {
    case: String,
    dim: usize,
    n: usize,
    re: f64,
    t_end: f64,
    dt: f64,
}

impl Default for NsReq {
    fn default() -> Self {
        Self { case: "cavity".into(), dim: 2, n: 32, re: 100.0, t_end: 5.0, dt: 0.01 }
    }
}

fn navier_stokes_mac(params: &str) -> Result<Value, String> {
    use crate::sim::ns_mac::{NsConfig, NsSolver};
    let r: NsReq = parse(params)?;
    check_dim(r.dim, 2, 3)?;
    let steps = (r.t_end / r.dt).round() as usize;
    let nu = 1.0 / r.re;
    let mut solver = match r.case.as_str() {
        "cavity" => {
            let mut cfg = NsConfig::new(&vec![r.n; r.dim], &vec![1.0; r.dim], nu)?;
            cfg.bcs = FluidFaceBcs::uniform(&FluidBoundary::wall()).with(Face::YMax, FluidBoundary::moving_wall([1.0, 0.0, 0.0]));
            NsSolver::new(cfg, r.dt, &|_| [0.0; 3])?
        },
        "taylor_green" => {
            let mut cfg = NsConfig::new(&vec![r.n; r.dim], &[vec![2.0 * PI; 2], if r.dim == 3 { vec![1.0] } else { vec![] }].concat(), nu)?;
            cfg.bcs = FluidFaceBcs::uniform(&FluidBoundary::Periodic);
            let init = |x: [f64; 3]| [-x[0].cos() * x[1].sin(), x[0].sin() * x[1].cos(), 0.0];
            NsSolver::new(cfg, r.dt, &init)?
        },
        "poiseuille" => {
            let mut cfg = NsConfig::new(&vec![r.n; r.dim], &vec![1.0; r.dim], nu)?;
            cfg.bcs = FluidFaceBcs::uniform(&FluidBoundary::wall()).with_axis(0, &FluidBoundary::Periodic);
            if r.dim == 3 {
                cfg.bcs = cfg.bcs.clone().with_axis(2, &FluidBoundary::Periodic);
            }
            cfg.force = Some(Box::new(|_, _| [1.0, 0.0, 0.0]));
            NsSolver::new(cfg, r.dt, &|_| [0.0; 3])?
        },
        other => return Err(format!("unknown case `{other}` (cavity, taylor_green, poiseuille)")),
    };
    solver.run(steps)?;
    let div = solver.divergence().iter().fold(0.0_f64, |m, d| m.max(d.abs()));
    let line: Vec<f64> = (0..r.n)
        .map(|k| {
            let y = (k as f64 + 0.5) / r.n as f64;
            let x = if r.case == "taylor_green" { PI } else { 0.5 };
            let yy = if r.case == "taylor_green" { 2.0 * PI * y } else { y };
            solver.sample_velocity(0, [x, yy, 0.5])
        })
        .collect();
    let decay = (r.case == "taylor_green").then(|| (-2.0 * nu * r.t_end).exp());
    Ok(json!({ "u_line": line, "max_divergence": div, "kinetic_energy": solver.kinetic_energy(), "steps": steps, "taylor_green_decay": decay }))
}

#[derive(Deserialize)]
#[serde(default)]
struct SplitReq {
    dim: usize,
    n: usize,
    len: f64,
    g: f64,
    mode: String,
    dt: f64,
    steps: usize,
}

impl Default for SplitReq {
    fn default() -> Self {
        Self { dim: 1, n: 64, len: 16.0, g: 0.0, mode: "ground_state".into(), dt: 0.05, steps: 2000 }
    }
}

fn split_step(params: &str) -> Result<Value, String> {
    use crate::sim::split_step::{SplitOrder, SplitStep};
    use num_complex::Complex;
    let r: SplitReq = parse(params)?;
    check_dim(r.dim, 1, 3)?;
    let s = SplitStep::centered(&vec![r.n; r.dim], &vec![r.len; r.dim], 1.0)?;
    let d = r.dim;
    let v = move |_t: f64, x: [f64; 3]| 0.5 * (0..d).map(|a| x[a] * x[a]).sum::<f64>();
    let guess = s.sample(|x| Complex::new((-(0..d).map(|a| x[a] * x[a]).sum::<f64>() / 8.0).exp(), 0.0));
    let (psi, energy, mu) = match r.mode.as_str() {
        "ground_state" => {
            let gs = s.ground_state(&guess, &v, r.g, r.dt, r.steps, 1e-12);
            (gs.psi, gs.energy, gs.mu)
        },
        "evolve" => {
            let shifted = s.sample(|x| Complex::new((-((x[0] - 1.0).powi(2) + (1..d).map(|a| x[a] * x[a]).sum::<f64>()) / 2.0).exp(), 0.0));
            let psi = s.evolve(&shifted, &v, r.g, 0.0, r.dt, r.steps, SplitOrder::Yoshida4);
            let e = s.energy(&psi, &v, r.dt * r.steps as f64, r.g);
            (psi, e, 0.0)
        },
        other => return Err(format!("unknown mode `{other}` (ground_state, evolve)")),
    };
    let mid = r.n / 2;
    let (j, k) = (if d > 1 { mid } else { 0 }, if d > 2 { mid } else { 0 });
    let line: Vec<f64> = (0..r.n).map(|i| psi[i + r.n * (j + r.n * k)].norm_sqr()).collect();
    Ok(json!({ "energy": energy, "mu": mu, "norm": s.norm2(&psi), "density_line": line }))
}

#[derive(Deserialize)]
#[serde(default)]
struct FemReq {
    dim: usize,
    degree: usize,
    n: usize,
}

impl Default for FemReq {
    fn default() -> Self {
        Self { dim: 2, degree: 1, n: 8 }
    }
}

fn fem_simplex(params: &str) -> Result<Value, String> {
    use crate::sim::fem_simplex::{FemProblem, FemSpace, SimplexMesh, solve_static};
    use crate::sim::structured::LinearSolver;
    let r: FemReq = parse(params)?;
    check_dim(r.dim, 1, 3)?;
    let mesh = match r.dim {
        1 => SimplexMesh::interval(r.n, 0.0, 1.0),
        2 => SimplexMesh::rectangle(r.n, r.n, [0.0; 2], [1.0; 2]),
        _ => SimplexMesh::boxed([r.n; 3], [0.0; 3], [1.0; 3]),
    };
    let space = FemSpace::new(mesh, r.degree)?;
    let d = r.dim;
    let exact = move |x: [f64; 3]| (0..d).map(|a| x[a]).sum::<f64>().exp();
    let bcs = MarkerBcs::with_default(Bc::Dirichlet(BcValue::space(exact)));
    let c = 1.0;
    let f = move |_t: f64, x: [f64; 3]| (c - d as f64) * exact(x);
    let prob = FemProblem { k: 1.0, c, f: &f, bcs: &bcs };
    let u = solve_static(&space, &prob, 0.0, LinearSolver::Auto)?;
    Ok(json!({ "dofs": space.size(), "l2_error": space.l2_error(&u, exact) }))
}

#[derive(Deserialize)]
#[serde(default)]
struct ElasticReq {
    dim: usize,
    degree: usize,
    length_cells: usize,
    young: f64,
    poisson: f64,
    load: f64,
}

impl Default for ElasticReq {
    fn default() -> Self {
        Self { dim: 2, degree: 2, length_cells: 20, young: 1000.0, poisson: 0.3, load: 1.0 }
    }
}

fn elasticity_fem(params: &str) -> Result<Value, String> {
    use crate::sim::boundary::ElasticBoundary;
    use crate::sim::elasticity_fem::{ElasticityProblem, Material, Model, solve};
    use crate::sim::fem_simplex::{FemSpace, SimplexMesh, marker};
    use crate::sim::structured::LinearSolver;
    let r: ElasticReq = parse(params)?;
    check_dim(r.dim, 2, 3)?;
    let (l, h) = (10.0, 1.0);
    let mat = Material { young: r.young, poisson: r.poisson, thickness: 1.0 };
    let (mesh, model, load) = if r.dim == 2 {
        (SimplexMesh::rectangle(r.length_cells, (r.length_cells / 10).max(2), [0.0; 2], [l, h]), Model::PlaneStress, [0.0, -r.load / h, 0.0])
    } else {
        (SimplexMesh::boxed([r.length_cells, 2, 2], [0.0; 3], [l, h, h]), Model::Solid3d, [0.0, 0.0, -r.load / (h * h)])
    };
    let bcs = MarkerBcs::new().set(marker::XMIN, ElasticBoundary::clamped()).set(marker::XMAX, ElasticBoundary::traction(load));
    let space = FemSpace::new(mesh, r.degree)?;
    let force = |_: [f64; 3]| [0.0; 3];
    let prob = ElasticityProblem { model, material: mat, force: &force, bcs: &bcs };
    let sol = solve(&space, &prob, LinearSolver::Auto)?;
    let tip = sol.eval([l, 0.5 * h, 0.5 * h]).ok_or("tip outside the mesh")?;
    let inertia = if r.dim == 2 { h.powi(3) / 12.0 } else { h.powi(4) / 12.0 };
    let beam = r.load * l.powi(3) / (3.0 * r.young * inertia);
    Ok(json!({ "tip_displacement": tip, "euler_bernoulli_deflection": beam, "dofs": space.size() }))
}

#[derive(Deserialize)]
#[serde(default)]
struct BemReq {
    subdivisions: usize,
    problem: String,
    method: String,
}

impl Default for BemReq {
    fn default() -> Self {
        Self { subdivisions: 2, problem: "interior".into(), method: "aca".into() }
    }
}

fn bem_laplace(params: &str) -> Result<Value, String> {
    use crate::sim::bem_aca::{AcaOptions, solve_laplace_bem_3d_aca, solve_laplace_bem_3d_exterior, solve_laplace_bem_3d_p1};
    use crate::sim::physics_bem::{BoundaryCondition as El, SurfaceMesh3D};
    let r: BemReq = parse(params)?;
    let mesh = SurfaceMesh3D::icosphere(1.0, [0.0; 3], r.subdivisions);
    let exterior = match r.problem.as_str() {
        "interior" => false,
        "exterior" => true,
        other => return Err(format!("unknown problem `{other}` (interior, exterior)")),
    };
    // interior: u = x, q = n_x ; exterior: u = 1 (capacitance), q = -1
    if r.method == "p1" {
        let bcs: Vec<El<f64>> = mesh.vertices.iter().map(|v| El::Potential(if exterior { 1.0 } else { v[0] })).collect();
        let sol = solve_laplace_bem_3d_p1(&mesh, &bcs, exterior)?;
        let err = mesh.vertices.iter().zip(&sol.q).map(|(v, q)| (q - if exterior { -1.0 } else { v[0] }).abs()).fold(0.0, f64::max);
        return Ok(json!({ "unknowns": mesh.vertices.len(), "max_flux_error": err }));
    }
    let bcs: Vec<El<f64>> = (0..mesh.triangles.len())
        .map(|t| El::Potential(if exterior { 1.0 } else { mesh.centroid(t).map_or(0.0, |c| c[0]) }))
        .collect();
    let opts = AcaOptions { leaf_size: 16, ..AcaOptions::default() };
    let (sol, stats) = if exterior { solve_laplace_bem_3d_exterior(&mesh, &bcs, &opts)? } else { solve_laplace_bem_3d_aca(&mesh, &bcs, &opts)? };
    let err = (0..mesh.triangles.len())
        .map(|t| (sol.q[t] - if exterior { -1.0 } else { mesh.normal(t).map_or(0.0, |n| n[0]) }).abs())
        .fold(0.0, f64::max);
    Ok(json!({
        "unknowns": mesh.triangles.len(),
        "max_flux_error": err,
        "gmres_iterations": stats.iterations,
        "compression": stats.compression(),
    }))
}

// ------------------------------------------------------------------ particles

#[derive(Deserialize)]
#[serde(default)]
struct SphReq {
    t_end: f64,
    grad_h: bool,
    gamma: f64,
}

impl Default for SphReq {
    fn default() -> Self {
        Self { t_end: 0.2, grad_h: false, gamma: 1.4 }
    }
}

fn sph(params: &str) -> Result<Value, String> {
    use crate::sim::sph::{Sph, SphBoundary};
    let r: SphReq = parse(params)?;
    let m = 0.001;
    let (mut pos, mut u) = (Vec::new(), Vec::new());
    let mut x = -0.5 + 0.0005;
    while x < 0.0 {
        pos.push([x]);
        u.push(1.0 / (r.gamma - 1.0));
        x += 0.001;
    }
    let mut x = 0.004;
    while x < 0.5 {
        pos.push([x]);
        u.push(0.1 / ((r.gamma - 1.0) * 0.125));
        x += 0.008;
    }
    let n = pos.len();
    let mut s = Sph::<1>::new(pos, vec![m; n], u, [-0.5], [0.5], [SphBoundary::Reflecting]);
    s.gamma = r.gamma;
    s.grad_h = r.grad_h;
    s.update_density();
    let e0 = s.total_energy();
    let mut t = 0.0;
    while t < r.t_end {
        let dt = s.suggest_dt(0.25).min(r.t_end - t);
        s.step(dt);
        t += dt;
    }
    Ok(json!({
        "x": s.pos.iter().map(|p| p[0]).collect::<Vec<_>>(),
        "rho": s.rho,
        "pressure": (0..n).map(|i| s.pressure(i)).collect::<Vec<_>>(),
        "energy_drift": (s.total_energy() - e0) / e0,
    }))
}

#[derive(Deserialize)]
#[serde(default)]
struct DpdReq {
    dim: usize,
    box_len: f64,
    density: f64,
    steps: usize,
    dt: f64,
    a: f64,
    seed: u64,
}

impl Default for DpdReq {
    fn default() -> Self {
        Self { dim: 3, box_len: 6.0, density: 3.0, steps: 400, dt: 0.04, a: 25.0, seed: 1 }
    }
}

fn dpd(params: &str) -> Result<Value, String> {
    use crate::sim::dpd::{Dpd, DpdParams};
    let r: DpdReq = parse(params)?;
    check_dim(r.dim, 2, 3)?;
    let n = (r.density * r.box_len.powi(r.dim as i32)).round() as usize;
    let p = DpdParams { dim: r.dim, box_len: [r.box_len; 3], a: r.a, gamma: 4.5, kt: 1.0, rc: 1.0, mass: 1.0 };
    let mut d = Dpd::new(n, p, r.seed)?;
    d.run(r.dt, r.steps);
    Ok(json!({ "particles": n, "temperature": d.temperature(), "pressure": d.pressure(), "momentum": d.momentum() }))
}

#[derive(Deserialize)]
#[serde(default)]
struct MdReq {
    dim: usize,
    cells: usize,
    density: f64,
    temperature: f64,
    steps: usize,
    dt: f64,
    thermostat: String,
    seed: u64,
}

impl Default for MdReq {
    fn default() -> Self {
        Self { dim: 3, cells: 6, density: 0.8, temperature: 1.0, steps: 200, dt: 0.004, thermostat: "nve".into(), seed: 1 }
    }
}

fn md_run<const D: usize>(r: &MdReq) -> Result<Value, String> {
    use crate::kernels::random::Rng;
    use crate::sim::md_engine::{ForceField, MdSimulation, MdSystem, PairTerm, SimBox, Thermostat, lattice_cubic};
    let spacing = (1.0 / r.density).powf(1.0 / D as f64);
    let len = spacing * r.cells as f64;
    let mut sys = MdSystem::new(lattice_cubic([r.cells; D], spacing), SimBox::periodic([len; D]));
    let mut rng = Rng::new(r.seed);
    sys.maxwell_boltzmann(r.temperature, &mut rng);
    let mut ff = ForceField::new(vec![PairTerm::LennardJones { eps: 1.0, sigma: 1.0 }], 2.5);
    ff.shift_force = true;
    let mut sim = MdSimulation::new(sys, ff, 0.3, r.dt, r.seed + 1)?;
    sim.thermostat = match r.thermostat.as_str() {
        "nve" => Thermostat::None,
        "langevin" => Thermostat::Langevin { t0: r.temperature, gamma: 1.0 },
        other => return Err(format!("unknown thermostat `{other}` (nve, langevin)")),
    };
    let e0 = sim.total_energy();
    sim.run(r.steps)?;
    Ok(json!({
        "particles": sim.sys.n(),
        "temperature": sim.temperature(),
        "total_energy": sim.total_energy(),
        "energy_drift": (sim.total_energy() - e0) / e0.abs().max(1e-12),
        "pressure": sim.pressure(),
    }))
}

fn md_lj(params: &str) -> Result<Value, String> {
    let r: MdReq = parse(params)?;
    check_dim(r.dim, 2, 3)?;
    if r.dim == 2 { md_run::<2>(&r) } else { md_run::<3>(&r) }
}

#[derive(Deserialize)]
#[serde(default)]
struct SpinReq {
    dim: usize,
    n: usize,
    beta: f64,
    sweeps: usize,
    model: String,
    update: String,
    seed: u64,
}

impl Default for SpinReq {
    fn default() -> Self {
        Self { dim: 2, n: 16, beta: 0.5, sweeps: 200, model: "ising".into(), update: "wolff".into(), seed: 1 }
    }
}

fn spin_run<const D: usize, const N: usize>(r: &SpinReq) -> Result<Value, String> {
    use crate::sim::spin_models::{OnModel, Update};
    let mut m = OnModel::<D, N>::new([r.n; D], true, r.beta, r.seed)?;
    let upd = match r.update.as_str() {
        "wolff" => Update::Wolff,
        "metropolis" => Update::Metropolis(0.5),
        "heat_bath" => Update::HeatBath,
        "swendsen_wang" => Update::SwendsenWang,
        other => return Err(format!("unknown update `{other}`")),
    };
    for _ in 0..r.sweeps {
        m.update(upd)?;
    }
    let sites = m.lattice.sites() as f64;
    Ok(json!({ "energy_per_site": m.energy() / sites, "abs_magnetisation_per_site": m.abs_magnetisation_per_site(), "sites": m.lattice.sites() }))
}

fn spin_model(params: &str) -> Result<Value, String> {
    let r: SpinReq = parse(params)?;
    check_dim(r.dim, 1, 3)?;
    match (r.model.as_str(), r.dim) {
        ("ising", 1) => spin_run::<1, 1>(&r),
        ("ising", 2) => spin_run::<2, 1>(&r),
        ("ising", _) => spin_run::<3, 1>(&r),
        ("xy", 1) => spin_run::<1, 2>(&r),
        ("xy", 2) => spin_run::<2, 2>(&r),
        ("xy", _) => spin_run::<3, 2>(&r),
        ("heisenberg", 1) => spin_run::<1, 3>(&r),
        ("heisenberg", 2) => spin_run::<2, 3>(&r),
        ("heisenberg", _) => spin_run::<3, 3>(&r),
        (other, _) => Err(format!("unknown model `{other}` (ising, xy, heisenberg)")),
    }
}

#[derive(Deserialize)]
struct NbodyReq {
    positions: Vec<[f64; 3]>,
    masses: Vec<f64>,
    #[serde(default = "one")]
    g: f64,
    #[serde(default = "small")]
    softening: f64,
    #[serde(default = "theta_default")]
    theta: f64,
}

const fn one() -> f64 {
    1.0
}
const fn small() -> f64 {
    1e-3
}
const fn theta_default() -> f64 {
    0.5
}

fn nbody_gravity(params: &str) -> Result<Value, String> {
    use crate::sim::nbody::gravity_accelerations;
    let r: NbodyReq = parse(params)?;
    let acc = gravity_accelerations(&r.positions, &r.masses, r.g, r.softening, r.theta)?;
    Ok(json!({ "accelerations": acc }))
}

#[derive(Deserialize)]
#[serde(default)]
struct KerrReq {
    mass: f64,
    spin: f64,
    charge: f64,
    radius: f64,
    prograde: bool,
    step: f64,
    steps: usize,
}

impl Default for KerrReq {
    fn default() -> Self {
        Self { mass: 1.0, spin: 0.5, charge: 0.0, radius: 10.0, prograde: true, step: 0.5, steps: 400 }
    }
}

fn kerr_geodesic(params: &str) -> Result<Value, String> {
    use crate::sim::relativity::{Geodesic, KerrNewman};
    let r: KerrReq = parse(params)?;
    let metric = KerrNewman { m: r.mass, a: r.spin, q: r.charge };
    let mut g = Geodesic::circular(metric, r.radius, r.prograde)?;
    let horizon = metric.horizon().unwrap_or(0.0);
    let end = g.integrate(r.step, r.steps, true, horizon * 1.001, 1.0e6)?;
    Ok(json!({
        "termination": format!("{end:?}"),
        "r": g.r,
        "phi": g.phi,
        "t": g.t,
        "constraint": g.constraint(),
        "energy": g.energy,
        "angular_momentum": g.lz,
    }))
}
