//! # Incompressible Navier–Stokes on a staggered (MAC) grid in 2D and 3D
//!
//! Second-order incremental pressure-correction (projection) solver for
//! `u_t + (u·∇)u = -∇p + ν Δu + f`, `∇·u = 0`:
//!
//! * staggered marker-and-cell grid (`u, v, w` on faces, `p` in cells);
//! * convection: conservative second-order central differences, advanced
//!   with Adams–Bashforth 2 (Euler in the first step);
//! * diffusion: Crank–Nicolson, solved implicitly for each component with a
//!   sparse LU factorisation computed once;
//! * pressure: incremental projection `Δφ = ∇·u*/Δt`, `u = u* - Δt ∇φ`,
//!   `p += φ` (optionally the rotational correction `- ν/2 ∇·u*`);
//! * boundary faces from [`FluidFaceBcs`]: no-slip (fixed or moving wall),
//!   slip / symmetry, prescribed inflow velocity (time dependent), outflow
//!   with prescribed pressure, periodic;
//! * obstacles: a Boolean cell mask (staircase no-slip, first-order at the
//!   body) or an **immersed body** given by a level set ([`ImmersedBody`]):
//!   Shortley–Weller viscous stencils with the exact wall distance, aperture
//!   weighted divergence and pressure Poisson on cut cells (second order for
//!   diffusion-dominated flow, see the Taylor–Couette test);
//! * **cut cells** ([`CutCellOptions`]): face apertures and cell volume
//!   fractions come from [`cutcell`](crate::sim::cutcell) (exact clipping of
//!   the level set, also in 3D). The convective term near the body is
//!   evaluated in conservative form on the face-centred dual cells with their
//!   own apertures and volume fractions (face values interpolated at the true
//!   face position between the node and the wall, moving rigid bodies
//!   included), and cut cells with a volume fraction below `beta_min` are
//!   stabilised by **flux redistribution**: the cell takes the stable
//!   non-conservative update and the conservation defect is spread over its
//!   neighbours in proportion to their volumes, so momentum stays conserved
//!   and the time step is not restricted by small cells;
//! * **convection schemes** ([`Convection`]): second-order central, QUICK or
//!   flux-limited TVD upwind (minmod, van Leer, superbee, MC) in flux form;
//! * **turbulence model** ([`LesModel`]): Smagorinsky or WALE eddy viscosity
//!   added explicitly as `∇·(ν_t (∇u + ∇uᵀ))`;
//! * **buoyancy** ([`Thermal`]): a cell-centred temperature advected and
//!   diffused (Crank–Nicolson) with a Boussinesq force
//!   `g β (T - T_ref)` on the momentum equation (Rayleigh–Bénard onset is
//!   reproduced, see the tests).
//!
//! | Feature | 2D | 3D |
//! |---|---|---|
//! | Navier–Stokes, constant density (Boussinesq buoyancy) | yes (`nz = 1`) | yes |
//! | no-slip / moving wall, slip, inflow, outflow, periodic | yes | yes |
//! | obstacles (cell mask) | yes | yes |
//! | immersed level-set body, exact apertures, flux redistribution | yes | yes |
//! | central / QUICK / TVD convection | yes | yes |
//! | Smagorinsky / WALE | yes | yes |
//! | body force `f(t, x)` | yes | yes |
//! | drag by momentum balance ([`NsSolver::force_on_region`]) | yes | yes |
//!
//! Limitations: uniform grids, no variable density (Boussinesq only), no
//! free surface, the thermal field needs a body-free domain, the
//! eddy-viscosity stress is not wall-modelled and is evaluated on the plain
//! grid stencils near immersed bodies, and the high-resolution convection
//! schemes fall back to central differences next to walls and bodies.
//!
//! Validation: the immersed-body solver reproduces the `Re = 100` flow past a
//! circular cylinder (`D = 1`, `U_inf = 1`).  A 40 x 20 box (5 % blockage)
//! with 8 cells per diameter gives `C_D = 1.399` and `St = 0.1653` after six
//! shedding periods (literature `C_D = 1.32-1.34`, `St = 0.164-0.166`; Braza,
//! Chassaing & Ha Minh 1986, Williamson 1989), i.e. +5.2 % and +0.2 %; a
//! 20 x 10 box at 12-16 cells per diameter gives `C_D = 1.49` and
//! `St = 0.178`, the difference being blockage.  See `tests/sim/ns_mac.rs`.

use crate::kernels::krylov::Csr;
use crate::sim::boundary::Face;
use crate::sim::boundary::FluidBoundary;
use crate::sim::boundary::FluidFaceBcs;
use crate::sim::cutcell;
use crate::sim::structured::LinearSolver;
use crate::sim::structured::LinearSystem;

/// Body-force closure `f(t, x)`.
pub type ForceFn = dyn Fn(f64, [f64; 3]) -> [f64; 3] + Send + Sync;

/// A rigid body immersed in the flow, given by a level set.
///
/// The surface is `level_set = 0`, with `level_set < 0` inside the body.
/// Velocity components inside the body take the (steady) rigid-body velocity
/// `velocity(x)`; the viscous operator uses Shortley–Weller stencils with the
/// exact wall distance along each axis (second order for diffusion), and the
/// pressure Poisson equation and the divergence use the open fraction
/// (aperture) of every cell face, which keeps the projection conservative on
/// cut cells.
pub struct ImmersedBody {
    /// Signed distance-like function, negative inside the body.
    pub level_set: Box<dyn Fn([f64; 3]) -> f64 + Send + Sync>,
    /// Rigid-body velocity at a point (zero for a fixed body).
    pub velocity: Box<dyn Fn([f64; 3]) -> [f64; 3] + Send + Sync>,
}

/// Flux limiter of the TVD convection scheme.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ConvLimiter {
    /// Minmod.
    Minmod,
    /// Van Leer.
    VanLeer,
    /// Superbee.
    Superbee,
    /// Monotonised central.
    Mc,
}

impl ConvLimiter {
    fn psi(
        self,
        r: f64,
    ) -> f64 {
        match self {
            | Self::Minmod => r.clamp(0.0, 1.0),
            | Self::VanLeer => (r + r.abs()) / (1.0 + r.abs()),
            | Self::Superbee => (2.0 * r).min(1.0).max(r.min(2.0)).max(0.0),
            | Self::Mc => f64::midpoint(1.0, r).min(2.0).min(2.0 * r).max(0.0),
        }
    }
}

/// Face reconstruction of the convective flux.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Convection {
    /// Second-order central (arithmetic mean of the two nodes).
    Central,
    /// QUICK (third-order upwind-biased quadratic interpolation).
    Quick,
    /// Flux-limited second-order upwind (total variation diminishing).
    Tvd(ConvLimiter),
}

/// Sub-grid-scale eddy-viscosity model.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum LesModel {
    /// No model.
    None,
    /// Smagorinsky: `ν_t = (C_s Δ)² |S|`, `|S| = sqrt(2 S_ij S_ij)`.
    Smagorinsky {
        /// Smagorinsky constant (about 0.1–0.17).
        cs: f64,
    },
    /// Wall-adapting local eddy viscosity (Nicoud–Ducros).
    Wale {
        /// WALE constant (about 0.325).
        cw: f64,
    },
}

/// Boundary condition of the temperature on one domain face.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ScalarBc {
    /// Zero normal gradient.
    Adiabatic,
    /// Fixed value.
    Dirichlet(f64),
    /// Periodic (both faces of the axis).
    Periodic,
}

/// Boussinesq temperature field: `T_t + u·∇T = κ ΔT`, buoyancy `g β (T - T_ref)`.
pub struct Thermal {
    /// Thermal diffusivity.
    pub kappa: f64,
    /// Buoyancy acceleration per unit temperature excess, `g β` times the
    /// direction pointing **against** gravity.
    pub gbeta: [f64; 3],
    /// Reference temperature of the Boussinesq force.
    pub t_ref: f64,
    /// Face conditions `XMin, XMax, YMin, YMax, ZMin, ZMax`.
    pub bcs: [ScalarBc; 6],
    /// Initial temperature `T(x)`.
    pub init: Box<dyn Fn([f64; 3]) -> f64 + Send + Sync>,
}

/// Cut-cell treatment of the convective term near an immersed body.
#[derive(Clone, Copy, Debug)]
pub struct CutCellOptions {
    /// Conservative dual-cell fluxes with apertures (otherwise the plain
    /// extrapolated-ghost form is used).
    pub conservative: bool,
    /// Cut cells with a dual volume fraction below this value are stabilised
    /// by flux redistribution (`0` disables it).
    pub beta_min: f64,
}

impl Default for CutCellOptions {
    fn default() -> Self {
        Self {
            conservative: true,
            beta_min: 0.5,
        }
    }
}

/// Solver configuration.
pub struct NsConfig {
    /// Cells per axis (`n[2] = 1` for 2D).
    pub n: [usize; 3],
    /// Dimension (2 or 3).
    pub dim: usize,
    /// Domain size per axis.
    pub len: [f64; 3],
    /// Kinematic viscosity.
    pub nu: f64,
    /// Face boundary conditions.
    pub bcs: FluidFaceBcs,
    /// Solid cells (`true` = solid), x fastest; `None` for no obstacle.
    pub solid: Option<Vec<bool>>,
    /// Optional body force.
    pub force: Option<Box<ForceFn>>,
    /// Optional immersed body (level-set description, second order for diffusion).
    pub body: Option<ImmersedBody>,
    /// Use the rotational pressure update.
    pub rotational: bool,
    /// Linear solver for the pressure system.
    pub solver: LinearSolver,
    /// Convective flux reconstruction.
    pub convection: Convection,
    /// Sub-grid-scale model.
    pub les: LesModel,
    /// Optional Boussinesq temperature field.
    pub thermal: Option<Thermal>,
    /// Cut-cell options (used with an immersed body).
    pub cut_cell: CutCellOptions,
}

impl NsConfig {
    /// Creates a configuration with all walls no-slip.
    ///
    /// # Errors
    /// For bad dimensions or sizes.
    pub fn new(
        n: &[usize],
        len: &[f64],
        nu: f64,
    ) -> Result<Self, String> {
        let dim = n.len();
        if !(2..=3).contains(&dim) || len.len() != dim || n.iter().any(|&m| m < 3) || nu <= 0.0 {
            return Err("invalid Navier-Stokes grid or viscosity".into());
        }
        let mut c = Self {
            n: [1; 3],
            dim,
            len: [1.0; 3],
            nu,
            bcs: FluidFaceBcs::uniform(&FluidBoundary::wall()),
            solid: None,
            force: None,
            body: None,
            rotational: true,
            solver: LinearSolver::Auto,
            convection: Convection::Central,
            les: LesModel::None,
            thermal: None,
            cut_cell: CutCellOptions::default(),
        };
        c.n[..dim].copy_from_slice(&n[..dim]);
        c.len[..dim].copy_from_slice(&len[..dim]);
        Ok(c)
    }
}

#[derive(Clone, Copy, PartialEq)]
enum Status {
    Unknown,
    Fixed,
}

enum Term {
    Const(f64),
    Fixed(usize),
    Wall { face: Face, pos: [f64; 3] },
}

/// Smallest wall-distance fraction kept in a stencil; closer faces are snapped onto the wall.
const THETA_MIN: f64 = 0.1;

struct Ib {
    inside: [Vec<bool>; 3],
    inside_val: [Vec<f64>; 3],
    cross: [Vec<bool>; 3],
    theta: [Vec<f64>; 3],
    wall: [Vec<f64>; 3],
    aperture: [Vec<f64>; 3],
    /// Face is inside the body or touches it (high-resolution stencils fall back there).
    rough: [Vec<bool>; 3],
    /// Fluid fraction of the face-centred dual cell.
    dual_beta: [Vec<f64>; 3],
    /// Apertures of the dual-cell faces, six slots per face (`2 a + side`).
    dual_ap: [Vec<f64>; 3],
    /// Fluid fraction of every pressure cell.
    vol_frac: Vec<f64>,
    /// The body surface carries a non-zero velocity.
    moving: bool,
}

struct Helm {
    sys: LinearSystem,
    rows: Vec<usize>, // face index of each unknown row
    terms: Vec<(usize, f64, Term)>,
}

/// Incompressible flow solver state.
pub struct NsSolver {
    cfg: NsConfig,
    /// Time step.
    pub dt: f64,
    /// Current time.
    pub t: f64,
    /// Staggered velocity components (`u`, `v`, `w`; `w` empty in 2D).
    pub u: [Vec<f64>; 3],
    /// Cell pressure.
    pub p: Vec<f64>,
    status: [Vec<Status>; 3],
    bvals: [Vec<f64>; 3],
    helm: Vec<Helm>,
    pmat: LinearSystem,
    pfluid: Vec<bool>,
    ib: Option<Ib>,
    conv_prev: Option<[Vec<f64>; 3]>,
    periodic: [bool; 3],
    steps: usize,
    th: Option<ThermalState>,
}

impl NsSolver {
    const fn ext(
        cfg: &NsConfig,
        c: usize,
    ) -> [usize; 3] {
        let mut e = cfg.n;
        if c < cfg.dim {
            e[c] += 1;
        }
        e
    }

    const fn fidx(
        e: [usize; 3],
        i: [usize; 3],
    ) -> usize {
        i[0] + e[0] * (i[1] + e[1] * i[2])
    }

    const fn fijk(
        e: [usize; 3],
        p: usize,
    ) -> [usize; 3] {
        [p % e[0], (p / e[0]) % e[1], p / (e[0] * e[1])]
    }

    fn cell_solid(
        cfg: &NsConfig,
        i: [isize; 3],
    ) -> bool {
        for a in 0..3 {
            if i[a] < 0 || i[a] >= cfg.n[a] as isize {
                return false;
            }
        }
        cfg.solid
            .as_ref()
            .is_some_and(|s| s[Self::fidx(cfg.n, [i[0] as usize, i[1] as usize, i[2] as usize])])
    }

    fn face_pos(
        cfg: &NsConfig,
        c: usize,
        i: [usize; 3],
    ) -> [f64; 3] {
        let mut x = [0.0; 3];
        for a in 0..cfg.dim {
            let h = cfg.len[a] / cfg.n[a] as f64;
            x[a] = (i[a] as f64 + if a == c { 0.0 } else { 0.5 }) * h;
        }
        x
    }

    fn build_ib(cfg: &NsConfig) -> Option<Ib> {
        let body = cfg.body.as_ref()?;
        let dim = cfg.dim;
        let hs: Vec<f64> = (0..dim).map(|a| cfg.len[a] / cfg.n[a] as f64).collect();
        let mut ib = Ib {
            inside: std::array::from_fn(|_| Vec::new()),
            inside_val: std::array::from_fn(|_| Vec::new()),
            cross: std::array::from_fn(|_| Vec::new()),
            theta: std::array::from_fn(|_| Vec::new()),
            wall: std::array::from_fn(|_| Vec::new()),
            aperture: std::array::from_fn(|_| Vec::new()),
            rough: std::array::from_fn(|_| Vec::new()),
            dual_beta: std::array::from_fn(|_| Vec::new()),
            dual_ap: std::array::from_fn(|_| Vec::new()),
            vol_frac: Vec::new(),
            moving: false,
        };
        let phi = |x: [f64; 3]| (body.level_set)(x);
        for c in 0..dim {
            let e = Self::ext(cfg, c);
            let total: usize = e.iter().product();
            let mut inside = vec![false; total];
            let mut inside_val = vec![0.0; total];
            let mut cross = vec![false; 6 * total];
            let mut theta = vec![1.0; 6 * total];
            let mut wall = vec![0.0; 6 * total];
            let mut aperture = vec![1.0; total];
            for p in 0..total {
                let i = Self::fijk(e, p);
                let x = Self::face_pos(cfg, c, i);
                if phi(x) <= 0.0 {
                    inside[p] = true;
                    inside_val[p] = (body.velocity)(x)[c];
                } else {
                    let mut snap: Option<f64> = None;
                    for a in 0..dim {
                        for s in 0..2 {
                            let dir = if s == 0 { -1.0 } else { 1.0 };
                            let mut xn = x;
                            xn[a] += dir * hs[a];
                            if xn[a] < -1e-12 || xn[a] > cfg.len[a] + 1e-12 || phi(xn) > 0.0 {
                                continue;
                            }
                            let (mut lo, mut hi) = (0.0_f64, 1.0_f64);
                            for _ in 0..50 {
                                let mid = f64::midpoint(lo, hi);
                                let mut xm = x;
                                xm[a] += dir * mid * hs[a];
                                if phi(xm) > 0.0 {
                                    lo = mid;
                                } else {
                                    hi = mid;
                                }
                            }
                            let th = f64::midpoint(lo, hi);
                            let mut xw = x;
                            xw[a] += dir * th * hs[a];
                            let uw = (body.velocity)(xw)[c];
                            if th < THETA_MIN {
                                snap = Some(uw);
                            } else {
                                cross[6 * p + 2 * a + s] = true;
                                theta[6 * p + 2 * a + s] = th;
                                wall[6 * p + 2 * a + s] = uw;
                            }
                        }
                    }
                    if let Some(v) = snap {
                        inside[p] = true;
                        inside_val[p] = v;
                        for k in 0..6 {
                            cross[6 * p + k] = false;
                            theta[6 * p + k] = 1.0;
                        }
                    }
                }
                // aperture of this face (cell face normal to axis c): exact clipping of the level set
                let others: Vec<usize> = (0..dim).filter(|&b| b != c).collect();
                let mut lo = x;
                let mut len = [0.0; 3];
                for &b in &others {
                    lo[b] -= 0.5 * hs[b];
                    len[b] = hs[b];
                }
                aperture[p] = cutcell::open_fraction(
                    &phi,
                    lo,
                    len,
                    &others,
                    cutcell::default_depth(others.len()),
                );
            }
            let rough: Vec<bool> = (0..total)
                .map(|p| inside[p] || (0..6).any(|k| cross[6 * p + k]))
                .collect();
            if inside_val
                .iter()
                .chain(wall.iter())
                .any(|v| v.abs() > 1e-14)
            {
                ib.moving = true;
            }
            let (mut dual_beta, mut dual_ap) = (vec![1.0; total], vec![1.0; 6 * total]);
            if cfg.cut_cell.conservative {
                let axes_all: Vec<usize> = (0..dim).collect();
                for p in 0..total {
                    let x = Self::face_pos(cfg, c, Self::fijk(e, p));
                    let mut lo = x;
                    let mut len = [0.0; 3];
                    for a in 0..dim {
                        lo[a] -= 0.5 * hs[a];
                    }
                    len[..dim].copy_from_slice(&hs[..dim]);
                    dual_beta[p] = cutcell::open_fraction(
                        &phi,
                        lo,
                        len,
                        &axes_all,
                        cutcell::default_depth(dim),
                    );
                    for a in 0..dim {
                        let others: Vec<usize> = (0..dim).filter(|&b| b != a).collect();
                        for s in 0..2 {
                            let mut lo = x;
                            let mut len = [0.0; 3];
                            for &b in &others {
                                lo[b] -= 0.5 * hs[b];
                                len[b] = hs[b];
                            }
                            lo[a] += if s == 0 {
                                -0.5 * hs[a]
                            } else {
                                0.5 * hs[a]
                            };
                            dual_ap[6 * p + 2 * a + s] = cutcell::open_fraction(
                                &phi,
                                lo,
                                len,
                                &others,
                                cutcell::default_depth(others.len()),
                            );
                        }
                    }
                }
            }
            ib.rough[c] = rough;
            ib.dual_beta[c] = dual_beta;
            ib.dual_ap[c] = dual_ap;
            ib.inside[c] = inside;
            ib.inside_val[c] = inside_val;
            ib.cross[c] = cross;
            ib.theta[c] = theta;
            ib.wall[c] = wall;
            ib.aperture[c] = aperture;
        }
        let ncell: usize = cfg.n.iter().product();
        let axes_all: Vec<usize> = (0..dim).collect();
        ib.vol_frac = (0..ncell)
            .map(|p| {
                let i = Self::fijk(cfg.n, p);
                let mut lo = [0.0; 3];
                let mut len = [0.0; 3];
                for a in 0..dim {
                    lo[a] = i[a] as f64 * hs[a];
                    len[a] = hs[a];
                }
                cutcell::open_fraction(&phi, lo, len, &axes_all, cutcell::default_depth(dim))
            })
            .collect();
        Some(ib)
    }

    fn ib_theta(
        &self,
        c: usize,
        p: usize,
        a: usize,
    ) -> (f64, f64) {
        self.ib.as_ref().map_or((1.0, 1.0), |ib| {
            (ib.theta[c][6 * p + 2 * a], ib.theta[c][6 * p + 2 * a + 1])
        })
    }

    fn ib_cross(
        &self,
        c: usize,
        p: usize,
        a: usize,
        s: usize,
    ) -> Option<f64> {
        let ib = self.ib.as_ref()?;
        ib.cross[c][6 * p + 2 * a + s].then(|| ib.wall[c][6 * p + 2 * a + s])
    }

    fn aperture(
        &self,
        a: usize,
        i: [usize; 3],
        side: isize,
    ) -> f64 {
        let Some(ib) = &self.ib else { return 1.0 };
        let mut f = i;
        if side > 0 {
            f[a] += 1;
        }
        ib.aperture[a][Self::fidx(Self::ext(&self.cfg, a), f)]
    }

    /// Creates a solver with the given initial velocity (faces) and time step.
    ///
    /// # Errors
    /// For inconsistent periodic faces or singular systems.
    pub fn new(
        cfg: NsConfig,
        dt: f64,
        init: &dyn Fn([f64; 3]) -> [f64; 3],
    ) -> Result<Self, String> {
        let dim = cfg.dim;
        let mut periodic = [false; 3];
        for a in 0..dim {
            let lo = matches!(
                cfg.bcs.get(Face::of_axis(a, false)),
                FluidBoundary::Periodic
            );
            let hi = matches!(cfg.bcs.get(Face::of_axis(a, true)), FluidBoundary::Periodic);
            if lo != hi {
                return Err("periodic must be set on both faces of an axis".into());
            }
            periodic[a] = lo;
        }
        if let Some(s) = &cfg.solid
            && s.len() != cfg.n.iter().product::<usize>()
        {
            return Err("solid mask has the wrong length".into());
        }
        let ib = Self::build_ib(&cfg);
        let mut status: [Vec<Status>; 3] = std::array::from_fn(|_| Vec::new());
        let mut bvals: [Vec<f64>; 3] = std::array::from_fn(|_| Vec::new());
        let mut u: [Vec<f64>; 3] = std::array::from_fn(|_| Vec::new());
        for c in 0..dim {
            let e = Self::ext(&cfg, c);
            let total: usize = e.iter().product();
            let mut st = vec![Status::Unknown; total];
            let mut uu = vec![0.0; total];
            for p in 0..total {
                let i = Self::fijk(e, p);
                let x = Self::face_pos(&cfg, c, i);
                uu[p] = init(x)[c];
                let mut fixed = false;
                // domain boundary normal faces
                if i[c] == 0 || i[c] == cfg.n[c] {
                    let max = i[c] == cfg.n[c];
                    match cfg.bcs.get(Face::of_axis(c, max)) {
                        | FluidBoundary::Outflow(_) => {},
                        | FluidBoundary::Periodic => {
                            if max {
                                fixed = true; // alias of face 0
                            }
                        },
                        | _ => fixed = true,
                    }
                }
                // solid cells on either side
                let mut ii = [i[0] as isize, i[1] as isize, i[2] as isize];
                let cell_a = ii;
                ii[c] -= 1;
                if Self::cell_solid(&cfg, cell_a) || Self::cell_solid(&cfg, ii) {
                    fixed = true;
                }
                if ib.as_ref().is_some_and(|b| b.inside[c][p]) {
                    fixed = true;
                }
                st[p] = if fixed {
                    Status::Fixed
                } else {
                    Status::Unknown
                };
            }
            status[c] = st;
            bvals[c] = vec![0.0; total];
            u[c] = uu;
        }
        let ncell: usize = cfg.n.iter().product();
        let mut pfluid = vec![true; ncell];
        if let Some(s) = &cfg.solid {
            for p in 0..ncell {
                pfluid[p] = !s[p];
            }
        }
        if let Some(b) = &ib {
            for p in 0..ncell {
                let i = Self::fijk(cfg.n, p);
                let open = (0..cfg.dim).any(|a| {
                    let e = Self::ext(&cfg, a);
                    let mut j = i;
                    j[a] += 1;
                    b.aperture[a][Self::fidx(e, i)] > 1e-12
                        || b.aperture[a][Self::fidx(e, j)] > 1e-12
                });
                if !open {
                    pfluid[p] = false;
                }
            }
        }
        let mut solver = Self {
            cfg,
            dt,
            t: 0.0,
            u,
            p: vec![0.0; ncell],
            status,
            bvals,
            helm: Vec::new(),
            pmat: LinearSystem::new(
                Csr::from_triplets(1, 1, &[(0, 0, 1.0)]),
                LinearSolver::Direct,
            )?,
            pfluid,
            ib,
            conv_prev: None,
            periodic,
            steps: 0,
            th: None,
        };
        solver.th = Self::build_thermal(&solver.cfg, dt)?;
        solver.update_boundary_values(0.0);
        solver.apply_fixed();
        solver.build_helmholtz()?;
        solver.build_poisson()?;
        Ok(solver)
    }

    fn h(
        &self,
        a: usize,
    ) -> f64 {
        self.cfg.len[a] / self.cfg.n[a] as f64
    }

    /// Boundary value of a fixed face at time `t`.
    fn bc_value(
        &self,
        c: usize,
        i: [usize; 3],
        t: f64,
    ) -> f64 {
        let cfg = &self.cfg;
        if let Some(ib) = &self.ib {
            let pf = Self::fidx(Self::ext(cfg, c), i);
            if ib.inside[c][pf] {
                return ib.inside_val[c][pf];
            }
        }
        let mut ii = [i[0] as isize, i[1] as isize, i[2] as isize];
        let cell_a = ii;
        ii[c] -= 1;
        if Self::cell_solid(cfg, cell_a) || Self::cell_solid(cfg, ii) {
            return 0.0;
        }
        if i[c] == 0 || i[c] == cfg.n[c] {
            let max = i[c] == cfg.n[c];
            let x = Self::face_pos(cfg, c, i);
            return match cfg.bcs.get(Face::of_axis(c, max)) {
                | FluidBoundary::NoSlip(v) | FluidBoundary::Inflow(v) => v[c].eval(t, x),
                | _ => 0.0,
            };
        }
        0.0
    }

    fn update_boundary_values(
        &mut self,
        t: f64,
    ) {
        for c in 0..self.cfg.dim {
            let e = Self::ext(&self.cfg, c);
            for p in 0..self.bvals[c].len() {
                if self.status[c][p] == Status::Fixed {
                    let i = Self::fijk(e, p);
                    let v = if self.periodic[c] && i[c] == self.cfg.n[c] {
                        0.0
                    } else {
                        self.bc_value(c, i, t)
                    };
                    self.bvals[c][p] = v;
                }
            }
        }
    }

    /// Writes fixed values into the velocity arrays and copies periodic aliases.
    fn apply_fixed(&mut self) {
        for c in 0..self.cfg.dim {
            let e = Self::ext(&self.cfg, c);
            for p in 0..self.u[c].len() {
                let i = Self::fijk(e, p);
                if self.periodic[c] && i[c] == self.cfg.n[c] {
                    let mut j = i;
                    j[c] = 0;
                    self.u[c][p] = self.u[c][Self::fidx(e, j)];
                } else if self.status[c][p] == Status::Fixed {
                    self.u[c][p] = self.bvals[c][p];
                }
            }
        }
    }

    fn build_helmholtz(&mut self) -> Result<(), String> {
        let dim = self.cfg.dim;
        let k = 0.5 * self.cfg.nu * self.dt;
        let mut out = Vec::new();
        for c in 0..dim {
            let e = Self::ext(&self.cfg, c);
            let total: usize = e.iter().product();
            let mut row_of: Vec<Option<usize>> = vec![None; total];
            let mut rows = Vec::new();
            for p in 0..total {
                if self.status[c][p] == Status::Unknown {
                    row_of[p] = Some(rows.len());
                    rows.push(p);
                }
            }
            let mut trip: Vec<(usize, usize, f64)> = Vec::new();
            let mut terms: Vec<(usize, f64, Term)> = Vec::new();
            for (r, &p) in rows.iter().enumerate() {
                let i = Self::fijk(e, p);
                trip.push((r, r, 1.0));
                for a in 0..dim {
                    let hh = self.h(a);
                    let (thm, thp) = self.ib_theta(c, p, a);
                    trip.push((r, r, 2.0 * k / (hh * hh * thm * thp)));
                    for (si, side) in [-1_isize, 1].into_iter().enumerate() {
                        let th_s = if si == 0 { thm } else { thp };
                        let w = 2.0 * k / (hh * hh * th_s * (thm + thp));
                        if let Some(wv) = self.ib_cross(c, p, a, si) {
                            terms.push((r, w, Term::Const(wv)));
                            continue;
                        }
                        let ni = i[a] as isize + side;
                        let (lo, hi): (isize, isize) = if a == c {
                            (0, self.cfg.n[a] as isize)
                        } else {
                            (0, self.cfg.n[a] as isize - 1)
                        };
                        let mut q = i;
                        if ni >= lo && ni <= hi {
                            q[a] = ni as usize;
                            let qp = Self::fidx(e, q);
                            let target = if self.periodic[a] && a == c && q[a] == self.cfg.n[a] {
                                let mut z = q;
                                z[a] = 0;
                                Self::fidx(e, z)
                            } else {
                                qp
                            };
                            if let Some(rr) = row_of[target] {
                                trip.push((r, rr, -w));
                            } else {
                                terms.push((r, w, Term::Fixed(target)));
                            }
                            continue;
                        }
                        // beyond the array
                        let n = self.cfg.n[a] as isize;
                        if self.periodic[a] {
                            let wrapped = if a == c {
                                if ni < lo { n - 1 } else { 1 }
                            } else if ni < lo {
                                n - 1
                            } else {
                                0
                            };
                            q[a] = wrapped as usize;
                            let qp = Self::fidx(e, q);
                            if let Some(rr) = row_of[qp] {
                                trip.push((r, rr, -w));
                            } else {
                                terms.push((r, w, Term::Fixed(qp)));
                            }
                        } else if a == c {
                            // outflow face: mirror about the boundary face
                            let mi = if ni < lo { 1 } else { n - 1 };
                            q[a] = mi as usize;
                            let qp = Self::fidx(e, q);
                            if let Some(rr) = row_of[qp] {
                                trip.push((r, rr, -w));
                            } else {
                                terms.push((r, w, Term::Fixed(qp)));
                            }
                        } else {
                            let max = ni > hi;
                            let face = Face::of_axis(a, max);
                            match self.cfg.bcs.get(face) {
                                | FluidBoundary::NoSlip(_) | FluidBoundary::Inflow(_) => {
                                    trip.push((r, r, w));
                                    let mut pos = Self::face_pos(&self.cfg, c, i);
                                    pos[a] = if max { self.cfg.len[a] } else { 0.0 };
                                    terms.push((r, 2.0 * w, Term::Wall { face, pos }));
                                },
                                | _ => trip.push((r, r, -w)),
                            }
                        }
                    }
                }
            }
            let a = Csr::from_triplets(rows.len(), rows.len(), &trip);
            out.push(Helm {
                sys: LinearSystem::new(a, LinearSolver::Auto)?,
                rows,
                terms,
            });
        }
        self.helm = out;
        Ok(())
    }

    fn build_poisson(&mut self) -> Result<(), String> {
        let cfg = &self.cfg;
        let n = cfg.n;
        let ncell: usize = n.iter().product();
        let mut trip: Vec<(usize, usize, f64)> = Vec::new();
        let mut has_dirichlet = false;
        for p in 0..ncell {
            if !self.pfluid[p] {
                trip.push((p, p, 1.0));
                continue;
            }
            let i = Self::fijk(n, p);
            for a in 0..cfg.dim {
                let w = 1.0 / (self.h(a) * self.h(a));
                for side in [-1_isize, 1] {
                    let ni = i[a] as isize + side;
                    let mut q = i;
                    let w = w * self.aperture(a, i, side);
                    if w <= 1e-12 {
                        continue;
                    }
                    if ni >= 0 && ni < n[a] as isize {
                        q[a] = ni as usize;
                        let qp = Self::fidx(n, q);
                        if self.pfluid[qp] {
                            trip.push((p, qp, w));
                            trip.push((p, p, -w));
                        }
                    } else if self.periodic[a] {
                        q[a] = if ni < 0 { n[a] - 1 } else { 0 };
                        let qp = Self::fidx(n, q);
                        if self.pfluid[qp] {
                            trip.push((p, qp, w));
                            trip.push((p, p, -w));
                        }
                    } else if matches!(
                        cfg.bcs.get(Face::of_axis(a, ni > 0)),
                        FluidBoundary::Outflow(_)
                    ) {
                        trip.push((p, p, -2.0 * w));
                        has_dirichlet = true;
                    }
                }
            }
        }
        if !has_dirichlet {
            // pin the first fluid cell
            if let Some(p0) = (0..ncell).find(|&p| self.pfluid[p]) {
                trip.retain(|&(r, _, _)| r != p0);
                trip.push((p0, p0, 1.0));
            }
        }
        self.pmat = LinearSystem::new(Csr::from_triplets(ncell, ncell, &trip), cfg.solver)?;
        Ok(())
    }

    /// Extended velocity component (one ghost layer on active axes) at time `t`.
    fn extended(
        &self,
        c: usize,
        t: f64,
    ) -> (Vec<f64>, [usize; 3]) {
        let cfg = &self.cfg;
        let e = Self::ext(cfg, c);
        let mut d = e;
        let mut off = [0usize; 3];
        for a in 0..cfg.dim {
            d[a] += 2;
            off[a] = 1;
        }
        let eid = |i: [usize; 3]| i[0] + d[0] * (i[1] + d[1] * i[2]);
        let mut out = vec![0.0; d.iter().product()];
        for p in 0..self.u[c].len() {
            let i = Self::fijk(e, p);
            out[eid([i[0] + off[0], i[1] + off[1], i[2] + off[2]])] = self.u[c][p];
        }
        for a in 0..cfg.dim {
            // iterate over positions along other axes within range
            let mut lim = e;
            lim[a] = 1;
            let total: usize = lim.iter().product();
            for q in 0..total {
                let mut i = Self::fijk(lim, q);
                for side in [false, true] {
                    // ghost position and the source positions
                    let (g_idx, in_idx, mir_idx, wrap_idx): (usize, usize, usize, usize);
                    let na = e[a];
                    if side {
                        g_idx = na + 1;
                        in_idx = na - 1;
                        mir_idx = na.saturating_sub(2);
                        wrap_idx = usize::from(a == c);
                    } else {
                        g_idx = 0;
                        in_idx = 0;
                        mir_idx = 1.min(na - 1);
                        wrap_idx = if a == c { na - 2 } else { na - 1 };
                    }
                    let mut src = i;

                    let value = if self.periodic[a] {
                        src[a] = wrap_idx;
                        self.u[c][Self::fidx(e, src)]
                    } else if a == c {
                        src[a] = mir_idx;
                        self.u[c][Self::fidx(e, src)]
                    } else {
                        src[a] = in_idx;
                        let uin = self.u[c][Self::fidx(e, src)];
                        let face = Face::of_axis(a, side);
                        match cfg.bcs.get(face) {
                            | FluidBoundary::NoSlip(v) | FluidBoundary::Inflow(v) => {
                                let mut pos = Self::face_pos(cfg, c, src);
                                pos[a] = if side { cfg.len[a] } else { 0.0 };
                                2.0 * v[c].eval(t, pos) - uin
                            },
                            | _ => uin,
                        }
                    };
                    i[a] = g_idx;
                    out[eid([
                        i[0] + if a == 0 { 0 } else { off[0] },
                        i[1] + if a == 1 { 0 } else { off[1] },
                        i[2] + if a == 2 { 0 } else { off[2] },
                    ])] = value;
                }
            }
        }
        (out, d)
    }

    /// Value of an extended component array at (possibly ghost) face indices.
    fn ext_get(
        &self,
        ext: &[(Vec<f64>, [usize; 3])],
        c: usize,
        i: [isize; 3],
    ) -> f64 {
        let (arr, d) = &ext[c];
        let mut ix = [0usize; 3];
        for a in 0..3 {
            ix[a] = (i[a] + isize::from(a < self.cfg.dim)) as usize;
        }
        arr[ix[0] + d[0] * (ix[1] + d[1] * ix[2])]
    }

    /// Like [`Self::ext_get`], `None` outside the one-cell ghost layer.
    fn ext_try(
        &self,
        ext: &[(Vec<f64>, [usize; 3])],
        c: usize,
        i: [isize; 3],
    ) -> Option<f64> {
        let e = Self::ext(&self.cfg, c);
        for a in 0..self.cfg.dim {
            if i[a] < -1 || i[a] > e[a] as isize {
                return None;
            }
        }
        Some(self.ext_get(ext, c, i))
    }

    /// Whether the face at `i` (array index) is inside or touches the body.
    fn is_rough(
        &self,
        c: usize,
        i: [isize; 3],
    ) -> bool {
        let Some(ib) = &self.ib else { return false };
        let e = Self::ext(&self.cfg, c);
        for a in 0..self.cfg.dim {
            if i[a] < 0 || i[a] >= e[a] as isize {
                return false;
            }
        }
        ib.rough[c][Self::fidx(e, [i[0] as usize, i[1] as usize, i[2] as usize])]
    }

    /// Value of the transported quantity on the dual-cell face between the
    /// nodes `nodes[1]` and `nodes[2]` (adjacent along an axis), with the
    /// upwind-biased neighbours `nodes[0]` (before) and `nodes[3]` (after);
    /// `v` is the transporting velocity.
    fn face_value(
        &self,
        ext: &[(Vec<f64>, [usize; 3])],
        c: usize,
        nodes: [[isize; 3]; 4],
        v: f64,
    ) -> f64 {
        let q0 = self.ext_get(ext, c, nodes[1]);
        let q1 = self.ext_get(ext, c, nodes[2]);
        let central = f64::midpoint(q0, q1);
        if self.cfg.convection == Convection::Central {
            return central;
        }
        // the stencil must lie in the ghosted array and away from walls of bodies
        let (up, upup, down) = if v >= 0.0 {
            (nodes[1], nodes[0], nodes[2])
        } else {
            (nodes[2], nodes[3], nodes[1])
        };
        let (Some(qu), Some(quu)) = (self.ext_try(ext, c, up), self.ext_try(ext, c, upup)) else {
            return central;
        };
        if [up, upup, down].iter().any(|&n| self.is_rough(c, n)) {
            return central;
        }
        let qd = self.ext_get(ext, c, down);
        match self.cfg.convection {
            | Convection::Central => central,
            | Convection::Quick => 0.75 * qu + 0.375 * qd - 0.125 * quu,
            | Convection::Tvd(lim) => {
                let den = qd - qu;
                if den.abs() < 1e-300 {
                    return qu;
                }
                let r = (qu - quu) / den;
                qu + 0.5 * lim.psi(r) * den
            },
        }
    }

    /// Transporting velocity `u_a` at the dual face (normal to axis `a`) on side `sgn` of face `ii` of component `c`.
    fn transport_velocity(
        &self,
        ext: &[(Vec<f64>, [usize; 3])],
        c: usize,
        a: usize,
        ii: [isize; 3],
        sgn: isize,
    ) -> f64 {
        let mut hi = ii;
        let mut lo = ii;
        if sgn > 0 {
            hi[a] += 1;
            lo[a] += 1;
        }
        lo[c] -= 1;
        f64::midpoint(self.ext_get(ext, a, lo), self.ext_get(ext, a, hi))
    }

    /// Convective term `∂_a (u_a u_c)` of face `p` of component `c`.
    ///
    /// With `cons` the flux form uses the dual-cell apertures and volume
    /// fraction (and the surface term of a moving rigid body), otherwise the
    /// plain form with a clamped linear ghost extrapolation at the wall.
    fn conv_dof(
        &self,
        ext: &[(Vec<f64>, [usize; 3])],
        c: usize,
        p: usize,
        cons: bool,
    ) -> f64 {
        let cfg = &self.cfg;
        let dim = cfg.dim;
        let e = Self::ext(cfg, c);
        let i = Self::fijk(e, p);
        let ii = [i[0] as isize, i[1] as isize, i[2] as isize];
        let uc = self.ext_get(ext, c, ii);
        let shift = |base: [isize; 3], a: usize, s: isize| {
            let mut b = base;
            b[a] += s;
            b
        };
        let mut total = 0.0;
        let mut sdiv = [0.0; 3];
        for a in 0..dim {
            let h = self.h(a);
            let mut fl = [0.0; 2];
            for (si, sgn) in [(0usize, -1_isize), (1, 1)] {
                let th = self.ib_theta(c, p, a);
                let th = if si == 0 { th.0 } else { th.1 };
                let wall = self.ib_cross(c, p, a, si);
                let (v, q) = wall.map_or_else(
                    || {
                        // nodes along a around the dual face: (before, n0, n1, after)
                        let nodes = if sgn > 0 {
                            [shift(ii, a, -1), ii, shift(ii, a, 1), shift(ii, a, 2)]
                        } else {
                            [shift(ii, a, -2), shift(ii, a, -1), ii, shift(ii, a, 1)]
                        };
                        let v = if a == c {
                            let nb = self.ext_get(ext, c, shift(ii, a, sgn));
                            f64::midpoint(uc, nb)
                        } else {
                            self.transport_velocity(ext, c, a, ii, sgn)
                        };
                        (v, self.face_value(ext, c, nodes, v))
                    },
                    |w| {
                        let val = if cons {
                            uc + (w - uc) * (0.5 / th).min(1.0)
                        } else {
                            f64::midpoint(uc, uc + (w - uc) / th.max(0.3))
                        };
                        (
                            if a == c {
                                val
                            } else {
                                self.transport_velocity(ext, c, a, ii, sgn)
                            },
                            val,
                        )
                    },
                );
                let ap = if cons {
                    self.ib
                        .as_ref()
                        .map_or(1.0, |ib| ib.dual_ap[c][6 * p + 2 * a + si])
                } else {
                    1.0
                };
                fl[si] = ap * v * q;
                sdiv[a] += if si == 1 { ap } else { -ap };
            }
            total += (fl[1] - fl[0]) / h;
        }
        if cons && let Some(ib) = &self.ib {
            let beta = ib.dual_beta[c][p].max(1e-12);
            if ib.moving
                && let Some(body) = &self.cfg.body
            {
                // surface momentum flux u_c (u.n) dA = w_c (w . S) with S = -sum (ap+ - ap-) A_a e_a
                let w = (body.velocity)(Self::face_pos(cfg, c, i));
                let ws: f64 = (0..dim).map(|a| -w[a] * sdiv[a] / self.h(a)).sum();
                total += w[c] * ws;
            }
            total /= beta;
        }
        total
    }

    /// Computes convection and diffusion (explicit part) at all unknown faces.
    fn convection_and_laplacian(
        &self,
        ext: &[(Vec<f64>, [usize; 3])],
    ) -> ([Vec<f64>; 3], [Vec<f64>; 3]) {
        let cfg = &self.cfg;
        let dim = cfg.dim;
        let mut conv: [Vec<f64>; 3] = std::array::from_fn(|_| Vec::new());
        let mut lap: [Vec<f64>; 3] = std::array::from_fn(|_| Vec::new());
        let cons = self.ib.is_some() && cfg.cut_cell.conservative;
        for c in 0..dim {
            let e = Self::ext(cfg, c);
            let total: usize = e.iter().product();
            conv[c] = vec![0.0; total];
            lap[c] = vec![0.0; total];
            let mut redist: Vec<(usize, f64)> = Vec::new();
            for p in 0..total {
                if self.status[c][p] != Status::Unknown {
                    continue;
                }
                let i = Self::fijk(e, p);
                let ii = [i[0] as isize, i[1] as isize, i[2] as isize];
                let uc = self.ext_get(ext, c, ii);
                let mut lp = 0.0;
                for a in 0..dim {
                    let h = self.h(a);
                    let mut ep = ii;
                    ep[a] += 1;
                    let mut em = ii;
                    em[a] -= 1;
                    let (thm, thp) = self.ib_theta(c, p, a);
                    let vp_raw = self.ext_get(ext, c, ep);
                    let vm_raw = self.ext_get(ext, c, em);
                    let (vm, vp) = (
                        self.ib_cross(c, p, a, 0).unwrap_or(vm_raw),
                        self.ib_cross(c, p, a, 1).unwrap_or(vp_raw),
                    );
                    lp += 2.0 / (h * h)
                        * (vm / (thm * (thm + thp)) + vp / (thp * (thm + thp)) - uc / (thm * thp));
                }
                lap[c][p] = lp;
                let cv = self.conv_dof(ext, c, p, cons);
                let beta = self.ib.as_ref().map_or(1.0, |ib| ib.dual_beta[c][p]);
                if cons && beta < cfg.cut_cell.beta_min {
                    // flux redistribution: take the stable non-conservative update, spread the defect
                    let nc = self.conv_dof(ext, c, p, false);
                    conv[c][p] = nc;
                    let hv: f64 = (0..dim).map(|a| self.h(a)).product();
                    let defect = beta * hv * (cv - nc);
                    let mut nbrs: Vec<(usize, f64)> = Vec::new();
                    for a in 0..dim {
                        for s in [-1_isize, 1] {
                            let ni = i[a] as isize + s;
                            if ni < 0 || ni >= e[a] as isize {
                                continue;
                            }
                            let mut q = i;
                            q[a] = ni as usize;
                            let qp = Self::fidx(e, q);
                            if self.status[c][qp] == Status::Unknown {
                                let vol =
                                    self.ib.as_ref().map_or(1.0, |ib| ib.dual_beta[c][qp]) * hv;
                                if vol > 0.0 {
                                    nbrs.push((qp, vol));
                                }
                            }
                        }
                    }
                    let vsum: f64 = nbrs.iter().map(|n| n.1).sum();
                    if vsum > 0.0 {
                        for (qp, _) in nbrs {
                            redist.push((qp, defect / vsum));
                        }
                    } else {
                        conv[c][p] = cv;
                    }
                } else {
                    conv[c][p] = cv;
                }
            }
            for (qp, add) in redist {
                conv[c][qp] += add;
            }
        }
        (conv, lap)
    }

    /// Pressure gradient component at an unknown face.
    fn pressure_gradient(
        &self,
        c: usize,
        i: [usize; 3],
        p: &[f64],
        bc_pressure: bool,
    ) -> f64 {
        let cfg = &self.cfg;
        let n = cfg.n;
        let h = self.h(c);
        let cell = |ix: [usize; 3]| p[Self::fidx(n, ix)];
        if i[c] == 0 {
            if self.periodic[c] {
                let mut q = i;
                q[c] = n[c] - 1;
                return (cell(i) - cell(q)) / h;
            }
            let pout = match cfg.bcs.get(Face::of_axis(c, false)) {
                | FluidBoundary::Outflow(v) if bc_pressure => {
                    v.eval(self.t, Self::face_pos(cfg, c, i))
                },
                | FluidBoundary::Outflow(_) => 0.0,
                | _ => cell(i),
            };
            return (cell(i) - pout) / (0.5 * h);
        }
        if i[c] == n[c] {
            let mut q = i;
            q[c] = n[c] - 1;
            let pout = match cfg.bcs.get(Face::of_axis(c, true)) {
                | FluidBoundary::Outflow(v) if bc_pressure => {
                    v.eval(self.t, Self::face_pos(cfg, c, i))
                },
                | FluidBoundary::Outflow(_) => 0.0,
                | _ => cell(q),
            };
            return (pout - cell(q)) / (0.5 * h);
        }
        let mut q = i;
        q[c] -= 1;
        (cell(i) - cell(q)) / h
    }

    /// Divergence at every cell.
    #[must_use]
    pub fn divergence(&self) -> Vec<f64> {
        let cfg = &self.cfg;
        let ncell: usize = cfg.n.iter().product();
        let mut d = vec![0.0; ncell];
        for p in 0..ncell {
            let i = Self::fijk(cfg.n, p);
            let mut s = 0.0;
            for c in 0..cfg.dim {
                let e = Self::ext(cfg, c);
                let mut j = i;
                j[c] += 1;
                let (ap_hi, ap_lo) = self.ib.as_ref().map_or((1.0, 1.0), |ib| {
                    (
                        ib.aperture[c][Self::fidx(e, j)],
                        ib.aperture[c][Self::fidx(e, i)],
                    )
                });
                s += (ap_hi * self.u[c][Self::fidx(e, j)] - ap_lo * self.u[c][Self::fidx(e, i)])
                    / self.h(c);
            }
            d[p] = s;
        }
        d
    }

    /// Advances one time step.
    ///
    /// # Errors
    /// If a linear solve fails.
    pub fn step(&mut self) -> Result<(), String> {
        let cfg_dim = self.cfg.dim;
        let dt = self.dt;
        let t_new = self.t + dt;
        let nu = self.cfg.nu;
        let ext: Vec<(Vec<f64>, [usize; 3])> =
            (0..cfg_dim).map(|c| self.extended(c, self.t)).collect();
        let (conv, lap) = self.convection_and_laplacian(&ext);
        let les = self.les_force(&ext);
        if self.th.is_some() {
            let adv = self
                .temperature_advection(&self.th.as_ref().map(|t| t.t.clone()).unwrap_or_default());
            if let Some(st) = &mut self.th {
                st.adv_now = adv;
            }
        }
        let adv: [Vec<f64>; 3] = std::array::from_fn(|c| {
            if c >= cfg_dim {
                return Vec::new();
            }
            self.conv_prev.as_ref().map_or_else(
                || conv[c].clone(),
                |prev| {
                    conv[c]
                        .iter()
                        .zip(&prev[c])
                        .map(|(a, b)| 1.5 * a - 0.5 * b)
                        .collect()
                },
            )
        });
        self.update_boundary_values(t_new);
        let mut ustar: [Vec<f64>; 3] = std::array::from_fn(|c| self.u[c].clone());
        for c in 0..cfg_dim {
            let e = Self::ext(&self.cfg, c);
            let helm = &self.helm[c];
            let mut rhs = vec![0.0; helm.rows.len()];
            for (r, &p) in helm.rows.iter().enumerate() {
                let i = Self::fijk(e, p);
                let x = Self::face_pos(&self.cfg, c, i);
                let f = self
                    .cfg
                    .force
                    .as_ref()
                    .map_or(0.0, |ff| ff(self.t + 0.5 * dt, x)[c])
                    + self.buoyancy(c, i)
                    + les.as_ref().map_or(0.0, |l| l[c][p]);
                let g = self.pressure_gradient(c, i, &self.p, true);
                rhs[r] = self.u[c][p] + dt * (-adv[c][p] - g + f) + 0.5 * nu * dt * lap[c][p];
            }
            for (r, coef, term) in &helm.terms {
                rhs[*r] += coef
                    * match term {
                        | Term::Const(v) => *v,
                        | Term::Fixed(idx) => self.bvals[c][*idx],
                        | Term::Wall { face, pos } => {
                            match self.cfg.bcs.get(*face) {
                                | FluidBoundary::NoSlip(v) | FluidBoundary::Inflow(v) => {
                                    v[c].eval(t_new, *pos)
                                },
                                | _ => 0.0,
                            }
                        },
                    };
            }
            let sol = helm.sys.solve(&rhs, None)?;
            for (r, &p) in helm.rows.iter().enumerate() {
                ustar[c][p] = sol[r];
            }
        }
        let old = std::mem::replace(&mut self.u, ustar);
        self.apply_fixed();
        // projection
        let div = self.divergence();
        let ncell = div.len();
        let rhs: Vec<f64> = (0..ncell)
            .map(|p| {
                if self.pfluid[p] {
                    div[p] / dt
                } else {
                    0.0
                }
            })
            .collect();
        let mut rhs = rhs;
        if !(0..self.cfg.dim).any(|a| {
            !self.periodic[a]
                && (matches!(
                    self.cfg.bcs.get(Face::of_axis(a, false)),
                    FluidBoundary::Outflow(_)
                ) || matches!(
                    self.cfg.bcs.get(Face::of_axis(a, true)),
                    FluidBoundary::Outflow(_)
                ))
        }) && let Some(p0) = (0..ncell).find(|&p| self.pfluid[p])
        {
            rhs[p0] = 0.0;
        }
        let phi = self.pmat.solve(&rhs, None)?;
        for c in 0..self.cfg.dim {
            let e = Self::ext(&self.cfg, c);
            for p in 0..self.u[c].len() {
                if self.status[c][p] == Status::Unknown {
                    let i = Self::fijk(e, p);
                    self.u[c][p] -= dt * self.pressure_gradient(c, i, &phi, false);
                }
            }
        }
        // periodic aliases
        self.apply_fixed();
        for p in 0..ncell {
            if self.pfluid[p] {
                self.p[p] += phi[p]
                    - if self.cfg.rotational {
                        0.5 * nu * div[p]
                    } else {
                        0.0
                    };
            }
        }
        let _ = old;
        self.advance_temperature()?;
        self.conv_prev = Some(std::array::from_fn(|c| conv[c].clone()));
        self.t = t_new;
        self.steps += 1;
        Ok(())
    }

    /// Runs `steps` steps.
    ///
    /// # Errors
    /// See [`NsSolver::step`].
    pub fn run(
        &mut self,
        steps: usize,
    ) -> Result<(), String> {
        for _ in 0..steps {
            self.step()?;
        }
        Ok(())
    }

    /// Number of steps taken.
    #[must_use]
    pub const fn steps(&self) -> usize {
        self.steps
    }

    /// Cells per axis.
    #[must_use]
    pub const fn n(&self) -> [usize; 3] {
        self.cfg.n
    }

    /// Cell size.
    #[must_use]
    pub fn cell_size(
        &self,
        a: usize,
    ) -> f64 {
        self.h(a)
    }

    /// Multilinear interpolation of velocity component `c` at `x`.
    #[must_use]
    pub fn sample_velocity(
        &self,
        c: usize,
        x: [f64; 3],
    ) -> f64 {
        let cfg = &self.cfg;
        let e = Self::ext(cfg, c);
        let mut lo = [0usize; 3];
        let mut fr = [0.0; 3];
        for a in 0..cfg.dim {
            let h = self.h(a);
            let s = x[a] / h - if a == c { 0.0 } else { 0.5 };
            let s = s.clamp(0.0, (e[a] - 1) as f64);
            let l = (s.floor() as usize).min(e[a].saturating_sub(2));
            lo[a] = l;
            fr[a] = s - l as f64;
        }
        let mut v = 0.0;
        let corners = 1 << cfg.dim;
        for m in 0..corners {
            let mut w = 1.0;
            let mut ix = lo;
            for a in 0..cfg.dim {
                if m >> a & 1 == 1 {
                    w *= fr[a];
                    ix[a] += 1;
                } else {
                    w *= 1.0 - fr[a];
                }
            }
            v += w * self.u[c][Self::fidx(e, ix)];
        }
        v
    }

    /// Cell-centred pressure at flat cell index.
    #[must_use]
    pub const fn cell_index(
        &self,
        i: [usize; 3],
    ) -> usize {
        Self::fidx(self.cfg.n, i)
    }

    /// Kinetic energy `½ ∫ |u|²` using face values.
    #[must_use]
    pub fn kinetic_energy(&self) -> f64 {
        let cfg = &self.cfg;
        let vol: f64 = (0..cfg.dim).map(|a| self.h(a)).product();
        let mut s = 0.0;
        for c in 0..cfg.dim {
            let e = Self::ext(cfg, c);
            for p in 0..self.u[c].len() {
                let i = Self::fijk(e, p);
                let wgt = if i[c] == 0 || i[c] == cfg.n[c] {
                    0.5
                } else {
                    1.0
                };
                s += wgt * self.u[c][p] * self.u[c][p];
            }
        }
        0.5 * s * vol
    }

    /// Volumetric flux through the plane `x_axis = index` (face index).
    #[must_use]
    pub fn flux_through(
        &self,
        axis: usize,
        face_index: usize,
    ) -> f64 {
        let cfg = &self.cfg;
        let e = Self::ext(cfg, axis);
        let mut s = 0.0;
        for p in 0..self.u[axis].len() {
            let i = Self::fijk(e, p);
            if i[axis] == face_index {
                let area: f64 = (0..cfg.dim)
                    .filter(|&a| a != axis)
                    .map(|a| self.h(a))
                    .product();
                s += self.u[axis][p] * area;
            }
        }
        s
    }

    /// Maximum speed over the faces (for CFL estimates).
    #[must_use]
    pub fn max_speed(&self) -> f64 {
        self.u
            .iter()
            .flat_map(|v| v.iter())
            .fold(0.0_f64, |m, x| m.max(x.abs()))
    }
}

/// Temperature field state of the Boussinesq model.
struct ThermalState {
    t: Vec<f64>,
    t_prev: Option<Vec<f64>>,
    adv_prev: Option<Vec<f64>>,
    adv_now: Vec<f64>,
    sys: LinearSystem,
    amat: Csr,
    bsrc: Vec<f64>,
}

impl NsSolver {
    fn build_thermal(
        cfg: &NsConfig,
        dt: f64,
    ) -> Result<Option<ThermalState>, String> {
        let Some(th) = &cfg.thermal else {
            return Ok(None);
        };
        if cfg.body.is_some() || cfg.solid.is_some() {
            return Err("the thermal field needs a domain without obstacles".into());
        }
        if th.kappa < 0.0 {
            return Err("thermal diffusivity must not be negative".into());
        }
        for a in 0..cfg.dim {
            if (th.bcs[2 * a] == ScalarBc::Periodic) != (th.bcs[2 * a + 1] == ScalarBc::Periodic) {
                return Err("periodic temperature must be set on both faces of an axis".into());
            }
        }
        let n = cfg.n;
        let ncell: usize = n.iter().product();
        let mut atrip: Vec<(usize, usize, f64)> = Vec::new();
        let mut mtrip: Vec<(usize, usize, f64)> = Vec::new();
        let mut bsrc = vec![0.0; ncell];
        let mut t = vec![0.0; ncell];
        for p in 0..ncell {
            let i = Self::fijk(n, p);
            let mut x = [0.0; 3];
            for a in 0..cfg.dim {
                x[a] = (i[a] as f64 + 0.5) * cfg.len[a] / n[a] as f64;
            }
            t[p] = (th.init)(x);
            mtrip.push((p, p, 1.0));
            for a in 0..cfg.dim {
                let h = cfg.len[a] / n[a] as f64;
                let w = 1.0 / (h * h);
                for (si, s) in [-1_isize, 1].into_iter().enumerate() {
                    let k = i[a] as isize + s;
                    let mut q = i;
                    let inside = k >= 0 && k < n[a] as isize;
                    let periodic = th.bcs[2 * a] == ScalarBc::Periodic;
                    if inside || periodic {
                        q[a] = k.rem_euclid(n[a] as isize) as usize;
                        let qp = Self::fidx(n, q);
                        atrip.push((p, qp, w));
                        atrip.push((p, p, -w));
                        mtrip.push((p, qp, -0.5 * th.kappa * dt * w));
                        mtrip.push((p, p, 0.5 * th.kappa * dt * w));
                    } else {
                        match th.bcs[2 * a + si] {
                            | ScalarBc::Adiabatic | ScalarBc::Periodic => {},
                            | ScalarBc::Dirichlet(tw) => {
                                atrip.push((p, p, -2.0 * w));
                                bsrc[p] += 2.0 * w * tw;
                                mtrip.push((p, p, 0.5 * th.kappa * dt * 2.0 * w));
                            },
                        }
                    }
                }
            }
        }
        let amat = Csr::from_triplets(ncell, ncell, &atrip);
        let sys = LinearSystem::new(Csr::from_triplets(ncell, ncell, &mtrip), LinearSolver::Auto)?;
        Ok(Some(ThermalState {
            t,
            t_prev: None,
            adv_prev: None,
            adv_now: vec![0.0; ncell],
            sys,
            amat,
            bsrc,
        }))
    }

    /// Face value of the temperature between cells `f - 1` and `f` along `a`
    /// (`f` in `0..=n`), for the transporting velocity `v`.
    fn temperature_face(
        &self,
        t: &[f64],
        i: [usize; 3],
        a: usize,
        f: usize,
        v: f64,
    ) -> f64 {
        let n = self.cfg.n;
        let th = self.cfg.thermal.as_ref();
        let cell = |k: isize| -> Option<f64> {
            let na = n[a] as isize;
            let mut q = i;
            if k >= 0 && k < na {
                q[a] = k as usize;
                Some(t[Self::fidx(n, q)])
            } else if th.is_some_and(|t| t.bcs[2 * a] == ScalarBc::Periodic) {
                q[a] = k.rem_euclid(na) as usize;
                Some(t[Self::fidx(n, q)])
            } else {
                None
            }
        };
        let fi = f as isize;
        let (c0, c1) = (cell(fi - 1), cell(fi));
        match (c0, c1) {
            | (Some(l), Some(r)) => {
                let central = f64::midpoint(l, r);
                if self.cfg.convection == Convection::Central {
                    return central;
                }
                let (up, upup, down) = if v >= 0.0 {
                    (l, cell(fi - 2), r)
                } else {
                    (r, cell(fi + 1), l)
                };
                let Some(uu) = upup else { return central };
                match self.cfg.convection {
                    | Convection::Central => central,
                    | Convection::Quick => 0.75 * up + 0.375 * down - 0.125 * uu,
                    | Convection::Tvd(lim) => {
                        let den = down - up;
                        if den.abs() < 1e-300 {
                            up
                        } else {
                            up + 0.5 * lim.psi((up - uu) / den) * den
                        }
                    },
                }
            },
            | (None, Some(r)) => {
                match th.map(|t| t.bcs[2 * a]) {
                    | Some(ScalarBc::Dirichlet(tw)) => tw,
                    | _ => r,
                }
            },
            | (Some(l), None) => {
                match th.map(|t| t.bcs[2 * a + 1]) {
                    | Some(ScalarBc::Dirichlet(tw)) => tw,
                    | _ => l,
                }
            },
            | (None, None) => 0.0,
        }
    }

    /// Conservative advection `∇·(u T)` of the temperature at every cell.
    fn temperature_advection(
        &self,
        t: &[f64],
    ) -> Vec<f64> {
        let cfg = &self.cfg;
        let n = cfg.n;
        let ncell: usize = n.iter().product();
        let mut out = vec![0.0; ncell];
        for p in 0..ncell {
            let i = Self::fijk(n, p);
            let mut s = 0.0;
            for a in 0..cfg.dim {
                let e = Self::ext(cfg, a);
                let mut fl = [0.0; 2];
                for (k, f) in [i[a], i[a] + 1].into_iter().enumerate() {
                    let mut fi = i;
                    fi[a] = f;
                    let v = self.u[a][Self::fidx(e, fi)];
                    fl[k] = v * self.temperature_face(t, i, a, f, v);
                }
                s += (fl[1] - fl[0]) / self.h(a);
            }
            out[p] = s;
        }
        out
    }

    /// Cell-centred temperature (`None` without a thermal model).
    #[must_use]
    pub fn temperature(&self) -> Option<&[f64]> {
        self.th.as_ref().map(|s| s.t.as_slice())
    }

    /// Overwrites the temperature field.
    ///
    /// # Errors
    /// Without a thermal model or for a wrong length.
    pub fn set_temperature(
        &mut self,
        t: &[f64],
    ) -> Result<(), String> {
        let st = self.th.as_mut().ok_or("no thermal model configured")?;
        if t.len() != st.t.len() {
            return Err("temperature field has the wrong length".into());
        }
        st.t.copy_from_slice(t);
        st.t_prev = None;
        st.adv_prev = None;
        Ok(())
    }

    /// Advances the temperature by one step with the velocity at the old time.
    fn advance_temperature(&mut self) -> Result<(), String> {
        let dt = self.dt;
        let Some(st) = &self.th else { return Ok(()) };
        let kappa = self.cfg.thermal.as_ref().map_or(0.0, |t| t.kappa);
        let adv = st.adv_now.clone();
        let adv_ab = st.adv_prev.as_ref().map_or_else(
            || adv.clone(),
            |prev| {
                adv.iter()
                    .zip(prev)
                    .map(|(a, b)| 1.5 * a - 0.5 * b)
                    .collect()
            },
        );
        let ncell = st.t.len();
        let mut lt = vec![0.0; ncell];
        st.amat.matvec(&st.t, &mut lt);
        let rhs: Vec<f64> = (0..ncell)
            .map(|p| st.t[p] - dt * adv_ab[p] + 0.5 * kappa * dt * lt[p] + kappa * dt * st.bsrc[p])
            .collect();
        let new = st.sys.solve(&rhs, Some(&st.t))?;
        if let Some(st) = &mut self.th {
            st.t_prev = Some(std::mem::replace(&mut st.t, new));
            st.adv_prev = Some(adv);
        }
        Ok(())
    }

    /// Boussinesq force on face `p` of component `c` from the extrapolated temperature.
    fn buoyancy(
        &self,
        c: usize,
        i: [usize; 3],
    ) -> f64 {
        let (Some(st), Some(th)) = (&self.th, &self.cfg.thermal) else {
            return 0.0;
        };
        let n = self.cfg.n;
        let tval = |q: [usize; 3]| {
            let k = Self::fidx(n, q);
            st.t_prev
                .as_ref()
                .map_or(st.t[k], |tp| 1.5 * st.t[k] - 0.5 * tp[k])
        };
        let mut hi = i;
        hi[c] = i[c].min(n[c] - 1);
        let mut lo = i;
        lo[c] = i[c].saturating_sub(1).min(n[c] - 1);
        th.gbeta[c] * (f64::midpoint(tval(hi), tval(lo)) - th.t_ref)
    }

    /// Velocity gradient `G[c][a] = ∂_a u_c` at the centre of cell `j`.
    fn centre_gradient(
        &self,
        ext: &[(Vec<f64>, [usize; 3])],
        j: [isize; 3],
    ) -> [[f64; 3]; 3] {
        let dim = self.cfg.dim;
        let mut g = [[0.0; 3]; 3];
        let ucen = |c: usize, cell: [isize; 3]| {
            let mut hi = cell;
            hi[c] += 1;
            f64::midpoint(self.ext_get(ext, c, cell), self.ext_get(ext, c, hi))
        };
        for c in 0..dim {
            for a in 0..dim {
                if a == c {
                    let mut hi = j;
                    hi[c] += 1;
                    g[c][c] = (self.ext_get(ext, c, hi) - self.ext_get(ext, c, j)) / self.h(c);
                } else {
                    let (mut jp, mut jm) = (j, j);
                    jp[a] += 1;
                    jm[a] -= 1;
                    g[c][a] = (ucen(c, jp) - ucen(c, jm)) / (2.0 * self.h(a));
                }
            }
        }
        g
    }

    fn eddy_viscosity_from(
        &self,
        ext: &[(Vec<f64>, [usize; 3])],
    ) -> Vec<f64> {
        let cfg = &self.cfg;
        let ncell: usize = cfg.n.iter().product();
        let delta = (0..cfg.dim)
            .map(|a| self.h(a))
            .product::<f64>()
            .powf(1.0 / cfg.dim as f64);
        let mut out = vec![0.0; ncell];
        for (p, nut) in out.iter_mut().enumerate() {
            let i = Self::fijk(cfg.n, p);
            let g = self.centre_gradient(ext, [i[0] as isize, i[1] as isize, i[2] as isize]);
            let mut sij = [[0.0; 3]; 3];
            for r in 0..3 {
                for c in 0..3 {
                    sij[r][c] = f64::midpoint(g[r][c], g[c][r]);
                }
            }
            let ss: f64 = sij.iter().flatten().map(|x| x * x).sum();
            *nut = match cfg.les {
                | LesModel::None => 0.0,
                | LesModel::Smagorinsky { cs } => (cs * delta).powi(2) * (2.0 * ss).sqrt(),
                | LesModel::Wale { cw } => {
                    let mut g2 = [[0.0; 3]; 3];
                    for r in 0..3 {
                        for c in 0..3 {
                            g2[r][c] = (0..3).map(|k| g[r][k] * g[k][c]).sum();
                        }
                    }
                    let tr = g2[0][0] + g2[1][1] + g2[2][2];
                    let mut sd2 = 0.0;
                    for r in 0..3 {
                        for c in 0..3 {
                            let v = f64::midpoint(g2[r][c], g2[c][r])
                                - if r == c { tr / 3.0 } else { 0.0 };
                            sd2 += v * v;
                        }
                    }
                    let den = ss.powf(2.5) + sd2.powf(1.25);
                    if den < 1e-300 {
                        0.0
                    } else {
                        (cw * delta).powi(2) * sd2.powf(1.5) / den
                    }
                },
            };
        }
        out
    }

    /// Eddy viscosity `ν_t` of the sub-grid-scale model at every cell centre.
    #[must_use]
    pub fn eddy_viscosity(&self) -> Vec<f64> {
        let ext: Vec<(Vec<f64>, [usize; 3])> = (0..self.cfg.dim)
            .map(|c| self.extended(c, self.t))
            .collect();
        self.eddy_viscosity_from(&ext)
    }

    /// Explicit sub-grid force `∇·(ν_t (∇u + ∇uᵀ))` at every face.
    fn les_force(
        &self,
        ext: &[(Vec<f64>, [usize; 3])],
    ) -> Option<[Vec<f64>; 3]> {
        if self.cfg.les == LesModel::None {
            return None;
        }
        let cfg = &self.cfg;
        let dim = cfg.dim;
        let nut = self.eddy_viscosity_from(ext);
        let nu_at = |c: [isize; 3]| {
            let mut q = [0usize; 3];
            for a in 0..3 {
                q[a] = c[a].clamp(0, cfg.n[a] as isize - 1) as usize;
            }
            nut[Self::fidx(cfg.n, q)]
        };
        let mut out: [Vec<f64>; 3] = std::array::from_fn(|_| Vec::new());
        for c in 0..dim {
            let e = Self::ext(cfg, c);
            let total: usize = e.iter().product();
            out[c] = vec![0.0; total];
            for p in 0..total {
                if self.status[c][p] != Status::Unknown {
                    continue;
                }
                let i = Self::fijk(e, p);
                let ii = [i[0] as isize, i[1] as isize, i[2] as isize];
                let uc = self.ext_get(ext, c, ii);
                let mut f = 0.0;
                for a in 0..dim {
                    let h = self.h(a);
                    if a == c {
                        let mut cl = ii;
                        cl[c] -= 1;
                        let tau = |cell: [isize; 3]| {
                            let mut cc = cell;
                            cc[c] = cc[c].clamp(0, cfg.n[c] as isize - 1);
                            let g = self.centre_gradient(ext, cc);
                            2.0 * nu_at(cc) * g[c][c]
                        };
                        f += (tau(ii) - tau(cl)) / h;
                    } else {
                        let mut tau = [0.0; 2];
                        for (si, sgn) in [(0usize, -1_isize), (1, 1)] {
                            let mut nb = ii;
                            let (dcu, dau);
                            if sgn > 0 {
                                nb[a] += 1;
                                dcu = (self.ext_get(ext, c, nb) - uc) / h;
                                let mut hi = ii;
                                hi[a] += 1;
                                let mut lo = hi;
                                lo[c] -= 1;
                                dau = (self.ext_get(ext, a, hi) - self.ext_get(ext, a, lo))
                                    / self.h(c);
                            } else {
                                nb[a] -= 1;
                                dcu = (uc - self.ext_get(ext, c, nb)) / h;
                                let mut lo = ii;
                                lo[c] -= 1;
                                dau = (self.ext_get(ext, a, ii) - self.ext_get(ext, a, lo))
                                    / self.h(c);
                            }
                            // viscosity at the edge: mean of the four cells around it
                            let ia = if sgn > 0 { ii[a] } else { ii[a] - 1 };
                            let mut nu_e = 0.0;
                            for dc in [-1_isize, 0] {
                                for da in [0_isize, 1] {
                                    let mut cell = ii;
                                    cell[c] += dc;
                                    cell[a] = ia + da;
                                    nu_e += 0.25 * nu_at(cell);
                                }
                            }
                            tau[si] = nu_e * (dcu + dau);
                        }
                        f += (tau[1] - tau[0]) / h;
                    }
                }
                out[c][p] = f;
            }
        }
        Some(out)
    }

    /// Multilinear interpolation of the cell-centred pressure at `x`.
    #[must_use]
    pub fn sample_pressure(
        &self,
        x: [f64; 3],
    ) -> f64 {
        let cfg = &self.cfg;
        let mut lo = [0usize; 3];
        let mut fr = [0.0; 3];
        for a in 0..cfg.dim {
            let s = (x[a] / self.h(a) - 0.5).clamp(0.0, (cfg.n[a] - 1) as f64);
            let l = (s.floor() as usize).min(cfg.n[a].saturating_sub(2));
            lo[a] = l;
            fr[a] = s - l as f64;
        }
        let mut v = 0.0;
        for m in 0..(1usize << cfg.dim) {
            let mut w = 1.0;
            let mut ix = lo;
            for a in 0..cfg.dim {
                if m >> a & 1 == 1 {
                    w *= fr[a];
                    ix[a] = (ix[a] + 1).min(cfg.n[a] - 1);
                } else {
                    w *= 1.0 - fr[a];
                }
            }
            v += w * self.p[Self::fidx(cfg.n, ix)];
        }
        v
    }

    /// Fluid volume fraction of pressure cell `i` (1 without an immersed body).
    #[must_use]
    pub fn volume_fraction(
        &self,
        i: [usize; 3],
    ) -> f64 {
        self.ib
            .as_ref()
            .map_or(1.0, |ib| ib.vol_frac[Self::fidx(self.cfg.n, i)])
    }

    /// Open fraction of the cell face normal to axis `axis` with face index `i`.
    #[must_use]
    pub fn face_aperture(
        &self,
        axis: usize,
        i: [usize; 3],
    ) -> f64 {
        self.ib.as_ref().map_or(1.0, |ib| {
            ib.aperture[axis][Self::fidx(Self::ext(&self.cfg, axis), i)]
        })
    }

    /// Force exerted by the fluid on whatever lies inside the box `[lo, hi]`,
    /// from the steady momentum balance on its surface:
    /// `F = -∮ [u (u·n) + p n - ν ∂_n u] dS` (densities are one).
    ///
    /// The box must lie in the fluid and the flow must be steady for the
    /// result to be the drag of the enclosed body. In 2D the force is per unit
    /// length.
    #[must_use]
    pub fn force_on_region(
        &self,
        lo: [f64; 3],
        hi: [f64; 3],
    ) -> [f64; 3] {
        let cfg = &self.cfg;
        let dim = cfg.dim;
        let mut force = [0.0; 3];
        let vel = |x: [f64; 3]| -> [f64; 3] {
            std::array::from_fn(|c| {
                if c < dim {
                    self.sample_velocity(c, x)
                } else {
                    0.0
                }
            })
        };
        for d in 0..dim {
            let others: Vec<usize> = (0..dim).filter(|&b| b != d).collect();
            let counts: Vec<usize> = others
                .iter()
                .map(|&b| (((hi[b] - lo[b]) / self.h(b)).round() as usize).max(1))
                .collect();
            let total: usize = counts.iter().product();
            let da: f64 = others
                .iter()
                .zip(&counts)
                .map(|(&b, &m)| (hi[b] - lo[b]) / m as f64)
                .product();
            for (side, nd) in [(lo[d], -1.0_f64), (hi[d], 1.0)] {
                for q in 0..total {
                    let mut x = [0.0; 3];
                    x[d] = side;
                    let mut rem = q;
                    for (k, &b) in others.iter().enumerate() {
                        let idx = rem % counts[k];
                        rem /= counts[k];
                        x[b] = lo[b] + (idx as f64 + 0.5) * (hi[b] - lo[b]) / counts[k] as f64;
                    }
                    let u = vel(x);
                    let pr = self.sample_pressure(x);
                    let dl = self.h(d);
                    let (mut xp, mut xm) = (x, x);
                    xp[d] += dl;
                    xm[d] -= dl;
                    let (up, um) = (vel(xp), vel(xm));
                    let un = nd * u[d];
                    for k in 0..dim {
                        let dn = nd * (up[k] - um[k]) / (2.0 * dl);
                        let normal = if k == d { nd } else { 0.0 };
                        force[k] -= da * (u[k] * un + pr * normal - cfg.nu * dn);
                    }
                }
            }
        }
        force
    }
}
