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
//!   diffusion-dominated flow, see the Taylor–Couette test; the convective term
//!   uses a clamped linear ghost extrapolation).
//!
//! | Feature | 2D | 3D |
//! |---|---|---|
//! | Navier–Stokes, constant density | yes (`nz = 1`) | yes |
//! | no-slip / moving wall, slip, inflow, outflow, periodic | yes | yes |
//! | obstacles (cell mask) | yes | yes |
//! | body force `f(t, x)` | yes | yes |
//!
//! Limitations: uniform grids, central convection (use a mesh Péclet number
//! `u h / ν` of a few units or less), no turbulence model, no free surface.

use crate::kernels::krylov::Csr;
use crate::sim::boundary::{Face, FluidBoundary, FluidFaceBcs};
use crate::sim::structured::{LinearSolver, LinearSystem};

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
}

impl NsConfig {
    /// Creates a configuration with all walls no-slip.
    ///
    /// # Errors
    /// For bad dimensions or sizes.
    pub fn new(n: &[usize], len: &[f64], nu: f64) -> Result<Self, String> {
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
}

struct Helm {
    sys: LinearSystem,
    rows: Vec<usize>,   // face index of each unknown row
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
}

impl NsSolver {
    const fn ext(cfg: &NsConfig, c: usize) -> [usize; 3] {
        let mut e = cfg.n;
        if c < cfg.dim {
            e[c] += 1;
        }
        e
    }

    const fn fidx(e: [usize; 3], i: [usize; 3]) -> usize {
        i[0] + e[0] * (i[1] + e[1] * i[2])
    }

    const fn fijk(e: [usize; 3], p: usize) -> [usize; 3] {
        [p % e[0], (p / e[0]) % e[1], p / (e[0] * e[1])]
    }

    fn cell_solid(cfg: &NsConfig, i: [isize; 3]) -> bool {
        for a in 0..3 {
            if i[a] < 0 || i[a] >= cfg.n[a] as isize {
                return false;
            }
        }
        cfg.solid.as_ref().is_some_and(|s| s[Self::fidx(cfg.n, [i[0] as usize, i[1] as usize, i[2] as usize])])
    }

    fn face_pos(cfg: &NsConfig, c: usize, i: [usize; 3]) -> [f64; 3] {
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
                // aperture of this face (cell face normal to axis c)
                let others: Vec<usize> = (0..dim).filter(|&b| b != c).collect();
                aperture[p] = if others.len() == 1 {
                    let b = others[0];
                    let (mut x1, mut x2) = (x, x);
                    x1[b] -= 0.5 * hs[b];
                    x2[b] += 0.5 * hs[b];
                    let (f1, f2) = (phi(x1), phi(x2));
                    if f1 > 0.0 && f2 > 0.0 {
                        1.0
                    } else if f1 <= 0.0 && f2 <= 0.0 {
                        0.0
                    } else {
                        let (pos, neg) = if f1 > 0.0 { (f1, f2) } else { (f2, f1) };
                        pos / (pos - neg)
                    }
                } else {
                    let ns = 8;
                    let mut open = 0;
                    for s1 in 0..ns {
                        for s2 in 0..ns {
                            let mut xs = x;
                            xs[others[0]] += ((f64::from(s1) + 0.5) / f64::from(ns) - 0.5) * hs[others[0]];
                            xs[others[1]] += ((f64::from(s2) + 0.5) / f64::from(ns) - 0.5) * hs[others[1]];
                            if phi(xs) > 0.0 {
                                open += 1;
                            }
                        }
                    }
                    f64::from(open) / f64::from(ns * ns)
                };
            }
            ib.inside[c] = inside;
            ib.inside_val[c] = inside_val;
            ib.cross[c] = cross;
            ib.theta[c] = theta;
            ib.wall[c] = wall;
            ib.aperture[c] = aperture;
        }
        Some(ib)
    }

    fn ib_theta(&self, c: usize, p: usize, a: usize) -> (f64, f64) {
        self.ib.as_ref().map_or((1.0, 1.0), |ib| (ib.theta[c][6 * p + 2 * a], ib.theta[c][6 * p + 2 * a + 1]))
    }

    fn ib_cross(&self, c: usize, p: usize, a: usize, s: usize) -> Option<f64> {
        let ib = self.ib.as_ref()?;
        ib.cross[c][6 * p + 2 * a + s].then(|| ib.wall[c][6 * p + 2 * a + s])
    }

    fn aperture(&self, a: usize, i: [usize; 3], side: isize) -> f64 {
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
    pub fn new(cfg: NsConfig, dt: f64, init: &dyn Fn([f64; 3]) -> [f64; 3]) -> Result<Self, String> {
        let dim = cfg.dim;
        let mut periodic = [false; 3];
        for a in 0..dim {
            let lo = matches!(cfg.bcs.get(Face::of_axis(a, false)), FluidBoundary::Periodic);
            let hi = matches!(cfg.bcs.get(Face::of_axis(a, true)), FluidBoundary::Periodic);
            if lo != hi {
                return Err("periodic must be set on both faces of an axis".into());
            }
            periodic[a] = lo;
        }
        if let Some(s) = &cfg.solid
            && s.len() != cfg.n.iter().product::<usize>() {
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
                        FluidBoundary::Outflow(_) => {},
                        FluidBoundary::Periodic => {
                            if max {
                                fixed = true; // alias of face 0
                            }
                        },
                        _ => fixed = true,
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
                st[p] = if fixed { Status::Fixed } else { Status::Unknown };
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
                    b.aperture[a][Self::fidx(e, i)] > 1e-12 || b.aperture[a][Self::fidx(e, j)] > 1e-12
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
            pmat: LinearSystem::new(Csr::from_triplets(1, 1, &[(0, 0, 1.0)]), LinearSolver::Direct)?,
            pfluid,
            ib,
            conv_prev: None,
            periodic,
            steps: 0,
        };
        solver.update_boundary_values(0.0);
        solver.apply_fixed();
        solver.build_helmholtz()?;
        solver.build_poisson()?;
        Ok(solver)
    }

    fn h(&self, a: usize) -> f64 {
        self.cfg.len[a] / self.cfg.n[a] as f64
    }

    /// Boundary value of a fixed face at time `t`.
    fn bc_value(&self, c: usize, i: [usize; 3], t: f64) -> f64 {
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
                FluidBoundary::NoSlip(v) | FluidBoundary::Inflow(v) => v[c].eval(t, x),
                _ => 0.0,
            };
        }
        0.0
    }

    fn update_boundary_values(&mut self, t: f64) {
        for c in 0..self.cfg.dim {
            let e = Self::ext(&self.cfg, c);
            for p in 0..self.bvals[c].len() {
                if self.status[c][p] == Status::Fixed {
                    let i = Self::fijk(e, p);
                    let v = if self.periodic[c] && i[c] == self.cfg.n[c] { 0.0 } else { self.bc_value(c, i, t) };
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
                        let (lo, hi): (isize, isize) = if a == c { (0, self.cfg.n[a] as isize) } else { (0, self.cfg.n[a] as isize - 1) };
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
                                FluidBoundary::NoSlip(_) | FluidBoundary::Inflow(_) => {
                                    trip.push((r, r, w));
                                    let mut pos = Self::face_pos(&self.cfg, c, i);
                                    pos[a] = if max { self.cfg.len[a] } else { 0.0 };
                                    terms.push((r, 2.0 * w, Term::Wall { face, pos }));
                                },
                                _ => trip.push((r, r, -w)),
                            }
                        }
                    }
                }
            }
            let a = Csr::from_triplets(rows.len(), rows.len(), &trip);
            out.push(Helm { sys: LinearSystem::new(a, LinearSolver::Auto)?, rows, terms });
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
                    } else if matches!(cfg.bcs.get(Face::of_axis(a, ni > 0)), FluidBoundary::Outflow(_)) {
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
    fn extended(&self, c: usize, t: f64) -> (Vec<f64>, [usize; 3]) {
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
                            FluidBoundary::NoSlip(v) | FluidBoundary::Inflow(v) => {
                                let mut pos = Self::face_pos(cfg, c, src);
                                pos[a] = if side { cfg.len[a] } else { 0.0 };
                                2.0 * v[c].eval(t, pos) - uin
                            },
                            _ => uin,
                        }
                    };
                    i[a] = g_idx;
                    out[eid([i[0] + if a == 0 { 0 } else { off[0] }, i[1] + if a == 1 { 0 } else { off[1] }, i[2] + if a == 2 { 0 } else { off[2] }])] = value;
                }
            }
        }
        (out, d)
    }

    /// Computes convection and diffusion (explicit part) at all unknown faces.
    fn convection_and_laplacian(&self, ext: &[(Vec<f64>, [usize; 3])]) -> ([Vec<f64>; 3], [Vec<f64>; 3]) {
        let cfg = &self.cfg;
        let dim = cfg.dim;
        let mut conv: [Vec<f64>; 3] = std::array::from_fn(|_| Vec::new());
        let mut lap: [Vec<f64>; 3] = std::array::from_fn(|_| Vec::new());
        let get = |c: usize, i: [isize; 3]| -> f64 {
            let (arr, d) = &ext[c];
            let mut ix = [0usize; 3];
            for a in 0..3 {
                ix[a] = (i[a] + isize::from(a < dim)) as usize;
            }
            arr[ix[0] + d[0] * (ix[1] + d[1] * ix[2])]
        };
        for c in 0..dim {
            let e = Self::ext(cfg, c);
            let total: usize = e.iter().product();
            conv[c] = vec![0.0; total];
            lap[c] = vec![0.0; total];
            for p in 0..total {
                if self.status[c][p] != Status::Unknown {
                    continue;
                }
                let i = Self::fijk(e, p);
                let ii = [i[0] as isize, i[1] as isize, i[2] as isize];
                let mut cv = 0.0;
                let mut lp = 0.0;
                let uc = get(c, ii);
                for a in 0..dim {
                    let h = self.h(a);
                    let mut ep = [0isize; 3];
                    ep[a] = 1;
                    let shift = |base: [isize; 3], d: [isize; 3], s: isize| [base[0] + s * d[0], base[1] + s * d[1], base[2] + s * d[2]];
                    let (thm, thp) = self.ib_theta(c, p, a);
                    let vp_raw = get(c, shift(ii, ep, 1));
                    let vm_raw = get(c, shift(ii, ep, -1));
                    let (vm, vp) = (self.ib_cross(c, p, a, 0).unwrap_or(vm_raw), self.ib_cross(c, p, a, 1).unwrap_or(vp_raw));
                    lp += 2.0 / (h * h) * (vm / (thm * (thm + thp)) + vp / (thp * (thm + thp)) - uc / (thm * thp));
                    let gp = self.ib_cross(c, p, a, 1).map_or(vp_raw, |w| uc + (w - uc) / thp.max(0.3));
                    let gm = self.ib_cross(c, p, a, 0).map_or(vm_raw, |w| uc + (w - uc) / thm.max(0.3));
                    if a == c {
                        let up = f64::midpoint(uc, gp);
                        let um = f64::midpoint(gm, uc);
                        cv += (up * up - um * um) / h;
                    } else {
                        let mut ec = [0isize; 3];
                        ec[c] = 1;
                        // u_a at the edges above/below the face
                        let top_a = f64::midpoint(get(a, shift(shift(ii, ep, 1), ec, -1)), get(a, shift(ii, ep, 1)));
                        let bot_a = f64::midpoint(get(a, shift(ii, ec, -1)), get(a, ii));
                        let top_c = f64::midpoint(uc, gp);
                        let bot_c = f64::midpoint(gm, uc);
                        cv += (top_a * top_c - bot_a * bot_c) / h;
                    }
                }
                conv[c][p] = cv;
                lap[c][p] = lp;
            }
        }
        (conv, lap)
    }

    /// Pressure gradient component at an unknown face.
    fn pressure_gradient(&self, c: usize, i: [usize; 3], p: &[f64], bc_pressure: bool) -> f64 {
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
                FluidBoundary::Outflow(v) if bc_pressure => v.eval(self.t, Self::face_pos(cfg, c, i)),
                FluidBoundary::Outflow(_) => 0.0,
                _ => cell(i),
            };
            return (cell(i) - pout) / (0.5 * h);
        }
        if i[c] == n[c] {
            let mut q = i;
            q[c] = n[c] - 1;
            let pout = match cfg.bcs.get(Face::of_axis(c, true)) {
                FluidBoundary::Outflow(v) if bc_pressure => v.eval(self.t, Self::face_pos(cfg, c, i)),
                FluidBoundary::Outflow(_) => 0.0,
                _ => cell(q),
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
                    (ib.aperture[c][Self::fidx(e, j)], ib.aperture[c][Self::fidx(e, i)])
                });
                s += (ap_hi * self.u[c][Self::fidx(e, j)] - ap_lo * self.u[c][Self::fidx(e, i)]) / self.h(c);
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
        let ext: Vec<(Vec<f64>, [usize; 3])> = (0..cfg_dim).map(|c| self.extended(c, self.t)).collect();
        let (conv, lap) = self.convection_and_laplacian(&ext);
        let adv: [Vec<f64>; 3] = std::array::from_fn(|c| {
            if c >= cfg_dim {
                return Vec::new();
            }
            self.conv_prev.as_ref().map_or_else(
                || conv[c].clone(),
                |prev| conv[c].iter().zip(&prev[c]).map(|(a, b)| 1.5 * a - 0.5 * b).collect(),
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
                let f = self.cfg.force.as_ref().map_or(0.0, |ff| ff(self.t + 0.5 * dt, x)[c]);
                let g = self.pressure_gradient(c, i, &self.p, true);
                rhs[r] = self.u[c][p] + dt * (-adv[c][p] - g + f) + 0.5 * nu * dt * lap[c][p];
            }
            for (r, coef, term) in &helm.terms {
                rhs[*r] += coef
                    * match term {
                        Term::Const(v) => *v,
                        Term::Fixed(idx) => self.bvals[c][*idx],
                        Term::Wall { face, pos } => match self.cfg.bcs.get(*face) {
                            FluidBoundary::NoSlip(v) | FluidBoundary::Inflow(v) => v[c].eval(t_new, *pos),
                            _ => 0.0,
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
        let rhs: Vec<f64> = (0..ncell).map(|p| if self.pfluid[p] { div[p] / dt } else { 0.0 }).collect();
        let mut rhs = rhs;
        if !(0..self.cfg.dim).any(|a| {
            !self.periodic[a]
                && (matches!(self.cfg.bcs.get(Face::of_axis(a, false)), FluidBoundary::Outflow(_))
                    || matches!(self.cfg.bcs.get(Face::of_axis(a, true)), FluidBoundary::Outflow(_)))
        })
            && let Some(p0) = (0..ncell).find(|&p| self.pfluid[p]) {
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
                self.p[p] += phi[p] - if self.cfg.rotational { 0.5 * nu * div[p] } else { 0.0 };
            }
        }
        let _ = old;
        self.conv_prev = Some(std::array::from_fn(|c| conv[c].clone()));
        self.t = t_new;
        self.steps += 1;
        Ok(())
    }

    /// Runs `steps` steps.
    ///
    /// # Errors
    /// See [`NsSolver::step`].
    pub fn run(&mut self, steps: usize) -> Result<(), String> {
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
    pub fn cell_size(&self, a: usize) -> f64 {
        self.h(a)
    }

    /// Multilinear interpolation of velocity component `c` at `x`.
    #[must_use]
    pub fn sample_velocity(&self, c: usize, x: [f64; 3]) -> f64 {
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
    pub const fn cell_index(&self, i: [usize; 3]) -> usize {
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
                let wgt = if i[c] == 0 || i[c] == cfg.n[c] { 0.5 } else { 1.0 };
                s += wgt * self.u[c][p] * self.u[c][p];
            }
        }
        0.5 * s * vol
    }

    /// Volumetric flux through the plane `x_axis = index` (face index).
    #[must_use]
    pub fn flux_through(&self, axis: usize, face_index: usize) -> f64 {
        let cfg = &self.cfg;
        let e = Self::ext(cfg, axis);
        let mut s = 0.0;
        for p in 0..self.u[axis].len() {
            let i = Self::fijk(e, p);
            if i[axis] == face_index {
                let area: f64 = (0..cfg.dim).filter(|&a| a != axis).map(|a| self.h(a)).product();
                s += self.u[axis][p] * area;
            }
        }
        s
    }

    /// Maximum speed over the faces (for CFL estimates).
    #[must_use]
    pub fn max_speed(&self) -> f64 {
        self.u.iter().flat_map(|v| v.iter()).fold(0.0_f64, |m, x| m.max(x.abs()))
    }
}
