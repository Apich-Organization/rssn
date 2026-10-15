//! # High-order finite volumes for Euler and shallow-water systems in 1D–3D
//!
//! Unsplit finite-volume solver on uniform Cartesian grids with
//!
//! * **reconstruction** of the primitive variables: first order, **MUSCL**
//!   with a minmod, van Leer, monotonised-central or superbee limiter, or
//!   fifth-order **WENO5** (Jiang–Shu);
//! * **Riemann solvers**: **HLLC**, **Roe** (Harten entropy fix) and
//!   Rusanov, for the Euler equations of an ideal gas and the shallow-water
//!   equations;
//! * **SSP-RK3** time stepping with a CFL condition summed over axes;
//! * boundary treatment per face: transmissive, reflective wall, periodic, a
//!   fixed (inflow) state or a linearised characteristic **non-reflecting**
//!   boundary (Euler only).
//!
//! | Equation | 1D | 2D | 3D | Riemann solvers | Reconstruction |
//! |---|---|---|---|---|---|
//! | Euler (`ρ, ρu, ρv, ρw, E`) | yes | yes | yes | HLLC, Roe, Rusanov | first order, MUSCL (4 limiters), WENO5 |
//! | Shallow water (`h, hu, hv`) | yes | yes | n/a | HLLC, Roe, Rusanov | same |
//!
//! Variables are stored as `[f64; 5]`: Euler uses all five slots, shallow
//! water the first `1 + dim`. Initial data are given in primitive variables
//! `[ρ, u, v, w, p]` (Euler) or `[h, u, v, 0, 0]` (shallow water).
//!
//! Limitations: no source terms (bathymetry, gravity), no dry-bed treatment
//! beyond a small depth floor, uniform grids only, no adaptive mesh.

/// Conserved or primitive state with up to five components.
pub type State = [f64; 5];

/// Slope limiter for MUSCL reconstruction.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SlopeLimiter {
    /// Minmod.
    Minmod,
    /// Van Leer.
    VanLeer,
    /// Monotonised central.
    Mc,
    /// Superbee.
    Superbee,
}

/// Spatial reconstruction.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Reconstruction {
    /// Piecewise constant.
    FirstOrder,
    /// MUSCL on primitive variables.
    Muscl(SlopeLimiter),
    /// Fifth-order WENO on primitive variables.
    Weno5,
}

/// Approximate Riemann solver.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RiemannSolver {
    /// HLLC (contact-resolving).
    Hllc,
    /// Roe linearisation with entropy fix.
    Roe,
    /// Local Lax–Friedrichs.
    Rusanov,
}

/// Conservation law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Law {
    /// Euler equations with adiabatic exponent `gamma`.
    Euler {
        /// Adiabatic exponent.
        gamma: f64,
    },
    /// Shallow water with gravity `g`.
    ShallowWater {
        /// Gravitational acceleration.
        g: f64,
    },
}

/// Boundary condition on one face of the box.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum GasBc {
    /// Zero gradient (outflow).
    Transmissive,
    /// Slip wall: mirrored normal velocity.
    Reflective,
    /// Periodic (must be set on both faces of the axis).
    Periodic,
    /// Fixed primitive state (inflow).
    Fixed(State),
    /// Linearised characteristic non-reflecting boundary about the given
    /// primitive far-field state (Euler only): the incoming acoustic
    /// characteristic is set to its far-field value, the outgoing ones are
    /// taken from the interior.
    NonReflecting(State),
}

/// Cartesian cell grid.
#[derive(Clone, Copy, Debug)]
pub struct FvGrid {
    /// Active axes.
    pub dim: usize,
    /// Cells per axis.
    pub n: [usize; 3],
    /// Lower corner.
    pub lo: [f64; 3],
    /// Upper corner.
    pub hi: [f64; 3],
}

impl FvGrid {
    /// Creates a grid.
    ///
    /// # Errors
    /// For a bad dimension or too few cells.
    pub fn new(n: &[usize], lo: &[f64], hi: &[f64]) -> Result<Self, String> {
        let dim = n.len();
        if !(1..=3).contains(&dim) || lo.len() != dim || hi.len() != dim || n.iter().any(|&c| c < 2) {
            return Err("invalid finite-volume grid".into());
        }
        let mut g = Self { dim, n: [1; 3], lo: [0.0; 3], hi: [1.0; 3] };
        for a in 0..dim {
            g.n[a] = n[a];
            g.lo[a] = lo[a];
            g.hi[a] = hi[a];
        }
        Ok(g)
    }

    /// Cell width along `axis`.
    #[must_use]
    pub fn h(&self, axis: usize) -> f64 {
        (self.hi[axis] - self.lo[axis]) / self.n[axis] as f64
    }

    /// Number of cells.
    #[must_use]
    pub fn len(&self) -> usize {
        self.n.iter().product()
    }

    /// Whether there are no cells.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Cell volume.
    #[must_use]
    pub fn cell_volume(&self) -> f64 {
        (0..self.dim).map(|a| self.h(a)).product()
    }

    /// Flat index of cell `(i, j, k)`.
    #[must_use]
    pub const fn index(&self, i: [usize; 3]) -> usize {
        i[0] + self.n[0] * (i[1] + self.n[1] * i[2])
    }

    /// Cell multi-index.
    #[must_use]
    pub const fn ijk(&self, p: usize) -> [usize; 3] {
        let r = p / self.n[0];
        [p % self.n[0], r % self.n[1], r / self.n[1]]
    }

    /// Cell centre.
    #[must_use]
    pub fn centre(&self, p: usize) -> [f64; 3] {
        let ijk = self.ijk(p);
        let mut x = [0.0; 3];
        for a in 0..self.dim {
            x[a] = self.lo[a] + (ijk[a] as f64 + 0.5) * self.h(a);
        }
        x
    }
}

/// Solver configuration.
#[derive(Clone, Copy, Debug)]
pub struct FvConfig {
    /// Equation.
    pub law: Law,
    /// Reconstruction.
    pub recon: Reconstruction,
    /// Riemann solver.
    pub riemann: RiemannSolver,
    /// CFL number (`< 1`).
    pub cfl: f64,
    /// Face boundary conditions in the order `XMin, XMax, YMin, YMax, ZMin, ZMax`.
    pub bcs: [GasBc; 6],
}

impl FvConfig {
    /// Euler solver with all faces transmissive, HLLC + WENO5, CFL 0.4.
    #[must_use]
    pub const fn euler(gamma: f64) -> Self {
        Self {
            law: Law::Euler { gamma },
            recon: Reconstruction::Weno5,
            riemann: RiemannSolver::Hllc,
            cfl: 0.4,
            bcs: [GasBc::Transmissive; 6],
        }
    }

    /// Shallow-water solver with all faces transmissive, HLLC + MUSCL(minmod).
    #[must_use]
    pub const fn shallow_water(g: f64) -> Self {
        Self {
            law: Law::ShallowWater { g },
            recon: Reconstruction::Muscl(SlopeLimiter::Minmod),
            riemann: RiemannSolver::Hllc,
            cfl: 0.4,
            bcs: [GasBc::Transmissive; 6],
        }
    }

    /// Sets both faces of an axis.
    #[must_use]
    pub const fn with_axis(mut self, axis: usize, bc: GasBc) -> Self {
        self.bcs[2 * axis] = bc;
        self.bcs[2 * axis + 1] = bc;
        self
    }

    /// Sets one face (`2 * axis + is_max`).
    #[must_use]
    pub const fn with_face(mut self, face: usize, bc: GasBc) -> Self {
        self.bcs[face] = bc;
        self
    }
}

/// Final state of a run.
#[derive(Clone, Debug)]
pub struct FvSolution {
    /// Grid.
    pub grid: FvGrid,
    /// Conserved variables per cell.
    pub q: Vec<State>,
    /// Final time.
    pub t: f64,
    /// Number of time steps.
    pub steps: usize,
    law: Law,
}

impl FvSolution {
    /// Primitive variables of cell `p`.
    #[must_use]
    pub fn primitive(&self, p: usize) -> State {
        to_prim(self.law, self.q[p])
    }

    /// Integrals of the conserved variables over the domain.
    #[must_use]
    pub fn totals(&self) -> State {
        let v = self.grid.cell_volume();
        let mut s = [0.0; 5];
        for q in &self.q {
            for c in 0..5 {
                s[c] += q[c] * v;
            }
        }
        s
    }

    /// Discrete L1 norm of `component(prim) - exact(x)`.
    pub fn l1_error(&self, component: usize, exact: impl Fn([f64; 3]) -> f64) -> f64 {
        let v = self.grid.cell_volume();
        (0..self.q.len()).map(|p| (self.primitive(p)[component] - exact(self.grid.centre(p))).abs() * v).sum()
    }
}

fn to_cons(law: Law, w: State) -> State {
    match law {
        Law::Euler { gamma } => {
            let ke = 0.5 * w[0] * (w[1] * w[1] + w[2] * w[2] + w[3] * w[3]);
            [w[0], w[0] * w[1], w[0] * w[2], w[0] * w[3], w[4] / (gamma - 1.0) + ke]
        },
        Law::ShallowWater { .. } => [w[0], w[0] * w[1], w[0] * w[2], 0.0, 0.0],
    }
}

fn to_prim(law: Law, q: State) -> State {
    match law {
        Law::Euler { gamma } => {
            let r = q[0];
            let (u, v, w) = (q[1] / r, q[2] / r, q[3] / r);
            [r, u, v, w, (gamma - 1.0) * (q[4] - 0.5 * r * (u * u + v * v + w * w))]
        },
        Law::ShallowWater { .. } => {
            let h = q[0];
            if h < 1e-10 { [h.max(0.0), 0.0, 0.0, 0.0, 0.0] } else { [h, q[1] / h, q[2] / h, 0.0, 0.0] }
        },
    }
}

fn admissible(law: Law, w: &State) -> bool {
    match law {
        Law::Euler { .. } => w[0] > 0.0 && w[4] > 0.0,
        Law::ShallowWater { .. } => w[0] > 0.0,
    }
}

fn sound_speed(law: Law, w: &State) -> f64 {
    match law {
        Law::Euler { gamma } => (gamma * w[4] / w[0]).max(0.0).sqrt(),
        Law::ShallowWater { g } => (g * w[0]).max(0.0).sqrt(),
    }
}

/// Physical flux in the first (normal) direction of a primitive state.
fn flux_x(law: Law, w: &State) -> State {
    match law {
        Law::Euler { gamma } => {
            let q = to_cons(Law::Euler { gamma }, *w);
            let u = w[1];
            [q[1], q[1] * u + w[4], q[2] * u, q[3] * u, (q[4] + w[4]) * u]
        },
        Law::ShallowWater { g } => {
            let (h, u, v) = (w[0], w[1], w[2]);
            [h * u, h * u * u + 0.5 * g * h * h, h * u * v, 0.0, 0.0]
        },
    }
}

const fn swap_dir(s: &mut State, d: usize) {
    if d > 0 {
        s.swap(1, 1 + d);
    }
}

/// Numerical flux across a face with normal along `x` (states already rotated).
fn riemann_flux(law: Law, solver: RiemannSolver, l: &State, r: &State) -> State {
    let fl = flux_x(law, l);
    let fr = flux_x(law, r);
    let (cl, cr) = (sound_speed(law, l), sound_speed(law, r));
    let (ql, qr) = (to_cons(law, *l), to_cons(law, *r));
    // Roe averages for wave-speed estimates
    let (sl_, sr_) = (l[0].sqrt(), r[0].sqrt());
    let ut = (sl_ * l[1] + sr_ * r[1]) / (sl_ + sr_);
    let ct = match law {
        Law::Euler { gamma } => {
            let (hl, hr) = ((ql[4] + l[4]) / l[0], (qr[4] + r[4]) / r[0]);
            let ht = (sl_ * hl + sr_ * hr) / (sl_ + sr_);
            let vt = (sl_ * l[2] + sr_ * r[2]) / (sl_ + sr_);
            let wt = (sl_ * l[3] + sr_ * r[3]) / (sl_ + sr_);
            ((gamma - 1.0) * (ht - 0.5 * (ut * ut + vt * vt + wt * wt))).max(1e-300).sqrt()
        },
        Law::ShallowWater { g } => (0.5 * g * (l[0] + r[0])).sqrt(),
    };
    let sl = (l[1] - cl).min(ut - ct);
    let sr = (r[1] + cr).max(ut + ct);
    let nv = match law {
        Law::Euler { .. } => 5,
        Law::ShallowWater { .. } => 3,
    };
    match solver {
        RiemannSolver::Rusanov => {
            let s = (l[1].abs() + cl).max(r[1].abs() + cr);
            let mut f = [0.0; 5];
            for c in 0..nv {
                f[c] = f64::midpoint(fl[c], fr[c]) - 0.5 * s * (qr[c] - ql[c]);
            }
            f
        },
        RiemannSolver::Hllc => {
            if sl >= 0.0 {
                return fl;
            }
            if sr <= 0.0 {
                return fr;
            }
            let (ul, ur) = (l[1], r[1]);
            let den = l[0] * (sl - ul) - r[0] * (sr - ur);
            let sm = match law {
                Law::Euler { .. } => (r[4] - l[4] + l[0] * ul * (sl - ul) - r[0] * ur * (sr - ur)) / den,
                Law::ShallowWater { .. } => (sl * r[0] * (ur - sr) - sr * l[0] * (ul - sl)) / (r[0] * (ur - sr) - l[0] * (ul - sl)),
            };
            let star = |w: &State, q: &State, s: f64| -> State {
                let k = w[0] * (s - w[1]) / (s - sm);
                match law {
                    Law::Euler { .. } => [
                        k,
                        k * sm,
                        k * w[2],
                        k * w[3],
                        k * (q[4] / w[0] + (sm - w[1]) * (sm + w[4] / (w[0] * (s - w[1])))),
                    ],
                    Law::ShallowWater { .. } => [k, k * sm, k * w[2], 0.0, 0.0],
                }
            };
            let (f, w, q, s) = if sm >= 0.0 { (&fl, l, &ql, sl) } else { (&fr, r, &qr, sr) };
            let qs = star(w, q, s);
            let mut out = [0.0; 5];
            for c in 0..nv {
                out[c] = f[c] + s * (qs[c] - q[c]);
            }
            out
        },
        RiemannSolver::Roe => {
            let mut diss = [0.0; 5];
            let fix = |lam: f64, delta: f64| {
                let a = lam.abs();
                if a < delta { f64::midpoint(lam * lam / delta, delta) } else { a }
            };
            match law {
                Law::Euler { gamma } => {
                    let rt = sl_ * sr_;
                    let (hl, hr) = ((ql[4] + l[4]) / l[0], (qr[4] + r[4]) / r[0]);
                    let ht = (sl_ * hl + sr_ * hr) / (sl_ + sr_);
                    let vt = (sl_ * l[2] + sr_ * r[2]) / (sl_ + sr_);
                    let wt = (sl_ * l[3] + sr_ * r[3]) / (sl_ + sr_);
                    let c = ct;
                    let _ = gamma;
                    let (dr, dp, du, dv, dw) = (r[0] - l[0], r[4] - l[4], r[1] - l[1], r[2] - l[2], r[3] - l[3]);
                    let a1 = (dp - rt * c * du) / (2.0 * c * c);
                    let a2 = dr - dp / (c * c);
                    let a3 = (dp + rt * c * du) / (2.0 * c * c);
                    let q2 = 0.5 * (ut * ut + vt * vt + wt * wt);
                    let delta = 0.1 * c;
                    let (l1, l2, l3) = (fix(ut - c, delta), ut.abs(), fix(ut + c, delta));
                    let r1 = [1.0, ut - c, vt, wt, ht - ut * c];
                    let r2 = [1.0, ut, vt, wt, q2];
                    let r3 = [1.0, ut + c, vt, wt, ht + ut * c];
                    let r4 = [0.0, 0.0, 1.0, 0.0, vt];
                    let r5 = [0.0, 0.0, 0.0, 1.0, wt];
                    for k in 0..5 {
                        diss[k] = l1 * a1 * r1[k]
                            + l2 * a2 * r2[k]
                            + l3 * a3 * r3[k]
                            + l2 * rt * dv * r4[k]
                            + l2 * rt * dw * r5[k];
                    }
                },
                Law::ShallowWater { .. } => {
                    let vt = (sl_ * l[2] + sr_ * r[2]) / (sl_ + sr_);
                    let ht = f64::midpoint(l[0], r[0]);
                    let c = ct;
                    let (dh, dm) = (r[0] - l[0], r[0] * r[1] - l[0] * l[1]);
                    let a1 = ((ut + c) * dh - dm) / (2.0 * c);
                    let a3 = (dm - (ut - c) * dh) / (2.0 * c);
                    let delta = 0.1 * c;
                    let (l1, l2, l3) = (fix(ut - c, delta), ut.abs(), fix(ut + c, delta));
                    let r1 = [1.0, ut - c, vt];
                    let r3 = [1.0, ut + c, vt];
                    let dv = r[2] - l[2];
                    for k in 0..3 {
                        diss[k] = l1 * a1 * r1[k] + l3 * a3 * r3[k];
                    }
                    diss[2] += l2 * ht * dv;
                },
            }
            let mut f = [0.0; 5];
            for c in 0..nv {
                f[c] = f64::midpoint(fl[c], fr[c]) - 0.5 * diss[c];
            }
            f
        },
    }
}

fn limit(lim: SlopeLimiter, a: f64, b: f64) -> f64 {
    // a = backward difference, b = forward difference
    if a * b <= 0.0 {
        return 0.0;
    }
    match lim {
        SlopeLimiter::Minmod => if a.abs() < b.abs() { a } else { b },
        SlopeLimiter::VanLeer => 2.0 * a * b / (a + b),
        SlopeLimiter::Mc => {
            let c = f64::midpoint(a, b);
            let m = (2.0 * a.abs()).min(2.0 * b.abs()).min(c.abs());
            m.copysign(a)
        },
        SlopeLimiter::Superbee => {
            let s1 = if a.abs() < 2.0 * b.abs() { a } else { 2.0 * b };
            let s2 = if 2.0 * a.abs() < b.abs() { 2.0 * a } else { b };
            if s1.abs() > s2.abs() { s1 } else { s2 }
        },
    }
}

/// Fifth-order WENO value at the right face of the cell holding `v2`
/// (stencil `v0 .. v4`, centred on `v2`).
fn weno5(v0: f64, v1: f64, v2: f64, v3: f64, v4: f64) -> f64 {
    let q0 = (2.0 * v0 - 7.0 * v1 + 11.0 * v2) / 6.0;
    let q1 = (-v1 + 5.0 * v2 + 2.0 * v3) / 6.0;
    let q2 = (2.0 * v2 + 5.0 * v3 - v4) / 6.0;
    let b0 = 13.0 / 12.0 * (v0 - 2.0 * v1 + v2).powi(2) + 0.25 * (v0 - 4.0 * v1 + 3.0 * v2).powi(2);
    let b1 = 13.0 / 12.0 * (v1 - 2.0 * v2 + v3).powi(2) + 0.25 * (v1 - v3).powi(2);
    let b2 = 13.0 / 12.0 * (v2 - 2.0 * v3 + v4).powi(2) + 0.25 * (3.0 * v2 - 4.0 * v3 + v4).powi(2);
    let eps = 1e-6;
    let a0 = 0.1 / (eps + b0).powi(2);
    let a1 = 0.6 / (eps + b1).powi(2);
    let a2 = 0.3 / (eps + b2).powi(2);
    (a0 * q0 + a1 * q1 + a2 * q2) / (a0 + a1 + a2)
}

/// Left and right states at the face between line cells `k - 1` and `k`.
fn face_states(recon: Reconstruction, line: &[State], k: usize, nv: usize) -> (State, State) {
    let mut l = line[k - 1];
    let mut r = line[k];
    match recon {
        Reconstruction::FirstOrder => {},
        Reconstruction::Muscl(lim) => {
            for c in 0..nv {
                let sl = limit(lim, line[k - 1][c] - line[k - 2][c], line[k][c] - line[k - 1][c]);
                let sr = limit(lim, line[k][c] - line[k - 1][c], line[k + 1][c] - line[k][c]);
                l[c] = line[k - 1][c] + 0.5 * sl;
                r[c] = line[k][c] - 0.5 * sr;
            }
        },
        Reconstruction::Weno5 => {
            for c in 0..nv {
                l[c] = weno5(line[k - 3][c], line[k - 2][c], line[k - 1][c], line[k][c], line[k + 1][c]);
                r[c] = weno5(line[k + 2][c], line[k + 1][c], line[k][c], line[k - 1][c], line[k - 2][c]);
            }
        },
    }
    (l, r)
}

/// Ghost state of the linearised characteristic non-reflecting boundary.
fn nrbc_ghost(law: Law, q: State, w_inf: &State, axis: usize, max: bool) -> State {
    let Law::Euler { gamma } = law else { return q };
    let wi = to_prim(law, q);
    let c = (gamma * w_inf[4] / w_inf[0]).sqrt();
    let rc = w_inf[0] * c;
    let sgn = if max { 1.0 } else { -1.0 };
    let (un_i, un_inf) = (sgn * wi[1 + axis], sgn * w_inf[1 + axis]);
    let rp = f64::midpoint(wi[4] - w_inf[4], rc * (un_i - un_inf));
    let p_g = (w_inf[4] + rp).max(1e-8);
    let un_g = un_inf + rp / rc;
    let mut wg = if un_inf >= 0.0 { wi } else { *w_inf };
    wg[0] = (wg[0] + (p_g - if un_inf >= 0.0 { wi[4] } else { w_inf[4] }) / (c * c)).max(1e-8);
    wg[4] = p_g;
    wg[1 + axis] = sgn * un_g;
    to_cons(law, wg)
}

const GHOST: usize = 3;

struct Padded {
    np: [usize; 3],
    g: [usize; 3],
}

impl Padded {
    fn new(grid: &FvGrid) -> Self {
        let mut np = [1; 3];
        let mut g = [0; 3];
        for a in 0..grid.dim {
            np[a] = grid.n[a] + 2 * GHOST;
            g[a] = GHOST;
        }
        Self { np, g }
    }

    const fn idx(&self, i: [usize; 3]) -> usize {
        i[0] + self.np[0] * (i[1] + self.np[1] * i[2])
    }

    fn len(&self) -> usize {
        self.np.iter().product()
    }
}

struct Stepper<'a> {
    grid: &'a FvGrid,
    cfg: &'a FvConfig,
    pad: Padded,
    nv: usize,
}

impl Stepper<'_> {
    fn fill_padded(&self, q: &[State]) -> Vec<State> {
        let g = self.grid;
        let pad = &self.pad;
        let law = self.cfg.law;
        let mut p = vec![[0.0; 5]; pad.len()];
        for c in 0..g.len() {
            let ijk = g.ijk(c);
            let pi = [ijk[0] + pad.g[0], ijk[1] + pad.g[1], ijk[2] + pad.g[2]];
            p[pad.idx(pi)] = q[c];
        }
        for a in 0..g.dim {
            let n = g.n[a] as isize;
            let gh = GHOST as isize;
            let ranges: Vec<std::ops::Range<usize>> = (0..3)
                .map(|b| {
                    if b == a {
                        0..0
                    } else if b < a || b >= g.dim {
                        0..pad.np[b]
                    } else {
                        pad.g[b]..pad.g[b] + g.n[b]
                    }
                })
                .collect();
            for i2 in if a == 2 { 0..1 } else { ranges[2].clone() } {
                for i1 in if a == 1 { 0..1 } else { ranges[1].clone() } {
                    for i0 in if a == 0 { 0..1 } else { ranges[0].clone() } {
                        for layer in 0..GHOST {
                            for max in [false, true] {
                                let mut dst = [i0, i1, i2];
                                let mut src = dst;
                                // ghost coordinate in interior index space
                                let gi: isize = if max { n + layer as isize } else { -1 - layer as isize };
                                dst[a] = (gi + gh) as usize;
                                let bc = self.cfg.bcs[2 * a + usize::from(max)];
                                let mut flip = false;
                                let s_i: isize = match bc {
                                    GasBc::Transmissive => if max { n - 1 } else { 0 },
                                    GasBc::Reflective => {
                                        flip = true;
                                        if max { n - 1 - layer as isize } else { layer as isize }
                                    },
                                    GasBc::Periodic => gi.rem_euclid(n),
                                    GasBc::Fixed(w) => {
                                        p[pad.idx(dst)] = to_cons(law, w);
                                        continue;
                                    },
                                    GasBc::NonReflecting(w_inf) => {
                                        let mut edge = dst;
                                        edge[a] = ((if max { n - 1 } else { 0 }) + gh) as usize;
                                        p[pad.idx(dst)] = nrbc_ghost(law, p[pad.idx(edge)], &w_inf, a, max);
                                        continue;
                                    },
                                };
                                src[a] = (s_i + gh) as usize;
                                let mut v = p[pad.idx(src)];
                                if flip {
                                    v[1 + a] = -v[1 + a];
                                }
                                p[pad.idx(dst)] = v;
                            }
                        }
                    }
                }
            }
        }
        p
    }

    fn max_rate(&self, q: &[State]) -> f64 {
        let g = self.grid;
        q.iter()
            .map(|&c| {
                let w = to_prim(self.cfg.law, c);
                let cs = sound_speed(self.cfg.law, &w);
                (0..g.dim).map(|a| (w[1 + a].abs() + cs) / g.h(a)).sum::<f64>()
            })
            .fold(0.0, f64::max)
    }

    fn residual(&self, q: &[State]) -> Vec<State> {
        let g = self.grid;
        let pad = &self.pad;
        let law = self.cfg.law;
        let qp = self.fill_padded(q);
        let wp: Vec<State> = qp.iter().map(|&c| to_prim(law, c)).collect();
        let mut res = vec![[0.0; 5]; g.len()];
        for d in 0..g.dim {
            let n = g.n[d];
            let hd = g.h(d);
            let o: Vec<usize> = (0..3).filter(|&b| b != d).collect();
            let mut line = vec![[0.0; 5]; pad.np[d]];
            let mut flux = vec![[0.0; 5]; n + 1];
            for ib in 0..g.n[o[1]].max(1) {
                for ia in 0..g.n[o[0]].max(1) {
                    let mut base = [0usize; 3];
                    base[o[0]] = ia + pad.g[o[0]];
                    base[o[1]] = ib + pad.g[o[1]];
                    for k in 0..pad.np[d] {
                        let mut i = base;
                        i[d] = k;
                        line[k] = wp[pad.idx(i)];
                    }
                    for f in 0..=n {
                        let (mut l, mut r) = face_states(self.cfg.recon, &line, GHOST + f, self.nv);
                        if !admissible(law, &l) || !admissible(law, &r) {
                            l = line[GHOST + f - 1];
                            r = line[GHOST + f];
                        }
                        swap_dir(&mut l, d);
                        swap_dir(&mut r, d);
                        let mut fl = riemann_flux(law, self.cfg.riemann, &l, &r);
                        swap_dir(&mut fl, d);
                        flux[f] = fl;
                    }
                    for c in 0..n {
                        let mut i = [0usize; 3];
                        i[o[0]] = ia;
                        i[o[1]] = ib;
                        i[d] = c;
                        let cell = g.index(i);
                        for v in 0..5 {
                            res[cell][v] -= (flux[c + 1][v] - flux[c][v]) / hd;
                        }
                    }
                }
            }
        }
        res
    }
}

/// Runs the solver from primitive initial data to `t_end`.
///
/// # Errors
/// For inconsistent periodic faces, a CFL number outside `(0, 1)`, an
/// unsupported dimension for shallow water, or non-admissible initial data.
pub fn solve(
    grid: &FvGrid,
    cfg: &FvConfig,
    init: &dyn Fn([f64; 3]) -> State,
    t_end: f64,
) -> Result<FvSolution, String> {
    if !(cfg.cfl > 0.0 && cfg.cfl < 1.0) {
        return Err("cfl must be in (0, 1)".into());
    }
    if matches!(cfg.law, Law::ShallowWater { .. }) && cfg.bcs.iter().any(|b| matches!(b, GasBc::NonReflecting(_))) {
        return Err("non-reflecting boundaries are available for the Euler equations only".into());
    }
    if matches!(cfg.law, Law::ShallowWater { .. }) && grid.dim == 3 {
        return Err("shallow water is two-dimensional at most".into());
    }
    for a in 0..grid.dim {
        if (cfg.bcs[2 * a] == GasBc::Periodic) != (cfg.bcs[2 * a + 1] == GasBc::Periodic) {
            return Err("periodic must be set on both faces of an axis".into());
        }
    }
    let nv = match cfg.law {
        Law::Euler { .. } => 5,
        Law::ShallowWater { .. } => 1 + grid.dim,
    };
    let mut q: Vec<State> = Vec::with_capacity(grid.len());
    for p in 0..grid.len() {
        let w = init(grid.centre(p));
        if !admissible(cfg.law, &w) {
            return Err("initial state is not admissible".into());
        }
        q.push(to_cons(cfg.law, w));
    }
    let st = Stepper { grid, cfg, pad: Padded::new(grid), nv };
    let mut t = 0.0;
    let mut steps = 0;
    let add = |a: &[State], b: &[State], wa: f64, wb: f64, dt: f64| -> Vec<State> {
        a.iter()
            .zip(b)
            .map(|(x, y)| {
                let mut o = [0.0; 5];
                for c in 0..5 {
                    o[c] = wa * x[c] + wb * (y[c] * dt);
                }
                o
            })
            .collect()
    };
    while t < t_end * (1.0 - 1e-14) {
        let rate = st.max_rate(&q);
        if !rate.is_finite() || rate <= 0.0 {
            return Err("non-finite wave speed".into());
        }
        let dt = (cfg.cfl / rate).min(t_end - t);
        let l0 = st.residual(&q);
        let q1: Vec<State> = add(&q, &l0, 1.0, 1.0, dt);
        let l1 = st.residual(&q1);
        let tmp = add(&q1, &l1, 1.0, 1.0, dt);
        let q2: Vec<State> = q.iter().zip(&tmp).map(|(a, b)| {
            let mut o = [0.0; 5];
            for c in 0..5 { o[c] = 0.75 * a[c] + 0.25 * b[c]; }
            o
        }).collect();
        let l2 = st.residual(&q2);
        let tmp = add(&q2, &l2, 1.0, 1.0, dt);
        q = q.iter().zip(&tmp).map(|(a, b)| {
            let mut o = [0.0; 5];
            for c in 0..5 { o[c] = a[c] / 3.0 + 2.0 * b[c] / 3.0; }
            o
        }).collect();
        t += dt;
        steps += 1;
        if steps > 2_000_000 {
            return Err("too many steps".into());
        }
    }
    Ok(FvSolution { grid: *grid, q, t, steps, law: cfg.law })
}
