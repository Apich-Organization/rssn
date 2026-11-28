//! High-order finite volumes for Euler and shallow water in 1D-3D.
//!
//! # High-order finite volumes for Euler and shallow-water systems in 1D–3D
//!
//! Unsplit finite-volume solver on uniform Cartesian grids with
//!
//! * **reconstruction** of the primitive variables: first order, **MUSCL**
//!   with a minmod, van Leer, monotonised-central or superbee limiter, or
//!   fifth-order **WENO5** (Jiang–Shu);
//! * **Riemann solvers**: **HLLC**, **HLLE**, **Roe** (Harten entropy fix)
//!   and Rusanov, for the Euler equations of an ideal gas and the shallow-water
//!   equations;
//! * **SSP-RK3** time stepping with a CFL condition summed over axes;
//! * boundary treatment per face: transmissive, reflective wall, periodic, a
//!   fixed (inflow) state or a linearised characteristic **non-reflecting**
//!   boundary (Euler only).
//!
//! | Equation | 1D | 2D | 3D | Riemann solvers | Reconstruction |
//! |---|---|---|---|---|---|
//! | Euler (`ρ, ρu, ρv, ρw, E`) | yes | yes | yes | HLLC, HLLE, Roe, Rusanov | first order, MUSCL (4 limiters), WENO5 |
//! | Shallow water (`h, hu, hv, hw`) | yes | yes | yes (barotropic) | HLLC, HLLE, Roe, Rusanov | same |
//!
//! Variables are stored as `[f64; 5]`: Euler uses all five slots, shallow
//! water the first `1 + dim`. Initial data are given in primitive variables
//! `[ρ, u, v, w, p]` (Euler) or `[h, u, v, w, 0]` (shallow water, `w = 0`
//! below three dimensions).
//!
//! The depth-averaged shallow-water model is planar. Its **three-dimensional
//! version** is the barotropic (isentropic, `γ = 2`) system with the pressure
//! law `p = g h² / 2` and state `(h, hu, hv, hw)` — the exact 3D analogue of
//! the shallow-water system (Euler with that pressure law and no energy
//! equation), with the wave speed `c = √(g h)`. A bed and Manning friction are
//! planar notions and are rejected in 3D; in 1D and 2D they work as before.
//!
//! Block-structured adaptive refinement with sub-cycling and refluxing is
//! provided by [`euler_amr`](crate::sim::euler_amr) for Euler and shallow
//! water in 1D–3D, built on the reconstruction and Riemann kernels of this
//! module.
//!
//! ## Source terms ([`solve_with_sources`])
//!
//! * **Gravity** (Euler): a uniform acceleration `g` adds `ρ g` to the
//!   momentum and `ρ u·g` to the energy, evaluated in the cell centres, and
//!   the momentum flux uses a **well-balanced** interface pressure: for a
//!   locally isothermal stratification (`θ = p/ρ` constant across a cell) the
//!   pressure `P = p_L a / (1 - e^{-a})`, `a = ln(p_R / p_L)`, satisfies
//!   `P_{i+1/2} - P_{i-1/2} = ρ_i g Δx` exactly, so an isothermal atmosphere
//!   at rest between two reflective walls is held to round-off (a reflective
//!   wall in the gravity direction uses the hydrostatic continuation of the
//!   boundary cell). A uniform gas in free fall is reproduced exactly (the
//!   interface pressure reduces to the common pressure and the time
//!   integrator is exact for the resulting polynomial in `t`). The
//!   well-balanced interface pressure is second order accurate, so gravity
//!   runs are second order in space even with the high-order reconstruction.
//! * **Bathymetry** (shallow water): a bed `b(x)` with the **hydrostatic
//!   reconstruction** of Audusse, Bouchut, Bristeau, Klein and Perthame:
//!   the free surface `η = h + b`, the velocity and `b` are reconstructed
//!   (so the configured MUSCL/WENO5 scheme applies), the face depths are
//!   `h* = max(0, η - max(b_L, b_R))`, the Riemann solver gets the starred
//!   states, each side's flux is corrected by `g (h² - h*²)/2` in the
//!   normal momentum, and the cell source is the centred
//!   `-g h̄ Δb / Δx`. The lake at rest `η = const, u = 0` is preserved to
//!   round-off by every reconstruction and Riemann solver, over smooth and
//!   discontinuous beds, in 1D and 2D (a bed is rejected in 3D). The initial
//!   data are the *depth* `h = η - b`.
//! * **Manning bed friction** (shallow water in 1D–2D, [`Sources::manning`]): the
//!   momentum source `-g n² h^{-7/3} |q| q` (equivalently
//!   `d(h u)/dt = -g n² h^{-1/3} |u| u`, `q = h u`) with a uniform Manning
//!   coefficient `n`. It is integrated in *operator-split (Strang) form*: a
//!   half step of the friction before and after the SSP-RK3 hyperbolic step.
//!   A friction half step is the exact solution of the source ODE with the
//!   depth held at its current value,
//!   `q ← q / (1 + dt g n² |q| h^{-7/3})`, which is unconditionally stable
//!   and leaves `|q|` monotone; `n = 0` is the frictionless solver and stays
//!   bit for bit unchanged. In a uniform flow over a flat bed this reduces to
//!   `u(t) = u₀ / (1 + t g n² u₀ / h^{4/3})`.
//!
//! Limitations: no dry-bed treatment beyond a small depth floor (the
//! well-balanced reconstruction is positivity-preserving in `h*`, but a
//! front over a dry bed is not resolved sharply); the Euler gravity source is
//! well-balanced only for a locally isothermal stratification (`θ = p/ρ`
//! constant across a cell, so not exactly for a polytropic atmosphere) and
//! only at reflective walls with gravity aligned with the grid axes; the
//! Manning coefficient is uniform (no per-cell roughness map) and friction is
//! not well-balanced against a wet/dry front; this solver itself uses a
//! uniform grid (adaptive refinement is [`euler_amr`](crate::sim::euler_amr)).
//!
//! Independent lines of a flux sweep can run on the active CPU pool. Each pass
//! updates a cell once, and the axis order is fixed across worker counts.

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
    /// Two-wave Harten–Lax–van Leer flux with Einfeldt wave-speed estimates.
    Hlle,
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
    /// Returns an error for an invalid dimension, extent, spacing or storage size.
    pub fn new(
        n: &[usize],
        lo: &[f64],
        hi: &[f64],
    ) -> Result<Self, String> {
        let dim = n.len();
        if !(1..=3).contains(&dim) || lo.len() != dim || hi.len() != dim {
            return Err("invalid finite-volume grid".into());
        }
        let mut g = Self {
            dim,
            n: [1; 3],
            lo: [0.0; 3],
            hi: [1.0; 3],
        };
        g.n[..dim].copy_from_slice(&n[..dim]);
        g.lo[..dim].copy_from_slice(&lo[..dim]);
        g.hi[..dim].copy_from_slice(&hi[..dim]);
        g.check()?;
        Ok(g)
    }

    pub(super) fn check(&self) -> Result<(), String> {
        if !(1..=3).contains(&self.dim)
            || self.n[self.dim..].iter().any(|&n| n != 1)
            || (0..self.dim).any(|axis| {
                let spacing = self.h(axis);
                self.n[axis] < 2 || !spacing.is_finite() || spacing <= 0.0
            })
        {
            return Err("invalid finite-volume grid or spacing".into());
        }
        let capacity = self.n[..self.dim].iter().try_fold(1usize, |size, n| {
            size.checked_mul(n.checked_add(2 * GHOST)?)
        });
        if capacity.is_none_or(|size| size > isize::MAX as usize / size_of::<CellUpdate>()) {
            return Err("finite-volume grid exceeds addressable storage".into());
        }
        Ok(())
    }

    /// Cell width along `axis`.
    #[must_use]
    pub fn h(
        &self,
        axis: usize,
    ) -> f64 {
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
    pub const fn index(
        &self,
        i: [usize; 3],
    ) -> usize {
        i[0] + self.n[0] * (i[1] + self.n[1] * i[2])
    }

    /// Cell multi-index.
    #[must_use]
    pub const fn ijk(
        &self,
        p: usize,
    ) -> [usize; 3] {
        let r = p / self.n[0];
        [p % self.n[0], r % self.n[1], r / self.n[1]]
    }

    /// Cell centre.
    #[must_use]
    pub fn centre(
        &self,
        p: usize,
    ) -> [f64; 3] {
        let ijk = self.ijk(p);
        let mut x = [0.0; 3];
        for a in 0..self.dim {
            x[a] = self.lo[a] + (ijk[a] as f64 + 0.5) * self.h(a);
        }
        x
    }
}

/// Source terms of [`solve_with_sources`].
#[derive(Clone, Copy)]
pub struct Sources<'a> {
    /// Uniform gravitational acceleration of the Euler equations (the
    /// components beyond the grid dimension are ignored).
    pub gravity: [f64; 3],
    /// Bed elevation `b(x)` of the shallow-water equations (hydrostatic
    /// reconstruction, 1D–2D; rejected in 3D); the initial data are then the
    /// depth `h`.
    pub bed: Option<&'a dyn Fn([f64; 3]) -> f64>,
    /// Uniform Manning bed-friction coefficient `n` (shallow water; zero is
    /// frictionless). The momentum source is `-g n² h^{-7/3} |q| q` with
    /// `q = h u`, integrated in Strang-split form with an exact, stable
    /// implicit half step.
    pub manning: f64,
}

impl Default for Sources<'_> {
    fn default() -> Self {
        Self {
            gravity: [0.0; 3],
            bed: None,
            manning: 0.0,
        }
    }
}

impl std::fmt::Debug for Sources<'_> {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        f.debug_struct("Sources")
            .field("gravity", &self.gravity)
            .field("bed", &self.bed.map(|_| "<function>"))
            .field("manning", &self.manning)
            .finish()
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
    pub const fn with_axis(
        mut self,
        axis: usize,
        bc: GasBc,
    ) -> Self {
        self.bcs[2 * axis] = bc;
        self.bcs[2 * axis + 1] = bc;
        self
    }

    /// Sets one face (`2 * axis + is_max`).
    #[must_use]
    pub const fn with_face(
        mut self,
        face: usize,
        bc: GasBc,
    ) -> Self {
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
    pub fn primitive(
        &self,
        p: usize,
    ) -> State {
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
    pub fn l1_error(
        &self,
        component: usize,
        exact: impl Fn([f64; 3]) -> f64,
    ) -> f64 {
        let v = self.grid.cell_volume();
        (0..self.q.len())
            .map(|p| (self.primitive(p)[component] - exact(self.grid.centre(p))).abs() * v)
            .sum()
    }
}

pub(crate) fn to_cons(
    law: Law,
    w: State,
) -> State {
    match law {
        | Law::Euler { gamma } => {
            let ke = 0.5 * w[0] * (w[1] * w[1] + w[2] * w[2] + w[3] * w[3]);
            [
                w[0],
                w[0] * w[1],
                w[0] * w[2],
                w[0] * w[3],
                w[4] / (gamma - 1.0) + ke,
            ]
        },
        | Law::ShallowWater { .. } => [w[0], w[0] * w[1], w[0] * w[2], w[0] * w[3], 0.0],
    }
}

pub(crate) fn to_prim(
    law: Law,
    q: State,
) -> State {
    match law {
        | Law::Euler { gamma } => {
            let r = q[0];
            let (u, v, w) = (q[1] / r, q[2] / r, q[3] / r);
            [
                r,
                u,
                v,
                w,
                (gamma - 1.0) * (q[4] - 0.5 * r * (u * u + v * v + w * w)),
            ]
        },
        | Law::ShallowWater { .. } => {
            let h = q[0];
            if h < 1e-10 {
                [h.max(0.0), 0.0, 0.0, 0.0, 0.0]
            } else {
                [h, q[1] / h, q[2] / h, q[3] / h, 0.0]
            }
        },
    }
}

pub(crate) fn admissible(
    law: Law,
    w: &State,
) -> bool {
    match law {
        | Law::Euler { .. } => w[0] > 0.0 && w[4] > 0.0,
        | Law::ShallowWater { .. } => w[0] > 0.0,
    }
}

fn finite_primitive(
    law: Law,
    state: &State,
) -> bool {
    let used = if matches!(law, Law::Euler { .. }) {
        5
    } else {
        4
    };
    admissible(law, state) && state[..used].iter().all(|value| value.is_finite())
}

pub(crate) fn sound_speed(
    law: Law,
    w: &State,
) -> f64 {
    match law {
        | Law::Euler { gamma } => (gamma * w[4] / w[0]).max(0.0).sqrt(),
        | Law::ShallowWater { g } => (g * w[0]).max(0.0).sqrt(),
    }
}

/// Well-balanced Euler interface pressure from the two adjacent cell-centre
/// pressures: with `a = ln(p_R / p_L)`, `P = p_L a / (1 - e^{-a})` (the limit
/// `(p_L + p_R)/2` at `p_L = p_R`). `P` reduces to the common pressure of a
/// uniform state — so a gas in free fall is untouched — and it is a second
/// order accurate approximation of the face pressure for a general smooth
/// state.
fn hydrostatic_face_pressure(
    p_left: f64,
    p_right: f64,
) -> f64 {
    if p_left > 0.0 && p_right > 0.0 {
        let a = (p_right / p_left).ln();
        if a.abs() > 1e-8 {
            return p_left * a / -(-a).exp_m1();
        }
    }
    f64::midpoint(p_left, p_right)
}

/// Physical flux in the first (normal) direction of a primitive state.
fn flux_x(
    law: Law,
    w: &State,
) -> State {
    match law {
        | Law::Euler { gamma } => {
            let q = to_cons(Law::Euler { gamma }, *w);
            let u = w[1];
            [q[1], q[1] * u + w[4], q[2] * u, q[3] * u, (q[4] + w[4]) * u]
        },
        | Law::ShallowWater { g } => {
            let (h, u, v, ww) = (w[0], w[1], w[2], w[3]);
            [
                h * u,
                h * u * u + 0.5 * g * h * h,
                h * u * v,
                h * u * ww,
                0.0,
            ]
        },
    }
}

pub(crate) const fn swap_dir(
    s: &mut State,
    d: usize,
) {
    if d > 0 {
        s.swap(1, 1 + d);
    }
}

/// Numerical flux across a face with normal along `x` (states already rotated).
pub(crate) fn riemann_flux(
    law: Law,
    solver: RiemannSolver,
    l: &State,
    r: &State,
) -> State {
    let fl = flux_x(law, l);
    let fr = flux_x(law, r);
    let (cl, cr) = (sound_speed(law, l), sound_speed(law, r));
    let (ql, qr) = (to_cons(law, *l), to_cons(law, *r));
    // Roe averages for wave-speed estimates
    let (sl_, sr_) = (l[0].sqrt(), r[0].sqrt());
    let ut = (sl_ * l[1] + sr_ * r[1]) / (sl_ + sr_);
    let ct = match law {
        | Law::Euler { gamma } => {
            let (hl, hr) = ((ql[4] + l[4]) / l[0], (qr[4] + r[4]) / r[0]);
            let ht = (sl_ * hl + sr_ * hr) / (sl_ + sr_);
            let vt = (sl_ * l[2] + sr_ * r[2]) / (sl_ + sr_);
            let wt = (sl_ * l[3] + sr_ * r[3]) / (sl_ + sr_);
            ((gamma - 1.0) * (ht - 0.5 * (ut * ut + vt * vt + wt * wt)))
                .max(1e-300)
                .sqrt()
        },
        | Law::ShallowWater { g } => (0.5 * g * (l[0] + r[0])).sqrt(),
    };
    let sl = (l[1] - cl).min(ut - ct);
    let sr = (r[1] + cr).max(ut + ct);
    let nv = match law {
        | Law::Euler { .. } => 5,
        | Law::ShallowWater { .. } => 4,
    };
    match solver {
        | RiemannSolver::Rusanov => {
            let s = (l[1].abs() + cl).max(r[1].abs() + cr);
            let mut f = [0.0; 5];
            for c in 0..nv {
                f[c] = f64::midpoint(fl[c], fr[c]) - 0.5 * s * (qr[c] - ql[c]);
            }
            f
        },
        | RiemannSolver::Hlle => {
            if sl >= 0.0 {
                return fl;
            }
            if sr <= 0.0 {
                return fr;
            }
            // Normalize the two wave speeds before forming their sum or product.
            let scale = sr.max(-sl);
            let left = -sl / scale;
            let right = sr / scale;
            let width = left + right;
            let diffusion = sr.min(-sl) / width;
            let mut f = [0.0; 5];
            for c in 0..nv {
                f[c] =
                    (right / width) * fl[c] + (left / width) * fr[c] - diffusion * (qr[c] - ql[c]);
            }
            f
        },
        | RiemannSolver::Hllc => {
            if sl >= 0.0 {
                return fl;
            }
            if sr <= 0.0 {
                return fr;
            }
            let (ul, ur) = (l[1], r[1]);
            let den = l[0] * (sl - ul) - r[0] * (sr - ur);
            let sm = match law {
                | Law::Euler { .. } => {
                    (r[4] - l[4] + l[0] * ul * (sl - ul) - r[0] * ur * (sr - ur)) / den
                },
                | Law::ShallowWater { .. } => {
                    (sl * r[0] * (ur - sr) - sr * l[0] * (ul - sl))
                        / (r[0] * (ur - sr) - l[0] * (ul - sl))
                },
            };
            let star = |w: &State, q: &State, s: f64| -> State {
                let k = w[0] * (s - w[1]) / (s - sm);
                match law {
                    | Law::Euler { .. } => {
                        [
                            k,
                            k * sm,
                            k * w[2],
                            k * w[3],
                            k * (q[4] / w[0] + (sm - w[1]) * (sm + w[4] / (w[0] * (s - w[1])))),
                        ]
                    },
                    | Law::ShallowWater { .. } => [k, k * sm, k * w[2], k * w[3], 0.0],
                }
            };
            let (f, w, q, s) = if sm >= 0.0 {
                (&fl, l, &ql, sl)
            } else {
                (&fr, r, &qr, sr)
            };
            let qs = star(w, q, s);
            let mut out = [0.0; 5];
            for c in 0..nv {
                out[c] = f[c] + s * (qs[c] - q[c]);
            }
            out
        },
        | RiemannSolver::Roe => {
            let mut diss = [0.0; 5];
            let fix = |lam: f64, delta: f64| {
                let a = lam.abs();
                if a < delta {
                    f64::midpoint(lam * lam / delta, delta)
                } else {
                    a
                }
            };
            match law {
                | Law::Euler { gamma } => {
                    let rt = sl_ * sr_;
                    let (hl, hr) = ((ql[4] + l[4]) / l[0], (qr[4] + r[4]) / r[0]);
                    let ht = (sl_ * hl + sr_ * hr) / (sl_ + sr_);
                    let vt = (sl_ * l[2] + sr_ * r[2]) / (sl_ + sr_);
                    let wt = (sl_ * l[3] + sr_ * r[3]) / (sl_ + sr_);
                    let c = ct;
                    let _ = gamma;
                    let (dr, dp, du, dv, dw) = (
                        r[0] - l[0],
                        r[4] - l[4],
                        r[1] - l[1],
                        r[2] - l[2],
                        r[3] - l[3],
                    );
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
                | Law::ShallowWater { .. } => {
                    let vt = (sl_ * l[2] + sr_ * r[2]) / (sl_ + sr_);
                    let wt = (sl_ * l[3] + sr_ * r[3]) / (sl_ + sr_);
                    let ht = f64::midpoint(l[0], r[0]);
                    let c = ct;
                    let (dh, dm) = (r[0] - l[0], r[0] * r[1] - l[0] * l[1]);
                    let a1 = ((ut + c) * dh - dm) / (2.0 * c);
                    let a3 = (dm - (ut - c) * dh) / (2.0 * c);
                    let delta = 0.1 * c;
                    let (l1, l2, l3) = (fix(ut - c, delta), ut.abs(), fix(ut + c, delta));
                    let r1 = [1.0, ut - c, vt, wt];
                    let r3 = [1.0, ut + c, vt, wt];
                    let (dv, dw) = (r[2] - l[2], r[3] - l[3]);
                    for k in 0..4 {
                        diss[k] = l1 * a1 * r1[k] + l3 * a3 * r3[k];
                    }
                    // The two transverse (shear) waves; in 1D and 2D the
                    // second one is identically zero.
                    diss[2] += l2 * ht * dv;
                    diss[3] += l2 * ht * dw;
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

fn limit(
    lim: SlopeLimiter,
    a: f64,
    b: f64,
) -> f64 {
    // a = backward difference, b = forward difference
    if a * b <= 0.0 {
        return 0.0;
    }
    match lim {
        | SlopeLimiter::Minmod => {
            if a.abs() < b.abs() {
                a
            } else {
                b
            }
        },
        | SlopeLimiter::VanLeer => 2.0 * a * b / (a + b),
        | SlopeLimiter::Mc => {
            let c = f64::midpoint(a, b);
            let m = (2.0 * a.abs()).min(2.0 * b.abs()).min(c.abs());
            m.copysign(a)
        },
        | SlopeLimiter::Superbee => {
            let s1 = if a.abs() < 2.0 * b.abs() {
                a
            } else {
                2.0 * b
            };
            let s2 = if 2.0 * a.abs() < b.abs() {
                2.0 * a
            } else {
                b
            };
            if s1.abs() > s2.abs() {
                s1
            } else {
                s2
            }
        },
    }
}

/// Fifth-order WENO value at the right face of the cell holding `v2`
/// (stencil `v0 .. v4`, centred on `v2`).
fn weno5(
    v0: f64,
    v1: f64,
    v2: f64,
    v3: f64,
    v4: f64,
) -> f64 {
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
pub(crate) fn face_states(
    recon: Reconstruction,
    line: &[State],
    k: usize,
    nv: usize,
) -> (State, State) {
    let mut l = line[k - 1];
    let mut r = line[k];
    match recon {
        | Reconstruction::FirstOrder => {},
        | Reconstruction::Muscl(lim) => {
            for c in 0..nv {
                let sl = limit(
                    lim,
                    line[k - 1][c] - line[k - 2][c],
                    line[k][c] - line[k - 1][c],
                );
                let sr = limit(
                    lim,
                    line[k][c] - line[k - 1][c],
                    line[k + 1][c] - line[k][c],
                );
                l[c] = line[k - 1][c] + 0.5 * sl;
                r[c] = line[k][c] - 0.5 * sr;
            }
        },
        | Reconstruction::Weno5 => {
            for c in 0..nv {
                l[c] = weno5(
                    line[k - 3][c],
                    line[k - 2][c],
                    line[k - 1][c],
                    line[k][c],
                    line[k + 1][c],
                );
                r[c] = weno5(
                    line[k + 2][c],
                    line[k + 1][c],
                    line[k][c],
                    line[k - 1][c],
                    line[k - 2][c],
                );
            }
        },
    }
    (l, r)
}

/// Ghost state of the linearised characteristic non-reflecting boundary.
fn nrbc_ghost(
    law: Law,
    q: State,
    w_inf: &State,
    axis: usize,
    max: bool,
) -> State {
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
    wg[0] = (wg[0]
        + (p_g
            - if un_inf >= 0.0 {
                wi[4]
            } else {
                w_inf[4]
            })
            / (c * c))
        .max(1e-8);
    wg[4] = p_g;
    wg[1 + axis] = sgn * un_g;
    to_cons(law, wg)
}

const GHOST: usize = 3;

/// One friction half step of the Manning source,
/// `q ← q / (1 + dt g n² |q| h^{-7/3})`, with the depth frozen at its current
/// value. This is the exact solution of `dq/dt = -g n² h^{-7/3} |q| q` over
/// `dt`: unconditionally stable for the source, monotone in `|q|` and
/// direction preserving.
fn manning_friction(
    law: Law,
    q: &mut State,
    dt: f64,
    n: f64,
) {
    let Law::ShallowWater { g } = law else { return };
    let h = q[0];
    let m = q[1].hypot(q[2]);
    if !(h > 0.0 && m > 0.0 && dt > 0.0 && n > 0.0 && h.is_finite() && m.is_finite()) {
        return;
    }
    let fac = 1.0 / (1.0 + dt * g * n * n * m / h.powf(7.0 / 3.0));
    q[1] *= fac;
    q[2] *= fac;
}

#[derive(Clone, Copy, Default)]
struct CellUpdate {
    delta: State,
    bed_slope: f64,
}

impl CellUpdate {
    fn subtract(
        self,
        out: &mut State,
        axis: usize,
        with_bed: bool,
    ) {
        for (r, delta) in out.iter_mut().zip(self.delta) {
            *r -= delta;
        }
        if with_bed {
            out[1 + axis] -= self.bed_slope;
        }
    }
}

struct BedFace {
    left: State,
    right: State,
    left_height: f64,
    left_bed: f64,
    right_height: f64,
    right_bed: f64,
}

#[derive(Default)]
struct ResidualWorkspace {
    padded: Vec<State>,
    updates: Vec<CellUpdate>,
    residual: Vec<State>,
}

const SWEEP_LINES_PER_CHUNK: usize = 8;
const SWEEP_MIN_CHUNKS: usize = 2;
const CELL_CHUNK: usize = 512;
pub(super) static LINE_STORAGE: crate::kernels::reuse::ScratchPool<State> =
    crate::kernels::reuse::ScratchPool::new(16, 256 * 1024);

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

    const fn idx(
        &self,
        i: [usize; 3],
    ) -> usize {
        i[0] + self.np[0] * (i[1] + self.np[1] * i[2])
    }

    const fn stride(
        &self,
        axis: usize,
    ) -> usize {
        [1, self.np[0], self.np[0] * self.np[1]][axis]
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
    gravity: [f64; 3],
    /// Bed elevation in the cell centres (shallow water).
    bed: Option<Vec<f64>>,
    /// Uniform Manning friction coefficient (zero = frictionless).
    manning: f64,
}

impl Stepper<'_> {
    fn fill_padded(
        &self,
        q: &[State],
        p: &mut [State],
    ) {
        let g = self.grid;
        let pad = &self.pad;
        let law = self.cfg.law;
        let width = g.n[0];
        for (row, states) in q.chunks_exact(width).enumerate() {
            let start = pad.idx([pad.g[0], row % g.n[1] + pad.g[1], row / g.n[1] + pad.g[2]]);
            let interior = &mut p[start..start + width];
            interior.copy_from_slice(states);
            if let Some(bed) = &self.bed {
                for (state, &height) in interior.iter_mut().zip(&bed[row * width..][..width]) {
                    state[4] = height;
                }
            }
        }
        for a in 0..g.dim {
            let n = g.n[a] as isize;
            let gh = GHOST as isize;
            let ranges: [std::ops::Range<usize>; 3] = std::array::from_fn(|b| {
                if b == a {
                    0..1
                } else if b < a || b >= g.dim {
                    0..pad.np[b]
                } else {
                    pad.g[b]..pad.g[b] + g.n[b]
                }
            });
            for i2 in ranges[2].clone() {
                for i1 in ranges[1].clone() {
                    for i0 in ranges[0].clone() {
                        for layer in 0..GHOST {
                            for max in [false, true] {
                                let mut dst = [i0, i1, i2];
                                let mut src = dst;
                                // ghost coordinate in interior index space
                                let gi: isize = if max {
                                    n + layer as isize
                                } else {
                                    -1 - layer as isize
                                };
                                dst[a] = (gi + gh) as usize;
                                let bc = self.cfg.bcs[2 * a + usize::from(max)];
                                let mut flip = false;
                                let s_i: isize = match bc {
                                    | GasBc::Transmissive => {
                                        if max {
                                            n - 1
                                        } else {
                                            0
                                        }
                                    },
                                    | GasBc::Reflective => {
                                        flip = true;
                                        if max {
                                            n - 1 - layer as isize
                                        } else {
                                            layer as isize
                                        }
                                    },
                                    | GasBc::Periodic => gi.rem_euclid(n),
                                    | GasBc::Fixed(w) => {
                                        let mut ghost = to_cons(law, w);
                                        if matches!(law, Law::ShallowWater { .. }) {
                                            // The bed (slot 4) is that of the nearest interior cell.
                                            let mut edge = dst;
                                            edge[a] = ((if max { n - 1 } else { 0 }) + gh) as usize;
                                            ghost[4] = p[pad.idx(edge)][4];
                                        }
                                        p[pad.idx(dst)] = ghost;
                                        continue;
                                    },
                                    | GasBc::NonReflecting(w_inf) => {
                                        let mut edge = dst;
                                        edge[a] = ((if max { n - 1 } else { 0 }) + gh) as usize;
                                        p[pad.idx(dst)] =
                                            nrbc_ghost(law, p[pad.idx(edge)], &w_inf, a, max);
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
    }

    fn max_rate(
        &self,
        q: &[State],
    ) -> f64 {
        let g = self.grid;
        q.iter()
            .map(|&c| {
                let w = to_prim(self.cfg.law, c);
                let cs = sound_speed(self.cfg.law, &w);
                (0..g.dim)
                    .map(|a| (w[1 + a].abs() + cs) / g.h(a))
                    .sum::<f64>()
            })
            .fold(0.0, f64::max)
    }

    /// Well-balanced pressure of face `f` of a direction-`d` line.
    ///
    /// Interior faces use [`hydrostatic_face_pressure`] with the two adjacent
    /// cell-centre pressures. A reflective face on an axis of gravity uses the
    /// hydrostatic continuation of the boundary cell across the wall, so the
    /// wall momentum flux is the pressure the cell-centred gravity source
    /// balances.
    fn well_balanced_face_pressure(
        &self,
        line: &[State],
        d: usize,
        n: usize,
        f: usize,
    ) -> f64 {
        let g = self.grid;
        let boundary = |k: usize, max: bool| -> f64 {
            let (p, rho) = (line[k][4], line[k][0]);
            if p > 0.0 && rho > 0.0 {
                let a = self.gravity[d] * g.h(d) / (p / rho);
                if max {
                    hydrostatic_face_pressure(p, p * a.exp())
                } else {
                    hydrostatic_face_pressure(p * (-a).exp(), p)
                }
            } else {
                hydrostatic_face_pressure(line[GHOST + f - 1][4], line[GHOST + f][4])
            }
        };
        if f == 0 && matches!(self.cfg.bcs[2 * d], GasBc::Reflective) {
            boundary(GHOST, false)
        } else if f == n && matches!(self.cfg.bcs[2 * d + 1], GasBc::Reflective) {
            boundary(GHOST + n - 1, true)
        } else {
            hydrostatic_face_pressure(line[GHOST + f - 1][4], line[GHOST + f][4])
        }
    }

    fn regular_face(
        &self,
        line: &[State],
        d: usize,
        f: usize,
        well_balanced: bool,
    ) -> State {
        let law = self.cfg.law;
        let (mut l, mut r) = face_states(self.cfg.recon, line, GHOST + f, self.nv);
        if !admissible(law, &l) || !admissible(law, &r) {
            l = line[GHOST + f - 1];
            r = line[GHOST + f];
        }
        swap_dir(&mut l, d);
        swap_dir(&mut r, d);
        if well_balanced {
            // Both Riemann states carry the face pressure balanced by the cell source.
            let pressure = self.well_balanced_face_pressure(line, d, self.grid.n[d], f);
            l[4] = pressure;
            r[4] = pressure;
        }
        let mut flux = riemann_flux(law, self.cfg.riemann, &l, &r);
        swap_dir(&mut flux, d);
        flux
    }

    fn bed_face(
        &self,
        line: &[State],
        d: usize,
        f: usize,
        gravity: f64,
    ) -> BedFace {
        // Reconstruct (eta, u, v, -, b), then use the higher bed on both sides.
        let law = self.cfg.law;
        let (mut l, mut r) = face_states(self.cfg.recon, line, GHOST + f, 5);
        let low = |s: &State| s[0] - s[4] < 0.0 || !s[0].is_finite();
        if low(&l) || low(&r) {
            l = line[GHOST + f - 1];
            r = line[GHOST + f];
        }
        let (h_l, h_r) = ((l[0] - l[4]).max(0.0), (r[0] - r[4]).max(0.0));
        let b_star = l[4].max(r[4]);
        let (hs_l, hs_r) = ((l[0] - b_star).max(0.0), (r[0] - b_star).max(0.0));
        let mut wl = [hs_l, l[1], l[2], l[3], 0.0];
        let mut wr = [hs_r, r[1], r[2], r[3], 0.0];
        swap_dir(&mut wl, d);
        swap_dir(&mut wr, d);
        let mut left = if hs_l + hs_r > 1e-12 && admissible(law, &wl) && admissible(law, &wr) {
            riemann_flux(law, self.cfg.riemann, &wl, &wr)
        } else if hs_l > 1e-12 && admissible(law, &wl) {
            flux_x(law, &wl)
        } else if hs_r > 1e-12 && admissible(law, &wr) {
            flux_x(law, &wr)
        } else {
            [0.0; 5]
        };
        swap_dir(&mut left, d);
        let mut right = left;
        left[1 + d] += 0.5 * gravity * (h_l * h_l - hs_l * hs_l);
        right[1 + d] += 0.5 * gravity * (h_r * h_r - hs_r * hs_r);
        BedFace {
            left,
            right,
            left_height: h_l,
            left_bed: l[4],
            right_height: h_r,
            right_bed: r[4],
        }
    }

    fn line_updates(
        &self,
        d: usize,
        line: &[State],
        mut store: impl FnMut(usize, CellUpdate),
    ) {
        let n = self.grid.n[d];
        let spacing = self.grid.h(d);
        if let (Some(_), Law::ShallowWater { g }) = (&self.bed, self.cfg.law) {
            let mut lower = self.bed_face(line, d, 0, g);
            for cell in 0..n {
                let upper = self.bed_face(line, d, cell + 1, g);
                let delta = std::array::from_fn(|v| (upper.left[v] - lower.right[v]) / spacing);
                let bed_slope = g
                    * 0.5
                    * (lower.right_height + upper.left_height)
                    * (upper.left_bed - lower.right_bed)
                    / spacing;
                store(cell, CellUpdate { delta, bed_slope });
                lower = upper;
            }
        } else {
            let well_balanced = matches!(self.cfg.law, Law::Euler { .. }) && self.gravity[d] != 0.0;
            let mut lower = self.regular_face(line, d, 0, well_balanced);
            for cell in 0..n {
                let upper = self.regular_face(line, d, cell + 1, well_balanced);
                let delta = std::array::from_fn(|v| (upper[v] - lower[v]) / spacing);
                store(
                    cell,
                    CellUpdate {
                        delta,
                        bed_slope: 0.0,
                    },
                );
                lower = upper;
            }
        }
    }

    fn residual(
        &self,
        q: &[State],
        work: &mut ResidualWorkspace,
    ) {
        use crate::kernels::par::for_chunks_mut;

        let grid = self.grid;
        let pad = &self.pad;
        let law = self.cfg.law;
        let with_bed = self.bed.is_some();
        work.padded.resize(pad.len(), [0.0; 5]);
        self.fill_padded(q, &mut work.padded);
        for_chunks_mut(&mut work.padded, CELL_CHUNK, 16, |_, cells| {
            for state in cells {
                let bed = state[4];
                *state = to_prim(law, *state);
                if with_bed {
                    state[0] += bed;
                    state[4] = bed;
                }
            }
        });
        work.residual.resize(q.len(), [0.0; 5]);
        work.residual.fill([0.0; 5]);

        let [nx, ny, _] = grid.n;
        let padded = &work.padded;
        // The first-axis lines are contiguous in both the padded states and the residual.
        for_chunks_mut(
            &mut work.residual,
            nx * SWEEP_LINES_PER_CHUNK,
            SWEEP_MIN_CHUNKS,
            |block, cells| {
                for (line, row) in cells.chunks_exact_mut(nx).enumerate() {
                    let flat = block * SWEEP_LINES_PER_CHUNK + line;
                    let start = pad.idx([0, flat % ny + pad.g[1], flat / ny + pad.g[2]]);
                    self.line_updates(0, &padded[start..start + pad.np[0]], |i, update| {
                        update.subtract(&mut row[i], 0, with_bed);
                    });
                }
            },
        );

        for d in 1..grid.dim {
            let n = grid.n[d];
            let outer_axis = 3 - d;
            let stride = pad.stride(d);
            let line_len = pad.np[d];
            work.updates.resize(q.len(), CellUpdate::default());
            // Pack one strided line per work block, and write its updates in line order.
            for_chunks_mut(
                &mut work.updates,
                n * SWEEP_LINES_PER_CHUNK,
                SWEEP_MIN_CHUNKS,
                |block, updates| {
                    let mut packed = LINE_STORAGE.checkout_for(line_len);
                    packed.resize(line_len, [0.0; 5]);
                    for (line, row) in updates.chunks_exact_mut(n).enumerate() {
                        let flat = block * SWEEP_LINES_PER_CHUNK + line;
                        let mut i = [0; 3];
                        i[0] = flat % nx + pad.g[0];
                        i[outer_axis] = flat / nx + pad.g[outer_axis];
                        let start = pad.idx(i);
                        for (k, value) in packed.iter_mut().enumerate() {
                            *value = padded[start + k * stride];
                        }
                        self.line_updates(d, &packed, |i, update| row[i] = update);
                    }
                },
            );
            // Each chunk owns complete x rows. Every cell sees x, then y, then z.
            let updates = &work.updates;
            for_chunks_mut(
                &mut work.residual,
                nx * SWEEP_LINES_PER_CHUNK,
                SWEEP_MIN_CHUNKS,
                |block, cells| {
                    for (line, row) in cells.chunks_exact_mut(nx).enumerate() {
                        let flat = block * SWEEP_LINES_PER_CHUNK + line;
                        let (y, z) = (flat % ny, flat / ny);
                        let start = if d == 1 {
                            z * nx * n + y
                        } else {
                            y * nx * n + z
                        };
                        for (i, out) in row.iter_mut().enumerate() {
                            updates[start + i * n].subtract(out, d, with_bed);
                        }
                    }
                },
            );
        }
        if matches!(law, Law::Euler { .. }) && self.gravity.iter().any(|&a| a != 0.0) {
            for (res, &state) in work.residual.iter_mut().zip(q) {
                let w = to_prim(law, state);
                for a in 0..grid.dim {
                    res[1 + a] += w[0] * self.gravity[a];
                    res[4] += w[0] * w[1 + a] * self.gravity[a];
                }
            }
        }
    }
}

/// Runs the solver from primitive initial data to `t_end`.
///
/// # Errors
/// Returns an error for an invalid grid, an invalid time or CFL number,
/// inconsistent periodic faces or inadmissible initial data.
pub fn solve(
    grid: &FvGrid,
    cfg: &FvConfig,
    init: &dyn Fn([f64; 3]) -> State,
    t_end: f64,
) -> Result<FvSolution, String> {
    solve_with_sources(grid, cfg, init, &Sources::default(), t_end)
}

/// Runs the solver with the source terms of [`Sources`].
///
/// Gravity (Euler) uses the well-balanced interface pressure described in the
/// module documentation; a bed (shallow water) uses the hydrostatic
/// reconstruction, which is well-balanced for the lake at rest.
///
/// # Errors
/// As [`solve`], and for gravity with shallow water, a bed or Manning friction
/// with Euler, or a bed or Manning friction with shallow water in 3D.
pub fn solve_with_sources(
    grid: &FvGrid,
    cfg: &FvConfig,
    init: &dyn Fn([f64; 3]) -> State,
    sources: &Sources<'_>,
    t_end: f64,
) -> Result<FvSolution, String> {
    grid.check()?;
    if !t_end.is_finite() || t_end < 0.0 {
        return Err("final time must be finite and nonnegative".into());
    }
    if !(cfg.cfl > 0.0 && cfg.cfl < 1.0) {
        return Err("cfl must be in (0, 1)".into());
    }
    let valid_law = match cfg.law {
        | Law::Euler { gamma } => gamma.is_finite() && gamma > 1.0,
        | Law::ShallowWater { g } => g.is_finite() && g >= 0.0,
    };
    if !valid_law || !sources.gravity.iter().all(|g| g.is_finite()) {
        return Err("invalid conservation law or gravity".into());
    }
    match cfg.law {
        | Law::ShallowWater { .. } if sources.gravity.iter().any(|&a| a != 0.0) => {
            return Err("shallow water has its own gravity: use the bed".into());
        },
        | Law::Euler { .. } if sources.bed.is_some() => {
            return Err("a bed is a shallow-water source".into());
        },
        | _ => {},
    }
    if !(sources.manning.is_finite() && sources.manning >= 0.0) {
        return Err("the Manning coefficient must be finite and non-negative".into());
    }
    if sources.manning > 0.0 && matches!(cfg.law, Law::Euler { .. }) {
        return Err("Manning friction is a shallow-water source".into());
    }
    if matches!(cfg.law, Law::ShallowWater { .. })
        && cfg.bcs.iter().any(|b| matches!(b, GasBc::NonReflecting(_)))
    {
        return Err("non-reflecting boundaries are available for the Euler equations only".into());
    }
    if matches!(cfg.law, Law::ShallowWater { .. }) && grid.dim == 3 {
        if sources.bed.is_some() {
            return Err(
                "a bed is a planar shallow-water notion: not available for the 3D barotropic system"
                    .into(),
            );
        }
        if sources.manning > 0.0 {
            return Err(
                "Manning bed friction is a planar shallow-water notion: not available for the 3D \
                 barotropic system"
                    .into(),
            );
        }
    }
    for a in 0..grid.dim {
        if (cfg.bcs[2 * a] == GasBc::Periodic) != (cfg.bcs[2 * a + 1] == GasBc::Periodic) {
            return Err("periodic must be set on both faces of an axis".into());
        }
    }
    for boundary in &cfg.bcs[..2 * grid.dim] {
        if let GasBc::Fixed(state) | GasBc::NonReflecting(state) = boundary
            && !finite_primitive(cfg.law, state)
        {
            return Err("boundary state is not admissible".into());
        }
    }
    let nv = match cfg.law {
        | Law::Euler { .. } => 5,
        | Law::ShallowWater { .. } => 1 + grid.dim,
    };
    let mut q: Vec<State> = Vec::with_capacity(grid.len());
    for p in 0..grid.len() {
        let mut w = init(grid.centre(p));
        if grid.dim < 3 && matches!(cfg.law, Law::ShallowWater { .. }) {
            // Slot 3 is the vertical-axis velocity of the 3D system only.
            w[3] = 0.0;
        }
        if !finite_primitive(cfg.law, &w) {
            return Err("initial state is not admissible".into());
        }
        let state = to_cons(cfg.law, w);
        if state.iter().any(|value| !value.is_finite()) {
            return Err("initial conserved state is not representable".into());
        }
        q.push(state);
    }
    let bed = sources
        .bed
        .map(|sample| {
            (0..grid.len())
                .map(|p| {
                    let height = sample(grid.centre(p));
                    if height.is_finite() {
                        Ok(height)
                    } else {
                        Err("bed elevation must be finite")
                    }
                })
                .collect::<Result<Vec<f64>, _>>()
        })
        .transpose()?;
    let st = Stepper {
        grid,
        cfg,
        pad: Padded::new(grid),
        nv,
        gravity: sources.gravity,
        bed,
        manning: sources.manning,
    };
    let mut t = 0.0;
    let mut steps = 0;
    let mut work = ResidualWorkspace::default();
    let mut stage = Vec::new();
    while t < t_end * (1.0 - 1e-14) {
        crate::sim::models_geodesic_relativity::native_charge();
        let rate = st.max_rate(&q);
        if !rate.is_finite() || rate <= 0.0 {
            return Err("non-finite wave speed".into());
        }
        let dt = (cfg.cfl / rate).min(t_end - t);
        let next_time = t + dt;
        if next_time <= t {
            return Err("finite-volume time step cannot advance time".into());
        }
        // Strang splitting of the Manning friction: the source is integrated
        // exactly over a half step before and after the hyperbolic update.
        let friction = st.manning > 0.0;
        if friction {
            for state in &mut q {
                manning_friction(cfg.law, state, 0.5 * dt, st.manning);
            }
        }
        // Reuse the same intermediate state and residual storage through SSP-RK3.
        st.residual(&q, &mut work);
        stage.resize(q.len(), [0.0; 5]);
        for ((out, &start), &rate) in stage.iter_mut().zip(&q).zip(&work.residual) {
            *out = std::array::from_fn(|c| start[c] + rate[c] * dt);
        }
        st.residual(&stage, &mut work);
        for ((out, &start), &rate) in stage.iter_mut().zip(&q).zip(&work.residual) {
            *out = std::array::from_fn(|c| 0.75 * start[c] + 0.25 * (out[c] + rate[c] * dt));
        }
        st.residual(&stage, &mut work);
        for ((out, &intermediate), &rate) in q.iter_mut().zip(&stage).zip(&work.residual) {
            *out = std::array::from_fn(|c| {
                out[c] / 3.0 + 2.0 * (intermediate[c] + rate[c] * dt) / 3.0
            });
        }
        if friction {
            for state in &mut q {
                manning_friction(cfg.law, state, 0.5 * dt, st.manning);
            }
        }
        t = next_time;
        steps += 1;
        if steps > 2_000_000 {
            return Err("too many steps".into());
        }
    }
    Ok(FvSolution {
        grid: *grid,
        q,
        t,
        steps,
        law: cfg.law,
    })
}
