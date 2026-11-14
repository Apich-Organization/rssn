//! # Scalar wave equation / linear acoustics with PML in 1D, 2D and 3D
//!
//! The scalar wave equation `p_tt = c² Δp` is solved in its first-order
//! (linear acoustics) form on a staggered grid, with `p` in the cells and the
//! velocity components on the faces. The medium is a wave speed `c(x)` and a
//! mass density `ρ(x)`; the scheme integrates
//!
//! ```text
//! p_t = -K ∇·v,   v_t = -ρ⁻¹ ∇p,   K = ρ c²,
//! ```
//!
//! which is the conservative form `∂ₜ(ρ⁻¹∂ₜp) = ∇·(ρ⁻¹∇p)` of the
//! variable-medium equation. It keeps the second order of the uniform scheme
//! across a jump in the medium: the density on a face is the harmonic mean of
//! the two cells it separates (the mean of the compliances, which is what
//! continuity of the normal velocity asks for) and the bulk modulus `ρc²` is
//! constant within a cell. The obvious `c(x)² ∇²p` is only first order at such
//! a jump. For constant `c` and `ρ = 1` every term collapses to the plain
//! `p_t = -c²∇·v`, `v_t = -∇p` scheme.
//!
//! [`AcousticsConfig::c`] and [`AcousticsConfig::rho`] describe a uniform
//! medium; [`Acoustics::set_medium`] replaces it with a spatially varying one
//! and re-derives [`Acoustics::dt`] from the fastest effective speed. The
//! invariant of the scheme is the discrete energy `½ Σ (p²/K + ρ v²) ΔV`
//! ([`Acoustics::energy`]); the naive `½ Σ (p²/c² + |v|²)` is not conserved
//! once `ρ` varies.
//!
//! Boundary faces use [`BoundaryCondition`]:
//!
//! | Face condition | Meaning |
//! |---|---|
//! | `Neumann` / `Symmetry` | rigid wall, `v·n = 0` |
//! | `Dirichlet(g)` | prescribed pressure `p = g(t, x)` at the face (`g = 0` is a pressure-release surface) |
//! | `Periodic` | wrap-around (both faces of the axis) |
//! | `Absorbing` | convolutional PML of [`AcousticsConfig::pml_cells`] cells ending in a rigid wall |
//!
//! The PML is the stretched-coordinate convolutional PML (`κ = 1`, `α = 0`)
//! with polynomial conductivity grading and one auxiliary variable per
//! derivative; because it acts on the first-order system it is exactly the
//! PML of the scalar wave equation for `p`.
//!
//! The Dirichlet datum is applied as a ghost value: the face sits half a cell
//! from the first cell centre, so the normal pressure gradient there is
//! `(p_cell - g) / (h/2)` on a low face and `(g - p_cell) / (h/2)` on a high
//! one, the outward normal flipping the sign between the two. A
//! time-dependent `g` is evaluated at the current time `t`; the leapfrog
//! scheme is second order, so the half-step offset of the velocity field does
//! not degrade the boundary condition.
//!
//! A Dirichlet face is a **pressure-release** surface and therefore a
//! reflector, not an absorber: it pins `p` and leaves `v` free, so a
//! transient keeps bouncing between it and the other faces instead of
//! radiating out. A boundary that absorbs needs [`BoundaryCondition::Absorbing`]
//! (the PML), and a driven problem with a steady state needs one of the two.
//!
//! Limitations: the medium is piecewise constant on the cells, so an
//! interface is only represented to second order when it lies on a face and is
//! smeared over one cell otherwise; `c` and `ρ` are time independent; there is
//! no attenuation (a complex, frequency-dependent bulk modulus) and no flow.
//! Robin faces are refused.

use crate::sim::boundary::BcValue;
use crate::sim::boundary::BoundaryCondition;
use crate::sim::boundary::Face;
use crate::sim::boundary::FaceBcs;

type SourceFn = Box<dyn Fn(f64) -> f64 + Send + Sync>;

/// Configuration of the acoustics solver.
pub struct AcousticsConfig {
    /// Cells per axis (inactive axes: 1).
    pub n: [usize; 3],
    /// Active axes.
    pub dim: usize,
    /// Cell size per axis.
    pub h: [f64; 3],
    /// Wave speed of the uniform medium.
    pub c: f64,
    /// Mass density of the uniform medium (`1` for a velocity-only medium).
    pub rho: f64,
    /// Courant number (`<= 1`).
    pub cfl: f64,
    /// Face conditions.
    pub bcs: FaceBcs,
    /// PML thickness in cells.
    pub pml_cells: usize,
    /// PML grading exponent.
    pub pml_order: f64,
}

impl AcousticsConfig {
    /// Rigid-wall configuration with unit density.
    ///
    /// # Errors
    /// For a bad dimension or sizes.
    pub fn new(
        n: &[usize],
        h: &[f64],
        c: f64,
    ) -> Result<Self, String> {
        let dim = n.len();
        if !(1..=3).contains(&dim) || h.len() != dim || n.iter().any(|&m| m < 4) || c <= 0.0 {
            return Err("invalid acoustics grid".into());
        }
        let mut cfg = Self {
            n: [1; 3],
            dim,
            h: [1.0; 3],
            c,
            rho: 1.0,
            cfl: 0.9,
            bcs: FaceBcs::uniform(&BoundaryCondition::neumann(0.0)),
            pml_cells: 10,
            pml_order: 3.0,
        };
        cfg.n[..dim].copy_from_slice(&n[..dim]);
        cfg.h[..dim].copy_from_slice(&h[..dim]);
        Ok(cfg)
    }
}

/// Staggered-grid acoustics solver with optional PML.
pub struct Acoustics {
    cfg: AcousticsConfig,
    /// Time step.
    pub dt: f64,
    /// Current time.
    pub t: f64,
    /// Pressure in the cells.
    pub p: Vec<f64>,
    /// Face velocities (`v[a]` has `n[a] + 1` entries along axis `a`).
    pub v: [Vec<f64>; 3],
    periodic: [bool; 3],
    /// Bulk modulus `ρc²` at the cell centres.
    bulk: Vec<f64>,
    /// Density at the faces, laid out like [`Acoustics::v`].
    rho_f: [Vec<f64>; 3],
    /// Fastest effective speed; sets [`Acoustics::dt`] and the PML strength.
    c_ref: f64,
    // (b, a) CPML coefficients per axis at cell centres and at faces
    pml_c: [Vec<(f64, f64)>; 3],
    pml_f: [Vec<(f64, f64)>; 3],
    psi_p: [Vec<f64>; 3],
    psi_v: [Vec<f64>; 3],
    use_pml: bool,
    /// Faces carrying a PML, per axis and face.
    absorbing: [[bool; 2]; 3],
    /// Dirichlet datum per axis and face; `Some(g)` replaces the rigid wall.
    dirichlet: [[Option<BcValue>; 2]; 3],
    sources: Vec<(usize, SourceFn)>,
}

impl Acoustics {
    /// Creates the solver.
    ///
    /// # Errors
    /// For unsupported boundary conditions or a PML that is too thick.
    pub fn new(cfg: AcousticsConfig) -> Result<Self, String> {
        if !(cfg.cfl > 0.0 && cfg.cfl <= 1.0) {
            return Err("cfl must be in (0, 1]".into());
        }
        if !(cfg.rho > 0.0 && cfg.rho.is_finite()) {
            return Err("density must be positive".into());
        }
        let inv: f64 = (0..cfg.dim).map(|a| 1.0 / (cfg.h[a] * cfg.h[a])).sum();
        let dt = cfg.cfl / (cfg.c * inv.sqrt());
        let mut periodic = [false; 3];
        let mut dirichlet: [[Option<BcValue>; 2]; 3] = std::array::from_fn(|_| [None, None]);
        let mut absorbing = [[false; 2]; 3];
        let mut use_pml = false;
        let pml_c: [Vec<(f64, f64)>; 3] = std::array::from_fn(|a| vec![(1.0, 0.0); cfg.n[a]]);
        let pml_f: [Vec<(f64, f64)>; 3] = std::array::from_fn(|a| vec![(1.0, 0.0); cfg.n[a] + 1]);
        for a in 0..cfg.dim {
            periodic[a] = cfg.bcs.axis_periodic(a)?;
            for (m, max) in [false, true].into_iter().enumerate() {
                match cfg.bcs.get(Face::of_axis(a, max)) {
                    | BoundaryCondition::Neumann(_)
                    | BoundaryCondition::Symmetry
                    | BoundaryCondition::Periodic => {},
                    | BoundaryCondition::Dirichlet(g) => {
                        dirichlet[a][m] = Some(g.clone());
                    },
                    | BoundaryCondition::Absorbing => {
                        if 2 * cfg.pml_cells >= cfg.n[a] {
                            return Err("PML too thick for the grid".into());
                        }
                        use_pml = true;
                        absorbing[a][m] = true;
                    },
                    | BoundaryCondition::Robin { .. } => {
                        return Err("Robin faces are not supported".into());
                    },
                }
            }
        }
        let nn: usize = cfg.n.iter().product();
        let v: [Vec<f64>; 3] = std::array::from_fn(|a| {
            if a < cfg.dim {
                let mut e = cfg.n;
                e[a] += 1;
                vec![0.0; e.iter().product()]
            } else {
                Vec::new()
            }
        });
        let rho_f: [Vec<f64>; 3] = std::array::from_fn(|a| {
            if a < cfg.dim {
                vec![cfg.rho; v[a].len()]
            } else {
                Vec::new()
            }
        });
        let psi_p: [Vec<f64>; 3] = std::array::from_fn(|a| {
            if use_pml && a < cfg.dim {
                vec![0.0; nn]
            } else {
                Vec::new()
            }
        });
        let psi_v: [Vec<f64>; 3] = std::array::from_fn(|a| {
            if use_pml && a < cfg.dim {
                vec![0.0; v[a].len()]
            } else {
                Vec::new()
            }
        });
        let bulk = vec![cfg.rho * cfg.c * cfg.c; nn];
        let c_ref = cfg.c;
        let mut solver = Self {
            cfg,
            dt,
            t: 0.0,
            p: vec![0.0; nn],
            v,
            periodic,
            bulk,
            rho_f,
            c_ref,
            pml_c,
            pml_f,
            psi_p,
            psi_v,
            use_pml,
            absorbing,
            dirichlet,
            sources: Vec::new(),
        };
        if use_pml {
            solver.set_pml();
        }
        Ok(solver)
    }

    /// Builds the convolutional-PML coefficients from `dt`, `c_ref` and the
    /// absorbing faces.
    fn set_pml(&mut self) {
        let d = self.cfg.pml_cells;
        let dt = self.dt;
        let order = self.cfg.pml_order;
        for a in 0..self.cfg.dim {
            let (n, h) = (self.cfg.n[a], self.cfg.h[a]);
            let smax = 0.8 * (order + 1.0) * self.c_ref / h;
            for (m, max) in [false, true].into_iter().enumerate() {
                if !self.absorbing[a][m] {
                    continue;
                }
                let prof = |pos: f64| -> Option<(f64, f64)> {
                    let depth = if max {
                        pos - (n - d) as f64
                    } else {
                        d as f64 - pos
                    };
                    (depth > 0.0).then(|| {
                        let sigma = smax * (depth / d as f64).powf(order);
                        let b = (-sigma * dt).exp();
                        (b, b - 1.0)
                    })
                };
                for i in 0..n {
                    if let Some(c) = prof(i as f64 + 0.5) {
                        self.pml_c[a][i] = c;
                    }
                }
                for i in 0..=n {
                    if let Some(c) = prof(i as f64) {
                        self.pml_f[a][i] = c;
                    }
                }
            }
        }
    }

    /// Sets a spatially varying medium from `x ↦ (c, ρ)` sampled at the cell
    /// centres, and re-derives [`Acoustics::dt`] and the PML from the result.
    ///
    /// The density on a face is the harmonic mean of the two cells it
    /// separates and the bulk modulus of a cell is `ρc²`, which is what keeps
    /// the scheme second order across an interface (an interface is exact when
    /// it lies on a face). Call this before stepping: the fields, the clock and
    /// the sources are left alone.
    ///
    /// # Errors
    /// A non-positive or non-finite `c` or `ρ` at a cell centre.
    pub fn set_medium(
        &mut self,
        f: impl Fn([f64; 3]) -> (f64, f64),
    ) -> Result<(), String> {
        let (n, dim, h) = (self.cfg.n, self.cfg.dim, self.cfg.h);
        let cells = self.p.len();
        let mut bulk = vec![0.0; cells];
        let mut speed = vec![0.0; cells];
        let mut rho = vec![0.0; cells];
        let mut rho_min = f64::INFINITY;
        for q in 0..cells {
            let ijk = [q % n[0], (q / n[0]) % n[1], q / (n[0] * n[1])];
            let mut x = [0.0; 3];
            for a in 0..dim {
                x[a] = (ijk[a] as f64 + 0.5) * h[a];
            }
            let (c, r) = f(x);
            if !(c > 0.0 && r > 0.0 && c.is_finite() && r.is_finite()) {
                return Err("invalid medium".into());
            }
            bulk[q] = r * c * c;
            speed[q] = c;
            rho[q] = r;
            rho_min = rho_min.min(r);
        }
        self.bulk = bulk;
        for a in 0..dim {
            let e = self.face_dims(a);
            let total = e.iter().product::<usize>();
            for q in 0..total {
                let ijk = [q % e[0], (q / e[0]) % e[1], q / (e[0] * e[1])];
                let i = ijk[a];
                let at = |ii: usize| {
                    let mut c = ijk;
                    c[a] = ii;
                    c[0] + n[0] * (c[1] + n[1] * c[2])
                };
                // Harmonic mean of the two cells: the mean of the compliances.
                // A boundary face has one cell and takes its value.
                let compliance = if i == 0 {
                    1.0 / rho[at(0)]
                } else if i == n[a] {
                    1.0 / rho[at(n[a] - 1)]
                } else {
                    f64::midpoint(1.0 / rho[at(i - 1)], 1.0 / rho[at(i)])
                };
                self.rho_f[a][q] = 1.0 / compliance;
            }
        }
        // The CFL bound of the coupled system is `dt <= 2 / sqrt(Σ 4 K/(ρ h²))`,
        // so the reference speed is the largest `sqrt(K/ρ_min)` over the cells.
        let mut c_ref: f64 = 0.0;
        for q in 0..cells {
            c_ref = c_ref.max(speed[q] * (rho[q] / rho_min).sqrt());
        }
        let inv: f64 = (0..dim).map(|a| 1.0 / (h[a] * h[a])).sum();
        self.c_ref = c_ref;
        self.dt = self.cfg.cfl / (c_ref * inv.sqrt());
        if self.use_pml {
            self.set_pml();
        }
        Ok(())
    }

    /// Cells per axis.
    #[must_use]
    pub const fn n(&self) -> [usize; 3] {
        self.cfg.n
    }

    /// Flat cell index.
    #[must_use]
    pub const fn idx(
        &self,
        i: usize,
        j: usize,
        k: usize,
    ) -> usize {
        i + self.cfg.n[0] * (j + self.cfg.n[1] * k)
    }

    /// Sets the pressure from a function of the cell-centre position (velocity zero).
    pub fn set_pressure(
        &mut self,
        f: impl Fn([f64; 3]) -> f64,
    ) {
        let n = self.cfg.n;
        for p in 0..self.p.len() {
            let ijk = [p % n[0], (p / n[0]) % n[1], p / (n[0] * n[1])];
            let mut x = [0.0; 3];
            for a in 0..self.cfg.dim {
                x[a] = (ijk[a] as f64 + 0.5) * self.cfg.h[a];
            }
            self.p[p] = f(x);
        }
    }

    /// Adds a soft source `s(t)` to the pressure equation of a cell.
    pub fn add_source(
        &mut self,
        cell: usize,
        s: impl Fn(f64) -> f64 + Send + Sync + 'static,
    ) {
        self.sources.push((cell, Box::new(s)));
    }

    /// Acoustic energy `½ Σ (p²/K + ρ|v|²) ΔV`, the invariant of the
    /// conservative scheme (`K = ρc²` in the cells, `ρ` the face density).
    #[must_use]
    pub fn energy(&self) -> f64 {
        let dv: f64 = (0..self.cfg.dim).map(|a| self.cfg.h[a]).product();
        let mut s: f64 = self.p.iter().zip(&self.bulk).map(|(x, k)| x * x / k).sum();
        for a in 0..self.cfg.dim {
            let mut e = self.cfg.n;
            e[a] += 1;
            for (q, x) in self.v[a].iter().enumerate() {
                let ia = (q / e[..a].iter().product::<usize>().max(1)) % e[a];
                let w = if ia == 0 || ia == self.cfg.n[a] {
                    0.5
                } else {
                    1.0
                };
                s += w * self.rho_f[a][q] * x * x;
            }
        }
        0.5 * s * dv
    }

    const fn face_dims(
        &self,
        a: usize,
    ) -> [usize; 3] {
        let mut e = self.cfg.n;
        e[a] += 1;
        e
    }

    /// Velocity half of the leapfrog step with step size `dt`.
    fn update_v(
        &mut self,
        dt: f64,
    ) {
        let cfg = &self.cfg;
        let (n, dim) = (cfg.n, cfg.dim);
        // velocity update at faces
        let t = self.t;
        for a in 0..dim {
            let e = self.face_dims(a);
            let total = e.iter().product::<usize>();
            for q in 0..total {
                let ijk = [q % e[0], (q / e[0]) % e[1], q / (e[0] * e[1])];
                let i = ijk[a];
                let h = cfg.h[a];
                let pidx = |ii: usize| {
                    let mut c = ijk;
                    c[a] = ii;
                    c[0] + n[0] * (c[1] + n[1] * c[2])
                };
                // Where a face sits, for a datum that depends on position.
                let face_at = |index: usize| -> [f64; 3] {
                    let mut x = [0.0; 3];
                    for b in 0..dim {
                        x[b] = if b == a {
                            index as f64 * cfg.h[b]
                        } else {
                            (ijk[b] as f64 + 0.5) * cfg.h[b]
                        };
                    }
                    x
                };
                let grad = if i == 0 {
                    if self.periodic[a] {
                        (self.p[pidx(0)] - self.p[pidx(n[a] - 1)]) / h
                    } else if let Some(g) = &self.dirichlet[a][0] {
                        // ghost value: the face is half a cell from the centre
                        (self.p[pidx(0)] - g.eval(t, face_at(0))) / (0.5 * h)
                    } else {
                        continue; // rigid wall: v stays zero
                    }
                } else if i == n[a] {
                    if self.periodic[a] {
                        // alias of face 0, set after the sweep
                        continue;
                    } else if let Some(g) = &self.dirichlet[a][1] {
                        // The outward normal is `-x` on this face, so the sign
                        // of the gradient is the mirror of the low face's.
                        (g.eval(t, face_at(n[a])) - self.p[pidx(n[a] - 1)]) / (0.5 * h)
                    } else {
                        continue;
                    }
                } else {
                    (self.p[pidx(i)] - self.p[pidx(i - 1)]) / h
                };
                let mut g = grad;
                if self.use_pml {
                    let (b, aa) = self.pml_f[a][i.min(n[a])];
                    let psi = &mut self.psi_v[a][q];
                    *psi = b * *psi + aa * g;
                    g += *psi;
                }
                let rho = self.rho_f[a][q];
                self.v[a][q] -= dt * g / rho;
            }
            if self.periodic[a] {
                for q in 0..total {
                    let ijk = [q % e[0], (q / e[0]) % e[1], q / (e[0] * e[1])];
                    if ijk[a] == n[a] {
                        let mut c = ijk;
                        c[a] = 0;
                        let src = c[0] + e[0] * (c[1] + e[1] * c[2]);
                        self.v[a][q] = self.v[a][src];
                    }
                }
            }
        }
    }

    fn update_p(
        &mut self,
        dt: f64,
    ) {
        let cfg = &self.cfg;
        let (n, dim) = (cfg.n, cfg.dim);
        // pressure update at cells
        for p in 0..self.p.len() {
            let ijk = [p % n[0], (p / n[0]) % n[1], p / (n[0] * n[1])];
            let mut div = 0.0;
            for a in 0..dim {
                let e = self.face_dims(a);
                let h = cfg.h[a];
                let mut hi = ijk;
                hi[a] += 1;
                let d = (self.v[a][hi[0] + e[0] * (hi[1] + e[1] * hi[2])]
                    - self.v[a][ijk[0] + e[0] * (ijk[1] + e[1] * ijk[2])])
                    / h;
                let mut dd = d;
                if self.use_pml {
                    let (b, aa) = self.pml_c[a][ijk[a]];
                    let psi = &mut self.psi_p[a][p];
                    *psi = b * *psi + aa * d;
                    dd += *psi;
                }
                div += dd;
            }
            let k = self.bulk[p];
            self.p[p] -= dt * k * div;
        }
    }

    /// Initialises `v^{-1/2}` from `p^0` for a start from rest (velocity zero at `t = 0`);
    /// call once after [`Acoustics::set_pressure`]. Not for PML runs.
    pub fn start_from_rest(&mut self) {
        self.update_v(-0.5 * self.dt);
    }

    /// One leapfrog step (velocity, then pressure).
    pub fn step(&mut self) {
        let dt = self.dt;
        self.update_v(dt);
        self.update_p(dt);
        let tm = self.t + 0.5 * dt;
        for (cell, s) in &self.sources {
            self.p[*cell] += dt * s(tm);
        }
        self.t += dt;
    }

    /// Runs `steps` steps.
    pub fn run(
        &mut self,
        steps: usize,
    ) {
        for _ in 0..steps {
            self.step();
        }
    }
}
