//! # Scalar wave equation / linear acoustics with PML in 1D, 2D and 3D
//!
//! The scalar wave equation `p_tt = c² Δp` is solved in its first-order
//! (linear acoustics) form on a staggered grid,
//! `p_t = -c² ∇·v`, `v_t = -∇p` (unit density), with `p` in the cells and the
//! velocity components on the faces. Boundary faces use
//! [`BoundaryCondition`]:
//!
//! | Face condition | Meaning |
//! |---|---|
//! | `Neumann` / `Symmetry` | rigid wall, `v·n = 0` |
//! | `Dirichlet(0)` | pressure-release surface, `p = 0` at the face |
//! | `Periodic` | wrap-around (both faces of the axis) |
//! | `Absorbing` | convolutional PML of [`AcousticsConfig::pml_cells`] cells ending in a rigid wall |
//!
//! The PML is the stretched-coordinate convolutional PML (`κ = 1`, `α = 0`)
//! with polynomial conductivity grading and one auxiliary variable per
//! derivative; because it acts on the first-order system it is exactly the
//! PML of the scalar wave equation for `p`.
//!
//! Limitations: homogeneous medium, no dispersive attenuation, Dirichlet data
//! must be zero (`p = 0`).

use crate::sim::boundary::{BoundaryCondition, Face, FaceBcs};

type SourceFn = Box<dyn Fn(f64) -> f64 + Send + Sync>;

/// Configuration of the acoustics solver.
pub struct AcousticsConfig {
    /// Cells per axis (inactive axes: 1).
    pub n: [usize; 3],
    /// Active axes.
    pub dim: usize,
    /// Cell size per axis.
    pub h: [f64; 3],
    /// Wave speed.
    pub c: f64,
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
    /// Rigid-wall configuration.
    ///
    /// # Errors
    /// For a bad dimension or sizes.
    pub fn new(n: &[usize], h: &[f64], c: f64) -> Result<Self, String> {
        let dim = n.len();
        if !(1..=3).contains(&dim) || h.len() != dim || n.iter().any(|&m| m < 4) || c <= 0.0 {
            return Err("invalid acoustics grid".into());
        }
        let mut cfg = Self {
            n: [1; 3],
            dim,
            h: [1.0; 3],
            c,
            cfl: 0.9,
            bcs: FaceBcs::uniform(&BoundaryCondition::neumann(0.0)),
            pml_cells: 10,
            pml_order: 3.0,
        };
        for a in 0..dim {
            cfg.n[a] = n[a];
            cfg.h[a] = h[a];
        }
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
    // (b, a) CPML coefficients per axis at cell centres and at faces
    pml_c: [Vec<(f64, f64)>; 3],
    pml_f: [Vec<(f64, f64)>; 3],
    psi_p: [Vec<f64>; 3],
    psi_v: [Vec<f64>; 3],
    use_pml: bool,
    release: [[bool; 2]; 3],
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
        let inv: f64 = (0..cfg.dim).map(|a| 1.0 / (cfg.h[a] * cfg.h[a])).sum();
        let dt = cfg.cfl / (cfg.c * inv.sqrt());
        let mut periodic = [false; 3];
        let mut release = [[false; 2]; 3];
        let mut use_pml = false;
        let mut pml_c: [Vec<(f64, f64)>; 3] = std::array::from_fn(|a| vec![(1.0, 0.0); cfg.n[a]]);
        let mut pml_f: [Vec<(f64, f64)>; 3] = std::array::from_fn(|a| vec![(1.0, 0.0); cfg.n[a] + 1]);
        for a in 0..cfg.dim {
            periodic[a] = cfg.bcs.axis_periodic(a)?;
            for (m, max) in [false, true].into_iter().enumerate() {
                match cfg.bcs.get(Face::of_axis(a, max)) {
                    BoundaryCondition::Neumann(_) | BoundaryCondition::Symmetry | BoundaryCondition::Periodic => {},
                    BoundaryCondition::Dirichlet(g) => {
                        if g.eval(0.0, [0.0; 3]) != 0.0 {
                            return Err("pressure-release faces need zero data".into());
                        }
                        release[a][m] = true;
                    },
                    BoundaryCondition::Absorbing => {
                        let d = cfg.pml_cells;
                        let n = cfg.n[a];
                        if 2 * d >= n {
                            return Err("PML too thick for the grid".into());
                        }
                        use_pml = true;
                        let smax = 0.8 * (cfg.pml_order + 1.0) * cfg.c / cfg.h[a];
                        let prof = |pos: f64| -> Option<(f64, f64)> {
                            let depth = if max { pos - (n - d) as f64 } else { d as f64 - pos };
                            (depth > 0.0).then(|| {
                                let sigma = smax * (depth / d as f64).powf(cfg.pml_order);
                                let b = (-sigma * dt).exp();
                                (b, b - 1.0)
                            })
                        };
                        for i in 0..n {
                            if let Some(c) = prof(i as f64 + 0.5) {
                                pml_c[a][i] = c;
                            }
                        }
                        for i in 0..=n {
                            if let Some(c) = prof(i as f64) {
                                pml_f[a][i] = c;
                            }
                        }
                    },
                    BoundaryCondition::Robin { .. } => return Err("Robin faces are not supported".into()),
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
        let psi_p: [Vec<f64>; 3] = std::array::from_fn(|a| if use_pml && a < cfg.dim { vec![0.0; nn] } else { Vec::new() });
        let psi_v: [Vec<f64>; 3] = std::array::from_fn(|a| if use_pml && a < cfg.dim { vec![0.0; v[a].len()] } else { Vec::new() });
        Ok(Self { cfg, dt, t: 0.0, p: vec![0.0; nn], v, periodic, pml_c, pml_f, psi_p, psi_v, use_pml, release, sources: Vec::new() })
    }

    /// Cells per axis.
    #[must_use]
    pub const fn n(&self) -> [usize; 3] {
        self.cfg.n
    }

    /// Flat cell index.
    #[must_use]
    pub const fn idx(&self, i: usize, j: usize, k: usize) -> usize {
        i + self.cfg.n[0] * (j + self.cfg.n[1] * k)
    }

    /// Sets the pressure from a function of the cell-centre position (velocity zero).
    pub fn set_pressure(&mut self, f: impl Fn([f64; 3]) -> f64) {
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
    pub fn add_source(&mut self, cell: usize, s: impl Fn(f64) -> f64 + Send + Sync + 'static) {
        self.sources.push((cell, Box::new(s)));
    }

    /// Acoustic energy `½ Σ (p²/c² + |v|²) ΔV`.
    #[must_use]
    pub fn energy(&self) -> f64 {
        let dv: f64 = (0..self.cfg.dim).map(|a| self.cfg.h[a]).product();
        let c2 = self.cfg.c * self.cfg.c;
        let mut s: f64 = self.p.iter().map(|x| x * x / c2).sum();
        for a in 0..self.cfg.dim {
            let mut e = self.cfg.n;
            e[a] += 1;
            for (q, x) in self.v[a].iter().enumerate() {
                let ia = (q / e[..a].iter().product::<usize>().max(1)) % e[a];
                let w = if ia == 0 || ia == self.cfg.n[a] { 0.5 } else { 1.0 };
                s += w * x * x;
            }
        }
        0.5 * s * dv
    }

    const fn face_dims(&self, a: usize) -> [usize; 3] {
        let mut e = self.cfg.n;
        e[a] += 1;
        e
    }

    /// Velocity half of the leapfrog step with step size `dt`.
    fn update_v(&mut self, dt: f64) {
        let cfg = &self.cfg;
        let (n, dim) = (cfg.n, cfg.dim);
        // velocity update at faces
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
                let grad = if i == 0 {
                    if self.periodic[a] {
                        (self.p[pidx(0)] - self.p[pidx(n[a] - 1)]) / h
                    } else if self.release[a][0] {
                        self.p[pidx(0)] / (0.5 * h)
                    } else {
                        continue; // rigid wall: v stays zero
                    }
                } else if i == n[a] {
                    if self.periodic[a] {
                        // alias of face 0, set after the sweep
                        continue;
                    } else if self.release[a][1] {
                        -self.p[pidx(n[a] - 1)] / (0.5 * h)
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
                self.v[a][q] -= dt * g;
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

    fn update_p(&mut self, dt: f64) {
        let cfg = &self.cfg;
        let (n, dim) = (cfg.n, cfg.dim);
        let c2 = cfg.c * cfg.c;
        // pressure update at cells
        for p in 0..self.p.len() {
            let ijk = [p % n[0], (p / n[0]) % n[1], p / (n[0] * n[1])];
            let mut div = 0.0;
            for a in 0..dim {
                let e = self.face_dims(a);
                let h = cfg.h[a];
                let mut hi = ijk;
                hi[a] += 1;
                let d = (self.v[a][hi[0] + e[0] * (hi[1] + e[1] * hi[2])] - self.v[a][ijk[0] + e[0] * (ijk[1] + e[1] * ijk[2])]) / h;
                let mut dd = d;
                if self.use_pml {
                    let (b, aa) = self.pml_c[a][ijk[a]];
                    let psi = &mut self.psi_p[a][p];
                    *psi = b * *psi + aa * d;
                    dd += *psi;
                }
                div += dd;
            }
            self.p[p] -= dt * c2 * div;
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
    pub fn run(&mut self, steps: usize) {
        for _ in 0..steps {
            self.step();
        }
    }
}
