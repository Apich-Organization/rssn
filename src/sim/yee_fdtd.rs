//! # Yee FDTD for Maxwell's equations in 1D, 2D and 3D
//!
//! Second-order staggered-grid finite-difference time-domain solver for
//! `ε E_t = ∇×H - J`, `μ H_t = -∇×E` (units `ε0 = μ0 = c = 1`) on a uniform
//! Cartesian Yee lattice. All six field components are always carried; lower
//! dimensions use `∂/∂x_a = 0` along the inactive axes, so one code path
//! covers 1D (`Ey/Hz` and `Ez/Hy` polarisations), 2D (TE and TM) and 3D.
//!
//! | Feature | 1D | 2D | 3D |
//! |---|---|---|---|
//! | Yee update, isotropic dielectric `ε_r(x)` | yes | yes | yes |
//! | PEC walls, periodic faces | yes | yes | yes |
//! | convolutional PML (graded, `κ = 1`, `α = 0`) per face | yes | yes | yes |
//! | soft current source (any component) | yes | yes | yes |
//! | TF/SF plane wave along `+x` with arbitrary `(y, z)` polarisation (1D auxiliary grid) | yes | yes | yes |
//! | TF/SF plane wave on an **arbitrary box with oblique incidence** ([`ObliquePlaneWave`]) | yes | yes | yes |
//!
//! Boundary faces use [`BoundaryCondition`]: `Dirichlet` is a PEC wall,
//! `Periodic` wraps, `Absorbing` is a PML of `pml_cells` cells terminated by
//! a PEC wall. PMC walls (`Neumann`/`Symmetry`) are rejected.
//!
//! The TF/SF box is a slab spanning the full transverse extent between the
//! `x` indices `ia` and `ib` (so it needs transverse periodic, PEC or PML
//! faces that the plane wave satisfies); the incident field comes from a
//! one-dimensional auxiliary Yee grid with the same dispersion, so the
//! scattered-field region is quiet to round-off in vacuum.
//!
//! The oblique TF/SF source evaluates the incident field analytically, so its
//! leakage is set by the grid dispersion (a few per cent at ten points per
//! wavelength; use the auxiliary-grid source for normal incidence when exact
//! cancellation is needed).
//!
//! **Extensions.**
//!
//! * *Dispersive media* ([`Pole`], [`YeeConfig::dispersive`]): any number of
//!   Lorentz / Drude poles per cell through auxiliary differential equations
//!   `P_tt + γ P_t + ω0² P = a E` (central differences, second order). A
//!   Drude pole has `ω0 = 0`; `ε(ω) = ε∞ + Σ a / (ω0² - ω² - i γ ω)`.
//! * *Complex-frequency-shifted PML*: `κ` stretching and `α` shift
//!   ([`YeeConfig::pml_kappa_max`], [`YeeConfig::pml_alpha_max`]) absorb
//!   evanescent and grazing waves that the plain PML reflects.
//! * *Subcell PEC* (Dey–Mittra, [`YeeConfig::pec_body`]): a PEC body given
//!   by a level set modifies the Faraday update with exact edge-length and
//!   face-area fractions (see [`cutcell`](crate::sim::cutcell)); faces with a
//!   free area below [`YeeConfig::dm_min_area`] are frozen and the time step is
//!   reduced to keep the scheme stable. Resonances of curved cavities
//!   converge with second order instead of the first order of a staircase.
//! * *Rotated TF/SF boxes*: [`Yee::set_oblique_plane_wave_rotated`] applies a
//!   rotation to the box of an [`ObliquePlaneWave`] (staircase boundary with
//!   the usual consistent total/scattered correction).
//!
//! Limitations: no anisotropic or magnetic media, no conductivity, the
//! Dey–Mittra treatment is not combined with the TF/SF source or periodic
//! faces, the auxiliary-grid plane wave and the analytic oblique wave assume
//! a vacuum incident medium, PMC boundaries not available.

use crate::sim::boundary::{BoundaryCondition, Face, FaceBcs};

/// Field component.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Comp {
    /// `E_x`.
    Ex,
    /// `E_y`.
    Ey,
    /// `E_z`.
    Ez,
    /// `H_x`.
    Hx,
    /// `H_y`.
    Hy,
    /// `H_z`.
    Hz,
}

/// Soft current source `J(t)` added to an electric component at a cell.
pub struct CurrentSource {
    /// Cell index `(i, j, k)`.
    pub cell: [usize; 3],
    /// Electric component driven.
    pub comp: Comp,
    /// Waveform `J(t)`.
    pub waveform: Box<dyn Fn(f64) -> f64 + Send + Sync>,
}

/// Plane-wave TF/SF source travelling along `+x`.
pub struct PlaneWave {
    /// First total-field `E` node along `x`.
    pub ia: usize,
    /// Last total-field `E` node along `x`.
    pub ib: usize,
    /// Polarisation angle in the `(y, z)` plane (`0` = `E` along `y`).
    pub psi: f64,
    /// Waveform of the incident electric field at the source plane.
    pub waveform: Box<dyn Fn(f64) -> f64 + Send + Sync>,
}

/// Plane wave with arbitrary direction on an arbitrary axis-aligned TF/SF box.
///
/// The incident field is evaluated analytically at the field positions
/// (`E = ê f(t - n̂·(x - x0))`, `H = (n̂ × ê) f(...)`, unit impedance), so
/// the box can have any position and size in 1D, 2D or 3D and the wave any
/// angle of incidence. Because the incident field is not propagated on the
/// grid, its numerical dispersion differs slightly from that of the host
/// lattice and a leakage of order `(k h)² · distance` remains in the
/// scattered-field region.
pub struct ObliquePlaneWave {
    /// Lower corner of the total-field box (physical coordinates).
    pub lo: [f64; 3],
    /// Upper corner of the total-field box.
    pub hi: [f64; 3],
    /// Propagation direction (normalised internally).
    pub direction: [f64; 3],
    /// Electric polarisation (projected orthogonal to `direction`).
    pub polarization: [f64; 3],
    /// Point where the waveform argument equals the time.
    pub origin: [f64; 3],
    /// Waveform `f(τ)`.
    pub waveform: Box<dyn Fn(f64) -> f64 + Send + Sync>,
}

struct ObliqueState {
    rot: Option<[[f64; 3]; 3]>,
    lo: [f64; 3],
    hi: [f64; 3],
    n: [f64; 3],
    e: [f64; 3],
    h: [f64; 3],
    origin: [f64; 3],
    waveform: Box<dyn Fn(f64) -> f64 + Send + Sync>,
}

/// One Lorentz / Drude pole of the electric susceptibility,
/// `P_tt + γ P_t + ω0² P = a E` (unit `ε0`).
#[derive(Clone, Copy, Debug)]
pub struct Pole {
    /// Resonance angular frequency `ω0` (zero for a Drude pole).
    pub omega0: f64,
    /// Damping rate `γ`.
    pub gamma: f64,
    /// Drive strength `a` (`Δε ω0²` for a Lorentz pole, `ω_p²` for a Drude pole).
    pub strength: f64,
}

impl Pole {
    /// Lorentz pole with static contribution `Δε`: `χ = Δε ω0² / (ω0² - ω² - i γ ω)`.
    #[must_use]
    pub fn lorentz(delta_eps: f64, omega0: f64, gamma: f64) -> Self {
        Self { omega0, gamma, strength: delta_eps * omega0 * omega0 }
    }

    /// Drude pole with plasma frequency `ω_p`: `χ = -ω_p² / (ω² + i γ ω)`.
    #[must_use]
    pub const fn drude(omega_p: f64, gamma: f64) -> Self {
        Self { omega0: 0.0, gamma, strength: omega_p * omega_p }
    }

    /// Susceptibility `χ(ω)` (time dependence `e^{-iωt}`) as `(re, im)`.
    #[must_use]
    pub fn susceptibility(&self, omega: f64) -> (f64, f64) {
        let (re, im) = (self.omega0 * self.omega0 - omega * omega, -self.gamma * omega);
        let den = re * re + im * im;
        (self.strength * re / den, -self.strength * im / den)
    }
}

/// Per-cell list of poles as a function of the position.
pub type PoleFn = dyn Fn([f64; 3]) -> Vec<Pole> + Send + Sync;

/// Solver configuration.
pub struct YeeConfig {
    /// Cells per axis (inactive axes: 1).
    pub n: [usize; 3],
    /// Number of active axes.
    pub dim: usize,
    /// Cell size per axis.
    pub h: [f64; 3],
    /// Courant number (`<= 1`); `dt = cfl / sqrt(Σ 1/h²)`.
    pub cfl: f64,
    /// Face conditions.
    pub bcs: FaceBcs,
    /// PML thickness in cells.
    pub pml_cells: usize,
    /// PML grading exponent.
    pub pml_order: f64,
    /// Relative permittivity at `E` nodes (vacuum when `None`); with poles this is `ε∞`.
    pub eps_r: Option<Box<dyn Fn([f64; 3]) -> f64 + Send + Sync>>,
    /// Maximum PML coordinate stretching `κ` (1 = plain PML).
    pub pml_kappa_max: f64,
    /// Maximum complex-frequency shift `α` of the PML (0 = plain PML), in
    /// inverse time units; it is graded to zero towards the PML interface.
    pub pml_alpha_max: f64,
    /// Dispersive poles per position (none when `None`).
    pub dispersive: Option<Box<PoleFn>>,
    /// Level set of a PEC body (`φ <= 0` inside) treated with the Dey–Mittra subcell scheme.
    pub pec_body: Option<Box<dyn Fn([f64; 3]) -> f64 + Send + Sync>>,
    /// Faces of the PEC scheme with a free-area fraction below this value are frozen.
    pub dm_min_area: f64,
}

impl YeeConfig {
    /// Vacuum configuration with PEC walls everywhere.
    ///
    /// # Errors
    /// For a bad dimension or sizes.
    pub fn new(n: &[usize], h: &[f64]) -> Result<Self, String> {
        let dim = n.len();
        if !(1..=3).contains(&dim) || h.len() != dim || n.iter().any(|&m| m < 4) {
            return Err("invalid Yee grid".into());
        }
        let mut c = Self {
            n: [1; 3],
            dim,
            h: [1.0; 3],
            cfl: 0.9,
            bcs: FaceBcs::dirichlet_zero(),
            pml_cells: 10,
            pml_order: 3.0,
            eps_r: None,
            pml_kappa_max: 1.0,
            pml_alpha_max: 0.0,
            dispersive: None,
            pec_body: None,
            dm_min_area: 0.05,
        };
        c.n[..dim].copy_from_slice(&n[..dim]);
        c.h[..dim].copy_from_slice(&h[..dim]);
        Ok(c)
    }
}

/// Yee solver state.
pub struct Yee {
    cfg: YeeConfig,
    /// Time step.
    pub dt: f64,
    /// Current time of the electric field.
    pub t: f64,
    /// Fields `[Ex, Ey, Ez, Hx, Hy, Hz]`.
    pub f: [Vec<f64>; 6],
    eps_inv: Vec<f64>,
    periodic: [bool; 3],
    pec_lo: [bool; 3],
    // PML coefficients per axis: (b, a, kappa) at integer and half positions
    pml_int: [Vec<(f64, f64, f64)>; 3],
    pml_half: [Vec<(f64, f64, f64)>; 3],
    layers: Vec<PoleLayer>,
    dm: Option<DeyMittra>,
    psi: Vec<Vec<f64>>,
    use_pml: bool,
    sources: Vec<CurrentSource>,
    pw: Option<PlaneWaveState>,
    obl: Option<ObliqueState>,
}

/// One pole slot of every cell with the polarisation of the last two steps.
struct PoleLayer {
    strength: Vec<f64>,
    omega2: Vec<f64>,
    gamma: Vec<f64>,
    p: [Vec<f64>; 3],
    pm: [Vec<f64>; 3],
}

/// Dey–Mittra geometry: free edge fractions of `E` and inverse free face areas of `H`.
struct DeyMittra {
    len: [Vec<f64>; 3],
    area_inv: [Vec<f64>; 3],
    scratch: [Vec<f64>; 3],
}

struct PlaneWaveState {
    spec: PlaneWave,
    e: Vec<f64>,
    h: Vec<f64>,
}

const fn comp_index(c: Comp) -> usize {
    match c {
        Comp::Ex => 0,
        Comp::Ey => 1,
        Comp::Ez => 2,
        Comp::Hx => 3,
        Comp::Hy => 4,
        Comp::Hz => 5,
    }
}

/// One term of a curl difference: field component `f`, differentiated along `axis`, with PML memory slot `psi`.
#[derive(Clone, Copy)]
struct CurlTerm {
    f: usize,
    axis: usize,
    psi: usize,
}

impl Yee {
    /// Creates a solver.
    ///
    /// # Errors
    /// For unsupported boundary conditions, inconsistent periodic faces, a
    /// CFL number above one, or a PML thicker than half the grid.
    pub fn new(cfg: YeeConfig) -> Result<Self, String> {
        if !(cfg.cfl > 0.0 && cfg.cfl <= 1.0) {
            return Err("cfl must be in (0, 1]".into());
        }
        let nn: usize = cfg.n.iter().product();
        let inv_sum: f64 = (0..cfg.dim).map(|a| 1.0 / (cfg.h[a] * cfg.h[a])).sum();
        let dt = cfg.cfl / inv_sum.sqrt();
        let mut periodic = [false; 3];
        let mut pec_lo = [false; 3];
        let mut use_pml = false;
        let mut pml_int: [Vec<(f64, f64, f64)>; 3] = std::array::from_fn(|a| vec![(1.0, 0.0, 1.0); cfg.n[a]]);
        let mut pml_half: [Vec<(f64, f64, f64)>; 3] = std::array::from_fn(|a| vec![(1.0, 0.0, 1.0); cfg.n[a]]);
        for a in 0..cfg.dim {
            periodic[a] = cfg.bcs.axis_periodic(a)?;
            pec_lo[a] = !periodic[a];
            for max in [false, true] {
                let bc = cfg.bcs.get(Face::of_axis(a, max));
                match bc {
                    BoundaryCondition::Dirichlet(_) | BoundaryCondition::Periodic => {},
                    BoundaryCondition::Absorbing => {
                        let d = cfg.pml_cells;
                        if 2 * d >= cfg.n[a] {
                            return Err("PML too thick for the grid".into());
                        }
                        use_pml = true;
                        let n = cfg.n[a];
                        let m = cfg.pml_order;
                        let smax = 0.8 * (m + 1.0) / cfg.h[a];
                        for i in 0..n {
                            for half in [false, true] {
                                let pos = i as f64 + if half { 0.5 } else { 0.0 };
                                let depth = if max { pos - (n - d) as f64 } else { d as f64 - pos };
                                if depth > 0.0 {
                                    let frac = (depth / d as f64).min(1.0);
                                    let sigma = smax * frac.powf(m);
                                    let kappa = 1.0 + (cfg.pml_kappa_max - 1.0) * frac.powf(m);
                                    let alpha = cfg.pml_alpha_max * (1.0 - frac);
                                    let b = (-(sigma / kappa + alpha) * dt).exp();
                                    let den = sigma * kappa + kappa * kappa * alpha;
                                    let coef = (b, if den > 0.0 { sigma / den * (b - 1.0) } else { 0.0 }, kappa);
                                    if half {
                                        pml_half[a][i] = coef;
                                    } else {
                                        pml_int[a][i] = coef;
                                    }
                                }
                            }
                        }
                    },
                    _ => return Err("YeeConfig supports PEC (Dirichlet), periodic and PML (Absorbing) faces".into()),
                }
            }
        }
        let eps_inv: Vec<f64> = (0..nn)
            .map(|p| {
                let i = [p % cfg.n[0], (p / cfg.n[0]) % cfg.n[1], p / (cfg.n[0] * cfg.n[1])];
                let x = [i[0] as f64 * cfg.h[0], i[1] as f64 * cfg.h[1], i[2] as f64 * cfg.h[2]];
                1.0 / cfg.eps_r.as_ref().map_or(1.0, |f| f(x))
            })
            .collect();
        let psi = if use_pml { vec![vec![0.0; nn]; 12] } else { Vec::new() };
        let mut layers: Vec<PoleLayer> = Vec::new();
        if let Some(pf) = &cfg.dispersive {
            let per_cell: Vec<Vec<Pole>> = (0..nn)
                .map(|p| {
                    let i = [p % cfg.n[0], (p / cfg.n[0]) % cfg.n[1], p / (cfg.n[0] * cfg.n[1])];
                    pf([i[0] as f64 * cfg.h[0], i[1] as f64 * cfg.h[1], i[2] as f64 * cfg.h[2]])
                })
                .collect();
            let kmax = per_cell.iter().map(Vec::len).max().unwrap_or(0);
            for k in 0..kmax {
                let mut l = PoleLayer {
                    strength: vec![0.0; nn],
                    omega2: vec![0.0; nn],
                    gamma: vec![0.0; nn],
                    p: std::array::from_fn(|_| vec![0.0; nn]),
                    pm: std::array::from_fn(|_| vec![0.0; nn]),
                };
                for (p, poles) in per_cell.iter().enumerate() {
                    if let Some(pole) = poles.get(k) {
                        l.strength[p] = pole.strength;
                        l.omega2[p] = pole.omega0 * pole.omega0;
                        l.gamma[p] = pole.gamma;
                    }
                }
                layers.push(l);
            }
        }
        let mut dt = dt;
        let mut dm = None;
        if let Some(body) = &cfg.pec_body {
            if use_pml || periodic.iter().any(|&b| b) {
                return Err("the subcell PEC scheme is available with PEC walls only".into());
            }
            let outside = |x: [f64; 3]| body(x);
            let mut len: [Vec<f64>; 3] = std::array::from_fn(|_| vec![1.0; nn]);
            let mut area_inv: [Vec<f64>; 3] = std::array::from_fn(|_| vec![1.0; nn]);
            let mut min_area: f64 = 1.0;
            for p in 0..nn {
                let i = [p % cfg.n[0], (p / cfg.n[0]) % cfg.n[1], p / (cfg.n[0] * cfg.n[1])];
                let x0 = [i[0] as f64 * cfg.h[0], i[1] as f64 * cfg.h[1], i[2] as f64 * cfg.h[2]];
                for c in 0..3 {
                    len[c][p] = if c < cfg.dim {
                        crate::sim::cutcell::open_length_fraction(&outside, x0, c, cfg.h[c])
                    } else {
                        f64::from(u8::from(outside(x0) > 0.0))
                    };
                    // face normal to c (H_c): spans the other active axes
                    let axes: Vec<usize> = (0..cfg.dim).filter(|&b| b != c).collect();
                    let mut hl = [0.0; 3];
                    for &b in &axes {
                        hl[b] = cfg.h[b];
                    }
                    let a = if axes.is_empty() {
                        1.0
                    } else {
                        crate::sim::cutcell::open_fraction(&outside, x0, hl, &axes, crate::sim::cutcell::default_depth(axes.len()))
                    };
                    if a >= cfg.dm_min_area {
                        area_inv[c][p] = 1.0 / a;
                        if a < 1.0 - 1e-12 {
                            min_area = min_area.min(a);
                        }
                    } else {
                        area_inv[c][p] = 0.0;
                    }
                }
            }
            dt *= min_area.sqrt();
            dm = Some(DeyMittra { len, area_inv, scratch: std::array::from_fn(|_| vec![0.0; nn]) });
        }
        Ok(Self {
            cfg,
            dt,
            layers,
            dm,
            t: 0.0,
            f: std::array::from_fn(|_| vec![0.0; nn]),
            eps_inv,
            periodic,
            pec_lo,
            pml_int,
            pml_half,
            psi,
            use_pml,
            sources: Vec::new(),
            pw: None,
            obl: None,
        })
    }

    /// Adds a soft current source.
    pub fn add_source(&mut self, s: CurrentSource) {
        self.sources.push(s);
    }

    /// Installs the TF/SF plane-wave source.
    ///
    /// # Errors
    /// If the slab is outside the grid or touches the PML.
    pub fn set_plane_wave(&mut self, pw: PlaneWave, total_steps: usize) -> Result<(), String> {
        let nx = self.cfg.n[0];
        if pw.ia < 2 || pw.ib + 2 >= nx || pw.ia >= pw.ib {
            return Err("TF/SF slab must lie inside the grid".into());
        }
        let n_aux = nx + (total_steps as f64 * self.cfg.cfl).ceil() as usize + 20;
        self.pw = Some(PlaneWaveState { spec: pw, e: vec![0.0; n_aux], h: vec![0.0; n_aux] });
        Ok(())
    }

    /// Installs an oblique plane wave on an arbitrary TF/SF box.
    ///
    /// # Errors
    /// For a zero direction or polarisation, a polarisation parallel to the
    /// direction, or components along inactive axes.
    pub fn set_oblique_plane_wave(&mut self, w: ObliquePlaneWave) -> Result<(), String> {
        let norm = |v: [f64; 3]| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        let nn = norm(w.direction);
        if nn == 0.0 {
            return Err("direction must be non-zero".into());
        }
        let n = [w.direction[0] / nn, w.direction[1] / nn, w.direction[2] / nn];
        let dotp = n[0] * w.polarization[0] + n[1] * w.polarization[1] + n[2] * w.polarization[2];
        let mut e = [w.polarization[0] - dotp * n[0], w.polarization[1] - dotp * n[1], w.polarization[2] - dotp * n[2]];
        let en = norm(e);
        if en < 1e-12 {
            return Err("polarisation must not be parallel to the direction".into());
        }
        for c in &mut e {
            *c /= en;
        }
        for a in self.cfg.dim..3 {
            if n[a].abs() > 1e-12 {
                return Err("direction has a component along an inactive axis".into());
            }
        }
        let h = [n[1] * e[2] - n[2] * e[1], n[2] * e[0] - n[0] * e[2], n[0] * e[1] - n[1] * e[0]];
        self.obl = Some(ObliqueState { rot: None, lo: w.lo, hi: w.hi, n, e, h, origin: w.origin, waveform: w.waveform });
        Ok(())
    }

    /// Installs an oblique plane wave on a **rotated** box: the box
    /// `[lo, hi]` of `w` is rotated about its centre by the orthonormal
    /// matrix `rotation` (columns are the images of the axes). The incident
    /// wave itself keeps the direction and polarisation given in `w`.
    ///
    /// # Errors
    /// As [`Yee::set_oblique_plane_wave`], and for a matrix that is not
    /// orthonormal or that mixes in an inactive axis.
    pub fn set_oblique_plane_wave_rotated(&mut self, w: ObliquePlaneWave, rotation: [[f64; 3]; 3]) -> Result<(), String> {
        for i in 0..3 {
            for j in 0..3 {
                let dot: f64 = (0..3).map(|k| rotation[k][i] * rotation[k][j]).sum();
                if (dot - f64::from(u8::from(i == j))).abs() > 1e-9 {
                    return Err("the rotation must be orthonormal".into());
                }
            }
        }
        for a in self.cfg.dim..3 {
            if (0..3).any(|b| b != a && (rotation[a][b].abs() > 1e-12 || rotation[b][a].abs() > 1e-12)) {
                return Err("the rotation must not mix in inactive axes".into());
            }
        }
        self.set_oblique_plane_wave(w)?;
        if let Some(o) = &mut self.obl {
            o.rot = Some(rotation);
        }
        Ok(())
    }

    /// Physical position of a field component at (possibly out-of-range) cell indices.
    fn pos_of(&self, comp: usize, ijk: [isize; 3]) -> [f64; 3] {
        let mut x = [0.0; 3];
        for a in 0..self.cfg.dim {
            let off = if comp < 3 { f64::from(u8::from(a == comp)) } else { f64::from(u8::from(a != comp - 3)) };
            x[a] = (ijk[a] as f64 + 0.5 * off) * self.cfg.h[a];
        }
        x
    }

    fn in_box(&self, pos: [f64; 3]) -> f64 {
        let Some(o) = &self.obl else { return 0.0 };
        let mut pos = pos;
        if let Some(r) = &o.rot {
            // local coordinates in the frame of the rotated box (rotation about its centre)
            let c: [f64; 3] = std::array::from_fn(|a| 0.5 * (o.lo[a] + o.hi[a]));
            let d: [f64; 3] = std::array::from_fn(|a| pos[a] - c[a]);
            pos = std::array::from_fn(|a| c[a] + (0..3).map(|b| r[b][a] * d[b]).sum::<f64>());
        }
        for a in 0..self.cfg.dim {
            let tol = 1e-9 * self.cfg.h[a];
            if pos[a] < o.lo[a] - tol || pos[a] > o.hi[a] + tol {
                return 0.0;
            }
        }
        1.0
    }

    fn incident(&self, comp: usize, pos: [f64; 3], t: f64) -> f64 {
        let Some(o) = &self.obl else { return 0.0 };
        let s: f64 = (0..3).map(|a| o.n[a] * (pos[a] - o.origin[a])).sum();
        let f = (o.waveform)(t - s);
        if comp < 3 { o.e[comp] * f } else { o.h[comp - 3] * f }
    }

    /// Correction of the difference quotient `∂_a G` used to update component
    /// `fc` so that total and scattered fields are mixed consistently.
    fn tf_corr(&self, fc: usize, g: usize, ijk: [usize; 3], a: usize, forward: bool) -> f64 {
        if self.obl.is_none() || self.cfg.n[a] == 1 {
            return 0.0;
        }
        let n = self.cfg.n[a];
        let cur = [ijk[0] as isize, ijk[1] as isize, ijk[2] as isize];
        let mut nb = cur;
        if forward {
            if ijk[a] + 1 >= n && !self.periodic[a] {
                return 0.0;
            }
            nb[a] += 1;
        } else {
            if ijk[a] == 0 && !self.periodic[a] {
                return 0.0;
            }
            nb[a] -= 1;
        }
        let mf = self.in_box(self.pos_of(fc, cur));
        let (pc, pn) = (self.pos_of(g, cur), self.pos_of(g, nb));
        let (mc, mn) = (self.in_box(pc), self.in_box(pn));
        if crate::kernels::exact_eq(mf, mc) && crate::kernels::exact_eq(mf, mn) {
            return 0.0;
        }
        let tg = if g < 3 { self.t } else { self.t + 0.5 * self.dt };
        let dcur = if crate::kernels::exact_eq(mf, mc) { 0.0 } else { (mf - mc) * self.incident(g, pc, tg) };
        let dnb = if crate::kernels::exact_eq(mf, mn) { 0.0 } else { (mf - mn) * self.incident(g, pn, tg) };
        let h = self.cfg.h[a];
        if forward { (dnb - dcur) / h } else { (dcur - dnb) / h }
    }

    /// Number of cells.
    #[must_use]
    pub fn size(&self) -> usize {
        self.cfg.n.iter().product()
    }

    /// Flat index of cell `(i, j, k)`.
    #[must_use]
    pub const fn idx(&self, i: usize, j: usize, k: usize) -> usize {
        i + self.cfg.n[0] * (j + self.cfg.n[1] * k)
    }

    /// Field values of one component.
    #[must_use]
    pub fn field(&self, c: Comp) -> &[f64] {
        &self.f[comp_index(c)]
    }

    /// Mutable field values of one component (for initial conditions).
    pub fn field_mut(&mut self, c: Comp) -> &mut [f64] {
        &mut self.f[comp_index(c)]
    }

    /// Cell size along an axis.
    #[must_use]
    pub const fn h(&self, a: usize) -> f64 {
        self.cfg.h[a]
    }

    /// Cells per axis.
    #[must_use]
    pub const fn n(&self) -> [usize; 3] {
        self.cfg.n
    }

    /// Electromagnetic energy `½ Σ (ε E² + H²) ΔV`.
    #[must_use]
    pub fn energy(&self) -> f64 {
        let dv: f64 = (0..self.cfg.dim).map(|a| self.cfg.h[a]).product();
        let mut s = 0.0;
        for p in 0..self.size() {
            for c in 0..3 {
                s += self.f[c][p] * self.f[c][p] / self.eps_inv[p] + self.f[3 + c][p] * self.f[3 + c][p];
            }
        }
        0.5 * s * dv
    }

    /// Forward difference `(f[i+1] - f[i]) / h` along `a` (zero beyond a PEC wall).
    fn dfw(&self, f: &[f64], ijk: [usize; 3], a: usize) -> f64 {
        let n = self.cfg.n[a];
        if n == 1 {
            return 0.0;
        }
        let mut q = ijk;
        let nb = if ijk[a] + 1 < n {
            q[a] += 1;
            f[self.idx(q[0], q[1], q[2])]
        } else if self.periodic[a] {
            q[a] = 0;
            f[self.idx(q[0], q[1], q[2])]
        } else {
            0.0
        };
        (nb - f[self.idx(ijk[0], ijk[1], ijk[2])]) / self.cfg.h[a]
    }

    /// Backward difference `(f[i] - f[i-1]) / h` along `a`.
    fn dbw(&self, f: &[f64], ijk: [usize; 3], a: usize) -> f64 {
        let n = self.cfg.n[a];
        if n == 1 {
            return 0.0;
        }
        let mut q = ijk;
        let nb = if ijk[a] > 0 {
            q[a] -= 1;
            f[self.idx(q[0], q[1], q[2])]
        } else if self.periodic[a] {
            q[a] = n - 1;
            f[self.idx(q[0], q[1], q[2])]
        } else {
            0.0
        };
        (f[self.idx(ijk[0], ijk[1], ijk[2])] - nb) / self.cfg.h[a]
    }

    /// Curl-difference `∂_a1 F1 - ∂_a2 F2` for the H update (forward differences) with PML memory.
    fn curl_h(&mut self, fc: usize, p: usize, ijk: [usize; 3], t1: CurlTerm, t2: CurlTerm) -> f64 {
        let CurlTerm { f: f1, axis: a1, psi: m1 } = t1;
        let CurlTerm { f: f2, axis: a2, psi: m2 } = t2;
        let mut d1 = self.dfw(&self.f[f1], ijk, a1) + self.tf_corr(fc, f1, ijk, a1, true);
        let mut d2 = self.dfw(&self.f[f2], ijk, a2) + self.tf_corr(fc, f2, ijk, a2, true);
        if self.use_pml {
            let (b, a, k) = self.pml_half[a1][ijk[a1]];
            self.psi[m1][p] = b * self.psi[m1][p] + a * d1;
            d1 = d1 / k + self.psi[m1][p];
            let (b, a, k) = self.pml_half[a2][ijk[a2]];
            self.psi[m2][p] = b * self.psi[m2][p] + a * d2;
            d2 = d2 / k + self.psi[m2][p];
        }
        d1 - d2
    }

    fn curl_e(&mut self, fc: usize, p: usize, ijk: [usize; 3], t1: CurlTerm, t2: CurlTerm) -> f64 {
        let CurlTerm { f: f1, axis: a1, psi: m1 } = t1;
        let CurlTerm { f: f2, axis: a2, psi: m2 } = t2;
        let mut d1 = self.dbw(&self.f[f1], ijk, a1) + self.tf_corr(fc, f1, ijk, a1, false);
        let mut d2 = self.dbw(&self.f[f2], ijk, a2) + self.tf_corr(fc, f2, ijk, a2, false);
        if self.use_pml {
            let (b, a, k) = self.pml_int[a1][ijk[a1]];
            self.psi[m1][p] = b * self.psi[m1][p] + a * d1;
            d1 = d1 / k + self.psi[m1][p];
            let (b, a, k) = self.pml_int[a2][ijk[a2]];
            self.psi[m2][p] = b * self.psi[m2][p] + a * d2;
            d2 = d2 / k + self.psi[m2][p];
        }
        d1 - d2
    }

    /// Advances `H` by `dt_h` using the current `E`.
    fn update_h(&mut self, dt_h: f64) {
        let n = self.cfg.n;
        if let Some(dm) = &mut self.dm {
            // the Faraday circulation uses the free parts of the edges: E -> L E
            for c in 0..3 {
                for p in 0..self.f[c].len() {
                    dm.scratch[c][p] = self.f[c][p] * dm.len[c][p];
                }
                std::mem::swap(&mut self.f[c], &mut dm.scratch[c]);
            }
        }
        for k in 0..n[2] {
            for j in 0..n[1] {
                for i in 0..n[0] {
                    let ijk = [i, j, k];
                    let p = self.idx(i, j, k);
                    // Hx: dy Ez - dz Ey ; Hy: dz Ex - dx Ez ; Hz: dx Ey - dy Ex
                    let cx = self.curl_h(3, p, ijk, CurlTerm { f: 2, axis: 1, psi: 0 }, CurlTerm { f: 1, axis: 2, psi: 1 });
                    let cy = self.curl_h(4, p, ijk, CurlTerm { f: 0, axis: 2, psi: 2 }, CurlTerm { f: 2, axis: 0, psi: 3 });
                    let cz = self.curl_h(5, p, ijk, CurlTerm { f: 1, axis: 0, psi: 4 }, CurlTerm { f: 0, axis: 1, psi: 5 });
                    let (sx, sy, sz) = self.dm.as_ref().map_or((1.0, 1.0, 1.0), |d| (d.area_inv[0][p], d.area_inv[1][p], d.area_inv[2][p]));
                    self.f[3][p] -= dt_h * cx * sx;
                    self.f[4][p] -= dt_h * cy * sy;
                    self.f[5][p] -= dt_h * cz * sz;
                }
            }
        }
        if let Some(dm) = &mut self.dm {
            for c in 0..3 {
                std::mem::swap(&mut self.f[c], &mut dm.scratch[c]);
            }
        }
    }

    fn update_e(&mut self) {
        let n = self.cfg.n;
        let dt = self.dt;
        for k in 0..n[2] {
            for j in 0..n[1] {
                for i in 0..n[0] {
                    let ijk = [i, j, k];
                    let p = self.idx(i, j, k);
                    // Ex: dy Hz - dz Hy ; Ey: dz Hx - dx Hz ; Ez: dx Hy - dy Hx
                    let cx = self.curl_e(0, p, ijk, CurlTerm { f: 5, axis: 1, psi: 6 }, CurlTerm { f: 4, axis: 2, psi: 7 });
                    let cy = self.curl_e(1, p, ijk, CurlTerm { f: 3, axis: 2, psi: 8 }, CurlTerm { f: 5, axis: 0, psi: 9 });
                    let cz = self.curl_e(2, p, ijk, CurlTerm { f: 4, axis: 0, psi: 10 }, CurlTerm { f: 3, axis: 1, psi: 11 });
                    let ei = self.eps_inv[p];
                    let mut dp = [0.0; 3];
                    for layer in &mut self.layers {
                        if layer.strength[p] == 0.0 && layer.omega2[p] == 0.0 && layer.gamma[p] == 0.0 {
                            continue;
                        }
                        let (g, o2, a) = (layer.gamma[p], layer.omega2[p], layer.strength[p]);
                        for c in 0..3 {
                            let (pn, pnm) = (layer.p[c][p], layer.pm[c][p]);
                            let pnew = ((2.0 - o2 * dt * dt) * pn - (1.0 - 0.5 * g * dt) * pnm + a * dt * dt * self.f[c][p]) / (1.0 + 0.5 * g * dt);
                            dp[c] += pnew - pn;
                            layer.pm[c][p] = pn;
                            layer.p[c][p] = pnew;
                        }
                    }
                    self.f[0][p] += ei * (dt * cx - dp[0]);
                    self.f[1][p] += ei * (dt * cy - dp[1]);
                    self.f[2][p] += ei * (dt * cz - dp[2]);
                    if let Some(dm) = &self.dm {
                        for c in 0..3 {
                            if dm.len[c][p] <= 0.0 {
                                self.f[c][p] = 0.0;
                            }
                        }
                    }
                }
            }
        }
        // PEC at the lower faces: tangential E vanishes
        for k in 0..n[2] {
            for j in 0..n[1] {
                for i in 0..n[0] {
                    let ijk = [i, j, k];
                    let p = self.idx(i, j, k);
                    for a in 0..self.cfg.dim {
                        if self.pec_lo[a] && ijk[a] == 0 {
                            for c in 0..3 {
                                if c != a {
                                    self.f[c][p] = 0.0;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /// Applies the TF/SF corrections of the H update (call after `update_h`).
    fn tfsf_h(&mut self, dt_h: f64) {
        let Some(pw) = &self.pw else { return };
        let (ia, ib) = (pw.spec.ia, pw.spec.ib);
        let (s, c) = pw.spec.psi.sin_cos();
        let (einc_a, einc_b) = (pw.e[ia], pw.e[ib]);
        let hx = self.cfg.h[0];
        let [_, ny, nz] = self.cfg.n;
        for k in 0..nz {
            for j in 0..ny {
                let pl = self.idx(ia - 1, j, k);
                let pr = self.idx(ib, j, k);
                // pair (Ey, Hz): amplitude c ; pair (Ez, Hy' = -Hy): amplitude s
                self.f[5][pl] += dt_h / hx * c * einc_a;
                self.f[5][pr] -= dt_h / hx * c * einc_b;
                self.f[4][pl] -= dt_h / hx * s * einc_a;
                self.f[4][pr] += dt_h / hx * s * einc_b;
            }
        }
    }

    fn tfsf_e(&mut self, hinc_at: &dyn Fn(usize) -> f64) {
        let Some(pw) = &self.pw else { return };
        let (ia, ib) = (pw.spec.ia, pw.spec.ib);
        let (s, c) = pw.spec.psi.sin_cos();
        let (ha, hb) = (hinc_at(ia - 1), hinc_at(ib));
        let hx = self.cfg.h[0];
        let dt = self.dt;
        let [_, ny, nz] = self.cfg.n;
        for k in 0..nz {
            for j in 0..ny {
                let pl = self.idx(ia, j, k);
                let pr = self.idx(ib, j, k);
                let (el, er) = (self.eps_inv[pl], self.eps_inv[pr]);
                self.f[1][pl] += dt / hx * el * c * ha;
                self.f[1][pr] -= dt / hx * er * c * hb;
                self.f[2][pl] += dt / hx * el * s * ha;
                self.f[2][pr] -= dt / hx * er * s * hb;
            }
        }
    }

    /// Sets `H^{-1/2}` from `E^0` assuming `H^0 = 0` (call once after setting
    /// `E` for a standing-wave initial condition; not for PML runs).
    pub fn half_step_h(&mut self) {
        self.update_h(-0.5 * self.dt);
    }

    /// One leapfrog step (`H` then `E`).
    pub fn step(&mut self) {
        let dt = self.dt;
        self.update_h(dt);
        self.tfsf_h(dt);
        // aux H update
        if let Some(pw) = &mut self.pw {
            let hx = self.cfg.h[0];
            let n = pw.e.len();
            for i in 0..n - 1 {
                pw.h[i] -= dt / hx * (pw.e[i + 1] - pw.e[i]);
            }
        }
        self.update_e();
        let h_snapshot: Option<Vec<f64>> = self.pw.as_ref().map(|p| p.h.clone());
        if let Some(hs) = &h_snapshot {
            self.tfsf_e(&|i| hs[i]);
        }
        if let Some(pw) = &mut self.pw {
            let hx = self.cfg.h[0];
            let n = pw.e.len();
            for i in 1..n {
                let hm = if i > 0 { pw.h[i - 1] } else { 0.0 };
                pw.e[i] -= dt / hx * (pw.h[i] - hm);
            }
            pw.e[0] = (pw.spec.waveform)(self.t + dt);
        }
        let t_now = self.t + 0.5 * dt;
        for s in &self.sources {
            let p = self.idx(s.cell[0], s.cell[1], s.cell[2]);
            let ci = comp_index(s.comp);
            self.f[ci][p] -= dt * self.eps_inv[p] * (s.waveform)(t_now);
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
