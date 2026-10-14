//! # Dimension-generic molecular dynamics engine
//!
//! A complete MD toolkit in reduced units (`k_B = 1`) for any dimension `D`
//! (const generic; 1, 2 and 3 are tested).
//!
//! | Feature | Support |
//! |---|---|
//! | Boundaries per axis ([`Boundary`]) | open, periodic (minimum image), reflecting (specular) walls, thermalising walls (flux-Maxwellian re-emission); Lees-Edwards shear (flow axis 0, gradient axis 1, `D >= 2`) |
//! | Initial conditions | simple-cubic / fcc (3-D) / triangular (2-D) lattices, random sequential (hard-core) placement, Maxwell-Boltzmann velocities with zero total momentum and exact temperature, user data by direct assignment |
//! | Pair potentials ([`PairTerm`]) | Lennard-Jones, Morse, Coulomb (bare or `erfc`-screened for Ewald real space), optional energy shift |
//! | Neighbour search ([`NeighbourList`]) | cell list + Verlet skin in any `D` (brute force under Lees-Edwards) |
//! | Bonded forces | harmonic bonds, harmonic angles (any `D`), periodic dihedrals (`D = 3`) |
//! | Thermostats ([`Thermostat`]) | Berendsen, Langevin (BAOAB), Nose-Hoover chains (Suzuki-Yoshida) |
//! | Barostats ([`Barostat`]) | Berendsen, isotropic MTK (Martyna-Tobias-Klein), and the flexible-cell Parrinello-Rahman / MTK barostat [`FlexibleBarostat`] (full triclinic cell matrix, orthorhombic or isotropic modes, target pressure tensor, cell Nose-Hoover chain = MTTK) |
//! | Cells | orthorhombic, and general triclinic cells ([`SimBox::triclinic`], cell lists in fractional coordinates, 3^D image refinement of the minimum image) |
//! | Constraints | SHAKE (positions) and RATTLE (velocities) for distance constraints; composable with Nose-Hoover chains, Langevin (constrained BAOAB) and the flexible barostat (bond flow `r.(v + v_g r) = 0`, constraint virial in the pressure) |
//! | Rigid molecules | [`crate::sim::md_rigid`]: centre-of-mass + quaternion integration (NO_SQUISH or exact symmetric-top flows, any even order), molecular virial |
//! | Observables | kinetic/potential energy, temperature, pressure tensor, [`Rdf`], [`msd`], [`vacf`], [`energy_drift_slope`] |
//!
//! Long-range electrostatics (Ewald, PME) live in [`crate::sim::ewald`] and plug
//! in through [`ForceField::extra`]. The older [`crate::sim::physics_md`] module
//! keeps its runtime-dimension `Vec<f64>` API and delegates its built-in force field,
//! neighbour search and RDF to this engine.
//!
//! Limitations: the legacy isotropic [`Barostat::Mtk`] is NPH-type and excludes thermostats
//! and constraints (use [`FlexibleBarostat`] with [`CellMode::Isotropic`] for those
//! combinations); the flexible cell needs a fully periodic box without Lees-Edwards shear and
//! does not compose with Langevin or Berendsen barostat coupling; Ewald/PME
//! ([`crate::sim::ewald`]) assumes an orthorhombic box; the constraint virial of the
//! flexible barostat is the SHAKE impulse (first-order lag in the time step); the
//! neighbour list of a triclinic cell is rebuilt when the affine deformation plus particle
//! displacement exceed the skin.
#![allow(
    clippy::arithmetic_side_effects,
    clippy::assigning_clones,
    clippy::cast_lossless,
    clippy::cast_possible_truncation,
    clippy::chunks_exact_to_as_chunks,
    clippy::cognitive_complexity,
    clippy::derive_partial_eq_without_eq,
    clippy::doc_markdown,
    clippy::excessive_precision,
    clippy::float_cmp,
    clippy::if_not_else,
    clippy::imprecise_flops,
    clippy::indexing_slicing,
    clippy::int_plus_one,
    clippy::items_after_statements,
    clippy::manual_let_else,
    clippy::manual_memcpy,
    clippy::manual_midpoint,
    clippy::many_single_char_names,
    clippy::missing_const_for_fn,
    clippy::module_name_repetitions,
    clippy::needless_pass_by_value,
    clippy::needless_range_loop,
    clippy::option_if_let_else,
    clippy::similar_names,
    clippy::single_call_fn,
    clippy::single_match_else,
    clippy::suboptimal_flops,
    clippy::suspicious_operation_groupings,
    clippy::too_long_first_doc_paragraph,
    clippy::too_many_arguments,
    clippy::too_many_lines,
    clippy::type_complexity,
    clippy::unreadable_literal,
    clippy::useless_let_if_seq,
    clippy::while_float
)]

use std::collections::HashSet;

use crate::kernels::random::Rng;
use crate::kernels::special::erfc_numerical;

#[path = "md_cell.rs"]
mod cell;

pub use cell::{mat_det, mat_exp_phi1, mat_inverse, mat_mul, mat_vec, CellMode, FlexibleBarostat};

/// Boundary condition of one box axis.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Boundary {
    /// No boundary: particles may leave.
    Open,
    /// Periodic with minimum image.
    Periodic,
    /// Specular reflecting walls at `0` and `L`.
    Reflecting,
    /// Walls that re-emit particles with a flux-Maxwellian at the given temperature.
    Thermal(f64),
}

/// Simulation box (orthorhombic).
#[derive(Clone, Debug)]
pub struct SimBox<const D: usize> {
    /// Edge lengths.
    pub len: [f64; D],
    /// Boundary per axis.
    pub bc: [Boundary; D],
    /// Lees-Edwards shear rate `du_x/dy` (requires `D >= 2`, axes 0 and 1 periodic).
    pub shear_rate: f64,
    /// Current Lees-Edwards displacement of the upper image along `x`.
    pub offset: f64,
    /// General (triclinic) cell: row `k` is the lattice vector `a_k`. `None` for an orthorhombic
    /// box. Set it through [`SimBox::set_cell`], which keeps `len` (the vector norms) and the
    /// cached inverse in sync. Triclinic boxes are fully periodic and have no Lees-Edwards shear.
    pub cell: Option<[[f64; D]; D]>,
    cell_inv: [[f64; D]; D],
}

impl<const D: usize> SimBox<D> {
    /// Fully periodic box.
    #[must_use]
    pub fn periodic(len: [f64; D]) -> Self {
        Self { len, bc: [Boundary::Periodic; D], shear_rate: 0.0, offset: 0.0, cell: None, cell_inv: [[0.0; D]; D] }
    }

    /// Fully periodic triclinic box with lattice vectors `h[k]` (rows).
    ///
    /// # Errors
    /// Fails if the lattice vectors are linearly dependent.
    pub fn triclinic(h: [[f64; D]; D]) -> Result<Self, String> {
        let mut b = Self::periodic([1.0; D]);
        b.set_cell(h)?;
        Ok(b)
    }

    /// Sets the general cell matrix (rows are lattice vectors), refreshing `len` and the inverse.
    ///
    /// # Errors
    /// Fails if the matrix is singular.
    pub fn set_cell(
        &mut self,
        h: [[f64; D]; D],
    ) -> Result<(), String> {
        let inv = crate::sim::md_engine::mat_inverse(&h).ok_or_else(|| "cell matrix is singular".to_string())?;
        self.cell = Some(h);
        self.cell_inv = inv;
        for k in 0..D {
            self.len[k] = h[k].iter().map(|x| x * x).sum::<f64>().sqrt();
        }
        self.bc = [Boundary::Periodic; D];
        self.shear_rate = 0.0;
        Ok(())
    }

    /// Cell matrix (rows are lattice vectors; diagonal for an orthorhombic box).
    #[must_use]
    pub fn matrix(&self) -> [[f64; D]; D] {
        self.cell.unwrap_or_else(|| {
            let mut h = [[0.0; D]; D];
            for k in 0..D {
                h[k][k] = self.len[k];
            }
            h
        })
    }

    /// Fractional coordinates `s` with `x = sum_k s_k a_k` (orthorhombic: `x / L`).
    #[must_use]
    pub fn to_fractional(
        &self,
        x: &[f64; D],
    ) -> [f64; D] {
        let mut s = [0.0; D];
        if self.cell.is_some() {
            for a in 0..D {
                for b in 0..D {
                    s[a] += x[b] * self.cell_inv[b][a];
                }
            }
        } else {
            for a in 0..D {
                s[a] = x[a] / self.len[a];
            }
        }
        s
    }

    /// Cartesian vector of fractional coordinates `s`.
    #[must_use]
    pub fn from_fractional(
        &self,
        s: &[f64; D],
    ) -> [f64; D] {
        let h = self.matrix();
        let mut x = [0.0; D];
        for a in 0..D {
            for b in 0..D {
                x[b] += s[a] * h[a][b];
            }
        }
        x
    }

    /// Perpendicular width of the cell along lattice direction `a` (distance between faces).
    #[must_use]
    pub fn height(
        &self,
        a: usize,
    ) -> f64 {
        if self.cell.is_some() {
            1.0 / (0..D).map(|b| self.cell_inv[b][a].powi(2)).sum::<f64>().sqrt()
        } else {
            self.len[a]
        }
    }

    /// Box with explicit per-axis boundaries.
    #[must_use]
    pub fn new(
        len: [f64; D],
        bc: [Boundary; D],
    ) -> Self {
        Self { len, bc, shear_rate: 0.0, offset: 0.0, cell: None, cell_inv: [[0.0; D]; D] }
    }

    /// Volume (length in 1-D, area in 2-D).
    #[must_use]
    pub fn volume(&self) -> f64 {
        if let Some(h) = &self.cell {
            return mat_det(h).abs();
        }
        self.len.iter().product()
    }

    /// Applies the minimum-image convention (including Lees-Edwards offset) to a separation.
    pub fn min_image(
        &self,
        d: &mut [f64; D],
    ) {
        if let Some(h) = &self.cell {
            let mut s = self.to_fractional(d);
            for c in &mut s {
                *c -= c.round();
            }
            let mut base = [0.0; D];
            for a in 0..D {
                for b in 0..D {
                    base[b] += s[a] * h[a][b];
                }
            }
            // a skewed cell can have a shorter image in a neighbouring cell
            let mut best = base;
            let mut best_n = dot(&base, &base);
            for r0 in 0..3usize.pow(D as u32) {
                let mut r = r0;
                let mut cand = base;
                for a in 0..D {
                    let o = (r % 3) as f64 - 1.0;
                    r /= 3;
                    for b in 0..D {
                        cand[b] += o * h[a][b];
                    }
                }
                let n = dot(&cand, &cand);
                if n < best_n {
                    best_n = n;
                    best = cand;
                }
            }
            *d = best;
            return;
        }
        if D >= 2 && self.shear_rate != 0.0 {
            let n = (d[1] / self.len[1]).round();
            d[1] -= n * self.len[1];
            d[0] -= n * self.offset;
        }
        for k in 0..D {
            if self.bc[k] == Boundary::Periodic {
                d[k] -= self.len[k] * (d[k] / self.len[k]).round();
            }
        }
    }
}

/// A particle system.
#[derive(Clone, Debug)]
pub struct MdSystem<const D: usize> {
    /// Positions (wrapped into the box on periodic axes).
    pub pos: Vec<[f64; D]>,
    /// Velocities (peculiar velocities are not separated under Lees-Edwards).
    pub vel: Vec<[f64; D]>,
    /// Continuous positions without periodic wrapping (for MSD).
    pub unwrapped: Vec<[f64; D]>,
    /// Masses.
    pub mass: Vec<f64>,
    /// Charges.
    pub charge: Vec<f64>,
    /// Box.
    pub bx: SimBox<D>,
    /// Elapsed time.
    pub time: f64,
}

impl<const D: usize> MdSystem<D> {
    /// Builds a system at rest from positions with unit masses and zero charges.
    #[must_use]
    pub fn new(
        pos: Vec<[f64; D]>,
        bx: SimBox<D>,
    ) -> Self {
        let n = pos.len();
        Self { unwrapped: pos.clone(), vel: vec![[0.0; D]; n], mass: vec![1.0; n], charge: vec![0.0; n], pos, bx, time: 0.0 }
    }

    /// Number of particles.
    #[must_use]
    pub fn n(&self) -> usize {
        self.pos.len()
    }

    /// Kinetic energy.
    #[must_use]
    pub fn kinetic_energy(&self) -> f64 {
        let mut k = 0.0;
        for i in 0..self.n() {
            for a in 0..D {
                k += 0.5 * self.mass[i] * self.vel[i][a] * self.vel[i][a];
            }
        }
        k
    }

    /// Total momentum.
    #[must_use]
    pub fn momentum(&self) -> [f64; D] {
        let mut p = [0.0; D];
        for i in 0..self.n() {
            for a in 0..D {
                p[a] += self.mass[i] * self.vel[i][a];
            }
        }
        p
    }

    /// Removes the centre-of-mass velocity.
    pub fn remove_com_velocity(&mut self) {
        let m: f64 = self.mass.iter().sum();
        let p = self.momentum();
        for i in 0..self.n() {
            for a in 0..D {
                self.vel[i][a] -= p[a] / m;
            }
        }
    }

    /// Temperature `2 K / n_dof` for the given number of degrees of freedom.
    #[must_use]
    pub fn temperature(
        &self,
        dof: f64,
    ) -> f64 {
        2.0 * self.kinetic_energy() / dof
    }

    /// Draws Maxwell-Boltzmann velocities at temperature `t`, removes the drift and rescales to exactly `t`.
    pub fn maxwell_boltzmann(
        &mut self,
        t: f64,
        rng: &mut Rng,
    ) {
        for i in 0..self.n() {
            let s = (t / self.mass[i]).sqrt();
            for a in 0..D {
                self.vel[i][a] = s * rng.normal();
            }
        }
        if self.n() > 1 {
            self.remove_com_velocity();
        }
        let dof = if self.n() > 1 { (D * (self.n() - 1)) as f64 } else { D as f64 };
        let cur = self.temperature(dof);
        if cur > 0.0 {
            let f = (t / cur).sqrt();
            for v in &mut self.vel {
                for c in v.iter_mut() {
                    *c *= f;
                }
            }
        }
    }

    /// Applies boundary conditions after a displacement `dx` of every particle already added to `pos`.
    /// Periodic axes wrap (with Lees-Edwards velocity shift), walls reflect or thermalise.
    pub fn apply_boundaries(
        &mut self,
        rng: &mut Rng,
    ) {
        if let Some(h) = self.bx.cell {
            for i in 0..self.n() {
                let mut s = self.bx.to_fractional(&self.pos[i]);
                for c in &mut s {
                    *c = c.floor();
                }
                for a in 0..D {
                    for b in 0..D {
                        self.pos[i][b] -= s[a] * h[a][b];
                    }
                }
            }
            return;
        }
        let l = self.bx.len;
        let shear = D >= 2 && self.bx.shear_rate != 0.0;
        for i in 0..self.n() {
            if shear {
                while self.pos[i][1] >= l[1] {
                    self.pos[i][1] -= l[1];
                    self.pos[i][0] -= self.bx.offset;
                    self.vel[i][0] -= self.bx.shear_rate * l[1];
                }
                while self.pos[i][1] < 0.0 {
                    self.pos[i][1] += l[1];
                    self.pos[i][0] += self.bx.offset;
                    self.vel[i][0] += self.bx.shear_rate * l[1];
                }
            }
            for a in 0..D {
                match self.bx.bc[a] {
                    Boundary::Open => {},
                    Boundary::Periodic => {
                        // the gradient axis was already wrapped above under shear
                        if !(a == 1 && shear) {
                            let w = (self.pos[i][a] / l[a]).floor();
                            self.pos[i][a] -= w * l[a];
                        }
                    },
                    Boundary::Reflecting | Boundary::Thermal(_) => {
                        let mut guard = 0;
                        while (self.pos[i][a] < 0.0 || self.pos[i][a] > l[a]) && guard < 8 {
                            let low = self.pos[i][a] < 0.0;
                            self.pos[i][a] = if low { -self.pos[i][a] } else { 2.0 * l[a] - self.pos[i][a] };
                            if let Boundary::Thermal(t) = self.bx.bc[a] {
                                let s = (t / self.mass[i]).sqrt();
                                let vn = s * (-2.0 * (1.0 - rng.uniform()).ln()).sqrt();
                                self.vel[i][a] = if low { vn } else { -vn };
                                for b in 0..D {
                                    if b != a {
                                        self.vel[i][b] = s * rng.normal();
                                    }
                                }
                            } else {
                                self.vel[i][a] = -self.vel[i][a];
                            }
                            guard += 1;
                        }
                        self.unwrapped[i][a] = self.pos[i][a];
                    },
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Initial configurations
// ---------------------------------------------------------------------------

/// Simple-cubic lattice `n[0] x ... x n[D-1]` with the given spacing (origin offset by half a spacing).
#[must_use]
pub fn lattice_cubic<const D: usize>(
    n: [usize; D],
    spacing: f64,
) -> Vec<[f64; D]> {
    let total: usize = n.iter().product();
    let mut out = Vec::with_capacity(total);
    for idx in 0..total {
        let mut r = idx;
        let mut p = [0.0; D];
        for a in 0..D {
            p[a] = ((r % n[a]) as f64 + 0.5) * spacing;
            r /= n[a];
        }
        out.push(p);
    }
    out
}

/// Face-centred cubic lattice with `n^3` conventional cells of edge `a` (4 atoms per cell).
#[must_use]
pub fn lattice_fcc(
    n: usize,
    a: f64,
) -> Vec<[f64; 3]> {
    let basis = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]];
    let mut out = Vec::new();
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                for b in &basis {
                    out.push([(i as f64 + b[0] + 0.25) * a, (j as f64 + b[1] + 0.25) * a, (k as f64 + b[2] + 0.25) * a]);
                }
            }
        }
    }
    out
}

/// Triangular (hexagonal close-packed) 2-D lattice, `nx` by `ny` sites with nearest-neighbour distance `a`
/// (use even `ny` and box `(nx a, ny a sqrt(3)/2)` for periodicity).
#[must_use]
pub fn lattice_triangular(
    nx: usize,
    ny: usize,
    a: f64,
) -> Vec<[f64; 2]> {
    let mut out = Vec::new();
    for j in 0..ny {
        for i in 0..nx {
            out.push([(i as f64 + 0.5 * (j % 2) as f64 + 0.25) * a, (j as f64 + 0.5) * a * 3.0_f64.sqrt() / 2.0]);
        }
    }
    out
}

/// Random sequential placement of `n` particles with minimum separation `dmin` (`None` if packing fails).
#[must_use]
pub fn random_positions<const D: usize>(
    n: usize,
    bx: &SimBox<D>,
    dmin: f64,
    rng: &mut Rng,
) -> Option<Vec<[f64; D]>> {
    let mut out: Vec<[f64; D]> = Vec::with_capacity(n);
    let mut tries = 0usize;
    while out.len() < n {
        tries += 1;
        if tries > 5000 * n.max(1) {
            return None;
        }
        let mut p = [0.0; D];
        for a in 0..D {
            p[a] = rng.uniform() * bx.len[a];
        }
        let ok = out.iter().all(|q| {
            let mut d = [0.0; D];
            for a in 0..D {
                d[a] = p[a] - q[a];
            }
            bx.min_image(&mut d);
            d.iter().map(|x| x * x).sum::<f64>() >= dmin * dmin
        });
        if ok {
            out.push(p);
        }
    }
    Some(out)
}

// ---------------------------------------------------------------------------
// Force field
// ---------------------------------------------------------------------------

/// One additive non-bonded pair interaction.
#[derive(Clone, Copy, Debug)]
pub enum PairTerm {
    /// `4 eps ((s/r)^12 - (s/r)^6)`.
    LennardJones {
        /// Well depth.
        eps: f64,
        /// Length scale.
        sigma: f64,
    },
    /// `D (1 - exp(-a (r - r0)))^2 - D`.
    Morse {
        /// Well depth.
        d: f64,
        /// Width parameter.
        alpha: f64,
        /// Equilibrium distance.
        r0: f64,
    },
    /// `k q_i q_j erfc(alpha r) / r` (`alpha = 0`: bare Coulomb).
    Coulomb {
        /// Coupling constant.
        k: f64,
        /// Ewald splitting parameter (0 for none).
        alpha: f64,
    },
}

impl PairTerm {
    /// Energy and `-dU/dr / r` at squared distance `r2` for charge product `qq`.
    #[must_use]
    pub fn eval(
        &self,
        r2: f64,
        qq: f64,
    ) -> (f64, f64) {
        match *self {
            Self::LennardJones { eps, sigma } => {
                let s2 = sigma * sigma / r2;
                let s6 = s2 * s2 * s2;
                (4.0 * eps * (s6 * s6 - s6), 24.0 * eps * (2.0 * s6 * s6 - s6) / r2)
            },
            Self::Morse { d, alpha, r0 } => {
                let r = r2.sqrt();
                let e = (-alpha * (r - r0)).exp();
                (d * ((1.0 - e) * (1.0 - e) - 1.0), -2.0 * d * alpha * e * (1.0 - e) / r)
            },
            Self::Coulomb { k, alpha } => {
                let r = r2.sqrt();
                if alpha == 0.0 {
                    (k * qq / r, k * qq / (r * r2))
                } else {
                    let ec = erfc_numerical(alpha * r);
                    let g = 2.0 * alpha / std::f64::consts::PI.sqrt() * (-alpha * alpha * r2).exp();
                    (k * qq * ec / r, k * qq * (ec / r + g) / r2)
                }
            },
        }
    }
}

/// Harmonic bond.
#[derive(Clone, Copy, Debug)]
pub struct Bond {
    /// First atom.
    pub i: usize,
    /// Second atom.
    pub j: usize,
    /// Force constant.
    pub k: f64,
    /// Rest length.
    pub r0: f64,
}

/// Harmonic angle with apex `j`.
#[derive(Clone, Copy, Debug)]
pub struct Angle {
    /// End atom.
    pub i: usize,
    /// Apex atom.
    pub j: usize,
    /// End atom.
    pub k: usize,
    /// Force constant (`U = k/2 (theta - theta0)^2`).
    pub ktheta: f64,
    /// Equilibrium angle.
    pub theta0: f64,
}

/// Periodic dihedral `U = k (1 + cos(n phi - delta))` (3-D only).
#[derive(Clone, Copy, Debug)]
pub struct Dihedral {
    /// Atoms `i-j-k-l`.
    pub atoms: [usize; 4],
    /// Barrier height.
    pub k: f64,
    /// Multiplicity.
    pub n: f64,
    /// Phase.
    pub delta: f64,
}

/// Callback adding extra (e.g. long-range) forces.
pub type ExtraForce<const D: usize> = Box<dyn Fn(&MdSystem<D>, &mut Forces<D>) + Send + Sync>;

/// Complete force field.
pub struct ForceField<const D: usize> {
    /// Pair interactions (summed).
    pub pair: Vec<PairTerm>,
    /// Pair cutoff.
    pub rcut: f64,
    /// Subtract the pair energy at the cutoff (energy shift).
    pub shift: bool,
    /// Shifted-force truncation: energy and force both vanish at the cutoff (smooth dynamics; overrides `shift`).
    pub shift_force: bool,
    /// Harmonic bonds.
    pub bonds: Vec<Bond>,
    /// Harmonic angles.
    pub angles: Vec<Angle>,
    /// Dihedrals.
    pub dihedrals: Vec<Dihedral>,
    /// Exclude bonded (1-2) and angle-ends (1-3) pairs from the non-bonded sum.
    pub exclude_bonded: bool,
    /// Additional force callbacks (Ewald, external fields, ...).
    pub extra: Vec<ExtraForce<D>>,
    /// Additional pairs `(i, j)` excluded from the non-bonded sum (e.g. sites of one rigid
    /// molecule); call [`ForceField::rebuild_exclusions`] after editing.
    pub extra_exclusions: Vec<(usize, usize)>,
    excl: HashSet<(usize, usize)>,
}

impl<const D: usize> ForceField<D> {
    /// New force field with the given pair terms and cutoff.
    #[must_use]
    pub fn new(
        pair: Vec<PairTerm>,
        rcut: f64,
    ) -> Self {
        Self {
            pair,
            rcut,
            shift: false,
            shift_force: false,
            bonds: Vec::new(),
            angles: Vec::new(),
            dihedrals: Vec::new(),
            exclude_bonded: true,
            extra: Vec::new(),
            extra_exclusions: Vec::new(),
            excl: HashSet::new(),
        }
    }

    /// Rebuilds the exclusion table (call after changing bonded terms).
    pub fn rebuild_exclusions(&mut self) {
        self.excl.clear();
        for &(i, j) in &self.extra_exclusions {
            self.excl.insert((i.min(j), i.max(j)));
        }
        if !self.exclude_bonded {
            return;
        }
        for b in &self.bonds {
            self.excl.insert((b.i.min(b.j), b.i.max(b.j)));
        }
        for a in &self.angles {
            self.excl.insert((a.i.min(a.k), a.i.max(a.k)));
        }
    }

    /// Validates the force field for dimension `D`.
    ///
    /// # Errors
    /// Dihedrals need `D = 3`; angles need `D >= 2`.
    pub fn validate(&self) -> Result<(), String> {
        if !self.dihedrals.is_empty() && D != 3 {
            return Err("dihedral terms require D = 3".into());
        }
        if !self.angles.is_empty() && D < 2 {
            return Err("angle terms require D >= 2".into());
        }
        Ok(())
    }

    fn pair_energy_shift(
        &self,
        qq: f64,
    ) -> f64 {
        if !self.shift {
            return 0.0;
        }
        self.pair.iter().map(|t| t.eval(self.rcut * self.rcut, qq).0).sum()
    }
}

/// Result of a force evaluation.
#[derive(Clone, Debug)]
pub struct Forces<const D: usize> {
    /// Force on each particle.
    pub f: Vec<[f64; D]>,
    /// Potential energy.
    pub epot: f64,
    /// Virial tensor `sum r_i (x) f_i` (translation invariant form).
    pub virial: [[f64; D]; D],
}

impl<const D: usize> Forces<D> {
    /// Zero forces for `n` particles.
    #[must_use]
    pub fn zeros(n: usize) -> Self {
        Self { f: vec![[0.0; D]; n], epot: 0.0, virial: [[0.0; D]; D] }
    }
}

/// Cell list + Verlet skin neighbour list.
#[derive(Clone, Debug)]
pub struct NeighbourList<const D: usize> {
    /// Interaction cutoff.
    pub rcut: f64,
    /// Verlet skin.
    pub skin: f64,
    /// Candidate pairs `i < j` within `rcut + skin` at build time.
    pub pairs: Vec<(usize, usize)>,
    ref_pos: Vec<[f64; D]>,
    ref_cell: [[f64; D]; D],
    /// Number of rebuilds so far.
    pub builds: usize,
}

impl<const D: usize> NeighbourList<D> {
    /// Creates an empty list.
    #[must_use]
    pub fn new(
        rcut: f64,
        skin: f64,
    ) -> Self {
        Self { rcut, skin, pairs: Vec::new(), ref_pos: Vec::new(), ref_cell: [[0.0; D]; D], builds: 0 }
    }

    /// Whether the list is stale: a particle moved (beyond the affine cell deformation) so that
    /// `2 max|dx| + |strain| (r_cut + skin) > skin`, the box shears, or the list was never built.
    #[must_use]
    pub fn needs_rebuild(
        &self,
        sys: &MdSystem<D>,
    ) -> bool {
        if self.ref_pos.len() != sys.n() || sys.bx.shear_rate != 0.0 || self.builds == 0 {
            return true;
        }
        let now = sys.bx.matrix();
        let mut strain2 = 0.0;
        let affine = now != self.ref_cell;
        let mut map = [[0.0; D]; D];
        if affine {
            let Some(inv) = mat_inverse(&self.ref_cell) else { return true };
            map = mat_mul(&inv, &now);
            for a in 0..D {
                for b in 0..D {
                    let id = if a == b { 1.0 } else { 0.0 };
                    strain2 += (map[a][b] - id).powi(2);
                }
            }
        }
        let rl = self.rcut + self.skin;
        let strain = strain2.sqrt();
        if strain * rl > 0.5 * self.skin {
            return true;
        }
        let tol = 0.5 * (self.skin - strain * rl);
        let lim = tol * tol;
        let ref_inv = if affine { mat_inverse(&self.ref_cell) } else { None };
        for i in 0..sys.n() {
            let mut pred = self.ref_pos[i];
            if let Some(inv) = &ref_inv {
                // x_pred = (x_ref h_ref^{-1}) h_now
                let mut s = [0.0; D];
                for a in 0..D {
                    for b in 0..D {
                        s[a] += self.ref_pos[i][b] * inv[b][a];
                    }
                }
                pred = [0.0; D];
                for a in 0..D {
                    for b in 0..D {
                        pred[b] += s[a] * now[a][b];
                    }
                }
            }
            let mut d = [0.0; D];
            for a in 0..D {
                d[a] = sys.pos[i][a] - pred[a];
            }
            sys.bx.min_image(&mut d);
            if d.iter().map(|x| x * x).sum::<f64>() > lim {
                return true;
            }
        }
        let _ = map;
        false
    }

    /// Rebuilds the list with a cell grid (brute force under Lees-Edwards shear).
    pub fn build(
        &mut self,
        sys: &MdSystem<D>,
    ) {
        self.builds += 1;
        let n = sys.n();
        let rl = self.rcut + self.skin;
        let rl2 = rl * rl;
        self.pairs.clear();
        self.ref_pos = sys.pos.clone();
        self.ref_cell = sys.bx.matrix();
        let close = |i: usize, j: usize| -> bool {
            let mut d = [0.0; D];
            for a in 0..D {
                d[a] = sys.pos[i][a] - sys.pos[j][a];
            }
            sys.bx.min_image(&mut d);
            d.iter().map(|x| x * x).sum::<f64>() < rl2
        };
        if sys.bx.shear_rate != 0.0 {
            for i in 0..n {
                for j in i + 1..n {
                    if close(i, j) {
                        self.pairs.push((i, j));
                    }
                }
            }
            return;
        }
        // grid; open axes extend over the particle bounding range
        let tri = sys.bx.cell.is_some();
        let pc: Vec<[f64; D]> = if tri {
            sys.pos
                .iter()
                .map(|p| {
                    let mut s = sys.bx.to_fractional(p);
                    for c in &mut s {
                        *c -= c.floor();
                    }
                    s
                })
                .collect()
        } else {
            sys.pos.clone()
        };
        let mut lo = [0.0; D];
        let mut ext = if tri { [1.0; D] } else { sys.bx.len };
        for a in 0..D {
            if sys.bx.bc[a] == Boundary::Open {
                let (mn, mx) = sys.pos.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(m, x), p| (m.min(p[a]), x.max(p[a])));
                lo[a] = mn;
                ext[a] = (mx - mn).max(1e-12);
            }
        }
        let mut nc = [1usize; D];
        for a in 0..D {
            let width = if tri { sys.bx.height(a) } else { ext[a] };
            let c = (width / rl).floor() as usize;
            nc[a] = if c >= 3 { c } else { 1 };
        }
        let total: usize = nc.iter().product();
        let cell_of = |p: &[f64; D]| -> usize {
            let mut idx = 0usize;
            let mut stride = 1usize;
            for a in 0..D {
                let mut c = (((p[a] - lo[a]) / ext[a]) * nc[a] as f64).floor() as i64;
                c = c.clamp(0, nc[a] as i64 - 1);
                idx += c as usize * stride;
                stride *= nc[a];
            }
            idx
        };
        let mut cells: Vec<Vec<usize>> = vec![Vec::new(); total];
        for i in 0..n {
            cells[cell_of(&pc[i])].push(i);
        }
        let offsets: Vec<[i64; D]> = {
            let count = 3usize.pow(D as u32);
            (0..count)
                .map(|mut r| {
                    let mut o = [0i64; D];
                    for a in 0..D {
                        o[a] = (r % 3) as i64 - 1;
                        r /= 3;
                    }
                    o
                })
                .collect()
        };
        for c in 0..total {
            let mut coords = [0i64; D];
            let mut r = c;
            for a in 0..D {
                coords[a] = (r % nc[a]) as i64;
                r /= nc[a];
            }
            for off in &offsets {
                let mut nb = 0usize;
                let mut stride = 1usize;
                let mut valid = true;
                for a in 0..D {
                    let mut x = coords[a] + off[a];
                    if nc[a] == 1 {
                        if off[a] != 0 {
                            valid = false;
                            break;
                        }
                        x = 0;
                    } else if x < 0 || x >= nc[a] as i64 {
                        if sys.bx.bc[a] == Boundary::Periodic {
                            x = x.rem_euclid(nc[a] as i64);
                        } else {
                            valid = false;
                            break;
                        }
                    }
                    nb += x as usize * stride;
                    stride *= nc[a];
                }
                if !valid || nb < c {
                    continue;
                }
                for &i in &cells[c] {
                    for &j in &cells[nb] {
                        if (nb == c && j <= i) || !close(i, j) {
                            continue;
                        }
                        self.pairs.push((i.min(j), i.max(j)));
                    }
                }
            }
        }
    }
}

fn dot<const D: usize>(
    a: &[f64; D],
    b: &[f64; D],
) -> f64 {
    let mut s = 0.0;
    for k in 0..D {
        s += a[k] * b[k];
    }
    s
}

fn add_virial<const D: usize>(
    out: &mut Forces<D>,
    r: &[f64; D],
    f: &[f64; D],
) {
    for a in 0..D {
        for b in 0..D {
            out.virial[a][b] += r[a] * f[b];
        }
    }
}

/// Evaluates all forces, rebuilding the neighbour list when needed.
///
/// # Errors
/// Propagates [`ForceField::validate`] errors.
pub fn compute_forces<const D: usize>(
    sys: &MdSystem<D>,
    ff: &ForceField<D>,
    nl: &mut NeighbourList<D>,
) -> Result<Forces<D>, String> {
    ff.validate()?;
    if nl.needs_rebuild(sys) {
        nl.build(sys);
    }
    let mut out = Forces::zeros(sys.n());
    let rc2 = ff.rcut * ff.rcut;
    for &(i, j) in &nl.pairs {
        if !ff.excl.is_empty() && ff.excl.contains(&(i, j)) {
            continue;
        }
        let mut d = [0.0; D];
        for a in 0..D {
            d[a] = sys.pos[i][a] - sys.pos[j][a];
        }
        sys.bx.min_image(&mut d);
        let r2 = dot(&d, &d);
        if r2 >= rc2 || r2 == 0.0 {
            continue;
        }
        let qq = sys.charge[i] * sys.charge[j];
        let mut fr = 0.0;
        let mut u = -ff.pair_energy_shift(qq);
        for t in &ff.pair {
            let (e, f) = t.eval(r2, qq);
            u += e;
            fr += f;
            if ff.shift_force {
                let (ec, fc) = t.eval(rc2, qq);
                let big_fc = fc * ff.rcut;
                let r = r2.sqrt();
                u += -ec + big_fc * (r - ff.rcut);
                fr -= big_fc / r;
            }
        }
        out.epot += u;
        let mut fv = [0.0; D];
        for a in 0..D {
            fv[a] = fr * d[a];
            out.f[i][a] += fv[a];
            out.f[j][a] -= fv[a];
        }
        add_virial(&mut out, &d, &fv);
    }
    bonded_forces(sys, ff, &mut out);
    for e in &ff.extra {
        e(sys, &mut out);
    }
    Ok(out)
}

fn sep<const D: usize>(
    sys: &MdSystem<D>,
    i: usize,
    j: usize,
) -> [f64; D] {
    let mut d = [0.0; D];
    for a in 0..D {
        d[a] = sys.pos[i][a] - sys.pos[j][a];
    }
    sys.bx.min_image(&mut d);
    d
}

fn bonded_forces<const D: usize>(
    sys: &MdSystem<D>,
    ff: &ForceField<D>,
    out: &mut Forces<D>,
) {
    for b in &ff.bonds {
        let d = sep(sys, b.i, b.j);
        let r = dot(&d, &d).sqrt();
        if r == 0.0 {
            continue;
        }
        out.epot += 0.5 * b.k * (r - b.r0).powi(2);
        let fr = -b.k * (r - b.r0) / r;
        let mut fv = [0.0; D];
        for a in 0..D {
            fv[a] = fr * d[a];
            out.f[b.i][a] += fv[a];
            out.f[b.j][a] -= fv[a];
        }
        add_virial(out, &d, &fv);
    }
    for g in &ff.angles {
        let r1 = sep(sys, g.i, g.j);
        let r2 = sep(sys, g.k, g.j);
        let (n1, n2) = (dot(&r1, &r1).sqrt(), dot(&r2, &r2).sqrt());
        if n1 == 0.0 || n2 == 0.0 {
            continue;
        }
        let c = (dot(&r1, &r2) / (n1 * n2)).clamp(-1.0, 1.0);
        let s = (1.0 - c * c).sqrt().max(1e-10);
        let th = c.acos();
        out.epot += 0.5 * g.ktheta * (th - g.theta0).powi(2);
        let pref = g.ktheta * (th - g.theta0) / s; // -dU/dtheta = -ktheta (th - th0); dtheta/dx = -1/s dc/dx
        let mut fi = [0.0; D];
        let mut fk = [0.0; D];
        let mut fj = [0.0; D];
        for a in 0..D {
            let dc_di = r2[a] / (n1 * n2) - c * r1[a] / (n1 * n1);
            let dc_dk = r1[a] / (n1 * n2) - c * r2[a] / (n2 * n2);
            fi[a] = pref * dc_di;
            fk[a] = pref * dc_dk;
            fj[a] = -fi[a] - fk[a];
            out.f[g.i][a] += fi[a];
            out.f[g.k][a] += fk[a];
            out.f[g.j][a] += fj[a];
        }
        add_virial(out, &r1, &fi);
        add_virial(out, &r2, &fk);
    }
    for dh in &ff.dihedrals {
        if D != 3 {
            continue;
        }
        let [i, j, k, l] = dh.atoms;
        let to3 = |v: [f64; D]| [v[0], v[1], v[2]];
        let rij = to3(sep(sys, i, j));
        let rkj = to3(sep(sys, k, j));
        let rkl = to3(sep(sys, k, l));
        let cross = |a: [f64; 3], b: [f64; 3]| [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
        let d3 = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        let m = cross(rij, rkj);
        let nn = cross(rkj, rkl);
        let (m2, n2) = (d3(m, m), d3(nn, nn));
        let rkj_len = d3(rkj, rkj).sqrt();
        if m2 < 1e-20 || n2 < 1e-20 {
            continue;
        }
        let cphi = (d3(m, nn) / (m2 * n2).sqrt()).clamp(-1.0, 1.0);
        let mut phi = cphi.acos();
        if d3(rij, nn) < 0.0 {
            phi = -phi;
        }
        out.epot += dh.k * (1.0 + (dh.n * phi - dh.delta).cos());
        let dudphi = -dh.k * dh.n * (dh.n * phi - dh.delta).sin();
        let mut dphi_di = [0.0; 3];
        let mut dphi_dl = [0.0; 3];
        for a in 0..3 {
            dphi_di[a] = rkj_len / m2 * m[a];
            dphi_dl[a] = -rkj_len / n2 * nn[a];
        }
        let p = d3(rij, rkj) / (rkj_len * rkj_len);
        let q = d3(rkl, rkj) / (rkj_len * rkj_len);
        let mut fi = [0.0; 3];
        let mut fj = [0.0; 3];
        let mut fk = [0.0; 3];
        let mut fl = [0.0; 3];
        for a in 0..3 {
            let dj = (p - 1.0) * dphi_di[a] - q * dphi_dl[a];
            let dk = (q - 1.0) * dphi_dl[a] - p * dphi_di[a];
            fi[a] = -dudphi * dphi_di[a];
            fj[a] = -dudphi * dj;
            fk[a] = -dudphi * dk;
            fl[a] = -dudphi * dphi_dl[a];
        }
        for a in 0..3 {
            out.f[i][a] += fi[a];
            out.f[j][a] += fj[a];
            out.f[k][a] += fk[a];
            out.f[l][a] += fl[a];
        }
        // virial relative to atom j (forces sum to zero)
        let mut ri = [0.0; D];
        let mut rk_ = [0.0; D];
        let mut rl_ = [0.0; D];
        for a in 0..3 {
            ri[a] = rij[a];
            rk_[a] = rkj[a];
            rl_[a] = rkj[a] - rkl[a];
        }
        let to_d = |v: [f64; 3]| {
            let mut o = [0.0; D];
            for a in 0..3 {
                o[a] = v[a];
            }
            o
        };
        add_virial(out, &ri, &to_d(fi));
        add_virial(out, &rk_, &to_d(fk));
        add_virial(out, &rl_, &to_d(fl));
    }
}

// ---------------------------------------------------------------------------
// Thermostats, barostats, constraints, simulation driver
// ---------------------------------------------------------------------------

/// Thermostat selection.
#[derive(Clone, Debug)]
pub enum Thermostat {
    /// Microcanonical.
    None,
    /// Berendsen weak coupling (`t0`, relaxation time `tau`).
    Berendsen {
        /// Target temperature.
        t0: f64,
        /// Coupling time.
        tau: f64,
    },
    /// Langevin dynamics integrated with BAOAB (friction `gamma`).
    Langevin {
        /// Target temperature.
        t0: f64,
        /// Friction coefficient.
        gamma: f64,
    },
    /// Nose-Hoover chain of the given length.
    NoseHooverChain(NoseHooverChain),
}

/// Nose-Hoover chain state (Martyna-Klein-Tuckerman).
#[derive(Clone, Debug)]
pub struct NoseHooverChain {
    /// Target temperature.
    pub t0: f64,
    /// Thermostat masses `Q_j`.
    pub q: Vec<f64>,
    /// Chain positions.
    pub xi: Vec<f64>,
    /// Chain velocities.
    pub vxi: Vec<f64>,
    /// Degrees of freedom coupled to the chain.
    pub nf: f64,
}

impl NoseHooverChain {
    /// Chain of `m` thermostats with relaxation time `tau`.
    #[must_use]
    pub fn new(
        t0: f64,
        tau: f64,
        m: usize,
        nf: f64,
    ) -> Self {
        let mut q = vec![t0 * tau * tau; m];
        q[0] *= nf;
        Self { t0, q, xi: vec![0.0; m], vxi: vec![0.0; m], nf }
    }

    /// Propagates the chain for `h / 2` and returns the velocity scale factor
    /// (`ke2` is `sum m v^2` before the call).
    pub fn propagate(
        &mut self,
        ke2_in: f64,
        h: f64,
    ) -> f64 {
        let m = self.q.len();
        let kt = self.t0;
        let w1 = 1.0 / (2.0 - 2.0_f64.cbrt());
        let ws = [w1, 1.0 - 2.0 * w1, w1];
        let mut ke2 = ke2_in;
        let mut scale = 1.0;
        let mut g = vec![0.0; m];
        g[0] = (ke2 - self.nf * kt) / self.q[0];
        for j in 1..m {
            g[j] = (self.q[j - 1] * self.vxi[j - 1] * self.vxi[j - 1] - kt) / self.q[j];
        }
        for w in ws {
            let wdt = w * h;
            self.vxi[m - 1] += g[m - 1] * wdt / 4.0;
            for j in (0..m - 1).rev() {
                let aa = (-self.vxi[j + 1] * wdt / 8.0).exp();
                self.vxi[j] = self.vxi[j] * aa * aa + wdt / 4.0 * g[j] * aa;
            }
            let aa = (-self.vxi[0] * wdt / 2.0).exp();
            scale *= aa;
            ke2 *= aa * aa;
            g[0] = (ke2 - self.nf * kt) / self.q[0];
            for j in 0..m {
                self.xi[j] += self.vxi[j] * wdt / 2.0;
            }
            for j in 0..m - 1 {
                let aa = (-self.vxi[j + 1] * wdt / 8.0).exp();
                self.vxi[j] = self.vxi[j] * aa * aa + wdt / 4.0 * g[j] * aa;
                g[j + 1] = (self.q[j] * self.vxi[j] * self.vxi[j] - kt) / self.q[j + 1];
            }
            self.vxi[m - 1] += g[m - 1] * wdt / 4.0;
        }
        scale
    }

    /// Thermostat contribution to the conserved quantity.
    #[must_use]
    pub fn energy(&self) -> f64 {
        let mut e = 0.0;
        for j in 0..self.q.len() {
            e += 0.5 * self.q[j] * self.vxi[j] * self.vxi[j];
            e += if j == 0 { self.nf * self.t0 * self.xi[0] } else { self.t0 * self.xi[j] };
        }
        e
    }
}

/// Barostat selection.
#[derive(Clone, Debug)]
pub enum Barostat {
    /// Constant volume.
    None,
    /// Berendsen weak coupling with isothermal compressibility `beta`.
    Berendsen {
        /// Target pressure.
        p0: f64,
        /// Coupling time.
        tau: f64,
        /// Compressibility.
        beta: f64,
    },
    /// Flexible-cell Parrinello-Rahman / Martyna-Tobias-Klein barostat (see [`FlexibleBarostat`]).
    Flexible(FlexibleBarostat),
    /// Isotropic MTK barostat with piston mass `w` (NPH-type dynamics).
    Mtk {
        /// Target pressure.
        p0: f64,
        /// Piston mass.
        w: f64,
        /// Piston momentum `p_eps`.
        p_eps: f64,
    },
}

/// Distance constraint between two atoms.
#[derive(Clone, Copy, Debug)]
pub struct Constraint {
    /// First atom.
    pub i: usize,
    /// Second atom.
    pub j: usize,
    /// Fixed distance.
    pub d: f64,
}

/// Molecular dynamics driver.
pub struct MdSimulation<const D: usize> {
    /// Particle system.
    pub sys: MdSystem<D>,
    /// Force field.
    pub ff: ForceField<D>,
    /// Neighbour list.
    pub nl: NeighbourList<D>,
    /// Current forces.
    pub forces: Forces<D>,
    /// Thermostat.
    pub thermostat: Thermostat,
    /// Barostat.
    pub barostat: Barostat,
    /// Distance constraints.
    pub constraints: Vec<Constraint>,
    /// Time step.
    pub dt: f64,
    /// Random number generator (Langevin, thermal walls).
    pub rng: Rng,
    /// Virial `sum r_ij (x) f_c` of the constraint forces from the last SHAKE sweep
    /// (added to the force virial in the pressure of constrained systems).
    pub constraint_virial: [[f64; D]; D],
}

impl<const D: usize> MdSimulation<D> {
    /// Creates a simulation and evaluates the initial forces.
    ///
    /// # Errors
    /// Propagates force-field validation errors.
    pub fn new(
        sys: MdSystem<D>,
        mut ff: ForceField<D>,
        skin: f64,
        dt: f64,
        seed: u64,
    ) -> Result<Self, String> {
        ff.rebuild_exclusions();
        let mut nl = NeighbourList::new(ff.rcut, skin);
        let forces = compute_forces(&sys, &ff, &mut nl)?;
        Ok(Self { sys, ff, nl, forces, thermostat: Thermostat::None, barostat: Barostat::None, constraints: Vec::new(), dt, rng: Rng::new(seed), constraint_virial: [[0.0; D]; D] })
    }

    /// Degrees of freedom (`D N - D - n_constraints`).
    #[must_use]
    pub fn dof(&self) -> f64 {
        (D * self.sys.n()).saturating_sub(D + self.constraints.len()) as f64
    }

    /// Instantaneous temperature.
    #[must_use]
    pub fn temperature(&self) -> f64 {
        self.sys.temperature(self.dof())
    }

    /// Total energy.
    #[must_use]
    pub fn total_energy(&self) -> f64 {
        self.sys.kinetic_energy() + self.forces.epot
    }

    /// Pressure tensor `(sum m v v + virial) / V`.
    #[must_use]
    pub fn pressure_tensor(&self) -> [[f64; D]; D] {
        let v = self.sys.bx.volume();
        let mut p = [[0.0; D]; D];
        for i in 0..self.sys.n() {
            for a in 0..D {
                for b in 0..D {
                    p[a][b] += self.sys.mass[i] * self.sys.vel[i][a] * self.sys.vel[i][b];
                }
            }
        }
        for a in 0..D {
            for b in 0..D {
                p[a][b] = (p[a][b] + self.forces.virial[a][b] + self.constraint_virial[a][b]) / v;
            }
        }
        p
    }

    /// Scalar pressure (trace of the pressure tensor over `D`).
    #[must_use]
    pub fn pressure(&self) -> f64 {
        let p = self.pressure_tensor();
        (0..D).map(|a| p[a][a]).sum::<f64>() / D as f64
    }

    /// Conserved quantity of the extended system (energy plus thermostat and MTK piston terms).
    #[must_use]
    pub fn conserved_energy(&self) -> f64 {
        let mut e = self.total_energy();
        if let Thermostat::NoseHooverChain(c) = &self.thermostat {
            e += c.energy();
        }
        if let Barostat::Mtk { p0, w, p_eps } = &self.barostat {
            e += p_eps * p_eps / (2.0 * w) + p0 * self.sys.bx.volume();
        }
        if let Barostat::Flexible(b) = &self.barostat {
            e += b.energy(self.sys.bx.volume());
        }
        e
    }

    fn recompute(&mut self) -> Result<(), String> {
        self.forces = compute_forces(&self.sys, &self.ff, &mut self.nl)?;
        Ok(())
    }

    fn half_kick(
        &mut self,
        h: f64,
    ) {
        for i in 0..self.sys.n() {
            for a in 0..D {
                self.sys.vel[i][a] += 0.5 * h * self.forces.f[i][a] / self.sys.mass[i];
            }
        }
    }

    fn drift(
        &mut self,
        h: f64,
    ) {
        let old = self.sys.pos.clone();
        for i in 0..self.sys.n() {
            for a in 0..D {
                let dx = h * self.sys.vel[i][a];
                self.sys.pos[i][a] += dx;
                self.sys.unwrapped[i][a] += dx;
            }
        }
        if !self.constraints.is_empty() {
            self.shake(&old, h);
        }
        if self.sys.bx.shear_rate != 0.0 {
            self.sys.bx.offset += self.sys.bx.shear_rate * self.sys.bx.len[1] * h;
            self.sys.bx.offset -= self.sys.bx.len[0] * (self.sys.bx.offset / self.sys.bx.len[0]).floor();
        }
        self.sys.apply_boundaries(&mut self.rng);
    }

    fn shake(
        &mut self,
        old: &[[f64; D]],
        h: f64,
    ) {
        self.constraint_virial = [[0.0; D]; D];
        for _ in 0..500 {
            let mut worst = 0.0_f64;
            for c in self.constraints.clone() {
                let (i, j) = (c.i, c.j);
                let mut s = [0.0; D];
                let mut r_old = [0.0; D];
                for a in 0..D {
                    s[a] = self.sys.pos[i][a] - self.sys.pos[j][a];
                    r_old[a] = old[i][a] - old[j][a];
                }
                self.sys.bx.min_image(&mut s);
                self.sys.bx.min_image(&mut r_old);
                let diff = dot(&s, &s) - c.d * c.d;
                worst = worst.max(diff.abs() / (c.d * c.d));
                let (mi, mj) = (self.sys.mass[i], self.sys.mass[j]);
                let g = diff / (2.0 * dot(&s, &r_old) * (1.0 / mi + 1.0 / mj));
                for a in 0..D {
                    for b in 0..D {
                        self.constraint_virial[a][b] -= 2.0 * g * r_old[a] * r_old[b] / (h * h);
                    }
                }
                for a in 0..D {
                    let di = -g * r_old[a] / mi;
                    let dj = g * r_old[a] / mj;
                    self.sys.pos[i][a] += di;
                    self.sys.pos[j][a] += dj;
                    self.sys.unwrapped[i][a] += di;
                    self.sys.unwrapped[j][a] += dj;
                    self.sys.vel[i][a] += di / h;
                    self.sys.vel[j][a] += dj / h;
                }
            }
            if worst < 1e-13 {
                break;
            }
        }
    }

    /// RATTLE velocity projection. With a cell velocity gradient `vg` the constraint
    /// condition is `r . (v_i - v_j) = -r . (vg r)` (the bond is carried by the cell flow).
    fn rattle_velocities(
        &mut self,
        vg: Option<&[[f64; D]; D]>,
    ) {
        for _ in 0..500 {
            let mut worst = 0.0_f64;
            for c in self.constraints.clone() {
                let (i, j) = (c.i, c.j);
                let mut r = [0.0; D];
                let mut v = [0.0; D];
                for a in 0..D {
                    r[a] = self.sys.pos[i][a] - self.sys.pos[j][a];
                    v[a] = self.sys.vel[i][a] - self.sys.vel[j][a];
                }
                self.sys.bx.min_image(&mut r);
                let (mi, mj) = (self.sys.mass[i], self.sys.mass[j]);
                let mut target = 0.0;
                if let Some(g) = vg {
                    for a in 0..D {
                        for b in 0..D {
                            target -= r[a] * g[a][b] * r[b];
                        }
                    }
                }
                let rv = dot(&r, &v) - target;
                worst = worst.max(rv.abs() / (c.d * c.d));
                let k = rv / (dot(&r, &r) * (1.0 / mi + 1.0 / mj));
                for a in 0..D {
                    self.sys.vel[i][a] -= k * r[a] / mi;
                    self.sys.vel[j][a] += k * r[a] / mj;
                }
            }
            if worst < 1e-13 {
                break;
            }
        }
    }

    fn scale_velocities(
        &mut self,
        f: f64,
    ) {
        for v in &mut self.sys.vel {
            for c in v.iter_mut() {
                *c *= f;
            }
        }
    }

    fn scale_box(
        &mut self,
        mu: f64,
    ) {
        for a in 0..D {
            self.sys.bx.len[a] *= mu;
        }
        for p in &mut self.sys.pos {
            for c in p.iter_mut() {
                *c *= mu;
            }
        }
        for p in &mut self.sys.unwrapped {
            for c in p.iter_mut() {
                *c *= mu;
            }
        }
    }

    /// Advances one time step.
    ///
    /// # Errors
    /// Fails for unsupported combinations (isotropic MTK with a thermostat or constraints, Langevin with a flexible cell) or force-field errors.
    pub fn step(&mut self) -> Result<(), String> {
        let h = self.dt;
        let n = self.sys.n();
        if matches!(self.barostat, Barostat::Flexible(_)) {
            return self.step_flexible();
        }
        if let Thermostat::Langevin { t0, gamma } = self.thermostat {
            // BAOAB; with constraints every velocity update (B, O) is followed by a RATTLE
            // projection and every drift (A) by SHAKE, so the constrained flow is sampled
            // with the tangent-space Gibbs measure.
            let constrained = !self.constraints.is_empty();
            self.half_kick(h);
            if constrained {
                self.rattle_velocities(None);
            }
            self.drift(0.5 * h);
            let c1 = (-gamma * h).exp();
            let c2 = (1.0 - c1 * c1).sqrt();
            for i in 0..n {
                let s = (t0 / self.sys.mass[i]).sqrt();
                for a in 0..D {
                    self.sys.vel[i][a] = c1 * self.sys.vel[i][a] + c2 * s * self.rng.normal();
                }
            }
            if constrained {
                self.rattle_velocities(None);
            }
            self.drift(0.5 * h);
            self.recompute()?;
            self.half_kick(h);
            if constrained {
                self.rattle_velocities(None);
            }
            self.sys.time += h;
            return Ok(());
        }
        if matches!(self.barostat, Barostat::Mtk { .. }) {
            if !matches!(self.thermostat, Thermostat::None) {
                return Err("the MTK barostat is NPH-type and cannot be combined with a thermostat".into());
            }
            return self.step_mtk();
        }
        if let Thermostat::NoseHooverChain(_) = &self.thermostat {
            let ke2 = 2.0 * self.sys.kinetic_energy();
            let s = if let Thermostat::NoseHooverChain(c) = &mut self.thermostat { c.propagate(ke2, h) } else { 1.0 };
            self.scale_velocities(s);
        }
        self.half_kick(h);
        self.drift(h);
        self.recompute()?;
        self.half_kick(h);
        if !self.constraints.is_empty() {
            self.rattle_velocities(None);
        }
        if let Thermostat::NoseHooverChain(_) = &self.thermostat {
            let ke2 = 2.0 * self.sys.kinetic_energy();
            let s = if let Thermostat::NoseHooverChain(c) = &mut self.thermostat { c.propagate(ke2, h) } else { 1.0 };
            self.scale_velocities(s);
        }
        if let Thermostat::Berendsen { t0, tau } = self.thermostat {
            let t = self.temperature();
            if t > 0.0 {
                let lam = (1.0 + h / tau * (t0 / t - 1.0)).max(0.0).sqrt();
                self.scale_velocities(lam);
            }
        }
        if let Barostat::Berendsen { p0, tau, beta } = self.barostat {
            let p = self.pressure();
            let mu = (1.0 - beta * h / tau * (p0 - p)).max(0.1).powf(1.0 / D as f64);
            self.scale_box(mu);
            self.recompute()?;
        }
        self.sys.time += h;
        Ok(())
    }

    fn step_mtk(&mut self) -> Result<(), String> {
        if !self.constraints.is_empty() {
            return Err("the MTK barostat cannot be combined with constraints".into());
        }
        let h = self.dt;
        let d = D as f64;
        let nf = (D * self.sys.n()) as f64;
        let Barostat::Mtk { p0, w, .. } = self.barostat else { return Ok(()) };
        let g_eps = |s: &Self| -> f64 {
            let ke2 = 2.0 * s.sys.kinetic_energy();
            let vir: f64 = (0..D).map(|a| s.forces.virial[a][a]).sum();
            ke2 * (1.0 + d / nf) + vir - d * s.sys.bx.volume() * p0
        };
        let set_p = |s: &mut Self, dp: f64| {
            if let Barostat::Mtk { p_eps, .. } = &mut s.barostat {
                *p_eps += dp;
            }
        };
        let get_p = |s: &Self| -> f64 {
            if let Barostat::Mtk { p_eps, .. } = &s.barostat { *p_eps } else { 0.0 }
        };
        // first half
        let g = g_eps(self);
        set_p(self, 0.5 * h * g);
        let ve = get_p(self) / w;
        let f = (-(1.0 + d / nf) * ve * 0.5 * h).exp();
        self.scale_velocities(f);
        self.half_kick(h);
        // drift with box scaling
        let ve = get_p(self) / w;
        let x = ve * 0.5 * h;
        let sinhc = if x.abs() < 1e-8 { 1.0 + x * x / 6.0 } else { x.sinh() / x };
        let (ex, ex2) = ((ve * h).exp(), x.exp());
        for i in 0..self.sys.n() {
            for a in 0..D {
                let nv = self.sys.pos[i][a] * ex + self.sys.vel[i][a] * ex2 * sinhc * h;
                let dx = nv - self.sys.pos[i][a] * ex;
                self.sys.unwrapped[i][a] = self.sys.unwrapped[i][a] * ex + dx;
                self.sys.pos[i][a] = nv;
            }
        }
        for a in 0..D {
            self.sys.bx.len[a] *= ex;
        }
        self.sys.apply_boundaries(&mut self.rng);
        self.recompute()?;
        self.half_kick(h);
        let ve = get_p(self) / w;
        let f = (-(1.0 + d / nf) * ve * 0.5 * h).exp();
        self.scale_velocities(f);
        let g = g_eps(self);
        set_p(self, 0.5 * h * g);
        self.sys.time += h;
        Ok(())
    }

    /// Runs `n` steps.
    ///
    /// # Errors
    /// Propagates step errors.
    pub fn run(
        &mut self,
        n: usize,
    ) -> Result<(), String> {
        for _ in 0..n {
            self.step()?;
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Observables
// ---------------------------------------------------------------------------

fn unit_ball_volume(d: usize) -> f64 {
    match d {
        1 => 2.0,
        2 => std::f64::consts::PI,
        3 => 4.0 / 3.0 * std::f64::consts::PI,
        _ => {
            let g = crate::kernels::special::gamma_numerical(d as f64 / 2.0 + 1.0);
            std::f64::consts::PI.powf(d as f64 / 2.0) / g
        },
    }
}

/// Radial distribution function accumulator (any dimension, periodic box, `rmax <= L/2`).
#[derive(Clone, Debug)]
pub struct Rdf<const D: usize> {
    /// Maximum distance.
    pub rmax: f64,
    /// Histogram of pair counts.
    pub hist: Vec<f64>,
    frames: usize,
    n: usize,
    volume: f64,
}

impl<const D: usize> Rdf<D> {
    /// New accumulator with `bins` bins up to `rmax`.
    #[must_use]
    pub fn new(
        bins: usize,
        rmax: f64,
    ) -> Self {
        Self { rmax, hist: vec![0.0; bins], frames: 0, n: 0, volume: 0.0 }
    }

    /// Adds one configuration.
    pub fn add_frame(
        &mut self,
        sys: &MdSystem<D>,
    ) {
        let nb = self.hist.len();
        let dr = self.rmax / nb as f64;
        for i in 0..sys.n() {
            for j in i + 1..sys.n() {
                let d = sep(sys, i, j);
                let r = dot(&d, &d).sqrt();
                if r < self.rmax {
                    self.hist[(r / dr) as usize] += 2.0;
                }
            }
        }
        self.frames += 1;
        self.n = sys.n();
        self.volume = sys.bx.volume();
    }

    /// Returns `(r, g(r))` pairs.
    #[must_use]
    pub fn result(&self) -> Vec<(f64, f64)> {
        let nb = self.hist.len();
        let dr = self.rmax / nb as f64;
        let rho = self.n as f64 / self.volume;
        let vb = unit_ball_volume(D);
        (0..nb)
            .map(|k| {
                let (r1, r2) = (k as f64 * dr, (k + 1) as f64 * dr);
                let shell = vb * (r2.powi(D as i32) - r1.powi(D as i32));
                (0.5 * (r1 + r2), self.hist[k] / (self.frames as f64 * self.n as f64 * rho * shell))
            })
            .collect()
    }
}

/// Mean square displacement from unwrapped positions relative to a reference snapshot.
#[must_use]
pub fn msd<const D: usize>(
    sys: &MdSystem<D>,
    reference: &[[f64; D]],
) -> f64 {
    let mut s = 0.0;
    for i in 0..sys.n() {
        for a in 0..D {
            s += (sys.unwrapped[i][a] - reference[i][a]).powi(2);
        }
    }
    s / sys.n() as f64
}

/// Velocity autocorrelation `<v(0) . v(k dt)>` averaged over time origins, for `lags` lags (frames of per-particle velocities).
#[must_use]
pub fn vacf<const D: usize>(
    frames: &[Vec<[f64; D]>],
    lags: usize,
) -> Vec<f64> {
    let mut out = vec![0.0; lags];
    for lag in 0..lags {
        let mut s = 0.0;
        let mut cnt = 0.0;
        for t0 in 0..frames.len().saturating_sub(lag) {
            for (a, b) in frames[t0].iter().zip(&frames[t0 + lag]) {
                s += dot(a, b);
                cnt += 1.0;
            }
        }
        out[lag] = s / cnt;
    }
    out
}

/// Least-squares slope of the energy against time (the energy drift rate).
#[must_use]
pub fn energy_drift_slope(
    times: &[f64],
    energies: &[f64],
) -> f64 {
    let n = times.len() as f64;
    let (mt, me) = (times.iter().sum::<f64>() / n, energies.iter().sum::<f64>() / n);
    let sxy: f64 = times.iter().zip(energies).map(|(t, e)| (t - mt) * (e - me)).sum();
    let sxx: f64 = times.iter().map(|t| (t - mt).powi(2)).sum();
    sxy / sxx
}
