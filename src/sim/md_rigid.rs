//! # Rigid-molecule molecular dynamics (3-D)
//!
//! Molecules made of several interaction sites move as rigid bodies: the centre of mass
//! translates under the total force, and the orientation (unit quaternion) and body-frame
//! angular momentum evolve under the total torque with the symplectic rotor splittings of
//! [`crate::sim::symplectic::RigidBody`] (NO_SQUISH or exact symmetric-top flows, composable to
//! any even order). No constraint iteration is needed: the rigid geometry is exact to round-off.
//!
//! Site-site forces come from the usual [`ForceField`] (intramolecular pairs are excluded
//! automatically), so Lennard-Jones, Coulomb (also with Ewald real-space screening) and
//! user-defined terms work unchanged. Linear molecules (diatomics, CO2) and single-site
//! particles are supported (zero principal moments are treated as non-rotating axes).
//! Pressure uses the molecular virial.


use crate::kernels::random::Rng;
use crate::sim::md_engine::{compute_forces, ForceField, Forces, MdSystem, NeighbourList};
use crate::sim::symplectic::{quat_to_matrix, RigidBody};

type V3 = [f64; 3];
type M3 = [[f64; 3]; 3];

/// One rigid molecule.
#[derive(Clone, Debug)]
pub struct RigidMolecule {
    /// Indices of the sites in the particle system.
    pub sites: Vec<usize>,
    /// Site coordinates in the principal body frame (about the centre of mass).
    pub body: Vec<V3>,
    /// Total mass.
    pub mass: f64,
    /// Centre-of-mass position (unwrapped).
    pub com: V3,
    /// Centre-of-mass velocity.
    pub vcom: V3,
    /// Orientation and angular momentum.
    pub rb: RigidBody,
}

/// Rigid-molecule MD driver.
pub struct RigidMdSimulation {
    /// Site system (positions and velocities are kept consistent with the rigid bodies).
    pub sys: MdSystem<3>,
    /// Force field acting between sites.
    pub ff: ForceField<3>,
    /// Neighbour list.
    pub nl: NeighbourList<3>,
    /// Current site forces.
    pub forces: Forces<3>,
    /// The rigid molecules.
    pub mols: Vec<RigidMolecule>,
    /// Time step.
    pub dt: f64,
    /// Even composition order of the rotational flow (2 = Strang).
    pub order: usize,
    /// Use exact symmetric-top flows (true) or `NO_SQUISH` axis rotations (false).
    pub top_split: bool,
    /// Use the closed-form Jacobi-elliptic free flow of the asymmetric top
    /// ([`RigidBody::free_flow_exact`]) instead of a splitting (overrides `top_split`).
    pub exact_rotor: bool,
    rng: Rng,
}

fn cross(
    a: V3,
    b: V3,
) -> V3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

fn sym_eigen(a: M3) -> ([f64; 3], M3) {
    let mut m = a;
    let mut v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for _ in 0..60 {
        let off = m[0][1].abs() + m[0][2].abs() + m[1][2].abs();
        if off < 1e-15 {
            break;
        }
        for (p, q) in [(0, 1), (0, 2), (1, 2)] {
            if m[p][q].abs() < 1e-300 {
                continue;
            }
            let theta = (m[q][q] - m[p][p]) / (2.0 * m[p][q]);
            let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
            let t = if theta == 0.0 { 1.0 } else { t };
            let c = 1.0 / (t * t + 1.0).sqrt();
            let s = t * c;
            for k in 0..3 {
                let (mkp, mkq) = (m[k][p], m[k][q]);
                m[k][p] = c * mkp - s * mkq;
                m[k][q] = s * mkp + c * mkq;
            }
            for k in 0..3 {
                let (mpk, mqk) = (m[p][k], m[q][k]);
                m[p][k] = c * mpk - s * mqk;
                m[q][k] = s * mpk + c * mqk;
            }
            for k in 0..3 {
                let (vkp, vkq) = (v[k][p], v[k][q]);
                v[k][p] = c * vkp - s * vkq;
                v[k][q] = s * vkp + c * vkq;
            }
        }
    }
    ([m[0][0], m[1][1], m[2][2]], v)
}

fn mat_to_quat(r: &M3) -> [f64; 4] {
    let tr = r[0][0] + r[1][1] + r[2][2];
    let q = if tr > 0.0 {
        let s = 2.0 * (tr + 1.0).sqrt();
        [0.25 * s, (r[2][1] - r[1][2]) / s, (r[0][2] - r[2][0]) / s, (r[1][0] - r[0][1]) / s]
    } else if r[0][0] > r[1][1] && r[0][0] > r[2][2] {
        let s = 2.0 * (1.0 + r[0][0] - r[1][1] - r[2][2]).sqrt();
        [(r[2][1] - r[1][2]) / s, 0.25 * s, (r[0][1] + r[1][0]) / s, (r[0][2] + r[2][0]) / s]
    } else if r[1][1] > r[2][2] {
        let s = 2.0 * (1.0 + r[1][1] - r[0][0] - r[2][2]).sqrt();
        [(r[0][2] - r[2][0]) / s, (r[0][1] + r[1][0]) / s, 0.25 * s, (r[1][2] + r[2][1]) / s]
    } else {
        let s = 2.0 * (1.0 + r[2][2] - r[0][0] - r[1][1]).sqrt();
        [(r[1][0] - r[0][1]) / s, (r[0][2] + r[2][0]) / s, (r[1][2] + r[2][1]) / s, 0.25 * s]
    };
    let n = q.iter().map(|x| x * x).sum::<f64>().sqrt();
    [q[0] / n, q[1] / n, q[2] / n, q[3] / n]
}

impl RigidMdSimulation {
    /// Builds the rigid-body state from the site positions and velocities in `sys`.
    ///
    /// `molecules[m]` lists the site indices of molecule `m`; the lists must partition
    /// the sites. The current site positions must already form rigid molecules (their
    /// shape defines the body frame); velocities are projected onto rigid motion.
    ///
    /// # Errors
    /// Fails if the molecules do not partition the sites or the box is not fully periodic.
    pub fn new(
        sys: MdSystem<3>,
        mut ff: ForceField<3>,
        skin: f64,
        dt: f64,
        molecules: &[Vec<usize>],
    ) -> Result<Self, String> {
        let n = sys.n();
        let mut seen = vec![false; n];
        for m in molecules {
            if m.is_empty() {
                return Err("empty molecule".into());
            }
            for &i in m {
                if i >= n || seen[i] {
                    return Err("molecules must partition the sites".into());
                }
                seen[i] = true;
            }
        }
        if seen.iter().any(|s| !s) {
            return Err("every site must belong to a molecule".into());
        }
        if sys.bx.bc.iter().any(|b| *b != crate::sim::md_engine::Boundary::Periodic) || sys.bx.shear_rate != 0.0 {
            return Err("rigid-molecule MD needs a fully periodic box without shear".into());
        }
        let mut mols = Vec::new();
        for sites in molecules {
            // unwrap the molecule around its first site
            let r0 = sys.pos[sites[0]];
            let mut x: Vec<V3> = Vec::new();
            for &i in sites {
                let mut d = [sys.pos[i][0] - r0[0], sys.pos[i][1] - r0[1], sys.pos[i][2] - r0[2]];
                sys.bx.min_image(&mut d);
                x.push([r0[0] + d[0], r0[1] + d[1], r0[2] + d[2]]);
            }
            let mass: f64 = sites.iter().map(|&i| sys.mass[i]).sum();
            let mut com = [0.0; 3];
            let mut vcom = [0.0; 3];
            for (k, &i) in sites.iter().enumerate() {
                for a in 0..3 {
                    com[a] += sys.mass[i] * x[k][a] / mass;
                    vcom[a] += sys.mass[i] * sys.vel[i][a] / mass;
                }
            }
            let rel: Vec<V3> = x.iter().map(|p| [p[0] - com[0], p[1] - com[1], p[2] - com[2]]).collect();
            let mut inertia = [[0.0; 3]; 3];
            let mut l = [0.0; 3];
            for (k, &i) in sites.iter().enumerate() {
                let m = sys.mass[i];
                let r = rel[k];
                let r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
                for a in 0..3 {
                    for b in 0..3 {
                        inertia[a][b] += m * (if a == b { r2 } else { 0.0 } - r[a] * r[b]);
                    }
                }
                let v = [sys.vel[i][0] - vcom[0], sys.vel[i][1] - vcom[1], sys.vel[i][2] - vcom[2]];
                let c = cross(r, v);
                for a in 0..3 {
                    l[a] += m * c[a];
                }
            }
            let (mut ev, mut vec) = sym_eigen(inertia);
            // sort ascending for reproducibility and make the frame right-handed
            let mut idx = [0usize, 1, 2];
            idx.sort_by(|&a, &b| ev[a].total_cmp(&ev[b]));
            let cols: Vec<V3> = idx.iter().map(|&c| [vec[0][c], vec[1][c], vec[2][c]]).collect();
            ev = [ev[idx[0]], ev[idx[1]], ev[idx[2]]];
            let c2 = cross(cols[0], cols[1]);
            let sign = if c2[0] * cols[2][0] + c2[1] * cols[2][1] + c2[2] * cols[2][2] < 0.0 { -1.0 } else { 1.0 };
            vec = [
                [cols[0][0], cols[1][0], sign * cols[2][0]],
                [cols[0][1], cols[1][1], sign * cols[2][1]],
                [cols[0][2], cols[1][2], sign * cols[2][2]],
            ];
            let imax = ev[2].max(0.0);
            let mut inertia_p = ev;
            for v in &mut inertia_p {
                if *v <= 1e-10 * imax || imax == 0.0 {
                    *v = f64::INFINITY;
                }
            }
            let body: Vec<V3> = rel
            .iter()
            .map(|r| {
                let mut b = [0.0; 3];
                for k in 0..3 {
                    b[k] = vec[0][k] * r[0] + vec[1][k] * r[1] + vec[2][k] * r[2];
                }
                b
            })
            .collect();
        let mut l_body = [0.0; 3];
            for k in 0..3 {
                l_body[k] = vec[0][k] * l[0] + vec[1][k] * l[1] + vec[2][k] * l[2];
                if inertia_p[k].is_infinite() {
                    l_body[k] = 0.0;
                }
            }
            let rb = RigidBody { inertia: inertia_p, q: mat_to_quat(&vec), l_body };
            mols.push(RigidMolecule { sites: sites.clone(), body, mass, com, vcom, rb });
        }
        for sites in molecules {
            for (a, &i) in sites.iter().enumerate() {
                for &j in &sites[a + 1..] {
                    ff.extra_exclusions.push((i, j));
                }
            }
        }
        ff.rebuild_exclusions();
        let mut nl = NeighbourList::new(ff.rcut, skin);
        let forces = compute_forces(&sys, &ff, &mut nl)?;
        let mut s = Self { sys, ff, nl, forces, mols, dt, order: 2, top_split: true, exact_rotor: false, rng: Rng::new(1) };
        s.sync_sites();
        s.forces = compute_forces(&s.sys, &s.ff, &mut s.nl)?;
        Ok(s)
    }

    /// Writes site positions and velocities from the rigid-body state.
    pub fn sync_sites(&mut self) {
        for m in &self.mols {
            let r = quat_to_matrix(&m.rb.q);
            let wb = [m.rb.l_body[0] / m.rb.inertia[0], m.rb.l_body[1] / m.rb.inertia[1], m.rb.l_body[2] / m.rb.inertia[2]];
            let mut ws = [0.0; 3];
            for i in 0..3 {
                for k in 0..3 {
                    ws[i] += r[i][k] * wb[k];
                }
            }
            for (k, &site) in m.sites.iter().enumerate() {
                let mut d = [0.0; 3];
                for i in 0..3 {
                    for j in 0..3 {
                        d[i] += r[i][j] * m.body[k][j];
                    }
                }
                let w = cross(ws, d);
                for a in 0..3 {
                    self.sys.pos[site][a] = m.com[a] + d[a];
                    self.sys.unwrapped[site][a] = m.com[a] + d[a];
                    self.sys.vel[site][a] = m.vcom[a] + w[a];
                }
            }
        }
        self.sys.apply_boundaries(&mut self.rng);
    }

    /// Lever arms `R b_s` of the sites of molecule `m`.
    fn arms(
        &self,
        m: usize,
    ) -> Vec<V3> {
        let mol = &self.mols[m];
        let r = quat_to_matrix(&mol.rb.q);
        mol.body
            .iter()
            .map(|b| {
                let mut d = [0.0; 3];
                for i in 0..3 {
                    for j in 0..3 {
                        d[i] += r[i][j] * b[j];
                    }
                }
                d
            })
            .collect()
    }

    fn half_kick(
        &mut self,
        h: f64,
    ) {
        for m in 0..self.mols.len() {
            let arms = self.arms(m);
            let mut f = [0.0; 3];
            let mut tau = [0.0; 3];
            for (k, &site) in self.mols[m].sites.iter().enumerate() {
                let fs = self.forces.f[site];
                let t = cross(arms[k], fs);
                for a in 0..3 {
                    f[a] += fs[a];
                    tau[a] += t[a];
                }
            }
            let mass = self.mols[m].mass;
            for a in 0..3 {
                self.mols[m].vcom[a] += 0.5 * h * f[a] / mass;
            }
            self.mols[m].rb.kick_space(0.5 * h, &tau);
            // keep non-rotating axes at rest
            for k in 0..3 {
                if self.mols[m].rb.inertia[k].is_infinite() {
                    self.mols[m].rb.l_body[k] = 0.0;
                }
            }
        }
    }

    /// Advances one step (velocity-Verlet for the centre of mass, symplectic rotor splitting).
    ///
    /// # Errors
    /// Propagates force-field errors or an invalid composition order.
    pub fn step(&mut self) -> Result<(), String> {
        let h = self.dt;
        self.half_kick(h);
        for m in &mut self.mols {
            for a in 0..3 {
                m.com[a] += h * m.vcom[a];
            }
            if self.exact_rotor {
                for c in crate::sim::symplectic::composition_weights(self.order)? {
                    m.rb.free_flow_exact(c * h);
                }
            } else {
                m.rb.free_step_composed(h, self.top_split, self.order)?;
            }
        }
        self.sync_sites();
        self.forces = compute_forces(&self.sys, &self.ff, &mut self.nl)?;
        self.half_kick(h);
        self.sync_sites();
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

    /// Translational plus rotational kinetic energy.
    #[must_use]
    pub fn kinetic_energy(&self) -> f64 {
        self.mols.iter().map(|m| 0.5 * m.mass * m.vcom.iter().map(|v| v * v).sum::<f64>() + m.rb.kinetic_energy()).sum()
    }

    /// Total energy.
    #[must_use]
    pub fn total_energy(&self) -> f64 {
        self.kinetic_energy() + self.forces.epot
    }

    /// Degrees of freedom: 3 translations plus the non-degenerate rotations, minus the 3 of the total momentum.
    #[must_use]
    pub fn dof(&self) -> f64 {
        let d: usize = self.mols.iter().map(|m| 3 + m.rb.inertia.iter().filter(|i| i.is_finite()).count()).sum();
        d.saturating_sub(3) as f64
    }

    /// Instantaneous temperature.
    #[must_use]
    pub fn temperature(&self) -> f64 {
        2.0 * self.kinetic_energy() / self.dof()
    }

    /// Total linear momentum.
    #[must_use]
    pub fn momentum(&self) -> V3 {
        let mut p = [0.0; 3];
        for m in &self.mols {
            for a in 0..3 {
                p[a] += m.mass * m.vcom[a];
            }
        }
        p
    }

    /// Molecular pressure tensor `(sum M v v + W_mol) / V` with the virial of the
    /// centre-of-mass separations (site virial minus the intramolecular lever arms).
    #[must_use]
    pub fn pressure_tensor(&self) -> M3 {
        let mut w = self.forces.virial;
        let mut kin = [[0.0; 3]; 3];
        for (mi, m) in self.mols.iter().enumerate() {
            let arms = self.arms(mi);
            for (k, &site) in m.sites.iter().enumerate() {
                for a in 0..3 {
                    for b in 0..3 {
                        w[a][b] -= arms[k][a] * self.forces.f[site][b];
                    }
                }
            }
            for a in 0..3 {
                for b in 0..3 {
                    kin[a][b] += m.mass * m.vcom[a] * m.vcom[b];
                }
            }
        }
        let v = self.sys.bx.volume();
        let mut p = [[0.0; 3]; 3];
        for a in 0..3 {
            for b in 0..3 {
                p[a][b] = (kin[a][b] + w[a][b]) / v;
            }
        }
        p
    }

    /// Maximum deviation of any intramolecular site distance from its rigid value (diagnostic).
    #[must_use]
    pub fn max_geometry_error(&self) -> f64 {
        let mut worst = 0.0_f64;
        for (mi, m) in self.mols.iter().enumerate() {
            let arms = self.arms(mi);
            for a in 0..m.sites.len() {
                for b in a + 1..m.sites.len() {
                    let d0: f64 = (0..3).map(|k| (m.body[a][k] - m.body[b][k]).powi(2)).sum::<f64>().sqrt();
                    let d1: f64 = (0..3).map(|k| (arms[a][k] - arms[b][k]).powi(2)).sum::<f64>().sqrt();
                    worst = worst.max((d0 - d1).abs());
                }
            }
        }
        worst
    }
}
