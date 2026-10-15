//! Flexible-cell (Parrinello-Rahman / Martyna-Tobias-Klein) barostat for `MdSimulation`.
//!
//! The cell matrix `h` (rows are lattice vectors) is a dynamical variable with
//! momentum `p_g` and mass `W`: with `v_g = p_g / W` the equations of motion are
//!
//! ```text
//! dr_i/dt = p_i/m_i + v_g r_i              dh/dt = v_g h
//! dp_i/dt = F_i - (v_g + Tr(v_g)/N_f) p_i  dp_g/dt = V (P_int - P_ext) + (2K/N_f) I
//! ```
//!
//! (`P_int V = sum p p^T / m + W_virial`). Particle and cell thermostats are
//! Nose-Hoover chains (MTTK), and holonomic distance constraints are carried
//! through SHAKE/RATTLE with the cell flow included in the velocity condition.
//! The Liouville splitting is symmetric, so the scheme is reversible and second order.

use super::{Barostat, MdSimulation, NoseHooverChain, Thermostat};

/// Square matrix helper type.
type Mat<const D: usize> = [[f64; D]; D];

/// Matrix product `a b`.
#[must_use]
pub fn mat_mul<const D: usize>(
    a: &Mat<D>,
    b: &Mat<D>,
) -> Mat<D> {
    let mut c = [[0.0; D]; D];
    for i in 0..D {
        for k in 0..D {
            for j in 0..D {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}

/// Matrix-vector product `a x`.
#[must_use]
pub fn mat_vec<const D: usize>(
    a: &Mat<D>,
    x: &[f64; D],
) -> [f64; D] {
    let mut y = [0.0; D];
    for i in 0..D {
        for j in 0..D {
            y[i] += a[i][j] * x[j];
        }
    }
    y
}

/// Determinant by Gaussian elimination with partial pivoting.
#[must_use]
pub fn mat_det<const D: usize>(a: &Mat<D>) -> f64 {
    let mut m = *a;
    let mut det = 1.0;
    for c in 0..D {
        let mut piv = c;
        for r in c + 1..D {
            if m[r][c].abs() > m[piv][c].abs() {
                piv = r;
            }
        }
        if m[piv][c] == 0.0 {
            return 0.0;
        }
        if piv != c {
            m.swap(piv, c);
            det = -det;
        }
        det *= m[c][c];
        for r in c + 1..D {
            let f = m[r][c] / m[c][c];
            for k in c..D {
                m[r][k] -= f * m[c][k];
            }
        }
    }
    det
}

/// Inverse by Gauss-Jordan elimination (`None` if singular).
#[must_use]
pub fn mat_inverse<const D: usize>(a: &Mat<D>) -> Option<Mat<D>> {
    let mut m = *a;
    let mut inv = [[0.0; D]; D];
    for i in 0..D {
        inv[i][i] = 1.0;
    }
    for c in 0..D {
        let mut piv = c;
        for r in c + 1..D {
            if m[r][c].abs() > m[piv][c].abs() {
                piv = r;
            }
        }
        if m[piv][c].abs() < 1e-300 {
            return None;
        }
        m.swap(piv, c);
        inv.swap(piv, c);
        let d = m[c][c];
        for k in 0..D {
            m[c][k] /= d;
            inv[c][k] /= d;
        }
        for r in 0..D {
            if r != c {
                let f = m[r][c];
                if f != 0.0 {
                    for k in 0..D {
                        m[r][k] -= f * m[c][k];
                        inv[r][k] -= f * inv[c][k];
                    }
                }
            }
        }
    }
    Some(inv)
}

/// Matrix exponential `e^A` and the first phi-function `phi1(A) = sum A^k/(k+1)!`
/// (so that `x(t) = e^{At} x0 + t phi1(At) b` solves `x' = A x + b`).
#[must_use]
pub fn mat_exp_phi1<const D: usize>(a: &Mat<D>) -> (Mat<D>, Mat<D>) {
    let norm = a.iter().map(|r| r.iter().map(|x| x.abs()).sum::<f64>()).fold(0.0_f64, f64::max);
    let mut s = 0u32;
    let mut scale = 1.0;
    while norm * scale > 0.25 && s < 60 {
        scale *= 0.5;
        s += 1;
    }
    let mut x = *a;
    for r in &mut x {
        for c in r.iter_mut() {
            *c *= scale;
        }
    }
    let mut id = [[0.0; D]; D];
    for i in 0..D {
        id[i][i] = 1.0;
    }
    let (mut e, mut phi) = (id, id);
    let mut term = id; // x^k / k!
    for k in 1..=16 {
        term = mat_mul(&term, &x);
        for r in &mut term {
            for c in r.iter_mut() {
                *c /= f64::from(k);
            }
        }
        for i in 0..D {
            for j in 0..D {
                e[i][j] += term[i][j];
                phi[i][j] += term[i][j] / f64::from(k + 1);
            }
        }
    }
    for _ in 0..s {
        let mut ep = e;
        for i in 0..D {
            ep[i][i] += 1.0;
        }
        let mut np = mat_mul(&ep, &phi);
        for r in &mut np {
            for c in r.iter_mut() {
                *c *= 0.5;
            }
        }
        phi = np;
        e = mat_mul(&e, &e);
    }
    (e, phi)
}

/// Which cell degrees of freedom fluctuate.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CellMode {
    /// All `D(D+1)/2` components of the (symmetrised) strain, including shears: the cell becomes triclinic.
    Full,
    /// Only the diagonal: the box stays orthorhombic, each edge fluctuates independently.
    Orthorhombic,
    /// A single isotropic degree of freedom (volume).
    Isotropic,
}

/// Flexible-cell barostat state (Parrinello-Rahman dynamics in the Martyna-Tobias-Klein form).
#[derive(Clone, Debug)]
pub struct FlexibleBarostat {
    /// Spatial dimension the barostat was built for.
    pub d: usize,
    /// Target (external) pressure tensor, row-major `d x d`; `p0 I` for hydrostatic pressure.
    pub p_ext: Vec<f64>,
    /// Piston (cell) mass `W`.
    pub w: f64,
    /// Cell momentum `p_g`, row-major `d x d`.
    pub pg: Vec<f64>,
    /// Active cell degrees of freedom.
    pub mode: CellMode,
    /// Optional Nose-Hoover chain thermostatting the cell momentum (MTTK).
    pub chain: Option<NoseHooverChain>,
}

impl FlexibleBarostat {
    /// Hydrostatic barostat at pressure `p0` with piston mass `w`.
    #[must_use]
    pub fn hydrostatic(
        d: usize,
        p0: f64,
        w: f64,
        mode: CellMode,
    ) -> Self {
        let mut p = vec![0.0; d * d];
        for i in 0..d {
            p[i * d + i] = p0;
        }
        Self { d, p_ext: p, w, pg: vec![0.0; d * d], mode, chain: None }
    }

    /// Same with a general target pressure tensor (row-major).
    #[must_use]
    pub fn with_tensor(
        d: usize,
        p_ext: Vec<f64>,
        w: f64,
        mode: CellMode,
    ) -> Self {
        Self { d, p_ext, w, pg: vec![0.0; d * d], mode, chain: None }
    }

    /// Number of independent cell degrees of freedom.
    #[must_use]
    pub const fn cell_dof(&self) -> usize {
        match self.mode {
            CellMode::Full => self.d * (self.d + 1) / 2,
            CellMode::Orthorhombic => self.d,
            CellMode::Isotropic => 1,
        }
    }

    /// Couples a Nose-Hoover chain of `m` thermostats (relaxation time `tau`, temperature `t0`) to the cell.
    #[must_use]
    pub fn with_chain(
        mut self,
        t0: f64,
        tau: f64,
        m: usize,
    ) -> Self {
        let nf = self.cell_dof() as f64;
        self.chain = Some(NoseHooverChain::new(t0, tau, m, nf));
        self
    }

    /// Cell contribution to the conserved quantity at volume `v`: `sum p_g^2/(2W) + P V (+ chain)`.
    #[must_use]
    pub fn energy(
        &self,
        v: f64,
    ) -> f64 {
        let d = self.d;
        let ke: f64 = self.pg.iter().map(|x| x * x).sum::<f64>() / (2.0 * self.w);
        let p0 = (0..d).map(|i| self.p_ext[i * d + i]).sum::<f64>() / d as f64;
        ke + p0 * v + self.chain.as_ref().map_or(0.0, NoseHooverChain::energy)
    }

    fn velocity_gradient<const D: usize>(&self) -> Mat<D> {
        let mut g = [[0.0; D]; D];
        for i in 0..D {
            for j in 0..D {
                g[i][j] = self.pg[i * D + j] / self.w;
            }
        }
        g
    }
}

impl<const D: usize> MdSimulation<D> {
    const fn constrained_dof(&self) -> f64 {
        (D * self.sys.n()).saturating_sub(self.constraints.len()) as f64
    }

    /// Full virial tensor `sum m v v^T + virial + constraint virial` (equal to `V P_int`).
    fn internal_stress_volume(&self) -> Mat<D> {
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
                p[a][b] += self.forces.virial[a][b] + self.constraint_virial[a][b];
            }
        }
        p
    }

    fn cell_kick(
        &self,
        b: &mut FlexibleBarostat,
        dt: f64,
    ) {
        let v = self.sys.bx.volume();
        let nf = self.constrained_dof();
        let w = self.internal_stress_volume();
        let ke2 = 2.0 * self.sys.kinetic_energy();
        let mut g = [[0.0; D]; D];
        for a in 0..D {
            for c in 0..D {
                g[a][c] = w[a][c] - v * b.p_ext[a * D + c];
            }
            g[a][a] += ke2 / nf;
        }
        match b.mode {
            CellMode::Full => {
                for a in 0..D {
                    for c in a + 1..D {
                        let m = f64::midpoint(g[a][c], g[c][a]);
                        g[a][c] = m;
                        g[c][a] = m;
                    }
                }
            },
            CellMode::Orthorhombic => {
                for a in 0..D {
                    for c in 0..D {
                        if a != c {
                            g[a][c] = 0.0;
                        }
                    }
                }
            },
            CellMode::Isotropic => {
                let tr = (0..D).map(|a| g[a][a]).sum::<f64>() / D as f64;
                for a in 0..D {
                    for c in 0..D {
                        g[a][c] = if a == c { tr } else { 0.0 };
                    }
                }
            },
        }
        for a in 0..D {
            for c in 0..D {
                b.pg[a * D + c] += dt * g[a][c];
            }
        }
    }

    fn thermostat_half(
        &mut self,
        b: &mut FlexibleBarostat,
    ) {
        let h = self.dt;
        if let Thermostat::NoseHooverChain(_) = &self.thermostat {
            let ke2 = 2.0 * self.sys.kinetic_energy();
            let s = if let Thermostat::NoseHooverChain(c) = &mut self.thermostat { c.propagate(ke2, h) } else { 1.0 };
            self.scale_velocities(s);
        }
        if let Some(c) = &mut b.chain {
            let ke2 = b.pg.iter().map(|x| x * x).sum::<f64>() / b.w;
            let s = c.propagate(ke2, h);
            for x in &mut b.pg {
                *x *= s;
            }
        }
    }

    pub(super) fn step_flexible(&mut self) -> Result<(), String> {
        let Barostat::Flexible(mut b) = std::mem::replace(&mut self.barostat, Barostat::None) else { return Ok(()) };
        let r = self.flexible_inner(&mut b);
        self.barostat = Barostat::Flexible(b);
        r
    }

    fn flexible_inner(
        &mut self,
        b: &mut FlexibleBarostat,
    ) -> Result<(), String> {
        if b.d != D {
            return Err("flexible barostat built for a different dimension".into());
        }
        if matches!(self.thermostat, Thermostat::Langevin { .. }) {
            return Err("the flexible-cell barostat composes with Nose-Hoover chains, Berendsen or no thermostat, not Langevin".into());
        }
        if self.sys.bx.shear_rate != 0.0 || self.sys.bx.bc.iter().any(|c| *c != super::Boundary::Periodic) {
            return Err("the flexible-cell barostat needs a fully periodic box without Lees-Edwards shear".into());
        }
        if b.mode == CellMode::Full && self.sys.bx.cell.is_none() {
            let h0 = self.sys.bx.matrix();
            self.sys.bx.set_cell(h0)?;
        }
        let h = self.dt;
        let nf = self.constrained_dof();
        // 1. thermostats, cell kick
        self.thermostat_half(b);
        self.cell_kick(b, 0.5 * h);
        // 2. particle half kick with cell damping
        let damp = |s: &mut Self, b: &FlexibleBarostat, t: f64| {
            let vg: Mat<D> = b.velocity_gradient();
            let alpha = (0..D).map(|a| vg[a][a]).sum::<f64>() / nf;
            let mut a_mat = vg;
            for a in 0..D {
                a_mat[a][a] += alpha;
            }
            for r in &mut a_mat {
                for c in r.iter_mut() {
                    *c *= -t;
                }
            }
            let (e, _) = mat_exp_phi1(&a_mat);
            for v in &mut s.sys.vel {
                *v = mat_vec(&e, v);
            }
        };
        damp(self, b, 0.25 * h);
        self.half_kick(h);
        damp(self, b, 0.25 * h);
        // 3. drift with the cell flow
        let vg: Mat<D> = b.velocity_gradient();
        let mut hv = vg;
        for r in &mut hv {
            for c in r.iter_mut() {
                *c *= h;
            }
        }
        let (e, phi) = mat_exp_phi1(&hv);
        let old = self.sys.pos.clone();
        for i in 0..self.sys.n() {
            let drift = mat_vec(&phi, &self.sys.vel[i]);
            let p = mat_vec(&e, &self.sys.pos[i]);
            let u = mat_vec(&e, &self.sys.unwrapped[i]);
            for a in 0..D {
                self.sys.pos[i][a] = p[a] + h * drift[a];
                self.sys.unwrapped[i][a] = u[a] + h * drift[a];
            }
        }
        if let Some(cell) = self.sys.bx.cell {
            // lattice vectors are columns of the deformation: a_k <- E a_k
            let mut nc = cell;
            for k in 0..D {
                nc[k] = mat_vec(&e, &cell[k]);
            }
            self.sys.bx.set_cell(nc)?;
        } else {
            for a in 0..D {
                self.sys.bx.len[a] *= e[a][a];
            }
        }
        if !self.constraints.is_empty() {
            self.shake(&old, h);
        }
        self.sys.apply_boundaries(&mut self.rng);
        // 4. forces, second half kick
        self.recompute()?;
        damp(self, b, 0.25 * h);
        self.half_kick(h);
        damp(self, b, 0.25 * h);
        if !self.constraints.is_empty() {
            let vg: Mat<D> = b.velocity_gradient();
            self.rattle_velocities(Some(&vg));
        }
        self.cell_kick(b, 0.5 * h);
        self.thermostat_half(b);
        if let Thermostat::Berendsen { t0, tau } = self.thermostat {
            let t = self.temperature();
            if t > 0.0 {
                let lam = (1.0 + h / tau * (t0 / t - 1.0)).max(0.0).sqrt();
                self.scale_velocities(lam);
            }
        }
        self.sys.time += h;
        Ok(())
    }
}
