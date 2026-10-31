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
//!
//! Extensions: a **Langevin piston** (Ornstein-Uhlenbeck step of the cell momentum with the
//! particle temperature, exact equipartition of the cell modes, checked against the ideal-gas
//! NPT distribution `p(V) ~ V^N exp(-PV/T)`), **SLLOD shear** on a deforming cell
//! ([`FlexibleBarostat::shear_flow`]: imposed velocity gradient, particle equations
//! `dr/dt = v + G r`, `dv/dt = F/m - G v` integrated exactly for the flow, lattice-equivalent
//! reduction of the sheared cell), **sliding-brick Lees-Edwards** boundaries combined with the
//! orthorhombic / isotropic cell (peculiar-velocity dynamics with the total velocity gradient
//! `G_cell + gamma e_x e_y^T`, offset advanced and scaled with the cell), and the constraint virial from the exact acceleration-level
//! constraint forces (including the flexible-cell flow terms) instead of the SHAKE impulse.

use super::Barostat;
use super::MdSimulation;
use super::NoseHooverChain;
use super::Thermostat;

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
    let norm = a
        .iter()
        .map(|r| r.iter().map(|x| x.abs()).sum::<f64>())
        .fold(0.0_f64, f64::max);
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
    /// Friction of the **Langevin piston**: when positive and the particle thermostat is
    /// [`Thermostat::Langevin`], the cell momentum receives an Ornstein-Uhlenbeck step of the same
    /// temperature (stochastic cell thermostat, exact equipartition of the cell degrees of freedom).
    pub cell_friction: f64,
    /// Prescribed constant velocity gradient `dr_a/dt = G_ab r_b` (row-major `d x d`). The cell
    /// then follows this flow instead of its own dynamics: with `G = gamma e_x e_y^T` this is
    /// SLLOD shear flow on a deforming (Lees-Edwards) cell, see [`FlexibleBarostat::shear_flow`].
    pub prescribed: Option<Vec<f64>>,
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
        Self {
            d,
            p_ext: p,
            w,
            pg: vec![0.0; d * d],
            mode,
            chain: None,
            cell_friction: 0.0,
            prescribed: None,
        }
    }

    /// Same with a general target pressure tensor (row-major).
    #[must_use]
    pub fn with_tensor(
        d: usize,
        p_ext: Vec<f64>,
        w: f64,
        mode: CellMode,
    ) -> Self {
        Self {
            d,
            p_ext,
            w,
            pg: vec![0.0; d * d],
            mode,
            chain: None,
            cell_friction: 0.0,
            prescribed: None,
        }
    }

    /// SLLOD shear flow on a deforming cell: the velocity gradient `u_x = rate * y` is imposed
    /// (no cell dynamics), particles obey `dr/dt = v + G r`, `dv/dt = F/m - G v`, and the lattice
    /// vectors are advected by the flow and reduced (lattice-equivalent re-wrapping) when the
    /// shear exceeds half a box. Combine with a Nose-Hoover chain or Langevin thermostat acting on
    /// the peculiar velocities. Requires `d >= 2`.
    #[must_use]
    pub fn shear_flow(
        d: usize,
        rate: f64,
    ) -> Self {
        let mut g = vec![0.0; d * d];
        if d >= 2 {
            g[1] = rate;
        }
        let mut b = Self::hydrostatic(d, 0.0, 1.0, CellMode::Full);
        b.prescribed = Some(g);
        b
    }

    /// Adds a Langevin piston with the given friction (see [`FlexibleBarostat::cell_friction`]).
    #[must_use]
    pub const fn with_langevin_cell(
        mut self,
        gamma: f64,
    ) -> Self {
        self.cell_friction = gamma;
        self
    }

    /// Number of independent cell degrees of freedom.
    #[must_use]
    pub const fn cell_dof(&self) -> usize {
        match self.mode {
            | CellMode::Full => self.d * (self.d + 1) / 2,
            | CellMode::Orthorhombic => self.d,
            | CellMode::Isotropic => 1,
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
                g[i][j] = self
                    .prescribed
                    .as_ref()
                    .map_or_else(|| self.pg[i * D + j] / self.w, |p| p[i * D + j]);
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

    /// Generalised force on the cell momentum, projected on the active cell modes.
    fn cell_force(
        &self,
        b: &FlexibleBarostat,
    ) -> Mat<D> {
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
            | CellMode::Full => {
                for a in 0..D {
                    for c in a + 1..D {
                        let m = f64::midpoint(g[a][c], g[c][a]);
                        g[a][c] = m;
                        g[c][a] = m;
                    }
                }
            },
            | CellMode::Orthorhombic => {
                for a in 0..D {
                    for c in 0..D {
                        if a != c {
                            g[a][c] = 0.0;
                        }
                    }
                }
            },
            | CellMode::Isotropic => {
                let tr = (0..D).map(|a| g[a][a]).sum::<f64>() / D as f64;
                for a in 0..D {
                    for c in 0..D {
                        g[a][c] = if a == c { tr } else { 0.0 };
                    }
                }
            },
        }
        g
    }

    fn cell_kick(
        &self,
        b: &mut FlexibleBarostat,
        dt: f64,
    ) {
        if b.prescribed.is_some() {
            return;
        }
        let g = self.cell_force(b);
        for a in 0..D {
            for c in 0..D {
                b.pg[a * D + c] += dt * g[a][c];
            }
        }
    }

    /// Instantaneous constraint-force virial `sum_c lambda_c r_c (x) r_c` from the **exact**
    /// acceleration-level constraint condition `d^2 |r_c|^2 / dt^2 = 0` (no lag in the time
    /// step): the multipliers solve `r_c . (a_i - a_j) + |dr_c/dt|^2 = 0` for the actual positions,
    /// velocities and non-constraint forces. With `flow = Some((vg, vg_dot, alpha))` the
    /// equations of motion of the flexible-cell dynamics are used (`dr/dt = v + vg r`,
    /// `dv/dt = F/m + G/m - (vg + alpha) v`, cell acceleration `vg_dot`).
    fn exact_constraint_virial(
        &self,
        flow: Option<(&Mat<D>, &Mat<D>, f64)>,
    ) -> Mat<D> {
        let nc = self.constraints.len();
        let zero = [[0.0; D]; D];
        let (vg, vgd, alpha) = flow.map_or((zero, zero, 0.0), |(a, b2, c)| (*a, *b2, c));
        let mut rs = Vec::with_capacity(nc);
        let mut base = Vec::with_capacity(nc);
        for c in &self.constraints {
            let (i, j) = (c.i, c.j);
            let mut r = [0.0; D];
            let mut u = [0.0; D];
            let mut df = [0.0; D];
            for a in 0..D {
                r[a] = self.sys.pos[i][a] - self.sys.pos[j][a];
                u[a] = self.sys.vel[i][a] - self.sys.vel[j][a];
                df[a] =
                    self.forces.f[i][a] / self.sys.mass[i] - self.forces.f[j][a] / self.sys.mass[j];
            }
            self.sys.bx.min_image(&mut r);
            let vg_r = mat_vec(&vg, &r);
            let rdot: [f64; D] = std::array::from_fn(|a| u[a] + vg_r[a]);
            let vg_u = mat_vec(&vg, &u);
            let vgd_r = mat_vec(&vgd, &r);
            let vg_rdot = mat_vec(&vg, &rdot);
            let mut acc = 0.0;
            for a in 0..D {
                acc += r[a] * (df[a] - vg_u[a] - alpha * u[a] + vgd_r[a] + vg_rdot[a])
                    + rdot[a] * rdot[a];
            }
            rs.push(r);
            base.push(acc);
        }
        let mut g = vec![[0.0; D]; self.sys.n()];
        let mut lam = vec![0.0; nc];
        for _ in 0..2000 {
            let mut worst = 0.0_f64;
            for (k, c) in self.constraints.iter().enumerate() {
                let (i, j) = (c.i, c.j);
                let (mi, mj) = (self.sys.mass[i], self.sys.mass[j]);
                let r = rs[k];
                let mut s = base[k];
                for a in 0..D {
                    s += r[a] * (g[i][a] / mi - g[j][a] / mj);
                }
                let rr: f64 = r.iter().map(|x| x * x).sum();
                let dl = -s / (rr * (1.0 / mi + 1.0 / mj));
                lam[k] += dl;
                for a in 0..D {
                    g[i][a] += dl * r[a];
                    g[j][a] -= dl * r[a];
                }
                worst = worst.max(s.abs() / (1.0 + base[k].abs()));
            }
            if worst < 1e-13 {
                break;
            }
        }
        let mut w = zero;
        for (k, r) in rs.iter().enumerate() {
            for a in 0..D {
                for b in 0..D {
                    w[a][b] += lam[k] * r[a] * r[b];
                }
            }
        }
        w
    }

    /// Replaces [`MdSimulation::constraint_virial`] by the exact instantaneous constraint-force
    /// virial of the current state (plain Newtonian flow, see the module notes).
    pub fn refresh_constraint_virial(&mut self) {
        if self.constraints.is_empty() {
            self.constraint_virial = [[0.0; D]; D];
        } else {
            self.constraint_virial = self.exact_constraint_virial(None);
        }
    }

    /// Exact constraint virial for the flexible-cell flow, with the cell acceleration iterated
    /// to self-consistency (the cell force contains the constraint virial).
    fn refresh_constraint_virial_flow(
        &mut self,
        b: &FlexibleBarostat,
    ) {
        if self.constraints.is_empty() {
            self.constraint_virial = [[0.0; D]; D];
            return;
        }
        let nf = self.constrained_dof();
        let vg: Mat<D> = self.total_gradient(b);
        let alpha = (0..D).map(|a| vg[a][a]).sum::<f64>() / nf;
        let mut vgd = [[0.0; D]; D];
        for _ in 0..if b.prescribed.is_some() {
            1
        } else {
            4
        } {
            self.constraint_virial = self.exact_constraint_virial(Some((&vg, &vgd, alpha)));
            if b.prescribed.is_none() {
                let g = self.cell_force(b);
                for a in 0..D {
                    for c in 0..D {
                        vgd[a][c] = g[a][c] / b.w;
                    }
                }
            }
        }
    }

    /// Ornstein-Uhlenbeck step of duration `dt` for the cell momentum (Langevin piston).
    fn cell_ou(
        &mut self,
        b: &mut FlexibleBarostat,
        t0: f64,
        dt: f64,
    ) {
        let c1 = (-b.cell_friction * dt).exp();
        let c2 = (1.0 - c1 * c1).sqrt();
        let kt_w = (b.w * t0).sqrt();
        match b.mode {
            | CellMode::Isotropic => {
                let xi = self.rng.normal();
                for a in 0..D {
                    let i = a * D + a;
                    if a == 0 {
                        b.pg[i] = c1 * b.pg[i] + c2 * kt_w / (D as f64).sqrt() * xi;
                    } else {
                        b.pg[i] = b.pg[0];
                    }
                }
            },
            | CellMode::Orthorhombic => {
                for a in 0..D {
                    b.pg[a * D + a] = c1 * b.pg[a * D + a] + c2 * kt_w * self.rng.normal();
                }
            },
            | CellMode::Full => {
                for a in 0..D {
                    b.pg[a * D + a] = c1 * b.pg[a * D + a] + c2 * kt_w * self.rng.normal();
                    for c in a + 1..D {
                        let v = c1 * b.pg[a * D + c]
                            + c2 * kt_w * std::f64::consts::FRAC_1_SQRT_2 * self.rng.normal();
                        b.pg[a * D + c] = v;
                        b.pg[c * D + a] = v;
                    }
                }
            },
        }
    }

    /// Lattice-equivalent reduction of a sheared lower-triangular cell: each row is shifted by
    /// integer multiples of the rows above it so that the off-diagonal elements stay within half
    /// the corresponding box length. Returns whether the cell changed.
    fn reduce_cell(&mut self) -> Result<bool, String> {
        let Some(mut cell) = self.sys.bx.cell else {
            return Ok(false);
        };
        let mut changed = false;
        for row in 1..D {
            for col in (0..row).rev() {
                let k = (cell[row][col] / cell[col][col]).round();
                if k != 0.0 {
                    for c in 0..D {
                        cell[row][c] -= k * cell[col][c];
                    }
                    changed = true;
                }
            }
        }
        if changed {
            self.sys.bx.set_cell(cell)?;
            self.nl.builds = 0;
        }
        Ok(changed)
    }

    fn thermostat_half(
        &mut self,
        b: &mut FlexibleBarostat,
    ) {
        let h = self.dt;
        if let Thermostat::NoseHooverChain(_) = &self.thermostat {
            let ke2 = 2.0 * self.sys.kinetic_energy();
            let s = if let Thermostat::NoseHooverChain(c) = &mut self.thermostat {
                c.propagate(ke2, h)
            } else {
                1.0
            };
            self.scale_velocities(s);
        }
        if let Thermostat::Langevin { t0, gamma } = self.thermostat {
            let c1 = (-gamma * 0.5 * h).exp();
            let c2 = (1.0 - c1 * c1).sqrt();
            for i in 0..self.sys.n() {
                let sd = (t0 / self.sys.mass[i]).sqrt();
                for a in 0..D {
                    self.sys.vel[i][a] = c1 * self.sys.vel[i][a] + c2 * sd * self.rng.normal();
                }
            }
            if b.cell_friction > 0.0 && b.prescribed.is_none() {
                self.cell_ou(b, t0, 0.5 * h);
            }
            if !self.constraints.is_empty() {
                let vg: Mat<D> = self.total_gradient(b);
                self.rattle_velocities(Some(&vg));
            }
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
        let Barostat::Flexible(mut b) = std::mem::replace(&mut self.barostat, Barostat::None)
        else {
            return Ok(());
        };
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
        if self
            .sys
            .bx
            .bc
            .iter()
            .any(|c| *c != super::Boundary::Periodic)
        {
            return Err("the flexible-cell barostat needs a fully periodic box".into());
        }
        let shear = self.sys.bx.shear_rate;
        if shear != 0.0
            && (D < 2
                || b.mode == CellMode::Full
                || b.prescribed.is_some()
                || self.sys.bx.cell.is_some())
        {
            return Err("sliding-brick Lees-Edwards needs an orthorhombic box and the Orthorhombic or Isotropic cell mode (use FlexibleBarostat::shear_flow for a deforming triclinic cell)".into());
        }
        if b.mode == CellMode::Full && self.sys.bx.cell.is_none() {
            let h0 = self.sys.bx.matrix();
            self.sys.bx.set_cell(h0)?;
        }
        let h = self.dt;
        let nf = self.constrained_dof();
        if shear != 0.0 {
            self.shear_profile_shift(-1.0);
        }
        if !self.constraints.is_empty()
            && self.constraint_virial.iter().flatten().all(|x| *x == 0.0)
        {
            self.refresh_constraint_virial_flow(b);
        }
        // 1. thermostats, cell kick
        self.thermostat_half(b);
        self.cell_kick(b, 0.5 * h);
        // 2. particle half kick with cell damping
        let damp = |s: &mut Self, b: &FlexibleBarostat, t: f64| {
            let vg: Mat<D> = s.total_gradient(b);
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
        // 3. drift with the cell flow (plus the shear flow of a sliding-brick box)
        let vg: Mat<D> = self.total_gradient(b);
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
            if shear != 0.0 {
                // the velocity profile is gamma (y - L_y/2): remove the uniform part of the shear drift
                let dx = 0.5 * shear * self.sys.bx.len[1] * h;
                self.sys.pos[i][0] -= dx;
                self.sys.unwrapped[i][0] -= dx;
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
        if b.prescribed.is_some() {
            self.reduce_cell()?;
        }
        if shear != 0.0 {
            let lx = self.sys.bx.len[0];
            let off = self.sys.bx.offset * e[0][0] + shear * self.sys.bx.len[1] * h;
            self.sys.bx.offset = off - lx * (off / lx).floor();
        }
        if !self.constraints.is_empty() {
            self.shake(&old, h);
        }
        self.wrap_peculiar(shear);
        // 4. forces, second half kick
        self.recompute()?;
        damp(self, b, 0.25 * h);
        self.half_kick(h);
        damp(self, b, 0.25 * h);
        if !self.constraints.is_empty() {
            let vg: Mat<D> = self.total_gradient(b);
            self.rattle_velocities(Some(&vg));
        }
        self.refresh_constraint_virial_flow(b);
        self.cell_kick(b, 0.5 * h);
        self.thermostat_half(b);
        if let Thermostat::Berendsen { t0, tau } = self.thermostat {
            let t = self.temperature();
            if t > 0.0 {
                let lam = (1.0 + h / tau * (t0 / t - 1.0)).max(0.0).sqrt();
                self.scale_velocities(lam);
            }
        }
        if shear != 0.0 {
            self.shear_profile_shift(1.0);
        }
        self.sys.time += h;
        Ok(())
    }

    /// Velocity gradient of the cell flow plus the Lees-Edwards shear (`du_x/dy`) of a sliding-brick box.
    fn total_gradient(
        &self,
        b: &FlexibleBarostat,
    ) -> Mat<D> {
        let mut g: Mat<D> = b.velocity_gradient();
        if D >= 2 {
            g[0][1] += self.sys.bx.shear_rate;
        }
        g
    }

    /// Converts between laboratory velocities of the sliding-brick convention
    /// `v_x = u_x + gamma (y - L_y/2)` and peculiar velocities `u` (`sign = -1`: lab to peculiar).
    fn shear_profile_shift(
        &mut self,
        sign: f64,
    ) {
        let (gamma, ly) = (self.sys.bx.shear_rate, self.sys.bx.len[1]);
        for i in 0..self.sys.n() {
            self.sys.vel[i][0] += sign * gamma * (self.sys.pos[i][1] - 0.5 * ly);
        }
    }

    /// Applies the periodic wrapping; for a sheared box the velocity jumps of the sliding-brick
    /// boundary are not applied because the velocities are peculiar (continuous) in the flexible step.
    fn wrap_peculiar(
        &mut self,
        shear: f64,
    ) {
        if shear == 0.0 {
            self.sys.apply_boundaries(&mut self.rng);
            return;
        }
        let vx: Vec<f64> = self.sys.vel.iter().map(|v| v[0]).collect();
        self.sys.apply_boundaries(&mut self.rng);
        for (v, x) in self.sys.vel.iter_mut().zip(vx) {
            v[0] = x;
        }
    }
}
