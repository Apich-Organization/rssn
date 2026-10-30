//! # Split-step Fourier solver for Schrödinger and Gross–Pitaevskii equations
//!
//! Solves `i ψ_t = -(1/2m) Δψ + V(t, x) ψ + g |ψ|² ψ` (units `ħ = 1`) in one,
//! two or three dimensions on a periodic box with power-of-two points per
//! axis, by operator splitting with FFTs ([`crate::sim::spectral`]):
//!
//! * real-time evolution with Strang splitting (second order) or the
//!   Yoshida triple-jump composition (fourth order), exactly unitary for the
//!   linear and the nonlinear sub-steps (the norm is conserved to round-off);
//! * **imaginary-time propagation** with renormalisation for ground states
//!   ([`SplitStep::ground_state`]), for the linear Schrödinger equation
//!   (`g = 0`) and the GPE;
//! * energy, chemical potential and norm diagnostics.
//!
//! | Equation | 1D | 2D | 3D | Boundary | Time stepping |
//! |---|---|---|---|---|---|
//! | linear Schrödinger | yes | yes | yes | periodic | Strang, Yoshida4, imaginary time |
//! | GPE (cubic nonlinearity) | yes | yes | yes | periodic | same |
//!
//! A **complex absorbing potential** `-i W(x)` ([`SplitStep::evolve_absorbing`]) damps the
//! wave function near the box edges, so outgoing waves in scattering and
//! expansion problems are absorbed instead of wrapping around the periodic box.
//!
//! Limitations: periodic boundaries (use the absorbing potential for open
//! problems), no rotating frame, no dipolar interactions.

use crate::sim::spectral::SpectralGrid;
use num_complex::Complex;

/// Splitting scheme.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SplitOrder {
    /// Strang splitting (order 2).
    Strang,
    /// Yoshida composition of three Strang steps (order 4).
    Yoshida4,
}

/// Ground state found by imaginary-time propagation.
#[derive(Clone, Debug)]
pub struct GroundState {
    /// Normalised wavefunction.
    pub psi: Vec<Complex<f64>>,
    /// Total energy per particle.
    pub energy: f64,
    /// Chemical potential.
    pub mu: f64,
    /// Steps performed.
    pub steps: usize,
    /// Whether the energy change fell below the tolerance.
    pub converged: bool,
}

/// Split-step solver on a periodic box.
#[derive(Clone, Debug)]
pub struct SplitStep {
    /// Spectral grid.
    pub grid: SpectralGrid,
    /// Coordinate of the grid corner.
    pub origin: [f64; 3],
    /// Particle mass.
    pub mass: f64,
}

impl SplitStep {
    /// Box `[-len/2, len/2]^dim` with `n` points per axis.
    ///
    /// # Errors
    /// See [`SpectralGrid::new`].
    pub fn centered(
        n: &[usize],
        len: &[f64],
        mass: f64,
    ) -> Result<Self, String> {
        let grid = SpectralGrid::new(n, len)?;
        let mut origin = [0.0; 3];
        for a in 0..grid.dim {
            origin[a] = -0.5 * len[a];
        }
        Ok(Self { grid, origin, mass })
    }

    /// Physical coordinate of point `p`.
    #[must_use]
    pub fn x(
        &self,
        p: usize,
    ) -> [f64; 3] {
        let c = self.grid.coord(p);
        [
            self.origin[0] + c[0],
            self.origin[1] + c[1],
            self.origin[2] + c[2],
        ]
    }

    fn cell(&self) -> f64 {
        (0..self.grid.dim).map(|a| self.grid.h(a)).product()
    }

    /// `∫ |ψ|²`.
    #[must_use]
    pub fn norm2(
        &self,
        psi: &[Complex<f64>],
    ) -> f64 {
        psi.iter().map(Complex::norm_sqr).sum::<f64>() * self.cell()
    }

    /// Samples a complex function on the grid.
    pub fn sample(
        &self,
        f: impl Fn([f64; 3]) -> Complex<f64>,
    ) -> Vec<Complex<f64>> {
        (0..self.grid.size()).map(|p| f(self.x(p))).collect()
    }

    /// Kinetic energy `∫ |∇ψ|² / (2m)`.
    #[must_use]
    pub fn kinetic_energy(
        &self,
        psi: &[Complex<f64>],
    ) -> f64 {
        let ph = self.grid.to_spec_c(psi);
        let n = self.grid.size() as f64;
        // Parseval: sum |psi|^2 = sum |psi_hat|^2 / N
        ph.iter()
            .enumerate()
            .map(|(p, v)| 0.5 / self.mass * self.grid.k2(p) * v.norm_sqr())
            .sum::<f64>()
            / n
            * self.cell()
    }

    /// Total energy `∫ [|∇ψ|²/2m + V|ψ|² + g|ψ|⁴/2]`.
    pub fn energy(
        &self,
        psi: &[Complex<f64>],
        v: &dyn Fn(f64, [f64; 3]) -> f64,
        t: f64,
        g: f64,
    ) -> f64 {
        let pot: f64 = psi
            .iter()
            .enumerate()
            .map(|(p, c)| {
                let d = c.norm_sqr();
                (v(t, self.x(p)) + 0.5 * g * d) * d
            })
            .sum::<f64>()
            * self.cell();
        self.kinetic_energy(psi) + pot
    }

    /// `∫ ψ* (H ψ)` with the nonlinear term `g|ψ|²` counted fully (chemical potential).
    pub fn chemical_potential(
        &self,
        psi: &[Complex<f64>],
        v: &dyn Fn(f64, [f64; 3]) -> f64,
        t: f64,
        g: f64,
    ) -> f64 {
        let pot: f64 = psi
            .iter()
            .enumerate()
            .map(|(p, c)| {
                let d = c.norm_sqr();
                (v(t, self.x(p)) + g * d) * d
            })
            .sum::<f64>()
            * self.cell();
        self.kinetic_energy(psi) + pot
    }

    #[allow(clippy::too_many_arguments)] // state, potential, nonlinearity, step, mode, kinetic factors and the optional absorber
    fn strang(
        &self,
        psi: &mut [Complex<f64>],
        vv: &[f64],
        g: f64,
        dt: f64,
        imaginary: bool,
        kin: &[Complex<f64>],
        damp: Option<&[f64]>,
    ) {
        let half = |psi: &mut [Complex<f64>]| {
            for (p, c) in psi.iter_mut().enumerate() {
                let phase = f64::midpoint(vv[p], g * c.norm_sqr()) * dt;
                *c *= if imaginary {
                    Complex::new((-phase).exp(), 0.0)
                } else {
                    Complex::new(0.0, -phase).exp()
                };
                if let Some(w) = damp {
                    *c *= (-0.5 * w[p] * dt).exp();
                }
            }
        };
        half(psi);
        let mut ph = self.grid.to_spec_c(psi);
        for (a, k) in ph.iter_mut().zip(kin) {
            *a *= k;
        }
        let back = self.grid.to_phys_c(&ph);
        psi.copy_from_slice(&back);
        half(psi);
    }

    fn kinetic_factors(
        &self,
        dt: f64,
        imaginary: bool,
    ) -> Vec<Complex<f64>> {
        (0..self.grid.size())
            .map(|p| {
                let e = 0.5 / self.mass * self.grid.k2(p) * dt;
                if imaginary {
                    Complex::new((-e).exp(), 0.0)
                } else {
                    Complex::new(0.0, -e).exp()
                }
            })
            .collect()
    }

    /// Real-time evolution for `steps` steps of size `dt` from time `t0`.
    #[allow(clippy::too_many_arguments)] // public evolution signature (psi0, v, g, t0, dt, steps, order) kept for API compatibility
    pub fn evolve(
        &self,
        psi0: &[Complex<f64>],
        v: &dyn Fn(f64, [f64; 3]) -> f64,
        g: f64,
        t0: f64,
        dt: f64,
        steps: usize,
        order: SplitOrder,
    ) -> Vec<Complex<f64>> {
        let mut psi = psi0.to_vec();
        let c3 = 2.0_f64.cbrt();
        let (w1, w0) = (1.0 / (2.0 - c3), -c3 / (2.0 - c3));
        let weights: Vec<f64> = match order {
            | SplitOrder::Strang => vec![1.0],
            | SplitOrder::Yoshida4 => vec![w1, w0, w1],
        };
        let kins: Vec<Vec<Complex<f64>>> = weights
            .iter()
            .map(|w| self.kinetic_factors(w * dt, false))
            .collect();
        let mut t = t0;
        for _ in 0..steps {
            let mut tt = t;
            for (w, kin) in weights.iter().zip(&kins) {
                let vv: Vec<f64> = (0..psi.len())
                    .map(|p| v(tt + 0.5 * w * dt, self.x(p)))
                    .collect();
                self.strang(&mut psi, &vv, g, w * dt, false, kin, None);
                tt += w * dt;
            }
            t += dt;
        }
        psi
    }

    /// Real-time evolution with a complex absorbing potential `-i W(x)` (`W >= 0`): after every
    /// potential half step the amplitude is multiplied by `exp(-W dt / 2)`.
    #[allow(clippy::too_many_arguments)] // same controls as `evolve` plus the absorber
    pub fn evolve_absorbing(
        &self,
        psi0: &[Complex<f64>],
        v: &dyn Fn(f64, [f64; 3]) -> f64,
        absorber: &dyn Fn([f64; 3]) -> f64,
        g: f64,
        t0: f64,
        dt: f64,
        steps: usize,
        order: SplitOrder,
    ) -> Vec<Complex<f64>> {
        let mut psi = psi0.to_vec();
        let c3 = 2.0_f64.cbrt();
        let (w1, w0) = (1.0 / (2.0 - c3), -c3 / (2.0 - c3));
        let weights: Vec<f64> = match order {
            | SplitOrder::Strang => vec![1.0],
            | SplitOrder::Yoshida4 => vec![w1, w0, w1],
        };
        let kins: Vec<Vec<Complex<f64>>> = weights
            .iter()
            .map(|w| self.kinetic_factors(w * dt, false))
            .collect();
        let damp: Vec<f64> = (0..psi.len())
            .map(|p| absorber(self.x(p)).max(0.0))
            .collect();
        let mut t = t0;
        for _ in 0..steps {
            let mut tt = t;
            for (w, kin) in weights.iter().zip(&kins) {
                let vv: Vec<f64> = (0..psi.len())
                    .map(|p| v(tt + 0.5 * w * dt, self.x(p)))
                    .collect();
                self.strang(&mut psi, &vv, g, w * dt, false, kin, Some(&damp));
                tt += w * dt;
            }
            t += dt;
        }
        psi
    }

    /// Ground state by imaginary-time propagation with renormalisation.
    ///
    /// Stops when the energy changes by less than `tol` over a step or after
    /// `max_steps` steps.
    pub fn ground_state(
        &self,
        guess: &[Complex<f64>],
        v: &dyn Fn(f64, [f64; 3]) -> f64,
        g: f64,
        dt: f64,
        max_steps: usize,
        tol: f64,
    ) -> GroundState {
        let mut psi = guess.to_vec();
        let renorm = |psi: &mut [Complex<f64>]| {
            let s = self.norm2(psi).sqrt();
            for c in psi.iter_mut() {
                *c /= s;
            }
        };
        renorm(&mut psi);
        let vv: Vec<f64> = (0..psi.len()).map(|p| v(0.0, self.x(p))).collect();
        let kin = self.kinetic_factors(dt, true);
        let mut e_old = self.energy(&psi, v, 0.0, g);
        let mut steps = 0;
        let mut converged = false;
        while steps < max_steps {
            self.strang(&mut psi, &vv, g, dt, true, &kin, None);
            renorm(&mut psi);
            steps += 1;
            if steps % 10 == 0 {
                let e = self.energy(&psi, v, 0.0, g);
                if (e - e_old).abs() < tol {
                    converged = true;
                    e_old = e;
                    break;
                }
                e_old = e;
            }
        }
        let mu = self.chemical_potential(&psi, v, 0.0, g);
        GroundState {
            psi,
            energy: e_old,
            mu,
            steps,
            converged,
        }
    }
}
