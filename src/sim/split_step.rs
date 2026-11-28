//! Split-step Fourier Schrodinger/GPE solver (1D-3D, imaginary time).
//!
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
//! [`LinearSplitPlan`] prepares the phases for repeated real-time propagation
//! in a fixed real potential. It accepts periodic grids with arbitrary positive
//! extents, including odd sizes.
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
//! A **rotating frame** `H - Ω L_z` (`L_z = x p_y - y p_x`, rotation in the `x`-`y` plane about the
//! box centre, 2D and 3D) is available through [`SplitStep::evolve_rotating`]: the operator
//! `exp(i Ω dt L_z)` is a rotation of the wave function, applied exactly (to the grid's spectral
//! accuracy) as three alternating *shears* (Paeth): each shear is a shift of every grid line by
//! an amount proportional to the other coordinate, i.e. a multiplication by
//! `exp(i k a y)` in the 1D Fourier space of the line. All three are unitary, so the norm is
//! conserved to round-off, and the rotation commutes with the kinetic term, so a free packet
//! or an isotropic trap is rotated exactly. The rotation is sandwiched symmetrically around the
//! Strang step (second order; Yoshida weights for fourth order of the non-rotating part).
//!
//! Limitations: periodic boundaries (use the absorbing potential for open problems); the rotating
//! frame needs a wave function that vanishes near the box edges (the shears wrap around the
//! periodic box, and a rotated corner would leave it) and each sub-step rotates by at most
//! `π/2`; no dipolar interactions.

use crate::kernels::fft::Norm;
use crate::kernels::fft::transform_axes_in_place;
use crate::kernels::numeric_common::scale_product;
use crate::sim::spectral::SpectralGrid;
use num_complex::Complex;
use std::f64::consts::TAU;

/// Splitting scheme.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SplitOrder {
    /// Strang splitting (order 2).
    Strang,
    /// Yoshida composition of three Strang steps (order 4).
    Yoshida4,
}

/// Prepared real-time split-step propagation for the linear Schrödinger equation.
///
/// The plan owns the phases for one fixed real potential and a signed time
/// step. Each evolution reuses those phases and the shared FFT plans. Wave
/// functions and potentials are flattened with the x axis varying fastest.
#[derive(Debug)]
pub struct LinearSplitPlan {
    shape: [usize; 3],
    dim: usize,
    size: usize,
    phases: Option<LinearPhases>,
}

#[derive(Debug)]
struct LinearPhases {
    kinetic: Vec<Complex<f64>>,
    half_potential: Vec<Complex<f64>>,
    step_potential: Vec<Complex<f64>>,
    middle: Option<YoshidaPhases>,
}

#[derive(Debug)]
struct YoshidaPhases {
    kinetic: Vec<Complex<f64>>,
    bridge_potential: Vec<Complex<f64>>,
}

impl LinearSplitPlan {
    /// Prepares a fixed real potential on a one- to three-dimensional periodic grid.
    ///
    /// Each positive extent in `n` corresponds to a finite positive box length
    /// and the mass must be finite and positive. The potential is sampled at
    /// every grid point. `dt` may be negative for inverse time evolution or zero
    /// for the identity.
    ///
    /// # Errors
    /// Invalid grid dimensions, lengths, mass, time step or potential, a phase
    /// angle outside the representable range, or failure to reserve storage.
    pub fn new(
        n: &[usize],
        lengths: &[f64],
        mass: f64,
        potential: &[f64],
        dt: f64,
        order: SplitOrder,
    ) -> Result<Self, String> {
        let dim = n.len();
        if !(1..=3).contains(&dim)
            || lengths.len() != dim
            || n.contains(&0)
            || lengths.iter().any(|&len| !len.is_finite() || len <= 0.0)
            || !mass.is_finite()
            || mass <= 0.0
            || !dt.is_finite()
        {
            return Err("invalid linear split-step grid, mass or time step".into());
        }
        let size = n
            .iter()
            .try_fold(1_usize, |size, &extent| size.checked_mul(extent))
            .ok_or("linear split-step grid size overflows")?;
        if potential.len() != size || potential.iter().any(|value| !value.is_finite()) {
            return Err("linear split-step needs one finite potential value per point".into());
        }

        let mut shape = [1; 3];
        shape[..dim].copy_from_slice(n);
        shape[..dim].reverse();
        if dt == 0.0 {
            return Ok(Self {
                shape,
                dim,
                size,
                phases: None,
            });
        }
        let (outer, middle) = match order {
            | SplitOrder::Strang => (1.0, None),
            | SplitOrder::Yoshida4 => {
                let root = 2.0_f64.cbrt();
                (1.0 / (2.0 - root), Some(-root / (2.0 - root)))
            },
        };
        let half_potential = Self::potential_phases(potential, dt, outer * 0.5)?;
        // Two adjoining half steps form a single complex multiplication.
        let mut step_potential = Self::reserve(size)?;
        step_potential.extend(half_potential.iter().map(|&phase| phase * phase));
        let phases = LinearPhases {
            kinetic: Self::kinetic_phases(n, lengths, size, mass, dt, outer)?,
            half_potential,
            step_potential,
            middle: if let Some(weight) = middle {
                Some(YoshidaPhases {
                    kinetic: Self::kinetic_phases(n, lengths, size, mass, dt, weight)?,
                    bridge_potential: Self::potential_phases(
                        potential,
                        dt,
                        f64::midpoint(outer, weight),
                    )?,
                })
            } else {
                None
            },
        };
        Ok(Self {
            shape,
            dim,
            size,
            phases: Some(phases),
        })
    }

    fn reserve<T>(size: usize) -> Result<Vec<T>, String> {
        let mut buffer = Vec::new();
        buffer
            .try_reserve_exact(size)
            .map_err(|_| "linear split-step phase storage allocation failed")?;
        Ok(buffer)
    }

    fn potential_phases(
        potential: &[f64],
        dt: f64,
        weight: f64,
    ) -> Result<Vec<Complex<f64>>, String> {
        let mut phases = Self::reserve(potential.len())?;
        let sign = -dt.signum() * weight.signum();
        for &value in potential {
            let angle = scale_product(value, &[dt.abs(), weight.abs()], &[])
                .map_err(|_| "linear split-step potential phase is not representable")?;
            phases.push(Complex::from_polar(1.0, sign * angle));
        }
        Ok(phases)
    }

    #[allow(
        clippy::single_call_fn,
        reason = "Build kinetic phases from separable Fourier axes."
    )]
    fn kinetic_phases(
        n: &[usize],
        lengths: &[f64],
        size: usize,
        mass: f64,
        dt: f64,
        weight: f64,
    ) -> Result<Vec<Complex<f64>>, String> {
        let mut axes = [vec![0.0], vec![0.0], vec![0.0]];
        let error = "linear split-step kinetic phase is not representable";
        for ((axis, &extent), &len) in axes.iter_mut().zip(n).zip(lengths) {
            *axis = Self::reserve(extent)?;
            for j in 0..extent {
                let angular_mode = TAU * j.min(extent - j) as f64;
                let angle = scale_product(
                    0.5 * angular_mode * angular_mode,
                    &[dt.abs(), weight.abs()],
                    &[mass, len, len],
                )
                .map_err(|_| error)?;
                axis.push(angle);
            }
        }
        let [x, y, z] = axes;
        let sign = -dt.signum() * weight.signum();
        let mut phases = Self::reserve(size)?;
        for z in z {
            for &y in &y {
                for &x in &x {
                    let angle = x + y + z;
                    if !angle.is_finite() {
                        return Err(error.into());
                    }
                    phases.push(Complex::from_polar(1.0, sign * angle));
                }
            }
        }
        Ok(phases)
    }

    /// Number of points in the periodic grid.
    #[must_use]
    pub const fn size(&self) -> usize {
        self.size
    }

    pub(crate) const fn heap_bytes(&self) -> usize {
        let Some(phases) = &self.phases else {
            return 0;
        };
        let middle = if let Some(middle) = &phases.middle {
            middle.kinetic.capacity() + middle.bridge_potential.capacity()
        } else {
            0
        };
        (phases.kinetic.capacity()
            + phases.half_potential.capacity()
            + phases.step_potential.capacity()
            + middle)
            * size_of::<Complex<f64>>()
    }

    fn apply_phase(
        psi: &mut [Complex<f64>],
        phases: &[Complex<f64>],
    ) {
        for (value, &phase) in psi.iter_mut().zip(phases) {
            *value *= phase;
        }
    }

    fn apply_kinetic(
        &self,
        psi: &mut [Complex<f64>],
        phases: &[Complex<f64>],
    ) {
        let shape = &self.shape[..self.dim];
        let axes = &[2, 1, 0][3 - self.dim..];
        transform_axes_in_place(psi, shape, axes, false, Norm::Backward);
        Self::apply_phase(psi, phases);
        transform_axes_in_place(psi, shape, axes, true, Norm::Backward);
    }

    /// Evolves one wave function for `steps` steps, with no resampling of the potential.
    ///
    /// A zero step count or a zero time step returns the original samples exactly.
    ///
    /// # Errors
    /// A wrong input length, a non-finite sample, or a nonrepresentable wave function.
    pub fn evolve(
        &self,
        psi0: &[Complex<f64>],
        steps: usize,
    ) -> Result<Vec<Complex<f64>>, String> {
        if psi0.len() != self.size() {
            return Err("linear split-step wave function does not match the grid".into());
        }
        let mut peak = 0.0_f64;
        for value in psi0 {
            if !value.re.is_finite() || !value.im.is_finite() {
                return Err("linear split-step wave function must be finite".into());
            }
            peak = peak.max(value.re.abs()).max(value.im.abs());
        }
        let mut psi = psi0.to_vec();
        let Some(phases) = &self.phases else {
            return Ok(psi);
        };
        if steps == 0 || peak == 0.0 {
            return Ok(psi);
        }

        // FFT intermediates can be much larger than the physical samples.
        // Power-of-two scaling also preserves headroom for subnormal inputs.
        let scale = (!(1e-140..=1e140).contains(&peak))
            .then(|| f64::from_bits(peak.to_bits() & 0x7ff0_0000_0000_0000).max(f64::MIN_POSITIVE));
        if let Some(scale) = scale {
            for value in &mut psi {
                *value /= scale;
            }
        }
        Self::apply_phase(&mut psi, &phases.half_potential);
        for step in 1..=steps {
            self.apply_kinetic(&mut psi, &phases.kinetic);
            if let Some(middle) = &phases.middle {
                Self::apply_phase(&mut psi, &middle.bridge_potential);
                self.apply_kinetic(&mut psi, &middle.kinetic);
                Self::apply_phase(&mut psi, &middle.bridge_potential);
                self.apply_kinetic(&mut psi, &phases.kinetic);
            }
            let following = if step == steps {
                &phases.half_potential
            } else {
                &phases.step_potential
            };
            Self::apply_phase(&mut psi, following);
        }
        if let Some(scale) = scale {
            for value in &mut psi {
                *value *= scale;
            }
        }
        if psi
            .iter()
            .any(|value| !value.re.is_finite() || !value.im.is_finite())
        {
            return Err("evolved linear split-step wave function is not representable".into());
        }
        Ok(psi)
    }
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

    /// Prepares a fixed real potential sampled at the points of this box.
    ///
    /// # Errors
    /// See [`LinearSplitPlan::new`].
    pub fn prepare_linear(
        &self,
        potential: &[f64],
        dt: f64,
        order: SplitOrder,
    ) -> Result<LinearSplitPlan, String> {
        let dim = self.grid.dim;
        if dim > self.grid.n.len() {
            return Err("invalid linear split-step grid dimension".into());
        }
        LinearSplitPlan::new(
            &self.grid.n[..dim],
            &self.grid.len[..dim],
            self.mass,
            potential,
            dt,
            order,
        )
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

    /// `∫ |ψ|²`, with rescaling when necessary to keep the integrated norm in range.
    #[must_use]
    pub fn norm2(
        &self,
        psi: &[Complex<f64>],
    ) -> f64 {
        let squared = psi.iter().map(Complex::norm_sqr).sum::<f64>();
        let cell = self.cell();
        let direct = squared * cell;
        if squared.is_normal() && cell.is_normal() && direct.is_normal() {
            return direct;
        }
        self.scaled_norm2(psi, direct)
    }

    #[allow(
        clippy::single_call_fn,
        reason = "Keep the ordinary norm's arithmetic in its fast path."
    )]
    fn scaled_norm2(
        &self,
        psi: &[Complex<f64>],
        direct: f64,
    ) -> f64 {
        let dim = self.grid.dim;
        let lengths = &self.grid.len[..dim];
        let n = &self.grid.n[..dim];
        if lengths.iter().any(|&len| !len.is_finite() || len <= 0.0) || n.contains(&0) {
            return direct;
        }
        let mut peak = 0.0_f64;
        for value in psi {
            if !value.re.is_finite() || !value.im.is_finite() {
                return direct;
            }
            peak = peak.max(value.re.abs()).max(value.im.abs());
        }
        if peak == 0.0 {
            return 0.0;
        }
        let squared = psi
            .iter()
            .map(|value| (*value / peak).norm_sqr())
            .sum::<f64>();
        let mut numerator = [peak; 5];
        numerator[2..2 + dim].copy_from_slice(lengths);
        let denominator = self.grid.n.map(|value| value as f64);
        scale_product(squared, &numerator[..2 + dim], &denominator[..dim]).unwrap_or(f64::INFINITY)
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

    /// Applies the shear `ψ(x_a) -> ψ(x_a + shift * x_b)` along axis `a` (physical coordinates
    /// relative to the box centre): every line along `a` is shifted by a Fourier phase.
    fn shear(
        &self,
        psi: &mut [Complex<f64>],
        a: usize,
        b: usize,
        shift: f64,
    ) {
        use crate::kernels::transforms::fft_slice;
        use crate::kernels::transforms::ifft_slice;
        let n = self.grid.n;
        let stride = [1, n[0], n[0] * n[1]];
        let len = n[a];
        let mut line = vec![Complex::new(0.0, 0.0); len];
        let phases_k: Vec<f64> = (0..len).map(|j| self.grid.wavenumber(a, j)).collect();
        for p in 0..psi.len() {
            let idx = [p % n[0], (p / n[0]) % n[1], p / (n[0] * n[1])];
            if idx[a] != 0 {
                continue;
            }
            let xb = self.origin[b] + idx[b] as f64 * self.grid.h(b);
            for (j, l) in line.iter_mut().enumerate() {
                *l = psi[p + j * stride[a]];
            }
            fft_slice(&mut line);
            for (l, k) in line.iter_mut().zip(&phases_k) {
                *l *= Complex::new(0.0, k * shift * xb).exp();
            }
            ifft_slice(&mut line);
            for (j, l) in line.iter().enumerate() {
                psi[p + j * stride[a]] = *l;
            }
        }
    }

    /// Rotates the pattern of `psi` counterclockwise by `angle` (`|angle| <= π/2`) in the
    /// `x`-`y` plane about the box centre, as three shears.
    ///
    /// # Panics
    /// For a one-dimensional grid.
    pub fn rotate(
        &self,
        psi: &mut [Complex<f64>],
        angle: f64,
    ) {
        assert!(
            self.grid.dim >= 2,
            "a rotation needs at least two dimensions"
        );
        let a = (0.5 * angle).tan();
        let b = -angle.sin();
        self.shear(psi, 0, 1, a);
        self.shear(psi, 1, 0, b);
        self.shear(psi, 0, 1, a);
    }

    /// Real-time evolution in a frame rotating with angular velocity `omega` about the `z` axis:
    /// `i ψ_t = (-Δ/2m + V + g|ψ|² - Ω L_z) ψ`. With `omega = 0` it is exactly [`SplitStep::evolve`].
    ///
    /// # Panics
    /// For a one-dimensional grid, or when a sub-step would rotate by more than `π/2`.
    #[allow(clippy::too_many_arguments)] // same controls as `evolve` plus the rotation rate
    pub fn evolve_rotating(
        &self,
        psi0: &[Complex<f64>],
        v: &dyn Fn(f64, [f64; 3]) -> f64,
        g: f64,
        omega: f64,
        t0: f64,
        dt: f64,
        steps: usize,
        order: SplitOrder,
    ) -> Vec<Complex<f64>> {
        if omega == 0.0 {
            return self.evolve(psi0, v, g, t0, dt, steps, order);
        }
        assert!(
            self.grid.dim >= 2,
            "a rotating frame needs at least two dimensions"
        );
        assert!(
            (omega * dt).abs() < std::f64::consts::FRAC_PI_2,
            "Omega dt must be below pi/2"
        );
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
                // exp(i Omega dt L_z) turns the pattern clockwise by Omega dt
                self.rotate(&mut psi, -0.5 * omega * w * dt);
                let vv: Vec<f64> = (0..psi.len())
                    .map(|p| v(tt + 0.5 * w * dt, self.x(p)))
                    .collect();
                self.strang(&mut psi, &vv, g, w * dt, false, kin, None);
                self.rotate(&mut psi, -0.5 * omega * w * dt);
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
