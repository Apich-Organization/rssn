//! Dimension-generic (1D/2D/3D) Poisson, transport, heat, Schrodinger and wave kernels.
//!
//! # Dimension-generic grid kernels (1D, 2D and 3D)
//!
//! The legacy solver modules (`physics_cfd`, `physics_fdm`, `physics_cnm`,
//! `classical`, ...) grew one function per space dimension. This module holds
//! the shared kernels that work on a flat, `x`-fastest array together with a
//! shape slice `dims = [nx]`, `[nx, ny]` or `[nx, ny, nz]` and the grid
//! spacings `h = [hx, hy, hz]` (as many as axes). The dimension-specific names
//! in the legacy modules are thin wrappers around these entry points.
//!
//! | Family | Entry point | Boundary treatment |
//! |---|---|---|
//! | Poisson `∇²u = f`, Jacobi / Gauss-Seidel / SOR | [`solve_poisson_iterative_nd`] | Dirichlet values taken from the initial guess |
//! | Explicit transport (advection, diffusion, advection-diffusion, Burgers) | [`solve_transport_nd`], [`solve_transport_sampled_nd`] | frozen boundary values; advection is periodic with period `n - 2` |
//! | Crank-Nicolson heat, factorised ADI (Beam-Warming) | [`solve_heat_adi_nd`] | zero Dirichlet |
//! | Schrodinger, Strang split with Cayley (Crank-Nicolson) kinetic steps | [`solve_schrodinger_cn_nd`] | zero Dirichlet |
//! | Wave, explicit leapfrog | [`solve_wave_leapfrog_nd`] | zero Dirichlet |
//!
//! Every axis needs at least three nodes. All schemes reduce, for a single
//! axis, to the corresponding 1D routine of the legacy modules.

use num_complex::Complex;
use std::num::NonZeroUsize;
use std::ops::Div;
use std::ops::Mul;
use std::ops::Sub;

use crate::kernels::numeric_common::trajectory_sample_count;

/// Shape bookkeeping for a flat array with `x` fastest.
#[derive(Clone, Copy, Debug)]
struct Lay {
    n: [usize; 3],
    s: [usize; 3],
    dim: usize,
    len: usize,
}

fn lay(
    dims: &[usize],
    h: &[f64],
) -> Result<Lay, String> {
    if dims.is_empty() || dims.len() > 3 {
        return Err("the grid must have 1, 2 or 3 axes".to_owned());
    }
    if h.len() != dims.len() {
        return Err("one grid spacing per axis is required".to_owned());
    }
    if dims.iter().any(|&n| n < 3) {
        return Err("every axis needs at least 3 nodes".to_owned());
    }
    if h.iter().any(|&v| !(v.is_finite() && v > 0.0)) {
        return Err("grid spacings must be finite and positive".to_owned());
    }
    let mut n = [1usize; 3];
    for (a, &d) in dims.iter().enumerate() {
        n[a] = d;
    }
    let plane = n[0]
        .checked_mul(n[1])
        .ok_or("the grid shape is too large")?;
    let len = plane
        .checked_mul(n[2])
        .ok_or("the grid shape is too large")?;
    Ok(Lay {
        n,
        s: [1, n[0], plane],
        dim: dims.len(),
        len,
    })
}

impl Lay {
    const fn ijk(
        &self,
        p: usize,
    ) -> [usize; 3] {
        [p % self.n[0], (p / self.s[1]) % self.n[1], p / self.s[2]]
    }

    fn is_interior(
        &self,
        q: [usize; 3],
    ) -> bool {
        (0..self.dim).all(|a| q[a] >= 1 && q[a] + 1 < self.n[a])
    }

    /// Flat indices of every node whose active coordinates are all interior.
    fn interior_nodes(&self) -> Vec<usize> {
        (0..self.len)
            .filter(|&p| self.is_interior(self.ijk(p)))
            .collect()
    }

    /// First node of every line along `axis` whose other coordinates are interior.
    fn lines(
        &self,
        axis: usize,
    ) -> Vec<usize> {
        let [x, y, z] = std::array::from_fn(|a| {
            if a == axis || a >= self.dim {
                0..1
            } else {
                1..self.n[a] - 1
            }
        });
        let mut bases = Vec::with_capacity(x.len() * y.len() * z.len());
        for k in z {
            for j in y.clone() {
                let base = j * self.s[1] + k * self.s[2];
                bases.extend(x.clone().map(|i| base + i));
            }
        }
        bases
    }

    fn zero_boundary<T: Copy + Default>(
        &self,
        u: &mut [T],
    ) {
        for p in 0..self.len {
            if !self.is_interior(self.ijk(p)) {
                u[p] = T::default();
            }
        }
    }

    fn fill_periodic_ghosts(
        &self,
        u: &mut [f64],
    ) {
        // Axis order also fills edges and corners from the wrapped interior.
        for a in 0..self.dim {
            let stride = self.s[a];
            let high = (self.n[a] - 1) * stride;
            for strip in u.chunks_exact_mut(self.n[a] * stride) {
                strip.copy_within(high - stride..high, 0);
                strip.copy_within(stride..2 * stride, high);
            }
        }
    }

    fn update_interior(
        &self,
        out: &mut [f64],
        update: impl Fn(usize) -> f64 + Send + Sync,
    ) {
        // Rows are independent. A long 1D row can also split into chunks.
        const PARALLEL_CELLS: usize = 32_768;
        let width = if self.dim == 1 {
            self.n[0].min(PARALLEL_CELLS / 4)
        } else {
            self.n[0]
        };
        let grain = PARALLEL_CELLS.div_ceil(width);
        crate::kernels::par::for_chunks_mut(out, width, grain, |chunk, row| {
            let base = chunk * width;
            let [x, y, z] = self.ijk(base);
            if (self.dim >= 2 && (y == 0 || y + 1 == self.n[1]))
                || (self.dim >= 3 && (z == 0 || z + 1 == self.n[2]))
            {
                return;
            }
            let start = usize::from(x == 0);
            let end = row.len() - usize::from(x + row.len() == self.n[0]);
            for (offset, result) in row.iter_mut().enumerate().take(end).skip(start) {
                *result = update(base + offset);
            }
        });
    }
}

// ---------------------------------------------------------------------------
// Poisson
// ---------------------------------------------------------------------------

/// Stationary iteration used by [`solve_poisson_iterative_nd`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum IterativeMethod {
    /// Jacobi (simultaneous update).
    Jacobi,
    /// Gauss-Seidel (lexicographic in-place update).
    GaussSeidel,
    /// Successive over-relaxation with the given factor `omega` (`1 < omega < 2`).
    Sor(f64),
}

/// Solves `∇²u = f` on a 1D, 2D or 3D grid by a stationary iteration.
///
/// `u0` is the initial guess; its boundary nodes are the (fixed) Dirichlet
/// data. Iteration stops when the largest update is below `tolerance` or after
/// `max_iter` sweeps. With a 1D grid this is the scalar tridiagonal iteration,
/// with `dims = [nx, ny]` it reproduces
/// [`crate::sim::physics_cfd::solve_poisson_2d_jacobi`] and friends.
///
/// # Errors
/// For inconsistent shapes, spacings or lengths.
pub fn solve_poisson_iterative_nd(
    f: &[f64],
    u0: &[f64],
    dims: &[usize],
    h: &[f64],
    method: IterativeMethod,
    max_iter: usize,
    tolerance: f64,
) -> Result<Vec<f64>, String> {
    let l = lay(dims, h)?;
    if f.len() != l.len || u0.len() != l.len {
        return Err("array length does not match the grid".to_owned());
    }
    let mut inv2 = [0.0; 3];
    for a in 0..l.dim {
        inv2[a] = 1.0 / (h[a] * h[a]);
    }
    let diag: f64 = (0..l.dim).map(|a| 2.0 * inv2[a]).sum();
    let nodes = l.interior_nodes();
    let mut u = u0.to_vec();
    let mut w = u0.to_vec();
    for _ in 0..max_iter {
        let mut max_diff = 0.0_f64;
        for &p in &nodes {
            let mut sum = 0.0;
            for a in 0..l.dim {
                sum += inv2[a] * (u[p + l.s[a]] + u[p - l.s[a]]);
            }
            let gs = (sum - f[p]) / diag;
            let old = u[p];
            match method {
                | IterativeMethod::Jacobi => {
                    w[p] = gs;
                    max_diff = max_diff.max((gs - old).abs());
                },
                | IterativeMethod::GaussSeidel => {
                    u[p] = gs;
                    max_diff = max_diff.max((gs - old).abs());
                },
                | IterativeMethod::Sor(omega) => {
                    let new = old + omega * (gs - old);
                    u[p] = new;
                    max_diff = max_diff.max((new - old).abs());
                },
            }
        }
        if matches!(method, IterativeMethod::Jacobi) {
            std::mem::swap(&mut u, &mut w);
        }
        if max_diff < tolerance {
            break;
        }
    }
    Ok(u)
}

// ---------------------------------------------------------------------------
// Transport
// ---------------------------------------------------------------------------

/// Explicit transport equation solved by [`solve_transport_nd`].
#[derive(Clone, Debug, PartialEq)]
pub enum Transport {
    /// `u_t + c·∇u = 0`, first-order upwind, periodic with period `n - 2`
    /// (node `0` is identified with node `n - 2`, node `n - 1` with node `1`).
    Advection {
        /// Velocity per axis.
        velocity: Vec<f64>,
    },
    /// `u_t = α Δu`, central differences, frozen boundary values.
    Diffusion {
        /// Diffusivity.
        alpha: f64,
    },
    /// `u_t + c·∇u = α Δu`, upwind + central, frozen boundary values.
    AdvectionDiffusion {
        /// Velocity per axis.
        velocity: Vec<f64>,
        /// Diffusivity.
        alpha: f64,
    },
    /// Scalar viscous Burgers `u_t + u (∂_1 + ... + ∂_d) u = ν Δu` with
    /// upwinded convection; in 1D this is `u_t + u u_x = ν u_xx`.
    Burgers {
        /// Viscosity.
        nu: f64,
    },
    /// As [`Transport::Burgers`] with central convection differences.
    BurgersCentral {
        /// Viscosity.
        nu: f64,
    },
}

/// Advances a [`Transport`] equation `steps` explicit Euler steps on a 1D, 2D
/// or 3D grid and returns all `steps + 1` snapshots (the first is `u0`).
///
/// The scheme is only conditionally stable: keep `Σ|c_a| dt/h_a <= 1` and
/// `2 α dt Σ 1/h_a² <= 1`.
///
/// # Errors
/// For inconsistent shapes, a velocity list of the wrong length, invalid time
/// controls, too many snapshots, or cancellation.
pub fn solve_transport_nd(
    u0: &[f64],
    dims: &[usize],
    h: &[f64],
    eq: &Transport,
    dt: f64,
    steps: usize,
) -> Result<Vec<Vec<f64>>, String> {
    solve_transport_sampled_nd(u0, dims, h, eq, dt, steps, 1)
}

/// Advances a [`Transport`] equation and retains snapshots at a fixed step cadence.
///
/// The result holds step 0, every positive multiple of `save_every` through
/// `steps`, and the final step if it is off cadence. The initial and final state
/// share one snapshot when `steps` is zero. Boundary and stability conventions
/// are those of [`solve_transport_nd`].
///
/// # Errors
/// For an inconsistent grid or equation, zero cadence, or an unrepresentable
/// snapshot count, or cancellation.
pub fn solve_transport_sampled_nd(
    u0: &[f64],
    dims: &[usize],
    h: &[f64],
    eq: &Transport,
    dt: f64,
    steps: usize,
    save_every: usize,
) -> Result<Vec<Vec<f64>>, String> {
    let count = transport_snapshot_count(steps, save_every)?;
    let mut history = Vec::new();
    history
        .try_reserve_exact(count)
        .map_err(|_| "too many transport snapshots")?;
    let final_state = run_transport(u0, dims, h, eq, dt, steps, |step, field| {
        if step < steps && step.is_multiple_of(save_every) {
            history.push(field.to_vec());
        }
    })?;
    history.push(final_state);
    Ok(history)
}

pub(crate) fn transport_snapshot_count(
    steps: usize,
    save_every: usize,
) -> Result<usize, &'static str> {
    let cadence = NonZeroUsize::new(save_every).ok_or("the snapshot cadence must be positive")?;
    trajectory_sample_count(steps, cadence).ok_or("too many transport snapshots")
}

/// Like [`solve_transport_nd`] but returns only the final state (no snapshot
/// history, so memory stays `O(n)`).
///
/// # Errors
/// As [`solve_transport_nd`].
pub fn advance_transport_nd(
    u0: &[f64],
    dims: &[usize],
    h: &[f64],
    eq: &Transport,
    dt: f64,
    steps: usize,
) -> Result<Vec<f64>, String> {
    run_transport(u0, dims, h, eq, dt, steps, |_, _| {})
}

fn run_transport(
    u0: &[f64],
    dims: &[usize],
    h: &[f64],
    eq: &Transport,
    dt: f64,
    steps: usize,
    mut snapshot: impl FnMut(usize, &[f64]),
) -> Result<Vec<f64>, String> {
    let l = lay(dims, h)?;
    if u0.len() != l.len {
        return Err("array length does not match the grid".to_owned());
    }
    if !(dt.is_finite() && dt > 0.0) {
        return Err("the time step must be positive".to_owned());
    }
    if let Transport::Advection { velocity } | Transport::AdvectionDiffusion { velocity, .. } = eq
        && velocity.len() != l.dim
    {
        return Err("one velocity component per axis is required".to_owned());
    }
    let periodic = matches!(eq, Transport::Advection { .. });
    let mut u = u0.to_vec();
    snapshot(0, &u);
    if steps == 0 {
        return Ok(u);
    }
    let mut next = u.clone();
    let spacing_squared: Vec<_> = h.iter().map(|spacing| spacing * spacing).collect();
    for step in 1..=steps {
        crate::cancellation::checkpoint().map_err(|error| error.to_string())?;
        l.update_interior(&mut next, |p| {
            let up = u[p];
            let mut conv = 0.0;
            let mut diff = 0.0;
            let mut visc = 0.0;
            for a in 0..l.dim {
                let (fwd, bwd) = (u[p + l.s[a]], u[p - l.s[a]]);
                let second = (fwd - 2.0 * up + bwd) / spacing_squared[a];
                let upwind = |v: f64| {
                    if v > 0.0 {
                        v * (up - bwd) / h[a]
                    } else {
                        v * (fwd - up) / h[a]
                    }
                };
                match eq {
                    | Transport::Advection { velocity } => conv += upwind(velocity[a]),
                    | Transport::Diffusion { alpha } => {
                        visc = *alpha;
                        diff += second;
                    },
                    | Transport::AdvectionDiffusion { velocity, alpha } => {
                        conv += upwind(velocity[a]);
                        visc = *alpha;
                        diff += second;
                    },
                    | Transport::Burgers { nu } => {
                        conv += upwind(up);
                        visc = *nu;
                        diff += second;
                    },
                    | Transport::BurgersCentral { nu } => {
                        conv += up * (fwd - bwd) / (2.0 * h[a]);
                        visc = *nu;
                        diff += second;
                    },
                }
            }
            up + dt * (visc * diff - conv)
        });
        if periodic {
            l.fill_periodic_ghosts(&mut next);
        }
        std::mem::swap(&mut u, &mut next);
        snapshot(step, &u);
    }
    Ok(u)
}

// ---------------------------------------------------------------------------
// Heat: factorised Crank-Nicolson ADI
// ---------------------------------------------------------------------------

/// Thomas factorization of a tridiagonal matrix with constant diagonals.
struct ThomasFactors<T> {
    off: T,
    rows: Vec<(T, T)>,
}

impl<T> ThomasFactors<T>
where
    T: Copy + Div<Output = T> + Mul<Output = T> + Sub<Output = T>,
{
    fn new(
        off: T,
        diag: T,
        len: usize,
    ) -> Self {
        let mut rows = Vec::with_capacity(len);
        let mut pivot = diag;
        for _ in 0..len {
            let upper = off / pivot;
            rows.push((pivot, upper));
            pivot = diag - off * upper;
        }
        Self { off, rows }
    }

    fn solve(
        &self,
        rhs: &mut [T],
    ) {
        rhs[0] = rhs[0] / self.rows[0].0;
        for i in 1..rhs.len() {
            rhs[i] = (rhs[i] - self.off * rhs[i - 1]) / self.rows[i].0;
        }
        for i in (0..rhs.len() - 1).rev() {
            rhs[i] = rhs[i] - self.rows[i].1 * rhs[i + 1];
        }
    }
}

/// Crank-Nicolson heat equation `u_t = D Δu` with zero Dirichlet boundaries on
/// a 1D, 2D or 3D grid, advanced `steps` steps of size `dt`.
///
/// The step is the factorised (Beam-Warming) approximate-factorisation form
/// `Π_a (I - ρ_a δ²_a) uⁿ⁺¹ = Π_a (I + ρ_a δ²_a) uⁿ` with `ρ_a = D dt / (2 h_a²)`:
/// unconditionally stable, second order in space and time, one tridiagonal solve
/// per line and axis. For data that is a product of functions of the single axes
/// it equals the product of the 1D Crank-Nicolson solutions. Boundary entries of
/// `u0` are ignored (treated as zero).
///
/// # Errors
/// For inconsistent shapes, non-positive `dt`, negative diffusivity, or an
/// unrepresentable diffusion coefficient at the given spacing and time step.
/// A cancelled computation also returns an error.
pub fn solve_heat_adi_nd(
    u0: &[f64],
    dims: &[usize],
    h: &[f64],
    d_coeff: f64,
    dt: f64,
    steps: usize,
) -> Result<Vec<f64>, String> {
    let l = lay(dims, h)?;
    if u0.len() != l.len {
        return Err("array length does not match the grid".to_owned());
    }
    if !(dt.is_finite() && dt > 0.0 && d_coeff.is_finite() && d_coeff >= 0.0) {
        return Err("dt must be positive and the diffusivity non-negative".to_owned());
    }
    let mut u = u0.to_vec();
    l.zero_boundary(&mut u);
    if steps == 0 || d_coeff == 0.0 {
        return Ok(u);
    }
    let rho: Vec<_> = h
        .iter()
        .map(|&spacing| {
            let numerator = d_coeff * dt;
            let denominator = 2.0 * spacing * spacing;
            if numerator.is_normal() && denominator.is_normal() {
                Ok(numerator / denominator)
            } else {
                crate::kernels::numeric_common::scale_product(
                    d_coeff,
                    &[dt],
                    &[2.0, spacing, spacing],
                )
                .map_err(|error| error.to_string())
            }
        })
        .collect::<Result<_, _>>()?;
    if rho.iter().any(|&r| !(1.0 + 2.0 * r).is_finite()) {
        return Err("heat ADI diagonal is not finite".to_owned());
    }
    let lines: Vec<Vec<usize>> = (0..l.dim).map(|a| l.lines(a)).collect();
    let factors: Vec<_> = rho
        .iter()
        .zip(dims)
        .map(|(&r, &n)| ThomasFactors::new(-r, 1.0 + 2.0 * r, n - 2))
        .collect();
    let mut buf = crate::kernels::reuse::f64_scratch(l.n.into_iter().max().unwrap_or(0));
    for _ in 0..steps {
        crate::cancellation::checkpoint().map_err(|error| error.to_string())?;
        for a in 0..l.dim {
            for &base in &lines[a] {
                buf.clear();
                buf.extend((0..l.n[a]).map(|t| u[base + t * l.s[a]]));
                for t in 1..l.n[a] - 1 {
                    u[base + t * l.s[a]] =
                        buf[t] + rho[a] * (buf[t + 1] - 2.0 * buf[t] + buf[t - 1]);
                }
            }
        }
        for a in 0..l.dim {
            for &base in &lines[a] {
                buf.clear();
                buf.extend((1..l.n[a] - 1).map(|t| u[base + t * l.s[a]]));
                factors[a].solve(&mut buf);
                for (t, &v) in buf.iter().enumerate() {
                    u[base + (t + 1) * l.s[a]] = v;
                }
            }
        }
    }
    Ok(u)
}

// ---------------------------------------------------------------------------
// Schrodinger
// ---------------------------------------------------------------------------

/// Schrodinger equation `i ψ_t = -½ Δψ + V ψ` (ħ = m = 1) with zero Dirichlet
/// boundaries on a 1D, 2D or 3D grid.
///
/// One step is the symmetric splitting `e^{-iV dt/2} · Π_a C_a · e^{-iV dt/2}`
/// where `C_a = (I + i dt/2 K_a)⁻¹ (I - i dt/2 K_a)` is the Crank-Nicolson
/// (Cayley) step of the one-dimensional kinetic operator `K_a = -½ δ²_a`. The
/// kinetic factors commute, every factor is unitary, so the scheme is second
/// order in space and time and conserves the norm to round-off. For `V = 0` and
/// separable data it equals the product of 1D Crank-Nicolson solutions.
///
/// # Errors
/// For inconsistent shapes, non-positive `dt`, or cancellation.
pub fn solve_schrodinger_cn_nd(
    psi0: &[Complex<f64>],
    v: &[f64],
    dims: &[usize],
    h: &[f64],
    dt: f64,
    steps: usize,
) -> Result<Vec<Complex<f64>>, String> {
    let l = lay(dims, h)?;
    if psi0.len() != l.len || v.len() != l.len {
        return Err("array length does not match the grid".to_owned());
    }
    if !(dt.is_finite() && dt > 0.0) {
        return Err("the time step must be positive".to_owned());
    }
    let mut psi = psi0.to_vec();
    l.zero_boundary(&mut psi);
    if steps == 0 {
        return Ok(psi);
    }
    let phase: Vec<Complex<f64>> = v
        .iter()
        .map(|&x| Complex::from_polar(1.0, -0.5 * dt * x))
        .collect();
    let lines: Vec<Vec<usize>> = (0..l.dim).map(|a| l.lines(a)).collect();
    let axes = h
        .iter()
        .zip(dims)
        .map(|(&spacing, &n)| {
            let q = Complex::new(0.0, dt / (4.0 * spacing * spacing));
            let e = Complex::new(1.0, -dt / (2.0 * spacing * spacing));
            let d = Complex::new(1.0, dt / (2.0 * spacing * spacing));
            (q, e, ThomasFactors::new(-q, d, n - 2))
        })
        .collect::<Vec<_>>();
    let mut buf = Vec::with_capacity(l.n.into_iter().max().unwrap_or(0));
    for _ in 0..steps {
        crate::cancellation::checkpoint().map_err(|error| error.to_string())?;
        for (x, ph) in psi.iter_mut().zip(&phase) {
            *x *= *ph;
        }
        for a in 0..l.dim {
            let (q, e, factors) = &axes[a];
            for &base in &lines[a] {
                buf.clear();
                buf.extend((1..l.n[a] - 1).map(|t| {
                    let c = psi[base + t * l.s[a]];
                    *e * c + *q * (psi[base + (t + 1) * l.s[a]] + psi[base + (t - 1) * l.s[a]])
                }));
                factors.solve(&mut buf);
                for (t, val) in buf.iter().enumerate() {
                    psi[base + (t + 1) * l.s[a]] = *val;
                }
            }
        }
        for (x, ph) in psi.iter_mut().zip(&phase) {
            *x *= *ph;
        }
    }
    Ok(psi)
}

// ---------------------------------------------------------------------------
// Wave
// ---------------------------------------------------------------------------

/// Explicit leapfrog for `u_tt = c² Δu` with zero Dirichlet boundaries on a
/// 1D, 2D or 3D grid. Returns the `steps + 1` snapshots (`u0`, `u1`, ...), the
/// first step using the initial velocity `ut0` (Taylor start).
///
/// # Errors
/// For inconsistent shapes, an unrepresentable snapshot count, a nonfinite
/// speed or invalid time step, cancellation, or a violated CFL condition
/// `|c| dt sqrt(Σ 1/h_a²) <= 1` (with a `10⁻¹²` rounding allowance).
pub fn solve_wave_leapfrog_nd(
    u0: &[f64],
    ut0: &[f64],
    dims: &[usize],
    h: &[f64],
    c: f64,
    dt: f64,
    steps: usize,
) -> Result<Vec<Vec<f64>>, String> {
    let l = lay(dims, h)?;
    if u0.len() != l.len || ut0.len() != l.len {
        return Err("array length does not match the grid".to_owned());
    }
    if !(c.is_finite() && dt.is_finite() && dt > 0.0) {
        return Err("wave speed must be finite and the time step positive".into());
    }
    let distance = c.abs() * dt;
    let mut weights = [0.0; 3];
    for (&spacing, weight) in h.iter().zip(&mut weights) {
        let courant = if distance.is_normal() {
            distance / spacing
        } else {
            crate::kernels::numeric_common::scale_product(c.abs(), &[dt], &[spacing])
                .map_err(|error| error.to_string())?
        };
        *weight = courant * courant;
    }
    let cfl = weights.into_iter().sum::<f64>().sqrt();
    if !cfl.is_finite() || cfl > 1.0 + 1e-12 {
        return Err(format!("CFL violation: |c| dt sqrt(sum 1/h^2) = {cfl} > 1"));
    }
    let stencil = |u: &[f64], p: usize| -> f64 {
        (0..l.dim)
            .map(|a| weights[a] * (u[p + l.s[a]] - 2.0 * u[p] + u[p - l.s[a]]))
            .sum()
    };
    let mut initial = u0.to_vec();
    l.zero_boundary(&mut initial);
    let mut out = Vec::new();
    out.try_reserve_exact(steps.checked_add(1).ok_or("too many wave snapshots")?)
        .map_err(|_| "too many wave snapshots")?;
    out.push(initial);
    for step in 1..=steps {
        crate::cancellation::checkpoint().map_err(|error| error.to_string())?;
        let current = &out[step - 1];
        let mut next = vec![0.0; l.len];
        if step == 1 {
            l.update_interior(&mut next, |p| {
                current[p] + dt * ut0[p] + 0.5 * stencil(current, p)
            });
        } else {
            let prior = &out[step - 2];
            l.update_interior(&mut next, |p| {
                2.0 * current[p] - prior[p] + stencil(current, p)
            });
        }
        out.push(next);
    }
    Ok(out)
}
