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
//! | Explicit transport (advection, diffusion, advection-diffusion, Burgers) | [`solve_transport_nd`] | frozen boundary values; advection is periodic with period `n - 2` |
//! | Crank-Nicolson heat, factorised ADI (Beam-Warming) | [`solve_heat_adi_nd`] | zero Dirichlet |
//! | Schrodinger, Strang split with Cayley (Crank-Nicolson) kinetic steps | [`solve_schrodinger_cn_nd`] | zero Dirichlet |
//! | Wave, explicit leapfrog | [`solve_wave_leapfrog_nd`] | zero Dirichlet |
//!
//! Every axis needs at least three nodes. All schemes reduce, for a single
//! axis, to the corresponding 1D routine of the legacy modules.

use num_complex::Complex;

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
    let s = [1, n[0], n[0] * n[1]];
    Ok(Lay {
        n,
        s,
        dim: dims.len(),
        len: n[0] * n[1] * n[2],
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
        (0..self.len)
            .filter(|&p| {
                let q = self.ijk(p);
                q[axis] == 0
                    && (0..self.dim)
                        .filter(|&a| a != axis)
                        .all(|a| q[a] >= 1 && q[a] + 1 < self.n[a])
            })
            .collect()
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
/// For inconsistent shapes or a velocity list of the wrong length.
pub fn solve_transport_nd(
    u0: &[f64],
    dims: &[usize],
    h: &[f64],
    eq: &Transport,
    dt: f64,
    steps: usize,
) -> Result<Vec<Vec<f64>>, String> {
    run_transport(u0, dims, h, eq, dt, steps, true)
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
    let mut all = run_transport(u0, dims, h, eq, dt, steps, false)?;
    all.pop().ok_or_else(|| "no state".to_owned())
}

fn run_transport(
    u0: &[f64],
    dims: &[usize],
    h: &[f64],
    eq: &Transport,
    dt: f64,
    steps: usize,
    keep: bool,
) -> Result<Vec<Vec<f64>>, String> {
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
    let nodes = l.interior_nodes();
    let periodic = matches!(eq, Transport::Advection { .. });
    let mut u = u0.to_vec();
    let mut out = Vec::with_capacity(steps + 1);
    out.push(u.clone());
    for _ in 0..steps {
        let mut next = if periodic {
            vec![0.0; l.len]
        } else {
            u.clone()
        };
        for &p in &nodes {
            let up = u[p];
            let mut conv = 0.0;
            let mut diff = 0.0;
            let mut visc = 0.0;
            for a in 0..l.dim {
                let (fwd, bwd) = (u[p + l.s[a]], u[p - l.s[a]]);
                let second = (fwd - 2.0 * up + bwd) / (h[a] * h[a]);
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
            next[p] = up + dt * (visc * diff - conv);
        }
        if periodic {
            for a in 0..l.dim {
                for p in 0..l.len {
                    let q = l.ijk(p);
                    if q[a] == 0 {
                        next[p] = next[p + (l.n[a] - 2) * l.s[a]];
                    } else if q[a] == l.n[a] - 1 {
                        next[p] = next[p - (l.n[a] - 2) * l.s[a]];
                    }
                }
            }
        }
        u = next;
        if keep {
            out.push(u.clone());
        }
    }
    if !keep {
        out.clear();
        out.push(u);
    }
    Ok(out)
}

// ---------------------------------------------------------------------------
// Heat: factorised Crank-Nicolson ADI
// ---------------------------------------------------------------------------

/// In-place Thomas solve for a constant tridiagonal matrix (`off`, `diag`, `off`).
fn thomas_const(
    off: f64,
    diag: f64,
    d: &mut [f64],
    scratch: &mut Vec<f64>,
) {
    let m = d.len();
    scratch.clear();
    scratch.resize(m, 0.0);
    let mut denom = diag;
    scratch[0] = off / denom;
    d[0] /= denom;
    for i in 1..m {
        denom = diag - off * scratch[i - 1];
        scratch[i] = off / denom;
        d[i] = (d[i] - off * d[i - 1]) / denom;
    }
    for i in (0..m - 1).rev() {
        d[i] -= scratch[i] * d[i + 1];
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
/// For inconsistent shapes or non-positive `dt`/`d_coeff`.
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
    let rho: Vec<f64> = (0..l.dim)
        .map(|a| d_coeff * dt / (2.0 * h[a] * h[a]))
        .collect();
    let lines: Vec<Vec<usize>> = (0..l.dim).map(|a| l.lines(a)).collect();
    let mut buf = Vec::new();
    let mut scratch = Vec::new();
    for _ in 0..steps {
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
                let mut d: Vec<f64> = (1..l.n[a] - 1).map(|t| u[base + t * l.s[a]]).collect();
                thomas_const(-rho[a], 1.0 + 2.0 * rho[a], &mut d, &mut scratch);
                for (t, v) in d.into_iter().enumerate() {
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

fn thomas_complex_const(
    off: Complex<f64>,
    diag: Complex<f64>,
    d: &mut [Complex<f64>],
) {
    let m = d.len();
    let mut cp = vec![Complex::new(0.0, 0.0); m];
    let mut denom = diag;
    cp[0] = off / denom;
    d[0] /= denom;
    for i in 1..m {
        denom = diag - off * cp[i - 1];
        cp[i] = off / denom;
        d[i] = (d[i] - off * d[i - 1]) / denom;
    }
    for i in (0..m - 1).rev() {
        let next = d[i + 1];
        d[i] -= cp[i] * next;
    }
}

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
/// For inconsistent shapes or non-positive `dt`.
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
    let phase: Vec<Complex<f64>> = v
        .iter()
        .map(|&x| Complex::from_polar(1.0, -0.5 * dt * x))
        .collect();
    let lines: Vec<Vec<usize>> = (0..l.dim).map(|a| l.lines(a)).collect();
    let mut buf: Vec<Complex<f64>> = Vec::new();
    for _ in 0..steps {
        for (x, ph) in psi.iter_mut().zip(&phase) {
            *x *= *ph;
        }
        for a in 0..l.dim {
            let q = Complex::new(0.0, dt / (4.0 * h[a] * h[a]));
            let e = Complex::new(1.0, -dt / (2.0 * h[a] * h[a]));
            let d = Complex::new(1.0, dt / (2.0 * h[a] * h[a]));
            for &base in &lines[a] {
                buf.clear();
                buf.extend((1..l.n[a] - 1).map(|t| {
                    let c = psi[base + t * l.s[a]];
                    e * c + q * (psi[base + (t + 1) * l.s[a]] + psi[base + (t - 1) * l.s[a]])
                }));
                thomas_complex_const(-q, d, &mut buf);
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
/// For inconsistent shapes or a violated CFL condition
/// `|c| dt sqrt(Σ 1/h_a²) <= 1`.
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
    let inv2: Vec<f64> = (0..l.dim).map(|a| 1.0 / (h[a] * h[a])).collect();
    let cfl = c.abs() * dt * inv2.iter().sum::<f64>().sqrt();
    if !(dt.is_finite() && dt > 0.0) || cfl > 1.0 {
        return Err(format!("CFL violation: |c| dt sqrt(sum 1/h^2) = {cfl} > 1"));
    }
    let nodes = l.interior_nodes();
    let laplacian = |u: &[f64], p: usize| -> f64 {
        (0..l.dim)
            .map(|a| inv2[a] * (u[p + l.s[a]] - 2.0 * u[p] + u[p - l.s[a]]))
            .sum()
    };
    let mut prev = u0.to_vec();
    l.zero_boundary(&mut prev);
    let mut cur = vec![0.0; l.len];
    for &p in &nodes {
        cur[p] = prev[p] + dt * ut0[p] + 0.5 * c * c * dt * dt * laplacian(&prev, p);
    }
    let mut out = Vec::with_capacity(steps + 1);
    out.push(prev.clone());
    if steps == 0 {
        return Ok(out);
    }
    out.push(cur.clone());
    for _ in 2..=steps {
        let mut next = vec![0.0; l.len];
        for &p in &nodes {
            next[p] = 2.0 * cur[p] - prev[p] + c * c * dt * dt * laplacian(&cur, p);
        }
        out.push(next.clone());
        prev = cur;
        cur = next;
    }
    Ok(out)
}
