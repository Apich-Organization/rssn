//! Riemann's method for the Euler–Poisson–Darboux equation of the hodograph
//! plane.
//!
//! The hodograph equation of a coupled affine system
//! `λ = a r + b s + c`, `μ = b r + a s + c` is the Euler–Poisson–Darboux
//! equation `E(N, N)`,
//!
//! ```text
//! (r - s) t_rs + N (t_s - t_r) = 0,      N = a / (a - b),
//! ```
//!
//! whose general solution is elementary for an integer `N` (Darboux,
//! `pde_hodograph_solution`) but not for the non-integer values that occur for
//! real gases (`N = 3/2` for shallow water, `γ = 2`). Riemann's method works
//! for every `N`. Written as `t_rs + a t_r + b t_s = 0` with
//! `a = -N/(r - s)`, `b = N/(r - s)`, the adjoint equation
//! `(ρ - σ)² R_ρσ - N (ρ - σ)(R_σ - R_ρ) - 2 N R = 0` has, for the field
//! point `(r, s)`, the **Riemann function**
//!
//! ```text
//! R(ρ, σ; r, s) = ((ρ - σ)/(r - s))^N  2F1(N, 1 - N; 1; z),
//! z = (ρ - r)(σ - s) / ((r - s)(ρ - σ)),
//! ```
//!
//! which equals `1` at `(ρ, σ) = (r, s)` and satisfies
//! `R_σ = a R` on `ρ = r` and `R_ρ = b R` on `σ = s`
//! (`R = ((r - σ)/(r - s))^N` and `((ρ - s)/(r - s))^N` there). For an
//! integer `N` the series terminates: `2F1(N, 1 - N; 1; z) = P_(N-1)(1 - 2z)`
//! is a Legendre polynomial. In the characteristic rectangle
//! `r0 <= ρ <= r`, `s0 <= σ <= s` with `r0 > s` (the region `ρ > σ` where the
//! system is hyperbolic) `0 <= z < 1`.
//!
//! [`epd_goursat`] evaluates the solution of the Goursat problem
//! `t(ρ, s0) = φ(ρ)`, `t(r0, σ) = ψ(σ)` by Green's identity,
//!
//! ```text
//! t(r, s) = R(r, s0) φ(r) + R(r0, s) ψ(s) - R(r0, s0) t(r0, s0)
//!         + ∫_{s0}^{s} ψ(σ) (a R - R_σ)(r0, σ) dσ
//!         + ∫_{r0}^{r} φ(ρ) (b R - R_ρ)(ρ, s0) dρ,
//! ```
//!
//! with the Riemann function and its derivatives in closed form
//! ([`epd_riemann`], [`epd_riemann_gradient`]).

use crate::kernels::hypergeometric::hyp2f1_cplx;
use crate::kernels::quadrature::integrate;
use num_complex::Complex64;

fn gauss(
    a: f64,
    b: f64,
    c: f64,
    z: f64,
) -> f64 {
    hyp2f1_cplx(
        Complex64::new(a, 0.0),
        Complex64::new(b, 0.0),
        Complex64::new(c, 0.0),
        Complex64::new(z, 0.0),
    )
    .re
}

/// The Riemann function `R(ρ, σ; r, s)` of `(r - s) t_rs + N (t_s - t_r) = 0`
/// for the field point `(r, s)` (`r > s`), evaluated at `(ρ, σ)` with
/// `ρ > σ`; `NaN` outside that region.
#[must_use]
pub fn epd_riemann(
    n: f64,
    r: f64,
    s: f64,
    rho: f64,
    sigma: f64,
) -> f64 {
    if !(r > s && rho > sigma) {
        return f64::NAN;
    }
    let a = (rho - sigma) / (r - s);
    let z = (rho - r) * (sigma - s) / ((r - s) * (rho - sigma));
    a.powf(n) * gauss(n, 1.0 - n, 1.0, z)
}

/// `(R_ρ, R_σ)` of [`epd_riemann`] in closed form (the derivative of
/// `2F1(N, 1-N; 1; z)` is `N (1-N) 2F1(N+1, 2-N; 2; z)`).
#[must_use]
pub fn epd_riemann_gradient(
    n: f64,
    r: f64,
    s: f64,
    rho: f64,
    sigma: f64,
) -> (f64, f64) {
    if !(r > s && rho > sigma) {
        return (f64::NAN, f64::NAN);
    }
    let a = (rho - sigma) / (r - s);
    let d = rho - sigma;
    let z = (rho - r) * (sigma - s) / ((r - s) * d);
    let f = gauss(n, 1.0 - n, 1.0, z);
    let fp = n * (1.0 - n) * gauss(n + 1.0, 2.0 - n, 2.0, z);
    let z_rho = (sigma - s) * (r - sigma) / ((r - s) * d * d);
    let z_sigma = (rho - r) * (rho - s) / ((r - s) * d * d);
    let scale = a.powf(n);
    (
        scale * (n / d * f + fp * z_rho),
        scale * (-n / d * f + fp * z_sigma),
    )
}

/// The solution `t(r, s)` of `(r - s) t_rs + N (t_s - t_r) = 0` of the Goursat
/// problem `t(ρ, s0) = φ(ρ)`, `t(r0, σ) = ψ(σ)` (`φ(r0) = ψ(s0)`), by
/// Riemann's method; see the module documentation.
///
/// The characteristic rectangle must lie in the hyperbolic region
/// `ρ > σ`, i.e. `r0 > s` (and `r >= r0`, `s >= s0`); `None` otherwise.
/// Accurate to the quadrature tolerance (about `1e-11`).
#[must_use]
pub fn epd_goursat(
    n: f64,
    (r0, s0): (f64, f64),
    phi: &dyn Fn(f64) -> f64,
    psi: &dyn Fn(f64) -> f64,
    (r, s): (f64, f64),
) -> Option<f64> {
    if !(r0 > s && r >= r0 && s >= s0 && n.is_finite()) {
        return None;
    }
    let boundary = f64::midpoint(phi(r0), psi(s0));
    let corner = epd_riemann(n, r, s, r0, s0);
    // R on the characteristics through (r, s) in closed form
    let along_rho = ((r - s0) / (r - s)).powf(n); // R(r, s0)
    let along_sigma = ((r0 - s) / (r - s)).powf(n); // R(r0, s)
    let mut total = along_rho * phi(r) + along_sigma * psi(s) - corner * boundary;
    if s > s0 {
        total += integrate(
            |sigma| {
                let a = -n / (r0 - sigma);
                let (_, r_sigma) = epd_riemann_gradient(n, r, s, r0, sigma);
                psi(sigma) * (a * epd_riemann(n, r, s, r0, sigma) - r_sigma)
            },
            s0,
            s,
            1e-12,
        )
        .value;
    }
    if r > r0 {
        total += integrate(
            |rho| {
                let b = n / (rho - s0);
                let (r_rho, _) = epd_riemann_gradient(n, r, s, rho, s0);
                phi(rho) * (b * epd_riemann(n, r, s, rho, s0) - r_rho)
            },
            r0,
            r,
            1e-12,
        )
        .value;
    }
    Some(total)
}
