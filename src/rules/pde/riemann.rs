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
//! ([`epd_riemann`], [`epd_riemann_gradient`]). Sampled characteristic data
//! can be evaluated at many field points with [`epd_goursat_samples`].

use crate::kernels::exact_eq;
use crate::kernels::hypergeometric::hyp2f1_cplx;
use crate::kernels::numeric_common::NumericError;
use crate::kernels::par::try_map_indexed;
use crate::kernels::quadrature::GkRule;
use crate::kernels::quadrature::integrate;
use crate::kernels::quadrature::integrate_gk_breakpoints;
use crate::kernels::special_ext::hyp2f1;
use num_complex::Complex64;

fn gauss(
    a: f64,
    b: f64,
    c: f64,
    z: f64,
) -> f64 {
    // The characteristic rectangle has real 0 <= z < 1. Use the real
    // kernel there, including its logarithmic continuation near one.
    if (0.0..1.0).contains(&z) && a.abs().max(b.abs()) <= 32.0 {
        let real = hyp2f1(a, b, c, z);
        if real.is_finite() {
            return real;
        }
    }
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
    let width = r - s;
    let d = rho - sigma;
    let a = d / width;
    let z = ((rho - r) / width) * ((sigma - s) / d);
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
    let width = r - s;
    let d = rho - sigma;
    let a = d / width;
    let z = ((rho - r) / width) * ((sigma - s) / d);
    let f = gauss(n, 1.0 - n, 1.0, z);
    let fp = n * (1.0 - n) * gauss(n + 1.0, 2.0 - n, 2.0, z);
    let z_rho = ((sigma - s) / d) * ((r - sigma) / width) / d;
    let z_sigma = ((rho - r) / width) * ((rho - s) / d) / d;
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
/// `ρ > σ`, i.e. `r0 > s` with `r >= r0` and `s >= s0`. Returns `None`
/// for invalid geometry or a nonfinite result. In regular regions the
/// estimate is accurate to the quadrature tolerance (about `1e-11`).
#[must_use]
pub fn epd_goursat(
    n: f64,
    (r0, s0): (f64, f64),
    phi: &dyn Fn(f64) -> f64,
    psi: &dyn Fn(f64) -> f64,
    (r, s): (f64, f64),
) -> Option<f64> {
    goursat(
        n,
        (r0, s0),
        (
            Characteristic {
                sample: phi,
                knots: &[],
            },
            Characteristic {
                sample: psi,
                knots: &[],
            },
        ),
        (r, s),
    )
}

struct Characteristic<'a> {
    sample: &'a dyn Fn(f64) -> f64,
    knots: &'a [f64],
}

impl Characteristic<'_> {
    fn weighted_integral(
        &self,
        start: f64,
        stop: f64,
        weight: impl Fn(f64) -> f64,
    ) -> Option<f64> {
        if stop <= start {
            return Some(0.0);
        }
        let integrand = |coordinate| (self.sample)(coordinate) * weight(coordinate);
        if self.knots.is_empty() {
            return Some(integrate(integrand, start, stop, 1e-12).value);
        }
        let end = self.knots.partition_point(|&knot| knot < stop);
        let interior = &self.knots[1..end];
        let quadrature = integrate_gk_breakpoints(
            integrand,
            start,
            stop,
            interior,
            GkRule::G7K15,
            1e-15,
            1e-12,
            2000.max(interior.len().saturating_mul(2)),
        )
        .ok()?;
        quadrature.converged.then_some(quadrature.value)
    }
}

fn goursat(
    n: f64,
    (r0, s0): (f64, f64),
    (phi, psi): (Characteristic<'_>, Characteristic<'_>),
    (r, s): (f64, f64),
) -> Option<f64> {
    if !([n, r0, s0, r, s].into_iter().all(f64::is_finite)
        && (r - s0).is_finite()
        && r0 > s
        && r >= r0
        && s >= s0)
    {
        return None;
    }
    let boundary = f64::midpoint((phi.sample)(r0), (psi.sample)(s0));
    let corner = epd_riemann(n, r, s, r0, s0);
    // R on the characteristics through (r, s) in closed form
    let along_rho = ((r - s0) / (r - s)).powf(n); // R(r, s0)
    let along_sigma = ((r0 - s) / (r - s)).powf(n); // R(r0, s)
    let mut total = along_rho * (phi.sample)(r) + along_sigma * (psi.sample)(s) - corner * boundary;
    let width = r - s;
    let differentiated = |rho: f64, sigma: f64, dz: f64| {
        let d = rho - sigma;
        let z = ((rho - r) / width) * ((sigma - s) / d);
        let derivative = n * (1.0 - n) * gauss(n + 1.0, 2.0 - n, 2.0, z);
        -(d / width).powf(n) * derivative * dz
    };
    // The two power-derivative terms cancel analytically in aR-R_σ and
    // bR-R_ρ. Evaluate only the derivative of the hypergeometric factor.
    if !exact_eq(n, 0.0) && !exact_eq(n, 1.0) {
        total += psi.weighted_integral(s0, s, |sigma| {
            let d = r0 - sigma;
            let dz = ((r0 - r) / width) * ((r0 - s) / d) / d;
            differentiated(r0, sigma, dz)
        })?;
        total += phi.weighted_integral(r0, r, |rho| {
            let d = rho - s0;
            let dz = ((s0 - s) / d) * ((r - s0) / width) / d;
            differentiated(rho, s0, dz)
        })?;
    }
    total.is_finite().then_some(total)
}

#[derive(Clone, Copy)]
struct Samples<'a> {
    coordinates: &'a [f64],
    values: &'a [f64],
}

impl Samples<'_> {
    fn validate(&self) -> Result<(), NumericError> {
        if self.coordinates.is_empty() || self.coordinates.len() != self.values.len() {
            return Err(NumericError::ShapeMismatch(
                "each characteristic needs matching nonempty coordinate and value arrays".into(),
            ));
        }
        if self
            .coordinates
            .iter()
            .chain(self.values)
            .any(|x| !x.is_finite())
            || self
                .coordinates
                .windows(2)
                .any(|p| p[0] >= p[1] || !(p[1] - p[0]).is_finite())
        {
            return Err(NumericError::Domain(
                "characteristic data must be finite and coordinates strictly increasing with finite spacings".into(),
            ));
        }
        Ok(())
    }

    fn at(
        &self,
        x: f64,
    ) -> f64 {
        let hi = self
            .coordinates
            .partition_point(|&point| point < x)
            .min(self.coordinates.len() - 1);
        if hi == 0 || exact_eq(self.coordinates[hi], x) {
            return self.values[hi];
        }
        let lo = hi - 1;
        let weight = (x - self.coordinates[lo]) / (self.coordinates[hi] - self.coordinates[lo]);
        self.values[lo].mul_add(1.0 - weight, self.values[hi] * weight)
    }
}

/// Evaluates a Goursat solution at ordered target pairs `[r, s]` from two
/// piecewise-linear sampled characteristics.
///
/// The first axis supplies
/// `t(ρ, s₀)` and the second supplies `t(r₀, σ)`. Their first coordinates
/// define `(r₀, s₀)`. Quadrature uses the sample knots as breakpoints.
///
/// Targets on a characteristic return the interpolated boundary data.
/// Independent interior targets can be evaluated concurrently and the output
/// order matches the input order.
///
/// # Errors
/// Returns a shape error for missing or misaligned samples. Domain errors
/// cover nonfinite or unordered data, inconsistent corner values, targets
/// outside `r ≥ r₀ > s ≥ s₀` or the sampled coordinate ranges, and quadrature
/// that cannot produce a finite, converged estimate.
pub fn epd_goursat_samples(
    n: f64,
    rho: &[f64],
    phi: &[f64],
    sigma: &[f64],
    psi: &[f64],
    targets: &[[f64; 2]],
) -> Result<Vec<f64>, NumericError> {
    let r_trace = Samples {
        coordinates: rho,
        values: phi,
    };
    let s_trace = Samples {
        coordinates: sigma,
        values: psi,
    };
    r_trace.validate()?;
    s_trace.validate()?;
    let (r0, s0) = (rho[0], sigma[0]);
    if !n.is_finite() || r0 <= s0 {
        return Err(NumericError::Domain(
            "Goursat order or corner geometry is invalid".into(),
        ));
    }
    let scale = phi[0].abs().max(psi[0].abs());
    if scale > 0.0 && (phi[0] / scale - psi[0] / scale).abs() > 16.0 * f64::EPSILON {
        return Err(NumericError::Domain(
            "characteristics disagree at their shared corner".into(),
        ));
    }
    let (r_max, s_max) = (rho[rho.len() - 1], sigma[sigma.len() - 1]);
    if targets
        .iter()
        .any(|&[r, s]| !(r0 <= r && r <= r_max && s0 <= s && s <= s_max && s < r0))
    {
        return Err(NumericError::Domain(
            "Goursat target is outside the sampled hyperbolic region".into(),
        ));
    }
    try_map_indexed(targets.len(), 16, |index| {
        let [r, s] = targets[index];
        if exact_eq(s, s0) {
            return Ok(r_trace.at(r));
        }
        if exact_eq(r, r0) {
            return Ok(s_trace.at(s));
        }
        goursat(
            n,
            (r0, s0),
            (
                Characteristic {
                    sample: &|x| r_trace.at(x),
                    knots: rho,
                },
                Characteristic {
                    sample: &|x| s_trace.at(x),
                    knots: sigma,
                },
            ),
            (r, s),
        )
        .ok_or_else(|| {
            NumericError::Domain(format!(
                "Goursat quadrature did not converge for target {index}"
            ))
        })
    })
}
