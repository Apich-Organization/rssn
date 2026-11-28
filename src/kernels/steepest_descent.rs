//! Saddle-point asymptotics and finite complex contour integration.
#![allow(
    clippy::arithmetic_side_effects,
    reason = "the floating-point phase formulas check finite values and the work counts are fixed"
)]

use std::f64::consts::FRAC_PI_2;
use std::f64::consts::PI;

use num_complex::Complex64;

use crate::kernels::complex::cauchy_fft_derivatives;
use crate::kernels::contour_deform::ContourSegment;
use crate::kernels::contour_deform::integrate_contour_checked;
use crate::kernels::quadrature::GkRule;
use crate::kernels::quadrature::gk_panel_complex;

const ZERO: Complex64 = Complex64::new(0.0, 0.0);
const NAN: Complex64 = Complex64::new(f64::NAN, f64::NAN);
const STEPS: usize = 64;

/// Critical point of the phase `phi` in `g(z) exp(lambda * phi(z))`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SaddlePoint {
    /// Saddle point location `z0`.
    pub z0: Complex64,
    /// Phase value `phi(z0)`.
    pub phi_val: Complex64,
    /// Second derivative `phi''(z0)`.
    pub phi_d2: Complex64,
    /// Directed tangent angle in radians for local steepest descent.
    pub steepest_angle: f64,
    /// Two for a simple saddle, three when the quadratic term vanishes.
    pub degeneracy: usize,
}

fn quotient(
    numerator: Complex64,
    denominator: Complex64,
) -> Complex64 {
    let scale = denominator.re.abs().max(denominator.im.abs());
    (numerator / scale) / (denominator / scale)
}

fn unit(z: Complex64) -> Complex64 {
    let scaled = z / z.re.abs().max(z.im.abs());
    scaled / scaled.norm()
}

/// Locates a saddle point by Newton iteration, verifying the phase derivative.
/// Returns `None` on invalid values, cancellation or failure to meet the tolerance.
#[must_use]
pub fn find_saddle_point<F, D1, D2>(
    phi: F,
    dphi: D1,
    d2phi: D2,
    z_guess: Complex64,
    tol: f64,
) -> Option<SaddlePoint>
where
    F: Fn(Complex64) -> Complex64,
    D1: Fn(Complex64) -> Complex64,
    D2: Fn(Complex64) -> Complex64,
{
    if !z_guess.is_finite() || !tol.is_finite() || tol <= 0.0 {
        return None;
    }
    let mut z = z_guess;
    for _ in 0..=80 {
        if crate::cancellation::is_cancelled() {
            return None;
        }
        let slope = dphi(z);
        let curvature = d2phi(z);
        if !slope.is_finite() || !curvature.is_finite() {
            return None;
        }
        if slope.norm() <= tol {
            let phase = phi(z);
            if !phase.is_finite() {
                return None;
            }
            let nondegenerate = curvature != ZERO;
            let mut angle = if nondegenerate {
                0.5 * (PI - curvature.arg())
            } else {
                0.0
            };
            if angle > FRAC_PI_2 {
                angle -= PI;
            }
            return Some(SaddlePoint {
                z0: z,
                phi_val: phase,
                phi_d2: curvature,
                steepest_angle: angle,
                degeneracy: if nondegenerate { 2 } else { 3 },
            });
        }
        if curvature == ZERO {
            return None;
        }
        let next = z - quotient(slope, curvature);
        if !next.is_finite() || next == z {
            return None;
        }
        z = next;
    }
    None
}

fn descending_tangent(
    saddle: &SaddlePoint,
    lambda: f64,
) -> Option<Complex64> {
    if !lambda.is_finite()
        || lambda <= 0.0
        || !saddle.z0.is_finite()
        || !saddle.phi_val.is_finite()
        || !saddle.phi_d2.is_finite()
        || !saddle.steepest_angle.is_finite()
        || saddle.degeneracy != 2
    {
        return None;
    }
    let direction = Complex64::from_polar(1.0, saddle.steepest_angle);
    let decay = -unit(saddle.phi_d2) * direction * direction;
    (decay.re > 0.0).then_some(direction)
}

fn gaussian(
    saddle: &SaddlePoint,
    lambda: f64,
) -> Complex64 {
    let Some(direction) = descending_tangent(saddle, lambda) else {
        return NAN;
    };
    let phase = -unit(saddle.phi_d2) * direction * direction;
    // Take square roots before combining the two positive scales. In
    // particular, lambda * |phi''| need not be representable by itself.
    let rho = saddle.phi_d2.re.abs().max(saddle.phi_d2.im.abs());
    let width =
        (2.0 * PI).sqrt() / lambda.sqrt() / rho.sqrt() / (saddle.phi_d2 / rho).norm().sqrt();
    direction * width / phase.sqrt()
}

/// Leading oriented Gaussian term at a simple saddle. Invalid geometry gives NaN.
#[must_use]
pub fn saddle_point_leading<G>(
    g: G,
    saddle: &SaddlePoint,
    lambda: f64,
) -> Complex64
where
    G: Fn(Complex64) -> Complex64,
{
    let factor = gaussian(saddle, lambda);
    if !factor.is_finite() || crate::cancellation::is_cancelled() {
        return NAN;
    }
    g(saddle.z0) * (lambda * saddle.phi_val).exp() * factor
}

/// Oriented Gaussian term and first inverse-`lambda` correction.
/// Uses circular Cauchy samples for the third and fourth phase derivatives
/// and first and second amplitude derivatives.
#[must_use]
pub fn saddle_point_asymptotic_order2<F, G>(
    phi: F,
    g: G,
    saddle: &SaddlePoint,
    lambda: f64,
) -> Complex64
where
    F: Fn(Complex64) -> Complex64,
    G: Fn(Complex64) -> Complex64,
{
    let factor = gaussian(saddle, lambda);
    if !factor.is_finite() || crate::cancellation::is_cancelled() {
        return NAN;
    }
    let radius = 0.1;
    if saddle.z0.re + radius <= saddle.z0.re || saddle.z0.im + radius <= saddle.z0.im {
        return NAN;
    }
    let p = cauchy_fft_derivatives(phi, saddle.z0, radius, 4, 32);
    let g_derivatives = cauchy_fft_derivatives(&g, saddle.z0, radius, 2, 32);
    let ([_, _, _, p3, p4], [_, g1, g2]) = (p.as_slice(), g_derivatives.as_slice()) else {
        return NAN;
    };
    let g0 = g(saddle.z0);
    // DLMF 2.4.16 under p = -phi. Keep the correction linear in g
    // and its derivatives, including the mixed term when g(z0) vanishes.
    let c3 = quotient(*p3, saddle.phi_d2);
    let c4 = quotient(*p4, saddle.phi_d2);
    let numerator = -0.5 * *g2 + 0.5 * *g1 * c3 + g0 * (0.125 * c4 - (5.0 / 24.0) * c3 * c3);
    let correction = quotient(numerator, saddle.phi_d2) / lambda;
    (g0 + correction) * (lambda * saddle.phi_val).exp() * factor
}

#[allow(
    clippy::single_call_fn,
    reason = "construct the two oriented phase branches independently of quadrature"
)]
fn trace_descent(
    phi: &impl Fn(Complex64) -> Complex64,
    saddle: &SaddlePoint,
    direction: Complex64,
    half_len: f64,
) -> Option<Vec<Complex64>> {
    let ds = half_len / STEPS as f64;
    let mut points = Vec::with_capacity(STEPS + 1);
    let mut point = saddle.z0;
    let mut tangent = direction;
    let mut last_real = saddle.phi_val.re;
    let phase_tol = 1e-12 * saddle.phi_val.norm().max(1.0);
    points.push(point);
    for _ in 0..STEPS {
        if crate::cancellation::is_cancelled() {
            return None;
        }
        let mut next = point + tangent * ds;
        if !next.is_finite() || next == point {
            return None;
        }
        let mut accepted = false;
        for _ in 0..8 {
            let phase = phi(next);
            if !phase.is_finite() {
                return None;
            }
            let h = f64::EPSILON.cbrt() * (1.0 + next.re.abs());
            let slope = (phi(next + h) - phi(next - h)) / (2.0 * h);
            if slope == ZERO || !slope.is_finite() {
                return None;
            }
            let residual = phase.im - saddle.phi_val.im;
            if residual.abs() <= phase_tol {
                if phase.re > last_real + phase_tol {
                    return None;
                }
                last_real = phase.re;
                tangent = -unit(slope).conj();
                accepted = true;
                break;
            }
            next -= quotient(Complex64::new(0.0, residual), slope);
            if !next.is_finite() {
                return None;
            }
        }
        if !accepted {
            return None;
        }
        point = next;
        points.push(point);
    }
    Some(points)
}

/// Integrates along two traced constant-phase branches of a simple saddle.
///
/// The finite path extends approximately `half_len` in each direction.
/// Returns NaN on invalid geometry, branch projection failure or nonconvergence.
#[must_use]
pub fn steepest_descent_integrate<F, G>(
    phi: F,
    g: G,
    saddle: &SaddlePoint,
    lambda: f64,
    half_len: f64,
    tol: f64,
) -> Complex64
where
    F: Fn(Complex64) -> Complex64,
    G: Fn(Complex64) -> Complex64,
{
    let Some(direction) = descending_tangent(saddle, lambda) else {
        return NAN;
    };
    if !half_len.is_finite()
        || half_len <= 0.0
        || !tol.is_finite()
        || tol <= 0.0
        || crate::cancellation::is_cancelled()
    {
        return NAN;
    }
    let actual = phi(saddle.z0);
    if !actual.is_finite()
        || (actual - saddle.phi_val).norm() > 1e-12 * saddle.phi_val.norm().max(1.0)
    {
        return NAN;
    }
    let Some(negative) = trace_descent(&phi, saddle, -direction, half_len) else {
        return NAN;
    };
    let Some(positive) = trace_descent(&phi, saddle, direction, half_len) else {
        return NAN;
    };
    let segments: Vec<_> = negative
        .windows(2)
        .rev()
        .filter_map(|pair| {
            let [end, start] = *pair else { return None };
            Some(ContourSegment::Line { start, end })
        })
        .chain(positive.windows(2).filter_map(|pair| {
            let [start, end] = *pair else { return None };
            Some(ContourSegment::Line { start, end })
        }))
        .collect();
    let evaluate = |z| g(z) * (lambda * phi(z)).exp();
    // For a resolved trace, reuse the shared complex panel's samples and
    // error estimates directly. Refine through the adaptive contour engine
    // only when this partition misses the requested absolute tolerance.
    let mut integral = ZERO;
    let mut rounding = ZERO;
    let mut error = 0.0;
    for segment in &segments {
        if crate::cancellation::is_cancelled() {
            return NAN;
        }
        let &ContourSegment::Line { start, end } = segment else {
            return NAN;
        };
        let delta = end - start;
        let (estimate, errors, _) = gk_panel_complex(0.0, 1.0, GkRule::G7K15, |t| {
            evaluate(start + t * delta) * delta
        });
        if !estimate.is_finite() || !errors[0].is_finite() || !errors[1].is_finite() {
            return NAN;
        }
        let adjusted = estimate - rounding;
        let next = integral + adjusted;
        rounding = (next - integral) - adjusted;
        integral = next;
        error += errors[0].hypot(errors[1]);
    }
    if error <= tol && integral.is_finite() {
        return integral;
    }
    integrate_contour_checked(evaluate, &segments, tol, 4096).map_or(NAN, |result| {
        if result.converged {
            result.value
        } else {
            NAN
        }
    })
}
