//! Steepest-descent and saddle-point asymptotic methods for oscillatory and
//! exponential contour integrals `I(lambda) = int_C g(z) exp(lambda * phi(z)) dz`.

use std::f64::consts::FRAC_PI_2;
use std::f64::consts::PI;

use num_complex::Complex64;

use crate::kernels::complex::cauchy_fft_derivatives;
use crate::kernels::integrate::gauss_kronrod_any;

/// Non-degenerate or higher-order saddle point `phi'(z0) = 0`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SaddlePoint {
    /// Saddle point location `z0`.
    pub z0: Complex64,
    /// Phase value `phi(z0)`.
    pub phi_val: Complex64,
    /// Second derivative `phi''(z0)`.
    pub phi_d2: Complex64,
    /// Steepest-descent direction angle `theta` where `Re(e^(2 i theta) phi''(z0)) < 0`.
    pub steepest_angle: f64,
    /// Degeneracy order (`2` for a simple non-degenerate saddle).
    pub degeneracy: usize,
}

fn orient_steepest_angle(f2: Complex64) -> f64 {
    if f2.norm() <= 1e-12 {
        return 0.0;
    }
    let mut theta = 0.5 * (PI - f2.arg());
    while theta > FRAC_PI_2 {
        theta -= PI;
    }
    while theta <= -FRAC_PI_2 {
        theta += PI;
    }
    theta
}

/// Locates a saddle point `phi'(z0) = 0` from `z_guess` via Newton's method.
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
    let mut z = z_guess;
    for _ in 0..80 {
        let f1 = dphi(z);
        if f1.norm() <= tol {
            let f2 = d2phi(z);
            return Some(SaddlePoint {
                z0: z,
                phi_val: phi(z),
                phi_d2: f2,
                steepest_angle: orient_steepest_angle(f2),
                degeneracy: if f2.norm() > 1e-10 { 2 } else { 3 },
            });
        }
        let f2 = d2phi(z);
        if f2.norm() < 1e-30 {
            return None;
        }
        let step = f1 / f2;
        z -= step;
        if step.norm() <= tol * (1.0 + z.norm()) {
            let f2_fin = d2phi(z);
            return Some(SaddlePoint {
                z0: z,
                phi_val: phi(z),
                phi_d2: f2_fin,
                steepest_angle: orient_steepest_angle(f2_fin),
                degeneracy: 2,
            });
        }
    }
    None
}

/// Leading-order saddle-point approximation:
/// `g(z0) * exp(lambda * phi(z0)) * sqrt(2*pi / (-lambda * phi''(z0)))`
/// with branch oriented along `saddle.steepest_angle`.
#[must_use]
pub fn saddle_point_leading<G>(
    g: G,
    saddle: &SaddlePoint,
    lambda: f64,
) -> Complex64
where
    G: Fn(Complex64) -> Complex64,
{
    let dir = Complex64::from_polar(1.0, saddle.steepest_angle);
    let a = -0.5 * lambda * saddle.phi_d2 * dir * dir; // Re(a) > 0
    let gauss = dir * (PI / a).sqrt();
    g(saddle.z0) * (lambda * saddle.phi_val).exp() * gauss
}

/// Second-order saddle-point approximation including the `O(1/lambda)` term,
/// using Cauchy-FFT derivatives of `phi` (to order 4) and `g` (to order 2) at `z0`.
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
    let r = 0.1_f64.min(0.25 * (1.0 + saddle.z0.norm()));
    let p_der = cauchy_fft_derivatives(phi, saddle.z0, r, 4, 32);
    let g_der = cauchy_fft_derivatives(g, saddle.z0, r, 2, 32);
    let g0 = g_der.first().copied().unwrap_or_default();
    let g1 = g_der.get(1).copied().unwrap_or_default();
    let g2 = g_der.get(2).copied().unwrap_or_default();
    let p2 = p_der.get(2).copied().unwrap_or(saddle.phi_d2);
    let p3 = p_der.get(3).copied().unwrap_or_default();
    let p4 = p_der.get(4).copied().unwrap_or_default();

    let dir = Complex64::from_polar(1.0, saddle.steepest_angle);
    let a = -0.5 * lambda * p2 * dir * dir;
    let gauss = dir * (PI / a).sqrt();

    // Standard Laplace/saddle-point O(1/lambda) correction:
    // 1 + (1/lambda) * [ -g2/(2*g0*p2) + g1*p3/(2*g0*p2^2) + p4/(8*p2^2) - 5*p3^2/(24*p2^3) ]
    let corr = if g0.norm() > 1e-14 {
        (-g2 / (2.0 * g0 * p2) + (g1 * p3) / (2.0 * g0 * p2 * p2) + p4 / (8.0 * p2 * p2)
            - (5.0 * p3 * p3) / (24.0 * p2 * p2 * p2))
            / lambda
    } else {
        Complex64::new(0.0, 0.0)
    };
    g0 * (lambda * saddle.phi_val).exp() * gauss * (Complex64::new(1.0, 0.0) + corr)
}

/// Numerical integration along the steepest-descent tangent ray
/// `z(s) = z0 + s * exp(i * steepest_angle)` for `s in [-half_len, half_len]`.
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
    let dir = Complex64::from_polar(1.0, saddle.steepest_angle);
    let integrand = |s: f64| -> Complex64 {
        let z = saddle.z0 + dir * s;
        g(z) * (lambda * phi(z)).exp() * dir
    };
    let re = gauss_kronrod_any(|s| integrand(s).re, -half_len, half_len, tol, 24).value;
    let im = gauss_kronrod_any(|s| integrand(s).im, -half_len, half_len, tol, 24).value;
    Complex64::new(re, im)
}
