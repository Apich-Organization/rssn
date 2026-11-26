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

    // Standard Laplace/saddle-point O(1/lambda) correction, or O(lambda^(-3/2)) when g(z0) == 0:
    if g0.norm() <= 1e-14 {
        let second_term = 0.5 * g2 / (-lambda * p2);
        return (lambda * saddle.phi_val).exp() * gauss * second_term;
    }
    let corr = (-g2 / (2.0 * g0 * p2) + (g1 * p3) / (2.0 * g0 * p2 * p2) + p4 / (8.0 * p2 * p2)
        - (5.0 * p3 * p3) / (24.0 * p2 * p2 * p2))
        / lambda;
    g0 * (lambda * saddle.phi_val).exp() * gauss * (Complex64::new(1.0, 0.0) + corr)
}

/// Numerical integration along the **curved** steepest-descent path
/// `Im(phi(z)) = Im(phi(z0))` with `Re(phi(z))` decreasing away from `z0`,
/// traced via predictor-corrector projection onto the constant-phase level curve.
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
    let im0 = saddle.phi_val.im;
    let trace_branch = |sign: f64| -> Complex64 {
        let n_steps = 64_usize;
        let ds = half_len / (n_steps as f64);
        let mut pts = Vec::with_capacity(n_steps + 1);
        pts.push(saddle.z0);
        let mut z = saddle.z0;
        let mut tangent = dir * sign;
        for _ in 0..n_steps {
            z += tangent * ds;
            // Project onto Im(phi(z)) = im0 via 3 Newton steps perpendicular to grad(Im phi).
            for _ in 0..4 {
                let h = 1e-6 * (1.0 + z.norm());
                let dp = (phi(z + h) - phi(z - h)) / (2.0 * h);
                if dp.norm() < 1e-14 {
                    break;
                }
                let err_im = phi(z).im - im0;
                if err_im.abs() < 1e-13 {
                    break;
                }
                // d/dt Im(phi(z + i * t / dp)) = Re(i * t) -> use z - i * err_im / dp
                z -= Complex64::new(0.0, err_im) / dp;
                tangent = -dp.conj() / dp.norm();
                if (tangent.re * (dir.re * sign) + tangent.im * (dir.im * sign)) < 0.0 {
                    tangent = -tangent;
                }
            }
            pts.push(z);
        }
        // Composite trapezoidal/Simpson on the traced steepest-descent curve
        let mut seg_sum = Complex64::new(0.0, 0.0);
        for w in pts.windows(2) {
            let a = w[0];
            let b = w[1];
            let re = gauss_kronrod_any(
                |t| {
                    let zt = a + (b - a) * t;
                    (g(zt) * (lambda * phi(zt)).exp() * (b - a)).re
                },
                0.0,
                1.0,
                tol,
                16,
            )
            .value;
            let im = gauss_kronrod_any(
                |t| {
                    let zt = a + (b - a) * t;
                    (g(zt) * (lambda * phi(zt)).exp() * (b - a)).im
                },
                0.0,
                1.0,
                tol,
                16,
            )
            .value;
            seg_sum += Complex64::new(re, im);
        }
        seg_sum
    };
    trace_branch(1.0) - trace_branch(-1.0)
}
