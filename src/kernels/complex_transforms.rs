//! Checked finite integral transforms and inverse Laplace quadrature.
#![allow(
    clippy::arithmetic_side_effects,
    reason = "floating-point formulas are checked for nonfinite estimates, with bounded order and panel counts"
)]

use std::f64::consts::PI;
use std::f64::consts::TAU;

use num_complex::Complex64 as C;

use crate::kernels::contour_deform::ContourIntegral;
use crate::kernels::contour_deform::ContourQuadratureError;
use crate::kernels::contour_deform::ContourSegment;
use crate::kernels::contour_deform::integrate_contour_checked;
use crate::kernels::contour_deform::integrate_segment_checked;
use crate::kernels::special_complex::bessel_j_complex;
use crate::kernels::special_complex::ln_gamma_complex;

const NAN: C = C::new(f64::NAN, f64::NAN);
const ONE: C = C::new(1.0, 0.0);
const ZERO: C = C::new(0.0, 0.0);
const DEFAULT_PANELS: usize = 400;
const FACTORIAL: [f64; 15] = [
    1.0,
    1.0,
    2.0,
    6.0,
    24.0,
    120.0,
    720.0,
    5040.0,
    40320.0,
    362880.0,
    3628800.0,
    39916800.0,
    479001600.0,
    6227020800.0,
    87178291200.0,
];
const BOREL_SCALE: f64 = 1.0 / (1_u64 << 32) as f64;

type CheckedIntegral = Result<ContourIntegral, ContourQuadratureError>;

fn on_interval(
    f: impl Fn(f64) -> C,
    a: f64,
    b: f64,
    tol: f64,
    panels: usize,
) -> CheckedIntegral {
    if !a.is_finite() || !b.is_finite() || a >= b {
        return Err(ContourQuadratureError::InvalidGeometry);
    }
    integrate_segment_checked(
        |t| f(t.re),
        ContourSegment::Line {
            start: C::new(a, 0.0),
            end: C::new(b, 0.0),
        },
        tol,
        panels,
    )
}

pub(crate) fn mellin_checked(
    f: impl Fn(f64) -> C,
    s: C,
    log_bound: f64,
    tol: f64,
    panels: usize,
) -> CheckedIntegral {
    if !(s.is_finite() && log_bound > 0.0 && log_bound.exp().is_finite()) {
        return Err(ContourQuadratureError::InvalidGeometry);
    }
    on_interval(
        |u| (s * u).exp() * f(u.exp()),
        -log_bound,
        log_bound,
        tol,
        panels,
    )
}

/// Estimates the Mellin transform on `[exp(-log_bound), exp(log_bound)]`.
pub fn mellin_transform<F: Fn(f64) -> C>(
    f: F,
    s: C,
    log_bound: f64,
    tol: f64,
) -> C {
    mellin_checked(f, s, log_bound, tol, DEFAULT_PANELS).map_or(NAN, |r| r.value)
}

pub(crate) fn inverse_mellin_checked(
    f: impl Fn(C) -> C,
    x: C,
    c: f64,
    t_max: f64,
    tol: f64,
    panels: usize,
) -> CheckedIntegral {
    if !x.is_finite() || x == ZERO || !c.is_finite() {
        return Err(ContourQuadratureError::InvalidGeometry);
    }
    let log_x = x.ln();
    on_interval(
        |t| {
            let s = C::new(c, t);
            (-s * log_x).exp() * f(s) / TAU
        },
        -t_max,
        t_max,
        tol,
        panels,
    )
}

/// Estimates inverse Mellin on the vertical line segment `c + i[-t_max, t_max]`.
pub fn inverse_mellin_transform<F: Fn(C) -> C>(
    f: F,
    x: C,
    c: f64,
    t_max: f64,
    tol: f64,
) -> C {
    inverse_mellin_checked(f, x, c, t_max, tol, DEFAULT_PANELS).map_or(NAN, |r| r.value)
}

#[derive(Clone, Copy)]
pub(crate) enum LaplaceContour {
    Talbot,
    Parabola,
}

impl LaplaceContour {
    pub(crate) const fn samples(
        self,
        n: usize,
    ) -> usize {
        match self {
            | Self::Talbot => n,
            | Self::Parabola => n / 2 * 2 + 1,
        }
    }
}

fn compensated_add<T>(
    sum: &mut T,
    correction: &mut T,
    term: T,
) where
    T: Copy + std::ops::Sub<Output = T> + std::ops::Add<Output = T>,
{
    let adjusted = term - *correction;
    let next = *sum + adjusted;
    *correction = (next - *sum) - adjusted;
    *sum = next;
}

pub(crate) fn invert_laplace(
    f: impl Fn(C) -> C,
    t: f64,
    sigma: f64,
    n: usize,
    contour: LaplaceContour,
) -> Result<C, String> {
    if !t.is_finite() || t <= 0.0 || !sigma.is_finite() || n < 4 {
        return Err("t must be finite and positive, sigma finite, and n_pts >= 4".into());
    }
    let shift = sigma * t;
    let n_f = n as f64;
    let (width, start, scale) = match contour {
        | LaplaceContour::Talbot => (TAU / n_f, -PI + PI / n_f, n_f),
        | LaplaceContour::Parabola => {
            let half = (n / 2) as f64;
            (3.0 / half, -3.0, PI * half / 12.0)
        },
    };
    let factor = C::new(0.0, -width / TAU);
    let mut sum = ZERO;
    let mut correction = ZERO;
    for k in 0..contour.samples(n) {
        if crate::cancellation::is_cancelled() {
            return Err("Laplace inversion cancelled".into());
        }
        let u = (k as f64).mul_add(width, start);
        let (position, tangent) = match contour {
            | LaplaceContour::Talbot => {
                const A: f64 = 0.6407;
                const B: f64 = 0.5017;
                let (ratio, slope) = if u.abs() < 1.0e-4 {
                    // Taylor limits of u cot(Au) and its derivative at zero.
                    let au = A * u;
                    let square = au * au;
                    (
                        (1.0 - square / 3.0 - square * square / 45.0) / A,
                        -u * A * (2.0 / 3.0 + 4.0 * square / 45.0),
                    )
                } else {
                    let (sin, cos) = (A * u).sin_cos();
                    let cot = cos / sin;
                    let sin_squared = sin * sin;
                    (u * cot, cot - (A * u) / sin_squared)
                };
                (
                    C::new(B.mul_add(ratio, -0.6122), 0.2645 * u),
                    C::new(B * slope, 0.2645),
                )
            },
            | LaplaceContour::Parabola => {
                let w = C::new(1.0, u);
                (w * w, C::new(0.0, 2.0) * w)
            },
        };
        let z = position * scale + shift;
        let s = z / t;
        if !s.is_finite() {
            return Err("Laplace contour leaves binary64 range".into());
        }
        let value = f(s);
        let term = z.exp() * (value * (scale / t)) * tangent;
        if !term.is_finite() {
            return Err("Laplace contour produced a nonfinite integrand".into());
        }
        compensated_add(&mut sum, &mut correction, term);
    }
    let value = sum * factor;
    if !value.is_finite() {
        return Err("Laplace inversion leaves binary64 range".into());
    }
    Ok(value)
}

/// Inverts Laplace along the optimized truncated Talbot contour shifted by `sigma`.
///
/// # Errors
/// Rejects invalid or nonfinite inputs, nonfinite evaluations and cancellation.
pub fn bromwich_talbot_shifted<F: Fn(C) -> C>(
    f: F,
    t: f64,
    sigma: f64,
    n_pts: usize,
) -> Result<C, String> {
    invert_laplace(f, t, sigma, n_pts, LaplaceContour::Talbot)
}

/// Inverts Laplace with the unshifted truncated Talbot contour.
///
/// # Errors
/// Rejects invalid or nonfinite inputs, nonfinite evaluations and cancellation.
pub fn bromwich_talbot<F: Fn(C) -> C>(
    f: F,
    t: f64,
    n_pts: usize,
) -> Result<C, String> {
    bromwich_talbot_shifted(f, t, 0.0, n_pts)
}

/// Inverts Laplace on Weideman's parabolic contour shifted by `sigma`.
///
/// # Errors
/// Rejects invalid or nonfinite inputs, nonfinite evaluations and cancellation.
pub fn bromwich_weideman_parabola_shifted<F: Fn(C) -> C>(
    f: F,
    t: f64,
    sigma: f64,
    n_pts: usize,
) -> Result<C, String> {
    invert_laplace(f, t, sigma, n_pts, LaplaceContour::Parabola)
}

/// Inverts Laplace on Weideman's unshifted parabolic contour.
///
/// # Errors
/// Rejects invalid or nonfinite inputs, nonfinite evaluations and cancellation.
pub fn bromwich_weideman_parabola<F: Fn(C) -> C>(
    f: F,
    t: f64,
    n_pts: usize,
) -> Result<C, String> {
    bromwich_weideman_parabola_shifted(f, t, 0.0, n_pts)
}

/// Gaver–Stehfest real-axis inverse Laplace transform with even order 2 through 14.
///
/// # Errors
/// Rejects invalid or nonfinite inputs, nonfinite evaluations and cancellation.
pub fn gaver_stehfest_inverse_laplace<F: Fn(f64) -> f64>(
    f: F,
    t: f64,
    order: usize,
) -> Result<f64, String> {
    if !t.is_finite() || t <= 0.0 || !(2..=14).contains(&order) || !order.is_multiple_of(2) {
        return Err(
            "t must be finite and positive and order must be an even integer in 2..=14".into(),
        );
    }

    let fact = |n: usize| FACTORIAL.get(n).copied().unwrap_or(f64::NAN);
    let half = order / 2;
    let ln2_t = std::f64::consts::LN_2 / t;
    let mut sum = 0.0;
    let mut correction = 0.0;
    for k in 1..=order {
        if crate::cancellation::is_cancelled() {
            return Err("Stehfest inversion cancelled".into());
        }
        let mut weight = 0.0;
        for j in k.div_ceil(2)..=k.min(half) {
            weight += (j as f64).powi(half as i32) * fact(2 * j)
                / (fact(half - j) * fact(j) * fact(j - 1) * fact(k - j) * fact(2 * j - k));
        }
        let s = (k as f64) * ln2_t;
        if !s.is_finite() {
            return Err("Stehfest sample leaves binary64 range".into());
        }
        let term = weight * f(s);
        if !term.is_finite() {
            return Err("Stehfest inversion produced a nonfinite integrand".into());
        }
        let term = if (k + half).is_multiple_of(2) {
            term
        } else {
            -term
        };
        compensated_add(&mut sum, &mut correction, term);
    }
    let value = ln2_t * sum;
    if !value.is_finite() {
        return Err("Stehfest inversion leaves binary64 range".into());
    }
    Ok(value)
}

#[derive(Clone, Copy)]
pub(crate) struct MellinBarnes<'a> {
    pub(crate) num_plus: &'a [C],
    pub(crate) num_minus: &'a [C],
    pub(crate) den_plus: &'a [C],
    pub(crate) den_minus: &'a [C],
}

fn gamma_pole(real: f64) -> bool {
    real <= 0.0 && real.fract() == 0.0
}

impl MellinBarnes<'_> {
    pub(crate) fn integrate(
        self,
        z: C,
        c: f64,
        t_max: f64,
        tol: f64,
        panels: usize,
    ) -> CheckedIntegral {
        let Self {
            num_plus,
            num_minus,
            den_plus,
            den_minus,
        } = self;
        if !z.is_finite()
            || z == ZERO
            || !c.is_finite()
            || num_plus
                .iter()
                .any(|x| !x.is_finite() || (x.im.abs() <= t_max && gamma_pole(x.re + c)))
            || num_minus
                .iter()
                .any(|x| !x.is_finite() || (x.im.abs() <= t_max && gamma_pole(x.re - c)))
            || den_plus.iter().chain(den_minus).any(|x| !x.is_finite())
        {
            return Err(ContourQuadratureError::InvalidGeometry);
        }
        let log_z = z.ln();
        on_interval(
            |t| {
                let s = C::new(c, t);
                let mut exponent = -s * log_z;
                for &a in num_plus {
                    exponent += ln_gamma_complex(a + s);
                }
                for &b in num_minus {
                    exponent += ln_gamma_complex(b - s);
                }
                for arg in den_plus
                    .iter()
                    .map(|x| x + s)
                    .chain(den_minus.iter().map(|x| x - s))
                {
                    // A denominator Gamma pole is a zero of the entire integrand.
                    if arg.im == 0.0 && gamma_pole(arg.re) {
                        return ZERO;
                    }
                    exponent -= ln_gamma_complex(arg);
                }
                exponent.exp() / TAU
            },
            -t_max,
            t_max,
            tol,
            panels,
        )
    }
}

/// Estimates a truncated Mellin–Barnes integral on a pole-free vertical segment.
#[allow(clippy::too_many_arguments)]
#[must_use]
pub fn mellin_barnes_integral(
    num_plus: &[C],
    num_minus: &[C],
    den_plus: &[C],
    den_minus: &[C],
    z: C,
    c: f64,
    t_max: f64,
    tol: f64,
) -> C {
    MellinBarnes {
        num_plus,
        num_minus,
        den_plus,
        den_minus,
    }
    .integrate(z, c, t_max, tol, DEFAULT_PANELS)
    .map_or(NAN, |estimate| estimate.value)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn meijer_checked(
    a: &[C],
    b: &[C],
    m: usize,
    n: usize,
    z: C,
    c: f64,
    t_max: f64,
    tol: f64,
    panels: usize,
) -> CheckedIntegral {
    let (b0, b1) = b
        .split_at_checked(m)
        .ok_or(ContourQuadratureError::InvalidGeometry)?;
    let (a0, a1) = a
        .split_at_checked(n)
        .ok_or(ContourQuadratureError::InvalidGeometry)?;
    let num_minus: Vec<_> = a0.iter().map(|a| ONE - a).collect();
    if b0.iter().any(|b| b.re + c <= 0.0) || num_minus.iter().any(|a| a.re - c <= 0.0) {
        return Err(ContourQuadratureError::InvalidGeometry);
    }
    let den_plus: Vec<_> = b1.iter().map(|b| ONE - b).collect();
    MellinBarnes {
        num_plus: b0,
        num_minus: &num_minus,
        den_plus: &den_plus,
        den_minus: a1,
    }
    .integrate(z, c, t_max, tol, panels)
}

/// Estimates Meijer G using its finite, pole-separating Mellin–Barnes contour.
///
/// # Errors
/// Rejects invalid indices or contours, failed quadrature and cancellation.
#[allow(clippy::too_many_arguments)]
pub fn meijer_g(
    a: &[C],
    b: &[C],
    m: usize,
    n: usize,
    z: C,
    c: f64,
    t_max: f64,
    tol: f64,
) -> Result<C, String> {
    meijer_checked(a, b, m, n, z, c, t_max, tol, DEFAULT_PANELS)
        .map(|estimate| estimate.value)
        .map_err(|error| error.to_string())
}

pub(crate) fn hilbert_checked(
    f: impl Fn(f64) -> f64,
    x: f64,
    bound: f64,
    tol: f64,
    panels: usize,
) -> CheckedIntegral {
    if !bound.is_finite() || !x.is_finite() || x.abs() >= bound {
        return Err(ContourQuadratureError::InvalidGeometry);
    }
    if crate::cancellation::is_cancelled() {
        return Err(crate::kernels::quadrature::GkError::Cancelled.into());
    }
    let x0 = x / bound;
    let fx = f(x);
    if !fx.is_finite() {
        return Err(ContourQuadratureError::InvalidGeometry);
    }
    let pv_const = (fx / PI) * (x0.ln_1p() - (-x0).ln_1p());
    if !pv_const.is_finite() {
        return Err(ContourQuadratureError::InvalidGeometry);
    }
    let split = C::new(x0, 0.0);
    let mut estimate = integrate_contour_checked(
        |u| {
            let sample = u.re * bound;
            let offset = x0 - u.re;
            let v = if offset == 0.0 {
                0.0
            } else {
                (f(sample) - fx) / offset / PI
            };
            C::new(v, 0.0)
        },
        &[
            ContourSegment::Line {
                start: -ONE,
                end: split,
            },
            ContourSegment::Line {
                start: split,
                end: ONE,
            },
        ],
        tol,
        panels,
    )?;
    estimate.value.re += pv_const;
    estimate.evaluations += 1;
    if !estimate.value.is_finite() {
        return Err(ContourQuadratureError::InvalidGeometry);
    }
    Ok(estimate)
}

/// Estimates a finite Hilbert transform by subtracting the removable singularity.
///
/// # Errors
/// Rejects invalid bounds or nonfinite estimates.
pub fn hilbert_transform_real<F: Fn(f64) -> f64>(
    f: F,
    x: f64,
    bound: f64,
    tol: f64,
) -> Result<f64, String> {
    hilbert_checked(f, x, bound, tol, DEFAULT_PANELS)
        .map(|estimate| estimate.value.re)
        .map_err(|error| error.to_string())
}

pub(crate) fn hankel_checked(
    f: impl Fn(f64) -> C,
    nu: f64,
    k: f64,
    r_max: f64,
    tol: f64,
    panels: usize,
) -> CheckedIntegral {
    if !nu.is_finite() || !k.is_finite() {
        return Err(ContourQuadratureError::InvalidGeometry);
    }
    on_interval(
        |r| f(r) * bessel_j_complex(nu, C::new(k * r, 0.0)) * r,
        0.0,
        r_max,
        tol,
        panels,
    )
}

/// Estimates the truncated Hankel transform on the finite interval `[0, r_max]`.
pub fn hankel_transform<F: Fn(f64) -> C>(
    f: F,
    nu: f64,
    k: f64,
    r_max: f64,
    tol: f64,
) -> C {
    hankel_checked(f, nu, k, r_max, tol, DEFAULT_PANELS).map_or(NAN, |estimate| estimate.value)
}

/// Evaluates the truncated Borel transform `sum(a_n t^n / n!)`.
#[must_use]
pub fn borel_transform_series(
    coeffs: &[C],
    t: C,
) -> C {
    let Some((last, coefficients)) = coeffs.split_last() else {
        return ZERO;
    };
    // In the backwards recurrence the divisions telescope to n!. The large
    // coefficients enter before division, avoiding an underflowed t^n/n!.
    let evaluate = |scale: f64| {
        coefficients
            .iter()
            .enumerate()
            .rev()
            .fold(*last * scale, |sum, (n, coeff)| {
                *coeff * scale + (t / ((n + 1) as f64)) * sum
            })
    };
    let direct = evaluate(1.0);
    if direct.is_finite() {
        return direct;
    }
    // An intermediate sum can overflow even when a subsequent multiplication
    // shrinks it. Retry with an exact power-of-two coefficient scale.
    evaluate(BOREL_SCALE) / BOREL_SCALE
}

pub(crate) fn stieltjes_checked(
    f: impl Fn(f64) -> C,
    z: C,
    r_max: f64,
    tol: f64,
    panels: usize,
) -> CheckedIntegral {
    if !z.is_finite() || (z.im == 0.0 && z.re <= 0.0) {
        return Err(ContourQuadratureError::InvalidGeometry);
    }
    on_interval(|t| f(t) / (t + z), 0.0, r_max, tol, panels)
}

/// Estimates the finite Stieltjes transform off the principal cut `(-∞, 0]`.
pub fn stieltjes_transform<F: Fn(f64) -> C>(
    f: F,
    z: C,
    r_max: f64,
    tol: f64,
) -> C {
    stieltjes_checked(f, z, r_max, tol, DEFAULT_PANELS).map_or(NAN, |estimate| estimate.value)
}
