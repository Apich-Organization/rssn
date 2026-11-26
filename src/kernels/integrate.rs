//! # Numerical Integration (Quadrature)
//!
//! This module provides various numerical integration (quadrature) methods for approximating
//! definite integrals of functions. It includes implementations of the Trapezoidal rule,
//! Simpson's rule, Adaptive Quadrature (Adaptive Simpson's), Romberg Integration, and
//! Gauss-Legendre Quadrature.
//!
//! # Methods
//!
//! * **Trapezoidal Rule**: Simple and robust, linear approximation.
//! * **Simpson's Rule**: Quadratic approximation, generally more accurate than Trapezoidal for smooth functions.
//! * **Adaptive Quadrature**: Recursively subdivides intervals to achieve a specified error tolerance.
//! * **Romberg Integration**: Uses Richardson extrapolation on Trapezoidal approximation to achieve high order accuracy.
//! * **Gauss-Legendre Quadrature**: Uses optimal sample points (roots of Legendre polynomials) for high accuracy with fewer function evaluations, ideal for smooth functions.

use serde::Deserialize;
use serde::Serialize;

/// Enum to select the numerical integration method.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum QuadratureMethod {
    /// Trapezoidal rule.
    Trapezoidal,
    /// Simpson's 1/3 rule.
    Simpson,
    /// Adaptive Simpson's quadrature.
    Adaptive,
    /// Romberg integration.
    Romberg,
    /// Gauss-Legendre quadrature (n=5).
    GaussLegendre,
}

/// # Pure Numerical Trapezoidal Rule
///
/// Performs numerical integration using the trapezoidal rule for a given function.
///
/// ## Arguments
/// * `f` - The function to integrate, as a closure.
/// * `range` - A tuple `(a, b)` representing the integration interval.
/// * `n_steps` - The number of steps to use.
///
/// ## Returns
/// The approximate value of the definite integral.
///
/// ## Example
/// ```
/// use rssn::kernels::integrate::trapezoidal_rule;
///
/// let f = |x: f64| x * x;
///
/// let res = trapezoidal_rule(f, (0.0, 1.0), 1000);
///
/// assert!((res - 1.0 / 3.0).abs() < 1e-4);
/// ```
pub fn trapezoidal_rule<F>(
    f: F,
    range: (f64, f64),
    n_steps: usize,
) -> f64
where
    F: Fn(f64) -> f64,
{
    let (a, b) = range;

    if !a.is_finite() || !b.is_finite() {
        return f64::NAN;
    }

    if n_steps == 0 {
        return 0.0;
    }

    // Correctly handle a == b (integral is 0)
    // and a > b (integral is negative of b to a)
    // The logic below works for a > b as h will be negative.
    if crate::kernels::exact_eq(a, b) {
        return 0.0;
    }

    let h = (b - a) / (n_steps as f64);

    let mut sum = f64::midpoint(f(a), f(b));

    for i in 1..n_steps {
        let x = (i as f64).mul_add(h, a);

        sum += f(x);
    }

    h * sum
}

/// # Pure Numerical Simpson's Rule
///
/// Performs numerical integration using Simpson's 1/3 rule.
///
/// ## Arguments
/// * `f` - The function to integrate.
/// * `range` - Interval `(a, b)`.
/// * `n_steps` - Number of steps (must be usually even, but we handle it).
///
/// ## Returns
/// `Result<f64, String>`
///
/// ## Errors
/// Returns an error if the number of steps is zero or if an internal calculation fails.
///
/// ## Example
/// ```
/// use rssn::kernels::integrate::simpson_rule;
///
/// let f = |x: f64| x * x;
///
/// let res = simpson_rule(f, (0.0, 1.0), 10).unwrap();
///
/// assert!((res - 1.0 / 3.0).abs() < 1e-10);
/// ```
pub fn simpson_rule<F>(
    f: F,
    range: (f64, f64),
    n_steps: usize,
) -> Result<f64, String>
where
    F: Fn(f64) -> f64,
{
    let (a, b) = range;

    if !a.is_finite() || !b.is_finite() {
        return Err("simpson_rule needs finite bounds".to_string());
    }

    if n_steps == 0 {
        return Ok(0.0);
    }

    if crate::kernels::exact_eq(a, b) {
        return Ok(0.0);
    }

    // Simpson's rule requires even number of intervals for the strict global formula.
    // If odd, we can warn or adjust. For now, enforce even.
    let steps = if n_steps.is_multiple_of(2) {
        n_steps
    } else {
        n_steps + 1
    };

    let h = (b - a) / (steps as f64);

    let mut sum = f(a) + f(b);

    for i in 1..steps {
        crate::cancellation::checkpoint().map_err(|error| error.to_string())?;
        let x = (i as f64).mul_add(h, a);

        let weight = if i % 2 == 0 { 2.0 } else { 4.0 };

        sum += weight * f(x);
    }

    Ok((h / 3.0) * sum)
}

// A panel's numerical estimate is independent of the refinement budget.
// Keep the ordinary arithmetic order, and rescale the operands when an
// intermediate width or sum would overflow or the divided width is subnormal.
fn simpson_panel_estimate(
    a: f64,
    b: f64,
    fa: f64,
    fm: f64,
    fb: f64,
) -> f64 {
    let width = b - a;
    let sixth = width / 6.0;
    let sum = 4.0_f64.mul_add(fm, fa) + fb;
    if sixth.is_normal() && sum.is_finite() {
        return sixth * sum;
    }
    scaled_simpson_estimate(a, b, width, [fa, fm, fb], sum)
}

#[cold]
fn scaled_simpson_estimate(
    a: f64,
    b: f64,
    width: f64,
    values: [f64; 3],
    sum: f64,
) -> f64 {
    if !sum.is_finite() {
        let [fa, fm, fb] = values;
        if !fa.is_finite() || !fm.is_finite() || !fb.is_finite() {
            return width * sum;
        }
        // The scaled sum is bounded by 3/4 of f64::MAX. Leave its inverse
        // factor outside the product, so the weighted mean need not overflow.
        let scaled = 4.0_f64.mul_add(0.125 * fm, 0.125 * fa) + 0.125 * fb;
        return if width.is_finite() {
            (width * (scaled / 3.0)) * 4.0
        } else {
            ((0.5 * b - 0.5 * a) * (scaled / 3.0)) * 8.0
        };
    }
    if width.is_finite() {
        width * (sum / 6.0)
    } else {
        ((0.5 * b - 0.5 * a) / 3.0) * sum
    }
}

// The individual f64 parameters let a recursive call pass the inherited
// samples and panel data directly. Each split evaluates only its new nodes.
#[allow(clippy::too_many_arguments)]
fn refine_simpson(
    f: &impl Fn(f64) -> f64,
    a: f64,
    mid: f64,
    b: f64,
    fa: f64,
    fm: f64,
    fb: f64,
    whole: f64,
    tolerance: f64,
    limit: usize,
    budget: &mut usize,
) -> f64 {
    if limit == 0 || *budget == 0 || !whole.is_finite() {
        return whole;
    }
    *budget -= 1;

    let left = f64::midpoint(a, mid);
    let right = f64::midpoint(mid, b);
    let exact_eq = crate::kernels::exact_eq;
    if exact_eq(a, left) || exact_eq(left, mid) || exact_eq(mid, right) || exact_eq(right, b) {
        return refine_simpson_at_resolution(
            f,
            [a, mid, b],
            [fa, fm, fb],
            whole,
            tolerance / 2.0,
            limit - 1,
            budget,
        );
    }

    let fl = f(left);
    let fr = f(right);
    let left_width = mid - a;
    let right_width = b - mid;
    let left_sum = 4.0_f64.mul_add(fl, fa) + fm;
    let right_sum = 4.0_f64.mul_add(fr, fm) + fb;
    let left_estimate = (left_width / 6.0) * left_sum;
    let right_estimate = (right_width / 6.0) * right_sum;
    let (left_estimate, right_estimate) = if (left_estimate + right_estimate).is_finite()
        && left_width.abs().min(right_width.abs()) >= 6.0 * f64::MIN_POSITIVE
    {
        (left_estimate, right_estimate)
    } else {
        (
            scaled_simpson_estimate(a, mid, left_width, [fa, fl, fm], left_sum),
            scaled_simpson_estimate(mid, b, right_width, [fm, fr, fb], right_sum),
        )
    };
    let sum = left_estimate + right_estimate;
    if !sum.is_finite() {
        return sum;
    }
    let correction = (sum - whole) / 15.0;
    if correction.abs() <= tolerance {
        return sum + correction;
    }

    let next_tolerance = tolerance / 2.0;
    refine_simpson(
        f,
        a,
        left,
        mid,
        fa,
        fl,
        fm,
        left_estimate,
        next_tolerance,
        limit - 1,
        budget,
    ) + refine_simpson(
        f,
        mid,
        right,
        b,
        fm,
        fr,
        fb,
        right_estimate,
        next_tolerance,
        limit - 1,
        budget,
    )
}

#[cold]
fn refine_simpson_at_resolution(
    f: &impl Fn(f64) -> f64,
    bounds: [f64; 3],
    values: [f64; 3],
    whole: f64,
    tolerance: f64,
    limit: usize,
    budget: &mut usize,
) -> f64 {
    let [a, mid, b] = bounds;
    let [fa, fm, fb] = values;
    let left = f64::midpoint(a, mid);
    let right = f64::midpoint(mid, b);
    let exact_eq = crate::kernels::exact_eq;
    let left_collapsed = exact_eq(a, left) || exact_eq(left, mid);
    let right_collapsed = exact_eq(mid, right) || exact_eq(right, b);
    if left_collapsed && right_collapsed {
        let left_width = (mid - a).abs();
        let right_width = (b - mid).abs();
        if exact_eq(left_width, right_width) {
            return whole;
        }
        // Adjacent spacings differ only at a binade boundary, where one
        // interval is twice the other. Integrating the interpolating quadratic
        // gives 3/4 of the middle value and 1/4 of the farther endpoint.
        let farther = if left_width > right_width {
            fa
        } else {
            fb
        };
        return (b - a) * f64::midpoint(fm, f64::midpoint(fm, farther));
    }
    // An adjacent pair uses the trapezoid. Richardson's estimate assumes
    // Simpson's rule in both children and does not apply in this case.
    if left_collapsed {
        let fr = f(right);
        let right_estimate = simpson_panel_estimate(mid, b, fm, fr, fb);
        (mid - a) * f64::midpoint(fa, fm)
            + refine_simpson(
                f,
                mid,
                right,
                b,
                fm,
                fr,
                fb,
                right_estimate,
                tolerance,
                limit,
                budget,
            )
    } else {
        let fl = f(left);
        let left_estimate = simpson_panel_estimate(a, mid, fa, fl, fm);
        refine_simpson(
            f,
            a,
            left,
            mid,
            fa,
            fl,
            fm,
            left_estimate,
            tolerance,
            limit,
            budget,
        ) + (b - mid) * f64::midpoint(fm, fb)
    }
}

/// # Adaptive Quadrature (Adaptive Simpson's Method)
///
/// Recursively refines the interval until the estimated error satisfies the tolerance
/// or the refinement limit is reached. The integrand must return a consistent result
/// for a given coordinate because refinement reuses the parent panel's samples.
/// A panel with consecutive floating-point bounds uses a trapezoid through its
/// two finite endpoint samples, or returns `NaN` if either sample is nonfinite.
///
/// ## Arguments
/// * `f` - Function to integrate.
/// * `range` - Interval `(a, b)`.
/// * `tolerance` - Error tolerance (e.g., 1e-6).
///
/// ## Example
/// ```
/// use rssn::kernels::integrate::adaptive_quadrature;
///
/// let f = |x: f64| x.sin();
///
/// let res = adaptive_quadrature(f, (0.0, std::f64::consts::PI), 1e-6);
///
/// assert!((res - 2.0).abs() < 1e-6);
/// ```
pub fn adaptive_quadrature<F>(
    f: F,
    range: (f64, f64),
    tolerance: f64,
) -> f64
where
    F: Fn(f64) -> f64,
{
    let (a, b) = range;

    if !a.is_finite() || !b.is_finite() || tolerance.is_nan() {
        return f64::NAN;
    }

    if crate::kernels::exact_eq(a, b) {
        return 0.0;
    }

    let mid = f64::midpoint(a, b);
    if crate::kernels::exact_eq(a, mid) || crate::kernels::exact_eq(mid, b) {
        let (fa, fb) = (f(a), f(b));
        return if fa.is_finite() && fb.is_finite() {
            (b - a) * f64::midpoint(fa, fb)
        } else {
            f64::NAN
        };
    }
    let fm = f(mid);
    let (fa, fb) = (f(a), f(b));

    if fa.is_nan() || fb.is_nan() || fm.is_nan() {
        return f64::NAN;
    }

    if !fa.is_finite() || !fb.is_finite() {
        // An integrable endpoint singularity such as 1/sqrt(x) on [0, 1]:
        // Simpson would evaluate the pole itself. The double-exponential rule
        // never samples the endpoints and converges for algebraic and
        // logarithmic singularities.
        return crate::kernels::quadrature::tanh_sinh(&f, a, b, tolerance.max(1e-14)).value;
    }

    let initial = simpson_panel_estimate(a, b, fa, fm, fb);

    // The tolerance floor and panel budget bound work even for zero tolerance
    // or a non-integrable singularity.
    let tolerance = tolerance.max(1e-15 * initial.abs());

    refine_simpson(
        &f,
        a,
        mid,
        b,
        fa,
        fm,
        fb,
        initial,
        tolerance,
        60,
        &mut 50_000,
    )
}

/// # Romberg Integration
///
/// Uses Richardson extrapolation on the Trapezoidal rule to improve accuracy.
///
/// ## Arguments
/// * `f` - Function to integrate.
/// * `range` - Interval.
/// * `max_steps` - Order of extrapolation (e.g., 5-10).
///
/// ## Example
/// ```
/// use rssn::kernels::integrate::romberg_integration;
///
/// let f = |x: f64| x.exp();
///
/// let res = romberg_integration(f, (0.0, 1.0), 6);
///
/// assert!((res - (std::f64::consts::E - 1.0)).abs() < 1e-10);
/// ```
pub fn romberg_integration<F>(
    f: F,
    range: (f64, f64),
    max_steps: usize,
) -> f64
where
    F: Fn(f64) -> f64,
{
    let (a, b) = range;

    if !a.is_finite() || !b.is_finite() {
        return f64::NAN;
    }

    if max_steps == 0 {
        return 0.0;
    }

    if crate::kernels::exact_eq(a, b) {
        return 0.0;
    }

    // 2^21 trapezoid panels is far beyond double precision; the cap also
    // keeps the shifts below and the table size (max_steps^2) sane.
    let max_steps = max_steps.min(22);

    let mut r = vec![vec![0.0; max_steps]; max_steps];

    // R[0][0]
    let h = b - a;

    r[0][0] = 0.5 * h * (f(a) + f(b));

    for i in 1..max_steps {
        // Calculate R[i][0] using Trapezoidal rule with 2^i segments
        // But we can update from R[i-1][0] efficiently
        let steps_prev = 1_u64 << (i - 1);

        let h_i = h / (1_u64 << i) as f64;

        let mut sum = 0.0;

        for k in 1..=steps_prev {
            let x = a + (2 * k - 1) as f64 * h_i;

            sum += f(x);
        }

        r[i][0] = 0.5f64.mul_add(r[i - 1][0], h_i * sum);

        // Richardson extrapolation
        for j in 1..=i {
            let k = 4.0_f64.powi(j as i32);

            r[i][j] = k.mul_add(r[i][j - 1], -r[i - 1][j - 1]) / (k - 1.0);
        }
    }

    r[max_steps - 1][max_steps - 1]
}

/// # Gauss-Legendre Quadrature
///
/// Uses standard weights and nodes for n = 5 (hardcoded for now as it's efficient for general purpose).
/// Can be extended to arbitrary n in the future.
///
/// ## Example
/// ```
/// use rssn::kernels::integrate::gauss_legendre_quadrature;
///
/// let f = |x: f64| x.powi(3);
///
/// let res = gauss_legendre_quadrature(f, (0.0, 1.0));
///
/// assert!((res - 0.25).abs() < 1e-10);
/// ```
pub fn gauss_legendre_quadrature<F>(
    f: F,
    range: (f64, f64),
) -> f64
where
    F: Fn(f64) -> f64,
{
    let (a, b) = range;

    if !a.is_finite() || !b.is_finite() {
        return f64::NAN;
    }

    if crate::kernels::exact_eq(a, b) {
        return 0.0;
    }

    crate::kernels::quadrature::gauss_legendre(f, a, b, 5)
}

/// Result of an adaptive quadrature: the estimate and an error bound
/// estimate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Quadrature {
    /// Estimate of the integral.
    pub value: f64,
    /// Estimated absolute error.
    pub error: f64,
    /// Number of integrand evaluations used.
    pub evaluations: usize,
}

// Gauss–Kronrod 7/15 nodes and weights on [-1, 1] (positive half).
const GK_NODES: [f64; 8] = [
    0.991_455_371_120_812_6,
    0.949_107_912_342_758_5,
    0.864_864_423_359_769_1,
    0.741_531_185_599_394_4,
    0.586_087_235_467_691_1,
    0.405_845_151_377_397_2,
    0.207_784_955_007_898_47,
    0.0,
];
const KRONROD_WEIGHTS: [f64; 8] = [
    0.022_935_322_010_529_225,
    0.063_092_092_629_978_55,
    0.104_790_010_322_250_18,
    0.140_653_259_715_525_92,
    0.169_004_726_639_267_9,
    0.190_350_578_064_785_4,
    0.204_432_940_075_298_9,
    0.209_482_141_084_727_82,
];
// Gauss weights for the odd-indexed nodes (1, 3, 5, 7).
const GAUSS_WEIGHTS: [f64; 4] = [
    0.129_484_966_168_869_7,
    0.279_705_391_489_276_64,
    0.381_830_050_505_118_9,
    0.417_959_183_673_469_4,
];

/// One Gauss–Kronrod 7/15 panel on `[a, b]`: `(estimate, error estimate)`.
fn gk15(
    f: &impl Fn(f64) -> f64,
    a: f64,
    b: f64,
) -> (f64, f64) {
    let centre = f64::midpoint(a, b);
    let half = 0.5 * b - 0.5 * a;
    let mut kronrod = 0.0;
    let mut gauss = 0.0;
    for (i, (&node, &weight)) in GK_NODES.iter().zip(&KRONROD_WEIGHTS).enumerate() {
        let sum = if node == 0.0 {
            f(centre)
        } else {
            f(centre - half * node) + f(centre + half * node)
        };
        kronrod += weight * sum;
        if i % 2 == 1 {
            gauss += GAUSS_WEIGHTS[i / 2] * sum;
        }
    }
    (kronrod * half, ((kronrod - gauss) * half).abs())
}

/// Globally adaptive Gauss–Kronrod quadrature of `f` over the finite
/// interval `[a, b]`.
///
/// The panel with the largest error estimate is bisected until the total
/// estimate drops below `tolerance` or `max_panels` panels exist. The
/// returned [`Quadrature::error`] is the sum of the panel estimates, so a
/// caller can tell a converged result from an exhausted budget.
/// Reversing the limits negates the integral. Non-finite integrand values
/// are propagated to the result rather than treated as zero.
#[allow(clippy::float_cmp)] // exact comparison against a sentinel / integer-valued input is intended
pub fn gauss_kronrod(
    f: impl Fn(f64) -> f64,
    a: f64,
    b: f64,
    tolerance: f64,
    max_panels: usize,
) -> Quadrature {
    if a > b {
        let flipped = gauss_kronrod(f, b, a, tolerance, max_panels);
        return Quadrature {
            value: -flipped.value,
            ..flipped
        };
    }
    if a == b {
        return Quadrature {
            value: 0.0,
            error: 0.0,
            evaluations: 0,
        };
    }
    let (value, error) = gk15(&f, a, b);
    let mut panels = vec![(a, b, value, error)];
    let mut evaluations = 15;
    loop {
        let total_error: f64 = panels.iter().map(|p| p.3).sum();
        let total: f64 = panels.iter().map(|p| p.2).sum();
        if total_error <= tolerance.max(f64::EPSILON * total.abs())
            || panels.len() >= max_panels
            || !total.is_finite()
        {
            return Quadrature {
                value: total,
                error: total_error,
                evaluations,
            };
        }
        let worst = panels
            .iter()
            .enumerate()
            .max_by(|x, y| x.1.3.total_cmp(&y.1.3))
            .map_or(0, |(i, _)| i);
        let (lo, hi, value, error) = panels.swap_remove(worst);
        let mid = f64::midpoint(lo, hi);
        if mid <= lo || mid >= hi {
            // The panel cannot be split any further in floating point.
            panels.push((lo, hi, value, error));
            let total: f64 = panels.iter().map(|p| p.2).sum();
            return Quadrature {
                value: total,
                error: f64::INFINITY,
                evaluations,
            };
        }
        for (x0, x1) in [(lo, mid), (mid, hi)] {
            let (v, e) = gk15(&f, x0, x1);
            panels.push((x0, x1, v, e));
        }
        evaluations += 30;
    }
}

/// Adaptive Gauss–Kronrod quadrature over an interval that may be
/// infinite at either end, by mapping it to a finite one.
pub fn gauss_kronrod_any(
    f: impl Fn(f64) -> f64,
    a: f64,
    b: f64,
    tolerance: f64,
    max_panels: usize,
) -> Quadrature {
    if a > b {
        let flipped = gauss_kronrod_any(f, b, a, tolerance, max_panels);
        return Quadrature {
            value: -flipped.value,
            ..flipped
        };
    }
    match (a.is_finite(), b.is_finite()) {
        | (true, true) => gauss_kronrod(f, a, b, tolerance, max_panels),
        | (true, false) => {
            // x = a + t/(1 - t), t in [0, 1)
            gauss_kronrod(
                |t| f(a + t / (1.0 - t)) / ((1.0 - t) * (1.0 - t)),
                0.0,
                1.0,
                tolerance,
                max_panels,
            )
        },
        | (false, true) => {
            // x = b - t/(1 - t), t in [0, 1)
            gauss_kronrod(
                |t| f(b - t / (1.0 - t)) / ((1.0 - t) * (1.0 - t)),
                0.0,
                1.0,
                tolerance,
                max_panels,
            )
        },
        | (false, false) => {
            // x = t/(1 - t^2), t in (-1, 1)
            gauss_kronrod(
                |t| {
                    let d = 1.0 - t * t;
                    f(t / d) * (1.0 + t * t) / (d * d)
                },
                -1.0,
                1.0,
                tolerance,
                max_panels,
            )
        },
    }
}

#[cfg(test)]
mod gauss_kronrod_tests {
    use super::*;

    #[test]
    fn smooth_integrands_converge_in_one_panel() {
        let q = gauss_kronrod(f64::sin, 0.0, std::f64::consts::PI, 1e-12, 100);
        assert!((q.value - 2.0).abs() < 1e-13, "{q:?}");
        assert!(q.evaluations <= 45, "{q:?}");
        // Degree 22 is integrated exactly by the Kronrod rule.
        let q = gauss_kronrod(|x| x.powi(22), -1.0, 1.0, 1e-14, 1);
        assert!((q.value - 2.0 / 23.0).abs() < 1e-15);
    }

    #[test]
    fn adaptivity_handles_peaks_and_endpoint_singularities() {
        let q = gauss_kronrod(|x| 1.0 / (1e-4 + x * x), -1.0, 1.0, 1e-9, 2_000);
        let exact = 2.0 * (1.0_f64 / 1e-2).atan() / 1e-2;
        assert!((q.value - exact).abs() < 1e-8 * exact, "{q:?} vs {exact}");
        let q = gauss_kronrod(f64::sqrt, 0.0, 1.0, 1e-10, 2_000);
        assert!((q.value - 2.0 / 3.0).abs() < 1e-9, "{q:?}");
        assert!(
            q.error >= (q.value - 2.0 / 3.0).abs() * 0.1,
            "the estimate is not wildly optimistic"
        );
    }

    #[test]
    fn infinite_intervals() {
        let q = gauss_kronrod_any(
            |x| (-x * x).exp(),
            f64::NEG_INFINITY,
            f64::INFINITY,
            1e-12,
            2_000,
        );
        assert!(
            (q.value - std::f64::consts::PI.sqrt()).abs() < 1e-10,
            "{q:?}"
        );
        let q = gauss_kronrod_any(|x| (-x).exp(), 0.0, f64::INFINITY, 1e-12, 2_000);
        assert!((q.value - 1.0).abs() < 1e-10, "{q:?}");
        let q = gauss_kronrod_any(
            |x| 1.0 / (1.0 + x * x),
            f64::NEG_INFINITY,
            0.0,
            1e-12,
            2_000,
        );
        assert!(
            (q.value - std::f64::consts::FRAC_PI_2).abs() < 1e-10,
            "{q:?}"
        );
        let q = gauss_kronrod_any(|x| x, 1.0, 0.0, 1e-12, 10);
        assert!(
            (q.value + 0.5).abs() < 1e-14,
            "reversed limits flip the sign"
        );
    }

    #[test]
    fn exhausted_budget_is_visible_in_the_error() {
        let q = gauss_kronrod(|x| (1.0 / x).sin(), 1e-6, 1.0, 1e-12, 8);
        assert!(q.error > 1e-6, "{q:?}");
    }

    #[test]
    fn reversed_finite_interval_retains_adaptivity() {
        let f = |x: f64| 1.0 / (1e-4 + x * x);
        let forward = gauss_kronrod(f, -1.0, 1.0, 1e-9, 2_000);
        let reverse = gauss_kronrod(f, 1.0, -1.0, 1e-9, 2_000);
        assert!(forward.evaluations > 15);
        assert_eq!(reverse.value, -forward.value);
        assert_eq!(reverse.error, forward.error);
        assert_eq!(reverse.evaluations, forward.evaluations);
    }

    #[test]
    fn non_finite_improper_integrand_is_not_silently_discarded() {
        let q = gauss_kronrod_any(|_| f64::NAN, 0.0, f64::INFINITY, 1e-9, 10);
        assert!(!q.value.is_finite(), "{q:?}");
    }

    #[test]
    fn wide_finite_interval_does_not_overflow_panel_width() {
        let q = gauss_kronrod(|_| 1e-308, -1e308, 1e308, 1e-12, 10);
        assert!((q.value - 2.0).abs() < 1e-14, "{q:?}");
    }

    #[test]
    fn unsplittable_panel_retains_its_estimate() {
        let a = 1.0_f64;
        let b = f64::from_bits(a.to_bits() + 1);
        let f = |x: f64| {
            if x.to_bits() == a.to_bits() {
                1.0
            } else {
                0.0
            }
        };
        let initial = gauss_kronrod(f, a, b, 0.0, 1);
        let refined = gauss_kronrod(f, a, b, 0.0, 10);
        assert_ne!(initial.value, 0.0);
        assert_eq!(refined.value, initial.value);
        assert_eq!(refined.error, f64::INFINITY);
    }
}
