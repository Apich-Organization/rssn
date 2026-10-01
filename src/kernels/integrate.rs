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

    if n_steps == 0 {
        return 0.0;
    }

    // Correctly handle a == b (integral is 0)
    // and a > b (integral is negative of b to a)
    // The logic below works for a > b as h will be negative.
    if (a - b).abs() < f64::EPSILON {
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

    if n_steps == 0 {
        return Ok(0.0);
    }

    if (a - b).abs() < f64::EPSILON {
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
        let x = (i as f64).mul_add(h, a);

        let weight = if i % 2 == 0 { 2.0 } else { 4.0 };

        sum += weight * f(x);
    }

    Ok((h / 3.0) * sum)
}

/// # Adaptive Quadrature (Adaptive Simpson's Method)
///
/// Recursively refines the interval until the estimated error satisfies the tolerance.
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
    // Inner recursive function
    fn adaptive_recursive<F>(
        f: &F,
        a: f64,
        b: f64,
        eps: f64,
        whole_simpson: f64,
        limit: usize,
    ) -> f64
    where
        F: Fn(f64) -> f64,
    {
        if limit == 0 {
            // Recursion limit reached, return current best guess
            return whole_simpson;
        }

        let mid = f64::midpoint(a, b);

        let sub_mid_left = f64::midpoint(a, mid);

        let sub_mid_right = f64::midpoint(mid, b);

        let fa = f(a);

        let fb = f(b);

        let fm = f(mid);

        let fml = f(sub_mid_left);

        let fmr = f(sub_mid_right);

        // Simp(a, b) = (b-a)/6 * (f(a) + 4f(m) + f(b))
        let left_simpson = (mid - a) / 6.0 * (4.0f64.mul_add(fml, fa) + fm);

        let right_simpson = (b - mid) / 6.0 * (4.0f64.mul_add(fmr, fm) + fb);

        let sum_halves = left_simpson + right_simpson;

        // Error estimate (1/15 rule)
        let error = (sum_halves - whole_simpson).abs() / 15.0;

        if error <= eps {
            // Richardson extrapolation: S + (S - S_whole)/15
            sum_halves + (sum_halves - whole_simpson) / 15.0
        } else {
            adaptive_recursive(f, a, mid, eps / 2.0, left_simpson, limit - 1)
                + adaptive_recursive(f, mid, b, eps / 2.0, right_simpson, limit - 1)
        }
    }

    let (a, b) = range;

    if (a - b).abs() < f64::EPSILON {
        return 0.0;
    }

    // Initial Simpson estimate
    let mid = f64::midpoint(a, b);

    let fm = f(mid);

    let initial_simpson = (b - a) / 6.0 * (4.0f64.mul_add(fm, f(a)) + f(b));

    adaptive_recursive(&f, a, b, tolerance, initial_simpson, 100) // limit depth to avoid stack overflow
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

    if max_steps == 0 {
        return 0.0;
    }

    if (a - b).abs() < f64::EPSILON {
        return 0.0;
    }

    let mut r = vec![vec![0.0; max_steps]; max_steps];

    // R[0][0]
    let h = b - a;

    r[0][0] = 0.5 * h * (f(a) + f(b));

    for i in 1..max_steps {
        // Calculate R[i][0] using Trapezoidal rule with 2^i segments
        // But we can update from R[i-1][0] efficiently
        let steps_prev = 1 << (i - 1);

        let h_i = h / f64::from(1 << i);

        let mut sum = 0.0;

        for k in 1..=steps_prev {
            let x = a + f64::from(2 * k - 1) * h_i;

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

    if (a - b).abs() < f64::EPSILON {
        return 0.0;
    }

    let mid = f64::midpoint(a, b);

    let half_len = (b - a) / 2.0;

    // Nodes and weights for n=5 (from standard tables)
    // x_i are for interval [-1, 1]
    let nodes = [
        0.0,
        0.538_469_310_105_683_1,
        -0.538_469_310_105_683_1,
        0.906_179_845_938_664,
        -0.906_179_845_938_664,
    ];

    let weights = [
        0.568_888_888_888_889,
        0.478_628_670_499_366_5,
        0.478_628_670_499_366_5,
        0.236_926_885_056_189_1,
        0.236_926_885_056_189_1,
    ];

    let mut sum = 0.0;

    for i in 0..5 {
        // Transform x from [-1, 1] to [a, b]
        let x = mid + half_len * nodes[i];

        sum += weights[i] * f(x);
    }

    half_len * sum
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
const GAUSS_WEIGHTS: [f64; 4] =
    [0.129_484_966_168_869_7, 0.279_705_391_489_276_64, 0.381_830_050_505_118_9, 0.417_959_183_673_469_4];

/// One Gauss–Kronrod 7/15 panel on `[a, b]`: `(estimate, error estimate)`.
fn gk15(
    f: &impl Fn(f64) -> f64,
    a: f64,
    b: f64,
) -> (f64, f64) {
    let centre = f64::midpoint(a, b);
    let half = 0.5 * (b - a);
    let mut kronrod = 0.0;
    let mut gauss = 0.0;
    for (i, (&node, &weight)) in GK_NODES.iter().zip(&KRONROD_WEIGHTS).enumerate() {
        let sum = if node == 0.0 { f(centre) } else { f(centre - half * node) + f(centre + half * node) };
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
#[allow(clippy::float_cmp)] // exact comparison against a sentinel / integer-valued input is intended
pub fn gauss_kronrod(
    f: impl Fn(f64) -> f64,
    a: f64,
    b: f64,
    tolerance: f64,
    max_panels: usize,
) -> Quadrature {
    if a == b {
        return Quadrature { value: 0.0, error: 0.0, evaluations: 0 };
    }
    let (value, error) = gk15(&f, a, b);
    let mut panels = vec![(a, b, value, error)];
    let mut evaluations = 15;
    loop {
        let total_error: f64 = panels.iter().map(|p| p.3).sum();
        let total: f64 = panels.iter().map(|p| p.2).sum();
        if total_error <= tolerance.max(f64::EPSILON * total.abs()) || panels.len() >= max_panels || !total.is_finite() {
            return Quadrature { value: total, error: total_error, evaluations };
        }
        let worst = panels
            .iter()
            .enumerate()
            .max_by(|x, y| x.1.3.total_cmp(&y.1.3))
            .map_or(0, |(i, _)| i);
        let (lo, hi, _, _) = panels.swap_remove(worst);
        let mid = f64::midpoint(lo, hi);
        if mid <= lo || mid >= hi {
            // The panel cannot be split any further in floating point.
            panels.push((lo, hi, 0.0, 0.0));
            let total: f64 = panels.iter().map(|p| p.2).sum();
            return Quadrature { value: total, error: f64::INFINITY, evaluations };
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
        return Quadrature { value: -flipped.value, ..flipped };
    }
    // Guard the transformed integrands against the end points, where the
    // Jacobian is infinite but the product is (for a convergent integral)
    // zero.
    let guard = |v: f64| if v.is_finite() { v } else { 0.0 };
    match (a.is_finite(), b.is_finite()) {
        | (true, true) => gauss_kronrod(f, a, b, tolerance, max_panels),
        | (true, false) => {
            // x = a + t/(1 - t), t in [0, 1)
            gauss_kronrod(
                |t| guard(f(a + t / (1.0 - t)) / ((1.0 - t) * (1.0 - t))),
                0.0,
                1.0,
                tolerance,
                max_panels,
            )
        },
        | (false, true) => {
            // x = b - t/(1 - t), t in [0, 1)
            gauss_kronrod(
                |t| guard(f(b - t / (1.0 - t)) / ((1.0 - t) * (1.0 - t))),
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
                    guard(f(t / d) * (1.0 + t * t) / (d * d))
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
        assert!(q.error >= (q.value - 2.0 / 3.0).abs() * 0.1, "the estimate is not wildly optimistic");
    }

    #[test]
    fn infinite_intervals() {
        let q = gauss_kronrod_any(|x| (-x * x).exp(), f64::NEG_INFINITY, f64::INFINITY, 1e-12, 2_000);
        assert!((q.value - std::f64::consts::PI.sqrt()).abs() < 1e-10, "{q:?}");
        let q = gauss_kronrod_any(|x| (-x).exp(), 0.0, f64::INFINITY, 1e-12, 2_000);
        assert!((q.value - 1.0).abs() < 1e-10, "{q:?}");
        let q = gauss_kronrod_any(|x| 1.0 / (1.0 + x * x), f64::NEG_INFINITY, 0.0, 1e-12, 2_000);
        assert!((q.value - std::f64::consts::FRAC_PI_2).abs() < 1e-10, "{q:?}");
        let q = gauss_kronrod_any(|x| x, 1.0, 0.0, 1e-12, 10);
        assert!((q.value + 0.5).abs() < 1e-14, "reversed limits flip the sign");
    }

    #[test]
    fn exhausted_budget_is_visible_in_the_error() {
        let q = gauss_kronrod(|x| (1.0 / x).sin(), 1e-6, 1.0, 1e-12, 8);
        assert!(q.error > 1e-6, "{q:?}");
    }
}
