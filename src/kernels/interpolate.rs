//! # Numerical Interpolation
//!
//! This module provides numerical interpolation methods for constructing new data points
//! within the range of a discrete set of known data points. It includes implementations
//! for Lagrange interpolation, piecewise-linear interpolation, cubic spline interpolation
//! with selectable end conditions, monotone cubic Hermite (PCHIP) interpolation, and
//! Bézier and B-spline curves.
//!
//! # Order
//!
//! [`linear_interpolation`] is second order in the node spacing `h`,
//! [`cubic_spline_with`] is fourth order when its end condition agrees with the data,
//! [`pchip_interpolation`] is third order away from extrema of the data and second order
//! at one, and [`lagrange_interpolation`] is exact for polynomials of degree `n - 1` but
//! is a global interpolant that does not improve locally with `h`. The natural end
//! condition (the default of the legacy [`cubic_spline_interpolation`]) forces a zero
//! second derivative at both ends, which a general smooth function does not have, so a
//! natural spline is only second order.
//!
//! # Extrapolation
//!
//! Every interpolant here **extrapolates**: a query outside `[x_0, x_n]` continues the
//! first or the last segment polynomial. Nothing is clamped, refused or reported, and no
//! error estimate is returned, so a caller that needs a defined range must check it.
//!
//! ## Limitations
//!
//! The spline degree is not configurable: [`cubic_spline_with`] selects the *end
//! condition*, not the order, and there is no quintic or higher-order spline. Monotone
//! interpolation is PCHIP only — `C^1` piecewise cubic — with no `C^2` monotone spline and
//! no monotonicity-constrained least squares. The abscissae must be strictly increasing
//! and are not sorted for the caller; the legacy [`cubic_spline_interpolation`] checks
//! only the count, so repeated abscissae reach a division by zero. The legacy closure
//! locates its interval by a linear scan, `O(n)` per evaluation, where the newer
//! interpolants binary-search. The parametric curves are not function interpolants:
//! [`b_spline`] requires a non-empty control polygon and a consistent knot vector, and
//! [`bezier_curve`] recurses once per control point.

use std::sync::Arc;

use crate::kernels::approx::ApproxError;
use crate::kernels::approx::CubicSpline;
use crate::kernels::approx::Pchip;
use crate::kernels::approx::SplineBc;
use crate::kernels::polynomial::Polynomial;

/// Constructs a Lagrange interpolating polynomial that passes through a given set of points.
///
/// Lagrange interpolation is a method of finding a polynomial that takes on certain values
/// at specified points. The resulting polynomial is unique and passes exactly through all given points.
///
/// # Arguments
/// * `points` - A slice of `(x, y)` tuples.
///
/// # Returns
/// A `Result` containing a `Polynomial` struct representing the interpolating polynomial.
///
/// # Errors
/// Returns an error if duplicate x-coordinates are found among the input points.
///
/// # Example
/// ```rust
/// use rssn::kernels::interpolate::lagrange_interpolation;
///
/// let points = vec![(0.0, 0.0), (1.0, 1.0), (2.0, 4.0)];
///
/// let poly = lagrange_interpolation(&points).unwrap();
///
/// assert!((poly.eval(1.5) - 2.25).abs() < 1e-9);
/// ```
pub fn lagrange_interpolation(points: &[(f64, f64)]) -> Result<Polynomial, String> {
    if points.is_empty() {
        return Ok(Polynomial { coeffs: vec![0.0] });
    }

    let mut total_poly = Polynomial { coeffs: vec![0.0] };

    for (j, (xj, yj)) in points.iter().enumerate() {
        let mut basis_poly = Polynomial { coeffs: vec![1.0] };

        for (i, (xi, _)) in points.iter().enumerate() {
            if i == j {
                continue;
            }

            let numerator = Polynomial {
                coeffs: vec![1.0, -xi],
            };

            let denominator = xj - xi;

            if denominator == 0.0 || denominator.is_nan() {
                return Err(format!("Duplicate or non-finite x-coordinates found: {xj}"));
            }

            basis_poly = basis_poly * (numerator / denominator);
        }

        total_poly = total_poly + (basis_poly * *yj);
    }

    Ok(total_poly)
}

/// Creates a cubic spline interpolator for a given set of points.
///
/// Cubic spline interpolation constructs a piecewise cubic polynomial that passes
/// through each data point, ensuring smoothness (continuity of first and second derivatives)
/// at the interior points. The input points must be sorted by their x-coordinates.
///
/// # Arguments
/// * `points` - A slice of `(x, y)` tuples. Must be sorted by x.
///
/// # Returns
/// A `Result` containing a closure `Arc<dyn Fn(f64) -> f64>` that can be used to evaluate
/// the spline at any point.
///
/// # Errors
/// Returns an error if fewer than two points are provided, or if the x-coordinates are not strictly increasing.
///
/// # Example
/// ```rust
/// use rssn::kernels::interpolate::cubic_spline_interpolation;
///
/// let points = vec![(0.0, 0.0), (1.0, 1.0), (2.0, 0.0)];
///
/// let spline = cubic_spline_interpolation(&points).unwrap();
///
/// // Check that interpolated value is reasonable (between the y values)
/// let val = spline(0.5);
///
/// assert!((spline(0.5) - 0.6875).abs() < 1e-9);
/// assert!(val > 0.0 && val < 1.0);
/// ```
pub fn cubic_spline_interpolation(
    points: &[(f64, f64)]
) -> Result<Arc<dyn Fn(f64) -> f64>, String> {
    let n = points.len();

    if n < 2 {
        return Err("At least two \
                    points are \
                    required."
            .to_string());
    }

    if points
        .iter()
        .any(|&(x, y)| !x.is_finite() || !y.is_finite())
    {
        return Err("The data must be finite.".to_string());
    }

    let mut h = vec![0.0; n - 1];

    for i in 0..n - 1 {
        h[i] = points[i + 1].0 - points[i].0;

        if h[i] <= 0.0 {
            return Err("The x-coordinates must be strictly increasing.".to_string());
        }
    }

    // 1. Calculate Alpha (size n)
    let mut alpha = vec![0.0; n];

    for i in 1..n - 1 {
        alpha[i] = (3.0f64 / h[i]).mul_add(
            points[i + 1].1 - points[i].1,
            -(3.0f64 / h[i - 1]) * (points[i].1 - points[i - 1].1),
        );
    }

    // 2. Solve Tridiagonal System
    let mut l = vec![1.0; n];

    let mut mu = vec![0.0; n];

    let mut z = vec![0.0; n];

    for i in 1..n - 1 {
        // The diagonal element is 2*(h[i-1] + h[i])
        l[i] = 2.0f64.mul_add(h[i - 1] + h[i], -(h[i - 1] * mu[i - 1]));

        mu[i] = h[i] / l[i];

        z[i] = h[i - 1].mul_add(-z[i - 1], alpha[i]) / l[i];
    }

    let mut c = vec![0.0; n];

    let mut b = vec![0.0; n - 1];

    let mut d = vec![0.0; n - 1];

    // 3. Back-substitution
    for j in (0..n - 1).rev() {
        c[j] = (mu[j] as f64).mul_add(-c[j + 1], z[j]);

        b[j] = (points[j + 1].1 - points[j].1) / h[j] - h[j] * 2.0f64.mul_add(c[j], c[j + 1]) / 3.0;

        d[j] = (c[j + 1] - c[j]) / (3.0 * h[j]);
    }

    let points_owned = points.to_vec();

    let spline = move |x: f64| -> f64 {
        // 4. Find the correct interval
        let mut i = 0;

        for j in 0..points_owned.len() - 1 {
            if x >= points_owned[j].0 && x <= points_owned[j + 1].0 {
                i = j;

                break;
            }
        }

        let dx = x - points_owned[i].0;

        dx.mul_add(dx.mul_add(dx.mul_add(d[i], c[i]), b[i]), points_owned[i].1)
    };

    Ok(Arc::new(spline))
}

/// End conditions of [`cubic_spline_with`].
///
/// All three conditions interpolate the data and make the spline `C^2` at the interior
/// knots; they differ in the two extra equations that close the system, hence in the
/// first and last intervals.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SplineEnd {
    /// Zero second derivative at both ends.
    Natural,
    /// Prescribed first derivatives `(left, right)` at the two ends.
    Clamped(f64, f64),
    /// Continuous third derivative across the second and the last-but-one knot.
    ///
    /// Needs at least four nodes and reproduces cubics exactly without derivative data.
    NotAKnot,
}

/// Message for a failure of the underlying approximation kernel.
fn approx_message(err: ApproxError) -> String {
    match err {
        | ApproxError::Invalid => {
            "the x-coordinates must be strictly increasing and \
                                   there must be at least two points (three for a \
                                   spline, four for not-a-knot)"
                .to_string()
        },
        | ApproxError::Singular => "the interpolation system is singular".to_string(),
    }
}

/// Piecewise-linear interpolation through `(x, y)`.
///
/// Second order in the node spacing for a smooth function: the error on a segment of
/// length `h` is `f''(xi) h^2 / 8`. Values outside the node range are extrapolated with
/// the first or the last segment.
///
/// # Arguments
/// * `points` - A slice of `(x, y)` tuples with strictly increasing `x`.
///
/// # Returns
/// A `Result` containing a closure `Arc<dyn Fn(f64) -> f64>` evaluating the interpolant.
///
/// # Errors
/// Returns an error if fewer than two points are given or the abscissae are not strictly
/// increasing.
///
/// # Example
/// ```rust
/// use rssn::kernels::interpolate::linear_interpolation;
///
/// let f = linear_interpolation(&[(0.0, 0.0), (1.0, 1.0), (2.0, 4.0)]).unwrap();
///
/// assert!((f(0.5) - 0.5).abs() < 1e-12);
///
/// assert!((f(1.5) - 2.5).abs() < 1e-12);
/// ```
pub fn linear_interpolation(points: &[(f64, f64)]) -> Result<Arc<dyn Fn(f64) -> f64>, String> {
    if points.len() < 2 {
        return Err("At least two points are required.".to_string());
    }

    if !points.windows(2).all(|w| {
        match w {
            | [a, b] => a.0 < b.0,
            | _ => false,
        }
    }) {
        return Err("The x-coordinates must be strictly increasing.".to_string());
    }

    let xs: Vec<f64> = points.iter().map(|p| p.0).collect();

    let ys: Vec<f64> = points.iter().map(|p| p.1).collect();

    let slopes: Vec<f64> = points
        .windows(2)
        .map(|w| {
            match w {
                | [a, b] => (b.1 - a.1) / (b.0 - a.0),
                | _ => f64::NAN,
            }
        })
        .collect();

    let eval = move |x: f64| -> f64 {
        let i = xs
            .partition_point(|&v| v <= x)
            .saturating_sub(1)
            .min(xs.len().saturating_sub(2));

        match (xs.get(i), ys.get(i), slopes.get(i)) {
            | (Some(&x0), Some(&y0), Some(&m)) => y0 + m * (x - x0),
            | _ => f64::NAN,
        }
    };

    Ok(Arc::new(eval))
}

/// Cubic spline interpolation with a selectable end condition.
///
/// Fourth order in the node spacing when the end condition matches the data — clamped
/// with the exact end slopes, or not-a-knot — and second order for the natural condition,
/// whose zero end curvature a general smooth function violates. Values outside the node
/// range are extrapolated with the first or the last cubic.
///
/// # Arguments
/// * `points` - A slice of `(x, y)` tuples with strictly increasing `x`; at least three
///   of them, four for [`SplineEnd::NotAKnot`].
/// * `ends` - The end condition.
///
/// # Returns
/// A `Result` containing a closure `Arc<dyn Fn(f64) -> f64>` evaluating the spline.
///
/// # Errors
/// Returns an error for too few points, non-increasing abscissae, or a singular system.
///
/// # Example
/// ```rust
/// use rssn::kernels::interpolate::SplineEnd;
/// use rssn::kernels::interpolate::cubic_spline_with;
///
/// let points = [(0.0, 0.0), (1.0, 1.0), (2.0, 0.0), (3.0, 1.0)];
///
/// let spline = cubic_spline_with(&points, SplineEnd::NotAKnot).unwrap();
///
/// assert!((spline(1.0) - 1.0).abs() < 1e-12);
/// ```
pub fn cubic_spline_with(
    points: &[(f64, f64)],
    ends: SplineEnd,
) -> Result<Arc<dyn Fn(f64) -> f64>, String> {
    let (xs, ys): (Vec<f64>, Vec<f64>) = points.iter().copied().unzip();

    let bc = match ends {
        | SplineEnd::Natural => SplineBc::Natural,
        | SplineEnd::Clamped(left, right) => SplineBc::Clamped(left, right),
        | SplineEnd::NotAKnot => SplineBc::NotAKnot,
    };

    let spline = CubicSpline::new(&xs, &ys, bc).map_err(approx_message)?;

    Ok(Arc::new(move |x: f64| spline.eval(x)))
}

/// Monotone piecewise cubic Hermite interpolation (PCHIP, Fritsch-Carlson).
///
/// The interpolant is `C^1` and, unlike a plain cubic spline, preserves the shape of the
/// data: on every interval its values stay inside the range of the two end ordinates, so
/// monotone data stay monotone and a step is not overshot. Third order in the node
/// spacing where the slope has one sign, second order at a data extremum, where the
/// slope estimate is clipped to zero — the price of shape preservation. Values outside
/// the node range are extrapolated with the first or the last Hermite cubic.
///
/// # Arguments
/// * `points` - A slice of `(x, y)` tuples with strictly increasing `x`.
///
/// # Returns
/// A `Result` containing a closure `Arc<dyn Fn(f64) -> f64>` evaluating the interpolant.
///
/// # Errors
/// Returns an error if fewer than two points are given or the abscissae are not strictly
/// increasing.
///
/// # Example
/// ```rust
/// use rssn::kernels::interpolate::pchip_interpolation;
///
/// // A step: the monotone interpolant stays within [0, 1].
/// let step = [(0.0, 0.0), (1.0, 0.0), (2.0, 1.0), (3.0, 1.0)];
///
/// let f = pchip_interpolation(&step).unwrap();
///
/// for k in 0..=30 {
///     let v = f(f64::from(k) / 10.0);
///
///     assert!((-1e-12..=1.0 + 1e-12).contains(&v));
/// }
/// ```
pub fn pchip_interpolation(points: &[(f64, f64)]) -> Result<Arc<dyn Fn(f64) -> f64>, String> {
    let (xs, ys): (Vec<f64>, Vec<f64>) = points.iter().copied().unzip();

    let pchip = Pchip::new(&xs, &ys).map_err(approx_message)?;

    Ok(Arc::new(move |x: f64| pchip.eval(x)))
}

/// Evaluates a point on a Bézier curve defined by a set of control points at parameter `t`.
///
/// This function uses De Casteljau's algorithm, a numerically stable and efficient method
/// for evaluating Bézier curves. The curve passes through the first and last control points
/// and is tangent to the segments connecting the first two and last two control points.
///
/// # Arguments
/// * `control_points` - A slice of points, where each point is a slice of `f64` coordinates.
/// * `t` - The parameter, typically in the range `[0, 1]`.
///
/// # Returns
/// The coordinates of the point on the curve.
///
/// # Example
/// ```rust
/// use rssn::kernels::interpolate::bezier_curve;
///
/// let control_points = vec![vec![0.0, 0.0], vec![1.0, 2.0], vec![2.0, 0.0]];
///
/// let p = bezier_curve(&control_points, 0.5);
///
/// assert!((p[0] - 1.0).abs() < 1e-9);
///
/// assert!((p[1] - 1.0).abs() < 1e-9);
/// ```
#[must_use]
pub fn bezier_curve(
    control_points: &[Vec<f64>],
    t: f64,
) -> Vec<f64> {
    if control_points.is_empty() {
        return vec![];
    }

    if control_points.len() == 1 {
        return control_points[0].clone();
    }

    let mut new_points = Vec::with_capacity(control_points.len() - 1);

    for i in 0..(control_points.len() - 1) {
        let p1 = &control_points[i];

        let p2 = &control_points[i + 1];

        let new_point: Vec<f64> = p1
            .iter()
            .zip(p2.iter())
            .map(|(&c1, &c2)| (1.0f64 - t).mul_add(c1, t * c2))
            .collect();

        new_points.push(new_point);
    }

    bezier_curve(&new_points, t)
}

/// Evaluates a point on a B-spline curve.
///
/// B-splines are generalizations of Bézier curves, offering more flexibility and local control.
/// They are defined by control points, a degree, and a knot vector. This function uses
/// the Cox-de Boor recursion formula to evaluate the curve at a given parameter `t`.
///
/// # Arguments
/// * `control_points` - The control points.
/// * `degree` - The degree of the spline, `p`.
/// * `knots` - The knot vector, `U`.
/// * `t` - The parameter.
///
/// # Returns
/// The coordinates of the point on the curve.
///
/// # Example
/// ```rust
/// use rssn::kernels::interpolate::b_spline;
///
/// let control_points = vec![vec![0.0], vec![1.0], vec![2.0]];
///
/// let knots = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
///
/// let p = b_spline(&control_points, 2, &knots, 0.5);
///
/// assert!((p.unwrap()[0] - 1.0).abs() < 1e-9);
/// ```
#[must_use]
pub fn b_spline(
    control_points: &[Vec<f64>],
    degree: usize,
    knots: &[f64],
    t: f64,
) -> Option<Vec<f64>> {
    let n = control_points.len() - 1;

    let m = knots.len() - 1;

    if degree > n || m != n + degree + 1 {
        return None;
    }

    let i = find_knot_span(n, degree, t, knots);

    let basis_vals = basis_functions(i, t, degree, knots);

    let mut point = vec![0.0; control_points[0].len()];

    for (j, _var) in basis_vals.iter().enumerate().take(degree + 1) {
        let pt_idx = i - degree + j;

        let p = &control_points[pt_idx];

        for d in 0..point.len() {
            point[d] += basis_vals[j] * p[d];
        }
    }

    Some(point)
}

/// Finds the knot span for a given parameter t.
pub(crate) fn find_knot_span(
    n: usize,
    p: usize,
    t: f64,
    knots: &[f64],
) -> usize {
    if t >= knots[n + 1] {
        return n;
    }

    if t < knots[p] {
        return p;
    }

    let mut low = p;

    let mut high = n + 1;

    let mut mid = usize::midpoint(low, high);

    while t < knots[mid] || t >= knots[mid + 1] {
        if t < knots[mid] {
            high = mid;
        } else {
            low = mid;
        }

        mid = usize::midpoint(low, high);
    }

    mid
}

/// Computes the B-spline basis functions using the Cox-de Boor formula.
pub(crate) fn basis_functions(
    i: usize,
    t: f64,
    p: usize,
    knots: &[f64],
) -> Vec<f64> {
    let mut n = vec![0.0; p + 1];

    let mut left = vec![0.0; p + 1];

    let mut right = vec![0.0; p + 1];

    n[0] = 1.0;

    for j in 1..=p {
        left[j] = t - knots[i + 1 - j];

        right[j] = knots[i + j] - t;

        let mut saved = 0.0;

        for r in 0..j {
            let temp = n[r] / (right[r + 1] + left[j - r]);

            n[r] = right[r + 1].mul_add(temp, saved);

            saved = left[j - r] * temp;
        }

        n[j] = saved;
    }

    n
}

#[cfg(test)]
mod tests {
    #![allow(
        clippy::unwrap_used,
        clippy::expect_used,
        clippy::indexing_slicing,
        clippy::arithmetic_side_effects,
        clippy::cast_precision_loss
    )]

    use super::*;

    /// Smooth test function; not a polynomial, so every order is visible.
    fn exact(x: f64) -> f64 {
        (3.0 * x).sin() + 0.5 * x
    }

    fn samples(
        a: f64,
        b: f64,
        count: usize,
    ) -> impl Iterator<Item = f64> {
        (0..=count).map(move |k| a + (b - a) * k as f64 / count as f64)
    }

    /// Maximum deviation of `f` from `g` on `[a, b]`.
    fn max_err(
        f: &Arc<dyn Fn(f64) -> f64>,
        g: &dyn Fn(f64) -> f64,
        a: f64,
        b: f64,
    ) -> f64 {
        samples(a, b, 4000)
            .map(|x| (f(x) - g(x)).abs())
            .fold(0.0_f64, f64::max)
    }

    /// Least-squares slope of `(ln h, ln e)`.
    fn slope(points: &[(f64, f64)]) -> f64 {
        let n = points.len() as f64;
        let mx = points.iter().map(|p| p.0).sum::<f64>() / n;
        let my = points.iter().map(|p| p.1).sum::<f64>() / n;
        let num: f64 = points.iter().map(|p| (p.0 - mx) * (p.1 - my)).sum();
        let den: f64 = points.iter().map(|p| (p.0 - mx).powi(2)).sum();
        num / den
    }

    /// Measures the convergence order of `build` on `n = 16 ... 512` nodes of `exact`.
    fn measured_order<B>(
        build: B,
        exact: &dyn Fn(f64) -> f64,
    ) -> f64
    where
        B: Fn(&[(f64, f64)]) -> Arc<dyn Fn(f64) -> f64>,
    {
        let points: Vec<(f64, f64)> = [16_usize, 32, 64, 128, 256, 512]
            .iter()
            .map(|&n| {
                let data: Vec<(f64, f64)> = (0..n)
                    .map(|i| {
                        let x = i as f64 / (n - 1) as f64;

                        (x, exact(x))
                    })
                    .collect();

                let f = build(&data);

                let e = max_err(&f, exact, 0.0, 1.0);

                (-((n - 1) as f64).ln(), e.ln())
            })
            .collect();

        slope(&points)
    }

    #[test]
    fn interp_linear_is_second_order() {
        let order = measured_order(|data| linear_interpolation(data).unwrap(), &exact);

        eprintln!("linear order = {order:.3}");

        assert!((1.9..2.1).contains(&order), "linear order {order}");
    }

    #[test]
    fn interp_cubic_is_fourth_order() {
        let order = measured_order(
            |data| cubic_spline_with(data, SplineEnd::NotAKnot).unwrap(),
            &exact,
        );

        eprintln!("cubic (not-a-knot) order = {order:.3}");

        assert!((3.7..4.5).contains(&order), "cubic order {order}");
    }

    #[test]
    fn interp_natural_is_only_second_order() {
        let order = measured_order(
            |data| cubic_spline_with(data, SplineEnd::Natural).unwrap(),
            &exact,
        );

        eprintln!("cubic (natural) order = {order:.3}");

        assert!((1.9..2.2).contains(&order), "natural order {order}");
    }

    #[test]
    fn interp_pchip_is_third_order() {
        // A monotone function: `exp` has no extremum, so the slope estimate is never
        // clipped to zero and the full third order shows. With `exact` (whose derivative
        // changes sign inside the range) the same test measures two, which is the
        // shape-preserving trade-off, not a defect.
        let order = measured_order(|data| pchip_interpolation(data).unwrap(), &|x: f64| x.exp());

        eprintln!("pchip order = {order:.3}");

        assert!((2.7..3.4).contains(&order), "pchip order {order}");
    }

    #[test]
    fn interp_pchip_preserves_monotone_data() {
        let data: Vec<(f64, f64)> = (0..7)
            .map(|i| (f64::from(i), f64::from(i).powi(3)))
            .collect();

        let f = pchip_interpolation(&data).unwrap();

        let mut previous = f(0.0);

        for x in samples(0.0, 6.0, 6000) {
            let v = f(x);

            assert!(v >= previous - 1e-9, "decrease at {x}: {previous} -> {v}");

            previous = v;
        }
    }

    #[test]
    fn interp_pchip_no_overshoot_on_a_step() {
        let data = [
            (0.0, 0.0),
            (1.0, 0.0),
            (2.0, 0.0),
            (3.0, 1.0),
            (4.0, 1.0),
            (5.0, 1.0),
        ];

        let f = pchip_interpolation(&data).unwrap();

        for pair in data.windows(2) {
            let (a, b) = (pair[0], pair[1]);

            let (lo, hi) = (a.1.min(b.1), a.1.max(b.1));

            for x in samples(a.0, b.0, 400) {
                let v = f(x);

                assert!(
                    v >= lo - 1e-12 && v <= hi + 1e-12,
                    "overshoot at {x}: {v} outside [{lo}, {hi}]"
                );
            }
        }
    }

    #[test]
    fn interp_plain_cubic_overshoots_the_step() {
        // Control: the same step under the natural spline leaves the data range.
        let data = [
            (0.0, 0.0),
            (1.0, 0.0),
            (2.0, 0.0),
            (3.0, 1.0),
            (4.0, 1.0),
            (5.0, 1.0),
        ];

        let f = cubic_spline_interpolation(&data).unwrap();

        let (lo, hi) = samples(0.0, 5.0, 4000)
            .map(|x| f(x))
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), v| {
                (l.min(v), h.max(v))
            });

        eprintln!("natural spline on the step: [{lo:.6}, {hi:.6}]");

        assert!(hi > 1.05, "expected an overshoot, got {hi}");

        assert!(lo < -0.05, "expected an undershoot, got {lo}");
    }

    #[test]
    fn interp_clamped_and_not_a_knot_reproduce_a_cubic() {
        let f = |t: f64| 2.0 * t.powi(3) - 3.0 * t.powi(2) + t + 1.0;

        let df = |t: f64| 6.0 * t * t - 6.0 * t + 1.0;

        let data: Vec<(f64, f64)> = (0..11)
            .map(|i| {
                let x = 3.0 * f64::from(i) / 10.0;

                (x, f(x))
            })
            .collect();

        let clamped = cubic_spline_with(&data, SplineEnd::Clamped(df(0.0), df(3.0))).unwrap();

        let not_a_knot = cubic_spline_with(&data, SplineEnd::NotAKnot).unwrap();

        let e_clamped = max_err(&clamped, &f, 0.0, 3.0);

        let e_nak = max_err(&not_a_knot, &f, 0.0, 3.0);

        eprintln!("cubic reproduction: clamped {e_clamped:.3e}, not-a-knot {e_nak:.3e}");

        assert!(e_clamped < 1e-12, "clamped error {e_clamped}");

        assert!(e_nak < 1e-10, "not-a-knot error {e_nak}");
    }

    #[test]
    fn interp_natural_does_not_reproduce_a_cubic() {
        // Natural sets f'' = 0 at both ends, a cubic does not have that.
        let f = |t: f64| 2.0 * t.powi(3) - 3.0 * t.powi(2) + t + 1.0;

        let data: Vec<(f64, f64)> = (0..11)
            .map(|i| {
                let x = 3.0 * f64::from(i) / 10.0;

                (x, f(x))
            })
            .collect();

        let natural = cubic_spline_with(&data, SplineEnd::Natural).unwrap();

        let e = max_err(&natural, &f, 0.0, 3.0);

        eprintln!("natural on a cubic: {e:.3e}");

        assert!(e > 1e-2, "natural error {e}");
    }

    #[test]
    fn interp_end_conditions_differ_only_near_the_ends() {
        let f = |t: f64| 2.0 * t.powi(3) - 3.0 * t.powi(2) + t + 1.0;

        let df = |t: f64| 6.0 * t * t - 6.0 * t + 1.0;

        let data: Vec<(f64, f64)> = (0..41)
            .map(|i| {
                let x = 3.0 * f64::from(i) / 40.0;

                (x, f(x))
            })
            .collect();

        let natural = cubic_spline_with(&data, SplineEnd::Natural).unwrap();

        let clamped = cubic_spline_with(&data, SplineEnd::Clamped(df(0.0), df(3.0))).unwrap();

        let diff = |a: f64, b: f64| {
            samples(a, b, 400)
                .map(|x| (natural(x) - clamped(x)).abs())
                .fold(0.0_f64, f64::max)
        };

        let left = diff(0.0, 3.0 / 40.0);

        let right = diff(117.0 / 40.0, 3.0);

        // The two conditions differ by a decaying end effect: the central half is far
        // below the end intervals, and the decay is geometric in the knot index.
        let middle = diff(0.75, 2.25);

        eprintln!("natural - clamped: left {left:.3e}, right {right:.3e}, middle {middle:.3e}");

        assert!(left > 1e-3, "left end differs by only {left}");

        assert!(right > 1e-3, "right end differs by only {right}");

        assert!(middle < 1e-6, "the middle differs by {middle}");

        assert!(middle * 1000.0 < left, "no decay: {middle} vs {left}");
    }

    #[test]
    fn interp_extrapolation_is_documented() {
        // Linear data: every interpolant continues the end segment exactly.
        let line = [(0.0, 1.0), (1.0, 3.0), (2.0, 5.0)];

        let linear = linear_interpolation(&line).unwrap();

        let natural = cubic_spline_with(&line, SplineEnd::Natural).unwrap();

        let pchip = pchip_interpolation(&line).unwrap();

        for f in [linear, natural, pchip] {
            assert!(
                (f(-2.0) + 3.0).abs() < 1e-12,
                "left extrapolation {}",
                f(-2.0)
            );

            assert!(
                (f(4.0) - 9.0).abs() < 1e-12,
                "right extrapolation {}",
                f(4.0)
            );
        }

        // Cubic data with clamped ends: the end cubic is the exact polynomial.
        let g = |t: f64| t.powi(3) - 2.0 * t + 1.0;

        let dg = |t: f64| 3.0 * t * t - 2.0;

        let data: Vec<(f64, f64)> = (0..8)
            .map(|i| {
                let x = 0.5 * f64::from(i);

                (x, g(x))
            })
            .collect();

        let clamped = cubic_spline_with(&data, SplineEnd::Clamped(dg(0.0), dg(3.5))).unwrap();

        assert!((clamped(-1.0) - g(-1.0)).abs() < 1e-9);

        assert!((clamped(4.0) - g(4.0)).abs() < 1e-9);

        // Natural: extrapolates, but with the end cubic that has the wrong curvature.
        let natural = cubic_spline_with(&data, SplineEnd::Natural).unwrap();

        assert!(natural(-1.0).is_finite());

        assert!((natural(-1.0) - g(-1.0)).abs() > 1e-3);
    }

    #[test]
    fn interp_rejects_degenerate_input() {
        assert!(linear_interpolation(&[(0.0, 0.0)]).is_err());

        assert!(linear_interpolation(&[(0.0, 0.0), (0.0, 1.0)]).is_err());

        assert!(pchip_interpolation(&[]).is_err());

        assert!(
            cubic_spline_with(&[(0.0, 0.0), (1.0, 1.0), (0.5, 2.0)], SplineEnd::Natural).is_err()
        );

        assert!(
            cubic_spline_with(&[(0.0, 0.0), (1.0, 1.0), (2.0, 2.0)], SplineEnd::NotAKnot).is_err()
        );
    }

    #[test]
    fn interp_legacy_spline_is_bit_identical() {
        // Regression guard: the natural spline of the legacy entry point must not move.
        let spline = cubic_spline_interpolation(&[(0.0, 0.0), (1.0, 1.0), (2.0, 0.0)]).unwrap();

        assert_eq!(spline(0.5).to_bits(), 0x3fe6_0000_0000_0000);

        assert_eq!(spline(0.25).to_bits(), 0x3fd7_8000_0000_0000);

        assert_eq!(spline(1.75).to_bits(), 0x3fd7_8000_0000_0000);
    }
}
