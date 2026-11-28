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
//! The one-dimensional function interpolants extrapolate beyond their nodes. The
//! rectilinear-grid routines require each query coordinate to lie inside its axis.
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
/// Returns an error for repeated coordinates, nonfinite samples, or overflowing
/// divided differences or power-basis coefficients. Empty input returns zero.
/// For nonempty input the polynomial retains one coefficient per sample.
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
    let coefficients =
        super::polynomial::lagrange_coefficients(|i| points[i].0, points.iter().map(|&(_, y)| y))?;
    Ok(Polynomial::new(coefficients))
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
        // 4. Find the correct interval in O(log n) via binary search.
        let last = points_owned.len() - 1;
        let i = if x >= points_owned[0].0 && x <= points_owned[last].0 {
            let idx = points_owned.partition_point(|&(px, _)| px <= x);
            idx.saturating_sub(1).min(last - 1)
        } else {
            0
        };

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
    // Empty inputs, mismatched knot counts and ragged control points have
    // no curve: `None` rather than an index or subtraction panic.
    let n = control_points.len().checked_sub(1)?;

    let m = knots.len().checked_sub(1)?;

    if degree > n || Some(m) != n.checked_add(degree)?.checked_add(1) {
        return None;
    }

    let dimension = control_points.first().map_or(0, Vec::len);

    if control_points.iter().any(|p| p.len() != dimension) {
        return None;
    }

    // The knot-span search is a bisection that assumes a sorted knot
    // vector; an unsorted or non-finite one could make it spin.
    if knots.iter().any(|k| !k.is_finite()) || knots.windows(2).any(|w| w[1] < w[0]) {
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

// ============================================================================
// Tensor-product interpolation on rectilinear grids (1D, 2D, 3D, ...)
// ============================================================================

/// Checks the axes and the value array of a rectilinear grid.
/// When supplied, the query must lie inside the grid.
fn check_rectilinear<'grid>(
    axes: impl ExactSizeIterator<Item = &'grid [f64]> + Clone,
    values: &[f64],
    point: Option<&[f64]>,
) -> Result<(), String> {
    if axes.len() == 0 {
        return Err("at least one axis is required".to_string());
    }

    if point.is_some_and(|point| point.len() != axes.len()) {
        return Err("the query point needs one coordinate per axis".to_string());
    }

    if axes.clone().any(|axis| axis.len() < 2) {
        return Err("every axis needs at least two nodes".to_string());
    }

    let size = axes
        .clone()
        .try_fold(1usize, |size, axis| size.checked_mul(axis.len()))
        .ok_or_else(|| "the grid size overflows usize".to_string())?;
    if values.len() != size {
        return Err("the number of values must equal the product of the axis lengths".to_string());
    }

    for (index, axis) in axes.enumerate() {
        if !axis.is_sorted_by(|left, right| left < right)
            || axis.iter().any(|node| !node.is_finite())
        {
            return Err("axes must be finite and strictly increasing".to_string());
        }

        if let Some(&coordinate) = point.and_then(|point| point.get(index)) {
            check_grid_coordinate(axis, coordinate)?;
        }
    }

    // Scan whole blocks so the finite checks can vectorize while an invalid
    // value still short-circuits the rest of the grid.
    if values.chunks(256).any(|block| {
        block
            .iter()
            .fold(false, |invalid, &value| invalid | !value.is_finite())
    }) {
        return Err("the values must be finite".to_string());
    }

    Ok(())
}

fn check_grid_coordinate(
    axis: &[f64],
    coordinate: f64,
) -> Result<(), String> {
    if coordinate.is_finite()
        && axis
            .first()
            .zip(axis.last())
            .is_some_and(|(&first, &last)| first <= coordinate && coordinate <= last)
    {
        Ok(())
    } else {
        Err("the query point lies outside the grid".to_string())
    }
}

#[derive(Clone, Copy)]
struct MultilinearCell {
    stride: usize,
    fraction: f64,
    starts_at_upper: bool,
}

/// A validated rectilinear grid for repeated multilinear interpolation.
///
/// The axes and the array of values are borrowed. Only per-axis metadata is
/// allocated on construction. Each query locates an interval on each axis
/// and interpolates its corners, using fixed-size stack storage.
/// Values are arranged with the first axis varying fastest.
#[derive(Debug)]
pub struct MultilinearGrid<'grid> {
    axes: Vec<&'grid [f64]>,
    values: &'grid [f64],
}

impl<'grid> MultilinearGrid<'grid> {
    /// Prepares a shared grid for calls to [`Self::evaluate`].
    ///
    /// # Errors
    /// Returns an error if there are no axes, any axis has fewer than two
    /// finite strictly increasing nodes, the size overflows, or the values
    /// are non-finite or have the wrong length.
    pub fn new(
        axes: impl IntoIterator<Item = &'grid [f64]>,
        values: &'grid [f64],
    ) -> Result<Self, String> {
        // At least two nodes on every axis makes this many axes impossible
        // for a grid addressable by usize, even before considering its values.
        let axes: Vec<_> = axes.into_iter().take(usize::BITS as usize).collect();
        check_rectilinear(axes.iter().copied(), values, None)?;
        Ok(Self { axes, values })
    }

    /// Interpolates an in-bounds point without allocating or rescanning the grid.
    ///
    /// # Errors
    /// Returns an error if the point has the wrong dimension or any coordinate
    /// is non-finite or outside its axis.
    pub fn evaluate(
        &self,
        point: &[f64],
    ) -> Result<f64, String> {
        evaluate_multilinear::<true>(self.axes.iter().copied(), self.values, point)
    }
}

/// Owned grid data for a session-local prepared interpolation handle.
#[derive(Debug)]
pub(crate) struct OwnedMultilinearGrid {
    axes: Vec<Vec<f64>>,
    values: Vec<f64>,
}

impl OwnedMultilinearGrid {
    pub(crate) fn new(
        axes: Vec<Vec<f64>>,
        values: Vec<f64>,
    ) -> Result<Self, String> {
        check_rectilinear(axes.iter().map(Vec::as_slice), &values, None)?;
        Ok(Self { axes, values })
    }

    pub(crate) fn evaluate(
        &self,
        point: &[f64],
    ) -> Result<f64, String> {
        evaluate_multilinear::<true>(self.axes.iter().map(Vec::as_slice), &self.values, point)
    }

    pub(crate) fn allocated_bytes(&self) -> usize {
        let entries = self
            .axes
            .iter()
            .map(Vec::capacity)
            .fold(self.values.capacity(), usize::saturating_add);
        size_of::<Self>()
            .saturating_add(self.axes.capacity().saturating_mul(size_of::<Vec<f64>>()))
            .saturating_add(entries.saturating_mul(size_of::<f64>()))
    }
}

/// Returns the scaled distances from each endpoint and the scaled width of
/// a finite interval whose width would otherwise overflow.
#[cold]
#[allow(
    clippy::arithmetic_side_effects,
    reason = "halving before subtraction keeps the distances representable"
)]
pub(super) fn scaled_grid_interval_distances(
    lower: f64,
    upper: f64,
    coordinate: f64,
) -> (f64, f64, f64) {
    let lower = lower * 0.5;
    let upper = upper * 0.5;
    let coordinate = coordinate * 0.5;
    (coordinate - lower, upper - coordinate, upper - lower)
}

#[allow(
    clippy::arithmetic_side_effects,
    reason = "validated dimensions bound flat cell offsets and the local coordinates are finite"
)]
fn evaluate_multilinear<'grid, const CHECK_QUERY: bool>(
    axes: impl ExactSizeIterator<Item = &'grid [f64]>,
    values: &[f64],
    point: &[f64],
) -> Result<f64, String> {
    const EMPTY_CELL: MultilinearCell = MultilinearCell {
        stride: 0,
        fraction: 0.0,
        starts_at_upper: false,
    };
    let dimension = axes.len();
    if CHECK_QUERY && point.len() != dimension {
        return Err("the query point needs one coordinate per axis".to_string());
    }
    let mut small_cells = [EMPTY_CELL; 3];
    let mut large_cells;
    let storage: &mut [MultilinearCell] = if dimension <= small_cells.len() {
        &mut small_cells
    } else {
        large_cells = [EMPTY_CELL; usize::BITS as usize];
        &mut large_cells
    };
    let cells = storage
        .get_mut(..dimension)
        .ok_or_else(|| "the grid size overflows usize".to_string())?;
    let mut base = 0usize;
    let mut stride = 1usize;
    for ((axis, &coordinate), cell) in axes.zip(point).zip(&mut *cells) {
        if CHECK_QUERY {
            check_grid_coordinate(axis, coordinate)?;
        }
        let (interval, located) = locate_grid_cell(axis, coordinate, stride)?;
        *cell = located;
        base += interval * stride;
        stride *= axis.len();
    }
    Ok(interpolate_cell(cells, values, base))
}

#[allow(
    clippy::arithmetic_side_effects,
    reason = "the checked grid shape bounds corner offsets and interpolation uses a convex combination"
)]
fn interpolate_cell(
    cells: &[MultilinearCell],
    values: &[f64],
    base: usize,
) -> f64 {
    let Some((cell, lower_axes)) = cells.split_last() else {
        return values.get(base).copied().unwrap_or(f64::NAN);
    };
    sample_grid_axis(cell, base, |index| {
        interpolate_cell(lower_axes, values, index)
    })
}

#[inline]
#[allow(
    clippy::arithmetic_side_effects,
    reason = "the checked shape bounds the upper index and the weights are convex"
)]
fn sample_grid_axis(
    cell: &MultilinearCell,
    base: usize,
    mut read: impl FnMut(usize) -> f64,
) -> f64 {
    let upper = base + cell.stride;
    let (near, far) = if cell.starts_at_upper {
        (upper, base)
    } else {
        (base, upper)
    };
    let anchor = read(near);
    if cell.fraction <= 0.0 {
        return anchor;
    }
    let other = read(far);
    let difference = other - anchor;
    if difference.is_finite() {
        cell.fraction.mul_add(difference, anchor)
    } else {
        (1.0 - cell.fraction).mul_add(anchor, cell.fraction * other)
    }
}

#[inline]
#[allow(
    clippy::arithmetic_side_effects,
    reason = "each checked axis has at least two nodes and the query was validated or clamped"
)]
fn locate_grid_cell(
    axis: &[f64],
    coordinate: f64,
    stride: usize,
) -> Result<(usize, MultilinearCell), String> {
    let interval = axis
        .partition_point(|&node| node <= coordinate)
        .saturating_sub(1)
        .min(axis.len() - 2);
    let &[lower, upper] = axis
        .get(interval..)
        .and_then(|tail| tail.first_chunk::<2>())
        .ok_or_else(|| "the query point lies outside the grid".to_string())?;
    let width = upper - lower;
    let (from_lower, from_upper, width) = if width.is_finite() {
        (coordinate - lower, upper - coordinate, width)
    } else {
        scaled_grid_interval_distances(lower, upper, coordinate)
    };
    let cell = MultilinearCell {
        stride,
        fraction: from_lower.min(from_upper) / width,
        starts_at_upper: from_upper < from_lower,
    };
    Ok((interval, cell))
}

/// Evaluates a validated two-dimensional grid with clamped coordinates.
#[allow(
    clippy::arithmetic_side_effects,
    reason = "checked grid dimensions bound the offsets"
)]
pub(super) fn evaluate_bilinear(
    xs: &[f64],
    ys: &[f64],
    values: &[f64],
    x: f64,
    y: f64,
) -> Result<f64, String> {
    let (column, x_cell) = locate_grid_cell(xs, x, 1)?;
    let (row, y_cell) = locate_grid_cell(ys, y, xs.len())?;
    let base = row * xs.len() + column;
    Ok(sample_grid_axis(&y_cell, base, |row| {
        sample_grid_axis(&x_cell, row, |index| {
            values.get(index).copied().unwrap_or(f64::NAN)
        })
    }))
}

fn multilinear_at<'grid>(
    axes: impl ExactSizeIterator<Item = &'grid [f64]> + Clone,
    values: &[f64],
    point: &[f64],
) -> Result<f64, String> {
    check_rectilinear(axes.clone(), values, Some(point))?;
    evaluate_multilinear::<false>(axes, values, point)
}

/// Multilinear (linear / bilinear / trilinear / ...) interpolation on a
/// rectilinear grid; the dimension is the number of axes.
///
/// `values` is flat with the **first axis fastest**
/// (`values[i0 + n0 * (i1 + n1 * i2)]`). The query point must lie inside the grid
/// (unlike the 1D interpolants of this module, nothing is extrapolated). The result
/// is exact for functions that are affine in every coordinate separately.
/// [`bilinear_interpolation`] and [`trilinear_interpolation`] are the 2D and 3D
/// names of this routine. Use [`MultilinearGrid`] to evaluate repeated queries
/// on the same validated grid.
///
/// # Errors
/// Returns an error for inconsistent sizes, non-increasing or non-finite axes, or a
/// point outside the grid.
pub fn multilinear_interpolation(
    axes: &[Vec<f64>],
    values: &[f64],
    point: &[f64],
) -> Result<f64, String> {
    multilinear_at(axes.iter().map(Vec::as_slice), values, point)
}

/// Bilinear interpolation on the grid `xs x ys` (`values[i + nx * j]`); the 2D case
/// of [`multilinear_interpolation`].
///
/// # Errors
/// As [`multilinear_interpolation`].
pub fn bilinear_interpolation(
    xs: &[f64],
    ys: &[f64],
    values: &[f64],
    x: f64,
    y: f64,
) -> Result<f64, String> {
    multilinear_at([xs, ys].into_iter(), values, &[x, y])
}

/// Trilinear interpolation on the grid `xs x ys x zs` (`values[i + nx * (j + ny * k)]`);
/// the 3D case of [`multilinear_interpolation`].
///
/// # Errors
/// As [`multilinear_interpolation`].
pub fn trilinear_interpolation(
    xs: &[f64],
    ys: &[f64],
    zs: &[f64],
    values: &[f64],
    x: f64,
    y: f64,
    z: f64,
) -> Result<f64, String> {
    multilinear_at([xs, ys, zs].into_iter(), values, &[x, y, z])
}

#[cfg(test)]
mod multilinear_grid_tests {
    #![allow(clippy::arithmetic_side_effects, clippy::indexing_slicing)]

    use super::MultilinearGrid;
    use super::OwnedMultilinearGrid;
    use super::bilinear_interpolation;
    use super::multilinear_interpolation;

    fn close(
        actual: f64,
        expected: f64,
    ) {
        assert!(
            (actual - expected).abs() <= 5e-13 * expected.abs().max(1.0),
            "{actual} != {expected}"
        );
    }

    #[test]
    fn owned_grid_counts_spare_capacity_and_uses_the_same_lookup() -> Result<(), String> {
        let mut x = Vec::with_capacity(12);
        x.extend([0.0, 1.0]);
        let mut y = Vec::with_capacity(6);
        y.extend([0.0, 1.0]);
        let mut values = Vec::with_capacity(20);
        values.extend([1.0, 2.0, 3.0, 4.0]);
        let entries = x.capacity() + y.capacity() + values.capacity();
        let mut axes = Vec::with_capacity(5);
        axes.extend([x, y]);
        let expected = size_of::<OwnedMultilinearGrid>()
            + axes.capacity() * size_of::<Vec<f64>>()
            + entries * size_of::<f64>();
        let grid = OwnedMultilinearGrid::new(axes, values)?;
        assert_eq!(grid.allocated_bytes(), expected);
        close(grid.evaluate(&[0.5, 0.5])?, 2.5);
        Ok(())
    }

    #[test]
    fn prepared_grid_reproduces_a_multiaffine_function_on_a_nonuniform_four_axis_grid()
    -> Result<(), String> {
        let axes = [
            &[-2.0, 1.0, 5.0][..],
            &[-3.0, 4.0][..],
            &[0.0, 0.5, 3.0, 7.5][..],
            &[1.0, 5.0][..],
        ];
        let reference = |x, y, z, w| 3.0 + x * (1.0 + y) - 2.0 * z * w + 0.2 * x * y * z * w;
        let mut values = Vec::new();
        for &w in axes[3] {
            for &z in axes[2] {
                for &y in axes[1] {
                    for &x in axes[0] {
                        values.push(reference(x, y, z, w));
                    }
                }
            }
        }
        let grid = MultilinearGrid::new(axes, &values)?;
        for [x, y, z, w] in [
            [-2.0, 4.0, 7.5, 1.0],
            [5.0, -3.0, 0.0, 5.0],
            [-0.25, 1.5, 0.125, 2.5],
            [3.25, -1.75, 5.5, 4.0],
        ] {
            close(grid.evaluate(&[x, y, z, w])?, reference(x, y, z, w));
        }
        close(
            bilinear_interpolation(
                &[0.0, 2.0, 4.0],
                &[1.0, 3.0],
                &[1.0, 3.0, 5.0, 3.0, 9.0, 15.0],
                3.0,
                2.0,
            )?,
            8.0,
        );
        Ok(())
    }

    #[test]
    fn wide_cells_subnormals_and_near_endpoints_stay_finite_and_accurate() -> Result<(), String> {
        let wide = [-f64::MAX, f64::MAX];
        let grid = MultilinearGrid::new([wide.as_slice()], &[1.0, 3.0])?;
        close(grid.evaluate(&[0.0])?, 2.0);
        close(grid.evaluate(&[-f64::MAX])?, 1.0);
        close(grid.evaluate(&[f64::MAX])?, 3.0);

        let extremes = MultilinearGrid::new([wide.as_slice()], &wide)?;
        close(extremes.evaluate(&[0.0])?, 0.0);
        close(extremes.evaluate(&[f64::MAX * 0.5])?, f64::MAX * 0.5);

        let adjacent = [0.0, f64::from_bits(4)];
        let tiny = MultilinearGrid::new([adjacent.as_slice()], &[-0.0, 8.0])?;
        close(tiny.evaluate(&[f64::from_bits(2)])?, 4.0);
        assert_eq!(tiny.evaluate(&[0.0])?.to_bits(), (-0.0_f64).to_bits());

        let near_endpoint = [-1e16, 1.0];
        let near = MultilinearGrid::new([near_endpoint.as_slice()], &[1e16, 0.0])?;
        close(near.evaluate(&[0.0])?, 1.0);
        Ok(())
    }

    #[test]
    fn grid_validation_rejects_overflow_and_invalid_nodes_values_or_queries() -> Result<(), String>
    {
        let axis = [0.0, 1.0];
        let count = usize::BITS as usize;
        let too_many = vec![axis.to_vec(); count];
        let point = vec![0.5; count];
        let overflow = "the grid size overflows usize";
        assert_eq!(
            multilinear_interpolation(&too_many, &[], &point),
            Err(overflow.to_string())
        );
        assert_eq!(
            MultilinearGrid::new(std::iter::repeat(axis.as_slice()), &[])
                .err()
                .as_deref(),
            Some(overflow)
        );
        assert!(MultilinearGrid::new(std::iter::empty::<&[f64]>(), &[]).is_err());
        assert!(MultilinearGrid::new([&[0.0][..]], &[0.0]).is_err());
        assert!(MultilinearGrid::new([&[0.0, -0.0][..]], &[0.0, 1.0]).is_err());
        assert!(MultilinearGrid::new([&[0.0, f64::INFINITY][..]], &[0.0, 1.0]).is_err());
        assert!(MultilinearGrid::new([axis.as_slice()], &[f64::NAN, 1.0]).is_err());
        assert!(MultilinearGrid::new([axis.as_slice()], &[1.0]).is_err());
        let grid = MultilinearGrid::new([axis.as_slice()], &[1.0, 2.0])?;
        for invalid in [&[][..], &[f64::NAN][..], &[-0.1][..], &[1.1][..]] {
            assert!(grid.evaluate(invalid).is_err());
        }
        Ok(())
    }
}

/// Tensor-product cubic spline interpolation on a rectilinear grid (1D, 2D, 3D,
/// ...) with not-a-knot end conditions along every axis.
///
/// The method is fourth order for smooth data. The first axis is reduced first, then the next, and so on, so the result
/// reproduces any function that is a polynomial of degree three or less in every
/// coordinate. `values` is flat with the first axis fastest; every axis needs at
/// least four nodes; the point must lie inside the grid.
///
/// # Errors
/// As [`multilinear_interpolation`], plus a failure of the 1D spline construction.
pub fn tensor_cubic_spline_interpolation(
    axes: &[Vec<f64>],
    values: &[f64],
    point: &[f64],
) -> Result<f64, String> {
    check_rectilinear(axes.iter().map(Vec::as_slice), values, Some(point))?;
    let dims: Vec<_> = axes.iter().map(Vec::len).collect();

    if dims.iter().any(|&n| n < 4) {
        return Err("the not-a-knot cubic spline needs at least four nodes per axis".to_string());
    }

    let mut current = values.to_vec();

    for a in 0..axes.len() {
        let n = dims[a];

        let reduced = {
            let lines = current.len() / n;
            let build = |line_index: usize| -> Result<f64, String> {
                let line = &current[line_index * n..(line_index + 1) * n];
                let points: Vec<(f64, f64)> =
                    axes[a].iter().copied().zip(line.iter().copied()).collect();

                let spline = cubic_spline_with(&points, SplineEnd::NotAKnot)?;

                Ok(spline(point[a]))
            };
            if lines >= 2 && crate::execution::should_split_work(lines, n.saturating_mul(4).max(1))
            {
                crate::kernels::par::try_map_indexed(lines, 64, build)?
            } else {
                (0..lines).map(build).collect::<Result<Vec<_>, String>>()?
            }
        };

        current = reduced;
    }

    Ok(current[0])
}

/// Dimension-generic interpolation entry point: the order of the tensor-product
/// interpolant used by [`interpolate_nd`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TensorOrder {
    /// Multilinear (second order).
    Linear,
    /// Tensor-product not-a-knot cubic spline (fourth order).
    CubicSpline,
}

/// Interpolates scalar grid data in any number of dimensions (1D, 2D, 3D, ...)
/// with the chosen [`TensorOrder`]; see [`multilinear_interpolation`] and
/// [`tensor_cubic_spline_interpolation`].
///
/// # Errors
/// As the underlying routine.
pub fn interpolate_nd(
    axes: &[Vec<f64>],
    values: &[f64],
    point: &[f64],
    order: TensorOrder,
) -> Result<f64, String> {
    match order {
        | TensorOrder::Linear => multilinear_interpolation(axes, values, point),
        | TensorOrder::CubicSpline => tensor_cubic_spline_interpolation(axes, values, point),
    }
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

    /// The per-line tensor spline reduction must be bit-identical for every
    /// thread count; the grid is wide enough to cross the parallel threshold.
    #[test]
    fn tensor_cubic_spline_is_thread_count_independent() {
        use crate::execution::with_scoped_pool;

        let nx = 32usize;
        let ny = 512usize;
        let ax: Vec<f64> = (0..nx).map(|i| i as f64 / (nx - 1) as f64).collect();
        let ay: Vec<f64> = (0..ny).map(|j| j as f64 / (ny - 1) as f64).collect();
        let f = |x: f64, y: f64| (2.0 * x).sin() * (3.0 * y).cos() + 0.25 * x * y;
        let mut values = vec![0.0; nx * ny];
        for j in 0..ny {
            for i in 0..nx {
                values[i + j * nx] = f(ax[i], ay[j]);
            }
        }
        let axes = [ax, ay];
        let point = [0.37, 0.61];

        let serial = with_scoped_pool(Some(1), || {
            tensor_cubic_spline_interpolation(&axes, &values, &point).unwrap()
        });
        for threads in [2usize, 4] {
            let parallel = with_scoped_pool(Some(threads), || {
                tensor_cubic_spline_interpolation(&axes, &values, &point).unwrap()
            });
            assert_eq!(parallel.to_bits(), serial.to_bits());
        }
    }
}
