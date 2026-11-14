//! Interpolation and approximation kernels.
//!
//! Splines (natural, clamped,
//! not-a-knot, PCHIP, Akima), B-splines, bicubic / N-D grid interpolation,
//! scattered-data interpolation (RBF, Shepard IDW, Delaunay), and Pade / rational.

use crate::kernels::dense::DenseError;
use crate::kernels::dense::Mat;
use crate::kernels::dense::{
    self,
};
use crate::kernels::geometry;

/// Error returned by extended interpolation kernels.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum InterpError {
    /// Input was empty.
    Empty,
    /// Dimensions or lengths did not match.
    ShapeMismatch(String),
    /// Invalid domain or non-monotonic knots/axes.
    Domain(String),
    /// Linear system was singular.
    Singular,
}

impl std::fmt::Display for InterpError {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        match self {
            | Self::Empty => f.write_str("interpolation input is empty"),
            | Self::ShapeMismatch(s) => write!(f, "shape mismatch: {s}"),
            | Self::Domain(s) => write!(f, "domain error: {s}"),
            | Self::Singular => f.write_str("singular system in interpolation"),
        }
    }
}

impl std::error::Error for InterpError {}

fn solve_flat(
    n: usize,
    a: &[f64],
    b: &[f64],
) -> Result<Vec<f64>, InterpError> {
    let mat = Mat {
        rows: n,
        cols: n,
        data: a.to_vec(),
    };
    dense::solve(&mat, b).map_err(|e| {
        match e {
            | DenseError::Singular => InterpError::Singular,
            | _ => InterpError::Domain(format!("dense solve failed: {e:?}")),
        }
    })
}

/// Boundary condition for [`cubic_spline`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SplineBc {
    /// Natural spline: `y''(x_0) = y''(x_{n-1}) = 0`.
    Natural,
    /// Clamped spline: prescribed first derivatives `(dy0, dyn)` at endpoints.
    Clamped(f64, f64),
    /// Not-a-knot spline: `y'''` is continuous across `x_1` and `x_{n-2}`.
    NotAKnot,
}

fn validate_1d(
    xs: &[f64],
    ys: &[f64],
) -> Result<(), InterpError> {
    if xs.len() != ys.len() {
        return Err(InterpError::ShapeMismatch(format!(
            "interpolation needs two samples of equal length, found {} and {}",
            xs.len(),
            ys.len()
        )));
    }
    if xs.is_empty() {
        return Err(InterpError::Empty);
    }
    if xs.iter().chain(ys).any(|v| !v.is_finite()) {
        return Err(InterpError::Domain(
            "interpolation data must be finite".into(),
        ));
    }
    for i in 1..xs.len() {
        if xs[i] <= xs[i - 1] {
            return Err(InterpError::Domain(
                "knots must be strictly increasing".to_string(),
            ));
        }
    }
    Ok(())
}

fn find_segment(
    xs: &[f64],
    x: f64,
) -> usize {
    let n = xs.len();
    if x <= xs[0] {
        return 0;
    }
    if x >= xs[n - 1] {
        return n - 2;
    }
    let mut lo = 0;
    let mut hi = n - 1;
    while hi - lo > 1 {
        let mid = usize::midpoint(lo, hi);
        if xs[mid] > x {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    lo
}

fn eval_hermite(
    x0: f64,
    x1: f64,
    y0: f64,
    y1: f64,
    d0: f64,
    d1: f64,
    x: f64,
) -> f64 {
    let h = x1 - x0;
    let t = (x - x0) / h;
    let t2 = t * t;
    let t3 = t2 * t;
    let h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
    let h10 = t3 - 2.0 * t2 + t;
    let h01 = -2.0 * t3 + 3.0 * t2;
    let h11 = t3 - t2;
    h00 * y0 + h10 * h * d0 + h01 * y1 + h11 * h * d1
}

/// Solves a tridiagonal system `a[i]*x[i-1] + b[i]*x[i] + c[i]*x[i+1] = d[i]`.
fn solve_tridiagonal(
    a: &[f64],
    b: &[f64],
    c: &[f64],
    d: &[f64],
) -> Result<Vec<f64>, InterpError> {
    let n = b.len();
    if n == 0 {
        return Ok(Vec::new());
    }
    let mut cp = vec![0.0; n];
    let mut dp = vec![0.0; n];
    if b[0] == 0.0 {
        return Err(InterpError::Singular);
    }
    cp[0] = c[0] / b[0];
    dp[0] = d[0] / b[0];
    for i in 1..n {
        let denom = b[i] - a[i] * cp[i - 1];
        if denom == 0.0 {
            return Err(InterpError::Singular);
        }
        cp[i] = if i + 1 < n {
            c[i] / denom
        } else {
            0.0
        };
        dp[i] = (d[i] - a[i] * dp[i - 1]) / denom;
    }
    let mut x = vec![0.0; n];
    x[n - 1] = dp[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = dp[i] - cp[i] * x[i + 1];
    }
    Ok(x)
}

/// Cubic spline interpolation with selectable boundary conditions
/// (`Natural`, `Clamped`, or `NotAKnot`).
///
/// # Errors
/// [`InterpError::ShapeMismatch`], [`InterpError::Empty`], or [`InterpError::Domain`]
/// when knots are not strictly increasing.
pub fn cubic_spline(
    xs: &[f64],
    ys: &[f64],
    x: f64,
    bc: SplineBc,
) -> Result<f64, InterpError> {
    if !x.is_finite() {
        return Err(InterpError::Domain("query must be finite".into()));
    }
    validate_1d(xs, ys)?;
    let n = xs.len();
    if n == 1 {
        return Ok(ys[0]);
    }
    if n == 2 {
        return match bc {
            | SplineBc::Clamped(d0, d1) => Ok(eval_hermite(xs[0], xs[1], ys[0], ys[1], d0, d1, x)),
            | _ => interp1(xs, ys, x),
        };
    }

    // Solve for first derivatives d[0..n] at the knots.
    let mut h = vec![0.0; n - 1];
    let mut delta = vec![0.0; n - 1];
    for i in 0..n - 1 {
        h[i] = xs[i + 1] - xs[i];
        delta[i] = (ys[i + 1] - ys[i]) / h[i];
    }

    let mut a = vec![0.0; n];
    let mut b = vec![0.0; n];
    let mut c = vec![0.0; n];
    let mut rhs = vec![0.0; n];

    for i in 1..n - 1 {
        a[i] = h[i];
        b[i] = 2.0 * (h[i - 1] + h[i]);
        c[i] = h[i - 1];
        rhs[i] = 3.0 * (h[i] * delta[i - 1] + h[i - 1] * delta[i]);
    }

    match bc {
        | SplineBc::Natural => {
            b[0] = 2.0;
            c[0] = 1.0;
            rhs[0] = 3.0 * delta[0];
            a[n - 1] = 1.0;
            b[n - 1] = 2.0;
            rhs[n - 1] = 3.0 * delta[n - 2];
        },
        | SplineBc::Clamped(d0, dn) => {
            b[0] = 1.0;
            c[0] = 0.0;
            rhs[0] = d0;
            a[n - 1] = 0.0;
            b[n - 1] = 1.0;
            rhs[n - 1] = dn;
        },
        | SplineBc::NotAKnot => {
            if n == 3 {
                // With 3 points, not-a-knot is the unique quadratic through all 3 points.
                let l0 = ((x - xs[1]) * (x - xs[2])) / ((xs[0] - xs[1]) * (xs[0] - xs[2]));
                let l1 = ((x - xs[0]) * (x - xs[2])) / ((xs[1] - xs[0]) * (xs[1] - xs[2]));
                let l2 = ((x - xs[0]) * (x - xs[1])) / ((xs[2] - xs[0]) * (xs[2] - xs[1]));
                return Ok(ys[0] * l0 + ys[1] * l1 + ys[2] * l2);
            }
            // Eliminate d[2] from row 0 using row 1:
            // h0^2 d0 + (h0 + h1)^2 d1 - ... -> standard not-a-knot condition:
            // h1^2 d0 + (h0 + h1) h1 d1 ... wait: using third derivative continuity at x_1:
            // y'''_0(x_1) = y'''_1(x_1) => (d0 + d1 - 2 delta0)/h0^2 = (d1 + d2 - 2 delta1)/h1^2
            // Combine with row 1: h1 d0 + 2(h0 + h1) d1 + h0 d2 = 3(h1 delta0 + h0 delta1)
            // Multiply row 1 by h0/h1^2? Easier: use dense or direct elimination:
            let r0 = h[1] * h[1];
            let r1 = h[0] * h[0];
            // r0 * (d0 + d1 - 2 delta0) = r1 * (d1 + d2 - 2 delta1)
            // => r0 d0 + (r0 - r1) d1 - r1 d2 = 2 r0 delta0 - 2 r1 delta1
            // Row 1: h1 d0 + 2(h0 + h1) d1 + h0 d2 = 3(h1 delta0 + h0 delta1)
            // Eliminate d2: multiply eq0 by h0, eq1 by r1 (= h0^2), add:
            // (h0 r0 + h1 r1) d0 + (h0(r0 - r1) + 2(h0+h1) r1) d1 = h0(2 r0 delta0 - 2 r1 delta1) + r1 * rhs[1]
            b[0] = h[0] * r0 + h[1] * r1;
            c[0] = h[0] * (r0 - r1) + 2.0 * (h[0] + h[1]) * r1;
            rhs[0] = h[0] * (2.0 * r0 * delta[0] - 2.0 * r1 * delta[1]) + r1 * rhs[1];

            // At x_{n-2}: (d_{n-3} + d_{n-2} - 2 delta_{n-3})/h_{n-3}^2 = (d_{n-2} + d_{n-1} - 2 delta_{n-2})/h_{n-2}^2
            let hn3 = h[n - 3];
            let hn2 = h[n - 2];
            let s0 = hn2 * hn2;
            let s1 = hn3 * hn3;
            // s0 d_{n-3} + (s0 - s1) d_{n-2} - s1 d_{n-1} = 2 s0 delta_{n-3} - 2 s1 delta_{n-2}
            // Row n-2: hn2 d_{n-3} + 2(hn3 + hn2) d_{n-2} + hn3 d_{n-1} = rhs[n-2]
            // Eliminate d_{n-3}: multiply eq_last by hn2, row n-2 by s0, subtract:
            a[n - 1] = hn2 * (s0 - s1) - 2.0 * (hn3 + hn2) * s0;
            b[n - 1] = -hn2 * s1 - hn3 * s0;
            rhs[n - 1] =
                hn2 * (2.0 * s0 * delta[n - 3] - 2.0 * s1 * delta[n - 2]) - s0 * rhs[n - 2];
        },
    }

    let d = solve_tridiagonal(&a, &b, &c, &rhs)?;
    let idx = find_segment(xs, x);
    Ok(eval_hermite(
        xs[idx],
        xs[idx + 1],
        ys[idx],
        ys[idx + 1],
        d[idx],
        d[idx + 1],
        x,
    ))
}

/// Piecewise Cubic Hermite Interpolating Polynomial (PCHIP, Fritsch–Carlson).
///
/// Preserves monotonicity of the data without overshoot.
///
/// # Errors
/// [`InterpError::ShapeMismatch`], [`InterpError::Empty`], or [`InterpError::Domain`].
pub fn pchip(
    xs: &[f64],
    ys: &[f64],
    x: f64,
) -> Result<f64, InterpError> {
    if !x.is_finite() {
        return Err(InterpError::Domain("query must be finite".into()));
    }
    validate_1d(xs, ys)?;
    let n = xs.len();
    if n == 1 {
        return Ok(ys[0]);
    }
    if n == 2 {
        return interp1(xs, ys, x);
    }

    let mut h = vec![0.0; n - 1];
    let mut delta = vec![0.0; n - 1];
    for i in 0..n - 1 {
        h[i] = xs[i + 1] - xs[i];
        delta[i] = (ys[i + 1] - ys[i]) / h[i];
    }

    let mut d = vec![0.0; n];
    for i in 1..n - 1 {
        if delta[i - 1] * delta[i] > 0.0 {
            let w1 = 2.0 * h[i] + h[i - 1];
            let w2 = h[i] + 2.0 * h[i - 1];
            d[i] = (w1 + w2) / (w1 / delta[i - 1] + w2 / delta[i]);
        } else {
            d[i] = 0.0;
        }
    }

    // Endpoint slopes (one-sided three-point formula with shape preservation).
    let endpoint_slope = |h0: f64, h1: f64, del0: f64, del1: f64| -> f64 {
        let mut slope = ((2.0 * h0 + h1) * del0 - h0 * del1) / (h0 + h1);
        if slope.signum() != del0.signum() {
            slope = 0.0;
        } else if del0.signum() != del1.signum() && slope.abs() > 3.0 * del0.abs() {
            slope = 3.0 * del0;
        }
        slope
    };
    d[0] = endpoint_slope(h[0], h[1], delta[0], delta[1]);
    d[n - 1] = endpoint_slope(h[n - 2], h[n - 3], delta[n - 2], delta[n - 3]);

    let idx = find_segment(xs, x);
    Ok(eval_hermite(
        xs[idx],
        xs[idx + 1],
        ys[idx],
        ys[idx + 1],
        d[idx],
        d[idx + 1],
        x,
    ))
}

/// Akima 1D cubic spline interpolation.
///
/// Reduces oscillation near sharp corners compared to global cubic splines.
///
/// # Errors
/// [`InterpError::ShapeMismatch`], [`InterpError::Empty`], or [`InterpError::Domain`].
pub fn akima(
    xs: &[f64],
    ys: &[f64],
    x: f64,
) -> Result<f64, InterpError> {
    if !x.is_finite() {
        return Err(InterpError::Domain("query must be finite".into()));
    }
    validate_1d(xs, ys)?;
    let n = xs.len();
    if n < 3 {
        return interp1(xs, ys, x);
    }

    // Secant slopes m[2..n+1] correspond to segments 0..n-1, padded by 2 on each side.
    let mut m = vec![0.0; n + 3];
    for i in 0..n - 1 {
        m[i + 2] = (ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i]);
    }
    m[1] = 2.0 * m[2] - m[3];
    m[0] = 2.0 * m[1] - m[2];
    m[n + 1] = 2.0 * m[n] - m[n - 1];
    m[n + 2] = 2.0 * m[n + 1] - m[n];

    let mut d = vec![0.0; n];
    for i in 0..n {
        let w1 = (m[i + 3] - m[i + 2]).abs();
        let w2 = (m[i + 1] - m[i]).abs();
        let denom = w1 + w2;
        d[i] = if denom == 0.0 {
            f64::midpoint(m[i + 1], m[i + 2])
        } else {
            (w1 * m[i + 1] + w2 * m[i + 2]) / denom
        };
    }

    let idx = find_segment(xs, x);
    Ok(eval_hermite(
        xs[idx],
        xs[idx + 1],
        ys[idx],
        ys[idx + 1],
        d[idx],
        d[idx + 1],
        x,
    ))
}

/// Evaluates a B-spline of given `degree` at `x` using the Cox–de Boor recursion.
///
/// Requires `knots.len() == coeffs.len() + degree + 1` and non-decreasing `knots`.
///
/// # Errors
/// [`InterpError::ShapeMismatch`] or [`InterpError::Domain`].
pub fn bspline_eval(
    knots: &[f64],
    coeffs: &[f64],
    degree: usize,
    x: f64,
) -> Result<f64, InterpError> {
    if !x.is_finite() {
        return Err(InterpError::Domain("query must be finite".into()));
    }
    if coeffs.is_empty() {
        return Err(InterpError::Empty);
    }
    if knots.len() != coeffs.len() + degree + 1 {
        return Err(InterpError::ShapeMismatch(format!(
            "B-spline requires knots.len() == coeffs.len() + degree + 1 ({}), found {}",
            coeffs.len() + degree + 1,
            knots.len()
        )));
    }
    for i in 1..knots.len() {
        if knots[i] < knots[i - 1] {
            return Err(InterpError::Domain(
                "B-spline knots must be non-decreasing".to_string(),
            ));
        }
    }
    let n = coeffs.len() - 1;
    let x_clamped = x.clamp(knots[degree], knots[n + 1]);

    // Find knot span k such that knots[k] <= x_clamped < knots[k+1]
    let k = if x_clamped >= knots[n + 1] {
        let mut idx = n;
        while idx > degree && crate::kernels::exact_eq(knots[idx], knots[n + 1]) {
            idx -= 1;
        }
        idx
    } else {
        let mut lo = degree;
        let mut hi = n + 1;
        while hi - lo > 1 {
            let mid = usize::midpoint(lo, hi);
            if x_clamped < knots[mid] {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        lo
    };

    // de Boor's algorithm
    let mut d: Vec<f64> = (0..=degree).map(|j| coeffs[k - degree + j]).collect();
    for r in 1..=degree {
        for j in (r..=degree).rev() {
            let i = k - degree + j;
            let denom = knots[i + degree + 1 - r] - knots[i];
            let alpha = if denom == 0.0 {
                0.0
            } else {
                (x_clamped - knots[i]) / denom
            };
            d[j] = (1.0 - alpha) * d[j - 1] + alpha * d[j];
        }
    }
    Ok(d[degree])
}

/// Fitted B-spline interpolant.
#[derive(Clone, Debug)]
pub struct BSpline {
    /// Clamped knot vector.
    pub knots: Vec<f64>,
    /// Control-point coefficients.
    pub coeffs: Vec<f64>,
    /// Spline degree.
    pub degree: usize,
}

impl BSpline {
    /// Fits an interpolating B-spline of `degree` through `(xs, ys)`.
    ///
    /// # Errors
    /// [`InterpError::ShapeMismatch`], [`InterpError::Domain`], or [`InterpError::Singular`].
    pub fn interpolate(
        xs: &[f64],
        ys: &[f64],
        degree: usize,
    ) -> Result<Self, InterpError> {
        validate_1d(xs, ys)?;
        let n = xs.len();
        if degree >= n {
            return Err(InterpError::Domain(format!(
                "B-spline degree ({degree}) must be < number of points ({n})"
            )));
        }
        // Clamped knot vector with averaging of interior data points.
        let m = n + degree + 1;
        let mut knots = vec![0.0; m];
        for i in 0..=degree {
            knots[i] = xs[0];
            knots[m - 1 - i] = xs[n - 1];
        }
        for j in 1..(n - degree) {
            let mut sum = 0.0;
            for i in j..(j + degree) {
                sum += xs[i];
            }
            knots[j + degree] = sum / degree as f64;
        }

        // Assemble collocation matrix N_{j,p}(xs[i])
        let mut mat = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                // Evaluate j-th basis function at xs[i]
                let mut basis = vec![0.0; n];
                basis[j] = 1.0;
                mat[i * n + j] = bspline_eval(&knots, &basis, degree, xs[i])?;
            }
        }
        let coeffs = solve_flat(n, &mat, ys)?;
        Ok(Self {
            knots,
            coeffs,
            degree,
        })
    }

    /// Evaluates the B-spline at `x`.
    ///
    /// # Errors
    /// [`InterpError::Domain`] when the query is nonfinite, or an error from
    /// [`bspline_eval`] if the stored spline representation is invalid.
    pub fn eval(
        &self,
        x: f64,
    ) -> Result<f64, InterpError> {
        bspline_eval(&self.knots, &self.coeffs, self.degree, x)
    }
}

/// Bicubic interpolation on a rectilinear 2D grid.
///
/// `grid` is row-major of length `ys.len() * xs.len()` where `grid[j * xs.len() + i]`
/// is the value at `(xs[i], ys[j])`.
///
/// # Errors
/// [`InterpError::ShapeMismatch`], [`InterpError::Empty`], or [`InterpError::Domain`].
pub fn bicubic(
    xs: &[f64],
    ys: &[f64],
    grid: &[f64],
    x: f64,
    y: f64,
) -> Result<f64, InterpError> {
    if xs.is_empty() || ys.is_empty() {
        return Err(InterpError::Empty);
    }
    if grid.len() != xs.len() * ys.len() {
        return Err(InterpError::ShapeMismatch(format!(
            "bicubic grid needs {} * {} = {} entries, found {}",
            ys.len(),
            xs.len(),
            ys.len() * xs.len(),
            grid.len()
        )));
    }
    let nx = xs.len();
    let ny = ys.len();
    let mut row_vals = Vec::with_capacity(ny);
    for j in 0..ny {
        let row = &grid[j * nx..(j + 1) * nx];
        row_vals.push(cubic_spline(xs, row, x, SplineBc::NotAKnot)?);
    }
    cubic_spline(ys, &row_vals, y, SplineBc::NotAKnot)
}

/// N-dimensional multilinear interpolation on a rectilinear grid.
///
/// `axes[d]` gives the strictly increasing coordinates along axis `d`,
/// `grid` is row-major of length `prod(axes[d].len())`, and `point` has length `axes.len()`.
///
/// # Errors
/// [`InterpError::ShapeMismatch`], [`InterpError::Empty`], or [`InterpError::Domain`].
pub fn interpn(
    axes: &[&[f64]],
    grid: &[f64],
    point: &[f64],
) -> Result<f64, InterpError> {
    let ndim = axes.len();
    if ndim > 20 {
        return Err(InterpError::Domain(
            "multilinear interpolation supports at most 20 axes".into(),
        ));
    }
    if point.iter().chain(grid).any(|v| !v.is_finite()) {
        return Err(InterpError::Domain("grid and query must be finite".into()));
    }
    if ndim == 0 {
        return Err(InterpError::Empty);
    }
    if point.len() != ndim {
        return Err(InterpError::ShapeMismatch(format!(
            "point has {} coordinates for {} axes",
            point.len(),
            ndim
        )));
    }
    let mut total_len = 1usize;
    for ax in axes {
        if ax.is_empty() {
            return Err(InterpError::Empty);
        }
        if ax.iter().any(|v| !v.is_finite()) {
            return Err(InterpError::Domain("grid axes must be finite".into()));
        }
        for i in 1..ax.len() {
            if ax[i] <= ax[i - 1] {
                return Err(InterpError::Domain(
                    "grid axes must be strictly increasing".to_string(),
                ));
            }
        }
        total_len = total_len
            .checked_mul(ax.len())
            .ok_or_else(|| InterpError::Domain("grid shape overflows usize".into()))?;
    }
    if grid.len() != total_len {
        return Err(InterpError::ShapeMismatch(format!(
            "grid has {} entries, expected {total_len}",
            grid.len()
        )));
    }

    let mut strides = vec![1usize; ndim];
    for d in (0..ndim - 1).rev() {
        strides[d] = strides[d + 1] * axes[d + 1].len();
    }

    let mut lo_idx = Vec::with_capacity(ndim);
    let mut hi_idx = Vec::with_capacity(ndim);
    let mut frac = Vec::with_capacity(ndim);
    for d in 0..ndim {
        let ax = axes[d];
        if ax.len() == 1 {
            lo_idx.push(0);
            hi_idx.push(0);
            frac.push(0.0);
        } else {
            let seg = find_segment(ax, point[d]);
            lo_idx.push(seg);
            hi_idx.push(seg + 1);
            frac.push((point[d] - ax[seg]) / (ax[seg + 1] - ax[seg]));
        }
    }

    let corners = 1usize << ndim;
    let mut acc = 0.0;
    for mask in 0..corners {
        let mut flat = 0usize;
        let mut weight = 1.0;
        for d in 0..ndim {
            if (mask >> d) & 1 == 1 {
                flat += hi_idx[d] * strides[d];
                weight *= frac[d];
            } else {
                flat += lo_idx[d] * strides[d];
                weight *= 1.0 - frac[d];
            }
        }
        acc += weight * grid[flat];
    }
    Ok(acc)
}

/// Radial basis function kernel for [`rbf_interp`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum RbfKernel {
    /// `exp(-(r / eps)^2)`
    Gaussian(f64),
    /// `sqrt(1 + (r / eps)^2)`
    Multiquadric(f64),
    /// `1 / sqrt(1 + (r / eps)^2)`
    InverseMultiquadric(f64),
    /// `r^2 ln(r)` (`0` at `r = 0`)
    ThinPlateSpline,
    /// `r^3`
    Cubic,
    /// `r`
    Linear,
}

impl RbfKernel {
    fn eval(
        self,
        r: f64,
    ) -> f64 {
        match self {
            | Self::Gaussian(eps) => {
                let s = r / eps;
                (-s * s).exp()
            },
            | Self::Multiquadric(eps) => {
                let s = r / eps;
                (1.0 + s * s).sqrt()
            },
            | Self::InverseMultiquadric(eps) => {
                let s = r / eps;
                1.0 / (1.0 + s * s).sqrt()
            },
            | Self::ThinPlateSpline => {
                if r == 0.0 {
                    0.0
                } else {
                    r * r * r.ln()
                }
            },
            | Self::Cubic => r * r * r,
            | Self::Linear => r,
        }
    }
}

fn euclid_dist(
    a: &[f64],
    b: &[f64],
) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| (x - y) * (x - y))
        .sum::<f64>()
        .sqrt()
}

/// Radial basis function (RBF) scattered-data interpolation in `R^d`.
///
/// # Errors
/// [`InterpError::Empty`], [`InterpError::ShapeMismatch`], or [`InterpError::Singular`].
pub fn rbf_interp(
    points: &[&[f64]],
    values: &[f64],
    query: &[f64],
    kernel: RbfKernel,
) -> Result<f64, InterpError> {
    if points.is_empty() {
        return Err(InterpError::Empty);
    }
    if points.len() != values.len() {
        return Err(InterpError::ShapeMismatch(format!(
            "rbf_interp got {} points and {} values",
            points.len(),
            values.len()
        )));
    }
    let dim = query.len();
    for p in points {
        if p.len() != dim {
            return Err(InterpError::ShapeMismatch(format!(
                "point dimension {} does not match query dimension {dim}",
                p.len()
            )));
        }
    }
    let n = points.len();
    let mut phi = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            phi[i * n + j] = kernel.eval(euclid_dist(points[i], points[j]));
        }
    }
    let weights = solve_flat(n, &phi, values)?;
    let mut sum = 0.0;
    for i in 0..n {
        sum += weights[i] * kernel.eval(euclid_dist(query, points[i]));
    }
    Ok(sum)
}

/// Shepard inverse-distance weighting (IDW) scattered-data interpolation in `R^d`.
///
/// # Errors
/// [`InterpError::Empty`], [`InterpError::ShapeMismatch`], or [`InterpError::Domain`]
/// when `power <= 0`.
pub fn shepard_idw(
    points: &[&[f64]],
    values: &[f64],
    query: &[f64],
    power: f64,
) -> Result<f64, InterpError> {
    if points.is_empty() {
        return Err(InterpError::Empty);
    }
    if points.len() != values.len() {
        return Err(InterpError::ShapeMismatch(format!(
            "shepard_idw got {} points and {} values",
            points.len(),
            values.len()
        )));
    }
    if !power.is_finite() || power <= 0.0 {
        return Err(InterpError::Domain(
            "IDW power must be positive".to_string(),
        ));
    }
    let dim = query.len();
    let mut num = 0.0;
    let mut den = 0.0;
    for (p, &v) in points.iter().zip(values.iter()) {
        if p.len() != dim {
            return Err(InterpError::ShapeMismatch(format!(
                "point dimension {} does not match query dimension {dim}",
                p.len()
            )));
        }
        let d = euclid_dist(p, query);
        if d == 0.0 {
            return Ok(v);
        }
        let w = d.powf(-power);
        num += w * v;
        den += w;
    }
    Ok(num / den)
}

/// 2D scattered-data linear interpolation over the Delaunay triangulation of `points`.
///
/// Uses barycentric coordinates inside the containing Delaunay triangle; if `query` lies
/// outside the convex hull, falls back to Shepard IDW (`power = 2`).
///
/// # Errors
/// [`InterpError::ShapeMismatch`] or [`InterpError::Domain`] when fewer than 3 points are given.
pub fn delaunay_linear_2d(
    points: &[(f64, f64)],
    values: &[f64],
    x: f64,
    y: f64,
) -> Result<f64, InterpError> {
    if points.len() != values.len() {
        return Err(InterpError::ShapeMismatch(format!(
            "delaunay_linear_2d got {} points and {} values",
            points.len(),
            values.len()
        )));
    }
    if points.len() < 3 {
        return Err(InterpError::Domain(
            "delaunay_linear_2d requires at least 3 points".to_string(),
        ));
    }
    let pts_arr: Vec<[f64; 2]> = points.iter().map(|&(px, py)| [px, py]).collect();
    let tris = geometry::delaunay_2d(&pts_arr);
    for &[i, j, k] in &tris {
        let (x0, y0) = points[i];
        let (x1, y1) = points[j];
        let (x2, y2) = points[k];
        let denom = (y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2);
        if denom.abs() < 1e-15 {
            continue;
        }
        let w0 = ((y1 - y2) * (x - x2) + (x2 - x1) * (y - y2)) / denom;
        let w1 = ((y2 - y0) * (x - x2) + (x0 - x2) * (y - y2)) / denom;
        let w2 = 1.0 - w0 - w1;
        if w0 >= -1e-12 && w1 >= -1e-12 && w2 >= -1e-12 {
            return Ok(w0 * values[i] + w1 * values[j] + w2 * values[k]);
        }
    }
    let pts: Vec<[f64; 2]> = points.iter().map(|&(px, py)| [px, py]).collect();
    let refs: Vec<&[f64]> = pts.iter().map(<[f64; 2]>::as_slice).collect();
    shepard_idw(&refs, values, &[x, y], 2.0)
}

/// Computes the `[m / n]` Padé approximant `P(x) / Q(x)` from Taylor series
/// coefficients `coeffs[0..=m+n]` (ascending power order).
///
/// Returns `(p, q)` with `p.len() == m + 1`, `q.len() == n + 1`, and `q[0] == 1.0`.
///
/// # Errors
/// [`InterpError::Domain`] when `coeffs.len() < m + n + 1`, or [`InterpError::Singular`]
/// when the Padé system is singular.
pub fn pade(
    coeffs: &[f64],
    m: usize,
    n: usize,
) -> Result<(Vec<f64>, Vec<f64>), InterpError> {
    if coeffs.len() < m + n + 1 {
        return Err(InterpError::Domain(format!(
            "Padé [{m}/{n}] needs at least {} coefficients, got {}",
            m + n + 1,
            coeffs.len()
        )));
    }
    let mut q = vec![0.0; n + 1];
    q[0] = 1.0;
    if n > 0 {
        let mut mat = vec![0.0; n * n];
        let mut rhs = vec![0.0; n];
        for row in 0..n {
            let k = m + 1 + row;
            rhs[row] = -coeffs[k];
            for col in 0..n {
                let idx = k as isize - (col as isize + 1);
                mat[row * n + col] = if idx >= 0 {
                    coeffs[idx as usize]
                } else {
                    0.0
                };
            }
        }
        let q_tail = solve_flat(n, &mat, &rhs)?;
        q[1..=n].copy_from_slice(&q_tail);
    }
    let mut p = vec![0.0; m + 1];
    for k in 0..=m {
        let mut s = 0.0;
        for j in 0..=k.min(n) {
            s += q[j] * coeffs[k - j];
        }
        p[k] = s;
    }
    Ok((p, q))
}

/// Evaluates a rational function `P(x) / Q(x)` with ascending-power coefficients.
///
/// # Errors
/// [`InterpError::Domain`] when `Q(x) == 0`.
pub fn rational_eval(
    p: &[f64],
    q: &[f64],
    x: f64,
) -> Result<f64, InterpError> {
    let den = polyval(q, x);
    if den == 0.0 {
        return Err(InterpError::Domain(
            "rational denominator is zero".to_string(),
        ));
    }
    Ok(polyval(p, x) / den)
}

/// Linear interpolation of `(xs, ys)` at `x`.
///
/// Outside the range the nearest endpoint is extended linearly from its own
/// segment, which is what a piecewise-linear interpolant does at its ends.
///
/// # Errors
///
/// [`InterpError::ShapeMismatch`] when the two samples differ in length, and
/// [`InterpError::Empty`] when there are no points.
pub fn interp1(
    xs: &[f64],
    ys: &[f64],
    x: f64,
) -> Result<f64, InterpError> {
    validate_1d(xs, ys)?;
    if !x.is_finite() {
        return Err(InterpError::Domain("query must be finite".into()));
    }
    if xs.len() != ys.len() {
        return Err(InterpError::ShapeMismatch(format!(
            "interp1 needs two samples of one length, found {} and {}",
            xs.len(),
            ys.len()
        )));
    }
    if xs.is_empty() {
        return Err(InterpError::Empty);
    }
    if xs.len() == 1 {
        return Ok(ys[0]);
    }
    // Locate the segment; the samples are assumed sorted, and a scan keeps the
    // function total when they are not.
    let mut index = 0;
    for candidate in 0..xs.len() - 1 {
        if x >= xs[candidate] && x <= xs[candidate + 1] {
            index = candidate;
            break;
        }
        if x < xs[0] {
            index = 0;
            break;
        }
        index = candidate;
    }
    let (x0, x1) = (xs[index], xs[index + 1]);
    let (y0, y1) = (ys[index], ys[index + 1]);
    if crate::kernels::exact_eq(x1, x0) {
        return Ok(y0);
    }
    Ok(y0 + (y1 - y0) * (x - x0) / (x1 - x0))
}

/// Natural cubic spline interpolation at `x`.
///
/// # Errors
///
/// [`InterpError::ShapeMismatch`] when the samples differ in length, and
/// [`InterpError::Domain`] when fewer than three points are given (a natural
/// spline needs a second derivative to pin down).
pub fn spline(
    xs: &[f64],
    ys: &[f64],
    x: f64,
) -> Result<f64, InterpError> {
    if xs.len() < 3 {
        return interp1(xs, ys, x);
    }
    cubic_spline(xs, ys, x, SplineBc::Natural)
}

/// Bilinear interpolation on a regular grid.
///
/// `grid` is `ny × nx` in row-major order, covering `xs × ys`.
///
/// # Errors
///
/// [`InterpError::ShapeMismatch`] when the grid does not match the axes.
pub fn interp2(
    xs: &[f64],
    ys: &[f64],
    grid: &[f64],
    x: f64,
    y: f64,
) -> Result<f64, InterpError> {
    validate_1d(xs, xs)?;
    validate_1d(ys, ys)?;
    if !x.is_finite() || !y.is_finite() || grid.iter().any(|v| !v.is_finite()) {
        return Err(InterpError::Domain(
            "grid and queries must be finite".into(),
        ));
    }

    let (nx, ny) = (xs.len(), ys.len());
    if nx.checked_mul(ny) != Some(grid.len()) {
        return Err(InterpError::ShapeMismatch(format!(
            "interp2 needs a {ny}x{nx} grid, found {} values",
            grid.len()
        )));
    }
    if nx < 2 || ny < 2 {
        return Err(InterpError::Domain(
            "interp2 needs at least two points per axis".to_string(),
        ));
    }
    let locate = |axis: &[f64], value: f64| -> (usize, f64) {
        let last = axis.len() - 1;
        if value <= axis[0] {
            return (0, 0.0);
        }
        if value >= axis[last] {
            return (last - 1, 1.0);
        }
        for index in 0..last {
            if value >= axis[index] && value <= axis[index + 1] {
                let width = axis[index + 1] - axis[index];
                let fraction = if width == 0.0 {
                    0.0
                } else {
                    (value - axis[index]) / width
                };
                return (index, fraction);
            }
        }
        (last - 1, 1.0)
    };
    let (column, fx) = locate(xs, x);
    let (row, fy) = locate(ys, y);
    let at = |row: usize, column: usize| grid[row * nx + column];
    let low = at(row, column) * (1.0 - fx) + at(row, column + 1) * fx;
    let high = at(row + 1, column) * (1.0 - fx) + at(row + 1, column + 1) * fx;
    Ok(low * (1.0 - fy) + high * fy)
}

/// The coefficients of the least-squares polynomial of `degree`, ascending.
///
/// # Errors
///
/// [`InterpError::ShapeMismatch`] when the samples differ in length, and
/// [`InterpError::Domain`] when there are fewer points than coefficients.
pub fn polyfit(
    xs: &[f64],
    ys: &[f64],
    degree: usize,
) -> Result<Vec<f64>, InterpError> {
    if xs.len() != ys.len() {
        return Err(InterpError::ShapeMismatch(format!(
            "polyfit needs two samples of one length, found {} and {}",
            xs.len(),
            ys.len()
        )));
    }
    let columns = degree
        .checked_add(1)
        .ok_or_else(|| InterpError::Domain("polynomial degree overflows".into()))?;
    if xs.len() < columns {
        return Err(InterpError::Domain(format!(
            "a degree-{degree} fit needs {columns} points, found {}",
            xs.len()
        )));
    }
    let size = xs
        .len()
        .checked_mul(columns)
        .ok_or_else(|| InterpError::Domain("polynomial design overflows".into()))?;
    let mut design = vec![0.0; size];
    for (row, x) in xs.iter().enumerate() {
        let mut power = 1.0;
        for column in 0..columns {
            design[row * columns + column] = power;
            power *= x;
        }
    }
    dense::lstsq(
        &Mat {
            rows: xs.len(),
            cols: columns,
            data: design,
        },
        ys,
    )
    .map_err(|e| {
        match e {
            | DenseError::Singular => InterpError::Singular,
            | _ => InterpError::Domain(format!("polynomial fit failed: {e:?}")),
        }
    })
}

/// Evaluates a polynomial given its coefficients in ascending order.
#[must_use]
pub fn polyval(
    coefficients: &[f64],
    x: f64,
) -> f64 {
    // Horner, from the highest power down.
    coefficients
        .iter()
        .rev()
        .fold(0.0, |accumulator, coefficient| {
            accumulator * x + coefficient
        })
}

/// The monomial coefficients of the polynomial interpolating `(xs, ys)` at
/// `degree + 1` Chebyshev nodes of the sample's range.
///
/// # Errors
///
/// [`InterpError::ShapeMismatch`] when the samples differ in length, and
/// [`InterpError::Empty`] when there are none.
pub fn chebyshev_fit(
    xs: &[f64],
    ys: &[f64],
    degree: usize,
) -> Result<Vec<f64>, InterpError> {
    if xs.len() != ys.len() {
        return Err(InterpError::ShapeMismatch(format!(
            "chebyshev_fit needs two samples of one length, found {} and {}",
            xs.len(),
            ys.len()
        )));
    }
    if xs.is_empty() {
        return Err(InterpError::Empty);
    }
    let (low, high) = xs
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(low, high), x| {
            (low.min(*x), high.max(*x))
        });
    let count = degree
        .checked_add(1)
        .ok_or_else(|| InterpError::Domain("polynomial degree overflows".into()))?;
    let nodes: Vec<f64> = (0..count)
        .map(|index| {
            let angle = std::f64::consts::PI * (index as f64 + 0.5) / count as f64;
            let mapped = angle.cos();
            if high > low {
                f64::midpoint(low, high) + 0.5 * (high - low) * mapped
            } else {
                low
            }
        })
        .collect();
    // Interpolate at the nodes in the Chebyshev basis, then change basis.
    let values: Vec<f64> = nodes
        .iter()
        .map(|node| interp1(xs, ys, *node))
        .collect::<Result<_, _>>()?;
    let scaled: Vec<f64> = nodes
        .iter()
        .map(|node| {
            if high > low {
                2.0 * (node - low) / (high - low) - 1.0
            } else {
                0.0
            }
        })
        .collect();
    // Solve the Vandermonde system in the Chebyshev basis by Gaussian
    // elimination; it is small (degree + 1) and well conditioned at these
    // nodes.
    let size = count
        .checked_mul(count)
        .ok_or_else(|| InterpError::Domain("Chebyshev design overflows".into()))?;
    let mut matrix = vec![0.0; size];
    for row in 0..count {
        for column in 0..count {
            matrix[row * count + column] =
                crate::kernels::special::chebyshev_t(column as u32, scaled[row]);
        }
    }
    let coefficients = solve_flat(count, &matrix, &values)?;
    // Convert sum c_k T_k(t) with t = (2x - low - high)/(high - low) into
    // monomial coefficients in x.
    let mut total = vec![0.0; count];
    let mut previous: Vec<f64> = vec![1.0]; // T_0
    let mut current: Vec<f64> = vec![0.0, 1.0]; // T_1 in t
    for degree_index in 0..count {
        let polynomial = if degree_index == 0 {
            previous.clone()
        } else if degree_index == 1 {
            current.clone()
        } else {
            // T_k = 2 t T_{k-1} - T_{k-2}
            let mut next = vec![0.0; degree_index + 1];
            for (power, coefficient) in current.iter().enumerate() {
                next[power + 1] += 2.0 * coefficient;
            }
            for (power, coefficient) in previous.iter().enumerate() {
                next[power] -= coefficient;
            }
            previous.clone_from(&current);
            current.clone_from(&next);
            next
        };
        for (power, coefficient) in polynomial.iter().enumerate() {
            total[power] += coefficients[degree_index] * coefficient;
        }
    }
    // `total` is in the scaled variable `t`; substitute t = a x + b by Horner's
    // rule, one degree at a time.
    if high > low {
        let a = 2.0 / (high - low);
        let b = -(high + low) / (high - low);
        let mut out: Vec<f64> = vec![0.0];
        for power in (0..count).rev() {
            let mut next = vec![0.0; out.len() + 1];
            for (index, coefficient) in out.iter().enumerate() {
                next[index] += coefficient * b;
                next[index + 1] += coefficient * a;
            }
            next[0] += total[power];
            out = next;
        }
        out.truncate(count);
        out.resize(count, 0.0);
        Ok(out)
    } else {
        let mut out = vec![0.0; count];
        out[0] = total[0];
        Ok(out)
    }
}
