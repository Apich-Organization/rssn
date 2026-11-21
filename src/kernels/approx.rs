//! # Interpolation and approximation
//!
//! * cubic splines (natural, clamped, not-a-knot) with derivatives and
//!   integrals;
//! * monotone shape-preserving PCHIP (Fritsch-Carlson / Fritsch-Butland);
//! * B-splines: Cox-de Boor basis, evaluation by de Boor's algorithm,
//!   interpolation with clamped knots;
//! * Chebyshev approximation: coefficients, Clenshaw evaluation,
//!   differentiation, integration, adaptive "chebfun-style" construction
//!   and roots (colleague matrix);
//! * AAA rational approximation (real data, barycentric form) and Pade
//!   approximants from Taylor coefficients.

use std::f64::consts::PI;

use num_integer::gcd;

use crate::kernels::dense::Mat;
use crate::kernels::dense::{
    self,
};
use crate::kernels::fft::Complex64;
use crate::kernels::fft::Norm;
use crate::kernels::fft::transform_in_place;
use crate::kernels::interp_ext::shape_preserving::Segment;
use crate::kernels::interp_ext::shape_preserving::pchip_polynomials;
use crate::kernels::interpolate::scaled_grid_interval_distances;
use crate::kernels::jacobian::BandLu;
use crate::kernels::par::map_indexed;
use crate::kernels::par::try_map_indexed;
use crate::kernels::reuse::f64_scratch;
use crate::kernels::reuse::fft_scratch_for;

/// Errors of the approximation routines.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ApproxError {
    /// Inputs have inconsistent lengths, are too short or not increasing.
    Invalid,
    /// An underlying linear solve failed.
    Singular,
}

impl std::fmt::Display for ApproxError {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        match self {
            | Self::Invalid => f.write_str("invalid approximation grid or parameters"),
            | Self::Singular => f.write_str("singular approximation system"),
        }
    }
}
impl std::error::Error for ApproxError {}

fn check_grid(
    x: &[f64],
    y: &[f64],
    min: usize,
) -> Result<(), ApproxError> {
    if x.len() != y.len() || x.len() < min || x.windows(2).any(|w| w[1] <= w[0]) {
        return Err(ApproxError::Invalid);
    }
    Ok(())
}

// LDLᵀ factorization with linear time and storage. Strictly increasing
// knots make the natural and clamped systems symmetric positive definite.
fn solve_symmetric_tridiagonal(
    mut diagonal: Vec<f64>,
    mut off: Vec<f64>,
    mut rhs: Vec<f64>,
) -> Result<Vec<f64>, ApproxError> {
    let n = diagonal.len();
    if diagonal
        .iter()
        .chain(&off)
        .chain(&rhs)
        .any(|v| !v.is_finite())
    {
        return Err(ApproxError::Singular);
    }
    for i in 1..n {
        if diagonal[i - 1] <= 0.0 {
            return Err(ApproxError::Singular);
        }
        let factor = off[i - 1] / diagonal[i - 1];
        diagonal[i] -= factor * off[i - 1];
        rhs[i] -= factor * rhs[i - 1];
        off[i - 1] = factor;
        if !diagonal[i].is_finite() || !rhs[i].is_finite() || !factor.is_finite() {
            return Err(ApproxError::Singular);
        }
    }
    if diagonal[n - 1] <= 0.0 {
        return Err(ApproxError::Singular);
    }
    for i in 0..n {
        rhs[i] /= diagonal[i];
    }
    for i in (0..n - 1).rev() {
        rhs[i] -= off[i] * rhs[i + 1];
    }
    if rhs.iter().any(|v| !v.is_finite()) {
        return Err(ApproxError::Singular);
    }
    Ok(rhs)
}

/// Boundary condition of a cubic spline.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SplineBc {
    /// Zero second derivative at both ends.
    Natural,
    /// Prescribed first derivatives `(left, right)`.
    Clamped(f64, f64),
    /// Third derivative continuous across the second and last-but-one knots.
    NotAKnot,
}

/// Piecewise cubic spline `y = a + b t + c t^2 + d t^3` with `t = x - x_i`.
/// Natural and clamped construction takes O(n) time and space; not-a-knot
/// currently uses a dense O(n³)-time, O(n²)-space solve.
#[derive(Debug, Clone)]
pub struct CubicSpline {
    x: Vec<f64>,
    a: Vec<f64>,
    b: Vec<f64>,
    c: Vec<f64>,
    d: Vec<f64>,
}

impl CubicSpline {
    /// Builds an interpolating spline through `(x, y)`.
    ///
    /// # Errors
    /// [`ApproxError::Invalid`] for fewer than 3 points (4 for not-a-knot)
    /// non-increasing abscissae, or non-finite data or clamped derivatives;
    /// [`ApproxError::Singular`] if the system overflows or cannot be solved.
    pub fn new(
        x: &[f64],
        y: &[f64],
        bc: SplineBc,
    ) -> Result<Self, ApproxError> {
        let n = x.len();
        check_grid(
            x,
            y,
            if bc == SplineBc::NotAKnot {
                4
            } else {
                3
            },
        )?;
        let h: Vec<f64> = x.windows(2).map(|w| w[1] - w[0]).collect();
        // Solve for second-derivative-like coefficients c_i (c = M/2).
        if x.iter().chain(y).any(|v| !v.is_finite()) || h.iter().any(|v| !v.is_finite()) {
            return Err(ApproxError::Invalid);
        }
        let mut r = vec![0.0; n];
        for i in 1..n - 1 {
            r[i] = 3.0 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]);
        }
        if r.iter().any(|v| !v.is_finite()) {
            return Err(ApproxError::Singular);
        }
        let c = if bc == SplineBc::NotAKnot {
            // This end condition has an extra off-diagonal in each end row.
            let mut m = Mat::zeros(n, n);
            for i in 1..n - 1 {
                m.set(i, i - 1, h[i - 1]);
                m.set(i, i, 2.0 * (h[i - 1] + h[i]));
                m.set(i, i + 1, h[i]);
            }
            m.set(0, 0, h[1]);
            m.set(0, 1, -(h[0] + h[1]));
            m.set(0, 2, h[0]);
            m.set(n - 1, n - 3, h[n - 2]);
            m.set(n - 1, n - 2, -(h[n - 3] + h[n - 2]));
            m.set(n - 1, n - 1, h[n - 3]);
            dense::solve(&m, &r).map_err(|_| ApproxError::Singular)?
        } else {
            // Natural end coefficients are zero, so their couplings may be
            // eliminated. Both resulting systems are symmetric positive definite.
            let mut diagonal = vec![1.0; n];
            let mut off = vec![0.0; n - 1];
            for i in 1..n - 1 {
                diagonal[i] = 2.0 * (h[i - 1] + h[i]);
            }
            off[1..n - 2].copy_from_slice(&h[1..n - 2]);
            if let SplineBc::Clamped(left, right) = bc {
                if !left.is_finite() || !right.is_finite() {
                    return Err(ApproxError::Invalid);
                }
                diagonal[0] = 2.0 * h[0];
                diagonal[n - 1] = 2.0 * h[n - 2];
                off[0] = h[0];
                off[n - 2] = h[n - 2];
                r[0] = 3.0 * ((y[1] - y[0]) / h[0] - left);
                r[n - 1] = 3.0 * (right - (y[n - 1] - y[n - 2]) / h[n - 2]);
            }
            solve_symmetric_tridiagonal(diagonal, off, r)?
        };
        let mut b = vec![0.0; n];
        let mut d = vec![0.0; n];
        for i in 0..n - 1 {
            b[i] = (y[i + 1] - y[i]) / h[i] - h[i] * (2.0 * c[i] + c[i + 1]) / 3.0;
            d[i] = (c[i + 1] - c[i]) / (3.0 * h[i]);
        }
        // last knot: derivative / third-derivative continuation
        let hl = h[n - 2];
        b[n - 1] = b[n - 2] + 2.0 * c[n - 2] * hl + 3.0 * d[n - 2] * hl * hl;
        if b.iter().chain(&c).chain(&d).any(|v| !v.is_finite()) {
            return Err(ApproxError::Singular);
        }
        Ok(Self {
            x: x.to_vec(),
            a: y.to_vec(),
            b,
            c,
            d,
        })
    }

    fn seg(
        &self,
        x: f64,
    ) -> usize {
        self.x
            .partition_point(|&v| v <= x)
            .saturating_sub(1)
            .min(self.x.len() - 2)
    }

    /// Spline value at `x` (extrapolates with the end cubics).
    #[must_use]
    pub fn eval(
        &self,
        x: f64,
    ) -> f64 {
        let i = self.seg(x);
        let t = x - self.x[i];
        self.a[i] + t * (self.b[i] + t * (self.c[i] + t * self.d[i]))
    }

    /// Derivative of order 1, 2 or 3 at `x` (order 0 is the value).
    #[must_use]
    pub fn derivative(
        &self,
        x: f64,
        order: u8,
    ) -> f64 {
        let i = self.seg(x);
        let t = x - self.x[i];
        match order {
            | 0 => self.eval(x),
            | 1 => self.b[i] + t * (2.0 * self.c[i] + 3.0 * self.d[i] * t),
            | 2 => 2.0 * self.c[i] + 6.0 * self.d[i] * t,
            | 3 => 6.0 * self.d[i],
            | _ => 0.0,
        }
    }

    /// Exact integral of the spline over `[lo, hi]`.
    #[must_use]
    pub fn integrate(
        &self,
        lo: f64,
        hi: f64,
    ) -> f64 {
        if lo > hi {
            return -self.integrate(hi, lo);
        }
        let anti = |i: usize, t: f64| {
            t * (self.a[i] + t * (self.b[i] / 2.0 + t * (self.c[i] / 3.0 + t * self.d[i] / 4.0)))
        };
        let (il, ih) = (self.seg(lo), self.seg(hi));
        if il == ih {
            return anti(il, hi - self.x[il]) - anti(il, lo - self.x[il]);
        }
        let mut s = anti(il, self.x[il + 1] - self.x[il]) - anti(il, lo - self.x[il]);
        for i in il + 1..ih {
            s += anti(i, self.x[i + 1] - self.x[i]);
        }
        s + anti(ih, hi - self.x[ih])
    }
}

/// Monotone piecewise cubic Hermite interpolant (PCHIP).
///
/// Preserves the monotonicity of the data: no overshoot between knots.
#[derive(Debug, Clone)]
pub struct Pchip {
    data: PchipData,
}

#[derive(Debug, Clone)]
enum PchipData {
    Ordinary {
        x: Vec<f64>,
        y: Vec<f64>,
        polynomial: Vec<[f64; 3]>,
        last_tangent: f64,
    },
    Scaled(Vec<Segment>),
}

impl Pchip {
    /// Builds the interpolant.
    ///
    /// # Errors
    /// [`ApproxError::Invalid`] for non-finite data, mismatched lengths,
    /// fewer than two samples, or non-increasing coordinates.
    pub fn new(
        x: &[f64],
        y: &[f64],
    ) -> Result<Self, ApproxError> {
        let data = if let Some((polynomial, last_tangent)) = pchip_polynomials(x, y) {
            PchipData::Ordinary {
                x: x.to_vec(),
                y: y.to_vec(),
                polynomial,
                last_tangent,
            }
        } else {
            Self::validate_samples(x, y)?;
            PchipData::Scaled(map_indexed(x.len() - 1, 4096, |i| Segment::pchip(x, y, i)))
        };
        Ok(Self { data })
    }

    pub(crate) fn validate_samples(
        x: &[f64],
        y: &[f64],
    ) -> Result<(), ApproxError> {
        check_grid(x, y, 2)?;
        if x.iter().chain(y).any(|v| !v.is_finite()) {
            return Err(ApproxError::Invalid);
        }
        Ok(())
    }

    /// A one-off query over checked samples, as in a portable model payload.
    #[allow(
        clippy::single_call_fn,
        reason = "lets the portable value adapter reuse the same segment arithmetic"
    )]
    pub(crate) fn query_checked_samples(
        x: &[f64],
        y: &[f64],
        query: f64,
        derivative: bool,
    ) -> f64 {
        if !query.is_finite() {
            return f64::NAN;
        }
        let i = pchip_segment_index(x, query);
        let segment = Segment::pchip(x, y, i);
        if derivative {
            segment.derivative(query)
        } else {
            segment.value(query)
        }
    }

    #[inline]
    fn evaluate(
        &self,
        query: f64,
        derivative: bool,
    ) -> f64 {
        if !query.is_finite() {
            return f64::NAN;
        }
        let eval = |segment: &Segment| {
            if derivative {
                segment.derivative(query)
            } else {
                segment.value(query)
            }
        };
        match &self.data {
            | PchipData::Ordinary {
                x,
                y,
                polynomial,
                last_tangent,
            } => {
                let i = pchip_segment_index(x, query);
                let width = x[i + 1] - x[i];
                let t = (query - x[i]) / width;
                let [left, quadratic, cubic] = polynomial[i];
                let right = crate::kernels::exact_eq(query, x[i + 1]);
                let direct = if derivative {
                    if right {
                        last_tangent / width
                    } else {
                        ((3.0 * cubic * t + 2.0 * quadratic) * t + left) / width
                    }
                } else if right {
                    y[i + 1]
                } else if t == 0.0 && crate::kernels::exact_eq(query, x[i]) {
                    y[i]
                } else {
                    ((cubic * t + quadratic) * t + left) * t + y[i]
                };
                if direct.is_finite() && (t != 0.0 || crate::kernels::exact_eq(query, x[i])) {
                    direct
                } else {
                    eval(&Segment::pchip(x, y, i))
                }
            },
            | PchipData::Scaled(segments) => {
                let i = segments
                    .partition_point(|segment| segment.left_coordinate() <= query)
                    .saturating_sub(1);
                eval(&segments[i])
            },
        }
    }

    /// Value at `x`.
    #[must_use]
    pub fn eval(
        &self,
        x: f64,
    ) -> f64 {
        self.evaluate(x, false)
    }

    /// Derivative at `x`.
    #[must_use]
    pub fn derivative(
        &self,
        x: f64,
    ) -> f64 {
        self.evaluate(x, true)
    }
}

fn pchip_segment_index(
    samples: &[f64],
    query: f64,
) -> usize {
    samples
        .partition_point(|&point| point <= query)
        .saturating_sub(1)
        .min(samples.len() - 2)
}

/// Value of the `i`-th B-spline basis function at `x`.
///
/// Uses the Cox–de Boor recurrence on non-decreasing knots, with
/// right-continuous spans and a closed last span. An invalid index or a
/// non-finite local knot gives NaN.
#[must_use]
pub fn bspline_basis(
    knots: &[f64],
    degree: usize,
    i: usize,
    x: f64,
) -> f64 {
    let Some(local) = degree
        .checked_add(2)
        .and_then(|len| i.checked_add(len))
        .and_then(|end| knots.get(i..end))
    else {
        return f64::NAN;
    };
    if x.is_nan() || local.iter().any(|knot| !knot.is_finite()) {
        return f64::NAN;
    }
    let closed_end = crate::kernels::exact_eq(x, knots[knots.len() - 1]);
    if x < local[0] || x > local[degree + 1] || (!closed_end && x >= local[degree + 1]) {
        return 0.0;
    }
    with_spline_scratch(degree + 1, |values| {
        for (value, endpoints) in values.iter_mut().zip(local.windows(2)) {
            let (lo, hi) = (endpoints[0], endpoints[1]);
            let inside = lo <= x && x < hi;
            let at_end = closed_end && lo < hi && crate::kernels::exact_eq(x, hi);
            *value = f64::from(u8::from(inside || at_end));
        }
        for order in 1..=degree {
            for j in 0..=degree - order {
                let mut weight = 0.0;
                if values[j] != 0.0 && local[j] < local[j + order] {
                    weight += spline_fraction(local[j], local[j + order], x, false) * values[j];
                }
                if values[j + 1] != 0.0 && local[j + 1] < local[j + order + 1] {
                    weight += spline_fraction(local[j + 1], local[j + order + 1], x, true)
                        * values[j + 1];
                }
                values[j] = weight;
            }
        }
        values[0]
    })
}

fn with_spline_scratch<T>(
    size: usize,
    compute: impl FnOnce(&mut [f64]) -> T,
) -> T {
    let mut small = [0.0; 16];
    if size <= small.len() {
        return compute(&mut small[..size]);
    }
    let mut scratch = f64_scratch(size);
    scratch.resize(size, 0.0);
    compute(&mut scratch)
}

fn spline_fraction(
    lo: f64,
    hi: f64,
    x: f64,
    from_right: bool,
) -> f64 {
    let width = hi - lo;
    let (left, right, width) = if width.is_finite() {
        (x - lo, hi - x, width)
    } else {
        scaled_grid_interval_distances(lo, hi, x)
    };
    (if from_right { right } else { left }) / width
}

fn spline_span(
    knots: &[f64],
    degree: usize,
    count: usize,
    x: f64,
) -> usize {
    let support = &knots[degree..=count];
    let offset = if x >= knots[count] {
        support.partition_point(|&knot| knot < x)
    } else {
        support.partition_point(|&knot| knot <= x)
    };
    degree + offset.saturating_sub(1).min(count - degree - 1)
}

/// Evaluates a B-spline curve `sum c_i B_{i,p}(x)` by de Boor's algorithm.
///
/// Requires finite, non-decreasing `knots`, `degree < coeffs.len()` and
/// `knots.len() == coeffs.len() + degree + 1`. The coordinate is clamped to
/// the valid range. Returns NaN for an invalid shape or endpoint.
#[must_use]
pub fn bspline_eval(
    knots: &[f64],
    coeffs: &[f64],
    degree: usize,
    x: f64,
) -> f64 {
    let count = coeffs.len();
    let expected = count.checked_add(degree).and_then(|n| n.checked_add(1));
    if degree >= count || expected != Some(knots.len()) || x.is_nan() {
        return f64::NAN;
    }
    let (lo, hi) = (knots[degree], knots[count]);
    if !lo.is_finite() || !hi.is_finite() || lo > hi {
        return f64::NAN;
    }
    bspline_eval_validated(knots, coeffs, degree, x)
}

/// Evaluates a shape-checked spline with a non-decreasing knot vector.
pub(super) fn bspline_eval_validated(
    knots: &[f64],
    coeffs: &[f64],
    degree: usize,
    x: f64,
) -> f64 {
    let count = coeffs.len();
    let x = x.clamp(knots[degree], knots[count]);
    let k = spline_span(knots, degree, count, x);
    if degree == 0 {
        return coeffs[k];
    }
    with_spline_scratch(degree + 1, |d| {
        d.copy_from_slice(&coeffs[k - degree..=k]);
        for r in 1..=degree {
            for j in (r..=degree).rev() {
                let i = j + k - degree;
                let (lo, hi) = (knots[i], knots[i + degree + 1 - r]);
                let alpha = if crate::kernels::exact_eq(lo, hi) {
                    0.0
                } else {
                    spline_fraction(lo, hi, x, false)
                };
                d[j] = (1.0 - alpha) * d[j - 1] + alpha * d[j];
            }
        }
        d[degree]
    })
}

#[allow(
    clippy::single_call_fn,
    reason = "keeps the basis recurrence separate from collocation row scheduling"
)]
fn spline_basis_on_span(
    knots: &[f64],
    span: usize,
    x: f64,
    values: &mut [f64],
) {
    values[0] = 1.0;
    for degree in 1..values.len() {
        let mut saved = 0.0;
        for (i, value) in values[..degree].iter_mut().enumerate() {
            let lo = knots[span + i + 1 - degree];
            let hi = knots[span + i + 1];
            let previous = *value;
            if previous == 0.0 || lo >= hi {
                *value = saved;
                saved = 0.0;
            } else {
                *value = saved + spline_fraction(lo, hi, x, true) * previous;
                saved = spline_fraction(lo, hi, x, false) * previous;
            }
        }
        values[degree] = saved;
    }
}

type SplineEntry = (usize, usize, f64);

#[allow(
    clippy::single_call_fn,
    reason = "builds the collocation matrix independently of the linear solver"
)]
fn spline_collocation(
    knots: &[f64],
    x: &[f64],
    degree: usize,
) -> Result<Vec<SplineEntry>, ApproxError> {
    const ROWS: usize = 64;
    let row_len = degree + 1;
    let work = ROWS.saturating_mul(row_len.saturating_mul(row_len));
    let grain = 16_384_usize.div_ceil(work);
    let count = x.len();
    let total = count.checked_mul(row_len).ok_or(ApproxError::Invalid)?;
    let blocks = try_map_indexed(count.div_ceil(ROWS), grain, |block| {
        let start = block * ROWS;
        let xs = &x[start..(start + ROWS).min(count)];
        let mut entries = Vec::new();
        entries
            .try_reserve_exact(xs.len() * row_len)
            .map_err(|_| ApproxError::Invalid)?;
        with_spline_scratch(row_len, |values| {
            for (offset, &coordinate) in xs.iter().enumerate() {
                let span = spline_span(knots, degree, count, coordinate);
                spline_basis_on_span(knots, span, coordinate, values);
                let row = start + offset;
                let first = span - degree;
                for (column, &value) in (first..).zip(values.iter()) {
                    if value != 0.0 {
                        entries.push((row, column, value));
                    }
                }
            }
        });
        Ok(entries)
    })?;
    let mut entries = Vec::new();
    entries
        .try_reserve_exact(total)
        .map_err(|_| ApproxError::Invalid)?;
    entries.extend(blocks.into_iter().flatten());
    Ok(entries)
}

#[allow(
    clippy::single_call_fn,
    reason = "keeps the overflow-safe average separate from knot assembly"
)]
fn spline_knot_average(values: &[f64]) -> f64 {
    let sum: f64 = values.iter().sum();
    if sum.is_finite() {
        return sum / values.len() as f64;
    }
    let (lo, hi) = (values[0], values[values.len() - 1]);
    let scale = lo.abs().max(hi.abs());
    let scaled: f64 = values.iter().map(|v| v / scale).sum();
    ((scaled / values.len() as f64) * scale).clamp(lo, hi)
}

/// Fits an interpolating B-spline and returns `(knots, coefficients)`.
///
/// Degree zero places boundaries at the midpoints between consecutive sites.
/// Higher degrees use a clamped knot vector with averaged interior knots.
///
/// # Errors
/// [`ApproxError::Invalid`] for inconsistent or non-finite samples, an
/// invalid degree, unordered sites, or an unrepresentable final degree-zero
/// interval. [`ApproxError::Singular`] if the collocation system is singular
/// or overflows.
pub fn bspline_interpolate(
    x: &[f64],
    y: &[f64],
    degree: usize,
) -> Result<(Vec<f64>, Vec<f64>), ApproxError> {
    let n = x.len();
    if degree >= n || x.iter().chain(y).any(|v| !v.is_finite()) {
        return Err(ApproxError::Invalid);
    }
    check_grid(x, y, degree + 1)?;
    if degree == 0 {
        let mut knots = Vec::with_capacity(n + 1);
        knots.push(x[0]);
        for pair in x.windows(2) {
            let midpoint = f64::midpoint(pair[0], pair[1]);
            knots.push(if midpoint > pair[0] {
                midpoint
            } else {
                pair[1]
            });
        }
        let upper = x[n - 1];
        let upper = if n > 1 && knots.last().is_some_and(|&knot| knot >= upper) {
            upper.next_up()
        } else {
            upper
        };
        if !upper.is_finite() {
            return Err(ApproxError::Invalid);
        }
        knots.push(upper);
        return Ok((knots, y.to_vec()));
    }
    let mut knots = vec![x[0]; degree + 1];
    for j in 1..n - degree {
        knots.push(spline_knot_average(&x[j..j + degree]));
    }
    knots.extend(std::iter::repeat_n(x[n - 1], degree + 1));
    if degree == 1 {
        return Ok((knots, y.to_vec()));
    }
    let entries = spline_collocation(&knots, x, degree)?;
    let coefficients = if degree.saturating_mul(3) < n {
        BandLu::factor(n, degree, degree, &entries)
            .ok_or(ApproxError::Singular)?
            .solve(y)
    } else {
        let mut a = Mat::zeros(n, n);
        for (i, j, value) in entries {
            a.set(i, j, value);
        }
        dense::solve(&a, y).map_err(|_| ApproxError::Singular)?
    };
    if coefficients.iter().any(|v| !v.is_finite()) {
        return Err(ApproxError::Singular);
    }
    Ok((knots, coefficients))
}

/// Chebyshev points of the second kind (extrema)
/// `cos(pi k / n)` mapped to `[a, b]`, `k = 0..=n`.
///
/// Degree zero gives the midpoint. For positive degree, the first and last
/// points are exactly `b` and `a`.
#[must_use]
pub fn chebyshev_points(
    n: usize,
    a: f64,
    b: f64,
) -> Vec<f64> {
    let interval = ChebInterval::new(a, b);
    (0..=n).map(|k| interval.node(k, n)).collect()
}

/// A shared affine map between the physical interval and `[-1, 1]`.
struct ChebInterval {
    a: f64,
    b: f64,
    midpoint: f64,
    width: f64,
    half_width: f64,
    use_endpoint_distances: bool,
}

impl ChebInterval {
    fn new(
        a: f64,
        b: f64,
    ) -> Self {
        let midpoint = f64::midpoint(a, b);
        let width = b - a;
        let half_width = if width.is_finite() {
            0.5 * width
        } else {
            0.5 * b - 0.5 * a
        };
        Self {
            a,
            b,
            midpoint,
            width,
            half_width,
            use_endpoint_distances: width.is_subnormal() || width.abs() < midpoint.abs(),
        }
    }

    fn node(
        &self,
        k: usize,
        degree: usize,
    ) -> f64 {
        if degree == 0 || degree.is_multiple_of(2) && k == degree / 2 {
            return self.midpoint;
        }
        if k == 0 {
            return self.b;
        }
        if k == degree {
            return self.a;
        }
        let t = (PI * (k as f64 / degree as f64)).cos();
        let x = if self.use_endpoint_distances {
            // Anchor at a real endpoint when the midpoint has been rounded
            // by much of a step in the physical coordinates.
            f64::midpoint(t, 1.0).mul_add(self.width, self.a)
        } else {
            self.midpoint + self.half_width * t
        };
        let lo = self.a.min(self.b);
        let hi = self.a.max(self.b);
        if x < lo {
            lo
        } else if x > hi {
            hi
        } else {
            x
        }
    }

    fn normalize(
        &self,
        x: f64,
    ) -> f64 {
        if self.width == 0.0 {
            return 0.0;
        }
        if crate::kernels::exact_eq(x, self.a) {
            return -1.0;
        }
        if crate::kernels::exact_eq(x, self.b) {
            return 1.0;
        }
        if self.use_endpoint_distances {
            ((x - self.a) - (self.b - x)) / self.width
        } else {
            (x - self.midpoint) / self.half_width
        }
    }
}

/// A function on `[a, b]` represented by Chebyshev coefficients.
#[derive(Debug, Clone, PartialEq)]
pub struct ChebSeries {
    /// Coefficients `c_0 .. c_n` of `sum c_k T_k(t)`.
    pub coeffs: Vec<f64>,
    /// Left end of the interval.
    pub a: f64,
    /// Right end of the interval.
    pub b: f64,
}

/// Clenshaw evaluation of `sum c_k T_k(t)` at `t in [-1, 1]`.
#[must_use]
pub fn clenshaw(
    c: &[f64],
    t: f64,
) -> f64 {
    if c.len() <= 1 {
        return c.first().copied().unwrap_or(0.0);
    }
    let value = clenshaw_scaled(c, t, 1.0);
    if !value.is_normal() && t.abs() <= 1.0 {
        return clenshaw_rescaled(c, t, value);
    }
    value
}

#[cold]
#[allow(clippy::single_call_fn)]
fn clenshaw_rescaled(
    c: &[f64],
    t: f64,
    direct_value: f64,
) -> f64 {
    let scale = chebyshev_scale(c.iter().fold(0.0_f64, |max, v| max.max(v.abs())));
    if crate::kernels::exact_eq(scale, 1.0) {
        direct_value
    } else {
        clenshaw_scaled(c, t, scale) / scale
    }
}

#[inline]
fn clenshaw_scaled(
    c: &[f64],
    t: f64,
    scale: f64,
) -> f64 {
    let (mut b1, mut b2) = (0.0, 0.0);
    for &ck in c.iter().skip(1).rev() {
        let b0 = ck * scale + 2.0 * t * b1 - b2;
        b2 = b1;
        b1 = b0;
    }
    c.first().copied().unwrap_or(0.0) * scale + t * b1 - b2
}

/// Keep linear transforms of very large or small finite values away from
/// range limits. Powers of two leave normal significands unchanged.
fn chebyshev_scale(maximum: f64) -> f64 {
    const LARGE: f64 = f64::from_bits(1535_u64 << 52); // 2^512
    const SMALL: f64 = f64::from_bits(511_u64 << 52); // 2^-512
    if maximum.is_finite() && maximum > LARGE {
        SMALL
    } else if maximum > 0.0 && maximum < SMALL {
        LARGE
    } else {
        1.0
    }
}

#[allow(clippy::single_call_fn)]
fn cheb_coefficients(values: &[f64]) -> Vec<f64> {
    let degree = values.len() - 1;
    let first = values[0];
    if first.is_finite() && values.iter().all(|&v| crate::kernels::exact_eq(v, first)) {
        let mut coefficients = vec![0.0; values.len()];
        coefficients[0] = first;
        return coefficients;
    }

    let scale = chebyshev_scale(values.iter().fold(0.0_f64, |max, v| max.max(v.abs())));
    let n = degree as f64;
    let normalize = |j: usize, dct: f64| {
        let endpoint_weight = if j == 0 || j == degree {
            0.5
        } else {
            1.0
        };
        (dct * endpoint_weight / n) / scale
    };
    if degree <= 8 {
        // Tiny transforms avoid the FFT plan lookup. A small table covers all
        // cosine factors once their index is reduced by periodicity.
        let period = 2 * degree;
        let mut cosine = [0.0; 9];
        cosine[0] = 1.0;
        cosine[degree] = -1.0;
        for (k, value) in cosine.iter_mut().enumerate().take(degree).skip(1) {
            if k != degree / 2 || !degree.is_multiple_of(2) {
                *value = (PI * (k as f64 / n)).cos();
            }
        }
        let mut weighted = [0.0; 9];
        for (k, (&value, weight)) in values.iter().zip(&mut weighted).enumerate() {
            *weight = if k == 0 || k == degree {
                value * scale
            } else {
                2.0 * (value * scale)
            };
        }
        return (0..=degree)
            .map(|j| {
                let (mut sum, mut phase) = (0.0, 0);
                for &value in &weighted[..=degree] {
                    sum += value * cosine[phase.min(period - phase)];
                    phase += j;
                    if phase >= period {
                        phase -= period;
                    }
                }
                normalize(j, sum)
            })
            .collect();
    }

    // DCT-I is the real part of the DFT of the even extension. The RSSN
    // FFT implementation shares the plan and borrows its own scratch space.
    let mut reflected = fft_scratch_for(2 * degree);
    reflected.extend(
        values
            .iter()
            .chain(values[1..degree].iter().rev())
            .map(|&value| Complex64::new(value * scale, 0.0)),
    );
    transform_in_place(&mut reflected, false, Norm::Backward);
    reflected
        .iter()
        .take(values.len())
        .enumerate()
        .map(|(j, value)| normalize(j, value.re))
        .collect()
}

impl ChebSeries {
    /// Interpolates `f` at the Chebyshev extrema in [`chebyshev_points`] order.
    ///
    /// Degree zero and equal endpoints use one function evaluation.
    pub fn from_fn<F: Fn(f64) -> f64>(
        f: F,
        a: f64,
        b: f64,
        n: usize,
    ) -> Self {
        let interval = ChebInterval::new(a, b);
        if interval.width == 0.0 || n == 0 {
            return Self::from_values(&[f(interval.midpoint)], a, b);
        }
        let values: Vec<_> = (0..=n).map(|k| f(interval.node(k, n))).collect();
        Self::from_values(&values, a, b)
    }

    /// Interpolates already sampled values in [`chebyshev_points`] order.
    ///
    /// `values[0]` is the sample at `b` for two or more values. One value
    /// represents a constant on the whole interval. An empty slice
    /// represents zero. With equal endpoints the first value is used.
    #[must_use]
    pub fn from_values(
        values: &[f64],
        a: f64,
        b: f64,
    ) -> Self {
        let coeffs = if values.len() <= 1 || crate::kernels::exact_eq(a, b) {
            vec![values.first().copied().unwrap_or(0.0)]
        } else {
            cheb_coefficients(values)
        };
        Self { coeffs, a, b }
    }

    /// Adaptive construction: doubles the degree (up to `max_degree`)
    /// until the trailing coefficients fall below `tol * max|c|`, then
    /// chops the negligible tail. Values at shared nodes are reused when
    /// refining the grid. The same coordinate should always give the same
    /// value when supplied to `f`.
    pub fn adaptive<F: Fn(f64) -> f64>(
        f: F,
        a: f64,
        b: f64,
        tol: f64,
        max_degree: usize,
    ) -> Self {
        let interval = ChebInterval::new(a, b);
        if interval.width == 0.0 || max_degree == 0 {
            return Self::from_values(&[f(interval.midpoint)], a, b);
        }
        let mut n = 8.min(max_degree);
        let mut values: Vec<_> = (0..=n).map(|k| f(interval.node(k, n))).collect();
        loop {
            let mut s = Self::from_values(&values, a, b);
            let big = s.coeffs.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
            let cut = tol * big.max(f64::MIN_POSITIVE);
            let tail_ok = s
                .coeffs
                .iter()
                .rev()
                .take(3.min(s.coeffs.len()))
                .all(|v| v.abs() <= cut);
            if tail_ok || n >= max_degree {
                while s.coeffs.len() > 1 && s.coeffs.last().is_some_and(|v| v.abs() <= cut) {
                    s.coeffs.pop();
                }
                return s;
            }
            let next = n.saturating_mul(2).min(max_degree);
            let shared = gcd(n, next);
            let (old_step, new_step) = (n / shared, next / shared);
            values.resize(next + 1, 0.0);
            // Move shared samples from the back so no source is overwritten.
            for k in (1..=shared).rev() {
                values[k * new_step] = values[k * old_step];
            }
            for (k, value) in values.iter_mut().enumerate().skip(1) {
                if !k.is_multiple_of(new_step) {
                    *value = f(interval.node(k, next));
                }
            }
            n = next;
        }
    }

    /// Degree of the series.
    #[must_use]
    pub const fn degree(&self) -> usize {
        self.coeffs.len().saturating_sub(1)
    }

    /// Value at `x in [a, b]`.
    #[must_use]
    pub fn eval(
        &self,
        x: f64,
    ) -> f64 {
        let t = if crate::kernels::exact_eq(self.a, -self.b) && self.b != 0.0 && self.b.is_finite()
        {
            x / self.b
        } else {
            ChebInterval::new(self.a, self.b).normalize(x)
        };
        clenshaw(&self.coeffs, t)
    }

    /// Derivative series.
    #[must_use]
    pub fn derivative(&self) -> Self {
        let n = self.coeffs.len();
        if n <= 1 {
            return Self::from_values(&[], self.a, self.b);
        }
        let mut d = vec![0.0; n.max(2) + 1];
        for k in (1..n).rev() {
            d[k - 1] = d[k + 1] + 2.0 * k as f64 * self.coeffs[k];
        }
        d.truncate(n.saturating_sub(1).max(1));
        d[0] *= 0.5;
        let s = 2.0 / (self.b - self.a);
        Self {
            coeffs: d.iter().map(|v| v * s).collect(),
            a: self.a,
            b: self.b,
        }
    }

    /// Antiderivative series vanishing at `a`.
    #[must_use]
    pub fn integral(&self) -> Self {
        if crate::kernels::exact_eq(self.a, self.b) {
            return Self::from_values(&[], self.a, self.b);
        }
        let n = self.coeffs.len();
        let scale = 0.5 * (self.b - self.a);
        let mut c = vec![0.0; n + 1];
        let get = |k: usize| self.coeffs.get(k).copied().unwrap_or(0.0);
        for k in 1..=n {
            let hi = get(k + 1);
            let lo = if k == 1 {
                2.0 * get(0)
            } else {
                get(k - 1)
            };
            c[k] = (lo - hi) / (2.0 * k as f64) * scale;
        }
        let at_minus_one: f64 = c
            .iter()
            .enumerate()
            .skip(1)
            .map(|(k, v)| if k % 2 == 0 { *v } else { -*v })
            .sum();
        c[0] = -at_minus_one;
        Self {
            coeffs: c,
            a: self.a,
            b: self.b,
        }
    }

    /// Definite integral over `[a, b]`.
    #[must_use]
    pub fn definite_integral(&self) -> f64 {
        if crate::kernels::exact_eq(self.a, self.b) {
            return 0.0;
        }
        let s = 0.5 * (self.b - self.a);
        let mut sum = 0.0;
        for (k, &c) in self.coeffs.iter().enumerate().step_by(2) {
            sum += c * 2.0 / (1.0 - (k * k) as f64);
        }
        s * sum
    }

    /// Real roots in `[a, b]` from the eigenvalues of the colleague matrix.
    ///
    /// # Errors
    /// [`ApproxError::Singular`] if the eigenvalue iteration fails.
    pub fn roots(&self) -> Result<Vec<f64>, ApproxError> {
        let mut c = self.coeffs.clone();
        while c.len() > 1 && c.last().is_some_and(|v| *v == 0.0) {
            c.pop();
        }
        let n = c.len() - 1;
        if n == 0 {
            return Ok(Vec::new());
        }
        let mut m = Mat::zeros(n, n);
        for i in 0..n - 1 {
            m.set(i, i + 1, 0.5);
            m.set(i + 1, i, 0.5);
        }
        if n > 1 {
            m.set(0, 1, 1.0);
        }
        if n == 1 {
            m.set(0, 0, -c[0] / c[1]);
        } else {
            for j in 0..n {
                let v = m.at(n - 1, j) - 0.5 * c[j] / c[n];
                m.set(n - 1, j, v);
            }
        }
        let ev = dense::eigenvalues(&m).map_err(|_| ApproxError::Singular)?;
        let mut out: Vec<f64> = ev
            .iter()
            .filter(|&&(re, im)| im.abs() < 1e-8 && re.abs() <= 1.0 + 1e-10)
            .map(|&(re, _)| f64::midpoint(self.a, self.b) + 0.5 * (self.b - self.a) * re)
            .collect();
        out.sort_by(f64::total_cmp);
        Ok(out)
    }
}

/// Barycentric rational approximant produced by [`aaa`].
#[derive(Debug, Clone)]
pub struct Aaa {
    /// Support points.
    pub z: Vec<f64>,
    /// Function values at the support points.
    pub f: Vec<f64>,
    /// Barycentric weights.
    pub w: Vec<f64>,
    /// Final maximal sample error.
    pub error: f64,
}

impl Aaa {
    /// Evaluates the rational function at `x`.
    #[must_use]
    pub fn eval(
        &self,
        x: f64,
    ) -> f64 {
        let (mut num, mut den) = (0.0, 0.0);
        for k in 0..self.z.len() {
            let d = x - self.z[k];
            if d == 0.0 {
                return self.f[k];
            }
            let t = self.w[k] / d;
            num += t * self.f[k];
            den += t;
        }
        num / den
    }

    /// Degree (`numerator = denominator = m - 1`).
    #[must_use]
    pub const fn degree(&self) -> usize {
        self.z.len().saturating_sub(1)
    }
}

/// AAA algorithm (Nakatsukasa, Sete, Trefethen) on real samples
/// `(z_i, f_i)`: greedy support-point selection and a least-squares
/// weight vector from the SVD of the Loewner matrix.
///
/// # Errors
/// [`ApproxError::Invalid`] for fewer than 3 samples.
pub fn aaa(
    z: &[f64],
    f: &[f64],
    tol: f64,
    max_terms: usize,
) -> Result<Aaa, ApproxError> {
    let m = z.len();
    if m < 3 || f.len() != m {
        return Err(ApproxError::Invalid);
    }
    let fmax = f.iter().fold(0.0_f64, |a, v| a.max(v.abs()));
    let mut rest: Vec<usize> = (0..m).collect();
    let mean = f.iter().sum::<f64>() / m as f64;
    let mut r = vec![mean; m];
    let (mut zj, mut fj): (Vec<f64>, Vec<f64>) = (Vec::new(), Vec::new());
    let mut w = Vec::new();
    let mut err = f64::INFINITY;
    for _ in 0..max_terms.min(m - 1) {
        let (pos, &j) = rest
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| (f[**a] - r[**a]).abs().total_cmp(&(f[**b] - r[**b]).abs()))
            .unwrap_or((0, &0));
        zj.push(z[j]);
        fj.push(f[j]);
        rest.remove(pos);
        let k = zj.len();
        if rest.len() < k {
            break;
        }
        let mut a = Mat::zeros(rest.len(), k);
        for (row, &i) in rest.iter().enumerate() {
            for c in 0..k {
                a.set(row, c, (f[i] - fj[c]) / (z[i] - zj[c]));
            }
        }
        let Ok(s) = dense::svd(&a) else { break };
        let last = k - 1;
        w = (0..k).map(|c| s.v.at(c, last)).collect();
        for i in 0..m {
            let (mut num, mut den) = (0.0, 0.0);
            let mut hit = None;
            for c in 0..k {
                let d = z[i] - zj[c];
                if d == 0.0 {
                    hit = Some(c);
                    break;
                }
                num += w[c] * fj[c] / d;
                den += w[c] / d;
            }
            r[i] = hit.map_or(num / den, |c| fj[c]);
        }
        err = rest
            .iter()
            .map(|&i| (f[i] - r[i]).abs())
            .fold(0.0, f64::max);
        if err <= tol * fmax {
            break;
        }
    }
    Ok(Aaa {
        z: zj,
        f: fj,
        w,
        error: err,
    })
}

/// Pade approximant `[m/n]` from Taylor coefficients `c_0, c_1, ...`
/// (needs `m + n + 1` coefficients). Returns `(numerator, denominator)`
/// with ascending coefficients and `denominator[0] = 1`.
///
/// # Errors
/// [`ApproxError::Invalid`] if too few coefficients, [`ApproxError::Singular`]
/// if the approximant does not exist.
pub fn pade(
    c: &[f64],
    m: usize,
    n: usize,
) -> Result<(Vec<f64>, Vec<f64>), ApproxError> {
    if c.len() < m + n + 1 {
        return Err(ApproxError::Invalid);
    }
    let cc = |i: isize| if i < 0 { 0.0 } else { c[i as usize] };
    let mut q = vec![1.0];
    if n > 0 {
        let mut a = Mat::zeros(n, n);
        let mut b = vec![0.0; n];
        for i in 0..n {
            for j in 0..n {
                a.set(i, j, cc((m + i) as isize - j as isize));
            }
            b[i] = -cc((m + i + 1) as isize);
        }
        q.extend(dense::solve(&a, &b).map_err(|_| ApproxError::Singular)?);
    }
    let p: Vec<f64> = (0..=m)
        .map(|k| (0..=k.min(n)).map(|j| q[j] * cc((k - j) as isize)).sum())
        .collect();
    Ok((p, q))
}

/// Evaluates a polynomial with ascending coefficients by Horner's scheme.
#[must_use]
pub fn horner(
    c: &[f64],
    x: f64,
) -> f64 {
    c.iter().rev().fold(0.0, |acc, &v| acc * x + v)
}
