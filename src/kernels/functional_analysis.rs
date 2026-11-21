//! # Numerical Functional Analysis
//!
//! This module provides numerical implementations for concepts from functional analysis.
//! It includes functions for calculating various norms (L1, L2, L-infinity) and inner products
//! for functions represented by discrete points. [`L2Grid`] reuses the
//! quadrature weights when several functions share the same coordinates.

use std::f64::consts::FRAC_1_SQRT_2;

/// Calculates the L1 norm of a function represented by discrete points.
///
/// The L1 norm (or Manhattan norm) is defined as `∫|f(x)|dx`. For discrete points,
/// it is approximated as `Σ|y_i|Δx_i`.
///
/// # Arguments
/// * `points` - A slice of `(x, y)` tuples representing the function's samples.
///
/// # Returns
/// The numerical value of the L1 norm.
#[must_use]
pub fn l1_norm(points: &[(f64, f64)]) -> f64 {
    points
        .windows(2)
        .map(|w| {
            let (x1, y1) = w[0];

            let (x2, y2) = w[1];

            f64::midpoint(y1.abs(), y2.abs()) * (x2 - x1)
        })
        .sum()
}

/// Calculates the L2 norm of a function represented by discrete points.
///
/// The L2 norm is `sqrt(∫|f(x)|²dx)`. For samples, the integral is the
/// trapezoidal rule applied to the squared ordinates. On a finite
/// nondecreasing grid, repeated coordinates give zero-width panels.
///
/// # Arguments
/// * `points` - A slice of `(x, y)` tuples representing the function's samples.
///
/// # Returns
/// The numerical value of the L2 norm.
#[must_use]
pub fn l2_norm(points: &[(f64, f64)]) -> f64 {
    l2_measure(points).value()
}

/// A trapezoidal L2 rule prepared for several functions on the same grid.
///
/// The cached roots of its quadrature weights depend only on the physical
/// coordinates. A call supplies the corresponding ordinates in grid order.
#[derive(Clone, Debug)]
pub struct L2Grid {
    root_weights: Vec<f64>,
}

impl L2Grid {
    /// Prepares the weights for finite, nondecreasing coordinates.
    ///
    /// # Errors
    /// Returns an error for a nonfinite coordinate or a decrease.
    pub fn new(coordinates: &[f64]) -> Result<Self, String> {
        if coordinates.iter().any(|x| !x.is_finite())
            || coordinates.windows(2).any(|window| window[0] > window[1])
        {
            return Err("L2 grid coordinates must be finite and nondecreasing".into());
        }
        Ok(Self {
            root_weights: trapezoid_weight_roots(coordinates.len(), |i| coordinates[i]).collect(),
        })
    }

    /// Computes the sampled L2 norm using the prepared weights.
    ///
    /// # Errors
    /// Returns an error for the wrong number of ordinates or a nonfinite one.
    pub fn norm(
        &self,
        ordinates: &[f64],
    ) -> Result<f64, String> {
        Ok(self.measure(ordinates)?.value())
    }

    /// Returns the L2 norm and normalized ordinates from one reduction.
    ///
    /// Zero measure leaves the ordinates as supplied. If the rounded norm
    /// overflows or underflows, normalization still uses its internal scale.
    ///
    /// # Errors
    /// Returns an error for the wrong number of ordinates or a nonfinite one.
    pub fn norm_and_normalize(
        &self,
        ordinates: &[f64],
    ) -> Result<(f64, Vec<f64>), String> {
        let (norm, normalized, _) = self.normalize_with_measure(ordinates)?;
        Ok((norm, normalized))
    }

    pub(crate) fn normalize_with_measure(
        &self,
        ordinates: &[f64],
    ) -> Result<(f64, Vec<f64>, bool), String> {
        let measure = self.measure(ordinates)?;
        let has_measure = measure.magnitude != 0.0;
        let normalized = if !has_measure {
            ordinates.to_vec()
        } else if measure.exponent == 0 {
            ordinates
                .iter()
                .map(|&value| value / measure.magnitude)
                .collect()
        } else {
            ordinates
                .iter()
                .map(|&value| measure.normalize_scaled(value))
                .collect()
        };
        Ok((measure.value(), normalized, has_measure))
    }

    pub(crate) const fn heap_bytes(&self) -> usize {
        self.root_weights.capacity() * size_of::<f64>()
    }

    fn measure(
        &self,
        ordinates: &[f64],
    ) -> Result<L2Measure, String> {
        if ordinates.len() != self.root_weights.len() {
            return Err("L2 ordinates must match the number of coordinates".into());
        }
        let samples = self
            .root_weights
            .iter()
            .copied()
            .zip(ordinates.iter().copied());
        let square: f64 = samples
            .clone()
            .map(|(root, value)| (root * value).powi(2))
            .sum();
        // A weighted sample that underflows before squaring cannot affect
        // a squared sum with this much headroom.
        if squared_sum_has_headroom(square, ordinates.len()) {
            Ok(L2Measure::direct(square.sqrt()))
        } else {
            checked_scaled_measure(samples)
        }
    }
}

/// The norm is `magnitude * 2^exponent`. Keeping the exponent separate
/// also supports normalization when the norm itself cannot fit in f64.
struct L2Measure {
    magnitude: f64,
    exponent: i32,
}

impl L2Measure {
    const fn direct(magnitude: f64) -> Self {
        Self {
            magnitude,
            exponent: 0,
        }
    }

    fn value(&self) -> f64 {
        scale_positive(self.magnitude, self.exponent)
    }

    fn normalize_scaled(
        &self,
        value: f64,
    ) -> f64 {
        if value == 0.0 {
            return value;
        }
        let (mantissa, exponent) = split_positive(value.abs());
        scale_positive(mantissa / self.magnitude, exponent - self.exponent).copysign(value)
    }
}

const fn squared_sum_has_headroom(
    square: f64,
    terms: usize,
) -> bool {
    // Each underflowing term can lose half a least subnormal unit. This
    // total bounds their combined loss to half an ulp of the sum.
    square >= (2.0 * f64::MIN_POSITIVE) * terms as f64 && square.is_finite()
}

#[inline]
fn window_square(
    left: f64,
    right: f64,
) -> f64 {
    right.mul_add(right, left * left) * 0.5
}

fn l2_measure(points: &[(f64, f64)]) -> L2Measure {
    if points.len() < 2 {
        return L2Measure::direct(0.0);
    }
    let mut square = 0.0;
    for &[(x0, y0), (x1, y1)] in points.array_windows::<2>() {
        let average = window_square(y0, y1);
        // Squaring may lose an ordinate before a large width restores it.
        if average < f64::MIN_POSITIVE && (y0 != 0.0 || y1 != 0.0) {
            return scaled_points_measure(points);
        }
        square += average * (x1 - x0);
    }
    if squared_sum_has_headroom(square, points.len() - 1) {
        L2Measure::direct(square.sqrt())
    } else if square == 0.0 && points.iter().all(|&(_, y)| y == 0.0) {
        L2Measure::direct(0.0)
    } else {
        scaled_points_measure(points)
    }
}

#[cold]
fn scaled_points_measure(points: &[(f64, f64)]) -> L2Measure {
    if points
        .iter()
        .any(|&(x, y)| !x.is_finite() || !y.is_finite())
        || points.windows(2).any(|window| window[0].0 > window[1].0)
    {
        // Preserve the signed trapezoidal arithmetic of the unchecked API
        // for inputs outside its finite, nondecreasing grid.
        let square: f64 = points
            .windows(2)
            .map(|window| window_square(window[0].1, window[1].1) * (window[1].0 - window[0].0))
            .sum();
        return L2Measure::direct(square.sqrt());
    }
    let roots = trapezoid_weight_roots(points.len(), |i| points[i].0);
    scaled_l2_measure(roots.zip(points.iter().map(|&(_, y)| y)))
}

/// A node's weight is half the distance between the adjacent outer
/// coordinates. At an endpoint, use the coordinate itself on that side.
fn trapezoid_weight_roots(
    len: usize,
    coordinate: impl Fn(usize) -> f64,
) -> impl Iterator<Item = f64> {
    (0..len).map(move |i| {
        let left = coordinate(i.saturating_sub(1));
        let right = coordinate((i + 1).min(len - 1));
        let width = right - left;
        if !width.is_finite() {
            (0.5 * right - 0.5 * left).sqrt()
        } else if width < 2.0 * f64::MIN_POSITIVE {
            // Halving a subnormal width first could discard it.
            width.sqrt() * FRAC_1_SQRT_2
        } else {
            (0.5 * width).sqrt()
        }
    })
}

#[cold]
#[allow(clippy::single_call_fn)]
fn checked_scaled_measure(
    samples: impl Iterator<Item = (f64, f64)> + Clone
) -> Result<L2Measure, String> {
    if samples.clone().any(|(_, value)| !value.is_finite()) {
        return Err("L2 ordinates must be finite".into());
    }
    Ok(scaled_l2_measure(samples))
}

/// Sum weighted squares with a common binary exponent. A weighted sample
/// may be outside the representable range even before it is squared.
fn scaled_l2_measure(samples: impl Iterator<Item = (f64, f64)>) -> L2Measure {
    let mut exponent = 0;
    let mut square = 0.0;
    for (root, value) in samples {
        if root == 0.0 || value == 0.0 {
            continue;
        }
        let (root_mantissa, root_exponent) = split_positive(root);
        let (value_mantissa, value_exponent) = split_positive(value.abs());
        let product = root_mantissa * value_mantissa;
        let term = product * product;
        let term_exponent = root_exponent + value_exponent;
        if square == 0.0 {
            square = term;
            exponent = term_exponent;
        } else if term_exponent > exponent {
            square = term + square * square_scale(exponent - term_exponent);
            exponent = term_exponent;
        } else {
            square += term * square_scale(term_exponent - exponent);
        }
    }
    L2Measure {
        magnitude: square.sqrt(),
        exponent,
    }
}

/// Scale of a squared term whose exponent is no greater than the accumulator's.
const fn square_scale(exponent_delta: i32) -> f64 {
    if exponent_delta < -511 {
        0.0
    } else {
        normal_power_of_two(2 * exponent_delta)
    }
}

/// Splits a positive finite number into a significand in [1, 2) and exponent.
fn split_positive(value: f64) -> (f64, i32) {
    let (normal, adjustment) = if value.is_subnormal() {
        (value * 18_014_398_509_481_984.0, -54)
    } else {
        (value, 0)
    };
    let bits = normal.to_bits();
    let mantissa = f64::from_bits((bits & 0x000f_ffff_ffff_ffff) | 0x3ff0_0000_0000_0000);
    (mantissa, ((bits >> 52) & 0x7ff) as i32 - 1023 + adjustment)
}

const fn normal_power_of_two(exponent: i32) -> f64 {
    f64::from_bits(((exponent + 1023) as u64) << 52)
}

fn scale_positive(
    value: f64,
    exponent: i32,
) -> f64 {
    if exponent == 0 || value == 0.0 {
        return value;
    }
    if (-1022..=1023).contains(&exponent) {
        return value * normal_power_of_two(exponent);
    }
    let (mantissa, current) = split_positive(value);
    let exponent = exponent + current;
    match exponent {
        | 1024.. => f64::INFINITY,
        | -1022..=1023 => mantissa * normal_power_of_two(exponent),
        | -1074..=-1023 => mantissa * f64::from_bits(1_u64 << (exponent + 1074)),
        | -1075 => (0.5 * mantissa) * f64::from_bits(1),
        | _ => 0.0,
    }
}

/// Calculates the L-infinity norm of a function represented by discrete points.
///
/// The L-infinity norm (or Chebyshev norm) is defined as `max(|f(x)|)`. For discrete points,
/// it is simply the maximum absolute value among the sampled points.
///
/// # Arguments
/// * `points` - A slice of `(x, y)` tuples representing the function's samples.
///
/// # Returns
/// The numerical value of the L-infinity norm.
pub fn infinity_norm(points: &[(f64, f64)]) -> f64 {
    points.iter().map(|(_, y)| y.abs()).fold(0.0, f64::max)
}

/// Calculates the inner product of two functions, `<f, g> = ∫f(x)g(x)dx`.
///
/// Both functions must be sampled at the same x-coordinates. For discrete points,
/// it is approximated as `Σ f(x_i)g(x_i)Δx_i`.
///
/// # Arguments
/// * `f_points` - A slice of `(x, y)` tuples representing the first function's samples.
/// * `g_points` - A slice of `(x, y)` tuples representing the second function's samples.
///
/// # Returns
/// A `Result` containing the numerical value of the inner product.
///
/// # Errors
/// Returns an error if the input functions have different numbers of sample points.
pub fn inner_product(
    f_points: &[(f64, f64)],
    g_points: &[(f64, f64)],
) -> Result<f64, String> {
    if f_points.len() != g_points.len() {
        return Err("Input functions \
                    must have the \
                    same number of \
                    sample points."
            .to_string());
    }

    if f_points
        .iter()
        .chain(g_points)
        .any(|(x, y)| !x.is_finite() || !y.is_finite())
        || f_points.windows(2).any(|w| w[0].0 >= w[1].0)
        || g_points.windows(2).any(|w| w[0].0 >= w[1].0)
        || f_points
            .iter()
            .zip(g_points)
            .any(|(f, g)| !crate::kernels::exact_eq(f.0, g.0))
    {
        return Err(
            "Inner products require finite samples on the same strictly increasing coordinate grid"
                .into(),
        );
    }

    let integral = f_points
        .windows(2)
        .enumerate()
        .map(|(i, w)| {
            let (x1, y1_f) = w[0];

            let (x2, y2_f) = w[1];

            let (_, y1_g) = g_points[i];

            let (_, y2_g) = g_points[i + 1];

            y1_f.mul_add(y1_g, y2_f * y2_g) / 2.0 * (x2 - x1)
        })
        .sum();

    Ok(integral)
}

/// Computes the projection of function `f` onto function `g`.
/// Both functions must be sampled at the same x-coordinates.
///
/// # Errors
/// Returns an error if the input functions have different numbers or inconsistent sample points.
pub fn project(
    f_points: &[(f64, f64)],
    g_points: &[(f64, f64)],
) -> Result<Vec<(f64, f64)>, String> {
    let ip_fg = inner_product(f_points, g_points)?;

    let ip_gg = inner_product(g_points, g_points)?;

    if !ip_gg.is_finite() || ip_gg < 0.0 {
        return Err("Projection requires a finite nonnegative target norm".into());
    }
    if crate::kernels::exact_eq(ip_gg, 0.0) {
        return Ok(g_points.iter().map(|&(x, _)| (x, 0.0)).collect());
    }

    let coeff = ip_fg / ip_gg;

    Ok(g_points.iter().map(|&(x, y)| (x, y * coeff)).collect())
}

/// Normalizes a function relative to its L2 norm.
///
/// Zero measure leaves the samples as supplied. Normalization uses the
/// norm's internal scale even if the rounded norm is zero or infinite.
#[must_use]
pub fn normalize(points: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let measure = l2_measure(points);
    if measure.magnitude == 0.0 {
        return points.to_vec();
    }
    if measure.exponent == 0 {
        points
            .iter()
            .map(|&(x, y)| (x, y / measure.magnitude))
            .collect()
    } else {
        points
            .iter()
            .map(|&(x, y)| (x, measure.normalize_scaled(y)))
            .collect()
    }
}

/// Performs the Gram-Schmidt process to orthogonalize a set of functions.
/// All functions must have the same sampling points.
///
/// # Errors
/// Returns an error if any pair of input functions has mismatched sample points.
pub fn gram_schmidt(basis: &[Vec<(f64, f64)>]) -> Result<Vec<Vec<(f64, f64)>>, String> {
    let mut orthogonal_basis: Vec<Vec<(f64, f64)>> = Vec::new();

    for b in basis {
        let mut v = b.clone();

        for u in &orthogonal_basis {
            let proj = project(b, u)?;

            for (j, (_, py)) in proj.into_iter().enumerate() {
                v[j].1 -= py;
            }
        }

        orthogonal_basis.push(v);
    }

    Ok(orthogonal_basis)
}

/// Performs the Gram-Schmidt process to orthonormalize a set of functions.
///
/// # Errors
/// Returns an error if any pair of input functions has mismatched sample points.
pub fn gram_schmidt_orthonormal(basis: &[Vec<(f64, f64)>]) -> Result<Vec<Vec<(f64, f64)>>, String> {
    let orthogonal_basis = gram_schmidt(basis)?;

    Ok(orthogonal_basis
        .into_iter()
        .map(|v| normalize(&v))
        .collect())
}
