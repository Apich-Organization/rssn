//! Local PCHIP and Akima Hermite segments.

#![allow(
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    reason = "stencil and query indices come from checked slices; finite results are checked at the boundary"
)]

use super::InterpError;
use super::find_segment;
use super::validate_1d;
use crate::kernels::exact_eq;
use crate::kernels::numeric_common::scale_product;
use crate::kernels::numeric_common::scaled_difference;
use crate::kernels::par;

#[derive(Clone, Copy)]
pub(super) enum Method {
    Pchip,
    Akima,
}

const BATCH_CHUNK: usize = 8192;

pub(super) fn evaluate(
    xs: &[f64],
    ys: &[f64],
    query: f64,
    method: Method,
) -> Result<f64, InterpError> {
    check_queries(&[query])?;
    validate_1d(xs, ys)?;
    if xs.len() == 1 {
        return Ok(ys[0]);
    }
    let index = find_segment(xs, query);
    let value = Segment::prepare(xs, ys, index, method).value(query);
    if value.is_finite() {
        Ok(value)
    } else {
        Err(unrepresentable())
    }
}

pub(super) fn evaluate_many(
    xs: &[f64],
    ys: &[f64],
    queries: &[f64],
    method: Method,
) -> Result<Vec<f64>, InterpError> {
    validate_1d(xs, ys)?;
    check_queries(queries)?;
    let mut output = Vec::new();
    output.try_reserve_exact(queries.len()).map_err(|_| {
        InterpError::Domain("interpolation batch output cannot be allocated".into())
    })?;
    if xs.len() == 1 {
        output.resize(queries.len(), ys[0]);
        return Ok(output);
    }
    output.resize(queries.len(), 0.0);
    par::for_chunks_mut(&mut output, BATCH_CHUNK, 8, |batch, values| {
        let start = batch * BATCH_CHUNK;
        let queries = &queries[start..start + values.len()];
        let mut index = find_segment(xs, queries[0]);
        let mut segment = Segment::prepare(xs, ys, index, method);
        for (&query, value) in queries.iter().zip(values) {
            let next = scan_from(xs, index, query);
            if index != next {
                index = next;
                segment = Segment::prepare(xs, ys, index, method);
            }
            *value = segment.value(query);
        }
    });
    if output.iter().all(|value| value.is_finite()) {
        Ok(output)
    } else {
        Err(unrepresentable())
    }
}

fn check_queries(queries: &[f64]) -> Result<(), InterpError> {
    if queries.iter().any(|query| !query.is_finite()) {
        Err(InterpError::Domain("query must be finite".into()))
    } else {
        Ok(())
    }
}

fn unrepresentable() -> InterpError {
    InterpError::Domain("local interpolation arithmetic exceeds the Float range".into())
}

// An ordered sweep normally moves no further than the next segment. Bound
// the forward scan so sparse or unordered queries still take logarithmic time.
fn scan_from(
    xs: &[f64],
    index: usize,
    query: f64,
) -> usize {
    let last = xs.len() - 2;
    if index > 0 && query < xs[index] {
        return find_segment(xs, query);
    }
    for candidate in index..(index + 4).min(last) {
        if query < xs[candidate + 1] {
            return candidate;
        }
    }
    let candidate = (index + 4).min(last);
    if candidate == last || query < xs[candidate + 1] {
        candidate
    } else {
        find_segment(xs, query)
    }
}

// Powers of two preserve the significant bits in local differences. Only
// rescale atypical magnitudes to keep ordinary segments on the direct path.
fn binary_scale(maximum: f64) -> f64 {
    if ordinary_magnitude(maximum) {
        1.0
    } else {
        f64::from_bits(maximum.to_bits() & 0x7ff0_0000_0000_0000).max(f64::MIN_POSITIVE)
    }
}

fn ordinary_magnitude(value: f64) -> bool {
    value == 0.0 || (1e-100..=1e100).contains(&value.abs())
}

/// Fits the cubic coefficients for ordinary scales in one pass through the
/// shared knot tangents. Returns the last tangent in segment coordinates.
#[allow(
    clippy::single_call_fn,
    reason = "keeps block scheduling with the shared tangent kernel"
)]
pub(in crate::kernels) fn pchip_polynomials(
    xs: &[f64],
    ys: &[f64],
) -> Option<(Vec<[f64; 3]>, f64)> {
    const BLOCK: usize = 1024;
    if xs.len() < 2 || xs.len() != ys.len() {
        return None;
    }
    let segments = xs.len() - 1;
    let count = segments.div_ceil(BLOCK);
    if !crate::execution::should_split(count, 4) {
        return pchip_polynomial_block(xs, ys);
    }
    // A block gets one neighboring span on either side so a boundary
    // between workers still uses the same interior knot tangent.
    let blocks = par::try_map_indexed(count, 4, |block| {
        let start = block * BLOCK;
        let end = (start + BLOCK).min(segments);
        let first = start.saturating_sub(1);
        let last = (end + 1).min(segments);
        let (polynomials, tangent) =
            pchip_polynomial_block(&xs[first..=last], &ys[first..=last]).ok_or(())?;
        Ok::<_, ()>((polynomials, tangent, start - first, end - start))
    })
    .ok()?;
    let last_tangent = blocks.last()?.1;
    let mut polynomials = Vec::with_capacity(segments);
    for (block, _, skip, keep) in blocks {
        polynomials.extend(block.into_iter().skip(skip).take(keep));
    }
    Some((polynomials, last_tangent))
}

fn pchip_polynomial_block(
    xs: &[f64],
    ys: &[f64],
) -> Option<(Vec<[f64; 3]>, f64)> {
    let mut points = xs.iter().copied().zip(ys.iter().copied());
    let mut previous = points.next()?;
    if !ordinary_magnitude(previous.0) || !ordinary_magnitude(previous.1) {
        return None;
    }
    let mut increments = points.map(|point| {
        let h = point.0 - previous.0;
        let rise = point.1 - previous.1;
        previous = point;
        let slope = rise / h;
        (ordinary_magnitude(point.0)
            && ordinary_magnitude(point.1)
            && h > 0.0
            && ordinary_magnitude(slope))
        .then_some((Secant { h, slope }, rise))
    });
    let (mut secant, mut rise) = increments.next()??;
    let Some(second) = increments.next() else {
        return Some((vec![[rise, 0.0, 0.0]], rise));
    };
    let (mut next, mut next_rise) = second?;
    let mut left = pchip_endpoint(secant, next);
    let coefficients = |rise: f64, left: f64, right: f64| {
        if exact_eq(left, rise) && exact_eq(right, rise) {
            [left, 0.0, 0.0]
        } else {
            [
                left,
                3.0 * rise - 2.0 * left - right,
                left + right - 2.0 * rise,
            ]
        }
    };
    let mut polynomials = Vec::with_capacity(xs.len() - 1);
    loop {
        let right = pchip_interior(secant, next);
        polynomials.push(coefficients(rise, left * secant.h, right * secant.h));
        let previous = secant;
        secant = next;
        rise = next_rise;
        left = right;
        if let Some(value) = increments.next() {
            (next, next_rise) = value?;
        } else {
            let last = pchip_endpoint(secant, previous) * secant.h;
            polynomials.push(coefficients(rise, left * secant.h, last));
            return Some((polynomials, last));
        }
    }
}

fn coordinate_scale(points: &[f64]) -> f64 {
    binary_scale(points.iter().fold(0.0_f64, |max, x| max.max(x.abs())))
}

#[derive(Clone, Copy)]
struct Secant {
    h: f64,
    slope: f64,
}

#[derive(Clone, Debug)]
pub(in crate::kernels) struct Segment {
    x: [f64; 2],
    y: [f64; 2],
    x_scale: f64,
    y_scale: f64,
    width: f64,
    tangent: [f64; 2],
    rise: f64,
}

impl Segment {
    pub(in crate::kernels) fn pchip(
        xs: &[f64],
        ys: &[f64],
        index: usize,
    ) -> Self {
        Self::prepare(xs, ys, index, Method::Pchip)
    }

    pub(in crate::kernels) const fn left_coordinate(&self) -> f64 {
        self.x[0]
    }

    fn prepare(
        xs: &[f64],
        ys: &[f64],
        index: usize,
        method: Method,
    ) -> Self {
        let radius = match method {
            | Method::Pchip => 1,
            | Method::Akima => 2,
        };
        let first = index.saturating_sub(radius);
        let end = (index + radius + 2).min(xs.len());
        let x_scale = coordinate_scale(&xs[first..end]);
        let y_scale = coordinate_scale(&ys[first..end]);
        let x = [xs[index], xs[index + 1]];
        let y = [ys[index], ys[index + 1]];
        let width = scaled_difference(x[1], x[0], x_scale);
        let rise = scaled_difference(y[1], y[0], y_scale);
        let tangent = if xs.len() == 2 {
            [rise, rise]
        } else {
            let mut secants = [Secant { h: 0.0, slope: 0.0 }; 5];
            let local = &mut secants[..end - first - 1];
            for (j, value) in (first..end - 1).zip(local.iter_mut()) {
                let h = scaled_difference(xs[j + 1], xs[j], x_scale);
                let dy = scaled_difference(ys[j + 1], ys[j], y_scale);
                *value = Secant { h, slope: dy / h };
            }
            let scale = binary_scale(local.iter().fold(0.0_f64, |max, s| max.max(s.slope.abs())));
            if local.iter().all(|s| s.h > 0.0 && s.slope.is_finite()) {
                for s in local.iter_mut() {
                    s.slope /= scale;
                }
                let slopes = match method {
                    | Method::Pchip => pchip_tangents(local, index - first),
                    | Method::Akima => akima_tangents(local, index - first),
                };
                slopes.map(|slope| {
                    let direct = slope * width * scale;
                    if direct.is_finite() && (direct != 0.0 || slope == 0.0) {
                        direct
                    } else {
                        scale_product(slope, &[width, scale], &[]).unwrap_or(f64::NAN)
                    }
                })
            } else {
                [f64::NAN; 2]
            }
        };
        Self {
            x,
            y,
            x_scale,
            y_scale,
            width,
            tangent,
            rise,
        }
    }

    #[inline]
    pub(in crate::kernels) fn value(
        &self,
        query: f64,
    ) -> f64 {
        if exact_eq(query, self.x[0]) {
            return self.y[0];
        }
        if exact_eq(query, self.x[1]) {
            return self.y[1];
        }
        let [left, right] = self.tangent;
        if exact_eq(self.y[0], self.y[1]) && left == 0.0 && right == 0.0 {
            return self.y[0];
        }
        let t = scaled_difference(query, self.x[0], self.x_scale) / self.width;
        let y0 = self.y[0] / self.y_scale;
        if exact_eq(left, self.rise) && exact_eq(right, self.rise) {
            let direct = (y0 + t * self.rise) * self.y_scale;
            return if direct.is_finite() && (t != 0.0 || self.rise == 0.0) {
                direct
            } else {
                self.extended_line(query)
            };
        }
        let u = 1.0 - t;
        let y1 = self.y[1] / self.y_scale;
        ((1.0 + 2.0 * t) * u * u * y0 + (3.0 - 2.0 * t) * t * t * y1 + t * u * u * left
            - t * t * u * right)
            * self.y_scale
    }

    #[inline]
    pub(in crate::kernels) fn derivative(
        &self,
        query: f64,
    ) -> f64 {
        let [left, right] = self.tangent;
        let scaled = if exact_eq(left, self.rise) && exact_eq(right, self.rise) {
            self.rise
        } else if exact_eq(query, self.x[0]) {
            left
        } else if exact_eq(query, self.x[1]) {
            right
        } else {
            let t = scaled_difference(query, self.x[0], self.x_scale) / self.width;
            let u = 1.0 - t;
            6.0 * t * u * self.rise + u * (1.0 - 3.0 * t) * left + t * (3.0 * t - 2.0) * right
        };
        if scaled == 0.0 {
            return scaled;
        }
        if exact_eq(self.x_scale, 1.0) && exact_eq(self.y_scale, 1.0) {
            return scaled / self.width;
        }
        let ordinate = scaled * self.y_scale;
        let relative = ordinate / self.width;
        let direct = relative / self.x_scale;
        if ordinate.is_normal() && relative.is_normal() && direct.is_finite() {
            direct
        } else {
            scale_product(scaled, &[self.y_scale], &[self.width, self.x_scale]).unwrap_or(f64::NAN)
        }
    }

    fn extended_line(
        &self,
        query: f64,
    ) -> f64 {
        // The relative coordinate can exceed the Float range even when the
        // product of the relative coordinate and a small ordinate change fits.
        let difference = query - self.x[0];
        let shift = if difference.is_finite() {
            1.0
        } else {
            coordinate_scale(&[query, self.x[0]])
        };
        let distance = scaled_difference(query, self.x[0], shift);
        let rise = if distance.is_sign_negative() {
            -self.rise
        } else {
            self.rise
        };
        let delta = scale_product(
            rise,
            &[distance.abs(), shift, self.y_scale],
            &[self.width, self.x_scale],
        );
        self.y[0] + delta.unwrap_or(f64::NAN)
    }
}

fn pchip_tangents(
    secants: &[Secant],
    current: usize,
) -> [f64; 2] {
    let center = secants[current];
    // A cubic stencil always contains a real secant on at least one side.
    if current == 0 {
        let next = secants[current + 1];
        [pchip_endpoint(center, next), pchip_interior(center, next)]
    } else if current == secants.len() - 1 {
        let previous = secants[current - 1];
        [
            pchip_interior(previous, center),
            pchip_endpoint(center, previous),
        ]
    } else {
        [
            pchip_interior(secants[current - 1], center),
            pchip_interior(center, secants[current + 1]),
        ]
    }
}

fn relative_width(
    left: f64,
    right: f64,
) -> f64 {
    let sum = left + right;
    if sum.is_finite() {
        return left / sum;
    }
    let max = left.max(right);
    let a = left / max;
    let b = right / max;
    a / (a + b)
}

fn pchip_interior(
    left: Secant,
    right: Secant,
) -> f64 {
    let a = left.slope;
    let b = right.slope;
    if a == 0.0 || b == 0.0 || a.is_sign_negative() != b.is_sign_negative() {
        return 0.0;
    }
    let fraction = relative_width(left.h, right.h);
    let w_left = (2.0 - fraction) / 3.0;
    let w_right = (1.0 + fraction) / 3.0;
    // Divide by the larger slope first, keeping the harmonic mean finite
    // even when computing either reciprocal individually would overflow.
    if a.abs() <= b.abs() {
        a / (w_left + w_right * (a / b))
    } else {
        b / (w_right + w_left * (b / a))
    }
}

fn pchip_endpoint(
    edge: Secant,
    neighbor: Secant,
) -> f64 {
    let fraction = relative_width(edge.h, neighbor.h);
    let slope = (1.0 + fraction) * edge.slope - fraction * neighbor.slope;
    if edge.slope == 0.0 || edge.slope.is_sign_negative() != slope.is_sign_negative() {
        0.0
    } else if edge.slope.is_sign_negative() != neighbor.slope.is_sign_negative()
        && slope.abs() > 3.0 * edge.slope.abs()
    {
        3.0 * edge.slope
    } else {
        slope
    }
}

fn akima_tangents(
    secants: &[Secant],
    current: usize,
) -> [f64; 2] {
    let center = secants[current].slope;
    let left = if current == 0 {
        2.0 * center - secants[1].slope
    } else {
        secants[current - 1].slope
    };
    let right = secants
        .get(current + 1)
        .map_or(2.0 * center - left, |s| s.slope);
    let leftmost = if current >= 2 {
        secants[current - 2].slope
    } else {
        2.0 * left - center
    };
    let rightmost = secants
        .get(current + 2)
        .map_or(2.0 * right - center, |s| s.slope);
    [
        akima_weighted(leftmost, left, center, right),
        akima_weighted(left, center, right, rightmost),
    ]
}

fn akima_weighted(
    a: f64,
    b: f64,
    c: f64,
    d: f64,
) -> f64 {
    let w1 = (d - c).abs();
    let w2 = (b - a).abs();
    let total = w1 + w2;
    if total > 0.0 {
        (w1 / total) * b + (w2 / total) * c
    } else {
        f64::midpoint(b, c)
    }
}
