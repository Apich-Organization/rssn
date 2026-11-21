//! # Robust statistics, quantiles and kernel density estimation
//!
//! Quantile definitions (Hyndman-Fan types 1 to 9), median absolute
//! deviation, trimmed/winsorised means, Huber and Tukey-biweight location
//! estimators, empirical CDF and weighted quantiles, bandwidth selectors,
//! Gaussian and Epanechnikov kernel density estimation, and seeded bootstrap
//! confidence intervals.
//!
//! # Limitations and coverage
//!
//! * Unweighted quantiles use Hyndman-Fan definitions 1–9. Weighted quantiles
//!   invert the weighted empirical CDF and ignore observations with zero weight.
//! * The location estimators are one-dimensional; the M-estimators use the
//!   normalised MAD as a fixed scale rather than estimating location and scale
//!   jointly, and there is no robust covariance, robust regression or outlier
//!   labelling.
//! * KDE supports only the Gaussian and Epanechnikov kernels with Silverman or
//!   Scott bandwidth; no adaptive-bandwidth, boundary-corrected or
//!   multivariate density is provided.
//! * Bootstrap confidence intervals use seeded IID resamples and either
//!   percentile or basic (reverse percentile) endpoints.

use std::f64::consts::PI;

use crate::kernels::descriptive;
use crate::kernels::numeric_common::scale_product;
use crate::kernels::numeric_common::scaled_difference;
use crate::kernels::par;
use crate::kernels::random::Rng;
use crate::kernels::reuse::f64_scratch;

fn sorted(data: &[f64]) -> Vec<f64> {
    let mut v: Vec<f64> = data.iter().copied().filter(|x| !x.is_nan()).collect();
    v.sort_unstable_by(f64::total_cmp);
    v
}

/// Sample quantile of `data` at probability `p` using the Hyndman-Fan
/// definition `kind` (1..=9; 7 is the NumPy/R default). NaNs are ignored.
/// Returns NaN for empty data or an invalid `p`/`kind`.
#[must_use]
pub fn quantile(
    data: &[f64],
    p: f64,
    kind: u8,
) -> f64 {
    if !(0.0..=1.0).contains(&p) || !(1..=9).contains(&kind) {
        return f64::NAN;
    }
    let mut sample: Vec<_> = data.iter().copied().filter(|x| !x.is_nan()).collect();
    quantile_select(&mut sample, p, kind)
}

/// Hyndman-Fan sample quantiles in query order, using one shared sort.
///
/// NaN observations are ignored. An empty sample or invalid definition yields
/// NaN at every query. An invalid probability yields NaN at that query.
#[must_use]
pub fn quantiles(
    data: &[f64],
    probabilities: &[f64],
    kind: u8,
) -> Vec<f64> {
    if probabilities.is_empty() {
        return Vec::new();
    }
    if !(1..=9).contains(&kind) {
        return vec![f64::NAN; probabilities.len()];
    }
    if let [p] = probabilities {
        return vec![quantile(data, *p, kind)];
    }
    let sample = sorted(data);
    probabilities
        .iter()
        .map(|&p| quantile_sorted(&sample, p, kind))
        .collect()
}

fn quantile_sorted(
    x: &[f64],
    p: f64,
    kind: u8,
) -> f64 {
    if x.is_empty() || !(0.0..=1.0).contains(&p) {
        return f64::NAN;
    }
    let (lower, upper, fraction) = quantile_position(x.len(), p, kind);
    interpolate(x[lower], x[upper], fraction)
}

pub(super) fn quantile_select(
    data: &mut [f64],
    p: f64,
    kind: u8,
) -> f64 {
    if data.is_empty() {
        return f64::NAN;
    }
    let (lower, upper, fraction) = quantile_position(data.len(), p, kind);
    let (_, &mut value, tail) = data.select_nth_unstable_by(lower, f64::total_cmp);
    if lower == upper || fraction == 0.0 {
        return value;
    }
    let next = tail.iter().copied().min_by(f64::total_cmp).unwrap_or(value);
    interpolate(value, next, fraction)
}

fn quantile_position(
    n: usize,
    p: f64,
    kind: u8,
) -> (usize, usize, f64) {
    let nf = n as f64;
    let at = |i: isize| (i.clamp(1, n as isize) - 1) as usize;
    let selected = |index| (index, index, 0.0);
    if n == 1 || p == 0.0 {
        return selected(0);
    }
    if p >= 1.0 {
        return selected(n - 1);
    }
    match kind {
        | 1 | 2 => {
            let np = nf * p;
            let j = np.floor();
            let g = np - j;
            let ji = j as isize;
            if g > 0.0 {
                selected(at(ji + 1))
            } else if kind == 1 {
                selected(at(ji))
            } else {
                (at(ji), at(ji + 1), 0.5)
            }
        },
        | 3 => {
            let h = nf * p - 0.5;
            let j = h.floor();
            if h.fract() == 0.0 && (j as isize) % 2 == 0 {
                selected(at(j as isize))
            } else {
                selected(at(j as isize + 1))
            }
        },
        | 7 => {
            let (lower, fraction) = descriptive::type7_position(n, p);
            (lower, (lower + 1).min(n - 1), fraction)
        },
        | _ => {
            let (a, b) = match kind {
                | 4 => (0.0, 1.0),
                | 5 => (0.5, 0.5),
                | 6 => (0.0, 0.0),
                | 8 => (1.0 / 3.0, 1.0 / 3.0),
                | _ => (0.375, 0.375),
            };
            // Plotting position m = a + p (1 - a - b), index h = n p + m.
            let m = a + p * (1.0 - a - b);
            let h = nf * p + m;
            let j = h.floor();
            (at(j as isize), at(j as isize + 1), h - j)
        },
    }
}

fn interpolate(
    left: f64,
    right: f64,
    fraction: f64,
) -> f64 {
    if fraction == 0.0 {
        return left;
    }
    if crate::kernels::exact_eq(left, right) {
        return f64::midpoint(left, right);
    }
    descriptive::interpolate(left, right, fraction)
}

fn weight_scale(
    data: &[f64],
    weights: &[f64],
) -> Option<f64> {
    if data.len() != weights.len() {
        return None;
    }
    let mut max_weight: f64 = 0.0;
    for (&sample, &weight) in data.iter().zip(weights) {
        if !weight.is_finite() || weight < 0.0 || (weight > 0.0 && sample.is_nan()) {
            return None;
        }
        max_weight = max_weight.max(weight);
    }
    (max_weight > 0.0).then_some(max_weight)
}

/// Smallest positive-weight observation whose weighted empirical CDF is at
/// least `p`. At `p = 0` and `p = 1` this gives the minimum and maximum
/// positive-weight observations.
///
/// Returns NaN for probabilities outside `[0, 1]`, unequal lengths, invalid
/// weights, NaN observations with positive weight or zero total weight.
#[must_use]
pub fn weighted_quantile(
    data: &[f64],
    weights: &[f64],
    p: f64,
) -> f64 {
    if !(0.0..=1.0).contains(&p) {
        return f64::NAN;
    }
    let Some(sample) = cumulative_weights(data, weights) else {
        return f64::NAN;
    };
    weighted_quantile_sorted(&sample, p)
}

fn cumulative_weights(
    data: &[f64],
    weights: &[f64],
) -> Option<Vec<(f64, f64)>> {
    let scale = weight_scale(data, weights)?;
    let mut sample: Vec<_> = data
        .iter()
        .copied()
        .zip(weights.iter().copied())
        .filter(|&(_, weight)| weight > 0.0)
        .collect();
    sample.sort_unstable_by(|(a, _), (b, _)| a.total_cmp(b));
    let mut cumulative = 0.0;
    for (_, weight) in &mut sample {
        cumulative += *weight / scale;
        *weight = cumulative;
    }
    Some(sample)
}

fn weighted_quantile_sorted(
    sample: &[(f64, f64)],
    p: f64,
) -> f64 {
    if !(0.0..=1.0).contains(&p) {
        return f64::NAN;
    }
    if p == 0.0 {
        return sample.first().map_or(f64::NAN, |&(value, _)| value);
    }
    let Some(&(last, total)) = sample.last() else {
        return f64::NAN;
    };
    if p >= 1.0 {
        return last;
    }
    let rank = sample.partition_point(|&(_, cumulative)| cumulative < p * total);
    sample.get(rank).map_or(last, |&(value, _)| value)
}

/// Weighted empirical quantiles, sharing one sort and cumulative weight scan.
///
/// Zero-weight observations are ignored. Invalid weights yield NaN at every
/// query and an invalid probability yields NaN at that query.
#[must_use]
pub fn weighted_quantiles(
    data: &[f64],
    weights: &[f64],
    probabilities: &[f64],
) -> Vec<f64> {
    if probabilities.is_empty() {
        return Vec::new();
    }
    let Some(sample) = cumulative_weights(data, weights) else {
        return vec![f64::NAN; probabilities.len()];
    };
    par::map_slice(probabilities, 1_024, |&p| {
        weighted_quantile_sorted(&sample, p)
    })
}

/// Fraction of total weight on observations at or below `x`.
///
/// Weights must be finite and nonnegative. Returns NaN for a NaN query,
/// unequal lengths, invalid weights, NaN observations with positive weight or
/// zero total weight.
#[must_use]
pub fn weighted_ecdf(
    data: &[f64],
    weights: &[f64],
    x: f64,
) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    let Some(scale) = weight_scale(data, weights) else {
        return f64::NAN;
    };
    let (mut cumulative_weight, mut total_weight) = (0.0, 0.0);
    for (&value, &weight) in data.iter().zip(weights) {
        let normalized = weight / scale;
        total_weight += normalized;
        if value <= x {
            cumulative_weight += normalized;
        }
    }
    cumulative_weight / total_weight
}

/// Median of `data` (NaNs ignored).
#[must_use]
pub fn median(data: &[f64]) -> f64 {
    let mut sample: Vec<_> = data.iter().copied().filter(|v| !v.is_nan()).collect();
    quantile_select(&mut sample, 0.5, 7)
}

/// Median absolute deviation; multiplied by 1.4826 when `normal` is true so
/// it estimates the standard deviation for Gaussian data.
#[must_use]
pub fn mad(
    data: &[f64],
    normal: bool,
) -> f64 {
    mad_about(data, median(data))
        * if normal {
            NORMAL_MAD_SCALE
        } else {
            1.0
        }
}

const NORMAL_MAD_SCALE: f64 = 1.482_602_218_505_602;

fn mad_about(
    data: &[f64],
    median_value: f64,
) -> f64 {
    let mut deviations: Vec<_> = data
        .iter()
        .map(|v| (v - median_value).abs())
        .filter(|v| !v.is_nan())
        .collect();
    quantile_select(&mut deviations, 0.5, 7)
}

/// Interquartile range (type 7 quantiles).
#[must_use]
pub fn iqr(data: &[f64]) -> f64 {
    let mut sample: Vec<_> = data.iter().copied().filter(|v| !v.is_nan()).collect();
    let (lower, upper) = quartiles(&mut sample);
    upper - lower
}

fn quartiles(sample: &mut [f64]) -> (f64, f64) {
    (
        quantile_select(sample, 0.25, 7),
        quantile_select(sample, 0.75, 7),
    )
}

/// Mean after discarding the fraction `trim` from each tail.
#[must_use]
pub fn trimmed_mean(
    data: &[f64],
    trim: f64,
) -> f64 {
    let x = sorted(data);
    if x.is_empty() || trim.is_nan() {
        return f64::NAN;
    }
    let k = (x.len() as f64 * trim.clamp(0.0, 0.499)).floor() as usize;
    descriptive::mean(&x[k..x.len() - k]).unwrap_or(f64::NAN)
}

/// Mean after replacing the lowest and highest fraction `trim` by the
/// nearest retained values.
#[must_use]
pub fn winsorized_mean(
    data: &[f64],
    trim: f64,
) -> f64 {
    let mut x = sorted(data);
    let n = x.len();
    let k = (n as f64 * trim.clamp(0.0, 0.499)).floor() as usize;
    if n == 0 || trim.is_nan() {
        return f64::NAN;
    }
    let (lo, hi) = (x[k], x[n - 1 - k]);
    for v in &mut x {
        *v = v.clamp(lo, hi);
    }
    descriptive::mean(&x).unwrap_or(f64::NAN)
}

/// Huber M-estimator of location with tuning constant `k` (1.345 is the
/// classical choice), scale fixed to the normalised MAD.
///
/// NaNs are ignored. An invalid tuning constant or an empty effective sample
/// yields NaN. An undefined or zero MAD gives the sample median.
#[must_use]
pub fn huber_location(
    data: &[f64],
    k: f64,
) -> f64 {
    if !k.is_finite() || k <= 0.0 {
        return f64::NAN;
    }
    let (mut samples, center) = huber_sample(data);
    if !center.is_finite() {
        return center;
    }
    for value in &mut samples {
        *value = (*value - center).abs();
    }
    let scale = quantile_select(&mut samples, 0.5, 7);
    if scale == 0.0 || !scale.is_finite() {
        return center;
    }
    huber_fit(data, &mut samples, center, [scale, NORMAL_MAD_SCALE, k])
}

/// Huber location using a caller-supplied scale in the units of the sample.
///
/// The clipping distance is `k * scale`. NaNs are ignored. The tuning
/// constant and scale must both be positive and finite. Returns NaN for
/// invalid parameters or an empty effective sample. For a flat interval
/// containing the median, returns the median.
#[must_use]
pub fn huber_location_with_scale(
    data: &[f64],
    k: f64,
    scale: f64,
) -> f64 {
    if !k.is_finite() || k <= 0.0 || !scale.is_finite() || scale <= 0.0 {
        return f64::NAN;
    }
    let (mut samples, center) = huber_sample(data);
    if !center.is_finite() {
        return center;
    }
    huber_fit(data, &mut samples, center, [scale, 1.0, k])
}

fn huber_sample(data: &[f64]) -> (Vec<f64>, f64) {
    let observations = data.iter().copied().filter(|v| !v.is_nan());
    let mut sample = Vec::with_capacity(observations.clone().count());
    sample.extend(observations);
    let center = quantile_select(&mut sample, 0.5, 7);
    (sample, center)
}

#[allow(
    clippy::float_cmp,
    reason = "An exact zero displacement preserves the original median."
)]
fn huber_fit(
    data: &[f64],
    samples: &mut [f64],
    center: f64,
    factors: [f64; 3],
) -> f64 {
    let scale = factors[0] * factors[1];
    let direct_cutoff = scale * factors[2];
    let cutoff = if scale.is_normal() && direct_cutoff.is_normal() {
        direct_cutoff
    } else {
        scale_product(1.0, &factors, &[]).unwrap_or(f64::INFINITY)
    };
    let observations = data.iter().copied().filter(|value| !value.is_nan());
    let mut fully_linear = true;
    let mut infinite_balance = 0_i64;
    for (normalized, value) in samples.iter_mut().zip(observations) {
        *normalized = if !value.is_finite() {
            infinite_balance += if value.is_sign_positive() {
                1
            } else {
                -1
            };
            value
        } else if cutoff.is_normal() {
            scaled_difference(value, center, cutoff)
        } else {
            let difference = value - center;
            let (numerator, scale) = if difference.is_finite() {
                (difference, 1.0)
            } else {
                (0.5 * value - 0.5 * center, 2.0)
            };
            scale_product(numerator, &[scale], &factors)
                .unwrap_or_else(|_| f64::INFINITY.copysign(numerator))
        };
        if value.is_finite() {
            fully_linear &= normalized.abs() <= 0.5;
        }
    }
    // A finite sample contained in one clipping window stays in the linear
    // part at its mean. Work in physical units if normalizing would lose
    // small observations. Opposite infinite observations cancel in the score.
    if fully_linear && infinite_balance == 0 {
        let mut count = 0;
        for (&value, slot) in data
            .iter()
            .filter(|value| value.is_finite())
            .zip(samples.iter_mut())
        {
            *slot = value;
            count += 1;
        }
        return descriptive::mean(&samples[..count]).unwrap_or(center);
    }

    let root = huber_offset(samples);
    if root == 0.0 {
        return center;
    }
    if cutoff.is_normal() {
        return center + cutoff * root;
    }
    if let Ok(displacement) = scale_product(root, &factors, &[]) {
        return center + displacement;
    }
    let half_displacement =
        scale_product(root, &factors, &[2.0]).unwrap_or_else(|_| f64::INFINITY.copysign(root));
    2.0 * (0.5 * center + half_displacement)
}

// Half the sample is on each side of the median. Moving by one clipping
// distance makes that half contribute its full signed weight, so the score
// has opposite weak signs at -1 and 1 in these coordinates.
#[allow(
    clippy::float_cmp,
    reason = "Detect exhaustion of representable points in the bracket."
)]
fn huber_offset(samples: &[f64]) -> f64 {
    let (mut lower, mut upper) = (-1.0, 1.0);
    let mut x = 0.0;
    for _ in 0..200 {
        let mut sum = 0.0;
        let mut clipped = 0_i64;
        let mut active = 0_usize;
        for &sample in samples {
            let difference = sample - x;
            if difference > 1.0 {
                clipped += 1;
            } else if difference < -1.0 {
                clipped -= 1;
            } else {
                sum += difference;
                active += 1;
            }
        }
        let score = sum + clipped as f64;
        if score == 0.0 {
            return x;
        }
        if score > 0.0 {
            lower = x;
        } else {
            upper = x;
        }
        let correction = score / active as f64;
        let newton = x + correction;
        if newton == x {
            return x;
        }
        let next = if lower < newton && newton < upper {
            // Keep the correction while avoiding repeated scans at the
            // rounding floor of the normalized coordinate.
            if correction.abs() <= 8.0 * f64::EPSILON {
                return newton;
            }
            newton
        } else {
            f64::midpoint(lower, upper)
        };
        if next == lower || next == upper {
            return next;
        }
        x = next;
    }
    x
}

/// Tukey biweight location estimator (tuning constant `c`, 4.685 typical).
#[must_use]
pub fn biweight_location(
    data: &[f64],
    c: f64,
) -> f64 {
    if !c.is_finite() || c <= 0.0 {
        return f64::NAN;
    }
    let mut mu = median(data);
    let s = NORMAL_MAD_SCALE * mad_about(data, mu);
    if s == 0.0 || !s.is_finite() {
        return mu;
    }
    for _ in 0..200 {
        let (mut num, mut den) = (0.0, 0.0);
        for &v in data {
            let u = (v - mu) / (c * s);
            if u.abs() < 1.0 {
                let w = (1.0 - u * u).powi(2);
                num += w * (v - mu);
                den += w;
            }
        }
        if den == 0.0 {
            break;
        }
        let d = num / den;
        mu += d;
        if d.abs() < 1e-13 * s {
            break;
        }
    }
    mu
}

/// Empirical CDF of `data` at `x`.
#[must_use]
pub fn ecdf(
    data: &[f64],
    x: f64,
) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    let (mut at_or_below, mut total) = (0_usize, 0_usize);
    for &value in data {
        if !value.is_nan() {
            total += 1;
            at_or_below += usize::from(value <= x);
        }
    }
    at_or_below as f64 / total as f64
}

/// Empirical CDF at several query points, in the supplied order.
///
/// NaN observations are ignored and NaN queries yield NaN. An empty effective
/// sample yields NaN at every query. The ordered sample is prepared once and
/// each query uses a binary search.
#[must_use]
pub fn ecdf_batch(
    data: &[f64],
    queries: &[f64],
) -> Vec<f64> {
    if queries.is_empty() {
        return Vec::new();
    }
    if let [query] = queries {
        return vec![ecdf(data, *query)];
    }
    let sample = sorted(data);
    par::map_slice(queries, 1_024, |&x| {
        if x.is_nan() {
            return f64::NAN;
        }
        let count = sample.partition_point(|&value| value <= x);
        count as f64 / sample.len() as f64
    })
}

/// Silverman's rule-of-thumb bandwidth `0.9 min(sd, IQR/1.34) n^(-1/5)`.
/// Uses `sd` when `IQR` is zero. Requires at least two finite observations.
#[must_use]
pub fn bandwidth_silverman(data: &[f64]) -> f64 {
    let Ok((sd, scale)) = descriptive::std_dev_parts(data, false) else {
        return f64::NAN;
    };
    if data.len() < 2 {
        return f64::NAN;
    }
    if scale == 0.0 {
        return 0.0;
    }
    let (lower, upper) = quartiles(&mut data.to_vec());
    let iqr = scaled_difference(upper, lower, scale);
    let spread = if iqr > 0.0 {
        sd.min(iqr / 1.34)
    } else {
        sd
    };
    scale * (0.9 * (data.len() as f64).powf(-0.2) * spread)
}

/// Scott's rule bandwidth `1.06 sd n^(-1/5)`.
/// Requires at least two finite observations.
#[must_use]
pub fn bandwidth_scott(data: &[f64]) -> f64 {
    let Ok((sd, scale)) = descriptive::std_dev_parts(data, false) else {
        return f64::NAN;
    };
    if data.len() < 2 {
        return f64::NAN;
    }
    scale * (1.06 * (data.len() as f64).powf(-0.2) * sd)
}

/// Kernel shape for [`kde`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kernel {
    /// Standard normal kernel.
    Gaussian,
    /// Epanechnikov kernel (compact support).
    Epanechnikov,
}

/// Kernel density estimate of `data` evaluated at `x` with bandwidth `h`.
#[must_use]
pub fn kde(
    data: &[f64],
    x: f64,
    h: f64,
    kernel: Kernel,
) -> f64 {
    if data.is_empty() || !h.is_finite() || h <= 0.0 || x.is_nan() {
        return f64::NAN;
    }
    let mut count = 0_usize;
    let s: f64 = data
        .iter()
        .filter(|d| !d.is_nan())
        .map(|&d| {
            count += 1;
            kernel_weight(scaled_difference(x, d, h), kernel)
        })
        .sum();
    (s / count as f64) / h
}

fn kernel_weight(
    u: f64,
    kernel: Kernel,
) -> f64 {
    match kernel {
        | Kernel::Gaussian => (-0.5 * u * u).exp() / (2.0 * PI).sqrt(),
        | Kernel::Epanechnikov if u.abs() < 1.0 => 0.75 * (1.0 - u * u),
        | Kernel::Epanechnikov => 0.0,
    }
}

pub(super) fn gaussian_kde_finite(
    data: &[f64],
    points: &[f64],
    h: f64,
) -> Vec<f64> {
    let grain = 32_768usize.div_ceil(data.len()).max(2);
    let count = data.len() as f64;
    par::map_slice(points, grain, |&point| {
        let sum: f64 = data
            .iter()
            .map(|&value| kernel_weight(scaled_difference(point, value, h), Kernel::Gaussian))
            .sum();
        (sum / count) / h
    })
}

/// Percentile bootstrap confidence interval for `statistic` at level
/// `1 - alpha`, with `resamples` seeded resamples.
///
/// Calls `statistic` once per
/// resample in draw order and excludes NaN replicates from the interval.
#[must_use]
pub fn bootstrap_ci<S: Fn(&[f64]) -> f64>(
    data: &[f64],
    statistic: S,
    resamples: usize,
    alpha: f64,
    seed: u64,
) -> (f64, f64) {
    bootstrap_serial(data, |sample| statistic(sample), resamples, alpha, seed)
}

fn bootstrap_serial(
    data: &[f64],
    statistic: impl Fn(&mut [f64]) -> f64,
    resamples: usize,
    alpha: f64,
    seed: u64,
) -> (f64, f64) {
    let Some(mut stats) = bootstrap_output(data, resamples, alpha) else {
        return (f64::NAN, f64::NAN);
    };
    let mut rng = Rng::new(seed);
    let mut sample = f64_scratch(data.len());
    sample.resize(data.len(), 0.0);
    for _ in 0..resamples {
        fill_resample(data, &mut rng, &mut sample);
        let value = statistic(&mut sample);
        if !value.is_nan() {
            stats.push(value);
        }
    }
    bootstrap_percentiles(&mut stats, alpha)
}

/// Percentile bootstrap interval with a statistic that can run concurrently.
///
/// Uses the same seed and samples as [`bootstrap_ci`]. For sufficiently large
/// workloads, evaluates the statistic on a bounded batch in the active Rayon
/// pool. A statistic that depends only on its input produces the same interval
/// in every pool size. NaN replicates are excluded.
#[must_use]
pub fn bootstrap_ci_parallel<S: Fn(&[f64]) -> f64 + Sync>(
    data: &[f64],
    statistic: S,
    resamples: usize,
    alpha: f64,
    seed: u64,
) -> (f64, f64) {
    bootstrap_sampled(data, |sample| statistic(sample), resamples, alpha, seed)
}

fn bootstrap_sampled(
    data: &[f64],
    statistic: impl Fn(&mut [f64]) -> f64 + Sync,
    resamples: usize,
    alpha: f64,
    seed: u64,
) -> (f64, f64) {
    const BATCH_DRAWS: usize = 262_144;
    const MIN_PARALLEL_DRAWS: usize = 65_536;
    let width = data.len();
    let batch_rows = BATCH_DRAWS.checked_div(width).unwrap_or(0).min(resamples);
    let min_rows = MIN_PARALLEL_DRAWS.div_ceil(width.max(1)).max(2);
    if !crate::execution::should_split(batch_rows, min_rows) {
        return bootstrap_serial(data, statistic, resamples, alpha, seed);
    }
    let Some(mut stats) = bootstrap_output(data, resamples, alpha) else {
        return (f64::NAN, f64::NAN);
    };
    let mut rng = Rng::new(seed);
    let capacity = batch_rows.saturating_mul(width);
    let mut samples = f64_scratch(capacity);
    samples.resize(capacity, 0.0);
    for start in (0..resamples).step_by(batch_rows) {
        let rows = batch_rows.min(resamples.saturating_sub(start));
        samples.truncate(rows.saturating_mul(width));
        fill_resample(data, &mut rng, &mut samples);
        par::for_chunks_mut(&mut samples, width, min_rows, |_, row| {
            let value = statistic(row);
            if let Some(first) = row.first_mut() {
                *first = value;
            }
        });
        stats.extend(
            samples
                .chunks_exact(width)
                .filter_map(|row| row.first().copied().filter(|value| !value.is_nan())),
        );
    }
    bootstrap_percentiles(&mut stats, alpha)
}

#[derive(Clone, Copy)]
pub(crate) enum BootstrapOrderStatistic {
    Median,
    Mad,
    Iqr,
}

impl BootstrapOrderStatistic {
    /// Uses mutable resamples directly. The bridge has already checked the input is finite.
    pub(crate) fn interval_finite(
        self,
        data: &[f64],
        resamples: usize,
        alpha: f64,
        seed: u64,
        basic: bool,
    ) -> (f64, f64) {
        let interval =
            bootstrap_sampled(data, |sample| self.evaluate(sample), resamples, alpha, seed);
        if !basic || !valid_bootstrap(data, resamples, alpha) {
            return interval;
        }
        let mut observed = f64_scratch(data.len());
        observed.extend_from_slice(data);
        reflect_bootstrap(self.evaluate(&mut observed), interval)
    }

    fn evaluate(
        self,
        sample: &mut [f64],
    ) -> f64 {
        match self {
            | Self::Median => quantile_select(sample, 0.5, 7),
            | Self::Mad => {
                let middle = quantile_select(sample, 0.5, 7);
                for value in &mut *sample {
                    *value = (*value - middle).abs();
                }
                quantile_select(sample, 0.5, 7)
            },
            | Self::Iqr => {
                let (lower, upper) = quartiles(sample);
                upper - lower
            },
        }
    }
}

/// Basic (reverse percentile) bootstrap confidence interval.
///
/// Reflects the seeded percentile endpoints about the statistic evaluated on
/// the observed sample. The statistic has the concurrency contract of
/// [`bootstrap_ci_parallel`]. Returns NaN endpoints if that estimate is not
/// finite. NaN replicates are excluded.
#[must_use]
pub fn bootstrap_basic_ci<S: Fn(&[f64]) -> f64 + Sync>(
    data: &[f64],
    statistic: S,
    resamples: usize,
    alpha: f64,
    seed: u64,
) -> (f64, f64) {
    if !valid_bootstrap(data, resamples, alpha) {
        return (f64::NAN, f64::NAN);
    }
    let estimate = statistic(data);
    if !estimate.is_finite() {
        return (f64::NAN, f64::NAN);
    }
    let interval = bootstrap_ci_parallel(data, statistic, resamples, alpha, seed);
    reflect_bootstrap(estimate, interval)
}

fn reflect_bootstrap(
    estimate: f64,
    (lower, upper): (f64, f64),
) -> (f64, f64) {
    if estimate.is_finite() {
        (estimate + (estimate - upper), estimate + (estimate - lower))
    } else {
        (f64::NAN, f64::NAN)
    }
}

fn valid_bootstrap(
    data: &[f64],
    resamples: usize,
    alpha: f64,
) -> bool {
    !data.is_empty() && resamples != 0 && (0.0..=1.0).contains(&alpha)
}

fn bootstrap_output(
    data: &[f64],
    resamples: usize,
    alpha: f64,
) -> Option<Vec<f64>> {
    if !valid_bootstrap(data, resamples, alpha) {
        return None;
    }
    let mut stats = Vec::new();
    stats.try_reserve_exact(resamples).ok()?;
    Some(stats)
}

fn fill_resample(
    data: &[f64],
    rng: &mut Rng,
    output: &mut [f64],
) {
    let len = data.len() as u64;
    for value in output {
        *value = data[rng.below(len) as usize];
    }
}

fn bootstrap_percentiles(
    stats: &mut [f64],
    alpha: f64,
) -> (f64, f64) {
    (
        quantile_select(stats, alpha / 2.0, 7),
        quantile_select(stats, 1.0 - alpha / 2.0, 7),
    )
}
