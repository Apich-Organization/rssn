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
//! * Bootstrap intervals are percentile only (`1 - alpha`, seeded resamples);
//!   there is no `BCa`, studentised or block bootstrap, and the resampling is
//!   single-threaded.

use std::f64::consts::PI;

use crate::kernels::random::Rng;

fn sorted(data: &[f64]) -> Vec<f64> {
    let mut v: Vec<f64> = data.iter().copied().filter(|x| !x.is_nan()).collect();
    v.sort_by(f64::total_cmp);
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
    quantile_sorted(&sorted(data), p, kind)
}

fn quantile_sorted(
    x: &[f64],
    p: f64,
    kind: u8,
) -> f64 {
    let Some(&first) = x.first() else {
        return f64::NAN;
    };
    let n = x.len();
    let nf = n as f64;
    let at = |i: isize| x[(i.clamp(1, n as isize) - 1) as usize];
    if n == 1 || p == 0.0 {
        return first;
    }
    if p >= 1.0 {
        return at(n as isize);
    }
    match kind {
        | 1 | 2 => {
            let np = nf * p;
            let j = np.floor();
            let g = np - j;
            let ji = j as isize;
            if g > 0.0 {
                at(ji + 1)
            } else if kind == 1 {
                at(ji)
            } else {
                interpolate(at(ji), at(ji + 1), 0.5)
            }
        },
        | 3 => {
            let h = nf * p - 0.5;
            let j = h.floor();
            if h.fract() == 0.0 && (j as isize) % 2 == 0 {
                at(j as isize)
            } else {
                at(j as isize + 1)
            }
        },
        | _ => {
            let (a, b) = match kind {
                | 4 => (0.0, 1.0),
                | 5 => (0.5, 0.5),
                | 6 => (0.0, 0.0),
                | 7 => (1.0, 1.0),
                | 8 => (1.0 / 3.0, 1.0 / 3.0),
                | _ => (0.375, 0.375),
            };
            // Plotting position m = a + p (1 - a - b), index h = n p + m.
            let m = a + p * (1.0 - a - b);
            let h = nf * p + m;
            let j = h.floor();
            interpolate(at(j as isize), at(j as isize + 1), h - j)
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
    let width = right - left;
    if width.is_finite() {
        left + fraction * width
    } else {
        // Avoid an overflowing difference and cancellation such as -inf + inf.
        (1.0 - fraction) * left + fraction * right
    }
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
    let Some(scale) = weight_scale(data, weights) else {
        return f64::NAN;
    };
    let mut sample: Vec<_> = data
        .iter()
        .copied()
        .zip(weights.iter().copied())
        .filter(|&(_, weight)| weight > 0.0)
        .collect();
    sample.sort_unstable_by(|(a, _), (b, _)| a.total_cmp(b));
    if p == 0.0 {
        return sample.first().map_or(f64::NAN, |&(value, _)| value);
    }
    let last = sample.last().map_or(f64::NAN, |&(value, _)| value);
    if p >= 1.0 {
        return last;
    }
    let total_weight: f64 = sample.iter().map(|&(_, weight)| weight / scale).sum();
    let target_weight = p * total_weight;
    let mut cumulative_weight = 0.0;
    for (value, weight) in sample {
        cumulative_weight += weight / scale;
        if cumulative_weight >= target_weight {
            return value;
        }
    }
    last
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
    let n = sample.len();
    if n == 0 {
        return f64::NAN;
    }
    let (lower, &mut middle, _) = sample.select_nth_unstable_by(n / 2, f64::total_cmp);
    if n.is_multiple_of(2) {
        let previous = lower
            .iter()
            .copied()
            .max_by(f64::total_cmp)
            .unwrap_or(middle);
        interpolate(previous, middle, 0.5)
    } else {
        middle
    }
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
    let deviations: Vec<f64> = data.iter().map(|v| (v - median_value).abs()).collect();
    median(&deviations)
}

/// Interquartile range (type 7 quantiles).
#[must_use]
pub fn iqr(data: &[f64]) -> f64 {
    let x = sorted(data);
    quantile_sorted(&x, 0.75, 7) - quantile_sorted(&x, 0.25, 7)
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
    crate::kernels::stats::mean(&x[k..x.len() - k])
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
    crate::kernels::stats::mean(&x)
}

/// Huber M-estimator of location with tuning constant `k` (1.345 is the
/// classical choice), scale fixed to the normalised MAD.
#[must_use]
pub fn huber_location(
    data: &[f64],
    k: f64,
) -> f64 {
    if !k.is_finite() || k <= 0.0 {
        return f64::NAN;
    }
    let mut mu = median(data);
    let s = NORMAL_MAD_SCALE * mad_about(data, mu);
    if s == 0.0 || !s.is_finite() {
        return mu;
    }
    let count = data.iter().filter(|v| !v.is_nan()).count() as f64;
    for _ in 0..200 {
        let step: f64 = data
            .iter()
            .filter(|v| !v.is_nan())
            .map(|&v| ((v - mu) / s).clamp(-k, k))
            .sum::<f64>()
            / count;
        mu += s * step;
        if step.abs() < 1e-13 {
            break;
        }
    }
    mu
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

/// Standard deviation with Bessel's correction.
fn std_dev(data: &[f64]) -> f64 {
    let n = data.len() as f64;
    let m = data.iter().sum::<f64>() / n;
    (data.iter().map(|v| (v - m).powi(2)).sum::<f64>() / (n - 1.0)).sqrt()
}

/// Silverman's rule-of-thumb bandwidth `0.9 min(sd, IQR/1.34) n^(-1/5)`.
#[must_use]
pub fn bandwidth_silverman(data: &[f64]) -> f64 {
    let sd = std_dev(data);
    let a = sd.min(iqr(data) / 1.34);
    let a = if a > 0.0 { a } else { sd };
    0.9 * a * (data.len() as f64).powf(-0.2)
}

/// Scott's rule bandwidth `1.06 sd n^(-1/5)`.
#[must_use]
pub fn bandwidth_scott(data: &[f64]) -> f64 {
    1.06 * std_dev(data) * (data.len() as f64).powf(-0.2)
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
            let difference = x - d;
            let u = if difference.is_infinite() && x.is_finite() && d.is_finite() {
                x / h - d / h
            } else {
                difference / h
            };
            match kernel {
                | Kernel::Gaussian => (-0.5 * u * u).exp() / (2.0 * PI).sqrt(),
                | Kernel::Epanechnikov => {
                    if u.abs() < 1.0 {
                        0.75 * (1.0 - u * u)
                    } else {
                        0.0
                    }
                },
            }
        })
        .sum();
    (s / count as f64) / h
}

/// Percentile bootstrap confidence interval for `statistic` at level
/// `1 - alpha`, with `resamples` seeded resamples.
pub fn bootstrap_ci<S: Fn(&[f64]) -> f64>(
    data: &[f64],
    statistic: S,
    resamples: usize,
    alpha: f64,
    seed: u64,
) -> (f64, f64) {
    if data.is_empty() || resamples == 0 || !(0.0..=1.0).contains(&alpha) {
        return (f64::NAN, f64::NAN);
    }
    let mut rng = Rng::new(seed);
    let n = data.len();
    let mut sample = vec![0.0; n];
    let mut stats = Vec::new();
    if stats.try_reserve_exact(resamples).is_err() {
        return (f64::NAN, f64::NAN);
    }
    for _ in 0..resamples {
        for value in &mut sample {
            *value = data[rng.below(n as u64) as usize];
        }
        let value = statistic(&sample);
        if !value.is_nan() {
            stats.push(value);
        }
    }
    stats.sort_unstable_by(f64::total_cmp);
    (
        quantile_sorted(&stats, alpha / 2.0, 7),
        quantile_sorted(&stats, 1.0 - alpha / 2.0, 7),
    )
}
