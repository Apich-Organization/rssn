//! # Robust statistics, quantiles and kernel density estimation
//!
//! Quantile definitions (Hyndman-Fan types 1 to 9), median absolute
//! deviation, trimmed/winsorised means, Huber and Tukey-biweight location
//! estimators, empirical CDF, bandwidth selectors, Gaussian /
//! Epanechnikov kernel density estimation, and seeded bootstrap
//! confidence intervals.
#![allow(
    clippy::manual_midpoint,
    clippy::missing_const_for_fn,
    clippy::struct_field_names,
    clippy::or_fun_call,
    clippy::manual_swap,
    clippy::if_not_else,
    clippy::unnecessary_sort_by,
    clippy::while_float,
    clippy::too_long_first_doc_paragraph,
    clippy::cast_sign_loss,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::needless_pass_by_value,
    clippy::manual_map,
    clippy::unnecessary_map_or,
    clippy::suboptimal_flops,
    clippy::similar_names,
    clippy::unreadable_literal,
    clippy::excessive_precision,
    clippy::needless_range_loop,
    clippy::float_cmp,
    clippy::too_many_lines,
    clippy::cognitive_complexity,
    clippy::option_if_let_else,
    clippy::many_single_char_names,
    clippy::type_complexity,
    clippy::too_many_arguments,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_precision_loss,
    clippy::map_unwrap_or,
    clippy::cast_lossless
)]

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
pub fn quantile(data: &[f64], p: f64, kind: u8) -> f64 {
    let x = sorted(data);
    let n = x.len();
    if n == 0 || !(0.0..=1.0).contains(&p) || !(1..=9).contains(&kind) {
        return f64::NAN;
    }
    let nf = n as f64;
    let at = |i: isize| x[(i.clamp(1, n as isize) - 1) as usize];
    if kind <= 3 {
        {
            let np = nf * p;
            let j = np.floor();
            let g = np - j;
            let ji = j as isize;
            match kind {
                1 => if g > 0.0 { at(ji + 1) } else { at(ji) },
                2 => if g > 0.0 { at(ji + 1) } else { 0.5 * (at(ji) + at(ji + 1)) },
                _ => {
                    let np = nf * p - 0.5;
                    let j = np.floor();
                    let g = np - j;
                    if g == 0.0 && (j as isize) % 2 == 0 { at(j as isize) } else { at(j as isize + 1) }
                }
            }
        }
    } else {
        {
            let (a, b) = match kind {
                4 => (0.0, 1.0),
                5 => (0.5, 0.5),
                6 => (0.0, 0.0),
                7 => (1.0, 1.0),
                8 => (1.0 / 3.0, 1.0 / 3.0),
                _ => (0.375, 0.375),
            };
            // plotting position m = a + p (1 - a - b), index h = n p + m
            let m = a + p * (1.0 - a - b);
            let h = nf * p + m;
            let j = h.floor();
            let g = h - j;
            let ji = j as isize;
            at(ji) + g * (at(ji + 1) - at(ji))
        }
    }
}

/// Median of `data` (NaNs ignored).
#[must_use]
pub fn median(data: &[f64]) -> f64 {
    quantile(data, 0.5, 7)
}

/// Median absolute deviation; multiplied by 1.4826 when `normal` is true so
/// it estimates the standard deviation for Gaussian data.
#[must_use]
pub fn mad(data: &[f64], normal: bool) -> f64 {
    let m = median(data);
    let dev: Vec<f64> = data.iter().map(|v| (v - m).abs()).collect();
    median(&dev) * if normal { 1.482_602_218_505_602 } else { 1.0 }
}

/// Interquartile range (type 7 quantiles).
#[must_use]
pub fn iqr(data: &[f64]) -> f64 {
    quantile(data, 0.75, 7) - quantile(data, 0.25, 7)
}

/// Mean after discarding the fraction `trim` from each tail.
#[must_use]
pub fn trimmed_mean(data: &[f64], trim: f64) -> f64 {
    let x = sorted(data);
    let k = (x.len() as f64 * trim.clamp(0.0, 0.499)).floor() as usize;
    let core = &x[k..x.len() - k];
    core.iter().sum::<f64>() / core.len().max(1) as f64
}

/// Mean after replacing the lowest and highest fraction `trim` by the
/// nearest retained values.
#[must_use]
pub fn winsorized_mean(data: &[f64], trim: f64) -> f64 {
    let mut x = sorted(data);
    let n = x.len();
    let k = (n as f64 * trim.clamp(0.0, 0.499)).floor() as usize;
    if n == 0 {
        return f64::NAN;
    }
    let (lo, hi) = (x[k], x[n - 1 - k]);
    for v in &mut x {
        *v = v.clamp(lo, hi);
    }
    x.iter().sum::<f64>() / n as f64
}

/// Huber M-estimator of location with tuning constant `k` (1.345 is the
/// classical choice), scale fixed to the normalised MAD.
#[must_use]
pub fn huber_location(data: &[f64], k: f64) -> f64 {
    let mut mu = median(data);
    let s = mad(data, true);
    if s == 0.0 {
        return mu;
    }
    for _ in 0..200 {
        let step: f64 = data.iter().map(|&v| ((v - mu) / s).clamp(-k, k)).sum::<f64>() / data.len() as f64;
        mu += s * step;
        if step.abs() < 1e-13 {
            break;
        }
    }
    mu
}

/// Tukey biweight location estimator (tuning constant `c`, 4.685 typical).
#[must_use]
pub fn biweight_location(data: &[f64], c: f64) -> f64 {
    let mut mu = median(data);
    let s = mad(data, true);
    if s == 0.0 {
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
pub fn ecdf(data: &[f64], x: f64) -> f64 {
    let v = sorted(data);
    v.partition_point(|&a| a <= x) as f64 / v.len().max(1) as f64
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
    let a = std_dev(data).min(iqr(data) / 1.34);
    let a = if a > 0.0 { a } else { std_dev(data) };
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
pub fn kde(data: &[f64], x: f64, h: f64, kernel: Kernel) -> f64 {
    let s: f64 = data
        .iter()
        .map(|&d| {
            let u = (x - d) / h;
            match kernel {
                Kernel::Gaussian => (-0.5 * u * u).exp() / (2.0 * PI).sqrt(),
                Kernel::Epanechnikov => if u.abs() < 1.0 { 0.75 * (1.0 - u * u) } else { 0.0 },
            }
        })
        .sum();
    s / (data.len() as f64 * h)
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
    let mut rng = Rng::new(seed);
    let n = data.len();
    let mut stats: Vec<f64> = (0..resamples)
        .map(|_| {
            let s: Vec<f64> = (0..n).map(|_| data[rng.below(n as u64) as usize]).collect();
            statistic(&s)
        })
        .collect();
    stats.sort_by(f64::total_cmp);
    (quantile(&stats, alpha / 2.0, 7), quantile(&stats, 1.0 - alpha / 2.0, 7))
}
