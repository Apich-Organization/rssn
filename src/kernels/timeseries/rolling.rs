//! Rolling mean and variance for finite real time series.

use std::fmt;

use crate::cancellation;
use crate::kernels::dd::Dd;
use crate::kernels::par::try_for_chunks_mut;
use crate::kernels::reuse::f64_scratch;
use crate::kernels::simd::min_max_f64;

const CHECK_INTERVAL: usize = 4096;
const TASK_SAMPLES: usize = 2048;
const PARALLEL_SAMPLES: usize = 16_384;
const SUMMARY_SIZE: usize = 4;
const SHORT_WINDOW: usize = 16;
const SHORT_MIN_MAGNITUDE: f64 = f64::from_bits((1023 - 400) << 52);
const SHORT_MAX_MAGNITUDE: f64 = f64::from_bits((1023 + 400) << 52);

/// Mean and variance of one full window.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct RollingMoment {
    /// Arithmetic mean.
    pub mean: f64,
    /// Sum of squared deviations divided by the window width minus `ddof`.
    pub variance: f64,
}

/// A rolling-moments calculation could not produce a result.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RollingError {
    /// The requested width was zero.
    InvalidWindow,
    /// The degrees-of-freedom correction was at least the window width.
    InvalidDegreesOfFreedom,
    /// A sample was NaN or infinite.
    NonFiniteInput {
        /// Zero-based sample index.
        index: usize,
    },
    /// The variance exceeded the range of `f64`.
    VarianceOverflow {
        /// Zero-based starting index of the first overflowing window.
        start: usize,
    },
    /// The computation was interrupted by its embedding host.
    Cancelled,
}

impl fmt::Display for RollingError {
    fn fmt(
        &self,
        formatter: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        match *self {
            | Self::InvalidWindow => formatter.write_str("rolling window width must be positive"),
            | Self::InvalidDegreesOfFreedom => {
                formatter.write_str("rolling variance ddof must be less than the window width")
            },
            | Self::NonFiniteInput { index } => {
                write!(formatter, "non-finite sample at index {index}")
            },
            | Self::VarianceOverflow { start } => {
                write!(
                    formatter,
                    "rolling variance overflow at window start {start}"
                )
            },
            | Self::Cancelled => cancellation::Cancelled.fmt(formatter),
        }
    }
}

impl std::error::Error for RollingError {}

#[derive(Clone, Copy, Default)]
struct Summary {
    count: usize,
    sum: Dd,
    // The unnormalised sum of squared deviations is scale² × squares.
    scale: f64,
    squares: f64,
}

impl Summary {
    const fn singleton(
        sample: f64,
        sum_scale: f64,
    ) -> Self {
        Self {
            count: 1,
            sum: Dd::from_f64(sample * sum_scale),
            scale: 0.0,
            squares: 0.0,
        }
    }

    fn add_scaled(
        &mut self,
        scale: f64,
        squares: f64,
    ) {
        if scale > self.scale {
            let ratio = self.scale / scale;
            self.squares = squares + self.squares * ratio * ratio;
            self.scale = scale;
        } else if scale > 0.0 {
            let ratio = scale / self.scale;
            self.squares += squares * ratio * ratio;
        }
    }

    fn merge(
        self,
        other: Self,
        sum_scale: f64,
    ) -> Self {
        let (base, other) = if self.count >= other.count {
            (self, other)
        } else {
            (other, self)
        };
        if other.count == 0 {
            return base;
        }
        let count = base.count + other.count;
        let mut combined = Self {
            count,
            sum: base.sum + other.sum,
            ..base
        };
        if !base.scale.is_finite() || !other.scale.is_finite() {
            combined.scale = f64::INFINITY;
            return combined;
        }

        combined.add_scaled(other.scale, other.squares);
        // Chan, Golub and LeVeque's pairwise variance update adds a
        // nonnegative correction for the difference between group means.
        // Form the common numerator before division to keep a small
        // difference even when the group means have a large common offset.
        let base_count = base.count as f64;
        let other_count = other.count as f64;
        let weights = base_count * other_count;
        let weighted_base = if other.count == 1 {
            base.sum
        } else {
            base.sum * Dd::from_f64(other_count)
        };
        let weighted_other = if base.count == 1 {
            other.sum
        } else {
            other.sum * Dd::from_f64(base_count)
        };
        let extent = ((weighted_other - weighted_base).to_f64() / (weights * sum_scale)).abs();
        if extent.is_finite() {
            combined.add_scaled(extent, weights / count as f64);
        } else {
            // The variance will overflow even if later observations dilute
            // the contribution of this difference between means.
            combined.scale = f64::INFINITY;
        }
        combined
    }

    const fn store(self) -> [f64; SUMMARY_SIZE] {
        [self.sum.hi, self.sum.lo, self.scale, self.squares]
    }

    const fn load(
        count: usize,
        [hi, lo, scale, squares]: [f64; SUMMARY_SIZE],
    ) -> Self {
        Self {
            count,
            sum: Dd { hi, lo },
            scale,
            squares,
        }
    }

    fn moment(
        self,
        divisor: f64,
        sum_scale: f64,
        start: usize,
    ) -> Result<RollingMoment, RollingError> {
        // Apply the divisor before the two scale multiplications so that an
        // overflowing sum of squares may still produce a finite variance.
        let variance = self.scale * (self.scale * (self.squares / divisor));
        if variance.is_finite() {
            // Sum the observations before dividing to preserve small residuals
            // when observations of opposite signs cancel.
            let count = self.count as f64;
            let leading = self.sum.hi / count;
            let remainder = (-leading).mul_add(count, self.sum.hi) + self.sum.lo;
            let mean = (leading + remainder / count) / sum_scale;
            Ok(RollingMoment { mean, variance })
        } else {
            Err(RollingError::VarianceOverflow { start })
        }
    }
}

fn checkpoint() -> Result<(), RollingError> {
    cancellation::checkpoint().map_err(|_| RollingError::Cancelled)
}

fn rolling_block(
    samples: &[f64],
    width: usize,
    divisor: f64,
    sum_scale: f64,
    start: usize,
    output: &mut [RollingMoment],
    suffixes: &mut [[f64; SUMMARY_SIZE]],
) -> Result<(), RollingError> {
    let left = &samples[start..start + width];
    let output_len = output.len();
    let mut suffix = Summary::default();
    // Saved suffixes are needed only at starts that have a full window.
    for values in left[output_len..].rchunks(CHECK_INTERVAL) {
        checkpoint()?;
        for &sample in values.iter().rev() {
            suffix = suffix.merge(Summary::singleton(sample, sum_scale), sum_scale);
        }
    }
    for (index, (saved, &sample)) in suffixes.iter_mut().zip(left).rev().enumerate() {
        if index.is_multiple_of(CHECK_INTERVAL) {
            checkpoint()?;
        }
        suffix = suffix.merge(Summary::singleton(sample, sum_scale), sum_scale);
        *saved = suffix.store();
    }

    let mut prefix = Summary::default();
    for (offset, (result, &saved)) in output.iter_mut().zip(suffixes.iter()).enumerate() {
        if offset.is_multiple_of(CHECK_INTERVAL) {
            checkpoint()?;
        }
        *result = Summary::load(width - offset, saved)
            .merge(prefix, sum_scale)
            .moment(divisor, sum_scale, start + offset)?;
        if offset + 1 < output_len {
            prefix = prefix.merge(
                Summary::singleton(samples[start + width + offset], sum_scale),
                sum_scale,
            );
        }
    }
    Ok(())
}

// Dispatch is limited to inputs whose differences from any anchor are exact
// and whose squared deviations are far from floating-point range limits.
fn short_windows(
    samples: &[f64],
    width: usize,
    divisor: f64,
    output: &mut [RollingMoment],
) -> Result<(), RollingError> {
    let mean_scale = 1.0 / width as f64;
    let variance_scale = 1.0 / divisor;
    try_for_chunks_mut(
        output,
        TASK_SAMPLES,
        PARALLEL_SAMPLES / TASK_SAMPLES,
        |task, result| {
            checkpoint()?;
            let starts = samples[task * TASK_SAMPLES..].windows(width);
            for (result, values) in result.iter_mut().zip(starts) {
                let anchor = values[0];
                let average_shift =
                    values.iter().map(|value| value - anchor).sum::<f64>() * mean_scale;
                let squares = values
                    .iter()
                    .map(|value| ((value - anchor) - average_shift).powi(2))
                    .sum::<f64>();
                *result = RollingMoment {
                    mean: anchor + average_shift,
                    variance: squares * variance_scale,
                };
            }
            Ok(())
        },
    )
}

/// Computes a mean and variance for every complete window in `samples`.
///
/// Result `i` describes `samples[i..i + window]`. The result is empty when
/// `window` exceeds the input length. Pass `ddof = 0` for population variance
/// or `ddof = 1` for sample variance. Samples are always required to be finite.
///
/// Windows are evaluated in linear time. Windows of at most 16 samples use
/// two passes around an anchor when the whole series has the same sign, its
/// magnitudes are within a factor of two, and their exponents leave enough
/// headroom to square the deviations. For other series and window widths,
/// fixed blocks of window starts combine suffix and prefix moments. Large
/// inputs use the configured worker pool, with identical arithmetic order
/// for any worker count.
///
/// # Errors
/// Returns [`RollingError`] for a zero window width, `ddof >= window`, a
/// non-finite sample, an overflowing variance, or host cancellation.
///
/// # Example
/// ```
/// use rssn::kernels::timeseries::rolling::RollingMoment;
/// use rssn::kernels::timeseries::rolling::rolling_moments;
///
/// assert_eq!(
///     rolling_moments(&[1.0, 3.0, 5.0], 2, 0)?,
///     vec![
///         RollingMoment {
///             mean: 2.0,
///             variance: 1.0
///         },
///         RollingMoment {
///             mean: 4.0,
///             variance: 1.0
///         },
///     ],
/// );
/// # Ok::<(), rssn::kernels::timeseries::rolling::RollingError>(())
/// ```
pub fn rolling_moments(
    samples: &[f64],
    window: usize,
    ddof: usize,
) -> Result<Vec<RollingMoment>, RollingError> {
    if window == 0 {
        return Err(RollingError::InvalidWindow);
    }
    if ddof >= window {
        return Err(RollingError::InvalidDegreesOfFreedom);
    }
    checkpoint()?;
    let mut minimum = f64::INFINITY;
    let mut maximum = f64::NEG_INFINITY;
    for (chunk, values) in samples.chunks(CHECK_INTERVAL).enumerate() {
        checkpoint()?;
        let (lo, hi) = min_max_f64(values).unwrap_or_default();
        if (!lo.is_finite() || !hi.is_finite())
            && let Some(offset) = values.iter().position(|x| !x.is_finite())
        {
            return Err(RollingError::NonFiniteInput {
                index: chunk * CHECK_INTERVAL + offset,
            });
        }
        minimum = minimum.min(lo);
        maximum = maximum.max(hi);
    }
    if window > samples.len() {
        return Ok(Vec::new());
    }
    if window == 1 {
        return Ok(samples
            .iter()
            .map(|&mean| RollingMoment { mean, variance: 0.0 })
            .collect());
    }

    let divisor = (window - ddof) as f64;
    let max_abs = minimum.abs().max(maximum.abs());
    let min_abs = if minimum > 0.0 {
        minimum
    } else {
        -maximum
    };
    let mut output = vec![RollingMoment::default(); samples.len() - window + 1];
    // Sterbenz's lemma guarantees exact subtraction for same-sign values
    // within a factor of two. The exponent margin also bounds squared differences.
    if window <= SHORT_WINDOW
        && min_abs >= max_abs * 0.5
        && (SHORT_MIN_MAGNITUDE..=SHORT_MAX_MAGNITUDE).contains(&max_abs)
    {
        short_windows(samples, window, divisor, &mut output)?;
        return Ok(output);
    }
    // A binary scale leaves headroom for weighted sums and
    // can be removed exactly from their average. Ordinary-size data needs no
    // scaling, so small observations retain their full exponent range.
    let width = window as f64;
    let sum_scale = if max_abs > (f64::MAX * 0.25) / (width * width) {
        let shift = window.next_power_of_two().trailing_zeros() * 2 + 2;
        f64::from_bits(u64::from(1023 - shift) << 52)
    } else {
        1.0
    };
    let task_length = TASK_SAMPLES.div_ceil(window) * window;
    let grain = PARALLEL_SAMPLES.div_ceil(task_length).max(2);
    try_for_chunks_mut(&mut output, task_length, grain, |task, results| {
        let scratch_len = results.len().min(window) * SUMMARY_SIZE;
        let mut scratch = f64_scratch(scratch_len);
        scratch.resize(scratch_len, 0.0);
        let (suffixes, _) = scratch.as_chunks_mut::<SUMMARY_SIZE>();
        for (block, results) in results.chunks_mut(window).enumerate() {
            let start = task * task_length + block * window;
            let suffixes = &mut suffixes[..results.len()];
            rolling_block(
                samples, window, divisor, sum_scale, start, results, suffixes,
            )?;
        }
        Ok(())
    })?;
    Ok(output)
}
