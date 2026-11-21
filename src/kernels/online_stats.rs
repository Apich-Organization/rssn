//! Mergeable online multivariate mean and covariance.
//!
//! Streaming and mergeable multivariate statistics.
//!
//! [`OnlineCovariance`] maintains a mean vector and the centred sum of products
//! using the multivariate Welford update. Independent accumulators can be
//! combined by Chan's parallel formula, without retaining their observations.
//! Storage is quadratic in the number of variables, independent of sample size.
//! Integer frequencies let a row represent several identical observations.

// Matrix storage is checked at construction. Row and packed offsets are
// derived from that checked dimension and the corresponding slice lengths.
#![allow(clippy::indexing_slicing, clippy::arithmetic_side_effects)]

/// Failure while constructing or updating an accumulator.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OnlineStatsError {
    /// The dimension is zero, too large, or could not be allocated.
    InvalidDimension,
    /// An observation or accumulator has a different dimension.
    DimensionMismatch,
    /// An observation contains NaN or infinity.
    NonFiniteInput,
    /// The count or an intermediate arithmetic result is not representable.
    Overflow,
}

impl std::fmt::Display for OnlineStatsError {
    fn fmt(
        &self,
        formatter: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        formatter.write_str(match self {
            | Self::InvalidDimension => "online covariance dimension is zero or unallocatable",
            | Self::DimensionMismatch => "online covariance dimensions disagree",
            | Self::NonFiniteInput => "online covariance input must be finite",
            | Self::Overflow => "online covariance count or arithmetic overflow",
        })
    }
}
impl std::error::Error for OnlineStatsError {}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(
        a: f64,
        b: f64,
        tol: f64,
    ) {
        assert!((a - b).abs() <= tol, "{a} != {b}");
    }

    #[test]
    fn multivariate_population_and_sample_covariance() {
        let mut s = OnlineCovariance::new(2).unwrap();
        assert!(s.mean().is_none());
        assert!(s.sample_covariance().is_none());
        for x in [[1.0, 2.0], [2.0, 4.0], [3.0, 6.0]] {
            s.push(&x).unwrap();
        }
        assert_eq!(s.count(), 3);
        assert_eq!(s.mean().unwrap(), &[2.0, 4.0]);
        assert_eq!(s.sample_covariance().unwrap(), vec![1.0, 2.0, 2.0, 4.0]);
        for (v, expected) in s.population_covariance().unwrap().iter().zip([
            2.0 / 3.0,
            4.0 / 3.0,
            4.0 / 3.0,
            8.0 / 3.0,
        ]) {
            close(*v, expected, 1e-15);
        }
    }

    #[test]
    fn merge_matches_sequential_and_handles_empty_shards() {
        let mut all = OnlineCovariance::new(2).unwrap();
        let mut left = all.clone();
        let mut right = all.clone();
        for i in 0..100 {
            let x = [1e12 + f64::from(i), -1e12 + 2.0 * f64::from(i)];
            all.push(&x).unwrap();
            if i < 43 {
                left.push(&x).unwrap();
            } else {
                right.push(&x).unwrap();
            }
        }
        left.merge(&right).unwrap();
        let empty = OnlineCovariance::new(2).unwrap();
        left.merge(&empty).unwrap();
        for (a, b) in left.mean().unwrap().iter().zip(all.mean().unwrap()) {
            close(*a, *b, 1e-4);
        }
        for (a, b) in left
            .sample_covariance()
            .unwrap()
            .iter()
            .zip(all.sample_covariance().unwrap())
        {
            close(*a, b, 1e-9);
        }
        let mut first = empty;
        first.merge(&all).unwrap();
        assert_eq!(first.count(), 100);
        assert_eq!(first.sample_covariance(), all.sample_covariance());
    }

    #[test]
    fn rejects_invalid_input_without_mutating_state() {
        assert_eq!(
            OnlineCovariance::new(0).unwrap_err(),
            OnlineStatsError::InvalidDimension
        );
        assert_eq!(
            OnlineCovariance::new(usize::MAX).unwrap_err(),
            OnlineStatsError::InvalidDimension
        );
        let mut s = OnlineCovariance::new(1).unwrap();
        s.push(&[1e308]).unwrap();
        for (x, err) in [
            (&[][..], OnlineStatsError::DimensionMismatch),
            (&[f64::NAN][..], OnlineStatsError::NonFiniteInput),
            (&[-1e308][..], OnlineStatsError::Overflow),
        ] {
            assert_eq!(s.push(x), Err(err));
            assert_eq!(s.count(), 1);
            assert_eq!(s.mean().unwrap(), &[1e308]);
            assert_eq!(s.population_covariance().unwrap(), vec![0.0]);
        }
        let mut bad = OnlineCovariance::new(1).unwrap();
        bad.push(&[-1e308]).unwrap();
        assert_eq!(s.merge(&bad), Err(OnlineStatsError::Overflow));
        assert_eq!(s.mean().unwrap(), &[1e308]);
        assert_eq!(s.count(), 1);
    }

    #[test]
    fn count_limit_is_checked_before_mutation() {
        let mut s = OnlineCovariance::new(1).unwrap();
        s.push(&[3.0]).unwrap();
        for _ in 0..53 {
            s.merge(&s.clone()).unwrap();
        }
        assert_eq!(s.count(), 1 << 53);
        assert_eq!(s.push(&[3.0]), Err(OnlineStatsError::Overflow));
        assert_eq!(
            s.push_batch(&[f64::NAN]),
            Err(OnlineStatsError::NonFiniteInput)
        );
        assert_eq!(s.merge(&s.clone()), Err(OnlineStatsError::Overflow));
        assert_eq!(s.count(), 1 << 53);
    }

    #[test]
    fn merge_anchors_the_mean_to_the_larger_shard() {
        let mut small = OnlineCovariance::new(1).unwrap();
        small.push(&[1e150]).unwrap();
        let mut large = OnlineCovariance::new(1).unwrap();
        large.push(&[0.0]).unwrap();
        for _ in 0..52 {
            large.merge(&large.clone()).unwrap();
        }
        let mut reversed = large.clone();
        reversed.merge(&small).unwrap();
        small.merge(&large).unwrap();
        let want = 1e150 / ((1_u64 << 52) + 1) as f64;
        let mean = small.mean().unwrap()[0];
        assert!((mean / want - 1.0).abs() < 4.0 * f64::EPSILON);
        assert_eq!(small.mean(), reversed.mean());
        assert_eq!(small.sample_covariance(), reversed.sample_covariance());
    }

    #[test]
    fn tiny_representable_correlation_survives_either_column_order() {
        for swapped in [false, true] {
            let mut dominant = OnlineCovariance::new(2).unwrap();
            let data = if swapped {
                [-1e150, 0.0, 1e150, 0.0]
            } else {
                [0.0, -1e150, 0.0, 1e150]
            };
            dominant.push_batch(&data).unwrap();
            let mut minor = OnlineCovariance::new(2).unwrap();
            minor
                .push_batch(&[-1e-150, -1e-150, 1e-150, 1e-150])
                .unwrap();
            dominant.merge(&minor).unwrap();
            let actual = dominant.correlation().unwrap();
            assert_eq!(actual[1], actual[2]);
            assert!((actual[1] / 1e-300 - 1.0).abs() < 4.0 * f64::EPSILON);
        }
    }

    #[test]
    fn a_subnormal_row_delta_keeps_representable_cross_covariance() {
        let tiny = f64::from_bits(1);
        let large = 1e150;
        for row in [[large, tiny], [tiny, large]] {
            let mut moments = OnlineCovariance::new(2).unwrap();
            moments.push(&[0.0, 0.0]).unwrap();
            moments.push(&row).unwrap();
            let covariance = moments.population_covariance().unwrap();
            assert_eq!(covariance[1], tiny * (0.25 * large));
            assert_eq!(covariance[1], covariance[2]);
        }
    }

    #[test]
    fn frequencies_match_explicit_repetitions_with_sample_and_population_divisors() {
        let mut compressed = OnlineCovariance::new(2).unwrap();
        compressed
            .push_frequency_batch(&[0.0, 0.0, 1e300, -1e300, 2.0, 6.0], &[2, 0, 2])
            .unwrap();
        let mut expanded = OnlineCovariance::new(2).unwrap();
        expanded
            .push_batch(&[0.0, 0.0, 0.0, 0.0, 2.0, 6.0, 2.0, 6.0])
            .unwrap();
        assert_eq!(compressed.count(), 4);
        assert_eq!(compressed.mean().unwrap(), [1.0, 3.0]);
        assert_eq!(
            compressed.population_covariance().unwrap(),
            [1.0, 3.0, 3.0, 9.0]
        );
        assert_eq!(compressed.sample_covariance(), expanded.sample_covariance());
        assert_eq!(
            compressed.sample_covariance().unwrap(),
            [4.0 / 3.0, 4.0, 4.0, 12.0]
        );
    }

    #[test]
    fn weighted_batches_check_shape_count_and_values_atomically() {
        let mut moments = OnlineCovariance::new(2).unwrap();
        moments.push_frequency_batch(&[], &[]).unwrap();
        moments.push_frequency_batch(&[4.0, 2.0], &[0]).unwrap();
        assert_eq!(
            moments.push_frequency_batch(&[4.0, 2.0], &[(1 << 53) + 1]),
            Err(OnlineStatsError::Overflow)
        );
        assert!(moments.mean().is_none());
        moments.push_frequency_batch(&[4.0, 2.0], &[2]).unwrap();
        for (data, frequencies, expected) in [
            (&[1.0][..], &[2][..], OnlineStatsError::DimensionMismatch),
            (&[1.0, 2.0], &[], OnlineStatsError::DimensionMismatch),
            (&[1.0, f64::NAN], &[0], OnlineStatsError::NonFiniteInput),
            (&[f64::MAX, 2.0], &[1], OnlineStatsError::Overflow),
            (&[0.0, 0.0], &[1 << 53], OnlineStatsError::Overflow),
            (&[0.0; 4], &[u64::MAX, 1], OnlineStatsError::Overflow),
        ] {
            assert_eq!(
                moments.push_frequency_batch(data, frequencies),
                Err(expected)
            );
            assert_eq!(moments.count(), 2);
            assert_eq!(moments.mean().unwrap(), [4.0, 2.0]);
            assert_eq!(moments.sample_covariance().unwrap(), [0.0; 4]);
        }
        moments
            .push_frequency_batch(&[4.0, 2.0], &[(1 << 53) - 2])
            .unwrap();
        assert_eq!(moments.count(), 1 << 53);
        assert_eq!(moments.push(&[4.0, 2.0]), Err(OnlineStatsError::Overflow));
    }

    #[test]
    fn weighted_batches_preserve_fixed_shards_and_rollback_across_workers() {
        let mut data: Vec<_> = (0_u64..8192)
            .flat_map(|i| [(i % 9) as f64, (i % 7) as f64, (i % 11) as f64])
            .collect();
        let counts: Vec<_> = (0_u64..8192).map(|i| i % 7).collect();
        let run = |workers, weights: Option<&[u64]>| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(workers)
                .build()
                .unwrap()
                .install(|| {
                    let mut moments = OnlineCovariance::new(3).unwrap();
                    if let Some(weights) = weights {
                        moments.push_frequency_batch(&data, weights).unwrap();
                    } else {
                        moments.push_batch(&data).unwrap();
                    }
                    moments
                })
        };
        let one = run(1, Some(&counts));
        let mut four = run(4, Some(&counts));
        let bits = |moments: &OnlineCovariance| {
            moments
                .mean()
                .unwrap()
                .iter()
                .copied()
                .chain(moments.population_covariance().unwrap())
                .map(f64::to_bits)
                .collect::<Vec<_>>()
        };
        assert_eq!(four.count(), counts.iter().sum::<u64>());
        assert_eq!(bits(&one), bits(&four));
        let unweighted = run(1, None);
        let ones = vec![1; counts.len()];
        let unit_weights = run(4, Some(&ones));
        assert_eq!(bits(&unweighted), bits(&unit_weights));

        *data.last_mut().unwrap() = f64::NAN;
        assert_eq!(
            four.push_frequency_batch(&data, &counts),
            Err(OnlineStatsError::NonFiniteInput)
        );
        assert_eq!(one.count(), four.count());
        assert_eq!(bits(&one), bits(&four));
    }

    #[test]
    fn weighted_shards_match_integer_reference_moments() {
        for seed in 0_i64..31 {
            let mut values = Vec::new();
            let mut frequencies = Vec::new();
            let mut count = 0;
            let mut sum = [0_i64; 3];
            let mut products = [[0_i64; 3]; 3];
            for row in 0_i64..113 {
                let sample = [
                    (row * 13 + seed * 11) % 41 - 20,
                    (row * 21 + seed * 17) % 59 - 29,
                    (row * 7 + seed * 23) % 29 - 14,
                ];
                let frequency = (row * 3 + seed * 5) % 11;
                count += frequency;
                values.extend(sample.map(|x| x as f64));
                frequencies.push(frequency as u64);
                for i in 0..3 {
                    sum[i] += frequency * sample[i];
                    for j in 0..3 {
                        products[i][j] += frequency * sample[i] * sample[j];
                    }
                }
            }
            let mut left = OnlineCovariance::new(3).unwrap();
            let mut right = OnlineCovariance::new(3).unwrap();
            let split = (7 * seed) as usize % frequencies.len();
            left.push_frequency_batch(&values[..split * 3], &frequencies[..split])
                .unwrap();
            right
                .push_frequency_batch(&values[split * 3..], &frequencies[split..])
                .unwrap();
            left.merge(&right).unwrap();
            assert_eq!(left.count(), count as u64);
            let mean = left.mean().unwrap();
            let sample = left.sample_covariance().unwrap();
            for (i, &total) in sum.iter().enumerate() {
                close(mean[i], total as f64 / count as f64, 1e-12);
                for (j, &other) in sum.iter().enumerate() {
                    let numerator = count * products[i][j] - total * other;
                    let denominator = count * (count - 1);
                    close(
                        sample[i * 3 + j],
                        numerator as f64 / denominator as f64,
                        1e-10,
                    );
                }
            }
        }
    }
}

/// Online mean and covariance of equally weighted vector observations.
///
/// Covariance matrices are returned in row-major order. Both [`Self::push`]
/// and [`Self::merge`] leave the accumulator unchanged on error. To avoid
/// silently rounding observation counts when converting them to `f64`, the
/// number of observations is limited to 2^53.
#[derive(Clone)]
pub struct OnlineCovariance {
    count: u64,
    mean: Vec<f64>,
    // Centred sum of products in row-major packed lower-triangle storage.
    m2: Vec<f64>,
    // Delta, staged mean and optional rescaled delta, reused across updates.
    scratch: Vec<f64>,
    // Staged triangle, swapped with m2 on commit.
    staged_m2: Vec<f64>,
}

impl std::fmt::Debug for OnlineCovariance {
    fn fmt(
        &self,
        formatter: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        formatter
            .debug_struct("OnlineCovariance")
            .field("count", &self.count)
            .field("mean", &self.mean)
            .field("m2", &self.m2)
            .finish_non_exhaustive()
    }
}

impl OnlineCovariance {
    const MAX_COUNT: u64 = 1_u64 << 53;

    /// Construct an empty accumulator for `dimension` variables.
    ///
    /// # Errors
    /// Returns [`OnlineStatsError::InvalidDimension`] for zero or an
    /// unallocatable matrix size.
    pub fn new(dimension: usize) -> Result<Self, OnlineStatsError> {
        let len = dimension
            .checked_mul(dimension)
            .filter(|&len| len != 0 && len <= isize::MAX as usize / std::mem::size_of::<f64>())
            .ok_or(OnlineStatsError::InvalidDimension)?;
        let packed = len / 2 + dimension.div_ceil(2);
        let zeroed = |len| -> Result<Vec<f64>, OnlineStatsError> {
            let mut values = Vec::new();
            values
                .try_reserve_exact(len)
                .map_err(|_| OnlineStatsError::InvalidDimension)?;
            values.resize(len, 0.0);
            Ok(values)
        };
        Ok(Self {
            count: 0,
            mean: zeroed(dimension)?,
            m2: zeroed(packed)?,
            scratch: zeroed(3 * dimension)?,
            staged_m2: zeroed(packed)?,
        })
    }

    /// Number of variables.
    #[must_use]
    pub const fn dimension(&self) -> usize {
        self.mean.len()
    }

    /// Number of observations incorporated so far.
    #[must_use]
    pub const fn count(&self) -> u64 {
        self.count
    }

    /// The mean, or `None` if no observation has been incorporated.
    #[must_use]
    pub fn mean(&self) -> Option<&[f64]> {
        (self.count != 0).then_some(&self.mean)
    }

    /// Incorporate one observation.
    ///
    /// # Errors
    /// Reports a dimension mismatch, non-finite input, or arithmetic overflow.
    pub fn push(
        &mut self,
        observation: &[f64],
    ) -> Result<(), OnlineStatsError> {
        if observation.len() != self.dimension() {
            return Err(OnlineStatsError::DimensionMismatch);
        }
        self.push_counted(observation, 1)
    }

    fn push_counted(
        &mut self,
        observation: &[f64],
        count: u64,
    ) -> Result<(), OnlineStatsError> {
        if observation.iter().any(|x| !x.is_finite()) {
            return Err(OnlineStatsError::NonFiniteInput);
        }
        if count == 0 {
            return Ok(());
        }
        if self.count == 0 {
            self.count = self.checked_count(count)?;
            self.mean.copy_from_slice(observation);
            return Ok(());
        }
        self.update(observation, count, None)
    }

    /// Incorporate a flat row-major batch atomically.
    ///
    /// Each row contains [`Self::dimension`] values. Large batches use bounded
    /// shards and combine them in row order, independent of worker count.
    /// An empty batch leaves the accumulator unchanged.
    ///
    /// # Errors
    /// Returns an error for an incomplete row, non-finite input, count or
    /// arithmetic overflow, or failure to allocate a shard.
    pub fn push_batch(
        &mut self,
        observations: &[f64],
    ) -> Result<(), OnlineStatsError> {
        self.push_rows(observations, None)
    }

    /// Incorporate a flat row-major batch with an integer frequency per row.
    ///
    /// Each row represents as many identical observations as its frequency,
    /// with no need to expand those observations. A zero frequency adds no
    /// observations. All row values must be finite. The resulting count is the
    /// sum of frequencies and the previous count. Sample covariance divides
    /// by the resulting count minus one. Floating-point results can differ
    /// from pushing each repetition individually. Shards and their merge order
    /// are independent of worker count.
    ///
    /// # Errors
    /// An incomplete row, a frequency vector of the wrong length, a non-finite
    /// value, count or arithmetic overflow, or a shard allocation failure
    /// leaves the accumulator unchanged.
    pub fn push_frequency_batch(
        &mut self,
        observations: &[f64],
        frequencies: &[u64],
    ) -> Result<(), OnlineStatsError> {
        self.push_rows(observations, Some(frequencies))
    }

    fn push_rows(
        &mut self,
        observations: &[f64],
        frequencies: Option<&[u64]>,
    ) -> Result<(), OnlineStatsError> {
        let dimension = self.dimension();
        if !observations.len().is_multiple_of(dimension) {
            return Err(OnlineStatsError::DimensionMismatch);
        }
        let rows = observations.len() / dimension;
        let count = if let Some(frequencies) = frequencies {
            if frequencies.len() != rows {
                return Err(OnlineStatsError::DimensionMismatch);
            }
            frequencies
                .iter()
                .copied()
                .try_fold(0_u64, u64::checked_add)
                .ok_or(OnlineStatsError::Overflow)?
        } else {
            u64::try_from(rows).map_err(|_| OnlineStatsError::Overflow)?
        };
        if rows == 1 {
            return self.push_counted(observations, count);
        }
        self.checked_count(count)?;
        if count == 0 {
            return if observations.iter().all(|value| value.is_finite()) {
                Ok(())
            } else {
                Err(OnlineStatsError::NonFiniteInput)
            };
        }

        // Shards depend on input shape only. Cap scratch retention and work
        // inflation when each partial covariance already occupies many bytes.
        let min_rows = (32_768 / self.m2.len()).clamp(64, 1_024);
        let max_shards = (262_144 / self.m2.len()).clamp(1, 8);
        let rows_per_shard = rows.div_ceil((rows / min_rows).clamp(1, max_shards));
        let shard_size = rows_per_shard * dimension;
        let shard_count = rows.div_ceil(rows_per_shard);
        let shards =
            super::par::map_indexed(shard_count, 2, |shard| -> Result<Self, OnlineStatsError> {
                let start = shard * shard_size;
                let end = start.saturating_add(shard_size).min(observations.len());
                let mut partial = Self::new(dimension)?;
                let rows = observations[start..end].chunks_exact(dimension);
                if let Some(frequencies) = frequencies {
                    let start_row = start / dimension;
                    for (row, &frequency) in rows.zip(&frequencies[start_row..]) {
                        partial.push_counted(row, frequency)?;
                    }
                } else {
                    for row in rows {
                        partial.push_counted(row, 1)?;
                    }
                }
                Ok(partial)
            });
        let mut shards = shards.into_iter();
        let mut batch = shards.next().ok_or(OnlineStatsError::DimensionMismatch)??;
        for shard in shards {
            batch.merge(&shard?)?;
        }
        self.merge(&batch)
    }

    /// Merge another independent accumulator into this one.
    ///
    /// # Errors
    /// Reports a dimension mismatch, count overflow, or arithmetic overflow.
    pub fn merge(
        &mut self,
        other: &Self,
    ) -> Result<(), OnlineStatsError> {
        if self.dimension() != other.dimension() {
            return Err(OnlineStatsError::DimensionMismatch);
        }
        if other.count == 0 {
            return Ok(());
        }
        if self.count == 0 {
            self.count = other.count;
            self.mean.copy_from_slice(&other.mean);
            self.m2.copy_from_slice(&other.m2);
            return Ok(());
        }
        self.update(&other.mean, other.count, Some(&other.m2))
    }

    fn checked_count(
        &self,
        added: u64,
    ) -> Result<u64, OnlineStatsError> {
        self.count
            .checked_add(added)
            .filter(|&n| n <= Self::MAX_COUNT)
            .ok_or(OnlineStatsError::Overflow)
    }

    fn update(
        &mut self,
        incoming_mean: &[f64],
        incoming_count: u64,
        incoming_m2: Option<&[f64]>,
    ) -> Result<(), OnlineStatsError> {
        let next = self.checked_count(incoming_count)?;
        // Start at the mean of the larger shard to avoid subtracting nearly
        // equal magnitudes when the smaller shard has a distant mean.
        let (base, target, larger, smaller) = if self.count >= incoming_count {
            (&self.mean[..], incoming_mean, self.count, incoming_count)
        } else {
            (incoming_mean, &self.mean[..], incoming_count, self.count)
        };
        let ratio = smaller as f64 / next as f64;
        let factor = larger as f64 * ratio;
        let n = self.dimension();
        let (delta, remaining) = self.scratch.split_at_mut(n);
        let (new_mean, scaled_delta) = remaining.split_at_mut(n);
        for i in 0..n {
            delta[i] = target[i] - base[i];
            new_mean[i] = base[i] + delta[i] * ratio;
        }
        if new_mean.iter().any(|x| !x.is_finite()) {
            return Err(OnlineStatsError::Overflow);
        }
        let mut start = 0;
        let mut finite = true;
        let mut scaled = false;
        for (i, &d) in delta.iter().enumerate() {
            let stop = start + i + 1;
            let area = d * factor;
            let (coefficient, source) = if area.abs() < f64::MIN_POSITIVE && d != 0.0 {
                // Scaling an already tiny row delta can underflow before
                // multiplication by a large column delta recovers the result.
                if !scaled {
                    for (result, value) in scaled_delta.iter_mut().zip(delta.iter()) {
                        *result = value * factor;
                    }
                    scaled = true;
                }
                (d, &scaled_delta[..])
            } else {
                (area, &delta[..])
            };
            let old_row = &self.m2[start..stop];
            let new_row = &mut self.staged_m2[start..stop];
            if let Some(other) = incoming_m2 {
                for (j, slot) in new_row.iter_mut().enumerate() {
                    *slot = old_row[j] + other[start + j] + coefficient * source[j];
                    finite &= slot.is_finite();
                }
            } else {
                for (j, slot) in new_row.iter_mut().enumerate() {
                    *slot = old_row[j] + coefficient * source[j];
                    finite &= slot.is_finite();
                }
            }
            start = stop;
        }
        if !finite {
            return Err(OnlineStatsError::Overflow);
        }
        self.mean.copy_from_slice(new_mean);
        std::mem::swap(&mut self.m2, &mut self.staged_m2);
        self.count = next;
        Ok(())
    }

    /// Population covariance (divisor `n`), or `None` when empty.
    #[must_use]
    pub fn population_covariance(&self) -> Option<Vec<f64>> {
        (self.count != 0).then(|| self.matrix(|_, _, v| v / self.count as f64))
    }

    /// Unbiased sample covariance (divisor `n - 1`), or `None` for fewer
    /// than two observations.
    #[must_use]
    pub fn sample_covariance(&self) -> Option<Vec<f64>> {
        (self.count >= 2).then(|| self.matrix(|_, _, v| v / (self.count - 1) as f64))
    }

    /// Row-major Pearson correlation matrix, or `None` when the sample has
    /// fewer than two observations or any column has zero variance.
    #[must_use]
    pub fn correlation(&self) -> Option<Vec<f64>> {
        if self.count < 2 {
            return None;
        }
        let mut diagonal = 0;
        let mut standard_deviations = Vec::with_capacity(self.dimension());
        for i in 0..self.dimension() {
            let deviation = self.m2[diagonal].sqrt();
            if deviation <= 0.0 {
                return None;
            }
            standard_deviations.push(deviation);
            diagonal += i + 2;
        }
        Some(self.matrix(|row, col, value| {
            if row == col {
                1.0
            } else {
                let a = standard_deviations[row];
                let b = standard_deviations[col];
                // Dividing by the smaller deviation first avoids losing a
                // representable correlation to an underflowed intermediate.
                ((value / a.min(b)) / a.max(b)).clamp(-1.0, 1.0)
            }
        }))
    }

    fn matrix(
        &self,
        transform: impl Fn(usize, usize, f64) -> f64,
    ) -> Vec<f64> {
        let n = self.dimension();
        let mut matrix = vec![0.0; n * n];
        let mut index = 0;
        for i in 0..n {
            for j in 0..=i {
                let value = transform(i, j, self.m2[index]);
                matrix[i * n + j] = value;
                matrix[j * n + i] = value;
                index += 1;
            }
        }
        matrix
    }
}
