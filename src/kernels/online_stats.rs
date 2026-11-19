//! Streaming and mergeable multivariate statistics.
//!
//! [`OnlineCovariance`] maintains a mean vector and the centred sum of products
//! using the multivariate Welford update. Independent accumulators can be
//! combined by Chan's parallel formula, without retaining their observations.
//! Storage is quadratic in the number of variables, independent of sample size.

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
        assert_eq!(s.merge(&s.clone()), Err(OnlineStatsError::Overflow));
        assert_eq!(s.count(), 1 << 53);
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
    // Delta, staged mean and staged triangle. Reused across updates.
    scratch: Vec<f64>,
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
        let scratch_len = packed
            .checked_add(dimension.saturating_mul(2))
            .ok_or(OnlineStatsError::InvalidDimension)?;
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
            scratch: zeroed(scratch_len)?,
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
        if observation.iter().any(|x| !x.is_finite()) {
            return Err(OnlineStatsError::NonFiniteInput);
        }
        if self.count == 0 {
            self.mean.copy_from_slice(observation);
            self.count = 1;
            return Ok(());
        }
        self.update(observation, 1, None)
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
        let dimension = self.dimension();
        if !observations.len().is_multiple_of(dimension) {
            return Err(OnlineStatsError::DimensionMismatch);
        }
        let rows = observations.len() / dimension;
        if rows == 0 {
            return Ok(());
        }
        if rows == 1 {
            return self.push(observations);
        }
        let count = u64::try_from(rows).map_err(|_| OnlineStatsError::Overflow)?;
        self.checked_count(count)?;

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
                for row in observations[start..end].chunks_exact(dimension) {
                    partial.push(row)?;
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
        let ratio = incoming_count as f64 / next as f64;
        let factor = self.count as f64 * ratio;
        let n = self.dimension();
        let (delta, stage) = self.scratch.split_at_mut(n);
        let (new_mean, new_m2) = stage.split_at_mut(n);
        for i in 0..n {
            delta[i] = incoming_mean[i] - self.mean[i];
            new_mean[i] = self.mean[i] + delta[i] * ratio;
        }
        let mut start = 0;
        for (i, &d) in delta.iter().enumerate() {
            let stop = start + i + 1;
            let area = d * factor;
            let old_row = &self.m2[start..stop];
            let new_row = &mut new_m2[start..stop];
            if let Some(other) = incoming_m2 {
                for (j, slot) in new_row.iter_mut().enumerate() {
                    *slot = old_row[j] + other[start + j] + area * delta[j];
                }
            } else {
                for (j, slot) in new_row.iter_mut().enumerate() {
                    *slot = old_row[j] + area * delta[j];
                }
            }
            start = stop;
        }
        if new_mean.iter().chain(new_m2.iter()).any(|x| !x.is_finite()) {
            return Err(OnlineStatsError::Overflow);
        }
        self.mean.copy_from_slice(new_mean);
        self.m2.copy_from_slice(new_m2);
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
                ((value / standard_deviations[row]) / standard_deviations[col]).clamp(-1.0, 1.0)
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
