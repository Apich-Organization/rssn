//! Streaming and mergeable multivariate statistics.
//!
//! [`OnlineCovariance`] maintains a mean vector and the centred sum of products
//! using the multivariate Welford update. Independent accumulators can be
//! combined by Chan's parallel formula, without retaining their observations.
//! Storage is quadratic in the number of variables, independent of sample size.

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
            let x = [1e12 + i as f64, -1e12 + 2.0 * i as f64];
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
#[derive(Debug, Clone)]
pub struct OnlineCovariance {
    count: u64,
    mean: Vec<f64>,
    // Symmetric centred sum of products, row-major.
    m2: Vec<f64>,
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
            .filter(|_| dimension != 0)
            .ok_or(OnlineStatsError::InvalidDimension)?;
        let mut m2 = Vec::new();
        m2.try_reserve_exact(len)
            .map_err(|_| OnlineStatsError::InvalidDimension)?;
        m2.resize(len, 0.0);
        let mut mean = Vec::new();
        mean.try_reserve_exact(dimension)
            .map_err(|_| OnlineStatsError::InvalidDimension)?;
        mean.resize(dimension, 0.0);
        Ok(Self { count: 0, mean, m2 })
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
        let next = self
            .count
            .checked_add(1)
            .filter(|&n| n <= Self::MAX_COUNT)
            .ok_or(OnlineStatsError::Overflow)?;
        let delta: Vec<_> = observation
            .iter()
            .zip(&self.mean)
            .map(|(x, m)| x - m)
            .collect();
        let ratio = 1.0 / next as f64;
        let mean: Vec<_> = self
            .mean
            .iter()
            .zip(&delta)
            .map(|(m, d)| m + d * ratio)
            .collect();
        let factor = self.count as f64 * ratio;
        let m2 = self.updated_m2(&delta, factor, None)?;
        if mean.iter().any(|x| !x.is_finite()) {
            return Err(OnlineStatsError::Overflow);
        }
        self.count = next;
        self.mean = mean;
        self.m2 = m2;
        Ok(())
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
            self.clone_from(other);
            return Ok(());
        }
        let next = self
            .count
            .checked_add(other.count)
            .filter(|&n| n <= Self::MAX_COUNT)
            .ok_or(OnlineStatsError::Overflow)?;
        let delta: Vec<_> = other
            .mean
            .iter()
            .zip(&self.mean)
            .map(|(b, a)| b - a)
            .collect();
        let ratio = other.count as f64 / next as f64;
        let mean: Vec<_> = self
            .mean
            .iter()
            .zip(&delta)
            .map(|(m, d)| m + d * ratio)
            .collect();
        let factor = self.count as f64 * ratio;
        let m2 = self.updated_m2(&delta, factor, Some(&other.m2))?;
        if mean.iter().any(|x| !x.is_finite()) {
            return Err(OnlineStatsError::Overflow);
        }
        self.count = next;
        self.mean = mean;
        self.m2 = m2;
        Ok(())
    }

    fn updated_m2(
        &self,
        delta: &[f64],
        factor: f64,
        other: Option<&[f64]>,
    ) -> Result<Vec<f64>, OnlineStatsError> {
        let n = self.dimension();
        let mut result = self.m2.clone();
        for i in 0..n {
            for j in 0..=i {
                let k = i * n + j;
                let value =
                    self.m2[k] + other.map_or(0.0, |m| m[k]) + (delta[i] * factor) * delta[j];
                if !value.is_finite() {
                    return Err(OnlineStatsError::Overflow);
                }
                result[k] = value;
                result[j * n + i] = value;
            }
        }
        Ok(result)
    }

    /// Population covariance (divisor `n`), or `None` when empty.
    #[must_use]
    pub fn population_covariance(&self) -> Option<Vec<f64>> {
        (self.count != 0).then(|| self.m2.iter().map(|v| v / self.count as f64).collect())
    }

    /// Unbiased sample covariance (divisor `n - 1`), or `None` for fewer
    /// than two observations.
    #[must_use]
    pub fn sample_covariance(&self) -> Option<Vec<f64>> {
        (self.count >= 2).then(|| {
            self.m2
                .iter()
                .map(|v| v / (self.count - 1) as f64)
                .collect()
        })
    }
}
