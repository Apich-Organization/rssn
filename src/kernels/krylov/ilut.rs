//! Threshold incomplete LU for validated CSR matrices.
use super::Csr;
use super::accurate_dot;
use std::collections::BTreeMap;

/// Controls for bounded-fill, unpivoted incomplete LU.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IlutOptions {
    /// Drop working entries whose magnitude is at most this fraction of the
    /// maximum absolute value in the original row, before elimination and in U.
    /// Must be finite and in `[0, 1]`; zero disables threshold dropping.
    pub drop_tolerance: f64,
    /// Largest number of off-diagonal entries retained in each of L and U per row.
    /// The largest magnitudes are kept, with lower column indices breaking ties.
    pub max_entries_per_side: usize,
    /// Reject a pivot whose magnitude divided by its original row maximum is
    /// at most this fraction. Zero accepts every nonzero finite pivot.
    pub pivot_tolerance: f64,
}
impl Default for IlutOptions {
    fn default() -> Self {
        Self {
            drop_tolerance: 1e-4,
            max_entries_per_side: 20,
            pivot_tolerance: 1e-12,
        }
    }
}
/// An invalid matrix, unsafe pivot or non-finite intermediate.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IlutError {
    /// Options are invalid, or CSR is not a nonempty square, sorted, unique matrix.
    InvalidInput,
    /// Input or a right-hand side contains a non-finite value.
    NonFiniteInput,
    /// A diagonal is absent, zero, or below the requested relative pivot threshold.
    UnsafePivot {
        /// Row containing the unsafe pivot.
        row: usize,
    },
    /// A finite-input intermediate or triangular solution was not representable.
    NumericalBreakdown {
        /// Row where arithmetic failed.
        row: usize,
    },
    /// The input and output slices do not match the factor dimension.
    Shape,
}
/// Per-factorisation storage and dropping statistics.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct IlutStats {
    /// Number of stored strict-lower entries.
    pub lower_entries: usize,
    /// Number of stored diagonal and upper entries.
    pub upper_entries: usize,
    /// Number of working or candidate entries dropped by threshold or fill cap.
    pub dropped_entries: usize,
}
/// Bounded-fill incomplete LU factors. No pivoting or global convergence
/// guarantee is implied; use a Krylov method and verify its original residual.
#[derive(Clone, Debug)]
pub struct Ilut {
    rows: Vec<Row>,
    stats: IlutStats,
}
#[derive(Clone, Debug)]
struct Row {
    lower: Vec<(usize, f64)>,
    diagonal: f64,
    upper: Vec<(usize, f64)>,
}

pub(super) fn valid_csr(a: &Csr) -> bool {
    let n = a.rows;
    n != 0
        && a.cols == n
        && n.checked_add(1) == Some(a.indptr.len())
        && a.indices.len() == a.values.len()
        && a.indptr.first() == Some(&0)
        && a.indptr.last() == Some(&a.values.len())
        && a.indptr
            .windows(2)
            .all(|w| w[0] <= w[1] && w[1] <= a.values.len())
        && a.indptr.windows(2).all(|w| {
            let row = &a.indices[w[0]..w[1]];
            row.iter().all(|&col| col < n) && row.windows(2).all(|w| w[0] < w[1])
        })
}
fn cap(
    entries: &mut Vec<(usize, f64)>,
    max: usize,
    dropped: &mut usize,
) {
    if entries.len() > max {
        entries.sort_by(|a, b| b.1.abs().total_cmp(&a.1.abs()).then_with(|| a.0.cmp(&b.0)));
        *dropped = dropped.saturating_add(entries.len() - max);
        entries.truncate(max);
    }
    entries.sort_by_key(|entry| entry.0);
}
impl Csr {
    /// Compute threshold incomplete LU with deterministic bounded row storage.
    ///
    /// # Errors
    /// Returns [`IlutError`] for malformed input, non-finite data, unsafe
    /// pivots or unrepresentable factorisation arithmetic. No pivoting occurs.
    pub fn ilut(
        &self,
        options: IlutOptions,
    ) -> Result<Ilut, IlutError> {
        if !valid_csr(self)
            || !options.drop_tolerance.is_finite()
            || !(0.0..=1.0).contains(&options.drop_tolerance)
            || !options.pivot_tolerance.is_finite()
            || !(0.0..=1.0).contains(&options.pivot_tolerance)
        {
            return Err(IlutError::InvalidInput);
        }
        if self.values.iter().any(|v| !v.is_finite()) {
            return Err(IlutError::NonFiniteInput);
        }
        let mut rows: Vec<Row> = Vec::with_capacity(self.rows);
        let mut stats = IlutStats {
            lower_entries: 0,
            upper_entries: 0,
            dropped_entries: 0,
        };
        for i in 0..self.rows {
            let range = self.indptr[i]..self.indptr[i + 1];
            let row_scale = self.values[range.clone()]
                .iter()
                .map(|v| v.abs())
                .fold(0.0_f64, f64::max);
            if row_scale == 0.0 {
                return Err(IlutError::UnsafePivot { row: i });
            }
            let mut work = BTreeMap::new();
            for k in range {
                work.insert(self.indices[k], self.values[k]);
            }
            if !work.contains_key(&i) {
                return Err(IlutError::UnsafePivot { row: i });
            }
            let mut lower = Vec::new();
            while let Some((&col, &value)) = work.first_key_value() {
                if col >= i {
                    break;
                }
                work.remove(&col);
                if value == 0.0
                    || (options.drop_tolerance > 0.0
                        && value.abs() / row_scale <= options.drop_tolerance)
                {
                    stats.dropped_entries = stats.dropped_entries.saturating_add(1);
                    continue;
                }
                let factor = value / rows[col].diagonal;
                if !factor.is_finite() || (factor == 0.0 && value != 0.0) {
                    return Err(IlutError::NumericalBreakdown { row: i });
                }
                lower.push((col, factor));
                for &(column, upper) in &rows[col].upper {
                    let old = work.get(&column).copied().unwrap_or(0.0);
                    let updated = (-factor).mul_add(upper, old);
                    if !updated.is_finite() {
                        return Err(IlutError::NumericalBreakdown { row: i });
                    }
                    if updated == 0.0
                        && accurate_dot::product_sum([(1.0, old), (-factor, upper)]).0 != 0.0
                    {
                        return Err(IlutError::NumericalBreakdown { row: i });
                    }
                    work.insert(column, updated);
                }
            }
            let diagonal = work.remove(&i).ok_or(IlutError::UnsafePivot { row: i })?;
            if diagonal == 0.0
                || !diagonal.is_finite()
                || (options.pivot_tolerance > 0.0
                    && diagonal.abs() / row_scale <= options.pivot_tolerance)
            {
                return Err(IlutError::UnsafePivot { row: i });
            }
            let mut upper = Vec::new();
            for (col, value) in work {
                if value == 0.0
                    || (options.drop_tolerance > 0.0
                        && value.abs() / row_scale <= options.drop_tolerance)
                {
                    stats.dropped_entries = stats.dropped_entries.saturating_add(1);
                } else {
                    upper.push((col, value));
                }
            }
            cap(
                &mut lower,
                options.max_entries_per_side,
                &mut stats.dropped_entries,
            );
            cap(
                &mut upper,
                options.max_entries_per_side,
                &mut stats.dropped_entries,
            );
            stats.lower_entries = stats.lower_entries.saturating_add(lower.len());
            stats.upper_entries = stats.upper_entries.saturating_add(upper.len() + 1);
            rows.push(Row {
                lower,
                diagonal,
                upper,
            });
        }
        Ok(Ilut { rows, stats })
    }
}
// Keep the exact accumulated numerator's exponent until after division by
// the pivot; reject a nonzero intermediate that cannot be represented.
fn scaled_divide(
    fraction: f64,
    exponent: i32,
    divisor: f64,
) -> Option<f64> {
    if fraction == 0.0 {
        return Some(0.0);
    }
    let (normal, adjust) = if divisor.abs() < f64::MIN_POSITIVE {
        (divisor * 18_014_398_509_481_984.0, -54)
    } else {
        (divisor, 0)
    };
    let bits = normal.to_bits();
    let divisor_exponent = ((bits >> 52) & 0x7ff) as i32 - 1023 + adjust;
    let divisor_fraction =
        f64::from_bits((bits & ((1u64 << 63) | ((1u64 << 52) - 1))) | (1023u64 << 52));
    let mut value = fraction / divisor_fraction;
    let mut shift = exponent - divisor_exponent;
    while shift > 512 {
        value *= 2.0_f64.powi(512);
        shift -= 512;
        if !value.is_finite() {
            return None;
        }
    }
    while shift < -512 {
        value *= 2.0_f64.powi(-512);
        shift += 512;
        if value == 0.0 {
            return None;
        }
    }
    value *= 2.0_f64.powi(shift);
    (value.is_finite() && value != 0.0).then_some(value)
}
impl Ilut {
    /// Factorisation work and retained storage counts.
    #[must_use]
    pub const fn stats(&self) -> IlutStats {
        self.stats
    }

    /// Apply `L^-1` and `U^-1` without changing `out` on failure.
    ///
    /// # Errors
    /// Rejects mismatched/non-finite inputs and non-finite triangular solves.
    pub fn apply(
        &self,
        rhs: &[f64],
        out: &mut [f64],
    ) -> Result<(), IlutError> {
        let n = self.rows.len();
        if rhs.len() != n || out.len() != n {
            return Err(IlutError::Shape);
        }
        if rhs.iter().any(|v| !v.is_finite()) {
            return Err(IlutError::NonFiniteInput);
        }
        let mut work = vec![0.0; n];
        for i in 0..n {
            let terms = std::iter::once((1.0, rhs[i])).chain(
                self.rows[i]
                    .lower
                    .iter()
                    .map(|&(col, value)| (-value, work[col])),
            );
            let (fraction, exponent) = accurate_dot::product_sum(terms);
            work[i] = scaled_divide(fraction, exponent, 1.0)
                .ok_or(IlutError::NumericalBreakdown { row: i })?;
        }
        for i in (0..n).rev() {
            let terms = std::iter::once((1.0, work[i])).chain(
                self.rows[i]
                    .upper
                    .iter()
                    .map(|&(col, value)| (-value, work[col])),
            );
            let (fraction, exponent) = accurate_dot::product_sum(terms);
            work[i] = scaled_divide(fraction, exponent, self.rows[i].diagonal)
                .ok_or(IlutError::NumericalBreakdown { row: i })?;
        }
        out.copy_from_slice(&work);
        Ok(())
    }
}
