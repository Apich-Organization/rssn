//! Exact products accumulated in fixed binary bins, rounded once at the boundary.
const LOW: i32 = -2148;
const BITS: usize = 32;
// Maximum product bit is 4195 relative to LOW. The remaining bins provide
// more than 64 carry bits for a slice with at most usize::MAX terms.
const BINS: usize = 136;
const BASE: i128 = 1i128 << BITS;

struct Sum {
    bins: [i128; BINS],
    exceptional: f64,
    first: usize,
    last: usize,
}
impl Sum {
    const fn new() -> Self {
        Self {
            bins: [0; BINS],
            exceptional: 0.0,
            first: BINS,
            last: 0,
        }
    }

    const fn parts(x: f64) -> (u64, i32) {
        let bits = x.to_bits();
        let fraction = bits & ((1 << 52) - 1);
        let exponent = ((bits >> 52) & 0x7ff) as i32;
        if exponent == 0 {
            (fraction, -1074)
        } else {
            (fraction | (1 << 52), exponent - 1075)
        }
    }

    fn add(
        &mut self,
        a: f64,
        b: f64,
    ) {
        if !a.is_finite() || !b.is_finite() {
            self.exceptional += a * b;
            return;
        }
        let (a_m, a_e) = Self::parts(a);
        let (b_m, b_e) = Self::parts(b);
        if a_m == 0 || b_m == 0 {
            return;
        }
        let product = u128::from(a_m) * u128::from(b_m);
        let offset = (a_e + b_e - LOW) as usize;
        let index = offset / BITS;
        let shift = offset % BITS;
        self.first = self.first.min(index);
        self.last = self.last.max(index + 3);
        let sign = if a.is_sign_negative() ^ b.is_sign_negative() {
            -1
        } else {
            1
        };
        // Splitting before shifting avoids overflowing u128 for 106-bit products.
        for limb in 0..4 {
            let word = ((product >> (limb * BITS)) & 0xffff_ffff) as i128 * sign;
            self.bins[index + limb] += word << shift;
        }
    }

    fn normalize(&mut self) -> bool {
        if self.first == BINS {
            return false;
        }
        // Normalize only the interval touched by products, extending it when
        // a carry reaches a previously untouched higher bin.
        let mut i = self.first;
        while i <= self.last {
            let carry = self.bins[i] / BASE;
            self.bins[i] %= BASE;
            if carry != 0 {
                self.bins[i + 1] += carry;
                self.last = self.last.max(i + 1);
            }
            i += 1;
        }
        let negative = self.bins[self.first..=self.last]
            .iter()
            .rev()
            .find(|&&x| x != 0)
            .is_some_and(|&x| x < 0);
        if negative {
            for bin in &mut self.bins[self.first..=self.last] {
                *bin = -*bin;
            }
        }
        let mut i = self.first;
        while i <= self.last {
            let carry = self.bins[i].div_euclid(BASE);
            self.bins[i] = self.bins[i].rem_euclid(BASE);
            if carry != 0 {
                self.bins[i + 1] += carry;
                self.last = self.last.max(i + 1);
            }
            i += 1;
        }
        negative
    }

    const fn bit(
        &self,
        pos: usize,
    ) -> bool {
        pos / BITS < BINS && ((self.bins[pos / BITS] >> (pos % BITS)) & 1) != 0
    }

    fn any_below(
        &self,
        pos: usize,
    ) -> bool {
        let full = (pos / BITS).min(BINS);
        (self.first < full && self.bins[self.first..full].iter().any(|&x| x != 0))
            || (full < BINS
                && !pos.is_multiple_of(BITS)
                && (self.bins[full] & ((1i128 << (pos % BITS)) - 1)) != 0)
    }

    fn rounded(
        &self,
        high: usize,
        low: usize,
    ) -> u64 {
        let mut n = if low <= high {
            let index = low / BITS;
            let word = |i: usize| self.bins.get(i).copied().unwrap_or(0) as u128;
            let chunk = word(index) | (word(index + 1) << 32) | (word(index + 2) << 64);
            let width = high - low + 1;
            ((chunk >> (low % BITS)) & ((1u128 << width) - 1)) as u64
        } else {
            0u64
        };
        if low > 0 && self.bit(low - 1) && (self.any_below(low - 1) || n & 1 != 0) {
            n += 1;
        }
        n
    }

    fn high(&self) -> Option<usize> {
        if self.first == BINS {
            return None;
        }
        let i = self.first
            + self.bins[self.first..=self.last]
                .iter()
                .rposition(|&v| v != 0)?;
        Some(i * BITS + (127 - self.bins[i].leading_zeros() as usize))
    }

    fn normalized(mut self) -> (f64, i32) {
        if self.exceptional != 0.0 || self.exceptional.is_nan() {
            return (self.exceptional, 0);
        }
        let neg = self.normalize();
        let Some(high) = self.high() else {
            return (0.0, 0);
        };
        let low = high.saturating_sub(52);
        let value = self.rounded(high, low) as f64;
        let sign = if neg { -1.0 } else { 1.0 };
        (
            sign * value / 2.0_f64.powi((high - low) as i32),
            LOW + high as i32,
        )
    }

    fn value(mut self) -> f64 {
        if self.exceptional != 0.0 || self.exceptional.is_nan() {
            return self.exceptional;
        }
        let neg = self.normalize();
        let Some(high) = self.high() else {
            return 0.0;
        };
        let exp = LOW + high as i32;
        let sign = if neg { 1u64 << 63 } else { 0 };
        if exp < -1022 {
            let low = (-1074 - LOW) as usize;
            let n = self.rounded(high, low);
            return f64::from_bits(sign | n);
        }
        let low = high.saturating_sub(52);
        let mut n = self.rounded(high, low);
        let mut exp = exp;
        if n == 1u64 << 53 {
            n >>= 1;
            exp += 1;
        }
        if exp > 1023 {
            return f64::from_bits(sign | (0x7ffu64 << 52));
        }
        f64::from_bits(sign | (((exp + 1023) as u64) << 52) | (n & ((1u64 << 52) - 1)))
    }
}
fn accumulate(pairs: impl IntoIterator<Item = (f64, f64)>) -> Sum {
    let mut sum = Sum::new();
    for (a, b) in pairs {
        sum.add(a, b);
    }
    sum
}
/// Return a signed normalized mantissa and unbounded binary exponent.
pub(super) fn product_sum<I>(pairs: I) -> (f64, i32)
where
    I: IntoIterator<Item = (f64, f64)>,
    I::IntoIter: Clone,
{
    let pairs = pairs.into_iter();
    if let Some(value) = filtered_bulk(pairs.clone())
        .or_else(|| filtered(pairs.clone()))
        .filter(|value| value.is_normal() && value.abs() > f64::MIN_POSITIVE)
    {
        let (mantissa, exponent) = Sum::parts(value);
        let fraction = mantissa as f64 / (1u64 << 52) as f64;
        return (
            if value.is_sign_negative() {
                -fraction
            } else {
                fraction
            },
            exponent + 52,
        );
    }
    accumulate(pairs).normalized()
}
// Certified common-range path for long dots. Every nonzero rounded product is
// between 2^-863 and 2^928: even usize::MAX terms cannot overflow an f64,
// and every exact product bit and compensation remains representable. The
// final outward error interval is checked exactly as in the general filter.
fn filtered_bulk<I>(pairs: I) -> Option<f64>
where
    I: Iterator<Item = (f64, f64)>,
{
    let (mut sum, mut correction, mut bound) = (0.0_f64, 0.0_f64, 0.0_f64);
    for (a, b) in pairs {
        if a == 0.0 || b == 0.0 {
            if !a.is_finite() || !b.is_finite() {
                return None;
            }
            continue;
        }
        let product = a * b;
        let exponent = (product.to_bits() >> 52) & 0x7ff;
        if !(160..=1950).contains(&exponent) {
            return None;
        }
        let residual = a.mul_add(b, -product);
        let next = sum + product;
        let back = next - sum;
        let error = (sum - (next - back)) + (product - back);
        sum = next;
        for item in [residual, error] {
            if item != 0.0 {
                let next = correction + item;
                let back = next - correction;
                let dropped = (correction - (next - back)) + (item - back);
                correction = next;
                if dropped != 0.0 {
                    bound = (bound + dropped.abs()).next_up();
                }
            }
        }
    }
    let result = sum + correction;
    if bound == 0.0 {
        return Some(if result == 0.0 { 0.0 } else { result });
    }
    if !result.is_normal() {
        return None;
    }
    let back = result - sum;
    let remainder = (sum - (result - back)) + (correction - back);
    let low = result.next_down();
    let high = result.next_up();
    if !low.is_finite() || !high.is_finite() {
        return None;
    }
    let error_low = (remainder - bound).next_down();
    let error_high = (remainder + bound).next_up();
    if error_low > -(result - low) * 0.5 && error_high < (high - result) * 0.5 {
        Some(result)
    } else {
        None
    }
}

// Error-free transformation for finite additions. Inputs and intermediates
// that would require special handling are rejected by the filter.
fn two_sum(
    a: f64,
    b: f64,
) -> Option<(f64, f64)> {
    let sum = a + b;
    if !sum.is_finite() {
        return None;
    }
    let back = sum - a;
    if !back.is_finite() {
        return None;
    }
    let error = (a - (sum - back)) + (b - back);
    if !error.is_finite() {
        return None;
    }
    Some((sum, error))
}

// The exact dot lies in s+c ± error_bound. Every discarded correction is
// accumulated with outward rounding; inconclusive cases use exact bins.
fn filtered_bounded(pairs: impl Iterator<Item = (f64, f64)>) -> Option<(f64, bool)> {
    let (mut sum, mut correction, mut error_bound) = (0.0, 0.0, 0.0_f64);
    for (a, b) in pairs {
        if !a.is_finite() || !b.is_finite() {
            return None;
        }
        if a == 0.0 || b == 0.0 {
            continue;
        }
        if !a.is_normal() || !b.is_normal() {
            return None;
        }
        let product = a * b;
        if !product.is_normal() {
            return None;
        }
        // Scaling a finite normal by an exact power of two preserves its
        // significand whenever the resulting product is normal.
        let power_of_two = |value: f64| value.to_bits() & ((1u64 << 52) - 1) == 0;
        let residue = if power_of_two(a) || power_of_two(b) {
            0.0
        } else {
            let (_, ae) = Sum::parts(a);
            let (_, be) = Sum::parts(b);
            // Every exact product bit, including the FMA residual, must be
            // representable even at the subnormal boundary.
            if ae + be < -1074 {
                return None;
            }
            a.mul_add(b, -product)
        };
        let (next, addition_error) = two_sum(sum, product)?;
        sum = next;
        for error in [residue, addition_error] {
            if error != 0.0 {
                let (next, dropped) = two_sum(correction, error)?;
                if dropped != 0.0 {
                    error_bound = (error_bound + dropped.abs()).next_up();
                    if !error_bound.is_finite() {
                        return None;
                    }
                }
                correction = next;
            }
        }
    }
    let (result, remainder) = two_sum(sum, correction)?;
    if error_bound == 0.0 {
        return Some((if result == 0.0 { 0.0 } else { result }, false));
    }
    if !result.is_normal() {
        return None;
    }
    let lower = result.next_down();
    let upper = result.next_up();
    if !lower.is_finite() || !upper.is_finite() {
        return None;
    }
    let lower_half = (result - lower) * 0.5;
    let upper_half = (upper - result) * 0.5;
    let error_low = (remainder - error_bound).next_down();
    let error_high = (remainder + error_bound).next_up();
    if error_low > -lower_half && error_high < upper_half {
        Some((result, true))
    } else {
        None
    }
}
fn filtered(pairs: impl Iterator<Item = (f64, f64)>) -> Option<f64> {
    filtered_bounded(pairs).map(|(value, _)| value)
}

// For short rows with exact power-of-two products, two error-free additions
// retain the entire sum whenever the accumulated correction itself is exact.
// Difficult rows are handed to the existing certified filter and exact bins.
fn short_exact_products<I>(pairs: I) -> Option<f64>
where
    I: Iterator<Item = (f64, f64)>,
{
    if pairs.size_hint().1.is_none_or(|length| length > 8) {
        return None;
    }
    let (mut sum, mut correction) = (0.0_f64, 0.0_f64);
    let (mut low, mut high) = (u64::MAX, 0u64);
    for (a, b) in pairs {
        if a == 0.0 || b == 0.0 {
            if !a.is_finite() || !b.is_finite() {
                return None;
            }
            continue;
        }
        let bits = a.to_bits();
        let exponent = (bits >> 52) & 0x7ff;
        if exponent == 0 || exponent == 0x7ff || bits & ((1u64 << 52) - 1) != 0 {
            return None;
        }
        let product = a * b;
        let product_exponent = (product.to_bits() >> 52) & 0x7ff;
        // Eight such products and all two-sum intermediates stay finite;
        // even the smallest correction cannot underflow.
        if !(65..=2040).contains(&product_exponent) {
            return None;
        }
        low = low.min(product_exponent);
        high = high.max(product_exponent);
        if high - low > 40 {
            return None;
        }
        let next = sum + product;
        let back = next - sum;
        let error = (sum - (next - back)) + (product - back);
        sum = next;
        // At most eight input floats span 40 exponent positions. Their
        // two-sum errors are multiples of the smallest input ulp, and their
        // total fits within 46 significant bits, so this addition is exact.
        correction += error;
    }
    let result = sum + correction;
    if result.is_finite() {
        Some(if result == 0.0 { 0.0 } else { result })
    } else {
        None
    }
}

/// Round the exact dot product directly to f64, including subnormal ties.
pub(super) fn dot_value<I>(pairs: I) -> f64
where
    I: IntoIterator<Item = (f64, f64)>,
    I::IntoIter: Clone,
{
    let pairs = pairs.into_iter();
    // A single finite multiplication is already correctly rounded by IEEE
    // binary64, including overflow and gradual underflow. Exact zero operands
    // are canonicalized as in the binary accumulator; non-finite cases retain
    // the general path's exceptional-value handling.
    if pairs.size_hint() == (1, Some(1))
        && let Some((a, b)) = pairs.clone().next()
        && a.is_finite()
        && b.is_finite()
    {
        return if a == 0.0 || b == 0.0 {
            0.0
        } else {
            a * b
        };
    }
    short_exact_products(pairs.clone())
        .or_else(|| filtered_bulk(pairs.clone()))
        .or_else(|| filtered(pairs.clone()))
        .unwrap_or_else(|| accumulate(pairs).value())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn mixed_scale_and_cancellation_are_exact() {
        assert_eq!(dot_value([(1e300, 1e-300), (1e-300, 1e300)]), 2.0);
        assert!((dot_value([(1e308, 2.0), (9e307, -2.0)]) / 2e307 - 1.0).abs() < 1e-15);
        assert_eq!(dot_value([(1e308, 2.0), (1e308, -2.0), (1.0, 1.0)]), 1.0);
    }
    #[test]
    fn filtered_rounding_matches_exact_bins_for_varied_products() {
        let mut seed = 0xa076_1d64_78bd_642fu64;
        let mut accepted = 0;
        for case in 0..40_000 {
            let mut pairs = Vec::new();
            for _ in 0..(case % 9 + 1) {
                seed ^= seed << 7;
                seed ^= seed >> 9;
                let a = f64::from_bits(
                    (seed & ((1u64 << 52) - 1))
                        | (((seed >> 52) % 2047) << 52)
                        | (seed & (1u64 << 63)),
                );
                seed ^= seed << 7;
                seed ^= seed >> 9;
                let b = f64::from_bits(
                    (seed & ((1u64 << 52) - 1))
                        | (((seed >> 52) % 2047) << 52)
                        | (seed & (1u64 << 63)),
                );
                pairs.push((a, b));
            }
            let exact = accumulate(pairs.iter().copied()).value();
            if let Some(value) = filtered(pairs.iter().copied()) {
                accepted += 1;
                assert_eq!(value.to_bits(), exact.to_bits(), "case {case}");
            }
            assert_eq!(dot_value(pairs).to_bits(), exact.to_bits(), "case {case}");
        }
        assert!(accepted > 100, "filter should handle ordinary finite cases");
        // Exercise the power-of-two shortcut over signs and the entire
        // normal exponent range, including overflow/subnormal fallbacks.
        for case in 0..20_000 {
            seed ^= seed << 7;
            seed ^= seed >> 9;
            let power = f64::from_bits((1 + seed % 2046) << 52 | (seed & (1u64 << 63)));
            seed ^= seed << 7;
            seed ^= seed >> 9;
            let operand = f64::from_bits(
                (seed & ((1u64 << 52) - 1))
                    | ((1 + ((seed >> 52) % 2046)) << 52)
                    | (seed & (1u64 << 63)),
            );
            let pairs = [(power, operand), (-power, operand), (1.0, operand)];
            assert_eq!(
                dot_value(pairs).to_bits(),
                accumulate(pairs).value().to_bits(),
                "power case {case}"
            );
        }

        // A half-ulp tie, exact cancellation, and the smallest subnormal.
        for pairs in [
            vec![(1.0, 1.0), (f64::EPSILON, 0.5)],
            vec![(1.0, 1.0), (-1.0, 1.0)],
            vec![(f64::from_bits(1), 1.0)],
        ] {
            assert_eq!(
                dot_value(pairs.clone()).to_bits(),
                accumulate(pairs).value().to_bits()
            );
        }
    }
    #[test]
    fn long_vectors_and_midpoints_match_exact_accumulation() {
        let mut seed = 0x9e37_79b9_7f4a_7c15u64;
        let mut bounded_accepts = 0;
        for case in 0..2_000 {
            let mut pairs = Vec::new();
            for _ in 0..(16 + case % 113) {
                seed ^= seed << 7;
                seed ^= seed >> 9;
                let a = f64::from_bits(
                    (seed & ((1u64 << 52) - 1)) | ((993 + seed % 61) << 52) | (seed & (1u64 << 63)),
                );
                seed ^= seed << 7;
                seed ^= seed >> 9;
                let b = f64::from_bits(
                    (seed & ((1u64 << 52) - 1)) | ((993 + seed % 61) << 52) | (seed & (1u64 << 63)),
                );
                pairs.push((a, b));
            }
            let exact = accumulate(pairs.iter().copied()).value();
            if let Some((value, bounded)) = filtered_bounded(pairs.iter().copied()) {
                bounded_accepts += usize::from(bounded);
                assert_eq!(value.to_bits(), exact.to_bits(), "bounded case {case}");
            }
            assert_eq!(
                dot_value(pairs.iter().copied()).to_bits(),
                exact.to_bits(),
                "case {case}"
            );
            let (actual, ae) = product_sum(pairs.iter().copied());
            let (mut expected, mut ee) = accumulate(pairs.iter().copied()).normalized();
            if expected.abs() == 2.0 {
                expected *= 0.5;
                ee += 1;
            }
            assert_eq!(
                (actual.to_bits(), ae),
                (expected.to_bits(), ee),
                "scaled case {case}"
            );
        }
        assert!(
            bounded_accepts > 0,
            "test must exercise nonzero error-bound acceptance"
        );
        for delta in [-f64::from_bits(1), 0.0, f64::from_bits(1)] {
            let pairs = [(1.0, 1.0), (f64::EPSILON, 0.5), (delta, 1.0)];
            assert_eq!(
                dot_value(pairs).to_bits(),
                accumulate(pairs).value().to_bits()
            );
        }
    }
    #[test]
    fn power_of_two_stencil_products_are_certified() {
        let pairs = [
            (4.0, 1.234),
            (-1.0, 0.42),
            (-1.0, -0.13),
            (-1.0, 0.83),
            (-1.0, 0.45),
        ];
        let fast =
            filtered(pairs.into_iter()).expect("power-of-two stencil should take the filter");
        assert_eq!(fast.to_bits(), accumulate(pairs).value().to_bits());
        let pairs = [(f64::MIN_POSITIVE, 0.5)];
        assert!(
            filtered(pairs.into_iter()).is_none(),
            "subnormal products need exact bins"
        );
        assert_eq!(
            dot_value(pairs).to_bits(),
            accumulate(pairs).value().to_bits()
        );
    }
    #[test]
    fn subnormal_and_overflow_rounding() {
        let tiny = f64::from_bits(1);
        assert_eq!(dot_value([(tiny, 0.5)]), 0.0);
        assert_eq!(dot_value([(tiny, 0.5), (tiny, 1.0)]), tiny * 2.0);
        assert!(dot_value([(f64::MAX, 2.0)]).is_infinite());
        assert!(dot_value([(f64::INFINITY, 1.0), (f64::NEG_INFINITY, 1.0)]).is_nan());
        let (m, e) = product_sum([(1e-300, 1e-300)]);
        assert!(m != 0.0 && e < -1074);
    }
    #[test]
    fn fast_paths_match_exact_bins_bitwise_across_scales_and_cancellation() {
        let mut seed = 0x3c79_ac49_2ba7_b653_u64;
        let mut short_accepted = 0;
        let mut bulk_accepted = 0;
        for case in 0..50_000 {
            let count = if case % 2 == 0 {
                1 + case % 8
            } else {
                9 + case % 64
            };
            let mut pairs = Vec::with_capacity(count);
            for index in 0..count {
                seed ^= seed << 7;
                seed ^= seed >> 9;
                seed ^= seed << 8;
                let sign = seed & (1 << 63);
                let exponent = 900 + (seed % 240);
                let b = f64::from_bits(sign | (exponent << 52) | (seed & ((1 << 52) - 1)));
                seed ^= seed << 7;
                seed ^= seed >> 9;
                seed ^= seed << 8;
                let a = if count <= 8 {
                    f64::from_bits((1020 + seed % 7) << 52 | (seed & (1 << 63)))
                } else {
                    f64::from_bits(
                        (1000 + seed % 47) << 52 | (seed & ((1 << 52) - 1)) | (seed & (1 << 63)),
                    )
                };
                pairs.push((a, b));
                if index > 0 && case % 13 == 0 && index % 3 == 0 {
                    let (x, y) = pairs[index - 1];
                    pairs[index] = (-x, y);
                }
            }
            let exact = accumulate(pairs.iter().copied()).value();
            if let Some(fast) = short_exact_products(pairs.iter().copied()) {
                short_accepted += 1;
                assert_eq!(fast.to_bits(), exact.to_bits(), "short case {case}");
            }
            if let Some(fast) = filtered_bulk(pairs.iter().copied()) {
                bulk_accepted += 1;
                assert_eq!(fast.to_bits(), exact.to_bits(), "bulk case {case}");
            }
            assert_eq!(
                dot_value(pairs.iter().copied()).to_bits(),
                exact.to_bits(),
                "result case {case}"
            );
        }
        assert!(short_accepted > 5_000 && bulk_accepted > 5_000);
        for pairs in [vec![(0.0, f64::INFINITY)], vec![(f64::NAN, 0.0)]] {
            assert!(short_exact_products(pairs.iter().copied()).is_none());
            assert!(dot_value(pairs.iter().copied()).is_nan());
        }
        let pairs = vec![(0.0, f64::INFINITY); 9];
        assert!(filtered_bulk(pairs.iter().copied()).is_none());
        // Exercise the accepted exponent envelope, including subnormal
        // operands with a normal product and severe mixed-sign cancellation.
        for case in 0..20_000u64 {
            let base = match case % 4 {
                | 0 => 65,
                | 1 => 160,
                | 2 => 1900,
                | _ => 2000,
            };
            let mut pairs = Vec::new();
            for index in 0..8u64 {
                seed ^= seed << 7;
                seed ^= seed >> 9;
                seed ^= seed << 8;
                let exp = base + (seed % 38);
                let b = f64::from_bits((exp << 52) | (seed & ((1 << 52) - 1)) | (seed & (1 << 63)));
                let a = if index % 3 == 0 { -1.0 } else { 1.0 };
                pairs.push((a, b));
                if index > 0 && case % 5 == 0 && index % 2 == 1 {
                    let (x, y) = pairs[index as usize - 1];
                    pairs[index as usize] = (-x, y);
                }
            }
            let exact = accumulate(pairs.iter().copied()).value();
            if let Some(fast) = short_exact_products(pairs.iter().copied()) {
                assert_eq!(fast.to_bits(), exact.to_bits(), "edge {case}");
            }
        }
        for case in 0..10_000u64 {
            let exp = if case % 2 == 0 { 160 } else { 1949 };
            let mut pairs = Vec::new();
            for index in 0..16u64 {
                seed ^= seed << 7;
                seed ^= seed >> 9;
                seed ^= seed << 8;
                let a = f64::from_bits((1023u64 << 52) | (seed & ((1 << 52) - 1)));
                seed ^= seed << 7;
                seed ^= seed >> 9;
                seed ^= seed << 8;
                let b = f64::from_bits((exp << 52) | (seed & ((1 << 52) - 1)) | (seed & (1 << 63)));
                pairs.push((a, b));
                if index > 0 && case % 7 == 0 && index % 2 == 1 {
                    let (x, y) = pairs[index as usize - 1];
                    pairs[index as usize] = (-x, y);
                }
            }
            let exact = accumulate(pairs.iter().copied()).value();
            if let Some(fast) = filtered_bulk(pairs.iter().copied()) {
                assert_eq!(fast.to_bits(), exact.to_bits(), "bulk edge {case}");
            }
        }
    }
    #[test]
    fn short_general_products_match_exact_bins_with_cancellation_and_range_edges() {
        let mut seed = 0x61c8_8646_80b5_83ebu64;
        let mut accepted = 0;
        for case in 0..80_000usize {
            let mut pairs = Vec::new();
            for index in 0..(1 + case % 8) {
                seed ^= seed << 7;
                seed ^= seed >> 9;
                seed ^= seed << 8;
                let a_exp = match case % 5 {
                    | 0 => 1 + seed % 200,
                    | 1 => 1800 + seed % 247,
                    | _ => 970 + seed % 107,
                };
                let a = f64::from_bits(
                    (a_exp << 52) | (seed & ((1u64 << 52) - 1)) | (seed & (1u64 << 63)),
                );
                seed ^= seed << 7;
                seed ^= seed >> 9;
                seed ^= seed << 8;
                let b_exp = match case % 5 {
                    | 0 => 1900 - seed % 200,
                    | 1 => 100 + seed % 147,
                    | _ => 970 + seed % 107,
                };
                let b = f64::from_bits(
                    (b_exp << 52) | (seed & ((1u64 << 52) - 1)) | (seed & (1u64 << 63)),
                );
                pairs.push((a, b));
                if case % 7 == 0 && index % 2 == 1 {
                    let (x, y) = pairs[index - 1];
                    pairs[index] = (-x, y);
                }
            }
            let exact = accumulate(pairs.iter().copied()).value();
            if let Some(fast) = filtered_bulk(pairs.iter().copied()) {
                accepted += 1;
                assert_eq!(
                    fast.to_bits(),
                    exact.to_bits(),
                    "short-general filter case {case}"
                );
            }
            assert_eq!(
                dot_value(pairs.iter().copied()).to_bits(),
                exact.to_bits(),
                "short-general dispatch case {case}"
            );
        }
        assert!(accepted > 10_000);
    }
    #[test]
    fn single_product_fast_path_matches_exact_bins_including_zero_and_underflow() {
        let mut seed = 0x0d77_1d88_4a5c_b389u64;
        for case in 0..150_000 {
            seed ^= seed << 7;
            seed ^= seed >> 9;
            seed ^= seed << 8;
            let a = f64::from_bits(seed);
            seed ^= seed << 7;
            seed ^= seed >> 9;
            seed ^= seed << 8;
            let b = f64::from_bits(seed);
            if a.is_finite() && b.is_finite() {
                let exact = accumulate([(a, b)]).value();
                assert_eq!(
                    dot_value([(a, b)]).to_bits(),
                    exact.to_bits(),
                    "single case {case}"
                );
            }
        }
        for (a, b) in [
            (0.0, -1.0),
            (-0.0, 1.0),
            (-f64::from_bits(1), 0.5),
            (f64::MAX, 2.0),
        ] {
            assert_eq!(
                dot_value([(a, b)]).to_bits(),
                accumulate([(a, b)]).value().to_bits()
            );
        }
        assert!(dot_value([(f64::INFINITY, 0.0)]).is_nan());
    }
}
