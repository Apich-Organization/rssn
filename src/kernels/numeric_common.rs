//! Shared failures and domain checks for numerical routines.
use std::num::NonZeroUsize;

/// Number of snapshots including the initial state and the final step.
pub(crate) const fn trajectory_sample_count(
    steps: usize,
    stride: NonZeroUsize,
) -> Option<usize> {
    steps.div_ceil(stride.get()).checked_add(1)
}

/// Why a numeric kernel refused.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum NumericError {
    /// A matrix was expected to be square and was not.
    NotSquare,
    /// The shapes of two operands do not agree.
    ShapeMismatch(String),
    /// The matrix has no inverse / no unique solution.
    Singular,
    /// A matrix that should be positive definite is not.
    NotPositiveDefinite,
    /// An iterative method ran out of iterations.
    NotConverged,
    /// The input is outside the function's domain.
    Domain(String),
    /// The input was empty where a value was required.
    Empty,
}

impl std::fmt::Display for NumericError {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        match self {
            | Self::NotSquare => f.write_str("the matrix is not square"),
            | Self::ShapeMismatch(detail) => write!(f, "shape mismatch: {detail}"),
            | Self::Singular => f.write_str("the matrix is singular"),
            | Self::NotPositiveDefinite => f.write_str("the matrix is not positive definite"),
            | Self::NotConverged => f.write_str("the iteration did not converge"),
            | Self::Domain(detail) => write!(f, "outside the domain: {detail}"),
            | Self::Empty => f.write_str("the input is empty"),
        }
    }
}

impl std::error::Error for NumericError {}

// Scaling before subtraction recovers finite results when the unscaled
// distance between two finite observations exceeds the float range.
pub(super) fn scaled_difference(
    left: f64,
    right: f64,
    scale: f64,
) -> f64 {
    let difference = left - right;
    if difference.is_finite() {
        difference / scale
    } else {
        left / scale - right / scale
    }
}

// Positive reference scales may individually be extreme while their product
// or quotient is representable. Keep mantissas normalized and accumulate the
// binary exponent separately. The final output is rounded from that accumulator.
fn parts(mut x: f64) -> (f64, i32) {
    let shift = if x < f64::MIN_POSITIVE {
        x *= 18_014_398_509_481_984.0;
        -54
    } else {
        0
    };
    let bits = x.to_bits();
    let e = ((bits >> 52) & 2047) as i32 - 1023 + shift;
    let m = f64::from_bits((bits & 0x000f_ffff_ffff_ffff) | (1023_u64 << 52));
    (m, e)
}

pub(crate) fn scale_product(
    value: f64,
    numerator: &[f64],
    denominator: &[f64],
) -> Result<f64, NumericError> {
    const OVERFLOW: &str = "reference-scaled result overflows Float";
    if !value.is_finite()
        || numerator.len() + denominator.len() > 256
        || numerator
            .iter()
            .chain(denominator)
            .any(|v| !v.is_finite() || *v <= 0.0)
    {
        return Err(NumericError::Domain(
            "invalid reference scale product".into(),
        ));
    }
    if value == 0.0 {
        return Ok(value);
    }
    // A single hardware operation can round directly to a subnormal. The
    // normalized accumulator would round first at a larger exponent.
    let single = match (numerator, denominator) {
        | ([], []) => return Ok(value),
        | ([factor], []) => Some(value * factor),
        | ([], [factor]) => Some(value / factor),
        | _ => None,
    };
    if let Some(out) = single {
        return if out.is_finite() {
            Ok(out)
        } else {
            Err(NumericError::Domain(OVERFLOW.into()))
        };
    }
    let (mut m, mut e) = parts(value.abs());
    for factor in numerator {
        let (a, b) = parts(*factor);
        m *= a;
        e += b;
        if m >= 2.0 {
            m *= 0.5;
            e += 1;
        }
    }
    for factor in denominator {
        let (a, b) = parts(*factor);
        m /= a;
        e -= b;
        if m < 1.0 {
            m *= 2.0;
            e -= 1;
        }
    }
    if e > 1023 {
        return Err(NumericError::Domain(OVERFLOW.into()));
    }
    let sign = value.to_bits() & (1_u64 << 63);
    let fraction = m.to_bits() & 0x000f_ffff_ffff_ffff;
    let bits = if e >= -1022 {
        (((e + 1023) as u64) << 52) | fraction
    } else {
        let shift = (-1022 - e) as u32;
        if shift > 53 {
            0
        } else {
            let mantissa = (1_u64 << 52) | fraction;
            let integer = mantissa >> shift;
            let remainder = mantissa & ((1_u64 << shift) - 1);
            let half = 1_u64 << (shift - 1);
            integer + u64::from(remainder > half || (remainder == half && integer & 1 != 0))
        }
    };
    Ok(f64::from_bits(sign | bits))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_single_scale_rounds_directly_at_the_subnormal_boundary() {
        let magnitude = f64::from_bits(0x0019_6791_79ab_83a0);
        let factor = f64::from_bits(0x3fc7_2982_a0f5_f2ab);
        for sign in [0, 1_u64 << 63] {
            let value = f64::from_bits(magnitude.to_bits() | sign);
            assert_eq!(
                scale_product(value, &[factor], &[]).map(f64::to_bits),
                Ok(0x0004_98d9_3d97_d067 | sign),
            );
        }
        let dividend = f64::from_bits(0x001c_4ac2_72ec_b6e0);
        let divisor = f64::from_bits(0x4014_df2c_80e4_b612);
        assert_eq!(
            scale_product(dividend, &[], &[divisor]).map(f64::to_bits),
            Ok(0x0005_6c0c_7b2b_6d33),
        );
    }

    #[test]
    fn scale_ranges_and_signed_zero_are_preserved() {
        let huge = f64::from_bits(1773_u64 << 52);
        let tiny = f64::from_bits(273_u64 << 52);
        assert_eq!(
            scale_product(huge, &[huge, tiny], &[]).map(f64::to_bits),
            Ok(huge.to_bits()),
        );
        assert_eq!(
            scale_product(-f64::from_bits(1), &[tiny], &[]).map(f64::to_bits),
            Ok((-0.0_f64).to_bits()),
        );
        assert_eq!(
            scale_product(-0.0, &[], &[tiny]).map(f64::to_bits),
            Ok((-0.0_f64).to_bits()),
        );
        assert!(scale_product(huge, &[huge], &[]).is_err());
        assert!(scale_product(-0.0, &[], &[0.0]).is_err());
        assert!(scale_product(1.0, &[f64::INFINITY], &[]).is_err());
    }
}
