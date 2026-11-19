//! Shared failures for numerical routines.
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
// binary exponent separately, including correctly rounded subnormal results.
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
        return Err(NumericError::Domain(
            "reference-scaled result overflows Float".into(),
        ));
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
