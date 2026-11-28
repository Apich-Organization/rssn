//! Complex scalar arithmetic shared by language adapters and numerical kernels.

/// A complex number, in the same layout `Complex<ρ,d>` has at run time.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Complex {
    /// The real part.
    pub re: f64,
    /// The imaginary part.
    pub im: f64,
}

impl Complex {
    /// A complex number.
    #[must_use]
    pub const fn new(
        re: f64,
        im: f64,
    ) -> Self {
        Self { re, im }
    }

    /// The modulus.
    #[must_use]
    pub fn abs(self) -> f64 {
        self.re.hypot(self.im)
    }

    /// The argument, in `(-π, π]`.
    #[must_use]
    pub fn arg(self) -> f64 {
        self.im.atan2(self.re)
    }

    /// The conjugate.
    #[must_use]
    pub fn conj(self) -> Self {
        Self::new(self.re, -self.im)
    }

    /// `self * other`.
    #[must_use]
    #[allow(clippy::should_implement_trait)] // Named scalar ABI operation rather than a trait implementation.
    pub fn mul(
        self,
        other: Self,
    ) -> Self {
        Self::new(
            self.re * other.re - self.im * other.im,
            self.re * other.im + self.im * other.re,
        )
    }

    /// `self / other`.
    ///
    /// Finite inputs with a nonzero divisor use binary scaling when their
    /// components are large or small enough to lose intermediate products.
    /// The zero divisor and nonfinite arguments follow the direct IEEE
    /// componentwise formula.
    #[must_use]
    #[inline]
    #[allow(clippy::should_implement_trait)] // Named scalar ABI operation rather than a trait implementation.
    pub fn div(
        self,
        other: Self,
    ) -> Self {
        let components = [self.re, self.im, other.re, other.im];
        let [a, b, c, d] = components.map(outside_product_range);
        if (a | b | c | d) != 0 {
            return divide_outside_product_range(self, other);
        }

        divide_direct(self, other)
    }

    /// The reciprocal `1 / self`, with the scaling and special points of
    /// [`Self::div`].
    #[must_use]
    #[inline]
    pub fn recip(self) -> Self {
        Self::new(1.0, 0.0).div(self)
    }

    /// `exp(self)`.
    #[must_use]
    pub fn exp(self) -> Self {
        let scale = self.re.exp();
        Self::new(scale * self.im.cos(), scale * self.im.sin())
    }

    /// `ln(self)`, the principal branch.
    #[must_use]
    pub fn ln(self) -> Self {
        Self::new(self.abs().ln(), self.arg())
    }

    /// The principal square root, with a nonnegative real part for finite inputs.
    ///
    /// The sign of the imaginary part selects the side of the negative real
    /// axis, including when it is signed zero. Binary scaling avoids an
    /// intermediate overflow when a finite input's modulus exceeds `f64::MAX`.
    /// Nonfinite inputs follow the componentwise IEEE midpoint formula.
    #[must_use]
    pub fn sqrt(self) -> Self {
        const DIRECT_MIN: u64 = normal_power_of_two(-500).to_bits();
        const DIRECT_MAX: u64 = normal_power_of_two(500).to_bits();
        if self.im == 0.0 && self.re.is_finite() {
            let magnitude = self.re.abs().sqrt();
            return if self.re < 0.0 {
                Self::new(0.0, magnitude.copysign(self.im))
            } else {
                Self::new(magnitude, self.im)
            };
        }
        let real_abs = self.re.abs();
        let max_bits = real_abs.to_bits().max(self.im.abs().to_bits());
        #[allow(
            clippy::imprecise_flops,
            clippy::manual_midpoint,
            reason = "the magnitude guard keeps squares and midpoint sums in range"
        )]
        let larger = if (DIRECT_MIN..DIRECT_MAX).contains(&max_bits) {
            // In this range, the sum of squares is normal and finite.
            let modulus = (self.re * self.re + self.im * self.im).sqrt();
            ((modulus + real_abs) * 0.5).sqrt()
        } else {
            if max_bits >= f64::INFINITY.to_bits() {
                let modulus = self.abs();
                let re = f64::midpoint(modulus, self.re).sqrt();
                let im = f64::midpoint(modulus, -self.re).sqrt().copysign(self.im);
                return Self::new(re, im);
            }
            scaled_sqrt_component(self, max_bits)
        };
        // Use 2*re*im = input.im to recover the smaller component without
        // subtracting nearly equal numbers in the two square-root formulas.
        let smaller = self.im / (2.0 * larger);
        if self.re < 0.0 {
            Self::new(smaller.abs(), larger.copysign(self.im))
        } else {
            Self::new(larger, smaller)
        }
    }
}

#[cold]
#[allow(
    clippy::single_call_fn,
    reason = "keeps scale-dependent branches off the ordinary path"
)]
fn scaled_sqrt_component(
    value: Complex,
    max_bits: u64,
) -> f64 {
    // Undo the square root of the binary scale after evaluating the larger part.
    let (scale, undo) = if max_bits > (f64::MAX * 0.5).to_bits() {
        (0.25, 2.0)
    } else if max_bits < f64::MIN_POSITIVE.to_bits() {
        (normal_power_of_two(54), normal_power_of_two(-27))
    } else {
        return f64::midpoint(value.abs(), value.re.abs()).sqrt();
    };
    let re = value.re * scale;
    let im = value.im * scale;
    f64::midpoint(re.hypot(im), re.abs()).sqrt() * undo
}

/// A nonzero result marks a component outside the direct-formula range.
#[allow(clippy::single_call_fn)]
fn outside_product_range(value: f64) -> u64 {
    let magnitude = value.to_bits() & 0x7fff_ffff_ffff_ffff;
    // The leading biased exponent bits 011 and 100 cover exponents -255
    // to 256. The first pair changes and the second pair stays the same.
    let transitions = (magnitude ^ (magnitude << 1)) & 0x6000_0000_0000_0000;
    let unexpected = transitions ^ 0x4000_0000_0000_0000;
    unexpected & 0_u64.wrapping_sub(u64::from(magnitude != 0))
}

fn divide_direct(
    numerator: Complex,
    divisor: Complex,
) -> Complex {
    let denominator = divisor.re * divisor.re + divisor.im * divisor.im;
    Complex::new(
        (numerator.re * divisor.re + numerator.im * divisor.im) / denominator,
        (numerator.im * divisor.re - numerator.re * divisor.im) / denominator,
    )
}

#[derive(Clone, Copy)]
struct BinaryScale {
    value: f64,
    exponent: i32,
}

impl BinaryScale {
    /// Write a finite nonzero number as a signed significand in [1, 2) in
    /// magnitude and an integer exponent. Signed zeros retain their sign.
    fn from_float(value: f64) -> Self {
        const FRACTION_AND_SIGN: u64 = 0x800f_ffff_ffff_ffff;
        const UNIT_EXPONENT: u64 = 0x3ff0_0000_0000_0000;
        if value == 0.0 {
            return Self { value, exponent: 0 };
        }
        let (normal, adjustment) = if value.is_subnormal() {
            (value * 18_014_398_509_481_984.0, -54)
        } else {
            (value, 0)
        };
        let bits = normal.to_bits();
        let exponent = ((bits >> 52) & 0x7ff) as i32 - 1023 + adjustment;
        Self {
            value: f64::from_bits((bits & FRACTION_AND_SIGN) | UNIT_EXPONENT),
            exponent,
        }
    }
}

/// `2^exponent` for an exponent in the normal f64 range.
const fn normal_power_of_two(exponent: i32) -> f64 {
    f64::from_bits(((exponent + 1023) as u64) << 52)
}

/// Adds two products of significands at a common exponent.
fn product_sum(
    a: BinaryScale,
    b: BinaryScale,
    c: BinaryScale,
    d: BinaryScale,
) -> BinaryScale {
    let left = a.value * b.value;
    let right = c.value * d.value;
    let left_exponent = a.exponent + b.exponent;
    let right_exponent = c.exponent + d.exponent;
    if left == 0.0 {
        return BinaryScale {
            value: left + right,
            exponent: right_exponent,
        };
    }
    // At a gap of 55 bits, the smaller term is below half an ulp.
    if right == 0.0 || left_exponent > right_exponent + 54 {
        return BinaryScale {
            value: left,
            exponent: left_exponent,
        };
    }
    if right_exponent > left_exponent + 54 {
        return BinaryScale {
            value: right,
            exponent: right_exponent,
        };
    }
    let exponent = left_exponent.max(right_exponent);
    let left_scale = normal_power_of_two(left_exponent - exponent);
    let right_scale = normal_power_of_two(right_exponent - exponent);
    let lhs = left * left_scale;
    let rhs = right * right_scale;
    let sum = lhs + rhs;
    let offset = sum - lhs;
    let sum_error = (lhs - (sum - offset)) + (rhs - offset);
    let left_error = a.value.mul_add(b.value, -left) * left_scale;
    let right_error = c.value.mul_add(d.value, -right) * right_scale;
    BinaryScale {
        value: sum + (sum_error + (left_error + right_error)),
        exponent,
    }
}

/// Restores a scaled quotient, rounding subnormal results at the final scale.
fn quotient(
    numerator: BinaryScale,
    denominator: BinaryScale,
) -> f64 {
    let value = numerator.value / denominator.value;
    let exponent = numerator.exponent - denominator.exponent;
    if value == 0.0 {
        return value;
    }
    if (-1022..=1023).contains(&exponent) {
        return value * normal_power_of_two(exponent);
    }
    let normalized = BinaryScale::from_float(value);
    let exponent = exponent + normalized.exponent;
    match exponent {
        | 1024.. => f64::INFINITY.copysign(value),
        | -1022..=1023 => normalized.value * normal_power_of_two(exponent),
        | -1074..=-1023 => normalized.value * f64::from_bits(1_u64 << (exponent + 1074)),
        | -1075 => (normalized.value * 0.5) * f64::from_bits(1),
        | _ => 0.0_f64.copysign(value),
    }
}

/// Computes the quotient without materializing the products at the input scale.
#[allow(clippy::single_call_fn)]
#[cold]
fn divide_outside_product_range(
    numerator: Complex,
    divisor: Complex,
) -> Complex {
    let components = [numerator.re, numerator.im, divisor.re, divisor.im];
    if !components.into_iter().all(f64::is_finite) || (divisor.re == 0.0 && divisor.im == 0.0) {
        return divide_direct(numerator, divisor);
    }
    let [a, b, c, d] = components.map(BinaryScale::from_float);
    let squared_norm = product_sum(c, c, d, d);
    let negative_a = BinaryScale { value: -a.value, ..a };
    Complex::new(
        quotient(product_sum(a, c, b, d), squared_norm),
        quotient(product_sum(b, c, negative_a, d), squared_norm),
    )
}

/// `re + im i`.
#[must_use]
pub const fn complex(
    re: f64,
    im: f64,
) -> Complex {
    Complex::new(re, im)
}

/// The modulus, `|z|`.
#[must_use]
pub fn magnitude(value: Complex) -> f64 {
    value.abs()
}

/// The argument, `arg z`, in `(-π, π]`.
#[must_use]
pub fn arg(value: Complex) -> f64 {
    value.arg()
}

/// `z * conj(z)` is real; this returns that real part, `|z|^2`.
#[must_use]
pub fn norm_squared(value: Complex) -> f64 {
    value.re * value.re + value.im * value.im
}

/// `z^n` for an integer exponent, by repeated squaring.
#[must_use]
pub fn pow(
    value: Complex,
    exponent: i64,
) -> Complex {
    if exponent == 0 {
        return Complex::new(1.0, 0.0);
    }
    if exponent == -1 {
        return value.recip();
    }
    let mut result = Complex::new(1.0, 0.0);
    let mut base = if exponent < 0 {
        value.recip()
    } else {
        value
    };
    let mut remaining = exponent.unsigned_abs();
    while remaining > 0 {
        if remaining & 1 == 1 {
            result = result.mul(base);
        }
        base = base.mul(base);
        remaining >>= 1;
    }
    result
}

/// `sin z`.
#[must_use]
pub fn sin(value: Complex) -> Complex {
    Complex::new(
        value.re.sin() * value.im.cosh(),
        value.re.cos() * value.im.sinh(),
    )
}

/// `cos z`.
#[must_use]
pub fn cos(value: Complex) -> Complex {
    Complex::new(
        value.re.cos() * value.im.cosh(),
        -value.re.sin() * value.im.sinh(),
    )
}

/// `z^w = exp(w ln z)`.
#[must_use]
pub fn powc(
    base: Complex,
    exponent: Complex,
) -> Complex {
    exponent.mul(base.ln()).exp()
}

#[cfg(test)]
mod tests {
    use super::Complex;
    use super::pow;

    #[test]
    fn square_root_scales_an_overflowing_modulus_and_a_subnormal_input() {
        let large = f64::MAX;
        let radius = large.sqrt();
        let root_half_diag = f64::midpoint(2.0_f64.sqrt(), 1.0).sqrt();
        let offset_half_diag = ((2.0_f64.sqrt() - 1.0) * 0.5).sqrt();
        let z = Complex::new(large, -large).sqrt();
        assert!(z.re.to_bits().abs_diff((radius * root_half_diag).to_bits()) <= 2);
        assert!(
            z.im.to_bits()
                .abs_diff((-radius * offset_half_diag).to_bits())
                <= 2
        );

        let tiny = f64::from_bits(1);
        // (2 + i)^2 = 3 + 4i, including after multiplication by the least subnormal.
        let z = Complex::new(3.0 * tiny, 4.0 * tiny).sqrt();
        assert_eq!(z, Complex::new(2.0 * tiny.sqrt(), tiny.sqrt()));
    }

    #[test]
    fn square_root_respects_signed_zero_and_existing_nonfinite_classes() {
        for imaginary_zero in [0.0_f64, -0.0] {
            let negative = Complex::new(-4.0, imaginary_zero).sqrt();
            assert_eq!(negative.re.to_bits(), 0);
            assert_eq!(
                negative.im.to_bits(),
                2.0_f64.copysign(imaginary_zero).to_bits()
            );
            for re in [4.0_f64, 0.0, -0.0] {
                let root = Complex::new(re, imaginary_zero).sqrt();
                assert_eq!(root.re.to_bits(), re.abs().sqrt().to_bits());
                assert_eq!(root.im.to_bits(), imaginary_zero.to_bits());
            }
        }
        let positive = Complex::new(f64::INFINITY, 0.0).sqrt();
        assert_eq!(positive.re, f64::INFINITY);
        assert!(positive.im.is_nan());
        let negative = Complex::new(f64::NEG_INFINITY, -0.0).sqrt();
        assert!(negative.re.is_nan());
        assert_eq!(negative.im, f64::NEG_INFINITY);
    }

    #[test]
    fn division_avoids_overflow_and_underflow_of_intermediates() {
        let large = Complex::new(f64::MAX, f64::MAX);
        assert_eq!(large.div(large), Complex::new(1.0, 0.0));
        for scale in [1e-200, 1e200] {
            let numerator = Complex::new(scale, -scale);
            let divisor = Complex::new(scale, scale);
            assert_eq!(numerator.div(divisor), Complex::new(0.0, -1.0));
        }

        let numerator = Complex::new(1e-200, 0.0);
        let divisor = Complex::new(1e-150, 0.0);
        assert!(
            numerator
                .div(divisor)
                .re
                .to_bits()
                .abs_diff(1e-50_f64.to_bits())
                <= 1
        );

        let least = f64::from_bits(1);
        let numerator = Complex::new(f64::MAX, least);
        let divisor = Complex::new(0.0, least);
        assert_eq!(numerator.div(divisor), Complex::new(1.0, f64::NEG_INFINITY));

        let divisor = Complex::new(2.0, 0.0);
        assert_eq!(
            Complex::new(-least, 0.0).div(divisor).re.to_bits(),
            (-0.0_f64).to_bits()
        );
        assert_eq!(Complex::new(3.0 * least, 0.0).div(divisor).re, 2.0 * least);
    }

    #[test]
    fn scaled_division_keeps_residual_after_product_cancellation() {
        let scale = 2.0_f64.powi(900);
        let epsilon = f64::EPSILON;
        let numerator = Complex::new((1.0 + epsilon) * scale, scale);
        let divisor = Complex::new(scale, (1.0 - epsilon) * scale);
        let expected = epsilon.powi(2) / (1.0 + (1.0 - epsilon).powi(2));
        assert!(
            numerator
                .div(divisor)
                .im
                .to_bits()
                .abs_diff(expected.to_bits())
                <= 1
        );
    }

    #[test]
    fn reciprocal_and_negative_first_power_preserve_subnormal_components() {
        let large = Complex::new(f64::MAX, f64::MAX);
        let expected = 0.5 / f64::MAX;
        let reciprocal = large.recip();
        assert!(reciprocal.re.to_bits().abs_diff(expected.to_bits()) <= 1);
        assert!(reciprocal.im.to_bits().abs_diff((-expected).to_bits()) <= 1);
        assert_eq!(pow(large, -1), reciprocal);
        assert_eq!(pow(Complex::new(2.0, 0.0), -3), Complex::new(0.125, 0.0));

        let tiny = Complex::new(f64::from_bits(1), -0.0);
        assert_eq!(pow(tiny, -1), Complex::new(f64::INFINITY, 0.0));
        assert_eq!(tiny.recip().im.to_bits(), 0);
        assert_eq!(Complex::new(1e200, 0.0).recip().im.to_bits(), 0);
    }

    #[test]
    fn division_keeps_ordinary_values_and_special_cases() {
        let ordinary = Complex::new(5.0, 3.0);
        assert_eq!(
            ordinary.div(Complex::new(1.0, -2.0)),
            Complex::new(-0.2, 2.6)
        );
        let division_by_zero = Complex::new(-0.0, 0.0).recip();
        assert!(division_by_zero.re.is_nan() && division_by_zero.im.is_nan());
        let infinite = Complex::new(f64::INFINITY, 1.0);
        assert_eq!(
            infinite.div(Complex::new(2.0, 3.0)),
            Complex::new(f64::INFINITY, f64::NEG_INFINITY)
        );
        assert_eq!(
            Complex::new(-2.0, 0.0).recip().im.to_bits(),
            (-0.0_f64).to_bits()
        );
    }
}
