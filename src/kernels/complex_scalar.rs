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
    #[must_use]
    #[allow(clippy::should_implement_trait)] // Named scalar ABI operation rather than a trait implementation.
    pub fn div(
        self,
        other: Self,
    ) -> Self {
        let denominator = other.re * other.re + other.im * other.im;
        Self::new(
            (self.re * other.re + self.im * other.im) / denominator,
            (self.im * other.re - self.re * other.im) / denominator,
        )
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

    /// The principal square root.
    #[must_use]
    pub fn sqrt(self) -> Self {
        let modulus = self.abs();
        if modulus == 0.0 {
            return Self::default();
        }
        let re = f64::midpoint(modulus, self.re).sqrt();
        let im = f64::midpoint(modulus, -self.re).sqrt().copysign(self.im);
        Self::new(re, im)
    }
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
    let mut result = Complex::new(1.0, 0.0);
    let mut base = if exponent < 0 {
        Complex::new(1.0, 0.0).div(value)
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
