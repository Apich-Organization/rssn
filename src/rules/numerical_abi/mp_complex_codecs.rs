// Portable dyadic rectangle codecs. Each coordinate is a checked real ball.
impl Arg for mp::MpComplex {
    fn ty() -> TypeExpr {
        <Self as Output>::ty()
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let record = native_fields(value, &["real", "imag"])?;
        let real = mp::Ball::decode(native_field(record, "real")?, index)?;
        let imag = mp::Ball::decode(native_field(record, "imag")?, index)?;
        Ok(Self::new(real, imag))
    }
}

impl Output for mp::MpComplex {
    fn ty() -> TypeExpr {
        TypeExpr::Record(vec![
            ("real".into(), <mp::Ball as Output>::ty()),
            ("imag".into(), <mp::Ball as Output>::ty()),
        ])
    }

    fn encode(self) -> Result<Value, OpError> {
        Ok(Value::Record(RecordValue {
            fields: vec![
                ("real".into(), self.re.encode()?),
                ("imag".into(), self.im.encode()?),
            ],
        }))
    }
}

binding!(
    mp_complex_new, "num_mp_complex_new", "numerical_special",
    "`num_mp_complex_new(real, imag)` constructs a rectangle from two closed dyadic intervals.",
    "crate::kernels::mp::complex::MpComplex::new", DimensionPolicy::Dimensionless,
    (real: super::mp::Ball, imag: super::mp::Ball) -> super::mp::MpComplex,
    |real: super::mp::Ball, imag: super::mp::Ball| { super::mp::MpComplex::new(real, imag) }
);

binding!(
    mp_complex_point, "num_mp_complex_point", "numerical_special",
    "`num_mp_complex_point(real, imag)` constructs an exact complex point from two dyadics.",
    "crate::kernels::mp::complex::MpComplex::point", DimensionPolicy::Dimensionless,
    (real: super::mp::Real, imag: super::mp::Real) -> super::mp::MpComplex,
    |real: super::mp::Real, imag: super::mp::Real| { super::mp::MpComplex::point(real, imag) }
);

binding!(
    mp_complex_from_real, "num_mp_complex_from_real", "numerical_special",
    "`num_mp_complex_from_real(real)` lifts the real interval into a rectangle with an exact zero imaginary part.",
    "crate::kernels::mp::complex::MpComplex::from_real", DimensionPolicy::Dimensionless,
    (real: super::mp::Ball) -> super::mp::MpComplex,
    |real: super::mp::Ball| { super::mp::MpComplex::from_real(real) }
);

binding!(
    mp_complex_zero, "num_mp_complex_zero", "numerical_special",
    "`num_mp_complex_zero(precision)` constructs the exact zero rectangle at the requested bit precision.",
    "crate::kernels::mp::complex::MpComplex::zero", DimensionPolicy::Dimensionless,
    (precision: super::MpRealPrecision) -> super::mp::MpComplex,
    |precision: super::MpRealPrecision| { super::mp::MpComplex::zero(precision.0) }
);

binding!(
    mp_complex_one, "num_mp_complex_one", "numerical_special",
    "`num_mp_complex_one(precision)` constructs the exact real unit at the requested bit precision.",
    "crate::kernels::mp::complex::MpComplex::one", DimensionPolicy::Dimensionless,
    (precision: super::MpRealPrecision) -> super::mp::MpComplex,
    |precision: super::MpRealPrecision| { super::mp::MpComplex::one(precision.0) }
);

binding!(
    mp_complex_i, "num_mp_complex_i", "numerical_special",
    "`num_mp_complex_i(precision)` constructs the exact imaginary unit at the requested bit precision.",
    "crate::kernels::mp::complex::MpComplex::i", DimensionPolicy::Dimensionless,
    (precision: super::MpRealPrecision) -> super::mp::MpComplex,
    |precision: super::MpRealPrecision| { super::mp::MpComplex::i(precision.0) }
);

binding!(
    mp_complex_neg, "num_mp_complex_neg", "numerical_special",
    "`num_mp_complex_neg(value)` negates both coordinate intervals exactly.",
    "crate::kernels::mp::complex::MpComplex::neg", DimensionPolicy::Dimensionless,
    (value: super::mp::MpComplex) -> super::mp::MpComplex,
    |value: super::mp::MpComplex| { value.neg() }
);

binding!(
    mp_complex_conj, "num_mp_complex_conj", "numerical_special",
    "`num_mp_complex_conj(value)` conjugates a rectangle by negating its imaginary interval.",
    "crate::kernels::mp::complex::MpComplex::conj", DimensionPolicy::Dimensionless,
    (value: super::mp::MpComplex) -> super::mp::MpComplex,
    |value: super::mp::MpComplex| { value.conj() }
);

binding!(
    mp_complex_add, "num_mp_complex_add", "numerical_special",
    "`num_mp_complex_add(left, right, precision)` encloses the pairwise sums at the requested bit precision.",
    "crate::kernels::mp::complex::MpComplex::add", DimensionPolicy::Dimensionless,
    (left: super::mp::MpComplex, right: super::mp::MpComplex, precision: super::MpRealPrecision) -> Result<super::mp::MpComplex, OpError>,
    |left: super::mp::MpComplex, right: super::mp::MpComplex, precision: super::MpRealPrecision| { left.add(&right, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_complex_sub, "num_mp_complex_sub", "numerical_special",
    "`num_mp_complex_sub(left, right, precision)` encloses the pairwise differences at the requested bit precision.",
    "crate::kernels::mp::complex::MpComplex::sub", DimensionPolicy::Dimensionless,
    (left: super::mp::MpComplex, right: super::mp::MpComplex, precision: super::MpRealPrecision) -> Result<super::mp::MpComplex, OpError>,
    |left: super::mp::MpComplex, right: super::mp::MpComplex, precision: super::MpRealPrecision| { left.sub(&right, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_complex_mul, "num_mp_complex_mul", "numerical_special",
    "`num_mp_complex_mul(left, right, precision)` encloses the pairwise products at the requested bit precision.",
    "crate::kernels::mp::complex::MpComplex::mul", DimensionPolicy::Dimensionless,
    (left: super::mp::MpComplex, right: super::mp::MpComplex, precision: super::MpRealPrecision) -> Result<super::mp::MpComplex, OpError>,
    |left: super::mp::MpComplex, right: super::mp::MpComplex, precision: super::MpRealPrecision| { left.mul(&right, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_complex_div, "num_mp_complex_div", "numerical_special",
    "`num_mp_complex_div(left, right, precision)` encloses the pairwise quotients at the requested bit precision. The divisor must be separated from zero by the squared-modulus enclosure.",
    "crate::kernels::mp::complex::MpComplex::div", DimensionPolicy::Dimensionless,
    (left: super::mp::MpComplex, right: super::mp::MpComplex, precision: super::MpRealPrecision) -> Result<super::mp::MpComplex, OpError>,
    |left: super::mp::MpComplex, right: super::mp::MpComplex, precision: super::MpRealPrecision| { left.div(&right, precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_complex_recip, "num_mp_complex_recip", "numerical_special",
    "`num_mp_complex_recip(value, precision)` encloses the reciprocals at the requested bit precision. The rectangle must be separated from zero by the squared-modulus enclosure.",
    "crate::kernels::mp::complex::MpComplex::recip", DimensionPolicy::Dimensionless,
    (value: super::mp::MpComplex, precision: super::MpRealPrecision) -> Result<super::mp::MpComplex, OpError>,
    |value: super::mp::MpComplex, precision: super::MpRealPrecision| { value.recip(precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_complex_abs_sq, "num_mp_complex_abs_sq", "numerical_special",
    "`num_mp_complex_abs_sq(value, precision)` encloses the squared magnitudes using the rectangle’s nearest and farthest coordinates from zero.",
    "crate::kernels::mp::complex::MpComplex::abs_sq", DimensionPolicy::Dimensionless,
    (value: super::mp::MpComplex, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |value: super::mp::MpComplex, precision: super::MpRealPrecision| { value.abs_sq(precision.0).map_err(OpError::Domain) }
);

binding!(
    mp_complex_powi, "num_mp_complex_powi", "numerical_special",
    "`num_mp_complex_powi(value, exponent, precision)` encloses the integer power at the requested bit precision. The exponent must fit a signed 64-bit integer. Negative powers require separation from zero.",
    "crate::kernels::mp::complex::MpComplex::powi", DimensionPolicy::Dimensionless,
    (value: super::mp::MpComplex, exponent: i64, precision: super::MpRealPrecision) -> Result<super::mp::MpComplex, OpError>,
    |value: super::mp::MpComplex, exponent: i64, precision: super::MpRealPrecision| { value.powi(exponent, precision.0).map_err(OpError::Domain) }
);
