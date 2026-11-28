// Portable dyadic midpoint, radius and bounds for multiprecision intervals.
const MP_BALL_FIELDS: [&str; 2] = ["midpoint", "radius"];

impl Arg for mp::Ball {
    fn ty() -> TypeExpr {
        <Self as Output>::ty()
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let record = native_fields(value, &MP_BALL_FIELDS)?;
        let mid = mp::Real::decode(native_field(record, "midpoint")?, index)?;
        let rad = mp::Real::decode(native_field(record, "radius")?, index)?;
        Self::new(mid, rad).map_err(OpError::Domain)
    }
}

impl Output for mp::Ball {
    fn ty() -> TypeExpr {
        TypeExpr::Record(
            MP_BALL_FIELDS
                .into_iter()
                .map(|name| (name.to_owned(), <mp::Real as Output>::ty()))
                .collect(),
        )
    }

    fn encode(self) -> Result<Value, OpError> {
        Ok(Value::Record(RecordValue {
            fields: MP_BALL_FIELDS
                .into_iter()
                .map(str::to_owned)
                .zip([self.mid.encode()?, self.rad.encode()?])
                .collect(),
        }))
    }
}

binding!(
    mp_ball_div, "num_mp_ball_div", "numerical_special",
    "Encloses all quotients of two closed dimensionless intervals. Both arguments and the result have canonical dyadic midpoint and radius fields. The radius is nonnegative and can have a different bit precision from the midpoint. The divisor interval must exclude zero, including at its endpoints. An out-of-range exponent is a domain error.",
    "crate::kernels::mp::ball::Ball::div", DimensionPolicy::Dimensionless,
    (left: crate::kernels::mp::Ball, right: crate::kernels::mp::Ball,
     precision: super::MpRealPrecision) -> Result<crate::kernels::mp::Ball, OpError>,
    |left: crate::kernels::mp::Ball, right: crate::kernels::mp::Ball,
     precision: super::MpRealPrecision| {
        left.div(&right, precision.0).map_err(OpError::Domain)
    }
);


struct MpBallBounds {
    lower: mp::Real,
    upper: mp::Real,
}

record!(MpBallBounds {
    lower: mp::Real,
    upper: mp::Real,
});


binding!(
    mp_ball_new, "num_mp_ball_new", "numerical_special",
    "Creates a closed interval from canonical dyadic midpoint and radius. The radius must be nonnegative.",
    "crate::kernels::mp::ball::Ball::new", DimensionPolicy::Dimensionless,
    (midpoint: crate::kernels::mp::Real, radius: crate::kernels::mp::Real) -> Result<crate::kernels::mp::Ball, OpError>,
    |midpoint: crate::kernels::mp::Real, radius: crate::kernels::mp::Real| {
        crate::kernels::mp::Ball::new(midpoint, radius).map_err(OpError::Domain)
    }
);


binding!(
    mp_ball_point, "num_mp_ball_point", "numerical_special",
    "Creates a closed interval containing exactly the given dyadic.",
    "crate::kernels::mp::ball::Ball::point", DimensionPolicy::Dimensionless,
    (value: crate::kernels::mp::Real) -> crate::kernels::mp::Ball,
    |value: crate::kernels::mp::Real| {
        crate::kernels::mp::Ball::point(value)
    }
);


binding!(
    mp_ball_contains_zero, "num_mp_ball_contains_zero", "numerical_special",
    "Tests whether the closed interval contains zero, including either endpoint.",
    "crate::kernels::mp::ball::Ball::contains_zero", DimensionPolicy::Dimensionless,
    (interval: crate::kernels::mp::Ball) -> bool,
    |interval: crate::kernels::mp::Ball| {
        interval.contains_zero()
    }
);


binding!(
    mp_ball_width, "num_mp_ball_width", "numerical_special",
    "Returns twice the radius as a dyadic. An out-of-range exponent is a domain error.",
    "crate::kernels::mp::ball::Ball::width", DimensionPolicy::Dimensionless,
    (interval: crate::kernels::mp::Ball) -> Result<crate::kernels::mp::Real, OpError>,
    |interval: crate::kernels::mp::Ball| {
        interval.width().map_err(OpError::Domain)
    }
);


binding!(
    mp_ball_neg, "num_mp_ball_neg", "numerical_special",
    "Reflects a closed interval about zero and retains its radius.",
    "crate::kernels::mp::ball::Ball::neg", DimensionPolicy::Dimensionless,
    (interval: crate::kernels::mp::Ball) -> crate::kernels::mp::Ball,
    |interval: crate::kernels::mp::Ball| {
        interval.neg()
    }
);


binding!(
    mp_ball_abs, "num_mp_ball_abs", "numerical_special",
    "Returns a dyadic enclosure of the absolute values in the closed interval at the requested precision.",
    "crate::kernels::mp::ball::Ball::abs", DimensionPolicy::Dimensionless,
    (interval: crate::kernels::mp::Ball, precision: super::MpRealPrecision) -> Result<crate::kernels::mp::Ball, OpError>,
    |interval: crate::kernels::mp::Ball, precision: super::MpRealPrecision| {
        interval.abs(precision.0).map_err(OpError::Domain)
    }
);


binding!(
    mp_ball_endpoints, "num_mp_ball_endpoints", "numerical_special",
    "Returns the lower endpoint rounded down and upper endpoint rounded up at the requested precision.",
    "crate::kernels::mp::ball::Ball::endpoints", DimensionPolicy::Dimensionless,
    (interval: crate::kernels::mp::Ball, precision: super::MpRealPrecision) -> Result<super::MpBallBounds, OpError>,
    |interval: crate::kernels::mp::Ball, precision: super::MpRealPrecision| {
        interval.endpoints(precision.0).map(|(lower, upper)| super::MpBallBounds { lower, upper }).map_err(OpError::Domain)
    }
);

binding!(
    mp_ball_add, "num_mp_ball_add", "numerical_special",
    "`num_mp_ball_add(left, right, precision)` encloses every sum of two dyadic intervals. Precision sets the midpoint bit width; the radius includes midpoint rounding.",
    "crate::kernels::mp::ball::Ball::add", DimensionPolicy::Dimensionless,
    (left: crate::kernels::mp::Ball, right: crate::kernels::mp::Ball, precision: super::MpRealPrecision)
        -> Result<crate::kernels::mp::Ball, OpError>,
    |left: crate::kernels::mp::Ball, right: crate::kernels::mp::Ball, precision: super::MpRealPrecision| {
        left.add(&right, precision.0).map_err(OpError::Domain)
    }
);

binding!(
    mp_ball_sub, "num_mp_ball_sub", "numerical_special",
    "`num_mp_ball_sub(left, right, precision)` encloses every difference of two dyadic intervals. Precision sets the midpoint bit width; the radius includes midpoint rounding.",
    "crate::kernels::mp::ball::Ball::sub", DimensionPolicy::Dimensionless,
    (left: crate::kernels::mp::Ball, right: crate::kernels::mp::Ball, precision: super::MpRealPrecision)
        -> Result<crate::kernels::mp::Ball, OpError>,
    |left: crate::kernels::mp::Ball, right: crate::kernels::mp::Ball, precision: super::MpRealPrecision| {
        left.sub(&right, precision.0).map_err(OpError::Domain)
    }
);

binding!(
    mp_ball_mul, "num_mp_ball_mul", "numerical_special",
    "`num_mp_ball_mul(left, right, precision)` encloses every product of two dyadic intervals. Precision sets the midpoint bit width; the radius includes midpoint rounding.",
    "crate::kernels::mp::ball::Ball::mul", DimensionPolicy::Dimensionless,
    (left: crate::kernels::mp::Ball, right: crate::kernels::mp::Ball, precision: super::MpRealPrecision)
        -> Result<crate::kernels::mp::Ball, OpError>,
    |left: crate::kernels::mp::Ball, right: crate::kernels::mp::Ball, precision: super::MpRealPrecision| {
        left.mul(&right, precision.0).map_err(OpError::Domain)
    }
);
