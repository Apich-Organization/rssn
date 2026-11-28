// Scoped arithmetic operations carrying exact dyadic enclosure records.

struct MpExpressionText(String);

impl Arg for MpExpressionText {
    fn ty() -> TypeExpr {
        TypeExpr::String
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let Value::String(text) = value else {
            return Err(OpError::Type {
                expected: Box::new(TypeExpr::String),
                got: Box::new(value.type_expr()),
                arg: index,
            });
        };
        if text.len() > mp::parse::MAX_EXPRESSION_TEXT {
            return Err(OpError::Domain("the expression is too long".into()));
        }
        Ok(Self(text.clone()))
    }
}

binding!(
    mp_expression_enclose, "num_mp_expression_enclose", "numerical_special",
    "Evaluates a bounded arithmetic expression to a dyadic enclosure. Invalid syntax, uncertifiable domains and out-of-range precision return errors.",
    "crate::kernels::mp::parse::evaluate", DimensionPolicy::Dimensionless,
    (text: super::MpExpressionText, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |text: super::MpExpressionText, precision: super::MpRealPrecision| {
        super::mp::parse::evaluate(&text.0, precision.0).map_err(OpError::Domain)
    }
);

binding!(
    mp_ball_lipschitz, "num_mp_ball_lipschitz", "numerical_special",
    "Given the exact mapped midpoint and a caller-supplied nonnegative Lipschitz bound, widens the radius by input radius times the bound and includes rounding of the midpoint.",
    "crate::kernels::mp::ball::Ball::lipschitz", DimensionPolicy::Dimensionless,
    (interval: super::mp::Ball, mapped_midpoint: super::mp::Real, bound: super::mp::Real, precision: super::MpRealPrecision) -> Result<super::mp::Ball, OpError>,
    |interval: super::mp::Ball, mapped_midpoint: super::mp::Real, bound: super::mp::Real, precision: super::MpRealPrecision| {
        interval.lipschitz(mapped_midpoint, &bound, precision.0).map_err(OpError::Domain)
    }
);
