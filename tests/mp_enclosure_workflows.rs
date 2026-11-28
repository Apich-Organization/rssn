//! Independent dyadic references for expression and endpoint enclosures.

use num_bigint::BigInt;
use num_rational::BigRational;
use rssn::kernels::mp::Ball;
use rssn::kernels::mp::Real;
use rssn::kernels::mp::Round;
use rssn::kernels::mp::funcs;
use rssn::kernels::mp::parse;
use std::cell::Cell;

fn int(n: i64) -> BigRational {
    BigRational::from_integer(n.into())
}

fn scale(exponent: u32) -> BigRational {
    BigRational::new(1.into(), BigInt::from(1) << exponent)
}

fn dyadic(n: i64) -> Real {
    Real::from_i64(n, 16).unwrap()
}

fn exact(x: &Real) -> BigRational {
    let shift = x.exponent().unsigned_abs();
    if x.exponent() < 0 {
        BigRational::new(x.mantissa().clone(), BigInt::from(1) << shift)
    } else {
        BigRational::from_integer(x.mantissa() << shift)
    }
}

fn bounds(ball: &Ball) -> (BigRational, BigRational) {
    let center = exact(&ball.mid);
    let radius = exact(&ball.rad);
    (&center - &radius, center + radius)
}

fn contains(
    ball: &Ball,
    point: &BigRational,
) {
    let (lo, hi) = bounds(ball);
    assert!(lo <= *point && hi >= *point, "{ball:?} missed {point}");
}

#[test]
fn parser_root_and_hypot_propagate_operand_uncertainty() {
    let root = parse::evaluate("sqrt((1 + 2^-96) - 1 + 2^-62)", 32).unwrap();
    let (below, above) = bounds(&root);
    let below = below.max(int(0));
    let square = scale(62) + scale(96);
    assert!(
        &below * &below <= square && &above * &above >= square,
        "{root:?}"
    );

    let cube = parse::evaluate("cbrt((1 + 2^-96) - 1)", 32).unwrap();
    contains(&cube, &scale(32));

    let norm = parse::evaluate("hypot((1 + 2^-96) - 1, 0)", 32).unwrap();
    contains(&norm, &scale(96));

    let exact = parse::evaluate("sqrt(81) + hypot(3, 4) + cbrt(-8)", 24).unwrap();
    assert_eq!(bounds(&exact), (int(12), int(12)));
    assert!(parse::evaluate("sqrt(-1)", 24).is_err());
}

#[test]
fn parser_borrows_utf8_spans_and_caps_recursive_grammar() {
    contains(
        &parse::evaluate(" 1e+3 + 2E-2 - 20e-3 ", 32).unwrap(),
        &int(1000),
    );
    assert!(parse::evaluate("π", 32).unwrap_err().contains("π"));
    assert!(parse::evaluate("7e+", 32).is_err());
    assert!(parse::evaluate(&format!("{}1{}", "sqrt(".repeat(100), ")".repeat(100)), 32).is_err());
    assert!(parse::evaluate(&format!("{}1", "1^".repeat(100)), 32).is_err());
    assert_eq!(
        bounds(&parse::evaluate(&format!("{}2", "-+".repeat(333)), 32).unwrap()),
        (int(-2), int(-2))
    );
    assert!(parse::evaluate("pow(1, 2, 3)", 32).is_err());
    assert_eq!(
        bounds(&parse::evaluate("2^3^2", 32).unwrap()),
        (int(512), int(512))
    );
}

#[test]
fn range_map_keeps_uncertain_endpoint_bounds_in_both_directions() {
    let interval = Ball::new(dyadic(1), dyadic(1)).unwrap();
    for increasing in [false, true] {
        let calls = Cell::new(0);
        let curve = |point: &Real, prec: u32| {
            calls.set(calls.get() + 1);
            let sign = if increasing { 1 } else { -1 };
            if point.cmp_value(&dyadic(2)) == 0 {
                return Ball::new(dyadic(sign), dyadic(7));
            }
            let mapped =
                point
                    .mul(point, prec, Round::NearestEven)?
                    .mul(point, prec, Round::NearestEven)?;
            let value = if increasing {
                mapped
            } else {
                mapped.neg()
            };
            Ok(Ball::point(value))
        };
        let range = funcs::map(&interval, &curve, increasing, 24).unwrap();
        assert_eq!(calls.get(), 2);
        let expected = if increasing {
            (int(0), int(8))
        } else {
            (int(-8), int(0))
        };
        assert_eq!(bounds(&range), expected);
    }
}

#[test]
fn directed_endpoint_map_and_lipschitz_include_midpoint_rounding() {
    let interval = Ball::new(dyadic(3), dyadic(3)).unwrap();
    let calls = Cell::new(0);
    let scale = |x: &Real, prec: u32, round| {
        calls.set(calls.get() + 1);
        x.mul(&dyadic(-2), prec, round)
    };
    let mapped = interval.map_monotone(&scale, false, 24).unwrap();
    assert_eq!(bounds(&mapped), (int(-12), int(0)));
    assert_eq!(calls.get(), 2);

    let widened = interval.lipschitz(dyadic(7), &dyadic(2), 2).unwrap();
    assert_eq!(widened.mid.precision(), 2);
    contains(&widened, &int(1));
    contains(&widened, &int(13));
    assert!(interval.lipschitz(dyadic(7), &dyadic(-1), 2).is_err());
}

#[test]
fn portable_enclosure_contracts_are_declared_and_reject_invalid_inputs() {
    use rssn::api::Session;
    use rssn::kernel::CallCtx;
    use rssn::kernel::Value;
    use rssn::kernel::manifest::SignatureSource;
    use rssn::kernel::signature::DimRuleInfo;
    use rssn::kernel::value::RecordValue;

    let record = |fields: &[(&str, Value)]| {
        Value::Record(RecordValue {
            fields: fields
                .iter()
                .map(|(name, value)| ((*name).into(), value.clone()))
                .collect(),
        })
    };
    let scalar = |mantissa: i64| {
        record(&[
            ("mantissa", Value::int(mantissa)),
            ("exponent", Value::int(0)),
            ("precision", Value::int(1)),
        ])
    };
    let session = Session::new();
    let context = CallCtx::new();
    let manifest = session.manifest();
    for name in ["num_mp_expression_enclose", "num_mp_ball_lipschitz"] {
        let contract = manifest.op(name).unwrap();
        assert_eq!(contract.signature_source, SignatureSource::Declared);
        assert_eq!(contract.dim_rule, Some(DimRuleInfo::Dimensionless));
    }
    let enclose = "num_mp_expression_enclose";
    let computed = session
        .call_value(
            enclose,
            &[Value::String("sqrt(1)".into()), Value::int(1)],
            &context,
        )
        .unwrap();
    let Value::Record(fields) = computed else {
        panic!("dyadic ball result");
    };
    assert!(
        fields
            .fields
            .iter()
            .any(|(name, value)| name == "midpoint" && *value == scalar(1))
    );
    assert!(
        session
            .call_value(
                enclose,
                &[Value::String("sqrt(-1)".into()), Value::int(32)],
                &context
            )
            .is_err()
    );
    assert!(
        session
            .call_value(
                enclose,
                &[Value::String("1".repeat(4097)), Value::int(1)],
                &context
            )
            .is_err()
    );
    assert!(
        session
            .call_value(
                enclose,
                &[Value::String("1".into()), Value::int(0)],
                &context
            )
            .is_err()
    );
    let interval = record(&[("midpoint", scalar(1)), ("radius", scalar(1))]);
    assert!(
        session
            .call_value(
                "num_mp_ball_lipschitz",
                &[interval, scalar(0), scalar(-1), Value::int(1)],
                &context
            )
            .is_err()
    );
}
