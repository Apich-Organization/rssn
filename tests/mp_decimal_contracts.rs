//! Exact reference checks for the multiprecision decimal interfaces.

use num_bigint::BigInt;
use num_rational::BigRational;
use rssn::kernels::mp::Ball;
use rssn::kernels::mp::MAX_EXP;
use rssn::kernels::mp::Real;
use rssn::kernels::mp::Round;
use rssn::kernels::mp::decimal;

fn dyadic(
    mantissa: i64,
    exponent: i64,
    precision: u32,
) -> Real {
    Real::rounded(mantissa.into(), exponent, precision, Round::NearestEven).unwrap()
}

fn exact(real: &Real) -> BigRational {
    let shift = real.exponent().unsigned_abs();
    if real.exponent() < 0 {
        BigRational::new(real.mantissa().clone(), BigInt::from(1_u32) << shift)
    } else {
        BigRational::from_integer(real.mantissa() << shift)
    }
}

#[test]
fn decimal_to_dyadic_rounds_once_on_both_sides_of_a_tie() {
    let above = "1.2500000000000000000000000000001";
    let below = "1.2499999999999999999999999999999";
    for (text, expected) in [
        ("1.25000000001", dyadic(3, -1, 2)),
        (above, dyadic(3, -1, 2)),
        ("1.25", dyadic(1, 0, 2)),
        (below, dyadic(1, 0, 2)),
    ] {
        let value = decimal::from_decimal(text, 2).unwrap();
        assert_eq!(value.cmp_value(&expected), 0, "{text}");
    }

    for (mode, positive, negative) in [
        (Round::NearestEven, 1, -1),
        (Round::TowardZero, 1, -1),
        (Round::AwayFromZero, 2, -2),
        (Round::TowardPositive, 2, -1),
        (Round::TowardNegative, 1, -2),
    ] {
        for (text, expected) in [("1.49", positive), ("-1.49", negative)] {
            let value = decimal::from_decimal_in(text, 1, mode).unwrap();
            assert_eq!(value.cmp_value(&dyadic(expected, 0, 1)), 0);
        }
    }
}

#[test]
fn formatting_preserves_decimal_rounding_and_exact_metadata() {
    let close_to_ten = dyadic(39, -2, 6);
    assert_eq!(decimal::to_decimal(&close_to_ten, 1).unwrap(), "10");
    assert_eq!(
        decimal::to_decimal_in(&close_to_ten.neg(), 1, Round::TowardPositive).unwrap(),
        "-9"
    );
    assert_eq!(
        decimal::to_decimal_in(&close_to_ten.neg(), 1, Round::TowardNegative).unwrap(),
        "-10"
    );
    assert_eq!(decimal::decimal_exponent(&close_to_ten), Some(0));
    assert_eq!(decimal::decimal_exponent(&dyadic(1, -4, 1)), Some(-2));
    assert_eq!(decimal::decimal_exponent(&Real::zero(16)), None);
    let detailed = dyadic(1025, -10, 11);
    assert_eq!(decimal::significant_digits(&detailed), 11);
    assert_eq!(decimal::to_decimal(&detailed, 11).unwrap(), "1.0009765625");
    assert_eq!(decimal::significant_digits(&dyadic(25, 2, 5)), 1);
    assert_eq!(decimal::significant_digits(&Real::zero(16)), 1);
}

#[test]
fn parsing_a_ball_contains_the_exact_decimal_and_certifies_rounding() {
    for (text, numerator) in [("0.1", 1), ("-0.1", -1)] {
        let ball = decimal::enclose_decimal(text, 12).unwrap();
        let expected = BigRational::new(numerator.into(), 10.into());
        let mid = exact(&ball.mid);
        let rad = exact(&ball.rad);
        assert!(&mid - &rad <= expected && expected <= mid + rad);
        assert_eq!(decimal::ball_to_decimal(&ball, 1).unwrap(), text);
    }
    let point = decimal::enclose_decimal("0.125", 12).unwrap();
    assert!(point.rad.is_zero());
    let wider = Ball::new(dyadic(1, 0, 8), dyadic(1, -1, 8)).unwrap();
    assert!(decimal::ball_to_decimal(&wider, 1).is_err());
    assert_eq!(
        rssn::kernels::mp::round_text("1.25000000000000000000000001", 2, Round::NearestEven)
            .unwrap(),
        "1.3"
    );
}

#[test]
fn formatted_values_round_trip_at_both_binary_exponent_limits() {
    for exponent in [-MAX_EXP, MAX_EXP] {
        let value = dyadic(1, exponent, 1);
        let text = decimal::to_decimal(&value, 24).unwrap();
        assert!((75_000..100_000).contains(&text.len()));
        let round = decimal::from_decimal(&text, 1).unwrap();
        assert_eq!(round.cmp_value(&value), 0, "{exponent}");
        let point = Ball::point(value);
        assert_eq!(decimal::ball_to_decimal(&point, 24).unwrap(), text);
    }
}

#[test]
fn malformed_or_oversized_decimal_input_is_rejected() {
    for text in [
        "",
        "NaN",
        "±1",
        "1.2.3",
        "1e",
        "2e100001",
        "2e-100001",
        "1e90000",
    ] {
        assert!(decimal::from_decimal(text, 128).is_err(), "{text}");
    }
    assert!(decimal::from_decimal(&"0".repeat(100_001), 4).is_err());
}

#[test]
fn decimal_session_routes_keep_dimensionless_signatures_and_errors() {
    use rssn::api::Session;
    use rssn::kernel::CallCtx;
    use rssn::kernel::Value;
    use rssn::kernel::manifest::SignatureSource;
    use rssn::kernel::signature::DimRuleInfo;

    let session = Session::new();
    let context = CallCtx::new();
    let call = |name, args: &[Value]| session.call_value(name, args, &context);
    let text = |s: &str| Value::String(s.into());
    let input = text("1.25000000001");
    let value = call("num_mp_decimal_from_text", &[input.clone(), Value::int(2)]).unwrap();
    assert_eq!(
        call("num_mp_decimal_to_text", &[value.clone(), Value::int(8)]).unwrap(),
        text("1.5")
    );
    assert_eq!(
        call(
            "num_mp_decimal_to_text_in",
            &[value.clone(), Value::int(1), text("up")],
        )
        .unwrap(),
        text("2")
    );
    let lower = call(
        "num_mp_decimal_from_text_in",
        &[input.clone(), Value::int(2), text("down")],
    )
    .unwrap();
    assert_ne!(value, lower);
    assert_eq!(
        call(
            "num_mp_decimal_significant_digits",
            std::slice::from_ref(&value)
        )
        .unwrap(),
        Value::int(2)
    );
    assert_eq!(
        call("num_mp_decimal_working_bits", &[Value::int(10)]).unwrap(),
        Value::int(50)
    );
    assert_eq!(
        call(
            "num_mp_decimal_round_text",
            &[input, Value::int(2), text("nearest")],
        )
        .unwrap(),
        text("1.3")
    );
    let ball = call("num_mp_decimal_enclose", &[text("0.125"), Value::int(12)]).unwrap();
    assert_eq!(
        call("num_mp_decimal_format_ball", &[ball, Value::int(3)]).unwrap(),
        text("0.125")
    );
    for (name, args) in [
        ("num_mp_decimal_working_bits", vec![Value::int(0)]),
        ("num_mp_decimal_working_bits", vec![Value::int(10_001)]),
        (
            "num_mp_decimal_from_text",
            vec![text("1e100001"), Value::int(3)],
        ),
        (
            "num_mp_decimal_from_text_in",
            vec![text("1.2"), Value::int(3), text("bad")],
        ),
    ] {
        assert!(call(name, &args).is_err(), "{name}");
    }
    let manifest = session.manifest();
    let ops: Vec<_> = manifest
        .ops
        .iter()
        .filter(|op| op.name.starts_with("num_mp_decimal_"))
        .collect();
    assert_eq!(ops.len(), 10);
    assert!(ops.iter().all(|op| {
        op.signature_source == SignatureSource::Declared
            && op.dim_rule == Some(DimRuleInfo::Dimensionless)
            && op.overloads.len() == 1
    }));
}
