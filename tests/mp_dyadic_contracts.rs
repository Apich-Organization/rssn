//! Exact multiprecision dyadic records through the public numerical ABI.

use num_bigint::BigInt;
use num_rational::BigRational;
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::manifest::SignatureSource;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::DimRuleInfo;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::RecordValue;
use rssn::kernels::mp::MAX_BITS;
use rssn::kernels::mp::MAX_EXP;

fn dyadic(
    mantissa: impl Into<BigInt>,
    exponent: i64,
    precision: u32,
) -> Value {
    Value::Record(RecordValue {
        fields: vec![
            ("mantissa".into(), Value::int(mantissa)),
            ("exponent".into(), Value::int(exponent)),
            ("precision".into(), Value::int(precision)),
        ],
    })
}

fn mode(name: &str) -> Value {
    Value::String(name.into())
}

#[test]
fn dyadic_construction_rounds_the_exact_input_once_in_each_mode() {
    let session = Session::new();
    let context = CallCtx::new();
    let round = |mantissa: BigInt, name| {
        session
            .call_value(
                "num_mp_real_round",
                &[
                    Value::int(mantissa),
                    Value::int(0),
                    Value::int(64),
                    mode(name),
                ],
                &context,
            )
            .unwrap()
    };

    let exact_half = (BigInt::from(1) << 120_usize) + (BigInt::from(1) << 56_usize);
    let kept = BigInt::from(1) << 63_usize;
    let next = &kept + 1_u32;
    for (rounding, positive_up, negative_down) in [
        ("nearest", false, false),
        ("zero", false, false),
        ("away", true, true),
        ("up", true, false),
        ("down", false, true),
    ] {
        let positive = if positive_up { &next } else { &kept };
        let negative = if negative_down {
            &next
        } else {
            &kept
        };
        assert_eq!(
            round(exact_half.clone(), rounding),
            dyadic(positive.clone(), 57, 64)
        );
        assert_eq!(round(-&exact_half, rounding), dyadic(-negative, 57, 64));
    }

    let odd_half = &exact_half + (BigInt::from(1) << 57_usize);
    assert_eq!(round(odd_half, "nearest"), dyadic(&kept + 2_u32, 57, 64));
}

#[test]
fn dyadic_arithmetic_keeps_bits_beyond_binary64() {
    let session = Session::new();
    let context = CallCtx::new();
    let binary = |name, left, right, precision, rounding| {
        session
            .call_value(
                name,
                &[left, right, Value::int(precision), mode(rounding)],
                &context,
            )
            .unwrap()
    };
    let compare = |left, right| {
        session
            .call_value("num_mp_real_compare", &[left, right], &context)
            .unwrap()
    };

    let power = BigInt::from(1) << 80_usize;
    let one_above = dyadic(&power + 1_u32, 0, 81);
    let base = dyadic(2, 79, 2);
    assert_eq!(compare(one_above.clone(), base.clone()), Value::int(1));
    assert_eq!(compare(base.clone(), one_above.clone()), Value::int(-1));
    assert_eq!(
        compare(dyadic(power.clone(), 0, 81), base.clone()),
        Value::int(0)
    );
    let Value::Record(mut reordered) = base.clone() else {
        panic!("the dyadic helper constructs a record");
    };
    reordered.fields.reverse();
    assert_eq!(
        compare(Value::Record(reordered), base.clone()),
        Value::int(0)
    );
    assert_eq!(
        binary(
            "num_mp_real_sub",
            one_above.clone(),
            base.clone(),
            4_u32,
            "nearest"
        ),
        dyadic(8, -3, 4)
    );
    assert_eq!(
        binary("num_mp_real_sub", base.clone(), base.clone(), 8_u32, "down"),
        dyadic(0, 0, 8)
    );
    assert_eq!(
        binary(
            "num_mp_real_add",
            one_above.clone(),
            base,
            82_u32,
            "nearest"
        ),
        dyadic((&power << 1_usize) + 1_u32, 0, 82)
    );
    assert_eq!(
        binary(
            "num_mp_real_mul",
            one_above,
            dyadic(3, 0, 2),
            82_u32,
            "nearest"
        ),
        dyadic((power + 1_u32) * 3_u32, 0, 82)
    );
    // An intermediate rounding to three bits would produce 20, which rounds to 16.
    assert_eq!(
        binary(
            "num_mp_real_mul",
            dyadic(7, 0, 3),
            dyadic(3, 0, 2),
            2,
            "nearest"
        ),
        dyadic(3, 3, 2)
    );
}

#[test]
fn dyadic_division_rounds_once_in_each_mode() {
    let session = Session::new();
    let context = CallCtx::new();
    let div = |left, right, precision: u32, rounding| {
        session.call_value(
            "num_mp_real_div",
            &[left, right, Value::int(precision), mode(rounding)],
            &context,
        )
    };
    let one = dyadic(1, 0, 1);
    let three = dyadic(3, 0, 2);
    for (rounding, positive, negative) in [
        ("nearest", 11, -11),
        ("zero", 10, -10),
        ("away", 11, -11),
        ("up", 11, -10),
        ("down", 10, -11),
    ] {
        assert_eq!(
            div(one.clone(), three.clone(), 4, rounding).unwrap(),
            dyadic(positive, -5, 4)
        );
        assert_eq!(
            div(dyadic(-1, 0, 1), three.clone(), 4, rounding).unwrap(),
            dyadic(negative, -5, 4)
        );
    }
    assert_eq!(
        div(dyadic(1025, 0, 11), one.clone(), 1, "up").unwrap(),
        dyadic(1, 11, 1)
    );
    assert_eq!(
        div(dyadic(5121, 0, 13), one.clone(), 2, "nearest").unwrap(),
        dyadic(3, 11, 2)
    );
    let leading = BigInt::from(1_u32) << (MAX_BITS - 1);
    let full_width = dyadic((&leading << 1_u32) - 1_u32, 0, MAX_BITS);
    assert_eq!(
        div(full_width.clone(), full_width, MAX_BITS, "nearest").unwrap(),
        dyadic(leading, 1 - i64::from(MAX_BITS), MAX_BITS)
    );
    assert_eq!(
        div(dyadic(0, 0, 1), three, 64, "away").unwrap(),
        dyadic(0, 0, 64)
    );
    assert!(matches!(
        div(one.clone(), dyadic(0, 0, 8), 64, "nearest"),
        Err(OpError::Domain(message)) if message.contains("division by zero")
    ));
    assert!(matches!(
        div(dyadic(1, MAX_EXP, 1), dyadic(1, -1, 1), 1, "nearest"),
        Err(OpError::Domain(message)) if message.contains("exponent out of range")
    ));
    assert!(div(one.clone(), one.clone(), 0, "nearest").is_err());
    assert!(div(one.clone(), one, 2, "invalid").is_err());
}

#[test]
fn dyadic_codec_rejects_invalid_and_dimensional_values() {
    let session = Session::new();
    let context = CallCtx::new();
    let valid = dyadic(1, 0, 1);
    for invalid in [
        dyadic(127, 0, 8),
        dyadic(0, 1, 8),
        dyadic(1, 0, 0),
        dyadic(1, 0, MAX_BITS + 1),
        dyadic(1, MAX_EXP + 1, 1),
        dyadic(1, -MAX_EXP - 1, 1),
        Value::Record(RecordValue {
            fields: vec![
                ("mantissa".into(), Value::int(1)),
                ("exponent".into(), Value::int(0)),
            ],
        }),
        Value::Record(RecordValue {
            fields: vec![
                ("mantissa".into(), Value::int(1)),
                ("exponent".into(), Value::int(0)),
                ("exponent".into(), Value::int(1)),
            ],
        }),
        Value::Record(RecordValue {
            fields: vec![
                (
                    "mantissa".into(),
                    Value::Dyn(DynValue {
                        value: Box::new(Value::int(1)),
                        dim: DimExpr::base("Length"),
                    }),
                ),
                ("exponent".into(), Value::int(0)),
                ("precision".into(), Value::int(1)),
            ],
        }),
    ] {
        assert!(
            session
                .call_value("num_mp_real_compare", &[invalid, valid.clone()], &context)
                .is_err()
        );
    }

    let round = |mantissa: BigInt, exponent, precision, rounding| {
        session.call_value(
            "num_mp_real_round",
            &[
                Value::int(mantissa),
                Value::int(exponent),
                Value::int(precision),
                mode(rounding),
            ],
            &context,
        )
    };
    assert!(round(BigInt::from(1) << (2 * MAX_BITS), 0_i64, 8_u32, "nearest").is_err());
    assert!(round(BigInt::from(0), MAX_EXP + 1, 8_u32, "nearest").is_err());
    assert!(round(BigInt::from(1), -MAX_EXP, 8_u32, "nearest").is_err());
    assert!(round(BigInt::from(1), 0, 8, "invalid").is_err());
}

fn ball(
    midpoint: Value,
    radius: Value,
) -> Value {
    Value::Record(RecordValue {
        fields: vec![("radius".into(), radius), ("midpoint".into(), midpoint)],
    })
}

fn field<'a>(
    value: &'a Value,
    key: &str,
) -> &'a Value {
    let Value::Record(record) = value else {
        panic!("expected a record");
    };
    &record
        .fields
        .iter()
        .find(|(name, _)| name == key)
        .unwrap()
        .1
}

fn integer(value: &Value) -> &BigInt {
    let Value::Int(number) = value else {
        panic!("expected an Int");
    };
    number
}

fn exact_dyadic(value: &Value) -> BigRational {
    let exponent = i64::try_from(integer(field(value, "exponent"))).unwrap();
    let mantissa = integer(field(value, "mantissa"));
    let shift = exponent.unsigned_abs();
    if exponent >= 0 {
        BigRational::from_integer(mantissa << shift)
    } else {
        BigRational::new(mantissa.clone(), BigInt::from(1_u32) << shift)
    }
}

#[test]
fn ball_division_encloses_exact_corners_and_its_result_is_composable() {
    let session = Session::new();
    let context = CallCtx::new();
    let div = |left, right, precision: u32| {
        session.call_value(
            "num_mp_ball_div",
            &[left, right, Value::int(precision)],
            &context,
        )
    };
    let one = dyadic(1, 0, 1);
    let three = dyadic(3, 0, 2);
    let zero = dyadic(0, 0, 1);
    let uncertainty = dyadic(1, -10, 1);
    let input = ball(one.clone(), uncertainty.clone());
    let negative_divisor = ball(dyadic(-3, 0, 2), uncertainty.clone());
    let positive_divisor = ball(three.clone(), uncertainty);
    let bound = |result: &Value, minimum, maximum| {
        let midpoint = exact_dyadic(field(result, "midpoint"));
        let radius = exact_dyadic(field(result, "radius"));
        assert!(&midpoint - &radius <= minimum);
        assert!(&midpoint + &radius >= maximum);
    };
    let rational = |n, d| BigRational::new(BigInt::from(n), BigInt::from(d));
    let positive = div(input.clone(), positive_divisor, 4).unwrap();
    bound(&positive, rational(1023, 3073), rational(1025, 3071));
    let composed = div(positive, ball(one.clone(), zero.clone()), 4).unwrap();
    bound(&composed, rational(1023, 3073), rational(1025, 3071));
    let negative = div(input, negative_divisor, 4).unwrap();
    bound(&negative, rational(-1025, 3071), rational(-1023, 3073));

    let point = div(ball(one, zero.clone()), ball(three, zero), 4).unwrap();
    assert_eq!(field(&point, "midpoint"), &dyadic(11, -5, 4));
    bound(&point, rational(1, 3), rational(1, 3));
}

#[test]
fn ball_division_validates_nested_records_and_rejects_zero_overlap() {
    let session = Session::new();
    let context = CallCtx::new();
    let one = dyadic(1, 0, 1);
    let zero = dyadic(0, 0, 1);
    let input = ball(one.clone(), zero.clone());
    let div = |right| {
        session.call_value(
            "num_mp_ball_div",
            &[input.clone(), right, Value::int(4)],
            &context,
        )
    };
    for midpoint in [zero.clone(), one.clone(), dyadic(-1, 0, 1)] {
        assert!(matches!(
            div(ball(midpoint, one.clone())),
            Err(OpError::Domain(message)) if message.contains("interval containing zero")
        ));
    }
    assert!(matches!(
        div(ball(one.clone(), dyadic(-1, 0, 1))),
        Err(OpError::Domain(message)) if message.contains("radius cannot be negative")
    ));
    assert!(div(ball(one.clone(), dyadic(1, 0, 2))).is_err());
    let measured = Value::Dyn(DynValue {
        value: Box::new(zero),
        dim: DimExpr::base("Length"),
    });
    assert!(div(ball(one, measured)).is_err());
}

#[test]
fn ball_division_declares_nested_dyadic_contract() {
    let manifest = Session::new().manifest();
    let dyadic = TypeExpr::Record(
        ["mantissa", "exponent", "precision"]
            .into_iter()
            .map(|name| (name.into(), TypeExpr::Int))
            .collect(),
    );
    let ball = TypeExpr::Record(vec![
        ("midpoint".into(), dyadic.clone()),
        ("radius".into(), dyadic),
    ]);
    let op = manifest.op("num_mp_ball_div").unwrap();
    let [signature] = op.overloads.as_slice() else {
        panic!("expected one overload");
    };
    assert_eq!(
        op.rust_path.as_deref(),
        Some("crate::kernels::mp::ball::Ball::div")
    );
    assert_eq!(op.signature_source, SignatureSource::Declared);
    assert_eq!(op.dim_rule, Some(DimRuleInfo::Dimensionless));
    assert_eq!(
        signature.params,
        vec![ball.clone(), ball.clone(), TypeExpr::Int]
    );
    assert_eq!(signature.ret, ball);
}

#[test]
fn dyadic_operators_declare_typed_dimensionless_contracts_and_native_paths() {
    let manifest = Session::new().manifest();
    let record = TypeExpr::Record(
        ["mantissa", "exponent", "precision"]
            .into_iter()
            .map(|name| (name.into(), TypeExpr::Int))
            .collect(),
    );
    for (name, method) in [
        ("round", "rounded"),
        ("compare", "cmp_value"),
        ("add", "add"),
        ("sub", "sub"),
        ("mul", "mul"),
        ("div", "div"),
    ] {
        let native = format!("crate::kernels::mp::real::Real::{method}");
        let op = manifest.op(&format!("num_mp_real_{name}")).unwrap();
        let [sig] = op.overloads.as_slice() else {
            panic!("expected a single declared overload for {name}");
        };
        assert_eq!(op.rust_path.as_deref(), Some(native.as_str()));
        assert_eq!(op.signature_source, SignatureSource::Declared);
        assert_eq!(op.dim_rule, Some(DimRuleInfo::Dimensionless));
        let params = match name {
            | "round" => {
                vec![
                    TypeExpr::Int,
                    TypeExpr::Int,
                    TypeExpr::Int,
                    TypeExpr::String,
                ]
            },
            | "compare" => vec![record.clone(), record.clone()],
            | _ => {
                vec![
                    record.clone(),
                    record.clone(),
                    TypeExpr::Int,
                    TypeExpr::String,
                ]
            },
        };
        assert_eq!(sig.params, params);
        assert_eq!(
            sig.ret,
            if name == "compare" {
                TypeExpr::Int
            } else {
                record.clone()
            }
        );
    }
}

#[test]
fn dyadic_conversions_keep_exact_bits_and_bound_inputs() {
    let session = Session::new();
    let context = CallCtx::new();
    let call = |name, args: &[Value]| session.call_value(name, args, &context);
    assert_eq!(
        call("num_mp_real_zero", &[Value::int(9)]).unwrap(),
        dyadic(0, 0, 9)
    );
    assert_eq!(
        call("num_mp_real_one", &[Value::int(9)]).unwrap(),
        dyadic(256, -8, 9)
    );
    let wide = (BigInt::from(1) << 80_usize) + 1_u32;
    assert_eq!(
        call(
            "num_mp_real_from_bigint",
            &[Value::int(wide.clone()), Value::int(81), mode("nearest")]
        )
        .unwrap(),
        dyadic(wide.clone(), 0, 81)
    );
    assert_eq!(
        call("num_mp_real_from_i64", &[Value::int(-11), Value::int(3)]).unwrap(),
        dyadic(-6, 1, 3)
    );
    for (mantissa, sign) in [(-1, -1), (0, 0), (1, 1)] {
        assert_eq!(
            call("num_mp_real_sign", &[dyadic(mantissa, 0, 1)]).unwrap(),
            Value::int(sign)
        );
    }
    for (name, args) in [
        ("num_mp_real_zero", vec![Value::int(0)]),
        ("num_mp_real_one", vec![Value::int(MAX_BITS + 1)]),
        (
            "num_mp_real_from_i64",
            vec![Value::int(wide), Value::int(81)],
        ),
        (
            "num_mp_real_from_bigint",
            vec![
                Value::int(BigInt::from(1) << (2 * MAX_BITS)),
                Value::int(81),
                mode("nearest"),
            ],
        ),
        (
            "num_mp_real_from_f64",
            vec![Value::Float(f64::INFINITY), Value::int(53)],
        ),
    ] {
        assert!(
            matches!(call(name, &args), Err(OpError::Domain(_))),
            "{name}"
        );
    }
}

#[test]
fn dyadic_binary64_conversion_rounds_subnormals_and_overflow_ties() {
    let session = Session::new();
    let context = CallCtx::new();
    let to_float = |value| session.call_value("num_mp_real_to_f64", &[value], &context);

    // Round trips across the binary64 exponent range, including the normal/subnormal boundary.
    for exponent in 0..2047_u64 {
        for fraction in [1, (1_u64 << 52) - 1] {
            for sign in [0, 1_u64 << 63] {
                let bits = sign | (exponent << 52) | fraction;
                let exact = session
                    .call_value(
                        "num_mp_real_from_f64",
                        &[Value::Float(f64::from_bits(bits)), Value::int(53)],
                        &context,
                    )
                    .unwrap();
                assert_eq!(to_float(exact).unwrap().as_f64().unwrap().to_bits(), bits);
            }
        }
    }
    for (mantissa, exponent, precision, expected) in [
        (-1_i64, -1075, 1, 1_u64 << 63),
        (3, -1075, 2, 2),
        ((1_i64 << 53) - 1, -1075, 53, f64::MIN_POSITIVE.to_bits()),
        ((1_i64 << 54) - 2, 970, 54, f64::MAX.to_bits()),
    ] {
        assert_eq!(
            to_float(dyadic(mantissa, exponent, precision))
                .unwrap()
                .as_f64()
                .unwrap()
                .to_bits(),
            expected
        );
    }
    assert!(matches!(
        to_float(dyadic((1_i64 << 54) - 1, 970, 54)),
        Err(OpError::Domain(_))
    ));
}

#[test]
fn ball_arithmetic_encloses_mixed_sign_corners_and_rejects_bad_records() {
    let session = Session::new();
    let context = CallCtx::new();
    let left = ball(dyadic(1, 0, 1), dyadic(1, 1, 1));
    let right = ball(dyadic(-3, 0, 2), dyadic(1, 0, 1));
    for (name, minimum, maximum) in [
        ("num_mp_ball_add", -5, 1),
        ("num_mp_ball_sub", 1, 7),
        ("num_mp_ball_mul", -12, 4),
    ] {
        let result = session
            .call_value(
                name,
                &[left.clone(), right.clone(), Value::int(3)],
                &context,
            )
            .unwrap();
        let midpoint = exact_dyadic(field(&result, "midpoint"));
        let radius = exact_dyadic(field(&result, "radius"));
        assert_eq!(
            field(field(&result, "midpoint"), "precision"),
            &Value::int(3)
        );
        assert!(
            &midpoint - &radius <= BigRational::from_integer(minimum.into()),
            "{name}"
        );
        assert!(
            &midpoint + &radius >= BigRational::from_integer(maximum.into()),
            "{name}"
        );
        assert!(
            session
                .call_value(name, &[result, right.clone(), Value::int(5)], &context)
                .is_ok(),
            "{name}"
        );
        let malformed = ball(dyadic(0, 0, 1), dyadic(-1, 0, 1));
        assert!(
            matches!(
                session.call_value(name, &[left.clone(), malformed, Value::int(3)], &context),
                Err(OpError::Domain(_))
            ),
            "{name}"
        );
    }
}
