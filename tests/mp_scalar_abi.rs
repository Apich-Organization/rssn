//! Canonical dyadic and interval operations through the public numerical ABI.
use num_bigint::BigInt;
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::value::RecordValue;
use rssn::kernel::value::VariantValue;

fn dyadic(
    mantissa: impl Into<BigInt>,
    exponent: i64,
    precision: u32,
) -> Value {
    record([
        ("mantissa", Value::int(mantissa)),
        ("exponent", Value::int(exponent)),
        ("precision", Value::int(precision)),
    ])
}

fn record(items: impl IntoIterator<Item = (&'static str, Value)>) -> Value {
    Value::Record(RecordValue {
        fields: items
            .into_iter()
            .map(|(key, value)| (key.into(), value))
            .collect(),
    })
}

fn option(payload: Option<i64>) -> Value {
    Value::Variant(VariantValue {
        tag: if payload.is_some() {
            "Some"
        } else {
            "None"
        }
        .into(),
        payload: payload.into_iter().map(Value::int).collect(),
    })
}

fn call(
    s: &Session,
    name: &str,
    args: &[Value],
) -> Value {
    s.call_value(name, args, &CallCtx::new())
        .unwrap_or_else(|e| panic!("{name}: {e:?}"))
}

#[test]
fn scalar_dyadics_preserve_exact_bits_and_round_once() {
    let s = Session::new();
    let negative = dyadic(-13, 0, 4);
    assert_eq!(
        call(&s, "num_mp_real_abs", std::slice::from_ref(&negative)),
        dyadic(13, 0, 4)
    );
    assert_eq!(call(&s, "num_mp_real_neg", &[dyadic(13, 0, 4)]), negative);
    assert_eq!(
        call(
            &s,
            "num_mp_real_is_integer",
            std::slice::from_ref(&negative)
        ),
        Value::Bool(true)
    );
    assert_eq!(
        call(&s, "num_mp_real_floor_log2", &[dyadic(-1, -1, 1)]),
        option(Some(-1))
    );
    assert_eq!(
        call(&s, "num_mp_real_floor_log2", &[dyadic(0, 0, 7)]),
        option(None)
    );
    assert_eq!(
        call(&s, "num_mp_real_is_zero", &[dyadic(0, 0, 7)]),
        Value::Bool(true)
    );
    assert_eq!(
        call(&s, "num_mp_real_is_integer", &[dyadic(7, -1, 3)]),
        Value::Bool(false)
    );
    assert_eq!(
        call(&s, "num_mp_real_trunc_to_bigint", &[dyadic(-7, -1, 3)]),
        Value::int(-3)
    );
    assert_eq!(
        call(&s, "num_mp_real_trunc_to_bigint", &[dyadic(1, 120, 1)]),
        Value::int(BigInt::from(1_u8) << 120_u32),
    );
    for (mode, root, requantized) in [("nearest", 11, 6), ("down", 11, 6), ("up", 12, 7)] {
        let mode = Value::String(mode.into());
        assert_eq!(
            call(
                &s,
                "num_mp_real_sqrt",
                &[dyadic(1, 1, 1), Value::int(4), mode.clone()]
            ),
            dyadic(root, -3, 4),
        );
        assert_eq!(
            call(
                &s,
                "num_mp_real_with_precision",
                &[dyadic(13, 0, 4), Value::int(3), mode]
            ),
            dyadic(requantized, 1, 3),
        );
    }
    for shift in [-262_144, -5, 0, 5, 262_144] {
        assert_eq!(
            call(
                &s,
                "num_mp_real_scale_pow2",
                &[dyadic(0, 0, 16), Value::int(shift)]
            ),
            dyadic(0, 0, 16),
        );
    }
    assert_eq!(
        call(
            &s,
            "num_mp_real_scale_pow2",
            &[dyadic(13, 0, 4), Value::int(-5)]
        ),
        dyadic(13, -5, 4),
    );
}

#[test]
fn ball_queries_keep_canonical_records_and_outward_bounds() {
    let s = Session::new();
    let point = call(&s, "num_mp_ball_point", &[dyadic(3, -1, 2)]);
    assert_eq!(
        point,
        record([("midpoint", dyadic(3, -1, 2)), ("radius", dyadic(0, 0, 2))])
    );
    let a = call(&s, "num_mp_ball_new", &[dyadic(3, -1, 2), dyadic(1, -1, 1)]);
    assert_eq!(
        call(&s, "num_mp_ball_width", std::slice::from_ref(&a)),
        dyadic(1, 0, 1),
    );
    assert_eq!(
        call(&s, "num_mp_ball_contains_zero", std::slice::from_ref(&a)),
        Value::Bool(false),
    );
    let reflected = call(&s, "num_mp_ball_neg", &[a]);
    assert_eq!(
        call(&s, "num_mp_ball_endpoints", &[reflected, Value::int(8)]),
        record([
            ("lower", dyadic(-128, -6, 8)),
            ("upper", dyadic(-128, -7, 8))
        ]),
    );
    let crossing = call(&s, "num_mp_ball_new", &[dyadic(1, 0, 1), dyadic(1, 1, 1)]);
    let magnitude = call(&s, "num_mp_ball_abs", &[crossing, Value::int(8)]);
    assert_eq!(
        call(
            &s,
            "num_mp_ball_contains_zero",
            std::slice::from_ref(&magnitude)
        ),
        Value::Bool(true),
    );
    assert_eq!(
        call(&s, "num_mp_ball_endpoints", &[magnitude, Value::int(8)]),
        record([("lower", dyadic(0, 0, 8)), ("upper", dyadic(192, -6, 8))]),
    );
}

#[test]
fn the_new_ops_reject_noncanonical_dyadics_and_invalid_requests() {
    let s = Session::new();
    let mode = || Value::String("nearest".into());
    let attempt = |op, args: &[Value]| s.call_value(op, args, &CallCtx::new());
    for (name, args) in [
        ("num_mp_real_is_zero", vec![dyadic(0, 1, 10)]),
        (
            "num_mp_real_sqrt",
            vec![dyadic(-1, 0, 1), Value::int(10), mode()],
        ),
        (
            "num_mp_real_sqrt",
            vec![dyadic(1, 0, 1), Value::int(0), mode()],
        ),
        (
            "num_mp_real_with_precision",
            vec![dyadic(1, 0, 1), Value::int(32_769), mode()],
        ),
        (
            "num_mp_real_scale_pow2",
            vec![dyadic(1, 262_144, 1), Value::int(1)],
        ),
        ("num_mp_ball_new", vec![dyadic(1, 0, 1), dyadic(-1, 0, 1)]),
    ] {
        assert!(
            matches!(attempt(name, &args), Err(OpError::Domain(_))),
            "{name}"
        );
    }
}
