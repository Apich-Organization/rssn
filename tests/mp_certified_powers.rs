//! Exact reference checks for the portable regular, integer and rational power routes.

use num_bigint::BigInt;
use num_rational::BigRational;
use rssn::kernels::mp::Ball;
use rssn::kernels::mp::MAX_EXP;
use rssn::kernels::mp::Real;
use rssn::kernels::mp::Round;
use rssn::kernels::mp::funcs;
use rssn::kernels::mp::parse;

fn int(n: i64) -> BigRational {
    BigRational::from_integer(n.into())
}

fn dyadic(
    mantissa: i64,
    exponent: i64,
    precision: u32,
) -> Real {
    Real::rounded(mantissa.into(), exponent, precision, Round::NearestEven).unwrap()
}

fn exact(x: &Real) -> BigRational {
    let shift = x.exponent().unsigned_abs();
    if x.exponent() < 0 {
        BigRational::new(x.mantissa().clone(), BigInt::from(1) << shift)
    } else {
        BigRational::from_integer(x.mantissa() << shift)
    }
}

fn endpoints(ball: &Ball) -> (BigRational, BigRational) {
    let midpoint = exact(&ball.mid);
    let radius = exact(&ball.rad);
    (&midpoint - &radius, midpoint + radius)
}

fn encloses(
    ball: &Ball,
    lo: &BigRational,
    hi: &BigRational,
) {
    let (lower, upper) = endpoints(ball);
    assert!(
        &lower <= lo && &upper >= hi,
        "{ball:?}: lower covered {}, upper covered {}",
        &lower <= lo,
        &upper >= hi
    );
}

#[test]
fn signed_integer_powers_enclose_exact_rational_values() {
    for input in [dyadic(3, 0, 2), dyadic(-3, -1, 2)] {
        for n in [-7_i32, -2, 0, 3, 4, 63] {
            let expected = exact(&input).pow(n);
            for p in [1, 4, 64] {
                let out = funcs::powi(&input, i64::from(n), p).unwrap();
                encloses(&out, &expected, &expected);
                assert_eq!(out.mid.precision(), p);
            }
        }
    }
    let value = funcs::powi(&Real::one(1).neg(), i64::MIN, 1).unwrap();
    assert_eq!(endpoints(&value), (int(1), int(1)));
    let odd = Real::from_i64((1 << 20) + 1, 21).unwrap();
    let value = funcs::pow(&Real::one(1).neg(), &odd, 4).unwrap();
    assert_eq!(endpoints(&value), (int(-1), int(-1)));
    assert!(funcs::powi(&dyadic(1, 1, 1), i64::MAX, 64).is_err());
    assert!(funcs::powi(&Real::zero(1), -1, 64).is_err());
}

#[test]
fn huge_integer_exponents_and_tiny_input_offsets_use_bounded_work() {
    let base = Real::rounded(
        (BigInt::from(1) << 50_u32) + 1_u32,
        -50,
        51,
        Round::NearestEven,
    )
    .unwrap();
    let actual = funcs::powi(&base, 1_i64 << 50, 32).unwrap();
    let mut factorial = BigInt::from(1);
    let mut e_lower = int(1);
    for n in 1..=40 {
        factorial *= n;
        e_lower += BigRational::new(1.into(), factorial.clone());
    }
    let e_upper = &e_lower + BigRational::new(2.into(), factorial * 41);
    let h = BigRational::new(1.into(), BigInt::from(1) << 50_u32);
    let limit_delta = BigRational::new(1.into(), BigInt::from(1) << 29_u32);
    let (lo, hi) = endpoints(&actual);
    // ln(1+h)/h is between 1-h/2 and 1, and exp(-t) >= 1-t.
    assert!(lo >= &e_lower * (int(1) - &h / int(2)) - &limit_delta);
    assert!(hi <= &e_upper + &limit_delta);
    let near_one = Real::rounded(
        (BigInt::from(1) << 63_u32) + 1_u32,
        -63,
        64,
        Round::NearestEven,
    )
    .unwrap();
    let (lo, hi) = endpoints(&funcs::powi(&near_one, i64::MIN, 32).unwrap());
    assert!(lo >= e_upper.recip() - &limit_delta);
    assert!(hi <= e_lower.recip() + limit_delta);
    assert!(funcs::powi(&dyadic(1, MAX_EXP, 1), 2, 64).is_err());
}

#[test]
fn rational_powers_reduce_first_and_enforce_the_real_domain() {
    let x = dyadic(-1, 3, 1);
    for (n, d, value) in [
        (2, 6, int(-2)),
        (2, 3, int(4)),
        (-1, 3, BigRational::new((-1).into(), 2.into())),
    ] {
        let actual = funcs::pow_ratio(&x, n, d, 64).unwrap();
        encloses(&actual, &value, &value);
    }
    assert!(funcs::pow_ratio(&x, 2, 4, 64).is_err());
    assert!(funcs::pow_ratio(&Real::one(1), 0, 0, 64).is_err());
    assert!(funcs::pow_ratio(&Real::zero(1), -1, 3, 64).is_err());
    let value = funcs::pow_ratio(&Real::zero(1), 0, 3, 64).unwrap();
    assert_eq!(endpoints(&value), (int(1), int(1)));

    let (lo, hi) = endpoints(&funcs::pow(&dyadic(1, 1, 1), &dyadic(1, -1, 1), 4).unwrap());
    assert!(lo > int(0) && lo.pow(2) <= int(2) && hi.pow(2) >= int(2));
    let parse = parse::evaluate("pow(2, 0.5)", 64).unwrap();
    let (lo, hi) = endpoints(&parse);
    assert!(lo > int(0) && lo.pow(2) <= int(2) && hi.pow(2) >= int(2));
}

#[test]
fn intervals_include_zero_minima_and_negative_exponent_endpoint_order() {
    let crosses = Ball::new(dyadic(-1, -1, 1), dyadic(3, -1, 2)).unwrap();
    encloses(
        &funcs::powi_ball(&crosses, 2, 32).unwrap(),
        &int(0),
        &int(4),
    );
    encloses(
        &funcs::powi_ball(&crosses, 3, 32).unwrap(),
        &int(-8),
        &int(1),
    );
    assert!(funcs::powi_ball(&crosses, -1, 32).is_err());
    assert!(funcs::pow_ratio_ball(&crosses, 1, 2, 32).is_err());
    let rational_crossing = Ball::new(dyadic(-7, -1, 3), dyadic(9, -1, 4)).unwrap();
    encloses(
        &funcs::pow_ratio_ball(&rational_crossing, 1, 3, 32).unwrap(),
        &int(-2),
        &int(1),
    );
    encloses(
        &funcs::pow_ratio_ball(&rational_crossing, 2, 3, 32).unwrap(),
        &int(0),
        &int(4),
    );
    let positive = Ball::new(dyadic(3, 0, 2), Real::one(1)).unwrap();
    encloses(
        &funcs::powi_ball(&positive, -2, 32).unwrap(),
        &BigRational::new(1.into(), 16.into()),
        &BigRational::new(1.into(), 4.into()),
    );
    let range = Ball::new(dyadic(3, -2, 2), dyadic(1, -2, 1)).unwrap();
    encloses(
        &funcs::pow_ball(&positive, &range, 32).unwrap(),
        &int(2),
        &int(4),
    );
}

#[test]
fn session_exposes_the_six_declared_power_contracts() {
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
                .map(|(n, v)| ((*n).into(), v.clone()))
                .collect(),
        })
    };
    let dyadic = |mantissa: i64, exponent: i64, precision: u32| {
        record(&[
            ("mantissa", Value::int(mantissa)),
            ("exponent", Value::int(exponent)),
            ("precision", Value::int(precision)),
        ])
    };
    let session = Session::new();
    let context = CallCtx::new();
    let call = |name: &str, args: &[Value]| session.call_value(name, args, &context);
    let zero = dyadic(0, 0, 1);
    let one = dyadic(1, 0, 1);
    let neg_one = dyadic(-1, 0, 1);
    let point_one = record(&[("midpoint", one.clone()), ("radius", zero.clone())]);
    let point_neg_eight = record(&[("midpoint", dyadic(-1, 3, 1)), ("radius", zero)]);
    assert_eq!(
        call(
            "num_mp_real_powi",
            &[neg_one, Value::int(i64::MIN), Value::int(1)]
        )
        .unwrap(),
        point_one,
    );
    assert!(
        call(
            "num_mp_real_pow_ratio",
            &[
                dyadic(-1, 3, 1),
                Value::int(2),
                Value::int(6),
                Value::int(16)
            ]
        )
        .is_ok()
    );
    assert!(
        call(
            "num_mp_ball_pow_ratio",
            &[
                point_neg_eight,
                Value::int(2),
                Value::int(4),
                Value::int(16)
            ]
        )
        .is_err()
    );
    assert!(
        call(
            "num_mp_real_pow_ratio",
            &[one, Value::int(1), Value::int(0), Value::int(16)]
        )
        .is_err()
    );
    let manifest = session.manifest();
    for (name, source) in [
        ("real_pow", "pow"),
        ("ball_pow", "pow_ball"),
        ("real_powi", "powi"),
        ("ball_powi", "powi_ball"),
        ("real_pow_ratio", "pow_ratio"),
        ("ball_pow_ratio", "pow_ratio_ball"),
    ] {
        let op = manifest.op(&format!("num_mp_{name}")).unwrap();
        assert_eq!(op.signature_source, SignatureSource::Declared);
        assert_eq!(op.dim_rule, Some(DimRuleInfo::Dimensionless));
        assert_eq!(
            op.rust_path,
            Some(format!("crate::kernels::mp::funcs::{source}"))
        );
    }
}
