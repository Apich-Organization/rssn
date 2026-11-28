//! Independent rational interval checks for certified real hyperbolic paths.

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
    let mid = exact(&ball.mid);
    let rad = exact(&ball.rad);
    (&mid - &rad, mid + rad)
}

fn encloses(
    value: &Ball,
    lower: &BigRational,
    upper: &BigRational,
) {
    let (lo, hi) = endpoints(value);
    assert!(
        &lo <= lower && &hi >= upper,
        "{value:?}: lower covered {}, upper covered {}",
        &lo <= lower,
        &hi >= upper
    );
}

fn exp_reference(x: BigRational) -> (BigRational, BigRational) {
    let mut sum = int(1);
    let mut term = sum.clone();
    for n in 1..=64 {
        term = term * &x / int(n);
        sum += &term;
    }
    // Beyond degree 64, the positive terms for x in [0, 2] decrease
    // with ratio less than one half.
    let upper = &sum + term * x * int(2) / int(65);
    (sum, upper)
}

fn sinh_from_exp(e: &BigRational) -> BigRational {
    (e - e.recip()) / int(2)
}

fn cosh_from_exp(e: &BigRational) -> BigRational {
    (e + e.recip()) / int(2)
}

fn tanh_from_exp(e: &BigRational) -> BigRational {
    let square = e * e;
    (&square - int(1)) / (square + int(1))
}

#[test]
fn point_functions_cover_independent_series_on_both_algorithm_paths() {
    for exact_x in [dyadic(1, -2, 1), dyadic(1, 0, 1), dyadic(1, 1, 1)] {
        let (lo, hi) = exp_reference(exact(&exact_x));
        for (evaluate, reference, odd) in [
            (
                funcs::sinh as fn(&Real, u32) -> Result<Ball, String>,
                sinh_from_exp as fn(&BigRational) -> BigRational,
                true,
            ),
            (funcs::cosh, cosh_from_exp, false),
            (funcs::tanh, tanh_from_exp, true),
        ] {
            let expected_low = reference(&lo);
            let expected_high = reference(&hi);
            for precision in [1, 4, 40, 128] {
                let value = evaluate(&exact_x, precision).unwrap();
                encloses(&value, &expected_low, &expected_high);
                assert_eq!(value.mid.precision(), precision);
                let negative = evaluate(&exact_x.neg(), precision).unwrap();
                if odd {
                    encloses(&negative, &(-&expected_high), &(-&expected_low));
                } else {
                    assert_eq!(endpoints(&value), endpoints(&negative));
                }
            }
        }
    }
}

#[test]
fn very_small_sinh_and_tanh_retain_relative_precision() {
    let input = dyadic(1, -120, 1);
    let x = exact(&input);
    let error = x.pow(3);
    let sinh = funcs::sinh(&input, 64).unwrap();
    let tanh = funcs::tanh(&input, 64).unwrap();
    encloses(&sinh, &x, &(&x + &error));
    encloses(&tanh, &(&x - &error), &x);
    assert!(exact(&sinh.rad) < &x / int(1 << 40));
    assert!(exact(&tanh.rad) < &x / int(1 << 40));
    let at_floor = dyadic(1, -MAX_EXP, 1);
    let sinh = funcs::sinh(&at_floor, 1).unwrap();
    let tanh = funcs::tanh(&at_floor, 1).unwrap();
    assert_eq!(sinh.mid.cmp_value(&at_floor), 0);
    assert_eq!(tanh.mid.cmp_value(&at_floor), 0);
}

#[test]
fn saturation_and_early_scaling_avoid_unnecessary_overflow() {
    let large = Real::from_i64(181_705, 18).unwrap();
    assert!(funcs::exp(&large.neg(), 17).is_err());
    for value in [
        funcs::sinh(&large, 1).unwrap(),
        funcs::cosh(&large, 1).unwrap(),
    ] {
        assert_eq!(value.mid.exponent(), MAX_EXP);
        assert!(value.rad.sign() > 0);
    }
    let extreme = dyadic(1, MAX_EXP, 1);
    let limit = BigRational::new(1.into(), BigInt::from(1) << 95_u32);
    for sign in [1, -1] {
        let input = if sign > 0 {
            extreme.clone()
        } else {
            extreme.neg()
        };
        let value = funcs::tanh(&input, 16).unwrap();
        assert_eq!(exact(&value.mid), int(sign));
        assert!(exact(&value.rad) <= limit);
    }
    assert!(funcs::sinh(&extreme, 64).is_err());
}

#[test]
fn interval_cosh_includes_the_interior_minimum() {
    let interval = Ball::new(dyadic(-1, -1, 1), dyadic(3, -1, 2)).unwrap();
    // The interval is [-2, 1].
    let (_, e2_hi) = exp_reference(int(2));
    let (_, e_hi) = exp_reference(int(1));
    encloses(
        &funcs::cosh_ball(&interval, 64).unwrap(),
        &int(1),
        &cosh_from_exp(&e2_hi),
    );
    encloses(
        &funcs::sinh_ball(&interval, 64).unwrap(),
        &(-sinh_from_exp(&e2_hi)),
        &sinh_from_exp(&e_hi),
    );
    encloses(
        &funcs::tanh_ball(&interval, 64).unwrap(),
        &(-tanh_from_exp(&e2_hi)),
        &tanh_from_exp(&e_hi),
    );
    let counted = parse::evaluate("cosh(1)", 64).unwrap();
    let (e_lo, _) = exp_reference(int(1));
    encloses(&counted, &cosh_from_exp(&e_lo), &cosh_from_exp(&e_hi));
}

#[test]
fn session_exposes_exact_hyperbolic_point_and_interval_signatures() {
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
    let expect_zero = record(&[("midpoint", zero.clone()), ("radius", zero.clone())]);
    let expect_one = record(&[("midpoint", one.clone()), ("radius", zero.clone())]);
    assert_eq!(
        call("num_mp_real_sinh", &[zero.clone(), Value::int(1)]).unwrap(),
        expect_zero
    );
    assert_eq!(
        call("num_mp_real_cosh", &[zero.clone(), Value::int(1)]).unwrap(),
        expect_one
    );
    assert!(call("num_mp_real_tanh", &[dyadic(1, MAX_EXP, 1), Value::int(16)]).is_ok());
    assert!(call("num_mp_real_sinh", &[zero.clone(), Value::int(0)]).is_err());
    let crosses_zero = record(&[("midpoint", zero), ("radius", one)]);
    for name in ["num_mp_ball_sinh", "num_mp_ball_cosh", "num_mp_ball_tanh"] {
        assert!(call(name, &[crosses_zero.clone(), Value::int(16)]).is_ok());
    }
    let manifest = session.manifest();
    for (key, source) in [
        ("real_sinh", "sinh"),
        ("ball_sinh", "sinh_ball"),
        ("real_cosh", "cosh"),
        ("ball_cosh", "cosh_ball"),
        ("real_tanh", "tanh"),
        ("ball_tanh", "tanh_ball"),
    ] {
        let op = manifest.op(&format!("num_mp_{key}")).unwrap();
        assert_eq!(op.signature_source, SignatureSource::Declared);
        assert_eq!(op.dim_rule, Some(DimRuleInfo::Dimensionless));
        assert_eq!(
            op.rust_path,
            Some(format!("crate::kernels::mp::funcs::{source}"))
        );
    }
}
