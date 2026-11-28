//! Independent rational checks for the certified MP real logarithm routes.

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
    if x.exponent() < 0 {
        BigRational::new(
            x.mantissa().clone(),
            BigInt::from(1) << x.exponent().unsigned_abs(),
        )
    } else {
        BigRational::from_integer(x.mantissa() << x.exponent().unsigned_abs())
    }
}

fn logarithm_reference(n: i64) -> (BigRational, BigRational) {
    // -ln(1 - 1/n) = sum 1/(k * n^k). The positive tail after the
    // first omitted term has ratio below 1/2 for every n >= 2.
    let n = BigInt::from(n);
    let mut power = n.clone();
    let mut sum = int(0);
    for k in 1..=170 {
        sum += BigRational::new(1.into(), &power * k);
        power *= &n;
    }
    let upper = &sum + BigRational::new(2.into(), power * 171);
    (sum, upper)
}

fn encloses(
    value: &Ball,
    lower: &BigRational,
    upper: &BigRational,
) {
    let mid = exact(&value.mid);
    let rad = exact(&value.rad);
    assert!(&mid - &rad <= *lower && mid + rad >= *upper, "{value:?}");
}

#[test]
fn real_logarithms_enclose_independent_series() {
    let (ln2_low, ln2_high) = logarithm_reference(2);
    let (ln15_low, ln15_high) = logarithm_reference(3);
    let (ln125_low, ln125_high) = logarithm_reference(5);
    let log2_low = ln15_low / &ln2_high;
    let log2_high = ln15_high / &ln2_low;
    let log10_low = &ln2_low / (ln125_high + int(3) * &ln2_high);
    let log10_high = &ln2_high / (ln125_low + int(3) * &ln2_low);

    for prec in [1, 4, 16, 64, 128] {
        let two = dyadic(1, 1, 1);
        encloses(&funcs::ln(&two, prec).unwrap(), &ln2_low, &ln2_high);
        encloses(
            &funcs::ln(&two.scale_pow2(-2).unwrap(), prec).unwrap(),
            &(-&ln2_high),
            &(-&ln2_low),
        );
        encloses(
            &funcs::log2(&dyadic(3, -1, 2), prec).unwrap(),
            &log2_low,
            &log2_high,
        );
        encloses(&funcs::log10(&two, prec).unwrap(), &log10_low, &log10_high);
    }
    let parsed = parse::evaluate("ln(2)", 64).unwrap();
    encloses(&parsed, &ln2_low, &ln2_high);
}

#[test]
fn logarithms_preserve_relative_precision_near_one_and_exact_binary_powers() {
    let base: BigInt = BigInt::from(1) << 100_u32;
    let distance = BigRational::new(1.into(), base.clone());
    for x in [
        Real::rounded(&base + 1_u32, -100, 101, Round::NearestEven).unwrap(),
        Real::rounded(&base - 1_u32, -100, 100, Round::NearestEven).unwrap(),
    ] {
        let v = funcs::ln(&x, 64).unwrap();
        if x.floor_log2() == Some(0) {
            encloses(&v, &(&distance - &distance * &distance / int(2)), &distance);
        } else {
            encloses(&v, &(-&distance - &distance * &distance), &(-&distance));
        }
        assert!(exact(&v.rad) < &distance / int(1 << 30));
    }
    let binary = funcs::log2(&dyadic(1, -MAX_EXP, 1), 20).unwrap();
    assert_eq!(exact(&binary.mid), int(-MAX_EXP));
    assert!(binary.rad.is_zero());
}

#[test]
fn binary_exponential_encloses_algebraic_values_and_large_exact_powers() {
    for (sign, square) in [(1, int(2)), (-1, BigRational::new(1.into(), 2.into()))] {
        for prec in [1, 4, 64] {
            let result = funcs::exp2(&dyadic(sign, -1, 1), prec).unwrap();
            let mid = exact(&result.mid);
            let rad = exact(&result.rad);
            assert!((&mid - &rad).pow(2) <= square && (mid + rad).pow(2) >= square);
        }
    }
    let result = funcs::exp2(&dyadic(1, 17, 1), 64).unwrap();
    assert_eq!(result.mid.cmp_value(&dyadic(1, 1 << 17, 1)), 0);
    assert!(result.rad.is_zero());

    let tiny = funcs::exp2(&dyadic(-1, -MAX_EXP, 1), 64).unwrap();
    assert_eq!(exact(&tiny.mid), int(1));
    assert!(tiny.rad.sign() > 0);
}

#[test]
fn intervals_handle_endpoints_and_invalid_logarithm_arguments() {
    let signs = Ball::new(Real::zero(1), Real::one(1)).unwrap();
    let half = BigRational::new(1.into(), 2.into());
    encloses(&funcs::exp2_ball(&signs, 32).unwrap(), &half, &int(2));
    for log in [funcs::ln_ball, funcs::log2_ball, funcs::log10_ball] {
        assert!(log(&signs, 64).is_err());
        let touches_zero = Ball::new(Real::one(1), Real::one(1)).unwrap();
        assert!(log(&touches_zero, 64).is_err());
    }
    let real_interval = Ball::new(dyadic(3, -1, 2), dyadic(1, -1, 1)).unwrap();
    let (_, ln2_high) = logarithm_reference(2);
    encloses(
        &funcs::ln_ball(&real_interval, 64).unwrap(),
        &int(0),
        &ln2_high,
    );
    encloses(
        &funcs::log2_ball(&real_interval, 64).unwrap(),
        &int(0),
        &int(1),
    );
}

#[test]
fn decimal_scaling_uses_one_rounding_and_validates_wide_exponents() {
    let one_and_quarter = dyadic(5, -2, 3);
    let scale = |n, mode| exact(&funcs::scale_decimal_in(&one_and_quarter, n, 1, mode).unwrap());
    assert_eq!(scale(1, Round::NearestEven), int(16));
    assert_eq!(scale(1, Round::TowardNegative), int(8));
    assert_eq!(
        scale(-1, Round::NearestEven),
        BigRational::new(1.into(), 8.into())
    );
    assert!(funcs::scale_decimal(&one_and_quarter, i64::MIN, 1).is_err());
    assert!(funcs::scale_decimal(&one_and_quarter, i64::MAX, 1).is_err());
    assert!(
        funcs::scale_decimal(&Real::zero(1), i64::MIN, 1)
            .unwrap()
            .is_zero()
    );
}

#[test]
fn session_has_preserved_typed_real_and_interval_paths() {
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
    let one = dyadic(1, 0, 1);
    let zero = dyadic(0, 0, 1);
    let point = record(&[("midpoint", one.clone()), ("radius", zero.clone())]);
    let touches_zero = record(&[("midpoint", one.clone()), ("radius", one)]);
    let session = Session::new();
    let context = CallCtx::new();
    let call = |name: &str, args: &[Value]| session.call_value(name, args, &context);
    let precision = Value::int(4);
    let ln_one = call("num_mp_ball_ln", &[point, precision.clone()]).unwrap();
    let canonical_zero = dyadic(0, 0, 4);
    assert_eq!(
        ln_one,
        record(&[
            ("midpoint", canonical_zero.clone()),
            ("radius", canonical_zero),
        ])
    );
    assert!(call("num_mp_ball_log10", &[touches_zero, precision]).is_err());
    assert!(call("num_mp_real_log2", &[zero, Value::int(12)]).is_err());
    assert_eq!(
        call(
            "num_mp_real_scale_decimal",
            &[dyadic(5, -2, 3), Value::int(1), Value::int(1)]
        )
        .unwrap(),
        dyadic(1, 4, 1)
    );
    assert_eq!(
        call(
            "num_mp_real_scale_decimal_in",
            &[
                dyadic(5, -2, 3),
                Value::int(1),
                Value::int(1),
                Value::String("down".into()),
            ]
        )
        .unwrap(),
        dyadic(1, 3, 1)
    );

    let manifest = session.manifest();
    for (name, native) in [
        ("real_ln", "ln"),
        ("ball_ln", "ln_ball"),
        ("real_log2", "log2"),
        ("ball_log2", "log2_ball"),
        ("real_log10", "log10"),
        ("ball_log10", "log10_ball"),
        ("real_exp2", "exp2"),
        ("ball_exp2", "exp2_ball"),
        ("real_scale_decimal", "scale_decimal"),
        ("real_scale_decimal_in", "scale_decimal_in"),
    ] {
        let op = manifest.op(&format!("num_mp_{name}")).unwrap();
        assert_eq!(op.signature_source, SignatureSource::Declared);
        assert_eq!(op.dim_rule, Some(DimRuleInfo::Dimensionless));
        assert_eq!(
            op.rust_path,
            Some(format!("crate::kernels::mp::funcs::{native}"))
        );
    }
}
