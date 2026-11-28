//! Exact rational checks for directed multiprecision roots.

use num_bigint::BigInt;
use num_rational::BigRational;
use rssn::kernels::mp::Ball;
use rssn::kernels::mp::Real;
use rssn::kernels::mp::Round;
use rssn::kernels::mp::funcs;

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

fn encloses_root(
    ball: &Ball,
    radicand: &BigRational,
    degree: i32,
) {
    let midpoint = exact(&ball.mid);
    let radius = exact(&ball.rad);
    let low = &midpoint - &radius;
    let high = midpoint + radius;
    assert!(degree == 3 || low >= BigRational::default());
    assert!(low.pow(degree) <= *radicand && high.pow(degree) >= *radicand);
}

#[test]
fn point_roots_enclose_their_exact_radicands() {
    let one = Real::one(1);
    let two = dyadic(1, 1, 1);
    for precision in [1, 4, 17, 96] {
        encloses_root(&funcs::cbrt(&two, precision).unwrap(), &exact(&two), 3);
        encloses_root(
            &funcs::cbrt(&two.neg(), precision).unwrap(),
            &exact(&two.neg()),
            3,
        );
        encloses_root(
            &funcs::hypot(&one, &one, precision).unwrap(),
            &exact(&two),
            2,
        );
    }
}

#[test]
fn directed_roots_keep_the_correct_side_and_ties() {
    let minus_two = dyadic(-1, 1, 1);
    for (mode, expected) in [
        (Round::TowardPositive, dyadic(-1, 0, 1)),
        (Round::TowardNegative, dyadic(-3, -1, 2)),
        (Round::NearestEven, dyadic(-3, -1, 2)),
    ] {
        assert_eq!(minus_two.cbrt(2, mode).unwrap().cmp_value(&expected), 0);
    }
    // The square root of 25/16 is exactly halfway between 1 and 1.5.
    let halfway = dyadic(25, -4, 5);
    let nearest = halfway.sqrt(2, Round::NearestEven).unwrap();
    let upper = halfway.sqrt(2, Round::TowardPositive).unwrap();
    assert_eq!(nearest.cmp_value(&Real::one(1)), 0);
    assert_eq!(upper.cmp_value(&dyadic(3, -1, 2)), 0);
}

#[test]
fn exact_roots_support_large_exponents_and_keep_small_addends_sticky() {
    let target = dyadic(1, 180_000, 1);
    let cube = funcs::cbrt(&target, 64).unwrap();
    assert_eq!(cube.mid.cmp_value(&dyadic(1, 60_000, 1)), 0);
    assert!(cube.rad.is_zero());
    let a = dyadic(3, 200_000, 2);
    let b = dyadic(1, 200_002, 1);
    let norm = funcs::hypot(&a, &b, 64).unwrap();
    assert_eq!(norm.mid.cmp_value(&dyadic(5, 200_000, 3)), 0);
    assert!(norm.rad.is_zero());

    let tiny = dyadic(1, -200_000, 1);
    let rounded = funcs::hypot(&b, &tiny, 64).unwrap();
    assert_eq!(rounded.mid.cmp_value(&b), 0);
    assert!(rounded.rad.sign() > 0);
    encloses_root(&rounded, &(exact(&b).pow(2) + exact(&tiny).pow(2)), 2);
}

#[test]
fn interval_square_root_checks_closed_boundaries() {
    let one = Real::one(1);
    let interval = Ball::new(one.clone(), one.clone()).unwrap();
    let enclosure = interval.sqrt(64).unwrap();
    let zero = BigRational::default();
    let mid = exact(&enclosure.mid);
    let rad = exact(&enclosure.rad);
    assert!(&mid - &rad <= zero);
    assert!((mid + rad).pow(2) >= BigRational::from_integer(2.into()));
    let negative = Ball::new(one, dyadic(1, 1, 1)).unwrap();
    assert!(negative.sqrt(64).is_err());

    let exact = Ball::point(dyadic(1, 4, 1)).sqrt(64).unwrap();
    assert_eq!(exact.mid.cmp_value(&dyadic(1, 2, 1)), 0);
    assert!(exact.rad.is_zero());
}

#[test]
fn session_exposes_declared_root_and_dyadic_inspection_contracts() {
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
    let dyadic = |mantissa: Value, exponent: i64, precision: u32| {
        record(&[
            ("mantissa", mantissa),
            ("exponent", Value::int(exponent)),
            ("precision", Value::int(precision)),
        ])
    };
    let text = |s: &str| Value::String(s.into());
    let one = dyadic(Value::int(1), 0, 1);
    let two = dyadic(Value::int(1), 1, 1);
    let m: BigInt = (BigInt::from(1) << 64_u32) + 1_u32;
    let wide = dyadic(Value::int(m.clone()), -3, 65);
    let interval = record(&[("midpoint", wide.clone()), ("radius", one.clone())]);
    let session = Session::new();
    let context = CallCtx::new();
    let call = |name, args: &[Value]| session.call_value(name, args, &context);
    for (name, arguments, expected) in [
        ("num_mp_real_mantissa", vec![wide.clone()], Value::int(m)),
        ("num_mp_real_exponent", vec![wide.clone()], Value::int(-3)),
        ("num_mp_real_precision", vec![wide.clone()], Value::int(65)),
        ("num_mp_ball_midpoint", vec![interval.clone()], wide),
        ("num_mp_ball_radius", vec![interval], one.clone()),
        ("num_mp_real_rounding_mode", vec![text("ceil")], text("up")),
        (
            "num_mp_real_cbrt_in",
            vec![two.clone(), Value::int(2), text("down")],
            dyadic(Value::int(2), -1, 2),
        ),
    ] {
        assert_eq!(call(name, &arguments).unwrap(), expected, "{name}");
    }
    let Value::Record(limits) = call("num_mp_real_limits", &[]).unwrap() else {
        panic!("expected a limits record");
    };
    assert_eq!(limits.get("max_precision"), Some(&Value::int(32_768)));
    assert_eq!(limits.get("min_exponent"), Some(&Value::int(-262_144)));
    assert!(call("num_mp_real_cbrt", &[two, Value::int(0)]).is_err());
    assert!(call("num_mp_real_rounding_mode", &[text("unknown")]).is_err());
    let norm = call(
        "num_mp_real_hypot",
        &[one.clone(), one.clone(), Value::int(4)],
    )
    .unwrap();
    assert_ne!(
        call("num_mp_ball_radius", &[norm]).unwrap(),
        dyadic(Value::int(0), 0, 4)
    );
    let real_interval = record(&[("midpoint", one.clone()), ("radius", one)]);
    assert!(call("num_mp_ball_sqrt", &[real_interval, Value::int(8)]).is_ok());

    let manifest = session.manifest();
    for name in [
        "num_mp_real_mantissa",
        "num_mp_real_exponent",
        "num_mp_real_precision",
        "num_mp_ball_midpoint",
        "num_mp_ball_radius",
        "num_mp_real_limits",
        "num_mp_real_rounding_mode",
        "num_mp_real_cbrt",
        "num_mp_real_cbrt_in",
        "num_mp_real_hypot",
        "num_mp_ball_sqrt",
    ] {
        let op = manifest.op(name).unwrap();
        assert_eq!(op.signature_source, SignatureSource::Declared);
        assert_eq!(op.dim_rule, Some(DimRuleInfo::Dimensionless));
        assert_eq!(op.overloads.len(), 1);
    }
}
