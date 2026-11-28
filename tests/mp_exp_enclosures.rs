//! Exact-rational checks for arbitrary-precision exponential enclosures.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::Signed;
use rssn::kernels::mp::Ball;
use rssn::kernels::mp::MAX_BITS;
use rssn::kernels::mp::MAX_EXP;
use rssn::kernels::mp::Real;
use rssn::kernels::mp::Round;
use rssn::kernels::mp::funcs;

fn integer(n: i64) -> BigRational {
    BigRational::from_integer(n.into())
}

fn exact(real: &Real) -> BigRational {
    let shift = real.exponent().unsigned_abs();
    if real.exponent() < 0 {
        BigRational::new(real.mantissa().clone(), BigInt::from(1) << shift)
    } else {
        BigRational::from_integer(real.mantissa() << shift)
    }
}

fn certified_unit_exponential(x: &BigRational) -> (BigRational, BigRational) {
    // For 0 <= x <= 1, the positive terms beyond degree 80 decrease
    // by a ratio below 1/2. Twice the first bounds the omitted tail.
    let mut term = integer(1);
    let mut sum = term.clone();
    for j in 1..=80 {
        term = term * x / integer(j);
        sum += &term;
    }
    let upper = &sum + term * x * integer(2) / integer(81);
    (sum, upper)
}

fn require_enclosure(
    actual: &Ball,
    expected_low: &BigRational,
    expected_high: &BigRational,
) {
    let mid = exact(&actual.mid);
    let rad = exact(&actual.rad);
    assert!(
        &mid - &rad <= *expected_low,
        "lower endpoint missed: {actual:?}"
    );
    assert!(
        mid + rad >= *expected_high,
        "upper endpoint missed: {actual:?}"
    );
}

#[test]
fn exponential_encloses_exact_rationals_before_and_after_rescaling() {
    let (lo_e, hi_e) = certified_unit_exponential(&integer(1));
    for (input, precision) in [(1_i32, 4), (100, 20), (-100, 20), (103, 128)] {
        let x = Real::from_i64(i64::from(input), 7).unwrap();
        let actual = funcs::exp(&x, precision).unwrap();
        let (low, high) = if input < 0 {
            (hi_e.pow(input), lo_e.pow(input))
        } else {
            (lo_e.pow(input), hi_e.pow(input))
        };
        require_enclosure(&actual, &low, &high);
        assert_eq!(actual.mid.precision(), precision);
        assert!(exact(&actual.rad) * integer(1 << 15) < high || precision == 4);
    }
}

#[test]
fn exponential_respects_dyadic_inputs_and_huge_precision() {
    for mantissa in [-7_i64, -3, 3, 7] {
        let x = Real::rounded(mantissa.into(), -3, 3, Round::NearestEven).unwrap();
        let (low, high) = certified_unit_exponential(&exact(&x).abs());
        let (low, high) = if mantissa < 0 {
            (high.recip(), low.recip())
        } else {
            (low, high)
        };
        for precision in [1, 17, 256] {
            require_enclosure(&funcs::exp(&x, precision).unwrap(), &low, &high);
        }
    }

    let tiny = Real::rounded(1.into(), -i64::from(MAX_BITS) - 12, 1, Round::NearestEven).unwrap();
    let x = exact(&tiny);
    let low = integer(1) + &x;
    let high = &low + &x * &x;
    let actual = funcs::exp(&tiny, MAX_BITS).unwrap();
    require_enclosure(&actual, &low, &high);
    assert_eq!(actual.mid.precision(), MAX_BITS);
}

#[test]
fn exponential_handles_zero_and_the_exponent_limit() {
    let zero = funcs::exp(&Real::zero(1), 64).unwrap();
    assert_eq!(exact(&zero.mid), integer(1));
    assert!(zero.rad.is_zero());
    assert_eq!(zero.mid.precision(), 64);

    let negative_near_floor = Real::from_i64(-181_704, 18).unwrap();
    let result = funcs::exp(&negative_near_floor, 1).unwrap();
    assert_eq!(result.mid.exponent(), -MAX_EXP);
    assert_eq!(result.rad.exponent(), -MAX_EXP);

    let huge = Real::rounded(1.into(), MAX_EXP, 1, Round::NearestEven).unwrap();
    assert!(funcs::exp(&huge, MAX_BITS).is_err());
    assert!(funcs::exp(&huge.neg(), MAX_BITS).is_err());
}

#[test]
fn interval_exponential_encloses_endpoints_and_the_decimal_route() {
    let zero = Real::zero(1);
    let one = Real::one(1);
    let interval = Ball::new(zero, one).unwrap();
    let (_, hi_e) = certified_unit_exponential(&integer(1));
    for precision in [4, 64, 256] {
        let result = funcs::exp_ball(&interval, precision).unwrap();
        require_enclosure(&result, &hi_e.recip(), &hi_e);
    }

    let point = Ball::point(Real::from_i64(1, 4).unwrap());
    let exact = funcs::exp(&point.mid, 16).unwrap();
    let ball = funcs::exp_ball(&point, 16).unwrap();
    assert_eq!(ball.mid.cmp_value(&exact.mid), 0);
    assert_eq!(ball.rad.cmp_value(&exact.rad), 0);

    let from_text = rssn::kernels::mp::parse::evaluate("exp(1/3)", 65).unwrap();
    let (low, high) = certified_unit_exponential(&BigRational::new(1.into(), 3.into()));
    require_enclosure(&from_text, &low, &high);
}

#[test]
fn exponential_numerical_abi_is_typed_and_preserves_enclosures() {
    use rssn::api::Session;
    use rssn::kernel::CallCtx;
    use rssn::kernel::Value;
    use rssn::kernel::manifest::SignatureSource;
    use rssn::kernel::signature::DimRuleInfo;
    use rssn::kernel::signature::TypeExpr;
    use rssn::kernel::value::RecordValue;

    let fields = |pairs: &[(&str, Value)]| {
        Value::Record(RecordValue {
            fields: pairs
                .iter()
                .map(|(name, value)| ((*name).into(), value.clone()))
                .collect(),
        })
    };
    let dyadic = |mantissa: i64, exponent: i64, precision: i64| {
        fields(&[
            ("mantissa", Value::int(mantissa)),
            ("exponent", Value::int(exponent)),
            ("precision", Value::int(precision)),
        ])
    };
    let one = dyadic(1, 0, 1);
    let interval = fields(&[("midpoint", one.clone()), ("radius", dyadic(0, 0, 1))]);
    let session = Session::new();
    let context = CallCtx::new();
    let call = |name, args: &[Value]| session.call_value(name, args, &context);
    let output = call("num_mp_real_exp", &[one.clone(), Value::int(4)]).unwrap();
    assert_eq!(
        output,
        call("num_mp_ball_exp", &[interval, Value::int(4)]).unwrap()
    );
    let Value::Record(output) = output else {
        panic!("expected a dyadic ball record")
    };
    assert_eq!(output.get("midpoint"), Some(&dyadic(11, -2, 4)));
    assert_eq!(output.get("radius"), Some(&dyadic(9, -8, 4)));
    assert!(call("num_mp_real_exp", &[one, Value::int(0)]).is_err());
    let invalid = fields(&[("midpoint", dyadic(0, 0, 1)), ("radius", dyadic(-1, 0, 1))]);
    assert!(call("num_mp_ball_exp", &[invalid, Value::int(4)]).is_err());

    let manifest = session.manifest();
    let real = TypeExpr::Record(
        ["mantissa", "exponent", "precision"]
            .map(|name| (name.into(), TypeExpr::Int))
            .into(),
    );
    let ball = TypeExpr::Record(vec![
        ("midpoint".into(), real.clone()),
        ("radius".into(), real.clone()),
    ]);
    for (name, input) in [("num_mp_real_exp", real), ("num_mp_ball_exp", ball.clone())] {
        let declared = manifest.op(name).unwrap();
        let [signature] = declared.overloads.as_slice() else {
            panic!("expected a declared signature for {name}");
        };
        assert_eq!(signature.params, [input, TypeExpr::Int]);
        assert_eq!(signature.ret, ball);
        assert_eq!(declared.signature_source, SignatureSource::Declared);
        assert_eq!(declared.dim_rule, Some(DimRuleInfo::Dimensionless));
    }
}
