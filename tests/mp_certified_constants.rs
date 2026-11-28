//! Independent interval checks for certified multiprecision constants.

use num_bigint::BigInt;
use num_rational::BigRational;
use rssn::kernels::mp::MAX_BITS;
use rssn::kernels::mp::Real;
use rssn::kernels::mp::Round;
use rssn::kernels::mp::consts;
use rssn::kernels::mp::parse;

fn int(n: i64) -> BigRational {
    BigRational::from_integer(n.into())
}

fn exact(x: &Real) -> BigRational {
    let shift = x.exponent().unsigned_abs();
    if x.exponent() < 0 {
        BigRational::new(x.mantissa().clone(), BigInt::from(1) << shift)
    } else {
        BigRational::from_integer(x.mantissa() << shift)
    }
}

fn atan_bounds(denominator: i64) -> (BigRational, BigRational) {
    const TERMS: i64 = 192;
    let z = int(denominator).recip();
    let square = &z * &z;
    let mut power = z;
    let mut sum = int(0);
    for k in 0..TERMS {
        let term = &power / int(2 * k + 1);
        if k % 2 == 0 {
            sum += term;
        } else {
            sum -= term;
        }
        power *= &square;
    }
    let upper = &sum + power / int(2 * TERMS + 1);
    (sum, upper)
}

fn pi_bounds() -> (BigRational, BigRational) {
    // arctan(1/2) + arctan(1/3) = arctan(1), independently of Machin's identity.
    let (a, b) = atan_bounds(2);
    let (c, d) = atan_bounds(3);
    (int(4) * (a + c), int(4) * (b + d))
}

fn ln2_bounds() -> (BigRational, BigRational) {
    const TERMS: u32 = 384;
    let mut sum = int(0);
    for n in 1..=TERMS {
        sum += BigRational::new(1.into(), BigInt::from(n) << n);
    }
    // The logarithmic series at 1/2 has positive terms with ratio below 1/2.
    let upper = &sum + BigRational::new(1.into(), BigInt::from(TERMS + 1) << TERMS);
    (sum, upper)
}

fn e_bounds() -> (BigRational, BigRational) {
    const TERMS: u32 = 100;
    let mut factorial = BigInt::from(1);
    let mut sum = int(1);
    for n in 1..=TERMS {
        factorial *= n;
        sum += BigRational::new(1.into(), factorial.clone());
    }
    let upper = &sum + BigRational::new(2.into(), factorial * (TERMS + 1));
    (sum, upper)
}

#[test]
fn certified_bounds_contain_independent_rational_intervals() {
    let (pi_low, pi_high) = pi_bounds();
    for (name, reference) in [
        ("pi", (pi_low.clone(), pi_high.clone())),
        ("tau", (int(2) * pi_low, int(2) * pi_high)),
        ("ln2", ln2_bounds()),
        ("e", e_bounds()),
    ] {
        // A high-precision cache entry can tighten the radius to the outward
        // rounding of the midpoint's error. The independent rational bounds
        // must resolve that much smaller rounding margin, not just 128 bits.
        consts::named(name, 512).unwrap();
        for prec in [1, 4, 31, 128] {
            let ball = consts::named_ball(name, prec).unwrap();
            let mid = exact(&ball.mid);
            let rad = exact(&ball.rad);
            assert!(
                &mid - &rad <= reference.0 && mid + rad >= reference.1,
                "{name} {prec}"
            );
            assert_eq!(ball.mid.cmp_value(&consts::named(name, prec).unwrap()), 0);
            assert_eq!(ball.mid.precision(), prec);
        }
    }
}

#[test]
fn algebraic_constants_bound_their_exact_roots() {
    for prec in [1, 7, 64, 128] {
        for (name, target) in [("sqrt2", int(2)), ("phi", int(1))] {
            let ball = consts::named_ball(name, prec).unwrap();
            let midpoint = exact(&ball.mid);
            let radius = exact(&ball.rad);
            let low = &midpoint - &radius;
            let high = midpoint + radius;
            let polynomial = |x: BigRational| {
                let square = &x * &x;
                if name == "phi" {
                    square - x
                } else {
                    square
                }
            };
            assert!(low >= BigRational::new(1.into(), 2.into()));
            assert!(polynomial(low) <= target && polynomial(high) >= target);
        }
    }
}

#[test]
fn fixed_point_pi_observes_the_final_bits_at_maximum_precision() {
    let (lo, hi) = pi_bounds();
    for bits in [0_u32, 8, 64, 128] {
        let factor = BigRational::from_integer(BigInt::from(1) << bits);
        let lower = (&lo * &factor).to_integer();
        let upper = (&hi * &factor).to_integer();
        assert_eq!(lower, upper);
        assert_eq!(consts::pi_scaled(bits).unwrap(), lower);
    }
    let pi = consts::pi_scaled(MAX_BITS).unwrap();
    // Independently certified by direct exact integer quotients with bounds
    // on each individual alternating Taylor term.
    assert_eq!(pi.iter_u64_digits().next(), Some(0xba38_209c_f746_ce76));
    assert!(consts::pi_scaled(MAX_BITS + 1).is_err());
}

#[test]
fn expression_literals_keep_the_constant_error_at_maximum_precision() {
    let value = parse::evaluate("pi", MAX_BITS).unwrap();
    assert!(!value.rad.is_zero());
    let (lo, hi) = pi_bounds();
    let v = parse::evaluate("pi", 64).unwrap();
    let mid = exact(&v.mid);
    let rad = exact(&v.rad);
    assert!(&mid - &rad <= lo && hi <= mid + rad);
    assert!(consts::named("sqrt3", 32).is_err());
    assert!(consts::named_ball("√2", 32).is_err());

    let e = consts::e(2).unwrap();
    let three = Real::from_i64(3, 2).unwrap();
    assert_eq!(e.cmp_value(&three), 0);
    assert_eq!(
        consts::ln2(1)
            .unwrap()
            .cmp_value(&Real::rounded(1.into(), -1, 1, Round::NearestEven).unwrap()),
        0
    );
}

#[test]
fn session_returns_typed_constants_enclosures_and_catalogs() {
    use rssn::api::Session;
    use rssn::kernel::CallCtx;
    use rssn::kernel::Value;
    use rssn::kernel::manifest::SignatureSource;
    use rssn::kernel::signature::DimRuleInfo;
    use rssn::kernel::signature::TypeExpr;

    let session = Session::new();
    let context = CallCtx::new();
    let call = |name: &str, args: &[Value]| session.call_value(name, args, &context);
    let text = |s: &str| Value::String(s.into());
    let precision = Value::int(80);
    for (key, name) in [
        ("pi", "pi"),
        ("e", "e"),
        ("ln2", "ln2"),
        ("sqrt2", "sqrt2"),
        ("phi", "phi"),
    ] {
        assert_eq!(
            call(
                &format!("num_mp_real_const_{key}"),
                std::slice::from_ref(&precision)
            )
            .unwrap(),
            call("num_mp_real_const_named", &[text(name), precision.clone()]).unwrap()
        );
    }
    assert_eq!(
        call("num_mp_real_const_pi_scaled", &[Value::int(8)]).unwrap(),
        Value::int(804_u32)
    );
    assert!(call("num_mp_real_const_pi_scaled", &[Value::int(32_769)]).is_err());
    assert!(
        call(
            "num_mp_real_const_named",
            &[text("sqrt3"), precision.clone()]
        )
        .is_err()
    );
    assert!(call("num_mp_real_const_e", &[Value::int(0)]).is_err());
    let Value::Record(ball) = call("num_mp_real_const_enclose", &[text("tau"), precision]).unwrap()
    else {
        panic!("expected a dyadic ball");
    };
    assert!(ball.get("midpoint").is_some() && ball.get("radius").is_some());

    for (name, expected) in [
        ("num_mp_real_const_names", consts::names()),
        ("num_mp_real_const_functions", parse::functions()),
    ] {
        let Value::Tensor(catalog) = call(name, &[]).unwrap() else {
            panic!("expected a name tensor for {name}");
        };
        assert_eq!(catalog.shape, [expected.len()]);
        assert_eq!(catalog.elem, TypeExpr::String);
        assert_eq!(
            catalog.data,
            expected.iter().map(|name| text(name)).collect::<Vec<_>>()
        );
    }
    let manifest = session.manifest();
    let ops: Vec<_> = manifest
        .ops
        .iter()
        .filter(|op| op.name.starts_with("num_mp_real_const_"))
        .collect();
    assert_eq!(ops.len(), 10);
    assert!(ops.iter().all(|op| {
        op.signature_source == SignatureSource::Declared
            && op.dim_rule == Some(DimRuleInfo::Dimensionless)
            && op.overloads.len() == 1
    }));
}
