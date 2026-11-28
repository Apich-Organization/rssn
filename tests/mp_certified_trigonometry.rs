//! Independent reference intervals for certified MP trigonometry.

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::Signed;
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
    ball: &Ball,
    lower: &BigRational,
    upper: &BigRational,
) {
    let (lo, hi) = endpoints(ball);
    assert!(
        &lo <= lower && &hi >= upper,
        "{ball:?}: lower covered {}, upper covered {}",
        &lo <= lower,
        &hi >= upper,
    );
}

fn taylor(
    x: &BigRational,
    sine: bool,
) -> (BigRational, BigRational) {
    let magnitude = x.abs();
    let square = &magnitude * &magnitude;
    let mut term = if sine { magnitude } else { int(1) };
    let mut sum = term.clone();
    let offset = i64::from(sine);
    for k in 1..=48 {
        term = term * &square / int((2 * k - 1 + offset) * (2 * k + offset));
        if k % 2 == 0 {
            sum += &term;
        } else {
            sum -= &term;
        }
    }
    // For |x| <= 6, by degree 96 the omitted terms decrease in magnitude.
    let next = term * square / int((97 + offset) * (98 + offset));
    let lower = &sum - next;
    if sine && x.is_negative() {
        (-sum, -lower)
    } else {
        (lower, sum)
    }
}

fn atan_bounds(denominator: i64) -> (BigRational, BigRational) {
    let ratio = int(denominator).recip();
    let square = &ratio * &ratio;
    let mut power = ratio;
    let mut sum = int(0);
    for k in 0..120 {
        let term = &power / int(2 * k + 1);
        if k % 2 == 0 {
            sum += &term;
        } else {
            sum -= &term;
        }
        power *= &square;
    }
    let upper = &sum + power / int(241);
    (sum, upper)
}

#[test]
fn direct_and_fused_trig_enclose_independent_taylor_series() {
    for n in [1, -2, 3, -5, 6] {
        let x = Real::from_i64(n, 4).unwrap();
        let sine_bounds = taylor(&int(n), true);
        let cosine_bounds = taylor(&int(n), false);
        for prec in [4, 24, 128] {
            let pair = funcs::sin_cos(&x, prec).unwrap();
            encloses(&pair.sin, &sine_bounds.0, &sine_bounds.1);
            encloses(&pair.cos, &cosine_bounds.0, &cosine_bounds.1);
            assert_eq!(
                endpoints(&pair.sin),
                endpoints(&funcs::sin(&x, prec).unwrap())
            );
            assert_eq!(
                endpoints(&pair.cos),
                endpoints(&funcs::cos(&x, prec).unwrap())
            );
        }
    }
    let (sine_low, sine_high) = taylor(&int(1), true);
    let (cosine_low, cosine_high) = taylor(&int(1), false);
    let tangent = funcs::tan(&Real::one(1), 32).unwrap();
    encloses(
        &tangent,
        &(sine_low / cosine_high),
        &(sine_high / cosine_low),
    );
}

#[test]
fn large_arguments_have_certified_quarter_turn_reduction() {
    // Identity independent of the production Machin formula.
    let (a, b) = atan_bounds(2);
    let (c, d) = atan_bounds(3);
    let (half_pi_low, half_pi_high) = (int(2) * (a + c), int(2) * (b + d));
    let x = int(1_000_000);
    let k = int(636_620);
    let lo = &x - &k * half_pi_high;
    let hi = &x - k * half_pi_low;
    let scale = BigInt::from(1) << 208_u32;
    let centre = BigRational::new((lo.numer() * &scale).div_floor(lo.denom()), scale);
    let uncertainty = hi - &centre;
    assert!(uncertainty < BigRational::new(1.into(), BigInt::from(1) << 206_u32));

    // k is divisible by four, so the reduced sine and cosine agree with
    // those of the original argument. Derivatives have magnitude ≤ 1.
    let result = funcs::sin_cos(&Real::from_i64(1_000_000, 20).unwrap(), 96).unwrap();
    for (ball, sine) in [(result.sin, true), (result.cos, false)] {
        let (lower, upper) = taylor(&centre, sine);
        encloses(&ball, &(lower - &uncertainty), &(upper + &uncertainty));
    }
}

#[test]
fn interval_tangent_rejects_a_pole_and_sine_cosine_keep_closed_bounds() {
    let crossing = Ball::new(dyadic(3, -1, 2), dyadic(1, -3, 1)).unwrap();
    assert!(funcs::tan_ball(&crossing, 64).is_err());

    let positive = Ball::new(dyadic(1, -1, 1), dyadic(1, -1, 1)).unwrap();
    let (sin_lo, sin_hi) = taylor(&int(1), true);
    let (cos_lo, cos_hi) = taylor(&int(1), false);
    encloses(&funcs::sin_ball(&positive, 32).unwrap(), &int(0), &sin_hi);
    encloses(&funcs::cos_ball(&positive, 32).unwrap(), &cos_lo, &int(1));
    encloses(
        &funcs::tan_ball(&positive, 32).unwrap(),
        &int(0),
        &(sin_hi.clone() / cos_lo),
    );
    let from_text = parse::evaluate("sin(1)", 64).unwrap();
    encloses(&from_text, &sin_lo, &sin_hi);
    assert!(cos_hi < int(1));
}

#[test]
fn tiny_and_boundary_arguments_follow_the_bounded_domain() {
    let tiny = dyadic(1, -MAX_EXP, 1);
    let trig = funcs::sin_cos(&tiny, 1).unwrap();
    assert_eq!(trig.sin.mid.cmp_value(&tiny), 0);
    assert_eq!(trig.cos.mid.cmp_value(&Real::one(1)), 0);
    assert!(trig.sin.rad.sign() > 0 && trig.cos.rad.sign() > 0);
    assert!(funcs::sin(&dyadic(1, 20, 1), 64).is_ok());
    assert!(funcs::cos(&dyadic(3, 19, 2), 64).is_err());
    assert!(funcs::tan(&dyadic(1, MAX_EXP, 1), 64).is_err());
}

#[test]
fn native_session_registers_the_fused_and_interval_contracts() {
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
    let args = [dyadic(1, 0, 1), Value::int(4)];
    let Value::Record(pair) = call("num_mp_real_sin_cos", &args).unwrap() else {
        panic!("expected a sine-cosine record")
    };
    assert_eq!(
        pair.get("sin"),
        Some(&call("num_mp_real_sin", &args).unwrap())
    );
    assert_eq!(
        pair.get("cos"),
        Some(&call("num_mp_real_cos", &args).unwrap())
    );
    let crossing = record(&[("midpoint", dyadic(3, -1, 2)), ("radius", dyadic(1, -3, 1))]);
    assert!(call("num_mp_ball_tan", &[crossing, Value::int(16)]).is_err());
    assert!(call("num_mp_real_cos", &[dyadic(3, 19, 2), Value::int(16)]).is_err());

    let manifest = session.manifest();
    for (name, source) in [
        ("real_sin", "sin"),
        ("ball_sin", "sin_ball"),
        ("real_cos", "cos"),
        ("ball_cos", "cos_ball"),
        ("real_tan", "tan"),
        ("ball_tan", "tan_ball"),
        ("real_sin_cos", "sin_cos"),
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
