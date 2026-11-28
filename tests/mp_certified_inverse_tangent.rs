//! Independent rational reference checks for certified real and rectangle argument.

use num_bigint::BigInt;
use num_rational::BigRational;
use rssn::kernels::mp::Ball;
use rssn::kernels::mp::MAX_EXP;
use rssn::kernels::mp::MpComplex;
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
    let radius = exact(&ball.rad);
    (&mid - &radius, mid + radius)
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
        &hi >= upper
    );
}

fn atan_reference(
    n: i64,
    d: i64,
) -> (BigRational, BigRational) {
    let ratio = BigRational::new(n.into(), d.into());
    let square = &ratio * &ratio;
    let mut power = ratio;
    let mut sum = int(0);
    for k in 0..160 {
        let term = &power / int(2 * k + 1);
        if k % 2 == 0 {
            sum += &term;
        } else {
            sum -= &term;
        }
        power *= &square;
    }
    // For 0 <= ratio < 1 the first omitted positive term bounds the error.
    let upper = &sum + power / int(321);
    (sum, upper)
}

fn quarter_pi() -> (BigRational, BigRational) {
    // atan(1/2) + atan(1/3) = atan(1), independent of the Machin formula.
    let (lo, hi) = atan_reference(1, 2);
    let (a, b) = atan_reference(1, 3);
    (lo + a, hi + b)
}

#[test]
fn real_inverse_tangent_and_axis_values_enclose_independent_bounds() {
    let (quarter_lo, quarter_hi) = quarter_pi();
    let (a, b) = atan_reference(1, 4);
    let (c, d) = atan_reference(2, 3);
    let zero = Real::zero(1);
    let one = Real::one(1);
    let two = dyadic(1, 1, 1);
    let three = dyadic(3, 0, 2);
    for precision in [1, 4, 40, 80] {
        encloses(
            &funcs::atan(&one, precision).unwrap(),
            &quarter_lo,
            &quarter_hi,
        );
        encloses(
            &funcs::atan(&dyadic(-1, -2, 1), precision).unwrap(),
            &(-&b),
            &(-&a),
        );
        encloses(
            &funcs::atan2(&one, &zero, precision).unwrap(),
            &(int(2) * &quarter_lo),
            &(int(2) * &quarter_hi),
        );
        encloses(
            &funcs::atan2(&zero, &one.neg(), precision).unwrap(),
            &(int(4) * &quarter_lo),
            &(int(4) * &quarter_hi),
        );
        encloses(
            &funcs::atan2(&three, &two.neg(), precision).unwrap(),
            &(int(2) * &quarter_lo + &c),
            &(int(2) * &quarter_hi + &d),
        );
        encloses(
            &funcs::atan2(&three.neg(), &two, precision).unwrap(),
            &(&c - int(2) * &quarter_hi),
            &(&d - int(2) * &quarter_lo),
        );
    }
    assert!(funcs::atan2(&zero, &zero, 64).is_err());
    let from_text = parse::evaluate("atan2(0,-1)", 64).unwrap();
    encloses(&from_text, &(int(4) * quarter_lo), &(int(4) * quarter_hi));
}

#[test]
fn exact_ratios_avoid_intermediate_overflow_and_keep_near_axis_precision() {
    let huge = dyadic(1, MAX_EXP, 1);
    let tiny = dyadic(1, -MAX_EXP, 1);
    let (quarter_lo, quarter_hi) = quarter_pi();
    let loose = BigRational::new(1.into(), BigInt::from(1) << 140_u32);
    encloses(
        &funcs::atan2(&huge, &tiny, 64).unwrap(),
        &(int(2) * quarter_lo - loose),
        &(int(2) * quarter_hi),
    );
    let close = dyadic(1, -130, 1);
    let q = exact(&close);
    let value = funcs::atan(&close, 64).unwrap();
    encloses(&value, &(&q - q.pow(3) / int(3)), &q);
    assert!(exact(&value.rad) < &q / int(1 << 30));
    assert_eq!(funcs::atan(&tiny, 1).unwrap().mid.cmp_value(&tiny), 0);
}

#[test]
fn rectangle_argument_encloses_corners_and_enforces_the_branch_cut() {
    let zero = Real::zero(1);
    let one = Real::one(1);
    let rectangle_x = Ball::new(zero.clone(), one.clone()).unwrap();
    let rectangle_y = Ball::point(one.clone());
    let angle = funcs::atan2_ball(&rectangle_y, &rectangle_x, 64).unwrap();
    let (quarter_lo, quarter_hi) = quarter_pi();
    encloses(&angle, &quarter_lo, &(int(3) * &quarter_hi));
    let complex = MpComplex::new(rectangle_x.clone(), rectangle_y.clone());
    assert_eq!(endpoints(&angle), endpoints(&complex.arg(64).unwrap()));
    assert!(funcs::atan2_ball(&rectangle_x, &rectangle_x, 64).is_err());
    assert!(funcs::atan2_ball(&rectangle_x, &Ball::point(one.neg()), 64).is_err());

    let negative_axis = Ball::new(one.neg(), dyadic(1, -2, 1)).unwrap();
    let point_y = Ball::point(zero);
    let axis = funcs::atan2_ball(&point_y, &negative_axis, 64).unwrap();
    encloses(&axis, &(int(4) * &quarter_lo), &(int(4) * &quarter_hi));
    assert!(MpComplex::new(negative_axis, point_y).arg(64).is_err());
    let atans = funcs::atan_ball(&rectangle_x, 64).unwrap();
    encloses(&atans, &(-&quarter_hi), &quarter_hi);
}

#[test]
fn session_keeps_the_declared_coordinate_order_and_cut_policy() {
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
    let point_y = record(&[("midpoint", one.clone()), ("radius", zero.clone())]);
    let point_zero = record(&[("midpoint", zero.clone()), ("radius", zero.clone())]);
    let point_negative = record(&[("midpoint", neg_one), ("radius", zero.clone())]);
    let interval_x = record(&[("midpoint", zero.clone()), ("radius", one.clone())]);
    let rectangle = record(&[("real", interval_x.clone()), ("imag", point_y.clone())]);
    let quadrant = call(
        "num_mp_real_atan2",
        &[one.clone(), zero.clone(), Value::int(8)],
    )
    .unwrap();
    let Value::Record(axis) = quadrant else {
        panic!("expected a real ball")
    };
    assert_ne!(axis.get("radius"), Some(&dyadic(0, 0, 8)));
    assert_eq!(
        call("num_mp_complex_arg", &[rectangle, Value::int(16)]).unwrap(),
        call(
            "num_mp_ball_atan2",
            &[point_y, interval_x.clone(), Value::int(16)]
        )
        .unwrap()
    );
    assert!(
        call(
            "num_mp_ball_atan2",
            &[interval_x, point_negative.clone(), Value::int(16)]
        )
        .is_err()
    );
    assert!(
        call(
            "num_mp_ball_atan2",
            &[point_zero, point_negative, Value::int(16)]
        )
        .is_ok()
    );
    assert!(call("num_mp_real_atan2", &[zero.clone(), zero, Value::int(16)]).is_err());
    assert!(call("num_mp_real_atan", &[one, Value::int(0)]).is_err());

    let manifest = session.manifest();
    for (key, source) in [
        ("real_atan", "funcs::atan"),
        ("ball_atan", "funcs::atan_ball"),
        ("real_atan2", "funcs::atan2"),
        ("ball_atan2", "funcs::atan2_ball"),
        ("complex_arg", "complex::MpComplex::arg"),
    ] {
        let op = manifest.op(&format!("num_mp_{key}")).unwrap();
        assert_eq!(op.signature_source, SignatureSource::Declared);
        assert_eq!(op.dim_rule, Some(DimRuleInfo::Dimensionless));
        assert_eq!(op.rust_path, Some(format!("crate::kernels::mp::{source}")));
    }
}
