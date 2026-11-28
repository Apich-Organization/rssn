//! Independent reference values and interval regressions for complex Bessel J.

use num_bigint::BigInt;
use num_rational::BigRational;
use rssn::kernels::mp::Ball;
use rssn::kernels::mp::MpComplex;
use rssn::kernels::mp::Real;
use rssn::kernels::mp::Round;

fn dyadic(
    mantissa: i64,
    exponent: i64,
) -> Real {
    Real::rounded(mantissa.into(), exponent, 24, Round::NearestEven).unwrap()
}

fn exact(x: &Real) -> BigRational {
    let shift = x.exponent().unsigned_abs();
    if x.exponent() < 0 {
        BigRational::new(x.mantissa().clone(), BigInt::from(1) << shift)
    } else {
        BigRational::from_integer(x.mantissa() << shift)
    }
}

fn contains(
    ball: &Ball,
    lower: &BigRational,
    upper: &BigRational,
) {
    let center = exact(&ball.mid);
    let radius = exact(&ball.rad);
    assert!(
        &center - &radius <= *lower && center + radius >= *upper,
        "{ball:?}"
    );
}

fn check_decimal(
    ball: &Ball,
    floor_text: &str,
    exponent: u32,
) {
    let floor = floor_text.parse::<BigInt>().unwrap();
    let denominator = BigInt::from(10).pow(exponent);
    let lower = BigRational::new(floor.clone(), denominator.clone());
    let upper = BigRational::new(floor + 1, denominator);
    contains(ball, &lower, &upper);
}

fn zero(ball: &Ball) {
    assert!(ball.mid.is_zero() && ball.rad.is_zero(), "{ball:?}");
}

fn point(
    re: i64,
    im: i64,
    exponent: i64,
) -> MpComplex {
    MpComplex::point(dyadic(re, exponent), dyadic(im, exponent))
}

#[test]
fn integer_order_is_entire_and_retains_precision_under_cancellation() {
    let negative = point(-3, 0, 0).bessel_j(&point(3, 0, 0), 128).unwrap();
    check_decimal(
        &negative.re,
        "-30906272225525164361826019494683314942913593599305679435447464972",
        65,
    );
    zero(&negative.im);
    let reflected_order = point(3, 0, 0).bessel_j(&point(-3, 0, 0), 128).unwrap();
    check_decimal(
        &reflected_order.re,
        "-30906272225525164361826019494683314942913593599305679435447464972",
        65,
    );
    let oscillatory = point(160, 0, 0).bessel_j(&point(0, 0, 0), 128).unwrap();
    check_decimal(
        &oscillatory.re,
        "-3368731249490915730199619610613892618717810799866040713570050113",
        65,
    );
    assert!(oscillatory.re.rad.floor_log2().is_some_and(|x| x < -120));
    zero(&oscillatory.im);
}

#[test]
fn noninteger_complex_order_contains_point_and_rectangle_references() {
    let point_order = point(2, 1, -2);
    let argument = point(10, 3, -3);
    let result = argument.bessel_j(&point_order, 128).unwrap();
    check_decimal(
        &result.re,
        "66388976813297097902144385878828367466745328309686315989955626219",
        65,
    );
    check_decimal(
        &result.im,
        "-6528846532667772599611860242632204769988236320936770152163053411",
        65,
    );
    assert!(result.re.rad.floor_log2().is_some_and(|x| x < -120));
    assert!(result.im.rad.floor_log2().is_some_and(|x| x < -120));

    let radius = dyadic(1, -5);
    let order = MpComplex::new(
        Ball::new(point_order.re.mid, radius.clone()).unwrap(),
        Ball::new(point_order.im.mid, radius.clone()).unwrap(),
    );
    let argument = MpComplex::new(
        Ball::new(argument.re.mid, radius.clone()).unwrap(),
        Ball::new(argument.im.mid, radius).unwrap(),
    );
    let enclosed = argument.bessel_j(&order, 128).unwrap();
    for (real, imag) in [
        (
            "67009132561933000302340643490592140164896924939439289802788528995",
            "-6164891719272042355879346134168748842761193861163290489903265839",
        ),
        (
            "66582092648572521756606823026865994513080441733878444254534915459",
            "-7040351133617125878497267685293441314942713003606260368857030057",
        ),
        (
            "66007716290131065076463558904176278950935193581577578678227874441",
            "-5556985524303773380980643174958936398544653343660947045998339648",
        ),
        (
            "65522194403743511314517884184092509888427449721158517265630380535",
            "-6730493294740641808619202989553643304481618252962698132570620291",
        ),
    ] {
        check_decimal(&enclosed.re, real, 65);
        check_decimal(&enclosed.im, imag, 65);
    }
}

#[test]
fn origin_and_excluded_inputs_are_explicit() {
    let origin = point(0, 0, 0);
    let order_zero = point(0, 0, 0);
    let one = origin.bessel_j(&order_zero, 1).unwrap();
    assert_eq!(exact(&one.re.mid), BigRational::from_integer(1.into()));
    assert!(one.re.rad.is_zero());
    zero(&one.im);
    let half = point(1, 0, -1);
    let z = origin.bessel_j(&half, 128).unwrap();
    zero(&z.re);
    zero(&z.im);
    assert!(point(-1, 0, 0).bessel_j(&half, 128).is_err());
    assert!(point(1, 0, 0).bessel_j(&half, 513).is_err());
    assert!(point(257, 0, 0).bessel_j(&order_zero, 128).is_err());
    assert!(point(1, 0, 0).bessel_j(&point(513, 0, 0), 128).is_err());
}

#[test]
fn session_declares_bessel_contract_and_checks_precision() {
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
                .map(|(name, value)| ((*name).into(), value.clone()))
                .collect(),
        })
    };
    let scalar = |n: i64| {
        record(&[
            ("mantissa", Value::int(n)),
            ("exponent", Value::int(0)),
            ("precision", Value::int(1)),
        ])
    };
    let ball = |n| record(&[("midpoint", scalar(n)), ("radius", scalar(0))]);
    let origin = record(&[("real", ball(0)), ("imag", ball(0))]);
    let one = record(&[("real", ball(1)), ("imag", ball(0))]);
    let session = Session::new();
    let ctx = CallCtx::new();
    let manifest = session.manifest();
    let name = "num_mp_complex_bessel_j";
    let op = manifest.op(name).unwrap();
    assert_eq!(op.signature_source, SignatureSource::Declared);
    assert_eq!(op.dim_rule, Some(DimRuleInfo::Dimensionless));
    assert_eq!(
        session
            .call_value(name, &[origin.clone(), origin.clone(), Value::int(1)], &ctx)
            .unwrap(),
        one
    );
    assert!(
        session
            .call_value(name, &[origin.clone(), origin, Value::int(0)], &ctx)
            .is_err()
    );
}
