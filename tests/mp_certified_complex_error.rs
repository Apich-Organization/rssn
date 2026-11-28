//! Independent reference values and interval regressions for complex erf and erfc.

use num_bigint::BigInt;
use num_rational::BigRational;
use rssn::kernels::mp::Ball;
use rssn::kernels::mp::MAX_EXP;
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

#[test]
fn taylor_range_contains_off_axis_references_and_operand_uncertainty() {
    let center = MpComplex::point(dyadic(1, -2), dyadic(1, -1));
    let point = center.erf(80).unwrap();
    check_decimal(
        &point.re,
        "35118665938138045452713964603416829251241528906974370858210345945027584441246181780958341286480258935",
        101,
    );
    check_decimal(
        &point.im,
        "57128662239433419331313105634450952245977715094887949576224592792780426208540798393464089927758871379",
        101,
    );

    let rect = MpComplex::new(
        Ball::new(center.re.mid.clone(), dyadic(1, -3)).unwrap(),
        Ball::new(center.im.mid.clone(), dyadic(1, -3)).unwrap(),
    );
    let image = rect.erf(80).unwrap();
    for (real, imag) in [
        (
            "17970269803383398219163629234155550409336003168927919466999430843295956087312398979796611123222842081",
            "60373536860991375686144923855025847833447266485361110186075606972766065434537278471455971935681285602",
        ),
        (
            "50726534151002978519768607262275407492843209693638470729843169815653448354017281191029117820280017412",
            "52099986231849488181328490784185454018867798440116645539477881336762328963969256302755533792071063111",
        ),
        (
            "31622022186173416421559443011603981941777286234685976025605755441789993621669305727923121738416062705",
            "41442012940442149554436765138651473025810568641764262457458795272377218022717807608696023224532477902",
        ),
        (
            "40187878974299296432199561961431277375184002902050396253784711534802606929630524084643667327935747395",
            "74626656943753562582865331596037672431104701330242810729290498774508407976793586234945713124710748900",
        ),
    ] {
        check_decimal(&image.re, real, 101);
        check_decimal(&image.im, imag, 101);
    }
}

#[test]
fn right_sector_complement_stays_informative_after_cancellation() {
    let real = MpComplex::point(dyadic(20, 0), Real::zero(1))
        .erfc(128)
        .unwrap();
    check_decimal(
        &real.re,
        "53958656116079009289349991679053456040882726709236052834701037849130187631566699333354307424023437103",
        276,
    );
    assert!(real.re.mid.sign() > 0);
    assert!(exact(&real.re.rad) * (BigInt::from(1_u32) << 100_u32) < exact(&real.re.mid));
    zero(&real.im);

    let moderate = MpComplex::point(dyadic(9, 0), Real::zero(1))
        .erfc(192)
        .unwrap();
    check_decimal(
        &moderate.re,
        "41370317465138102380539034673625245957101919859479624598059099019076203547681843597048990538589086393",
        137,
    );
    assert!(moderate.re.mid.sign() > 0);
    assert!(exact(&moderate.re.rad) * (BigInt::from(1_u32) << 150_u32) < exact(&moderate.re.mid));

    let center = MpComplex::point(dyadic(20, 0), dyadic(1, -2));
    let off_axis = center.erfc(128).unwrap();
    check_decimal(
        &off_axis.re,
        "-47798119656980263543056547338898896728947034995202394271223668131808059984826196212960434841861930735",
        276,
    );
    check_decimal(
        &off_axis.im,
        "31843858472046427947751324965763663317790445563926414645616432001199510722714059954415891846186526221",
        276,
    );

    let rect = MpComplex::new(
        Ball::new(center.re.mid, dyadic(1, -3)).unwrap(),
        Ball::new(center.im.mid, dyadic(1, -3)).unwrap(),
    );
    let image = rect.erfc(80).unwrap();
    for (real, imag, exponent) in [
        (
            "-73060768379983390327271070023813094218737171107967692743418730150028709298045003980604247886545469413",
            "42344452714734916869028037170804110673756516466931959536713363601030859196838545384254323359000683190",
            274,
        ),
        (
            "-30139439776152907708144448934191100685460151949562512628384001840776205723308646134580163445676007693",
            "22917522759008903669518836592110906790266710773173654757886348343313139886390666113994383232460575672",
            278,
        ),
        (
            "-47919831388281279500244165066246916737952206501259167590976415585473159836118498032374850925138110405",
            "-39490549910736213454489070131689054665626637596207394051808998088248432501940172856252707405571154619",
            276,
        ),
    ] {
        check_decimal(&image.re, real, exponent);
        check_decimal(&image.im, imag, exponent);
    }
}

#[test]
fn axes_and_extreme_exponents_preserve_analytic_limits() {
    let origin = MpComplex::zero(1);
    let ordinary = origin.erf(1).unwrap();
    zero(&ordinary.re);
    zero(&ordinary.im);
    let complement = origin.erfc(1).unwrap();
    assert_eq!(
        exact(&complement.re.mid),
        BigRational::from_integer(1.into())
    );
    zero(&complement.im);

    let wide_positive = MpComplex::new(
        Ball::new(dyadic(600, 0), dyadic(50, 0)).unwrap(),
        Ball::point(Real::zero(1)),
    );
    let complement = wide_positive.erfc(256).unwrap();
    assert!(complement.re.mid.is_zero());
    assert_eq!(complement.re.rad.exponent(), -MAX_EXP);
    zero(&complement.im);
    let negative = wide_positive.neg().erfc(256).unwrap();
    assert_eq!(exact(&negative.re.mid), BigRational::from_integer(2.into()));
    zero(&negative.im);
    let ordinary = wide_positive.neg().erf(256).unwrap();
    assert_eq!(
        exact(&ordinary.re.mid),
        BigRational::from_integer((-1).into())
    );
    zero(&ordinary.im);

    let tiny = MpComplex::point(dyadic(1, -200000), Real::zero(1));
    let bound = tiny.erfc(128).unwrap();
    assert_eq!(exact(&bound.re.mid), BigRational::from_integer(1.into()));
    assert!(!bound.re.rad.is_zero());
    zero(&bound.im);

    let imag_axis = MpComplex::new(
        Ball::point(Real::zero(1)),
        Ball::new(dyadic(1, -1), dyadic(1, -3)).unwrap(),
    );
    let imaginary_erf = imag_axis.erf(80).unwrap();
    zero(&imaginary_erf.re);
    assert!(imaginary_erf.im.mid.sign() > 0);
    assert!(
        MpComplex::point(Real::zero(1), dyadic(64, 0))
            .erf(64)
            .is_err()
    );
}

#[test]
fn session_declares_error_function_contracts_and_check_precision() {
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
    for (name, expected) in [
        ("num_mp_complex_erf", &origin),
        ("num_mp_complex_erfc", &one),
    ] {
        let op = manifest.op(name).unwrap();
        assert_eq!(op.signature_source, SignatureSource::Declared);
        assert_eq!(op.dim_rule, Some(DimRuleInfo::Dimensionless));
        assert_eq!(
            &session
                .call_value(name, &[origin.clone(), Value::int(1)], &ctx)
                .unwrap(),
            expected
        );
        assert!(
            session
                .call_value(name, &[origin.clone(), Value::int(0)], &ctx)
                .is_err()
        );
    }
}
