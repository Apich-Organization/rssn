//! Independent reference values and interval regressions for the complex Dirichlet functions.

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
fn zeta_continues_into_the_critical_strip_and_left_half_plane() {
    let zeta_two = point(2, 0, 0).riemann_zeta(160).unwrap();
    check_decimal(
        &zeta_two.re,
        "164493406684822643647241516664602518921894990120679843773555822937000747040320087",
        80,
    );
    zero(&zeta_two.im);
    for (order, (real, imag)) in [
        (
            point(4, 113, -3),
            (
                "124172918521375103466730858028533506471402210663611630156008204327317574687330",
                "-760673133561305551064304749353703391143979315521862591503906558871352278050427",
            ),
        ),
        (
            point(-1, 4, -1),
            (
                "22809497171652632980496113661822963199718695465599006634157462073864918465426889",
                "-14452917173371359641989033762685346144289816409075083960150524279989383786076658",
            ),
        ),
    ] {
        let value = order.riemann_zeta(160).unwrap();
        check_decimal(&value.re, real, 80);
        check_decimal(&value.im, imag, 80);
        assert!(value.re.rad.floor_log2().is_some_and(|bits| bits < -140));
        assert!(value.im.rad.floor_log2().is_some_and(|bits| bits < -140));
    }
    let center = point(-1, 4, -1);
    let rad = dyadic(1, -6);
    let order = MpComplex::new(
        Ball::new(center.re.mid, rad.clone()).unwrap(),
        Ball::new(center.im.mid, rad).unwrap(),
    );
    let value = order.riemann_zeta(128).unwrap();
    for (real, imag) in [
        (
            "22845922041307172222779176851372736982124639153003838346114484824765109278112668",
            "-13977520705938844195102410384407433284348081055562415178797425122236976295683821",
        ),
        (
            "22770260180509406837430059649667067692275357883582044201571093508810524671185529",
            "-14933660699137279994953887870303176854539840656445290934141192087843622526220896",
        ),
    ] {
        check_decimal(&value.re, real, 80);
        check_decimal(&value.im, imag, 80);
    }
}

#[test]
fn polylog_complex_orders_and_near_unit_arguments_are_enclosed() {
    let argument = point(14, 1, -4);
    let order = point(2, -1, -3);
    let value = argument.polylog(&order, 160).unwrap();
    check_decimal(
        &value.re,
        "341537565172495773404281018268908266179836528171352822508535118324240923484084326",
        80,
    );
    check_decimal(
        &value.im,
        "235596764887093743577207190556772255432864722733961241887725087575775497264806686",
        80,
    );
    assert!(value.re.rad.floor_log2().is_some_and(|bits| bits < -145));
    let value = point(7, 0, -3).polylog(&point(-5, 0, 0), 128).unwrap();
    check_decimal(
        &value.re,
        "2116805600000000000000000000000000000000000000000000000000000000000000000000000000000000",
        80,
    );
    zero(&value.im);

    let rad = dyadic(1, -6);
    let argument = MpComplex::new(
        Ball::new(argument.re.mid, rad.clone()).unwrap(),
        Ball::new(argument.im.mid, rad.clone()).unwrap(),
    );
    let order = MpComplex::new(
        Ball::new(order.re.mid, rad.clone()).unwrap(),
        Ball::new(order.im.mid, rad).unwrap(),
    );
    let value = argument.polylog(&order, 128).unwrap();
    for (real, imag) in [
        (
            "393345880217287386512870936297790446622135270734712392398644976264878876037196489",
            "251629712494692291900665801233161603203914443570315487426473209450452885785599347",
        ),
        (
            "301283284210320306329877848895051592361522018640095398262868307521923774183172366",
            "221426719603437168858660862749839297895388265030920411993419417531011468132297800",
        ),
    ] {
        check_decimal(&value.re, real, 80);
        check_decimal(&value.im, imag, 80);
    }
}

#[test]
fn exact_orders_and_domain_errors_are_explicit() {
    let zero_order = point(0, 0, 0);
    let one = point(1, 0, 0);
    let half = point(1, 0, -1);
    let geom = half.polylog(&zero_order, 128).unwrap();
    assert_eq!(exact(&geom.re.mid), BigRational::from_integer(1.into()));
    assert!(geom.re.rad.is_zero());
    zero(&geom.im);
    let trivial = point(-129, 0, 0).riemann_zeta(128).unwrap_err();
    assert!(trivial.contains("real order"));
    let trivial = point(-128, 0, 0).riemann_zeta(128).unwrap();
    zero(&trivial.re);
    zero(&trivial.im);
    assert!(one.riemann_zeta(128).is_err());
    assert!(one.polylog(&half, 64).is_err());
    assert!(half.polylog(&half, 513).is_err());
}

#[test]
fn session_declares_dirichlet_contracts_and_checks_precision() {
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
    let session = Session::new();
    let ctx = CallCtx::new();
    let manifest = session.manifest();
    let one = record(&[("real", ball(1)), ("imag", ball(0))]);
    for name in ["num_mp_complex_riemann_zeta", "num_mp_complex_polylog"] {
        let op = manifest.op(name).unwrap();
        assert_eq!(op.signature_source, SignatureSource::Declared);
        assert_eq!(op.dim_rule, Some(DimRuleInfo::Dimensionless));
    }
    assert!(
        session
            .call_value(
                "num_mp_complex_riemann_zeta",
                &[one.clone(), Value::int(96)],
                &ctx
            )
            .is_err()
    );
    assert!(
        session
            .call_value(
                "num_mp_complex_polylog",
                &[one, origin.clone(), Value::int(96)],
                &ctx
            )
            .is_err()
    );
    assert_eq!(
        session
            .call_value(
                "num_mp_complex_polylog",
                &[origin.clone(), origin.clone(), Value::int(1)],
                &ctx
            )
            .unwrap(),
        origin
    );
}
