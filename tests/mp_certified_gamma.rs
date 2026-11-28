//! Independent reference values and interval regressions for the complex gamma family.

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

type GammaOperation = fn(&MpComplex, u32) -> Result<MpComplex, String>;
const OPERATIONS: [GammaOperation; 3] = [MpComplex::ln_gamma, MpComplex::gamma, MpComplex::digamma];

#[test]
fn point_values_enclose_independent_references_and_log_gamma_windings() {
    for (z, references) in [
        (
            point(10, 3, -3),
            [
                (
                    "-18028447798999403118405971657283029841861806566945193492089926514",
                    "-7411331731361256105328496041325531975749488300285230306988587321",
                ),
                (
                    "83274035190189125575451550651379096169896451222498816001677798250",
                    "-6183039842048889927037681115758467118638124667225075038936795137",
                ),
                (
                    "-14045274813250388032920330500380056633703414785282235781062618909",
                    "42639313090147813228914288817045562753895783149756411387723762930",
                ),
            ],
        ),
        (
            point(-18, 3, -3),
            [
                (
                    "-25526939408942423909331582279896905570684135670662034478299663835",
                    "-835138557485022757503380027331955307731779588937691612369273393099",
                ),
                (
                    "-36964845662284093280735400412080038517555952388288364677718748917",
                    "-68083195030300489181787322508712077936162124575715879378452886690",
                ),
                (
                    "161614424011940823027465457581391449613250805443213388301195804888",
                    "295151393444439550185258087684769546679364715891157545086619709319",
                ),
            ],
        ),
        (
            point(-18, -3, -3),
            [
                (
                    "-25526939408942423909331582279896905570684135670662034478299663835",
                    "835138557485022757503380027331955307731779588937691612369273393098",
                ),
                (
                    "-36964845662284093280735400412080038517555952388288364677718748917",
                    "68083195030300489181787322508712077936162124575715879378452886689",
                ),
                (
                    "161614424011940823027465457581391449613250805443213388301195804888",
                    "-295151393444439550185258087684769546679364715891157545086619709320",
                ),
            ],
        ),
    ] {
        for (method, (re, im)) in OPERATIONS.into_iter().zip(references) {
            let result = method(&z, 160).unwrap();
            check_decimal(&result.re, re, 65);
            check_decimal(&result.im, im, 65);
            assert!(result.re.rad.floor_log2().is_some_and(|bits| bits < -144));
            assert!(result.im.rad.floor_log2().is_some_and(|bits| bits < -144));
        }
    }
}

#[test]
fn negative_nonpoles_keep_real_axis_and_log_gamma_rejects_its_cut() {
    let z = point(-1, 0, -1);
    assert!(z.ln_gamma(128).is_err());
    for (method, (re, _)) in OPERATIONS[1..].iter().zip([
        (
            "-354490770181103205459633496668229036559509891224477425642761557971",
            "0",
        ),
        (
            "3648997397857652055902366700124443280684039533956589295287274612",
            "0",
        ),
    ]) {
        let result = method(&z, 128).unwrap();
        check_decimal(&result.re, re, 65);
        zero(&result.im);
    }
    let z = point(-2, 0, 0);
    for method in OPERATIONS {
        assert!(method(&z, 128).is_err());
    }
}

#[test]
fn input_rectangles_and_exact_positive_integer_values() {
    let center = point(-18, 3, -3);
    let rad = dyadic(1, -5);
    let value = MpComplex::new(
        Ball::new(center.re.mid, rad.clone()).unwrap(),
        Ball::new(center.im.mid, rad).unwrap(),
    );
    let outputs = OPERATIONS.map(|method| method(&value, 128).unwrap());
    for sample in [
        [
            (
                "-21716639792306600902229190393662749581242947059751396318812042303",
                "-849346410562912365235149476202776540923909044579556071723823258381",
            ),
            (
                "-48028638632028542365340061502824154259716121543214512466462717933",
                "-64577213935615569123223577846770209080029789573241194398485473223",
            ),
            (
                "170716169842153283296621848947688291902880592469798112139515100232",
                "280727898014888742265013237875232842922104418801238856169365761173",
            ),
        ],
        [
            (
                "-30036120111731826387379260671297417457680779381112690535199060533",
                "-820860101004569066615275528657907735896762104500955209136088049091",
            ),
            (
                "-25714396866242101331964461546676464887915603341894564720921628794",
                "-69447267468963054693468820578046589194576790233675004623793991499",
            ),
            (
                "150830273421634991341493923906898046882501588245113782770243665640",
                "305359693758751791845795068278159723245318841722371668629837297126",
            ),
        ],
    ] {
        for (output, (re, im)) in outputs.iter().zip(sample) {
            check_decimal(&output.re, re, 65);
            check_decimal(&output.im, im, 65);
        }
    }
    let factorial = point(5, 0, 0).gamma(8).unwrap();
    assert_eq!(
        exact(&factorial.re.mid),
        BigRational::from_integer(24.into())
    );
    assert!(factorial.re.rad.is_zero());
    zero(&factorial.im);
    let log_one = point(2, 0, 0).ln_gamma(128).unwrap();
    zero(&log_one.re);
    zero(&log_one.im);
    assert!(point(5, 0, 0).gamma(1025).is_err());
}

#[test]
fn session_declares_gamma_family_contracts_and_checks_precision() {
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
    for name in [
        "num_mp_complex_log_gamma",
        "num_mp_complex_gamma",
        "num_mp_complex_digamma",
    ] {
        let op = manifest.op(name).unwrap();
        assert_eq!(op.signature_source, SignatureSource::Declared);
        assert_eq!(op.dim_rule, Some(DimRuleInfo::Dimensionless));
        assert!(
            session
                .call_value(name, &[origin.clone(), Value::int(96)], &ctx)
                .is_err()
        );
        assert!(
            session
                .call_value(name, &[origin.clone(), Value::int(0)], &ctx)
                .is_err()
        );
    }
    let one = record(&[("real", ball(1)), ("imag", ball(0))]);
    assert_eq!(
        session
            .call_value("num_mp_complex_gamma", &[one.clone(), Value::int(1)], &ctx)
            .unwrap(),
        one
    );
}
