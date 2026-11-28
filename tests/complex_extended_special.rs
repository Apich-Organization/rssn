use num_complex::Complex64 as C;
use rssn::api::Session;
use rssn::kernel::call::CallCtx;
use rssn::kernel::value::ComplexValue;
use rssn::kernel::value::Value;

fn complex(z: C) -> Value {
    Value::Complex(ComplexValue::Float { re: z.re, im: z.im })
}

fn rssn_complex(value: &Value) -> C {
    let Value::Complex(ComplexValue::Float { re, im }) = value else {
        panic!("expected complex value, got {value:?}");
    };
    C::new(*re, *im)
}

fn near(
    actual: C,
    expected: C,
    tolerance: f64,
) {
    let relative = (actual - expected).norm() / expected.norm().max(f64::MIN_POSITIVE);
    assert!(
        relative <= tolerance,
        "relative error {relative}: {actual} versus {expected}"
    );
}

#[test]
fn complex_bessel_orders_against_high_precision_references() {
    // Real orders passed through the complex-order adapter use the same cut.
    let order: f64 = 0.7;
    let point = C::new(20.0, 0.0);
    let j = rssn::kernels::special_complex::bessel_j_complex(order, point);
    let phase = C::from_polar(1.0, order * std::f64::consts::PI);
    for (side, phase) in [(0.0, phase), (-0.0, phase.conj())] {
        near(
            rssn::kernels::special_complex::bessel_j_complex(order, C::new(-20.0, side)),
            phase * j,
            1e-13,
        );
    }
    let session = Session::new();
    let ctx = CallCtx::new();
    for (order, z, j, i) in [
        (
            C::new(0.75, 0.3),
            C::new(1.2, 0.8),
            C::new(0.6356840418226485, 0.0023860259500369865),
            C::new(0.7610651360302318, 0.4153095121776765),
        ),
        (
            C::new(-4.0, 1e-20),
            C::new(1.2, 0.4),
            C::new(0.0020397442508668286, 0.005910380262010307),
            C::new(0.0016636022704707946, 0.006908728825839641),
        ),
        (
            C::new(65.0, 0.3),
            C::new(0.5, 0.75),
            C::new(2.5933624298821363e-114, -1.242569404980434e-114),
            C::new(2.5943118255943474e-114, -1.2249151765260754e-114),
        ),
        (
            C::new(0.25, -0.8),
            C::new(-2.25, 1e-10),
            C::new(10.101731267573044, -0.3685019267239774),
            C::new(22.949149896537097, 31.892273590667358),
        ),
        (
            C::new(2.0, -0.2),
            C::new(15.0, 3.0),
            C::new(0.7718374653711628, 2.5613426258120757),
            C::new(-291310.5323139816, 55576.85716285055),
        ),
        (
            C::new(2.0, -0.2),
            C::new(21.0, 3.0),
            C::new(-0.11609924278751739, 2.307376885413188),
            C::new(-103138000.67277376, 18820727.269086324),
        ),
    ] {
        for (operator, expected) in [("bessel_j_order", j), ("bessel_i_order", i)] {
            let actual = session
                .call_value(
                    &format!("num_special_complex_{operator}"),
                    &[complex(order), complex(z)],
                    &ctx,
                )
                .unwrap();
            near(rssn_complex(&actual), expected, 5e-10);
        }
    }
}

#[test]
fn legendre_q_minimal_solution_and_exponential_integral_against_references() {
    let session = Session::new();
    let ctx = CallCtx::new();
    for (order, z, expected) in [
        (0, C::new(2.5, 0.0), C::new(0.4236489301936018, 0.0)),
        (1, C::new(1e+170, 1e+169), C::new(0.0, 0.0)),
        (
            8,
            C::new(2.1, -0.4),
            C::new(-7.203678299514355e-07, 1.9757426679857226e-06),
        ),
        (
            50,
            C::new(1.15, 0.35),
            C::new(-5.1487215664372866e-18, 1.60291491136409e-17),
        ),
        (85, C::new(1.01, 0.0), C::new(2.0184051911820726e-06, 0.0)),
        (
            52,
            C::new(-1.03, -0.12),
            C::new(2.268873740382543e-10, -1.1200033647077126e-10),
        ),
        (
            112,
            C::new(0.87, 0.008),
            C::new(-0.017059561053198984, -0.021091452954572727),
        ),
        (
            6,
            C::new(0.2, -0.013),
            C::new(-0.43984625568025737, -0.11553464312458438),
        ),
        (
            19,
            C::new(4.01, 3.0),
            C::new(3.6249983790363264e-21, -1.874514327554399e-21),
        ),
    ] {
        let actual = session
            .call_value(
                "num_special_complex_legendre_q",
                &[Value::int(order), complex(z)],
                &ctx,
            )
            .unwrap();
        near(rssn_complex(&actual), expected, 2e-11);
    }
    for (order, z, expected) in [
        (
            0,
            C::new(0.4, 0.5),
            C::new(0.18200048478208858, -1.030921978747187),
        ),
        (
            2,
            C::new(-3.0, 0.25),
            C::new(-10.292739079644885, -6.952919202954863),
        ),
        (
            4,
            C::new(-11.5, -0.25),
            C::new(-14412.32501011523, -2284.4250260478975),
        ),
        (
            7,
            C::new(4.25, 5.0),
            C::new(0.0008022015292875177, 0.000881738418016608),
        ),
        (
            31,
            C::new(0.01, 0.002),
            C::new(0.032990214612926845, -6.825487357100607e-05),
        ),
        (
            32,
            C::new(0.01, 0.002),
            C::new(0.031926381128847874, -6.598052348136447e-05),
        ),
        (
            64,
            C::new(14.0, 0.6),
            C::new(8.84334849275412e-09, -6.152606712988533e-09),
        ),
        (
            128,
            C::new(-0.1, 0.25),
            C::new(0.008433989997176204, -0.0021713998058296327),
        ),
        (
            4096,
            C::new(1.254, 0.125),
            C::new(6.912009784944528e-05, -8.687438416978307e-06),
        ),
    ] {
        let actual = session
            .call_value(
                "num_special_complex_exp_integral_en",
                &[Value::int(order), complex(z)],
                &ctx,
            )
            .unwrap();
        near(rssn_complex(&actual), expected, 5e-11);
    }
    let huge = C::new(0.0, 1e308);
    let exponential = (-huge).exp();
    for order in [0, 1, 4096] {
        let v = session
            .call_value(
                "num_special_complex_exp_integral_en",
                &[Value::int(order), complex(huge)],
                &ctx,
            )
            .unwrap();
        // Compare after rescaling since 1/z is subnormal here.
        near(
            rssn_complex(&v) * 1e308,
            exponential * C::new(0.0, -1.0),
            1e-14,
        );
    }
    for side in [0.0_f64, -0.0] {
        let e1 = session
            .call_value(
                "num_special_complex_exp_integral_en",
                &[Value::int(1), complex(C::new(-32.0, side))],
                &ctx,
            )
            .unwrap();
        let q0 = session
            .call_value(
                "num_special_complex_legendre_q",
                &[Value::int(0), complex(C::new(0.2, side))],
                &ctx,
            )
            .unwrap();
        assert!((rssn_complex(&e1).im + std::f64::consts::PI.copysign(side)).abs() < 1e-14);
        assert!((rssn_complex(&q0).im + std::f64::consts::FRAC_PI_2.copysign(side)).abs() < 1e-14);
    }
    assert!(
        session
            .call_value(
                "num_special_complex_legendre_q",
                &[Value::int(4097), complex(C::new(2.0, 0.0))],
                &ctx
            )
            .is_err()
    );
    assert!(
        session
            .call_value(
                "num_special_complex_legendre_q",
                &[Value::int(1), complex(C::new(1.0, 0.0))],
                &ctx
            )
            .is_err()
    );
}

#[test]
fn complex_anger_weber_and_kelvin_records_match_references() {
    let session = Session::new();
    let ctx = CallCtx::new();
    for (order, z, anger, weber, kelvin) in [
        (
            0.6,
            C::new(1.0, 0.5),
            C::new(0.853950368272301, 0.038591795070701485),
            C::new(0.21476741029097676, -0.2919447593850454),
            [
                C::new(0.0648568328986426, -0.10954035136499028),
                C::new(0.7633357592497063, 0.22891768937506454),
                C::new(-0.3551084864807993, 0.014740903895644381),
                C::new(-0.4081290972509181, 0.3995917313621866),
            ],
        ),
        (
            0.75,
            C::new(2.1, 0.0),
            C::new(0.6211215727693332, 0.0),
            C::new(-0.22247174785862772, 0.0),
            [
                C::new(-0.8762323671111644, 0.0),
                C::new(0.8319394937020183, 0.0),
                C::new(-0.20572840569322026, 0.0),
                C::new(-0.009541708642578434, 0.0),
            ],
        ),
        (
            -0.75,
            C::new(0.88, -0.45),
            C::new(-0.19352643947702056, 0.23904838795077207),
            C::new(-0.7833619894084602, -0.08219264463237712),
            [
                C::new(0.288282072886818, -0.2891999483362099),
                C::new(-0.5031861319180683, -0.15361621225858998),
                C::new(0.6516647991631219, 0.4125358035976994),
                C::new(-0.062145903201084644, 0.2279428941094235),
            ],
        ),
        (
            0.0,
            C::new(8.5, 1.3),
            C::new(0.047364224939895345, -0.4623580643169684),
            C::new(-0.6067175054232399, -0.01974451069537564),
            [
                C::new(14.687463728401518, 72.1012345965978),
                C::new(-68.4302812260195, -7.310778673258873),
                C::new(0.0007174639548815963, -0.0013233344045417708),
                C::new(-0.001003904022163363, -0.0004550543000787671),
            ],
        ),
        (
            2.0,
            C::new(-1.7, 0.35),
            C::new(0.2839973117720119, -0.08805735876145432),
            C::new(-0.19520496158179354, -0.027678643912068555),
            [
                C::new(0.06505448737764234, -0.06816495066010944),
                C::new(-0.3428586057108238, 0.14041045795362633),
                C::new(0.09426096647732875, -0.14748581865292876),
                C::new(0.8667811746868056, 1.3587818858023435),
            ],
        ),
        (
            3.2,
            C::new(11.0, -0.15),
            C::new(0.18616059531480908, -0.019829507718699405),
            C::new(0.12611832176406573, 0.028982701876684825),
            [
                C::new(204.7667642764725, -17.38104147256608),
                C::new(36.417114200448054, -24.615644208844436),
                C::new(0.00012803473672019845, 3.3140047665664585e-05),
                C::new(-0.000176285057004927, -7.319923709713618e-06),
            ],
        ),
    ] {
        let args = [Value::Float(order), complex(z)];
        for (operator, expected) in [("anger_j", anger), ("weber_e", weber)] {
            let actual = session
                .call_value(&format!("num_special_complex_{operator}"), &args, &ctx)
                .unwrap();
            near(rssn_complex(&actual), expected, 5e-11);
        }
        let Value::Record(record) = session
            .call_value("num_special_complex_kelvin", &args, &ctx)
            .unwrap()
        else {
            panic!("expected named Kelvin components");
        };
        for (label, expected) in ["ber", "bei", "ker", "kei"].into_iter().zip(kelvin) {
            let actual = record
                .fields
                .iter()
                .find(|(name, _)| name == label)
                .unwrap();
            near(rssn_complex(&actual.1), expected, 1e-7);
        }
    }
}
