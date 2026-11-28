//! Closed-form and typed-ABI regressions for complex orthogonal polynomials.
use num_complex::Complex64;
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::value::ComplexValue;
use rssn::kernel::value::Value;
use rssn::kernels::special_complex_ext::orthopoly_hyper as kernel;
use std::f64::consts::FRAC_PI_2;
use std::f64::consts::PI;

fn value(z: Complex64) -> Value {
    Value::Complex(ComplexValue::Float { re: z.re, im: z.im })
}

fn assert_close(
    computed: Complex64,
    expected: Complex64,
) {
    let scale = expected.norm().max(f64::MIN_POSITIVE);
    assert!(
        (computed - expected).norm() <= 2e-12 * scale,
        "expected {expected}, got {computed}"
    );
}

#[test]
fn jacobi_handles_vanishing_recurrence_denominators() {
    let minus_one = Complex64::new(-1.0, 0.0);
    let minus_two = Complex64::new(-2.0, 0.0);
    for z in [
        Complex64::new(0.25, 0.5),
        Complex64::new(-0.25, 0.5),
        Complex64::new(-1.0, 0.0),
        Complex64::new(1.0, 0.0),
    ] {
        assert_close(
            kernel::jacobi_p_complex(2, minus_one, minus_one, z),
            (z * z - 1.0) / 4.0,
        );
        assert_close(
            kernel::jacobi_p_complex(3, minus_one, minus_two, z),
            (z + 1.0) * (z + 1.0) * (z - 1.0) / 8.0,
        );
    }
}

#[test]
fn harmonics_remain_normalized_at_high_and_negative_orders() {
    let theta = Complex64::new(FRAC_PI_2, 0.0);
    let phi = Complex64::new(0.0, 0.0);
    // Closed form from the ratio of adjacent central binomial coefficients.
    let sum: f64 = (1..=256).map(|k| (0.5 / f64::from(k)).ln_1p()).sum();
    let expected = Complex64::new((0.5 * sum).exp() / (4.0 * PI).sqrt(), 0.0);
    assert_close(
        kernel::spherical_harmonic_complex(256, 256, theta, phi),
        expected,
    );
    assert_close(
        kernel::spherical_harmonic_complex(256, -256, theta, phi),
        expected,
    );
    let phase = Complex64::from_polar((3.0 / (8.0 * PI)).sqrt(), 0.3);
    let phi = Complex64::new(0.3, 0.0);
    assert_close(kernel::spherical_harmonic_complex(1, 1, theta, phi), -phase);
    assert_close(
        kernel::spherical_harmonic_complex(1, -1, theta, phi),
        phase.conj(),
    );
}

#[test]
fn all_declared_polynomial_abi_calls_match_closed_forms() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let z = Complex64::new(0.25, -0.5);
    let zero = Complex64::new(0.0, 0.0);
    let a = Complex64::new(0.5, 0.25);
    let b = Complex64::new(-0.25, 0.125);
    let cases = [
        (
            "legendre_p",
            vec![Value::int(2), value(z)],
            (3.0 * z * z - 1.0) / 2.0,
        ),
        (
            "associated_legendre_p",
            vec![Value::int(1), Value::int(0), value(z)],
            z,
        ),
        (
            "chebyshev_t",
            vec![Value::int(2), value(z)],
            2.0 * z * z - 1.0,
        ),
        (
            "chebyshev_u",
            vec![Value::int(2), value(z)],
            4.0 * z * z - 1.0,
        ),
        (
            "hermite_h",
            vec![Value::int(2), value(z)],
            4.0 * z * z - 2.0,
        ),
        (
            "laguerre_l",
            vec![Value::int(1), value(a), value(z)],
            1.0 + a - z,
        ),
        (
            "jacobi_p",
            vec![Value::int(1), value(a), value(b), value(z)],
            (a - b + (a + b + 2.0) * z) / 2.0,
        ),
        (
            "gegenbauer_c",
            vec![Value::int(1), value(a), value(z)],
            2.0 * a * z,
        ),
        (
            "spherical_harmonic",
            vec![Value::int(0), Value::int(0), value(zero), value(zero)],
            Complex64::new(1.0 / (4.0 * PI).sqrt(), 0.0),
        ),
    ];
    for (op, args, expected) in cases {
        let name = format!("num_special_complex_{op}");
        let Value::Complex(actual) = session
            .call_value(&name, &args, &ctx)
            .unwrap_or_else(|err| panic!("{name}: {err}"))
        else {
            panic!("{name} returned a non-complex result");
        };
        let (re, im) = actual.to_f64();
        assert_close(Complex64::new(re, im), expected);
    }
    for bad_degree in [Value::int(4097), Value::int(-1), Value::Float(2.5)] {
        assert!(
            session
                .call_value(
                    "num_special_complex_chebyshev_t",
                    &[bad_degree, value(z)],
                    &ctx
                )
                .is_err()
        );
    }
}

#[test]
fn shared_bernoulli_coefficients_and_related_abi_routes() {
    use rssn::kernels::special_complex_ext::zeta_polylog;

    let session = Session::new();
    let ctx = CallCtx::new();
    let z = Complex64::new(0.25, 0.125);
    let zero = Complex64::new(0.0, 0.0);
    let reference = Complex64::new(-174_611.0 / 330.0, 0.0);
    assert_close(zeta_polylog::bernoulli_poly_complex(20, zero), reference);
    for (kind, args, expected) in [
        ("bernoulli", vec![Value::int(20), value(zero)], reference),
        ("euler", vec![Value::int(2), value(z)], z * z - z),
        (
            "clausen",
            vec![Value::int(2), value(Complex64::new(PI / 2.0, 0.0))],
            Complex64::new(0.915_965_594_177_219, 0.0),
        ),
    ] {
        let name = format!("num_special_complex_{kind}");
        let Value::Complex(result) = session.call_value(&name, &args, &ctx).unwrap() else {
            panic!("{name} returned a non-complex result");
        };
        let (re, im) = result.to_f64();
        assert_close(Complex64::new(re, im), expected);
    }
    let Value::Record(pair) = session
        .call_value(
            "num_special_complex_fresnel_auxiliary",
            &[value(zero)],
            &ctx,
        )
        .unwrap()
    else {
        panic!("expected named Fresnel components");
    };
    assert_eq!(
        pair.fields,
        [
            ("f".into(), value(Complex64::new(0.5, 0.0))),
            ("g".into(), value(Complex64::new(0.5, 0.0)))
        ]
    );
    let name = "num_special_complex_log_barnes_g";
    assert_eq!(
        session
            .call_value(name, &[value(Complex64::new(1.0, 0.0))], &ctx)
            .unwrap(),
        value(zero)
    );
    assert!(session.call_value(name, &[value(zero)], &ctx).is_err());
}
