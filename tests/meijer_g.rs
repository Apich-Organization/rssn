//! Integration tests for the Meijer G-function through the public API.
//!
//! The symbolic layer (`rules::special::meijer`) and the numeric kernel
//! (`kernels::meijer`) are exercised together: every closed form the engine
//! returns is evaluated and compared with the kernel, and the numeric phase
//! is compared with values of the special functions the G-function reduces
//! to.

use rssn::api::Config;
use rssn::api::Session;
use rssn::kernels::meijer::meijer_g_real;
use rssn::kernels::special;

/// A symbolic answer and whether it is free of heavy operators.
fn symbolic(
    session: &Session,
    text: &str,
) -> (String, bool) {
    let term = session
        .parse(text)
        .unwrap_or_else(|e| panic!("cannot parse `{text}`: {e}"));
    let answer = session
        .compute(term, &Config::new())
        .unwrap_or_else(|e| panic!("{text}: {e}"));
    (answer.term.to_string(), answer.reduced)
}

/// The numeric value of a request under bindings.
fn numeric(
    session: &Session,
    text: &str,
    bindings: &[(&str, f64)],
) -> Result<f64, rssn::api::ComputeError> {
    let term = session
        .parse(text)
        .unwrap_or_else(|e| panic!("cannot parse `{text}`: {e}"));
    let mut config = Config::new().numeric(1e-12);
    for &(name, value) in bindings {
        config = config.bind(name, value);
    }
    session
        .compute(term, &config)
        .map(|answer| answer.value.unwrap_or(f64::NAN))
}

fn close(
    got: f64,
    want: f64,
    tolerance: f64,
    what: &str,
) {
    assert!(
        (got - want).abs() <= tolerance * (1.0 + want.abs()),
        "{what}: got {got}, want {want}"
    );
}

/// `G` for a literal argument, as the kernel sees it.
fn kernel(
    m: usize,
    n: usize,
    a: &[f64],
    b: &[f64],
    z: f64,
) -> f64 {
    meijer_g_real(m, n, a, b, z).unwrap_or(f64::NAN)
}

#[test]
fn exact_reductions_are_identities() {
    // (request, parameter values, expected value)
    let cases: &[(&str, &[(&str, f64)], f64)] = &[
        (
            "meijerg(1, 0, list(), list(0), z)",
            &[("z", 0.37)],
            (-0.37_f64).exp(),
        ),
        (
            "meijerg(1, 0, list(), list(0), z)",
            &[("z", 1.9)],
            (-1.9_f64).exp(),
        ),
        (
            "meijerg(1, 0, list(), list(1/2), z)",
            &[("z", 1.44)],
            1.2 * (-1.44_f64).exp(),
        ),
        (
            "meijerg(1, 0, list(), list(2/3), z)",
            &[("z", 0.49)],
            0.49_f64.powf(2.0 / 3.0) * (-0.49_f64).exp(),
        ),
        (
            "meijerg(1, 1, list(1/2), list(1/5), z)",
            &[("z", 0.8)],
            special::gamma_numerical(0.7) * 0.8_f64.powf(0.2) * 1.8_f64.powf(-0.7),
        ),
        (
            "meijerg(1, 0, list(), list(0, 1/2), z)",
            &[("z", 0.49)],
            (2.0 * 0.7_f64).cos() / std::f64::consts::PI.sqrt(),
        ),
        (
            "meijerg(1, 0, list(), list(1/2, 0), z)",
            &[("z", 0.81)],
            (1.8_f64).sin() / std::f64::consts::PI.sqrt(),
        ),
        (
            "meijerg(1, 1, list(1), list(1/2, 0), z)",
            &[("z", 0.64)],
            std::f64::consts::PI.sqrt() * special::erf_numerical(0.8),
        ),
        (
            "meijerg(2, 0, list(1), list(1/2, 0), z)",
            &[("z", 0.64)],
            std::f64::consts::PI.sqrt() * special::erfc_numerical(0.8),
        ),
        (
            "meijerg(1, 2, list(1, 1), list(1, 0), z)",
            &[("z", 0.5)],
            1.5_f64.ln(),
        ),
        (
            "meijerg(1, 2, list(1, 1), list(1, 0), z)",
            &[("z", 2.0)],
            3.0_f64.ln(),
        ),
        (
            "meijerg(1, 2, list(1/2, 1), list(1/2, 0), z)",
            &[("z", 1.21)],
            2.0 * 1.1_f64.atan(),
        ),
        (
            "meijerg(1, 0, list(), list(3/2, 1), z)",
            &[("z", 0.25)],
            0.25 * (2.0 * 0.5_f64).sin() / std::f64::consts::PI.sqrt(),
        ),
    ];
    let session = Session::new();
    for &(request, bindings, want) in cases {
        let (text, reduced) = symbolic(&session, request);
        assert!(reduced, "`{request}` was not reduced: {text}");
        assert!(!text.contains("meijerg"), "`{request}` -> {text}");
        let term = session
            .parse(&text)
            .unwrap_or_else(|e| panic!("cannot read answer `{text}`: {e}"));
        let got = term
            .eval(bindings)
            .unwrap_or_else(|| panic!("`{text}` is not evaluable"));
        close(got, want, 1e-10, request);
    }
}

#[test]
fn bessel_and_bessel_y_reductions() {
    let session = Session::new();
    for nu in [0.0, 0.5, 1.3, 2.0] {
        let (text, reduced) = symbolic(
            &session,
            &format!("meijerg(1, 0, list(), list({nu}/2, -{nu}/2), z)"),
        );
        assert!(reduced, "J_nu reduction failed: {text}");
        // At nu = 1/2 the elementary sine form is preferred to besselj.
        assert!(text.contains("besselj") || text.contains("sin"), "{text}");
        let term = session.parse(&text).unwrap_or_else(|e| panic!("{e}"));
        for z in [0.3_f64, 1.7, 4.0] {
            let got = term.eval(&[("nu", nu), ("z", z)]).unwrap_or(f64::NAN);
            close(got, special::bessel_j(nu, 2.0 * z.sqrt()), 1e-9, "J_nu");
        }
    }
    let (text, _) = symbolic(&session, "meijerg(2, 0, list(), list(nu/2, -nu/2), z)");
    assert!(text.contains("besselk"), "{text}");
    let term = session.parse(&text).unwrap_or_else(|e| panic!("{e}"));
    close(
        term.eval(&[("nu", 0.7), ("z", 1.4)]).unwrap_or(f64::NAN),
        2.0 * special::bessel_k(0.7, 2.0 * 1.4_f64.sqrt()),
        1e-8,
        "K_nu",
    );
    // A symbolic negative argument uses the global J form; a literal negative
    // argument is in the `w > 0` domain of the modified-Bessel I form.
    let (text, reduced) = symbolic(&session, "meijerg(1, 0, list(), list(0, 0), -z)");
    assert!(reduced, "J reduction failed: {text}");
    assert!(text.contains("besselj"), "{text}");
    let term = session.parse(&text).unwrap_or_else(|e| panic!("{e}"));
    for z in [-0.3_f64, -1.6] {
        close(
            term.eval(&[("z", z)]).unwrap_or(f64::NAN),
            special::bessel_j(0.0, 2.0 * (-z).sqrt()),
            1e-9,
            "J_0",
        );
    }
    let (text, reduced) = symbolic(&session, "meijerg(1, 0, list(), list(0, 0), -2)");
    assert!(reduced, "I_0 reduction failed: {text}");
    let term = session.parse(&text).unwrap_or_else(|e| panic!("{e}"));
    close(
        term.eval(&[]).unwrap_or(f64::NAN),
        special::bessel_i(0.0, 2.0 * 2.0_f64.sqrt()),
        1e-9,
        "I_0 literal",
    );
    // General shifts for J and K.
    let (text, _) = symbolic(&session, "meijerg(1, 0, list(), list(1/3, 1/4), z)");
    assert!(text.contains("besselj"), "{text}");
    let term = session.parse(&text).unwrap_or_else(|e| panic!("{e}"));
    close(
        term.eval(&[("z", 0.7)]).unwrap_or(f64::NAN),
        0.421_144_699_126_796_1,
        1e-10,
        "J shift",
    );
    let (text, _) = symbolic(&session, "meijerg(2, 0, list(), list(1/3, 1/4), z)");
    assert!(text.contains("besselk"), "{text}");
    let term = session.parse(&text).unwrap_or_else(|e| panic!("{e}"));
    close(
        term.eval(&[("z", 0.7)]).unwrap_or(f64::NAN),
        0.309_062_261_385_866_1,
        1e-10,
        "K shift",
    );
    // The Y form needs `a1 = b2 - 1/2`; without it the request stays
    // unevaluated instead of returning `Y_2`.
    let (text, reduced) = symbolic(&session, "meijerg(2, 0, list(2), list(1, -1, 2), z)");
    assert!(!reduced, "invalid Y reduction installed: {text}");
    close(
        numeric(&session, "meijerg(2, 0, list(2), list(1, -1, 2), 0.5)", &[]).unwrap_or(f64::NAN),
        -0.210_852_477_325_523_64,
        1e-10,
        "Y invalid relation",
    );
    // Integer order: the numeric kernel takes a confluent-pole limit.
    let (text, _) = symbolic(
        &session,
        "meijerg(2, 0, list((-nu - 1)/2), list(nu/2, -nu/2, (-nu - 1)/2), z)",
    );
    assert!(text.contains("bessely"), "{text}");
    let term = session.parse(&text).unwrap_or_else(|e| panic!("{e}"));
    for nu in [0.0, 1.0, 2.0] {
        for z in [0.4_f64, 2.5] {
            close(
                term.eval(&[("nu", nu), ("z", z)]).unwrap_or(f64::NAN),
                special::bessel_y(nu, 2.0 * z.sqrt()),
                1e-7,
                "Y_nu",
            );
        }
    }
}

#[test]
fn terminating_slater_series_reduces_to_a_monomial() {
    // A = 1 + b_1 - a_i contains 0 and -1, so the 1F2 terminates at k = 0:
    // the whole G is a constant times z^{b_1}.
    let session = Session::new();
    let request = "meijerg(1, 0, list(3, 2), list(1, 5/2, 7/10), z)";
    let (text, reduced) = symbolic(&session, request);
    assert!(reduced, "{request} -> {text}");
    assert!(!text.contains("meijerg"), "{text}");
    let term = session.parse(&text).unwrap_or_else(|e| panic!("{e}"));
    for z in [0.4_f64, 1.3] {
        let got = term.eval(&[("z", z)]).unwrap_or(f64::NAN);
        let want = kernel(1, 0, &[3.0, 2.0], &[1.0, 2.5, 0.7], z);
        close(got, want, 1e-7, "terminating series");
    }
}

#[test]
fn derivative_and_antiderivative() {
    let session = Session::new();
    // d/dx e^{-x} = -e^{-x}.
    let (text, reduced) = symbolic(&session, "diff(meijerg(1, 0, list(), list(0), x), x)");
    assert!(reduced, "derivative not reduced: {text}");
    assert!(!text.contains("diff"), "{text}");
    let term = session.parse(&text).unwrap_or_else(|e| panic!("{e}"));
    close(
        term.eval(&[("x", 0.6)]).unwrap_or(f64::NAN),
        -((-0.6_f64).exp()),
        1e-10,
        "d/dx exp",
    );

    // d/dx x^{1/3} e^{-x}.
    let (text, _) = symbolic(&session, "diff(meijerg(1, 0, list(), list(1/3), x), x)");
    let term = session.parse(&text).unwrap_or_else(|e| panic!("{e}"));
    let h: f64 = 1e-6;
    let up = (1.2 + h).powf(1.0 / 3.0) * (-(1.2 + h)).exp();
    let down = (1.2 - h).powf(1.0 / 3.0) * (-(1.2 - h)).exp();
    close(
        term.eval(&[("x", 1.2)]).unwrap_or(f64::NAN),
        (up - down) / (2.0 * h),
        1e-6,
        "d/dx x^(1/3)",
    );

    // A general parameter set keeps the parameter-shift form: the derivative
    // is a shifted G divided by z, which is the DLMF identity.
    let (text, _) = symbolic(
        &session,
        "diff(meijerg(1, 0, list(5/2), list(2/5, 13/10), x), x)",
    );
    assert!(text.contains("meijerg"), "expected a shifted G: {text}");
    assert!(
        text.contains("list(0, )") || text.contains("list(0"),
        "{text}"
    );

    // Antiderivative of a G that reduces to an elementary function; check it
    // differentiates back to the integrand.
    let (primitive, reduced) = symbolic(&session, "integral(meijerg(1, 0, list(), list(0), x), x)");
    assert!(reduced, "antiderivative not reduced: {primitive}");
    let term = session.parse(&primitive).unwrap_or_else(|e| panic!("{e}"));
    let h: f64 = 1e-6;
    let up = term.eval(&[("x", 0.7 + h)]).unwrap_or(f64::NAN);
    let down = term.eval(&[("x", 0.7 - h)]).unwrap_or(f64::NAN);
    close(
        (up - down) / (2.0 * h),
        (-0.7_f64).exp(),
        1e-5,
        "integral e^{-x}",
    );

    // A linear argument: integral(e^{-2x}) = -e^{-2x}/2.
    let (primitive, _) = symbolic(&session, "integral(meijerg(1, 0, list(), list(0), 2*x), x)");
    let term = session.parse(&primitive).unwrap_or_else(|e| panic!("{e}"));
    let up = term.eval(&[("x", 0.4 + h)]).unwrap_or(f64::NAN);
    let down = term.eval(&[("x", 0.4 - h)]).unwrap_or(f64::NAN);
    close(
        (up - down) / (2.0 * h),
        (-0.8_f64).exp(),
        1e-5,
        "integral e^{-2x}",
    );
}

#[test]
fn mellin_transform_is_the_gamma_product() {
    let session = Session::new();
    let (text, reduced) = symbolic(&session, "meijer_mellin(1, 0, list(), list(0), s)");
    assert!(reduced, "Mellin not reduced: {text}");
    assert!(text.contains("gamma"), "{text}");
    let term = session.parse(&text).unwrap_or_else(|e| panic!("{e}"));
    for s in [0.4, 1.3, 2.7, 3.5] {
        close(
            term.eval(&[("s", s)]).unwrap_or(f64::NAN),
            special::gamma_numerical(s),
            1e-10,
            "Gamma(s)",
        );
    }
    // The general product: int_0^inf x^{s-1} G^{1,0}_{0,1}(x|b) dx =
    // Gamma(b+s).
    let (text, _) = symbolic(&session, "meijer_mellin(1, 0, list(), list(1/3), s)");
    let term = session.parse(&text).unwrap_or_else(|e| panic!("{e}"));
    close(
        term.eval(&[("s", 0.8)]).unwrap_or(f64::NAN),
        special::gamma_numerical(0.8 + 1.0 / 3.0),
        1e-10,
        "Gamma(b+s)",
    );
}

#[test]
fn residue_series_at_zero_and_infinity() {
    let session = Session::new();
    let (text, reduced) = symbolic(&session, "meijer_series(1, 0, list(), list(0), z, 5)");
    assert!(reduced, "series not reduced: {text}");
    let term = session.parse(&text).unwrap_or_else(|e| panic!("{e}"));
    // The truncation is good for small z and carries a big_o term.
    close(
        term.eval(&[("z", 1e-3)]).unwrap_or(f64::NAN),
        (-1e-3_f64).exp(),
        1e-8,
        "series",
    );
    assert!(text.contains("big_o"), "{text}");
    // At infinity the expansion of 1/(1+z) is 1/z - 1/z^2 + 1/z^3 - ...
    let (text, reduced) = symbolic(&session, "meijer_series(1, 1, list(0), list(0), z, 3, oo)");
    assert!(reduced, "asymptotic series not reduced: {text}");
    let term = session.parse(&text).unwrap_or_else(|e| panic!("{e}"));
    close(
        term.eval(&[("z", 50.0)]).unwrap_or(f64::NAN),
        1.0 / 51.0,
        1e-3,
        "series at oo",
    );
}

#[test]
fn numeric_phase_agrees_with_the_kernel() {
    // (m, n, a, b, z)
    let cases: &[(usize, usize, &[f64], &[f64], f64)] = &[
        (1, 0, &[], &[0.0], 0.37),
        (1, 0, &[], &[0.5], 1.9),
        (1, 0, &[], &[0.0, 0.5], 0.49),
        (1, 0, &[], &[0.5, 0.0], 0.81),
        (1, 1, &[1.0], &[0.5, 0.0], 0.64),
        (2, 0, &[1.0], &[0.5, 0.0], 0.64),
        (1, 2, &[1.0, 1.0], &[1.0, 0.0], 2.4),
        (1, 2, &[0.5, 1.0], &[0.5, 0.0], 1.21),
        (1, 0, &[], &[0.4], 5.0),
        (1, 1, &[0.3], &[0.7], 0.2),
        (2, 1, &[0.4, 1.2], &[0.8, 2.0], 1.5),
        (2, 2, &[0.9, 2.3], &[1.4, 0.2], 0.6),
        (2, 0, &[], &[0.35, -0.35], 1.4),
    ];
    let session = Session::new();
    for &(m, n, a, b, z) in cases {
        let a_list = a
            .iter()
            .map(|v| v.to_string())
            .collect::<Vec<_>>()
            .join(", ");
        let b_list = b
            .iter()
            .map(|v| v.to_string())
            .collect::<Vec<_>>()
            .join(", ");
        let request = format!("meijerg({m}, {n}, list({a_list}), list({b_list}), {z})");
        let got = numeric(&session, &request, &[]).unwrap_or_else(|e| panic!("{request}: {e}"));
        let want = kernel(m, n, a, b, z);
        close(got, want, 1e-9, &request);
    }
    // Bound parameter symbols work too.
    let got = numeric(
        &session,
        "meijerg(1, 1, list(a), list(b), z)",
        &[("a", 0.5), ("b", 0.7), ("z", 1.2)],
    )
    .unwrap_or_else(|e| panic!("{e}"));
    close(
        got,
        kernel(1, 1, &[0.5], &[0.7], 1.2),
        1e-9,
        "bound parameters",
    );
}

#[test]
fn numeric_values_match_the_special_functions() {
    // A moderately large parameter table: the numeric phase must reproduce
    // the defining special function at several points.
    let session = Session::new();
    for z in [0.2_f64, 0.8, 2.5] {
        close(
            numeric(
                &session,
                &format!("meijerg(1, 0, list(), list(0), {z})"),
                &[],
            )
            .unwrap_or(f64::NAN),
            (-z).exp(),
            1e-10,
            "exp",
        );
        close(
            numeric(
                &session,
                &format!("meijerg(1, 2, list(1, 1), list(1, 0), {z})"),
                &[],
            )
            .unwrap_or(f64::NAN),
            (1.0 + z).ln(),
            1e-8,
            "log",
        );
        close(
            numeric(
                &session,
                &format!("meijerg(2, 0, list(1), list(1/2, 0), {z})"),
                &[],
            )
            .unwrap_or(f64::NAN),
            std::f64::consts::PI.sqrt() * special::erfc_numerical(z.sqrt()),
            1e-9,
            "erfc",
        );
        close(
            numeric(
                &session,
                &format!("meijerg(1, 0, list(), list(1/2, 0), {z})"),
                &[],
            )
            .unwrap_or(f64::NAN),
            (2.0 * z.sqrt()).sin() / std::f64::consts::PI.sqrt(),
            1e-10,
            "sin",
        );
    }
    for (nu, z) in [(0.0_f64, 0.5_f64), (1.0, 1.5), (2.0, 3.0)] {
        close(
            numeric(
                &session,
                &format!("meijerg(2, 0, list(), list({nu}/2, -{nu}/2), {z})"),
                &[],
            )
            .unwrap_or(f64::NAN),
            2.0 * special::bessel_k(nu, 2.0 * z.sqrt()),
            1e-7,
            "K",
        );
        close(
            numeric(
                &session,
                &format!(
                    "meijerg(2, 0, list((-{nu} - 1)/2), list({nu}/2, -{nu}/2, (-{nu} - 1)/2), {z})"
                ),
                &[],
            )
            .unwrap_or(f64::NAN),
            special::bessel_y(nu, 2.0 * z.sqrt()),
            1e-6,
            "Y",
        );
    }
}

#[test]
fn invalid_arguments_stay_unevaluated() {
    let session = Session::new();
    // m > q.
    let (text, reduced) = symbolic(&session, "meijerg(3, 0, list(0), list(0), z)");
    assert!(!reduced, "invalid order was reduced to {text}");
    assert!(text.contains("meijerg"), "{text}");
    // n > p.
    let (text, reduced) = symbolic(&session, "meijerg(0, 2, list(0), list(0), z)");
    assert!(!reduced, "invalid order was reduced to {text}");
    assert!(text.contains("meijerg"), "{text}");
    // A numeric target on an invalid request fails cleanly.
    assert!(
        numeric(&session, "meijerg(3, 0, list(0), list(0), 1)", &[]).is_err(),
        "invalid order must not evaluate"
    );
    // A numeric target without a binding for the argument fails cleanly.
    assert!(numeric(&session, "meijerg(1, 0, list(), list(0), x)", &[]).is_err());
    // And with the binding it succeeds.
    close(
        numeric(&session, "meijerg(1, 0, list(), list(0), x)", &[("x", 1.0)]).unwrap_or(f64::NAN),
        (-1.0_f64).exp(),
        1e-10,
        "bound exp",
    );
}

#[test]
fn m_zero_and_n_zero_are_zero_in_their_contours() {
    let session = Session::new();
    // m = 0 with p < q: the contour encloses no pole.
    let (_, reduced) = symbolic(&session, "meijerg(0, 1, list(1/2), list(1/3, 2), z)");
    let value =
        numeric(&session, "meijerg(0, 1, list(1/2), list(1/3, 2), 0.7)", &[]).unwrap_or(f64::NAN);
    assert!(value == 0.0, "G^{{0,1}} must vanish, got {value}");
    // n = 0 with p > q.
    let value =
        numeric(&session, "meijerg(1, 0, list(1/2, 3), list(1/3), 2.0)", &[]).unwrap_or(f64::NAN);
    assert!(value == 0.0, "G^{{1,0}} must vanish, got {value}");
    let _ = reduced;
    // z = 0: e^{-0} = 1.
    close(
        numeric(&session, "meijerg(1, 0, list(), list(0), 0)", &[]).unwrap_or(f64::NAN),
        1.0,
        1e-12,
        "at zero",
    );
}
