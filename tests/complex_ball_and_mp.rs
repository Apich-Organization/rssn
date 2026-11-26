//! Phase 1 integration tests: ComplexBall, Cball elementary enclosures,
//! Krawczyk root certification, certified Rouché, MpComplex, and Session operators.

use num_complex::Complex64;
use rssn::api::Config;
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::value::TensorValue;
use rssn::kernel::value::Value;
use rssn::kernels::ball::Cball;
use rssn::kernels::complex_ball::ComplexBall;
use rssn::kernels::complex_ball::ComplexBallJet;
use rssn::kernels::complex_ball::certify_root_1d;
use rssn::kernels::complex_ball::certify_rouche_circle;
use rssn::kernels::complex_ball::solve_linear_ball;
use rssn::kernels::mp;

fn sample_disk(b: ComplexBall) -> Vec<Complex64> {
    let mut pts = vec![b.mid];
    for r_frac in [0.25, 0.5, 0.99] {
        for k in 0..8 {
            let theta = std::f64::consts::TAU * (k as f64) / 8.0;
            let (s, c) = theta.sin_cos();
            pts.push(Complex64::new(
                b.mid.re + b.rad * r_frac * c,
                b.mid.im + b.rad * r_frac * s,
            ));
        }
    }
    pts
}

fn sample_rect(b: Cball) -> Vec<Complex64> {
    let rs = [b.re.lo, f64::midpoint(b.re.lo, b.re.hi), b.re.hi];
    let is = [b.im.lo, f64::midpoint(b.im.lo, b.im.hi), b.im.hi];
    let mut pts = Vec::new();
    for r in rs {
        for i in is {
            pts.push(Complex64::new(r, i));
        }
    }
    pts
}

#[test]
fn complex_ball_encloses_all_elementary_functions() {
    let b = ComplexBall::new(Complex64::new(0.4, 0.5), 0.05).unwrap();
    let pts = sample_disk(b);

    let b_exp = b.exp();
    let b_ln = b.ln().unwrap();
    let b_sqrt = b.sqrt().unwrap();
    let b_sin = b.sin();
    let b_cos = b.cos();
    let b_tan = b.tan().unwrap();
    let b_sinh = b.sinh();
    let b_cosh = b.cosh();
    let b_tanh = b.tanh().unwrap();
    let b_asin = b.asin().unwrap();
    let b_acos = b.acos().unwrap();
    let b_atan = b.atan().unwrap();
    let b_asinh = b.asinh().unwrap();
    let b_atanh = b.atanh().unwrap();

    for z in pts {
        assert!(b_exp.contains(z.exp()), "exp failed at {z}");
        assert!(b_ln.contains(z.ln()), "ln failed at {z}");
        assert!(b_sqrt.contains(z.sqrt()), "sqrt failed at {z}");
        assert!(b_sin.contains(z.sin()), "sin failed at {z}");
        assert!(b_cos.contains(z.cos()), "cos failed at {z}");
        assert!(b_tan.contains(z.tan()), "tan failed at {z}");
        assert!(b_sinh.contains(z.sinh()), "sinh failed at {z}");
        assert!(b_cosh.contains(z.cosh()), "cosh failed at {z}");
        assert!(b_tanh.contains(z.tanh()), "tanh failed at {z}");
        assert!(b_asin.contains(z.asin()), "asin failed at {z}");
        assert!(b_acos.contains(z.acos()), "acos failed at {z}");
        assert!(b_atan.contains(z.atan()), "atan failed at {z}");
        assert!(b_asinh.contains(z.asinh()), "asinh failed at {z}");
        assert!(b_atanh.contains(z.atanh()), "atanh failed at {z}");
    }

    // acosh away from (-inf, 1]
    let b2 = ComplexBall::new(Complex64::new(2.0, 0.5), 0.05).unwrap();
    let b_acosh = b2.acosh().unwrap();
    for z in sample_disk(b2) {
        assert!(b_acosh.contains(z.acosh()), "acosh failed at {z}");
    }
}

#[test]
fn complex_ball_branch_cut_rejection() {
    let on_neg_axis = ComplexBall::new(Complex64::new(-1.0, 0.01), 0.05).unwrap();
    assert!(on_neg_axis.ln().is_none());
    assert!(on_neg_axis.sqrt().is_none());
    assert!(on_neg_axis.acosh().is_none());

    let near_origin = ComplexBall::new(Complex64::new(0.01, 0.0), 0.02).unwrap();
    assert!(near_origin.ln().is_none());
    assert!(near_origin.recip().is_none());
    let s = near_origin.sqrt().unwrap();
    assert!(s.contains(Complex64::new(0.01, 0.0).sqrt()));
}

#[test]
fn cball_elementary_enclosures_hold() {
    let rect = Cball::rectangle((0.5, 0.6), (0.3, 0.4));
    let pts = sample_rect(rect);

    let r_exp = rect.exp();
    let r_ln = rect.ln().unwrap();
    let r_sqrt = rect.sqrt().unwrap();
    let r_sin = rect.sin();
    let r_cos = rect.cos();
    let r_tan = rect.tan().unwrap();
    let r_sinh = rect.sinh();
    let r_cosh = rect.cosh();
    let r_tanh = rect.tanh().unwrap();
    let w = Cball::rectangle((1.2, 1.3), (-0.1, 0.1));
    let r_pow = rect.powc(&w).unwrap();

    for z in pts {
        assert!(r_exp.contains(z.exp()));
        assert!(r_ln.contains(z.ln()));
        assert!(r_sqrt.contains(z.sqrt()));
        assert!(r_sin.contains(z.sin()));
        assert!(r_cos.contains(z.cos()));
        assert!(r_tan.contains(z.tan()));
        assert!(r_sinh.contains(z.sinh()));
        assert!(r_cosh.contains(z.cosh()));
        assert!(r_tanh.contains(z.tanh()));
        assert!(r_pow.contains(z.powc(Complex64::new(1.25, 0.0))));
    }
}

#[test]
fn wrapping_effect_advantage_over_cball() {
    let angle = std::f64::consts::TAU / 100.0;
    let rot = Complex64::from_polar(1.0, angle);
    let mut ball = ComplexBall::new(Complex64::new(1.0, 0.0), 1e-12).unwrap();
    let mut rect = Cball::square(1.0, 0.0, 1e-12);
    let rot_ball = ComplexBall::of(rot);
    let rot_rect = Cball::of(rot);
    for _ in 0..100 {
        ball = ball * rot_ball;
        rect = rect * rot_rect;
    }
    assert!(ball.rad < rect.half_width() * 0.01);
    assert!(ball.contains(Complex64::new(1.0, 0.0)));
}

#[test]
fn krawczyk_certifies_transcendental_and_rejects_empty_box() {
    // exp(z) - z - 2 = 0 has a root near -1.8414056604369606 + 0i (and complex roots)
    let root_approx = Complex64::new(-1.841_405_660_437, 0.0);
    let ball = ComplexBall::new(root_approx, 1e-3).unwrap();
    let cert = certify_root_1d(
        |z| {
            let ez = z.exp();
            let two = ComplexBallJet::constant(ComplexBall::point(2.0, 0.0));
            Some(ez.sub(&z).sub(&two))
        },
        ball,
    );
    assert!(cert.is_some());

    // sin(z) - 2 = 0 has a root at pi/2 + i*acosh(2)
    let sin_root = Complex64::new(std::f64::consts::FRAC_PI_2, 2.0_f64.acosh());
    let ball_sin = ComplexBall::new(sin_root + Complex64::new(1e-4, -1e-4), 1e-2).unwrap();
    let cert_sin = certify_root_1d(
        |z| {
            let two = ComplexBallJet::constant(ComplexBall::point(2.0, 0.0));
            Some(z.sin().sub(&two))
        },
        ball_sin,
    );
    assert!(cert_sin.is_some());
    assert!(cert_sin.unwrap().contains(sin_root));

    // Empty box far from any root
    let empty = ComplexBall::new(Complex64::new(0.0, 0.0), 1e-2).unwrap();
    let cert_empty = certify_root_1d(
        |z| {
            let two = ComplexBallJet::constant(ComplexBall::point(2.0, 0.0));
            Some(z.sin().sub(&two))
        },
        empty,
    );
    assert!(cert_empty.is_none());
}

#[test]
fn rouche_circle_certification() {
    // z^5 + 3z + 1 on |z|=1: |3z| = 3 > |z^5 + 1| <= 2 -> 1 root
    assert!(certify_rouche_circle(
        |z| Some(z.scale_real(3.0)),
        |z| Some(z.powi(5)? + ComplexBall::one()),
        Complex64::new(0.0, 0.0),
        1.0,
        64,
    ));
    // z^5 + 3z + 1 on |z|=2: |z^5| = 32 > |3z + 1| <= 7 -> 5 roots
    assert!(certify_rouche_circle(
        |z| z.powi(5),
        |z| Some(z.scale_real(3.0) + ComplexBall::one()),
        Complex64::new(0.0, 0.0),
        2.0,
        64,
    ));
    // exp(z) - 3z^2 on |z|=1: |3z^2| = 3 > |exp(z)| <= e < 3 -> 2 roots
    assert!(certify_rouche_circle(
        |z| Some(z.powi(2)?.scale_real(3.0)),
        |z| Some(z.exp()),
        Complex64::new(0.0, 0.0),
        1.0,
        128,
    ));
}

#[test]
fn complex_ball_linear_system_enclosure() {
    let a = vec![
        vec![
            ComplexBall::new(Complex64::new(2.0, 1.0), 1e-10).unwrap(),
            ComplexBall::new(Complex64::new(0.0, -1.0), 1e-10).unwrap(),
        ],
        vec![
            ComplexBall::new(Complex64::new(1.0, 0.0), 1e-10).unwrap(),
            ComplexBall::new(Complex64::new(3.0, 1.0), 1e-10).unwrap(),
        ],
    ];
    let x_true = [Complex64::new(1.0, 2.0), Complex64::new(-1.0, 0.5)];
    let b = vec![
        ComplexBall::of(a[0][0].mid * x_true[0] + a[0][1].mid * x_true[1]),
        ComplexBall::of(a[1][0].mid * x_true[0] + a[1][1].mid * x_true[1]),
    ];
    let sol = solve_linear_ball(&a, &b).unwrap();
    assert!(sol[0].contains(x_true[0]));
    assert!(sol[1].contains(x_true[1]));
    assert!(sol[0].rad < 1e-8);
    assert!(sol[1].rad < 1e-8);
}

#[test]
fn mp_complex_50_digit_identities() {
    let (re, im) = mp::eval_complex_text("asin(sin(0.3 + 0.4*I))", 50).unwrap();
    assert_eq!(re, "0.3");
    assert_eq!(im, "0.4");

    let (re2, im2) = mp::eval_complex_text("atan(tan(0.2 - 0.3*I))", 50).unwrap();
    assert_eq!(re2, "0.2");
    assert_eq!(im2, "-0.3");
}

fn float_vec(xs: Vec<f64>) -> Value {
    Value::Tensor(TensorValue::vector(
        xs.into_iter().map(Value::Float).collect(),
    ))
}

#[test]
fn session_phase1_operators_and_diagnostics() {
    let session = Session::new();
    let ctx = CallCtx::new();

    // mp_complex_eval & mp_complex_ball (return list of decimal strings -> Value::Term or Tensor)
    let val = session
        .call_value(
            "mp_complex_eval",
            &[Value::String("exp(I*pi/3)".to_string()), Value::int(25)],
            &ctx,
        )
        .unwrap();
    assert!(matches!(val, Value::Term(_) | Value::Tensor(_)));

    let ball_val = session
        .call_value(
            "mp_complex_ball",
            &[Value::String("1 + 2*I".to_string()), Value::int(20)],
            &ctx,
        )
        .unwrap();
    assert!(matches!(ball_val, Value::Term(_) | Value::Tensor(_)));

    // cball_mul & cball_exp
    let b1 = float_vec(vec![1.0, 2.0, 0.01]);
    let b2 = float_vec(vec![3.0, -1.0, 0.01]);
    let prod = session
        .call_value("cball_mul", &[b1.clone(), b2], &ctx)
        .unwrap();
    assert!(matches!(prod, Value::Tensor(_)));
    let exp_b = session.call_value("cball_exp", &[b1], &ctx).unwrap();
    assert!(matches!(exp_b, Value::Tensor(_)));

    // cinterval_add via term API
    let t = session
        .parse("cinterval_add(list(0.9, 1.1, -0.1, 0.1), list(1.9, 2.1, 0.9, 1.1))")
        .unwrap();
    let sum = session.compute(t, &Config::default()).unwrap();
    assert!(sum.reduced);

    // Diagnostics on invalid domain: cball_ln on (-1, 0, 0.1)
    let bad_ball = float_vec(vec![-1.0, 0.0, 0.1]);
    let err = session.apply_diag("cball_ln", &[bad_ball], &ctx);
    assert!(err.is_err());
}
