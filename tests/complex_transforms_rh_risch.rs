use num_complex::Complex;
use rssn::api::Config;
use rssn::api::Session;
use rssn::kernels::complex_transforms::borel_transform_series;
use rssn::kernels::complex_transforms::bromwich_talbot;
use rssn::kernels::complex_transforms::bromwich_weideman_parabola;
use rssn::kernels::complex_transforms::gaver_stehfest_inverse_laplace;
use rssn::kernels::complex_transforms::hankel_transform;
use rssn::kernels::complex_transforms::hilbert_transform_real;
use rssn::kernels::complex_transforms::inverse_mellin_transform;
use rssn::kernels::complex_transforms::meijer_g;
use rssn::kernels::complex_transforms::mellin_barnes_integral;
use rssn::kernels::complex_transforms::mellin_transform;
use rssn::kernels::riemann_hilbert::Matrix2;
use rssn::kernels::riemann_hilbert::cauchy_transform_arc;
use rssn::kernels::riemann_hilbert::plemelj_boundary_values;
use rssn::kernels::riemann_hilbert::scalar_rh_factorize_arc;
use rssn::kernels::riemann_hilbert::solve_matrix_rh_chebyshev;
use rssn::kernels::special_complex::gamma_complex;

fn eval_sym(src: &str) -> String {
    let session = Session::new();
    let term = session.parse(src).unwrap_or_else(|e| panic!("{src}: {e}"));
    let ans = session
        .compute(term, &Config::default())
        .unwrap_or_else(|e| panic!("{src}: {e}"));
    ans.term.to_string()
}

#[test]
fn mellin_bromwich_meijer_g_hilbert_hankel_borel() {
    // Mellin[exp(-x)](s) = Gamma(s)
    let s = Complex::new(1.5, 0.4);
    let m_exp = mellin_transform(|x| Complex::new((-x).exp(), 0.0), s, 12.0, 1e-10);
    let g_s = gamma_complex(s);
    assert!((m_exp - g_s).norm() < 1e-8, "m_exp={m_exp} vs g_s={g_s}");

    // Inverse Mellin of Gamma(s) at x = 0.7 -> exp(-0.7)
    let inv_m = inverse_mellin_transform(gamma_complex, Complex::new(0.7, 0.0), 1.0, 25.0, 1e-9);
    assert!(
        (inv_m - Complex::new((-0.7_f64).exp(), 0.0)).norm() < 1e-7,
        "inv_m={inv_m}"
    );

    // Bromwich Laplace inversions for F(s) = 1 / (s^2 + 1) -> sin(t)
    let f_lap = |z: Complex<f64>| Complex::new(1.0, 0.0) / (z * z + Complex::new(1.0, 0.0));
    let t = 1.7;
    let tal = bromwich_talbot(f_lap, t, 32).unwrap();
    let par = bromwich_weideman_parabola(f_lap, t, 32).unwrap();
    let gs = gaver_stehfest_inverse_laplace(|r| 1.0 / (r + 1.5), t, 12).unwrap();
    assert!(
        (tal.re - t.sin()).abs() < 1e-9 && tal.im.abs() < 1e-9,
        "tal={tal}, want={}",
        t.sin()
    );
    assert!((par.re - t.sin()).abs() < 1e-8 && par.im.abs() < 1e-8);
    assert!((gs - (-1.5 * t).exp()).abs() < 1e-4, "gs={gs}");

    // Meijer G_{0,1}^{1,0}(z | 0) = exp(-z)
    let z = Complex::new(1.3, 0.2);
    let mg_exp = meijer_g(&[], &[Complex::new(0.0, 0.0)], 1, 0, z, 0.5, 25.0, 1e-9).unwrap();
    assert!((mg_exp - (-z).exp()).norm() < 1e-7, "mg_exp={mg_exp}");

    // Meijer G_{1,1}^{1,1}(z | 0; 0) = 1 / (1 + z)
    let mg_rat = meijer_g(
        &[Complex::new(0.0, 0.0)],
        &[Complex::new(0.0, 0.0)],
        1,
        1,
        z,
        0.5,
        25.0,
        1e-9,
    )
    .unwrap();
    let exact_rat = Complex::new(1.0, 0.0) / (Complex::new(1.0, 0.0) + z);
    assert!((mg_rat - exact_rat).norm() < 1e-7, "mg_rat={mg_rat}");

    // Mellin-Barnes directly for exp(-z)
    let mb = mellin_barnes_integral(&[Complex::new(0.0, 0.0)], &[], &[], &[], z, 0.5, 25.0, 1e-9);
    assert!((mb - (-z).exp()).norm() < 1e-7);

    // Hilbert transform of 1/(1 + t^2) is x/(1 + x^2)
    let x = 0.6;
    let ht = hilbert_transform_real(|u| 1.0 / (1.0 + u * u), x, 200.0, 1e-7).unwrap();
    assert!((ht - x / (1.0 + x * x)).abs() < 1e-4, "ht={ht}");

    // Hankel transform order 0 of exp(-r^2 / 2) is exp(-k^2 / 2)
    let k = 1.2;
    let hk = hankel_transform(
        |r| Complex::new((-0.5 * r * r).exp(), 0.0),
        0.0,
        k,
        12.0,
        1e-9,
    );
    assert!((hk.re - (-0.5 * k * k).exp()).abs() < 1e-7, "hk={hk}");

    // Borel transform of sum n! t^n is geometric series sum t^n = 1/(1 - t)
    let facts: Vec<Complex<f64>> = (0..30)
        .map(|n| {
            let mut p = 1.0_f64;
            for j in 2..=n {
                p *= j as f64;
            }
            Complex::new(p, 0.0)
        })
        .collect();
    let bt = borel_transform_series(&facts, Complex::new(0.3, 0.2));
    let exact_bt = Complex::new(1.0, 0.0) / (Complex::new(1.0, 0.0) - Complex::new(0.3, 0.2));
    assert!((bt - exact_bt).norm() < 1e-12);
}

#[test]
fn riemann_hilbert_plemelj_scalar_and_matrix_chebyshev() {
    // Plemelj jump C_+ - C_- = phi
    let phi = |t: f64| Complex::new(t * t + 1.0, (2.0 * t).sin());
    let x = 0.35;
    let (cp, cm) = plemelj_boundary_values(&phi, x, 1e-10).unwrap();
    assert!((cp - cm - phi(x)).norm() < 1e-9);

    // Off-arc Cauchy transform agrees with Plemelj limit as eps -> 0+
    let c_above = cauchy_transform_arc(&phi, Complex::new(x, 1e-5), 1e-9);
    let c_below = cauchy_transform_arc(&phi, Complex::new(x, -1e-5), 1e-9);
    assert!((c_above - cp).norm() < 1e-3);
    assert!((c_below - cm).norm() < 1e-3);

    // Scalar RH: X^+(t) = G(t) X^-(t) with G(t) = exp((1 - t^2) * (0.4 + 0.3i))
    let g_scalar = |t: f64| ((1.0 - t * t) * Complex::new(0.4, 0.3)).exp();
    let xp = scalar_rh_factorize_arc(&g_scalar, Complex::new(x, 1e-6), 1e-9);
    let xm = scalar_rh_factorize_arc(&g_scalar, Complex::new(x, -1e-6), 1e-9);
    assert!((xp / xm - g_scalar(x)).norm() < 1e-4);

    // 2x2 Matrix RH with non-commuting jump G(t) = [[1, a(t)], [b(t), 1 + a(t)*b(t)]]
    // where a(t) = 0.3*(1 - t^2), b(t) = 0.2*i*(1 - t^2)^2 (G(+-1) = I, det G = 1)
    let jump_g = |t: f64| {
        let s = 1.0 - t * t;
        let a = Complex::new(0.3 * s, 0.0);
        let b = Complex::new(0.0, 0.2 * s * s);
        Matrix2::new(Complex::new(1.0, 0.0), a, b, Complex::new(1.0, 0.0) + a * b)
    };
    let sol = solve_matrix_rh_chebyshev(jump_g, 24).unwrap();
    // Asymptotic normalization Y(z) -> I as |z| -> oo
    let y_far = sol.eval(Complex::new(0.0, 50.0));
    assert!((y_far.sub(Matrix2::identity())).norm() < 0.05);

    // Jump condition Y^+(x) ≈ Y^-(x) * G(x) near x = 0.25
    let x_test = 0.25;
    let y_plus = sol.eval(Complex::new(x_test, 1e-4));
    let y_minus = sol.eval(Complex::new(x_test, -1e-4));
    let resid = y_plus.sub(y_minus.mul(jump_g(x_test))).norm();
    assert!(resid < 5e-3, "RH matrix jump residual={resid}");
}

#[test]
fn symbolic_phase9_transforms_and_complex_risch() {
    assert_eq!(eval_sym("mellin_transform(exp(-x), x, s)"), "gamma(s)");
    assert_eq!(
        eval_sym("inverse_mellin_transform(gamma(s), s, x)"),
        "exp(-x)"
    );

    let me1 = eval_sym("meijer_g(list(), list(0), 0, 1, z)");
    assert_eq!(me1, "exp(-z)");
    let me2 = eval_sym("meijer_g(list(0), list(0), 1, 1, z)");
    assert_eq!(me2, "1/(z + 1)");

    let mb = eval_sym("mellin_barnes_to_meijer_g(list(0), list(), list(), list(), z)");
    assert_eq!(mb, "exp(-z)");

    assert_eq!(eval_sym("hilbert_transform(cos(x), x, y)"), "sin(y)");
    assert_eq!(eval_sym("hilbert_transform(sin(x), x, y)"), "-cos(y)");

    let lap = eval_sym("bromwich_inverse_laplace(1/(s + 2), s, t)");
    assert_eq!(lap, "exp(-2*t)");

    let r1 = eval_sym("complex_risch_integrate(1/(z^2 + 1), z)");
    assert!(r1.contains("ln("), "expected complex log form, got {r1}");

    let r2 = eval_sym("complex_risch_integrate(exp(I*z)*sin(z), z)");
    assert!(
        !r2.contains("complex_risch_integrate") && !r2.contains("integral"),
        "got {r2}"
    );

    // Regression tests for Phase 9 fixes
    // 1. Abscissa-shifted Bromwich for RHP pole 1/(s - 5) -> exp(5*t)
    let f_rhp = |s: Complex<f64>| Complex::new(1.0, 0.0) / (s - Complex::new(5.0, 0.0));
    let t_val = 0.6_f64;
    let inv_t =
        rssn::kernels::complex_transforms::bromwich_talbot_shifted(f_rhp, t_val, 6.0, 32).unwrap();
    let inv_w = rssn::kernels::complex_transforms::bromwich_weideman_parabola_shifted(
        f_rhp, t_val, 6.0, 40,
    )
    .unwrap();
    let want_e5 = (5.0 * t_val).exp();
    assert!(
        (inv_t.re - want_e5).abs() < 1e-8,
        "talbot shifted: {inv_t} vs {want_e5}"
    );
    assert!(
        (inv_w.re - want_e5).abs() < 1e-6,
        "weideman shifted: {inv_w} vs {want_e5}"
    );

    // 2. Gaver-Stehfest caps m_even > 14
    assert!(gaver_stehfest_inverse_laplace(|s| 1.0 / s, 1.0, 16).is_err());

    // 3. Stieltjes transform of exp(-t) at z = 1: int_0^oo exp(-t)/(t+1) dt = e * E_1(1)
    let stj = rssn::kernels::complex_transforms::stieltjes_transform(
        |t| Complex::new((-t).exp(), 0.0),
        Complex::new(1.0, 0.0),
        30.0,
        1e-9,
    );
    let want_stj = 1.0_f64.exp() * rssn::kernels::special::exp_integral_e1(1.0);
    assert!(
        (stj.re - want_stj).abs() < 1e-6,
        "stieltjes: {stj} vs {want_stj}"
    );

    // 4. Symbolic Bromwich 1/sqrt(s) and repeated-pole complex_risch_integrate
    let br_sqrt = eval_sym("bromwich_inverse_laplace(1/s^(1/2), s, t)");
    assert!(
        br_sqrt.contains("pi") && br_sqrt.contains("t"),
        "bromwich 1/sqrt(s): {br_sqrt}"
    );
    let r_rep = eval_sym("complex_risch_integrate(1/(z^2 + 1)^2, z)");
    assert!(
        !r_rep.contains("complex_risch_integrate") && !r_rep.contains("integral"),
        "repeated pole: {r_rep}"
    );
}
