use num_complex::Complex;
use rssn::api::Config;
use rssn::api::Session;
use rssn::kernels::complex_dynamics::FixedPointKind;
use rssn::kernels::complex_dynamics::bottcher_coordinate;
use rssn::kernels::complex_dynamics::classify_periodic_point;
use rssn::kernels::complex_dynamics::external_ray_point;
use rssn::kernels::complex_dynamics::find_periodic_points_quadratic;
use rssn::kernels::complex_dynamics::julia_escape_time;
use rssn::kernels::complex_dynamics::koenigs_coordinate;
use rssn::kernels::complex_dynamics::lyapunov_exponent;
use rssn::kernels::complex_dynamics::mandelbrot_escape_time;
use rssn::kernels::complex_dynamics::orbit;
use std::f64::consts::PI;

fn eval_sym(src: &str) -> String {
    let session = Session::new();
    let term = session.parse(src).unwrap_or_else(|e| panic!("{src}: {e}"));
    let ans = session
        .compute(term, &Config::default())
        .unwrap_or_else(|e| panic!("{src}: {e}"));
    ans.term.to_string()
}

#[test]
fn orbits_periodic_points_and_classification() {
    let c = Complex::new(-1.0, 0.0);
    let f = |z: Complex<f64>| z * z + c;
    let df = |z: Complex<f64>| z * 2.0;

    let orb = orbit(f, Complex::new(0.0, 0.0), 4);
    assert_eq!(orb[0], Complex::new(0.0, 0.0));
    assert_eq!(orb[1], Complex::new(-1.0, 0.0));
    assert_eq!(orb[2], Complex::new(0.0, 0.0));

    let (mult2, kind2) = classify_periodic_point(f, df, Complex::new(0.0, 0.0), 2);
    assert_eq!(kind2, FixedPointKind::Superattracting);
    assert!(mult2.norm() < 1e-12);

    let phi_golden = 0.5 * (1.0 + 5.0_f64.sqrt());
    let (mult_fp, kind_fp) = classify_periodic_point(f, df, Complex::new(phi_golden, 0.0), 1);
    assert_eq!(kind_fp, FixedPointKind::Repelling);
    assert!((mult_fp.re - 2.0 * phi_golden).abs() < 1e-12);

    let (_, kind_att) = classify_periodic_point(
        |z| z * 0.5 + Complex::new(0.0, 0.1),
        |_| Complex::new(0.5, 0.0),
        Complex::new(0.0, 0.2),
        1,
    );
    assert_eq!(kind_att, FixedPointKind::Attracting);

    let (_, kind_neu) = classify_periodic_point(
        |z| z * Complex::from_polar(1.0, 1.0),
        |_| Complex::from_polar(1.0, 1.0),
        Complex::new(0.0, 0.0),
        1,
    );
    assert_eq!(kind_neu, FixedPointKind::Neutral);

    let seeds = vec![
        Complex::new(0.1, 0.1),
        Complex::new(-0.9, -0.1),
        Complex::new(1.5, 0.1),
        Complex::new(-0.6, 0.1),
    ];
    let pts = find_periodic_points_quadratic(c, 2, &seeds, 40, 1e-12);
    assert!(pts.iter().any(|&p| p.norm() < 1e-10));
    assert!(
        pts.iter()
            .any(|&p| (p - Complex::new(-1.0, 0.0)).norm() < 1e-10)
    );
}

#[test]
fn green_potential_bottcher_koenigs_lyapunov_and_external_ray() {
    // Mandelbrot: c = 0 is inside, c = 2 escapes
    let (m_in, g_in) = mandelbrot_escape_time(Complex::new(0.0, 0.0), 100, 1e4);
    assert_eq!(m_in, 100);
    assert_eq!(g_in, 0.0);
    let (m_out, g_out) = mandelbrot_escape_time(Complex::new(2.0, 0.0), 100, 1e4);
    assert!(m_out < 100 && g_out > 0.0);

    // Green potential functional equation: G_c(z^2 + c) = 2 * G_c(z)
    let c = Complex::new(-0.12, 0.74);
    let z = Complex::new(3.2, 1.5);
    let (_, gz) = julia_escape_time(z, c, 200, 1e12);
    let (_, gfz) = julia_escape_time(z * z + c, c, 200, 1e12);
    assert!(
        (gfz - 2.0 * gz).abs() < 1e-10,
        "gfz={gfz}, 2*gz={}",
        2.0 * gz
    );

    // Bottcher functional equation: phi_c(z^2 + c) = phi_c(z)^2
    let phi_z = bottcher_coordinate(z, c, 32).unwrap();
    let phi_fz = bottcher_coordinate(z * z + c, c, 32).unwrap();
    assert!((phi_fz - phi_z * phi_z).norm() < 1e-10);

    // Koenigs functional equation: psi(f(z)) = lambda * psi(z)
    let lam = Complex::new(0.4, 0.2);
    let f_k = |w: Complex<f64>| lam * w + w * w * 0.3;
    let zk = Complex::new(0.15, -0.1);
    let psi_z = koenigs_coordinate(f_k, Complex::new(0.0, 0.0), lam, zk, 45).unwrap();
    let psi_fz = koenigs_coordinate(f_k, Complex::new(0.0, 0.0), lam, f_k(zk), 45).unwrap();
    assert!((psi_fz - lam * psi_z).norm() < 1e-10);

    // Lyapunov exponent of z^2 - 2 on [-2, 2] is ln(2)
    let lyap = lyapunov_exponent(
        |w| w * w - Complex::new(2.0, 0.0),
        |w| w * 2.0,
        Complex::new(0.31415926, 0.0),
        4000,
        200,
    )
    .unwrap();
    assert!((lyap - 2.0_f64.ln()).abs() < 0.05, "lyap={lyap}");

    // External ray at potential r = 0.5, angle 1/7
    let r = 0.5_f64;
    let theta = 1.0 / 7.0;
    let ray_pt = external_ray_point(c, theta, r, 10).unwrap();
    let phi_ray = bottcher_coordinate(ray_pt, c, 24).unwrap();
    let target = Complex::from_polar(r.exp(), 2.0 * PI * theta);
    assert!(
        (phi_ray - target).norm() < 1e-9,
        "phi_ray={phi_ray} vs target={target}"
    );
}

#[test]
fn symbolic_phase10_dynamics_operators() {
    assert_eq!(eval_sym("orbit_iterate(z^2 + 1, z, 0, 3)"), "5");
    assert_eq!(eval_sym("multiplier(z^2 - 1, z, 0, 2)"), "0");
    assert_eq!(eval_sym("is_attracting_fixed_point(z/2, z, 0)"), "1");
    assert_eq!(eval_sym("is_attracting_fixed_point(2*z, z, 0)"), "0");

    let bs = eval_sym("bottcher_series(c, z, 3)");
    assert!(bs.contains("c/(2*z)") || bs.contains("1/2*c/z"), "bs={bs}");

    // Regression tests for Phase 10 fixes
    // 1. Mandelbrot Green's potential G_M(c) = G_c(c) (fixes factor-of-2 bug)
    let c_test = Complex::new(2.0, 0.0);
    let (_, g_m) = mandelbrot_escape_time(c_test, 100, 1e10);
    let (_, g_j) = julia_escape_time(c_test, c_test, 100, 1e10);
    assert!((g_m - g_j).abs() < 1e-10, "g_m={g_m} vs g_j={g_j}");

    // 2. external_ray_point does not overflow at n_steps = 15, r = 0.5
    let ray15 = external_ray_point(Complex::new(-0.12, 0.74), 1.0 / 7.0, 0.5, 15).unwrap();
    assert!(ray15.is_finite());

    // 3. Period-2 points exclude period-1 fixed points
    let p2 = find_periodic_points_quadratic(
        Complex::new(-1.0, 0.0),
        2,
        &[Complex::new(0.05, 0.01), Complex::new(1.6, 0.01)],
        40,
        1e-12,
    );
    let phi_g = 0.5 * (1.0 + 5.0_f64.sqrt());
    assert!(
        p2.iter()
            .all(|&p| (p - Complex::new(phi_g, 0.0)).norm() > 1e-3)
    );

    // 4. Newton fractal basin for z^3 - 1
    let (root_n, _) = rssn::kernels::complex_dynamics::newton_fractal_basin(
        |w| w.powi(3) - Complex::new(1.0, 0.0),
        |w| 3.0 * w * w,
        Complex::new(0.8, 0.2),
        50,
        1e-12,
    )
    .unwrap();
    assert!((root_n - Complex::new(1.0, 0.0)).norm() < 1e-10);

    // 5. Symbolic multiplier at period 10 (chain rule, no swell)
    assert_eq!(eval_sym("multiplier(z^2 - 1, z, 0, 10)"), "0");
}
