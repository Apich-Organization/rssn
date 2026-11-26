//! Phase 5 & Phase 6 integration tests: contours, Cauchy PV, adaptive
//! singularity/branch-cut deformation, Laurent FFT, Padé, Borel-Padé,
//! holonomic D-finite continuation & monodromy, and steepest descent.

use std::f64::consts::PI;

use num_complex::Complex64;
use rssn::api::Config;
use rssn::api::Session;
use rssn::kernels::analytic_continuation as ac;
use rssn::kernels::contour_deform as cd;
use rssn::kernels::steepest_descent as sd;

fn eval_sym(src: &str) -> String {
    let session = Session::new();
    let term = session.parse(src).unwrap_or_else(|e| panic!("{src}: {e}"));
    let ans = session
        .compute(term, &Config::default())
        .unwrap_or_else(|e| panic!("{src}: {e}"));
    ans.term.to_string()
}

#[test]
fn contours_keyhole_dogbone_indented_and_deformation() {
    let i_unit = Complex64::new(0.0, 1.0);
    // Keyhole around z^(1/2) / (z + 1) at branch point 0, radius 2:
    // By residue theorem at z = -1, pole is on cut if cut_angle = pi; put pole at z = 0.5*i instead:
    // Residue at 0.5*i of sqrt(z) / (z - 0.5*i) is sqrt(0.5*i) => 2*pi*i * sqrt(0.5*i).
    // Wait: keyhole also has cut jump along (-oo, 0], so test pure pole inside keyhole of 1/(z - 0.5i):
    let kh = cd::keyhole_contour(Complex64::new(0.0, 0.0), 0.1, 2.0, PI, 1e-4);
    let val = cd::integrate_contour(
        &|z| Complex64::new(1.0, 0.0) / (z - 0.5 * i_unit),
        &kh,
        1e-9,
    );
    assert!((val - 2.0 * PI * i_unit).norm() < 1e-6, "keyhole: {val}");

    // Dogbone around [-1, 1] for 1/sqrt(z^2 - 1) (standard branch):
    let db = cd::dogbone_contour(Complex64::new(-1.0, 0.0), Complex64::new(1.0, 0.0), 0.3);
    let val_db = cd::integrate_contour(&|z| Complex64::new(1.0, 0.0) / z, &db, 1e-9);
    assert!(
        (val_db - 2.0 * PI * i_unit).norm() < 1e-6,
        "dogbone: {val_db}"
    );

    // Indented semicircle for sin(x)/x on [-30, 30]
    let ind = cd::indented_semicircle_contour(30.0, &[0.0], 0.2, true);
    let val_ind = cd::integrate_contour(&|z| (i_unit * z).exp() / z, &ind, 1e-7);
    // Closed contour in UHP with pole 0 outside => integral = 0
    assert!(val_ind.norm() < 1e-5, "indented semicircle: {val_ind}");

    // Adaptive deformation around pole at z = 0.5 on [0, 1] for 1/(z - 0.5)
    let val_def = cd::integrate_deformed(
        &|z| Complex64::new(1.0, 0.0) / (z - 0.5),
        Complex64::new(0.0, 0.0),
        Complex64::new(1.0, 0.0),
        &[Complex64::new(0.5, 0.0)],
        &[],
        0.1,
        1e-9,
    );
    // P.V. = 0, plus semicircle detour +- pi*i
    assert!((val_def.re).abs() < 1e-7 && (val_def.im.abs() - PI).abs() < 1e-6);
}

#[test]
fn cauchy_pv_and_residue_at_infinity() {
    // P.V. int_{-1}^2 1/x dx = ln(2)
    let pv = cd::cauchy_principal_value(|x| Complex64::new(1.0 / x, 0.0), -1.0, 2.0, &[0.0], 1e-10);
    assert!((pv.re - 2.0_f64.ln()).abs() < 1e-9);

    // Res_{z=oo} (1/z) = -1
    let r_inf = cd::residue_at_infinity_numeric(|z| Complex64::new(1.0, 0.0) / z, 5.0, 64);
    assert!((r_inf - Complex64::new(-1.0, 0.0)).norm() < 1e-12);
}

#[test]
fn laurent_fft_pade_wynn_and_borel_pade() {
    // Laurent of exp(z)/z^3 at 0: c_{-3}=1, c_{-2}=1, c_{-1}=1/2, c_0=1/6
    let coeffs = ac::laurent_coeffs_fft(
        |z| z.exp() / z.powi(3),
        Complex64::new(0.0, 0.0),
        0.8,
        3,
        2,
        64,
    );
    let get_c = |k: i32| coeffs.iter().find(|x| x.0 == k).unwrap().1;
    assert!((get_c(-3) - Complex64::new(1.0, 0.0)).norm() < 1e-12);
    assert!((get_c(-2) - Complex64::new(1.0, 0.0)).norm() < 1e-12);
    assert!((get_c(-1) - Complex64::new(0.5, 0.0)).norm() < 1e-12);
    assert!((get_c(0) - Complex64::new(1.0 / 6.0, 0.0)).norm() < 1e-12);

    // Padé [6/6] of ln(1 + z) evaluated outside unit disk at z = 2.0 + 1.0i
    let mut taylor = Vec::with_capacity(13);
    taylor.push(Complex64::new(0.0, 0.0));
    for k in 1..=12 {
        let sign = if k % 2 == 1 { 1.0 } else { -1.0 };
        taylor.push(Complex64::new(sign / (k as f64), 0.0));
    }
    let (p, q) = ac::pade_approximant_complex(&taylor, 6, 6).unwrap();
    let z_far = Complex64::new(2.0, 1.0);
    let approx = ac::eval_pade_complex(&p, &q, z_far);
    assert!(
        (approx - (1.0 + z_far).ln()).norm() < 1e-5,
        "{approx} vs {}",
        (1.0 + z_far).ln()
    );

    // Wynn epsilon on alternating geometric series sum (-0.8)^n = 1/1.8
    let mut psums = Vec::new();
    let mut s = Complex64::new(0.0, 0.0);
    let mut r = Complex64::new(1.0, 0.0);
    for _ in 0..11 {
        s += r;
        psums.push(s);
        r *= Complex64::new(-0.8, 0.0);
    }
    let w = ac::wynn_epsilon_complex(&psums);
    assert!((w - Complex64::new(1.0 / 1.8, 0.0)).norm() < 1e-10);

    // Borel-Padé resummation of Euler's divergent series sum (-1)^n n! z^n at z = 0.5
    // Exact value = (e^2) * E_1(2)
    let mut euler_c = Vec::with_capacity(15);
    let mut fact = 1.0_f64;
    for k in 0..15 {
        if k > 0 {
            fact *= k as f64;
        }
        let sign = if k % 2 == 0 { 1.0 } else { -1.0 };
        euler_c.push(Complex64::new(sign * fact, 0.0));
    }
    let resum = ac::borel_pade_resum(&euler_c, Complex64::new(0.5, 0.0), 6, 6).unwrap();
    let exact = 2.0_f64.exp() * rssn::kernels::special::exp_integral_e1(2.0) / 0.5;
    assert!((resum.re - exact).abs() < 1e-6, "{resum} vs {exact}");
}

#[test]
fn holonomic_continuation_and_monodromy() {
    // 1. First-order ODE: z * y' - 0.5 * y = 0 => y = sqrt(z), monodromy around 0 is [-1]
    let ode1 = ac::HolonomicOde {
        coeffs: vec![
            vec![Complex64::new(-0.5, 0.0)],
            vec![Complex64::new(0.0, 0.0), Complex64::new(1.0, 0.0)],
        ],
    };
    let loop_sq = [
        Complex64::new(1.0, 0.0),
        Complex64::new(0.0, 1.0),
        Complex64::new(-1.0, 0.0),
        Complex64::new(0.0, -1.0),
    ];
    let mono1 = ac::holonomic_monodromy_matrix(&ode1, &loop_sq, 24).unwrap();
    assert!((mono1[0][0] - Complex64::new(-1.0, 0.0)).norm() < 1e-10);

    // 2. Airy ODE: y'' - z * y = 0 (entire, monodromy = I_2)
    let airy_ode = ac::HolonomicOde {
        coeffs: vec![
            vec![Complex64::new(0.0, 0.0), Complex64::new(-1.0, 0.0)],
            vec![Complex64::new(0.0, 0.0)],
            vec![Complex64::new(1.0, 0.0)],
        ],
    };
    let mono_airy = ac::holonomic_monodromy_matrix(&airy_ode, &loop_sq, 24).unwrap();
    assert!((mono_airy[0][0] - Complex64::new(1.0, 0.0)).norm() < 1e-10);
    assert!((mono_airy[1][1] - Complex64::new(1.0, 0.0)).norm() < 1e-10);
    assert!(mono_airy[0][1].norm() < 1e-10);
    assert!(mono_airy[1][0].norm() < 1e-10);
}

#[test]
fn steepest_descent_stirling_asymptotics() {
    // int_0^oo exp(lambda * (ln(z) - z)) dz = Gamma(lambda + 1) / lambda^(lambda + 1)
    // Saddle point at z0 = 1, phi(1) = -1, phi''(1) = -1
    let phi = |z: Complex64| z.ln() - z;
    let dphi = |z: Complex64| Complex64::new(1.0, 0.0) / z - 1.0;
    let d2phi = |z: Complex64| -Complex64::new(1.0, 0.0) / (z * z);
    let g = |_z: Complex64| Complex64::new(1.0, 0.0);

    let saddle = sd::find_saddle_point(phi, dphi, d2phi, Complex64::new(1.2, 0.1), 1e-12).unwrap();
    assert!((saddle.z0 - Complex64::new(1.0, 0.0)).norm() < 1e-11);

    let lambda = 50.0;
    let exact = rssn::kernels::special::gamma_numerical(lambda + 1.0) / lambda.powf(lambda + 1.0);
    let lead = sd::saddle_point_leading(g, &saddle, lambda);
    let ord2 = sd::saddle_point_asymptotic_order2(phi, g, &saddle, lambda);
    let num_sd = sd::steepest_descent_integrate(phi, g, &saddle, lambda, 0.8, 1e-10);

    // O(1/lambda) Stirling correction: 1 + 1/(12*lambda)
    assert!(
        (ord2.re - exact).abs() < 1e-4 * exact,
        "ord2={ord2}, exact={exact}"
    );
    assert!((lead.re - exact).abs() < 0.01 * exact);
    assert!((num_sd.re - exact).abs() < 1e-4 * exact);
}

#[test]
fn symbolic_phase5_and_phase6_operators() {
    assert_eq!(eval_sym("residue_at_infinity(1/z, z)"), "-1");
    assert_eq!(eval_sym("residue_real_integral(1/(x^2 + 1), x)"), "pi");
    assert_eq!(eval_sym("cauchy_pv(1/x, x, -1, 2, 0)"), "ln(2)");
    assert_eq!(
        eval_sym("laurent_principal_part(1/z + 2 + 3*z, z, 0, 2)"),
        "1/z"
    );
    assert_eq!(
        eval_sym("laurent_analytic_part(1/z + 2 + 3*z, z, 0, 2)"),
        "3*z + 2"
    );

    let trig_res = eval_sym("residue_trig_integral(1/(2 + cos(t)), t)");
    assert!(
        trig_res.contains("pi"),
        "unexpected trig integral: {trig_res}"
    );

    let ml = eval_sym("mittag_leffler_expansion(1/(z^2 - 1), z)");
    assert!(ml.contains("z"), "unexpected ML: {ml}");

    let sad = eval_sym("saddle_point_asymptotic(1, -z^2, z, 0, lam)");
    assert!(sad.contains("pi"), "unexpected saddle: {sad}");
}
