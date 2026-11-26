//! Phase 2 & Phase 3 integration tests: coordinates, Riemann sphere,
//! custom branch cuts, Riemann surface topology, full symbolic real/imag split,
//! Wirtinger derivatives, Cauchy-Riemann, harmonic conjugates, and Cauchy-FFT.

use std::f64::consts::PI;

use num_complex::Complex64;
use rssn::api::Session;
use rssn::kernels::complex as num_complex_k;
use rssn::kernels::complex_coords;
use rssn::prelude::Config;

fn eval_sym(src: &str) -> String {
    let session = Session::new();
    let term = session.parse(src).unwrap_or_else(|e| panic!("{src}: {e}"));
    let ans = session
        .compute(term, &Config::default())
        .unwrap_or_else(|e| panic!("{src}: {e}"));
    ans.term.to_string()
}

#[test]
fn polar_and_riemann_sphere_round_trips() {
    let z = Complex64::new(-3.0, 4.0);
    let (r, theta) = complex_coords::to_polar(z);
    assert!((r - 5.0).abs() < 1e-14);
    assert!((complex_coords::from_polar(r, theta) - z).norm() < 1e-13);

    for pt in [
        Complex64::new(0.0, 0.0),
        Complex64::new(1.0, -2.0),
        Complex64::new(1e150, -1e150),
    ] {
        let (x, y, z3) = complex_coords::to_riemann_sphere(pt);
        assert!((x * x + y * y + z3 * z3 - 1.0).abs() < 1e-12);
        if pt.norm() < 1e100 {
            let back = complex_coords::from_riemann_sphere(x, y, z3).unwrap();
            assert!((back - pt).norm() < 1e-12);
        }
    }
    assert!(complex_coords::from_riemann_sphere(0.0, 0.0, 1.0).is_none());

    // Chordal distance invariance under z -> 1/z
    let a = Complex64::new(1.2, -0.7);
    let b = Complex64::new(-0.5, 2.3);
    let d1 = complex_coords::chordal_distance(a, b);
    let d2 = complex_coords::chordal_distance(1.0 / a, 1.0 / b);
    assert!((d1 - d2).abs() < 1e-13);
}

#[test]
fn hyperbolic_uhp_disk_isometry() {
    let z1 = Complex64::new(0.5, 1.2);
    let z2 = Complex64::new(-1.3, 2.7);
    let w1 = complex_coords::uhp_to_disk(z1).unwrap();
    let w2 = complex_coords::uhp_to_disk(z2).unwrap();
    assert!((complex_coords::disk_to_uhp(w1).unwrap() - z1).norm() < 1e-13);
    let d_uhp = complex_coords::uhp_distance(z1, z2).unwrap();
    let d_disk = complex_coords::poincare_disk_distance(w1, w2).unwrap();
    assert!((d_uhp - d_disk).abs() < 1e-12);
}

#[test]
fn custom_branch_cuts_and_analytic_continuation() {
    // Cut along positive real axis (theta0 = 0): arg in (0, 2*pi]
    let z_above = Complex64::new(1.0, 1e-6);
    let z_below = Complex64::new(1.0, -1e-6);
    assert!((complex_coords::arg_branch(z_above, 0.0) - 1e-6).abs() < 1e-8);
    assert!((complex_coords::arg_branch(z_below, 0.0) - (2.0 * PI - 1e-6)).abs() < 1e-8);
    let jump = complex_coords::log_cut(z_below, 0.0) - complex_coords::log_cut(z_above, 0.0);
    assert!((jump - Complex64::new(0.0, 2.0 * PI)).norm() < 1e-5);

    // Analytic continuation of z^(1/3) around unit circle picks up exp(2*pi*i/3) per loop
    let n_pts = 128;
    let mut path = Vec::with_capacity(n_pts * 3 + 1);
    for k in 0..=(n_pts * 3) {
        let theta = 2.0 * PI * (k as f64) / (n_pts as f64);
        path.push(Complex64::from_polar(1.0, theta));
    }
    let vals = complex_coords::analytic_continuation_root(&path, 3, 0);
    let omega = Complex64::from_polar(1.0, 2.0 * PI / 3.0);
    assert!((vals[0] - Complex64::new(1.0, 0.0)).norm() < 1e-12);
    assert!((vals[n_pts] - omega).norm() < 1e-12);
    assert!((vals[3 * n_pts] - Complex64::new(1.0, 0.0)).norm() < 1e-12);
}

#[test]
fn riemann_surface_genus_and_monodromy() {
    // Elliptic curve w^2 = (z-1)(z-2)(z-3): genus 1, branched at infinity
    let ell = complex_coords::algebraic_riemann_surface_hyperelliptic(&[
        Complex64::new(1.0, 0.0),
        Complex64::new(2.0, 0.0),
        Complex64::new(3.0, 0.0),
    ]);
    assert_eq!(ell.genus, 1);
    assert!(ell.has_branch_at_infinity);

    // Genus 2 hyperelliptic curve: 6 roots, not branched at infinity
    let g2 = complex_coords::algebraic_riemann_surface_hyperelliptic(&[
        Complex64::new(1.0, 0.0),
        Complex64::new(2.0, 0.0),
        Complex64::new(3.0, 0.0),
        Complex64::new(4.0, 0.0),
        Complex64::new(5.0, 0.0),
        Complex64::new(6.0, 0.0),
    ]);
    assert_eq!(g2.genus, 2);
    assert!(!g2.has_branch_at_infinity);

    // Superelliptic w^3 = z^3 - 1 (3 roots, n=3 -> gcd(3,3)=3 -> unbranched at inf, g=1)
    let fermat = complex_coords::algebraic_riemann_surface_superelliptic(
        &[
            Complex64::new(1.0, 0.0),
            Complex64::new(-0.5, 0.75_f64.sqrt()),
            Complex64::new(-0.5, -0.75_f64.sqrt()),
        ],
        3,
    );
    assert_eq!(fermat.genus, 1);
    assert!(!fermat.has_branch_at_infinity);

    // Monodromy of w^2 - z = 0 around z = 0 is the transposition [1, 0]
    let poly = vec![
        vec![Complex64::new(0.0, 0.0), Complex64::new(-1.0, 0.0)], // -z * w^0
        vec![Complex64::new(0.0, 0.0)],                            // 0 * w^1
        vec![Complex64::new(1.0, 0.0)],                            // 1 * w^2
    ];
    let loop_pts: Vec<Complex64> = (0..64)
        .map(|k| Complex64::from_polar(1.0, 2.0 * PI * (k as f64) / 64.0))
        .collect();
    let perm = complex_coords::monodromy_around_loop(&poly, &loop_pts).unwrap();
    assert_eq!(perm, vec![1, 0]);
}

#[test]
fn symbolic_phase2_and_phase3_operators() {
    assert_eq!(eval_sym("riemann_genus(2, 3)"), "1");
    assert_eq!(eval_sym("riemann_genus(2, 6)"), "2");
    assert_eq!(eval_sym("riemann_genus(3, 3)"), "1");
    assert_eq!(eval_sym("abs_sq(3 + 4*I)"), "25");
    assert_eq!(eval_sym("cis(0)"), "1");
    assert_eq!(eval_sym("from_polar(2, pi/2)"), "2*I");
    assert_eq!(eval_sym("branch_cut_jump(ln(z), z, -1)"), "2*I*pi");

    // Wirtinger & Cauchy-Riemann & harmonic conjugates
    assert_eq!(eval_sym("wirtinger_dzbar(z^3 + exp(z), z)"), "0");
    assert_eq!(eval_sym("wirtinger_dz(z^2, z)"), "2*z");
    assert_eq!(eval_sym("complex_laplacian(z*conj(z), z)"), "4");
    assert_eq!(eval_sym("cauchy_riemann_cartesian(z^2, z)"), "list(0, 0)");
    assert_eq!(
        eval_sym("cauchy_riemann_polar(z^2, r, theta)"),
        "list(0, 0)"
    );
    assert_eq!(eval_sym("harmonic_conjugate(x^2 - y^2, x, y)"), "2*x*y");
    assert_eq!(eval_sym("is_entire(exp(z) + z^2, z)"), "true");
    assert_eq!(eval_sym("is_meromorphic(1/(z^2 + 1), z)"), "true");
    assert_eq!(eval_sym("is_entire(1/(z^2 + 1), z)"), "false");
    assert_eq!(eval_sym("is_meromorphic(conj(z), z)"), "false");
    assert_eq!(eval_sym("complex_diff(z^4, z, 2)"), "12*z^2");
}

#[test]
fn numeric_wirtinger_complex_step_and_cauchy_fft() {
    let z = Complex64::new(0.4, -0.3);
    // f(z) = z * conj(z) = |z|^2 => df/dz = conj(z), df/dzbar = z
    let f = |w: Complex64| w * w.conj();
    assert!((num_complex_k::wirtinger_dz_numeric(f, z, 1e-6) - z.conj()).norm() < 1e-9);
    assert!((num_complex_k::wirtinger_dzbar_numeric(f, z, 1e-6) - z).norm() < 1e-9);
    assert!(!num_complex_k::check_cauchy_riemann_numeric(
        f, z, 1e-6, 1e-8
    ));
    assert!(num_complex_k::check_cauchy_riemann_numeric(
        |w| w.sin(),
        z,
        1e-6,
        1e-8
    ));

    // Complex-step derivative of exp(x)*sin(x) at x = 0.7
    let x: f64 = 0.7;
    let exact = x.exp() * (x.sin() + x.cos());
    let cs = num_complex_k::complex_step_derivative(|w| w.exp() * w.sin(), x, 1e-20);
    assert!((cs - exact).abs() < 1e-14);

    // Cauchy-FFT derivatives of 1/(1 - z) at z0 = 0.2: f^(k)(0.2) = k! / (0.8)^(k+1)
    let z0 = Complex64::new(0.2, 0.0);
    let derivs = num_complex_k::cauchy_fft_derivatives(
        |w| Complex64::new(1.0, 0.0) / (Complex64::new(1.0, 0.0) - w),
        z0,
        0.3,
        8,
        64,
    );
    let mut fact = 1.0_f64;
    for (k, &got) in derivs.iter().enumerate() {
        if k > 0 {
            fact *= k as f64;
        }
        let want = Complex64::new(fact / (0.8_f64.powi(k as i32 + 1)), 0.0);
        assert!(
            (got - want).norm() < 1e-11 * (1.0 + want.norm()),
            "order {k}: {got} vs {want}"
        );
    }
}
