use num_complex::Complex;
use rssn::api::Config;
use rssn::api::Session;
use rssn::kernels::complex::Mobius;
use rssn::kernels::complex_roots::aberth_ehrlich;
use rssn::kernels::complex_roots::argument_principle_moments;
use rssn::kernels::complex_roots::argument_principle_winding;
use rssn::kernels::complex_roots::cauchy_root_bound;
use rssn::kernels::complex_roots::cauchy_sharp_root_bound;
use rssn::kernels::complex_roots::delves_lyness_roots;
use rssn::kernels::complex_roots::durand_kerner;
use rssn::kernels::complex_roots::enestrom_kakeya_bounds;
use rssn::kernels::complex_roots::fujiwara_root_bound;
use rssn::kernels::complex_roots::halley_root;
use rssn::kernels::complex_roots::homotopy_track_1d;
use rssn::kernels::complex_roots::laguerre_root;
use rssn::kernels::complex_roots::muller_root;
use rssn::kernels::complex_roots::rouche_check_circle;
use rssn::kernels::conformal::MobiusClass;
use rssn::kernels::conformal::bergman_kernel_disk;
use rssn::kernels::conformal::bergman_kernel_polynomial_basis;
use rssn::kernels::conformal::classify_mobius;
use rssn::kernels::conformal::conformal_scale_and_rotation;
use rssn::kernels::conformal::cross_ratio;
use rssn::kernels::conformal::disk_automorphism;
use rssn::kernels::conformal::is_conformal_numeric;
use rssn::kernels::conformal::joukowsky;
use rssn::kernels::conformal::joukowsky_inverse;
use rssn::kernels::conformal::karman_trefftz;
use rssn::kernels::conformal::koebe_map;
use rssn::kernels::conformal::mobius_fixed_points;
use rssn::kernels::conformal::mobius_from_three_points;
use rssn::kernels::conformal::solve_sc_uhp;
use rssn::kernels::conformal::szego_kerzman_stein_map;
use rssn::kernels::univalent::S_RADIUS_OF_CONVEXITY;
use rssn::kernels::univalent::S_RADIUS_OF_STARLIKENESS;
use rssn::kernels::univalent::check_bieberbach_de_branges;
use rssn::kernels::univalent::check_grunsky_inequality;
use rssn::kernels::univalent::check_nehari_univalence;
use rssn::kernels::univalent::fekete_szego_functional;
use rssn::kernels::univalent::grunsky_matrix;
use rssn::kernels::univalent::hankel_determinant;
use rssn::kernels::univalent::koebe_distortion_bounds;
use rssn::kernels::univalent::koebe_growth_bounds;
use rssn::kernels::univalent::koebe_rotation_bound;
use rssn::kernels::univalent::loewner_slit_step;
use rssn::kernels::univalent::radius_of_convexity_numeric;
use rssn::kernels::univalent::radius_of_starlikeness_numeric;
use rssn::kernels::univalent::schwarzian_derivative_numeric;
use std::f64::consts::PI;

fn eval_sym(src: &str) -> String {
    let session = Session::new();
    let term = session.parse(src).unwrap_or_else(|e| panic!("{src}: {e}"));
    let ans = session
        .compute(term, &Config::default())
        .unwrap_or_else(|e| panic!("{src}: {e}"));
    ans.term.to_string()
}

fn eval_poly(
    coeffs: &[Complex<f64>],
    z: Complex<f64>,
) -> Complex<f64> {
    let mut acc = Complex::new(0.0, 0.0);
    for &c in coeffs.iter().rev() {
        acc = acc * z + c;
    }
    acc
}

#[test]
fn root_bounds_and_polynomial_solvers() {
    // p(z) = (z - 1)(z - 2i)(z + 1 + i) = -2 - 2i + (3i)z + (-i)z^2 + z^3
    let r1 = Complex::new(1.0, 0.0);
    let r2 = Complex::new(0.0, 2.0);
    let r3 = Complex::new(-1.0, -1.0);
    let c0 = -r1 * r2 * r3;
    let c1 = r1 * r2 + r2 * r3 + r3 * r1;
    let c2 = -(r1 + r2 + r3);
    let c3 = Complex::new(1.0, 0.0);
    let coeffs = vec![c0, c1, c2, c3];

    let cb = cauchy_root_bound(&coeffs).unwrap();
    let cs = cauchy_sharp_root_bound(&coeffs).unwrap();
    let fb = fujiwara_root_bound(&coeffs).unwrap();
    let max_root = r1.norm().max(r2.norm()).max(r3.norm());
    assert!(cb >= max_root && cs >= max_root && fb >= max_root);
    assert!(cs <= cb + 1e-12);

    let (ek_lo, ek_hi) = enestrom_kakeya_bounds(&[1.0, 2.0, 3.0, 4.0]).unwrap();
    assert!((ek_lo - 0.5).abs() < 1e-12 && (ek_hi - 0.75).abs() < 1e-12);

    let roots_ae = aberth_ehrlich(&coeffs, 200, 1e-13).unwrap();
    let roots_dk = durand_kerner(&coeffs, 300, 1e-12).unwrap();
    assert_eq!(roots_ae.len(), 3);
    assert_eq!(roots_dk.len(), 3);
    for z in &roots_ae {
        assert!(eval_poly(&coeffs, *z).norm() < 1e-11);
    }
    for z in &roots_dk {
        assert!(eval_poly(&coeffs, *z).norm() < 1e-10);
    }

    let r_lag = laguerre_root(&coeffs, Complex::new(0.2, 2.1), 50, 1e-13).unwrap();
    assert!(eval_poly(&coeffs, r_lag).norm() < 1e-12);

    let r_mul = muller_root(
        |z| z.cos() - z,
        Complex::new(0.5, 0.1),
        Complex::new(0.7, 0.2),
        Complex::new(0.8, 0.0),
        100,
        1e-13,
    )
    .unwrap();
    assert!((r_mul.cos() - r_mul).norm() < 1e-12);

    let r_hal = halley_root(
        |z| z.exp() - Complex::new(2.0, 1.0),
        |z| z.exp(),
        |z| z.exp(),
        Complex::new(0.8, 0.4),
        50,
        1e-13,
    )
    .unwrap();
    assert!((r_hal.exp() - Complex::new(2.0, 1.0)).norm() < 1e-13);
}

#[test]
fn argument_principle_delves_lyness_rouche_and_homotopy() {
    // f(z) = (z - (0.3 + 0.4i)) * (z - (-0.5 + 0.2i)) * exp(z)
    let z1 = Complex::new(0.3, 0.4);
    let z2 = Complex::new(-0.5, 0.2);
    let f = |z: Complex<f64>| (z - z1) * (z - z2) * z.exp();
    let df = |z: Complex<f64>| ((z - z1) + (z - z2) + (z - z1) * (z - z2)) * z.exp();

    let w = argument_principle_winding(f, df, Complex::new(0.0, 0.0), 1.5, 128).unwrap();
    assert_eq!(w, 2);

    let moms = argument_principle_moments(f, df, Complex::new(0.0, 0.0), 1.5, 2, 128).unwrap();
    assert!((moms[0] - Complex::new(2.0, 0.0)).norm() < 1e-10);
    assert!((moms[1] - (z1 + z2)).norm() < 1e-10);
    assert!((moms[2] - (z1 * z1 + z2 * z2)).norm() < 1e-10);

    let found = delves_lyness_roots(f, df, Complex::new(0.0, 0.0), 1.5, 128).unwrap();
    assert_eq!(found.len(), 2);
    let d1 = found
        .iter()
        .map(|&r| (r - z1).norm())
        .fold(f64::INFINITY, f64::min);
    let d2 = found
        .iter()
        .map(|&r| (r - z2).norm())
        .fold(f64::INFINITY, f64::min);
    assert!(d1 < 1e-10 && d2 < 1e-10);

    // Rouche: on |z| = 1, |3z| = 3 > |z^5 + 1| <= 2
    assert!(
        rouche_check_circle(
            |z| z * 3.0,
            |z| z.powi(5) + Complex::new(1.0, 0.0),
            Complex::new(0.0, 0.0),
            1.0,
            128
        )
        .unwrap()
    );

    // Homotopy: H(z, t) = (1-t)(z^3 - 1) + t(z^3 + 2z - 5) = z^3 - 1 + t*(2z - 4)
    let z_end = homotopy_track_1d(
        |z, t| z.powi(3) - Complex::new(1.0, 0.0) + (z * 2.0 - Complex::new(4.0, 0.0)) * t,
        |z, t| z * z * 3.0 + Complex::new(2.0 * t, 0.0),
        |z, _t| z * 2.0 - Complex::new(4.0, 0.0),
        Complex::new(1.0, 0.0),
        50,
        10,
        1e-13,
    )
    .unwrap();
    assert!((z_end.powi(3) + z_end * 2.0 - Complex::new(5.0, 0.0)).norm() < 1e-12);
}

#[test]
fn conformal_maps_sc_szego_and_bergman() {
    let z1 = Complex::new(0.0, 0.0);
    let z2 = Complex::new(1.0, 0.0);
    let z3 = Complex::new(0.0, 1.0);
    let w1 = Complex::new(1.0, 1.0);
    let w2 = Complex::new(2.0, -1.0);
    let w3 = Complex::new(-1.0, 0.5);
    let m = mobius_from_three_points(z1, z2, z3, w1, w2, w3).unwrap();
    assert!((m.apply(z1) - w1).norm() < 1e-12);
    assert!((m.apply(z2) - w2).norm() < 1e-12);
    assert!((m.apply(z3) - w3).norm() < 1e-12);

    let z4 = Complex::new(2.0, 3.0);
    assert!(
        (cross_ratio(z1, z2, z3, z4)
            - cross_ratio(m.apply(z1), m.apply(z2), m.apply(z3), m.apply(z4)))
        .norm()
            < 1e-11
    );

    let ell = Mobius {
        a: Complex::from_polar(1.0, 0.5),
        b: Complex::new(0.0, 0.0),
        c: Complex::new(0.0, 0.0),
        d: Complex::from_polar(1.0, -0.5),
    };
    assert_eq!(classify_mobius(&ell).unwrap(), MobiusClass::Elliptic);
    assert_eq!(mobius_fixed_points(&ell), vec![Complex::new(0.0, 0.0)]);

    let auto = disk_automorphism(Complex::new(0.3, -0.4), 0.7).unwrap();
    assert!((auto.apply(Complex::new(0.3, -0.4))).norm() < 1e-12);
    assert!((auto.apply(Complex::from_polar(1.0, 1.2)).norm() - 1.0).abs() < 1e-12);

    let z = Complex::new(1.5, 0.4);
    let w = joukowsky(z);
    let (r1, r2) = joukowsky_inverse(w);
    assert!((r1 - z).norm() < 1e-12 || (r2 - z).norm() < 1e-12);
    let _ = karman_trefftz(Complex::new(1.5, 0.2), 1.95);
    assert!((koebe_map(Complex::new(0.5, 0.0)) - Complex::new(2.0, 0.0)).norm() < 1e-12);

    assert!(is_conformal_numeric(
        &|x| x.exp(),
        Complex::new(0.3, 0.4),
        1e-5,
        1e-6
    ));
    assert!(!is_conformal_numeric(
        &|x| x * x,
        Complex::new(0.0, 0.0),
        1e-5,
        1e-6
    ));
    let (sc, rot) = conformal_scale_and_rotation(&|x| x * x, Complex::new(1.0, 1.0), 1e-5);
    assert!((sc - 2.0 * 2.0_f64.sqrt()).abs() < 1e-7);
    assert!((rot - PI / 4.0).abs() < 1e-7);

    // Schwarz-Christoffel on unit square [0, 1] x [0, 1]
    let square = [
        Complex::new(0.0, 0.0),
        Complex::new(1.0, 0.0),
        Complex::new(1.0, 1.0),
        Complex::new(0.0, 1.0),
    ];
    let sc_map = solve_sc_uhp(&square, 40, 1e-7).unwrap();
    let img = sc_map.eval(Complex::new(0.5, 0.5)).unwrap();
    assert!(
        img.re > 0.0 && img.re < 1.0 && img.im > 0.0 && img.im < 1.0,
        "img={img}"
    );

    // Kerzman-Stein on unit disk centered at 0, basepoint a = 0.2 -> exact Riemann map is (z - a)/(1 - conj(a)*z)
    let a = Complex::new(0.2, 0.0);
    let pt = Complex::new(0.4, 0.3);
    let (phi_z, s_za) = szego_kerzman_stein_map(
        |t| Complex::from_polar(1.0, t),
        |t| Complex::new(-t.sin(), t.cos()),
        64,
        a,
        pt,
    )
    .unwrap();
    let exact_phi = (pt - a) / (Complex::new(1.0, 0.0) - a.conj() * pt);
    let exact_szego =
        Complex::new(1.0 / (2.0 * PI), 0.0) / (Complex::new(1.0, 0.0) - pt * a.conj());
    assert!((phi_z - exact_phi).norm() < 1e-5, "{phi_z} vs {exact_phi}");
    assert!(
        (s_za - exact_szego).norm() < 1e-5,
        "{s_za} vs {exact_szego}"
    );

    // Bergman kernel on disk vs orthonormal polynomial series
    let z_b = Complex::new(0.2, 0.1);
    let w_b = Complex::new(0.1, -0.2);
    let exact_b = bergman_kernel_disk(z_b, w_b);
    let basis: Vec<Vec<Complex<f64>>> = (0..25)
        .map(|k| {
            let mut v = vec![Complex::new(0.0, 0.0); k + 1];
            v[k] = Complex::new(((k + 1) as f64 / PI).sqrt(), 0.0);
            v
        })
        .collect();
    let approx_b = bergman_kernel_polynomial_basis(&basis, z_b, w_b);
    assert!((exact_b - approx_b).norm() < 1e-10);
}

#[test]
fn univalent_functions_theorems_and_loewner() {
    // Koebe coefficients a_k = k
    let koebe_coeffs: Vec<Complex<f64>> = (0..=8).map(|k| Complex::new(k as f64, 0.0)).collect();
    assert!(check_bieberbach_de_branges(&koebe_coeffs, 1e-12));
    assert!(
        (fekete_szego_functional(
            Complex::new(2.0, 0.0),
            Complex::new(3.0, 0.0),
            Complex::new(1.0, 0.0)
        ) - 1.0)
            .abs()
            < 1e-12
    );

    // Hankel determinant H_2(1) for Koebe = det [[1, 2], [2, 3]] = -1
    let h21 = hankel_determinant(&koebe_coeffs, 2, 1).unwrap();
    assert!((h21 - Complex::new(-1.0, 0.0)).norm() < 1e-12);

    let (g_lo, g_hi) = koebe_growth_bounds(0.5).unwrap();
    assert!((g_lo - 0.5 / 2.25).abs() < 1e-12 && (g_hi - 2.0).abs() < 1e-12);
    let (d_lo, d_hi) = koebe_distortion_bounds(0.5).unwrap();
    assert!((d_lo - 0.5 / 3.375).abs() < 1e-12 && (d_hi - 1.5 / 0.125).abs() < 1e-12);
    assert!(koebe_rotation_bound(0.8).unwrap() > 4.0 * 0.5_f64.asin());

    // Schwarzian of any Mobius map is 0
    let mob = |z: Complex<f64>| (z * 2.0 + Complex::new(1.0, 0.0)) / (z + Complex::new(3.0, -1.0));
    let s_mob = schwarzian_derivative_numeric(&mob, Complex::new(0.2, 0.1), 0.05).unwrap();
    assert!(s_mob.norm() < 1e-9, "s_mob={s_mob}");

    assert!(check_nehari_univalence(&mob, &[0.2, 0.5, 0.8], 8, 0.05).unwrap());

    // Koebe starlikeness radius is 1 (starlike on all of D), convexity radius is 2 - sqrt(3)
    let r_star = radius_of_starlikeness_numeric(
        &koebe_map,
        &|z| (Complex::new(1.0, 0.0) + z) / (Complex::new(1.0, 0.0) - z).powi(3),
        0.99,
        64,
    );
    assert!(r_star > 0.98 && S_RADIUS_OF_STARLIKENESS > 0.65);

    let r_conv = radius_of_convexity_numeric(
        &|z| (Complex::new(1.0, 0.0) + z) / (Complex::new(1.0, 0.0) - z).powi(3),
        &|z| (Complex::new(4.0, 0.0) + z * 2.0) / (Complex::new(1.0, 0.0) - z).powi(4),
        0.99,
        64,
    );
    assert!(
        (r_conv - S_RADIUS_OF_CONVEXITY).abs() < 1e-3,
        "r_conv={r_conv}"
    );

    // Grunsky inequality for Koebe: G_{1,1} = +1, |G x x| <= ||x||^2
    let g_mat = grunsky_matrix(&koebe_coeffs, 3).unwrap();
    assert!((g_mat[0][0] - Complex::new(1.0, 0.0)).norm() < 1e-12);
    assert!(
        check_grunsky_inequality(
            &g_mat,
            &[
                Complex::new(0.5, 0.1),
                Complex::new(-0.3, 0.2),
                Complex::new(0.1, -0.4)
            ],
            1e-9
        )
        .unwrap()
    );

    // Loewner step keeps points inside unit disk moving inward
    let w0 = Complex::new(0.3, 0.4);
    let w1 = loewner_slit_step(w0, 0.0, 0.01);
    assert!(w1.norm() < 1.0 && (w1 - w0).norm() > 0.0);
}

#[test]
fn symbolic_phase7_and_phase8_operators() {
    assert_eq!(eval_sym("joukowsky(2)"), "5/4");
    assert_eq!(eval_sym("koebe(1/2)"), "2");
    assert_eq!(eval_sym("cayley_transform(I)"), "0");
    assert_eq!(eval_sym("cross_ratio(0, 1, 2, 3)"), "4/3");

    // Schwarzian of Mobius (2*z + 1)/(z + 3) is 0
    let sw = eval_sym("schwarzian_derivative((2*z + 1)/(z + 3), z)");
    assert_eq!(sw, "0", "sw={sw}");

    assert_eq!(eval_sym("is_conformal_at(z^2, z, 1)"), "1");
    assert_eq!(eval_sym("is_conformal_at(z^2, z, 0)"), "0");
    assert_eq!(eval_sym("rouche_compare(3*z, z^5 + 1, z, 0, 1)"), "1");
    assert_eq!(eval_sym("cauchy_root_bound(z^3 - 2*z + 5, z)"), "6");

    // Regression tests for Phase 7 fixes
    // 1. Verified quadtree + Krawczyk root isolation for z^3 - 1 = 0
    let balls = rssn::kernels::complex_roots::isolate_roots_verified(
        &[
            Complex::new(-1.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(1.0, 0.0),
        ],
        1e-6,
    );
    assert_eq!(balls.len(), 3, "expected 3 isolated root balls");

    // 2. Delves-Lyness rejects poles inside disk
    let pole_err = delves_lyness_roots(
        |z| Complex::new(1.0, 0.0) / z,
        |z| -Complex::new(1.0, 0.0) / (z * z),
        Complex::new(0.0, 0.0),
        1.0,
        64,
    );
    assert!(pole_err.is_err());

    // 3. solve_complex with complex coefficients: z^2 - I = 0
    let sol_c = eval_sym("solve_complex(z^2 - I, z)");
    assert!(
        sol_c.contains("0.7071") && sol_c.contains("I"),
        "solve_complex(z^2-I)={sol_c}"
    );
}
