//! Coverage tests for `rssn::kernels::fractal_geometry_and_chaos`.

use rssn::kernels::fractal_geometry_and_chaos::FractalData;
use rssn::kernels::fractal_geometry_and_chaos::barnsley_fern_ifs;
use rssn::kernels::fractal_geometry_and_chaos::box_counting_dimension;
use rssn::kernels::fractal_geometry_and_chaos::correlation_dimension;
use rssn::kernels::fractal_geometry_and_chaos::generate_burning_ship;
use rssn::kernels::fractal_geometry_and_chaos::generate_henon_map;
use rssn::kernels::fractal_geometry_and_chaos::generate_ifs_fractal;
use rssn::kernels::fractal_geometry_and_chaos::generate_julia_set;
use rssn::kernels::fractal_geometry_and_chaos::generate_lorenz_attractor;
use rssn::kernels::fractal_geometry_and_chaos::generate_mandelbrot_set;
use rssn::kernels::fractal_geometry_and_chaos::generate_multibrot;
use rssn::kernels::fractal_geometry_and_chaos::generate_rossler_attractor;
use rssn::kernels::fractal_geometry_and_chaos::logistic_bifurcation;
use rssn::kernels::fractal_geometry_and_chaos::logistic_map_iterate;
use rssn::kernels::fractal_geometry_and_chaos::lyapunov_exponent_logistic;
use rssn::kernels::fractal_geometry_and_chaos::orbit_density;
use rssn::kernels::fractal_geometry_and_chaos::orbit_entropy;
use rssn::kernels::fractal_geometry_and_chaos::sierpinski_triangle_ifs;

#[test]
fn fractal_escape_sets_and_data_accessors() {
    let mut fd = FractalData::new(4, 3, 50);
    fd.set(2, 1, 17);
    assert_eq!(fd.get(2, 1), Some(17));
    assert_eq!(fd.get(4, 1), None);

    // Mandelbrot: origin is bounded (max_iter), far corner escapes immediately.
    let m = generate_mandelbrot_set(9, 9, (-2.0, 2.0), (-2.0, 2.0), 40);
    assert_eq!(m.len(), 9);
    assert_eq!(m[4][4], 40);
    assert!(m[0][0] < 5);

    let j = generate_julia_set(7, 7, (-1.5, 1.5), (-1.5, 1.5), (0.0, 0.0), 30);
    assert_eq!(j[3][3], 30);

    let bs = generate_burning_ship(5, 5, (-2.0, 1.0), (-2.0, 1.0), 25);
    assert_eq!(bs.len(), 5);
    assert_eq!(bs[0].len(), 5);

    let mb = generate_multibrot(5, 5, (-1.5, 1.5), (-1.5, 1.5), 3.0, 25);
    assert_eq!(mb[2][2], 25);
}

#[test]
fn attractors_chaotic_maps_and_lyapunov_exponents() {
    let lorenz = generate_lorenz_attractor((1.0, 1.0, 1.0), 0.01, 200);
    assert_eq!(lorenz.len(), 200);
    assert!(
        lorenz
            .iter()
            .all(|&(x, y, z)| x.is_finite() && y.is_finite() && z.is_finite())
    );

    let rossler = generate_rossler_attractor((1.0, 0.0, 0.0), 0.01, 200, 0.2, 0.2, 5.7);
    assert_eq!(rossler.len(), 200);
    assert!(
        rossler
            .iter()
            .all(|&(x, y, z)| x.is_finite() && y.is_finite() && z.is_finite())
    );

    let henon = generate_henon_map((0.1, 0.1), 100, 1.4, 0.3);
    assert_eq!(henon.len(), 100);
    assert!(henon.iter().all(|&(x, y)| x.is_finite() && y.is_finite()));

    // Logistic map fixed point at r = 2.0 is x* = 0.5.
    let fixed = logistic_map_iterate(0.2, 2.0, 50);
    assert!((fixed.last().copied().unwrap() - 0.5).abs() < 1e-10);

    let bif = logistic_bifurcation((2.5, 3.8), 10, 50, 10, 0.3);
    assert_eq!(bif.len(), 100);

    // Stable regime has negative Lyapunov exponent; r = 4.0 has ~ ln(2).
    let lam_stable = lyapunov_exponent_logistic(2.5, 0.3, 200, 500);
    let lam_chaos = lyapunov_exponent_logistic(4.0, 0.314_159, 200, 1500);
    assert!(lam_stable < -0.1, "lam_stable = {lam_stable}");
    assert!(
        (lam_chaos - std::f64::consts::LN_2).abs() < 0.1,
        "lam_chaos = {lam_chaos}"
    );
}

#[test]
fn dimension_estimation_density_and_ifs() {
    // 1D line segment in 2D has box-counting and correlation dimension near 1.0.
    let line: Vec<(f64, f64)> = (0..500)
        .map(|i| {
            let t = i as f64 / 499.0;
            (t, 0.5 * t)
        })
        .collect();
    let box_dim = box_counting_dimension(&line, 6);
    let corr_dim = correlation_dimension(&line, 6);
    assert!((box_dim - 1.0).abs() < 0.25, "box_dim = {box_dim}");
    assert!((corr_dim - 1.0).abs() < 0.25, "corr_dim = {corr_dim}");

    let density = orbit_density(&line, 10, 10, (0.0, 1.0), (0.0, 0.5));
    let entropy = orbit_entropy(&density);
    assert!(entropy > 0.0);

    let (s_t, s_p) = sierpinski_triangle_ifs();
    let ifs_pts = generate_ifs_fractal(&s_t, &s_p, (0.0, 0.0), 300, 50);
    assert_eq!(ifs_pts.len(), 300);
    let (b_t, b_p) = barnsley_fern_ifs();
    assert_eq!(b_t.len(), b_p.len());
}
