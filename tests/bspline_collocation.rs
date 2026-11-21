//! Shared B-spline basis, curve evaluation and collocation.

use rssn::kernels::approx::bspline_basis;
use rssn::kernels::approx::bspline_eval;
use rssn::kernels::approx::bspline_interpolate;
use rssn::kernels::interp_ext;

fn close(
    actual: f64,
    expected: f64,
    tolerance: f64,
) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} vs {expected}"
    );
}

#[test]
fn degree_zero_interpolates_each_sample_and_supports_singletons() {
    let (knots, coefficients) = bspline_interpolate(&[0.0, 1.0, 3.0], &[2.0, 4.0, 7.0], 0).unwrap();
    assert_eq!(knots, [0.0, 0.5, 2.0, 3.0]);
    for (x, expected) in [(-2.0, 2.0), (0.5, 4.0), (1.99, 4.0), (2.0, 7.0), (4.0, 7.0)] {
        close(bspline_eval(&knots, &coefficients, 0, x), expected, 0.0);
    }
    let singleton = interp_ext::BSpline::interpolate(&[4.0], &[7.0], 0).unwrap();
    close(singleton.eval(50.0).unwrap(), 7.0, 0.0);

    let next = 1.0_f64.next_up();
    let (knots, coefficients) = bspline_interpolate(&[1.0, next], &[2.0, 4.0], 0).unwrap();
    close(bspline_eval(&knots, &coefficients, 0, 1.0), 2.0, 0.0);
    close(bspline_eval(&knots, &coefficients, 0, next), 4.0, 0.0);
}

#[test]
fn basis_partition_and_curve_agree_at_repeated_knots() {
    for degree in 0..=8 {
        let mut knots = vec![0.0; degree + 1];
        knots.extend([0.25, 0.5, 0.5, 0.75]);
        knots.extend(std::iter::repeat_n(1.0, degree + 1));
        let count = knots.len() - degree - 1;
        let coefficients: Vec<_> = (0..count).map(|i| ((i as f64) * 1.7).sin()).collect();
        for x in [
            0.0,
            0.25,
            0.5_f64.next_down(),
            0.5,
            0.5_f64.next_up(),
            0.75,
            1.0,
        ] {
            let weights: Vec<_> = (0..count)
                .map(|i| bspline_basis(&knots, degree, i, x))
                .collect();
            assert!(weights.iter().all(|weight| (0.0..=1.0).contains(weight)));
            close(weights.iter().sum(), 1.0, 4e-15);
            let expected: f64 = coefficients.iter().zip(weights).map(|(c, w)| c * w).sum();
            let raw = bspline_eval(&knots, &coefficients, degree, x);
            let checked = interp_ext::bspline_eval(&knots, &coefficients, degree, x).unwrap();
            assert_eq!(raw.to_bits(), checked.to_bits());
            close(raw, expected, 4e-15);
        }
    }

    let excess = [0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
    let coefficients = [2.0, 4.0, 300.0, 400.0];
    close(bspline_eval(&excess, &coefficients, 1, 1.0), 4.0, 0.0);
    close(
        interp_ext::bspline_eval(&excess, &coefficients, 1, 1.0).unwrap(),
        4.0,
        0.0,
    );
}

#[test]
fn clamped_spline_fit_recovers_polynomials_in_both_solver_regimes() {
    let cubic = |x: f64| 2.0 * x.powi(3) - 3.0 * x * x + 0.5 * x + 1.0;
    for (count, degree) in [(8, 3), (153, 3), (330, 7)] {
        let x: Vec<_> = (0..count)
            .map(|i| (f64::from(i) / f64::from(count - 1)).powf(1.15))
            .collect();
        let y: Vec<_> = x.iter().copied().map(cubic).collect();
        let (knots, coeffs) = bspline_interpolate(&x, &y, degree).unwrap();
        let checked = interp_ext::BSpline::interpolate(&x, &y, degree).unwrap();
        assert_eq!(knots, checked.knots);
        assert_eq!(coeffs, checked.coeffs);
        for i in 0..=100 {
            let point = f64::from(i) / 100.0;
            close(
                bspline_eval(&knots, &coeffs, degree, point),
                cubic(point),
                3e-11,
            );
        }
    }
}

#[test]
fn spline_intermediate_distances_allow_finite_extremes() {
    let knots = [-1e308, -1e308, 1e308, 1e308];
    close(bspline_basis(&knots, 1, 0, 0.0), 0.5, 0.0);
    close(bspline_eval(&knots, &[2.0, 4.0], 1, 0.0), 3.0, 0.0);

    let x = [1e308, 1.1e308, 1.2e308, 1.3e308, 1.4e308];
    let y = x.map(|coordinate| (coordinate / 1e308_f64).powi(2));
    let (knots, coeffs) = bspline_interpolate(&x, &y, 2).unwrap();
    assert!(knots.iter().all(|v| v.is_finite()));
    for (&x, &y) in x.iter().zip(&y) {
        close(bspline_eval(&knots, &coeffs, 2, x), y, 3e-14);
    }

    let x = [0.0, 1e-310, 2e-310, 3e-310];
    let y = x.map(|coordinate| (coordinate / 1e-310_f64).powi(2));
    let (knots, coeffs) = bspline_interpolate(&x, &y, 3).unwrap();
    for (&x, &y) in x.iter().zip(&y) {
        close(bspline_eval(&knots, &coeffs, 3, x), y, 3e-14);
    }
}

#[test]
fn spline_shapes_and_nonfinite_controls_are_rejected() {
    assert!(bspline_basis(&[0.0, 1.0], 0, usize::MAX, 0.5).is_nan());
    assert!(bspline_basis(&[0.0, f64::NAN], 0, 0, 0.5).is_nan());
    assert!(bspline_eval(&[], &[], 0, 0.0).is_nan());
    assert!(bspline_eval(&[1.0, 0.0], &[2.0], 0, 0.0).is_nan());
    assert!(bspline_interpolate(&[0.0], &[1.0], usize::MAX).is_err());
    assert!(bspline_interpolate(&[0.0, f64::NAN], &[1.0, 2.0], 1).is_err());
    assert!(bspline_interpolate(&[0.0, 1.0], &[1.0, f64::INFINITY], 1).is_err());
    assert!(bspline_interpolate(&[f64::MAX.next_down(), f64::MAX], &[1.0, 2.0], 0).is_err());
    assert!(interp_ext::bspline_eval(&[0.0, 1.0], &[2.0], 0, f64::NAN).is_err());
}

#[test]
fn collocation_rows_preserve_values_across_worker_counts() {
    let x: Vec<_> = (0..2049).map(|i| f64::from(i) / 2048.0).collect();
    let y: Vec<_> = x.iter().map(|x| (10.0 * x).cos()).collect();
    let run = |threads| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| bspline_interpolate(&x, &y, 4).unwrap())
    };
    let (knots, coefficients) = run(1);
    let (parallel_knots, parallel_coefficients) = run(4);
    assert_eq!(knots, parallel_knots);
    let bits = |values: Vec<f64>| values.into_iter().map(f64::to_bits).collect::<Vec<_>>();
    assert_eq!(bits(coefficients), bits(parallel_coefficients));
}
