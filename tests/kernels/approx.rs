//! Tests for the interpolation / approximation kernels.

use rssn::kernels::approx::ChebSeries;
use rssn::kernels::approx::CubicSpline;
use rssn::kernels::approx::Pchip;
use rssn::kernels::approx::SplineBc;
use rssn::kernels::approx::aaa;
use rssn::kernels::approx::bspline_basis;
use rssn::kernels::approx::bspline_eval;
use rssn::kernels::approx::bspline_interpolate;
use rssn::kernels::approx::chebyshev_points;
use rssn::kernels::approx::clenshaw;
use rssn::kernels::approx::horner;
use rssn::kernels::approx::pade;

#[test]
fn cubic_spline_variants() {
    let xs: Vec<f64> = (0..=10).map(|i| f64::from(i) * 0.3).collect();
    let ys: Vec<f64> = xs.iter().map(|x| x.sin()).collect();
    let nat = CubicSpline::new(&xs, &ys, SplineBc::Natural).unwrap();
    let clamped =
        CubicSpline::new(&xs, &ys, SplineBc::Clamped(0.0_f64.cos(), 3.0_f64.cos())).unwrap();
    let nak = CubicSpline::new(&xs, &ys, SplineBc::NotAKnot).unwrap();
    for x in [0.5, 1.4, 2.2] {
        assert!((clamped.eval(x) - x.sin()).abs() < 1e-4);
        assert!((nak.eval(x) - x.sin()).abs() < 1e-4);
        assert!((nat.eval(x) - x.sin()).abs() < 1e-2);
    }
    // interpolation property and derivative
    assert!((clamped.eval(xs[4]) - ys[4]).abs() < 1e-14);
    assert!((clamped.derivative(1.5, 1) - 1.5_f64.cos()).abs() < 1e-3);
    assert!((clamped.integrate(0.0, 3.0) - (1.0 - 3.0_f64.cos())).abs() < 1e-4);
    // not-a-knot reproduces cubics exactly
    let cx: Vec<f64> = (0..6).map(f64::from).collect();
    let cy: Vec<f64> = cx.iter().map(|x| x * x * x - 2.0 * x).collect();
    let s = CubicSpline::new(&cx, &cy, SplineBc::NotAKnot).unwrap();
    assert!((s.eval(2.5) - (2.5_f64.powi(3) - 5.0)).abs() < 1e-11);
    assert!(CubicSpline::new(&[0.0, 1.0], &[0.0, 1.0], SplineBc::Natural).is_err());
}

#[test]
fn cubic_spline_large_irregular_grid_and_boundary_conditions() {
    let xs = [0.0, 0.1, 0.7, 1.8, 3.0];
    let f = |x: f64| x.powi(3) - 2.0 * x + 1.0;
    let df = |x: f64| 3.0 * x * x - 2.0;
    let ys: Vec<_> = xs.iter().map(|&x| f(x)).collect();
    let clamped = CubicSpline::new(&xs, &ys, SplineBc::Clamped(df(0.0), df(3.0))).unwrap();
    for x in [0.05, 0.4, 1.3, 2.9] {
        assert!((clamped.eval(x) - f(x)).abs() < 1e-12);
        assert!((clamped.derivative(x, 1) - df(x)).abs() < 1e-11);
    }
    let natural = CubicSpline::new(&xs, &ys, SplineBc::Natural).unwrap();
    assert!(natural.derivative(xs[0], 2).abs() < 1e-12);
    assert!(natural.derivative(xs[4], 2).abs() < 1e-12);

    // This problem requires only linear storage; a dense matrix for the same
    // grid would require tens of gigabytes.
    let xs: Vec<_> = (0..50_000).map(|i| f64::from(i) * 0.01).collect();
    let ys: Vec<_> = xs.iter().map(|x| 2.0 * x + 1.0).collect();
    for bc in [SplineBc::Natural, SplineBc::Clamped(2.0, 2.0)] {
        let spline = CubicSpline::new(&xs, &ys, bc).unwrap();
        assert!((spline.eval(250.125) - 501.25).abs() < 1e-10);
    }
    assert!(CubicSpline::new(&[0.0, 1.0, 2.0], &[0.0, f64::NAN, 2.0], SplineBc::Natural).is_err());
    assert!(
        CubicSpline::new(
            &[0.0, 1.0, 2.0],
            &[0.0, 1.0, 2.0],
            SplineBc::Clamped(f64::NAN, 1.0)
        )
        .is_err()
    );
}

#[test]
fn pchip_is_monotone() {
    let x = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
    let y = [0.0, 0.0, 0.1, 5.0, 5.0, 5.1];
    let p = Pchip::new(&x, &y).unwrap();
    let mut prev = f64::NEG_INFINITY;
    for i in 0..=500 {
        let v = p.eval(f64::from(i) * 0.01);
        assert!(v >= prev - 1e-12);
        assert!((-1e-12..=5.1 + 1e-12).contains(&v));
        prev = v;
    }
    assert!((p.eval(2.0) - 0.1).abs() < 1e-14);
    assert!(p.derivative(0.5) >= 0.0);
}

#[test]
fn bsplines() {
    // partition of unity
    let knots = [0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0, 3.0];
    for x in [0.0, 0.4, 1.5, 2.9, 3.0] {
        let s: f64 = (0..6).map(|i| bspline_basis(&knots, 3, i, x)).sum();
        assert!((s - 1.0).abs() < 1e-14, "x={x}");
    }
    let coeffs = [1.0, 2.0, 0.5, -1.0, 3.0, 2.0];
    for x in [0.3, 1.2, 2.5] {
        let direct: f64 = (0..6)
            .map(|i| coeffs[i] * bspline_basis(&knots, 3, i, x))
            .sum();
        assert!((bspline_eval(&knots, &coeffs, 3, x) - direct).abs() < 1e-13);
    }
    let xs: Vec<f64> = (0..9).map(|i| f64::from(i) * 0.5).collect();
    let ys: Vec<f64> = xs.iter().map(|x| x.cos()).collect();
    let (k, c) = bspline_interpolate(&xs, &ys, 3).unwrap();
    for (x, y) in xs.iter().zip(&ys) {
        assert!((bspline_eval(&k, &c, 3, *x) - y).abs() < 1e-12);
    }
    assert!((bspline_eval(&k, &c, 3, 1.25) - 1.25_f64.cos()).abs() < 1e-3);
}

#[test]
fn chebyshev_approximation() {
    let s = ChebSeries::from_fn(f64::exp, -1.0, 1.0, 16);
    for x in [-0.9, 0.1, 0.77] {
        assert!((s.eval(x) - x.exp()).abs() < 1e-14);
    }
    // known coefficient: c0 of e^x on [-1,1] is I0(1)
    assert!((s.coeffs[0] - 1.266_065_877_752_008_4).abs() < 1e-14);
    let d = s.derivative();
    assert!((d.eval(0.3) - 0.3_f64.exp()).abs() < 1e-12);
    let i = s.integral();
    assert!((i.eval(1.0) - (1.0_f64.exp() - (-1.0_f64).exp())).abs() < 1e-13);
    assert!((s.definite_integral() - (1.0_f64.exp() - (-1.0_f64).exp())).abs() < 1e-13);
    // adaptive on a shifted interval chops the tail
    let a = ChebSeries::adaptive(f64::sin, 0.0, 6.0, 1e-14, 256);
    assert!(a.degree() < 40);
    assert!((a.eval(2.2) - 2.2_f64.sin()).abs() < 1e-13);
    // roots of sin on [0, 10]
    let r = ChebSeries::adaptive(f64::sin, 0.0, 10.0, 1e-14, 256)
        .roots()
        .unwrap();
    let expect = [
        0.0,
        std::f64::consts::PI,
        2.0 * std::f64::consts::PI,
        3.0 * std::f64::consts::PI,
    ];
    assert_eq!(r.len(), 4, "{r:?}");
    for (v, e) in r.iter().zip(expect) {
        assert!((v - e).abs() < 1e-9);
    }
    // Clenshaw directly: T_3(x) = 4x^3 - 3x
    assert!((clenshaw(&[0.0, 0.0, 0.0, 1.0], 0.5) - (4.0 * 0.125 - 1.5)).abs() < 1e-15);
    let p = chebyshev_points(4, -1.0, 1.0);
    assert!((p[0] - 1.0).abs() < 1e-15 && (p[4] + 1.0).abs() < 1e-15 && p[2].abs() < 1e-15);
}

#[test]
fn chebyshev_samples_reconstruct_known_chebyshev_modes() {
    for degree in [7_usize, 9, 17, 64, 257] {
        let samples: Vec<_> = (0..=degree)
            .map(|k| {
                let angle = std::f64::consts::PI * (k as f64 / degree as f64);
                let nyquist = if k.is_multiple_of(2) {
                    0.125
                } else {
                    -0.125
                };
                0.75 - 0.2 * angle.cos() + 0.18 * (3.0 * angle).cos() + nyquist
            })
            .collect();
        let series = ChebSeries::from_values(&samples, 1.0, 4.0);
        assert_eq!(series.degree(), degree);
        for (j, &c) in series.coeffs.iter().enumerate() {
            let expected = match j {
                | 0 => 0.75,
                | 1 => -0.2,
                | 3 => 0.18,
                | _ if j == degree => 0.125,
                | _ => 0.0,
            };
            assert!(
                (c - expected).abs() < 3e-15,
                "degree {degree}, mode {j}: {c}"
            );
        }
    }
}

#[test]
fn adaptive_chebyshev_reuses_nested_and_shared_values_within_degree_cap() {
    use std::cell::Cell;
    for (cap, expected_calls) in [(2, 3), (18, 33), (64, 65)] {
        let calls = Cell::new(0);
        let target = |x: f64| (101.0 * x).sin() + (4.27 * x).cos();
        let adaptive = ChebSeries::adaptive(
            |x| {
                calls.set(calls.get() + 1);
                target(x)
            },
            -1.0,
            1.0,
            0.0,
            cap,
        );
        assert_eq!(adaptive.degree(), cap);
        assert_eq!(calls.get(), expected_calls);
        assert_eq!(adaptive, ChebSeries::from_fn(target, -1.0, 1.0, cap));
    }
}

#[test]
fn chebyshev_points_and_evaluation_handle_degenerate_and_extreme_intervals() {
    assert_eq!(chebyshev_points(0, 2.0, 4.0), vec![3.0]);
    assert_eq!(chebyshev_points(4, 2.0, 2.0), vec![2.0; 5]);
    assert_eq!(ChebSeries::from_fn(|x| x, 2.0, 4.0, 0).eval(2.0), 3.0);
    assert_eq!(ChebSeries::adaptive(|x| x, 2.0, 4.0, 0.0, 0).eval(4.0), 3.0);
    assert_eq!(ChebSeries::from_values(&[], 1.0, 2.0).eval(1.5), 0.0);
    assert_eq!(
        ChebSeries::from_values(&[4.0, 9.0], 3.0, 3.0).eval(3.0),
        4.0
    );
    let calls = std::cell::Cell::new(0);
    let degenerate = ChebSeries::from_fn(
        |x| {
            calls.set(calls.get() + 1);
            x * x
        },
        3.0,
        3.0,
        64,
    );
    assert_eq!(calls.get(), 1);
    assert_eq!(degenerate.eval(3.0), 9.0);

    let max = f64::MAX;
    let large_constant = ChebSeries::from_values(&[max], 3.0, 3.0);
    assert_eq!(large_constant.derivative().eval(3.0), 0.0);
    assert_eq!(large_constant.integral().eval(3.0), 0.0);
    assert_eq!(large_constant.definite_integral(), 0.0);
    let points = chebyshev_points(4, -max, max);
    assert!(points.iter().all(|x| x.is_finite()));
    assert_eq!((points[0], points[2], points[4]), (max, 0.0, -max));
    let linear = ChebSeries {
        coeffs: vec![0.0, 1.0],
        a: -max,
        b: max,
    };
    assert_eq!(
        [linear.eval(-max), linear.eval(0.0), linear.eval(max)],
        [-1.0, 0.0, 1.0]
    );
    assert_eq!(linear.eval(0.5 * max), 0.5);

    let a = max.next_down().next_down().next_down();
    let narrow = ChebSeries { a, b: max, ..linear };
    assert_eq!(narrow.eval(a.next_up()), -1.0 / 3.0);
    let least = f64::from_bits(1);
    let tiny = ChebSeries {
        a: -2.0 * least,
        b: least,
        ..narrow
    };
    assert_eq!(tiny.eval(0.0), 1.0 / 3.0);
}

#[test]
fn chebyshev_transform_and_clenshaw_keep_extreme_finite_samples() {
    let max = f64::MAX;
    let large = ChebSeries::from_fn(|x| x, -0.5 * max, 0.5 * max, 64);
    assert!((large.eval(-0.5 * max) / max + 0.5).abs() < 2e-15);
    assert!((large.eval(0.5 * max) / max - 0.5).abs() < 2e-15);
    assert_eq!(
        ChebSeries::from_values(&[max; 20], -1.0, 1.0).eval(1.0),
        max
    );
    assert_eq!(clenshaw(&[0.0, 0.0, max], 1.0), max);

    let tiny = 2.0 * f64::from_bits(1);
    let samples: Vec<_> = (0..=17)
        .map(|k| if k % 2 == 0 { tiny } else { -tiny })
        .collect();
    let series = ChebSeries::from_values(&samples, -1.0, 1.0);
    assert_eq!(series.coeffs[17], tiny);
    assert_eq!(series.eval(-1.0), -tiny);
    assert_eq!(series.eval(1.0), tiny);
}

#[test]
fn aaa_rational_approximation() {
    // recovers a rational function essentially exactly
    let z: Vec<f64> = (0..200)
        .map(|i| -1.0 + 2.0 * f64::from(i) / 199.0)
        .collect();
    let f: Vec<f64> = z.iter().map(|x| (1.0 + x) / (2.0 + x * x)).collect();
    let r = aaa(&z, &f, 1e-13, 20).unwrap();
    assert!(r.degree() <= 4, "degree {}", r.degree());
    for x in [-0.95, 0.123, 0.8] {
        assert!((r.eval(x) - (1.0 + x) / (2.0 + x * x)).abs() < 1e-11);
    }
    // approximates |x| much better than a polynomial of the same degree could
    let f: Vec<f64> = z.iter().map(|x| x.abs()).collect();
    let r = aaa(&z, &f, 1e-6, 40).unwrap();
    assert!((r.eval(0.5) - 0.5).abs() < 1e-5);
    assert!(aaa(&[0.0, 1.0], &[0.0, 1.0], 1e-8, 5).is_err());
}

#[test]
fn pade_approximants() {
    // [2/2] of exp: (1 + x/2 + x^2/12) / (1 - x/2 + x^2/12)
    let c = [1.0, 1.0, 0.5, 1.0 / 6.0, 1.0 / 24.0];
    let (p, q) = pade(&c, 2, 2).unwrap();
    let expect_p = [1.0, 0.5, 1.0 / 12.0];
    let expect_q = [1.0, -0.5, 1.0 / 12.0];
    for i in 0..3 {
        assert!((p[i] - expect_p[i]).abs() < 1e-14);
        assert!((q[i] - expect_q[i]).abs() < 1e-14);
    }
    let x = 0.5;
    assert!((horner(&p, x) / horner(&q, x) - x.exp()).abs() < 1e-4);
    assert!(pade(&c, 3, 3).is_err());
}
