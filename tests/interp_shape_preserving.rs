//! PCHIP and Akima interpolation with shared local segments.

use rssn::kernels::interp_ext::InterpError;
use rssn::kernels::interp_ext::akima;
use rssn::kernels::interp_ext::akima_many;
use rssn::kernels::interp_ext::pchip;
use rssn::kernels::interp_ext::pchip_many;

type Scalar = fn(&[f64], &[f64], f64) -> Result<f64, InterpError>;
type Batch = fn(&[f64], &[f64], &[f64]) -> Result<Vec<f64>, InterpError>;

const METHODS: [(Scalar, Batch); 2] = [(pchip, pchip_many), (akima, akima_many)];

#[test]
fn pchip_uses_nonuniform_secant_weights_and_preserves_local_ranges() {
    let x = [0.0, 1.0, 3.0];
    let y = [0.0, 1.0, 5.0];
    // The first and second knot slopes are 2/3 and 9/7, respectively.
    let expected = 0.5 + (2.0 / 3.0 - 9.0 / 7.0) / 8.0;
    assert!((pchip(&x, &y, 0.5).unwrap() - expected).abs() < 5e-16);

    let x = [0.0, 0.2, 1.0, 1.01, 2.0, 4.0];
    let y = [3.0, 2.0, 2.0, 2.001, 5.0, -1.0];
    for (x_seg, y_seg) in x.windows(2).zip(y.windows(2)) {
        let lower = f64::min(y_seg[0], y_seg[1]);
        let upper = f64::max(y_seg[0], y_seg[1]);
        let queries: Vec<_> = (0..=100)
            .map(|i| x_seg[0] + (x_seg[1] - x_seg[0]) * f64::from(i) * 0.01)
            .collect();
        for value in pchip_many(&x, &y, &queries).unwrap() {
            assert!((lower - 1e-14..=upper + 1e-14).contains(&value));
        }
    }
}

#[test]
fn akima_reproduces_a_quadratic_including_edge_segments() {
    let knots = [0.0, 1.0, 2.0, 3.0, 4.0];
    let values = knots.map(|x| x * x);
    for query in [-1.25, 0.125, 0.625, 1.0, 1.5, 2.5, 3.5, 4.0, 4.75] {
        let actual = akima(&knots, &values, query).unwrap();
        assert!((actual - query * query).abs() < 4e-14);
    }
}

#[test]
fn akima_matches_exact_reference_values_for_nonuniform_knots() {
    let x = [0.0, 0.125, 2.0, 2.25, 4.0, 9.0, 10.0];
    let y = [-2.0, 1.0, 1.5, 0.0, 12.0, 4.0, -1.0];
    let q = [-0.125, 0.03125, 1.0, 2.125, 3.0, 5.0, 9.5, 10.25];
    // Independent rational evaluations of the Akima Hermite segments.
    let exact = [
        -21_058.0 / 3_375.0,
        -268_241.0 / 288_000.0,
        9_780_780_589.0 / 2_917_518_750.0,
        1_923_897.0 / 2_969_866.0,
        688_360_068.0 / 150_864_091.0,
        351_840_812.0 / 29_516_875.0,
        60_901.0 / 33_200.0,
        -150_749.0 / 53_120.0,
    ];
    for (actual, expected) in akima_many(&x, &y, &q).unwrap().into_iter().zip(exact) {
        assert!((actual - expected).abs() < 4e-14);
    }
}

#[test]
fn adjoining_segments_have_matching_first_derivatives_at_knots() {
    let knots = [0.0, 0.5, 1.25, 1.3125, 3.0, 6.0];
    let values = [2.0, 5.0, 1.0, -4.0, 6.0, -12.0];
    for (scalar, _) in METHODS {
        for (&x, &y) in knots[1..knots.len() - 1].iter().zip(&values[1..]) {
            let distance = 1e-7;
            let before = scalar(&knots, &values, x - distance).unwrap();
            let after = scalar(&knots, &values, x + distance).unwrap();
            let left = (y - before) / distance;
            let right = (after - y) / distance;
            assert!((left - right).abs() < 5e-4);
        }
    }
}

#[test]
fn batch_matches_scalar_for_unsorted_queries_and_duplicate_queries_at_knots() {
    let x = [-8.0, -2.0, -0.25, 0.0, 0.5, 0.501, 2.0, 4.0, 7.0];
    let y = [0.0, 2.0, -3.0, 2.0, 4.0, 1.0, -3.0, 1.0, 5.0];
    let queries = [-10.0, 0.2, 1.0, 7.0, 0.501, 0.501, -2.0, 0.0, 8.0];
    for (scalar, batch) in METHODS {
        let result = batch(&x, &y, &queries).unwrap();
        let expected: Vec<_> = queries
            .iter()
            .map(|q| scalar(&x, &y, *q).unwrap())
            .collect();
        assert_eq!(
            result.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
            expected.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
        );
        assert_eq!(batch(&x, &y, &x).unwrap(), y);
    }
}

#[test]
fn one_or_two_knots_represent_a_constant_or_line() {
    for (scalar, batch) in METHODS {
        assert_eq!(
            batch(&[1.0], &[3.0], &[-f64::MAX, 1.0, f64::MAX]).unwrap(),
            [3.0; 3]
        );
        assert_eq!(scalar(&[-2.0, 2.0], &[-5.0, 7.0], 3.0).unwrap(), 10.0);
        assert_eq!(
            scalar(&[-f64::MAX, f64::MAX], &[-f64::MAX, f64::MAX], 0.0).unwrap(),
            0.0
        );
        assert_eq!(
            scalar(&[0.0, 1.0, 2.0], &[f64::MAX; 3], -f64::MAX).unwrap(),
            f64::MAX
        );
        let distant = scalar(&[0.0, 1e-300], &[0.0, 1e-300], 1e308).unwrap();
        assert!((distant / 1e308 - 1.0).abs() < 4e-16);
        let adjacent = scalar(&[0.0, f64::MAX], &[0.0, f64::MAX], f64::from_bits(1)).unwrap();
        assert_eq!(adjacent.to_bits(), 1);
    }
}

#[test]
fn affine_rescaling_keeps_values_near_both_ends_of_float_range() {
    for (scalar, _) in METHODS {
        let at_midpoint = scalar(
            &[-f64::MAX, -f64::MAX / 2.0, f64::MAX / 2.0, f64::MAX],
            &[-2.0, -1.0, 1.0, 2.0],
            0.0,
        )
        .unwrap();
        assert!(at_midpoint.abs() < 2e-15);
        let big = scalar(&[0.0, 1.0, 2.0], &[-f64::MAX, 0.0, f64::MAX], 0.5).unwrap();
        assert!((big / f64::MAX + 0.5).abs() < 1e-15);

        let step = f64::from_bits(128);
        let tiny = scalar(
            &[0.0, step, 2.0 * step],
            &[0.0, step, 2.0 * step],
            step * 1.25,
        )
        .unwrap();
        assert_eq!(tiny.to_bits(), 160);

        let small = scalar(&[1e-300, 2e-300, 3e-300], &[1e300, 2e300, 3e300], 2.5e-300).unwrap();
        assert!((small / 1e300 - 2.5).abs() < 2e-15);
    }
}

#[test]
fn affine_rescaling_of_nonuniform_nonlinear_data_preserves_shape() {
    let knots = [0.0, 0.0625, 1.25, 1.5, 1.75, 7.0];
    let values = [3.0, 2.0, 2.25, -0.5, 11.0, 8.75];
    let points = [-0.08, 0.05, 0.13, 0.5, 1.35, 1.75, 3.125, 7.15];
    for (_, batch) in METHODS {
        let expected = batch(&knots, &values, &points).unwrap();
        for (dx, dy) in [(1e-300, 1e300), (1e300, 1e-300), (1e-160, 1e-160)] {
            let x = knots.map(|x| x * dx);
            let y = values.map(|y| y * dy);
            let queries = points.map(|q| q * dx);
            let actual = batch(&x, &y, &queries).unwrap();
            for (observed, want) in actual.iter().zip(&expected) {
                assert!((observed / dy - want).abs() <= 5e-14 * (1.0 + want.abs()));
            }
        }
    }
}

#[test]
fn validation_applies_to_scalar_and_batch_even_for_empty_queries() {
    for (scalar, batch) in METHODS {
        assert!(batch(&[0.0], &[2.0], &[]).unwrap().is_empty());
        assert_eq!(batch(&[], &[], &[]), Err(InterpError::Empty));
        assert!(matches!(
            batch(&[0.0], &[], &[]),
            Err(InterpError::ShapeMismatch(_))
        ));
        assert!(matches!(
            batch(&[-0.0, 0.0], &[1.0, 2.0], &[]),
            Err(InterpError::Domain(_))
        ));
        assert!(scalar(&[1.0, 0.0], &[1.0, 2.0], 0.5).is_err());
        assert!(batch(&[0.0], &[f64::INFINITY], &[]).is_err());
        assert!(batch(&[0.0], &[0.0], &[0.0, f64::NAN]).is_err());
        assert!(scalar(&[0.0], &[0.0], f64::INFINITY).is_err());
        assert!(scalar(&[0.0, 1.0], &[0.0, f64::MAX], 2.0).is_err());
    }
}

#[test]
fn worker_count_does_not_change_the_order_or_bits() {
    let points: Vec<_> = (0..251).map(f64::from).collect();
    let values: Vec<_> = points.iter().map(|x| (x * 0.125).sin()).collect();
    let queries: Vec<_> = (0..65_540)
        .map(|i| f64::from(i % 2_500) * 0.1 + 0.011)
        .collect();
    let serial = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap();
    let parallel = rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build()
        .unwrap();
    for (_, batch) in METHODS {
        let left = serial
            .install(|| batch(&points, &values, &queries))
            .unwrap();
        let right = parallel
            .install(|| batch(&points, &values, &queries))
            .unwrap();
        assert!(
            left.iter()
                .zip(right)
                .all(|(a, b)| a.to_bits() == b.to_bits())
        );
    }
}
