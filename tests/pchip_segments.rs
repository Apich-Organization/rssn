//! Reusable PCHIP segments, scale covariance and derivative continuity.

use rssn::kernels::approx::Pchip;
use rssn::kernels::interp_ext::pchip;
use rssn::kernels::interp_ext::pchip_many;

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
fn fitted_pchip_is_covariant_under_small_and_large_ordinate_scales() {
    let x = [0.0, 1.0, 2.0, 3.0];
    for scale in [1e-200, 1.0, 1e200] {
        let y = [0.0, 1.0, 4.0, 9.0].map(|value| value * scale);
        let model = Pchip::new(&x, &y).unwrap();
        close(model.eval(1.5) / scale, 2.21875, 2e-14);
        close(model.derivative(1.5) / scale, 3.1875, 2e-14);
        close(model.derivative(1.0_f64.next_down()) / scale, 1.5, 2e-14);
        close(model.derivative(1.0) / scale, 1.5, 2e-14);
        close(model.derivative(1.0_f64.next_up()) / scale, 1.5, 2e-14);
        assert_eq!(
            model.eval(1.5).to_bits(),
            pchip(&x, &y, 1.5).unwrap().to_bits()
        );
    }
}

#[test]
fn fitted_pchip_avoids_overflow_and_underflow_in_physical_slopes() {
    let extreme = [-1e308, 1e308];
    let line = Pchip::new(&extreme, &extreme).unwrap();
    close(line.eval(0.0), 0.0, 0.0);
    close(line.derivative(0.0), 1.0, 0.0);
    close(line.derivative(extreme[1]), 1.0, 0.0);

    let scale = 1e-310;
    let x = [0.0, 1.0, 2.0, 3.0].map(|value| value * scale);
    let y = [0.0, 1.0, 4.0, 9.0].map(|value| value * scale);
    let model = Pchip::new(&x, &y).unwrap();
    close(model.eval(1.5 * scale) / scale, 2.21875, 3e-12);
    close(model.derivative(1.5 * scale), 3.1875, 3e-12);
}

#[test]
fn prepared_and_single_query_paths_share_segment_boundaries() {
    let x = [-2.1, -0.83, -0.22, 0.78, 2.7, 4.0];
    let y = [0.8, 0.8, 1.0, 1.5, 1.6, 3.9];
    let model = Pchip::new(&x, &y).unwrap();
    let queries = [4.7, -2.1, -1.3, -0.83, -3.0, 0.78, 0.78, 2.7, 3.7, 4.0];
    let batch = pchip_many(&x, &y, &queries).unwrap();
    for (query, expected) in queries.into_iter().zip(batch) {
        close(model.eval(query), expected, 2e-14 * expected.abs().max(1.0));
    }
    for (point, expected) in x.into_iter().zip(y) {
        assert_eq!(model.eval(point).to_bits(), expected.to_bits());
    }
}

#[test]
fn fitted_pchip_rejects_invalid_samples_and_nonfinite_queries() {
    assert!(Pchip::new(&[0.0], &[1.0]).is_err());
    assert!(Pchip::new(&[0.0, 1.0], &[1.0]).is_err());
    assert!(Pchip::new(&[0.0, 0.0], &[1.0, 2.0]).is_err());
    assert!(Pchip::new(&[0.0, f64::NAN], &[1.0, 2.0]).is_err());
    assert!(Pchip::new(&[0.0, 1.0], &[1.0, f64::INFINITY]).is_err());
    let line = Pchip::new(&[0.0, 1.0], &[2.0, 4.0]).unwrap();
    assert!(line.eval(f64::NAN).is_nan());
    assert!(line.derivative(f64::INFINITY).is_nan());
}

#[test]
fn prepared_segments_are_worker_independent() {
    let x: Vec<_> = (0..4600).map(|i| f64::from(i) / 30.0).collect();
    for scale in [1.0, 1e-200] {
        let y: Vec<_> = x.iter().map(|x| (x.sin() + x / 100.0) * scale).collect();
        let build = |threads| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(|| Pchip::new(&x, &y).unwrap())
        };
        let single = build(1);
        let parallel = build(4);
        for query in (0..=1000).map(|i| f64::from(i) * 0.151) {
            assert_eq!(single.eval(query).to_bits(), parallel.eval(query).to_bits());
            assert_eq!(
                single.derivative(query).to_bits(),
                parallel.derivative(query).to_bits()
            );
        }
    }
}
