//! Numerical and API coverage for complete-window rolling moments.

use std::sync::Arc;
use std::sync::atomic::AtomicBool;

use rssn::kernels::timeseries::rolling::RollingError;
use rssn::kernels::timeseries::rolling::RollingMoment;
use rssn::kernels::timeseries::rolling::rolling_moments;

fn relative_eq(
    actual: f64,
    expected: f64,
    tolerance: f64,
) {
    let error = (actual - expected).abs();
    assert!(
        error <= tolerance * expected.abs(),
        "{actual} != {expected} (error {error})"
    );
}

#[test]
fn full_windows_and_degrees_of_freedom() {
    let series = [1.0, 2.0, 3.0, 6.0];
    let population = rolling_moments(&series, 3, 0).unwrap();
    let sample = rolling_moments(&series, 3, 1).unwrap();
    assert_eq!(population.len(), 2);
    assert_eq!(sample.len(), 2);

    for (result, expected_mean, expected_variance) in [
        (population[0], 2.0, 2.0 / 3.0),
        (population[1], 11.0 / 3.0, 26.0 / 9.0),
        (sample[0], 2.0, 1.0),
        (sample[1], 11.0 / 3.0, 13.0 / 3.0),
    ] {
        relative_eq(result.mean, expected_mean, 1e-15);
        relative_eq(result.variance, expected_variance, 1e-15);
    }

    assert!(
        rolling_moments(&series, series.len() + 1, 0)
            .unwrap()
            .is_empty()
    );
    assert!(rolling_moments(&[], 1, 0).unwrap().is_empty());
    assert_eq!(
        rolling_moments(&series, 1, 0).unwrap(),
        series.map(|mean| RollingMoment { mean, variance: 0.0 }),
    );
    assert_eq!(rolling_moments(&series, series.len(), 0).unwrap().len(), 1);
}

#[test]
fn agrees_with_anchored_two_pass_reference_at_large_offsets() {
    let data: Vec<_> = (0..130)
        .map(|i| 1e16 + f64::from((i * 41 + i * i * 13) % 59 - 29) * 2.0)
        .collect();
    for width in 1..=data.len() {
        let ddof = width / 2;
        let roll = rolling_moments(&data, width, ddof).unwrap();
        assert_eq!(roll.len(), data.len() - width + 1);
        for (actual, window) in roll.iter().zip(data.windows(width)) {
            let anchor = window[0];
            let mean = window.iter().map(|x| x - anchor).sum::<f64>() / width as f64;
            let variance = window
                .iter()
                .map(|x| ((x - anchor) - mean).powi(2))
                .sum::<f64>()
                / (width - ddof) as f64;
            assert!((actual.mean - (anchor + mean)).abs() <= 2.0);
            relative_eq(actual.variance, variance, 3e-14);
        }
    }
}

#[test]
fn preserves_extreme_means_and_representable_variances() {
    for mean in [f64::MAX, -f64::MAX] {
        assert_eq!(
            rolling_moments(&[mean; 10], 4, 1).unwrap(),
            vec![RollingMoment { mean, variance: 0.0 }; 7],
        );
        let series = [mean; 8193];
        for width in [3, 63, 1023, 1024, 1025, 4095, 8193] {
            let expected = RollingMoment { mean, variance: 0.0 };
            assert!(
                rolling_moments(&series, width, 1)
                    .unwrap()
                    .iter()
                    .all(|result| *result == expected)
            );
        }
    }

    let data = [0.0, 1.6e154, 0.0, 1.6e154, 0.0, 1.6e154];
    for result in rolling_moments(&data, 4, 0).unwrap() {
        relative_eq(result.mean, 8e153, 1e-15);
        relative_eq(result.variance, 6.4e307, 1e-15);
    }

    let amplitude = 1e-160;
    let tiny = rolling_moments(&[-amplitude, amplitude, -amplitude, amplitude], 4, 0).unwrap();
    let expected = amplitude * amplitude;
    assert!(tiny[0].variance.is_subnormal());
    assert!(tiny[0].variance.to_bits().abs_diff(expected.to_bits()) <= 1);
}

#[test]
fn preserves_small_mean_residuals_after_mixed_sign_cancellation() {
    for (large, residual) in [(1e100, 1.0), (1e150, 1e-160), (1e150, -1e-160)] {
        for series in [
            [large, residual, -large],
            [large, -large, residual],
            [residual, large, -large],
        ] {
            let result = rolling_moments(&series, 3, 0).unwrap()[0];
            relative_eq(result.mean, residual / 3.0, 1e-15);
            relative_eq(result.variance, (large * (2.0 / 3.0)) * large, 1e-15);
        }
    }
    let result = rolling_moments(&[1e100, 1e100, -2e100], 3, 0).unwrap()[0];
    assert_eq!(result.mean.to_bits(), 0.0_f64.to_bits());
    assert_eq!(
        rolling_moments(&[f64::MAX, f64::MAX, -f64::MAX], 3, 0),
        Err(RollingError::VarianceOverflow { start: 0 }),
    );
}

#[test]
fn reports_invalid_requests_and_first_bad_window() {
    assert_eq!(
        rolling_moments(&[0.0], 0, 0),
        Err(RollingError::InvalidWindow)
    );
    assert_eq!(
        rolling_moments(&[0.0], 1, 1),
        Err(RollingError::InvalidDegreesOfFreedom)
    );
    assert_eq!(
        rolling_moments(&[0.0, f64::NAN, f64::INFINITY], 10, 0),
        Err(RollingError::NonFiniteInput { index: 1 }),
    );
    assert_eq!(
        rolling_moments(&[0.0, 0.0, f64::MAX, -f64::MAX], 2, 0),
        Err(RollingError::VarianceOverflow { start: 1 }),
    );
}

#[test]
fn worker_count_preserves_values_and_first_error() {
    let mut data: Vec<_> = (0..65_537)
        .map(|i| f64::from((i * 421) % 1979 - 989) + 1e12)
        .collect();
    let pool = |threads| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
    };
    let serial = pool(1);
    let parallel = pool(4);
    for width in [16, 47] {
        let result =
            |pool: &rayon::ThreadPool| pool.install(|| rolling_moments(&data, width, 1).unwrap());
        let serial_moments = result(&serial);
        let parallel_moments = result(&parallel);
        assert_eq!(serial_moments.len(), data.len() - width + 1);
        assert_eq!(serial_moments.len(), parallel_moments.len());
        for (a, b) in serial_moments.iter().zip(parallel_moments) {
            assert_eq!(a.mean.to_bits(), b.mean.to_bits());
            assert_eq!(a.variance.to_bits(), b.variance.to_bits());
        }
    }

    data[41_000] = f64::MAX;
    data[62_000] = -f64::MAX;
    assert_eq!(
        parallel.install(|| rolling_moments(&data, 47, 1)),
        Err(RollingError::VarianceOverflow { start: 41_000 - 46 }),
    );
}

#[test]
fn honors_current_host_cancellation() {
    let token = Arc::new(AtomicBool::new(true));
    let _guard = rssn::cancellation::install(token);
    assert_eq!(
        rolling_moments(&[1.0, 2.0, 3.0], 2, 0),
        Err(RollingError::Cancelled)
    );
}

#[test]
fn short_windows_at_both_guard_exponents_and_signs() {
    for exp in [-400, 400] {
        let anchor = f64::from_bits(((1023 + exp) as u64) << 52);
        let neighbour = if exp < 0 {
            anchor.next_up()
        } else {
            anchor.next_down()
        };
        for sign in [-1.0, 1.0] {
            let [a, b] = [anchor * sign, neighbour * sign];
            let result = rolling_moments(&[a, b], 2, 0).unwrap()[0];
            assert_eq!(result.mean.to_bits(), a.midpoint(b).to_bits());
            assert_eq!(result.variance.to_bits(), ((a - b) * 0.5).powi(2).to_bits());
        }
    }
}
