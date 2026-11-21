//! Numerical and scheduling regressions for batch statistics.

use rssn::kernels::dense::Mat;
use rssn::kernels::descriptive;
use rssn::kernels::distributions::DistError;
use rssn::kernels::distributions::Distribution;
use rssn::kernels::distributions::multivariate_normal_log_pdf;
use rssn::kernels::distributions::multivariate_normal_log_pdfs;
use rssn::kernels::online_stats::OnlineCovariance;
use rssn::kernels::online_stats::OnlineStatsError;
use rssn::kernels::robust;
use rssn::kernels::stats_ext;

fn close(
    actual: f64,
    expected: f64,
) {
    assert!(
        (actual - expected).abs() <= 3e-13 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

fn pool(threads: usize) -> rayon::ThreadPool {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap()
}

#[test]
fn batch_quantiles_match_scalar_algorithms_and_handle_nonfinite_input() {
    let data = [10.0, -0.0, -5.0, 10.0, 0.0, 4.0];
    let probabilities = [0.75, 0.0, 1.0, 0.25, 0.5, 0.25];
    let expected = [8.5, -5.0, 10.0, 0.0, 2.0, 0.0];
    let batch = descriptive::quantiles(&data, &probabilities).unwrap();
    assert_eq!(batch, expected);
    for (&p, &value) in probabilities.iter().zip(&batch) {
        assert_eq!(
            descriptive::quantile(&data, p).unwrap().to_bits(),
            value.to_bits()
        );
    }
    assert_eq!(
        descriptive::quantile(&[0.0, -0.0], 0.0).unwrap().to_bits(),
        (-0.0_f64).to_bits()
    );
    assert_eq!(
        descriptive::quantile(&[-0.0, 0.0], 1.0).unwrap().to_bits(),
        0.0_f64.to_bits()
    );
    assert!(descriptive::quantiles(&data, &[f64::NAN]).is_err());
    assert!(descriptive::quantiles(&[f64::INFINITY], &[]).is_err());
    assert_eq!(
        descriptive::quantiles(&data, &[]).unwrap(),
        Vec::<f64>::new()
    );

    let with_nan = [10.0, -0.0, -5.0, 10.0, 0.0, 4.0, f64::NAN];
    for kind in 1..=9 {
        let batch = robust::quantiles(&with_nan, &probabilities, kind);
        for (&p, &value) in probabilities.iter().zip(&batch) {
            assert_eq!(robust::quantile(&data, p, kind).to_bits(), value.to_bits());
        }
    }
    let batch = robust::quantiles(&data, &[0.0, f64::NAN, -1.0, 1.0], 7);
    assert_eq!(batch[0], -5.0);
    assert!(batch[1].is_nan() && batch[2].is_nan());
    assert_eq!(batch[3], 10.0);
    assert!(
        robust::quantiles(&data, &probabilities, 0)
            .iter()
            .all(|p| p.is_nan())
    );
}

#[test]
fn weighted_quantile_and_ecdf_preserve_order_ties_and_endpoints() {
    let sample = [10.0, 2.0, 6.0, -50.0, f64::NAN, 2.0];
    let weights = [1.0, 0.5, 1.0, 0.0, 0.0, 0.5];
    let queries = [1.0, 0.0, 1.0 / 3.0, 0.5, f64::NAN];
    let q = robust::weighted_quantiles(&sample, &weights, &queries);
    assert_eq!(&q[..4], &[10.0, 2.0, 2.0, 6.0]);
    assert!(q[4].is_nan());
    assert!(
        robust::weighted_quantiles(&sample, &[0.0; 6], &queries)
            .iter()
            .all(|p| p.is_nan())
    );

    let cdf = robust::ecdf_batch(&sample, &[10.0, f64::NEG_INFINITY, 2.0, f64::NAN, 5.0]);
    assert_eq!(&cdf[..3], &[1.0, 0.0, 0.6]);
    assert!(cdf[3].is_nan());
    assert_eq!(cdf[4], 0.6);
    assert!(
        robust::ecdf_batch(&[f64::NAN], &[1.0, 2.0])
            .iter()
            .all(|p| p.is_nan())
    );
}

#[test]
fn adjusted_p_values_control_the_expected_error_rates() {
    let probabilities = [0.04, 0.01, 0.03, 0.04, 1.0, 0.0];
    let cases = [
        (
            stats_ext::adjust_p_bh as fn(&[f64]) -> _,
            [0.048, 0.03, 0.048, 0.048, 1.0, 0.0],
        ),
        (
            stats_ext::adjust_p_by,
            [0.1176, 0.0735, 0.1176, 0.1176, 1.0, 0.0],
        ),
        (stats_ext::adjust_p_holm, [0.12, 0.05, 0.12, 0.12, 1.0, 0.0]),
    ];
    for (adjust, expected) in cases {
        for (actual, expected) in adjust(&probabilities).unwrap().into_iter().zip(expected) {
            close(actual, expected);
        }
        assert!(adjust(&[]).unwrap().is_empty());
        assert_eq!(adjust(&[0.03]).unwrap(), vec![0.03]);
        for invalid in [f64::NAN, f64::INFINITY, -0.01, 1.01] {
            assert!(adjust(&[0.5, invalid]).is_err());
        }
    }
}

#[test]
fn online_covariance_batch_is_transactional_and_independent_of_worker_count() {
    let samples: Vec<f64> = (0..5_000)
        .flat_map(|i| {
            let x = f64::from(i % 17);
            [x, 2.0 * x, -4.0 * x]
        })
        .collect();
    let prepare = || {
        let mut covariance = OnlineCovariance::new(3).unwrap();
        covariance.push_batch(&samples).unwrap();
        covariance
    };
    let sequential = pool(1).install(prepare);
    let mut parallel = pool(4).install(prepare);
    assert_eq!(parallel.count(), 5_000);
    assert_eq!(sequential.sample_covariance(), parallel.sample_covariance());
    assert_eq!(sequential.mean(), parallel.mean());

    let moments = parallel.sample_covariance().unwrap();
    let one_variable: Vec<_> = samples
        .as_chunks::<3>()
        .0
        .iter()
        .map(|row| row[0])
        .collect();
    close(moments[0], descriptive::var(&one_variable, false).unwrap());
    for (&value, coefficient) in moments
        .iter()
        .zip([1., 2., -4., 2., 4., -8., -4., -8., 16.])
    {
        close(value, coefficient * moments[0]);
    }
    for (value, expected) in parallel
        .correlation()
        .unwrap()
        .into_iter()
        .zip([1., 1., -1., 1., 1., -1., -1., -1., 1.])
    {
        close(value, expected);
    }

    assert_eq!(
        parallel.push_batch(&[2.0, 3.0]),
        Err(OnlineStatsError::DimensionMismatch)
    );
    assert_eq!(
        parallel.push_batch(&[1.0, 2.0, 3.0, f64::NAN, 2.0, 3.0]),
        Err(OnlineStatsError::NonFiniteInput)
    );
    assert_eq!(
        parallel.push_batch(&[f64::MAX; 3]),
        Err(OnlineStatsError::Overflow)
    );
    parallel.push_batch(&[]).unwrap();
    assert_eq!(parallel.count(), 5_000);
    assert_eq!(parallel.sample_covariance().unwrap(), moments);
    let mut constant = OnlineCovariance::new(2).unwrap();
    constant.push_batch(&[2.0, 1.0, 2.0, 3.0]).unwrap();
    assert!(constant.correlation().is_none());
}

#[test]
fn parallel_gaussian_kde_matches_scalar_and_stays_finite_at_extreme_scales() {
    let data: Vec<_> = (0..128).map(|i| f64::from(i % 13) - 6.0).collect();
    let points: Vec<_> = (0..1024).map(|i| f64::from(i) / 100.0 - 5.0).collect();
    let bandwidth = robust::bandwidth_silverman(&data);
    let evaluate = || stats_ext::kde(&data, &points, None).unwrap();
    let (one, four) = (pool(1).install(evaluate), pool(4).install(evaluate));
    assert_eq!(one, four);
    for (&value, &point) in four.iter().zip(&points) {
        close(
            value,
            robust::kde(&data, point, bandwidth, robust::Kernel::Gaussian),
        );
    }
    let large_sample = [-f64::MAX, f64::MAX];
    assert!(robust::bandwidth_silverman(&large_sample).is_finite());
    assert!(stats_ext::kde(&large_sample, &[0.0], None).unwrap()[0] > 0.0);
}

#[test]
fn incremental_jackknife_buffer_preserves_callback_row_order() {
    let data = [10.0, 20.0, 30.0, 40.0];
    let positional_sum = |sample: &[f64]| {
        Ok(sample
            .iter()
            .zip([1.0, 2.0, 3.0, 4.0])
            .map(|(x, i)| x * i)
            .sum())
    };
    let result = stats_ext::jackknife(&data, &positional_sum).unwrap();
    assert_eq!(result.estimate, 300.0);
    assert_eq!(result.jackknife_estimates, vec![200.0, 190.0, 170.0, 140.0]);
}

#[test]
fn batch_normal_density_reuses_factor_without_changing_query_order() {
    let mean = [3.0, 4.0];
    let covariance = Mat {
        rows: 2,
        cols: 2,
        data: vec![2.0, 0.5, 0.5, 1.0],
    };
    let points = Mat {
        rows: 12_288,
        cols: 2,
        data: (0..12_288)
            .flat_map(|i| [f64::from(i % 20), f64::from(i % 30)])
            .collect(),
    };
    let evaluate = || multivariate_normal_log_pdfs(&mean, &covariance, &points).unwrap();
    let (one, four) = (pool(1).install(evaluate), pool(4).install(evaluate));
    assert_eq!(one, four);
    for (&value, point) in four.iter().zip(points.data.as_chunks::<2>().0) {
        assert_eq!(
            value,
            multivariate_normal_log_pdf(&mean, &covariance, point).unwrap()
        );
    }
    let empty = Mat {
        rows: 0,
        cols: 2,
        data: Vec::new(),
    };
    assert!(
        multivariate_normal_log_pdfs(&mean, &covariance, &empty)
            .unwrap()
            .is_empty()
    );
    let asymmetric = Mat {
        rows: 2,
        cols: 2,
        data: vec![1.0, 0.0, 0.3, 1.0],
    };
    assert_eq!(
        multivariate_normal_log_pdfs(&mean, &asymmetric, &points),
        Err(DistError::NotPositiveDefinite)
    );
    let malformed = Mat {
        rows: usize::MAX,
        cols: 2,
        data: vec![0.0; 2],
    };
    assert!(matches!(
        multivariate_normal_log_pdfs(&mean, &covariance, &malformed),
        Err(DistError::ShapeMismatch(_))
    ));
}

#[test]
fn distribution_evaluation_avoids_overflow_and_tail_cancellation() {
    let normal = Distribution::Normal {
        mu: -f64::MAX,
        sigma: f64::MAX,
    };
    let standard = Distribution::Normal { mu: 0.0, sigma: 1.0 };
    let density = normal.pdf(f64::MAX);
    assert!(density > 0.0);
    close(density * f64::MAX, standard.pdf(2.0));
    close(normal.cdf(f64::MAX), standard.cdf(2.0));
    assert!(
        Distribution::LogNormal {
            mu: 0.0,
            sigma: 1e300
        }
        .pdf(1e20)
            > 0.0
    );

    let geometric = Distribution::Geometric { p: 1e-20 };
    close(geometric.pdf(1e20) / 1e-20, (-1.0_f64).exp());
    let cauchy = Distribution::Cauchy { loc: 0.0, scale: 1.0 };
    close(cauchy.cdf(-1e50) * 1e50, 1.0 / std::f64::consts::PI);
    assert_eq!(Distribution::F { df1: 2.0, df2: 2.0 }.cdf(f64::MAX), 1.0);

    let covariance = Mat {
        rows: 1,
        cols: 1,
        data: vec![f64::MAX],
    };
    let result = multivariate_normal_log_pdf(&[-1e308], &covariance, &[1e308]).unwrap();
    assert!(result.is_finite() && result < -1e307);
}
