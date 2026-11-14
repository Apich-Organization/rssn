//! Numerical conditioning and scale invariance of the CP least-squares core.
use ndarray::ArrayD;
use ndarray::IxDyn;
use rssn::kernels::tensor::CpOptions;
use rssn::kernels::tensor::cp_als;
use rssn::kernels::tensor::norm;

fn tensor(scale: f64) -> ArrayD<f64> {
    ArrayD::from_shape_fn(IxDyn(&[3, 4, 5]), |i| {
        scale * ((i[0] as f64 - 0.4) * (0.3 + i[1] as f64) * (0.7 - i[2] as f64))
    })
}

#[test]
fn cp_rank_one_is_scale_invariant_and_reports_measured_error() {
    for scale in [1e-250, 1.0, 1e250] {
        let input = tensor(scale);
        let mut options = CpOptions::new(1);
        options.max_sweeps = 50;
        let result = cp_als(&input, &options).unwrap();
        assert!(result.converged, "scale {scale:e}: {result:?}");
        let error = input
            .iter()
            .zip(&result.reconstruction)
            .fold(0.0_f64, |n, (a, b)| n.hypot((a - b) / scale));
        let measured = error / (norm(&input) / scale);
        assert!(measured < 1e-12, "scale {scale:e}: {measured:e}");
        assert!((measured - result.relative_residual).abs() < 1e-14);
        assert!(result.weights.iter().all(|x| x.is_finite()));
        assert!(norm(&input).is_finite() && norm(&input) > 0.0);
    }
}

#[test]
fn cp_overcomplete_zero_tensor_reports_rank_deficiency() {
    let input = ArrayD::zeros(IxDyn(&[2, 2, 2]));
    let mut options = CpOptions::new(5);
    options.max_sweeps = 5;
    let result = cp_als(&input, &options).unwrap();
    assert!(result.converged && result.rank_deficient);
    assert_eq!(result.residual, 0.0);
    assert!(result.reconstruction.iter().all(|&x| x == 0.0));
}

#[test]
fn cp_correlated_components_descend_without_normal_equations() {
    let input = ArrayD::from_shape_fn(IxDyn(&[4, 5, 3]), |i| {
        let [a, b, c] = [i[0] as f64, i[1] as f64, i[2] as f64];
        (a + 0.3) * (b - 0.2) * (c + 0.4) + 0.7 * (a + 0.31) * (b - 0.18) * (c - 0.5)
    });
    let mut previous = f64::INFINITY;
    for sweeps in [1, 2, 5, 20, 100] {
        let mut options = CpOptions::new(2);
        options.max_sweeps = sweeps;
        options.tol = 0.0;
        let result = cp_als(&input, &options).unwrap();
        assert!(result.relative_residual <= previous + 1e-13);
        assert!(!result.diverged);
        previous = result.relative_residual;
    }
    assert!(previous < 1e-3, "{previous:e}");
}

#[test]
fn cp_rejects_invalid_numeric_options() {
    let input = tensor(1.0);
    for tol in [f64::NAN, f64::INFINITY, -1.0] {
        let mut options = CpOptions::new(1);
        options.tol = tol;
        assert!(cp_als(&input, &options).is_err());
    }
    for limit in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let mut options = CpOptions::new(1);
        options.divergence_limit = limit;
        assert!(cp_als(&input, &options).is_err());
    }
}

#[test]
fn cp_recovers_three_distinct_components_across_seeds() {
    let mut input = ArrayD::zeros(IxDyn(&[3, 4, 5]));
    input[[0, 0, 0]] = 1.0;
    input[[1, 1, 1]] = -0.7;
    input[[2, 2, 2]] = 0.4;
    for seed in [2, 7] {
        let mut options = CpOptions::new(3);
        options.seed = seed;
        options.max_sweeps = 3000;
        let result = cp_als(&input, &options).unwrap();
        assert!(
            result.converged && result.relative_residual < 1e-10,
            "{result:?}"
        );
        assert!(!result.diverged);
    }
}

#[test]
fn cp_border_rank_degeneracy_is_reported_honestly() {
    let mut input = ArrayD::zeros(IxDyn(&[2, 3, 4]));
    input[[1, 0, 0]] = 1.0;
    input[[0, 2, 0]] = 1.0;
    input[[0, 0, 3]] = 1.0;
    let mut options = CpOptions::new(2);
    options.seed = 3;
    options.max_sweeps = 400;
    options.divergence_limit = 1.0;
    let result = cp_als(&input, &options).unwrap();
    assert!(!result.converged);
    assert!(result.diverged || result.stagnated);
    let measured = input
        .iter()
        .zip(&result.reconstruction)
        .fold(0.0_f64, |s, (x, y)| s.hypot(x - y))
        / norm(&input);
    assert!((measured - result.relative_residual).abs() < 1e-12);
}
