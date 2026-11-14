use ndarray::ArrayD;
use ndarray::IxDyn;
use rssn::kernels::tensor::hooi;
use rssn::kernels::tensor::hosvd;

fn relative_error(
    a: &ArrayD<f64>,
    b: &ArrayD<f64>,
) -> f64 {
    a.iter().zip(b).fold(0.0_f64, |s, (x, y)| s.hypot(x - y))
        / a.iter().fold(0.0_f64, |s, &x| s.hypot(x))
}

#[test]
fn tucker_refinement_improves_a_zero_hosvd_projection() {
    let mut tensor = ArrayD::zeros(IxDyn(&[2, 2, 2]));
    tensor[IxDyn(&[0, 1, 1])] = 1.0;
    tensor[IxDyn(&[1, 0, 1])] = 1.0;
    tensor[IxDyn(&[1, 1, 0])] = 1.0;
    let initial = hosvd(&tensor, Some(&[1, 1, 1])).unwrap();
    let initial_error = relative_error(&tensor, &initial.reconstruct().unwrap());
    let fit = hooi(&tensor, &[1, 1, 1], 50, 1e-10).unwrap();
    assert!(initial_error > 0.99);
    assert!(fit.relative_error < 0.9, "{fit:?}");
    assert!(fit.converged);
    assert!(
        (relative_error(&tensor, &fit.reconstruct().unwrap()) - fit.relative_error).abs() < 1e-12
    );
    for pair in fit.error_history.windows(2) {
        assert!(pair[1] <= pair[0] + 1e-12);
    }
    for u in &fit.factors {
        for i in 0..u.cols {
            for j in 0..u.cols {
                let dot: f64 = (0..u.rows).map(|k| u.at(k, i) * u.at(k, j)).sum();
                assert!((dot - if i == j { 1.0 } else { 0.0 }).abs() < 1e-12);
            }
        }
    }
}

#[test]
fn tucker_refines_dense_data_and_reports_a_finite_budget() {
    let data = (0..120)
        .map(|i| (i as f64 * 0.37).sin() + ((i * i % 113) as f64 * 0.21).cos())
        .collect();
    let tensor = ArrayD::from_shape_vec(IxDyn(&[4, 5, 6]), data).unwrap();
    let initial = hosvd(&tensor, Some(&[2, 2, 2])).unwrap();
    let initial_error = relative_error(&tensor, &initial.reconstruct().unwrap());
    let fit = hooi(&tensor, &[2, 2, 2], 40, 1e-8).unwrap();
    assert!(fit.relative_error <= initial_error + 1e-12);
    assert!(
        (fit.relative_error - relative_error(&tensor, &fit.reconstruct().unwrap())).abs() < 1e-12
    );
    let budget = hooi(&tensor, &[2, 2, 2], 1, 0.0).unwrap();
    assert_eq!(budget.sweeps, 1);
    assert!(budget.relative_error <= initial_error + 1e-12);
}

#[test]
fn tucker_handles_scale_zero_scalar_and_rejects_invalid_inputs() {
    for scale in [1e-250, 1.0, 1e250] {
        let data = (0..24)
            .map(|i| scale * ((i / 12 + 1) * (i / 4 % 3 + 1) * (i % 4 + 1)) as f64)
            .collect();
        let tensor = ArrayD::from_shape_vec(IxDyn(&[2, 3, 4]), data).unwrap();
        let fit = hooi(&tensor, &[1, 1, 1], 10, 1e-10).unwrap();
        let normed = tensor.mapv(|x| x / scale);
        let back = fit.reconstruct().unwrap().mapv(|x| x / scale);
        assert!(relative_error(&normed, &back) < 1e-12);
        assert!(fit.converged);
    }
    let zeros = ArrayD::zeros(IxDyn(&[2, 3]));
    let fit = hooi(&zeros, &[2, 2], 5, 1e-8).unwrap();
    assert_eq!(fit.reconstruct().unwrap(), zeros);
    assert_eq!(fit.relative_error, 0.0);
    let scalar = ArrayD::from_elem(IxDyn(&[]), 3.0);
    assert_eq!(
        hooi(&scalar, &[], 1, 0.0).unwrap().reconstruct().unwrap(),
        scalar
    );
    assert!(hooi(&zeros, &[0, 1], 10, 1e-8).is_err());
    assert!(hooi(&zeros, &[1], 10, 1e-8).is_err());
    assert!(hooi(&zeros, &[1, 1], 0, 1e-8).is_err());
    assert!(hooi(&zeros, &[1, 1], 1, f64::NAN).is_err());
    assert!(hooi(&ArrayD::from_elem(IxDyn(&[1]), f64::NAN), &[1], 1, 0.0).is_err());
}
