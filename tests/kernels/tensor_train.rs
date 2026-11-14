use rssn::kernels::tensor_train::TensorTrain;
use rssn::kernels::tensor_train::TrainCore;
use rssn::kernels::tensor_train::TrainError;
use rssn::kernels::tensor_train::tt_svd;

#[test]
fn compressed_core_import_round_trips_and_checks_untrusted_shapes() {
    let original = tt_svd(&[1.0, 2.0, 3.0, 4.0], &[2, 2], 0.0, None).unwrap();
    let cores = original
        .cores()
        .iter()
        .map(|c| TrainCore::new(c.shape(), c.data().to_vec()).unwrap())
        .collect();
    let imported = TensorTrain::from_cores(cores).unwrap();
    for (a, b) in original
        .to_dense()
        .unwrap()
        .iter()
        .zip(imported.to_dense().unwrap())
    {
        assert!((a - b).abs() < 1e-12);
    }
    assert!(imported.to_dense_with_limit(3).is_err());
    assert!(TrainCore::new([0, 2, 1], vec![]).is_err());
    assert!(TrainCore::new([1, 2, 1], vec![1.0]).is_err());
    assert!(TrainCore::new([1, 1, 1], vec![f64::NAN]).is_err());
    let a = TrainCore::new([1, 2, 2], vec![1.0; 4]).unwrap();
    let b = TrainCore::new([3, 2, 1], vec![1.0; 6]).unwrap();
    assert!(TensorTrain::from_cores(vec![a, b]).is_err());
    let high_order = TensorTrain::from_cores(
        (0..100)
            .map(|_| TrainCore::new([1, 2, 1], vec![1.0; 2]).unwrap())
            .collect(),
    )
    .unwrap();
    assert!((high_order.entry(&vec![1; 100]).unwrap() - 1.0).abs() < 1e-12);
    assert!(high_order.to_dense().is_err());
    assert!(high_order.storage_len() < 500);
}

#[test]
fn compressed_algebra_and_rounding_do_not_need_dense_intermediates() {
    let shape = [3, 4, 5];
    let a: Vec<_> = (0..60).map(|i| ((i * 17 % 23) as f64) - 11.0).collect();
    let b: Vec<_> = (0..60).map(|i| ((i * 7 % 19) as f64) - 9.0).collect();
    let ta = tt_svd(&a, &shape, 0.0, None).unwrap();
    let tb = tt_svd(&b, &shape, 0.0, None).unwrap();
    for (train, expected) in [
        (
            ta.add(&tb).unwrap(),
            a.iter().zip(&b).map(|(a, b)| a + b).collect::<Vec<_>>(),
        ),
        (
            ta.hadamard(&tb).unwrap(),
            a.iter().zip(&b).map(|(a, b)| a * b).collect::<Vec<_>>(),
        ),
    ] {
        let rounded = train.round(1e-12, None).unwrap();
        assert!(rounded.tolerance_met());
        let reconstructed = rounded.to_dense().unwrap();
        let error = expected
            .iter()
            .zip(&reconstructed)
            .fold(0.0_f64, |s, (a, b)| s.hypot(a - b));
        let norm = expected.iter().fold(0.0_f64, |s, &v| s.hypot(v));
        assert!(error / norm < 1e-11, "{error}/{norm}");
        assert!((rounded.frobenius_norm() / norm - 1.0).abs() < 1e-11);
        assert!(rounded.storage_len() <= train.storage_len());
    }
    let capped = ta.round(0.0, Some(1)).unwrap();
    assert!(!capped.tolerance_met());
    let error = a
        .iter()
        .zip(capped.to_dense().unwrap())
        .fold(0.0_f64, |s, (a, b)| s.hypot(a - b));
    assert!((error - capped.truncation_error()).abs() < 1e-10);
    let vector = tt_svd(&[1.0, 2.0], &[2], 0.0, None).unwrap();
    assert_eq!(
        vector.add(&vector).unwrap().to_dense().unwrap(),
        vec![2.0, 4.0]
    );
    assert_eq!(
        vector.hadamard(&vector).unwrap().to_dense().unwrap(),
        vec![1.0, 4.0]
    );
}

#[test]
fn tensor_train_exact_and_compressed_contractions() {
    let shape = [4, 5, 6, 3];
    let mut values = Vec::new();
    for a in 0..4 {
        for b in 0..5 {
            for c in 0..6 {
                for d in 0..3 {
                    values.push((a + 1) as f64 * (b + 2) as f64 * (c + 3) as f64 * (d + 1) as f64);
                }
            }
        }
    }
    let train = tt_svd(&values, &shape, 1e-12, None).unwrap();
    assert!(train.tolerance_met());
    assert!(train.storage_len() < values.len() / 3);
    let reconstructed = train.to_dense().unwrap();
    let error = values
        .iter()
        .zip(&reconstructed)
        .fold(0.0_f64, |s, (a, b)| s.hypot(a - b));
    let norm = values.iter().fold(0.0_f64, |s, &v| s.hypot(v));
    assert!((train.frobenius_norm() / norm - 1.0).abs() < 1e-12);
    assert!(error / norm < 1e-12, "{error} / {norm}");
    let dense_inner: f64 = values.iter().map(|x| x * x).sum();
    assert!((train.inner(&train).unwrap() / dense_inner - 1.0).abs() < 1e-12);
    assert!((train.entry(&[2, 3, 4, 1]).unwrap() - 3.0 * 5.0 * 7.0 * 2.0).abs() < 1e-10);
}

#[test]
fn tensor_train_rank_cap_reports_unmet_tolerance() {
    let values = [1.0, 0.0, 0.0, 2.0];
    let capped = tt_svd(&values, &[2, 2], 1e-12, Some(1)).unwrap();
    assert!(!capped.tolerance_met());
    assert!((capped.truncation_error() - 1.0).abs() < 1e-12);
    let dense = capped.to_dense().unwrap();
    let error = values
        .iter()
        .zip(dense)
        .fold(0.0_f64, |s, (a, b)| s.hypot(a - b));
    assert!((error - capped.truncation_error()).abs() < 1e-12);
}

#[test]
fn tensor_train_scale_zero_vector_and_invalid_inputs() {
    for scale in [1e-250, 1.0, 1e250] {
        let values = [scale, 2.0 * scale, 3.0 * scale, 6.0 * scale];
        let train = tt_svd(&values, &[2, 2], 1e-12, None).unwrap();
        for (a, b) in values.iter().zip(train.to_dense().unwrap()) {
            assert!((b / a - 1.0).abs() < 1e-12);
        }
    }
    let zero = tt_svd(&[0.0; 24], &[2, 3, 4], 0.0, None).unwrap();
    assert_eq!(zero.to_dense().unwrap(), vec![0.0; 24]);
    assert_eq!(zero.truncation_error(), 0.0);
    let vector = tt_svd(&[1.0, 2.0], &[2], 0.0, None).unwrap();
    assert_eq!(vector.to_dense().unwrap(), vec![1.0, 2.0]);
    assert!(matches!(
        tt_svd(&[f64::NAN], &[1], 0.0, None),
        Err(TrainError::NonFinite)
    ));
    assert!(tt_svd(&[], &[0], 0.0, None).is_err());
    assert!(tt_svd(&[1.0], &[1], 0.0, Some(0)).is_err());
    assert!(vector.entry(&[2]).is_err());
}

#[test]
fn dense_svd_rank_deficiency_orthogonality_and_scale() {
    use rssn::kernels::dense::Mat;
    use rssn::kernels::dense::svd;
    for scale in [1e-250, 1.0, 1e250] {
        let a = Mat {
            rows: 3,
            cols: 3,
            data: vec![
                scale,
                2.0 * scale,
                3.0 * scale,
                2.0 * scale,
                4.0 * scale,
                6.0 * scale,
                0.0,
                0.0,
                0.0,
            ],
        };
        let d = svd(&a).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                let back: f64 = (0..3).map(|k| d.u.at(i, k) * d.s[k] * d.v.at(j, k)).sum();
                assert!(((back - a.at(i, j)) / scale).abs() < 1e-12);
                let uu: f64 = (0..3).map(|k| d.u.at(k, i) * d.u.at(k, j)).sum();
                assert!((uu - if i == j { 1.0 } else { 0.0 }).abs() < 1e-12);
            }
        }
    }
}
