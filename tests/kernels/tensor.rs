//! N-dimensional tensors (ported from `numerical_tensor_test.rs`).

use assert_approx_eq::assert_approx_eq;
use ndarray::ArrayD;
use ndarray::array;
use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::kernels::tensor::TensorData;
use rssn::kernels::tensor::contract;
use rssn::kernels::tensor::hosvd;
use rssn::kernels::tensor::inner_product;
use rssn::kernels::tensor::norm;
use rssn::kernels::tensor::outer_product;
use rssn::kernels::tensor::tensor_vec_mul;
use rssn::kernels::tensor::tensordot;

fn cfg() -> ProptestConfig {
    ProptestConfig {
        rng_seed: RngSeed::Fixed(0x5EED),
        failure_persistence: None,
        ..ProptestConfig::default()
    }
}

#[test]
fn frobenius_norm() {
    let a = array![[1.0, 2.0], [3.0, 4.0]].into_dyn();
    assert_approx_eq!(norm(&a), 30f64.sqrt());
}

#[test]
fn inner_and_outer_products() {
    let a = array![1.0, 2.0].into_dyn();
    let b = array![3.0, 4.0].into_dyn();
    assert_eq!(inner_product(&a, &b), Ok(11.0));
    let outer = outer_product(&a, &b).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(outer.shape(), &[2, 2]);
    assert_eq!(outer, array![[3.0, 4.0], [6.0, 8.0]].into_dyn());
    let c = array![1.0, 2.0, 3.0].into_dyn();
    assert!(inner_product(&a, &c).is_err());
    assert_eq!(
        outer_product(&a, &c).map(|t| t.shape().to_vec()),
        Ok(vec![2, 3])
    );
}

#[test]
fn tensor_vector_product() {
    let a = array![[1.0, 2.0], [3.0, 4.0]].into_dyn();
    let res = tensor_vec_mul(&a, &[1.0, 2.0]).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(res.shape(), &[2]);
    assert_eq!(res[[0]], 5.0);
    assert_eq!(res[[1]], 11.0);
    assert!(tensor_vec_mul(&a, &[1.0, 2.0, 3.0]).is_err());
}

#[test]
fn contraction_of_a_matrix_is_its_trace() {
    let a = array![[1.0, 2.0], [3.0, 4.0]].into_dyn();
    let t = contract(&a, 0, 1).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(t.iter().copied().sum::<f64>(), 5.0);
    assert!(contract(&a, 0, 0).is_err());
    let rect = array![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]].into_dyn();
    assert!(contract(&rect, 0, 1).is_err());
}

#[test]
fn contraction_of_higher_rank_tensors() {
    // a[i][j][k] = 100 i + 10 j + k; contracting axes 0 and 2 leaves
    // r[j] = sum_t a[t][j][t].
    let a = ndarray::Array3::from_shape_fn((3, 2, 3), |(i, j, k)| (100 * i + 10 * j + k) as f64)
        .into_dyn();
    let r = contract(&a, 2, 0).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(r.shape(), &[2]);
    for j in 0..2 {
        let want: f64 = (0..3).map(|t| (100 * t + 10 * j + t) as f64).sum();
        assert_eq!(r[[j]], want);
    }
    assert!(contract(&a, 0, 1).is_err());
    assert!(contract(&a, 0, 5).is_err());
}

#[test]
fn tensordot_of_matrices_is_matrix_product() {
    let a = array![[1.0, 2.0], [3.0, 4.0]].into_dyn();
    let b = array![[5.0, 6.0], [7.0, 8.0]].into_dyn();
    let p = tensordot(&a, &b, &[1], &[0]).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(p, array![[19.0, 22.0], [43.0, 50.0]].into_dyn());
    assert!(tensordot(&a, &b, &[1], &[0, 1]).is_err());
}

#[test]
fn serde_round_trip() {
    let a = array![[1.0, 2.0], [3.0, 4.0]].into_dyn();
    let json = serde_json::to_string(&TensorData::from(&a)).unwrap_or_else(|e| panic!("{e}"));
    let decoded: TensorData = serde_json::from_str(&json).unwrap_or_else(|e| panic!("{e}"));
    let back = decoded.to_arrayd().unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(back.shape(), &[2, 2]);
    assert_eq!(back, a);
    let bad = TensorData {
        shape: vec![3, 3],
        data: vec![1.0],
    };
    assert!(bad.to_arrayd().is_err());
}

/// The inner product of two slices of `core` along `mode`.
fn slice_inner(
    core: &ArrayD<f64>,
    mode: usize,
    first: usize,
    second: usize,
) -> f64 {
    let ndim = core.ndim();
    let mut permutation = vec![mode];
    permutation.extend((0..ndim).filter(|&axis| axis != mode));
    let standard = core
        .clone()
        .permuted_axes(permutation)
        .as_standard_layout()
        .into_owned();
    let data = standard.into_raw_vec_and_offset().0;
    let slices = core.shape()[mode];
    let block = data.len() / slices;
    (0..block)
        .map(|offset| data[first * block + offset] * data[second * block + offset])
        .sum()
}

#[test]
fn hosvd_is_exact_full_rank_with_orthonormal_factors() {
    let tensor = ndarray::Array3::from_shape_fn((3, 4, 5), |(i, j, k)| {
        0.5 + ((i * 7 + j * 3 + k) % 5) as f64 + 0.25 * ((i + 2 * j + 3 * k) as f64).sin()
    })
    .into_dyn();
    let decomposition = hosvd(&tensor, None).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(decomposition.ranks, vec![3, 4, 5]);
    assert_eq!(decomposition.core.shape(), &[3, 4, 5]);
    assert_eq!(decomposition.factors.len(), 3);
    for (mode, factor) in decomposition.factors.iter().enumerate() {
        assert_eq!(factor.cols, decomposition.ranks[mode]);
        assert_eq!(factor.rows, tensor.shape()[mode]);
        // Orthonormal columns: U^T U = I.
        let gram = factor
            .transpose()
            .matmul(factor)
            .unwrap_or_else(|e| panic!("{e}"));
        for i in 0..gram.rows {
            for j in 0..gram.cols {
                let want = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (gram.at(i, j) - want).abs() < 1e-10,
                    "mode {mode}: U^T U[{i},{j}] = {}",
                    gram.at(i, j)
                );
            }
        }
        // All-orthogonality of the core along this mode.
        let slices = decomposition.core.shape()[mode];
        for a in 0..slices {
            for b in (a + 1)..slices {
                let inner = slice_inner(&decomposition.core, mode, a, b);
                assert!(inner.abs() < 1e-12, "mode {mode}: <C_{a}, C_{b}> = {inner}");
            }
        }
    }
    let back = decomposition
        .reconstruct()
        .unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(back.shape(), tensor.shape());
    let worst = tensor
        .iter()
        .zip(back.iter())
        .map(|(left, right)| (left - right).abs())
        .fold(0.0, f64::max);
    assert!(worst < 1e-10, "reconstruction error {worst}");
}

#[test]
fn hosvd_truncation_is_exact_for_a_rank_one_tensor() {
    let u = [1.0, -2.0, 0.5];
    let v = [0.3, 1.5, -2.5, 0.7];
    let w = [2.0, -1.0, 0.5, 3.0, -0.5];
    let tensor =
        ndarray::Array3::from_shape_fn((3, 4, 5), |(i, j, k)| u[i] * v[j] * w[k]).into_dyn();
    let decomposition = hosvd(&tensor, Some(&[1, 1, 1])).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(decomposition.ranks, vec![1, 1, 1]);
    assert_eq!(decomposition.core.shape(), &[1, 1, 1]);
    assert!(decomposition.singular_values[0][0] > 0.0);
    let back = decomposition
        .reconstruct()
        .unwrap_or_else(|e| panic!("{e}"));
    let worst = tensor
        .iter()
        .zip(back.iter())
        .map(|(left, right)| (left - right).abs())
        .fold(0.0, f64::max);
    assert!(worst < 1e-10, "rank-one truncation error {worst}");
}

#[test]
fn hosvd_checks_ranks_and_clamps_them() {
    let tensor = ndarray::Array2::from_shape_fn((3, 4), |(i, j)| (i + j) as f64).into_dyn();
    assert!(hosvd(&tensor, Some(&[2])).is_err(), "one rank per mode");
    assert!(hosvd(&tensor, Some(&[0, 2])).is_err(), "a zero rank");
    // `i + j` has numerical rank two, so a larger request is clamped to it.
    let clamped = hosvd(&tensor, Some(&[9, 2])).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(
        clamped.ranks,
        vec![2, 2],
        "ranks are clamped to the numerical rank of the unfolding"
    );
    assert_eq!(clamped.core.shape(), &[2, 2]);
    let back = clamped.reconstruct().unwrap_or_else(|e| panic!("{e}"));
    let error: f64 = tensor
        .iter()
        .zip(back.iter())
        .map(|(left, right)| (left - right).powi(2))
        .sum::<f64>()
        .sqrt();
    assert!(
        error < 1e-12,
        "rank-two reconstruction of a rank-two tensor: {error}"
    );
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn prop_norm_scaling(s in -100.0..100.0f64) {
        let a: ArrayD<f64> = array![1.0, 2.0, 3.0].into_dyn();
        let sa = (&a * s).into_dyn();
        prop_assert!((norm(&sa) - s.abs() * norm(&a)).abs() < 1e-9 * (1.0 + s.abs()));
    }

    #[test]
    fn prop_outer_product_norm_factorises(
        u in prop::collection::vec(-10.0..10.0f64, 1..5), v in prop::collection::vec(-10.0..10.0f64, 1..5),
    ) {
        let a = ArrayD::from_shape_vec(vec![u.len()], u).map_err(|e| TestCaseError::fail(e.to_string()))?;
        let b = ArrayD::from_shape_vec(vec![v.len()], v).map_err(|e| TestCaseError::fail(e.to_string()))?;
        let o = outer_product(&a, &b).map_err(TestCaseError::fail)?;
        prop_assert!((norm(&o) - norm(&a) * norm(&b)).abs() < 1e-9 * (1.0 + norm(&o)));
    }

    #[test]
    fn prop_inner_product_with_self_is_squared_norm(u in prop::collection::vec(-10.0..10.0f64, 1..8)) {
        let a = ArrayD::from_shape_vec(vec![u.len()], u).map_err(|e| TestCaseError::fail(e.to_string()))?;
        let ip = inner_product(&a, &a).map_err(TestCaseError::fail)?;
        prop_assert!((ip - norm(&a).powi(2)).abs() < 1e-9 * (1.0 + ip));
    }
}
