//! N-dimensional tensors (ported from `numerical_tensor_test.rs`).

use assert_approx_eq::assert_approx_eq;
use ndarray::{ArrayD, array};
use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::kernels::tensor::{
    TensorData, contract, inner_product, norm, outer_product, tensor_vec_mul, tensordot,
};

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
    assert_eq!(outer_product(&a, &c).map(|t| t.shape().to_vec()), Ok(vec![2, 3]));
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
    let bad = TensorData { shape: vec![3, 3], data: vec![1.0] };
    assert!(bad.to_arrayd().is_err());
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
