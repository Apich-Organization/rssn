//! Dense vector helpers (ported from `numerical_vector_test.rs`).

use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::kernels::vector::{
    angle, cosine_similarity, cross_product, distance, dot_product, is_orthogonal, is_parallel,
    l1_norm, lerp, linf_norm, lp_norm, norm, normalize, project, reflect, scalar_mul, vec_add,
    vec_sub,
};

fn cfg() -> ProptestConfig {
    ProptestConfig {
        rng_seed: RngSeed::Fixed(0x5EED),
        failure_persistence: None,
        ..ProptestConfig::default()
    }
}

#[test]
fn add_and_sub() {
    let (a, b) = ([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
    assert_eq!(vec_add(&a, &b), Ok(vec![5.0, 7.0, 9.0]));
    assert_eq!(vec_sub(&a, &b), Ok(vec![-3.0, -3.0, -3.0]));
    assert!(vec_add(&a, &[1.0]).is_err());
    assert!(vec_sub(&a, &[1.0]).is_err());
}

#[test]
fn scalar_multiplication() {
    assert_eq!(scalar_mul(&[1.0, 2.0, 3.0], 2.0), vec![2.0, 4.0, 6.0]);
}

#[test]
fn dot_product_value_and_mismatch() {
    assert_eq!(dot_product(&[1.0, 2.0, 3.0], &[4.0, 5.0, 6.0]), Ok(32.0));
    assert!(dot_product(&[1.0], &[1.0, 2.0]).is_err());
}

#[test]
fn norms() {
    assert_eq!(norm(&[3.0, 4.0]), 5.0);
    assert_eq!(l1_norm(&[3.0, -4.0]), 7.0);
    assert_eq!(linf_norm(&[3.0, -4.0]), 4.0);
    assert_approx_eq!(lp_norm(&[3.0, 4.0], 2.0), 5.0);
    assert_approx_eq!(lp_norm(&[1.0, 1.0, 1.0], 3.0), 3f64.powf(1.0 / 3.0));
    assert_eq!(lp_norm(&[3.0, -4.0], f64::INFINITY), 4.0);
    assert!(lp_norm(&[1.0], 0.0).is_nan());
    assert_eq!(norm(&[]), 0.0);
}

#[test]
fn cross_product_of_basis_vectors() {
    assert_eq!(cross_product(&[1.0, 0.0, 0.0], &[0.0, 1.0, 0.0]), Ok(vec![0.0, 0.0, 1.0]));
    assert_eq!(cross_product(&[0.0, 1.0, 0.0], &[1.0, 0.0, 0.0]), Ok(vec![0.0, 0.0, -1.0]));
    assert_eq!(cross_product(&[1.0, 2.0, 3.0], &[4.0, 5.0, 6.0]), Ok(vec![-3.0, 6.0, -3.0]));
    assert!(cross_product(&[1.0, 2.0], &[3.0, 4.0]).is_err());
}

#[test]
fn normalize_vector() {
    let r = normalize(&[3.0, 4.0]).unwrap_or_else(|e| panic!("{e}"));
    assert_approx_eq!(r[0], 0.6);
    assert_approx_eq!(r[1], 0.8);
    assert_approx_eq!(norm(&r), 1.0);
    assert!(normalize(&[0.0, 0.0]).is_err());
}

#[test]
fn projection() {
    assert_eq!(project(&[3.0, 4.0], &[1.0, 0.0]), Ok(vec![3.0, 0.0]));
    assert_eq!(project(&[3.0, 4.0], &[0.0, 0.0]), Ok(vec![0.0, 0.0]));
    assert!(project(&[3.0, 4.0], &[1.0]).is_err());
}

#[test]
fn reflection_about_a_unit_normal() {
    let r = reflect(&[1.0, -1.0], &[0.0, 1.0]).unwrap_or_else(|e| panic!("{e}"));
    assert_approx_eq!(r[0], 1.0);
    assert_approx_eq!(r[1], 1.0);
}

#[test]
fn angles() {
    assert_approx_eq!(angle(&[1.0, 0.0], &[0.0, 1.0]).unwrap_or(f64::NAN), std::f64::consts::FRAC_PI_2);
    assert_approx_eq!(angle(&[1.0, 0.0], &[-1.0, 0.0]).unwrap_or(f64::NAN), std::f64::consts::PI);
    assert_approx_eq!(angle(&[1.0, 0.0], &[1.0, 1.0]).unwrap_or(f64::NAN), std::f64::consts::FRAC_PI_4);
    assert_eq!(angle(&[0.0, 0.0], &[1.0, 1.0]), Ok(0.0));
}

#[test]
fn distance_lerp_and_predicates() {
    assert_approx_eq!(distance(&[0.0, 0.0], &[3.0, 4.0]).unwrap_or(f64::NAN), 5.0);
    assert_eq!(lerp(&[0.0, 10.0], &[10.0, 20.0], 0.25), Ok(vec![2.5, 12.5]));
    assert_eq!(is_orthogonal(&[1.0, 0.0], &[0.0, 3.0], 1e-9), Ok(true));
    assert_eq!(is_orthogonal(&[1.0, 1.0], &[0.0, 3.0], 1e-9), Ok(false));
    assert_eq!(is_parallel(&[1.0, 2.0], &[-2.0, -4.0], 1e-9), Ok(true));
    assert_eq!(is_parallel(&[1.0, 2.0], &[2.0, 1.0], 1e-9), Ok(false));
}

#[test]
fn cosine_similarity_values() {
    assert_approx_eq!(cosine_similarity(&[1.0, 0.0], &[1.0, 1.0]).unwrap_or(f64::NAN), std::f64::consts::FRAC_1_SQRT_2);
    assert!(cosine_similarity(&[0.0, 0.0], &[1.0, 1.0]).is_err());
}

fn same_len_pair() -> impl Strategy<Value = (Vec<f64>, Vec<f64>)> {
    (1..10usize).prop_flat_map(|n| {
        (
            prop::collection::vec(-100.0..100.0f64, n),
            prop::collection::vec(-100.0..100.0f64, n),
        )
    })
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn prop_add_commutes((a, b) in same_len_pair()) {
        prop_assert_eq!(vec_add(&a, &b), vec_add(&b, &a));
    }

    #[test]
    fn prop_add_then_sub_recovers((a, b) in same_len_pair()) {
        let sum = vec_add(&a, &b).map_err(TestCaseError::fail)?;
        let back = vec_sub(&sum, &b).map_err(TestCaseError::fail)?;
        for (x, y) in back.iter().zip(&a) {
            prop_assert!((x - y).abs() < 1e-9);
        }
    }

    #[test]
    fn prop_dot_commutes((a, b) in same_len_pair()) {
        prop_assert_eq!(dot_product(&a, &b), dot_product(&b, &a));
    }

    #[test]
    fn prop_cauchy_schwarz_and_triangle((a, b) in same_len_pair()) {
        let dot = dot_product(&a, &b).map_err(TestCaseError::fail)?;
        prop_assert!(dot.abs() <= norm(&a) * norm(&b) + 1e-9);
        let sum = vec_add(&a, &b).map_err(TestCaseError::fail)?;
        prop_assert!(norm(&sum) <= norm(&a) + norm(&b) + 1e-9);
    }

    #[test]
    fn prop_norm_scaling(v in prop::collection::vec(-100.0..100.0f64, 1..10), s in -10.0..10.0f64) {
        let scaled = norm(&scalar_mul(&v, s));
        prop_assert!((scaled - norm(&v) * s.abs()).abs() < 1e-9 * scaled.max(1.0));
    }

    #[test]
    fn prop_normalize_yields_unit_vector(v in prop::collection::vec(-100.0..100.0f64, 1..10)) {
        prop_assume!(norm(&v) > 1e-9);
        let unit = normalize(&v).map_err(TestCaseError::fail)?;
        prop_assert!((norm(&unit) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn prop_cross_product_is_orthogonal_to_operands(
        a in prop::collection::vec(-10.0..10.0f64, 3), b in prop::collection::vec(-10.0..10.0f64, 3),
    ) {
        let c = cross_product(&a, &b).map_err(TestCaseError::fail)?;
        prop_assert!(dot_product(&c, &a).map_err(TestCaseError::fail)?.abs() < 1e-8);
        prop_assert!(dot_product(&c, &b).map_err(TestCaseError::fail)?.abs() < 1e-8);
    }

    #[test]
    fn prop_reflection_preserves_norm(v in prop::collection::vec(-10.0..10.0f64, 3), n in prop::collection::vec(-10.0..10.0f64, 3)) {
        prop_assume!(norm(&n) > 1e-3);
        let n = normalize(&n).map_err(TestCaseError::fail)?;
        let r = reflect(&v, &n).map_err(TestCaseError::fail)?;
        prop_assert!((norm(&r) - norm(&v)).abs() < 1e-9 * norm(&v).max(1.0));
    }
}
