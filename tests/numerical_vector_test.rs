use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;
use rssn::prelude::numerical::*;
use rssn::prelude::*;

#[test]

fn test_vec_add() {

    let v1 = vec![1.0, 2.0, 3.0];

    let v2 = vec![4.0, 5.0, 6.0];

    let res = numerical_vec_add(&v1, &v2).unwrap();

    assert_eq!(res, vec![5.0, 7.0, 9.0]);
}

#[test]

fn test_vec_sub() {

    let v1 = vec![1.0, 2.0, 3.0];

    let v2 = vec![4.0, 5.0, 6.0];

    let res = numerical_vec_sub(&v1, &v2).unwrap();

    assert_eq!(res, vec![-3.0, -3.0, -3.0]);
}

#[test]

fn test_scalar_mul() {

    let v = vec![1.0, 2.0, 3.0];

    let res = numerical_scalar_mul(&v, 2.0);

    assert_eq!(res, vec![2.0, 4.0, 6.0]);
}

#[test]

fn test_dot_product() {

    let v1 = vec![1.0, 2.0, 3.0];

    let v2 = vec![4.0, 5.0, 6.0];

    let res = numerical_dot_product(&v1, &v2).unwrap();

    assert_eq!(res, 32.0); // 1*4 + 2*5 + 3*6 = 4 + 10 + 18 = 32
}

#[test]

fn test_norm() {

    let v = vec![3.0, 4.0];

    let res = numerical_norm(&v);

    assert_eq!(res, 5.0);
}

#[test]

fn test_cross_product() {

    let v1 = vec![1.0, 0.0, 0.0];

    let v2 = vec![0.0, 1.0, 0.0];

    let res = numerical_cross_product(&v1, &v2).unwrap();

    assert_eq!(res, vec![0.0, 0.0, 1.0]);
}

#[test]

fn test_normalize() {

    let v = vec![3.0, 4.0];

    let res = numerical_normalize(&v).unwrap();

    assert_approx_eq!(res[0], 0.6);

    assert_approx_eq!(res[1], 0.8);

    assert_approx_eq!(numerical_norm(&res), 1.0);
}

#[test]

fn test_project() {

    let v1 = vec![3.0, 4.0];

    let v2 = vec![1.0, 0.0];

    let res = numerical_project(&v1, &v2).unwrap();

    assert_eq!(res, vec![3.0, 0.0]);
}

#[test]

fn test_reflect() {

    let v = vec![1.0, -1.0];

    let n = vec![0.0, 1.0]; // Normal pointing up
    let res = numerical_reflect(&v, &n).unwrap();

    assert_approx_eq!(res[0], 1.0);

    assert_approx_eq!(res[1], 1.0);
}

#[test]

fn test_angle() {

    let v1 = vec![1.0, 0.0];

    let v2 = vec![0.0, 1.0];

    let res = numerical_angle(&v1, &v2).unwrap();

    assert_approx_eq!(res, std::f64::consts::FRAC_PI_2);
}

// Property Tests
proptest! {
    #[test]
    fn prop_vec_add_commutes(v1 in prop::collection::vec(-1000.0..1000.0, 1..100),
                             v2 in prop::collection::vec(-1000.0..1000.0, 1..100)) {
        if v1.len() == v2.len() {
            let res1 = numerical_vec_add(&v1, &v2).unwrap();
            let res2 = numerical_vec_add(&v2, &v1).unwrap();
            for (a, b) in res1.iter().zip(res2.iter()) {
                  assert_approx_eq!(a, b);
            }
        }
    }

    #[test]
    fn prop_dot_product_commutes(v1 in prop::collection::vec(-100.0..100.0, 1..10),
                                 v2 in prop::collection::vec(-100.0..100.0, 1..10)) {
        if v1.len() == v2.len() {
            let res1 = numerical_dot_product(&v1, &v2).unwrap();
            let res2 = numerical_dot_product(&v2, &v1).unwrap();
            assert_approx_eq!(res1, res2);
        }
    }

    #[test]
    fn prop_triangle_inequality(v1 in prop::collection::vec(-100.0..100.0, 1..10),
                                v2 in prop::collection::vec(-100.0..100.0, 1..10)) {
        if v1.len() == v2.len() {
            let norm_v1 = numerical_norm(&v1);
            let norm_v2 = numerical_norm(&v2);
            let sum = numerical_vec_add(&v1, &v2).unwrap();
            let norm_sum = numerical_norm(&sum);
            assert!(norm_sum <= norm_v1 + norm_v2 + 1e-9);
        }
    }

    #[test]
    fn prop_norm_scaling(v in prop::collection::vec(-100.0..100.0, 1..10),
                         s in -10.0..10.0) {
        let n_v = numerical_norm(&v);
        let scaled_v = numerical_scalar_mul(&v, s);
        let n_scaled = numerical_norm(&scaled_v);
        assert_approx_eq!(n_scaled, n_v * s.abs());
    }

    #[test]
    fn prop_normalize_yields_unit_vector(v in prop::collection::vec(-100.0..100.0, 1..10)) {
        let n = numerical_norm(&v);
        if n > 1e-9 {
            let unit = numerical_normalize(&v).unwrap();
            assert_approx_eq!(numerical_norm(&unit), 1.0);
        }
    }
}
