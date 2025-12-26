use ::ndarray::array;
use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;
use rssn::prelude::numerical::*;
use rssn::prelude::*;

#[test]

fn test_tensor_basic() {

    let a = array![
        [1.0, 2.0],
        [3.0, 4.0]
    ]
    .into_dyn();

    assert_eq!(
        numerical_tensor_norm(&a),
        (1.0 + 4.0 + 9.0 + 16.0f64)
            .sqrt()
    );
}

#[test]

fn test_tensor_inner_outer() {

    let a = array![1.0, 2.0].into_dyn();

    let b = array![3.0, 4.0].into_dyn();

    assert_eq!(
        numerical_tensor_inner_product(
            &a, &b
        )
        .unwrap(),
        11.0
    );

    let outer =
        numerical_outer_product(&a, &b)
            .unwrap();

    assert_eq!(
        outer.shape(),
        &[2, 2]
    );

    assert_eq!(outer[[0, 0]], 3.0);

    assert_eq!(outer[[1, 1]], 8.0);
}

#[test]

fn test_tensor_vec_mul() {

    let a = array![
        [1.0, 2.0],
        [3.0, 4.0]
    ]
    .into_dyn();

    let v = vec![1.0, 2.0];

    let res = numerical_tensor_vec_mul(
        &a, &v,
    )
    .unwrap();

    assert_eq!(res.shape(), &[2]);

    assert_eq!(res[[0]], 5.0);

    assert_eq!(res[[1]], 11.0);
}

#[test]

fn test_tensor_serde() {

    let a = array![
        [1.0, 2.0],
        [3.0, 4.0]
    ]
    .into_dyn();

    let data =
        numerical_TensorData::from(&a);

    let json =
        serde_json::to_string(&data)
            .unwrap();

    let decoded : numerical_TensorData =
        serde_json::from_str(&json)
            .unwrap();

    let arr_back = decoded
        .to_arrayd()
        .unwrap();

    assert_eq!(
        arr_back.shape(),
        &[2, 2]
    );

    assert_eq!(
        arr_back[[1, 1]],
        4.0
    );
}

proptest! {
    #[test]
    fn prop_tensor_norm_scaling(s in -100.0..100.0f64) {
        let a = array![1.0, 2.0, 3.0].into_dyn();
        let sa = &a * s;
        assert_approx_eq!(numerical_tensor_norm(&sa.into_dyn()), s.abs() * numerical_tensor_norm(&a), 1e-9);
    }
}
