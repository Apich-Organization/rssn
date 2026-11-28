//! Tensor contraction against a scalar multi-index definition.

use ndarray::Array2;
use ndarray::Array3;
use ndarray::Array4;
use ndarray::ArrayD;
use ndarray::ArrayViewD;
use ndarray::IxDyn;
use ndarray::ShapeBuilder;
use ndarray::s;
use rssn::kernels::tensor::tensor_vec_mul;
use rssn::kernels::tensor::tensordot;
use rssn::kernels::tensor::tensordot_views;

fn scalar_reference(
    a: ArrayViewD<'_, f64>,
    b: ArrayViewD<'_, f64>,
    axes_a: &[usize],
    axes_b: &[usize],
) -> ArrayD<f64> {
    let free_a: Vec<_> = (0..a.ndim()).filter(|i| !axes_a.contains(i)).collect();
    let free_b: Vec<_> = (0..b.ndim()).filter(|i| !axes_b.contains(i)).collect();
    let output: Vec<_> = free_a
        .iter()
        .map(|&i| a.shape()[i])
        .chain(free_b.iter().map(|&j| b.shape()[j]))
        .collect();
    let contracted: Vec<_> = axes_a.iter().map(|&i| a.shape()[i]).collect();
    ArrayD::from_shape_fn(IxDyn(&output), |index| {
        let mut ai = vec![0; a.ndim()];
        let mut bi = vec![0; b.ndim()];
        for (j, &axis) in free_a.iter().enumerate() {
            ai[axis] = index[j];
        }
        for (j, &axis) in free_b.iter().enumerate() {
            bi[axis] = index[free_a.len() + j];
        }
        ndarray::indices(IxDyn(&contracted))
            .into_iter()
            .map(|index| {
                for (j, (&axis_a, &axis_b)) in axes_a.iter().zip(axes_b).enumerate() {
                    ai[axis_a] = index[j];
                    bi[axis_b] = index[j];
                }
                a[IxDyn(&ai)] * b[IxDyn(&bi)]
            })
            .sum()
    })
}

#[test]
fn reversed_strided_unequal_ranks_keep_free_and_paired_axis_orders() {
    let left = Array3::from_shape_fn((6, 5, 4), |(i, j, k)| (i * 11 + j * 3 + k) as f64 - 35.0);
    let right = Array4::from_shape_fn((4, 2, 3, 6), |(i, j, k, l)| {
        (i * 9 + j * 3 + k * 2 + l) as f64 - 17.0
    });
    let a = left.slice(s![..;-2, .., ..;-1]).into_dyn();
    let b = right.slice(s![..;-1, .., .., ..;2]).into_dyn();
    assert_eq!(a.shape(), &[3, 5, 4]);
    assert_eq!(b.shape(), &[4, 2, 3, 3]);
    for (aa, bb) in [([2, 0], [0, 2]), ([0, 2], [2, 0])] {
        let expected = scalar_reference(a.view(), b.view(), &aa, &bb);
        let got = tensordot_views(a.view(), b.view(), &aa, &bb).unwrap();
        assert_eq!(got.shape(), &[5, 2, 3]);
        assert_eq!(got, expected);
        assert_eq!(
            tensordot(&a.to_owned(), &b.to_owned(), &aa, &bb).unwrap(),
            expected
        );
    }
    assert!(tensordot_views(a.view(), b.view(), &[2, 2], &[0, 2]).is_err());
    assert!(tensordot_views(a.view(), b.view(), &[2, 1], &[0, 2]).is_err());
    assert!(tensordot_views(a.view(), b.view(), &[2], &[4]).is_err());
}

#[test]
fn column_major_and_broadcast_views_match_logical_indexing() {
    let a = Array2::from_shape_fn((4, 3).f(), |(i, j)| (i * 4 + j) as f64 - 5.0).into_dyn();
    let b = Array2::from_shape_fn((3, 2).f(), |(i, j)| (i * 3 + j) as f64 + 1.0).into_dyn();
    let got = tensordot_views(a.view(), b.view(), &[1], &[0]).unwrap();
    assert_eq!(got, scalar_reference(a.view(), b.view(), &[1], &[0]));

    let singleton = Array3::from_shape_fn((1, 3, 1), |(_, i, _)| i as f64 - 1.0).into_dyn();
    let broadcast = singleton.broadcast(IxDyn(&[2, 3, 4])).unwrap();
    let got = tensordot_views(broadcast.view(), b.view(), &[1], &[0]).unwrap();
    assert_eq!(got.shape(), &[2, 4, 2]);
    assert_eq!(
        got,
        scalar_reference(broadcast.view(), b.view(), &[1], &[0])
    );

    let source = a.clone().slice_move(s![..;-1, ..;-1]).into_dyn();
    let v = [1.0, -2.0, 3.0];
    let expected: Vec<_> = source
        .outer_iter()
        .map(|row| row.iter().zip(v).map(|(x, y)| x * y).sum::<f64>())
        .collect();
    assert_eq!(
        tensor_vec_mul(&source, &v).unwrap().as_slice().unwrap(),
        expected
    );
}

#[test]
fn distinguishes_empty_free_axes_empty_sums_and_rank_zero() {
    let scalar = ArrayD::from_elem(IxDyn(&[]), -2.0);
    let vector = ArrayD::from_shape_vec(IxDyn(&[3]), vec![1.0, 2.0, -3.0]).unwrap();
    let outer = tensordot(&scalar, &vector, &[], &[]).unwrap();
    assert_eq!(outer.shape(), &[3]);
    assert_eq!(outer.as_slice().unwrap(), &[-2.0, -4.0, 6.0]);
    let scalar_result = tensordot(&vector, &vector, &[0], &[0]).unwrap();
    assert_eq!(scalar_result.ndim(), 0);
    assert_eq!(scalar_result[IxDyn(&[])], 14.0);

    let left = ArrayD::zeros(IxDyn(&[2, 0, 4]));
    let right = ArrayD::zeros(IxDyn(&[3, 0]));
    let result = tensordot(&left, &right, &[1], &[1]).unwrap();
    assert_eq!(result.shape(), &[2, 4, 3]);
    assert_eq!(
        result,
        scalar_reference(left.view(), right.view(), &[1], &[1])
    );
    assert!(result.iter().all(|v| v.to_bits() == 0));

    let no_values = ArrayD::zeros(IxDyn(&[0]));
    let scalar_zero = tensordot(&no_values, &no_values, &[0], &[0]).unwrap();
    assert_eq!(scalar_zero.ndim(), 0);
    assert_eq!(scalar_zero[IxDyn(&[])].to_bits(), 0);

    let left = ArrayD::zeros(IxDyn(&[0, 3]));
    let result = tensordot(&left, &vector, &[1], &[0]).unwrap();
    assert_eq!(result.shape(), &[0]);
    assert!(result.is_empty());
    let result = tensor_vec_mul(&ArrayD::zeros(IxDyn(&[2, 0])), &[]).unwrap();
    assert_eq!(result.as_slice().unwrap(), &[0.0, 0.0]);
}

#[test]
fn output_shape_overflow_returns_an_error_even_without_input_storage() {
    let wide = (isize::MAX as usize) / 16 + 1;
    let a = ArrayD::<f64>::zeros(IxDyn(&[wide, 0]));
    let b = ArrayD::<f64>::zeros(IxDyn(&[0, 3]));
    assert!(tensordot(&a, &b, &[1], &[0]).is_err());

    let a = ArrayD::<f64>::zeros(IxDyn(&[0, wide, 0]));
    let b = ArrayD::<f64>::zeros(IxDyn(&[0, wide]));
    assert!(tensordot(&a, &b, &[2], &[0]).is_err());
}
