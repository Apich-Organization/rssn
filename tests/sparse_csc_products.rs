//! Reference cases for direct forward and transposed CSC products.

use faer::sparse::SymbolicSparseColMat;
use rssn::kernels::sparse::SparseMat;
use rssn::kernels::sparse::csr_from_triplets;
use rssn::kernels::sparse::sp_mat_transpose_vec_mul;
use rssn::kernels::sparse::sp_mat_transpose_vec_mul_into;
use rssn::kernels::sparse::sp_mat_vec_mul;
use rssn::kernels::sparse::sp_mat_vec_mul_into;

#[test]
fn rectangular_products_match_an_independently_assembled_dense_matrix() {
    let triplets = [
        (2, 1, -3.0),
        (0, 4, 5.25),
        (1, 0, 2.0),
        (0, 0, -0.25),
        (2, 1, 1.5),
        (1, 2, 0.5),
        (0, 1, -3.0),
    ];
    let mut dense = [[0.0; 5]; 3];
    for (row, column, value) in triplets {
        dense[row][column] += value;
    }
    let matrix = csr_from_triplets(dense.len(), dense[0].len(), &triplets);
    let x = [1.5, -2.0, 7.0, 3.0, -5.0];
    let y = [-1.25, 2.75, 0.5];
    let expected_forward = dense.map(|row| row.into_iter().zip(x).map(|(a, x)| a * x).sum::<f64>());
    let expected_transposed: Vec<_> = (0..x.len())
        .map(|column| {
            dense
                .iter()
                .zip(y)
                .map(|(row, y)| row[column] * y)
                .sum::<f64>()
        })
        .collect();
    assert_eq!(sp_mat_vec_mul(&matrix, &x).unwrap(), expected_forward);
    assert_eq!(
        sp_mat_transpose_vec_mul(&matrix, &y).unwrap(),
        expected_transposed
    );
    let mut forward = [f64::NAN; 3];
    let mut reverse = [f64::NAN; 5];
    sp_mat_vec_mul_into(&matrix, &x, &mut forward).unwrap();
    sp_mat_transpose_vec_mul_into(&matrix, &y, &mut reverse).unwrap();
    assert_eq!(forward, expected_forward);
    assert_eq!(reverse.as_slice(), expected_transposed);
    assert_eq!(
        forward.into_iter().zip(y).map(|(a, b)| a * b).sum::<f64>(),
        reverse.into_iter().zip(x).map(|(a, b)| a * b).sum::<f64>()
    );

    sp_mat_vec_mul_into(&matrix, &[0.0; 5], &mut forward).unwrap();
    sp_mat_transpose_vec_mul_into(&matrix, &[0.0; 3], &mut reverse).unwrap();
    assert_eq!(forward, [0.0; 3]);
    assert_eq!(reverse, [0.0; 5]);
}

#[test]
fn both_products_skip_unused_column_capacity_and_outside_slots() {
    let padding = usize::MAX;
    let symbolic = SymbolicSparseColMat::new_checked(
        3,
        4,
        vec![1, 4, 6, 9, 11],
        Some(vec![2, 0, 1, 2]),
        vec![
            padding, 0, 2, padding, padding, padding, 1, padding, 0, 0, 2, padding,
        ],
    );
    let matrix = SparseMat::new(
        symbolic,
        vec![
            f64::NAN,
            1.5,
            -2.0,
            f64::NAN,
            f64::NAN,
            f64::NAN,
            4.0,
            f64::NAN,
            -777.0,
            5.0,
            -1.0,
            f64::NAN,
        ],
    );
    assert_eq!(matrix.compute_nnz(), 5);
    let x = [2.0, f64::NAN, -3.0, 0.75];
    let y = [0.5, -1.25, 3.0];
    assert_eq!(sp_mat_vec_mul(&matrix, &x).unwrap(), [6.75, -12.0, -4.75]);
    assert_eq!(
        sp_mat_transpose_vec_mul(&matrix, &y).unwrap(),
        [-5.25, 0.0, -5.0, -0.5]
    );
    let mut forward = [5.0; 3];
    let mut reverse = [5.0; 4];
    sp_mat_vec_mul_into(&matrix, &x, &mut forward).unwrap();
    sp_mat_transpose_vec_mul_into(&matrix, &y, &mut reverse).unwrap();
    assert_eq!(forward, [6.75, -12.0, -4.75]);
    assert_eq!(reverse, [-5.25, 0.0, -5.0, -0.5]);
}

#[test]
fn empty_rows_columns_and_dimensions_are_zeroed() {
    for (rows, cols) in [(0, 0), (3, 0), (0, 4), (3, 4)] {
        let matrix = csr_from_triplets(rows, cols, &[]);
        let x = vec![f64::NAN; cols];
        let y = vec![f64::NAN; rows];
        let mut forward = vec![f64::INFINITY; rows];
        let mut reverse = vec![f64::INFINITY; cols];
        sp_mat_vec_mul_into(&matrix, &x, &mut forward).unwrap();
        sp_mat_transpose_vec_mul_into(&matrix, &y, &mut reverse).unwrap();
        assert_eq!(forward, vec![0.0; rows]);
        assert_eq!(reverse, vec![0.0; cols]);
        assert_eq!(sp_mat_vec_mul(&matrix, &x).unwrap(), forward);
        assert_eq!(sp_mat_transpose_vec_mul(&matrix, &y).unwrap(), reverse);
    }
}

#[test]
fn shape_errors_leave_caller_owned_slices_untouched() {
    let matrix = csr_from_triplets(3, 4, &[(2, 1, 2.0)]);
    let mut tmp = [42.0; 5];
    assert!(sp_mat_vec_mul(&matrix, &[0.0; 3]).is_err());
    assert!(sp_mat_transpose_vec_mul(&matrix, &[0.0; 4]).is_err());
    assert!(sp_mat_vec_mul_into(&matrix, &[1.0; 4], &mut tmp[..2]).is_err());
    assert!(sp_mat_vec_mul_into(&matrix, &[1.0; 3], &mut tmp[..3]).is_err());
    assert!(sp_mat_transpose_vec_mul_into(&matrix, &[1.0; 3], &mut tmp[..3]).is_err());
    assert!(sp_mat_transpose_vec_mul_into(&matrix, &[1.0; 4], &mut tmp[..4]).is_err());
    assert_eq!(tmp, [42.0; 5]);
}

#[test]
fn products_keep_ieee_behavior_for_explicit_zeros_and_stored_order() {
    let symbolic = SymbolicSparseColMat::new_checked(2, 2, vec![0, 1, 1], None, vec![0]);
    let matrix = SparseMat::new(symbolic, vec![0.0]);
    let forward = sp_mat_vec_mul(&matrix, &[f64::NAN, 1.0]).unwrap();
    let reverse = sp_mat_transpose_vec_mul(&matrix, &[f64::NAN, 1.0]).unwrap();
    for product in [forward, reverse] {
        assert!(product[0].is_nan());
        assert_eq!(product[1].to_bits(), 0.0_f64.to_bits());
    }

    let symbolic =
        SymbolicSparseColMat::new_unsorted_checked(3, 1, vec![0, 3], None, vec![0, 2, 1]);
    let matrix = SparseMat::new(symbolic, vec![1.0e16, -1.0e16, 1.0]);
    assert_eq!(sp_mat_transpose_vec_mul(&matrix, &[1.0; 3]).unwrap(), [1.0]);
}
