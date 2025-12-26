use ::ndarray::Array1;
use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;
use rssn::prelude::numerical::*;
use rssn::prelude::*;

proptest! {
    #[test]
    fn prop_sparse_transpose_transpose(
        rows in 1..20usize,
        cols in 1..20usize,
        elements in proptest::collection::vec((0..20usize, 0..20usize, -100.0..100.0f64), 0..20)
    ) {
        let triplets: Vec<(usize, usize, f64)> = elements.into_iter()
            .map(|(r, c, v)| (r % rows, c % cols, v))
            .collect();
        let mat = numerical_csr_from_triplets(rows, cols, &triplets);
        let tt = numerical_transpose(&numerical_transpose(&mat));
        assert_eq!(mat.rows(), tt.rows());
        assert_eq!(mat.cols(), tt.cols());
        assert_eq!(mat.nnz(), tt.nnz());
        // Check data
        for (val, (r, c)) in mat.iter() {
            assert_approx_eq!(val, tt.get(r, c).cloned().unwrap_or(0.0), 1e-9);
        }
    }

    #[test]
    fn prop_sparse_trace_transpose(
        size in 1..20usize,
        elements in proptest::collection::vec((0..20usize, 0..20usize, -100.0..100.0f64), 0..20)
    ) {
        let triplets: Vec<(usize, usize, f64)> = elements.into_iter()
            .map(|(r, c, v)| (r % size, c % size, v))
            .collect();
        let mat = numerical_csr_from_triplets(size, size, &triplets);
        let t1 = numerical_trace(&mat).unwrap();
        let t2 = numerical_trace(&numerical_transpose(&mat)).unwrap();
        assert_approx_eq!(t1, t2, 1e-9);
    }
}

#[test]

fn test_sparse_basic() {

    let triplets = vec![
        (0, 0, 1.0),
        (1, 1, 2.0),
        (2, 2, 3.0),
    ];

    let mat = numerical_csr_from_triplets(3, 3, &triplets);

    assert_eq!(mat.rows(), 3);

    assert_eq!(mat.cols(), 3);

    assert_eq!(mat.nnz(), 3);
}

#[test]

fn test_sparse_spmv() {

    let triplets = vec![
        (0, 0, 1.0),
        (0, 2, 2.0),
        (1, 1, 3.0),
    ];

    let mat = numerical_csr_from_triplets(3, 3, &triplets);

    let v = vec![1.0, 2.0, 3.0];

    let res = numerical_sp_mat_vec_mul(&mat, &v).unwrap();

    assert_eq!(res, vec![7.0, 6.0, 0.0]);
}

#[test]

fn test_sparse_solve_cg() {

    let triplets = vec![
        (0, 0, 4.0),
        (0, 1, 1.0),
        (1, 0, 1.0),
        (1, 1, 3.0),
    ];

    let a = numerical_csr_from_triplets(2, 2, &triplets);

    let b = Array1::from_vec(vec![1.0, 2.0]);

    let res = numerical_solve_conjugate_gradient(&a, &b, None, 100, 1e-10).unwrap();

    assert_approx_eq!(res[0], 1.0 / 11.0, 1e-9);

    assert_approx_eq!(res[1], 7.0 / 11.0, 1e-9);
}

#[test]

fn test_sparse_trace_norm() {

    let triplets = vec![
        (0, 0, 1.0),
        (1, 1, -2.0),
        (2, 2, 3.0),
    ];

    let mat = numerical_csr_from_triplets(3, 3, &triplets);

    assert_eq!(numerical_trace(&mat).unwrap(), 2.0);

    assert_approx_eq!(numerical_frobenius_norm(&mat), 14.0f64.sqrt());

    assert_eq!(numerical_l1_norm(&mat), 3.0);

    assert_eq!(numerical_linf_norm(&mat), 3.0);
}

#[test]

fn test_sparse_predicates() {

    let triplets = vec![
        (0, 0, 1.0),
        (1, 1, 2.0),
    ];

    let mat = numerical_csr_from_triplets(2, 2, &triplets);

    assert!(numerical_is_symmetric(&mat, 1e-9));

    assert!(numerical_is_diagonal(&mat));

    let t2 = vec![(0, 1, 1.0)];

    let m2 = numerical_csr_from_triplets(2, 2, &t2);

    assert!(!numerical_is_symmetric(&m2, 1e-9));

    assert!(!numerical_is_diagonal(&m2));
}

#[test]

fn test_sparse_data_serde() {

    let triplets = vec![
        (0, 0, 1.0),
        (1, 2, 2.0),
    ];

    let mat = numerical_csr_from_triplets(3, 3, &triplets);

    let data = numerical_SparseMatrixData::from(&mat);

    let json = serde_json::to_string(&data).unwrap();

    let decoded: numerical_SparseMatrixData = serde_json::from_str(&json).unwrap();

    let mat_back = decoded.to_csmat();

    assert_eq!(mat_back.rows(), 3);

    assert_eq!(mat_back.cols(), 3);

    assert_eq!(mat_back.nnz(), 2);

    assert_eq!(mat_back.get(1, 2), Some(&2.0));
}
