//! CSR sparse matrices (ported from `numerical_sparse_test.rs`).

use assert_approx_eq::assert_approx_eq;
use ndarray::{Array1, array};
use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::kernels::sparse::{
    SparseMatrixData, csr_from_triplets, frobenius_norm, is_diagonal, is_symmetric, l1_norm,
    linf_norm, rank, solve_conjugate_gradient, sp_mat_vec_mul, to_csr, to_dense, trace, transpose,
};

fn cfg() -> ProptestConfig {
    ProptestConfig {
        rng_seed: RngSeed::Fixed(0x5EED),
        failure_persistence: None,
        ..ProptestConfig::default()
    }
}

#[test]
fn construction() {
    let mat = csr_from_triplets(3, 3, &[(0, 0, 1.0), (1, 1, 2.0), (2, 2, 3.0)]);
    assert_eq!((mat.rows(), mat.cols(), mat.nnz()), (3, 3, 3));
    assert_eq!(mat.get(1, 1), Some(&2.0));
    assert_eq!(mat.get(0, 1), None);
}

#[test]
fn matrix_vector_product() {
    let mat = csr_from_triplets(3, 3, &[(0, 0, 1.0), (0, 2, 2.0), (1, 1, 3.0)]);
    assert_eq!(sp_mat_vec_mul(&mat, &[1.0, 2.0, 3.0]), Ok(vec![7.0, 6.0, 0.0]));
    assert!(sp_mat_vec_mul(&mat, &[1.0]).is_err());
}

#[test]
fn conjugate_gradient_solves_spd_system() {
    let a = csr_from_triplets(2, 2, &[(0, 0, 4.0), (0, 1, 1.0), (1, 0, 1.0), (1, 1, 3.0)]);
    let b = Array1::from_vec(vec![1.0, 2.0]);
    let x = solve_conjugate_gradient(&a, &b, None, 100, 1e-10).unwrap_or_else(|e| panic!("{e}"));
    assert_approx_eq!(x[0], 1.0 / 11.0, 1e-9);
    assert_approx_eq!(x[1], 7.0 / 11.0, 1e-9);
}

#[test]
fn conjugate_gradient_rejects_bad_dimensions() {
    let a = csr_from_triplets(2, 2, &[(0, 0, 1.0), (1, 1, 1.0)]);
    assert!(solve_conjugate_gradient(&a, &Array1::from_vec(vec![1.0]), None, 10, 1e-9).is_err());
}

#[test]
fn trace_and_norms() {
    let mat = csr_from_triplets(3, 3, &[(0, 0, 1.0), (1, 1, -2.0), (2, 2, 3.0)]);
    assert_eq!(trace(&mat), Ok(2.0));
    assert_approx_eq!(frobenius_norm(&mat), 14f64.sqrt());
    assert_eq!(l1_norm(&mat), 3.0);
    assert_eq!(linf_norm(&mat), 3.0);
    assert!(trace(&csr_from_triplets(2, 3, &[])).is_err());
}

#[test]
fn predicates() {
    let d = csr_from_triplets(2, 2, &[(0, 0, 1.0), (1, 1, 2.0)]);
    assert!(is_symmetric(&d, 1e-9) && is_diagonal(&d));
    let u = csr_from_triplets(2, 2, &[(0, 1, 1.0)]);
    assert!(!is_symmetric(&u, 1e-9) && !is_diagonal(&u));
    assert!(!is_symmetric(&csr_from_triplets(2, 3, &[]), 1e-9));
}

#[test]
fn dense_round_trip_and_rank() {
    let dense = array![[1.0, 0.0, 2.0], [0.0, 0.0, 0.0], [0.0, 3.0, 0.0]].into_dyn();
    let sp = to_csr(&dense);
    assert_eq!(sp.nnz(), 3);
    let back = to_dense(&sp);
    for i in 0..3 {
        for j in 0..3 {
            assert_eq!(back[[i, j]], dense[[i, j]]);
        }
    }
    assert_eq!(rank(&sp), 2);
}

#[test]
fn transpose_moves_entries() {
    let mat = csr_from_triplets(2, 3, &[(0, 2, 5.0), (1, 0, 7.0)]);
    let t = transpose(&mat);
    assert_eq!((t.rows(), t.cols()), (3, 2));
    assert_eq!(t.get(2, 0), Some(&5.0));
    assert_eq!(t.get(0, 1), Some(&7.0));
}

#[test]
fn serde_round_trip() {
    let mat = csr_from_triplets(3, 3, &[(0, 0, 1.0), (1, 2, 2.0)]);
    let json = serde_json::to_string(&SparseMatrixData::from(&mat)).unwrap_or_else(|e| panic!("{e}"));
    let decoded: SparseMatrixData = serde_json::from_str(&json).unwrap_or_else(|e| panic!("{e}"));
    let back = decoded.to_csmat();
    assert_eq!((back.rows(), back.cols(), back.nnz()), (3, 3, 2));
    assert_eq!(back.get(1, 2), Some(&2.0));
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn prop_transpose_is_an_involution(
        rows in 1..20usize, cols in 1..20usize,
        elements in prop::collection::vec((0..20usize, 0..20usize, -100.0..100.0f64), 0..20),
    ) {
        let triplets: Vec<_> = elements.into_iter().map(|(r, c, v)| (r % rows, c % cols, v)).collect();
        let mat = csr_from_triplets(rows, cols, &triplets);
        let tt = transpose(&transpose(&mat));
        prop_assert_eq!((mat.rows(), mat.cols(), mat.nnz()), (tt.rows(), tt.cols(), tt.nnz()));
        for (val, (r, c)) in mat.iter() {
            prop_assert!((val - tt.get(r, c).copied().unwrap_or(0.0)).abs() < 1e-9);
        }
    }

    #[test]
    fn prop_trace_is_transpose_invariant(
        size in 1..20usize,
        elements in prop::collection::vec((0..20usize, 0..20usize, -100.0..100.0f64), 0..20),
    ) {
        let triplets: Vec<_> = elements.into_iter().map(|(r, c, v)| (r % size, c % size, v)).collect();
        let mat = csr_from_triplets(size, size, &triplets);
        let t1 = trace(&mat).map_err(TestCaseError::fail)?;
        let t2 = trace(&transpose(&mat)).map_err(TestCaseError::fail)?;
        prop_assert!((t1 - t2).abs() < 1e-9);
    }

    #[test]
    fn prop_spmv_matches_dense(
        n in 1..8usize,
        elements in prop::collection::vec((0..8usize, 0..8usize, -10.0..10.0f64), 0..20),
        v in prop::collection::vec(-10.0..10.0f64, 8),
    ) {
        let triplets: Vec<_> = elements.into_iter().map(|(r, c, x)| (r % n, c % n, x)).collect();
        let mat = csr_from_triplets(n, n, &triplets);
        let dense = to_dense(&mat);
        let got = sp_mat_vec_mul(&mat, &v[..n]).map_err(TestCaseError::fail)?;
        for i in 0..n {
            let want: f64 = (0..n).map(|j| dense[[i, j]] * v[j]).sum();
            prop_assert!((got[i] - want).abs() < 1e-9);
        }
    }

    #[test]
    fn prop_cg_solves_diagonally_dominant_spd_systems(n in 2..8usize, seed in 0..1000u64) {
        // Tridiagonal SPD matrix with a planted solution.
        let mut t = vec![];
        for i in 0..n {
            t.push((i, i, 4.0));
            if i + 1 < n { t.push((i, i + 1, -1.0)); t.push((i + 1, i, -1.0)); }
        }
        let a = csr_from_triplets(n, n, &t);
        let x_true: Vec<f64> = (0..n).map(|i| ((i as u64 * 7 + seed) % 13) as f64 - 6.0).collect();
        let b = Array1::from_vec(sp_mat_vec_mul(&a, &x_true).map_err(TestCaseError::fail)?);
        let x = solve_conjugate_gradient(&a, &b, None, 200, 1e-12).map_err(TestCaseError::fail)?;
        for (g, w) in x.iter().zip(&x_true) {
            prop_assert!((g - w).abs() < 1e-7, "{g} vs {w}");
        }
    }
}
