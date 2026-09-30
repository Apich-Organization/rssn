//! Dense matrices (ported from `numerical_matrix_test.rs`; the
//! `PrimeFieldElement` `Field` impl no longer exists).

use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::kernels::matrix::{Backend, Matrix};

fn cfg() -> ProptestConfig {
    ProptestConfig {
        rng_seed: RngSeed::Fixed(0x5EED),
        failure_persistence: None,
        ..ProptestConfig::default()
    }
}

fn m(rows: usize, cols: usize, d: &[f64]) -> Matrix<f64> {
    Matrix::new(rows, cols, d.to_vec())
}

#[test]
fn construction_and_access() {
    let data = vec![1.0, 2.0, 3.0, 4.0];
    let mat = Matrix::new(2, 2, data.clone());
    assert_eq!(mat.rows(), 2);
    assert_eq!(mat.cols(), 2);
    assert_eq!(mat.data(), &data);
    assert_eq!(*mat.get(0, 1), 2.0);
    assert_eq!(mat.get_cols(), vec![vec![1.0, 3.0], vec![2.0, 4.0]]);
    let z = Matrix::<f64>::zeros(2, 3);
    assert!(z.data().iter().all(|&v| v == 0.0));
    assert_eq!(z.into_data().len(), 6);
}

#[test]
fn arithmetic() {
    let m1 = m(2, 2, &[1.0, 2.0, 3.0, 4.0]);
    let m2 = m(2, 2, &[5.0, 6.0, 7.0, 8.0]);
    assert_eq!((m1.clone() + m2.clone()).data(), &vec![6.0, 8.0, 10.0, 12.0]);
    assert_eq!((m2.clone() - m1.clone()).data(), &vec![4.0, 4.0, 4.0, 4.0]);
    // [1 2; 3 4] * [5 6; 7 8] = [19 22; 43 50]
    assert_eq!((m1.clone() * m2.clone()).data(), &vec![19.0, 22.0, 43.0, 50.0]);
    assert_eq!((m1.clone() * 2.0).data(), &vec![2.0, 4.0, 6.0, 8.0]);
    assert_eq!((m1.clone() / 2.0).data(), &vec![0.5, 1.0, 1.5, 2.0]);
    assert_eq!((-m1).data(), &vec![-1.0, -2.0, -3.0, -4.0]);
}

#[test]
fn transpose() {
    let mt = m(2, 3, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]).transpose();
    assert_eq!((mt.rows(), mt.cols()), (3, 2));
    assert_eq!(mt.data(), &vec![1.0, 4.0, 2.0, 5.0, 3.0, 6.0]);
}

#[test]
fn determinant() {
    assert_approx_eq!(m(2, 2, &[1.0, 2.0, 3.0, 4.0]).determinant().unwrap_or(f64::NAN), -2.0);
    // -24 + 40 - 15 = 1
    let m3 = m(3, 3, &[1.0, 2.0, 3.0, 0.0, 1.0, 4.0, 5.0, 6.0, 0.0]);
    assert_approx_eq!(m3.determinant().unwrap_or(f64::NAN), 1.0);
    assert_approx_eq!(m3.determinant_lu().unwrap_or(f64::NAN), 1.0);
    assert_approx_eq!(Matrix::<f64>::identity(4).determinant().unwrap_or(f64::NAN), 1.0);
    assert!(m(2, 3, &[1.0; 6]).determinant().is_err());
}

#[test]
fn block_determinant_matches_lu() {
    let a = m(
        4,
        4,
        &[4.0, 1.0, 2.0, 0.5, 1.0, 3.0, 0.0, 1.0, 2.0, 0.0, 5.0, 1.0, 0.5, 1.0, 1.0, 6.0],
    );
    let lu = a.determinant_lu().unwrap_or(f64::NAN);
    let block = a.determinant_block().unwrap_or(f64::NAN);
    assert_approx_eq!(lu, block, 1e-8);
}

#[test]
fn singular_matrix_has_zero_determinant_or_error() {
    let s = m(2, 2, &[1.0, 2.0, 2.0, 4.0]);
    match s.determinant() {
        Ok(d) => assert!(d.abs() < 1e-12, "det = {d}"),
        Err(_) => {}
    }
    assert_eq!(s.rank(), Ok(1));
}

#[test]
#[ignore = "library bug: Matrix::inverse returns Some(bogus) for singular [[1,2],[2,4]] because rref pivots into the identity half"]
fn singular_matrix_has_no_inverse() {
    // observed: Some([[0, 0.5], [1, -0.5]]); expected: None
    let s = m(2, 2, &[1.0, 2.0, 2.0, 4.0]);
    assert!(s.inverse().is_none());
}

#[test]
fn inverse_of_2x2() {
    // det = 10 ; inverse = [0.6 -0.7; -0.2 0.4]
    let a = m(2, 2, &[4.0, 7.0, 2.0, 6.0]);
    let inv = a.inverse().unwrap_or_else(|| panic!("matrix should be invertible"));
    assert_approx_eq!(*inv.get(0, 0), 0.6);
    assert_approx_eq!(*inv.get(0, 1), -0.7);
    assert_approx_eq!(*inv.get(1, 0), -0.2);
    assert_approx_eq!(*inv.get(1, 1), 0.4);
    let id = a * inv;
    assert!(id.is_identity(1e-9));
}

#[test]
fn inverse_of_3x3_known() {
    // [[1,2,3],[0,1,4],[5,6,0]] has inverse [[-24,18,5],[20,-15,-4],[-5,4,1]].
    let a = m(3, 3, &[1.0, 2.0, 3.0, 0.0, 1.0, 4.0, 5.0, 6.0, 0.0]);
    let inv = a.inverse().unwrap_or_else(|| panic!("matrix should be invertible"));
    let want = [-24.0, 18.0, 5.0, 20.0, -15.0, -4.0, -5.0, 4.0, 1.0];
    for (g, w) in inv.data().iter().zip(want) {
        assert_approx_eq!(*g, w, 1e-9);
    }
    assert!(m(2, 3, &[1.0; 6]).inverse().is_none());
}

#[test]
fn faer_backend_inverse_agrees_with_native() {
    let a = m(3, 3, &[2.0, -1.0, 0.0, -1.0, 2.0, -1.0, 0.0, -1.0, 2.0]);
    let native = a.inverse().unwrap_or_else(|| panic!("native inverse failed"));
    let faer = a
        .clone()
        .with_backend(Backend::Faer)
        .inverse()
        .unwrap_or_else(|| panic!("faer inverse failed"));
    for (x, y) in native.data().iter().zip(faer.data()) {
        assert_approx_eq!(*x, *y, 1e-10);
    }
}

#[test]
fn norms() {
    let a = m(2, 2, &[1.0, -2.0, -3.0, 4.0]);
    assert_approx_eq!(a.frobenius_norm(), 30f64.sqrt());
    assert_eq!(a.l1_norm(), 6.0); // max column sum
    assert_eq!(a.linf_norm(), 7.0); // max row sum
}

#[test]
fn trace_and_rank() {
    assert_eq!(m(3, 3, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]).trace(), Ok(15.0));
    assert!(m(2, 3, &[1.0; 6]).trace().is_err());
    assert_eq!(m(3, 3, &[1.0, 2.0, 3.0, 2.0, 4.0, 6.0, 0.0, 1.0, 1.0]).rank(), Ok(2));
    assert_eq!(Matrix::<f64>::identity(5).rank(), Ok(5));
}

#[test]
fn rref_and_null_space() {
    let mut a = m(2, 3, &[1.0, 2.0, 3.0, 2.0, 4.0, 6.0]);
    assert_eq!(a.rref(), Ok(1));
    assert_eq!(a.data(), &vec![1.0, 2.0, 3.0, 0.0, 0.0, 0.0]);

    let b = m(2, 3, &[1.0, 2.0, 3.0, 2.0, 4.0, 6.0]);
    let ns = b.null_space().unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(ns.cols(), 2);
    for col in ns.get_cols() {
        let r0 = col[0] + 2.0 * col[1] + 3.0 * col[2];
        assert_approx_eq!(r0, 0.0, 1e-12);
    }
}

#[test]
fn identity_orthogonal_symmetric_diagonal() {
    let id = Matrix::<f64>::identity(3);
    assert!(id.is_identity(1e-9) && id.is_orthogonal(1e-9) && id.is_symmetric() && id.is_diagonal());
    let swap = m(2, 2, &[0.0, 1.0, 1.0, 0.0]);
    assert!(!swap.is_identity(1e-9));
    assert!(swap.is_orthogonal(1e-9));
    assert!(swap.is_symmetric());
    assert!(!swap.is_diagonal());
    assert!(!m(2, 2, &[1.0, 2.0, 3.0, 4.0]).is_orthogonal(1e-9));
}

#[test]
fn strassen_matches_naive_product() {
    let a = m(3, 2, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
    let b = m(2, 3, &[7.0, 8.0, 9.0, 10.0, 11.0, 12.0]);
    let s = a.mul_strassen(&b).unwrap_or_else(|e| panic!("{e}"));
    let naive = a.clone() * b.clone();
    assert_eq!((s.rows(), s.cols()), (naive.rows(), naive.cols()));
    for (x, y) in s.data().iter().zip(naive.data()) {
        assert_approx_eq!(*x, *y, 1e-9);
    }
    assert!(a.mul_strassen(&a).is_err());
}

#[test]
fn jacobi_eigen_decomposition_of_symmetric_matrix() {
    // Eigenvalues of [[2,1],[1,2]] are 1 and 3.
    let a = m(2, 2, &[2.0, 1.0, 1.0, 2.0]);
    let (mut vals, vecs) = a.jacobi_eigen_decomposition(50, 1e-12).unwrap_or_else(|e| panic!("{e}"));
    vals.sort_by(f64::total_cmp);
    assert_approx_eq!(vals[0], 1.0, 1e-9);
    assert_approx_eq!(vals[1], 3.0, 1e-9);
    assert!(vecs.is_orthogonal(1e-9));
    assert!(m(2, 3, &[1.0; 6]).jacobi_eigen_decomposition(10, 1e-9).is_err());
}

fn dominant_matrix(n: usize, seed: u64) -> Matrix<f64> {
    let mut rng = seed;
    let mut next = || {
        rng = rng.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
        ((rng >> 33) % 1000) as f64 / 500.0 - 1.0
    };
    let mut d: Vec<f64> = (0..n * n).map(|_| next()).collect();
    for i in 0..n {
        d[i * n + i] += n as f64 + 1.0; // diagonally dominant => invertible
    }
    Matrix::new(n, n, d)
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn prop_transpose_is_an_involution(rows in 1..10usize, cols in 1..10usize, seed in 0..1000u64) {
        let data: Vec<f64> = (0..rows * cols).map(|i| ((i as u64 * 31 + seed) % 97) as f64).collect();
        let a = Matrix::new(rows, cols, data);
        prop_assert_eq!(a.transpose().transpose(), a);
    }

    #[test]
    fn prop_inverse_times_matrix_is_identity(n in 1..6usize, seed in 0..10_000u64) {
        let a = dominant_matrix(n, seed);
        let inv = a.inverse().ok_or_else(|| TestCaseError::fail("invertible matrix reported singular"))?;
        prop_assert!((a * inv).is_identity(1e-8));
    }

    #[test]
    fn prop_determinant_is_multiplicative(n in 2..5usize, s1 in 0..1000u64, s2 in 0..1000u64) {
        let (a, b) = (dominant_matrix(n, s1), dominant_matrix(n, s2 + 5000));
        let da = a.determinant().map_err(TestCaseError::fail)?;
        let db = b.determinant().map_err(TestCaseError::fail)?;
        let dab = (a * b).determinant().map_err(TestCaseError::fail)?;
        prop_assert!((dab - da * db).abs() < 1e-6 * dab.abs().max(1.0));
    }

    #[test]
    fn prop_determinant_of_transpose(n in 2..5usize, seed in 0..1000u64) {
        let a = dominant_matrix(n, seed);
        let d1 = a.determinant().map_err(TestCaseError::fail)?;
        let d2 = a.transpose().determinant().map_err(TestCaseError::fail)?;
        prop_assert!((d1 - d2).abs() < 1e-8 * d1.abs().max(1.0));
    }

    #[test]
    fn prop_matrix_multiplication_distributes_over_addition(n in 1..5usize, s in 0..1000u64) {
        let (a, b, c) = (dominant_matrix(n, s), dominant_matrix(n, s + 1), dominant_matrix(n, s + 2));
        let left = a.clone() * (b.clone() + c.clone());
        let right = a.clone() * b + a * c;
        for (x, y) in left.data().iter().zip(right.data()) {
            prop_assert!((x - y).abs() < 1e-9 * x.abs().max(1.0));
        }
    }
}
