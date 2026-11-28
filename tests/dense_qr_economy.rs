//! Economy QR retains orthogonality, shape, scale, and reproducibility.

use rssn::kernels::dense::DenseError;
use rssn::kernels::dense::Mat;
use rssn::kernels::dense::qr;
use rssn::kernels::dense::qr_thin;

fn assert_qr(
    matrix: &Mat,
    q: &Mat,
    r: &Mat,
) {
    let product = q.matmul(r).unwrap();
    for j in 0..matrix.cols {
        let scale = (0..matrix.rows)
            .map(|i| matrix.at(i, j).abs())
            .fold(0.0_f64, f64::max);
        let scale = if scale > 0.0 { scale } else { 1.0 };
        for i in 0..matrix.rows {
            assert!(
                (product.at(i, j) / scale - matrix.at(i, j) / scale).abs() < 2e-13,
                "column {j}, row {i}, scale {scale}",
            );
        }
        for i in 0..r.rows {
            if i > j {
                assert!((r.at(i, j) / scale).abs() < 2e-13);
            }
        }
    }
    let gram = q.transpose().matmul(q).unwrap();
    for i in 0..gram.rows {
        for j in 0..gram.cols {
            assert!((gram.at(i, j) - f64::from(i == j)).abs() < 2e-13);
        }
    }
}

#[test]
fn tall_square_wide_and_rank_deficient_factors_reconstruct_input() {
    for (m, n) in [(0, 4), (4, 0), (1, 5), (5, 1), (7, 3), (5, 5), (3, 7)] {
        let mut matrix = Mat::zeros(m, n);
        for i in 0..m {
            for j in 0..n {
                let j_as_float = f64::from(u32::try_from(j).unwrap());
                let input = f64::from(u32::try_from(i + j % 2).unwrap());
                matrix.set(i, j, input.sin() * (j_as_float + 1.0));
            }
        }
        let (q, r) = qr_thin(&matrix).unwrap();
        assert_eq!((q.rows, q.cols, r.rows, r.cols), (m, m.min(n), m.min(n), n));
        assert_qr(&matrix, &q, &r);
        let (q, r) = qr(&matrix);
        assert_eq!((q.rows, q.cols, r.rows, r.cols), (m, m, m, n));
        assert_qr(&matrix, &q, &r);
    }
}

#[test]
fn independent_small_columns_survive_beside_huge_columns() {
    let matrix = Mat::from_rows(&[
        vec![3e250, 0.0, 1e-250],
        vec![4e250, 2e-280, 0.0],
        vec![0.0, 3e-280, 3e-250],
        vec![1e250, 0.0, 2e-250],
    ])
    .unwrap();
    let (q, r) = qr_thin(&matrix).unwrap();
    assert_qr(&matrix, &q, &r);
    assert!(r.at(1, 1).abs() > 1e-280);
    let (q, r) = qr(&matrix);
    assert_qr(&matrix, &q, &r);
}

#[test]
fn checked_qr_reports_nonfinite_input_and_unrepresentable_factors() {
    for rows in [
        vec![vec![1.0], vec![f64::NAN]],
        vec![vec![f64::MAX], vec![f64::MAX]],
    ] {
        let matrix = Mat::from_rows(&rows).unwrap();
        assert!(matches!(qr_thin(&matrix), Err(DenseError::NonFinite)));
    }
}
