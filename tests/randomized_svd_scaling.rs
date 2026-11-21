//! Matrix-free randomized SVD against a known rotated rectangular spectrum.

use rssn::kernels::dense::Mat;
use rssn::kernels::dense::Svd;
use rssn::kernels::lowrank::LowRankError;
use rssn::kernels::lowrank::RandomizedOptions;
use rssn::kernels::lowrank::randomized_svd;

const ROWS: usize = 9;
const COLS: usize = 7;
const SINGULAR: [f64; 3] = [7.0, 3.0, 1.0];

fn known_matrix() -> Mat {
    // U and V have three orthonormal columns with disjoint two-row support.
    let mut matrix = Mat::zeros(ROWS, COLS);
    for (block, singular) in SINGULAR.into_iter().enumerate() {
        for (i, u) in [3.0 / 5.0, -4.0 / 5.0].into_iter().enumerate() {
            for (j, v) in [5.0 / 13.0, 12.0 / 13.0].into_iter().enumerate() {
                matrix.set(2 * block + i, 2 * block + j, u * singular * v);
            }
        }
    }
    matrix
}

fn compute(
    matrix: &Mat,
    scale: f64,
    rank: usize,
    options: &RandomizedOptions,
) -> Svd {
    randomized_svd(
        &|x, y| {
            for (out, row) in y.iter_mut().zip(matrix.data.chunks_exact(COLS)) {
                *out = scale * row.iter().zip(x).map(|(a, x)| a * x).sum::<f64>();
            }
        },
        &|x, y| {
            for (column, out) in y.iter_mut().enumerate() {
                *out = scale
                    * x.iter()
                        .enumerate()
                        .map(|(i, x)| x * matrix.at(i, column))
                        .sum::<f64>();
            }
        },
        ROWS,
        COLS,
        rank,
        options,
    )
    .unwrap()
}

fn assert_orthogonal(factor: &Mat) {
    for i in 0..factor.cols {
        for j in 0..factor.cols {
            let inner: f64 = factor
                .data
                .chunks_exact(factor.cols)
                .map(|row| row[i] * row[j])
                .sum();
            let exact = f64::from(u8::from(i == j));
            assert!((inner - exact).abs() < 5e-14);
        }
    }
}

#[test]
fn preserves_rotated_singular_values_at_both_scale_extremes() {
    let matrix = known_matrix();
    for exponent in [-600, 0, 600] {
        let scale = 2.0_f64.powi(exponent);
        for power_iterations in [0, 3] {
            let options = RandomizedOptions {
                oversample: 4,
                power_iterations,
                ..RandomizedOptions::default()
            };
            let got = compute(&matrix, scale, 3, &options);
            assert_eq!(got.s.len(), 3);
            for (&actual, expected) in got.s.iter().zip(SINGULAR) {
                assert!((actual / scale - expected).abs() < 5e-13);
            }
            assert_eq!((got.u.rows, got.u.cols), (ROWS, 3));
            assert_eq!((got.v.rows, got.v.cols), (COLS, 3));
            assert_orthogonal(&got.u);
            assert_orthogonal(&got.v);
            for i in 0..ROWS {
                for j in 0..COLS {
                    let reconstructed: f64 = (0..3)
                        .map(|k| got.u.at(i, k) * (got.s[k] / scale) * got.v.at(j, k))
                        .sum();
                    assert!((reconstructed - matrix.at(i, j)).abs() < 5e-13);
                }
            }
        }
    }
}

#[test]
fn rank_deficient_range_keeps_its_rank_and_seed_determinism() {
    let matrix = known_matrix();
    let options = RandomizedOptions {
        oversample: usize::MAX,
        ..RandomizedOptions::default()
    };
    let large = 2.0_f64.powi(600);
    let first = compute(&matrix, large, 4, &options);
    let second = compute(&matrix, large, 4, &options);
    assert_eq!(first.s.len(), 3);
    assert_eq!(first.s, second.s);
    assert_eq!(first.u.data, second.u.data);
    assert_eq!(first.v.data, second.v.data);
}

#[test]
fn still_rejects_no_range_and_nonfinite_adjoint_products() {
    let options = RandomizedOptions {
        oversample: 0,
        power_iterations: 0,
        ..RandomizedOptions::default()
    };
    let identity = |x: &[f64], y: &mut [f64]| y.copy_from_slice(x);
    let zero = |_: &[f64], y: &mut [f64]| y.fill(0.0);
    assert!(matches!(
        randomized_svd(&zero, &zero, 4, 4, 2, &options),
        Err(LowRankError::Shape(_))
    ));
    let invalid_adjoint = |x: &[f64], y: &mut [f64]| {
        y.copy_from_slice(x);
        y[1] = f64::NAN;
    };
    assert!(matches!(
        randomized_svd(&identity, &invalid_adjoint, 4, 4, 2, &options),
        Err(LowRankError::Shape(_))
    ));
    assert!(matches!(
        randomized_svd(
            &identity,
            &invalid_adjoint,
            4,
            4,
            2,
            &RandomizedOptions {
                power_iterations: 1,
                ..options
            }
        ),
        Err(LowRankError::Shape(_))
    ));
}
