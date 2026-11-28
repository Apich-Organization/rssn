//! Numerical reference checks for live-entry CSC norms.

use faer::sparse::SymbolicSparseColMat;
use rssn::kernels::sparse::SparseMat;
use rssn::kernels::sparse::SparseNorms;
use rssn::kernels::sparse::csr_from_triplets;
use rssn::kernels::sparse::frobenius_norm;
use rssn::kernels::sparse::l1_norm;
use rssn::kernels::sparse::linf_norm;
use rssn::kernels::sparse::matrix_norms;

fn assert_same_norms(
    matrix: &SparseMat,
    expected: SparseNorms,
) {
    let joint = matrix_norms(matrix);
    for (individual, joined, reference) in [
        (frobenius_norm(matrix), joint.frobenius, expected.frobenius),
        (l1_norm(matrix), joint.l1, expected.l1),
        (linf_norm(matrix), joint.linf, expected.linf),
    ] {
        if reference.is_nan() {
            assert!(individual.is_nan() && joined.is_nan());
        } else {
            assert_eq!(individual.to_bits(), reference.to_bits());
            assert_eq!(joined.to_bits(), reference.to_bits());
        }
    }
}

#[test]
fn rectangular_norms_match_an_independent_dense_reference() {
    let triplets = [
        (1, 0, 1.25),
        (0, 3, -3.0),
        (2, 1, 1.5),
        (2, 4, -2.0),
        (1, 0, 0.25),
        (0, 1, 5.0),
        (2, 1, -0.5),
    ];
    let matrix = csr_from_triplets(3, 5, &triplets);
    let mut dense = [[0.0_f64; 5]; 3];
    for (row, col, value) in triplets {
        dense[row][col] += value;
    }
    let l1 = (0..5)
        .map(|col| dense.iter().map(|row| row[col].abs()).sum::<f64>())
        .fold(0.0, f64::max);
    let linf = dense
        .iter()
        .map(|row| row.iter().map(|value| value.abs()).sum::<f64>())
        .fold(0.0, f64::max);
    let frobenius = dense
        .iter()
        .flatten()
        .map(|value| value.powi(2))
        .sum::<f64>()
        .sqrt();
    assert_same_norms(&matrix, SparseNorms { frobenius, l1, linf });
}

#[test]
fn padding_and_values_outside_the_live_column_offsets_are_ignored() {
    let pad = usize::MAX;
    let symbolic = SymbolicSparseColMat::new_checked(
        3,
        4,
        vec![1, 4, 6, 8, 9],
        Some(vec![2, 0, 1, 0]),
        vec![pad, 0, 2, pad, pad, 0, 1, pad, pad, pad],
    );
    let matrix = SparseMat::new(
        symbolic,
        vec![
            f64::NAN,
            3.0,
            -12.0,
            f64::MAX,
            f64::NAN,
            99.0,
            -4.0,
            f64::INFINITY,
            -999.0,
            f64::NAN,
        ],
    );
    assert_eq!(matrix.compute_nnz(), 3);
    assert_same_norms(
        &matrix,
        SparseNorms {
            frobenius: 13.0,
            l1: 15.0,
            linf: 12.0,
        },
    );
}

#[test]
fn the_contiguous_fast_path_respects_nonzero_offsets() {
    let symbolic = SymbolicSparseColMat::new_checked(
        2,
        2,
        vec![1, 3, 4],
        None,
        vec![usize::MAX, 0, 1, 0, usize::MAX],
    );
    let matrix = SparseMat::new(symbolic, vec![f64::NAN, 3.0, -4.0, 12.0, f64::INFINITY]);
    assert_same_norms(
        &matrix,
        SparseNorms {
            frobenius: 13.0,
            l1: 12.0,
            linf: 15.0,
        },
    );
}

#[test]
fn scaled_frobenius_recovers_a_subnormal_and_an_extremely_large_finite_norm() {
    let microscopic = csr_from_triplets(
        2,
        1,
        &[(0, 0, f64::from_bits(3)), (1, 0, -f64::from_bits(4))],
    );
    assert_same_norms(
        &microscopic,
        SparseNorms {
            frobenius: f64::from_bits(5),
            l1: f64::from_bits(7),
            linf: f64::from_bits(4),
        },
    );

    let scale = 2.0_f64.powi(700);
    let huge = csr_from_triplets(2, 1, &[(0, 0, 3.0 * scale), (1, 0, -4.0 * scale)]);
    assert_same_norms(
        &huge,
        SparseNorms {
            frobenius: 5.0 * scale,
            l1: 7.0 * scale,
            linf: 4.0 * scale,
        },
    );

    let huge = csr_from_triplets(2, 1, &[(0, 0, f64::MAX), (1, 0, -f64::MAX)]);
    assert_same_norms(
        &huge,
        SparseNorms {
            frobenius: f64::INFINITY,
            l1: f64::INFINITY,
            linf: f64::MAX,
        },
    );
}

#[test]
fn half_scale_recovery_handles_both_sides_of_the_overflow_threshold() {
    let ulp = 2.0_f64.powi(971);
    let cases = [
        (
            [f64::MAX, ulp * 0.25, ulp * 0.25, ulp * 0.25],
            f64::INFINITY,
        ),
        (
            [f64::MAX - 2.0 * ulp, 0.625 * ulp, 0.625 * ulp, 0.625 * ulp],
            f64::MAX,
        ),
    ];
    for (entries, expected) in cases {
        let vertical: Vec<_> = entries
            .into_iter()
            .enumerate()
            .map(|(i, value)| (i, 0, value))
            .collect();
        let horizontal: Vec<_> = entries
            .into_iter()
            .enumerate()
            .map(|(i, value)| (0, i, -value))
            .collect();
        let vertical = csr_from_triplets(4, 1, &vertical);
        let horizontal = csr_from_triplets(1, 4, &horizontal);
        let max_entry = entries[0];
        assert_same_norms(
            &vertical,
            SparseNorms {
                frobenius: max_entry,
                l1: expected,
                linf: max_entry,
            },
        );
        assert_same_norms(
            &horizontal,
            SparseNorms {
                frobenius: max_entry,
                l1: max_entry,
                linf: expected,
            },
        );
    }
}

#[test]
fn empty_axes_and_explicit_zeros_have_zero_norms() {
    for (rows, cols, entries) in [
        (0, 0, &[][..]),
        (1_000_000, 0, &[][..]),
        (0, 7, &[][..]),
        (3, 5, &[(0, 1, -0.0)][..]),
    ] {
        assert_same_norms(
            &csr_from_triplets(rows, cols, entries),
            SparseNorms::default(),
        );
    }
}

#[test]
fn live_nan_and_infinity_follow_the_existing_norm_conventions() {
    let nan = csr_from_triplets(3, 2, &[(0, 0, f64::NAN), (1, 0, -4.0), (2, 1, 3.0)]);
    assert_same_norms(
        &nan,
        SparseNorms {
            frobenius: f64::NAN,
            l1: 3.0,
            linf: 4.0,
        },
    );

    let inf = csr_from_triplets(2, 2, &[(0, 0, f64::NEG_INFINITY), (1, 1, 1.0)]);
    assert_same_norms(
        &inf,
        SparseNorms {
            frobenius: f64::INFINITY,
            l1: f64::INFINITY,
            linf: f64::INFINITY,
        },
    );
}
