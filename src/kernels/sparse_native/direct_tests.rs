#![allow(
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    clippy::unwrap_used,
    reason = "small fixture arrays and tolerance comparisons"
)]

use super::Csr;
use super::NumericError;
use super::cholesky_solve;
use super::cholesky_solve_many;
use super::lu_solve;
use super::lu_solve_many;
use super::solve_auto_many;

fn check_values(
    actual: &[f64],
    expected: &[f64],
) {
    assert_eq!(actual.len(), expected.len());
    for (&actual, &expected) in actual.iter().zip(expected) {
        assert!((actual - expected).abs() <= 1e-13 * expected.abs().max(1.0));
    }
}

#[test]
fn one_factor_solves_multiple_rhs_for_unsorted_duplicate_csr() {
    let a = Csr {
        rows: 3,
        cols: 3,
        row_ptr: &[0, 4, 8, 11],
        col_idx: &[1, 0, 1, 0, 2, 1, 0, 2, 2, 1, 2],
        values: &[0.25, 1.5, 0.75, 2.5, 0.25, 3.0, 1.0, 0.25, 0.5, 0.5, 1.5],
    };
    let x = [1.0, 2.0, -1.0, -3.0, 4.0, 0.0, 2.0, 0.5, -2.0];
    let b = [1.0, 12.0, -4.0, -7.0, 14.25, -2.0, 2.5, 3.0, -4.0];
    let lu = lu_solve_many(&a, &b, 3).unwrap();
    let cholesky = cholesky_solve_many(&a, &b, 3).unwrap();
    check_values(&lu, &x);
    check_values(&cholesky, &x);
    let (auto, method) = solve_auto_many(&a, &b, 3).unwrap();
    assert_eq!(method, "Cholesky");
    check_values(&auto, &x);
    for col in 0..3 {
        let bc = [b[col], b[col + 3], b[col + 6]];
        let lc = [lu[col], lu[col + 3], lu[col + 6]];
        let cc = [cholesky[col], cholesky[col + 3], cholesky[col + 6]];
        check_values(&lu_solve(&a, &bc).unwrap(), &lc);
        check_values(&cholesky_solve(&a, &bc).unwrap(), &cc);
    }
}

#[test]
fn duplicate_entries_cannot_hide_a_nonsymmetric_matrix_from_auto() {
    let a = Csr {
        rows: 2,
        cols: 2,
        row_ptr: &[0, 2, 5],
        col_idx: &[0, 1, 0, 0, 1],
        values: &[2.0, 1.0, 1.0, 1.0, 3.0],
    };
    let b = [3.0, 1.0, 5.0, -1.0];
    assert_eq!(
        cholesky_solve_many(&a, &b, 2),
        Err(NumericError::NotPositiveDefinite)
    );
    let (solution, algorithm) = solve_auto_many(&a, &b, 2).unwrap();
    assert_eq!(algorithm, "LU");
    check_values(&solution, &[1.0, 1.0, 1.0, -1.0]);
}

#[test]
fn automatic_lu_after_cholesky_breakdown_preserves_the_rhs() {
    let a = Csr {
        rows: 2,
        cols: 2,
        row_ptr: &[0, 2, 4],
        col_idx: &[0, 1, 0, 1],
        values: &[1.0, 2.0, 2.0, 1.0],
    };
    let b = [7.0, 4.0, 5.0, 5.0];
    let (solution, algorithm) = solve_auto_many(&a, &b, 2).unwrap();
    assert_eq!(algorithm, "LU");
    check_values(&solution, &[1.0, 2.0, 3.0, 1.0]);
}

#[test]
fn dense_rhs_shapes_and_empty_matrices_are_checked_before_factoring() {
    let empty = Csr {
        rows: 0,
        cols: 0,
        row_ptr: &[0],
        col_idx: &[],
        values: &[],
    };
    assert!(lu_solve_many(&empty, &[], 5).unwrap().is_empty());
    assert_eq!(solve_auto_many(&empty, &[], 5), Ok((Vec::new(), "LU")));
    assert_eq!(
        cholesky_solve_many(&empty, &[], 5),
        Err(NumericError::NotPositiveDefinite)
    );
    assert_eq!(
        lu_solve_many(&Csr { cols: 5, ..empty }, &[], 5),
        Err(NumericError::NotSquare)
    );

    let a = Csr {
        rows: 2,
        cols: 2,
        row_ptr: &[0, 1, 2],
        col_idx: &[0, 1],
        values: &[2.0, 3.0],
    };
    assert_eq!(solve_auto_many(&a, &[], 0), Ok((Vec::new(), "LU")));
    assert_eq!(cholesky_solve_many(&a, &[], 0), Ok(Vec::new()));
    assert_eq!(
        lu_solve_many(
            &Csr {
                values: &[2.0, 0.0],
                ..a
            },
            &[],
            0
        ),
        Ok(Vec::new())
    );
    for (invalid, rhs_cols) in [(&[][..], usize::MAX), (&[0.0; 3][..], 2)] {
        assert!(matches!(
            solve_auto_many(&a, invalid, rhs_cols),
            Err(NumericError::ShapeMismatch(_))
        ));
    }
    let broken = Csr {
        row_ptr: &[0, usize::MAX, 2],
        ..a
    };
    assert!(matches!(
        lu_solve_many(&broken, &[], 0),
        Err(NumericError::ShapeMismatch(_))
    ));
    let broken = Csr {
        col_idx: &[0, 2],
        ..a
    };
    assert!(matches!(
        lu_solve_many(&broken, &[], 0),
        Err(NumericError::ShapeMismatch(_))
    ));
    let broken = Csr {
        values: &[2.0, f64::INFINITY],
        ..a
    };
    assert!(matches!(
        lu_solve_many(&broken, &[], 0),
        Err(NumericError::Domain(_))
    ));
    assert_eq!(
        lu_solve_many(&a, &[1.0, f64::NAN], 1),
        Err(NumericError::Singular)
    );
}

#[test]
fn direct_solvers_preserve_widely_separated_diagonal_scales() {
    let subnormal = f64::from_bits(1);
    let a = Csr {
        rows: 2,
        cols: 2,
        row_ptr: &[0, 1, 2],
        col_idx: &[0, 1],
        values: &[subnormal, f64::MAX],
    };
    let rhs = [subnormal, -2.0 * subnormal, f64::MAX, -f64::MAX];
    let expected = [1.0, -2.0, 1.0, -1.0];
    check_values(&lu_solve_many(&a, &rhs, 2).unwrap(), &expected);
    check_values(&cholesky_solve_many(&a, &rhs, 2).unwrap(), &expected);
    let singular = Csr {
        values: &[0.0, f64::MAX],
        ..a
    };
    assert_eq!(
        lu_solve_many(&singular, &rhs, 2),
        Err(NumericError::Singular)
    );
}
