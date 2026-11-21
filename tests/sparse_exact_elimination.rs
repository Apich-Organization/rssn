//! Regression coverage for exact sparse row updates and augmented solves.

#![allow(clippy::unwrap_used)]

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::Zero;
use rssn::kernels::sparse_exact::Domain;
use rssn::kernels::sparse_exact::Rationals;
use rssn::kernels::sparse_exact::Solution;
use rssn::kernels::sparse_exact::SparseMat;
use rssn::kernels::sparse_exact::echelon;
use rssn::kernels::sparse_exact::solve;

fn q(n: i64) -> BigRational {
    BigRational::from_integer(BigInt::from(n))
}

fn frac(
    n: i64,
    d: i64,
) -> BigRational {
    BigRational::new(BigInt::from(n), BigInt::from(d))
}

#[test]
fn records_fill_and_cancellation_in_sparse_pivot_order() {
    let original = [
        [0, 2, 1, -1, 2, 0],
        [1, 1, 2, 1, 2, 0],
        [0, -1, 0, -1, 0, -1],
        [-1, 0, 0, 1, 0, 0],
        [2, 1, 0, 2, 2, 1],
    ]
    .map(|row| row.map(q).to_vec());
    let matrix = SparseMat::from_dense(&mut Rationals, &original);
    let result = echelon(&mut Rationals, &matrix, 4, true);

    // The first update of row 1 adds column 5 and cancels column 3. The
    // second reintroduces column 3, which must count exactly once.
    assert_eq!(result.pivots, [(2, 1), (0, 2), (1, 0), (3, 3)]);
    assert_eq!(
        result.ops,
        [
            (0, 0, q(-2)),
            (1, 0, q(-1)),
            (4, 0, q(-1)),
            (1, 1, q(2)),
            (3, 2, q(-1)),
            (4, 2, q(2)),
            (4, 3, frac(-11, 7)),
        ]
    );
    assert_eq!(result.rest, [(4, vec![(4, frac(20, 7)), (5, frac(-9, 7))])]);

    let mut dense = original;
    for (target, pivot, factor) in &result.ops {
        let source = dense[result.pivots[*pivot].0].clone();
        for (value, entry) in dense[*target].iter_mut().zip(&source) {
            *value -= factor * entry;
        }
    }
    for (&(i, _), row) in result.pivots.iter().zip(&result.rows) {
        let expected: Vec<_> = dense[i]
            .iter()
            .enumerate()
            .filter(|(_, value)| !value.is_zero())
            .map(|(j, value)| (j, value.clone()))
            .collect();
        assert_eq!(*row, expected);
    }
    assert_eq!(dense[4], [q(0), q(0), q(0), q(0), frac(20, 7), frac(-9, 7)]);

    let unrecorded = echelon(&mut Rationals, &matrix, 4, false);
    assert!(unrecorded.ops.is_empty());
    assert_eq!(unrecorded.pivots, result.pivots);
    assert_eq!(unrecorded.rows, result.rows);
    assert_eq!(unrecorded.rest, result.rest);
}

#[test]
fn augmented_solve_matches_an_exact_rational_family() {
    // Two independently specified constraints: x + (2/3) z = 5/7 and
    // y - (1/4) z = -2/5. The first and last rows are linear combinations.
    let coefficients = [
        vec![q(2), q(-3), frac(25, 12)],
        vec![q(1), q(0), frac(2, 3)],
        vec![q(0), q(1), frac(-1, 4)],
        vec![q(-2), q(3), frac(-25, 12)],
    ];
    let mut rhs = [frac(92, 35), frac(5, 7), frac(-2, 5), frac(-92, 35)];
    let matrix = SparseMat::from_dense(&mut Rationals, &coefficients);
    let found = solve(&mut Rationals, &matrix, &rhs).unwrap();
    let Solution::Family {
        particular,
        nullspace,
    } = found
    else {
        panic!("consistent system returned {found:?}");
    };
    assert_eq!(particular[0], frac(5, 7) - frac(2, 3) * &particular[2]);
    assert_eq!(particular[1], frac(-2, 5) + frac(1, 4) * &particular[2]);
    assert_eq!(nullspace.len(), 1);
    let direction = &nullspace[0];
    assert!(!direction[2].is_zero());
    assert_eq!(direction[0], frac(-2, 3) * &direction[2]);
    assert_eq!(direction[1], frac(1, 4) * &direction[2]);

    rhs[3] += frac(1, 97);
    assert_eq!(
        solve(&mut Rationals, &matrix, &rhs),
        Some(Solution::Inconsistent)
    );
    assert_eq!(solve(&mut Rationals, &matrix, &rhs[..3]), None);
}

#[test]
fn checks_rhs_and_augmented_dimensions_before_elimination() {
    // The first shape would overflow the column count and the second
    // would overflow the column workspace, despite having no rows.
    for cols in [usize::MAX, usize::MAX - 1] {
        let unaugmentable = SparseMat::zeros(0, cols);
        assert_eq!(solve(&mut Rationals, &unaugmentable, &[]), None);
    }

    let mut invalid_row_count = SparseMat::zeros(2, 1);
    invalid_row_count.data.truncate(1);
    assert_eq!(
        solve(&mut Rationals, &invalid_row_count, &[q(0), q(0)]),
        None
    );

    let zero_columns = SparseMat::zeros(2, 0);
    assert_eq!(
        solve(&mut Rationals, &zero_columns, &[q(0), q(0)]),
        Some(Solution::Family {
            particular: vec![],
            nullspace: vec![]
        })
    );
    assert_eq!(
        solve(&mut Rationals, &zero_columns, &[q(0), q(1)]),
        Some(Solution::Inconsistent)
    );
}

#[derive(Default)]
struct LoggingIntegerDivision {
    divisions: Vec<(i64, i64)>,
}

impl Domain for LoggingIntegerDivision {
    type E = i64;

    fn zero(&mut self) -> i64 {
        0
    }

    fn one(&mut self) -> i64 {
        1
    }

    fn is_zero(
        &mut self,
        a: &i64,
    ) -> bool {
        *a == 0
    }

    fn add(
        &mut self,
        a: &i64,
        b: &i64,
    ) -> i64 {
        a + b
    }

    fn sub(
        &mut self,
        a: &i64,
        b: &i64,
    ) -> i64 {
        a - b
    }

    fn mul(
        &mut self,
        a: &i64,
        b: &i64,
    ) -> i64 {
        a * b
    }

    fn div(
        &mut self,
        a: &i64,
        b: &i64,
    ) -> Option<i64> {
        self.divisions.push((*a, *b));
        if a.checked_rem(*b) == Some(0) {
            a.checked_div(*b)
        } else {
            None
        }
    }
}

#[test]
fn pivot_order_retains_rows_after_failed_division() {
    for padding in [0, 128] {
        let mut dom = LoggingIntegerDivision::default();
        let mut rows = vec![
            vec![2, 0, 0, 0, 0],
            vec![3, 1, 0, 0, 0],
            vec![4, 0, 1, 0, 0],
            vec![0, 1, 0, 0, 0],
            vec![0, 0, 0, 0, 0],
            vec![0, 0, 0, 0, 9],
        ];
        rows.resize(6 + padding, vec![0; 5]);
        let matrix = SparseMat::from_dense(&mut dom, &rows);
        let result = echelon(&mut dom, &matrix, 4, true);

        // 3/2 fails, so the row stays eligible and retains its early entry.
        // Updating row 2 brings its width down to one, ahead of row 3 in the tie.
        assert_eq!(dom.divisions, [(3, 2), (4, 2), (1, 1)]);
        assert_eq!(result.pivots, [(0, 0), (2, 2), (3, 1), (1, 0)]);
        assert_eq!(result.ops, [(2, 0, 2), (1, 2, 1)]);
        assert_eq!(&result.rest[..2], [(4, vec![]), (5, vec![(4, 9)])]);
        assert_eq!(result.rest.len(), 2 + padding);
        assert!(result.rest[2..].iter().all(|(_, row)| row.is_empty()));
    }
}
