//! Exact Lie-algebra coordinates against independent matrix formulas.

use num_bigint::BigInt;
use num_traits::Zero;
use rssn::kernels::lie_structure::Lie;
use rssn::kernels::qlinalg::Q;
use rssn::kernels::qlinalg::QMat;
use rssn::kernels::qlinalg::q;

fn rat(
    num: i64,
    den: i64,
) -> Q {
    Q::new(BigInt::from(num), BigInt::from(den))
}

fn matrix2(entries: [Q; 4]) -> QMat {
    let [a, b, c, d] = entries;
    vec![vec![a, b], vec![c, d]]
}

fn matrix_unit(
    order: usize,
    row: usize,
    column: usize,
) -> QMat {
    let mut matrix = vec![vec![q(0); order]; order];
    matrix[row][column] = q(1);
    matrix
}

fn commutator(
    a: &QMat,
    b: &QMat,
) -> QMat {
    let order = a.len();
    (0..order)
        .map(|i| {
            (0..order)
                .map(|j| {
                    (0..order)
                        .map(|k| &a[i][k] * &b[k][j] - &b[i][k] * &a[k][j])
                        .sum()
                })
                .collect()
        })
        .collect()
}

fn represent(
    basis: &[QMat],
    coordinates: &[Q],
) -> QMat {
    let order = basis[0].len();
    (0..order)
        .map(|i| {
            (0..order)
                .map(|j| {
                    basis
                        .iter()
                        .zip(coordinates)
                        .map(|(matrix, coefficient)| coefficient * &matrix[i][j])
                        .sum()
                })
                .collect()
        })
        .collect()
}

fn trace(matrix: &QMat) -> Q {
    matrix.iter().enumerate().map(|(i, row)| &row[i]).sum()
}

#[test]
fn rational_gl2_in_an_unordered_noncanonical_basis() {
    let basis = [
        matrix2([q(0), rat(3, 2), q(0), rat(-1, 3)]),
        matrix2([rat(-1, 2), q(1), q(0), q(0)]),
        matrix2([q(0), q(0), q(0), rat(2, 3)]),
        matrix2([q(4), q(-1), q(2), rat(1, 5)]),
    ];
    let lie = Lie::from_matrices(&basis).unwrap();
    assert_eq!(lie.n, basis.len());
    for (i, a) in basis.iter().enumerate() {
        assert!(lie.c[i][i].iter().all(Q::is_zero));
        for (j, b) in basis.iter().enumerate() {
            assert_eq!(represent(&basis, &lie.c[i][j]), commutator(a, b));
            assert!(lie.c[i][j].iter().zip(&lie.c[j][i]).all(|(x, y)| x == &-y));
        }
    }

    let center = lie.center();
    assert_eq!(center.len(), 1);
    let center_matrix = represent(&basis, &center[0]);
    assert!(center_matrix[0][1].is_zero() && center_matrix[1][0].is_zero());
    assert!(!center_matrix[0][0].is_zero());
    assert_eq!(center_matrix[0][0], center_matrix[1][1]);

    // On gl_d, the exact Killing form is 2d tr(AB) - 2 tr(A) tr(B).
    let killing = lie.killing();
    for (i, a) in basis.iter().enumerate() {
        for (j, b) in basis.iter().enumerate() {
            let trace_ab: Q = (0..2)
                .flat_map(|i| (0..2).map(move |j| &a[i][j] * &b[j][i]))
                .sum();
            assert_eq!(killing[i][j], q(4) * trace_ab - q(2) * trace(a) * trace(b));
        }
    }
}

#[test]
fn detects_closure_beyond_pivot_entries_in_a_sparse_matrix_embedding() {
    let a = matrix_unit(6, 1, 3);
    let b = matrix_unit(6, 3, 5);
    let z = matrix_unit(6, 1, 5);
    assert!(Lie::from_matrices(&[a.clone(), b.clone()]).is_none());
    let lie = Lie::from_matrices(&[z, b, a]).unwrap();
    assert_eq!(lie.c[2][1], [q(1), q(0), q(0)]);
    assert_eq!(lie.c[1][2], [q(-1), q(0), q(0)]);
    assert!(lie.c[0].iter().flatten().all(Q::is_zero));
    assert_eq!(
        lie.lower_central_series()
            .iter()
            .map(Vec::len)
            .collect::<Vec<_>>(),
        [3, 1, 0]
    );
}

#[test]
fn rejects_empty_malformed_and_dependent_matrix_bases() {
    let ragged = vec![vec![q(1), q(1)], vec![q(1)]];
    let rectangular = vec![vec![q(1)], vec![q(2)]];
    let e = matrix_unit(2, 0, 1);
    for basis in [
        vec![],
        vec![vec![]],
        vec![ragged],
        vec![rectangular],
        vec![e.clone(), matrix_unit(3, 1, 1)],
        vec![e.clone(), e.clone()],
        vec![vec![vec![q(0)]]],
    ] {
        assert!(Lie::from_matrices(&basis).is_none());
    }
    let abelian = Lie::from_matrices(&[e]).unwrap();
    assert_eq!(abelian.n, 1);
    assert_eq!(abelian.c[0][0], [q(0)]);
}
