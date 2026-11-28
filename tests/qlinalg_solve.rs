//! Exact rational reference cases for the dense particular solver.

#![allow(clippy::unwrap_used)]

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Zero;
use rssn::kernels::qlinalg::solve;

type Q = BigRational;

fn q(
    n: i64,
    d: i64,
) -> Q {
    Q::new(BigInt::from(n), BigInt::from(d))
}

#[test]
fn validates_rhs_shape_before_reducing_or_returning_zero() {
    let matrix = [vec![Q::one()], vec![Q::one()]];
    assert_eq!(solve(&matrix, 1, &[q(2, 1)]), None);
    assert_eq!(solve(&matrix, 1, &[q(2, 1), q(2, 1), q(2, 1)]), None);
    assert_eq!(solve(&[vec![Q::one(), Q::zero()]], 1, &[q(2, 1)]), None);
    assert_eq!(solve(&[vec![]], 1, &[q(2, 1)]), None);
    assert_eq!(solve(&[vec![]], 1, &[Q::zero()]), None);
    assert_eq!(solve(&[], usize::MAX, &[]), None);
    assert_eq!(solve(&[], usize::MAX - 1, &[]), None);
}

#[test]
fn rectangular_system_keeps_earliest_pivots_and_zero_free_coordinates() {
    let a = [
        vec![Q::zero(), Q::zero(), Q::zero(), q(1, 3)],
        vec![Q::zero(), q(2, 3), q(1, 2), q(-4, 5)],
        vec![Q::zero(), q(4, 3), Q::one(), q(-8, 5)],
    ];
    let mut b = vec![q(1, 2), q(-22, 15), q(-44, 15)];
    let expected = vec![Q::zero(), q(-2, 5), Q::zero(), q(3, 2)];
    assert_eq!(solve(&a, 4, &b), Some(expected));
    b[2] += q(1, 11);
    assert!(solve(&a, 4, &b).is_none());
}

#[test]
fn homogeneous_and_zero_dimension_systems_match_their_definitions() {
    assert_eq!(solve(&[], 3, &[]), Some(vec![Q::zero(); 3]));
    assert_eq!(
        solve(&[vec![], vec![]], 0, &[Q::zero(), Q::zero()]),
        Some(vec![])
    );
    assert_eq!(solve(&[vec![]], 0, &[Q::one()]), None);
    let matrix = [vec![q(17, 23), q(-3, 11)], vec![q(34, 23), q(-6, 11)]];
    assert_eq!(
        solve(&matrix, 2, &[Q::zero(), Q::zero()]),
        Some(vec![Q::zero(); 2])
    );
}

// Independent Laplace expansion for the small Cramer's rule references.
fn determinant(matrix: &[Vec<Q>]) -> Q {
    if matrix.is_empty() {
        return Q::one();
    }
    matrix[0]
        .iter()
        .enumerate()
        .filter(|(_, v)| !v.is_zero())
        .map(|(column, value)| {
            let minor: Vec<_> = matrix[1..]
                .iter()
                .map(|row| {
                    row.iter()
                        .enumerate()
                        .filter(|(j, _)| *j != column)
                        .map(|(_, v)| v.clone())
                        .collect()
                })
                .collect();
            let term = value * determinant(&minor);
            if column % 2 == 0 {
                term
            } else {
                -term
            }
        })
        .sum()
}

#[test]
fn small_dense_systems_match_cramers_rule() {
    let mut seed = 515721_u64;
    let mut next = |modulus| {
        seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
        ((seed >> 32) % modulus) as i64
    };
    for size in [1, 2, 3, 4, 5] {
        for _ in 0..4 {
            let mut a: Vec<Vec<_>> = (0..size)
                .map(|i| {
                    (0..size)
                        .map(|j| {
                            if i == j {
                                q(32 + next(5), 2)
                            } else {
                                q(next(5) - 2, next(5) + 1)
                            }
                        })
                        .collect()
                })
                .collect();
            a.rotate_left(size - 1);
            let b: Vec<_> = (0..size).map(|_| q(next(15) - 7, next(9) + 1)).collect();
            let divisor = determinant(&a);
            assert!(!divisor.is_zero());
            let expected: Vec<_> = (0..size)
                .map(|column| {
                    let mut numerator = a.clone();
                    for (row, rhs) in numerator.iter_mut().zip(&b) {
                        row[column] = rhs.clone();
                    }
                    determinant(&numerator) / &divisor
                })
                .collect();
            assert_eq!(solve(&a, size, &b), Some(expected));
        }
    }
}

#[test]
fn exact_large_rational_coefficients_are_independent_of_worker_count() {
    let num = (BigInt::one() << 131_usize) + 49_u8;
    let den = (BigInt::one() << 109_usize) + 69_u8;
    let scale = Q::new(num, den);
    let cols = 6;
    let a: Vec<Vec<_>> = (0..72)
        .map(|i| {
            (0..cols)
                .map(|j| {
                    if i < cols {
                        if i == j {
                            &scale * q((i + 1) as i64, 7)
                        } else {
                            Q::zero()
                        }
                    } else {
                        q(((i * j + j + 2) % 9) as i64 - 4, ((i + j) % 7 + 1) as i64)
                    }
                })
                .collect()
        })
        .collect();
    let expected: Vec<_> = (0..cols)
        .map(|i| &scale * q(i as i64 - 2, i as i64 + 3))
        .collect();
    let b: Vec<_> = a
        .iter()
        .map(|row| row.iter().zip(&expected).map(|(x, y)| x * y).sum())
        .collect();
    for threads in [1, 3] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        assert_eq!(pool.install(|| solve(&a, cols, &b)), Some(expected.clone()));
    }
}
