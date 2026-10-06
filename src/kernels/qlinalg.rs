//! # Exact linear algebra over the rationals
//!
//! Dense matrices of [`BigRational`] stored as rows: reduced row echelon
//! forms, ranks, null spaces, solutions of linear systems, inverses and
//! determinants. Used by the Lie-algebra structure theory and by the
//! homology bases of simplicial complexes.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Zero;

/// An exact rational number.
pub type Q = BigRational;

/// A matrix stored as a list of rows.
pub type QMat = Vec<Vec<Q>>;

/// The integer `n` as a rational.
#[must_use]
pub fn q(n: i64) -> Q {
    Q::from_integer(BigInt::from(n))
}

/// The zero vector of length `n`.
#[must_use]
pub fn zeros(n: usize) -> Vec<Q> {
    vec![Q::zero(); n]
}

/// The identity matrix of size `n`.
#[must_use]
pub fn identity(n: usize) -> QMat {
    (0..n)
        .map(|i| {
            let mut row = zeros(n);
            row[i] = Q::one();
            row
        })
        .collect()
}

/// Brings `m` (rows of length `cols`) to reduced row echelon form in place,
/// drops the zero rows and returns the pivot columns.
pub fn rref(
    m: &mut QMat,
    cols: usize,
) -> Vec<usize> {
    let mut pivots = Vec::new();
    let mut r = 0;
    for c in 0..cols {
        if r == m.len() {
            break;
        }
        let Some(p) = (r..m.len()).find(|&i| !m[i][c].is_zero()) else {
            continue;
        };
        m.swap(r, p);
        let inv = m[r][c].recip();
        for x in &mut m[r] {
            *x = &*x * &inv;
        }
        let pivot_row = m[r].clone();
        for (i, row) in m.iter_mut().enumerate() {
            if i == r || row[c].is_zero() {
                continue;
            }
            let f = row[c].clone();
            for (x, y) in row.iter_mut().zip(&pivot_row) {
                *x -= &f * y;
            }
        }
        pivots.push(c);
        r += 1;
    }
    m.truncate(r);
    pivots
}

/// The rank of a matrix whose rows have length `cols`.
#[must_use]
pub fn rank(
    m: &[Vec<Q>],
    cols: usize,
) -> usize {
    let mut a = m.to_vec();
    rref(&mut a, cols).len()
}

/// A basis (reduced echelon rows) of the span of `vectors` of length `dim`.
#[must_use]
pub fn span_basis(
    vectors: &[Vec<Q>],
    dim: usize,
) -> QMat {
    let mut a = vectors.to_vec();
    rref(&mut a, dim);
    a
}

/// A basis of `{x : m x = 0}` for `m` with `cols` columns.
#[must_use]
pub fn nullspace(
    m: &[Vec<Q>],
    cols: usize,
) -> QMat {
    let mut a = m.to_vec();
    let pivots = rref(&mut a, cols);
    let mut out = Vec::new();
    for free in (0..cols).filter(|c| !pivots.contains(c)) {
        let mut x = zeros(cols);
        x[free] = Q::one();
        for (row, &p) in a.iter().zip(&pivots) {
            x[p] = -row[free].clone();
        }
        out.push(x);
    }
    out
}

/// A particular solution of `m x = b` (free variables set to zero), if any.
#[must_use]
pub fn solve(
    m: &[Vec<Q>],
    cols: usize,
    b: &[Q],
) -> Option<Vec<Q>> {
    let mut aug: QMat = m.iter().zip(b).map(|(row, bi)| row.iter().cloned().chain(std::iter::once(bi.clone())).collect()).collect();
    let pivots = rref(&mut aug, cols + 1);
    if pivots.last() == Some(&cols) {
        return None;
    }
    let mut x = zeros(cols);
    for (row, &p) in aug.iter().zip(&pivots) {
        x[p] = row[cols].clone();
    }
    Some(x)
}

/// The product of two matrices (`a` is `r x m`, `b` is `m x c`).
#[must_use]
pub fn matmul(
    a: &[Vec<Q>],
    b: &[Vec<Q>],
) -> QMat {
    let cols = b.first().map_or(0, Vec::len);
    a.iter()
        .map(|row| {
            (0..cols)
                .map(|j| {
                    let mut s = Q::zero();
                    for (x, brow) in row.iter().zip(b) {
                        if !x.is_zero() {
                            s += x * &brow[j];
                        }
                    }
                    s
                })
                .collect()
        })
        .collect()
}

/// The product of a matrix and a vector.
#[must_use]
pub fn matvec(
    a: &[Vec<Q>],
    v: &[Q],
) -> Vec<Q> {
    a.iter()
        .map(|row| {
            let mut s = Q::zero();
            for (x, y) in row.iter().zip(v) {
                if !x.is_zero() && !y.is_zero() {
                    s += x * y;
                }
            }
            s
        })
        .collect()
}

/// The transpose of an `r x c` matrix.
#[must_use]
pub fn transpose(
    a: &[Vec<Q>],
    cols: usize,
) -> QMat {
    (0..cols).map(|j| a.iter().map(|row| row[j].clone()).collect()).collect()
}

/// The inverse of a square matrix.
#[must_use]
pub fn inverse(m: &[Vec<Q>]) -> Option<QMat> {
    let n = m.len();
    let mut aug: QMat = m.iter().enumerate().map(|(i, row)| row.iter().cloned().chain((0..n).map(|j| q(i64::from(i == j)))).collect()).collect();
    let pivots = rref(&mut aug, 2 * n);
    if pivots.len() != n || pivots.last().is_some_and(|&p| p >= n) {
        return None;
    }
    Some(aug.into_iter().map(|row| row[n..].to_vec()).collect())
}

/// The determinant of a square matrix.
#[must_use]
pub fn det(m: &[Vec<Q>]) -> Q {
    let n = m.len();
    let mut a = m.to_vec();
    let mut d = Q::one();
    for c in 0..n {
        let Some(p) = (c..n).find(|&i| !a[i][c].is_zero()) else {
            return Q::zero();
        };
        if p != c {
            a.swap(p, c);
            d = -d;
        }
        d *= a[c][c].clone();
        let pivot_row = a[c].clone();
        for row in a.iter_mut().skip(c + 1) {
            if row[c].is_zero() {
                continue;
            }
            let f = &row[c] / &pivot_row[c];
            for (x, y) in row.iter_mut().zip(&pivot_row).skip(c) {
                *x -= &f * y;
            }
        }
    }
    d
}

/// Reduces `v` modulo the span of the echelon rows `basis` (with their
/// `pivots`, as returned by [`rref`]): the remainder vanishes exactly on
/// the span.
#[must_use]
pub fn reduce_mod(
    v: &[Q],
    basis: &[Vec<Q>],
    pivots: &[usize],
) -> Vec<Q> {
    let mut v = v.to_vec();
    for (row, &p) in basis.iter().zip(pivots) {
        if v[p].is_zero() {
            continue;
        }
        let f = v[p].clone();
        for (x, y) in v.iter_mut().zip(row) {
            *x -= &f * y;
        }
    }
    v
}

/// Whether a vector is zero.
#[must_use]
pub fn is_zero_vec(v: &[Q]) -> bool {
    v.iter().all(Zero::is_zero)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn m(rows: &[&[i64]]) -> QMat {
        rows.iter().map(|r| r.iter().map(|&x| q(x)).collect()).collect()
    }

    #[test]
    fn linear_algebra() {
        let a = m(&[&[1, 2, 3], &[2, 4, 6], &[1, 0, 1]]);
        assert_eq!(rank(&a, 3), 2);
        assert_eq!(nullspace(&a, 3).len(), 1);
        assert_eq!(det(&m(&[&[2, 1], &[1, 1]])), q(1));
        let inv = inverse(&m(&[&[2, 1], &[1, 1]])).expect("invertible");
        assert_eq!(inv, m(&[&[1, -1], &[-1, 2]]));
        assert!(inverse(&a).is_none());
        let x = solve(&m(&[&[1, 1], &[1, -1]]), 2, &[q(3), q(1)]).expect("solvable");
        assert_eq!(x, vec![q(2), q(1)]);
        assert!(solve(&m(&[&[1, 1], &[1, 1]]), 2, &[q(1), q(2)]).is_none());
        assert_eq!(matmul(&m(&[&[1, 2]]), &m(&[&[3], &[4]])), m(&[&[11]]));
        let mut b = m(&[&[0, 2], &[1, 1]]);
        let pivots = rref(&mut b, 2);
        assert_eq!(pivots, vec![0, 1]);
        assert!(is_zero_vec(&reduce_mod(&[q(5), q(7)], &b, &pivots)));
    }
}
