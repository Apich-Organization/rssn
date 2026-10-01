//! Partial fraction decomposition over the rationals.
//!
//! `n / d = q + sum_i sum_{j <= m_i} a_ij / f_i^j`, where `d = c prod f_i^m_i`
//! is the factorisation of the denominator into irreducible factors over `Q`
//! and every numerator `a_ij` has degree below `deg f_i`. The numerators are
//! found by solving one linear system over `Q` (the classical
//! undetermined-coefficients method), which is exact and handles repeated
//! and irreducible higher-degree factors uniformly.

use num_rational::BigRational;
use num_traits::One;
use num_traits::Zero;

use super::univariate;
use super::univariate::QPoly;

/// One term `numerator / factor^power` of a decomposition.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Piece {
    /// An irreducible factor of the denominator, primitive over `Z`.
    pub factor: QPoly,
    /// Its power in this term.
    pub power: u32,
    /// A numerator of degree below `deg factor` (never zero).
    pub numerator: QPoly,
}

/// A complete decomposition of a rational function.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Apart {
    /// The polynomial part.
    pub quotient: QPoly,
    /// The proper fractions.
    pub pieces: Vec<Piece>,
}

/// Decomposes `numer / denom`. Returns `None` for a zero denominator.
#[must_use]
pub fn apart(
    numer: &[BigRational],
    denom: &[BigRational],
) -> Option<Apart> {
    let trim = |p: &[BigRational]| {
        let mut p = p.to_vec();
        while p.last().is_some_and(Zero::is_zero) {
            p.pop();
        }
        p
    };
    let (numer, denom) = (trim(numer), trim(denom));
    if denom.is_empty() {
        return None;
    }
    let (quotient, remainder) = univariate::divrem(&numer, &denom)?;
    if remainder.is_empty() || denom.len() < 2 {
        return Some(Apart { quotient, pieces: Vec::new() });
    }
    let (_, factors) = univariate::factor(&denom);
    // Unknown numerators: for factor i and power j, a polynomial of degree
    // below deg(f_i). Column = x^k * denom / f_i^j.
    let n = denom.len() - 1;
    let mut columns: Vec<QPoly> = Vec::with_capacity(n);
    let mut labels: Vec<(usize, u32)> = Vec::with_capacity(n);
    let rational_factors: Vec<QPoly> =
        factors.iter().map(|(f, _)| f.iter().cloned().map(BigRational::from_integer).collect()).collect();
    for (i, ((factor, multiplicity), fq)) in factors.iter().zip(&rational_factors).enumerate() {
        let mut power: QPoly = vec![BigRational::one()];
        for j in 1..=*multiplicity {
            power = univariate::mul(&power, fq);
            let (cofactor, _) = univariate::divrem(&denom, &power)?;
            for k in 0..factor.len() - 1 {
                let mut shifted = vec![BigRational::zero(); k];
                shifted.extend(cofactor.iter().cloned());
                columns.push(shifted);
                labels.push((i, j));
            }
        }
    }
    if columns.len() != n {
        return None;
    }
    let solution = solve_columns(&columns, &remainder, n)?;
    let mut pieces = Vec::new();
    for (i, ((_, multiplicity), fq)) in factors.iter().zip(&rational_factors).enumerate() {
        for j in 1..=*multiplicity {
            let mut numerator: QPoly =
                labels.iter().zip(&solution).filter(|((fi, fj), _)| *fi == i && *fj == j).map(|(_, v)| v.clone()).collect();
            while numerator.last().is_some_and(Zero::is_zero) {
                numerator.pop();
            }
            if numerator.is_empty() {
                continue;
            }
            pieces.push(Piece { factor: fq.clone(), power: j, numerator });
        }
    }
    Some(Apart { quotient, pieces })
}

/// Solves `sum_k a_k * column_k = target` over the rationals, where every
/// polynomial has degree below `n` and there are `n` columns.
fn solve_columns(
    columns: &[QPoly],
    target: &[BigRational],
    n: usize,
) -> Option<Vec<BigRational>> {
    // Augmented matrix: row = coefficient of x^row.
    let mut m: Vec<Vec<BigRational>> = (0..n)
        .map(|row| {
            let mut line: Vec<BigRational> =
                columns.iter().map(|c| c.get(row).cloned().unwrap_or_else(BigRational::zero)).collect();
            line.push(target.get(row).cloned().unwrap_or_else(BigRational::zero));
            line
        })
        .collect();
    for col in 0..n {
        let pivot = (col..n).find(|&r| !m[r][col].is_zero())?;
        m.swap(col, pivot);
        let lead = m[col][col].clone();
        for value in &mut m[col] {
            *value = &*value / &lead;
        }
        for row in 0..n {
            if row != col && !m[row][col].is_zero() {
                let factor = m[row][col].clone();
                let pivot_row = m[col].clone();
                for (value, p) in m[row].iter_mut().zip(&pivot_row) {
                    *value = &*value - &factor * p;
                }
            }
        }
    }
    Some(m.into_iter().map(|row| row.last().cloned().unwrap_or_else(BigRational::zero)).collect())
}

#[cfg(test)]
mod tests {
    use num_bigint::BigInt;

    use super::*;

    fn q(values: &[i64]) -> QPoly {
        values.iter().map(|&v| BigRational::from_integer(BigInt::from(v))).collect()
    }

    /// Recombines a decomposition over the original denominator.
    fn recombine(
        parts: &Apart,
        denom: &[BigRational],
    ) -> QPoly {
        let mut total = univariate::mul(&parts.quotient, denom);
        for piece in &parts.pieces {
            let mut power: QPoly = vec![BigRational::one()];
            for _ in 0..piece.power {
                power = univariate::mul(&power, &piece.factor);
            }
            let (cofactor, rest) = univariate::divrem(denom, &power).unwrap_or_default();
            assert!(rest.is_empty());
            total = univariate::add(&total, &univariate::mul(&cofactor, &piece.numerator));
        }
        total
    }

    #[test]
    fn decompositions_recombine() {
        let cases: [(&[i64], &[i64]); 5] = [
            (&[1], &[-1, 0, 1]),                 // 1/(x^2 - 1)
            (&[0, 0, 0, 1], &[-1, 0, 1]),        // x^3/(x^2 - 1)
            (&[1, 2], &[0, 0, 1, 1]),            // (2x + 1)/(x^2 (x + 1))
            (&[3, 0, 1], &[1, 0, 2, 0, 1]),      // (x^2 + 3)/(x^2 + 1)^2
            (&[1], &[-1, 0, 0, 1]),              // 1/(x^3 - 1)
        ];
        for (n, d) in cases {
            let (n, d) = (q(n), q(d));
            let parts = apart(&n, &d).unwrap_or_else(|| panic!("no decomposition"));
            assert_eq!(recombine(&parts, &d), n, "{parts:?}");
            for piece in &parts.pieces {
                assert!(piece.numerator.len() < piece.factor.len());
            }
        }
    }

    #[test]
    fn repeated_irreducible_quadratic() {
        let parts = apart(&q(&[3, 0, 1]), &q(&[1, 0, 2, 0, 1])).unwrap_or_else(|| panic!("none"));
        // (x^2 + 3)/(x^2 + 1)^2 = 1/(x^2 + 1) + 2/(x^2 + 1)^2
        assert_eq!(parts.pieces.len(), 2);
        assert_eq!(parts.pieces[0].numerator, q(&[1]));
        assert_eq!(parts.pieces[1].numerator, q(&[2]));
    }
}
