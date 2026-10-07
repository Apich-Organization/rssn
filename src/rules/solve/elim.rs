//! Determinants and resultants over polynomial rings, for eliminating an
//! unknown from a pair of polynomial equations.

use super::cap_default;
use crate::graph::Number;
use crate::rules::poly::repr::Poly;

/// Determinant of a square matrix of polynomials by the fraction-free
/// Bareiss algorithm (every division is exact).
pub(super) fn determinant(mut m: Vec<Vec<Poly>>) -> Option<Poly> {
    let n = m.len();
    if n == 0 {
        return Some(Poly::constant(Number::from(1)));
    }
    let cap = cap_default();
    let mut negate = false;
    let mut previous = Poly::constant(Number::from(1));
    for k in 0..n.saturating_sub(1) {
        if m.get(k)?.get(k)?.is_zero() {
            let swap = (k + 1..n).find(|&i| m.get(i).and_then(|r| r.get(k)).is_some_and(|p| !p.is_zero()));
            match swap {
                | Some(i) => {
                    m.swap(k, i);
                    negate = !negate;
                },
                | None => return Some(Poly::zero()),
            }
        }
        let pivot = m.get(k)?.get(k)?.clone();
        for i in k + 1..n {
            let lead = m.get(i)?.get(k)?.clone();
            for j in k + 1..n {
                let kj = m.get(k)?.get(j)?.clone();
                let ij = m.get(i)?.get(j)?.clone();
                let numerator = ij.mul(&pivot, cap)?.sub(&lead.mul(&kj, cap)?);
                let value = if previous.as_constant().is_some_and(|c| c == Number::from(1)) {
                    numerator
                } else {
                    numerator.div_exact(&previous, cap)?
                };
                *m.get_mut(i)?.get_mut(j)? = value;
            }
        }
        previous = pivot;
    }
    let last = m.get(n - 1)?.get(n - 1)?.clone();
    Some(if negate { last.neg() } else { last })
}

/// The resultant of `f` and `g` with respect to the generator `v`.
pub(super) fn resultant(
    f: &Poly,
    g: &Poly,
    v: u32,
) -> Option<Poly> {
    let (df, dg) = (usize::try_from(f.degree_in(v)).ok()?, usize::try_from(g.degree_in(v)).ok()?);
    if df == 0 || dg == 0 {
        return None;
    }
    // Coefficients, highest degree first.
    let mut fc = f.coefficients_in(v);
    fc.reverse();
    let mut gc = g.coefficients_in(v);
    gc.reverse();
    let size = df + dg;
    let mut matrix: Vec<Vec<Poly>> = Vec::with_capacity(size);
    for i in 0..dg {
        let mut row = vec![Poly::zero(); size];
        for (j, c) in fc.iter().enumerate() {
            *row.get_mut(i + j)? = c.clone();
        }
        matrix.push(row);
    }
    for i in 0..df {
        let mut row = vec![Poly::zero(); size];
        for (j, c) in gc.iter().enumerate() {
            *row.get_mut(i + j)? = c.clone();
        }
        matrix.push(row);
    }
    determinant(matrix)
}
