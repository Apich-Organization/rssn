//! Linear algebra over a coefficient field.
#![allow(
    // The tower arithmetic indexes coefficient vectors and shifts degrees by
    // hand; every index comes from a length that the surrounding loop just
    // checked, and a tower of two or three levels keeps the vectors short.
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::pedantic,
    clippy::nursery
)]

use super::elem::Coeff;

/// Gauss–Jordan elimination. Returns one solution with the free variables
/// set to zero, or `None` when the system is inconsistent.
pub fn solve_field<F: Coeff>(
    rows: &mut [Vec<F>],
    rhs: &mut [F],
    columns: usize,
) -> Option<Vec<F>> {
    let m = rows.len();
    let mut pivots: Vec<usize> = Vec::new();
    let mut r = 0;
    for c in 0..columns {
        let Some(p) = (r..m).find(|&i| !rows[i][c].is_zero()) else {
            continue;
        };
        rows.swap(r, p);
        rhs.swap(r, p);
        let lead = rows[r][c].clone();
        let inv = lead.inv()?;
        for v in &mut rows[r] {
            *v = v.mul(&inv);
        }
        rhs[r] = rhs[r].mul(&inv);
        for i in 0..m {
            if i != r && !rows[i][c].is_zero() {
                let factor = rows[i][c].clone();
                let pivot = rows[r].clone();
                for (v, pv) in rows[i].iter_mut().zip(&pivot) {
                    *v = v.sub(&factor.mul(pv));
                }
                let pr = rhs[r].clone();
                rhs[i] = rhs[i].sub(&factor.mul(&pr));
            }
        }
        pivots.push(c);
        r += 1;
        if r == m {
            break;
        }
    }
    for i in r..m {
        if !rhs[i].is_zero() {
            return None;
        }
    }
    let mut out = vec![F::zero(); columns];
    for (i, &c) in pivots.iter().enumerate() {
        out[c] = rhs[i].clone();
    }
    Some(out)
}
