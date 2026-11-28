//! Sparse multivariate polynomials over `Q` in the flat variable set
//! `(x, theta_1, ..., theta_n)`, indexed `0, 1, ..., n`.
#![allow(
    // The tower arithmetic indexes coefficient vectors and shifts degrees by
    // hand; every index comes from a length that the surrounding loop just
    // checked, and a tower of two or three levels keeps the vectors short.
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::pedantic,
    clippy::nursery
)]

use std::collections::BTreeMap;

use num_bigint::BigInt;
use num_traits::Zero;

use super::Q;

/// A monomial: `(variable, exponent)` pairs, sorted by variable, positive
/// exponents only.
pub type Mono = Vec<(u16, u32)>;

#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct Flat {
    terms: BTreeMap<Mono, Q>,
}

impl Flat {
    pub fn zero() -> Self {
        Self::default()
    }

    pub fn one() -> Self {
        Self::constant(Q::from_integer(BigInt::from(1)))
    }

    pub fn constant(value: Q) -> Self {
        let mut terms = BTreeMap::new();
        if !value.is_zero() {
            terms.insert(Vec::new(), value);
        }
        Self { terms }
    }

    pub fn monomial(
        var: u16,
        exp: u32,
    ) -> Self {
        let mut terms = BTreeMap::new();
        if exp > 0 {
            terms.insert(vec![(var, exp)], Q::from_integer(BigInt::from(1)));
        } else {
            terms.insert(Vec::new(), Q::from_integer(BigInt::from(1)));
        }
        Self { terms }
    }

    pub fn len(&self) -> usize {
        self.terms.len()
    }

    pub fn terms(&self) -> impl Iterator<Item = (&Mono, &Q)> {
        self.terms.iter()
    }

    pub fn add(
        &self,
        other: &Self,
    ) -> Self {
        let mut out = self.clone();
        for (mono, c) in &other.terms {
            let slot = out.terms.entry(mono.clone()).or_insert_with(Q::zero);
            *slot += c;
            if slot.is_zero() {
                out.terms.remove(mono);
            }
        }
        out
    }

    pub fn mul(
        &self,
        other: &Self,
    ) -> Self {
        if self.terms.is_empty() || other.terms.is_empty() {
            return Self::zero();
        }
        let mut out: BTreeMap<Mono, Q> = BTreeMap::new();
        for (ma, ca) in &self.terms {
            for (mb, cb) in &other.terms {
                let mono = mono_mul(ma, mb);
                let slot = out.entry(mono).or_insert_with(Q::zero);
                *slot += ca * cb;
            }
        }
        out.retain(|_, c| !c.is_zero());
        Self { terms: out }
    }

    /// Shifts every monomial down by `mins` (the common factor).
    fn shift(
        &self,
        mins: &BTreeMap<u16, u32>,
    ) -> Self {
        if mins.values().all(|&e| e == 0) {
            return self.clone();
        }
        let terms = self
            .terms
            .iter()
            .map(|(mono, c)| {
                let shifted: Mono = mono
                    .iter()
                    .map(|&(v, e)| (v, e - mins.get(&v).copied().unwrap_or(0)))
                    .filter(|&(_, e)| e > 0)
                    .collect();
                (shifted, c.clone())
            })
            .collect();
        Self { terms }
    }
}

/// Divides the common monomial factor out of a fraction `n/d`.
pub fn reduce_pair(
    n: &Flat,
    d: &Flat,
) -> (Flat, Flat) {
    let mut mins: BTreeMap<u16, u32> = BTreeMap::new();
    let first = n
        .terms()
        .map(|(m, _)| m)
        .chain(d.terms().map(|(m, _)| m))
        .next();
    let Some(first) = first else {
        return (Flat::zero(), Flat::one());
    };
    for &(v, e) in first {
        mins.insert(v, e);
    }
    for mono in n.terms().map(|(m, _)| m).chain(d.terms().map(|(m, _)| m)) {
        for (v, slot) in mins.iter_mut() {
            let e = mono.iter().find(|&&(w, _)| w == *v).map_or(0, |&(_, e)| e);
            *slot = (*slot).min(e);
        }
    }
    (n.shift(&mins), d.shift(&mins))
}

fn mono_mul(
    a: &Mono,
    b: &Mono,
) -> Mono {
    let mut out: Mono = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        match a[i].0.cmp(&b[j].0) {
            | std::cmp::Ordering::Less => {
                out.push(a[i]);
                i += 1;
            },
            | std::cmp::Ordering::Greater => {
                out.push(b[j]);
                j += 1;
            },
            | std::cmp::Ordering::Equal => {
                out.push((a[i].0, a[i].1 + b[j].1));
                i += 1;
                j += 1;
            },
        }
    }
    out.extend_from_slice(&a[i..]);
    out.extend_from_slice(&b[j..]);
    out
}

/// Gaussian elimination over `Q`. Returns one solution with the free
/// variables set to zero, or `None` when the system is inconsistent.
pub fn solve_q(
    rows: &mut [Vec<Q>],
    rhs: &mut [Q],
    columns: usize,
) -> Option<Vec<Q>> {
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
        for v in &mut rows[r] {
            *v /= &lead;
        }
        rhs[r] /= &lead;
        for i in 0..m {
            if i != r && !rows[i][c].is_zero() {
                let factor = rows[i][c].clone();
                let pivot = rows[r].clone();
                for (v, pv) in rows[i].iter_mut().zip(&pivot) {
                    *v -= &factor * pv;
                }
                let pr = rhs[r].clone();
                rhs[i] -= &factor * pr;
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
    let mut out = vec![Q::zero(); columns];
    for (i, &c) in pivots.iter().enumerate() {
        out[c] = rhs[i].clone();
    }
    Some(out)
}
