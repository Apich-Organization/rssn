//! Gröbner bases over the rationals (Buchberger's algorithm).
//!
//! Polynomials are stored densely sorted by the chosen monomial order, with
//! exponent vectors over a fixed list of variables. The implementation uses
//! the normal selection strategy (smallest least common multiple first) and
//! Buchberger's two criteria, and returns the *reduced* basis, which is
//! unique for a given ideal and order — so two descriptions of the same
//! ideal reduce to literally the same list.

use std::cmp::Ordering;

use num_rational::BigRational;
use num_traits::One;
use num_traits::Zero;

/// A monomial order.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Order {
    /// Lexicographic: eliminates variables in the given order.
    Lex,
    /// Total degree, ties broken lexicographically.
    GradedLex,
    /// Total degree, ties broken reverse-lexicographically. Usually the
    /// fastest.
    GradedRevLex,
}

impl Order {
    /// Compares two exponent vectors of equal length.
    #[must_use]
    pub fn cmp(
        self,
        a: &[u32],
        b: &[u32],
    ) -> Ordering {
        let total = |m: &[u32]| m.iter().map(|&e| u64::from(e)).sum::<u64>();
        match self {
            | Self::Lex => a.cmp(b),
            | Self::GradedLex => total(a).cmp(&total(b)).then_with(|| a.cmp(b)),
            | Self::GradedRevLex => total(a).cmp(&total(b)).then_with(|| {
                // The monomial with the smaller exponent in the last
                // variable where they differ is the larger one.
                a.iter().zip(b).rev().find(|(x, y)| x != y).map_or(Ordering::Equal, |(x, y)| y.cmp(x))
            }),
        }
    }
}

/// A term: exponent vector and non-zero coefficient.
pub type Term = (Vec<u32>, BigRational);

/// A polynomial as terms in strictly descending monomial order.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GPoly {
    terms: Vec<Term>,
}

impl GPoly {
    /// Builds a polynomial from arbitrary terms: sorts them by `order`,
    /// merges equal monomials and drops zeros.
    #[must_use]
    pub fn new(
        mut terms: Vec<Term>,
        order: Order,
    ) -> Self {
        terms.sort_by(|a, b| order.cmp(&b.0, &a.0));
        let mut merged: Vec<Term> = Vec::with_capacity(terms.len());
        for (mono, coeff) in terms {
            match merged.last_mut() {
                | Some((last, acc)) if *last == mono => *acc += coeff,
                | _ => merged.push((mono, coeff)),
            }
        }
        merged.retain(|(_, c)| !c.is_zero());
        Self { terms: merged }
    }

    /// The terms in descending order.
    #[must_use]
    pub fn terms(&self) -> &[Term] {
        &self.terms
    }

    /// Whether this is the zero polynomial.
    #[must_use]
    pub const fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    fn leading(&self) -> Option<&Term> {
        self.terms.first()
    }

    fn monic(&self) -> Self {
        match self.leading() {
            | Some((_, lead)) => {
                Self { terms: self.terms.iter().map(|(m, c)| (m.clone(), c / lead)).collect() }
            },
            | None => self.clone(),
        }
    }

    /// `self - factor * x^shift * other`.
    fn sub_scaled(
        &self,
        other: &Self,
        shift: &[u32],
        factor: &BigRational,
        order: Order,
    ) -> Self {
        let mut out = Vec::with_capacity(self.terms.len() + other.terms.len());
        let (mut i, mut j) = (0, 0);
        let shifted = |t: &Term| -> Term {
            (t.0.iter().zip(shift).map(|(a, b)| a + b).collect(), -(&t.1 * factor))
        };
        let mut pending = other.terms.first().map(&shifted);
        loop {
            match (self.terms.get(i), &pending) {
                | (Some(a), Some(b)) => match order.cmp(&a.0, &b.0) {
                    | Ordering::Greater => {
                        out.push(a.clone());
                        i += 1;
                    },
                    | Ordering::Less => {
                        out.push(b.clone());
                        j += 1;
                        pending = other.terms.get(j).map(&shifted);
                    },
                    | Ordering::Equal => {
                        let sum = &a.1 + &b.1;
                        if !sum.is_zero() {
                            out.push((a.0.clone(), sum));
                        }
                        i += 1;
                        j += 1;
                        pending = other.terms.get(j).map(&shifted);
                    },
                },
                | (Some(a), None) => {
                    out.push(a.clone());
                    i += 1;
                },
                | (None, Some(b)) => {
                    out.push(b.clone());
                    j += 1;
                    pending = other.terms.get(j).map(&shifted);
                },
                | (None, None) => return Self { terms: out },
            }
        }
    }
}

fn divides(
    a: &[u32],
    b: &[u32],
) -> bool {
    a.iter().zip(b).all(|(x, y)| x <= y)
}

fn lcm(
    a: &[u32],
    b: &[u32],
) -> Vec<u32> {
    a.iter().zip(b).map(|(&x, &y)| x.max(y)).collect()
}

/// The normal form of `p` modulo `basis`: no term of the result is
/// divisible by a leading monomial of the basis.
#[must_use]
pub fn reduce(
    p: &GPoly,
    basis: &[GPoly],
    order: Order,
) -> GPoly {
    let mut work = p.clone();
    let mut done: Vec<Term> = Vec::new();
    while let Some((mono, coeff)) = work.leading().cloned() {
        let reducer = basis.iter().find_map(|g| {
            let (lead_mono, lead_coeff) = g.leading()?;
            divides(lead_mono, &mono).then_some((g, lead_mono, lead_coeff))
        });
        match reducer {
            | Some((g, lead_mono, lead_coeff)) => {
                let shift: Vec<u32> = mono.iter().zip(lead_mono).map(|(a, b)| a - b).collect();
                work = work.sub_scaled(g, &shift, &(&coeff / lead_coeff), order);
            },
            | None => {
                done.push((mono, coeff));
                work.terms.remove(0);
            },
        }
    }
    GPoly { terms: done }
}

fn s_polynomial(
    f: &GPoly,
    g: &GPoly,
    order: Order,
) -> GPoly {
    let (Some((fm, fc)), Some((gm, gc))) = (f.leading(), g.leading()) else {
        return GPoly { terms: Vec::new() };
    };
    let common = lcm(fm, gm);
    let shift_f: Vec<u32> = common.iter().zip(fm).map(|(a, b)| a - b).collect();
    let shift_g: Vec<u32> = common.iter().zip(gm).map(|(a, b)| a - b).collect();
    let zero = GPoly { terms: Vec::new() };
    // (1/fc) x^sf f - (1/gc) x^sg g
    let left = zero.sub_scaled(f, &shift_f, &(-BigRational::one() / fc), order);
    left.sub_scaled(g, &shift_g, &(BigRational::one() / gc), order)
}

/// Limits for [`groebner`].
#[derive(Copy, Clone, Debug)]
pub struct GroebnerLimits {
    /// Maximum number of S-polynomial reductions.
    pub reductions: usize,
    /// Maximum size of the intermediate basis.
    pub basis: usize,
}

impl Default for GroebnerLimits {
    fn default() -> Self {
        Self { reductions: 20_000, basis: 500 }
    }
}

/// The reduced Gröbner basis of the ideal generated by `generators`.
///
/// Returns `None` when a limit is exceeded. The result is sorted by leading
/// monomial, ascending, and every element is monic.
#[must_use]
pub fn groebner(
    generators: &[GPoly],
    order: Order,
    limits: GroebnerLimits,
) -> Option<Vec<GPoly>> {
    let mut basis: Vec<GPoly> = generators.iter().filter(|g| !g.is_zero()).map(GPoly::monic).collect();
    let mut pairs: Vec<(usize, usize)> =
        (0..basis.len()).flat_map(|i| (0..i).map(move |j| (j, i))).collect();
    let mut reductions = 0_usize;
    while !pairs.is_empty() {
        // Normal strategy: the pair whose leading monomials have the
        // smallest least common multiple.
        let pick = (0..pairs.len()).min_by(|&a, &b| {
            let key = |(i, j): (usize, usize)| match (basis[i].leading(), basis[j].leading()) {
                | (Some((x, _)), Some((y, _))) => lcm(x, y),
                | _ => Vec::new(),
            };
            order.cmp(&key(pairs[a]), &key(pairs[b]))
        })?;
        let (i, j) = pairs.swap_remove(pick);
        let (Some((lead_i, _)), Some((lead_j, _))) = (basis[i].leading(), basis[j].leading()) else {
            continue;
        };
        // First criterion: coprime leading monomials reduce to zero.
        if lead_i.iter().zip(lead_j).all(|(&a, &b)| a == 0 || b == 0) {
            continue;
        }
        // Second criterion: some k whose leading monomial divides the lcm
        // and whose pairs with i and j have already been handled.
        let common = lcm(lead_i, lead_j);
        let covered = (0..basis.len()).any(|k| {
            k != i
                && k != j
                && basis[k].leading().is_some_and(|(m, _)| divides(m, &common))
                && !pairs.contains(&(i.min(k), i.max(k)))
                && !pairs.contains(&(j.min(k), j.max(k)))
        });
        if covered {
            continue;
        }
        reductions += 1;
        if reductions > limits.reductions {
            return None;
        }
        let remainder = reduce(&s_polynomial(&basis[i], &basis[j], order), &basis, order);
        if !remainder.is_zero() {
            if basis.len() >= limits.basis {
                return None;
            }
            let new = basis.len();
            pairs.extend((0..new).map(|k| (k, new)));
            basis.push(remainder.monic());
        }
    }

    // Minimal: drop elements whose leading monomial is divisible by another's.
    let mut minimal: Vec<GPoly> = Vec::new();
    for (index, g) in basis.iter().enumerate() {
        let Some((lead, _)) = g.leading() else {
            continue;
        };
        let redundant = basis.iter().enumerate().any(|(other, h)| {
            other != index
                && h.leading().is_some_and(|(m, _)| divides(m, lead) && (m != lead || other < index))
        });
        if !redundant {
            minimal.push(g.clone());
        }
    }
    // Reduced: each element fully reduced modulo the others.
    let mut reduced = Vec::with_capacity(minimal.len());
    for (index, g) in minimal.iter().enumerate() {
        let others: Vec<GPoly> =
            minimal.iter().enumerate().filter(|&(k, _)| k != index).map(|(_, h)| h.clone()).collect();
        reduced.push(reduce(g, &others, order).monic());
    }
    reduced.sort_by(|a, b| match (a.leading(), b.leading()) {
        | (Some((x, _)), Some((y, _))) => order.cmp(x, y),
        | _ => Ordering::Equal,
    });
    Some(reduced)
}

#[cfg(test)]
mod tests {
    use num_bigint::BigInt;

    use super::*;

    fn q(n: i64) -> BigRational {
        BigRational::from_integer(BigInt::from(n))
    }

    /// Builds a polynomial in (x, y, z) from `(coeff, [ex, ey, ez])`.
    fn p(
        terms: &[(i64, [u32; 3])],
        order: Order,
    ) -> GPoly {
        GPoly::new(terms.iter().map(|(c, e)| (e.to_vec(), q(*c))).collect(), order)
    }

    #[test]
    fn orders() {
        let (a, b): (&[u32], &[u32]) = (&[1, 2, 0], &[0, 3, 4]);
        assert_eq!(Order::Lex.cmp(a, b), Ordering::Greater, "x*y^2 > y^3*z^4 in lex");
        assert_eq!(Order::GradedLex.cmp(a, b), Ordering::Less);
        let (c, d): (&[u32], &[u32]) = (&[1, 1, 1], &[0, 3, 0]);
        assert_eq!(Order::GradedLex.cmp(c, d), Ordering::Greater, "x*y*z > y^3 in grlex");
        assert_eq!(Order::GradedRevLex.cmp(c, d), Ordering::Less, "y^3 > x*y*z in grevlex");
    }

    #[test]
    fn construction_merges_and_sorts() {
        let poly = p(&[(1, [0, 0, 0]), (2, [1, 0, 0]), (3, [1, 0, 0]), (4, [0, 1, 0]), (-4, [0, 1, 0])], Order::Lex);
        assert_eq!(poly.terms(), &[(vec![1, 0, 0], q(5)), (vec![0, 0, 0], q(1))]);
    }

    #[test]
    fn linear_system() {
        // x + y - 3, x - y - 1  ->  x - 2, y - 1
        let order = Order::Lex;
        let basis = groebner(
            &[
                p(&[(1, [1, 0, 0]), (1, [0, 1, 0]), (-3, [0, 0, 0])], order),
                p(&[(1, [1, 0, 0]), (-1, [0, 1, 0]), (-1, [0, 0, 0])], order),
            ],
            order,
            GroebnerLimits::default(),
        )
        .unwrap_or_default();
        assert_eq!(basis, vec![
            p(&[(1, [0, 1, 0]), (-1, [0, 0, 0])], order),
            p(&[(1, [1, 0, 0]), (-2, [0, 0, 0])], order)
        ]);
    }

    #[test]
    fn circle_and_line_eliminate_to_a_univariate_polynomial() {
        // x^2 + y^2 - 1, x - y  ->  lex basis {y^2 - 1/2, x - y}
        let order = Order::Lex;
        let basis = groebner(
            &[
                p(&[(1, [2, 0, 0]), (1, [0, 2, 0]), (-1, [0, 0, 0])], order),
                p(&[(1, [1, 0, 0]), (-1, [0, 1, 0])], order),
            ],
            order,
            GroebnerLimits::default(),
        )
        .unwrap_or_default();
        let half = BigRational::new(BigInt::from(-1), BigInt::from(2));
        assert_eq!(basis.len(), 2);
        assert_eq!(basis[0].terms(), &[(vec![0, 2, 0], q(1)), (vec![0, 0, 0], half)]);
        assert_eq!(basis[1], p(&[(1, [1, 0, 0]), (-1, [0, 1, 0])], order));
    }

    #[test]
    fn inconsistent_system_gives_the_unit_ideal() {
        let order = Order::GradedRevLex;
        let basis = groebner(
            &[p(&[(1, [1, 0, 0])], order), p(&[(1, [1, 0, 0]), (1, [0, 0, 0])], order)],
            order,
            GroebnerLimits::default(),
        )
        .unwrap_or_default();
        assert_eq!(basis, vec![p(&[(1, [0, 0, 0])], order)]);
    }

    #[test]
    fn the_reduced_basis_does_not_depend_on_the_generators() {
        // Cox, Little, O'Shea: <x^3 - 2xy, x^2y - 2y^2 + x> in grlex.
        let order = Order::GradedLex;
        let f1 = p(&[(1, [3, 0, 0]), (-2, [1, 1, 0])], order);
        let f2 = p(&[(1, [2, 1, 0]), (-2, [0, 2, 0]), (1, [1, 0, 0])], order);
        let basis = groebner(&[f1.clone(), f2.clone()], order, GroebnerLimits::default()).unwrap_or_default();
        // Known reduced basis: {x^2, xy, y^2 - x/2}.
        let half = BigRational::new(BigInt::from(-1), BigInt::from(2));
        let expected = vec![
            GPoly::new(vec![(vec![0, 2, 0], q(1)), (vec![1, 0, 0], half)], order),
            p(&[(1, [1, 1, 0])], order),
            p(&[(1, [2, 0, 0])], order),
        ];
        assert_eq!(basis, expected);
        // Same ideal, different generators and order of input.
        let sum = GPoly::new([f1.terms(), f2.terms()].concat(), order);
        let again = groebner(&[f2, sum, f1], order, GroebnerLimits::default()).unwrap_or_default();
        assert_eq!(again, expected);
        // Every generator reduces to zero modulo the basis.
        for g in &expected {
            assert!(reduce(g, &basis, order).is_zero());
        }
    }

    #[test]
    fn limits_abort_instead_of_running_away() {
        let order = Order::Lex;
        let gens = [
            p(&[(1, [5, 1, 0]), (1, [0, 3, 2]), (-1, [0, 0, 0])], order),
            p(&[(1, [2, 4, 1]), (1, [1, 0, 5]), (3, [0, 1, 0])], order),
            p(&[(1, [3, 0, 3]), (1, [0, 5, 1]), (-2, [1, 1, 1])], order),
        ];
        assert_eq!(groebner(&gens, order, GroebnerLimits { reductions: 3, basis: 500 }), None);
    }
}
