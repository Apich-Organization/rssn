//! Gröbner bases over the rationals (Buchberger's algorithm).
//!
//! Polynomials are stored densely sorted by the chosen monomial order, with
//! exponent vectors over a fixed list of variables. The implementation uses
//! the normal selection strategy (smallest least common multiple first) and
//! Buchberger's two criteria, and returns the *reduced* basis, which is
//! unique for a given ideal and order — so two descriptions of the same
//! ideal reduce to literally the same list.
//!
//! [`fglm`] converts a reduced basis of a zero-dimensional ideal from one
//! order to another by linear algebra in the finite-dimensional quotient
//! ring (Faugère–Gianni–Lazard–Mora), and [`groebner_via`] uses it for the
//! lex and graded-lex targets: the basis is computed in the cheap graded
//! reverse-lexicographic order and converted, which is how lexicographic
//! elimination (triangular systems) is done for zero-dimensional ideals.
//! Ideals that are not zero-dimensional, and quotients beyond
//! [`FGLM_MAX_DIMENSION`], fall back to Buchberger in the requested order.
//!
//! Parallelism: with several rayon threads the S-pairs of the smallest total degree are reduced in
//! batches against one basis snapshot (the remainders are inserted afterwards in pair order, after
//! a further reduction against the elements added in the batch). The reduced basis is unique, so it
//! is the same list for every thread count; with one thread the batch size is one. [`fglm`] stays
//! sequential (its multiplication columns are not the bottleneck).

use std::cmp::Ordering;
use std::collections::HashMap;
use std::collections::HashSet;

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
            | Self::GradedRevLex => {
                total(a).cmp(&total(b)).then_with(|| {
                    // The monomial with the smaller exponent in the last
                    // variable where they differ is the larger one.
                    a.iter()
                        .zip(b)
                        .rev()
                        .find(|(x, y)| x != y)
                        .map_or(Ordering::Equal, |(x, y)| y.cmp(x))
                })
            },
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

    /// The leading term, if any.
    #[must_use]
    pub fn leading(&self) -> Option<&Term> {
        self.terms.first()
    }

    fn monic(&self) -> Self {
        match self.leading() {
            | Some((_, lead)) => Self {
                terms: self
                    .terms
                    .iter()
                    .map(|(m, c)| (m.clone(), c / lead))
                    .collect(),
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
            (
                t.0.iter().zip(shift).map(|(a, b)| a + b).collect(),
                -(&t.1 * factor),
            )
        };
        let mut pending = other.terms.first().map(shifted);
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
                        pending = other.terms.get(j).map(shifted);
                    },
                    | Ordering::Equal => {
                        let sum = &a.1 + &b.1;
                        if !sum.is_zero() {
                            out.push((a.0.clone(), sum));
                        }
                        i += 1;
                        j += 1;
                        pending = other.terms.get(j).map(shifted);
                    },
                },
                | (Some(a), None) => {
                    out.push(a.clone());
                    i += 1;
                },
                | (None, Some(b)) => {
                    out.push(b.clone());
                    j += 1;
                    pending = other.terms.get(j).map(shifted);
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

/// S-pairs per thread reduced together in one batch of [`groebner`].
const BATCH_PER_THREAD: usize = 2;

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
        Self {
            reductions: 20_000,
            basis: 500,
        }
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
    let mut basis: Vec<GPoly> = generators
        .iter()
        .filter(|g| !g.is_zero())
        .map(GPoly::monic)
        .collect();
    let mut pairs: Vec<(usize, usize)> = (0..basis.len())
        .flat_map(|i| (0..i).map(move |j| (j, i)))
        .collect();
    let mut reductions = 0_usize;
    // With several threads the pairs of the smallest total degree are reduced
    // in batches against the same basis snapshot; with one thread the batch
    // size is 1 and this is the classic one-pair-at-a-time loop. The reduced
    // basis is unique, so it is the same either way.
    let cap = if rayon::current_num_threads() > 1 {
        BATCH_PER_THREAD * rayon::current_num_threads()
    } else {
        1
    };
    while !pairs.is_empty() {
        let mut batch: Vec<(usize, usize)> = Vec::new();
        let mut batch_degree: Option<u64> = None;
        while !pairs.is_empty() && batch.len() < cap {
            // Normal strategy: the pair whose leading monomials have the
            // smallest least common multiple.
            let pick = (0..pairs.len()).min_by(|&a, &b| {
                let key = |(i, j): (usize, usize)| match (basis[i].leading(), basis[j].leading()) {
                    | (Some((x, _)), Some((y, _))) => lcm(x, y),
                    | _ => Vec::new(),
                };
                order.cmp(&key(pairs[a]), &key(pairs[b]))
            })?;
            let (i, j) = pairs[pick];
            let (Some((lead_i, _)), Some((lead_j, _))) = (basis[i].leading(), basis[j].leading())
            else {
                pairs.swap_remove(pick);
                continue;
            };
            let common = lcm(lead_i, lead_j);
            let degree: u64 = common.iter().map(|&e| u64::from(e)).sum();
            if batch_degree.is_some_and(|d| d != degree) {
                break;
            }
            pairs.swap_remove(pick);
            // First criterion: coprime leading monomials reduce to zero.
            if lead_i.iter().zip(lead_j).all(|(&a, &b)| a == 0 || b == 0) {
                continue;
            }
            // Second criterion: some k whose leading monomial divides the lcm
            // and whose pairs with i and j have already been handled.
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
            batch_degree = Some(degree);
            batch.push((i, j));
        }
        if batch.is_empty() {
            continue;
        }
        let snapshot = basis.len();
        let remainders: Vec<GPoly> = if batch.len() >= 2 {
            let basis_ref = &basis;
            crate::kernels::par::map_slice(&batch, 2, |&(i, j)| {
                reduce(
                    &s_polynomial(&basis_ref[i], &basis_ref[j], order),
                    basis_ref,
                    order,
                )
            })
        } else {
            batch
                .iter()
                .map(|&(i, j)| reduce(&s_polynomial(&basis[i], &basis[j], order), &basis, order))
                .collect()
        };
        for remainder in remainders {
            // elements added earlier in this batch may reduce it further
            let remainder = if basis.len() > snapshot {
                reduce(&remainder, &basis, order)
            } else {
                remainder
            };
            if !remainder.is_zero() {
                if basis.len() >= limits.basis {
                    return None;
                }
                let new = basis.len();
                pairs.extend((0..new).map(|k| (k, new)));
                basis.push(remainder.monic());
            }
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
                && h.leading()
                    .is_some_and(|(m, _)| divides(m, lead) && (m != lead || other < index))
        });
        if !redundant {
            minimal.push(g.clone());
        }
    }
    // Reduced: each element fully reduced modulo the others.
    let mut reduced = Vec::with_capacity(minimal.len());
    for (index, g) in minimal.iter().enumerate() {
        let others: Vec<GPoly> = minimal
            .iter()
            .enumerate()
            .filter(|&(k, _)| k != index)
            .map(|(_, h)| h.clone())
            .collect();
        reduced.push(reduce(g, &others, order).monic());
    }
    reduced.sort_by(|a, b| match (a.leading(), b.leading()) {
        | (Some((x, _)), Some((y, _))) => order.cmp(x, y),
        | _ => Ordering::Equal,
    });
    Some(reduced)
}

/// Largest dimension of the quotient ring [`fglm`] works in.
pub const FGLM_MAX_DIMENSION: usize = 400;

/// A row of the incremental echelon form used by [`fglm`]: the reduced
/// vector (pivot entry `1`) and the combination of the staircase vectors
/// that equals it.
struct EchelonRow {
    vector: Vec<BigRational>,
    combination: Vec<BigRational>,
    pivot: usize,
}

/// The reduced Gröbner basis in the order `to` of a zero-dimensional ideal
/// given by its reduced Gröbner basis `basis` in the order `from`
/// (FGLM).
///
/// The standard monomials of `basis` span the quotient ring; the matrices
/// of multiplication by the variables are read off by reducing, and the
/// monomials are then walked in the order `to`: a monomial whose normal
/// form depends linearly on the smaller ones gives an element of the new
/// basis, an independent one joins the staircase and contributes its
/// multiples. `None` when the ideal is not zero-dimensional, `basis` is
/// empty, or the quotient is larger than [`FGLM_MAX_DIMENSION`].
#[must_use]
pub fn fglm(
    basis: &[GPoly],
    from: Order,
    to: Order,
) -> Option<Vec<GPoly>> {
    let leads: Vec<&Vec<u32>> = basis
        .iter()
        .filter_map(|g| g.leading().map(|(m, _)| m))
        .collect();
    let n = leads.first()?.len();
    if leads.iter().any(|m| m.iter().all(|&e| e == 0)) {
        return Some(vec![GPoly::new(vec![(vec![0; n], BigRational::one())], to)]);
    }
    // Zero-dimensional: a pure power of every variable leads some element.
    let pure_power = |i: usize| {
        leads.iter().any(|m| {
            m.get(i).is_some_and(|&e| e > 0) && m.iter().enumerate().all(|(j, &e)| j == i || e == 0)
        })
    };
    if !(0..n).all(pure_power) {
        return None;
    }
    // The standard monomials, found from 1 by multiplying with the variables.
    let unit_monomial = vec![0_u32; n];
    let mut standard: Vec<Vec<u32>> = Vec::new();
    let mut index: HashMap<Vec<u32>, usize> = HashMap::new();
    let mut frontier = vec![unit_monomial.clone()];
    while let Some(m) = frontier.pop() {
        if index.contains_key(&m) || leads.iter().any(|l| divides(l, &m)) {
            continue;
        }
        if standard.len() >= FGLM_MAX_DIMENSION {
            return None;
        }
        index.insert(m.clone(), standard.len());
        standard.push(m.clone());
        for i in 0..n {
            let mut next = m.clone();
            *next.get_mut(i)? += 1;
            frontier.push(next);
        }
    }
    let dimension = standard.len();
    // Multiplication by x_i on the standard monomials, as sparse columns.
    let mut multiplication: Vec<Vec<Vec<(usize, BigRational)>>> = Vec::with_capacity(n);
    for i in 0..n {
        let mut columns = Vec::with_capacity(dimension);
        for m in &standard {
            let mut product = m.clone();
            *product.get_mut(i)? += 1;
            let column = match index.get(&product) {
                | Some(&position) => vec![(position, BigRational::one())],
                | None => {
                    let normal = reduce(
                        &GPoly {
                            terms: vec![(product, BigRational::one())],
                        },
                        basis,
                        from,
                    );
                    let mut entries = Vec::with_capacity(normal.terms.len());
                    for (monomial, coefficient) in &normal.terms {
                        entries.push((*index.get(monomial)?, coefficient.clone()));
                    }
                    entries
                },
            };
            columns.push(column);
        }
        multiplication.push(columns);
    }
    let multiply = |i: usize, v: &[BigRational]| -> Option<Vec<BigRational>> {
        let mut out = vec![BigRational::zero(); dimension];
        for (j, value) in v.iter().enumerate() {
            if value.is_zero() {
                continue;
            }
            for (r, coefficient) in multiplication.get(i)?.get(j)? {
                *out.get_mut(*r)? += coefficient * value;
            }
        }
        Some(out)
    };
    // Walk the monomials in the target order.
    let mut unit_vector = vec![BigRational::zero(); dimension];
    *unit_vector.get_mut(*index.get(&unit_monomial)?)? = BigRational::one();
    let mut queue: Vec<(Vec<u32>, Vec<BigRational>)> = vec![(unit_monomial.clone(), unit_vector)];
    let mut seen: HashSet<Vec<u32>> = HashSet::from([unit_monomial]);
    let mut staircase: Vec<Vec<u32>> = Vec::new();
    let mut rows: Vec<EchelonRow> = Vec::new();
    let mut new_leads: Vec<Vec<u32>> = Vec::new();
    let mut result: Vec<GPoly> = Vec::new();
    while !queue.is_empty() {
        let smallest = (0..queue.len()).min_by(|&a, &b| {
            let key = |k: usize| queue.get(k).map(|e| e.0.as_slice());
            match (key(a), key(b)) {
                | (Some(x), Some(y)) => to.cmp(x, y),
                | _ => Ordering::Equal,
            }
        })?;
        let (monomial, vector) = queue.swap_remove(smallest);
        if new_leads.iter().any(|l| divides(l, &monomial)) {
            continue;
        }
        // Reduce against the echelon rows, tracking the combination.
        let mut rest = vector.clone();
        let mut coefficients = vec![BigRational::zero(); staircase.len()];
        for row in &rows {
            let c = rest.get(row.pivot)?.clone();
            if c.is_zero() {
                continue;
            }
            for (slot, value) in rest.iter_mut().zip(&row.vector) {
                *slot -= &c * value;
            }
            for (slot, value) in coefficients.iter_mut().zip(&row.combination) {
                *slot += &c * value;
            }
        }
        if let Some(pivot) = rest.iter().position(|v| !v.is_zero()) {
            // Independent: it joins the staircase and contributes its multiples.
            let scale = rest.get(pivot)?.clone();
            let mut combination: Vec<BigRational> =
                coefficients.iter().map(|c| -c / &scale).collect();
            combination.push(BigRational::one() / &scale);
            rows.push(EchelonRow {
                vector: rest.iter().map(|v| v / &scale).collect(),
                combination,
                pivot,
            });
            for i in 0..n {
                let mut next = monomial.clone();
                *next.get_mut(i)? += 1;
                if seen.insert(next.clone()) {
                    queue.push((next, multiply(i, &vector)?));
                }
            }
            staircase.push(monomial);
        } else {
            // monomial = Σ coefficients_j staircase_j modulo the ideal.
            let mut terms: Vec<Term> = vec![(monomial.clone(), BigRational::one())];
            for (s, c) in staircase.iter().zip(&coefficients) {
                if !c.is_zero() {
                    terms.push((s.clone(), -c.clone()));
                }
            }
            result.push(GPoly::new(terms, to));
            new_leads.push(monomial);
        }
    }
    result.sort_by(|a, b| match (a.leading(), b.leading()) {
        | (Some((x, _)), Some((y, _))) => to.cmp(x, y),
        | _ => Ordering::Equal,
    });
    Some(result)
}

/// The reduced Gröbner basis of `generators` in `order`.
///
/// For the lex and graded-lex orders the basis is computed in the graded reverse-lex order
/// and converted by [`fglm`] when the ideal is zero-dimensional (usually
/// far less work than Buchberger in the target order); otherwise, and for
/// graded reverse-lex itself, this is [`groebner`]. The result is the same
/// reduced basis either way.
#[must_use]
pub fn groebner_via(
    generators: &[GPoly],
    order: Order,
    limits: GroebnerLimits,
) -> Option<Vec<GPoly>> {
    if order != Order::GradedRevLex {
        let graded: Vec<GPoly> = generators
            .iter()
            .map(|g| GPoly::new(g.terms().to_vec(), Order::GradedRevLex))
            .collect();
        if let Some(converted) = groebner(&graded, Order::GradedRevLex, limits)
            .and_then(|basis| fglm(&basis, Order::GradedRevLex, order))
        {
            return Some(converted);
        }
    }
    groebner(generators, order, limits)
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
        GPoly::new(
            terms.iter().map(|(c, e)| (e.to_vec(), q(*c))).collect(),
            order,
        )
    }

    #[test]
    fn orders() {
        let (a, b): (&[u32], &[u32]) = (&[1, 2, 0], &[0, 3, 4]);
        assert_eq!(
            Order::Lex.cmp(a, b),
            Ordering::Greater,
            "x*y^2 > y^3*z^4 in lex"
        );
        assert_eq!(Order::GradedLex.cmp(a, b), Ordering::Less);
        let (c, d): (&[u32], &[u32]) = (&[1, 1, 1], &[0, 3, 0]);
        assert_eq!(
            Order::GradedLex.cmp(c, d),
            Ordering::Greater,
            "x*y*z > y^3 in grlex"
        );
        assert_eq!(
            Order::GradedRevLex.cmp(c, d),
            Ordering::Less,
            "y^3 > x*y*z in grevlex"
        );
    }

    #[test]
    fn construction_merges_and_sorts() {
        let poly = p(
            &[
                (1, [0, 0, 0]),
                (2, [1, 0, 0]),
                (3, [1, 0, 0]),
                (4, [0, 1, 0]),
                (-4, [0, 1, 0]),
            ],
            Order::Lex,
        );
        assert_eq!(
            poly.terms(),
            &[(vec![1, 0, 0], q(5)), (vec![0, 0, 0], q(1))]
        );
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
        assert_eq!(
            basis,
            vec![
                p(&[(1, [0, 1, 0]), (-1, [0, 0, 0])], order),
                p(&[(1, [1, 0, 0]), (-2, [0, 0, 0])], order)
            ]
        );
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
        assert_eq!(
            basis[0].terms(),
            &[(vec![0, 2, 0], q(1)), (vec![0, 0, 0], half)]
        );
        assert_eq!(basis[1], p(&[(1, [1, 0, 0]), (-1, [0, 1, 0])], order));
    }

    #[test]
    fn inconsistent_system_gives_the_unit_ideal() {
        let order = Order::GradedRevLex;
        let basis = groebner(
            &[
                p(&[(1, [1, 0, 0])], order),
                p(&[(1, [1, 0, 0]), (1, [0, 0, 0])], order),
            ],
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
        let basis = groebner(&[f1.clone(), f2.clone()], order, GroebnerLimits::default())
            .unwrap_or_default();
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
        assert_eq!(
            groebner(
                &gens,
                order,
                GroebnerLimits {
                    reductions: 3,
                    basis: 500
                }
            ),
            None
        );
    }

    /// A polynomial in `n` variables from `(coeff, exponents)`.
    fn poly_n(
        terms: &[(i64, Vec<u32>)],
        order: Order,
    ) -> GPoly {
        GPoly::new(
            terms.iter().map(|(c, e)| (e.clone(), q(*c))).collect(),
            order,
        )
    }

    /// `x_i`-power vector of length `n`.
    fn mono(
        n: usize,
        powers: &[(usize, u32)],
    ) -> Vec<u32> {
        let mut m = vec![0; n];
        for &(i, e) in powers {
            m[i] += e;
        }
        m
    }

    /// Katsura-`m` in `m + 1` variables: `x_0 + 2 Σ x_l = 1` and
    /// `Σ_{l=-m}^{m} x_{|l|} x_{|k-l|} = x_k` for `k = 0 .. m-1`.
    fn katsura(
        m: usize,
        order: Order,
    ) -> Vec<GPoly> {
        let n = m + 1;
        let mut out = Vec::new();
        let mut linear = vec![(1, mono(n, &[(0, 1)])), (-1, mono(n, &[]))];
        for l in 1..=m {
            linear.push((2, mono(n, &[(l, 1)])));
        }
        out.push(poly_n(&linear, order));
        for k in 0..m {
            let mut terms = vec![(-1, mono(n, &[(k, 1)]))];
            for l in -(m as i64)..=(m as i64) {
                let (a, b) = (
                    l.unsigned_abs() as usize,
                    (k as i64 - l).unsigned_abs() as usize,
                );
                if a <= m && b <= m {
                    terms.push((1, mono(n, &[(a, 1), (b, 1)])));
                }
            }
            out.push(poly_n(&terms, order));
        }
        out
    }

    /// Cyclic-`n`: the elementary cyclic sums minus one for the last.
    fn cyclic(
        n: usize,
        order: Order,
    ) -> Vec<GPoly> {
        let mut out = Vec::new();
        for d in 1..n {
            let terms: Vec<(i64, Vec<u32>)> = (0..n)
                .map(|start| {
                    (
                        1,
                        mono(n, &(0..d).map(|j| ((start + j) % n, 1)).collect::<Vec<_>>()),
                    )
                })
                .collect();
            out.push(poly_n(&terms, order));
        }
        out.push(poly_n(&[(1, vec![1; n]), (-1, mono(n, &[]))], order));
        out
    }

    #[test]
    fn fglm_agrees_with_buchberger_in_the_target_order() {
        for (name, gens) in [
            (
                "circle and line",
                vec![
                    poly_n(
                        &[(1, vec![2, 0]), (1, vec![0, 2]), (-1, vec![0, 0])],
                        Order::GradedRevLex,
                    ),
                    poly_n(&[(1, vec![1, 0]), (-1, vec![0, 1])], Order::GradedRevLex),
                ],
            ),
            ("cyclic-3", cyclic(3, Order::GradedRevLex)),
            ("katsura-2", katsura(2, Order::GradedRevLex)),
            ("katsura-3", katsura(3, Order::GradedRevLex)),
            (
                "multiple root",
                vec![
                    poly_n(&[(1, vec![2, 0]), (1, vec![0, 1])], Order::GradedRevLex),
                    poly_n(&[(1, vec![0, 2])], Order::GradedRevLex),
                ],
            ),
            (
                "sphere, plane and a quadric",
                vec![
                    poly_n(
                        &[
                            (1, vec![2, 0, 0]),
                            (1, vec![0, 2, 0]),
                            (1, vec![0, 0, 2]),
                            (-1, vec![0, 0, 0]),
                        ],
                        Order::GradedRevLex,
                    ),
                    poly_n(
                        &[
                            (1, vec![1, 0, 0]),
                            (1, vec![0, 1, 0]),
                            (1, vec![0, 0, 1]),
                            (-1, vec![0, 0, 0]),
                        ],
                        Order::GradedRevLex,
                    ),
                    poly_n(
                        &[(1, vec![1, 1, 0]), (-1, vec![0, 0, 1])],
                        Order::GradedRevLex,
                    ),
                ],
            ),
        ] {
            let grevlex = groebner(&gens, Order::GradedRevLex, GroebnerLimits::default()).unwrap();
            for target in [Order::Lex, Order::GradedLex] {
                let converted = fglm(&grevlex, Order::GradedRevLex, target)
                    .unwrap_or_else(|| panic!("{name}: FGLM declined"));
                let regenerated: Vec<GPoly> = gens
                    .iter()
                    .map(|g| GPoly::new(g.terms().to_vec(), target))
                    .collect();
                let direct = groebner(&regenerated, target, GroebnerLimits::default()).unwrap();
                assert_eq!(converted, direct, "{name}, {target:?}");
                assert_eq!(
                    groebner_via(&regenerated, target, GroebnerLimits::default()).unwrap(),
                    direct,
                    "{name} via"
                );
            }
        }
    }

    #[test]
    fn fglm_declines_what_is_not_zero_dimensional_and_handles_the_unit_ideal() {
        // x*y - 1 has a hyperbola of solutions; (x) in two variables a line.
        let order = Order::GradedRevLex;
        let hyperbola = groebner(
            &[poly_n(&[(1, vec![1, 1]), (-1, vec![0, 0])], order)],
            order,
            GroebnerLimits::default(),
        )
        .unwrap();
        assert_eq!(fglm(&hyperbola, order, Order::Lex), None);
        let line = groebner(
            &[poly_n(&[(1, vec![1, 0])], order)],
            order,
            GroebnerLimits::default(),
        )
        .unwrap();
        assert_eq!(fglm(&line, order, Order::Lex), None);
        assert_eq!(fglm(&[], order, Order::Lex), None, "the zero ideal");
        // The fallback gives the Buchberger basis.
        let gens = [poly_n(&[(1, vec![1, 1]), (-1, vec![0, 0])], Order::Lex)];
        assert_eq!(
            groebner_via(&gens, Order::Lex, GroebnerLimits::default()),
            groebner(&gens, Order::Lex, GroebnerLimits::default())
        );
        // The unit ideal.
        let unit = groebner(
            &[
                poly_n(&[(1, vec![1, 0])], order),
                poly_n(&[(1, vec![1, 0]), (1, vec![0, 0])], order),
            ],
            order,
            GroebnerLimits::default(),
        )
        .unwrap();
        assert_eq!(
            fglm(&unit, order, Order::Lex),
            Some(vec![poly_n(&[(1, vec![0, 0])], Order::Lex)])
        );
    }

    #[test]
    fn fglm_converts_cyclic_5_with_70_solutions_and_reduces_the_generators_to_zero() {
        // Cyclic-5 has 70 solutions; Buchberger in lex is out of reach for
        // this implementation, the grevlex basis plus FGLM is immediate.
        let start = std::time::Instant::now();
        let gens = cyclic(5, Order::GradedRevLex);
        let grevlex = groebner(&gens, Order::GradedRevLex, GroebnerLimits::default()).unwrap();
        let lex = fglm(&grevlex, Order::GradedRevLex, Order::Lex).unwrap();
        eprintln!("cyclic-5: grevlex + FGLM in {:?}", start.elapsed());
        // The lex basis has 11 elements (it is not in shape position) and a
        // univariate element in the last variable.
        assert_eq!(lex.len(), 11);
        let univariate = lex
            .iter()
            .find(|g| {
                g.leading()
                    .is_some_and(|(m, _)| m[..4].iter().all(|&e| e == 0))
            })
            .unwrap();
        assert_eq!(univariate.leading().unwrap().0[4], 15);
        // The converted list is a Gröbner basis (Buchberger changes nothing)
        // and every generator lies in the ideal (it reduces to zero in lex).
        assert_eq!(
            groebner(&lex, Order::Lex, GroebnerLimits::default()).as_deref(),
            Some(lex.as_slice())
        );
        for g in &gens {
            let in_lex = GPoly::new(g.terms().to_vec(), Order::Lex);
            assert!(reduce(&in_lex, &lex, Order::Lex).is_zero());
        }
    }

    #[test]
    fn fglm_scales_to_katsura_4() {
        // Katsura-4: 16 solutions in five variables.
        let gens = katsura(4, Order::GradedRevLex);
        let start = std::time::Instant::now();
        let grevlex = groebner(&gens, Order::GradedRevLex, GroebnerLimits::default()).unwrap();
        let lex = fglm(&grevlex, Order::GradedRevLex, Order::Lex).unwrap();
        eprintln!("katsura-4: grevlex + FGLM in {:?}", start.elapsed());
        let univariate = lex
            .iter()
            .find(|g| {
                g.leading()
                    .is_some_and(|(m, _)| m[..4].iter().all(|&e| e == 0))
            })
            .unwrap();
        assert_eq!(univariate.leading().unwrap().0[4], 16);
        assert_eq!(
            groebner(&lex, Order::Lex, GroebnerLimits::default()).as_deref(),
            Some(lex.as_slice())
        );
        for g in &gens {
            let in_lex = GPoly::new(g.terms().to_vec(), Order::Lex);
            assert!(reduce(&in_lex, &lex, Order::Lex).is_zero());
        }
    }
}
