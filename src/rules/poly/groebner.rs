//! Gröbner bases over the rationals and over prime fields (Buchberger's
//! algorithm).
//!
//! Polynomials are stored densely sorted by the chosen monomial order, with
//! exponent vectors over a fixed list of variables. The implementation uses
//! the normal selection strategy (smallest least common multiple first) and
//! Buchberger's two criteria, and returns the *reduced* basis, which is
//! unique for a given ideal and order — so two descriptions of the same
//! ideal reduce to literally the same list.
//!
//! [`groebner`], [`reduce`], [`fglm`] and [`groebner_via`] run over `Q`.
//! Their `_in` counterparts ([`groebner_in`], [`reduce_in`], [`fglm_in`],
//! [`groebner_via_in`]) take a [`Domain`] and run the same exact algorithms
//! over a prime field `GF(p)`: every addition, multiplication and division is
//! modular, the pivot of the elimination is a modular inverse, and the result
//! is monic, reduced and sorted with the canonical representatives `0 .. p`.
//! The F4 and F5 operators use these entry points through
//! `groebner_f4(polys, vars[, order][, p])` and `groebner_f5(...)`; see
//! [`crate::rules::poly::field`] for the coefficient domain and its limits
//! (prime fields only — extension fields are out of scope).
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

use crate::rules::poly::field::Domain;

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
    /// Builds a rational polynomial from arbitrary terms: sorts them by
    /// `order`, merges equal monomials and drops zeros.
    #[must_use]
    pub fn new(
        terms: Vec<Term>,
        order: Order,
    ) -> Self {
        Self::new_in(&Domain::Rational, terms, order)
    }

    /// Builds a polynomial over `domain`: sorts the terms by `order`, merges
    /// equal monomials with the field addition and drops zeros.
    ///
    /// The coefficients must already be canonical elements of `domain` — for a
    /// prime field, integers in `0 .. p`. Use [`Self::in_domain`] to convert a
    /// polynomial with arbitrary rational coefficients; the engines only call
    /// this constructor with coefficients their own arithmetic produced.
    #[must_use]
    pub fn new_in(
        domain: &Domain,
        mut terms: Vec<Term>,
        order: Order,
    ) -> Self {
        terms.sort_by(|a, b| order.cmp(&b.0, &a.0));
        let mut merged: Vec<Term> = Vec::with_capacity(terms.len());
        for (mono, coeff) in terms {
            match merged.last_mut() {
                | Some((last, acc)) if *last == mono => *acc = domain.add(acc, &coeff),
                | _ => merged.push((mono, coeff)),
            }
        }
        merged.retain(|(_, c)| !domain.is_zero(c));
        Self { terms: merged }
    }

    /// This polynomial with every coefficient mapped into `domain`.
    ///
    /// For the rationals nothing changes. For a prime field the coefficient
    /// `a/b` becomes `a b^-1 mod p`, as the canonical integer in `0 .. p`;
    /// `None` when a denominator is divisible by `p`, i.e. when the
    /// coefficient has no image in the prime field.
    #[must_use]
    pub fn in_domain(
        &self,
        domain: &Domain,
        order: Order,
    ) -> Option<Self> {
        let mut terms = Vec::with_capacity(self.terms.len());
        for (monomial, coefficient) in &self.terms {
            terms.push((monomial.clone(), domain.normalize(coefficient)?));
        }
        Some(Self::new_in(domain, terms, order))
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

    fn monic_in(
        &self,
        domain: &Domain,
    ) -> Self {
        match self.leading() {
            | Some((_, lead)) => {
                let inverse = domain.inv(lead);
                match inverse {
                    | Some(inverse) => {
                        Self {
                            terms: self
                                .terms
                                .iter()
                                .map(|(m, c)| (m.clone(), domain.mul(c, &inverse)))
                                .collect(),
                        }
                    },
                    | None => Self { terms: Vec::new() },
                }
            },
            | None => self.clone(),
        }
    }

    /// `self - factor * x^shift * other`.
    fn sub_scaled(
        &self,
        domain: &Domain,
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
                domain.neg(&domain.mul(&t.1, factor)),
            )
        };
        let mut pending = other.terms.first().map(shifted);
        loop {
            match (self.terms.get(i), &pending) {
                | (Some(a), Some(b)) => {
                    match order.cmp(&a.0, &b.0) {
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
                            let sum = domain.add(&a.1, &b.1);
                            if !domain.is_zero(&sum) {
                                out.push((a.0.clone(), sum));
                            }
                            i += 1;
                            j += 1;
                            pending = other.terms.get(j).map(shifted);
                        },
                    }
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

/// The normal form of `p` modulo `basis` over the rationals: no term of the
/// result is divisible by a leading monomial of the basis.
#[must_use]
pub fn reduce(
    p: &GPoly,
    basis: &[GPoly],
    order: Order,
) -> GPoly {
    reduce_in(&Domain::Rational, p, basis, order)
}

/// [`reduce`] over `domain`.
#[must_use]
pub fn reduce_in(
    domain: &Domain,
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
                let factor = domain
                    .div(&coeff, lead_coeff)
                    .unwrap_or_else(|| domain.zero());
                work = work.sub_scaled(domain, g, &shift, &factor, order);
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
    domain: &Domain,
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
    let left_factor = domain.neg(
        &domain
            .div(&domain.one(), fc)
            .unwrap_or_else(|| domain.zero()),
    );
    let right_factor = domain
        .div(&domain.one(), gc)
        .unwrap_or_else(|| domain.zero());
    let left = zero.sub_scaled(domain, f, &shift_f, &left_factor, order);
    left.sub_scaled(domain, g, &shift_g, &right_factor, order)
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

/// The reduced Gröbner basis of the ideal generated by `generators` over the
/// rationals.
///
/// Returns `None` when a limit is exceeded. The result is sorted by leading
/// monomial, ascending, and every element is monic.
#[must_use]
pub fn groebner(
    generators: &[GPoly],
    order: Order,
    limits: GroebnerLimits,
) -> Option<Vec<GPoly>> {
    groebner_in(&Domain::Rational, generators, order, limits)
}

/// [`groebner`] over `domain`.
///
/// The generators must have canonical coefficients in `domain` (see
/// [`GPoly::in_domain`]); for a prime field the result is monic with the
/// canonical representatives `0 .. p`, reduced and sorted exactly like the
/// rational result.
#[must_use]
pub fn groebner_in(
    domain: &Domain,
    generators: &[GPoly],
    order: Order,
    limits: GroebnerLimits,
) -> Option<Vec<GPoly>> {
    let mut basis: Vec<GPoly> = generators
        .iter()
        .filter(|g| !g.is_zero())
        .map(|g| g.monic_in(domain))
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
                let key = |(i, j): (usize, usize)| {
                    match (basis[i].leading(), basis[j].leading()) {
                        | (Some((x, _)), Some((y, _))) => lcm(x, y),
                        | _ => Vec::new(),
                    }
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
                reduce_in(
                    domain,
                    &s_polynomial(domain, &basis_ref[i], &basis_ref[j], order),
                    basis_ref,
                    order,
                )
            })
        } else {
            batch
                .iter()
                .map(|&(i, j)| {
                    reduce_in(
                        domain,
                        &s_polynomial(domain, &basis[i], &basis[j], order),
                        &basis,
                        order,
                    )
                })
                .collect()
        };
        for remainder in remainders {
            // elements added earlier in this batch may reduce it further
            let remainder = if basis.len() > snapshot {
                reduce_in(domain, &remainder, &basis, order)
            } else {
                remainder
            };
            if !remainder.is_zero() {
                if basis.len() >= limits.basis {
                    return None;
                }
                let new = basis.len();
                pairs.extend((0..new).map(|k| (k, new)));
                basis.push(remainder.monic_in(domain));
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
        reduced.push(reduce_in(domain, g, &others, order).monic_in(domain));
    }
    reduced.sort_by(|a, b| {
        match (a.leading(), b.leading()) {
            | (Some((x, _)), Some((y, _))) => order.cmp(x, y),
            | _ => Ordering::Equal,
        }
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
/// given by its reduced Gröbner basis `basis` in the order `from` (FGLM),
/// over the rationals.
///
/// See [`fglm_in`] for the general coefficient domain. `None` when the ideal
/// is not zero-dimensional, `basis` is empty, or the quotient is larger than
/// [`FGLM_MAX_DIMENSION`].
#[must_use]
pub fn fglm(
    basis: &[GPoly],
    from: Order,
    to: Order,
) -> Option<Vec<GPoly>> {
    fglm_in(&Domain::Rational, basis, from, to)
}

/// [`fglm`] over `domain`.
///
/// The standard monomials of `basis` span the quotient ring; the matrices of
/// multiplication by the variables are read off by reducing, and the monomials
/// are then walked in the order `to`: a monomial whose normal form depends
/// linearly on the smaller ones gives an element of the new basis, an
/// independent one joins the staircase and contributes its multiples. Over a
/// prime field the echelon division is a modular inverse, so the conversion is
/// exact in `GF(p)`.
#[must_use]
pub fn fglm_in(
    domain: &Domain,
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
        return Some(vec![GPoly::new_in(
            domain,
            vec![(vec![0; n], domain.one())],
            to,
        )]);
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
                | Some(&position) => vec![(position, domain.one())],
                | None => {
                    let normal = reduce_in(
                        domain,
                        &GPoly {
                            terms: vec![(product, domain.one())],
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
        let mut out = vec![domain.zero(); dimension];
        for (j, value) in v.iter().enumerate() {
            if domain.is_zero(value) {
                continue;
            }
            for (r, coefficient) in multiplication.get(i)?.get(j)? {
                let slot = out.get_mut(*r)?;
                *slot = domain.add(slot, &domain.mul(coefficient, value));
            }
        }
        Some(out)
    };
    // Walk the monomials in the target order.
    let mut unit_vector = vec![domain.zero(); dimension];
    *unit_vector.get_mut(*index.get(&unit_monomial)?)? = domain.one();
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
        let mut coefficients = vec![domain.zero(); staircase.len()];
        for row in &rows {
            let c = rest.get(row.pivot)?.clone();
            if domain.is_zero(&c) {
                continue;
            }
            for (slot, value) in rest.iter_mut().zip(&row.vector) {
                *slot = domain.sub(slot, &domain.mul(&c, value));
            }
            for (slot, value) in coefficients.iter_mut().zip(&row.combination) {
                *slot = domain.add(slot, &domain.mul(&c, value));
            }
        }
        if let Some(pivot) = rest.iter().position(|v| !domain.is_zero(v)) {
            // Independent: it joins the staircase and contributes its multiples.
            let scale = rest.get(pivot)?.clone();
            let mut combination: Vec<BigRational> = coefficients
                .iter()
                .map(|c| {
                    domain
                        .div(&domain.neg(c), &scale)
                        .unwrap_or_else(|| domain.zero())
                })
                .collect();
            combination.push(
                domain
                    .div(&domain.one(), &scale)
                    .unwrap_or_else(|| domain.zero()),
            );
            rows.push(EchelonRow {
                vector: rest
                    .iter()
                    .map(|v| domain.div(v, &scale).unwrap_or_else(|| domain.zero()))
                    .collect(),
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
            let mut terms: Vec<Term> = vec![(monomial.clone(), domain.one())];
            for (s, c) in staircase.iter().zip(&coefficients) {
                if !domain.is_zero(c) {
                    terms.push((s.clone(), domain.neg(c)));
                }
            }
            result.push(GPoly::new_in(domain, terms, to));
            new_leads.push(monomial);
        }
    }
    result.sort_by(|a, b| {
        match (a.leading(), b.leading()) {
            | (Some((x, _)), Some((y, _))) => to.cmp(x, y),
            | _ => Ordering::Equal,
        }
    });
    Some(result)
}

/// The reduced Gröbner basis of `generators` in `order`, over the rationals.
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
    groebner_via_in(&Domain::Rational, generators, order, limits)
}

/// [`groebner_via`] over `domain`.
#[must_use]
pub fn groebner_via_in(
    domain: &Domain,
    generators: &[GPoly],
    order: Order,
    limits: GroebnerLimits,
) -> Option<Vec<GPoly>> {
    if order != Order::GradedRevLex {
        let graded: Vec<GPoly> = generators
            .iter()
            .map(|g| GPoly::new_in(domain, g.terms().to_vec(), Order::GradedRevLex))
            .collect();
        if let Some(converted) = groebner_in(domain, &graded, Order::GradedRevLex, limits)
            .and_then(|basis| fglm_in(domain, &basis, Order::GradedRevLex, order))
        {
            return Some(converted);
        }
    }
    groebner_in(domain, generators, order, limits)
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

    // -----------------------------------------------------------------
    // Prime fields GF(p)
    // -----------------------------------------------------------------

    use crate::rules::poly::field::Domain;
    use crate::rules::poly::field::oracle;

    use num_traits::One;

    /// The canonical representative of every coefficient must lie in
    /// `0 .. p`.
    fn assert_canonical(
        modulus: i64,
        basis: &[GPoly],
    ) {
        for g in basis {
            assert_eq!(
                g.leading().map(|(_, c)| c),
                Some(&BigRational::from_integer(BigInt::from(1))),
                "a GF(p) basis element is not monic"
            );
            for (_, coefficient) in g.terms() {
                let value = coefficient.numer();
                assert!(
                    *value >= BigInt::from(0) && *value < BigInt::from(modulus),
                    "coefficient {coefficient} is not the canonical representative mod {modulus}"
                );
                assert!(
                    coefficient.denom().is_one(),
                    "a GF(p) coefficient is not integral"
                );
            }
        }
    }

    /// The reference is a reduced basis: Buchberger on it changes nothing.
    fn assert_reference_quality(
        domain: &Domain,
        modulus: i64,
        basis: &[GPoly],
        order: Order,
    ) {
        assert_canonical(modulus, basis);
        // Idempotent: the basis is a reduced Gröbner basis.
        assert_eq!(
            groebner_in(domain, basis, order, GroebnerLimits::default()).as_deref(),
            Some(basis),
            "Buchberger changed the computed basis"
        );
        // Every S-polynomial of the basis reduces to zero.
        for (i, f) in basis.iter().enumerate() {
            for g in basis.iter().skip(i + 1) {
                let s = s_polynomial(domain, f, g, order);
                assert!(
                    reduce_in(domain, &s, basis, order).is_zero(),
                    "an S-polynomial does not reduce to zero"
                );
            }
        }
    }

    #[test]
    fn groebner_over_prime_fields_matches_the_independent_reference() {
        let mut cases = 0_usize;
        for modulus in oracle::PRIMES {
            let domain = Domain::prime(BigInt::from(modulus)).unwrap_or_else(|| {
                panic!("{modulus} must be accepted as a prime");
            });
            for order in [Order::GradedRevLex, Order::GradedLex, Order::Lex] {
                for (name, system) in oracle::sample_systems(modulus) {
                    let generators = oracle::generators(modulus, &system, order);
                    let computed =
                        groebner_via_in(&domain, &generators, order, GroebnerLimits::default())
                            .unwrap_or_else(|| panic!("{name} over GF({modulus}) in {order:?}"));
                    // The reference is a naive Buchberger implementation: in
                    // lex it is hopeless on the random dense systems, so the
                    // random cases are cross-checked in grevlex (where the
                    // engine computes the basis directly) and certified
                    // structurally in the other orders. The structured systems
                    // and every grevlex run are compared literally.
                    let random = name.starts_with("random");
                    if !random || order == Order::GradedRevLex {
                        let expected = oracle::buchberger_fp(modulus, &system, order);
                        assert_eq!(computed, expected, "{name} over GF({modulus}) in {order:?}");
                    }
                    assert_reference_quality(&domain, modulus, &computed, order);
                    // Every generator is in the ideal.
                    for generator in &generators {
                        assert!(
                            reduce_in(&domain, generator, &computed, order).is_zero(),
                            "{name}: a generator is not in the computed ideal"
                        );
                    }
                    cases = cases.saturating_add(1);
                }
            }
        }
        eprintln!("groebner over GF(p): {cases} systems cross-checked");
        assert!(cases >= 15, "the finite-field suite is too small");
    }

    #[test]
    fn direct_buchberger_over_prime_fields_in_every_order() {
        // Direct Buchberger in the target order (no FGLM) on small systems,
        // so the actual elimination is exercised in lex and grlex too.
        let systems: Vec<(&str, Vec<Vec<(i64, Vec<u32>)>>)> = vec![
            ("katsura-2", oracle::katsura(2)),
            ("cyclic-3", oracle::cyclic(3)),
            (
                "linear system",
                vec![
                    vec![(1, vec![1, 0, 0]), (1, vec![0, 1, 0]), (-3, vec![0, 0, 0])],
                    vec![(1, vec![1, 0, 0]), (-1, vec![0, 1, 0]), (-1, vec![0, 0, 0])],
                ],
            ),
            (
                "unit ideal",
                vec![
                    vec![(1, vec![1, 0, 0])],
                    vec![(1, vec![1, 0, 0]), (1, vec![0, 0, 0])],
                ],
            ),
        ];
        for modulus in oracle::PRIMES {
            let domain = Domain::prime(BigInt::from(modulus)).unwrap_or_else(|| {
                panic!("{modulus} must be accepted as a prime");
            });
            for order in [Order::GradedRevLex, Order::GradedLex, Order::Lex] {
                for (name, system) in &systems {
                    let generators = oracle::generators(modulus, system, order);
                    let computed =
                        groebner_in(&domain, &generators, order, GroebnerLimits::default())
                            .unwrap_or_else(|| panic!("{name} over GF({modulus}) in {order:?}"));
                    assert_eq!(
                        computed,
                        oracle::buchberger_fp(modulus, system, order),
                        "{name} over GF({modulus}) in {order:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn groebner_over_prime_fields_is_independent_of_the_generators() {
        // The reduced basis of an ideal does not depend on the presentation;
        // over GF(p) the same must hold, including with redundant and
        // non-monic generators and different input orders.
        for modulus in oracle::PRIMES {
            let domain = Domain::prime(BigInt::from(modulus)).unwrap_or_else(|| {
                panic!("{modulus} must be accepted as a prime");
            });
            for order in [Order::GradedRevLex, Order::GradedLex, Order::Lex] {
                let system = oracle::cyclic(3);
                let generators = oracle::generators(modulus, &system, order);
                let expected = groebner_in(&domain, &generators, order, GroebnerLimits::default())
                    .unwrap_or_default();
                // Add a non-trivial multiple of the first generator, and
                // reorder the list: the reduced basis must not change.
                let two = domain.add(&domain.one(), &domain.one());
                let mut extended = vec![generators[0].clone()];
                extended.extend(generators.iter().cloned());
                let multiple = GPoly::new_in(
                    &domain,
                    generators[0]
                        .terms()
                        .iter()
                        .map(|(m, c)| (m.clone(), domain.mul(c, &two)))
                        .collect(),
                    order,
                );
                extended.push(multiple);
                let computed = groebner_in(&domain, &extended, order, GroebnerLimits::default())
                    .unwrap_or_default();
                assert_eq!(computed, expected, "GF({modulus}) in {order:?}");
            }
        }
    }

    #[test]
    fn coefficients_without_a_prime_field_image_are_rejected() {
        // Over GF(5) the coefficient 1/6 has no image because 6 = 1 mod 5 is
        // invertible... but 1/5 does not, and the request must not silently
        // drop it.
        let domain = Domain::prime(BigInt::from(5)).unwrap();
        let null = GPoly::new(
            vec![(
                vec![1, 0],
                BigRational::new(BigInt::from(1), BigInt::from(5)),
            )],
            Order::Lex,
        );
        assert_eq!(null.in_domain(&domain, Order::Lex), None);
        // 1/2 has an image: 3.
        let half = GPoly::new(
            vec![(
                vec![1, 0],
                BigRational::new(BigInt::from(1), BigInt::from(2)),
            )],
            Order::Lex,
        );
        let converted = half
            .in_domain(&domain, Order::Lex)
            .unwrap_or_else(|| panic!("1/2 is invertible mod 5"));
        assert_eq!(converted.terms(), &[(vec![1, 0], q(3))]);
        // The same coefficient over Q is untouched.
        assert_eq!(
            half.in_domain(&Domain::Rational, Order::Lex),
            Some(half.clone())
        );
    }

    #[test]
    fn fglm_over_prime_fields_matches_direct_buchberger() {
        let mut conversions = 0_usize;
        for modulus in oracle::PRIMES {
            let domain = Domain::prime(BigInt::from(modulus)).unwrap_or_else(|| {
                panic!("{modulus} must be accepted as a prime");
            });
            for system in [
                oracle::katsura(2),
                oracle::katsura(3),
                vec![
                    vec![(1, vec![2, 0]), (1, vec![0, 2]), (-1, vec![0, 0])],
                    vec![(1, vec![1, 0]), (-1, vec![0, 1])],
                ],
                oracle::cyclic(3),
            ] {
                let generators = oracle::generators(modulus, &system, Order::GradedRevLex);
                let graded = groebner_in(
                    &domain,
                    &generators,
                    Order::GradedRevLex,
                    GroebnerLimits::default(),
                )
                .unwrap_or_default();
                for target in [Order::Lex, Order::GradedLex] {
                    // Direct Buchberger in lex is too slow on katsura: the F4
                    // matrix engine (a different algorithm, no FGLM) is the
                    // reference for the target order.
                    let direct = crate::rules::poly::f4::f4_basis_in(
                        &domain,
                        &oracle::generators(modulus, &system, target),
                        target,
                        crate::rules::poly::f4::F4Limits::default(),
                    )
                    .unwrap_or_default();
                    assert_eq!(
                        oracle::buchberger_fp(modulus, &system, Order::GradedRevLex),
                        graded,
                        "the reference disagrees with the grevlex basis over GF({modulus})"
                    );
                    match fglm_in(&domain, &graded, Order::GradedRevLex, target) {
                        | Some(converted) => {
                            assert_eq!(
                                converted, direct,
                                "FGLM over GF({modulus}) into {target:?}"
                            );
                            conversions = conversions.saturating_add(1);
                        },
                        | None => {
                            // Not zero-dimensional in this characteristic (for
                            // example the katsura linear equation degenerates
                            // over GF(2)); the direct basis is still checked.
                        },
                    }
                }
            }
        }
        assert!(
            conversions >= 15,
            "FGLM succeeded only {conversions} times over the prime fields"
        );
    }
}
